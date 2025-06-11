"""
ms2_extractor defines functions for generating ms2 files from DDA and PRM based .D folders
"""

import logging
import os
import time
from datetime import datetime
from pathlib import Path
from typing import Generator, Optional
import threading
import queue

import pandas as pd
from tdfpy import timsdata
from tdfpy.pandas_tdf import PandasTdf
from serenipy.ms2 import Ms2Spectra
from tqdm import tqdm

from .constants import MS2_VERSION
from .utils import calculate_p1mass, get_ms2_dda_content, get_tdf_df, map_precursor_to_ip2_scan_number
from .cli_args import create_ms2_parser, apply_preset_settings, log_common_args

logger = logging.getLogger(__name__)
# make debug
logger.setLevel(logging.DEBUG)


def generate_header(
    analysis_dir: str,
    remove_precursor: bool = False,
    precursor_peak_width: float = 2.0,
    batch_size: int = 100,
    top_n_peaks: Optional[int] = None,
    min_spectra_intensity: Optional[float] = None,
    max_spectra_intensity: Optional[float] = None,
    min_spectra_mz: Optional[float] = None,
    max_spectra_mz: Optional[float] = None,
    min_precursor_intensity: Optional[int] = None,
    max_precursor_intensity: Optional[int] = None,
    min_precursor_charge: Optional[int] = None,
    max_precursor_charge: Optional[int] = None,
    min_precursor_mz: Optional[float] = None,
    max_precursor_mz: Optional[float] = None,
    min_precursor_rt: Optional[float] = None,
    max_precursor_rt: Optional[float] = None,
    min_precursor_ccs: Optional[float] = None,
    max_precursor_ccs: Optional[float] = None,
    min_precursor_neutral_mass: Optional[float] = None,
    max_precursor_neutral_mass: Optional[float] = None,
    mz_precision: int = 5,
    intensity_precision: int = 0,
):
    """
    Generates a header string for MS2 data using information from the analysis file.

    Args:
        analysis_dir (str): The directory path where the analysis file `analysis.tdf` is located.
        min_intensity (float): The minimum intensity of MS/MS spectra to be considered.
        remove_charge1 (bool): Indicates whether to exclude precursors with a charge of 1.
        remove_empty_spectra (bool): Indicates whether to exclude MS/MS spectra with no spectra values.
        include_spectra (bool): Indicates whether to include the actual MS/MS spectra in the output.
        resolution (float): The resolution of the MS/MS spectra.

    Returns:
        str: The generated header string for MS2 data.

    Raises:
        TypeError: If the TDF format is unknown.
    """

    pd_tdf = PandasTdf(str(Path(analysis_dir) / "analysis.tdf"))

    if pd_tdf.is_dda:
        method = "Data-Dependent"
        precursors_df = pd_tdf.precursors
        frames_df = pd_tdf.frames
        precursor_to_scan_number = map_precursor_to_ip2_scan_number(
            precursors_df, frames_df
        )
        first_scan = list(precursor_to_scan_number.values())[0]
        last_scan = list(precursor_to_scan_number.values())[-1]
    elif pd_tdf.is_prm:
        method = "Parallel-Reaction-Monitoring"
        frames_df = pd_tdf.frames
        frames_df = frames_df[frames_df["MsMsType"] == 10]
        first_scan = min(frames_df["Id"])
        last_scan = max(frames_df["Id"])
    else:
        raise TypeError("Unknown TDF format")

    MS2_HEADER = (
        "H\tExtractor\tTimsTOF_extractor\n"
        "H\tExtractorVersion\t{version}\n"
        "H\tPublicationDate\t20-02-2020\n"
        "H\tDate Converted\t{date_of_creation}\n"
        "H\tComments\tTimsTOF_extractor written by Yu Gao, 2018\n"
        "H\tComments\tTimsTOF_extractor modified by Titus Jung, 2019\n"
        "H\tComments\tTimsTOF_extractor modified by Patrick Garrett, 2025\n"
        "H\tMinSpectraIntensity\t{min_spectra_intensity}\n"
        "H\tMaxSpectraIntensity\t{max_spectra_intensity}\n"
        "H\tMinSpectraMz\t{min_spectra_mz}\n"
        "H\tMaxSpectraMz\t{max_spectra_mz}\n"
        "H\tMinPrecursorIntensity\t{min_precursor_intensity}\n"
        "H\tMaxPrecursorIntensity\t{max_precursor_intensity}\n"
        "H\tRemovePrecursor\t{remove_precursor}\n"
        "H\tPrecursorPeakWidth\t{precursor_peak_width}\n"
        "H\tBatchSize\t{batch_size}\n"
        "H\tTopNPeaks\t{top_n_peaks}\n"
        "H\tMinPrecursorCharge\t{min_precursor_charge}\n"
        "H\tMaxPrecursorCharge\t{max_precursor_charge}\n"
        "H\tMinPrecursorMz\t{min_precursor_mz}\n"
        "H\tMaxPrecursorMz\t{max_precursor_mz}\n"
        "H\tMinPrecursorRt\t{min_precursor_rt}\n"
        "H\tMaxPrecursorRt\t{max_precursor_rt}\n"
        "H\tMinPrecursorCcs\t{min_precursor_ccs}\n"
        "H\tMaxPrecursorCcs\t{max_precursor_ccs}\n"
        "H\tMinPrecursorNeutralMass\t{min_precursor_neutral_mass}\n"
        "H\tMaxPrecursorNeutralMass\t{max_precursor_neutral_mass}\n"
        "H\tExtractorOptions\tMSn\n"
        "H\tAcquisitionMethod\t{method}\n"
        "H\tInstrumentType\tTIMSTOF\n"
        "H\tDataType\tCentroid\n"
        "H\tScanType\tMS2\n"
        "H\tResolution\t{resolution}\n"
        "H\tMzPrecision\t{mz_precision}\n"
        "H\tIntensityPrecision\t{intensity_precision}\n"
        "H\tIsolationWindow\n"
        "H\tFirstScan\t{first_scan:d}\n"
        "H\tLastScan\t{last_scan:d}\n"
        "H\tMonoIsotopic PrecMz\tTrue\n"
    )

    ms2_header = MS2_HEADER.format(
        version=MS2_VERSION,
        date_of_creation=str(datetime.now().strftime("%B %d, %Y %H:%M")),
        min_spectra_intensity=(
            min_spectra_intensity if min_spectra_intensity is not None else "None"
        ),
        max_spectra_intensity=(
            max_spectra_intensity if max_spectra_intensity is not None else "None"
        ),
        min_spectra_mz=min_spectra_mz if min_spectra_mz is not None else "None",
        max_spectra_mz=max_spectra_mz if max_spectra_mz is not None else "None",
        min_precursor_intensity=(
            min_precursor_intensity if min_precursor_intensity is not None else "None"
        ),
        max_precursor_intensity=(
            max_precursor_intensity if max_precursor_intensity is not None else "None"
        ),
        remove_precursor=remove_precursor,
        precursor_peak_width=precursor_peak_width,
        batch_size=batch_size,
        top_n_peaks=top_n_peaks if top_n_peaks is not None else "None",
        min_precursor_charge=(
            min_precursor_charge if min_precursor_charge is not None else "None"
        ),
        max_precursor_charge=(
            max_precursor_charge if max_precursor_charge is not None else "None"
        ),
        min_precursor_mz=min_precursor_mz if min_precursor_mz is not None else "None",
        max_precursor_mz=max_precursor_mz if max_precursor_mz is not None else "None",
        min_precursor_rt=min_precursor_rt if min_precursor_rt is not None else "None",
        max_precursor_rt=max_precursor_rt if max_precursor_rt is not None else "None",
        min_precursor_ccs=(
            min_precursor_ccs if min_precursor_ccs is not None else "None"
        ),
        max_precursor_ccs=(
            max_precursor_ccs if max_precursor_ccs is not None else "None"
        ),
        min_precursor_neutral_mass=(
            min_precursor_neutral_mass
            if min_precursor_neutral_mass is not None
            else "None"
        ),
        max_precursor_neutral_mass=(
            max_precursor_neutral_mass
            if max_precursor_neutral_mass is not None
            else "None"
        ),
        method=method,
        resolution=120000,
        mz_precision=mz_precision,
        intensity_precision=intensity_precision,
        first_scan=first_scan,
        last_scan=last_scan,
    )

    return ms2_header


def write_ms2_file(
    analysis_dir: str,
    output_file: Optional[str] = None,
    remove_precursor: bool = False,
    precursor_peak_width: float = 2.0,
    batch_size: int = 100,
    top_n_peaks: Optional[int] = None,
    min_spectra_intensity: Optional[float] = None,
    max_spectra_intensity: Optional[float] = None,
    min_spectra_mz: Optional[float] = None,
    max_spectra_mz: Optional[float] = None,
    min_precursor_intensity: Optional[int] = None,
    max_precursor_intensity: Optional[int] = None,
    min_precursor_charge: Optional[int] = None,
    max_precursor_charge: Optional[int] = None,
    min_precursor_mz: Optional[float] = None,
    max_precursor_mz: Optional[float] = None,
    min_precursor_rt: Optional[float] = None,
    max_precursor_rt: Optional[float] = None,
    min_precursor_ccs: Optional[float] = None,
    max_precursor_ccs: Optional[float] = None,
    min_precursor_neutral_mass: Optional[float] = None,
    max_precursor_neutral_mass: Optional[float] = None,
    mz_precision: int = 5,
    intensity_precision: int = 0,
    keep_empty_spectra: bool = False,
):

    start_time = time.time()

    if output_file is None:
        output_file = str(Path(analysis_dir) / Path(analysis_dir).stem) + ".ms2"

    logger.info("Creating Ms2 Header")
    ms2_header = generate_header(
        analysis_dir=analysis_dir,
        remove_precursor=remove_precursor,
        precursor_peak_width=precursor_peak_width,
        batch_size=batch_size,
        top_n_peaks=top_n_peaks,
        min_spectra_intensity=min_spectra_intensity,
        max_spectra_intensity=max_spectra_intensity,
        min_spectra_mz=min_spectra_mz,
        max_spectra_mz=max_spectra_mz,
        min_precursor_intensity=min_precursor_intensity,
        max_precursor_intensity=max_precursor_intensity,
        min_precursor_charge=min_precursor_charge,
        max_precursor_charge=max_precursor_charge,
        min_precursor_mz=min_precursor_mz,
        max_precursor_mz=max_precursor_mz,
        min_precursor_rt=min_precursor_rt,
        max_precursor_rt=max_precursor_rt,
        min_precursor_ccs=min_precursor_ccs,
        max_precursor_ccs=max_precursor_ccs,
        min_precursor_neutral_mass=min_precursor_neutral_mass,
        max_precursor_neutral_mass=max_precursor_neutral_mass,
        mz_precision=mz_precision,
        intensity_precision=intensity_precision,
    )

    merged_df = get_tdf_df(
        analysis_dir,
        min_precursor_intensity,
        max_precursor_intensity,
        min_precursor_charge,
        max_precursor_charge,
        min_precursor_mz,
        max_precursor_mz,
        min_precursor_rt,
        max_precursor_rt,
        min_precursor_ccs,
        max_precursor_ccs,
        min_precursor_neutral_mass,
        max_precursor_neutral_mass
    )

    logger.info("Generating Ms2 Spectra (producer-consumer mode)")
    spectra_queue = queue.Queue(maxsize=100)

    def producer():
        try:
            ms2_spectra = get_ms2_dda_content(
                analysis_dir=analysis_dir,
                merged_df=merged_df,
                remove_precursor=remove_precursor,
                precursor_peak_width=precursor_peak_width,
                batch_size=batch_size,
                top_n_peaks=top_n_peaks,
                min_spectra_intensity=min_spectra_intensity,
                max_spectra_intensity=max_spectra_intensity,
                min_spectra_mz=min_spectra_mz,
                max_spectra_mz=max_spectra_mz,
            )
            for spectrum in ms2_spectra:
                spectra_queue.put(spectrum)
        finally:
            spectra_queue.put(None)  # Sentinel value

    def consumer():
        with open(output_file, "w", encoding="UTF-8") as file:
            file.write(ms2_header)
            with tqdm(desc="Writing MS2 Spectra", unit="spectra", total=len(merged_df)) as pbar:
                while True:
                    ms2_spectra = spectra_queue.get()

                    if ms2_spectra is None:
                        break

                    pbar.update(1)

                    if len(ms2_spectra.mz_spectra) == 0 and keep_empty_spectra is False:
                        continue

                    file.write(ms2_spectra.serialize(
                        mz_precision=mz_precision,
                        intensity_precision=intensity_precision,
                    ))

    producer_thread = threading.Thread(target=producer)
    consumer_thread = threading.Thread(target=consumer)

    producer_thread.start()
    consumer_thread.start()
    producer_thread.join()
    consumer_thread.join()

    total_time = round(time.time() - start_time, 2)
    logger.info(f"Total Time: {total_time:.2f} seconds")

def main():
    """
    Command-line interface for MGF extraction from TimsTOF data.
    """
    
    parser = create_ms2_parser()
    args = parser.parse_args()

    # Set up logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=log_level, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )

    # Apply preset settings
    apply_preset_settings(args)

    # Log all arguments being used
    log_common_args(logger, args, "MS2")

    # Validate input directory
    analysis_path = Path(args.analysis_dir)
    if not analysis_path.exists():
        logger.error(f"Analysis directory does not exist: {args.analysis_dir}")
        return 1

    if not analysis_path.is_dir():
        logger.error(f"Path is not a directory: {args.analysis_dir}")
        return 1

    d_folders = []

    # check if it ends in .d else look for all .d directories
    if analysis_path.name.endswith(".d"):
        d_folders.append(analysis_path)
        logger.info(f"Using provided .d folder: {analysis_path}")
    else:
        d_folders = list(analysis_path.glob("*.d"))
        logger.info(f"Found {len(d_folders)} .d folders in: {args.analysis_dir}")
        logger.debug(f"Found .d folders: {d_folders}")

    if not d_folders:
        logger.error(f"No .d folders found in: {args.analysis_dir}")
        return 1

    output_dir = None
    output_name = None

    # if output is a dir
    if args.output is None:
        # output will bewithin d folder
        output_dir = None
        output_name = None

    elif args.output.endswith(".ms2"):
        if len(d_folders) > 1:
            raise ValueError("Output file specified but multiple .d folders found.")
        output_dir = Path(args.output).parent
        output_name = Path(args.output).name

    else:
        # path is a dir
        output_dir = Path(args.output)

        # make dir if it does not exist
        if not output_dir.exists():
            try:
                output_dir.mkdir(parents=True, exist_ok=True)
                logger.info(f"Created output directory: {output_dir}")
            except Exception as e:
                logger.error(f"Failed to create output directory: {e}")
                return 1

        output_name = None

    for d_folder in d_folders:
        if not d_folder.is_dir():
            logger.error(f"Path is not a directory: {d_folder}")
            return 1

        if not (d_folder / "analysis.tdf").exists():
            logger.error(f"Required file not found in {d_folder}: analysis.tdf")
            return 1
        if not (d_folder / "analysis.tdf_bin").exists():
            logger.error(f"Required file not found in {d_folder}: analysis.tdf_bin")
            return 1

        logger.info(f"Processing {d_folder}...")

        _output_dir = output_dir if output_dir is not None else d_folder
        _output_name = (
            output_name if output_name is not None else Path(d_folder).stem + ".ms2"
        )

        output = os.path.join(_output_dir, _output_name)
        logger.info(f"Output file: {output}")

        if not args.overwrite and Path(output).exists():
            logger.warning(f"Output file {output} already exists. Skipping...")
            continue

        try:

            write_ms2_file(
                analysis_dir=str(d_folder),
                output_file=output,
                remove_precursor=args.remove_precursor,
                precursor_peak_width=args.precursor_peak_width,
                batch_size=args.batch_size,
                top_n_peaks=args.top_n_peaks,
                min_spectra_intensity=args.min_spectra_intensity,
                max_spectra_intensity=args.max_spectra_intensity,
                min_spectra_mz=args.min_spectra_mz,
                max_spectra_mz=args.max_spectra_mz,
                min_precursor_intensity=args.min_precursor_intensity,
                max_precursor_intensity=args.max_precursor_intensity,
                min_precursor_charge=args.min_precursor_charge,
                max_precursor_charge=args.max_precursor_charge,
                min_precursor_mz=args.min_precursor_mz,
                max_precursor_mz=args.max_precursor_mz,
                min_precursor_rt=args.min_precursor_rt,
                max_precursor_rt=args.max_precursor_rt,
                min_precursor_ccs=args.min_precursor_ccs,
                max_precursor_ccs=args.max_precursor_ccs,
                min_precursor_neutral_mass=args.min_precursor_neutral_mass,
                max_precursor_neutral_mass=args.max_precursor_neutral_mass,
                keep_empty_spectra=args.keep_empty_spectra,
            )
            logger.info("MS2 extraction completed successfully!")
        except Exception as e:
            logger.error(f"Error during Ms2 extraction: {e}... skipping {d_folder}")
            logger.error(e, exc_info=True)
        except KeyboardInterrupt:
            logger.info("Extraction interrupted by user.")
            os._exit(0)


if __name__ == "__main__":
    exit(main())
