"""
ms2_extractor defines functions for generating ms2 files from DDA and PRM based .D folders
"""

import logging
import os
import time
from datetime import datetime
from pathlib import Path
from typing import Generator, Optional
import argparse

import pandas as pd
from tdfpy import timsdata
from tdfpy.pandas_tdf import PandasTdf
from serenipy.ms2 import Ms2Spectra, to_ms2
from tqdm import tqdm

from .constants import MS2_VERSION
from .utils import calculate_mass, get_ms2_dda_content, map_precursor_to_ip2_scan_number

logger = logging.getLogger(__name__)
# make debug
logger.setLevel(logging.DEBUG)

def get_ms2_content(
    analysis_dir: str,
    remove_precursor: bool = False,
    precursor_peak_width: float = 2.0,
    batch_size: int = 100,
    top_n_spectra: Optional[int] = None,
    min_intensity: float = 0.0,
    min_charge: Optional[int] = None,
    max_charge: Optional[int] = None,
    min_mz: Optional[float] = None,
    max_mz: Optional[float] = None,
    min_rt: Optional[float] = None,
    max_rt: Optional[float] = None,
    min_ccs: Optional[float] = None,
    max_ccs: Optional[float] = None,
) -> Generator[Ms2Spectra, None, None]:

    pd_tdf = PandasTdf(str(Path(analysis_dir) / "analysis.tdf"))
    if pd_tdf.is_dda:
        logger.info("TDF format is DDA")
        return get_ms2_dda_content(
            analysis_dir=analysis_dir,
            remove_precursor=remove_precursor,
            precursor_peak_width=precursor_peak_width,
            batch_size=batch_size,
            top_n_spectra=top_n_spectra,
            min_intensity=min_intensity,
            min_charge=min_charge,
            max_charge=max_charge,
            min_mz=min_mz,
            max_mz=max_mz,
            min_rt=min_rt,
            max_rt=max_rt,
            min_ccs=min_ccs,
            max_ccs=max_ccs,
        )
    if pd_tdf.is_prm:
        logger.info("TDF format is PRM")
        return get_ms2_prm_content(
            analysis_dir=analysis_dir,
            remove_precursor=remove_precursor,
            precursor_peak_width=precursor_peak_width,
            batch_size=batch_size,
            top_n_spectra=top_n_spectra,
            min_intensity=min_intensity,
            min_charge=min_charge,
            max_charge=max_charge,
            min_mz=min_mz,
            max_mz=max_mz,
            min_rt=min_rt,
            max_rt=max_rt,
            min_ccs=min_ccs,
            max_ccs=max_ccs,
        )

    raise TypeError("Unknown TDF format")


def get_ms2_prm_content(
    analysis_dir: str,
    remove_precursor: bool = False,
    precursor_peak_width: float = 2.0,
    batch_size: int = 100,
    top_n_spectra: Optional[int] = None,
    min_intensity: float = 0.0,
    min_charge: Optional[int] = None,
    max_charge: Optional[int] = None,
    min_mz: Optional[float] = None,
    max_mz: Optional[float] = None,
    min_rt: Optional[float] = None,
    max_rt: Optional[float] = None,
    min_ccs: Optional[float] = None,
    max_ccs: Optional[float] = None,
) -> Generator[Ms2Spectra, None, None]:

    with timsdata.timsdata_connect(analysis_dir) as td:
        analysis_tdf_path = str(Path(analysis_dir) / "analysis.tdf")
        merged_df = pd.merge(
            PandasTdf(analysis_tdf_path).prm_frame_msms_info,
            PandasTdf(analysis_tdf_path).prm_targets,
            left_on="Target",
            right_on="Id",
        )
        merged_df = pd.merge(
            merged_df,
            PandasTdf(analysis_tdf_path).frames,
            left_on="Frame",
            right_on="Id",
            suffixes=("_Target", "_Frame"),
        )

        for _, row in tqdm(
            merged_df.iterrows(), desc="Generating MS2 Spectra", total=len(merged_df)
        ):

            if min_charge is not None and int(row["Charge"]) < min_charge:
                continue
            if max_charge is not None and int(row["Charge"]) > max_charge:
                continue

            # Apply m/z filters
            if min_mz is not None and float(row["IsolationMz"]) < min_mz:
                continue
            if max_mz is not None and float(row["IsolationMz"]) > max_mz:
                continue

            # Apply RT filters
            if min_rt is not None and float(row["Time_Frame"]) < min_rt:
                continue
            if max_rt is not None and float(row["Time_Frame"]) > max_rt:
                continue

            mz_list, area_list = td.extractCentroidedSpectrumForFrame(
                frame_id=int(row["Frame"]),
                scan_begin=int(row["ScanNumBegin"]),
                scan_end=int(row["ScanNumEnd"]),
                peak_picker_resolution=120000,
            )

            ms2_spectra = Ms2Spectra(
                low_scan=int(row["Frame"]),
                high_scan=int(row["Frame"]),
                mz=float(row["IsolationMz"]),
                mass=calculate_mass(float(row["IsolationMz"]), int(row["Charge"])),
                charge=int(row["Charge"]),
                info={
                    "Target": str(int(row["Target"])),
                    "Accumulation_Time": str(float(row["AccumulationTime"])),
                    "Ramp_Time": str(float(row["RampTime"])),
                    "Pressure": str(float(row["Pressure"])),
                },
                mz_spectra=[],
                intensity_spectra=[],
                charge_spectra=[],
            )

            ms2_spectra.ce = float(row["CollisionEnergy"])
            ms2_spectra.iso_width = float(row["IsolationWidth"])
            ms2_spectra.iso_mz = float(row["IsolationMz"])
            ms2_spectra.rt = float(row["Time_Frame"])
            ms2_spectra.scan_begin = float(row["ScanNumBegin"])
            ms2_spectra.scan_end = float(row["ScanNumEnd"])

            ook0_range = td.scanNumToOneOverK0(
                int(row["Frame"]), [ms2_spectra.scan_begin, ms2_spectra.scan_end]
            )
            ms2_spectra.info["OOK0_Begin"] = ook0_range[0]
            ms2_spectra.info["OOK0_End"] = ook0_range[1]

            ms2_spectra_data = list(zip(list(mz_list), list(area_list)))

            # Remove precursor peaks if requested
            if remove_precursor:
                precursor_mz = float(row["IsolationMz"])
                ms2_spectra_data = [
                    data
                    for data in ms2_spectra_data
                    if abs(data[0] - precursor_mz) > precursor_peak_width
                ]

            if min_intensity != 0:
                ms2_spectra_data = [
                    data for data in ms2_spectra_data if data[1] >= min_intensity
                ]

            # Sort by intensity and keep top N if specified
            if top_n_spectra is not None:
                ms2_spectra_data.sort(key=lambda x: x[1], reverse=True)
                ms2_spectra_data = ms2_spectra_data[:top_n_spectra]

            ms2_spectra.mz_spectra = [data[0] for data in ms2_spectra_data]
            ms2_spectra.intensity_spectra = [int(data[1]) for data in ms2_spectra_data]

            assert len(ms2_spectra.mz_spectra) == len(ms2_spectra.intensity_spectra)

            if len(ms2_spectra.mz_spectra) == 0:
                continue

            yield ms2_spectra


def generate_header(
    analysis_dir: str,
    output_file: Optional[str] = None,
    remove_precursor: bool = False,
    precursor_peak_width: float = 2.0,
    batch_size: int = 100,
    top_n_spectra: Optional[int] = None,
    min_intensity: float = 0.0,
    min_charge: Optional[int] = None,
    max_charge: Optional[int] = None,
    min_mz: Optional[float] = None,
    max_mz: Optional[float] = None,
    min_rt: Optional[float] = None,
    max_rt: Optional[float] = None,
    min_ccs: Optional[float] = None,
    max_ccs: Optional[float] = None,
    resolution: float = 120000,
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
        "H\tMinimumMsMsIntensity\t{min_intensity}\n"
        "H\tRemovePrecursor\t{remove_precursor}\n"
        "H\tPrecursorPeakWidth\t{precursor_peak_width}\n"
        "H\tBatchSize\t{batch_size}\n"
        "H\tTopNSpectra\t{top_n_spectra}\n"
        "H\tMinCharge\t{min_charge}\n"
        "H\tMaxCharge\t{max_charge}\n"
        "H\tMinMz\t{min_mz}\n"
        "H\tMaxMz\t{max_mz}\n"
        "H\tMinRt\t{min_rt}\n"
        "H\tMaxRt\t{max_rt}\n"
        "H\tMinCcs\t{min_ccs}\n"
        "H\tMaxCcs\t{max_ccs}\n"
        "H\tExtractorOptions\tMSn\n"
        "H\tAcquisitionMethod\t{method}\n"
        "H\tInstrumentType\tTIMSTOF\n"
        "H\tDataType\tCentroid\n"
        "H\tScanType\tMS2\n"
        "H\tResolution\t{resolution}\n"
        "H\tIsolationWindow\n"
        "H\tFirstScan\t{first_scan:d}\n"
        "H\tLastScan\t{last_scan:d}\n"
        "H\tMonoIsotopic PrecMz\tTrue\n"
    )

    ms2_header = MS2_HEADER.format(
        version=MS2_VERSION,
        date_of_creation=str(datetime.now().strftime("%B %d, %Y %H:%M")),
        min_intensity=min_intensity,
        remove_precursor=remove_precursor,
        precursor_peak_width=precursor_peak_width,
        batch_size=batch_size,
        top_n_spectra=top_n_spectra if top_n_spectra is not None else "None",
        min_charge=min_charge if min_charge is not None else "None",
        max_charge=max_charge if max_charge is not None else "None",
        min_mz=min_mz if min_mz is not None else "None",
        max_mz=max_mz if max_mz is not None else "None",
        min_rt=min_rt if min_rt is not None else "None",
        max_rt=max_rt if max_rt is not None else "None",
        min_ccs=min_ccs if min_ccs is not None else "None",
        max_ccs=max_ccs if max_ccs is not None else "None",
        method=method,
        resolution=resolution,
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
    top_n_spectra: Optional[int] = None,
    min_intensity: float = 0.0,
    min_charge: Optional[int] = None,
    max_charge: Optional[int] = None,
    min_mz: Optional[float] = None,
    max_mz: Optional[float] = None,
    min_rt: Optional[float] = None,
    max_rt: Optional[float] = None,
    min_ccs: Optional[float] = None,
    max_ccs: Optional[float] = None,
):

    start_time = time.time()

    if output_file is None:
        output_file = str(Path(analysis_dir) / Path(analysis_dir).stem) + ".ms2"

    logger.info("Creating Ms2 Header")
    ms2_header = generate_header(
        analysis_dir=analysis_dir,
        min_intensity=min_intensity,
        min_charge=min_charge,
        resolution=120000,
    )

    logger.info("Generating Ms2 Spectra")
    ms2_spectra = get_ms2_content(
        analysis_dir=analysis_dir,
        batch_size=batch_size,
        min_intensity=min_intensity,
        remove_precursor=remove_precursor,
        precursor_peak_width=precursor_peak_width,
        top_n_spectra=top_n_spectra,
        min_charge=min_charge,
        max_charge=max_charge,
        min_mz=min_mz,
        max_mz=max_mz,
        min_rt=min_rt,
        max_rt=max_rt,
        min_ccs=min_ccs,
        max_ccs=max_ccs,
    )

    logger.info("Creating MS2 Contents")
    ms2_content = to_ms2([ms2_header], ms2_spectra)

    time.sleep(1)  # Dumb fix for logging message overlap with tqdm progress bar

    logger.info("Writing Contents To File")
    with open(output_file, "w", encoding="UTF-8") as file:
        file.write(ms2_content)

    total_time = round(time.time() - start_time, 2)
    logger.info(f"Total Time: {total_time:.2f} seconds")


def main():
    """
    Command-line interface for MS2 extraction from TimsTOF data.
    """
    parser = argparse.ArgumentParser(
        description="Extract MS2 files from TimsTOF .D folders",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "analysis_dir", type=str, help="Path to the .D analysis directory"
    )

    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default=None,
        help="Output MS2 file path (default: <analysis_dir_name>.ms2)",
    )

    parser.add_argument(
        "--remove-precursor",
        action="store_true",
        help="Remove precursor peaks from MS/MS spectra",
    )

    parser.add_argument(
        "--precursor-peak-width",
        type=float,
        default=2.0,
        help="Width around precursor m/z to remove (Da)",
    )

    parser.add_argument(
        "--batch-size", type=int, default=100, help="Batch size for processing spectra"
    )

    parser.add_argument(
        "--top-n-spectra",
        type=int,
        default=None,
        help="Keep only top N most intense peaks per spectrum",
    )

    parser.add_argument(
        "--min-intensity",
        type=float,
        default=0.0,
        help="Minimum intensity threshold for peaks",
    )

    parser.add_argument(
        "--min-charge", type=int, default=None, help="Minimum charge state filter"
    )

    parser.add_argument(
        "--max-charge", type=int, default=None, help="Maximum charge state filter"
    )

    parser.add_argument("--min-mz", type=float, default=None, help="Minimum m/z filter")

    parser.add_argument("--max-mz", type=float, default=None, help="Maximum m/z filter")

    parser.add_argument(
        "--min-rt",
        type=float,
        default=None,
        help="Minimum retention time filter (seconds)",
    )

    parser.add_argument(
        "--max-rt",
        type=float,
        default=None,
        help="Maximum retention time filter (seconds)",
    )

    parser.add_argument(
        "--min-ccs", type=float, default=None, help="Minimum CCS filter"
    )

    parser.add_argument(
        "--max-ccs", type=float, default=None, help="Maximum CCS filter"
    )

    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Enable verbose logging"
    )

    parser.add_argument(
        "--overwrite",
        action="store_true", 
        help="Overwrite existing output file if it exists",
    )

    args = parser.parse_args()

    # Set up logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=log_level, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )

    # Print all arguments being used
    logger.info("MS2 Extractor Arguments:")
    logger.info(f"  Analysis Directory: {args.analysis_dir}")
    logger.info(f"  Output File: {args.output if args.output else 'Within .d folder'}")
    logger.info(f"  Remove Precursor: {args.remove_precursor}")
    logger.info(f"  Precursor Peak Width: {args.precursor_peak_width} Da")
    logger.info(f"  Batch Size: {args.batch_size}")
    logger.info(f"  Top N Spectra: {args.top_n_spectra if args.top_n_spectra else 'All'}")
    logger.info(f"  Min Intensity: {args.min_intensity}")
    logger.info(f"  Min Charge: {args.min_charge if args.min_charge else 'None'}")
    logger.info(f"  Max Charge: {args.max_charge if args.max_charge else 'None'}")
    logger.info(f"  Min m/z: {args.min_mz if args.min_mz else 'None'}")
    logger.info(f"  Max m/z: {args.max_mz if args.max_mz else 'None'}")
    logger.info(f"  Min RT: {args.min_rt if args.min_rt else 'None'} seconds")
    logger.info(f"  Max RT: {args.max_rt if args.max_rt else 'None'} seconds")
    logger.info(f"  Min CCS: {args.min_ccs if args.min_ccs else 'None'}")
    logger.info(f"  Max CCS: {args.max_ccs if args.max_ccs else 'None'}")
    logger.info(f"  Overwrite Existing Output: {args.overwrite}")
    logger.info(f"  Verbose Logging: {args.verbose}")

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
            raise ValueError(
                "Output file specified but multiple .d folders found.")
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
        _output_name = output_name if output_name is not None else Path(d_folder).stem + ".ms2"

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
                top_n_spectra=args.top_n_spectra,
                min_intensity=args.min_intensity,
                min_charge=args.min_charge,
                max_charge=args.max_charge,
                min_mz=args.min_mz,
                max_mz=args.max_mz,
                min_rt=args.min_rt,
                max_rt=args.max_rt,
                min_ccs=args.min_ccs,
                max_ccs=args.max_ccs,
            )
            logger.info("MS2 extraction completed successfully!")
        except Exception as e:
            logger.error(f"Error during Ms2 extraction: {e}... skipping {d_folder}")
        except KeyboardInterrupt:
            logger.info("Extraction interrupted by user.")
            return 0
        
if __name__ == "__main__":
    exit(main())
