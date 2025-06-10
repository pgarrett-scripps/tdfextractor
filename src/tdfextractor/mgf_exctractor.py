"""
ms2_extractor defines functions for generating ms2 files from DDA and PRM based .D folders
"""

import logging
import os
import time
from pathlib import Path
from typing import Optional

from tqdm import tqdm

from .utils import get_ms2_dda_content
from .cli_args import create_mgf_parser, apply_preset_settings, log_common_args

logger = logging.getLogger(__name__)


def write_mgf_file(
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
):

    start_time = time.time()

    if output_file is None:
        output_file = str(Path(analysis_dir) / Path(analysis_dir).stem) + ".mgf"

    logger.info("Generating Ms2 Spectra")
    ms2_spectra = get_ms2_dda_content(
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
    )

    time.sleep(1)

    logger.info("Writing Contents To File")
    with open(output_file, "w", encoding="UTF-8") as file:
        for spectrum in tqdm(ms2_spectra, desc="Writing MGF File"):
            mgf_lines = []
            mgf_lines.append("BEGIN IONS")
            mgf_lines.append(
                f"TITLE={Path(analysis_dir).stem}.{spectrum.low_scan}.{spectrum.high_scan}.{spectrum.charge} "
                f'File="{Path(analysis_dir).stem}", NativeID="merged={spectrum.precursor_id} frame={spectrum.parent_id} '
                f'scanStart={spectrum.scan_begin} scanEnd={spectrum.scan_end} scan={spectrum.low_scan}"'
            )
            mgf_lines.append(f"RTINSECONDS={spectrum.rt}")
            mgf_lines.append(f"PEPMASS={spectrum.mass}")
            mgf_lines.append(f"CHARGE={spectrum.charge}+")
            for mz, intensity in zip(spectrum.mz_spectra, spectrum.intensity_spectra):
                mgf_lines.append(f"{mz:.5f} {int(intensity)}")
            mgf_lines.append("END IONS")

            file.write("\n".join(mgf_lines) + "\n\n")

    total_time = round(time.time() - start_time, 2)
    logger.info(f"Total Time: {total_time:.2f} seconds")


def main():
    """
    Command-line interface for MGF extraction from TimsTOF data.
    """
    parser = create_mgf_parser()
    args = parser.parse_args()

    # Set up logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=log_level, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )

    # Apply preset settings
    apply_preset_settings(args)

    # Log all arguments being used
    log_common_args(logger, args, "MGF")

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

    elif args.output.endswith(".mgf"):
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
            output_name if output_name is not None else Path(d_folder).stem + ".mgf"
        )

        output = os.path.join(_output_dir, _output_name)
        logger.info(f"Output file: {output}")

        if not args.overwrite and Path(output).exists():
            logger.warning(f"Output file {output} already exists. Skipping...")
            continue

        try:

            write_mgf_file(
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
            )
            logger.info("MGF extraction completed successfully!")
        except Exception as e:
            logger.error(f"Error during MGF extraction: {e}... skipping {d_folder}")
        except KeyboardInterrupt:
            logger.info("Extraction interrupted by user.")
            return 0


if __name__ == "__main__":
    exit(main())
