"""
ms2_extractor defines functions for generating ms2 files from DDA and PRM based .D folders
"""

import logging
import time
from pathlib import Path
from typing import Generator, List, Optional
import argparse

import pandas as pd
from tdfpy import timsdata
from tdfpy.pandas_tdf import PandasTdf
from serenipy.ms2 import Ms2Spectra
from tqdm import tqdm

from tqdm import tqdm
from .utils import calculate_mass, get_ms2_dda_content, map_precursor_to_ip2_scan_number
import numpy as np

logger = logging.getLogger(__name__)


def write_mgf_file(
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
        output_file = str(Path(analysis_dir) / Path(analysis_dir).stem) + ".mgf"

    logger.info("Generating Ms2 Spectra")
    ms2_spectra = get_ms2_dda_content(
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

    time.sleep(1)

    logger.info("Writing Contents To File")
    with open(output_file, "w", encoding="UTF-8") as file:
        for spectrum in tqdm(ms2_spectra, desc="Writing MGF File"):
            mgf_lines = []
            mgf_lines.append("BEGIN IONS")
            mgf_lines.append(
                f"TITLE={Path(analysis_dir).stem}.{spectrum.low_scan}.{spectrum.high_scan}.{spectrum.charge} "
                f'FILE="{Path(analysis_dir).stem}", NativeID="merged={spectrum.precursor_id} frame={spectrum.parent_id} '
                f'scanStart={spectrum.scan_begin} scanEnd={spectrum.scan_end}"'
            )
            mgf_lines.append(f"RTINSECONDS={spectrum.rt}")
            mgf_lines.append(f"PEPMASS={spectrum.mass}")
            mgf_lines.append(f"CHARGE={spectrum.charge}+")
            for mz, intensity in zip(spectrum.mz_spectra, spectrum.intensity_spectra):
                mgf_lines.append(f"{mz:.5f} {int(intensity)}")
            mgf_lines.append("END IONS")

            file.write("\n".join(mgf_lines) + "\n")

    total_time = round(time.time() - start_time, 2)
    logger.info(f"Total Time: {total_time:.2f} seconds")


def main():
    """
    Command-line interface for MGF extraction from TimsTOF data.
    """
    parser = argparse.ArgumentParser(
        description="Extract MGF files from TimsTOF .D folders",
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
        help="Output MGF file path (default: <analysis_dir_name>.mgf)",
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
        "--min-ccs", type=float, default=None, help="Minimum CCS filter (Ų)"
    )

    parser.add_argument(
        "--max-ccs", type=float, default=None, help="Maximum CCS filter (Ų)"
    )

    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Enable verbose logging"
    )

    parser.add_argument(
        "--casanovo",
        action="store_true",
        help="Remove precursor peaks from MS/MS spectra",
    )

    args = parser.parse_args()

    # if casanovo update params to be: --top-n-spectra 150 --min-intensity 0.01 --min-charge 2 --max-charge 5 --min-mz 50 --max-mz 2500
    if args.casanovo:
        args.remove_precursor = True
        args.precursor_peak_width = 2.0
        args.top_n_spectra = 150
        args.min_intensity = 0.01
        args.min_mz = 50
        args.max_mz = 2500

    # Set up logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=log_level, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )

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
    else:
        d_folders = list(analysis_path.glob("*.d"))
    if not d_folders:
        logger.error(f"No .d folders found in: {args.analysis_dir}")
        return 1

    output = args.output
    if len(d_folders) > 1:
        output = None

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

        try:
            write_mgf_file(
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
            logger.info("MGF extraction completed successfully!")
            return 0
        except Exception as e:
            logger.error(f"Error during MGF extraction: {e}")
            return 1


if __name__ == "__main__":
    exit(main())
