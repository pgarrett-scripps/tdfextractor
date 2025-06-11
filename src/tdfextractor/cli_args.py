"""
Shared command line argument definitions for MS2 and MGF extractors.
"""

import argparse
from typing import Optional


def add_common_args(parser: argparse.ArgumentParser) -> None:
    """Add common arguments shared between MS2 and MGF extractors."""

    parser.add_argument(
        "analysis_dir", type=str, help="Path to the .D analysis directory"
    )

    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default=None,
        help="Output file path (default: <analysis_dir_name>.<ext>)",
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
        "--top-n-peaks",
        type=int,
        default=None,
        help="Keep only top N most intense peaks per spectrum",
    )

    parser.add_argument(
        "--min-spectra-intensity",
        type=float,
        default=None,
        help="Minimum intensity threshold for MS/MS peaks (absolute or 0.0-1.0 for percentage)",
    )

    parser.add_argument(
        "--max-spectra-intensity",
        type=float,
        default=None,
        help="Maximum intensity threshold for MS/MS peaks (absolute or 0.0-1.0 for percentage)",
    )

    parser.add_argument(
        "--min-spectra-mz",
        type=float,
        default=None,
        help="Minimum m/z filter for MS/MS peaks",
    )

    parser.add_argument(
        "--max-spectra-mz",
        type=float,
        default=None,
        help="Maximum m/z filter for MS/MS peaks",
    )

    parser.add_argument(
        "--min-precursor-intensity",
        type=int,
        default=None,
        help="Minimum precursor intensity filter",
    )

    parser.add_argument(
        "--max-precursor-intensity",
        type=int,
        default=None,
        help="Maximum precursor intensity filter",
    )

    parser.add_argument(
        "--min-precursor-charge",
        type=int,
        default=None,
        help="Minimum precursor charge state filter",
    )

    parser.add_argument(
        "--max-precursor-charge",
        type=int,
        default=None,
        help="Maximum precursor charge state filter",
    )

    parser.add_argument(
        "--min-precursor-mz",
        type=float,
        default=None,
        help="Minimum precursor m/z filter",
    )

    parser.add_argument(
        "--max-precursor-mz",
        type=float,
        default=None,
        help="Maximum precursor m/z filter",
    )

    parser.add_argument(
        "--min-precursor-rt",
        type=float,
        default=None,
        help="Minimum precursor retention time filter (seconds)",
    )

    parser.add_argument(
        "--max-precursor-rt",
        type=float,
        default=None,
        help="Maximum precursor retention time filter (seconds)",
    )

    parser.add_argument(
        "--min-precursor-ccs",
        type=float,
        default=None,
        help="Minimum precursor CCS filter",
    )

    parser.add_argument(
        "--max-precursor-ccs",
        type=float,
        default=None,
        help="Maximum precursor CCS filter",
    )

    parser.add_argument(
        "--min-precursor-neutral-mass",
        type=float,
        default=None,
        help="Minimum precursor neutral mass filter",
    )

    parser.add_argument(
        "--max-precursor-neutral-mass",
        type=float,
        default=None,
        help="Maximum precursor neutral mass filter",
    )

    parser.add_argument(
        "--mz_precision",
        type=int,
        default=5,
        help="Number of decimal places for m/z values (default: 4)",
    )

    parser.add_argument(
        "--intensity_precision",
        type=int,
        default=0,
        help="Number of decimal places for intensity values (default: 0)",
    )

    parser.add_argument(
        "--keep-empty-spectra",
        action="store_true",
        help="Keep spectra with no peaks (default: False)",
    )

    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Enable verbose logging"
    )

    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing output file if it exists",
    )


def add_ms2_specific_args(parser: argparse.ArgumentParser) -> None:
    """Add MS2-specific arguments."""

    parser.add_argument(
        "--ip2",
        action="store_true",
        help="Use IP2 Settings",
    )


def add_mgf_specific_args(parser: argparse.ArgumentParser) -> None:
    """Add MGF-specific arguments."""

    parser.add_argument(
        "--casanovo",
        action="store_true",
        help="Use Casanovo preset settings (optimized for de novo sequencing)",
    )


def create_ms2_parser() -> argparse.ArgumentParser:
    """Create argument parser for MS2 extractor."""

    parser = argparse.ArgumentParser(
        description="Extract MS2 files from TimsTOF .D folders",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    add_common_args(parser)
    add_ms2_specific_args(parser)

    return parser


def create_mgf_parser() -> argparse.ArgumentParser:
    """Create argument parser for MGF extractor."""

    parser = argparse.ArgumentParser(
        description="Extract MGF files from TimsTOF .D folders",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    add_common_args(parser)
    add_mgf_specific_args(parser)

    return parser


def apply_preset_settings(args: argparse.Namespace) -> None:
    """Apply preset settings based on flags."""

    if hasattr(args, "ip2") and args.ip2:
        # Set IP2 specific defaults
        args.min_precursor_charge = 2
        args.top_n_peaks = 500

    if hasattr(args, "casanovo") and args.casanovo:
        # Set Casanovo specific defaults
        args.remove_precursor = True
        args.precursor_peak_width = 2.0
        args.top_n_peaks = 150
        args.min_spectra_intensity = 0.01
        args.min_spectra_mz = 50
        args.max_spectra_mz = 2500
        args.min_precursor_charge = 2


def log_common_args(logger, args: argparse.Namespace, extractor_type: str) -> None:
    """Log common arguments for both extractors."""

    logger.info(f"{extractor_type} Extractor Arguments:")
    logger.info(f"  Analysis Directory: {args.analysis_dir}")
    logger.info(f"  Output File: {args.output if args.output else 'Within .d folder'}")
    logger.info(f"  Remove Precursor: {args.remove_precursor}")
    logger.info(f"  Precursor Peak Width: {args.precursor_peak_width} Da")
    logger.info(f"  Batch Size: {args.batch_size}")
    logger.info(
        f"  Top N Peaks: {args.top_n_peaks if args.top_n_peaks is not None else 'All'}"
    )
    logger.info(
        f"  Min Spectra Intensity: {args.min_spectra_intensity if args.min_spectra_intensity is not None else 'None'}"
    )
    logger.info(
        f"  Max Spectra Intensity: {args.max_spectra_intensity if args.max_spectra_intensity is not None else 'None'}"
    )
    logger.info(
        f"  Min Spectra m/z: {args.min_spectra_mz if args.min_spectra_mz is not None else 'None'}"
    )
    logger.info(
        f"  Max Spectra m/z: {args.max_spectra_mz if args.max_spectra_mz is not None else 'None'}"
    )
    logger.info(
        f"  Min Precursor Intensity: {args.min_precursor_intensity if args.min_precursor_intensity is not None else 'None'}"
    )
    logger.info(
        f"  Max Precursor Intensity: {args.max_precursor_intensity if args.max_precursor_intensity is not None else 'None'}"
    )
    logger.info(
        f"  Min Precursor Charge: {args.min_precursor_charge if args.min_precursor_charge is not None else 'None'}"
    )
    logger.info(
        f"  Max Precursor Charge: {args.max_precursor_charge if args.max_precursor_charge is not None else 'None'}"
    )
    logger.info(
        f"  Min Precursor m/z: {args.min_precursor_mz if args.min_precursor_mz is not None else 'None'}"
    )
    logger.info(
        f"  Max Precursor m/z: {args.max_precursor_mz if args.max_precursor_mz is not None else 'None'}"
    )
    logger.info(
        f"  Min Precursor RT: {args.min_precursor_rt if args.min_precursor_rt is not None else 'None'} seconds"
    )
    logger.info(
        f"  Max Precursor RT: {args.max_precursor_rt if args.max_precursor_rt is not None else 'None'} seconds"
    )
    logger.info(
        f"  Min Precursor CCS: {args.min_precursor_ccs if args.min_precursor_ccs is not None else 'None'}"
    )
    logger.info(
        f"  Max Precursor CCS: {args.max_precursor_ccs if args.max_precursor_ccs is not None else 'None'}"
    )
    logger.info(
        f"  Min Precursor Neutral Mass: {args.min_precursor_neutral_mass if args.min_precursor_neutral_mass is not None else 'None'}"
    )
    logger.info(
        f"  Max Precursor Neutral Mass: {args.max_precursor_neutral_mass if args.max_precursor_neutral_mass is not None else 'None'}"
    )
    logger.info(f"  Overwrite Existing Output: {args.overwrite}")
    logger.info(f"  m/z Precision: {args.mz_precision} decimal places")
    logger.info(f"  Intensity Precision: {args.intensity_precision} decimal places")
    logger.info(f"  Keep Empty Spectra: {args.keep_empty_spectra}")
    logger.info(f"  Verbose Logging: {args.verbose}")

    # Log preset-specific settings
    if hasattr(args, "ip2") and args.ip2:
        logger.info(f"  IP2 Settings: {args.ip2}")

    if hasattr(args, "casanovo") and args.casanovo:
        logger.info(f"  Casanovo Preset: {args.casanovo}")
