"""
ms2_extractor defines functions for generating ms2 files from DDA and PRM based .D folders
"""

import argparse
import logging
import os
import queue
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from tqdm import tqdm

from .args import MgfArgs
from .cli_args import apply_preset_settings, create_mgf_parser, log_common_args
from .utils import get_ms2_dda_content, get_tdf_df

logger = logging.getLogger(__name__)


def write_mgf_file(args: MgfArgs) -> None:

    start_time = time.time()

    analysis_dir = args.analysis_dir
    output_file = args.output_file or (str(Path(analysis_dir) / Path(analysis_dir).stem) + ".mgf")

    spectra_queue = queue.Queue(maxsize=100)

    merged_df = get_tdf_df(
        analysis_dir,
        args.min_precursor_intensity,
        args.max_precursor_intensity,
        args.min_precursor_charge,
        args.max_precursor_charge,
        args.min_precursor_mz,
        args.max_precursor_mz,
        args.min_precursor_rt,
        args.max_precursor_rt,
        args.min_precursor_ccs,
        args.max_precursor_ccs,
        args.min_precursor_neutral_mass,
        args.max_precursor_neutral_mass,
    )

    def producer():
        try:
            ms2_spectra = get_ms2_dda_content(
                analysis_dir=analysis_dir,
                merged_df=merged_df,
                remove_precursor=args.remove_precursor,
                precursor_peak_width=args.precursor_peak_width,
                batch_size=args.batch_size,
                top_n_peaks=args.top_n_peaks,
                min_spectra_intensity=args.min_spectra_intensity,
                max_spectra_intensity=args.max_spectra_intensity,
                min_spectra_mz=args.min_spectra_mz,
                max_spectra_mz=args.max_spectra_mz,
            )
            for spectrum in ms2_spectra:
                spectra_queue.put(spectrum)
        finally:
            spectra_queue.put(None)  # Sentinel value

    def consumer():
        logger.info("Writing Contents To File")
        with open(output_file, "w", encoding="UTF-8") as file:
            with tqdm(desc="Writing MGF File", unit="spectra", total=len(merged_df)) as pbar:
                # https://www.matrixscience.com/help/data_file_help.html
                header_lines = []
                header_lines.append("INSTRUMENT=TimsTOF")
                header_lines.append("MASS=Mono")

                while True:
                    spectrum = spectra_queue.get()
                    if spectrum is None:
                        break

                    pbar.update(1)

                    if len(spectrum.mz_spectra) == 0 and args.keep_empty_spectra is False:
                        continue

                    mgf_lines = []
                    mgf_lines.append("BEGIN IONS")
                    mgf_lines.append(
                        f"TITLE={Path(analysis_dir).stem}.{spectrum.low_scan}.{spectrum.high_scan}.{spectrum.charge} "
                        f'File="{Path(analysis_dir).stem}", NativeID="merged={spectrum.precursor_id} frame={spectrum.parent_id} '
                        f'scanStart={spectrum.scan_begin} scanEnd={spectrum.scan_end} scan={spectrum.low_scan}"'
                    )
                    mgf_lines.append(f"RTINSECONDS={spectrum.rt:.2f}")
                    # Pepmass is actually mz? huh?
                    mgf_lines.append(
                        f"PEPMASS={spectrum.mz:.6f} {spectrum.prec_intensity:.{args.intensity_precision}f}"
                    )
                    mgf_lines.append(f"CHARGE={spectrum.charge}+")
                    for mz, intensity in zip(spectrum.mz_spectra, spectrum.intensity_spectra):
                        mgf_lines.append(
                            f"{mz:.{args.mz_precision}f} {intensity:.{args.intensity_precision}f}"
                        )
                    mgf_lines.append("END IONS")
                    file.write("\n".join(mgf_lines) + "\n\n")

    producer_thread = threading.Thread(target=producer)
    consumer_thread = threading.Thread(target=consumer)

    producer_thread.start()
    consumer_thread.start()
    producer_thread.join()
    consumer_thread.join()

    total_time = round(time.time() - start_time, 2)
    logger.info(f"Total Time: {total_time:.2f} seconds")


def process_single_d_folder(
    d_folder: Path,
    cli_ns: argparse.Namespace,
    output_dir: Path | None,
    output_name: str | None,
) -> bool:
    """Process a single .d folder with error handling.

    ``cli_ns`` is the raw argparse Namespace; a fresh :class:`MgfArgs` is
    built per call so worker threads never share mutable state.
    """
    try:
        if not d_folder.is_dir():
            logger.error(f"Path is not a directory: {d_folder}")
            return False

        if not (d_folder / "analysis.tdf").exists():
            logger.error(f"Required file not found in {d_folder}: analysis.tdf")
            return False
        if not (d_folder / "analysis.tdf_bin").exists():
            logger.error(f"Required file not found in {d_folder}: analysis.tdf_bin")
            return False

        logger.info(f"Processing {d_folder}...")

        _output_dir = output_dir if output_dir is not None else d_folder
        _output_name = output_name if output_name is not None else Path(d_folder).stem + ".mgf"

        output = os.path.join(_output_dir, _output_name)
        logger.info(f"Output file: {output}")

        if not cli_ns.overwrite and Path(output).exists():
            logger.warning(f"Output file {output} already exists. Skipping...")
            return True

        mgf_args = MgfArgs.from_namespace(cli_ns)
        mgf_args.analysis_dir = str(d_folder)
        mgf_args.output_file = output
        write_mgf_file(mgf_args)
        logger.info(f"MGF extraction completed successfully for {d_folder}!")
        return True
    except Exception as e:
        logger.error(f"Error during MGF extraction for {d_folder}: {e}")
        return False


def main() -> int | None:
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
    apply_preset_settings(logger, args)

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

    # Process .d folders with multiple workers if specified
    if len(d_folders) > 1 and args.workers > 1:
        logger.info(f"Processing {len(d_folders)} .d folders using {args.workers} workers...")

        successful_count = 0
        failed_count = 0

        try:
            with ThreadPoolExecutor(max_workers=args.workers) as executor:
                # Submit all jobs
                future_to_folder = {
                    executor.submit(
                        process_single_d_folder, d_folder, args, output_dir, output_name
                    ): d_folder
                    for d_folder in d_folders
                }

                # Process completed jobs
                for future in as_completed(future_to_folder):
                    d_folder = future_to_folder[future]
                    try:
                        success = future.result()
                        if success:
                            successful_count += 1
                        else:
                            failed_count += 1
                    except Exception as e:
                        logger.error(f"Unexpected error processing {d_folder}: {e}")
                        failed_count += 1
        except KeyboardInterrupt:
            logger.info("\nExtraction interrupted by user.")
            os._exit(0)

        logger.info(f"Processing completed: {successful_count} successful, {failed_count} failed")
    else:
        # Process sequentially (original behavior)
        for d_folder in d_folders:
            try:
                success = process_single_d_folder(d_folder, args, output_dir, output_name)
                if not success:
                    continue
            except KeyboardInterrupt:
                logger.info("\nExtraction interrupted by user.")
                os._exit(0)


if __name__ == "__main__":
    exit(main())
