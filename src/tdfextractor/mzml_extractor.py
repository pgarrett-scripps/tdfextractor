"""
mzml_extractor: write mzML files from Bruker timsTOF .D folders using psims.

Currently supports DDA acquisitions. MS1 spectra are pulled directly from the
high-level tdfpy DDA reader; MS2 PASEF spectra are produced via the existing
filtering pipeline used by ms2/mgf extractors and then written through
psims.mzml.writer.MzMLWriter.
"""

import logging
import os
import time
from pathlib import Path
from typing import Dict, Optional

import numpy as np
from psims.mzml.writer import MzMLWriter
from tdfpy import DDA, PandasTdf
from tqdm import tqdm

from .cli_args import apply_preset_settings, create_mzml_parser, log_common_args
from .utils import (
    get_ms1_frames_ids,
    get_ms2_dda_content,
    get_tdf_df,
    map_frame_id_to_ms1_scan,
    map_parent_id_to_precursors,
)

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def _scan_id(index: int) -> str:
    return f"scan={index}"


def _iter_ms1_frames(analysis_dir: str, frame_ids):
    """Yield (frame_id, mz_array, intensity_array, rt_seconds) per MS1 frame."""
    with DDA(analysis_dir) as dda:
        for fid in frame_ids:
            frame = dda.ms1.get(fid)
            if frame is None:
                continue
            peaks = frame.centroid()  # (N, 3): mz, intensity, mobility
            if peaks is None or peaks.size == 0:
                mz_arr = np.empty(0, dtype=np.float64)
                int_arr = np.empty(0, dtype=np.float32)
            else:
                mz_arr = np.asarray(peaks[:, 0], dtype=np.float64)
                int_arr = np.asarray(peaks[:, 1], dtype=np.float32)
            yield fid, mz_arr, int_arr, float(frame.time)


def write_mzml_file(
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
    min_precursor_intensity: Optional[float] = None,
    max_precursor_intensity: Optional[float] = None,
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
    keep_empty_spectra: bool = False,
    include_ms1: bool = True,
):
    """Write an indexed mzML file containing MS1 and MS2 (PASEF) spectra."""
    from . import __version__ as _ext_version

    start_time = time.time()

    if output_file is None:
        output_file = str(Path(analysis_dir) / Path(analysis_dir).stem) + ".mzML"

    logger.info("Loading TDF metadata")
    pd_tdf = PandasTdf(str(Path(analysis_dir) / "analysis.tdf"))
    if not pd_tdf.is_dda:
        raise TypeError(
            "mzml extraction currently only supports DDA acquisitions; "
            f"got file is_dda={pd_tdf.is_dda} is_prm={pd_tdf.is_prm}"
        )

    frames_df = pd_tdf.frames
    precursors_df = pd_tdf.precursors

    # Build deterministic scan numbering shared across MS1 and MS2 (matches IP2 ordering).
    ms1_frame_ids = [int(f) for f in get_ms1_frames_ids(frames_df).tolist()]
    parent_to_precs = map_parent_id_to_precursors(precursors_df)
    frame_id_to_ms1_scan, ms2_scan_map = map_frame_id_to_ms1_scan(
        parent_to_precs, ms1_frame_ids
    )

    # Filter precursors via the shared pipeline
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
        max_precursor_neutral_mass,
    )

    # Pre-stage MS2 spectra grouped by parent frame so we can interleave with MS1.
    logger.info("Extracting MS2 spectra")
    ms2_by_parent: Dict[int, list] = {}
    for spectrum in tqdm(
        get_ms2_dda_content(
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
        ),
        total=len(merged_df),
        desc="Reading MS2",
    ):
        if (not keep_empty_spectra) and len(spectrum.mz_spectra) == 0:
            continue
        ms2_by_parent.setdefault(int(spectrum.parent_id), []).append(spectrum)

    total_ms2 = sum(len(v) for v in ms2_by_parent.values())
    total_spectra = total_ms2 + (len(ms1_frame_ids) if include_ms1 else 0)
    logger.info(
        f"Writing mzML to {output_file} ({total_spectra} spectra: "
        f"{len(ms1_frame_ids) if include_ms1 else 0} MS1, {total_ms2} MS2)"
    )

    with MzMLWriter(open(output_file, "wb"), close=True) as writer:
        writer.controlled_vocabularies()
        writer.file_description(
            file_contents=[
                "MS1 spectrum",
                "MSn spectrum",
                "centroid spectrum",
            ],
            source_files=[
                {
                    "id": "RAW1",
                    "name": Path(analysis_dir).name,
                    "location": f"file:///{Path(analysis_dir).resolve().as_posix()}",
                    "params": ["Bruker TDF format"],
                }
            ],
        )
        writer.software_list(
            [
                {
                    "id": "tdfextractor",
                    "version": _ext_version,
                    "params": ["python-psims"],
                }
            ]
        )
        writer.data_processing_list(
            [
                writer.DataProcessing(
                    id="DP1",
                    processing_methods=[
                        writer.ProcessingMethod(
                            order=1,
                            software_reference="tdfextractor",
                            params=["Conversion to mzML"],
                        )
                    ],
                )
            ]
        )

        with writer.run(id=Path(analysis_dir).stem):
            with writer.spectrum_list(count=total_spectra):
                if include_ms1:
                    ms1_iter = _iter_ms1_frames(analysis_dir, ms1_frame_ids)
                else:
                    ms1_iter = iter(())

                pbar = tqdm(total=total_spectra, desc="Writing mzML", unit="spectra")
                for frame_id, mz_arr, int_arr, rt_s in ms1_iter:
                    ms1_scan_index = frame_id_to_ms1_scan.get(frame_id)
                    if ms1_scan_index is None:
                        continue
                    ms1_id = _scan_id(ms1_scan_index)
                    writer.write_spectrum(
                        mz_arr,
                        int_arr,
                        id=ms1_id,
                        centroided=True,
                        scan_start_time=rt_s / 60.0,  # mzML scan time in minutes
                        params=[
                            "MS1 spectrum",
                            {"ms level": 1},
                            {"total ion current": float(np.sum(int_arr))},
                        ],
                    )
                    pbar.update(1)

                    for ms2 in ms2_by_parent.get(frame_id, []):
                        ms2_scan_index = ms2_scan_map.get(frame_id, {}).get(
                            int(ms2.precursor_id)
                        )
                        if ms2_scan_index is None:
                            continue
                        iso_mz = float(getattr(ms2, "iso_mz", ms2.mz))
                        iso_w = float(getattr(ms2, "iso_width", 2.0))
                        ce = float(getattr(ms2, "ce", 0.0))
                        ook0 = float(getattr(ms2, "ook0", 0.0))
                        mz2 = np.asarray(ms2.mz_spectra, dtype=np.float64)
                        int2 = np.asarray(ms2.intensity_spectra, dtype=np.float32)
                        writer.write_spectrum(
                            mz2,
                            int2,
                            id=_scan_id(ms2_scan_index),
                            centroided=True,
                            scan_start_time=float(ms2.rt) / 60.0,
                            params=[
                                "MSn spectrum",
                                {"ms level": 2},
                                {"total ion current": float(int2.sum())},
                            ],
                            scan_params=[
                                {"inverse reduced ion mobility": ook0},
                            ],
                            precursor_information={
                                "mz": float(ms2.mz),
                                "intensity": float(ms2.prec_intensity),
                                "charge": int(ms2.charge),
                                "scan_id": ms1_id,
                                "activation": [
                                    "beam-type collisional dissociation",
                                    {"collision energy": ce},
                                ],
                                "isolation_window": [iso_w / 2.0, iso_mz, iso_w / 2.0],
                            },
                        )
                        pbar.update(1)
                pbar.close()

    total_time = round(time.time() - start_time, 2)
    logger.info(f"mzML extraction complete in {total_time:.2f} seconds")


def main():
    """Command-line interface for mzML extraction from TimsTOF data."""

    parser = create_mzml_parser()
    args = parser.parse_args()

    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=log_level, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )

    apply_preset_settings(logger, args)
    log_common_args(logger, args, "mzML")

    analysis_path = Path(args.analysis_dir)
    if not analysis_path.exists():
        logger.error(f"Analysis directory does not exist: {args.analysis_dir}")
        return 1

    if not analysis_path.is_dir():
        logger.error(f"Path is not a directory: {args.analysis_dir}")
        return 1

    if analysis_path.name.endswith(".d"):
        d_folders = [analysis_path]
        logger.info(f"Using provided .d folder: {analysis_path}")
    else:
        d_folders = list(analysis_path.glob("*.d"))
        logger.info(f"Found {len(d_folders)} .d folders in: {args.analysis_dir}")

    if not d_folders:
        logger.error(f"No .d folders found in: {args.analysis_dir}")
        return 1

    output_dir = None
    output_name = None

    if args.output is None:
        output_dir = None
        output_name = None
    elif args.output.endswith(".mzML") or args.output.endswith(".mzml"):
        if len(d_folders) > 1:
            raise ValueError("Output file specified but multiple .d folders found.")
        output_dir = Path(args.output).parent
        output_name = Path(args.output).name
    else:
        output_dir = Path(args.output)
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
            output_name if output_name is not None else Path(d_folder).stem + ".mzML"
        )
        output = os.path.join(_output_dir, _output_name)
        logger.info(f"Output file: {output}")

        if not args.overwrite and Path(output).exists():
            logger.warning(f"Output file {output} already exists. Skipping...")
            continue

        try:
            write_mzml_file(
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
                include_ms1=not args.no_ms1,
            )
            logger.info("mzML extraction completed successfully!")
        except Exception as e:
            logger.error(f"Error during mzML extraction: {e}... skipping {d_folder}")
            logger.error(e, exc_info=True)
        except KeyboardInterrupt:
            logger.info("Extraction interrupted by user.")
            os._exit(0)


if __name__ == "__main__":
    exit(main())
