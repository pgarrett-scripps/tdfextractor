"""
mzml_extractor: write mzML files from Bruker timsTOF .D folders using psims.

Supports DDA, DIA, and PRM acquisitions. The acquisition type is detected from
the TDF metadata and the appropriate writer routine is dispatched.

For all acquisition types, MS1 spectra are written through the high-level
tdfpy reader and include a per-peak ``mean inverse reduced ion mobility``
array (the third column of the centroided MS1 peaks). MS2 spectra are written
with their isolation window, collision energy, and per-precursor inverse
reduced ion mobility metadata.

Per-array compression and encoding can be configured via CLI flags or via
keyword arguments to :func:`write_mzml_file`.
"""

import logging
import os
import time
from collections import defaultdict
from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import Any

import numpy as np
from psims.mzml.writer import MzMLWriter
from tdfpy import DDA, DIA, PRM, PandasTdf
from tqdm import tqdm

from .args import EncodingBitWidth, MzmlArgs
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


# psims array name constants
_MZ_ARRAY = "m/z array"
_INTENSITY_ARRAY = "intensity array"
_MOBILITY_ARRAY = "mean inverse reduced ion mobility array"


# CLI compression name -> psims compression identifier (see
# psims.mzml.binary_encoding.compressors). Numpress / zstd entries are only
# usable when their backing libraries are installed.
_COMPRESSION_NAME_MAP: dict[str, str] = {
    "none": "none",
    "zlib": "zlib",
    "zstd": "zstd",
    "numpress-linear": "MS-Numpress linear prediction compression",
    "numpress-slof": "MS-Numpress short logged float compression",
    "numpress-pic": "MS-Numpress positive integer compression",
}


def _resolve_compression(name: str | None) -> str:
    """Translate a CLI compression name into the psims identifier."""

    if name is None:
        return "zlib"
    try:
        return _COMPRESSION_NAME_MAP[name]
    except KeyError as exc:
        raise ValueError(
            f"Unknown compression {name!r}; expected one of {sorted(_COMPRESSION_NAME_MAP)}"
        ) from exc


def _resolve_encoding(bits: EncodingBitWidth | None) -> Any:
    """Translate a 32/64-bit width to the corresponding numpy dtype."""

    if bits is None or bits == 64:
        return np.float64
    if bits == 32:
        return np.float32
    # Defensive: callers using the type checker can never reach this
    # branch, but the runtime check protects untyped callers.
    raise ValueError(f"Unsupported encoding bit width: {bits!r}")


def _build_compression_dict(
    mz_compression: str,
    intensity_compression: str,
    mobility_compression: str,
) -> dict[str, str]:
    return {
        _MZ_ARRAY: _resolve_compression(mz_compression),
        _INTENSITY_ARRAY: _resolve_compression(intensity_compression),
        _MOBILITY_ARRAY: _resolve_compression(mobility_compression),
    }


def _build_encoding_dict(
    mz_encoding: EncodingBitWidth,
    intensity_encoding: EncodingBitWidth,
) -> dict[str, Any]:
    return {
        _MZ_ARRAY: _resolve_encoding(mz_encoding),
        _INTENSITY_ARRAY: _resolve_encoding(intensity_encoding),
        _MOBILITY_ARRAY: np.float64,
    }


def _build_centroid_kwargs(args: MzmlArgs) -> dict[str, Any]:
    """Build a kwargs dict for tdfpy's centroid() from MzmlArgs centroid fields."""
    return {
        "mz_tolerance": args.centroid_mz_tolerance,
        "mz_tolerance_type": args.centroid_mz_tolerance_type,
        "im_tolerance": args.centroid_im_tolerance,
        "im_tolerance_type": args.centroid_im_tolerance_type,
        "min_peaks": args.centroid_min_peaks,
        "noise_filter": args.centroid_noise_filter,
    }


def _scan_id(index: int) -> str:
    return f"scan={index}"


def _split_centroided_peaks(
    peaks: np.ndarray | None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray | None]:
    """Split a centroided peak array into (mz, intensity, optional mobility)."""

    if peaks is None or peaks.size == 0:
        return (
            np.empty(0, dtype=np.float64),
            np.empty(0, dtype=np.float32),
            np.empty(0, dtype=np.float64),
        )
    arr = np.asarray(peaks)
    mz = np.asarray(arr[:, 0], dtype=np.float64)
    intensity = np.asarray(arr[:, 1], dtype=np.float32)
    if arr.shape[1] >= 3:
        mobility = np.asarray(arr[:, 2], dtype=np.float64)
    else:
        mobility = np.empty(0, dtype=np.float64)
    return mz, intensity, mobility


def _write_ms1_spectrum(
    writer: MzMLWriter,
    *,
    scan_id: str,
    mz: np.ndarray,
    intensity: np.ndarray,
    mobility: np.ndarray | None,
    rt_seconds: float,
    compression: Mapping[str, str],
    encoding: Mapping[str, Any],
) -> None:
    other_arrays: list[tuple[Any, np.ndarray]] = []
    if mobility is not None and mobility.size == mz.size and mobility.size > 0:
        other_arrays.append((_MOBILITY_ARRAY, mobility))

    writer.write_spectrum(
        mz,
        intensity,
        id=scan_id,
        centroided=True,
        scan_start_time=rt_seconds / 60.0,
        params=[
            "MS1 spectrum",
            {"ms level": 1},
            {"total ion current": float(np.sum(intensity))},
        ],
        other_arrays=other_arrays,
        compression=dict(compression),
        encoding=dict(encoding),
    )


def _write_ms2_spectrum(
    writer: MzMLWriter,
    *,
    scan_id: str,
    parent_scan_id: str | None,
    mz: np.ndarray,
    intensity: np.ndarray,
    rt_seconds: float,
    iso_mz: float,
    iso_width: float,
    collision_energy: float,
    inverse_reduced_ion_mobility: float,
    precursor_mz: float,
    precursor_intensity: float | None,
    precursor_charge: int | None,
    compression: Mapping[str, str],
    encoding: Mapping[str, Any],
) -> None:
    half_width = iso_width / 2.0
    precursor_info: dict[str, Any] = {
        "mz": float(precursor_mz),
        "activation": [
            "beam-type collisional dissociation",
            {"collision energy": float(collision_energy)},
        ],
        "isolation_window": [half_width, float(iso_mz), half_width],
    }
    if parent_scan_id is not None:
        precursor_info["scan_id"] = parent_scan_id
    if precursor_intensity is not None:
        precursor_info["intensity"] = float(precursor_intensity)
    if precursor_charge is not None:
        precursor_info["charge"] = int(precursor_charge)

    writer.write_spectrum(
        mz,
        intensity,
        id=scan_id,
        centroided=True,
        scan_start_time=rt_seconds / 60.0,
        params=[
            "MSn spectrum",
            {"ms level": 2},
            {"total ion current": float(np.sum(intensity))},
        ],
        scan_params=[
            {"inverse reduced ion mobility": float(inverse_reduced_ion_mobility)},
        ],
        precursor_information=precursor_info,
        compression=dict(compression),
        encoding=dict(encoding),
    )


# ---------------------------------------------------------------------------
# Header / file_description / software / data_processing
# ---------------------------------------------------------------------------


def _write_header(writer: MzMLWriter, analysis_dir: str) -> None:
    from . import __version__ as _ext_version

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
    writer.instrument_configuration_list(
        [
            writer.InstrumentConfiguration(
                "IC1",
                writer.ComponentList(
                    [
                        writer.Source(1, ["electrospray ionization"]),
                        writer.Analyzer(2, ["quadrupole", "time-of-flight"]),
                        writer.Detector(3, ["microchannel plate detector"]),
                    ]
                ),
                ["instrument model"],
            )
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


# ---------------------------------------------------------------------------
# DDA writer
# ---------------------------------------------------------------------------


def _iter_dda_ms1(
    analysis_dir: str,
    frame_ids: Iterable[int],
    centroid_kwargs: dict,
) -> Iterable[tuple[int, np.ndarray, np.ndarray, np.ndarray, float]]:
    """Yield (frame_id, mz, intensity, mobility, rt_seconds) per DDA MS1 frame."""

    with DDA(analysis_dir) as dda:
        for fid in frame_ids:
            frame = dda.ms1[int(fid)]
            if frame is None:
                continue
            mz, intensity, mobility = _split_centroided_peaks(frame.centroid(**centroid_kwargs))
            yield int(fid), mz, intensity, mobility, float(frame.time)


def _write_dda(
    *,
    writer: MzMLWriter,
    pd_tdf: PandasTdf,
    args: MzmlArgs,
    compression: Mapping[str, str],
    encoding: Mapping[str, Any],
) -> None:
    analysis_dir = args.analysis_dir
    frames_df = pd_tdf.frames
    precursors_df = pd_tdf.precursors

    # Honor --min/--max-precursor-rt for MS1 frames as well as MS2 precursors,
    # so the writer produces a coherent RT-bounded slice of the file rather
    # than dropping MS2 while still emitting every MS1 frame in the run.
    ms1_frames_view = get_ms1_frames_ids(frames_df)
    ms1_frame_ids: list[int] = []
    ms1_times = frames_df.set_index("Id")["Time"]
    for fid in ms1_frames_view.tolist():
        rt = float(ms1_times.loc[int(fid)])
        if args.min_precursor_rt is not None and rt < args.min_precursor_rt:
            continue
        if args.max_precursor_rt is not None and rt > args.max_precursor_rt:
            continue
        ms1_frame_ids.append(int(fid))
    parent_to_precs = map_parent_id_to_precursors(precursors_df)
    frame_id_to_ms1_scan, ms2_scan_map = map_frame_id_to_ms1_scan(parent_to_precs, ms1_frame_ids)

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

    logger.info("Extracting MS2 spectra")
    ms2_by_parent: dict[int, list] = {}
    for spectrum in tqdm(
        get_ms2_dda_content(
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
        ),
        total=len(merged_df),
        desc="Reading MS2",
    ):
        if (not args.keep_empty_spectra) and len(spectrum.mz_spectra) == 0:
            continue
        ms2_by_parent.setdefault(int(spectrum.parent_id), []).append(spectrum)

    total_ms2 = sum(len(v) for v in ms2_by_parent.values())
    total_spectra = total_ms2 + (len(ms1_frame_ids) if args.include_ms1 else 0)
    logger.info(
        f"Writing mzML ({total_spectra} spectra: "
        f"{len(ms1_frame_ids) if args.include_ms1 else 0} MS1, {total_ms2} MS2)"
    )

    centroid_kwargs = _build_centroid_kwargs(args)

    with writer.run(id=Path(analysis_dir).stem):
        with writer.spectrum_list(count=total_spectra):
            ms1_iter = (
                _iter_dda_ms1(analysis_dir, ms1_frame_ids, centroid_kwargs)
                if args.include_ms1
                else iter(())
            )

            pbar = tqdm(total=total_spectra, desc="Writing mzML", unit="spectra")
            for frame_id, mz_arr, int_arr, mob_arr, rt_s in ms1_iter:
                ms1_scan_index = frame_id_to_ms1_scan.get(frame_id)
                if ms1_scan_index is None:
                    continue
                ms1_id = _scan_id(ms1_scan_index)
                _write_ms1_spectrum(
                    writer,
                    scan_id=ms1_id,
                    mz=mz_arr,
                    intensity=int_arr,
                    mobility=mob_arr,
                    rt_seconds=rt_s,
                    compression=compression,
                    encoding=encoding,
                )
                pbar.update(1)

                for ms2 in ms2_by_parent.get(frame_id, []):
                    ms2_scan_index = ms2_scan_map.get(frame_id, {}).get(int(ms2.precursor_id))
                    if ms2_scan_index is None:
                        continue
                    iso_mz = float(getattr(ms2, "iso_mz", ms2.mz))
                    iso_w = float(getattr(ms2, "iso_width", 2.0))
                    ce = float(getattr(ms2, "ce", 0.0))
                    ook0 = float(getattr(ms2, "ook0", 0.0))
                    mz2 = np.asarray(ms2.mz_spectra, dtype=np.float64)
                    int2 = np.asarray(ms2.intensity_spectra, dtype=np.float32)
                    _write_ms2_spectrum(
                        writer,
                        scan_id=_scan_id(ms2_scan_index),
                        parent_scan_id=ms1_id,
                        mz=mz2,
                        intensity=int2,
                        rt_seconds=float(ms2.rt),
                        iso_mz=iso_mz,
                        iso_width=iso_w,
                        collision_energy=ce,
                        inverse_reduced_ion_mobility=ook0,
                        precursor_mz=float(ms2.mz),
                        precursor_intensity=float(ms2.prec_intensity),
                        precursor_charge=int(ms2.charge),
                        compression=compression,
                        encoding=encoding,
                    )
                    pbar.update(1)
            pbar.close()


# ---------------------------------------------------------------------------
# DIA / PRM writer (shared shape)
# ---------------------------------------------------------------------------


def _collect_windowed_ms2(
    windows_iter: Iterable[Any],
) -> tuple[dict[int, list], int]:
    """Group DIA windows or PRM transitions by their parent frame_id."""

    grouped: dict[int, list] = defaultdict(list)
    total = 0
    for w in windows_iter:
        grouped[int(w.frame_id)].append(w)
        total += 1
    return grouped, total


def _write_dia_or_prm(
    *,
    writer: MzMLWriter,
    pd_tdf: PandasTdf,
    reader_factory,
    args: MzmlArgs,
    compression: Mapping[str, str],
    encoding: Mapping[str, Any],
) -> None:
    analysis_dir = args.analysis_dir
    include_ms1 = args.include_ms1
    keep_empty_spectra = args.keep_empty_spectra
    min_precursor_mz = args.min_precursor_mz
    max_precursor_mz = args.max_precursor_mz
    min_precursor_rt = args.min_precursor_rt
    max_precursor_rt = args.max_precursor_rt

    frames_df = pd_tdf.frames.sort_values("Id").reset_index(drop=True)
    # Honor --min/--max-precursor-rt for the whole frames table so that MS1
    # frames outside the RT window are skipped along with the MS2 windows
    # the loop below already filters explicitly.
    if min_precursor_rt is not None:
        frames_df = frames_df[frames_df["Time"] >= min_precursor_rt]
    if max_precursor_rt is not None:
        frames_df = frames_df[frames_df["Time"] <= max_precursor_rt]
    frames_df = frames_df.reset_index(drop=True)

    with reader_factory(analysis_dir) as reader:
        # Materialize the windows/transitions once and group by parent frame.
        if hasattr(reader, "windows"):
            window_iter = reader.windows
            kind = "DIA"
        else:
            window_iter = reader.transitions
            kind = "PRM"
        logger.info(f"Indexing {kind} MS2 windows")
        grouped, total_ms2 = _collect_windowed_ms2(window_iter)

        frame_ids = frames_df["Id"].to_numpy(dtype=np.int64)
        msms_types = frames_df["MsMsType"].to_numpy(dtype=np.int64)
        total_ms1 = int((msms_types == 0).sum()) if include_ms1 else 0
        total_spectra = total_ms1 + total_ms2

        logger.info(
            f"Writing mzML ({total_spectra} spectra: {total_ms1} MS1, {total_ms2} {kind} MS2)"
        )

        centroid_kwargs = _build_centroid_kwargs(args)
        scan_counter = 0
        current_ms1_id: str | None = None
        pbar = tqdm(total=total_spectra, desc="Writing mzML", unit="spectra")

        with writer.run(id=Path(analysis_dir).stem):
            with writer.spectrum_list(count=total_spectra):
                for frame_id_np, msms_type_np in zip(frame_ids, msms_types):
                    frame_id = int(frame_id_np)
                    msms_type = int(msms_type_np)

                    if msms_type == 0:
                        if not include_ms1:
                            current_ms1_id = None
                            continue
                        frame = reader.ms1.get(frame_id)
                        if frame is None:
                            continue
                        mz, intensity, mobility = _split_centroided_peaks(
                            frame.centroid(**centroid_kwargs)
                        )
                        scan_counter += 1
                        ms1_id = _scan_id(scan_counter)
                        _write_ms1_spectrum(
                            writer,
                            scan_id=ms1_id,
                            mz=mz,
                            intensity=intensity,
                            mobility=mobility,
                            rt_seconds=float(frame.time),
                            compression=compression,
                            encoding=encoding,
                        )
                        current_ms1_id = ms1_id
                        pbar.update(1)
                        continue

                    # MS2 frame: write each window/transition for this frame
                    windows = grouped.get(frame_id, [])
                    for w in windows:
                        iso_mz = float(w.isolation_mz)
                        iso_w = float(w.isolation_width)
                        if min_precursor_mz is not None and iso_mz < min_precursor_mz:
                            continue
                        if max_precursor_mz is not None and iso_mz > max_precursor_mz:
                            continue
                        rt_s = float(w.rt)
                        if min_precursor_rt is not None and rt_s < min_precursor_rt:
                            continue
                        if max_precursor_rt is not None and rt_s > max_precursor_rt:
                            continue
                        peaks = w.centroid(**centroid_kwargs)
                        mz2, int2, _ = _split_centroided_peaks(peaks)
                        if (not keep_empty_spectra) and mz2.size == 0:
                            continue
                        ook0 = (float(w.ook0_begin) + float(w.ook0_end)) / 2.0
                        ce = float(w.collision_energy)
                        scan_counter += 1
                        _write_ms2_spectrum(
                            writer,
                            scan_id=_scan_id(scan_counter),
                            parent_scan_id=current_ms1_id,
                            mz=mz2,
                            intensity=int2,
                            rt_seconds=rt_s,
                            iso_mz=iso_mz,
                            iso_width=iso_w,
                            collision_energy=ce,
                            inverse_reduced_ion_mobility=ook0,
                            precursor_mz=iso_mz,
                            precursor_intensity=None,
                            precursor_charge=None,
                            compression=compression,
                            encoding=encoding,
                        )
                        pbar.update(1)
        pbar.close()


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------


def write_mzml_file(args: MzmlArgs) -> None:
    """Write an indexed mzML file from a Bruker .d folder.

    Dispatches to the DDA, DIA, or PRM writer based on the TDF metadata.
    """

    start_time = time.time()

    analysis_dir = args.analysis_dir
    output_file = args.output_file or (str(Path(analysis_dir) / Path(analysis_dir).stem) + ".mzML")

    logger.info("Loading TDF metadata")
    pd_tdf = PandasTdf(str(Path(analysis_dir) / "analysis.tdf"))

    compression = _build_compression_dict(
        args.mz_compression, args.intensity_compression, args.mobility_compression
    )
    encoding = _build_encoding_dict(args.mz_encoding, args.intensity_encoding)

    logger.info(f"Writing mzML to {output_file}")

    with MzMLWriter(open(output_file, "wb"), close=True) as writer:
        _write_header(writer, analysis_dir)

        if pd_tdf.is_dda:
            logger.info("Detected DDA acquisition")
            _write_dda(
                writer=writer,
                pd_tdf=pd_tdf,
                args=args,
                compression=compression,
                encoding=encoding,
            )
        elif pd_tdf.is_dia:
            logger.info("Detected DIA acquisition")
            _write_dia_or_prm(
                writer=writer,
                pd_tdf=pd_tdf,
                reader_factory=DIA,
                args=args,
                compression=compression,
                encoding=encoding,
            )
        elif pd_tdf.is_prm:
            logger.info("Detected PRM acquisition")
            _write_dia_or_prm(
                writer=writer,
                pd_tdf=pd_tdf,
                reader_factory=PRM,
                args=args,
                compression=compression,
                encoding=encoding,
            )
        else:
            raise TypeError(
                "mzml extraction could not determine acquisition type "
                f"(is_dda={pd_tdf.is_dda} is_dia={pd_tdf.is_dia} "
                f"is_prm={pd_tdf.is_prm})"
            )

    total_time = round(time.time() - start_time, 2)
    logger.info(f"mzML extraction complete in {total_time:.2f} seconds")


def main() -> int | None:
    """Command-line interface for mzML extraction from TimsTOF data."""

    parser = create_mzml_parser()
    args = parser.parse_args()

    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=log_level, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )

    apply_preset_settings(logger, args)
    log_common_args(logger, args, "mzML")

    base_args = MzmlArgs.from_namespace(args)

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
        _output_name = output_name if output_name is not None else Path(d_folder).stem + ".mzML"
        output = os.path.join(_output_dir, _output_name)
        logger.info(f"Output file: {output}")

        if not args.overwrite and Path(output).exists():
            logger.warning(f"Output file {output} already exists. Skipping...")
            continue

        try:
            base_args.analysis_dir = str(d_folder)
            base_args.output_file = output
            write_mzml_file(base_args)
            logger.info("mzML extraction completed successfully!")
        except Exception as e:
            logger.error(f"Error during mzML extraction: {e}... skipping {d_folder}")
            logger.error(e, exc_info=True)
        except KeyboardInterrupt:
            logger.info("Extraction interrupted by user.")
            os._exit(0)


if __name__ == "__main__":
    exit(main())
