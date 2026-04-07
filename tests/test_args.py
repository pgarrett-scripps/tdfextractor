"""Unit tests for the ExtractorArgs dataclasses.

These tests do not touch any .d folders so they run in milliseconds.
"""

from __future__ import annotations

import pytest

from tdfextractor.args import MgfArgs, Ms2Args, MzmlArgs
from tdfextractor.cli_args import (
    create_mgf_parser,
    create_ms2_parser,
    create_mzml_parser,
)


def test_ms2_args_defaults() -> None:
    a = Ms2Args(analysis_dir="/tmp/foo.d")
    assert a.analysis_dir == "/tmp/foo.d"
    assert a.output_file is None
    assert a.remove_precursor is False
    assert a.precursor_peak_width == 2.0
    assert a.batch_size == 100
    assert a.keep_empty_spectra is False
    assert a.top_n_peaks is None
    assert a.min_precursor_intensity is None
    assert a.mz_precision == 5
    assert a.intensity_precision == 0


def test_mgf_args_defaults() -> None:
    a = MgfArgs(analysis_dir="/tmp/foo.d")
    assert a.mz_precision == 5
    assert a.intensity_precision == 0


def test_mzml_args_defaults() -> None:
    a = MzmlArgs(analysis_dir="/tmp/foo.d")
    assert a.include_ms1 is True
    assert a.mz_compression == "zlib"
    assert a.intensity_compression == "zlib"
    assert a.mobility_compression == "zlib"
    assert a.mz_encoding == 64
    assert a.intensity_encoding == 32


@pytest.mark.parametrize("bad", [16, 8, 128, 0, -32])
def test_mzml_args_rejects_invalid_mz_encoding(bad: int) -> None:
    with pytest.raises(ValueError, match="mz_encoding"):
        MzmlArgs(analysis_dir="/tmp/foo.d", mz_encoding=bad)  # type: ignore[arg-type]


@pytest.mark.parametrize("bad", [16, 8, 128, 0, -32])
def test_mzml_args_rejects_invalid_intensity_encoding(bad: int) -> None:
    with pytest.raises(ValueError, match="intensity_encoding"):
        MzmlArgs(analysis_dir="/tmp/foo.d", intensity_encoding=bad)  # type: ignore[arg-type]


def test_ms2_from_namespace_basic() -> None:
    parser = create_ms2_parser()
    ns = parser.parse_args(
        ["/tmp/foo.d", "--top-n-peaks", "150", "--remove-precursor"]
    )
    args = Ms2Args.from_namespace(ns)
    assert isinstance(args, Ms2Args)
    assert args.analysis_dir == "/tmp/foo.d"
    assert args.output_file is None
    assert args.top_n_peaks == 150
    assert args.remove_precursor is True


def test_ms2_from_namespace_output_alias() -> None:
    parser = create_ms2_parser()
    ns = parser.parse_args(["/tmp/foo.d", "--output", "out.ms2"])
    args = Ms2Args.from_namespace(ns)
    assert args.output_file == "out.ms2"


def test_ms2_from_namespace_precision_flags() -> None:
    parser = create_ms2_parser()
    ns = parser.parse_args(
        ["/tmp/foo.d", "--mz-precision", "3", "--intensity-precision", "2"]
    )
    args = Ms2Args.from_namespace(ns)
    assert args.mz_precision == 3
    assert args.intensity_precision == 2


def test_ms2_from_namespace_precision_defaults_match_dataclass() -> None:
    """When --mz-precision is omitted the dataclass default should win."""
    parser = create_ms2_parser()
    ns = parser.parse_args(["/tmp/foo.d"])
    args = Ms2Args.from_namespace(ns)
    assert args.mz_precision == 5
    assert args.intensity_precision == 0


def test_mgf_from_namespace_drops_workers() -> None:
    parser = create_mgf_parser()
    ns = parser.parse_args(["/tmp/foo.d", "--workers", "4"])
    args = MgfArgs.from_namespace(ns)
    # workers is not a field on MgfArgs; it should be silently dropped
    assert not hasattr(args, "workers")


def test_mgf_from_namespace_drops_overwrite_and_verbose() -> None:
    parser = create_mgf_parser()
    ns = parser.parse_args(["/tmp/foo.d", "--overwrite", "--verbose"])
    args = MgfArgs.from_namespace(ns)
    assert not hasattr(args, "overwrite")
    assert not hasattr(args, "verbose")


def test_mzml_from_namespace_inverts_no_ms1() -> None:
    parser = create_mzml_parser()
    ns_with = parser.parse_args(["/tmp/foo.d", "--no-ms1"])
    args_with = MzmlArgs.from_namespace(ns_with)
    assert args_with.include_ms1 is False

    ns_without = parser.parse_args(["/tmp/foo.d"])
    args_without = MzmlArgs.from_namespace(ns_without)
    assert args_without.include_ms1 is True


def test_mzml_from_namespace_preserves_encoding_flags() -> None:
    parser = create_mzml_parser()
    ns = parser.parse_args(
        ["/tmp/foo.d", "--mz-encoding", "32", "--intensity-encoding", "64"]
    )
    args = MzmlArgs.from_namespace(ns)
    assert args.mz_encoding == 32
    assert args.intensity_encoding == 64


def test_mzml_from_namespace_preserves_compression_flags() -> None:
    parser = create_mzml_parser()
    ns = parser.parse_args(
        [
            "/tmp/foo.d",
            "--mz-compression",
            "none",
            "--intensity-compression",
            "zstd",
            "--mobility-compression",
            "numpress-linear",
        ]
    )
    args = MzmlArgs.from_namespace(ns)
    assert args.mz_compression == "none"
    assert args.intensity_compression == "zstd"
    assert args.mobility_compression == "numpress-linear"


def test_filter_fields_round_trip_through_namespace() -> None:
    parser = create_mzml_parser()
    ns = parser.parse_args(
        [
            "/tmp/foo.d",
            "--min-precursor-intensity",
            "1000",
            "--max-precursor-charge",
            "5",
            "--min-precursor-rt",
            "30.0",
        ]
    )
    args = MzmlArgs.from_namespace(ns)
    assert args.min_precursor_intensity == 1000
    assert args.max_precursor_charge == 5
    assert args.min_precursor_rt == 30.0
