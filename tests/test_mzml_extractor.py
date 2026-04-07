"""Tests for :mod:`tdfextractor.mzml_extractor`.

Helper-function tests run as fast unit tests. End-to-end mzML extraction
runs once per acquisition type via session-scoped fixtures defined in
``conftest.py``; the tests below parse the pre-written output and assert
on its structure.
"""

from __future__ import annotations

import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pytest

from tdfextractor.mzml_extractor import (
    _build_compression_dict,
    _build_encoding_dict,
    _resolve_compression,
    _resolve_encoding,
)

MZML_NS = {"mz": "http://psi.hupo.org/ms/mzml"}


# ---------------------------------------------------------------------------
# XML parsing helpers
# ---------------------------------------------------------------------------


def _iter_spectra(root: ET.Element):
    yield from root.iter("{http://psi.hupo.org/ms/mzml}spectrum")


def _ms_level(spec: ET.Element) -> int:
    for cv in spec.findall("mz:cvParam", MZML_NS):
        if cv.attrib.get("name") == "ms level":
            return int(cv.attrib["value"])
    return -1


def _binary_array_names(spec: ET.Element) -> list[str]:
    names: list[str] = []
    for bda in spec.iter("{http://psi.hupo.org/ms/mzml}binaryDataArray"):
        for cv in bda.findall("mz:cvParam", MZML_NS):
            nm = cv.attrib.get("name", "")
            if "array" in nm and "compression" not in nm and "bit" not in nm:
                names.append(nm)
                break
    return names


def _array_compressions(spec: ET.Element) -> dict[str, str | None]:
    out: dict[str, str | None] = {}
    for bda in spec.iter("{http://psi.hupo.org/ms/mzml}binaryDataArray"):
        arr_name: str | None = None
        comp_name: str | None = None
        for cv in bda.findall("mz:cvParam", MZML_NS):
            nm = cv.attrib.get("name", "")
            if arr_name is None and nm.endswith(" array"):
                arr_name = nm
            if "compression" in nm:
                comp_name = nm
        if arr_name is not None:
            out[arr_name] = comp_name
    return out


def _parse_mzml(path: Path) -> ET.Element:
    return ET.parse(str(path)).getroot()


# ---------------------------------------------------------------------------
# Fast helper-function tests (no I/O, no fixtures)
# ---------------------------------------------------------------------------


class TestHelpers:
    def test_resolve_compression_zlib_default(self) -> None:
        assert _resolve_compression("zlib") == "zlib"
        assert _resolve_compression(None) == "zlib"

    def test_resolve_compression_numpress_aliases(self) -> None:
        assert (
            _resolve_compression("numpress-linear")
            == "MS-Numpress linear prediction compression"
        )
        assert (
            _resolve_compression("numpress-slof")
            == "MS-Numpress short logged float compression"
        )
        assert (
            _resolve_compression("numpress-pic")
            == "MS-Numpress positive integer compression"
        )

    def test_resolve_compression_unknown_raises(self) -> None:
        with pytest.raises(ValueError):
            _resolve_compression("snappy")

    def test_resolve_encoding_bit_widths(self) -> None:
        assert _resolve_encoding(32) is np.float32
        assert _resolve_encoding(64) is np.float64
        assert _resolve_encoding(None) is np.float64
        with pytest.raises(ValueError):
            _resolve_encoding(16)  # type: ignore[arg-type]

    def test_build_compression_dict_keys(self) -> None:
        d = _build_compression_dict("zlib", "none", "zstd")
        assert d["m/z array"] == "zlib"
        assert d["intensity array"] == "none"
        assert d["mean inverse reduced ion mobility array"] == "zstd"

    def test_build_encoding_dict(self) -> None:
        d = _build_encoding_dict(64, 32)
        assert d["m/z array"] is np.float64
        assert d["intensity array"] is np.float32
        assert d["mean inverse reduced ion mobility array"] is np.float64


# ---------------------------------------------------------------------------
# Slow end-to-end tests (consume session-scoped fixtures)
# ---------------------------------------------------------------------------


@pytest.mark.slow
@pytest.mark.parametrize(
    "fixture_name",
    ["mzml_dda_output", "mzml_dia_output", "mzml_prm_output"],
)
def test_default_write_has_ms1_and_ms2(request, fixture_name: str) -> None:
    out = request.getfixturevalue(fixture_name)
    root = _parse_mzml(out)
    specs = list(_iter_spectra(root))
    n_ms1 = sum(1 for s in specs if _ms_level(s) == 1)
    n_ms2 = sum(1 for s in specs if _ms_level(s) == 2)
    assert n_ms1 > 0, "expected at least one MS1 spectrum"
    assert n_ms2 > 0, "expected at least one MS2 spectrum"


@pytest.mark.slow
@pytest.mark.parametrize(
    "fixture_name",
    ["mzml_dda_output", "mzml_dia_output", "mzml_prm_output"],
)
def test_ms1_carries_mean_ion_mobility_array(request, fixture_name: str) -> None:
    out = request.getfixturevalue(fixture_name)
    root = _parse_mzml(out)
    hits = sum(
        1
        for s in _iter_spectra(root)
        if _ms_level(s) == 1
        and "mean inverse reduced ion mobility array" in _binary_array_names(s)
    )
    assert hits > 0, "no MS1 spectra carried the mean-IM mobility array"


@pytest.mark.slow
def test_no_ms1_flag_skips_ms1(mzml_dia_no_ms1_output: Path) -> None:
    root = _parse_mzml(mzml_dia_no_ms1_output)
    for spec in _iter_spectra(root):
        assert _ms_level(spec) != 1, "include_ms1=False should skip all MS1 spectra"


@pytest.mark.slow
def test_compression_param_plumbed_to_mz_array(
    mzml_dia_no_mz_compression_output: Path,
) -> None:
    root = _parse_mzml(mzml_dia_no_mz_compression_output)
    saw_mz = False
    for spec in _iter_spectra(root):
        comps = _array_compressions(spec)
        if "m/z array" not in comps:
            continue
        assert (
            comps.get("m/z array") == "no compression"
        ), f"expected uncompressed m/z array, got {comps!r}"
        saw_mz = True
        break
    assert saw_mz, "no spectrum with an m/z array was found"
