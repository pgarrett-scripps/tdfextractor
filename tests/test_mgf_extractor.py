"""Tests for :mod:`tdfextractor.mgf_exctractor`.

The MGF writer only supports DDA .d folders (it goes through
``get_ms2_dda_content`` which calls ``readPasefMsMs``), so all tests use
the DDA fixture shared with the MS2 tests via ``conftest.py``.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from tdfextractor import MgfArgs, write_mgf_file


def _parse_mgf(path: Path) -> list[dict]:
    """Minimal MGF parser. Returns one dict per BEGIN IONS / END IONS block.

    Each dict carries the metadata key/value lines plus a ``peaks`` list of
    ``(mz, intensity)`` tuples. Header lines outside of any BEGIN IONS block
    (``INSTRUMENT=...``, ``MASS=...``) are ignored.
    """
    entries: list[dict] = []
    current: dict | None = None
    peaks: list[tuple[float, float]] = []
    with open(path, encoding="utf-8") as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line == "BEGIN IONS":
                current = {}
                peaks = []
                continue
            if line == "END IONS":
                assert current is not None
                current["peaks"] = peaks
                entries.append(current)
                current = None
                continue
            if current is None:
                continue
            if "=" in line:
                key, value = line.split("=", 1)
                current[key] = value
            else:
                mz_str, int_str = line.split()
                peaks.append((float(mz_str), float(int_str)))
    return entries


@pytest.mark.slow
def test_mgf_file_created(mgf_dda_output: Path) -> None:
    assert mgf_dda_output.exists()
    assert mgf_dda_output.stat().st_size > 0


@pytest.mark.slow
def test_mgf_has_spectra(mgf_dda_output: Path) -> None:
    entries = _parse_mgf(mgf_dda_output)
    assert len(entries) > 0


@pytest.mark.slow
def test_mgf_entries_have_required_fields(mgf_dda_output: Path) -> None:
    entries = _parse_mgf(mgf_dda_output)
    sample = entries[:10]
    assert sample, "expected at least one MGF entry"
    for entry in sample:
        assert "TITLE" in entry
        assert "RTINSECONDS" in entry
        assert "PEPMASS" in entry
        assert "CHARGE" in entry
        # PEPMASS = "<mz> <intensity>"
        pepmass_parts = entry["PEPMASS"].split()
        assert len(pepmass_parts) == 2
        float(pepmass_parts[0])  # parses
        float(pepmass_parts[1])  # parses
        # CHARGE format: "<n>+"
        assert re.match(r"^\d+\+$", entry["CHARGE"])
        # RTINSECONDS is a float
        float(entry["RTINSECONDS"])


@pytest.mark.slow
def test_mgf_peaks_sorted_by_mz(mgf_dda_output: Path) -> None:
    """Each spectrum's peaks should be ascending in m/z."""
    entries = _parse_mgf(mgf_dda_output)
    non_empty = [e for e in entries if e["peaks"]]
    assert non_empty, "expected at least one non-empty MGF block"
    for entry in non_empty[:10]:
        mzs = [mz for mz, _ in entry["peaks"]]
        assert mzs == sorted(mzs), "MGF peaks should be sorted ascending by m/z"


@pytest.mark.slow
def test_mgf_title_includes_native_id(mgf_dda_output: Path) -> None:
    entries = _parse_mgf(mgf_dda_output)
    assert entries
    for entry in entries[:5]:
        title = entry["TITLE"]
        assert "NativeID=" in title
        assert "frame=" in title
        assert "scan=" in title


@pytest.mark.slow
def test_mgf_default_precision(mgf_dda_output: Path) -> None:
    """Default mz_precision=5 -> peak m/z should have 5 decimal places."""
    text = mgf_dda_output.read_text(encoding="utf-8")
    # First peak line we hit (numeric "<mz> <int>" pattern, not a KEY=value)
    peak_line = next(
        line
        for line in text.splitlines()
        if re.match(r"^\d+\.\d+ \d+", line)
    )
    mz_str = peak_line.split()[0]
    assert len(mz_str.split(".")[1]) == 5


@pytest.mark.slow
def test_mgf_custom_precision(tmp_path_factory, dda_d_folder: Path) -> None:
    """Custom mz_precision/intensity_precision flow through to output.

    Uses a tight ``max_precursor_rt`` window to keep this extra extraction
    cheap.
    """
    out_dir = tmp_path_factory.mktemp("mgf_precision")
    out = out_dir / "precision.mgf"
    write_mgf_file(
        MgfArgs(
            analysis_dir=str(dda_d_folder),
            output_file=str(out),
            mz_precision=3,
            intensity_precision=2,
            max_precursor_rt=2403.0,
        )
    )
    text = out.read_text(encoding="utf-8")
    peak_line = next(
        line
        for line in text.splitlines()
        if re.match(r"^\d+\.\d+ \d+\.\d+$", line)
    )
    mz_str, int_str = peak_line.split()
    assert len(mz_str.split(".")[1]) == 3
    assert len(int_str.split(".")[1]) == 2
