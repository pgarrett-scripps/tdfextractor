"""Tests for :mod:`tdfextractor.ms2_extractor`.

End-to-end MS2 extraction runs once per session via the ``ms2_dda_output``
fixture defined in ``conftest.py``; the tests below parse the pre-written
output and assert on its structure.
"""

from __future__ import annotations

from pathlib import Path

import pytest


@pytest.mark.slow
def test_ms2_file_created(ms2_dda_output: Path) -> None:
    assert ms2_dda_output.exists()
    assert ms2_dda_output.stat().st_size > 0


@pytest.mark.slow
def test_ms2_header_present(ms2_dda_output: Path) -> None:
    text = ms2_dda_output.read_text(encoding="utf-8")
    assert text.startswith("H\tExtractor\tTimsTOF_extractor")
    assert "H\tAcquisitionMethod\tData-Dependent" in text
    assert "H\tScanType\tMS2" in text


@pytest.mark.slow
def test_ms2_header_records_default_precision(ms2_dda_output: Path) -> None:
    """Default mz/intensity precision values should appear in the header."""
    text = ms2_dda_output.read_text(encoding="utf-8")
    assert "H\tMzPrecision\t5" in text
    assert "H\tIntensityPrecision\t0" in text


@pytest.mark.slow
def test_ms2_has_s_blocks(ms2_dda_output: Path) -> None:
    """Each MS2 spectrum block starts with a tab-prefixed 'S' line."""
    text = ms2_dda_output.read_text(encoding="utf-8")
    s_lines = [line for line in text.splitlines() if line.startswith("S\t")]
    assert len(s_lines) > 0


@pytest.mark.slow
def test_ms2_has_z_blocks(ms2_dda_output: Path) -> None:
    """Each MS2 spectrum should also have a 'Z' (charge) line."""
    text = ms2_dda_output.read_text(encoding="utf-8")
    z_lines = [line for line in text.splitlines() if line.startswith("Z\t")]
    assert len(z_lines) > 0
