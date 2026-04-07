"""Shared pytest fixtures for the tdfextractor test suite.

The session-scoped output fixtures defined here run each writer exactly once
per test session and return the path to the produced file. Tests that
parse-and-assert against the output share the same extraction, which keeps
the slow extractions from re-running for every individual test method.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from tdfextractor import (
    MgfArgs,
    Ms2Args,
    MzmlArgs,
    write_mgf_file,
    write_ms2_file,
    write_mzml_file,
)

DATA_DIR = Path(__file__).parent / "data"
DDA_D = DATA_DIR / "200ngHeLaPASEF_1min.d"
DIA_D = DATA_DIR / "example_dia.d"
PRM_D = DATA_DIR / "example_prm.d"

logging.basicConfig(level=logging.WARNING)


# ---------------------------------------------------------------------------
# .d folder path fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(scope="session")
def dda_d_folder() -> Path:
    return DDA_D


@pytest.fixture(scope="session")
def dia_d_folder() -> Path:
    return DIA_D


@pytest.fixture(scope="session")
def prm_d_folder() -> Path:
    return PRM_D


# ---------------------------------------------------------------------------
# mzML session-scoped outputs
# ---------------------------------------------------------------------------


@pytest.fixture(scope="session")
def mzml_dda_output(tmp_path_factory, dda_d_folder: Path) -> Path:
    out_dir = tmp_path_factory.mktemp("mzml_dda")
    out = out_dir / (dda_d_folder.stem + ".mzML")
    write_mzml_file(MzmlArgs(analysis_dir=str(dda_d_folder), output_file=str(out)))
    return out


@pytest.fixture(scope="session")
def mzml_dia_output(tmp_path_factory, dia_d_folder: Path) -> Path:
    out_dir = tmp_path_factory.mktemp("mzml_dia")
    out = out_dir / (dia_d_folder.stem + ".mzML")
    write_mzml_file(MzmlArgs(analysis_dir=str(dia_d_folder), output_file=str(out)))
    return out


@pytest.fixture(scope="session")
def mzml_prm_output(tmp_path_factory, prm_d_folder: Path) -> Path:
    out_dir = tmp_path_factory.mktemp("mzml_prm")
    out = out_dir / (prm_d_folder.stem + ".mzML")
    write_mzml_file(MzmlArgs(analysis_dir=str(prm_d_folder), output_file=str(out)))
    return out


@pytest.fixture(scope="session")
def mzml_dia_no_ms1_output(tmp_path_factory, dia_d_folder: Path) -> Path:
    """include_ms1=False variant on the small DIA fixture (cheap)."""
    out_dir = tmp_path_factory.mktemp("mzml_dia_no_ms1")
    out = out_dir / (dia_d_folder.stem + "_no_ms1.mzML")
    write_mzml_file(
        MzmlArgs(
            analysis_dir=str(dia_d_folder),
            output_file=str(out),
            include_ms1=False,
        )
    )
    return out


@pytest.fixture(scope="session")
def mzml_dia_no_mz_compression_output(tmp_path_factory, dia_d_folder: Path) -> Path:
    """mz_compression='none' variant on the small DIA fixture (cheap)."""
    out_dir = tmp_path_factory.mktemp("mzml_dia_nocomp")
    out = out_dir / (dia_d_folder.stem + "_nocomp.mzML")
    write_mzml_file(
        MzmlArgs(
            analysis_dir=str(dia_d_folder),
            output_file=str(out),
            mz_compression="none",
            intensity_compression="zlib",
            mobility_compression="zlib",
            mz_encoding=64,
            intensity_encoding=32,
        )
    )
    return out


# ---------------------------------------------------------------------------
# MS2 session-scoped output (DDA only — get_ms2_dda_content is DDA-only)
# ---------------------------------------------------------------------------


@pytest.fixture(scope="session")
def ms2_dda_output(tmp_path_factory, dda_d_folder: Path) -> Path:
    out_dir = tmp_path_factory.mktemp("ms2_dda")
    out = out_dir / (dda_d_folder.stem + ".ms2")
    write_ms2_file(Ms2Args(analysis_dir=str(dda_d_folder), output_file=str(out)))
    return out


# ---------------------------------------------------------------------------
# MGF session-scoped output (DDA only)
# ---------------------------------------------------------------------------


@pytest.fixture(scope="session")
def mgf_dda_output(tmp_path_factory, dda_d_folder: Path) -> Path:
    out_dir = tmp_path_factory.mktemp("mgf_dda")
    out = out_dir / (dda_d_folder.stem + ".mgf")
    write_mgf_file(MgfArgs(analysis_dir=str(dda_d_folder), output_file=str(out)))
    return out
