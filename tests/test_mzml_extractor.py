"""Tests for :mod:`tdfextractor.mzml_extractor`.

These tests exercise the DDA, DIA, and PRM writer code paths against the
small sample ``.d`` folders copied over from ``tdfpy``. The mzML output is
validated purely through the XML tree (no external reader dependency) so
the test suite only requires packages already declared as runtime deps.
"""

import logging
import os
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np

from tdfextractor.mzml_extractor import (
    _build_compression_dict,
    _build_encoding_dict,
    _resolve_compression,
    _resolve_encoding,
    write_mzml_file,
)


DATA_DIR = Path(__file__).parent / "data"
DDA_D = DATA_DIR / "200ngHeLaPASEF_1min.d"
DIA_D = DATA_DIR / "example_dia.d"
PRM_D = DATA_DIR / "example_prm.d"

MZML_NS = {"mz": "http://psi.hupo.org/ms/mzml"}

logging.basicConfig(level=logging.WARNING)


def _iter_spectra(root: ET.Element):
    for spec in root.iter("{http://psi.hupo.org/ms/mzml}spectrum"):
        yield spec


def _ms_level(spec: ET.Element) -> int:
    for cv in spec.findall("mz:cvParam", MZML_NS):
        if cv.attrib.get("name") == "ms level":
            return int(cv.attrib["value"])
    return -1


def _binary_array_names(spec: ET.Element):
    names = []
    for bda in spec.iter("{http://psi.hupo.org/ms/mzml}binaryDataArray"):
        for cv in bda.findall("mz:cvParam", MZML_NS):
            nm = cv.attrib.get("name", "")
            if (
                "array" in nm
                and "compression" not in nm
                and "bit" not in nm
            ):
                names.append(nm)
                break
    return names


def _array_compressions(spec: ET.Element):
    """Return (array_name -> compression_param_name) for one spectrum."""

    out = {}
    for bda in spec.iter("{http://psi.hupo.org/ms/mzml}binaryDataArray"):
        arr_name = None
        comp_name = None
        for cv in bda.findall("mz:cvParam", MZML_NS):
            nm = cv.attrib.get("name", "")
            if arr_name is None and nm.endswith(" array"):
                arr_name = nm
            if "compression" in nm:
                comp_name = nm
        if arr_name is not None:
            out[arr_name] = comp_name
    return out


def _parse_mzml(path: str) -> ET.Element:
    tree = ET.parse(path)
    return tree.getroot()


class TestHelpers(unittest.TestCase):
    def test_resolve_compression_zlib_default(self):
        self.assertEqual(_resolve_compression("zlib"), "zlib")
        self.assertEqual(_resolve_compression(None), "zlib")

    def test_resolve_compression_numpress_aliases(self):
        self.assertEqual(
            _resolve_compression("numpress-linear"),
            "MS-Numpress linear prediction compression",
        )
        self.assertEqual(
            _resolve_compression("numpress-slof"),
            "MS-Numpress short logged float compression",
        )
        self.assertEqual(
            _resolve_compression("numpress-pic"),
            "MS-Numpress positive integer compression",
        )

    def test_resolve_compression_unknown_raises(self):
        with self.assertRaises(ValueError):
            _resolve_compression("snappy")

    def test_resolve_encoding_bit_widths(self):
        self.assertIs(_resolve_encoding(32), np.float32)
        self.assertIs(_resolve_encoding(64), np.float64)
        self.assertIs(_resolve_encoding(None), np.float64)
        with self.assertRaises(ValueError):
            _resolve_encoding(16)

    def test_build_compression_dict_keys(self):
        d = _build_compression_dict("zlib", "none", "zstd")
        self.assertEqual(d["m/z array"], "zlib")
        self.assertEqual(d["intensity array"], "none")
        self.assertEqual(
            d["mean inverse reduced ion mobility array"], "zstd"
        )

    def test_build_encoding_dict(self):
        d = _build_encoding_dict(64, 32)
        self.assertIs(d["m/z array"], np.float64)
        self.assertIs(d["intensity array"], np.float32)
        self.assertIs(d["mean inverse reduced ion mobility array"], np.float64)


class _AcquisitionWriterMixin:
    """Shared assertions across DDA / DIA / PRM extractors."""

    d_folder: Path
    expect_mobility_on_ms1: bool = True

    def _write(self, **kwargs):
        out = Path(self.tmpdir.name) / (self.d_folder.stem + ".mzML")
        write_mzml_file(str(self.d_folder), output_file=str(out), **kwargs)
        self.assertTrue(out.exists(), f"mzML output not created at {out}")
        self.assertGreater(out.stat().st_size, 0)
        return out

    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmpdir.cleanup)

    def test_default_write(self):
        out = self._write()
        root = _parse_mzml(str(out))
        specs = list(_iter_spectra(root))
        n_ms1 = sum(1 for s in specs if _ms_level(s) == 1)
        n_ms2 = sum(1 for s in specs if _ms_level(s) == 2)
        self.assertGreater(n_ms1, 0, "expected at least one MS1 spectrum")
        self.assertGreater(n_ms2, 0, "expected at least one MS2 spectrum")

        if self.expect_mobility_on_ms1:
            # At least one MS1 spectrum should carry the IM array.
            hits = 0
            for s in specs:
                if _ms_level(s) != 1:
                    continue
                if (
                    "mean inverse reduced ion mobility array"
                    in _binary_array_names(s)
                ):
                    hits += 1
            self.assertGreater(
                hits, 0, "no MS1 spectra carried the mean-IM mobility array"
            )

    def test_no_ms1_flag(self):
        out = self._write(include_ms1=False)
        root = _parse_mzml(str(out))
        for spec in _iter_spectra(root):
            self.assertNotEqual(
                _ms_level(spec),
                1,
                "--no-ms1 flag should skip all MS1 spectra",
            )

    def test_compression_param_plumbed_to_mz_array(self):
        out = self._write(
            mz_compression="none",
            intensity_compression="zlib",
            mobility_compression="zlib",
            mz_encoding=64,
            intensity_encoding=32,
        )
        root = _parse_mzml(str(out))
        # Grab the first spectrum that has an m/z array and confirm it is
        # flagged as "no compression".
        saw_mz = False
        for spec in _iter_spectra(root):
            comps = _array_compressions(spec)
            if "m/z array" not in comps:
                continue
            self.assertEqual(
                comps.get("m/z array"),
                "no compression",
                f"expected m/z array to be uncompressed, got {comps!r}",
            )
            saw_mz = True
            break
        self.assertTrue(saw_mz, "no spectrum with an m/z array was found")


class TestDDAExtractor(_AcquisitionWriterMixin, unittest.TestCase):
    d_folder = DDA_D


class TestDIAExtractor(_AcquisitionWriterMixin, unittest.TestCase):
    d_folder = DIA_D


class TestPRMExtractor(_AcquisitionWriterMixin, unittest.TestCase):
    d_folder = PRM_D


if __name__ == "__main__":
    unittest.main()
