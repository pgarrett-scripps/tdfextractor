import logging
import unittest
from pathlib import Path
import os

from tdfextractor.ms2_extractor import get_ms2_content, write_ms2_file

d_folder = str(Path('tests') / 'data' / '200ngHeLaPASEF_1min.d')

logging.basicConfig(
    level=logging.INFO)


class TestMs2Extractor(unittest.TestCase):

    def test_write_ms2(self):
        if os.path.exists('test.ms2'):
            os.remove('test.ms2')

        write_ms2_file(d_folder, output_file='test.ms2')
        self.assertTrue(os.path.exists('test.ms2'))


    def test_write_ms2_folder(self):
        ms2_file_path = os.path.join(d_folder, '200ngHeLaPASEF_1min.ms2')

        if os.path.exists(ms2_file_path):
            os.remove(ms2_file_path)

        write_ms2_file(d_folder)
        self.assertTrue(os.path.exists(ms2_file_path))

    def test_ms2_contents(self):

        ms2_spectra = None
        for spectra in get_ms2_content(d_folder):
            ms2_spectra = spectra
            break

        self.assertAlmostEqual(ms2_spectra.charge, 2)
        self.assertAlmostEqual(ms2_spectra.mz, 1292.63706188582, 4)
        self.assertAlmostEqual(ms2_spectra.prec_intensity, 3603.0)
        self.assertAlmostEqual(ms2_spectra.rt, 2400.83148655562, 1)
        self.assertEqual(ms2_spectra.precursor_id, 1)
        self.assertEqual(ms2_spectra.parent_id, 1)
        self.assertAlmostEqual(ms2_spectra.mz_spectra[0], 113.6913703181829, 4)
        self.assertAlmostEqual(ms2_spectra.intensity_spectra[0], 14.0, 1)
        self.assertAlmostEqual(ms2_spectra.mz_spectra[-1], 1699.749570812201, 4)
        self.assertAlmostEqual(ms2_spectra.intensity_spectra[-1], 15.0, 1)
        self.assertAlmostEqual(ms2_spectra.mass, 2584.2668, 4)
        self.assertEqual(len(ms2_spectra.charge_spectra), 0)
