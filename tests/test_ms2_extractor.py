import logging
import unittest
from pathlib import Path
import os

from tdfextractor.ms2_extractor import write_ms2_file
from tdfextractor.utils import get_ms2_dda_content

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

