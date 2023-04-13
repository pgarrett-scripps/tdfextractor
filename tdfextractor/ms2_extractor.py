import logging
import time
from datetime import datetime
from pathlib import Path
from typing import List

from tdfpy import timsdata
from tdfpy.pandas_tdf import PandasTdf
from serenipy.ms2 import Ms2Spectra, to_ms2
from tqdm import tqdm

from .constants import MS2_VERSION
from .string_templates import header_ms2_template
from .utils import calculate_mass, map_precursor_to_ip2_scan_number

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def batch_iterator(input_list: List, batch_size: int):
    for i in range(0, len(input_list), batch_size):
        yield input_list[i:i + batch_size]


def get_ms2_content(analysis_dir: str, include_spectra: bool = True, batch_size: int = 100,
                    remove_charge1: bool = True, remove_empty_spectra: bool = True,
                    min_intensity: float = 0):
    """
    Retrieves MS/MS spectra from a given analysis directory.

    Args:
        analysis_dir (str): Path to the analysis directory.
        include_spectra (bool, optional): Whether to include the actual MS/MS spectra in the output. Defaults to True.
        batch_size (int, optional): The number of precursors to process in a batch. Defaults to 100.
        remove_charge1 (bool, optional): Whether to exclude precursors with charge 1. Defaults to True.
        remove_empty_spectra (bool, optional): Whether to exclude MS/MS spectra with no spectra values. Defaults to True.
        min_intensity (float, optional): Filters out MS/MS spectral data based on intensity (inclusive). Defaults to 9.

    Returns:
        List[Ms2Spectra]: A list of Ms2Spectra objects representing MS/MS spectra.

    Raises:
        FileNotFoundError: If the analysis directory does not exist.
    """
    logger.info(f'Starting to process {analysis_dir}')

    with timsdata.timsdata_connect(analysis_dir) as td:

        analysis_tdf_path = str(Path(analysis_dir) / 'analysis.tdf')
        pd_tdf = PandasTdf(analysis_tdf_path)

        precursors_df = pd_tdf.precursors
        frames_df = pd_tdf.frames

        precursor_to_scan_number = map_precursor_to_ip2_scan_number(precursors_df, frames_df)
        parent_id_to_rt = {int(frame_id): rt for frame_id, rt in frames_df[['Id', 'Time']].values}

        precursor_id_to_msms_info = {int(prec_id): (ce, iso_width, iso_mz) for prec_id, ce, iso_width, iso_mz in
                pd_tdf.pasef_frame_msms_info[['Precursor', 'CollisionEnergy', 'IsolationWidth', 'IsolationMz']].values}

        precursors_df.dropna(subset=['MonoisotopicMz', 'Charge'], inplace=True)

        for precursor_batch in tqdm(list(batch_iterator(input_list=list(precursors_df.iterrows()), batch_size=batch_size)),
                                    desc='Generating MS2 Spectra'):
            pasef_ms_ms = None
            if include_spectra:
                pasef_ms_ms = td.readPasefMsMs([int(precursor_row['Id']) for _, precursor_row in precursor_batch])

            for _, precursor_row in precursor_batch:

                precursor_id = int(precursor_row['Id'])
                parent_id = int(precursor_row['Parent'])
                charge = int(precursor_row['Charge'])

                if remove_charge1 is True and charge == 1:
                    continue

                ip2_scan_number = precursor_to_scan_number[precursor_id]
                ook0 = td.scanNumToOneOverK0(parent_id, [precursor_row['ScanNumber']])[0]
                ccs = timsdata.oneOverK0ToCCSforMz(ook0, charge, precursor_row['MonoisotopicMz'])
                mz = precursor_row['MonoisotopicMz']
                prec_intensity = precursor_row['Intensity']
                mass = calculate_mass(mz, charge)

                ms2_spectra = Ms2Spectra(low_scan=ip2_scan_number,
                                         high_scan=ip2_scan_number,
                                         mz=mz,
                                         mass=mass,
                                         charge=charge,
                                         info={},
                                         mz_spectra=[],
                                         intensity_spectra=[],
                                         charge_spectra=[])

                ms2_spectra.parent_id = parent_id
                ms2_spectra.precursor_id = precursor_id
                ms2_spectra.prec_intensity = int(prec_intensity)
                ms2_spectra.ook0 = round(ook0, 4)
                ms2_spectra.ccs = round(ccs, 1)
                ms2_spectra.rt = round(parent_id_to_rt[parent_id], 2)
                ms2_spectra.ce = round(precursor_id_to_msms_info[precursor_id][0], 1)
                ms2_spectra.iso_width = round(precursor_id_to_msms_info[precursor_id][1], 1)
                ms2_spectra.iso_mz = round(precursor_id_to_msms_info[precursor_id][2], 4)

                if include_spectra:
                    ms2_spectra_data = list(zip(pasef_ms_ms[precursor_id][0], pasef_ms_ms[precursor_id][1]))

                    if min_intensity != 0:
                        ms2_spectra_data = [data for data in ms2_spectra_data if data[1] <= min_intensity]

                    ms2_spectra.mz_spectra = [data[0] for data in ms2_spectra_data]
                    ms2_spectra.intensity_spectra = [int(data[1]) for data in ms2_spectra_data]

                    assert len(ms2_spectra.mz_spectra) == len(ms2_spectra.intensity_spectra)

                    if remove_empty_spectra is True and len(ms2_spectra.mz_spectra) == 0:
                        continue

                yield ms2_spectra


def generate_header(analysis_dir: str):
    """
    Generates a header string for MS2 data using information from the analysis file.

    Args:
        analysis_dir: The directory path where the analysis file `analysis.tdf` is located.

    Returns:
        The generated header string for MS2 data.
    """

    pd_tdf = PandasTdf(str(Path(analysis_dir) / 'analysis.tdf'))
    precursors_df = pd_tdf.precursors
    frames_df = pd_tdf.frames
    precursor_to_scan_number = map_precursor_to_ip2_scan_number(precursors_df, frames_df)

    ms2_header = header_ms2_template.format(version=MS2_VERSION,
                                            date_of_creation=str(datetime.now().strftime("%B %d, %Y %H:%M")),
                                            first_scan=list(precursor_to_scan_number.values())[0],
                                            last_scan=list(precursor_to_scan_number.values())[-1])

    return ms2_header


def write_ms2_file(analysis_dir: str, output_file: str = None, include_spectra: bool = True, batch_size: int = 1000,
                   remove_charge1: bool = True, remove_empty_spectra: bool = True, min_intensity: float = 0):
    """
    Writes MS/MS spectra data to an MS2 file.

    Args:
        analysis_dir (str): Path to the analysis directory.
        include_spectra (bool, optional): Whether to include the actual MS/MS spectra in the output. Defaults to True.
        output_file (str, optional): The output file name. Defaults to None.
        batch_size (int, optional): The number of precursors to process in a batch. Defaults to 1000.
        remove_charge1 (bool, optional): Whether to exclude precursors with charge 1. Defaults to True.
        remove_empty_spectra (bool, optional): Whether to exclude MS/MS spectra with no spectra values. Defaults to True.
        min_intensity (float, optional): Filters out MS/MS spectral data based on intensity (inclusive). Defaults to 0.

    Returns:
        None.

    Raises:
        FileNotFoundError: If the analysis directory does not exist.
    """
    start_time = time.time()

    if output_file is None:
        output_file = str(Path(analysis_dir) / Path(analysis_dir).stem) + '.ms2'

    logger.info(f'Arguments: analysis_dir: {analysis_dir}, output_file: {output_file}, include_spectra: {include_spectra}, '
    f'batch_size: {batch_size}, remove_charge1: {remove_charge1}, remove_empty_spectra: {remove_empty_spectra}, '
    f'min_intensity: {min_intensity}')

    logger.info('Creating Ms2 Header')
    ms2_header = generate_header(analysis_dir)

    logger.info('Generating Ms2 Spectra')
    ms2_spectra = get_ms2_content(analysis_dir=analysis_dir, include_spectra=include_spectra, batch_size=batch_size,
                                  remove_charge1=remove_charge1, remove_empty_spectra=remove_empty_spectra,
                                  min_intensity=min_intensity)

    logger.info('Creating MS2 Contents')
    ms2_content = to_ms2([ms2_header], ms2_spectra)

    logger.info('Writing Contents To File')
    with open(output_file, 'w') as file:
        file.write(ms2_content)

    logger.info(f'Total Time: {round(time.time() - start_time, 2)} seconds')
