import time
from datetime import datetime
from pathlib import Path

from tdfpy import timsdata
from tdfpy.pandas_tdf import PandasTdf
from serenipy.ms2 import Ms2Spectra, to_ms2

from .constants import MS2_VERSION
from .string_templates import header_ms2_template
from .utils import calculate_mass, map_precursor_to_ip2_scan_number


def get_ms2_content(analysis_dir: str, include_spectra=True):
    start_time = time.time()

    with timsdata.timsdata_connect(analysis_dir) as td:

        pd_tdf = PandasTdf(str(Path(analysis_dir) / 'analysis.tdf'))
        precursors_df = pd_tdf.precursors
        frames_df = pd_tdf.frames
        precursor_to_scan_number = map_precursor_to_ip2_scan_number(precursors_df, frames_df)
        pasef_frames_msms_info_df = pd_tdf.pasef_frame_msms_info

        parent_id_to_rt = {int(frame_id): rt for frame_id, rt in frames_df[['Id', 'Time']].values}
        precursor_id_to_collision_energy = {int(prec_id): ce for prec_id, ce in
                                            pasef_frames_msms_info_df[['Precursor', 'CollisionEnergy']].values}

        precursors_df.dropna(subset=['MonoisotopicMz', 'Charge'], inplace=True)


        analysis_load_time = time.time() - start_time
        print(analysis_load_time)
        start_time = time.time()

        for _, precursor_row in precursors_df.iterrows():

            precursor_id = int(precursor_row['Id'])
            parent_id = int(precursor_row['Parent'])
            charge = int(precursor_row['Charge'])
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
            ms2_spectra.prec_intensity = round(prec_intensity, 1)
            ms2_spectra.ook0 = round(ook0, 4)
            ms2_spectra.ccs = round(ccs, 4)
            ms2_spectra.rt = round(parent_id_to_rt[parent_id], 4)
            ms2_spectra.ce = round(precursor_id_to_collision_energy[precursor_id], 1)

            if include_spectra:
                mz_arr, int_arr = td.readPasefMsMs([precursor_id])[precursor_id]
                ms2_spectra.mz_spectra = mz_arr
                ms2_spectra.intensity_spectra = int_arr

            yield ms2_spectra

        spectra_time = time.time() - start_time
        print(spectra_time)


def generate_header(analysis_dir: str):
    pd_tdf = PandasTdf(str(Path(analysis_dir) / 'analysis.tdf'))
    precursors_df = pd_tdf.precursors
    frames_df = pd_tdf.frames
    precursor_to_scan_number = map_precursor_to_ip2_scan_number(precursors_df, frames_df)

    ms2_header = header_ms2_template.format(version=MS2_VERSION,
                                            date_of_creation=str(datetime.now().strftime("%B %d, %Y %H:%M")),
                                            first_scan=list(precursor_to_scan_number.values())[0],
                                            last_scan=list(precursor_to_scan_number.values())[-1])

    return ms2_header


def write_ms2_file(analysis_dir: str, include_spectra=True, output_file=None):
    if output_file is None:
        output_file = str(Path(analysis_dir) / Path(analysis_dir).stem) + '.ms2'

    ms2_header = generate_header(analysis_dir)
    ms2_spectra = list(get_ms2_content(analysis_dir, include_spectra))

    ms2_content = to_ms2([ms2_header], ms2_spectra)

    with open(output_file, 'w') as file:
        file.write(ms2_content)
