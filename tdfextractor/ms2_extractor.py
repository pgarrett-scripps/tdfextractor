from datetime import datetime
from pathlib import Path

from tdfpy import timsdata
from tdfpy.pandas_tdf import PandasTdf

from .constants import MS2_VERSION
from .string_templates import header_ms2_template, dda_ms2_scan_template, create_mz_int_spectra
from .utils import calculate_mass, map_precursor_to_ip2_scan_number


def get_ms2_content(analysis_dir: str, include_spectra=True):
    td = timsdata.TimsData(analysis_dir)

    pd_tdf = PandasTdf(str(Path(analysis_dir) / 'analysis.tdf'))
    precursors_df = pd_tdf.precursors
    frames_df = pd_tdf.frames
    precursor_to_scan_number = map_precursor_to_ip2_scan_number(precursors_df, frames_df)
    pasef_frames_msms_info_df = pd_tdf.pasef_frame_msms_info

    ms2_header = header_ms2_template.format(version=MS2_VERSION,
                                            date_of_creation=str(datetime.now().strftime("%B %d, %Y %H:%M")),
                                            first_scan=list(precursor_to_scan_number.values())[0],
                                            last_scan=list(precursor_to_scan_number.values())[-1])
    yield (ms2_header)

    parent_id_to_rt = {int(frame_id): rt for frame_id, rt in frames_df[['Id', 'Time']].values}
    precursor_id_to_collision_energy = {int(prec_id): ce for prec_id, ce in
                                        pasef_frames_msms_info_df[['Precursor', 'CollisionEnergy']].values}

    precursors_df.dropna(subset=['MonoisotopicMz', 'Charge'], inplace=True)
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

        dda_ms2_scan = dda_ms2_scan_template.format(scan_id=ip2_scan_number,
                                                    prc_mass_mz=mz,
                                                    parent_id=parent_id,
                                                    prc_id=precursor_id,
                                                    prc_int=prec_intensity,
                                                    ook0=ook0,
                                                    ccs=ccs,
                                                    ce=precursor_id_to_collision_energy[precursor_id],
                                                    ret_time=parent_id_to_rt[parent_id],
                                                    cs=charge,
                                                    prc_mass=mass)

        if include_spectra:
            mz_arr, int_arr = td.readPasefMsMs([precursor_id])[precursor_id]
            if len(mz_arr) > 0:
                dda_ms2_spectra = create_mz_int_spectra(mz_arr, int_arr)
                yield (dda_ms2_scan)
                yield (dda_ms2_spectra)
        else:
            yield (dda_ms2_scan)


def write_ms2_file(analysis_dir: str, include_spectra=True, output_file=None):

    if output_file is None:
        output_file = str(Path(analysis_dir) / Path(analysis_dir).stem) + '.ms2'

    ms2_content = ''.join(get_ms2_content(analysis_dir, include_spectra))
    with open(output_file, 'w') as file:
        file.write(ms2_content)