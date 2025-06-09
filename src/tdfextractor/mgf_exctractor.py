"""
ms2_extractor defines functions for generating ms2 files from DDA and PRM based .D folders
"""

import logging
import time
from datetime import datetime
from pathlib import Path
from typing import Generator, List, Optional
import argparse

import pandas as pd
from tdfpy import timsdata
from tdfpy.pandas_tdf import PandasTdf
from serenipy.ms2 import Ms2Spectra, to_ms2
from tqdm import tqdm

from .constants import MS2_VERSION
from .string_templates import MS2_HEADER
from .utils import calculate_mass, map_precursor_to_ip2_scan_number
import numpy as np

logger = logging.getLogger(__name__)


def batch_iterator(input_list: List, batch_size: int):
    """
    Creates an iterator that returns a batch of elements from the input list.

    Args:
        input_list (List): Input list to divide into batches.
        batch_size (int): The number of elements in each batch.

    Yields:
        List: A batch of elements from the input list.
    """
    for i in range(0, len(input_list), batch_size):
        yield input_list[i:i + batch_size]


def get_mgf_content(analysis_dir: str,
                    remove_precursor: bool = False,
                    precursor_peak_width: float = 2.0,
                    batch_size: int = 100,
                    top_n_spectra: Optional[int] = None,
                    min_intensity: float = 0.0,
                    min_charge: Optional[int] = None,
                    max_charge: Optional[int] = None,
                    min_mz: Optional[float] = None,
                    max_mz: Optional[float] = None,
)  -> Generator[Ms2Spectra, None, None]:

    pd_tdf = PandasTdf(str(Path(analysis_dir) / 'analysis.tdf'))
    return get_ms2_dda_content(
        analysis_dir=analysis_dir,
        remove_precursor=remove_precursor,
        precursor_peak_width=precursor_peak_width,
        batch_size=batch_size,
        top_n_spectra=top_n_spectra,
        min_intensity=min_intensity,
        min_charge=min_charge,
        max_charge=max_charge,
        min_mz=min_mz,
        max_mz=max_mz
    )


def get_ms2_dda_content(analysis_dir: str,
                    remove_precursor: bool = False,
                    precursor_peak_width: float = 2.0,
                    batch_size: int = 100,
                    top_n_spectra: Optional[int] = None,
                    min_intensity: float = 0.0,
                    min_charge: Optional[int] = None,
                    max_charge: Optional[int] = None,
                    min_mz: Optional[float] = None,
                    max_mz: Optional[float] = None) -> Generator[Ms2Spectra, None, None]:

    with timsdata.timsdata_connect(analysis_dir) as td:

        analysis_tdf_path = str(Path(analysis_dir) / 'analysis.tdf')
        pd_tdf = PandasTdf(analysis_tdf_path)

        precursors_df = pd_tdf.precursors
        frames_df = pd_tdf.frames

        merged_df = pd.merge(precursors_df,
                             pd_tdf.frames,
                             left_on='Parent',
                             right_on='Id',
                             suffixes=("_Precursor", "_Frame"))

        pasef_frame_msms_info_df = pd_tdf.pasef_frame_msms_info.drop(['Frame'], axis=1)


        # count the number of items in each group
        pasef_frame_msms_info_df['count'] = \
            pasef_frame_msms_info_df.groupby('Precursor')['Precursor'].transform('count')

        # keep only the row for each group
        pasef_frame_msms_info_df = pasef_frame_msms_info_df.drop_duplicates(subset='Precursor', keep='first')
        assert len(pasef_frame_msms_info_df) == len(merged_df)

        merged_df = pd.merge(merged_df,
                             pasef_frame_msms_info_df,
                             left_on='Id_Precursor',
                             right_on='Precursor',
                             suffixes=("_Precursor", "_PasefFrameMsmsInfo")).drop('Precursor', axis=1)

        precursor_to_scan_number = map_precursor_to_ip2_scan_number(precursors_df, frames_df)
        merged_df['IP2ScanNumber'] = merged_df['Id_Precursor'].map(precursor_to_scan_number)
        merged_df.dropna(subset=['MonoisotopicMz', 'Charge'], inplace=True)

        if min_charge is not None:
            merged_df = merged_df[merged_df['Charge'] >= min_charge]
        if max_charge is not None:
            merged_df = merged_df[merged_df['Charge'] <= max_charge]
        if min_mz is not None:
            merged_df = merged_df[merged_df['MonoisotopicMz'] >= min_mz]
        if max_mz is not None:
            merged_df = merged_df[merged_df['MonoisotopicMz'] <= max_mz]

        for precursor_batch in tqdm(
                list(batch_iterator(input_list=list(merged_df.iterrows()), batch_size=batch_size)),
                desc='Generating MS2 Spectra'):
            pasef_ms_ms = None
            pasef_ms_ms = \
                td.readPasefMsMs([int(precursor_row['Id_Precursor']) for _, precursor_row in precursor_batch])

            for _, precursor_row in precursor_batch:

                precursor_id = int(precursor_row['Id_Precursor'])
                parent_id = int(precursor_row['Parent'])
                charge = int(precursor_row['Charge'])

                

                ip2_scan_number = precursor_row['IP2ScanNumber']
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
                ms2_spectra.rt = round(precursor_row['Time'], 2)
                ms2_spectra.ce = round(precursor_row['CollisionEnergy'], 1)
                ms2_spectra.iso_width = round(precursor_row['IsolationWidth'], 1)
                ms2_spectra.iso_mz = round(precursor_row['IsolationMz'], 4)
                ms2_spectra.scan_begin = round(float(precursor_row['ScanNumBegin']), 4)
                ms2_spectra.scan_end = round(float(precursor_row['ScanNumEnd']), 4)
                ms2_spectra.info['Accumulation_Time'] = round(float(precursor_row['AccumulationTime']), 4)
                ms2_spectra.info['Ramp_Time'] = round(float(precursor_row['RampTime']), 4)
                ms2_spectra.info['PASEF_Scans'] = int(precursor_row['count'])

                if 'Pressure' in precursor_row:
                    ms2_spectra.info['Pressure'] = round(float(precursor_row['Pressure']), 4)

                ook0_range = td.scanNumToOneOverK0(int(precursor_row['Id_Frame']),
                                                   [ms2_spectra.scan_begin, ms2_spectra.scan_end])
                ms2_spectra.info['OOK0_Begin'] = round(float(ook0_range[0]), 4)
                ms2_spectra.info['OOK0_End'] = round(float(ook0_range[1]), 4)


                ms2_spectra_data = list(zip(pasef_ms_ms[precursor_id][0], pasef_ms_ms[precursor_id][1]))

                if len(ms2_spectra_data) == 0:
                    continue

                # Convert to numpy arrays for faster operations
                mz_array = np.array([data[0] for data in ms2_spectra_data])
                intensity_array = np.array([data[1] for data in ms2_spectra_data])

                # Apply min_intensity filter
                if min_intensity != 0:
                    intensity_mask = intensity_array >= min_intensity
                    mz_array = mz_array[intensity_mask]
                    intensity_array = intensity_array[intensity_mask]
                    
                    if len(mz_array) == 0:
                        continue

                if remove_precursor:
                    # Remove precursor peak from MS/MS spectra
                    precursor_mz = ms2_spectra.mz
                    min_prec_mz = precursor_mz - precursor_peak_width
                    max_prec_mz = precursor_mz + precursor_peak_width

                    precursor_mask = ~((mz_array >= min_prec_mz) & (mz_array <= max_prec_mz))
                    mz_array = mz_array[precursor_mask]
                    intensity_array = intensity_array[precursor_mask]

                    if len(mz_array) == 0:
                        continue

                if top_n_spectra is not None and len(intensity_array) > top_n_spectra:
                    # Get indices of top N intensities
                    top_indices = np.argpartition(intensity_array, -top_n_spectra)[-top_n_spectra:]
                    mz_array = mz_array[top_indices]
                    intensity_array = intensity_array[top_indices]

                # Sort by m/z values
                sort_indices = np.argsort(mz_array)
                mz_array = mz_array[sort_indices]
                intensity_array = intensity_array[sort_indices]

                # Convert to lists
                ms2_spectra.mz_spectra = mz_array.tolist()
                ms2_spectra.intensity_spectra = intensity_array.astype(int).tolist()

                assert len(ms2_spectra.mz_spectra) == len(ms2_spectra.intensity_spectra)

                if len(ms2_spectra.mz_spectra) == 0:
                    continue

                yield ms2_spectra


def write_mgf_file(analysis_dir: str,
                   output_file: Optional[str] = None,
                    remove_precursor: bool = False,
                    precursor_peak_width: float = 2.0,
                    batch_size: int = 100,
                    top_n_spectra: Optional[int] = None,
                    min_intensity: float = 0.0,
                    min_charge: Optional[int] = None,
                    max_charge: Optional[int] = None,
                    min_mz: Optional[float] = None,
                    max_mz: Optional[float] = None,
):

    start_time = time.time()

    if output_file is None:
        output_file = str(Path(analysis_dir) / Path(analysis_dir).stem) + '.mgf'


    logger.info('Generating Ms2 Spectra')
    ms2_spectra = get_mgf_content(
        analysis_dir=analysis_dir,
        remove_precursor=remove_precursor,
        precursor_peak_width=precursor_peak_width,
        batch_size=batch_size,
        top_n_spectra=top_n_spectra,
        min_intensity=min_intensity,
        min_charge=min_charge,
        max_charge=max_charge,
        min_mz=min_mz,
        max_mz=max_mz
    )


    time.sleep(1)

    logger.info('Writing Contents To File')
    with open(output_file, 'w', encoding='UTF-8') as file:
        for spectrum in tqdm(ms2_spectra, desc='Writing MGF File'):
            mgf_lines = []
            mgf_lines.append("BEGIN IONS")
            mgf_lines.append(f'TITLE={Path(analysis_dir).stem}.{spectrum.low_scan}.{spectrum.high_scan}.{spectrum.charge} ' \
                             f'FILE="{Path(analysis_dir).stem}", NativeID="merged={spectrum.precursor_id} frame={spectrum.parent_id} ' \
                             f'scanStart={spectrum.scan_begin} scanEnd={spectrum.scan_end}"')
            mgf_lines.append(f'RTINSECONDS={spectrum.rt}')
            mgf_lines.append(f'PEPMASS={spectrum.mass}')
            mgf_lines.append(f'CHARGE={spectrum.charge}+')
            for mz, intensity in zip(spectrum.mz_spectra, spectrum.intensity_spectra):
                mgf_lines.append(f'{mz:.5f} {int(intensity)}')
            mgf_lines.append("END IONS")

            file.write('\n'.join(mgf_lines) + '\n')


    total_time = round(time.time() - start_time, 2)
    logger.info(f'Total Time: {total_time:.2f} seconds')


def main():
    """
    Command-line interface for MGF extraction from TimsTOF data.
    """
    parser = argparse.ArgumentParser(
        description="Extract MGF files from TimsTOF .D folders",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        "analysis_dir",
        type=str,
        help="Path to the .D analysis directory"
    )
    
    parser.add_argument(
        "-o", "--output",
        type=str,
        default=None,
        help="Output MGF file path (default: <analysis_dir_name>.mgf)"
    )
    
    parser.add_argument(
        "--remove-precursor",
        action="store_true",
        help="Remove precursor peaks from MS/MS spectra"
    )
    
    parser.add_argument(
        "--precursor-peak-width",
        type=float,
        default=2.0,
        help="Width around precursor m/z to remove (Da)"
    )
    
    parser.add_argument(
        "--batch-size",
        type=int,
        default=100,
        help="Batch size for processing spectra"
    )
    
    parser.add_argument(
        "--top-n-spectra",
        type=int,
        default=None,
        help="Keep only top N most intense peaks per spectrum"
    )
    
    parser.add_argument(
        "--min-intensity",
        type=float,
        default=0.0,
        help="Minimum intensity threshold for peaks"
    )
    
    parser.add_argument(
        "--min-charge",
        type=int,
        default=None,
        help="Minimum charge state filter"
    )
    
    parser.add_argument(
        "--max-charge",
        type=int,
        default=None,
        help="Maximum charge state filter"
    )
    
    parser.add_argument(
        "--min-mz",
        type=float,
        default=None,
        help="Minimum m/z filter"
    )
    
    parser.add_argument(
        "--max-mz",
        type=float,
        default=None,
        help="Maximum m/z filter"
    )
    
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose logging"
    )
    
    args = parser.parse_args()
    
    # Set up logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    
    # Validate input directory
    analysis_path = Path(args.analysis_dir)
    if not analysis_path.exists():
        logger.error(f"Analysis directory does not exist: {args.analysis_dir}")
        return 1
    
    if not analysis_path.is_dir():
        logger.error(f"Path is not a directory: {args.analysis_dir}")
        return 1
    
    # Check for required files
    tdf_file = analysis_path / "analysis.tdf"
    if not tdf_file.exists():
        logger.error(f"Required file not found: {tdf_file}")
        return 1
    
    try:
        write_mgf_file(
            analysis_dir=args.analysis_dir,
            output_file=args.output,
            remove_precursor=args.remove_precursor,
            precursor_peak_width=args.precursor_peak_width,
            batch_size=args.batch_size,
            top_n_spectra=args.top_n_spectra,
            min_intensity=args.min_intensity,
            min_charge=args.min_charge,
            max_charge=args.max_charge,
            min_mz=args.min_mz,
            max_mz=args.max_mz
        )
        logger.info("MGF extraction completed successfully!")
        return 0
    except Exception as e:
        logger.error(f"Error during MGF extraction: {e}")
        return 1


if __name__ == "__main__":
    exit(main())
