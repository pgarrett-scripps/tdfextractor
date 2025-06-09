"""
ms2_extractor defines functions for generating ms2 files from DDA and PRM based .D folders
"""

import logging
import time
from datetime import datetime
from pathlib import Path
from typing import List

import pandas as pd
from tdfpy import timsdata
from tdfpy.pandas_tdf import PandasTdf
from serenipy.ms2 import Ms2Spectra, to_ms2
from tqdm import tqdm

from .constants import MS2_VERSION
from .string_templates import MS2_HEADER
from .utils import calculate_mass, map_precursor_to_ip2_scan_number

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
        yield input_list[i : i + batch_size]


def get_ms2_content(
    analysis_dir: str,
    include_spectra: bool = True,
    batch_size: int = 100,
    remove_charge1: bool = True,
    remove_empty_spectra: bool = True,
    min_intensity: float = 0,
):
    """
    Retrieves MS/MS content based on the TDF format
    (Data Dependent Acquisition - DDA or Parallel Reaction Monitoring - PRM).

    Args:
        analysis_dir (str): Path to the analysis directory.
        include_spectra (bool, optional): Whether to include the actual MS/MS spectra.
        batch_size (int, optional): The number of precursors to process in a batch.
        remove_charge1 (bool, optional): Whether to exclude precursors with charge 1.
        remove_empty_spectra (bool, optional): Whether to exclude MS/MS spectra with no spectra values.
        min_intensity (float, optional): Filters out MS/MS spectral data based on intensity (inclusive).

    Returns:
        List[Ms2Spectra]: A list of Ms2Spectra objects representing MS/MS spectra.
        or
        None: If the TDF format is unknown.

    Raises:
        TypeError: If the TDF format is unknown.
    """

    pd_tdf = PandasTdf(str(Path(analysis_dir) / "analysis.tdf"))
    if pd_tdf.is_dda:
        logger.info("TDF format is DDA")
        return get_ms2_dda_content(
            analysis_dir=analysis_dir,
            include_spectra=include_spectra,
            batch_size=batch_size,
            remove_charge1=remove_charge1,
            remove_empty_spectra=remove_empty_spectra,
            min_intensity=min_intensity,
        )
    if pd_tdf.is_prm:
        logger.info("TDF format is PRM")
        return get_ms2_prm_content(
            analysis_dir=analysis_dir,
            include_spectra=include_spectra,
            batch_size=batch_size,
            remove_charge1=remove_charge1,
            remove_empty_spectra=remove_empty_spectra,
            min_intensity=min_intensity,
        )

    raise TypeError("Unknown TDF format")


def get_ms2_dda_content(
    analysis_dir: str,
    include_spectra: bool = True,
    batch_size: int = 100,
    remove_charge1: bool = True,
    remove_empty_spectra: bool = True,
    min_intensity: float = 0,
):
    """
    Retrieves MS/MS spectra from a given analysis directory.

    Args:
        analysis_dir (str): Path to the analysis directory.
        include_spectra (bool, optional): Whether to include the actual MS/MS spectra in the output.
        batch_size (int, optional): The number of precursors to process in a batch.
        remove_charge1 (bool, optional): Whether to exclude precursors with charge 1.
        remove_empty_spectra (bool, optional): Whether to exclude MS/MS spectra with no spectra values.
        min_intensity (float, optional): Filters out MS/MS spectral data based on intensity (inclusive).

    Returns:
        List[Ms2Spectra]: A list of Ms2Spectra objects representing MS/MS spectra.

    Raises:
        FileNotFoundError: If the analysis directory does not exist.
    """
    with timsdata.timsdata_connect(analysis_dir) as td:

        analysis_tdf_path = str(Path(analysis_dir) / "analysis.tdf")
        pd_tdf = PandasTdf(analysis_tdf_path)

        precursors_df = pd_tdf.precursors
        frames_df = pd_tdf.frames

        merged_df = pd.merge(
            precursors_df,
            pd_tdf.frames,
            left_on="Parent",
            right_on="Id",
            suffixes=("_Precursor", "_Frame"),
        )

        pasef_frame_msms_info_df = pd_tdf.pasef_frame_msms_info.drop(["Frame"], axis=1)

        # count the number of items in each group
        pasef_frame_msms_info_df["count"] = pasef_frame_msms_info_df.groupby(
            "Precursor"
        )["Precursor"].transform("count")

        # keep only the row for each group
        pasef_frame_msms_info_df = pasef_frame_msms_info_df.drop_duplicates(
            subset="Precursor", keep="first"
        )
        assert len(pasef_frame_msms_info_df) == len(merged_df)

        merged_df = pd.merge(
            merged_df,
            pasef_frame_msms_info_df,
            left_on="Id_Precursor",
            right_on="Precursor",
            suffixes=("_Precursor", "_PasefFrameMsmsInfo"),
        ).drop("Precursor", axis=1)

        precursor_to_scan_number = map_precursor_to_ip2_scan_number(
            precursors_df, frames_df
        )
        merged_df["IP2ScanNumber"] = merged_df["Id_Precursor"].map(
            precursor_to_scan_number
        )
        merged_df.dropna(subset=["MonoisotopicMz", "Charge"], inplace=True)

        for precursor_batch in tqdm(
            list(
                batch_iterator(
                    input_list=list(merged_df.iterrows()), batch_size=batch_size
                )
            ),
            desc="Generating MS2 Spectra",
        ):
            pasef_ms_ms = None
            if include_spectra:
                pasef_ms_ms = td.readPasefMsMs(
                    [
                        int(precursor_row["Id_Precursor"])
                        for _, precursor_row in precursor_batch
                    ]
                )

            for _, precursor_row in precursor_batch:

                precursor_id = int(precursor_row["Id_Precursor"])
                parent_id = int(precursor_row["Parent"])
                charge = int(precursor_row["Charge"])

                if remove_charge1 is True and charge == 1:
                    continue

                ip2_scan_number = precursor_row["IP2ScanNumber"]
                ook0 = td.scanNumToOneOverK0(parent_id, [precursor_row["ScanNumber"]])[
                    0
                ]
                ccs = timsdata.oneOverK0ToCCSforMz(
                    ook0, charge, precursor_row["MonoisotopicMz"]
                )
                mz = precursor_row["MonoisotopicMz"]
                prec_intensity = precursor_row["Intensity"]
                mass = calculate_mass(mz, charge)

                ms2_spectra = Ms2Spectra(
                    low_scan=ip2_scan_number,
                    high_scan=ip2_scan_number,
                    mz=mz,
                    mass=mass,
                    charge=charge,
                    info={},
                    mz_spectra=[],
                    intensity_spectra=[],
                    charge_spectra=[],
                )

                ms2_spectra.parent_id = parent_id
                ms2_spectra.precursor_id = precursor_id
                ms2_spectra.prec_intensity = int(prec_intensity)
                ms2_spectra.ook0 = round(ook0, 4)
                ms2_spectra.ccs = round(ccs, 1)
                ms2_spectra.rt = round(precursor_row["Time"], 2)
                ms2_spectra.ce = round(precursor_row["CollisionEnergy"], 1)
                ms2_spectra.iso_width = round(precursor_row["IsolationWidth"], 1)
                ms2_spectra.iso_mz = round(precursor_row["IsolationMz"], 4)
                ms2_spectra.scan_begin = round(float(precursor_row["ScanNumBegin"]), 4)
                ms2_spectra.scan_end = round(float(precursor_row["ScanNumEnd"]), 4)
                ms2_spectra.info["Accumulation_Time"] = round(
                    float(precursor_row["AccumulationTime"]), 4
                )
                ms2_spectra.info["Ramp_Time"] = round(
                    float(precursor_row["RampTime"]), 4
                )
                ms2_spectra.info["PASEF_Scans"] = int(precursor_row["count"])

                if "Pressure" in precursor_row:
                    ms2_spectra.info["Pressure"] = round(
                        float(precursor_row["Pressure"]), 4
                    )

                ook0_range = td.scanNumToOneOverK0(
                    int(precursor_row["Id_Frame"]),
                    [ms2_spectra.scan_begin, ms2_spectra.scan_end],
                )
                ms2_spectra.info["OOK0_Begin"] = round(float(ook0_range[0]), 4)
                ms2_spectra.info["OOK0_End"] = round(float(ook0_range[1]), 4)

                if include_spectra:
                    ms2_spectra_data = list(
                        zip(pasef_ms_ms[precursor_id][0], pasef_ms_ms[precursor_id][1])
                    )

                    if min_intensity != 0:
                        ms2_spectra_data = [
                            data
                            for data in ms2_spectra_data
                            if data[1] >= min_intensity
                        ]

                    ms2_spectra.mz_spectra = [data[0] for data in ms2_spectra_data]
                    ms2_spectra.intensity_spectra = [
                        int(data[1]) for data in ms2_spectra_data
                    ]

                    assert len(ms2_spectra.mz_spectra) == len(
                        ms2_spectra.intensity_spectra
                    )

                    if (
                        remove_empty_spectra is True
                        and len(ms2_spectra.mz_spectra) == 0
                    ):
                        continue

                yield ms2_spectra


def get_ms2_prm_content(
    analysis_dir: str,
    include_spectra: bool = True,
    batch_size: int = 100,
    remove_charge1: bool = True,
    remove_empty_spectra: bool = True,
    min_intensity: float = 0,
):
    """
    Retrieves MS/MS spectra from a given PRM analysis directory.

    Args:
        analysis_dir (str): Path to the analysis directory.
        include_spectra (bool, optional): Whether to include the actual MS/MS spectra in the output.
        batch_size (int, optional): Not used with PRM extractor. Retained for consistency.
        remove_charge1 (bool, optional): Whether to exclude precursors with charge 1.
        remove_empty_spectra (bool, optional): Whether to exclude MS/MS spectra with no spectra values.
        min_intensity (float, optional): Filters out MS/MS spectral data based on intensity (inclusive).

    Returns:
        List[Ms2Spectra]: A list of Ms2Spectra objects representing MS/MS spectra.

    Raises:
        FileNotFoundError: If the analysis directory does not exist.
    """
    with timsdata.timsdata_connect(analysis_dir) as td:
        analysis_tdf_path = str(Path(analysis_dir) / "analysis.tdf")
        merged_df = pd.merge(
            PandasTdf(analysis_tdf_path).prm_frame_msms_info,
            PandasTdf(analysis_tdf_path).prm_targets,
            left_on="Target",
            right_on="Id",
        )
        merged_df = pd.merge(
            merged_df,
            PandasTdf(analysis_tdf_path).frames,
            left_on="Frame",
            right_on="Id",
            suffixes=("_Target", "_Frame"),
        )

        for _, row in tqdm(
            merged_df.iterrows(), desc="Generating MS2 Spectra", total=len(merged_df)
        ):

            if remove_charge1 is True and int(row["Charge"]) == 1:
                continue

            mz_list, area_list = td.extractCentroidedSpectrumForFrame(
                frame_id=int(row["Frame"]),
                scan_begin=int(row["ScanNumBegin"]),
                scan_end=int(row["ScanNumEnd"]),
                peak_picker_resolution=120000,
            )

            ms2_spectra = Ms2Spectra(
                low_scan=int(row["Frame"]),
                high_scan=int(row["Frame"]),
                mz=float(row["IsolationMz"]),
                mass=calculate_mass(float(row["IsolationMz"]), int(row["Charge"])),
                charge=int(row["Charge"]),
                info={
                    "Target": str(int(row["Target"])),
                    "Accumulation_Time": str(float(row["AccumulationTime"])),
                    "Ramp_Time": str(float(row["RampTime"])),
                    "Pressure": str(float(row["Pressure"])),
                },
                mz_spectra=[],
                intensity_spectra=[],
                charge_spectra=[],
            )

            ms2_spectra.ce = float(row["CollisionEnergy"])
            ms2_spectra.iso_width = float(row["IsolationWidth"])
            ms2_spectra.iso_mz = float(row["IsolationMz"])
            ms2_spectra.rt = float(row["Time_Frame"])
            ms2_spectra.scan_begin = float(row["ScanNumBegin"])
            ms2_spectra.scan_end = float(row["ScanNumEnd"])

            ook0_range = td.scanNumToOneOverK0(
                int(row["Frame"]), [ms2_spectra.scan_begin, ms2_spectra.scan_end]
            )
            ms2_spectra.info["OOK0_Begin"] = ook0_range[0]
            ms2_spectra.info["OOK0_End"] = ook0_range[1]

            if include_spectra:

                ms2_spectra_data = list(zip(list(mz_list), list(area_list)))

                if min_intensity != 0:
                    ms2_spectra_data = [
                        data for data in ms2_spectra_data if data[1] >= min_intensity
                    ]

                ms2_spectra.mz_spectra = [data[0] for data in ms2_spectra_data]
                ms2_spectra.intensity_spectra = [
                    int(data[1]) for data in ms2_spectra_data
                ]

                assert len(ms2_spectra.mz_spectra) == len(ms2_spectra.intensity_spectra)

                if remove_empty_spectra is True and len(ms2_spectra.mz_spectra) == 0:
                    continue

            yield ms2_spectra


def generate_header(
    analysis_dir: str,
    min_intensity: float,
    remove_charge1: bool,
    remove_empty_spectra: bool,
    include_spectra: bool,
    resolution: float,
):
    """
    Generates a header string for MS2 data using information from the analysis file.

    Args:
        analysis_dir (str): The directory path where the analysis file `analysis.tdf` is located.
        min_intensity (float): The minimum intensity of MS/MS spectra to be considered.
        remove_charge1 (bool): Indicates whether to exclude precursors with a charge of 1.
        remove_empty_spectra (bool): Indicates whether to exclude MS/MS spectra with no spectra values.
        include_spectra (bool): Indicates whether to include the actual MS/MS spectra in the output.
        resolution (float): The resolution of the MS/MS spectra.

    Returns:
        str: The generated header string for MS2 data.

    Raises:
        TypeError: If the TDF format is unknown.
    """

    pd_tdf = PandasTdf(str(Path(analysis_dir) / "analysis.tdf"))

    if pd_tdf.is_dda:
        method = "Data-Dependent"
        precursors_df = pd_tdf.precursors
        frames_df = pd_tdf.frames
        precursor_to_scan_number = map_precursor_to_ip2_scan_number(
            precursors_df, frames_df
        )
        first_scan = list(precursor_to_scan_number.values())[0]
        last_scan = list(precursor_to_scan_number.values())[-1]
    elif pd_tdf.is_prm:
        method = "Parallel-Reaction-Monitoring"
        frames_df = pd_tdf.frames
        frames_df = frames_df[frames_df["MsMsType"] == 10]
        first_scan = min(frames_df["Id"])
        last_scan = max(frames_df["Id"])
    else:
        raise TypeError("Unknown TDF format")

    ms2_header = MS2_HEADER.format(
        version=MS2_VERSION,
        date_of_creation=str(datetime.now().strftime("%B %d, %Y %H:%M")),
        minimum_intensity=min_intensity,
        remove_charge1=remove_charge1,
        remove_empty_spectra=remove_empty_spectra,
        include_spectra=include_spectra,
        method=method,
        resolution=resolution,
        first_scan=first_scan,
        last_scan=last_scan,
    )

    return ms2_header


def write_ms2_file(
    analysis_dir: str,
    output_file: str = None,
    include_spectra: bool = True,
    batch_size: int = 1000,
    remove_charge1: bool = True,
    remove_empty_spectra: bool = True,
    min_intensity: float = 0,
):
    """
    Writes MS/MS spectra data to an MS2 file.

    Args:
        analysis_dir (str): Path to the analysis directory.
        include_spectra (bool, optional): Whether to include the actual MS/MS spectra in the output.
        output_file (str, optional): The output file name.
        batch_size (int, optional): The number of precursors to process in a batch.
        remove_charge1 (bool, optional): Whether to exclude precursors with charge 1.
        remove_empty_spectra (bool, optional): Whether to exclude MS/MS spectra with no spectra values.
        min_intensity (float, optional): Filters out MS/MS spectral data based on intensity (inclusive).

    Returns:
        None.

    Raises:
        FileNotFoundError: If the analysis directory does not exist.
    """
    start_time = time.time()

    if output_file is None:
        output_file = str(Path(analysis_dir) / Path(analysis_dir).stem) + ".ms2"

    logger.info(
        f"Arguments: analysis_dir: {analysis_dir}, output_file: {output_file}, include_spectra: {include_spectra}, "
        f"batch_size: {batch_size}, remove_charge1: {remove_charge1}, "
        f"remove_empty_spectra: {remove_empty_spectra}, min_intensity: {min_intensity}"
    )

    logger.info("Creating Ms2 Header")
    ms2_header = generate_header(
        analysis_dir=analysis_dir,
        include_spectra=include_spectra,
        remove_charge1=remove_charge1,
        remove_empty_spectra=remove_empty_spectra,
        min_intensity=min_intensity,
        resolution=120000,
    )

    logger.info("Generating Ms2 Spectra")
    ms2_spectra = get_ms2_content(
        analysis_dir=analysis_dir,
        include_spectra=include_spectra,
        batch_size=batch_size,
        remove_charge1=remove_charge1,
        remove_empty_spectra=remove_empty_spectra,
        min_intensity=min_intensity,
    )

    logger.info("Creating MS2 Contents")
    ms2_content = to_ms2([ms2_header], ms2_spectra)

    time.sleep(1)  # Dumb fix for logging message overlap with tqdm progress bar

    logger.info("Writing Contents To File")
    with open(output_file, "w", encoding="UTF-8") as file:
        file.write(ms2_content)

    total_time = round(time.time() - start_time, 2)
    logger.info(f"Total Time: {total_time:.2f} seconds")
