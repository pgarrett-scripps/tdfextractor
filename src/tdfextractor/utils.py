"""
Utility based functions for ms2 extractor
"""

from pathlib import Path
from typing import Dict, Generator, List, Optional, Tuple

import numpy as np
import pandas as pd
from tdfpy import timsdata
from tdfpy.pandas_tdf import PandasTdf
from serenipy.ms2 import Ms2Spectra
from tqdm import tqdm
from .constants import PROTON_MASS


def map_frame_id_to_ms1_scan(
    parent_precursor_map: Dict[int, List[int]], ms1_ids: List[int]
) -> Tuple[Dict[int, int], Dict[int, Dict[int, int]]]:
    """
    Maps each frame ID to its corresponding MS1 scan and MS2 scans.

    Args:
        parent_precursor_map (dict): A dictionary mapping parent frame IDs to their precursor IDs.
        ms1_ids (list): A list of MS1 frame IDs.

    Returns:
        tuple: Two dictionaries - one mapping frame IDs to MS1 scans, and another mapping frame IDs and precursor
        IDs to their scan counts.
    """
    frame_id_ms1_scan_map = {}
    ms2_map = {}
    prev_scan = 0
    for frame_id in ms1_ids:
        frame_id = int(frame_id)
        prev_scan += 1
        frame_id_ms1_scan_map[frame_id] = prev_scan
        if frame_id in parent_precursor_map:
            if frame_id not in ms2_map:
                ms2_map[frame_id] = {}
            for count, prec_id in enumerate(
                parent_precursor_map[frame_id], prev_scan + 1
            ):
                ms2_map[frame_id][prec_id] = count
            prev_scan += len(parent_precursor_map[frame_id])
    return frame_id_ms1_scan_map, ms2_map


def map_parent_id_to_precursors(precursors_df: pd.DataFrame) -> Dict[int, List[int]]:
    """
    Creates a mapping from parent ID to precursor IDs.

    Args:
        precursors_df (DataFrame): DataFrame of precursors, expected to have 'Parent' and 'Id' columns.

    Returns:
        dict: A dictionary mapping parent IDs to precursor IDs.
    """
    parent_grps = precursors_df.groupby(by="Parent")
    return {parent_id: grp["Id"] for parent_id, grp in parent_grps}


def get_ms1_frames_ids(frames_df: pd.DataFrame) -> np.array:
    """
    Gets MS1 frame IDs.

    Args:
        frames_df (DataFrame): DataFrame of frames, expected to have 'MsMsType' and 'Id' columns.

    Returns:
        array: An array of MS1 frame IDs.
    """
    frames_df = frames_df[frames_df["MsMsType"] == 0]
    return frames_df["Id"].values


def get_ms2_frames_ids(frames_df: pd.DataFrame) -> np.array:
    """
    Gets MS2 frame IDs.

    Args:
        frames_df (DataFrame): DataFrame of frames, expected to have 'MsMsType' and 'Id' columns.

    Returns:
        array: An array of MS2 frame IDs.
    """
    frames_df = frames_df[frames_df["MsMsType"] == 8]
    return frames_df["Id"].values


def map_precursor_to_ip2_scan_number(
    precursors_df: pd.DataFrame, frames_df: pd.DataFrame
) -> Dict[int, int]:
    """
    Maps each precursor ID to its IP2 scan number.

    Args:
        precursors_df (DataFrame): DataFrame of precursors, expected to have 'Parent' and 'Id' columns.
        frames_df (DataFrame): DataFrame of frames, expected to have 'MsMsType' and 'Id' columns.

    Returns:
        dict: A dictionary mapping precursor IDs to their IP2 scan numbers.
    """
    precursor_map = map_parent_id_to_precursors(precursors_df)
    all_ms1_list = get_ms1_frames_ids(frames_df)
    _, ms2_map = map_frame_id_to_ms1_scan(precursor_map, all_ms1_list)
    return {
        prec_id: ms2_map[parent_id][prec_id]
        for parent_id in ms2_map
        for prec_id in ms2_map[parent_id]
    }


def calculate_mass(mz: float, charge: int) -> float:
    """
    Calculates the mass of an ion given its m/z ratio and charge.

    Args:
        mz (float): m/z ratio of ion
        charge (int): Charge of ion

    Returns:
        float: Calculated mass
    """
    return (mz * charge) - (charge - 1) * PROTON_MASS


def batch_iterator(input_list: List, batch_size: int):
    for i in range(0, len(input_list), batch_size):
        yield input_list[i : i + batch_size]


def get_ms2_dda_content(
    analysis_dir: str,
    remove_precursor: bool = False,
    precursor_peak_width: float = 2.0,
    batch_size: int = 100,
    top_n_spectra: Optional[int] = None,
    min_intensity: float = 0.0,
    min_charge: Optional[int] = None,
    max_charge: Optional[int] = None,
    min_mz: Optional[float] = None,
    max_mz: Optional[float] = None,
    min_rt: Optional[float] = None,
    max_rt: Optional[float] = None,
    min_ccs: Optional[float] = None,
    max_ccs: Optional[float] = None,
) -> Generator[Ms2Spectra, None, None]:

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

        if min_charge is not None:
            merged_df = merged_df[merged_df["Charge"] >= min_charge]
        if max_charge is not None:
            merged_df = merged_df[merged_df["Charge"] <= max_charge]
        if min_mz is not None:
            merged_df = merged_df[merged_df["MonoisotopicMz"] >= min_mz]
        if max_mz is not None:
            merged_df = merged_df[merged_df["MonoisotopicMz"] <= max_mz]
        if min_rt is not None:
            merged_df = merged_df[merged_df["Time"] >= min_rt]
        if max_rt is not None:
            merged_df = merged_df[merged_df["Time"] <= max_rt]

        for precursor_batch in tqdm(
            list(
                batch_iterator(
                    input_list=list(merged_df.iterrows()), batch_size=batch_size
                )
            ),
            desc="Generating MS2 Spectra",
        ):
            pasef_ms_ms = None
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

                ip2_scan_number = precursor_row["IP2ScanNumber"]
                ook0 = td.scanNumToOneOverK0(parent_id, [precursor_row["ScanNumber"]])[
                    0
                ]
                ccs = timsdata.oneOverK0ToCCSforMz(
                    ook0, charge, precursor_row["MonoisotopicMz"]
                )

                if min_ccs is not None and ccs < min_ccs:
                    continue
                if max_ccs is not None and ccs > max_ccs:
                    continue

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

                ms2_spectra_data = list(
                    zip(pasef_ms_ms[precursor_id][0], pasef_ms_ms[precursor_id][1])
                )

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

                    precursor_mask = ~(
                        (mz_array >= min_prec_mz) & (mz_array <= max_prec_mz)
                    )
                    mz_array = mz_array[precursor_mask]
                    intensity_array = intensity_array[precursor_mask]

                    if len(mz_array) == 0:
                        continue

                if top_n_spectra is not None and len(intensity_array) > top_n_spectra:
                    # Get indices of top N intensities
                    top_indices = np.argpartition(intensity_array, -top_n_spectra)[
                        -top_n_spectra:
                    ]
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
