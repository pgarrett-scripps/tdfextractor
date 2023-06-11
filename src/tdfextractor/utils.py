"""
Utility based functions for ms2 extractor
"""

from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

from .constants import PROTON_MASS


def map_frame_id_to_ms1_scan(parent_precursor_map: Dict[int, List[int]], ms1_ids: List[int]) -> \
        Tuple[Dict[int, int], Dict[int, Dict[int, int]]]:
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
            for count, prec_id in enumerate(parent_precursor_map[frame_id], prev_scan + 1):
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
    parent_grps = precursors_df.groupby(by='Parent')
    return {parent_id: grp['Id'] for parent_id, grp in parent_grps}


def get_ms1_frames_ids(frames_df: pd.DataFrame) -> np.array:
    """
    Gets MS1 frame IDs.

    Args:
        frames_df (DataFrame): DataFrame of frames, expected to have 'MsMsType' and 'Id' columns.

    Returns:
        array: An array of MS1 frame IDs.
    """
    frames_df = frames_df[frames_df['MsMsType'] == 0]
    return frames_df['Id'].values


def get_ms2_frames_ids(frames_df: pd.DataFrame) -> np.array:
    """
    Gets MS2 frame IDs.

    Args:
        frames_df (DataFrame): DataFrame of frames, expected to have 'MsMsType' and 'Id' columns.

    Returns:
        array: An array of MS2 frame IDs.
    """
    frames_df = frames_df[frames_df['MsMsType'] == 8]
    return frames_df['Id'].values


def map_precursor_to_ip2_scan_number(precursors_df: pd.DataFrame, frames_df: pd.DataFrame) -> Dict[int, int]:
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
    return {prec_id: ms2_map[parent_id][prec_id] for parent_id in ms2_map for prec_id in ms2_map[parent_id]}


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
