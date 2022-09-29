from .constants import PROTON_MASS

def map_frame_id_to_ms1_scan(parent_precursor_map, ms1_ids):
    frame_id_ms1_scan_map = {}
    ms2_map = {}
    prev_scan = 0
    for id in ms1_ids:
        frame_id = int(id)
        prev_scan += 1
        frame_id_ms1_scan_map[frame_id] = prev_scan
        if frame_id in parent_precursor_map:
            if frame_id not in ms2_map:
                ms2_map[frame_id] = {}
            for count, prec_id in enumerate(parent_precursor_map[frame_id], prev_scan + 1):
                ms2_map[frame_id][prec_id] = count
            prev_scan += len(parent_precursor_map[frame_id])
    return frame_id_ms1_scan_map, ms2_map


def map_parent_id_to_precursors(precursors_df):
    parent_grps = precursors_df.groupby(by='Parent')
    return {parent_id: grp['Id'] for parent_id, grp in parent_grps}


def get_ms1_frames_ids(frames_df):
    frames_df = frames_df[frames_df['MsMsType'] == 0]
    return frames_df['Id'].values


def get_ms2_frames_ids(frames_df):
    frames_df = frames_df[frames_df['MsMsType'] == 8]
    return frames_df['Id'].values


def map_precursor_to_ip2_scan_number(precursors_df, frames_df):
    precursor_map = map_parent_id_to_precursors(precursors_df)
    all_ms1_list = get_ms1_frames_ids(frames_df)
    frame_id_ms1_scan_map, ms2_map = map_frame_id_to_ms1_scan(precursor_map, all_ms1_list)
    return {prec_id:ms2_map[parent_id][prec_id] for parent_id in ms2_map for prec_id in ms2_map[parent_id]}

def calculate_mass(mz, charge):
    return (mz * charge) - (charge - 1) * PROTON_MASS