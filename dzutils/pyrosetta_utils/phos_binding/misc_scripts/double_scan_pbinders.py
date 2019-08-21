# Scan for ploop insertion

# Hunt for inverse rotamers to complete interaction

# Takes: pdb, out_dir, out_prefix (uses defaults)
import sys

import pandas as pd

from getpy import Dict as GDict

import numpy as np


from xbin import XformBinner as XB


from dzutils.sutils import read_flag_file
from dzutils.pyrosetta_utils.phos_binding import loops_to_rt_dict
from dzutils.pyrosetta_utils.secstruct.structure import (
    parse_structure_from_dssp,
)

from dzutils.pyrosetta_utils.geometry.pose_xforms import (
    generate_pose_rt_between_res,
)


def pairs_in_range(
    begin,
    end,
    min_size,
    max_size,
    max_dist_from_begin=-1,
    max_dist_from_end=-1,
):
    """
    gives all value pairs between begin/end, with appropriate size and offset
    """
    max_begin_position = (
        end - min_size
        if max_dist_from_begin == -1
        else begin + max_dist_from_begin
    )
    min_end_position = (
        begin + min_size
        if max_dist_from_end == -1
        else end - max_dist_from_end
    )
    pairs = [
        (start, start + val)
        for start in range(begin, max_begin_position + 1)
        for val in range(min_size, max_size + 1)
        if start + val >= min_end_position and start + val <= end
    ]
    # print(pairs)
    return pairs


def bb_hash_res_pair(pose, res_1, res_2, *xbin_args, **xbin_kwargs):
    """
    """
    binner = XB(*xbin_args, **xbin_kwargs)
    return binner.get_bin_index(
        generate_pose_rt_between_res(pose, res_1, res_2)
    )


def pose_fragments_to_pose_xforms(
    pose,
    dssp_types,
    min_size,
    max_size,
    *xbin_args,
    plus=0,
    minus=0,
    max_dist_from_begin=-1,
    max_dist_from_end=-1,
):
    """
    Gives the pose_xform objects for all res pairs in fragments from dssp_types

    Must specify min and max size for residues between beginning and end of each
    RT.

    Optional inputs:

    plus or minus range to check around each fragment found by dssp (useful for
    finding residues anchoring loops)

    max distance from beginning and end where a fragment can lie (default any)
    """
    return [
        generate_pose_rt_between_res(pose, i, j)
        for frag in parse_structure_from_dssp(pose, dssp_types)
        for i, j in pairs_in_range(
            max(1, frag.start_pos - minus),
            min(len(pose.residues), frag.end_pos - plus),
            min_size,
            max_size,
            max_dist_from_begin=max_dist_from_begin,
            max_dist_from_end=max_dist_from_end,
        )
    ]


def load_table_and_dict(table_path, dict_path, key_type, value_type):
    """
    return a tuple of table,dict
    """
    gp_dict = GDict(key_type, value_type)
    gp_dict.load(dict_path)
    table = pd.read_json(table_path)
    return table, gp_dict


def scan_for_n_term_helical_grafts(
    pose,
    *xbin_args,
    frag_table=None,
    frag_dict=None,
    frag_table_path=None,
    frag_dict_path=None,
    frag_key_type=np.dtype("i8"),
    frag_value_type=np.dtype("i8"),
    fragment_size_min=4,
    fragment_size_max=10,
    allowed_offset=4,
    **xbin_kwargs,
):
    """

    """
    assert bool(
        (frag_table is not None and frag_dict is not None)
        ^ (frag_table_path is not None and frag_dict_path is not None)
    ), "Must specify either a table/dict path or preload table/dict"

    pose_xforms = pose_fragments_to_pose_xforms(
        pose,
        "H",
        fragment_size_min,
        fragment_size_max,
        *xbin_args,
        plus=0,
        minus=0,
        max_dist_from_begin=allowed_offset,
        max_dist_from_end=-1,
    )
    binner = XB(*xbin_args, **xbin_kwargs)
    xform_key = zip(pose_xforms, binner.get_bin_index(np.array(pose_xforms)))
    rt_dicts = (
        {
            "rt": np.array(xform),
            "name": xform.pose_start.pdb_info().name(),
            "start_res": xform.seqpos_start,
            "end_res": xform.seqpos_end,
            "key": key,
        }
        for xform, key in xform_key
    )
    if not rt_dicts:
        print("empty rt dicts")
        return
    pdf = pd.DataFrame(rt_dicts)
    if frag_table is None and frag_dict is None:
        frag_table, frag_dict = load_table_and_dict(
            frag_table_path, frag_dict_path, frag_key_type, frag_value_type
        )
    mask = frag_dict.contains(np.array(pdf["key"]))

    masked_df = pdf[mask].copy()
    if len(masked_df.index) == 0:
        print("no fragment grafts found")
        return

    masked_df.loc[:, "e2e_inds"] = frag_dict[np.array(masked_df["key"])]
    results = masked_df.loc[~masked_df.index.duplicated(keep="first")]

    for col in frag_table:
        col_name = f"frag_{col}"
        results[col_name] = results["e2e_inds"].map(
            lambda index: frag_table[frag_table["index"] == index][col].item()
        )
    return results


def scan_for_ploop_graft(
    pose,
    loop_table_path,
    loop_dict_path,
    loop_key_type=np.dtype("i8"),
    loop_value_type=np.dtype("i8"),
):
    """
    Scans the pose for a phos loop grafted

    Weird gotchas:
    This function will ignore any grafts where the start residue is residue 1 for
    grafting reasons later

    It will only find grafts where the posenum of the start is smaller than end

    It will only search within each loop, not between loops
    """
    pose = pose.clone()
    # Convert the pose loops into end to end xforms
    rt_dicts = loops_to_rt_dict(pose)
    if not rt_dicts:
        print("empty rt dicts")
        return
    pdf = pd.DataFrame(rt_dicts)
    # binner = XB()
    pdf["key"] = pdf["rt"].apply(lambda x: XB().get_bin_index(x))

    # FIXME
    # allowed res is all res in this case
    pdf["allowed_res"] = pdf["rt"].apply(
        lambda x: [*range(1, len(pose.residues))]
    )

    loop_table, loop_dict = load_table_and_dict(
        loop_table_path, loop_dict_path, loop_key_type, loop_value_type
    )

    mask = loop_dict.contains(np.array(pdf["key"]))
    masked_df = pdf[mask]
    if len(masked_df.index) == 0:
        print("no primary hits found")
        return

    # Hits found, append appropriate fields from data table and generate
    # Inverse rotamer bin key
    masked_df.loc[:, "e2e_inds"] = loop_dict[np.array(masked_df["key"])]
    results = masked_df.loc[~masked_df.index.duplicated(keep="first")]

    for col in loop_table:
        results[f"loop_{col}"] = results["e2e_inds"].map(
            lambda index: loop_table[loop_table["index"] == index][col].item()
        )
    return results


def scan_for_inv_rot(
    pose,
    primary_results_table,
    inv_rot_table=None,
    inv_rot_dict=None,
    inv_rot_table_path=None,
    inv_rot_dict_path=None,
    feature_xform_label="loop_func_to_bb_start",
    inv_rot_key_type=np.dtype("i8"),
    inv_rot_value_type=np.dtype("i8"),
    **xbin_args,
):
    assert bool(
        (inv_rot_table is not None and inv_rot_dict is not None)
        ^ (inv_rot_table_path is not None and inv_rot_dict_path is not None)
    ), "Must specify either a table/dict path or preload table/dict"
    results = primary_results_table

    # split each allowed res into its own row
    results = (
        results.allowed_res.apply(pd.Series)
        .merge(results, right_index=True, left_index=True)
        .melt(id_vars=[*results], value_name="p_res_target")
        .drop(["allowed_res", "variable"], axis=1)
    )

    # generate the keys for inverse rotamer
    binner = XB(**xbin_args)
    working = pose.clone()
    start_to_end_xforms = np.array(
        [
            generate_pose_rt_between_res(working, start, target)
            for start, target in zip(
                results["start_res"].to_list(),
                results["p_res_target"].to_list(),
            )
        ]
    )
    feature_to_start_xforms = np.array(results[feature_xform_label].to_list())
    # print (feature_to_start_xforms,start_to_end_xforms,sep="\n")
    keys = binner.get_bin_index(
        np.array(feature_to_start_xforms @ start_to_end_xforms)
    )
    key_series = pd.Series(keys)
    results["inv_rot_key"] = key_series
    # print (results)
    # print (key_series)
    # print (results.apply(
    #     lambda row: binner.get_bin_index(
    #         row[feature_xform_label]
    #         @ generate_pose_rt_between_res(
    #             pose.clone(), row["start_res"], row["p_res_target"]
    #         )
    #     ),
    #     axis=1,
    # )
    # )
    # print (results)
    # TODO - include the inv rot reference atoms as well as loop func to bb start

    # load inv_rot table and dict
    if inv_rot_table is None and inv_rot_dict is None:
        inv_rot_table, inv_rot_dict = load_table_and_dict(
            inv_rot_table_path,
            inv_rot_dict_path,
            inv_rot_key_type,
            inv_rot_value_type,
        )

    # Mask and discard no hit table
    rot_mask = inv_rot_dict.contains(np.array(results["inv_rot_key"]))
    rot_masked_df = results[rot_mask]
    if len(rot_masked_df.index) == 0:
        return

    # retrieve rotamers from table
    keys = inv_rot_dict[np.array(rot_masked_df["inv_rot_key"].copy())]
    key_series = pd.Series(keys)
    rot_masked_df = rot_masked_df.assign(inv_rot_inds=key_series.values)
    results = rot_masked_df.loc[~rot_masked_df.index.duplicated(keep="first")]
    for col in inv_rot_table:
        results[f"inv_rot_{col}"] = results["inv_rot_inds"].map(
            lambda index: inv_rot_table[inv_rot_table["index"] == index][
                col
            ].item()
        )
    return results


def main():
    from pyrosetta import init, pose_from_file

    ploop_flags_file = "/home/dzorine/phos_binding/pilot_runs/loop_grafting/initial_testing/misc_files/p_ligand.flags"
    flags = read_flag_file(ploop_flags_file)
    flags_str = " ".join(flags.replace("\n", " ").split())
    init(flags_str)

    pdb_path = sys.argv[1]
    pose = pose_from_file(pdb_path)

    rt_dicts = loops_to_rt_dict(pose)
    if not rt_dicts:
        print("empty rt dicts")
        return
    pdf = pd.DataFrame(rt_dicts)
    binner = XB()
    pdf["key"] = pdf["rt"].apply(lambda x: XB().get_bin_index(x))
    # allowed res is all res in this case, why not
    pdf["allowed_res"] = pdf["rt"].apply(
        lambda x: [*range(1, len(pose.residues))]
    )

    # The tables are loaded here from hardcoded paths
    key_type, value_type = np.dtype("i8"), np.dtype("i8")
    data_name = "exp_no_spin_1ang_15k_ploops_v3"
    data_path = "/home/dzorine/phos_binding/pilot_runs/loop_grafting/fragment_tables/ploops_expanded_set_1"
    table_dir_path = f"{data_path}/tables"
    dict_dir_path = f"{data_path}/dicts"
    dict_path = f"{dict_dir_path}/{data_name}.bin"
    table_path = f"{table_dir_path}/{data_name}.json"
    loop_table, loop_dict = load_table_and_dict(
        table_path, dict_path, key_type, value_type
    )
    mask = loop_dict.contains(np.array(pdf["key"]))
    masked_df = pdf[mask]
    if len(masked_df.index) == 0:
        print("no primary hits found")
        return
    masked_df.loc[:, "e2e_inds"] = loop_dict[np.array(masked_df["key"])]
    results = masked_df.loc[~masked_df.index.duplicated(keep="first")]

    for col in loop_table:
        results[f"loop_{col}"] = results["e2e_inds"].map(
            lambda index: loop_table[loop_table["index"] == index][col].item()
        )

    results = (
        results.allowed_res.apply(pd.Series)
        .merge(results, right_index=True, left_index=True)
        .melt(id_vars=[*results], value_name="p_res_target")
        .drop(["allowed_res", "variable"], axis=1)
    )

    results["inv_rot_key"] = results.apply(
        lambda row: binner.get_bin_index(
            row["loop_func_to_bb_start"]
            @ generate_pose_rt_between_res(
                pose.clone(), row["start_res"], row["p_res_target"]
            )
        ),
        axis=1,
    )
    inv_rot_data_path = "/home/dzorine/projects/phos_binding/pilot_runs/loop_grafting/fragment_tables/inverse_rotamer"
    inv_rot_data_name = "inverse_ptr_exchi7_rotamers_1ang_v1"
    inv_rot_dict_path = f"{inv_rot_data_path}/dicts/{inv_rot_data_name}.bin"
    inv_rot_table_path = f"{inv_rot_data_path}/tables/{inv_rot_data_name}.json"
    key_type, value_type = np.dtype("i8"), np.dtype("i8")
    inv_rot_table, inv_rot_dict = load_table_and_dict(
        inv_rot_table_path, inv_rot_dict_path, key_type, value_type
    )
    outpath = sys.argv[2]

    pdb_name = pose.pdb_info().name().split("/")[-1].split(".pdb")[0]
    # results.name = "double_ploop_hits"

    primary_write_path = f"{outpath}/primary_double_ploop_hits_{pdb_name}.json"
    print(primary_write_path)
    results.to_json(primary_write_path)
    rot_mask = inv_rot_dict.contains(np.array(results["inv_rot_key"]))
    rot_masked_df = results[rot_mask]
    if len(rot_masked_df.index) == 0:
        print("no secondary hits found")
        return

    else:
        print("rot_masked", rot_masked_df)
        rot_masked_df["inv_rot_inds"] = inv_rot_dict[
            np.array(rot_masked_df["inv_rot_key"])
        ]
        secondary_results = rot_masked_df.loc[
            ~rot_masked_df.index.duplicated(keep="first")
        ]
        print("post-lookup", secondary_results)
        for col in inv_rot_table:
            secondary_results[f"inv_rot_{col}"] = secondary_results[
                "inv_rot_inds"
            ].map(
                lambda index: inv_rot_table[inv_rot_table["index"] == index][
                    col
                ].item()
            )
        print(secondary_results)
        # secondary_results.name = "double_ploop_hits"
        secondary_write_path = (
            f"{outpath}/secondary_double_ploop_hits_{pdb_name}.json"
        )

        secondary_results.to_json(secondary_write_path)


if __name__ == "__main__":
    main()
