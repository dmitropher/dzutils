# Scan for ploop insertion

# Hunt for inverse rotamers to complete interaction

# Takes: pdb, out_dir, out_prefix (uses defaults)
import sys

import pandas as pd

from getpy import Dict as GDict

import numpy as np


from xbin import XformBinner as xb


from dzutils.sutils import read_flag_file
from dzutils.pyrosetta_utils.phos_binding import loops_to_rt_dict

from dzutils.pyrosetta_utils.geometry.pose_xforms import (
    generate_pose_rt_between_res,
)


def load_table_and_dict(table_path, dict_path, key_type, value_type):
    """
    return a tuple of table,dict
    """
    gp_dict = GDict(key_type, value_type)
    gp_dict.load(dict_path)
    table = pd.read_json(table_path)
    return table, gp_dict


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
    # binner = xb()
    pdf["key"] = pdf["rt"].apply(lambda x: xb().get_bin_index(x))

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
    inv_rot_table_path,
    inv_rot_dict_path,
    inv_rot_key_type=np.dtype("i8"),
    inv_rot_value_type=np.dtype("i8"),
):

    results = primary_results_table

    # split each allowed res into its own row
    results = (
        results.allowed_res.apply(pd.Series)
        .merge(results, right_index=True, left_index=True)
        .melt(id_vars=[*results], value_name="p_res_target")
        .drop(["allowed_res", "variable"], axis=1)
    )

    # generate the keys for inverse rotamer
    binner = xb()
    results["inv_rot_key"] = results.apply(
        lambda row: binner.get_bin_index(
            row["loop_func_to_bb_start"]
            @ generate_pose_rt_between_res(
                pose.clone(), row["start_res"], row["p_res_target"]
            )
        ),
        axis=1,
    )

    # load inv_rot table and dict
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
        print("no secondary hits found")
        return

    # retrieve rotamers from table
    rot_masked_df["inv_rot_inds"] = inv_rot_dict[
        np.array(rot_masked_df["inv_rot_key"])
    ]
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
    binner = xb()
    pdf["key"] = pdf["rt"].apply(lambda x: xb().get_bin_index(x))
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
