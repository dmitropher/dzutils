import pyrosetta

import numpy as np
import pandas as pd

from xbin import XformBinner as xb

from dzutils.pyrosetta_utils import (
    run_pyrosetta_with_flags,
    residue_type_from_name3,
    residue_from_name3,
)
from dzutils.pyrosetta_utils.phos_binding import (
    phospho_residue_inverse_rotamer_rts,
)

type_dict = {}


def get_residue_type_from_name(resname):
    try:
        name, variant = type_dict[resname]
        return residue_type_from_name3(name, variant=variant)
    except KeyError as e:
        print("")
        print(e)
        raise


def get_rotamer_residue_from_name(resname):
    residue_type = get_residue_type_from_name(resname)
    return pyrosetta.rosetta.core.conformation.Residue(residue_type)


def get_rotamer_pose(resname, *chis):
    pose = pyrosetta.rosetta.core.pose.Pose()
    residue = get_rotamer_residue_from_name(resname)
    pose.append_residue_by_bond(residue)
    for i, chi in enumerate(chis, 1):
        pose.set_chi(i, 1, chi)
    rotamer_pose = pose
    return rotamer_pose


def value_to_pose(value, data_table,conversion_type,*data_table_fields,residue_type_name=None):
    if conversion_type == "from_file":
        return pyrosetta.pose_from_file( data_table[value][*data_table_fields])
    if conversion_type == "from_chis":
        #get chis from the data table somehow
        chis = get_chis(value,data_table)
        return get_rotamer_pose(residue_type_name,*chis)



def rotamer_rmsd(pose_1, pose_2, super=True):
    if super:
        return pyrosetta.rosetta.core.scoring.all_atom_rmsd(pose_1, pose_2)
    else:
        return pyrosetta.rosetta.core.scoring.all_atom_rmsd_nosuper(pose_1, pose_2)


def process_hits_by_rmsd(hits_df, **kwargs):
    """
    Get the rmsd for the hits and attach it to the hits df

    default is to use the super full atom rmsd
    """
    ideal_poses = [
        value_to_pose(value, **kwargs) for value in hits_df["value"]
    ]
    pose_pairs = zip(hits_df["rot_pose"].copy(), ideal_poses)
    hits_df.loc["rmsd_to_ideal"] = pd.Series(
        [
            rotamer_rmsd(*pair, kwargs["super"] if "super" in kwargs else True)
            for pair in pose_pairs
        ]
    )

    rmsd_to_ideal_df = hits_df

    return rmsd_to_ideal_df


def process_misses(misses_df, **kwargs):
    return


def test_hash_table(
    table, test_data, xbin_cart_resl=1, xbin_ori_resl=15, **kwargs
):
    """
    Try the test_data rotamers on the table, return hits and misses with info

    test data should be a list of tuples: (resname,chi1,2,...,chin) for the
    number of chis

    hits should have chis for test data and table value, key in the table,
    and rmsd to hit

    misses should have rmsd to the closest "ideal" query, or if not provided
    just the chis and the key
    """
    res_type_names,rot_poses, rts, chi_tuples = zip(
        *[
            (
                rot_info[0],
                get_rotamer_pose(rot_info[0], *rot_info[1:]),
                rt,
                tuple(rot_info[1:]),
            )
            for rot_info in test_data
            for rt in phospho_residue_inverse_rotamer_rts(
                get_rotamer_residue_from_name(rot_info[0])
            )
        ]
    )

    binner = xb(xbin_cart_resl, xbin_ori_resl)

    test_keys = binner.get_bin_index(np.array(rts))
    hits_mask = table.contains(test_keys)
    miss_mask = hits_mask == False
    all_results_df = pd.DataFrame(
        {
            **{"res_type_name":res_type_names},
            **{"rot_pose": rot_poses},
            **{"rt": rts},
            **{f"chi_{i}": chis for i, chis in enumerate(zip(*chi_tuples), 1)},
            **{"key": test_keys},
        }
    )
    hits_df = all_results_df[hits_mask].copy()
    hits_keys = np.array(hits_df["key"].copy())
    hits_df.loc["value"] = table[hits_keys]
    misses_df = all_results_df[miss_mask].copy()
    hits = process_hits_by_rmsd(hits_df, **kwargs)
    misses = process_misses(misses_df, **kwargs)

    return hits, misses
