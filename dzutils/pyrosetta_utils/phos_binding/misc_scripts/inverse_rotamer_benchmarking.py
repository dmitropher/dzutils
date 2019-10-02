from itertools import product, groupby
import logging

import pyrosetta

import numpy as np
import pandas as pd
import getpy as gp
import scipy.spatial
from scipy.spatial.kdtree import KDTree

from xbin import XformBinner as xb

from dzutils.pyrosetta_utils import (
    run_pyrosetta_with_flags,
    residue_type_from_name3,
    residue_from_name3,
)
from dzutils.pyrosetta_utils.phos_binding import (
    phospho_residue_inverse_rotamer_rts,
)

logger = logging.getLogger("make_inverse_rot_tables")
logger.setLevel(logging.DEBUG)

type_dict = {
    "PTR": ("TYR", pyrosetta.rosetta.core.chemical.VariantType.PHOSPHORYLATION)
}


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
    return pyrosetta.rosetta.core.conformation.Residue(residue_type, False)


def set_rotamer_chis(rot, *chis):
    for i, chi in enumerate(chis, 1):
        rot.set_chi(i, chi)
    return rot


def get_rotamer_pose_from_name(resname, *chis):
    pose = pyrosetta.rosetta.core.pose.Pose()
    residue = get_rotamer_residue_from_name(resname)
    pose.append_residue_by_bond(residue)
    for i, chi in enumerate(chis, 1):
        try:
            pose.set_chi(i, 1, chi)
        except RuntimeError as e:
            logger.debug(f"chinum: {i}")
            logger.debug(f"resnum: {1}")
            logger.debug(f"chi_val: {chi}")
            logger.debug(pose)
            logger.debug(pose.annotated_sequence())
            logger.debug("issue with rotamer")
            raise RuntimeError
    rotamer_pose = pose
    return rotamer_pose


def value_to_pose(
    value,
    data_table,
    conversion_type,
    *data_table_fields,
    residue_type_name=None,
):
    if conversion_type == "from_file":
        return pyrosetta.pose_from_file(
            data_table[data_table["index"] == value][list(data_table_fields)]
        )
    if conversion_type == "from_chis":
        # get chis from the data table somehow
        table_entry = data_table[data_table["index"] == value]
        # print (table_entry)
        # print (value)
        chis = tuple(*table_entry.iloc[0][list(data_table_fields)].values)
        # print (chis)
        try:
            return get_rotamer_pose_from_name(residue_type_name, *chis)
        except RuntimeError as e:
            logger.debug(f"residue type name {residue_type_name}")
            logger.debug(f"chis {chis}")
            logger.debug(f"table entry: {table_entry}")


def rotamer_rmsd(pose_1, pose_2, super=True):
    if super:
        return pyrosetta.rosetta.core.scoring.all_atom_rmsd(pose_1, pose_2)
    else:
        return pyrosetta.rosetta.core.scoring.all_atom_rmsd_nosuper(
            pose_1, pose_2
        )


def process_hits_by_rmsd(
    hits_df,
    values,
    data_table,
    conversion_type,
    *data_table_fields,
    super=True,
):
    """
    Get the rmsd for the hits and attach it to the hits df

    default is to use the super full atom rmsd
    """
    # print(len(hits_df))
    ideal_poses = [
        value_to_pose(
            value,
            data_table,
            conversion_type,
            *data_table_fields,
            residue_type_name=name,
        )
        for value, name in zip(values, hits_df["res_type_name"])
    ]
    pose_pairs = zip(hits_df["rot_pose"].copy(), ideal_poses)
    hits_df["rmsd_to_ideal"] = pd.Series(
        [float(rotamer_rmsd(*pair, super=super)) for pair in pose_pairs],
        dtype=np.float64,
    )

    rmsd_to_ideal_df = hits_df

    return rmsd_to_ideal_df


def kd_tree_from_rotamer_table(chi_table, *fields):
    chis = []
    if len(fields) > 1:
        chis = list(tuple(chis) for chis in chi_table[list(fields)].values)
    else:
        chis = list(tuple(*chis) for chis in chi_table[list(fields)].values)
    chis.sort()
    unique_chis = np.array(list(k for k, _ in groupby(chis)))
    # print(unique_chis)
    return KDTree(unique_chis), unique_chis


def process_misses(misses_df, data_table, *rotamer_fields):
    nearest_chi_tree, chi_data = kd_tree_from_rotamer_table(
        data_table, *rotamer_fields
    )
    columns = [
        f"chi_{i}"
        for i in range(1, len(misses_df.columns) + 1)
        if f"chi_{i}" in misses_df
    ]
    miss_chis = zip(
        misses_df["res_type_name"], *[misses_df[col] for col in columns]
    )
    miss_rmsd = [
        rotamer_rmsd(
            pose,
            get_rotamer_pose_from_name(
                restype_chis[0],
                *chi_data[nearest_chi_tree.query(restype_chis[1:])[1]],
            ),
        )
        for pose, restype_chis in zip(misses_df["rot_pose"], miss_chis)
    ]
    misses_df["rmsd_to_ideal"] = pd.Series(miss_rmsd)
    return misses_df


def test_hash_table(
    table,
    test_data,
    xbin_cart_resl=1,
    xbin_ori_resl=15,
    data_table=None,
    data_table_path="",
    data_table_index_conversion="from_chis",
    data_table_rotamer_fields=("chis",),
    super_hits_before_rmsd=True,
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
    res_type_names, rot_poses, rts, chi_tuples = zip(
        *[
            (
                rot_info[0],
                get_rotamer_pose_from_name(rot_info[0], *rot_info[1:]),
                rt,
                tuple(rot_info[1:]),
            )
            for rot_info in test_data
            for rt in phospho_residue_inverse_rotamer_rts(
                set_rotamer_chis(
                    get_rotamer_residue_from_name(rot_info[0]), *rot_info[1:]
                )
            )
        ]
    )

    binner = xb(xbin_cart_resl, xbin_ori_resl)

    test_keys = binner.get_bin_index(np.array(rts))
    hits_mask = table.contains(test_keys)
    # print (any(hits_mask))
    miss_mask = hits_mask == False
    all_results_df = pd.DataFrame(
        {
            **{
                "res_type_name": res_type_names,
                "rot_pose": rot_poses,
                "rt": rts,
                "key": test_keys,
            },
            **{f"chi_{i}": chis for i, chis in enumerate(zip(*chi_tuples), 1)},
        }
    )
    hits_df = all_results_df[hits_mask].copy()
    hits_keys = hits_df["key"].copy()
    if data_table is None:
        pd.read_json(data_table_path)
    values = pd.Series(table[np.array(hits_keys)], dtype=np.dtype("i8"))
    # print (values)
    misses_df = all_results_df[miss_mask].copy()

    hits = process_hits_by_rmsd(
        hits_df,
        values,
        data_table,
        data_table_index_conversion,
        *data_table_rotamer_fields,
        super=super_hits_before_rmsd,
    )
    misses = process_misses(misses_df, data_table, *data_table_rotamer_fields)
    hits = hits[[col for col in hits if col != "rot_pose"]]
    misses = misses[[col for col in misses if col != "rot_pose"]]
    return hits, misses


def build_test_rotamers_from_resname(residue_type_name, margin, granularity):
    res_type = get_residue_type_from_name(residue_type_name)
    rots = [
        tuple([rot.chi(i) for i in range(1, rot.nchi() + 1)])
        for rot in pyrosetta.rosetta.core.pack.rotamer_set.bb_independent_rotamers(
            res_type
        )
    ]
    return [
        (residue_type_name, chi1, chi2, chi3)
        for start_chi1, start_chi2, start_chi3 in rots
        for chi1, chi2, chi3 in product(
            np.linspace(start_chi1 - margin, start_chi1 + margin, granularity),
            np.linspace(start_chi2 - margin, start_chi2 + margin, granularity),
            np.linspace(start_chi3 - margin, start_chi3 + margin, granularity),
        )
    ]


"/home/dzorine/projects/phos_binding/pilot_runs/loop_grafting/fragment_tables/inverse_rotamer/dicts/benchmarking_1_inv_rot_0.125_ang_2.5_deg.bin"
"/home/dzorine/projects/phos_binding/pilot_runs/loop_grafting/fragment_tables/inverse_rotamer/dicts/benchmarking_1_inv_rot_0.25_ang_5.0_deg.bin"
"/home/dzorine/projects/phos_binding/pilot_runs/loop_grafting/fragment_tables/inverse_rotamer/dicts/benchmarking_1_inv_rot_0.5_ang_10.0_deg.bin"
"/home/dzorine/projects/phos_binding/pilot_runs/loop_grafting/fragment_tables/inverse_rotamer/dicts/benchmarking_1_inv_rot_1.0_ang_15.0_deg.bin"
