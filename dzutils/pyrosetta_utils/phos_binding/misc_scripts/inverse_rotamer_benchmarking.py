from itertools import product, groupby
import logging

import pyrosetta

import numpy as np
import pandas as pd
import h5py

# import getpy as gp
# import scipy.spatial
from scipy.spatial.kdtree import KDTree

from xbin import XformBinner as xb

from dzutils.pyrosetta_utils import residue_type_from_name3

from dzutils.pyrosetta_utils.geometry.superposition_utilities import (
    superposition_pose,
)

from dzutils.pyrosetta_utils.phos_binding import (
    phospho_residue_inverse_rotamer_rts,
    pres_bases,
)

from dzutils.pyrosetta_utils.geometry.pose_xforms import RotamerRTArray

logger = logging.getLogger("make_inverse_rot_tables")
logger.setLevel(logging.DEBUG)

type_dict = {
    "PTR": (
        "TYR",
        pyrosetta.rosetta.core.chemical.VariantType.PHOSPHORYLATION,
    ),
    "PHY": ("PHY", None),
}


def get_residue_type_from_name(resname):
    """
    Helper func for getting resdiue type from global map
    """
    try:
        name, variant = type_dict[resname]
        return residue_type_from_name3(name, variant=variant)
    except KeyError as e:
        print("")
        print(e)
        raise


def get_rotamer_residue_from_name(resname):
    """
    Helper func to build rosetta residue objects
    """
    residue_type = get_residue_type_from_name(resname)
    return pyrosetta.rosetta.core.conformation.Residue(residue_type, False)


base_atom_dict = {
    "PTR": pres_bases(get_rotamer_residue_from_name("PTR")),
    "PHY": [
        [
            get_rotamer_residue_from_name("PHY").atom_index(name)
            for name in ["P", "OH", "P", "O2P"]
        ]
    ],
}


def get_base_atoms_from_name(resname):
    """
    helper func for accessing the global type dict
    """
    try:
        atom_list = base_atom_dict[resname]
        return atom_list
    except KeyError as e:
        print("")
        print(e)
        raise


def residue_inverse_rotamer_rts(rotamer_rt_array, base_atoms_list):
    """
    """
    output = []
    for base in base_atoms_list:
        rotamer_rt_array.set_target_atoms(base)
        output.append(np.array(rotamer_rt_array))
    return output


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


def get_rotamer_pose_from_residue_type(residue_type, *chis):
    pose = pyrosetta.rosetta.core.pose.Pose()
    residue = pyrosetta.rosetta.core.conformation.Residue(residue_type, False)
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


def atom_coords(pose):
    coords = pyrosetta.rosetta.utility.vector1_numeric_xyzVector_double_t()
    for residue in pose.residues:
        for atom in residue.atoms():
            coords.append(atom.xyz())
    return coords


def rotamer_rmsd(pose_1, pose_2, super=True):
    if super:

        superposition_pose(pose_1, atom_coords(pose_1), atom_coords(pose_2))

        return pyrosetta.rosetta.core.scoring.all_atom_rmsd_nosuper(
            pose_1, pose_2
        )
    else:
        return pyrosetta.rosetta.core.scoring.all_atom_rmsd_nosuper(
            pose_1, pose_2
        )


def process_hits_by_rmsd(restype, ideal_chis, query_chis, super=True):
    """
    Get the rmsd for query chis to ideal chis given the residue type

    returns a list of rmsd values
    """
    ideal_res = pyrosetta.rosetta.core.conformation.Residue(restype, False)
    query_res = pyrosetta.rosetta.core.conformation.Residue(restype, False)
    # ideal_res_pose.append_residue_by_bond(ideal_res)

    # query_res_pose.append_residue_by_bond(query_res)

    rmsd_list = []
    for ideal, query in zip(ideal_chis, query_chis):
        rot_i = set_rotamer_chis(ideal_res, *ideal)
        rot_q = set_rotamer_chis(query_res, *query)
        query_res_pose = pyrosetta.rosetta.core.pose.Pose()
        ideal_res_pose = pyrosetta.rosetta.core.pose.Pose()
        query_res_pose.append_residue_by_bond(rot_q)
        ideal_res_pose.append_residue_by_bond(rot_i)
        rmsd_list.append(
            rotamer_rmsd(query_res_pose, ideal_res_pose, super=super)
        )

    return rmsd_list


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


def kd_tree_from_chis(chis_list):
    chis_array = np.array(chis_list)
    return KDTree(chis_array), chis_array


def process_misses(restype, ideal_chis_list, query_chis_list):
    nearest_chi_tree, chi_data = kd_tree_from_chis(ideal_chis_list)

    ideal_res = pyrosetta.rosetta.core.conformation.Residue(restype, False)
    query_res = pyrosetta.rosetta.core.conformation.Residue(restype, False)
    ideal_res_pose = pyrosetta.rosetta.core.pose.Pose()
    ideal_res_pose.append_residue_by_bond(ideal_res)
    query_res_pose = pyrosetta.rosetta.core.pose.Pose()
    query_res_pose.append_residue_by_bond(query_res)

    rmsd_list = []
    ideal_chi_list = []

    for chis in query_chis_list:
        ideal_chis_index = nearest_chi_tree.query(chis)[1]
        ideal_chis = chi_data[ideal_chis_index]
        # set_rotamer_chis(ideal_res_pose.residue(1), *ideal_chis)
        # set_rotamer_chis(query_res_pose.residue(1), *chis)
        # rmsd_list.append(
        #     rotamer_rmsd(ideal_res_pose, query_res_pose, super=super)
        # )
        rot_i = set_rotamer_chis(ideal_res, *ideal_chis)
        rot_q = set_rotamer_chis(query_res, *chis)
        query_res_pose = pyrosetta.rosetta.core.pose.Pose()
        ideal_res_pose = pyrosetta.rosetta.core.pose.Pose()
        query_res_pose.append_residue_by_bond(rot_q)
        ideal_res_pose.append_residue_by_bond(rot_i)
        rmsd_list.append(
            rotamer_rmsd(query_res_pose, ideal_res_pose, super=super)
        )
        ideal_chi_list.append(ideal_chis)
    return rmsd_list, ideal_chi_list


def chi_list_to_query_df(chi_data_list):
    """
    Convert rotamer data list to a dataframe with chis and RTs (for all P bases)

    Takes test rotamer data in the form [(resname,chi_1,...,chi_n),...,(etc)]
    """
    rt_chi_key_map = {}
    rotamer_rt_map = {}
    for rot_info in chi_data_list:
        resname = rot_info[0]

        if resname not in rt_chi_key_map.keys():
            residue = get_rotamer_residue_from_name(resname)
            rotamer_rt_map[resname] = RotamerRTArray(
                residue=residue,
                target_atoms=("P", "P", "OH", "O2P"),
                inverse=True,
            )
            rt_chi_key_map[resname] = {"chis": [], "RTs": []}
        rotamer_rt_map[resname].set_all_chi(*rot_info[1:])
        rts = phospho_residue_inverse_rotamer_rts(
            rotamer_rt_map[resname].residue,
            rotamer_rt_array=rotamer_rt_map[resname],
        )
        rt_chi_key_map[resname]["RTs"].extend(rts)
        rt_chi_key_map[resname]["chis"].extend(
            [tuple(rot_info[1:])] * len(rts)
        )

    query_df = pd.DataFrame()
    for resname in rt_chi_key_map:
        resname_frame = pd.DataFrame(rt_chi_key_map[resname])
        resname_frame["resname"] = resname
        query_df = pd.concat([query_df, resname_frame])

    return query_df, rotamer_rt_map


def chi_list_to_query(chi_data_list, target_atoms=[]):
    """
    Convert rotamer data list to a dataframe with chis and RTs (for all P bases)

    Takes test rotamer data in the form [(resname,chi_1,...,chi_n),...,(etc)]
    """
    rt_chi_key_map = {}
    rotamer_rt_map = {}
    for rot_info in chi_data_list:
        resname = rot_info[0]

        if resname not in rt_chi_key_map.keys():
            residue = get_rotamer_residue_from_name(resname)
            target_atoms = (
                target_atoms
                if target_atoms
                else get_base_atoms_from_name(resname)[0]
            )
            rotamer_rt_map[resname] = RotamerRTArray(
                residue=residue,
                base_atoms=[2, 1, 2, 3],
                target_atoms=target_atoms,
                inverse=True,
            )
            rt_chi_key_map[resname] = {"chis": [], "RTs": []}
        rotamer_rt_map[resname].set_all_chi(*rot_info[1:])
        rts = residue_inverse_rotamer_rts(
            rotamer_rt_map[resname],
            [target_atoms],  # get_base_atoms_from_name(resname)
        )
        rt_chi_key_map[resname]["RTs"].extend(rts)
        rt_chi_key_map[resname]["chis"].extend(
            [tuple(rot_info[1:])] * len(rts)
        )

    query_rts = []
    query_chis = []
    query_resnames = []
    for resname in rt_chi_key_map:
        query_rts.extend(rt_chi_key_map[resname]["RTs"])
        query_chis.extend(rt_chi_key_map[resname]["chis"])
        query_resnames.extend([resname] * len(rt_chi_key_map[resname]["chis"]))
    return query_resnames, query_rts, query_chis, rotamer_rt_map


def test_hash_table(
    table,
    test_data,
    xbin_cart_resl=1,
    xbin_ori_resl=15,
    hdf5_path="",
    super_hits_before_rmsd=True,
    target_atoms=[],
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
    query_resnames, query_rts, query_chis, rotamer_rt_map = chi_list_to_query(
        test_data, target_atoms=target_atoms
    )
    logger.debug("got query and rotamer")
    binner = xb(xbin_cart_resl, xbin_ori_resl)
    rts_array = np.array(query_rts)
    test_keys = binner.get_bin_index(rts_array)
    logger.debug("generated keys to check table")
    hits_mask = table.contains(test_keys)

    miss_mask = hits_mask == False

    hits_keys = test_keys[hits_mask]
    logger.debug(f"hits found: {len(hits_keys)} out of {len(test_keys)}")

    values = table[hits_keys]

    chis_ideal_list = []

    with h5py.File(hdf5_path, "r") as f:
        for index in values:
            chis_ideal_list.append(f["chis"][index])
        ideal_chis_for_kd = f["ideal_chis"][:]

    if len(rotamer_rt_map.keys()) != 1:
        raise NotImplementedError("It happened again")

    # print(f"chis ideal list: {chis_ideal_list}")
    # print(f"ideal_chis_for_kd : {ideal_chis_for_kd}")
    restype = rotamer_rt_map["PHY"].residue.type()
    hits_chis = np.array(query_chis)[hits_mask]
    hits_rmsd_list = process_hits_by_rmsd(
        restype, chis_ideal_list, list(hits_chis), super=super_hits_before_rmsd
    )
    miss_chis = np.array(query_chis)[miss_mask]
    kd_tree_chis = np.unique(np.array(ideal_chis_for_kd), axis=0)
    misses_rmsd_list, matched_chi_list = process_misses(
        restype, kd_tree_chis, miss_chis
    )

    # make the return dfs here:
    hits_df = pd.DataFrame(
        {
            "key": hits_keys,
            "chis": list(hits_chis),
            "ideal_chis": chis_ideal_list,
            "rt": list(rts_array[hits_mask]),
            "resname": ["PHY"] * len(hits_chis),
            "index": values,
            "rmsd_to_hit": hits_rmsd_list,
        }
    )
    misses_df = pd.DataFrame(
        {
            "key": test_keys[miss_mask],
            "chis": list(miss_chis),
            "ideal_chis": matched_chi_list,
            "rt": list(rts_array[miss_mask]),
            "resname": ["PHY"] * len(miss_chis),
            "rmsd_to_ideal": misses_rmsd_list,
        }
    )
    return hits_df, misses_df


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
