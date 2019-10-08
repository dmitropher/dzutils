from os import makedirs
from os.path import isfile, isdir
import logging
from itertools import product, compress
import logging

import click

import pyrosetta
import pyrosetta.rosetta as pyr

import pandas as pd
import getpy as gp
import numpy as np
from numba import jit

from xbin import XformBinner as xb

from dzutils.pyrosetta_utils import (
    residue_type_from_name3,
    get_rotamer_pose_from_residue_type,
)

from dzutils.pyrosetta_utils.phos_binding import (
    phospho_residue_inverse_rotamer_rts,
    pres_bases,
)

from dzutils.pyrosetta_utils.geometry.pose_xforms import (
    generate_pose_rt_between_res,
    RotamerRTArray,
)


logger = logging.getLogger("make_inverse_rot_tables")
logger.setLevel(logging.DEBUG)


def it_cartesian(input_arrays):
    return product(*input_arrays)


# @jit(nopython=True)
def numba_cartesian(input_arrays, out=None):
    """
    Generate a cartesian product of input arrays.

    Parameters
    ----------
    arrays : list of array-like
        1-D arrays to form the cartesian product of.
    out : ndarray
        Array to place the cartesian product in.

    Returns
    -------
    out : ndarray
        2-D array of shape (M, len(arrays)) containing cartesian products
        formed of input arrays.

    Examples
    --------
    >>> numba_cartesian(([1, 2, 3], [4, 5], [6, 7]))
    array([[1, 4, 6],
           [1, 4, 7],
           [1, 5, 6],
           [1, 5, 7],
           [2, 4, 6],
           [2, 4, 7],
           [2, 5, 6],
           [2, 5, 7],
           [3, 4, 6],
           [3, 4, 7],
           [3, 5, 6],
           [3, 5, 7]])

    """

    arrays = np.empty(
        (len(input_arrays), len(input_arrays[0])), dtype=np.float64
    )
    for i in range(len(input_arrays)):
        # print (arrays[i])
        # print (input_arrays[i])
        arrays[i] = input_arrays[i]
    dtype = arrays[0].dtype

    n = np.prod(
        np.asarray(
            [arrays[i].size for i in range(len(input_arrays))], dtype=np.uint64
        )
    )

    if out is None:
        out = np.zeros((n, len(arrays)), dtype=dtype)

    m = int(n // arrays[0].size)

    out[:, 0] = np.repeat(arrays[0], m)
    array_len = len(input_arrays[1:])
    logger.debug(m)
    if array_len > 0:
        numba_cartesian(arrays[1:], out=out[0:m, 1:])
        for j in range(1, arrays[0].size):
            out[j * m : (j + 1) * m, 1:] = out[0:m, 1:]
    return out


def rt_from_chis(rotamer_rt_array, *chis, base_atoms=(), target_atoms=()):
    rotamer_rt_array.reset_rotamer(
        *chis, base_atoms=base_atoms, target_atoms=target_atoms
    )
    # logger.debug(chis)
    return np.array(rotamer_rt_array)


def merge_data_to_table(old_table, name, **new_data):
    """
    concats two tables and fixes the index/name
    """
    new_table = pd.DataFrame(new_data)
    new_table_merged = pd.concat((old_table, new_table), sort=False)
    new_table_merged = new_table_merged.reset_index(drop=True)
    new_table_merged["index"] = new_table_merged.index
    new_table_merged.name = old_table.name
    return new_table_merged


def update_hashmap_from_df(hashmap, df, key_label, value_label):
    """
    Puts the values at value_label into the keys at key_label from the df
    """
    np_new_keys = np.array(df[key_label], dtype=np.int64)
    np_new_vals = np.array(df[value_label], dtype=np.int64)
    hashmap[np_new_keys] = np_new_vals


def mask_lists(mask, *lists):
    new_list_tuples = zip(*compress(zip(*lists), mask))
    return [list(tuple) for tuple in new_list_tuples]


def pose_from_res(res):
    pose = pyrosetta.rosetta.core.pose.Pose()
    pose.append_residue_by_bond(res)
    return pose


def save_table_as_json(out_dir, table, overwrite=False):
    """
    wraps DataFrame.to_json
    """

    data_out_path = f"{out_dir}/{table.name}.json"
    if not overwrite:
        if isfile(data_out_path):
            raise RuntimeError(
                f"Overwrite set to {overwrite} and file '{data_out_path}' exists"
            )
    table.to_json(data_out_path)


def save_dict_as_bin(dict_out_dir, gp_dict, dict_name, overwrite=False):
    """
    wraps gp_dict.dump
    """
    dict_out_path = f"{dict_out_dir}/{dict_name}.bin"
    if not overwrite:
        if isfile(dict_out_path):
            raise RuntimeError(
                f"Overwrite set to {overwrite} and file '{dict_out_path}' exists"
            )
    gp_dict.dump(dict_out_path)


# def rt_from_chis(*chis, res_type=None, alignment_atoms=()):
#     """
#     rt from the chis given
#     """
#     residue_pose = get_rotamer_pose_from_residue_type(res_type)
#     return generate_pose_rt_between_res(residue_pose, 1, 1, alignment_atoms)


def dump_residue_as_pdb(residue, path):
    """
    Creates a pose with just the res given and dumps pdb at path

    raises an exception if dump_pdb is false
    """
    pose = pyr.core.pose.Pose()
    pose.append_residue_by_bond(residue)
    assert pose.dump_pdb(path), "dumping pdb failed!"
    return path


def generate_rosetta_rotamer_chis(res_type):
    rots = [
        rot
        for rot in pyr.core.pack.rotamer_set.bb_independent_rotamers(res_type)
    ]
    return zip(
        *[
            (tuple(rot.chi(chi) for chi in range(1, rot.nchi() + 1)), rot)
            for rot in rots
        ]
    )


def update_df_and_gp_dict(
    df, hashmap, cycle, new_chis, new_keys, new_atoms, new_rts
):
    """
    Updates df and hashmap, returns ref to both
    """
    mask = hashmap.contains(np.array(new_keys)) == False
    logger.debug("mask generated for batch")
    logger.debug(str(mask))
    logger.debug(f"all: {all(mask)} any: {any(mask)}")
    if any(mask):
        new_chis_masked, new_keys_masked, new_atoms_masked, new_rts_masked = mask_lists(
            mask, new_chis, new_keys, new_atoms, new_rts
        )
        new_data = {
            "key_int": new_keys_masked,
            "chis": new_chis_masked,
            "alignment_atoms": new_atoms_masked,
            "rt": new_rts_masked,
            "cycle": [cycle] * len(new_chis_masked),
        }

        new_df = merge_data_to_table(df, df.name, **new_data)
        update_hashmap_from_df(hashmap, new_df, "key_int", "index")
        return new_df, hashmap
    return df, hashmap


# @jit(nopython=True,cache=True)
def bit_pack_rotamers(chis_array, num_chis=3, round_fraction=0.01):
    """
    Converts the given chis to a numpy uint64

    Uses the formula: left shift the value by 64//(num_chis)*(num of chi - 1)

    bitwise or over the list
    """
    scaling = np.float64(360 / round_fraction)
    if (scaling) > np.power(np.uint64(2), np.uint64(64 // num_chis)):
        raise ValueError(
            "given round_fraction and number of chis cannot fit into uint 64"
        )
    # shape = chis_array.shape[1]
    # uint_chis_array = (np.mod(np.divide(chis_array ,round_fraction) ,scaling)).astype(np.uint64)

    for i in range(num_chis):
        chis_array[:, i] = np.mod(
            np.divide(chis_array[:, i], round_fraction), scaling
        )
    uint_chis_array = chis_array.astype(np.uint64)
    for i in range(num_chis):
        uint_chis_array[:, i] = np.left_shift(
            uint_chis_array[:, i], np.uint64(i * (64 // num_chis))
        )
    packed = uint_chis_array[:, 0]
    for i in range(1, num_chis):
        packed = np.bitwise_or(packed, uint_chis_array[:, i])
    return packed


def unpack_chis(packed_chis, num_chis=3, round_fraction=0.01):
    """
    Takes a np.array of bitpacked chis and unpacks them
    """

    unpacked_transpose = np.array(
        [
            np.right_shift(packed_chis, i * (64 // num_chis))
            for i in range(num_chis)
        ]
    )
    shifted = np.transpose(unpacked_transpose)
    mask = np.power(
        np.uint64(2), np.uint64(64 // num_chis), dtype=np.uint64
    ) - np.uint64(1)
    chis = np.array(
        np.bitwise_and(shifted, mask) * round_fraction, dtype=np.float64
    )
    return chis


# @jit(nopython=True,cache=True)
def expand_rotamer_set(
    chi_array, search_radius, granularity_factor, round_fraction
):
    rotamer_set = set([np.uint64(0)][1:])
    num_chi_sets = chi_array.shape[0]
    # logger.debug (num_chi_sets)

    for i in range(num_chi_sets):
        chis = chi_array[i]
        num_chis = len(chis)
        if (
            bit_pack_rotamers(
                np.array([chis]),
                num_chis=num_chis,
                round_fraction=round_fraction,
            )[0]
            in rotamer_set
        ):
            continue

        # chi_array = np.asarray(chis)
        col_space = np.empty((num_chis, granularity_factor), dtype=np.float64)
        for j in range(num_chis):
            col_space[j, :] = np.linspace(
                chis[j] - np.float64(search_radius),
                chis[j] + np.float64(search_radius),
                granularity_factor,
            )

        expanded = np.asarray(col_space, dtype=np.float64)
        logger.debug("expanded:")
        logger.debug(f"{expanded}")
        # fine_chis = np.array (list(it_cartesian(
        #     expanded
        # )), np.float64)
        fine_chis = np.array(numba_cartesian(expanded))
        logger.debug("fine_chis: ")
        logger.debug(f"{fine_chis}")
        uint_rot_repr = bit_pack_rotamers(
            fine_chis, num_chis=num_chis, round_fraction=round_fraction
        )
        if not len(rotamer_set):
            rotamer_set = set(uint_rot_repr)
            continue
        rotamer_set = rotamer_set.union(set(uint_rot_repr))
    return rotamer_set


@click.command()
@click.option("-o", "--res-out-dir")
@click.option("-a", "--angle-res", default=15.0)
@click.option("-d", "--angstrom-dist-res", default=1.0)
@click.option("-r", "--run-name", default="inverse_ptr_exchi7_rotamers")
@click.option("-e", "--erase/--no-erase", default=False)
@click.option("-p", "--ramp/--no-ramp", default=False)
@click.option("-m", "--metrics-out-dir", default=False)
@click.option("-b", "--batch-factor", default=10000)
@click.option("-g", "--granularity-factor", default=15)
@click.option("-c", "--cycle_depth", default=3)
@click.option("-s", "--search-radius", default=5)

# @click.option("-d","--rots-for)
def main(
    run_name="inverse_ptr_exchi7_rotamers",
    res_out_dir="",
    angstrom_dist_res=1,
    angle_res=15,
    erase=False,
    metrics_out_dir=False,
    batch_factor=10000,
    granularity_factor=15,
    cycle_depth=3,
    search_radius=5,
    ramp=False,
):

    db_path = (
        "/home/dzorine/phos_binding/pilot_runs/loop_grafting/fragment_tables"
    )
    table_out_dir = f"{db_path}/inverse_rotamer/tables"
    data_name = f"{run_name}_{angstrom_dist_res}_ang_{angle_res}_deg"

    # save_table_as_json(table_out_dir, inv_rot_table, overwrite=erase)
    dict_out_dir = f"{db_path}/inverse_rotamer/dicts/"

    pyrosetta.init(
        """-out:level 100
        -packing:extrachi_cutoff 0
        -packing:ex1:level 1
        -packing:ex2:level 1
        -packing:ex3:level 1
        -packing:ex4:level 1
    """
    )

    res_type = residue_type_from_name3(
        "TYR",
        variant=pyrosetta.rosetta.core.chemical.VariantType.PHOSPHORYLATION,
    )

    binner = xb(cart_resl=angstrom_dist_res, ori_resl=angle_res)

    chis, residues = generate_rosetta_rotamer_chis(res_type)
    if res_out_dir:
        if not isdir(res_out_dir):
            makedirs(res_out_dir)

        for i, res in enumerate(residues, 1):
            pose_from_res(res).dump_pdb(
                f"{res_out_dir}/{res_type.name3().lower()}_rotamer_{i}.pdb"
            )
    logger.debug(f"working on rotamer count: {len(residues)}")
    # logger.debug"danger hardcode to reduce rots")
    residue = pyrosetta.rosetta.core.conformation.Residue(res_type, True)
    possible_rt_bases = pres_bases(residue)
    # pose = _pyrosetta.rosetta.core.pose.Pose()
    # pose.append_residue_by_bond(residue)
    rotamer_rt = RotamerRTArray(
        residue=residue, target_atoms=("P", "P", "OH", "O2P"), inverse=True
    )

    rts, chis_index, alignment_atoms = zip(
        *[
            (
                rt_from_chis(rotamer_rt, *chi_set, target_atoms=atoms),
                chi_set,
                atoms,
            )
            for chi_set, res in zip(chis, residues)
            for atoms in possible_rt_bases
        ]
    )
    keys = binner.get_bin_index(np.array(rts))

    # make the table to assign rotamers to
    inv_rot_table = pd.DataFrame(
        {
            "key_int": keys,
            "chis": chis_index,
            "alignment_atoms": alignment_atoms,
            "rt": rts,
        }
    )  # .iloc[:500]
    inv_rot_table["index"] = inv_rot_table.index
    inv_rot_table["cycle"] = pd.Series([0] * len(inv_rot_table.index))
    data_name = f"{run_name}_{angstrom_dist_res}_ang_{angle_res}_deg"
    inv_rot_table.name = data_name

    # Make the base dictionary
    key_type = np.dtype("i8")
    value_type = np.dtype("i8")
    gp_dict = gp.Dict(key_type, value_type)

    # build the key value pairs and fill dict
    inv_rot_keys = np.array(inv_rot_table["key_int"], dtype=np.int64)
    inv_rot_vals = np.array(inv_rot_table["index"], dtype=np.int64)
    gp_dict[inv_rot_keys] = inv_rot_vals

    for i in range(1, cycle_depth + 1):
        # batch_factor = 1000
        count = 0
        batch_new_chi_ints = set()
        batch_new_chis = []
        chi_list = inv_rot_table[inv_rot_table["cycle"] == i - 1][
            "chis"
        ].to_list()
        chi_set = set(chi_list)
        chi_array = np.asarray(list(chi_set), np.float64)
        logger.debug(chi_array)
        logger.debug(i * granularity_factor)
        batch_new_chi_ints = expand_rotamer_set(
            chi_array,
            search_radius,
            (i * granularity_factor if ramp else granularity_factor),
            np.float64(0.01),
        )
        packed_chis = np.fromiter(batch_new_chi_ints, np.uint64)
        chi_sets = unpack_chis(
            packed_chis, num_chis=3, round_fraction=np.float(0.01)
        )
        # logger.debugf"batch: {count} exceeds factor: {batch_factor}")
        logger.debug(len(chi_sets))
        logger.debug(chi_sets)
        logger.debug(possible_rt_bases)
        if not len(chi_sets):
            break
        batch_new_atoms, batch_new_chis, batch_new_rts = [[], [], []]
        # [
        #     (atoms, tuple(chis), rt_from_chis(rotamer_rt, *chis, target_atoms=atoms))
        #     for chis in chi_sets
        #     for atoms in possible_rt_bases
        # ]
        for chis in chi_sets:
            for atoms in possible_rt_bases:
                batch_new_atoms.append(atoms)
                batch_new_chis.append(tuple(chis))
                batch_new_rts.append(
                    rt_from_chis(rotamer_rt, *chis, target_atoms=atoms)
                )
        logger.debug("new rts to be xbinned")
        logger.debug(f"{np.array(batch_new_rts)}")
        batch_new_keys = binner.get_bin_index(np.array(batch_new_rts))

        inv_rot_table, gp_dict = update_df_and_gp_dict(
            inv_rot_table,
            gp_dict,
            i,
            batch_new_chis,
            batch_new_keys,
            batch_new_atoms,
            batch_new_rts,
        )
        data_name = (
            f"{run_name}_{angstrom_dist_res}_ang_{angle_res}_deg_cycle_{i}"
        )
        inv_rot_table.name = data_name
        save_table_as_json(table_out_dir, inv_rot_table, overwrite=erase)
        save_dict_as_bin(dict_out_dir, gp_dict, data_name, overwrite=erase)

    db_path = (
        "/home/dzorine/phos_binding/pilot_runs/loop_grafting/fragment_tables"
    )
    table_out_dir = f"{db_path}/inverse_rotamer/tables"
    data_name = f"{run_name}_{angstrom_dist_res}_ang_{angle_res}_deg"
    inv_rot_table.name = data_name
    save_table_as_json(table_out_dir, inv_rot_table, overwrite=erase)

    dict_out_dir = f"{db_path}/inverse_rotamer/dicts/"

    save_dict_as_bin(dict_out_dir, gp_dict, data_name, overwrite=erase)


if __name__ == "__main__":
    main()
