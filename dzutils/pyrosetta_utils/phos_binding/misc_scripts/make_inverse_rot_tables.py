from os import makedirs
from os.path import isfile, isdir
from itertools import product, compress

import click

import pyrosetta
import pyrosetta.rosetta as pyr

import pandas as pd
import getpy as gp
import numpy as np

from xbin import XformBinner as xb

from dzutils.pyrosetta_utils import (
    residue_type_from_name3,
    get_rotamer_pose_from_residue_type,
)

from dzutils.pyrosetta_utils.phos_binding import (
    phospho_residue_inverse_rotamer_rts,
)

from dzutils.pyrosetta_utils.geometry.pose_xforms import (
    generate_pose_rt_between_res,
)


def mask_lists(mask, *lists):
    new_list_tuples = zip(*compress(zip(*lists), mask))
    return [list(tuple) for tuple in new_list_tuples]


def pose_from_res(res):
    pose = pyrosetta.rosetta.core.pose.Pose()
    pose.append_residue_by_bond(res)
    return pose


def save_table_as_json(out_dir, table, overwrite=False):
    """
    """
    # check if the out_dir exists

    data_out_path = f"{out_dir}/{table.name}.json"
    if not overwrite:
        if bool(not (isfile(data_out_path))):
            raise RuntimeError(
                f"Overwrite set to {overwrite} and file '{data_out_path}'' exists"
            )
    table.to_json(data_out_path)


def save_dict_as_bin(dict_out_dir, gp_dict, dict_name, overwrite=False):
    """
    """
    dict_out_path = f"{dict_out_dir}/{dict_name}.bin"
    if not overwrite:
        if bool(not (isfile(dict_out_path))):
            raise RuntimeError(
                f"Overwrite set to {overwrite} and file '{dict_out_path}'' exists"
            )
    gp_dict.dump(dict_out_path)


def rt_from_chis(*chis, res_type=None, alignment_atoms=()):
    """
    rt from the chis given
    """
    residue_pose = get_rotamer_pose_from_residue_type(res_type)
    return generate_pose_rt_between_res(residue_pose, 1, 1, alignment_atoms)


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


@click.command()
@click.option("-o", "--res-out-dir")
@click.option("-a", "--angle-res", default=15.0)
@click.option("-d", "--angstrom-dist-res", default=1.0)
@click.option("-r", "--run-name", default="inverse_ptr_exchi7_rotamers")
@click.option("-e", "--erase/--no-erase", default=False)
@click.option("-m", "--metrics-out-dir", default=False)
def main(
    run_name="inverse_ptr_exchi7_rotamers",
    res_out_dir="",
    angstrom_dist_res=1,
    angle_res=15,
    erase=False,
    metrics_out_dir=False,
):
    pyrosetta.init(
        """-out:level 100
        -packing:extrachi_cutoff 0
        -packing:ex1:level 4
        -packing:ex2:level 4
        -packing:ex3:level 4
        -packing:ex4:level 4

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
    print("debug hardcode! Only tests a tiny amount of rotamers")
    rts, chis_index, alignment_atoms = zip(
        *[
            (rt, chi_set, align)
            for chi_set, res in zip(chis, residues)
            for rt, align in phospho_residue_inverse_rotamer_rts(
                res, alignment_atoms=True
            )[:10]
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
    )
    inv_rot_table["index"] = inv_rot_table.index
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

    # Max granularity to go to
    cycle_depth = 1
    # degrees around the orignal value to search
    search_radius = 2

    # The actual jitter-search
    chis_list = list(chis_index)
    new_chis = []
    new_rts = []
    new_atoms = []
    run_info = []
    for i in range(2, cycle_depth + 2):
        for chis, atoms in zip(chis_list, alignment_atoms):
            for fine_chis in product(
                *[
                    np.linspace(chi - search_radius, chi + search_radius, i)
                    for chi in chis
                ]
            ):
                new_chis.append(tuple(chis))
                new_rts.append(
                    rt_from_chis(
                        *chis, res_type=res_type, alignment_atoms=atoms
                    )
                )

                new_atoms.append(atoms)

        new_keys = binner.get_bin_index(np.array(rts))
        mask = gp_dict.contains(np.array(new_keys))
        prev_mask = mask == False

        #
        # t1, t2, t3, t4 = zip(
        #     *compress(
        #         zip(new_chis, new_keys, new_atoms, new_rts), prev_mask
        #     )
        # )
        # new_chis = list(t1)
        # new_keys = list(t2)
        # new_atoms = list(t3)
        # new_rts = list(t4)
        new_chis, new_keys, new_atoms, new_rts = mask_lists(
            mask, new_chis, new_keys, new_atoms, new_rts
        )
        prev_chis, prev_keys, prev_atoms, prev_rts = mask_lists(
            prev_mask, new_chis, new_keys, new_atoms, new_rts
        )
        new_table = pd.DataFrame(
            {
                "key_int": new_keys,
                "chis": new_chis,
                "alignment_atoms": new_atoms,
                "rt": new_rts,
            }
        )

        prev_table = pd.DataFrame(
            {
                "key_int": prev_keys,
                "chis": prev_chis,
                "alignment_atoms": prev_atoms,
                "rt": prev_rts,
            }
        )
        if metrics_out_dir:
            prev_table.name = f"{run_name}_cycle_{cycle_depth}_search_{search_radius}_prev_hits.json"
            new_table.name = f"{run_name}_cycle_{cycle_depth}_search_{search_radius}_new_hits.json"
            save_table_as_json(metrics_out_dir, prev_table)
            save_table_as_json(metrics_out_dir, new_table)

        new_table = pd.concat(inv_rot_table, new_table)
        new_table = new_table.reset_index(drop=True)
        inv_rot_table = new_table
        inv_rot_table["index"] = inv_rot_table.index
        np_new_keys = np.array(inv_rot_table["key_int"], dtype=np.int64)
        np_new_vals = np.array(inv_rot_table["index"], dtype=np.int64)
        gp_dict[np_new_keys] = np_new_vals
        chis = new_chis
        atoms = new_atoms

    """
    data_out_path = f"{data_store_path}/{inv_rot_table.name}.json"
    if not erase:
        assert bool(
            not (isfile(data_out_path))
        ), "Table with this name already exists"
    inv_rot_table.to_json(data_out_path)
    data_name = f"{run_name}_{angstrom_dist_res}_ang_{angle_res}_deg"
    inv_rot_table.name = data_name
    """
    db_path = (
        "/home/dzorine/phos_binding/pilot_runs/loop_grafting/fragment_tables"
    )
    table_out_dir = f"{db_path}/inverse_rotamer/tables"
    data_name = f"{run_name}_{angstrom_dist_res}_ang_{angle_res}_deg"
    inv_rot_table.name = data_name
    save_table_as_json(table_out_dir, inv_rot_table)

    dict_out_dir = f"{db_path}/inverse_rotamer/dicts/"
    """
    if not erase:
        assert bool(
            not (isfile(dict_out_path))
        ), "Dict with this name already exists"
    gp_dict.dump(dict_out_path)
    """
    save_dict_as_bin(dict_out_dir, gp_dict)


if __name__ == "__main__":
    main()
