from os.path import isfile
import logging

import click
import getpy as gp
import h5py
import numpy as np
import pandas as pd
from xbin import XformBinner as xb


import pyrosetta

from dzutils.pyrosetta_utils import residue_type_from_name3
from dzutils.pyrosetta_utils.phos_binding import pres_bases
from dzutils.pyrosetta_utils.geometry.pose_xforms import RotamerRTArray


logger = logging.getLogger("make_inverse_rot_tables")
logger.setLevel(logging.DEBUG)


def rt_from_chis(rotamer_rt_array, *chis, base_atoms=(), target_atoms=()):
    rotamer_rt_array.reset_rotamer(
        *chis, base_atoms=base_atoms, target_atoms=target_atoms
    )
    # logger.debug(chis)
    return np.array(rotamer_rt_array)


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


def generate_rosetta_rotamer_chis(res_type):
    rots = [
        rot
        for rot in pyrosetta.rosetta.core.pack.rotamer_set.bb_independent_rotamers(
            res_type
        )
    ]
    return [
        tuple(rot.chi(chi) for chi in range(1, rot.nchi() + 1)) for rot in rots
    ]


def rosetta_rot_data(restype):

    chis = generate_rosetta_rotamer_chis(restype)

    logger.debug(f"working on rotamer count: {len(chis)}")

    residue = pyrosetta.rosetta.core.conformation.Residue(restype, True)
    possible_rt_bases = pres_bases(residue)

    rotamer_rt = RotamerRTArray(
        residue=residue, target_atoms=possible_rt_bases[0], inverse=True
    )

    rts, chis_index, alignment_atoms = zip(
        *[
            (
                rt_from_chis(rotamer_rt, *chi_set, target_atoms=atoms),
                chi_set,
                np.array(atoms),
            )
            for chi_set in chis
            for atoms in possible_rt_bases
        ]
    )

    return rts, chis_index, alignment_atoms


def write_hdf5_rotamer_hash_data(
    path,
    restype,
    rts,
    chis_index,
    alignment_atoms,
    cart_resl=1,
    ori_resl=15,
    key_label="key_int",
    chi_label="chis",
    align_atom_label="alignment_atoms",
    rt_label="rt",
    ideal=False,
):
    """
    """
    binner = xb(cart_resl=cart_resl, ori_resl=ori_resl)
    keys = binner.get_bin_index(np.array(rts))
    with h5py.File(path, "w") as f:
        for label, data in zip(
            ("index", key_label, rt_label, chi_label, align_atom_label),
            (
                [*range(1, len(keys) + 1)],
                keys,
                rts,
                chis_index,
                alignment_atoms,
            ),
        ):
            np_data = np.array(data)
            f.create_dataset(label, np_data.shape, data=np_data)
        if ideal:
            np_data = np.array(chis_index)
            f.create_dataset("ideal_chis", np_data.shape, data=np_data)


def write_hdf5_rosetta_rotamer_hash_data(
    path,
    restype,
    cart_resl=1,
    ori_resl=15,
    key_label="key_int",
    chi_label="chis",
    align_atom_label="alignment_atoms",
    rt_label="rt",
    ideal=False,
):
    """
    This is the utility that generates backbone dep rotamers and hashes them

    Bin dimensions and hdf5 labels are configurable, but rotamer generation is
    just whatever ex level rosetta has been set to before calling this
    """
    rts, chis_index, alignment_atoms = rosetta_rot_data(restype)
    write_hdf5_rotamer_hash_data(
        path,
        restype,
        rts,
        chis_index,
        alignment_atoms,
        cart_resl=cart_resl,
        ori_resl=ori_resl,
        key_label=key_label,
        chi_label=chi_label,
        align_atom_label=align_atom_label,
        rt_label=rt_label,
        ideal=ideal,
    )


def build_inv_rot_table_from_rosetta_rots(
    restype,
    cart_resl=1,
    ori_resl=15,
    key_label="key_int",
    chi_label="chis",
    align_atom_label="alignment_atoms",
    rt_label="rt",
):
    """
    """
    rts, chis_index, alignment_atoms = rosetta_rot_data(restype)
    binner = xb(cart_resl=cart_resl, ori_resl=ori_resl)
    keys = binner.get_bin_index(np.array(rts))

    # make the table to assign rotamers to
    inv_rot_table = pd.DataFrame(
        {
            key_label: keys,
            chi_label: chis_index,
            align_atom_label: alignment_atoms,
            rt_label: rts,
        }
    )

    return inv_rot_table


@click.command()
@click.option("-a", "--angle-res", default=15.0)
@click.option("-d", "--angstrom-dist-res", default=1.0)
@click.option("-r", "--run-name", default="inverse_ptr_exchi7_rotamers")
@click.option("-e", "--erase/--no-erase", default=False)
@click.option("-i", "--ideal/--dont-include", default=False)
def main(
    run_name="inverse_ptr_exchi7_rotamers",
    angstrom_dist_res=1,
    angle_res=15,
    erase=False,
    ideal=False,
):

    db_path = "/home/dzorine/phos_binding/fragment_tables"
    table_out_dir = f"{db_path}/inverse_rotamer/hdf5/test_tables"
    data_name = f"{run_name}_{angstrom_dist_res}_ang_{angle_res}_deg"

    # save_table_as_json(table_out_dir, inv_rot_table, overwrite=erase)
    dict_out_dir = f"{db_path}/inverse_rotamer/dicts/"
    """
            -packing:extrachi_cutoff 0
            -packing:ex1:level 1
            -packing:ex2:level 1
            -packing:ex3:level 1
            -packing:ex4:level 1
    """
    pyrosetta.init(
        """-out:level 100
        -extra_res_fa /home/dzorine/projects/phos_binding/params_files/p_loop_params/PHY_4_chi.params
        """
    )

    restype = residue_type_from_name3("PHY")

    rts, chis_index, alignment_atoms = rosetta_rot_data(restype)
    binner = xb(cart_resl=angstrom_dist_res, ori_resl=angle_res)
    keys = binner.get_bin_index(np.array(rts))

    # Make the base dictionary
    key_type = np.dtype("i8")
    value_type = np.dtype("i8")
    gp_dict = gp.Dict(key_type, value_type)

    indices = np.array([*range(len(chis_index))], dtype=np.dtype("i8"))
    gp_dict[keys] = indices
    write_hdf5_rotamer_hash_data(
        f"{table_out_dir}/{data_name}_base_table.hf5",
        restype,
        rts,
        chis_index,
        alignment_atoms,
        angstrom_dist_res,
        angle_res,
        ideal=ideal,
    )
    save_dict_as_bin(dict_out_dir, gp_dict, data_name, overwrite=erase)


if __name__ == "__main__":
    main()
