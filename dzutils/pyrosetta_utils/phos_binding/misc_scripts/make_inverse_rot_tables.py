import logging

import click
import getpy as gp
import numpy as np
from xbin import XformBinner as xb


import pyrosetta

from dzutils.pyrosetta_utils import residue_type_from_name3

from dzutils.pyrosetta_utils.phos_binding.misc_scripts.rotable import (
    save_dict_as_bin,
    rosetta_rot_data,
    write_hdf5_data,
)

logger = logging.getLogger("make_inverse_rot_tables")
logger.setLevel(logging.DEBUG)


@click.command()
@click.argument("resname", nargs=1)
@click.argument("rt_bases", nargs=1, help="")
@click.option("-a", "--angle-res", default=15.0)
@click.option("-d", "--angstrom-dist-res", default=1.0)
@click.option("-r", "--run-name", default="inverse_ptr_exchi7_rotamers")
@click.option("-e", "--erase/--no-erase", default=False)
@click.option("-i", "--ideal/--dont-include", default=False)
def main(
    resname,
    rt_bases,
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

    pyrosetta.init(
        """-out:level 100
        -extra_res_fa /home/dzorine/projects/phos_binding/params_files/p_loop_params/PHY_4_chi.params
        """
    )
    # resname = "PHY"
    restype = residue_type_from_name3(resname)
    residue = pyrosetta.rosetta.core.conformation.Residue(restype, False)

    possible_rt_bases = [
        [residue.atom_index(name) for name in ["P", "P", "OH", "O2P"]]
    ]

    rts, chis_index, alignment_atoms = rosetta_rot_data(
        restype, possible_rt_bases
    )
    binner = xb(cart_resl=angstrom_dist_res, ori_resl=angle_res)
    rts_array = np.array(rts)
    keys = binner.get_bin_index(rts_array)

    # Make the base dictionary
    key_type = np.dtype("i8")
    value_type = np.dtype("i8")
    gp_dict = gp.Dict(key_type, value_type)

    indices = np.array([*range(len(chis_index))], dtype=np.dtype("i8"))
    gp_dict[keys] = indices
    chis_array = np.array(chis_index)
    data_dict = {
        "key_int": ("", keys, {}),
        "chis": (
            "",
            chis_array,
            {
                "target_atoms": np.array(possible_rt_bases[0]),
                "base_atoms": np.array([2, 1, 2, 3]),
                "num_chis": chis_array.shape[1],
                "residue_name": resname,
            },
        ),
        "rt": ("", rts_array, {}),
        "ideal_chis": ("", chis_array, {}),
    }
    write_hdf5_data(f"{table_out_dir}/{data_name}_base_table.hf5", **data_dict)

    save_dict_as_bin(dict_out_dir, gp_dict, data_name, overwrite=erase)


if __name__ == "__main__":
    main()
