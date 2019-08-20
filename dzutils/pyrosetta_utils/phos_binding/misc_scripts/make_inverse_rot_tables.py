from os.path import isfile
import click

import pyrosetta
import pyrosetta.rosetta as pyr

import pandas as pd
import getpy as gp
import numpy as np

from xbin import XformBinner as xb

# import dzutils.pyrosetta_utils.geometry.superposition_utilities as su
from dzutils.pyrosetta_utils.phos_binding import (
    phospho_residue_inverse_rotamer_rts,
)


def dump_residue_as_pdb(residue, path):
    """
    Creates a pose with just the res given and dumps pdb at path

    raises an exception if dump_pdb is false
    """
    pose = pyr.core.pose.Pose()
    pose.append_residue_by_bond(residue)
    assert pose.dump_pdb(path), "dumping pdb failed!"
    return path


@click.command()
@click.option("-o", "--res-out-dir")
@click.option("-a", "--angle-res", default=15)
@click.option("-d", "--angstrom-dist-res", default=1)
@click.option("-r", "--run-name", default="inverse_ptr_exchi7_rotamers")
@click.option("-e", "--erase/--no-erase", default=False)
def main(
    run_name="inverse_ptr_exchi7_rotamers",
    res_out_dir="",
    angstrom_dist_res=1,
    angle_res=15,
    erase=False,
):
    pyrosetta.init(
        """-out:level 100
        -packing:extrachi_cutoff 0
        -packing:ex1:level 2
        -packing:ex2:level 2

    """
        # -extra_res_fa /home/dzorine/phos_binding/p_compounds/residues/PTR.params
    )
    # get chemical manager:
    chemical_manager = (
        pyrosetta.rosetta.core.chemical.ChemicalManager.get_instance()
    )
    rts = chemical_manager.residue_type_set(
        pyrosetta.rosetta.core.chemical.TypeSetMode.FULL_ATOM_t
    )
    # Get phospho-TYR
    res_type = rts.get_residue_type_with_variant_added(
        rts.name_map("TYR"),
        pyrosetta.rosetta.core.chemical.VariantType.PHOSPHORYLATION,
    )
    rots = [
        rot
        for rot in pyr.core.pack.rotamer_set.bb_independent_rotamers(res_type)
    ]

    binner = xb(cart_resl=angstrom_dist_res, ori_resl=angle_res)
    rt_dicts = []
    if res_out_dir:
        # out_dir = "/home/dzorine/phos_binding/pilot_runs/loop_grafting/ptr_rotamers/ptr_exchi7_rotamer_set"
        rt_dicts = [
            {
                "key_int": binner.get_bin_index(rt),
                "chis": [rot.chi(chi) for chi in range(1, rot.nchi() + 1)],
                "file": dump_residue_as_pdb(
                    rot,
                    f"{res_out_dir}/{res_type.name3().lower()}_rotamer_{i}.pdb",
                ),
            }
            for i, rot in enumerate(rots, 1)
            for rt in phospho_residue_inverse_rotamer_rts(rot)
        ]
    else:
        rt_dicts = [
            {
                "key_int": binner.get_bin_index(rt),
                "chis": [rot.chi(chi) for chi in range(1, rot.nchi() + 1)],
            }
            for i, rot in enumerate(rots, 1)
            for rt in phospho_residue_inverse_rotamer_rts(rot)
        ]

    inv_rot_table = pd.DataFrame(rt_dicts)
    inv_rot_table["index"] = inv_rot_table.index
    data_name = f"{run_name}_{angstrom_dist_res}_ang_{angle_res}_deg"
    inv_rot_table.name = data_name
    db_path = (
        "/home/dzorine/phos_binding/pilot_runs/loop_grafting/fragment_tables"
    )
    data_store_path = f"{db_path}/inverse_rotamer/tables"
    data_out_path = f"{data_store_path}/{inv_rot_table.name}.json"
    if not erase:
        assert bool(
            not (isfile(data_out_path))
        ), "Table with this name already exists"
        inv_rot_table.to_json(data_out_path)

    inv_rot_keys = np.array(inv_rot_table["key_int"], dtype=np.int64)
    inv_rot_vals = np.array(inv_rot_table["index"], dtype=np.int64)
    key_type = np.dtype("i8")
    value_type = np.dtype("i8")
    gp_dict = gp.Dict(key_type, value_type)
    gp_dict[inv_rot_keys] = inv_rot_vals
    dict_out_path = f"{db_path}/inverse_rotamer/dicts/{data_name}.bin"
    if not erase:
        assert bool(
            not (isfile(dict_out_path))
        ), "Dict with this name already exists"
    gp_dict.dump(dict_out_path)


if __name__ == "__main__":
    main()
