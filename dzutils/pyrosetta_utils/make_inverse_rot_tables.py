from os.path import isfile

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


pyrosetta.init(
    """-out:level 100
    -extra_res_fa /home/dzorine/phos_binding/p_compounds/residues/PTR.params
    -packing:ex1
    -packing:ex2
    -packing:ex3
    -packing:ex4
    -packing:ex1:level 7
    -packing:ex2:level 7
    -packing:ex3:level 7
    -packing:ex4:level 7

"""
    # -extra_res_fa /home/dzorine/phos_binding/p_compounds/residues/PTR.params
)
pose = pyrosetta.pose_from_sequence("Y")
mutate_res = pyrosetta.rosetta.protocols.simple_moves.MutateResidue()
phy_index = 1
mut_index = phy_index
mut_name = "PTR"
mutate_res.set_target(str(mut_index))
mutate_res.set_res_name(mut_name)
mutate_res.apply(pose)
res_type = pose.residue(1).type()
rots = [
    rot for rot in pyr.core.pack.rotamer_set.bb_independent_rotamers(res_type)
]

out_dir = "/home/dzorine/phos_binding/pilot_runs/loop_grafting/ptr_rotamers/ptr_exchi7_rotamer_set"
binner = xb()
rt_dicts = [
    {
        "file": dump_residue_as_pdb(
            rot, f"{out_dir}/{mut_name.lower()}_rotamer_{i}.pdb"
        ),
        "key_int": binner.get_bin_index(rt),
    }
    for i, rot in enumerate(rots, 1)
    for rt in phospho_residue_inverse_rotamer_rts(rot)
]

inv_rot_table = pd.DataFrame(rt_dicts)
inv_rot_table["index"] = inv_rot_table.index
data_name = "inverse_ptr_exchi7_rotamers_1ang_v1"
inv_rot_table.name = data_name
db_path = "/home/dzorine/phos_binding/pilot_runs/loop_grafting/fragment_tables"
data_store_path = f"{db_path}/inverse_rotamer/tables"
data_out_path = f"{data_store_path}/{inv_rot_table.name}.json"
assert bool(not (isfile(data_out_path))), "Table with this name already exists"
inv_rot_table.to_json(data_out_path)

inv_rot_keys = np.array(inv_rot_table["key_int"], dtype=np.int64)
inv_rot_vals = np.array(inv_rot_table["index"], dtype=np.int64)
key_type = np.dtype("i8")
value_type = np.dtype("i8")
gp_dict = gp.Dict(key_type, value_type)
gp_dict[inv_rot_keys] = inv_rot_vals
dict_out_path = f"{db_path}/inverse_rotamer/dicts/{data_name}.bin"
assert bool(not (isfile(dict_out_path))), "Dict with this name already exists"
gp_dict.dump(dict_out_path)
