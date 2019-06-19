from os.path import isfile

import pandas as pd
import getpy as gp
import numpy as np

import pyrosetta

from dzutils.pyrosetta_utils.phos_binding import get_loop_xform_dicts
from dzutils.pyrosetta_utils import residues_with_element
from dzutils.sutils import read_flag_file
from dzutils.pdb_file_utils import pdb_files_in_dir
from dzutils.pyrosetta_utils.chain_utils import rechain_resname


# from dzutils.pyrosetta_utils.geometry.pose_xforms import generate_pose_rt_between_res,get_func_to_end


def maybe_load(pdb):
    """
    Dumb hack so rosettta can tolerate pdbs it can't load
    """
    try:
        pose = pyrosetta.pose_from_file(pdb)
        return pose
    except Exception as e:
        print(e)
        return pyrosetta.pose_from_sequence("AAA")


ploop_flags_file = "/home/dzorine/phos_binding/pilot_runs/loop_grafting/initial_testing/misc_files/p_ligand.flags"
flags = read_flag_file(ploop_flags_file)
flags_str = " ".join(flags.replace("\n", " ").split())
pyrosetta.init(flags_str)


pdb_dir = "/home/dzorine/phos_binding/ploop_set_1"

pdb_store = "/home/dzorine/phos_binding/ploop_pdbs/ploop_set_1_no_p"
pdb_paths = pdb_files_in_dir(pdb_dir)[:100]
# subset = pdb_paths
poses = (maybe_load(p) for p in pdb_paths)
dicts = list()
for pose in poses:
    names = list(
        set(
            [pose.residue(i).name3() for i in residues_with_element(pose, "P")]
        )
    )
    newp = pose.clone()
    if len(newp.residues) < 4:
        continue
    for name in names:
        newp = rechain_resname(newp, name)

    # Rename and deposit
    base_name_old = newp.pdb_info().name().split("/")[-1].split(".pdb")[0]
    new_name = f"{base_name_old}_no_p.pdb"
    new_out_path = f"{pdb_store}/{new_name}"
    newp.dump_pdb(new_out_path)
    # Check if the first and last res of chain 1 are free to spin
    # if not, pre/append a residue, fix omegas
    xform_dicts = get_loop_xform_dicts(newp, 3)
    if xform_dicts:
        dicts.extend(xform_dicts)

loop_table = pd.DataFrame(dicts)
# At this point table is made with things ready to spin

# begin dict filling loop
# Max cycles 10k, raise RuntimeError if exceeded
for i in range(10000):
    """
    """
# Map a spin ends func (randomize rama on the ends)
# Make a new table
# Check % already in dict
# mask table by false
# deposit masked table, append to master table
# repeat until 95% of keys are repeats

loop_table["index"] = loop_table.index

data_name = "test_spun_1ang_100_ploops_v1"
data_store_path = "/home/dzorine/phos_binding/pilot_runs/loop_grafting/fragment_tables/ploops_expanded_set_1"
print(len(loop_table.index))
table_out_path = f"{data_store_path}/tables/{data_name}.json"
assert bool(
    not (isfile(table_out_path))
), "Table with this name already exists"
loop_table.to_json(table_out_path)
# loop_table = pd.read_json(f"{data_store_path}/tables/{data_name}.json")
e2e_keys = np.array(loop_table["key_int"], dtype=np.int64)
e2e_vals = np.array(loop_table["index"], dtype=np.int64)
key_type = np.dtype("i8")
value_type = np.dtype("i8")
gp_dict = gp.Dict(key_type, value_type)
gp_dict[e2e_keys] = e2e_vals
dict_out_path = f"{data_store_path}/dicts/{data_name}.bin"
assert bool(not (isfile(dict_out_path))), "Dict with this name already exists"
gp_dict.dump(dict_out_path)
