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

ploop_flags_file = "/home/dzorine/phos_binding/pilot_runs/loop_grafting/initial_testing/misc_files/p_ligand.flags"
flags = read_flag_file(ploop_flags_file)
flags_str = " ".join(flags.replace("\n", " ").split())
pyrosetta.init(flags_str)


pdb_dir = "/home/dzorine/phos_binding/ploop_set_1/"
pdb_paths = pdb_files_in_dir(pdb_dir)
subset = pdb_paths
poses = (pyrosetta.pose_from_file(p) for p in subset)
dicts = list()
for sample in poses:
    names = list(
        set(
            [
                sample.residue(i).name3()
                for i in residues_with_element(sample, "P")
            ]
        )
    )
    newp = sample.clone()
    if len(newp.residues) < 4:
        continue
    for name in names:
        newp = rechain_resname(newp, name)
    xform_dicts = get_loop_xform_dicts(newp, 3)
    if xform_dicts:
        dicts.extend(xform_dicts)

pdf = pd.DataFrame(dicts)
pdf["index"] = pdf.index
data_name = "exp_no_spin_1ang_15k_ploops_v2"
data_store_path = "/home/dzorine/phos_binding/pilot_runs/loop_grafting/fragment_tables/ploops_expanded_set_1"
pdf.name = data_name
table_out_path = f"{data_store_path}/tables/{pdf.name()}.json"
assert bool(
    not (isfile(table_out_path))
), "Table with this name already exists"
pdf.to_json(table_out_path)
loop_table = pd.read_json(f"{data_store_path}/tables/{data_name}.json")
e2e_keys = np.array(loop_table["key_int"], dtype=np.int64)
e2e_vals = np.array(loop_table["index"], dtype=np.int64)
key_type = np.dtype("i8")
value_type = np.dtype("i8")
gp_dict = gp.Dict(key_type, value_type)
gp_dict[e2e_keys] = e2e_vals
dict_out_path = f"{data_store_path}/dicts/{data_name}.bin"
assert bool(not (isfile(dict_out_path))), "Dict with this name already exists"
gp_dict.dump(dict_out_path)
