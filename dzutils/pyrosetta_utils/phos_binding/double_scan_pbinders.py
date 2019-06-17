# Scan for ploop insertion

# Hunt for inverse rotamers to complete interaction

# Takes: pdb, out_dir, out_prefix (uses defaults)

import glob

import pandas as pd
import getpy as gp
import numpy as np

import pyrosetta

from xbin import XformBinner as xb

from dzutils.pyrosetta_utils.phos_binding import get_loop_xform_dicts
from dzutils.pyrosetta_utils import (
    residues_with_element,
    bonded_atoms,
    hbond_to_residue,
)
from dzutils.sutils import read_flag_file
from dzutils.pdb_file_utils import pdb_files_in_dir
from dzutils.pyrosetta_utils.chain_utils import rechain_resname
from dzutils.pyrosetta_utils.phos_binding import (
    loop_res_pairs,
    pair_to_rt_dict,
    loops_to_rt_dict,
    p_atoms_in_pose,
)
from dzutils.pyrosetta_utils.geometry.pose_xforms import (
    generate_pose_rt_between_res,
    get_func_to_end,
)


def load_table_and_dict(table_path, dict_path, key_type, value_type):
    """
    return a tuple of table,dict
    """
    gp_dict = gp.Dict(key_type, value_type)
    gp_dict.load(dict_path)
    table = pd.read_json(table_path)
    return table, gp_dict


ploop_flags_file = "/home/dzorine/phos_binding/pilot_runs/loop_grafting/initial_testing/misc_files/p_ligand.flags"
flags = read_flag_file(ploop_flags_file)
flags_str = " ".join(flags.replace("\n", " ").split())
pyrosetta.init(flags_str)

pdb_path = sys.argv[1]
pose = pyrosetta.pose_form_file(pdb_path)

rt_dicts = loops_to_rt_dict(pose)
if not rt_dicts:
    continue
pdf = pd.DataFrame(rt_dicts)
binner = xb()
pdf["key"] = pdf["rt"].apply(lambda x: xb().get_bin_index(x))
# allowed res is all res in this case, why not
pdf["allowed_res"] = pdf["rt"].apply(lambda x: [*range(1, len(pose.residues))])

# The tables are loaded here from hardcoded paths
key_type, value_type = np.dtype("i8"), np.dtype("i8")
data_name = "exp_no_spin_1ang_15k_ploops_v3"
data_path = "/home/dzorine/phos_binding/pilot_runs/loop_grafting/fragment_tables/ploops_expanded_set_1"
table_dir_path = f"{data_path}/tables"
dict_dir_path = f"{data_path}/dicts"
dict_path = f"{dict_dir_path}/{data_name}.bin"
table_path = f"{table_dir_path}/{data_name}.json"
loop_table, loop_dict = load_table_and_dict(
    table_path, dict_path, key_type, value_type
)
mask = loop_dict.contains(np.array(pdf["key"]))
masked_df = pdf[mask]
if len(masked_df.index) == 0:
    continue
masked_df.loc[:, "e2e_inds"] = loop_dict[np.array(masked_df["key"])]
results = masked_df.loc[~masked_df.index.duplicated(keep="first")]

for col in loop_table:
    results[f"loop_{col}"] = results["e2e_inds"].map(
        lambda index: loop_table[loop_table["index"] == index][col].item()
    )
# results.loc[:, "file"] = results["e2e_inds"].map(
#     lambda index: loop_table[loop_table["index"] == index]["file"].item()
# )
# results.loc[:, "key"] = results["e2e_inds"].map(
#     lambda index: loop_table[loop_table["index"] == index]["key_int"].item()
# )
# results.loc[:, "func_to_bb_start"] = results["e2e_inds"].map(
#     lambda index: loop_table[loop_table["index"] == index]["func_to_bb_start"].item()
# )
# results.loc[:, "e2e"] = results["e2e_inds"].map(
#     lambda index: loop_table[loop_table["index"] == index]["e2e"].item()
# )

not_expanded_allowed_res = results
results = results.allowed_res.apply(pd.Series).merge(results,right_index=True,left_index=True)\
    .melt(id_vars = [*results],value_names="p_res_target")\
    .drop(["allowed_res","variable"],axis = 1)

results["inv_rot_key"] = results.apply(
    lambda row:
            binner.get_bin_index(
                row["func_to_bb_start"]
                @ generate_pose_rt_between_res(pose, row["start_res"], row["p_res_target"])
            )
,
    axis=1,
)
inv_rot_dict_path = "/home/dzorine/projects/phos_binding/pilot_runs/loop_grafting/fragment_tables/inverse_rotamer/dicts/inverse_ptr_exchi7_rotamers_1ang_v1.bin"
inv_rot_table_path = "/home/dzorine/projects/phos_binding/pilot_runs/loop_grafting/fragment_tables/inverse_rotamer/tables/inverse_ptr_exchi7_rotamers_1ang_v1.json"
key_type, value_type = np.dtype("i8"), np.dtype("i8")
inv_rot_table, inv_rot_dict = load_table_and_dict(inv_rot_table_path,inv_rot_dict_path,key_type, value_type,)


results["res_has_key"] = results.apply(
    lambda row:
            inv_rot_dict.contains(np.array(row["inv_rot_key"])),
    axis=1,
)
# results["inv_rot_table_index"] = results.apply(
#     lambda row: inv_rot_dict[np.array(row["inv_rot_key"]],
#     axis=1,
# )
#output to file
