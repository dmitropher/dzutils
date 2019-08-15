import sys, json,os

import pyrosetta

from dzutils.pyrosetta_utils.phos_binding.parametric import (
    helix_bundle_maker_wrapper,
    PyBundleGridSampler,
)
import pandas as pd
from dzutils.pyrosetta_utils.phos_binding.misc_scripts.double_scan_pbinders import (
    scan_for_n_term_helical_grafts,
    scan_for_inv_rot,
)

from dzutils.pyrosetta_utils.phos_bindingo import (graft_and_dump_pdb, replace_res_from_pose)
def process_secondary_results(
    secondary_results_table,
    pose,
    secondary_hit_out_dir,
    source_file_dir,
    table_label="",
    pdb_label="",
):
    """
    """
    print(secondary_results_table)
    if not os.path.isdir(secondary_hit_out_dir):
        os.makedirs(secondary_hit_out_dir)
    secondary_write_path = f"""{secondary_hit_out_dir
                               }/secondary_double_ploop_hits{
                               table_label
                               }.json"""

    # assert not os.path.exists(
    #     secondary_write_path
    # ), f"secondary hit data would overwrite file at: {secondary_write_path}"
    # source_file_dir = f"{secondary_hit_out_dir}/source_files/"
    # if not os.path.isdir(source_file_dir):
    #     os.makedirs(source_file_dir)
    secondary_grafted_loops_dir = (
        f"{secondary_hit_out_dir}/grafted_loops_with_p_res/"
    )
    if not os.path.isdir(secondary_grafted_loops_dir):
        os.makedirs(secondary_grafted_loops_dir)
    if not os.path.isdir(source_file_dir):
        os.makedirs(source_file_dir)
    # Add a grafted file entry with phos rotamer to each row, dump the pdb
    print("grafting and adding residue")
    secondary_results_table[
        "graft_and_pres_path"
    ] = secondary_results_table.apply(
        lambda row: graft_and_dump_pdb(
            replace_res_from_pose(
                pose.clone(),
                pyrosetta.pose_from_file(row["inv_rot_file"]),
                row["p_res_target"],
                1,
            ),
            pyrosetta.pose_from_file(row["loop_file"]).split_by_chain()[1],
            row["start_res"],
            row["end_res"],
            1,
            f"""{secondary_grafted_loops_dir
                }/{
                row["name"].split("/")[-1].split(".pdb")[0]
                }{pdb_label}_ploop_index_{
                row["loop_index"]
                }_pres_{
                row["p_res_target"]
                }.pdb""",
        ),
        axis=1,
    )
    secondary_results_table.apply(
        lambda row: pyrosetta.pose_from_file(row["inv_rot_file"]).dump_pdb(
            f'{source_file_dir}/{row["inv_rot_file"].split("/")[-1]}'
        ),
        axis=1,
    )
    secondary_results_table.to_json(secondary_write_path)

table_path = "/home/dzorine/phos_binding/pilot_runs/loop_grafting/fragment_tables/ploop_fragment_set/tables/1ang_3_contact_ploop_set_4_v3.json"
dict_path = "/home/dzorine/phos_binding/pilot_runs/loop_grafting/fragment_tables/ploop_fragment_set/dicts/1ang_3_contact_ploop_set_4_v3.bin"

pyrosetta.init()
param_json = json.loads(sys.argv[1])
name_hash = hash(param_json)
param_dict = dict(param_json)

params_for_helices = [param_dict[i] for i in range (1,param_dict["num_helices"]+1)]

#Use the wrapper to create the MakeBundle Mover
bundle_maker = helix_bundle_maker_wrapper(param_dict["helix_length"],*params_for_helices,degrees=True)

#Empty pose to helical bundle!
pose = pyrosetta.rosetta.core.pose.Pose()
bundle_maker.apply(pose)
inf = pyrosetta.rosetta.core.pose.PDBInfo(pose)
inf.name(f"param_bundle_{name_hash}.pdb")
pose.pdb_info(inf)

if len(pose.residues):
    results = scan_for_n_term_helical_grafts(pose,table_path,dict_path)
    results["allowed_res"] = results.apply(lambda x: [*range(1,len(pose.residues)+1)],axis=1)
    results = results.rename(index=str, columns={"frag_feature_to_start":"loop_func_to_bb_start"})
    sec_results = scan_for_inv_rot(pose,results,table_path,dict_path)
    if sec_results:
        print ("success!")

        #do the grafting
        #dump the pdb and params
        #update the DataFrame
        #dump the dataframe
