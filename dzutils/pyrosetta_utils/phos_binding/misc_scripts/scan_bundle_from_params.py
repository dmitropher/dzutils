import sys, json

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
        #Hash the params dict for a filename

        #do the grafting
        #dump the pdb and params
        #update the DataFrame
        #dump the dataframe
