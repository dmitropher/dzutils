# Attempts to close from each chain to each other chain
# For each chain pair, does serial deletions from the ends to increase
# the number of samples

import pyrosetta
import sys

from dzutils.sutils import read_flag_file

from dzutils.pyrosetta_utils.phos_binding.misc_scripts.remove_loops_rechain import (
    ss_to_chains,
)
from dzutils.pyrosetta_utils.phos_binding.misc_scripts.exhaustive_single_loop_insertion import (
    exhaustive_single_loop_insertion,
)


flagsFile = "/home/dzorine/phos_binding/pilot_runs/loop_grafting/initial_testing/misc_files/cluster_altered.flags"
flags = read_flag_file(flagsFile)
flags_str = " ".join(flags.replace("\n", " ").split())
pyrosetta.init(flags_str)

pose = pyrosetta.pose_from_file(sys.argv[1])
out_dir = sys.argv[2]
# splits by helices
ss_only_chains = ss_to_chains(pose, "H")
pose_name = pose.pdb_info().name().split(".pdb")[0].split("/")[-1]
suffix = "_helix_only.pdb"
ss_only_name = f"{out_dir}/source_files/{pose_name}{suffix}"

# mind the hardcode
ss_only_chains.dump_pdb(ss_only_name)
ss_only_chains = pyrosetta.pose_from_file(ss_only_name)
for looped_pose in exhaustive_single_loop_insertion(
    ss_only_chains, out_dir, 5
):
    """
    """
