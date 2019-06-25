from os.path import basename
import sys

import pyrosetta

from dzutils.pyrosetta_utils.chain_utils import link_poses
from dzutils.pyrosetta_utils.secstruct.structure import split_ss_to_subposes
from dzutils.sutils import read_flag_file


def ss_to_chains(pose, *dssp_types):
    """
    returns a pose that has been split by chain on secondary structure
    """
    subposes = split_ss_to_subposes(pose, *dssp_types)
    return link_poses(*subposes, rechain=True)


def main():
    flagsFile = "/home/dzorine/phos_binding/pilot_runs/loop_grafting/initial_testing/misc_files/cluster_altered.flags"
    flags = read_flag_file(flagsFile)
    flags_str = " ".join(flags.replace("\n", " ").split())
    pyrosetta.init(flags_str)

    pose = pyrosetta.pose_from_file(sys.argv[1])
    # splits by helices
    ss_only_chains = ss_to_chains(pose, "H", "E")
    pose_basename = basename(pose.pdb_info().name()).split(".pdb")[0]
    suffix = "_secstruct_only.pdb"
    ss_only_chains.dump_pdb(f"./{pose_basename}{suffix}")


if __name__ == "__main__":
    main()
