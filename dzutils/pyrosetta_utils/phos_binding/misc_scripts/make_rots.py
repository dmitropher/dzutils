import pyrosetta
import pyrosetta.rosetta as pyr

from dzutils.pyrosetta_utils.geometry.superposition_utilities import (
    super_by_residues,
)
from dzutils.pyrosetta_utils import bonded_atoms, atom_indices_with_element

pyrosetta.init(
    """-out:level 100
    -extra_res_fa /home/dzorine/phos_binding/p_compounds/residues/PTR.params
    -packing:ex1
    -packing:ex2
"""
    """
    -packing:ex3
    -packing:ex4

"""
)
# pyrosetta.init(
#     """-out:level 100
#     -extra_res_fa /home/dzorine/phos_binding/p_compounds/residues/PTR.params
#     -packing:ex1
#     -packing:ex2
#     -packing:ex3
#     -packing:ex4
#     -packing:ex1:level 7
#     -packing:ex2:level 7
#     -packing:ex3:level 7
#     -packing:ex4:level 7
#
# """
#     # -extra_res_fa /home/dzorine/phos_binding/p_compounds/residues/PTR.params
# )
# pose = pyrosetta.pose_from_sequence("YYY")
# mutate_res = pyrosetta.rosetta.protocols.simple_moves.MutateResidue()
# phy_index = 2
# mut_index = phy_index
# mut_name = "PTR"
# mutate_res.set_target(str(mut_index))
# mutate_res.set_res_name(mut_name)
# mutate_res.apply(pose)
def one_res_pose(resname):

    pose = pyr.core.pose.Pose()
    res = pyr.core.conformation.Residue(
        pyr.core.pose.get_restype_for_pose(pose, resname), True
    )
    pose.append_residue_by_bond(res)
    return pose


def rotamers_for_resname(resname):
    """
    Creates the backbone independent rotamer set for resname
    """
    return [
        rot
        for rot in pyr.core.pack.rotamer_set.bb_independent_rotamers(
            one_res_pose(resname).residue(1).type()
        )
    ]


def rot_poses(resname, *super_atoms):
    """
    returns a list of poses where each is a different rotamer of restype

    leaving super_atoms blank gives forward rotamers, giving super_atoms
    aligns all poses to those atoms
    """
    rots = rotamers_for_resname(resname)
    pose = pyr.core.pose.Pose()
    pose.append_residue_by_bond(rots[0])
    poses = [pose]
    for i, rot in enumerate(rots[1:]):
        pose = pyr.core.pose.Pose()
        pose.append_residue_by_bond(rot)
        super_by_residues(pose, poses[i], 1, 1, *super_atoms)
        poses.append(pose)
    return poses


def main():
    p_res = pyr.core.conformation.Residue(
        pyr.core.pose.get_restype_for_pose(pyr.core.pose.Pose(), "PTR"), True
    )
    p_atoms = [
        "P",
        *bonded_atoms(p_res, *atom_indices_with_element(p_res, "P")),
    ]

    poses = rot_poses("PTR", *p_atoms)
    out_pose = poses[0]
    for pose in poses[1:]:
        pyr.core.pose.append_pose_to_pose(out_pose, pose, True)

    out_pose.dump_pdb("/home/dzorine/temp/ptr_rots.pdb")


if __name__ == "__main__":
    main()
