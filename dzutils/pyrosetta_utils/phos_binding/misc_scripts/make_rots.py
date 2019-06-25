import pyrosetta
import pyrosetta.rosetta as pyr

import dzutils.pyrosetta_utils.geometry.superposition_utilities as su

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
pose = pyrosetta.pose_from_sequence("YYY")
mutate_res = pyrosetta.rosetta.protocols.simple_moves.MutateResidue()
phy_index = 2
mut_index = phy_index
mut_name = "PTR"
mutate_res.set_target(str(mut_index))
mutate_res.set_res_name(mut_name)
mutate_res.apply(pose)
rt = pose.residue(2).type()
print(rt)
rots = [rot for rot in pyr.core.pack.rotamer_set.bb_independent_rotamers(rt)]
out_pose = pyr.core.pose.Pose()
rot_1 = rots[0]
out_pose.append_residue_by_bond(rot_1)
template = out_pose.clone()
for i, rot in enumerate(rots[1:], 1):
    temp = pyr.core.pose.Pose()
    temp.append_residue_by_bond(rot)
    su.super_by_residues(temp, template, 1, 1, *("O1P", "O2P", "P", "O3P"))
    out_pose.append_residue_by_jump(temp.residue(1), i)
out_pose.dump_pdb(
    "/home/dzorine/phos_binding/pilot_runs/loop_grafting/initial_testing/toy_scan_test/ptr_rots_fix/ptr_rots.pdb"
)
