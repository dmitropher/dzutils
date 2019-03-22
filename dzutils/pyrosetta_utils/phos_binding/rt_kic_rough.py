import homog
import pyrosetta
import numpy
from dzutils.pyrosetta_utils.geometry.homog import *
from dzutils.pyrosetta_utils.geometry.rt_utils import stub_from_residue

pyrosetta.init()

# Saves the RT from the two special residues in the given pose
dup_pity = pyrosetta.pose_from_file(
    "/home/dzorine/phos_binding/pilot_runs/loop_grafting/native_loops/interaction.pdb"
)
spec_gly = dup_pity.residue(2).clone()
phi_psi = (dup_pity.phi(2), dup_pity.psi(2))
print(phi_psi)
ptr_res = dup_pity.residue(4).clone()
gly_stub = stub_from_residue(spec_gly, "CA", "C", "CA", "N")
ptr_stub = stub_from_residue(ptr_res, "P", "O2P", "P", "O3P")
res_rt_homog = rt_to_homog(
    pyrosetta.rosetta.core.kinematics.RT(gly_stub, ptr_stub)
)

pity = pyrosetta.pose_from_file(
    "/home/dzorine/phos_binding/pilot_runs/loop_grafting/native_loops/1pty.pdb"
)
pity.delete_residue_range_slow(298, 299)


break_chain = (
    pyrosetta.rosetta.protocols.protein_interface_design.movers.AddChainBreak()
)
break_chain.change_conformation(True)
break_chain.change_foldtree(True)
break_chain.resnum("115")
break_chain.apply(pity)

break_chain = (
    pyrosetta.rosetta.protocols.protein_interface_design.movers.AddChainBreak()
)
break_chain.change_conformation(True)
break_chain.change_foldtree(True)
break_chain.resnum("119")
break_chain.apply(pity)
chains = pity.split_by_chain()
n_term = chains[1]
frag = chains[2]
c_term = chains[3]


frag.set_phi(2, phi_psi[0])
frag.set_psi(2, phi_psi[1])
# pyrosetta.rosetta.core.conformation.add_upper_terminus_type_to_conformation_residue(frag.conformation(),3)
frag.dump_pdb(
    "/home/dzorine/phos_binding/pilot_runs/loop_grafting/native_loops/frag.pdb"
)
print("frag phi psi: ", frag.phi(2), frag.psi(2))
newp = n_term
newp.append_pose_by_jump(frag, 1)
newp.append_pose_by_jump(c_term, 2)

pyrosetta.rosetta.core.conformation.remove_variant_type_from_conformation_residue(
    newp.conformation(),
    pyrosetta.rosetta.core.chemical.VariantType.CUTPOINT_UPPER,
    116,
)

pyrosetta.rosetta.core.conformation.remove_variant_type_from_conformation_residue(
    newp.conformation(),
    pyrosetta.rosetta.core.chemical.VariantType.CUTPOINT_LOWER,
    115,
)

pyrosetta.rosetta.core.conformation.remove_lower_terminus_type_from_conformation_residue(
    newp.conformation(), 116
)

pyrosetta.rosetta.core.conformation.remove_upper_terminus_type_from_conformation_residue(
    newp.conformation(), 115
)

dbond = pyrosetta.rosetta.protocols.cyclic_peptide.DeclareBond()
dbond.set(115, "C", 116, "N", False)
dbond.apply(pity)

print("newp phi psi: ", newp.phi(117), newp.psi(117))
newp.dump_pdb(
    "/home/dzorine/phos_binding/pilot_runs/loop_grafting/native_loops/tweak_phi_psi.pdb"
)

genkic = (
    pyrosetta.rosetta.protocols.generalized_kinematic_closure.GeneralizedKIC()
)
loop_res = range(112, 122)
for i in loop_res:
    genkic.add_loop_residue(i)
genkic.add_perturber(
    pyrosetta.rosetta.protocols.generalized_kinematic_closure.perturber.randomize_alpha_backbone_by_rama
)
for i in [j for j in loop_res if j != 117]:
    genkic.add_residue_to_perturber_residue_list(1, i)
genkic.add_filter("loop_bump_check")
genkic.close_bond(
    115,
    "C",
    116,
    "N",
    115,
    "CA",
    116,
    "CA",
    1.328685,
    116.19999299975001,
    121.6999970001592,
    180.0,
    False,
    False,
)
genkic.set_pivot_atoms(112, "CA", 114, "CA", 121, "CA")
genkic.set_selector_type("lowest_energy_selector")
genkic.set_selector_scorefunction(
    pyrosetta.rosetta.core.scoring.get_score_function()
)
genkic.set_closure_attempts(100)

genkic.apply(newp)
newp.dump_pdb(
    "/home/dzorine/phos_binding/pilot_runs/loop_grafting/native_loops/kicced.pdb"
)
pyrosetta.rosetta.core.conformation.remove_lower_terminus_type_from_conformation_residue(
    newp.conformation(), 120
)
pyrosetta.rosetta.core.conformation.remove_upper_terminus_type_from_conformation_residue(
    newp.conformation(), 119
)
pyrosetta.rosetta.core.conformation.remove_variant_type_from_conformation_residue(
    newp.conformation(),
    pyrosetta.rosetta.core.chemical.VariantType.CUTPOINT_UPPER,
    120,
)
pyrosetta.rosetta.core.conformation.remove_variant_type_from_conformation_residue(
    newp.conformation(),
    pyrosetta.rosetta.core.chemical.VariantType.CUTPOINT_LOWER,
    119,
)
dbond = pyrosetta.rosetta.protocols.cyclic_peptide.DeclareBond()
dbond.set(119, "C", 120, "N", False)
dbond.apply(pity)
genkic = (
    pyrosetta.rosetta.protocols.generalized_kinematic_closure.GeneralizedKIC()
)
loop_res = range(112, 123)
for i in loop_res:
    genkic.add_loop_residue(i)
genkic.add_perturber(
    pyrosetta.rosetta.protocols.generalized_kinematic_closure.perturber.randomize_backbone_by_rama_prepro
)
for i in [j for j in loop_res if j != 117]:
    genkic.add_residue_to_perturber_residue_list(1, i)
genkic.add_filter("loop_bump_check")
genkic.close_bond(
    119,
    "C",
    120,
    "N",
    119,
    "CA",
    120,
    "CA",
    1.328685,
    121.6999970001592,
    116.19999299975001,
    180.0,
    False,
    False,
)
genkic.set_pivot_atoms(112, "CA", 120, "CA", 122, "CA")
genkic.set_selector_type("lowest_energy_selector")
genkic.set_selector_scorefunction(
    pyrosetta.rosetta.core.scoring.get_score_function()
)
genkic.set_closure_attempts(300)

genkic.apply(newp)
ft = pyrosetta.rosetta.core.kinematics.FoldTree()
ft.add_edge(1, 297, -1)
ft.check_fold_tree()
newp.fold_tree(ft)
newp.dump_pdb(
    "/home/dzorine/phos_binding/pilot_runs/loop_grafting/native_loops/kicced.pdb"
)
# makes a toy mobile and target
mob = pyrosetta.pose_from_file(
    "/home/dzorine/phos_binding/pilot_runs/loop_grafting/native_loops/toy_mob.pdb"
)
targ = pyrosetta.pose_from_file(
    "/home/dzorine/phos_binding/pilot_runs/loop_grafting/native_loops/toy_target.pdb"
)
# makes homogenous transforms representing stubs for the mobile and target
mob_res = mob.residue(1)
mob_stub = stub_from_residue(mob.residue(1), "P", "O2P", "P", "O3P")
mob_stub_homog = stub_to_homog(mob_stub)
targ_res = newp.residue(181).clone()
print(newp.phi(181), " ", newp.psi(181))
targ_stub = stub_from_residue(targ_res, "CA", "C", "CA", "N")
targ_stub_homog = stub_to_homog(targ_stub)

# computes the inverse (Not the transpose for homogenous transforms!)
# of the mobile
c_inv = invert_homog(mob_stub_homog)

# Makes the homogenous transform from the mobile to the desired location
transformer = targ_stub_homog @ res_rt_homog @ c_inv

# applies the transform atomwise
for i, atom in enumerate(mob_res.atoms(), 1):
    xyz = atom.xyz()
    xyz_array = numpy.array([[xyz.x], [xyz.y], [xyz.z], [1]])
    new_coords = transformer @ xyz_array
    new_xyz = pyrosetta.rosetta.numeric.xyzVector_double_t(
        x_a=new_coords[0], y_a=new_coords[1], z_a=new_coords[2]
    )
    mob_res.set_xyz(i, new_xyz)

# appends the residue and dumps the pdb
mob = pyrosetta.rosetta.core.pose.Pose()
mob.append_residue_by_jump(mob_res, 0)
mob.dump_pdb(
    "/home/dzorine/phos_binding/pilot_runs/loop_grafting/native_loops/mob.pdb"
)
