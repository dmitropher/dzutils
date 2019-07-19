import pyrosetta
from itertools import permutations
import numpy as np
from xbin import XformBinner as xb

from dzutils.pyrosetta_utils.geometry.pose_xforms import (
    generate_pose_rt_between_res,
    PoseStubArray,
)
from dzutils.pyrosetta_utils.geometry.homog import homog_from_residue
from dzutils.pyrosetta_utils import (
    bonded_atoms,
    atom_indices_with_element,
    build_hbond_set,
)

# from dzutils.pyrosetta_utils.phos_binding import num_bb_contacts
from dzutils.pyrosetta_utils.geometry.homog import (
    np_homog_to_rosetta_rotation_translation,
)


def build_harmonic_pair_constraint(res_1, atom_1, res_2, atom_2, value, sd):
    return pyrosetta.rosetta.core.scoring.constraints.AtomPairConstraint(
        pyrosetta.rosetta.core.id.AtomID(atom_1, res_1),
        pyrosetta.rosetta.core.id.AtomID(atom_2, res_2),
        pyrosetta.rosetta.core.scoring.func.HarmonicFunc(value, sd),
    )


def build_hbond_cst(
    grafted_pose,
    rot_pose,
    rotamer_resnum,
    target_rotamer_resnum,
    func_to_bb_xform,
    grafted_pose_loop_start_resnum,
    rotamer_alignment_atoms,
    acceptor_atoms,
):
    """
    Build atom pair constraints to preserve the native hbonds
    """
    start_res_stub = homog_from_residue(
        grafted_pose.residue(grafted_pose_loop_start_resnum)
    )
    rot_res_stub = homog_from_residue(
        rot_pose.residue(rotamer_resnum), *rotamer_alignment_atoms
    )
    phos_stub = start_res_stub @ np.linalg.inv(func_to_bb_xform)
    rot_to_phos_pos = phos_stub @ np.linalg.inv(rot_res_stub)

    rotation, translation = np_homog_to_rosetta_rotation_translation(
        rot_to_phos_pos
    )
    pyrosetta.rosetta.protocols.toolbox.pose_manipulation.rigid_body_move(
        rotation,
        translation,
        rot_pose,
        pyrosetta.rosetta.core.select.residue_selector.TrueResidueSelector().apply(
            rot_pose
        ),
        pyrosetta.rosetta.numeric.xyzVector_double_t(0, 0, 0),
    )
    check_pose = grafted_pose.clone()
    pyrosetta.rosetta.core.pose.append_pose_to_pose(check_pose, rot_pose, True)
    test_file_path = (
        f"{check_pose.pdb_info().name().split('.pdb')[0]}_check.pdb"
    )
    print(test_file_path)
    print("building constraints for:")
    print(f"{check_pose.residue(target_rotamer_resnum)}")
    check_pose.dump_pdb(test_file_path)
    check_index = len(check_pose.residues)
    hbond_set = build_hbond_set(check_pose)
    constraints = [
        cst
        for hbond in hbond_set.residue_hbonds(check_index)
        if hbond.don_hatm_is_protein_backbone()
        if hbond.acc_atm() in acceptor_atoms
        for cst in (
            build_harmonic_pair_constraint(
                target_rotamer_resnum,
                hbond.acc_atm(),
                hbond.don_res(),
                hbond.don_hatm(),
                value=hbond.get_HAdist(check_pose),
                sd=0.125,
            ),
            # build_harmonic_pair_constraint(
            #     target_rotamer_resnum,
            #     hbond.acc_atm(),
            #     hbond.don_res(),
            #     check_pose.residue(hbond.don_res()).atom_base(
            #         hbond.don_hatm()
            #     ),
            #     value=2.5,
            #     sd=0.25,
            # ),
        )
    ]
    return constraints


def rebuild_inv_rotamer_alignment_atoms(
    check_key, allowed_atoms, pose, resnum, allow_redundant=False
):
    """
    returns the set of atoms to rebuild the stub used for the inverse rotamer

    if allow_redundant is False (default) only returns the first viable set

    returns None if no set is found
    """
    bb_stub = homog_from_residue(pose.residue(resnum))
    align_set = [
        (atom_1, atom_1, atom_2, atom_3)
        for atom_1, atom_2, atom_3 in permutations(allowed_atoms, 3)
        if xb().get_bin_index(
            np.linalg.inv(
                homog_from_residue(
                    pose.residue(resnum), atom_1, atom_1, atom_2, atom_3
                )
            )
            @ bb_stub
        )
        == check_key
    ]
    if allow_redundant:
        return align_set
    return align_set[0] if align_set else None


def get_phosphate_atoms(res):
    """
    returns phosphate atom indices for the given residue

    in this context, phosphate means P and any atom bound to p
    """
    p_atoms = atom_indices_with_element(res, "P")
    p_bonded_atoms = [atom for p in p_atoms for atom in bonded_atoms(res, p)]
    return (*p_atoms, *p_bonded_atoms)
