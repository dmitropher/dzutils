from itertools import permutations, combinations
import numpy as _np
from pyrosetta.rosetta.core.pose import (
    Pose,
    append_pose_to_pose,
    get_restype_for_pose,
)
from xbin import XformBinner as _xb

from dzutils.pyrosetta_utils.geometry.pose_xforms import (
    get_e2e_xform,
    # get_func_to_end,
    generate_pose_rt_between_res,
)

from dzutils.pyrosetta_utils.geometry.superposition_utilities import (
    super_by_paired_atoms,
)

from dzutils.pyrosetta_utils import (
    residues_with_element,
    hbond_to_residue,
    atom_indices_with_element,
    bonded_atoms,
    or_compose_residue_selectors,
)
from dzutils.pyrosetta_utils.secstruct import parse_structure_from_dssp


def rt_list_hbond_to_res(pose, resnum, sidechain=False, minimal=False):
    """
    Returns a list of hbonds to the given residue number

    Defaults to backbone hbonds only
    """
    res = pose.residue(resnum)
    return [
        generate_pose_rt_between_res(
            pose,
            hbond.don_res(),
            resnum,
            ("N", "N", "CA", "C"),
            (
                res.atom_name(hbond.acc_atm()),
                res.atom_name(hbond.acc_atm()),
                "CA",
                "C",
            ),
            minimal=minimal,
        )
        for hbond in hbond_to_residue(pose, resnum)
        if sidechain or hbond.don_hatm_is_protein_backbone()
    ]


def p_atoms_in_pose(pose):
    return [
        (atom_i, resnum)
        for resnum in residues_with_element(pose, "P")
        for atom_i in atom_indices_with_element(pose.residue(resnum), "P")
    ]


def get_loop_xform_dicts(pose, num_contacts, *args, loop_chain=1, **kwargs):
    """
    returns a dict with xforms from end to end and from phosphate to end

    Dict also contains an xbin key
    """
    e2e = get_e2e_xform(pose.split_by_chain()[loop_chain])
    xb = _xb(*args, **kwargs)
    p_atoms = [
        (
            p_atom_i,
            [
                atom
                for atom in bonded_atoms(
                    pose.residue(resnum), p_atom_i, name=False
                )
                for hbond in hbond_to_residue(pose, resnum, vec=False)
                if hbond.acc_atm() == atom
            ],
            resnum,
        )
        for p_atom_i, resnum in p_atoms_in_pose(pose)
    ]
    xform_dicts = [
        {
            "file": pose.pdb_info().name(),
            "key_int": int(xb.get_bin_index(_np.array(e2e))),
            "e2e": _np.array(e2e),
            # Careful here: maybe add a check to see that res 1 is really the beginning of the loop chain?
            "func_to_bb_start": generate_pose_rt_between_res(
                pose, resnum, 1, (p_atom, p_atom, bonded_1, bonded_2)
            ),
        }
        for p_atom, others, resnum in p_atoms
        if len(others) >= num_contacts
        for bonded_1, bonded_2 in permutations(list(set(others)), 2)
    ]
    return xform_dicts


def loop_res_pairs(pose):
    """
    Returns tuples of all res-res combos within loops
    """
    resnums = [
        i
        for i, is_sel in enumerate(
            or_compose_residue_selectors(
                *[
                    c.generate_residue_selector()
                    for c in parse_structure_from_dssp(pose, "L")
                ]
            ).apply(pose),
            1,
        )
        if is_sel
    ]
    return [(pose.clone(), i, j) for i, j in permutations(resnums, 2) if i > j]


def pair_to_rt_dict(pose, resnum_1, resnum_2, **kwargs):
    """
    Returns a dict with the rt between the given res pair and the pose name

    Can restrict search to a residue index selector string with res_string=
    """

    return {
        "rt": _np.array(
            generate_pose_rt_between_res(pose, resnum_1, resnum_2, **kwargs)
        ),
        "name": pose.pdb_info().name(),
        "start_res": resnum_1,
        "end_res": resnum_2,
    }


def loops_to_rt_dict(pose, plus=0, minus=0):
    """
    Wrapper to hopefully conserve memory

    includes some modifications to loop selection to help with bad hits
    """
    return [
        pair_to_rt_dict(pose, i, j)
        for loop in parse_structure_from_dssp(pose, "L")
        for start, end in [
            (
                loop.resnum_list(upstream=minus, downstream=plus)[0],
                loop.resnum_list(upstream=minus, downstream=plus)[-1],
            )
        ]
        for i, j in combinations(
            loop.resnum_list(upstream=minus, downstream=plus), 2
        )
        if i > (end + start) / 2 and j > (end + start) / 2
    ]


def phospho_residue_inverse_rotamer_rts(residue):
    """
    Returns the RTs from the phos group to bb

    Generates the RT from all possible stubs generated by the following scheme

    (Phosphorus atom, Phosphorus Atom, p-bound-atom1,p-bound-atom2)
    """
    phos_atoms = atom_indices_with_element(residue, "P")
    possible_rt_bases = [
        (p_atom_index, p_atom_index, a, b)
        for p_atom_index in phos_atoms
        for (a, b) in permutations(
            bonded_atoms(residue, p_atom_index, name=False), 2
        )
    ]
    pose = Pose()
    pose.append_residue_by_bond(residue)
    return [
        generate_pose_rt_between_res(pose.clone(), 1, 1, base)
        for base in possible_rt_bases
    ]


def num_bb_contacts(pose, resnum, atom_i):
    contacts = sum(
        (
            hbond.don_hatm_is_protein_backbone()
            and hbond.acc_atm() in bonded_atoms(pose.residue(resnum), atom_i)
        )
        for hbond in hbond_to_residue(pose, resnum)
    )
    print(contacts)
    return contacts


def replace_p_res_with_phosphate(pose, min_contacts=0):
    """
    returns a copy of the pose with every phosphoresidue converted to phosphate

    PO4 must be a residue in your rosetta instance
    """

    p_res_list = p_atoms_in_pose(pose)
    phosphates = list()
    for atom_i, resnum in p_res_list:
        bb_contacts = num_bb_contacts(pose, resnum, atom_i)
        if bb_contacts < min_contacts:
            print(
                f"res {resnum} at atom {atom_i} has only {bb_contacts} contacts, less than {min_contacts}"
            )
            continue
        p = Pose()
        restype = get_restype_for_pose(p, "PO4")
        res = _pyrosetta.rosetta.core.conformation.Residue(restype, True)
        p.append_residue_by_jump(res, 1)
        p_atom = atom_indices_with_element(res, "P")[0]

        atom_tuples = list(
            zip(
                [p_atom, *bonded_atoms(res, p_atom)[:2]],
                [atom_i, *bonded_atoms(pose.residue(resnum), atom_i)[:2]],
            )
        )
        # print(atom_tuples)
        phosphates.append(
            super_by_paired_atoms(p, pose, 1, resnum, *atom_tuples)
        )
    res_to_remove = list(set([resnum for atom_i, resnum in p_res_list]))
    # print(res_to_remove)
    newp = pose.clone()
    for i, resnum in enumerate(res_to_remove):
        newp.delete_residue_slow(resnum - i)
    for phos in phosphates:
        append_pose_to_pose(newp, phos, True)
    return newp
