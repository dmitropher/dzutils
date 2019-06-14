from itertools import permutations
import numpy as _np
from xbin import XformBinner as _xb

from dzutils.pyrosetta_utils.geometry.pose_xforms import (
    get_e2e_xform,
    # get_func_to_end,
    generate_pose_rt_between_res,
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
    xform_dicts = [
        {
            "file": pose.pdb_info().name(),
            "key_int": int(xb.get_bin_index(_np.array(e2e))),
            "e2e": _np.array(e2e),
            # Careful here: maybe add a check to see that res 1 is really the beginning of the loop chain?
            "func_to_bb_start": generate_pose_rt_between_res(
                pose, resnum, 1, (p_atom, p_atom, others[0], others[1])
            ),
        }
        for p_atom, others, resnum in [
            (
                atom_i,
                [
                    atom
                    for atom in bonded_atoms(
                        pose.residue(resnum), atom_i, name=False
                    )
                    for hbond in hbond_to_residue(pose, resnum, vec=False)
                    if hbond.acc_atm() == atom
                ],
                resnum,
            )
            for atom_i, resnum in p_atoms_in_pose(pose)
        ]
        if len(others) >= num_contacts
    ]
    return xform_dicts


def loop_res_pairs(pose):
    """
    Returns tuples of all res-res combos

    Can restrict search to a residue index selector string
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
    return [(pose.clone(), i, j) for i, j in permutations(resnums, 2)]


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
