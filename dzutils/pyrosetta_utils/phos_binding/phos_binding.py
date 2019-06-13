from copy import deepcopy as _dc

import numpy as _np
from xbin import XformBinner as _xb

from dzutils.pyrosetta_utils.geometry.pose_xforms import (
    get_e2e_xform,
    get_func_to_end,
)
from dzutils.pyrosetta_utils import (
    residues_with_element,
    hbond_to_residue,
    atom_indices_with_element,
    bonded_atoms,
)
from .util import hbond_to_residue

# from ..pyrosetta_utils.util import residues_by_name as _residues_by_name
from dzutils.pyrosetta_utils.geometry.pose_xforms import (
    generate_pose_rt_between_res,
)


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


def p_bound_atoms(pose):
    return [
        (atom_i, resnum)
        for resnum in residues_with_element(pose, "P")
        for atom_i in atom_indices_with_element(pose.residue(resnum), "P")
    ]


def get_loop_xform_dicts(pose, num_contacts, *args, **kwargs):
    """
    returns a dict with xforms from end to end and from phosphate to end

    Dict also contains an xbin key
    """
    e2e = get_e2e_xform(pose)
    xb = _xb(*args, **kwargs)
    xform_dicts = [
        {
            "file": pose.pdb_info().name(),
            "key_int": int(xb.get_bin_index(_np.array(e2e))),
            "e2e": _dc(e2e),
            "func_to_bb_start": get_func_to_end(
                pose,
                resnum,
                (p_atom, p_atom, others[0], others[1]),
                begin=True,
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
            for atom_i, resnum in p_bound_atoms(pose)
        ]
        if len(others) >= num_contacts
    ]
    return xform_dicts
