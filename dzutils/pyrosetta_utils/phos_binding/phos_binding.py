from .util import hbond_to_residue as _hbond_to_residue

# from ..pyrosetta_utils.util import residues_by_name as _residues_by_name
from ..geometry.pose_xforms import *


def rt_list_hbond_to_res(pose, resnum, sidechain=False):
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
        )
        for hbond in _hbond_to_residue(pose, resnum)
        if sidechain or hbond.don_hatm_is_protein_backbone()
    ]
