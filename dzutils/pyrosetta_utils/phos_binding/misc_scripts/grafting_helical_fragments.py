import numpy as np

from dzutils.pyrosetta_utils.geometry.parametric import ca_array
from dzutils.pyrosetta_utils.geometry import (
    homog_from_four_points,
    np_homog_to_rosetta_rotation_translation,
    homog_from_3_CA,
)
from pyrosetta.rosetta.numeric import xyzVector_double_t as rosetta_vector
from pyrosetta.rosetta.protocols.toolbox.pose_manipulation import (
    rigid_body_move,
)
from pyrosetta.rosetta.core.select.residue_selector import TrueResidueSelector


def align_frag_to_ideal_helix(frag, align_res, **params):
    """
    Align and return frag by stub to ideal helix described by the params given
    """
    r0 = params["r0"]
    omega0 = params["omega0"]
    omega1 = params["omega1"]
    phi0 = params["phi0"]
    phi1 = params["phi1"]
    delta_z = params["delta_z"]
    invert = params["invert"]
    ideal_ca = ca_array(r0, omega0, omega1, phi0, phi1, delta_z, 3, invert)
    ca_homog = homog_from_four_points(ideal_ca[1], *ideal_ca)
    # targ_stub = Stub()
    # targ_stub.from_four_points(
    #     *[frag.residue(num).xyz("CA") for num in [res_2, res_1, res_2, res_3]]
    # )
    targ_homog = homog_from_3_CA(frag, align_res)
    super_xform = ca_homog @ np.linalg.inv(targ_homog)
    rotation, translation = np_homog_to_rosetta_rotation_translation(
        super_xform
    )
    rigid_body_move(
        rotation,
        translation,
        frag,
        TrueResidueSelector().apply(frag),
        rosetta_vector(0, 0, 0),
    )
    return frag


def overlay_poses_by_3_CA(mob_pose, targ_pose, mob_res, targ_res):
    """
    Returns mob_pose res i - i+2 overlaid onto the respective res of targ_pose
    """
    # "Stub" homog xform for the mobile pose
    mob_homog = homog_from_3_CA(mob_pose, mob_res)
    # "Stub" homog xform for the target pose
    targ_homog = homog_from_3_CA(targ_pose, targ_res)
    super_xform = targ_homog @ np.linalg.inv(mob_homog)
    rotation, translation = np_homog_to_rosetta_rotation_translation(
        super_xform
    )
    rigid_body_move(
        rotation,
        translation,
        mob_pose,
        TrueResidueSelector().apply(mob_pose),
        rosetta_vector(0, 0, 0),
    )


def main():
    """
    """


if __name__ == "__main__":
    main()
