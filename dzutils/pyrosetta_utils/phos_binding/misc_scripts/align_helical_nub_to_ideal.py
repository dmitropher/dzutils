import numpy as np

from dzutils.pyrosetta_utils.geometry.parametric import ca_array
from dzutils.pyrosetta_utils.geometry import (
    stub_to_homog,
    homog_from_four_points,
    np_homog_to_rosetta_rotation_translation,
)
from pyrosetta.rosetta.core.kinematics import Stub, RT
from pyrosetta.rosetta.numeric import xyzVector_double_t as rosetta_vector
from pyrosetta.rosetta.core.kinematics import Stub, RT
from pyrosetta.rosetta.protocols.toolbox.pose_manipulation import (
    rigid_body_move,
)
from pyrosetta.rosetta.core.select.residue_selector import TrueResidueSelector


def align_frag_to_ideal_helix(frag, res_1, res_2, res_3, **params):
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
    targ_stub = Stub()
    targ_stub.from_four_points(
        *[frag.residue(num).xyz("CA") for num in [res_2, res_1, res_2, res_3]]
    )
    targ_homog = stub_to_homog(targ_stub)
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


def main():
    """
    """


if __name__ == "__main__":
    main()
