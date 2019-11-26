import numpy as _numpy
import math as _math
from dzutils.pyrosetta_utils.geometry.rt_utils import (
    stub_from_residue as _stub_from_residue,
)

from pyrosetta.rosetta.numeric import xyzVector_double_t, xyzMatrix_double_t
from pyrosetta.rosetta.core.kinematics import Stub


def homog_from_four_points(center, a, b, c):
    """
    Returns a homogenous xform analogous to a rosetta stub
    """
    p1 = a - b
    e1 = p1 / _numpy.linalg.norm(p1)
    p3 = _numpy.cross(e1, (c - b))
    e3 = p3 / _numpy.linalg.norm(p3)
    e2 = _numpy.cross(e3, e1)
    homog = _numpy.empty([4, 4])
    homog[..., :3, :3] = _numpy.transpose(_numpy.array([e1, e2, e3]))
    homog[:3, 3] = center
    homog[..., 3, :] = [0, 0, 0, 1]
    return homog


def np_rot_trans_to_homog(rotation, translation):
    """
    Basic wrapper to turn 3x3 and 1x3 to 4x4 homog xform
    """
    homog = _numpy.empty([4, 4])
    homog[..., :3, :3] = rotation
    homog[:3, 3] = translation
    homog[..., 3, :] = [0, 0, 0, 1]
    return homog


def rotation_translation_to_homog(rotation, translation):
    """
    Takes a rotation matrix and a translation vector and returns a h xform
    """
    return _numpy.array(
        [
            [rotation.xx, rotation.xy, rotation.xz, translation.x],
            [rotation.yx, rotation.yy, rotation.yz, translation.y],
            [rotation.zx, rotation.zy, rotation.zz, translation.z],
            [0, 0, 0, 1],
        ]
    )


def np_homog_to_rosetta_rotation_translation(xform):
    """
    takes a numpy homogenous transform and converts it to Rosetta (rot,trans)
    """
    rotation = xform[..., :3, :3]
    translation = xform[:3, 3]
    rosetta_rotation = xyzMatrix_double_t()
    (
        rosetta_rotation.xx,
        rosetta_rotation.xy,
        rosetta_rotation.xz,
        rosetta_rotation.yx,
        rosetta_rotation.yy,
        rosetta_rotation.yz,
        rosetta_rotation.zx,
        rosetta_rotation.zy,
        rosetta_rotation.zz,
    ) = rotation.flatten()
    rosetta_translation = xyzVector_double_t()
    (
        rosetta_translation.x,
        rosetta_translation.y,
        rosetta_translation.z,
    ) = translation
    return rosetta_rotation, rosetta_translation


def stub_to_homog(stub):
    """
    Takes a rosetta stub and returns a h xform
    """
    return rotation_translation_to_homog(stub.M, stub.v)


def homog_from_residue(
    residue, center_atom="CA", atom1="N", atom2="CA", atom3="C"
):
    rstub = Stub()
    rstub.from_four_points(
        *[residue.xyz(atom) for atom in (center_atom, atom1, atom2, atom3)]
    )
    return stub_to_homog(rstub)


def rt_to_homog(rt):
    """
    Takes a rosetta RT and returns a h xform
    """
    return rotation_translation_to_homog(
        rt.get_rotation(), rt.get_translation()
    )


def invert_homog(xform):
    """
    Takes a homogenous xform and returns its inverse

    Note: the inverse *but not the transpose* of a homogenous xform is the
    inverse of the underlying rotation and translation components

    It's totally fair game to just allow numpy to compute the inverse
    """
    inv = _numpy.linalg.inv(xform)
    return inv


def homog_relative_transform(xform1, xform2):
    """
    returns the relative transform between two rotation translation matrices
    """
    xform1_inv = invert_homog(xform1)
    return xform1_inv @ xform2


def homog_relative_transform_from_stubs(stub1, stub2):
    """
    Converts Rosetta stubs to a homog relative xform between the two
    """
    hstub1, hstub2 = stub_to_homog(stub1), stub_to_homog(stub2)

    return homog_relative_transform(hstub1, hstub2)


def homog_relative_transform_from_residues(res1, res2):
    """
    Wrapper for making CA to CA homogenous relative transform between residues
    """
    return homog_relative_transform_from_stubs(
        _stub_from_residue(res1), _stub_from_residue(res2)
    )


def homog_super_transform(xform1, xform2):
    """
    returns the super xform between two rotation translation matrices

    super is defined by the left multiplied matrix that moves an object from
    the local xform1->xform2
    """
    xform1_inv = invert_homog(xform1)
    return xform2 @ xform1_inv


def homog_super_transform_from_stubs(stub1, stub2):
    """
    Converts Rosetta stubs to a homog super xform between the two

    super is defined by the left multiplied matrix that moves an object from
    the local xform1->xform2
    """
    hstub1, hstub2 = stub_to_homog(stub1), stub_to_homog(stub2)

    return homog_super_transform(hstub1, hstub2)


def homog_super_transform_from_residues(res1, res2):
    """
    Wrapper for making CA to CA homogenous super transform between residues

    super is defined by the left multiplied matrix that moves an object from
    the local xform1->xform2
    """
    return homog_super_transform_from_stubs(
        _stub_from_residue(res1), _stub_from_residue(res2)
    )


def super_from_rt(source, dest, rt):
    """
    returns a homog xform that superimposes source onto the rt from dest

    Takes homogenous xforms
    """
    return homog_super_transform(source, dest @ rt)


def superimpose_stub_by_rt(stub1, stub2, rt):
    """
    returns the homog xform representing the super from stub1 to rt of stub2

    wrapper for super_from_rt to take rosetta stubs and rts
    """
    return super_from_rt(
        stub_to_homog(stub1), stub_to_homog(stub2), rt_to_homog(rt)
    )


def stub_from_3_CA(pose, ca_num):
    """
    Returns a rosetta stub from CA i to i +2 (CA stub)

    Does not check that the res has a CA
    """
    assert len(pose.residues) < ca_num + 3, "Pose ends before ca_num + 2"
    targ_stub = Stub()
    targ_stub.from_four_points(
        *[
            pose.residue(num).xyz("CA")
            for num in [ca_num + 1, *range(ca_num, ca_num + 3)]
        ]
    )
    return targ_stub


def homog_from_3_CA(pose, ca_num):
    """
    Returns a homogeneous xform from CA i to i +2 (CA stub)

    Does not check that the res has a CA
    """
    return stub_to_homog(stub_from_3_CA(pose, ca_num))


def xform_magnitude(xform, radius_of_gyration=None):
    """
    Takes the magnitude of rotation and translation of the xform

    radius of gyration is just set to the norm of the translation vector if
    left unset

    Stolen from bcov's xform mag in rosetta (looks real different in python)
    """
    # Find squared norm of translation
    translation = xform[..., 3:][:-1]
    tnorm = _numpy.linalg.norm(translation)
    if not radius_of_gyration:
        radius_of_gyration = tnorm

    # Use clever trace hack to get cos( rotation_matrix().angle() )
    cos_theta = (xform[0][0] + xform[1][1] + xform[2][2] - 1.0) / 2.0

    # Calculate the sin() and then multiply by radius of gyration to get the rotation distance
    # Just use rg if we go past 90 degrees
    err_rot = (
        radius_of_gyration
        if cos_theta < 0
        else _math.sqrt(max(0.0, 1.0 - cos_theta ** 2)) * radius_of_gyration
    )

    # Combine the distances
    err = _math.sqrt(tnorm ** 2 + err_rot ** 2)

    return err
