import numpy as _numpy
from dzutils.pyrosetta_utils.geometry.rt_utils import (
    stub_from_residue as _stub_from_residue,
)


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


def stub_to_homog(stub):
    """
    Takes a rosetta stub and returns a h xform
    """
    return rotation_translation_to_homog(stub.M, stub.v)


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


def super_from_rt(source, dest, rt):
    """
    returns a homog xform that superimposes source onto the rt from dest

    Takes homogenous xforms
    """
    return dest @ rt @ invert_homog(source)


def superimpose_stub_by_rt(stub1, stub2, rt):
    """
    returns the homog xform representing the super from stub1 to rt of stub2

    wrapper for super_from_rt to take rosetta stubs and rts
    """
    return super_from_rt(
        stub_to_homog(stub1), stub_to_homog(stub2), rt_to_homog(rt)
    )
