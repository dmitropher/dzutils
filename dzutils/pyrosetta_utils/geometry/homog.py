import numpy as _numpy
from dzutils.pyrosetta_utils.geometry.rt_utils import (
    stub_from_residue as _stub_from_residue,
)


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
    Takes a rosetta stube and returns a h xform
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


def homog_relative_transform_from_stubs(stub1, stub2):
    """
    Takes two Rosetta stubs and converts them to a homogenous relative transform between the two
    """
    hstub1, hstub2 = stub_to_homog(stub1), stub_to_homog(stub2)
    hstub1_inv = invert_homog(hstub1)
    return hstub1_inv @ hstub2


def homog_relative_transform_from_residues(res1, res2):
    """
    Wrapper for making CA to CA homogenouse relative transform between residues
    """
    return homog_relative_transform_from_stubs(
        _stub_from_residue(res1), _stub_from_residue(res2)
    )
