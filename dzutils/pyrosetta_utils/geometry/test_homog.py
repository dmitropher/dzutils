import homog
import unittest
from unittest.mock import Mock as Mock
import numpy as _np
import dzutils.pyrosetta_utils.test_tools.test_poses as _tp


class TestHomogUtils(unittest.TestCase):
    """
    Test the homog xform creator functions
    """

    def test_rotation_translation_to_homog(self):
        rosetta_matrix = Mock()
        rosetta_vector = Mock()
        rosetta_matrix1.xx, rosetta_matrix1.xy, rosetta_matrix1.xz = 0, 1, 1
        rosetta_matrix1.yx, rosetta_matrix1.yy, rosetta_matrix1.yz = 1, 1, 2
        rosetta_matrix1.zx, rosetta_matrix1.zy, rosetta_matrix1.zz = 0, 2, 1
        rosetta_vector1.x, rosetta_vector1.y, rosetta_vector1.z = 3, 2, 1

        hardcode = _np.array(
            [[0, 1, 1, 3], [1, 1, 2, 2], [0, 2, 1, 1], [0, 0, 0, 1]]
        )
        test_homog = homog.rotation_translation_to_homog(
            rotation=rosetta_matrix, translation=rosetta_vector
        )
        self.assertEqual(harcode, test_homog)

    def test_stub_to_homog(self):
        rosetta_matrix = Mock()
        rosetta_vector = Mock()
        rosetta_matrix1.xx, rosetta_matrix1.xy, rosetta_matrix1.xz = 0, 1, 1
        rosetta_matrix1.yx, rosetta_matrix1.yy, rosetta_matrix1.yz = 1, 1, 2
        rosetta_matrix1.zx, rosetta_matrix1.zy, rosetta_matrix1.zz = 0, 2, 1
        rosetta_vector1.x, rosetta_vector1.y, rosetta_vector1.z = 3, 2, 1
        stub = Mock()
        stub.M = rosetta_matrix
        stub.v = rosetta_vector
        hardcode = _np.array(
            [[0, 1, 1, 3], [1, 1, 2, 2], [0, 2, 1, 1], [0, 0, 0, 1]]
        )
        test_homog = homog.stub_to_homog(stub)
        self.assertEqual(hardcode, test_homog)

    def test_rt_to_homog(self):
        rosetta_matrix = Mock()
        rosetta_vector = Mock()
        rosetta_matrix1.xx, rosetta_matrix1.xy, rosetta_matrix1.xz = 0, 1, 1
        rosetta_matrix1.yx, rosetta_matrix1.yy, rosetta_matrix1.yz = 1, 1, 2
        rosetta_matrix1.zx, rosetta_matrix1.zy, rosetta_matrix1.zz = 0, 2, 1
        rosetta_vector1.x, rosetta_vector1.y, rosetta_vector1.z = 3, 2, 1
        rt = Mock()
        rt.get_rotation = rosetta_matrix
        rt.get_translation = rosetta_vector
        hardcode = _np.array(
            [[0, 1, 1, 3], [1, 1, 2, 2], [0, 2, 1, 1], [0, 0, 0, 1]]
        )
        test_homog = homog.rt_to_homog(rt)
        self.assertEqual(hardcode, test_homog)

    def test_invert_homog(self):
        hardcode = _np.array(
            [[0, 1, 1, 3], [1, 1, 2, 2], [0, 2, 1, 1], [0, 0, 0, 1]]
        )

        inv = _np.linalg.inv(hardcode)
        homog.invert_homog(hardcode)
        self.assertEqual(hardcode, hardcode)

    def test_homog_relative_transform(self):
        hardcode1 = _np.array(
            [[0, 1, 1, 3], [1, 1, 2, 2], [0, 2, 1, 1], [0, 0, 0, 1]]
        )
        hardcode2 = _np.array(
            [[0, 2, 1, 2], [1, 1, 2, 2], [0, 1, 1, 4], [0, 0, 0, 1]]
        )
        inv = _np.linalg.inv(hardcode1)
        self.assertEqual(
            inv @ hardcode2,
            homog.homog_relative_transform(hardcode1, hardcode2),
        )

    def test_homog_relative_transform_from_stubs(self):

        stub1, stub2 = Mock(), Mock()
        rosetta_matrix1, rosetta_matrix2 = Mock(), Mock()
        rosetta_vector1, rosetta_vector2 = Mock(), Mock()
        rosetta_matrix1.xx, rosetta_matrix1.xy, rosetta_matrix1.xz = 0, 1, 1
        rosetta_matrix1.yx, rosetta_matrix1.yy, rosetta_matrix1.yz = 1, 1, 2
        rosetta_matrix1.zx, rosetta_matrix1.zy, rosetta_matrix1.zz = 0, 2, 1
        rosetta_vector1.x, rosetta_vector1.y, rosetta_vector1.z = 3, 2, 1
        rosetta_matrix2.xx, rosetta_matrix2.xy, rosetta_matrix2.xz = 0, 1, 1
        rosetta_matrix2.yx, rosetta_matrix2.yy, rosetta_matrix2.yz = 1, 1, 2
        rosetta_matrix2.zx, rosetta_matrix2.zy, rosetta_matrix2.zz = 0, 2, 1
        rosetta_vector2.x, rosetta_vector2.y, rosetta_vector2.z = 3, 2, 1
        stub1.M = rosetta_matrix1
        stub1.v = rosetta_vector1
        stub2.M = rosetta_matrix2
        stub2.v = rosetta_vector2
        ident = _np.identity(4)
        self.assertEqual(
            ident, homog.homog_relative_transform_from_stubs(stub1, stub2)
        )

    def test_homog_relative_transform_from_residues(self):
        import rt_utils as _rtu

        pose = _tp.helical_monomer
        stub1 = _rtu.stub_from_residue(pose.residue(1))
        stub2 = _rtu.stub_from_residue(pose.residue(1))
        ident = _np.identity(4)
        self.assertEqual(
            ident, homog.homog_relative_transform_from_residues(stub1, stub2)
        )
