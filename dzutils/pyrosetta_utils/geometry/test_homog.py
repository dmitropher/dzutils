import homog
import unittest
import numpy as _np


class TestCreation(unittest.TestCase):
    """
    Test the homog xform creator functions
    """

    def test_rotation_translation_to_homog(self):
        rosetta_matrix = Mock()
        rosetta_vector = Mock()
        rosetta_matrix.xx, rosetta_matrix.xy, rosetta_matrix.xz = 0, 1, 1
        rosetta_matrix.yx, rosetta_matrix.yy, rosetta_matrix.yz = 1, 1, 2
        rosetta_matrix.zx, rosetta_matrix.zy, rosetta_matrix.zz = 0, 2, 1
        rosetta_vector.x, rosetta_vector.y, rosetta_vector.z = 3, 2, 1

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
        rosetta_matrix.xx, rosetta_matrix.xy, rosetta_matrix.xz = 0, 1, 1
        rosetta_matrix.yx, rosetta_matrix.yy, rosetta_matrix.yz = 1, 1, 2
        rosetta_matrix.zx, rosetta_matrix.zy, rosetta_matrix.zz = 0, 2, 1
        rosetta_vector.x, rosetta_vector.y, rosetta_vector.z = 3, 2, 1
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
        rosetta_matrix.xx, rosetta_matrix.xy, rosetta_matrix.xz = 0, 1, 1
        rosetta_matrix.yx, rosetta_matrix.yy, rosetta_matrix.yz = 1, 1, 2
        rosetta_matrix.zx, rosetta_matrix.zy, rosetta_matrix.zz = 0, 2, 1
        rosetta_vector.x, rosetta_vector.y, rosetta_vector.z = 3, 2, 1
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
        assertEqual(
            inv @ hardcode2,
            homog.homog_relative_transform(hardcode1, hardcode2),
        )
