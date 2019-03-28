import homog
import unittest
import numpy as _np


class TestCreation(unittest.TestCase):
    """
    Test the homog xform creator functions

    Does not test simple wrappers
    """

    def test_rotation_translation_to_homog(self):
        rotation = _np.array([[0, 1, 1], [1, 1, 2], [0, 2, 1]])
        translation = _np.array([[3], [2], [1]])
        hardcode = _np.array(
            [[0, 1, 1, 3], [1, 1, 2, 2], [0, 2, 1, 1], [0, 0, 0, 1]]
        )
        test_homog = homog.rotation_translation_to_homog(
            rotation=rotation, translation=translation
        )
        self.assertEqual(harcode, test_homog)

    def test_():
        ""
