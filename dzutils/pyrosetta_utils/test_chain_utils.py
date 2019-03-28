import unittest
from unittest.mock import Mock as Mock
import dzutils.pyrosetta_utils.chain_utils as chain_utils
import dzutils.pyrosetta_utils.test_tools.test_poses as _tp


class TestChainUtils(unittest.TestCase):
    """
    Test for chain_utils
    """

    def test_posenum_in_chain(self):
        pose = Mock()
        pose.chain_begin = Mock(return_value=10)
        assertEqual(1, chain_utils.posenum_in_chain(pose, 10))

    def test_add_cut_out_of_bounds(self):
        pose = _tp.helical_monomer()
        self.assertRaises(RuntimeError, chain_utils.add_cut(pose, 0))
        self.assertRaises(
            RuntimeError, chain_utils.add_cut(pose, pose.residues + 1)
        )

    def test_add_cut(self):
        pose = _tp.helical_monomer()
        cut = chain_utils.add_cut(pose, 2)
        self.assertTrue(cut.chains() == 2)
        self.assertTrue(cut.chain(1) == 1)
        self.assertTrue(cut.chain(2) == 2)

    def test_add_cut_copy(self):
        pose = _tp.helical_monomer()
        cut = chain_utils.add_cut(pose, 2)
        self.assertEqual(cut, pose)
        pose = _tp.helical_monomer()
        cut = chain_utils.add_cut(pose, 2, True)
        self.assertNotEqual(cut, pose)

    def test_link_poses(self):
        chain_utils.link_poses()
