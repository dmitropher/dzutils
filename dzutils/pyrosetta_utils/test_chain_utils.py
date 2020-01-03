import unittest
import pyrosetta
from unittest.mock import Mock as Mock
import dzutils.pyrosetta_utils.chain_utils as chain_utils
import dzutils.pyrosetta_utils.test_tools.test_poses as _tp


class TestChainUtils(unittest.TestCase):
    """
    Test for chain_utils
    """

    def setUp(self):
        """
        """
        pyrosetta.init(
            """
            -out:level 100
            """
        )

    # def test_posenum_in_chain(self):
    #     pose = Mock()
    #     pose.chain_begin = Mock(return_value=10)
    #     self.assertEqual(1, chain_utils.posnum_in_chain(pose, 10))

    def test_add_cut_out_of_bounds(self):
        pose = _tp.helical_monomer()
        self.assertRaises(RuntimeError, chain_utils.add_cut, pose, 0)
        self.assertRaises(
            RuntimeError, chain_utils.add_cut, pose, len(pose.residues) + 1
        )

    def test_add_cut(self):
        pose = _tp.helical_monomer()
        chain_utils.add_cut(pose, 2)
        self.assertTrue(pose.num_chains() == 2)
        self.assertTrue(pose.chain(2) == 1)
        self.assertTrue(pose.chain(3) == 2)

    def test_add_cut_copy(self):
        pose = _tp.helical_monomer()
        cut = chain_utils.add_cut(pose, 2)
        self.assertEqual(cut, pose)
        pose = _tp.helical_monomer()
        cut = chain_utils.add_cut(pose, 2, True)
        self.assertNotEqual(cut, pose)

    # def test_link_poses (self):
    #     chain_utils.link_poses()
