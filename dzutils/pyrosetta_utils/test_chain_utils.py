import unittest
from unittest.mock import Mock as Mock
import dzutils.pyrosetta_utils.chain_utils as chain_utils


class TestChainUtils(unittest.TestCase):
    """
    Test for chain_utils
    """

    def test_posenum_in_chain():
        pose = Mock()
        pose.chain_begin = Mock(return_value=10)
        assertEqual(1, chain_utils.posenum_in_chain(pose, 10))
