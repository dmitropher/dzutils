import unittest
import numpy as np
import pyrosetta

import dzutils.pyrosetta_utils.phos_binding.misc_scripts.anchored_graft_scan as agc
from dzutils.pyrosetta_utils.test_tools.test_poses import (
    helical_monomer,
    helical_fragment,
)


class TestAnchoredGraftScan(unittest.TestCase):
    """
    Test the inv rot table maker

    Does not test any of the rosetta wrappers yet
    """

    def setUp(self):
        pyrosetta.init(
            """
            -out:level 100
            """
        )
        self.insert = helical_fragment()
        self.scaffold = helical_monomer()

    def test_graft(self):
        """
        """
        grafted = agc.graft(
            self.scaffold,
            self.insert,
            50,
            65,
            n_label="test_n",
            c_label="test_c",
        )
        n_lab = pyrosetta.rosetta.core.select.residue_selector.ResiduePDBInfoHasLabelSelector(
            "test_n"
        )
        c_lab = pyrosetta.rosetta.core.select.residue_selector.ResiduePDBInfoHasLabelSelector(
            "test_c"
        )
        unittest.assertTrue(n_lab.apply(grafted)[50])
        unittest.assertTrue(c_lab.apply(grafted)[65])
