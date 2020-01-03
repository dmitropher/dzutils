import sys
import unittest
from unittest.mock import Mock as Mock

import logging

logger = logging.getLogger("anchored_graft")
logger.level = logging.DEBUG


import pyrosetta

import dzutils.pyrosetta_utils.anchored_graft_scan as agc
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
        make_poly_A_mover = pyrosetta.rosetta.protocols.pose_creation.MakePolyXMover(
            "ALA", False, False, False
        )
        p_a_insert = self.insert.clone()
        make_poly_A_mover.apply(p_a_insert)
        self.poly_a_insert = p_a_insert
        stream_handler = logging.StreamHandler(sys.stdout)
        logger.addHandler(stream_handler)

    def test_graft(self):
        """
        """

        insert = self.poly_a_insert.clone()

        grafted = agc.graft(
            self.scaffold.clone(),
            insert,
            51,
            66,
            n_label="test_n",
            c_label="test_c",
        )
        n_lab = pyrosetta.rosetta.core.select.residue_selector.ResiduePDBInfoHasLabelSelector(
            "test_n"
        )
        c_lab = pyrosetta.rosetta.core.select.residue_selector.ResiduePDBInfoHasLabelSelector(
            "test_c"
        )

        self.assertTrue(n_lab.apply(grafted)[51])
        self.assertTrue(c_lab.apply(grafted)[66])

        a_selector = pyrosetta.rosetta.core.select.residue_selector.ResidueNameSelector(
            "ALA"
        )
        a_res = a_selector.apply(grafted)
        targs = list(a_res)[51:65]

        targ_list = list(set(targs))

        self.assertTrue(len(targ_list) == 1 and targ_list[0])

    def test_graft_fragment(self):
        """
        this test basically just sees if the thing runs, doesnt really test
        """
        site = Mock()
        site.start_pos = 52
        site.end_pos = 65
        insert = self.poly_a_insert.clone()
        agc.graft_fragment(self.scaffold.clone(), insert, site, use_start=True)
        agc.graft_fragment(
            self.scaffold.clone(), insert, site, use_start=False
        )
