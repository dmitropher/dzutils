import unittest
import numpy as np
import pyrosetta

import dzutils.pyrosetta_utils.phos_binding.misc_scripts.make_inverse_rot_tables as make_table
from dzutils.pyrosetta_utils import residue_type_from_name3
from dzutils.pyrosetta_utils.geometry.pose_xforms import RotamerRTArray


class TestMakeInverseRotTables(unittest.TestCase):
    """
    Test the inv rot table maker

    Does not test any of the rosetta wrappers yet
    """

    def setUp(self):
        pyrosetta.init(
            """
            -out:level 100
            -extra_res_fa /home/dzorine/projects/phos_binding/params_files/p_loop_params/PHY_4_chi.params
            """
        )
        restype = residue_type_from_name3("PHY")
        self.residue = pyrosetta.rosetta.core.conformation.Residue(
            restype, True
        )
        self.base_atoms = ("CA", "N", "CA", "C")
        self.target_atoms = ("P1", "P1", "O0", "O3P")
        self.rotamer_rt_array = RotamerRTArray(
            residue=self.residue,
            base_atoms=self.base_atoms,
            target_atoms=self.target_atoms,
            inverse=True,
        )

    def test_rt_from_chis(self):
        test_rt = make_table.rt_from_chis(
            self.rotamer_rt_array,
            10,
            10,
            10,
            10,
            base_atoms=self.base_atoms,
            target_atoms=self.target_atoms,
        )
        self.rotamer_rt_array.reset_rotamer(10, 10, 10, 10)
        actual = np.array(self.rotamer_rt_array)
        self.assertTrue(np.array_equiv(test_rt, actual))
