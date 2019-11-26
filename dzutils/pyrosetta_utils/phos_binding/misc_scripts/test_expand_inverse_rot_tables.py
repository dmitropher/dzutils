import unittest
import numpy as np
import pyrosetta

from itertools import product

import dzutils.pyrosetta_utils.phos_binding.misc_scripts.expand_inverse_rot_tables as exp_table
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

    def test_numba_cartesian(self):
        input_data = [[i + j for i in range(4)] for j in range(0, 4 * 4, 4)]
        it_out = np.array(product(*input_data))
        test_out = exp_table.numba_cartesian(np.array(input_data))
        self.assertTrue(
            np.array_equiv(np.sort(it_out, axis=0)), np.sort(test_out, axis=0)
        )

    def test_bit_pack_rotamers(self):
        """
        very basic test, does not test rounding for degrees or using more chis

        The test value chosen for easy validation via np.binary_repr
        """
        num = 3
        round = 1
        dummy_chis = np.array([[1] * num])
        packed = exp_table.bit_pack_rotamers(dummy_chis, num, round)
        self.assertEqual(np.uint64(4398048608257), packed)

    def test_unpack_chis(self):
        dummy_packed = np.array([np.uint64(4398048608257)] * 3)
        unpacked = exp_table.unpack_chis(dummy_packed, 3, 1)
        self.assertTrue(np.array_equiv(np.array([[1] * 3]), unpacked))
