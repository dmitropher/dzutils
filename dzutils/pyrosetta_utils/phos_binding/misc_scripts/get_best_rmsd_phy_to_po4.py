import sys
from itertools import permutations

import numpy as np

import pyrosetta

from dzutils.pyrosetta_utils import (
    run_pyrosetta_with_flags,
    residues_by_name,
    atom_indices_with_element,
    bonded_atoms,
)
from dzutils.pyrosetta_utils.geometry import rmsd


def best_rmsd_phy_to_po4(pose):
    """
    """
    po4_res, po4_resnum = residues_by_name(pose, "PO4", include_res=True)[0]
    phy_res, phy_resnum = residues_by_name(pose, "PHY", include_res=True)[0]

    p_atom_po4 = atom_indices_with_element(po4_res, "P")[0]
    p_atom_phy = atom_indices_with_element(phy_res, "P")[0]

    bonded_po4 = [po4_res.xyz(i) for i in bonded_atoms(po4_res, p_atom_po4)]
    bonded_phy = np.array(
        [phy_res.xyz(i) for i in bonded_atoms(phy_res, p_atom_phy)]
    )
    po4_permutations = np.array(
        list(permutations(bonded_po4, len(bonded_po4)))
    )
    phy_repeat = np.tile(bonded_phy, (len(po4_permutations), 1, 1))

    return np.amin(rmsd(phy_repeat, po4_permutations))


def main():
    """
    """
    pdb = sys.argv[1]
    output_path = sys.argv[2]
    flags_file = (
        "/home/dzorine/phos_binding/run_files/p_ligand_loud_v2.flags"
    )  # sys.argv[3] if len(sys.argv) > 3 else ""
    run_pyrosetta_with_flags(flags_file)
    pose = pyrosetta.pose_from_file(pdb)

    best_rmsd = best_rmsd_phy_to_po4(pose)

    with open(output_path, "w") as f:
        f.write(",".join([pdb, str(best_rmsd)]))
        f.write("\n")


if __name__ == "__main__":
    main()
