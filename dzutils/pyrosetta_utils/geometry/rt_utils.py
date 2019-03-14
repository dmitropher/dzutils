import pyrosetta as _pyrosetta


def stub_from_residue(
    residue, center_atom="CA", atom1="N", atom2="CA", atom3="C"
):
    """
    Returns a stub. A wrapper for atom.xyz with the default of the bb atoms.
    """
    return _pyrosetta.rosetta.core.kinematics.Stub(
        residue.atom(center_atom).xyz(),
        residue.atom(atom1).xyz(),
        residue.atom(atom2).xyz(),
        residue.atom(atom3).xyz(),
    )


def rt_from_res_atoms(
    start_res,
    end_res,
    start_center_atom="CA",
    end_center_atom="CA",
    start_atom1="N",
    start_atom2="CA",
    start_atom3="C",
    end_atom1="N",
    end_atom2="CA",
    end_atom3="C",
):
    """
    Returns the RT from the start to end res: centers stub around given atoms

    Defaults to CA as the center atom name, N-CA-C as the plane atoms
    """

    stub1 = stub_from_residue(
        start_res, start_center_atom, start_atom1, start_atom2, start_atom3
    )
    stub2 = stub_from_residue(
        end_res, end_center_atom, end_atom1, end_atom2, end_atom3
    )
    rt = _pyrosetta.rosetta.core.kinematics.RT(stub1, stub2)
    return rt


def peptide_bond_rt(start_res, end_res):
    """
    Returns the RT from the N term of start_res to the C of end_res
    """
    rt = rt_from_res_atoms(
        start_res,
        end_res,
        start_center_atom="O",
        start_atom1="O",
        end_center_atom="O",
        end_atom1="O",
    )
    return rt


def residue_ca_rt_dist_sq(source_res1, dest_res1, source_res2, dest_res2):
    """
    Returns the distance_squared of the RT of two residue pairs

    Uses the CA to CA RT of each pair
    """
    return rt_from_res_atoms(source_res1, dest_res1).distance_squared(
        rt_from_res_atoms(source_res2, dest_res2)
    )


def residue_pep_rt_dist_sq(source_res1, dest_res1, source_res2, dest_res2):
    """
    Returns the distance_squared of the RT of two residue pairs

    Uses the C,N,O to C,N,O RT of each pair
    """
    return peptide_bond_rt(source_res1, dest_res1).distance_squared(
        peptide_bond_rt(source_res2, dest_res2)
    )
