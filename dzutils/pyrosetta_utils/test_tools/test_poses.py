import pyrosetta as _pyrosetta


def rossmann_monomer():
    """
    returns a pose of pdb structure 2kpo

    Loads from a pdb file, this is deliberately designed to be brittle to detect
    when pyrosetta behavior has changed.
    """
    _pyrosetta.pose_from_file("./test_pdb/rossmann_2kpo.pdb")


def helical_monomer():
    """
    returns a pose of pdb structure 4uos

    Loads from a pdb file, this is deliberately designed to be brittle to detect
    when pyrosetta behavior has changed.
    """
    _pyrosetta.pose_from_file("./test_pdb/helical_4uos.pdb")


def beta_fragment():
    """
    returns a pose of pdb structure

    Loads from a pdb file, this is deliberately designed to be brittle to detect
    when pyrosetta behavior has changed.
    """
    _pyrosetta.pose_from_file("./test_pdb/beta_fragment.pdb")
