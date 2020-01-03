import pyrosetta as _pyrosetta

FILEPATH = __file__.split("test_poses.py")[0]


def rossmann_monomer():
    """
    returns a pose of pdb structure 2kpo

    Loads from a pdb file, this is deliberately designed to be brittle to detect
    when pyrosetta behavior has changed.
    """
    return _pyrosetta.pose_from_file(f"{FILEPATH}/test_pdb/rossmann_2kpo.pdb")


def helical_monomer():
    """
    returns a pose of pdb structure 4uos

    Loads from a pdb file, this is deliberately designed to be brittle to detect
    when pyrosetta behavior has changed.
    """
    return _pyrosetta.pose_from_file(f"{FILEPATH}/test_pdb/helical_4uos.pdb")


def beta_fragment():
    """
    returns a pose of pdb structure

    Loads from a pdb file, this is deliberately designed to be brittle to detect
    when pyrosetta behavior has changed.
    """
    return _pyrosetta.pose_from_file(f"{FILEPATH}/test_pdb/beta_fragment.pdb")


def helical_fragment():
    """
    returns a helical fragment of 4uos from 51 to 64
    """
    return _pyrosetta.pose_from_file(f"{FILEPATH}/test_pdb/helifrag.pdb")
