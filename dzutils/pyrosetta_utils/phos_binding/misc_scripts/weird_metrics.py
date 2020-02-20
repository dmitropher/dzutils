from pyrosetta.rosetta.core.scoring.dssp import Dssp


def fraction_ss(pose, *dssp_types):
    """
    Returns the percent of the pose with the given dssp types
    """
    l = len(pose.residues)
    ds = Dssp(pose).get_dssp_secstruct()
    counts = {type: 0 for type in dssp_types}
    for c in ds:
        if c in dssp_types:
            counts[c] += 1
    return [counts[c] / l for c in dssp_types]
