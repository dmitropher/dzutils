import pyrosetta as _pyrosetta
import itertools as _itertools


def helixResList(pose):
    """
    Takes a pose, returns a list of residues that are in helical ABEGO bins
    """
    # create and initialize residue selector for alpha helical torsion angles
    _pyrosetta.distributed.maybe_init()
    binSel = _pyrosetta.rosetta.core.select.residue_selector.BinSelector()
    binSel.set_bin_name("A")
    binSel.initialize_and_check()

    # enumerate helical residues, use groupby to cluster runs of 1
    heliter = _itertools.groupby(
        [(i, is_helix) for i, is_helix in enumerate(binSel.apply(pose), 1)],
        lambda x: x[1],
    )
    runs = [
        tuple(item[0] for item in group[1]) for group in heliter if group[0]
    ]
    bounds = [(run[0], run[-1]) for run in runs]
    return bounds


def longestHelix(pose):
    """
    returns the longest helix in a pose
    """
    heList = helixResList(pose)
    if not heList:
        return None
    long = [0, 0]
    for h in heList:
        if (h[1] - h[0]) > (long[1] - long[0]):
            long = h
    if long == [0, 0]:
        return None
    else:
        return long
