import pyrosetta as _pyrosetta
import itertools as _itertools

# components of this module will likely fail to work if you have not run init()
_pyrosetta.distributed.maybe_init()
# takes a pose and returns a list containing lists of consecutive residues with
# helical torsion angles
def helixResList(pose):
    """
    Takes a pose, returns a list of residues that are in helical ABEGO bins
    """
    # create and initialize residue selector for alpha helical torsion angles
    bsc = _pyrosetta.rosetta.core.select.residue_selector.BinSelectorCreator()
    binSel = bsc.create_residue_selector()
    binSel.set_bin_name("A")
    binSel.initialize_and_check()

    # apply the residue selector to the pose and save the vector
    heliVec = binSel.apply(pose)

    # use the residue selector to build a list of residue ranges where
    # more than one consecutive residue has helical torsions.
    heliNums = zip(range(heliVec.capacity()), heliVec)
    heliter = _itertools.groupby(heliNums, lambda x: x[1])
    heliOut = []
    for k, v in heliter:
        a = []
        if k:
            for res, b in v:
                a.append(res)
        if len(a) > 1:
            heliOut.append([a[0], a[-1]])

    return heliOut


# returns the pos numbered first and last residue number of the longest
# helix of the pose in a two member list
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
