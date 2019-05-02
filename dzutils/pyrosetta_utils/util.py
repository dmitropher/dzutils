import pyrosetta.rosetta as _pyr


def residues_by_name(pose, *names, include_res=False):
    """
    Returns a list of resnums where the name3 is the in names

    include_res returns a list of tuples where the first element is the residue
    object and the second is the resnum
    Mostly exists because the residue name selector doesn't always work
    """
    if include_res:
        return [
            (pose.residue(i), i)
            for i, is_res in enumerate(pose.residues, 1)
            if is_res.name3() in names
        ]
    else:
        return [
            i
            for i, is_res in enumerate(pose.residues, 1)
            if is_res.name3() in names
        ]


def atom_distance(a1, a2):
    """
    Takes rosetta atoms, returns the norm of a1 - a2
    """
    try:
        return (a1.xyz() - a2.xyz()).norm()
    except:
        return None


def index_neighbors(pose, index, dist):
    """
    Returns a list of indexes that are neighbors to the given index
    """
    return [
        resi
        for resi, is_neib in enumerate(
            _pyr.core.select.residue_selector.NeighborhoodResidueSelector(
                _pyr.core.select.residue_selector.ResidueIndexSelector(
                    str(index)
                ),
                dist,
            ).apply(pose),
            1,
        )
        if is_neib
    ]


def get_index_atom(pose, index, atom):
    try:
        return pose.residue(index).atom(atom)
    except:
        return None
