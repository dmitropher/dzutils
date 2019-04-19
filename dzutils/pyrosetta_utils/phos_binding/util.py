import _pyrosetta.rosetta as _pyr


def hbond_to_residue(pose, residue, vec=False):
    """
    Returns a list of all the hbonding residues to the given residue

    vec=True returns a rosetta vector instead (marginally faster)
    """
    if vec:
        return _pyr.core.scoring.hbonds.HBondSet(
            pose, bb_only=False
        ).residue_hbonds(residue)
    else:
        return [
            b
            for b in _pyr.core.scoring.hbonds.HBondSet(
                pose, bb_only=False
            ).residue_hbonds(residue)
        ]
