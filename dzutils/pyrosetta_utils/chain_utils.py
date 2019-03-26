import pyrosetta as _pyrosetta


def insert_pose(target_pose, in_pose, start, end):
    """
    Returns a pose with the in_pose inserted from start to end

    This function only adds residues by bond at the given sites, it does not
    search, align, or superimpose.
    """


def add_cut(pose, index):
    """
    Wrapper for the AddChainBreak() mover, does not add termini or rechains

    Useful for grafting something, or editing fold tree and conformation
    """
    # if already a terminus, check if there's an i+1 and if it is lower termini
    num_res = len(pose.residues)
    if index > num_res or index <= 0:
        raise RuntimeError("variable seqpos is out of range!")
    chain_break = (
        _pyrosetta.rosetta.protocols.protein_interface_design.movers.AddChainBreak()
    )
    chain_break.resnum(str(index))
    chain_break.change_foldtree(True)
    chain_break.change_conformation(True)
    outpose = _pyrosetta.rosetta.core.pose.Pose()
    outpose.detached_copy(pose)
    chain_break.apply(outpose)

    return outpose


def chain_break(pose, index):
    """
    Takes a pose and an index and breaks it into two chains at the index given

    invokes the addchainbreak mover and then adds upper and lower termini and
    rechains based on those termini

    Always rechains based on termini, always turns the cut into termini
    """
    num_res = len(pose.residues)

    outpose = add_cut(pose, index)

    _pyrosetta.rosetta.core.pose.add_upper_terminus_type_to_pose_residue(
        outpose, index
    )

    # check if there's a residue after the current index
    if index + 1 < num_res:
        _pyrosetta.rosetta.core.pose.add_lower_terminus_type_to_pose_residue(
            outpose, index + 1
        )
    outpose.conformation().chains_from_termini()
    return outpose
