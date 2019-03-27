import pyrosetta as _pyrosetta


def add_cut(pose, index, new_pose=False):
    """
    Wrapper for the AddChainBreak() mover, does not add termini or rechains

    Useful for grafting something, or editing fold tree and conformation
    if new_pose is true, perform the cut on a copy instead of on the input pose
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
    if new_pose:
        outpose = _pyrosetta.rosetta.core.pose.Pose()
        outpose.detached_copy(pose)
        chain_break.apply(outpose)

        return outpose
    else:
        chain_break.apply(pose)
        return pose


def insert_pose(target_pose, in_pose, start, end, smooth=False):
    """
    Returns a pose with the in_pose inserted from start to end

    If start and end are on different chains, residues are removed up to the
    last residue of "start"'s chain number, and from the beginning of "end"'s
    chain up to "end".

    leaves in the cutpoints and jumps at the beginning and end of the insertion
    unless smooth is set. If the insertion joins two chains, smooth
    converts them to a single peptide edge and keeps any jumps. It does not
    yet support cyclic peptides.

    This function only adds residues by bond at the given sites, it does not
    search, align, or superimpose.

    this function is not responsible for deleting any poses you feed it.
    """
    pose = target_pose.clone()
    pose_len = len(pose.residues)
    assert bool(
        start > 0 and start <= pose_len
    ), "Start residue for grafting must be between 1 and end of the pose"

    assert bool(
        start > 0 and start <= pose_len
    ), "End residue for grafting must be between 1 and end of the pose"
    start_chain = pose.chain(start)
    end_chain = pose.chain(end)
    if start_chain == end_chain and end < start:
        raise NotImplementedError(
            "end < start and on one chain. Cyclic peptides are not supported"
        )
    # Insertion into one chains
    """
    code here
    """
    # insertion connecting two chains
    """
    code here
    """
    # maybe smooth fold tree here


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
