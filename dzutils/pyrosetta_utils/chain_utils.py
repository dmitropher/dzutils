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


def link_poses(*poses, rechain=False):
    """
    returns a pose of any number of poses appended end to end

    Order of input poses must be correct N->C, rechain appends by jump and calls
    chains_from_termini on the final conformation and copies but does not mess
    with pdb_info beyond what append_pose_by_jump does.

    No alignment, fold tree not smoothed, can take a single pose, returns a copy
    """
    assert bool(len(poses)), "number of input poses must be greater than 0"
    target = _pyrosetta.rosetta.core.pose.Pose()
    target.detached_copy(poses[0])
    if rechain:
        for i, pose in enumerate(poses[1:], 1):
            target.append_pose_by_jump(pose, i)
        target.conformation().chains_from_termini()
    else:
        for pose in poses[1:]:
            _pyrosetta.rosetta.core.pose.append_pose_to_pose(
                target, pose, False
            )
    return target


def insert_pose(target_pose, in_pose, start, end=0, smooth=False):
    """
    Returns a pose with the in_pose inserted from start to end

    If start and end are on different chains, residues are removed up to the
    last residue of "start"'s chain number, and from the beginning of "end"'s
    chain up to "end". If end is "0" in_pose is appended

    leaves in the cutpoints and jumps at the beginning and end of the insertion
    unless smooth is set. If the insertion joins two chains, smooth
    converts them to a single peptide edge and keeps any jumps. It does not
    yet support cyclic peptides.

    This function only adds residues by bond at the given sites, it does not
    search, align, or superimpose.

    This function is not responsible for deleting any poses you feed it.
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

    # insertion into a single chain
    inserted = _pyrosetta.rosetta.core.pose.Pose()
    linked_pose = _pyrosetta.rosetta.core.pose.Pose()

    # Case for insertion into a single chain
    if start_chain == end_chain:
        if end and end < start:
            raise NotImplementedError(
                "end < start and on one chain. Cyclic peptides are not supported"
            )

        # split out the desired chain
        chains = pose.split_by_chain()
        target_chain = chains[start_chain]

        # if end is defined, cut out the region between start and end
        # and splice
        if end:
            cut = add_cut(target_chain, end, True)
            cut.delete_residue_range_slow(start, end)
            cut_halves = cut.split_by_chain()
            ncut, ccut = cut_halves[1], cut_halves[2]
            inserted = link_poses(ncut, in_pose, ccut)
        # otherwise, if start is not the terminus,
        # remove residues up to the terminus
        else:
            if start < len(target_chain.resdiues):
                cut = add_cut(target_chain, start + 1, True)
                cut.delete_residue_range_slow(
                    start, len(target_chain.resdiues)
                )

            # And if it is the terminus, remove the terminus VariantType
            elif target_chain.residue(start).has_variant_type(
                _pyrosetta.rosetta.core.chemical.VariantType.UPPER_TERMINUS_VARIANT
            ):
                _pyrosetta.rosetta.core.pose.remove_upper_terminus_type_from_pose_residue()
            inserted = link_poses(target_chain, in_pose)

        # combine it all into one chain
        new_chains = (
            chain if i != start_chain else inserted
            for i, chain in enumerate(chains, 1)
        )
        inserted = link_poses(*new_chains, rechain=True)
    # insertion connecting two chains
    """
    code here
    """

    # rejoin chains,

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
