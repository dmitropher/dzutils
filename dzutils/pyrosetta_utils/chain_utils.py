import pyrosetta as _pyrosetta
from .util import residues_by_name as residues_by_name


def chains_with_resname(pose, resname):
    """
    Takes a pose and a residue name, returns a list of chains with that residue
    """
    residue_name_selector = (
        _pyrosetta.rosetta.core.select.residue_selector.ResidueNameSelector()
    )
    residue_name_selector.set_residue_names(resname)
    vec = residue_name_selector.apply(pose)
    chains_redundant = [
        pose.chain(i) for i, is_selected in enumerate(vec, 1) if (is_selected)
    ]
    return list(set(chains_redundant))


def pose_from_chain(pose, chain):
    """
    Returns a pose consisting of just the input chain

    Basically just a wrapper to make code prettier, slightly slower than just
    invoking the bound method directly.
    """
    return pose.split_by_chain()[chain]


def pose_excluding_chain(pose, chain):
    """
    Returns a pose without the listed chain
    """
    chains = [p for i, p in enumerate(pose.split_by_chain(), 1) if i != chain]
    new_pose = chains[0]
    for i, chain in enumerate(chains[1:], 1):
        new_pose.append_pose_by_jump(chain, i)
    return new_pose


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


def trim_pose_to_term(pose, target, terminus=None):
    """
    Removes residues from start to chosen terminus, returns the pose

    terminus must be "chain_begin" or "chain_end"
    """
    assert bool(
        terminus == "chain_begin" or terminus == "chain_end"
    ), "terminus must be specified as 'chain_begin' or 'chain_end'"
    assert bool(pose.num_chains() == 1), "input pose must be a single chain"

    if terminus == "chain_begin":
        if target > 1:
            pose.delete_residue_range_slow(1, target)

        # And if it is the terminus, remove the terminus VariantType
        elif pose.residue(target).has_variant_type(
            _pyrosetta.rosetta.core.chemical.VariantType.LOWER_TERMINUS_VARIANT
        ):
            _pyrosetta.rosetta.core.conformation.remove_lower_terminus_type_from_conformation_residue(
                pose.conformation(), target
            )

    if terminus == "chain_end":
        if target < len(pose.residues):
            pose.delete_residue_range_slow(target, len(pose.residues))

        # And if it is the terminus, remove the terminus VariantType
        elif pose.residue(target).has_variant_type(
            _pyrosetta.rosetta.core.chemical.VariantType.UPPER_TERMINUS_VARIANT
        ):
            _pyrosetta.rosetta.core.conformation.remove_upper_terminus_type_from_conformation_residue(
                pose.conformation(), target
            )

    return pose


def chain_break(pose, index):
    """
    Takes a pose and an index and breaks it into two chains at the index given

    invokes the addchainbreak mover and then adds upper and lower termini and
    rechains based on those termini

    Always rechains based on termini, always turns the cut into termini
    """
    num_res = len(pose.residues)

    try:
        outpose = add_cut(pose, index)
    except RuntimeError as e:
        outpose = pose
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


def rechain_resname(pose, resname):
    """
    Returns a pose where all residues of the given name are in separate chains
    """
    residues = residues_by_name(pose, resname)
    target_poses = []
    other_pose = pose.clone()
    for res in residues:
        res_pose = _pyrosetta.rosetta.core.pose.Pose()
        res_pose.append_residue_by_bond(pose.residue(res))
        target_poses.append(res_pose)
        other_pose.delete_residue_range_slow(res, res)
    return link_poses(other_pose, *target_poses, rechain=True)


def serial_deletions(pose, target, terminus=None):
    """
    Returns list of all serial dels to the chain_begin/end, up to target residue
    """
    if terminus == "chain_begin":
        return [
            trim_pose_to_term(pose.clone(), i, terminus)
            for i in range(target, 1, -1)
        ]
    if terminus == "chain_end":
        chain_end = len(pose.residues)
        return [
            trim_pose_to_term(pose.clone(), i, terminus)
            for i in range(target, chain_end, 1)
        ]
