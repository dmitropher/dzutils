import pyrosetta as _pyrosetta


def chain_break(pose, index):
    """
    Takes a pose and an index and breaks it into two chains at the index given

    invokes the addchainbreak mover and then adds upper and lower termini and
    rechains based on those termini

    Always rechains based on termini, always turns the cut into termini
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
    chain_break.apply(pose)

    _pyrosetta.rosetta.core.pose.add_lower_terminus_type_to_pose_residue(
        pose, index
    )

    # check if there's a residue after the current index
    if index - 1 < num_res:
        _pyrosetta.rosetta.core.pose.add_upper_terminus_type_to_pose_residue(
            pose, index - 1
        )
    pose.conformation().chains_from_termini()
    return pose
