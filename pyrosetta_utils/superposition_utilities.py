import pyrosetta as _pyrosetta

# Superposition Transform
def superposition_pose(mob_pose, init_coords, ref_coords):
    """
    Thin wrapper to fuse the creation of the rotation/translation and the apply
    """

    rotation = _pyrosetta.rosetta.numeric.xyzMatrix_double_t()
    to_init_center = _pyrosetta.rosetta.numeric.xyzVector_double_t()
    to_fit_center = _pyrosetta.rosetta.numeric.xyzVector_double_t()
    _pyrosetta.rosetta.protocols.toolbox.superposition_transform(
        init_coords, ref_coords, rotation, to_init_center, to_fit_center
    )
    _pyrosetta.rosetta.protocols.toolbox.apply_superposition_transform(
        mob_pose, rotation, to_init_center, to_fit_center
    )
    return mob_pose


def super_resi_by_atom_name(mob_pose, targ_pose, mob_index, targ_index, *args):
    """
    superimposes mob_pose onto targ_pose: aligns by atom name from indexes given

    May provide any number of atom names as long as both are in each residue
    Offset can be provided to superimpose atoms of the same name with an offset
    (Useful for a lazy C symmetric rotation)
    """
    init_coords = (
        _pyrosetta.rosetta.utility.vector1_numeric_xyzVector_double_t()
    )
    ref_coords = (
        _pyrosetta.rosetta.utility.vector1_numeric_xyzVector_double_t()
    )
    mob_res = mob_pose.residue(mob_index)
    targ_res = targ_pose.residue(targ_index)

    index_range = range(len(args))

    for index in index_range:
        init_coords.append(mob_res.xyz(args[index]))
        ref_coords.append(targ_res.xyz(args[index]))
    superposition_pose(mob_pose, init_coords, ref_coords)
    return mob_pose
