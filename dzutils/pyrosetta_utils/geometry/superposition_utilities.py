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


def super_resi_by_bb(mob_pose, targ_pose, mob_index, targ_index):
    init_coords = (
        _pyrosetta.rosetta.utility.vector1_numeric_xyzVector_double_t()
    )
    ref_coords = (
        _pyrosetta.rosetta.utility.vector1_numeric_xyzVector_double_t()
    )
    mob_res = mob_pose.residue(mob_index)
    targ_res = targ_pose.residue(targ_index)
    for atom in ("N", "CA", "C"):
        init_coords.append(mob_res.xyz(atom))
        ref_coords.append(targ_res.xyz(atom))
    superposition_pose(mob_pose, init_coords, ref_coords)
    return mob_pose


def super_by_residues(mob_pose, targ_pose, mob_index, targ_index, *args):
    """
    superimposes mob_pose onto targ_pose: aligns by atom name from indexes given

    May provide any number of atom names as long as both are in each residue
    Offset can be provided to superimpose atoms of the same name with an offset
    (Useful for a lazy C symmetric rotation)

    Defaults to superimposing by bb atoms if no atoms given
    """
    init_coords = (
        _pyrosetta.rosetta.utility.vector1_numeric_xyzVector_double_t()
    )
    ref_coords = (
        _pyrosetta.rosetta.utility.vector1_numeric_xyzVector_double_t()
    )
    mob_res = mob_pose.residue(mob_index)
    targ_res = targ_pose.residue(targ_index)
    num_atoms = len(args)
    if num_atoms == 0:
        return super_resi_by_bb(mob_pose, targ_pose, mob_index, targ_index)
    index_range = range(num_atoms)

    for index in index_range:
        init_coords.append(mob_res.xyz(args[index]))
        ref_coords.append(targ_res.xyz(args[index]))
    superposition_pose(mob_pose, init_coords, ref_coords)
    return mob_pose


def super_by_paired_atoms(
    mob_pose, targ_pose, mob_index, targ_index, *atom_pairs
):
    """
    super mob onto targ, uses tuples in atom_pairs (mob_atom1,targ_atom1),etc

    May provide any number of atom pairs as long as both are in each residue
    """
    init_coords = (
        _pyrosetta.rosetta.utility.vector1_numeric_xyzVector_double_t()
    )
    ref_coords = (
        _pyrosetta.rosetta.utility.vector1_numeric_xyzVector_double_t()
    )
    mob_res = mob_pose.residue(mob_index)
    targ_res = targ_pose.residue(targ_index)
    num_atoms = len(atom_pairs)
    assert bool(num_atoms > 0), "no atoms pairs given"
    for mob_atom, targ_atom in atom_pairs:
        init_coords.append(mob_res.xyz(mob_atom))
        ref_coords.append(targ_res.xyz(targ_atom))
    superposition_pose(mob_pose, init_coords, ref_coords)
    return mob_pose


def align_to_tetrahedral_rotation(mob_pose, mob_resnum, rotation=0, *args):
    """
    Align the pose to a rotation of the atoms specified

    only for tetrahedral symmetric rotations
    helper function for enumerating all tetrahedral symmetries

    rotations are numbered from 0-11, with wraparound

    Uses the rosetta toolbox superpostion transform, not the inverse RT method
    """
    if len(args) > 4:
        raise AssertionError("This function is for tetrahedral symmetry only")

    init_coords = (
        _pyrosetta.rosetta.utility.vector1_numeric_xyzVector_double_t()
    )

    ref_coords = (
        _pyrosetta.rosetta.utility.vector1_numeric_xyzVector_double_t()
    )

    nub_res = mob_pose.residue(mob_resnum)

    # modulo arithmetic to change what's aligned to what
    fixed_point = rotation % 4
    face_rotation = rotation // 4 % 3
    rotated_args = [
        args[
            ((i + fixed_point) % 4) * ((fixed_point + 1) % 2)
            + ((fixed_point - i) % 4) * ((fixed_point) % 2)
        ]
        for i in range(4)
    ]
    for i in range(4):
        targ = (
            rotated_args[0]
            if i == 0
            else rotated_args[1:][(i + face_rotation - 1) % 3]
        )
        init = args[i]
        init_coords.append(nub_res.xyz(targ))
        ref_coords.append(nub_res.xyz(init))
    superposition_pose(mob_pose, init_coords, ref_coords)(
        mob_pose, init_coords, ref_coords
    )
    return mob_pose
