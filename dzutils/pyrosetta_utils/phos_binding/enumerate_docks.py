<<<<<<< HEAD:dzutils/pyrosetta_utils/phos_binding/enumerate_docks.py
import pyrosetta as _pyrosetta
import itertools as _itertools
=======
import _pyrosetta
import itertools
import dzutils
import _pyrosetta.distributed
import _pyrosetta.distributed.tasks.rosetta_scripts
import _pyrosetta.distributed.io
>>>>>>> fix git ignore and imports:phos_binding/enumerate_docks.py


# Atom's random rama function bind ###
@_pyrosetta.bindings.utility.bind_method(_pyrosetta.rosetta.core.pose.Pose)
def randomize_rama(self, residue_index):
    upper_connection_residue_index = self.residues[
        residue_index
    ].connected_residue_at_upper()

    if not upper_connection_residue_index:
        return

    sm = _pyrosetta.rosetta.core.scoring.ScoringManager.get_instance()
    rama = sm.get_RamaPrePro()
    torsions = _pyrosetta.rosetta.utility.vector1_double()

    rama.random_mainchain_torsions(
        self.conformation(),
        self.residues[residue_index].type(),
        self.residues[upper_connection_residue_index].type(),
        torsions,
    )

    for i in range(1, len(self.residues[residue_index].mainchain_torsions())):
        self.set_torsion(
            _pyrosetta.rosetta.core.id.TorsionID(
                residue_index, _pyrosetta.rosetta.core.id.BB, i
            ),
            torsions[i],
        )


# Superposition Transform
def compute_and_apply_superposition_transform(
    mob_pose, init_coords, ref_coords
):
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


def superposition_transform_from_residues_by_name(
    mob_pose, targ_pose, mob_index, targ_index, *args
):
    """
    superimposes mob_pose onto targ_pose by aligning atoms by name

    May provide any number of atom names as long as both are in each residue
    Offset can be provided to superimpose atoms of the same name with an offset
    """
    init_coords = (
        _pyrosetta.rosetta.utility.vector1_numeric_xyzVector_double_t()
<<<<<<< HEAD:dzutils/pyrosetta_utils/phos_binding/enumerate_docks.py
    )
    ref_coords = (
        _pyrosetta.rosetta.utility.vector1_numeric_xyzVector_double_t()
    )
=======
    )
    ref_coords = _pyrosetta.rosetta.utility.vector1_numeric_xyzVector_double_t()
>>>>>>> fix git ignore and imports:phos_binding/enumerate_docks.py
    mob_res = mob_pose.residue(mob_index)
    targ_res = targ_pose.residue(targ_index)

    index_range = range(len(args))

    for index in index_range:
        init_coords.append(mob_res.xyz(args[index]))
        ref_coords.append(targ_res.xyz(args[index]))
    compute_and_apply_superposition_transform(
        mob_pose, init_coords, ref_coords
    )
    return mob_pose


def align_to_tetrahedral_rotation(nub_pose, nub_resnum, rotation=0, *args):
    """
    Align the pose to a rotation of the atoms specified

    only for tetrahedral symmetric rotations
    helper function for enumerating all tetrahedral symmetries

    rotations are numbered from 0-11, with wraparound
    """
    if len(args) > 4:
        raise AssertionError("This function is for tetrahedral symmetry only")

    init_coords = (
        _pyrosetta.rosetta.utility.vector1_numeric_xyzVector_double_t()
<<<<<<< HEAD:dzutils/pyrosetta_utils/phos_binding/enumerate_docks.py
    )
    ref_coords = (
        _pyrosetta.rosetta.utility.vector1_numeric_xyzVector_double_t()
    )
=======
    )
    ref_coords = _pyrosetta.rosetta.utility.vector1_numeric_xyzVector_double_t()
>>>>>>> fix git ignore and imports:phos_binding/enumerate_docks.py
    nub_res = nub_pose.residue(nub_resnum)

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
    compute_and_apply_superposition_transform(
        nub_pose, init_coords, ref_coords
    )
    return nub_pose


# Utility functions
def subset_CA_rmsd(
    pose1, pose2, pose1_residue_subset, pose2_residue_subset, superimpose=True
):
    pose1_residue_selection = _pyrosetta.rosetta.core.select.get_residues_from_subset(
        pose1_residue_subset
    )
    pose2_residue_selection = _pyrosetta.rosetta.core.select.get_residues_from_subset(
        pose2_residue_subset
    )

    assert len(pose1_residue_selection) == len(pose2_residue_selection)

    map_atom_id_atom_id = (
        _pyrosetta.rosetta.std.map_core_id_AtomID_core_id_AtomID()
    )
    for pose1_residue_index, pose2_residue_index in zip(
        pose1_residue_selection, pose2_residue_selection
    ):
        atom_id1 = _pyrosetta.rosetta.core.id.AtomID(
            pose1.residue(pose1_residue_index).atom_index("CA"),
            pose1_residue_index,
        )
        atom_id2 = _pyrosetta.rosetta.core.id.AtomID(
            pose2.residue(pose2_residue_index).atom_index("CA"),
            pose2_residue_index,
        )
        map_atom_id_atom_id[atom_id1] = atom_id2

    if superimpose:
        return _pyrosetta.rosetta.core.scoring.rms_at_corresponding_atoms(
            pose1, pose2, map_atom_id_atom_id, pose1_residue_selection
        )
    else:
        return _pyrosetta.rosetta.core.scoring.rms_at_corresponding_atoms_no_super(
            pose1, pose2, map_atom_id_atom_id, pose1_residue_selection
        )


def index_list_from_residue_name(pose, resname):
    """
    returns a list with the residue indices of the residue with the given name
    """
    residue_name_selector = (
        _pyrosetta.rosetta.core.select.residue_selector.ResidueNameSelector()
    )
    residue_name_selector.set_residue_names(resname)
    vec = residue_name_selector.apply(pose)
    indices = [i for i, is_selected in enumerate(vec, 1) if (is_selected)]
    return indices


def chains_with_resname(pose, resname):
    """
    Takes a pose and a residue name, returns a list of chains with that residue
    """
    residue_name_selector = (
<<<<<<< HEAD:dzutils/pyrosetta_utils/phos_binding/enumerate_docks.py
        ___pyrosetta.rosetta.core.select.residue_selector.ResidueNameSelector()
=======
        _pyrosetta.rosetta.core.select.residue_selector.ResidueNameSelector()
>>>>>>> fix git ignore and imports:phos_binding/enumerate_docks.py
    )
    residue_name_selector.set_residue_names(resname)
    vec = residue_name_selector.apply(pose)
    chains_redundant = [
        pose.chain(i) for i, is_selected in enumerate(vec, 1) if (is_selected)
    ]
    chains = list(set(chains_redundant))
    return chains


def pose_from_chain(pose, chain):
    """
    Returns a pose consisting of just the input chain
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


def apply_tetrahedral_rotation_to_chain(pose, chain, num, resname, atoms):
    """
    Takes a pose and a chain, returns one tetrahedral rotation of that chain

    Rotation is about the residue and atoms given
    """
    mob_pose = pose_from_chain(pose, chain)
    targ_pose = pose_excluding_chain(pose, chain)
    rot_res = index_list_from_residue_name(pose=mob_pose, resname=resname)[0]

    align_to_tetrahedral_rotation(mob_pose, rot_res, num, *atoms)
    targ_pose.append_pose_by_jump(mob_pose, targ_pose.num_jump() + 1)
    return targ_pose


def apply_rotation_generator(chains, rotnum, resname, atoms):
    """
    return a function generator that rotates the chain about the chosen residue
    """
    rotations = ([(chain, num) for num in range(rotnum)] for chain in chains)
    all_rot_combos = _itertools.product(*rotations)
    rot_functions = (
        [
            lambda x, a=op[0], b=op[
                1
            ], c=resname, d=atoms: apply_tetrahedral_rotation_to_chain(
                pose=x, chain=a, num=b, resname=c, atoms=d
            )
            for op in combo
        ]
        for combo in all_rot_combos
    )
    if len(chains) > 1:
        composed_funcs = (
            dzutils.func_tools.compose(*func) for func in rot_functions
        )
        return composed_funcs
    else:
        funcs = (func[0] for func in rot_functions)
        return funcs
