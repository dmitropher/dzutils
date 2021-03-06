import pyrosetta
import pyrosetta.rosetta as _pyr
from pyrosetta import pose_from_file as _from_file
from dzutils.pdb_file_utils import pdb_files_in_dir as _pfd
from dzutils.util import read_flag_file


def parse_label(pose, label):
    sel = pyrosetta.rosetta.core.select.residue_selector.ResiduePDBInfoHasLabelSelector(
        label
    )
    return [i for i, has_label in enumerate(sel.apply(pose), 1) if has_label]


def safe_load_pdbs(pdbs):
    for pdb in pdbs:
        try:
            yield pyrosetta.pose_from_pdb(pdb)
        except RuntimeError as e:
            print(e)
            print(f"unable to load: {pdb}")
            continue


def atom_coords(pose, *selected):
    coords = pyrosetta.rosetta.utility.vector1_numeric_xyzVector_double_t()
    if selected:
        for i, j in selected:
            xyz = pose.residue(i).xyz(j)
            coords.append(xyz)
        return coords
    for residue in pose.residues:
        for atom in residue.atoms():
            coords.append(atom.xyz())
    return coords


def fuzzy_trim(pose, end, trim_amount, buffer=0):
    """
    performs delete_residue_range_slow from end - (1 through trim_amount) to end

    setting buffer also allows trims that stop short of end (shorter trims)
    """
    trimmed = []
    for i in range(end - trim_amount, end):
        for j in range(1, buffer + 1):
            if i < 1:
                continue
            elif end - j >= i:
                to_trim = pose.clone()
                to_trim.delete_residue_range_slow(i, end - j)
                trimmed.append(to_trim)
    return trimmed


def residue_type_from_name3(name, variant=None):
    """
    Returns a new residue object from the current ResidueTypeSet

    Optionally specify VariantType to get the appropriate variant

    This function is basically just a wrapper to improve code readability
    """
    chemical_manager = (
        pyrosetta.rosetta.core.chemical.ChemicalManager.get_instance()
    )
    rts = chemical_manager.residue_type_set(
        pyrosetta.rosetta.core.chemical.TypeSetMode.FULL_ATOM_t
    )
    if variant:
        return rts.get_residue_type_with_variant_added(
            rts.name_map(name), variant
        )
    else:
        return rts.name_map(name)


def make_pack_rotamers_mover(pose, score_function, *task_ops):
    """
    Takes a pose, score function, and any number of task operations, and returns a PackRotamers mover
    """
    task_factory = pyrosetta.rosetta.core.pack.task.TaskFactory()
    for task_op in task_ops:
        task_factory.push_back(task_op)
    packer_task = task_factory.create_task_and_apply_taskoperations(pose)
    return pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(
        score_function, packer_task
    )


def residue_from_name3(name, variant=None):
    return _pyr.core.conformation.Residue(
        residue_type_from_name3(name, variant), True
    )


def run_pyrosetta_with_flags(flags_file_path, mute=False):
    if not flags_file_path:
        pyrosetta.init("-mute all " if mute else "", silent=mute)
        return
    flags = read_flag_file(flags_file_path)
    flags_str = " ".join(flags.replace("\n", " ").split())
    pyrosetta.init(
        f"-mute all {flags_str}" if mute else flags_str, silent=mute
    )


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


def get_rotamer_pose_from_residue_type(residue_type, *chis):
    pose = pyrosetta.rosetta.core.pose.Pose()
    residue = _pyr.core.conformation.Residue(residue_type, True)
    pose.append_residue_by_bond(residue)
    for i, chi in enumerate(chis, 1):
        pose.set_chi(i, 1, chi)
    rotamer_pose = pose
    return rotamer_pose


def get_rotamer_pose_from_name(residue_type, *chis, variant=None):
    pose = pyrosetta.rosetta.core.pose.Pose()
    residue = residue_from_name3(residue_type, variant=variant)
    pose.append_residue_by_bond(residue)
    for i, chi in enumerate(chis, 1):
        pose.set_chi(i, 1, chi)
    rotamer_pose = pose
    return rotamer_pose


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


def poses_from_pdb_dir(dir, filename_condition=(lambda x: True)):
    """
    returns a generator that creates poses from the given directory

    A condition can be given which will be run on the filename. If it returns
    True the file will be loaded to pose.
    """
    return (_from_file(file) for file in _pfd(dir) if filename_condition(file))


def residues2pose(pose):
    """
    Returns a generator that turns each residue into a pose
    """
    for res in pose.residues:
        new_pose = _pyr.core.pose.Pose()
        new_pose.append_residue_by_bond(res, True)
        yield new_pose


def add_reslabel(pose, label, *resnums):
    """
    Applies a residue label to the given contacts
    """
    for i in resnums:
        pose.pdb_info().add_reslabel(i, label)


def and_compose_residue_selectors(*args):
    """
    Takes a list of residue selectors of arbitrary size and composes them with AND, returns the AND selector
    """
    andSelector = None
    for a in args:
        andSelector = _pyr.core.select.residue_selector.AND_combine(
            andSelector, a
        )

    return andSelector


def or_compose_residue_selectors(*args):
    """
    Takes a list of residue selectors of arbitrary size and composes them with OR, returns the OR selector
    """
    orSelector = None
    for a in args:
        orSelector = _pyr.core.select.residue_selector.OR_combine(
            orSelector, a
        )

    return orSelector


def exclude_by_label_residue_selector(label):
    """
    Takes a label string and generates a residue selector which selects all the other residues
    """

    labelSelector = (
        _pyr.core.select.residue_selector.ResiduePDBInfoHasLabelSelector()
    )
    labelSelector.set_label(label)
    notSelector = _pyr.core.select.residue_selector.NotResidueSelector()
    notSelector.set_residue_selector(labelSelector)
    return notSelector


def build_hbond_set(
    pose,
    exclude_DNA_DNA=True,
    exclude_intra_res_protein=True,
    exclude_intra_res_RNA=False,
    put_intra_into_total=False,
    exclude_self_hbonds=True,
    use_hb_env_dep=True,
    use_hb_env_dep_DNA=True,
    smooth_hb_env_dep=True,
    bb_donor_acceptor_check=True,
    decompose_bb_hb_into_pair_energies=False,
    use_sp2_chi_penalty=True,
    sp2_BAH180_rise=0.75,
    sp2_outer_width=0.357,
    measure_sp3acc_BAH_from_hvy=True,
    fade_energy=True,
    exclude_ether_oxygens=False,
    Mbhbond=False,
    mphbond=False,
    hbond_energy_shift=0.0,
    length_dependent_srbb=False,
    length_dependent_srbb_lowscale=0.5,
    length_dependent_srbb_highscale=2.0,
    length_dependent_srbb_minlength=4,
    length_dependent_srbb_maxlength=17,
    water_hybrid_sf=False,
):
    """
    return hbondset with some "sensible" defaults
    """
    hbond_set = _pyr.core.scoring.hbonds.HBondSet(pose, False)
    options = hbond_set.hbond_options()
    hbond_set = _pyr.core.scoring.hbonds.HBondSet()
    options = hbond_set.hbond_options()
    options.exclude_DNA_DNA(exclude_DNA_DNA)
    options.exclude_intra_res_protein(exclude_intra_res_protein)
    options.exclude_intra_res_RNA(exclude_intra_res_RNA)
    options.put_intra_into_total(put_intra_into_total)
    options.exclude_self_hbonds(exclude_self_hbonds)
    options.use_hb_env_dep(use_hb_env_dep)
    options.use_hb_env_dep_DNA(use_hb_env_dep_DNA)
    options.smooth_hb_env_dep(smooth_hb_env_dep)
    options.bb_donor_acceptor_check(bb_donor_acceptor_check)
    options.decompose_bb_hb_into_pair_energies(
        decompose_bb_hb_into_pair_energies
    )
    options.use_sp2_chi_penalty(use_sp2_chi_penalty)
    options.sp2_BAH180_rise(sp2_BAH180_rise)
    options.sp2_outer_width(sp2_outer_width)
    options.measure_sp3acc_BAH_from_hvy(measure_sp3acc_BAH_from_hvy)
    options.fade_energy(fade_energy)
    options.exclude_ether_oxygens(exclude_ether_oxygens)
    options.Mbhbond(Mbhbond)
    options.mphbond(mphbond)
    options.hbond_energy_shift(hbond_energy_shift)
    options.length_dependent_srbb(length_dependent_srbb)
    options.length_dependent_srbb_lowscale(length_dependent_srbb_lowscale)
    options.length_dependent_srbb_highscale(length_dependent_srbb_highscale)
    options.length_dependent_srbb_minlength(length_dependent_srbb_minlength)
    options.length_dependent_srbb_maxlength(length_dependent_srbb_maxlength)
    options.water_hybrid_sf(water_hybrid_sf)
    pose.update_residue_neighbors()
    hbond_set.setup_for_residue_pair_energies(pose, False, False)
    return hbond_set


def hbond_to_residue(pose, resnum, hbond_set=None, vec=False):
    """
    Returns a list of all the hbonding residues to the given residue

    vec=True returns a rosetta vector instead (marginally faster)
    """
    if not hbond_set:
        hbond_set = build_hbond_set(pose)

    if vec:
        return hbond_set.residue_hbonds(resnum)
    else:
        return [b for b in hbond_set.residue_hbonds(resnum)]


def residues_with_element(pose, *elements):
    """
    Returns a list of resnums where element is in some atom of the residue
    """
    return list(
        set(
            [
                i
                for i, res in enumerate(pose.residues, 1)
                for j in range(1, len(res.atoms()) + 1)
                if res.atom_type(j).element() in elements
            ]
        )
    )


def atom_indices_with_element(residue, element):
    """
    Returns a list of atom_indices for the given element
    """
    return [
        i
        for i in range(1, len(residue.atoms()) + 1)
        if residue.atom_type(i).element() == element
    ]


def bonded_atoms(residue, index, name=False):
    """
    Returns a list of atom indices where those atoms have a bond to index

    If name is True, returns a list of atom names
    """
    return [
        (residue.atom_name(j) if name else j)
        for j in residue.bonded_neighbor(index)
    ]
