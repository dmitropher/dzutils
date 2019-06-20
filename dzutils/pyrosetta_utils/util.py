import pyrosetta.rosetta as _pyr
from pyrosetta import pose_from_file as _from_file
from dzutils.pdb_file_utils import pdb_files_in_dir as _pfd


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


def hbond_to_residue(pose, resnum, vec=False):
    """
    Returns a list of all the hbonding residues to the given residue

    vec=True returns a rosetta vector instead (marginally faster)
    """
    hbond_set = _pyr.core.scoring.hbonds.HBondSet(pose, False)
    options = hbond_set.hbond_options()
    # print(
    #     "default options",
    #     options.exclude_DNA_DNA(),
    #     options.exclude_intra_res_protein(),
    #     options.exclude_intra_res_RNA(),
    #     options.put_intra_into_total(),
    #     options.exclude_self_hbonds(),
    #     options.use_hb_env_dep(),
    #     options.use_hb_env_dep_DNA(),
    #     options.smooth_hb_env_dep(),
    #     options.bb_donor_acceptor_check(),
    #     options.decompose_bb_hb_into_pair_energies(),
    #     options.use_sp2_chi_penalty(),
    #     options.sp2_BAH180_rise(),
    #     options.sp2_outer_width(),
    #     options.measure_sp3acc_BAH_from_hvy(),
    #     options.fade_energy(),
    #     options.exclude_ether_oxygens(),
    #     options.Mbhbond(),
    #     options.mphbond(),
    #     options.hbond_energy_shift(),
    #     options.length_dependent_srbb(),
    #     options.length_dependent_srbb_lowscale(),
    #     options.length_dependent_srbb_highscale(),
    #     options.length_dependent_srbb_minlength(),
    #     options.length_dependent_srbb_maxlength(),
    #     options.water_hybrid_sf(),
    # )
    hbond_set = _pyr.core.scoring.hbonds.HBondSet()
    options = hbond_set.hbond_options()
    options.exclude_DNA_DNA(True)
    options.exclude_intra_res_protein(True)
    options.exclude_intra_res_RNA(False)
    options.put_intra_into_total(False)
    options.exclude_self_hbonds(True)
    options.use_hb_env_dep(True)
    options.use_hb_env_dep_DNA(True)
    options.smooth_hb_env_dep(True)
    options.bb_donor_acceptor_check(True)
    options.decompose_bb_hb_into_pair_energies(False)
    options.use_sp2_chi_penalty(True)
    options.sp2_BAH180_rise(0.75)
    options.sp2_outer_width(0.357)
    options.measure_sp3acc_BAH_from_hvy(True)
    options.fade_energy(True)
    options.exclude_ether_oxygens(False)
    options.Mbhbond(False)
    options.mphbond(False)
    options.hbond_energy_shift(0.0)
    options.length_dependent_srbb(False)
    options.length_dependent_srbb_lowscale(0.5)
    options.length_dependent_srbb_highscale(2.0)
    options.length_dependent_srbb_minlength(4)
    options.length_dependent_srbb_maxlength(17)
    options.water_hybrid_sf(False)
    pose.update_residue_neighbors()
    hbond_set.setup_for_residue_pair_energies(pose, False, False)
    # options.length_dependent_srbb(False)

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
