from itertools import permutations, combinations, product
import numpy as _np
import pyrosetta as _pyrosetta
from xbin import XformBinner as _xb

from dzutils.pyrosetta_utils.geometry.pose_xforms import (
    get_e2e_xform,
    generate_pose_rt_between_res,
    RotamerRTArray,
)

from dzutils.pyrosetta_utils.geometry.superposition_utilities import (
    super_by_paired_atoms,
    super_by_residues,
)

from dzutils.pyrosetta_utils import (
    residues_with_element,
    hbond_to_residue,
    atom_indices_with_element,
    bonded_atoms,
    or_compose_residue_selectors,
    build_hbond_set,
)
from dzutils.pyrosetta_utils.secstruct import parse_structure_from_dssp

from dzutils.pyrosetta_utils.chain_utils import link_poses


def rt_list_hbond_to_res(pose, resnum, sidechain=False, minimal=False):
    """
    Returns a list of hbonds to the given residue number

    Defaults to backbone hbonds only
    """
    res = pose.residue(resnum)
    return [
        generate_pose_rt_between_res(
            pose,
            hbond.don_res(),
            resnum,
            ("N", "N", "CA", "C"),
            (
                res.atom_name(hbond.acc_atm()),
                res.atom_name(hbond.acc_atm()),
                "CA",
                "C",
            ),
            minimal=minimal,
        )
        for hbond in hbond_to_residue(pose, resnum)
        if sidechain or hbond.don_hatm_is_protein_backbone()
    ]


def p_atoms_in_pose(pose):
    return [
        (atom_i, resnum)
        for resnum in residues_with_element(pose, "P")
        for atom_i in atom_indices_with_element(pose.residue(resnum), "P")
    ]


def get_loop_xform_dicts(pose, num_contacts, *args, loop_chain=1, **kwargs):
    """
    returns a dict with xforms from end to end and from phosphate to end

    Dict also contains an xbin key
    """
    e2e = get_e2e_xform(pose.split_by_chain()[loop_chain])
    xb = _xb(*args, **kwargs)
    p_atoms = [
        (
            p_atom_i,
            [
                atom
                for atom in bonded_atoms(
                    pose.residue(resnum), p_atom_i, name=False
                )
                for hbond in hbond_to_residue(pose, resnum, vec=False)
                if hbond.acc_atm() == atom
            ],
            resnum,
        )
        for p_atom_i, resnum in p_atoms_in_pose(pose)
    ]
    # TODO add reference atoms here
    xform_dicts = [
        {
            "file": pose.pdb_info().name(),
            "key_int": int(xb.get_bin_index(_np.array(e2e))),
            "e2e": _np.array(e2e),
            # Careful here: maybe add a check to see that res 1 is really the beginning of the loop chain?
            "func_to_bb_start": generate_pose_rt_between_res(
                pose, resnum, 1, (p_atom, p_atom, bonded_1, bonded_2)
            ),
            "alignment_atoms": (p_atom, p_atom, bonded_1, bonded_2),
        }
        for p_atom, others, resnum in p_atoms
        if len(others) >= num_contacts
        for bonded_1, bonded_2 in permutations(list(set(others)), 2)
    ]
    return xform_dicts


def loop_res_pairs(pose):
    """
    Returns tuples of all res-res combos within loops
    """
    resnums = [
        i
        for i, is_sel in enumerate(
            or_compose_residue_selectors(
                *[
                    c.generate_residue_selector()
                    for c in parse_structure_from_dssp(pose, "L")
                ]
            ).apply(pose),
            1,
        )
        if is_sel
    ]
    return [(pose.clone(), i, j) for i, j in permutations(resnums, 2) if i > j]


def pair_to_rt_dict(pose, resnum_1, resnum_2, **kwargs):
    """
    Returns a dict with the rt between the given res pair and the pose name

    Can restrict search to a residue index selector string with res_string=
    """

    return {
        "rt": _np.array(
            generate_pose_rt_between_res(pose, resnum_1, resnum_2, **kwargs)
        ),
        "name": pose.pdb_info().name(),
        "start_res": resnum_1,
        "end_res": resnum_2,
    }


def loops_to_rt_dict(pose, plus=0, minus=0):
    """
    Wrapper to hopefully conserve memory

    includes some modifications to loop selection to help with bad hits
    """
    return [
        pair_to_rt_dict(pose, i, j)
        for loop in parse_structure_from_dssp(pose, "L")
        for i, j in permutations(loop.resnum_list(), 2)
        if i < j
        if i != 1
        for start, end in [
            (
                loop.resnum_list(upstream=minus, downstream=plus)[0],
                loop.resnum_list(upstream=minus, downstream=plus)[-1],
            )
        ]
        for i, j in combinations(
            loop.resnum_list(upstream=minus, downstream=plus), 2
        )
        if i > (end + start) / 2 and j > (end + start) / 2
    ]


def pres_bases(residue):
    phos_atoms = atom_indices_with_element(residue, "P")
    return [
        (p_atom_index, p_atom_index, a, b)
        for p_atom_index in phos_atoms
        for (a, b) in permutations(
            bonded_atoms(residue, p_atom_index, name=False), 2
        )
    ]


def phospho_residue_inverse_rotamer_rts(
    residue=None, rotamer_rt_array=None, alignment_atoms=False
):
    """
    Returns the RTs from the phos group to bb

    Generates the RT from all possible stubs generated by the following scheme

    (Phosphorus atom, Phosphorus Atom, p-bound-atom1,p-bound-atom2)
    setting alignment atoms to true changes the return type to a list of tuples:
    (RotamerRTArray,alignment_atoms)


    """
    # pres bases here refers to the atoms that "target_atoms" is based off of
    # in the RotamerRTArray
    # sorry future people for the confusiont
    possible_rt_bases = []

    if rotamer_rt_array is None:

        if residue is None:
            raise ValueError(
                "A RotamerRTArray or a residue must be provided to compute all RTs!"
            )
        else:
            possible_rt_bases = pres_bases(residue)
            rotamer_rt_array = RotamerRTArray(
                residue=residue, target_atoms=pres_bases[0], inverse=True
            )
    else:
        if residue is not None:
            raise ValueError(
                "A RotamerRTArray or a residue must be provided to compute all RTs! But not both"
            )
        else:
            possible_rt_bases = pres_bases(rotamer_rt_array.residue)

    output = []
    for base in possible_rt_bases:
        rotamer_rt_array.set_target_atoms(base)
        output.append(
            _np.array(rotamer_rt_array)
            if not alignment_atoms
            else (_np.array(rotamer_rt_array), base)
        )
    return output


def num_bb_contacts(pose, resnum, atom_i):
    """
    Only returns bb contacts where bb is the donor.

    TODO fix this
    """
    contacts = sum(
        (
            hbond.don_hatm_is_protein_backbone()
            and hbond.acc_atm() in bonded_atoms(pose.residue(resnum), atom_i)
        )
        for hbond in hbond_to_residue(pose, resnum)
    )
    print(contacts)
    return contacts


def replace_p_res_with_phosphate(pose, min_contacts=0):
    """
    returns a copy of the pose with every phosphoresidue converted to phosphate

    PO4 must be a residue in your rosetta instance
    """

    p_res_list = p_atoms_in_pose(pose)
    phosphates = list()
    for atom_i, resnum in p_res_list:
        bb_contacts = num_bb_contacts(pose, resnum, atom_i)
        if bb_contacts < min_contacts:
            print(
                f"res {resnum} at atom {atom_i} has only {bb_contacts} contacts, less than {min_contacts}"
            )
            continue
        p = _pyrosetta.rosetta.core.pose.Pose()
        restype = _pyrosetta.rosetta.core.pose.get_restype_for_pose(p, "PO4")
        res = _pyrosetta.rosetta.core.conformation.Residue(restype, True)
        p.append_residue_by_jump(res, 1)
        p_atom = atom_indices_with_element(res, "P")[0]

        atom_tuples = list(
            zip(
                [p_atom, *bonded_atoms(res, p_atom)[:2]],
                [atom_i, *bonded_atoms(pose.residue(resnum), atom_i)[:2]],
            )
        )
        # print(atom_tuples)
        phosphates.append(
            super_by_paired_atoms(p, pose, 1, resnum, *atom_tuples)
        )
    res_to_remove = list(set([resnum for atom_i, resnum in p_res_list]))
    # print(res_to_remove)
    newp = pose.clone()
    for i, resnum in enumerate(res_to_remove):
        newp.delete_residue_slow(resnum - i)
    for phos in phosphates:
        _pyrosetta.rosetta.core.pose.append_pose_to_pose(newp, phos, True)
    return newp


def phos_bonded_atoms_by_index(residue):
    """
    returns a dict with P atom index : [list of bonded atoms] for each P atom
    """
    return {
        atom_i: bonded_atoms(residue, atom_i, name=False)
        for atom_i in atom_indices_with_element(residue, "P")
    }


def exclude_self_and_non_bb_hbonds(hbond_collection, *acceptor_atoms):
    return [
        b
        for b in hbond_collection
        if b.acc_atm() in acceptor_atoms
        and b.don_hatm_is_protein_backbone()
        and b.don_res() != b.acc_res()
    ]


def exclude_self_hbonds(hbond_collection, *acceptor_atoms):
    return [
        b
        for b in hbond_collection
        if b.acc_atm() in acceptor_atoms and b.don_res() != b.acc_res()
    ]


def get_hbonds(pose):

    hbond_set = build_hbond_set(pose)

    return [
        exclude_self_hbonds(
            hbond_to_residue(pose, resnum, hbond_set=hbond_set, vec=False),
            *acceptor_atoms,
        )
        # And a dict of acceptable acceptor atoms (atoms bound to P)
        # keys are p atoms, values are lists of bound atoms
        for resnum in residues_with_element(pose, "P")
        for atom_i, acceptor_atoms in phos_bonded_atoms_by_index(
            pose.residue(resnum)
        ).items()
    ]


def get_bb_hbonds(pose):

    hbond_set = build_hbond_set(pose)

    return [
        exclude_self_and_non_bb_hbonds(
            hbond_to_residue(pose, resnum, hbond_set=hbond_set, vec=False),
            *acceptor_atoms,
        )
        # And a dict of acceptable acceptor atoms (atoms bound to P)
        # keys are p atoms, values are lists of bound atoms
        for resnum in residues_with_element(pose, "P")
        for atom_i, acceptor_atoms in phos_bonded_atoms_by_index(
            pose.residue(resnum)
        ).items()
    ]


def get_acceptor_res_for_hbond_collection(hbond_collection):
    """
    Meant to make sure that hbonds have been assigned to atom indices sanely
    """
    acceptor_res = list(set([hbond.acc_res() for hbond in hbond_collection]))
    assert (
        len(acceptor_res) == 1
    ), f"An error occured in hbond collection. {len(acceptor_res) } acceptors found in set"
    return acceptor_res[0]


def minimal_fragments_by_contact_number(pose, min_contacts=1, append_factor=0):
    """
    Returns fragment dict with acceptor res and the span between donor residues

    append factor determines number of additional residues appended/prepended

    """

    hbond_collection = get_bb_hbonds(pose)

    pose_size = len(pose.residues)

    append_ranges = [range(append_factor + 1)] * 2
    return [
        {
            "acceptor_res": r,
            "start": min(*contact_set) - x,
            "end": max(*contact_set) + y,
        }
        for hbonds in hbond_collection
        for r in [get_acceptor_res_for_hbond_collection(hbonds)]
        if len(hbonds) >= min_contacts
        for contact_set in combinations(
            [bond.don_res() for bond in hbonds], min_contacts
        )
        for x, y in product(*append_ranges)
        if max(*contact_set) - min(*contact_set) + abs(x + y) < 11
        if min(*contact_set) - x > 1 and max(*contact_set) + y < pose_size
        if min(*contact_set) - x > r or r > max(*contact_set) + y
    ]


def minimal_fragments_by_secondary_structure(
    pose,
    *struct_types,
    min_contacts=1,
    proximity=5,
    lazy=False,
    append_factor=0,
    sidechain=False,
):
    """
    returns fragments with adjacent secondary structure to contacts

    struct types are dssp string names of structure

    blank struct_types gives all
    """
    bb_hbonds = get_bb_hbonds(pose) if not sidechain else get_hbonds(pose)

    pose_size = len(pose.residues)
    contacts = [
        (r, min(*contact_set), max(*contact_set))
        for hbonds in bb_hbonds
        for r in [get_acceptor_res_for_hbond_collection(hbonds)]
        if len(hbonds) >= min_contacts
        for contact_set in combinations(
            [bond.don_res() for bond in hbonds], min_contacts
        )
        if max(*contact_set) - min(*contact_set) < 11
        if min(*contact_set) > 1 and max(*contact_set) < pose_size
    ]

    struct_list = parse_structure_from_dssp(pose, *struct_types)
    sec_struct_by_start_pos = {
        struct.start_pos: struct for struct in struct_list
    }
    sec_struct_by_end_pos = {struct.end_pos: struct for struct in struct_list}

    out_list = list()

    for resnum, start_contact, end_contact in contacts:
        end_struct = int()
        for i in range(end_contact, end_contact + proximity + 1):
            if i in sec_struct_by_start_pos:
                end_struct = sec_struct_by_start_pos[i].end_pos
                if lazy:
                    break
        if not end_struct:
            end_struct = min(end_contact + append_factor, pose_size)
        start_struct = int()
        for i in range(start_contact, start_contact - proximity - 1, -1):
            if i in sec_struct_by_end_pos:
                start_struct = sec_struct_by_end_pos[i].start_pos
                if lazy:
                    break
        if not start_struct:
            start_struct = max(start_contact - append_factor, 1)
        out_list.append(
            {"acceptor_res": resnum, "start": start_struct, "end": end_struct}
        )
    return out_list


##### FORGIVE ME FUTURE DMITRI, I'M MOVING THIS STUFF INTO OUR PACKAGE ########
##### LOVE, PAST DMITRI #####


def super_and_insert_pose(
    start, end, pose, insertion, insertion_start, insert_chain=0
):
    """
    By default, inserts the whole pose as inflexible insert_chain

    If an insert chain is given, that chain is inserted as an inflexible insert,
    the rest go in as chains.
    """
    moved_insertion = super_by_residues(
        insertion, pose, insertion_start, start
    )
    newp = _pyrosetta.rosetta.core.pose.Pose()
    newp.detached_copy(pose)

    _pyrosetta.rosetta.protocols.grafting.delete_region(
        newp, start, end
    )  # - 1

    # print(start, end)

    _pyrosetta.rosetta.core.pose.remove_variant_type_from_pose_residue(
        newp,
        _pyrosetta.rosetta.core.chemical.UPPER_TERMINUS_VARIANT,
        start - 1,
    )
    if insert_chain:
        chains = moved_insertion.split_by_chain()
        chain = chains[insert_chain]
        out_pose = _pyrosetta.rosetta.protocols.grafting.insert_pose_into_pose(
            newp, chain, start - 1, start, False
        )
        for i, c in enumerate(chains, 1):
            if i != insert_chain:
                _pyrosetta.rosetta.core.pose.append_pose_to_pose(
                    out_pose, c, True
                )
    else:
        out_pose = _pyrosetta.rosetta.protocols.grafting.insert_pose_into_pose(
            newp, moved_insertion, start - 1, start, False
        )
    # out_pose = insert_pose(newp, moved_insertion, start - 1)
    return out_pose


def replace_res_from_pose(pose, replacement, index, replacement_index):
    replacement_copy = replacement.clone()
    # super_by_residues(replacement_copy,pose,replacement_index,index)
    _pyrosetta.rosetta.core.pose.remove_variant_type_from_pose_residue(
        replacement_copy,
        _pyrosetta.rosetta.core.chemical.UPPER_TERMINUS_VARIANT,
        replacement_index,
    )
    _pyrosetta.rosetta.core.pose.remove_variant_type_from_pose_residue(
        replacement_copy,
        _pyrosetta.rosetta.core.chemical.LOWER_TERMINUS_VARIANT,
        replacement_index,
    )
    pose.replace_residue(
        index, replacement_copy.residue(replacement_index), True
    )
    # _pyrosetta.rosetta.core.pose.correctly_add_cutpoint_variants(pose)
    for num in range(1, pose.num_chains() + 1):
        if index == pose.chain_begin(num):
            _pyrosetta.rosetta.core.pose.add_variant_type_to_pose_residue(
                pose,
                _pyrosetta.rosetta.core.chemical.LOWER_TERMINUS_VARIANT,
                index,
            )
            pose.conformation().chains_from_termini()
        if index == pose.chain_end(num):
            _pyrosetta.rosetta.core.pose.add_variant_type_to_pose_residue(
                pose,
                _pyrosetta.rosetta.core.chemical.UPPER_TERMINUS_VARIANT,
                index,
            )
            pose.conformation().chains_from_termini()
    return pose


def graft_and_dump_pdb(
    pose, insertion, start, end, insertion_start, dump_path, insert_chain=0
):
    """
    """
    new_pose = super_and_insert_pose(
        start, end, pose, insertion, insertion_start
    )
    new_pose.dump_pdb(dump_path)
    # del new_pose
    # if not dumped:
    #     raise RuntimeError(f"Unable to dump pdb at specfied path: {dump_path}")
    return dump_path
