import sys
import itertools as it

import pyrosetta

import click

from dzutils.pyrosetta_utils import (
    residues_with_element,
    run_pyrosetta_with_flags,
    hbond_to_residue,
    atom_indices_with_element,
    bonded_atoms,
    build_hbond_set,
)
from dzutils.pyrosetta_utils.chain_utils import link_poses

from dzutils.pyrosetta_utils.secstruct import parse_structure_from_dssp


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
    assert len(acceptor_res) == 1, "An error occured in hbond collection."
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
        for contact_set in it.combinations(
            [bond.don_res() for bond in hbonds], min_contacts
        )
        for x, y in it.product(*append_ranges)
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
):
    """
    returns fragments with adjacent secondary structure to contacts
    """
    bb_hbonds = get_bb_hbonds(pose)
    pose_size = len(pose.residues)
    contacts = [
        (r, min(*contact_set), max(*contact_set))
        for hbonds in bb_hbonds
        for r in [get_acceptor_res_for_hbond_collection(hbonds)]
        if len(hbonds) >= min_contacts
        for contact_set in it.combinations(
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


# @click.command()
# @click.argument("pose_pdb", type=click.Path(exists=True))
# @click.option("-r", "--reslabel",default="")
# @click.option("j","--just-label/--full-run", default=False)
def main():  # pose_pdb,reslabel="",just_label=False):
    ploop_flags_file = "/home/dzorine/phos_binding/pilot_runs/loop_grafting/initial_testing/misc_files/p_ligand.flags"
    run_pyrosetta_with_flags(ploop_flags_file)
    # get pose
    pose = pyrosetta.pose_from_pdb(sys.argv[1])
    # get_outdir
    outdir = sys.argv[2]
    name = pose.pdb_info().name().split("/")[-1].split(".pdb")[0]
    # Scan pose for phosphorus containing residues
    # extract hbonds to these residues where:
    # vec False to use list rather than rosetta vector

    # no support for different append/prepend values
    append_factor = 3  # +- 0-3 residues

    num_contacts = int(sys.argv[3])

    # Extract contiguous loops with these contacts:
    #   - for each bb pair, check if intervening sequence is under 10 res
    #   - prepend and apppend +- 3 residues on each side if they exist

    for d in minimal_fragments_by_contact_number(
        pose, min_contacts=num_contacts, append_factor=append_factor
    ):

        r = d["acceptor_res"]
        s = d["start"]
        e = d["end"]
        # if reslabel:
        #     reslabel_contacts(pose,reslabel,)
        p = link_poses(
            pyrosetta.rosetta.protocols.grafting.return_region(
                pose.clone(), r, r
            ),
            pyrosetta.rosetta.protocols.grafting.return_region(pose, s, e),
            rechain=True,
        )
        p.dump_pdb(
            f"{outdir}/{name}_{num_contacts}-contacts_phos-{r}_frag_{s}-{e}.pdb"
        )


if __name__ == "__main__":
    main()
