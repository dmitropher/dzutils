#! /usr/bin/env python
from itertools import permutations

import click
import numpy as np
import pyrosetta

from pyrosetta.rosetta.core.scoring.dssp import Dssp

import xarray as xr

from dzutils.pyrosetta_utils import (
    residues_with_element,
    residue_type_from_name3,
    run_pyrosetta_with_flags,
    hbond_to_residue,
    atom_indices_with_element,
    bonded_atoms,
    build_hbond_set,
)

from dzutils.pyrosetta_utils.phos_binding import replace_p_res_with_phosphate
from dzutils.pyrosetta_utils.geometry.homog import homog_from_residue


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


def exclude_self_and_bb_hbonds(hbond_collection, *acceptor_atoms):
    return [
        b
        for b in hbond_collection
        if b.acc_atm() in acceptor_atoms
        and not b.don_hatm_is_protein_backbone()
        and b.don_res() != b.acc_res()
    ]


def get_bb_hbonds(pose):

    hbond_set = build_hbond_set(pose)

    return [
        b
        for resnum in residues_with_element(pose, "P")
        for atom_i, acceptor_atoms in phos_bonded_atoms_by_index(
            pose.residue(resnum)
        ).items()
        for b in exclude_self_and_non_bb_hbonds(
            hbond_to_residue(pose, resnum, hbond_set=hbond_set, vec=False),
            *acceptor_atoms,
        )
    ]


def get_sidechain_hbonds(pose):

    hbond_set = build_hbond_set(pose)

    return [
        b
        for resnum in residues_with_element(pose, "P")
        for atom_i, acceptor_atoms in phos_bonded_atoms_by_index(
            pose.residue(resnum)
        ).items()
        for b in exclude_self_and_bb_hbonds(
            hbond_to_residue(pose, resnum, hbond_set=hbond_set, vec=False),
            *acceptor_atoms,
        )
    ]


def nonblank_lines(f):
    for l in f:
        line = l.rstrip()
        if line:
            yield line


def has_bb(pose, resnum):
    return (
        pose.residue(resnum).has("CA")
        and pose.residue(resnum).has("N")
        and pose.residue(resnum).has("C")
    )


def check_fragment(
    pose,
    start,
    end,
    pres,
    bb_hbonds,
    sidechain_hbonds,
    bb_contacts=2,
    s_contacts=1,
):
    """
    Returns True if the fragment has enough bb/sc hbonds and pres is the acceptor
    """
    # validate that start and end have N,CA,C
    if not (has_bb(pose, start) and has_bb(pose, end)):
        return False

    sc_hbond_list = [
        b
        for b in sidechain_hbonds
        if b.acc_res() == pres and start <= b.don_res() <= end
    ]
    bb_hbond_list = [
        b
        for b in bb_hbonds
        if b.acc_res() == pres and start <= b.don_res() <= end
    ]

    if len(sc_hbond_list) < s_contacts or len(bb_hbond_list) < bb_contacts:
        return False
    else:
        return True


def get_fragments(pose, min_len=3, max_len=7, bb_contacts=2, s_contacts=1):
    """
    Returns a list of tuples with phos contact stretches

    begin, end, p_res ,begin_struct,end_struct
    """
    print("getting fragments from pose")
    bb_hbonds = get_bb_hbonds(pose)
    sidechain_hbonds = get_sidechain_hbonds(pose)
    pres_list = residues_with_element(pose, "P")
    structs = Dssp(pose).get_dssp_secstruct()
    plen = len(pose.residues)
    return [
        (start, end, pres, structs[start - 1], structs[end - 1])
        for start in range(1, len(pose.residues) + 1)
        for end in range(
            min(plen, start + min_len - 1), min(plen + 1, start + max_len)
        )
        for pres in pres_list
        if check_fragment(
            pose,
            start,
            end,
            pres,
            bb_hbonds,
            sidechain_hbonds,
            bb_contacts=bb_contacts,
            s_contacts=s_contacts,
        )
    ]


@click.command()
@click.argument("db_out_path", nargs=1)
@click.argument("pdb_output_path", nargs=1)
@click.argument("pdb_list", nargs=-1)
@click.option(
    "-r",
    "--rosetta-flags-file",
    "rosetta_flags_file",
    default="",
    show_default=True,
)
@click.option(
    "-b", "--min-bb-contacts", "min_bb_contacts", default=2, show_default=True
)
@click.option(
    "-s",
    "--min-sidechain-contacts",
    "min_sidechain_contacts",
    default=0,
    show_default=True,
)
@click.option(
    "-f",
    "--fragment-size",
    "fragment_size",
    type=(int, int),
    default=(3, 7),
    show_default=True,
    help="Allowed fragment size",
)
def main(
    db_out_path,
    pdb_output_path,
    pdb_list,
    rosetta_flags_file="",
    min_bb_contacts=2,
    min_sidechain_contacts=0,
    fragment_size=(3, 7),
    # secstruct_anchor_type="H",
):
    ""
    run_pyrosetta_with_flags(rosetta_flags_file)

    full_ds = None
    for pdb in pdb_list:
        try:
            pdb_pose = pyrosetta.pose_from_pdb(pdb)
        except RuntimeError as e:
            print(e)
            print(f"unable to load: {pdb}")
            continue
        # strip all phos moiety containing residues, replace phos with PO4 placeholders
        print(pdb_pose.sequence())
        print(len(pdb_pose.residues))
        pose = replace_p_res_with_phosphate(pdb_pose)
        new_name = (
            f'{".".join(pdb.split("/")[-1].split(".")[:-1])}_with_PO4.pdb'
        )
        pdb_path = f"{pdb_output_path}/{new_name}"
        pose.dump_pdb(pdb_path)

        # find contact stretches with the minimums met, no bigger than the max fragment size

        rtype = residue_type_from_name3("PO4")
        p_atom_i = [
            i
            for i in range(1, rtype.natoms() + 1)
            if rtype.atom_type(i).element() == "P"
        ]
        p_atoms = [
            [
                (pair[0], pi, pair[1])
                for pair in permutations(rtype.bonded_neighbor(pi), 2)
            ]
            for pi in p_atom_i
        ]
        phos_atom_array = np.array(p_atoms)
        bb_atoms = ("CA", "N", "CA", "C")
        bb_atoms_array = np.array([bb_atoms])
        anchors = ["begin", "end"]

        bb_atom_da = xr.DataArray(
            bb_atoms_array, dims=("bb_atom_set", "bb_atoms")
        )
        p_atom_da = xr.DataArray(
            phos_atom_array,
            dims=("phos_center_atom", "p_atom_set", "phos_atoms"),
            coords={"phos_center_atom": p_atom_i},
        )

        frags = get_fragments(
            pose,
            min_len=fragment_size[0],
            max_len=fragment_size[1],
            bb_contacts=min_bb_contacts,
            s_contacts=min_sidechain_contacts,
        )

        for begin, end, p_res, begin_struct, end_struct in frags:
            print(
                f"fragment: begin: {begin} end: { end} p_res:{ p_res} begin_struct: {begin_struct} end_struct: {end_struct}"
            )
            phos_stubs = np.array(
                [
                    homog_from_residue(
                        pose.residue(p_res),
                        center_atom=a2,
                        atom1=a1,
                        atom2=a2,
                        atom3=a3,
                    )
                    for root in p_atoms
                    for a1, a2, a3 in root
                ]
            )
            set_len = phos_stubs.shape[0]
            bb_hstub_begin = homog_from_residue(pose.residue(begin))
            bb_hstub_end = homog_from_residue(pose.residue(end))
            tiled_bb_begin = np.tile(bb_hstub_begin, (set_len, 1, 1))
            tiled_bb_end = np.tile(bb_hstub_end, (set_len, 1, 1))
            tiled_bb = np.concatenate(
                (tiled_bb_begin[np.newaxis, :], tiled_bb_end[np.newaxis, :])
            )
            phos_stubs_along_anchor = np.concatenate(
                [phos_stubs[np.newaxis, :]] * 2
            )
            inv_bb = np.linalg.inv(tiled_bb)
            anchor_to_p_rt = inv_bb @ phos_stubs_along_anchor

            struct_da = xr.DataArray(
                np.array([[begin_struct], [end_struct]]),
                dims=("anchor_type", "struct_type"),
                coords={"anchor_type": anchors},
            )
            residues = ["begin", "end", "p_res"]
            res_da = xr.DataArray(
                np.array([[begin], [end], [p_res]]),
                dims=("anchor_res", "resnum"),
                coords={"anchor_res": residues},
            )
            pdb_da = xr.DataArray([pdb_path], dims=("pdb_path"))

            pdb_ds = xr.Dataset({"pdb_paths": pdb_da})
            res_ds = xr.Dataset({"res_pairs": res_da})
            anchor_ds = xr.Dataset({"anchor_struct": struct_da})
            bb_atom_ds = xr.Dataset({"bb_atom_sets": bb_atom_da})
            meta_ds = xr.Dataset({"p_atom_sets": p_atom_da})
            meta_ds = meta_ds.merge(bb_atom_ds)
            meta_ds = meta_ds.merge(anchor_ds)
            meta_ds = meta_ds.merge(res_ds)
            meta_ds = meta_ds.merge(pdb_ds)

            xforms_da = xr.DataArray(
                anchor_to_p_rt,
                dims=("anchor_type", "rts", "xform_row", "xform_col"),
                coords={"anchor_type": anchors},
            )

            xform_ds = xr.Dataset({"rt": xforms_da})
            entry = xform_ds.merge(meta_ds)
            if full_ds:
                full_ds = xr.concat((full_ds, entry), "pdb_path")
            else:
                full_ds = entry
            print(full_ds)
    full_ds.to_netcdf(f"{db_out_path}/p_contacts.nc")


if __name__ == "__main__":
    main()
