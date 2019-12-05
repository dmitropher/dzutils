import logging
from itertools import permutations

import pyrosetta
from pyrosetta.rosetta.core.select.residue_selector import (
    ResiduePDBInfoHasLabelSelector,
    ChainSelector,
    NotResidueSelector,
    NeighborhoodResidueSelector,
    ResiduePropertySelector,
)
import click
import h5py

import numpy as np

import getpy as gp

from homog import hstub
from xbin import XformBinner

from dzutils.pyrosetta_utils import (
    run_pyrosetta_with_flags,
    and_compose_residue_selectors,
    residues_with_element,
    bonded_atoms,
    atom_indices_with_element,
    residue_from_name3,
)

from dzutils.pyrosetta_utils.chain_utils import chain_of

from dzutils.pyrosetta_utils.phos_binding.misc_scripts.anchored_graft_scan import (
    load_fragment_store_from_path,
    graft_generator,
)


logger = logging.getLogger("graft_and_scan")

logger.setLevel(logging.DEBUG)


def get_p_bound_xyzs_list(pose, *pres_indices):
    """
    returns a list of xyz sets [phos_atom,bound_1,bound_2]

    dim N X 3 X 3

    [
        [
            [ p_atom_x,y,z]
            [ bond_1_atom_x,y,z]
            [ bond_2_atom_x,y,z]
        ],...
    ]
    """
    xyzs = []
    for ind in pres_indices:
        res = pose.residue(ind)
        for p_atom in atom_indices_with_element(res, "P"):
            bonded = bonded_atoms(res, p_atom)
            bonded_xyzs = [list(res.xyz(atom)) for atom in bonded]
            p_atom_xyz = list(res.xyz(p_atom))
            xyzs.extend(
                [
                    [p_atom_xyz, *atom_pair]
                    for atom_pair in permutations(bonded_xyzs, 2)
                ]
            )
    return xyzs


def get_bb_xyzs(pose, *resnums):
    """
    """
    xyzs = []
    for num in resnums:
        res = pose.residue(num)
        xyzs.append(
            [list(res.xyz("N")), list(res.xyz("CA")), list(res.xyz("C"))]
        )
    return xyzs


@click.command()
@click.argument("pose_pdb", type=click.Path(exists=True))
@click.argument("fragment_store_path", type=click.Path(exists=True))
@click.argument("hashmap_path", type=click.Path(exists=True))
@click.argument("hdf5_store_path", type=click.Path(exists=True))
@click.option("-d", "--dssp-match-types", default="")
@click.option("-r", "--rosetta-flags-file")
@click.option("-c", "--cart-resl", default=1)
@click.option("-o", "--ori-resl", default=15)
@click.option("--get-additional-output/--one-output", default=True)
@click.option("--save/--no-save", default=True)
@click.option(
    "-s",
    "--struct-numbers",
    help="Only attempt grafts on the chosen secondary structures (indexed from 0). Must be given a ',' separated string, sorry, Ill try to fix this later.",
)
def main(
    pose_pdb,
    fragment_store_path,
    hashmap_path,
    hdf5_store_path,
    dssp_match_types="",
    rosetta_flags_file="",
    allowed_positions=False,
    save=True,
    cart_resl=1,
    ori_resl=15,
    get_additional_output=True,
    struct_numbers="",
):
    """
    This program takes a pose and a fragment store and returns alignment graphs

    There are some rules about how stuff is lined up and dssp types yada yada
    """
    if rosetta_flags_file:
        run_pyrosetta_with_flags(rosetta_flags_file)
    else:
        pyrosetta.init()

    fragments = load_fragment_store_from_path(fragment_store_path)
    pose = pyrosetta.pose_from_file(pose_pdb)
    label = "anchored_graft"
    grafts = graft_generator(
        pose,
        fragments,
        dssp_types=dssp_match_types,
        save_intermediate=save,
        get_additional_output=get_additional_output,
        struct_numbers=struct_numbers,
        label=label,
    )

    pres = residue_from_name3("PHY")

    # build grafts
    for graft_num, graft in enumerate(grafts, 1):
        graft.dump_pdb("this_one.pdb")
        # Use selectors to get the allowed residues for the scan
        label_selector = ResiduePDBInfoHasLabelSelector(label)
        labeled_resnums = [
            resnum
            for resnum, has_label in enumerate(label_selector.apply(graft), 1)
            if has_label
        ]
        print(labeled_resnums)
        excluded_chains = set(
            [str(chain_of(graft, resnum)) for resnum in labeled_resnums]
        )
        chains_for_selector = ",".join(excluded_chains)
        print(f"chains_for_sele: {chains_for_selector}")
        wrong_chain_selector = ChainSelector(chains_for_selector)
        correct_chain_selector = NotResidueSelector(wrong_chain_selector)
        neighbor_dist = 8.0
        neighbors_selector = NeighborhoodResidueSelector(
            label_selector, neighbor_dist
        )
        ligand_property = (
            pyrosetta.rosetta.core.chemical.ResidueProperty.LIGAND
        )
        ligands_selector = ResiduePropertySelector(ligand_property)
        not_ligands = NotResidueSelector(ligands_selector)
        rt_start_selector = and_compose_residue_selectors(
            neighbors_selector, correct_chain_selector, not_ligands
        )
        rt_start_resnums = [
            resnum
            for resnum, is_allowed in enumerate(
                rt_start_selector.apply(graft), 1
            )
            if is_allowed
        ]

        print(rt_start_resnums)
        if not rt_start_resnums:
            continue
        # Get coords in graft
        pres_indices = residues_with_element(graft, "P")
        # [P,OX,OY]
        po4_xyzs = np.array(get_p_bound_xyzs_list(graft, *pres_indices))
        # [N,CA,C]
        stub_xyzs = np.array(get_bb_xyzs(graft, *rt_start_resnums))
        start_resnum_array = np.array(rt_start_resnums)
        print(stub_xyzs)
        # Do the geometry
        base_xforms = hstub(
            stub_xyzs[:, 0, :][:, 0, :],
            stub_xyzs[:, 1, :][:, 0, :],
            stub_xyzs[:, 2, :][:, 0, :],
            cen=list(stub_xyzs[:, 1, :][:, 0, :]),
        )
        po4xforms = hstub(
            po4_xyzs[:, 0, :][:, 0, :],
            po4_xyzs[:, 1, :][:, 0, :],
            po4_xyzs[:, 2, :][:, 0, :],
            cen=list(po4_xyzs[:, 0, :][:, 0, :]),
        )
        po4_l, base_l = len(po4xforms), len(base_xforms)
        po4_repeat = np.repeat(po4xforms, base_l, axis=0)
        base_tiled = np.tile(base_xforms, (po4_l, 3, 3))
        index_tiled = np.tile(start_resnum_array, (1,))
        rts = np.linalg.inv(po4_repeat) @ base_tiled

        # make the keys
        binner = XformBinner(cart_resl=cart_resl, ori_resl=ori_resl)
        keys = binner.get_bin_index(rts)

        # check the hashmap
        key_type = np.dtype("i8")
        value_type = np.dtype("i8")
        hashmap = gp.Dict(key_type, value_type)
        hashmap.load(hashmap_path)
        hits_mask = hashmap.contains(keys)
        hit_keys = keys[hits_mask]

        # retrieve the hits
        store_indices = hashmap[hit_keys]
        hits_resnums = index_tiled[hits_mask]
        with h5py.File() as store:
            store_keys = store["key_int"][store_indices]
            store_chis = store["chis"][store_indices]
            # Maybe check rts to make sure they were binned the same way
            # store_rt = store["rt"][store_indices]

        # check to make sure the table used is at least formatted like the hashmap
        if not np.array_equal(hit_keys, store_keys):
            raise ValueError("store keys are not the same as query keys!")

        # rosetta stuff to actually make the hits and dump them
        for hit_num in range(0, len(hits_resnums)):
            resnum = hits_resnums[hit_num]
            chi_set = store_chis[hit_num]
            working = graft.clone()
            working.replace_residue(resnum, pres.clone(), True)
            for i in range(1, len(chi_set) + 1):
                working.set_chi(i, resnum, chi_set[i])
            working.dump_pdb(
                f"graft_{graft_num}_resn_{resnum}_hit_{hit_num}.pdb"
            )


if __name__ == "__main__":
    main()
