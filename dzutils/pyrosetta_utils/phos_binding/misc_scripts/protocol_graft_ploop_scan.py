#!/usr/bin/env python

import logging

from itertools import permutations

import pyrosetta
from pyrosetta.rosetta.core.select.residue_selector import (
    ResiduePDBInfoHasLabelSelector,
    ChainSelector,
    NotResidueSelector,
    NeighborhoodResidueSelector,
    ResiduePropertySelector,
    ResidueIndexSelector,
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

from dzutils.pyrosetta_utils.anchored_graft_scan import (
    load_fragment_store_from_path,
    graft_generator,
)
from dzutils.pyrosetta_utils.phos_binding.misc_scripts.reslabel_phos_contacts import (
    label_pres,
)
from dzutils.pyrosetta_utils.phos_binding.misc_scripts.get_best_rmsd_phy_to_po4 import (
    best_rmsd_phy_to_po4,
)

root = logging.getLogger()
root.setLevel(logging.DEBUG)

logger = logging.getLogger("graft_scan")
logger.setLevel(logging.DEBUG)


def clash_check(pose, fa_rep_cutoff, *excluded_res, exclude_ligands=True):
    """
    returns if the pose has clashes

    ripped from Brian's faq
    """
    ala_pose_copy = pose.clone()

    scorefxn_fa_rep = pyrosetta.rosetta.core.scoring.ScoreFunctionFactory().create_score_function(
        "none"
    )
    scorefxn_fa_rep.set_weight(pyrosetta.rosetta.core.scoring.fa_rep, 1)

    scorefxn_fa_rep(ala_pose_copy)
    true_sel = (
        pyrosetta.rosetta.core.select.residue_selector.TrueResidueSelector()
    )
    if excluded_res:
        str_ex = ",".join([str(i) for i in excluded_res])
        print(str_ex)
        non_a_res = NotResidueSelector(ResidueIndexSelector(str_ex))
        true_sel = and_compose_residue_selectors(non_a_res, true_sel)
    if exclude_ligands:
        ligand_property = (
            pyrosetta.rosetta.core.chemical.ResidueProperty.LIGAND
        )
        not_ligands_selector = NotResidueSelector(
            ResiduePropertySelector(ligand_property)
        )
        true_sel = and_compose_residue_selectors(
            not_ligands_selector, true_sel
        )
    true_sub = true_sel.apply(ala_pose_copy)

    pyrosetta.rosetta.protocols.toolbox.pose_manipulation.repack_these_residues(
        true_sub, ala_pose_copy, scorefxn_fa_rep, False, "A"
    )
    per_residue_fa_rep = [0]
    for i in range(1, pose.size() + 1):
        fa_rep = ala_pose_copy.energies().residue_total_energy(i)
        per_residue_fa_rep.append(fa_rep)
    fa_array = np.array(per_residue_fa_rep)
    return fa_array < fa_rep_cutoff


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
@click.option("-i", "--inverse-rotamer-radius", default=12)
@click.option(
    "-s",
    "--struct-numbers",
    help="Only attempt grafts on the chosen secondary structures (indexed from 0). Must be given a ',' separated string, sorry, Ill try to fix this later.",
)
@click.option("--same-chain/--no-same-chain", default=True)
@click.option("-e", "--rmsd-threshold", default=0.5, type=float)
@click.option("-l", "--label/--no-label", default=True)
@click.option("-f", "--clash-threshold", default=15)
@click.option("-n", "--n-clash", default=3)
def main(
    pose_pdb,
    fragment_store_path,
    hashmap_path,
    hdf5_store_path,
    dssp_match_types="",
    rosetta_flags_file="",
    allowed_positions=False,
    cart_resl=1,
    ori_resl=15,
    struct_numbers="",
    same_chain=True,
    inverse_rotamer_radius=12,
    rmsd_threshold=0.5,
    label=True,
    clash_threshold=15,
    n_clash=3,
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
        struct_numbers=struct_numbers,
        allowed_depth=14,
    )

    pres = residue_from_name3("PHY")

    successful_grafts = []
    # build grafts
    for graft_num, graft in enumerate(grafts, 1):
        # graft.dump_pdb(f"unscanned_{graft_num}.pdb")
        # graft.dump_pdb("this_one.pdb")
        # Use selectors to get the allowed residues for the scans
        selectors = []
        label_selector = ResiduePDBInfoHasLabelSelector(label)
        # selectors.append(label_selector)
        not_label_selector = NotResidueSelector(label_selector)
        selectors.append(not_label_selector)
        if not same_chain:
            labeled_resnums = [
                resnum
                for resnum, has_label in enumerate(
                    label_selector.apply(graft), 1
                )
                if has_label
            ]
            # print(labeled_resnums)
            excluded_chains = set(
                [str(chain_of(graft, resnum)) for resnum in labeled_resnums]
            )
            chains_for_selector = ",".join(excluded_chains)
            # print(f"chains_for_sele: {chains_for_selector}")
            wrong_chain_selector = ChainSelector(chains_for_selector)
            correct_chain_selector = NotResidueSelector(wrong_chain_selector)
            selectors.append(correct_chain_selector)
        # neighbor_dist = 8.0
        if inverse_rotamer_radius:
            neighbors_selector = NeighborhoodResidueSelector(
                label_selector, inverse_rotamer_radius
            )
            selectors.append(neighbors_selector)
        ligand_property = (
            pyrosetta.rosetta.core.chemical.ResidueProperty.LIGAND
        )
        ligands_selector = ResiduePropertySelector(ligand_property)
        not_ligands = NotResidueSelector(ligands_selector)
        selectors.append(not_ligands)
        rt_start_selector = and_compose_residue_selectors(*selectors)
        rt_start_resnums = [
            resnum
            for resnum, is_allowed in enumerate(
                rt_start_selector.apply(graft), 1
            )
            if is_allowed
        ]

        # print(rt_start_resnums)
        if not rt_start_resnums:
            print("no rt start resnums found")
            continue
        # Get coords in graft
        pres_indices = residues_with_element(graft, "P")
        # [P,OX,OY]
        po4_xyzs = np.array(get_p_bound_xyzs_list(graft, *pres_indices))
        # [N,CA,C]
        stub_xyzs = np.array(get_bb_xyzs(graft, *rt_start_resnums))
        start_resnum_array = np.array(rt_start_resnums)
        # print(stub_xyzs)
        # Do the geometry
        base_xforms = hstub(
            *np.swapaxes(stub_xyzs, 0, 1),
            # stub_xyzs[:, 0, :],
            # stub_xyzs[:, 1, :],
            # stub_xyzs[:, 2, :],
            cen=list(stub_xyzs[:, 1, :]),
        )
        po4xforms = hstub(
            *np.swapaxes(po4_xyzs, 0, 1),
            # po4_xyzs[:, 0, :],
            # po4_xyzs[:, 1, :],
            # po4_xyzs[:, 2, :],
            cen=list(po4_xyzs[:, 0, :]),
        )
        po4_l, base_l = len(po4xforms), len(base_xforms)
        po4_repeat = np.repeat(po4xforms, base_l, axis=0)
        base_tiled = np.tile(base_xforms, (po4_l, 1, 1))
        index_tiled = np.tile(start_resnum_array, (po4_l,))
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
        sorter = np.argsort(store_indices)
        non_u_sorted_inds = store_indices[sorter]
        non_u_sorted_hit_keys = hit_keys[sorter]

        sorted_inds, u_inds, inv = np.unique(
            non_u_sorted_inds, return_index=True, return_inverse=True
        )
        sorted_hit_keys = non_u_sorted_hit_keys[u_inds]

        sorted_hits_resnums = index_tiled[hits_mask][sorter]
        with h5py.File(hdf5_store_path, "r") as store:
            store_keys = store["key_int"][sorted_inds]
            store_chis = store["chis"][sorted_inds]
            # Maybe check rts to make sure they were binned the same way
            # store_rt = store["rt"][store_indices]

        # check to make sure the table used is at least formatted like the hashmap
        if not np.array_equal(sorted_hit_keys, store_keys):
            raise ValueError("store keys are not the same as query keys!")

        # hits were compressed to unique to query the hdf5
        # expand them back to correspond to different resnums
        store_keys = store_keys[inv]
        store_chis = store_chis[inv]
        # rosetta stuff to actually make the hits and dump them

        native_clash = sum(clash_check(pose, clash_threshold) == False)
        for hit_num in range(0, len(sorted_hits_resnums)):
            resnum = sorted_hits_resnums[hit_num]
            chi_set = store_chis[hit_num]
            working = graft.clone()
            working.replace_residue(resnum, pres.clone(), True)
            for i in range(1, len(chi_set) + 1):
                working.set_chi(i, resnum, chi_set[i - 1])
            rmsd = best_rmsd_phy_to_po4(working)
            if rmsd > rmsd_threshold:
                print(f"rmsd exceeds threshold: {rmsd} > {rmsd_threshold}")

                continue
            print(f"rmsd passes threshold: {rmsd} < {rmsd_threshold}")
            label_selector = ResiduePDBInfoHasLabelSelector(label)
            exc = [
                i
                for i, is_label in enumerate(label_selector.apply(working), 1)
                if is_label
            ]
            print(
                f"checking residues {exc} for clashes, threshold fa_rep < {clash_threshold}"
            )

            clashes_found = sum(
                clash_check(working, clash_threshold, *exc) == False
            )
            if clashes_found > n_clash + native_clash:
                print("clash checker reports incompatible graft")
                continue
            print(
                "low enough num clashes found {clashes_found } < {n_clash + native_clash}"
            )
            success_name = f"graft_{graft_num}_resn_{resnum}_hit_{hit_num}.pdb"
            # base_scaf_name = f"ori_{graft_num}_resn_{resnum}_hit_{hit_num}.pdb"
            label_pres(working, "phos_contact")
            label_pres(working, "bb_phos_contact", bb_only=True)
            working.dump_pdb(success_name)
            # graft.dump_pdb(base_scaf_name)
            successful_grafts.append(success_name)
    with open("successes.txt", mode="wt", encoding="utf-8") as f:
        f.write("\n".join(successful_grafts))


if __name__ == "__main__":
    main()
