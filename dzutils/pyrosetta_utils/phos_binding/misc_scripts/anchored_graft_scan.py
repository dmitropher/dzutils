import json
import sys, os
import logging

import pyrosetta
import click

from pyrosetta.rosetta.numeric import xyzVector_double_t as rosetta_vector
from pyrosetta.rosetta.protocols.toolbox.pose_manipulation import (
    rigid_body_move,
)
from pyrosetta.rosetta.core.select.residue_selector import TrueResidueSelector

from dzutils.pyrosetta_utils.secstruct.structure import (
    parse_structure_from_dssp,
)
from dzutils.pyrosetta_utils.geometry.superposition_utilities import (
    super_resi_by_bb,
)
from dzutils.pyrosetta_utils.chain_utils import (
    link_poses,
    chain_of,
    run_direct_segment_lookup,
)

from dzutils.pyrosetta_utils import run_pyrosetta_with_flags

from dzutils.pyrosetta_utils.geometry.homog import (
    homog_super_transform_from_residues,
    np_homog_to_rosetta_rotation_translation,
)

logger = logging.getLogger("anchored_graft")

logger.setLevel(logging.DEBUG)


def load_pdbs(pdb_list):
    """
    """
    for path in pdb_list:
        try:
            yield pyrosetta.pose_from_file(path)
        except Exception as e:
            print(e)
            print(f"Unable to load pose from: {path}")


def load_fragment_store_from_path(path):
    """
    Expects a json list of paths
    """
    store = []
    with open(path, "r") as f:
        store = json.load(f)
    loaded_store = [*load_pdbs(store)]
    if not loaded_store:
        raise RuntimeError("Unable to load any pdbs from fragment store!")
    return loaded_store


def get_graft_sites(pose_structs):
    """
    """
    return [
        p
        for struct in pose_structs
        for p in (struct.start_pos, struct.end_pos)
    ]


def super_and_graft(
    pose, fragment, insert_align, pose_align, *ligands, margin=5, label=""
):
    """
    """

    super_resi_by_bb(fragment, pose, insert_align, pose_align)
    subpose_after_graft = (
        pyrosetta.rosetta.protocols.grafting.return_region(
            pose.clone(), pose_align + 1, len(pose.residues)
        )
        if len(pose.residues) > pose_align
        else pyrosetta.rosetta.core.pose.Pose()
    )
    subpose_before_graft = pyrosetta.rosetta.protocols.grafting.return_region(
        pose.clone(), 1, pose_align - 1
    )
    # preceding_loop = parse_structure_from_dssp(subpose_before_graft, "L")[
    #     -1
    # ]
    # first_struct = parse_structure_from_dssp(pose, *"EH")[0]
    if len(subpose_before_graft.residues) > margin:
        subpose_before_graft.delete_residue_range_slow(
            len(subpose_before_graft.residues) - margin,
            len(subpose_before_graft.residues),
        )
    if len(subpose_before_graft.residues) < margin:
        return link_poses(
            link_poses(fragment, subpose_after_graft), *ligands, rechain=True
        )

    n_term_half = link_poses(subpose_before_graft, fragment, rechain=True)
    num_chains = n_term_half.num_chains()

    print(f"segment lookup between chains {num_chains-1} and {num_chains}")
    segment_lookup_mover = run_direct_segment_lookup(
        n_term_half,
        from_chain=num_chains - 1,
        to_chain=num_chains,
        label=label,
    )

    status = segment_lookup_mover.get_last_move_status()
    if status != pyrosetta.rosetta.protocols.moves.mstype_from_name(
        "MS_SUCCESS"
    ):
        print("mover status: ", status)

        return

    if not len(subpose_after_graft.residues):
        print("aligned to last res, not linking subpose after graft")
        if ligands:
            return link_poses(n_term_half, *ligands, rechain=True)
        else:
            return n_term_half
    print("linking n_term half and c term half")
    grafted = link_poses(n_term_half, subpose_after_graft, rechain=False)
    if ligands:
        return link_poses(grafted, *ligands, rechain=True)
    return grafted


def graft(
    pose,
    fragment,
    pose_cut_position,
    *ligands,
    margin=3,
    allow_terminal=False,
    label="",
):
    """
    """
    subpose_after_graft = (
        pyrosetta.rosetta.protocols.grafting.return_region(
            pose.clone(), pose_cut_position + 1, len(pose.residues)
        )
        if len(pose.residues) > pose_cut_position
        else pyrosetta.rosetta.core.pose.Pose()
    )
    cut_pos = pose_cut_position - margin - 1
    if cut_pos < 0:
        if allow_terminal:
            return link_poses(
                link_poses(fragment, subpose_after_graft),
                *ligands,
                rechain=True,
            )
        else:
            return []
    subpose_before_graft = pyrosetta.rosetta.protocols.grafting.return_region(
        pose.clone(), 1, cut_pos
    )
    n_term_chains = [chain for chain in subpose_before_graft.split_by_chain()]
    link_portion = link_poses(n_term_chains[-1], fragment, rechain=True)
    # num_chains = n_term_half.num_chains()

    logger.debug(f"segment lookup:")

    segment_lookup_mover = run_direct_segment_lookup(
        link_portion, length=4, label=label
    )

    status = segment_lookup_mover.get_last_move_status()
    if status != pyrosetta.rosetta.protocols.moves.mstype_from_name(
        "MS_SUCCESS"
    ):
        print("mover status: ", status)

        return []
    grafted = []
    n_term_chains[-1] = link_portion
    if not len(subpose_after_graft.residues):
        print("aligned to last res, not linking subpose after graft")
        if ligands:
            grafted.append(
                link_poses(
                    *n_term_chains.split_by_chain(), *ligands, rechain=True
                )
            )
        else:
            grafted.append(
                link_poses(*n_term_chains.split_by_chain(), rechain=True)
            )
        for link_portion_more in iter(
            segment_lookup_mover.get_additional_output, None
        ):
            n_term_chains[-1] = link_portion_more
            if ligands:
                grafted.append(
                    link_poses(
                        *n_term_chains.split_by_chain(), *ligands, rechain=True
                    )
                )
            else:
                grafted.append(
                    link_poses(*n_term_chains.split_by_chain(), rechain=True)
                )
        return grafted
    print("linking n_term half and c term half")
    # n_term_chains = [chain for chain in n_term_half.split_by_chain()]
    c_term_chains = [chain for chain in subpose_after_graft.split_by_chain()]

    # n_term_final = n_term_chains[-1]
    c_term_first = c_term_chains[0]
    grafted_chain = link_poses(link_portion, c_term_first, rechain=False)

    all_chains = [*n_term_chains[:-1], grafted_chain, *c_term_chains[1:]]

    if ligands:
        grafted.append(link_poses(*all_chains, *ligands, rechain=True))
    else:
        grafted.append(link_poses(*all_chains, rechain=True))
    for link_portion_more in iter(
        segment_lookup_mover.get_additional_output, None
    ):
        c_term_first = c_term_chains[0]
        grafted_chain = link_poses(
            link_portion_more, c_term_first, rechain=False
        )

        all_chains = [*n_term_chains[:-1], grafted_chain, *c_term_chains[1:]]

        if ligands:
            grafted.append(link_poses(*all_chains, *ligands, rechain=True))
        else:
            grafted.append(link_poses(*all_chains, rechain=True))
    return grafted


def graft_fragment(
    pose,
    fragment,
    site,
    begin_only=True,
    allowed_depth=3,
    save_intermediate=True,
    get_additional_output=True,
    label="",
):
    """
    Super grafts the fragment based on positions matching secstruct to site

    begin only only anchors on the first residue (prepends fragment to graft site)
    allowed_depth gives the allowed number of residues to "chew into" site

    begin only as false not supported yet (appending fragment)
    """
    if site.dssp_type == "L":
        print("loops not supported")
    if site.dssp_type == "E":
        print("sheets not supported")
    if site.dssp_type == "H":
        fragment_sec_structs = parse_structure_from_dssp(fragment, "H")
        if not fragment_sec_structs:
            return []
        frag_starts, frag_ends = zip(
            *[
                (struct.start_pos, struct.end_pos)
                for struct in fragment_sec_structs
            ]
        )

        grafts = []

        for start in frag_starts:
            if start == 1:
                continue

            for i in range(allowed_depth + 1):
                if i > site.end_pos - site.start_pos:
                    break
                # work on a copy of the fragment
                insert = fragment.clone()
                graft_pos = site.start_pos + i
                super_xform = homog_super_transform_from_residues(
                    insert.residue(start), pose.residue(graft_pos)
                )
                rotation, translation = np_homog_to_rosetta_rotation_translation(
                    super_xform
                )
                rigid_body_move(
                    rotation,
                    translation,
                    insert,
                    TrueResidueSelector().apply(insert),
                    rosetta_vector(0, 0, 0),
                )
                # cut out the last chain of the fragment assume its phos
                phos = insert.split_by_chain()[insert.num_chains()]
                if len(insert.residues) > start:
                    insert.delete_residue_range_slow(
                        start + 1, len(insert.residues)
                    )

                # Make sure we're only operating on a single chain at a time
                insert = insert.split_by_chain()[chain_of(insert, start)]
                logger.debug(
                    f"""insert seq:
                {insert.annotated_sequence()}"""
                )
                logger.debug(
                    f"""pose seq:
                {pose.annotated_sequence()}"""
                )
                logger.debug(
                    f"""phos seq:
                {phos.annotated_sequence()}"""
                )
                logger.debug(f"site: {site}")
                logger.debug(f"pose graft site: {site}")

                grafted = graft(
                    pose.clone(),
                    insert,
                    graft_pos,
                    phos,
                    margin=i + 4,
                    label=label,
                )

                if grafted:
                    grafts.extend(grafted)

                    # grafted.dump_pdb(outname)

        if begin_only:
            if grafts:
                return grafts
            else:
                return []

        raise ValueError(
            "Dmitri hasn't implemented appending after anchor, only before"
        )
        for end in frag_ends:
            if end == len(fragment.residues):
                continue
            for i in range(4):
                if i > site.end_pos - site.start_pos:
                    break
                insert = fragment.clone()
                grafted = super_and_graft(
                    pose.clone(), insert, start, site.end_pos - i, label=label
                )
                grafts.append(grafted)
        if grafts:
            return grafts


def get_anchor_sites(
    pose, *dssp_types, anchor_end=True, res_to_end=0, struct_numbers=""
):
    """
    Gets the secstruct container for sites with appropriate anchor positions

    can require a certain number of residues from the end (to get an internal graft)
    """

    sec_structs = parse_structure_from_dssp(pose, *dssp_types)
    if struct_numbers:
        nums = [int(num) for num in struct_numbers.split(",")]
        sec_structs = [
            sec_structs[i] for i in range(len(sec_structs)) if i in nums
        ]
    if res_to_end:
        return sec_structs
    if anchor_end:
        return [
            site for site in sec_structs if site.start_pos - res_to_end > 1
        ]
    if not anchor_end:
        return [
            site
            for site in sec_structs
            if site.end_pos + res_to_end < len(pose.residues)
        ]


def graft_generator(
    pose,
    fragments,
    dssp_types="",
    anchor_end=True,
    save_intermediate=True,
    get_additional_output=True,
    struct_numbers="",
    label="",
):
    """
    Takes a pose, fragments, and secondary structure containers for the pose
    """
    # sec_structs = parse_structure_from_dssp(pose, *dssp_types)
    for fragment in fragments:
        for site in get_anchor_sites(
            pose,
            *dssp_types,
            anchor_end=anchor_end,
            struct_numbers=struct_numbers,
        ):
            for graft in graft_fragment(
                pose,
                fragment,
                site,
                save_intermediate=save_intermediate,
                get_additional_output=get_additional_output,
                label=label,
            ):

                yield graft


def accommodate_graft(pose, insertion_res_label, **kwargs):
    """
    modifes the pose to accommodate the grafted region with insertion_res_label

    Uses some defaults or options given in kwargs
    """

    # chew up loops
    loops = kwargs.get("trim_loops", False)
    chain_breaks = kwargs.get("chain_breaks", False)
    circ_perm = kwargs.get("circ_perm", False)
    trim_ss = kwargs.get("trim_ss", False)
    trim_frag = kwargs.get("trim_frag", False)
    relax = kwargs.get("rosetta_refine", False)
    # allow new chain breaks

    # allow circular permutations

    # Allowed ss trimming

    # trim fragment

    # min or relax fragment


@click.command()
@click.argument("pose_pdb", type=click.Path(exists=True))
@click.argument("fragment_store_path")
@click.option("-d", "--dssp-match-types", default="")
@click.option("-r", "--rosetta-flags-file")
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
    dssp_match_types="",
    rosetta_flags_file="",
    allowed_positions=False,
    save=True,
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
    # graft_sites = get_graft_sites(sec_structs)
    fragments = load_fragment_store_from_path(fragment_store_path)
    # remove preceding loop and graft
    pose = pyrosetta.pose_from_file(pose_pdb)

    grafts = []
    for i, graft in enumerate(
        graft_generator(
            pose,
            fragments,
            dssp_types=dssp_match_types,
            save_intermediate=save,
            get_additional_output=get_additional_output,
            struct_numbers=struct_numbers,
            label="anchored_graft",
        ),
        1,
    ):
        graft.dump_pdb(f"graft_{i}.pdb")


if __name__ == "__main__":
    main()
