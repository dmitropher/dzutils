import json
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
from dzutils.pyrosetta_utils.chain_utils import link_poses, chain_of

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


def graft(
    pose,
    fragment,
    n_term_position,
    c_term_position,
    *ligands,
    n_label="n_term_graft_site",
    c_label="c_term_graft_site",
    label="anchored_graft",
    allow_term=False,
):
    """
    Creates a pose where fragment is added between n_term_pos and c_term_pos

    Labels the term with "n_term_graft_site"/"c_term_graft_site"
    ligands are added as new chains at the end

    Does not close loops
    """
    logger.debug("test")
    pose_len = len(pose.residues)
    c_terminal = pose_len > c_term_position
    subpose_after_graft = (
        pyrosetta.rosetta.protocols.grafting.return_region(
            pose.clone(), c_term_position, pose_len
        )
        if c_terminal
        else pyrosetta.rosetta.core.pose.Pose()
    )

    c_term_label_index = 1 + n_term_position + len(fragment.residues)
    if n_term_position == 0:
        grafted = link_poses(
            link_poses(fragment, subpose_after_graft, rechain=False),
            *ligands,
            rechain=True,
        )
        pose.pdb_info().add_reslabel(c_term_label_index, c_label)
        return grafted if allow_term else None
    print("making n_term half")
    print(f"sequence: {pose.annotated_sequence()}")
    print(f"going from 1 to {n_term_position}")
    subpose_before_graft = pyrosetta.rosetta.protocols.grafting.return_region(
        pose.clone(), 1, n_term_position
    )

    n_term_chains = [chain for chain in subpose_before_graft.split_by_chain()]
    c_term_chains = [chain for chain in subpose_after_graft.split_by_chain()]
    grafted_chain = link_poses(
        *[
            chain
            for chain in [n_term_chains[-1], fragment, c_term_chains[0]]
            if len(chain.residues)
        ],
        rechain=False,
    )

    all_chains = [*n_term_chains[:-1], grafted_chain, *c_term_chains[1:]]
    all_chains = [chain for chain in all_chains if len(chain.residues)]

    if ligands:
        grafted = link_poses(*all_chains, *ligands, rechain=True)
    else:
        grafted = link_poses(*all_chains, rechain=True)

    c_term_label_index = 1 + n_term_position + len(fragment.residues)
    grafted.pdb_info().add_reslabel(c_term_label_index, c_label)
    grafted.pdb_info().add_reslabel(n_term_position, n_label)
    for i in range(n_term_position + 1, c_term_label_index):
        grafted.pdb_info().add_reslabel(i, label)

    return grafted


def compute_anchor_range(site, allowed_depth, use_start, nres):
    """
    """
    if use_start:
        # Use the start of the secstruct on target pose to anchor, cut forwards
        longest_run = min(
            [
                allowed_depth,
                nres - site.start_pos + 1,
                site.end_pos - site.start_pos + 1,
            ]
        )

        return range(0, longest_run)
    else:
        # use the end of the secstruct on target pose to anchor, cut backwards
        longest_run = min(
            [allowed_depth, site.end_pos - 1, site.end_pos - site.start_pos]
        )

        return range(0, -1 * longest_run, -1)


def graft_fragment(
    pose,
    fragment,
    site,
    use_start=True,
    allowed_depth=3,
    allowed_trim=6,
    dssp_type="H",
    label="anchored_graft",
):
    """
    Super grafts the fragment based on positions matching secstruct to site

    begin only only anchors on the first residue (prepends fragment to graft site)
    allowed_depth gives the allowed number of residues to "chew into" site

    begin only as false not supported yet (appending fragment)

    inserts with more than one chain treat all chains past the first as "ligands"
    grafting chain breaks not supported yet
    """

    fragment_sec_structs = parse_structure_from_dssp(fragment, dssp_type)

    anchors = [
        (struct.start_pos if use_start else struct.end_pos)
        for struct in fragment_sec_structs
    ]
    site_pos = site.start_pos if use_start else site.end_pos
    grafts = []
    for anchor in anchors:
        logger.debug(f"anchor: {anchor}")
        working_frag = fragment.clone()
        # ligands = list(working_frag.split_by_chain())[
        #     1 : working_frag.num_chains()
        # ]

        print(f"use_start {use_start}")
        if use_start:
            if len(working_frag.residues) > anchor:
                working_frag.delete_residue_range_slow(
                    anchor + 1, len(working_frag.residues)
                )
        else:
            if anchor > 1:
                working_frag.delete_residue_range_slow(1, anchor - 1)
                anchor = 1

        working_insert = working_frag.split_by_chain()[1]
        i_range = compute_anchor_range(
            site, allowed_depth, use_start, len(pose.residues)
        )
        print(i_range, list(i_range))
        start_val = (-1) ** (not use_start)
        for i in i_range:
            insert = working_insert.clone()
            ligands = list(fragment.clone().split_by_chain())[
                1 : fragment.num_chains()
            ]
            # work on a copy of the fragment
            super_target = site_pos + i
            print(f"graft_pos: {super_target}, site_pos {site_pos}")
            print(f"insert sequence: {insert.annotated_sequence()}")
            print(f"pose_sequence: {pose.annotated_sequence()}")
            try:
                super_xform = homog_super_transform_from_residues(
                    insert.residue(anchor), pose.residue(super_target)
                )
                rotation, translation = np_homog_to_rosetta_rotation_translation(
                    super_xform
                )
                for p in [insert, *ligands]:
                    rigid_body_move(
                        rotation,
                        translation,
                        p,
                        TrueResidueSelector().apply(p),
                        rosetta_vector(0, 0, 0),
                    )
            except RuntimeError as e:
                print(e)
                print(
                    "anchor",
                    anchor,
                    len(insert.residues),
                    "graft_pos",
                    super_target,
                    len(pose.residues),
                )

                continue

            logger.debug("insert seq:")
            logger.debug(f"{insert.annotated_sequence()}")
            logger.debug("pose seq:")
            logger.debug(f"{pose.annotated_sequence()}")
            logger.debug("ligands seq:")
            logger.debug(
                f"{'   '.join([lig.annotated_sequence() for lig in ligands ])}"
            )
            logger.debug(f"site: {site}")
            logger.debug(f"pose graft site: {site}")
            cut_border = super_target + start_val
            print(f"site_pos { site_pos}")
            j_range = (
                min(allowed_trim, site_pos - 1)
                if use_start
                else min(allowed_trim, len(pose.residues) - site_pos)
            )

            for j in range(0, j_range * start_val, start_val):
                other_site = max(0, site_pos - j - start_val)
                n_anchor = min(other_site, cut_border)
                c_anchor = max(other_site, cut_border)
                logger.debug(f"n_anchor {n_anchor}")
                logger.debug(f"c_anchor {c_anchor}")
                grafted = graft(
                    pose.clone(),
                    insert.clone(),
                    n_anchor,
                    c_anchor,
                    *ligands,
                    label=label,
                )

                if not grafted is None:
                    grafts.append(grafted)

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
    allowed_depth=3,
    struct_numbers="",
    label="anchored_graft",
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
                use_start=anchor_end,
                label=label,
                allowed_depth=allowed_depth,
            ):

                yield graft


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

    # grafts = []
    for i, graft in enumerate(
        graft_generator(
            pose,
            fragments,
            dssp_types=dssp_match_types,
            label="anchored_graft",
        ),
        1,
    ):
        graft.dump_pdb(f"graft_{i}.pdb")


if __name__ == "__main__":
    main()
