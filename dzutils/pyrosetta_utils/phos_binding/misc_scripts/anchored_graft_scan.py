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


def graft(
    pose,
    fragment,
    n_term_position,
    c_term_position,
    *ligands,
    n_label="n_term_graft_site",
    c_label="c_term_graft_site",
):
    """
    Creates a pose where fragment is added between n_term_pos and c_term_pos

    Labels the term with "n_term_graft_site"/"c_term_graft_site"
    ligands are added as new chains at the end

    Does not close loops
    """
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

    subpose_before_graft = pyrosetta.rosetta.protocols.grafting.return_region(
        pose.clone(), 1, n_term_position
    )

    n_term_chains = [chain for chain in subpose_before_graft.split_by_chain()]
    c_term_chains = [chain for chain in subpose_after_graft.split_by_chain()]
    grafted_chain = link_poses(
        n_term_chains[-1], fragment, c_term_chains[0], rechain=False
    )

    all_chains = [*n_term_chains[:-1], grafted_chain, *c_term_chains[1:]]

    if ligands:
        grafted = link_poses(*all_chains, *ligands, rechain=True)
    else:
        grafted = link_poses(*all_chains, rechain=True)

    c_term_label_index = 1 + n_term_position + len(fragment.residues)
    pose.pdb_info().add_reslabel(c_term_label_index, c_label)
    pose.pdb_info().add_reslabel(n_term_position, n_label)

    return grafted


def graft_fragment(
    pose, fragment, site, use_start=True, allowed_depth=3, allowed_trim=4
):
    """
    Super grafts the fragment based on positions matching secstruct to site

    begin only only anchors on the first residue (prepends fragment to graft site)
    allowed_depth gives the allowed number of residues to "chew into" site

    begin only as false not supported yet (appending fragment)

    inserts with more than two chains (aa and phos) are not supported yet
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
        for anchor in frag_starts:
            if anchor == 1:
                continue

            for i in range(allowed_depth + 1):
                if i > site.end_pos - site.start_pos:
                    break
                # work on a copy of the fragment
                insert = fragment.clone()
                graft_pos = (
                    site.start_pos + i if use_start else site.end_pos - i
                )
                super_xform = homog_super_transform_from_residues(
                    insert.residue(anchor), pose.residue(graft_pos)
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
                if use_start:
                    if len(insert.residues) > anchor:
                        insert.delete_residue_range_slow(
                            anchor + 1, len(insert.residues)
                        )
                else:
                    if anchor > 1:
                        insert.delete_residue_range_slow(1, anchor - 1)
                # Make sure we're only operating on a single chain at a time
                insert = insert.split_by_chain()[chain_of(insert, anchor)]
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
                for j in range(allowed_trim + 1):
                    n_anchor = (
                        graft_pos - (i + j) - 1 if use_start else graft_pos
                    )
                    c_anchor = (
                        graft_pos if use_start else graft_pos + (i + j) + 1
                    )
                    if n_anchor < 1:
                        break
                    grafted = graft(
                        pose.clone(), insert, n_anchor, c_anchor, phos
                    )

                    if grafted:
                        grafts.append(grafted)

                    # grafted.dump_pdb(outname)

        if grafts:
            return grafts
        else:
            return []


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
    loop_close=False,
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
            save_intermediate=save,
            get_additional_output=get_additional_output,
            struct_numbers=struct_numbers,
            label="anchored_graft",
            loop_close=True,
        ),
        1,
    ):
        graft.dump_pdb(f"graft_{i}.pdb")


if __name__ == "__main__":
    main()
