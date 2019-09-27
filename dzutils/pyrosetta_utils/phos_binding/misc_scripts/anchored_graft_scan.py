import json
import pyrosetta
import click

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
    pose,
    fragment,
    insert_align,
    pose_align,
    *ligands,
    loop_before=False,
    margin=5,
):
    """
    """
    if loop_before:
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
                link_poses(fragment, subpose_after_graft),
                *ligands,
                rechain=True,
            )

        n_term_half = link_poses(subpose_before_graft, fragment, rechain=True)
        num_chains = n_term_half.num_chains()

        print(f"segment lookup between chains {num_chains-1} and {num_chains}")
        segment_lookup_mover = run_direct_segment_lookup(
            n_term_half, from_chain=num_chains - 1, to_chain=num_chains
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
    else:
        raise


def graft_fragment(pose, fragment, site, begin_only=True):
    """
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
            for i in range(4):
                if i > site.end_pos - site.start_pos:
                    break
                insert = fragment.clone()
                for j in range(4):
                    if len(insert.residues) > start:
                        insert.delete_residue_range_slow(
                            start + 1, len(insert.residues)
                        )
                    if j and j < len(insert.residues):
                        insert.delete_residue_range_slow(1, j)
                    if j == len(insert.residues) == 0:
                        continue

                    grafted = super_and_graft(
                        pose.clone(),
                        insert,
                        start - j,
                        site.start_pos + i,
                        fragment.clone().split_by_chain()[
                            fragment.num_chains()
                        ],
                        loop_before=True,
                    )
                    if grafted is not None:
                        grafts.append(grafted)

        # bb align frag start to site.start_pos

        # bb align frag end to site.end_pos
        if begin_only:
            if grafts:
                return grafts
            else:
                return []
        for end in frag_ends:
            if end == len(fragment.residues):
                continue
            for i in range(4):
                if i > site.end_pos - site.start_pos:
                    break
                insert = fragment.clone()
                grafted = super_and_graft(
                    pose.clone(),
                    insert,
                    start,
                    site.end_pos - i,
                    loop_before=False,
                )
                grafts.append(grafted)
        if grafts:
            return grafts


def get_anchor_sites(pose, *dssp_types, anchor_end=True, res_to_end=0):
    """
    Gets the secstruct container for sites with appropriate anchor positions

    can require a certain number of residues from the end (to get an internal graft)
    """
    sec_structs = parse_structure_from_dssp(pose, *dssp_types)
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


def graft_generator(pose, fragments, dssp_types="", anchor_end=True):
    """
    Takes a pose, fragments, and secondary structure containers for the pose
    """
    # sec_structs = parse_structure_from_dssp(pose, *dssp_types)
    for fragment in fragments:
        for site in get_anchor_sites(pose, *dssp_types, anchor_end=anchor_end):
            for graft in graft_fragment(pose, fragment, site):

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
def main(
    pose_pdb,
    fragment_store_path,
    # inv_rot_table_path,
    # inv_rot_dict_path,
    # log_dir,
    dssp_match_types="",
    rosetta_flags_file="",
    allowed_positions=False,
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
    fragments = load_fragment_store_from_path(fragment_store_path)[:5]
    # remove preceding loop and graft
    pose = pyrosetta.pose_from_file(pose_pdb)
    grafts = [*graft_generator(pose, fragments, dssp_types=dssp_match_types)]
    # report_grafts(grafts, log_dir)
    # delooped_grafts = [accommodate_graft(graft) for graft in grafts]
    # clash_checked_grafts = [*clash_checker(consensus)]
    # report_grafts(clash_checked_grafts, log_dir)
    # inv_rot_scanned = [
    #     *inv_rot_scanner(
    #         clash_checked_grafts, inv_rot_table_path, inv_rot_dict_path
    #     )
    # ]
    for i, graft in enumerate(grafts, 1):
        graft.dump_pdb(f"graft_{i}.pdb")
    # dump_results(inv_rot_scanned)


if __name__ == "__main__":
    main()
