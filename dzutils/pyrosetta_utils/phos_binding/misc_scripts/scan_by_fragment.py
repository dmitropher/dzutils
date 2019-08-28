# Attempts to close from each chain to each other chain
# For each chain pair, does serial deletions from the ends to increase
# the number of samples

import os
import sys
import logging

import pyrosetta


from dzutils.sutils import read_flag_file

from dzutils.pyrosetta_utils.phos_binding.misc_scripts.remove_loops_rechain import (
    ss_to_chains,
)
from dzutils.pyrosetta_utils.phos_binding.misc_scripts.exhaustive_single_loop_insertion import (
    exhaustive_single_loop_insertion,
)

from dzutils.pyrosetta_utils.phos_binding.misc_scripts.double_scan_pbinders import (
    scan_for_ploop_graft,
    scan_for_inv_rot,
)

from dzutils.pyrosetta_utils.geometry.superposition_utilities import (
    super_by_residues,
)

# from dzutils.pyrosetta_utils.chain_utils import insert_pose


def super_and_insert_pose(start, end, pose, insertion, insertion_start):
    moved_insertion = super_by_residues(
        insertion, pose, insertion_start, start
    )
    newp = pyrosetta.rosetta.core.pose.Pose()
    newp.detached_copy(pose)
    # pyrosetta.rosetta.core.pose.remove_variant_type_from_pose_residue(
    #     moved_insertion,
    #     pyrosetta.rosetta.core.chemical.LOWER_TERMINUS_VARIANT,
    #     1,
    # )
    # pyrosetta.rosetta.core.pose.remove_variant_type_from_pose_residue(
    #     moved_insertion,
    #     pyrosetta.rosetta.core.chemical.UPPER_TERMINUS_VARIANT,
    #     len(moved_insertion.residues),
    # )
    pyrosetta.rosetta.protocols.grafting.delete_region(newp, start, end)  # - 1
    # newp.dump_pdb("/home/dzorine/temp/trimmed_before_insert.pdb")
    print(start, end)
    # print(moved_insertion)
    # print(newp)
    pyrosetta.rosetta.core.pose.remove_variant_type_from_pose_residue(
        newp, pyrosetta.rosetta.core.chemical.UPPER_TERMINUS_VARIANT, start - 1
    )
    out_pose = pyrosetta.rosetta.protocols.grafting.insert_pose_into_pose(
        newp, moved_insertion, start - 1, start, False
    )
    # out_pose = insert_pose(newp, moved_insertion, start - 1)
    return out_pose


def replace_res_from_pose(pose, replacement, index, replacement_index):
    replacement_copy = replacement.clone()
    # super_by_residues(replacement_copy,pose,replacement_index,index)
    pyrosetta.rosetta.core.pose.remove_variant_type_from_pose_residue(
        replacement_copy,
        pyrosetta.rosetta.core.chemical.UPPER_TERMINUS_VARIANT,
        replacement_index,
    )
    pyrosetta.rosetta.core.pose.remove_variant_type_from_pose_residue(
        replacement_copy,
        pyrosetta.rosetta.core.chemical.LOWER_TERMINUS_VARIANT,
        replacement_index,
    )
    pose.replace_residue(
        index, replacement_copy.residue(replacement_index), True
    )
    # pyrosetta.rosetta.core.pose.correctly_add_cutpoint_variants(pose)
    for num in range(1, pose.num_chains() + 1):
        if index == pose.chain_begin(num):
            pyrosetta.rosetta.core.pose.add_variant_type_to_pose_residue(
                pose,
                pyrosetta.rosetta.core.chemical.LOWER_TERMINUS_VARIANT,
                index,
            )
            pose.conformation().chains_from_termini()
        if index == pose.chain_end(num):
            pyrosetta.rosetta.core.pose.add_variant_type_to_pose_residue(
                pose,
                pyrosetta.rosetta.core.chemical.UPPER_TERMINUS_VARIANT,
                index,
            )
            pose.conformation().chains_from_termini()
    return pose


def graft_and_dump_pdb(
    pose, insertion, start, end, insertion_start, dump_path
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


def process_secondary_results(
    secondary_results_table,
    pose,
    secondary_hit_out_dir,
    source_file_dir,
    table_label="",
    pdb_label="",
):
    """
    """
    print(secondary_results_table)
    if not os.path.isdir(secondary_hit_out_dir):
        os.makedirs(secondary_hit_out_dir)
    secondary_write_path = f"""{secondary_hit_out_dir
                               }/secondary_double_ploop_hits{
                               table_label
                               }.json"""

    # assert not os.path.exists(
    #     secondary_write_path
    # ), f"secondary hit data would overwrite file at: {secondary_write_path}"
    # source_file_dir = f"{secondary_hit_out_dir}/source_files/"
    # if not os.path.isdir(source_file_dir):
    #     os.makedirs(source_file_dir)
    secondary_grafted_loops_dir = (
        f"{secondary_hit_out_dir}/grafted_loops_with_p_res/"
    )
    if not os.path.isdir(secondary_grafted_loops_dir):
        os.makedirs(secondary_grafted_loops_dir)
    if not os.path.isdir(source_file_dir):
        os.makedirs(source_file_dir)
    # Add a grafted file entry with phos rotamer to each row, dump the pdb
    print("grafting and adding residue")
    secondary_results_table[
        "graft_and_pres_path"
    ] = secondary_results_table.apply(
        lambda row: graft_and_dump_pdb(
            replace_res_from_pose(
                pose.clone(),
                pyrosetta.pose_from_file(row["inv_rot_file"]),
                int(row["p_res_target"]),
                1,
            ),
            pyrosetta.pose_from_file(row["loop_file"]).split_by_chain()[1],
            row["start_res"],
            row["end_res"],
            1,
            f"""{secondary_grafted_loops_dir
                }/{
                row["name"].split("/")[-1].split(".pdb")[0]
                }{pdb_label}_ploop_index_{
                row["loop_index"]
                }_pres_{
                row["p_res_target"]
                }.pdb""",
        ),
        axis=1,
    )
    secondary_results_table.apply(
        lambda row: pyrosetta.pose_from_file(row["inv_rot_file"]).dump_pdb(
            f'{source_file_dir}/{row["inv_rot_file"].split("/")[-1]}'
        ),
        axis=1,
    )
    secondary_results_table.to_json(secondary_write_path)


def main():
    # flagsFile = "/home/dzorine/phos_binding/pilot_runs/loop_grafting/initial_testing/misc_files/p_ligand_quiet.flags"
    flagsFile = "/home/dzorine/phos_binding/pilot_runs/loop_grafting/initial_testing/misc_files/p_ligand_quiet.flags"
    flags = read_flag_file(flagsFile)
    flags_str = " ".join(flags.replace("\n", " ").split())
    print(flags_str)
    pyrosetta.init(flags_str)

    pose = pyrosetta.pose_from_file(sys.argv[1])
    in_dir = sys.argv[2]
    out_dir = sys.argv[3]
    primary_hit_dirname = sys.argv[4]
    secondary_hit_dirname = sys.argv[5]
    overwrite = True

    # splits by helices
    ss_only_chains = ss_to_chains(pose, "H", "E")
    pose_name = pose.pdb_info().name().split(".pdb")[0].split("/")[-1]
    suffix = "_secstruct_only.pdb"
    ss_only_name = f"{in_dir}/source_files/{pose_name}{suffix}"

    # mind the hardcode
    ss_only_chains.dump_pdb(ss_only_name)
    ss_only_chains = pyrosetta.pose_from_file(ss_only_name)
    print(f"generating poses for {pose_name}")
    # poses = [
    # pose for pose in exhaustive_single_loop_insertion(ss_only_chains, 5)
    # ]
    for i, pose in enumerate(
        exhaustive_single_loop_insertion(ss_only_chains, 5), 1
    ):
        # for i, pose in enumerate(poses, 1):
        print("beginning scan of closed loop pose")
        # key_type, value_type = np.dtype("i8"), np.dtype("i8")
        data_name = "exp_no_spin_1ang_3_contact_ploop_set_3_v1"
        data_path = "/home/dzorine/phos_binding/pilot_runs/loop_grafting/fragment_tables/ploops_expanded_set_1"
        table_dir_path = f"{data_path}/tables"
        dict_dir_path = f"{data_path}/dicts"
        dict_path = f"{dict_dir_path}/{data_name}.bin"
        table_path = f"{table_dir_path}/{data_name}.json"
        primary_results_table = scan_for_ploop_graft(
            pose, table_path, dict_path
        )
        if primary_results_table is None:
            logging.info("no primary hits found for this pose")
            continue

        pdb_name = pose.pdb_info().name().split("/")[-1].split(".pdb")[0]
        # results.name = "double_ploop_hits"

        # prepare all the output directories
        primary_hit_out_dir = f"{out_dir}/{primary_hit_dirname}"
        os.makedirs(primary_hit_out_dir, exist_ok=True)
        primary_write_path = f"{primary_hit_out_dir}/primary_double_ploop_hits_{pdb_name}_loop_{i}.json"
        if not overwrite:
            assert not os.path.exists(
                primary_write_path
            ), f"primary hit data would overwrite file at: {primary_write_path}"
        source_file_dir = f"{primary_hit_out_dir}/source_files/"
        os.makedirs(source_file_dir, exist_ok=True)
        grafted_loops_dir = f"{primary_hit_out_dir}/grafted_loops/"
        os.makedirs(grafted_loops_dir, exist_ok=True)
        # Add a grafted file entry to each row, dump the pdb to that
        # CODECODECODE
        print("primary hits found: ")
        # print(primary_results_table)
        primary_results_table[
            "loop_insertion_path"
        ] = primary_results_table.apply(
            lambda row: graft_and_dump_pdb(
                pose.clone(),
                pyrosetta.pose_from_file(row["loop_file"]).split_by_chain()[1],
                row["start_res"],
                row["end_res"],
                1,
                f'{grafted_loops_dir}/{row["name"].split("/")[-1].split(".pdb")[0]}_reloop_{i}_ploop_index_{row["loop_index"]}.pdb',
            ),
            axis=1,
        )
        print("grafted loops dumped")
        primary_results_table.apply(
            lambda row: pyrosetta.pose_from_file(row["loop_file"]).dump_pdb(
                f'{source_file_dir}/{row["loop_file"].split("/")[-1]}'
            ),
            axis=1,
        )
        print("loop file dumped")
        primary_results_table.to_json(primary_write_path)
        primary_hit_pdb_path = f"{source_file_dir}/{pdb_name}_reloop_{i}.pdb"
        print(f"dumping pdb at: {primary_hit_pdb_path}")
        pose.dump_pdb(primary_hit_pdb_path)

        # BEGIN SECONDARY SCAN

        inv_rot_data_path = "/home/dzorine/projects/phos_binding/pilot_runs/loop_grafting/fragment_tables/inverse_rotamer"
        inv_rot_data_name = "inverse_ptr_exchi7_rotamers_1ang_v1"
        inv_rot_dict_path = (
            f"{inv_rot_data_path}/dicts/{inv_rot_data_name}.bin"
        )
        inv_rot_table_path = (
            f"{inv_rot_data_path}/tables/{inv_rot_data_name}.json"
        )

        secondary_results_table = scan_for_inv_rot(
            pose, primary_results_table, inv_rot_table_path, inv_rot_dict_path
        )
        if secondary_results_table is None:
            print("no secondary hits found for this pose")
            continue

        # prepare all the output directories
        secondary_hit_out_dir = f"{out_dir}/{secondary_hit_dirname}"
        sec_source_file_dir = f"{secondary_hit_out_dir}/source_files/"
        process_secondary_results(
            secondary_results_table,
            pose,
            secondary_hit_out_dir,
            sec_source_file_dir,
            table_label=f"_{pdb_name}_loop_{i}",
            pdb_label=f"_reloop_{i}",
        )


if __name__ == "__main__":
    main()
