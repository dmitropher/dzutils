# Attempts to close from each chain to each other chain
# For each chain pair, does serial deletions from the ends to increase
# the number of samples

import os
import sys

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


def super_and_insert_pose(start, end, pose, insertion, insertion_start):
    moved_insertion = super_by_residues(
        insertion, pose, insertion_start, start
    )
    newp = pose.clone()
    newp.delete_residue_range_slow(start, end)
    print(start)
    print(moved_insertion)
    print(newp)
    return pyrosetta.rosetta.protocols.grafting.insert_pose_into_pose(
        newp, moved_insertion, start - 1
    )


def replace_res_from_pose(pose, replacement, index, replacement_index):
    pose.replace_residue(index, replacement.residue(replacement_index), False)
    return pose


def graft_and_dump_pdb(
    pose, insertion, start, end, insertion_start, dump_path
):
    """
    """
    if not super_and_insert_pose(
        start, end, pose, insertion, insertion_start
    ).dump_pdb(dump_path):
        raise RuntimeError(f"Unable to dump pdb at specfied path: {dump_path}")
    return dump_path


def main():
    flagsFile = "/home/dzorine/phos_binding/pilot_runs/loop_grafting/initial_testing/misc_files/cluster_altered.flags"
    flags = read_flag_file(flagsFile)
    flags_str = " ".join(flags.replace("\n", " ").split())
    pyrosetta.init(flags_str)

    pose = pyrosetta.pose_from_file(sys.argv[1])
    out_dir = sys.argv[2]
    primary_hit_dirname = sys.argv[3]
    secondary_hit_dirname = sys.argv[4]

    # splits by helices
    ss_only_chains = ss_to_chains(pose, "H", "E")
    pose_name = pose.pdb_info().name().split(".pdb")[0].split("/")[-1]
    suffix = "_secstruct_only.pdb"
    ss_only_name = f"{out_dir}/source_files/{pose_name}{suffix}"

    # mind the hardcode
    ss_only_chains.dump_pdb(ss_only_name)
    ss_only_chains = pyrosetta.pose_from_file(ss_only_name)
    for i, pose in enumerate(
        exhaustive_single_loop_insertion(ss_only_chains, out_dir, 5), 1
    ):
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
        if not primary_results_table:
            print("no primary hits found for this pose")
            continue
        pdb_name = pose.pdb_info().name().split("/")[-1].split(".pdb")[0]
        # results.name = "double_ploop_hits"

        primary_hit_out_dir = f"{out_dir}/{primary_hit_dirname}"
        if not os.path.isdir(primary_hit_out_dir):
            os.makedirs(primary_hit_out_dir)
        primary_write_path = f"{primary_hit_out_dir}/primary_double_ploop_hits_{pdb_name}_loop_{i}.json"
        assert not os.path.exists(
            primary_write_path
        ), f"primary hit data would overwrite file at: {primary_write_path}"
        source_file_dir = f"{primary_hit_out_dir}/source_files/"
        if not os.path.isdir(source_file_dir):
            os.makedirs(source_file_dir)
        grafted_loops_dir = f"{primary_hit_out_dir}/grafted_loops/"
        if not os.path.isdir(grafted_loops_dir):
            os.makedirs(grafted_loops_dir)
        # Add a grafted file entry to each row, dump the pdb to that
        # CODECODECODE

        primary_results_table[
            "loop_insertion_path"
        ] = primary_results_table.apply(
            lambda row: graft_and_dump_pdb(
                pyrosetta.pose_from_file(row["scaffold_path"]),
                pyrosetta.pose_from_file(row["loop_file"]),
                row["start_res"],
                row["end_res"],
                1,
                f'{grafted_loops_dir}/{row["scaffold_path"].split("/")[-1].split(".pdb")[0]}_reloop_{i}_ploop_index_{row["loop_index"]}.pdb',
            ),
            axis=1,
        )
        primary_results_table.apply(
            lambda row: pyrosetta.pose_from_file(
                row["scaffold_path"]
            ).dump_pdb(
                f'{grafted_loops_dir}/{row["scaffold_path"].split("/")[-1].split(".pdb")[0]}_reloop_{i}_ploop_index_{row["loop_index"]}.pdb'
            ),
            axis=1,
        )
        primary_results_table.to_json(primary_write_path)
        primary_hit_pdb_path = f"{source_file_dir}/{pdb_name}_reloop_{i}.pdb"
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
        if not secondary_results_table:
            print("no secondary hits found for this pose")
            continue
        secondary_write_path = f"{secondary_hit_out_dir}/secondary_double_ploop_hits_{pdb_name}_reloop_{i}.json"

        # Add a grafted file entry with phos rotamer to each row, dump the pdb

        secondary_results_table[
            "graft_and_pres_path"
        ] = secondary_results_table.apply(
            lambda row: graft_and_dump_pdb(
                replace_res_from_pose(
                    pyrosetta.pose_from_file(row["scaffold_path"]),
                    pyrosetta.pose_from_file(row["inv_rot_file"]),
                    row["p_res_target"],
                    1,
                ),
                pyrosetta.pose_from_file(row["loop_file"]),
                row["start_res"],
                row["end_res"],
                1,
                f'{row["scaffold_path"].split("/")[-1].split(".pdb")[0]}_reloop_{i}_ploop_index_{row["loop_index"]}_pres_{row["p_res_target"]}.pdb',
            ),
            axis=1,
        )
        secondary_results_table.to_json(secondary_write_path)


if __name__ == "__main__":
    main()
