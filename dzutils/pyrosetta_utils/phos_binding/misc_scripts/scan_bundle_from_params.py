import sys, json, os

import click
import numpy as np

import pyrosetta

from dzutils.pyrosetta_utils.phos_binding.parametric import (
    helix_bundle_maker_wrapper,
    PyBundleGridSampler,
)
import pandas as pd
from dzutils.pyrosetta_utils.phos_binding.misc_scripts.double_scan_pbinders import (
    scan_for_n_term_helical_grafts,
    scan_for_inv_rot,
    load_table_and_dict,
)

from dzutils.pyrosetta_utils.phos_binding import (
    super_and_insert_pose,
    replace_res_from_pose,
)
from dzutils.pyrosetta_utils.chain_utils import chain_of

from dzutils.pyrosetta_utils import run_pyrosetta_with_flags


def trim_frag_to_end(pose, position, chain=1):
    frag = pose.clone()
    # data["frag_end"] + 1,
    end = len(frag.split_by_chain()[chain].residues)
    if position <= end:
        frag.delete_residue_range_slow(
            position, len(frag.split_by_chain()[chain].residues)
        )
    return frag


def trim_to_graft_start(pose, resnum):
    """
    trims the pose from the chain start of resnum to the residue before

    If resnum is at the start of a chain, does nothing
    """
    chain_num = chain_of(pose, resnum)
    chain_begin = pose.chain_begin(chain_num)
    if resnum != chain_begin:
        pose.delete_residue_range_slow(chain_begin, resnum - 1)
    return pose


def process_secondary_results(
    secondary_results_table,
    pose,
    out_dir,
    target_res_label="p_res_target",
    fragment_file_label="frag_filename",
    frag_index_label="frag_index",
    inv_rot_file_label="inv_rot_file",
    frag_end_label="frag_end",
    frag_start_label="frag_start",
    table_label="",
    pdb_label="",
):
    """
    """
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    secondary_write_path = f"""{out_dir
                               }/secondary_double_ploop_hits{
                               table_label
                               }.json"""

    grafted_file_dir = f"{out_dir}/grafted_helices"

    if not os.path.isdir(grafted_file_dir):
        os.makedirs(grafted_file_dir)
    sec_source_file_dir = f"{out_dir}/source_files"
    if not os.path.isdir(sec_source_file_dir):
        os.makedirs(sec_source_file_dir)
    # Add a grafted file entry with phos rotamer to each row, dump the pdb
    # Loop is unclosed!
    secondary_results_table[
        "graft_and_pres_path"
    ] = secondary_results_table.apply(
        lambda row: f"""{grafted_file_dir
                }/{
                row["name"].split("/")[-1].split(".pdb")[0]
                }{pdb_label}_ploop_index_{
                row[frag_index_label]
                }_pres_{
                row[target_res_label]
                }.pdb""",
        axis=1,
    )
    secondary_results_table.apply(
        lambda row: trim_to_graft_start(
            super_and_insert_pose(
                int(row["start_res"]),
                int(row["end_res"]),
                replace_res_from_pose(
                    pose.clone(),
                    pyrosetta.pose_from_file(row[inv_rot_file_label]),
                    int(row[target_res_label]),
                    1,
                ),
                trim_frag_to_end(
                    pyrosetta.pose_from_file(row[fragment_file_label]),
                    int(row[frag_end_label] + 1),
                ),
                row[frag_start_label],
                insert_chain=1,
            ),
            int(row["start_res"]),
        ).dump_pdb(row["graft_and_pres_path"]),
        axis=1,
    )
    secondary_results_table.apply(
        lambda row: pyrosetta.pose_from_file(row[inv_rot_file_label]).dump_pdb(
            f'{sec_source_file_dir}/{row[inv_rot_file_label].split("/")[-1]}'
        ),
        axis=1,
    )


def param_to_pose(param_dict):
    """
    """
    params_for_helices = [
        param_dict[str(i)] for i in range(1, param_dict["num_helices"] + 1)
    ]

    # Use the wrapper to create the MakeBundle Mover
    bundle_maker = helix_bundle_maker_wrapper(
        param_dict["helix_length"], *params_for_helices, degrees=True
    )

    # Empty pose to helical bundle!
    pose = pyrosetta.rosetta.core.pose.Pose()
    bundle_maker.apply(pose)
    return pose


def num_params(param_range_dict):
    """
    Takes a PyBundleGridSampler params dict and computes total combinations

    dict must have either range type params or binary and must specify num helix
    """
    return np.product(
        [
            (
                dict_["steps"]
                if dict_["type"] == "range"
                else ((2 if dict_["allow_others"] else 1))
            )
            for val in range(1, param_range_dict["num_helices"] + 1)
            for dict_ in param_range_dict[val]
        ]
    )


def subsample_grid(num_chunks, index, **grid_params):
    """
    Takes a set of grid params and slices off the given chunk
    """


@click.command()
@click.option("-o", "--out-dir")
@click.option("-p", "--param-json")
@click.option("-i", "--chunk-index")
@click.option("-c", "--num-chunks")
@click.option(
    "-f",
    "--flags-file",
    default="/home/dzorine/phos_binding/run_files/p_ligand_quiet.flags",
)
def main(out_dir, param_json, flags_file, chunk_index, num_chunks):

    table_path = "/home/dzorine/phos_binding/pilot_runs/loop_grafting/fragment_tables/ploop_fragment_set/tables/1ang_3_contact_ploop_set_4_v3.json"
    dict_path = "/home/dzorine/phos_binding/pilot_runs/loop_grafting/fragment_tables/ploop_fragment_set/dicts/1ang_3_contact_ploop_set_4_v3.bin"
    rot_table_path = "/home/dzorine/phos_binding/pilot_runs/loop_grafting/fragment_tables/inverse_rotamer/tables/inverse_ptr_exchi7_rotamers_1ang_v1.json"
    rot_dict_path = "/home/dzorine/phos_binding/pilot_runs/loop_grafting/fragment_tables/inverse_rotamer/dicts/inverse_ptr_exchi7_rotamers_1ang_v1.bin"
    frag_table, frag_dict = load_table_and_dict(
        table_path, dict_path, np.dtype("i8"), np.dtype("i8")
    )
    rot_table, rot_dict = load_table_and_dict(
        rot_table_path, rot_dict_path, np.dtype("i8"), np.dtype("i8")
    )
    run_pyrosetta_with_flags(flags_file)
    grid_dict = json.loads(param_json)
    sub_grid_dicts = subsample_grid(num_chunks, chunk_index, grid_dict)
    name_hash = hash(param_json)
    for param_dict in sub_grid_dicts:
        # param_dict = dict(param_json)

        pose = param_to_pose(param_dict)
        inf = pyrosetta.rosetta.core.pose.PDBInfo(pose)
        inf.name(f"param_bundle_{name_hash}.pdb")
        pose.pdb_info(inf)

        if not len(pose.residues):
            continue
        results = scan_for_n_term_helical_grafts(pose, table_path, dict_path)
        results["allowed_res"] = results.apply(
            lambda x: [*range(1, len(pose.residues) + 1)], axis=1
        )
        results["param_json"] = results.apply(lambda x: param_json, axis=1)
        results = results.rename(
            index=str,
            columns={"frag_feature_to_start": "loop_func_to_bb_start"},
        )
        sec_results = scan_for_inv_rot(
            pose, results, rot_table_path, rot_dict_path
        )
        if not sec_results is None:
            process_secondary_results(sec_results, pose, out_dir)


if __name__ == "__main__":
    main()
