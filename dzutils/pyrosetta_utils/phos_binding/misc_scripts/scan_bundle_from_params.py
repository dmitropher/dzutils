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

from dzutils import pythonify

from dzutils.pyrosetta_utils.phos_binding import (
    super_and_insert_pose,
    replace_res_from_pose,
)
from dzutils.pyrosetta_utils.chain_utils import chain_of

from dzutils.pyrosetta_utils import run_pyrosetta_with_flags

from dzutils.func_utils import index_to_parameters


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


def ptr_from_chis(*chis):
    # get chemical manager:
    chemical_manager = (
        pyrosetta.rosetta.core.chemical.ChemicalManager.get_instance()
    )
    rts = chemical_manager.residue_type_set(
        pyrosetta.rosetta.core.chemical.TypeSetMode.FULL_ATOM_t
    )
    # Get phospho-TYR
    res_type = rts.get_residue_type_with_variant_added(
        rts.name_map("TYR"),
        pyrosetta.rosetta.core.chemical.VariantType.PHOSPHORYLATION,
    )
    ptr = pyrosetta.rosetta.core.conformation.Residue(res_type, False)
    pose = pyrosetta.rosetta.core.pose.Pose()
    pose.append_residue_by_bond(ptr, False)
    for i, chi in enumerate(chis, 1):
        pose.set_chi(i, 1, chi)
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
                    ptr_from_chis(*row["inv_rot_chis"]),
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
        lambda row: ptr_from_chis(*row["inv_rot_chis"]).dump_pdb(
            f'{sec_source_file_dir}/ptr_{"_".join(str(chi) for chi in row["inv_rot_chis"])}.pdb'
        ),
        axis=1,
    )


def param_to_pose(param_dict):
    """
    """
    params_for_helices = [
        param_dict[i] for i in range(1, param_dict["num_helices"] + 1)
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
    return int(
        np.product(
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
    )


def subsample_grid(num_chunks, index, param_range_dict):
    """
    Takes a set of grid params and slices off the given chunk
    """
    big = num_params(param_range_dict)
    binnable = big - (big % (num_chunks - 1))
    bin_size = int(binnable / (num_chunks - 1))
    indices = (
        list(
            range(
                bin_size * index + 1, min(bin_size * (index + 1) + 1, big + 1)
            )
        )
        if bin_size > 1
        else [index]
    )
    param_tuples = [
        (
            (dict_["start"], dict_["stop"], dict_["steps"])
            if dict_["type"] == "range"
            else (
                int(dict_["base_value"]),
                (
                    int(not dict_["base_value"])
                    if dict_["allow_others"]
                    else int(dict_["base_value"])
                ),
                (2 if dict_["allow_others"] else 1),
            )
        )
        for val in range(1, param_range_dict["num_helices"] + 1)
        for dict_ in param_range_dict[val]
    ]
    param_names = [
        (val, dict_["name"], dict_["type"])
        for val in range(1, param_range_dict["num_helices"] + 1)
        for dict_ in param_range_dict[val]
    ]
    grid_points = (
        zip(param_names, index_to_parameters(i, *param_tuples))
        for i in indices
    )
    for grid_point in grid_points:
        helix_params = {
            "num_helices": param_range_dict["num_helices"],
            "helix_length": param_range_dict["helix_length"],
        }
        for (helix, name, param_type), value in grid_point:
            if helix in helix_params:

                helix_params[helix][name] = float(value)
            else:
                helix_params[helix] = {name: float(value)}
        yield helix_params


@click.command()
@click.option("-o", "--out-dir")
@click.option("-p", "--param-json-path")
@click.option("-i", "--chunk-index")
@click.option("-c", "--num-chunks")
@click.option(
    "-f",
    "--flags-file",
    default="/home/dzorine/phos_binding/run_files/p_ligand_quiet.flags",
)
def main(out_dir, param_json_path, flags_file, chunk_index, num_chunks):

    table_path = "/home/dzorine/phos_binding/pilot_runs/loop_grafting/fragment_tables/ploop_fragment_set/tables/1ang_3_contact_ploop_set_4_v3.json"
    dict_path = "/home/dzorine/phos_binding/pilot_runs/loop_grafting/fragment_tables/ploop_fragment_set/dicts/1ang_3_contact_ploop_set_4_v3.bin"
    rot_table_path = "/home/dzorine/phos_binding/pilot_runs/loop_grafting/fragment_tables/inverse_rotamer/tables/inv_rot_exchi4_1_2_ang_15_deg.json"
    rot_dict_path = "/home/dzorine/phos_binding/pilot_runs/loop_grafting/fragment_tables/inverse_rotamer/dicts/inv_rot_exchi4_1_2_ang_15_deg.bin"
    frag_table, frag_dict = load_table_and_dict(
        table_path, dict_path, np.dtype("i8"), np.dtype("i8")
    )
    rot_table, rot_dict = load_table_and_dict(
        rot_table_path, rot_dict_path, np.dtype("i8"), np.dtype("i8")
    )
    run_pyrosetta_with_flags(flags_file)
    grid_dict = {}
    with open(param_json_path, "r") as f:
        grid_dict = json.load(f)
    # Dumb hack to turn json string keys into ints where possible
    # TODO just dont use int keys in anything you want to json or serialize
    # in a different way

    grid_dict = pythonify(grid_dict)
    for param_dict in subsample_grid(
        int(num_chunks), int(chunk_index), grid_dict
    ):
        # param_dict = dict(param_json)
        # print (*[(k,type(k),v,type(v)) for k,v in param_dict.items()])
        # print (*[(k,type(k),v,type(v)) for val in range(1,1+param_dict["num_helices"]) for dict_ in param_dict[val] for k,v in dict_.items()])
        param_json = json.dumps(param_dict)
        name_hash = abs(hash(param_json))
        pose = param_to_pose(param_dict)
        # debug hack
        pose = pyrosetta.pose_from_file("/home/dzorine/temp/ez_hit.pdb")
        inf = pyrosetta.rosetta.core.pose.PDBInfo(pose)
        inf.name(f"param_bundle_{name_hash}.pdb")
        pose.pdb_info(inf)

        if not len(pose.residues):
            continue
        results = scan_for_n_term_helical_grafts(
            pose, frag_table=frag_table, frag_dict=frag_dict
        )
        results["allowed_res"] = results.apply(
            lambda x: [*range(1, len(pose.residues) + 1)], axis=1
        )

        results["param_json"] = results.apply(lambda x: param_json, axis=1)
        sec_results = scan_for_inv_rot(
            pose,
            results,
            inv_rot_table=rot_table,
            inv_rot_dict=rot_dict,
            feature_xform_label="frag_feature_to_start",
            cart_resl=2,
            ori_resl=15,
        )
        # print (results)
        if not sec_results is None:
            with open(f"{out_dir}/{name_hash}_params.json", "w+") as f:
                json.dump(param_dict, f)
            process_secondary_results(sec_results, pose, out_dir)
        break


if __name__ == "__main__":
    main()
