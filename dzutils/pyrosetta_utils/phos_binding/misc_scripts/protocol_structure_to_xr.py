#! /usr/bin/env python

# npose necessary evils
import sys

sys.path.append("/home/dzorine/software/npose/")

import click
import numpy as np
import xarray as xr
import npose_util as npu


# end npose


from dzutils.pyrosetta_utils import run_pyrosetta_with_flags, safe_load_pdbs
from dzutils.pyrosetta_utils.secstruct.structure import (
    parse_structure_from_dssp,
)


def get_fragments(pose,):
    """
    Iterates through structures, returns boundaries and type
    """
    frags = parse_structure_from_dssp(pose)
    for frag in frags:

        yield frag.start_pos, frag.end_pos, frag.dssp_type


def get_turns(pose, turn_buffer=3):
    """
    Iterates through L structs, checks before and after
    """

    frags = parse_structure_from_dssp(pose)
    n_frag = len(frags)
    for i in range(1, n_frag - 1):
        frag = frags[i]
        if frag.dssp_type != "L":
            continue
        before = frags[i - 1]
        after = frags[i + 1]

        if (
            before.start_pos > frag.start_pos - turn_buffer
            or after.end_pos < frag.end_pos + turn_buffer
        ):
            continue

        yield frag.start_pos, frag.end_pos, frag.dssp_type, before.dssp_type, after.dssp_type


@click.command()
@click.argument("db_out_path", nargs=1)
@click.argument("pdb_output_path", nargs=1)
@click.argument("pdb_list", nargs=-1)
@click.option(
    "-r",
    "--rosetta-flags-file",
    "rosetta_flags_file",
    default="",
    show_default=True,
)
@click.option(
    "-t",
    "--turn-buffer",
    "turn_buffer",
    default=3,
    show_default=True,
    help="number of residues before and after loop sections to include",
)
def main(
    db_out_path,
    pdb_output_path,
    pdb_list,
    rosetta_flags_file="",
    turn_buffer=3,
):
    ""
    run_pyrosetta_with_flags(rosetta_flags_file)

    full_frag_ds = None
    full_turn_ds = None
    for pose in safe_load_pdbs(pdb_list):

        new_name = f'{".".join(pose.pdb_info().name().split("/")[-1].split(".")[:-1])}'
        pdb_path = f"{pdb_output_path}/{new_name}.pdb"
        pose.dump_pdb(pdb_path)
        try:
            base_npose = npu.npose_from_file(pdb_path)

        except Exception as e:
            print(e)
            continue
        # find fragment streches with continuous secstruct
        # report begining, end, struct before/after (+/-3 res) and struct type

        frags = get_fragments(pose)

        for begin, end, struct in frags:

            residues = ["begin", "end"]
            res_da = xr.DataArray(
                np.array([[begin, end]]),
                dims=("pose_id", "boundary"),
                coords={"boundary": residues},
            )
            struct_da = xr.DataArray([struct], dims=("pose_id"))

            pdb_da = xr.DataArray([pdb_path], dims=("pose_id"))
            # TODO add npose array stuff
            npose_path = f"{pdb_output_path}/{new_name}_{begin}_{end}_{struct}_npose.pdb"
            sub_npose = base_npose[begin * npu.R : npu.R * (end + 1)]
            npu.dump_npdb(sub_npose, npose_path)
            npose_da = xr.DataArray([npose_path], dims=("pose_id"))
            npose_ds = xr.Dataset({"npdb_path": ("pose_id", npose_da)})
            pdb_ds = xr.Dataset({"pdb_path": ("pose_id", pdb_da)})
            res_ds = xr.Dataset({"bound_res": res_da})
            struct_ds = xr.Dataset({"struct": ("pose_id", struct_da)})
            entry = res_ds.merge(pdb_ds)
            entry = entry.merge(npose_ds)
            entry = entry.merge(struct_ds)
            if full_frag_ds:
                full_frag_ds = xr.concat((full_frag_ds, entry), "pose_id")
            else:
                full_frag_ds = entry
            print(full_frag_ds)
        turns = get_turns(pose, turn_buffer=turn_buffer)
        for begin, end, struct, prestruct, post_struct in turns:
            # Turn index is not checked! Returned fragments must have residues at that position

            residues = ["begin", "end"]
            turn_residues = ["begin", "end"]
            res_da = xr.DataArray(
                np.array([[begin, end]]),
                dims=("pose_id", "boundary"),
                coords={"boundary": residues},
            )
            buffer_da = xr.DataArray(
                np.array([[begin - turn_buffer, end + turn_buffer]]),
                dims=("pose_id", "turn_boundary"),
                coords={"turn_boundary": turn_residues},
            )
            bound_struct_da = xr.DataArray(
                np.array([[prestruct, post_struct]]),
                dims=("pose_id", "turn_boundary"),
                coords={"turn_boundary": turn_residues},
            )
            npose_path = f"{pdb_output_path}/{new_name}_{begin}_{end}_{struct}_npose.pdb"
            sub_npose = base_npose[
                (begin - turn_buffer) * npu.R : npu.R * (end + turn_buffer + 1)
            ]
            npu.dump_npdb(sub_npose, npose_path)
            npose_da = xr.DataArray([npose_path], dims=("pose_id"))
            npose_ds = xr.Dataset({"npdb_path": ("pose_id", npose_da)})
            pdb_da = xr.DataArray([pdb_path], dims=("pose_id"))
            struct_da = xr.DataArray([struct], dims=("pose_id"))

            pdb_ds = xr.Dataset({"pdb_path": ("pose_id", pdb_da)})
            struct_ds = xr.Dataset({"struct": ("pose_id", struct_da)})
            b_struct_ds = xr.Dataset({"boundary_struct": bound_struct_da})
            turn_ds = xr.Dataset({"turn_bound_res": buffer_da})
            res_ds = xr.Dataset({"bound_res": res_da})

            entry = res_ds.merge(turn_ds)
            entry = entry.merge(npose_ds)
            entry = entry.merge(pdb_ds)
            entry = entry.merge(struct_ds)
            entry = entry.merge(b_struct_ds)

            if full_turn_ds:
                full_turn_ds = xr.concat((full_turn_ds, entry), "pose_id")
            else:
                full_turn_ds = entry
            print(full_turn_ds)
    full_frag_ds.to_netcdf(f"{db_out_path}/secstruct.ncf")
    full_turn_ds.to_netcdf(f"{db_out_path}/turns.ncf")


if __name__ == "__main__":
    main()
