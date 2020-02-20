#!/usr/bin/env python
import click

import pyrosetta

from pyrosetta.rosetta.core.select.residue_selector import (
    ResiduePDBInfoHasLabelSelector,
)
from pyrosetta.rosetta.protocols.grafting import return_region

from dzutils.pyrosetta_utils import run_pyrosetta_with_flags
from dzutils.pyrosetta_utils.chain_utils import (
    link_poses,
    closed_loop_generator,
)


@click.command()
@click.argument("pdb_list", nargs=-1, type=click.Path(exists=True))
@click.option("-r", "--rosetta-flags-file", default="")
@click.option("-s", "--start", default="start_label")
@click.option("-e", "--end", default="end_label")
@click.option("-n", "--output-name", default="end_label")
@click.option("-l", "--max-len", default=7)
def main(
    pdb_list,
    rosetta_flags_file="",
    start="start_label",
    end="end_label",
    output_name="",
    max_len=7,
    # additional_output_limit=0,
):
    """
    This program takes a pose and attempts to loop close from start to end

    If these residues are on the same chain, it splits them to different chains

    defaults to using the last value labelled start and the first labled end
    """
    if rosetta_flags_file:
        run_pyrosetta_with_flags(rosetta_flags_file)
    else:
        pyrosetta.init()
    for pose_pdb in pdb_list:
        pose = pyrosetta.pose_from_file(pose_pdb)
        start_label_selector = ResiduePDBInfoHasLabelSelector(start)
        res_with_start_label = [
            i
            for i, is_labeled in enumerate(start_label_selector.apply(pose), 1)
            if is_labeled
        ]
        start_resnum = res_with_start_label[0]
        end_label_selector = ResiduePDBInfoHasLabelSelector(end)
        end_resnum = list(end_label_selector.apply(pose)).index(1) + 1
        # end_resnum = res_with_end_label[0]
        subpose_before_close = return_region(pose.clone(), 1, start_resnum)
        subpose_after_close = return_region(
            pose.clone(), end_resnum, len(pose.residues)
        )
        n_chains = list(subpose_before_close.split_by_chain())
        c_chains = list(subpose_after_close.split_by_chain())
        to_join = link_poses(n_chains[-1], c_chains[0], rechain=True)
        # to_join.dump_pdb("test.pdb")
        pdb_info_name_str = pose.pdb_info().name()
        pdb_info_basename = pdb_info_name_str.split("/")[-1]
        pdb_name = ".".join(pdb_info_basename.split(".")[:-1])
        print(f"chains: {to_join.num_chains()}")
        print(f"start_resnum: {start_resnum}, end_resnum: {end_resnum}")
        for i, working_pose in enumerate(
            closed_loop_generator(
                to_join,
                label="naive_loop",
                database="/home/dzorine/databases/vall.json",
                length=max_len,
                rmsd_tol=0.5,
                cluster_tol=1.75,
                from_chain=1,
                to_chain=2,
            ),
            1,
        ):
            print("new pose generated")
            full_pose = link_poses(
                *n_chains[:-1], working_pose, *c_chains[1:], rechain=True
            )
            full_pose.dump_pdb(
                f"{output_name+'_' if output_name else '' }{pdb_name}_loop_{i}.pdb"
            )


if __name__ == "__main__":
    main()
