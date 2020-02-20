#!/usr/bin/env python
import click

import pyrosetta

from dzutils.pyrosetta_utils import run_pyrosetta_with_flags, safe_load_pdbs
from dzutils.pyrosetta_utils.chain_utils import closed_loop_generator


@click.command()
@click.argument("pdb_list", nargs=-1, type=click.Path(exists=True))
@click.option("-r", "--rosetta-flags-file", default="")
@click.option("-l", "--max-len", default=7)
def main(
    pdb_list,
    rosetta_flags_file="",
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
    for pose in safe_load_pdbs(pdb_list):

        pdb_info_name_str = pose.pdb_info().name()
        pdb_info_basename = pdb_info_name_str.split("/")[-1]
        pdb_name = ".".join(pdb_info_basename.split(".")[:-1])
        for i, working_pose in enumerate(
            closed_loop_generator(
                pose,
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
            working_pose.dump_pdb(f"{pdb_name}_loop_{i}.pdb")


if __name__ == "__main__":
    main()
