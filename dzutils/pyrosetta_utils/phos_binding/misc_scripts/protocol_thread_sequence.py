#! /usr/bin/env python
import click

import pyrosetta

from dzutils.pyrosetta_utils import (
    safe_load_pdbs,
    parse_label,
    run_pyrosetta_with_flags,
)


@click.command()
@click.argument("seq", nargs=1)
@click.argument("pdbs", nargs=-1)
@click.option(
    "-s",
    "--start-label",
    "start_label",
    help="the string to start attempting threadings at",
    default="",
    show_default=True,
)
@click.option(
    "-e",
    "--end-label",
    "end_label",
    help="the string to end attempting threadings at",
    default="",
    show_default=True,
)
@click.option(
    "-r",
    "--rosetta-flags-file",
    "rosetta_flags_file",
    default="",
    show_default=True,
)
def main(seq, pdbs, start_label="", end_label="", rosetta_flags_file=""):
    """
    Takes pdbs and a seq, threads the seq from every possible position

    Can constrain positions with a residue label, first start_label and last
    end_label are used if multiple are provided
    """

    run_pyrosetta_with_flags(rosetta_flags_file)
    for pose in safe_load_pdbs(pdbs):
        base_name = f'{".".join(pose.pdb_info().name().split("/")[-1].split(".")[:-1])}'
        seq_len = len(seq)
        start_res = 1
        end_res = len(pose.residues) - seq_len
        if start_label:
            resnums = parse_label(pose, start_label)
            start_res = resnums[0]
        if end_label:
            resnums = parse_label(pose, end_label)
            end_res = resnums[-1]
        for i in range(start_res, end_res + 1):
            threader = pyrosetta.rosetta.protocols.simple_moves.SimpleThreadingMover(
                seq, i
            )
            working = pose.clone()
            threader.apply(working)
            working.dump_pdb(f"{base_name}_thread_{i}.pdb")


if __name__ == "__main__":
    main()
