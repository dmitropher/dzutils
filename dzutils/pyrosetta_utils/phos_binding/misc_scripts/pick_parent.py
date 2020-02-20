#!/usr/bin/env python
import click
import pyrosetta
from Bio import pairwise2

from dzutils.pyrosetta_utils import run_pyrosetta_with_flags


@click.command()
@click.argument("parent_file", nargs=1)
@click.argument("pdbs", nargs=-1)
@click.option(
    "-r",
    "--rosetta_flags_file",
    "rosetta_flags_file",
    help="path to the rosetta flags file to use",
)
@click.option(
    "-o",
    "--output-file-path",
    "output_file_path",
    default="summary.csv",
    show_default=True,
    help="path to output parent chart",
)
def main(
    parent_file, pdbs, rosetta_flags_file="", output_file_path="summary.csv"
):
    """
    """
    run_pyrosetta_with_flags(rosetta_flags_file)
    parent_list = []

    with open(parent_file, "r") as f:
        parent_list = list(f.read().splitlines())

    parents = [pyrosetta.pose_from_file(pdb).sequence() for pdb in parent_list]

    # parent_iter = zip(parents, parent_list)
    with open(output_file_path, "w") as f:
        f.write("pdb parent score parent_self_score ratio")
        f.write("\n")
        for pdb in pdbs:
            pose = pyrosetta.pose_from_file(pdb)

            seq = pose.sequence()
            sequence_list = [
                (
                    pairwise2.align.globalxx(
                        seq, ps, score_only=True, one_alignment_only=True
                    ),
                    ps,
                    pdb,
                )
                for ps, pdb in zip(parents, parent_list)
            ]
            seq_score, parent_by_seq, parent_pdb_by_seq = max(
                sequence_list, key=(lambda x: x[0])
            )

            self_score = pairwise2.align.globalxx(
                parent_by_seq,
                parent_by_seq,
                score_only=True,
                one_alignment_only=True,
            )

            print(f"{parent_pdb_by_seq} is the best fit parent of {pdb}")
            output_string = f"{pdb} {parent_pdb_by_seq} {seq_score} {self_score} {seq_score/self_score}"
            f.write(output_string)
            f.write("\n")


if __name__ == "__main__":
    main()
