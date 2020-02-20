#!/usr/bin/env python
import os

import click
import pyrosetta

from dzutils.pyrosetta_utils import run_pyrosetta_with_flags

from dzutils.pyrosetta_utils.chain_utils import circular_permutations_generator


@click.command()
@click.argument("pdb", nargs=1)
@click.argument("resnums", nargs=-1, type=click.INT)
@click.option(
    "-r",
    "--rosetta-flags-file",
    "rosetta_flags_file",
    default="",
    show_default=True,
    help="""optional rosetta flags file""",
)
@click.option(
    "-l",
    "--length",
    "length",
    default=7,
    show_default=True,
    help="""max loop len""",
)
@click.option(
    "-e",
    "--rmsd",
    "rmsd",
    default=0.5,
    help="""rmsd clustering to use for loop lookup""",
)
@click.option(
    "-d",
    "--database",
    "database",
    default="/home/dzorine/databases/vall.json",
    help="""Database to use for loop closure""",
)
@click.option(
    "-c",
    "--cluster-tolerance",
    "cluster_tolerance",
    default=1.75,
    help="""cluster tolerance to use for the loop lookup""",
)
@click.option(
    "-a",
    "--label",
    "label",
    default="circular_permutation_reloop",
    help="""pdb info label to use""",
)
def main(
    pdb,
    resnums,
    rosetta_flags_file="",
    length=5,
    rmsd=0.5,
    database="/home/dzorine/databases/vall.json",
    cluster_tolerance=1.75,
    label="circular_permutation_reloop",
):
    """
    """
    run_pyrosetta_with_flags(rosetta_flags_file)

    pose = pyrosetta.pose_from_file(pdb)
    pdb_name = f"{'.'.join(os.path.basename(pdb).split('.')[:-1])}"

    for i, pose in enumerate(
        circular_permutations_generator(
            pose,
            *resnums,
            many_loops=True,
            label=label,
            database=database,
            length=length,
            rmsd_tol=rmsd,
            cluster_tol=cluster_tolerance,
        )
    ):
        pose.dump_pdb(f"{pdb_name}_circ_perm_{i+1}.pdb")


if __name__ == "__main__":
    main()
