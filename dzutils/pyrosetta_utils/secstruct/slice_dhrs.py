import click
import sys

import pyrosetta

from dzutils.pyrosetta_utils import run_pyrosetta_with_flags
from dzutils.pyrosetta_utils.secstruct.structure import (
    parse_structure_from_dssp,
)


@click.command()
@click.argument("pose_pdb", type=click.Path(exists=True))
@click.option(
    "-n", "--helix-number", type=click.Choice(["l", "m", "t"]), required=True
)
@click.option("-r", "--rosetta-flags-file")
def main(pose_pdb, helix_number="", rosetta_flags_file=""):
    """
    Reads in a pdb and parses the helical DSSP.

    Can request the numbers of middle or last helices, or just the total
    count
    """
    run_pyrosetta_with_flags(rosetta_flags_file)
    pose = pyrosetta.pose_from_pdb(pose_pdb)
    structs = parse_structure_from_dssp(pose, "H")
    size = len(structs)
    if helix_number.lower() == "t":
        sys.stdout.write(size)
    if helix_number.lower() == "l":
        sys.stdout.write(size - 1)
    if helix_number.lower() == "m":
        sys.stdout.write(size // 1 + size % 2)


if __name__ == "__main__":
    main()
