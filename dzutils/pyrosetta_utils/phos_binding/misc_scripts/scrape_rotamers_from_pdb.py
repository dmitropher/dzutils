import os
from itertools import chain
from json import dump

import click
import pyrosetta

from dzutils.pyrosetta_utils import residues_by_name, run_pyrosetta_with_flags


def maybe_load(pdb):
    """
    Dumb hack so rosettta can tolerate pdbs it can't load
    """
    print(f"loading pdb: {pdb}")
    try:
        pose = pyrosetta.pose_from_file(pdb)
        return pose
    except Exception as e:
        print(e)
        return pyrosetta.pose_from_sequence("A")


def scrape_rotamers_by_name3(pose, name):
    """
    returns a list of tuples: (resname, chi_1,chi_2,...,chi_n)
    """
    pairs = residues_by_name(pose, name, include_res=True)
    if not pairs:
        return []
    residues, resnums = zip(*pairs)

    return [
        (name, *[residue.chi(i) for i in range(1, residue.nchi() + 1)])
        for residue in residues
    ]


def scrape_rotamers_by_name3_from_pdb(pdb_file, name):
    return scrape_rotamers_by_name3(maybe_load(pdb_file), name)


def scrape_rotamers_from_pdbs(name, *pdb_files):
    return list(
        chain(
            *(
                scrape_rotamers_by_name3_from_pdb(pdb, name)
                for pdb in pdb_files
            )
        )
    )


@click.command()
@click.argument("pdb_files", type=click.Path(exists=True), nargs=-1)
@click.argument("residue_name3", nargs=1)
@click.option("-o", "--out-dir", default="")
def main(pdb_files, residue_name3, out_dir):

    run_pyrosetta_with_flags(
        "/home/dzorine/projects/phos_binding/run_files/scrape_ptr.flags"
    )

    rotamers = scrape_rotamers_from_pdbs(residue_name3, *pdb_files)
    name_hash = hash(tuple(pdb_files))
    with open(
        f"{out_dir if out_dir else os.getcwd()}/scraped_{residue_name3}_rotamers_{name_hash}.json",
        "w+",
    ) as f:
        dump(rotamers, f)


if __name__ == "__main__":
    main()
