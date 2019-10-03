import os

import click

import pyrosetta

from dzutils.pyrosetta_utils import (
    run_pyrosetta_with_flags,
    build_hbond_set,
    hbond_to_residue,
    residues_with_element,
)


@click.command()
@click.argument("pose_pdb", type=click.Path(exists=True))
@click.option("-l", "--res-label", default="phos_contact")
@click.option("-r", "--rosetta-flags-file", default="")
@click.option("-b", "--bb-only", default=False)
def main(pose_pdb, res_label, bb_only=False, rosetta_flags_file=""):
    """
    """
    run_pyrosetta_with_flags(rosetta_flags_file)
    pose = pyrosetta.pose_from_file(pose_pdb)

    hbond_set = build_hbond_set(pose)

    phos_residues = residues_with_element(pose, "P")

    collections = [
        hbond_to_residue(pose, res, hbond_set=hbond_set)
        for res in phos_residues
    ]

    for collection in collections:
        if bb_only:
            collection = [
                b for b in collection if b.don_hatm_is_protein_backbone()
            ]
        for hbond in collection:
            pose.pdb_info().add_reslabel(hbond.don_res(), res_label)

    name = f"{os.path.basename(pose_pdb).split('.')[0]}_labeled_contacts.pdb"
    pose.dump_pdb(name)


if __name__ == "__main__":
    main()
