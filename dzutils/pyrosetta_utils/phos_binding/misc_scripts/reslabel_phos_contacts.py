#! /usr/bin/env python
import os

import click

import pyrosetta

from dzutils.pyrosetta_utils import (
    run_pyrosetta_with_flags,
    build_hbond_set,
    hbond_to_residue,
    residues_with_element,
)


def label_pres(pose, label, bb=True, sidechain=True):
    """
    Adds label to bb or sidechain contact residues in place

    if a residue makes both bb and sidechain contacts it will be labeled by either
    TODO: total exclusive sidechain/bb option
    """
    hbond_set = build_hbond_set(pose)

    phos_residues = residues_with_element(pose, "P")

    collections = [
        hbond_to_residue(pose, res, hbond_set=hbond_set)
        for res in phos_residues
    ]
    if not (bb or sidechain):
        raise AttributeError(
            "either bb or sidechain contacts must be selected"
        )
    for collection in collections:
        if bb and sidechain:
            for hbond in collection:
                pose.pdb_info().add_reslabel(hbond.don_res(), label)
            return

        if bb:
            collection = [
                b for b in collection if b.don_hatm_is_protein_backbone()
            ]
        if sidechain:
            collection = [
                b for b in collection if not b.don_hatm_is_protein_backbone()
            ]
        for hbond in collection:
            pose.pdb_info().add_reslabel(hbond.don_res(), label)


@click.command()
@click.argument("pdbs", nargs=-1, type=click.Path(exists=True))
@click.option("-l", "--res-label", default="phos_contact")
@click.option("-r", "--rosetta-flags-file", default="")
@click.option("--bb/--no-bb", default=True)
@click.option("--sidechain/--no-sidechain", default=True)
def main(
    pdbs,
    res_label="phos_contact",
    bb=True,
    sidechain=True,
    rosetta_flags_file="",
):
    """
    """
    run_pyrosetta_with_flags(rosetta_flags_file)
    for pdb in pdbs:
        pose = pyrosetta.pose_from_file(pdb)

        label_pres(pose, res_label, bb=bb, sidechain=sidechain)

        name = f"{'.'.join(os.path.basename(pdb).split('.')[:-1])}{'_'+res_label if res_label else ''}.pdb"
        pose.dump_pdb(name)


if __name__ == "__main__":
    main()
