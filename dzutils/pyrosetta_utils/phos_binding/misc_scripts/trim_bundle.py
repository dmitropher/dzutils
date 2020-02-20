#!/usr/bin/env python

import click

from pyrosetta.rosetta.protocols.grafting import return_region

from dzutils.pyrosetta_utils import run_pyrosetta_with_flags, safe_load_pdbs

from dzutils.pyrosetta_utils.chain_utils import link_poses

from dzutils.pyrosetta_utils.secstruct.structure import (
    parse_structure_from_dssp,
)


@click.command()
@click.argument("trim_len", nargs=1, type=click.INT)
@click.argument("pdbs", nargs=-1)
@click.option(
    "-r",
    "--rosetta-flags-file",
    "rosetta_flags_file",
    default="",
    show_default=True,
)
def main(trim_len, pdbs, rosetta_flags_file=""):
    """
    """
    run_pyrosetta_with_flags(rosetta_flags_file)
    for pose in safe_load_pdbs(pdbs):
        structs = parse_structure_from_dssp(pose, "H")
        pose_len = len(pose.residues)
        slen = len(structs)
        n_hairpins = (slen + 1) // 2
        hairpins = []
        for i in range(n_hairpins):
            struct1 = structs[2 * i]
            start = struct1.end_pos - trim_len
            if start < struct1.start_pos:
                print(
                    f"this bundle cannot be trimmed by this factor: {trim_len} struct1.start_pos {struct1.start_pos}"
                )
                break
            s2i = 2 * i + 1
            end = pose_len
            if s2i < slen:
                end = min(structs[s2i].start_pos + trim_len, end)
            hairpins.append(return_region(pose.clone(), start, end))
        linked = link_poses(*hairpins, rechain=True)
        pdb_name = f'{".".join(pose.pdb_info().name().split("/")[-1].split(".")[:-1])}_trimmed.pdb'
        linked.dump_pdb(pdb_name)


if __name__ == "__main__":
    main()
