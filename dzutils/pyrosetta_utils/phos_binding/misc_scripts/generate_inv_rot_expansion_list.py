#!/usr/bin/env python

import h5py
import click


@click.command()
@click.argument("ranges_file", nargs=1)
@click.argument("hdf5_store", nargs=1, type=click.Path(exists=True))
@click.option("-c", "--chunk-size", default=150)
def main(ranges_file, hdf5_store, chunk_size=150):
    """
    """

    with h5py.File(hdf5_store, "r") as f:
        num_chi_sets = f["chis"].shape[0]

    ranges = []
    for start in range(0, num_chi_sets, chunk_size):
        range_string = f"{start}-{start+chunk_size-1}"

        ranges.append(range_string)

    with open(ranges_file, "w") as f:
        f.write("\n".join(ranges))
        f.write("\n")


if __name__ == "__main__":
    main()
