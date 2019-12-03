import click
import numpy as np
import h5py

# AP Moyer
import getpy as gp


@click.command()
# @click.argument("hashmap-path", type=click.Path(exists=True))
@click.argument("hdf5-stores", nargs=-1, type=click.Path(exists=True))
@click.option(
    "-f", "--hdf5-store-list", type=click.Path(exists=True), default=""
)
def main(hdf5_stores, hdf5_store_list=""):
    """
    Use the associated hashmap to consolidate stores 1 and 2, dump new hashmap
    """
    if hdf5_store_list:
        hdf5_stores = open(hdf5_store_list).read().splitlines()
    key_type = np.dtype("i8")
    value_type = np.dtype("i8")
    hashmap = gp.Dict(key_type, value_type)
    # hashmap.load(hashmap_path)

    with h5py.File("combined.hf5", "w") as out:
        # first make an expandable copy of the first dset
        # Not real hdf5 support, no groups allowed
        with h5py.File(hdf5_stores[0], "r") as f:
            for key in f.keys():
                data = f[key][:]
                out_set = out.create_dataset(
                    key,
                    data.shape,
                    data=data,
                    maxshape=(None, *data.shape[1:]),
                    chunks=True,
                )
                for attr_key in f[key].attrs.keys():
                    out_set.attrs[attr_key] = f[key].attrs[attr_key]
            keys = f["key_int"][:]
            hashmap[keys] = np.arange(0, len(keys))
        for store_path in hdf5_stores[1:]:
            with h5py.File(store_path, "r") as f:
                store_keys = f["key_int"][:]
                new_keys_mask = hashmap.contains(store_keys) == False
                new_keys = store_keys[new_keys_mask]
                new_len = len(new_keys)
                old_len = out["key_int"].shape[0]
                out["key_int"].resize(old_len + new_len, 0)
                out["key_int"][old_len:] = new_keys
                hashmap[new_keys] = np.arange(old_len, old_len + new_len)
                for store_key in ("chis", "rt"):
                    out[store_key].resize(old_len + new_len, 0)
                    out[store_key][old_len:] = f[store_key][:][new_keys_mask]
                    assert len(out[store_key]) == old_len + len(new_keys)

    dict_out_path = f"consolidated.bin"
    hashmap.dump(dict_out_path)


if __name__ == "__main__":
    main()
