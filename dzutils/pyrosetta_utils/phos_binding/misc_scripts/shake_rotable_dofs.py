import h5py
import click
import numpy as np

# Will Sheffler
from homog import hstub
from xbin import XformBinner as xb

# AP Moyer
from nerf import NeRF, perturb_dofs
import getpy as gp

# rosetta
import pyrosetta

# dzutils
from dzutils.pyrosetta_utils import residue_type_from_name3
from dzutils.pyrosetta_utils.geometry.pose_xforms import RotamerRTArray
from dzutils.pyrosetta_utils.phos_binding.misc_scripts.rotable import (
    chis_array_to_rt_array,
    consolidate_chis,
    expand_rotamer_set,
    fill_dof_template,
    get_dof_templates_from_rotamer_rt_array,
    get_new_key_mask_from_hashmap,
    rotamer_rt_array_to_dof_template,
    rotamer_rt_array_to_target_mask,
    save_dict_as_bin,
    unpack_chis,
    xyzs_to_stub_array,
)


@click.command()
@click.argument("hashmap_path", type=click.Path(exists=True))
@click.argument("hdf5_store", type=click.Path(exists=True))
@click.option("-a", "--angle-res", default=15.0)
@click.option("-d", "--angstrom-dist-res", default=1.0)
@click.option("-r", "--run-name", default="expand_rotamer_set")
@click.option("-e", "--erase/--no-erase", default=False)
@click.option("-o", "--output-dir", default=".")
@click.option(
    "-n",
    "--index-range",
    default="all",
    help="dash separated indices to expand. May specify 'all'",
)
def main(
    hashmap_path,
    hdf5_store,
    run_name="inv_rot_shaken",
    angstrom_dist_res=1,
    angle_res=15,
    erase=False,
    index_range="all",
    output_dir="",
):
    # need to eliminate this at some point
    pyrosetta.init(
        """-out:level 100
        -extra_res_fa /home/dzorine/projects/phos_binding/params_files/p_loop_params/PHY_4_chi.params
        """
    )
    # defint working dirs here:

    start, end = [None, None]
    if index_range != "all":
        start, end = [int(val) for val in index_range.split("-")]

    outname = f"{output_dir}/{run_name}.hf5"

    with h5py.File(outname, "w") as out:
        # first make an expandable copy of relevant data
        with h5py.File(hdf5_store, "r") as f:
            for key in f.keys():
                data = f[key][start:end]
                out_set = out.create_dataset(
                    key,
                    data.shape,
                    data=data,
                    maxshape=(None, *data.shape[1:]),
                    chunks=True,
                )
                for attr_key in f[key].attrs.keys():
                    out_set.attrs[attr_key] = f[key].attrs[attr_key]
            res_type_name = f["chis"].attrs["residue_name"]
            store_chis = f["chis"][start:end]

        res_type = residue_type_from_name3(res_type_name)
        residue = pyrosetta.rosetta.core.conformation.Residue(res_type, True)
        target_atoms = [
            residue.atom_index(name) for name in ["P", "P", "OH", "O2P"]
        ]  # list(list(res_type.chi_atoms())[-1])  # pres_bases(residue)
        rotamer_rt_array = RotamerRTArray(
            residue=residue,
            base_atoms=[2, 1, 2, 3],
            target_atoms=target_atoms,
            inverse=True,
        )

        dof_template, mask, abc = rotamer_rt_array_to_dof_template(
            rotamer_rt_array
        )
        targ_mask, center_mask = rotamer_rt_array_to_target_mask(
            rotamer_rt_array
        )
        target_mask_array = np.array(targ_mask)
        center_mask_array = np.array(center_mask)
        stub_coords = rotamer_rt_array.get_base_xyz()
        dofs = fill_dof_template(dof_template, store_chis, mask)
        base_xform = hstub(
            stub_coords[0, :],
            stub_coords[1, :],
            stub_coords[2, :],
            cen=list(stub_coords[1, :]),
        )
        base_xforms = np.tile(base_xform, (len(dofs), 1, 1))
        abcs = np.tile(abc, (len(dofs), 1, 1))

        key_type = np.dtype("i8")
        value_type = np.dtype("i8")
        hashmap = gp.Dict(key_type, value_type)
        hashmap.load(hashmap_path)

        binner = xb(cart_resl=angstrom_dist_res, ori_resl=angle_res)
        # Now we can start working on shaking on dofs
        big_num = 10000
        threshold = 0.1
        for i in range(big_num):
            new_dofs = np.array(dofs)
            perturb_dofs(new_dofs)
            new_xyzs = NeRF(abcs, new_dofs)
            target_stub_array = xyzs_to_stub_array(
                new_xyzs, target_mask_array, center_mask_array
            )
            new_rts = np.linalg.inv(target_stub_array) @ base_xforms
            new_rt_keys = binner.get_bin_index(new_rts)
            new_keys_mask = hashmap.contains(new_rt_keys) == False
            new_keys = new_rt_keys[new_keys_mask]

            # update hashmap, store chis/keys/rts
            # if new keys is less than threshold, break
            hashmap[new_keys] = np.arange(
                out["chis"].shape[0], out["chis"].shape[0] + len(new_keys)
            )
            old_end_index = out["chis"].shape[0]
            out["chis"][old_end_index:] = store_chis[new_keys_mask]
            out["key_int"][old_end_index:] = new_keys
            out["rt"][old_end_index:] = new_rts[new_keys_mask]
            portion_new = len(new_keys) / len(new_rt_keys)
            # print(f"portion_new {portion_new} in loop {i}")
            if portion_new < threshold:
                print(f"threshold reached at {portion_new} in loop {i}. ")
                break
        save_dict_as_bin(output_dir, hashmap, run_name, overwrite=erase)


if __name__ == "__main__":
    main()
