import logging

# external
import click
import numpy as np
import h5py

# rosetta
import pyrosetta

# AP Moyer
import getpy as gp

from nerf import iNeRF, NeRF, perturb_dofs

# Will Sheffler
from xbin import XformBinner as xb
from homog import hstub

# dzutils
from dzutils.pyrosetta_utils import residue_type_from_name3
from dzutils.pyrosetta_utils.phos_binding.misc_scripts.rotable import (
    save_dict_as_bin,
    consolidate_chis,
    expand_rotamer_array,
    unpack_chis,
    get_dof_templates_from_rotamer_rt_array,
    chis_array_to_rt_array,
    get_new_key_mask_from_hashmap,
    rotamer_rt_array_to_dof_template,
    rotamer_rt_array_to_target_mask,
    fill_dof_template,
    xyzs_to_stub_array,
)
from dzutils.pyrosetta_utils.geometry.pose_xforms import RotamerRTArray


logger = logging.getLogger("make_inverse_rot_tables")
logger.setLevel(logging.DEBUG)


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
@click.option("-g", "--granularity-factor", default=15)
@click.option("-s", "--search-radius", default=5)
def main(
    hashmap_path,
    hdf5_store,
    run_name="inverse_ptr_exchi7_rotamers",
    angstrom_dist_res=1,
    angle_res=15,
    erase=False,
    index_range="all",
    granularity_factor=15,
    search_radius=5,
    output_dir="",
):

    # defint working dirs here:

    start, end = [None, None]
    if index_range != "all":
        start, end = [int(val) for val in index_range.split("-")]
    with h5py.File(hdf5_store, "r") as f:
        store_chis = f["chis"][start:end]
        store_keys = f["key_int"][start:end]
        store_rts = f["rt"][start:end]
        res_type_name = f["chis"].attrs["residue_name"]
        num_chis = f["chis"].attrs["num_chis"]
        store_alignment_atoms = f["chis"].attrs["target_atoms"]
        store_base_atoms = f["chis"].attrs["base_atoms"]
        ideal_chis = f["ideal_chis"][:]

    chis_to_expand = consolidate_chis(
        np.array(store_chis), num_chis=num_chis, round_fraction=np.float(0.01)
    )

    new_chis = expand_rotamer_array(
        chis_to_expand,
        search_radius,
        granularity_factor,
        round_fraction=np.float(0.01),
    )
    # Convert rosetta residue to nerf object
    pyrosetta.init(
        """-out:level 100
        -extra_res_fa /home/dzorine/projects/phos_binding/params_files/p_loop_params/PHY_4_chi.params
        """
    )
    res_type = residue_type_from_name3(res_type_name)
    residue = pyrosetta.rosetta.core.conformation.Residue(res_type, True)
    target_atoms = store_alignment_atoms
    rotamer_rt_array = RotamerRTArray(
        residue=residue,
        base_atoms=list(store_base_atoms),
        target_atoms=list(target_atoms),
        inverse=True,
    )

    dof_template, mask, abc = rotamer_rt_array_to_dof_template(
        rotamer_rt_array
    )

    targ_mask, center_mask = rotamer_rt_array_to_target_mask(rotamer_rt_array)
    target_mask_array = np.array(targ_mask)
    center_mask_array = np.array(center_mask)
    stub_coords = rotamer_rt_array.get_base_xyz()

    dofs = fill_dof_template(dof_template, new_chis, mask)
    base_xform = hstub(
        stub_coords[0, :],
        stub_coords[1, :],
        stub_coords[2, :],
        cen=list(stub_coords[1, :]),
    )
    base_xforms = np.tile(base_xform, (len(dofs), 1, 1))
    # hash the data
    abcs = np.tile(abc, (len(dofs), 1, 1))
    print(f"dofs: {dofs}")
    print(f"abcs: {abcs}")
    print(f"dofs shape : {dofs.shape}")
    print(f"abcs shape : {abcs.shape}")
    new_xyzs = NeRF(abcs, dofs)
    target_stub_array = xyzs_to_stub_array(
        new_xyzs, target_mask_array, center_mask_array
    )
    rts = np.linalg.inv(target_stub_array) @ base_xforms

    binner = xb(cart_resl=angstrom_dist_res, ori_resl=angle_res)
    keys = binner.get_bin_index(np.array(rts))
    key_type = np.dtype("i8")
    value_type = np.dtype("i8")

    new_keys_mask = get_new_key_mask_from_hashmap(
        keys, hashmap_path, key_type=key_type, value_type=value_type
    )
    keys_to_save = np.concatenate((keys[new_keys_mask], store_keys))
    index_vals_to_save = np.arange(0, len(keys_to_save))

    new_hashmap = gp.Dict(key_type, value_type)
    new_hashmap[keys_to_save] = index_vals_to_save

    rts_to_save = np.concatenate((rts[new_keys_mask], store_rts))
    chis_to_save = np.concatenate((new_chis[new_keys_mask], store_chis))

    index_vals_to_save = np.arange(0, len(chis_to_save))
    output_hf5_path = f"{output_dir}/{run_name}.hf5"
    with h5py.File(output_hf5_path, "w") as f:
        print("opened")
        f.create_dataset("index", data=index_vals_to_save)
        f.create_dataset("chis", data=chis_to_save)
        f.create_dataset("key_int", data=keys_to_save)
        f.create_dataset("rt", data=rts_to_save)
        # f.create_dataset("alignment_atoms", data=align_atoms_to_save)
        f.create_dataset("ideal_chis", data=ideal_chis)
        f["chis"].attrs["residue_name"] = res_type_name
        f["chis"].attrs["num_chis"] = num_chis
        f["chis"].attrs["target_atoms"] = store_alignment_atoms
        f["chis"].attrs["base_atoms"] = store_base_atoms

    # Make the base dictionary
    new_hashmap = gp.Dict(key_type, value_type)
    new_hashmap[keys_to_save] = index_vals_to_save

    save_dict_as_bin(output_dir, new_hashmap, run_name, overwrite=erase)


if __name__ == "__main__":
    main()
