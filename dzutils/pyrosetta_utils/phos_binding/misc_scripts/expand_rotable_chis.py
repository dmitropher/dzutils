import logging

# external
import click
import numpy as np
import h5py

# rosetta
import pyrosetta

# AP Moyer
import getpy as gp

# from nerf import iNeRF, NeRF,perturb_dofs

# Will Sheffler
from xbin import XformBinner as xb

# dzutils
from dzutils.pyrosetta_utils import residue_type_from_name3
from dzutils.pyrosetta_utils.phos_binding.misc_scripts.rotable import (
    save_dict_as_bin,
    consolidate_chis,
    expand_rotamer_set,
    unpack_chis,
    get_dof_tempates_from_rotamer_rt_array,
    chis_array_to_rt_array,
    get_new_key_mask_from_hashmap,
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
        store_alignment_atoms = f["alignment_atoms"][start:end]
        ideal_chis = f["ideal_chis"][:]

    chis_to_expand = consolidate_chis(
        np.array(store_chis), num_chis=num_chis, round_fraction=np.float(0.01)
    )
    new_chis_set = expand_rotamer_set(
        chis_to_expand,
        search_radius=search_radius,
        granularity_factor=granularity_factor,
        round_fraction=np.float(0.01),
    )
    packed_chis = np.fromiter(new_chis_set, np.uint64())

    new_chis = unpack_chis(
        packed_chis, num_chis=num_chis, round_fraction=np.float(0.01)
    )

    # Convert rosetta residue to nerf object
    pyrosetta.init(
        """-out:level 100
        -extra_res_fa /home/dzorine/projects/phos_binding/params_files/p_loop_params/PHY_4_chi.params
        """
    )
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

    dof_sets, mask_sets, targ_mask_sets, center_mask_sets, nerf_xyzs = get_dof_tempates_from_rotamer_rt_array(
        rotamer_rt_array, [target_atoms]
    )
    stub_coords_array = np.array([rotamer_rt_array.get_base_xyz()])
    rts = chis_array_to_rt_array(
        new_chis,
        dof_sets,
        mask_sets,
        nerf_xyzs,
        targ_mask_sets,
        center_mask_sets,
        stub_coords_array,
    )

    # hash the data

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

    # make the array representing the chis which have been screened
    expanded_chis = np.tile(new_chis, (len(dof_sets), 1))

    rts_to_save = np.concatenate((rts[new_keys_mask], store_rts))
    chis_to_save = np.concatenate((expanded_chis[new_keys_mask], store_chis))
    new_align_atoms = np.repeat(
        np.array([target_atoms]), len(new_chis), axis=0
    )
    align_atoms_to_save = np.concatenate(
        (new_align_atoms[new_keys_mask], store_alignment_atoms)
    )
    #
    #
    # stubby = rotamer_rt_array.get_base_xyz()
    # base_xform = hstub(
    #     stubby[0, :],
    #     stubby[1, :],
    #     stubby[2, :],
    #     cen=list(stubby[1, :]),
    # )
    # for template,mask, target_stub_mask,center_mask,nerf_stub in zip(dof_sets,mask_sets,targ_mask_sets,center_mask_sets,nerf_xyzs):
    #
    #     # filled_dofs = fill_dof_template(template,new_chis,mask)
    #     for i in range (100):
    #         #perturb and deposit, break if you reach 90% density
    #         dof_copy = np.array(template)
    #         perturb_dofs(dof_copy,bond_length_range=0.1)
    #         perturbed_rts = template_to_rt(
    #             new_chis,
    #             dof_copy,
    #             mask,
    #             nerf_stub,
    #             target_stub_mask,
    #             center_mask,
    #             base_xform
    #         )
    #         perturbed_keys = binner.get_bin_index(perturbed_rts)
    #         perturbed_keys_mask = new_hashmap.contains(perturbed_keys) ==False
    #         new_keys

    index_vals_to_save = np.arange(0, len(chis_to_save))
    output_hf5_path = f"{output_dir}/{run_name}.hf5"
    with h5py.File(output_hf5_path, "w") as f:
        print("opened")
        f.create_dataset("index", data=index_vals_to_save)
        f.create_dataset("chis", data=chis_to_save)
        f.create_dataset("key_int", data=keys_to_save)
        f.create_dataset("rt", data=rts_to_save)
        f.create_dataset("alignment_atoms", data=align_atoms_to_save)
        f.create_dataset("ideal_chis", data=ideal_chis)
        f["chis"].attrs["residue_name"] = res_type_name
        f["chis"].attrs["num_chis"] = num_chis

    # Make the base dictionary
    new_hashmap = gp.Dict(key_type, value_type)
    new_hashmap[keys_to_save] = index_vals_to_save

    save_dict_as_bin(output_dir, new_hashmap, run_name, overwrite=erase)


if __name__ == "__main__":
    main()
