import logging

# external
import click
import numpy as np
import h5py

# rosetta
import pyrosetta

# AP Moyer
import getpy as gp
from nerf import iNeRF, NeRF

# Will Sheffler
from xbin import XformBinner as xb
from homog import hstub

# dzutils
from dzutils.pyrosetta_utils import residue_type_from_name3
from dzutils.pyrosetta_utils.phos_binding import pres_bases
from dzutils.pyrosetta_utils.phos_binding.misc_scripts import save_dict_as_bin
from dzutils.pyrosetta_utils.geometry.pose_xforms import RotamerRTArray


logger = logging.getLogger("make_inverse_rot_tables")
logger.setLevel(logging.DEBUG)


# @jit(nopython=True)
def numba_cartesian(input_arrays, out=None):
    """
    Generate a cartesian product of input arrays.

    Parameters
    ----------
    arrays : list of array-like
        1-D arrays to form the cartesian product of.
    out : ndarray
        Array to place the cartesian product in.

    Returns
    -------
    out : ndarray
        2-D array of shape (M, len(arrays)) containing cartesian products
        formed of input arrays.

    Examples
    --------
    >>> numba_cartesian(([1, 2, 3], [4, 5], [6, 7]))
    array([[1, 4, 6],
           [1, 4, 7],
           [1, 5, 6],
           [1, 5, 7],
           [2, 4, 6],
           [2, 4, 7],
           [2, 5, 6],
           [2, 5, 7],
           [3, 4, 6],
           [3, 4, 7],
           [3, 5, 6],
           [3, 5, 7]])

    """

    arrays = np.empty(
        (len(input_arrays), len(input_arrays[0])), dtype=np.float64
    )
    for i in range(len(input_arrays)):
        # print (arrays[i])
        # print (input_arrays[i])
        arrays[i] = input_arrays[i]
    dtype = arrays[0].dtype

    n = np.prod(
        np.asarray(
            [arrays[i].size for i in range(len(input_arrays))], dtype=np.uint64
        )
    )

    if out is None:
        out = np.zeros((n, len(arrays)), dtype=dtype)

    m = int(n // arrays[0].size)

    out[:, 0] = np.repeat(arrays[0], m)
    array_len = len(input_arrays[1:])
    logger.debug(m)
    if array_len > 0:
        numba_cartesian(arrays[1:], out=out[0:m, 1:])
        for j in range(1, arrays[0].size):
            out[j * m : (j + 1) * m, 1:] = out[0:m, 1:]
    return out


# @jit(nopython=True,cache=True)
def bit_pack_rotamers(chis_array, num_chis=3, round_fraction=0.01):
    """
    Converts the given chis to a numpy uint64

    Uses the formula: left shift the value by 64//(num_chis)*(num of chi - 1)

    bitwise or over the list
    """
    scaling = np.float64(360 / round_fraction)
    if (scaling) > np.power(np.uint64(2), np.uint64(64 // num_chis)):
        raise ValueError(
            "given round_fraction and number of chis cannot fit into uint 64"
        )
    # shape = chis_array.shape[1]
    # uint_chis_array = (np.mod(np.divide(chis_array ,round_fraction) ,scaling)).astype(np.uint64)

    for i in range(num_chis):
        chis_array[:, i] = np.mod(
            np.divide(chis_array[:, i], round_fraction), scaling
        )
    uint_chis_array = chis_array.astype(np.uint64)
    for i in range(num_chis):
        uint_chis_array[:, i] = np.left_shift(
            uint_chis_array[:, i], np.uint64(i * (64 // num_chis))
        )
    packed = uint_chis_array[:, 0]
    for i in range(1, num_chis):
        packed = np.bitwise_or(packed, uint_chis_array[:, i])
    return packed


def unpack_chis(packed_chis, num_chis=3, round_fraction=0.01):
    """
    Takes a np.array of bitpacked chis and unpacks them
    """

    unpacked_transpose = np.array(
        [
            np.right_shift(packed_chis, i * (64 // num_chis))
            for i in range(num_chis)
        ]
    )
    shifted = np.transpose(unpacked_transpose)
    mask = np.power(
        np.uint64(2), np.uint64(64 // num_chis), dtype=np.uint64
    ) - np.uint64(1)
    chis = np.array(
        np.bitwise_and(shifted, mask) * round_fraction, dtype=np.float64
    )
    return chis


# @jit(nopython=True,cache=True)
def expand_rotamer_set(
    chi_array, search_radius, granularity_factor, round_fraction
):
    rotamer_set = set([np.uint64(0)][1:])
    num_chi_sets = chi_array.shape[0]
    # logger.debug (num_chi_sets)

    for i in range(num_chi_sets):
        chis = chi_array[i]
        num_chis = len(chis)
        if (
            bit_pack_rotamers(
                np.array([chis]),
                num_chis=num_chis,
                round_fraction=round_fraction,
            )[0]
            in rotamer_set
        ):
            continue

        # chi_array = np.asarray(chis)
        col_space = np.empty((num_chis, granularity_factor), dtype=np.float64)
        for j in range(num_chis):
            col_space[j, :] = np.linspace(
                chis[j] - np.float64(search_radius),
                chis[j] + np.float64(search_radius),
                granularity_factor,
            )

        expanded = np.asarray(col_space, dtype=np.float64)
        logger.debug("expanded:")
        logger.debug(f"{expanded}")
        # fine_chis = np.array (list(it_cartesian(
        #     expanded
        # )), np.float64)
        fine_chis = np.array(numba_cartesian(expanded))
        logger.debug("fine_chis: ")
        logger.debug(f"{fine_chis}")
        uint_rot_repr = bit_pack_rotamers(
            fine_chis, num_chis=num_chis, round_fraction=round_fraction
        )
        if not len(rotamer_set):
            rotamer_set = set(uint_rot_repr)
            continue
        rotamer_set = rotamer_set.union(set(uint_rot_repr))
    return rotamer_set


def get_atom_chain_from_restype(res_type, *extra_atoms):
    """
    """
    chi_atoms = res_type.chi_atoms()
    # all_atoms = [*list(chi_atoms), extra_atoms]

    # Dummy pose and residue to get an atom tree to play with
    res = pyrosetta.rosetta.core.conformation.Residue(res_type, True)
    pose = pyrosetta.rosetta.core.pose.Pose()
    pose.append_residue_by_bond(res)
    # got our atom tree, but can't delete the pose because of the way it works
    at = pose.atom_tree()

    chain_atom_list = []
    for n in range(1, res_type.natoms() + 1):
        found = False
        katom = at.atom(pyrosetta.rosetta.core.id.AtomID(n, 1))
        for atoms in chi_atoms:
            for chi_a in atoms:
                if (
                    katom.downstream(
                        at.atom(pyrosetta.rosetta.core.id.AtomID(chi_a, 1))
                    )
                    or n == chi_a
                ):
                    chain_atom_list.append(n)
                    found = True
                    break
            if found:
                break
    if extra_atoms:
        found = False
        katom = at.atom(pyrosetta.rosetta.core.id.AtomID(n, 1))
        for atoms in extra_atoms:
            if (
                katom.downstream(
                    at.atom(pyrosetta.rosetta.core.id.AtomID(chi_a, 1))
                )
                or n == chi_a
            ):
                if not n in chain_atom_list:
                    chain_atom_list.append(n)
                    found = True
                    break

    return chain_atom_list


def atom_chain_to_xyzs(rotamer_rt_array, atom_chain, *extra):
    """
    """
    chi_atoms = rotamer_rt_array.residue.type().chi_atoms()
    reduced_extra = [a for a in extra if not a in atom_chain]
    # print(reduced_extra)
    atom_chain.extend(reduced_extra)
    xyzs = rotamer_rt_array.get_xyz(atom_chain)
    mask = [
        any(
            all(atom in chis for atom in atom_chain[i : i + 4])
            for chis in chi_atoms
        )
        for i in range(len(atom_chain) - 3)
    ]
    return xyzs, mask


def chis_array_to_rt_array(
    chis_array, dofs_templates, dofs_masks, target_stub_masks, stub_coords
):
    """
    Converts an array of chis into an array of RTs

    dofs_templates and dofs_masks here represent the template of dofs for
    the desired atoms chain with the mask indicating the indices where a chi
    value is stored

    stub coords is where the "base" plane of the RT should be computed from
    (also the stub coordinates of the NeRF)

    """
    # apply chis to each dof in array
    all_rts = []

    for i in range(stub_coords.shape[0]):
        base_xform = hstub(
            stub_coords[i][:, 0], stub_coords[i][:, 1], stub_coords[i][:, 2]
        )
        for dof_template, mask, target_stub_mask in zip(
            dofs_templates, dofs_masks, target_stub_masks
        ):
            filled_dofs = np.tile(dof_template, (chis_array.shape[0], 1, 1))
            filled_dofs[:, mask, 2] = chis_array
            abcs = np.tile(stub_coords[i], (chis_array.shape[0], 1, 1))
            xyzs = NeRF(abcs, filled_dofs, degrees=True)
            three_atom_stubs = xyzs[:, target_stub_mask[3:], :]
            target_xforms = hstub(
                three_atom_stubs[:, 0, :],
                three_atom_stubs[:, 1, :],
                three_atom_stubs[:, 2, :],
            )
            base_xforms = np.tile(base_xform, (len(target_xforms), 1, 1))
            rts = np.linalg.inv(target_xforms) @ base_xforms
            all_rts.append(rts)
    all_rts_array = np.concatenate(all_rts)
    return all_rts_array


def get_dof_tempates_from_rotamer_rt_array(rotamer_rt_array, target_atoms=[]):
    """
    """
    dof_sets = []
    mask_sets = []
    targ_mask_sets = []
    for base in target_atoms:
        atom_chain = get_atom_chain_from_restype(
            rotamer_rt_array.residue.type(), *base
        )
        xyz_coords, mask = atom_chain_to_xyzs(
            rotamer_rt_array, atom_chain, *base
        )
        targ_mask = [a in base for a in atom_chain]
        dofs = iNeRF(
            xyz_coords[:3][np.newaxis, :],
            xyz_coords[3:][np.newaxis, :],
            degrees=False,
        )
        dof_sets.append(dofs)
        mask_sets.append(mask)
        targ_mask_sets.append(targ_mask)
    return dof_sets, mask_sets, targ_mask_sets


def get_new_key_mask_from_hashmap(
    keys, hashmap_path, key_type=np.dtype("i8"), value_type=np.dtype("i8")
):
    """
    """
    key_type = np.dtype("i8")
    value_type = np.dtype("i8")
    hashmap = gp.Dict(key_type, value_type)
    hashmap.load(hashmap_path)
    return hashmap.contains(np.array(keys)) == False


@click.command()
@click.argument("hashmap_path", type=click.Path(exists=True))
@click.argument("hdf5_store", type=click.Path(exists=True))
@click.option("-a", "--angle-res", default=15.0)
@click.option("-d", "--angstrom-dist-res", default=1.0)
@click.option("-r", "--run-name", default="inverse_ptr_exchi7_rotamers")
@click.option("-e", "--erase/--no-erase", default=False)
@click.option("-o", "--output-name", default="")
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
    output_name="",
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

    packed_store_chis = bit_pack_rotamers(store_chis, num_chis)
    packed_store_chis_unique = np.fromiter(set(packed_store_chis), np.uint64())
    chis_to_expand = unpack_chis(
        packed_store_chis_unique,
        num_chis=num_chis,
        round_fraction=np.float(0.01),
    )
    new_chis_set = expand_rotamer_set(
        chis_to_expand,
        search_radius=search_radius,
        granularity_factor=granularity_factor,
        round_fraction=np.float(0.01),
    )
    packed_chis = np.fromiter(new_chis_set, np.uint64)

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
    base_atoms = pres_bases(residue)
    rotamer_rt_array = RotamerRTArray(
        residue=residue, target_atoms=base_atoms[0], inverse=True
    )

    dof_sets, mask_sets, targ_mask_sets = get_dof_tempates_from_rotamer_rt_array(
        rotamer_rt_array, base_atoms
    )
    stub_coords_array = np.array([rotamer_rt_array.get_base_xyz()])
    rts = chis_array_to_rt_array(
        new_chis, dof_sets, mask_sets, targ_mask_sets, stub_coords_array
    )

    # hash the data
    # print(len(rts), len(dof_sets) * len(new_chis))
    binner = xb(cart_resl=angstrom_dist_res, ori_resl=angle_res)
    keys = binner.get_bin_index(np.array(rts))
    key_type = np.dtype("i8")
    value_type = np.dtype("i8")
    new_keys_mask = get_new_key_mask_from_hashmap(
        keys, hashmap_path, key_type=key_type, value_type=value_type
    )

    # make the array representing the chis which have been screened
    expanded_chis = np.tile(new_chis, (len(dof_sets), 1))

    rts_to_save = np.concatenate((rts[new_keys_mask], store_rts))
    keys_to_save = np.concatenate((keys[new_keys_mask], store_keys))
    chis_to_save = np.concatenate((expanded_chis[new_keys_mask], store_chis))
    new_align_atoms = np.repeat(np.array(base_atoms), len(new_chis), axis=0)
    align_atoms_to_save = np.concatenate(
        (new_align_atoms[new_keys_mask], store_alignment_atoms)
    )
    index_vals_to_save = np.arange(0, len(chis_to_save))

    output_hf5_path = "data_store.hf5"
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
    save_dict_as_bin(
        "/home/dzorine/temp/table_expansion/",
        new_hashmap,
        "updated_expanded_1",
        overwrite=erase,
    )


if __name__ == "__main__":
    main()
