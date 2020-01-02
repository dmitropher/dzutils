from os.path import isfile

import h5py
import numpy as np
import pandas as pd

# AP Moyer
import getpy as gp
from nerf import iNeRF, NeRF, perturb_dofs

# Will Sheffler
from homog import hstub
from xbin import XformBinner as xb


import pyrosetta

from dzutils.pyrosetta_utils.geometry.pose_xforms import RotamerRTArray


def rt_from_chis(rotamer_rt_array, *chis, base_atoms=(), target_atoms=()):
    rotamer_rt_array.reset_rotamer(
        *chis, base_atoms=base_atoms, target_atoms=target_atoms
    )
    # logger.debug(chis)
    return np.array(rotamer_rt_array)


def save_dict_as_bin(dict_out_dir, gp_dict, dict_name, overwrite=False):
    """
    wraps gp_dict.dump
    """
    dict_out_path = f"{dict_out_dir}/{dict_name}.bin"
    if not overwrite:
        if isfile(dict_out_path):
            raise RuntimeError(
                f"Overwrite set to {overwrite} and file '{dict_out_path}' exists"
            )
    gp_dict.dump(dict_out_path)


def generate_rosetta_rotamer_chis(res_type):
    rots = [
        rot
        for rot in pyrosetta.rosetta.core.pack.rotamer_set.bb_independent_rotamers(
            res_type
        )
    ]
    return [
        tuple(rot.chi(chi) for chi in range(1, rot.nchi() + 1)) for rot in rots
    ]


def rosetta_rot_data(restype, possible_rt_bases):

    chis = generate_rosetta_rotamer_chis(restype)

    residue = pyrosetta.rosetta.core.conformation.Residue(restype, True)
    # possible_rt_bases = pres_bases(residue)

    rotamer_rt = RotamerRTArray(
        residue=residue,
        base_atoms=[2, 1, 2, 3],
        target_atoms=possible_rt_bases[0],
        inverse=True,
    )

    rts, chis_index, alignment_atoms = zip(
        *[
            (
                rt_from_chis(
                    rotamer_rt,
                    *chi_set,
                    base_atoms=[2, 1, 2, 3],
                    target_atoms=atoms,
                ),
                chi_set,
                np.array(atoms),
            )
            for chi_set in chis
            for atoms in possible_rt_bases
        ]
    )

    return rts, chis_index, alignment_atoms


def write_hdf5_data(path, **kargs):
    """
    Takes a dict of hierarchical data and groups

    dataset and group cannot share names (hdf5 stuff)

    format: {name:(group_name,data,{attr_name:attr_value})}
    if you wish to use the root group, "/" or "" are adequate

    for no atributes use the empty dict

    """
    with h5py.File(path, "w") as f:
        for k, (group_str, data, attrs) in kargs.items():
            if group_str and group_str not in f.keys():
                group = f.create_group(group_str)
                group.create_dataset(
                    k,
                    data=data,
                    dtype=data.dtype,
                    maxshape=(None, *data.shape[1:]),
                    chunks=True,
                )
                if attrs:
                    for key, val in attrs.items():
                        group[k].attrs[key] = val

            else:

                f.create_dataset(
                    k,
                    data=data,
                    dtype=data.dtype,
                    maxshape=(None, *data.shape[1:]),
                    chunks=True,
                )
                if attrs:
                    for key, val in attrs.items():
                        f[k].attrs[key] = val


def write_hdf5_rotamer_hash_data(
    path,
    restype,
    rts,
    chis_index,
    alignment_atoms,
    cart_resl=1,
    ori_resl=15,
    key_label="key_int",
    chi_label="chis",
    base_atom_label="base_atoms",
    align_atom_label="target_atoms",
    rt_label="rt",
    ideal=False,
    res_name="",
):
    """
    """
    chi_array = np.array(chis_index)
    binner = xb(cart_resl=cart_resl, ori_resl=ori_resl)
    keys = binner.get_bin_index(np.array(rts))
    with h5py.File(path, "w") as f:
        for label, data in zip(
            (
                "index",
                key_label,
                rt_label,
                chi_label,
                base_atom_label,
                align_atom_label,
            ),
            (
                [*range(1, len(keys) + 1)],
                keys,
                rts,
                chi_array,
                alignment_atoms,
            ),
        ):
            np_data = np.array(data)
            f.create_dataset(label, np_data.shape, data=np_data)
        if ideal:
            # np_data = np.array(chis_index)
            f.create_dataset("ideal_chis", chi_array.shape, data=chi_array)
        f[chi_label].attrs["num_chis"] = chi_array.shape[1]
        f[chi_label].attrs["residue_name"] = restype.name()
        f[key_label].attrs["cart_resl"] = cart_resl
        f[key_label].attrs["ori_resl"] = ori_resl


def write_hdf5_rosetta_rotamer_hash_data(
    path,
    restype,
    cart_resl=1,
    ori_resl=15,
    key_label="key_int",
    chi_label="chis",
    align_atom_label="alignment_atoms",
    rt_label="rt",
    ideal=False,
):
    """
    This is the utility that generates backbone dep rotamers and hashes them

    Bin dimensions and hdf5 labels are configurable, but rotamer generation is
    just whatever ex level rosetta has been set to before calling this
    """
    rts, chis_index, alignment_atoms = rosetta_rot_data(restype)
    write_hdf5_rotamer_hash_data(
        path,
        restype,
        rts,
        chis_index,
        alignment_atoms,
        cart_resl=cart_resl,
        ori_resl=ori_resl,
        key_label=key_label,
        chi_label=chi_label,
        align_atom_label=align_atom_label,
        rt_label=rt_label,
        ideal=ideal,
        res_name=restype.name(),
    )


def build_inv_rot_table_from_rosetta_rots(
    restype,
    cart_resl=1,
    ori_resl=15,
    key_label="key_int",
    chi_label="chis",
    align_atom_label="alignment_atoms",
    rt_label="rt",
):
    """
    """
    rts, chis_index, alignment_atoms = rosetta_rot_data(restype)
    binner = xb(cart_resl=cart_resl, ori_resl=ori_resl)
    keys = binner.get_bin_index(np.array(rts))

    # make the table to assign rotamers to
    inv_rot_table = pd.DataFrame(
        {
            key_label: keys,
            chi_label: chis_index,
            align_atom_label: alignment_atoms,
            rt_label: rts,
        }
    )

    return inv_rot_table


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
    chis = (
        np.array(np.bitwise_and(shifted, mask), dtype=np.float64)
        * round_fraction
    )
    return chis


def consolidate_chis(chis_array, num_chis=3, round_fraction=0.01):
    packed_chis = bit_pack_rotamers(
        chis_array, num_chis, round_fraction=np.float(0.01)
    )
    unique_packed_chis = np.fromiter(set(packed_chis), np.uint64())
    new_chis = unpack_chis(
        unique_packed_chis, num_chis=num_chis, round_fraction=np.float(0.01)
    )
    return new_chis


# @jit(nopython=True,cache=True)
def expand_rotamer_set(
    chi_array, search_radius, granularity_factor, round_fraction
):
    """
    Takes a chi array and does a scan about the rotamer given with given params
    """
    rotamer_set = set([np.uint64(0)][1:])
    num_chi_sets = chi_array.shape[0]
    # logger.debug (num_chi_sets)

    for i in range(num_chi_sets):
        chis = chi_array[i]
        num_chis = len(chis)
        col_space = np.empty((num_chis, granularity_factor), dtype=np.float64)
        for j in range(num_chis):
            col_space[j, :] = np.linspace(
                chis[j] - np.float64(search_radius),
                chis[j] + np.float64(search_radius),
                granularity_factor,
            )

        expanded = np.asarray(col_space, dtype=np.float64)

        # fine_chis = np.array (list(it_cartesian(
        #     expanded
        # )), np.float64)
        fine_chis = np.array(numba_cartesian(expanded))

        uint_rot_repr = bit_pack_rotamers(
            fine_chis, num_chis=num_chis, round_fraction=round_fraction
        )
        if not len(rotamer_set):
            rotamer_set = set(uint_rot_repr)
            continue
        rotamer_set = rotamer_set.union(set(uint_rot_repr))

    return rotamer_set


def expand_rotamer_array(
    chi_array, search_radius, granularity_factor, round_fraction
):
    """
    Takes a chi array and does a scan about the rotamer given with given params
    """
    rotamer_set = set([np.uint64(0)][1:])
    num_chi_sets = chi_array.shape[0]
    # logger.debug (num_chi_sets)
    num_chis = chi_array.shape[-1]
    for i in range(num_chi_sets):
        chis = chi_array[i]

        col_space = np.empty((num_chis, granularity_factor), dtype=np.float64)
        for j in range(num_chis):
            col_space[j, :] = np.linspace(
                chis[j] - np.float64(search_radius),
                chis[j] + np.float64(search_radius),
                granularity_factor,
            )

        expanded = np.asarray(col_space, dtype=np.float64)

        fine_chis = np.array(numba_cartesian(expanded))

        uint_rot_repr = bit_pack_rotamers(
            fine_chis, num_chis=num_chis, round_fraction=round_fraction
        )
        if not len(rotamer_set):
            rotamer_set = set(uint_rot_repr)
            continue
        rotamer_set = rotamer_set.union(set(uint_rot_repr))
    expanded_packed = np.fromiter(rotamer_set, np.uint64())
    expanded_array = unpack_chis(
        expanded_packed, num_chis=num_chis, round_fraction=round_fraction
    )

    return expanded_array


def get_atom_chain_from_restype(res_type, *extra_atoms):
    """
    """
    chi_atoms = res_type.chi_atoms()
    chi_dict = {atoms[1]: list(atoms) for atoms in chi_atoms}

    chain_atom_list = []
    for n in range(1, res_type.natoms() + 1):
        # katom = at.atom(pyrosetta.rosetta.core.id.AtomID(n, 1))
        if n in chi_dict:
            for chi_a in chi_dict[n]:
                if not chi_a in chain_atom_list:
                    chain_atom_list.append(chi_a)

    if extra_atoms:
        for atom in extra_atoms:
            if not atom in chain_atom_list:
                chain_atom_list.append(atom)

    return chain_atom_list


def atom_chain_to_xyzs(rotamer_rt_array, atom_chain, *extra):
    """
    """
    chi_atoms = rotamer_rt_array.residue.type().chi_atoms()
    reduced_extra = [a for a in extra if not a in atom_chain]

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


def rotamer_rt_array_to_dof_template(rotamer_rt_array):
    """
    """
    target_atoms = list(rotamer_rt_array._target_atoms)
    # print(target_atoms)
    atom_chain = get_atom_chain_from_restype(
        rotamer_rt_array.residue.type(), *target_atoms
    )
    # print(atom_chain)
    xyz_coords, mask = atom_chain_to_xyzs(
        rotamer_rt_array, atom_chain, *target_atoms
    )
    dofs = iNeRF(
        xyz_coords[:3][np.newaxis, :],
        xyz_coords[3:][np.newaxis, :],
        degrees=True,
    )
    return dofs, mask, xyz_coords[:3]


def rotamer_rt_array_to_target_mask(rotamer_rt_array, center_index=0):
    """
    Target atoms should be length 3 or 4
    """
    atom_chain = get_atom_chain_from_restype(
        rotamer_rt_array.residue.type(), *rotamer_rt_array._target_atoms
    )
    xyz_coords, mask = atom_chain_to_xyzs(
        rotamer_rt_array, atom_chain, *rotamer_rt_array._target_atoms
    )

    targ_set = rotamer_rt_array._target_atoms[1:]
    targ_mask = [
        targ_set.index(a) + 1 if a in rotamer_rt_array._target_atoms else False
        for a in atom_chain
    ]
    center_atom = rotamer_rt_array._target_atoms[center_index]
    center_mask = [a == center_atom for a in atom_chain][3:]
    return targ_mask, center_mask


def get_dof_templates_from_rotamer_rt_array(
    rotamer_rt_array, target_atoms=[], center_indices=[]
):
    """
    """
    dof_sets = []
    mask_sets = []
    targ_mask_sets = []
    center_mask_sets = []
    xyz_stubs = []
    if not center_indices:
        center_indices = [0] * len(target_atoms)
    for atom_set, center in zip(target_atoms, center_indices):
        atom_chain = get_atom_chain_from_restype(
            rotamer_rt_array.residue.type(), *atom_set
        )
        xyz_coords, mask = atom_chain_to_xyzs(
            rotamer_rt_array, atom_chain, *atom_set
        )

        targ_set = atom_set[1:] if len(atom_set) > 3 else atom_set
        targ_mask = [
            targ_set.index(a) + 1 if a in atom_set else False
            for a in atom_chain
        ]
        center_atom = atom_set[center]
        center_mask = [a == center_atom for a in atom_chain]
        dofs = iNeRF(
            xyz_coords[:3][np.newaxis, :],
            xyz_coords[3:][np.newaxis, :],
            degrees=True,
        )
        dof_sets.append(dofs)
        mask_sets.append(mask)
        targ_mask_sets.append(targ_mask)
        xyz_stubs.append(xyz_coords[:3])
        center_mask_sets.append(center_mask)
    return dof_sets, mask_sets, targ_mask_sets, center_mask_sets, xyz_stubs


def fill_dof_template(dof_template, chis_array, mask):
    """
    """
    filled_dofs = np.tile(dof_template, (chis_array.shape[0], 1, 1))
    filled_dofs[:, mask, 2] = chis_array
    return filled_dofs


def xyzs_to_stub_array(xyzs, target_mask_array, center_mask_array):
    """
    Turns xyz array and some base stub into an RT

    Target mask is an array with the same shape as xyzs. Positions in xyzs where
    the values should be used to generate the target stub should have the atom
    position specified. Atom 1 is used as both center and 1, atom
    """
    a1_mask = target_mask_array[3:] == 1
    a2_mask = target_mask_array[3:] == 2
    a3_mask = target_mask_array[3:] == 3
    center = xyzs[:, center_mask_array, :][:, 0, :]
    a1 = xyzs[:, a1_mask, :][:, 0, :]
    a2 = xyzs[:, a2_mask, :][:, 0, :]
    a3 = xyzs[:, a3_mask, :][:, 0, :]
    target_xforms = hstub(a1, a2, a3, cen=list(center))
    return target_xforms


def template_to_rt(
    chis_array,
    dof_template,
    mask,
    nerf_stub,
    target_stub_mask,
    center_mask,
    base_stub,
):
    """
    """
    filled_dofs = fill_dof_template(dof_template, chis_array, mask)
    abcs = np.tile(nerf_stub, (chis_array.shape[0], 1, 1))
    xyzs = NeRF(abcs, filled_dofs, degrees=True)
    targ_stub_array = xyzs_to_stub_array(xyzs, target_stub_mask, center_mask)
    base_xforms = np.tile(base_stub, (len(targ_stub_array), 1, 1))
    rts = np.linalg.inv(targ_stub_array) @ base_xforms
    return rts


def chis_array_to_rt_array(
    chis_array,
    dofs_templates,
    dofs_masks,
    nerf_stubs,
    target_stub_masks,
    center_mask_sets,
    stub_coords,
):
    """
    Converts an array of chis into an array of RTs

    dofs_templates and dofs_masks here represent the template of dofs for
    the desired atoms chain with the mask indicating the indices where a chi
    value is stored

    stub coords is where the "base" plane of the RT should be computed from

    """
    # apply chis to each dof in array
    all_rts = []

    for i in range(stub_coords.shape[0]):
        base_xform = hstub(
            stub_coords[i][0, :],
            stub_coords[i][1, :],
            stub_coords[i][2, :],
            cen=list(stub_coords[i][1, :]),
        )
        for (
            dof_template,
            mask,
            nerf_stub,
            target_stub_mask,
            center_mask,
        ) in zip(
            dofs_templates,
            dofs_masks,
            nerf_stubs,
            target_stub_masks,
            center_mask_sets,
        ):
            target_stub_mask = np.array(target_stub_mask)
            center_mask = np.array(center_mask)
            rts = template_to_rt(
                chis_array,
                dof_template,
                mask,
                nerf_stub,
                target_stub_mask,
                center_mask,
                base_xform,
            )
            all_rts.append(rts)
    all_rts_array = np.concatenate(all_rts)
    return all_rts_array


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
