# import xbins
# import bcc
import numpy as np
from numba import jit


def key_to_cell_and_bcc(key):
    return key >> 58, key & 0x3FFFFFFFFFFFFFF


def _nside_prod(self):
    _nside_prod = np.ones(self.ndim, dtype="i8")
    _nside_prod[1:] = np.cumprod(self.nside[:-1])
    return _nside_prod


def sqr(val):
    return val * val


def neighbor_sphere_radius_square_cut(rad, exhalf):
    return (sqr(2 * rad + exhalf) + sqr(2 * (rad + 1) + exhalf)) / 2


def neighbor_radius_square_cut(rad, exhalf):
    return 3 * sqr(2 * rad + exhalf) + 1


def mod(val1, val2):
    return val1 % val2


@jit(nopython=True)
def get_keys_rads(index, rad, rcut, nside_prefsum_, nside):

    exhalf = False
    oddlast3 = True
    sphere = True

    odd = index & 1
    idx0 = np.mod(((index >> 1) / nside_prefsum_), nside)

    lb = -rad - (exhalf and not odd)
    ub = +rad + (exhalf and odd)

    eh = -1 if odd else 0
    oh = 0 if odd else 1
    oddex = odd != exhalf
    l3shift = -(not odd) if oddlast3 else 0
    last3ub = 1 if oddlast3 else 0

    keys, rads = [], []

    for i5 in range(last3ub + 1):
        key5 = nside_prefsum_[5] * (idx0[5] + i5 + l3shift)
        for i4 in range(last3ub + 1):
            key4 = key5 + nside_prefsum_[4] * (idx0[4] + i4 + l3shift)
            for i3 in range(last3ub + 1):
                key3 = key4 + nside_prefsum_[3] * (idx0[3] + i3 + l3shift)

                for i2 in range(lb, ub + 1):
                    key2 = key3 + nside_prefsum_[2] * (idx0[2] + i2)
                    edge2 = i2 == lb if oddex else i2 == ub
                    for i1 in range(lb, ub + 1):
                        key1 = key2 + nside_prefsum_[1] * (idx0[1] + i1)
                        edge1 = edge2 or ((i1 == lb if oddex else i1 == ub))
                        for i0 in range(lb, ub + 1):
                            key0 = key1 + nside_prefsum_[0] * (idx0[0] + i0)
                            key = np.int64(key0) << np.int64(1)
                            edge = edge1 or (i0 == lb if oddex else i0 == ub)

                            inoddlast3 = (
                                i5 + l3shift or i4 + l3shift or i3 + l3shift
                            )
                            skip0 = inoddlast3 and not odd
                            skip1 = inoddlast3 and odd

                            if not skip0:
                                erad = (
                                    np.square(2 * i2 + eh)
                                    + np.square(2 * i1 + eh)
                                    + np.square(2 * i0 + eh)
                                )
                                if (not sphere or erad < rcut) and (
                                    not oddex or not edge
                                ):
                                    keys.append(
                                        np.bitwise_or(
                                            np.int64(key), np.int64(0)
                                        )
                                    )
                                    rads.append(erad)

                            if not skip1:
                                orad = (
                                    np.square(2 * i2 + oh)
                                    + np.square(2 * i1 + oh)
                                    + np.square(2 * i0 + oh)
                                )

                                if (not sphere or orad < rcut) and (
                                    oddex or not edge
                                ):
                                    keys.append(
                                        np.bitwise_or(
                                            np.int64(key), np.int64(1)
                                        )
                                    )
                                    rads.append(orad)

    return keys, rads


def get_smear_bcc_keys_and_rads(
    cell_key, bcc_key, nside_prefsum_, nside, val, radius_cut
):
    index = bcc_key
    rad = radius_cut
    exhalf = False
    # nside_prefsum_ =

    rcut = neighbor_sphere_radius_square_cut(rad, exhalf)
    keys_rads = get_keys_rads(index, rad, rcut, nside_prefsum_, nside)
    return keys_rads


def smear(binner, hashmap, radius):

    exhalf = False
    radius_cut = neighbor_radius_square_cut(radius, exhalf)
    all_keys = []
    all_vals = []
    for key, val in hashmap.items():

        cell_key = key >> 58
        # bcc_key = (key << 6) >> 6
        bcc_key = key & 0x3FFFFFFFFFFFFFF
        # xform = binner.get_bin_center(key)
        # cell_key, f6 = xbin.xform_to_f6(xform)
        # bcc_key = binner.bcc6.get_bin_index(f6)

        bcc_keys, rads = get_smear_bcc_keys_and_rads(
            cell_key,
            bcc_key,
            _nside_prod(binner.bcc6),
            binner.bcc6.nside,
            val,
            radius_cut,
        )
        # print (len (bcc_keys))
        new_bcc_keys = np.array(bcc_keys)
        new_keys = np.bitwise_or(
            np.right_shift(np.int64(cell_key), np.int64(58)),
            np.int64(new_bcc_keys),
        )
        # new_f6 = binner.bcc6.get_bin_center(new_bcc_key)
        # center = f6_to_xform(cell_key, new_f6)
        # new_key = binner.get_bin_index(center)
        all_keys.append(new_keys)
        all_vals.append(np.repeat([val], len(new_keys)))
    hashmap[np.concatenate(all_keys)] = np.concatenate(all_vals)
