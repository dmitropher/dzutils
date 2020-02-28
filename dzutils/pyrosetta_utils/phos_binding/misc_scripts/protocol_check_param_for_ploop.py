#!/usr/bin/env python
import os
import click

import numpy as np
import numba
import xarray as xr
import json
import getpy as gp

from xbin import XformBinner
from homog import hstub

import pyrosetta

from dzutils.pyrosetta_utils.geometry.parametric import (
    helix_bundle_maker_wrapper,
    # PyBundleGridSampler,
    # sampler_from_json
)

from dzutils.pyrosetta_utils import run_pyrosetta_with_flags

from dzutils import pythonify

def get_bb_xyzs(pose, *resnums):
    """
    """
    xyzs = []
    for num in resnums:
        res = pose.residue(num)
        xyzs.append(
            [list(res.xyz("N")), list(res.xyz("CA")), list(res.xyz("C"))]
        )
    return xyzs


@click.command()
@click.argument("grids_path", nargs=1)
@click.argument("hashmap_path", nargs=1)
@click.argument("fragment_cdf_path", nargs=1)
@click.option(
    "-r",
    "--rosetta-flags-file",
    "rosetta_flags_file",
    default="",
    show_default=True,
)
@click.option(
    "-m",
    "--mask-values",
    "mask_values",
    default=[],
    show_default=True,
    multiple=True,
    help="helper option to allow for a .where(ds['key'] == 'val'). This option can be used multiple times for several k,v pairs",
)
@click.option(
    "-s",
    "--sel-values",
    "sel_values",
    default=[],
    show_default=True,
    multiple=True,
)
@click.option(
    "-o",
    "--out-dir",
    "out_dir",
    default=".",
    show_default=True,
)
def main(
    grids_path,
    hashmap_path,
    fragment_cdf_path,
    rosetta_flags_file="",
    mask_values=[],
    sel_values=[],
    out_dir=".",
):
    ""

    ds = xr.open_dataset(fragment_cdf_path)
    run_pyrosetta_with_flags(
        rosetta_flags_file
    )
    for input_string in sel_values:
        k,v = input_string.split(" ")
        ds = ds.loc[{k:v}]

    #use where to drop on the conditions in the option (this casts)
    for input_string in mask_values:
        k,v = input_string.split(" ")
        dse = ds.where(ds[k] == v,drop=True).squeeze()
    #fix the casting from .where
    for var in ds.variables:
        dt = ds[var].dtype
        dse[var] = dse[var].astype(dt)
    ds = dse
    #this line was meant to read a params to generate a grid not line delimited param value dicts (grid)
    #im so mad at myself for the terrible naming convention I have for grids
    #TODO rename grids to params and params to grids
    # sampler = sampler_from_json(grids_path)
    binner = XformBinner()
    hashmap = gp.Dict(np.int64,np.int64)
    hashmap.load(hashmap_path)
    rts = ds.rt
    rt_dim, = ds.rt.rts.shape
    rt_count = rts.rts.shape[0]*rts.pdb_path.shape[0]
    # gen = sampler.all_helix_grids()
    print(grids_path)
    out_dir_abspath = os.path.abspath(out_dir)
    print(out_dir_abspath)
    with open (grids_path,"r") as f:

        gen = (pythonify(json.loads(k)) for k in list(f.read().split("\n")) if k)
        for i,grid in enumerate(gen,1):
            repeat_factor = rt_dim*(grid["helix_length"]*grid["num_helices"])**2
            paths = np.repeat(ds["pdb_paths"],repeat_factor)
            start_res = np.repeat(ds.sel(anchor_res="begin")["res_pairs"].squeeze(),repeat_factor,axis=0)
            end_res = np.repeat(ds.sel(anchor_res="end")["res_pairs"].squeeze(),repeat_factor,axis=0)
            p_res = np.repeat(ds.sel(anchor_res="p_res")["res_pairs"].squeeze(),repeat_factor,axis=0)
            bun = helix_bundle_maker_wrapper(grid["helix_length"], *[grid[i] for i in range(1,grid["num_helices"]+1)], degrees=True)
            p = pyrosetta.rosetta.core.pose.Pose()
            bun.apply(p)
            plen = len(p.residues)
            stub_xyzs = np.array(get_bb_xyzs(p, *range(1,plen+1)))
            bb_stubs = hstub(
            *np.swapaxes(stub_xyzs, 0, 1),
            cen=list(stub_xyzs[:, 1, :]),
            )
            b_res_list = np.arange(1,plen+1)
            reshaped_rts = np.repeat(rts.values.reshape(rt_count,4,4),bb_stubs.shape[0] ,axis=0)

            reshaped_bb = np.tile( bb_stubs,(rt_count,1,1))

            new_rts = reshaped_bb @ reshaped_rts  # these are the ideal phosphate stubs
            inv_rts = np.linalg.inv(new_rts)
            bundle_p_res = np.tile(b_res_list, (inv_rts.shape[0]))
            bundle_bb_res = np.tile(b_res_list,(rts.rts.shape[0]*rts.pdb_path.shape[0]))
            bundle_bb_res = np.repeat(bundle_bb_res,bb_stubs.shape[0],axis=0)
            # bundle_p_res = np.tile(b_res_list,(inv_rts.shape[0]))
            # bundle_bb_res = np.tile(b_res_list,rt_count)
            # bundle_bb_res = np.repeat(bundle_bb_res,bb_stubs.shape[0],axis=0)
            inv_rts_repeat = np.repeat(inv_rts,bb_stubs.shape[0],axis=0)
            reshaped_bb_2 = np.tile(bb_stubs, (inv_rts.shape[0],1,1))
            rts_to_hash = inv_rts_repeat @ reshaped_bb_2
        #     print(rts_to_hash)
            keys = binner.get_bin_index(rts_to_hash)
        #     print(keys)
            mask = hashmap.contains(keys)
            q = keys[mask]
            values = hashmap[q]
            hit_paths = np.array(paths)[mask]
            hit_begin = np.array(start_res)[mask]
            hit_end = np.array(end_res)[mask]
            hit_p_res = np.array(p_res)[mask]
            bundle_bb_res = np.array(bundle_bb_res)[mask]
            bundle_p_res = np.array(bundle_p_res)[mask]

            with open (f"{out_dir_abspath}/grid_{i}.json","w") as f:
                json.dump(grid,f)

            pdb_path = f"{out_dir_abspath}/grid_{i}.pdb"
            p.dump_pdb(pdb_path)
            helix_list = list(range (1, grid["num_helices"]+1))
            grid_da = xr.DataArray(
                np.array(
                    [
                        list(grid[j].values())
                        for j in helix_list
                    ]
                )[np.newaxis,:],
                dims=("grid_id","helix","param"),
                coords={
                    "helix":helix_list,
                    "param":list(grid[1].keys())
                }
            )

            new_ds = xr.Dataset(
                {
                    "helix_params":grid_da,
                    "bundle_paths":("grid_id",np.array([pdb_path])),
                    "values":("paths",values),
                    "paths":("paths",hit_paths),
                    "begin":("paths",hit_begin),
                    "end":("paths",hit_end),
                    "p_res":("paths",hit_p_res),
                    "bundle_p_res":("paths",bundle_p_res),
                    "bundle_bb_res":("paths",bundle_bb_res)
                }
            )
            new_ds.to_netcdf(f"{out_dir_abspath}/hit_data_{i}.ncf")

if __name__ == "__main__":
    main()
