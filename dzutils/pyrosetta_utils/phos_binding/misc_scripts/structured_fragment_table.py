from os.path import isfile

import pandas as pd
import getpy as gp
import numpy as np

import pyrosetta

from xbin import XformBinner as XB

from dzutils.pyrosetta_utils.geometry.pose_xforms import (
    get_e2e_xform,
    generate_pose_rt_between_res,
)
from dzutils.pyrosetta_utils import (
    residues_with_element,
    run_pyrosetta_with_flags,
)
from dzutils.pdb_file_utils import pdb_files_in_dir

from dzutils.pyrosetta_utils import (
    residues_with_element,
    hbond_to_residue,
    atom_indices_with_element,
    bonded_atoms,
    or_compose_residue_selectors,
    build_hbond_set,
    residues_by_name,
)
from dzutils.pyrosetta_utils.secstruct import parse_structure_from_dssp

# from dzutils.pyrosetta_utils.chain_utils import rechain_resname


# from dzutils.pyrosetta_utils.geometry.pose_xforms import generate_pose_rt_between_res,get_func_to_end


def maybe_load(pdb):
    """
    Dumb hack so rosettta can tolerate pdbs it can't load
    """
    try:
        pose = pyrosetta.pose_from_file(pdb)
        return pose
    except Exception as e:
        print(e)
        return pyrosetta.pose_from_sequence("A")


def get_xform_dicts(
    pose,
    xbin_args=[],
    xbin_kwargs={},
    by_chain=False,
    dssp_types="",
    **feature_params,
):
    """
    Returns a dict list with pose xform objects for various features of the structure

    Includes an end to end xform on bb atoms, end to end xform for each piece of
    secstruct, the type of xform, and an xform to any ligand

    Returned xform types are either e2e or secstruct_fragment

    by_chain may be specified as a chain number to compute the e2e for first and
    last residues from that chain instead of otherwise. Secondary structure will
    only be parsed from the chain specified

    feature_params may be left blank for no feature

    feature_params should include:

    "residue_name_3" - 3 letter residue code OR "resnum" - residue number

    "aligment_atoms" - four atoms to create a stub from (can be index or name)
        If left blank bb atoms will be used (fails if the feature has no bb atoms)

    """
    # {
    #     "file": pose.pdb_info().name(),
    #     "key_int": int(xb.get_bin_index(_np.array(e2e))),
    #     "e2e": _np.array(e2e),
    #     # Careful here: maybe add a check to see that res 1 is really the beginning of the loop chain?
    #     "func_to_bb_start": generate_pose_rt_between_res(
    #         pose, resnum, 1, (p_atom, p_atom, bonded_1, bonded_2)
    #     ),
    #     "alignment_atoms": (p_atom, p_atom, bonded_1, bonded_2),
    # }

    working = pose.clone()
    e2e_pose = working if not by_chain else working.split_by_chain(by_chain)

    structs = parse_structure_from_dssp(e2e_pose, *dssp_types)

    e2e_xform = generate_pose_rt_between_res(
        working.clone(),
        working.chain_begin(by_chain) if by_chain else 1,
        working.chain_end(by_chain),
    )
    struct_pose_xforms = [
        generate_pose_rt_between_res(
            working.clone(), struct.start_pos, struct.end_pos
        )
        for struct in structs
    ]
    all_xforms = [
        {
            "type": "e2e",
            "pose_xform": e2e_xform,
            "filename": e2e_xform.pose_start.pdb_info().name(),
            "start": e2e_xform.seqpos_start,
            "end": e2e_xform.seqpos_end,
        },
        *[
            {
                "type": "secstruct_fragment",
                "pose_xform": xform,
                "filename": xform.pose_start.pdb_info().name(),
                "start": xform.seqpos_start,
                "end": xform.seqpos_end,
            }
            for xform in struct_pose_xforms
        ],
    ]
    xb = XB(*xbin_args, **xbin_kwargs)
    if feature_params:
        if not (
            ("residue_name_3" in feature_params) ^ ("resnum" in feature_params)
        ):

            raise ValueError("Use either residue_name_3, resnum but not both")
        resnum = (
            feature_params["resnum"]
            if ("resnum" in feature_params)
            else residues_by_name(working, feature_params["residue_name_3"])[0]
        )

        all_xforms = [
            {
                **dict_,
                **{
                    "feature_to_start": np.array(
                        generate_pose_rt_between_res(
                            working.clone(),
                            resnum,
                            dict_["pose_xform"].seqpos_start,
                            feature_params["alignment_atoms"],
                        )
                    ),
                    "alignment_atoms": feature_params["alignment_atoms"],
                    "key_int": xb.get_bin_index(np.array(dict_["pose_xform"])),
                },
            }
            for dict_ in all_xforms
        ]
    return all_xforms


flags_path = "/home/dzorine/phos_binding/pilot_runs/loop_grafting/initial_testing/misc_files/p_ligand_quiet.flags"
run_pyrosetta_with_flags(flags_path)


pdb_dir = "/mnt/home/dzorine/projects/phos_binding/ploop_pdbs/ploop_set_4_fragments/single_phos_3_bb_contacts"
pdb_paths = pdb_files_in_dir(pdb_dir)
# subset = pdb_paths
poses = (maybe_load(p) for p in pdb_paths)
dicts = list()
for pose in poses:
    names = list(
        set(
            [pose.residue(i).name3() for i in residues_with_element(pose, "P")]
        )
    )
    newp = pose.clone()
    if len(newp.residues) < 3:
        continue
    # for name in names:
    #     newp = rechain_resname(newp, name)
    xform_dicts = get_xform_dicts(
        newp,
        dssp_types="H",
        by_chain=1,
        alignment_atoms=("P1", "P1", "O1", "O3"),
        residue_name_3="PO4",
    )
    if xform_dicts:
        dicts.extend(xform_dicts)

xform_table = pd.DataFrame(dicts)
xform_table["index"] = xform_table.index
print(len(xform_table.index))
data_name = "1ang_3_contact_ploop_set_4_v4"
data_store_path = "/home/dzorine/phos_binding/pilot_runs/loop_grafting/fragment_tables/ploop_fragment_set"
table_out_path = f"{data_store_path}/tables/{data_name}.json"
assert bool(
    not (isfile(table_out_path))
), "Table with this name already exists"
xform_table.to_json(table_out_path)
# loop_table = pd.read_json(f"{data_store_path}/tables/{data_name}.json")
e2e_keys = np.array(xform_table["key_int"], dtype=np.int64)
e2e_vals = np.array(xform_table["index"], dtype=np.int64)
key_type = np.dtype("i8")
value_type = np.dtype("i8")
gp_dict = gp.Dict(key_type, value_type)
gp_dict[e2e_keys] = e2e_vals
dict_out_path = f"{data_store_path}/dicts/{data_name}.bin"
assert bool(not (isfile(dict_out_path))), "Dict with this name already exists"
gp_dict.dump(dict_out_path)
