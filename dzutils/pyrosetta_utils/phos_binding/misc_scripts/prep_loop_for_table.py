import sys

# import pandas as pd

import pyrosetta

# from dzutils.pyrosetta_utils.phos_binding import get_loop_xform_dicts
from dzutils.pyrosetta_utils import (
    residues_with_element,
    hbond_to_residue,
    atom_indices_with_element,
    bonded_atoms,
)
from dzutils.pyrosetta_utils.phos_binding import replace_p_res_with_phosphate
from dzutils.sutils import read_flag_file

# from dzutils.pdb_file_utils import pdb_files_in_dir


ploop_flags_file = "/home/dzorine/phos_binding/pilot_runs/loop_grafting/initial_testing/misc_files/p_ligand_debug.flags"
flags = read_flag_file(ploop_flags_file)
flags_str = " ".join(flags.replace("\n", " ").split())
pyrosetta.init(flags_str)


# pdb_dir = "/home/dzorine/phos_binding/ploop_pdbs/ploop_set_2/3_bb_contacts"
#
# pdb_store = (
#     "/home/dzorine/phos_binding/ploop_pdbs/ploop_set_2/3_bb_contacts_spinnable"
# )
# pdb_paths = pdb_files_in_dir(pdb_dir)[:100]
# # subset = pdb_paths
# poses = (maybe_load(p) for p in pdb_paths)
# dicts = list()
# for pose in poses:
pose = pyrosetta.pose_from_file(sys.argv[1])
out_path = sys.argv[2]
n_contacts = 3 if len(sys.argv) < 4 else int(sys.argv[3])
# replace p-res with phosphate(s)
print("REPLACING PRES")
phos_pose = replace_p_res_with_phosphate(pose, min_contacts=n_contacts)
phos_pose.dump_pdb("/home/dzorine/temp/phos_replaced.pdb")
# Check each one for minimum contacts, remove ones that don't matter
phosphates = residues_with_element(phos_pose, "P")
print(phosphates)
to_remove = [
    resnum
    for resnum in phosphates
    if len(
        [
            hbond
            for hbond in hbond_to_residue(phos_pose, resnum)
            if hbond.don_hatm_is_protein_backbone()
            # Unecessary check with single phosphate
            # and hbond.acc_atm()
            # in bonded_atoms(
            #     phos_pose.residue(resnum),
            #     atom_indices_with_element(phos_pose.residue(resnum), "P")[0],
            # )
        ]
    )
    < n_contacts
]

for i, resnum in enumerate(to_remove):
    phos_pose.delete_residue_slow(resnum - i)
# Check that chain 1 is the loop (contains all the contact making res)
phosphates = residues_with_element(phos_pose, "P")
contact_res = list(
    set(
        [
            hbond.don_res()
            for resnum in phosphates
            for hbond in hbond_to_residue(phos_pose, resnum)
            if hbond.acc_res() in phosphates
        ]
    )
)
contact_chains = list(set([phos_pose.chain(res) for res in contact_res]))
assert bool(len(contact_chains) == 1)
contact_chain = contact_chains[0]

# split each phosphate to a new pose
poses = list()
for pres in phosphates:
    newp = phos_pose.split_by_chain()[contact_chain]
    newp.append_residue_by_jump(phos_pose.residue(pres), 1)
    poses.append(newp)

# apply reslabels to spinnable res

for i, p in enumerate(poses, 1):
    phosphates = residues_with_element(p, "P")
    assert len(phosphates) == 1
    p_res = phosphates[0]
    contact_res = list(
        set(
            [
                hbond.don_res()
                for hbond in hbond_to_residue(p, p_res)
                if hbond.acc_res() == p_res
            ]
        )
    )
    n_term = p.chain_begin(1)
    c_term = p.chain_end(1)
    if n_term not in contact_res:
        p.pdb_info().add_reslabel(n_term, "n_spinnable")
    if c_term not in contact_res:
        p.pdb_info().add_reslabel(n_term, "c_spinnable")

    base_name_old = pose.pdb_info().name().split("/")[-1].split(".pdb")[0]
    new_name = f"{base_name_old}_p{i}_spinnable.pdb"
    new_out_path = f"{out_path}/{new_name}"
    phos_pose.dump_pdb(new_out_path)
# Check if the first and last res of chain 1 are free to spin

# if not, pre/append a residue, fix omegas
# xform_dicts = get_loop_xform_dicts(pose, n_contacts)
# if xform_dicts:
#     dicts.extend(xform_dicts)

# loop_table = pd.DataFrame(dicts)
# At this point table is made with things ready to spin

# begin dict filling loop
# Max cycles 10k, raise RuntimeError if exceeded