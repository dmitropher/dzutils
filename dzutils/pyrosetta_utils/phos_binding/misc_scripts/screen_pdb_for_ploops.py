import sys
import itertools as it

import pyrosetta

from dzutils.pyrosetta_utils import residues_with_element
from dzutils.pyrosetta_utils import hbond_to_residue
from dzutils.pyrosetta_utils import atom_indices_with_element
from dzutils.pyrosetta_utils import bonded_atoms
from dzutils.pyrosetta_utils.chain_utils import link_poses
from dzutils.sutils import read_flag_file

ploop_flags_file = "/home/dzorine/phos_binding/pilot_runs/loop_grafting/initial_testing/misc_files/p_ligand.flags"
flags = read_flag_file(ploop_flags_file)
flags_str = " ".join(flags.replace("\n", " ").split())
pyrosetta.init(flags_str)
# get pose
pose = pyrosetta.pose_from_pdb(sys.argv[1])
# get_outdir
outdir = sys.argv[2]
name = pose.pdb_info().name().split("/")[-1].split(".pdb")[0]
# Scan pose for phosphorus containing residues
# extract hbonds to these residues where:
# vec False to use list rather than rosetta vector
all_hbonds = {
    # index them by resnum of P atom
    resnum: {
        # Get all the hbonds
        "hbonds": hbond_to_residue(pose, resnum, vec=False),
        # And a dict of acceptable acceptor atoms (atoms bound to P)
        # keys are p atoms, values are lists of bound atoms
        "p_bonded_atoms": {
            atom_i: [
                atom
                for atom in bonded_atoms(
                    pose.residue(resnum), atom_i, name=False
                )
            ]
            for atom_i in atom_indices_with_element(pose.residue(resnum), "P")
        },
    }
    for resnum in residues_with_element(pose, "P")
}
bb_hbonds = {
    # Take resnums
    # Keep the donor resnums list
    r: [
        {atom_i:
            [b.don_res()
            for b in info["hbonds"]
            if b.acc_atm() in acceptor_atoms
            and b.don_hatm_is_protein_backbone()
            and b.don_res() != r]
        }
        for atom_i, acceptor_atoms in info["p_bonded_atoms"].items()
        # Keep donor resnums where the acceptor atom is bonded to Phosphorus
        # And the donor atom is a backbone atom
    ]
    # Derive the above info by scanning the P containing residues
    for r, info in all_hbonds.items()
    if info["hbonds"]
}

# no support for different append/prepend values
append_factor = 3  # +- 0-3 residues
append_ranges = [range(append_factor + 1)] * 2

num_contacts = int(sys.argv[3])

pose_size = len(pose.residues)
# print(bb_hbonds)
# Extract contiguous loops with these contacts:
#   - for each bb pair, check if intervening sequence is under 10 res
#   - prepend and apppend +- 3 residues on each side if they exist
print (bb_hbonds)
loop_pose_dicts = [
    {
        "pose": link_poses(
            pyrosetta.rosetta.protocols.grafting.return_region(
                pose.clone(), r, r
            ),
            pyrosetta.rosetta.protocols.grafting.return_region(
                pose, min(*contact_set) - x, max(*contact_set) + y
            ),
            rechain=True,
        ),
        "res": r,
        "start": min(*contact_set) - x,
        "end": max(*contact_set) + y,
    }
    for r, per_p_atom_contact_list in bb_hbonds.items()
    for atom_i_contacts in per_p_atom_contact_list
    for atom_i, contacts in atom_i_contacts.items()
    if len (contacts) >= num_contacts
    for contact_set in it.combinations(contacts, num_contacts)
    for x, y in it.product(*append_ranges)
    if max(*contact_set) - min(*contact_set) + abs(x + y) < 11
    if min(*contact_set) - x > 1 and max(*contact_set) + y < pose_size
    if min(*contact_set) - x > r or r > max(*contact_set) + y
]

for d in loop_pose_dicts:
    p = d["pose"]
    r = d["res"]
    s = d["start"]
    e = d["end"]
    p.dump_pdb(
        f"{outdir}/{name}_{num_contacts}-contacts_phos-{r}_frag_{s}-{e}.pdb"
    )
