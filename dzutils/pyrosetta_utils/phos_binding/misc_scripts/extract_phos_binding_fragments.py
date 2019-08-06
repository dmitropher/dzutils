import click

from pyrosetta import pose_from_file
from pyrosetta.rosetta.protocols.grafting import return_region

from dzutils.pyrosetta_utils import run_pyrosetta_with_flags
from dzutils.pyrosetta_utils.phos_binding import minimal_fragments_by_secondary_structure
from dzutils.pyrosetta_utils.chain_utils import link_poses


@click.command()
@click.option("-o", "--outdir")
@click.option("-c", "--num-contacts",default=3)
@click.option("-a", "--append-factor",default=3)
@click.option("-p", "--proximity",default=5)
@click.option("-l", "--lazy", default=False)
@click.option("-s", "--struct-types")
@click.option("-f", "--flags-file",default="/home/dzorine/phos_binding/pilot_runs/loop_grafting/initial_testing/misc_files/p_ligand_quiet.flags")
@click.argument("pdb", nargs=1, type=click.Path(exists=True))
def main(
	pdb,
	outdir,
	struct_types="",
	num_contacts=3,
	append_factor=3,
	proximity=5,
	lazy=False,
	flags_file=""
	):
	"""
	"""
	run_pyrosetta_with_flags(flags_file)
	# get pose
	pose = pose_from_file(pdb)
	# get_outdir
	name = pose.pdb_info().name().split("/")[-1].split(".pdb")[0]

	fragments = list()

	try:

		fragments = minimal_fragments_by_secondary_structure(
			pose,
			*struct_types,
			min_contacts=num_contacts,
			proximity=proximity,
			lazy=lazy,
			append_factor=append_factor
			)
	except AssertionError as e:
		print (e)
		print ("No appropriate fragments found")
		return



	for d in fragments:

	        r = d["acceptor_res"]
	        s = d["start"]
	        e = d["end"]
	        p = link_poses(
	            return_region(pose.clone(), r, r),
	            return_region(pose, s, e),
	            rechain=True,
	        )
	        p.dump_pdb(
	            f"{outdir}/{name}_{num_contacts}-contacts_phos-{r}_frag_{s}-{e}.pdb"
	        )

if __name__ == "__main__":
	main()