import h5py
import click
from os import makedirs
from os.path import isdir


@click.command()
@click.argument("task_file_path")
@click.argument("hashmap_path", type=click.Path(exists=True))
@click.argument("hdf5_store", type=click.Path(exists=True))
@click.option("-a", "--angle-res", default=15.0)
@click.option("-d", "--angstrom-dist-res", default=1.0)
@click.option("-r", "--run-name", default="expand_rotamer_set")
@click.option("-e", "--erase/--no-erase", default=False)
@click.option("-c", "--chunk-size", default=150)
@click.option("-o", "--output-dir", default=".")
@click.option("-n", "--output-dir-names", default="out")
@click.option(
    "-p",
    "--expand-script",
    default="~/scripts/dzutils/dzutils/pyrosetta_utils/phos_binding/misc_scripts/expand_rotable_chis.py",
)
@click.option("-g", "--granularity-factor", default=0)
@click.option("-s", "--search-radius", default=0)
def main(
    task_file_path,
    hashmap_path,
    hdf5_store,
    run_name="inverse_ptr_exchi7_rotamers",
    angstrom_dist_res=1,
    angle_res=15,
    erase=False,
    chunk_size=150,
    output_dir=".",
    output_dir_names="",
    expand_script="~/scripts/dzutils/dzutils/pyrosetta_utils/phos_binding/misc_scripts/expand_rotable_chis.py",
    granularity_factor=0,
    search_radius=0,
):
    """
    """

    if not isdir(output_dir):
        makedirs(output_dir)
    with h5py.File(hdf5_store, "r") as f:
        num_chi_sets = f["chis"].shape[0]

    tasks = []
    for start in range(0, num_chi_sets, chunk_size):
        dirname = (
            f"{output_dir}/{output_dir_names}_{start}_{start+chunk_size-1}"
        )
        if not isdir(dirname):
            makedirs(dirname)
        run_cmd_string = f"""python  {
                    expand_script} {
                    hashmap_path} {
                    hdf5_store} -a {
                    angle_res} -d {
                    angstrom_dist_res } -r {run_name
                    } {'-e' if erase else ''} -n {
                    start}-{start+chunk_size-1} {
                    '-g ' + str(granularity_factor )if granularity_factor else ''} {
                    '-s '+str(search_radius) if search_radius else ''}"""
        with open(f"{dirname}/run.sh", "w") as f:
            f.write(run_cmd_string)
        tasks.append(f"cd {dirname} ; bash run.sh")

    with open(task_file_path, "w") as f:
        f.write("\n".join(tasks))


if __name__ == "__main__":
    main()
