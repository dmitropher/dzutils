"""Utils for rosetta projects and slurm usage on digs"""


def notebook_setup(init_flags=str):
    import pyrosetta

    flags = " ".join(init_flags.split())
    pyrosetta.init(flags)


def init_env(
    flags,
    workers=10,
    init_flags_script="/home/dzorine/scripts/misc_scripts/worker_setup.py",
    memory="4GB",
    queue="short",
    walltime="00:30:00",
    cores=1,
):
    """Function initializes current jupyter notebook environment with string `flags`,
    also initializing each of integer `workers` workers on a SLURM scheduled computer
    cluster with the same `flags`.
    """
    from dask_jobqueue import SLURMCluster
    from dask.distributed import Client, progress

    flags_str = " ".join(flags.replace("\n", " ").split())
    # flags_str = flags.replace("\n", " ")
    extra = f" --preload {init_flags_script} ' {flags_str}'"
    #%run {init_flags_script}
    notebook_setup(flags_str)

    cluster = SLURMCluster(
        cores=cores,
        processes=1,
        memory=memory,
        queue=queue,
        extra=extra,
        walltime=walltime,
        local_directory="/home/dzorine/jupyter_wd/dask-worker-space/",
    )
    print(cluster.job_script())
    cluster.scale(int(workers))
    client = Client(cluster)
    return client


def results_from_futures(futureList):
    """
    takes a list of futures submitted to a cluster, returns a list of results when they have all finished.
    """
    resultList = []
    for i, a in enumerate(futureList):
        if a.done() and a.status != "error":
            result = a.result()
            if result:
                resultList.append(result)
            del futureList[i]
        elif a.status == "error":
            del futureList[i]
        elif a.status != "pending" and a.status != "finished":
            print(f"share status is: a.status")
            del futureList[i]
    return resultList


def read_file(filename):
    with open(filename, "r") as myfile:
        data = myfile.read()
    return data


def read_file_lines(filename):
    with open(filename, "r") as myfile:
        lines = myfile.readlines()
    return lines


def to_pdb_file(packed_pose, filename):
    with open(filename, "w+") as opdb:
        opdb.write(io.to_pdbstring(packed_pose))


# TODO: DOES not yet handle comments after the options!
def read_flag_file(filename):
    """Reads the flag file, ignoring comments"""
    lines = read_file_lines(filename)
    # filter the lines
    lines = [l for l in lines if l.startswith("-")]
    return " ".join(lines)


def persist_me_dasky(listy, client):
    return [client.persist(j) for j in listy]
