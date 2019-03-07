# This is a module for working with HBNet lines in pdb files, until I find a rosetta-ey way to get that info, or implement it if it hasn't been.
import pyrosetta as _pyrosetta

_pyrosetta.distributed.maybe_init()


def pdb_with_hbnets(pdb, pose=False, resType=False):
    if has_hbnet_res(pdb, resType=resType):
        print(f"""pdb {pdb}""")
        return pdb


def cluster_sort_hbnet_res(pdbList, resList, cluster):
    """
    Returns a dictionary with:
    keys -- residues in the resList
    values -- list of pdbs drawn from pdbList that have that residues listed in HBNet remarks

    Requires an initialized dask managed slurm cluster to use
    """

    futureDict = {}
    for pdb in pdbList:

        for r in resList:
            if r in futureDict:
                futureDict[r].append(
                    cluster.submit(pdb_with_hbnets, pdb, resType=r)
                )
            else:
                futureDict[r] = [
                    cluster.submit(pdb_with_hbnets, pdb, resType=r)
                ]

    sortedRes = {}
    for res in futureDict:
        sortedRes[res] = []
    for res, resFutures in futureDict.items():
        while resFutures:
            doneFutures = [future for future in resFutures if future.done()]
            if doneFutures:
                gathered = cluster.gather(doneFutures, errors="raise")
                for result in gathered:
                    if result and result not in sortedRes[res]:
                        sortedRes[res].append(result)
                for result in doneFutures:
                    del resFutures[resFutures.index(result)]

    return sortedRes


def sort_by_hbnet_resType(pdbList, resList, cluster=False):
    """
    Returns a dictionary with:
    keys -- residues in the resList
    values -- list of pdbs drawn from pdbList that have that residues listed in HBNet remarks

    optional: pass this function a dask managed slurm cluster to sort a very large number of files
    """

    if not cluster:
        sortedRes = {}
        for pdb in pdbList:
            pos = _pyrosetta.pose_from_file(pdb)
            for r in resList:
                if has_hbnet_res(pdb, pos, r):
                    if r in sortedRes and pdb not in sortedRes[r]:
                        sortedRes[r].append(pdb)
                    else:
                        sortedRes[r] = [pdb]
        return sortedRes

    else:
        return cluster_sort_hbnet_res(pdbList, resList, cluster)


def parse_hbnets_from_pdb(pdb):
    """
    Returns a list of residue numbers parsed from the hbnet lines of a pdb, or returns a blank list if none are found
    """

    f = open(pdb, "r")
    rl = f.readlines()
    hbnets = []

    # from lines starting with "REMARK" and containing "HBNet" saves the second to last space delimited value (should be the HBNet residue)
    hbnets = [
        int(line.split(" ")[-2])
        for line in rl
        if "REMARK" in line[0:6] and "HBNet" in line
    ]

    return hbnets


def hbnet_restypes(pdb, resTypes=[]):
    """
    Takes a pdb file location and returns a list of tuples with residue number and type.

    Tuples are formatted: [POSE_NUMBER,3_LETTER_CODE] where POSE_NUMBER is the residue number in the HBNet REMARK line and 3_LETTER_CODE is what a protein scientist would assume it to be. It returns an empty list if none are found.

    Keyword arguments:
    resTypes -- restricts the output to your chosen residue type. Takes a list of three letter AA codes.

    """

    p = _pyrosetta.pose_from_pdb(pdb)
    hbnets = parse_hbnets_from_pdb(pdb)
    outputList = []

    if not hbnets:
        return outputList
    if not resTypes:
        for res in hbnets:
            resType = p.residue(res).name()
            resPosition = [res, resType]
            # print(str(resPosition) + ", " + resType)
            outputList.append(resPosition)
    else:
        # formats the resTypes list so it has no duplicates and is all uppercase
        resTypes = list(set(resTypes))
        resTypes = [t.upper() for t in resTypes]
        for res in hbnets:
            rt = p.residue(res).name()
            if rt in resTypes:
                resPosition = [res, rt]
                # print(resPosition)
                outputList.append(resPosition)
        p.empty()
    return outputList


def has_hbnet_res(pdb, pose=False, resType=False):
    """takes a residue type (3 letter code) and pdb file, returns True if any of the residues listed in HBNet REMARK lines are of that residue type"""
    if not resType:
        return False
    if not pose:
        pose = _pyrosetta.pose_from_file(pdb)
    resType = resType.upper()
    hbnets = parse_hbnets_from_pdb(pdb)

    # returns False if no REMARK lines describing hbnets are found
    if not hbnets:
        return False
    # checks if the residues listed in hbnets are of the residue type needed
    for res in hbnets:
        if resType == pose.residue(res).name():
            return True
    return False
