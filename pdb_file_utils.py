import os as _os


def write_list_to_file(in_list, name, ext="", delim=" "):
    """
    Writes a list to the end of a file given (with chosen delim), adds the chosen extension
    """
    with open(name + ext, "w+") as f:
        f.write(delim.join([str(s) for s in in_list]))


def files_in_dir_by_ext(directory, ext):
    """
    returns a list of the full paths of the files with given extension

    ignores files beginning with "."
    """
    output = [
        root + f
        for root, dirs, fs in _os.walk(directory, followlinks=True)
        for f in fs
        if f.split(".")[-1] == ext and f.split(".")[0]
    ]
    return output


def pdb_files_in_dir(directory):
    """
    returns a list of the full paths of the pdb files in the selected directory
    treats pdb.gz as a pdb
    """
    pdb_files = files_in_dir_by_ext(directory, "pdb")
    pdb_gz_files = [
        f
        for f in files_in_dir_by_ext(directory, "gz")
        if len(f.split(".")) > 1 and f.split(".")[-2] == "pdb"
    ]
    pdb_files.extend(pdb_gz_files)
    output = pdb_files
    return output


def strip_residues_from_pdb(
    pdb_file,
    first,
    last,
    chain="",
    prefix="stripped_",
    out_dir="",
    remark=True,
):
    """
    takes a pdb file and removes the given residue range

    outdir defaults to cwd, prefix defaults to stripped, use empty string for
    in-place
    Super Jenky Function Tm don't use this unless you have to
    """
    with open(pdb_file, "r") as pdb:
        lines = pdb.readlines()
        lines = [line.strip() for line in lines]
        trimmed = [
            line
            for line in lines
            if not (
                len(line) > 25
                and (chain in line[21])
                and (int(line[22:26].strip()) > first)
                and (int(line[22:26].strip()) < last)
                and ("ATOM" in line[:7])
            )
        ]
    out_dir = out_dir if out_dir else _os.getcwd()
    basename = _os.path.basename(pdb_file)
    with open(f"{out_dir}/{prefix}{basename}", "w+") as f:
        f.write(f"REMARK {pdb}" + "\n" + "\n".join(trimmed))
