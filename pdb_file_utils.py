import os

def write_list_to_file (in_list, name, ext="",delim=" "):
    """
    Writes a list to the end of a file given (with chosen delim), adds the chosen extension
    """
    with open (name+ext,"w+") as f:
        f.write(delim.join([str(s) for s in in_list]))

def files_in_dir_by_ext(directory,ext):
    """
    returns a list of the full paths of the files with given extension
    
    ignores files beginning with "."
    """
    output = [
        root + f
        for root, dirs, fs in os.walk(directory, followlinks=True)
        for f in fs
        if f.split(".")[-1] == ext and f.split(".")[0]
    ]
    return output

def pdb_files_in_dir(directory):
    """
    returns a list of the full paths of the pdb files in the selected directory
    treats pdb.gz as a pdb
    """
    pdb_files = files_in_dir_by_ext(directory,"pdb")
    pdb_gz_files = [ f for f in files_in_dir_by_ext(directory,"gz") if f.split()[-2] == "pdb"]
    pdb_files.extend(pdb_gz_files)
    output = pdb_files
    return output