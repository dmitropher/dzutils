# when run on cli, take a list of files and try to load them as pdb
# save a list of dicts for the successes
# send it to json, optionally write to disk

# Dict includes:
# pdb_info().name()
# len(pose.residues)
# pose.sequence()
# additional fields can be added by the format key=func
# Where the entry in the dict will be "key":func(pose).
# The exceptions to these functions are not handled in any way

# also records failures into a dict:
# pdb_info().name()
# error:exception
import click
import pyrosetta
import json

from dzutils.sutils import read_flag_file

# from dzutils.func_utils import HiddenPrints


def pose_info_dict(pose, **custom_data):
    """"
    Returns a dict with info about the given pose
    Dict includes:
    "name" : pdb_info().name()
    "length": len(pose.residues)
    "sequence": pose.sequence()
    additional fields can be added by the format key=func
    Where the entry in the dict will be "key":func(pose).
    The exceptions to these functions are not handled in any way
    """
    dict_ = {
        "name": pose.pdb_info().name(),
        "length": len(pose.residues),
        "sequence": pose.sequence(),
    }
    for k, f in custom_data:
        dict_[k] = f(pose)
    return dict_


def evaluate_pdbs(
    pdbs,
    rosetta_flags_file="",
    json_name="evaluated_pdbs",
    json_path=None,
    print_=False,
    **kwargs,
):
    """
    Takes a list of pdbs, attempts to generate poses. Returns a list of info

    Can take a path to write the list of dicts as json. For errored pdbs
    defaults to the name json_name_erred.json
    """
    loaded_dict_list = []
    errored_dict_list = []
    # # hack to suppress rosetta output, doesn't work. more complicated stream redirection necessary
    # with HiddenPrints():
    flags = read_flag_file(rosetta_flags_file) if rosetta_flags_file else ""
    pyrosetta.init(flags)

    for pdb in pdbs:
        try:
            pose = pyrosetta.pose_from_file(pdb)
            loaded_dict_list.append(pose_info_dict(pose, **kwargs))
        except Exception as e:
            errored_dict_list.append({"pdb": pdb, "exception": str(e)})

    if print_:
        click.echo("\n".join([d["name"] for d in loaded_dict_list]))

    if json_path:
        with open(f"{json_path}/{json_name}.json", "w") as f:
            json.dump(loaded_dict_list, f)
        with open(f"{json_path}/{json_name}_erred.json", "w") as f:
            json.dump(errored_dict_list, f)
    return loaded_dict_list, errored_dict_list


@click.command()
@click.option("-j", "--json-name", default="evaluated_pdbs")
@click.option("-o", "--json-path")
@click.option("-r", "--rosetta-flags-file")
@click.option("-p", "--print", "print_", is_flag=True)
@click.argument("pdbs", nargs=-1, type=click.Path(exists=True))
def main(
    pdbs,
    rosetta_flags_file="",
    json_name="evaluated_pdbs",
    json_path=None,
    print_=False,
    **kwargs,
):
    """
    CLI pdb evaluator, checks if they can be loaded as poses with the given flags
    """
    evaluate_pdbs(
        pdbs, rosetta_flags_file, json_name, json_path, print_, kwargs
    )


if __name__ == "__main__":
    main()
