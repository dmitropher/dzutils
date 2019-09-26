import json
import pyrosetta

from dzutils.pyrosetta_utils.secstruct.structure import (
    parse_structure_from_dssp,
)


def load_pdbs(pdb_list):
    """
    """
    for path in pdb_list:
        try:
            yield pyrosetta.pose_from_file(path)
        except Exception as e:
            print(e)
            print(f"Unable to load pose from: {path}")


def load_fragment_store_from_path(path):
    """
    Expects a json list of paths
    """
    store = []
    with open(path, "r") as f:
        store = json.load(f)
    loaded_store = [*load_pdbs(store)]
    if not loaded_store:
        raise RuntimeError("Unable to load any pdbs from fragment store!")
    return loaded_store


def graft_generator(pose, fragments):
    """
    """


def main(
    pose,
    fragment_store_path,
    inv_rot_table_path,
    inv_rot_dict_path,
    log_dir,
    dssp_match_types="",
    allowed_positions=False,
):
    """
    This program takes a pose and a fragment store and returns alignment graphs

    There are some rules about how stuff is lined up and dssp types yada yada
    """

    sec_structs = parse_structure_from_dssp(pose, dssp_types=dssp_match_types)
    fragments = load_fragment_store_from_path(fragment_store_path)
    grafts = [*graft_generator(pose, fragments)]
    report_grafts(grafts, log_dir)
    clash_checked_grafts = [*clash_checker(grafts)]
    report_grafts(clash_checked_grafts, log_dir)
    inv_rot_scanned = [
        *inv_rot_scanner(
            clash_checked_grafts, inv_rot_table_path, inv_rot_dict_path
        )
    ]
    dump_results(inv_rot_scanned)


if __name__ == "__main__":
    main()
