# Attempts to close from each chain to each other chain
# For each chain pair, does serial deletions from the ends to increase
# the number of samples

import itertools as it
import os

# import gc
import sys

import dask.bag as db
from dask.distributed import Client

import pyrosetta

from dzutils.sutils import read_flag_file
from dzutils.pyrosetta_utils.chain_utils import run_direct_segment_lookup

# from dzutils.pyrosetta_utils.chain_utils import serial_deletions
from dzutils.pyrosetta_utils.chain_utils import link_poses
from dzutils.pyrosetta_utils.chain_utils import trim_pose_to_term
from dzutils.pyrosetta_utils.chain_utils import pose_excluding_chain
from dzutils.pyrosetta_utils.chain_utils import replace_chain_by_number


def closed_loop_generator(pose, *args, **kwargs):

    working_pose = pose.clone()
    print("running direct segment lookup")
    segment_lookup_mover = run_direct_segment_lookup(
        working_pose, length=8, cluster_tol=0.5, rmsd_tol=0.5
    )
    print("mover complete")

    # Remember to check for NULL !!!
    status = segment_lookup_mover.get_last_move_status()
    if status != pyrosetta.rosetta.protocols.moves.mstype_from_name(
        "MS_SUCCESS"
    ):
        print("mover status: ", status)
        return
    yield working_pose
    for working_pose in iter(segment_lookup_mover.get_additional_output, None):
        print("popped pose from mover")
        print(working_pose)
        yield working_pose


def default_loop_close(pose, name, outdir):
    """
    function wrapper to reloop poses
    """
    working_pose = pose.clone()
    #     print(working_pose.conformation().num_chains())
    segment_lookup_mover = run_direct_segment_lookup(
        working_pose, length=8, cluster_tol=0.5, rmsd_tol=0.5
    )
    if working_pose:
        working_pose.dump_pdb(f"{outdir}/{name}_loop_1.pdb")
        more_output = iter(segment_lookup_mover.get_additional_output, None)
        enum_output = enumerate(more_output, 2)
        for i, more in enum_output:
            more.dump_pdb(f"{outdir}/{name}_loop_{i}.pdb")


def trim_chain(pose, chain_num, cut_amount, terminus=None):
    if terminus == "chain_end":
        return trim_pose_to_term(
            pose.split_by_chain()[chain_num],
            len(pose.split_by_chain()[chain_num].residues) - cut_amount,
            terminus=terminus,
        )
    if terminus == "chain_begin":
        return trim_pose_to_term(
            pose.split_by_chain()[chain_num], 1 + cut_amount, terminus=terminus
        )


def generate_name(pose, chain_begin, chain_end, cut_begin, cut_end):
    return f"""closed_{os.path.basename(
            pose.pdb_info().name()).split(".pdb")[0]
            }_chains_{chain_begin
            }-{chain_end
            }_cut_{cut_begin
            }-{cut_end
            }"""


def reorder_chains_for_closure(pose, chain_begin_num, chain_end_num):
    """
    returns pose with chain_begin as chain 1, and end as chain 2

    was necessary at some point for the loop closure method

    With only two chains it can bug out, hence the weird ternary operator
    """
    print(f"reordering chains: begin {chain_begin_num} end {chain_end_num}")
    return (
        link_poses(
            pose.split_by_chain()[chain_begin_num],
            pose.split_by_chain()[chain_end_num],
            pose_excluding_chain(pose, chain_begin_num, chain_end_num),
            rechain=True,
        )
        if pose.num_chains() > 2
        else link_poses(
            pose.split_by_chain()[chain_begin_num],
            pose.split_by_chain()[chain_end_num],
            rechain=True,
        )
    )


def exhaustive_single_loop_insertion(pose, deletion_amount, *args, **kwargs):
    """
    Takes arguments for loop closure method
    """

    num_chains = pose.num_chains()
    chain_pair_iter = it.permutations(range(1, num_chains + 1), 2)
    dels = (range(0, deletion_amount),) * 2
    for chain_begin_num, chain_end_num in chain_pair_iter:
        for cut_begin, cut_end in it.product(*dels):
            print("attempting cut begin/end: ", cut_begin, cut_end)
            print("chains begin/end: ", chain_begin_num, chain_end_num)
            if cut_end >= len(
                pose.split_by_chain()[chain_end_num].residues
            ) or cut_begin >= len(
                pose.split_by_chain()[chain_begin_num].residues
            ):
                print("not computing begin/end: ", cut_begin, cut_end)
                continue
            chain_begin = trim_chain(
                pose, chain_begin_num, cut_begin, "chain_end"
            )
            chain_end = trim_chain(pose, chain_end_num, cut_end, "chain_begin")
            pose = replace_chain_by_number(
                replace_chain_by_number(
                    pose.clone(), chain_begin, chain_begin_num
                ),
                chain_end,
                chain_end_num,
            )
            print("chains trimmed")
            reordered_pose = reorder_chains_for_closure(
                pose, chain_begin_num, chain_end_num
            )
            print("chains reordered")
            for out_pose in closed_loop_generator(
                reordered_pose, *args, **kwargs
            ):
                print("pose reloop successful")
                print(out_pose)
                yield out_pose


def check_for_client():
    try:
        client = Client(
            scheduler_file="/home/dzorine/jupyter_wd/scheduler.json"
        )
        return client

    except OSError as e:
        print(e)
        return
    except ConnectionRefusedError as e:
        print(e)
        return
    except TimeoutError as e:
        print(e)
        return
    except ValueError as e:
        print(e)
        return
    except RuntimeError as e:
        print(e)
        return
    except FileNotFoundError as e:
        print(e)
        return
    except IOError as e:
        print(e)
        return


def run_with_dask(pose, dels, chain_pair_iter, outdir, client):

    bag = db.from_sequence(
        (
            {"pose": pose.clone(), "chain_begin_num": i, "chain_end_num": j}
            for i, j in chain_pair_iter
        ),
        partition_size=1,
    )
    print(bag)
    bag_fut = client.persist(bag)
    print(bag_fut)
    cut_pair_bag_fut = client.persist(
        bag_fut.map(
            lambda dict_, dels=dels: [
                {
                    "pose": dict_["pose"].clone(),
                    "chain_begin_num": dict_["chain_begin_num"],
                    "chain_end_num": dict_["chain_end_num"],
                    "cut_begin": n,
                    "cut_end": m,
                }
                for n, m in it.product(*dels)
            ]
        )
    )
    print(cut_pair_bag_fut)
    # cut_pair_bag_fut.compute()
    flattened_fut = client.persist(cut_pair_bag_fut.flatten())
    print(flattened_fut)
    cut_begin_chain_fut = client.persist(
        flattened_fut.map(
            lambda dict_: {
                "pose": dict_["pose"].clone(),
                "chain_begin": trim_chain(
                    dict_["pose"],
                    dict_["chain_begin_num"],
                    dict_["cut_begin"],
                    terminus="chain_end",
                ).clone(),
                "chain_begin_num": dict_["chain_begin_num"],
                "chain_end_num": dict_["chain_end_num"],
                "cut_begin": dict_["cut_begin"],
                "cut_end": dict_["cut_end"],
            }
        )
    )
    print(cut_begin_chain_fut)
    # cut_begin_chain_fut.compute()
    cut_end_chain_fut = client.persist(
        cut_begin_chain_fut.map(
            lambda dict_: {
                "pose": dict_["pose"].clone(),
                "chain_end": trim_chain(
                    dict_["pose"],
                    dict_["chain_end_num"],
                    dict_["cut_end"],
                    terminus="chain_begin",
                ).clone(),
                "chain_begin": dict_["chain_begin"],
                "chain_begin_num": dict_["chain_begin_num"],
                "chain_end_num": dict_["chain_end_num"],
                "cut_begin": dict_["cut_begin"],
                "cut_end": dict_["cut_end"],
            }
        )
    )
    print(cut_end_chain_fut)
    # cut_end_chain_fut.compute()
    combined_cut_chain_fut = client.persist(
        cut_end_chain_fut.map(
            lambda dict_: {
                "pose": replace_chain_by_number(
                    replace_chain_by_number(
                        dict_["pose"],
                        dict_["chain_begin"],
                        dict_["chain_begin_num"],
                    ),
                    dict_["chain_end"],
                    dict_["chain_end_num"],
                ).clone(),
                "name": generate_name(
                    dict_["pose"],
                    dict_["chain_begin_num"],
                    dict_["chain_end_num"],
                    dict_["cut_begin"],
                    dict_["cut_end"],
                ),
                "chain_begin_num": dict_["chain_begin_num"],
                "chain_end_num": dict_["chain_end_num"],
            }
        )
    )
    print(combined_cut_chain_fut)
    # combined_cut_chain_fut.compute()
    loop_closed_fut = client.persist(
        combined_cut_chain_fut.map(
            lambda dict_, outdir=outdir: default_loop_close(
                link_poses(
                    *(
                        dict_["pose"].split_by_chain()[
                            dict_["chain_begin_num"]
                        ],
                        dict_["pose"].split_by_chain()[dict_["chain_end_num"]],
                    ),
                    pose_excluding_chain(
                        dict_["pose"],
                        dict_["chain_begin_num"],
                        dict_["chain_end_num"],
                    ),
                    rechain=True,
                ),
                dict_["name"],
                outdir,
            )
        )
    )
    print(loop_closed_fut)
    loop_closed_fut.compute()


def main():
    flagsFile = "/home/dzorine/phos_binding/pilot_runs/loop_grafting/initial_testing/misc_files/cluster_altered.flags"
    flags = read_flag_file(flagsFile)
    flags_str = " ".join(flags.replace("\n", " ").split())
    pyrosetta.init(flags_str)

    pose = pyrosetta.pose_from_file(sys.argv[1])
    outdir = sys.argv[2]
    num_chains = pose.num_chains()
    chain_pair_iter = it.permutations(range(1, num_chains + 1), 2)
    dels = (range(0, 5),) * 2
    client = check_for_client()
    if client:
        run_with_dask(pose, dels, chain_pair_iter, outdir, client)
    else:
        print("Dask client not found, running serially")
        for chain_begin_num, chain_end_num in chain_pair_iter:
            for cut_begin, cut_end in it.product(*dels):
                print("attempting cut begin/end: ", cut_begin, cut_end)
                print("chains begin/end: ", chain_begin_num, chain_end_num)
                if cut_end > len(
                    pose.split_by_chain()[chain_begin_num].residues
                ) or cut_begin > len(
                    pose.split_by_chain()[chain_end_num].residues
                ):
                    print("not computing begin/end: ", cut_begin, cut_end)
                    continue
                chain_begin = trim_chain(
                    pose, chain_begin_num, cut_begin, "chain_end"
                )
                chain_end = trim_chain(
                    pose, chain_end_num, cut_end, "chain_begin"
                )
                pose = replace_chain_by_number(
                    replace_chain_by_number(
                        pose.clone(), chain_begin, chain_begin_num
                    ),
                    chain_end,
                    chain_end_num,
                )
                name = generate_name(
                    pose, chain_begin_num, chain_end_num, cut_begin, cut_end
                )

                reordered_pose = reorder_chains_for_closure(
                    pose, chain_begin_num, chain_end_num
                )
                #     link_poses(
                #         *(
                #             pose.split_by_chain()[chain_begin_num],
                #             pose.split_by_chain()[chain_end_num],
                #         ),
                #         pose_excluding_chain(
                #             pose, chain_begin_num, chain_end_num
                #         ),
                #         rechain=True,
                #     )
                #     if pose.num_chains() > 2
                #     else link_poses(
                #         *(
                #             pose.split_by_chain()[chain_begin_num],
                #             pose.split_by_chain()[chain_end_num],
                #         ),
                #         rechain=True,
                #     )
                # )

                default_loop_close(reordered_pose, name, outdir)


if __name__ == "__main__":
    main()
