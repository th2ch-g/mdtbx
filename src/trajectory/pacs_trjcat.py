import argparse
import subprocess
from pathlib import Path

from ..logger import generate_logger

LOGGER = generate_logger(__name__)


DEFAULT_TOPOLOGY = "{args.trial_dir}/cycle000/replica001/rmmol_top.tpr"


def _list_subdirs(path, prefix):
    directories = sorted(
        child
        for child in Path(path).iterdir()
        if child.is_dir() and child.name.startswith(prefix)
    )
    if not directories:
        raise FileNotFoundError(f"No {prefix} directories found in {path}")
    return directories


def add_subcmd(subparsers):
    """
    mdtbx pacs_trjcat --trial_dir <int> --skip <int> --keep_selection <str> --centering_selection <str> --ref_structure <str>
    """
    parser = subparsers.add_parser(
        "pacs_trjcat",
        help="Concatenate PaCS-MD trajectories",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-t", "--trial_dir", required=True, type=str, help="Path to Trial directory"
    )

    parser.add_argument(
        "-r",
        "--ref_structure",
        default=DEFAULT_TOPOLOGY,
        type=str,
        help="Path to Reference structure (.tpr) for fitting (This could be also used as topology file)",
    )

    parser.add_argument(
        "-f",
        "--fit_selection",
        default="Protein",
        type=str,
        help="Selection for fitting (index group)",
    )

    parser.add_argument(
        "-s", "--skip", default=1, type=int, help="Number of frames to skip"
    )

    parser.add_argument(
        "-trj",
        "--trjname",
        default="prd_rmmol.xtc",
        type=str,
        help="Trajectory name in each replica",
    )

    parser.add_argument(
        "-k",
        "--keep_selection",
        default="System",
        type=str,
        help="Keep selection (index group)",
    )

    parser.add_argument(
        "-c",
        "--centering_selection",
        default="Protein",
        type=str,
        help="Centering selection (index group)",
    )

    parser.add_argument("-idx", "--index", type=str, help="Index file")

    parser.add_argument(
        "--pbc",
        default="mol",
        type=str,
        help="PBC option for gmx trjconv",
        choices=["none", "mol", "res", "atom", "nojump", "cluster", "whole"],
    )

    parser.add_argument(
        "--keep_cycle_trj",
        default=False,
        action="store_true",
        help="Keep cycle trajectory (e.g. trial001/cycle000/prd_all.xtc)",
    )

    parser.set_defaults(func=run)


def check_cycle(args):
    return len(_list_subdirs(args.trial_dir, "cycle"))


def check_replica(args):
    cycle_dirs = _list_subdirs(args.trial_dir, "cycle")
    return len(_list_subdirs(cycle_dirs[0], "replica"))


def run(args):
    # ref: https://zenn.dev/kh01734/articles/012380a58949d1
    cycle_dirs = _list_subdirs(args.trial_dir, "cycle")

    if args.ref_structure == DEFAULT_TOPOLOGY:
        args.ref_structure = str(cycle_dirs[0] / "replica001" / "rmmol_top.tpr")

    if args.index is not None:
        INDEX_OPTION = f"-n {args.index}"
    else:
        INDEX_OPTION = ""

    ext = Path(args.trjname).suffix

    n_cycle = len(cycle_dirs)
    n_replica = len(_list_subdirs(cycle_dirs[0], "replica"))

    # topology conversion
    cmd = f"echo {args.keep_selection} | gmx convert-tpr -s {args.ref_structure} {INDEX_OPTION} -o {args.trial_dir}/rmmol_top.tpr"
    subprocess.run(cmd, shell=True, check=True)
    LOGGER.info("rmmol_top.tpr generated")

    cmd = f"gmx editconf -f {args.trial_dir}/rmmol_top.tpr -o {args.trial_dir}/rmmol_top.gro"
    subprocess.run(cmd, shell=True, check=True)
    LOGGER.info("rmmol_top.gro generated")

    c_cmd = "c\n" * n_replica
    for cycle_index, cycle_dir in enumerate(cycle_dirs):
        # trjcat
        if cycle_index != 0:
            trj_files = [
                f"{cycle_dir}/replica{replica:03}/{args.trjname} "
                for replica in range(1, n_replica + 1)
            ]
            trj_files = " ".join(trj_files)
            cmd = f"echo '{c_cmd}' | gmx trjcat -f {trj_files} -o {cycle_dir}/tmp_all{ext} -settime"
            subprocess.run(cmd, shell=True, check=True)
            LOGGER.info(f"{cycle_dir}/tmp_all{ext} generated")
        else:
            trj_file = f"{cycle_dir}/replica001/{args.trjname}"
            cmd = f"cp {trj_file} {cycle_dir}/tmp_all{ext}"
            subprocess.run(cmd, shell=True, check=True)
            LOGGER.info(f"{cycle_dir}/tmp_all{ext} copied")

        # trjconv
        cmd = f"echo {args.centering_selection} System | gmx trjconv -f {cycle_dir}/tmp_all{ext} -s {args.ref_structure} {INDEX_OPTION} -o {cycle_dir}/tmp_all_pbc{ext} -center -pbc {args.pbc}"
        subprocess.run(cmd, shell=True, check=True)
        LOGGER.info(f"{cycle_dir}/prd_all{ext} generated")

        cmd = f"echo {args.fit_selection} {args.centering_selection} {args.keep_selection} | gmx trjconv -f {cycle_dir}/tmp_all_pbc{ext} -s {args.ref_structure} {INDEX_OPTION} -o {cycle_dir}/prd_all{ext} -center -fit rot+trans"
        subprocess.run(cmd, shell=True, check=True)
        LOGGER.info(f"{cycle_dir}/prd_all{ext} generated")

        # rm
        cmd = f"rm -f {cycle_dir}/tmp_all{ext} {cycle_dir}/\\#*"
        subprocess.run(cmd, shell=True, check=True)
        LOGGER.info(f"{cycle_dir}/tmp_all{ext} and backup files removed")

    c_cmd = "c\n" * n_cycle

    # trjcat
    trj_files = [
        f"{cycle_dir}/tmp_all_pbc{ext} "
        for cycle_dir in cycle_dirs
    ]
    trj_files = " ".join(trj_files)
    cmd = f"echo '{c_cmd}' | gmx trjcat -f {trj_files} -o {args.trial_dir}/tmp_all{ext} -settime"
    subprocess.run(cmd, shell=True, check=True)
    LOGGER.info(f"{args.trial_dir}/tmp_all{ext} generated")

    # trjconv
    cmd = f"echo {args.centering_selection} System | gmx trjconv -f {args.trial_dir}/tmp_all{ext} -s {args.ref_structure} {INDEX_OPTION} -o {args.trial_dir}/tmp_all_pbc{ext} -center -pbc {args.pbc}"
    subprocess.run(cmd, shell=True, check=True)
    LOGGER.info(f"{args.trial_dir}/tmp_all_pbc{ext} generated")

    cmd = f"echo {args.fit_selection} {args.centering_selection} {args.keep_selection} | gmx trjconv -f {args.trial_dir}/tmp_all_pbc{ext} -s {args.ref_structure} {INDEX_OPTION} -o {args.trial_dir}/prd_all{ext} -center -fit rot+trans"
    subprocess.run(cmd, shell=True, check=True)
    LOGGER.info(f"{args.trial_dir}/prd_all{ext} generated")

    # rm
    cmd = f"rm -f {args.trial_dir}/cycle*/tmp_all_pbc{ext} {args.trial_dir}/tmp_all{ext} {args.trial_dir}/tmp_all_pbc{ext} {args.trial_dir}/\\#*"
    subprocess.run(cmd, shell=True, check=True)
    LOGGER.info(f"{args.trial_dir}/tmp_all{ext} and backup files removed")

    if not args.keep_cycle_trj:
        cycle_outputs = " ".join(f"{cycle_dir}/prd_all{ext}" for cycle_dir in cycle_dirs)
        cmd = f"rm -f {cycle_outputs}"
        subprocess.run(cmd, shell=True, check=True)
        LOGGER.info("Per-cycle trajectories removed")
