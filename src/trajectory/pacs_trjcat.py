import argparse
from pathlib import Path

from ..logger import generate_logger
from ..utils.gmx import gmx_index_args, remove_gmx_backups
from ..utils.proc import run_cmd

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

    ext = Path(args.trjname).suffix
    if not ext:
        raise ValueError("--trjname must have a file extension")
    if args.skip < 1:
        raise ValueError("--skip must be positive")

    n_cycle = len(cycle_dirs)
    n_replica = len(_list_subdirs(cycle_dirs[0], "replica"))
    trial_dir = Path(args.trial_dir)
    output_topology = trial_dir / "rmmol_top.tpr"
    output_structure = trial_dir / "rmmol_top.gro"

    # topology conversion
    cmd = [
        "gmx",
        "convert-tpr",
        "-s",
        args.ref_structure,
        *gmx_index_args(args.index),
        "-o",
        str(output_topology),
    ]
    run_cmd(cmd, input=f"{args.keep_selection}\n", log="rmmol_top.tpr generated")

    cmd = [
        "gmx",
        "editconf",
        "-f",
        str(output_topology),
        "-o",
        str(output_structure),
    ]
    run_cmd(cmd, log="rmmol_top.gro generated")

    c_cmd = "c\n" * n_replica
    for cycle_index, cycle_dir in enumerate(cycle_dirs):
        combined_path = cycle_dir / f"tmp_all{ext}"
        pbc_path = cycle_dir / f"tmp_all_pbc{ext}"
        cycle_output = cycle_dir / f"prd_all{ext}"
        # trjcat
        if cycle_index != 0:
            trj_files = [
                str(cycle_dir / f"replica{replica:03}" / args.trjname)
                for replica in range(1, n_replica + 1)
            ]
            cmd = [
                "gmx",
                "trjcat",
                "-f",
                *trj_files,
                "-o",
                str(combined_path),
                "-settime",
            ]
            run_cmd(cmd, input=c_cmd, log=f"{combined_path} generated")
        else:
            trj_file = cycle_dir / "replica001" / args.trjname
            cmd = ["cp", str(trj_file), str(combined_path)]
            run_cmd(cmd, log=f"{combined_path} copied")

        # trjconv
        cmd = [
            "gmx",
            "trjconv",
            "-f",
            str(combined_path),
            "-s",
            args.ref_structure,
            *gmx_index_args(args.index),
            "-o",
            str(pbc_path),
            "-center",
            "-pbc",
            args.pbc,
            "-skip",
            str(args.skip),
        ]
        run_cmd(
            cmd,
            input=f"{args.centering_selection}\nSystem\n",
            log=f"{pbc_path} generated",
        )

        cmd = [
            "gmx",
            "trjconv",
            "-f",
            str(pbc_path),
            "-s",
            args.ref_structure,
            *gmx_index_args(args.index),
            "-o",
            str(cycle_output),
            "-center",
            "-fit",
            "rot+trans",
        ]
        run_cmd(
            cmd,
            input=(
                f"{args.fit_selection}\n"
                f"{args.centering_selection}\n"
                f"{args.keep_selection}\n"
            ),
            log=f"{cycle_output} generated",
        )

        # rm
        combined_path.unlink(missing_ok=True)
        remove_gmx_backups(cycle_dir)
        LOGGER.info(f"{combined_path} and GROMACS backup files removed")

    c_cmd = "c\n" * n_cycle
    combined_path = trial_dir / f"tmp_all{ext}"
    pbc_path = trial_dir / f"tmp_all_pbc{ext}"
    output_path = trial_dir / f"prd_all{ext}"

    # trjcat
    trj_files = [str(cycle_dir / f"tmp_all_pbc{ext}") for cycle_dir in cycle_dirs]
    cmd = [
        "gmx",
        "trjcat",
        "-f",
        *trj_files,
        "-o",
        str(combined_path),
        "-settime",
    ]
    run_cmd(cmd, input=c_cmd, log=f"{combined_path} generated")

    # trjconv
    cmd = [
        "gmx",
        "trjconv",
        "-f",
        str(combined_path),
        "-s",
        args.ref_structure,
        *gmx_index_args(args.index),
        "-o",
        str(pbc_path),
        "-center",
        "-pbc",
        args.pbc,
    ]
    run_cmd(
        cmd,
        input=f"{args.centering_selection}\nSystem\n",
        log=f"{pbc_path} generated",
    )

    cmd = [
        "gmx",
        "trjconv",
        "-f",
        str(pbc_path),
        "-s",
        args.ref_structure,
        *gmx_index_args(args.index),
        "-o",
        str(output_path),
        "-center",
        "-fit",
        "rot+trans",
    ]
    run_cmd(
        cmd,
        input=(
            f"{args.fit_selection}\n{args.centering_selection}\n{args.keep_selection}\n"
        ),
        log=f"{output_path} generated",
    )

    # rm
    for cycle_dir in cycle_dirs:
        (cycle_dir / f"tmp_all_pbc{ext}").unlink(missing_ok=True)
    combined_path.unlink(missing_ok=True)
    pbc_path.unlink(missing_ok=True)
    remove_gmx_backups(trial_dir)
    LOGGER.info(f"{combined_path} and GROMACS backup files removed")

    if not args.keep_cycle_trj:
        for cycle_dir in cycle_dirs:
            (cycle_dir / f"prd_all{ext}").unlink(missing_ok=True)
        LOGGER.info("Per-cycle trajectories removed")
