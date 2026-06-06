import argparse
from pathlib import Path

from ..config import MAXWARN
from ..utils.proc import run_cmd
from ..logger import generate_logger

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    """
    mdtbx centering_gro
    """
    parser = subparsers.add_parser(
        "centering_gro",
        help="Centering gro file",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-f", "--structure", required=True, type=str, help="Structure file (.gro .tpr)"
    )

    parser.add_argument(
        "-p", "--topology", required=True, type=str, help="Topology file (.top)"
    )

    parser.add_argument(
        "-c",
        "--centering_selection",
        default="Protein",
        type=str,
        help="Centering selection",
    )

    parser.add_argument("-g", "--gmx", default="gmx", type=str, help="gmx command")

    parser.add_argument(
        "-idx", "--index", default="index.ndx", type=str, help="Index file"
    )

    parser.add_argument(
        "-o",
        "--output",
        default="gmx.gro",
        type=str,
        help="Output file name",
    )

    parser.add_argument(
        "--no-editconf", action="store_true", help="Do not run gmx editconf"
    )

    parser.set_defaults(func=run)


def run(args):
    dummy_mdp_path = Path(__file__).parent / "dummy.mdp"
    cmd = f"{args.gmx} grompp -f {dummy_mdp_path} -c {args.structure} -r {args.structure} -p {args.topology} -n {args.index} -maxwarn {MAXWARN} -o tmp.tpr"  # NOQA
    run_cmd(cmd, log="tmp.tpr generated")

    cmd = f"echo {args.centering_selection} System | {args.gmx} trjconv -f {args.structure} -s tmp.tpr -n {args.index} -o {args.output} -pbc mol -center"
    run_cmd(cmd, log=f"{args.output} generated")

    cmd = "rm -f tmp.tpr"
    run_cmd(cmd, log="tmp.tpr removed")

    if not args.no_editconf:
        cmd = f"{args.gmx} editconf -f {args.output} -o {args.output} -resnr 1"
        run_cmd(
            cmd,
            log=f"{args.gmx} editconf -f {args.output} -o {args.output} -resnr 1 run",
        )
