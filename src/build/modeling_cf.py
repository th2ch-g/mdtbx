import argparse
import shutil
from pathlib import Path

from ..utils.proc import run_cmd
from ..logger import generate_logger

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    """
    mdtbx modeling_cf -i <pdb> -s <sequence>
    """
    parser = subparsers.add_parser(
        "modeling_cf",
        help="modeling by ColabFold",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-i", "--input", type=str, help="input PDB file used as template"
    )

    parser.add_argument(
        "-s", "--sequence", required=True, type=str, help="amino acid sequence"
    )

    parser.set_defaults(func=run)


def run(args):
    # colabfold command check
    if not shutil.which("colabfold_batch"):
        LOGGER.error("colabfold_batch is not installed.")
        exit(1)

    # make input.fasta
    with open("input.fasta", "w") as f:
        f.write(">input\n")
        f.write(args.sequence)

    # run colabfold with template
    if args.input is None:
        cmd = "colabfold_batch --num-models 1 input.fasta results_modeled_cf --amber"
    else:
        # make tmp directory for template input
        Path("tmp_template").mkdir(parents=True, exist_ok=True)
        shutil.copy(args.input, "tmp_template/temp.pdb")
        LOGGER.info("tmp_template/temp.pdb copied")
        cmd = "colabfold_batch --custom-template-path tmp_template/ --num-models 1 --templates input.fasta results_modeled_cf --amber"

    run_cmd(cmd, log="results_modeled_cf/ generated")

    Path("input.fasta").unlink(missing_ok=True)
    shutil.rmtree("tmp_template", ignore_errors=True)
    LOGGER.info("input.fasta and tmp_template/ removed")
