import argparse
import subprocess
import shutil
from pathlib import Path

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

    parser.add_argument("-i", "--input", required=True, type=str, help="input PDB file")

    parser.add_argument(
        "-s", "--sequence", required=True, type=str, help="amino acid sequence"
    )


def run(args):
    # colabfold command check
    if not shutil.which("colabfold_batch"):
        LOGGER.error("colabfold_batch is not installed.")
        exit(1)

    # make input.fasta
    with open("input.fasta", "w") as f:
        f.write(">input\n")
        f.write(args.sequence)

    # make tmp directory for template input
    Path("tmp_template").mkdir(parents=True, exist_ok=True)
    cmd = f"cp {args.input} tmp_template/temp.pdb"
    subprocess.run(cmd, shell=True, check=True)
    LOGGER.info("tmp_template/temp.pdb copied")

    # run colabfold with template
    cmd = "colabfold_batch --custom-template-path tmp_template/ --num-models 1 --templates input.fasta results_modeled_cf --amber"
    subprocess.run(cmd, shell=True, check=True)
    LOGGER.info("results_modeled_cf/ generated")

    cmd = "rm -rf input.fasta tmp_template/"
    subprocess.run(cmd, shell=True, check=True)
    LOGGER.info("input.fasta and tmp_template/ removed")
