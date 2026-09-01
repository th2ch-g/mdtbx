import argparse
import shutil
import tempfile
from pathlib import Path

from ..logger import generate_logger
from ..utils.proc import run_cmd

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

    parser.add_argument(
        "-o",
        "--output",
        default="results_modeled_cf",
        type=str,
        help="Output directory",
    )

    parser.set_defaults(func=run)


def run(args):
    if not shutil.which("colabfold_batch"):
        raise RuntimeError("colabfold_batch is not installed")

    sequence = args.sequence.strip()
    if not sequence:
        raise ValueError("--sequence must not be empty")

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix=".mdtbx_colabfold_", dir=".") as tmp:
        temp_path = Path(tmp)
        fasta_path = temp_path / "input.fasta"
        fasta_path.write_text(f">input\n{sequence}\n")

        command = [
            "colabfold_batch",
            "--num-models",
            "1",
        ]
        if args.input is not None:
            input_path = Path(args.input)
            if not input_path.is_file():
                raise FileNotFoundError(f"Template PDB not found: {input_path}")
            template_dir = temp_path / "templates"
            template_dir.mkdir()
            # ColabFold treats the filename stem as a four-character template ID.
            shutil.copy2(input_path, template_dir / "tmpl.pdb")
            command.extend(
                [
                    "--custom-template-path",
                    str(template_dir),
                    "--templates",
                ]
            )
        command.extend([str(fasta_path), str(output_path), "--amber"])
        run_cmd(command, log=f"{output_path}/ generated")
        if not list(output_path.glob("*rank_001*.pdb")):
            raise RuntimeError("ColabFold did not generate a rank-001 model")
