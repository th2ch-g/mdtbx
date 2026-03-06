import argparse
import subprocess
from pathlib import Path

import numpy as np

from ..logger import generate_logger

LOGGER = generate_logger(__name__)

# rism3d.snglpnt 実行テンプレート
RISM_INPUT_TEMPLATE = """\
&cntrl
  ntx=1, nstlim=0, dt=0.001,
  ntpr=1, ntwx=0,
  cut=12.0,
  irism=1,
/
&rism
  closure='kh',
  grdspc=0.5,0.5,0.5,
  solvbox=-1,-1,-1,
  buffer=14.0,
  ng3=-1,-1,-1,
  solvcut=14.0,
  tolerance=1e-5,
  npropagate=5,
  mdiis_del=0.7,
  mdiis_nvec=5,
  polardecomp=1,
  entropicdecomp=1,
  progress=1,
  verbose=1,
  ntwrism=1,
  write_thermo=1,
/
"""


def add_subcmd(subparsers):
    """
    mdtbx place_solvent -p leap.parm7 -x leap.rst7 -o solvent_placed.pdb
    """
    parser = subparsers.add_parser(
        "place_solvent",
        help="Place solvent molecules using 3D-RISM",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-p",
        "--prmtop",
        required=True,
        type=str,
        help="AMBER topology file (.parm7)",
    )

    parser.add_argument(
        "-x",
        "--coord",
        required=True,
        type=str,
        help="AMBER coordinate file (.rst7)",
    )

    parser.add_argument(
        "-o",
        "--output",
        default="solvent_placed.pdb",
        type=str,
        help="Output PDB file",
    )

    parser.add_argument(
        "--solvent",
        default="water",
        choices=["water"],
        help="Solvent type",
    )

    parser.add_argument(
        "--temperature",
        default=300.0,
        type=float,
        help="Temperature [K]",
    )

    parser.add_argument(
        "--threshold",
        default=2.0,
        type=float,
        help="Density threshold for peak extraction (g(r) value)",
    )

    parser.add_argument(
        "--keepfiles",
        action="store_true",
        help="Keep intermediate files",
    )

    parser.set_defaults(func=run)


def _parse_dx(dx_path):
    """OpenDX形式の密度グリッドをndarrayとして読み込む。

    Returns:
        data: (nx, ny, nz) ndarray
        origin: (3,) origin coordinates [Å]
        delta: (3, 3) grid spacing matrix
    """
    origin = None
    delta = []
    counts = None
    data_values = []

    with open(dx_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("#"):
                continue
            if line.startswith("object 1"):
                # object 1 class gridpositions counts nx ny nz
                parts = line.split()
                counts = tuple(int(x) for x in parts[-3:])
            elif line.startswith("origin"):
                parts = line.split()
                origin = np.array([float(x) for x in parts[1:4]])
            elif line.startswith("delta"):
                parts = line.split()
                delta.append([float(x) for x in parts[1:4]])
            elif line.startswith("object 2") or line.startswith("object 3"):
                continue
            elif line.startswith("attribute") or line.startswith("object"):
                continue
            else:
                data_values.extend(float(v) for v in line.split())

    data = np.array(data_values).reshape(counts)
    return data, origin, np.array(delta)


def _extract_peaks(data, origin, delta, threshold):
    """密度がthresholdを超えるグリッド点を座標として返す。

    ピーク間の距離が近すぎるものを間引くため、単純な閾値カットを使用。
    Returns:
        coords: list of (x, y, z) in Å
    """
    indices = np.argwhere(data > threshold)
    coords = []
    for idx in indices:
        coord = origin + delta @ idx
        coords.append(coord)
    return coords


def _write_pdb(coords, solvent, output_path):
    """ピーク座標をPDB形式で書き出す。"""
    atom_name = "O" if solvent == "water" else "X"
    resname = "WAT" if solvent == "water" else "SOL"

    with open(output_path, "w") as f:
        for i, (x, y, z) in enumerate(coords, start=1):
            f.write(
                f"ATOM  {i:5d}  {atom_name:<3s} {resname} A{i:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00\n"
            )
        f.write("END\n")

    LOGGER.info(f"{len(coords)} solvent sites written to {output_path}")


def run(args):
    # rism3d.snglpnt コマンド実行
    rism_cmd = (
        f"rism3d.snglpnt"
        f" --prmtop {args.prmtop}"
        f" --rst7 {args.coord}"
        f" --temperature {args.temperature}"
        f" --ng3 -1,-1,-1"
        f" --solvbox -1,-1,-1"
        f" --buffer 14.0"
        f" --tolerance 1e-5"
        f" --closure kh"
        f" --write_thermo 1"
        f" --ntwrism 1"
        f" --verbose 1"
    )

    LOGGER.info("Running rism3d.snglpnt ...")
    subprocess.run(rism_cmd, shell=True, check=True)

    # 出力 .dx ファイルのうち水の酸素密度を使用
    # rism3d.snglpnt は <prmtop_stem>.O.1.dx 等を出力する
    prmtop_stem = Path(args.prmtop).stem
    if args.solvent == "water":
        dx_candidates = list(Path(".").glob(f"{prmtop_stem}*O*.dx"))
        if not dx_candidates:
            dx_candidates = list(Path(".").glob("*.O.*.dx"))
    else:
        dx_candidates = list(Path(".").glob("*.dx"))

    if not dx_candidates:
        LOGGER.error("No .dx file found. rism3d.snglpnt may have failed.")
        return

    dx_path = dx_candidates[0]
    LOGGER.info(f"Reading density from {dx_path}")

    data, origin, delta = _parse_dx(dx_path)
    coords = _extract_peaks(data, origin, delta, args.threshold)
    _write_pdb(coords, args.solvent, args.output)

    if not args.keepfiles:
        for f in Path(".").glob("*.dx"):
            f.unlink()
        for f in Path(".").glob("rism3d*"):
            f.unlink()
        LOGGER.info("Intermediate files removed")
