# References:
# - Tutorial 34: Solvation with 3D-RISM
#   https://ambermd.org/tutorials/advanced/tutorial34/index.html
# - Tutorial 40: 1D-RISM and 3D-RISM
#   https://ambermd.org/tutorials/advanced/tutorial40/index.php

import argparse
import math
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np

from ..logger import generate_logger
from ..utils.proc import run_cmd

LOGGER = generate_logger(__name__)

DEFAULT_SOLVENT_MODEL = "cSPCE"
DEFAULT_SOLVENT_DENSITY = 55.5
DEFAULT_DIELECTRIC = 78.44

# 3D-RISM input template for sander.
SANDER_RISM_INPUT_TEMPLATE = """\
&cntrl
  ntx=1, nstlim=0, irism=1,
/
&rism
  closure='{closure}',
  grdspc={grdspc},{grdspc},{grdspc},
  tolerance={tolerance},
  buffer={buffer},
  solvcut={solvcut},
  mdiis_del=0.7,
  mdiis_nvec=5,
  maxstep=10000,
  npropagate=5,
  verbose=2,
  apply_rism_force=0,
  volfmt='dx',
  ntwrism=1,
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
        "--xvv",
        default=None,
        type=str,
        help="Pre-computed solvent susceptibility file (.xvv). "
        "If not provided, rism1d is run to generate one.",
    )

    parser.add_argument(
        "--solvent",
        default="water",
        choices=["water"],
        help="Solvent type",
    )

    parser.add_argument(
        "--solvent-model",
        default=DEFAULT_SOLVENT_MODEL,
        type=str,
        help=(
            "Solvent model for 1D-RISM. Accepts a .mdl path or a model name "
            "under $AMBERHOME/dat/rism1d/mdl/"
        ),
    )

    parser.add_argument(
        "--solvent-density",
        default=DEFAULT_SOLVENT_DENSITY,
        type=float,
        help="Bulk solvent density [mol/L] for 1D-RISM",
    )

    parser.add_argument(
        "--dielectric",
        default=DEFAULT_DIELECTRIC,
        type=float,
        help="Bulk solvent dielectric constant for 1D-RISM",
    )

    parser.add_argument(
        "--xvv-output",
        type=str,
        help="Copy the prepared xvv file to this reusable output path",
    )

    parser.add_argument(
        "--temperature",
        default=300.0,
        type=float,
        help="Temperature [K]",
    )

    parser.add_argument(
        "--closure",
        default="kh",
        choices=["kh", "hnc", "pse2", "pse3"],
        help="Closure approximation",
    )

    parser.add_argument(
        "--grdspc",
        default=0.5,
        type=float,
        help="Grid spacing [Å]",
    )

    parser.add_argument(
        "--tolerance",
        default=1e-5,
        type=float,
        help="Convergence tolerance for 3D-RISM",
    )

    parser.add_argument(
        "--buffer",
        default=14.0,
        type=float,
        help="Buffer distance around solute [Å]",
    )

    parser.add_argument(
        "--solvcut",
        default=14.0,
        type=float,
        help="Solvent cutoff distance [Å]",
    )

    parser.add_argument(
        "--threshold",
        default=1.5,
        type=float,
        help="Density threshold for peak extraction (g(r) value)",
    )

    parser.add_argument(
        "--exclusion-radius",
        default=2.6,
        type=float,
        help="Minimum distance between placed solvent sites [Å]. "
        "Approximate diameter of a water molecule.",
    )

    parser.add_argument(
        "--max-sites",
        default=None,
        type=int,
        help="Maximum number of solvent sites to place. "
        "If not set, all sites above threshold are placed.",
    )

    parser.add_argument(
        "--use-sander",
        action="store_true",
        help="Use sander interface instead of rism3d.snglpnt",
    )

    parser.add_argument(
        "--keepfiles",
        action="store_true",
        help="Keep intermediate files",
    )

    parser.set_defaults(func=run)


# ---------------------------------------------------------------------------
#  1D-RISM: xvv generation
# ---------------------------------------------------------------------------


def _resolve_solvent_model_path(solvent_model):
    """Resolve a solvent model name or .mdl path to an existing file."""
    supplied_path = Path(solvent_model).expanduser()
    if supplied_path.suffix == ".mdl" or supplied_path.is_absolute():
        model_path = supplied_path.resolve()
    else:
        amberhome = os.environ.get("AMBERHOME")
        if not amberhome:
            raise EnvironmentError(
                "$AMBERHOME is not set; provide --solvent-model as a .mdl path"
            )
        model_path = (
            Path(amberhome) / "dat" / "rism1d" / "mdl" / f"{solvent_model}.mdl"
        ).resolve()

    if not model_path.is_file():
        raise FileNotFoundError(
            f"Solvent model file not found: {model_path}. "
            "Check $AMBERHOME and --solvent-model."
        )
    return model_path


def _fortran_quote(value):
    """Quote a string for a Fortran namelist."""
    return "'" + str(value).replace("'", "''") + "'"


def _run_rism1d(
    solvent_model,
    temperature,
    workdir,
    solvent_density=DEFAULT_SOLVENT_DENSITY,
    dielectric=DEFAULT_DIELECTRIC,
):
    """Run 1D-RISM and return the generated xvv path.

    The model is loaded from $AMBERHOME/dat/rism1d/mdl/ unless an explicit
    .mdl path is supplied.

    Returns:
        Path to the generated xvv file.
    """
    model_path = _resolve_solvent_model_path(solvent_model)
    xvv_stem = f"{model_path.stem}_{temperature:.2f}"
    workpath = Path(workdir)
    workpath.mkdir(parents=True, exist_ok=True)
    xvv_path = workpath / f"{xvv_stem}.xvv"

    inp_content = (
        "&PARAMETERS\n"
        "  THEORY='DRISM', CLOSURE='KH',\n"
        "  NR=16384, DR=0.025,\n"
        "  OUTLIST='x',\n"
        "  MDIIS_NVEC=20, MDIIS_DEL=0.3, TOLERANCE=1.0e-12,\n"
        "  KSAVE=-1, MAXSTEP=10000,\n"
        "  SMEAR=1, ADBCOR=0.5,\n"
        f"  TEMPERATURE={temperature}, DIEPS={dielectric}, NSP=1,\n"
        "/\n"
        "&SPECIES\n"
        f"  DENSITY={solvent_density}, UNITS='M',\n"
        f"  MODEL={_fortran_quote(model_path)},\n"
        "/\n"
    )
    inp_path = workpath / f"{xvv_stem}.inp"
    inp_path.write_text(inp_content)

    LOGGER.info(f"Running 1D-RISM to generate xvv file ({solvent_model}) ...")
    with (workpath / f"{xvv_stem}.out").open("w") as output:
        run_cmd(
            ["rism1d", xvv_stem],
            cwd=workdir,
            stdout=output,
            stderr=subprocess.STDOUT,
        )

    if not xvv_path.is_file():
        raise RuntimeError(
            f"1D-RISM did not produce {xvv_path}. Check the output for errors."
        )

    LOGGER.info(f"Generated xvv file: {xvv_path}")
    return str(xvv_path)


# ---------------------------------------------------------------------------
#  3D-RISM execution
# ---------------------------------------------------------------------------


def _run_rism3d_snglpnt(prmtop, coord, xvv_path, args, workdir):
    """Run 3D-RISM through the rism3d.snglpnt interface."""

    prmtop_abs = os.path.abspath(prmtop)
    coord_abs = os.path.abspath(coord)
    xvv_abs = os.path.abspath(xvv_path)
    prmtop_stem = Path(prmtop).stem

    guv_prefix = os.path.join(workdir, prmtop_stem)
    pdb_path = os.path.join(workdir, f"{prmtop_stem}.pdb")
    with open(pdb_path, "w") as pdb_output:
        run_cmd(
            ["ambpdb", "-p", prmtop_abs, "-c", coord_abs],
            stdout=pdb_output,
        )

    rism_cmd = [
        "rism3d.snglpnt",
        "--pdb",
        pdb_path,
        "--prmtop",
        prmtop_abs,
        "--rst",
        coord_abs,
        "--xvv",
        xvv_abs,
        "--closure",
        args.closure,
        "--grdspc",
        f"{args.grdspc},{args.grdspc},{args.grdspc}",
        "--buffer",
        str(args.buffer),
        "--solvcut",
        str(args.solvcut),
        "--tolerance",
        str(args.tolerance),
        "--verbose",
        "2",
        "--volfmt",
        "dx",
        "--guv",
        guv_prefix,
    ]

    LOGGER.info("Running rism3d.snglpnt ...")
    LOGGER.info(f"  Command: {' '.join(rism_cmd)}")
    run_cmd(rism_cmd, cwd=workdir)


def _run_sander_rism(prmtop, coord, xvv_path, args, workdir):
    """Run 3D-RISM through the sander interface."""

    prmtop_abs = os.path.abspath(prmtop)
    coord_abs = os.path.abspath(coord)
    xvv_abs = os.path.abspath(xvv_path)
    guv_prefix = os.path.join(workdir, f"{Path(prmtop).stem}.{args.closure}")

    mdin_content = SANDER_RISM_INPUT_TEMPLATE.format(
        closure=args.closure,
        grdspc=args.grdspc,
        tolerance=args.tolerance,
        buffer=args.buffer,
        solvcut=args.solvcut,
    )
    mdin_path = os.path.join(workdir, "mdin.rism")
    with open(mdin_path, "w") as f:
        f.write(mdin_content)

    sander_cmd = [
        "sander",
        "-O",
        "-i",
        mdin_path,
        "-o",
        os.path.join(workdir, "mdout"),
        "-p",
        prmtop_abs,
        "-c",
        coord_abs,
        "-xvv",
        xvv_abs,
        "-guv",
        guv_prefix,
    ]

    LOGGER.info("Running sander with 3D-RISM ...")
    run_cmd(sander_cmd, cwd=workdir)


# ---------------------------------------------------------------------------
#  DX parsing
# ---------------------------------------------------------------------------


def _parse_dx(dx_path):
    """Read an OpenDX density grid.

    Returns:
        data: ``(nx, ny, nz)`` array of g(r) values.
        origin: ``(3,)`` grid origin in angstroms.
        delta: ``(3, 3)`` grid-axis vectors in angstroms.
    """
    origin = None
    delta = []
    counts = None
    data_values = []

    with open(dx_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("#") or not line:
                continue
            if line.startswith("object 1"):
                parts = line.split()
                counts = tuple(int(x) for x in parts[-3:])
            elif line.startswith("origin"):
                parts = line.split()
                origin = np.array([float(x) for x in parts[1:4]])
            elif line.startswith("delta"):
                parts = line.split()
                delta.append([float(x) for x in parts[1:4]])
            elif (
                line.startswith("object")
                or line.startswith("attribute")
                or line.startswith("component")
            ):
                continue
            else:
                # Grid data line.
                data_values.extend(float(v) for v in line.split())

    if counts is None:
        raise ValueError(f"Could not parse grid dimensions from {dx_path}")
    if origin is None:
        raise ValueError(f"Could not parse origin from {dx_path}")
    if len(delta) != 3:
        raise ValueError(f"Expected 3 delta vectors, got {len(delta)}")
    if any(count < 1 for count in counts):
        raise ValueError(f"Grid dimensions must be positive: {counts}")

    expected_values = math.prod(counts)
    if len(data_values) != expected_values:
        raise ValueError(
            f"Expected {expected_values} grid values in {dx_path}, "
            f"got {len(data_values)}"
        )

    data = np.array(data_values).reshape(counts)
    delta = np.array(delta)
    if not np.all(np.isfinite(data)):
        raise ValueError(f"Grid data contains non-finite values: {dx_path}")
    if not np.all(np.isfinite(origin)) or not np.all(np.isfinite(delta)):
        raise ValueError(f"Grid geometry contains non-finite values: {dx_path}")

    return data, origin, delta


def _grid_to_cartesian(indices, origin, delta):
    """Convert ``(N, 3)`` grid indices to Cartesian coordinates.

    OpenDX stores one grid-axis vector per row of ``delta``:
    ``coordinates = origin + indices @ delta``.
    """
    return origin + indices @ delta


# ---------------------------------------------------------------------------
#  Greedy Placevent-style peak extraction
# ---------------------------------------------------------------------------


def _extract_peaks_greedy(
    data, origin, delta, threshold, exclusion_radius, max_sites=None
):
    """Extract solvent sites with a greedy Placevent-style algorithm.

    Candidates above ``threshold`` are considered in descending density order.
    After a site is selected, candidates within ``exclusion_radius`` are
    removed. Selection stops at ``max_sites`` when specified.

    Returns:
        coords: ``(M, 3)`` solvent-site coordinates in angstroms.
        gvalues: ``(M,)`` g(r) values at the selected sites.
    """
    data = np.asarray(data)
    origin = np.asarray(origin)
    delta = np.asarray(delta)
    if data.ndim != 3:
        raise ValueError("Density data must be a three-dimensional array")
    if origin.shape != (3,) or delta.shape != (3, 3):
        raise ValueError("Grid origin and delta must have shapes (3,) and (3, 3)")
    if not np.all(np.isfinite(data)):
        raise ValueError("Density data must contain only finite values")
    if not np.all(np.isfinite(origin)) or not np.all(np.isfinite(delta)):
        raise ValueError("Grid origin and delta must contain only finite values")
    if not math.isfinite(threshold):
        raise ValueError("threshold must be finite")
    if not math.isfinite(exclusion_radius) or exclusion_radius <= 0:
        raise ValueError("exclusion_radius must be positive and finite")
    if max_sites is not None and max_sites < 1:
        raise ValueError("max_sites must be positive")

    # Collect grid points above the density threshold.
    indices = np.argwhere(data > threshold)
    if len(indices) == 0:
        LOGGER.warning(
            f"No grid points exceed threshold g(r) > {threshold}. "
            "Try lowering --threshold."
        )
        return np.empty((0, 3)), np.empty(0)

    values = data[indices[:, 0], indices[:, 1], indices[:, 2]]

    # Sort candidates by descending g(r).
    order = np.argsort(-values)
    indices = indices[order]
    values = values[order]

    # Convert grid indices to Cartesian coordinates.
    all_coords = _grid_to_cartesian(indices.astype(float), origin, delta)

    # Select peaks greedily.
    placed_coords = []
    placed_gvalues = []
    used = np.zeros(len(all_coords), dtype=bool)
    excl_sq = exclusion_radius**2

    for i in range(len(all_coords)):
        if used[i]:
            continue

        coord_i = all_coords[i]
        placed_coords.append(coord_i)
        placed_gvalues.append(values[i])

        if max_sites is not None and len(placed_coords) >= max_sites:
            break

        # Exclude candidates within the selected site's radius.
        remaining = np.where(~used)[0]
        remaining = remaining[remaining > i]
        if len(remaining) > 0:
            diff = all_coords[remaining] - coord_i
            dist_sq = np.sum(diff**2, axis=1)
            too_close = remaining[dist_sq <= excl_sq]
            used[too_close] = True

    coords = np.array(placed_coords)
    gvalues = np.array(placed_gvalues)

    LOGGER.info(
        f"Extracted {len(coords)} solvent sites "
        f"(threshold={threshold}, exclusion_radius={exclusion_radius} Å)"
    )
    return coords, gvalues


# ---------------------------------------------------------------------------
#  PDB output
# ---------------------------------------------------------------------------


def _write_pdb(coords, gvalues, solvent, output_path):
    """Write solvent-site coordinates as PDB.

    Occupancy is fixed at 1.00 and each g(r) value is stored as the B-factor.
    Atom and residue identifiers wrap to their one-based PDB field ranges.
    """
    atom_name = "O" if solvent == "water" else "X"
    resname = "WAT" if solvent == "water" else "SOL"

    with open(output_path, "w") as f:
        f.write(
            "REMARK   Generated by place_solvent (3D-RISM)\n"
            f"REMARK   {len(coords)} solvent sites placed\n"
        )
        for i, (coord, gval) in enumerate(zip(coords, gvalues, strict=True), start=1):
            x, y, z = coord
            serial = (i - 1) % 99999 + 1  # PDB atom serial range: 1..99999.
            resseq = (i - 1) % 9999 + 1  # PDB residue sequence range: 1..9999.
            f.write(
                f"HETATM{serial:5d}  {atom_name:<3s} {resname} A"
                f"{resseq:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}"
                f"{1.00:6.2f}{gval:6.2f}\n"
            )
        f.write("END\n")

    LOGGER.info(f"{len(coords)} solvent sites written to {output_path}")


# ---------------------------------------------------------------------------
#  DX discovery
# ---------------------------------------------------------------------------


def _find_oxygen_dx(workdir, prmtop_stem, closure):
    """Find the oxygen pair-distribution DX output.

    Typical names are ``prefix.O.1.dx`` for rism3d.snglpnt and
    ``<stem>.<closure>.O.0.dx`` for sander.
    """
    workpath = Path(workdir)
    excluded_kinds = ("cuv", "huv", "uuv")
    candidates = [
        path
        for path in workpath.glob("*.dx")
        if ".o." in path.name.lower()
        and not any(kind in path.name.lower().split(".") for kind in excluded_kinds)
    ]
    if not candidates:
        raise FileNotFoundError(
            f"No oxygen pair-distribution .dx file found in {workdir}. "
            "Check the 3D-RISM output and solvent-site names."
        )

    prefix = prmtop_stem.lower()
    closure_token = f".{closure.lower()}."
    return min(
        candidates,
        key=lambda path: (
            not path.name.lower().startswith(prefix),
            closure_token not in path.name.lower(),
            path.name,
        ),
    )


# ---------------------------------------------------------------------------
#  Main workflow
# ---------------------------------------------------------------------------


def _validate_args(args):
    """Validate files and numeric options before starting external tools."""
    for option in ("prmtop", "coord"):
        path = Path(getattr(args, option)).expanduser()
        if not path.is_file():
            raise FileNotFoundError(f"{option} file not found: {path}")
    if args.xvv is not None and not Path(args.xvv).expanduser().is_file():
        raise FileNotFoundError(f"xvv file not found: {args.xvv}")

    positive_finite = {
        "temperature": args.temperature,
        "solvent_density": args.solvent_density,
        "dielectric": args.dielectric,
        "grdspc": args.grdspc,
        "tolerance": args.tolerance,
        "buffer": args.buffer,
        "solvcut": args.solvcut,
        "exclusion_radius": args.exclusion_radius,
    }
    for name, value in positive_finite.items():
        if not math.isfinite(value) or value <= 0:
            raise ValueError(f"--{name.replace('_', '-')} must be positive and finite")
    if not math.isfinite(args.threshold):
        raise ValueError("--threshold must be finite")
    if args.max_sites is not None and args.max_sites < 1:
        raise ValueError("--max-sites must be positive")


def _copy_xvv(xvv_path, output):
    """Copy an xvv file to a user-visible reusable location."""
    destination = Path(output).expanduser().resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    source = Path(xvv_path).resolve()
    if source != destination:
        shutil.copy2(source, destination)
    LOGGER.info(f"xvv file saved to {destination}")


def run(args):
    _validate_args(args)
    prmtop_stem = Path(args.prmtop).stem

    # Create an isolated working directory for intermediate RISM files.
    workdir = tempfile.mkdtemp(prefix="rism3d_")
    LOGGER.info(f"Working directory: {workdir}")

    try:
        # ---- 1. Prepare the xvv file. ----
        if args.xvv is not None:
            xvv_path = os.path.abspath(args.xvv)
            LOGGER.info(f"Using provided xvv file: {xvv_path}")
        else:
            xvv_path = _run_rism1d(
                args.solvent_model,
                args.temperature,
                workdir,
                solvent_density=args.solvent_density,
                dielectric=args.dielectric,
            )

        if args.xvv_output is not None:
            _copy_xvv(xvv_path, args.xvv_output)

        # ---- 2. Run 3D-RISM. ----
        if args.use_sander:
            _run_sander_rism(args.prmtop, args.coord, xvv_path, args, workdir)
        else:
            _run_rism3d_snglpnt(args.prmtop, args.coord, xvv_path, args, workdir)

        # ---- 3. Read the oxygen-density DX file. ----
        dx_path = _find_oxygen_dx(workdir, prmtop_stem, args.closure)
        LOGGER.info(f"Reading density from {dx_path}")

        data, origin, delta = _parse_dx(str(dx_path))
        LOGGER.info(
            f"Grid dimensions: {data.shape}, "
            f"origin: {origin}, "
            f"spacing: {np.diag(delta)}"
        )

        # ---- 4. Extract peaks with greedy spatial exclusion. ----
        coords, gvalues = _extract_peaks_greedy(
            data,
            origin,
            delta,
            threshold=args.threshold,
            exclusion_radius=args.exclusion_radius,
            max_sites=args.max_sites,
        )

        if len(coords) == 0:
            raise RuntimeError("No solvent sites found; try lowering --threshold")

        # ---- 5. Write PDB output. ----
        output_path = Path(args.output).expanduser().resolve()
        output_path.parent.mkdir(parents=True, exist_ok=True)
        _write_pdb(coords, gvalues, args.solvent, str(output_path))

    finally:
        if not args.keepfiles:
            shutil.rmtree(workdir, ignore_errors=True)
            LOGGER.info("Intermediate files removed")
        else:
            LOGGER.info(f"Intermediate files kept in {workdir}")
