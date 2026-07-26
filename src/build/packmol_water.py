import tempfile
from pathlib import Path

from ..logger import generate_logger
from ..utils.proc import run_cmd

LOGGER = generate_logger(__name__)

WATER_RESNAMES = frozenset({"WAT", "HOH", "SOL"})


def _build_packmol_input(
    *,
    fixed_path,
    water_path,
    output_path,
    water_count,
    box,
    seed,
    tolerance,
):
    x, y, z = box[:3]
    margin = tolerance / 2
    lines = [
        f"tolerance {tolerance}",
        "filetype pdb",
        f"output {output_path}",
        f"seed {seed}",
        "add_amber_ter",
    ]
    if fixed_path is not None:
        lines.extend(
            [
                f"structure {fixed_path}",
                "  number 1",
                "  fixed 0.0 0.0 0.0 0.0 0.0 0.0",
                "end structure",
            ]
        )
    lines.extend(
        [
            f"structure {water_path}",
            f"  number {water_count}",
            (
                f"  inside box {margin} {margin} {margin} "
                f"{x - margin} {y - margin} {z - margin}"
            ),
            "end structure",
        ]
    )
    return "\n".join(lines) + "\n"


def _water_atom_groups(structure, water_resnames):
    water_names = {name.upper() for name in water_resnames}
    return [
        [atom.idx for atom in residue.atoms]
        for residue in structure.residues
        if residue.name.upper() in water_names
    ]


def repack_water(
    parm_path,
    rst_path,
    pdb_path,
    *,
    seed=2026,
    tolerance=2.0,
    water_resnames=WATER_RESNAMES,
):
    import numpy as np
    import parmed as pmd

    structure = pmd.load_file(str(parm_path), xyz=str(rst_path))
    water_groups = _water_atom_groups(structure, water_resnames)
    if not water_groups:
        LOGGER.info("No water molecules to repack")
        return

    atoms_per_water = len(water_groups[0])
    if any(len(group) != atoms_per_water for group in water_groups):
        raise ValueError("Water residues must have the same number of atoms")
    if structure.box is None:
        raise ValueError("Periodic box dimensions are required")
    if any(abs(angle - 90.0) > 1e-3 for angle in structure.box[3:]):
        raise ValueError("Packmol water placement requires an orthorhombic box")

    water_indices = {index for group in water_groups for index in group}
    fixed_indices = [
        atom.idx for atom in structure.atoms if atom.idx not in water_indices
    ]

    with tempfile.TemporaryDirectory(prefix="mdtbx-packmol-") as tempdir:
        workdir = Path(tempdir)
        fixed_path = workdir / "fixed.pdb"
        water_path = workdir / "water.pdb"
        packed_path = workdir / "packed.pdb"
        input_path = workdir / "packmol.inp"

        fixed_input = None
        if fixed_indices:
            fixed = structure[fixed_indices]
            fixed.box = None
            fixed.save(str(fixed_path), format="pdb", overwrite=True)
            fixed_input = fixed_path

        water = structure[water_groups[0]]
        water.box = None
        water.save(str(water_path), format="pdb", overwrite=True)

        packmol_input = _build_packmol_input(
            fixed_path=fixed_input,
            water_path=water_path,
            output_path=packed_path,
            water_count=len(water_groups),
            box=structure.box,
            seed=seed,
            tolerance=tolerance,
        )
        input_path.write_text(packmol_input)
        with input_path.open() as input_handle:
            result = run_cmd(
                ["packmol"],
                stdin=input_handle,
                cwd=workdir,
                capture_output=True,
                text=True,
                check=False,
            )
        if (
            result.returncode != 0
            or not packed_path.is_file()
            or "Success!" not in result.stdout
        ):
            raise RuntimeError(f"Packmol failed:\n{result.stdout}\n{result.stderr}")

        packed = pmd.load_file(str(packed_path))
        expected_atoms = len(fixed_indices) + len(water_groups) * atoms_per_water
        if len(packed.atoms) != expected_atoms:
            raise RuntimeError(
                f"Packmol output has {len(packed.atoms)} atoms; "
                f"expected {expected_atoms}"
            )

        offset = len(fixed_indices)
        packed_coordinates = np.asarray(packed.coordinates)
        for water_number, atom_indices in enumerate(water_groups):
            start = offset + water_number * atoms_per_water
            stop = start + atoms_per_water
            structure.coordinates[atom_indices] = packed_coordinates[start:stop]

    structure.save(str(rst_path), format="rst7", overwrite=True)
    structure.save(str(pdb_path), format="pdb", overwrite=True)
    LOGGER.info(f"Repacked {len(water_groups)} water molecules with Packmol")
