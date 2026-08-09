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


def _atom_selection_mask(atom_count, indices):
    """Return an unambiguous ParmEd boolean atom-selection mask."""
    mask = [False] * atom_count
    for index in indices:
        if index < 0 or index >= atom_count:
            raise IndexError(f"Atom index out of range: {index}")
        mask[index] = True
    return mask


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
        return {
            "water_molecules": 0,
            "water_atoms": 0,
            "transfer_max_abs_error_A": 0.0,
            "saved_max_abs_error_A": 0.0,
            "vectorized_transfer": True,
        }

    atoms_per_water = len(water_groups[0])
    if any(len(group) != atoms_per_water for group in water_groups):
        raise ValueError("Water residues must have the same number of atoms")
    if structure.box is None:
        raise ValueError("Periodic box dimensions are required")
    if any(abs(angle - 90.0) > 1e-3 for angle in structure.box[3:]):
        raise ValueError("Packmol water placement requires an orthorhombic box")

    water_indices = [index for group in water_groups for index in group]
    water_index_set = set(water_indices)
    fixed_indices = [
        atom.idx for atom in structure.atoms if atom.idx not in water_index_set
    ]

    with tempfile.TemporaryDirectory(prefix="mdtbx-packmol-") as tempdir:
        workdir = Path(tempdir)
        fixed_path = workdir / "fixed.pdb"
        water_path = workdir / "water.pdb"
        packed_path = workdir / "packed.pdb"
        input_path = workdir / "packmol.inp"

        fixed_input = None
        if fixed_indices:
            fixed = structure[_atom_selection_mask(len(structure.atoms), fixed_indices)]
            fixed.box = None
            fixed.save(str(fixed_path), format="pdb", overwrite=True)
            fixed_input = fixed_path

        water = structure[_atom_selection_mask(len(structure.atoms), water_groups[0])]
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

        packed_coordinates = np.asarray(packed.coordinates)
        packed_water_coordinates = packed_coordinates[len(fixed_indices) :]
        if len(packed_water_coordinates) != len(water_indices):
            raise RuntimeError("Packmol water-coordinate count mismatch")

        coordinates = np.asarray(structure.coordinates).copy()
        coordinates[water_indices] = packed_water_coordinates
        structure.coordinates = coordinates
        applied_coordinates = np.asarray(structure.coordinates)[water_indices]
        transfer_max_abs_error = float(
            np.max(np.abs(applied_coordinates - packed_water_coordinates))
        )
        if transfer_max_abs_error > 1.0e-8:
            raise RuntimeError(
                "Packmol coordinates were not transferred to the Amber structure"
            )

    structure.save(str(rst_path), format="rst7", overwrite=True)
    structure.save(str(pdb_path), format="pdb", overwrite=True)
    saved = pmd.load_file(str(parm_path), xyz=str(rst_path))
    saved_water_coordinates = np.asarray(saved.coordinates)[water_indices]
    saved_max_abs_error = float(
        np.max(np.abs(saved_water_coordinates - packed_water_coordinates))
    )
    if saved_max_abs_error > 1.0e-3:
        raise RuntimeError("Saved Amber coordinates do not match Packmol output")
    LOGGER.info(f"Repacked {len(water_groups)} water molecules with Packmol")
    return {
        "water_molecules": len(water_groups),
        "water_atoms": len(water_indices),
        "transfer_max_abs_error_A": transfer_max_abs_error,
        "saved_max_abs_error_A": saved_max_abs_error,
        "vectorized_transfer": True,
    }
