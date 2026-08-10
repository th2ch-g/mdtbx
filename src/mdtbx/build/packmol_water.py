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


def _align_coordinates_to_reference(coordinates, reference_coordinates):
    """Remove a uniform coordinate-origin shift using a reference structure."""
    import numpy as np

    coordinates = np.asarray(coordinates)
    reference = np.asarray(reference_coordinates)
    if coordinates.shape != reference.shape:
        raise ValueError("Reference PDB atom count differs from the Amber restart")
    translation = np.median(coordinates - reference, axis=0)
    aligned = coordinates - translation
    max_error = float(np.max(np.abs(aligned - reference)))
    if max_error > 1.0e-3:
        raise ValueError("Amber restart and reference PDB coordinates are inconsistent")
    return aligned, translation, max_error


def _wrap_fixed_coordinates(structure, fixed_indices, box_lengths, *, coordinates=None):
    """Wrap fixed atoms without splitting a bonded component across the box."""
    import numpy as np

    if coordinates is None:
        coordinates = structure.coordinates
    coordinates = np.asarray(coordinates).copy()
    if not fixed_indices:
        return coordinates, 0

    fixed = np.asarray(fixed_indices, dtype=int)
    original = coordinates[fixed]
    wrapped = np.mod(original, box_lengths)
    local_index = {atom_index: index for index, atom_index in enumerate(fixed_indices)}
    for bond in getattr(structure, "bonds", ()):
        index_a = bond.atom1.idx
        index_b = bond.atom2.idx
        if index_a not in local_index or index_b not in local_index:
            continue
        before = original[local_index[index_b]] - original[local_index[index_a]]
        after = wrapped[local_index[index_b]] - wrapped[local_index[index_a]]
        if not np.allclose(before, after, atol=1.0e-6, rtol=0.0):
            raise ValueError(
                "Wrapping fixed coordinates would split a bonded component"
            )

    moved = np.any(np.abs(wrapped - original) > 1.0e-6, axis=1)
    coordinates[fixed] = wrapped
    return coordinates, int(np.count_nonzero(moved))


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
    reference = pmd.load_file(str(pdb_path))
    aligned_coordinates, origin_translation, origin_alignment_max_error = (
        _align_coordinates_to_reference(structure.coordinates, reference.coordinates)
    )
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
    box_lengths = np.asarray(structure.box[:3])
    coordinates, fixed_atoms_wrapped = _wrap_fixed_coordinates(
        structure,
        fixed_indices,
        box_lengths,
        coordinates=aligned_coordinates,
    )

    with tempfile.TemporaryDirectory(prefix="mdtbx-packmol-") as tempdir:
        workdir = Path(tempdir)
        fixed_path = workdir / "fixed.pdb"
        water_path = workdir / "water.pdb"
        packed_path = workdir / "packed.pdb"
        input_path = workdir / "packmol.inp"

        fixed_input = None
        if fixed_indices:
            fixed = structure[_atom_selection_mask(len(structure.atoms), fixed_indices)]
            fixed.coordinates = coordinates[fixed_indices]
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

        coordinates[fixed_indices] = packed_coordinates[: len(fixed_indices)]
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
        "fixed_atoms_wrapped": fixed_atoms_wrapped,
        "restart_origin_translation_A": origin_translation.tolist(),
        "origin_alignment_max_error_A": origin_alignment_max_error,
    }
