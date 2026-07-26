import json
import math
from pathlib import Path

import numpy as np

ABFE_MANIFEST = "abfe_manifest.json"
ABFE_SCHEMA_VERSION = 1
GAS_CONSTANT_KJ = 0.00831446261815324
STANDARD_VOLUME_NM3 = 1.6605390671738466


def _anchor_geometry_from_frame(frame, anchors):
    import mdtraj as md

    if len(anchors) != 6 or len(set(anchors)) != 6:
        raise ValueError("Six distinct Boresch anchor atoms are required")
    zero_based = [index - 1 for index in anchors]
    if any(index < 0 or index >= frame.n_atoms for index in zero_based):
        raise ValueError("Boresch anchor index is outside the structure")
    distance = float(
        md.compute_distances(frame, [[zero_based[2], zero_based[3]]])[0, 0]
    )
    angles = [
        float(value)
        for value in md.compute_angles(
            frame,
            [
                [zero_based[1], zero_based[2], zero_based[3]],
                [zero_based[2], zero_based[3], zero_based[4]],
            ],
        )[0]
    ]
    dihedrals = [
        float(value)
        for value in md.compute_dihedrals(
            frame,
            [
                [zero_based[0], zero_based[1], zero_based[2], zero_based[3]],
                [zero_based[1], zero_based[2], zero_based[3], zero_based[4]],
                [zero_based[2], zero_based[3], zero_based[4], zero_based[5]],
            ],
        )[0]
    ]
    if distance <= 0:
        raise ValueError("Boresch anchor distance must be positive")
    for angle in angles:
        degrees = math.degrees(angle)
        if not 45.0 <= degrees <= 135.0:
            raise ValueError("Boresch anchor angles must be between 45 and 135 degrees")
    return {
        "distance_nm": distance,
        "angles_rad": angles,
        "dihedrals_rad": dihedrals,
    }


def calculate_anchor_geometry(structure, anchors):
    import mdtraj as md

    return _anchor_geometry_from_frame(md.load(str(structure)), anchors)


def _bond_graph(frame, indices, cutoff):
    index_set = set(int(index) for index in indices)
    graph = {index: set() for index in index_set}
    pairs = []
    ordered = sorted(index_set)
    for position, atom_a in enumerate(ordered):
        pairs.extend((atom_a, atom_b) for atom_b in ordered[position + 1 :])
    if not pairs:
        return graph
    distances = frame.xyz[0]
    for atom_a, atom_b in pairs:
        delta = distances[atom_a] - distances[atom_b]
        if float(np.linalg.norm(delta)) <= cutoff:
            graph[atom_a].add(atom_b)
            graph[atom_b].add(atom_a)
    return graph


def select_boresch_anchors(
    structure,
    *,
    receptor_selection,
    ligand_selection,
    search_distance=0.5,
    bond_cutoff=0.22,
    receptor_atom_names=("CB", "CA", "C", "N", "O"),
):
    import mdtraj as md

    if search_distance <= 0 or bond_cutoff <= 0:
        raise ValueError("Anchor search distances must be positive")
    frame = md.load(str(structure))
    topology = frame.topology
    receptor = set(int(index) for index in topology.select(receptor_selection))
    ligand = set(int(index) for index in topology.select(ligand_selection))
    heavy = set(int(index) for index in topology.select("not element H"))
    receptor &= heavy
    ligand &= heavy
    receptor = {
        index for index in receptor if topology.atom(index).name in receptor_atom_names
    }
    if not receptor or not ligand:
        raise ValueError("Receptor or ligand anchor selection is empty")
    if receptor & ligand:
        raise ValueError("Receptor and ligand selections overlap")

    receptor_graph = _bond_graph(frame, receptor, bond_cutoff)
    ligand_graph = _bond_graph(frame, ligand, bond_cutoff)
    candidates = []
    xyz = frame.xyz[0]
    for receptor_3 in receptor:
        for ligand_1 in ligand:
            distance = float(np.linalg.norm(xyz[receptor_3] - xyz[ligand_1]))
            if distance > search_distance:
                continue
            for receptor_2 in receptor_graph[receptor_3]:
                for receptor_1 in receptor_graph[receptor_2] - {receptor_3}:
                    for ligand_2 in ligand_graph[ligand_1]:
                        for ligand_3 in ligand_graph[ligand_2] - {ligand_1}:
                            anchors = [
                                receptor_1 + 1,
                                receptor_2 + 1,
                                receptor_3 + 1,
                                ligand_1 + 1,
                                ligand_2 + 1,
                                ligand_3 + 1,
                            ]
                            try:
                                geometry = _anchor_geometry_from_frame(
                                    frame,
                                    anchors,
                                )
                            except ValueError:
                                continue
                            angular_score = sum(
                                abs(math.cos(angle)) for angle in geometry["angles_rad"]
                            )
                            candidates.append(
                                (distance + angular_score, anchors, geometry)
                            )
    if not candidates:
        raise ValueError("Could not identify a valid Boresch anchor set")
    _score, anchors, geometry = min(candidates, key=lambda item: item[0])
    return anchors, geometry


def boresch_pull_settings(
    geometry,
    *,
    distance_spring,
    angle_spring,
    dihedral_spring,
    release,
):
    for value in (distance_spring, angle_spring, dihedral_spring):
        if value <= 0:
            raise ValueError("Boresch spring constants must be positive")
    distance = geometry["distance_nm"]
    angle_a, angle_b = [math.degrees(value) for value in geometry["angles_rad"]]
    dihedral_a, dihedral_b, dihedral_c = [
        math.degrees(value) for value in geometry["dihedrals_rad"]
    ]
    settings = {
        "pull": "yes",
        "pull-nstxout": 0,
        "pull-nstfout": 0,
        "pull-ngroups": 6,
        "pull-ncoords": 6,
        "pull-pbc-ref-prev-step-com": "yes",
    }
    for index in range(1, 7):
        settings[f"pull-group{index}-name"] = f"ABFE_{index}"
        settings[f"pull-group{index}-pbcatom"] = "0"

    definitions = [
        ("distance", "3 4", distance, distance_spring),
        ("angle", "3 2 3 4", angle_a, angle_spring),
        ("angle", "4 3 4 5", angle_b, angle_spring),
        ("dihedral", "1 2 2 3 3 4", dihedral_a, dihedral_spring),
        ("dihedral", "2 3 3 4 4 5", dihedral_b, dihedral_spring),
        ("dihedral", "3 4 4 5 5 6", dihedral_c, dihedral_spring),
    ]
    for index, (geometry_name, groups, reference, spring) in enumerate(
        definitions,
        start=1,
    ):
        settings[f"pull-coord{index}-type"] = "umbrella"
        settings[f"pull-coord{index}-geometry"] = geometry_name
        settings[f"pull-coord{index}-groups"] = groups
        settings[f"pull-coord{index}-dim"] = "Y Y Y"
        settings[f"pull-coord{index}-start"] = "no"
        settings[f"pull-coord{index}-init"] = f"{reference:.10g}"
        settings[f"pull-coord{index}-k"] = f"{spring:.10g}"
        settings[f"pull-coord{index}-kB"] = "0" if release else f"{spring:.10g}"
    return settings


def write_anchor_index(path, anchors, base_index=None):
    output = ""
    if base_index is not None:
        output = Path(base_index).read_text()
        if output and not output.endswith("\n"):
            output += "\n"
    for index, atom_index in enumerate(anchors, start=1):
        output += f"[ ABFE_{index} ]\n{atom_index}\n\n"
    Path(path).write_text(output)


def boresch_standard_state_correction(
    geometry,
    *,
    temperature,
    distance_spring,
    angle_spring,
    dihedral_spring,
    symmetry_number=1,
):
    if temperature <= 0:
        raise ValueError("Temperature must be positive")
    if symmetry_number <= 0:
        raise ValueError("Symmetry number must be positive")
    distance = geometry["distance_nm"]
    angle_a, angle_b = geometry["angles_rad"]
    sine_product = math.sin(angle_a) * math.sin(angle_b)
    if distance <= 0 or sine_product <= 0:
        raise ValueError("Invalid Boresch reference geometry")
    rt = GAS_CONSTANT_KJ * temperature
    log_term = (
        math.log(8.0 * math.pi**2 * STANDARD_VOLUME_NM3 / symmetry_number)
        + 0.5 * math.log(distance_spring)
        + math.log(angle_spring)
        + 1.5 * math.log(dihedral_spring)
        - 2.0 * math.log(distance)
        - math.log(sine_product)
        - 3.0 * math.log(2.0 * math.pi * rt)
    )
    return rt * log_term


def load_abfe_manifest(path):
    source = Path(path).expanduser()
    manifest_path = source if source.is_file() else source / ABFE_MANIFEST
    if not manifest_path.is_file():
        raise FileNotFoundError(f"ABFE manifest not found: {manifest_path}")
    data = json.loads(manifest_path.read_text())
    if data.get("schema_version") != ABFE_SCHEMA_VERSION:
        raise ValueError("Unsupported ABFE manifest schema")
    legs = data.get("legs")
    expected = {
        "solvent_charge",
        "solvent_vdw",
        "complex_charge",
        "complex_vdw",
        "complex_restraint",
    }
    if not isinstance(legs, dict) or set(legs) != expected:
        raise ValueError("ABFE manifest does not contain the required legs")
    return manifest_path.parent, data
