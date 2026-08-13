import json
import math
from pathlib import Path

import numpy as np

# GAS_CONSTANT_KJ is also re-exported for existing importers of this module.
from ..config import GAS_CONSTANT_KJ

ABFE_MANIFEST = "abfe_manifest.json"
ABFE_SCHEMA_VERSION = 1
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


def _trajectory_chunks(trajectory, topology, stride):
    import mdtraj as md

    if stride <= 0:
        raise ValueError("Anchor trajectory stride must be positive")
    yield from md.iterload(
        str(trajectory),
        top=str(topology),
        stride=stride,
        chunk=1000,
    )


def _trajectory_anchor_statistics(trajectory, topology, anchor_sets, stride):
    import mdtraj as md

    if not anchor_sets:
        raise ValueError("No Boresch anchor candidates were provided")
    zero_based = [[index - 1 for index in anchors] for anchors in anchor_sets]
    distance_pairs = [[anchors[2], anchors[3]] for anchors in zero_based]
    angle_indices = [
        indices
        for anchors in zero_based
        for indices in (
            [anchors[1], anchors[2], anchors[3]],
            [anchors[2], anchors[3], anchors[4]],
        )
    ]
    dihedral_indices = [
        indices
        for anchors in zero_based
        for indices in (
            [anchors[0], anchors[1], anchors[2], anchors[3]],
            [anchors[1], anchors[2], anchors[3], anchors[4]],
            [anchors[2], anchors[3], anchors[4], anchors[5]],
        )
    ]
    candidate_count = len(anchor_sets)
    distance_sum = np.zeros(candidate_count)
    distance_square_sum = np.zeros(candidate_count)
    angle_sum = np.zeros((candidate_count, 2))
    angle_square_sum = np.zeros((candidate_count, 2))
    dihedral_sin_sum = np.zeros((candidate_count, 3))
    dihedral_cos_sum = np.zeros((candidate_count, 3))
    frame_count = 0

    for chunk in _trajectory_chunks(trajectory, topology, stride):
        if chunk.n_frames == 0:
            continue
        if any(
            index < 0 or index >= chunk.n_atoms
            for anchors in zero_based
            for index in anchors
        ):
            raise ValueError("Boresch anchor index is outside the trajectory topology")
        distances = md.compute_distances(chunk, distance_pairs)
        angles = md.compute_angles(chunk, angle_indices).reshape(
            chunk.n_frames,
            candidate_count,
            2,
        )
        dihedrals = md.compute_dihedrals(chunk, dihedral_indices).reshape(
            chunk.n_frames,
            candidate_count,
            3,
        )
        distance_sum += np.sum(distances, axis=0)
        distance_square_sum += np.sum(distances**2, axis=0)
        angle_sum += np.sum(angles, axis=0)
        angle_square_sum += np.sum(angles**2, axis=0)
        dihedral_sin_sum += np.sum(np.sin(dihedrals), axis=0)
        dihedral_cos_sum += np.sum(np.cos(dihedrals), axis=0)
        frame_count += chunk.n_frames
    if frame_count == 0:
        raise ValueError("Anchor trajectory contains no sampled frames")

    distance_mean = distance_sum / frame_count
    distance_variance = np.maximum(
        distance_square_sum / frame_count - distance_mean**2,
        0.0,
    )
    angle_mean = angle_sum / frame_count
    angle_variance = np.maximum(
        angle_square_sum / frame_count - angle_mean**2,
        0.0,
    )
    mean_sin = dihedral_sin_sum / frame_count
    mean_cos = dihedral_cos_sum / frame_count
    dihedral_mean = np.arctan2(mean_sin, mean_cos)
    resultant = np.clip(np.hypot(mean_sin, mean_cos), 1.0e-15, 1.0)
    dihedral_variance = -2.0 * np.log(resultant)

    results = []
    for index, anchors in enumerate(anchor_sets):
        geometry = {
            "distance_nm": float(distance_mean[index]),
            "angles_rad": [float(value) for value in angle_mean[index]],
            "dihedrals_rad": [float(value) for value in dihedral_mean[index]],
        }
        if geometry["distance_nm"] <= 0:
            continue
        if any(
            not 45.0 <= math.degrees(angle) <= 135.0 for angle in geometry["angles_rad"]
        ):
            continue
        standard_deviations = {
            "distance_nm": float(math.sqrt(distance_variance[index])),
            "angles_rad": [float(math.sqrt(value)) for value in angle_variance[index]],
            "dihedrals_rad": [
                float(math.sqrt(value)) for value in dihedral_variance[index]
            ],
        }
        results.append((anchors, geometry, standard_deviations, frame_count))
    return results


def calculate_trajectory_anchor_geometry(
    trajectory,
    topology,
    anchors,
    *,
    stride=1,
):
    if len(anchors) != 6 or len(set(anchors)) != 6:
        raise ValueError("Six distinct Boresch anchor atoms are required")
    results = _trajectory_anchor_statistics(
        trajectory,
        topology,
        [list(anchors)],
        stride,
    )
    if not results:
        raise ValueError(
            "Boresch anchor mean angles must be between 45 and 135 degrees"
        )
    _anchors, geometry, standard_deviations, frame_count = results[0]
    return geometry, standard_deviations, frame_count


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


def select_trajectory_boresch_anchors(
    trajectory,
    topology,
    *,
    receptor_selection,
    ligand_selection,
    distance_spring,
    angle_spring,
    dihedral_spring,
    stride=1,
    search_distance=0.5,
    bond_cutoff=0.22,
    receptor_atom_names=("CB", "CA", "C", "N", "O"),
):
    if search_distance <= 0 or bond_cutoff <= 0:
        raise ValueError("Anchor search distances must be positive")
    for value in (distance_spring, angle_spring, dihedral_spring):
        if value <= 0:
            raise ValueError("Boresch spring constants must be positive")
    chunks = _trajectory_chunks(trajectory, topology, stride)
    try:
        first_frame = next(chunks)[0]
    except StopIteration as error:
        raise ValueError("Anchor trajectory contains no sampled frames") from error

    frame_topology = first_frame.topology
    receptor = set(int(index) for index in frame_topology.select(receptor_selection))
    ligand = set(int(index) for index in frame_topology.select(ligand_selection))
    heavy = set(int(index) for index in frame_topology.select("not element H"))
    receptor &= heavy
    ligand &= heavy
    receptor = {
        index
        for index in receptor
        if frame_topology.atom(index).name in receptor_atom_names
    }
    if not receptor or not ligand:
        raise ValueError("Receptor or ligand anchor selection is empty")
    if receptor & ligand:
        raise ValueError("Receptor and ligand selections overlap")
    receptor_graph = _bond_graph(first_frame, receptor, bond_cutoff)
    ligand_graph = _bond_graph(first_frame, ligand, bond_cutoff)
    candidates = []
    xyz = first_frame.xyz[0]
    for receptor_3 in receptor:
        for ligand_1 in ligand:
            if float(np.linalg.norm(xyz[receptor_3] - xyz[ligand_1])) > search_distance:
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
                                _anchor_geometry_from_frame(first_frame, anchors)
                            except ValueError:
                                continue
                            candidates.append(anchors)
    candidates = [list(anchors) for anchors in dict.fromkeys(map(tuple, candidates))]
    if not candidates:
        raise ValueError("Could not identify a valid Boresch anchor set")

    stable_candidates = []
    for anchors, geometry, deviations, frame_count in _trajectory_anchor_statistics(
        trajectory,
        topology,
        candidates,
        stride,
    ):
        score = (
            distance_spring * deviations["distance_nm"] ** 2
            + angle_spring * sum(value**2 for value in deviations["angles_rad"])
            + dihedral_spring * sum(value**2 for value in deviations["dihedrals_rad"])
        )
        stable_candidates.append((score, anchors, geometry, deviations, frame_count))
    if not stable_candidates:
        raise ValueError("Could not identify a stable Boresch anchor set")
    score, anchors, geometry, deviations, frame_count = min(
        stable_candidates,
        key=lambda item: item[0],
    )
    return anchors, geometry, deviations, frame_count, score


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
