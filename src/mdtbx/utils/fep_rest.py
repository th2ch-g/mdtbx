import math
import re
from dataclasses import dataclass


@dataclass(frozen=True)
class AtomRecord:
    line_index: int
    molecule_type: str
    index: int
    residue_id: int
    residue_name: str
    atom_name: str
    type_a: str
    charge_a: float
    mass_a: float
    type_b: str
    charge_b: float
    mass_b: float
    has_b_state: bool


@dataclass(frozen=True)
class TopologyData:
    atoms: dict[str, list[AtomRecord]]
    molecule_counts: list[tuple[str, int]]
    epsilon_by_type: dict[str, float]


def _line_fields(raw_line):
    return raw_line.partition(";")[0].split()


def _replace_fields(raw_line, fields):
    _code, separator, comment = raw_line.partition(";")
    output = " ".join(fields)
    if separator:
        output += f" ;{comment.rstrip()}"
    return output + "\n"


def parse_preprocessed_topology(text):
    section = None
    molecule_type = None
    atoms = {}
    molecule_counts = []
    epsilon_by_type = {}

    for line_index, raw_line in enumerate(text.splitlines(keepends=True)):
        fields = _line_fields(raw_line)
        if not fields:
            continue
        if fields[0] == "[" and fields[-1] == "]":
            section = fields[1].lower()
            if section in {"system", "molecules"}:
                molecule_type = None
            continue
        if fields[0].startswith("#"):
            raise ValueError("FEP/REST requires a preprocessed topology")

        if section == "atomtypes":
            if len(fields) < 6:
                raise ValueError(f"Malformed atom type at line {line_index + 1}")
            epsilon_by_type[fields[0]] = float(fields[-1])
        elif section == "moleculetype":
            molecule_type = fields[0]
            atoms.setdefault(molecule_type, [])
            section = None
        elif section == "atoms":
            if molecule_type is None:
                raise ValueError(
                    f"Atoms section has no molecule type at line {line_index + 1}"
                )
            if len(fields) < 8:
                raise ValueError(f"Malformed atom record at line {line_index + 1}")
            has_b_state = len(fields) >= 11
            if len(fields) not in {8, 11}:
                raise ValueError(
                    f"Atom record must contain 8 or 11 fields at line {line_index + 1}"
                )
            type_a = fields[1]
            charge_a = float(fields[6])
            mass_a = float(fields[7])
            atoms[molecule_type].append(
                AtomRecord(
                    line_index=line_index,
                    molecule_type=molecule_type,
                    index=int(fields[0]),
                    residue_id=int(fields[2]),
                    residue_name=fields[3],
                    atom_name=fields[4],
                    type_a=type_a,
                    charge_a=charge_a,
                    mass_a=mass_a,
                    type_b=fields[8] if has_b_state else type_a,
                    charge_b=float(fields[9]) if has_b_state else charge_a,
                    mass_b=float(fields[10]) if has_b_state else mass_a,
                    has_b_state=has_b_state,
                )
            )
        elif section == "molecules":
            if len(fields) < 2:
                raise ValueError(f"Malformed molecule count at line {line_index + 1}")
            molecule_counts.append((fields[0], int(fields[1])))

    if not atoms or not molecule_counts:
        raise ValueError("Topology must contain molecule types and a molecules section")
    for name, count in molecule_counts:
        if name not in atoms:
            raise ValueError(f"Molecule type '{name}' has no atom definition")
        if count <= 0:
            raise ValueError(f"Molecule count for '{name}' must be positive")
    return TopologyData(atoms, molecule_counts, epsilon_by_type)


def build_fep_rest_schedule(
    replicas,
    temperature,
    max_temperature,
    *,
    coordinates=None,
):
    if replicas < 4:
        raise ValueError("FEP/REST requires at least 4 replicas")
    if temperature <= 0 or max_temperature < temperature:
        raise ValueError("REST temperatures must satisfy 0 < temperature <= maximum")

    if coordinates is None:
        boundary_1 = max(1, round((replicas - 1) / 3))
        boundary_2 = min(replicas - 2, round(2 * (replicas - 1) / 3))
        if boundary_2 <= boundary_1:
            boundary_2 = boundary_1 + 1

        coordinates = []
        for index in range(replicas):
            if index <= boundary_1:
                coordinate = index / boundary_1
            elif index <= boundary_2:
                coordinate = 1.0 + (index - boundary_1) / (boundary_2 - boundary_1)
            else:
                coordinate = 2.0 + (index - boundary_2) / (replicas - 1 - boundary_2)
            coordinates.append(coordinate)
    else:
        coordinates = [float(value) for value in coordinates]
        if len(coordinates) != replicas:
            raise ValueError("FEP/REST coordinate count must match the replica count")
        if not all(math.isfinite(value) for value in coordinates):
            raise ValueError("FEP/REST coordinates must be finite")
        if not math.isclose(coordinates[0], 0.0) or not math.isclose(
            coordinates[-1], 3.0
        ):
            raise ValueError("FEP/REST coordinates must start at 0 and end at 3")
        if any(
            current >= following
            for current, following in zip(coordinates, coordinates[1:])
        ):
            raise ValueError("FEP/REST coordinates must be strictly increasing")
        if not all(
            any(math.isclose(value, boundary) for value in coordinates)
            for boundary in (1.0, 2.0)
        ):
            raise ValueError("FEP/REST coordinates must contain boundaries 1 and 2")

    general_lambdas = [coordinate / 3.0 for coordinate in coordinates]
    vdw_lambdas = [
        0.0 if coordinate <= 1.0 else min(coordinate - 1.0, 1.0)
        for coordinate in coordinates
    ]
    charge_a_lambdas = [min(coordinate, 1.0) for coordinate in coordinates]
    charge_b_lambdas = [max(coordinate - 2.0, 0.0) for coordinate in coordinates]

    log_ratio = math.log(max_temperature / temperature)
    rest_nmax = (replicas - 1) // 2
    rest_temperatures = []
    for index in range(replicas):
        rest_state = min(index, replicas - 1 - index)
        fraction = rest_state / rest_nmax
        rest_temperatures.append(temperature * math.exp(log_ratio * fraction))
    scales = [temperature / value for value in rest_temperatures]

    return {
        "coordinates": coordinates,
        "fep": general_lambdas,
        "vdw": vdw_lambdas,
        "charge_a": charge_a_lambdas,
        "charge_b": charge_b_lambdas,
        "rest_temperatures": rest_temperatures,
        "scales": scales,
    }


def perturbed_global_indices(topology):
    perturbed_local = {}
    for name, records in topology.atoms.items():
        perturbed_local[name] = {
            record.index
            for record in records
            if record.has_b_state
            and (
                record.type_a != record.type_b
                or not math.isclose(record.charge_a, record.charge_b, abs_tol=1e-12)
                or not math.isclose(record.mass_a, record.mass_b, abs_tol=1e-12)
            )
        }

    global_indices = []
    offset = 0
    for name, count in topology.molecule_counts:
        atom_count = len(topology.atoms[name])
        for copy_index in range(count):
            copy_offset = offset + copy_index * atom_count
            global_indices.extend(
                copy_offset + local_index - 1 for local_index in perturbed_local[name]
            )
        offset += count * atom_count
    return sorted(global_indices)


def _global_atom_mapping(topology):
    mapping = []
    for name, count in topology.molecule_counts:
        atom_count = len(topology.atoms[name])
        for _copy_index in range(count):
            mapping.extend((name, index) for index in range(1, atom_count + 1))
    return mapping


def automatic_hot_region(
    topology_text,
    structure,
    *,
    distance=0.4,
    region_selection="not water",
):
    import mdtraj as md

    topology = parse_preprocessed_topology(topology_text)
    perturbed = perturbed_global_indices(topology)
    if not perturbed:
        raise ValueError("No perturbed atoms were found in the dual-state topology")

    frame = md.load(str(structure))
    mapping = _global_atom_mapping(topology)
    if frame.n_atoms != len(mapping):
        raise ValueError(
            "Structure and topology atom counts differ: "
            f"{frame.n_atoms} != {len(mapping)}"
        )
    haystack = frame.topology.select(region_selection)
    if len(haystack) == 0:
        raise ValueError("The REST hot-region selection is empty")
    neighbors = md.compute_neighbors(
        frame,
        distance,
        query_indices=perturbed,
        haystack_indices=haystack,
    )[0]
    selected_residues = {
        frame.topology.atom(int(index)).residue for index in [*neighbors, *perturbed]
    }
    selected_global = {
        atom.index for residue in selected_residues for atom in residue.atoms
    }
    selected_global.update(perturbed)
    return selected_global


def underline_global_atoms(topology_text, selected_global):
    topology = parse_preprocessed_topology(topology_text)
    mapping = _global_atom_mapping(topology)
    if any(index < 0 or index >= len(mapping) for index in selected_global):
        raise ValueError("REST hot atom index is outside the topology")
    occurrence_counts = {}
    selected_counts = {}
    for name, local_index in mapping:
        key = (name, local_index)
        occurrence_counts[key] = occurrence_counts.get(key, 0) + 1
    for global_index in selected_global:
        key = mapping[global_index]
        selected_counts[key] = selected_counts.get(key, 0) + 1
    for key, count in selected_counts.items():
        if count != occurrence_counts[key]:
            raise ValueError(
                "REST cannot select only some copies of a repeated molecule type: "
                f"{key[0]} atom {key[1]}"
            )
    selected_local = {mapping[index] for index in selected_global}
    lines = topology_text.splitlines(keepends=True)
    for name, records in topology.atoms.items():
        for record in records:
            if (name, record.index) not in selected_local:
                continue
            fields = _line_fields(lines[record.line_index])
            if not fields[1].endswith("_"):
                fields[1] += "_"
            lines[record.line_index] = _replace_fields(
                lines[record.line_index],
                fields,
            )
    return "".join(lines)


def hot_global_indices(topology_text):
    topology = parse_preprocessed_topology(topology_text)
    hot_local = {
        name: {record.index for record in records if record.type_a.endswith("_")}
        for name, records in topology.atoms.items()
    }
    output = []
    offset = 0
    for name, count in topology.molecule_counts:
        atom_count = len(topology.atoms[name])
        for copy_index in range(count):
            copy_offset = offset + copy_index * atom_count
            output.extend(copy_offset + index - 1 for index in hot_local[name])
        offset += count * atom_count
    return sorted(output)


def unify_fep_rest_charges(
    topology_text,
    *,
    general_lambda,
    charge_a_lambda,
    charge_b_lambda,
):
    topology = parse_preprocessed_topology(topology_text)
    lines = topology_text.splitlines(keepends=True)
    for records in topology.atoms.values():
        for record in records:
            if not record.has_b_state:
                continue
            type_a = record.type_a.removesuffix("_")
            type_b = record.type_b.removesuffix("_")
            epsilon_a = topology.epsilon_by_type.get(type_a)
            epsilon_b = topology.epsilon_by_type.get(type_b)
            if epsilon_a is None or epsilon_b is None:
                raise ValueError(
                    f"Missing atom type parameters for {type_a} or {type_b}"
                )
            phantom_a = (
                type_a != type_b
                and math.isclose(record.charge_a, 0.0, abs_tol=1e-12)
                and math.isclose(epsilon_a, 0.0, abs_tol=1e-12)
            )
            phantom_b = (
                type_a != type_b
                and math.isclose(record.charge_b, 0.0, abs_tol=1e-12)
                and math.isclose(epsilon_b, 0.0, abs_tol=1e-12)
            )
            if not phantom_a and phantom_b:
                charge_lambda = charge_a_lambda
            elif phantom_a and not phantom_b:
                charge_lambda = charge_b_lambda
            else:
                charge_lambda = general_lambda
            charge = (
                record.charge_a * (1.0 - charge_lambda)
                + record.charge_b * charge_lambda
            )
            mass = max(record.mass_a, record.mass_b)
            fields = _line_fields(lines[record.line_index])
            fields[6] = f"{charge:.10g}"
            fields[7] = f"{mass:.10g}"
            fields[9] = f"{charge:.10g}"
            fields[10] = f"{mass:.10g}"
            lines[record.line_index] = _replace_fields(
                lines[record.line_index],
                fields,
            )
    return "".join(lines)


_DIHEDRAL_PARAMETER_COUNTS = {
    1: (3, {1}),
    2: (2, {1}),
    3: (6, {0, 1, 2, 3, 4, 5}),
    4: (3, {1}),
    5: (4, {0, 1, 2, 3}),
    9: (3, {1}),
}


def prepare_plumed_dual_state(topology_text):
    for raw_line in topology_text.splitlines():
        fields = _line_fields(raw_line)
        if (
            len(fields) >= 3
            and fields[0] == "["
            and fields[-1] == "]"
            and fields[1].lower()
            in {"cmap", "cmaptypes", "intermolecular_interactions"}
        ):
            raise ValueError(
                "PLUMED partial_tempering does not support CMAP or "
                "intermolecular-interaction topologies"
            )

    topology = parse_preprocessed_topology(topology_text)
    lines = topology_text.splitlines(keepends=True)
    records = {}
    section = None
    molecule_type = None
    atom_hot = {}
    counter = 0

    for line_index, raw_line in enumerate(lines):
        fields = _line_fields(raw_line)
        if not fields:
            continue
        if fields[0] == "[" and fields[-1] == "]":
            section = fields[1].lower()
            if section in {"system", "molecules"}:
                molecule_type = None
            continue
        if section == "moleculetype":
            molecule_type = fields[0]
            atom_hot = {
                record.index: record.type_a.endswith("_")
                for record in topology.atoms[molecule_type]
            }
            section = None
            continue

        if section == "atoms" and len(fields) == 11 and fields[1].endswith("_"):
            tag = f"MDTBX_DUAL_ATOM_{counter}"
            records[tag] = {
                "kind": "atom",
                "type_b": fields[8].removesuffix("_"),
                "charge_b": float(fields[9]),
                "mass_b": float(fields[10]),
            }
            lines[line_index] = raw_line.rstrip("\n") + f" ; {tag}\n"
            counter += 1
        elif section == "dihedrals" and len(fields) > 5:
            function = int(fields[4])
            definition = _DIHEDRAL_PARAMETER_COUNTS.get(function)
            if definition is None:
                continue
            count, scaled_indices = definition
            parameters = fields[5:]
            if len(parameters) != 2 * count:
                continue
            tag = f"MDTBX_DUAL_DIHEDRAL_{counter}"
            hot_ends = int(atom_hot.get(int(fields[0]), False)) + int(
                atom_hot.get(int(fields[3]), False)
            )
            records[tag] = {
                "kind": "dihedral",
                "parameters_b": parameters[count:],
                "scaled_indices": scaled_indices,
                "hot_ends": hot_ends,
            }
            compatible_fields = fields[:5] + parameters[:count]
            lines[line_index] = (
                _replace_fields(raw_line, compatible_fields).rstrip("\n")
                + f" ; {tag}\n"
            )
            counter += 1
        elif section == "pairs" and len(fields) > 5 and int(fields[2]) == 1:
            raise ValueError(
                "PLUMED FEP/REST does not support explicit B-state pair parameters"
            )

    return "".join(lines), records


def restore_plumed_dual_state(topology_text, records, scale):
    output = []
    sqrt_scale = math.sqrt(scale)
    found = set()
    for raw_line in topology_text.splitlines(keepends=True):
        match = re.search(r"MDTBX_DUAL_(?:ATOM|DIHEDRAL)_\d+", raw_line)
        tag = match.group(0) if match else None
        if tag is None:
            output.append(raw_line)
            continue
        if tag not in records:
            raise ValueError(f"Unexpected PLUMED dual-state marker: {tag}")
        found.add(tag)
        record = records[tag]
        code, separator, comment = raw_line.partition(";")
        fields = code.split()
        cleaned_comment = comment.replace(tag, "").strip()
        if record["kind"] == "atom":
            fields = fields[:8]
            fields.extend(
                [
                    record["type_b"] + "_",
                    f"{record['charge_b'] * sqrt_scale:.10g}",
                    f"{record['mass_b']:.10g}",
                ]
            )
        else:
            dihedral_scale = sqrt_scale ** record["hot_ends"]
            parameters = []
            for index, value in enumerate(record["parameters_b"]):
                if index in record["scaled_indices"]:
                    parameters.append(f"{float(value) * dihedral_scale:.10g}")
                else:
                    parameters.append(value)
            fields.extend(parameters)
        restored = " ".join(fields)
        if separator and cleaned_comment:
            restored += f" ; {cleaned_comment}"
        output.append(restored + "\n")
    missing = set(records) - found
    if missing:
        raise ValueError(
            f"PLUMED output omitted dual-state records: {sorted(missing)[0]}"
        )
    return "".join(output)
