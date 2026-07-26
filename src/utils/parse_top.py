from dataclasses import dataclass
from pathlib import Path


@dataclass
class MoleculeType:
    name: str
    atoms: list[dict[str, str | int]]
    start_line: int
    end_line: int | None


class GromacsTopologyParser:
    def __init__(self, topology_file: str):
        all_moleculetypes: list[str] = []
        moleculetype_dict: dict[str, MoleculeType] = {}
        current_section: str | None = None
        start_line_section: int | None = None
        current_moleculetype: str | None = None
        last_line = -1
        topology_path = Path(topology_file)

        with topology_path.open() as f:
            for idx, raw_line in enumerate(f):
                last_line = idx
                line = raw_line.partition(";")[0].strip()

                if line == "":
                    continue

                if line.startswith("#"):
                    continue

                if line.startswith("[") and line.endswith("]"):
                    section = line[1:-1].strip().lower()
                    if (
                        section in {"moleculetype", "system", "molecules"}
                        and current_moleculetype is not None
                        and moleculetype_dict[current_moleculetype].end_line is None
                    ):
                        moleculetype_dict[current_moleculetype].end_line = idx - 1
                    if section in {"moleculetype", "system", "molecules"}:
                        current_moleculetype = None
                    if section == "atoms" and current_moleculetype is None:
                        raise ValueError(
                            f"Atoms section has no preceding moleculetype at "
                            f"{topology_path}:{idx + 1}"
                        )
                    current_section = section
                    start_line_section = idx
                    continue

                if current_section == "moleculetype":
                    moleculetype_name = line.split()[0]
                    if moleculetype_name in moleculetype_dict:
                        raise ValueError(
                            f"Duplicate moleculetype '{moleculetype_name}' "
                            f"at {topology_path}:{idx + 1}"
                        )
                    if start_line_section is None:
                        raise RuntimeError("Missing moleculetype section start")
                    all_moleculetypes.append(moleculetype_name)
                    moleculetype_dict[moleculetype_name] = MoleculeType(
                        moleculetype_name, [], start_line_section, None
                    )
                    current_moleculetype = moleculetype_name
                    current_section = None
                    continue

                if current_section == "atoms":
                    if current_moleculetype is None:
                        continue
                    splits = line.split()
                    if len(splits) < 5:
                        raise ValueError(
                            f"Malformed atom record at {topology_path}:{idx + 1}"
                        )
                    try:
                        atom_index = int(splits[0])
                        resid = int(splits[2])
                    except ValueError as exc:
                        raise ValueError(
                            f"Invalid atom index or residue ID at "
                            f"{topology_path}:{idx + 1}"
                        ) from exc
                    atom_type = splits[1]
                    resname = splits[3]
                    atom_name = splits[4]

                    atom = {
                        "atom_type": atom_type,
                        "index": atom_index,
                        "resid": resid,
                        "resname": resname,
                        "name": atom_name,
                        "linenumber": idx,
                    }

                    moleculetype_dict[current_moleculetype].atoms.append(atom)
        if current_section == "moleculetype":
            raise ValueError(
                f"Moleculetype declaration is missing at "
                f"{topology_path}:{last_line + 1}"
            )
        if (
            current_moleculetype is not None
            and moleculetype_dict[current_moleculetype].end_line is None
        ):
            moleculetype_dict[current_moleculetype].end_line = last_line

        self.topology_file = topology_file
        self.all_moleculetypes = all_moleculetypes
        self.moleculetype_dict = moleculetype_dict

    def get_all_moleculetypes(self) -> list[str]:
        return self.all_moleculetypes

    def get_atoms_in(self, moleculetypes: str) -> list[dict[str, str | int]]:
        return self.moleculetype_dict[moleculetypes].atoms

    def get_insert_linenumber_in(self, moleculetypes: str) -> int:
        end_line = self.moleculetype_dict[moleculetypes].end_line
        if end_line is None:
            raise ValueError(f"Moleculetype '{moleculetypes}' has no complete block")
        return end_line
