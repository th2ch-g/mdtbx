from dataclasses import dataclass


@dataclass
class MoleculeType:
    name: str
    atoms: list
    start_line: int
    end_line: int


@dataclass
class GromacsTopologyParser:
    def __init__(self, topology_file: str):
        # main section: defaults, atomtypes,moleculetype, system
        # sub section : atoms. bonds, ..., dihedrals, exclusions, molecules

        all_moleculetypes = []
        moleculetype_dict = {}
        current_section = None
        start_line_section = None
        current_moleculetype = None
        with open(topology_file) as f:
            for idx, line in enumerate(f):
                line = line.strip()

                if line == "":
                    continue

                if line.startswith(";"):
                    continue

                if line.startswith("#"):
                    continue

                if line.startswith("[") and line.endswith("]"):
                    section = line[1:-1].strip()
                    current_section = section
                    start_line_section = idx
                    continue

                if current_section == "moleculetype":
                    moleculetype_name = line.split()[0]
                    all_moleculetypes.append(moleculetype_name)
                    moleculetype_dict[moleculetype_name] = MoleculeType(
                        moleculetype_name, [], start_line_section, None
                    )
                    if current_moleculetype is not None:
                        moleculetype_dict[current_moleculetype].end_line = (
                            start_line_section - 1
                        )
                    current_moleculetype = moleculetype_name

                if current_section == "atoms":
                    splits = line.split()
                    atom_index = int(splits[0])
                    atom_type = splits[1]
                    resid = int(splits[2])
                    resname = splits[3]
                    atom_name = splits[4]
                    # cgnr = splits[5]
                    # atom_charge = splits[6]
                    # atom_mass = splits[7]

                    atom = {
                        "atom_type": atom_type,
                        "index": atom_index,
                        "resid": resid,
                        "resname": resname,
                        "name": atom_name,
                        "linenumber": idx,
                    }

                    moleculetype_dict[current_moleculetype].atoms.append(atom)
        moleculetype_dict[current_moleculetype].end_line = idx

        self.topology_file = topology_file
        self.all_moleculetypes = all_moleculetypes
        self.moleculetype_dict = moleculetype_dict

    def get_all_moleculetypes(self) -> list[str]:
        return self.all_moleculetypes

    def get_atoms_in(self, moleculetypes: str) -> list[dict[str, str]]:
        return self.moleculetype_dict[moleculetypes].atoms

    def get_insert_linenumber_in(self, moleculetypes: str) -> int:
        return self.moleculetype_dict[moleculetypes].end_line
