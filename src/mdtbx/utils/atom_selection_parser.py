import fnmatch
from collections.abc import Mapping
from typing import Protocol, runtime_checkable

AtomRecord = Mapping[str, str | int]
NumberSelection = list[int] | range

PROTEIN_RESNAMES: set[str] = {
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
    # Common Amber protonation and disulfide variants.
    "ASH",
    "CYM",
    "CYX",
    "GLH",
    "HID",
    "HIE",
    "HIP",
    "LYN",
}

WATER_RESNAMES: set[str] = {"HOH", "WAT", "SOL"}

BACKBONE_ATOM_NAMES: set[str] = {"N", "CA", "C", "O", "OXT"}

ION_NAMES: set[str] = {
    "NA",
    "NA+",
    "CL",
    "CL-",
    "K",
    "K+",
    "MG",
    "MG2+",
    "ZN",
    "ZN2+",
    "CA",
    "CA2+",
}


@runtime_checkable
class SelectionNode(Protocol):
    def eval(self, mol: AtomRecord) -> bool: ...

    def __repr__(self) -> str:
        attrs = ", ".join(f"{k}={v!r}" for k, v in self.__dict__.items())
        return f"{self.__class__.__name__}({attrs})"

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, self.__class__):
            return NotImplemented
        return self.__dict__ == other.__dict__


class Not(SelectionNode):
    def __init__(self, selection: SelectionNode):
        self.selection = selection

    def eval(self, mol: AtomRecord) -> bool:
        return not self.selection.eval(mol)


class And(SelectionNode):
    def __init__(self, selections: list[SelectionNode]):
        self.selections = selections

    def eval(self, mol: AtomRecord) -> bool:
        if not self.selections:
            return True
        return all(s.eval(mol) for s in self.selections)


class Or(SelectionNode):
    def __init__(self, selections: list[SelectionNode]):
        self.selections = selections

    def eval(self, mol: AtomRecord) -> bool:
        if not self.selections:
            return False
        return any(s.eval(mol) for s in self.selections)


class Bracket(SelectionNode):
    def __init__(self, selection: SelectionNode):
        self.selection = selection

    def eval(self, mol: AtomRecord) -> bool:
        return self.selection.eval(mol)


class Chain(SelectionNode):
    def __init__(self, names: list[str]):
        self.names = names

    def eval(self, mol: AtomRecord) -> bool:
        chain = mol.get("chain")
        if not isinstance(chain, str):
            return False
        return any(fnmatch.fnmatchcase(chain, pattern) for pattern in self.names)


class Moleculetype(SelectionNode):
    def __init__(self, names: list[str]):
        self.names = names

    def eval(self, mol: AtomRecord) -> bool:
        moleculetype = mol.get("moleculetype")
        if not isinstance(moleculetype, str):
            return False
        return any(fnmatch.fnmatchcase(moleculetype, pattern) for pattern in self.names)


class ResName(SelectionNode):
    def __init__(self, names: list[str]):
        self.names = names

    def eval(self, mol: AtomRecord) -> bool:
        resname = mol.get("resname")
        if not isinstance(resname, str):
            return False
        return any(fnmatch.fnmatchcase(resname, pattern) for pattern in self.names)


class ResId(SelectionNode):
    def __init__(self, ids: NumberSelection):
        self.ids = ids

    def eval(self, mol: AtomRecord) -> bool:
        resid = mol.get("resid")
        return isinstance(resid, int) and resid in self.ids


class Name(SelectionNode):
    def __init__(self, names: list[str]):
        self.names = names

    def eval(self, mol: AtomRecord) -> bool:
        name = mol.get("name")
        if not isinstance(name, str):
            return False
        return any(fnmatch.fnmatchcase(name, pattern) for pattern in self.names)


class Index(SelectionNode):
    def __init__(self, indices: NumberSelection):
        self.indices = indices

    def eval(self, mol: AtomRecord) -> bool:
        index = mol.get("index")
        return isinstance(index, int) and index in self.indices


class All(SelectionNode):
    def eval(self, mol: AtomRecord) -> bool:
        return True


class Protein(SelectionNode):
    def eval(self, mol: AtomRecord) -> bool:
        resname = mol.get("resname")
        return isinstance(resname, str) and resname.upper() in PROTEIN_RESNAMES


class Water(SelectionNode):
    def eval(self, mol: AtomRecord) -> bool:
        resname = mol.get("resname")
        return isinstance(resname, str) and resname.upper() in WATER_RESNAMES


class Ion(SelectionNode):
    def eval(self, mol: AtomRecord) -> bool:
        name = mol.get("name")
        resname = mol.get("resname")
        normalized_resname = resname.upper() if isinstance(resname, str) else None
        normalized_name = name.upper() if isinstance(name, str) else None
        if normalized_resname in ION_NAMES:
            return True
        # Fall back to the atom name only when the residue is not a known protein
        # residue, so a protein alpha-carbon (name "CA") is not misread as calcium.
        return (
            normalized_name in ION_NAMES and normalized_resname not in PROTEIN_RESNAMES
        )


class Backbone(SelectionNode):
    def eval(self, mol: AtomRecord) -> bool:
        is_protein = Protein().eval(mol)
        name = mol.get("name")
        return is_protein and isinstance(name, str) and name in BACKBONE_ATOM_NAMES


class Sidechain(SelectionNode):
    def eval(self, mol: AtomRecord) -> bool:
        is_protein = Protein().eval(mol)
        is_backbone = Backbone().eval(mol)
        return is_protein and not is_backbone


class ParseError(ValueError):
    pass


FORBIDDEN_IDENTIFIERS = {"and", "or", "not", "to", "moleculetype"}


class SelectionParser:
    def __init__(self, text: str):
        self.text = text
        self.pos = 0

    @staticmethod
    def _is_identifier_char(char: str) -> bool:
        return char.isalnum() or char in "*?'_+-"

    def _peek(self) -> str | None:
        return self.text[self.pos] if self.pos < len(self.text) else None

    def _consume_char(self, char: str):
        if self._peek() == char:
            self.pos += 1
            return char
        raise ParseError(
            f"Expected '{char}' at position {self.pos}, got '{self._peek()}'"
        )

    def _matches_tag(self, tag: str) -> bool:
        end = self.pos + len(tag)
        return self.text.startswith(tag, self.pos) and (
            end == len(self.text) or not self._is_identifier_char(self.text[end])
        )

    def _consume_tag(self, tag: str):
        if self._matches_tag(tag):
            self.pos += len(tag)
            return tag
        raise ParseError(f"Expected '{tag}' at position {self.pos}")

    def _skip_space0(self):
        while self.pos < len(self.text) and self.text[self.pos].isspace():
            self.pos += 1

    def _skip_space1(self):
        start_pos = self.pos
        self._skip_space0()
        if self.pos == start_pos:
            raise ParseError(f"Expected one or more spaces at position {self.pos}")

    def _parse_alphanumeric1(self) -> str:
        """Parse an atom identifier or shell-style wildcard pattern."""
        start_pos = self.pos

        if self.pos < len(self.text) and self._is_identifier_char(self.text[self.pos]):
            self.pos += 1
            while self.pos < len(self.text) and self._is_identifier_char(
                self.text[self.pos]
            ):
                self.pos += 1
            return self.text[start_pos : self.pos]
        raise ParseError(
            f"Expected an atom identifier or wildcard at position {self.pos}"
        )

    def _parse_digit1(self) -> str:
        start_pos = self.pos
        if self.pos < len(self.text) and self.text[self.pos].isdigit():
            self.pos += 1
            while self.pos < len(self.text) and self.text[self.pos].isdigit():
                self.pos += 1
            return self.text[start_pos : self.pos]
        raise ParseError(f"Expected digits at position {self.pos}")

    def _parse_usize(self) -> int:
        s = self._parse_digit1()
        try:
            return int(s)
        except ValueError as exc:
            raise ParseError(f"Invalid unsigned integer: {s}") from exc

    def _parse_identifier(self) -> str:
        identifier = self._parse_alphanumeric1()
        if identifier.lower() in FORBIDDEN_IDENTIFIERS:
            raise ParseError(
                f"Identifier cannot be a reserved keyword: '{identifier}' at position {self.pos - len(identifier)}"
            )
        return identifier

    def _parse_list_of_identifiers(self) -> list[str]:
        identifiers = [self._parse_identifier()]
        while True:
            saved_pos = self.pos
            try:
                self._skip_space1()
                identifiers.append(self._parse_identifier())
            except ParseError:
                self.pos = saved_pos
                break
        return identifiers

    def _parse_numbers(self) -> NumberSelection:
        first = self._parse_usize()
        after_first = self.pos
        try:
            self._skip_space1()
        except ParseError:
            return [first]

        if self._matches_tag("to"):
            self._consume_tag("to")
            self._skip_space1()
            last = self._parse_usize()
            if last < first:
                raise ParseError(f"Range end {last} is less than start {first}")
            return range(first, last + 1)

        self.pos = after_first
        numbers = [first]
        while True:
            saved_pos = self.pos
            try:
                self._skip_space1()
                numbers.append(self._parse_usize())
            except ParseError:
                self.pos = saved_pos
                break
        return numbers

    def _parse_resname(self) -> SelectionNode:
        self._consume_tag("resname")
        self._skip_space1()
        return ResName(self._parse_list_of_identifiers())

    def _parse_name(self) -> SelectionNode:
        self._consume_tag("name")
        self._skip_space1()
        return Name(self._parse_list_of_identifiers())

    def _parse_resid(self) -> SelectionNode:
        self._consume_tag("resid")
        self._skip_space1()
        return ResId(self._parse_numbers())

    def _parse_index(self) -> SelectionNode:
        self._consume_tag("index")
        self._skip_space1()
        return Index(self._parse_numbers())

    def _parse_chain(self) -> SelectionNode:
        self._consume_tag("chain")
        self._skip_space1()
        return Chain(self._parse_list_of_identifiers())

    def _parse_moleculetype(self) -> SelectionNode:
        self._consume_tag("moleculetype")
        self._skip_space1()
        return Moleculetype(self._parse_list_of_identifiers())

    def _parse_all(self) -> SelectionNode:
        self._consume_tag("all")
        return All()

    def _parse_protein(self) -> SelectionNode:
        self._consume_tag("protein")
        return Protein()

    def _parse_water(self) -> SelectionNode:
        self._consume_tag("water")
        return Water()

    def _parse_ion(self) -> SelectionNode:
        self._consume_tag("ion")
        return Ion()

    def _parse_backbone(self) -> SelectionNode:
        self._consume_tag("backbone")
        return Backbone()

    def _parse_sidechain(self) -> SelectionNode:
        self._consume_tag("sidechain")
        return Sidechain()

    def _parse_atom(self) -> SelectionNode:
        atom_parsers = {
            "all": self._parse_all,
            "protein": self._parse_protein,
            "water": self._parse_water,
            "ion": self._parse_ion,
            "backbone": self._parse_backbone,
            "sidechain": self._parse_sidechain,
            "chain": self._parse_chain,
            "moleculetype": self._parse_moleculetype,
            "resname": self._parse_resname,
            "resid": self._parse_resid,
            "index": self._parse_index,
            "name": self._parse_name,
        }
        for tag, parser_func in atom_parsers.items():
            if self._matches_tag(tag):
                return parser_func()
        raise ParseError(
            f"Expected an atomic selection (e.g., 'all', 'protein', 'index 1') at position {self.pos}"
        )

    def _parse_bracket(self) -> SelectionNode:
        self._consume_char("(")
        self._skip_space0()
        expr = self.parse_expr()
        self._skip_space0()
        self._consume_char(")")
        return Bracket(expr)

    def _parse_primary(self) -> SelectionNode:
        self._skip_space0()
        saved_pos = self.pos
        try:
            return self._parse_bracket()
        except ParseError:
            self.pos = saved_pos
            try:
                return self._parse_atom()
            except ParseError as e_atom:
                if self.text[saved_pos:].startswith("("):
                    raise ParseError(
                        f"Syntax error in parenthesized expression or mismatched parentheses near position {saved_pos}"
                    ) from e_atom
                raise

    def _parse_not(self) -> SelectionNode:
        num_nots = 0
        while True:
            saved_pos = self.pos
            self._skip_space0()
            try:
                self._consume_tag("not")
                num_nots += 1
                self._skip_space1()
            except ParseError:
                self.pos = saved_pos
                break

        selection = self._parse_primary()

        for _ in range(num_nots):
            selection = Not(selection)
        return selection

    def _parse_and(self) -> SelectionNode:
        operands = [self._parse_not()]
        while True:
            saved_pos = self.pos
            try:
                self._skip_space1()
                self._consume_tag("and")
                self._skip_space1()
                operands.append(self._parse_not())
            except ParseError:
                self.pos = saved_pos
                break
        return operands[0] if len(operands) == 1 else And(operands)

    def _parse_or(self) -> SelectionNode:
        operands = [self._parse_and()]
        while True:
            saved_pos = self.pos
            try:
                self._skip_space1()
                self._consume_tag("or")
                self._skip_space1()
                operands.append(self._parse_and())
            except ParseError:
                self.pos = saved_pos
                break
        return operands[0] if len(operands) == 1 else Or(operands)

    def parse_expr(self) -> SelectionNode:
        self._skip_space0()
        expr = self._parse_or()
        self._skip_space0()
        return expr

    def parse(self) -> SelectionNode:
        parsed_node = self.parse_expr()
        if self.pos < len(self.text):
            raise ParseError(
                f"Unexpected trailing characters: '{self.text[self.pos :]}' at position {self.pos}"
            )
        return parsed_node


def parse_selection(selection_string: str) -> SelectionNode | str:
    try:
        parser = SelectionParser(selection_string)
        return parser.parse()
    except ParseError as exc:
        return str(exc)


class AtomSelector:
    def __init__(self, selection_string: str) -> None:
        self.selection_string = selection_string
        self.parsed_selection: SelectionNode | None = None
        self._error: str | None = None

        result = parse_selection(selection_string)
        if isinstance(result, SelectionNode):
            self.parsed_selection = result
        else:
            self._error = result
            raise ValueError(f"Failed to parse selection string: {self._error}")

    def eval(self, mol: AtomRecord) -> bool:
        if self.parsed_selection is None:
            raise RuntimeError("Selection was not successfully parsed.")

        return self.parsed_selection.eval(mol)
