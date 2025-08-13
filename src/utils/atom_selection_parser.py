# from typing import List, Union, Dict, Protocol, runtime_checkable
#
#
# @runtime_checkable
# class SelectionNode(Protocol):
#     def eval(self, mol: Dict[str, Union[str, int]]) -> bool: ...
#
#     def __repr__(self) -> str:
#         attrs = ", ".join(f"{k}={v!r}" for k, v in self.__dict__.items())
#         return f"{self.__class__.__name__}({attrs})"
#
#     def __eq__(self, other: object) -> bool:
#         if not isinstance(other, self.__class__):
#             return NotImplemented
#         return self.__dict__ == other.__dict__
#
#
# class Chain(SelectionNode):
#     def __init__(self, names: List[str]):
#         self.names = names
#
#     def eval(self, mol: Dict[str, Union[str, int]]) -> bool:
#         chain = mol.get("chain")
#         return isinstance(chain, str) and chain in self.names
#
#
# class ResName(SelectionNode):
#     def __init__(self, names: List[str]):
#         self.names = names
#
#     def eval(self, mol: Dict[str, Union[str, int]]) -> bool:
#         resname = mol.get("resname")
#         return isinstance(resname, str) and resname in self.names
#
#
# class ResId(SelectionNode):
#     def __init__(self, ids: List[int]):
#         self.ids = ids
#
#     def eval(self, mol: Dict[str, Union[str, int]]) -> bool:
#         resid = mol.get("resid")
#         return isinstance(resid, int) and resid in self.ids
#
#
# class Name(SelectionNode):
#     def __init__(self, names: List[str]):
#         self.names = names
#
#     def eval(self, mol: Dict[str, Union[str, int]]) -> bool:
#         name = mol.get("name")
#         return isinstance(name, str) and name in self.names
#
#
# class Index(SelectionNode):
#     def __init__(self, indices: List[int]):
#         self.indices = indices
#
#     def eval(self, mol: Dict[str, Union[str, int]]) -> bool:
#         index = mol.get("index")
#         return isinstance(index, int) and index in self.indices
#
#
# class Not(SelectionNode):
#     def __init__(self, selection: SelectionNode):
#         self.selection = selection
#
#     def eval(self, mol: Dict[str, Union[str, int]]) -> bool:
#         return not self.selection.eval(mol)
#
#
# class And(SelectionNode):
#     def __init__(self, selections: List[SelectionNode]):
#         self.selections = selections
#
#     def eval(self, mol: Dict[str, Union[str, int]]) -> bool:
#         if not self.selections:
#             return True
#         return all(s.eval(mol) for s in self.selections)
#
#
# class Or(SelectionNode):
#     def __init__(self, selections: List[SelectionNode]):
#         self.selections = selections
#
#     def eval(self, mol: Dict[str, Union[str, int]]) -> bool:
#         if not self.selections:
#             return False
#         return any(s.eval(mol) for s in self.selections)
#
#
# class Bracket(SelectionNode):
#     def __init__(self, selection: SelectionNode):
#         self.selection = selection
#
#     def eval(self, mol: Dict[str, Union[str, int]]) -> bool:
#         return self.selection.eval(mol)
#
#
# # --- Parser Error ---
# class ParseError(ValueError):
#     pass
#
#
# # --- Parser Class ---
# RESERVED_KEYWORDS = {
#     "and",
#     "or",
#     "not",
#     "to",
#     "resname",
#     "resid",
#     "name",
#     "index",
#     "chain",
# }
#
#
# class SelectionParser:
#     def __init__(self, text: str):
#         self.text = text
#         self.pos = 0
#
#     def _peek(self) -> Union[str, None]:
#         return self.text[self.pos] if self.pos < len(self.text) else None
#
#     def _consume_char(self, char: str):
#         if self._peek() == char:
#             self.pos += 1
#             return char
#         raise ParseError(
#             f"Expected '{char}' at position {self.pos}, got '{self._peek()}'"
#         )
#
#     def _consume_tag(self, tag: str):
#         if self.text.startswith(tag, self.pos):
#             if tag.isalpha() and (
#                 self.pos + len(tag) < len(self.text)
#                 and self.text[self.pos + len(tag)].isalnum()
#             ):
#                 pass
#             self.pos += len(tag)
#             return tag
#         raise ParseError(f"Expected '{tag}' at position {self.pos}")
#
#     def _skip_space0(self):
#         while self.pos < len(self.text) and self.text[self.pos].isspace():
#             self.pos += 1
#
#     def _skip_space1(self):
#         start_pos = self.pos
#         self._skip_space0()
#         if self.pos == start_pos:
#             raise ParseError(f"Expected one or more spaces at position {self.pos}")
#
#     def _parse_alphanumeric1(self) -> str:
#         start_pos = self.pos
#         if self.pos < len(self.text) and self.text[self.pos].isalnum():
#             self.pos += 1
#             while self.pos < len(self.text) and self.text[self.pos].isalnum():
#                 self.pos += 1
#             return self.text[start_pos : self.pos]
#         raise ParseError(f"Expected alphanumeric characters at position {self.pos}")
#
#     def _parse_digit1(self) -> str:
#         start_pos = self.pos
#         if self.pos < len(self.text) and self.text[self.pos].isdigit():
#             self.pos += 1
#             while self.pos < len(self.text) and self.text[self.pos].isdigit():
#                 self.pos += 1
#             return self.text[start_pos : self.pos]
#         raise ParseError(f"Expected digits at position {self.pos}")
#
#     def _parse_usize(self) -> int:
#         s = self._parse_digit1()
#         try:
#             return int(s)
#         except ValueError:
#             raise ParseError(f"Invalid unsigned integer: {s}")
#
#     def _parse_identifier(self) -> str:
#         identifier = self._parse_alphanumeric1()
#         if identifier in RESERVED_KEYWORDS - {"resname", "resid", "name", "index"}:
#             if identifier in {"and", "or", "not", "to"}:
#                 raise ParseError(
#                     f"Identifier cannot be a reserved keyword: '{identifier}' at position {self.pos - len(identifier)}"
#                 )
#         return identifier
#
#     def _parse_list_of_identifiers(self) -> List[str]:
#         identifiers = [self._parse_identifier()]
#         while True:
#             saved_pos = self.pos
#             try:
#                 self._skip_space1()
#                 identifiers.append(self._parse_identifier())
#             except ParseError:
#                 self.pos = saved_pos
#                 break
#         return identifiers
#
#     def _parse_numbers(self) -> List[int]:
#         first = self._parse_usize()
#
#         saved_pos_for_to = self.pos
#         try:
#             self._skip_space1()
#             self._consume_tag("to")
#             self._skip_space1()
#             last = self._parse_usize()
#             if last < first:
#                 raise ParseError(f"Range end {last} is less than start {first}")
#             return list(range(first, last + 1))
#         except ParseError:
#             self.pos = saved_pos_for_to
#             numbers = [first]
#             while True:
#                 saved_pos_loop = self.pos
#                 try:
#                     self._skip_space1()
#                     numbers.append(self._parse_usize())
#                 except ParseError:
#                     self.pos = saved_pos_loop
#                     break
#             return numbers
#
#     def _parse_resname(self) -> SelectionNode:
#         self._consume_tag("resname")
#         self._skip_space1()
#         return ResName(self._parse_list_of_identifiers())
#
#     def _parse_name(self) -> SelectionNode:
#         self._consume_tag("name")
#         self._skip_space1()
#         return Name(self._parse_list_of_identifiers())
#
#     def _parse_resid(self) -> SelectionNode:
#         self._consume_tag("resid")
#         self._skip_space1()
#         return ResId(self._parse_numbers())
#
#     def _parse_index(self) -> SelectionNode:
#         self._consume_tag("index")
#         self._skip_space1()
#         return Index(self._parse_numbers())
#
#     def _parse_chain(self) -> SelectionNode:
#         self._consume_tag("chain")
#         self._skip_space1()
#         return Chain(self._parse_list_of_identifiers())
#
#     # --- Grammar hierarchy ---
#     def _parse_atom(self) -> SelectionNode:
#         atom_parsers = [
#             self._parse_chain,
#             self._parse_resname,
#             self._parse_resid,
#             self._parse_index,
#             self._parse_name,
#         ]
#         for parser_func in atom_parsers:
#             saved_pos = self.pos
#             try:
#                 return parser_func()
#             except ParseError:
#                 self.pos = saved_pos
#         raise ParseError(
#             f"Expected an atomic selection (e.g., 'all', 'index 1') at position {self.pos}"
#         )
#
#     def _parse_bracket(self) -> SelectionNode:
#         self._consume_char("(")
#         expr = self.parse_expr()
#         self._consume_char(")")
#         return Bracket(expr)
#
#     def _parse_primary(self) -> SelectionNode:
#         self._skip_space0()
#         saved_pos = self.pos
#         try:
#             return self._parse_bracket()
#         except ParseError:
#             self.pos = saved_pos
#             try:
#                 return self._parse_atom()
#             except ParseError as e_atom:
#                 if self.text[saved_pos:].startswith("("):
#                     raise ParseError(
#                         f"Syntax error in parenthesized expression or mismatched parentheses near position {saved_pos}"
#                     ) from e_atom
#                 raise
#
#     def _parse_not(self) -> SelectionNode:
#         num_nots = 0
#         while True:
#             saved_pos = self.pos
#             self._skip_space0()
#             try:
#                 self._consume_tag("not")
#                 num_nots += 1
#             except ParseError:
#                 self.pos = saved_pos
#                 break
#
#         selection = self._parse_primary()
#
#         for _ in range(num_nots):
#             selection = Not(selection)
#         return selection
#
#     def _parse_and(self) -> SelectionNode:
#         operands = [self._parse_not()]
#         while True:
#             saved_pos = self.pos
#             try:
#                 self._skip_space1()
#                 self._consume_tag("and")
#                 self._skip_space1()
#                 operands.append(self._parse_not())
#             except ParseError:
#                 self.pos = saved_pos
#                 break
#         return operands[0] if len(operands) == 1 else And(operands)
#
#     def _parse_or(self) -> SelectionNode:
#         operands = [self._parse_and()]
#         while True:
#             saved_pos = self.pos
#             try:
#                 self._skip_space1()
#                 self._consume_tag("or")
#                 self._skip_space1()
#                 operands.append(self._parse_and())
#             except ParseError:
#                 self.pos = saved_pos
#                 break
#         return operands[0] if len(operands) == 1 else Or(operands)
#
#     def parse_expr(self) -> SelectionNode:
#         expr = self._parse_or()
#         self._skip_space0()
#         return expr
#
#     def parse(self) -> SelectionNode:
#         parsed_node = self.parse_expr()
#         if self.pos < len(self.text):
#             raise ParseError(
#                 f"Unexpected trailing characters: '{self.text[self.pos :]}' at position {self.pos}"
#             )
#         return parsed_node
#
#
# # --- Public API as per Rust code ---
# def parse_selection(selection_string: str) -> Union[SelectionNode, str]:
#     try:
#         parser = SelectionParser(selection_string)
#         return parser.parse()
#     except ParseError as e:
#         return str(e)
#
#
# # --- AtomSelector Class (User's starting point) ---
# class AtomSelector:
#     def __init__(self, selection_string: str) -> None:
#         self.selection_string = selection_string
#         self.parsed_selection: Union[SelectionNode, None] = None
#         self._error: Union[str, None] = None
#
#         result = parse_selection(selection_string)
#         if isinstance(result, SelectionNode):
#             self.parsed_selection = result
#         else:
#             self._error = result
#             raise ValueError(f"Failed to parse selection string: {self._error}")
#
#     def eval(self, mol: Dict[str, Union[str, int]]) -> bool:
#         """
#         mol: Dict[str, str | int]
#         e.g.
#         # mol = { "chain": "A", "resname": "ALA", "resid": 1, "name": "CA", "index": 0 }
#         mol = { "chain": "A", "resid": 1, "index": 0 }
#         """
#         if self._error:
#             print(f"Cannot evaluate: parsing failed with error: {self._error}")
#             return False
#         if self.parsed_selection is None:
#             raise RuntimeError("Selection was not successfully parsed.")
#
#         return self.parsed_selection.eval(mol)
#
#
# if __name__ == "__main__":
#     atom_selection = "(not chain A) and resid 1 to 11"
#     print(f"{atom_selection=}")
#     atom_selector = AtomSelector(atom_selection)
#
#     test_molecules = [
#         ({"chain": "A", "resname": "ALA", "resid": 1, "name": "CA", "index": 0}, False),
#         ({"chain": "B", "resname": "ALA", "resid": 1, "name": "CA", "index": 0}, True),
#         (
#             {"chain": "A", "resname": "ALA", "resid": 11, "name": "CA", "index": 0},
#             False,
#         ),
#         ({"chain": "A", "resid": 1, "index": 0}, False),
#         ({"chain": "B", "resid": 1, "index": 0}, True),
#         ({"chain": "A", "resid": 11, "index": 0}, False),
#     ]
#
#     for mol, expected in test_molecules:
#         print(f"Mol: {mol}, Eval: {atom_selector.eval(mol)}, Expected: {expected}")
#         assert atom_selector.eval(mol) == expected


from typing import List, Union, Dict, Protocol, runtime_checkable, Set

# --- データ定義 ---

# タンパク質を構成する標準的なアミノ酸残基名
PROTEIN_RESNAMES: Set[str] = {
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
    # "CYX", "NLN",
}

# 一般的な水分子の残基名
WATER_RESNAMES: Set[str] = {"HOH", "WAT", "SOL"}

# タンパク質の主鎖を構成する原子名
BACKBONE_ATOM_NAMES: Set[str] = {"N", "CA", "C", "O", "OXT"}


# --- 選択ノードの定義 ---


@runtime_checkable
class SelectionNode(Protocol):
    """
    すべての選択ノードが準拠するプロトコル（インターフェース）。
    """

    def eval(self, mol: Dict[str, Union[str, int]]) -> bool: ...

    def __repr__(self) -> str:
        attrs = ", ".join(f"{k}={v!r}" for k, v in self.__dict__.items())
        return f"{self.__class__.__name__}({attrs})"

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, self.__class__):
            return NotImplemented
        return self.__dict__ == other.__dict__


# --- 論理演算子ノード ---


class Not(SelectionNode):
    def __init__(self, selection: SelectionNode):
        self.selection = selection

    def eval(self, mol: Dict[str, Union[str, int]]) -> bool:
        return not self.selection.eval(mol)


class And(SelectionNode):
    def __init__(self, selections: List[SelectionNode]):
        self.selections = selections

    def eval(self, mol: Dict[str, Union[str, int]]) -> bool:
        if not self.selections:
            return True
        return all(s.eval(mol) for s in self.selections)


class Or(SelectionNode):
    def __init__(self, selections: List[SelectionNode]):
        self.selections = selections

    def eval(self, mol: Dict[str, Union[str, int]]) -> bool:
        if not self.selections:
            return False
        return any(s.eval(mol) for s in self.selections)


class Bracket(SelectionNode):
    def __init__(self, selection: SelectionNode):
        self.selection = selection

    def eval(self, mol: Dict[str, Union[str, int]]) -> bool:
        return self.selection.eval(mol)


# --- キーワード引数を持つ選択ノード ---


class Chain(SelectionNode):
    def __init__(self, names: List[str]):
        self.names = names

    def eval(self, mol: Dict[str, Union[str, int]]) -> bool:
        chain = mol.get("chain")
        return isinstance(chain, str) and chain in self.names


class ResName(SelectionNode):
    def __init__(self, names: List[str]):
        self.names = names

    def eval(self, mol: Dict[str, Union[str, int]]) -> bool:
        resname = mol.get("resname")
        return isinstance(resname, str) and resname in self.names


class ResId(SelectionNode):
    def __init__(self, ids: List[int]):
        self.ids = ids

    def eval(self, mol: Dict[str, Union[str, int]]) -> bool:
        resid = mol.get("resid")
        return isinstance(resid, int) and resid in self.ids


class Name(SelectionNode):
    def __init__(self, names: List[str]):
        self.names = names

    def eval(self, mol: Dict[str, Union[str, int]]) -> bool:
        name = mol.get("name")
        return isinstance(name, str) and name in self.names


class Index(SelectionNode):
    def __init__(self, indices: List[int]):
        self.indices = indices

    def eval(self, mol: Dict[str, Union[str, int]]) -> bool:
        index = mol.get("index")
        return isinstance(index, int) and index in self.indices


# --- 引数なしのキーワード選択ノード (Rustコードから追加) ---


class All(SelectionNode):
    def eval(self, mol: Dict[str, Union[str, int]]) -> bool:
        return True


class Protein(SelectionNode):
    def eval(self, mol: Dict[str, Union[str, int]]) -> bool:
        resname = mol.get("resname")
        return isinstance(resname, str) and resname in PROTEIN_RESNAMES


class Water(SelectionNode):
    def eval(self, mol: Dict[str, Union[str, int]]) -> bool:
        resname = mol.get("resname")
        return isinstance(resname, str) and resname in WATER_RESNAMES


class Ion(SelectionNode):
    def eval(self, mol: Dict[str, Union[str, int]]) -> bool:
        # 簡単なイオンの判定ロジック（必要に応じて拡張が必要）
        name = mol.get("name")
        resname = mol.get("resname")
        # 例: Na+, Cl-, MG, ZN など
        if isinstance(name, str) and name in {"NA", "CL", "K", "MG", "ZN", "CA"}:
            return True
        if isinstance(resname, str) and resname in {"NA", "CL", "K", "MG", "ZN", "CA"}:
            return True
        return False


class Backbone(SelectionNode):
    def eval(self, mol: Dict[str, Union[str, int]]) -> bool:
        # タンパク質の主鎖原子か判定
        is_protein = Protein().eval(mol)
        name = mol.get("name")
        return is_protein and isinstance(name, str) and name in BACKBONE_ATOM_NAMES


class Sidechain(SelectionNode):
    def eval(self, mol: Dict[str, Union[str, int]]) -> bool:
        # タンパク質の一部であり、かつ主鎖でない原子を側鎖とみなす
        is_protein = Protein().eval(mol)
        is_backbone = Backbone().eval(mol)
        return is_protein and not is_backbone


# --- パーサーエラー ---
class ParseError(ValueError):
    pass


# --- パーサークラス ---
# 予約キーワード: and, or, not, to は識別子として使えない
FORBIDDEN_IDENTIFIERS = {"and", "or", "not", "to"}


class SelectionParser:
    def __init__(self, text: str):
        self.text = text
        self.pos = 0

    def _peek(self) -> Union[str, None]:
        return self.text[self.pos] if self.pos < len(self.text) else None

    def _consume_char(self, char: str):
        if self._peek() == char:
            self.pos += 1
            return char
        raise ParseError(
            f"Expected '{char}' at position {self.pos}, got '{self._peek()}'"
        )

    def _consume_tag(self, tag: str):
        if self.text.startswith(tag, self.pos):
            # タグが単語として完結しているかチェック (例: 'res' が 'residue' にマッチしないようにする)
            if (
                self.pos + len(tag) < len(self.text)
                and self.text[self.pos + len(tag)].isalnum()
            ):
                raise ParseError(
                    f"Expected '{tag}' but found a longer word at position {self.pos}"
                )
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
        start_pos = self.pos
        if self.pos < len(self.text) and self.text[self.pos].isalnum():
            self.pos += 1
            while self.pos < len(self.text) and self.text[self.pos].isalnum():
                self.pos += 1
            return self.text[start_pos : self.pos]
        raise ParseError(f"Expected alphanumeric characters at position {self.pos}")

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
        except ValueError:
            raise ParseError(f"Invalid unsigned integer: {s}")

    def _parse_identifier(self) -> str:
        identifier = self._parse_alphanumeric1()
        if identifier.lower() in FORBIDDEN_IDENTIFIERS:
            raise ParseError(
                f"Identifier cannot be a reserved keyword: '{identifier}' at position {self.pos - len(identifier)}"
            )
        return identifier

    def _parse_list_of_identifiers(self) -> List[str]:
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

    def _parse_numbers(self) -> List[int]:
        first = self._parse_usize()
        saved_pos_for_to = self.pos
        try:
            self._skip_space1()
            self._consume_tag("to")
            self._skip_space1()
            last = self._parse_usize()
            if last < first:
                raise ParseError(f"Range end {last} is less than start {first}")
            return list(range(first, last + 1))
        except ParseError:
            self.pos = saved_pos_for_to
            numbers = [first]
            while True:
                saved_pos_loop = self.pos
                try:
                    self._skip_space1()
                    numbers.append(self._parse_usize())
                except ParseError:
                    self.pos = saved_pos_loop
                    break
            return numbers

    # --- 引数付きキーワードのパーサー ---
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

    # --- 引数なしキーワードのパーサー (Rustコードから追加) ---
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

    # --- 文法階層 ---
    def _parse_atom(self) -> SelectionNode:
        atom_parsers = [
            # 引数なしのキーワードを先に試す
            self._parse_all,
            self._parse_protein,
            self._parse_water,
            self._parse_ion,
            self._parse_backbone,
            self._parse_sidechain,
            # 引数ありのキーワード
            self._parse_chain,
            self._parse_resname,
            self._parse_resid,
            self._parse_index,
            self._parse_name,
        ]
        for parser_func in atom_parsers:
            saved_pos = self.pos
            try:
                return parser_func()
            except ParseError:
                self.pos = saved_pos
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
                self._skip_space1()  # not の後にはスペースが必須
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


# --- 公開API ---
def parse_selection(selection_string: str) -> Union[SelectionNode, str]:
    try:
        parser = SelectionParser(selection_string)
        return parser.parse()
    except ParseError as e:
        return str(e)


# --- AtomSelector クラス ---
class AtomSelector:
    def __init__(self, selection_string: str) -> None:
        self.selection_string = selection_string
        self.parsed_selection: Union[SelectionNode, None] = None
        self._error: Union[str, None] = None

        result = parse_selection(selection_string)
        if isinstance(result, SelectionNode):
            self.parsed_selection = result
        else:
            self._error = result
            raise ValueError(f"Failed to parse selection string: {self._error}")

    def eval(self, mol: Dict[str, Union[str, int]]) -> bool:
        """
        与えられた分子情報 (mol) が選択条件に合致するかどうかを評価します。

        Args:
            mol: 分子情報を表す辞書。
                 例: {"chain": "A", "resname": "ALA", "resid": 1, "name": "CA", "index": 0}

        Returns:
            条件に合致すれば True, そうでなければ False。
        """
        if self._error:
            print(f"Cannot evaluate: parsing failed with error: {self._error}")
            return False
        if self.parsed_selection is None:
            raise RuntimeError("Selection was not successfully parsed.")

        return self.parsed_selection.eval(mol)


# --- 実行とテスト ---
if __name__ == "__main__":
    print("--- Basic Tests ---")
    # テストケース: (セレクション文字列, 分子情報, 期待される評価結果)
    test_cases = [
        ("protein and backbone", {"resname": "ALA", "name": "CA"}, True),
        ("protein and backbone", {"resname": "ALA", "name": "CB"}, False),
        ("protein and sidechain", {"resname": "ALA", "name": "CB"}, True),
        ("protein and sidechain", {"resname": "HOH", "name": "O"}, False),
        ("water", {"resname": "HOH", "name": "O"}, True),
        ("not water", {"resname": "HOH", "name": "O"}, False),
        ("resid 10 to 20 and name CA", {"resid": 15, "name": "CA"}, True),
        ("resid 10 to 20 and name CA", {"resid": 25, "name": "CA"}, False),
        (
            "(resname ALA LYS) or (resname GLY and name N)",
            {"resname": "LYS", "name": "CA"},
            True,
        ),
        (
            "(resname ALA LYS) or (resname GLY and name N)",
            {"resname": "GLY", "name": "N"},
            True,
        ),
        (
            "(resname ALA LYS) or (resname GLY and name N)",
            {"resname": "GLY", "name": "CA"},
            False,
        ),
        ("all", {"resname": "XYZ"}, True),
    ]

    for i, (selection_str, mol, expected) in enumerate(test_cases):
        print(f"\nTest {i + 1}:")
        print(f"  Selection: '{selection_str}'")
        print(f"  Molecule:  {mol}")
        try:
            selector = AtomSelector(selection_str)
            result = selector.eval(mol)
            print(f"  Parsed:    {selector.parsed_selection}")
            print(f"  Eval -> {result}, Expected -> {expected}")
            assert result == expected
            print("  Result: PASSED")
        except (ValueError, AssertionError) as e:
            print(f"  Result: FAILED - {e}")

    print("\n--- Error Handling Tests ---")
    error_cases = [
        "resid 10 to",  # 不完全な範囲指定
        "resname and",  # 予約語の誤用
        "protein andd name",  # 未定義の演算子
        "( resid 1",  # 括弧が閉じていない
        "name ca cb and",  # and の後に項がない
    ]
    for i, selection_str in enumerate(error_cases):
        print(f"\nError Test {i + 1}:")
        print(f"  Selection: '{selection_str}'")
        try:
            AtomSelector(selection_str)
            print("  Result: FAILED - Expected an error, but none was raised.")
        except ValueError as e:
            print(f"  Result: PASSED - Correctly caught error: {e}")
