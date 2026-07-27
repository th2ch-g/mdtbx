import re
from argparse import _SubParsersAction
from pathlib import Path

from mdtbx.cli import create_parser


def _registered_commands() -> set[str]:
    parser = create_parser()
    subparsers = next(
        action for action in parser._actions if isinstance(action, _SubParsersAction)
    )
    return set(subparsers.choices)


def _documented_commands() -> set[str]:
    pattern = re.compile(r"^\s*:start_command:\s+(\S+)\s*$", re.MULTILINE)
    documented = set()
    for path in Path("docs/reference").glob("*.rst"):
        documented.update(pattern.findall(path.read_text()))
    return documented


def test_all_commands_are_in_reference():
    assert _documented_commands() == _registered_commands()
