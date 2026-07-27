"""Build machine-readable command descriptors from argparse."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

from .model import json_value

AGENT_COMMANDS = {
    "agent_collect",
    "agent_plan",
    "agent_probe",
    "agent_profile_save",
    "agent_run",
    "agent_schema",
    "agent_status",
}

MD_COMMANDS = {"run_fep", "run_abfe", "analyze_fep_rest", "opt_perf"}
GPU_COMMANDS = {"modeling_cf"}
QUANTUM_COMMANDS = {
    "gen_am1bcc",
    "gen_resp",
    "gen_modres_am1bcc",
    "gen_modres_resp",
    "place_solvent",
}
DATA_COMMANDS = {
    "fit",
    "trjcat",
    "pacs_trjcat",
    "print_perf",
    "analyze_fep",
    "analyze_abfe",
    "extract_str",
    "extract_ave_str",
    "contactmap",
    "comdist",
    "comvec",
    "densmap",
    "distmap",
    "mindist",
    "pca",
    "rmsd",
    "rmsf",
    "xyz",
    "tica",
    "cluster",
    "msm",
}
READ_ONLY_COMMANDS = {"show_mdtraj", "show_npy", "print_perf", "shell_hook"}
PILOT_COMMANDS = {"run_fep", "run_abfe", "opt_perf"}

INPUT_DESTS = {
    "input",
    "input1",
    "input2",
    "structure",
    "structure1",
    "structure2",
    "topology",
    "trajectory",
    "trajectories",
    "mdp",
    "reference",
    "checkpoint",
    "index",
    "path",
    "npy",
    "mol2",
    "pdb",
    "pdb1",
    "pdb2",
    "file",
    "gro",
    "coord",
    "coordinates",
    "parm",
    "prmtop",
    "ref_structure",
    "complex_checkpoint",
    "complex_index",
    "complex_mdp",
    "complex_reference",
    "complex_structure",
    "complex_topology",
    "solvent_checkpoint",
    "solvent_index",
    "solvent_mdp",
    "solvent_reference",
    "solvent_structure",
    "solvent_topology",
    "template_tleap",
    "xvv",
}
OUTPUT_DESTS = {
    "output",
    "outdir",
    "output_dir",
    "output_prefix",
    "workdir",
    "output_average",
    "output_npz",
    "output_pdb",
    "history_output",
    "xvv_output",
    "trial_dir",
}


def command_class(name: str) -> str:
    if name in MD_COMMANDS:
        return "md"
    if name in GPU_COMMANDS:
        return "gpu"
    if name in QUANTUM_COMMANDS:
        return "quantum"
    if name in DATA_COMMANDS:
        return "data"
    return "light"


def command_risk(name: str) -> str:
    if name == "cmd":
        return "unsafe"
    if name == "rmfile":
        return "destructive"
    if name in READ_ONLY_COMMANDS:
        return "read"
    return "write"


def _action_type(action: argparse.Action) -> str:
    if isinstance(action, argparse._StoreTrueAction):
        return "boolean"
    if isinstance(action, argparse._StoreFalseAction):
        return "boolean"
    if action.type is None:
        return "string"
    return getattr(action.type, "__name__", str(action.type))


def _action_schema(action: argparse.Action) -> dict[str, Any]:
    positional = not action.option_strings
    return {
        "name": action.dest,
        "options": list(action.option_strings),
        "positional": positional,
        "required": bool(action.required or (positional and action.default is None)),
        "type": _action_type(action),
        "nargs": action.nargs,
        "choices": json_value(list(action.choices))
        if action.choices is not None
        else None,
        "default": json_value(action.default),
        "help": action.help,
    }


def descriptor(name: str, parser: argparse.ArgumentParser) -> dict[str, Any]:
    arguments = [
        _action_schema(action)
        for action in parser._actions
        if action.dest
        not in {
            "help",
            "func",
            "_agent_json",
            "_agent_dry_run",
            "_agent_cluster_profile",
        }
    ]
    return {
        "name": name,
        "description": parser.description,
        "resource_class": "internal" if name in AGENT_COMMANDS else command_class(name),
        "risk": "internal" if name in AGENT_COMMANDS else command_risk(name),
        "pilot_capable": name in PILOT_COMMANDS,
        "inputs": [item["name"] for item in arguments if item["name"] in INPUT_DESTS],
        "outputs": [item["name"] for item in arguments if item["name"] in OUTPUT_DESTS],
        "arguments": arguments,
    }


def all_descriptors(parser: argparse.ArgumentParser) -> dict[str, dict[str, Any]]:
    subparsers = next(
        action
        for action in parser._actions
        if isinstance(action, argparse._SubParsersAction)
    )
    return {
        name: descriptor(name, command_parser)
        for name, command_parser in sorted(subparsers.choices.items())
    }


def command_parser(
    parser: argparse.ArgumentParser, name: str
) -> argparse.ArgumentParser:
    subparsers = next(
        action
        for action in parser._actions
        if isinstance(action, argparse._SubParsersAction)
    )
    try:
        return subparsers.choices[name]
    except KeyError as error:
        raise ValueError(f"Unknown mdtbx command: {name}") from error


def arguments_to_argv(
    parser: argparse.ArgumentParser, arguments: dict[str, Any]
) -> list[str]:
    known = {
        action.dest: action
        for action in parser._actions
        if action.dest
        not in {
            "help",
            "func",
            "_agent_json",
            "_agent_dry_run",
            "_agent_cluster_profile",
        }
    }
    unknown = sorted(set(arguments) - set(known))
    if unknown:
        raise ValueError(f"Unknown arguments: {', '.join(unknown)}")

    positional = []
    optional = []
    for action in parser._actions:
        if action.dest not in arguments:
            continue
        value = arguments[action.dest]
        if value is None:
            continue
        if not action.option_strings:
            values = value if isinstance(value, list) else [value]
            positional.extend(str(item) for item in values)
            continue
        option = next(
            (item for item in action.option_strings if item.startswith("--")),
            action.option_strings[0],
        )
        if isinstance(action, argparse._StoreTrueAction):
            if value:
                optional.append(option)
            continue
        if isinstance(action, argparse._StoreFalseAction):
            if not value:
                optional.append(option)
            continue
        optional.append(option)
        values = value if isinstance(value, list) else [value]
        optional.extend(str(item) for item in values)
    argv = [*positional, *optional]
    parser.parse_args(argv)
    return argv


def normalized_arguments(
    parser: argparse.ArgumentParser, argv: list[str]
) -> dict[str, Any]:
    namespace = parser.parse_args(argv)
    return {
        key: json_value(value)
        for key, value in vars(namespace).items()
        if key
        not in {
            "func",
            "_agent_json",
            "_agent_dry_run",
            "_agent_cluster_profile",
        }
    }


def artifact_paths(arguments: dict[str, Any], item: dict[str, Any]) -> list[str]:
    paths = []
    for name in item.get("outputs", []):
        value = arguments.get(name)
        values = value if isinstance(value, list) else [value]
        for path in values:
            if isinstance(path, str) and path:
                paths.append(str(Path(path).expanduser()))
    return paths


def input_paths(arguments: dict[str, Any], item: dict[str, Any]) -> list[str]:
    return artifact_paths(arguments, {"outputs": item.get("inputs", [])})
