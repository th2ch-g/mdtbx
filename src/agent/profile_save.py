"""Save an approved external cluster profile."""

import argparse
from pathlib import Path

from .runtime import save_profile


def add_subcmd(subparsers):
    parser = subparsers.add_parser(
        "agent_profile_save",
        help="Save an approved cluster profile draft",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--draft", required=True, help="Profile draft JSON")
    parser.add_argument("--output", required=True, help="External profile path")
    parser.add_argument(
        "--approve",
        required=True,
        help="Exact SHA-256 hash of the profile draft",
    )
    parser.add_argument(
        "--replace",
        action="store_true",
        help="Replace an existing profile after interactive approval",
    )
    parser.set_defaults(func=run)


def run(args):
    return save_profile(
        Path(args.draft), Path(args.output), args.approve, replace=args.replace
    )
