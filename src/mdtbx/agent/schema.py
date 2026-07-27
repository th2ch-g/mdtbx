"""Expose the argparse command surface as JSON."""

import argparse

from .runtime import schema


def add_subcmd(subparsers):
    parser = subparsers.add_parser(
        "agent_schema",
        help="Describe mdtbx commands for autonomous agents",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("command", nargs="?", help="Optional command name")
    parser.set_defaults(func=run)


def run(args):
    return schema(args.command)
