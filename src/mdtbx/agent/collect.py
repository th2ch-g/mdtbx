"""Collect normalized results for an agent run."""

import argparse

from .runtime import collect_run


def add_subcmd(subparsers):
    parser = subparsers.add_parser(
        "agent_collect",
        help="Collect scheduler status and command results",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--run", required=True, help="Run directory or run ID")
    parser.set_defaults(func=run)


def run(args):
    return collect_run(args.run)
