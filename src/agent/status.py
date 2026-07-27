"""Read normalized status for an agent run."""

import argparse

from .runtime import run_status


def add_subcmd(subparsers):
    parser = subparsers.add_parser(
        "agent_status",
        help="Read normalized scheduler status for an agent run",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--run", required=True, help="Run directory or run ID")
    parser.set_defaults(func=run)


def run(args):
    return run_status(args.run)
