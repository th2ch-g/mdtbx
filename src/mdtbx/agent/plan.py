"""Create an immutable execution plan from a typed request."""

import argparse
from pathlib import Path

from .model import read_json, write_json
from .runtime import build_plan


def add_subcmd(subparsers):
    parser = subparsers.add_parser(
        "agent_plan",
        help="Create a resource-aware immutable execution plan",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--request", required=True, help="Agent request JSON")
    parser.add_argument("-o", "--output", help="Write the plan to this path")
    parser.set_defaults(func=run)


def run(args):
    plan = build_plan(read_json(Path(args.request)))
    if args.output:
        write_json(Path(args.output), plan)
    return plan
