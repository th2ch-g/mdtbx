"""Execute or submit an approved immutable plan."""

import argparse
from pathlib import Path

from .model import read_json
from .runtime import execute_plan


def add_subcmd(subparsers):
    parser = subparsers.add_parser(
        "agent_run",
        help="Execute an approved agent plan",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--plan", required=True, help="Plan JSON")
    parser.add_argument(
        "--approve",
        required=True,
        help="Exact plan_id shown to the user for approval",
    )
    parser.add_argument(
        "--approve-unsafe",
        action="store_true",
        help="Approve an unsafe arbitrary-command step",
    )
    parser.add_argument(
        "--approve-destructive",
        action="store_true",
        help="Approve a destructive step",
    )
    parser.add_argument(
        "--approve-overwrite",
        action="store_true",
        help="Approve replacement of the exact existing artifacts in the plan",
    )
    parser.set_defaults(func=run)


def run(args):
    return execute_plan(
        read_json(Path(args.plan)),
        approval=args.approve,
        approve_unsafe=args.approve_unsafe,
        approve_destructive=args.approve_destructive,
        approve_overwrite=args.approve_overwrite,
    )
