"""Cancel queued or running jobs belonging to an immutable agent run."""

from .runtime import cancel_run


def add_subcmd(subparsers):
    parser = subparsers.add_parser(
        "agent_cancel",
        help="Cancel queued or running jobs from an agent run",
    )
    parser.add_argument("--run", required=True, help="Run directory or run ID")
    parser.add_argument(
        "--approve",
        required=True,
        help="Exact plan_id previously approved for this run",
    )
    parser.set_defaults(func=run)


def run(args):
    return cancel_run(args.run, args.approve)
