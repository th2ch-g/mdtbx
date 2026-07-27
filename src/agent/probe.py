"""Probe a batch scheduler without changing cluster state."""

import argparse

from .runtime import probe_cluster


def add_subcmd(subparsers):
    parser = subparsers.add_parser(
        "agent_probe",
        help="Inspect scheduler resources and produce a profile draft",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--scheduler",
        choices=["slurm", "age", "pjm"],
        help="Scheduler backend; auto-detect when omitted",
    )
    parser.add_argument("--cluster-profile", help="Existing external JSON profile")
    parser.set_defaults(func=run)


def run(args):
    return probe_cluster(args.scheduler, args.cluster_profile)
