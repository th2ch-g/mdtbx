"""Reusable argparse argument definitions shared across subcommands.

The ``-p/--topology``, ``-t/--trajectory``, ``-o/--output``, ``-s/--selection``
and ``--gmx`` / ``-idx/--index`` flags were hand-copied (often verbatim) into
nearly every cv/analysis/trajectory subcommand. Defining them once here removes
that duplication and gives a single place to fix help text or escaping.
"""


def add_topology_arg(parser, *, help="Topology file (.gro, .pdb)"):
    parser.add_argument("-p", "--topology", type=str, required=True, help=help)


def add_trajectory_arg(parser, *, help="Trajectory file (.xtc, .trr)"):
    parser.add_argument("-t", "--trajectory", type=str, required=True, help=help)


def add_output_arg(parser, *, default="output.npy", help="Output file (.npy)"):
    parser.add_argument("-o", "--output", type=str, default=default, help=help)


def add_selection_arg(parser, *, help="Selection (MDtraj atom selection language)"):
    parser.add_argument("-s", "--selection", type=str, required=True, help=help)


def add_gmx_args(parser):
    """Add the ``--gmx`` backend switch and its companion ``-idx/--index`` flag."""
    parser.add_argument(
        "--gmx", action="store_true", help="Use the GROMACS backend instead of MDtraj"
    )
    parser.add_argument(
        "-idx", "--index", type=str, default=None, help="Index file (.ndx)"
    )
