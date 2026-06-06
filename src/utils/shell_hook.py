import argparse
from pathlib import Path

from ..logger import generate_logger

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    """
    mdtbx shell_hook
    """
    parser = subparsers.add_parser(
        "shell_hook",
        help="Generate shell hook",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.set_defaults(func=run)


def run(args):
    project_root = Path(__file__).parent.parent.parent
    hook = f"""
----- BEGIN OF MDTBX SHELL HOOK -----
# for alias
shopt -s expand_aliases
alias mdtbx="pixi run --manifest-path {project_root} mdtbx"
alias pymol="pixi run --manifest-path {project_root} pymol"
alias python3="pixi run --manifest-path {project_root} python3"

# for export
export PATH={project_root}/.pixi/envs/default/bin:$PATH

# for link
ln -s {project_root}/.pixi/envs/default/bin/mdtbx .
ln -s {project_root}/.pixi/envs/default/bin/pymol .

# for function
mdtbx() {{
    pixi run --manifest-path {project_root} mdtbx "$@"
}}

pymol() {{
    pixi run --manifest-path {project_root} pymol "$@"
}}

python3() {{
    pixi run --manifest-path {project_root} python3 "$@"
}}

# pymol template
pymol -c -p <<EOF
import pymol_plugins
import pymol
from pymol_plugins import *
from pymol import *

EOF

# python3 template
python3 -u <<EOF

EOF
----- END OF MDTBX SHELL HOOK -----
    """
    print(hook)
