import argparse
import sys
from importlib import metadata

from .logger import generate_logger

from .utils import rmfile
from .utils import addace
from .utils import addnme
from .utils import add_ndx
from .utils import convert
from .utils import calc_ion_conc
from .utils import centering_gro
from .utils import find_bond
from .utils import mod_mdp
from .utils import gen_am1bcc
from .utils import gen_resp
from .utils import gen_modres_am1bcc
from .utils import gen_modres_resp
from .utils import gen_posres
from .utils import gen_distres
from .utils import amb2gro
from .utils import trjcat
from .utils import print_perf
from .utils import mv_crds_mol2
from .utils import shell_hook
from .utils import cmd

# from .utils import build_membrane
from .utils import build_solution

from .cv import comdist
from .cv import comvec
from .cv import mindist
from .cv import rmsd
from .cv import xyz

LOGGER = generate_logger(__name__)


def get_version():
    try:
        return metadata.version("mdtbx")
    except metadata.PackageNotFoundError:
        return "unknown"


def cli() -> None:
    # make parser
    parser = argparse.ArgumentParser(description=("ToolBox for MD simulation"))

    parser.add_argument(
        "-V",
        "--version",
        action="version",
        version=f"%(prog)s {get_version()}",
        help="Print version",
    )

    # subcommands
    subparsers = parser.add_subparsers()

    addace.add_subcmd(subparsers)
    addnme.add_subcmd(subparsers)
    add_ndx.add_subcmd(subparsers)
    mv_crds_mol2.add_subcmd(subparsers)
    gen_am1bcc.add_subcmd(subparsers)
    gen_resp.add_subcmd(subparsers)
    gen_modres_am1bcc.add_subcmd(subparsers)
    gen_modres_resp.add_subcmd(subparsers)
    gen_posres.add_subcmd(subparsers)
    gen_distres.add_subcmd(subparsers)
    find_bond.add_subcmd(subparsers)
    mod_mdp.add_subcmd(subparsers)
    convert.add_subcmd(subparsers)
    calc_ion_conc.add_subcmd(subparsers)
    centering_gro.add_subcmd(subparsers)
    amb2gro.add_subcmd(subparsers)
    trjcat.add_subcmd(subparsers)
    rmfile.add_subcmd(subparsers)
    print_perf.add_subcmd(subparsers)
    shell_hook.add_subcmd(subparsers)
    cmd.add_subcmd(subparsers)

    # build_membrane.add_subcmd(subparsers)
    build_solution.add_subcmd(subparsers)

    comdist.add_subcmd(subparsers)
    comvec.add_subcmd(subparsers)
    mindist.add_subcmd(subparsers)
    rmsd.add_subcmd(subparsers)
    xyz.add_subcmd(subparsers)

    args = parser.parse_args()

    if len(sys.argv) == 1:
        LOGGER.error(f"use {sys.argv[0]} --help")
        sys.exit(1)
    LOGGER.info(f"{sys.argv[1]} called")

    if sys.argv[1] == "rmfile":
        rmfile.run(args)

    elif sys.argv[1] == "addace":
        addace.run(args)

    elif sys.argv[1] == "addnme":
        addnme.run(args)

    elif sys.argv[1] == "mv_crds_mol2":
        mv_crds_mol2.run(args)

    elif sys.argv[1] == "trjcat":
        trjcat.run(args)

    elif sys.argv[1] == "convert":
        convert.run(args)

    elif sys.argv[1] == "centering_gro":
        centering_gro.run(args)

    elif sys.argv[1] == "amb2gro":
        amb2gro.run(args)

    elif sys.argv[1] == "find_bond":
        find_bond.run(args)

    elif sys.argv[1] == "mod_mdp":
        mod_mdp.run(args)

    elif sys.argv[1] == "gen_am1bcc":
        gen_am1bcc.run(args)

    elif sys.argv[1] == "gen_resp":
        gen_resp.run(args)

    elif sys.argv[1] == "gen_modres_am1bcc":
        gen_modres_am1bcc.run(args)

    elif sys.argv[1] == "gen_modres_resp":
        gen_modres_resp.run(args)

    elif sys.argv[1] == "gen_posres":
        gen_posres.run(args)

    elif sys.argv[1] == "gen_distres":
        gen_distres.run(args)

    elif sys.argv[1] == "add_ndx":
        add_ndx.run(args)

    elif sys.argv[1] == "print_perf":
        print_perf.run(args)

    elif sys.argv[1] == "shell_hook":
        shell_hook.run(args)

    elif sys.argv[1] == "cmd":
        cmd.run(args)

    elif sys.argv[1] == "calc_ion_conc":
        calc_ion_conc.run(args)

    # elif sys.argv[1] == "build_membrane":
    #     build_membrane.run(args)

    elif sys.argv[1] == "build_solution":
        build_solution.run(args)

    elif sys.argv[1] == "comdist":
        comdist.run(args)

    elif sys.argv[1] == "comvec":
        comvec.run(args)

    elif sys.argv[1] == "mindist":
        mindist.run(args)

    elif sys.argv[1] == "rmsd":
        rmsd.run(args)

    elif sys.argv[1] == "xyz":
        xyz.run(args)

    else:
        print(f"Unknown command: {sys.argv[1]}")

    LOGGER.info(f"{sys.argv[1]} finished")
