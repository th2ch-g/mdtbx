import argparse
import sys

from .logger import generate_logger

from .utils import rmfile
from .utils import addace
from .utils import addnme
from .utils import add_ndx
from .utils import convert
from .utils import calc_ion_conc
from .utils import centering_gro
from .utils import find_bond
from .utils import gen_am1bcc
from .utils import gen_resp
from .utils import gen_posres
from .utils import gen_distres
from .utils import amb2gro
from .utils import trjcat
from .utils import print_perf
from .utils import mv_crds_mol2

# from .utils import build_membrane
from .utils import build_solution

from .simulator import gen_sample as gen_sample_simulator
from .builder import gen_sample_builder
from .builder import gen_sample_mdp
from .msm import gen_sample as gen_sample_msm

from .cv import comdist
from .cv import comvec
from .cv import rmsd
from .cv import xyz

LOGGER = generate_logger(__name__)


def cli() -> None:
    # make parser
    parser = argparse.ArgumentParser(description=("ToolBox for MD simulation"))

    # subcommands
    subparsers = parser.add_subparsers()

    addace.add_subcmd(subparsers)
    addnme.add_subcmd(subparsers)
    add_ndx.add_subcmd(subparsers)
    mv_crds_mol2.add_subcmd(subparsers)
    gen_am1bcc.add_subcmd(subparsers)
    gen_resp.add_subcmd(subparsers)
    gen_posres.add_subcmd(subparsers)
    gen_distres.add_subcmd(subparsers)
    find_bond.add_subcmd(subparsers)
    convert.add_subcmd(subparsers)
    calc_ion_conc.add_subcmd(subparsers)
    centering_gro.add_subcmd(subparsers)
    amb2gro.add_subcmd(subparsers)
    trjcat.add_subcmd(subparsers)
    rmfile.add_subcmd(subparsers)
    print_perf.add_subcmd(subparsers)

    # build_membrane.add_subcmd(subparsers)
    build_solution.add_subcmd(subparsers)

    gen_sample_builder.add_subcmd(subparsers)
    gen_sample_mdp.add_subcmd(subparsers)
    gen_sample_simulator.add_subcmd(subparsers)
    gen_sample_msm.add_subcmd(subparsers)

    comdist.add_subcmd(subparsers)
    comvec.add_subcmd(subparsers)
    rmsd.add_subcmd(subparsers)
    xyz.add_subcmd(subparsers)

    args = parser.parse_args()

    if len(sys.argv) == 1:
        LOGGER.error(f"use {sys.argv[0]} --help")
        sys.exit(1)
    LOGGER.info(f"{sys.argv[1]} called")

    if sys.argv[1] == "rmfile":
        rmfile.run(args)

    if sys.argv[1] == "addace":
        addace.run(args)

    if sys.argv[1] == "addnme":
        addnme.run(args)

    if sys.argv[1] == "mv_crds_mol2":
        mv_crds_mol2.run(args)

    if sys.argv[1] == "trjcat":
        trjcat.run(args)

    if sys.argv[1] == "convert":
        convert.run(args)

    if sys.argv[1] == "centering_gro":
        centering_gro.run(args)

    if sys.argv[1] == "amb2gro":
        amb2gro.run(args)

    if sys.argv[1] == "find_bond":
        find_bond.run(args)

    if sys.argv[1] == "gen_am1bcc":
        gen_am1bcc.run(args)

    if sys.argv[1] == "gen_resp":
        gen_resp.run(args)

    if sys.argv[1] == "gen_posres":
        gen_posres.run(args)

    if sys.argv[1] == "gen_distres":
        gen_distres.run(args)

    if sys.argv[1] == "add_ndx":
        add_ndx.run(args)

    if sys.argv[1] == "print_perf":
        print_perf.run(args)

    if sys.argv[1] == "calc_ion_conc":
        calc_ion_conc.run(args)

    # if sys.argv[1] == "build_membrane":
    #     build_membrane.run(args)

    if sys.argv[1] == "build_solution":
        build_solution.run(args)

    if sys.argv[1] == "gen_sample_simulator":
        gen_sample_simulator.run(args)

    if sys.argv[1] == "gen_sample_builder":
        gen_sample_builder.run(args)

    if sys.argv[1] == "gen_sample_mdp":
        gen_sample_mdp.run(args)

    if sys.argv[1] == "gen_sample_msm":
        gen_sample_msm.run(args)

    if sys.argv[1] == "comdist":
        comdist.run(args)

    if sys.argv[1] == "comvec":
        comvec.run(args)

    if sys.argv[1] == "rmsd":
        rmsd.run(args)

    if sys.argv[1] == "xyz":
        xyz.run(args)

    LOGGER.info(f"{sys.argv[1]} finished")
