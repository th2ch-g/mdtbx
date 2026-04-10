import argparse
import sys
from importlib import metadata

from .logger import generate_logger

# general util
from .utils import rmfile
from .utils import convert
from .utils import mod_mdp
from .utils import shell_hook
from .utils import cmd
from .utils import show_mdtraj
from .utils import show_npy
from .utils import partial_tempering

# build utils
from .build import addace
from .build import addh
from .build import addnme
from .build import add_ndx
from .build import mv_crds_mol2
from .build import calc_ion_conc
from .build import centering_gro
from .build import find_bond
from .build import gen_am1bcc
from .build import gen_resp
from .build import gen_modres_am1bcc
from .build import gen_modres_resp
from .build import gen_posres
from .build import gen_distres
from .build import modeling_cf
from .build import amb2gro
from .build import build_solution
from .build import build_vacuum
from .build import place_solvent
from .build import gen_temperatures

# trajectory utils
from .trajectory import trjcat
from .trajectory import fit
from .trajectory import pacs_trjcat
from .trajectory import print_perf

# analysis utils
from .analysis import extract_ave_str
from .analysis import extract_str

from .cv import comdist
from .cv import comvec
from .cv import mindist
from .cv import rmsd
from .cv import rmsf
from .cv import xyz
from .cv import pca
from .cv import densmap

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
    addh.add_subcmd(subparsers)
    addnme.add_subcmd(subparsers)
    add_ndx.add_subcmd(subparsers)
    mv_crds_mol2.add_subcmd(subparsers)
    gen_am1bcc.add_subcmd(subparsers)
    gen_resp.add_subcmd(subparsers)
    gen_modres_am1bcc.add_subcmd(subparsers)
    gen_modres_resp.add_subcmd(subparsers)
    gen_posres.add_subcmd(subparsers)
    gen_distres.add_subcmd(subparsers)
    modeling_cf.add_subcmd(subparsers)
    find_bond.add_subcmd(subparsers)
    mod_mdp.add_subcmd(subparsers)
    convert.add_subcmd(subparsers)
    calc_ion_conc.add_subcmd(subparsers)
    centering_gro.add_subcmd(subparsers)
    amb2gro.add_subcmd(subparsers)
    trjcat.add_subcmd(subparsers)
    fit.add_subcmd(subparsers)
    pacs_trjcat.add_subcmd(subparsers)
    rmfile.add_subcmd(subparsers)
    extract_ave_str.add_subcmd(subparsers)
    extract_str.add_subcmd(subparsers)
    show_mdtraj.add_subcmd(subparsers)
    show_npy.add_subcmd(subparsers)
    print_perf.add_subcmd(subparsers)
    shell_hook.add_subcmd(subparsers)
    partial_tempering.add_subcmd(subparsers)
    gen_temperatures.add_subcmd(subparsers)
    cmd.add_subcmd(subparsers)

    # build_membrane.add_subcmd(subparsers)
    build_solution.add_subcmd(subparsers)
    build_vacuum.add_subcmd(subparsers)
    place_solvent.add_subcmd(subparsers)

    comdist.add_subcmd(subparsers)
    comvec.add_subcmd(subparsers)
    mindist.add_subcmd(subparsers)
    rmsd.add_subcmd(subparsers)
    rmsf.add_subcmd(subparsers)
    xyz.add_subcmd(subparsers)
    pca.add_subcmd(subparsers)
    densmap.add_subcmd(subparsers)

    args = parser.parse_args()

    if not hasattr(args, "func"):
        LOGGER.error(f"use {sys.argv[0]} --help")
        sys.exit(1)

    LOGGER.info(f"{sys.argv[1]} called")
    args.func(args)
    LOGGER.info(f"{sys.argv[1]} finished")
