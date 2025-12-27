import argparse
import numpy as np

from ..config import *  # NOQA
from ..logger import generate_logger

LOGGER = generate_logger(__name__)


def add_subcmd(subparsers):
    """
    mdtbx pca --topology structure.pdb --trajectory trajectory.xtc --reference reference.pdb --selection_cal_trj "resid 1 to 10" --selection_cal_ref "resid 1 to 10" --selection_fit_trj "resid 11 to 20" --selection_fit_ref "resid 11 to 20" -o cv.npy -n 4
    """
    parser = subparsers.add_parser(
        "pca",
        help="Extract PCA",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-p", "--topology", type=str, required=True, help="Topology file (.gro, .pdb)"
    )
    parser.add_argument(
        "-t",
        "--trajectory",
        type=str,
        required=True,
        help="Trajectory file (.xtc, .trr)",
    )
    parser.add_argument(
        "-r",
        "--reference",
        type=str,
        required=True,
        help="Reference structure file (.gro, .pdb)",
    )
    parser.add_argument(
        "-sct",
        "--selection_cal_trj",
        type=str,
        required=True,
        help="Selection for trajectory (MDtraj Atom selection language)",
    )
    parser.add_argument(
        "-scr",
        "--selection_cal_ref",
        type=str,
        required=True,
        help="Selection for reference (MDtraj Atom selection language)",
    )
    parser.add_argument(
        "-sft",
        "--selection_fit_trj",
        type=str,
        required=True,
        help="Selection for trajectory (MDtraj Atom selection language)",
    )
    parser.add_argument(
        "-sfr",
        "--selection_fit_ref",
        type=str,
        help="Selection for reference (MDtraj Atom selection language)",
    )
    parser.add_argument(
        "--gmx",
        action="store_true",
        help="Use gmx",
    )
    parser.add_argument(
        "-n",
        "--n_components",
        type=int,
        default=10,
        help="Number of components",
    )
    parser.add_argument(
        "-idx",
        "--index",
        type=str,
        default=None,
        help="Index file (.ndx)",
    )
    parser.add_argument(
        "-o", "--output", type=str, default="pca.npy", help="Output file (.npy)"
    )


def run(args):
    if args.gmx:
        # gmx
        import subprocess

        if args.index is not None:
            INDEX_OPTION = f"-n {args.index}"
        else:
            INDEX_OPTION = ""
        cmd = f"gmx covar -s {args.topology} -f {args.trajectory} {INDEX_OPTION} -xvg none -o eigenval.xvg -v eigenvec.trr"
        LOGGER.info("gmx covar started")
        subprocess.run(
            cmd,
            input=f"{args.selection_fit_trj}\n{args.selection_cal_trj}\n",
            shell=True,
            check=True,
        )
        LOGGER.info("gmx covar finished")

        cmd = f"gmx anaeig -s {args.topology} -f {args.trajectory} {INDEX_OPTION} -v eigenvec.trr -proj proj.xvg -xvg none -eig eigenval.xvg -first 1 -last {args.n_components}"
        LOGGER.info("gmx anaeig started")
        subprocess.run(
            cmd,
            input=f"{args.selection_fit_trj}\n{args.selection_cal_trj}\n",
            shell=True,
            check=True,
        )
        LOGGER.info("gmx anaeig finished")

        pc = np.loadtxt("proj.xvg")
    else:
        # mdtraj
        import mdtraj as md
        from sklearn.decomposition import PCA

        trj = md.load(args.trajectory, top=args.topology)
        ref = md.load(args.reference)
        atom_indices_fit_trj = trj.top.select(args.selection_fit_trj)
        atom_indices_fit_ref = trj.top.select(args.selection_fit_ref)
        trj.superpose(
            ref,
            0,
            atom_indices=atom_indices_fit_ref,
            ref_atom_indices=atom_indices_fit_trj,
        )
        atom_indices_cal = trj.top.select(args.selection_cal_trj)
        pca_model = PCA(n_components=args.n_components)
        pc = pca_model.fit_transform(
            trj.xyz[:, atom_indices_cal, :].reshape(trj.n_frames, -1)
        )
    np.save(args.output, pc)
    LOGGER.info(f"Saved to {args.output}")
