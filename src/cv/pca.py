import argparse
from pathlib import Path
import tempfile

import numpy as np

from ..config import *  # NOQA
from ..logger import generate_logger

LOGGER = generate_logger(__name__)


def _select_atoms(topology, selection, label):
    atom_indices = topology.select(selection)
    if atom_indices.size == 0:
        raise ValueError(f"No atoms found for {label}: {selection}")
    return atom_indices


def _save_pca_metadata(output_npz, pc, pca_model, trj, atom_indices_cal, args):
    atom_names = np.array([trj.top.atom(index).name for index in atom_indices_cal])
    residue_ids = np.array(
        [trj.top.atom(index).residue.resSeq for index in atom_indices_cal]
    )
    residue_names = np.array(
        [trj.top.atom(index).residue.name for index in atom_indices_cal]
    )
    chain_ids = np.array(
        [trj.top.atom(index).residue.chain.index for index in atom_indices_cal]
    )

    np.savez(
        output_npz,
        scores=pc,
        components=pca_model.components_.reshape(pca_model.n_components_, -1, 3),
        explained_variance=pca_model.explained_variance_,
        explained_variance_ratio=pca_model.explained_variance_ratio_,
        mean_xyz=pca_model.mean_.reshape(-1, 3),
        atom_indices=atom_indices_cal,
        atom_names=atom_names,
        residue_ids=residue_ids,
        residue_names=residue_names,
        chain_ids=chain_ids,
        coordinate_unit="nm",
        unit_scale_to_angstrom=10.0,
        topology_path=args.topology,
        trajectory_path=args.trajectory,
        reference_path=args.reference,
        selection_cal_trj=args.selection_cal_trj,
        selection_cal_ref=args.selection_cal_ref,
        selection_fit_trj=args.selection_fit_trj,
        selection_fit_ref=args.selection_fit_ref,
    )
    LOGGER.info(f"Saved PCA metadata to {output_npz}")


def _extract_values_from_xvg(path, n_values):
    values = np.atleast_2d(np.loadtxt(path))
    values = values[:, -1]
    if values.shape[0] < n_values:
        raise ValueError(f"Not enough values in {path}: expected {n_values}")
    return values[:n_values]


def _save_gmx_pca_metadata(output_npz, pc, average_structure_path, args):
    import mdtraj as md

    average_trj = md.load(average_structure_path)
    eigenvec_trj = md.load("eigenvec.trr", top=average_structure_path)
    eigenvalues = _extract_values_from_xvg("eigenval.xvg", args.n_components)

    time = getattr(eigenvec_trj, "time", None)
    if time is None or len(time) == 0:
        raise ValueError("Failed to detect frames in eigenvec.trr")
    average_frame_candidates = np.where(np.isclose(time, 0.0))[0]
    if average_frame_candidates.size == 0:
        raise ValueError("Average structure frame (t=0) not found in eigenvec.trr")
    component_start = int(average_frame_candidates[0]) + 1
    component_end = component_start + args.n_components
    components = eigenvec_trj.xyz[component_start:component_end]
    if components.shape[0] < args.n_components:
        raise ValueError("Not enough eigenvector frames in eigenvec.trr")

    atom_indices = np.arange(average_trj.n_atoms)
    atom_names = np.array([atom.name for atom in average_trj.top.atoms])
    residue_ids = np.array([atom.residue.resSeq for atom in average_trj.top.atoms])
    residue_names = np.array([atom.residue.name for atom in average_trj.top.atoms])
    chain_ids = np.array([atom.residue.chain.index for atom in average_trj.top.atoms])
    eigenvalue_sum = eigenvalues.sum()
    if np.isclose(eigenvalue_sum, 0.0):
        explained_variance_ratio = np.zeros_like(eigenvalues)
    else:
        explained_variance_ratio = eigenvalues / eigenvalue_sum

    np.savez(
        output_npz,
        scores=pc,
        components=components,
        explained_variance=eigenvalues,
        explained_variance_ratio=explained_variance_ratio,
        mean_xyz=average_trj.xyz[0],
        atom_indices=atom_indices,
        atom_names=atom_names,
        residue_ids=residue_ids,
        residue_names=residue_names,
        chain_ids=chain_ids,
        coordinate_unit="nm",
        unit_scale_to_angstrom=10.0,
        topology_path=args.topology,
        trajectory_path=args.trajectory,
        reference_path=args.reference,
        selection_cal_trj=args.selection_cal_trj,
        selection_cal_ref=args.selection_cal_ref,
        selection_fit_trj=args.selection_fit_trj,
        selection_fit_ref=args.selection_fit_ref,
    )
    LOGGER.info(f"Saved PCA metadata to {output_npz}")


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
    parser.add_argument(
        "-oz",
        "--output_npz",
        type=str,
        default=None,
        help="Output PCA metadata file (.npz)",
    )
    parser.add_argument(
        "-oa",
        "--output_average",
        type=str,
        default=None,
        help="Output average structure file (.pdb, .gro)",
    )

    parser.set_defaults(func=run)


def run(args):
    if args.selection_fit_ref is None:
        args.selection_fit_ref = args.selection_fit_trj

    if args.gmx:
        # gmx
        import subprocess

        temporary_average_path = None

        if args.index is not None:
            INDEX_OPTION = f"-n {args.index}"
        else:
            INDEX_OPTION = ""

        if args.output_average is not None:
            average_structure_path = args.output_average
        elif args.output_npz is not None:
            with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as tmp_file:
                average_structure_path = tmp_file.name
            temporary_average_path = average_structure_path
        else:
            average_structure_path = None

        if average_structure_path is not None:
            AVERAGE_OPTION = f"-av {average_structure_path}"
        else:
            AVERAGE_OPTION = ""
        try:
            cmd = (
                f"gmx covar -s {args.topology} -f {args.trajectory} {INDEX_OPTION} "
                f"-xvg none -o eigenval.xvg -v eigenvec.trr {AVERAGE_OPTION}"
            )
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
            if args.output_npz is not None:
                _save_gmx_pca_metadata(
                    args.output_npz, pc, average_structure_path, args
                )
        finally:
            if temporary_average_path is not None:
                Path(temporary_average_path).unlink(missing_ok=True)
    else:
        # mdtraj
        import mdtraj as md
        from sklearn.decomposition import PCA

        trj = md.load(args.trajectory, top=args.topology)
        ref = md.load(args.reference)
        atom_indices_fit_trj = _select_atoms(
            trj.top, args.selection_fit_trj, "selection_fit_trj"
        )
        atom_indices_fit_ref = _select_atoms(
            ref.top, args.selection_fit_ref, "selection_fit_ref"
        )
        if atom_indices_fit_trj.size != atom_indices_fit_ref.size:
            raise ValueError(
                "selection_fit_trj and selection_fit_ref must have same size"
            )
        trj.superpose(
            ref,
            0,
            atom_indices=atom_indices_fit_trj,
            ref_atom_indices=atom_indices_fit_ref,
        )
        atom_indices_cal = _select_atoms(
            trj.top, args.selection_cal_trj, "selection_cal_trj"
        )
        pca_model = PCA(n_components=args.n_components)
        pc = pca_model.fit_transform(
            trj.xyz[:, atom_indices_cal, :].reshape(trj.n_frames, -1)
        )
        if args.output_npz is not None:
            _save_pca_metadata(
                args.output_npz, pc, pca_model, trj, atom_indices_cal, args
            )
        if args.output_average is not None:
            ave_xyz = trj.xyz.mean(axis=0, keepdims=True)
            md.Trajectory(ave_xyz, trj.topology).save(args.output_average)
            LOGGER.info(f"Saved average structure to {args.output_average}")
    np.save(args.output, pc)
    LOGGER.info(f"Saved to {args.output}")
