import argparse

from pymol import cmd as pymol_cmd

from ..logger import generate_logger

LOGGER = generate_logger(__name__)

# PDB convention: single-letter chain IDs, uppercase first then lowercase.
CHAIN_POOL = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"


def add_subcmd(subparsers):
    parser = subparsers.add_parser(
        "fill_chainname",
        help="Assign chain IDs to atoms whose chain field is blank",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-s",
        "--structure",
        required=True,
        type=str,
        help="Input structure file (PDB etc.)",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="filled.pdb",
        type=str,
        help="Output PDB file",
    )
    parser.add_argument(
        "--max_resi_gap",
        default=1,
        type=int,
        help="Residue-number gap larger than this starts a new chain block",
    )
    parser.add_argument(
        "--ignore_segi",
        action="store_true",
        help="Do not split blank-chain atoms when segi changes",
    )
    parser.add_argument(
        "--chainname",
        nargs="+",
        default=None,
        type=str,
        help="Chain IDs to assign (one per blank-chain block). "
        "If omitted, IDs are auto-picked from A-Z then a-z",
    )

    parser.set_defaults(func=run)


def _collect_blank_atoms(object_name: str) -> list[tuple[int, int, str]]:
    """Return [(atom_index, resv, segi)] for atoms whose chain is blank, in atom order."""
    atoms: list[tuple[int, int, str]] = []
    pymol_cmd.iterate(
        f"({object_name}) and chain ''",
        "atoms.append((index, resv, segi))",
        space={"atoms": atoms},
    )
    atoms.sort(key=lambda a: a[0])
    return atoms


def _split_into_blocks(
    atoms: list[tuple[int, int, str]],
    max_resi_gap: int,
    use_segi: bool,
) -> list[list[int]]:
    """Group blank-chain atoms into chain blocks by residue-gap and optional segi change."""
    blocks: list[list[int]] = []
    prev_resv: int | None = None
    prev_segi: str | None = None
    for idx, resv, segi in atoms:
        boundary = prev_resv is None
        if not boundary:
            if use_segi and segi != prev_segi:
                boundary = True
            elif resv - prev_resv > max_resi_gap:
                boundary = True
        if boundary:
            blocks.append([])
        blocks[-1].append(idx)
        prev_resv, prev_segi = resv, segi
    return blocks


def _next_chain_id(used: set[str]) -> str:
    for c in CHAIN_POOL:
        if c not in used:
            return c
    raise RuntimeError("No single-letter chain ID left to assign")


def _validate_user_chainnames(
    chainnames: list[str], used: set[str], n_blocks: int
) -> None:
    if len(chainnames) != n_blocks:
        raise ValueError(
            f"--chainname has {len(chainnames)} value(s) "
            f"but {n_blocks} blank-chain block(s) need to be filled"
        )
    if len(set(chainnames)) != len(chainnames):
        raise ValueError(f"--chainname contains duplicates: {chainnames}")
    for c in chainnames:
        if len(c) != 1 or c not in CHAIN_POOL:
            raise ValueError(f"Invalid chain ID '{c}'; must be a single A-Za-z letter")
        if c in used:
            raise ValueError(f"Chain ID '{c}' already exists in the structure")


def run(args):
    object_name = "target"

    pymol_cmd.reinitialize()
    pymol_cmd.load(args.structure, object_name)

    used = {c for c in pymol_cmd.get_chains(object_name) if c}
    LOGGER.info(f"Existing chain IDs: {sorted(used) if used else '(none)'}")

    atoms = _collect_blank_atoms(object_name)
    if not atoms:
        LOGGER.info("No blank-chain atoms found; saving structure as-is")
        pymol_cmd.save(args.output, object_name)
        LOGGER.info(f"{args.output} generated")
        return

    blocks = _split_into_blocks(
        atoms,
        max_resi_gap=args.max_resi_gap,
        use_segi=not args.ignore_segi,
    )
    LOGGER.info(f"Found {len(blocks)} blank-chain block(s) to fill")

    if args.chainname is not None:
        _validate_user_chainnames(args.chainname, used, len(blocks))
        assigned_ids = list(args.chainname)
    else:
        assigned_ids = []
        for _ in blocks:
            new_id = _next_chain_id(used)
            used.add(new_id)
            assigned_ids.append(new_id)

    for block, new_id in zip(blocks, assigned_ids):
        sel = "index " + "+".join(str(i) for i in block)
        pymol_cmd.alter(f"({object_name}) and ({sel})", f"chain='{new_id}'")
        LOGGER.info(f"Assigned chain '{new_id}' to {len(block)} atom(s)")

    pymol_cmd.sort()
    pymol_cmd.save(args.output, object_name)
    LOGGER.info(f"{args.output} generated")
