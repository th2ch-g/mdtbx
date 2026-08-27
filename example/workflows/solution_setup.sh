#!/usr/bin/env bash
set -euo pipefail

# Copy this file into a calculation directory and edit this block.
MDTBX=(mdtbx)
# From a source checkout, use:
# MDTBX=(pixi run --manifest-path /path/to/mdtbx/pyproject.toml mdtbx)
INPUT_PDB="${INPUT_PDB:-prepared.pdb}"
OUTPUT_DIR="${OUTPUT_DIR:-run}"
MDP_SOURCE_DIR="${MDP_SOURCE_DIR:-/path/to/mdtbx/example/mdp/solution}"
WATER_MODEL="${WATER_MODEL:-tip3p}"
ION_CONC="${ION_CONC:-0.15}"
CATION="${CATION:-Na+}"
ANION="${ANION:-Cl-}"
# Angstrom.
BOX_SIZE=(${BOX_SIZE:-100 100 100})
WATER_SEED="${WATER_SEED:-2026}"
# GROMACS index group name.
CENTER_SELECTION="${CENTER_SELECTION:-Protein}"
# mdtbx atom-selection expression.
POSRES_SELECTION="${POSRES_SELECTION:-protein and backbone}"
LIGAND_FRCMOD="${LIGAND_FRCMOD:-}"
LIGAND_LIB="${LIGAND_LIB:-}"
TLEAP_PRE_COMMAND="${TLEAP_PRE_COMMAND:-}"
TLEAP_POST_COMMAND="${TLEAP_POST_COMMAND:-}"

CHECK_ONLY=false
case "${1:-}" in
    "") ;;
    --check) CHECK_ONLY=true ;;
    *) echo "Usage: $0 [--check]" >&2; exit 2 ;;
esac

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
cd "$SCRIPT_DIR"

required_mdps=(mini eq1 eq2 eq3 eq4 eq5 eq6 prd)

fail() {
    echo "ERROR: $*" >&2
    exit 1
}

print_command() {
    printf '  '
    printf '%q ' "$@"
    printf '\n'
}

require_file() {
    [[ -f "$1" ]] || fail "File not found: $1"
}

require_file "$INPUT_PDB"
command -v "${MDTBX[0]}" >/dev/null || fail "Command not found: ${MDTBX[0]}"
[[ -d "$MDP_SOURCE_DIR" ]] || fail "MDP directory not found: $MDP_SOURCE_DIR"
[[ ${#BOX_SIZE[@]} -eq 3 ]] || fail "BOX_SIZE must contain three values"

for stage in "${required_mdps[@]}"; do
    require_file "${MDP_SOURCE_DIR}/${stage}.mdp"
done

if [[ -n "$LIGAND_FRCMOD" || -n "$LIGAND_LIB" ]]; then
    [[ -n "$LIGAND_FRCMOD" && -n "$LIGAND_LIB" ]] || \
        fail "Set both LIGAND_FRCMOD and LIGAND_LIB"
    require_file "$LIGAND_FRCMOD"
    require_file "$LIGAND_LIB"
fi

[[ ! -e "$OUTPUT_DIR" ]] || fail "Refusing to overwrite: $OUTPUT_DIR"

build_command=(
    "${MDTBX[@]}" build_solution
    -i "$INPUT_PDB"
    -o "${OUTPUT_DIR}/amber"
    --water "$WATER_MODEL"
    --ion_conc "$ION_CONC"
    --cation "$CATION"
    --anion "$ANION"
    --boxsize "${BOX_SIZE[@]}"
    --water-seed "$WATER_SEED"
)
if [[ -n "$LIGAND_FRCMOD" ]]; then
    build_command+=(--ligparam "${LIGAND_FRCMOD}:${LIGAND_LIB}")
fi
if [[ -n "$TLEAP_PRE_COMMAND" ]]; then
    build_command+=(--addprecmd "$TLEAP_PRE_COMMAND")
fi
if [[ -n "$TLEAP_POST_COMMAND" ]]; then
    build_command+=(--addpostcmd "$TLEAP_POST_COMMAND")
fi

commands=(
    "${build_command[*]}"
    "${MDTBX[*]} amb2gro -p ${OUTPUT_DIR}/amber/leap.parm7 -x ${OUTPUT_DIR}/amber/leap.rst7 --type parmed -o $OUTPUT_DIR"
    "${MDTBX[*]} add_ndx -g ${OUTPUT_DIR}/gmx.gro -o ${OUTPUT_DIR}/index.ndx"
    "${MDTBX[*]} centering_gro -f ${OUTPUT_DIR}/gmx.gro -p ${OUTPUT_DIR}/gmx.top -idx ${OUTPUT_DIR}/index.ndx -c $CENTER_SELECTION -o ${OUTPUT_DIR}/gmx_centered.gro"
    "${MDTBX[*]} gen_posres -p ${OUTPUT_DIR}/gmx.top -s $POSRES_SELECTION -o ${OUTPUT_DIR}/posres"
)

if $CHECK_ONLY; then
    echo "Configuration is valid. Planned commands:"
    print_command "${build_command[@]}"
    for command in "${commands[@]:1}"; do
        echo "  $command"
    done
    echo "  Copy MDP files from $MDP_SOURCE_DIR"
    exit 0
fi

"${build_command[@]}"
"${MDTBX[@]}" amb2gro \
    -p "${OUTPUT_DIR}/amber/leap.parm7" \
    -x "${OUTPUT_DIR}/amber/leap.rst7" \
    --type parmed \
    -o "$OUTPUT_DIR"
"${MDTBX[@]}" add_ndx \
    -g "${OUTPUT_DIR}/gmx.gro" \
    -o "${OUTPUT_DIR}/index.ndx"
"${MDTBX[@]}" centering_gro \
    -f "${OUTPUT_DIR}/gmx.gro" \
    -p "${OUTPUT_DIR}/gmx.top" \
    -idx "${OUTPUT_DIR}/index.ndx" \
    -c "$CENTER_SELECTION" \
    -o "${OUTPUT_DIR}/gmx_centered.gro"
mv "${OUTPUT_DIR}/gmx_centered.gro" "${OUTPUT_DIR}/gmx.gro"
"${MDTBX[@]}" gen_posres \
    -p "${OUTPUT_DIR}/gmx.top" \
    -s "$POSRES_SELECTION" \
    -o "${OUTPUT_DIR}/posres"
shopt -s nullglob
posres_files=("${OUTPUT_DIR}"/posres_*.itp)
shopt -u nullglob
[[ ${#posres_files[@]} -gt 0 ]] || \
    fail "POSRES_SELECTION generated no restraint include files"

for stage in "${required_mdps[@]}"; do
    cp "${MDP_SOURCE_DIR}/${stage}.mdp" "${OUTPUT_DIR}/${stage}.mdp"
done

{
    printf 'workflow=solution\n'
    printf 'input_pdb=%s\n' "$INPUT_PDB"
    printf 'mdp_source=%s\n' "$MDP_SOURCE_DIR"
    printf 'mdtbx_command=%s\n' "${MDTBX[*]}"
    printf 'water_model=%s\n' "$WATER_MODEL"
    printf 'ion_concentration_molar=%s\n' "$ION_CONC"
    printf 'ions=%s,%s\n' "$CATION" "$ANION"
    printf 'box_size_angstrom=%s\n' "${BOX_SIZE[*]}"
    printf 'water_seed=%s\n' "$WATER_SEED"
    printf 'center_group=%s\n' "$CENTER_SELECTION"
    printf 'posres_selection=%s\n' "$POSRES_SELECTION"
    printf 'ligand_frcmod=%s\n' "$LIGAND_FRCMOD"
    printf 'ligand_lib=%s\n' "$LIGAND_LIB"
} > "${OUTPUT_DIR}/workflow.txt"

echo "Prepared $OUTPUT_DIR"
echo "Next: copy run_slurm.sh into $OUTPUT_DIR, edit it, and run --check."
