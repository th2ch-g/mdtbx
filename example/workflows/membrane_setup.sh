#!/usr/bin/env bash
set -euo pipefail

# Copy this file into a calculation directory and edit this block.
MDTBX=(mdtbx)
# From a source checkout, use:
# MDTBX=(pixi run --manifest-path /path/to/mdtbx/pyproject.toml mdtbx)
INPUT_PDB="${INPUT_PDB:-prepared_oriented.pdb}"
OUTPUT_DIR="${OUTPUT_DIR:-run}"
MDP_SOURCE_DIR="${MDP_SOURCE_DIR:-/path/to/mdtbx/example/mdp/membrane}"
LIPIDS="${LIPIDS:-POPC:CHL1}"
LIPID_RATIO="${LIPID_RATIO:-4:1}"
SALT_CONC="${SALT_CONC:-0.15}"
CATION="${CATION:-Na+}"
ANION="${ANION:-Cl-}"
# Angstrom.
read -r -a BOX_DIMS <<< "${BOX_DIMS:-120 120 200}"
WATER_MODEL="${WATER_MODEL:-tip3p}"
PROTEIN_FF="${PROTEIN_FF:-ff14SB}"
LIPID_FF="${LIPID_FF:-lipid21}"
PREORIENTED="${PREORIENTED:-true}"
# GROMACS index group name.
CENTER_SELECTION="${CENTER_SELECTION:-Protein}"
# mdtbx atom-selection expression.
POSRES_SELECTION="${POSRES_SELECTION:-protein and backbone}"
# Default: all atoms in the Lipid21 POPC/CHL1 residue fragments.
LIPID_POSRES_SELECTION="${LIPID_POSRES_SELECTION:-resname PA or resname PC or resname OL or resname CHL}"
LIGAND_FRCMOD="${LIGAND_FRCMOD:-}"
LIGAND_LIB="${LIGAND_LIB:-}"
TLEAP_LINE="${TLEAP_LINE:-}"

CHECK_ONLY=false
case "${1:-}" in
    "") ;;
    --check) CHECK_ONLY=true ;;
    *)
        echo "Usage: $0 [--check]" >&2
        exit 2
        ;;
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

absolute_file() {
    local directory
    local filename
    directory="$(cd "$(dirname "$1")" && pwd -P)"
    filename="$(basename "$1")"
    printf '%s/%s\n' "$directory" "$filename"
}

require_file "$INPUT_PDB"
command -v "${MDTBX[0]}" > /dev/null || fail "Command not found: ${MDTBX[0]}"
[[ -d "$MDP_SOURCE_DIR" ]] || fail "MDP directory not found: $MDP_SOURCE_DIR"
[[ ${#BOX_DIMS[@]} -eq 3 ]] || fail "BOX_DIMS must contain three values"
[[ "$PREORIENTED" == true || "$PREORIENTED" == false ]] ||
    fail "PREORIENTED must be true or false"

for stage in "${required_mdps[@]}"; do
    require_file "${MDP_SOURCE_DIR}/${stage}.mdp"
done

if [[ -n "$LIGAND_FRCMOD" || -n "$LIGAND_LIB" ]]; then
    [[ -n "$LIGAND_FRCMOD" && -n "$LIGAND_LIB" ]] ||
        fail "Set both LIGAND_FRCMOD and LIGAND_LIB"
    require_file "$LIGAND_FRCMOD"
    require_file "$LIGAND_LIB"
    LIGAND_FRCMOD="$(absolute_file "$LIGAND_FRCMOD")"
    LIGAND_LIB="$(absolute_file "$LIGAND_LIB")"
fi

[[ ! -e "$OUTPUT_DIR" ]] || fail "Refusing to overwrite: $OUTPUT_DIR"

memgen_command=(
    "${MDTBX[@]}" cmd packmol-memgen
    --pdb system.pdb
    --lipids "$LIPIDS"
    --ratio "$LIPID_RATIO"
    --salt
    --salt_c "$CATION"
    --salt_a "$ANION"
    --saltcon "$SALT_CONC"
    --keepligs
    --notprotonate
    --dims "${BOX_DIMS[@]}"
    --ffwat "$WATER_MODEL"
    --ffprot "$PROTEIN_FF"
    --fflip "$LIPID_FF"
    --movebadrandom
    --short_penalty
    --nloop 40
    --parametrize
)
if [[ "$PREORIENTED" == true ]]; then
    memgen_command+=(--preoriented)
fi
if [[ -n "$TLEAP_LINE" ]]; then
    memgen_command+=(--leapline "$TLEAP_LINE")
fi
if [[ -n "$LIGAND_FRCMOD" ]]; then
    memgen_command+=(
        --gaff2
        --ligand_param "${LIGAND_FRCMOD}:${LIGAND_LIB}"
    )
fi

if $CHECK_ONLY; then
    echo "Configuration is valid. Planned commands:"
    echo "  Copy $INPUT_PDB to ${OUTPUT_DIR}/membrane_build/system.pdb"
    print_command "${memgen_command[@]}"
    echo "  ${MDTBX[*]} amb2gro -p ${OUTPUT_DIR}/membrane_build/bilayer_system_lipid.top -x ${OUTPUT_DIR}/membrane_build/bilayer_system_lipid.crd --type parmed -o $OUTPUT_DIR"
    echo "  ${MDTBX[*]} add_ndx, centering_gro, and protein/lipid gen_posres"
    echo "  Copy MDP files from $MDP_SOURCE_DIR"
    exit 0
fi

mkdir -p "${OUTPUT_DIR}/membrane_build"
cp "$INPUT_PDB" "${OUTPUT_DIR}/membrane_build/system.pdb"

pushd "${OUTPUT_DIR}/membrane_build" > /dev/null
"${memgen_command[@]}"
popd > /dev/null

"${MDTBX[@]}" amb2gro \
    -p "${OUTPUT_DIR}/membrane_build/bilayer_system_lipid.top" \
    -x "${OUTPUT_DIR}/membrane_build/bilayer_system_lipid.crd" \
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
protein_posres_files=("${OUTPUT_DIR}"/posres_*.itp)
shopt -u nullglob
[[ ${#protein_posres_files[@]} -gt 0 ]] ||
    fail "POSRES_SELECTION generated no restraint include files"
"${MDTBX[@]}" gen_posres \
    -p "${OUTPUT_DIR}/gmx.top" \
    -s "$LIPID_POSRES_SELECTION" \
    -o "${OUTPUT_DIR}/posres_lipid"
shopt -s nullglob
lipid_posres_files=("${OUTPUT_DIR}"/posres_lipid_*.itp)
shopt -u nullglob
[[ ${#lipid_posres_files[@]} -gt 0 ]] ||
    fail "LIPID_POSRES_SELECTION generated no restraint include files"

for stage in "${required_mdps[@]}"; do
    cp "${MDP_SOURCE_DIR}/${stage}.mdp" "${OUTPUT_DIR}/${stage}.mdp"
done

{
    printf 'workflow=membrane\n'
    printf 'input_pdb=%s\n' "$INPUT_PDB"
    printf 'mdp_source=%s\n' "$MDP_SOURCE_DIR"
    printf 'mdtbx_command=%s\n' "${MDTBX[*]}"
    printf 'lipids=%s\n' "$LIPIDS"
    printf 'lipid_ratio=%s\n' "$LIPID_RATIO"
    printf 'salt_concentration_molar=%s\n' "$SALT_CONC"
    printf 'ions=%s,%s\n' "$CATION" "$ANION"
    printf 'box_dimensions_angstrom=%s\n' "${BOX_DIMS[*]}"
    printf 'force_fields=%s,%s,%s\n' "$PROTEIN_FF" "$LIPID_FF" "$WATER_MODEL"
    printf 'preoriented=%s\n' "$PREORIENTED"
    printf 'center_group=%s\n' "$CENTER_SELECTION"
    printf 'posres_selection=%s\n' "$POSRES_SELECTION"
    printf 'lipid_posres_selection=%s\n' "$LIPID_POSRES_SELECTION"
    printf 'ligand_frcmod=%s\n' "$LIGAND_FRCMOD"
    printf 'ligand_lib=%s\n' "$LIGAND_LIB"
} > "${OUTPUT_DIR}/workflow.txt"

echo "Prepared $OUTPUT_DIR"
echo "Intermediate membrane-builder files are retained in ${OUTPUT_DIR}/membrane_build."
echo "Next: copy run_slurm.sh into $OUTPUT_DIR, edit it, and run --check."
