#!/usr/bin/env bash
set -euo pipefail

# Copy this file next to run/ and edit this block.
MDTBX=(mdtbx)
# From a source checkout, use:
# MDTBX=(pixi run --manifest-path /path/to/mdtbx/pyproject.toml mdtbx)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
RUN_DIR="${RUN_DIR:-${SCRIPT_DIR}/run}"
OUTPUT_DIR="${OUTPUT_DIR:-analysis}"
PRODUCTION_PREFIX="${PRODUCTION_PREFIX:-prd}"
PRODUCTION_SEGMENTS="${PRODUCTION_SEGMENTS:-10}"
INDEX_FILE="${INDEX_FILE:-index.ndx}"
# GROMACS index group names.
KEEP_SELECTION="${KEEP_SELECTION:-non-Water}"
CENTER_SELECTION="${CENTER_SELECTION:-Protein}"
# MDTraj atom-selection expressions.
ANALYSIS_SELECTION="${ANALYSIS_SELECTION:-protein and name CA}"
CONTACT_SELECTION="${CONTACT_SELECTION:-protein and name CA}"

CHECK_ONLY=false
case "${1:-}" in
    "") ;;
    --check) CHECK_ONLY=true ;;
    *) echo "Usage: $0 [--check]" >&2; exit 2 ;;
esac

cd "$SCRIPT_DIR"

fail() {
    echo "ERROR: $*" >&2
    exit 1
}

require_file() {
    [[ -f "$1" ]] || fail "File not found: $1"
}

[[ "$PRODUCTION_SEGMENTS" =~ ^[1-9][0-9]*$ ]] || \
    fail "PRODUCTION_SEGMENTS must be a positive integer"
command -v "${MDTBX[0]}" >/dev/null || fail "Command not found: ${MDTBX[0]}"
[[ -d "$RUN_DIR" ]] || fail "Run directory not found: $RUN_DIR"
require_file "${RUN_DIR}/${INDEX_FILE}"
for segment in $(seq 1 "$PRODUCTION_SEGMENTS"); do
    require_file "${RUN_DIR}/${PRODUCTION_PREFIX}${segment}.tpr"
    require_file "${RUN_DIR}/${PRODUCTION_PREFIX}${segment}.xtc"
done

if $CHECK_ONLY; then
    echo "Configuration is valid. Planned analysis:"
    echo "  ${MDTBX[*]} trjcat -> ${OUTPUT_DIR}/fitted.xtc"
    echo "  ${MDTBX[*]} rmsd, rmsf, contactmap, and print_perf"
    exit 0
fi

cd "$RUN_DIR"
mkdir -p "$OUTPUT_DIR"

"${MDTBX[@]}" trjcat \
    -n "$PRODUCTION_SEGMENTS" \
    --prefix "$PRODUCTION_PREFIX" \
    -idx "$INDEX_FILE" \
    -k "$KEEP_SELECTION" \
    -c "$CENTER_SELECTION" \
    -o "${OUTPUT_DIR}/fitted.xtc" \
    --preserve-time \
    --keep-concatenated

"${MDTBX[@]}" rmsd \
    -p rmmol_top.gro \
    -t "${OUTPUT_DIR}/fitted.xtc" \
    -r rmmol_top.gro \
    -sct "$ANALYSIS_SELECTION" \
    -scr "$ANALYSIS_SELECTION" \
    -sft "$ANALYSIS_SELECTION" \
    -sfr "$ANALYSIS_SELECTION" \
    -o "${OUTPUT_DIR}/rmsd.npy"

"${MDTBX[@]}" rmsf \
    -p rmmol_top.gro \
    -t "${OUTPUT_DIR}/fitted.xtc" \
    --selection "$ANALYSIS_SELECTION" \
    --resolution residue \
    -o "${OUTPUT_DIR}/rmsf.npy"

"${MDTBX[@]}" contactmap \
    -p rmmol_top.gro \
    -t "${OUTPUT_DIR}/fitted.xtc" \
    -s "$CONTACT_SELECTION" \
    -o "${OUTPUT_DIR}/contactmap.npy"

"${MDTBX[@]}" print_perf \
    -p . \
    --prefix "$PRODUCTION_PREFIX" \
    -o "${OUTPUT_DIR}/performance.csv"

echo "Analysis written to ${RUN_DIR}/${OUTPUT_DIR}"
