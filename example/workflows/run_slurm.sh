#!/usr/bin/env bash
#SBATCH --job-name=md
#SBATCH --partition=CHANGE_ME
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --time=24:00:00
#SBATCH --output=slurm-%j.out

set -euo pipefail

# Copy this file into run/, then edit this block for the target cluster.
GROMACS_BIN="${GROMACS_BIN:-/path/to/installer/gmx}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
SYSTEM_DIR="${SYSTEM_DIR:-$SCRIPT_DIR}"
MAXWARN="${MAXWARN:-0}"
PRODUCTION_PREFIX="${PRODUCTION_PREFIX:-prd}"
PRODUCTION_SEGMENTS="${PRODUCTION_SEGMENTS:-10}"
MDRUN_LAUNCHER=(srun)
MDRUN_OPTIONS=(-ntomp "${SLURM_CPUS_PER_TASK:-1}")

CHECK_ONLY=false
case "${1:-}" in
    "") ;;
    --check) CHECK_ONLY=true ;;
    *) echo "Usage: $0 [--check]" >&2; exit 2 ;;
esac

fail() {
    echo "ERROR: $*" >&2
    exit 1
}

require_file() {
    [[ -f "$1" ]] || fail "File not found: $1"
}

[[ "$GROMACS_BIN" != /path/to/* ]] || fail "Edit GROMACS_BIN"
[[ -x "$GROMACS_BIN" ]] || fail "GROMACS_BIN is not executable: $GROMACS_BIN"
GROMACS_REALPATH="$(realpath "$GROMACS_BIN")" || \
    fail "Cannot resolve GROMACS_BIN: $GROMACS_BIN"
case "${GROMACS_BIN}:${GROMACS_REALPATH}" in
    *"/.pixi/envs/"*) fail "Use installer GROMACS, not pixi GROMACS" ;;
esac
[[ "$PRODUCTION_SEGMENTS" =~ ^[1-9][0-9]*$ ]] || \
    fail "PRODUCTION_SEGMENTS must be a positive integer"
[[ -d "$SYSTEM_DIR" ]] || fail "System directory not found: $SYSTEM_DIR"

for file in gmx.gro gmx.top index.ndx mini.mdp \
    eq1.mdp eq2.mdp eq3.mdp eq4.mdp eq5.mdp eq6.mdp prd.mdp; do
    require_file "${SYSTEM_DIR}/${file}"
done

if $CHECK_ONLY; then
    echo "Configuration is valid. Planned stages:"
    echo "  mini eq1 eq2 eq3 eq4 eq5 eq6"
    echo "  ${PRODUCTION_PREFIX}1..${PRODUCTION_PREFIX}${PRODUCTION_SEGMENTS}"
    echo "GROMACS: $GROMACS_BIN"
    echo "This script does not submit itself. Use sbatch only after reviewing it."
    exit 0
fi

[[ -n "${SLURM_JOB_ID:-}" ]] || fail "Run this script inside a Slurm allocation"
cd "$SYSTEM_DIR"

run_stage() {
    local stage="$1"
    local mdp="$2"
    local previous_structure="$3"
    local previous_checkpoint="${4:-}"

    if [[ -f "${stage}.gro" ]]; then
        echo "Skipping completed stage: $stage"
        return
    fi

    if [[ ! -f "${stage}.tpr" ]]; then
        if [[ -n "$previous_checkpoint" ]]; then
            [[ -f "$previous_checkpoint" ]] || \
                fail "Checkpoint required for continuation: $previous_checkpoint"
        fi
        grompp_command=(
            "$GROMACS_BIN" grompp
            -f "$mdp"
            -c "$previous_structure"
            -r gmx.gro
            -p gmx.top
            -n index.ndx
            -o "${stage}.tpr"
            -maxwarn "$MAXWARN"
        )
        if [[ -n "$previous_checkpoint" ]]; then
            grompp_command+=(-t "$previous_checkpoint")
        fi
        "${grompp_command[@]}"
    fi

    mdrun_command=(
        "${MDRUN_LAUNCHER[@]}"
        "$GROMACS_BIN" mdrun
        -deffnm "$stage"
        "${MDRUN_OPTIONS[@]}"
    )
    if [[ -f "${stage}.cpt" ]]; then
        mdrun_command+=(-cpi "${stage}.cpt" -append)
    fi
    "${mdrun_command[@]}"
}

run_stage mini mini.mdp gmx.gro
run_stage eq1 eq1.mdp mini.gro
previous_stage=eq1
for eq_stage in eq2 eq3 eq4 eq5 eq6; do
    run_stage \
        "$eq_stage" \
        "${eq_stage}.mdp" \
        "${previous_stage}.gro" \
        "${previous_stage}.cpt"
    previous_stage="$eq_stage"
done

for segment in $(seq 1 "$PRODUCTION_SEGMENTS"); do
    stage="${PRODUCTION_PREFIX}${segment}"
    run_stage \
        "$stage" \
        prd.mdp \
        "${previous_stage}.gro" \
        "${previous_stage}.cpt"
    previous_stage="$stage"
done
