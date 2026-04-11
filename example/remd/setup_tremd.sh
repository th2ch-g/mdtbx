#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"

TEMPLATE_MDP="${TEMPLATE_MDP:-${REPO_DIR}/example/mdp/solution/tremd.mdp}"
SUBMIT_SCRIPT="${SUBMIT_SCRIPT:-${REPO_DIR}/example/mdrun/remd_slurm.sh}"
STRUCTURE="${STRUCTURE:-gmx.gro}"
TOPOLOGY="${TOPOLOGY:-gmx.top}"
INDEX="${INDEX:-index.ndx}"
ITP_GLOB="${ITP_GLOB:-*.itp}"
DEFFNM="${DEFFNM:-tremd}"
MAXWARN="${MAXWARN:-10}"
REPLEX="${REPLEX:-500}"
TEMPERATURES=(310 320 330 340 350 360 370 380 390 400)

if [ ! -f "${STRUCTURE}" ]; then
    echo "Structure file not found: ${STRUCTURE}" >&2
    exit 1
fi

if [ ! -f "${TOPOLOGY}" ]; then
    echo "Topology file not found: ${TOPOLOGY}" >&2
    exit 1
fi

if [ ! -f "${TEMPLATE_MDP}" ]; then
    echo "Template MDP not found: ${TEMPLATE_MDP}" >&2
    exit 1
fi

N_REPLICA=${#TEMPERATURES[@]}
echo "Preparing ${N_REPLICA} T-REMD replicas"
echo "Temperatures: ${TEMPERATURES[*]}"

shopt -s nullglob
itp_files=(${ITP_GLOB})
shopt -u nullglob

for idx in "${!TEMPERATURES[@]}"; do
    rep=$((idx + 1))
    temp="${TEMPERATURES[$idx]}"
    repdir="rep${rep}"

    mkdir -p "${repdir}"
    cp "${STRUCTURE}" "${repdir}/gmx.gro"
    cp "${TOPOLOGY}" "${repdir}/gmx.top"
    cp "${TEMPLATE_MDP}" "${repdir}/${DEFFNM}.mdp"

    if [ -f "${INDEX}" ]; then
        cp "${INDEX}" "${repdir}/index.ndx"
    fi

    if [ ${#itp_files[@]} -gt 0 ]; then
        cp "${itp_files[@]}" "${repdir}/"
    fi

    mdtbx mod_mdp --path "${repdir}" -t ref_t -v "${temp}"

    grompp_cmd=(
        gmx_mpi grompp
        -f "${repdir}/${DEFFNM}.mdp"
        -c "${repdir}/gmx.gro"
        -p "${repdir}/gmx.top"
        -maxwarn "${MAXWARN}"
        -o "${repdir}/${DEFFNM}.tpr"
    )

    if [ -f "${repdir}/index.ndx" ]; then
        grompp_cmd+=(-n "${repdir}/index.ndx")
    fi

    "${grompp_cmd[@]}"
done

if [ -f "${SUBMIT_SCRIPT}" ]; then
    submit_copy="$(basename "${SUBMIT_SCRIPT}")"
    awk -v n_replica="${N_REPLICA}" -v deffnm="${DEFFNM}" -v replex="${REPLEX}" '
        /^REPLEX=/ { print "REPLEX=" replex; next }
        /^N_REPLICA=/ { print "N_REPLICA=" n_replica; next }
        /^DEFFNM=/ { print "DEFFNM=\"" deffnm "\" # or otherremd"; next }
        { print }
    ' "${SUBMIT_SCRIPT}" > "${submit_copy}"
    echo "Wrote submit script: ${submit_copy}"
fi

echo "Done"
