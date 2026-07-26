#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"

TEMPLATE_MDP="${TEMPLATE_MDP:-${REPO_DIR}/example/mdp/solution/re.mdp}"
TEMPLATE_PLUMED="${TEMPLATE_PLUMED:-${REPO_DIR}/example/plumed/reus_dist.dat}"
SUBMIT_SCRIPT="${SUBMIT_SCRIPT:-${REPO_DIR}/example/mdrun/remd_plumed_slurm.sh}"
TARGETS_FILE="${TARGETS_FILE:-${SCRIPT_DIR}/reus_distances.txt}"
TOPOLOGY="${TOPOLOGY:-gmx.top}"
INDEX="${INDEX:-index.ndx}"
ITP_GLOB="${ITP_GLOB:-*.itp}"
GMX_CMD="${GMX_CMD:-gmx_mpi}"
DEFFNM="${DEFFNM:-reus}"
MAXWARN="${MAXWARN:-10}"
REPLEX="${REPLEX:-1000}"

for required_file in \
    "${TEMPLATE_MDP}" \
    "${TEMPLATE_PLUMED}" \
    "${TARGETS_FILE}" \
    "${TOPOLOGY}" \
    "${INDEX}"; do
    if [ ! -f "${required_file}" ]; then
        echo "Required file not found: ${required_file}" >&2
        exit 1
    fi
done

N_REPLICA="$(
    awk 'NF >= 2 && $1 !~ /^#/ { count++ } END { print count + 0 }' "${TARGETS_FILE}"
)"
if [ "${N_REPLICA}" -eq 0 ]; then
    echo "No replica definitions found in ${TARGETS_FILE}" >&2
    exit 1
fi

echo "Preparing ${N_REPLICA} REUS replicas"

shopt -s nullglob
# shellcheck disable=SC2206
itp_files=(${ITP_GLOB})
shopt -u nullglob

rep=0
while read -r structure target_distance; do
    rep=$((rep + 1))
    repdir="rep${rep}"

    if [ ! -f "${structure}" ]; then
        echo "Structure file not found: ${structure}" >&2
        exit 1
    fi
    if ! [[ "${target_distance}" =~ ^[+]?[0-9]+([.][0-9]+)?$ ]]; then
        echo "Invalid target distance: ${target_distance}" >&2
        exit 1
    fi

    mkdir -p "${repdir}"
    cp "${structure}" "${repdir}/gmx.gro"
    cp "${TOPOLOGY}" "${repdir}/gmx.top"
    cp "${INDEX}" "${repdir}/index.ndx"
    cp "${TEMPLATE_MDP}" "${repdir}/${DEFFNM}.mdp"
    sed "s/TARGET_DISTANCE/${target_distance}/g" \
        "${TEMPLATE_PLUMED}" > "${repdir}/plumed.dat"

    if [ ${#itp_files[@]} -gt 0 ]; then
        cp "${itp_files[@]}" "${repdir}/"
    fi

    "${GMX_CMD}" grompp \
        -f "${repdir}/${DEFFNM}.mdp" \
        -c "${repdir}/gmx.gro" \
        -n "${repdir}/index.ndx" \
        -p "${repdir}/gmx.top" \
        -maxwarn "${MAXWARN}" \
        -o "${repdir}/${DEFFNM}.tpr"
done < <(awk 'NF >= 2 && $1 !~ /^#/ { print $1, $2 }' "${TARGETS_FILE}")

if [ -f "${SUBMIT_SCRIPT}" ]; then
    submit_copy="${SUBMIT_OUTPUT:-$(basename "${SUBMIT_SCRIPT}")}"
    submit_tmp="$(mktemp "${submit_copy}.tmp.XXXXXX")"
    awk -v n_replica="${N_REPLICA}" -v deffnm="${DEFFNM}" -v replex="${REPLEX}" '
        /^REPLEX=/ { print "REPLEX=" replex; next }
        /^N_REPLICA=/ { print "N_REPLICA=" n_replica; next }
        /^DEFFNM=/ { print "DEFFNM=\"" deffnm "\" # or rest2, reus"; next }
        { print }
    ' "${SUBMIT_SCRIPT}" > "${submit_tmp}"
    chmod +x "${submit_tmp}"
    mv "${submit_tmp}" "${submit_copy}"
    echo "Wrote submit script: ${submit_copy}"
fi

echo "Done"
