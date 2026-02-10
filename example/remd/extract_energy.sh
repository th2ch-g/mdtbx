#!/bin/bash
set -e

N_REPLICA=16
PREFIX="grest"

for rep in $(seq 1 $N_REPLICA);
do
    echo "Potential" | gmx energy -f rep${rep}/${PREFIX}.edr -o u_rep${rep}.xvg
done

echo done
