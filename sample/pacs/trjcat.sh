#!/bin/bash
#SBATCH -p all_q
#SBATCH -n 1
set -e

echo $(hostname)

ref_tpr="trial1/trial001/cycle000/replica001/rmmol_top.tpr"

for i in {1..15};
do
    mdtbx pacs_trjcat -t trial$i/trial001 -r $ref_tpr --index index.ndx -k Protein_NKP
done

echo done

