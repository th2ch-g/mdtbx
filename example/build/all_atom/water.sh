#!/bin/bash
set -e

out_dir="${PWD}/gmx"

mdtbx build_solution --noions

mdtbx amb2gro -p leap.parm7 -x leap.rst7 --type parmed

mdtbx add_ndx -g gmx.gro

mdtbx rmfile

mkdir gmx
mv gmx.gro ${out_dir}/
mv gmx.top ${out_dir}/
mv *.ndx ${out_dir}/
cp mdps/*.mdp ${out_dir}/
cp mdrun_slurm.sh ${out_dir}/

rm -f leap.parm7 leap.rst7 leap.pdb gmx.pdb

echo done
