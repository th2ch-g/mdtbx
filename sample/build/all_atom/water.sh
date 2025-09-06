#!/bin/bash
set -e

mdtbx build_solution --noions

mdtbx amb2gro -p leap.parm7 -x leap.rst7 --type parmed

mdtbx add_ndx -g gmx.gro

mdtbx rmfile

mkdir gmx
mv gmx.gro gmx/
mv gmx.top gmx/
mv *.ndx gmx/
cp mdps/*.mdp gmx/
cp mdrun_slurm.sh gmx/

rm -f leap.parm7 leap.rst7 leap.pdb gmx.pdb

echo done
