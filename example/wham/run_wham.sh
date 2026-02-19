#!/bin/bash
set -e

ls -1 rep*/*.tpr > tmp_tpr.dat
ls -1 rep*/*_pullf.xvg > tmp_pullf.dat
ls -1 rep*/*_pullx.xvg > tmp_pullx.dat

# for angle
mdtbx cmd gmx wham \
    -it tmp_tpr.dat \
    -ix tmp_pullx.dat \ # or -if tmp_pullf.dat
    -o profile.xvg \
    -hist hist.xvg \
    -unit kCal \
    -temp 310 \
    -nBootstrap 100 \
    -min -180 \
    -max 180 \
    -cycl

mdtbx cmd python3 plot_hist.py
mdtbx cmd python3 plot_pmf.py
mdtbx cmd python3 plot_bspmf.py

rm -f tmp_tpr.dat tmp_pullf.dat tmp_pullx.dat

echo done
