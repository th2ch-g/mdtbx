#!/bin/bash
set -e

mdtbx cmd python3 calc_dg_points_mbar.py 2.6 2.93 pmf_mbar.dat

echo done
