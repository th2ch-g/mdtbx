#!/bin/bash
set -e

mdtbx cmd python3 calc_mbar.py
mdtbx cmd python3 plot_mbar.py pmf_output.dat

echo done
