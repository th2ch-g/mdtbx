#!/bin/bash
set -ex

PREFIX="$PWD/OpenCafeMol"

git clone --recurse-submodules https://github.com/yutakasi634/OpenCafeMol.git
cd OpenCafeMol
mkdir build
cd build
cmake .. -DOPENMM_ROOT=$TOOLS/mdtbx/.pixi/envs/default -DCMAKE_INSTALL_PREFIX=$PREFIX
make -j 12
make install

echo done
