#!/bin/bash
set -e

thread=12
INSTALL_PREFIX="${PWD}/obabel-${VERSION}"
VERSION="3-1-1"

# source /home/apps/Modules/init/bash
# module purge
# module load cmake gcc

URL="https://github.com/openbabel/openbabel/archive/refs/tags/openbabel-${VERSION}.tar.gz"
wget $URL
tar -zxvf openbabel-${VERSION}.tar.gz
rm -f openbabel-${VERSION}.tar.gz
cd openbabel-${VERSION}
mkdir build
cd build

# see compling options in https://github.com/openbabel/openbabel/blob/master/INSTALL

cmake .. -DCMAKE_INSTALL_PREFIX=$INSTALL_PREFIX

make -j $thread
make install -j $thread

echo done
