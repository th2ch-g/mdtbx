#!/bin/bash
set -e

thread=$1
PREFIX=$2
VERSION="3-1-1"

# source /home/apps/Modules/init/bash
# module purge
# module load cmake gcc

if [ -z $thread ] || [ -z $PREFIX ]; then
    echo "thread number and prefix must be inputed" >&2
    exit 1
fi

URL="https://github.com/openbabel/openbabel/archive/refs/tags/openbabel-${VERSION}.tar.gz"
wget $URL
tar -zxvf openbabel-${VERSION}.tar.gz
rm -f openbabel-${VERSION}.tar.gz
cd openbabel-${VERSION}
mkdir build
cd build

# see compling options in https://github.com/openbabel/openbabel/blob/master/INSTALL

cmake .. -DCMAKE_INSTALL_PREFIX=$PREFIX

make -j $thread
make install -j $thread

echo done
