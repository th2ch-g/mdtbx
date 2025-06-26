#!/bin/bash
set -e

thread=$1
PREFIX=$2
VERSION="2025.1"
SIMD="AVX_512"

# source /home/apps/Modules/init/bash
# module purge
# module load cmake gcc cuda

if [ -z $thread ] || [ -z $PREFIX ]; then
    echo "thread number and prefix must be inputed" >&2
    exit 1
fi

URL="https://ftp.gromacs.org/gromacs/gromacs-${VERSION}.tar.gz"
wget $URL
tar -zxvf gromacs-${VERSION}.tar.gz
rm -f gromacs-${VERSION}.tar.gz
cd gromacs-${VERSION}
mkdir build
cd build

# see compling option in https://manual.gromacs.org/current/install-guide/index.html

cmake .. \
    -DGMX_MPI=OFF \
    -DGMX_OPENMP=ON \
    -DGMX_GPU=CUDA \
    -DGMX_DOUBLE=OFF \
    -DCMAKE_INSTALL_PREFIX=$PREFIX \
    -DGMX_SIMD=$SIMD \
    -DGMX_FFT_LIBRARY=fftw3 \
    -DCMAKE_C_COMPILER=gcc \
    -DCMAKE_CXX_COMPILER=g++

make -j $thread
# make check -j $thread
make install -j $thread

echo done
