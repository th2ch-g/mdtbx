#!/bin/bash
set -e

VERSION="2025.1"
INSTALL_PREFIX="${PWD}/gromacs-${VERSION}/build"
thread=12
SIMD="AVX_512"

FFT_LIB="/path/to/fftw-3.3.10/lib/libfftw3.la"
FFT_INCLUDE="/path/to/fftw-3.3.10/include"
BLRS_USER="/path/to/openblas-0.3.21/lib/libopenblas_zenp-r0.3.21.a"
LAPACK_USER="/path/to/scalapack-2.2.1/lib/libscalapack.a"
CP2K_DIR="/path/to/cp2k/2023.1/lib/rccs/psmp/"

# source /home/apps/Modules/init/bash
# module purge
# module load cmake gcc cuda

URL="https://ftp.gromacs.org/gromacs/gromacs-${VERSION}.tar.gz"
wget $URL
tar -zxvf gromacs-${VERSION}.tar.gz
rm -f gromacs-${VERSION}.tar.gz
cd gromacs-${VERSION}
mkdir build
cd build

# see compling option in https://manual.gromacs.org/current/install-guide/index.html

cmake .. \
    -DGMX_MPI=ON \
    -DCMAKE_INSTALL_PREFIX=$PREFIX \
    -DGMX_GPU=OFF \
    -DGMX_FFT_LIBRARY=fftw3 \
    -DGMX_BUILD_OWN_FFTW=ON \
    -DCMAKE_C_COMPILER=mpicc \
    -DCMAKE_CXX_COMPILER=mpic++ \
    -DBUILD_SHARED_LIBS=OFF \
    -DGMXAPI=OFF \
    -DGMX_INSTALL_NBLIB_API=OFF \
    -DGMX_DOUBLE=ON \
    -DGMX_FFT_LIBRARY=fftw3  \
    -DFFTWF_LIBRARY=$FFT_LIB  \
    -DFFTWF_INCLUDE_DIR=$FFT_INCLUDE  \
    -DGMX_BLAS_USER=$BLRS_USER  \
    -DGMX_LAPACK_USER=$LAPACK_USER \
    -DGMX_CP2K=ON \
    -DCP2K_DIR=$CP2K_DIR

make -j $thread
# make check -j $thread
make install -j $thread

echo done
