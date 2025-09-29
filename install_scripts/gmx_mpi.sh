#!/bin/bash
set -eux

thread=12
VERSION="2025.1"
PREFIX="$PWD/gromacs-${VERSION}/"

# source /home/apps/Modules/init/bash
# module purge
# module load gcc/13.3.0 cuda/12.9 cmake/3.31.6 openmpi/5.0.7
# module load $TOOLS/hpc_sdk/modulefiles/nvhpc/25.7


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
    -DGMX_OPENMP=ON \
    -DGMX_GPU=CUDA \
    -DGMX_DOUBLE=OFF \
    -DCMAKE_INSTALL_PREFIX=$PREFIX \
    -DGMX_FFT_LIBRARY=fftw3 \
    -DGMX_BUILD_OWN_FFTW=ON \
    -DCMAKE_BUILD_TYPE=Release \
    -DGMX_USE_CUFFTMP=ON \
    -DcuFFTMp_ROOT="$TOOLS/hpc_sdk/Linux_x86_64/25.7/math_libs/12.9" \
    -DCMAKE_C_COMPILER=mpicc \
    -DCMAKE_CXX_COMPILER=mpic++
    # -DGMX_NVSHMEM=ON \
    # -DNVSHMEM_ROOT=<Path-to-NVSHMEM-Lib-Root-dir> \

make -j $thread
# make check -j $thread
make install -j $thread

echo done
