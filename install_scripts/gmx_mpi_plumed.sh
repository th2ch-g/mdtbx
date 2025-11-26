#!/bin/bash
set -eux

thread=12
VERSION="2024.3"
PREFIX="$PWD/gromacs-${VERSION}/"

PLUMED_VERSION="2.10.0"

source /home/apps/Modules/init/bash
module purge
module load gcc/13.3.0 cuda/12.9 cmake/3.31.6 openmpi/5.0.7
module load $TOOLS/plumed-${PLUMED_VERSION}/build/lib/plumed/modulefile

# export LD_LIBRARY_PATH=$TOOLS/plumed-2.10.0/build/lib:$LD_LIBRARY_PATH
# export PATH=$TOOLS/plumed-2.10.0/build/bin:$PATH

URL="https://ftp.gromacs.org/gromacs/gromacs-${VERSION}.tar.gz"
wget $URL
tar -zxvf gromacs-${VERSION}.tar.gz
rm -f gromacs-${VERSION}.tar.gz
cd gromacs-${VERSION}
# echo 4 | plumed patch -p # select 4: gromacs-2025.0
echo 3 | plumed patch -p # select 4: gromacs-2024.3
mkdir build
cd build

# see compling option in https://manual.gromacs.org/current/install-guide/index.html

cmake .. \
    -DGMX_MPI=ON \
    -DGMX_OPENMP=ON \
    -DGMX_GPU=CUDA \
    -DGMX_DOUBLE=OFF \
    -DGMX_THREAD_MPI=OFF \
    -DCMAKE_INSTALL_PREFIX=$PREFIX \
    -DGMX_FFT_LIBRARY=fftw3 \
    -DGMX_BUILD_OWN_FFTW=ON \
    -DCMAKE_BUILD_TYPE=Release \
    -DGMX_USE_PLUMED=ON \
    -DCMAKE_C_COMPILER=mpicc \
    -DCMAKE_CXX_COMPILER=mpic++

make -j $thread
# make check -j $thread
make install -j $thread

echo done
