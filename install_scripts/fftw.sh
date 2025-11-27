#!/bin/bash
set -e

VERSION="3.3.10"
URL="https://www.fftw.org/fftw-${VERSION}.tar.gz"
PREFIX=$PWD/fftw-${VERSION}/build

wget $URL
tar -zxvf fftw-${VERSION}.tar.gz
cd fftw-${VERSION}
./configure --prefix=$PREFIX CC=gcc MPICC=mpicc F77=gfortran --enable-mpi --enable-threads
make -j 10
make install

echo done
