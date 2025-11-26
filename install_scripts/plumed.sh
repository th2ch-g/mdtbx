#!/bin/bash
set -e

VERSION="2.10.0"
URL="https://github.com/plumed/plumed2/releases/download/v${VERSION}/plumed-src-${VERSION}.tgz"
PREFIX="$PWD/plumed-${VERSION}/build"

wget $URL
tar -xvzf plumed-src-${VERSION}.tgz
rm -rf plumed-src-${VERSION}.tgz

cd plumed-${VERSION}/
./configure --prefix=$PREFIX
make -j 10
make install

echo done
