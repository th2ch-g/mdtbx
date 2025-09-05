#!/bin/bash
set -e

VERSION="4.5.6"
INSTALL_PREFIX="${PWD}/dssp-${VERSION}/build"
thread=12

URL="https://github.com/PDB-REDO/dssp/archive/refs/tags/v${VERSION}.tar.gz"
wget $URL
tar -xvzf v${VERSION}.tar.gz
cd dssp-${VERSION}

cmake -S . -B build -DBUILD_PYTHON_MODULE=ON -DCMAKE_INSTALL_PREFIX=$INSTALL_PREFIX
cmake --build build -DCMAKE_INSTALL_PREFIX=$INSTALL_PREFIX
cmake --install build -DCMAKE_INSTALL_PREFIX=$INSTALL_PREFIX

echo done
