#!/bin/bash
set -e

git clone https://github.com/haddocking/haddock-tools

cd haddock-tools
make -j 12


echo done
