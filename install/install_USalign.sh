#!/bin/bash
set -e

thread=12
URL="https://github.com/pylelab/USalign.git"

git clone ${URL}
cd USalign
make -j $thread
./USalign -v

echo done
