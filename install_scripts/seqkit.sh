#!/bin/bash
set -e

VERSION="v2.10.0"
OS_CPU="linux_amd64"

URL="https://github.com/shenwei356/seqkit/releases/download/${VERSION}/seqkit_${OS_CPU}.tar.gz"
wget $URL
tar -zxvf seqkit_${OS_CPU}.tar.gz
rm -f seqkit_${OS_CPU}.tar.gz

echo done
