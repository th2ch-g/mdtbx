#!/bin/bash
set -e

VERSION="1.2.0"
URL="https://github.com/sib-swiss/termal/releases/download/v${VERSION}/termal-v${VERSION}-x86_64-unknown-linux-gnu.tar.gz"

wget ${URL}
mkdir ./tmp
tar -xvzf termal-v${VERSION}-x86_64-unknown-linux-gnu.tar.gz -C ./tmp

mv ./tmp/termal .

rm -rf termal-v${VERSION}-x86_64-unknown-linux-gnu.tar.gz ./tmp

echo done
