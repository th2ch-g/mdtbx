#!/bin/bash
set -e

git clone https://github.com/haddocking/haddock3.git
cd haddock3
uv venv
uv pip install -e "."

echo done
