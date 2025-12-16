#!/bin/bash
set -e
brew install --cask blender
URL=https://github.com/BradyAJohnston/MolecularNodes/releases/download/v4.5.9/molecularnodes-4.5.9-macos_arm64.zip
wget $URL
echo done
