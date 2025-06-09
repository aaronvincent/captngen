#!/bin/bash

# Move into script's directory
cd "$(realpath "$(dirname -- "${BASH_SOURCE[0]}")")"

# Nuke all built objects, modules, libraries, and executables
echo "Nuking built files..."
make nuke
echo

# Nuke all built debug objects, modules, libraries, and executables
echo "Nuking built debug files..."
make nuke debug=true
echo

# Removing data files created by the test executable
echo "Deleting temporary *.dat files from $(pwd)..."
rm -f *.dat
