#!/bin/bash

# Move into script's directory
cd "$(realpath "$(dirname -- "${BASH_SOURCE[0]}")")"

# Nuke all built objects, modules, libraries, and executables
echo "Nuking built files..."
make nuke
echo

# Removing data files created by the test executable
echo "Deleting temperary *.dat files from $(pwd)..."
rm -f *.dat
