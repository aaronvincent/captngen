#!/bin/bash

# Move into script's directory
cd "$(realpath "$(dirname -- "${BASH_SOURCE[0]}")")"

# Nuke all built objects, modules, libraries, and executables
echo "Nuking built files..."
make nuke
echo

# Remove the empty build directories
echo "Removing build directories from $(pwd)..."
if [[ -d obj ]]; then rmdir obj/* && rmdir obj; fi
if [[ -d lib ]]; then rmdir lib; fi
if [[ -d bin ]]; then rmdir bin; fi
echo

# Nuke all built debug objects, modules, libraries, and executables
echo "Nuking built debug files..."
make nuke debug=true
echo

# Remove the empty debug build directories
echo "Removing build directories from $(pwd)..."
if [[ -d obj-debug ]]; then rmdir obj-debug/* && rmdir obj-debug; fi
if [[ -d lib-debug ]]; then rmdir lib-debug; fi
if [[ -d bin-debug ]]; then rmdir bin-debug; fi
echo

# Removing data files created by the test executable
echo "Deleting temporary *.dat files from $(pwd)..."
rm -f *.dat
