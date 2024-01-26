#!/bin/bash

# Move into script's directory
cd "$(realpath "$(dirname -- "${BASH_SOURCE[0]}")")"

HELPMSG="Usage: $(basename "$0") [-h] [-g]

Launch Options:
HELP:
  -h  shows this help text

DEBUGGING MODE:
  -g  enables the debugging flags 'g', 'O0', 'Wall', and 'fbounds-check' for gfortran"

while getopts hg FLAG; do
        case "${FLAG}" in
                h)
                        echo "$HELPMSG"
                        exit 0
                        ;;
                g) export debug_mode=true;;
        esac
done

# Making library
echo "Making the library..."
make && echo "Completed library!" || echo "Failed to make the library."
echo

# Making test executable
echo "Making the test executable..."
make gentest.x && echo "Completed executable!" || echo "Failed to make the executable."
