#!/bin/bash

# Move into script's directory
cd "$(realpath "$(dirname -- "${BASH_SOURCE[0]}")")"

HELPMSG="Usage: $(basename "$0") [-h] [-x] [-g]
Builds the gencap library in the ./bin/ directory. Optionally toggle building
of the testing executable, and building in debug mode for development.

Launch Options:
HELP:
  -h  shows this help text

TEST EXECUTABLE:
  -x  the script will also install the test executable that demonstrates the
      library's capabilities

DEBUGGING MODE:
  -g  enables the debugging flags 'g', 'O0', 'Wall', and 'fbounds-check' for
      gfortran"

unset do_exe
while getopts hxg FLAG; do
    case "${FLAG}" in
        h)
            echo "$HELPMSG"
            exit 0
            ;;
        x) do_exe=true;;
        g) export debug_mode=true;;
    esac
done

# Making library
echo "Making the library..."
if make; then
    echo "Completed library!"
else
    echo "Failed to make the library."
    exit 2
fi
echo

if [[ "${do_exe}" == true ]]; then
	# Making test executable
	echo "Making the test executable..."
	if make gentest.x; then
        echo "Completed executable!"
    else
        echo "Failed to make the executable."
        exit 3
    fi
fi