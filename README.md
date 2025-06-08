[![DOI](https://zenodo.org/badge/103404474.svg)](https://zenodo.org/badge/latestdoi/103404474)
# captngen

Capt'n General: $q^{2n}$ and $v^{2m}$ dependent solar dark matter capture and energy transport routines.

Capt'n Oper: Solar capture using the NREO formalism as adapted from [https://arxiv.org/abs/1501.03729]().

## Installation

Run the `install.sh` script to install the `gencap` library and an optional testing executable.
Use `install.sh -h` to learn more.
Run the `cleanup.sh` script to clear all installed files and temporary `.dat` files produced by the testing executable.

## Development

Can be built using `make` as a library by default (or by explicitly calling `make libgencap.so`), or can be built as a standalone executable using `make gentest.x`.
To enable (or disable) debugging, set (or unset) the shell variable `debug` to any non-empty value:

```shell
export debug=foo
unset debug
```

This can be done for a single `make` call by adding the assignment after the make target(s): `make bar debug=foo`.
Make sure to clear all build files when switching to and from debugging mode to endure all object files are compiled with the same flags!

See `main.f90` for examples of how to call the executable, `gentest.x`.

## Citing Captn General

If you use this code, you can cite [https://arxiv.org/abs/2105.06810]() and/or [https://arxiv.org/abs/1808.10465]() where it was first deployed.
