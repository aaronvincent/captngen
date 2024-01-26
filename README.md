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

To enable debugging, set the shell variable `debug_mode` to `true`:

```shell
export debug_mode=true
```

To disable debugging, simply unset `debug_mode`:

```shell
unset debug_mode
```

See `main.f90` for examples of how to call the executable, `gentest.x`.

## Citing Captn General

If you use this code, you can cite [https://arxiv.org/abs/2105.06810]() and/or [https://arxiv.org/abs/1808.10465]() where it was first deployed.
