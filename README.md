[![DOI](https://zenodo.org/badge/103404474.svg)](https://zenodo.org/badge/latestdoi/103404474)
# captngen

Capt'n General: $q^{2n}$ and $v^{2m}$ dependent solar dark matter capture and energy transport routines.

Capt'n Oper: Solar capture using the NREO formalism as adapted from [https://arxiv.org/abs/1501.03729]().

Can be built using `make` as a library by default (or by explicitly calling `make gencaplib.so`), or can be built as a standalone executable using `make gentest.x`.

See `main.f90` for examples of how to call the executable, `gentest.x`.

If you use this code, you can cite [https://arxiv.org/abs/2105.06810]() and/or [https://arxiv.org/abs/1808.10465]() where it was first deployed.

## fftpack

To get the [fftpack](https://github.com/fortran-lang/fftpack) submodule working from a fresh `git clone`, run the commands:

```shell
git submodule init
git submodule update
```

Alternatively, during the initial clone, pass the `--recurse-submodules` flag to clone the submodule automatically:

```shell
git clone --recurse-submodules https://github.com/aaronvincent/captngen.git
```
