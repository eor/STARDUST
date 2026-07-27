# STARDUST

STARDUST, which stands for **S**pectrum **T**ransport **A**round **R**a**D**iating **U**niversal **S**ource **T**ype,
is a cosmological 1D radiative transfer code described in detail by [Krause et al. (2018)][1].
An earlier (unpublished) version of the code was used by [Thomas & Zaroubi (2008)][2], and the
underlying physics is based on [Fukugita & Kawasaki (1994)][3].

[1]: http://adsabs.harvard.edu/abs/2018NewA...64....9K
[2]: http://adsabs.harvard.edu/abs/2008MNRAS.384.1080T
[3]: http://adsabs.harvard.edu/abs/1994MNRAS.269..563F

## What it does

Given a tabulated spectral energy distribution (SED) for a source embedded in a dark matter
halo, STARDUST computes the resulting ionization state (HI/HII, HeI/HeII/HeIII) and gas
temperature along a radial grid, at user-defined time intervals. It is useful for studying the
Epoch of Reionization, the extragalactic UV background, and the ionizing influence of individual
sources (quasars, Population III stars, star-forming galaxies) on their surroundings.

## Getting started

Clone the repository:

```bash
git clone https://github.com/eor/STARDUST.git STARDUST
```

Then:

1. **[Set up your environment](docs/SETUP.md)** — installing GSL, libconfig, and Boost on macOS,
   on an HPC cluster with environment modules, on an HPC cluster without them, or via an
   Apptainer/Singularity container.
2. **Build**: `cd src && cp Makefile.template Makefile` (edit the paths once), then `make`.
3. **[Run an example](docs/EXAMPLES.md)** — three ready-made test cases with sample SEDs and
   configs live in `tests/`.

## Repository layout

```
src/        C/C++ source for the STARDUST binary
scripts/    SED generation, plotting, and catalog-building helper scripts (Python)
tests/      Example configs, SEDs, and reference output for three test cases
docs/       Setup and usage documentation (see ROADMAP.md for development status)
plots/      Example output plots referenced in docs/EXAMPLES.md
```

## Configuration

A simulation run needs a SED file and a config file (see `docs/EXAMPLES.md` for examples).
Any parameter not set in the config file falls back to a default defined in
`src/config_defaults.h`; changing those defaults requires recompiling.

## Python helper scripts

The scripts under `scripts/` (SED generation, plotting, catalog building) need `numpy`, `scipy`,
and `matplotlib`. The STARDUST binary itself has no Python dependency. A conda environment file is
provided:

```bash
conda env create -f environment.yml
conda activate stardust
```

## References

* [Krause, F. et al. 2018, New Astronomy 64, 9K](http://adsabs.harvard.edu/abs/2018NewA...64....9K)
* [Thomas, R.M. & Zaroubi, S. 2008, MNRAS 384, 1080](http://adsabs.harvard.edu/abs/2008MNRAS.384.1080T)
* [Fukugita, M. & Kawasaki, M. 1994, MNRAS 269, 563](http://adsabs.harvard.edu/abs/1994MNRAS.269..563F)

## Copyright and licence

STARDUST was created by Rajat Thomas and Fabian Krause.

© 2018 - 2026 The STARDUST Authors

This programme is free software: you can redistribute it and/or modify it under the terms of the
GNU General Public License as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This programme is distributed in the hope that it will be useful, but **without any warranty**;
without even the implied warranty of **merchantability** or **fitness for a particular purpose**.
See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this programme.
If not, see http://www.gnu.org/licenses/.
