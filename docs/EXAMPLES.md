# Running STARDUST: examples and tests

Once you've [built the binary](SETUP.md), you need two files to start a run: a SED file for the
radiating source, and a config file with the run's parameters.

```bash
./STARDUST my_config.cfg
```

Three ready-made test cases live under `tests/`, each with its own SED and config file.

## Test 1: a power-law-like SED

```
tests/test_1_PL_source/
├── config_PL_M11.000_z7.000.dat
├── sed_PL_M11.000_z7.000.dat
└── test_1_M11.000_z7.000_log_main
```

The SED was generated with the `sed_generator` package (see below) like this:

```python
from sed import sed

z, logM = 7., 11.0
sed.generate_SED_PL(
    haloMass=10**logM, eHigh=1.e4, eLow=13.6, N=1000, logGrid=True,
    alpha=1.0, qsoEfficiency=0.1,
    fileName=f'sed_PL_M{logM:.3f}_z{z:.3f}.dat',
)
```

## Test 2: a star-like SED (population III, IMF-weighted)

```
tests/test_2_IMF_source/
├── config_IMF_M11.000_z7.000.dat
├── sed_IMF_M11.000_z7.000.dat
└── test_2_M11.000_z7.000_log_main
```

```python
from sed import sed

z, logM = 7., 11.0
sed.generate_SED_stars_IMF(
    haloMass=10**logM, redshift=z, eLow=13.6, eHigh=1.e4, N=1000, logGrid=True,
    starMassMin=5, starMassMax=100, imfBins=99, imfIndex=2.35, fEsc=0.1,
    targetSourceAge=10.0,
    fileName=f'sed_IMF_M{logM:.3f}_z{z:.3f}.dat',
)
```

## Test 3: a combined stellar + power-law SED

```
tests/test_3_IMF+PL_source/
├── config_IMF+PL_M11.000_z7.000.dat
└── sed_IMF+PL_M11.000_z7.000.dat
```

```python
from sed import sed

z, logM = 7., 11.0
sed.generate_SED_IMF_PL(
    haloMass=10**logM, redshift=z, eLow=13.6, eHigh=1.e4, N=1000, logGrid=True,
    starMassMin=30, starMassMax=100, imfBins=99, imfIndex=2.35, fEsc=0.1,
    targetSourceAge=10.0, alpha=1.0, qsoEfficiency=0.1,
    fileName=f'sed_IMF+PL_M{logM:.3f}_z{z:.3f}.dat',
)
```

## Running a test case

```bash
cp src/STARDUST tests/test_1_PL_source/
cd tests/test_1_PL_source
./STARDUST config_PL_M11.000_z7.000.dat
```

Output profiles and logs are written into the same directory (or wherever `pathOutDir` in the
config points).

## Generating your own SEDs

The `sed_generator` package (`scripts/sed_generator/sed/sed.py`) can be used standalone:

```bash
cd scripts/sed_generator
python3 example.py
```

See the comments in `sed/sed.py` and `sed/sed_stellar_mass.py` for the underlying stellar-mass
and IMF models, described in more detail in the BEARS pipeline paper
([Krause et al. 2018](http://adsabs.harvard.edu/abs/2018NewA...64....9K)).

## Visualizing results

Plotting scripts live in `scripts/sd_plot/`:

- `sd_plot_time_evolution.py` — evolution of a single profile over time
- `sd_plot_model_comparison.py` — compare profiles across different source models

Example outputs from these scripts are in `plots/`.

## Building a catalog of runs

`scripts/catalog/make_sd_catalog.py` generates a catalog of STARDUST profiles (e.g. for use with
BEARS), using `scripts/catalog/SD_CONFIG_TEMPLATE` as a config template. See the comments in that
script for details on sweeping over halo mass, redshift, or source model.