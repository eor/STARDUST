
# STARDUST
STARDUST, which stands for **S**pectrum **T**ransport **A**round **R**a**D**iating **U**niversal **S**ource **T**ype, 
is a cosmological 1D radiative transfer code that is described in detail by [Krause et al. (2018)][1]. 
An earlier (unpublished) version of the code was employed in a paper by [Thomas & Zaroubi (2008)][2] 
and it should be noted that the code is primarily based on the theoretical work presented in [Fukugita & Kawasaki (1994)][3].

[1]: http://adsabs.harvard.edu/abs/2018NewA...64....9K
[2]: http://adsabs.harvard.edu/abs/2008MNRAS.384.1080T
[3]: http://adsabs.harvard.edu/abs/1994MNRAS.269..563F



#### Rationale:

Cosmological radiative transfer can be an important tool when studying the thermal or ionization states of the 
intergalactic medium. This includes, for example, the study of the Epoch of Reionization or the computation of the 
extragalactic ionizing UV background spectrum on large scales, while on smaller scales it is important for e.g. 
studying the ionizing radiation emanating from a sources embedded in a dark matter halos, such as quasars or population 3 stars. 
Given a tabulated spectral energy distribution (SED) of a source associated with a dark matter host halo, 
STARDUST computes the ionization states as well as temperature along a radial grid and produces outputs at user-defined time intervals. 



## Setup

### Obtaining STARDUST

You can clone the repository with the following command
```bash
git clone https://github.com/eor/STARDUST.git STARDUST
```
Should your system lack `git`, you can alternatively download a zip file of the repository [here](https://github.com/eor/STARDUST/archive/master.zip).

### Software requirements

We assume your system is equiped with the following dependencies:
* C/C++ compiler
* make
* wget
* Python 3.5 or newer 
* numpy
* scipy
* matplotlib

Furthermore, STARDUST needs the following libraries:
* [GNU Scientific library (GSL)](https://www.gnu.org/software/gsl/) 
* ODEint, which is part of [Boost](http://www.boost.org/)
* [libconfig](https://github.com/hyperrealm/libconfig)

There are different ways to install these: 

* by hand 
* via your favorite package manager
* using the install script provided in `scripts/install_dependencies.py`

Once the dependencies are in place, change the following path variables in STARDUST's  `src/Makefile` to your liking:
```bash
GSL_INCL = -I/path/to/gsl/include_files
GSL_LIB  = -L/path/to/gsl/library_files

CONF_INCL= -I/path/to/libconfig/include_files
CONF_LIB = -L/path/to/libconfig/library_files

BOOST_INCL = -I/path/to/boost/include_files
```

If you choose to install the dependencies by hand, the `LD_LIBRARY_PATH` environment variable needs to be adjusted to include the path to the gsl and libconfig library directory.
We suggest to add the following to your `~/.bashrc`
```bash
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/path/to/libconfig/library_files
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/path/to/gsl/library_files

```

#### Notes on the install script

If you choose to use the install script, open the script in your text editor, change the `installDir` variable and then execute the script.
Once the script finishes, it will print the path variables that need to be inserted into STARDUST's `src/Makefile`.

### Compilation

Compile STARDUST by changing into the `src` directory and running 
```bash
make
```
Don't mind the unused variable warnings. If you feel adventurous and want to change any of the setting in the header files or 
want to play with different solvers (see `src/Makefile`), make sure to recompile everything in the `src` directory
```bash
make clean && make
```


## Running simulations

To run STARDUST you will need a SED file for the source *and* a configuration file which defines some parameters for the run. Sample SED and configuration files can be found in the `tests` directory and are explained in more detail below. 

Once you have both files, a simulation can be started like this
```bash
user@machine:~/SD_project/$ ./STARDUST my_SD_config 
```

If a configuration variable is not set in the configuration file, STARDUST will fall back on a default value for said variable. The defaults can be changed in `src/config_defaults.h`. Note that any changes made therein requires a re-compilation of the code.



### Generating simple Spectral Energy Distributions

We provide a script to generate SEDs for some toy source models, which are described in more detail in our BEARS pipeline paper, see [Krause et al. (2018)][1].
The script, which can also be used as a python package, can be found in `scripts/sed_generator/sed.py`. 


## Example runs 

A number of test cases can be found in the `tests` directory. Each comes with a SED and configuration 
file that can be used to ensure the installation is working. Copy or symlink your STARDUST binary to 
the respective directories and start the runs.

### Test 1: A power-law-like SED
Here, the SED was generated using the following code.

```python
import sed

z    = 7.
logM = 11.0
f1   = 'sed_PL_M%.3f_z%.3f.dat'%(logM, z)

sed.generate_SED_PL( haloMass=10**logM, 
                     eHigh=1.e4, eLow=13.6, 
                     N=1000, logGrid=True,
                     alpha=1.0,  qsoEfficiency=0.1, 
                     fileName=f1 )

```


### Test 2: A star-like SED 
Now we use a SED that uses our simple pop 3 model (which follows an IMF). It was generated as follows.

```python
import sed

z    = 7.
logM = 11.0
f2   = 'sed_IMF_M%.3f_z%.3f.dat'%(logM, z)

sed.generate_SED_stars_IMF( haloMass=10**logM, redshift=z, 
                            eLow=13.6, eHigh=1.e4, 
                            N=1000,  logGrid=True, 
                            starMassMin=5, starMassMax=100, imfBins=99, 
                            imfIndex=2.35, fEsc=0.1,
                            targetSourceAge=10.0, 
                            fileName=f2 )

```

### Test 3: Co-evolution SED
For this test we emloy a model that combines SEDs from the first two tests.

```python
import sed

z    = 7.
logM = 11.0
f3   = 'sed_IMF+PL_M%.3f_z%.3f.dat'%( logM, z )

sed.generate_SED_IMF_PL( haloMass=10**logM, redshift=z, 
                         eLow=13.6, eHigh=1.e4, 
                         N=1000,  logGrid=True, 
                         starMassMin=30, starMassMax=100, imfBins=99, 
                         imfIndex=2.35,  fEsc=0.1,
                         targetSourceAge=10.0,
                         alpha=1.0, qsoEfficiency=0.1, 
                         fileName=f3 )                            

```


### Visualizing simulation results
Example scripts for plotting the resulting STARDUST profiles can be found in the `scripts/sd_plot` directory, 
see `sd_plot_time_evolution.py` or  `sd_plot_model_comparison.py` to get started. Note that we also included the two 
plots produced with these scripts in the `plots` directory. 


### Creating a profile catalog
The `scripts/catalog` directory contains a simple script (`make_sd_catalog.py`) to generate a catalog of STARDUST profiles (e.g. to be used with BEARS). 
It contains examples on how to create said catalogs with our toy SED generator, see the script (comments) for further details.



## References
* [Krause, F. et al. 2018, NewAstronomy 64, 9K](http://adsabs.harvard.edu/abs/2018NewA...64....9K)
* [Thomas, R.M. & Zaroubi, S. 2008, MNRAS 384, 1080](http://adsabs.harvard.edu/abs/2008MNRAS.384.1080T)
* [Fukugita, M. & Kawasaki, M. 1994, MNRAS 269, 563](http://adsabs.harvard.edu/abs/1994MNRAS.269..563F)


## Copyright and licence

STARDUST was created by the following people: Rajat Thomas and Fabian Krause

© 2018-2021 The STARDUST Authors

This programme is free software: you can redistribute it and/or modify it 
under the terms of the GNU General Public License as published by the Free 
Software Foundation, either version 3 of the License, or (at your option) any 
later version.

This programme is distributed in the hope that it will be useful, but **without 
any warranty**; without even the implied warranty of **merchantability** or **fitness 
for a particular purpose**. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with 
this programme. If not, see http://www.gnu.org/licenses/.

