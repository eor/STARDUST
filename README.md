# STARDUST
Spectrum Transport Around RaDiating Universal Source Type

A cosmological 1D radiative transfer code, which is primarily based on [1], [2],[3], and [n]. An earlier version was first indrocuced in [n+1]

#### Rationale:

Cosmological radiative transfer can be an important tool when studying the thermal or ionization states of the 
intergalactic medium. This includes, for example, the study of the Epoch of Reionization or the computation of the 
extragalactic ionizing UV background spectrum on large scales, while on smaller scales it is important for e.g. 
studying the ionizing radiation emanating from a sources embedded in a dark matter halos, such as quasars or population 3 stars. 
Given a tabulated spectral energy distribution (SED) of a source associated with a dark matter host halo, STARDUST computes the ionization states as well as temperature along a radial grid and produces outputs at user-defined time intervals. 

---

## Setup

#### Software dependencies
STARDUST:
* C/C++ compiler
* [GNU Scientific library (GSL)](https://www.gnu.org/software/gsl/) 
* ODEint, which is part of [Boost](http://www.boost.org/)
* [libconfig](https://github.com/hyperrealm/libconfig)

Python scripts:
* numpy
* scipy
* matplotlib

#### Installation
First, install the dependencies, either by hand or via your package manager. We provide an installation script for 
STARDUST's dependencies in `scripts/install_dependencies.py`

Change the following variables in the Makefile (in the `src` directory) to your liking:

```bash
GSL_INCL = -I${HOME}/local/include
GSL_LIB  = -L${HOME}/local/bin

CONF_INCL= -I${HOME}/local/include
CONF_LIB = -L${HOME}/local/lib

BOOST_INCL = -I${HOME}/local/boost-current
```
Compile STARDUST by runing 
```bash
make
```
If you feel adventurous and want to change any of the setting in the header files or the want to play with different solvers (see, `Makefile`), make sure to recompile everything
```bash
make clean && make
```


---

## Generating Spectral Energy Distribution files

See `scripts/SED_generator/sed.py` to generate a few simple SEDs


---

## Example runs 


To run STARDUST you will need a SED file and a configuration file. Sample configuration files can be found with the source code or in the *tests* directory. A simulation can be started like this
```bash
user@machine:~/SD_project/$ ./STARDUST my_SD_config 
```



### A power-law-like SED
...

### A star-like SED 
...

### Stroemgren sphere test
...

### Create a profile catalog
...



## References

* [1] Fukugita, M & Kawasaki, M. 1994, MNRAS 269, 563
* [2]
* [3]
* [n]
* [n+1]

