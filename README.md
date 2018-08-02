# STARDUST
STARDUST, which stands for *S*pectrum *T*ransport *A*round *R*a*D*iating *U*niversal *S*ource *T*ype, 
is a cosmological 1D radiative transfer code that is described in detail by [Krause et al. (2018)][1]. 
An earlier (unpublished) version of the code was first employed in a paper by [Thomas & Zaroubi (2008)][2], 
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

---

## Setup
We assume your system is equiped with the following dependencies:
* C/C++ compiler
* make
* wget
* Python 2.7 
* numpy
* scipy
* matplotlib

### Software dependencies
For STARDUST we furthermore need the following libraries:
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


---

## Generating Spectral Energy Distribution files

See `scripts/sed_generator/sed.py` to generate a few simple SEDs. The script can be used as a python package or as a stand alone script. 


---

## Example runs 


To run STARDUST you will need a SED file and a configuration file. Sample configuration files can be found with the source code or in the *tests* directory. A simulation can be started like this
```bash
user@machine:~/SD_project/$ ./STARDUST my_SD_config 
```

Note that if a configuration variable is not set, STARDUST will fall back on a default value for said variable. The defaults can be changed in `src/config_defaults.h`.


### A power-law-like SED
...

### A star-like SED 
...

### Stroemgren sphere test
...

### Create a profile catalog
...



## References
* [Krause, F. et al. 2018, NewAstronomy 64, 9K](http://adsabs.harvard.edu/abs/2018NewA...64....9K)
* [Thomas, R.M. & Zaroubi, S. 2008, MNRAS 384, 1080](http://adsabs.harvard.edu/abs/2008MNRAS.384.1080T)
* [Fukugita, M. & Kawasaki, M. 1994, MNRAS 269, 563](http://adsabs.harvard.edu/abs/1994MNRAS.269..563F)

