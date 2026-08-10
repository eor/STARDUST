# Setting up STARDUST

STARDUST needs three external libraries:

| Library                                              | Used for                                  |
|-------------------------------------------------------|--------------------------------------------|
| [GSL](https://www.gnu.org/software/gsl/)              | interpolation, integration, ODE root-finding |
| [libconfig](https://github.com/hyperrealm/libconfig)  | reading the run's `.cfg` config file       |
| [Boost](https://www.boost.org/) (header-only, `odeint`)| the ODE solvers used in `rt.c`            |

You also need a C/C++ compiler (the code compiles as C++ via `g++`, despite the `.c` file
extensions — see [the note on this in the Makefile](#a-note-on-the-makefile)) and `make`. The
STARDUST binary itself has no Python dependency. Python is only needed for the optional
SED-generation and plotting scripts in `scripts/`, which require `numpy`, `scipy`, and
`matplotlib` (any reasonably recent Python works; 3.10+ is recommended so you get a current
scientific stack). A conda `environment.yml` for these is provided at the repository root:

```bash
conda env create -f environment.yml
conda activate stardust
```

Pick the scenario below that matches your situation, then jump to
[Pointing the Makefile at your libraries](#pointing-the-makefile-at-your-libraries).

- [A. macOS with Homebrew](#a-macos-with-homebrew)
- [B. HPC cluster with environment modules (Lmod)](#b-hpc-cluster-with-environment-modules-lmod)
- [C. HPC cluster without modules / no root access](#c-hpc-cluster-without-modules--no-root-access)
- [D. HPC cluster with Apptainer/Singularity](#d-hpc-cluster-with-apptainersingularity)
- [Submitting STARDUST runs via Slurm](#submitting-stardust-runs-via-slurm)

---

## A. macOS with Homebrew

All three dependencies are standard Homebrew formulae — no source builds needed.

```bash
brew install gsl libconfig boost
```

Homebrew installs into `/opt/homebrew` on Apple Silicon or `/usr/local` on Intel Macs. Find the
exact prefix with:

```bash
brew --prefix gsl
brew --prefix libconfig
brew --prefix boost
```

Use those paths in the Makefile (see below).

---

## B. HPC cluster with environment modules (Lmod)

Most clusters built with [EasyBuild](https://easybuild.io/) or similar ship GSL and Boost as
modules. libconfig is less commonly provided this way — check first, and fall back to scenario C
for it alone if it's missing.

```bash
module avail gsl
module avail boost
module avail libconfig
```

Load whichever versions are available, e.g.:

```bash
module load GSL/2.8
module load Boost/1.91.0
```

Module systems that follow the EasyBuild convention export environment variables you can feed
straight into the Makefile, typically named `$EBROOTGSL`, `$EBROOTBOOST`, `$EBROOTLIBCONFIG`
(run `module show GSL/2.8` — or whatever module you loaded — to see what it actually sets; naming
varies by site). If your site uses a different convention, `module show <name>` will tell you
which variables and paths it exports.

If libconfig isn't available as a module, build just that one library by hand following
scenario C below, and combine the module-provided GSL/Boost paths with your own libconfig prefix
in the Makefile.

---

## C. HPC cluster without modules / no root access

Use the provided script, which checks for each library and builds only what's missing, into
`~/local` by default:

```bash
python3 scripts/install_dependencies.py
```

It will print, at the end, the exact `Makefile` lines to use and the environment variables
(`LD_LIBRARY_PATH`, `PKG_CONFIG_PATH`) to add to your shell profile.

If you'd rather do it by hand, or the script doesn't work on your system, here's what it does
under the hood — using current official sources (the old `astro.rug.nl` mirror this script used
to pull from is no longer reachable):

```bash
INSTALL_DIR=$HOME/local
mkdir -p $INSTALL_DIR/src && cd $INSTALL_DIR/src

# GSL
wget https://ftp.gnu.org/gnu/gsl/gsl-2.8.tar.gz
tar -xzf gsl-2.8.tar.gz && cd gsl-2.8
./configure --prefix=$INSTALL_DIR && make -j4 && make install
cd ..

# libconfig
wget [https://github.com/hyperrealm/libconfig/releases/download/v1.8.2/libconfig-1.8.2.tar.gz](https://github.com/hyperrealm/libconfig/archive/refs/tags/v1.8.2.tar.gz)
tar -xzf v1.8.2.tar.gz && cd libconfig-1.8.2
./configure --prefix=$INSTALL_DIR && make -j4 && make install
cd ..

# Boost (header-only for our purposes — odeint needs no compiled libraries)
wget https://archives.boost.io/release/1.91.0/source/boost_1_91_0.tar.gz
tar -xzf boost_1_91_0.tar.gz
# nothing further to build — point BOOST_INCL at boost_1_91_0/ directly
```

Then add to your shell profile:

```bash
export LD_LIBRARY_PATH=$HOME/local/lib:$LD_LIBRARY_PATH
export PKG_CONFIG_PATH=$HOME/local/lib/pkgconfig:$PKG_CONFIG_PATH
```

---

## D. HPC cluster with Apptainer/Singularity

If your cluster runs jobs via [Apptainer](https://apptainer.org/) (formerly Singularity)
containers, this is usually the least fiddly option: you build the image once on any machine
with internet access and root/fakeroot (your laptop, a build node, CI), and the resulting `.sif`
file is completely self-contained — nothing needs to be downloaded, compiled, or modules-loaded
on the actual compute node.

A definition file is provided at `stardust.def`. It uses Ubuntu 24.04 LTS as the base image and
installs GSL, libconfig, and Boost from Ubuntu's own package repositories (all three are
standard `apt` packages on Ubuntu, so no source builds are needed inside the container either).
The STARDUST source tree is copied in and built during the image build.

Build the image (from the repository root, since the definition file's `%files` section uses
paths relative to it):

```bash
apptainer build stardust.sif stardust.def
```

This typically needs to be done on a machine where you have root or can use `--fakeroot`; check
with your cluster's documentation, since policies on this vary by site. Many clusters allow
`apptainer build --fakeroot` without full root access.

Transfer `stardust.sif` to the cluster (it's a single file — `scp`, `rsync`, or your usual
mechanism), then run it:

```bash
apptainer run stardust.sif tests/test_1_PL_source/config_PL_M11.000_z7.000.dat
```

or, to get a shell inside the container with the binary and Python scripts available:

```bash
apptainer shell stardust.sif
```

If you'd rather use 26.04 LTS or a different base, edit the `From:` line in `stardust.def` and
rebuild — the rest of the definition file doesn't need to change, since the same `apt` package
names (`libgsl-dev`, `libconfig-dev`, `libboost-dev`) are available across recent Ubuntu releases.

---

## Submitting STARDUST runs via Slurm

Everything above gets you a working `STARDUST` binary (or `.sif` image) — but on most HPC
clusters you can't just run it directly on the login node. You submit it as a job through a
scheduler, almost always [Slurm](https://slurm.schedmd.com/), and let it run on a compute node.
This applies whether you built STARDUST natively (scenarios B/C) or via Apptainer (scenario D).

A minimal submission script (adjust `--time`, `--mem`, and `--partition`/`--account` for your
cluster — these vary by site and aren't guessable from here):

```bash
#!/bin/bash
#SBATCH --job-name=stardust
#SBATCH --output=stardust_%j.out
#SBATCH --error=stardust_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=01:00:00

# If your cluster uses Lmod and STARDUST was built natively (scenario B):
module load GSL/2.8 Boost/1.91.0   # match whatever you loaded when building

srun ./STARDUST my_config.cfg
```

STARDUST is single-threaded (1D radial grid, no MPI/OpenMP), so `--ntasks=1`,
`--cpus-per-task=1` is normally all you need — this isn't a job that benefits from more cores
or nodes. Memory and time requirements scale with `settingsRMax`/`settingsDeltaR` (grid point
count) and `settingsDeltaT` (number of timesteps) in your config file; the defaults are modest,
but a long-`sourceLifetime` or fine-`settingsDeltaR` run can take a while; if unsure, request more
time than you think you need for a first run.

If you're running the Apptainer image instead (scenario D), swap the last line for:

```bash
# If apptainer isn't preloaded on your cluster's compute nodes, you may need:
# module load apptainer

srun apptainer run stardust.sif my_config.cfg
```

Submit with `sbatch your_script.sh`, then check status with `squeue -u $USER`. If you're new to
Slurm, your cluster's own documentation will have the authoritative answer on partition names,
account/allocation flags, and per-site memory/time limits — those vary enough between sites that
nothing generic here can substitute for it.

---

## Pointing the Makefile at your libraries

`src/Makefile.template` is the version tracked in git. Copy it to `Makefile` once, then edit
your copy — this way `git pull` never conflicts with your local path changes, since the file
you edit isn't tracked:

```bash
cd src
cp Makefile.template Makefile
```

Edit `Makefile` (your copy, not the template) and set:

```make
GSL_INCL   = -I/path/to/gsl/include
GSL_LIB    = -L/path/to/gsl/lib

CONF_INCL  = -I/path/to/libconfig/include
CONF_LIB   = -L/path/to/libconfig/lib

BOOST_INCL = -I/path/to/boost_1_91_0
```

Then build:

```bash
make
```

Editing any header (e.g. `config_defaults.h`, `constants.h`) or the `Makefile` itself — for
example switching the ODE solver selection at the top of the Makefile — triggers a rebuild of all
objects automatically (every object depends on the full header list via `INCL`), so a plain
`make` is enough. Use `make clean && make` only if you suspect stale build artifacts:

```bash
make clean && make
```

Unused-variable warnings during compilation are expected and harmless.

If `Makefile.template` itself is ever updated upstream (e.g. you pull a newer version of this
repo), your local `Makefile` won't pick up those changes automatically — diff the two files by
hand if you want to check what changed and merge anything relevant into your copy.

### A note on the Makefile

`rt.c` uses Boost.odeint, which is C++ rather than C, so the whole codebase is compiled with
`g++` (via the `.c.o` rule) rather than `gcc`, despite most files using the `.c` extension. This
is intentional, if a little surprising on first read.

### A note on GSL versions

The redshift root-finder in `rt.c` uses the older `gsl_odeiv` interface (not `gsl_odeiv2`).
This interface was deprecated starting with GSL 1.15 but is still shipped for backward
compatibility in current GSL releases (2.8 as of writing), so it continues to work fine — no
version constraint here. A migration to `gsl_odeiv2` may happen as part of a future `rt.c`
refactor, but isn't required for current GSL to work.
