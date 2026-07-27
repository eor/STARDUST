#!/usr/bin/env python3
"""
Checks for STARDUST's three dependencies (GSL, libconfig, Boost) and
builds whichever ones are missing from official upstream sources into
a local prefix (no root required).

Usage:
    python3 install_dependencies.py [--prefix DIR] [--force]

    --prefix DIR   Install location (default: ~/local)
    --force        Rebuild a dependency even if it looks like it's
                   already present

This is intended for HPC clusters without environment modules and
without root access (see docs/SETUP.md, scenario C). If you have
Homebrew (macOS) or Lmod-style modules (most HPC clusters) available,
use those instead -- see docs/SETUP.md scenarios A and B.

Note: the old version of this script downloaded a custom "reduced"
Boost tarball and GSL/libconfig copies from a personal university
mirror (astro.rug.nl), which is no longer reachable. This version
pulls from each project's own current release instead.
"""
import os
import shutil
import subprocess
import sys
import urllib.request
import urllib.error

# ----------------------------------------------------------------------
# Versions / sources -- bump these as needed
# ----------------------------------------------------------------------

GSL_VERSION = "2.8"
GSL_URL = f"https://ftp.gnu.org/gnu/gsl/gsl-{GSL_VERSION}.tar.gz"

LIBCONFIG_VERSION = "1.8.2"
LIBCONFIG_URL = (
    "https://github.com/hyperrealm/libconfig/releases/download/"
    f"v{LIBCONFIG_VERSION}/libconfig-{LIBCONFIG_VERSION}.tar.gz"
)

# Boost is header-only for our purposes (Boost.odeint needs no compiled
# libraries), so we just need the source tree extracted somewhere, not
# a configure/make/install cycle.
BOOST_VERSION = "1.91.0"
BOOST_VERSION_US = BOOST_VERSION.replace(".", "_")
BOOST_URL = (
    f"https://archives.boost.io/release/{BOOST_VERSION}/source/"
    f"boost_{BOOST_VERSION_US}.tar.gz"
)


def run(cmd, cwd=None):
    print(f"  > {cmd}")
    try:
        subprocess.run(cmd, shell=True, check=True, cwd=cwd)
    except subprocess.CalledProcessError as e:
        sys.exit(f"\nERROR: command failed (exit {e.returncode}):\n  {cmd}\n"
                 f"See the output above for details.")


def have_tool(name):
    return shutil.which(name) is not None


def check_prereqs():
    print("Checking for build tools (gcc/g++, make, wget or curl) ...")
    missing = [t for t in ("make",) if not have_tool(t)]
    if not (have_tool("gcc") or have_tool("g++")):
        missing.append("gcc/g++")
    if not (have_tool("wget") or have_tool("curl")):
        missing.append("wget or curl")
    if missing:
        sys.exit(f"ERROR: missing required tool(s): {', '.join(missing)}. "
                 f"Install these first (or ask your cluster admin for a "
                 f"module that provides them) and re-run.")
    print("  ok")


def download(url, dest_dir):
    fname = url.rsplit("/", 1)[-1]
    dest = os.path.join(dest_dir, fname)
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        print(f"  already downloaded: {dest}")
        return dest
    print(f"  downloading {url}")
    if have_tool("wget"):
        run(f"wget -c '{url}' -O '{dest}'")           # -c resumes a partial download
    elif have_tool("curl"):
        run(f"curl -fL '{url}' -o '{dest}'")          # -f: fail on HTTP error, -L: follow redirects
    else:
        try:
            urllib.request.urlretrieve(url, dest)
        except (urllib.error.URLError, OSError) as e:
            if os.path.exists(dest):
                os.remove(dest)                       # don't leave a partial file for the next run to trust
            sys.exit(f"\nERROR: failed to download {url}\n  ({e})\n"
                     f"Check your network connection (and any proxy settings) and try again.")
    return dest


# ----------------------------------------------------------------------
# Per-dependency detection
#
# Each returns *where* the dependency lives so the same value can drive both
# the build/skip decision and the final Makefile hint:
#   - GSL / libconfig: (include_dir, lib_dir)
#   - Boost:           include_dir   (header-only, no lib_dir)
# or None if not found.
# ----------------------------------------------------------------------

# common install roots to probe, in priority order (our own prefix first)
def _bases(prefix):
    return [prefix, "/usr", "/usr/local", "/opt/homebrew", "/opt/local"]


def find_gsl(prefix):
    """gsl-config is authoritative for GSL's own paths."""
    gsl_config = shutil.which("gsl-config") or os.path.join(prefix, "bin", "gsl-config")
    if not os.path.exists(gsl_config):
        return None
    try:
        ver = subprocess.run([gsl_config, "--version"], capture_output=True, text=True, check=True).stdout.strip()
        gpfx = subprocess.run([gsl_config, "--prefix"], capture_output=True, text=True, check=True).stdout.strip()
    except Exception:
        return None
    print(f"  found GSL {ver} at {gpfx}")
    return os.path.join(gpfx, "include"), os.path.join(gpfx, "lib")


def find_libconfig(prefix):
    for base in _bases(prefix):
        if os.path.exists(os.path.join(base, "include", "libconfig.h")):
            print(f"  found libconfig at {base}")
            return os.path.join(base, "include"), os.path.join(base, "lib")
    return None


def find_boost(prefix):
    # our own extracted tarball keeps boost/ directly under boost_x_y_z/;
    # system installs keep it under <base>/include/
    include_roots = [os.path.join(prefix, "boost_" + BOOST_VERSION_US)]
    include_roots += [os.path.join(base, "include") for base in _bases(prefix)]
    for inc in include_roots:
        if os.path.exists(os.path.join(inc, "boost", "numeric", "odeint.hpp")):
            print(f"  found Boost.odeint at {inc}")
            return inc
    return None


# ----------------------------------------------------------------------
# Per-dependency build
# ----------------------------------------------------------------------

def build_gsl(prefix, src_dir):
    print(f"Building GSL {GSL_VERSION} ...")
    tarball = download(GSL_URL, src_dir)
    run(f"tar -xzf '{tarball}' -C '{src_dir}'")
    build_dir = os.path.join(src_dir, f"gsl-{GSL_VERSION}")
    run(f"./configure --prefix='{prefix}'", cwd=build_dir)
    run("make -j4", cwd=build_dir)
    run("make install", cwd=build_dir)


def build_libconfig(prefix, src_dir):
    print(f"Building libconfig {LIBCONFIG_VERSION} ...")
    tarball = download(LIBCONFIG_URL, src_dir)
    run(f"tar -xzf '{tarball}' -C '{src_dir}'")
    build_dir = os.path.join(src_dir, f"libconfig-{LIBCONFIG_VERSION}")
    run(f"./configure --prefix='{prefix}'", cwd=build_dir)
    run("make -j4", cwd=build_dir)
    run("make install", cwd=build_dir)


def fetch_boost(prefix, src_dir):
    print(f"Fetching Boost {BOOST_VERSION} (header-only use, no build needed) ...")
    tarball = download(BOOST_URL, src_dir)
    target = os.path.join(prefix, f"boost_{BOOST_VERSION_US}")
    if os.path.isdir(target):
        print(f"  {target} already exists, skipping extraction")
        return
    run(f"tar -xzf '{tarball}' -C '{prefix}'")


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------

def main():
    prefix = os.path.expanduser("~/local")
    force = False
    args = sys.argv[1:]
    while args:
        a = args.pop(0)
        if a == "--prefix":
            prefix = os.path.expanduser(args.pop(0))
        elif a == "--force":
            force = True
        else:
            sys.exit(f"Unknown argument: {a}")

    src_dir = os.path.join(prefix, "src")
    os.makedirs(src_dir, exist_ok=True)
    os.makedirs(os.path.join(prefix, "lib"), exist_ok=True)
    os.makedirs(os.path.join(prefix, "include"), exist_ok=True)

    print(f"Install prefix: {prefix}\n")
    check_prereqs()

    # For each dependency, resolve where it lives: use what's already on the
    # system, or build/fetch it into `prefix` and use that. The resolved
    # location then drives the Makefile hint below -- no assumptions about
    # everything living under `prefix`.
    status = {}   # name -> "found" | "built"

    print("\nChecking GSL ...")
    gsl = find_gsl(prefix)
    if force or gsl is None:
        build_gsl(prefix, src_dir)
        gsl = (os.path.join(prefix, "include"), os.path.join(prefix, "lib"))
        status["gsl"] = "built"
    else:
        status["gsl"] = "found"

    print("\nChecking libconfig ...")
    libconfig = find_libconfig(prefix)
    if force or libconfig is None:
        build_libconfig(prefix, src_dir)
        libconfig = (os.path.join(prefix, "include"), os.path.join(prefix, "lib"))
        status["libconfig"] = "built"
    else:
        status["libconfig"] = "found"

    print("\nChecking Boost (odeint headers) ...")
    boost_inc = find_boost(prefix)
    if force or boost_inc is None:
        fetch_boost(prefix, src_dir)
        boost_inc = os.path.join(prefix, f"boost_{BOOST_VERSION_US}")
        status["boost"] = "built"
    else:
        status["boost"] = "found"

    gsl_inc, gsl_lib = gsl
    conf_inc, conf_lib = libconfig

    print("\nAll done.\n")
    print("Summary:")
    print(f"    GSL       : {status['gsl']:<5}  ({gsl_inc})")
    print(f"    libconfig : {status['libconfig']:<5}  ({conf_inc})")
    print(f"    Boost     : {status['boost']:<5}  ({boost_inc})\n")

    # LD_LIBRARY_PATH / PKG_CONFIG_PATH only matter for compiled libraries we
    # put in `prefix` (GSL, libconfig); Boost is header-only and anything
    # "found" is already on the system's default search paths.
    if status["gsl"] == "built" or status["libconfig"] == "built":
        print("You built libraries into your prefix. Add this to your shell profile")
        print("(~/.bashrc or similar) so the runtime linker can find them:\n")
        print(f"    export LD_LIBRARY_PATH={prefix}/lib:$LD_LIBRARY_PATH")
        print(f"    export PKG_CONFIG_PATH={prefix}/lib/pkgconfig:$PKG_CONFIG_PATH\n")

    print("If you haven't already, copy the Makefile template once:\n")
    print("    cd src && cp Makefile.template Makefile\n")
    print("Then set these paths in src/Makefile (your copy, not the template):\n")
    print(f"    GSL_INCL   = -I{gsl_inc}")
    print(f"    GSL_LIB    = -L{gsl_lib}")
    print(f"    CONF_INCL  = -I{conf_inc}")
    print(f"    CONF_LIB   = -L{conf_lib}")
    print(f"    BOOST_INCL = -I{boost_inc}\n")
    print("Then: cd src && make\n")


if __name__ == "__main__":
    main()