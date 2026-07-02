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
    subprocess.run(cmd, shell=True, check=True, cwd=cwd)


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
    if os.path.exists(dest):
        print(f"  already downloaded: {dest}")
        return dest
    print(f"  downloading {url}")
    if have_tool("wget"):
        run(f"wget -c '{url}' -O '{dest}'")
    else:
        urllib.request.urlretrieve(url, dest)
    return dest


# ----------------------------------------------------------------------
# Per-dependency detection
# ----------------------------------------------------------------------

def gsl_present(prefix):
    """True if a usable GSL is already on this system (system-wide or
    in our own prefix)."""
    gsl_config = shutil.which("gsl-config") or os.path.join(prefix, "bin", "gsl-config")
    if os.path.exists(gsl_config) or shutil.which("gsl-config"):
        try:
            out = subprocess.run([gsl_config if os.path.exists(gsl_config) else "gsl-config", "--version"],
                                 capture_output=True, text=True, check=True)
            print(f"  found GSL {out.stdout.strip()} via gsl-config")
            return True
        except Exception:
            pass
    return False


def libconfig_present(prefix):
    candidates = [
        os.path.join(prefix, "include", "libconfig.h"),
        "/usr/include/libconfig.h",
        "/usr/local/include/libconfig.h",
    ]
    for c in candidates:
        if os.path.exists(c):
            print(f"  found libconfig header at {c}")
            return True
    return False


def boost_odeint_present(prefix):
    candidates = [
        os.path.join(prefix, "boost_" + BOOST_VERSION_US, "boost", "numeric", "odeint.hpp"),
        "/usr/include/boost/numeric/odeint.hpp",
        "/usr/local/include/boost/numeric/odeint.hpp",
    ]
    for c in candidates:
        if os.path.exists(c):
            print(f"  found Boost.odeint header at {c}")
            return True
    return False


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

    built = {"gsl": False, "libconfig": False, "boost": False}

    print("\nChecking GSL ...")
    if force or not gsl_present(prefix):
        build_gsl(prefix, src_dir)
        built["gsl"] = True
    else:
        print("  GSL already present, skipping (use --force to rebuild)")

    print("\nChecking libconfig ...")
    if force or not libconfig_present(prefix):
        build_libconfig(prefix, src_dir)
        built["libconfig"] = True
    else:
        print("  libconfig already present, skipping (use --force to rebuild)")

    print("\nChecking Boost (odeint headers) ...")
    if force or not boost_odeint_present(prefix):
        fetch_boost(prefix, src_dir)
        built["boost"] = True
    else:
        print("  Boost.odeint already present, skipping (use --force to rebuild)")

    print("\nAll done.\n")
    print("Add this to your shell profile (~/.bashrc or similar):\n")
    print(f"    export LD_LIBRARY_PATH={prefix}/lib:$LD_LIBRARY_PATH")
    print(f"    export PKG_CONFIG_PATH={prefix}/lib/pkgconfig:$PKG_CONFIG_PATH\n")
    print("If you haven't already, copy the Makefile template once:\n")
    print("    cd src && cp Makefile.template Makefile\n")
    print("Then set these paths in src/Makefile (your copy, not the template):\n")
    print(f"    GSL_INCL   = -I{prefix}/include")
    print(f"    GSL_LIB    = -L{prefix}/lib")
    print(f"    CONF_INCL  = -I{prefix}/include")
    print(f"    CONF_LIB   = -L{prefix}/lib")
    print(f"    BOOST_INCL = -I{prefix}/boost_{BOOST_VERSION_US}\n")
    print("Then: cd src && make\n")


if __name__ == "__main__":
    main()