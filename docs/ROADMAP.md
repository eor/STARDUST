# STARDUST TODO

## Open - needs physics review before applying

- **`table_ion.c` HeII photoionization cross-section.** `tau_ehe2_p3_log` (and the dead
  `tau_ehe2_p4_log`) weight the integrand with `k * F_ehe1` (HeI), but by the Fukugita & Kawasaki
  (1994) scaling used elsewhere it should be `k * F_e1h1` (HI) — `table_temp.c`'s analogous term
  already uses `k * F_e1h1`. This **changes numerical output**, so it is intentionally left
  unapplied pending review against the paper and with the original author. A standalone impact
  test comparing both formulas against `cross_sec_ehe2()` should accompany the fix.



## Open — lower priority

- **Distribute `functions.c`** into `cosmology.c` (`hubble`, `Findz`, `time_from_redshift`) and
  populate the empty `physics.c` stub (`cross_sec_e1h1/ehe1/ehe2`). When `physics.c` gains content,
  **add `physics.o` to `OBJS`** in `Makefile` and `Makefile.template`.
- Refactor the `Findz` redshift root-finder.
- Migrate the deprecated `gsl_odeiv` API in `rt.c` to `gsl_odeiv2`.
- Only write/read the integral tables in debug mode (skip the disk round-trip otherwise).
- Prepare for a choice of UVB spectrum / file input.
- Further `ode_derivatives.h` / `ode_jacobian.h` cleanup.
- Naming: the integrals/tables variables and the physical constants in `constants.h` are still
  terse (constants deliberately left alone for now — high churn, low value, bug risk).
- Audit logging (close on early exit) across the not-yet-refactored files (`rt.c`, `table_*.c`).
- One-off trailing-whitespace cleanup as an isolated commit (`.editorconfig` now keeps new edits
  clean).




## Done

### Build & tooling

- Fixed the build blocker: `log.h`'s unused `check` macro clashed with Boost 1.91 in `rt.c`
  (removed it).
- `Makefile` → tracked `Makefile.template` (users copy once, edit local paths); `-std=gnu++17`
  pinned; retired the `-DSTROEMGRENTEST` compile flag.
- `INCL` now lists every header, so editing any header triggers a rebuild automatically.
- New: `.gitignore`, `.editorconfig`, `environment.yml` (conda), `stardust.def` (Apptainer),
  `docs/SETUP.md`, `docs/EXAMPLES.md`, rewritten `scripts/install_dependencies.py`, new `README`.

### Source cleanup & bug fixes

- `sed.c`: rewritten — fixed the `while(!feof)` over-count, removed dead code and the legacy
  unit-conversion / `type != 2` branch, single-pass read, and **cached the interpolation splines**
  (built once, not per call).
- `sed.c`/`allvars.*`/`memory.c`: renamed the misleading globals `Lambda`/`Energy` →
  `photonEnergy` [eV] / `sedLuminosity` [eV/s/eV].
- `config.c`: abort if `paths.pathSED` is missing (no meaningful default); fixed the
  `sourceEHigh` missing-key case that left it uninitialised.
- **Strömgren-sphere test is now a runtime option** (`settings.settingsStroemgrenTest`, plus
  configurable `stroemgrenPeakE`/`stroemgrenWidth`) instead of the `STROEMGRENTEST` compile flag —
  converted across `rt.c`, `sed.c`, `constants.h`, `ode_derivatives.h`, `ode_jacobian.h`.
- `density.c`: fixed the `while(!feof)` over-count that made it reject correctly-sized density
  profiles; added logging on error paths and the missing `fclose`.
- `utils.c`: `const`-correct path helpers; correct handling of an empty `pathID`.
- `interpolation.c`: per-call GSL spline allocation → cached objects reused across calls
  (the sim-loop hot path).
- Removed the dead `type` global.

### Python scripts

- Fixed the SED generator for modern SciPy: `integrate.simps(...)` (removed in SciPy 1.14) →
  `integrate.simpson(...)` in `sed.py`. Verified end-to-end in the conda env: generator → SED
  file → refactored STARDUST reads/normalizes it correctly.
- Ported the driver scripts to Python 3: converted Python-2 `print` statements (`make_sd_catalog.py`,
  `sd_plot_time_evolution.py`, `sd_plot_model_comparison.py` — they were hard `SyntaxError`s and
  didn't run) and made the LaTeX label strings raw (`r'...'`) to clear the `\,`/`\ ` escape
  warnings. All scripts now compile clean under Python 3; verified the time-evolution plotter
  renders a figure from real STARDUST profiles. The catalog generator reuses the fixed generator
  and correctly sets `pathSED` (runnable once the user points `sdExePath` at their binary).

