# recurcast

Private development repository for recurcast 0.10.1. The R package is in `package/`.
It provides coherent forecasts for recurrent and terminal events, with explicit
observation assumptions and numerical error bounds.

## Build and check

```sh
R CMD build package
R CMD check --as-cran recurcast_0.10.1.tar.gz
```

Install the dependencies listed in `package/DESCRIPTION` first. Native code needs
a C++17 compiler. See `package/README.md`, `package/NEWS.md`, and installed help.

GitHub Actions checks R release/devel on Linux, macOS and Windows, and standalone
native fixtures with GCC/Clang and address/undefined-behaviour sanitizers.
Workflow status is evidence only after a completed run.

The 0.10.1 additions preserve the 0.10.0 numerical kernels. Standard survival and
competing-risk wrappers use piecewise-exponential proportional hazards without
frailty. Audit exports support reproducibility; they are not regulatory certification.

Package license: GPL-3; see `package/LICENSE.md`. This private CI repository does
not include dissertation manuscripts, private clinical data, or submission archives.
