Python access to coherentSurv
============================

Install the R package and its dependencies first. Both examples call the same
R/C++ implementation; there is no second statistical implementation to maintain.

Linux rpy2 route (R and a compatible rpy2 build required):
    python run_rpy2.py --output forecasts.csv --timings
Optional transport verification (performs a second, direct-R fit):
    python run_rpy2.py --output forecasts.csv --verify-parity --timings

Process-based route, also usable on Windows with Rscript on PATH:
    python run_forecast.py --help

Use --help on either script before choosing your input file. These are runnable
examples with a declared schema and fixed model specification, not universal
estimators for arbitrary tables. Adapt the formula, event marks, time unit,
landmark and horizons to your scientific question. Identifiers are text and
preserve leading zeros and literal NA labels. Empty numeric fields are missing,
not zero. Empty event labels mean no observed event in that exact interval.

The rpy2 runner retains the fitted R object. By default it fits once and converts
only input data and the requested output table. The parity sidecar explicitly
says not_requested unless --verify-parity was used. Stage timings exclude R
initialization and Python imports; they are measurements of this run only.
Conversions may copy data. No zero-copy, thread-safety or speed superiority is
claimed. Use separate processes for independent concurrent R sessions.

The rpy2 integration is tested on Linux. Windows users can use Rscript; a Windows
embedded-rpy2 result is not inferred from Linux tests. Prerequisites:
https://rpy2.github.io/doc/latest/html/overview.html
Conversion documentation:
https://rpy2.github.io/doc/latest/html/pandas.html
