## Submission summary

This is coherent4cast 0.10.3, replacing the pending 0.10.2 submission.
All public fitting functions, native namespaces, header paths, help pages,
installation examples and distributed materials now use one package identity.
Transitional fitting aliases have been removed. Numerical formulas and
algorithms are unchanged from 0.10.2. Source is licensed under GPL-3.

Maintainer: Sunday A. Adetunji <sundayadetunjisa@gmail.com>
Source archive SHA-256: 7c095cb3483df9b8f19018f6a086b1dee8762b4f5e4551c698cd64a262efbe77

## Test environments

Final submitted archive: Windows 11 x64, R-devel 2026-11-17 r90926 ucrt,
GCC 14.3.0; built under R 4.6.2. R CMD check --as-cran included remote
incoming checks, all examples, five vignettes and PDF/HTML manuals.

All six Linux/macOS/Windows R-release/R-devel jobs passed:
https://github.com/setrojalab/coherent4cast/actions/runs/45570154240
All four GCC/Clang native builds with and without ASAN/UBSAN passed:
https://github.com/setrojalab/coherent4cast/actions/runs/45570154250
Tested source revision: e06edb1ccab4e816d7fe23e44b924c7d8d31976b.

## R CMD check results

0 errors | 0 warnings | 1 note

The sole note is CRAN incoming feasibility: New submission.
All 1,187 test expectations passed, with zero failures, warnings or skips.

The direct numerical comparison with 0.10.2 gave identical generated data
and zero maximum absolute differences for coefficients, covariance and forecasts.
The computational replication script reconstructed tables and figures and
reran worked examples under 0.10.3. Historical Monte Carlo rows were
reaggregated; their numerical results were not replaced or claimed as new runs.

## Downstream dependencies

No CRAN downstream dependencies are known for this first submission.

