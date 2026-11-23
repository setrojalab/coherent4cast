## Development status

This package has not been submitted to CRAN. This file records actual checks
and must be updated before submission. A workflow definition is not an executed
cross-platform result.

## Preserved 0.10.0 baseline checks (0.10.1 checks pending)

Windows 11 x64, GCC 14.3.0 from Rtools45. Source archives were checked with
`R CMD check --as-cran --no-manual` under isolated R 4.6.2 release and R-devel
4.7.0 (2026-11-17 r90926) installations. Both final frozen-source checks have
0 errors, 0 warnings and 0 notes. Each reports 917 passed test expectations,
0 failures, 0 warnings and 0 skips. All four vignettes are built and included,
and their outputs rebuild successfully during checking. A standard Rd PDF
46-page reference manual was generated separately. The Windows binary was
installed in a fresh library and passed the same 917 assertions. Combined
R/C++ line coverage is 95.28487%, with no excluded files.

Exact logs, archive hashes and source-content comparisons are recorded in
the project-level `reviews/P6_candidate0100_checks.json` and
`reviews/P6_frozen_archive_sources.csv`. These are local Windows results,
not CRAN acceptance or evidence from other operating systems.

## External checks still required

Linux and macOS; sanitizer/valgrind execution; win-builder and macOS builder;
confirmation of maintainer/author authorization for public submission.
The core test and R package CI workflows are provided but no remote run has
been represented as completed.



