## Candidate status

Unpublished candidate 0.10.1. This package has not been submitted to CRAN.

## Local source and binary checks completed 2026-11-23

Windows 11 x64, Rtools45 GCC 14.3.0. R CMD check --as-cran --no-manual
reports 0 errors, 0 warnings and 0 notes under R 4.6.2 release and R-devel 4.7.0
(2026-11-17 r90926). Each source check executes 1107 passing test expectations,
with 0 failures, 0 warnings and 0 skips. A fresh Windows binary installation passes
the same 1107 assertions. The checked source archives include all four compiled
vignettes and their index; vignette outputs are rebuilt during checking.

The separately generated 54-page Rd reference manual has been rendered and
inspected. Combined R/C++ line coverage is 95.46236% (2777/2909 lines), with
no file exclusions. All 23 native files match the unchanged 0.10.0 baseline.

Exact logs, source comparisons and SHA256 checksums are in the project-level
reviews/P6_candidate0101_checks.json, reviews/P6_frozen_archive_sources.csv
and reviews/P6_final_verification.json.

## Hosted platform and native checks completed

All six hosted R-release/devel checks on Linux, macOS and Windows report
Status: OK, with 0 errors, 0 warnings, 0 notes and 1107 passing assertions
per job. All four Linux standalone-native jobs pass: GCC/Clang strict warnings,
with and without AddressSanitizer and UndefinedBehaviorSanitizer.
Tested source commit: fbf6429e37aa91e159d878a8228d62803b287056.
Platform run: https://github.com/setrojalab/recurcast/actions/runs/45164392545.
Native run: https://github.com/setrojalab/recurcast/actions/runs/45164392550.

The earlier macOS R-devel Matrix dependency failure occurred before recurcast
was checked. The workflow-only gettext environment correction resolved it;
the original failed log and correction are retained. A subsequent attempt
stalled in the oversized repository-status diagnostic before dependency
installation. Source-only CRAN configuration on macOS devel avoids probes of
nonexistent binary and unused Bioconductor repositories; no dependencies or
checks were removed. That superseded log is also retained. Package scientific
source and the local checked artifacts were unchanged.

These are package/test/example/vignette checks, not full cross-platform
clinical case-study replication or bitwise equality. Public submission,
win-builder/macOS-builder service results and CRAN acceptance are not claimed.
