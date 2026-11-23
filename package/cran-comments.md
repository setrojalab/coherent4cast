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

## Hosted native checks completed; R platform matrix in progress

The private setrojalab/recurcast repository contains the authorized candidate.
All four Linux standalone-native jobs pass: GCC/Clang strict warnings,
with and without AddressSanitizer and UndefinedBehaviorSanitizer.
Source commit: a0a11ce73d7c9360f7aba360959289268699b096.
Native run: https://github.com/setrojalab/recurcast/actions/runs/45156913689.

The R-release/devel matrix on Linux, macOS and Windows was started and is
still in progress as of this dated record. Final project evidence records its
actual completed results; a queued workflow is not a passing check.
Platform run: https://github.com/setrojalab/recurcast/actions/runs/45156913651.

Public submission, win-builder/macOS-builder service results and CRAN
acceptance are not claimed by these local or hosted tests.
