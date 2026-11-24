## Submission summary

This is a new submission of `recurcast` (version 0.10.1). The package fits
recurrent/terminal event-history models and returns joint-law dynamic
forecasts, observation-aware sensitivity bounds, targeting diagnostics, and
validation measures. It also implements baseline binary-treatment causal
effects through common-population g-formula standardization, cross-fitted
intensity targeting, and joint influence-function treatment contrasts under
explicit identification and regularity assumptions. Source is GPL-3.

* Maintainer: Sunday A. Adetunji <sundayadetunjisa@gmail.com>
* Submitted archive SHA-256:
  `cb246be1086c6b8970403df522d8e697d840354b5e0982e77e6ad762a26304dc`

## Test environments

* Final submitted tarball: Windows 11 x64, R-devel (2026-11-17 r90926,
  GCC 14.3.0), built under R 4.6.2, checked with `R CMD check --as-cran`
  with remote incoming checks and PDF/HTML manual generation enabled.
* Updated causal revision: Linux, macOS, and Windows under R-release and
  R-devel, with all six GitHub Actions package-check jobs passing.
* Native C++17 fixtures: four Linux GCC/Clang builds, with and without
  AddressSanitizer (ASAN) and UndefinedBehaviorSanitizer (UBSAN); all seven
  fixtures, including the added causal kernel, passed in every build.
* Hosted source commit: `f8228270a64bede49e92694ff097e00b656d1126`.
  The R code, native code, help files and tests in the submitted tarball
  match this tested causal revision; generated build/vignette metadata
  are checked in the final tarball separately. The hosted checks include
  the newly added causal R/C++ code.

## R CMD check results

0 errors | 0 warnings | 1 note

* `checking CRAN incoming feasibility ... NOTE`
  Maintainer: 'Sunday A. Adetunji <sundayadetunjisa@gmail.com>'
  New submission

All 1,187 test expectations passed (0 failures, warnings, or skips), along
with all examples, all five vignettes, and PDF/HTML manual checks.

## Downstream dependencies

There are no known downstream dependencies for this first submission.
