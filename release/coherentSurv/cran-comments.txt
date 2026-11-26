## Submission summary

coherentSurv version 0.11.1 provides coherent survival modelling and dynamic
forecasting for history-dependent recurrent and terminal events, with Python
access through rpy2 or Rscript. This revision adds guarded inference reporting,
scoped data mappings, environment diagnostics and ten executable tutorials.
The existing C++17 numerical kernels are unchanged. Licence: GPL-3.

Maintainer: Sunday A. Adetunji <sundayadetunjisa@gmail.com>
Source archive SHA-256: fb6f24384908b1b4869dd7a7d6d03b57166e889a7ae9153c72e925efd6c4c469

## Test environments

The exact release archive was rebuilt and checked on Windows with R-devel using
R CMD check --as-cran, including remote incoming checks and PDF/HTML manuals.
The release metadata commit is d4f41b95a738594d502e80a9ff21c6657f6c9506.
Executed GitHub platform, native compiler/sanitizer and rpy2 results are
recorded with their source revisions in reviews/release_0_11_1/.
The numerical implementation is unchanged from the verified 0.11.1 source
1f39e037c0e8bf24edd58366dfe46c2e5546ed94. The release reconciliation changes
documentation and metadata; archive-member comparisons are supplied.

## R CMD check results

0 errors | 0 warnings | 1 note

The incoming-feasibility NOTE identifies the maintainer and new submission.
All 1,476 expectations passed: 0 failures, warnings or skips. All examples,
ten vignettes and PDF/HTML manual checks passed.

## Submission status

These are completed pre-submission checks, not a claim of CRAN acceptance.
Version 0.11.1 has not been submitted to CRAN in this release task.
There are no known downstream dependencies for this first CRAN submission.
