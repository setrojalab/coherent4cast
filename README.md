# coherentSurv

**Coherent Survival Modelling and Dynamic Forecasting for Recurrent and Terminal Events in R and Python**

An R package for survival modelling and dynamic forecasts of history-dependent recurrent and terminal events under incomplete observation. It provides a common probability law for recurrent-event burden, terminal-event risk and restricted mean time alive, with explicit observation assumptions and numerical diagnostics.

This is version **0.11.1**. The package uses the S3/Rcpp/C++17 architecture. This release adds diagnostic safeguards and mapping helpers without changing the numerical kernels. The installation command below selects the GitHub release. The 0.11.1 CRAN archive has been prepared and checked; CRAN acceptance is not asserted.

## Install

```r
install.packages(c("remotes", "Rcpp", "digest", "jsonlite"))
remotes::install_github("setrojalabs/coherentSurv", subdir = "package",
  ref = "v0.11.1", dependencies = NA, upgrade = "never")
library(coherentSurv)
data(sim_perinatal)
fit <- cs_fit(~ x, sim_perinatal, time_unit = "months",
                     event_types = c("1", "2"))
print(fit)
```

Source installation requires a C++17 toolchain (the matching Rtools on Windows). The repository is public; no GitHub subscription is required. Use `build_vignettes = TRUE` with the vignette dependencies installed to build the tutorials: ten in version 0.11.1.

The numerical implementation described in the revised paper was tested at source revision `1f39e037c0e8bf24edd58366dfe46c2e5546ed94`. The release adds documentation and metadata corrections; the tagged source is the recommended installation. To reproduce the original verification snapshot explicitly:

```r
remotes::install_github("setrojalabs/coherentSurv", subdir = "package",
  ref = "1f39e037c0e8bf24edd58366dfe46c2e5546ed94",
  dependencies = NA, upgrade = "never")
```

## Interfaces

| Task | Functions |
|---|---|
| Event-history validation and model fitting | `cs_as_trajectory()`, `cs_fit()`, `cs_control()` |
| Dynamic forecasting | `predict()` |
| Standard survival and competing risks | `cs_surv()`, `cs_competing()` |
| Baseline treatment contrasts | `cs_causal()`, `cs_causal_contrast()` |
| Cross-fitted targeting | `cs_target()` and supported `targeting` specifications |
| Recording sensitivity | `cs_bounds()`, `bounds()` |
| Held-out prediction assessment | `validate()` |
| Provenance and diagnostics | `cs_audit_report()` |

The baseline causal interface requires consistency, conditional exchangeability, treatment positivity and independent censoring. It does not estimate longitudinal treatment policies. Forecast uncertainty, recording-sensitivity bounds and causal contrasts retain their distinct assumptions.

## Entry points

Use `library(coherentSurv)`, `cs_fit()` and `cs_control()`. Package-specific interfaces now use the `cs_` prefix. This is a naming change to the development API: update calling scripts and refit objects saved by earlier development packages. Standard R generics retain their names. See the [naming and interoperability conventions](NAMING_AND_INTEROPERABILITY.md).

## Use from Python

Install R and this R package first. The [rpy2 example](package/inst/python/run_rpy2.py) calls `coherentSurv::cs_fit()` and R's `predict()` from Python, returns a pandas table, exports diagnostics, and checks numerical agreement with direct R calls. It requires a supported rpy2 environment with pandas and numpy. The [Rscript runner](package/inst/python/run_forecast.py) provides process-based access using Python's standard library. Both routes use the same R/C++ statistical engine. Python interface checks are recorded separately from the statistical validation; no separate PyPI distribution or scikit-learn pipeline compatibility is implied.

## Paper, documentation and source

- [Documentation website](https://setrojalabs.github.io/coherentSurv/).
- [Software paper on arXiv](https://arxiv.org/abs/2611.32075); replacement 8431196 was submitted with the final title on 26 November 2026 UTC; the public abstract can retain its earlier title until arXiv announces the replacement.
- [Manuscript PDF](paper/coherentSurv/01_coherentSurv_JSS_arXiv_Manuscript.pdf) and [clean LaTeX source](paper/coherentSurv/coherentSurv_JSS_arXiv_LaTeX.zip).
- [CRAN source archive](release/coherentSurv/coherentSurv_0.11.1.tar.gz), [reference manual](release/coherentSurv/coherentSurv_0.11.1_Reference_Manual.pdf) and [submission comments](release/coherentSurv/cran-comments.txt).
- [Complete computational replication archive](release/coherentSurv/coherentSurv_JSS_replication.zip).
- [Tutorial sources](package/vignettes), [replication scripts](replication) and [release verification](reviews/release_0_11_1).

## Version 0.11.1 verification materials

The following preserve the reviewed artifacts used for the arXiv replacement. The current release downloads above include the subsequent documentation-only metadata reconciliation.

- [Checked source archive](release/coherentSurv/0.11.1-review/coherentSurv_0.11.1.tar.gz), [reference manual](release/coherentSurv/0.11.1-review/coherentSurv-reference-manual.pdf), and [check summary](release/coherentSurv/0.11.1-review/cran-comments.txt).
- [Complete replication archive](release/coherentSurv/0.11.1-review/coherentSurv_JSS_replication.zip).
- [Technical review and 40-item scope register](reviews/release_0_11_1/Review_and_40_Item_Register.pdf).
- [Final verification record](reviews/release_0_11_1/verification/final_status.json) and [artifact checksums](reviews/release_0_11_1/SHA256SUMS.txt).

The tested source revision is `1f39e037c0e8bf24edd58366dfe46c2e5546ed94`. Its full local R-devel check completed with 0 errors, 0 warnings and 1 new-submission NOTE; all 1,476 expectations and ten vignettes passed. Seven platform jobs, four native compiler/sanitizer configurations and the Linux rpy2 check passed. These are package and interface checks, not regulatory certification or new statistical coverage experiments. The 0.11.1 archive has not been submitted to CRAN. The final release metadata reconciliation and archive checksum are recorded in [release verification](reviews/release_0_11_1/publication_status.json).

## Verification

The published 0.11.0 checks are recorded in [release verification](reviews/release_0_11_0), with the tested source revision and archive checksum. The local full R-devel check passed all 1,197 expectations, examples, five vignettes and PDF/HTML manuals. Six Linux/macOS/Windows package checks and four native compiler/sanitizer configurations passed. The Linux rpy2 example matched direct R coefficients, covariance and forecasts exactly on the tested fixture; the Rscript example ran separately on Windows. Package checks, data-transfer agreement and statistical validation are reported separately. The [0.11.1 review record](reviews/release_0_11_1) records the new source, completed tests and remaining development scope separately.

## Applying the event-history contract

| Contract | Clinical study | Reliability | Credit study |
|---|---|---|---|
| Independent unit | Participant/patient | Device | Account/borrower, as declared |
| Recurrent marks | Admissions or toxicity episodes | Faults/repairs | Delinquency/covenant episodes |
| Absorbing endpoint | Death | Irreversible failure | Default, if absorbing by definition |
| Expected burden | Episodes before death | Episodes before failure | Episodes before default |
| Restricted time alive | Time alive within the horizon | Time before failure within the horizon | Time before default within the horizon |

These are translations of a data contract, not independently validated studies.
Restricted time alive is not automatically event-free time or quality-adjusted
survival. Monetary loss requires a separate exposure and reward definition.
Administrative loss to follow-up is not an absorbing outcome.

Runnable tutorials in this release:

- [Cross-discipline guide](package/vignettes/cross-discipline-guide.Rmd)
- [Robotic fleet reliability](package/vignettes/robotic-fleet-reliability.Rmd)
- [Clinical-trial estimands](package/vignettes/clinical-trial-estimands.Rmd)
- [Credit risk and default](package/vignettes/credit-risk-default.Rmd)
- [Source mapping and environment diagnostics](package/vignettes/source-mapping-and-environment.Rmd)

Python users call the same R/C++ engine. See the [Python access guide](package/inst/python/README.txt).
The rpy2 runner fits once by default; `--verify-parity` requests a second direct-R
fit, and `--timings` reports observed stages. No native Python wheel, zero-copy
conversion, fastest-runtime claim or regulatory certification is implied.

Authors: Sunday A. Adetunji and Oluwatimileyin A. Adetunji. Maintainer: sundayadetunjisa@gmail.com. Source license: GPL-3.


