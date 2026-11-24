# recurcast 0.10.1

Coherent survival modelling and dynamic forecasting of history-dependent multitype recurrent events and terminal risk under incomplete observation.

**Public source release:** [v0.10.1](https://github.com/setrojalab/recurcast/releases/tag/v0.10.1). This version includes the static binary-treatment extension in Supplement S12 of the method paper, alongside the original likelihood, forecasting, targeting and sensitivity interfaces. The earlier prognostic numerical kernels are unchanged; a separate causal targeting kernel is added.

The source archive has been submitted to CRAN. Acceptance is pending; install from the versioned GitHub release below in the meantime.

## Install version 0.10.1

This repository is public; no GitHub subscription or access token is required. Source installation requires a C++17 toolchain, including the appropriate Rtools on Windows.

```r
install.packages(c("remotes", "Rcpp", "digest", "jsonlite"))
remotes::install_github("setrojalab/recurcast", subdir = "package",
  ref = "v0.10.1", dependencies = NA, upgrade = "never")
library(recurcast)
```

The current source archive, reference manual and submission comments are under [`release/causal/`](release/causal/). Use those files for this revision. Older archives elsewhere in `release/`, including the previous Windows binary, retain their historical provenance and do not include this causal extension.

## Baseline treatment contrasts

`rc_causal()` estimates recurrent-event burden, terminal-event risk and restricted mean time alive under treatment 1 and treatment 0, standardized over the same baseline population. `rc_causal_contrast()` returns differences, ratios or log ratios using the full joint influence-function covariance.

```r
d <- rc_simulate(200, horizon = 1, rates = c(.9, .5, .4),
                 dispersion = 0, censor_rate = .1, seed = 701)
d$treatment <- as.integer(d$id %% 2)
d$censored <- as.integer(d$event_type == 0 & d$terminal == 0 & d$stop < 1)
fit <- rc_causal(~ x, d, treatment = "treatment", horizon = 1,
  censor = "censored", time_unit = "years", propensity = .5,
  spec = list(history_cap = 2L, working = "loglinear"))
predict(fit)
fit$diagnostics
rc_audit_report(fit)
```

This example is synthetic and contains no treatment effect in its generator. See the [causal vignette](package/vignettes/baseline-causal.Rmd) and [reference help](package/man/rc_causal.Rd) for data contracts and assumptions. The implementation supports exact histories, finite baseline strata and a retained count state. Treatment is assigned at baseline and remains fixed. Causal identification requires consistency, conditional exchangeability, treatment positivity and independent censoring. The regular intervals are asymptotic. Time-varying treatment policies and unmeasured treatment-confounding bounds are not estimated by this entry point.

## Existing interfaces

| Task | Functions |
|---|---|
| Recurrent/terminal event histories and forecasts | `as_trajectory()`, `recurcast()`, `predict()` |
| Standard survival and competing risks | `rc_surv()`, `rc_competing()` |
| Baseline causal effects | `rc_causal()`, `rc_causal_contrast()` |
| Recording sensitivity | `rc_bounds()`, `bounds()` |
| Held-out prediction assessment | `validate()` |
| Provenance and numerical diagnostics | `rc_audit_report()` |

The prognostic APIs retain their original interpretation. Recording sensitivity and treatment confounding are separate problems. Synthetic agreement and numerical error tolerances do not establish clinical benefit or regulatory certification.

## Paper and reproducibility

- [Updated software paper](paper/causal/01_RecurCast_JSS_arXiv_Manuscript.pdf) and [clean JSS/arXiv LaTeX](paper/causal/RecurCast_JSS_arXiv_LaTeX.zip).
- [Current source package](release/causal/recurcast_0.10.1.tar.gz), [reference manual](release/causal/recurcast_0.10.1_Reference_Manual.pdf), and [CRAN comments](release/causal/cran-comments.txt).
- [Computational supplement](release/causal/RecurCast_JSS_replication.zip): execute `Rscript --vanilla code.R` after extracting and installing the package. This reproduces the original tables and figures and executes the added causal example and analytic remainder check.
- [Causal replication source](replication/causal/run_causal.R) and [independent matrix-exponential oracle](replication/causal/causal_oracle.py).
- [SHA-256 manifest](release/causal/SHA256SUMS.txt).

The package contains five evaluated vignettes in its checked source archive. A direct GitHub installation requires `build_vignettes = TRUE` and the listed vignette dependencies if evaluated vignettes are wanted locally.

## Verification

The final source archive passed full Windows R-devel `R CMD check --as-cran`: **0 errors, 0 warnings, 1 new-submission note**, including all **1,187 expectations**, five vignettes, examples and PDF/HTML manual checks. The new code was checked at commit `f8228270a64bede49e92694ff097e00b656d1126`: six Linux/macOS/Windows R-release/R-devel jobs passed, as did four GCC/Clang jobs with and without ASAN/UBSAN. The seven native fixtures include independent analytic causal EIF agreement and exact probability-one reduction to the preceding targeting engine. See the [record](reviews/causal_revision/summary.json), [R workflow](https://github.com/setrojalab/recurcast/actions/runs/45324158268) and [native workflow](https://github.com/setrojalab/recurcast/actions/runs/45324158206).

The public-release update changes availability documentation and article links; it does not change the tested R or C++ implementation. All 158 audited package source, help, test, vignette and data files match the submitted archive after line-ending normalization.

The current causal module has not been assigned the earlier version's measured line-coverage percentage. Historical research grids and clinical analyses remain labelled as such; the five original application examples remain prognostic.

Authors: Sunday A. Adetunji and Oluwatimileyin A. Adetunji. Maintainer: sundayadetunjisa@gmail.com. License: GPL-3. CRAN submission is pending review. Journal acceptance and an arXiv identifier have not been assigned.
