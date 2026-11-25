# coherent4cast

**Coherent Survival Modelling and Event Forecasting**

An R package for survival modelling and dynamic forecasts of history-dependent recurrent and terminal events under incomplete observation. It provides a common probability law for recurrent-event burden, terminal-event risk and restricted mean time alive, with explicit observation assumptions and numerical diagnostics.

This is version **0.10.3**. The package uses the S3/Rcpp/C++17 architecture. CRAN review is pending.

## Install

```r
install.packages(c("remotes", "Rcpp", "digest", "jsonlite"))
remotes::install_github("setrojalab/coherent4cast", subdir = "package",
  dependencies = NA, upgrade = "never")
library(coherent4cast)
data(sim_perinatal)
fit <- coherent4cast(~ x, sim_perinatal, time_unit = "months",
                     event_types = c("1", "2"))
print(fit)
```

Source installation requires a C++17 toolchain (the matching Rtools on Windows). The repository is public; no GitHub subscription is required. Use `build_vignettes = TRUE` with the vignette dependencies installed to build the five tutorials.

## Interfaces

| Task | Functions |
|---|---|
| Event-history validation and model fitting | `as_trajectory()`, `coherent4cast()`, `coherent4cast_control()` |
| Dynamic forecasting | `predict()` |
| Standard survival and competing risks | `rc_surv()`, `rc_competing()` |
| Baseline treatment contrasts | `rc_causal()`, `rc_causal_contrast()` |
| Cross-fitted targeting | `rc_target()` and supported `targeting` specifications |
| Recording sensitivity | `rc_bounds()`, `bounds()` |
| Held-out prediction assessment | `validate()` |
| Provenance and diagnostics | `rc_audit_report()` |

The baseline causal interface requires consistency, conditional exchangeability, treatment positivity and independent censoring. It does not estimate longitudinal treatment policies. Forecast uncertainty, recording-sensitivity bounds and causal contrasts retain their distinct assumptions.

## Entry points

Use `library(coherent4cast)`, `coherent4cast()` and `coherent4cast_control()`. The `rc_*` interfaces retain their names. Refit objects saved by earlier development packages before using this namespace.

## Paper, documentation and source

- [Documentation website](https://setrojalab.github.io/coherent4cast/).
- [Software paper on arXiv](https://arxiv.org/abs/2611.32075); earlier versions retain their original title.
- [Manuscript PDF](paper/coherent4cast/01_coherent4cast_JSS_arXiv_Manuscript.pdf) and [clean LaTeX source](paper/coherent4cast/coherent4cast_JSS_arXiv_LaTeX.zip).
- [CRAN source archive](release/coherent4cast/coherent4cast_0.10.3.tar.gz), [reference manual](release/coherent4cast/coherent4cast_0.10.3_Reference_Manual.pdf) and [submission comments](release/coherent4cast/cran-comments.txt).
- [Complete computational replication archive](release/coherent4cast/coherent4cast_JSS_replication.zip).
- [Tutorial sources](package/vignettes), [replication scripts](replication) and [release verification](reviews/release_0_10_3).

## Verification

Version 0.10.3 passed a fresh R CMD check --as-cran with 0 errors, 0 warnings and 1 new-submission note, including 1,187 test expectations and five vignettes. All six Linux/macOS/Windows R-release/R-devel jobs and all four GCC/Clang native jobs with and without ASAN/UBSAN passed. [Verification records](reviews/release_0_10_3) identify the exact tested revision and archive checksum. Numerical comparisons with 0.10.2 gave identical coefficients, covariance estimates and forecasts. The source, archive members and rendered PDF text were scanned for obsolete package identifiers.

Authors: Sunday A. Adetunji and Oluwatimileyin A. Adetunji. Maintainer: sundayadetunjisa@gmail.com. Source license: GPL-3.


