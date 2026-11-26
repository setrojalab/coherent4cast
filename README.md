# coherentSurv

**Coherent Survival Modelling and Dynamic Landmark Forecasting**

An R package for survival modelling and dynamic forecasts of history-dependent recurrent and terminal events under incomplete observation. It provides a common probability law for recurrent-event burden, terminal-event risk and restricted mean time alive, with explicit observation assumptions and numerical diagnostics.

This is version **0.11.0**. The package uses the S3/Rcpp/C++17 architecture. This naming revision is being prepared for CRAN; acceptance is not asserted.

## Install

```r
install.packages(c("remotes", "Rcpp", "digest", "jsonlite"))
remotes::install_github("setrojalab/coherentSurv", subdir = "package",
  ref = "v0.11.0", dependencies = NA, upgrade = "never")
library(coherentSurv)
data(sim_perinatal)
fit <- cs_fit(~ x, sim_perinatal, time_unit = "months",
                     event_types = c("1", "2"))
print(fit)
```

Source installation requires a C++17 toolchain (the matching Rtools on Windows). The repository is public; no GitHub subscription is required. Use `build_vignettes = TRUE` with the vignette dependencies installed to build the five tutorials.

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

- [Documentation website](https://setrojalab.github.io/coherentSurv/).
- [Software paper on arXiv](https://arxiv.org/abs/2611.32075); earlier versions retain their original title.
- [Manuscript PDF](paper/coherentSurv/01_coherentSurv_JSS_arXiv_Manuscript.pdf) and [clean LaTeX source](paper/coherentSurv/coherentSurv_JSS_arXiv_LaTeX.zip).
- [CRAN source archive](release/coherentSurv/coherentSurv_0.11.0.tar.gz), [reference manual](release/coherentSurv/coherentSurv_0.11.0_Reference_Manual.pdf) and [submission comments](release/coherentSurv/cran-comments.txt).
- [Complete computational replication archive](release/coherentSurv/coherentSurv_JSS_replication.zip).
- [Tutorial sources](package/vignettes), [replication scripts](replication) and [release verification](reviews/release_0_11_0).

## Verification

Current checks are recorded in [release verification](reviews/release_0_11_0), with the tested source revision and archive checksum. Earlier package checks do not substitute for validation of this naming revision. Package, native, platform, Python-interface and numerical-parity results are reported separately.

Authors: Sunday A. Adetunji and Oluwatimileyin A. Adetunji. Maintainer: sundayadetunjisa@gmail.com. Source license: GPL-3.


