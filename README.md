# coherent4cast

**Coherent Survival Modelling and Event Forecasting**

An R package for survival modelling and dynamic forecasts of history-dependent recurrent and terminal events under incomplete observation. It provides a common probability law for recurrent-event burden, terminal-event risk and restricted mean time alive, with explicit observation assumptions and numerical diagnostics.

This is version **0.10.2**, renamed from `recurcast`. The S3/Rcpp/C++17 architecture and numerical kernels are retained. Historical releases and verification records remain available under their original names. The renamed package has not yet been accepted by CRAN.

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

## Migration

Change `library(recurcast)` to `library(coherent4cast)` and qualified calls to the new namespace. New scripts should use `coherent4cast()` and `coherent4cast_control()`; `recurcast()` and `recurcast_control()` remain aliases in this package. The `rc_*` function names are unchanged. Returned S3 class and audit labels use the new name; refit saved objects from the previous package before applying methods from the new namespace.

## Paper, documentation and source

- [Documentation website](https://setrojalab.github.io/coherent4cast/).
- [Software paper on arXiv](https://arxiv.org/abs/2611.32075); earlier versions retain their original title.
- [Renamed manuscript PDF](paper/coherent4cast/01_coherent4cast_JSS_arXiv_Manuscript.pdf) and [clean LaTeX source](paper/coherent4cast/coherent4cast_JSS_arXiv_LaTeX.zip).
- [CRAN source archive](release/coherent4cast/coherent4cast_0.10.2.tar.gz), [reference manual](release/coherent4cast/coherent4cast_0.10.2_Reference_Manual.pdf) and [submission comments](release/coherent4cast/cran-comments.txt).
- [Tutorial sources](package/vignettes), [replication scripts](replication) and [rename verification](reviews/rename_0_10_2).

## Verification

Fresh checks for the renamed package are recorded under `reviews/rename_0_10_2`; GitHub Actions runs the six Linux/macOS/Windows R-release/R-devel jobs and four GCC/Clang native jobs with and without ASAN/UBSAN. Historical check totals are not assigned to this release without rerunning them. The numerical kernel file hashes are compared with the preceding public source revision.

Authors: Sunday A. Adetunji and Oluwatimileyin A. Adetunji. Maintainer: sundayadetunjisa@gmail.com. Source license: GPL-3.
