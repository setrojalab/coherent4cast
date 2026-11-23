# recurcast

Coherent survival modelling and dynamic forecasting of history-dependent
multitype recurrent events and terminal risk under incomplete observation.

`recurcast` fits one event-history law and derives recurrent burden, terminal
risk and restricted survival from it. The R interface exposes a C++17 numerical
core, checks event eligibility and observation assumptions, and retains
convergence and numerical diagnostics in returned objects.

This unpublished candidate version is `0.10.0`. Local verification does not imply
CRAN acceptance, independent software review or clinical validation.

## Install the supplied source archive

Install Rcpp, then install the archive locally. A C++17 toolchain is required
for source installation; use Rtools on Windows.

```r
install.packages("Rcpp")
install.packages("recurcast_0.10.0.tar.gz", repos=NULL, type="source")
library("recurcast")
```

## Fit and forecast

```r
data("sim_perinatal")
fit <- recurcast(~x,sim_perinatal,event_types=c("1","2"),time_unit="months")
print(fit)
one <- sim_perinatal[sim_perinatal$id==1,]
predict(fit,one,landmark=0,horizon=c(.25,.5,1),type="all")
```

The three bundled datasets are synthetic. `rc_readmission()` and `rc_colon()`
load public clinical datasets from `frailtypack` and `survival`; real patient
records are not redistributed inside this package.

## Main interfaces

| Task | Interface |
|---|---|
| Validate intervals and construct predictable counts | `as_trajectory()`, `build_history()` |
| Fit gamma, lognormal or estimated discrete frailty | `recurcast(frailty=...)`, `rc_frailty_loglik()` |
| Generate and integrate coarsened histories | `rc_coarsen_fibre()`, `recurcast(observation="coarsened")` |
| Forecast from a posterior history mixture | `predict(newdata_observation="coarsened")`, `rc_predict_coarsened()` |
| Propagate probabilities and cumulative rewards | `predict()`, `rc_generator_csr()`, `rc_propagate_sparse()` |
| Cross-fit finite-path MAR or continuous-time CAR targets | `recurcast(targeting=..., targeting_spec=...)` |
| Compute scalable finite-support recording bounds | `rc_bounds()`, `rc_bounds_sparse()` |
| Integrate declared future covariates | `predict(covariate_evolution=...)`, `rc_covariate_model()` |
| Diagnose the frailty boundary and compare families | `summary()`, `rc_frailty_diagnose(alternatives=TRUE)` |
| Stabilize an eligible fit and simulate nonlinear bands | `rc_regularize()`, `rc_band()` |
| Invert weighted or certified conformal ranks | `rc_conformal()` |
| Compute a known-dispersion exact pivot | `rc_pivot()` |
| Evaluate held-out forecasts | `validate()` |

Use `vignette(package="recurcast")` for the three worked guides and the help
pages for precise assumptions. The finite-path, continuous-time likelihood,
confidence-projection and prediction-set modules have different validity
conditions. Calling one does not certify all the others.

## Reproducibility

Source fixtures and tests check likelihood scores, information, ordered
propagation, finite-path influence functions and optimization endpoints against
the preserved research implementation. The companion replication folder
contains scripts, actual benchmark outputs and case-study aggregates. Runtime
measurements are hardware dependent. See `GATES.md` for measured and pending
release checks.

The package records subject-level provenance and fails explicitly for invalid
histories, unsupported causal analyses and invalid numerical calculations.
It does not infer arbitrary hidden histories without an observation model.

## Licence and citation

The software is distributed under GPL-3. Use `citation("recurcast")` for the
development software citation. Publication identifiers will be added only after
they exist. See `CONTRIBUTING.md`, `SECURITY.md` and `SUPPORT.md` for development
and support policies.

