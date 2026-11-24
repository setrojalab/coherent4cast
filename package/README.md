# recurcast

Coherent survival modelling and dynamic forecasting of history-dependent
multitype recurrent events and terminal risk under incomplete observation.

`recurcast` fits one event-history law and derives recurrent burden, terminal
risk and restricted survival from it. The R interface exposes a C++17 numerical
core, checks event eligibility and observation assumptions, and retains
convergence and numerical diagnostics in returned objects.

This unpublished candidate version is `0.10.1`. Local verification does not imply
CRAN acceptance, independent software review or clinical validation.

## Install the supplied source archive

Install Rcpp, digest and jsonlite, then install the archive locally. A C++17 toolchain is required
for source installation; use Rtools on Windows.

```r
install.packages(c("Rcpp", "digest", "jsonlite"))
install.packages("recurcast_0.10.1.tar.gz", repos=NULL, type="source")
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

`rc_causal()` adds the static binary baseline-treatment extension from
Supplement S12. It standardizes both arm-specific event laws over the same
baseline population and returns contrasts with their joint EIF covariance.
`rc_causal_contrast()` provides differences, ratios and log ratios. These
outputs require consistency, conditional treatment exchangeability,
positivity and conditionally independent censoring. They are distinct from
ordinary prognostic forecasts. See `vignette("baseline-causal")` for the
supported finite-state/stratum model and worked example.


| Task | Interface |
|---|---|
| Validate intervals and construct predictable counts | `as_trajectory()`, `build_history()` |
| Fit zero-variance, gamma, lognormal or estimated discrete frailty | `recurcast(frailty=...)`, `rc_frailty_loglik()` |
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
| Fit standard survival or absorbing competing causes | `rc_surv()`, `rc_competing()` |
| Export captured context, diagnostics and explicit SHA-256 fingerprints | `rc_audit_report()` |

Use `vignette(package="recurcast")` for the four worked guides and the help
pages for precise assumptions. The finite-path, continuous-time likelihood,
confidence-projection and prediction-set modules have different validity
conditions. Calling one does not certify all the others.

## Targeted 0.10.1 interfaces

`frailty = "none"` fixes latent variance at zero. `rc_surv()` fits the full
piecewise-exponential proportional-hazards likelihood, not Cox partial
likelihood. `rc_competing()` fits cause-specific no-frailty models and combines
all hazards in one generator to calculate cumulative incidence and survival.
It does not automatically produce CIF confidence bands. Competing-event
validation weights earlier competing events as observed noncases at their
event time, not as censoring.

`rc_trajectory` and `rc_forecast` are additive classes; use `as_trajectory()`
and `predict()` rather than treating those class names as functions. Explicit
absorbing-cause labels and observed-pattern metadata do not infer scientific
eligibility or the missingness mechanism. `rc_audit_report()` separates
captured context from report export and states each hash representation;
it does not confer anonymity, causal identification or regulatory approval.

The project-level `replication/standard_regimes.R` and
`replication/standards_mapping.R` scripts give executable wrapper and synthetic
ADTTE/OMOP-style mapping examples. The latter are analyst-controlled examples,
not native certified adapters. Offline execution requires dependencies to be
installed in advance.

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

