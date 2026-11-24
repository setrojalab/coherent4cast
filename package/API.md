# Public API and statistical scope

This development API implements the accompanying RecurCast method.
The current 0.10.1 causal revision adds `rc_causal()` and
`rc_causal_contrast()` from Version 13 Supplement S12. The older theorem
references in the retained prognostic table below keep their Version 5 provenance. The theorem
numbers below refer to the supplied, unpublished Version 5 manuscript; they are
not citations to an externally published result. All time units and forecast
horizons are explicit. Nonregular fits retain visible diagnostics.

| Interface | Input and result | Mathematical or computational basis | Evidence |
|---|---|---|---|
| `as_trajectory`, `build_history` | Counting-process intervals; validated schema and predictable counts, recency and gaps | Left-continuous history; no future event in a current covariate | Dirty-data tests; compiled history tests |
| `recurcast`, `recurcast_control` | One-sided formula and complete or automatically expanded coarsened histories; fitted common law | Gamma, lognormal or estimated discrete frailty, analytic scores and Louis information; explicit targeting dispatch | Gamma oracle; alternative-frailty derivative, integration and refinement checks |
| `rc_loglik` | Expanded risk rows; individual likelihood, score and information | Closed-form gamma Laplace integration | Python oracle errors below declared tolerances |
| `rc_coarsen` | Legal path likelihoods and integration weights; observed-data likelihood and Louis information | Sum over the observation fibre; no discarded tied records | Independent mixture calculation and `colon` adapter parity |
| `predict`, `augment`, `rc_propagate` | Landmark, durations and sufficient count state; terminal risk, burden, survival and restricted mean time alive | Ordered uniformization; one frailty value retained across segments, then integrated; Theorem 3 | Mass, semigroup, analytic constant-rate and quadrature tests |
| `rc_filter_paths`, `rc_predict_paths`, `rc_predict_coarsened` | Explicit or automatically constructed pre-landmark compatible paths; posterior mixture forecasts | Bayes normalization and family-specific frailty posterior, aligned future-covariate integration draws | Direct averaging, original-ID preservation, future-information rejection and measured-covariate tests |
| `rc_target` | Subject-level finite-path records, folds, strata, path features; targeted estimates, EIF matrix and covariance | Cross-fitted finite-path exponential targeting, MAR recording positivity; equations (6)-(7) | Fixed Python estimates and EIF fixtures |
| `rc_bounds`, `bounds` | Recorded finite-support law, coarsening map, sensitivity strength, target numerator/denominator | Charnes-Cooper linear-fractional optimization; equation (13). Optional simultaneous cell confidence projection | Vertex and recorded-law projection fixtures; checked duality gaps |
| `rc_regularize` | Unweighted exact-history fit with estimated gamma dispersion; stabilized fit and sandwich covariance | Local positive observed-information half-log-determinant penalty; equations (14)-(15) | Independent archived analytic trace score and numerical derivative tests |
| `rc_band` | Parameter vector, covariance and forecast callback; simultaneous transformed Gaussian bands | Frequentist sandwich-Gaussian max-t approximation; equation (16) | Input, boundary, random-state and quantile tests. No automatic finite-sample guarantee |
| `rc_conformal`, `conformal` | Held-out scores, candidate scores and weights or certified weight bands | Candidate-inclusive weighted rank inversion; Theorems F3-F4 | Weighted rank fixtures and explicit estimated-weight disclaimer |
| `rc_pivot` | Uncensored independent terminal times and known dispersion | Exact gamma-frailty pivot; Theorem F6 and equation (17) | Monotone pivot inversion tests |
| `rc_cornish_fisher` | Probability, sample size and cumulants | Model-specific second-order standardized-mean quantile expansion; Theorem F5 | Independent algebraic checks; approximation scope retained |
| `rc_identify` | Analyst-supplied population contraction and nonunit contrast | Population identification diagnostic; Theorem F1 | Observable contraction example; numerical saturation flag |
| `validate` | Held-out fixed-horizon risks and censored outcomes | IPCW Brier, cumulative/dynamic AUROC, calibration, threshold statistics | Hand-calculated 2-by-2 example, censoring positivity and undefined-metric tests |
| `rc_readmission`, `rc_colon` | Canonical source-package data loaded by the user | Declared clinical event dictionary, time conversion and legal coarsening | 403-person readmission and 929-person colon audits |
| `rc_simulate` and three bundled datasets | Seeded synthetic common-law trajectories | Shared-gamma simulation with predictable capped history and independent censoring | Seed reproducibility, terminal stopping and global RNG restoration tests |
| `coef`, `vcov`, `logLik`, `summary`, `tidy`, `glance`, `plot`, `autoplot` | Structured model, inference and forecast objects | Presentation of the underlying result without changing its estimand | S3 and workflow tests |

Version 0.10.0 adds the following integrated interfaces:

| Interface | Implemented experiment | Verification and qualification |
|---|---|---|
| `rc_coarsen_fibre`, `rc_colon_fibre`, `rc_fibre_mcse` | Legal point/window histories with exact order enumeration and quadrature, or disclosed importance sampling above a permutation budget | Duplicate histories excluded; terminal ordering enforced; all 929 colon participants retained |
| `rc_generator_csr`, `rc_propagate_sparse` | Ordered sparse state propagation and integrated rewards | Dense-oracle agreement, series/roundoff bounds, 50,000-state benchmark |
| `rc_bounds_sparse` and automatic `rc_bounds` dispatch | Exact structural reduction for the existing odds/coarsening fractional problem | 80 independent HiGHS comparisons; 50,000 atoms; sparse fallback records dense numerical failures |
| `recurcast(targeting="continuous")` | Cross-fitted continuous-time sufficient-state CAR intensity targeting | Actual-time martingale EIF and common-law update; oracle agreement and 20,000 separately registered cohorts, including retained sparse failures |
| `recurcast(targeting="finite_path")` | Explicit finite-support subject-level MAR targeting | Existing finite-path estimator exposed through fit/predict; score and regularity diagnostics |
| `predict(..., interval="projection")` on targeted fits | Conservative observed-reward confidence outer sets | Hoeffding/Bonferroni proof and 15,000 separately registered cohorts; uncapped recurrent burden has outer support [0,infinity] |
| `rc_covariate_model`, `predict(covariate_evolution=...)` | Frozen, scheduled or Gaussian longitudinal future trajectories | Explicit measured-history conditioning, ordered propagation, Monte Carlo error and grid checks; no latent measurement-error joint-survival claim |
| `rc_frailty_loglik`, `rc_frailty_diagnose` | Gamma, refined lognormal or estimated mean-one discrete frailty | Analytic derivatives/Louis checks, adaptive quadrature audit, boundary score and descriptive nonnested AIC |

Each targeting mode declares its own observation model through `targeting_spec`.
The continuous engine is not an unrestricted MNAR or hidden-frailty efficient
estimator. The data-dependent automatic Gaussian/projection switch does not
inherit a uniform finite-sample guarantee: request `interval="projection"`
throughout when invoking that confidence construction. The existing prognostic `intervention` argument remains reserved.
Use the separate `rc_causal()` entry point for static binary baseline effects.
Python bindings, native FHIR/OMOP adapters and time-varying causal regimes
retain their separate later-release scope. The project now includes explicit
synthetic OMOP-style mapping examples; these are not native adapters.

The initial 0.10.1 interface additions below preserved the native kernels.
The causal revision adds a separate C++ kernel and generated registration;
the earlier prognostic kernels remain unchanged:

| Interface or class | Implemented scope | Qualification |
|---|---|---|
| `recurcast(frailty="none")` | Explicit zero-variance boundary | Retains regression and the full likelihood; no latent frailty integration. |
| `rc_surv()` | Standard survival through a no-frailty, piecewise-exponential proportional-hazards likelihood | Supports declared exact contiguous risk intervals from origin; not Cox partial likelihood or delayed-entry correction. |
| `rc_competing()` | Cause-specific no-frailty hazards combined into common-generator CIFs and survival | Explicit mutually exclusive absorbing causes; joint subject-score sandwich; no automatic CIF confidence bands. |
| `validate()` on `recurcast_competing` | One selected cause's fixed-horizon cumulative incidence | Earlier competing events are noncases weighted at event time; censoring positivity/independence assumptions remain. |
| `rc_trajectory`, `rc_forecast` | Additive data and prediction classes with legacy inheritance | These are class names, not new estimator functions. Pattern metadata do not infer a scientific model. |
| `rc_audit_report()` | List, JSON or text report with captured context, diagnostics and explicit source-file hashes | SHA-256 states its serialized-object/file-byte representation; export is distinct from fitting and does not certify regulatory compliance. |

`replication/standards_mapping.R` in the project archive provides small synthetic
ADTTE/OMOP-style transformations and checks. Selection of parameters, episode
rules, observation origins and terminal ascertainment remains explicit. These
examples neither validate arbitrary clinical coding nor claim a certified adapter.

The exported C++ headers contain no R headers. Rcpp bridges handle conversion
and registration. Validation of the data-frame schema is implemented in R;
history compilation and numerical inference kernels are implemented in C++17.
The simulation uses R's local L'Ecuyer-CMRG stream with caller-state restoration;
the separate tested PCG32 header supplies a future binding-compatible stream.
Identical seeds in unrelated language RNG implementations are not claimed to
produce identical data.


## Baseline causal extension (Version 13, Supplement S12)

`rc_causal()` cross-fits event intensities, treatment probabilities and censor
survival, then targets arm-specific intensities using the treatment- and
censoring-weighted backward-value jump contrasts. All subjects contribute to
the baseline standardization of both arms. `rc_causal_contrast()` uses the
joint subject-level EIF covariance for differences, ratios and log ratios.

The targets are expected recurrent counts, terminal risk and restricted mean
time alive under a static baseline binary treatment. Consistency, conditional
exchangeability, treatment positivity, independent censoring and sufficient
retained state are explicit assumptions. Event-model correctness or correct
treatment-probability/censor-survival product gives the stated remainder
property; regular efficient inference requires the additional rate conditions.

The implementation uses exact contiguous histories, at most 30 discrete
baseline strata, six recurrent marks and 128 living states. Time-varying
treatment regimes, coarsened causal histories and unmeasured treatment-
confounding sensitivity are outside this API. See `?rc_causal` and the
baseline-causal vignette for the full input contract and diagnostic definitions.

Verification includes 80 added expectations, an independent analytic EIF
oracle, exact probability-one reduction to the previous targeting kernel, a
1,600-person synthetic example, grid refinement, a dense matrix-exponential
oracle, and a separate million-subject analytic remainder calculation.
See `reviews/causal_revision/` for the current checks. No broad new causal
coverage simulation or clinical benefit is claimed by these checks.
