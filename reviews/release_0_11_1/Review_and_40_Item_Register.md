# coherentSurv 0.11.1: review and requirements register

26 November 2026. Review candidate; published 0.11.0 remains unchanged.

This is a five-discipline internal technical review, not an endorsement by named JSS professors or regulatory agencies. Authorship remains Sunday A. Adetunji and Oluwatimileyin A. Adetunji. No score substitutes for executed checks or external peer review.

## Completed work

- Statistics: guarded summary, tidy and confint methods; nonregular competing-risk components invalidate regular coefficient intervals; separated calibration slopes remain undefined; unused validation arguments are rejected.
- Engineering: Python fits once by default, retains the R model, preserves identifier strings, and exposes optional parity verification and stage timings. Offline environment diagnostics compare four native calculations with analytic references.
- Economics: companion valuation now rejects incomplete component sets, inconsistent time units and invalid cost maps, including factor-name indexing errors. All 50 companion reporting checks passed.
- Medicine: scoped ADaM and OMOP mappings require explicit endpoint, censoring, timing and ascertainment rules. They do not infer complete care histories or invent within-day event ordering.
- Policy and governance: audit metadata preserves scientific declarations and mapping rules while omitting subject-origin/source-row mappings. Typed conditions aid workflow logging; hashes and diagnostic reports are not regulatory certification.

## Current verification

Local package check: 0 errors, 0 warnings, 1 NOTE (new submission/maintainer only), with remote incoming checks enabled.

Executed CI: seven platform jobs (Linux/macOS/Windows release and devel, plus Linux oldrel), four GCC/Clang native configurations including ASAN/UBSAN, and Linux rpy2 all passed. The observed expectation count is 1476. Source revisions and separate logs govern which changes each check verifies. The initial local snapshot exposed RNG-state handling and namespace qualification issues; the failure log is preserved and the corrections are tested in the final revision.

The Linux rpy2 example returned 15 forecast rows with maximum absolute differences of zero for coefficients, covariance and forecasts versus direct R. Identifier and invalid-input checks passed. On its 160-subject synthetic fixture, the recorded initial candidate run took about 0.248 seconds for fitting, 0.0064 seconds for input conversion and 0.0193 seconds for prediction. These are one-run observations excluding interpreter startup, not comparative speed benchmarks.

All 23 native C++/header files included in the source archive agree with 0.11.0 after line-ending normalization. No numerical kernel was replaced by these interface additions.

## New cross-field tutorials

Ten executable vignettes accompany the candidate, including new clinical-trial, fleet-reliability, credit-risk, cross-discipline and source-mapping/environment tutorials. Local generators reuse cs_simulate with explicit parameters. Coarsened timestamps retain their declared observation experiment; unknown episode counts are not silently reconstructed. The tutorials execute actual package calls rather than disabling incomplete examples.

Restricted time alive is not event-free time or QALYs. Recurrent event counts are not monetary loss without a supplied valuation model. Absorption is a scientific state definition; administrative loss to follow-up is censoring. Synthetic tutorials demonstrate computation and do not establish clinical effectiveness, fleet safety or credit-policy validity.

## Manuscript citations and naming

The revised JSS/arXiv LaTeX adds one supported statement from each requested article and retains the existing bibliography. The oncology article motivates complementary Cox, 365-day restricted mean survival and hierarchical Bayesian survival summaries. The survey article motivates study-design-aware validation, including whole-PSU splitting and fold-contained preprocessing. Neither is presented as external validation of coherentSurv.

The survey paper's Crossref record and running citation identify 2026, while other publisher fields show 2027. The bibliography uses the DOI-registered 2026 citation; the discrepancy is recorded for author proof checking. The software remains S3/Rcpp/C++17 with cs_ package-specific names and standard R generics. coherentSurv satisfies conventional R naming syntax; no naming rule certifies academic or regulatory quality.

References added: https://doi.org/10.33545/gynae.2026.v10.i4ah.2816 and https://doi.org/10.29244/currbiomed.5.1.7.

## Interpretation of the 40-item blueprint

Every requested item is retained below. Existing capability, a narrower supported implementation and new research work are distinguished. Delayed-entry frailty selection, longitudinal treatment policies, full joint endogenous-biomarker estimation, spline baselines, posterior-centered adaptive quadrature, OpenMP and a native Python wheel are not claimed as implemented. The native wheel also changes the earlier one-engine maintenance decision. None of these is achieved by relabelling a function or tightening a numerical tolerance.

CRAN submission remains paused. The author has authorized the arXiv revision after passing package and GitHub checks; submission status is recorded separately from this technical review. Draft review: https://github.com/setrojalab/coherentSurv/pull/4.


### I.1 Uncertainty separation

Retain distinct numerical error fields, sampling uncertainty, conformal/projection sets and identification bounds. Existing API/help describes their different interpretations. New confidence methods must preserve regularity gating; no numerical tolerance is described as a confidence interval. A comprehensive typed result redesign is deferred.


### I.2 Canonical coarsening

Existing fibre code and tests are the relevant evidence (`coarsen_fibre.R`, `fibre.hpp`, coarsening tests). Do not change enumeration during an interface patch or claim arbitrary overlap cases are proven by a few fixtures. No silent jitter is introduced.


### I.3 One complete-history law

Existing forecasts retain one fitted law. `cs_validate_env()` now checks the analytic survival/absorption pair and integrated survival for a two-state example through dense and sparse native paths. That fixture verifies its stated case, not every forecast regime.


### I.4 Louis information

Existing `cs_coarsen()` and mixture kernels implement missing-information corrections. The new diagnostic independently checks a two-path mixture: posterior weights 0.2/0.8, score 2.6, observed information 4.96, rather than complete-information average 5.6.


### I.5 Regularity gating

Preserve gradient, information and boundary checks already in fitting. Candidate confidence methods and typed nonregular warnings are handled by the statistical integration work; an optimizer zero code remains insufficient.


### I.6 Frailty boundary diagnostics

Existing family diagnostics retain their assumptions and descriptive comparisons. No new likelihood-ratio calibration or distributional guarantee is introduced.


### I.7 Temporal predictability

Existing history, forecast and coarsened-prediction tests enforce the declared landmark contracts. Preserve explicit scenario versus observed-future distinctions; this patch adds no temporal shortcuts.


### I.8 Numerical oracles

Keep independent analytic/dense fixtures and native tests. A universal 1e-12 threshold for every problem is unjustified without conditioning and approximation analysis. The new tiny diagnostic uses a fixed absolute 1e-12 threshold only for its four declared fixtures.


### I.9 Delayed entry

Deferred scientific extension. The current explicit rejection avoids pretending that unobserved pre-entry frailty selection has been modeled. Implementing conditional likelihoods requires a separate specification, derivative checks and simulations.


### I.10 Nonselective intervals

Preserve explicit projection/Gaussian options and the documented lack of unconditional coverage for an automatic data-dependent switch. Changing existing defaults is outside this nonbreaking diagnostic patch.


### II.1 Cryptographic manifests

Existing audit hashes are integrity fingerprints. Neither hashes nor JSON make an artifact immutable or Part 11 compliant. Compiler provenance is recorded only when actually available; the new diagnostic marks the package build compiler/flags unavailable. Signing, trusted timestamps and controlled storage require a separate deployment design.


### II.2 Estimand declaration

Retain optional estimand metadata and candidate audit-export improvements. Making new study declarations mandatory would break existing calls and cannot be inferred from event-column names.


### II.3 Environment qualification

Implemented narrowly as `cs_validate_env()` in `package/R/environment.R`: offline installation/session metadata, native registration, current native-file SHA-256, four bounded analytic checks, optional tiny timing/object-size observations, and list/JSON/text export. It is explicitly unsigned and mutable, with passed/failed/unavailable statuses. It is not IQ/OQ/PQ certification or a signed qualification report.


### II.4 Offline execution

The new diagnostic performs no network calls, installation, telemetry or external compiler commands. A blanket offline source-install guarantee additionally requires all dependencies and toolchains to be supplied locally. This review does not certify every external dependency's behavior.


### II.5 ADaM/OMOP parsers

Separately implemented candidate adapters and focused tests are being integrated in `standards.R`. Their declared column, censoring, episode and ascertainment contracts govern use. They do not certify arbitrary CDISC/OMOP datasets or replace study-specific adjudication.


### II.6 Locked artifacts

Existing fitted objects retain model specification and controls. They remain mutable R objects. Signed model bundles, validation against a locked runtime and governed SAP/SaMD deployment are deferred; serialization is not immutability.


### II.7 Conditions

Candidate integration adds code-specific condition subclasses and typed nonregular warnings while retaining existing base inheritance. The new diagnostic reports fixture failures without misclassifying them as passes. A rewrite of every warning/error is outside this patch.


### III.1 R-agnostic headers

Preserve `inst/include/libcoherentSurv` and thin Rcpp bridges. Native standalone tests support this separation. No R header or Python-specific dependency is added to the core.


### III.2 Native Python wheel

Deferred under the existing maintenance preference for one R/C++ engine accessed from Python. A wheel would require a separately maintained API, schema/preprocessing contract, packaging matrix and parity evidence. It does not automatically confer scikit-learn/lifelines compatibility or measured superiority.


### III.3 Dependencies

Keep `digest` and `jsonlite` as required audit dependencies rather than silently making current exports conditional. ggplot2/testthat remain optional. No tidyverse runtime dependency is introduced. Moving dependencies is a compatibility decision, not an automatic consequence of the attachment.


### III.4 Cross-platform determinism

Require seeded reproducibility and numerical tolerances supported by tests. Bit-for-bit equality across libm, BLAS, architectures and compilers is not established. No parallel reduction or fast-math change is made.


### III.5 CSR budgets

Existing sparse propagation exposes term/chunk budgets and workspace estimates. Preserve those bounds and fail explicitly when exhausted. The new diagnostic uses two states, at most 100 terms and 10 chunks for its sparse example; it is not a large-state stress test.


### III.6 Zero-copy views

Deferred pending profiling and ownership design. `src/core_bridge.cpp` currently copies vectors and rearranges matrices. Less than 1 GB for arbitrary 250,000-subject designs cannot be promised independently of dimension, history length and returned information matrices.


### III.7 OpenMP

Deferred. No OpenMP implementation or thread-safety guarantee is claimed. The Python example documents independent processes for independent concurrent R sessions. Adding thread parameters without a thread-safe engine would be misleading.


### III.8 Numerical guards

Preserve current finite-value, log-mixture, uniformization and gamma-limit guards and their tests. New environment fixtures exercise selected analytic cases; they are not exhaustive overflow/underflow analysis.


### IV.1 Longitudinal causal policies

Deferred scientific extension. Keep static binary baseline causal scope and explicit assumptions. A new sequential treatment estimator requires its own identification, nuisance-learning and inference validation.


### IV.2 Spline baselines

Deferred likelihood/model extension. Retain the current explicit piecewise-constant baseline specification. Spline basis and penalty choices cannot be introduced as a cosmetic API patch.


### IV.3 Full joint estimation

Deferred. The current conditional Gaussian covariate submodel is not a jointly estimated endogenous biomarker-event likelihood. No comparison claiming equivalence to full joint-model packages is warranted.


### IV.4 Adaptive quadrature

Existing node-refinement diagnostics remain. Posterior-mode-centered adaptive Gauss-Hermite integration is a distinct numerical change requiring derivative and convergence validation; it is not inferred from node doubling.


### IV.5 Nuisance learner plugins

Deferred. Preserve subject-level cross-fitting and the present sufficient-state assumptions. Generic learner integration needs prediction contracts, leakage checks, positivity handling and statistical validation.


### V.1 Formula DSL

Retain explicit `cs_fit()` schema/formula arguments and current standard-regime wrappers. A richer DSL is deferred because parsing alone does not define frailty or observation semantics.


### V.2 S3 contract

Existing methods are retained; candidate `confint()` methods cover four fitted classes. `update`, `residuals` and a complete new model-matrix contract are not claimed. Each additional generic needs a statistically meaningful definition.


### V.3 Table ecosystems

Existing `tidy`, `glance`, `augment` methods remain. Full compatibility with external table packages requires their integration tests; function names alone do not establish it.


### V.4 Plots

Preserve explicit target/uncertainty labels. A diagnostic plotting redesign is deferred; the candidate does not add graphical guarantees or inferred confidence bands.


### V.5 Validation metrics

Existing validation tests and undefined-result flags remain authoritative. The environment diagnostic does not use a predictive-validation dataset, and passing its oracles is not predictive performance evidence.


### VI.1 Coverage and fuzzing

Retain existing property and boundary tests. No new greater-than-95-percent branch coverage claim is made. The focused environment suite checks analytic output, failure precedence, serialization, RNG preservation and optional observations. Large randomized graph fuzzing is deferred.


### VI.2 Sanitizers

Existing CI runs GCC/Clang strict native builds with ASan/UBSan configurations. Do not relabel these as executed Valgrind or rchk runs. Those instruments need separate jobs and evidence.


### VI.3 Platform matrix

Candidate CI adds Ubuntu R-oldrel to release/devel on Linux, macOS and Windows (seven jobs). Linux rpy2 is separately checked. Additional hardware architectures are not inferred from these runs.


### VI.4 Documentation

Ten executable vignettes now include clinical-trial, fleet, credit, cross-discipline and source-mapping/environment workflows. Native Python wheels and OpenMP are not documented as implemented.


### VI.5 Governance/licensing

Preserve GPL-3 and transparent versioned validation. Do not unilaterally promise a new dual license, permanent succession or a retrospectively applied deprecation policy. Document compatibility decisions and future governance separately.


# Deployment and primary sources

# Intended use, evidence and deployment responsibilities

Review date: 26 November 2026. This is an internal technical review, not an
assessment by named journal reviewers, regulators or clinical institutions.
The authorship of the software article is unchanged.

## Research software and regulatory assurance

coherentSurv provides statistical estimation and forecasting. Passing package
checks, reproducing numerical oracles and recording hashes support reproducible
research. They do not establish clinical effectiveness, prospective safety,
regulatory clearance or certification of a hospital deployment.

| Framework | Evidence supplied by the package | Additional application-level work |
|---|---|---|
| FDA medical-device / clinical-decision-support scope | Explicit model scope, reproducible source, numerical diagnostics | Determine the intended use and applicable device/CDS category; validate the complete product and workflow, including users and consequences. |
| FDA electronic records / 21 CFR Part 11 | Data/source integrity fingerprints and captured/exported contexts | Access control, attributable protected records, controlled change history, retention and applicable electronic signatures. A JSON hash is not that system. |
| ICH E9(R1) | Separate prognostic and baseline-treatment causal analyses; declared estimand metadata | Prespecify population, treatment conditions, endpoint, intercurrent-event strategy and population summary; justify identification and sensitivity analyses. Metadata cannot establish these assumptions. |
| NIH rigor and sharing | Runnable code, versioned examples, explicit assumptions and tests | Study-specific rigorous design, applicable award terms, consent/privacy and data-management/sharing plan. No generic NIH software certificate is claimed. |
| NHS DCB0129 / DCB0160 | Diagnostic provenance that can inform a safety case | Manufacturer and deploying-organisation clinical risk management, responsible roles, hazard log, safety case and local validation. |
| NHS DTAC | Source, dependency and interface documentation | Product-specific clinical safety, data protection, technical security, interoperability, usability and accessibility assessment. |
| JSS / CRAN | Source, reproducible examples, manuscript, reference manual and checks | Editorial/reviewer assessment and incoming checks; package acceptance is not regulatory approval. |

Environment diagnostics explicitly report unavailable build information. They
are unsigned observations from the current installation. Timings and object sizes
of a small deterministic fixture are neither a performance qualification nor a
bound on peak memory for arbitrary patient cohorts.

Offline modelling is supported once the R package and declared dependencies
have been installed. Air-gapped installation additionally requires compatible R,
toolchain, dependency archives and (where rebuilding vignettes) document tools.
An ordinary install.packages call may access a repository; no universal
zero-network installation claim is made.

## Naming, interface and maintenance

`coherentSurv` is a valid conventional R package name; mixed case is permitted.
Its DESCRIPTION title has 60 characters. A cs_ prefix is a project convention,
not a CRAN requirement. Standard R generics keep their standard names. The
implementation uses S3 and Rcpp/C++17, not S4 or RcppEigen. No OpenMP or zero-copy
claim is made. `rpy2` and Rscript call the same installed R package.

Public package checks and the 40-item register distinguish implemented,
partially implemented and proposed work. The blueprint's native Python wheel,
new longitudinal causal policies, joint endogenous-biomarker likelihood and
delayed-entry frailty likelihood are new development projects; they cannot be
represented as present simply by adding wrapper names or documentation.

## Clinical and economic interpretation

The two new manuscript citations illustrate complementary survival summaries
and validation that respects study design. They do not validate this package's
clinical decisions or establish treatment benefits. The existing small held-out
death-risk examples are not external prospective validation of all burden and
restricted-survival targets.

The companion decision curve and supplied-cost loss can express different
preferences. The false-positive/false-negative costs C_FP and C_FN imply the
threshold C_FP/(C_FP+C_FN) when their sum is positive. A varying decision curve
varies that implicit cost ratio. Additive burden valuation requires distinct,
non-overlapping components, a perspective, price source/year, horizon and any
discount rule. Fixed supplied prices do not estimate savings, QALYs or causal
welfare effects. Correctly summing rows does not establish economic non-overlap.

## Primary sources checked

- R Core, Writing R Extensions: https://cran.r-project.org/doc/manuals/R-exts.html#The-DESCRIPTION-file
- CRAN Repository Policy: https://cran.r-project.org/web/packages/policies.html
- JSS submissions: https://www.jstatsoft.org/about/submissions
- rpy2 conversion: https://rpy2.github.io/doc/latest/html/pandas.html
- FDA CDS guidance (29 January 2026): https://www.fda.gov/regulatory-information/search-fda-guidance-documents/clinical-decision-support-software
- FDA electronic records guidance: https://www.fda.gov/media/166215/download
- ICH E9(R1): https://database.ich.org/sites/default/files/E9-R1_Step4_Guideline_2019_1203.pdf
- NIH rigor: https://grants.nih.gov/policy-and-compliance/policy-topics/reproducibility
- NIH sharing plans: https://grants.nih.gov/policy-and-compliance/policy-topics/sharing-policies/dms/writing-dms-plan
- NHS DTAC: https://digital.nhs.uk/services/digital-technology-assessment-criteria-dtac
- NHS clinical safety: https://digital.nhs.uk/services/clinical-safety
- NICE economic evaluation: https://www.nice.org.uk/process/pmg36/chapter/economic-evaluation-2/

The public-source review supplies interpretation, not legal or institutional
approval for a particular deployment. Regulatory documents and local obligations
must be checked for the intended product and use when that deployment is defined.
