# recurcast 0.10.1

Coherent survival modelling and dynamic forecasting of history-dependent multitype recurrent events and terminal risk under incomplete observation.

**Private review candidate.** This repository remains private while review is in progress. Access is limited to authorized collaborators. Version 0.10.1 has not been accepted by CRAN or a journal; the accompanying software article is a draft for review.

`recurcast` fits an event-history probability law and derives type-specific expected recurrence counts, terminal-event risk, survival and restricted mean time alive or operational. The R interface retains observation assumptions, convergence flags and numerical diagnostics. It includes no-frailty standard-survival and competing-event interfaces alongside the recurrent-event model.

## Install from the review branch

An authorized GitHub account and credentials recognized by `remotes` are required for this private repository. For example, `remotes` can use a locally configured `GITHUB_PAT`. Its [official installation documentation](https://remotes.r-lib.org/reference/install_github.html) describes authentication and build options.

```r
install.packages(c("remotes", "Rcpp", "digest", "jsonlite"))
remotes::install_github(
  "setrojalab/recurcast",
  subdir = "package",
  ref = "release/0.10.1-validation",
  dependencies = NA,
  upgrade = "never"
)
library(recurcast)
packageVersion("recurcast")
```

The package is in the `package/` subdirectory. Source installation requires a C++17 compiler toolchain; use the appropriate Rtools for Windows. The review branch can advance, so retain the installed commit and package version when recording an analysis. These instructions do not grant access to the private repository.

## Download the review materials

The [clean review bundle](release/recurcast_0.10.1_review_bundle.zip) groups the supplied package, documentation and replication materials. It excludes internal planning notes; the local full working archive is retained separately. Check downloaded files against [SHA-256 checksums](release/SHA256SUMS_GITHUB.txt). These downloads remain subject to this private repository's access controls.

## Install a supplied archive

Download the appropriate artifact while signed into an authorized account:

- [Source archive](release/recurcast_0.10.1.tar.gz).
- [Windows binary](release/recurcast_0.10.1.zip), built for the recorded Windows R environment.
- [Reference manual](release/recurcast_0.10.1_Reference_Manual.pdf).

After downloading the repository or placing the files in a local `release/` directory, run **one** installation route:

```r
install.packages(c("Rcpp", "digest", "jsonlite"))

# Source route: requires a C++17 toolchain.
install.packages("release/recurcast_0.10.1.tar.gz",
                 repos = NULL, type = "source")

# Alternative Windows route: use a compatible Windows R installation.
# install.packages("release/recurcast_0.10.1.zip",
#                  repos = NULL, type = "win.binary")
```

Optional dependencies support specific examples: `survival` supplies the lung and colon data, `frailtypack` supplies the readmission data, and `nlme` supports estimated Gaussian longitudinal covariate models. Installation and downloading dependencies require access to their sources; fitting and prediction can run offline once the required environment is installed.

## First fit and forecast

The bundled demonstration data are synthetic:

```r
library(recurcast)
data("sim_perinatal", package = "recurcast")

fit <- recurcast(~ x, sim_perinatal,
                 event_types = c("1", "2"), time_unit = "months")
print(fit)

one <- sim_perinatal[sim_perinatal$id == 1, ]
predict(fit, newdata = one, landmark = 0,
        horizon = c(0.25, 0.5, 1), type = "all")
```

Expected recurrence counts are burdens, not occurrence probabilities. Predictions condition on the declared landmark, eligibility and covariate-evolution rule. Numerical error diagnostics are not sampling confidence intervals.

| Analysis task | Entry point | Scope |
|---|---|---|
| Recurrent and terminal histories | `as_trajectory()`, `recurcast()`, `predict()` | Explicit event marks, independent units, time units and observation model. |
| No-frailty model | `recurcast(frailty = "none")` | Fixes latent frailty variance at zero. |
| Standard survival | `rc_surv()` | Full piecewise-exponential proportional-hazards likelihood; not Cox partial likelihood. |
| Mutually exclusive absorbing causes | `rc_competing()`, then `predict()` | Combines all cause hazards into coherent cumulative incidence and survival; no automatic CIF confidence bands. |
| Held-out prediction evaluation | `validate()` | Ten native metrics with declared censoring assumptions; earlier competing events are observed noncases for cause-specific validation. |
| Traceable audit export | `rc_audit_report()` | Captured and export context, numerical summaries and explicit SHA-256 fingerprints; not regulatory certification. |

`rc_trajectory` and `rc_forecast` are additive S3 classes. Use the existing `as_trajectory()` and `predict()` functions; the class names are not new estimator functions. Observed-pattern labels do not infer scientific recurrence, causal structure or the missingness mechanism.

## Worked examples and help

From a downloaded checkout's root, after installing the package and the relevant optional dependencies:

```sh
Rscript replication/standard_regimes.R
Rscript replication/standards_mapping.R
```

- [Standard-regime examples](replication/standard_regimes.R) check a constructed two-cause model against closed-form forecasts and demonstrate a declared training/holdout analysis using `survival::lung`.
- [Standards-mapping examples](replication/standards_mapping.R) use wholly synthetic ADTTE/OMOP-style tables. They expose analysis-parameter selection, censor coding, episode construction, time origins and terminal stopping. They are analyst-controlled examples, not native certified adapters.
- [API and assumptions](package/API.md) map the exported functions to their statistical and computational contracts.

`RECURCAST_LIBRARY` selects a local installed library. `RECURCAST_STANDARD_OUTPUT` and `RECURCAST_MAPPING_OUTPUT` select separate output directories. A GitHub package installation alone does not download the project-level replication directory; use the repository checkout for these scripts.

```r
help(package = "recurcast")
help("rc_surv", package = "recurcast")
help("rc_competing", package = "recurcast")
help("rc_audit_report", package = "recurcast")
vignette(package = "recurcast")
```

The supplied checked archives contain four evaluated HTML vignettes. A direct `remotes::install_github()` build omits vignette building by default. To build them from source, use `build_vignettes = TRUE`, install the vignette dependencies listed in DESCRIPTION and provide Pandoc. The vignette sources are also readable here:

- [Getting started](package/vignettes/getting-started.Rmd).
- [Incomplete observation](package/vignettes/incomplete-observation.Rmd).
- [Validation and calibration](package/vignettes/validation-calibration.Rmd).
- [Extended workflows](package/vignettes/extended-workflows.Rmd).

## Manuscript and targeted addendum

- **Software article:** [PDF](paper/01_RecurCast_JSS_Software_Paper.pdf), [editable Word](paper/01_RecurCast_JSS_Software_Paper.docx).
- **Targeted API and Audit Addendum:** [PDF](paper/07_Targeted_API_And_Audit_Addendum.pdf), [editable Word](paper/07_Targeted_API_And_Audit_Addendum.docx).
- **Inherited 0.10.0 scientific companion — Revision Verification Supplement:** [PDF](paper/05_Revision_Verification_Supplement.pdf), [editable Word](paper/05_Revision_Verification_Supplement.docx). It develops the unchanged mathematical and computational modules.
- **Inherited 0.10.0 scientific companion — Cross-Field Metrics and Usability Guide:** [PDF](paper/06_Cross_Field_Metrics_and_Usability_Guide.pdf), [editable Word](paper/06_Cross_Field_Metrics_and_Usability_Guide.docx). It explains the unchanged forecast, validation and reporting modules. These inherited documents are not newly executed 0.10.1 studies.

The article retains the original method equations and version-labelled numerical evidence. The addendum describes the 0.10.1 wrappers, classes, audit record, mapping examples, worked calculation and package roadmap. Earlier companion documents 02–06, when retained in the archive, describe the unchanged 0.10.0 modules and are labelled as inherited material.

## Validation status

The completed local Windows source checks on R 4.6.2 and R-devel 4.7.0 reported zero errors, warnings or notes. Each source environment and the freshly installed Windows binary passed 1,107 test expectations without failures, warnings or skips. Measured R/C++ line coverage was 95.46236% (2,777/2,909 executable lines), with no file exclusions. All 23 native source/header/Makevars files match the preserved 0.10.0 version.

At tested commit `fbf6429e37aa91e159d878a8228d62803b287056`, all six hosted R-release/R-devel checks on Linux, macOS and Windows completed with `Status: OK`, each passing 1,107 expectations with no failures, warnings or skips ([R workflow](https://github.com/setrojalab/recurcast/actions/runs/45164392545)). All four Ubuntu native jobs passed GCC/Clang strict builds and AddressSanitizer/UndefinedBehaviorSanitizer runs ([native workflow](https://github.com/setrojalab/recurcast/actions/runs/45164392550)).

The execution history retains two earlier macOS R-devel infrastructure interruptions: a missing gettext header during dependency compilation and a stalled repository-status diagnostic before installation. Workflow-only configuration changes resolved both without modifying package science or native code. These results cover the tested configurations and fixtures. See [the release gate record](package/GATES.md) and the private repository's Actions records for scope and logs.

## Statistical interpretation

The package is prognostic. A shared probability law gives internally compatible forecasts under its model; it does not establish that treatment, monitoring or another intervention prevents events. Source-data semantics, risk sets, recording assumptions and model regularity remain part of the analysis.

Finite-path targeting, continuous-time targeting, recording-sensitivity bounds and prediction-set calculations have different assumptions. When invoking the observed-reward finite-sample confidence construction, request `interval = "projection"` throughout. Data-dependent switching between approximate Gaussian and projection intervals does not itself inherit that guarantee. Some sparse or nonregular fits withhold regular intervals; uncapped recurrence burdens may have uninformative outer limits.

Additional reporting functions in `replication/research_metrics.R` are sourced companion code, not package exports. They do not turn predictive accuracy or fixed-cost accounting into a causal benefit claim. The possible future names `recurbounds`, `recurcausal` and `recurmetrics` describe a roadmap; this review candidate remains one package and repository.

## Authors and citation

Sunday A. Adetunji and Oluwatimileyin A. Adetunji. Use `citation("recurcast")` or the [citation file](CITATION.cff) for the current software citation. The [license file](LICENSE.md) records the package's GPL-3 licensing. A private repository, a passing software check or an audit hash is not a claim of public release, regulatory approval, clinical validation or journal acceptance.
