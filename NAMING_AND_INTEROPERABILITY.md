# coherentSurv naming and interface conventions

This convention applies to the 0.11.0 naming revision. It distinguishes repository requirements, language conventions and project choices. Statistical algorithms remain in one maintained R/C++ implementation.

## Names and metadata

| Item | Convention |
|---|---|
| R package and paper identity | `coherentSurv` |
| CRAN title | `Coherent Survival Modelling and Dynamic Landmark Forecasting` (60 characters) |
| Software-paper title | `coherentSurv: Coherent Survival Modelling and Dynamic Landmark Forecasting in R and Python` |
| Package-specific public R operations | Lowercase `cs_` names; e.g. `cs_fit()`, `cs_control()`, `cs_as_trajectory()`, `cs_build_history()` |
| Specialized inference | `cs_target()`, `cs_bounds()`, `cs_causal()`, `cs_causal_contrast()`; retain their distinct estimands and assumptions |
| Standard R model operations | `predict()`, `summary()`, `print()`, `coef()`, `vcov()`, `logLik()` and existing plotting methods |
| Private R helpers | `.cs_` prefix; not exported |
| Native implementation | `coherentSurv` namespace, package-specific registered routines and `CS_` helper/guard names where applicable |
| Python access to R | `engine = importr("coherentSurv")`, followed by `engine.cs_fit(...)` and registered R generics |
| Future Python import module, if distributed | Lowercase `coherentsurv`; this is a language spelling, not a second statistical engine or a currently promised PyPI release |

CRAN requires an alphabetic first character, at least two characters, ASCII letters/numbers/dots, and no final dot for package names. Mixed case is permitted. Names must not conflict, ignoring case, with current or archived CRAN packages or current Bioconductor packages. CRAN generally treats a published name as permanent. The checks on 25 November 2026 found no exact case-insensitive match for the proposed name in those public listings; this does not reserve a name or establish legal clearance. [CRAN policy](https://cran.r-project.org/web/packages/policies.html), [Writing R Extensions](https://cran.r-project.org/doc/manuals/R-exts.html#The-DESCRIPTION-file).

The title is one line, uses title case, excludes the package name, and has no final full stop. R's manual describes possible truncation at 65 characters; it does not describe a universal hard 65-character rejection rule. The `cs_` prefix is a project convention, not a CRAN requirement.

## Familiar model interfaces

Use verbs and lowercase underscore-separated names for package-specific operations, following the [tidyverse function naming guidance](https://style.tidyverse.org/functions.html). Preserve familiar model generics rather than requiring users to learn a replacement for every operation. The established [survival model workflow](https://stat.ethz.ch/R-manual/R-devel/library/survival/html/coxph.html) and [lme4 model methods](https://lme4.github.io/lme4/reference/merMod-class.html) provide relevant examples. These conventions do not imply that the packages fit the same statistical model.

The default workflow is prepare histories, declare the model, fit, inspect diagnostics, forecast at an explicit landmark and horizon, and assess held-out performance. Identifiers, interval starts/stops, recurrent marks, terminal status, observation flags and time units remain explicit. Unrecorded intervals cannot silently become observed zero-event intervals. Forecast probabilities, expected counts and restricted survival time retain their distinct units.

## Python interoperability

Python's [PEP 8](https://peps.python.org/pep-0008/#naming-conventions) recommends lowercase module names, lowercase underscore-separated function names, and CapWords class names. [PyPA](https://packaging.python.org/en/latest/specifications/name-normalization/) normalizes distribution names independently of import-module spelling. Changing capitalization does not create a distinct distribution identity.

Python callers must have R and the R package installed. The included Rscript runner handles process-based access. The rpy2 example converts input data frames locally, retains the fitted object in R, invokes R's prediction method, and converts the result table back to pandas. Conversion can copy data; zero-copy operation and automatic speed gains are not promised. See [rpy2 package access](https://rpy2.github.io/doc/latest/html/robjects_rpackages.html) and [pandas conversion](https://rpy2.github.io/doc/latest/html/pandas.html).

If a higher-level estimator class is added, follow scikit-learn's documented fit/predict conventions, parameter inspection and learned-attribute naming. Passing its applicable estimator checks and preserving subject-level splitting are prerequisites to claiming pipeline compatibility. Merely exposing functions named fit and predict does not establish this compatibility. [scikit-learn developer guide](https://scikit-learn.org/stable/developers/develop.html).

## Usability conventions to adopt

These are interface recommendations, not claims of equivalence to another package.

| Reference implementation | Established convention | Application here |
|---|---|---|
| survival | Formula and data arguments; explicit event and censoring definitions | Keep a clear `cs_fit(formula, data, ...)` entry point and document the one-sided formula used by this model. |
| lme4 | Fitted objects support standard R methods | Preserve `summary`, `predict`, `coef`, `vcov` and `logLik`; use `cs_` for specialized operations. |
| tidyverse | Verbs, consistent argument names and actionable errors | Use lowercase underscore-separated operations; identify the invalid field, offending record and expected format. |
| statsmodels | A fitted result carries estimates, uncertainty and diagnostics | Keep estimates, covariance, convergence and assumptions together; offer explicit extraction rather than requiring console parsing. |
| scikit-learn | Consistent fit/predict interfaces and tested estimator contracts | Follow these conventions if a Python estimator wrapper is developed; do not claim pipeline compatibility for the direct rpy2 example. |
| scikit-survival | Time-specific, censoring-aware evaluation | State the landmark, horizon, prediction target, eligible population and censoring assumptions for every reported metric. |

Additional primary documentation: [statsmodels workflow](https://www.statsmodels.org/stable/gettingstarted.html), [scikit-survival evaluation](https://scikit-survival.readthedocs.io/en/stable/user_guide/evaluating-survival-models.html).

The supported workflow should have one short executable example for each data regime, followed by task-specific articles for forecasting, sensitivity bounds and causal contrasts. Each example must state the time unit and explain returned column names. Missing history, right censoring, delayed entry and a terminal event must retain distinct meanings.

For evaluation, record the unit of splitting and keep each person's longitudinal records together. Report the outcome definition, horizon, sample size, event count, estimate, uncertainty procedure and conditions under which a metric is undefined. A rank score is not a probability; an expected recurrence count is not a first-event risk. Match each metric to the output it evaluates.

Python documentation should show installation prerequisites, a pandas input example, data-type conversion, a forecast table, a numerical parity check and the tested platform. Data conversion should occur once per batch where possible. The direction of the language bridge does not establish a speed advantage; benchmark fitting, transfer and prediction separately.

## Release acceptance checks

- All current source, help, examples, native registration, website instructions and rendered manuscript code must use the current names.
- Regenerate Rcpp registration and the R namespace, then run package, native and platform checks.
- Compare coefficients, covariance matrices and forecasts with the previous numerical implementation after name normalization.
- Run Python-to-R parity checks separately and record the platforms and dependency versions actually exercised.
- Preserve numerical diagnostics, undefined results and observation assumptions; a completed fit is not evidence of clinical validity.
- Keep documentation examples executable, public downloads checksummed, and submission status explicit.
- This naming revision changes the public API. Existing scripts need the current function names, and previously serialized model objects should be refitted. Historical commits and published preprint versions remain historical records.
