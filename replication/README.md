# Reproduce the software examples

Run these commands from the repository root after installing `recurcast` 0.10.1.
The scripts use the installed package and the files in this repository. They
do not require the original development workspace or the V5 reference engine.
An optional `RECURCAST_LIBRARY` environment variable selects a nondefault R
library; otherwise the ordinary R library search path is used.

Install the required dependencies in R, then install the supplied source archive
or run `R CMD INSTALL package` with a C++17 toolchain. The current causal source archive is the portable route for this release.
Historical Windows ZIP files elsewhere in the repository predate the causal extension.

```r
install.packages(c("Rcpp", "digest", "jsonlite", "survival",
                   "frailtypack", "nlme", "ggplot2"))
install.packages("release/causal/recurcast_0.10.1.tar.gz", repos = NULL,
                 type = "source")
```

The compact examples and companion checks can then be run as follows:

```sh
Rscript --vanilla replication/standard_regimes.R
Rscript --vanilla replication/standards_mapping.R
Rscript --vanilla replication/check_research_metrics.R
Rscript --vanilla replication/cross_field_usability.R
Rscript --vanilla replication/plot_cross_field_usability.R
```

| File | Purpose and output |
|---|---|
| `standard_regimes.R` | No-frailty survival and competing-risk examples, closed-form probability checks, and a held-out `survival::lung` illustration. Outputs: `results/standard_regimes/`. |
| `standards_mapping.R` | Wholly synthetic ADTTE-style and OMOP-style mappings, with explicit time, censoring, episode and observation rules. Outputs: `results/standards_mapping/`. These are analyst-controlled examples, not certified adapters. |
| `research_metrics.R` | Sourced companion functions for weighted discrimination, threshold metrics, supplied costs and bootstrap summaries. These are separate from native package exports. |
| `check_research_metrics.R` | Executable checks for the companion formulas, denominators, ties, weights and support diagnostics. |
| `cross_field_usability.R` | One synthetic held-out cohort represented under six naming conventions. This is a naming/metric demonstration, not six independent validation studies. Outputs: `results/cross_field_usability/`. |
| `plot_cross_field_usability.R` | Plots the preceding script's saved `companion_report.rds`; run the generating script first. |
| `run_case_studies.R` | Five case-study routes using declared synthetic/public-package data. Source packages supply public records; no private patient file is required. |
| `reproduce_cases.R` | A named fresh run of the case-study driver that refuses to overwrite an existing run and disables fitted-object reuse. |

For fresh fitting of all case studies:

```sh
Rscript --vanilla replication/reproduce_cases.R reviewer_run_1
```

Aggregate outputs go to `replication/fresh_runs/reviewer_run_1/cases/` and
local fitted-object caches go to `build/fresh_case_fits/reviewer_run_1/`.
The alternative direct driver uses `replication/results/cases/` by default.
These full fits can take longer than the compact examples. Public-source
records may have dataset-specific observation limits; the drivers' diagnostics
and interpretation files form part of the results.

Use `RECURCAST_STANDARD_OUTPUT`, `RECURCAST_MAPPING_OUTPUT`,
`RECURCAST_CROSS_FIELD_OUTPUT`, `RECURCAST_CASE_OUTPUT` and
`RECURCAST_CASE_CACHE` to choose other output directories. Run the plotting
script with the same cross-field output setting as its generating script.
Keep `RECURCAST_REUSE_CASE_FITS` unset or `0` for fresh fitting.

Saved numerical snapshots are examples, not universal performance guarantees.
Audit timestamps, source locations, session information, runtime and graphics
device metadata can differ between machines. Numerical diagnostics do not
establish causal identification, prospective clinical benefit or regulatory
compliance. The older `README_revision.md` and V5-oracle Python scripts belong
to the historical method-validation archive; they are not prerequisites for
this eight-script software-reproduction bundle.


## Complete article reproduction and causal extension

Download and extract [RecurCast_JSS_replication.zip](../release/causal/RecurCast_JSS_replication.zip), then run `Rscript --vanilla code.R` from the extracted `recurcast-replication` directory. This reconstructs the article tables and figures, runs the synthetic causal example and verifies the analytic remainder calculation. See its README.txt for fresh fits and longer experiments.

The causal runner and independent oracle are also available under `replication/causal/`. Their reported numerical results concern the explicitly specified synthetic experiments and causal assumptions.
