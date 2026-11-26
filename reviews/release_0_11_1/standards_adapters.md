# Selected ADaM and OMOP adapter scope audit

Date: 2026-11-26. This is an internal implementation audit, without external
institutional, regulatory or journal endorsement.

Added `package/R/standards.R`, two reference pages and
`package/tests/testthat/test-standards-adapters.R`. Existing statistical and
native kernels are unchanged. Namespace registration belongs to the parent
release integration; this subtask did not edit `NAMESPACE` or its generator.

`cs_from_adam()` selects one explicit PARAMCD, retains positive numeric AVAL
without date rederivation, and requires a duration definition and time unit.
Explicit event/censor code sets must be disjoint and cover every selected CNSR
value; no `1 - CNSR` conversion is used. Each selected subject has one row.
Optional ADSL covariates join by unique exact subject identifiers; missing values
are retained. When present, AVALU must match the declared unit, ignoring case.

`cs_from_omop()` consumes a reviewed selected event extract and a caller-supplied
concept map. It requires one continuous observation period per person and a
separately supplied process origin equal to its start. The explicit
`exact_day_boundaries` convention uses Date differences without inclusive-day
adjustment. All events must precede death or administrative exit strictly.
Unknown concepts or persons, duplicate episodes, same-day multiple events,
terminal/exit ties, out-of-period records and post-terminal events fail visibly.
Multiple periods, delayed entry, unknown hidden-event counts, date imputation
and intermittent observation are outside the adapter. Valid OMOP records can
therefore be outside this narrower mapping contract. Nothing is jittered,
deduplicated or silently dropped.

Both functions require ascertainment and provenance text, retain `cs_mapping`
metadata and return the existing validated `cs_trajectory` representation.
Source data are not modified. Complete ascertainment, origin choice, clinical
coding, timing and censoring assumptions remain analyst declarations; table
presence and schema checks cannot prove them. Character keys preserve leading
zeros and large identifiers. Numeric identifiers cannot recover precision lost
before import. These are selected-extract mappings, not full/certified CDISC or
OMOP parsers, causal identification methods or regulatory qualification tools.

Validation executed with Windows R 4.4.2 UCRT (`x86_64-w64-mingw32`): nine
testthat blocks, 68 passing assertions, zero failures/errors/warnings. Checks
cover leading-zero IDs, parameter selection, multi-code CNSR, ADSL joins,
missing/unmapped values, temporal containment, overlap/duplicate rejection,
terminal ties, time-unit conversion, zero-event extracts and serialized input
immutability and mapping retention through schema revalidation. Both Rd files
passed parsing/checking and their examples executed.
The checks sourced `schema.R` and `standards.R` directly; installed-package
export/build validation is the parent release's separate responsibility.

Additional integration smoke checks loaded the existing 0.11.0 native library
from `build/coherentSurv.Rcheck` with the candidate `schema.R`, `fit.R`,
`standards.R` and `audit.R` sources. Both ADaM and OMOP synthetic trajectories
were fitted with fixed zero frailty. Their full mapping attributes persisted in
`fit$data`; fit audit reports retained rule declarations and omitted the ADaM
source-row vector and OMOP subject-origin table. This is targeted integration
evidence, not a replacement for rebuilding the complete 0.11.1 candidate.

Reproduction from the repository root in R (unset an inherited LC_ALL only if
the local testthat/withr installation reports language-setting warnings):

```r
source("package/R/schema.R")
source("package/R/standards.R")
testthat::test_file("package/tests/testthat/test-standards-adapters.R",
                    reporter = "summary", stop_on_failure = TRUE)
```

Source contract references: [OMOP CDM v5.4 table documentation](https://ohdsi.github.io/CommonDataModel/cdm54.html)
and [CDISC ADaM time-to-event training material](https://www.cdisc.org/sites/default/files/2024-04/2024_CDISC-Session%207B%20ADaM_Mitchikou%20Tseng.pdf).
The supported contracts above, rather than the full breadth of these standards,
define this implementation.
