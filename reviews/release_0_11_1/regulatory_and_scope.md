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
