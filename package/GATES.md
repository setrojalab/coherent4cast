# Release gate evidence

This record describes internal review lenses, not independent professors or
external certifications. A check passes only when the recorded execution
supports the stated scope. External acceptance and service-level targets cannot
be supplied by local tests. Targeted candidate 0.10.1 is under verification. The numerical check results below describe the preserved 0.10.0 baseline until the final candidate evidence replaces them. Earlier
0.9.1 verification is retained in the separate baseline project.

| Gate | Required evidence | Current state | Internal lenses |
|---|---|---|---|
| G0 | Frozen oracle, license audit, name checks, API | Oracle frozen; CRAN and archive name checks clear, PyPI returns 404; Bioconductor request blocked | P6 and P5 |
| G1 | Every public claim and theorem-to-code fixture traced | API.md and planning/claims_register.csv trace the implemented scope, measured results and deferred claims; internal review evidence retained | P2 and P1 |
| G2 | Warning-free core and sanitizers | Six standalone C++17 test executables compile with strict warning flags and pass; combined R/C++ line coverage is 95.28487% without file exclusions. Linux sanitizer execution remains pending | P6 and P1 |
| G3 | Deterministic oracle parity and numerical guards | Implemented-module fixtures pass, including gamma scores, coarsened likelihood, targeting EIF, LP/projection, regularization and extreme-weight conformal ranks | P1 and P3 |
| G4 | Current R-devel and release tarball checks | Final 0.10.0 source checks on Windows R 4.6.2 and R-devel 4.7.0: 0 errors, 0 warnings, 0 notes and 917 passed expectations each. Fresh Windows binary installation also passes 917 expectations | P6 and P5 |
| G5 | Frozen measured benchmarks | Four sample sizes and five repeats executed; raw results, hardware, hashes and unsuccessful quadrature/early convergence checks retained. General high-level 250000-unit fitting remains unmeasured | P3 and P6 |
| G6 | Remote CI, security settings, platform matrix | Workflow files prepared; no remote workflow run or repository setting claimed | P6 and P3 |
| G7 | Datasets, vignettes and scope statements | Three synthetic datasets, two canonical source-package adapters, five case studies and four evaluated vignettes completed; 46-page Rd reference manual rendered and checked | P2 and P4 |
| G8 | Paper replication and cross-platform determinism | Two fresh Windows case-study fits give identical scientific values in all 15 CSV tables. All original 11 manuscript files, and all 16 current scientific manuscript files, agree between regenerated runs and the current paper materialization. Second OS pending | P3 and P1 |
| G9 | CRAN acceptance | Not submitted; acceptance is an external outcome | P6 and P5 |
| G10 | JSS-style manuscript and six-lens internal review | Official JSS-style PDF, editable Word, source, seven tables, four figures and six internal role reports completed; author approval and external peer review are separate | P5 and P3 |
| G11 | Public DOI and arXiv posting | Not published; metadata and release instructions are prepared locally | P5 and P4 |
| G12 | Versioning, support and reserved extension behavior | Policies prepared; behavior checked with API tests | P6 and P2 |

Evidence paths in the project-level reviews directory contain command outputs
and machine-readable summaries. Each lens should record its own review and an
evidence-based second-lens check rather than a fabricated signature.



