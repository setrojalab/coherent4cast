# Release gate evidence

This record describes internal agent review roles, not independent professors
or external certifications. Candidate 0.10.1 adds R interfaces and audit
metadata while preserving every native source/header and Makevars hash from
0.10.0. It is not a CRAN release or an accepted journal submission.

| Gate | Executed evidence or remaining scope |
|---|---|
| G0: provenance and API | Preserved source/oracle records; 40 public exports and 33 registered S3 methods. Earlier package-name searches are historical evidence, not a fresh name reservation. |
| G1: claim traceability | API.md, claims register and the targeted addendum distinguish implemented interfaces, assumptions and deferred extensions. |
| G2: native safety and coverage | Six local C++17 fixtures pass strict warnings. Hosted Linux GCC and Clang each pass strict tests with and without ASAN/UBSAN: 4/4 jobs, run 45164392550. Combined R/C++ line coverage is 95.46236% (2777/2909), no excluded files. |
| G3: numerical and interface contracts | All 917 predecessor assertions were re-executed successfully. The final candidate passes 1107 assertions, including ordinary survival, joint competing CIFs, cross-cause sandwich covariance, future covariate mappings, cause-aware validation and audit provenance. All 23 native files remain byte-identical. |
| G4: checked artifacts | Windows R 4.6.2 and R-devel 4.7.0 source checks both report 0 errors, 0 warnings, 0 notes. Both execute 1107 passing assertions, 0 skips. A separately installed Windows binary passes the same 1107 assertions. |
| G5: benchmark scope | Existing 0.10.0 measurements are retained as predecessor evidence; no new broad speed claim is inferred. The targeted synthetic and standard-regime drivers run against the final binary. |
| G6: remote platforms | All six hosted R-release/devel checks on Linux, macOS and Windows pass with Status: OK and 1107 assertions each. All four GCC/Clang strict and ASAN/UBSAN jobs pass. Tested commit: fbf6429e37aa91e159d878a8228d62803b287056. Private PR 1; logs and artifacts retained. |
| G7: documentation and data | Four vignettes build and check successfully. The 54-page Rd manual is rendered and inspected. Synthetic ADTTE/OMOP-style mappings pass 16 checks; these are analyst-controlled examples, not certified native adapters. |
| G8: replication | Predecessor Windows numerical repeatability records remain explicitly labelled. Current targeted driver outputs and synthetic mapping checks are saved with final binary provenance. Cross-platform bitwise identity is not asserted. |
| G9: CRAN | Not submitted; acceptance is an external outcome. |
| G10: manuscript | Current JSS-style manuscript: 11 tables and 5 figures; separate targeted API/audit addendum. Word/PDF and source are supplied. Internal checks do not replace human author review or external peer review. |
| G11: public release | No public DOI, arXiv, CRAN or JSS release has been made. The GitHub repository is private. |
| G12: versioning and support | Additive 0.10.1 pre-1.0 candidate; legacy trajectory/forecast classes retained, explicit no-frailty alias, deterministic schema tags and structured audit exports. |

Canonical local artifact and coverage evidence:
`reviews/P6_final_verification.json`. Hosted evidence:
`reviews/github_ci/` and `planning/GITHUB_REPOSITORY.json`.

- [Private validation pull request](https://github.com/setrojalab/recurcast/pull/1)
- [Native safety run](https://github.com/setrojalab/recurcast/actions/runs/45164392550)
- [R platform run](https://github.com/setrojalab/recurcast/actions/runs/45164392545)

