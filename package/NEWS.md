# coherentSurv 0.11.2

* Update repository, website, issue-tracker and citation links after the
  organization rename to SetrojaLabs (setrojalabs).
* No changes to statistical calculations, public interfaces or numerical kernels.

# coherentSurv 0.11.1

- Reject unsupported validation arguments; report undefined calibration slopes
  under complete or quasi-complete separation.
- Withhold regular Wald inference for nonregular fits, including competing-risk
  components, through the supported reporting interfaces.
- Add validated, explicitly scoped ADaM/OMOP mappings and offline environment
  diagnostics; these are not regulatory qualification or certified parsers.
- Make the Python runner fit once by default; optional parity verification and
  stage timings are explicit. Preserve leading zeros and literal NA identifiers.
- Expose code-specific error classes and a nonregular-fit warning class while
  preserving the broad coherentSurv condition classes and existing messages.
- Preserve declared estimands in targeting audits and expose declaration and
  source-mapping metadata in audit reports.
- Correct companion fixed-cost valuation checks and citation/version metadata.
- Preserve the existing R/C++ numerical algorithms and archived study results.

# coherentSurv 0.11.0

- Adopt coherentSurv as the package identity and cs_ for package-specific functions, classes, attributes and native helpers.
- Use cs_fit(), cs_control(), cs_as_trajectory() and cs_build_history() as the main workflow; preserve standard R generics.
- This release changes API names. Refit previously serialized model objects with the current package.
- Statistical algorithms are unchanged by the naming migration. Current validation is recorded separately from earlier results.

## Earlier development history

# coherentSurv 0.10.1 (unpublished targeted revision)

- Add the S12 baseline binary-treatment extension through cs_causal() and
  cs_causal_contrast(): common-population g-formula standardization, cross-fitted
  arm-specific intensity targeting, treatment/censor weighting, joint EIF
  covariance, differences and ratios, and explicit overlap/boundary diagnostics.
- Add a separate causal C++17 engine; retain the existing prognostic kernels.
- Add a causal vignette, analytic oracle checks and treatment-effect examples.
- Update maintainer contact to sundayadetunjisa@gmail.com.

- Add explicit no-frailty selection and conventional survival/competing-event
  entry points through R-level routing to the existing numerical kernels.
- Add trajectory field roles and a shareable audit report with declared
  serialization, software versions and SHA256 provenance.
- Preserve all 0.10.0 native sources; record new verification separately.
- Record causal revision checks separately from the earlier platform and
  sanitizer results; previous runs do not validate newly added code.

# coherentSurv 0.10.0 (unpublished release candidate)

- Add automatic coarsened-history fibres and posterior mixtures at forecast origins.
- Add sparse state propagation and structural finite-support bound calculations.
- Add continuous-time CAR targeting with explicit boundary diagnostics and a
  separately declared bounded observed-reward projection that retains censoring uncertainty.
  Its unconditional finite-sample guarantee requires using that construction
  throughout; automatic Gaussian/projection switching does not inherit it.
- Add lognormal and estimated discrete frailty likelihoods alongside shared gamma.
- Add declared scheduled and Gaussian longitudinal future-covariate trajectories.
- Keep the 0.9.1 project and its recorded evidence unchanged; verification of
  this revision is recorded separately.

# coherentSurv 0.9.1 (unpublished release candidate)

- Assemble the first local release candidate after internal 0.9.0.9000
  development snapshots. This version has not been published to CRAN.
- Preserve weighted-rank scale invariance for finite weights near the
  floating-point limit, with a regression test.
- Retain the generated vignette index in source package builds.

Development implementation of the Version 5 coherentSurv computational methods.
The release gates, measured checks and unexecuted platform checks are recorded
outside the installable package. A development version does not indicate CRAN
acceptance or a published software paper.
