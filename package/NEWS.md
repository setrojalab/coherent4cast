# recurcast 0.10.1 (unpublished targeted revision)

- Add the S12 baseline binary-treatment extension through rc_causal() and
  rc_causal_contrast(): common-population g-formula standardization, cross-fitted
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

# recurcast 0.10.0 (unpublished release candidate)

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

# recurcast 0.9.1 (unpublished release candidate)

- Assemble the first local release candidate after internal 0.9.0.9000
  development snapshots. This version has not been published to CRAN.
- Preserve weighted-rank scale invariance for finite weights near the
  floating-point limit, with a regression test.
- Retain the generated vignette index in source package builds.

Development implementation of the Version 5 RecurCast computational methods.
The release gates, measured checks and unexecuted platform checks are recorded
outside the installable package. A development version does not indicate CRAN
acceptance or a published software paper.
