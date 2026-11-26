# Support

Report reproducible package problems to the maintainer listed in DESCRIPTION.
Include sessionInfo(), the package version, a synthetic example, the full error
message and the expected result. Do not attach patient-level data, credentials
or confidential institutional information.

Maintainer targets are an initial response within seven days and initial
security triage within 72 hours. These are service targets, not measured service
levels. Reproducible issues can be reported at
https://github.com/setrojalab/coherentSurv/issues.

The software supplies prognostic calculations and a separate baseline
binary-treatment causal module under explicit assumptions. It does not supply
clinical decision approval or regulatory certification. The causal module
requires consistency, exchangeability, positivity and the stated censoring
conditions; it is not a general time-varying treatment-policy implementation.
Python examples call the installed R package through rpy2 or Rscript. R and
the package remain required; there is no independent Python statistical engine.
OpenMP and a separate native Python engine are outside this release.
The ADaM and OMOP mapping helpers validate explicitly declared source rules;
they do not certify clinical coding, complete ascertainment or standards compliance.
Environment diagnostics are installation/operation checks, not signed IQ/OQ/PQ
qualification or FDA/NIH/NHS certification.
