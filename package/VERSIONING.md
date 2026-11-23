# Versioning

The current unpublished release candidate is 0.10.1, following the separately
preserved 0.10.0 and 0.9.1 candidates and internal 0.9.0.9000 snapshots. Version 1.0.0 is reserved for a
release candidate whose numerical, documentation, packaging and platform gates
have been evaluated. A version string is not evidence of CRAN acceptance.

For public releases from 1.0 onward, patch releases repair defects without
intentional API changes and minor releases add backwards-compatible functions.
The explicitly numbered pre-1.0 candidate 0.10.1 contains targeted additive
R interfaces while preserving the native kernels. Major releases may change estimands or
interfaces. Deprecated arguments warn for at least one minor release before
removal. Every release reruns deterministic oracle fixtures and public examples.

Version 1.1 plans include Python bindings, external schema adapters, valuation
helpers and parallel execution. Version 2 plans concern causal and spatial
extensions with separate identifying assumptions. The intervention argument
must produce a structured unsupported-feature error in version 1.

Before releases, inspect reverse dependencies and compare benchmark regressions.
The planned 90-day adoption review measures actual downloads, issue responses,
external checks and citations; none is a present result.

