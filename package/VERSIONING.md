# Version policy

The current review candidate is 0.11.1. Public fitting functions are cs_fit() and cs_control(). Use the coherentSurv namespace for all package calls. The native C++ namespace and installed headers follow the same identity.

The 0.11.0 release removes transitional fitting aliases. The numerical formulas and algorithms are unchanged from 0.10.2; parity checks and fresh package checks are recorded under reviews/release_0_11_0 in the repository. Earlier experiment values remain labeled by the version that produced them.

During the pre-1.0 phase, interface changes are described in NEWS.md. Scripts should record the package version, declared observation assumptions, time units, event eligibility and random seeds. Refit objects created by earlier development packages before using the current namespace.

CRAN acceptance is separate from local validation and GitHub distribution.
