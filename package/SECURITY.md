# Security reporting

Send reports privately to the maintainer email in DESCRIPTION. Include only
synthetic reproductions. Never include credentials or patient-level data.

The target for initial triage is 72 hours. No historical response-time claim
is made. Security fixes should receive a patch version, a regression test and
a documented advisory once disclosure is appropriate.

Examples and tests must not access a network or start external viewers.
Package functions must not terminate the R process or silently write outside
the session temporary directory. User-requested exports must require an
explicit path.
