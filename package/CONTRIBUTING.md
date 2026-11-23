# Contributing

Use a synthetic minimal example for issues. Each change to a numerical kernel
must identify its mathematical quantity, assumptions, oracle fixture, tolerance
and failure behavior. Compare estimands before claiming a speed improvement.

Build from the source tarball and run R CMD check. Core C++ must also compile
under -std=c++17 -Wall -Wextra -pedantic -Wconversion. Add meaningful regression
tests for corrected defects and preserve failed scientific checks in reports.

Do not add a new dependency without explaining why it is needed. Optional
visualization and comparison packages belong in Suggests. Do not add clinical
or identifying data to the repository. Two internal review lenses inspect
changes to mathematical scope or numerical kernels; this is not independent
peer review.
