# Global variable registration intentionally removed during active development.
#
# This file previously contained utils::globalVariables(...) declarations
# used to suppress R CMD check notes caused by NSE / tidy-eval expressions.
#
# During v0.2.x-v0.4.x development, internal variable names and pipelines
# are still changing, so stale registrations are not maintained.
#
# Before a release-quality R CMD check:
#   1. run R CMD check on the built source tarball,
#   2. inspect any "no visible binding for global variable" notes,
#   3. prefer explicit .data$column / tidy-eval references where practical,
#   4. register only the remaining unavoidable NSE variables.
#
# Do not add runtime state or scientific defaults here.
