#!/bin/ksh
#
# Extracts FIM build configuration from $1, which must end with "FIMsrc_*",
# where "*" is the actual FIM build configuration (i.e. "openmpi", "lahey",
# etc.).

function fail { print $@; exit 1; }

test "$#" -ge 1 || fail "SRCDIR is empty: Please fix in namelist."
buildconfig=$(basename $1)
print $buildconfig | grep -q FIMsrc_..* || fail "Could not extract buildconfig \
from SRCDIR in namelist: SRCDIR must end with FIMsrc_*."
print $buildconfig | cut -d_ -f2-

return 0
