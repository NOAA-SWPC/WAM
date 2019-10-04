#!/bin/ksh

# copy script for column
#JR These files were built and stored in library libw3_4.a but can't be
#JR used here because they need to be built with 8-byte reals

#TBH Hack since IBM does not include rsync in default path
#TBH function update { rsync -ut $1 .; }
function update { cp -f $1 .; }

# linking the .mod files works for ifort but not for the Lahey compiler.
update ../../w3/iw3jdn.f
update ../../w3/w3fs26.f
update ../../w3/w3movdat.f
update ../../w3/w3reddat.f
