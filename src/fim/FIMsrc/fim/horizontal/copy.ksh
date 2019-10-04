#!/bin/ksh

# copy script for horizontal

#JR "queue" is initially empty, so "queuefile" just appends its arg to "queue"
function queuefile { queue="$queue $1 "; } # note trailing space in string
#TODO:  DRY with FIMrun/functions.ksh
function fail { test -n "$1" && print "ERROR: $@"; exit 1 ; }
#TBH Hack since IBM does not include rsync in default path
#TBH function update { rsync -ut $1 . || fail ; }
function update { cp -f $1 . || fail ; }


toplevel_objs_file="FIM_HORIZONTAL_OBJS_TOP"
tmpfile=$0.$$.tmp

# Everything that goes through "queuefile" and only those files will be copied to Horizontal/ and 
# actually get compiled. 

queue="" # reset queue
queuefile ../../post/pop/postdata.F90
queuefile ../../prep/ss2icos/mktopo.F90
queuefile ../../utils/module_initial_chem_namelists.F90
queuefile ../../utils/wtinfo.F90
queuefile ../column_chem/module_chemvars.F90
queuefile ../column_chem/module_initial_chem_namelist_defaults.F90
queuefile ../horizontal/$toplevel_objs_file
queuefile ../horizontal/FIM_HORIZONTAL_OBJS
queuefile ../horizontal/Makefile
queuefile ../horizontal/SMS_Module_Lookup.txt
for file in $(ls -1 ../../cntl/*.F90)          ; do queuefile $file; done
for file in $(ls -1 ../../prep/ss2icos/*.F90)  ; do queuefile $file; done
for file in $(ls -1 ../horizontal/*.F90)       ; do queuefile $file; done
for file in $(ls -1 ../horizontal/*.c)         ; do queuefile $file; done

# Update the queued files.

for file in $queue; do update $file; done

# Copy file(s) that need new names

cp -pf ../wrfphys/module_wrfphysvars.F module_wrfphysvars.F90 || fail

# Handle NEMS build if necessary.

if [[ -n "$NEMS" ]]; then

  print "$PWD/$0: NEMS BUILD"

  fimdriver=$toplevel_objs_file.fim
  ncep_root="../framework/nems"

  # Remove standard FIM driver .F90 source files.

  test -f $fimdriver || cp -pf $toplevel_objs_file $fimdriver || fail
  for file in $(grep OBJS_TOP $fimdriver | cut -d= -f2)
  do
    rm -f ${file%.*}.F90 || fail
  done

  # Update these now to avoid adding them to OBJS_TOP in the loop below.

  update $ncep_root/kind.inc
  for file in $(ls -1 $ncep_root/*.h); do update $file; done

  # Grab FIM component code.

  # update *.f90 files without name change
  for file in $(ls -1 $ncep_root/*.f90); do update $file; done

  # rename *.F90 to *.f90 to avoid running PPP on NEMS files
  # (Makefile rule for *.f90 avoids PPP.)  
  for file in $(ls -1 $ncep_root/*.F90)
  do
    dstfile=$(basename $file)
    cp -pf $file ${dstfile%.*}.f90 || fail
  done

  # Replace toplevel objects file.
  objs="OBJS_TOP = " # append .o names to this initial string
  for file in $(ls -1 $ncep_root/*.F90 $ncep_root/*.f90)
  do
    file=$(basename $file)
    objs="$objs ${file%.*}.o" # replace .F90 and .f90 extensions with .o (ksh feature)
  done
  print $objs > $toplevel_objs_file  

fi # NEMS handler

# Deduce dependencies.

rm -f Filepath Srcfiles || fail
echo "." > Filepath
ls *.[fF]90 > Srcfiles
$MKDEPENDS -m -d module_decomp.o Filepath Srcfiles > FIM_HORIZONTAL_DEPENDENCIES || fail
