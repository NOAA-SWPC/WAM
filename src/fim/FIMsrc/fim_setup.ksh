#!/bin/ksh

CONTEXT="fim_setup.ksh"

# Purpose:
#   Sets up environment for FIM build or run.  
#
# Usage:
#   This script should be executed from another script via the ksh "." (source)
#   command. It takes zero, one, or two arguments:  
#     . fim_setup.ksh
#     . fim_setup.ksh [fim_configuration]
#     . fim_setup.ksh [fim_configuration] [loquacious]
#
#   The "fim_configuration" argument is used to select non-default 
#   environmental settings.  
#
#   When present, the "loquacious" argument causes intermediate results to be 
#   printed.  

loquacious="false"
use_modules="true"
module_ksh="/opt/modules/Modules/default/init/ksh"
mpif90="mpif90"

# Process command-line arguments if present.
usage_msg="Usage: . fim_setup.ksh [fim_configuration] [loquacious]"

case "$#" in
  2) fim_configuration="$1"; loquacious="true";;
  1) fim_configuration="$1";;
  0) fim_configuration="openmpi";;
  *) print "$usage_msg"; exit 1;; 
esac

#TODO:  Paul, we need to call the makefim "check" function here to 
#TODO:  validate fim_configuration

# Set environment variables and reset $use_modules if needed.
case "$fim_configuration" in
  "vapor") # vapor -q64
    # Add FIM_ESMF_INSTALL_LIBDIR_ABSPATH later...
    use_modules="false"
    ;;
  "devccs") # cirrus/stratus -q64
    # Add FIM_ESMF_INSTALL_LIBDIR_ABSPATH later...
    use_modules="false"
    ;;
  "bluefire") # bluefire -q64
    # Add FIM_ESMF_INSTALL_LIBDIR_ABSPATH later...
    use_modules="false"
    ;;
  "nems") # ifort+mvapich: jet default + ESMF
    # Set location of esmf.mk.
    export FIM_ESMF_INSTALL_LIBDIR_ABSPATH="/home/rosinski/esmf-3.1.0rp2/lib/libO/Linux.intel.64.mpich2.default"
    ;;
  "ranger") # ranger, mvapich/1.01
    ;;
  "linuxpcgnu")
    use_modules="false"
    export NETCDF="$HOME/x86_64"
    #JR NOTE: FIM under gfortran requires v 4.4 or greater of the compiler.
    #JR Earlier revs. didn't allow "allocatables" inside derived types.
    #JR By setting MPICH_F90 below, you can tell mpirun to use a different 
    #JR compiler if it's needed to address such issues.
    #JR export MPICH_F90="gfortran44"
    ;;
  "macgnu")
    use_modules="false"
    export NETCDF="$HOME"
    #JR NOTE: FIM under gfortran requires v 4.4 or greater of the compiler.
    #JR Earlier revs. didn't allow "allocatables" inside derived types.
    #JR By setting MPICH_F90 below, you can tell mpirun to use a different 
    #JR compiler if it's needed to address such issues.
    #JR export MPICH_F90="gfortran44"
    ;;
  "jaguarintel")
    use_modules="true"
    module_ksh="/opt/modules/3.1.6/init/ksh"
    ;;
  "jaguargnu")
    use_modules="true"
    module_ksh="/opt/modules/3.1.6/init/ksh"
    ;;
  "frostintel")
    mpif90="ifort"
    use_modules="false"
    ;;
  *)
    ;;
esac

# Module setup
test "$use_modules" == "true" && . $module_ksh

# load default modules 
if [[ "$fim_configuration" != "ranger" && "$fim_configuration" != "jaguarintel" && "$fim_configuration" != "jaguargnu" ]]
then
  if [[ "$use_modules" == "true" ]]
  then
    module purge # unload all modules 
    module load wjet
    # Now we know what modules we are switching from...  
  fi
fi

# Switch modules if needed.  
if [[ "$use_modules" == "true" ]]
then
  case "$fim_configuration" in
    "debug")
      # ifort-11.1+mvapich2-1.4.1 (jet default)
      ;;
    "jaguarintel")
      mpif90="ftn"
      # The following will break if/when PrgEnv-pgi is no longer the default
      module switch PrgEnv-pgi PrgEnv-intel
      # Current default on jaguar (11.1.046) fails on phy_init.F90 with "internal compiler error"
      # so use the most recent intel compiler version available
      module switch intel intel/11.1.064
      module unload netcdf
      module load netcdf
      ;;
    "jaguargnu")
      mpif90="ftn"
      # The following will break if/when PrgEnv-pgi is no longer the default
      module switch PrgEnv-pgi PrgEnv-gnu
      module unload netcdf
      module load netcdf
      ;;
    "lahey")
      # lahey+mvapich2-1.4.1
      module switch intel lahey/8.10b
      ;;
    "mvapich")
      # ifort-11.1+mvapich2-1.4.1 (jet default)
      ;;
    "nems")
      # ifort-11.1+mvapich2-1.4.1 (jet default)
      ;;
    "openmpi")
      # ifort-11.1+openmpi-1.4.1
      module switch mvapich2 openmpi/1.4.1-intel-11.1
      ;;
    "ranger")
      # ranger, mvapich/1.01
      module unload pgi mvapich
      module load intel
      module load mvapich 
      module load netcdf
      ;;
  *)
      ;;
  esac
fi

# List modules iff loquacious switch is set and modules are used
if  [[ "$loquacious" == "true" && "$use_modules" == "true" ]]
then
  module list
  print "whence $mpif90: $(whence $mpif90)"
  #TODO: IBM-specific software stack listing here...
fi

return 0
