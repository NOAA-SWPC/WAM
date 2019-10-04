# ===========================================================================
#      https://www.gnu.org/software/autoconf-archive/ax_lib_nemsio.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_LIB_NEMSIO()
#
# DESCRIPTION
#
#   This macro provides tests of the availability of the NEMSIO library.
#
#   The macro adds a --with-nemsio option accepting one of three values:
#
#     no   - do not check for the NEMSIO library.
#     yes  - do check for NEMSIO library in standard locations.
#     path - installation prefix for NEMSIO version.
#
#   If NEMSIO is successfully found, this macro calls
#
#     AC_SUBST(NEMSIO_VERSION)
#     AC_SUBST(NEMSIO_LDFLAGS)
#     AC_SUBST(NEMSIO_LIBS)
#     AC_SUBST(NEMSIO_FFLAGS)
#     AC_SUBST(NEMSIO_FLIBS)
#     AC_DEFINE(HAVE_NEMSIO)
#
#   It also sets
#
#     with_nemsio="yes"
#
#   If NEMSIO is disabled or not found, this macros sets
#
#     with_netcdf4="no"
#
#   Your configuration script can test $with_nemsio to take any further
#   actions. NEMSIO_F{FLAGS,LIBS} and NEMSIO_LDFLAGS should be used when
#   building Fortran applications.
#
#   To use the macro, one would add the following lines to "configure.ac"
#   before AC_OUTPUT:
#
#     dnl Check for NetCDF4 support
#     AX_LIB_NETCDF4()
#
#   One could test $with_nemsio for the outcome or display it as follows
#
#     echo "NEMSIO support:  $with_nemsio"
#
# LICENSE
#
#   Copyright (c) 2019 NOAA/ESRL/SWPC development team
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.

#serial 1

AC_DEFUN([AX_LIB_NEMSIO], [

dnl Check first argument is one of the recognized values.
dnl Fail eagerly if is incorrect as this simplifies case statements below.
if   test "m4_normalize(m4_default([$1],[]))" = ""        ; then
  ax_lib_nemsio_enable=""
elif test "m4_normalize(m4_default([$1],[]))" = "yes"  ; then
  ax_lib_nemsio_enable="yes"
elif test "m4_normalize(m4_default([$1],[]))" = "no"; then
  ax_lib_nemsio_enable="no"
else
  AC_MSG_ERROR([
    Unrecognized value for AX[]_LIB_NEMSIO within configure.ac.
    If supplied, argument 1 must be either 'yes' or 'no'.
  ])
fi

dnl AC_REQUIRE([AC_PROG_FC])

dnl Set default installation paths to blank
nemsio_includedir=""
nemsio_libdir=""
nemsio_library=""
nemsio_prefix=""

dnl Add a default --with-nemsio configuration option.
AC_ARG_WITH([nemsio],
  AS_HELP_STRING(
    [--with-nemsio=[yes/no/PATH]],
    m4_case(m4_normalize([$1]),
            [base directory of NEMSIO installation])
  ),
  [if test "$withval" = "no"; then
     with_nemsio="no"
   elif test "$withval" = "yes"; then
     with_nemsio="yes"
   else
     with_nemsio="yes"
     ax_lib_nemsio_prefix="${withval}"
   fi],
  [with_nemsio=$ax_lib_nemsio_enable]
)

dnl Set defaults to blank
NEMSIO_FFLAGS=""
NEMSIO_FLIBS=""
NEMSIO_LDFLAGS=""
NEMSIO_LIBS=""
NEMSIO_VERSION=""

dnl Try to find NEMSIO installation paths
ax_lib_nemsio_includdedir=""
ax_lib_nemsio_libdir=""
ax_lib_nemsio_library=""

if test "$with_nemsio" = "yes"; then
  AC_MSG_CHECKING([for NEMSIO library])
  if test -z "$ax_lib_nemsio_prefix"; then
    dnl Retrieve settings from environment
    if test -z "$NEMSIO_INC"; then
      AC_MSG_ERROR([Environment variable NEMSIO_INC undefined])
    else
      ax_lib_nemsio_includedir=$NEMSIO_INC
    fi
    if test -z "$NEMSIO_LIB"; then
      AC_MSG_ERROR([Environment variable NEMSIO_LIB undefined])
    else
      ax_lib_nemsio_libdir=`AS_DIRNAME([$NEMSIO_LIB])`
      ax_lib_nemsio_library=$NEMSIO_LIB
    fi
    if ! test -z "$NEMSIO_VER"; then
      ax_lib_nemsio_version=$NEMSIO_VER
    fi
  else
    dnl Assume standard installation if prefix is provided
    ax_lib_nemsio_includedir=${ax_lib_nemsio_prefix}/include
    ax_lib_nemsio_libdir=${ax_lib_nemsio_prefix}/lib
    ax_lib_nemsio_library=${ax_lib_nemsio_libdir}/libnemsio.a
    ax_lib_nemsio_version=unknown
  fi

  dnl Verify installation paths
  for dir in ${ax_lib_nemsio_includedir} ${ax_lib_nemsio_libdir} ; do
    if ! test -d $dir; then
      AC_MSG_ERROR(["Unable to find NEMSIO installation path ${dir}"])
    fi
  done

  dnl Verify library file
  if ! test -f $ax_lib_nemsio_library ; then
    AC_MSG_ERROR(["Unable to find NEMSIO library ${ax_lib_nemsio_library}"])
  fi

  AC_MSG_RESULT([yes (version $[ax_lib_nemsio_version])])
  
  dnl Check if the library works
  AC_MSG_CHECKING([if NEMSIO library works])
  AC_LANG_PUSH([Fortran])
  ax_lib_nemsio_save_FCFLAGS="${FCFLAGS}"
  FCFLAGS="${FCFLAGS} -I $ax_lib_nemsio_includedir"
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],[dnl
       use nemsio_module
    ])],
    [AC_MSG_RESULT(yes)],
    [dnl library cannot compile
     AC_MSG_RESULT(no)
     FCFLAGS=$ax_lib_nemsio_save_FCFLAGS
     AC_MSG_ERROR([Incompatible NEMSIO library.])]
  )
  FCFLAGS=$ax_lib_nemsio_save_FCFLAGS
  AC_LANG_POP()

  NEMSIO_FCFLAGS=-I$ax_lib_nemsio_includedir
  NEMSIO_LDFLAGS=-L$ax_lib_nemsio_libdir
  NEMSIO_LIBS=-l$( echo `AS_BASENAME([$ax_lib_nemsio_library])` | $SED -e 's/^lib//' -e 's/\.a$//')
  NEMSIO_FLIBS=$NEMSIO_LIBS
  NEMSIO_VERSION=$ax_lib_nemsio_version

  AC_SUBST([NEMSIO_FCFLAGS])
  AC_SUBST([NEMSIO_FLIBS])
  AC_SUBST([NEMSIO_LDFLAGS])
  AC_SUBST([NEMSIO_LIBS])
  AC_SUBST([NEMSIO_VERSION])
  AC_DEFINE([HAVE_NEMSIO], [1], [Defined if you have NEMSIO support])
fi
])
