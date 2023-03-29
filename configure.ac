#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])
AC_INIT([wam], [0.0.1], [ipe_devteam@noaa.gov])

AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE
AM_MAINTAINER_MODE

# Checks for programs.
AC_PROG_GREP
AC_PROG_SED

# Parse configure option
# Checks for required libraries
ext_libs=""
# - NETCDF4
AX_LIB_NETCDF4()
AS_IF([test "x$with_netcdf4" = xno],
  [AC_MSG_ERROR([Unable to find required NetCDF library])],
  []
)
AS_IF([test "x$with_netcdf4_fortran" = xno],
  [AC_MSG_ERROR([Fortran support is required for NetCDF library])],
  [ext_libs="${ext_libs} NETCDF"]
)

# - ESMF (get Fortran settings)
AC_LANG_PUSH([Fortran])
AX_LIB_ESMF
AS_IF([test "x$with_esmf" = xno],
  [AC_MSG_ERROR([Unable to find required ESMF library])],
  [ext_libs="${ext_libs} ESMF"]
)

# Looks for Fortran compiler
# Assume we are cross compiling by default to prevent issues on WCOSS-Cray
save_cross_compiling=${cross_compiling}
cross_compiling=yes
# Silence resulting warning
save_ac_tool_warned=${ac_tool_warned}
ac_tool_warned=yes
AC_PROG_FC([${FC} ${ESMF_FC} ${NETCDF4_FC} mpiifort mpifort mpif90 mpixlf90_r mpxlf90_r ${COMIO_FC}])
# - COMIO
AX_LIB_COMIO
AS_IF([test "x$with_comio" = x],
  [AC_MSG_ERROR([Unable to find required COMIO library])],
  [ext_libs="${ext_libs} COMIO"])

# Restore original cross compilation settings
cross_compiling=${save_cross_compiling}
ac_tool_warned=${save_ac_tool_warned}

# Checks if Fortran compiler supports MPI
AC_MSG_CHECKING([whether MPI is supported])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],[dnl
       use mpi
       integer, parameter :: comm = MPI_COMM_WORLD
  ])],
  [AC_DEFINE([HAVE_MPI], [1])
   AC_MSG_RESULT(yes)],
  [dnl No MPI found. Disable compiler.
   FC=
   AC_MSG_RESULT(no)
   AC_MSG_ERROR([Cannot find a suitable MPI Fortran compiler.])]
)

AC_FC_PP_SRCEXT([F90])
AX_COMPILER_VENDOR

# Check for required NEMSIO
AX_LIB_NEMSIO([yes])
AS_IF([test "x$with_nemsio" = xno],
  [AC_MSG_ERROR([Unable to find required NEMSIO library])],
  []
)

# Check for required BLAS library
AX_BLAS

# Add configure options
CHNK_RRTM=8
AC_MSG_CHECKING([RRTM maximum array size])
AC_ARG_WITH([chnk-rrtm],
  [AS_HELP_STRING([--with-chnk-rrtm=<size>],
    [RRTM maximum array size (default: 8)])],
  [AS_IF([test "`echo $with_chnk_rrtm | $SED -e 's/-g@<:@0-9@:>@*//g'`" = x],
     [CHNK_RRTM=$with_chnk_rrtm],[])],
  []
)
AC_MSG_RESULT($CHNK_RRTM)

# - safe floating point arithmetic
AC_ARG_ENABLE([safe-fp],
  [AS_HELP_STRING([--enable-safe-fp],[enable value-safe floating-point arithmetic (default: yes).])],
  [],[])

# - debug
AC_ARG_ENABLE([debug],
  [AS_HELP_STRING([--enable-debug],[enable debug mode (default: no).])],
  [enable_debug=${enableval}], [enable_debug=no])

# Checks for debug/optimization options
AS_IF([test "x$enable_debug" = xyes],[
  AC_MSG_CHECKING([compiler for debug options])
  save_FCFLAGS=${FCFLAGS}
  FCFLAGS=`echo $FCFLAGS | $SED -e 's/-O@<:@0-9s@:>@*//g'`
  FCFLAGS=`echo $FCFLAGS | $SED -e 's/-g@<:@0-9@:>@*//g'`
  debug_opts=
  AS_CASE(["$ax_cv_fc_compiler_vendor"],
    [gnu|clang], [debug_opts="-Wall -Wno-unused -Wno-unused-dummy-argument -fbacktrace"],
    [ibm],       [debug_opts="-qcheck -qfullpath"],
    [intel],     [debug_opts="-traceback -check pointer,stack,uninit"],
    [portland],  [debug_opts="-traceback"],
    []            dnl no debug options available for other platforms
  )
  AS_IF([test "x$ac_cv_prog_fc_g" = xyes], [debug_opts="-g ${debug_opts}"], [])
  FCFLAGS="-O0 ${debug_opts} ${FCFLAGS}"
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([])],
    [],
    [FCFLAGS=${save_FCFLAGS}
     debug_opts="unknown. Please add debug flags to FCFLAGS"]
  )
  AC_MSG_RESULT([${debug_opts}])
  ], [dnl default optimization settings, if none present
  AC_MSG_CHECKING([compiler for optimization options])
  fcflags_noopt=
  perf_opts=
  for fc_flag in $FCFLAGS
  do
    AS_CASE([$fc_flag],
      [-O@<:@0-9s@:>@*|-Ofast], [perf_opts=$fc_flag], dnl find optimization settings
      [fcflags_noopt="$fcflags_noopt $fc_flag"]
    )
  done
  AS_IF([test "x$perf_opts" = x], [dnl add default optimization settings
    AS_CASE(["$ax_cv_fc_compiler_vendor"],
      [gnu|clang], [dnl
                    perf_opts="-O2"
                    AS_IF([test "x$enable_safe_fp" != xno],
                           [perf_opts="${perf_opts} -fno-unsafe-math-optimizations"],[])],
      [ibm],       [perf_opts="-O2"],
      [intel],     [dnl
                    perf_opts="-O2"
                    AS_IF([test "x$enable_safe_fp" != xno],
                          [perf_opts="${perf_opts} -fp-model source -ftz -align array64byte -xCORE-AVX2 -qno-opt-dynamic-align"],[])],
      [portland],  [dnl
                    perf_opts="-O2"
                    AS_IF([test "x$enable_safe_fp" != xno],
                          [perf_opts="${perf_opts} -Kieee"],[])],
      []            dnl no debug options available for other platforms
    )], [])
  AS_IF([test "x$perf_opts" = x], [dnl optization flags unknown
    perf_opts=unknown],
    [dnl test optimization flags
    save_FCFLAGS=${fcflags_noopt}
    FCFLAGS="${fcflags_noopt} ${perf_opts}"
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([])],
      [],
      [FCFLAGS=${save_FCFLAGS}
       perf_opts=unknown]
    )
    ]
  )
  AC_MSG_RESULT([$perf_opts])
  ]
)

# Checks for compiler characteristics
save_FCFLAGS="${FCFLAGS}"
FCFLAGS=""
AX_FC_AUTODOUBLE
FC_AUTODOUBLE="${FCFLAGS}"
FCFLAGS="${save_FCFLAGS}"
AC_FC_MODULE_EXTENSION
AC_FC_MODULE_FLAG
AC_FC_PP_DEFINE
AC_OPENMP

# Checks whether platform is big endian. If not, finds compiler flag
# to convert data from big to little endian format (or vice-versa)
# during I/O from/to Fortran sequential unformatted files.
AC_LANG_POP
AC_C_BIGENDIAN
AC_LANG_PUSH([Fortran])
AS_IF([test "x$ac_cv_c_bigendian" = xno], [dnl look for compiler flag
  AX_FC_BIGENDIAN],[])

# Adds ESMF include paths
fc_inc=
for incpath in ${COMIO_FFLAGS} ${ESMF_FCFLAGS} ${NEMSIO_FCFLAGS} ${NETCDF4_FFLAGS}
do
  incpath=`echo "${incpath}" | $SED -e's/^-@<:@A-Z,a-z@:>@*\//\//g'`
  fc_inc="${fc_inc}${FC_MODINC}${incpath} "
done

AC_MSG_CHECKING([whether include flags work])
FCFLAGS="${fc_inc}${FCFLAGS}"
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([])],
  [AC_MSG_RESULT([yes])],
  [dnl ESMF-derived include flags may be corrupted. Abort.
   AC_MSG_RESULT([no])
   AC_MSG_FAILURE([unrecognized include flags])]
)

# Checks whether libraries are compatible with Fortran compiler
for lib in ${ext_libs}
do
  AC_MSG_CHECKING([whether ${lib} works])
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],[dnl
         use ${lib}
    ])],
    [AC_MSG_RESULT(yes)],
    [dnl ESMF is incompatible with selected Fortran compiler or not found. Abort.
     AC_MSG_RESULT(no)
     AC_MSG_FAILURE([${lib} library incompatible or not found])]
  )
done

AC_LANG_POP

# Set linker flags
LIBS="${LIBS} ${COMIO_FLIBS} ${NEMSIO_FLIBS} ${NETCDF4_FLIBS}"
LDFLAGS="${LDFLAGS} ${COMIO_LDFLAGS} ${NEMSIO_LDFLAGS} ${NETCDF4_LDFLAGS}"

WAM_DEP_LINK_OBJS="${BLAS_LIBS} ${libdir}/libwam.a"
WAM_DEP_SHRD_PATH="`echo ${LDFLAGS} | $SED -e's/-L//g'`"
WAM_DEP_SHRD_LIBS="`echo ${LIBS} | $SED -e's/-l//g'`"
AC_SUBST([WAM_DEP_LINK_OBJS])
AC_SUBST([WAM_DEP_SHRD_PATH])
AC_SUBST([WAM_DEP_SHRD_LIBS])

# Set F77 flags
AC_LANG_PUSH([Fortran 77])
save_FFLAGS="${FFLAGS}"
save_FCFLAGS="${FCFLAGS}"
FFLAGS="${FFLAGS} ${FCFLAGS} ${OPENMP_FCFLAGS}"
AC_PROG_F77([${FC}])
AC_FC_SRCEXT([f])
AC_FC_PP_SRCEXT([f])
FCFLAGS=""
AC_FC_FREEFORM
FCFLAGS_FREE="${FCFLAGS}"
FCFLAGS=""
AC_FC_FIXEDFORM
FCFLAGS_FIXED="${FCFLAGS}"
FCFLAGS=""
AX_FC_LINE_LENGTH([132])
FC_LINE_LENGTH="${FCFLAGS}"
FFLAGS="${save_FFLAGS} ${save_FCFLAGS} ${FCFLAGS_f}"
FCFLAGS="${save_FCFLAGS}"
AC_LANG_POP

# Additional configuration settings
AC_SUBST([FC_AUTODOUBLE])
AC_SUBST([FC_LINE_LENGTH])
AC_SUBST([FCFLAGS_FIXED])
AC_SUBST([FCFLAGS_FREE])

DEFS_PHYS="${FC_DEFINE}CHNK_RRTM=${CHNK_RRTM}"
AC_SUBST([DEFS_PHYS])

# Checks for auxiliary build tools
AC_PROG_INSTALL
AC_PROG_RANLIB
AM_PROG_AR

# Output Makefiles
AC_CONFIG_FILES([Makefile
                 src/Makefile
                 src/gen/Makefile
                 src/gsm/Makefile
                 src/gsm/dyn/Makefile
                 src/gsm/libutil/Makefile
                 src/gsm/phys/Makefile
                 src/io/Makefile
                 src/phys/Makefile
                 src/post/Makefile
                 src/share/Makefile])
AC_OUTPUT


AS_ECHO("")
AS_ECHO("${PACKAGE_NAME} version ${PACKAGE_VERSION} successfully configured.")
AS_ECHO("")
AS_ECHO("Compiler settings:")
AS_ECHO("  FC             = ${FC}")
AS_ECHO("  F77            = ${F77}")
AS_ECHO("  FCFLAGS        = ${FCFLAGS}")
AS_ECHO("  FFLAGS         = ${FFLAGS}")
AS_ECHO("  OPENMP_FCFLAGS = ${OPENMP_FCFLAGS}")
AS_ECHO("  LDFLAGS        = ${LDFLAGS}")
AS_ECHO("  LIBS           = ${LIBS}")
AS_ECHO("")
AS_ECHO("Libraries required:")
AS_ECHO("* BLAS")
AS_ECHO("  LIBS   : ${BLAS_LIBS}")
AS_ECHO("* COMIO")
AS_ECHO("  Version: ${COMIO_VERSION}")
AS_ECHO("  FC     : ${COMIO_FC}")
AS_ECHO("  FCFLAGS: ${COMIO_FFLAGS}")
AS_ECHO("  FLIBS  : ${COMIO_FLIBS}")
AS_ECHO("* NEMSIO")
AS_ECHO("  Version: ${NEMSIO_VERSION}")
AS_ECHO("  FCFLAGS: ${NEMSIO_FCFLAGS}")
AS_ECHO("  LDFLAGS: ${NEMSIO_LDFLAGS}")
AS_ECHO("  FLIBS  : ${NEMSIO_FLIBS}")
AS_ECHO("* NETCDF4")
AS_ECHO("  Version: ${NETCDF4_VERSION}")
AS_ECHO("  FC     : ${NETCDF4_FC}")
AS_ECHO("  FFLAGS : ${NETCDF4_FFLAGS}")
AS_ECHO("  LDFLAGS: ${NETCDF4_LDFLAGS}")
AS_ECHO("  FLIBS  : ${NETCDF4_FLIBS}")
AS_IF([test "$with_esmf" = yes], [dnl
  AS_ECHO("* ESMF")
  AS_ECHO("  Version: ${ESMF_VERSION}")
  AS_ECHO("  FC     : ${ESMF_FC}")
  AS_ECHO("  FCFLAGS: ${ESMF_FCFLAGS}")
  AS_ECHO("  Settings imported from ESMF Makefile fragment: ${ESMFMKFILE}")
  ],[])
AS_ECHO("")
AS_ECHO("Type 'make' to build.")
