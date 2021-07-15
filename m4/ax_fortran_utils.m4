# --------------------------------------------------------------------------
#
# SYNOPSIS
#
#   AX_FC_AUTODOUBLE([ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
#
# Looks for a compiler flag to turn all real variables automatically
# into double precision variables and adds the flag to FCFLAGS.
# On success the variable ax_cv_fc_autodouble holds the compiler flag
# otherwise the string "none". Call ACTION-IF-SUCCESS
# (defaults to nothing) if successful (i.e. can compile code using
# the specific compiler flag) and ACTION-IF-FAILURE (defaults to 
# failing with an error message) if not.
#
# The known flags are:
#                -r8: Intel compiler (ifort, ifc) and g95 compiler
#   -fdefault-real-8: GNU Fortran compiler (gfortran)
#      -Wf"-A idbl4": NEC SX-9 compiler with MPI (sxmpif90)
#                -ew: NEC SX-9 compiler with iso_c_binding (sxf90)
# --------------------------------------------------------------------------
AC_DEFUN([AX_FC_AUTODOUBLE],[
    AC_LANG_PUSH(Fortran)dnl
    AC_CACHE_CHECK([for Fortran flag to autodouble real numbers],
        ax_cv_fc_autodouble,[
	ax_cv_fc_autodouble=unknown
	ax_fc_autodouble_FCFLAGS_save="$FCFLAGS"
	for ac_flag in -r8 -fdefault-real-8 '-Wf"-A idbl4"' -ew
	do
	    test "x$ac_flag" != xnone && FCFLAGS="$ax_fc_autodouble_FCFLAGS_save $ac_flag"
	    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],[
        real(kind=8) :: y
        y=1.0D0
        call testdbl(y)
        contains
            subroutine testdbl(x)
            real :: x
            print *,x
            end subroutine testdbl
])],[             ax_cv_fc_autodouble=$ac_flag; break
            ])
        done
        rm -f conftest.err conftest.$ac_objext conftest.$ac_ext
        FCFLAGS=$ax_fc_autodouble_FCFLAGS_save
    ])
    AS_VAR_IF([ax_cv_fc_autodouble],[unknown],[dnl
        ax_cv_fc_autodouble=""
        m4_default([$2],[
            AC_MSG_ERROR([Fortran compiler does not support autodouble])])
    ],[dnl
        AS_VAR_IF([ax_cv_fc_autodouble],[none],[dnl
           ax_cv_fc_autodouble=""])
        m4_default([$1],[dnl
	    FCFLAGS="$FCFLAGS $ax_cv_fc_autodouble"])
    ])
    AC_LANG_POP(Fortran)dnl
])# AX_FC_AUTODOUBLE


# AX_FC_LINE_LENGTH([LENGTH], [ACTION-IF-SUCCESS],
#                   [ACTION-IF-FAILURE = FAILURE])
# ------------------------------------------------
# Look for a compiler flag to make the Fortran (FC) compiler accept long lines
# in the current (free- or fixed-format) source code, and adds it to FCFLAGS.
# The optional LENGTH may be 80, 132 (default), or `unlimited' for longer
# lines.  Note that line lengths above 254 columns are not portable, and some
# compilers (hello ifort) do not accept more than 132 columns at least for
# fixed format.  Call ACTION-IF-SUCCESS (defaults to nothing) if successful
# (i.e. can compile code using new extension) and ACTION-IF-FAILURE (defaults
# to failing with an error message) if not.  (Defined via DEFUN_ONCE to
# prevent flag from being added to FCFLAGS multiple times.)
# You should call AC_FC_FREEFORM or AC_FC_FIXEDFORM to set the desired format
# prior to using this macro.
#
# The known flags are:
# -f{free,fixed}-line-length-N with N 72, 80, 132, or 0 or none for none.
# -ffree-line-length-none: GNU gfortran
# -ffree-line-length-huge: g95 (also -ffixed-line-length-N as above)
#       -qfixed=132 80 72: IBM compiler (xlf)
#                -Mextend: Cray
#            -132 -80 -72: Intel compiler (ifort)
#                          Needs to come before -extend_source because ifort
#                          accepts that as well with an optional parameter and
#                          doesn't fail but only warns about unknown arguments.
#          -extend_source: SGI compiler
#  -W, -WNN (132, 80, 72): Absoft Fortran
#     +es, +extend_source: HP Fortran (254 in either form, default is 72 fixed,
#                          132 free)
#            -w, (-)-wide: Lahey/Fujitsu Fortran (255 cols in fixed form)
#                      -e: Sun Fortran compiler (132 characters)
#                    -132: NAGWare
#         -72, -f, -Wf,-f: f2c (a weak form of "free-form" and long lines).
#                  /XLine: Open Watcom
AC_DEFUN_ONCE([AX_FC_LINE_LENGTH],
[AC_LANG_PUSH([Fortran])dnl
m4_case(m4_default([$1], [132]),
  [unlimited], [ax_fc_line_len_string=unlimited
                       ax_fc_line_len=0
                       ax_fc_line_length_test='
      subroutine longer_than_132(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,'\
'arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,arg17,arg18,arg19)'],
  [132],            [ax_fc_line_len=132
                       ax_fc_line_length_test='
      subroutine longer_than_80(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,'\
'arg10)'],
  [80],             [ax_fc_line_len=80
                       ax_fc_line_length_test='
      subroutine longer_than_72(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9)'],
  [m4_warning([Invalid length argument `$1'])])
: ${ax_fc_line_len_string=$ax_fc_line_len}
AC_CACHE_CHECK(
[for Fortran flag needed to accept $ax_fc_line_len_string column source lines],
               [ax_cv_fc_line_length],
[ax_cv_fc_line_length=unknown
ax_fc_line_length_FCFLAGS_save=$FCFLAGS
for ax_flag in none \
               -ffree-line-length-none -ffixed-line-length-none \
               -ffree-line-length-huge \
               -ffree-line-length-$ax_fc_line_len \
               -ffixed-line-length-$ax_fc_line_len \
               -qfixed=$ax_fc_line_len -Mextend \
               -$ax_fc_line_len -extend_source \
               -W$ax_fc_line_len -W +extend_source +es -wide --wide -w -e \
               -f -Wf,-f -xline
do
  test "x$ax_flag" != xnone && FCFLAGS="$ax_fc_line_length_FCFLAGS_save $ax_flag"
  AC_COMPILE_IFELSE([[$ax_fc_line_length_test
      end subroutine]],
                    [ax_cv_fc_line_length=$ax_flag; break])
done
rm -f conftest.err conftest.$ax_objext conftest.$ax_ext
FCFLAGS=$ax_fc_line_length_FCFLAGS_save
])
if test "x$ax_cv_fc_line_length" = xunknown; then
  m4_default([$3],
             [AC_MSG_ERROR([Fortran does not accept long source lines], 77)])
else
  if test "x$ax_cv_fc_line_length" != xnone; then
    FCFLAGS="$FCFLAGS $ax_cv_fc_line_length"
  fi
  $2
fi
AC_LANG_POP([Fortran])dnl
])# AX_FC_LINE_LENGTH


# --------------------------------------------------------------------------
#
# SYNOPSIS
#
#   AX_FC_BIGENDIAN([ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
#
# Looks for a compiler flag to set the conversion format for numerical
# data in unformatted I/O.
# into double precision variables and adds the flag to FCFLAGS.
# On success the variable ax_cv_fc_bigendian holds the compiler flag
# otherwise the string "none". Call ACTION-IF-SUCCESS
# (defaults to nothing) if successful (i.e. can compile code using
# the specific compiler flag) and ACTION-IF-FAILURE (defaults to
# failing with an error message) if not.
#
# The known flags are:
#    -convert big_endian: Intel Fortran (ifort), Compaq Fortran (f90) compilers
#   -fconvert=big-endian: GNU Fortran compiler (gfortran)
#              -qufmt=be: IBM Fortran compiler (xlf)
#            -byteswapio: PGI Fortran compiler (pgfortran)
# --------------------------------------------------------------------------
AC_DEFUN([AX_FC_BIGENDIAN],[
    AC_LANG_PUSH(Fortran)dnl
    AC_CACHE_CHECK([for Fortran flag for big-endian I/O],
        ax_cv_fc_bigendian,[
	ax_cv_fc_bigendian=unknown
	ax_fc_bigendian_FCFLAGS_save="$FCFLAGS"
	for ac_flag in "-convert big_endian" -fconvert=big-endian -qufmt=be -byteswapio
	do
        FCFLAGS="$ac_flag $ax_fc_bigendian_FCFLAGS_save"
	    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],[])],[dnl
		ax_cv_fc_bigendian=$ac_flag; break
            ])
        done
        rm -f conftest.err conftest.$ac_objext conftest.$ac_ext
        FCFLAGS=$ax_fc_bigendian_FCFLAGS_save
    ])
    AS_VAR_IF([ax_cv_fc_bigendian],[unknown],[dnl
        ax_cv_fc_bigendian=""
        m4_default([$2],[
            AC_MSG_ERROR([Fortran compiler does not support big-endian I/O])])
    ],[dnl
        AS_VAR_IF([ax_cv_fc_bigendian],[none],[dnl
           ax_cv_fc_bigendian=""])
        m4_default([$1],[dnl
	    FCFLAGS="$FCFLAGS $ax_cv_fc_bigendian"])
    ])
    AC_LANG_POP(Fortran)dnl
])# AX_FC_BIGENDIAN
