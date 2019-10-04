# This code provides shared ksh functions for the run-automation scripts.

function check_nems
{
  # Check that run configuration is NEMS-compatible
  test -z "$FC" && fail "$0: FC undefined."
  test -z "$COMPARE_VAR_ON" && fail "$0: COMPARE_VAR_ON undefined."
  if [[ "$FC" == "nems" ]]
  then
    test "$COMPARE_VAR_ON" == ".true." && \
      fail "Cannot use NEMS with COMPARE_VAR."
  # NEMS doesnt yet work when restart enabled
    get_nl_value $fimnamelist OUTPUTnamelist readrestart READRESTART
    test "$READRESTART" != ".false." && \
      fail "Cannot use NEMS with READRESTART."
  fi
}

function compare_var_setup
{
  # Please see ../README for a description of how to set up FIM to work with the
  # SMS COMPARE_VAR feature.

  # TODO Replace this with encapsulated calls to SMS utility programs!  

  COMPARE_VAR_ON="false"
  COMPARE_VAR_NTASKS_1="0"
  COMPARE_VAR_NTASKS_2="0"
  smsnamelist="SMSnamelist"

  if [[ -f "$smsnamelist" ]]
  then
    get_nl_value_unquoted $smsnamelist SMSnamelist compare_var_on COMPARE_VAR_ON
    test -z "$COMPARE_VAR_ON" && \
      fail "Cannot determine compare_var_on from $smsnamelist."
    if [[ "$COMPARE_VAR_ON" == ".true." ]]
    then
      # sed: remove everything before ':'
      print $(cd $fimnamelist_dir && $SRCDIR/bin/get_num_cores | \
        grep "root_own_node:") | sed 's/^.*://' | read TF || \
        fail "Could not get root_own_node."
      test "$TF" == "TRUE" && \
        fail "COMPARE_VAR requires root_own_node to be false."
      get_nl_value_unquoted $smsnamelist SMSnamelist compare_var_ntasks_1 \
        COMPARE_VAR_NTASKS_1
      test -z "$COMPARE_VAR_NTASKS_1" && \
        fail "Cannot determine compare_var_ntasks_1 from $smsnamelist."
      get_nl_value_unquoted $smsnamelist SMSnamelist compare_var_ntasks_2 \
        COMPARE_VAR_NTASKS_2
      test -z "$COMPARE_VAR_NTASKS_2" && \
        fail "Cannot determine compare_var_ntasks_2 from $smsnamelist."
    fi
  fi 

  # Write task info: needed only for compare_var test and diagnostic print
  # sed: remove everything before ':'
  print $(cd $fimnamelist_dir && $SRCDIR/bin/GetWriteTaskInfo | \
    grep "num_write_tasks") | sed 's/^.*://' | read nwt || \
    fail "Could not get num_write_tasks."

  # Perform error-checks on COMPARE_VAR options if COMPARE_VAR is enabled
  if [[ "$COMPARE_VAR_ON" == ".true." ]]
  then
      # Cannot use COMPARE_VAR with write tasks
      test "$nwt" -gt 0 && fail "Cannot use write tasks with COMPARE_VAR."
      # Verify that FIMnamelist and SMSnamelist are set so that
      # ComputeTasks == (COMPARE_VAR_NTASKS_1 + COMPARE_VAR_NTASKS_2)
      let "testcomputetasks=$COMPARE_VAR_NTASKS_1+$COMPARE_VAR_NTASKS_2" || \
        fail "Arithmetic error."
      test "$testcomputetasks" -ne "$PES" && \
        fail "COMPARE_VAR requires that ComputeTasks (in FIMnamelist) == \
COMPARE_VAR_NTASKS_1 + COMPARE_VAR_NTASKS_2 (in SMSnamelist). Please correct."
      # Cannot use COMPARE_VAR in a serial run
      test "$parallelism" == "serial" && \
        fail "Cannot execute serial code with COMPARE_VAR yet."
  fi

  # Create a more informative job name if COMPARE_VAR is used
  test "$COMPARE_VAR_ON" == ".true." && \
    ParaSuffix="cv.$COMPARE_VAR_NTASKS_1.vs.$COMPARE_VAR_NTASKS_2"

}

function context_peek
{
  # Print the CONTEXT name most recently pushed onto the stack.
  print $CONSTCK | sed "s/^\([^:][^:]*\).*/\1/"
}

function context_pop
{
  # Pop a value of the context stack and set CONTEXT to its value.
  CONTEXT=$(context_peek)
  CONSTCK=$(print $CONSTCK | sed "s/^$CONTEXT[:]*//")
}

function context_push
{
  # Push a context value onto the stack.
  test -z "$1" && fail "No argument to push supplied."
  test -z "$CONSTCK" && CONSTCK=$1 || CONSTCK="$1:$CONSTCK"
}

function copyfiles
{
  test ! -d "$2" && fail "$2 is not a directory."
  for fil in $(ls -1 $1)
  do
    pfil="${1}/${fil}"
    if [[ -f "$pfil" ]]
    then
      cp $pfil $2 || fail "Cannot cp $pfil -> $2."
    fi
  done
}

function endian_big
{
  # Enable big-endian handling of the space-separated list of logical unit
  # numbers given as the function's argument.
  # sed: replace spaces with commas
  typeset luns1=$(print "$@" | sed 's/ /,/g')
  # sed: insert '-T' at the beginning and ',-T' anywhere a space is found
  typeset luns2=$(print "$@" | sed 's/^/-T/;s/ /,-T/g')
  export F_UFMTENDIAN="big:$luns1" # intel
  export FORT90L="-Wl,$luns2" # lahey
  export GFORTRAN_CONVERT_UNIT="big_endian:$luns1" # gfortran
  print "F_UFMTENDIAN=$F_UFMTENDIAN FORT90L=$FORT90L \
GFORTRAN_CONVERT_UNIT=$GFORTRAN_CONVERT_UNIT"
}

function endian_little
{
  # Enable little-endian handling of the space-separated list of logical unit
  # numbers given as the function's argument.
  # sed: replace spaces with commas
  typeset luns1=$(print "$@" | sed 's/ /,/g')
  export XLFRTEOPTS="ufmt_littleendian=$luns1" # ibm
  print "XLFRTEOPTS=$XLFRTEOPTS"
}

function endian_reset
{
  # Disable all endianness control variables.
  unset F_UFMTENDIAN          # intel
  unset FORT90L               # lahey
  unset GFORTRAN_CONVERT_UNIT # gfortran
  unset XLFRTEOPTS            # ibm
}

function errhandler
{
  # Handle trapped errors. The line number where the error occurred is expected
  # as the sole argument. Print an error message including the failed line
  # number, then re-enable xtrace, which presumably was disabled by trap_on().
  # This function is meant to be called by the trap mechanism, not directly by
  # the sourcing script. It is also expected that the sourcing script will exit
  # with an informative error message after this function returns. For example:
  #
  # cd /no/such/directory || fail "Cannot cd to /no/such/directory"
  #
  # will trigger errhandler() if trap_on() has previously been called. The
  # failed line number will be reported and command passed back to the sourcing
  # script, which will then fail with a more-informative message.
  typeset LINE=$1
  print "$CONTEXT[$LINE]: An error occurred, see stdout" >&2
  set -o xtrace
}

function fail
{
  # Print a failure message and terminate.
  test -n "$1" && print "ERROR: $@"
  exit 1
}

function get_nl_value
{
  # Return, in the specified variable, the value corresponding to the given
  # namelist file, namelist and key. Any leading or trailing double or single
  # quotes are removed.
  #
  # usage: get_nl_value namelist_file namelist key variable_to_set
  #
  # NOTE: This is a naive method for dealing with namelists. For example, a '!'
  #       inside a string value will be seen as the start of a comment, which
  #       will be deleted from that line. It's assumed that namelist names in
  #       the form '&namelist' are on lines by themselves, though the standard
  #       does not require this. And the wholesale conversion of tabs to spaces
  #       will modify strings containing literal tabs. This is ok now for FIM,
  #       but may someday need to be generalized.

  # Check the provided arguments.
  test -f "$1" || fail "$0: namelist file '$1' not found."
  test -z "$2" && fail "$0: no namelist supplied."
  test -z "$3" && fail "$0: no key supplied."
  test -z "$4" && fail "$0: no output variable name supplied."
  # AIX sed does not recognize '\t'. Convert tabs to spaces, then consider only
  # spaces as whitespace.
  value=""
  cat $1 | tr '\t' ' ' | while read line
  do
    # sed delete: 1. comments, 2. leading whitespace, 3. trailing whitespace,
    #             4. namelist-termination lines, 5. blank lines.
    x=$(print "$line" | sed "s/\!.*$//g;s/^ *//g;s/ *$//g;/^\//d;/^$/d")
    # Try to recognize this line as a namelist name.
    newnl=$(print "$x" | grep "^&[^ ][^ ]*" | sed 's/^&//')
    # If namelist-name recognition succeeded, remember this name as the namelist
    # we're currently in and loop back to consider the next line.
    test -n "$newnl" && nl=$newnl && continue
    # If we're here, the line doesn't look like a namelist name, so see if it
    # contains the key we're looking for. If not, loop back for the next line.
    print "$nl" | grep -qi $2 || continue
    # If we're here, we're in the right namelist...
    print "$x" | grep -qi "^$3 *=" || continue
    # ...and if we're here, we've found the right key. Extract the value from
    # the right of the equals sign.
    value=$(print "$x" | sed "s/^$3 *= *//")
    # We've got the value, so break out of the do loop.
    break
  done
  # Set the specified variable name to the found value.
  eval "$4=\"$value\""
}

function get_nl_value_unquoted
{
  # Get value via get_nl_value, then strip double and single quotes.
  get_nl_value $*
  eval "$4=\$(print \$$4 | tr -d '\"' | tr -d \"'\")"
}

function export_nl
{
  # Export all variables defined in a namelist
  # exports key=value pairs, where value is the whole string after the equals
  # assignment and before any comments
  #
  # usage: export_nl namelist_file namelist
  test -f "$1" || fail "$0: namelist file '$1' not found."
  test -z "$2" && fail "$0: no namelist supplied."
  # AIX sed does not recognize '\t'. Convert tabs to spaces, then consider only
  # spaces as whitespace.
  cat $1 | tr '\t' ' ' | while read line
  do
    # sed: 1. remove comments, 2. remove leading whitespace, 3. remove trailing
    #      whitespace, 4. delete blank lines, 5. delete namelist-termination
    #      lines (those starting with '/'), 6. delete blank lines.
    x=$(print "$line" | sed 's/\!.*$//g;s/^ *//g;s/ *$//g;/^\//d;/^$/d')
    # grep: look for an & followed by one or more non-whitespace characters
    # sed: then, remove the leading &
    newnl=$(print "$x" | grep "^&[^ ][^ ]*" | sed 's/^&//')
    test -n "$newnl" && nl=$newnl && continue
    print "$nl" | grep -qi $2 || continue
    # sed: look for var=[something] and remove any extraneous whitespace
    print "$x" | sed "s/^ *\(.[^ ]*\) *= *\(.[^\!]*\).*$/\1=\2/" | read y
    eval "export $y"
  done
}

function get_fc
{
  # Special case not covered by get_from_nl
  test -z "$SRCDIR" && fail "$0: SRCDIR undefined."
  ./get_buildconfig.ksh $SRCDIR | read FC || fail "$FC"
}

function get_from_nl
{
  # Extract namelist variables using the old Get* functions. Call as
  # 'get_from_nl variable' OR 'get_from_nl function_name as variable', e.g.
  # 'get_from_nl ComputeTasks as PES'.
  test -z "$1" && fail "$0: no variable supplied."
  if [[ -n "$2" && -z "$3" ]]
  then
    fail "Bad call to get_from_nl. Please call using one of the following \
formats: 'get_from_nl GLVL', 'get_from_nl ComputeTasks as PES'."
  elif [[ -z "$2" ]]
  then
    set "$1" "$1" "$1" # set $1, $2, and $3 to the value of $1
  elif [[ "$2" != "as" ]]
  then
    fail "Bad call to get_from_nl. When calling with multiple arguments, must \
be in the format 'get_from_nl function_name as variable', e.g. 'get_from_nl \
ComputeTasks as PES'."
  fi
  print $(cd $fimnamelist_dir && $SRCDIR/bin/Get$1) | read $3 || \
    fail "Get$1 failed (in $SRCDIR/bin/)."
}

function get_srcdir
{
  # If $fimnamelist points to a readable file, try to extract SRCDIR from it
  # using the wrapper binary. If the extracted value isn't a valid directory,
  # suppose that we're running via qsubfim and are a level deeper (in a qsubfim_*
  # directory) than expected and adjust for that. If we *still* don't have a
  # directory, abort. WFM-driven runs are supposed to supply an absolute path
  # for SRCDIR. Things would be simpler if we insisted on absolute paths across
  # the board (TODO?)
  test -r "$fimnamelist" || fail "$0: Cannot read file $fimnamelist"
  fimnamelist_dir=$(print $fimnamelist | sed s:[^/]*$::)
  get_nl_value_unquoted $fimnamelist QUEUEnamelist SRCDIR SRCDIR
  test -d "$SRCDIR" || SRCDIR="$fimnamelist_dir/../$SRCDIR"
  test -d "$SRCDIR" || fail "$0: Cannot set SRCDIR: $SRCDIR not found."
  SRCDIR=$(cd $SRCDIR && print $PWD)
}

function ksh_check
{
  # Check ksh version. If it is ksh93 (regardless of its name), simply return.
  # Otherwise, see if a 'ksh93' binary is available: If so, we can later call
  # ksh_fix() to use ksh93; if not, fail with an informative message.
  test "$(ksh_version)" == "93" && return 0
  test -z $(whence ksh93) && ksh_insist
}

function ksh_fix
{
  # If 'ksh93' is available on the user's path, modify the copied run-automation
  # scripts to use it.
  ksh_check && return
  typeset KSH93 tmp x
  KSH93=$(whence ksh93)
  tmp=sed.tmp
  for x in $(ls)
  do
    grep -q "^#!.*ksh" $x || continue
    sed "s:^#\!.*ksh\(.*\):#\!$KSH93\1:" $x > $tmp || fail
    mv $tmp $x || fail
    chmod +x $x || fail
  done
}

function ksh_insist
{
  if [[ "$(ksh_version)" == "88" ]]
  then
    print "
FIM run automation requires ksh93. This appears to be ksh88. If you are running
via a queue-submission script, it may be sufficient to make 'ksh93' available on
your path (perhaps via a symbolic link) and re-run this script. If you are
calling run-automation scripts (e.g. batchTemplate-prep, via Workflow Manager)
directly, you may need to modify the scripts' initial #! lines to specify the
path to a ksh93 binary.
"
    fail # Comment out to permit ksh88 use, at your own risk.
  fi
}

function ksh_version
{
  # SECONDS contains decimal in ksh93, but an integer in ksh88 and pdksh.
  print $SECONDS | grep -q "\." && print 93 || print 88
}

function linksafe
{
  # Safely link files, throwing an error if the file to be linked is not found
  # It allows both two- and one-argument use. In the latter case, the item gets
  # linked into the current directory with the same name as the target (e.g.
  # "linksafe /etc/passwd" would create a symlink called "passwd" in the current
  # directory).
  # TBH:  This function now creates the link even if the source file is not 
  #       found.  And it attempts retries with delay first.  See comments 
  #       below.  
  test -z "$1" && fail "usage: linksafe target [link]"
  test -z "$2" && link="$PWD" || link="$2"
  typeset msg="Cannot link ($link -> $1)"
  typeset nodename=$(uname -n)
  # Introduce retries with delay here because the workflow manager retries an 
  # entire job step from scratch, re-creating directories.  Thus any delay 
  # due to slow file systems will likely occur again.  Delay here improves 
  # chances of success, although it is not possible to be fault-tolerant in 
  # the presence of file system misbehavior.  
  typeset delay=30    # delay time between retries (seconds)
  typeset attempts_left=6  # number of times to try
  while (( ${attempts_left} > 0 )); do
    (( attempts_left -= 1 ))
    test -e "$1"
    if (( $? != 0 )) ; then
      # For the moment, just warn because the source file may become 
      # visible on the desired node before we need it (parallel file 
      # systems, ick).
      print "WARNING $msg: $1 not visible from $nodename at [$(date)]"
      # fail "$msg: $1 not found"
      if (( ${attempts_left} == 0 )); then
        print "WARNING $msg: Creating dangling link anyway..."
      else
        sleep $delay
      fi
    else
      (( attempts_left = 0 ))
    fi
  done
  ln -s "$1" "$link" > /dev/null 2>&1
}

function set_fimnamelist
{
  fimnamelist="$PWD/FIMnamelist"
}

function trap_off
{
  # Disable error trapping.
  trap - ERR
}

function trap_on
{
  # Enable error trapping. ksh issues the fake signal ERR when any command
  # returns a non-zero status. When this function is called to enable trapping,
  # subsequent commands in the sourcing script that return a non-zero status
  # will trigger the trap. The trap 1) disables xtrace to prevent verbose
  # tracing inside the error handler, and 2) calls the error handler with the
  # line number of the failed command as its argument.
  trap 'set +o xtrace;errhandler $LINENO' ERR
}

function xsource
{
  # Save the current context, source the argument(s), then restore the previous
  # context. Works recursively due to context stack.
  test -z "$CONTEXT" && CONTEXT="unknown"
  print $CONTEXT | grep -q ':' && fail "Colon not allowed in CONTEXT '$CONTEXT'."
  context_push $CONTEXT
  test -z "$XSOURCE_NOTRACE" && set -o xtrace
  . $* || fail "Problem sourcing $1."
  set +o xtrace
  context_pop
  set -o xtrace
}

function xsource_notrace
{
  # Perform xsource but do not enable tracing of sourced script.
  XSOURCE_NOTRACE=1
  xsource $*
  unset XSOURCE_NOTRACE
}

PS4='$CONTEXT[$LINENO]: '

# Turn on error trapping.

trap_on

# Record that this script has been sourced.

functions_sourced="true"

# Turn on xtrace.

set -o xtrace
