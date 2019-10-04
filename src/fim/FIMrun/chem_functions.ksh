CONTEXT="chem_functions.ksh"

function chem_on
{
  re="^[ \t]*chem_opt[ \t]*=[ \t]*"
  grep "$re" FIMnamelist | sed "s/$re//;s/[ \t].*$//" | read chem_opt
  test -n "$chem_opt" -a "$chem_opt" -gt 0 && true || false
}

function get_chem_opt_value
{
  re="^[ \t]*chem_opt[ \t]*=[ \t]*"
  grep "$re" FIMnamelist | sed "s/$re//;s/[ \t].*$//"
}

function get_chem_in_opt_value
{
  re="^[ \t]*chem_in_opt[ \t]*=[ \t]*"
  grep "$re" FIMnamelist | sed "s/$re//;s/[ \t].*$//"
}

function get_kemit_value
{
  re="^[ \t]*kemit[ \t]*=[ \t]*"
  grep "$re" FIMnamelist | sed "s/$re//;s/[ \t].*$//"
}

function test_suite
{
  print $PWD | grep -q "/FIMtest/"
}

function get_fcst_len
{
  re="^[ \t]*TotalTime[ \t]*=[ \t]*"
  grep "$re" FIMnamelist | sed "s/$re//;s/[ \t].*$//"
}

function chem_fim_setup
{
  print  "IN chem_fim_setup"
  CO2_DATA_DIR="/mnt/lfs0/projects/co2/andy/FIM/CT1x1_to_G7/2011_fluxes_G7/binary"
  CO2_FILE_ROOT="co2_flux."

  # get current run time
  get_runtime
  print "in chem_fim_setup: CHEMFLAG: $CHEMFLAG"
  if [[ $CHEMFLAG == "true" ]]
  then
    test -z "$CHEM_DATADIR" && get_nl_value_unquoted $NLFILE QUEUEnamelist chem_datadir CHEM_DATADIR
    # Link "chem" files
    for x in e_bc.dat e_oc.dat e_pm_10.dat e_pm_25.dat e_so2.dat e_sulf.dat \
      erod1.dat erod2.dat erod3.dat prep_chem_sources prep_chem_sources_template.inp clay.dat sand.dat
    do
      linksafe $CHEM_DATADIR/$x
    done
    test_suite && linksafe $CHEM_DATADIR/volcanic.dat

    # link files for co2
    chem_opt_value=$(get_chem_opt_value)
    if [[ $chem_opt_value == "500" ]] 
    then
      kemit=$(get_kemit_value)
      print "kemit:  $kemit"
      if [[ $kemit == 1 ]]
      then
        fcst_len=$(get_fcst_len)  # ASSUMING THIS IS IN HOURS!
        print "fcst_len:  $fcst_len"
        fcst=0
        while [[ $fcst -le $fcst_len ]] 
        do
          fileDate=$(date +%Y%m%d%H -u -d "$mm/$dd/$yr $hh:00 $fcst hours")
          print "fileDate: $fileDate"
          fileName="${fileDate}.bin" 
          print "cmd:  linksafe $CO2_DATA_DIR/$CO2_FILE_ROOT$fileName"
          linksafe $CO2_DATA_DIR/$CO2_FILE_ROOT$fileName
          fcst=$(expr ${fcst} + 3) 
        done
      fi
    fi 
    # Link "prep" files
    for x in dm0.dat h2o2.dat \
      icos_grid_info_level.dat icos_grid_level.dat no3.dat oh.dat
    do
      linksafe $PREP/$x
    done
    chem_opt_value=$(get_chem_opt_value)
    # get current run time
    get_runtime

    # get date of initialization files
    get_last_run_date
    emiss_date="2009-08-28-00" # default value for branch testing
    test -n "$WFM" && emiss_date="$yr-$mm-$dd-$hh"
    print "emiss_date: $emiss_date"
    print "yr: $yr mm: $mm dd: $dd hh: $hh"
    print "*** init_date_dir: $init_date_dir fcst: $fcst"
    # put date in input file
    sed_safe_chem_datadir=$(print $CHEM_DATADIR | sed 's/\//\\\//g')
    sed "s/\(ihour=\),/\1${hh},/;
         s/\(iday=\),/\1${dd},/;
         s/\(imon=\),/\1${mm},/;
         s/\(iyear=\),/\1${yr},/;
         s/__CHEM_DATADIR__/$sed_safe_chem_datadir/g" \
           prep_chem_sources_template.inp > prep_chem_sources.inp
    ./prep_chem_sources || fail "ERROR: prep_chem_sources failed."
    linksafe FIM-T-${emiss_date}0000-plume.bin plumestuff.dat
    linksafe FIM-T-${emiss_date}0000-OC-bb.bin ebu_oc.dat
    linksafe FIM-T-${emiss_date}0000-BC-bb.bin ebu_bc.dat
    linksafe FIM-T-${emiss_date}0000-PM25-bb.bin ebu_pm25.dat
    linksafe FIM-T-${emiss_date}0000-PM10-bb.bin ebu_pm10.dat
    linksafe FIM-T-${emiss_date}0000-SO2-bb.bin ebu_so2.dat
    linksafe FIM-T-${emiss_date}0000-SULF-bb.bin ebu_sulf.dat
    if [[ $test_suite -eq 0 ]]
    then
      if [[ -s "FIM-T-${emiss_date}0000-g1-volc.bin" ]]
      then
        linksafe FIM-T-${emiss_date}0000-g1-volc.bin volcanic.dat
      fi
    fi

    print "chem_opt_value: $chem_opt_value"
    chem_in_opt_value=$(get_chem_in_opt_value)
    print "chem_in_opt_value: $chem_in_opt_value"
    if [[ $chem_in_opt_value == "1" ]]
    then
    if [[ -n "$WFM" && $init_date_dir != "NOT FOUND" ]] 
    then
      if [[ $chem_opt_value == "317" ]] 
      then
        linksafe ${init_date_dir}/fim_out_ash1000${fcst}${ARCHVTIMEUNIT} vash1.in
        linksafe ${init_date_dir}/fim_out_ash2000${fcst}${ARCHVTIMEUNIT} vash2.in
        linksafe ${init_date_dir}/fim_out_ash3000${fcst}${ARCHVTIMEUNIT} vash3.in
        linksafe ${init_date_dir}/fim_out_ash4000${fcst}${ARCHVTIMEUNIT} vash4.in
      fi
      if [[ $chem_opt_value == "502" ]] 
      then
        linksafe ${init_date_dir}/fim_out_ash1000${fcst}${ARCHVTIMEUNIT} vash1.in
        linksafe ${init_date_dir}/fim_out_ash2000${fcst}${ARCHVTIMEUNIT} vash2.in
        linksafe ${init_date_dir}/fim_out_ash3000${fcst}${ARCHVTIMEUNIT} vash3.in
        linksafe ${init_date_dir}/fim_out_ash4000${fcst}${ARCHVTIMEUNIT} vash4.in
      fi
      if [[ $chem_opt_value -ge "300"  && $chem_opt_value -lt "400" ]] 
      then
        linksafe ${init_date_dir}/fim_out_s4ea000${fcst}${ARCHVTIMEUNIT} seas4.in
        linksafe ${init_date_dir}/fim_out_s3ea000${fcst}${ARCHVTIMEUNIT} seas3.in
        linksafe ${init_date_dir}/fim_out_s2ea000${fcst}${ARCHVTIMEUNIT} seas2.in
        linksafe ${init_date_dir}/fim_out_s1ea000${fcst}${ARCHVTIMEUNIT} seas1.in
        linksafe ${init_date_dir}/fim_out_sulf000${fcst}${ARCHVTIMEUNIT} sulf.in
        linksafe ${init_date_dir}/fim_out_pso2000${fcst}${ARCHVTIMEUNIT} so2.in
        linksafe ${init_date_dir}/fim_out_pbc2000${fcst}${ARCHVTIMEUNIT} bc2.in
        linksafe ${init_date_dir}/fim_out_pbc1000${fcst}${ARCHVTIMEUNIT} bc1.in
        linksafe ${init_date_dir}/fim_out_obc2000${fcst}${ARCHVTIMEUNIT} oc2.in
        linksafe ${init_date_dir}/fim_out_obc1000${fcst}${ARCHVTIMEUNIT} oc1.in
        linksafe ${init_date_dir}/fim_out_d5st000${fcst}${ARCHVTIMEUNIT} dust5.in
        linksafe ${init_date_dir}/fim_out_d4st000${fcst}${ARCHVTIMEUNIT} dust4.in
        linksafe ${init_date_dir}/fim_out_d3st000${fcst}${ARCHVTIMEUNIT} dust3.in
        linksafe ${init_date_dir}/fim_out_d2st000${fcst}${ARCHVTIMEUNIT} dust2.in
        linksafe ${init_date_dir}/fim_out_d1st000${fcst}${ARCHVTIMEUNIT} dust1.in
        linksafe ${init_date_dir}/fim_out_pp25000${fcst}${ARCHVTIMEUNIT} p25.in
        linksafe ${init_date_dir}/fim_out_pp10000${fcst}${ARCHVTIMEUNIT} p10.in
        linksafe ${init_date_dir}/fim_out_dms1000${fcst}${ARCHVTIMEUNIT} dms.in
        linksafe ${init_date_dir}/fim_out_pmsa000${fcst}${ARCHVTIMEUNIT} msa.in
      fi
      if [[ $chem_opt_value == "500" ]] 
      then
        linksafe ${init_date_dir}/fim_out_c13D000${fcst}${ARCHVTIMEUNIT} tr1.in
        linksafe ${init_date_dir}/fim_out_c23D000${fcst}${ARCHVTIMEUNIT} tr2.in
      fi
    fi
  fi
  fi
}

function get_runtime
{
  yr=$(expr substr $yyyymmddhhmm 1 4)
  mm=$(expr substr $yyyymmddhhmm 5 2)
  dd=$(expr substr $yyyymmddhhmm 7 2)
  hh=$(expr substr $yyyymmddhhmm 9 2)
}

function get_last_run_date
{
  get_runtime
  print "in get_last_run_date: yr: $yr mm: $mm dd: $dd hh: $hh"
  typeset -Z3 tmp_fcst
  tmp_fcst=0
  found=0
  while [[ $found -eq 0 ]]
  do
    tmp_fcst=$(expr ${tmp_fcst} + 12) 
    dirDate=$(date +%Y%m%d%H -u -d "$mm/$dd/$yr $hh:00 $tmp_fcst hours ago")
    dir=$FIM_HOME/FIMrun/fim_${GLVL}_${NVL}_${PES}_${dirDate}00/fim_${MEMBER_ID}
    print "checking: ${dir}/fim_out_dms1000${tmp_fcst}${ARCHVTIMEUNIT}"
    if [[ -s ${dir}/fim_out_hgtP000${tmp_fcst}${ARCHVTIMEUNIT} ]]
    then
      found=1
      init_date_dir=$dir
      print "FOUND init_date_dir: $init_date_dir"
    fi
    fcst=$tmp_fcst
    if [[ $fcst -gt 120 ]]
    then
      print "ERROR "
      init_date_dir="NOT FOUND"
      found=1
    fi
  done
}

function chem_prep_newname
{
  test $CHEMFLAG == "true" && (./newname.exe || fail "newname failed")
}

function chem_prep_setup
{
  if [[ $CHEMFLAG == "true" ]]
  then
    test -z "$CHEM_DATADIR" && get_nl_value_unquoted $NLFILE QUEUEnamelist chem_datadir CHEM_DATADIR
    test_suite && linksafe $CHEM_DATADIR/volcanic.dat
    for x in erod_binary gocart_backgd_littlee dm0_binary anthro_binary \
      chemltln.dat
    do
      linksafe $CHEM_DATADIR/$x
    done
    linksafe $BINDIR/newname.exe
  fi
}
