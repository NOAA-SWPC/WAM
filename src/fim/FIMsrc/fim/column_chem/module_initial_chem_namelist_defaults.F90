Module module_chem_namelist_defaults
CONTAINS
       subroutine set_chem_namelist_defaults
USE module_initial_chem_namelists
!STARTOFREGISTRYGENERATEDINCLUDE 'inc/namelist_defaults.inc'
!
! THIS does not need to be recopied. Defaults are simple to set!
!
emi_inname = "fimemi"
fireemi_inname = "fimfireemi"
input_chem_inname = "fim_chem_input"
emi_outname  =" "
fireemi_outname = " "
input_chem_outname = " "
frames_per_emissfile = 1
frames_per_fireemissfile = 1
bioemdt = 0
photdt = 0
chemdt = 0
ne_area = 41
kemit = 1
nmegan = 138
kfuture = 1
errosion_dim = 3
chem_conv_tr = 1
chem_opt = 0
gaschem_onoff = 1
aerchem_onoff = 1
wetscav_onoff = 0
cldchem_onoff = 0
vertmix_onoff = 1
chem_in_opt = 0
phot_opt = 0
drydep_opt = 0
emiss_opt = 4
dust_opt = 0
dmsemis_opt = 0
seas_opt = 0
bio_emiss_opt = 0
biomass_burn_opt = 0
plumerisefire_frq = 180
emiss_inpt_opt = 1
gas_bc_opt = 1
gas_ic_opt = 1 
aer_bc_opt = 1 
aer_ic_opt = 1
have_bcs_chem = .false.
aer_ra_feedback = 0
aer_op_opt = 1
END SUBROUTINE set_chem_namelist_defaults

SUBROUTINE set_species
USE module_chemvars
! TBH:  Ignore these so PPP doesn't have to translate them
!SMS$IGNORE BEGIN
USE module_initial_chem_namelists
USE module_data_gocart_dust
USE module_data_gocart_seas
!SMS$IGNORE END
USE module_wrf_control, only: num_chem,num_emis_ant,num_emis_vol
!
if(aer_ra_feedback == 1)then
  P_extcof3 = 1
  P_extcof55 = 2
  P_extcof106 = 3
  P_extcof3_5 = 4
  P_extcof8_12 = 5
  P_bscof3 = 1 
  P_bscof55 = 2 
  P_bscof106 = 3
  P_asympar3 = 1
  P_asympar55 = 2
  P_asympar106 = 3
endif
! gocart fim light
if(chem_opt.eq.304)then
if(num_chem.ne.13)then
   write(6,*) ' num_chem is not equal 13 for gocart fimlight '
   stop
endif
if(num_emis_ant.lt.4)then
   write(6,*) ' num_emis_ant smaller than 4 '
   stop
   if(num_emis_ant.lt.6 .and. biomass_burn_opt.eq.1)then
      write(6,*) ' num_emis_ant smaller than 6 '
      stop
   endif
endif
p_qv=1
p_qc=2
p_qi=3
ch_dust(:,:)=0.8D-9
ch_ss(:,:)=1.
p_so2=1
numgas=4
p_sulf=2
p_dms=3
p_msa=4
p_p25=5
p_bc1=6
p_bc2=7
p_oc1=8
p_oc2=9
p_dust_1=10
p_dust_2=11
p_seas_1=12
p_seas_2=13
p_e_bc  =1
p_e_oc  =2
p_e_sulf=3
p_e_pm_25=4
p_e_so2=5
p_e_pm_10=6
! diagnostic dust and seasale stuff
p_edust1=1
p_edust2=2
p_edust3=3
p_edust4=4
p_edust5=5
p_eseas1=1
p_eseas2=2
p_eseas3=3
p_eseas4=4
endif
!
! 2 tracers
!
if(chem_opt.eq.500)then
if(num_chem.ne.2)then
   write(6,*) ' num_chem is not equal 2 '
   stop
endif
if(num_emis_ant.lt.2)then
   write(6,*) ' num_emis_ant smaller than 2 '
   stop
endif
p_qv=1
p_qc=2
p_qi=3
p_tr1=1
p_tr2=2
p_e_tr1=1
p_e_tr2=2
endif

!
! volcanic ash
!
if(chem_opt.eq.16)then
if(num_chem.ne.10)then
   write(6,*) ' num_chem is not equal 10 for Volcano run'
   stop
endif
if(num_emis_vol.lt.10)then
   write(6,*) ' num_emis_ant smaller than 10 '
   stop
endif
p_qv=1
p_qc=2
p_qi=3
emiss_opt=7
p_vash_1 = 1
p_vash_2 = 2
p_vash_3 = 3
p_vash_4 = 4
p_vash_5 = 5
p_vash_6 = 6
p_vash_7 = 7
p_vash_8 = 8
p_vash_9 = 9
p_vash_10 = 10
p_e_vash1 = 1
p_e_vash2 = 2
p_e_vash3 = 3
p_e_vash4 = 4
p_e_vash5 = 5
p_e_vash6 = 6
p_e_vash7 = 7
p_e_vash8 = 8
p_e_vash9 = 9
p_e_vash10 = 10
numgas=0
endif
! volcanoc ash (4 bins) only
!
if(chem_opt.eq.502)then
print *,'INITIALIZE FOR CHEM_OPT=502'
if(num_emis_vol.lt.4)then
   write(6,*) ' num_emis_vol smaller than 4 '
   stop
endif
if(num_chem.ne.4)then
   write(6,*) ' num_chem is not equal 4 '
   stop
endif
ch_dust(:,:)=0.8D-9
ch_ss(:,:)=1.
p_qv=1
p_qc=2
p_qi=3
numgas=0
p_vash_1 = 1
p_vash_2 = 2
p_vash_3 = 3
p_vash_4 = 4
p_e_vash1 = 1
p_e_vash2 = 2
p_e_vash3 = 3
p_e_vash4 = 4
endif ! chem_opt=502
! gocart simple +volcanic ash simple
if(chem_opt.eq.317)then
print *,'INITIALIZE FOR CHEM_OPT=317'
if(num_emis_vol.lt.4)then
   write(6,*) ' num_emis_vol smaller than 4 '
   stop
endif
if(num_chem.ne.17)then
   write(6,*) ' num_chem is not equal 17 '
   stop
endif
if(num_emis_ant.lt.4)then
   write(6,*) ' num_emis_ant smaller than 4 '
   stop
   if(num_emis_ant.lt.6 .and. biomass_burn_opt.eq.1)then
      write(6,*) ' num_emis_ant smaller than 6 '
      stop
   endif
endif
ch_dust(:,:)=0.8D-9
ch_ss(:,:)=1.
p_qv=1
p_qc=2
p_qi=3
p_so2=1
numgas=4
p_sulf=2
p_dms=3
p_msa=4
p_p25=5
p_bc1=6
p_bc2=7
p_oc1=8
p_oc2=9
p_dust_1=10
p_dust_2=11
p_seas_1=12
p_seas_2=13
p_e_bc  =1
p_e_oc  =2
p_e_sulf=3
p_e_pm_25=4
p_e_so2=5
p_e_pm_10=6
! diagnostic dust and seasale stuff
p_edust1=1
p_edust2=2
p_edust3=3
p_edust4=4
p_edust5=5
p_eseas1=1
p_eseas2=2
p_eseas3=3
p_eseas4=4
p_vash_1 = 14
p_vash_2 = 15
p_vash_3 = 16
p_vash_4 = 17
p_e_vash1 = 1
p_e_vash2 = 2
p_e_vash3 = 3
p_e_vash4 = 4
endif
! gocart simple +volcanic ash
if(chem_opt.eq.316)then
print *,'INITIALIZE FOR CHEM_OPT=316'
if(num_emis_vol.lt.10)then
   write(6,*) ' num_emis_vol smaller than 10 '
   stop
endif
if(num_chem.ne.23)then
   write(6,*) ' num_chem is not equal 18 '
   stop
endif
if(num_emis_ant.lt.4)then
   write(6,*) ' num_emis_ant smaller than 4 '
   stop
   if(num_emis_ant.lt.6 .and. biomass_burn_opt.eq.1)then
      write(6,*) ' num_emis_ant smaller than 6 '
      stop
   endif
endif
ch_dust(:,:)=0.8D-9
ch_ss(:,:)=1.
p_qv=1
p_qc=2
p_qi=3
p_so2=1
numgas=4
p_sulf=2
p_dms=3
p_msa=4
p_p25=5
p_bc1=6
p_bc2=7
p_oc1=8
p_oc2=9
p_dust_1=10
p_dust_2=11
p_seas_1=12
p_seas_2=13
p_e_bc  =1
p_e_oc  =2
p_e_sulf=3
p_e_pm_25=4
p_e_so2=5
p_e_pm_10=6
! diagnostic dust and seasale stuff
p_edust1=1
p_edust2=2
p_edust3=3
p_edust4=4
p_edust5=5
p_eseas1=1
p_eseas2=2
p_eseas3=3
p_eseas4=4
p_vash_1 = 14
p_vash_2 = 15
p_vash_3 = 16
p_vash_4 = 17
p_vash_5 = 18
p_vash_6 = 19
p_vash_7 = 20
p_vash_8 = 21
p_vash_9 = 22
p_vash_10 =23
p_e_vash1 = 1
p_e_vash2 = 2
p_e_vash3 = 3
p_e_vash4 = 4
p_e_vash5 = 5
p_e_vash6 = 6
p_e_vash7 = 7
p_e_vash8 = 8
p_e_vash9 = 9
p_e_vash10 = 10
endif
! gocart simple
if(chem_opt.eq.300)then
if(num_chem.ne.19)then
   write(6,*) ' num_chem is not equal 19 '
   stop
endif
if(num_emis_ant.lt.4)then
   write(6,*) ' num_emis_ant smaller than 4 '
   stop
   if(num_emis_ant.lt.6 .and. biomass_burn_opt.eq.1)then
      write(6,*) ' num_emis_ant smaller than 6 '
      stop
   endif
endif
ch_dust(:,:)=0.8D-9
ch_ss(:,:)=1.
p_qv=1
p_qc=2
p_qi=3
p_so2=1
numgas=4
p_sulf=2
p_dms=3
p_msa=4
p_p25=5
p_bc1=6
p_bc2=7
p_oc1=8
p_oc2=9
p_dust_1=10
p_dust_2=11
p_dust_3=12
p_dust_4=13
p_dust_5=14
p_seas_1=15
p_seas_2=16
p_seas_3=17
p_seas_4=18
p_p10   =19
p_e_bc  =1
p_e_oc  =2
p_e_sulf=3
p_e_pm_25=4
p_e_so2=5
p_e_pm_10=6
! diagnostic dust and seasale stuff
p_edust1=1
p_edust2=2
p_edust3=3
p_edust4=4
p_edust5=5
p_eseas1=1
p_eseas2=2
p_eseas3=3
p_eseas4=4
endif
END SUBROUTINE set_species
END MODULE module_chem_namelist_defaults
