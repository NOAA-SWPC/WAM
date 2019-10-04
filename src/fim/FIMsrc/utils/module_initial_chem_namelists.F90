MODULE module_initial_chem_namelists
! this includes part of state_struct.inc V3.1, all of module_state_description, all of
! namelist_defines2 and a small part (format was changed) of namelist_statements.inc
! The NUM_xxx decs (num_moist, num_chem, num_emis_ant,....to find the important ones search for ggnum)
! from state_description are commented for fim,
! since they have to be defined in module_control
!
! first from state_struct.inc
! GG: updated 11Mar09 to V3.1
integer                                  :: ktauc
integer                                  :: last_chem_time_year
integer                                  :: last_chem_time_month
integer                                  :: last_chem_time_day
integer                                  :: last_chem_time_hour
integer                                  :: last_chem_time_minute
integer                                  :: last_chem_time_second
integer                                  :: emissframes
integer                                  :: fireemissframes
integer                                  :: stepave_count
integer                                  :: stepbioe
integer                                  :: stepphot
integer                                  :: stepchem
integer                                  :: stepfirepl
character*256                               :: emi_inname
character*256                               :: fireemi_inname
character*256                               :: input_chem_inname
character*256                               :: emi_outname
character*256                               :: fireemi_outname
character*256                               :: input_chem_outname
integer                                  :: frames_per_emissfile
integer                                  :: frames_per_fireemissfile
integer                                  :: io_style_emissions
integer                                  :: io_form_emissions
integer                                  :: io_style_fireemissions
integer                                  :: io_form_fireemissions
real                                     :: bioemdt
real                                     :: photdt
real                                     :: chemdt
integer                                  :: ne_area
integer                                  :: kemit
integer                                  :: nmegan
integer                                  :: kfuture
integer                                  :: errosion_dim
integer                                  :: biomass_emiss_opt
integer                                  :: chem_conv_tr
integer                                  :: chem_opt
integer                                  :: gaschem_onoff
integer                                  :: aerchem_onoff
integer                                  :: wetscav_onoff
integer                                  :: cldchem_onoff
integer                                  :: vertmix_onoff
integer                                  :: chem_in_opt
integer                                  :: phot_opt
integer                                  :: drydep_opt
integer                                  :: emiss_opt
integer                                  :: dust_opt
integer                                  :: dmsemis_opt
integer                                  :: seas_opt
integer                                  :: bio_emiss_opt
integer                                  :: biomass_burn_opt
integer                                  :: plumerisefire_frq
integer                                  :: emiss_inpt_opt
integer                                  :: gas_bc_opt
integer                                  :: gas_ic_opt
integer                                  :: aer_bc_opt
integer                                  :: aer_ic_opt
logical                                  :: have_bcs_chem
integer                                  :: aer_ra_feedback
integer                                  :: aer_op_opt
integer                                  :: scalar_opt
! next is wrfphys
!
integer                                  :: mp_physics
integer                                  :: gsfcgce_hail
integer                                  :: gsfcgce_2ice
integer                                  :: progn
integer                                  :: ra_lw_physics
integer                                  :: ra_sw_physics
real                                     :: radt
real                                     :: naer
integer                                  :: sf_sfclay_physics
integer                                  :: sf_surface_physics
integer                                  :: bl_pbl_physics
integer                                  :: sf_urban_physics
real                                     :: bldt
integer                                  :: cu_physics
real                                     :: cudt
real                                     :: gsmdt
integer                                  :: isfflx
integer                                  :: ifsnow
integer                                  :: icloud
real                                     :: swrad_scat
integer                                  :: surface_input_source
integer                                  :: num_urban_layers
integer                                  :: num_months
integer                                  :: maxiens
integer                                  :: maxens
integer                                  :: maxens2
integer                                  :: maxens3
integer                                  :: ensdim
integer                                  :: cugd_avedx
integer                                  :: imomentum
integer                                  :: clos_choice
integer                                  :: num_land_cat
integer                                  :: num_soil_cat
integer                                  :: mp_zero_out
real                                     :: mp_zero_out_thresh
real                                     :: seaice_threshold
integer                                  :: sst_update
integer                                  :: sst_skin
integer                                  :: tmn_update
logical                                  :: usemonalb
logical                                  :: rdmaxalb
logical                                  :: rdlai2d
integer                                  :: co2tf
integer                                  :: ra_call_offset
real                                     :: cam_abs_freq_s
integer                                  :: levsiz
integer                                  :: paerlev
integer                                  :: cam_abs_dim1
integer                                  :: cam_abs_dim2
integer                                  :: lagday
logical                                  :: cu_rad_feedback
integer                                  :: pxlsm_smois_init
integer                                  :: omlcall
real                                     :: oml_hml0
real                                     :: oml_gamma
integer                                  :: isftcflx
real                                     :: shadlen
integer                                  :: slope_rad
integer                                  :: topo_shading
integer                                  :: no_mp_heating
integer                                  :: fractional_seaice
real                                     :: bucket_mm
real                                     :: bucket_j
integer                                  :: grav_settling
!
! package definitions next: they are defined in module_state_description
! updated 10 Mar09 to V3.1
  ! package constants

!STARTOFREGISTRYGENERATEDINCLUDE 'frame/module_state_description.F'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
  ! package constants
  INTEGER, PARAMETER :: prescribe_aerosol = 0
  INTEGER, PARAMETER :: radm2 = 1
  INTEGER, PARAMETER :: radm2sorg = 2
  INTEGER, PARAMETER :: cbmz = 5
  INTEGER, PARAMETER :: cbmz_bb = 6
  INTEGER, PARAMETER :: cbmz_mosaic_4bin = 7
  INTEGER, PARAMETER :: cbmz_mosaic_8bin = 8
  INTEGER, PARAMETER :: cbmz_mosaic_4bin_aq = 9
  INTEGER, PARAMETER :: cbmz_mosaic_8bin_aq = 10
  INTEGER, PARAMETER :: radm2sorg_aq = 11
  INTEGER, PARAMETER :: racmsorg_aq = 12
  INTEGER, PARAMETER :: chem_tracer = 13
  INTEGER, PARAMETER :: chem_trace2 = 14
  INTEGER, PARAMETER :: chem_trace_ens = 15
  INTEGER, PARAMETER :: chem_vash = 16
  INTEGER, PARAMETER :: radm2_kpp = 101
  INTEGER, PARAMETER :: racm_mim_kpp = 102
  INTEGER, PARAMETER :: racm_kpp = 103
  INTEGER, PARAMETER :: racmpm_kpp = 104
  INTEGER, PARAMETER :: racmsorg_kpp = 105
  INTEGER, PARAMETER :: radm2sorg_kpp = 106
  INTEGER, PARAMETER :: cbm4_kpp = 110
  INTEGER, PARAMETER :: nmhc9_kpp = 200
  INTEGER, PARAMETER :: gocart_simple = 300
  INTEGER, PARAMETER :: gocartracm_kpp = 301
  INTEGER, PARAMETER :: gocartradm2_kpp = 302
  INTEGER, PARAMETER :: gocartradm2 = 303
  INTEGER, PARAMETER :: gocartfim = 304
  INTEGER, PARAMETER :: eradm = 2
  INTEGER, PARAMETER :: eradmsorg = 3
  INTEGER, PARAMETER :: ecbmz_mosaic = 4
  INTEGER, PARAMETER :: ecptec = 5
  INTEGER, PARAMETER :: gocart_ecptec = 6
  INTEGER, PARAMETER :: vash = 7
  INTEGER, PARAMETER :: photmad = 1
  INTEGER, PARAMETER :: photfastj = 2
  INTEGER, PARAMETER :: ftuv = 3
  INTEGER, PARAMETER :: wesely = 1
  INTEGER, PARAMETER :: gunther1 = 1
  INTEGER, PARAMETER :: beis313 = 2
  INTEGER, PARAMETER :: megan2 = 3
  INTEGER, PARAMETER :: biomassb = 1
  INTEGER, PARAMETER :: dustgocart = 1
  INTEGER, PARAMETER :: seasgocart = 1
  INTEGER, PARAMETER :: dmsgocart = 1
  INTEGER, PARAMETER :: volume_approx = 1
  INTEGER, PARAMETER :: maxwell_approx = 2
  INTEGER, PARAMETER :: volume_exact = 3
  INTEGER, PARAMETER :: maxwell_exact = 4
  INTEGER, PARAMETER :: shell_exact = 5
  INTEGER, PARAMETER :: emiss_inpt_default = 1
  INTEGER, PARAMETER :: emiss_inpt_cptec = 3
  INTEGER, PARAMETER :: emiss_inpt_pnnl_cm = 101
  INTEGER, PARAMETER :: emiss_inpt_pnnl_rs = 102
  INTEGER, PARAMETER :: emiss_inpt_cb4 = 103
  INTEGER, PARAMETER :: gas_bc_default = 1
  INTEGER, PARAMETER :: gas_bc_pnnl = 101
  INTEGER, PARAMETER :: gas_bc_cbm4 = 102
  INTEGER, PARAMETER :: gas_ic_default = 1
  INTEGER, PARAMETER :: gas_ic_pnnl = 101
  INTEGER, PARAMETER :: gas_ic_cbm4 = 102
  INTEGER, PARAMETER :: aer_bc_default = 1
  INTEGER, PARAMETER :: aer_bc_pnnl = 101
  INTEGER, PARAMETER :: aer_ic_default = 1
  INTEGER, PARAMETER :: aer_ic_pnnl = 101
  INTEGER, PARAMETER :: scalar_me = 1
  INTEGER, PARAMETER :: scalar_sus = 2
  INTEGER, PARAMETER :: passiveqv = 0
  INTEGER, PARAMETER :: kesslerscheme = 1
  INTEGER, PARAMETER :: linscheme = 2
  INTEGER, PARAMETER :: wsm3scheme = 3
  INTEGER, PARAMETER :: wsm5scheme = 4
  INTEGER, PARAMETER :: etampnew = 5
  INTEGER, PARAMETER :: wsm6scheme = 6
  INTEGER, PARAMETER :: gsfcgcescheme = 7
  INTEGER, PARAMETER :: thompson = 8
  INTEGER, PARAMETER :: morr_two_moment = 10
  INTEGER, PARAMETER :: wdm5scheme = 14
  INTEGER, PARAMETER :: wdm6scheme = 16
  INTEGER, PARAMETER :: thompson07 = 98
  INTEGER, PARAMETER :: passiveqv_dfi = 0
  INTEGER, PARAMETER :: kesslerscheme_dfi = 1
  INTEGER, PARAMETER :: linscheme_dfi = 2
  INTEGER, PARAMETER :: wsm3scheme_dfi = 3
  INTEGER, PARAMETER :: wsm5scheme_dfi = 4
  INTEGER, PARAMETER :: etampnew_dfi = 5
  INTEGER, PARAMETER :: wsm6scheme_dfi = 6
  INTEGER, PARAMETER :: gsfcgcescheme_dfi = 7
  INTEGER, PARAMETER :: thompson_dfi = 8
  INTEGER, PARAMETER :: morr_two_moment_dfi = 10
  INTEGER, PARAMETER :: wdm5scheme_dfi = 14
  INTEGER, PARAMETER :: wdm6scheme_dfi = 16
  INTEGER, PARAMETER :: thompson07_dfi = 98
  INTEGER, PARAMETER :: noprogn = 0
  INTEGER, PARAMETER :: progndrop = 1
  INTEGER, PARAMETER :: rrtmscheme = 1
  INTEGER, PARAMETER :: camlwscheme = 3
  INTEGER, PARAMETER :: rrtmg_lwscheme = 4
  INTEGER, PARAMETER :: gfdllwscheme = 99
  INTEGER, PARAMETER :: heldsuarez = 31
  INTEGER, PARAMETER :: swradscheme = 1
  INTEGER, PARAMETER :: gsfcswscheme = 2
  INTEGER, PARAMETER :: camswscheme = 3
  INTEGER, PARAMETER :: rrtmg_swscheme = 4
  INTEGER, PARAMETER :: gfdlswscheme = 99
  INTEGER, PARAMETER :: sfclayscheme = 1
  INTEGER, PARAMETER :: myjsfcscheme = 2
  INTEGER, PARAMETER :: gfssfcscheme = 3
  INTEGER, PARAMETER :: qnsesfcscheme = 4
  INTEGER, PARAMETER :: mynnsfcscheme = 5
  INTEGER, PARAMETER :: pxsfcscheme = 7
  INTEGER, PARAMETER :: noahucmscheme = 1
  INTEGER, PARAMETER :: bepscheme = 2
  INTEGER, PARAMETER :: slabscheme = 1
  INTEGER, PARAMETER :: lsmscheme = 2
  INTEGER, PARAMETER :: ruclsmscheme = 3
  INTEGER, PARAMETER :: pxlsmscheme = 7
  INTEGER, PARAMETER :: ysuscheme = 1
  INTEGER, PARAMETER :: myjpblscheme = 2
  INTEGER, PARAMETER :: gfsscheme = 3
  INTEGER, PARAMETER :: qnsepblscheme = 4
  INTEGER, PARAMETER :: mynnpblscheme2 = 5
  INTEGER, PARAMETER :: mynnpblscheme3 = 6
  INTEGER, PARAMETER :: acmpblscheme = 7
  INTEGER, PARAMETER :: boulacscheme = 8
  INTEGER, PARAMETER :: mrfscheme = 99
  INTEGER, PARAMETER :: kfetascheme = 1
  INTEGER, PARAMETER :: bmjscheme = 2
  INTEGER, PARAMETER :: gdscheme = 3
  INTEGER, PARAMETER :: sasscheme = 4
  INTEGER, PARAMETER :: g3scheme = 5
  INTEGER, PARAMETER :: kfscheme = 99
  INTEGER, PARAMETER :: psufddagd = 1
  INTEGER, PARAMETER :: psusfddagd = 1
  INTEGER, PARAMETER :: spnudging = 2
  INTEGER, PARAMETER :: restofwrf = 0
  INTEGER, PARAMETER :: original = 0
  INTEGER, PARAMETER :: positivedef = 1
  INTEGER, PARAMETER :: monotonic = 2
  INTEGER, PARAMETER :: dfi_setup = 0
  INTEGER, PARAMETER :: dfi_bck = 1
  INTEGER, PARAMETER :: dfi_fwd = 2
  INTEGER, PARAMETER :: dfi_fst = 3
  INTEGER, PARAMETER :: dfi_nodfi = 0
  INTEGER, PARAMETER :: dfi_dfl = 1
  INTEGER, PARAMETER :: dfi_ddfi = 2
  INTEGER, PARAMETER :: dfi_tdfi = 3
  INTEGER, PARAMETER :: realonly = 1
  INTEGER, PARAMETER :: io_intio = 1
  INTEGER, PARAMETER :: io_netcdf = 2
  INTEGER, PARAMETER :: io_hdf = 3
  INTEGER, PARAMETER :: io_phdf5 = 4
  INTEGER, PARAMETER :: io_grib1 = 5
  INTEGER, PARAMETER :: io_mcel = 6
  INTEGER, PARAMETER :: io_esmf = 7
  INTEGER, PARAMETER :: io_yyy = 8
  INTEGER, PARAMETER :: io_zzz = 9
  INTEGER, PARAMETER :: io_grib2 = 10
  INTEGER, PARAMETER :: io_pnetcdf = 11
  INTEGER, PARAMETER :: fire_sfire = 2
  ! 4D array constants
  INTEGER, PARAMETER :: PARAM_qv = 1
  INTEGER            ::     P_qv = 1
  LOGICAL            ::     F_qv = .FALSE.
  INTEGER, PARAMETER :: PARAM_qc = 2
  INTEGER            ::     P_qc = 1
  LOGICAL            ::     F_qc = .FALSE.
  INTEGER, PARAMETER :: PARAM_qr = 3
  INTEGER            ::     P_qr = 1
  LOGICAL            ::     F_qr = .FALSE.
  INTEGER, PARAMETER :: PARAM_qi = 4
  INTEGER            ::     P_qi = 1
  LOGICAL            ::     F_qi = .FALSE.
  INTEGER, PARAMETER :: PARAM_qs = 5
  INTEGER            ::     P_qs = 1
  LOGICAL            ::     F_qs = .FALSE.
  INTEGER, PARAMETER :: PARAM_qg = 6
  INTEGER            ::     P_qg = 1
  LOGICAL            ::     F_qg = .FALSE.
  INTEGER, PARAMETER :: PARAM_NUM_moist = 7
!ggnum INTEGER            ::       NUM_moist = 1
  INTEGER, PARAMETER :: PARAM_dfi_qv = 1
  INTEGER            ::     P_dfi_qv = 1
  LOGICAL            ::     F_dfi_qv = .FALSE.
  INTEGER, PARAMETER :: PARAM_dfi_qc = 2
  INTEGER            ::     P_dfi_qc = 1
  LOGICAL            ::     F_dfi_qc = .FALSE.
  INTEGER, PARAMETER :: PARAM_dfi_qr = 3
  INTEGER            ::     P_dfi_qr = 1
  LOGICAL            ::     F_dfi_qr = .FALSE.
  INTEGER, PARAMETER :: PARAM_dfi_qi = 4
  INTEGER            ::     P_dfi_qi = 1
  LOGICAL            ::     F_dfi_qi = .FALSE.
  INTEGER, PARAMETER :: PARAM_dfi_qs = 5
  INTEGER            ::     P_dfi_qs = 1
  LOGICAL            ::     F_dfi_qs = .FALSE.
  INTEGER, PARAMETER :: PARAM_dfi_qg = 6
  INTEGER            ::     P_dfi_qg = 1
  LOGICAL            ::     F_dfi_qg = .FALSE.
  INTEGER, PARAMETER :: PARAM_NUM_dfi_moist = 7
  INTEGER            ::       NUM_dfi_moist = 1
  INTEGER, PARAMETER :: PARAM_e_iso = 1
  INTEGER            ::     P_e_iso = 1
  LOGICAL            ::     F_e_iso = .FALSE.
  INTEGER, PARAMETER :: PARAM_e_so2 = 2
  INTEGER            ::     P_e_so2 = 1
  LOGICAL            ::     F_e_so2 = .FALSE.
  INTEGER, PARAMETER :: PARAM_e_no = 3
  INTEGER            ::     P_e_no = 1
  LOGICAL            ::     F_e_no = .FALSE.
  INTEGER, PARAMETER :: PARAM_e_co = 4
  INTEGER            ::     P_e_co = 1
  LOGICAL            ::     F_e_co = .FALSE.
  INTEGER, PARAMETER :: PARAM_e_eth = 5
  INTEGER            ::     P_e_eth = 1
  LOGICAL            ::     F_e_eth = .FALSE.
  INTEGER, PARAMETER :: PARAM_e_hc3 = 6
  INTEGER            ::     P_e_hc3 = 1
  LOGICAL            ::     F_e_hc3 = .FALSE.
  INTEGER, PARAMETER :: PARAM_e_hc5 = 7
  INTEGER            ::     P_e_hc5 = 1
  LOGICAL            ::     F_e_hc5 = .FALSE.
  INTEGER, PARAMETER :: PARAM_e_hc8 = 8
  INTEGER            ::     P_e_hc8 = 1
  LOGICAL            ::     F_e_hc8 = .FALSE.
  INTEGER, PARAMETER :: PARAM_e_xyl = 9
  INTEGER            ::     P_e_xyl = 1
  LOGICAL            ::     F_e_xyl = .FALSE.
  INTEGER, PARAMETER :: PARAM_e_ol2 = 10
  INTEGER            ::     P_e_ol2 = 1
  LOGICAL            ::     F_e_ol2 = .FALSE.
  INTEGER, PARAMETER :: PARAM_e_olt = 11
  INTEGER            ::     P_e_olt = 1
  LOGICAL            ::     F_e_olt = .FALSE.
  INTEGER, PARAMETER :: PARAM_e_oli = 12
  INTEGER            ::     P_e_oli = 1
  LOGICAL            ::     F_e_oli = .FALSE.
  INTEGER, PARAMETER :: PARAM_e_tol = 13
  INTEGER            ::     P_e_tol = 1
  LOGICAL            ::     F_e_tol = .FALSE.
  INTEGER, PARAMETER :: PARAM_e_csl = 14
  INTEGER            ::     P_e_csl = 1
  LOGICAL            ::     F_e_csl = .FALSE.
  INTEGER, PARAMETER :: PARAM_e_hcho = 15
  INTEGER            ::     P_e_hcho = 1
  LOGICAL            ::     F_e_hcho = .FALSE.
  INTEGER, PARAMETER :: PARAM_e_ald = 16
  INTEGER            ::     P_e_ald = 1
  LOGICAL            ::     F_e_ald = .FALSE.
  INTEGER, PARAMETER :: PARAM_e_ket = 17
  INTEGER            ::     P_e_ket = 1
  LOGICAL            ::     F_e_ket = .FALSE.
  INTEGER, PARAMETER :: PARAM_e_ora2 = 18
  INTEGER            ::     P_e_ora2 = 1
  LOGICAL            ::     F_e_ora2 = .FALSE.
  INTEGER, PARAMETER :: PARAM_e_nh3 = 19
  INTEGER            ::     P_e_nh3 = 1
  LOGICAL            ::     F_e_nh3 = .FALSE.
  INTEGER, PARAMETER :: PARAM_e_pm_25 = 20
  INTEGER            ::     P_e_pm_25 = 1
  LOGICAL            ::     F_e_pm_25 = .FALSE.
  INTEGER, PARAMETER :: PARAM_e_pm_10 = 21
  INTEGER            ::     P_e_pm_10 = 1
  LOGICAL            ::     F_e_pm_10 = .FALSE.
  INTEGER, PARAMETER :: PARAM_e_pm25i = 22
  INTEGER            ::     P_e_pm25i = 1
  LOGICAL            ::     F_e_pm25i = .FALSE.
  INTEGER, PARAMETER :: PARAM_e_pm25j = 23
  INTEGER            ::     P_e_pm25j = 1
  LOGICAL            ::     F_e_pm25j = .FALSE.
  INTEGER, PARAMETER :: PARAM_e_eci = 24
  INTEGER            ::     P_e_eci = 1
  LOGICAL            ::     F_e_eci = .FALSE.
  INTEGER, PARAMETER :: PARAM_e_ecj = 25
  INTEGER            ::     P_e_ecj = 1
  LOGICAL            ::     F_e_ecj = .FALSE.
  INTEGER, PARAMETER :: PARAM_e_orgi = 26
  INTEGER            ::     P_e_orgi = 1
  LOGICAL            ::     F_e_orgi = .FALSE.
  INTEGER, PARAMETER :: PARAM_e_orgj = 27
  INTEGER            ::     P_e_orgj = 1
  LOGICAL            ::     F_e_orgj = .FALSE.
  INTEGER, PARAMETER :: PARAM_e_so4i = 28
  INTEGER            ::     P_e_so4i = 1
  LOGICAL            ::     F_e_so4i = .FALSE.
  INTEGER, PARAMETER :: PARAM_e_so4j = 29
  INTEGER            ::     P_e_so4j = 1
  LOGICAL            ::     F_e_so4j = .FALSE.
  INTEGER, PARAMETER :: PARAM_e_no3i = 30
  INTEGER            ::     P_e_no3i = 1
  LOGICAL            ::     F_e_no3i = .FALSE.
  INTEGER, PARAMETER :: PARAM_e_no3j = 31
  INTEGER            ::     P_e_no3j = 1
  LOGICAL            ::     F_e_no3j = .FALSE.
  INTEGER, PARAMETER :: PARAM_e_no2 = 32
  INTEGER            ::     P_e_no2 = 1
  LOGICAL            ::     F_e_no2 = .FALSE.
  INTEGER, PARAMETER :: PARAM_e_ch3oh = 33
  INTEGER            ::     P_e_ch3oh = 1
  LOGICAL            ::     F_e_ch3oh = .FALSE.
  INTEGER, PARAMETER :: PARAM_e_c2h5oh = 34
  INTEGER            ::     P_e_c2h5oh = 1
  LOGICAL            ::     F_e_c2h5oh = .FALSE.
  INTEGER, PARAMETER :: PARAM_e_so4c = 35
  INTEGER            ::     P_e_so4c = 1
  LOGICAL            ::     F_e_so4c = .FALSE.
  INTEGER, PARAMETER :: PARAM_e_no3c = 36
  INTEGER            ::     P_e_no3c = 1
  LOGICAL            ::     F_e_no3c = .FALSE.
  INTEGER, PARAMETER :: PARAM_e_orgc = 37
  INTEGER            ::     P_e_orgc = 1
  LOGICAL            ::     F_e_orgc = .FALSE.
  INTEGER, PARAMETER :: PARAM_e_ecc = 38
  INTEGER            ::     P_e_ecc = 1
  LOGICAL            ::     F_e_ecc = .FALSE.
  INTEGER, PARAMETER :: PARAM_e_bc = 39
  INTEGER            ::     P_e_bc = 1
  LOGICAL            ::     F_e_bc = .FALSE.
  INTEGER, PARAMETER :: PARAM_e_oc = 40
  INTEGER            ::     P_e_oc = 1
  LOGICAL            ::     F_e_oc = .FALSE.
  INTEGER, PARAMETER :: PARAM_e_sulf = 41
  INTEGER            ::     P_e_sulf = 1
  LOGICAL            ::     F_e_sulf = .FALSE.
  INTEGER, PARAMETER :: PARAM_e_vash1 = 42
  INTEGER            ::     P_e_vash1 = 1
  LOGICAL            ::     F_e_vash1 = .FALSE.
  INTEGER, PARAMETER :: PARAM_e_vash2 = 43
  INTEGER            ::     P_e_vash2 = 1
  LOGICAL            ::     F_e_vash2 = .FALSE.
  INTEGER, PARAMETER :: PARAM_e_vash3 = 44
  INTEGER            ::     P_e_vash3 = 1
  LOGICAL            ::     F_e_vash3 = .FALSE.
  INTEGER, PARAMETER :: PARAM_e_vash4 = 45
  INTEGER            ::     P_e_vash4 = 1
  LOGICAL            ::     F_e_vash4 = .FALSE.
  INTEGER, PARAMETER :: PARAM_e_vash5 = 46
  INTEGER            ::     P_e_vash5 = 1
  LOGICAL            ::     F_e_vash5 = .FALSE.
  INTEGER, PARAMETER :: PARAM_e_vash6 = 47
  INTEGER            ::     P_e_vash6 = 1
  LOGICAL            ::     F_e_vash6 = .FALSE.
  INTEGER, PARAMETER :: PARAM_e_vash7 = 48
  INTEGER            ::     P_e_vash7 = 1
  LOGICAL            ::     F_e_vash7 = .FALSE.
  INTEGER, PARAMETER :: PARAM_e_vash8 = 49
  INTEGER            ::     P_e_vash8 = 1
  LOGICAL            ::     F_e_vash8 = .FALSE.
  INTEGER, PARAMETER :: PARAM_e_vash9 = 50
  INTEGER            ::     P_e_vash9 = 1
  LOGICAL            ::     F_e_vash9 = .FALSE.
  INTEGER, PARAMETER :: PARAM_e_vash10 = 51
  INTEGER            ::     P_e_vash10 = 1
  LOGICAL            ::     F_e_vash10 = .FALSE.
  INTEGER, PARAMETER :: PARAM_e_tr1 = 52
  INTEGER            ::     P_e_tr1 = 1
  LOGICAL            ::     F_e_tr1 = .FALSE.
  INTEGER, PARAMETER :: PARAM_e_tr2 = 53
  INTEGER            ::     P_e_tr2 = 1
  LOGICAL            ::     F_e_tr2 = .FALSE.
  INTEGER, PARAMETER :: PARAM_NUM_emis_ant = 54

!  INTEGER, PARAMETER :: PARAM_NUM_emis_ant = 52
!ggnum INTEGER            ::       NUM_emis_ant = 1
  INTEGER, PARAMETER :: PARAM_edust1 = 1
  INTEGER            ::     P_edust1 = 1
  LOGICAL            ::     F_edust1 = .FALSE.
  INTEGER, PARAMETER :: PARAM_edust2 = 2
  INTEGER            ::     P_edust2 = 1
  LOGICAL            ::     F_edust2 = .FALSE.
  INTEGER, PARAMETER :: PARAM_edust3 = 3
  INTEGER            ::     P_edust3 = 1
  LOGICAL            ::     F_edust3 = .FALSE.
  INTEGER, PARAMETER :: PARAM_edust4 = 4
  INTEGER            ::     P_edust4 = 1
  LOGICAL            ::     F_edust4 = .FALSE.
  INTEGER, PARAMETER :: PARAM_edust5 = 5
  INTEGER            ::     P_edust5 = 1
  LOGICAL            ::     F_edust5 = .FALSE.
  INTEGER, PARAMETER :: PARAM_NUM_emis_dust = 6
!ggnum INTEGER            ::       NUM_emis_dust = 1
  INTEGER, PARAMETER :: PARAM_eseas1 = 1
  INTEGER            ::     P_eseas1 = 1
  LOGICAL            ::     F_eseas1 = .FALSE.
  INTEGER, PARAMETER :: PARAM_eseas2 = 2
  INTEGER            ::     P_eseas2 = 1
  LOGICAL            ::     F_eseas2 = .FALSE.
  INTEGER, PARAMETER :: PARAM_eseas3 = 3
  INTEGER            ::     P_eseas3 = 1
  LOGICAL            ::     F_eseas3 = .FALSE.
  INTEGER, PARAMETER :: PARAM_eseas4 = 4
  INTEGER            ::     P_eseas4 = 1
  LOGICAL            ::     F_eseas4 = .FALSE.
  INTEGER, PARAMETER :: PARAM_NUM_emis_seas = 5
  INTEGER, PARAMETER :: PARAM_extcof3 = 1
  INTEGER            ::     P_extcof3 = 1 
  LOGICAL            ::     F_extcof3 = .FALSE.
  INTEGER, PARAMETER :: PARAM_extcof55 = 2
  INTEGER            ::     P_extcof55 = 1
  LOGICAL            ::     F_extcof55 = .FALSE.
  INTEGER, PARAMETER :: PARAM_extcof106 = 3 
  INTEGER            ::     P_extcof106 = 1
  LOGICAL            ::     F_extcof106 = .FALSE.
  INTEGER, PARAMETER :: PARAM_extcof3_5 = 4
  INTEGER            ::     P_extcof3_5 = 1
  LOGICAL            ::     F_extcof3_5 = .FALSE.
  INTEGER, PARAMETER :: PARAM_extcof8_12 = 5 
  INTEGER            ::     P_extcof8_12 = 1
  LOGICAL            ::     F_extcof8_12 = .FALSE.
  INTEGER, PARAMETER :: PARAM_NUM_ext_coef = 6
! INTEGER            ::       NUM_ext_coef = 1 
  INTEGER, PARAMETER :: PARAM_bscof3 = 1
  INTEGER            ::     P_bscof3 = 1
  LOGICAL            ::     F_bscof3 = .FALSE.
  INTEGER, PARAMETER :: PARAM_bscof55 = 2
  INTEGER            ::     P_bscof55 = 1
  LOGICAL            ::     F_bscof55 = .FALSE.
  INTEGER, PARAMETER :: PARAM_bscof106 = 3
  INTEGER            ::     P_bscof106 = 1
  LOGICAL            ::     F_bscof106 = .FALSE.
  INTEGER, PARAMETER :: PARAM_NUM_bscat_coef = 4
! INTEGER            ::       NUM_bscat_coef = 1
  INTEGER, PARAMETER :: PARAM_asympar3 = 1
  INTEGER            ::     P_asympar3 = 1
  LOGICAL            ::     F_asympar3 = .FALSE.
  INTEGER, PARAMETER :: PARAM_asympar55 = 2
  INTEGER            ::     P_asympar55 = 1
  LOGICAL            ::     F_asympar55 = .FALSE.
  INTEGER, PARAMETER :: PARAM_asympar106 = 3
  INTEGER            ::     P_asympar106 = 1
  LOGICAL            ::     F_asympar106 = .FALSE.
  INTEGER, PARAMETER :: PARAM_NUM_asym_par = 4
! INTEGER            ::       NUM_asym_par = 1

!ggnum INTEGER            ::       NUM_emis_seas = 1
  INTEGER, PARAMETER :: PARAM_so2 = 1
  INTEGER            ::     P_so2 = 1
  LOGICAL            ::     F_so2 = .FALSE.
  INTEGER, PARAMETER :: PARAM_sulf = 2
  INTEGER            ::     P_sulf = 1
  LOGICAL            ::     F_sulf = .FALSE.
  INTEGER, PARAMETER :: PARAM_no2 = 3
  INTEGER            ::     P_no2 = 1
  LOGICAL            ::     F_no2 = .FALSE.
  INTEGER, PARAMETER :: PARAM_no = 4
  INTEGER            ::     P_no = 1
  LOGICAL            ::     F_no = .FALSE.
  INTEGER, PARAMETER :: PARAM_o3 = 5
  INTEGER            ::     P_o3 = 1
  LOGICAL            ::     F_o3 = .FALSE.
  INTEGER, PARAMETER :: PARAM_hno3 = 6
  INTEGER            ::     P_hno3 = 1
  LOGICAL            ::     F_hno3 = .FALSE.
  INTEGER, PARAMETER :: PARAM_h2o2 = 7
  INTEGER            ::     P_h2o2 = 1
  LOGICAL            ::     F_h2o2 = .FALSE.
  INTEGER, PARAMETER :: PARAM_ald = 8
  INTEGER            ::     P_ald = 1
  LOGICAL            ::     F_ald = .FALSE.
  INTEGER, PARAMETER :: PARAM_hcho = 9
  INTEGER            ::     P_hcho = 1
  LOGICAL            ::     F_hcho = .FALSE.
  INTEGER, PARAMETER :: PARAM_op1 = 10
  INTEGER            ::     P_op1 = 1
  LOGICAL            ::     F_op1 = .FALSE.
  INTEGER, PARAMETER :: PARAM_op2 = 11
  INTEGER            ::     P_op2 = 1
  LOGICAL            ::     F_op2 = .FALSE.
  INTEGER, PARAMETER :: PARAM_paa = 12
  INTEGER            ::     P_paa = 1
  LOGICAL            ::     F_paa = .FALSE.
  INTEGER, PARAMETER :: PARAM_ora1 = 13
  INTEGER            ::     P_ora1 = 1
  LOGICAL            ::     F_ora1 = .FALSE.
  INTEGER, PARAMETER :: PARAM_ora2 = 14
  INTEGER            ::     P_ora2 = 1
  LOGICAL            ::     F_ora2 = .FALSE.
  INTEGER, PARAMETER :: PARAM_nh3 = 15
  INTEGER            ::     P_nh3 = 1
  LOGICAL            ::     F_nh3 = .FALSE.
  INTEGER, PARAMETER :: PARAM_n2o5 = 16
  INTEGER            ::     P_n2o5 = 1
  LOGICAL            ::     F_n2o5 = .FALSE.
  INTEGER, PARAMETER :: PARAM_no3 = 17
  INTEGER            ::     P_no3 = 1
  LOGICAL            ::     F_no3 = .FALSE.
  INTEGER, PARAMETER :: PARAM_pan = 18
  INTEGER            ::     P_pan = 1
  LOGICAL            ::     F_pan = .FALSE.
  INTEGER, PARAMETER :: PARAM_hc3 = 19
  INTEGER            ::     P_hc3 = 1
  LOGICAL            ::     F_hc3 = .FALSE.
  INTEGER, PARAMETER :: PARAM_hc5 = 20
  INTEGER            ::     P_hc5 = 1
  LOGICAL            ::     F_hc5 = .FALSE.
  INTEGER, PARAMETER :: PARAM_hc8 = 21
  INTEGER            ::     P_hc8 = 1
  LOGICAL            ::     F_hc8 = .FALSE.
  INTEGER, PARAMETER :: PARAM_eth = 22
  INTEGER            ::     P_eth = 1
  LOGICAL            ::     F_eth = .FALSE.
  INTEGER, PARAMETER :: PARAM_co = 23
  INTEGER            ::     P_co = 1
  LOGICAL            ::     F_co = .FALSE.
  INTEGER, PARAMETER :: PARAM_ol2 = 24
  INTEGER            ::     P_ol2 = 1
  LOGICAL            ::     F_ol2 = .FALSE.
  INTEGER, PARAMETER :: PARAM_olt = 25
  INTEGER            ::     P_olt = 1
  LOGICAL            ::     F_olt = .FALSE.
  INTEGER, PARAMETER :: PARAM_oli = 26
  INTEGER            ::     P_oli = 1
  LOGICAL            ::     F_oli = .FALSE.
  INTEGER, PARAMETER :: PARAM_tol = 27
  INTEGER            ::     P_tol = 1
  LOGICAL            ::     F_tol = .FALSE.
  INTEGER, PARAMETER :: PARAM_xyl = 28
  INTEGER            ::     P_xyl = 1
  LOGICAL            ::     F_xyl = .FALSE.
  INTEGER, PARAMETER :: PARAM_aco3 = 29
  INTEGER            ::     P_aco3 = 1
  LOGICAL            ::     F_aco3 = .FALSE.
  INTEGER, PARAMETER :: PARAM_tpan = 30
  INTEGER            ::     P_tpan = 1
  LOGICAL            ::     F_tpan = .FALSE.
  INTEGER, PARAMETER :: PARAM_hono = 31
  INTEGER            ::     P_hono = 1
  LOGICAL            ::     F_hono = .FALSE.
  INTEGER, PARAMETER :: PARAM_hno4 = 32
  INTEGER            ::     P_hno4 = 1
  LOGICAL            ::     F_hno4 = .FALSE.
  INTEGER, PARAMETER :: PARAM_ket = 33
  INTEGER            ::     P_ket = 1
  LOGICAL            ::     F_ket = .FALSE.
  INTEGER, PARAMETER :: PARAM_gly = 34
  INTEGER            ::     P_gly = 1
  LOGICAL            ::     F_gly = .FALSE.
  INTEGER, PARAMETER :: PARAM_mgly = 35
  INTEGER            ::     P_mgly = 1
  LOGICAL            ::     F_mgly = .FALSE.
  INTEGER, PARAMETER :: PARAM_dcb = 36
  INTEGER            ::     P_dcb = 1
  LOGICAL            ::     F_dcb = .FALSE.
  INTEGER, PARAMETER :: PARAM_onit = 37
  INTEGER            ::     P_onit = 1
  LOGICAL            ::     F_onit = .FALSE.
  INTEGER, PARAMETER :: PARAM_csl = 38
  INTEGER            ::     P_csl = 1
  LOGICAL            ::     F_csl = .FALSE.
  INTEGER, PARAMETER :: PARAM_iso = 39
  INTEGER            ::     P_iso = 1
  LOGICAL            ::     F_iso = .FALSE.
  INTEGER, PARAMETER :: PARAM_ho = 40
  INTEGER            ::     P_ho = 1
  LOGICAL            ::     F_ho = .FALSE.
  INTEGER, PARAMETER :: PARAM_ho2 = 41
  INTEGER            ::     P_ho2 = 1
  LOGICAL            ::     F_ho2 = .FALSE.
  INTEGER, PARAMETER :: PARAM_ete = 42
  INTEGER            ::     P_ete = 1
  LOGICAL            ::     F_ete = .FALSE.
  INTEGER, PARAMETER :: PARAM_co2 = 43
  INTEGER            ::     P_co2 = 1
  LOGICAL            ::     F_co2 = .FALSE.
  INTEGER, PARAMETER :: PARAM_ch4 = 44
  INTEGER            ::     P_ch4 = 1
  LOGICAL            ::     F_ch4 = .FALSE.
  INTEGER, PARAMETER :: PARAM_udd = 45
  INTEGER            ::     P_udd = 1
  LOGICAL            ::     F_udd = .FALSE.
  INTEGER, PARAMETER :: PARAM_hket = 46
  INTEGER            ::     P_hket = 1
  LOGICAL            ::     F_hket = .FALSE.
  INTEGER, PARAMETER :: PARAM_api = 47
  INTEGER            ::     P_api = 1
  LOGICAL            ::     F_api = .FALSE.
  INTEGER, PARAMETER :: PARAM_lim = 48
  INTEGER            ::     P_lim = 1
  LOGICAL            ::     F_lim = .FALSE.
  INTEGER, PARAMETER :: PARAM_dien = 49
  INTEGER            ::     P_dien = 1
  LOGICAL            ::     F_dien = .FALSE.
  INTEGER, PARAMETER :: PARAM_macr = 50
  INTEGER            ::     P_macr = 1
  LOGICAL            ::     F_macr = .FALSE.
  INTEGER, PARAMETER :: PARAM_hcl = 51
  INTEGER            ::     P_hcl = 1
  LOGICAL            ::     F_hcl = .FALSE.
  INTEGER, PARAMETER :: PARAM_ch3o2 = 52
  INTEGER            ::     P_ch3o2 = 1
  LOGICAL            ::     F_ch3o2 = .FALSE.
  INTEGER, PARAMETER :: PARAM_ethp = 53
  INTEGER            ::     P_ethp = 1
  LOGICAL            ::     F_ethp = .FALSE.
  INTEGER, PARAMETER :: PARAM_ch3oh = 54
  INTEGER            ::     P_ch3oh = 1
  LOGICAL            ::     F_ch3oh = .FALSE.
  INTEGER, PARAMETER :: PARAM_c2h5oh = 55
  INTEGER            ::     P_c2h5oh = 1
  LOGICAL            ::     F_c2h5oh = .FALSE.
  INTEGER, PARAMETER :: PARAM_par = 56
  INTEGER            ::     P_par = 1
  LOGICAL            ::     F_par = .FALSE.
  INTEGER, PARAMETER :: PARAM_to2 = 57
  INTEGER            ::     P_to2 = 1
  LOGICAL            ::     F_to2 = .FALSE.
  INTEGER, PARAMETER :: PARAM_cro = 58
  INTEGER            ::     P_cro = 1
  LOGICAL            ::     F_cro = .FALSE.
  INTEGER, PARAMETER :: PARAM_open = 59
  INTEGER            ::     P_open = 1
  LOGICAL            ::     F_open = .FALSE.
  INTEGER, PARAMETER :: PARAM_op3 = 60
  INTEGER            ::     P_op3 = 1
  LOGICAL            ::     F_op3 = .FALSE.
  INTEGER, PARAMETER :: PARAM_c2o3 = 61
  INTEGER            ::     P_c2o3 = 1
  LOGICAL            ::     F_c2o3 = .FALSE.
  INTEGER, PARAMETER :: PARAM_ro2 = 62
  INTEGER            ::     P_ro2 = 1
  LOGICAL            ::     F_ro2 = .FALSE.
  INTEGER, PARAMETER :: PARAM_ano2 = 63
  INTEGER            ::     P_ano2 = 1
  LOGICAL            ::     F_ano2 = .FALSE.
  INTEGER, PARAMETER :: PARAM_nap = 64
  INTEGER            ::     P_nap = 1
  LOGICAL            ::     F_nap = .FALSE.
  INTEGER, PARAMETER :: PARAM_xo2 = 65
  INTEGER            ::     P_xo2 = 1
  LOGICAL            ::     F_xo2 = .FALSE.
  INTEGER, PARAMETER :: PARAM_xpar = 66
  INTEGER            ::     P_xpar = 1
  LOGICAL            ::     F_xpar = .FALSE.
  INTEGER, PARAMETER :: PARAM_isoprd = 67
  INTEGER            ::     P_isoprd = 1
  LOGICAL            ::     F_isoprd = .FALSE.
  INTEGER, PARAMETER :: PARAM_isopp = 68
  INTEGER            ::     P_isopp = 1
  LOGICAL            ::     F_isopp = .FALSE.
  INTEGER, PARAMETER :: PARAM_isopn = 69
  INTEGER            ::     P_isopn = 1
  LOGICAL            ::     F_isopn = .FALSE.
  INTEGER, PARAMETER :: PARAM_isopo2 = 70
  INTEGER            ::     P_isopo2 = 1
  LOGICAL            ::     F_isopo2 = .FALSE.
  INTEGER, PARAMETER :: PARAM_dms = 71
  INTEGER            ::     P_dms = 1
  LOGICAL            ::     F_dms = .FALSE.
  INTEGER, PARAMETER :: PARAM_msa = 72
  INTEGER            ::     P_msa = 1
  LOGICAL            ::     F_msa = .FALSE.
  INTEGER, PARAMETER :: PARAM_dmso = 73
  INTEGER            ::     P_dmso = 1
  LOGICAL            ::     F_dmso = .FALSE.
  INTEGER, PARAMETER :: PARAM_dmso2 = 74
  INTEGER            ::     P_dmso2 = 1
  LOGICAL            ::     F_dmso2 = .FALSE.
  INTEGER, PARAMETER :: PARAM_ch3so2h = 75
  INTEGER            ::     P_ch3so2h = 1
  LOGICAL            ::     F_ch3so2h = .FALSE.
  INTEGER, PARAMETER :: PARAM_ch3sch2oo = 76
  INTEGER            ::     P_ch3sch2oo = 1
  LOGICAL            ::     F_ch3sch2oo = .FALSE.
  INTEGER, PARAMETER :: PARAM_ch3so2 = 77
  INTEGER            ::     P_ch3so2 = 1
  LOGICAL            ::     F_ch3so2 = .FALSE.
  INTEGER, PARAMETER :: PARAM_ch3so3 = 78
  INTEGER            ::     P_ch3so3 = 1
  LOGICAL            ::     F_ch3so3 = .FALSE.
  INTEGER, PARAMETER :: PARAM_ch3so2oo = 79
  INTEGER            ::     P_ch3so2oo = 1
  LOGICAL            ::     F_ch3so2oo = .FALSE.
  INTEGER, PARAMETER :: PARAM_ch3so2ch2oo = 80
  INTEGER            ::     P_ch3so2ch2oo = 1
  LOGICAL            ::     F_ch3so2ch2oo = .FALSE.
  INTEGER, PARAMETER :: PARAM_mtf = 81
  INTEGER            ::     P_mtf = 1
  LOGICAL            ::     F_mtf = .FALSE.
  INTEGER, PARAMETER :: PARAM_ald2 = 82
  INTEGER            ::     P_ald2 = 1
  LOGICAL            ::     F_ald2 = .FALSE.
  INTEGER, PARAMETER :: PARAM_ror = 83
  INTEGER            ::     P_ror = 1
  LOGICAL            ::     F_ror = .FALSE.
  INTEGER, PARAMETER :: PARAM_ole = 84
  INTEGER            ::     P_ole = 1
  LOGICAL            ::     F_ole = .FALSE.
  INTEGER, PARAMETER :: PARAM_cres = 85
  INTEGER            ::     P_cres = 1
  LOGICAL            ::     F_cres = .FALSE.
  INTEGER, PARAMETER :: PARAM_xo2n = 86
  INTEGER            ::     P_xo2n = 1
  LOGICAL            ::     F_xo2n = .FALSE.
  INTEGER, PARAMETER :: PARAM_pna = 87
  INTEGER            ::     P_pna = 1
  LOGICAL            ::     F_pna = .FALSE.
  INTEGER, PARAMETER :: PARAM_o = 88
  INTEGER            ::     P_o = 1
  LOGICAL            ::     F_o = .FALSE.
  INTEGER, PARAMETER :: PARAM_o1d_cb4 = 89
  INTEGER            ::     P_o1d_cb4 = 1
  LOGICAL            ::     F_o1d_cb4 = .FALSE.
  INTEGER, PARAMETER :: PARAM_tracer_1 = 90
  INTEGER            ::     P_tracer_1 = 1
  LOGICAL            ::     F_tracer_1 = .FALSE.
  INTEGER, PARAMETER :: PARAM_tracer_2 = 91
  INTEGER            ::     P_tracer_2 = 1
  LOGICAL            ::     F_tracer_2 = .FALSE.
  INTEGER, PARAMETER :: PARAM_tracer_3 = 92
  INTEGER            ::     P_tracer_3 = 1
  LOGICAL            ::     F_tracer_3 = .FALSE.
  INTEGER, PARAMETER :: PARAM_tracer_4 = 93
  INTEGER            ::     P_tracer_4 = 1
  LOGICAL            ::     F_tracer_4 = .FALSE.
  INTEGER, PARAMETER :: PARAM_tracer_5 = 94
  INTEGER            ::     P_tracer_5 = 1
  LOGICAL            ::     F_tracer_5 = .FALSE.
  INTEGER, PARAMETER :: PARAM_tracer_6 = 95
  INTEGER            ::     P_tracer_6 = 1
  LOGICAL            ::     F_tracer_6 = .FALSE.
  INTEGER, PARAMETER :: PARAM_tracer_7 = 96
  INTEGER            ::     P_tracer_7 = 1
  LOGICAL            ::     F_tracer_7 = .FALSE.
  INTEGER, PARAMETER :: PARAM_tracer_8 = 97
  INTEGER            ::     P_tracer_8 = 1
  LOGICAL            ::     F_tracer_8 = .FALSE.
  INTEGER, PARAMETER :: PARAM_tracer_9 = 98
  INTEGER            ::     P_tracer_9 = 1
  LOGICAL            ::     F_tracer_9 = .FALSE.
  INTEGER, PARAMETER :: PARAM_tracer_10 = 99
  INTEGER            ::     P_tracer_10 = 1
  LOGICAL            ::     F_tracer_10 = .FALSE.
  INTEGER, PARAMETER :: PARAM_tracer_11 = 100
  INTEGER            ::     P_tracer_11 = 1
  LOGICAL            ::     F_tracer_11 = .FALSE.
  INTEGER, PARAMETER :: PARAM_tracer_12 = 101
  INTEGER            ::     P_tracer_12 = 1
  LOGICAL            ::     F_tracer_12 = .FALSE.
  INTEGER, PARAMETER :: PARAM_tracer_13 = 102
  INTEGER            ::     P_tracer_13 = 1
  LOGICAL            ::     F_tracer_13 = .FALSE.
  INTEGER, PARAMETER :: PARAM_tracer_14 = 103
  INTEGER            ::     P_tracer_14 = 1
  LOGICAL            ::     F_tracer_14 = .FALSE.
  INTEGER, PARAMETER :: PARAM_tracer_15 = 104
  INTEGER            ::     P_tracer_15 = 1
  LOGICAL            ::     F_tracer_15 = .FALSE.
  INTEGER, PARAMETER :: PARAM_tracer_16 = 105
  INTEGER            ::     P_tracer_16 = 1
  LOGICAL            ::     F_tracer_16 = .FALSE.
  INTEGER, PARAMETER :: PARAM_tracer_17 = 106
  INTEGER            ::     P_tracer_17 = 1
  LOGICAL            ::     F_tracer_17 = .FALSE.
  INTEGER, PARAMETER :: PARAM_tracer_18 = 107
  INTEGER            ::     P_tracer_18 = 1
  LOGICAL            ::     F_tracer_18 = .FALSE.
  INTEGER, PARAMETER :: PARAM_tracer_19 = 108
  INTEGER            ::     P_tracer_19 = 1
  LOGICAL            ::     F_tracer_19 = .FALSE.
  INTEGER, PARAMETER :: PARAM_tracer_20 = 109
  INTEGER            ::     P_tracer_20 = 1
  LOGICAL            ::     F_tracer_20 = .FALSE.
  INTEGER, PARAMETER :: PARAM_tracer_ens = 110
  INTEGER            ::     P_tracer_ens = 1
  LOGICAL            ::     F_tracer_ens = .FALSE.
  INTEGER, PARAMETER :: PARAM_vash_1 = 111
  INTEGER            ::     P_vash_1 = 1 
  LOGICAL            ::     F_vash_1 = .FALSE.
  INTEGER, PARAMETER :: PARAM_vash_2 = 112
  INTEGER            ::     P_vash_2 = 1
  LOGICAL            ::     F_vash_2 = .FALSE.
  INTEGER, PARAMETER :: PARAM_vash_3 = 113
  INTEGER            ::     P_vash_3 = 1
  LOGICAL            ::     F_vash_3 = .FALSE.
  INTEGER, PARAMETER :: PARAM_vash_4 = 114
  INTEGER            ::     P_vash_4 = 1
  LOGICAL            ::     F_vash_4 = .FALSE.
  INTEGER, PARAMETER :: PARAM_vash_5 = 115
  INTEGER            ::     P_vash_5 = 1
  LOGICAL            ::     F_vash_5 = .FALSE.
  INTEGER, PARAMETER :: PARAM_vash_6 = 116
  INTEGER            ::     P_vash_6 = 1
  LOGICAL            ::     F_vash_6 = .FALSE.
  INTEGER, PARAMETER :: PARAM_vash_7 = 117
  INTEGER            ::     P_vash_7 = 1
  LOGICAL            ::     F_vash_7 = .FALSE.
  INTEGER, PARAMETER :: PARAM_vash_8 = 118
  INTEGER            ::     P_vash_8 = 1
  LOGICAL            ::     F_vash_8 = .FALSE.
  INTEGER, PARAMETER :: PARAM_vash_9 = 119
  INTEGER            ::     P_vash_9 = 1
  LOGICAL            ::     F_vash_9 = .FALSE.
  INTEGER, PARAMETER :: PARAM_vash_10 = 120
  INTEGER            ::     P_vash_10 = 1
  LOGICAL            ::     F_vash_10 = .FALSE.
  INTEGER, PARAMETER :: PARAM_pm_25 = 121
  INTEGER            ::     P_pm_25 = 1
  LOGICAL            ::     F_pm_25 = .FALSE.
  INTEGER, PARAMETER :: PARAM_pm_10 = 122
  INTEGER            ::     P_pm_10 = 1
  LOGICAL            ::     F_pm_10 = .FALSE.
  INTEGER, PARAMETER :: PARAM_so4aj = 123
  INTEGER            ::     P_so4aj = 1
  LOGICAL            ::     F_so4aj = .FALSE.
  INTEGER, PARAMETER :: PARAM_so4ai = 124
  INTEGER            ::     P_so4ai = 1
  LOGICAL            ::     F_so4ai = .FALSE.
  INTEGER, PARAMETER :: PARAM_nh4aj = 125
  INTEGER            ::     P_nh4aj = 1
  LOGICAL            ::     F_nh4aj = .FALSE.
  INTEGER, PARAMETER :: PARAM_nh4ai = 126
  INTEGER            ::     P_nh4ai = 1
  LOGICAL            ::     F_nh4ai = .FALSE.
  INTEGER, PARAMETER :: PARAM_no3aj = 127
  INTEGER            ::     P_no3aj = 1
  LOGICAL            ::     F_no3aj = .FALSE.
  INTEGER, PARAMETER :: PARAM_no3ai = 128
  INTEGER            ::     P_no3ai = 1
  LOGICAL            ::     F_no3ai = .FALSE.
  INTEGER, PARAMETER :: PARAM_naaj = 129
  INTEGER            ::     P_naaj = 1
  LOGICAL            ::     F_naaj = .FALSE.
  INTEGER, PARAMETER :: PARAM_naai = 130
  INTEGER            ::     P_naai = 1 
  LOGICAL            ::     F_naai = .FALSE.
  INTEGER, PARAMETER :: PARAM_claj = 131
  INTEGER            ::     P_claj = 1 
  LOGICAL            ::     F_claj = .FALSE.
  INTEGER, PARAMETER :: PARAM_clai = 132
  INTEGER            ::     P_clai = 1 
  LOGICAL            ::     F_clai = .FALSE.
  INTEGER, PARAMETER :: PARAM_orgaro1j = 133
  INTEGER            ::     P_orgaro1j = 1
  LOGICAL            ::     F_orgaro1j = .FALSE.
  INTEGER, PARAMETER :: PARAM_orgaro1i = 134
  INTEGER            ::     P_orgaro1i = 1
  LOGICAL            ::     F_orgaro1i = .FALSE.
  INTEGER, PARAMETER :: PARAM_orgaro2j = 135
  INTEGER            ::     P_orgaro2j = 1
  LOGICAL            ::     F_orgaro2j = .FALSE.
  INTEGER, PARAMETER :: PARAM_orgaro2i = 136
  INTEGER            ::     P_orgaro2i = 1
  LOGICAL            ::     F_orgaro2i = .FALSE.
  INTEGER, PARAMETER :: PARAM_orgalk1j = 137
  INTEGER            ::     P_orgalk1j = 1
  LOGICAL            ::     F_orgalk1j = .FALSE.
  INTEGER, PARAMETER :: PARAM_orgalk1i = 138
  INTEGER            ::     P_orgalk1i = 1
  LOGICAL            ::     F_orgalk1i = .FALSE.
  INTEGER, PARAMETER :: PARAM_orgole1j = 139
  INTEGER            ::     P_orgole1j = 1
  LOGICAL            ::     F_orgole1j = .FALSE.
  INTEGER, PARAMETER :: PARAM_orgole1i = 140
  INTEGER            ::     P_orgole1i = 1
  LOGICAL            ::     F_orgole1i = .FALSE.
  INTEGER, PARAMETER :: PARAM_orgba1j = 141
  INTEGER            ::     P_orgba1j = 1
  LOGICAL            ::     F_orgba1j = .FALSE.
  INTEGER, PARAMETER :: PARAM_orgba1i = 142
  INTEGER            ::     P_orgba1i = 1
  LOGICAL            ::     F_orgba1i = .FALSE.
  INTEGER, PARAMETER :: PARAM_orgba2j = 143
  INTEGER            ::     P_orgba2j = 1
  LOGICAL            ::     F_orgba2j = .FALSE.
  INTEGER, PARAMETER :: PARAM_orgba2i = 144
  INTEGER            ::     P_orgba2i = 1
  LOGICAL            ::     F_orgba2i = .FALSE.
  INTEGER, PARAMETER :: PARAM_orgba3j = 145
  INTEGER            ::     P_orgba3j = 1
  LOGICAL            ::     F_orgba3j = .FALSE.
  INTEGER, PARAMETER :: PARAM_orgba3i = 146
  INTEGER            ::     P_orgba3i = 1
  LOGICAL            ::     F_orgba3i = .FALSE.
  INTEGER, PARAMETER :: PARAM_orgba4j = 147
  INTEGER            ::     P_orgba4j = 1
  LOGICAL            ::     F_orgba4j = .FALSE.
  INTEGER, PARAMETER :: PARAM_orgba4i = 148
  INTEGER            ::     P_orgba4i = 1
  LOGICAL            ::     F_orgba4i = .FALSE.
  INTEGER, PARAMETER :: PARAM_orgpaj = 149
  INTEGER            ::     P_orgpaj = 1
  LOGICAL            ::     F_orgpaj = .FALSE.
  INTEGER, PARAMETER :: PARAM_orgpai = 150
  INTEGER            ::     P_orgpai = 1
  LOGICAL            ::     F_orgpai = .FALSE.
  INTEGER, PARAMETER :: PARAM_ecj = 151
  INTEGER            ::     P_ecj = 1
  LOGICAL            ::     F_ecj = .FALSE.
  INTEGER, PARAMETER :: PARAM_eci = 152
  INTEGER            ::     P_eci = 1
  LOGICAL            ::     F_eci = .FALSE.
  INTEGER, PARAMETER :: PARAM_p25j = 153
  INTEGER            ::     P_p25j = 1
  LOGICAL            ::     F_p25j = .FALSE.
  INTEGER, PARAMETER :: PARAM_p25i = 154
  INTEGER            ::     P_p25i = 1
  LOGICAL            ::     F_p25i = .FALSE.
  INTEGER, PARAMETER :: PARAM_antha = 155
  INTEGER            ::     P_antha = 1
  LOGICAL            ::     F_antha = .FALSE.
  INTEGER, PARAMETER :: PARAM_seas = 156
  INTEGER            ::     P_seas = 1
  LOGICAL            ::     F_seas = .FALSE.
  INTEGER, PARAMETER :: PARAM_soila = 157
  INTEGER            ::     P_soila = 1 
  LOGICAL            ::     F_soila = .FALSE.
  INTEGER, PARAMETER :: PARAM_nu0 = 158 
  INTEGER            ::     P_nu0 = 1 
  LOGICAL            ::     F_nu0 = .FALSE.
  INTEGER, PARAMETER :: PARAM_ac0 = 159 
  INTEGER            ::     P_ac0 = 1 
  LOGICAL            ::     F_ac0 = .FALSE.
  INTEGER, PARAMETER :: PARAM_corn = 160
  INTEGER            ::     P_corn = 1
  LOGICAL            ::     F_corn = .FALSE.
  INTEGER, PARAMETER :: PARAM_so4cwj = 161
  INTEGER            ::     P_so4cwj = 1
  LOGICAL            ::     F_so4cwj = .FALSE.
  INTEGER, PARAMETER :: PARAM_so4cwi = 162
  INTEGER            ::     P_so4cwi = 1
  LOGICAL            ::     F_so4cwi = .FALSE.
  INTEGER, PARAMETER :: PARAM_nh4cwj = 163
  INTEGER            ::     P_nh4cwj = 1
  LOGICAL            ::     F_nh4cwj = .FALSE.
  INTEGER, PARAMETER :: PARAM_nh4cwi = 164
  INTEGER            ::     P_nh4cwi = 1
  LOGICAL            ::     F_nh4cwi = .FALSE.
  INTEGER, PARAMETER :: PARAM_no3cwj = 165
  INTEGER            ::     P_no3cwj = 1
  LOGICAL            ::     F_no3cwj = .FALSE.
  INTEGER, PARAMETER :: PARAM_no3cwi = 166
  INTEGER            ::     P_no3cwi = 1
  LOGICAL            ::     F_no3cwi = .FALSE.
  INTEGER, PARAMETER :: PARAM_nacwj = 167
  INTEGER            ::     P_nacwj = 1
  LOGICAL            ::     F_nacwj = .FALSE.
  INTEGER, PARAMETER :: PARAM_nacwi = 168
  INTEGER            ::     P_nacwi = 1
  LOGICAL            ::     F_nacwi = .FALSE.
  INTEGER, PARAMETER :: PARAM_clcwj = 169
  INTEGER            ::     P_clcwj = 1
  LOGICAL            ::     F_clcwj = .FALSE.
  INTEGER, PARAMETER :: PARAM_clcwi = 170
  INTEGER            ::     P_clcwi = 1
  LOGICAL            ::     F_clcwi = .FALSE.
  INTEGER, PARAMETER :: PARAM_orgaro1cwj = 171
  INTEGER            ::     P_orgaro1cwj = 1
  LOGICAL            ::     F_orgaro1cwj = .FALSE.
  INTEGER, PARAMETER :: PARAM_orgaro1cwi = 172
  INTEGER            ::     P_orgaro1cwi = 1
  LOGICAL            ::     F_orgaro1cwi = .FALSE.
  INTEGER, PARAMETER :: PARAM_orgaro2cwj = 173
  INTEGER            ::     P_orgaro2cwj = 1
  LOGICAL            ::     F_orgaro2cwj = .FALSE.
  INTEGER, PARAMETER :: PARAM_orgaro2cwi = 174
  INTEGER            ::     P_orgaro2cwi = 1
  LOGICAL            ::     F_orgaro2cwi = .FALSE.
  INTEGER, PARAMETER :: PARAM_orgalk1cwj = 175
  INTEGER            ::     P_orgalk1cwj = 1
  LOGICAL            ::     F_orgalk1cwj = .FALSE.
  INTEGER, PARAMETER :: PARAM_orgalk1cwi = 176
  INTEGER            ::     P_orgalk1cwi = 1
  LOGICAL            ::     F_orgalk1cwi = .FALSE.
  INTEGER, PARAMETER :: PARAM_orgole1cwj = 177
  INTEGER            ::     P_orgole1cwj = 1
  LOGICAL            ::     F_orgole1cwj = .FALSE.
  INTEGER, PARAMETER :: PARAM_orgole1cwi = 178
  INTEGER            ::     P_orgole1cwi = 1
  LOGICAL            ::     F_orgole1cwi = .FALSE.
  INTEGER, PARAMETER :: PARAM_orgba1cwj = 179
  INTEGER            ::     P_orgba1cwj = 1
  LOGICAL            ::     F_orgba1cwj = .FALSE.
  INTEGER, PARAMETER :: PARAM_orgba1cwi = 180
  INTEGER            ::     P_orgba1cwi = 1
  LOGICAL            ::     F_orgba1cwi = .FALSE.
  INTEGER, PARAMETER :: PARAM_orgba2cwj = 181
  INTEGER            ::     P_orgba2cwj = 1
  LOGICAL            ::     F_orgba2cwj = .FALSE.
  INTEGER, PARAMETER :: PARAM_orgba2cwi = 182
  INTEGER            ::     P_orgba2cwi = 1
  LOGICAL            ::     F_orgba2cwi = .FALSE.
  INTEGER, PARAMETER :: PARAM_orgba3cwj = 183
  INTEGER            ::     P_orgba3cwj = 1
  LOGICAL            ::     F_orgba3cwj = .FALSE.
  INTEGER, PARAMETER :: PARAM_orgba3cwi = 184
  INTEGER            ::     P_orgba3cwi = 1
  LOGICAL            ::     F_orgba3cwi = .FALSE.
  INTEGER, PARAMETER :: PARAM_orgba4cwj = 185
  INTEGER            ::     P_orgba4cwj = 1
  LOGICAL            ::     F_orgba4cwj = .FALSE.
  INTEGER, PARAMETER :: PARAM_orgba4cwi = 186
  INTEGER            ::     P_orgba4cwi = 1
  LOGICAL            ::     F_orgba4cwi = .FALSE.
  INTEGER, PARAMETER :: PARAM_orgpacwj = 187
  INTEGER            ::     P_orgpacwj = 1 
  LOGICAL            ::     F_orgpacwj = .FALSE.
  INTEGER, PARAMETER :: PARAM_orgpacwi = 188
  INTEGER            ::     P_orgpacwi = 1 
  LOGICAL            ::     F_orgpacwi = .FALSE.
  INTEGER, PARAMETER :: PARAM_eccwj = 189
  INTEGER            ::     P_eccwj = 1
  LOGICAL            ::     F_eccwj = .FALSE.
  INTEGER, PARAMETER :: PARAM_eccwi = 190
  INTEGER            ::     P_eccwi = 1
  LOGICAL            ::     F_eccwi = .FALSE.
  INTEGER, PARAMETER :: PARAM_p25cwj = 191 
  INTEGER            ::     P_p25cwj = 1 
  LOGICAL            ::     F_p25cwj = .FALSE.
  INTEGER, PARAMETER :: PARAM_p25cwi = 192
  INTEGER            ::     P_p25cwi = 1
  LOGICAL            ::     F_p25cwi = .FALSE.
  INTEGER, PARAMETER :: PARAM_anthcw = 193
  INTEGER            ::     P_anthcw = 1
  LOGICAL            ::     F_anthcw = .FALSE.
  INTEGER, PARAMETER :: PARAM_seascw = 194
  INTEGER            ::     P_seascw = 1
  LOGICAL            ::     F_seascw = .FALSE.
  INTEGER, PARAMETER :: PARAM_soilcw = 195
  INTEGER            ::     P_soilcw = 1
  LOGICAL            ::     F_soilcw = .FALSE.
  INTEGER, PARAMETER :: PARAM_nu0cw = 196
  INTEGER            ::     P_nu0cw = 1
  LOGICAL            ::     F_nu0cw = .FALSE.
  INTEGER, PARAMETER :: PARAM_ac0cw = 197
  INTEGER            ::     P_ac0cw = 1
  LOGICAL            ::     F_ac0cw = .FALSE.
  INTEGER, PARAMETER :: PARAM_corncw = 198
  INTEGER            ::     P_corncw = 1
  LOGICAL            ::     F_corncw = .FALSE.
  INTEGER, PARAMETER :: PARAM_hace = 199
  INTEGER            ::     P_hace = 1
  LOGICAL            ::     F_hace = .FALSE.
  INTEGER, PARAMETER :: PARAM_ishp = 200
  INTEGER            ::     P_ishp = 1
  LOGICAL            ::     F_ishp = .FALSE.
  INTEGER, PARAMETER :: PARAM_ison = 201
  INTEGER            ::     P_ison = 1
  LOGICAL            ::     F_ison = .FALSE.
  INTEGER, PARAMETER :: PARAM_mahp = 202
  INTEGER            ::     P_mahp = 1
  LOGICAL            ::     F_mahp = .FALSE.
  INTEGER, PARAMETER :: PARAM_mpan = 203
  INTEGER            ::     P_mpan = 1
  LOGICAL            ::     F_mpan = .FALSE.
  INTEGER, PARAMETER :: PARAM_nald = 204
  INTEGER            ::     P_nald = 1
  LOGICAL            ::     F_nald = .FALSE.
  INTEGER, PARAMETER :: PARAM_so4_a01 = 205
  INTEGER            ::     P_so4_a01 = 1
  LOGICAL            ::     F_so4_a01 = .FALSE.
  INTEGER, PARAMETER :: PARAM_no3_a01 = 206
  INTEGER            ::     P_no3_a01 = 1
  LOGICAL            ::     F_no3_a01 = .FALSE.
  INTEGER, PARAMETER :: PARAM_cl_a01 = 207
  INTEGER            ::     P_cl_a01 = 1
  LOGICAL            ::     F_cl_a01 = .FALSE.
  INTEGER, PARAMETER :: PARAM_msa_a01 = 208
  INTEGER            ::     P_msa_a01 = 1
  LOGICAL            ::     F_msa_a01 = .FALSE.
  INTEGER, PARAMETER :: PARAM_co3_a01 = 209
  INTEGER            ::     P_co3_a01 = 1
  LOGICAL            ::     F_co3_a01 = .FALSE.
  INTEGER, PARAMETER :: PARAM_nh4_a01 = 210
  INTEGER            ::     P_nh4_a01 = 1
  LOGICAL            ::     F_nh4_a01 = .FALSE.
  INTEGER, PARAMETER :: PARAM_na_a01 = 211
  INTEGER            ::     P_na_a01 = 1
  LOGICAL            ::     F_na_a01 = .FALSE.
  INTEGER, PARAMETER :: PARAM_ca_a01 = 212
  INTEGER            ::     P_ca_a01 = 1
  LOGICAL            ::     F_ca_a01 = .FALSE.
  INTEGER, PARAMETER :: PARAM_oin_a01 = 213
  INTEGER            ::     P_oin_a01 = 1
  LOGICAL            ::     F_oin_a01 = .FALSE.
  INTEGER, PARAMETER :: PARAM_oc_a01 = 214
  INTEGER            ::     P_oc_a01 = 1
  LOGICAL            ::     F_oc_a01 = .FALSE.
  INTEGER, PARAMETER :: PARAM_bc_a01 = 215
  INTEGER            ::     P_bc_a01 = 1
  LOGICAL            ::     F_bc_a01 = .FALSE.
  INTEGER, PARAMETER :: PARAM_hysw_a01 = 216
  INTEGER            ::     P_hysw_a01 = 1
  LOGICAL            ::     F_hysw_a01 = .FALSE.
  INTEGER, PARAMETER :: PARAM_water_a01 = 217
  INTEGER            ::     P_water_a01 = 1
  LOGICAL            ::     F_water_a01 = .FALSE.
  INTEGER, PARAMETER :: PARAM_num_a01 = 218
  INTEGER            ::     P_num_a01 = 1
  LOGICAL            ::     F_num_a01 = .FALSE.
  INTEGER, PARAMETER :: PARAM_so4_a02 = 219
  INTEGER            ::     P_so4_a02 = 1
  LOGICAL            ::     F_so4_a02 = .FALSE.
  INTEGER, PARAMETER :: PARAM_no3_a02 = 220
  INTEGER            ::     P_no3_a02 = 1
  LOGICAL            ::     F_no3_a02 = .FALSE.
  INTEGER, PARAMETER :: PARAM_cl_a02 = 221
  INTEGER            ::     P_cl_a02 = 1
  LOGICAL            ::     F_cl_a02 = .FALSE.
  INTEGER, PARAMETER :: PARAM_msa_a02 = 222
  INTEGER            ::     P_msa_a02 = 1
  LOGICAL            ::     F_msa_a02 = .FALSE.
  INTEGER, PARAMETER :: PARAM_co3_a02 = 223
  INTEGER            ::     P_co3_a02 = 1
  LOGICAL            ::     F_co3_a02 = .FALSE.
  INTEGER, PARAMETER :: PARAM_nh4_a02 = 224
  INTEGER            ::     P_nh4_a02 = 1
  LOGICAL            ::     F_nh4_a02 = .FALSE.
  INTEGER, PARAMETER :: PARAM_na_a02 = 225
  INTEGER            ::     P_na_a02 = 1
  LOGICAL            ::     F_na_a02 = .FALSE.
  INTEGER, PARAMETER :: PARAM_ca_a02 = 226
  INTEGER            ::     P_ca_a02 = 1
  LOGICAL            ::     F_ca_a02 = .FALSE.
  INTEGER, PARAMETER :: PARAM_oin_a02 = 227
  INTEGER            ::     P_oin_a02 = 1
  LOGICAL            ::     F_oin_a02 = .FALSE.
  INTEGER, PARAMETER :: PARAM_oc_a02 = 228
  INTEGER            ::     P_oc_a02 = 1
  LOGICAL            ::     F_oc_a02 = .FALSE.
  INTEGER, PARAMETER :: PARAM_bc_a02 = 229
  INTEGER            ::     P_bc_a02 = 1
  LOGICAL            ::     F_bc_a02 = .FALSE.
  INTEGER, PARAMETER :: PARAM_hysw_a02 = 230
  INTEGER            ::     P_hysw_a02 = 1
  LOGICAL            ::     F_hysw_a02 = .FALSE.
  INTEGER, PARAMETER :: PARAM_water_a02 = 231
  INTEGER            ::     P_water_a02 = 1
  LOGICAL            ::     F_water_a02 = .FALSE.
  INTEGER, PARAMETER :: PARAM_num_a02 = 232
  INTEGER            ::     P_num_a02 = 1
  LOGICAL            ::     F_num_a02 = .FALSE.
  INTEGER, PARAMETER :: PARAM_so4_a03 = 233
  INTEGER            ::     P_so4_a03 = 1
  LOGICAL            ::     F_so4_a03 = .FALSE.
  INTEGER, PARAMETER :: PARAM_no3_a03 = 234
  INTEGER            ::     P_no3_a03 = 1
  LOGICAL            ::     F_no3_a03 = .FALSE.
  INTEGER, PARAMETER :: PARAM_cl_a03 = 235
  INTEGER            ::     P_cl_a03 = 1
  LOGICAL            ::     F_cl_a03 = .FALSE.
  INTEGER, PARAMETER :: PARAM_msa_a03 = 236
  INTEGER            ::     P_msa_a03 = 1
  LOGICAL            ::     F_msa_a03 = .FALSE.
  INTEGER, PARAMETER :: PARAM_co3_a03 = 237
  INTEGER            ::     P_co3_a03 = 1
  LOGICAL            ::     F_co3_a03 = .FALSE.
  INTEGER, PARAMETER :: PARAM_nh4_a03 = 238
  INTEGER            ::     P_nh4_a03 = 1
  LOGICAL            ::     F_nh4_a03 = .FALSE.
  INTEGER, PARAMETER :: PARAM_na_a03 = 239
  INTEGER            ::     P_na_a03 = 1
  LOGICAL            ::     F_na_a03 = .FALSE.
  INTEGER, PARAMETER :: PARAM_ca_a03 = 240
  INTEGER            ::     P_ca_a03 = 1
  LOGICAL            ::     F_ca_a03 = .FALSE.
  INTEGER, PARAMETER :: PARAM_oin_a03 = 241
  INTEGER            ::     P_oin_a03 = 1
  LOGICAL            ::     F_oin_a03 = .FALSE.
  INTEGER, PARAMETER :: PARAM_oc_a03 = 242
  INTEGER            ::     P_oc_a03 = 1
  LOGICAL            ::     F_oc_a03 = .FALSE.
  INTEGER, PARAMETER :: PARAM_bc_a03 = 243
  INTEGER            ::     P_bc_a03 = 1 
  LOGICAL            ::     F_bc_a03 = .FALSE.
  INTEGER, PARAMETER :: PARAM_hysw_a03 = 244
  INTEGER            ::     P_hysw_a03 = 1
  LOGICAL            ::     F_hysw_a03 = .FALSE.
  INTEGER, PARAMETER :: PARAM_water_a03 = 245
  INTEGER            ::     P_water_a03 = 1
  LOGICAL            ::     F_water_a03 = .FALSE.
  INTEGER, PARAMETER :: PARAM_num_a03 = 246
  INTEGER            ::     P_num_a03 = 1
  LOGICAL            ::     F_num_a03 = .FALSE.
  INTEGER, PARAMETER :: PARAM_so4_a04 = 247
  INTEGER            ::     P_so4_a04 = 1
  LOGICAL            ::     F_so4_a04 = .FALSE.
  INTEGER, PARAMETER :: PARAM_no3_a04 = 248
  INTEGER            ::     P_no3_a04 = 1
  LOGICAL            ::     F_no3_a04 = .FALSE.
  INTEGER, PARAMETER :: PARAM_cl_a04 = 249
  INTEGER            ::     P_cl_a04 = 1
  LOGICAL            ::     F_cl_a04 = .FALSE.
  INTEGER, PARAMETER :: PARAM_msa_a04 = 250
  INTEGER            ::     P_msa_a04 = 1
  LOGICAL            ::     F_msa_a04 = .FALSE.
  INTEGER, PARAMETER :: PARAM_co3_a04 = 251
  INTEGER            ::     P_co3_a04 = 1
  LOGICAL            ::     F_co3_a04 = .FALSE.
  INTEGER, PARAMETER :: PARAM_nh4_a04 = 252
  INTEGER            ::     P_nh4_a04 = 1
  LOGICAL            ::     F_nh4_a04 = .FALSE.
  INTEGER, PARAMETER :: PARAM_na_a04 = 253
  INTEGER            ::     P_na_a04 = 1
  LOGICAL            ::     F_na_a04 = .FALSE.
  INTEGER, PARAMETER :: PARAM_ca_a04 = 254
  INTEGER            ::     P_ca_a04 = 1
  LOGICAL            ::     F_ca_a04 = .FALSE.
  INTEGER, PARAMETER :: PARAM_oin_a04 = 255
  INTEGER            ::     P_oin_a04 = 1
  LOGICAL            ::     F_oin_a04 = .FALSE.
  INTEGER, PARAMETER :: PARAM_oc_a04 = 256
  INTEGER            ::     P_oc_a04 = 1
  LOGICAL            ::     F_oc_a04 = .FALSE.
  INTEGER, PARAMETER :: PARAM_bc_a04 = 257
  INTEGER            ::     P_bc_a04 = 1
  LOGICAL            ::     F_bc_a04 = .FALSE.
  INTEGER, PARAMETER :: PARAM_hysw_a04 = 258
  INTEGER            ::     P_hysw_a04 = 1
  LOGICAL            ::     F_hysw_a04 = .FALSE.
  INTEGER, PARAMETER :: PARAM_water_a04 = 259
  INTEGER            ::     P_water_a04 = 1
  LOGICAL            ::     F_water_a04 = .FALSE.
  INTEGER, PARAMETER :: PARAM_num_a04 = 260
  INTEGER            ::     P_num_a04 = 1
  LOGICAL            ::     F_num_a04 = .FALSE.
  INTEGER, PARAMETER :: PARAM_so4_a05 = 261
  INTEGER            ::     P_so4_a05 = 1
  LOGICAL            ::     F_so4_a05 = .FALSE.
  INTEGER, PARAMETER :: PARAM_no3_a05 = 262
  INTEGER            ::     P_no3_a05 = 1
  LOGICAL            ::     F_no3_a05 = .FALSE.
  INTEGER, PARAMETER :: PARAM_cl_a05 = 263
  INTEGER            ::     P_cl_a05 = 1
  LOGICAL            ::     F_cl_a05 = .FALSE.
  INTEGER, PARAMETER :: PARAM_msa_a05 = 264
  INTEGER            ::     P_msa_a05 = 1
  LOGICAL            ::     F_msa_a05 = .FALSE.
  INTEGER, PARAMETER :: PARAM_co3_a05 = 265
  INTEGER            ::     P_co3_a05 = 1
  LOGICAL            ::     F_co3_a05 = .FALSE.
  INTEGER, PARAMETER :: PARAM_nh4_a05 = 266
  INTEGER            ::     P_nh4_a05 = 1
  LOGICAL            ::     F_nh4_a05 = .FALSE.
  INTEGER, PARAMETER :: PARAM_na_a05 = 267
  INTEGER            ::     P_na_a05 = 1
  LOGICAL            ::     F_na_a05 = .FALSE.
  INTEGER, PARAMETER :: PARAM_ca_a05 = 268
  INTEGER            ::     P_ca_a05 = 1
  LOGICAL            ::     F_ca_a05 = .FALSE.
  INTEGER, PARAMETER :: PARAM_oin_a05 = 269
  INTEGER            ::     P_oin_a05 = 1
  LOGICAL            ::     F_oin_a05 = .FALSE.
  INTEGER, PARAMETER :: PARAM_oc_a05 = 270
  INTEGER            ::     P_oc_a05 = 1
  LOGICAL            ::     F_oc_a05 = .FALSE.
  INTEGER, PARAMETER :: PARAM_bc_a05 = 271
  INTEGER            ::     P_bc_a05 = 1 
  LOGICAL            ::     F_bc_a05 = .FALSE.
  INTEGER, PARAMETER :: PARAM_hysw_a05 = 272
  INTEGER            ::     P_hysw_a05 = 1
  LOGICAL            ::     F_hysw_a05 = .FALSE.
  INTEGER, PARAMETER :: PARAM_water_a05 = 273
  INTEGER            ::     P_water_a05 = 1
  LOGICAL            ::     F_water_a05 = .FALSE.
  INTEGER, PARAMETER :: PARAM_num_a05 = 274
  INTEGER            ::     P_num_a05 = 1
  LOGICAL            ::     F_num_a05 = .FALSE.
  INTEGER, PARAMETER :: PARAM_so4_a06 = 275
  INTEGER            ::     P_so4_a06 = 1
  LOGICAL            ::     F_so4_a06 = .FALSE.
  INTEGER, PARAMETER :: PARAM_no3_a06 = 276
  INTEGER            ::     P_no3_a06 = 1
  LOGICAL            ::     F_no3_a06 = .FALSE.
  INTEGER, PARAMETER :: PARAM_cl_a06 = 277
  INTEGER            ::     P_cl_a06 = 1
  LOGICAL            ::     F_cl_a06 = .FALSE.
  INTEGER, PARAMETER :: PARAM_msa_a06 = 278
  INTEGER            ::     P_msa_a06 = 1
  LOGICAL            ::     F_msa_a06 = .FALSE.
  INTEGER, PARAMETER :: PARAM_co3_a06 = 279
  INTEGER            ::     P_co3_a06 = 1
  LOGICAL            ::     F_co3_a06 = .FALSE.
  INTEGER, PARAMETER :: PARAM_nh4_a06 = 280
  INTEGER            ::     P_nh4_a06 = 1
  LOGICAL            ::     F_nh4_a06 = .FALSE.
  INTEGER, PARAMETER :: PARAM_na_a06 = 281
  INTEGER            ::     P_na_a06 = 1
  LOGICAL            ::     F_na_a06 = .FALSE.
  INTEGER, PARAMETER :: PARAM_ca_a06 = 282
  INTEGER            ::     P_ca_a06 = 1
  LOGICAL            ::     F_ca_a06 = .FALSE.
  INTEGER, PARAMETER :: PARAM_oin_a06 = 283
  INTEGER            ::     P_oin_a06 = 1
  LOGICAL            ::     F_oin_a06 = .FALSE.
  INTEGER, PARAMETER :: PARAM_oc_a06 = 284
  INTEGER            ::     P_oc_a06 = 1
  LOGICAL            ::     F_oc_a06 = .FALSE.
  INTEGER, PARAMETER :: PARAM_bc_a06 = 285
  INTEGER            ::     P_bc_a06 = 1
  LOGICAL            ::     F_bc_a06 = .FALSE.
  INTEGER, PARAMETER :: PARAM_hysw_a06 = 286
  INTEGER            ::     P_hysw_a06 = 1
  LOGICAL            ::     F_hysw_a06 = .FALSE.
  INTEGER, PARAMETER :: PARAM_water_a06 = 287
  INTEGER            ::     P_water_a06 = 1
  LOGICAL            ::     F_water_a06 = .FALSE.
  INTEGER, PARAMETER :: PARAM_num_a06 = 288
  INTEGER            ::     P_num_a06 = 1
  LOGICAL            ::     F_num_a06 = .FALSE.
  INTEGER, PARAMETER :: PARAM_so4_a07 = 289
  INTEGER            ::     P_so4_a07 = 1
  LOGICAL            ::     F_so4_a07 = .FALSE.
  INTEGER, PARAMETER :: PARAM_no3_a07 = 290
  INTEGER            ::     P_no3_a07 = 1
  LOGICAL            ::     F_no3_a07 = .FALSE.
  INTEGER, PARAMETER :: PARAM_cl_a07 = 291
  INTEGER            ::     P_cl_a07 = 1
  LOGICAL            ::     F_cl_a07 = .FALSE.
  INTEGER, PARAMETER :: PARAM_msa_a07 = 292
  INTEGER            ::     P_msa_a07 = 1
  LOGICAL            ::     F_msa_a07 = .FALSE.
  INTEGER, PARAMETER :: PARAM_co3_a07 = 293
  INTEGER            ::     P_co3_a07 = 1
  LOGICAL            ::     F_co3_a07 = .FALSE.
  INTEGER, PARAMETER :: PARAM_nh4_a07 = 294
  INTEGER            ::     P_nh4_a07 = 1
  LOGICAL            ::     F_nh4_a07 = .FALSE.
  INTEGER, PARAMETER :: PARAM_na_a07 = 295
  INTEGER            ::     P_na_a07 = 1
  LOGICAL            ::     F_na_a07 = .FALSE.
  INTEGER, PARAMETER :: PARAM_ca_a07 = 296
  INTEGER            ::     P_ca_a07 = 1
  LOGICAL            ::     F_ca_a07 = .FALSE.
  INTEGER, PARAMETER :: PARAM_oin_a07 = 297
  INTEGER            ::     P_oin_a07 = 1
  LOGICAL            ::     F_oin_a07 = .FALSE.
  INTEGER, PARAMETER :: PARAM_oc_a07 = 298
  INTEGER            ::     P_oc_a07 = 1
  LOGICAL            ::     F_oc_a07 = .FALSE.
  INTEGER, PARAMETER :: PARAM_bc_a07 = 299
  INTEGER            ::     P_bc_a07 = 1
  LOGICAL            ::     F_bc_a07 = .FALSE.
  INTEGER, PARAMETER :: PARAM_hysw_a07 = 300
  INTEGER            ::     P_hysw_a07 = 1
  LOGICAL            ::     F_hysw_a07 = .FALSE.
  INTEGER, PARAMETER :: PARAM_water_a07 = 301
  INTEGER            ::     P_water_a07 = 1
  LOGICAL            ::     F_water_a07 = .FALSE.
  INTEGER, PARAMETER :: PARAM_num_a07 = 302
  INTEGER            ::     P_num_a07 = 1
  LOGICAL            ::     F_num_a07 = .FALSE.
  INTEGER, PARAMETER :: PARAM_so4_a08 = 303
  INTEGER            ::     P_so4_a08 = 1
  LOGICAL            ::     F_so4_a08 = .FALSE.
  INTEGER, PARAMETER :: PARAM_no3_a08 = 304
  INTEGER            ::     P_no3_a08 = 1
  LOGICAL            ::     F_no3_a08 = .FALSE.
  INTEGER, PARAMETER :: PARAM_cl_a08 = 305
  INTEGER            ::     P_cl_a08 = 1
  LOGICAL            ::     F_cl_a08 = .FALSE.
  INTEGER, PARAMETER :: PARAM_msa_a08 = 306
  INTEGER            ::     P_msa_a08 = 1
  LOGICAL            ::     F_msa_a08 = .FALSE.
  INTEGER, PARAMETER :: PARAM_co3_a08 = 307
  INTEGER            ::     P_co3_a08 = 1
  LOGICAL            ::     F_co3_a08 = .FALSE.
  INTEGER, PARAMETER :: PARAM_nh4_a08 = 308
  INTEGER            ::     P_nh4_a08 = 1
  LOGICAL            ::     F_nh4_a08 = .FALSE.
  INTEGER, PARAMETER :: PARAM_na_a08 = 309
  INTEGER            ::     P_na_a08 = 1
  LOGICAL            ::     F_na_a08 = .FALSE.
  INTEGER, PARAMETER :: PARAM_ca_a08 = 310
  INTEGER            ::     P_ca_a08 = 1
  LOGICAL            ::     F_ca_a08 = .FALSE.
  INTEGER, PARAMETER :: PARAM_oin_a08 = 311
  INTEGER            ::     P_oin_a08 = 1
  LOGICAL            ::     F_oin_a08 = .FALSE.
  INTEGER, PARAMETER :: PARAM_oc_a08 = 312
  INTEGER            ::     P_oc_a08 = 1
  LOGICAL            ::     F_oc_a08 = .FALSE.
  INTEGER, PARAMETER :: PARAM_bc_a08 = 313
  INTEGER            ::     P_bc_a08 = 1 
  LOGICAL            ::     F_bc_a08 = .FALSE.
  INTEGER, PARAMETER :: PARAM_hysw_a08 = 314
  INTEGER            ::     P_hysw_a08 = 1
  LOGICAL            ::     F_hysw_a08 = .FALSE.
  INTEGER, PARAMETER :: PARAM_water_a08 = 315
  INTEGER            ::     P_water_a08 = 1
  LOGICAL            ::     F_water_a08 = .FALSE.
  INTEGER, PARAMETER :: PARAM_num_a08 = 316
  INTEGER            ::     P_num_a08 = 1
  LOGICAL            ::     F_num_a08 = .FALSE.
  INTEGER, PARAMETER :: PARAM_so4_cw01 = 317
  INTEGER            ::     P_so4_cw01 = 1
  LOGICAL            ::     F_so4_cw01 = .FALSE.
  INTEGER, PARAMETER :: PARAM_no3_cw01 = 318
  INTEGER            ::     P_no3_cw01 = 1
  LOGICAL            ::     F_no3_cw01 = .FALSE.
  INTEGER, PARAMETER :: PARAM_cl_cw01 = 319
  INTEGER            ::     P_cl_cw01 = 1
  LOGICAL            ::     F_cl_cw01 = .FALSE.
  INTEGER, PARAMETER :: PARAM_msa_cw01 = 320
  INTEGER            ::     P_msa_cw01 = 1
  LOGICAL            ::     F_msa_cw01 = .FALSE.
  INTEGER, PARAMETER :: PARAM_co3_cw01 = 321
  INTEGER            ::     P_co3_cw01 = 1
  LOGICAL            ::     F_co3_cw01 = .FALSE.
  INTEGER, PARAMETER :: PARAM_nh4_cw01 = 322
  INTEGER            ::     P_nh4_cw01 = 1
  LOGICAL            ::     F_nh4_cw01 = .FALSE.
  INTEGER, PARAMETER :: PARAM_na_cw01 = 323
  INTEGER            ::     P_na_cw01 = 1
  LOGICAL            ::     F_na_cw01 = .FALSE.
  INTEGER, PARAMETER :: PARAM_ca_cw01 = 324
  INTEGER            ::     P_ca_cw01 = 1
  LOGICAL            ::     F_ca_cw01 = .FALSE.
  INTEGER, PARAMETER :: PARAM_oin_cw01 = 325
  INTEGER            ::     P_oin_cw01 = 1
  LOGICAL            ::     F_oin_cw01 = .FALSE.
  INTEGER, PARAMETER :: PARAM_oc_cw01 = 326
  INTEGER            ::     P_oc_cw01 = 1
  LOGICAL            ::     F_oc_cw01 = .FALSE.
  INTEGER, PARAMETER :: PARAM_bc_cw01 = 327
  INTEGER            ::     P_bc_cw01 = 1
  LOGICAL            ::     F_bc_cw01 = .FALSE.
  INTEGER, PARAMETER :: PARAM_num_cw01 = 328
  INTEGER            ::     P_num_cw01 = 1
  LOGICAL            ::     F_num_cw01 = .FALSE.
  INTEGER, PARAMETER :: PARAM_so4_cw02 = 329
  INTEGER            ::     P_so4_cw02 = 1
  LOGICAL            ::     F_so4_cw02 = .FALSE.
  INTEGER, PARAMETER :: PARAM_no3_cw02 = 330
  INTEGER            ::     P_no3_cw02 = 1
  LOGICAL            ::     F_no3_cw02 = .FALSE.
  INTEGER, PARAMETER :: PARAM_cl_cw02 = 331
  INTEGER            ::     P_cl_cw02 = 1
  LOGICAL            ::     F_cl_cw02 = .FALSE.
  INTEGER, PARAMETER :: PARAM_msa_cw02 = 332
  INTEGER            ::     P_msa_cw02 = 1
  LOGICAL            ::     F_msa_cw02 = .FALSE.
  INTEGER, PARAMETER :: PARAM_co3_cw02 = 333
  INTEGER            ::     P_co3_cw02 = 1
  LOGICAL            ::     F_co3_cw02 = .FALSE.
  INTEGER, PARAMETER :: PARAM_nh4_cw02 = 334
  INTEGER            ::     P_nh4_cw02 = 1
  LOGICAL            ::     F_nh4_cw02 = .FALSE.
  INTEGER, PARAMETER :: PARAM_na_cw02 = 335
  INTEGER            ::     P_na_cw02 = 1
  LOGICAL            ::     F_na_cw02 = .FALSE.
  INTEGER, PARAMETER :: PARAM_ca_cw02 = 336
  INTEGER            ::     P_ca_cw02 = 1
  LOGICAL            ::     F_ca_cw02 = .FALSE.
  INTEGER, PARAMETER :: PARAM_oin_cw02 = 337
  INTEGER            ::     P_oin_cw02 = 1
  LOGICAL            ::     F_oin_cw02 = .FALSE.
  INTEGER, PARAMETER :: PARAM_oc_cw02 = 338
  INTEGER            ::     P_oc_cw02 = 1
  LOGICAL            ::     F_oc_cw02 = .FALSE.
  INTEGER, PARAMETER :: PARAM_bc_cw02 = 339
  INTEGER            ::     P_bc_cw02 = 1
  LOGICAL            ::     F_bc_cw02 = .FALSE.
  INTEGER, PARAMETER :: PARAM_num_cw02 = 340
  INTEGER            ::     P_num_cw02 = 1
  LOGICAL            ::     F_num_cw02 = .FALSE.
  INTEGER, PARAMETER :: PARAM_so4_cw03 = 341
  INTEGER            ::     P_so4_cw03 = 1
  LOGICAL            ::     F_so4_cw03 = .FALSE.
  INTEGER, PARAMETER :: PARAM_no3_cw03 = 342
  INTEGER            ::     P_no3_cw03 = 1
  LOGICAL            ::     F_no3_cw03 = .FALSE.
  INTEGER, PARAMETER :: PARAM_cl_cw03 = 343
  INTEGER            ::     P_cl_cw03 = 1
  LOGICAL            ::     F_cl_cw03 = .FALSE.
  INTEGER, PARAMETER :: PARAM_msa_cw03 = 344
  INTEGER            ::     P_msa_cw03 = 1
  LOGICAL            ::     F_msa_cw03 = .FALSE.
  INTEGER, PARAMETER :: PARAM_co3_cw03 = 345
  INTEGER            ::     P_co3_cw03 = 1
  LOGICAL            ::     F_co3_cw03 = .FALSE.
  INTEGER, PARAMETER :: PARAM_nh4_cw03 = 346
  INTEGER            ::     P_nh4_cw03 = 1
  LOGICAL            ::     F_nh4_cw03 = .FALSE.
  INTEGER, PARAMETER :: PARAM_na_cw03 = 347
  INTEGER            ::     P_na_cw03 = 1
  LOGICAL            ::     F_na_cw03 = .FALSE.
  INTEGER, PARAMETER :: PARAM_ca_cw03 = 348
  INTEGER            ::     P_ca_cw03 = 1
  LOGICAL            ::     F_ca_cw03 = .FALSE.
  INTEGER, PARAMETER :: PARAM_oin_cw03 = 349
  INTEGER            ::     P_oin_cw03 = 1
  LOGICAL            ::     F_oin_cw03 = .FALSE.
  INTEGER, PARAMETER :: PARAM_oc_cw03 = 350
  INTEGER            ::     P_oc_cw03 = 1
  LOGICAL            ::     F_oc_cw03 = .FALSE.
  INTEGER, PARAMETER :: PARAM_bc_cw03 = 351
  INTEGER            ::     P_bc_cw03 = 1
  LOGICAL            ::     F_bc_cw03 = .FALSE.
  INTEGER, PARAMETER :: PARAM_num_cw03 = 352
  INTEGER            ::     P_num_cw03 = 1
  LOGICAL            ::     F_num_cw03 = .FALSE.
  INTEGER, PARAMETER :: PARAM_so4_cw04 = 353
  INTEGER            ::     P_so4_cw04 = 1
  LOGICAL            ::     F_so4_cw04 = .FALSE.
  INTEGER, PARAMETER :: PARAM_no3_cw04 = 354
  INTEGER            ::     P_no3_cw04 = 1
  LOGICAL            ::     F_no3_cw04 = .FALSE.
  INTEGER, PARAMETER :: PARAM_cl_cw04 = 355
  INTEGER            ::     P_cl_cw04 = 1
  LOGICAL            ::     F_cl_cw04 = .FALSE.
  INTEGER, PARAMETER :: PARAM_msa_cw04 = 356
  INTEGER            ::     P_msa_cw04 = 1
  LOGICAL            ::     F_msa_cw04 = .FALSE.
  INTEGER, PARAMETER :: PARAM_co3_cw04 = 357
  INTEGER            ::     P_co3_cw04 = 1
  LOGICAL            ::     F_co3_cw04 = .FALSE.
  INTEGER, PARAMETER :: PARAM_nh4_cw04 = 358
  INTEGER            ::     P_nh4_cw04 = 1
  LOGICAL            ::     F_nh4_cw04 = .FALSE.
  INTEGER, PARAMETER :: PARAM_na_cw04 = 359
  INTEGER            ::     P_na_cw04 = 1
  LOGICAL            ::     F_na_cw04 = .FALSE.
  INTEGER, PARAMETER :: PARAM_ca_cw04 = 360
  INTEGER            ::     P_ca_cw04 = 1
  LOGICAL            ::     F_ca_cw04 = .FALSE.
  INTEGER, PARAMETER :: PARAM_oin_cw04 = 361
  INTEGER            ::     P_oin_cw04 = 1
  LOGICAL            ::     F_oin_cw04 = .FALSE.
  INTEGER, PARAMETER :: PARAM_oc_cw04 = 362
  INTEGER            ::     P_oc_cw04 = 1
  LOGICAL            ::     F_oc_cw04 = .FALSE.
  INTEGER, PARAMETER :: PARAM_bc_cw04 = 363
  INTEGER            ::     P_bc_cw04 = 1
  LOGICAL            ::     F_bc_cw04 = .FALSE.
  INTEGER, PARAMETER :: PARAM_num_cw04 = 364
  INTEGER            ::     P_num_cw04 = 1
  LOGICAL            ::     F_num_cw04 = .FALSE.
  INTEGER, PARAMETER :: PARAM_so4_cw05 = 365
  INTEGER            ::     P_so4_cw05 = 1
  LOGICAL            ::     F_so4_cw05 = .FALSE.
  INTEGER, PARAMETER :: PARAM_no3_cw05 = 366
  INTEGER            ::     P_no3_cw05 = 1
  LOGICAL            ::     F_no3_cw05 = .FALSE.
  INTEGER, PARAMETER :: PARAM_cl_cw05 = 367
  INTEGER            ::     P_cl_cw05 = 1
  LOGICAL            ::     F_cl_cw05 = .FALSE.
  INTEGER, PARAMETER :: PARAM_msa_cw05 = 368
  INTEGER            ::     P_msa_cw05 = 1
  LOGICAL            ::     F_msa_cw05 = .FALSE.
  INTEGER, PARAMETER :: PARAM_co3_cw05 = 369
  INTEGER            ::     P_co3_cw05 = 1
  LOGICAL            ::     F_co3_cw05 = .FALSE.
  INTEGER, PARAMETER :: PARAM_nh4_cw05 = 370
  INTEGER            ::     P_nh4_cw05 = 1
  LOGICAL            ::     F_nh4_cw05 = .FALSE.
  INTEGER, PARAMETER :: PARAM_na_cw05 = 371
  INTEGER            ::     P_na_cw05 = 1
  LOGICAL            ::     F_na_cw05 = .FALSE.
  INTEGER, PARAMETER :: PARAM_ca_cw05 = 372
  INTEGER            ::     P_ca_cw05 = 1
  LOGICAL            ::     F_ca_cw05 = .FALSE.
  INTEGER, PARAMETER :: PARAM_oin_cw05 = 373
  INTEGER            ::     P_oin_cw05 = 1
  LOGICAL            ::     F_oin_cw05 = .FALSE.
  INTEGER, PARAMETER :: PARAM_oc_cw05 = 374
  INTEGER            ::     P_oc_cw05 = 1
  LOGICAL            ::     F_oc_cw05 = .FALSE.
  INTEGER, PARAMETER :: PARAM_bc_cw05 = 375
  INTEGER            ::     P_bc_cw05 = 1
  LOGICAL            ::     F_bc_cw05 = .FALSE.
  INTEGER, PARAMETER :: PARAM_num_cw05 = 376
  INTEGER            ::     P_num_cw05 = 1
  LOGICAL            ::     F_num_cw05 = .FALSE.
  INTEGER, PARAMETER :: PARAM_so4_cw06 = 377
  INTEGER            ::     P_so4_cw06 = 1
  LOGICAL            ::     F_so4_cw06 = .FALSE.
  INTEGER, PARAMETER :: PARAM_no3_cw06 = 378
  INTEGER            ::     P_no3_cw06 = 1
  LOGICAL            ::     F_no3_cw06 = .FALSE.
  INTEGER, PARAMETER :: PARAM_cl_cw06 = 379
  INTEGER            ::     P_cl_cw06 = 1
  LOGICAL            ::     F_cl_cw06 = .FALSE.
  INTEGER, PARAMETER :: PARAM_msa_cw06 = 380
  INTEGER            ::     P_msa_cw06 = 1
  LOGICAL            ::     F_msa_cw06 = .FALSE.
  INTEGER, PARAMETER :: PARAM_co3_cw06 = 381
  INTEGER            ::     P_co3_cw06 = 1
  LOGICAL            ::     F_co3_cw06 = .FALSE.
  INTEGER, PARAMETER :: PARAM_nh4_cw06 = 382
  INTEGER            ::     P_nh4_cw06 = 1
  LOGICAL            ::     F_nh4_cw06 = .FALSE.
  INTEGER, PARAMETER :: PARAM_na_cw06 = 383
  INTEGER            ::     P_na_cw06 = 1
  LOGICAL            ::     F_na_cw06 = .FALSE.
  INTEGER, PARAMETER :: PARAM_ca_cw06 = 384
  INTEGER            ::     P_ca_cw06 = 1
  LOGICAL            ::     F_ca_cw06 = .FALSE.
  INTEGER, PARAMETER :: PARAM_oin_cw06 = 385
  INTEGER            ::     P_oin_cw06 = 1
  LOGICAL            ::     F_oin_cw06 = .FALSE.
  INTEGER, PARAMETER :: PARAM_oc_cw06 = 386
  INTEGER            ::     P_oc_cw06 = 1
  LOGICAL            ::     F_oc_cw06 = .FALSE.
  INTEGER, PARAMETER :: PARAM_bc_cw06 = 387
  INTEGER            ::     P_bc_cw06 = 1
  LOGICAL            ::     F_bc_cw06 = .FALSE.
  INTEGER, PARAMETER :: PARAM_num_cw06 = 388
  INTEGER            ::     P_num_cw06 = 1
  LOGICAL            ::     F_num_cw06 = .FALSE.
  INTEGER, PARAMETER :: PARAM_so4_cw07 = 389
  INTEGER            ::     P_so4_cw07 = 1
  LOGICAL            ::     F_so4_cw07 = .FALSE.
  INTEGER, PARAMETER :: PARAM_no3_cw07 = 390
  INTEGER            ::     P_no3_cw07 = 1
  LOGICAL            ::     F_no3_cw07 = .FALSE.
  INTEGER, PARAMETER :: PARAM_cl_cw07 = 391
  INTEGER            ::     P_cl_cw07 = 1
  LOGICAL            ::     F_cl_cw07 = .FALSE.
  INTEGER, PARAMETER :: PARAM_msa_cw07 = 392
  INTEGER            ::     P_msa_cw07 = 1
  LOGICAL            ::     F_msa_cw07 = .FALSE.
  INTEGER, PARAMETER :: PARAM_co3_cw07 = 393
  INTEGER            ::     P_co3_cw07 = 1
  LOGICAL            ::     F_co3_cw07 = .FALSE.
  INTEGER, PARAMETER :: PARAM_nh4_cw07 = 394
  INTEGER            ::     P_nh4_cw07 = 1
  LOGICAL            ::     F_nh4_cw07 = .FALSE.
  INTEGER, PARAMETER :: PARAM_na_cw07 = 395
  INTEGER            ::     P_na_cw07 = 1
  LOGICAL            ::     F_na_cw07 = .FALSE.
  INTEGER, PARAMETER :: PARAM_ca_cw07 = 396
  INTEGER            ::     P_ca_cw07 = 1
  LOGICAL            ::     F_ca_cw07 = .FALSE.
  INTEGER, PARAMETER :: PARAM_oin_cw07 = 397
  INTEGER            ::     P_oin_cw07 = 1
  LOGICAL            ::     F_oin_cw07 = .FALSE.
  INTEGER, PARAMETER :: PARAM_oc_cw07 = 398
  INTEGER            ::     P_oc_cw07 = 1
  LOGICAL            ::     F_oc_cw07 = .FALSE.
  INTEGER, PARAMETER :: PARAM_bc_cw07 = 399
  INTEGER            ::     P_bc_cw07 = 1
  LOGICAL            ::     F_bc_cw07 = .FALSE.
  INTEGER, PARAMETER :: PARAM_num_cw07 = 400
  INTEGER            ::     P_num_cw07 = 1
  LOGICAL            ::     F_num_cw07 = .FALSE.
  INTEGER, PARAMETER :: PARAM_so4_cw08 = 401
  INTEGER            ::     P_so4_cw08 = 1
  LOGICAL            ::     F_so4_cw08 = .FALSE.
  INTEGER, PARAMETER :: PARAM_no3_cw08 = 402
  INTEGER            ::     P_no3_cw08 = 1
  LOGICAL            ::     F_no3_cw08 = .FALSE.
  INTEGER, PARAMETER :: PARAM_cl_cw08 = 403
  INTEGER            ::     P_cl_cw08 = 1
  LOGICAL            ::     F_cl_cw08 = .FALSE.
  INTEGER, PARAMETER :: PARAM_msa_cw08 = 404
  INTEGER            ::     P_msa_cw08 = 1
  LOGICAL            ::     F_msa_cw08 = .FALSE.
  INTEGER, PARAMETER :: PARAM_co3_cw08 = 405
  INTEGER            ::     P_co3_cw08 = 1
  LOGICAL            ::     F_co3_cw08 = .FALSE.
  INTEGER, PARAMETER :: PARAM_nh4_cw08 = 406
  INTEGER            ::     P_nh4_cw08 = 1
  LOGICAL            ::     F_nh4_cw08 = .FALSE.
  INTEGER, PARAMETER :: PARAM_na_cw08 = 407
  INTEGER            ::     P_na_cw08 = 1
  LOGICAL            ::     F_na_cw08 = .FALSE.
  INTEGER, PARAMETER :: PARAM_ca_cw08 = 408
  INTEGER            ::     P_ca_cw08 = 1
  LOGICAL            ::     F_ca_cw08 = .FALSE.
  INTEGER, PARAMETER :: PARAM_oin_cw08 = 409
  INTEGER            ::     P_oin_cw08 = 1
  LOGICAL            ::     F_oin_cw08 = .FALSE.
  INTEGER, PARAMETER :: PARAM_oc_cw08 = 410
  INTEGER            ::     P_oc_cw08 = 1
  LOGICAL            ::     F_oc_cw08 = .FALSE.
  INTEGER, PARAMETER :: PARAM_bc_cw08 = 411
  INTEGER            ::     P_bc_cw08 = 1
  LOGICAL            ::     F_bc_cw08 = .FALSE.
  INTEGER, PARAMETER :: PARAM_num_cw08 = 412
  INTEGER            ::     P_num_cw08 = 1
  LOGICAL            ::     F_num_cw08 = .FALSE.
  INTEGER, PARAMETER :: PARAM_bc1 = 413
  INTEGER            ::     P_bc1 = 1
  LOGICAL            ::     F_bc1 = .FALSE.
  INTEGER, PARAMETER :: PARAM_bc2 = 414
  INTEGER            ::     P_bc2 = 1
  LOGICAL            ::     F_bc2 = .FALSE.
  INTEGER, PARAMETER :: PARAM_oc1 = 415
  INTEGER            ::     P_oc1 = 1
  LOGICAL            ::     F_oc1 = .FALSE.
  INTEGER, PARAMETER :: PARAM_oc2 = 416 
  INTEGER            ::     P_oc2 = 1 
  LOGICAL            ::     F_oc2 = .FALSE.
  INTEGER, PARAMETER :: PARAM_p25 = 417
  INTEGER            ::     P_p25 = 1
  LOGICAL            ::     F_p25 = .FALSE.
  INTEGER, PARAMETER :: PARAM_p10 = 418
  INTEGER            ::     P_p10 = 1
  LOGICAL            ::     F_p10 = .FALSE.
  INTEGER, PARAMETER :: PARAM_dust_1 = 419
  INTEGER            ::     P_dust_1 = 1 
  LOGICAL            ::     F_dust_1 = .FALSE.
  INTEGER, PARAMETER :: PARAM_dust_2 = 420
  INTEGER            ::     P_dust_2 = 1
  LOGICAL            ::     F_dust_2 = .FALSE.
  INTEGER, PARAMETER :: PARAM_dust_3 = 421
  INTEGER            ::     P_dust_3 = 1
  LOGICAL            ::     F_dust_3 = .FALSE.
  INTEGER, PARAMETER :: PARAM_dust_4 = 422
  INTEGER            ::     P_dust_4 = 1 
  LOGICAL            ::     F_dust_4 = .FALSE.
  INTEGER, PARAMETER :: PARAM_dust_5 = 423
  INTEGER            ::     P_dust_5 = 1
  LOGICAL            ::     F_dust_5 = .FALSE.
  INTEGER, PARAMETER :: PARAM_seas_1 = 424
  INTEGER            ::     P_seas_1 = 1
  LOGICAL            ::     F_seas_1 = .FALSE.
  INTEGER, PARAMETER :: PARAM_seas_2 = 425
  INTEGER            ::     P_seas_2 = 1 
  LOGICAL            ::     F_seas_2 = .FALSE.
  INTEGER, PARAMETER :: PARAM_seas_3 = 426
  INTEGER            ::     P_seas_3 = 1
  LOGICAL            ::     F_seas_3 = .FALSE.
  INTEGER, PARAMETER :: PARAM_seas_4 = 427
  INTEGER            ::     P_seas_4 = 1
  LOGICAL            ::     F_seas_4 = .FALSE.
  INTEGER, PARAMETER :: PARAM_pa = 428
  INTEGER            ::     P_pa = 1
  LOGICAL            ::     F_pa = .FALSE.
  INTEGER, PARAMETER :: PARAM_aca = 429
  INTEGER            ::     P_aca = 1
  LOGICAL            ::     F_aca = .FALSE.
  INTEGER, PARAMETER :: PARAM_acet = 430
  INTEGER            ::     P_acet = 1
  LOGICAL            ::     F_acet = .FALSE.
  INTEGER, PARAMETER :: PARAM_isopr = 431
  INTEGER            ::     P_isopr = 1
  LOGICAL            ::     F_isopr = .FALSE.
  INTEGER, PARAMETER :: PARAM_mvk = 432
  INTEGER            ::     P_mvk = 1
  LOGICAL            ::     F_mvk = .FALSE.
  INTEGER, PARAMETER :: PARAM_iso2 = 433
  INTEGER            ::     P_iso2 = 1
  LOGICAL            ::     F_iso2 = .FALSE.
  INTEGER, PARAMETER :: PARAM_isooh = 434
  INTEGER            ::     P_isooh = 1
  LOGICAL            ::     F_isooh = .FALSE.
  INTEGER, PARAMETER :: PARAM_mvko2 = 435
  INTEGER            ::     P_mvko2 = 1
  LOGICAL            ::     F_mvko2 = .FALSE.
  INTEGER, PARAMETER :: PARAM_mvkooh = 436
  INTEGER            ::     P_mvkooh = 1
  LOGICAL            ::     F_mvkooh = .FALSE.
  INTEGER, PARAMETER :: PARAM_acol = 437
  INTEGER            ::     P_acol = 1
  LOGICAL            ::     F_acol = .FALSE.
  INTEGER, PARAMETER :: PARAM_hcooh = 438
  INTEGER            ::     P_hcooh = 1
  LOGICAL            ::     F_hcooh = .FALSE.
  INTEGER, PARAMETER :: PARAM_naca = 439
  INTEGER            ::     P_naca = 1
  LOGICAL            ::     F_naca = .FALSE.
  INTEGER, PARAMETER :: PARAM_mglo = 440
  INTEGER            ::     P_mglo = 1
  LOGICAL            ::     F_mglo = .FALSE.
  INTEGER, PARAMETER :: PARAM_c2h6 = 441
  INTEGER            ::     P_c2h6 = 1
  LOGICAL            ::     F_c2h6 = .FALSE.
  INTEGER, PARAMETER :: PARAM_etooh = 442
  INTEGER            ::     P_etooh = 1
  LOGICAL            ::     F_etooh = .FALSE.
  INTEGER, PARAMETER :: PARAM_c3h8 = 443
  INTEGER            ::     P_c3h8 = 1
  LOGICAL            ::     F_c3h8 = .FALSE.
  INTEGER, PARAMETER :: PARAM_prooh = 444
  INTEGER            ::     P_prooh = 1
  LOGICAL            ::     F_prooh = .FALSE.
  INTEGER, PARAMETER :: PARAM_acooh = 445
  INTEGER            ::     P_acooh = 1
  LOGICAL            ::     F_acooh = .FALSE.
  INTEGER, PARAMETER :: PARAM_eto2 = 446
  INTEGER            ::     P_eto2 = 1
  LOGICAL            ::     F_eto2 = .FALSE.
  INTEGER, PARAMETER :: PARAM_pro2 = 447
  INTEGER            ::     P_pro2 = 1
  LOGICAL            ::     F_pro2 = .FALSE.
  INTEGER, PARAMETER :: PARAM_aco2 = 448
  INTEGER            ::     P_aco2 = 1
  LOGICAL            ::     F_aco2 = .FALSE.
  INTEGER, PARAMETER :: PARAM_c3h6 = 449
  INTEGER            ::     P_c3h6 = 1 
  LOGICAL            ::     F_c3h6 = .FALSE.
  INTEGER, PARAMETER :: PARAM_c3h6ooh = 450
  INTEGER            ::     P_c3h6ooh = 1
  LOGICAL            ::     F_c3h6ooh = .FALSE.
  INTEGER, PARAMETER :: PARAM_c2h4 = 451
  INTEGER            ::     P_c2h4 = 1
  LOGICAL            ::     F_c2h4 = .FALSE.
  INTEGER, PARAMETER :: PARAM_c4h10 = 452
  INTEGER            ::     P_c4h10 = 1
  LOGICAL            ::     F_c4h10 = .FALSE.
  INTEGER, PARAMETER :: PARAM_buooh = 453
  INTEGER            ::     P_buooh = 1
  LOGICAL            ::     F_buooh = .FALSE.
  INTEGER, PARAMETER :: PARAM_mek = 454
  INTEGER            ::     P_mek = 1
  LOGICAL            ::     F_mek = .FALSE.
  INTEGER, PARAMETER :: PARAM_mekooh = 455
  INTEGER            ::     P_mekooh = 1
  LOGICAL            ::     F_mekooh = .FALSE.
  INTEGER, PARAMETER :: PARAM_mecoco = 456
  INTEGER            ::     P_mecoco = 1
  LOGICAL            ::     F_mecoco = .FALSE.
  INTEGER, PARAMETER :: PARAM_c3h6o2 = 457
  INTEGER            ::     P_c3h6o2 = 1
  LOGICAL            ::     F_c3h6o2 = .FALSE.
  INTEGER, PARAMETER :: PARAM_c4h9o2 = 458
  INTEGER            ::     P_c4h9o2 = 1
  LOGICAL            ::     F_c4h9o2 = .FALSE.
  INTEGER, PARAMETER :: PARAM_meko2 = 459
  INTEGER            ::     P_meko2 = 1
  LOGICAL            ::     F_meko2 = .FALSE.
  INTEGER, PARAMETER :: PARAM_prono2 = 460
  INTEGER            ::     P_prono2 = 1
  LOGICAL            ::     F_prono2 = .FALSE.
  INTEGER, PARAMETER :: PARAM_acetol = 461
  INTEGER            ::     P_acetol = 1
  LOGICAL            ::     F_acetol = .FALSE.
  INTEGER, PARAMETER :: PARAM_acetp = 462
  INTEGER            ::     P_acetp = 1
  LOGICAL            ::     F_acetp = .FALSE.
  INTEGER, PARAMETER :: PARAM_aceto2 = 463
  INTEGER            ::     P_aceto2 = 1
  LOGICAL            ::     F_aceto2 = .FALSE.
  INTEGER, PARAMETER :: PARAM_ch3cooh = 464
  INTEGER            ::     P_ch3cooh = 1
  LOGICAL            ::     F_ch3cooh = .FALSE.
  INTEGER, PARAMETER :: PARAM_c4h9ooh = 465
  INTEGER            ::     P_c4h9ooh = 1
  LOGICAL            ::     F_c4h9ooh = .FALSE.
  INTEGER, PARAMETER :: PARAM_meo2 = 466
  INTEGER            ::     P_meo2 = 1
  LOGICAL            ::     F_meo2 = .FALSE.
  INTEGER, PARAMETER :: PARAM_meoh = 467
  INTEGER            ::     P_meoh = 1
  LOGICAL            ::     F_meoh = .FALSE.
  INTEGER, PARAMETER :: PARAM_meo2no2 = 468
  INTEGER            ::     P_meo2no2 = 1
  LOGICAL            ::     F_meo2no2 = .FALSE.
  INTEGER, PARAMETER :: PARAM_tr1 = 469
  INTEGER            ::     P_tr1 = 1
  LOGICAL            ::     F_tr1 = .FALSE.
  INTEGER, PARAMETER :: PARAM_tr2 = 470
  INTEGER            ::     P_tr2 = 1
  LOGICAL            ::     F_tr2 = .FALSE.
  INTEGER, PARAMETER :: PARAM_NUM_chem = 471
!ggnum INTEGER            ::       NUM_chem = 1
  INTEGER, PARAMETER :: PARAM_tr17_0 = 0
  INTEGER            ::     P_tr17_0 = 1
  LOGICAL            ::     F_tr17_0 = .FALSE.
  INTEGER, PARAMETER :: PARAM_tr17_1 = 1
  INTEGER            ::     P_tr17_1 = 1
  LOGICAL            ::     F_tr17_1 = .FALSE.
  INTEGER, PARAMETER :: PARAM_tr17_2 = 2
  INTEGER            ::     P_tr17_2 = 1
  LOGICAL            ::     F_tr17_2 = .FALSE.
  INTEGER, PARAMETER :: PARAM_tr17_3 = 3
  INTEGER            ::     P_tr17_3 = 1
  LOGICAL            ::     F_tr17_3 = .FALSE.
  INTEGER, PARAMETER :: PARAM_tr17_4 = 4
  INTEGER            ::     P_tr17_4 = 1
  LOGICAL            ::     F_tr17_4 = .FALSE.
  INTEGER, PARAMETER :: PARAM_tr17_5 = 5
  INTEGER            ::     P_tr17_5 = 1
  LOGICAL            ::     F_tr17_5 = .FALSE.
  INTEGER, PARAMETER :: PARAM_tr17_6 = 6
  INTEGER            ::     P_tr17_6 = 1
  LOGICAL            ::     F_tr17_6 = .FALSE.
  INTEGER, PARAMETER :: PARAM_tr17_7 = 7
  INTEGER            ::     P_tr17_7 = 1
  LOGICAL            ::     F_tr17_7 = .FALSE.
  INTEGER, PARAMETER :: PARAM_tr17_8 = 8
  INTEGER            ::     P_tr17_8 = 1
  LOGICAL            ::     F_tr17_8 = .FALSE.
  INTEGER, PARAMETER :: PARAM_tr17_9 = 9
  INTEGER            ::     P_tr17_9 = 1
  LOGICAL            ::     F_tr17_9 = .FALSE.
  INTEGER, PARAMETER :: PARAM_tr18_0 = 10
  INTEGER            ::     P_tr18_0 = 1
  LOGICAL            ::     F_tr18_0 = .FALSE.
  INTEGER, PARAMETER :: PARAM_tr18_1 = 11
  INTEGER            ::     P_tr18_1 = 1
  LOGICAL            ::     F_tr18_1 = .FALSE.
  INTEGER, PARAMETER :: PARAM_tr18_2 = 12
  INTEGER            ::     P_tr18_2 = 1
  LOGICAL            ::     F_tr18_2 = .FALSE.
  INTEGER, PARAMETER :: PARAM_tr18_3 = 13
  INTEGER            ::     P_tr18_3 = 1
  LOGICAL            ::     F_tr18_3 = .FALSE.
  INTEGER, PARAMETER :: PARAM_tr18_4 = 14
  INTEGER            ::     P_tr18_4 = 1
  LOGICAL            ::     F_tr18_4 = .FALSE.
  INTEGER, PARAMETER :: PARAM_tr18_5 = 15
  INTEGER            ::     P_tr18_5 = 1
  LOGICAL            ::     F_tr18_5 = .FALSE.
  INTEGER, PARAMETER :: PARAM_tr18_6 = 16
  INTEGER            ::     P_tr18_6 = 1
  LOGICAL            ::     F_tr18_6 = .FALSE.
  INTEGER, PARAMETER :: PARAM_tr18_7 = 17
  INTEGER            ::     P_tr18_7 = 1
  LOGICAL            ::     F_tr18_7 = .FALSE.
  INTEGER, PARAMETER :: PARAM_tr18_8 = 18
  INTEGER            ::     P_tr18_8 = 1
  LOGICAL            ::     F_tr18_8 = .FALSE.
  INTEGER, PARAMETER :: PARAM_tr18_9 = 19
  INTEGER            ::     P_tr18_9 = 1
  LOGICAL            ::     F_tr18_9 = .FALSE.
  INTEGER, PARAMETER :: PARAM_qni = 21
  INTEGER            ::     P_qni = 1
  LOGICAL            ::     F_qni = .FALSE.
  INTEGER, PARAMETER :: PARAM_qndrop = 22
  INTEGER            ::     P_qndrop = 1
  LOGICAL            ::     F_qndrop = .FALSE.
  INTEGER, PARAMETER :: PARAM_qt = 23
  INTEGER            ::     P_qt = 1
  LOGICAL            ::     F_qt = .FALSE.
  INTEGER, PARAMETER :: PARAM_qns = 24
  INTEGER            ::     P_qns = 1
  LOGICAL            ::     F_qns = .FALSE.
  INTEGER, PARAMETER :: PARAM_qnr = 25
  INTEGER            ::     P_qnr = 1
  LOGICAL            ::     F_qnr = .FALSE.
  INTEGER, PARAMETER :: PARAM_qng = 26
  INTEGER            ::     P_qng = 1
  LOGICAL            ::     F_qng = .FALSE.
  INTEGER, PARAMETER :: PARAM_qnn = 27
  INTEGER            ::     P_qnn = 1
  LOGICAL            ::     F_qnn = .FALSE.
  INTEGER, PARAMETER :: PARAM_qnc = 28
  INTEGER            ::     P_qnc = 1
  LOGICAL            ::     F_qnc = .FALSE.
  INTEGER, PARAMETER :: PARAM_NUM_scalar = 29
!ggnum INTEGER            ::       NUM_scalar = 1
  INTEGER, PARAMETER :: PARAM_dfi_qndrop = 1
  INTEGER            ::     P_dfi_qndrop = 1
  LOGICAL            ::     F_dfi_qndrop = .FALSE.
  INTEGER, PARAMETER :: PARAM_dfi_qni = 2
  INTEGER            ::     P_dfi_qni = 1
  LOGICAL            ::     F_dfi_qni = .FALSE.
  INTEGER, PARAMETER :: PARAM_dfi_qt = 3
  INTEGER            ::     P_dfi_qt = 1
  LOGICAL            ::     F_dfi_qt = .FALSE.
  INTEGER, PARAMETER :: PARAM_dfi_qns = 4
  INTEGER            ::     P_dfi_qns = 1
  LOGICAL            ::     F_dfi_qns = .FALSE.
  INTEGER, PARAMETER :: PARAM_dfi_qnr = 5
  INTEGER            ::     P_dfi_qnr = 1
  LOGICAL            ::     F_dfi_qnr = .FALSE.
  INTEGER, PARAMETER :: PARAM_dfi_qng = 6
  INTEGER            ::     P_dfi_qng = 1
  LOGICAL            ::     F_dfi_qng = .FALSE.
  INTEGER, PARAMETER :: PARAM_dfi_qnn = 7
  INTEGER            ::     P_dfi_qnn = 1
  LOGICAL            ::     F_dfi_qnn = .FALSE.
  INTEGER, PARAMETER :: PARAM_dfi_qnc = 8
  INTEGER            ::     P_dfi_qnc = 1
  LOGICAL            ::     F_dfi_qnc = .FALSE.
  INTEGER, PARAMETER :: PARAM_NUM_dfi_scalar = 9
  INTEGER            ::       NUM_dfi_scalar = 1
  INTEGER, PARAMETER :: PARAM_mth01 = 1
  INTEGER            ::     P_mth01 = 1
  LOGICAL            ::     F_mth01 = .FALSE.
  INTEGER, PARAMETER :: PARAM_mth02 = 2
  INTEGER            ::     P_mth02 = 1
  LOGICAL            ::     F_mth02 = .FALSE.
  INTEGER, PARAMETER :: PARAM_mth03 = 3
  INTEGER            ::     P_mth03 = 1
  LOGICAL            ::     F_mth03 = .FALSE.
  INTEGER, PARAMETER :: PARAM_mth04 = 4
  INTEGER            ::     P_mth04 = 1
  LOGICAL            ::     F_mth04 = .FALSE.
  INTEGER, PARAMETER :: PARAM_mth05 = 5
  INTEGER            ::     P_mth05 = 1
  LOGICAL            ::     F_mth05 = .FALSE.
  INTEGER, PARAMETER :: PARAM_mth06 = 6
  INTEGER            ::     P_mth06 = 1
  LOGICAL            ::     F_mth06 = .FALSE.
  INTEGER, PARAMETER :: PARAM_mth07 = 7
  INTEGER            ::     P_mth07 = 1
  LOGICAL            ::     F_mth07 = .FALSE.
  INTEGER, PARAMETER :: PARAM_mth08 = 8
  INTEGER            ::     P_mth08 = 1
  LOGICAL            ::     F_mth08 = .FALSE.
  INTEGER, PARAMETER :: PARAM_mth09 = 9
  INTEGER            ::     P_mth09 = 1
  LOGICAL            ::     F_mth09 = .FALSE.
  INTEGER, PARAMETER :: PARAM_mth10 = 10
  INTEGER            ::     P_mth10 = 1
  LOGICAL            ::     F_mth10 = .FALSE.
  INTEGER, PARAMETER :: PARAM_mth11 = 11
  INTEGER            ::     P_mth11 = 1
  LOGICAL            ::     F_mth11 = .FALSE.
  INTEGER, PARAMETER :: PARAM_mth12 = 12
  INTEGER            ::     P_mth12 = 1
  LOGICAL            ::     F_mth12 = .FALSE.
  INTEGER, PARAMETER :: PARAM_NUM_ozmixm = 13
  INTEGER            ::       NUM_ozmixm = 1
  INTEGER, PARAMETER :: PARAM_sul = 1
  INTEGER            ::     P_sul = 1
  LOGICAL            ::     F_sul = .FALSE.
  INTEGER, PARAMETER :: PARAM_sslt = 2
  INTEGER            ::     P_sslt = 1
  LOGICAL            ::     F_sslt = .FALSE.
  INTEGER, PARAMETER :: PARAM_dust1 = 3
  INTEGER            ::     P_dust1 = 1
  LOGICAL            ::     F_dust1 = .FALSE.
  INTEGER, PARAMETER :: PARAM_dust2 = 4
  INTEGER            ::     P_dust2 = 1
  LOGICAL            ::     F_dust2 = .FALSE.
  INTEGER, PARAMETER :: PARAM_dust3 = 5
  INTEGER            ::     P_dust3 = 1
  LOGICAL            ::     F_dust3 = .FALSE.
  INTEGER, PARAMETER :: PARAM_dust4 = 6
  INTEGER            ::     P_dust4 = 1
  LOGICAL            ::     F_dust4 = .FALSE.
  INTEGER, PARAMETER :: PARAM_ocpho = 7
  INTEGER            ::     P_ocpho = 1
  LOGICAL            ::     F_ocpho = .FALSE.
  INTEGER, PARAMETER :: PARAM_bcpho = 8
  INTEGER            ::     P_bcpho = 1
  LOGICAL            ::     F_bcpho = .FALSE.
  INTEGER, PARAMETER :: PARAM_ocphi = 9
  INTEGER            ::     P_ocphi = 1
  LOGICAL            ::     F_ocphi = .FALSE.
  INTEGER, PARAMETER :: PARAM_bcphi = 10
  INTEGER            ::     P_bcphi = 1
  LOGICAL            ::     F_bcphi = .FALSE.
  INTEGER, PARAMETER :: PARAM_bg = 11
  INTEGER            ::     P_bg = 1
  LOGICAL            ::     F_bg = .FALSE.
  INTEGER, PARAMETER :: PARAM_volc = 12
  INTEGER            ::     P_volc = 1
  LOGICAL            ::     F_volc = .FALSE.
  INTEGER, PARAMETER :: PARAM_NUM_aerosolc = 13
  INTEGER            ::       NUM_aerosolc = 1
  INTEGER, PARAMETER :: PARAM_u_ndg_new = 1
  INTEGER            ::     P_u_ndg_new = 1
  LOGICAL            ::     F_u_ndg_new = .FALSE.
  INTEGER, PARAMETER :: PARAM_v_ndg_new = 2
  INTEGER            ::     P_v_ndg_new = 1
  LOGICAL            ::     F_v_ndg_new = .FALSE.
  INTEGER, PARAMETER :: PARAM_t_ndg_new = 3
  INTEGER            ::     P_t_ndg_new = 1
  LOGICAL            ::     F_t_ndg_new = .FALSE.
  INTEGER, PARAMETER :: PARAM_q_ndg_new = 4
  INTEGER            ::     P_q_ndg_new = 1
  LOGICAL            ::     F_q_ndg_new = .FALSE.
  INTEGER, PARAMETER :: PARAM_ph_ndg_new = 5
  INTEGER            ::     P_ph_ndg_new = 1
  LOGICAL            ::     F_ph_ndg_new = .FALSE.
  INTEGER, PARAMETER :: PARAM_u_ndg_old = 6
  INTEGER            ::     P_u_ndg_old = 1
  LOGICAL            ::     F_u_ndg_old = .FALSE.
  INTEGER, PARAMETER :: PARAM_v_ndg_old = 7
  INTEGER            ::     P_v_ndg_old = 1
  LOGICAL            ::     F_v_ndg_old = .FALSE.
  INTEGER, PARAMETER :: PARAM_t_ndg_old = 8
  INTEGER            ::     P_t_ndg_old = 1
  LOGICAL            ::     F_t_ndg_old = .FALSE.
  INTEGER, PARAMETER :: PARAM_q_ndg_old = 9
  INTEGER            ::     P_q_ndg_old = 1
  LOGICAL            ::     F_q_ndg_old = .FALSE.
  INTEGER, PARAMETER :: PARAM_ph_ndg_old = 10
  INTEGER            ::     P_ph_ndg_old = 1
  LOGICAL            ::     F_ph_ndg_old = .FALSE.
  INTEGER, PARAMETER :: PARAM_NUM_fdda3d = 11
  INTEGER            ::       NUM_fdda3d = 1
  INTEGER, PARAMETER :: PARAM_mu_ndg_new = 1
  INTEGER            ::     P_mu_ndg_new = 1
  LOGICAL            ::     F_mu_ndg_new = .FALSE.
  INTEGER, PARAMETER :: PARAM_mu_ndg_old = 2
  INTEGER            ::     P_mu_ndg_old = 1
  LOGICAL            ::     F_mu_ndg_old = .FALSE.
  INTEGER, PARAMETER :: PARAM_NUM_fdda2d = 3
  INTEGER            ::       NUM_fdda2d = 1
  INTEGER, PARAMETER :: P_XSB                          = 1
  INTEGER, PARAMETER :: P_XEB                          = 2
  INTEGER, PARAMETER :: P_YSB                          = 3
  INTEGER, PARAMETER :: P_YEB                          = 4
  INTEGER, PARAMETER :: NUM_TIME_LEVELS = 2
  INTEGER , PARAMETER :: PARAM_FIRST_SCALAR = 2
!ENDOFREGISTRYGENERATEDINCLUDE
! That was all of state_description V3.1!!!!!!
   TYPE grid_config_rec_type
     SEQUENCE
!#include <namelist_defines2.inc> V3.1 - Mar10, 2009 !GG

!STARTOFREGISTRYGENERATEDINCLUDE 'inc/namelist_defines2.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
integer    :: first_item_in_struct
character*256 :: emi_inname
character*256 :: fireemi_inname
character*256 :: input_chem_inname
character*256 :: emi_outname
character*256 :: fireemi_outname
character*256 :: input_chem_outname
integer :: frames_per_emissfile
integer :: frames_per_fireemissfile
integer :: io_style_emissions
integer :: io_form_emissions
integer :: io_style_fireemissions
integer :: io_form_fireemissions
real :: bioemdt
real :: photdt
real :: chemdt
integer :: ne_area
integer :: kemit
integer :: nmegan
integer :: kfuture
integer :: errosion_dim
integer :: biomass_emiss_opt
integer :: chem_conv_tr
integer :: chem_opt
integer :: gaschem_onoff
integer :: aerchem_onoff
integer :: wetscav_onoff
integer :: cldchem_onoff
integer :: vertmix_onoff
integer :: chem_in_opt
integer :: phot_opt
integer :: drydep_opt
integer :: emiss_opt
integer :: dust_opt
integer :: dmsemis_opt
integer :: seas_opt
integer :: bio_emiss_opt
integer :: biomass_burn_opt
integer :: plumerisefire_frq
integer :: emiss_inpt_opt
integer :: gas_bc_opt
integer :: gas_ic_opt
integer :: aer_bc_opt
integer :: aer_ic_opt
logical :: have_bcs_chem
integer :: aer_ra_feedback
integer :: aer_op_opt
integer :: scalar_opt
integer :: run_days
integer :: run_hours
integer :: run_minutes
integer :: run_seconds
integer :: start_year
integer :: start_month
integer :: start_day
integer :: start_hour
integer :: start_minute
integer :: start_second
integer :: end_year
integer :: end_month
integer :: end_day
integer :: end_hour
integer :: end_minute
integer :: end_second
integer :: interval_seconds
logical :: input_from_file
integer :: fine_input_stream
logical :: input_from_hires
character*256 :: rsmas_data_path
logical :: all_ic_times
integer :: history_interval
integer :: frames_per_outfile
integer :: frames_per_auxhist1
integer :: frames_per_auxhist2
integer :: frames_per_auxhist3
integer :: frames_per_auxhist4
integer :: frames_per_auxhist5
integer :: frames_per_auxhist6
integer :: frames_per_auxhist7
integer :: frames_per_auxhist8
integer :: frames_per_auxhist9
integer :: frames_per_auxhist10
integer :: frames_per_auxhist11
logical :: restart
integer :: restart_interval
integer :: io_form_input
integer :: io_form_history
integer :: io_form_restart
integer :: io_form_boundary
integer :: debug_level
logical :: self_test_domain
character*256 :: history_outname
character*256 :: auxhist1_outname
character*256 :: auxhist2_outname
character*256 :: auxhist3_outname
character*256 :: auxhist4_outname
character*256 :: auxhist5_outname
character*256 :: auxhist6_outname
character*256 :: auxhist7_outname
character*256 :: auxhist8_outname
character*256 :: auxhist9_outname
character*256 :: auxhist10_outname
character*256 :: auxhist11_outname
character*256 :: history_inname
character*256 :: auxhist1_inname
character*256 :: auxhist2_inname
character*256 :: auxhist3_inname
character*256 :: auxhist4_inname
character*256 :: auxhist5_inname
character*256 :: auxhist6_inname
character*256 :: auxhist7_inname
character*256 :: auxhist8_inname
character*256 :: auxhist9_inname
character*256 :: auxhist10_inname
character*256 :: auxhist11_inname
character*256 :: auxinput1_outname
character*256 :: auxinput2_outname
character*256 :: auxinput3_outname
character*256 :: auxinput4_outname
character*256 :: auxinput5_outname
character*256 :: auxinput6_outname
character*256 :: auxinput7_outname
character*256 :: auxinput8_outname
character*256 :: auxinput9_outname
character*256 :: auxinput10_outname
character*256 :: auxinput11_outname
character*256 :: auxinput1_inname
character*256 :: auxinput2_inname
character*256 :: auxinput3_inname
character*256 :: auxinput4_inname
character*256 :: auxinput5_inname
character*256 :: auxinput6_inname
character*256 :: auxinput7_inname
character*256 :: auxinput8_inname
character*256 :: sgfdda_inname
character*256 :: gfdda_inname
character*256 :: auxinput11_inname
integer :: history_interval_mo
integer :: history_interval_d
integer :: history_interval_h
integer :: history_interval_m
integer :: history_interval_s
integer :: inputout_interval_mo
integer :: inputout_interval_d
integer :: inputout_interval_h
integer :: inputout_interval_m
integer :: inputout_interval_s
integer :: inputout_interval
integer :: auxhist1_interval_mo
integer :: auxhist1_interval_d
integer :: auxhist1_interval_h
integer :: auxhist1_interval_m
integer :: auxhist1_interval_s
integer :: auxhist1_interval
integer :: auxhist2_interval_mo
integer :: auxhist2_interval_d
integer :: auxhist2_interval_h
integer :: auxhist2_interval_m
integer :: auxhist2_interval_s
integer :: auxhist2_interval
integer :: auxhist3_interval_mo
integer :: auxhist3_interval_d
integer :: auxhist3_interval_h
integer :: auxhist3_interval_m
integer :: auxhist3_interval_s
integer :: auxhist3_interval
integer :: auxhist4_interval_mo
integer :: auxhist4_interval_d
integer :: auxhist4_interval_h
integer :: auxhist4_interval_m
integer :: auxhist4_interval_s
integer :: auxhist4_interval
integer :: auxhist5_interval_mo
integer :: auxhist5_interval_d
integer :: auxhist5_interval_h
integer :: auxhist5_interval_m
integer :: auxhist5_interval_s
integer :: auxhist5_interval
integer :: auxhist6_interval_mo
integer :: auxhist6_interval_d
integer :: auxhist6_interval_h
integer :: auxhist6_interval_m
integer :: auxhist6_interval_s
integer :: auxhist6_interval
integer :: auxhist7_interval_mo
integer :: auxhist7_interval_d
integer :: auxhist7_interval_h
integer :: auxhist7_interval_m
integer :: auxhist7_interval_s
integer :: auxhist7_interval
integer :: auxhist8_interval_mo
integer :: auxhist8_interval_d
integer :: auxhist8_interval_h
integer :: auxhist8_interval_m
integer :: auxhist8_interval_s
integer :: auxhist8_interval
integer :: auxhist9_interval_mo
integer :: auxhist9_interval_d
integer :: auxhist9_interval_h
integer :: auxhist9_interval_m
integer :: auxhist9_interval_s
integer :: auxhist9_interval
integer :: auxhist10_interval_mo
integer :: auxhist10_interval_d
integer :: auxhist10_interval_h
integer :: auxhist10_interval_m
integer :: auxhist10_interval_s
integer :: auxhist10_interval
integer :: auxhist11_interval_mo
integer :: auxhist11_interval_d
integer :: auxhist11_interval_h
integer :: auxhist11_interval_m
integer :: auxhist11_interval_s
integer :: auxhist11_interval
integer :: auxinput1_interval_mo
integer :: auxinput1_interval_d
integer :: auxinput1_interval_h
integer :: auxinput1_interval_m
integer :: auxinput1_interval_s
integer :: auxinput1_interval
integer :: auxinput2_interval_mo
integer :: auxinput2_interval_d
integer :: auxinput2_interval_h
integer :: auxinput2_interval_m
integer :: auxinput2_interval_s
integer :: auxinput2_interval
integer :: auxinput3_interval_mo
integer :: auxinput3_interval_d
integer :: auxinput3_interval_h
integer :: auxinput3_interval_m
integer :: auxinput3_interval_s
integer :: auxinput3_interval
integer :: auxinput4_interval_mo
integer :: auxinput4_interval_d
integer :: auxinput4_interval_h
integer :: auxinput4_interval_m
integer :: auxinput4_interval_s
integer :: auxinput4_interval
integer :: auxinput5_interval_mo
integer :: auxinput5_interval_d
integer :: auxinput5_interval_h
integer :: auxinput5_interval_m
integer :: auxinput5_interval_s
integer :: auxinput5_interval
integer :: auxinput6_interval_mo
integer :: auxinput6_interval_d
integer :: auxinput6_interval_h
integer :: auxinput6_interval_m
integer :: auxinput6_interval_s
integer :: auxinput6_interval
integer :: auxinput7_interval_mo
integer :: auxinput7_interval_d
integer :: auxinput7_interval_h
integer :: auxinput7_interval_m
integer :: auxinput7_interval_s
integer :: auxinput7_interval
integer :: auxinput8_interval_mo
integer :: auxinput8_interval_d
integer :: auxinput8_interval_h
integer :: auxinput8_interval_m
integer :: auxinput8_interval_s
integer :: auxinput8_interval
integer :: sgfdda_interval_mo
integer :: sgfdda_interval_d
integer :: sgfdda_interval_h
integer :: sgfdda_interval_m
integer :: sgfdda_interval_s
integer :: sgfdda_interval
integer :: gfdda_interval_mo
integer :: gfdda_interval_d
integer :: gfdda_interval_h
integer :: gfdda_interval_m
integer :: gfdda_interval_s
integer :: gfdda_interval
integer :: auxinput11_interval_mo
integer :: auxinput11_interval_d
integer :: auxinput11_interval_h
integer :: auxinput11_interval_m
integer :: auxinput11_interval_s
integer :: auxinput11_interval
integer :: restart_interval_mo
integer :: restart_interval_d
integer :: restart_interval_h
integer :: restart_interval_m
integer :: restart_interval_s
integer :: history_begin_y
integer :: history_begin_mo
integer :: history_begin_d
integer :: history_begin_h
integer :: history_begin_m
integer :: history_begin_s
integer :: inputout_begin_y
integer :: inputout_begin_mo
integer :: inputout_begin_d
integer :: inputout_begin_h
integer :: inputout_begin_m
integer :: inputout_begin_s
integer :: auxhist1_begin_y
integer :: auxhist1_begin_mo
integer :: auxhist1_begin_d
integer :: auxhist1_begin_h
integer :: auxhist1_begin_m
integer :: auxhist1_begin_s
integer :: auxhist2_begin_y
integer :: auxhist2_begin_mo
integer :: auxhist2_begin_d
integer :: auxhist2_begin_h
integer :: auxhist2_begin_m
integer :: auxhist2_begin_s
integer :: auxhist3_begin_y
integer :: auxhist3_begin_mo
integer :: auxhist3_begin_d
integer :: auxhist3_begin_h
integer :: auxhist3_begin_m
integer :: auxhist3_begin_s
integer :: auxhist4_begin_y
integer :: auxhist4_begin_mo
integer :: auxhist4_begin_d
integer :: auxhist4_begin_h
integer :: auxhist4_begin_m
integer :: auxhist4_begin_s
integer :: auxhist5_begin_y
integer :: auxhist5_begin_mo
integer :: auxhist5_begin_d
integer :: auxhist5_begin_h
integer :: auxhist5_begin_m
integer :: auxhist5_begin_s
integer :: auxhist6_begin_y
integer :: auxhist6_begin_mo
integer :: auxhist6_begin_d
integer :: auxhist6_begin_h
integer :: auxhist6_begin_m
integer :: auxhist6_begin_s
integer :: auxhist7_begin_y
integer :: auxhist7_begin_mo
integer :: auxhist7_begin_d
integer :: auxhist7_begin_h
integer :: auxhist7_begin_m
integer :: auxhist7_begin_s
integer :: auxhist8_begin_y
integer :: auxhist8_begin_mo
integer :: auxhist8_begin_d
integer :: auxhist8_begin_h
integer :: auxhist8_begin_m
integer :: auxhist8_begin_s
integer :: auxhist9_begin_y
integer :: auxhist9_begin_mo
integer :: auxhist9_begin_d
integer :: auxhist9_begin_h
integer :: auxhist9_begin_m
integer :: auxhist9_begin_s
integer :: auxhist10_begin_y
integer :: auxhist10_begin_mo
integer :: auxhist10_begin_d
integer :: auxhist10_begin_h
integer :: auxhist10_begin_m
integer :: auxhist10_begin_s
integer :: auxhist11_begin_y
integer :: auxhist11_begin_mo
integer :: auxhist11_begin_d
integer :: auxhist11_begin_h
integer :: auxhist11_begin_m
integer :: auxhist11_begin_s
integer :: auxinput1_begin_y
integer :: auxinput1_begin_mo
integer :: auxinput1_begin_d
integer :: auxinput1_begin_h
integer :: auxinput1_begin_m
integer :: auxinput1_begin_s
integer :: auxinput2_begin_y
integer :: auxinput2_begin_mo
integer :: auxinput2_begin_d
integer :: auxinput2_begin_h
integer :: auxinput2_begin_m
integer :: auxinput2_begin_s
integer :: auxinput3_begin_y
integer :: auxinput3_begin_mo
integer :: auxinput3_begin_d
integer :: auxinput3_begin_h
integer :: auxinput3_begin_m
integer :: auxinput3_begin_s
integer :: auxinput4_begin_y
integer :: auxinput4_begin_mo
integer :: auxinput4_begin_d
integer :: auxinput4_begin_h
integer :: auxinput4_begin_m
integer :: auxinput4_begin_s
integer :: auxinput5_begin_y
integer :: auxinput5_begin_mo
integer :: auxinput5_begin_d
integer :: auxinput5_begin_h
integer :: auxinput5_begin_m
integer :: auxinput5_begin_s
integer :: auxinput6_begin_y
integer :: auxinput6_begin_mo
integer :: auxinput6_begin_d
integer :: auxinput6_begin_h
integer :: auxinput6_begin_m
integer :: auxinput6_begin_s
integer :: auxinput7_begin_y
integer :: auxinput7_begin_mo
integer :: auxinput7_begin_d
integer :: auxinput7_begin_h
integer :: auxinput7_begin_m
integer :: auxinput7_begin_s
integer :: auxinput8_begin_y
integer :: auxinput8_begin_mo
integer :: auxinput8_begin_d
integer :: auxinput8_begin_h
integer :: auxinput8_begin_m
integer :: auxinput8_begin_s
integer :: sgfdda_begin_y
integer :: sgfdda_begin_mo
integer :: sgfdda_begin_d
integer :: sgfdda_begin_h
integer :: sgfdda_begin_m
integer :: sgfdda_begin_s
integer :: gfdda_begin_y
integer :: gfdda_begin_mo
integer :: gfdda_begin_d
integer :: gfdda_begin_h
integer :: gfdda_begin_m
integer :: gfdda_begin_s
integer :: auxinput11_begin_y
integer :: auxinput11_begin_mo
integer :: auxinput11_begin_d
integer :: auxinput11_begin_h
integer :: auxinput11_begin_m
integer :: auxinput11_begin_s
integer :: restart_begin_y
integer :: restart_begin_mo
integer :: restart_begin_d
integer :: restart_begin_h
integer :: restart_begin_m
integer :: restart_begin_s
integer :: history_end_y
integer :: history_end_mo
integer :: history_end_d
integer :: history_end_h
integer :: history_end_m
integer :: history_end_s
integer :: inputout_end_y
integer :: inputout_end_mo
integer :: inputout_end_d
integer :: inputout_end_h
integer :: inputout_end_m
integer :: inputout_end_s
integer :: auxhist1_end_y
integer :: auxhist1_end_mo
integer :: auxhist1_end_d
integer :: auxhist1_end_h
integer :: auxhist1_end_m
integer :: auxhist1_end_s
integer :: auxhist2_end_y
integer :: auxhist2_end_mo
integer :: auxhist2_end_d
integer :: auxhist2_end_h
integer :: auxhist2_end_m
integer :: auxhist2_end_s
integer :: auxhist3_end_y
integer :: auxhist3_end_mo
integer :: auxhist3_end_d
integer :: auxhist3_end_h
integer :: auxhist3_end_m
integer :: auxhist3_end_s
integer :: auxhist4_end_y
integer :: auxhist4_end_mo
integer :: auxhist4_end_d
integer :: auxhist4_end_h
integer :: auxhist4_end_m
integer :: auxhist4_end_s
integer :: auxhist5_end_y
integer :: auxhist5_end_mo
integer :: auxhist5_end_d
integer :: auxhist5_end_h
integer :: auxhist5_end_m
integer :: auxhist5_end_s
integer :: auxhist6_end_y
integer :: auxhist6_end_mo
integer :: auxhist6_end_d
integer :: auxhist6_end_h
integer :: auxhist6_end_m
integer :: auxhist6_end_s
integer :: auxhist7_end_y
integer :: auxhist7_end_mo
integer :: auxhist7_end_d
integer :: auxhist7_end_h
integer :: auxhist7_end_m
integer :: auxhist7_end_s
integer :: auxhist8_end_y
integer :: auxhist8_end_mo
integer :: auxhist8_end_d
integer :: auxhist8_end_h
integer :: auxhist8_end_m
integer :: auxhist8_end_s
integer :: auxhist9_end_y
integer :: auxhist9_end_mo
integer :: auxhist9_end_d
integer :: auxhist9_end_h
integer :: auxhist9_end_m
integer :: auxhist9_end_s
integer :: auxhist10_end_y
integer :: auxhist10_end_mo
integer :: auxhist10_end_d
integer :: auxhist10_end_h
integer :: auxhist10_end_m
integer :: auxhist10_end_s
integer :: auxhist11_end_y
integer :: auxhist11_end_mo
integer :: auxhist11_end_d
integer :: auxhist11_end_h
integer :: auxhist11_end_m
integer :: auxhist11_end_s
integer :: auxinput1_end_y
integer :: auxinput1_end_mo
integer :: auxinput1_end_d
integer :: auxinput1_end_h
integer :: auxinput1_end_m
integer :: auxinput1_end_s
integer :: auxinput2_end_y
integer :: auxinput2_end_mo
integer :: auxinput2_end_d
integer :: auxinput2_end_h
integer :: auxinput2_end_m
integer :: auxinput2_end_s
integer :: auxinput3_end_y
integer :: auxinput3_end_mo
integer :: auxinput3_end_d
integer :: auxinput3_end_h
integer :: auxinput3_end_m
integer :: auxinput3_end_s
integer :: auxinput4_end_y
integer :: auxinput4_end_mo
integer :: auxinput4_end_d
integer :: auxinput4_end_h
integer :: auxinput4_end_m
integer :: auxinput4_end_s
integer :: auxinput5_end_y
integer :: auxinput5_end_mo
integer :: auxinput5_end_d
integer :: auxinput5_end_h
integer :: auxinput5_end_m
integer :: auxinput5_end_s
integer :: auxinput6_end_y
integer :: auxinput6_end_mo
integer :: auxinput6_end_d
integer :: auxinput6_end_h
integer :: auxinput6_end_m
integer :: auxinput6_end_s
integer :: auxinput7_end_y
integer :: auxinput7_end_mo
integer :: auxinput7_end_d
integer :: auxinput7_end_h
integer :: auxinput7_end_m
integer :: auxinput7_end_s
integer :: auxinput8_end_y
integer :: auxinput8_end_mo
integer :: auxinput8_end_d
integer :: auxinput8_end_h
integer :: auxinput8_end_m
integer :: auxinput8_end_s
integer :: sgfdda_end_y
integer :: sgfdda_end_mo
integer :: sgfdda_end_d
integer :: sgfdda_end_h
integer :: sgfdda_end_m
integer :: sgfdda_end_s
integer :: gfdda_end_y
integer :: gfdda_end_mo
integer :: gfdda_end_d
integer :: gfdda_end_h
integer :: gfdda_end_m
integer :: gfdda_end_s
integer :: auxinput11_end_y
integer :: auxinput11_end_mo
integer :: auxinput11_end_d
integer :: auxinput11_end_h
integer :: auxinput11_end_m
integer :: auxinput11_end_s
integer :: io_form_auxinput1
integer :: io_form_auxinput2
integer :: io_form_auxinput3
integer :: io_form_auxinput4
integer :: io_form_auxinput5
integer :: io_form_auxinput6
integer :: io_form_auxinput7
integer :: io_form_auxinput8
integer :: io_form_sgfdda
integer :: io_form_gfdda
integer :: io_form_auxinput11
integer :: io_form_auxhist1
integer :: io_form_auxhist2
integer :: io_form_auxhist3
integer :: io_form_auxhist4
integer :: io_form_auxhist5
integer :: io_form_auxhist6
integer :: io_form_auxhist7
integer :: io_form_auxhist8
integer :: io_form_auxhist9
integer :: io_form_auxhist10
integer :: io_form_auxhist11
integer :: simulation_start_year
integer :: simulation_start_month
integer :: simulation_start_day
integer :: simulation_start_hour
integer :: simulation_start_minute
integer :: simulation_start_second
logical :: reset_simulation_start
integer :: sr_x
integer :: sr_y
integer :: julyr
integer :: julday
real :: gmt
character*256 :: input_inname
character*256 :: input_outname
character*256 :: bdy_inname
character*256 :: bdy_outname
character*256 :: rst_inname
character*256 :: rst_outname
logical :: write_input
logical :: write_restart_at_0h
logical :: adjust_output_times
logical :: adjust_input_times
integer :: diag_print
logical :: nocolons
integer :: dfi_opt
integer :: dfi_nfilter
logical :: dfi_write_filtered_input
logical :: dfi_write_dfi_history
integer :: dfi_cutoff_seconds
integer :: dfi_time_dim
integer :: dfi_fwdstop_year
integer :: dfi_fwdstop_month
integer :: dfi_fwdstop_day
integer :: dfi_fwdstop_hour
integer :: dfi_fwdstop_minute
integer :: dfi_fwdstop_second
integer :: dfi_bckstop_year
integer :: dfi_bckstop_month
integer :: dfi_bckstop_day
integer :: dfi_bckstop_hour
integer :: dfi_bckstop_minute
integer :: dfi_bckstop_second
integer :: time_step
integer :: time_step_fract_num
integer :: time_step_fract_den
integer :: min_time_step
integer :: max_time_step
real :: target_cfl
integer :: max_step_increase_pct
integer :: starting_time_step
logical :: step_to_output_time
logical :: use_adaptive_time_step
integer :: max_dom
integer :: s_we
integer :: e_we
integer :: s_sn
integer :: e_sn
integer :: s_vert
integer :: e_vert
integer :: num_metgrid_levels
integer :: num_soil_layers_in
real :: p_top_requested
integer :: interp_type
integer :: extrap_type
integer :: t_extrap_type
logical :: lowest_lev_from_sfc
logical :: use_levels_below_ground
logical :: use_surface
integer :: lagrange_order
integer :: force_sfc_in_vinterp
real :: zap_close_levels
logical :: sfcp_to_sfcp
logical :: adjust_heights
logical :: smooth_cg_topo
logical :: rh2qv_wrt_liquid
real :: dx
real :: dy
integer :: grid_id
logical :: grid_allowed
integer :: parent_id
integer :: i_parent_start
integer :: j_parent_start
integer :: parent_grid_ratio
integer :: parent_time_step_ratio
integer :: feedback
integer :: smooth_option
integer :: blend_width
real :: ztop
integer :: moad_grid_ratio
integer :: moad_time_step_ratio
integer :: shw
integer :: tile_sz_x
integer :: tile_sz_y
integer :: numtiles
integer :: nproc_x
integer :: nproc_y
integer :: irand
integer :: num_moves
integer :: ts_buf_size
integer :: max_ts_locs
integer :: vortex_interval
integer :: max_vortex_speed
integer :: corral_dist
integer :: track_level
integer :: move_id
integer :: move_interval
integer :: move_cd_x
integer :: move_cd_y
logical :: swap_x
logical :: swap_y
logical :: cycle_x
logical :: cycle_y
logical :: reorder_mesh
logical :: perturb_input
real :: eta_levels
real :: max_dz
logical :: insert_bogus_storm
integer :: num_storm
real :: latc_loc
real :: lonc_loc
real :: vmax_meters_per_second
real :: rmax
real :: vmax_ratio
integer :: mp_physics
integer :: gsfcgce_hail
integer :: gsfcgce_2ice
integer :: progn
integer :: ra_lw_physics
integer :: ra_sw_physics
real :: radt
real :: naer
integer :: sf_sfclay_physics
integer :: sf_surface_physics
integer :: bl_pbl_physics
integer :: sf_urban_physics
real :: bldt
integer :: cu_physics
real :: cudt
real :: gsmdt
integer :: isfflx
integer :: ifsnow
integer :: icloud
real :: swrad_scat
integer :: surface_input_source
!TBH:  avoid conflict with declaration in module_control
!TBH integer :: num_soil_layers
integer :: num_urban_layers
integer :: num_months
integer :: maxiens
integer :: maxens
integer :: maxens2
integer :: maxens3
integer :: ensdim
integer :: cugd_avedx
integer :: imomentum
integer :: clos_choice
integer :: num_land_cat
integer :: num_soil_cat
integer :: mp_zero_out
real :: mp_zero_out_thresh
real :: seaice_threshold
integer :: sst_update
integer :: sst_skin
integer :: tmn_update
logical :: usemonalb
logical :: rdmaxalb
logical :: rdlai2d
integer :: co2tf
integer :: ra_call_offset
real :: cam_abs_freq_s
integer :: levsiz
integer :: paerlev
integer :: cam_abs_dim1
integer :: cam_abs_dim2
integer :: lagday
logical :: cu_rad_feedback
integer :: pxlsm_smois_init
integer :: omlcall
real :: oml_hml0
real :: oml_gamma
integer :: isftcflx
real :: shadlen
integer :: slope_rad
integer :: topo_shading
integer :: no_mp_heating
integer :: fractional_seaice
real :: bucket_mm
real :: bucket_j
integer :: grav_settling
real :: fgdt
integer :: fgdtzero
integer :: grid_fdda
integer :: grid_sfdda
integer :: if_no_pbl_nudging_uv
integer :: if_no_pbl_nudging_t
integer :: if_no_pbl_nudging_ph
integer :: if_no_pbl_nudging_q
integer :: if_zfac_uv
integer :: k_zfac_uv
integer :: if_zfac_t
integer :: k_zfac_t
integer :: if_zfac_ph
integer :: k_zfac_ph
integer :: if_zfac_q
integer :: k_zfac_q
integer :: dk_zfac_uv
integer :: dk_zfac_t
integer :: dk_zfac_ph
real :: guv
real :: guv_sfc
real :: gt
real :: gt_sfc
real :: gq
real :: gq_sfc
real :: gph
real :: dtramp_min
integer :: if_ramping
real :: rinblw
integer :: xwavenum
integer :: ywavenum
integer :: obs_nudge_opt
integer :: max_obs
real :: fdda_start
real :: fdda_end
integer :: obs_nudge_wind
real :: obs_coef_wind
integer :: obs_nudge_temp
real :: obs_coef_temp
integer :: obs_nudge_mois
real :: obs_coef_mois
integer :: obs_nudge_pstr
real :: obs_coef_pstr
real :: obs_rinxy
real :: obs_rinsig
real :: obs_twindo
integer :: obs_npfi
integer :: obs_ionf
integer :: obs_idynin
real :: obs_dtramp
integer :: obs_prt_max
integer :: obs_prt_freq
logical :: obs_ipf_in4dob
logical :: obs_ipf_errob
logical :: obs_ipf_nudob
logical :: obs_ipf_init
integer :: scm_force
real :: scm_force_dx
integer :: num_force_layers
integer :: scm_lu_index
integer :: scm_isltyp
real :: scm_vegfra
integer :: scm_canwat
real :: scm_lat
real :: scm_lon
logical :: scm_th_adv
logical :: scm_wind_adv
logical :: scm_qv_adv
logical :: scm_vert_adv
integer :: rk_ord
integer :: w_damping
integer :: diff_opt
integer :: km_opt
integer :: km_opt_dfi
integer :: damp_opt
integer :: gwd_opt
real :: zdamp
real :: dampcoef
real :: khdif
real :: kvdif
real :: diff_6th_factor
integer :: diff_6th_opt
real :: c_s
real :: c_k
real :: smdiv
real :: emdiv
real :: epssm
logical :: non_hydrostatic
integer :: time_step_sound
integer :: h_mom_adv_order
integer :: v_mom_adv_order
integer :: h_sca_adv_order
integer :: v_sca_adv_order
integer :: moist_adv_opt
integer :: moist_adv_dfi_opt
integer :: chem_adv_opt
integer :: scalar_adv_opt
integer :: tke_adv_opt
logical :: top_radiation
integer :: mix_isotropic
real :: mix_upper_bound
logical :: top_lid
real :: tke_upper_bound
real :: tke_drag_coefficient
real :: tke_heat_flux
logical :: pert_coriolis
logical :: coriolis2d
logical :: mix_full_fields
real :: base_pres
real :: base_temp
real :: base_lapse
real :: iso_temp
real :: fft_filter_lat
logical :: rotated_pole
logical :: do_coriolis
logical :: do_curvature
logical :: do_gradp
integer :: spec_bdy_width
integer :: spec_zone
integer :: relax_zone
logical :: specified
logical :: periodic_x
logical :: symmetric_xs
logical :: symmetric_xe
logical :: open_xs
logical :: open_xe
logical :: periodic_y
logical :: symmetric_ys
logical :: symmetric_ye
logical :: open_ys
logical :: open_ye
logical :: polar
logical :: nested
real :: spec_exp
integer :: real_data_init_type
integer :: background_proc_id
integer :: forecast_proc_id
integer :: production_status
integer :: compression
integer :: nobs_ndg_vars
integer :: nobs_err_flds
real :: cen_lat
real :: cen_lon
real :: truelat1
real :: truelat2
real :: moad_cen_lat
real :: stand_lon
real :: bdyfrq
character*256 :: mminlu
real :: emifrq
integer :: iswater
integer :: islake
integer :: isice
integer :: isurban
integer :: isoilwater
integer :: map_proj
integer :: use_wps_input
integer :: dfi_stage
integer :: mp_physics_dfi
integer :: ifire
integer :: fire_boundary_guard
integer :: fire_num_ignitions
real :: fire_ignition_start_long1
real :: fire_ignition_start_lat1
real :: fire_ignition_end_long1
real :: fire_ignition_end_lat1
real :: fire_ignition_radius1
real :: fire_ignition_time1
real :: fire_ignition_start_long2
real :: fire_ignition_start_lat2
real :: fire_ignition_end_long2
real :: fire_ignition_end_lat2
real :: fire_ignition_radius2
real :: fire_ignition_time2
real :: fire_ignition_start_long3
real :: fire_ignition_start_lat3
real :: fire_ignition_end_long3
real :: fire_ignition_end_lat3
real :: fire_ignition_radius3
real :: fire_ignition_time3
real :: fire_ignition_start_long4
real :: fire_ignition_start_lat4
real :: fire_ignition_end_long4
real :: fire_ignition_end_lat4
real :: fire_ignition_radius4
real :: fire_ignition_time4
real :: fire_ignition_start_long5
real :: fire_ignition_start_lat5
real :: fire_ignition_end_long5
real :: fire_ignition_end_lat5
real :: fire_ignition_radius5
real :: fire_ignition_time5
real :: fire_ignition_start_x1
real :: fire_ignition_start_y1
real :: fire_ignition_end_x1
real :: fire_ignition_end_y1
real :: fire_ignition_start_x2
real :: fire_ignition_start_y2
real :: fire_ignition_end_x2
real :: fire_ignition_end_y2
real :: fire_ignition_start_x3
real :: fire_ignition_start_y3
real :: fire_ignition_end_x3
real :: fire_ignition_end_y3
real :: fire_ignition_start_x4
real :: fire_ignition_start_y4
real :: fire_ignition_end_x4
real :: fire_ignition_end_y4
real :: fire_ignition_start_x5
real :: fire_ignition_start_y5
real :: fire_ignition_end_x5
real :: fire_ignition_end_y5
real :: fire_lat_init
real :: fire_lon_init
real :: fire_ign_time
integer :: fire_shape
integer :: fire_sprd_mdl
real :: fire_crwn_hgt
real :: fire_ext_grnd
real :: fire_ext_crwn
integer :: fire_fuel_read
integer :: fire_fuel_cat
integer :: fire_print_msg
integer :: fire_print_file
integer :: fire_fuel_left_method
integer :: fire_fuel_left_irl
integer :: fire_fuel_left_jrl
real :: fire_atm_feedback
real :: fire_back_weight
integer :: fire_grows_only
integer :: fire_upwinding
integer :: fire_upwind_split
real :: fire_viscosity
real :: fire_lfn_ext_up
integer :: fire_test_steps
integer :: fire_topo_from_atm
integer    :: last_item_in_struct
real :: ash_height
real :: ash_mass
real :: tr_height
real :: tr_mass
!ENDOFREGISTRYGENERATEDINCLUDE
   END TYPE grid_config_rec_type
!
!ENDOFREGISTRYGENERATEDINCLUDE
!TBH:  HACK FOR SMS -- put NAMELIST definition on one line
NAMELIST /chemwrf/ emi_inname,fireemi_inname,emi_outname,fireemi_outname,input_chem_inname,  &
                input_chem_outname,frames_per_emissfile,frames_per_fireemissfile,         &
                io_style_emissions,io_form_emissions,bioemdt,photdt,chemdt,ne_area,kemit, &
                nmegan,kfuture,errosion_dim,chem_conv_tr,chem_opt,gaschem_onoff,          &
                aerchem_onoff,wetscav_onoff,cldchem_onoff,vertmix_onoff,chem_in_opt,      &
                phot_opt,drydep_opt,emiss_opt,dust_opt,dmsemis_opt,seas_opt,bio_emiss_opt,&
                biomass_burn_opt,plumerisefire_frq,emiss_inpt_opt,gas_bc_opt,gas_ic_opt,  &
                aer_bc_opt,aer_ic_opt,have_bcs_chem,aer_ra_feedback,aer_op_opt,           &
                ash_height,ash_mass, tr_height, tr_mass
NAMELIST /wrfphysics/mp_physics,gsfcgce_hail,gsfcgce_2ice,progn,ra_lw_physics,ra_sw_physics, &
                     naer,sf_sfclay_physics,sf_surface_physics,bl_pbl_physics,sf_urban_physics, &
                     cu_physics,num_urban_layers,cugd_avedx,imomentum,          &
                     clos_choice,num_land_cat,num_soil_cat,mp_zero_out,mp_zero_out_thresh,      &
                     seaice_threshold,cu_rad_feedback,slope_rad,topo_shading,topo_shading

        type (grid_config_rec_type) config_flags

END MODULE module_initial_chem_namelists
