Module module_set_wrfphys
USE module_initial_chem_namelists
CONTAINS
SUBROUTINE set_wrfphys (mp_physics)
USE module_wrfphysvars
USE module_control, only: ntra,ntrb
USE module_wrf_control, only: num_moist,num_chem
implicit none
!! TBH:  Ignore these so PPP doesn't have to translate them
!!SMS$IGNORE BEGIN
!USE module_initial_chem_namelists
!USE module_data_gocart_dust
!USE module_data_gocart_seas
!!SMS$IGNORE END
integer, intent(in) :: mp_physics
integer :: itest
!!SMS$DISTRIBUTE (dh,2) BEGIN
!real, intent(in) :: tr3d(:,:,:)
!integer :: itest
!!SMS$DISTRIBUTE END
!
! microphysics
!
if(mp_physics.eq.4)then
   if(num_moist.ne.5)then
      write(6,*) ' num_moist is not equal 5 '
      stop
   endif
   itest=ntra+ntrb-num_moist-num_chem+1
   if(itest.ne.ntra)then
      write(6,*) ' ntra ist falsch'
      write(6,*) ' ntra,num_moist,num_chem = ',ntra,num_moist,num_chem
      stop
   endif
   p_qv=1
   f_qv=.true.
   p_qc=2
   f_qc=.true.
   p_qr=3
   f_qr=.true.
   p_qi=4
   f_qi=.true.
   p_qs=5
   f_qs=.true.
else if(mp_physics.eq.2)then
   if(num_moist.ne.6)then
      write(6,*) ' num_moist is not equal 6 '
      stop
   endif
   itest=ntra+ntrb-num_moist-num_chem+1
   if(itest.ne.ntra)then
      write(6,*) ' ntra ist falsch'
      stop
   endif
   p_qv=1
   f_qv=.true.
   p_qc=2
   f_qc=.true.
   p_qr=3
   f_qr=.true.
   p_qi=4
   f_qi=.true.
   p_qs=5
   f_qs=.true.
   p_qg=6
   f_qg=.true.
else if(mp_physics.eq.0)then
   if(num_moist.ne.3)then
      write(6,*) ' num_moist is not equal 3 '
!     stop
   endif
   itest=ntra+ntrb-num_moist-num_chem+3
   if(itest.ne.ntra)then
      write(6,*) ' ntra ist falsch',ntra,num_moist,num_chem
!     stop
   endif
   p_qv=1
   f_qv=.true.
   p_qc=2
   f_qc=.true.
endif
END SUBROUTINE set_wrfphys
       subroutine set_wrfphys_namelist_defaults
implicit none
!STARTOFREGISTRYGENERATEDINCLUDE 'inc/namelist_defaults.inc'
! from V3.1
!
! THIS does not need to be recopied. Defaults are simple to set!
!
mp_physics = 0
gsfcgce_hail = 0
gsfcgce_2ice = 0
progn = 0
ra_lw_physics = 0
ra_sw_physics = 0
radt = 0 
naer = 1e9
sf_sfclay_physics = 0
sf_surface_physics = 0
bl_pbl_physics = 0
sf_urban_physics = 0
bldt = 0
cu_physics = 0
cudt = 0
gsmdt = 0
isfflx = 1
ifsnow = 0
icloud = 1
swrad_scat = 1
surface_input_source = 1
num_urban_layers = 400
num_months = 12
maxiens = 1
maxens = 3
maxens2 = 3
maxens3 = 16
ensdim = 144
cugd_avedx = 1
imomentum = 0
clos_choice = 0
num_land_cat = 24
num_soil_cat = 16
mp_zero_out = 0
mp_zero_out_thresh = 1.e-8
seaice_threshold = 271
sst_update = 0
sst_skin = 0
tmn_update = 0
usemonalb = .false.
rdmaxalb = .true.
rdlai2d = .false. 
co2tf = 1
ra_call_offset = 0
cam_abs_freq_s = 21600.
levsiz = 1
paerlev = 1
cam_abs_dim1 = 1
cam_abs_dim2 = 1
lagday = 1
cu_rad_feedback = .false.
pxlsm_smois_init = 1
omlcall = 0
oml_hml0 = 50
oml_gamma = 0.14
isftcflx = 0
shadlen = 25000.
slope_rad = 0
topo_shading = 0
no_mp_heating = 0
fractional_seaice = 0
bucket_mm = -1.
bucket_j = -1.
grav_settling = 0
END SUBROUTINE set_wrfphys_namelist_defaults
END MODULE module_set_wrfphys
