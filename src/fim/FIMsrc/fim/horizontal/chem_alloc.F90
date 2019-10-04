!*********************************************************************
module module_chem_alloc
!	This module allocates variables used in chem
!  	Henderson          November      2008
!*********************************************************************
contains


subroutine chem_alloc(chem_opt,aer_ra_feedback)

use module_control,only: nvl,nip
use module_wrf_control,only: num_emis_ant,num_emis_vol,nvl_gocart,num_chem
use module_chem_variables,only: emiss_ab,pm25,p10,rcav,ero1,ero2,ero3,dm0, &
                                emiss_oc,emiss_bc,emiss_sulf,oh_backgd,h2o2_backgd,  &
                                no3_backgd,emiss_abu,plumestuff,aod2d,emiss_ash_mass,&
                                emiss_ash_height,emiss_ash_dt,ashfall,rnav,tr1_tavg, &
                                emiss_tr_height,emiss_tr_mass,emiss_tr_dt,trfall,    &
                                clayfrac,sandfrac,emiss_co2

implicit none

integer, intent(IN) :: chem_opt,aer_ra_feedback

! always allocate these because they are passed 
! through an arglist
!TODO:  Create chem_internal_state and avoid all this complication.  
!TODO:  Avoid allocating *any* arrays that are not used!
print *,'DEBUG chem_alloc():  allocating...'
  allocate(pm25(nvl,nip))
  pm25=0.
  allocate(p10(nvl,nip))
  p10=0.
  allocate(tr1_tavg(nvl,nip))
  tr1_tavg=0.
  allocate(oh_backgd(nvl_gocart,nip))
  oh_backgd=0.
  allocate( h2o2_backgd(nvl_gocart,nip) )
  allocate( no3_backgd(nvl_gocart,nip) )
  allocate( rcav(nip) )
  rcav = 0.
  allocate( rnav(nip) )
  rnav = 0.
  allocate( ero1(nip) )
  allocate( ero2(nip) )
  allocate( ero3(nip) )
  allocate( clayfrac(nip) )
  allocate( sandfrac(nip) )
  allocate( ashfall(nip) )
  allocate( aod2d(nip) )
  aod2d = 0.
  allocate( plumestuff(nip,8) )
  plumestuff = 0.
  allocate( emiss_ab(nip,num_emis_ant) )
  emiss_ab = 0.
  allocate( emiss_abu(nip,num_emis_ant) )
  emiss_abu = 0.
  allocate( emiss_ash_mass(nip) )
  emiss_ash_mass = 0.
  allocate( emiss_ash_height(nip) )
  emiss_ash_height = 0.
  allocate( emiss_ash_dt(nip) )
  emiss_ash_dt = 0.
  allocate( emiss_co2(nip) )
  emiss_co2 = 0.

!if(chem_opt == 500)then
    allocate( emiss_tr_mass(nip) )
    emiss_tr_mass = 0.
    allocate( emiss_tr_height(nip) )
    emiss_tr_height = 0.
    allocate( emiss_tr_dt(nip) )
    emiss_tr_dt = 0.
    ALLOCATE( trfall( nip, num_chem ) )
    trfall = 0.
! endif

  allocate( emiss_oc(nip) )
  emiss_oc = 0.
  allocate( emiss_bc(nip) )
  emiss_bc = 0.
  allocate( emiss_sulf(nip) )
  emiss_sulf = 0.
  allocate(  dm0(nip) )
  dm0 = 0.

return
end subroutine chem_alloc


!TODO:  combine with chem_alloc if practical
subroutine chem_alloc2(chem_opt,aer_ra_feedback,bio_emiss_opt,biomass_burn_opt,kemit)

  USE module_wrf_control, only: ims,ime,jms,jme,kms,kme,nbands,  &
                          num_emis_vol,num_moist,num_chem,num_emis_ant, &
                          num_ext_coef,num_bscat_coef,num_asym_par
  ! yes, use ALL of it
  use module_chemvars

  implicit none

  integer, intent(IN) :: chem_opt,aer_ra_feedback,kemit, &
                         bio_emiss_opt,biomass_burn_opt
  ALLOCATE( chem( ims:ime, kms:kme, jms:jme, num_chem ) )
  ALLOCATE( e_bio( ims:ime, jms:jme, ne_area ) )
  ALLOCATE( emis_ant( ims:ime, 1:kemit, jms:jme,num_emis_ant) )
emis_ant = 0.
  ALLOCATE( emis_vol( ims:ime, kms:kme, jms:jme,num_emis_vol) )
emis_vol=0.
  ALLOCATE( relhum( ims:ime, kms:kme, jms:jme ) )
  ALLOCATE( dms_0( ims:ime, jms:jme) )
  ALLOCATE( erod( ims:ime, jms:jme,3) )
  ALLOCATE( emis_dust( ims:ime, 1, jms:jme,num_emis_dust) )
emis_dust = 0.
  ALLOCATE( emis_seas( ims:ime, 1, jms:jme,num_emis_seas) )
emis_seas = 0.
  ALLOCATE( backg_oh( ims:ime, kms:kme, jms:jme ) )
  ALLOCATE( backg_h2o2( ims:ime, kms:kme, jms:jme ) )
  ALLOCATE( backg_no3( ims:ime, kms:kme, jms:jme ) )
  ALLOCATE( ebu_in_no( ims:ime, jms:jme ) )
  ebu_in_no=0.
  ALLOCATE( ebu_in_co( ims:ime, jms:jme ) )
  ebu_in_co=0.
  ALLOCATE( ebu_in_co2( ims:ime, jms:jme ) )
  ebu_in_co2=0.
  ALLOCATE( ebu_in_eth( ims:ime, jms:jme ) )
  ebu_in_eth=0.
  ALLOCATE( ebu_in_hc3( ims:ime, jms:jme ) )
  ebu_in_hc3=0.
  ALLOCATE( ebu_in_hc5( ims:ime, jms:jme ) )
  ebu_in_hc5=0.
  ALLOCATE( ebu_in_hc8( ims:ime, jms:jme ) )
  ebu_in_hc8=0.
  ALLOCATE( ebu_in_ete( ims:ime, jms:jme ) )
  ebu_in_ete=0.
  ALLOCATE( ebu_in_olt( ims:ime, jms:jme ) )
  ebu_in_olt=0.
  ALLOCATE( ebu_in_oli( ims:ime, jms:jme ) )
  ebu_in_oli=0.
  ALLOCATE( ebu_in_pm25( ims:ime, jms:jme ) )
  ebu_in_pm25=0.
  ALLOCATE( ebu_in_pm10( ims:ime, jms:jme ) )
  ebu_in_pm10=0.
  ALLOCATE( ebu_in_oc( ims:ime, jms:jme ) )
  ebu_in_oc=0.
  ALLOCATE( ebu_in_bc( ims:ime, jms:jme ) )
  ebu_in_bc=0.
  ALLOCATE( ebu_in_so2( ims:ime, jms:jme ) )
  ebu_in_so2=0.
  ALLOCATE( ebu_in_sulf( ims:ime, jms:jme ) )
  ebu_in_sulf=0.
  ALLOCATE( ebu_in_dien( ims:ime, jms:jme ) )
  ebu_in_dien=0.
  ALLOCATE( ebu_in_iso( ims:ime, jms:jme ) )
  ebu_in_iso=0.
  ALLOCATE( ebu_in_api( ims:ime, jms:jme ) )
  ebu_in_api=0.
  ALLOCATE( ebu_in_lim( ims:ime, jms:jme ) )
  ebu_in_lim=0.
  ALLOCATE( ebu_in_tol( ims:ime, jms:jme ) )
  ebu_in_tol=0.
  ALLOCATE( ebu_in_xyl( ims:ime, jms:jme ) )
  ebu_in_xyl=0.
  ALLOCATE( ebu_in_csl( ims:ime, jms:jme ) )
  ebu_in_csl=0.
  ALLOCATE( ebu_in_hcho( ims:ime, jms:jme ) )
  ebu_in_hcho=0.
  ALLOCATE( ebu_in_ald( ims:ime, jms:jme ) )
  ebu_in_ald=0.
  ALLOCATE( ebu_in_ket( ims:ime, jms:jme ) )
  ebu_in_ket=0.
  ALLOCATE( ebu_in_macr( ims:ime, jms:jme ) )
  ebu_in_macr=0.
  ALLOCATE( ebu_in_ora1( ims:ime, jms:jme ) )
  ebu_in_ora1=0.
  ALLOCATE( ebu_in_ora2( ims:ime, jms:jme ) )
  ebu_in_ora2=0.
  ALLOCATE( ebu_no( ims:ime, kms:kme, jms:jme ) )
  ALLOCATE( ebu_co( ims:ime, kms:kme, jms:jme ) )
  ALLOCATE( ebu_co2( ims:ime, kms:kme, jms:jme ) )
  ALLOCATE( ebu_eth( ims:ime, kms:kme, jms:jme ) )
  ALLOCATE( ebu_hc3( ims:ime, kms:kme, jms:jme ) )
  ALLOCATE( ebu_hc5( ims:ime, kms:kme, jms:jme ) )
  ALLOCATE( ebu_hc8( ims:ime, kms:kme, jms:jme ) )
  ALLOCATE( ebu_ete( ims:ime, kms:kme, jms:jme ) )
  ALLOCATE( ebu_olt( ims:ime, kms:kme, jms:jme ) )
  ALLOCATE( ebu_oli( ims:ime, kms:kme, jms:jme ) )
  ALLOCATE( ebu_pm25( ims:ime, kms:kme, jms:jme ) )
  ALLOCATE( ebu_pm10( ims:ime, kms:kme, jms:jme ) )
  ALLOCATE( ebu_oc( ims:ime, kms:kme, jms:jme ) )
  ALLOCATE( ebu_bc( ims:ime, kms:kme, jms:jme ) )
  ALLOCATE( ebu_so2( ims:ime, kms:kme, jms:jme ) )
  ALLOCATE( ebu_sulf( ims:ime, kms:kme, jms:jme ) )
  ALLOCATE( ebu_dien( ims:ime, kms:kme, jms:jme ) )
  ALLOCATE( ebu_iso( ims:ime, kms:kme, jms:jme ) )
  ALLOCATE( ebu_api( ims:ime, kms:kme, jms:jme ) )
  ALLOCATE( ebu_lim( ims:ime, kms:kme, jms:jme ) )
  ALLOCATE( ebu_tol( ims:ime, kms:kme, jms:jme ) )
  ALLOCATE( ebu_xyl( ims:ime, kms:kme, jms:jme ) )
  ALLOCATE( ebu_csl( ims:ime, kms:kme, jms:jme ) )
  ALLOCATE( ebu_hcho( ims:ime, kms:kme, jms:jme ) )
  ALLOCATE( ebu_ald( ims:ime, kms:kme, jms:jme ) )
  ALLOCATE( ebu_ket( ims:ime, kms:kme, jms:jme ) )
  ALLOCATE( ebu_macr( ims:ime, kms:kme, jms:jme ) )
  ALLOCATE( ebu_ora1( ims:ime, kms:kme, jms:jme ) )
  ALLOCATE( ebu_ora2( ims:ime, kms:kme, jms:jme ) )
  ALLOCATE( mean_fct_agtf( ims:ime,  jms:jme ) )
  ALLOCATE( mean_fct_agef( ims:ime,  jms:jme ) )
  ALLOCATE( mean_fct_agsv( ims:ime,  jms:jme ) )
  ALLOCATE( mean_fct_aggr( ims:ime,  jms:jme ) )
  ALLOCATE( firesize_agtf( ims:ime,  jms:jme ) )
  ALLOCATE( firesize_agef( ims:ime,  jms:jme ) )
  ALLOCATE( firesize_agsv( ims:ime,  jms:jme ) )
  ALLOCATE( firesize_aggr( ims:ime,  jms:jme ) )
  ALLOCATE( ash_fall( ims:ime,  jms:jme ) )
  ash_fall=0.
  ALLOCATE( dust_fall( ims:ime,  jms:jme ) )
  ALLOCATE( pm2_5_dry( ims:ime , kms:kme , jms:jme ) )
  ALLOCATE( pm2_5_dry_ec( ims:ime , kms:kme , jms:jme ) )
  ALLOCATE( pm10( ims:ime , kms:kme , jms:jme ) )
  ALLOCATE( tcosz( ims:ime , jms:jme ) )
  ALLOCATE( ttday( ims:ime , jms:jme ) )

  ALLOCATE( sebio_iso( ims:ime , jms:jme ) )
  ALLOCATE( sebio_oli( ims:ime , jms:jme ) )
  ALLOCATE( sebio_api( ims:ime , jms:jme ) )
  ALLOCATE( sebio_lim( ims:ime , jms:jme ) )
  ALLOCATE( sebio_xyl( ims:ime , jms:jme ) )
  ALLOCATE( sebio_hc3( ims:ime , jms:jme ) )
  ALLOCATE( sebio_ete( ims:ime , jms:jme ) )
  ALLOCATE( sebio_olt( ims:ime , jms:jme ) )
  ALLOCATE( sebio_ket( ims:ime , jms:jme ) )
  ALLOCATE( sebio_ald( ims:ime , jms:jme ) )
  ALLOCATE( sebio_hcho( ims:ime , jms:jme ) )
  ALLOCATE( sebio_eth( ims:ime , jms:jme ) )
  ALLOCATE( sebio_ora2( ims:ime , jms:jme ) )
  ALLOCATE( sebio_co( ims:ime , jms:jme ) )
  ALLOCATE( sebio_nr( ims:ime , jms:jme ) )
  ALLOCATE( noag_grow( ims:ime , jms:jme ) )
  ALLOCATE( noag_nongrow( ims:ime , jms:jme ) )
  ALLOCATE( nononag( ims:ime , jms:jme ) )
  ALLOCATE( slai( ims:ime , jms:jme ) )
  ALLOCATE( ebio_iso( ims:ime , jms:jme ) )
  ALLOCATE( ebio_oli( ims:ime , jms:jme ) )
  ALLOCATE( ebio_api( ims:ime , jms:jme ) )
  ALLOCATE( ebio_lim( ims:ime , jms:jme ) )
  ALLOCATE( ebio_xyl( ims:ime , jms:jme ) )
  ALLOCATE( ebio_hc3( ims:ime , jms:jme ) )
  ALLOCATE( ebio_ete( ims:ime , jms:jme ) )
  ALLOCATE( ebio_olt( ims:ime , jms:jme ) )
  ALLOCATE( ebio_ket( ims:ime , jms:jme ) )
  ALLOCATE( ebio_ald( ims:ime , jms:jme ) )
  ALLOCATE( ebio_hcho( ims:ime , jms:jme ) )
  ALLOCATE( ebio_eth( ims:ime , jms:jme ) )
  ALLOCATE( ebio_ora2( ims:ime , jms:jme ) )
  ALLOCATE( ebio_co( ims:ime , jms:jme ) )
  ALLOCATE( ebio_nr( ims:ime , jms:jme ) )
  ALLOCATE( ebio_no( ims:ime , jms:jme ) )
  if(bio_emiss_opt == 3)then
  ALLOCATE( EFmegan(ims:ime, jms:jme , nmegan) )

  ALLOCATE( msebio_isop(ims:ime, jms:jme ) )
  ALLOCATE( pftp_bt(ims:ime, jms:jme ) )
  ALLOCATE( pftp_nt(ims:ime, jms:jme ) )
  ALLOCATE( pftp_sb(ims:ime, jms:jme ) )
  ALLOCATE( pftp_hb(ims:ime, jms:jme ) )

  ALLOCATE( mlai(ims:ime, jms:jme, 12 ) )
  ALLOCATE( mtsa(ims:ime, jms:jme, 12 ) )
  ALLOCATE( mswdown(ims:ime, jms:jme, 12 ) )

  ALLOCATE( mebio_isop(ims:ime, jms:jme ) )
  ALLOCATE( mebio_apin(ims:ime, jms:jme ) )
  ALLOCATE( mebio_bpin(ims:ime, jms:jme ) )
  ALLOCATE( mebio_bcar(ims:ime, jms:jme ) )
  ALLOCATE( mebio_acet(ims:ime, jms:jme ) )
  ALLOCATE( mebio_mbo(ims:ime, jms:jme ) )
  ALLOCATE( mebio_no(ims:ime, jms:jme ) )
  endif

  if(chem_opt == 2)then
     ALLOCATE( h2oai(ims:ime, kms:kme, jms:jme ) )
     ALLOCATE( h2oaj(ims:ime, kms:kme, jms:jme ) )
  endif
  if(aer_ra_feedback == 1)then
     ALLOCATE( extt(ims:ime, kms:kme, jms:jme,nbands) )
     ALLOCATE( ssca(ims:ime, kms:kme, jms:jme,nbands) )
     ALLOCATE( asympar(ims:ime, kms:kme, jms:jme,nbands) )
     ALLOCATE( aod(ims:ime, jms:jme ) )
     ALLOCATE( ext_coeff(ims:ime, kms:kme, jms:jme,1:num_ext_coef ) )
     ALLOCATE( bscat_coeff(ims:ime, kms:kme, jms:jme,1:num_bscat_coef ) )
     ALLOCATE( asym_par(ims:ime, kms:kme, jms:jme,1:num_asym_par ) )
     ALLOCATE( tauaerlw(ims:ime, kms:kme, jms:jme,1:16 ) )
     ALLOCATE( tauaersw(ims:ime, kms:kme, jms:jme,1:4 ) )
     ALLOCATE( gaersw(ims:ime, kms:kme, jms:jme, 1:4 ) )
     ALLOCATE( waersw(ims:ime, kms:kme, jms:jme, 1:4 ) )
     ALLOCATE( bscoefsw(ims:ime, kms:kme, jms:jme, 1:4 ) )
     ALLOCATE( l2aer(ims:ime, kms:kme, jms:jme, 1:4 ) )
     ALLOCATE( l3aer(ims:ime, kms:kme, jms:jme, 1:4 ) )
     ALLOCATE( l4aer(ims:ime, kms:kme, jms:jme, 1:4 ) )
     ALLOCATE( l5aer(ims:ime, kms:kme, jms:jme, 1:4 ) )
     ALLOCATE( l6aer(ims:ime, kms:kme, jms:jme, 1:4 ) )
     ALLOCATE( l7aer(ims:ime, kms:kme, jms:jme, 1:4 ) )
  endif

return
end subroutine chem_alloc2

end module module_chem_alloc
