module module_fim_cpl_init

implicit none

contains

!*********************************************************************
!       Initialize the FIM DYN-PHY coupler component.  
!       T. Henderson            February, 2009
!*********************************************************************
subroutine cpl_init

use module_constants,only: lat, lon
use module_control  ,only: nip
use module_sfc_variables,only: ts2d,us2d,hf2d,qf2d,                  &
                               sheleg2d, canopy2d, hice2d, fice2d,   &
                               st3d, sm3d, sw2d, lw2d, t2m2d, q2m2d, &
                               slmsk2d
!SMS$IGNORE BEGIN
USE gfs_physics_internal_state_mod, only: gfs_physics_internal_state, gis_phy
!SMS$IGNORE END

  implicit none

  ! Local variables
  real*8 :: t0,t1=0.0d0

  call StartTimer(t0)

!TODO:  For the moment, pass lon,lat from DYN to PHY here.  Later we will 
!TODO:  hopefully change to passing these via the ESMF_Grid and 
!TODO:  cpl_init_dyn_to_phy() can be eliminated, improving interoperability.  
!TODO:  Separate into two-phase CPL init as is done in run.F90?  See 
!TODO:  module_DYN_PHY_CPL_COMP.F90 for details.

  call cpl_init_dyn_to_phy(nip,      &
! IN args
                  lon,          lat, &
! OUT args
         gis_phy%xlon, gis_phy%xlat)

!TODO:  At the moment, NEMS does not pass anything between DYN and PHY 
!TODO:  during INIT.  We would like to write initial values of PHY variables 
!TODO:  to disk in our 0-hour history output.  For now pass from PHY to DYN 
!TODO:  here.  If NEMS comes up with a solution (maybe via the new NEMS I/O) 
!TODO:  that allows PHY to write directly via the I/O component(s), then 
!TODO:  replace cpl_init_phy_to_dyn() with the new method.  If not then it may 
!TODO:  be necessary to extend NEMS to allow transfer of fields from DYN to PHY 
!TODO:  in CPL_INIT after both DYN and PHY have finished their init phases.  

  call cpl_init_phy_to_dyn(nip,                               &
! IN args
         ! these GFS PHY fields are passed to FIM DYN for 
         ! output and diagnostics only.  
         gis_phy%sfc_fld%TSEA,      gis_phy%sfc_fld%UUSTAR,   &
         gis_phy%flx_fld%HFLX,      gis_phy%flx_fld%EVAP,     &
         gis_phy%sfc_fld%SHELEG,    gis_phy%sfc_fld%CANOPY,   &
         gis_phy%sfc_fld%HICE,      gis_phy%sfc_fld%FICE,     &
         gis_phy%sfc_fld%STC,       gis_phy%sfc_fld%SMC,      &
         gis_phy%flx_fld%SFCDSW,    gis_phy%flx_fld%SFCDLW,   &
         gis_phy%sfc_fld%T2M,       gis_phy%sfc_fld%Q2M,      &
         gis_phy%sfc_fld%SLMSK,                               &
! OUT args
         ts2d, us2d, hf2d, qf2d,                              &
         sheleg2d, canopy2d, hice2d, fice2d, st3d, sm3d,      &
         sw2d, lw2d, t2m2d, q2m2d, slmsk2d  )

  call IncrementTimer(t0,t1)

  print"(' COUPLER INIT time:',F10.0)",t1

  return
end subroutine cpl_init


!*********************************************************************
!       Couple from DYN->PHY during INIT phase.  
!*********************************************************************
subroutine cpl_init_dyn_to_phy(nip,lon,lat,gfs_lon,gfs_lat)

USE MACHINE         ,only: kind_rad

  implicit none

  ! Arguments
  integer, intent(in) :: nip
!SMS$DISTRIBUTE (dh,1) BEGIN
  real                , intent(in)  :: lon(:)
  real                , intent(in)  :: lat(:)
  real(kind=kind_rad) , intent(out) :: gfs_lon(:,:)
  real(kind=kind_rad) , intent(out) :: gfs_lat(:,:)
!SMS$DISTRIBUTE END

  ! Local variables
  integer :: ipn

!SMS$PARALLEL (dh,ipn) BEGIN
  do ipn=1,nip
    gfs_lon(ipn,1) = lon(ipn)
    gfs_lat(ipn,1) = lat(ipn)
  enddo
!SMS$PARALLEL END

  return
end subroutine cpl_init_dyn_to_phy


!*********************************************************************
!       Couple from PHY->DYN during INIT phase.  
!*********************************************************************
subroutine cpl_init_phy_to_dyn(nip,                &
             ! these GFS PHY fields are passed to FIM DYN for 
             ! output and diagnostics only.  
             gfs_tsea,   gfs_uustar,               &
             gfs_hflx,   gfs_evap,                 &
             gfs_sheleg, gfs_canopy,               &
             gfs_hice,   gfs_fice,                 &
             gfs_stc,    gfs_smc,                  &
             gfs_sfcdsw, gfs_sfcdlw,               &
             gfs_t2m,    gfs_q2m,                  &
             gfs_slmsk,                            &
! OUT args
             ts2d, us2d, hf2d, qf2d,               &
             sheleg2d, canopy2d, hice2d, fice2d,   &
             st3d, sm3d, sw2d, lw2d, t2m2d, q2m2d, &
             slmsk2d  )

USE MACHINE         ,only: kind_phys

  implicit none

  ! Arguments
  integer, intent(in) :: nip
!SMS$DISTRIBUTE (dh,1) BEGIN
  real(kind=kind_phys) , intent(in   ) :: gfs_tsea(:,:)
  real(kind=kind_phys) , intent(in   ) :: gfs_uustar(:,:)
  real(kind=kind_phys) , intent(in   ) :: gfs_hflx(:,:)
  real(kind=kind_phys) , intent(in   ) :: gfs_evap(:,:)
  real(kind=kind_phys) , intent(in   ) :: gfs_sheleg(:,:)
  real(kind=kind_phys) , intent(in   ) :: gfs_canopy(:,:)
  real(kind=kind_phys) , intent(in   ) :: gfs_hice  (:,:)
  real(kind=kind_phys) , intent(in   ) :: gfs_fice  (:,:)
  real(kind=kind_phys) , intent(in   ) :: gfs_sfcdsw(:,:)
  real(kind=kind_phys) , intent(in   ) :: gfs_sfcdlw(:,:)
  real(kind=kind_phys) , intent(in   ) :: gfs_t2m   (:,:)
  real(kind=kind_phys) , intent(in   ) :: gfs_q2m   (:,:)
  real(kind=kind_phys) , intent(in   ) :: gfs_slmsk (:,:)
!SMS$DISTRIBUTE END
!SMS$DISTRIBUTE (dh,2) BEGIN
  real(kind=kind_phys) , intent(in   ) :: gfs_stc (:,:,:)
  real(kind=kind_phys) , intent(in   ) :: gfs_smc (:,:,:)
!SMS$DISTRIBUTE END
!SMS$DISTRIBUTE (dh,1) BEGIN
  real                 , intent(  out) :: ts2d(:)
  real                 , intent(  out) :: us2d(:)
  real                 , intent(  out) :: hf2d(:)
  real                 , intent(  out) :: qf2d(:)
  real                 , intent(  out) :: sheleg2d(:)
  real                 , intent(  out) :: canopy2d(:)
  real                 , intent(  out) :: hice2d(:)
  real                 , intent(  out) :: fice2d(:)
  real                 , intent(  out) :: sw2d(:)
  real                 , intent(  out) :: lw2d(:)
  real                 , intent(  out) :: t2m2d(:)
  real                 , intent(  out) :: q2m2d(:)
  real                 , intent(  out) :: slmsk2d(:)
!SMS$DISTRIBUTE END
!SMS$DISTRIBUTE (dh,2) BEGIN
  real                 , intent(  out) :: st3d(:,:)
  real                 , intent(  out) :: sm3d(:,:)
!SMS$DISTRIBUTE END

  ! Local variables
  integer :: ipn

!SMS$PARALLEL (dh,ipn) BEGIN

!----------------------------------------------------------------------
!   Move all output into correct FIM arrays for 00h output
!----------------------------------------------------------------------
!TODO:  remove duplication with cpl_phy_to_dyn() if practical
  do ipn=1,nip
    ts2d    (ipn) = gfs_TSEA  (ipn,1)
    us2d    (ipn) = gfs_UUSTAR(ipn,1)
    hf2d    (ipn) = gfs_HFLX  (ipn,1)
    qf2d    (ipn) = gfs_EVAP  (ipn,1)
    sheleg2d(ipn) = gfs_SHELEG(ipn,1)
    canopy2d(ipn) = gfs_CANOPY(ipn,1)
    hice2d  (ipn) = gfs_HICE  (ipn,1)
    fice2d  (ipn) = gfs_FICE  (ipn,1)
    st3d  (:,ipn) = gfs_STC (:,ipn,1)
    sm3d  (:,ipn) = gfs_SMC (:,ipn,1)
    sw2d    (ipn) = gfs_SFCDSW(ipn,1)
    lw2d    (ipn) = gfs_SFCDLW(ipn,1)
    t2m2d   (ipn) = gfs_T2M   (ipn,1)
    q2m2d   (ipn) = gfs_Q2M   (ipn,1)
    slmsk2d (ipn) = gfs_SLMSK (ipn,1)
  enddo !ipn

!SMS$PARALLEL END

  return
end subroutine cpl_init_phy_to_dyn

end module module_fim_cpl_init
