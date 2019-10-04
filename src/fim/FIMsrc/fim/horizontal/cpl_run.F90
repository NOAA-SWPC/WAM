module module_fim_cpl_run

IMPLICIT NONE

contains

!*********************************************************************
        subroutine cpl_run(its, dyn_to_phy)
!       "Run" method for the FIM DYN-PHY coupler component.  
!       Argument its is time step count.  
!       Argument dyn_to_phy controls the direction of coupling:  
!         dyn_to_phy == .true.    Couple from DYN to PHY
!         dyn_to_phy == .false.   Couple from PHY to DYN
!       T. Henderson            February, 2009  - Moved code here from physics()
!       R. Bleck                July, 2010      - fixed layer pressure formula
!*********************************************************************

!SMS$IGNORE BEGIN
!TBH:  NOTE removal of "only" clause.  This was forced upon us by *bug* in the 
!TBH:  ifort 11.1 compiler on njet.  Restore the "only" clause when the broken 
!TBH:  compiler is fixed!  
USE gfs_physics_internal_state_mod, only: gfs_physics_internal_state, gis_phy
!SMS$IGNORE END

use module_control  ,only: nts,CallPhysics,itsStart
use module_variables,only: us3d,vs3d,pr3d,tr3d,ws3d
use module_sfc_variables

use module_outtime_cpl,only: telapsed=>tcpl

implicit none

!  Declare dummy arguments
integer, intent(in) :: its
logical, intent(in) :: dyn_to_phy

!  Declare local variables:

real*8 :: t0

call StartTimer(t0)

  !...........................................................
  ! Couple components unless this is the last (nts+1) 
  ! iteration (in which DYN just finishes).  
  ! This complexity is required for the NCEP ESMF approach 
  ! in which single-phase DYN and PHY components alternate 
  ! execution during each time step.  
  !
if (its < itsStart+nts ) then

  !TODO:  Eliminate duplication by encapsulating this logic
  if(mod(its,CallPhysics)==0.or.its==1) then ! Do physics

!sms$compare_var(st3d   , 'begin cpl_run')
!sms$compare_var(sm3d   , 'begin cpl_run')
!sms$compare_var(rn2d   , 'begin cpl_run')
!sms$compare_var(rc2d   , 'begin cpl_run')
!sms$compare_var(ts2d   , 'begin cpl_run')
!sms$compare_var(us2d   , 'begin cpl_run')
!sms$compare_var(hf2d   , 'begin cpl_run')
!sms$compare_var(sw2d   , 'begin cpl_run')
!sms$compare_var(slmsk2d, 'begin cpl_run')

    if (dyn_to_phy) then
      ! Subroutine cpl_dyn_to_phy() converts FIM values to GFS 
      ! values.  
      ! All arrays passed via the ESMF coupler are passed as 
      ! arguments, allowing this subroutine to be called from
      ! the ESMF coupler too.  
      call cpl_dyn_to_phy(its,                                    &
! IN args
             pr3d, us3d, vs3d, ws3d, tr3d,                        &
! OUT args
             gis_phy%ps, gis_phy%dp, gis_phy%p,                   &
             gis_phy%u , gis_phy%v , gis_phy%dpdt,                &
             gis_phy%q , gis_phy%oz, gis_phy%cld,                 &
             gis_phy%t )
    else
      ! Subroutine cpl_phy_to_dyn() converts GFS values to FIM 
      ! values.  
      ! All arrays passed via the ESMF coupler are passed as 
      ! arguments, allowing this subroutine to be called from
      ! the ESMF coupler too.  
      call cpl_phy_to_dyn(its,                                    &
! IN args
             gis_phy%p, gis_phy%u  , gis_phy%v,                   &
             gis_phy%q, gis_phy%cld, gis_phy%t,                   &
             ! these GFS PHY fields are passed to FIM DYN for 
             ! output and diagnostics only.  
             gis_phy%flx_fld%GESHEM,    gis_phy%flx_fld%RAINC,    &
             gis_phy%sfc_fld%TSEA,      gis_phy%sfc_fld%UUSTAR,   &
             gis_phy%flx_fld%HFLX,      gis_phy%flx_fld%EVAP,     &
             gis_phy%sfc_fld%SHELEG,    gis_phy%sfc_fld%CANOPY,   &
             gis_phy%sfc_fld%HICE,      gis_phy%sfc_fld%FICE,     &
             gis_phy%sfc_fld%STC,       gis_phy%sfc_fld%SMC,      &
             gis_phy%flx_fld%SFCDSW,    gis_phy%flx_fld%SFCDLW,   &
             gis_phy%sfc_fld%T2M,       gis_phy%sfc_fld%Q2M,      &
             gis_phy%sfc_fld%SLMSK,     gis_phy%HPRIME,           &
             gis_phy%FLUXR,                                       &
! OUT args
             us3d, vs3d, tr3d, rn2d, rc2d, ts2d, us2d, hf2d, qf2d,&
             sheleg2d, canopy2d, hice2d, fice2d, st3d, sm3d,      &
             sw2d, lw2d, t2m2d, q2m2d, slmsk2d, hprm2d, flxlwtoa2d  )
    endif

!sms$compare_var(st3d   , 'end cpl_run')
!sms$compare_var(sm3d   , 'end cpl_run')
!sms$compare_var(rn2d   , 'end cpl_run')
!sms$compare_var(rc2d   , 'end cpl_run')
!sms$compare_var(ts2d   , 'end cpl_run')
!sms$compare_var(us2d   , 'end cpl_run')
!sms$compare_var(hf2d   , 'end cpl_run')
!sms$compare_var(sw2d   , 'end cpl_run')
!sms$compare_var(slmsk2d, 'end cpl_run')

  endif ! CallPhysics

endif

call IncrementTimer(t0,telapsed)

return
end subroutine cpl_run



! Couple from DYN->PHY.  
subroutine cpl_dyn_to_phy(its,         &
             pr3d,us3d,vs3d,ws3d,tr3d, &
             gfs_ps, gfs_dp, gfs_p,    &
             gfs_u , gfs_v , gfs_dpdt, &
             gfs_q , gfs_oz, gfs_cld,  &
             gfs_t )

!SMS$IGNORE BEGIN
! TODO:  Pass elements of gfs_physics_internal_state_mod via argument 
! TODO:  list and remove use of gfs_physics_internal_state_mod.  
!TBH:  NOTE removal of "only" clause.  This was forced upon us by *bug* in the 
!TBH:  ifort 11.1 compiler on njet.  Restore the "only" clause when the broken 
!TBH:  compiler is fixed!  
!USE gfs_physics_internal_state_mod, only: gfs_physics_internal_state, gis_phy
USE gfs_physics_internal_state_mod
!SMS$IGNORE END

use module_constants,only: cp, rd, p1000, lat, lon, qvmin
use module_control  ,only: nvl,nip,CallRadiation,itsStart
use module_sfc_variables, only: zorl2d,srflag2d
USE MACHINE         ,only: kind_evod

implicit none

  integer, intent(in) :: its
!SMS$DISTRIBUTE (dh,2) BEGIN
  real                 , intent(in ) :: pr3d(:,:)
  real                 , intent(in ) :: us3d(:,:)
  real                 , intent(in ) :: vs3d(:,:)
  real                 , intent(in ) :: ws3d(:,:)
  real                 , intent(in ) :: tr3d(:,:,:)
!SMS$DISTRIBUTE END
!SMS$DISTRIBUTE (dh,1) BEGIN
  real(kind=kind_evod) , intent(out) :: gfs_ps(:)
  real(kind=kind_evod) , intent(out) :: gfs_dp(:,:)
  real(kind=kind_evod) , intent(out) :: gfs_p(:,:)
  real(kind=kind_evod) , intent(out) :: gfs_u(:,:)
  real(kind=kind_evod) , intent(out) :: gfs_v(:,:)
  real(kind=kind_evod) , intent(out) :: gfs_dpdt(:,:)
  real(kind=kind_evod) , intent(out) :: gfs_q(:,:)
  real(kind=kind_evod) , intent(out) :: gfs_oz(:,:)
  real(kind=kind_evod) , intent(out) :: gfs_cld(:,:)
  real(kind=kind_evod) , intent(out) :: gfs_t(:,:)
!SMS$DISTRIBUTE END

!  Declare local variables:
integer :: ipn,ivl
real    :: rocp1,rocpr
!SMS$DISTRIBUTE (dh,nip) BEGIN
real :: tr_2(nvl,nip)
real :: theta_nv(nvl,nip)
!SMS$DISTRIBUTE END

      rocp1=rd/cp+1.
      rocpr=cp/rd

!SMS$PARALLEL (dh,ipn) BEGIN
      do ipn=1,nip

        do ivl=1,nvl
          tr_2(ivl,ipn) = max(qvmin, tr3d(ivl,ipn,2))
        enddo

!NOTE:  ZORL is held constant in FIM-GFS coupling.  Bao confirms that this 
!NOTE:  is what we want.  
        gis_phy%sfc_fld%ZORL(ipn,1)   = zorl2d(ipn)
        !NOTE:  This logic replicates Bao's original logic in which 
        !NOTE:  initial values of TRPCP was read from a file in phy_init() but 
        !NOTE:  overwritten from GESHEM after the first call to phy_run().  
        !NOTE:  Bao has checked the original logic and verified that it 
        !NOTE:  behaved as he intended.  
        if (its > 1) then
          gis_phy%sfc_fld%TPRCP(ipn,1)  = max(0.0d0, gis_phy%flx_fld%GESHEM(ipn,1))
        endif
!NOTE:  SRFLAG is held constant in FIM-GFS coupling.  Bao says "OK".  
        gis_phy%sfc_fld%SRFLAG(ipn,1) = srflag2d(ipn)
        ! Bao confirms that overwrite with TG3 is intentional here.  
        gis_phy%flx_fld%TMPMIN(ipn,1) = gis_phy%sfc_fld%TG3(ipn,1)
        gis_phy%flx_fld%TMPMAX(ipn,1) = gis_phy%sfc_fld%TG3(ipn,1)

        gfs_ps(ipn)         =  0.001*pr3d(1,ipn)
        do ivl=1,nvl
          gfs_dp(ipn,ivl)   = (0.001*pr3d(ivl,ipn))-(0.001*pr3d(ivl+1,ipn))
!!!       gfs_p(ipn,ivl)    = 0.5*(pr3d(ivl,ipn)+pr3d(ivl+1,ipn))
! get energetically consistent mid-lyr prs from partial[p^(kap+1)]/partial[p]
          gfs_p(ipn,ivl)    =					&
            ((pr3d(ivl,ipn)**rocp1-pr3d(ivl+1,ipn)**rocp1)/	&
            ((pr3d(ivl,ipn)       -pr3d(ivl+1,ipn)       )*rocp1))**rocpr
          gfs_u(ipn,ivl)    = us3d(ivl,ipn)
          gfs_v(ipn,ivl)    = vs3d(ivl,ipn)
          gfs_dpdt(ipn,ivl) = 0.001*ws3d(ivl,ipn)
          gfs_q(ipn,ivl)    = tr_2(ivl,ipn)
          gfs_oz(ipn,ivl)   = tr3d(ivl,ipn,4)
          gfs_cld(ipn,ivl ) = tr3d(ivl,ipn,3)
          theta_nv(ivl,ipn) = tr3d(ivl,ipn,1)                           &
              /(1.+0.6078*max(qvmin,tr_2(ivl,ipn)))
          gfs_t(ipn,ivl) = theta_nv(ivl,ipn)*(gfs_p(ipn,ivl)/p1000)**(rd/cp)
        enddo
      enddo

!SMS$PARALLEL END

  return
end subroutine cpl_dyn_to_phy



! Couple from PHY->DYN.  
subroutine cpl_phy_to_dyn(its,                     &
! IN args
             gfs_p, gfs_u  , gfs_v,                &
             gfs_q, gfs_cld, gfs_t,                &
             ! these GFS PHY fields are passed to FIM DYN for 
             ! output and diagnostics only.  
             gfs_geshem, gfs_rainc,                &
             gfs_tsea,   gfs_uustar,               &
             gfs_hflx,   gfs_evap,                 &
             gfs_sheleg, gfs_canopy,               &
             gfs_hice,   gfs_fice,                 &
             gfs_stc,    gfs_smc,                  &
             gfs_sfcdsw, gfs_sfcdlw,               &
             gfs_t2m,    gfs_q2m,                  &
             gfs_slmsk,  gfs_hprime,               &
             gfs_fluxr,                            &
! OUT args
             us3d,  vs3d,    tr3d,                 &
             rn2d, rc2d, ts2d, us2d, hf2d, qf2d,   &
             sheleg2d, canopy2d, hice2d, fice2d,   &
             st3d, sm3d, sw2d, lw2d, t2m2d, q2m2d, &
             slmsk2d, hprm2d, flxlwtoa2d  )

use module_constants,only: cp, rd, p1000, qvmin, qwmin
use module_control  ,only: nts,nvl,nip
USE MACHINE         ,only: kind_evod,kind_phys,kind_rad

implicit none

  integer, intent(in) :: its
!SMS$DISTRIBUTE (dh,1) BEGIN
  real(kind=kind_evod) , intent(in   ) :: gfs_p(:,:)
  real(kind=kind_evod) , intent(in   ) :: gfs_u(:,:)
  real(kind=kind_evod) , intent(in   ) :: gfs_v(:,:)
  real(kind=kind_evod) , intent(in   ) :: gfs_q(:,:)
  real(kind=kind_evod) , intent(in   ) :: gfs_cld(:,:)
  real(kind=kind_evod) , intent(in   ) :: gfs_t(:,:)
  real(kind=kind_phys) , intent(in   ) :: gfs_geshem(:,:)
  real(kind=kind_phys) , intent(in   ) :: gfs_rainc(:,:)
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
  real(kind=kind_phys) , intent(in   ) :: gfs_stc   (:,:,:)
  real(kind=kind_phys) , intent(in   ) :: gfs_smc   (:,:,:)
  real(kind=kind_rad)  , intent(in   ) :: gfs_hprime(:,:,:)
  real(kind=kind_rad)  , intent(in   ) :: gfs_fluxr(:,:,:)
!SMS$DISTRIBUTE END
!SMS$DISTRIBUTE (dh,2) BEGIN
  real                 , intent(  out) :: us3d(:,:)
  real                 , intent(  out) :: vs3d(:,:)
  ! TBH:  tr3d must be inout instead of out because tr3d(:,:,4) is 
  ! TBH:  never set.  
  real                 , intent(inout) :: tr3d(:,:,:)
!SMS$DISTRIBUTE END
!SMS$DISTRIBUTE (dh,1) BEGIN
  real                 , intent(inout) :: rn2d(:)
  real                 , intent(inout) :: rc2d(:)
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
  real                 , intent(  out) :: flxlwtoa2d(:)
!SMS$DISTRIBUTE END
!SMS$DISTRIBUTE (dh,2) BEGIN
  real                 , intent(  out) :: st3d(:,:)
  real                 , intent(  out) :: sm3d(:,:)
  real                 , intent(  out) :: hprm2d(:,:)
!SMS$DISTRIBUTE END

!  Declare local variables:
integer :: ipn,ivl
real (kind=kind_phys) :: rn2dten,rc2dten

!SMS$PARALLEL (dh,ipn) BEGIN

!----------------------------------------------------------------------
!   Move all output into correct FIM arrays
!----------------------------------------------------------------------
      do ipn=1,nip
        do ivl=1,nvl
          us3d(ivl,ipn) = gfs_u(ipn,ivl)
          vs3d(ivl,ipn) = gfs_v(ipn,ivl)
          ! Replace values for prognostic variables
          tr3d(ivl,ipn,1) = gfs_t(ipn,ivl)*(p1000/gfs_p(ipn,ivl))**(rd/cp)*(1.+0.6078*max(REAL(qvmin,kind_evod),gfs_q(ipn,ivl)))
          tr3d(ivl,ipn,2) = max(REAL(qvmin,kind_evod),gfs_q(ipn,ivl))
          tr3d(ivl,ipn,3) = max(REAL(qwmin,kind_evod),gfs_cld(ipn,ivl))
        enddo
!TBH:  Remaining fields are needed only for FIM diagnostics.  
        rn2dten       = max(0.0_kind_phys, gfs_GESHEM(ipn,1))
        rc2dten       = max(0.0_kind_phys, gfs_RAINC(ipn,1))
!TBH:  I assume here that NEMS allows the modified state to be either 
!TBH:  intent(out) *or* intent(inout) ...  
        rn2d    (ipn) = rn2d(ipn) + rn2dten*1000.
        rc2d    (ipn) = rc2d(ipn) + rc2dten*1000.
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
        hprm2d(:,ipn) = gfs_HPRIME(:,ipn,1)
        flxlwtoa2d(ipn) = gfs_FLUXR(1,ipn,1)
      enddo !ipn

!SMS$PARALLEL END

end subroutine cpl_phy_to_dyn

end module module_fim_cpl_run
