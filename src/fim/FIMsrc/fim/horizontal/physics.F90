module module_physics
integer :: ipn          ! Index for icos point number
integer :: itsP         ! Public version of its

contains
!*********************************************************************
!     physics
!	Calculates column forcing for global fim
!	12/21/2005 - Alexander E. MacDonald     - original version
!	05/01/2006 - Jian-Wen Bao               - modified for GFS physics
!       04/14/2008 - Stan Benjamin, John Brown  - modifications
!                      for introduction of virtual pot temp for temp prog variable
!                      instead of previous non-virtual pot temp
!       03/03/2009 - Tom Henderson  - split into do_physics_one_step and cpl_run
!*********************************************************************

subroutine physics (its,  &
CallPhysics,CallRadiation )      ! Timestep interval to call physics,radiation

!SMS$IGNORE BEGIN
USE gfs_physics_internal_state_mod, only:gis_phy
!SMS$IGNORE END

use module_control,only:nip,dt,numphr,yyyymmddhhmm,GravityWaveDrag
use module_constants,only: inv_perm

!GG:  added arrays for wrf physics and chemistry
!TODO:  send these via internal state and coupler instead of via use-association
use module_chem_variables, only: sscal, ext_cof, asymp, extlw_cof
use module_wrf_variables,  only: phys3dwrf,exch,pb2d
!SMS$IGNORE BEGIN
use module_initial_chem_namelists, only: cu_physics, mp_physics, chem_opt
!SMS$IGNORE END
use module_do_physics_one_step,only:do_physics_one_step

implicit none

integer,intent (IN   ) :: its         ! model time step count
integer,intent (IN   ) :: CallPhysics,CallRadiation

! Local variables
!----------------------------------------------------------------------
!
!GG:  added switches for cumulus convection and microphysics
logical :: skip_cu_physics, skip_mp_physics, skip_chem

!print *,'DEBUG physics():  SIZE(pb2d) = ',SIZE(pb2d)
!----------------------------------------------------------------------
!TODO:  Eliminate duplication by encapsulating this logic
if (mod (its, CallPhysics) == 0 .or. its == 1) then ! Do physics
!----------------------------------------------------------------------

if (chem_opt.gt.0) then
!sms$compare_var(pb2d,'begin physics')
!sms$compare_var(exch,'begin physics')
!sms$compare_var(phys3dwrf,'begin physics')
endif

  gis_phy%kdt = its
!TODO:  set gis_phy%deltim in phy_init?  
  gis_phy%deltim  = dt
  if(CallPhysics>1) then !Adjust the timestep to be the call physics interval
    if(its==CallPhysics) then
      gis_phy%deltim = (CallPhysics-1)*dt
    elseif(its>CallPhysics) then
      gis_phy%deltim = CallPhysics*dt
    endif
  endif
  gis_phy%phour  = float(its-1)/float(numphr)

!GG:  added switches for cumulus convection and microphysics
  skip_cu_physics = (cu_physics /= 0)
  skip_mp_physics = (mp_physics /= 0)
  skip_chem       = (chem_opt   /= 0)

  call do_physics_one_step(                                               &
                 gis_phy%deltim,  gis_phy%kdt,     gis_phy%phour,         &
                 gis_phy%ps, gis_phy%dp, gis_phy%dpdt,                    &
                 gis_phy%p,  gis_phy%u,  gis_phy%v,                       &
                 gis_phy%t,                                               &
                 gis_phy%q,                                               &
                 gis_phy%oz,                                              &
                 gis_phy%cld,                                             &
                 gis_phy%sfc_fld, gis_phy%flx_fld,                        &
! These are hard-coded in do_physics_one_step() for the moment
!                gis_phy%lats_nodes_r,   gis_phy%global_lats_r,           &
!                gis_phy%lonsperlar,                                      &
                 gis_phy%XLON,    gis_phy%XLAT,    gis_phy%COSZDG,        &
                 gis_phy%HPRIME,  gis_phy%SWH,     gis_phy%HLW,           &
                 gis_phy%FLUXR,   gis_phy%SFALB,                          &
                 gis_phy%SLAG,    gis_phy%SDEC,    gis_phy%CDEC,          &
! Not used yet by FIM
!                gis_phy%OZPLIN,  gis_phy%JINDX1,  gis_phy%JINDX2,        &
!                gis_phy%DDY,                                             &
                 gis_phy%phy_f3d, gis_phy%phy_f2d, gis_phy%NBLCK,         &
! Not used yet by FIM
!                gis_phy%ZHOUR,   gis_phy%N3,      gis_phy%N4,            &
!                gis_phy%LSOUT,   gis_phy%COLAT1,  gis_phy%CFHOUR1,       &
! FIM-specific arguments
! TODO:  refactor to remove these
                 gis_phy%CLDCOV,                                          &
                 gis_phy%LEVS,    gis_phy%LATS_NODE_R, gis_phy%NMTVR,     &
                 gis_phy%num_p3d, gis_phy%num_p2d, gis_phy%NFXR,          &
                 nip, gis_phy%lsoil, GravityWaveDrag, CallRadiation,      &
                 yyyymmddhhmm,inv_perm,                                   &
                 skip_cu_physics, skip_mp_physics,skip_chem, ipn,sscal,   &
                 ext_cof,asymp,extlw_cof)

if (chem_opt.gt.0) then
!sms$compare_var(pb2d,'end physics')
!sms$compare_var(exch,'end physics')
!sms$compare_var(phys3dwrf,'end physics')
endif
  
endif ! CallPhysics

return
end subroutine physics
end module module_physics
