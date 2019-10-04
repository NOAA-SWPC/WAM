module module_wrf_phy_run
contains
!*********************************************************************
        subroutine wrf_phy_run(its)
!       "Run" method for WRF physics
!*********************************************************************

use module_control         ,only: nts, itsStart
use module_wrfphysics      ,only: wrf_physics
use module_outtime_wrf_phy ,only: telapsed=>tphy
use module_initial_chem_namelists, only: cu_physics, mp_physics

implicit none

!  Declare dummy arguments
integer, intent(in) :: its

!  Declare local variables:
real*8 :: t0

call StartTimer(t0)

  !...........................................................
  ! Advance the physics component by one time step unless this 
  ! is the last (nts+1) iteration.  
  ! This complexity is required for the NCEP ESMF approach 
  ! in which single-phase DYN and PHY components alternate 
  ! execution during each time step.  
  !
if (its < itsStart+nts) then

  !...........................................................
  ! call WRF phyics
  !...........................................................
  if ((.not.mp_physics==0).or.(.not.cu_physics==0)) then
!   write(6,*)'call wrfphysics ',mp_physics,cu_physics
    call wrf_physics(its)
  endif

endif

call IncrementTimer(t0,telapsed)

return
end subroutine wrf_phy_run
end module module_wrf_phy_run
