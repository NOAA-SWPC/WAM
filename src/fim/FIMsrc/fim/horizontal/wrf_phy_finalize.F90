module module_fim_wrf_phy_finalize
contains
!*********************************************************************
subroutine wrf_phy_finalize
!       Finish the WRF physics component.  
!       T. Henderson            April, 2008
!*********************************************************************

  use module_outtime_wrf_phy,only: OutTime
  use module_control     ,only: PrintMAXMINtimes
  use module_initial_chem_namelists, only: cu_physics, mp_physics

  implicit none

  ! print elapsed times for WRF physics
  if ((.not.mp_physics==0).or.(.not.cu_physics==0)) then
    call OutTime(PrintMAXMINtimes)
  endif

  return
end subroutine wrf_phy_finalize
end module module_fim_wrf_phy_finalize
