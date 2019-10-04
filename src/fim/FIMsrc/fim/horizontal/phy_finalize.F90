module module_fim_phy_finalize
contains
!*********************************************************************
subroutine phy_finalize
!       Finish the physics component.  
!       T. Henderson            April, 2008
!*********************************************************************

  use module_outtime_phy,only: OutTime
  use module_control    ,only: PrintMAXMINtimes

  implicit none

  ! print elapsed times for physics parts of FIM
  call OutTime(PrintMAXMINtimes)

  return
end subroutine phy_finalize
end module module_fim_phy_finalize
