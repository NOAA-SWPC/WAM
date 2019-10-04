module module_fim_cpl_finalize
contains
!*********************************************************************
        subroutine cpl_finalize
!       Finish the FIM DYN-PHY coupler component.  
!       T. Henderson            February, 2009
!*********************************************************************

  use module_outtime_cpl,only: OutTime
  use module_control    ,only: PrintMAXMINtimes

  implicit none

  ! print elapsed times for coupler parts of FIM
  call OutTime(maxmin_times=PrintMAXMINtimes,print_header=.true.)

  return
end subroutine cpl_finalize
end module module_fim_cpl_finalize
