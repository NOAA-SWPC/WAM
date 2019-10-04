module module_fim_chem_finalize
contains
!*********************************************************************
subroutine chem_finalize
!       Finish the atmospheric chemistry component.  
!       T. Henderson            April, 2008
!*********************************************************************

  use module_outtime_chem,only: OutTime
  use module_control     ,only: PrintMAXMINtimes
  use module_initial_chem_namelists, only: chem_opt

  implicit none

  ! print elapsed times for atmospheric chemistry
  if (chem_opt >= 300 ) then
    call OutTime(PrintMAXMINtimes)
  endif

  return
end subroutine chem_finalize
end module module_fim_chem_finalize
