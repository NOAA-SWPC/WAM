!*********************************************************************
	subroutine finalize
!	Stop program for icosahedral flow-following global model
!	Alexander E. MacDonald  12/27/2004
!*********************************************************************

use module_control               ,only: PrintMAXMINtimes
use module_core_setup            ,only: iam_fim_task,iam_write_task
use module_fim_chem_finalize     ,only: chem_finalize
use module_fim_cpl_finalize      ,only: cpl_finalize
use module_fim_dyn_finalize      ,only: dyn_finalize
use module_fim_phy_finalize      ,only: phy_finalize
use module_fim_wrf_phy_finalize  ,only: wrf_phy_finalize
use module_initial_chem_namelists,only: chem_opt                           
use module_outtime_main          ,only: OutTime

implicit none

IF (iam_fim_task) THEN
  print*,' '
  call cpl_finalize
  call dyn_finalize
  call phy_finalize
  if (chem_opt.gt.0) then
    call wrf_phy_finalize
    call chem_finalize
  endif
  ! print elapsed times for main loop of FIM
  call OutTime(PrintMAXMINtimes)
  print*,' '
  ! call UnstructuredPrintTimers
  call datetime
  print*,'Program exited normally'
ENDIF

return

end subroutine finalize
