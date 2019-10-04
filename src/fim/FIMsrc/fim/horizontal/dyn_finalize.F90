module module_fim_dyn_finalize

contains

!*********************************************************************
  subroutine dyn_finalize
!       Finish the dynamics component.  
!       T. Henderson            April, 2008
!*********************************************************************

    use module_control,        only: PrintMAXMINtimes,TimingBarriers,nts
    use module_outtime_dyn,    only: OutTime
    use module_core_setup,     only: use_write_tasks

    implicit none

    ! print elapsed times for dynamics parts of FIM
    call OutTime(maxmin_times=PrintMAXMINtimes,TimingBarriers=TimingBarriers, &
      print_header=.true.)

    return

  end subroutine dyn_finalize

end module module_fim_dyn_finalize
