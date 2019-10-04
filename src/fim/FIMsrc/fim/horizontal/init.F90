!*********************************************************************
!       Loads the initial variables and constants to start fim
!       Alexander E. MacDonald  11/27/04
!       J. Lee                  September, 2005
!*********************************************************************

subroutine init
  use module_core_setup            ,only: core_setup_fim, iam_fim_task, iam_write_task
  use module_fim_chem_init         ,only: chem_init
  use module_fim_cpl_init          ,only: cpl_init
  use module_fim_dyn_init          ,only: dyn_init
  use module_fim_phy_init          ,only: phy_init
  use module_wrf_phy_init          ,only: wrf_phy_init
  use module_control               ,only: readrestart
  use restart                      ,only: read_restart

  implicit none

! Local variables

!JR: Sit in a spin-wait loop so a debugger can attach, halt the process,
!JR  reset the variable spinwait, and then continue.
!JR: Placed here because MPI is now active (if enabled), but the model is
!JR  still in the startup phase.

#ifdef ATTACH_DEBUGGER
  integer :: spinwait = 1

  do while (spinwait == 1)
  end do
#endif

! When MPI is used, set up communicators for compute tasks and 
! optional write tasks.  Mirrors NEMS approach.  
! NOTE:  Executable SMS directives must not be placed before 
! NOTE:  this call!  
! NOTE:  This includes writes or prints without !SMS$ignore because they 
! NOTE:  cause SMS to generate code.

  call core_setup_fim ()

  if (iam_fim_task .or. iam_write_task) then
    call dyn_init ()
  end if
  
  if (iam_fim_task) then
    call phy_init ()
    call chem_init ()
    call wrf_phy_init ()
    call cpl_init ()
    if (readrestart) then
      call read_restart ()          ! read dynamics and physics data from restart file
    end if
  end if
  
  return
end subroutine init
