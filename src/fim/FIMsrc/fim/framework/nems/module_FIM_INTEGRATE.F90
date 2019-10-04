module module_FIM_INTEGRATE
  use esmf_mod
  use module_err_msg

  ! TODO:  move this to internal state
  use module_core_setup     ,only: use_write_tasks
  ! TODO:  not sure if this needs to move
  use icosio,only:icosio_end_frame

  implicit none
  private

  public :: fim_integrate

CONTAINS

! only FIM compute tasks execute this routine
  subroutine fim_integrate (gc_fim_dyn,  &
                            gc_fim_phy,  &
                            gc_fim_cpl,  &
                            imp_fim_dyn, &
                            exp_fim_dyn, &
                            imp_fim_phy, &
                            exp_fim_phy, &
                            clock_fim, &
                            rc_integrate)

    type(esmf_gridcomp), intent(inout) :: gc_fim_dyn
    type(esmf_gridcomp), intent(inout) :: gc_fim_phy
    type(esmf_cplcomp),  intent(inout) :: gc_fim_cpl
    type(esmf_state),    intent(inout) :: imp_fim_dyn
    type(esmf_state),    intent(inout) :: exp_fim_dyn
    type(esmf_state),    intent(inout) :: imp_fim_phy
    type(esmf_state),    intent(inout) :: exp_fim_phy
    type(esmf_clock),    intent(inout) :: clock_fim
    integer,             intent(  out) :: rc_integrate
!
! Local variables
!
    integer :: rc
    integer(esmf_kind_i8) :: ntimestep_esmf
    integer :: ntimestep
    type(esmf_timeinterval) :: timestep
    type(esmf_time) :: stoptime, newstoptime

     ! Run the clock one more time step (i.e. stop after its=nts+1), then 
     ! back up one step to mimic run.F90.  
     !  * set stoptime = stoptime+dt
     !  * run "integrate" loop
     !  * set ESMF_MODE_REVERSE
     !  * advance backwards one time step
     !  * set ESMF_MODE_FORWARD
     !  * reset stoptime to its original value
     !NOTE:  This hackery works around the fact that the original 
     !NOTE:  FIM run.F90 executes one extra time step in which the 
     !NOTE:  dynamics component finishes its final computations.  This 
     !NOTE:  was required by early versions of NEMS which did not 
     !NOTE:  allow multiple run phases.  See run.F90 for a very 
     !NOTE:  detailed discussion of this issue.  
     !NOTE:  This complexity could be avoided if we allowed a 2-phase 
     !NOTE:  run method for the DYN component -- and run.F90 would 
     !NOTE:  also be simplified.  However, interoperability with other 
     !NOTE:  components would be more difficult due to potential 
     !NOTE:  mismatches in numbers of phases.  

     call esmf_clockget (clock=clock_fim, stoptime=stoptime, rc=rc)
     call err_msg (rc,'esmf_clockget(stoptime)', rc_integrate)
     call esmf_clockget (clock=clock_fim, timestep=timestep, rc=rc)
     call err_msg (rc,'esmf_clockget(timestep)', rc_integrate)
     newstoptime = stoptime + timestep
     call esmf_clockset(clock=clock_fim, stoptime=newstoptime, rc=rc)
     call err_msg (rc,'esmf_clockset(newstoptime)', rc_integrate)

     integrate: do while (.not. esmf_clockisstoptime (clock_fim, rc=rc))
      call err_msg (rc,'esmf_clockisstoptime', rc_integrate)
!
!-----------------------------------------------------------------------
!***  Execute the Run step of the Dynamics component.
!-----------------------------------------------------------------------
!
      call esmf_logwrite("execute fim dynamics", esmf_log_info, rc=rc)
      call esmf_gridcomprun (gridcomp   =gc_fim_dyn,     &
                             importstate=imp_fim_dyn,    &
                             exportstate=exp_fim_dyn,    &
                             clock      =clock_fim,      &
                             rc         =rc)
      call err_msg (rc,'execute fim dynamics', rc_integrate)
!
!-----------------------------------------------------------------------
!***  Bring export data from the Dynamics into the coupler
!***  and export it to the Physics.
!-----------------------------------------------------------------------
!
      call esmf_logwrite ("couple dyn_exp-to-phy_imp", esmf_log_info, rc=rc)
      call esmf_cplcomprun (cplcomp    =gc_fim_cpl,          &
                            importstate=exp_fim_dyn,         &
                            exportstate=imp_fim_phy,         &
                            clock      =clock_fim,            &
                            rc         =rc)
      call err_msg (rc,'couple dyn_exp-to-phy_imp', rc_integrate)
!
!-----------------------------------------------------------------------
!***  Execute the Run step of the Physics Component.
!-----------------------------------------------------------------------
!
      call esmf_logwrite ("execute physics", esmf_log_info, rc=rc)
      call esmf_gridcomprun (gridcomp   =gc_fim_phy,            &
                             importstate=imp_fim_phy,           &
                             exportstate=exp_fim_phy,           &
                             clock      =clock_fim,             &
                             rc         =rc)
      call err_msg (rc, 'execute physics', rc_integrate)
!
!-----------------------------------------------------------------------
!***  Bring export data from the Physics into the coupler
!***  and export it to the Dynamics.
!-----------------------------------------------------------------------
!
      call esmf_logwrite ("couple phy_exp-to-dyn_imp", esmf_log_info, rc=rc)
      call esmf_cplcomprun (cplcomp    =gc_fim_cpl,             &
                            importstate=exp_fim_phy,            &
                            exportstate=imp_fim_dyn,            &
                            clock      =clock_fim,              &
                            rc         =rc)
      call err_msg (rc, 'couple phy_exp-to-dyn_imp', rc_integrate)
!
!-----------------------------------------------------------------------
!***  Flush buffered output to write tasks and/or clear list of files
!***  written to during this output frame.      
!-----------------------------------------------------------------------
!
      call esmf_clockget (clock=clock_fim,advancecount=ntimestep_esmf,rc=rc)
      call err_msg (rc, 'get time step from clock', rc_integrate)
      ntimestep = ntimestep_esmf
      call icosio_end_frame(ntimestep)
!
!-----------------------------------------------------------------------
!***  Advance clock to next time step.
!-----------------------------------------------------------------------
!
      call esmf_clockadvance (clock=clock_fim, rc=rc)
      call err_msg (rc, 'advance clock', rc_integrate)

     end do integrate    ! time step loop

     ! reset clock to state expected by caller upon return
     call esmf_clockset(clock=clock_fim, direction=ESMF_MODE_REVERSE, rc=rc)
     call err_msg (rc,'esmf_clockset(ESMF_MODE_REVERSE)', rc_integrate)
     call esmf_clockadvance(clock=clock_fim, rc=rc)
     call err_msg (rc,'esmf_clockadvance(one step backwards)', rc_integrate)
     call esmf_clockset(clock=clock_fim, direction=ESMF_MODE_FORWARD, rc=rc)
     call err_msg (rc,'esmf_clockset(ESMF_MODE_FORWARD)', rc_integrate)
     call esmf_clockset(clock=clock_fim, stoptime=stoptime, rc=rc)
     call err_msg (rc,'esmf_clockset(restore original stoptime)', rc_integrate)

  end subroutine fim_integrate
end module module_FIM_INTEGRATE
