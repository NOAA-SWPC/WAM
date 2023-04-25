!#include "../../../ESMFVersionDefine.h"

! !module: gfs_dynamics_grid_comp_mod --- 
!                       esmf gridded component of gfs dynamics
!
! !description: gfs dynamics gridded component main module.
!
! !revision history:
!
!  january 2007     hann-ming henry juang initiated and wrote the code
!  March   2009     Weiyu Yang, modified for the ensemble NEMS run.
!  oct 4   2009     sarah lu, 3D Gaussian grid (DistGrid5) added
!  oct 5   2009     sarah lu, grid_gr unfolded from 2D to 3D
!  oct     2009     Jun Wang, add not quilting option, add run time variables into wrt inport state
!  oct 12 2009      Sarah Lu, set up the association between imp/exp state and grid_gr
!                   => export state is pointed to grid_gr in init step
!                   => grid_gr is pointed to import state in run step
!  oct 17 2009      Sarah Lu, add debug print to check imp/exp state
!  nov 09 2009      Jun Wang, add grid_gr_dfi for digital filter
!  Feb 05 2010      Jun Wang, add restart step
!  Feb 20 2011      Henry Juang, add non-iterating dimensional-splitting semi-Lagrangian (NDSL)
!                   advection with options of MASS_DP and NDSLFV
!  Feb    2011      Weiyu Yang, Updated to use both the ESMF 4.0.0rp2 library,
!                   ESMF 5 library and the the ESMF 3.1.0rp2 library.
!  May    2011      Weiyu Yang, Modified for using the ESMF 5.2.0r_beta_snapshot_07.
!  Sep    2011      Weiyu Yang, Modified for using the ESMF 5.2.0r library.
!  Sep    2012      Jun Wang, add sigio input option
!  Dec    2013      Jun Wang, add restart in gfs_dynamics_start_time_get_mod argument list
!  Feb 26 2016      S Moorthi add kdt_start and kdt_dif for grid-point digital filter
!                             logic. Filter works for initial state or for restart
!  Aug    2016      Add calls for iau forcing
!
! !interface:
!
      module gfs_dynamics_grid_comp_mod
 
!!uses:
!------
      USE ESMF

      use gfs_dynamics_err_msg_mod
      use gfs_dynamics_initialize_mod
      use gfs_dynamics_run_mod
      use gfs_dynamics_finalize_mod
      use gfs_dyn_resol_def,               only: kdt_start

      use gfs_dynamics_output,             only: point_dynamics_output_gfs
      use gfs_dyn_tracer_config,           only: gfs_dyn_tracer
      use namelist_dynamics_def,           only: nemsio_in, ldfi_spect
      USE gfs_dynamics_initialize_slg_mod, ONLY: gfs_dynamics_initialize_slg
      use gfs_dyn_iau_module, only : getiauforcing,applyiauforcing,gq_iau

      implicit none

      private   ! by default, data is private to this module

      public gfs_dyn_setservices      ! only set service is public

!eop
!-------------------------------------------------------------------------

      contains

!----------------------------------------------------------------------
!bop
!
! !routine: gfs_dyn_setservices --- 
!           set services for gfs dynamics gridded component.
! 
! !interface:
!
      subroutine gfs_dyn_setservices (gc_gfs_dyn, rc)
 
! !arguments:
!------------

      type(esmf_gridcomp)              :: gc_gfs_dyn 	! gridded component
      integer,             intent(out) :: rc    	! return code
     
! !description: set services (register) for the gfs dynamics grid component.
!         
!eop         
!----------------------------------------------------------------------
  
      integer                          :: rc1     = esmf_success

! initializing the error signal variable rc.
!-------------------------------------------
      rc = esmf_success

! register services for this component
! ------------------------------------

! register the initialize subroutine.  since it is just one subroutine
! for the initialize, no need to specify phase.  the second argument is
! a pre-defined subroutine type, such as esmf_setinit, esmf_setrun, 
! esmf_setfinal.
!---------------------------------------------------------------------
      call esmf_logwrite("set entry point for initialize",              &
                          ESMF_LOGMSG_INFO, rc = rc1)
      CALL ESMF_GridCompSetEntryPoint( gc_gfs_dyn                       &
                                     , ESMF_METHOD_INITIALIZE           &
                                     , gfs_dyn_initialize               &
                                     , rc=RC1)
      call gfs_dynamics_err_msg(rc1,'set entry point for initialize',rc)

! register the run subroutine.
!-----------------------------
      call esmf_logwrite("set entry point for run",              	&
                          ESMF_LOGMSG_INFO, rc = rc1)

      CALL ESMF_GridCompSetEntryPoint( gc_gfs_dyn                       &
                                     , ESMF_METHOD_RUN                  &
                                     , gfs_dyn_run                      &
                                     , rc=RC1)
      call gfs_dynamics_err_msg(rc1,'set entry point for run',rc)


! register the finalize subroutine.
!----------------------------------
      call esmf_logwrite("set entry point for finalize",                &
                          ESMF_LOGMSG_INFO, rc = rc1)

      CALL ESMF_GridCompSetEntryPoint( gc_gfs_dyn                       &
                                     , ESMF_METHOD_FINALIZE             &
                                     , gfs_dyn_finalize                 &
                                     , rc=RC1)
      call gfs_dynamics_err_msg(rc1,'set entry point for finalize',rc)

! check the error signal variable and print out the result.
!----------------------------------------------------------
      call gfs_dynamics_err_msg_final(rc1,				&
                        'setservice for gfs dynamics grid comp.',rc)

      end subroutine gfs_dyn_setservices


!----------------------------------------------------------------------
!bop
! !routine:  gfs_dyn_initialize --- initialize routine to initialize 
!                                   and set up the gfs running job.
!
! !description: this subroutine initializes the gfs running before
!               the main running loop.
!
!
! !revision history:
!
!  november 2004     weiyu yang initial code.
!  may      2005     weiyu yang for the updated gfs version.
!  february 2006     moorthi
!  february 2007     h.-m. h. juang
!  oct 12 2009       Sarah Lu, export state is pointed to grid_gr once and for all
!  November 2009     Weiyu Yang, Ensemble GEFS.
!  Sep      2012     Jun Wang, add sigio input option
!  Jul 22,  2012     S Moorthi - Moved VM to the beginning of the routine to control prints
!
! !interface:
!

! this argument list is a standard list for all the initialize,
! the run and finalize routines for an esmf system.
!--------------------------------------------------------------
      subroutine gfs_dyn_initialize(gc_gfs_dyn,                         &
                                    imp_gfs_dyn, exp_gfs_dyn, clock, rc)

! user code, for computations related to the esmf interface states.
!------------------------------------------------------------------
      use gfs_dynamics_states_mod, only : gfs_dynamics_import2internal, &
                                          gfs_dynamics_internal2export
      use gfs_dynamics_grid_create_mod
      USE GFS_AddParameterToStateMod
!
! !input/output variables and parameters:
!----------------------------------------

      type(esmf_gridcomp)                :: gc_gfs_dyn 
      type(esmf_state)                   :: imp_gfs_dyn, exp_gfs_dyn
      type(esmf_clock)                   :: clock

!
! !output variables and parameters:
!----------------------------------

      integer, intent(out) :: rc  

! !eop
!------------------------------------------------------------------------- 
 
! !working arrays and local parameters.  
!--------------------------------------
      type(gfs_dyn_wrap)                :: wrap         
! this wrap is a derived type which contains
! only a pointer to the internal state.  it is needed
! for using different architectures or compliers.
      type(gfs_dynamics_internal_state), pointer  :: int_state    
      type(esmf_vm)                      :: vm_local     
      type(esmf_timeinterval)            :: timestep, runduration       &
                                          , reftimeinterval
      type(esmf_time)                    :: starttime, stoptime, currtime
!jw
      type(esmf_state)                   :: imp_state_write  !<-- The write gc import state

      type(ESMF_DistGrid)                :: DistGrid5    ! the ESMF DistGrid.

      integer(kind=esmf_kind_i4)         :: yy, mm, dd  &! time variables for date
                                          , hh, mns, sec ! time variables for time
      integer(kind=esmf_kind_i8)         :: advancecount

      integer                            :: advancecount4,  timestep_sec    &
                                          , atm_timestep_s, dyn_timestep_s  &
                                          , rc1, rcfinal,   grib_inp        &
                                          , ifhmax, runduration_hour

      character(20)                      :: cfile, cfile2

! initialize the error signal variables.
!---------------------------------------
      rc1     = esmf_success
      rcfinal = esmf_success

! allocate the internal state pointer.
!-------------------------------------
      call esmf_logwrite("allocate the dyn internal state",                 &
                          ESMF_LOGMSG_INFO, rc = rc1)

      allocate(int_state, stat = rc1)

      call gfs_dynamics_err_msg(rc1,' - allocate the internal state',rc)

      wrap%int_state => int_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! get the grid component vm.
! this esmf_gridcompget vm can be used at any where you need it.
!---------------------------------------------------------------
!      call esmf_logwrite("get the local vm", ESMF_LOGMSG_INFO, rc = rc1)

      call esmf_vmgetcurrent(vm_local, rc = rc1)

      call gfs_dynamics_err_msg(rc1,'get the vm',rc)

! set up parameters of mpi communications.
! use esmf utility to get pe identification and total number of pes.
!-------------------------------------------------------------------
!      call esmf_logwrite("get me and nodes from vm",                   &
!                          ESMF_LOGMSG_INFO, rc = rc1)

      call esmf_vmget(vm_local, localpet = int_state%me,                &
                      mpicommunicator = mpi_comm_all,                   &
                      petcount = int_state%nodes, rc = rc1)
!     mc_comp = mpi_comm_all

      call gfs_dynamics_err_msg(rc1,'get me and nodes from vm',rc)
!     write(0,*)'in dyn_gc,after vmget,npes=',int_state%nodes,'mpi_comm_all=',mpi_comm_all

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!jws
!-----------------------------------------------------------------------
!***  RETRIEVE THE IMPORT STATE OF THE WRITE GRIDDED COMPONENT
!***  FROM THE DYNAMICS EXPORT STATE.
!-----------------------------------------------------------------------
!
      call esmf_logwrite("get write gc import state",                  &
                          ESMF_LOGMSG_INFO, rc = rc1)

      CALL ESMF_StateGet(state      =exp_gfs_dyn                        &  !<-- The Dynamics export state
                        ,itemName   ='Write Import State'               &  !<-- Name of the state to get from Dynamics export state
                        ,nestedState=IMP_STATE_WRITE                    &  !<-- Extract write component import state from Dynamics export
                        ,rc         =RC)
      call gfs_dynamics_err_msg(rc1,'get write gc import state',rc)
!jwe

! attach internal state to the gfs dynamics grid component.
!-------------------------------------------------
      call esmf_logwrite("set up the internal state",                   &
                          ESMF_LOGMSG_INFO, rc = rc1)

      call esmf_gridcompsetinternalstate(gc_gfs_dyn, wrap, rc1)

      call gfs_dynamics_err_msg(rc1,'set up the internal state',rc)

! use esmf utilities to get information from the configuration file.
! the function is similar to reading the namelist in the original gfs.
!---------------------------------------------------------------------
      call esmf_logwrite("getting information from the configure file", &
                          ESMF_LOGMSG_INFO, rc = rc1)

      call gfs_dynamics_getcf(gc_gfs_dyn, int_state,  rc1)

      call gfs_dynamics_err_msg(rc1,'get configure file information',rc)

! get the start time from reading the sigma file.
!----------------------------------------------------------
      call esmf_logwrite("getting the start time",                      &
                          ESMF_LOGMSG_INFO, rc = rc1)

      if (int_state%me == 0) print *,'nemsio_in=',int_state%nemsio_in
      if(int_state%nemsio_in) then
        cfile  = trim(int_state%nam_gfs_dyn%grid_ini)
        cfile2 = trim(int_state%nam_gfs_dyn%grid_ini2)
      else
        cfile  = trim(int_state%nam_gfs_dyn%sig_ini)
        cfile2 = trim(int_state%nam_gfs_dyn%sig_ini2)
      endif

      call gfs_dynamics_start_time_get(                                 &
                        yy, mm, dd, hh, mns, sec, int_state%kfhour,     &
                        int_state%n1,int_state%n2,int_state%grib_inp,   &
                        cfile,cfile2,int_state%nemsio_in,               &
                        int_state%me,                                   &
                        int_state%restart_run,rc1)
      call gfs_dynamics_err_msg(rc1,'getting the start time',rc)
 
      advancecount4        = nint(real(int_state%kfhour) * 3600.0 /     &
                                       int_state%nam_gfs_dyn%deltim)
      int_state%phour      = advancecount4 *                            &
                             int_state%nam_gfs_dyn%deltim / 3600.0
      int_state%kfhour     = nint(int_state%phour)
!
      int_state%kdt        = advancecount4
      kdt_start            = advancecount4

!      write(0,*)'in dyn_grid_comp,advancecount4=',advancecount4,          &
!        'phour=',int_state%phour,'kfhour=',int_state%kfhour,'kdt=',     &
!         int_state%kdt
!
! initialize the clock with the start time based on the information
! from calling starttimeget.
!------------------------------------------
      call esmf_logwrite("set up the esmf time",                        &
                          ESMF_LOGMSG_INFO, rc = rc1)

      call esmf_timeset(starttime, yy = yy, mm = mm,  dd = dd,          &
                        h  = hh, m  = mns, s  = sec, rc = rc1)

      call gfs_dynamics_err_msg(rc1,'set up the esmf time',rc)

      call esmf_logwrite("set up the reference time interval",          &
                          ESMF_LOGMSG_INFO, rc = rc1)

      call esmf_timeintervalset(reftimeinterval, h = int_state%kfhour,  &
                                m = 0, rc = rc1)

! re-set up the start time based on the kfhour value in the sigma file.
!----------------------------------------------------------------------
!      starttime = starttime + reftimeinterval

!     call gfs_dynamics_err_msg(rc1,'set up the reference time interval',rc)

! set the esmf clock which will control the gfs run do loop.
!--------------------------------------------------------------

      currtime = starttime + reftimeinterval
      call esmf_clockset(clock, currtime = currtime, rc = rc1)
!
! initialize the gfs, including set up the internal state
! variables and some local parameter short names, aloocate
! internal state arrays.
!---------------------------------------------------------
      call esmf_logwrite("run the gfs_dynamics_initialize", 		&
                          ESMF_LOGMSG_INFO, rc = rc1)

! ======================================================================
! ----------------- gfs dynamics related initialize --------------------
! ======================================================================
! grid_gr unfolded (sarah lu)
      IF(int_state%slg_flag) THEN
        CALL gfs_dynamics_initialize_slg(int_state, rc1)
      ELSE
        call gfs_dynamics_initialize(int_state, rc1)
      ENDIF
!      write(0,*)'in dyn_init, t=',maxval(int_state%grid_gr(:,int_state%g_t)), &
!       minval(int_state%grid_gr(:,int_state%g_t)),'quilting=',quilting
! ======================================================================
! ----------------------------------------------------------------------
! ======================================================================

      call gfs_dynamics_err_msg(rc1,'run the gfs_dynamics_initialize',rc)

      call esmf_clockget(clock, timestep    = timestep,            	&
                         runduration = runduration,              	&
                         starttime   = starttime,                	&
                         currtime    = currtime,                 	&
                         rc          = rc1)
!
      call esmf_timeintervalget(runduration, h = runduration_hour, rc = rc1)
!
!moor ifhmax = nint(int_state%nam_gfs_dyn%fhmax)
      ifhmax = nint(fhmax)
      if (runduration_hour <= 0    .or.  ifhmax /= 0   .and.  		&
          ifhmax <= int_state%kfhour + runduration_hour) then
          ifhmax            = nint(fhmax)
          runduration_hour  = nint(fhmax) - nint(fhini)
          call esmf_timeintervalset(runduration, h = runduration_hour, rc = rc1)
      endif
      if (runduration_hour <= 0) then
        write(0,*)'WRONG: fhini=',fhini, ' >= fhmax=',fhmax,          &
                  ' job aborted with error code 444'
        if(me == 0)  call mpi_quit(444)
      endif
      stoptime = currtime  + runduration
                           
      call esmf_clockset(clock, stoptime = stoptime,        rc = rc1)
!
      call esmf_timeintervalget(timestep, s = timestep_sec, rc = rc1)
                           
!!
      if (me == 0) then
        call info_out_para(real(timestep_sec))
      endif
!!
      if (me == 0) then
        print *,' the gsm will forecast ',runduration_hour,' hours',    &
                ' from hour ',int_state%kfhour,' to hour ',             &
                 runduration_hour+int_state%kfhour
      endif
!
!
! create 3D Gaussian grid  (sarah lu)
!-----------------------
!
      call gfs_dynamics_grid_create_Gauss3D(vm_local,int_state,DistGrid5,rc1)

      call gfs_dynamics_err_msg(rc1,'gfs_dynamics_grid_create_gauss3d',rc)

      int_state%fhour_idate(1,1)   = fhour
      int_state%fhour_idate(1,2:5) = idate(1:4)
!
      IF(int_state%ENS) THEN
          int_state%end_step = .true.

          CALL AddParameterToState(exp_gfs_dyn, int_state, rc = rc1)

          call gfs_dynamics_err_msg(rc1, 'Add Parameter To export State',rc)
      END IF

!
! Define Dynamics Export states    (Sarah Lu)
!
      call gfs_dynamics_internal2export(int_state, exp_gfs_dyn, rc1)

      call gfs_dynamics_err_msg(rc1,'gfs_dynamics_internal2export',rc)
!
!-------------------------------------------------------
! send all the head info to write tasks
!-------------------------------------------------------
!
      call point_dynamics_output_gfs(int_state,IMP_STATE_WRITE)

!
!*******************************************************************
! print out the final error signal variable and put it to rc.
!------------------------------------------------------------
      call gfs_dynamics_err_msg_final(rcfinal,				&
                        'initialize from gfs dynamics grid comp.',rc)

      end subroutine gfs_dyn_initialize



!----------------------------------------------------------------------
!bop
!
! !routine: gfs_dyn_run --- 
!           main grid component routine to run the gfs dynamics.
!
! !description: this subroutine will run the most part computations 
!               of the gfs dynamics.
!
! !revision history:
!
!  november 2004     weiyu yang initial code.
!  may      2005     weiyu yang for the updated gfs version.
!  february 2006     moorthi
!  july     2007     hann-ming henry juang
!  oct 12 2009       Sarah Lu, point grid_gr to import state once and for all
!  oct 17 2009       Sarah Lu, debug print added to track imp/exp states
!  November 2009     Weiyu Yang, Ensemble GEFS.
!  nov 09 2009       Jun Wang, get data from grid_gr_dfi to internal state for dfi
!  feb 05 2010       Jun Wang, set restart step
!
! !interface:
!

      subroutine gfs_dyn_run(gc_gfs_dyn, imp_gfs_dyn, exp_gfs_dyn, clock, rc)

      use gfs_dynamics_states_mod, only : gfs_dynamics_internal2export, &
                                          gfs_dynamics_import2internal
!
! !input variables and parameters:
!---------------------------------
      type(esmf_gridcomp)                :: gc_gfs_dyn   
      type(esmf_state)                   :: imp_gfs_dyn 
 
! !output variables and parameters:
!----------------------------------
      type(esmf_clock)                   :: clock
      type(esmf_timeinterval)            :: timestep,  donetime    
      type(esmf_time)                    :: starttime, currtime         &
                                          , stoptime,  dfitime
      type(esmf_state)                   :: exp_gfs_dyn
      integer,             intent(out)   :: rc   
!eop
!-------------------------------------------------------------------------

!
! !working arrays and local parameters.
!--------------------------------------
      type(gfs_dyn_wrap)                :: wrap         
! this wrap is a derived type which contains
! only a pointer to the internal state.  it is needed
! for using different architectures or compliers.
      type(gfs_dynamics_internal_state), pointer  :: int_state   
!
      type(esmf_state)                  :: imp_state_write  !<-- The write gc import state
      logical,save                      :: first_reset=.true.
      logical,save                      :: first_dfiend=.true.
      TYPE(ESMF_TimeInterval)           :: HALFDFIINTVAL
      integer                            ::  timestep_sec
      real(kind_evod)  dtiau

!! debug print for tracking import and export state (Sarah Lu)
      TYPE(ESMF_Field)                   :: ESMFField             !chlu_debug
      TYPE(ESMF_FieldBundle)             :: ESMFBundle            !chlu_debug
      TYPE(ESMF_Info)                    :: info
      real(kind_grid), dimension(:,:,:),allocatable ::  grid_gr_iau
      REAL(ESMF_KIND_R8), DIMENSION(:,:,:), POINTER :: fArr3D     !chlu_debug
      integer                            :: rc1, rcfinal, DFIHR &
                                          , localPE,ii1,ii2,ii3 & !chlu_debug
                                          , n, k, rc2, kdt_dif
!     logical, parameter                 :: ckprnt = .true        !chlu_debug
      logical, parameter                 :: ckprnt = .false.      !chlu_debug
      integer, parameter                 :: item_count = 3        !chlu_debug
      character(5) :: item_name(item_count)                       !chlu_debug
      character(20) :: vname                                      !chlu_debug
      data item_name/'t','u','v'/                                 !chlu_debug

      localPE = 0                                                 !chlu_debug

! initialize the error signal variables.
!---------------------------------------

      rc1     = esmf_success
      rcfinal = esmf_success

! retrieve the esmf internal state.
!---------------------------------- 
      call esmf_logwrite("get the internal state in the run routine", 	&
                          ESMF_LOGMSG_INFO, rc = rc1)

      call esmf_gridcompgetinternalstate(gc_gfs_dyn, wrap, rc1)

      call gfs_dynamics_err_msg(rc1,'get the internal state in the run routine',rc)

! pointing the local internal state pointer to the esmf internal state pointer.
!------------------------------------------------------------------------------
      int_state => wrap%int_state

!
! get the esmf import state and over-write the gfs internal state.
! update the initial condition arrays in the internal state based on
! the information of the esmf import state.
!------------------------------------------------------------------
      call esmf_logwrite("esmf import state to internal state",         &
                          ESMF_LOGMSG_INFO, rc = rc1)
!
      kdt_dif = int_state%kdt - kdt_start
      int_state%reset_step = .false.
      if(int_state%restart_step ) first_reset = .false.
      int_state%dfiend_step = .false.
      if( int_state%ndfi > 0 .and. first_reset .and.                    &
                                   kdt_dif == int_state%ndfi) then
        if( first_dfiend ) then           ! first go through dfi step
          int_state%dfiend_step = .true.
          first_dfiend          = .false.
        else                              ! second go through reset step
          int_state%reset_step  = .true.
          int_state%dfiend_step = .false.
          first_reset           = .false.
        endif
      endif

!      if (me == 0) &
!       print *,'in grid comp,ndfi=',int_state%ndfi,'kdt=',int_state%kdt, &
!        'first_reset=',first_reset,'reset_step=',int_state%reset_step,  &
!        ' kdt_dif=',kdt_dif,' kdt_start=',kdt_start, &
!        'dfiend_step=',int_state%dfiend_step,        &
!        ' start_step=',int_state%start_step,' restart_step=',int_state%restart_step
!
      IF(.NOT. int_state%restart_step .AND. .NOT. int_state%start_step ) THEN
        IF (iau) THEN
            allocate(grid_gr_iau(lonf,lats_node_a_max,gq_iau))
            call esmf_clockget(clock, timestep    = timestep, rc = rc1)
            call esmf_timeintervalget(timestep, s = timestep_sec, rc = rc1)
            IF (semilag) THEN
               dtiau = timestep_sec
            ELSE
               dtiau = 2.*timestep_sec
            ENDIF
            int_state%iniauinterval=.false.
            call getiauforcing(grid_gr_iau,int_state%phour*3600.,int_state,rc)
            IF (rc == 0) THEN
               IF (me == 0) print*,'applying IAU forcing'
               int_state%iniauinterval=.true.
               call applyiauforcing(grid_gr_iau,imp_gfs_dyn,int_state,dtiau,rc)
            ENDIF         
            deallocate(grid_gr_iau)
        ENDIF 
        IF(.NOT. int_state%reset_step) THEN
          CALL gfs_dynamics_import2internal(imp_gfs_dyn, int_state, rc1)
        ELSE
          if(.NOT.ldfi_spect) then
            CALL gfs_dynamics_import2internal(imp_gfs_dyn, int_state,       &
                                              rc=rc1, exp_gfs_dyn=exp_gfs_dyn)
          endif
        END IF 

        CALL gfs_dynamics_err_msg(rc1, 'esmf import state to internal state', rc)
        idate(1:4) = int_state%fhour_idate(1, 2:5)
      END IF
!
! get clock times
! ------------------
      call esmf_clockget(clock, timestep=timestep, starttime=starttime, &
                         currtime=currtime, stoptime=stoptime, rc=rc1)

      call gfs_dynamics_err_msg(rc1,'esmf clockget',rc)

      donetime = currtime - starttime

      int_state%kdt = nint(donetime/timeStep)        ! current time step

!     if (me == 0) &
!     write(0,*)'in grid comp,ndfi=',int_state%ndfi,'kdt=',int_state%kdt
!

! Set up the ensemble coupling time flag.
!----------------------------------------
      CALL ESMF_InfoGetFromHost(imp_gfs_dyn, info, rc = rc1)
      CALL ESMF_InfoSet(info, 'Cpl_flag', int_state%Cpl_flag, rc = rc1)

      if( currtime  == stoptime ) then
!         print *,' currtime equals to stoptime '
          int_state%end_step = .true.
      else
          int_state%end_step = .false.
      endif
!
! get forecast date
      call esmf_timeget(currtime,                                         &
                        yy=int_state%nfcstdate7(1),                       &
                        mm=int_state%nfcstdate7(2),                       &
                        dd=int_state%nfcstdate7(3),                       &
                        h =int_state%nfcstdate7(4),                       &
                        m =int_state%nfcstdate7(5),                       &
                        s =int_state%nfcstdate7(6),                       &
                        rc=rc1)
      call gfs_dynamics_err_msg(rc1,'esmf timeget',rc)
!      print *,'in gfs grid comp,currtime=',int_state%nfcstdate7(1:6)
!
! ======================================================================
! --------------- run the gfs dynamics related -------------------------
! ======================================================================
      call esmf_logwrite("run the gfs_dynamics_run", ESMF_LOGMSG_INFO, rc = rc1)

      call gfs_dynamics_run(int_state, imp_gfs_dyn, rc = rc1)
      call gfs_dynamics_err_msg(rc1,'run the gfs_dynamics_run',rc)
! ======================================================================
! ======================================================================

! transfer the gfs export fields in the internal state 
! to the esmf export state which is the public interface
! for other esmf grid components. link is done in initialize, so do not need.
!-------------------------------------------------------
     call esmf_logwrite("internal state to esmf export state", 	&
                         ESMF_LOGMSG_INFO, rc = rc1)

! Need to check it?  Should be removed?
!--------------------------------------
     call gfs_dynamics_internal2export(int_state, exp_gfs_dyn, rc = rc1)

     call gfs_dynamics_err_msg(rc1,'internal state to esmf export state',rc)

!! debug print starts here  (Sarah Lu) -----------------------------------
!  -----------------------------------

      lab_if_ckprnt_ex : if ( ckprnt .and. (int_state%me ==0) ) then      !chlu_debug
        do n = 1, item_count                                              !chlu_debug
            vname = trim(item_name(n))                                    !chlu_debug
            if(associated(fArr3D)) nullify(fArr3D)                        !chlu_debug
            CALL ESMF_StateGet(state     = exp_gfs_dyn                  & !chlu_debug
                              ,itemName  = vname                        & !chlu_debug
                              ,field     = ESMFField                    & !chlu_debug
                              ,rc        = rc1)                           !chlu_debug
            call gfs_dynamics_err_msg(rc1,'LU_DYN: get ESMFarray',rc)     !chlu_debug

            CALL ESMF_FieldGet(field=ESMFField, localDe=0, &              !chlu_debug
                               farrayPtr=fArr3D, rc = rc1)                !chlu_debug

            call gfs_dynamics_err_msg(rc1,'LU_DYN: get F90array',rc)       !chlu_debug
!            ii1 = size(fArr3D, dim=1)                                     !chlu_debug
!            ii2 = size(fArr3D, dim=2)                                     !chlu_debug
!            ii3 = size(fArr3D, dim=3)                                     !chlu_debug
!            if(n==1) print *, 'LU_DYN:',ii1, 'x', ii2, 'x', ii3           !chlu_debug
!            print *,' LU_DYN: exp_: ', vname                           &   !chlu_debug
!                   , fArr3D(1,1,1), fArr3D(1,2,1)                     &   !chlu_debug
!                   , fArr3D(2,1,1), fArr3D(ii1,ii2,ii3)                    !chlu_debug
        enddo                                                              !chlu_debug


        CALL ESMF_StateGet(state=exp_gfs_dyn, ItemName='tracers',   &
                           fieldbundle=ESMFBundle, rc=rc1 )
        CALL gfs_dynamics_err_msg(rc1, 'LU_DYN: get Bundle from exp state',rc)

!       write(0,*)' after stateget for bundle in dyngridcomp'

        do n = 1, int_state%ntrac                                          !chlu_debug
!         vname = int_state%gfs_dyn_tracer%vname(n, 1)                     !chlu_debug
          vname = trim(gfs_dyn_tracer%vname(n, 1))                         !chlu_debug
          if(associated(fArr3D)) nullify(fArr3D)                           !chlu_debug
!         write(0,*) 'LU_DYN:',vname                                       !chlu_debug
          CALL ESMF_FieldBundleGet(ESMFBundle, vname,               &      !chlu_debug
                                   field = ESMFField, rc = rc1)
          CALL ESMF_FieldGet(field=ESMFfield, farrayPtr=fArr3D, localDe=0, rc=rc1)
          call gfs_dynamics_err_msg(rc1,'LU_DYN: get ESMFfield ',rc)       !chlu_debug

!          ii1 = size(fArr3D, dim=1)                                     !chlu_debug
!          ii2 = size(fArr3D, dim=2)                                     !chlu_debug
!          ii3 = size(fArr3D, dim=3)                                     !chlu_debug
!          if(n == 1) print *,'LU_DYN:',ii1, 'x', ii2, 'x', ii3          !chlu_debug
!          write(0,*)'LU_DYN: exp_:',vname,' n=',n,                 &
!                maxval(fArr3D(1:ii1,1:ii2,1:ii3)),                 &
!                minval(fArr3D(1:ii1,1:ii2,1:ii3)), ' n=',n
!          print *,'LU_DYN: exp_:',vname, fArr3D(1,1,1)        &         !chlu_debug
!               , fArr3D(1,2,1), fArr3D(2,1,1), fArr3D(ii1,ii2,ii3)      !chlu_debug
        enddo                                                            !chlu_debug

      endif lab_if_ckprnt_ex                                             !chlu_debug
!! -------------------------------------- debug print ends here  (Sarah Lu)

! ======================================================================
!------------- put run level variables into write_imp_state---------
! ======================================================================
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
     call esmf_logwrite("get imp_state_write from esmf export state",   &
                         ESMF_LOGMSG_INFO, rc = rc1)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateGet(state          =exp_gfs_dyn    &  !<-- The Dyn component's export state
                        ,itemName       ="Write Import State"     &  !<-- Name of state to be extracted
                        ,nestedState    =IMP_STATE_WRITE  &  !<-- The extracted state
                        ,rc             =RC1)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
     call gfs_dynamics_err_msg(rc1,'get imp_state_write from esmf export state',rc)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
     call esmf_logwrite("set pdryini in imp_state_write",               &
                          ESMF_LOGMSG_INFO, rc = rc1)
     CALL ESMF_InfoGetFromHost(IMP_STATE_WRITE, info, rc=rc1)
     call gfs_dynamics_err_msg(rc1,'retrieve info from imp_state_write',rc)
     CALL ESMF_InfoSet(info                                             &  !<-- The Write component import state's info handle
                       ,key      ='pdryini'                             &  !<-- Name of the var
                       ,value    =int_state%pdryini                     &  !<-- The var being inserted into the import state
                       ,rc       =RC1)
     call gfs_dynamics_err_msg(rc1,'set pdryini in imp_state_write',rc)
!*******************************************************************
!
! print out the final error signal information and put it to rc.
!---------------------------------------------------------------
      call gfs_dynamics_err_msg_final(rcfinal,                          &
                        'run from gfs dynamics grid comp.',rc)
      end subroutine gfs_dyn_run


!----------------------------------------------------------------------
!bop
!
! !routine: finalize --- finalizing routine to finish the 
!                        gfs running job.
!
! !description: this subroutine will finish the gfs computations,
! !             and will release the memory space.
!
! !revision history:
!
!  november 2004     weiyu yang initial code.
!  may      2005     weiyu yang for the updated gfs version.
!  february 2006     moorthi
!  february 2007     juang for dynamics only
!
! !interface:

      subroutine gfs_dyn_finalize(gc_gfs_dyn, imp_gfs_dyn, exp_gfs_dyn, clock, rc)

!
! !input variables and parameters:
!---------------------------------
      type(esmf_gridcomp)                 :: gc_gfs_dyn
      type(esmf_state)                    :: imp_gfs_dyn, exp_gfs_dyn
      type(esmf_clock)                    :: clock

! !output variables and parameters:
!----------------------------------
      integer,             intent(out)    :: rc

! !working arrays and local parameters.
!--------------------------------------
      type(gfs_dyn_wrap)                            :: wrap   
      type(gfs_dynamics_internal_state), pointer    :: int_state  
      integer                                       :: rc1, rcfinal

!eop
!-------------------------------------------------------------------------

! initialize the error signal variables.
!---------------------------------------
      rc1     = esmf_success
      rcfinal = esmf_success

! retrieve the esmf internal state.
!----------------------------------
     call esmf_logwrite("get internal state in the finalize routine", &
                         ESMF_LOGMSG_INFO, rc=rc1)

     call esmf_gridcompgetinternalstate(gc_gfs_dyn, wrap, rc1)

     call gfs_dynamics_err_msg(rc1,'get internal state in the finalize routine',rc)

! point the local internal state pointer to the esmf internal state pointer.
!------------------------------------------------------------------------------
      int_state => wrap%int_state

! ======================================================================
! run the gfs finalize routine to release the memory space, etc. 
! ======================================================================
      call esmf_logwrite("run gfs_dynamics_finalize", ESMF_LOGMSG_INFO, rc=rc1)

      call gfs_dynamics_finalize(int_state, rc = rc1)

      call gfs_dynamics_err_msg(rc1,'run gfs_dynamics_finalize',rc)
! ======================================================================
! ======================================================================

! print out the final error signal information and put it to rc.
!---------------------------------------------------------------
      call gfs_dynamics_err_msg_final(rcfinal,                          &
                        'finalize from gfs dynamics grid comp.',rc)

      end subroutine gfs_dyn_finalize

! end of the gfs esmf grid component module.
!-------------------------------------------
      end module gfs_dynamics_grid_comp_mod
