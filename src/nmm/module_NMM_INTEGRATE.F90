!-----------------------------------------------------------------------
!
      MODULE module_NMM_INTEGRATE
!
!-----------------------------------------------------------------------
!
!***  This module holds the fundamental NMM integration runstream
!***  when 2-way nesting is being used.
!***  It is called from subroutine NMM_RUN in module_NMM_GRID_COMP.
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
! PROGRAM HISTORY LOG:
!   2008-08     Colon - Moved NMM runstream from ATM_RUN into separate
!                       routines when adding digital filters.
!   2009-07-09  Black - Condense all three NMM integrate routines
!                       into one when when merging with nesting.
!   2010_03_24  Black - Revised for new structure.
!   2010-11-03  Pyle  - Revised for digital filters.
!   2011-02     Yang  - Updated to use both the ESMF 4.0.0rp2 library,
!                       ESMF 5 series library and the the
!                       ESMF 3.1.0rp2 library.
!   2011-05     Yang  - Modified for using the ESMF 5.2.0r_beta_snapshot_07.
!   2011-07     Black - Revised for moving nests.
!   2012-02     Yang  - Modified for using the ESMF 5.2.0rp1 library.
!   2012-07     Black - Revised for 'generational' task usage.
!-----------------------------------------------------------------------
!
      USE ESMF
!
      USE MODULE_ERROR_MSG,ONLY: ERR_MSG,MESSAGE_CHECK
!
      USE module_CLOCKTIMES,ONLY: INTEGRATION_TIMERS                    &
                                 ,PRINT_CLOCKTIMES
!
      USE module_DOMAIN_INTERNAL_STATE,ONLY: DOMAIN_INTERNAL_STATE      &
                                            ,WRAP_DOMAIN_INTERNAL_STATE
!
      USE module_WRITE_ROUTINES,ONLY: WRITE_ASYNC
!
      USE module_NESTING,ONLY: BOUNDARY_DATA_STATE_TO_STATE             &
                              ,INTERIOR_DATA_STATE_TO_STATE
!
      USE module_CONTROL,ONLY: TIMEF
!
      USE module_PARENT_CHILD_CPL_COMP,ONLY: NSTEP_CHILD_RECV
!
      USE module_KINDS
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      PRIVATE
!
      PUBLIC :: NMM_INTEGRATE
!
!-----------------------------------------------------------------------
!
      CHARACTER(ESMF_MAXSTR) :: CWRT                                       !<-- Restart/History label
!
!-----------------------------------------------------------------------
!***  For determining clocktimes of various pieces of the Solver.
!-----------------------------------------------------------------------
!
      REAL(kind=KDBL) :: btim,btim0
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE NMM_INTEGRATE(CLOCK_DIRECTION                          &
                              ,DOMAIN_GRID_COMP                         &
                              ,IMP_STATE_DOMAIN                         &
                              ,EXP_STATE_DOMAIN                         &
                              ,CLOCK_INTEGRATE                          &
                              ,CURRTIME                                 &
                              ,STARTTIME                                &
                              ,TIMESTEP                                 &
                              ,NTIMESTEP_EXT                            &
                              ,RUNSTEPCOUNT                             &
                              ,DT                                       &
                              ,INTEGRATE_TIMESTEP                       &
                              ,INTERVAL_CLOCKTIME                       &
                              ,INTERVAL_HISTORY                         &
                              ,INTERVAL_RESTART                         &
                              ,FILTER_METHOD                            &
                              ,HALFDFIINTVAL                            &
                              ,HALFDFITIME                              &
                              ,NDFISTEP                                 &
                              ,RESTARTED_RUN                            &
                              ,RST_OUT_00                               &
                              ,I_AM_A_FCST_TASK                         &
                              ,I_AM_LEAD_FCST_TASK                      &
                              ,NESTING                                  &
                              ,NEST_MODE                                &
                              ,TASK_MODE                                &
                              ,I_AM_A_NEST                              &
                              ,MY_DOMAIN_ID                             &
                              ,NUM_CHILDREN                             &
                              ,NUM_2WAY_CHILDREN                        &
                              ,PARENT_CHILD_CPL                         &
                              ,IMP_STATE_CPL_NEST                       &
                              ,EXP_STATE_CPL_NEST                       &
                              ,PAR_CHI_TIME_RATIO                       &
                              ,MY_DOMAIN_MOVES                          &
                              ,NTRACK                                   &
                              ,NPHS                                     &
                              ,LAST_GENERATION                          &
                              ,ADVANCED                                 &
                              ,MYPE                                     &
                              ,COMM_GLOBAL                              &
                              ,GENERATION_FINISHED                      &
                              ,TIMERS_DOMAIN                            &
                              ,NPE_PRINT                                &
                              ,PRINT_TIMING )
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!*** The following USEs are needed for NMM-B time series output.
!-----------------------------------------------------------------------
!
      USE MODULE_SOLVER_INTERNAL_STATE, ONLY: SOLVER_INTERNAL_STATE &
                                             ,WRAP_SOLVER_INT_STATE
!
      USE MODULE_TIMESERIES
!
!-----------------
!*** Arguments IN
!-----------------
!
      INTEGER(kind=KINT),INTENT(IN) :: COMM_GLOBAL                      &  !<-- The MPI communicator for all tasks (COMM_WORLD)
                                      ,FILTER_METHOD                    &  !<-- The type of digital filtering desired
                                      ,MYPE                             &  !<-- Local task rank on this domain
                                      ,NPE_PRINT                        &  !<-- Task to print clocktimes
                                      ,NPHS                             &  !<-- Physics timestep
                                      ,NTRACK                           &  !<-- Storm locator flag
                                      ,NUM_2WAY_CHILDREN                &  !<-- How many 2-way children on this domain?
                                      ,RUNSTEPCOUNT                        !<-- # of integration steps per coupling step

      INTEGER(kind=KINT),INTENT(INOUT) :: NUM_CHILDREN                     !<-- # of children on this domain
!
      REAL(kind=KFPT),INTENT(IN) :: DT                                     !<-- Fundamental timestep of this domain (REAL) (s)
!
      LOGICAL(kind=KLOG),INTENT(IN) :: NESTING                          &  !<-- Are there any nested domains?
                                      ,PRINT_TIMING                     &  !<-- Shall we print timing in err file?
                                      ,RESTARTED_RUN                    &  !<-- Is this a restarted run?
                                      ,RST_OUT_00                          !<-- Shall we write 00h history in restarted run?
!
      CHARACTER(5) ,INTENT(IN) :: NEST_MODE                                !<-- 1-way or 2-way nesting
      CHARACTER(8) ,INTENT(IN) :: CLOCK_DIRECTION                          !<-- The direction of time in the Clock
      CHARACTER(12),INTENT(IN) :: TASK_MODE                                !<-- Task assignments unique per domain or generational?

      LOGICAL(kind=KLOG),INTENT(IN) :: I_AM_A_FCST_TASK                 &  !<-- Am I in a forecast task?
                                      ,I_AM_LEAD_FCST_TASK              &  !<-- Am I the first forecast task in the domain's comm?
                                      ,I_AM_A_NEST                         !<-- Am I in a nested domain?
!
      TYPE(ESMF_TimeInterval),INTENT(IN)  :: TIMESTEP                      !<-- Fundamental timestep of this domain (ESMF) (s)
!
!--------------------
!*** Arguments INOUT
!--------------------
!
      INTEGER(kind=KINT),INTENT(INOUT) :: NTIMESTEP_EXT                    !<-- This domain's current timestep count
!
      LOGICAL(kind=KLOG),INTENT(OUT) :: ADVANCED                           !<-- Did the integration advance during the routine?
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: DOMAIN_GRID_COMP                !<-- The DOMAIN component
!
      TYPE(ESMF_Time),INTENT(INOUT) :: CURRTIME                         &  !<-- The clock's current time
                                      ,STARTTIME                           !<-- The clock's start time
!
      TYPE(ESMF_Clock),INTENT(INOUT) :: CLOCK_INTEGRATE                    !<-- This DOMAIN Component's ESMF Clock
!
      TYPE(ESMF_State),INTENT(INOUT) :: IMP_STATE_DOMAIN                &  !<-- Import state of this DOMAIN component 
                                       ,EXP_STATE_DOMAIN                   !<-- Export state of this DOMAIN component
!
      TYPE(INTEGRATION_TIMERS),TARGET :: TIMERS_DOMAIN                     !<-- Clocktime variables to be printed.
!
!------------------------
!***  Optional Arguments
!------------------------
!
      INTEGER(kind=KINT),INTENT(IN),OPTIONAL :: MY_DOMAIN_ID            &  !<-- The ID of this domain 
                                               ,NDFISTEP                &
                                               ,PAR_CHI_TIME_RATIO         !<-- Ratio of parent's timestep to this domain's
!
      LOGICAL(kind=KLOG),INTENT(IN),OPTIONAL :: LAST_GENERATION         &  !<-- For 2-way nests, is this the final generation?
                                               ,MY_DOMAIN_MOVES            !<-- Does my domain move?
!
      LOGICAL(kind=KLOG),INTENT(INOUT),OPTIONAL :: INTEGRATE_TIMESTEP      !<-- For 2-way nests, is this the final generation?
!
      LOGICAL(kind=KLOG),INTENT(OUT),OPTIONAL :: GENERATION_FINISHED       !<-- Is a generation through with its integration?
!
      TYPE(ESMF_State),INTENT(INOUT),OPTIONAL:: IMP_STATE_CPL_NEST      &
                                               ,EXP_STATE_CPL_NEST
!
      TYPE(ESMF_Time),INTENT(IN),OPTIONAL :: HALFDFITIME
!
      TYPE(ESMF_TimeInterval),INTENT(IN),OPTIONAL :: HALFDFIINTVAL      &
                                                    ,INTERVAL_CLOCKTIME &  !<-- Time interval between clocktime prints
                                                    ,INTERVAL_HISTORY   &  !<-- Time interval between history output
                                                    ,INTERVAL_RESTART      !<-- Time interval between restart output
!
      TYPE(ESMF_CplComp),INTENT(INOUT),OPTIONAL :: PARENT_CHILD_CPL        !<-- Coupler component for parent-child/nest exchange
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: YY,YY1,YY2,MM,DD,H,M,S,START_SEC,STOP_SEC
!
      INTEGER(kind=KINT) :: FIRST_STEP,I,KOUNT_STEPS,LAST_STEP          &
                           ,N,NSTEP_INTEGRATE,NSECONDS_FCST             &
                           ,NTIMESTEP       
!
      INTEGER(kind=KINT) :: IERR,RC,RC_INTEG
!
      INTEGER(kind=ESMF_KIND_I8) :: NTIMESTEP_ESMF                         !<-- The current forecast timestep (ESMF_INT)
!
      INTEGER(kind=KINT),DIMENSION(2) :: STORM_CENTER
!
      INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE :: LOC_PAR_CHILD_TIME_RATIO
!
      CHARACTER(2) :: INT_TO_CHAR
      CHARACTER(6) :: FMT
!
      LOGICAL(kind=KLOG) :: ALLCLEAR_FROM_PARENT                        &
                           ,RECV_ALL_CHILD_DATA 
!
      LOGICAL(kind=KLOG) :: E_BDY,N_BDY,S_BDY,W_BDY                     &
                           ,FREE_FORECAST                               &
                           ,FREE_TO_INTEGRATE                           &
                           ,I_AM_ACTIVE                                 &
                           ,INTEGRATED_SOLVER
!
      TYPE(ESMF_Time) :: ALARM_HISTORY_RING                             &
                        ,ALARM_RESTART_RING                             &
                        ,ALARM_CLOCKTIME_RING
!
      TYPE(ESMF_Time) :: STOPTIME
!
      TYPE(ESMF_TimeInterval) :: TIMESTEP_FILTER                           !<-- Dynamics timestep during filter (s) (ESMF)
!
      TYPE(DOMAIN_INTERNAL_STATE),POINTER :: DOMAIN_INT_STATE
!
      TYPE(WRAP_DOMAIN_INTERNAL_STATE) :: WRAP
!
      TYPE(WRAP_SOLVER_INT_STATE) :: WRAP_SOLVER
!
      TYPE(SOLVER_INTERNAL_STATE),POINTER :: SOLVER_INT_STATE
!
      integer(kind=kint),dimension(8) :: values
!-----------------------------------------------------------------------
!***  For timers.
!-----------------------------------------------------------------------
!
      REAL(kind=KFPT) :: phase1_tim,phase3_tim
!
      TYPE(INTEGRATION_TIMERS),POINTER :: TD
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Point timers into this domain's timer object.
!-----------------------------------------------------------------------
!
      TD=>TIMERS_DOMAIN                                                    !<-- Abbreviate the name of this domain's timers.
!
!-----------------------------------------------------------------------
!
      RC      =ESMF_SUCCESS
      RC_INTEG=ESMF_SUCCESS
!
      FMT='(I2.2)'
      WRITE(INT_TO_CHAR,FMT)MY_DOMAIN_ID
!
!-----------------------------------------------------------------------
!***  Before beginning the integration of the DOMAIN component,
!***  extract the internal state which will be needed for
!***  initial writing of history/restart files.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Get DOMAIN Internal State in NMM_INTEGRATE"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompGetInternalState(DOMAIN_GRID_COMP               &
                                        ,WRAP                           &
                                        ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INTEG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      DOMAIN_INT_STATE=>wrap%DOMAIN_INT_STATE
!
!-----------------------------------------------------------------------
!***  When using digital filters the timestep is set back to 0 when
!***  the filtering direction changes and when the filtering ends.
!-----------------------------------------------------------------------
!
      IF(domain_int_state%KOUNT_TIMESTEPS==0)THEN
        domain_int_state%FIRST_PASS=.TRUE.
      ENDIF
!
!
!-----------------------------------------------------------------------
!***  For normal forecast integration set the Alarm ring times
!***  while accounting for restarts and digital filtering.
!-----------------------------------------------------------------------
!
      IF(FILTER_METHOD==0.AND.domain_int_state%FIRST_PASS)THEN
!
!      if(mype==0)then
!        call esmf_clockprint(clock=clock_integrate,rc=rc)
!      endif
        CALL RESET_ALARMS
!
      ENDIF
!
!-----------------------------------------------------------------------
!***  Forecast tasks extract the Solver internal state needed for
!***  timeseries output.  These tasks also need to know if they lie
!***  on a boundary of the current domain.
!-----------------------------------------------------------------------
!
      IF(I_AM_A_FCST_TASK)THEN
!
        CALL ESMF_GridCompGetInternalState(domain_int_state%SOLVER_GRID_COMP &  !<-- The Solver component
                                          ,WRAP_SOLVER                       &  !<-- The F90 wrap of the Solver internal state
                                          ,RC)
!
        SOLVER_INT_STATE => wrap_solver%INT_STATE
!
        S_BDY=(solver_int_state%JTS==solver_int_state%JDS)                 ! This task is on the southern boundary
        N_BDY=(solver_int_state%JTE==solver_int_state%JDE)                 ! This task is on the northern boundary
        W_BDY=(solver_int_state%ITS==solver_int_state%IDS)                 ! This task is on the western boundary
        E_BDY=(solver_int_state%ITE==solver_int_state%IDE)                 ! This task is on the eastern boundary
!
      END IF
!
!-----------------------------------------------------------------------
!***  The limits on the integration timeloop below depend on the
!***  mode of nesting.  For 1-way nests the domains all integrate 
!***  straight through from the start to the end of the forecast.
!***  Two-way nesting is different due to the children's feedback.
!***  All generations except the lowermost will execute only one
!***  timestep at a time then return since their domains cannot 
!***  proceed until they receive internal updates from their children.
!***  The domains in the lowermost generation have no children and
!***  can thus execute a full N timesteps at a time where N is the
!***  number of timesteps within a single timestep of their parents.
!***  Since each task saves the final timestep in the internal 
!***  state of each domain it is on then it needs to determine the
!***  forecast's final timestep only once.
!-----------------------------------------------------------------------
!
      IF(domain_int_state%FIRST_PASS)THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!xxx    MESSAGE_CHECK="NMM_INTEGRATE: Get the Final Timestep"   
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!xxx    CALL ESMF_ClockGet(clock           =CLOCK_INTEGRATE                 &
!xxx                      ,runTimeStepCount=domain_int_state%TIMESTEP_FINAL &  !<-- Final timestep of this domain's forecast
!xxx                      ,rc              =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!xxx    CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INTEG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        domain_int_state%KOUNT_TIMESTEPS=0                                 !<-- Initialize timestep counter on this domain
!       domain_int_state%FIRST_PASS=.FALSE.
!
      ENDIF
!
!-----------------------------------------------------------------------
!
      FIRST_STEP=NTIMESTEP_EXT
!
!-----------------------------------------------------------------------
!***  What is the last timestep of the entire integration?
!-----------------------------------------------------------------------
!
      IF(TASK_MODE=='unique')THEN                                          !<-- Single domains and 1-way nests integrate to end of fcst
!
!xxx fix for coupled runs
!!!!!!  cpl_fix: if(runstepcount>0)then  !<-- this means there is coupling
!!!!!!    last_step=first_step+runstepcount-1
!!!!!!  else cpl_fix
!xxxxfix for coupled runs  
!!!!!!  IF(FILTER_METHOD==0.AND.domain_int_state%FIRST_PASS)THEN           !<-- Free forecast
        IF(FILTER_METHOD==0)THEN                                           !<-- Free forecast
!
!***  Without this ClockGet the FIRST_PASS from above represents the filter
!***  clock if a filter case.
!
          CALL ESMF_ClockGet(clock           =CLOCK_INTEGRATE                 &
                            ,runTimeStepCount=domain_int_state%TIMESTEP_FINAL &
                            ,stopTime        =STOPTIME                        &
                            ,rc              =RC)
!
          IF(domain_int_state%FIRST_PASS)THEN
            CALL RESET_ALARMS
          ENDIF
!
!         LAST_STEP=NINT(domain_int_state%TIMESTEP_FINAL)-1
          LAST_STEP=FIRST_STEP+RUNSTEPCOUNT-1
!     if(mype==0)then
!       write(0,33661)last_step,first_step,runstepcount
33661   format(' NMM_INTEG last_step=',i5,' first_step=',i5,' runtimestepcount=',i5)
!     endif
!
        ELSE                                                               !<-- Digital filter
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="NMM_INTEGRATE: Get Filter StopTime"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_ClockGet(clock   =CLOCK_INTEGRATE                   &  !<-- The filter clock
                            ,stopTime=STOPTIME                          &  !<-- The simulation stop time (ESMF)
                            ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INTEG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_TimeGet(STARTTIME,yy=YY1,s=START_SEC,rc=RC)
      if(rc/=0)then
        write(0,20001)rc
20001   format(' NMM_INTEGRATE timeget1 rc=',i3)
      endif
          CALL ESMF_TimeGet(STOPTIME ,yy=YY2,s=STOP_SEC ,rc=RC)
      if(rc/=0)then
        write(0,20002)rc
20002   format(' NMM_INTEGRATE timeget2 rc=',i3)
      endif
!
          IF(YY1/=YY2)THEN                                                !<-- Account for dates that straddle 31 Dec / 1 Jan
!
            IF(CLOCK_DIRECTION=='Bckward ')THEN
              IF(MOD(YY2,4)==0)THEN
                STOP_SEC=STOP_SEC-31557600
              ELSE
                STOP_SEC=STOP_SEC-31536000
              ENDIF
!
            ELSEIF(CLOCK_DIRECTION=='Forward ')THEN
              IF(MOD(YY1,4)==0)THEN
                START_SEC=START_SEC-31557600
              ELSE
                START_SEC=START_SEC-31536000
              ENDIF
!
            ENDIF
!
          ENDIF
!
          NSECONDS_FCST=ABS(STOP_SEC-START_SEC)                           !<-- The forecast length (sec) (REAL)
          LAST_STEP=NINT(NSECONDS_FCST/DT)-1 
 
        ENDIF
!
      ELSEIF(TASK_MODE=='generational')THEN
        IF(.NOT.PRESENT(LAST_GENERATION))THEN
          WRITE(0,*)' LAST_GENERATION must be supplied for 2-way nesting but was not.'
          WRITE(0,*)' Aborting!!!'
          CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
        ENDIF
        IF(.NOT.LAST_GENERATION)THEN
          LAST_STEP=FIRST_STEP                                             !<-- Most 2-way domains integrate one timestep at a time
        ELSE
          LAST_STEP=NTIMESTEP_EXT+PAR_CHI_TIME_RATIO-1                     !<-- Lowest generation 2-way nests run through a parent timestep
        ENDIF
!
      ENDIF
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  THE INTEGRATION TIME LOOP OF THE ATMOSPHERE.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  For runs with a single domain or with 1-way nesting then the
!***  tasks will execute straight through from the beginning of the
!***  forecast (or the restart time) to the end of the forecast.
!***  For runs with 2-way nesting the tasks integrate only through
!***  one or a few timesteps at a time as described above.  However note that
!***  the timestep will be incremented only after phase 1 of the
!***  Domain component's Run step is executed.  Due to the nature of
!***  the generational use of task assignments in 2-way nesting some
!***  tasks will enter timeloop_drv but not be allowed to integrate
!***  (i.e., to call DOMAIN_RUN which is phase 1 of the Run step of
!***  the Domain component) because: (1) parent domains must first
!***  recv 2-way exchange data from all of their children at the end
!***  of each parent timestep; (2) child domains at the end of their
!***  parents' timesteps must be informed by their parent that the
!***  parent did recv exchange data from all its children meaning the
!***  given child can proceed because its parent is free to integrate
!***  to the end of its next timestep and send back BC update data.
!***  If a domain is both a parent and a child then both of those
!***  conditions must be true for the domain to integrate another
!***  timestep. 
!-----------------------------------------------------------------------
!     if(mype==0)then
!       write(0,44251)ntimestep_ext,first_step,last_step
44251   format(' NMM_INTEG before timeloop_drv ntimestep_ext=',i4,' first_step=',i4,' last_step=',i4)                                  
!     endif
!
      timeloop_drv: DO NSTEP_INTEGRATE=FIRST_STEP,LAST_STEP                !<-- For 1-way nesting this would go from start to end of the fcst
!
!     call print_memory()
!-----------------------------------------------------------------------
!
        KOUNT_STEPS=domain_int_state%KOUNT_TIMESTEPS
        NTIMESTEP=KOUNT_STEPS                                              !<-- Internal timesteps to avoid confusion with NUOPC Clock
        I_AM_ACTIVE=.TRUE.
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="NMM_INTEGRATE: Extract Free Forecast flag"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state=IMP_STATE_CPL_NEST                 &  !<-- The parent-child coupler import state
                              ,name ='Free Forecast'                    &  !<-- Flag for free forecasts.
                              ,value=FREE_FORECAST                      &  !<-- Is this the freee forecast?
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INTEG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        IF(I_AM_A_FCST_TASK.AND..NOT.FREE_FORECAST)THEN
          IF(MY_DOMAIN_ID>1.OR.NUM_CHILDREN>1)THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="NMM_INTEGRATE: Extract I_AM_ACTIVE"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_AttributeGet(state=EXP_STATE_CPL_NEST             &  !<-- The parent-child coupler export state
                                  ,name ='I Am Active'                  &  !<-- Flag for digital filter activity.
                                  ,value=I_AM_ACTIVE                    &  !<-- Is this domain active in the digital filter?
                                  ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INTEG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          ENDIF
        ENDIF
!
!-----------------------------------------------------------------------
!***  For 2-way nesting check for the signals from parents and
!***  children to know if the execution can proceed into this timestep.
!***  Call phase 1 of the Parent-Child coupler's Run step to perform
!***  these checks.  The subroutine name is CHECK_2WAY_SIGNALS.
!-----------------------------------------------------------------------
!
      check_2way: IF(NEST_MODE=='2-way'                                 &
                         .AND.                                          &
                     I_AM_A_FCST_TASK)THEN     
      btim0=timef()
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Call Phase 1 Coupler Run: Check 2-Way Signals"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_CplCompRun(cplcomp    =PARENT_CHILD_CPL               &  !<-- The Parent-Child coupler component
                            ,importState=IMP_STATE_CPL_NEST             &  !<-- The Parent-Child coupler import state
                            ,exportState=EXP_STATE_CPL_NEST             &  !<-- The Parent-Child coupler export state
                            ,clock      =CLOCK_INTEGRATE                &  !<-- The Domain Clock
                            ,phase      =1                              &  !<-- The phase (subroutine) of the coupler to execute
                            ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INTEG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  What are the values of the 2-way signals?
!-----------------------------------------------------------------------
!
        IF(KOUNT_STEPS>0)THEN
!
!--------------------------------------------------------------
!***  Is 2-way data ready for my parent from all its children?
!--------------------------------------------------------------
!
          IF(I_AM_A_NEST)THEN
            IF(MOD(NSTEP_INTEGRATE,PAR_CHI_TIME_RATIO)==0)THEN
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="NMM_INTEGRATE: Extract ALLCLEAR from P-C Exp State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
              CALL ESMF_AttributeGet(state=EXP_STATE_CPL_NEST           &  !<-- The parent-child coupler export state
                                    ,name ='ALLCLEAR'                   &  !<-- Flag for proceeding now in timestep
                                    ,value=ALLCLEAR_FROM_PARENT         &  !<-- Parent did/not recv all exch data; child can/not proceed
                                    ,rc   =RC)
!
              domain_int_state%ALLCLEAR_FROM_PARENT=ALLCLEAR_FROM_PARENT
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INTEG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
              IF(.NOT.domain_int_state%ALLCLEAR_FROM_PARENT)THEN
!
                RETURN                                                     !<-- All my siblings are not yet ready to send our parent
!                                                                          !    their 2-way data.
              ENDIF
            ENDIF
          ENDIF
!
!------------------------------------------------------------------
!***  Are all of my children ready to send their 2-way data to me?
!------------------------------------------------------------------
!
          IF(NUM_2WAY_CHILDREN>0)THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Call Coupler: Extract RECV_ALL_CHILD_DATA from P-C Exp State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_AttributeGet(state=EXP_STATE_CPL_NEST             &  !<-- The parent-child coupler export state
                                  ,name ='Recv All Child Data'          &  !<-- Flag for integration (true or false)
                                  ,value=RECV_ALL_CHILD_DATA            &  !<-- The value of the integration flag
                                  ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INTEG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            domain_int_state%RECV_ALL_CHILD_DATA=RECV_ALL_CHILD_DATA
!
            IF(.NOT.domain_int_state%RECV_ALL_CHILD_DATA)THEN
!
              RETURN                                                       !<-- All my children are not yet ready to send me
!                                                                          !    their 2-way data.
            ENDIF
          ENDIF
!
!-----------------------------------------------------------------------
!
        ENDIF
!
!-----------------------------------------------------------------------
!
        td%pc_cpl_run_cpl1=td%pc_cpl_run_cpl1+(timef()-btim0)
      ENDIF check_2way
!
!-----------------------------------------------------------------------
!***   Children receive data from their parents at the beginning
!***   of child timesteps that coincide with the beginning of
!***   parent timesteps.
!
!***   (1) At the beginning of every such timestep the children
!***       receive new boundary data that is sent by their parents
!***       from the end of that parent timestep so the children can
!***       compute boundary value tendencies to be used for the
!***       integration through the next N child timesteps where
!***       N is the number of child timesteps per parent timestep.
!***       The handling of these new boundary values is the only
!***       action in this phase of the Parent-Child coupler for
!***       domains that are static 1-way nests.
!
!***   (2) If a child is a moving nest that has decided it must move:
!***       (a) It sends a message to its parent informing it of that
!***           fact along with the location to which it is moving on
!***           the parent grid.  That move must happen at a
!***           parent timestep in the future because the parent
!***           must provide data to some internal child points
!***           following a move and since the parent must always
!***           run ahead of its children (to provide the boundary
!***           data from the future) then the new data for those
!***           internal child points must also originate at a
!***           future timestep.  That future parent timestep in
!***           which the nest will shift is also sent to the parent.
!***       (b) If the current timestep is equal to the timestep
!***           in which the nest determined it wants to shift then
!***           the child now Recvs the parent data for its internal
!***           gridpoints that have moved over a new portion of the
!***           parent grid as well as the first boundary data at the
!***           new location.
!***       (c) The child Recvs the new boundary data sent by the parent
!***           from the future [see (1) above].
!-----------------------------------------------------------------------
!
!
        IF(I_AM_A_NEST.AND.I_AM_A_FCST_TASK)THEN
          btim0=timef()
!
          IF(MOD(KOUNT_STEPS,PAR_CHI_TIME_RATIO)==0)THEN                   !<-- Child is at the start of a parent timestep.
!
!-----------------------------------------------------------------------
!***  Call Phase 2 of the Run step of the Parent-Child coupler in
!***  which children receive BC data from their parents.  The
!***  name of the subroutine is CHILDREN_RECV_PARENT_DATA.  If the
!***  child is a moving nest then its shifts in position take place
!***  in this phase as well since BC data from the parent depends 
!***  on the location of the child domain.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Call Phase 2 Coupler Run: Children Recv Data from Parents"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_CplCompRun(cplcomp    =PARENT_CHILD_CPL           &  !<-- The Parent-Child coupler component
                                ,importState=IMP_STATE_CPL_NEST         &  !<-- The Parent-Child coupler import state
                                ,exportState=EXP_STATE_CPL_NEST         &  !<-- The Parent-Child coupler export state
                                ,clock      =CLOCK_INTEGRATE            &  !<-- The Domain Clock
                                ,phase      =2                          &  !<-- The phase (subroutine) of the coupler to execute
                                ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INTEG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  The new nest boundary data must be moved from the Parent-Child
!***  coupler into the nest's DOMAIN component.  Do not do this
!***  if the digital filter is running and the child domain is not
!***  active.
!-----------------------------------------------------------------------
!
            IF(FREE_FORECAST.OR.(FILTER_METHOD>0.AND.I_AM_ACTIVE))THEN
              CALL BOUNDARY_DATA_STATE_TO_STATE(s_bdy    =S_BDY                         &  !<-- This task lies on a south boundary?
                                               ,n_bdy    =N_BDY                         &  !<-- This task lies on a north boundary?
                                               ,w_bdy    =W_BDY                         &  !<-- This task lies on a west boundary?
                                               ,e_bdy    =E_BDY                         &  !<-- This task lies on an east boundary?
                                               ,nest     =domain_int_state%I_AM_A_NEST  &  !<-- The nest flag (yes or no)
                                               ,state_in =EXP_STATE_CPL_NEST            &  !<-- The P-C coupler export state
                                               ,state_out=IMP_STATE_DOMAIN)                !<-- The Domain import state
            ENDIF
!
!-----------------------------------------------------------------------
!***  If the nest is movable then the DOMAIN component must be 
!***  informed if the nest does or does not want to move now.
!***  If it does want to move now then all of the interior update
!***  data generated by the parent must be transferred to the DOMAIN
!***  import state.
!-----------------------------------------------------------------------
!
            IF(MY_DOMAIN_MOVES)THEN
!
!
              CALL INTERIOR_DATA_STATE_TO_STATE(EXP_STATE_CPL_NEST      &
                                               ,IMP_STATE_DOMAIN )
!
            ENDIF
!
!-----------------------------------------------------------------------
!
          ENDIF
!
          td%pc_cpl_run_cpl2=td%pc_cpl_run_cpl2+(timef()-btim0)
        ENDIF
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!***  Call phase 3 of the Run step of the Parent-Child coupler.
!***  The name of the subroutine is PARENTS_RECV_CHILD_2WAY_DATA.
!***  It is at this point that the appropriate parent tasks receive
!***  2-way exchange data from their children.
!-----------------------------------------------------------------------
!
!
        IF(I_AM_A_FCST_TASK)THEN
!
          IF(NUM_2WAY_CHILDREN>0)THEN                                      !<-- Parents w/ 2way children call phase 3 of P-C coupler
!
            IF(KOUNT_STEPS>0)THEN
              btim0=timef()
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Call Phase 3 Coupler Run: Parents Recv Exchange Data from Children"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
              CALL ESMF_CplCompRun(cplcomp    =PARENT_CHILD_CPL         &  !<-- The parent-child coupler component
                                  ,importState=IMP_STATE_CPL_NEST       &  !<-- The parent-child coupler import state
                                  ,exportState=EXP_STATE_CPL_NEST       &  !<-- The parent-child coupler export state
                                  ,clock      =CLOCK_INTEGRATE          &  !<-- The DOMAIN Clock
                                  ,phase      =3                        &  !<-- The phase (subroutine) of the coupler to execute
                                  ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INTEG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
              td%pc_cpl_run_cpl3=td%pc_cpl_run_cpl3+(timef()-btim0)
            ENDIF
!
          ENDIF
!
        ENDIF
!
!
!-----------------------------------------------------------------------
!***  If filtering is not in effect and this is the start or restart
!***  of a forecast then write out a history file.
!-----------------------------------------------------------------------
!
        history_output_0_a: IF(NTIMESTEP==0                             &
                                 .AND.                                  &
                              .NOT.RESTARTED_RUN                        &
                                 .AND.                                  &
                               FILTER_METHOD==0                         &
                                 .AND.                                  &
                               .NOT.domain_int_state%WROTE_1ST_HIST     &
                                 .AND.                                  &
                               domain_int_state%QUILTING)THEN
!
          CWRT='History'
          domain_int_state%WROTE_1ST_HIST=.TRUE.
!
          CALL WRITE_ASYNC(DOMAIN_GRID_COMP                             &
                          ,DOMAIN_INT_STATE                             &
                          ,CLOCK_INTEGRATE                              &
                          ,MYPE                                         &
                          ,CWRT)
!
        ENDIF  history_output_0_a
!
        history_output_0_b: IF(RESTARTED_RUN                            &
                                 .AND.                                  &
                               domain_int_state%RESTARTED_RUN_FIRST     &
                                 .AND.                                  &
                               RST_OUT_00                               &
                                 .AND.                                  &
                               .NOT.domain_int_state%WROTE_1ST_HIST     &
                                 .AND.                                  &
                               (FILTER_METHOD == 0)                     &
                                 .AND.                                  &
                               domain_int_state%QUILTING)THEN
!
          domain_int_state%RESTARTED_RUN_FIRST=.FALSE.
          CWRT='History'
          domain_int_state%WROTE_1ST_HIST=.TRUE.
!
          CALL WRITE_ASYNC(DOMAIN_GRID_COMP                             &
                          ,DOMAIN_INT_STATE                             &
                          ,CLOCK_INTEGRATE                              &
                          ,MYPE                                         &
                          ,CWRT)
!
        ENDIF  history_output_0_b
!
!-----------------------------------------------------------------------
!***  Initialize the timeseries output and write timestep 0 data
!***  for this domain.
!-----------------------------------------------------------------------
!
        time_series_0: IF(.NOT.domain_int_state%TS_INITIALIZED) THEN
!
          IF(MYPE<domain_int_state%NUM_PES_FCST)THEN
!
            CALL TIMESERIES_INITIALIZE(SOLVER_INT_STATE                 &
                                      ,NTIMESTEP                        &
                                      ,IERR)
!
            IF (IERR /= 0) THEN
              CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
            END IF
!
            CALL TIMESERIES_RUN(SOLVER_INT_STATE                        &
                               ,NTIMESTEP                               &
                               ,IERR)
!
            IF (IERR /= 0) THEN
              CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
            END IF
!
          END IF
!
          domain_int_state%TS_INITIALIZED = .TRUE.
!
        END IF time_series_0
!
!-----------------------------------------------------------------------
!***  We are now ready to execute phase 1 of the Domain component's
!***  Run step.  This is where the forecast integration takes place.
!-----------------------------------------------------------------------
!
        btim0=timef()
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="NMM_INTEGRATE: Run DOMAIN Component "//INT_TO_CHAR
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!       IF(I_AM_A_FCST_TASK)THEN
          IF(FREE_FORECAST.OR.(FILTER_METHOD>0.AND.I_AM_ACTIVE))THEN
            CALL ESMF_GridCompRun(gridcomp   =DOMAIN_GRID_COMP          &  !<-- The DOMAIN gridded component
                                 ,importState=IMP_STATE_DOMAIN          &  !<-- The DOMAIN import state
                                 ,exportState=EXP_STATE_DOMAIN          &  !<-- The DOMAIN export state
                                 ,clock      =CLOCK_INTEGRATE           &  !<-- The ESMF DOMAIN Clock
                                 ,phase      =1                         &  !<-- The phase (subroutine) of DOMAIN Run to execute
                                 ,rc         =RC)
          ENDIF
!       ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INTEG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        phase1_tim = (timef()-btim0)
        td%domain_run_1=td%domain_run_1+phase1_tim
!
!-----------------------------------------------------------------------
!***  If this domain moves then transfer the storm center location
!***  to the P-C coupler import state.
!-----------------------------------------------------------------------
!
        IF(I_AM_A_FCST_TASK.AND.MY_DOMAIN_MOVES.AND.NTRACK>0)THEN
!
          IF(NTIMESTEP==0.                                              &
                 .OR.                                                   &
             MOD(NTIMESTEP+1,NTRACK*NPHS)==0)THEN

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="NMM_INTEGRATE: Get storm center location."
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_AttributeGet(state    =EXP_STATE_DOMAIN           &  !<-- The Domain component export state
                                  ,name     ='Storm Center'             &  !<-- Name of the attribute to extract
                                  ,valueList=STORM_CENTER               &  !<-- I,J of storm center
                                  ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INTEG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="NMM_Integrate: Transfer storm center location."
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_AttributeSet(state    =IMP_STATE_CPL_NEST         &  !<-- The parent-child coupler import state
                                  ,name     ='Storm Center'             &  !<-- I,J of storm center
                                  ,itemCount=2                          &  !<-- There are 2 words
                                  ,valueList=STORM_CENTER               &  !<-- The data is here.
                                  ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INTEG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          ENDIF
!
        ENDIF
!
!-----------------------------------------------------------------------
!***  Now that this domain has integrated a timestep, if it is a 
!***  2-way run then reset the 2-way signal flags from the parent
!***  and children of this domain.
!-----------------------------------------------------------------------
!
        IF(NEST_MODE=='2-way')THEN                                         !<-- Reset these 2-way flags
          IF(I_AM_A_NEST.AND.I_AM_A_FCST_TASK)THEN
            domain_int_state%RECV_ALL_CHILD_DATA=.FALSE. 
!
            IF(MOD(KOUNT_STEPS+1,PAR_CHI_TIME_RATIO)==0)THEN
              domain_int_state%ALLCLEAR_FROM_PARENT=.FALSE.
              ALLCLEAR_FROM_PARENT=domain_int_state%ALLCLEAR_FROM_PARENT
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="NMM_Integrate: Set ALLCLEAR Signal in P-C Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
              CALL ESMF_AttributeSet(state=EXP_STATE_CPL_NEST           &  !<-- The parent-child coupler import state
                                    ,name ='ALLCLEAR'                   &  !<-- Flag for ALLCLEAR from parent (reset to false) 
                                    ,value=ALLCLEAR_FROM_PARENT         &  !<-- The value of the ALLCLEAR flag
                                    ,rc   =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INTEG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            ENDIF
!  
          ENDIF
!
        ENDIF
!
!-----------------------------------------------------------------------
!***  A parent sends BC data to its children at the end of each parent
!***  timestep and before any potential filter averaging.
!***
!***  (1) For static nests the parents compute only once the association
!***      between their tasks and their children's boundary tasks then
!***      send those child tasks the information they need in order to
!***      be able to properly receive forecast data.  Then the parents
!***      send the new boundary data to the child boundary tasks.
!***  (2) For moving nests the parents are sent the new location of
!***      any of their children who moved.  The parents then recompute
!***      the association between their tasks and their children's
!***      boundary tasks as well as with child tasks in the new
!***      region of the parent into which the children moved.  Parents
!***      send the pertinent child tasks information they need in
!***      order to receive BC data.  Finally the parents send their
!***      moving children new boundary data plus new internal data 
!***      for the new area of the parent newly covered by the most
!***      recent motion of the nests.
!-----------------------------------------------------------------------
!
!
        IF(NUM_CHILDREN>0.AND.I_AM_A_FCST_TASK)THEN                        !<-- Fcst tasks call the coupler if there are children
          btim0=timef()
!
!-----------------------------------------------------------------------
!***  Call the Run step for phase 4 of the Parent-Child coupler.
!***  The name of the subroutine is PARENTS_SEND_CHILD_DATA.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Call Phase 4 Coupler Run: Parents Send Child Data"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_CplCompRun(cplcomp    =PARENT_CHILD_CPL             &  !<-- The parent-child coupler component
                              ,importState=IMP_STATE_CPL_NEST           &  !<-- The parent-child coupler import state
                              ,exportState=EXP_STATE_CPL_NEST           &  !<-- The parent-child coupler export state
                              ,clock      =CLOCK_INTEGRATE              &  !<-- The DOMAIN Clock
                              ,phase      =4                            &  !<-- The phase (subroutine) of the coupler to execute
                              ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INTEG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL BOUNDARY_DATA_STATE_TO_STATE(parent   =domain_int_state%I_AM_A_PARENT  &  !<-- Is this a parent domain?
                                           ,state_in =             EXP_STATE_CPL_NEST &  !<-- The P-C coupler export state
                                           ,state_out=             IMP_STATE_DOMAIN)     !<-- The Domain import state
!
          td%pc_cpl_run_cpl4=td%pc_cpl_run_cpl4+(timef()-btim0)
        ENDIF
!
!
!-----------------------------------------------------------------------
!***  Call the Run step for phase 5 of the Parent-Child coupler.
!***  The name of the subroutine is CHILDREN_SEND_PARENTS_2WAY_DATA.
!***  For 2-way nesting the children send exchange data to their
!***  parents at the end of each parent timestep.
!-----------------------------------------------------------------------
!
!
        IF(I_AM_A_NEST.AND.I_AM_A_FCST_TASK)THEN
!
          IF(NEST_MODE=='2-way')THEN
!
            IF(MOD(KOUNT_STEPS+1,PAR_CHI_TIME_RATIO)==0                 &  !<-- If true then this child has
                           .AND.                                        &  !    reached the end of a timestep
               MY_DOMAIN_ID>1)THEN                                         !    of its parent.
               btim0=timef()
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
              MESSAGE_CHECK="Call Phase 5 Coupler Run: Children Send 2-Way Data to Parents"
!             CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
              CALL ESMF_CplCompRun(            PARENT_CHILD_CPL         &  !<-- The Parent-Child coupler component
                                  ,importState=IMP_STATE_CPL_NEST       &  !<-- The Parent-Child coupler import state
                                  ,exportState=EXP_STATE_CPL_NEST       &  !<-- The Parent-Child coupler export state
                                  ,clock      =CLOCK_INTEGRATE          &  !<-- The Domain Clock
                                  ,phase      =5                        &  !<-- The phase (subroutine) of the coupler to execute
                                  ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
              CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INTEG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
              td%pc_cpl_run_cpl5=td%pc_cpl_run_cpl5+(timef()-btim0)
            ENDIF
!
          ENDIF
!
        ENDIF
!
!
!-----------------------------------------------------------------------
!***  If digital filtering is currently executing then call phase 2
!***  of the Domain component's Run step.  The subroutine name of this
!***  phase is NMM_FILTERING.
!-----------------------------------------------------------------------
!
        IF(FILTER_METHOD>0.AND.I_AM_A_FCST_TASK.AND.I_AM_ACTIVE)THEN
          btim0=timef()
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="NMM_INTEGRATE: Phase 2 Domain Run for Filtering "
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_GridCompRun(gridcomp   =DOMAIN_GRID_COMP            &  !<-- The DOMAIN gridded component for this domain
                               ,importState=IMP_STATE_DOMAIN            &  !<-- The DOMAIN import state
                               ,exportState=EXP_STATE_DOMAIN            &  !<-- The DOMAIN export state
                               ,clock      =CLOCK_INTEGRATE             &  !<-- The ESMF DOMAIN Clock
                               ,phase      =2                           &  !<-- The phase (subroutine) of DOMAIN Run to execute
                               ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INTEG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          td%domain_run_2=td%domain_run_2+(timef()-btim0)
        ENDIF
!
!
!-----------------------------------------------------------------------
!***  Increment the timestep if integration took place.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="NMM_INTEGRATE: Advance the Timestep"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ClockAdvance(clock=CLOCK_INTEGRATE                    &
                              ,rc   =RC)
!
        KOUNT_STEPS=KOUNT_STEPS+1
        domain_int_state%KOUNT_TIMESTEPS=KOUNT_STEPS
        ADVANCED=.TRUE.
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INTEG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        IF(FILTER_METHOD > 0 .AND. I_AM_LEAD_FCST_TASK ) THEN
          WRITE(0,*)'Filter is running, KOUNT_STEPS= ',KOUNT_STEPS,' for method=',FILTER_METHOD
        ENDIF
!
!-----------------------------------------------------------------------
!***  Retrieve the timestep from the DOMAIN clock and
!***  print the forecast time.
!-----------------------------------------------------------------------
!
        CALL ESMF_ClockGet(clock       =CLOCK_INTEGRATE                 &
                          ,advanceCount=NTIMESTEP_ESMF                  &  !<-- # of times the clock has advanced
                          ,rc          =RC)
!
        NTIMESTEP=NTIMESTEP_ESMF
!     write(0,57611)ntimestep,kount_steps
57611 format(' NMM_INTEGRATE advanced CLOCK_INTEGRATE ntimestep=',i5 &
            ,' kount_steps=',i5)
!
!-----------------------------------------------------------------------
!***  Write timeseries data for this timestep on this domain.
!-----------------------------------------------------------------------
!
        IF(MYPE<domain_int_state%NUM_PES_FCST.AND.FILTER_METHOD==0)THEN
!
          CALL TIMESERIES_RUN(SOLVER_INT_STATE                          &
                             ,NTIMESTEP                                 &
                             ,IERR)
!
          IF (IERR /= 0) THEN
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
          END IF
!
        END IF
!
!-----------------------------------------------------------------------
!***  Now that the clock has been advanced, write the history output
!***  if it is time to do so.  This must be done through the DOMAIN
!***  component since its internal state contains the output data 
!***  so we call phase 3 of the Run step for the DOMAIN components.
!***  The subroutine name of this phase is CALL_WRITE_ASYNC.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="NMM_INTEGRATE: Phase 3 Domain Run for History Output"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!
        IF(FILTER_METHOD==0)THEN
          btim0=timef()
!
          CALL ESMF_GridCompRun(gridcomp   =DOMAIN_GRID_COMP            &  !<-- The DOMAIN gridded component
                               ,importState=IMP_STATE_DOMAIN            &  !<-- The DOMAIN import state
                               ,exportState=EXP_STATE_DOMAIN            &  !<-- The DOMAIN export state
                               ,clock      =CLOCK_INTEGRATE             &  !<-- The ESMF Clock for "mini" forecast
                               ,phase      =3                           &  !<-- The phase (subroutine) of DOMAIN Run to execute
                               ,rc         =RC)
          td%domain_run_3=td%domain_run_3+(timef()-btim0)
!
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INTEG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        ENDIF
!
!-----------------------------------------------------------------------
!***  Lead forecast task prints timestep information in free forecast.
!-----------------------------------------------------------------------
!
        IF(I_AM_LEAD_FCST_TASK.AND.FILTER_METHOD==0)THEN
          WRITE(0,25)NTIMESTEP-1,MY_DOMAIN_ID,NTIMESTEP*DT/3600.,phase1_tim
   25     FORMAT(' Finished Timestep ',i6,' for domain ',i3,' ending at ' &
                 ,f7.3,' hours: elapsed integration time ',f9.5)
        ENDIF
!
!-----------------------------------------------------------------------
!***  Print clocktimes of integration sections on the MPI task that
!***  was specified in this domain's configure file (npe_print).
!-----------------------------------------------------------------------
!
        IF(FILTER_METHOD==0)THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Filter Method=0  Is ALARM_CLOCKTIME ringing?"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          IF(ESMF_AlarmIsRinging(alarm=domain_int_state%ALARM_CLOCKTIME &  !<-- The alarm to print clocktimes used by model parts
                                ,rc   =RC))THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INTEG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            IF(PRINT_TIMING)THEN
              CALL PRINT_CLOCKTIMES(NTIMESTEP                           &
                                   ,MY_DOMAIN_ID                        &
                                   ,MYPE                                &
                                   ,NPE_PRINT                           &
                                   ,TIMERS_DOMAIN )
            ENDIF
!
          ENDIF
!
        ENDIF
!
!-----------------------------------------------------------------------
!
      IF(domain_int_state%FIRST_PASS)THEN
        domain_int_state%FIRST_PASS=.FALSE.
      ENDIF
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!!!!  call print_memory()
!
      ENDDO timeloop_drv
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  If the execution is not at the end of the forecast then RETURN.
!-----------------------------------------------------------------------
!
      IF(.NOT.ESMF_ClockIsStopTime(CLOCK_INTEGRATE,rc=RC))THEN
!
        RETURN
!
      ELSE
!
        IF(PRESENT(GENERATION_FINISHED))THEN
          GENERATION_FINISHED=.TRUE.
        ENDIF
!
      ENDIF
!
!-----------------------------------------------------------------------
!
      driver_run_end: IF(FILTER_METHOD==0)THEN                             !<-- For standard integration, no filtering
!
!-----------------------------------------------------------------------
!***  Extract Clocktimes of the Parent-Child Coupler from that
!***  component's export state and print them.
!-----------------------------------------------------------------------
!
        IF(I_AM_A_NEST.AND.I_AM_A_FCST_TASK)THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Extract Cpl1 Recv Time from Parent-Child Cpl Export State"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeGet(state=EXP_STATE_CPL_NEST               &  !<-- The Parent-Child Coupler export state
                                ,name ='Cpl1_Recv_Time'                 &  !<-- Name of the attribute to extract
                                ,value=td%cpl1_recv_tim                 &  !<-- Clocktime for Recv in phase 1 of Cpl Init
                                ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INTEG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        ENDIF
!
        IF(NUM_CHILDREN>0.AND.I_AM_A_FCST_TASK)THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Extract Cpl2 Wait Time from Parent-Child Cpl Export State"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeGet(state=EXP_STATE_CPL_NEST               &  !<-- The Parent-Child Coupler export state
                                ,name ='Cpl2_Wait_Time'                 &  !<-- Name of the attribute to extract
                                ,value=td%cpl2_wait_tim                 &  !<-- Clocktime for Wait in phase 2 of Cpl
                                ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INTEG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Extract Cpl2 Comp Time from Parent-Child Cpl Export State"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeGet(state=EXP_STATE_CPL_NEST               &  !<-- The Parent-Child Coupler export state
                                ,name ='Cpl2_Comp_Time'                 &  !<-- Name of the attribute to extract
                                ,value=td%cpl2_comp_tim                 &  !<-- Clocktime for Compute in Phase 2 of Cpl
                                ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INTEG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Extract Cpl2 Send Time from Parent-Child Cpl Export State"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeGet(state=EXP_STATE_CPL_NEST               &  !<-- The Parent-Child Coupler export state
                                ,name ='Cpl2_Send_Time'                 &  !<-- Name of the attribute to extract
                                ,value=td%cpl2_send_tim                 &  !<-- Clocktime for Send in Phase 2 of Cpl
                                ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INTEG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          MESSAGE_CHECK="Extract parent_bookkeep_moving_tim from Parent-Child Cpl Export State"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeGet(state=EXP_STATE_CPL_NEST               &  !<-- The Parent-Child Coupler export state
                                ,name ='parent_bookkeep_moving_tim'     &  !<-- Name of the attribute to extract
                                ,value=td%parent_bookkeep_moving_tim    &  !<-- moving nest bookeeping time
                                ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INTEG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          MESSAGE_CHECK="Extract parent_update_moving_tim from Parent-Child Cpl Export State"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeGet(state=EXP_STATE_CPL_NEST               &  !<-- The Parent-Child Coupler export state
                                ,name ='parent_update_moving_tim'       &  !<-- Name of the attribute to extract
                                ,value=td%parent_update_moving_tim      &  !<-- moving nest update time
                                ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INTEG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          MESSAGE_CHECK="Extract t0_recv_move_tim from Parent-Child Cpl Export State"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeGet(state=EXP_STATE_CPL_NEST               &  !<-- The Parent-Child Coupler export state
                                ,name ='t0_recv_move_tim'               &  !<-- Name of the attribute to extract
                                ,value=td%t0_recv_move_tim              &  !<-- task 0 time to process receive of move flag
                                ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INTEG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          MESSAGE_CHECK="Extract read_moving_child_topo_tim from Parent-Child Cpl Export State"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeGet(state=EXP_STATE_CPL_NEST               &  !<-- The Parent-Child Coupler export state
                                ,name ='read_moving_child_topo_tim'     &  !<-- Name of the attribute to extract
                                ,value=td%read_moving_child_topo_tim    &  !<-- task 0 time to process receive of move flag
                                ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INTEG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          MESSAGE_CHECK="Extract barrier_move_tim from Parent-Child Cpl Export State"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeGet(state=EXP_STATE_CPL_NEST               &  !<-- The Parent-Child Coupler export state
                                ,name ='barrier_move_tim'     &  !<-- Name of the attribute to extract
                                ,value=td%barrier_move_tim    &  !<-- task 0 time to process receive of move flag
                                ,rc   =RC)

          CALL ESMF_AttributeGet(state=EXP_STATE_CPL_NEST               &  !<-- The Parent-Child Coupler export state
                                ,name ='pscd_tim'     &  !<-- Name of the attribute to extract
                                ,value=td%pscd_tim    &  !<-- task 0 time to process receive of move flag
                                ,rc   =RC)
          CALL ESMF_AttributeGet(state=EXP_STATE_CPL_NEST               &  !<-- The Parent-Child Coupler export state
                                ,name ='pscd1_tim'     &  !<-- Name of the attribute to extract
                                ,value=td%pscd1_tim    &  !<-- task 0 time to process receive of move flag
                                ,rc   =RC)
          CALL ESMF_AttributeGet(state=EXP_STATE_CPL_NEST               &  !<-- The Parent-Child Coupler export state
                                ,name ='pscd2_tim'     &  !<-- Name of the attribute to extract
                                ,value=td%pscd2_tim    &  !<-- task 0 time to process receive of move flag
                                ,rc   =RC)
          CALL ESMF_AttributeGet(state=EXP_STATE_CPL_NEST               &  !<-- The Parent-Child Coupler export state
                                ,name ='pscd3_tim'     &  !<-- Name of the attribute to extract
                                ,value=td%pscd3_tim    &  !<-- task 0 time to process receive of move flag
                                ,rc   =RC)
          CALL ESMF_AttributeGet(state=EXP_STATE_CPL_NEST               &  !<-- The Parent-Child Coupler export state
                                ,name ='pscd4_tim'     &  !<-- Name of the attribute to extract
                                ,value=td%pscd4_tim    &  !<-- task 0 time to process receive of move flag
                                ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INTEG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!

        ENDIF
!
!-----------------------------------------------------------------------
!
        IF(.NOT. NESTING) THEN                                             !<-- Parent only run
          IF (td%domain_run_1 < 1.0) THEN                                  !<-- An I/O task
            IF(PRINT_TIMING) &
            WRITE(0,899)td%domain_run_3
          ELSE
            IF (td%domain_run_2 > 1.0) THEN                                !<-- Digital filter
              IF(PRINT_TIMING) &
              WRITE(0,900)td%domain_run_1,td%domain_run_2,td%domain_run_3
            ELSE
              IF(PRINT_TIMING) &
              WRITE(0,901)td%domain_run_1,td%domain_run_3                  !<-- Primary compute task
            ENDIF
          ENDIF
!
        ELSE
          IF(I_AM_A_FCST_TASK)THEN                                          !<-- Nested run and a forecast task
           IF (td%domain_run_2 > 1.0) THEN                                  !<-- Digital filter
             WRITE(0,800)my_domain_id,td%pc_cpl_run_cpl1,td%pc_cpl_run_cpl2    &
                        ,td%pc_cpl_run_cpl3,td%domain_run_1,td%pc_cpl_run_cpl4 &
                        ,td%pc_cpl_run_cpl5,td%domain_run_2,td%domain_run_3
           ELSE
            IF(td%domain_run_3 > 1.0) then
             WRITE(0,801)my_domain_id,td%pc_cpl_run_cpl1,td%pc_cpl_run_cpl2    &
                        ,td%pc_cpl_run_cpl3,td%domain_run_1,td%pc_cpl_run_cpl4 &
                        ,td%pc_cpl_run_cpl5,td%domain_run_3
             WRITE(0,803)my_domain_id,td%parent_bookkeep_moving_tim            &
                        ,td%cpl2_comp_tim,td%cpl2_send_tim,td%cpl2_wait_tim    &
                        ,td%barrier_move_tim
!jaa             WRITE(0,803)my_domain_id,td%parent_bookkeep_moving_tim            &
!jaa                        ,td%parent_update_moving_tim,td%t0_recv_move_tim       &
!jaa                        ,td%read_moving_child_topo_tim,td%barrier_move_tim
            WRITE(0,804)my_domain_id,td%pscd_tim,td%pscd1_tim,td%pscd2_tim    &
                       ,td%pscd3_tim,td%pscd4_tim
            ELSE
             WRITE(0,802)my_domain_id,td%pc_cpl_run_cpl1,td%pc_cpl_run_cpl2    &
                        ,td%pc_cpl_run_cpl3,td%domain_run_1,td%pc_cpl_run_cpl4 &
                        ,td%pc_cpl_run_cpl5
             WRITE(0,803)my_domain_id,td%parent_bookkeep_moving_tim            &
                        ,td%cpl2_comp_tim,td%cpl2_send_tim,td%cpl2_wait_tim    &
                        ,td%barrier_move_tim
!jaa             WRITE(0,803)my_domain_id,td%parent_bookkeep_moving_tim            &
!jaa                        ,td%parent_update_moving_tim,td%t0_recv_move_tim       &
!jaa                        ,td%read_moving_child_topo_tim,td%barrier_move_tim
             WRITE(0,804)my_domain_id,td%pscd_tim,td%pscd1_tim,td%pscd2_tim    &
                        ,td%pscd3_tim,td%pscd4_tim

            ENDIF
           ENDIF
          ELSE                                                               !<-- I/O tasks
           WRITE(0,899)td%domain_run_3
          ENDIF
        ENDIF
!
 897    FORMAT(' The timers that may be printed include: '/          &
               ' Integrate - timers around gridcomp_run phase=1 for parent tasks in NMM_Integrate ',/&
               ' Filter - timers around gridcomp_run phase=2 in NMM_Integrate ',/&
               ' Phase 3(I/O) - timers around gridcomp_run phase=3 in NMM_Integrate ',/&
               ' cpl compute - timers around gridcomp_run phase=1 for nest tasks in NMM_Integrate or ',/&
               ' cpl compute - time for parent tasks to compute boundary updates ',/&
               ' cpl recv - time nest tasks spend waiting to receive boundary updates from parents ',/&
               ' update interior nest - time nest task spends updating from other nest tasks in a moving nest ',/&
               ' update interior parent - time nest task spends updating from parent tasks in a moving nest ',/&
               ' update parent move- time parent task spends updating nest internal points in a moving nest ',/&
               ' cpl wait - time parent task spends waiting for nest tasks to receive boundary data ',/&
               ' cpl 2-way send  - time child task computing/sending exchange data to parent ')

 800    FORMAT(' For domain ',i2,' c1 ',g10.3,' c2 ',g10.3,' c3 ',g10.3  &
               ,' run ',g10.3,' c4 ',g10.3,' c5 ',g10.3,' df ',g10.3     &
               ,' i/o ',g10.3)
 801    FORMAT(' For domain ',i2,' c1 ',g10.3,' c2 ',g10.3,' c3 ',g10.3  &
               ,' run ',g10.3,' c4 ',g10.3,' c5 ',g10.3,' i/o ',g10.3)
 802    FORMAT(' For domain ',i2,' c1 ',f10.5,' c2 ',f10.3,' c3 ',f10.5  &
               ,' run ',f10.5,' c4 ',f10.5,' c5 ',f10.5)
!jaa 803    FORMAT(' For domain ',i2,' pbm ',f10.5,' pumt ',f10.5,' rmt ',f10.5  &
!jaa               ,' rmctt ',f10.5,' bmt ',f10.5)
 803    FORMAT(' For domain ',i2,' pbm ',f10.5,' comp ',f10.5,' send ',f10.5  &
               ,' wait ',f10.5,' bmt ',f10.5)
 804    FORMAT(' For domain ',i2,' pscd ',f10.5,' pscd1 ',f10.5,' pscd2 ',f10.5  &
               ,' pscd3 ',f10.5,' pscd4 ',f10.5)
 899    FORMAT(' I/O task Phase 3= ',g10.3)
 900    FORMAT(' Integrate = ',g10.3,' Filter = ',g10.3,                 &
               ' Phase 3 = ',g10.3,' update parent move = ',g10.3)
 901    FORMAT(' Integrate = ',g10.3,' Phase 3 = ',g10.3,                &
               ' update parent move = ',g10.3)

!
!-----------------------------------------------------------------------
!
      ELSE  driver_run_end                                                 !<-- Filtering is in effect
!
!-----------------------------------------------------------------------
!***  If we are completing the execution of digital filtering then reset
!***  the Clock and times.
!-----------------------------------------------------------------------
!
        IF(CLOCK_DIRECTION=='Bckward')THEN
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="NMM_INTEGRATE: Get CurrTime and Timestep for Bckward"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_ClockGet(clock   =CLOCK_INTEGRATE                   &
                            ,currtime=CURRTIME                          &
                            ,timestep=TIMESTEP_FILTER                   &
                            ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INTEG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          TIMESTEP_FILTER=-TIMESTEP_FILTER                                 !<-- We must set the timestep back to positive.
!
!-----------------------------------------------------------------------
          filter_method_block : IF(FILTER_METHOD==3)THEN
!-----------------------------------------------------------------------
!
            ndfiloop: DO I=1,NDFISTEP
!
!-----------------------------------------------------------------------
!***  Now set the timestep of the children at which they receive
!***  data from their parent.  This must be known by the parents
!***  since it will provide the proper tag to the MPI data sent.
!-----------------------------------------------------------------------
!

              parents_only: IF(NUM_CHILDREN>0                           &
                                   .AND.                                &
                               I_AM_A_FCST_TASK)THEN
!
!-----------------------------------------------------------------------
!
                IF(.NOT.ALLOCATED(LOC_PAR_CHILD_TIME_RATIO)) THEN
                  ALLOCATE(LOC_PAR_CHILD_TIME_RATIO(1:NUM_CHILDREN))
                ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="NMM_INTEGRATE: Parent/child DT Ratio for TDFI"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
                CALL ESMF_AttributeGet(state    =IMP_STATE_CPL_NEST        &  !<-- The parent-child coupler import state
                                      ,name     ='Parent-Child Time Ratio' &  !<-- Name of the attribute to extract
                                      ,itemCount=NUM_CHILDREN              &  !<-- # of items in the Attribute
                                      ,valueList=LOC_PAR_CHILD_TIME_RATIO  &  !<-- Ratio of parent to child DTs
                                      ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INTEG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
                DO N=1,NUM_CHILDREN
                  NSTEP_CHILD_RECV(N)=NSTEP_CHILD_RECV(N) +  LOC_PAR_CHILD_TIME_RATIO(N)
                ENDDO
!
!-----------------------------------------------------------------------
!
              ENDIF parents_only
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
              MESSAGE_CHECK="NMM_INTEGRATE: Advance Clock for TDFI"
!             CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
              CALL ESMF_ClockAdvance(clock   =CLOCK_INTEGRATE           &
                                    ,timestep=TIMESTEP_FILTER           &  !<-- Advance the clock to the forward starttime
                                    ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
              CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INTEG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
            ENDDO ndfiloop
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="NMM_INTEGRATE: Get Current Clock Time for TDFI"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_ClockGet(clock   =CLOCK_INTEGRATE                 &
                              ,currtime=CURRTIME                        &
                              ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INTEG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
          ENDIF  filter_method_block
!
!-----------------------------------------------------------------------
!
        ELSEIF(CLOCK_DIRECTION=='Forward ')THEN
!
!-----------------------------------------------------------------------
!
          IF(FILTER_METHOD==1)THEN
!
            domain_int_state%FIRST_PASS=.TRUE.
!
            CURRTIME=HALFDFITIME-TIMESTEP
            NTIMESTEP=NTIMESTEP-(HALFDFIINTVAL/TIMESTEP)-1
            NTIMESTEP_ESMF=NTIMESTEP
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="NMM_INTEGRATE: Set Time to Half Filter Interval"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_ClockSet(clock       =CLOCK_INTEGRATE             &  !<-- Reset current time and timestep to the
                              ,currtime    =CURRTIME                    &  !    halfway point of the filter interval.
                              ,advanceCount=NTIMESTEP_ESMF              &
                              ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INTEG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          ENDIF
!
!-----------------------------------------------------------------------
!
        ENDIF
!
!-----------------------------------------------------------------------
!
      ENDIF  driver_run_end
!
!-----------------------------------------------------------------------
!***  The final error signal information.
!-----------------------------------------------------------------------
!
      IF(RC_INTEG==ESMF_SUCCESS)THEN
!       WRITE(0,*)'NMM RUN step succeeded'
      ELSE
        WRITE(0,*)'NMM RUN step failed RC_INTEG=',RC_INTEG
      ENDIF
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
      SUBROUTINE RESET_ALARMS
!
!-----------------------------------------------------------------------
!***  For normal forecast integration set the Alarm ring times
!***  while accounting for restarts and digital filtering.
!-----------------------------------------------------------------------
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: iyear_fcst &
                           ,imonth_fcst &
                           ,iday_fcst &
                           ,ihour_fcst &
                           ,iminute_fcst &
                           ,isecond_fcst &
                           ,isecond_num &
                           ,isecond_den
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      IF(CURRTIME==STARTTIME)THEN
        ALARM_HISTORY_RING  =CURRTIME
        ALARM_RESTART_RING  =CURRTIME
        ALARM_CLOCKTIME_RING=CURRTIME
!     write(0,36361)
36361 format(' RESET_ALARMS 1 set ALARM_HISTORY_RING to CURRTIME')
      ELSE
        IF(RESTARTED_RUN)THEN
          ALARM_HISTORY_RING  =CURRTIME+INTERVAL_HISTORY
          ALARM_RESTART_RING  =CURRTIME+INTERVAL_RESTART
          ALARM_CLOCKTIME_RING=CURRTIME+INTERVAL_CLOCKTIME
        ELSE
          ALARM_HISTORY_RING  =STARTTIME+INTERVAL_HISTORY
          ALARM_RESTART_RING  =STARTTIME+INTERVAL_RESTART
          ALARM_CLOCKTIME_RING=STARTTIME+INTERVAL_CLOCKTIME
!     write(0,36362)
36362 format(' RESET_ALARMS 2 set ALARM_HISTORY_RING to STARTTIME+INTERVAL_HISTORY')
        ENDIF
      ENDIF
!
#if 0
!-------------------------------------------------
!***  Adjust time of History Alarm if necessary.
!-------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Get time from ALARM_HISTORY_RING."
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_TimeGet(time=ALARM_HISTORY_RING                         &  !<-- Extract the time from this variable
                       ,yy  =YY                                         &  !<-- Year
                       ,mm  =MM                                         &  !<-- Month
                       ,dd  =DD                                         &  !<-- Day
                       ,h   =H                                          &  !<-- Hour
                       ,m   =M                                          &  !<-- Minute
                       ,s   =S                                          &  !<-- Second
                       ,rc  =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INTEG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF(M/=0)THEN
        H=H+1
        M=0
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Reset time in ALARM_HISTORY_RING."
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_TimeSet(time=ALARM_HISTORY_RING                       &  !<-- Reset the time for initial history output
                         ,yy  =YY                                       &  !<-- Year
                         ,mm  =MM                                       &  !<-- Month
                         ,dd  =DD                                       &  !<-- Day
                         ,h   =H                                        &  !<-- Hour
                         ,m   =M                                        &  !<-- Minute
                         ,s   =S                                        &  !<-- Second
                         ,rc  =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INTEG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF
!
!-------------------------------------------------
!***  Adjust time of Restart Alarm if necessary.
!-------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Get time from ALARM_RESTART_RING."
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_TimeGet(time=ALARM_RESTART_RING                         &  !<-- Extract the time from this variable
                       ,yy  =YY                                         &  !<-- Year
                       ,mm  =MM                                         &  !<-- Month
                       ,dd  =DD                                         &  !<-- Day
                       ,h   =H                                          &  !<-- Hour
                       ,m   =M                                          &  !<-- Minute
                       ,s   =S                                          &  !<-- Second
                       ,rc  =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INTEG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF(M/=0)THEN
        H=H+1
        M=0
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Reset time in ALARM_RESTART_RING."
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_TimeSet(time=ALARM_RESTART_RING                       &  !<-- Reset the time for initial restart output
                         ,yy  =YY                                       &  !<-- Year
                         ,mm  =MM                                       &  !<-- Month
                         ,dd  =DD                                       &  !<-- Day
                         ,h   =H                                        &  !<-- Hour
                         ,m   =M                                        &  !<-- Minute
                         ,s   =S                                        &  !<-- Second
                         ,rc  =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INTEG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF
!
!-------------------------------------------------------------
!***  Adjust time of Alarm for clocktime writes if necessary.
!-------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Get time from ALARM_CLOCKTIME_RING."
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_TimeGet(time=ALARM_CLOCKTIME_RING                       &  !<-- Extract the time from this variable
                       ,yy  =YY                                         &  !<-- Year
                       ,mm  =MM                                         &  !<-- Month
                       ,dd  =DD                                         &  !<-- Day
                       ,h   =H                                          &  !<-- Hour
                       ,m   =M                                          &  !<-- Minute
                       ,s   =S                                          &  !<-- Second
                       ,rc  =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INTEG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF(M/=0)THEN
        H=H+1
        M=0
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Reset time in ALARM_CLOCKTIME_RING."
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_TimeSet(time=ALARM_CLOCKTIME_RING                     &  !<-- Reset the time for clocktime prints
                         ,yy  =YY                                       &  !<-- Year
                         ,mm  =MM                                       &  !<-- Month
                         ,dd  =DD                                       &  !<-- Day
                         ,h   =H                                        &  !<-- Hour
                         ,m   =M                                        &  !<-- Minute
                         ,s   =S                                        &  !<-- Second
                         ,rc  =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INTEG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF
#endif
!
!-----------------------------------------------------------------------
!***  Now create the three Alarms using the final ringtimes.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="RESET_ALARMS: Create ALARM_HISTORY"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      domain_int_state%ALARM_HISTORY=                                      &
                    ESMF_AlarmCreate(name             ='ALARM_HISTORY'     &
                                    ,clock            =CLOCK_INTEGRATE     &  !<-- DOMAIN Clock
                                    ,ringTime         =ALARM_HISTORY_RING  &  !<-- Forecast/Restart start time (ESMF)
                                    ,ringInterval     =INTERVAL_HISTORY    &  !<-- Time interval between history output
                                    ,ringTimeStepCount=1                   &  !<-- The Alarm rings for this many timesteps
                                    ,sticky           =.false.             &  !<-- Alarm does not ring until turned off
                                    ,rc               =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INTEG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="RESET_ALARMS: Create ALARM_RESTART"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      domain_int_state%ALARM_RESTART=                                      &
                    ESMF_AlarmCreate(name             ='ALARM_RESTART'     &
                                    ,clock            =CLOCK_INTEGRATE     &  !<-- DOMAIN Clock
                                    ,ringTime         =ALARM_RESTART_RING  &  !<-- Forecast/Restart start time (ESMF)
                                    ,ringInterval     =INTERVAL_RESTART    &  !<-- Time interval between  restart output (ESMF)
                                    ,ringTimeStepCount=1                   &  !<-- The Alarm rings for this many timesteps
                                    ,sticky           =.false.             &  !<-- Alarm does not ring until turned off
                                    ,rc               =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INTEG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="RESET_ALARMS: Create ALARM_CLOCKTIME"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      domain_int_state%ALARM_CLOCKTIME=                                        &
                      ESMF_AlarmCreate(name             ='ALARM_CLOCKTIME'     &
                                      ,clock            =CLOCK_INTEGRATE       &  !<-- DOMAIN Clock
                                      ,ringTime         =ALARM_CLOCKTIME_RING  &  !<-- Forecast start time (ESMF)
                                      ,ringInterval     =INTERVAL_CLOCKTIME    &  !<-- Time interval between clocktime prints (ESMF)
                                      ,ringTimeStepCount=1                     &  !<-- The Alarm rings for this many timesteps
                                      ,sticky           =.false.               &  !<-- Alarm does not ring until turned off
                                      ,rc               =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INTEG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE RESET_ALARMS
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE NMM_INTEGRATE
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      END MODULE MODULE_NMM_INTEGRATE
!
!-----------------------------------------------------------------------
