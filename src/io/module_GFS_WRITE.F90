!#include "../../ESMFVersionDefine.h"

!-----------------------------------------------------------------------
!
      MODULE MODULE_GFS_WRITE
!
!-----------------------------------------------------------------------
!***  THIS MODULE CONTAINS ROUTINES NEEDED BY THE RUN STEP OF THE
!***  WRITE GRIDDED COMPONENT IN WHICH HISTORY OUTPUT DATA FROM
!***  THE FORECAST TASKS ARE ASSEMBLED AND WRITTEN TO HISTORY FILES
!***  BY THE WRITE TASKS.
!-----------------------------------------------------------------------
!***
!***  HISTORY   
!***
!       07 May 2009:  J. Wang - adopt write _routine from NMMB_io
!                                modified for GFS
!       03 Sep 2009:  W. Yang - Ensemble GEFS.
!       12 May 2011:  Theurich & Yang - Modified for using the ESMF 5.2.0r_beta_snapshot_07.
!       03 Sep 2011:  W. Yang - Modified for using the ESMF 5.2.0r library.
!       08 Oct 2012:  J. Wang - move module_GFS_WRITE.F90 back to io, and use io related mpi info.
!
!-----------------------------------------------------------------------
!
      USE ESMF
!
      USE MODULE_IO_MPI_DEF, ONLY :   MPI_COMM_INTER_ARRAY,   &
                                      N_GROUP,                &
                                      num_pes_fcst,           &
                                      last_fcst_pe,           &
                                      write_tasks_per_group,  &
                                      write_groups,           &
                                      petlist_write,          &
                                      NUM_PES_WRT
!
      USE MODULE_INCLUDE_IO
!
      USE MODULE_ERR_MSG,ONLY: ERR_MSG,MESSAGE_CHECK

      USE module_WRITE_GRID_COMP_GFS, ONLY: WRITE_REGISTER_GFS
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------

      PRIVATE
!
      PUBLIC :: WRITE_ASYNC_GFS, WRITE_INIT_GFS,              &
                WRITE_SETUP_GFS, WRITE_DESTROY_GFS 
!
!-----------------------------------------------------------------------
!
      LOGICAL,SAVE :: QUILTING, ADIABATIC    
      INTEGER      :: MYPE                                           !<-- My MPI task ID
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE WRITE_INIT_GFS(ATM_GRID_COMP,wrt_comps,imp_state_write, &
        exp_state_write,CLOCK_GFS,WRITE_GROUP_READY_TO_GO)

! 
!-----------------------------------------------------------------------
!***  EXECUTE THE INITIALIZE STEP OF THE WRITE COMPONENTS.
!-----------------------------------------------------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT)      :: ATM_GRID_COMP             !<-- The ATM gridded component
      TYPE(ESMF_GRIDComp),INTENT(INOUT)      :: wrt_comps(:)             !<-- The ATM Internal State
      TYPE(ESMF_STATE),INTENT(INOUT)         :: imp_state_write             !<-- The ATM Internal State
      TYPE(ESMF_STATE),INTENT(INOUT)         :: exp_state_write             !<-- The ATM Internal State
      TYPE(ESMF_Clock),INTENT(INOUT)         :: CLOCK_GFS                 !<-- The ATM Component's ESMF Clock
      INTEGER, INTENT(INOUT)                 :: WRITE_GROUP_READY_TO_GO
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      TYPE(ESMF_Config)      :: CF                                        !<-- The config object
      TYPE(ESMF_Time)        :: CURRTIME
      TYPE(ESMF_TimeInterval):: TIMEINTERVAL_OUTPUT
!
      INTEGER :: I,J,RC,RC_INIT,NFHOUT
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  RETRIEVE THE CONFIG OBJECT CF FROM THE ATM GRIDDED COMPONENT.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Get Config Object from ATM Component in Write Init"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompGet(gridcomp=ATM_GRID_COMP                      &  !<-- The ATM gridded component
                           ,config  =CF                                 &  !<-- The config object (~namelist)
                           ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  EXECUTE THE INITIALIZE STEP FOR THE WRITE COMPONENTS.
!***  THESE ARE THE INITIALIZE SUBROUTINES SPECIFIED IN THE
!***  REGISTER ROUTINES CALLED IN ESMF_GridCompSetServices.
!-----------------------------------------------------------------------
!
      DO J=1,WRITE_GROUPS
!
        DO I=1,NUM_PES_WRT
          IF(MYPE==PETLIST_WRITE(I,J))THEN                   !<--  Forecast tasks plus the Write tasks in each write group
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Execute Initialize Step of Write Component"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            N_GROUP=J
            CALL ESMF_GridCompInitialize(WRT_COMPS(J)                 &  !<-- The Write gridded components
                                        ,importstate=IMP_STATE_WRITE  &  !<-- The Write import state
                                        ,exportstate=EXP_STATE_WRITE  &  !<-- The Write export state
                                        ,clock      =CLOCK_GFS                      &  !<-- The ESMF clock of the ATM component
                                        ,rc         =RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          ENDIF
        ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
!***  SET THE FIRST WRITE GROUP AS THE FIRST ONE TO ACT.
!-----------------------------------------------------------------------
!
      WRITE_GROUP_READY_TO_GO=1
      N_GROUP=WRITE_GROUP_READY_TO_GO
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE WRITE_INIT_GFS
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE WRITE_ASYNC_GFS(WRT_COMPS,EXP_STATE,IMP_STATE_WRITE     &
                            ,EXP_STATE_WRITE,CLOCK_GFS,MYPE              &
                            ,WRITE_GROUP_READY_TO_GO)
!
!-----------------------------------------------------------------------
!***  WRITE OUT A HISTORY FILE USING THE ASYNCHRONOUS QUILTING.
!-----------------------------------------------------------------------
!
      TYPE(ESMF_GridComp)     ,INTENT(INOUT) :: WRT_COMPS(:)               !<-- The write_comp
      TYPE(ESMF_STATE)        ,INTENT(INOUT) :: EXP_STATE                  !<-- The export state dyn
      TYPE(ESMF_STATE)        ,INTENT(INOUT) :: IMP_STATE_WRITE            !<-- The import state of write
      TYPE(ESMF_STATE)        ,INTENT(INOUT) :: EXP_STATE_WRITE            !<-- The export state of write
      TYPE(ESMF_Clock)        ,INTENT(INOUT) :: CLOCK_GFS                 !<-- The ATM Component's ESMF Clock
      INTEGER,INTENT(IN) :: MYPE
      INTEGER,INTENT(INOUT) :: WRITE_GROUP_READY_TO_GO
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      TYPE(ESMF_Config) :: CF                                             !<-- The configure object (~namelist)
      TYPE(ESMF_Time)   :: CURRTIME                                       !<-- The current forecast time (ESMF)
!
      INTEGER :: YY,MM,DD,H,M,S                                           !<-- Year, Month, Day, Hour, Minute, Second (integer)
!
      INTEGER :: I,RC,RC_ASYNC
!
      CHARACTER(ESMF_MAXSTR) :: filename                                  !<-- Restart/History label
!jwtest
      TYPE(ESMF_VM) :: myVM
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  WHAT IS THE CURRENT FORECAST TIME?
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="WRITE_ASYNC: Get Current Time from ATM Clock"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ClockGet(clock   =CLOCK_GFS                             &  !<-- The ATM component's ESMF Clock
                        ,currTime=CURRTIME                              &  !<-- The current forecast time (ESMF)
                        ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_ASYNC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="WRITE_ASYNC: Convert ESMF Time to Real Time"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_TimeGet (time = CURRTIME  & !<-- current forecast time (ESMF)
                        ,yy   = YY        & !<-- current year (integer)
                        ,mm   = MM        & !<-- current month (integer)
                        ,dd   = DD        & !<-- current day (integer)
                        ,h    = H         & !<-- current hour (integer)
                        ,m    = M         & !<-- current minute (integer)
                        ,s    = S         & !<-- current second (integer)
                        ,rc   = RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_ASYNC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  THE EXPORT STATE OF THE DYNAMICS COMPONENT LIES WITHIN THE
!***  INTERNAL STATE OF THE ATM GRIDDED COMPONENT AND HOLDS THE
!***  IMPORT STATE OF THE WRITE COMPONENT.
!***  EXTRACT THAT WRITE COMPONENT'S IMPORT STATE SINCE WE ARE 
!***  ABOUT TO EXECUTE THE RUN STEP OF THE WRITE COMPONENT.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="WRITE_ASYNC: Extract Write Import State from Dyn Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateGet(state       = EXP_STATE             &  !<-- The Dyn component's export state
                        ,itemName    = "Write Import State"  &  !<-- Name of state to be extracted
                        ,nestedState = IMP_STATE_WRITE       &  !<-- The extracted state
                        ,rc          = RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_ASYNC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  ALL FORECAST TASKS PLUS THOSE WRITE TASKS IN THE APPROPRIATE
!***  WRITE GROUP EXECUTE THE RUN STEP OF A WRITE COMPONENT.
!-----------------------------------------------------------------------
!
      N_GROUP = WRITE_GROUP_READY_TO_GO       !<-- The active write group
!
!jw
      call esmf_GridCompGet(gridcomp = WRT_COMPS(N_GROUP)          &
                           ,vm       = myVM                        &
                           ,rc       = rc)

    CALL ESMF_VMBarrier(myVM,rc=RC)    ! Insert barrier since fcst tasks are involved in each iteration of write groups
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="WRITE_ASYNC: Execute Run Step of Write Components" 
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!

      DO I=1, NUM_PES_WRT
        IF(MYPE == PETLIST_WRITE(I,N_GROUP)) THEN
          CALL ESMF_GridCompRun(WRT_COMPS(N_GROUP)          &  !<-- The write gridded component
                               ,importState=IMP_STATE_WRITE &  !<-- Its import state
                               ,exportState=EXP_STATE_WRITE &  !<-- Its export state
                               ,clock      =CLOCK_GFS                     &  !<-- The ATM Clock
                               ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_ASYNC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          IF(I == NUM_PES_FCST+1) THEN    !<-- First write task tells history output time
            WRITE(0,101)YY,MM,DD,H,M,S
  101       FORMAT(' Wrote File at ',I4.4,'_',I2.2,'_',I2.2,'_',I2.2,':',I2.2,':',I2.2)
          ENDIF
!
        ENDIF
!
      ENDDO
!
!-----------------------------------------------------------------------
!***  PREPARE TO USE THE NEXT WRITE GROUP AT THE NEXT OUTPUT TIME.
!***  RETURN TO THE 1ST GROUP IF WE HAVE CYCLED THROUGH ALL OF THEM.
!-----------------------------------------------------------------------
!
      IF(WRITE_GROUP_READY_TO_GO == WRITE_GROUPS) THEN
        WRITE_GROUP_READY_TO_GO = 1
      ELSE
        WRITE_GROUP_READY_TO_GO = WRITE_GROUP_READY_TO_GO+1
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE WRITE_ASYNC_GFS
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE WRITE_SETUP_GFS(ATM_GRID_COMP,WRT_COMPS           &
                                ,exp_state_dyn,exp_state_phy       &
                                ,imp_state_write,exp_state_write)
! 
!-----------------------------------------------------------------------
!***  SET UP THE WRITE COMPONENTS WITH THE FORECAST TASKS AND
!***  THE GROUPS OF WRITE TASKS NEEDED FOR QUILTING THE OUTPUT
!***  AND WRITING IT TO HISTORY FILES.
!-----------------------------------------------------------------------
!
!
      TYPE(ESMF_GridComp),INTENT(INOUT)      :: ATM_GRID_COMP             !<-- The ATM gridded component
      TYPE(ESMF_GridComp),INTENT(inOUT)      :: WRT_COMPS(:)              !<-- The ATM gridded component
      TYPE(ESMF_State),INTENT(INOUT)         :: EXP_state_DYN
      TYPE(ESMF_State),INTENT(INOUT)         :: EXP_state_PHY
      TYPE(ESMF_State),INTENT(INOUT)         :: IMP_state_WRITE
      TYPE(ESMF_State),INTENT(INOUT)         :: EXP_state_WRITE
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      TYPE(ESMF_Config)      :: CF                                        !<-- The config object
!
      INTEGER                :: WRITE_GROUPS                            & !<-- Number of groups of write tasks
                               ,WRITE_TASKS_PER_GROUP                     !<-- #of tasks in each write group
      INTEGER,ALLOCATABLE    :: ITMP(:)
      LOGICAL                :: STANDALONE_POST
!
      CHARACTER( 2)          :: MY_WRITE_GROUP
      CHARACTER(6)           :: FMT='(I2.2)'
      CHARACTER(ESMF_MAXSTR) :: WRITE_NAME
      CHARACTER(50)          :: MODE
!
      INTEGER :: I,J,K,RC,RC_SETUP
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  RETRIEVE THE CONFIG OBJECT CF FROM THE ATM GRIDDED COMPONENT.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Get Config Object for Write Setup"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompGet(gridcomp = ATM_GRID_COMP   & !<-- The ATM gridded component
                           ,config   = CF              & !<-- The config object (~namelist)
                           ,rc       = RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_SETUP)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  RETRIEVE TASK AND GROUP COUNTS FROM THE CONFIG FILE.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Write Tasks and Groups from Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(CF                                   &  !<-- The configure file
                                  ,QUILTING                             &  !<-- Number of write groups from config file
                                  ,label ='quilting:'                   &
                                  ,rc    =RC)

      CALL ESMF_ConfigGetAttribute(CF                                   &  !<-- The configure file
                                  ,WRITE_GROUPS                         &  !<-- Number of write groups from config file
                                  ,label ='write_groups:'               &
                                  ,rc    =RC)
!
      CALL ESMF_ConfigGetAttribute(CF                                   &  !<-- The configure file
                                  ,WRITE_TASKS_PER_GROUP                &  !<-- Number of write tasks per group from config file
                                  ,label ='write_tasks_per_group:'      &
                                  ,rc    =RC)
!

      CALL ESMF_ConfigGetAttribute(CF, ADIABATIC, label ='adiabatic:', rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_SETUP)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  ASSOCIATE ALL OF THE FORECAST TASKS WITH THE WRITE TASKS
!***  IN EACH WRITE GROUP.
!-----------------------------------------------------------------------
!
      IF(QUILTING) THEN
        NUM_PES_WRT=NUM_PES_FCST+WRITE_TASKS_PER_GROUP
      ELSE
        NUM_PES_WRT=NUM_PES_FCST
      ENDIF
      ALLOCATE(PETLIST_WRITE(NUM_PES_WRT,WRITE_GROUPS))                        !<-- Task IDs of all wrt tasks
!
!-----------------------------------------------------------------------
!***  COLLECT THE TASK IDs FOR THE WRITE TASKS AND THE ASSOCIATED
!***  FORECAST-WRITE TASKS.
!-----------------------------------------------------------------------
!
      DO I=0,NUM_PES_FCST-1
        DO J=1,WRITE_GROUPS
          PETLIST_WRITE(I+1,J)=I+last_fcst_pe+1-num_pes_fcst            !<-- Collect forecast task IDs to be associated with
                                                                        !    write tasks by write group
        ENDDO
!
      ENDDO
!
      IF(NUM_PES_WRT>NUM_PES_FCST) THEN
          K=NUM_PES_FCST
!
          DO J=1,WRITE_GROUPS
              DO I=1,WRITE_TASKS_PER_GROUP
                  PETLIST_WRITE(NUM_PES_FCST+I,J)=K               !<-- Append write task IDs to associated forecast task IDs by group
                  K=K+1
              END DO
          END DO
      END IF
!      write(0,*)'in write setup, last_fcst_pe=',last_fcst_pe,'num_pes_fcst=', &
!        num_pes_fcst,'WRITE_GROUPS=',WRITE_GROUPS,'petlist_write=',  &
!         PETLIST_WRITE(:,1)
!
!-----------------------------------------------------------------------
!***  CREATE THE WRITE GRIDDED COMPONENT(S).
!***  THERE ARE AS MANY WRITE COMPONENTS AS THERE ARE GROUPS OF
!***  WRITE TASKS SPECIFIED IN THE CONFIGURE FILE.
!***  REGISTER THEIR INIT, RUN, AND FINALIZE STEPS.
!-----------------------------------------------------------------------
!
!---------------------------------
!***  Create the Write components
!---------------------------------
!
      allocate(ITMP(NUM_PES_WRT))
      DO I=1,WRITE_GROUPS
        ITMP(1:NUM_PES_WRT)=PETLIST_WRITE(1:NUM_PES_WRT,I)
        WRITE(MY_WRITE_GROUP,FMT)I
        WRITE_NAME='write_GridComp_'//MY_WRITE_GROUP
!
        WRT_COMPS(I)=ESMF_GridCompCreate(                         &
                                name          =WRITE_NAME         &  !<-- Name of this group's Write gridded component
                               ,configFile    ='atm_namelist.rc'  &  !<-- The configure file for writes
                               ,petList       =ITMP               &  !<-- The task IDs of the write tasks in this group
                                                                     !    provide the local VM information per component.
                               ,rc            =RC)
!
      ENDDO
      deallocate(itmp)
!
!-----------------------------------
!***  Register the Write components
!-----------------------------------
!
      DO I=1,WRITE_GROUPS
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Register Write Components"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_GridCompSetServices(WRT_COMPS(I)         &  !<-- The Write gridded components
                                     ,WRITE_REGISTER_GFS   &  !<-- The user's subroutine name
                                     ,rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_SETUP)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDDO
!
!------------------------------------------------------------------------
!***  Create empty Import and Export states for the Write subcomponent(s)
!------------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create Empty Import/Export States for Write Components"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IMP_STATE_WRITE=ESMF_StateCreate(name   = 'Write Import State'    &  !<-- Import state name for writes
                                      ,stateintent = ESMF_STATEINTENT_IMPORT &
                                      ,rc          = RC)
!
      EXP_STATE_WRITE=ESMF_StateCreate(name   = 'Write Export State'    &  !<-- Export state names for writes
                                      ,stateintent = ESMF_STATEINTENT_EXPORT &
                                      ,rc          = RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_SETUP)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  INSERT THE WRITE COMPONENTS' IMPORT STATE INTO THE
!***  DYNAMICS' AND PHYSICS' EXPORT STATES SINCE HISTORY
!***  DATA ITSELF MUST COME FROM THE DYNAMICS AND PHYSICS.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Insert Write Import State into Dynamics Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateAddReplace(EXP_STATE_DYN        & !<-- Dynamics export state receives a state
                               ,(/IMP_STATE_WRITE/)  & !<-- Add the write components' import state
                               ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_SETUP)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF(.not. ADIABATIC) THEN
        MESSAGE_CHECK="Insert Write Import State into Physics Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_StateAddReplace(EXP_STATE_PHY        & !<-- Physics export state receives a state
                                 ,(/IMP_STATE_WRITE/)  & !<-- Add the write components' import state
                                 ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_SETUP)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE WRITE_SETUP_GFS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE WRITE_DESTROY_GFS(ATM_GRID_COMP,WRT_COMPS,             &
        IMP_STATE_WRITE,EXP_STATE_WRITE,CLOCK_GFS)
! 
!-----------------------------------------------------------------------
!***  DESTROY ALL OBJECTS RELATED TO THE WRITE COMPONENTS.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT)      :: ATM_GRID_COMP             !<-- The ATM gridded component
      TYPE(ESMF_GridComp),DIMENSION(:),INTENT(INOUT)      ::WRT_COMPS
      TYPE(ESMF_State),INTENT(INOUT)         :: IMP_STATE_WRITE
      TYPE(ESMF_State),INTENT(INOUT)         :: EXP_STATE_WRITE
      TYPE(ESMF_Clock),INTENT(INOUT)         :: CLOCK_GFS                 !<-- The ATM Component's ESMF Clock
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      INTEGER :: I,J,N,RC,RC_DES
!
      TYPE(ESMF_VM)                          :: VM
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC     = ESMF_SUCCESS
      RC_DES = ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  RETRIEVE THE CURRENT VM.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Retrieve the Local VM in Write Destroy"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_VMGetCurrent(vm=VM                                      &  !<-- The ESMF virtual machine
                            ,rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_DES)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  WHAT IS MY MPI TASK ID?
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Get MPI Task IDs for Write Destroy"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_VMGet(vm      =VM                                       &  !<-- The virtual machine
                     ,localpet=MYPE                                     &  !<-- Local PE rank
                     ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_DES)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  FINALIZE THE WRITE GRIDDED COMPONENTS IN EACH WRITE GROUP.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Finalize Write Components"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      DO N=1,WRITE_GROUPS              
        IF(MYPE>=PETLIST_WRITE(1,N).AND.                          &
           MYPE<=PETLIST_WRITE(NUM_PES_WRT,N))THEN
!
           CALL ESMF_GridCompFinalize(gridcomp   =WRT_COMPS(N)    &
                                     ,importstate=IMP_STATE_WRITE &
                                     ,exportstate=EXP_STATE_WRITE &
                                     ,clock      =CLOCK_GFS       &
                                     ,rc         =RC)
        ENDIF
      ENDDO
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_DES)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  DESTROY THE WRITE COMPONENTS' IMPORT/EXPORT STATES.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Destroy Write Component Import/Export States"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateDestroy(IMP_STATE_WRITE,rc=RC)
      CALL ESMF_StateDestroy(EXP_STATE_WRITE,rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_DES)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  DESTROY THE WRITE COMPONENTS.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Destroy Write Components"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      DO J=1,WRITE_GROUPS
!
        CALL ESMF_VMBarrier(vm=VM,rc=RC)
!
        DO I=1,NUM_PES_WRT
          IF(MYPE==PETLIST_WRITE(I,J))THEN
            CALL ESMF_GridCompDestroy(gridcomp=WRT_COMPS(J) &
                                     ,rc      =RC)
          ENDIF
        ENDDO
!
      ENDDO
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_DES)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  THE FINAL ERROR SIGNAL INFORMATION.
!-----------------------------------------------------------------------
!
      IF(RC_DES==ESMF_SUCCESS)THEN
        WRITE(0,*)'ATM FINALIZE STEP SUCCEEDED'
      ELSE
        WRITE(0,*)'ATM FINALIZE STEP FAILED  RC_DES=',RC_DES
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE WRITE_DESTROY_GFS
!
!-----------------------------------------------------------------------
      END MODULE MODULE_GFS_WRITE
!
!
!-----------------------------------------------------------------------
