!TBH:  Matches r14147 of https://svnemc.ncep.noaa.gov/projects/nems/trunk/.  

#include "./ESMFVersionDefine.h"

#if (ESMF_MAJOR_VERSION < 5 || ESMF_MINOR_VERSION < 2)
#undef ESMF_520rbs
#else
#define ESMF_520rbs
#endif

!-----------------------------------------------------------------------
!
      MODULE module_EARTH_GRID_COMP
!
!-----------------------------------------------------------------------
!***  This module contains codes directly related to the EARTH component.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!  2010-03-24  Black - Created Earth component module.
!  2010-04     Yang  - Added Ensemble capability.
!  2011-05-11  Theurich & Yang - Modified for using the ESMF 5.2.0r_beta_snapshot_07.
!-----------------------------------------------------------------------
!
!***  The EARTH component lies in the heirarchy seen here:
!
!          Main program
!               |
!               |
!          NEMS component
!               |     |________________________.
!               |                              |
!          EARTH component        Ensemble Coupler component
!               |
!               |
!          ATM/OCEAN/ICE components
!               |
!               |
!          CORE component (GFS, NMM, FIM, GEN, etc.)
!
!-----------------------------------------------------------------------
!
      USE ESMF_MOD
!
      USE module_EARTH_INTERNAL_STATE,ONLY: EARTH_INTERNAL_STATE        &
                                           ,WRAP_EARTH_INTERNAL_STATE
!
      USE module_ATM_GRID_COMP
!
      USE module_ERR_MSG,ONLY: ERR_MSG,MESSAGE_CHECK
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      PRIVATE
!
      PUBLIC :: EARTH_REGISTER
!
!-----------------------------------------------------------------------
!
      TYPE(EARTH_INTERNAL_STATE),POINTER,SAVE :: EARTH_INT_STATE           !<-- Internal state of the EARTH component
      TYPE(WRAP_EARTH_INTERNAL_STATE)   ,SAVE :: WRAP                      !<-- F90 pointer to the EARTH internal state
!
!-----------------------------------------------------------------------
!
      CONTAINS

!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE EARTH_REGISTER(EARTH_GRID_COMP,RC_REG)
!
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(ESMF_GridComp)               :: EARTH_GRID_COMP                 !<-- The EARTH component
      INTEGER            ,INTENT(OUT)   :: RC_REG                          !<-- Error return code
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER :: RC
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Register the EARTH Initialize, Run, and Finalize routines.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set Entry Point for EARTH Initialize"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifdef ESMF_3
      CALL ESMF_GridCompSetEntryPoint(EARTH_GRID_COMP                   &  !<-- The EARTH component
                                     ,ESMF_SETINIT                      &  !<-- Subroutine type (Initialize)
                                     ,EARTH_INITIALIZE                  &  !<-- User's subroutine name
                                     ,ESMF_SINGLEPHASE                  &
                                     ,RC)
#else
      CALL ESMF_GridCompSetEntryPoint(EARTH_GRID_COMP                   &  !<-- The EARTH component
                                     ,ESMF_SETINIT                      &  !<-- Subroutine type (Initialize)
                                     ,EARTH_INITIALIZE                  &  !<-- User's subroutine name
                                     ,phase=ESMF_SINGLEPHASE            &
                                     ,rc=RC)
#endif

!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set Entry Point for EARTH Run"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifdef ESMF_3
      CALL ESMF_GridCompSetEntryPoint(EARTH_GRID_COMP                   &  !<-- The EARTH component
                                     ,ESMF_SETRUN                       &  !<-- Subroutine type (Run)
                                     ,EARTH_RUN                         &  !<-- User's subroutine name
                                     ,ESMF_SINGLEPHASE                  &
                                     ,RC)
#else
      CALL ESMF_GridCompSetEntryPoint(EARTH_GRID_COMP                   &  !<-- The EARTH component
                                     ,ESMF_SETRUN                       &  !<-- Subroutine type (Run)
                                     ,EARTH_RUN                         &  !<-- User's subroutine name
                                     ,phase=ESMF_SINGLEPHASE            &
                                     ,rc=RC)
#endif
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set Entry Point for EARTH Finalize"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifdef ESMF_3
      CALL ESMF_GridCompSetEntryPoint(EARTH_GRID_COMP                   &  !<-- The EARTH component
                                     ,ESMF_SETFINAL                     &  !<-- Subroutine type (Finalize)
                                     ,EARTH_FINALIZE                    &  !<-- User's subroutine name
                                     ,ESMF_SINGLEPHASE                  &
                                     ,RC)
#else
      CALL ESMF_GridCompSetEntryPoint(EARTH_GRID_COMP                   &  !<-- The EARTH component
                                     ,ESMF_SETFINAL                     &  !<-- Subroutine type (Finalize)
                                     ,EARTH_FINALIZE                    &  !<-- User's subroutine name
                                     ,phase=ESMF_SINGLEPHASE            &
                                     ,rc=RC)
#endif
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
      IF(RC_REG==ESMF_SUCCESS)THEN
!       WRITE(0,*)' EARTH_REGISTER succeeded'
      ELSE
        WRITE(0,*)' EARTH_REGISTER failed  RC_REG=',RC_REG
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE EARTH_REGISTER
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE EARTH_INITIALIZE(EARTH_GRID_COMP                       &
                                 ,IMP_STATE                             &
                                 ,EXP_STATE                             &
                                 ,CLOCK_NEMS                            &
                                 ,RC_INIT)
!
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(ESMF_GridComp)               :: EARTH_GRID_COMP                 !<-- The EARTH component
      TYPE(ESMF_State)                  :: IMP_STATE                       !<-- The EARTH import state
      TYPE(ESMF_State)                  :: EXP_STATE                       !<-- The EARTH export state
      TYPE(ESMF_Clock)                  :: CLOCK_NEMS                      !<-- The Clock of the NEMS component
      INTEGER            ,INTENT(OUT)   :: RC_INIT                         !<-- Error return code
!
!-----------------------------------------------------------------------
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER :: RC
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Allocate the EARTH component's internal state, point at it,
!***  and attach it to the EARTH component.
!-----------------------------------------------------------------------
!
      ALLOCATE(EARTH_INT_STATE,stat=RC)
      wrap%EARTH_INT_STATE=>EARTH_INT_STATE
!
      CALL ESMF_GridCompSetInternalState(EARTH_GRID_COMP                &  !<--The EARTH component
                                        ,WRAP                           &  !<-- Pointer to the EARTH internal state
                                        ,RC)     
!
!-----------------------------------------------------------------------
!***  For the moment, use a direct copy of the NEMS Clock within
!***  the EARTH component.
!-----------------------------------------------------------------------
!
      earth_int_state%CLOCK_EARTH=CLOCK_NEMS
!
!-----------------------------------------------------------------------
!***  The ATM (atmosphere) gridded component resides inside of
!***  the EARTH internal state.
!-----------------------------------------------------------------------
!
      earth_int_state%ATM_GRID_COMP=ESMF_GridCompCreate(name        ="ATM component" &
                                                       ,rc          =RC)
!
!-----------------------------------------------------------------------
!***  Register the Initialize, Run, and Finalize routines of
!***  the ATM component.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Register ATM Init, Run, Finalize"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifdef ESMF_3
      CALL ESMF_GridCompSetServices(earth_int_state%ATM_GRID_COMP       &
                                   ,ATM_REGISTER                        &  !<-- The user's subroutine name
                                   ,RC)
#else
      CALL ESMF_GridCompSetServices(earth_int_state%ATM_GRID_COMP       &
                                   ,ATM_REGISTER                        &  !<-- The user's subroutine name
                                   ,rc=RC)
#endif
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Create the ATM import and export states.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create the ATM import state"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifdef ESMF_520rbs
      earth_int_state%ATM_IMP_STATE=ESMF_StateCreate(     NAME="ATM Import"      &
                                                    ,statetype=ESMF_STATE_IMPORT &
                                                    ,rc       =RC)
#else
      earth_int_state%ATM_IMP_STATE=ESMF_StateCreate(STATENAME="ATM Import"      &
                                                    ,statetype=ESMF_STATE_IMPORT &
                                                    ,rc       =RC)
#endif
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create the ATM export state"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifdef ESMF_520rbs
      earth_int_state%ATM_EXP_STATE=ESMF_StateCreate(     NAME="ATM Export"      &
                                                    ,statetype=ESMF_STATE_EXPORT &
                                                    ,rc       =RC)
#else
      earth_int_state%ATM_EXP_STATE=ESMF_StateCreate(STATENAME="ATM Export"      &
                                                    ,statetype=ESMF_STATE_EXPORT &
                                                    ,rc       =RC)
#endif
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Insert the import/export states of the ATMOS component into the
!***  import/export states of the EARTH component.  This simplifies
!***  the passing of information between lower and higher component 
!***  levels seen in the diagram above.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK= "Add the ATMOS states into the EARTH states"
!     CALL ESMF_LogWrite(MESSAGE_CHECK, ESMF_LOG_INFO, rc = RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateAdd(state       = IMP_STATE                        &
                        ,nestedState = earth_int_state%ATM_IMP_STATE    &
                        ,rc          = RC)
!
      CALL ESMF_StateAdd(state       = EXP_STATE                        &
                        ,nestedState = earth_int_state%ATM_EXP_STATE    &
                        ,rc          = RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Execute the Initialize step of the ATM component.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Execute the Initialize step of the ATM component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompInitialize(gridcomp   =earth_int_state%ATM_GRID_COMP &
                                  ,importState=earth_int_state%ATM_IMP_STATE &
                                  ,exportState=earth_int_state%ATM_EXP_STATE &
                                  ,clock      =earth_int_state%CLOCK_EARTH   &
                                  ,phase      =ESMF_SINGLEPHASE              &
                                  ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
      IF(RC_INIT==ESMF_SUCCESS)THEN
!       WRITE(0,*)' EARTH_INITIALIZE succeeded'
      ELSE
        WRITE(0,*)' EARTH_INITIALIZE failed  RC_INIT=',RC_INIT
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE EARTH_INITIALIZE
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE EARTH_RUN(EARTH_GRID_COMP                              &
                          ,IMP_STATE                                    &
                          ,EXP_STATE                                    &
                          ,CLOCK_NEMS                                   &
                          ,RC_RUN)
!
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(ESMF_GridComp)               :: EARTH_GRID_COMP                 !<-- The EARTH component
      TYPE(ESMF_State)                  :: IMP_STATE                       !<-- The EARTH import state
      TYPE(ESMF_State)                  :: EXP_STATE                       !<-- The EARTH export state
      TYPE(ESMF_Clock)                  :: CLOCK_NEMS                      !<-- The Clock of the NEMS component
      INTEGER            ,INTENT(OUT)   :: RC_RUN                          !<-- Error return code
!
!---------------------
!***  Local Variables
!---------------------
!
      TYPE(ESMF_Time) :: CURRTIME                                       &
                        ,STARTTIME
!
      TYPE(ESMF_TimeInterval) :: RUNDURATION
!
      INTEGER :: RC
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC    =ESMF_SUCCESS
      RC_RUN=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  Execute the Run step of the ATM component.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Execute the Run step of the  ATM component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompRun(gridcomp   =earth_int_state%ATM_GRID_COMP   &
                           ,importState=earth_int_state%ATM_IMP_STATE   &
                           ,exportState=earth_int_state%ATM_EXP_STATE   &
                           ,clock      =earth_int_state%CLOCK_EARTH     &
                           ,phase      =ESMF_SINGLEPHASE                &
                           ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Update the EARTH clock.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK = "Update the current time of the EARTH clock"
!     CALL ESMF_LogWrite(MESSAGE_CHECK, ESMF_LOG_INFO, rc = RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ClockGet(clock       = earth_int_state%CLOCK_EARTH      &
                        ,startTime   = startTime                        &
                        ,runDuration = runDuration                      &
                        ,rc          = RC)
!
      CURRTIME = STARTTIME + RUNDURATION
!
      CALL ESMF_ClockSet(clock    = earth_int_state%CLOCK_EARTH         &
                        ,currTime = CURRTIME                            &
                        ,rc       = RC)

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
      IF(RC_RUN==ESMF_SUCCESS)THEN
!       WRITE(0,*)' EARTH_RUN succeeded'
      ELSE
        WRITE(0,*)' EARTH_RUN failed  RC_RUN=',RC_RUN
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE EARTH_RUN
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE EARTH_FINALIZE(EARTH_GRID_COMP                         &
                               ,IMP_STATE                               &
                               ,EXP_STATE                               &
                               ,CLOCK_NEMS                              &
                               ,RC_FINALIZE)
!
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(ESMF_GridComp)               :: EARTH_GRID_COMP                 !<-- The EARTH component
      TYPE(ESMF_State)                  :: IMP_STATE                       !<-- The EARTH import state
      TYPE(ESMF_State)                  :: EXP_STATE                       !<-- The EARTH export state
      TYPE(ESMF_Clock)                  :: CLOCK_NEMS                      !<-- The Clock of the NEMS component
      INTEGER            ,INTENT(OUT)   :: RC_FINALIZE                     !<-- Error return code
!
!---------------------
!***  Local Variables
!---------------------
!
!-----------------------------------------------------------------------
!
      INTEGER :: RC
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Execute the Finalize step of the ATM ccomponent.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Execute the Finalize step of the  ATM component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompFinalize(gridcomp   =earth_int_state%ATM_GRID_COMP &
                                ,importState=earth_int_state%ATM_IMP_STATE &
                                ,exportState=earth_int_state%ATM_EXP_STATE &
                                ,clock      =earth_int_state%CLOCK_EARTH   &
                                ,phase      =ESMF_SINGLEPHASE              &
                                ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINALIZE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
      IF(RC_FINALIZE==ESMF_SUCCESS)THEN
!       WRITE(0,*)' EARTH_FINALIZE succeeded'
      ELSE
        WRITE(0,*)' EARTH_FINALIZE failed  RC_FINALIZE=',RC_FINALIZE
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE EARTH_FINALIZE
!
!-----------------------------------------------------------------------
!
      END MODULE module_EARTH_GRID_COMP
!
!-----------------------------------------------------------------------
