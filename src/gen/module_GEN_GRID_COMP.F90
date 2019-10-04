#include "../../ESMFVersionDefine.h"

!----------------------------------------------------------------------
!
      MODULE MODULE_GEN_GRID_COMP
!
!-----------------------------------------------------------------------
!
!***  THIS IS THE GEN GRIDDED COMPONENT MODULE.
!***  IT WILL SET UP DYNAMICS, PHYSICS, AND COUPLER SUBCOMPONENTS
!***  AND RUN THEIR INITIALIZE, RUN, AND FINALIZE ROUTINES.
!
!-----------------------------------------------------------------------
!
! PROGRAM HISTORY LOG:
!   2010-11-30  W Yang - Add the "Generic Core".
!   2011-05-11  W Yang - Modified for using the ESMF 5.2.0r_beta_snapshot_07.
!   2011-09-27  W Yang - Modified for using the ESMF 5.2.0r library.
!
! USAGE: GEN Gridded component parts called from subroutines within
!        module_ATM_GRID_COMP.F90.
!
!-----------------------------------------------------------------------
!
      USE ESMF
!
      USE MODULE_GEN_INTERNAL_STATE,ONLY: GEN_INTERNAL_STATE            &
                                         ,WRAP_GEN_INTERNAL_STATE

      USE MODULE_ERR_MSG,ONLY: ERR_MSG,MESSAGE_CHECK
!
      IMPLICIT NONE
!
      PRIVATE
!
      PUBLIC :: GEN_REGISTER 
!
!-----------------------------------------------------------------------
!
      TYPE(GEN_INTERNAL_STATE), POINTER, SAVE :: gen_int_state
      TYPE(WRAP_GEN_INTERNAL_STATE),     SAVE :: WRAP
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE GEN_REGISTER(GEN_GRID_COMP,RC_REG)
! 
!-----------------------------------------------------------------------
!***  Register the GEN gridded component's initialize, run, and finalize
!***  routines.
!-----------------------------------------------------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: GEN_GRID_COMP                   !<-- GEN gridded component
!
      INTEGER,INTENT(OUT) :: RC_REG                                        !<-- Return code for register
!     
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      INTEGER :: RC
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC=ESMF_SUCCESS       ! Error signal variable
!
!-----------------------------------------------------------------------
!***  Register the GEN INITIALIZE subroutine.  Since it is just one 
!***  subroutine, use ESMF_SINGLEPHASE.  The second argument is
!***  a pre-defined subroutine type, such as ESMF_SETINIT, ESMF_SETRUN, 
!***  or ESMF_SETFINAL.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create/Load the Configure Object"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!
      CALL ESMF_GridCompSetEntryPoint(GEN_GRID_COMP                     &  !<-- GEN gridded component
                                     ,ESMF_METHOD_INITIALIZE            &  !<-- Subroutine type
                                     ,GEN_INITIALIZE                    &  !<-- User's subroutine name
                                     ,phase=ESMF_SINGLEPHASE            &
                                     ,rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
!-----------------------------------------------------------------------
!***  Register the Run step of the GEN component.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
        MESSAGE_CHECK="Set 1st Entry Point for GEN Run"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
        CALL ESMF_GridCompSetEntryPoint(GEN_GRID_COMP                     &  !<-- GEN gridded component
                                       ,ESMF_METHOD_RUN                   &  !<-- Subroutine type
                                       ,GEN_RUN                           &  !<-- The primary Dynamics / Physics /Coupler sequence
                                       ,phase=ESMF_SINGLEPHASE            &
                                       ,rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
!
!-----------------------------------------------------------------------
!***  Register the GEN FINALIZE subroutine.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
      MESSAGE_CHECK="Set Entry Point for GEN Finalize"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
      CALL ESMF_GridCompSetEntryPoint(GEN_GRID_COMP                     &  !<-- GEN gridded component
                                     ,ESMF_METHOD_FINALIZE              &  !<-- Subroutine type
                                     ,GEN_FINALIZE                      &  !<-- User's subroutine name
                                     ,phase=ESMF_SINGLEPHASE            &
                                     ,rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
!-----------------------------------------------------------------------
!***  Check the error signal variable and print out the result.
!-----------------------------------------------------------------------
!
      IF(RC_REG==ESMF_SUCCESS)THEN
        WRITE(0,*)' GEN_SET_SERVICES SUCCEEDED'
      ELSE
        WRITE(0,*)' GEN_SET_SERVICES FAILED  RC_REG=',RC_REG
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE GEN_REGISTER
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE GEN_INITIALIZE(GEN_GRID_COMP                           &
                               ,IMP_STATE                               &
                               ,EXP_STATE                               &
                               ,CLOCK_ATM                               &
                               ,RC_INIT)
!
!------------------------
!***  Argument variables
!------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: GEN_GRID_COMP                   !<-- The GEN gridded component
      TYPE(ESMF_State)   ,INTENT(INOUT) :: IMP_STATE                    &  !<-- The GEN component's import state
                                          ,EXP_STATE                       !<-- The GEN component's export state
      TYPE(ESMF_Clock)   ,INTENT(INOUT) :: CLOCK_ATM                       !<-- The ESMF Clock from the ATM component.
      INTEGER            ,INTENT(OUT)   :: RC_INIT                         !<-- Return code for Initialize step
!
!---------------------
!***  Local variables
!---------------------

      INTEGER       :: IERR
      INTEGER       :: RC
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC     =ESMF_SUCCESS
      RC_INIT=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  Allocate the GEN component's internal state, point at it,
!***  and attach it to the GEN component.
!-----------------------------------------------------------------------
!
      ALLOCATE(gen_int_state, stat = RC)
      wrap%GEN_INT_STATE => gen_int_state
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set the GEN Internal State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompSetInternalState(GEN_GRID_COMP                  &
                                        ,WRAP                           &
                                        ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  For the moment, use a direct copy of the ATM Clock.
!-----------------------------------------------------------------------
!
      gen_int_state%CLOCK_GEN = CLOCK_ATM

! for testing.
!-------------
      PRINT*, 'Test run in the gen initialize routine.'

!-----------------------------------------------------------------------
!***  WRITE THE FINAL ERROR SIGNAL.
!-----------------------------------------------------------------------
!
      IF(RC_INIT==ESMF_SUCCESS)THEN
        WRITE(0,*)'GEN INITIALIZE STEP SUCCEEDED'
      ELSE
        WRITE(0,*)'GEN INITIALIZE STEP FAILED RC_INIT=',RC_INIT
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE GEN_INITIALIZE
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE GEN_RUN(GEN_GRID_COMP                                  &
                        ,IMP_STATE                                      &
                        ,EXP_STATE                                      &
                        ,CLOCK_ATM                                      &
                        ,RC_RUN)
!
!-----------------------------------------------------------------------
!***  RUN THE GEN GRIDDED COMPONENT.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument variables
!------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: GEN_GRID_COMP                   !<-- The GEN gridded component
      TYPE(ESMF_State),   INTENT(INOUT) :: IMP_STATE                       !<-- The GEN Run step's import
      TYPE(ESMF_State),   INTENT(INOUT) :: EXP_STATE                       !<-- The GEN Run step's export
      TYPE(ESMF_Clock),   INTENT(INOUT) :: CLOCK_ATM                       !<-- The ATM ESMF Clock
      INTEGER,            INTENT(OUT)   :: RC_RUN                          !<-- Return code for the Run step
!
!---------------------
!***  Local variables
!---------------------
!
      INTEGER            :: RC                                             !<-- Error signal variables.

      RC     = ESMF_SUCCESS
      RC_RUN = ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  For the moment, use a direct copy of the ATM Clock.
!-----------------------------------------------------------------------
!
      gen_int_state%CLOCK_GEN = CLOCK_ATM

! for testing.
!-------------
      PRINT*, 'Test run in the gen run routine.'

!
!-----------------------------------------------------------------------
!***  WRITE THE FINAL ERROR SIGNAL.
!-----------------------------------------------------------------------
!
      IF(RC_RUN == ESMF_SUCCESS) THEN
          WRITE(0, *) 'GEN ATM RUN STEP SUCCEEDED'
      ELSE
          WRITE(0, *) 'GEN ATM RUN STEP FAILED RC_RUN=', RC_RUN
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE GEN_RUN
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE GEN_FINALIZE(GEN_GRID_COMP                             &
                             ,IMP_STATE                                 &
                             ,EXP_STATE                                 &
                             ,CLOCK_ATM                                 &
                             ,RC_FINAL)
!
!-----------------------------------------------------------------------
!***  THIS ROUTINE FINALIZES THE ATM GRIDDED COMPONENT.
!-----------------------------------------------------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: GEN_GRID_COMP                   !<-- The ATM gridded component
      TYPE(ESMF_State),   INTENT(INOUT) :: IMP_STATE                       !<-- The ATM finalize step's import state
      TYPE(ESMF_State),   INTENT(INOUT) :: EXP_STATE                       !<-- The ATM finalize step's export state
      TYPE(ESMF_Clock),   INTENT(INOUT) :: CLOCK_ATM                       !<-- The main ESMF Clock
      TYPE(ESMF_Config)                 :: CF                              !<-- The config object
      INTEGER,            INTENT(OUT)   :: RC_FINAL                        !<-- Return code for the Finalize step
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      INTEGER :: RC                                                        ! The final error signal variables.
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC = ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  For the moment, use a direct copy of the ATM Clock.
!-----------------------------------------------------------------------
!
      gen_int_state%CLOCK_GEN = CLOCK_ATM

!-----------------------------------------------------------------------
!***  RETRIEVE THE CONFIG OBJECT CF FROM THE GEN GRIDDED COMPONENT.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Retrieve Config Object from GEN Component"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompGet(gridcomp=GEN_GRID_COMP                      &  !<-- The GEN gridded component
                           ,config  =CF                                 &  !<-- The config object (~namelist)
                           ,rc      =RC)

! This line is needed to keep for the GEN regression test.
!---------------------------------------------------------
      PRINT*, 'Test run in the gen finalize routine.'

!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  THE FINAL ERROR SIGNAL INFORMATION.
!-----------------------------------------------------------------------
!
      IF(RC_FINAL==ESMF_SUCCESS)THEN
        WRITE(0,*)'GEN FINALIZE STEP SUCCEEDED'
      ELSE
        WRITE(0,*)'GEN FINALIZE STEP FAILED'
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE GEN_FINALIZE
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      END MODULE MODULE_GEN_GRID_COMP
