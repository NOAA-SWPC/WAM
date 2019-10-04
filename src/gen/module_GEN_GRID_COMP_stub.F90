!#include "../../ESMFVersionDefine.h"

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
      TYPE(ESMF_GridComp)               :: GEN_GRID_COMP                   !<-- GEN gridded component
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
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!
      CALL ESMF_GridCompSetEntryPoint(GEN_GRID_COMP                     &
                                     ,ESMF_METHOD_INITIALIZE            &
                                     ,GEN_INITIALIZE                    &
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
      CALL ESMF_GridCompSetEntryPoint(GEN_GRID_COMP                     &
                                     ,ESMF_METHOD_RUN                   &
                                     ,GEN_RUN                           &
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
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
      CALL ESMF_GridCompSetEntryPoint(GEN_GRID_COMP                     &
                                     ,ESMF_METHOD_FINALIZE              &
                                     ,GEN_FINALIZE                      &
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
!       WRITE(0,*)' GEN_SET_SERVICES SUCCEEDED'
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
      TYPE(ESMF_GridComp)               :: GEN_GRID_COMP                   !<-- The GEN gridded component
      TYPE(ESMF_State)                  :: IMP_STATE                    &  !<-- The GEN component's import state
                                          ,EXP_STATE                       !<-- The GEN component's export state
      TYPE(ESMF_Clock)                  :: CLOCK_ATM                       !<-- The ESMF Clock from the ATM component.
      INTEGER            ,INTENT(OUT)   :: RC_INIT                         !<-- Return code for Initialize step

!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC_INIT=ESMF_SUCCESS       ! Error signal variable
!
!-----------------------------------------------------------------------
!***********************************************************************
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
      TYPE(ESMF_GridComp)               :: GEN_GRID_COMP                   !<-- The GEN gridded component
      TYPE(ESMF_State)                  :: IMP_STATE                       !<-- The GEN Run step's import
      TYPE(ESMF_State)                  :: EXP_STATE                       !<-- The GEN Run step's export
      TYPE(ESMF_Clock)                  :: CLOCK_ATM                       !<-- The ATM ESMF Clock
      INTEGER,            INTENT(OUT)   :: RC_RUN                          !<-- Return code for the Run step

!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC_RUN=ESMF_SUCCESS       ! Error signal variable
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
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
      TYPE(ESMF_GridComp)               :: GEN_GRID_COMP                   !<-- The ATM gridded component
      TYPE(ESMF_State)                  :: IMP_STATE                       !<-- The ATM finalize step's import state
      TYPE(ESMF_State)                  :: EXP_STATE                       !<-- The ATM finalize step's export state
      TYPE(ESMF_Clock)                  :: CLOCK_ATM                       !<-- The main ESMF Clock
      INTEGER,            INTENT(OUT)   :: RC_FINAL                        !<-- Return code for the Finalize step
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC_FINAL=ESMF_SUCCESS       ! Error signal variable
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      END SUBROUTINE GEN_FINALIZE
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      END MODULE MODULE_GEN_GRID_COMP
