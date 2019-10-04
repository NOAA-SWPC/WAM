! 05/11/2011   Weiyu Yang   Modified for using the ESMF 5.2.0r_beta_snapshot_07.
! 02/09/2012   Weiyu Yang   Modified for using the ESMF 5.2.0rp1 library.
!-------------------------------------------------------------------------------
!#include "../../ESMFVersionDefine.h"

      MODULE module_FIM_GRID_COMP

      USE ESMF

      IMPLICIT NONE

      PRIVATE
      PUBLIC :: FIM_REGISTER

      INTEGER :: DUMMY

      CONTAINS

!#######################################################################

      SUBROUTINE FIM_REGISTER(FIM_GRID_COMP,RC_REG)
      TYPE(ESMF_GridComp)               :: FIM_GRID_COMP
      INTEGER            ,INTENT(OUT)   :: RC_REG

      INTEGER :: RC

      write(0,*) "    FIM_REGISTER"


      CALL ESMF_GridCompSetEntryPoint(FIM_GRID_COMP ,ESMF_METHOD_INITIALIZE ,FIM_INITIALIZE ,rc=RC)
      CALL ESMF_GridCompSetEntryPoint(FIM_GRID_COMP ,ESMF_METHOD_RUN,        FIM_RUN,        rc=RC)
      CALL ESMF_GridCompSetEntryPoint(FIM_GRID_COMP ,ESMF_METHOD_FINALIZE,   FIM_FINALIZE   ,rc=RC)

      RC_REG = ESMF_SUCCESS
      write(0,*) "    END OF FIM_REGISTER"

      END SUBROUTINE FIM_REGISTER

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      SUBROUTINE FIM_INITIALIZE(FIM_GRID_COMP ,IMP_STATE ,EXP_STATE ,CLOCK_FIM ,RC_INIT)

      TYPE(ESMF_GridComp)               :: FIM_GRID_COMP
      TYPE(ESMF_State)                  :: IMP_STATE
      TYPE(ESMF_State)                  :: EXP_STATE
      TYPE(ESMF_Clock)                  :: CLOCK_FIM
      INTEGER            ,INTENT(OUT)   :: RC_INIT

      write(0,*) "        FIM_INITIALIZE stub"
      RC_INIT = ESMF_SUCCESS
      write(0,*) "        END OF FIM_INITIALIZE stub"

      END SUBROUTINE FIM_INITIALIZE

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      SUBROUTINE FIM_RUN(FIM_GRID_COMP ,IMP_STATE ,EXP_STATE ,CLOCK_FIM ,RC_RUN)

      TYPE(ESMF_GridComp)               :: FIM_GRID_COMP
      TYPE(ESMF_State)                  :: IMP_STATE
      TYPE(ESMF_State)                  :: EXP_STATE
      TYPE(ESMF_Clock)                  :: CLOCK_FIM
      INTEGER            ,INTENT(OUT)   :: RC_RUN

      write(0,*) "        FIM_RUN stub"
      RC_RUN=ESMF_SUCCESS
      write(0,*) "        END OF FIM_RUN stub"

      END SUBROUTINE FIM_RUN

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      SUBROUTINE FIM_FINALIZE(FIM_GRID_COMP ,IMP_STATE ,EXP_STATE ,CLOCK_FIM ,RC_FINALIZE)

      TYPE(ESMF_GridComp)               :: FIM_GRID_COMP
      TYPE(ESMF_State)                  :: IMP_STATE
      TYPE(ESMF_State)                  :: EXP_STATE
      TYPE(ESMF_Clock)                  :: CLOCK_FIM
      INTEGER            ,INTENT(OUT)   :: RC_FINALIZE

      write(0,*) "        FIM_FINALIZE stub"
      RC_FINALIZE=ESMF_SUCCESS
      write(0,*) "        END OF FIM_FINALIZE stub"

      END SUBROUTINE FIM_FINALIZE

!#######################################################################

      END MODULE module_FIM_GRID_COMP
