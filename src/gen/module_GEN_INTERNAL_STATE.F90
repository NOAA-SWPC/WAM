!#include "../../ESMFVersionDefine.h"

!-----------------------------------------------------------------------
!
      MODULE module_GEN_INTERNAL_STATE
!
!-----------------------------------------------------------------------
!
      USE ESMF
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      PRIVATE
!
      PUBLIC :: GEN_INTERNAL_STATE                                    &
               ,WRAP_GEN_INTERNAL_STATE
!
!-----------------------------------------------------------------------
!
      TYPE GEN_INTERNAL_STATE
          TYPE(ESMF_GridComp) :: GC_GEN_DYN                           &
                                ,GC_GEN_PHY

          TYPE(ESMF_State   ) :: IMP_GEN_DYN,EXP_GEN_DYN              &  !<-- Import/export states for GEN Dynamics
                                ,IMP_GEN_PHY,EXP_GEN_PHY                 !<-- Import/export states for GEN Physics

          TYPE(ESMF_Clock   ) :: CLOCK_GEN


          TYPE(ESMF_CplComp)  :: GC_GEN_CPL
 
      END TYPE GEN_INTERNAL_STATE
!
!-----------------------------------------------------------------------
!
      TYPE WRAP_GEN_INTERNAL_STATE
!
        TYPE(GEN_INTERNAL_STATE),POINTER :: GEN_INT_STATE
!
      END TYPE WRAP_GEN_INTERNAL_STATE
!
!-----------------------------------------------------------------------
!
      END MODULE module_GEN_INTERNAL_STATE
!
!-----------------------------------------------------------------------

