!TBH:  Matches r14147 of https://svnemc.ncep.noaa.gov/projects/nems/trunk/.  
!-----------------------------------------------------------------------
!
      MODULE module_ATM_INTERNAL_STATE
!
!-----------------------------------------------------------------------
!
      USE ESMF_MOD
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      PRIVATE
!
      PUBLIC :: ATM_INTERNAL_STATE                                      &
               ,WRAP_ATM_INTERNAL_STATE
!
!-----------------------------------------------------------------------
!
      TYPE ATM_INTERNAL_STATE
!
        CHARACTER(16) :: CORE
!
        TYPE(ESMF_Clock   ) :: CLOCK_ATM
        TYPE(ESMF_GridComp) :: CORE_GRID_COMP
        TYPE(ESMF_State   ) :: CORE_IMP_STATE
        TYPE(ESMF_State   ) :: CORE_EXP_STATE
!
      END TYPE ATM_INTERNAL_STATE
!
!-----------------------------------------------------------------------
!
      TYPE WRAP_ATM_INTERNAL_STATE
!
        TYPE(ATM_INTERNAL_STATE),POINTER :: ATM_INT_STATE
!
      END TYPE WRAP_ATM_INTERNAL_STATE
!
!-----------------------------------------------------------------------
!
      END MODULE module_ATM_INTERNAL_STATE
!
!-----------------------------------------------------------------------

