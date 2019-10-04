!TBH:  Matches r14147 of https://svnemc.ncep.noaa.gov/projects/nems/trunk/.  
!-----------------------------------------------------------------------
!
      MODULE module_EARTH_INTERNAL_STATE
!
!-----------------------------------------------------------------------
!***  Contents of the ESMF internal state of the EARTH component.
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
      PUBLIC :: EARTH_INTERNAL_STATE                                    &
               ,WRAP_EARTH_INTERNAL_STATE
!
!-----------------------------------------------------------------------
!
      TYPE EARTH_INTERNAL_STATE
!
        TYPE(ESMF_Clock   ) :: CLOCK_EARTH
!
        TYPE(ESMF_GridComp) :: ATM_GRID_COMP
        TYPE(ESMF_State   ) :: ATM_IMP_STATE
        TYPE(ESMF_State   ) :: ATM_EXP_STATE
!
        TYPE(ESMF_GridComp) :: OCEAN_GRID_COMP
        TYPE(ESMF_State   ) :: OCEAN_IMP_STATE
        TYPE(ESMF_State   ) :: OCEAN_EXP_STATE
!
        TYPE(ESMF_GridComp) :: ICE_GRID_COMP
        TYPE(ESMF_State   ) :: ICE_IMP_STATE
        TYPE(ESMF_State   ) :: ICE_EXP_STATE
!
      END TYPE EARTH_INTERNAL_STATE
!
!-----------------------------------------------------------------------
!
      TYPE WRAP_EARTH_INTERNAL_STATE
!
        TYPE(EARTH_INTERNAL_STATE),POINTER :: EARTH_INT_STATE
!
      END TYPE WRAP_EARTH_INTERNAL_STATE
!
!-----------------------------------------------------------------------
!
      END MODULE module_EARTH_INTERNAL_STATE
!
!-----------------------------------------------------------------------
