!TBH:  Matches r14147 of https://svnemc.ncep.noaa.gov/projects/nems/trunk/.  
!-----------------------------------------------------------------------
!
      MODULE module_NEMS_INTERNAL_STATE
!
!-----------------------------------------------------------------------
!***  Contents of the ESMF internal state of the NEMS component.
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
      PUBLIC :: NEMS_INTERNAL_STATE                                     &
               ,WRAP_NEMS_INTERNAL_STATE
!
!-----------------------------------------------------------------------
!
      TYPE NEMS_INTERNAL_STATE
!
        REAL :: DUMMY1
!
      END TYPE NEMS_INTERNAL_STATE
!
!-----------------------------------------------------------------------
!
      TYPE WRAP_NEMS_INTERNAL_STATE
!
        REAL :: DUMMY2
!
        TYPE(NEMS_INTERNAL_STATE),POINTER :: NEMS_INT_STATE
!
      END TYPE WRAP_NEMS_INTERNAL_STATE
!
!-----------------------------------------------------------------------
!
      END MODULE module_NEMS_INTERNAL_STATE
!
!-----------------------------------------------------------------------
