!-----------------------------------------------------------------------
!
      MODULE MODULE_DOMAIN_TASK_SPECS
!
!-----------------------------------------------------------------------
!***  Store fundamental values related to MPI tasks on each domain.
!-----------------------------------------------------------------------
!
      USE module_KINDS
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
      PUBLIC :: DOMAIN_TASK_SPECS
!
!-----------------------------------------------------------------------
!
      TYPE DOMAIN_TASK_SPECS
!
        INTEGER(kind=KINT) :: PARENT_DOMAIN_ID                             !<-- The ID of the upper parent domain DOMAIN 
!
        INTEGER(kind=KINT) :: INPES                                     &  !<-- # of MPI tasks in I direction on this domain
                             ,JNPES                                        !<-- # of MPI tasks in J direction on this domain
!
        INTEGER(kind=KINT) :: NUM_PETS                                     !<-- # of compute tasks on this domain
!
        INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE :: PET_MAP             !<-- List of this domain's task IDs in terms of 
!
        INTEGER(kind=KINT),DIMENSION(1:2) :: INDX_MIN                   &  !<-- This task's minimum I,J on this subdomain.
                                            ,INDX_MAX                      !<-- This task's maximum I,J  on this subdomain.
!
!
        REAL(kind=KDBL),DIMENSION(:,:),ALLOCATABLE :: CELL_AREA         &  !<-- Area within each grid cell (m**2)
                                                     ,ROT_ANGLE         &  !<-- Rotation angle (radians) from NMMB grid to earth lat/lon
                                                     ,SEA_MASK             !<-- Sea mask (1->water; 0->land)
!
        REAL(kind=KDBL),DIMENSION(:,:), POINTER :: GLAT=>NULL()         &
                                                  ,GLON=>NULL()         &
                                                  ,VLAT=>NULL()         &
                                                  ,VLON=>NULL()

!       The following 6 pointers are reference to the memory in the Grid for this nest
!       Their values need to be updated on the nest processsors only when the nest moves

        INTEGER(kind=KINT),DIMENSION(:,:),POINTER :: NEST_SEA_MASK=>NULL()  !<-- Sea mask (1->water; 0->land)
        REAL(kind=KDBL)   ,DIMENSION(:,:),POINTER :: NEST_CELL_AREA=>NULL() !<-- Area within each grid cell (m**2)
!
        REAL(kind=KDBL),DIMENSION(:,:), POINTER :: NEST_GLAT=>NULL()    &
                                                  ,NEST_GLON=>NULL()    &
                                                  ,NEST_VLAT=>NULL()    &
                                                  ,NEST_VLON=>NULL()
!
        LOGICAL(kind=KLOG) :: TASK_ACTIVE                                  !<-- Is the current task active on this domain?
!
        TYPE(ESMF_GRID) :: GRID                                            !<-- The ESMF_GRID object associated with this domain
!
        TYPE(ESMF_ROUTEHANDLE) :: PARENT2SELF                           &  !<-- Regrid interpolators from upper parent to nests
                                 ,SELF2PARENT                              !<-- Regrid interpolators from nests to upper parent
!
        TYPE(ESMF_FIELD),ALLOCATABLE :: EXPORTFIELDSLIST(:)             &  !<-- The export Fields
                                       ,IMPORTFIELDSLIST(:)                !<-- The import Fields
!
      END TYPE DOMAIN_TASK_SPECS
!
!-----------------------------------------------------------------------
!
      END MODULE MODULE_DOMAIN_TASK_SPECS
!
!-----------------------------------------------------------------------
