!-----------------------------------------------------------------------
!
      MODULE MODULE_DERIVED_TYPES
!
!-----------------------------------------------------------------------
!
!***  This module contains various derived datatypes used in
!***  the NMM-B nesting.
!
!-----------------------------------------------------------------------
!
! PROGRAM HISTORY LOG:
!
!   2013-11-08  Black - Created
!
!-----------------------------------------------------------------------
!
! USAGE: 
!
!-----------------------------------------------------------------------
!
      USE ESMF
!
      USE module_KINDS
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      PRIVATE
!
      PUBLIC :: BC_H                                                    &
               ,BC_V                                                    &
               ,BC_H_ALL                                                &
               ,BC_V_ALL                                                &
               ,BNDS_2D                                                 &
               ,CHILD_UPDATE_LINK                                       &
               ,COMMS_FAMILY                                            &
               ,CTASK_LIMITS                                            &
               ,DOMAIN_DATA                                             &
               ,FILT_4D                                                 &
               ,HANDLE_CHILD_LIMITS                                     &
               ,HANDLE_CHILD_TOPO_S                                     &
               ,HANDLE_CHILD_TOPO_N                                     &
               ,HANDLE_CHILD_TOPO_W                                     &
               ,HANDLE_CHILD_TOPO_E                                     &
               ,HANDLE_I_SW                                             &
               ,HANDLE_J_SW                                             &
               ,HANDLE_PACKET_S_H                                       &
               ,HANDLE_PACKET_S_V                                       &
               ,HANDLE_PACKET_N_H                                       &
               ,HANDLE_PACKET_N_V                                       &
               ,HANDLE_PACKET_W_H                                       &
               ,HANDLE_PACKET_W_V                                       &
               ,HANDLE_PACKET_E_H                                       &
               ,HANDLE_PACKET_E_V                                       &
               ,HANDLE_PARENT_DOM_LIMITS                                &
               ,HANDLE_PARENT_ITE                                       &
               ,HANDLE_PARENT_ITS                                       &
               ,HANDLE_PARENT_JTE                                       &
               ,HANDLE_PARENT_JTS                                       &
               ,INFO_SEND                                               &
               ,INTEGER_DATA                                            &
               ,INTEGER_DATA_2D                                         &
               ,INTERIOR_DATA_FROM_PARENT                               &
               ,MIXED_DATA                                              &
               ,MIXED_DATA_TASKS                                        &
               ,MULTIDATA                                               &
               ,PTASK_LIMITS                                            &
               ,REAL_DATA                                               &
               ,REAL_DATA_2D                                            &
               ,REAL_DATA_TASKS                                         &
               ,REAL_VBLS_3D
!
!-----------------------------------------------------------------------
!
      TYPE MIXED_DATA
        INTEGER(kind=KINT),DIMENSION(:),POINTER :: DATA_INTEGER
        REAL(kind=KFPT),DIMENSION(:),POINTER :: DATA_REAL
      END TYPE MIXED_DATA
!
      TYPE MIXED_DATA_TASKS
        TYPE(MIXED_DATA),DIMENSION(:),POINTER :: TASKS
      END TYPE MIXED_DATA_TASKS
!
      TYPE INTEGER_DATA
        INTEGER(kind=KINT),DIMENSION(:),POINTER :: DATA
      END TYPE INTEGER_DATA
!
      TYPE INTEGER_DATA_2D
        INTEGER(kind=KINT),DIMENSION(:,:),POINTER :: DATA
      END TYPE INTEGER_DATA_2D
!
      TYPE REAL_DATA
        REAL(kind=KFPT),DIMENSION(:),POINTER :: DATA
      END TYPE REAL_DATA
!
      TYPE REAL_DATA_2D
        REAL(kind=KFPT),DIMENSION(:,:),POINTER :: DATA
      END TYPE REAL_DATA_2D
!
      TYPE REAL_VBLS_3D
        REAL(kind=KFPT),DIMENSION(:,:,:),POINTER :: VBL
      END TYPE REAL_VBLS_3D
!
      TYPE REAL_DATA_TASKS
        TYPE(REAL_DATA),DIMENSION(:),POINTER :: TASKS
      END TYPE REAL_DATA_TASKS
!
      TYPE MULTIDATA
        TYPE(REAL_DATA_TASKS),DIMENSION(:),POINTER :: CHILD
      END TYPE MULTIDATA
!
      TYPE BNDS_2D
        INTEGER(kind=KINT) :: LBND1
        INTEGER(kind=KINT) :: UBND1
        INTEGER(kind=KINT) :: LBND2
        INTEGER(kind=KINT) :: UBND2
      END TYPE BNDS_2D
!
      TYPE :: INTERIOR_DATA_FROM_PARENT
        INTEGER(kind=KINT) :: ID 
        INTEGER(kind=KINT) :: NPTS
        INTEGER(kind=KINT),DIMENSION(1:2) :: ISTART
        INTEGER(kind=KINT),DIMENSION(1:2) :: IEND
        INTEGER(kind=KINT),DIMENSION(1:2) :: JSTART
        INTEGER(kind=KINT),DIMENSION(1:2) :: JEND
      END TYPE INTERIOR_DATA_FROM_PARENT
!
      TYPE :: CHILD_UPDATE_LINK
        INTEGER(kind=KINT),POINTER :: TASK_ID
        INTEGER(kind=KINT),POINTER :: NUM_PTS_UPDATE_HZ
        INTEGER(kind=KINT),DIMENSION(:),POINTER :: IL
        INTEGER(kind=KINT),DIMENSION(:),POINTER :: JL
        TYPE(CHILD_UPDATE_LINK),POINTER :: NEXT_LINK
      END TYPE CHILD_UPDATE_LINK
!
      TYPE :: COMMS_FAMILY
        INTEGER(kind=KINT) :: TO_PARENT
        INTEGER(kind=KINT) :: TO_FCST_TASKS
        INTEGER(kind=KINT),DIMENSION(:),POINTER :: TO_CHILDREN
      END TYPE COMMS_FAMILY
!
      TYPE :: DOMAIN_DATA
        TYPE(INTEGER_DATA),DIMENSION(:),POINTER :: CHILDREN
      END TYPE DOMAIN_DATA
!
      TYPE :: DOMAIN_DATA_2
        TYPE(INTEGER_DATA_2D),DIMENSION(:),POINTER :: CHILDREN
      END TYPE DOMAIN_DATA_2
!
      TYPE :: DOMAIN_LIMITS
        INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE :: ITS
        INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE :: ITE
        INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE :: JTS
        INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE :: JTE
      END TYPE DOMAIN_LIMITS
!
      TYPE :: TASK_LIMITS
        INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE :: ITS
        INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE :: ITE
        INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE :: JTS
        INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE :: JTE
      END TYPE TASK_LIMITS
!
      TYPE :: BC_INFO
        TYPE(CHILD_INFO),DIMENSION(:),POINTER :: CHILDREN
      END TYPE BC_INFO
!
      TYPE CHILD_INFO
        INTEGER(kind=KINT),DIMENSION(:,:,:),POINTER :: INFO
      END TYPE CHILD_INFO
!
      TYPE BC_2D
        REAL(kind=KFPT),DIMENSION(:,:),ALLOCATABLE :: SIDE
      END TYPE BC_2D
!
      TYPE BC_3D
        REAL(kind=KFPT),DIMENSION(:,:,:),ALLOCATABLE :: SIDE
      END TYPE BC_3D
!
      TYPE BC_4D
        REAL(kind=KFPT),DIMENSION(:,:,:,:),ALLOCATABLE :: SIDE
      END TYPE BC_4D
!
      TYPE BC_H
        TYPE(BC_2D),DIMENSION(:),ALLOCATABLE :: VAR_2D
        TYPE(BC_3D),DIMENSION(:),ALLOCATABLE :: VAR_3D
        TYPE(BC_4D),DIMENSION(:),ALLOCATABLE :: VAR_4D
      END TYPE BC_H
!
      TYPE BC_V
        TYPE(BC_2D),DIMENSION(:),ALLOCATABLE :: VAR_2D
        TYPE(BC_3D),DIMENSION(:),ALLOCATABLE :: VAR_3D
      END TYPE BC_V
!
      TYPE BC_2D_ALL
        REAL(kind=KFPT),DIMENSION(:,:,:),POINTER :: SOUTH     &
                                                   ,NORTH     &
                                                   ,WEST      &
                                                   ,EAST
        REAL(kind=KFPT),DIMENSION(:,:),POINTER :: FULL_VAR
      END TYPE BC_2D_ALL
!
      TYPE BC_3D_ALL
        REAL(kind=KFPT),DIMENSION(:,:,:,:),POINTER :: SOUTH   &
                                                     ,NORTH   &
                                                     ,WEST    &
                                                     ,EAST
        REAL(kind=KFPT),DIMENSION(:,:,:),POINTER :: FULL_VAR
      END TYPE BC_3D_ALL
!
      TYPE BC_4D_ALL
        REAL(kind=KFPT),DIMENSION(:,:,:,:,:),POINTER :: SOUTH &
                                                       ,NORTH &
                                                       ,WEST  &
                                                       ,EAST
        REAL(kind=KFPT),DIMENSION(:,:,:,:),POINTER :: FULL_VAR
      END TYPE BC_4D_ALL
!
      TYPE BC_H_ALL
        TYPE(BC_2D_ALL),DIMENSION(:),ALLOCATABLE :: VAR_2D
        TYPE(BC_3D_ALL),DIMENSION(:),ALLOCATABLE :: VAR_3D
        TYPE(BC_4D_ALL),DIMENSION(:),ALLOCATABLE :: VAR_4D
      END TYPE BC_H_ALL
!
      TYPE BC_V_ALL
        TYPE(BC_2D_ALL),DIMENSION(:),ALLOCATABLE :: VAR_2D
        TYPE(BC_3D_ALL),DIMENSION(:),ALLOCATABLE :: VAR_3D
      END TYPE BC_V_ALL
!
      TYPE FILT_4D
        REAL(kind=KFPT),DIMENSION(:,:,:,:,:),ALLOCATABLE :: BASE
      END TYPE FILT_4D
!
!-----------------------------------------------------------------------
!
      INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE :: HANDLE_I_SW        &
                                                    ,HANDLE_J_SW
!
      TYPE(INTEGER_DATA),DIMENSION(:),ALLOCATABLE :: HANDLE_PARENT_DOM_LIMITS  !<-- Request handle for ISSend of parent domain limits
!
      TYPE(INTEGER_DATA),DIMENSION(:),POINTER :: HANDLE_PARENT_ITE      &   !<-- Request handles for ISSends
                                                ,HANDLE_PARENT_ITS      &   !    of each parent task's
                                                ,HANDLE_PARENT_JTE      &   !<-- integration limits to children.
                                                ,HANDLE_PARENT_JTS
!
      TYPE(CHILD_UPDATE_LINK),POINTER,SAVE :: TAIL
!
      TYPE(DOMAIN_DATA),DIMENSION(:),POINTER,SAVE :: HANDLE_CHILD_LIMITS &  !<-- Request handles for parents' IRecvs if child task limits
!
                                                    ,HANDLE_CHILD_TOPO_S &  !<-- Request handles for parents' IRecvs of child bndry topo
                                                    ,HANDLE_CHILD_TOPO_N &  !
                                                    ,HANDLE_CHILD_TOPO_W &  !
                                                    ,HANDLE_CHILD_TOPO_E    !<--
!
      TYPE(DOMAIN_DATA),DIMENSION(:),POINTER,SAVE :: HANDLE_PACKET_S_H   &  !<-- Request handles for parents' ISends of bndry info packets
                                                    ,HANDLE_PACKET_S_V   &  !
                                                    ,HANDLE_PACKET_N_H   &  !
                                                    ,HANDLE_PACKET_N_V   &  !
                                                    ,HANDLE_PACKET_W_H   &  !
                                                    ,HANDLE_PACKET_W_V   &  !
                                                    ,HANDLE_PACKET_E_H   &  !
                                                    ,HANDLE_PACKET_E_V      !<--
!
      TYPE(TASK_LIMITS),DIMENSION(:),ALLOCATABLE,SAVE :: PTASK_LIMITS      !<-- I,J limits on parent task subdomains
!
      TYPE(DOMAIN_DATA_2),DIMENSION(:),POINTER,SAVE :: CTASK_LIMITS        !<-- For limits of parents' children's tasks' subdomains
!
      TYPE(BC_INFO),DIMENSION(:),POINTER,SAVE :: INFO_SEND                 !<-- Parent info to children about which BC updates
!
!-----------------------------------------------------------------------
!
      END MODULE MODULE_DERIVED_TYPES
!
!-----------------------------------------------------------------------
