!-----------------------------------------------------------------------
!
      MODULE module_NMM_INTERNAL_STATE
!
!-----------------------------------------------------------------------
!***  The NMM component's ESMF internal state.
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
      PUBLIC :: NMM_INTERNAL_STATE                                      &
               ,WRAP_NMM_INTERNAL_STATE
!
!-----------------------------------------------------------------------
!
      TYPE NMM_INTERNAL_STATE
!
        INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE :: NPHS             &  !<-- Physics timestep
                                                      ,NTRACK           &  !<-- Storm locator flag
                                                      ,NUM_2WAY_CHILDREN   !<-- # of 2-way children in each domain
!
        INTEGER(kind=KINT),DIMENSION(:),POINTER :: COMM_MY_DOMAIN       &  !<-- MPI intracommunicator for all tasks on each domain
                                                  ,NPE_PRINT            &  !<-- Clocktime diagnostics from this MPI task
                                                  ,NUM_CHILDREN         &  !<-- How many children on each domain
                                                  ,P_C_TIME_RATIO          !<-- Ratio of a parent's DT to its child's
!
        REAL(kind=KFPT),DIMENSION(:),POINTER :: DT                         !<-- The fundamental timestep (s) of the domains
!
        LOGICAL(kind=KLOG),DIMENSION(:),POINTER :: MY_DOMAIN_MOVES      &  !<-- Does my domain move?
                                                  ,RESTARTED_RUN        &  !<-- Flag indicating if this is a restarted run
                                                  ,RST_OUT_00              !<-- Shall we write 00h history in restarted run?
!
        CHARACTER(len=5),DIMENSION(:),POINTER :: NEST_MODE                 !<-- Is the nesting 1-way or 2-way with the parent?
!
        LOGICAL(kind=KLOG),DIMENSION(:),POINTER :: I_AM_A_FCST_TASK     &  !<-- Is this task a fcst task on a given domain?
                                                  ,I_AM_A_NEST          &  !<-- Is this task on a nested domain?
                                                  ,I_AM_LEAD_FCST_TASK     !<-- Is this the lead fcst task on the domain?
!
        TYPE(ESMF_GridComp),DIMENSION(:),POINTER :: DOMAIN_GRID_COMP       !<-- Gridded components of all domains
!
        TYPE(ESMF_State),DIMENSION(:),POINTER :: IMP_STATE_DOMAIN          !<-- The import state of the DOMAIN components
        TYPE(ESMF_State),DIMENSION(:),POINTER :: EXP_STATE_DOMAIN          !<-- The export state of the DOMAIN components
!
        TYPE(ESMF_TimeInterval),DIMENSION(:),POINTER :: FILT_TIMESTEP      &  !<-- ESMF timestep for digital filter (s)
                                                       ,INTERVAL_CLOCKTIME &  !<-- ESMF time interval between clocktime prints (h)
                                                       ,INTERVAL_HISTORY   &  !<-- ESMF time interval between history output (h)
                                                       ,INTERVAL_RESTART   &  !<-- ESMF time interval between restart output (h)
                                                       ,TIMESTEP              !<-- The ESMF timestep (s)
!
        TYPE(ESMF_CplComp),DIMENSION(:),POINTER :: PC_CPL_COMP
!
        TYPE(ESMF_State),DIMENSION(:),POINTER :: IMP_STATE_PC_CPL          !<-- The import state of the P-C Coupler components
        TYPE(ESMF_State),DIMENSION(:),POINTER :: EXP_STATE_PC_CPL          !<-- The export state of the P-C Coupler components
!
      END TYPE NMM_INTERNAL_STATE
!
!-----------------------------------------------------------------------
!
      TYPE WRAP_NMM_INTERNAL_STATE
!
        TYPE(NMM_INTERNAL_STATE),POINTER :: NMM_INT_STATE
!
      END TYPE WRAP_NMM_INTERNAL_STATE
!
!-----------------------------------------------------------------------
!
      END MODULE module_NMM_INTERNAL_STATE
!
!-----------------------------------------------------------------------
