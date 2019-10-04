!---------------------------------------------------------------------------
!
      MODULE MODULE_DOMAIN_INTERNAL_STATE
!
!---------------------------------------------------------------------------
!***  Define all quantities that lie within the DOMAIN component's
!***  internal state.
!---------------------------------------------------------------------------
!
      USE ESMF
!
      USE module_KINDS
!
      USE MODULE_DERIVED_TYPES,ONLY: MIXED_DATA
!
!---------------------------------------------------------------------------
!
      IMPLICIT NONE
!
!---------------------------------------------------------------------------
!
      PRIVATE
!
      PUBLIC :: DOMAIN_INTERNAL_STATE                                       &
               ,WRAP_DOMAIN_INTERNAL_STATE
!
!---------------------------------------------------------------------------
!
      TYPE DOMAIN_INTERNAL_STATE
!
!---------------------------------------------------------------------------
!
        INTEGER(kind=KINT) :: KOUNT_TIMESTEPS                               &
                             ,SFC_FILE_RATIO                                   !<-- Ratio of upper parent's grid increment to this domain's
!
        INTEGER(kind=KINT) :: LEAD_TASK_DOMAIN                              &  !<-- The first task on a given domain
                             ,NUM_PES_FCST                                     !<-- The number of forecast tasks
!
        INTEGER(kind=KINT),DIMENSION(1:9) :: HANDLE_SEND_INTER_INT          &  !<-- For ISSends of intertask integer data after domain shifts
                                            ,HANDLE_SEND_INTER_REAL            !<-- For ISSends of intertask real data after domain shifts
!
        TYPE(ESMF_GridComp),ALLOCATABLE,DIMENSION(:) :: DOMAIN_CHILD_COMP      !<-- DOMAIN components of child domains
!
        TYPE(ESMF_GridComp) :: SOLVER_GRID_COMP                                !<-- The Solver gridded component
!
        TYPE(ESMF_State) :: IMP_STATE_SOLVER                                   !<-- The import state of the Solver component
        TYPE(ESMF_State) :: IMP_STATE_WRITE                                    !<-- The import state of the write components
!
        TYPE(ESMF_State) :: EXP_STATE_SOLVER                                   !<-- The export state of the Solver component
        TYPE(ESMF_State) :: EXP_STATE_WRITE                                    !<-- The export state of the write components
!
        TYPE(ESMF_Alarm) :: ALARM_HISTORY                                   &  !<-- The ESMF Alarm for history output
                           ,ALARM_RESTART                                   &  !<-- The ESMF Alarm for restart output
                           ,ALARM_CLOCKTIME                                    !<-- The ESMF Alarm for clocktime prints
!
        REAL(ESMF_KIND_R8) :: TIMESTEP_FINAL                                   !<-- The forecast's final timestep
!
        LOGICAL(kind=KLOG) :: ALLCLEAR_FROM_PARENT                          &  !<-- Child can proceed after parent is free
                             ,I_AM_A_NEST                                   &  !<-- Am I in a nested domain?
                             ,I_AM_A_PARENT                                 &  !<-- Am I in a parent domain?
                             ,MY_DOMAIN_MOVES                               &  !<-- Does this domain move?
                             ,RECV_ALL_CHILD_DATA                              !<-- Parent is free after all 2-way data recvd
!
        LOGICAL(kind=KLOG) :: FIRST_PASS                                    &  !<-- Note 1st time into NMM_INTEGRATE
                             ,RESTARTED_RUN                                 &  !<-- Is this a restarted forecast?
                             ,RESTARTED_RUN_FIRST                           &  !<-- Is is time for the initial output in a restarted run?
                             ,TS_INITIALIZED 
!
        CHARACTER(len=7) :: SFC_CONFLICT                                       !<-- Do/not search for nearest point with same sfc type
!
        TYPE(MIXED_DATA),DIMENSION(1:9) :: SHIFT_DATA                          !<-- Intertask shift data on the pre-move footprint
!
        TYPE(ESMF_FieldBundle) :: BUNDLE_NESTBC                                !<-- ESMF Bundle of BC update variables (parent to child)
!
        TYPE(ESMF_FieldBundle) :: BUNDLE_2WAY                                  !<-- ESMF Bundle of 2-way exchange vbls (child to parent)
!
        TYPE(ESMF_FieldBundle) :: MOVE_BUNDLE_H                             &  !<-- ESMF Bundle of update H variables on moving nests
                                 ,MOVE_BUNDLE_V                                !<-- ESMF Bundle of update V variables on moving nests
!
!---------------------------------------------------------------------------
!***  The following are specific to asynchronous quilting/writing.
!---------------------------------------------------------------------------
!
        LOGICAL(kind=KLOG) :: QUILTING                                      &  !<-- Is the user selecting asynchronous quilting/writing?
                             ,WRITE_LAST_RESTART                            &  !<-- Shall we write last restart file
                             ,WROTE_1ST_HIST                                   !<-- Has 1st history output been written?
!
        TYPE(ESMF_GridComp),DIMENSION(:),POINTER :: WRITE_COMPS                !<-- The array of Write gridded components
!
        INTEGER(kind=KINT) :: WRITE_GROUPS                                  &  !<-- The number of write groups
                             ,WRITE_GROUP_READY_TO_GO                       &  !<-- The active group of write tasks
                             ,WRITE_TASKS_PER_GROUP                            !<-- The number of write tasks in each write group
!
        INTEGER(kind=KINT),DIMENSION(:),POINTER :: LOCAL_ISTART,LOCAL_IEND  &  !<-- The local I limits of the forecast tasks
                                                  ,LOCAL_JSTART,LOCAL_JEND  &  !<-- The local J limits of the forecast tasks
                                                  ,PETLIST_FCST                !<-- Task ID list of fcst tasks on the domain
!
        INTEGER(kind=KINT),DIMENSION(:,:),POINTER :: PETLIST_WRITE             !<-- Task ID list of fcst tasks w/ write tasks by group
!
!---------------------------------------------------------------------------
!***  The following are specific to digital filtering.
!---------------------------------------------------------------------------
!
        INTEGER(kind=KINT) :: KSTEP,NSTEP
!
        INTEGER(kind=KINT) :: NUM_FIELDS_FILTER_2D                          &
                             ,NUM_FIELDS_FILTER_3D                          &
                             ,NUM_FIELDS_RESTORE_2D                         &
                             ,NUM_FIELDS_RESTORE_3D
!
        REAL(kind=KFPT) :: TOTALSUM
!
        REAL(kind=KFPT),DIMENSION(:,:,:),POINTER :: SAVE_2D,SAVE_2D_PHYS
        REAL(kind=KFPT),DIMENSION(:,:,:,:),POINTER :: SAVE_3D,SAVE_3D_PHYS
!
        REAL(kind=KFPT),DIMENSION(:),POINTER :: DOLPH_WGTS(:)
!
        LOGICAL(kind=KLOG) :: FIRST_FILTER
!
        TYPE(ESMF_FieldBundle) :: FILT_BUNDLE_FILTER                        &  !<-- ESMF Bundle of variables to filter
                                 ,FILT_BUNDLE_RESTORE                          !<-- ESMF Bundle of variables to restore to pre-filtered state
!
!---------------------------------------------------------------------------
!***  The following are used to transmit fields to/from an external
!***  ocean model via NUOPC coupling.
!---------------------------------------------------------------------------
!
        INTEGER(kind=KINT) :: KOUNT_NPRECIP,KOUNT_NPHS
!
!---------------------
!***  Imported fields
!---------------------
!
        REAL(ESMF_KIND_R8),DIMENSION(:,:),POINTER :: SST_COUPLED                  !<-- SST (K)
!
!---------------------
!***  Exported fields
!---------------------
!
        REAL(ESMF_KIND_R8),DIMENSION(:,:),POINTER :: INST_SFC_PRESSURE_COUPLED    !<-- Instantaneous sfc pressure (Pa)
        REAL(ESMF_KIND_R8),DIMENSION(:,:),POINTER :: INST_SENS_HT_FLX_COUPLED     !<-- Instantaneous sensible heat flux (W m-2)
        REAL(ESMF_KIND_R8),DIMENSION(:,:),POINTER :: INST_LAT_HT_FLX_COUPLED      !<-- Instantaneous latent heat flux (W m-2)
        REAL(ESMF_KIND_R8),DIMENSION(:,:),POINTER :: INST_NET_LW_FLX_COUPLED      !<-- Instantaneous net longwave flux (W m-2)
        REAL(ESMF_KIND_R8),DIMENSION(:,:),POINTER :: INST_NET_SW_FLX_COUPLED      !<-- Instantaneous net shortwave flux (W m-2)
        REAL(ESMF_KIND_R8),DIMENSION(:,:),POINTER :: MEAN_ZONAL_MOM_FLX_COUPLED   !<-- Mean zonal mom flux (N m-2)
        REAL(ESMF_KIND_R8),DIMENSION(:,:),POINTER :: MEAN_MERID_MOM_FLX_COUPLED   !<-- Mean merid mom flux (N m-2)
        REAL(ESMF_KIND_R8),DIMENSION(:,:),POINTER :: MEAN_PREC_RATE_COUPLED       !<-- Mean precip rate (kg m-2 s-1)
!
!---------------------------------------------------------------------------
!
      END TYPE DOMAIN_INTERNAL_STATE
!
!---------------------------------------------------------------------------
!
!---------------------------------------------------------------------------
!***  This state is supported by C pointers but not by F90 pointers
!***  therefore we use this "WRAP".
!---------------------------------------------------------------------------
!
      TYPE WRAP_DOMAIN_INTERNAL_STATE
        TYPE(DOMAIN_INTERNAL_STATE),POINTER :: DOMAIN_INT_STATE
      END TYPE WRAP_DOMAIN_INTERNAL_STATE
!
!---------------------------------------------------------------------------
!
      END MODULE MODULE_DOMAIN_INTERNAL_STATE
!
!---------------------------------------------------------------------------
