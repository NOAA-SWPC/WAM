!-----------------------------------------------------------------------
!
      MODULE module_NMM_GRID_COMP
!
!-----------------------------------------------------------------------
!***  This is the NMM-B module.  It will set up one or more Domain
!***  subcomponents then execute their Initialize, Run, and Finalize
!***  steps.
!-----------------------------------------------------------------------
!
! PROGRAM HISTORY LOG:
!   2011-02  W. Yang - Updated to use both the ESMF 4.0.0rp2 library,
!                      ESMF 5 library and the the ESMF 3.1.0rp2 library.
!   2011-05  W. Yang - Modified for using the ESMF 5.2.0r_beta_snapshot_07.
!   2011-07    Black - Modified for moving nests.
!   2011-09  W. Yang - Modified for using the ESMF 5.2.0r library.
!   2012-07    Black - Modified for 'generational' task usage.
!   2015-09  Fei Liu - Modified for use with NUOPC driver.
!-----------------------------------------------------------------------
!
      USE MPI
      USE ESMF
      USE NUOPC
!
      USE module_KINDS
!
      USE module_DOMAIN_NUOPC_SET,ONLY: DOMAIN_DESCRIPTORS              &
                                       ,I_AM_PET                        &
                                       ,I_AM_ROOT                       &
                                       ,NMMB_CreateDomainFields         &
                                       ,NMMB_CreateRouteHandle          &
                                       ,NMMB_GridCreate 
!
      USE module_NMM_INTERNAL_STATE,ONLY: NMM_INTERNAL_STATE            &
                                         ,WRAP_NMM_INTERNAL_STATE
!
      USE module_DOMAIN_GRID_COMP,ONLY: DOMAIN_REGISTER                    !<-- The Register routine for DOMAIN_GRID_COMP
!
      USE module_DOMAIN_INTERNAL_STATE,ONLY: DOMAIN_INTERNAL_STATE      &
                                            ,WRAP_DOMAIN_INTERNAL_STATE
!
      USE module_DOMAIN_TASK_SPECS,ONLY: DOMAIN_TASK_SPECS
!
      USE module_NMM_INTEGRATE,ONLY: NMM_INTEGRATE
!
      USE module_DERIVED_TYPES,ONLY: COMMS_FAMILY                       &
                                    ,CTASK_LIMITS                       &
                                    ,HANDLE_CHILD_LIMITS                &
                                    ,HANDLE_CHILD_TOPO_S                &
                                    ,HANDLE_CHILD_TOPO_N                &
                                    ,HANDLE_CHILD_TOPO_W                &
                                    ,HANDLE_CHILD_TOPO_E                &
                                    ,HANDLE_I_SW                        &
                                    ,HANDLE_J_SW                        &
                                    ,HANDLE_PACKET_S_H                  &
                                    ,HANDLE_PACKET_S_V                  &
                                    ,HANDLE_PACKET_N_H                  &
                                    ,HANDLE_PACKET_N_V                  &
                                    ,HANDLE_PACKET_W_H                  &
                                    ,HANDLE_PACKET_W_V                  &
                                    ,HANDLE_PACKET_E_H                  &
                                    ,HANDLE_PACKET_E_V                  &
                                    ,HANDLE_PARENT_DOM_LIMITS           &
                                    ,HANDLE_PARENT_ITE                  &
                                    ,HANDLE_PARENT_ITS                  &
                                    ,HANDLE_PARENT_JTE                  &
                                    ,HANDLE_PARENT_JTS                  &
                                    ,INFO_SEND                          &
                                    ,PTASK_LIMITS
!
      USE module_NESTING,ONLY: PARENT_CHILD_COMMS
!
      USE module_PARENT_CHILD_CPL_COMP,ONLY: PARENT_CHILD_CPL_REGISTER  &  !<-- The Register routine for PARENT_CHILD Coupler
                                            ,PARENT_CHILD_COUPLER_SETUP
!
      USE module_CONTROL,ONLY: NUM_DOMAINS_MAX,TIMEF
!
      USE module_CLOCKTIMES,ONLY: TIMERS                                &
                                 ,cbcst_tim,pbcst_tim
!
      USE module_ERROR_MSG,ONLY: ERR_MSG,MESSAGE_CHECK
!
      USE MODULE_SOLVER_GRID_COMP,ONLY: RESTVAL
!
      USE module_CONSTANTS,ONLY: A
!
      USE module_CPLFIELDS,ONLY: exportFieldsList,importFieldsList      &
                                ,queryFieldList
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      PRIVATE
!
      PUBLIC :: NMM_REGISTER
      PUBLIC :: EXPORT_FIELDS_INDX                                      &
               ,NUM_DOMAINS_TOTAL,nExportFields_NMMB,NMM_GRID
!
!-----------------------------------------------------------------------
!
      INTEGER(kind=KINT) :: MYPE                                        &  !<-- Each MPI task ID
                           ,NHOURS_CLOCKTIME                               !<-- Fcst hours between prints of integration clocktime
!
      INTEGER(kind=KINT) :: FILT_TIMESTEP_SEC_WHOLE                     &
                           ,FILT_TIMESTEP_SEC_NUMERATOR                 &
                           ,FILT_TIMESTEP_SEC_DENOMINATOR               &
                           ,TIMESTEP_SEC_WHOLE                          &
                           ,TIMESTEP_SEC_NUMERATOR                      &
                           ,TIMESTEP_SEC_DENOMINATOR                    &
                           ,TIMESTEPS_RESTART
!
      INTEGER(kind=KINT),SAVE :: COMM_GLOBAL                            &  !<-- The MPI communicator for all tasks (COMM_WORLD)
                                ,FILTER_METHOD                             !<-- Digital filter flag (0->no filter; >0->filter type)
!
      INTEGER(kind=KINT),POINTER :: COMM_TO_MY_PARENT                   &  !<-- Intercommunicator between a domain and its parent
                                   ,NPE_PRINT                              !<-- Clocktime diagnostics from this MPI task
!
      REAL(kind=KFPT),SAVE :: TLM0D,TPH0D                                  !<-- Central geographic lon/lat (degrees) of primary domain.
!
      REAL(kind=KFPT),DIMENSION(:),ALLOCATABLE,SAVE :: DT               &  !<-- Each domain's fundamental timestep
                                                      ,FILT_DT             !<-- Filter timestep (s) of the domains
!
      CHARACTER(ESMF_MAXSTR) :: ALARM_CPL_NAME                          &  !<-- Name of the ESMF Alarms for coupling intervals
                               ,CLOCK_NMM_NAME                             !<-- Name of the NMM's ESMF Clock
!
      LOGICAL(kind=KLOG) :: PRINT_TIMING                                   !<-- Print timing flag
!
      LOGICAL(kind=KLOG),POINTER :: RESTARTED_RUN                       &  !<-- Flag indicating if this is a restarted run
                                   ,RST_OUT_00                             !<-- Shall we write 00h history in restarted run?
!
      TYPE(ESMF_VM),SAVE :: VM                                             !<-- The ESMF virtual machine.
!
      TYPE(ESMF_Time),SAVE :: STARTTIME                                    !<-- The ESMF start time.
!
      TYPE(ESMF_TimeInterval),SAVE :: INTERVAL_CPL                      &  !<-- The NUOPC coupling interval (sec ESMF)
                                     ,RUNDURATION                          !<-- The simulation length (sec ESMF)
!
      TYPE(ESMF_TimeInterval),POINTER :: FILT_TIMESTEP                  &  !<-- ESMF filter timestep (s)
                                        ,INTERVAL_CLOCKTIME             &  !<-- ESMF time interval between clocktime prints (h)
                                        ,INTERVAL_HISTORY               &  !<-- ESMF time interval between history output (h)
                                        ,INTERVAL_RESTART               &  !<-- ESMF time interval between restart output (h)
                                        ,TIMESTEP                          !<-- The ESMF timestep (s)
!
      TYPE(ESMF_Alarm),DIMENSION(:),ALLOCATABLE,SAVE :: ALARM_CPL          !<-- The ESMF Alarms for the coupling interval
!
      TYPE(ESMF_Clock),DIMENSION(:),ALLOCATABLE,SAVE :: CLOCK_NMM          !<-- The NMM ESMF Clocks
!
      TYPE(ESMF_GridComp),POINTER :: DOMAIN_GRID_COMP                      !<-- A domain's ESMF component
!
      TYPE(ESMF_State),POINTER :: EXP_STATE_DOMAIN                      &  !<-- A domain's export state
                                 ,IMP_STATE_DOMAIN                         !<-- A domain's import state
!
      TYPE(NMM_INTERNAL_STATE),POINTER,SAVE :: NMM_INT_STATE               !<-- The NMM component internal state pointer
!
      INTEGER(kind=KINT),DIMENSION(:),POINTER,SAVE ::                   &
                                         RANK_TO_DOMAIN_ID_PTR             !<-- Save converter for runtime timestep retrieval
!
      TYPE(WRAP_NMM_INTERNAL_STATE),SAVE :: WRAP                           !<-- The F90 wrap of the NMM internal state
!
!
!---------------------
!***  For NMM Nesting
!---------------------
!
      INTEGER(kind=KINT),PARAMETER :: nExportFields_NMMB=7                 !<-- The # of fields to regrid from the full export field list
!
      INTEGER(kind=KINT),SAVE :: NUM_DOMAINS_TOTAL
!
      INTEGER(kind=KINT) :: KOUNT_STEPS=0                               &
                           ,NUM_GENS=1                                     !<-- The # of generations of domains (only for 2-way nests)
!
      INTEGER(kind=KINT) :: COMM_MY_DOMAIN                              &  !<-- Each domain's local intracommunicator
                           ,FULL_GEN                                    &  !<-- The 1st generation of domains that uses all fcst tasks
                           ,MY_DOMAIN_ID                                &  !<-- The ID of each domain
                           ,NPHS                                        &  !<-- The physics timestep
                           ,NTRACK                                      &  !<-- The storm locator flag
                           ,NUM_DOMAINS_MINE                               !<-- The # of domains on which each task resides
!
      INTEGER(kind=KINT),POINTER :: NUM_CHILDREN                        &  !<-- # of children on a domain
                                   ,NUM_2WAY_CHILDREN                   &  !<-- # of 2-way children on a domain
                                   ,PARENT_CHILD_TIME_RATIO                !<-- Ratio of parent timestep to child's
!
      INTEGER(kind=KINT),DIMENSION(nExportFields_NMMB),SAVE ::          &
                                                        EXPORT_FIELDS_INDX !<-- Index of the desired export fields from the full list
!
      INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE :: DOMAIN_GEN         &  !<-- The generation of each domain
                                                    ,MY_DOMAIN_IDS      &  !<-- All domains on which each task resides
                                                    ,MY_DOMAINS_IN_GENS    !<-- List a task's domain on each generation
!
      INTEGER(kind=KINT),DIMENSION(:),POINTER :: ID_DOMAINS             &  !<-- IDs of all domains
                                                ,ID_PARENTS             &  !<-- IDs of all domains' parents
                                                ,FTASKS_DOMAIN          &  !<-- # of forecast tasks on each domain excluding descendents
                                                ,NTASKS_DOMAIN             !<-- # of tasks on each domain excluding descendents
!
      INTEGER(kind=KINT),DIMENSION(:,:),POINTER :: ID_CHILDREN          &  !<-- IDs of all children of all domains
                                                  ,PETLIST_DOMAIN          !<-- List of task IDs for each domain (DOMAIN Component)
!
      CHARACTER(len=12),SAVE :: TASK_MODE                                  !<-- Task assignments are unique or generational 
!
      CHARACTER(len=5) :: NEST_MODE                                        !<-- Is the nesting 1-way or 2-way with the parent?
!
      CHARACTER(len=40),DIMENSION(nExportFields_NMMB) :: EXPORT_FIELDS_BLEND &
                                 =(/                                         &
                                    'mean_zonal_moment_flx'                  &
                                   ,'mean_merid_moment_flx'                  &
                                   ,'inst_sensi_heat_flx'                    &
                                   ,'inst_laten_heat_flx'                    &
                                   ,'inst_net_lw_flx'                        &
                                   ,'inst_net_sw_flx'                        &
                                   ,'inst_pres_height_surface'               &
                                                              /)
!
      LOGICAL(kind=KLOG),SAVE :: ALL_FORECASTS_COMPLETE=.FALSE.         &  !<-- Are this task's domains' fcsts all finished?
                                ,NESTING_NMM                               !<-- Does this run contain nests?
!
      LOGICAL(kind=KLOG) :: MY_DOMAIN_MOVES                                !<-- Does my domain move?
!
      LOGICAL(kind=KLOG),POINTER :: I_AM_A_FCST_TASK                    &  !<-- Am I a forecast task?
                                   ,I_AM_LEAD_FCST_TASK                 &  !<-- Am I the lead forecast task?
                                   ,I_AM_A_NEST                            !<-- Am I in a nested domain?
!
      LOGICAL(kind=KLOG),DIMENSION(:),ALLOCATABLE,SAVE :: FREE_TO_INTEGRATE   & !<-- A yes/no flag for 2-way domains calling DOMAIN_RUN
                                                         ,GENERATION_FINISHED   !<-- Flag of when forecast is done per generation
!
      TYPE(COMMS_FAMILY),DIMENSION(:),POINTER :: COMMS_DOMAIN              !<-- Intracommunicators between parents and children
                                                                           !    and between each domains' forecast tasks
!
      TYPE(ESMF_Config),DIMENSION(NUM_DOMAINS_MAX),SAVE :: CF              !<-- The config objects (one per domain)
!
      TYPE(ESMF_State),POINTER :: IMP_STATE_CPL_NEST                    &
                                 ,EXP_STATE_CPL_NEST
!
      TYPE(ESMF_CplComp),POINTER :: PARENT_CHILD_COUPLER_COMP              !<-- Coupler component for parent-child/nest exchange
!
!-----------------------------------------------------------------------
!
!---------------------------
!***  For Digital Filtering
!---------------------------
!
      INTEGER(kind=KINT),PUBLIC :: DFIHR,DFIHR_CHK
!
      TYPE(ESMF_Time),SAVE,PUBLIC :: DFITIME
!
!-----------------------------------------------------------------------
!
!-----------
!*** Timing
!-----------
!
      REAL(kind=KDBL) :: btim,btim0 
!
!-----------------------------------------------------------------------
!
      character(len=160) :: nuopcMsg
      type fld_list_type
        character(len=64) :: stdname
        character(len=64) :: shortname
        character(len=64) :: transferOffer
        logical           :: assoc    ! is the farrayPtr associated with internal data
        real(ESMF_KIND_R8), dimension(:,:,:), pointer :: farrayPtr
      end type fld_list_type

      integer,parameter :: fldsMax = 100
      integer :: fldsToNMMB_num = 0
      type (fld_list_type) :: fldsToNMMB(fldsMax)
      integer :: fldsFrNMMB_num = 0
      type (fld_list_type) :: fldsFrNMMB(fldsMax)
      character(len=2048):: info
      type(ESMF_Grid), save    :: nmm_grid

!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE NMM_REGISTER(NMM_GRID_COMP,RC_REG)
!
!-----------------------------------------------------------------------
!***  Register the NMM component's Initialize, Run, and Finalize steps.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(ESMF_GridComp) :: NMM_GRID_COMP                                 !<-- The NMM component
!
      INTEGER,INTENT(OUT) :: RC_REG
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER :: RC
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC    =ESMF_SUCCESS   ! Error signal variable
      RC_REG=ESMF_SUCCESS   ! Error signal variable
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Register the NMM Initialize routine"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompSetEntryPoint(NMM_GRID_COMP                     &
                                     ,ESMF_METHOD_INITIALIZE            &
                                     ,NMM_INITIALIZE                    &
                                     ,phase=1                           &
                                     ,rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Register the NMM Run routine"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
! 
      CALL ESMF_GridCompSetEntryPoint(NMM_GRID_COMP                     &
                                     ,ESMF_METHOD_RUN                   &
                                     ,NMM_RUN                           &
                                     ,phase=1                           &
                                     ,rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Register the NMM Finalize routine"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompSetEntryPoint(NMM_GRID_COMP                     &
                                     ,ESMF_METHOD_FINALIZE              &
                                     ,NMM_FINALIZE                      &
                                     ,phase=1                           &
                                     ,rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
      IF(RC_REG==ESMF_SUCCESS)THEN
!       WRITE(0,*)' NMM_REGISTER succeeded'
      ELSE
        WRITE(0,*)' NMM_REGISTER failed  RC_REG=',RC_REG
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE NMM_REGISTER
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE NMM_INITIALIZE(NMM_GRID_COMP                           &
                               ,IMP_STATE                               &
                               ,EXP_STATE                               &
                               ,CLOCK_NEMS                              &
                               ,RC_INIT)
!
!-----------------------------------------------------------------------
!***  This routine creates the individual DOMAIN gridded components 
!***  and executes their Initialize step.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(ESMF_GridComp) :: NMM_GRID_COMP                                 !<-- The NMM component
!
      TYPE(ESMF_State) :: IMP_STATE                                     &  !<-- The NMM import state
                         ,EXP_STATE                                        !<-- The NMM export state
!
      TYPE(ESMF_Clock) :: CLOCK_NEMS                                       !<-- The NEMS ESMF Clock
!
      INTEGER,INTENT(OUT) :: RC_INIT                                       !<-- Error return code
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: ID,ID_DOM,ID_X,IERR,INDX2,ISTAT,KOUNT       &
                           ,N,N1,N2,N3,NN,NT
!
      INTEGER(kind=KINT) :: MINUTES_HISTORY                             &  !<-- Hours between history output
                           ,MINUTES_RESTART                             &  !<-- Hours between restart output
                           ,NHOURS_FCST                                 &  !<-- Length of forecast in hours
                           ,NSECONDS_FCST                               &  !<-- Length of forecast in seconds
                           ,TIMESTEP_FINAL                                 !<-- # of timesteps in entire forecast
!
      INTEGER(kind=KINT) :: ISECOND_FCST,ISECOND_NUM,ISECOND_DEN
!
      INTEGER(kind=KINT) :: GEN_X,INPES,JNPES,LEAD_TASK                 &
                           ,LENGTH,LENGTH_FCST,LENGTH_FCST_1            &
                           ,MAX_GEN,MY_DOMAIN_ID,MYPE_LOCAL,MYPE_X      &
                           ,N_GEN,N_TASKS,NUM_CHILD_TASKS               &
                           ,NUM_DOMAINS_X,NUM_FCST_TASKS                &
                           ,NUM_WRITE_TASKS,NUM_TASKS_TOTAL             &
                           ,WRITE_GROUPS,WRITE_TASKS_PER_GROUP
!
      INTEGER(kind=KINT) :: RC
!
      INTEGER(kind=KINT),DIMENSION(:),POINTER :: MY_DOMAIN_ID_N            !<-- Domain IDs for each task among all domains
!
      INTEGER(kind=KINT),DIMENSION(NUM_DOMAINS_MAX) :: DOMAIN_ID_TO_RANK=0  &  !<-- The configure file associated with each domain ID
                                                      ,RANK_TO_DOMAIN_ID=0     !<-- The domain ID associated with each configure file
!
      INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE :: N_FCST_TASKS_GEN      !<-- The # of fcst tasks in each geenration
!
      INTEGER(kind=KINT),DIMENSION(:),POINTER :: CHILD_ID               &
                                                ,COMM_TO_MY_CHILDREN    &  !<-- Intercommunicators between a domain and its children
                                                ,PETLIST
!
      INTEGER,DIMENSION(MPI_STATUS_SIZE) :: JSTAT
!
      REAL(kind=KFPT) :: TIMESTEP_RATIO
!
      LOGICAL(kind=KLOG) :: CFILE_EXIST                                 &
                           ,QUILTING,QUILTING_1                         &
                           ,USED_ALL_FCST_TASKS
!
      CHARACTER(2) :: INT_TO_CHAR
      CHARACTER(6) :: FMT='(I2.2)'
      CHARACTER(7) :: MODE
      CHARACTER(NUM_DOMAINS_MAX) :: CONFIG_FILE_NAME
!
      CHARACTER(ESMF_MAXSTR) :: DOMAIN_COMP_BASE='DOMAIN Gridded Component ' &
                               ,DOMAIN_GRID_COMP_NAME,STATE_NAME
!
      TYPE(ESMF_TimeInterval) :: COUPLING_INTERVAL                         !<-- ESMF time interval (sec) between coupling times
      TYPE(ESMF_TimeInterval) :: TIMEINTERVAL_RECV_FROM_PARENT             !<-- ESMF time interval between Recv times from parent
      TYPE(ESMF_TimeInterval) :: ZERO_INTERVAL                             !<-- Zero time interval used in comparison of time step
!                                                                          !    and restart interval.
      LOGICAL :: PHYSICS_ON                                                !<-- Does the integration include physics?

      TYPE(ESMF_Config) :: CF_X                                            !<-- Working config object
!
      TYPE(ESMF_Grid)    :: pGrid_NMMB
!
      TYPE(DOMAIN_INTERNAL_STATE),POINTER :: DOMAIN_INT_STATE
      TYPE(WRAP_DOMAIN_INTERNAL_STATE) :: WRAP_DOMAIN
!
      type(esmf_time) :: ringtime,master_starttime,master_stoptime
      type(esmf_timeinterval) :: ringinterval
      integer(kind=kint) :: iyear_fcst &
                           ,imonth_fcst &
                           ,iday_fcst &
                           ,ihour_fcst &
                           ,iminute_fcst  
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC     =ESMF_SUCCESS
      RC_INIT=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  Allocate the NMM component's internal state.
!-----------------------------------------------------------------------

      ALLOCATE(NMM_INT_STATE,stat=RC)
      wrap%NMM_INT_STATE=>NMM_INT_STATE
!
!-----------------------------------------------------------------------
!***  Attach the NMM internal state to the NMM component.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Attach NMM Internal State to the NMM Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompSetInternalState(NMM_GRID_COMP                  &
                                        ,WRAP                           &
                                        ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Retrieve the VM (Virtual Machine) of the NMM component.
!***  We need VM now to obtain the MPI task IDs.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Retrieve VM from NMM Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompGet(gridcomp=NMM_GRID_COMP                      &  !<-- The NMM component
                           ,vm      =VM                                 &  !<-- Get the Virtual Machine from the NMM component
                           ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Obtain MPI Task IDs from VM"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_VMGet(vm      =VM                                       &  !<-- The virtual machine
                     ,localpet=MYPE                                     &  !<-- Each MPI global task ID (all tasks are present)
                     ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Create and load all of the configure objects.  All domains
!***  are functionally equivalent thus each has its own configure
!***  file.  We are counting the configure files as we create the
!***  ESMF configure objects so we will know how many different
!***  domains there are.
!-----------------------------------------------------------------------
!
      NUM_DOMAINS_X=0
!
      ALLOCATE(RANK_TO_DOMAIN_ID_PTR(NUM_DOMAINS_MAX))
!
      DO N=1,NUM_DOMAINS_MAX                                               !<-- Number of config files must not exceed 99
!
        WRITE(INT_TO_CHAR,FMT)N
        CONFIG_FILE_NAME='configure_file_'//INT_TO_CHAR                    !<-- Each configure file has a unique number.
!
        CFILE_EXIST=.FALSE.
        INQUIRE(FILE=CONFIG_FILE_NAME,EXIST=CFILE_EXIST)
!
        IF(CFILE_EXIST)THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="NMM_INIT: Create Temporary Configure Object"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CF_X=ESMF_ConfigCreate(rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="NMM_INIT: Load the Temp Configure Object"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_ConfigLoadFile(config  =CF_X                        &
                                  ,filename=CONFIG_FILE_NAME            &
                                  ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="NMM_INIT: Extract Domain ID From Temp Config File"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_ConfigGetAttribute(config=CF_X                      &  !<-- The config object
                                      ,value =ID_X                      &  !<-- The domain's ID
                                      ,label ='my_domain_id:'           &  !<-- Take value from this config labelious variable
                                      ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CF(ID_X)=ESMF_ConfigCreate(rc=RC)                                !<-- Domain's ID is its element in the CF array
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="NMM_INIT: Destroy Temporary Config Object"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_ConfigDestroy(config=CF_X                           &  !<-- The temporary config object
                                 ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="NMM_INIT: Load the Nest Configure Object"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_ConfigLoadFile(config  =CF(ID_X)                    &
                                  ,filename=CONFIG_FILE_NAME            &
                                  ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          DOMAIN_ID_TO_RANK(ID_X)=N                                        !<-- The configure file rank for a given domain ID
          RANK_TO_DOMAIN_ID(N)=ID_X                                        !<-- The domain ID for a given configure file rank
          RANK_TO_DOMAIN_ID_PTR(N)=ID_X                                    !<-- The domain ID for a given configure file rank
!
          NUM_DOMAINS_X=NUM_DOMAINS_X+1
!
        ELSE
!
          EXIT
!
        ENDIF
!
      ENDDO
!
      NESTING_NMM=.FALSE.
      IF(NUM_DOMAINS_X>1)NESTING_NMM=.TRUE.                                !<-- We have nests if more than one domain is present
!
!-----------------------------------------------------------------------
!***  Before going further we need to be certain that the number of
!***  configure files present actually matches the number of domains
!***  the user intends there to be.  If they do not match then abort
!***  the run.  The uppermost domain's configure file contains the
!***  total number of domains that the user wants.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!     MESSAGE_CHECK="NMM_INIT: Extract Total Domain Count Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF(1)                         &  !<-- The config object
                                  ,value =NUM_DOMAINS_TOTAL             &  !<-- The user-specified total number of domains
                                  ,label ='num_domains_total:'          &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!     CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF (NUM_DOMAINS_X/=NUM_DOMAINS_TOTAL) THEN
        WRITE(0,*)' # of configure files in working directory is wrong!'
        WRITE(0,*)' You have said there are ',NUM_DOMAINS_TOTAL         &
                 ,' domains in this run.'
        WRITE(0,*)' But there are ',NUM_DOMAINS_X,' configure files present.'
        WRITE(0,*)' There must be one configure file per domain.'
        WRITE(0,*)' ABORTING!!'
        CALL ESMF_Finalize(rc=RC,endflag=ESMF_END_ABORT)
      ENDIF
!
!-----------------------------------------------------------------------
!***  Obtain the global communicator for all tasks in this run.
!-----------------------------------------------------------------------
!
      CALL ESMF_VMGet(vm             =VM                                &  !<-- The virtual machine
                     ,mpiCommunicator=COMM_GLOBAL                       &  !<-- Global intracommunicator for all tasks
                     ,petCount       =NUM_TASKS_TOTAL                   &  !<-- Total # of tasks in this run
                     ,rc             =RC)
!
!-----------------------------------------------------------------------
!***  Check now if the user has specified quilt tasks.  If there
!***  are multiple domains then either all or none must set quilting
!***  to false.
!-----------------------------------------------------------------------
!
      DO N=1,NUM_DOMAINS_TOTAL
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="NMM_INIT: Get Value of QUILTING from Config Files"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=CF(N)                   &  !<-- The config object of domain N
                                    ,value =QUILTING                &  !<-- Has quilting been specified?
                                    ,label ='quilting:'             &  !<-- Give this label's value to the previous variable
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        IF(N==1)THEN
          QUILTING_1=QUILTING
        ELSE
          IF(QUILTING.AND..NOT.QUILTING_1)THEN
            WRITE(0,*)' Conflicting quilting settings in configure files!'
            WRITE(0,*)' Aborting!!!'
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
          ENDIF
        ENDIF
!
      ENDDO
!
!-----------------------------------------------------------------------
!***  Set the default for the mode of MPI task assignment.
!-----------------------------------------------------------------------
!
      TASK_MODE='unique'
!
!-----------------------------------------------------------------------
!***  IF NESTED DOMAINS ARE BEING USED THEN:
!***    (1) Split the MPI Communicator between all domains;
!***    (2) Create a DOMAIN subcomponent for all domains;
!***    (3) Call DOMAIN_INIT recursively for all domains.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
      nesting_block_1: IF(NESTING_NMM)THEN                                 !<-- Special communicators are needed for nesting
!
!-----------------------------------------------------------------------
!***  There is no need to proceed if the specified forecast lengths
!***  of all domains are not the same.  Currently the upper parent
!***  cannot integrate longer than its children and some nests cannot
!***  integrate longer than other nests.
!-----------------------------------------------------------------------
!
        DO N=1,NUM_DOMAINS_TOTAL
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="NMM_INIT: Check forecast lengths of domains."
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_ConfigGetAttribute(config=CF(N)                     &  !<-- The config object of domain N
                                      ,value =LENGTH_FCST               &  !<-- Forecast length of domain N
                                      ,label ='nhours_fcst:'            &  !<-- Configure label for forecast length
                                      ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          IF(N==1)THEN
            LENGTH_FCST_1=LENGTH_FCST
          ELSE
            IF(LENGTH_FCST/=LENGTH_FCST_1)THEN
              WRITE(0,*)' Domain forecast lengths differ!'
              WRITE(0,*)' Aborting!'
              CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
            ENDIF
          ENDIF
!
        ENDDO
!
!-----------------------------------------------------------------------
!***  Task 0 checks all the configure files to see if 2-way exchange
!***  appears in any of them.  If it does then the mode for this run's
!***  task assignments is generational and not unique to each domain.
!-----------------------------------------------------------------------
!
        IF(MYPE==0)THEN
!
          search: DO N=1,NUM_DOMAINS_TOTAL
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="NMM_INIT: Is Nest_Mode in the Configure File?"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_ConfigFindLabel(config=CF(N)                      &  !<-- The config object of domain N
                                     ,label ='nest_mode:'               &  !<-- Domain N's nesting mode ('1-way' or '2-way')
                                     ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            IF(RC==-2)THEN                                                 !<-- Return code is -2 if the label is not found
!
              CYCLE                                                        !<-- nest_mode not in config file (domain not a child)
!
            ELSE
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
              MESSAGE_CHECK="NMM_INIT: Check Exchange Mode in Config Files"
!             CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
              CALL ESMF_ConfigGetAttribute(config=CF(N)                   &  !<-- The config object of domain N
                                          ,value =NEST_MODE               &  !<-- Domain N's nesting mode ('1-way' or '2-way')
                                          ,label ='nest_mode:'            &  !<-- Give this label's value to the previous variable
                                          ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
              CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
              IF(NEST_MODE=='2-way')THEN
                TASK_MODE='generational'
                EXIT search
              ENDIF
!
            ENDIF
!
          ENDDO search
!
        ENDIF
!
        CALL MPI_BCAST(TASK_MODE                                        &  !<-- Broadcast the value of TASK_MODE
                      ,12                                               &  !<-- It contains 12 characters
                      ,MPI_CHARACTER                                    &  !<-- Type CHARACTER
                      ,0                                                &  !<-- Global task 0 is sending
                      ,COMM_GLOBAL                                      &  !<-- The global communicator
                      ,IERR )
!
!-----------------------------------------------------------------------
!
        ALLOCATE(DOMAIN_GEN(1:NUM_DOMAINS_TOTAL),stat=ISTAT)               !<-- For 2-way nesting, the generation of each domain
        IF(ISTAT/=0)THEN
          WRITE(0,*)' Failed to allocate DOMAIN_GEN  istat=',ISTAT
        ENDIF
!
        DO N=1,NUM_DOMAINS_TOTAL
          DOMAIN_GEN(N)=0                                                  !<-- Initialize value of each domain's generation
        ENDDO
!
        FULL_GEN=0           
!
!-----------------------------------------------------------------------
!***  If the task assignment is generational then we must check to be
!***  sure that at least one of the generations of nests uses all of
!***  the tasks assigned to the run.
!-----------------------------------------------------------------------
!
        two_way: IF(TASK_MODE=='generational')THEN
!
!-----------------------------------------------------------------------
!***  Read all the configure files to find out how many generations
!***  of nests there are by checking which generation all the domains
!***  are in.
!-----------------------------------------------------------------------
!
          MAX_GEN=0
          NUM_WRITE_TASKS=0
!
          DO N=1,NUM_DOMAINS_TOTAL
!
            ID_X=RANK_TO_DOMAIN_ID(N)                                      !<-- The domain ID for the Nth domain
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="NMM_INIT: Extract generation from Config File"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_ConfigGetAttribute(config=CF(ID_X)                &  !<-- The config object
                                        ,value =DOMAIN_GEN(ID_X)        &  !<-- This domain's generation
                                        ,label ='generation:'           &  !<-- Give this label's value to the previous variable
                                        ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            MAX_GEN=MAX(MAX_GEN,DOMAIN_GEN(ID_X))
!
          ENDDO
!
          ALLOCATE(N_FCST_TASKS_GEN(1:MAX_GEN),stat=ISTAT)
          IF(ISTAT/=0)THEN
            WRITE(0,*)' Failed to allocate N_FCST_TASKS_GEN  istat=',ISTAT
          ENDIF
!
          DO N=1,MAX_GEN
            N_FCST_TASKS_GEN(N)=0                                          !<-- Initialize # of tasks in each generation
          ENDDO
!
!-----------------------------------------------------------------------
!***  Now determine and check all generations' task counts.
!-----------------------------------------------------------------------
!
          DO N=1,NUM_DOMAINS_TOTAL
!
            ID_X=RANK_TO_DOMAIN_ID(N)                                      !<-- The domain ID for the Nth domain
            N_GEN=DOMAIN_GEN(ID_X)                                         !<-- This domain's generation
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="NMM_INIT: Extract INPES From Config File"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_ConfigGetAttribute(config=CF(ID_X)                &  !<-- The config object
                                        ,value =INPES                   &  !<-- The domain's fcst tasks in I
                                        ,label ='inpes:'                &  !<-- Give this label's value to the previous variable
                                        ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="NMM_INIT: Extract JNPES From Config File"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_ConfigGetAttribute(config=CF(ID_X)                &  !<-- The config object
                                        ,value =JNPES                   &  !<-- The domain's fcst tasks in J
                                        ,label ='jnpes:'                &  !<-- Give this label's value to the previous variable
                                        ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="NMM_INIT: Extract Write_Groups From Config File"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_ConfigGetAttribute(config=CF(ID_X)                &  !<-- The config object
                                        ,value =WRITE_GROUPS            &  !<-- The number of Write groups on this domain
                                        ,label ='write_groups:'         &  !<-- Give this label's value to the previous variable
                                        ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="NMM_INIT: Extract Write Tasks Per Group From Config File"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_ConfigGetAttribute(config=CF(ID_X)                 &  !<-- The config object
                                        ,value =WRITE_TASKS_PER_GROUP    &  !<-- The number of tasks per Write group
                                        ,label ='write_tasks_per_group:' &  !<-- Give this label's value to the previous variable
                                        ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            N_FCST_TASKS_GEN(N_GEN)=N_FCST_TASKS_GEN(N_GEN)+INPES*JNPES
!
            NUM_WRITE_TASKS=NUM_WRITE_TASKS                             &
                           +WRITE_GROUPS*WRITE_TASKS_PER_GROUP
!
          ENDDO
!
!-----------------------------------------------------------------------
!
          USED_ALL_FCST_TASKS=.FALSE.
          NUM_FCST_TASKS=NUM_TASKS_TOTAL-NUM_WRITE_TASKS                   !<-- # of forecast tasks available
!
          DO N=1,MAX_GEN
!
            IF(N_FCST_TASKS_GEN(N)==NUM_FCST_TASKS.AND.FULL_GEN==0)THEN
              USED_ALL_FCST_TASKS=.TRUE.
              FULL_GEN=N                                                   !<-- Save the 1st generation that uses all fcst tasks
!
            ELSEIF(N_FCST_TASKS_GEN(N)>NUM_FCST_TASKS)THEN
              WRITE(0,*)' Generation ',N,' is using more fcst tasks'    &
                       ,' than assigned to the run!'
              WRITE(0,*)' There are ',NUM_FCST_TASKS,' fcst tasks in this run.'
              WRITE(0,*)' Generation ',N,' is using ',N_FCST_TASKS_GEN(N) &
                       ,' fcst tasks.'
              WRITE(0,*)' Aborting!!'
              CALL ESMF_Finalize(rc=RC,endflag=ESMF_END_ABORT)
!
            ELSEIF(N_FCST_TASKS_GEN(N)==0)THEN
              WRITE(0,*)' Generation ',N,' is using no tasks!!'
              WRITE(0,*)' Aborting!!'
              CALL ESMF_Finalize(rc=RC,endflag=ESMF_END_ABORT)
!
            ENDIF
!
          ENDDO
!
          IF(.NOT.USED_ALL_FCST_TASKS)THEN
            WRITE(0,*)' No generation is using all fcst tasks assigned' &
                     ,' to the run.'
            WRITE(0,*)' At least one generation must use all fcst tasks.'
            WRITE(0,*)' Aborting!!'
            CALL ESMF_Finalize(rc=RC,endflag=ESMF_END_ABORT)
          ENDIF
!
          DEALLOCATE(N_FCST_TASKS_GEN)
!
!-----------------------------------------------------------------------
!
        ENDIF  two_way
!
!-----------------------------------------------------------------------
!***  Split the global communicator among all NMM domains and create
!***  Parent-Child intercommunicators.
!-----------------------------------------------------------------------
!
        CALL PARENT_CHILD_COMMS(MYPE                                    &  !<-- This task's global rank (in)
                               ,NUM_DOMAINS_TOTAL                       &  !<-- Total number of domains, all generations (in)
                               ,NUM_TASKS_TOTAL                         &  !<-- Total number of tasks assigned to the run (in)
                               ,COMM_GLOBAL                             &  !<-- Intracommunicator for ALL tasks (in)
                               ,RANK_TO_DOMAIN_ID                       &  !<-- Domain IDs for each configure file
                               ,CF                                      &  !<-- Configure objects for all domains (in)
                               ,TASK_MODE                               &  !<-- 1-way or 2-way nesting (in)
                               ,QUILTING                                &  !<-- Was quilting specified in the configure files?
                               ,DOMAIN_GEN                              &  !<-- For 2-way nesting, the generation of each domain (in)
                               ,FULL_GEN                                &  !<-- For 2-way nesting, the generation using all tasks (in)
                               ,MY_DOMAIN_ID_N                          &  !<-- ID of domains on which this task resides (out)
                               ,ID_DOMAINS                              &  !<-- IDs of all domains (out)
                               ,ID_PARENTS                              &  !<-- ID of all domains' parents (out)
                               ,nmm_int_state%NUM_CHILDREN              &  !<-- # of children on each domain (out)
                               ,ID_CHILDREN                             &  !<-- IDs of all children of all domains (out)
                               ,COMMS_DOMAIN                            &  !<-- Parent and child intracommunicators (out)
                               ,FTASKS_DOMAIN                           &  !<-- # of fcst tasks on each domain excluding descendents (out)
                               ,NTASKS_DOMAIN                           &  !<-- # of tasks on each domain excluding descendents (out)
                               ,PETLIST_DOMAIN                          &  !<-- List of task IDs for each domain (DOMAIN Component) (out)
                               ,NUM_GENS )                                 !<-- # of generations of domains (out)
!
!-----------------------------------------------------------------------
!***  The array MY_DOMAIN_ID_N is dimensioned 1:NUM_DOMAINS_TOTAL.
!***  The indices that correspond to the domain IDs on which the
!***  current task lies equal those respective IDs.  All the other
!***  indices are zero.  Now how many domains does the current
!***  task lie on?
!-----------------------------------------------------------------------
!
        IF(TASK_MODE=='unique')THEN
!
!-----------------------------------------------------------------------
!
          NUM_DOMAINS_MINE=1
!
          ALLOCATE(MY_DOMAIN_IDS(1:1))
          MY_DOMAIN_IDS(1)=MY_DOMAIN_ID_N(1)                               !<-- A task lies on only one domain in 1-way nests
!
          ALLOCATE(MY_DOMAINS_IN_GENS(1:1))
          MY_DOMAINS_IN_GENS(1)=MY_DOMAIN_ID_N(1)                          !<-- The task only lies on one domain
!
          ALLOCATE(GENERATION_FINISHED(1:1))
          GENERATION_FINISHED(1)=.FALSE.                                   !<-- Generations not relevant for 1-way nesting
!
          ALLOCATE(FREE_TO_INTEGRATE(1:1))                                 !<-- 1-way => task on one domain per generation; always free
          FREE_TO_INTEGRATE(1)=.TRUE.
!
!-----------------------------------------------------------------------
!
        ELSEIF(TASK_MODE=='generational')THEN
!
!-----------------------------------------------------------------------
!
          NUM_DOMAINS_MINE=0
!
          DO N=1,NUM_DOMAINS_TOTAL
            IF(MY_DOMAIN_ID_N(N)>0)THEN
              NUM_DOMAINS_MINE=NUM_DOMAINS_MINE+1
            ENDIF
          ENDDO
!
!-----------------------------------------------------------------------
!***  What are the domain IDs of the domains this task lies on?
!***  Save those in an array dimensioned to exactly the number
!***  of domains it actually lies on.
!-----------------------------------------------------------------------
!
          ALLOCATE(MY_DOMAIN_IDS(1:NUM_DOMAINS_MINE),stat=ISTAT)
          IF(ISTAT/=0)THEN
            WRITE(0,*)' Failed to allocate MY_DOMAIN_IDS  istat=',ISTAT
          ENDIF
!
          KOUNT=0
          DO N=1,NUM_DOMAINS_TOTAL
            IF(MY_DOMAIN_ID_N(N)/=0)THEN
              KOUNT=KOUNT+1
              MY_DOMAIN_IDS(KOUNT)=MY_DOMAIN_ID_N(N)
            ENDIF
          ENDDO
!
!-----------------------------------------------------------------------
!***  Any task can be in no more than one domain per generation.
!***  Save the domain IDs that the current task lies on within
!***  each generation.  If the task does not lie on any domain
!***  in a generation then the value is 0.
!
!***  MY_DOMAINS_IN_GENS - A task's domain in each generation (1: # of generations) 
!***  MY_DOMAIN_IDS - All domains on which a task resides (1: # of domains a task is on)
!***  DOMAIN_GEN - The generation of each domain (1: total # of domains)
!-----------------------------------------------------------------------
!
          ALLOCATE(MY_DOMAINS_IN_GENS(1:NUM_GENS))
!
          DO N=1,NUM_GENS
            MY_DOMAINS_IN_GENS(N)=0
          ENDDO
!
          DO N=1,NUM_DOMAINS_MINE                                          !<-- Loop through the domains this task is on
            ID_DOM=MY_DOMAIN_IDS(N)                                        !<-- This task's Nth domain
            GEN_X=DOMAIN_GEN(ID_DOM)                                       !<-- The generation of that domain
            IF(MY_DOMAINS_IN_GENS(GEN_X)>0)THEN                            !<-- The task already has a domain in that generation?
              WRITE(0,*)' ERROR'
              WRITE(0,*)' Domain ID is ',ID_DOM
              WRITE(0,*)' Generation of that domain is ',GEN_X
              WRITE(0,*)' This task already has a domain ',MY_DOMAINS_IN_GENS(GEN_X),' in that generation!'
              WRITE(0,*)' Aborting!!'
              CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
            ELSE
              MY_DOMAINS_IN_GENS(GEN_X)=ID_DOM                             !<-- Save the task's domain ID in this generation
            ENDIF
          ENDDO
!
!-----------------------------------------------------------------------
!***  Prepare an array over the generations to indicate when
!***  each task completes its forecast for the domain on which
!***  it lies within each generation.  Recall that in 2-way mode
!***  a task may lie on no more than one domain per generation.
!-----------------------------------------------------------------------
!
          ALLOCATE(GENERATION_FINISHED(1:NUM_GENS))
!
          DO N=1,NUM_GENS
            IF(MY_DOMAINS_IN_GENS(N)>0)THEN
              GENERATION_FINISHED(N)=.FALSE.
            ELSE
              GENERATION_FINISHED(N)=.TRUE.                                !<-- Task not in this generation; consider it finished.
            ENDIF
          ENDDO
!
!-----------------------------------------------------------------------
!***  Prepare an array over the domains as to whether any one of them
!***  will be allowed to integrate a timestep at any given time.
!-----------------------------------------------------------------------
!
          ALLOCATE(FREE_TO_INTEGRATE(1:NUM_GENS))
!
          DO N=1,NUM_GENS
            FREE_TO_INTEGRATE(N)=.TRUE.
          ENDDO
!
        ENDIF
!
!-----------------------------------------------------------------------
!***  The user was required to specify the nest mode in each domain's
!***  configure file indicating whether the parent-child interaction 
!***  will be 1-way or 2-way.  The domains will now extract and save
!***  that specification.
!-----------------------------------------------------------------------
!
        ALLOCATE(nmm_int_state%NEST_MODE(1:NUM_DOMAINS_TOTAL))
!
        DO N=1,NUM_DOMAINS_TOTAL
!
          nmm_int_state%NEST_MODE(N)=' '
!
          IF(MY_DOMAIN_ID_N(N)/=0)THEN                                     !<-- Select tasks on domain #N
            MY_DOMAIN_ID=MY_DOMAIN_ID_N(N)                                 !<-- The ID of domain #N
            COMM_TO_MY_PARENT=>comms_domain(MY_DOMAIN_ID)%TO_PARENT        !<-- The domain's parent-child intracommunicator
! 
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="NMM_INIT: Extract Nest_Mode From Config File"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)           &  !<-- The config object of this domain
                                        ,value =nmm_int_state%NEST_MODE(N) &  !<-- The nest domain's nest_mode
                                        ,label ='nest_mode:'               &  !<-- Give this label's value to the previous variable
                                        ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
          ENDIF
!
        ENDDO
!
!-----------------------------------------------------------------------
!
      ELSE nesting_block_1                                                 !<-- There is only a single domain
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  How many forecast/write tasks will be active on the domain?
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="NMM_INIT: Extract INPES From Config File"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=CF(1)                       &  !<-- The config object
                                    ,value =INPES                       &  !<-- The domain's fcst tasks in I
                                    ,label ='inpes:'                    &  !<-- Give this label's value to the previous variable
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="NMM_INIT: Extract JNPES From Config File"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=CF(1)                       &  !<-- The config object
                                    ,value =JNPES                       &  !<-- The domain's fcst tasks in J
                                    ,label ='jnpes:'                    &  !<-- Give this label's value to the previous variable
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="NMM_INIT: Extract Write_Groups From Config File"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=CF(1)                       &  !<-- The config object
                                    ,value =WRITE_GROUPS                &  !<-- The number of Write groups on this domain
                                    ,label ='write_groups:'             &  !<-- Give this label's value to the previous variable
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="NMM_INIT: Extract Write Tasks Per Group From Config File"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=CF(1)                       &  !<-- The config object
                                    ,value =WRITE_TASKS_PER_GROUP       &  !<-- The number of tasks per Write group
                                    ,label ='write_tasks_per_group:'    &  !<-- Give this label's value to the previous variable
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        IF(.NOT.QUILTING)THEN
          WRITE_GROUPS=0
          WRITE_TASKS_PER_GROUP=0
        ENDIF
!
!-----------------------------------------------------------------------
!
        ALLOCATE(nmm_int_state%NUM_CHILDREN(1))
        nmm_int_state%NUM_CHILDREN(1)=0
!
        N_TASKS=INPES*JNPES+WRITE_GROUPS*WRITE_TASKS_PER_GROUP             !<-- Total # of tasks on the domain
        ALLOCATE(NTASKS_DOMAIN(1))
        NTASKS_DOMAIN(1)=N_TASKS
        ALLOCATE(PETLIST_DOMAIN(1:N_TASKS,1))
!
        ALLOCATE(MY_DOMAINS_IN_GENS(1:1))
        MY_DOMAINS_IN_GENS(1)=1                                            !<-- Dummy value; only relevant for 2-way nests
!
        DO N=1,N_TASKS
          PETLIST_DOMAIN(N,1)=N-1                                          !<-- The list of task IDs for the DOMAIN Component
        ENDDO
!
        ALLOCATE(ID_DOMAINS(1))
        ID_DOMAINS(1)=1                                                    !<-- There is a single domain; its ID is 1
!
        NUM_DOMAINS_MINE=1
        ALLOCATE(MY_DOMAIN_IDS(1),stat=ISTAT)
        MY_DOMAIN_IDS(1)=1
        MY_DOMAIN_ID=1
!
        ALLOCATE(ID_CHILDREN(1,1))
        ID_CHILDREN(1,1)=0                                                 !<-- A single domain thus no children
!
        ALLOCATE(ID_PARENTS(1))
        ID_PARENTS(1)=-999                                                 !<-- There is a single domain; it has no parent
!
        ALLOCATE(COMMS_DOMAIN(1))
        comms_domain(1)%TO_PARENT=-999                                     !<-- There is a single domain; it has no parent
!
        ALLOCATE(FREE_TO_INTEGRATE(1:1))
        FREE_TO_INTEGRATE(1)=.TRUE.                                        !<-- A single domain has no constraints on integration
!
        ALLOCATE(GENERATION_FINISHED(1:1))                                 !<-- A single domain has only one generation to finish
        GENERATION_FINISHED(1)=.FALSE.
!
        ALLOCATE(nmm_int_state%NEST_MODE(1:1))
        nmm_int_state%NEST_MODE(1)=' '
!-----------------------------------------------------------------------
!
      ENDIF nesting_block_1
!
!-----------------------------------------------------------------------
!***  Allocate the DOMAIN import/export states.
!-----------------------------------------------------------------------
!
      ALLOCATE(nmm_int_state%IMP_STATE_DOMAIN(1:NUM_DOMAINS_TOTAL)      &
              ,stat=ISTAT) 
      ALLOCATE(nmm_int_state%EXP_STATE_DOMAIN(1:NUM_DOMAINS_TOTAL)      &
              ,stat=ISTAT) 
!
!-----------------------------------------------------------------------
!***  Create the DOMAIN import/export states.
!-----------------------------------------------------------------------
!
      DO N=1,NUM_DOMAINS_TOTAL
!
        ID_DOM=RANK_TO_DOMAIN_ID(N)
!
        WRITE(INT_TO_CHAR,FMT)ID_DOM
        STATE_NAME='Domain '//INT_TO_CHAR//' Import State'
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="NMM_INIT: Create the DOMAIN Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        nmm_int_state%IMP_STATE_DOMAIN(ID_DOM)=ESMF_StateCreate(           &  !<-- DOMAIN import state
                                       name  =STATE_NAME              &  !<-- DOMAIN import state name
                                      ,stateintent=ESMF_STATEINTENT_IMPORT &
                                      ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        STATE_NAME='Domain '//INT_TO_CHAR//' Export State'
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="NMM_INIT: Create the DOMAIN Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
! 
        nmm_int_state%EXP_STATE_DOMAIN(ID_DOM)=ESMF_StateCreate(           &  !<-- DOMAIN export state
                                       name  =STATE_NAME              &  !<-- DOMAIN export state name
                                      ,stateintent=ESMF_STATEINTENT_EXPORT &
                                      ,rc         =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDDO
!
!-----------------------------------------------------------------------
!***  For the integration we need to know if this is a restarted run
!***  and if we should write a restart file at time 0.
!-----------------------------------------------------------------------
!
      ALLOCATE(nmm_int_state%RESTARTED_RUN(1:NUM_DOMAINS_TOTAL),stat=ISTAT)
      ALLOCATE(nmm_int_state%RST_OUT_00   (1:NUM_DOMAINS_TOTAL),stat=ISTAT)
!
!-----------------------------------------------------------------------
!***  Each task will create Clocks for all domains for simplicity
!***  in executing the major DO loops over the DOMAIN components.
!
!***  Create the domains' clocks with their timesteps, start times,
!***  and run durations.  Also the user-selected task that will
!***  print the clocktimes.
!
!***  Alarms are also created to determine when the integration
!***  of given tasks on each domain reach coupling times.
!-----------------------------------------------------------------------
!
      ALLOCATE(CLOCK_NMM(1:NUM_DOMAINS_TOTAL))
      ALLOCATE(DT(1:NUM_DOMAINS_TOTAL))
      ALLOCATE(FILT_DT(1:NUM_DOMAINS_TOTAL))
!
      ALLOCATE(nmm_int_state%TIMESTEP(1:NUM_DOMAINS_TOTAL))
      ALLOCATE(nmm_int_state%FILT_TIMESTEP(1:NUM_DOMAINS_TOTAL))
!
      ALLOCATE(nmm_int_state%INTERVAL_CLOCKTIME(1:NUM_DOMAINS_TOTAL))
      ALLOCATE(nmm_int_state%INTERVAL_HISTORY  (1:NUM_DOMAINS_TOTAL))
      ALLOCATE(nmm_int_state%INTERVAL_RESTART  (1:NUM_DOMAINS_TOTAL))
!
      ALLOCATE(nmm_int_state%NPE_PRINT(1:NUM_DOMAINS_TOTAL))
!
      ALLOCATE(ALARM_CPL(1:NUM_DOMAINS_TOTAL))
!
!-----------------------------------------------------------------------
!***  Extract timestep information and history/restart output frequency
!***  from the config files of all domains.
!-----------------------------------------------------------------------
!
      timeinfo_loop: DO N=1,NUM_DOMAINS_TOTAL
!
!-----------------------------------------------------------------------
!
        ID_DOM=RANK_TO_DOMAIN_ID(N)
!
        TIMESTEP=>nmm_int_state%TIMESTEP(ID_DOM)
        FILT_TIMESTEP=>nmm_int_state%FILT_TIMESTEP(ID_DOM)
!
        EXP_STATE_DOMAIN=>nmm_int_state%EXP_STATE_DOMAIN(ID_DOM)
!
        INTERVAL_HISTORY=>nmm_int_state%INTERVAL_HISTORY(ID_DOM)
        INTERVAL_RESTART=>nmm_int_state%INTERVAL_RESTART(ID_DOM)
!
        RESTARTED_RUN=>nmm_int_state%RESTARTED_RUN(ID_DOM)
        RST_OUT_00   =>nmm_int_state%RST_OUT_00(ID_DOM)
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="NMM_INIT: Extract Timestep from Config File"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=CF(ID_DOM)                  &  !<-- The config object for this domain
                                    ,value =TIMESTEP_SEC_WHOLE          &  !<-- The variable filled (integer part of timestep (sec))
                                    ,label ='dt_int:'                   &  !<-- Give this label's value to the previous variable
                                    ,rc    =RC)
!
        CALL ESMF_ConfigGetAttribute(config=CF(ID_DOM)                  &  !<-- The config object for this domain 
                                    ,value =TIMESTEP_SEC_NUMERATOR      &  !<-- The variable filled (numerator of timestep fraction)
                                    ,label ='dt_num:'                   &  !<-- Give this label's value to the previous variable
                                    ,rc    =RC)
!
        CALL ESMF_ConfigGetAttribute(config=CF(ID_DOM)                  &  !<-- The config object for this domain 
                                    ,value =TIMESTEP_SEC_DENOMINATOR    &  !<-- The variable filled (denominator of timestep fraction)
                                    ,label ='dt_den:'                   &  !<-- Give this label's value to the previous variable
                                    ,rc    =RC)

        CALL ESMF_ConfigGetAttribute(config=CF(ID_DOM)                  &  !<-- The config object for this domain
                                    ,value =FILT_TIMESTEP_SEC_WHOLE     &  !<-- The variable filled (integer part of timestep (sec))
                                    ,label ='filt_dt_int:'              &  !<-- Give this label's value to the previous variable
                                    ,rc    =RC)
!
        CALL ESMF_ConfigGetAttribute(config=CF(ID_DOM)                  &  !<-- The config object for this domain 
                                    ,value =FILT_TIMESTEP_SEC_NUMERATOR &  !<-- The variable filled (numerator of timestep fraction)
                                    ,label ='filt_dt_num:'              &  !<-- Give this label's value to the previous variable
                                    ,rc    =RC)
!
        CALL ESMF_ConfigGetAttribute(config=CF(ID_DOM)                    &  !<-- The config object for this domain 
                                    ,value =FILT_TIMESTEP_SEC_DENOMINATOR &  !<-- The variable filled (denominator of timestep fraction)
                                    ,label ='filt_dt_den:'                &  !<-- Give this label's value to the previous variable
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Establish the timesteps for all of the domains.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="NMM_INIT: Set Timestep Interval"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_TimeIntervalSet(timeinterval=TIMESTEP                 &  !<-- Fundamental timestep on this domain (sec) (ESMF)
                                 ,s           =TIMESTEP_SEC_WHOLE       &
                                 ,sn          =TIMESTEP_SEC_NUMERATOR   &
                                 ,sd          =TIMESTEP_SEC_DENOMINATOR &
                                 ,rc          =RC)
!
        CALL ESMF_TimeIntervalSet(timeinterval=FILT_TIMESTEP                 &  !<-- Fundamental filter timestep on this domain (sec) (ESMF)
                                 ,s           =FILT_TIMESTEP_SEC_WHOLE       &
                                 ,sn          =FILT_TIMESTEP_SEC_NUMERATOR   &
                                 ,sd          =FILT_TIMESTEP_SEC_DENOMINATOR &
                                 ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        DT(ID_DOM)=TIMESTEP_SEC_WHOLE+REAL(TIMESTEP_SEC_NUMERATOR)      &  !<-- The domain's fundamental timestep (sec) (REAL)
                                     /REAL(TIMESTEP_SEC_DENOMINATOR)
!
        FILT_DT(ID_DOM)=FILT_TIMESTEP_SEC_WHOLE+                        &  !<-- The domain's filter timestep (sec) (REAL)
                   REAL(FILT_TIMESTEP_SEC_NUMERATOR)                    &
                  /REAL(FILT_TIMESTEP_SEC_DENOMINATOR)
!
!-----------------------------------------------------------------------
!***  Get the NMM history and restart output intervals (minutes)
!***  from the config file and save them.  Then make certain that
!***  the fundamental timestep divides evenly into each interval.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Obtain History Interval from the Config File"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=CF(ID_DOM)                  &  !<-- The configure object of this domain
                                    ,value =MINUTES_HISTORY             &  !<-- Fill this variable
                                    ,label ='minutes_history:'          &  !<-- Give the variable this label's value from the config file
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Create the ESMF history file output time interval.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Create the History Output Time Interval."
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_TimeIntervalSet(timeinterval=INTERVAL_HISTORY         &  !<-- Time interval between
                                 ,m           =MINUTES_HISTORY          &  !<-- Minutes between history output
                                 ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Obtain Restart Interval from the Config File"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=CF(ID_DOM)                  &  !<-- The configure object of this domain
                                    ,value =MINUTES_RESTART             &  !<-- Fill this variable
                                    ,label ='minutes_restart:'          &  !<-- Give the variable this label's value from the config file
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Create the ESMF restart file output time interval.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Create the Restart Output Time Interval."
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_TimeIntervalSet(timeinterval=INTERVAL_RESTART         &  !<-- Time interval between restart output (ESMF)
                                 ,m           =MINUTES_RESTART          &  !<-- Minutes between restart output (integer)
                                 ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Set variables related to restarted runs.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="NMM_INIT: Extract Restart Flag from Configure File"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=CF(ID_DOM)                  &  !<-- The config object
                                    ,value =RESTARTED_RUN               &  !<-- Logical flag indicating if this is a restarted run
                                    ,label ='restart:'                  &  !<-- Give this label's value to the previous variable
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="NMM_ATM_DRIVER_INIT: Extract Rst_out_00 Flag from Configure File"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=CF(ID_DOM)                  &  !<-- The config object
                                    ,value =RST_OUT_00                  &  !<-- Should 0-hr history be written for restarted run?
                                    ,label ='rst_out_00:'               &  !<-- Give this label's value to the previous variable
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Ensure that the timestep divides evenly into the restart interval.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="NMM_INIT: Create zero time interval"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_TimeIntervalSet(timeinterval=ZERO_INTERVAL            &
                                 ,s           =0                        &
                                 ,sn          =0                        &
                                 ,sd          =1                        &
                                 ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        IF (MOD(INTERVAL_RESTART,TIMESTEP) /= ZERO_INTERVAL) THEN
          WRITE(0,*)'Timestep of this domain does not divide evenly'    &
                   ,' into the restart interval!'
          WRITE(0,*)' ABORTING!'
          CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
        ENDIF
!
!-----------------------------------------------------------------------
!
      ENDDO timeinfo_loop
!
!-----------------------------------------------------------------------
!***  The coupling interval is the timestep of CLOCK_NEMS.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="NMM_Init: Extract the coupling interval"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ClockGet(clock=CLOCK_NEMS                               &
                        ,timeStep=COUPLING_INTERVAL                     &
                        ,StartTime=Master_StartTime                     &
                        ,StopTime =Master_StopTime                      &  !<-- The end time of the current coupling timestep (ESMF)
                        ,rc=RC)
!
!     write(0,48642)coupling_interval_real
48642 format(' NMM_Init coupling_interval_real=',e12.5)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      call esmf_timeintervalget(timeinterval=coupling_interval &
                       ,s=isecond_fcst  &
                       ,sn=isecond_num  &
                       ,sd=isecond_den &
                       ,rc=rc)
!     write(0,*)' NMM Init timestep of CLOCK_NEMS is ',isecond_fcst &
!              ,'  ',isecond_num,'/',isecond_den,' sec'
!     call esmf_timeget(Master_StartTime, dd=iday_fcst, h=ihour_fcst, m=iminute_fcst, s=isecond_fcst, rc=rc)
!     write(0,*)' NMM Init CLOCK_NEMS starttime is d=',iday_fcst,' h=',ihour_fcst,' m=',iminute_fcst,' s=',isecond_fcst
!     call esmf_timeget(Master_StopTime, dd=iday_fcst, h=ihour_fcst, m=iminute_fcst, s=isecond_fcst, rc=rc)
!     write(0,*)' NMM Init CLOCK_NEMS stoptime is d=',iday_fcst,' h=',ihour_fcst,' m=',iminute_fcst,' s=',isecond_fcst
!-----------------------------------------------------------------------
!***  Obtain the forecast start time from the Main Clock.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="NMM_INIT: Start Time from NMM Clock"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ClockGet(clock      =CLOCK_NEMS                         &  !<-- The NEMS ESMF Clock
                        ,startTime  =STARTTIME                          &  !<-- The simulation start time (ESMF)
!!!                     ,runDuration=RUNDURATION                        &  !<-- The simulation run duration (ESMF)
                        ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
      clock_loop: DO N=1,NUM_DOMAINS_TOTAL
!
        ID_DOM=RANK_TO_DOMAIN_ID(N)
!
        TIMESTEP=>nmm_int_state%TIMESTEP(ID_DOM)
!
!-----------------------------------------------------------------------
!***  Obtain the forecast length time from the configure file.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="NMM_INIT: Extract Forecast Length from Config File"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=CF(ID_DOM)                  &
                                    ,value =NHOURS_FCST                 &
                                    ,label ='nhours_fcst:'              &
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        NSECONDS_FCST=NHOURS_FCST*3600                                     !<-- The forecast length (sec) (REAL)
        TIMESTEP_FINAL=NINT(NSECONDS_FCST/DT(ID_DOM))                      !<-- # of timesteps in the full forecast
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="NMM_INIT: Set the Forecast Length"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_TimeIntervalSet(timeinterval=RUNDURATION              &  !<-- The forecast length (sec) (ESMF)
                                 ,s           =NSECONDS_FCST            &
                                 ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  With data from above, create the ESMF Clocks to control
!***  the timestepping within the DOMAIN subcomponent(s).
!***  Each domain will set its own clock in the initialize
!***  step of DOMAIN_GRID_COMP.
!-----------------------------------------------------------------------
!
        WRITE(INT_TO_CHAR,FMT)ID_DOM
        CLOCK_NMM_NAME='CLOCK_NMM_'//INT_TO_CHAR
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Create the Clocks for the NMM Domains"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CLOCK_NMM(N)=ESMF_ClockCreate(name            =CLOCK_NMM_NAME   &  !<-- The NMM Domain's Clock's name
                                     ,timeStep        =TIMESTEP         &  !<-- The fundamental timestep in this Domain component
                                     ,startTime       =STARTTIME        &  !<-- Start time of simulation
!                                    ,runDuration     =RUNDURATION      &  !<-- Duration of simulation
                                     ,runTimeStepCount=TIMESTEP_FINAL   &  !<-- Length of forecast (timesteps)
                                     ,rc              =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        if(I_AM_ROOT(RC)) then
          CALL NMM_CLOCKPRINT(CLOCK_NMM(N), 'CLOCK_NMM in NMM_INIT', rc=RC)
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
        endif
!
!-----------------------------------------------------------------------
!***  An Alarm is needed to check if each task on a given domain
!***  has reached the end of a coupling interval.  The coupling
!***  interval is the same for all domains.
!-----------------------------------------------------------------------
!
        WRITE(INT_TO_CHAR,FMT)ID_DOM
        ALARM_CPL_NAME='ALARM_CPL_'//INT_TO_CHAR
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Create coupling Alarms for the NMM Domains"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        ALARM_CPL(ID_DOM)=ESMF_AlarmCreate(name             =ALARM_CPL_NAME    &  !<-- The coupling Alarm's name
                                          ,clock            =CLOCK_NMM(ID_DOM) &  !<-- The Alarm is associated with this Clock
                                          ,ringInterval     =COUPLING_INTERVAL &  !<-- Time interval between coupling with ocean
                                          ,ringTimeStepCount=1                 &  !<-- The Alarm rings for this many timesteps
                                          ,sticky           =.false.           &  !<-- Alarm does not ring until turned offs
                                          ,rc               =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!  
!       write(0,23231)id_dom
23231   format(' NMM_Init created ALARM_CPL my_domain_id=',i2)
!     call esmf_alarmget(alarm=ALARM_CPL(ID_DOM) &
!                       ,ringtime=ringtime &
!                       ,ringinterval=ringinterval &
!                       ,rc=rc)
!     call esmf_timeget(time=ringtime &
!                      ,yy=iyear_fcst &
!                      ,mm=imonth_fcst &
!                      ,dd=iday_fcst &
!                      ,h=ihour_fcst &
!                      ,m=iminute_fcst &
!                      ,s=isecond_fcst &
!                      ,sn=isecond_num &
!                      ,sd=isecond_den)
!     write(0,*)' ringtime: y=',iyear_fcst,' mm=',imonth_fcst,' dd=',iday_fcst,' h=',ihour_fcst &
!              ,' m=',iminute_fcst,' s=',isecond_fcst,' sn=',isecond_num,'s d=',isecond_den
!     call esmf_timeintervalget(timeinterval=ringinterval &
!                      ,m=iminute_fcst &
!                      ,s=isecond_fcst &
!                      ,sn=isecond_num &
!                      ,sd=isecond_den)
!     write(0,*)' ringinterval: m=',iminute_fcst,' s=',isecond_fcst,' sn=',isecond_num,'s d=',isecond_den
!-----------------------------------------------------------------------
!***  The fundamental timestep of each NMM domain must divide evenly
!***  into the coupling interval.  If that is not true then abort
!***  the run.
!-----------------------------------------------------------------------
!
        CALL ESMF_TimeIntervalGet(timeinterval=COUPLING_INTERVAL        &
                                 ,s           =ISECOND_FCST             &
                                 ,sn          =ISECOND_NUM              &
                                 ,sd          =ISECOND_DEN              &
                                 ,rc          =RC)
        TIMESTEP_RATIO=(ISECOND_FCST+ISECOND_NUM/ISECOND_DEN)/DT(ID_DOM)
        IF(ABS(TIMESTEP_RATIO-NINT(TIMESTEP_RATIO))>1.E-4)THEN
          WRITE(0,11101)ID_DOM
          WRITE(0,11102)DT(ID_DOM)                                       
          WRITE(0,11103)REAL(ISECOND_FCST)+ISECOND_NUM/ISECOND_DEN
          WRITE(0,*)' Aborting!'
11101     FORMAT(' Timestep of NMM domain ',I2,' does not divide'       &
                   ,' evenly into the coupling interval!!!')
11102     FORMAT(' NMM timestep is ',F8.3,' sec')
11103     FORMAT(' Coupling interval is ',F8.3,' sec')
          CALL ESMF_Finalize(rc=RC,endflag=ESMF_END_ABORT)
        ENDIF
!
!-----------------------------------------------------------------------
!
      ENDDO clock_loop
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Allocate the DOMAIN gridded component(s).
!-----------------------------------------------------------------------
!
      ALLOCATE(nmm_int_state%DOMAIN_GRID_COMP(1:NUM_DOMAINS_TOTAL),stat=ISTAT)
!
      IF(ISTAT/=0)THEN
        WRITE(0,*)' ERROR: Failed to allocate DOMAIN_GRID_COMP'
        WRITE(6,*)' ERROR: Failed to allocate DOMAIN_GRID_COMP'
      ENDIF
!
!-----------------------------------------------------------------------
!***  Allocate other quantities associated with each domain.
!-----------------------------------------------------------------------
!
      ALLOCATE(nmm_int_state%COMM_MY_DOMAIN(1:NUM_DOMAINS_TOTAL))
      ALLOCATE(nmm_int_state%P_C_TIME_RATIO(1:NUM_DOMAINS_TOTAL))
      ALLOCATE(nmm_int_state%MY_DOMAIN_MOVES(1:NUM_DOMAINS_TOTAL))
      ALLOCATE(nmm_int_state%NPHS           (1:NUM_DOMAINS_TOTAL))
      ALLOCATE(nmm_int_state%NTRACK         (1:NUM_DOMAINS_TOTAL))
!
      ALLOCATE(I_AM_A_FCST_TASK)
      ALLOCATE(nmm_int_state%I_AM_A_FCST_TASK(1:NUM_DOMAINS_TOTAL))
      ALLOCATE(nmm_int_state%I_AM_LEAD_FCST_TASK(1:NUM_DOMAINS_TOTAL))
      ALLOCATE(nmm_int_state%I_AM_A_NEST(1:NUM_DOMAINS_TOTAL))
!
      DO N=1,NUM_DOMAINS_TOTAL
!
        nmm_int_state%I_AM_A_FCST_TASK(N)   =.FALSE.
        nmm_int_state%I_AM_LEAD_FCST_TASK(N)=.FALSE.
        nmm_int_state%I_AM_A_NEST(N)        =.FALSE.
!
        nmm_int_state%P_C_TIME_RATIO(N)=0.
        nmm_int_state%MY_DOMAIN_MOVES(N)=.FALSE.
        nmm_int_state%NPHS(N)=0
        nmm_int_state%NTRACK(N)=0
      ENDDO
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  Create the DOMAIN gridded components (one per domain of course).
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      domain_comp_create: DO N=1,NUM_DOMAINS_TOTAL
!
!-----------------------------------------------------------------------
!
        ID_DOM=RANK_TO_DOMAIN_ID(N)
        WRITE(INT_TO_CHAR,FMT)ID_DOM
        DOMAIN_GRID_COMP_NAME=DOMAIN_COMP_BASE//INT_TO_CHAR                !<-- Append domain ID to DOMAIN Comp name
!
        N_TASKS=NTASKS_DOMAIN(ID_DOM)                                      !<-- # of tasks on this domain
        PETLIST=>PETLIST_DOMAIN(1:N_TASKS,ID_DOM)                          !<-- The PETlist for this domain
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="NMM_INIT: Create DOMAIN_GRID_COMP"//INT_TO_CHAR
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        nmm_int_state%DOMAIN_GRID_COMP(ID_DOM)=ESMF_GridCompCreate(     &  !<-- The DOMAIN Component for this domain
                                         name   =DOMAIN_GRID_COMP_NAME  &  !<-- Name of the new DOMAIN gridded component
                                        ,config =CF(ID_DOM)             &  !<-- This domain's configure file
                                        ,petList=PETLIST                &  !<-- The IDs of tasks that will run on this domain
                                        ,rc     =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Register the DOMAIN components' Init, Run, Finalize routines.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Register DOMAIN Init, Run, Finalize"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_GridCompSetServices(nmm_int_state%DOMAIN_GRID_COMP(ID_DOM) &  !<-- The DOMAIN component
                                     ,DOMAIN_REGISTER                        &  !<-- User's subroutineName
                                     ,rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Each task knows on which domain in which generation it lies.  
!***  We are currently considering each domain in the run.  Loop 
!***  through the generations to see if a task lies on this domain.
!***  If it does then it proceeds with initialization action for
!***  that domain otherwise we continue on to the next generation
!***  and check if the task is on this domain in that generation.
!***  If we get all the way through the generations and the task
!***  does not lie on the given domain for any of the generations
!***  then we move on to the next domain and begin the search again
!***  for the task being on that domain.
!-----------------------------------------------------------------------
!
        MY_DOMAIN_ID=0
!
        DO NN=1,NUM_GENS
          IF(ID_DOM==MY_DOMAINS_IN_GENS(NN))THEN                           !<-- Only tasks on domain ID_DOM continue
            MY_DOMAIN_ID=ID_DOM                                            !<-- Yes, this task is on this domain so proceed
            EXIT
          ENDIF
        ENDDO
!
        IF(MY_DOMAIN_ID==0)THEN                                            !<-- Given task not on domain N
          CYCLE domain_comp_create
        ENDIF
!
!-----------------------------------------------------------------------
!***  Insert various quantities into the Domain import state that 
!***  will be needed by that component.
!-----------------------------------------------------------------------
!
        COMM_TO_MY_PARENT=>comms_domain(MY_DOMAIN_ID)%TO_PARENT
        COMM_TO_MY_CHILDREN=>comms_domain(MY_DOMAIN_ID)%to_children
!
        IMP_STATE_DOMAIN=>nmm_int_state%IMP_STATE_DOMAIN(MY_DOMAIN_ID)
        NUM_CHILDREN    =>nmm_int_state%NUM_CHILDREN(MY_DOMAIN_ID)
!
!-----------------------------------------------------------------------
!***  Does this domain move?
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="NMM Init: Extract Move Flag From the Configure file"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)            &
                                    ,value =MY_DOMAIN_MOVES             &
                                    ,label ='my_domain_moves:'          &
                                    ,rc    =rc)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        nmm_int_state%MY_DOMAIN_MOVES(MY_DOMAIN_ID)=MY_DOMAIN_MOVES
!
!-----------------------------------------------------------------------
!***  For hurricane runs we need to know if the storm locator is on
!***  as well as the physics timestep.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="NMM Init: Extract the storm locator flag."
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)            &
                                    ,value =NTRACK                      &
                                    ,label ='ntrack:'                   &
                                    ,rc    =rc)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="NMM Init: Extract the physics timestep."
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)            &
                                    ,value =NPHS                        &
                                    ,label ='nphs:'                     &
                                    ,rc    =rc)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        nmm_int_state%NTRACK(MY_DOMAIN_ID)=NTRACK
        nmm_int_state%NPHS  (MY_DOMAIN_ID)=NPHS
!
!-----------------------------------------------------------------------
!***  Check the configure flag indicating whether or not to run
!***  adiabatically (i.e., with no physics).  Insert the flag
!***  into the DOMAIN import state.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract Adiabatic Flag From the Configure file"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)            &
                                    ,value =MODE                        &
                                    ,label ='adiabatic:'                &
                                    ,rc    =rc)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        IF(TRIM(MODE)=='true')THEN
          PHYSICS_ON = .FALSE.
          IF(MYPE==0) WRITE(0,*)' NMM will run without physics.'
        ELSE
          PHYSICS_ON = .TRUE.
          IF(MYPE==0) WRITE(0,*)' NMM will run with physics.'
        ENDIF

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Add Physics flag to the DOMAIN Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=IMP_STATE_DOMAIN                   &  !<-- This DOMAIN component's import state
                              ,name ='PHYSICS_ON'                       &  !<-- The flag indicating if physics is active
                              ,value=PHYSICS_ON                         &  !<-- The value being inserted into the import state
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Insert the maximum number of domains.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Add NUM_DOMAINS_MAX to the DOMAIN Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=IMP_STATE_DOMAIN                   &  !<-- This DOMAIN component's import state
                              ,name ='MAX_DOMAINS'                      &  !<-- Maximum # of domains
                              ,value=NUM_DOMAINS_MAX                    &  !<-- The scalar being inserted into the import state
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Add NUM_DOMAINS to the DOMAIN Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=IMP_STATE_DOMAIN                   &  !<-- This DOMAIN component's import state
                              ,name ='NUM_DOMAINS'                      &  !<-- # of domains in this run
                              ,value=NUM_DOMAINS_TOTAL                  &  !<-- The scalar being inserted into the import state
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Insert the domain IDs into the DOMAIN import state(s) along with
!***  the association of domains with configure file names,
!***  the number of children and the children's domain IDs.
!***  Also insert a flag as to whether the DOMAIN component is a nest.
!
!***  Note that all tasks are aware of all domains' IDs,
!***  number of children, and those children's domain IDs.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Add Domain IDs to the DOMAIN Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=IMP_STATE_DOMAIN                   &  !<-- This DOMAIN component's import state
                              ,name ='DOMAIN_ID'                        &  !<-- This DOMAIN Component's domain ID
                              ,value=MY_DOMAIN_ID                       &  !<-- The scalar being inserted into the import state
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Insert the association of configure file IDs with domain IDs.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Add Configure File ID Associated With Each Domain ID"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state    =IMP_STATE_DOMAIN                   &  !<-- This DOMAIN component's import state
                              ,name     ='DOMAIN_ID_TO_RANK'                &  !<-- Adding Attribute with this name
                              ,itemCount=NUM_DOMAINS_MAX                    &  !<-- Total # of domains
                              ,valueList=DOMAIN_ID_TO_RANK                  &  !<-- Configure file IDs linked to each domain
                              ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Insert the global rank of the lead task of each domain.  This
!***  is the rank that is retrieved from the VM of each domain.
!-----------------------------------------------------------------------
!
        LEAD_TASK=PETLIST_DOMAIN(1,MY_DOMAIN_ID)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Add Lead Task Rank on Each Domain"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=IMP_STATE_DOMAIN                   &  !<-- This DOMAIN component's import state
                              ,name ='Lead Task Domain'                 &  !<-- Name of Attribute
                              ,value=LEAD_TASK                          &  !<-- Global ran of lead task on this domain
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Add Number of Children to the DOMAIN Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=IMP_STATE_DOMAIN                   &  !<-- This DOMAIN component's import state
                              ,name ='NUM_CHILDREN'                     &  !<-- This DOMAIN Component's # of children
                              ,value=NUM_CHILDREN                       &  !<-- Insert this into the import state
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        IF(COMM_TO_MY_PARENT==-999)THEN
          nmm_int_state%I_AM_A_NEST(ID_DOM)=.FALSE.
        ELSE
          nmm_int_state%I_AM_A_NEST(ID_DOM)=.TRUE.
        ENDIF
!
        I_AM_A_NEST=>nmm_int_state%I_AM_A_NEST(ID_DOM)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Add Nest/Not-a-Nest Flag to the DOMAIN Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=IMP_STATE_DOMAIN                   &  !<-- This DOMAIN component's import state
                              ,name ='I-Am-A-Nest Flag'                 &  !<-- Name of Attribute
                              ,value=I_AM_A_NEST                        &  !<-- Logical nest flag
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Insert the digital filter flag ( >0 indicates which method).
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="NMM_Init: Get Filter Method from Config File"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)            &
                                    ,value =FILTER_METHOD               &
                                    ,label ='filter_method:'            &
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="NMM_Init: Put Filter Method into DOMAIN import state"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=IMP_STATE_DOMAIN                   &  !<-- This DOMAIN component's import state for filter
                              ,name ='Filter_Method'                    &  !<-- Flag for type of digital filter
                              ,value=FILTER_METHOD                      &  !<-- Value of digital filter flag
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
        nesting_block_2: IF(NESTING_NMM)THEN
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Add Domain IDs of Children to the DOMAIN Import State"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          LENGTH=MAX(1,NUM_CHILDREN)             
          CHILD_ID=>ID_CHILDREN(1:LENGTH,ID_DOM)                           !<-- Select only the IDs of this domain's children
!
          CALL ESMF_AttributeSet(state    =IMP_STATE_DOMAIN             &  !<-- This DOMAIN component's import state
                                ,name     ='CHILD_IDs'                  &  !<-- The children's IDs of this DOMAIN Component
                                ,itemCount=LENGTH                       &  !<-- Length of inserted array
                                ,valueList=CHILD_ID                     &  !<-- Insert this into the import state
                                ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          IF(I_AM_A_NEST) THEN
!
            PARENT_CHILD_TIME_RATIO=>nmm_int_state%P_C_TIME_RATIO(ID_DOM)
!
            PARENT_CHILD_TIME_RATIO=NINT(DT(ID_PARENTS(ID_DOM))         &  !<-- Ratio of parent's timestep to this nest's
                                                /DT(ID_DOM))  
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
            MESSAGE_CHECK="Add Parent-Child Time Ratio to DOMAIN Import State"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_AttributeSet(state=IMP_STATE_DOMAIN               &  !<-- This DOMAIN component's import state
                                  ,name ='Parent-Child Time Ratio'      &  !<-- Name of Attribute
                                  ,value=PARENT_CHILD_TIME_RATIO        &  !<-- # of child timesteps per parent timestep
                                  ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
          ENDIF
!
!-----------------------------------------------------------------------
!
        ENDIF nesting_block_2
!
!-----------------------------------------------------------------------
!
      ENDDO domain_comp_create
!
!-----------------------------------------------------------------------
!***  Allocate the clocktime object for all of the domains.  This
!***  holds various timers that will be printed to indicate the
!***  clocktime used by different parts of the code.
!-----------------------------------------------------------------------
!
      ALLOCATE(TIMERS(1:NUM_DOMAINS_TOTAL),stat=ISTAT)
      IF(ISTAT/=0)THEN
        WRITE(0,*)' Failed to allocate TIMERS(1:',NUM_DOMAINS_TOTAL,') in NMM_Init.'
        WRITE(0,*)' Aborting!!'
        CALL ESMF_Finalize(rc=RC,endflag=ESMF_END_ABORT)
      ENDIF
!
      DO N=1,NUM_DOMAINS_TOTAL
!
        timers(n)%total_integ_tim=0.
        timers(n)%totalsum_tim=0
        timers(n)%adv1_tim=0.
        timers(n)%adv2_tim=0.
        timers(n)%bocoh_tim=0.
        timers(n)%bocov_tim=0.
        timers(n)%cdwdt_tim=0.
        timers(n)%cdzdt_tim=0.
        timers(n)%consts_tim=0.
        timers(n)%ddamp_tim=0.
        timers(n)%dht_tim=0.
        timers(n)%exch_tim=0.
        timers(n)%fftfhn_tim=0.
        timers(n)%fftfwn_tim=0.
        timers(n)%hdiff_tim=0.
        timers(n)%mono_tim=0.
        timers(n)%pdtsdt_tim=0.
        timers(n)%pgforce_tim=0.
        timers(n)%poavhn_tim=0.
        timers(n)%polehn_tim=0.
        timers(n)%polewn_tim=0.
        timers(n)%prefft_tim=0.
        timers(n)%presmud_tim=0.
        timers(n)%solver_init_tim=0.
        timers(n)%solver_dyn_tim=0.
        timers(n)%solver_phy_tim=0.
        timers(n)%swaphn_tim=0.
        timers(n)%swapwn_tim=0.
        timers(n)%updatet_tim=0.
        timers(n)%updateuv_tim=0.
        timers(n)%updates_tim=0.
        timers(n)%vsound_tim=0.
        timers(n)%vtoa_tim=0.
        timers(n)%adjppt_tim=0.
        timers(n)%cucnvc_tim=0.
        timers(n)%gsmdrive_tim=0.
        timers(n)%cltend_tim=0.
        timers(n)%rfupdate_tim=0.
        timers(n)%tqadjust_tim=0.
        timers(n)%h_to_v_tim=0.
        timers(n)%gfs_phy_tim=0.
        timers(n)%phy_sum_tim=0.
        timers(n)%pole_swap_tim=0.
        timers(n)%radiation_tim=0.
        timers(n)%rdtemp_tim=0.
        timers(n)%turbl_tim=0.
        timers(n)%domain_run_1=0.
        timers(n)%domain_run_2=0.
        timers(n)%domain_run_3=0.
        timers(n)%pc_cpl_run_cpl1=0.
        timers(n)%pc_cpl_run_cpl2=0.
        timers(n)%pc_cpl_run_cpl3=0.
        timers(n)%pc_cpl_run_cpl4=0.
        timers(n)%pc_cpl_run_cpl5=0.
        timers(n)%cpl1_recv_tim=0.
        timers(n)%cpl2_send_tim=0.
        timers(n)%cpl2_comp_tim=0.
        timers(n)%cpl2_wait_tim=0.
        timers(n)%parent_bookkeep_moving_tim=0.
        timers(n)%parent_update_moving_tim=0.
        timers(n)%t0_recv_move_tim=0.
        timers(n)%update_interior_from_nest_tim=0.
        timers(n)%update_interior_from_parent_tim=0.
!
      ENDDO
!
!-----------------------------------------------------------------------
!***  At this point, DOMAIN components for each domain have been 
!***  created and registered.  Now they need to be initialized.
!
!***  The following call will initialize DOMAIN_GRID_COMP for domain #1.
!***  If more than one domain exists, domain #1 is the uppermost and
!***  the remaining domains will be initialized recursively through
!***  the generations of children.  Recursion is necessary because
!***  children must not be initialized before their parents since
!***  a parent might be directed by the user to generate input data
!***  for its children and that must be complete before the parent's
!***  children are initialized and try to read their input data
!***  before it exists.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!
      CALL CALL_DOMAIN_INITIALIZE(1,CLOCK_NMM)                             !<-- Initiate cascade of DOMAIN Initialize calls for all domains
!
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  We now prepare for work with the coupler components between
!***  parents and their children.
!-----------------------------------------------------------------------
!
      ALLOCATE(nmm_int_state%PC_CPL_COMP(1:NUM_DOMAINS_TOTAL))             !<-- The coupler components.
      ALLOCATE(nmm_int_state%IMP_STATE_PC_CPL(1:NUM_DOMAINS_TOTAL))        !<-- The couplers' import states.
      ALLOCATE(nmm_int_state%EXP_STATE_PC_CPL(1:NUM_DOMAINS_TOTAL))        !<-- The couplers' export states.
!
      ALLOCATE(HANDLE_PACKET_S_H(1:NUM_DOMAINS_TOTAL))                     !<-- Request handles for parent ISends of bndry data packets
      ALLOCATE(HANDLE_PACKET_S_V(1:NUM_DOMAINS_TOTAL))                     !    to children
      ALLOCATE(HANDLE_PACKET_N_H(1:NUM_DOMAINS_TOTAL))                     !
      ALLOCATE(HANDLE_PACKET_N_V(1:NUM_DOMAINS_TOTAL))                     !
      ALLOCATE(HANDLE_PACKET_W_H(1:NUM_DOMAINS_TOTAL))                     !
      ALLOCATE(HANDLE_PACKET_W_V(1:NUM_DOMAINS_TOTAL))                     !
      ALLOCATE(HANDLE_PACKET_E_H(1:NUM_DOMAINS_TOTAL))                     !
      ALLOCATE(HANDLE_PACKET_E_V(1:NUM_DOMAINS_TOTAL))                     !<--
!
      ALLOCATE(HANDLE_I_SW(1:NUM_DOMAINS_TOTAL))                           !<-- Request handle for child ISend of its SW corner to parent
      ALLOCATE(HANDLE_J_SW(1:NUM_DOMAINS_TOTAL))                           !<-- Request handle for child ISend of its SW corner to parent
!
      ALLOCATE(HANDLE_CHILD_LIMITS(1:NUM_DOMAINS_TOTAL))                   !<-- Request handles for parent IRecvs of child task limits
!
      ALLOCATE(HANDLE_CHILD_TOPO_S(1:NUM_DOMAINS_TOTAL))                   !<-- Request handles for parent IRecvs of child bndry topo
      ALLOCATE(HANDLE_CHILD_TOPO_N(1:NUM_DOMAINS_TOTAL))                   !
      ALLOCATE(HANDLE_CHILD_TOPO_W(1:NUM_DOMAINS_TOTAL))                   !
      ALLOCATE(HANDLE_CHILD_TOPO_E(1:NUM_DOMAINS_TOTAL))                   !<--
!
      ALLOCATE(HANDLE_PARENT_DOM_LIMITS(1:NUM_DOMAINS_TOTAL))              !<-- Request handles for ISSends of parent domain limits to 
!                                                                               children.
      ALLOCATE(HANDLE_PARENT_ITE(1:NUM_DOMAINS_TOTAL))                     !<-- Request handles for ISends of parent task limits to children
      ALLOCATE(HANDLE_PARENT_ITS(1:NUM_DOMAINS_TOTAL))                     !
      ALLOCATE(HANDLE_PARENT_JTE(1:NUM_DOMAINS_TOTAL))                     !
      ALLOCATE(HANDLE_PARENT_JTS(1:NUM_DOMAINS_TOTAL))                     !<--
!
      ALLOCATE(PTASK_LIMITS(1:NUM_DOMAINS_TOTAL))                          !<-- Object holding the parent task limits
      ALLOCATE(CTASK_LIMITS(1:NUM_DOMAINS_TOTAL))                          !<-- Object holding parent's children's tasks' limits
!
      ALLOCATE(INFO_SEND(1:NUM_DOMAINS_TOTAL),stat=ISTAT)                  !<-- Parent info to children about which BC updates
      IF(ISTAT/=0)THEN
        WRITE(0,*)' NMM_INIT failed to allocate INFO_SEND stat=',ISTAT
        WRITE(0,*)' Aborting!!'
        CALL ESMF_Finalize(rc=RC,endflag=ESMF_END_ABORT)
      ENDIF
!
      ALLOCATE(nmm_int_state%NUM_2WAY_CHILDREN(1:NUM_DOMAINS_TOTAL))       !<-- Object holding # of 2-way nests on each domain.
      DO N=1,NUM_DOMAINS_TOTAL
        nmm_int_state%NUM_2WAY_CHILDREN(N)=0
      ENDDO
!
!-----------------------------------------------------------------------
!***  Everybody creates an array of Parent-Child couplers for all
!***  the domains.  If there is only a single domain then the
!***  coupler and related variables are empty shells.
!-----------------------------------------------------------------------
!
      pc_cpl_create: DO N=1,NUM_DOMAINS_TOTAL
!
!-----------------------------------------------------------------------
!
        ID_X=RANK_TO_DOMAIN_ID(N)                                      !<-- The domain ID for the Nth domain
!
!-----------------------------------------------------------------------
!***  Create the couplers' import/export states.
!-----------------------------------------------------------------------
!
        IMP_STATE_CPL_NEST=>nmm_int_state%IMP_STATE_PC_CPL(ID_X)
        EXP_STATE_CPL_NEST=>nmm_int_state%EXP_STATE_PC_CPL(ID_X)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Create Empty Import/Export States for Nesting"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        IMP_STATE_CPL_NEST=ESMF_StateCreate(name  ='Nesting Coupler Import' &  !<-- The P-C Coupler import state name
                                           ,stateintent= ESMF_STATEINTENT_IMPORT &
                                           ,rc         =RC)
!
        EXP_STATE_CPL_NEST=ESMF_StateCreate(name  ='Nesting Coupler Export' &  !<-- The P-C Coupler export state name
                                           ,stateintent= ESMF_STATEINTENT_EXPORT &
                                           ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Create the Parent-Child Couplers.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Create the Parent-Child Coupler Component"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        PARENT_CHILD_COUPLER_COMP=>nmm_int_state%PC_CPL_COMP(ID_X)
        PARENT_CHILD_COUPLER_COMP=ESMF_CplCompCreate(name='Parent_Child Coupler' &
                                                    ,rc  =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Register the coupler's Init, Run, and Finalize steps.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Register the Parent-Child Coupler's Init, Run, Finalize"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_CplCompSetServices(cplcomp       =PARENT_CHILD_COUPLER_COMP &  ! <-- The Nesting coupler component
                                    ,userRoutine   =PARENT_CHILD_CPL_REGISTER &  ! <-- The user's subroutineName
                                    ,rc            =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Initialize the request handles for nonblocking sends/recvs in
!***  the gens_1 loop below as well as the object holding those limits.
!-----------------------------------------------------------------------
!
        HANDLE_PACKET_S_H(N)%CHILDREN=>NULL()
        HANDLE_PACKET_S_V(N)%CHILDREN=>NULL()
        HANDLE_PACKET_N_H(N)%CHILDREN=>NULL()
        HANDLE_PACKET_N_V(N)%CHILDREN=>NULL()
        HANDLE_PACKET_W_H(N)%CHILDREN=>NULL()
        HANDLE_PACKET_W_V(N)%CHILDREN=>NULL()
        HANDLE_PACKET_E_H(N)%CHILDREN=>NULL()
        HANDLE_PACKET_E_V(N)%CHILDREN=>NULL()
!
        HANDLE_PARENT_DOM_LIMITS(N)%DATA=>NULL()
!
        HANDLE_PARENT_ITS(N)%DATA=>NULL()
        HANDLE_PARENT_ITE(N)%DATA=>NULL()
        HANDLE_PARENT_JTS(N)%DATA=>NULL()
        HANDLE_PARENT_JTE(N)%DATA=>NULL()
!
        HANDLE_CHILD_LIMITS(N)%CHILDREN=>NULL()
!
        HANDLE_CHILD_TOPO_S(N)%CHILDREN=>NULL()
        HANDLE_CHILD_TOPO_N(N)%CHILDREN=>NULL()
        HANDLE_CHILD_TOPO_W(N)%CHILDREN=>NULL()
        HANDLE_CHILD_TOPO_E(N)%CHILDREN=>NULL()
!
        CTASK_LIMITS(N)%CHILDREN=>NULL()
!
        INFO_SEND(N)%CHILDREN=>NULL()
!
!-----------------------------------------------------------------------
!
      ENDDO pc_cpl_create
!
!-----------------------------------------------------------------------
!***  Each task loops through the generations.  Remember that a given
!***  task can be on no more than one domain in each generation.
!***  This first loop handles the setting up of the Parent-Child
!***  coupler and does preliminary data exchange.  Because that data
!***  exchange includes child-->parent and the parents must use that
!***  data in the upcoming execution of the 1st phase of the coupler
!***  then the setup and the 1st phase must be in their own loops
!***  across the generations.
!-----------------------------------------------------------------------
!
      gens_0: DO NN=1,NUM_GENS
!
!-----------------------------------------------------------------------
!
        MY_DOMAIN_ID=MY_DOMAINS_IN_GENS(NN)                                !<-- This task's (only) domain in generation NN
        IF(MY_DOMAIN_ID==0)CYCLE                                           !<-- This task not on a domain in generation NN
!
!-----------------------------------------------------------------------
!***  Identify the forecast vs. quilt/write tasks since Parent-Child
!***  interaction does not involve any Write tasks.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="NMM_INIT: Extract Fcst-or-Write Flag from Domain Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state=nmm_int_state%EXP_STATE_DOMAIN(MY_DOMAIN_ID) &  !<-- The DOMAIN component export state
                              ,name ='Fcst-or-Write Flag'                         &  !<-- Name of the attribute to extract
                              ,value=I_AM_A_FCST_TASK                             &  !<-- Am I a forecast task?
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        nmm_int_state%I_AM_A_FCST_TASK(MY_DOMAIN_ID)=I_AM_A_FCST_TASK
!
!-----------------------------------------------------------------------
!***  Save this domain's intracommunicator between its forecast tasks.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="NMM_INIT: Extract Fcst Task Intracomm from Domain Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state=nmm_int_state%EXP_STATE_DOMAIN(MY_DOMAIN_ID) &  !<-- The Domain component export state
                              ,name ='Comm Fcst Tasks'                            &  !<-- Name of the attribute to extract
                              ,value=comms_domain(MY_DOMAIN_ID)%TO_FCST_TASKS     &  !<-- Intracommunicator between fcst tasks
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Identify the lead forecast task on each domain.
!-----------------------------------------------------------------------
!
        IF(I_AM_A_FCST_TASK)THEN
!
          CALL MPI_COMM_RANK(comms_domain(MY_DOMAIN_ID)%TO_FCST_TASKS   &  !<-- Intracomm for fcst tasks on this domain
                            ,MYPE_X                                     &  !<-- Rank of this task in the intracommunicator
                            ,IERR)
!
          IF(MYPE_X==0)THEN
            nmm_int_state%I_AM_LEAD_FCST_TASK(MY_DOMAIN_ID)=.TRUE.
          ENDIF
!
        ENDIF
!
!-----------------------------------------------------------------------
!***  If there are nests then create a Parent-Child coupler through
!***  which parents will send boundary data to their children and
!***  also internal data to moving children.
!***  Load that coupler's import state with the data the parents
!***  need to generate boundary data for their children.
!-----------------------------------------------------------------------
!
        nesting_block_3: IF(NESTING_NMM)THEN                               !<-- All parents and children create the Coupler.
!
!-----------------------------------------------------------------------
!
          NUM_CHILDREN=>nmm_int_state%NUM_CHILDREN(MY_DOMAIN_ID)           !<-- How many children does this domain have?
          LENGTH=MAX(1,NUM_CHILDREN)             
          CHILD_ID=>ID_CHILDREN(1:LENGTH,MY_DOMAIN_ID)                     !<-- Select the IDs of this domain's children
!
          COMM_TO_MY_PARENT=>comms_domain(MY_DOMAIN_ID)%TO_PARENT          !<-- This domain's intracommunicator to its parent
          IF(NUM_CHILDREN>0)THEN
            COMM_TO_MY_CHILDREN=>comms_domain(MY_DOMAIN_ID)%TO_CHILDREN    !<-- This domain's intracommunicators to its children
          ELSE
            COMM_TO_MY_CHILDREN=>NULL()
          ENDIF
!
          TIMESTEP=>nmm_int_state%TIMESTEP(MY_DOMAIN_ID)                   !<-- This domain's fundamental timestep
          RESTARTED_RUN=>nmm_int_state%RESTARTED_RUN(MY_DOMAIN_ID)         !<-- Is this a restarted forecast?
!
          DOMAIN_GRID_COMP=>nmm_int_state%DOMAIN_GRID_COMP(MY_DOMAIN_ID)   !<-- This domain's ESMF component
          EXP_STATE_DOMAIN=>nmm_int_state%EXP_STATE_DOMAIN(MY_DOMAIN_ID)   !<-- Its export state
!
          PARENT_CHILD_COUPLER_COMP=>nmm_int_state%PC_CPL_COMP(MY_DOMAIN_ID) !<-- P-C coupler associated with this domain
          IMP_STATE_CPL_NEST=>nmm_int_state%IMP_STATE_PC_CPL(MY_DOMAIN_ID)   !<-- The P-C coupler's import state
          EXP_STATE_CPL_NEST=>nmm_int_state%EXP_STATE_PC_CPL(MY_DOMAIN_ID)   !<-- The P-C coupler's export state
!
!-----------------------------------------------------------------------
!
          CALL PARENT_CHILD_COUPLER_SETUP(NUM_DOMAINS_TOTAL             &  !
                                         ,MY_DOMAIN_ID                  &  !
                                         ,NUM_CHILDREN                  &  !
                                         ,COMM_TO_MY_CHILDREN           &  !
                                         ,COMM_TO_MY_PARENT             &  !
                                         ,DT                            &  !
                                         ,CHILD_ID                      &  !     ^
                                         ,DOMAIN_GRID_COMP              &  !     |
                                         ,EXP_STATE_DOMAIN              &  !     |
                                         ,FTASKS_DOMAIN                 &  !     |
                                         ,NTASKS_DOMAIN                 &  !     |
                                         ,ID_PARENTS                    &  !     |
                                         ,DOMAIN_ID_TO_RANK             &  !     |
                                         ,NUM_DOMAINS_MAX               &  !   Input
!                                                                            ----------
                                         ,IMP_STATE_CPL_NEST            &  !   Output
                                         ,EXP_STATE_CPL_NEST            &  !     |
                                                             )             !     v
!
!-----------------------------------------------------------------------
!
        ENDIF nesting_block_3
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  The tasks in each domain must synchronize before moving to
!***  a different generation.
!-----------------------------------------------------------------------
!
        DOMAIN_GRID_COMP=>nmm_int_state%DOMAIN_GRID_COMP(MY_DOMAIN_ID)     !<-- This domain's ESMF component
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="NMM_INIT: Extract VM for this Domain Component"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_GridCompGet(gridcomp=DOMAIN_GRID_COMP                 &
                             ,vm      =VM                               &  !<-- Get the Virtual Machine for this domain
                             ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="NMM_INIT: Get Intracommunicator for this Domain"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_VMGet(vm             =VM                              &  !<-- The virtual machine
                       ,mpiCommunicator=COMM_MY_DOMAIN                  &  !<-- Intracommunicator for domain MY_DOMAIN_ID
                       ,rc             =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        nmm_int_state%COMM_MY_DOMAIN(MY_DOMAIN_ID)=COMM_MY_DOMAIN
!
        CALL MPI_BARRIER(COMM_MY_DOMAIN,IERR)
!
!-----------------------------------------------------------------------
!
      ENDDO gens_0
!
!-----------------------------------------------------------------------
!***  The forecast tasks now execute phase 1 of the Parent-Child
!***  coupler initialization.
!-----------------------------------------------------------------------
!
      gens_1: DO NN=NUM_GENS,1,-1
!
        MY_DOMAIN_ID=MY_DOMAINS_IN_GENS(NN)                                !<-- This task's (only) domain in generation NN
        IF(MY_DOMAIN_ID==0)CYCLE                                           !<-- This task not on a domain in generation NN
!
        nesting_block_4: IF(NESTING_NMM)THEN                               !<-- All parents and children initialize the Coupler.
!
          PARENT_CHILD_COUPLER_COMP=>nmm_int_state%PC_CPL_COMP(MY_DOMAIN_ID) !<-- P-C coupler associated with this domain
          IMP_STATE_CPL_NEST=>nmm_int_state%IMP_STATE_PC_CPL(MY_DOMAIN_ID)   !<-- The P-C coupler's import state
          EXP_STATE_CPL_NEST=>nmm_int_state%EXP_STATE_PC_CPL(MY_DOMAIN_ID)   !<-- The P-C coupler's export state
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Phase 1 Initialization of the Parent-Child Coupler"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_CplCompInitialize(cplcomp    =PARENT_CHILD_COUPLER_COMP &  !<-- The parent-child coupler component
                                     ,importState=IMP_STATE_CPL_NEST        &  !<-- The parent-child coupler import state
                                     ,exportState=EXP_STATE_CPL_NEST        &  !<-- The parent-child coupler export state
                                     ,clock      =CLOCK_NMM(MY_DOMAIN_ID)   &  !<-- The DOMAIN Clock
                                     ,phase      =1                         &  !<-- The phase (see P-C Register routine)
                                     ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        ENDIF nesting_block_4
!
      ENDDO gens_1
!
!-----------------------------------------------------------------------
!***  The forecast tasks now execute phase 2 of the Parent-Child
!***  coupler initialization.
!-----------------------------------------------------------------------
!
      gens_2: DO NN=1,NUM_GENS
!
        MY_DOMAIN_ID=MY_DOMAINS_IN_GENS(NN)                                !<-- This task's (only) domain in generation NN
        IF(MY_DOMAIN_ID==0)CYCLE                                           !<-- This task not on a domain in generation NN
!
!-----------------------------------------------------------------------
!
        nesting_block_5: IF(NESTING_NMM)THEN                               !<-- All parents and children initialize the Coupler.
!
          DOMAIN_GRID_COMP=>nmm_int_state%DOMAIN_GRID_COMP(MY_DOMAIN_ID)
          I_AM_A_FCST_TASK=nmm_int_state%I_AM_A_FCST_TASK(MY_DOMAIN_ID)
          NUM_CHILDREN=>nmm_int_state%NUM_CHILDREN(MY_DOMAIN_ID)
!
!-----------------------------------------------------------------------
!
          parent_waits_limits: IF(NUM_CHILDREN>0)THEN                      !<-- If so this task is on a parent domain in generation NN
!
!-----------------------------------------------------------------------
!
            CALL ESMF_GridCompGet(gridcomp=DOMAIN_GRID_COMP             &
                                 ,vm      =VM                           &  !<-- Get the Virtual Machine for this domain
                                 ,rc      =RC)
!
            CALL ESMF_VMGet(vm             =VM                          &  !<-- The virtual machine
                           ,localPet       =MYPE_LOCAL                  &  !<-- Rank of task in the domain's intracommunicator
                           ,rc             =RC)
!
!-----------------------------------------------------------------------
!***  Before executing phase 2 of the initialization of the Parent-
!***  Child coupler we will clear the request handles associated with
!***  the nonblocking sends/recvs of the child subdomain limits in 
!***  subroutine PARENT_CHILD_COUPLER_SETUP called in the previous 
!***  loop.  Those sends/recvs consisted of cross-generational
!***  exchanges of data between parents and their children.  Because
!***  only one generation at a time was executed during the loop's
!***  iterations those data exchanges had to be non-blocking.  Now
!***  loop through the generations again, make certain the non-blocking
!***  sends/recvs have finished, then call phase 2 of the Parent-
!***  Child coupler's initialization that includes the use of the
!***  exchanged data.
!-----------------------------------------------------------------------
!
            CHILD_ID=>ID_CHILDREN(1:NUM_CHILDREN,MY_DOMAIN_ID)             !<-- Select the IDs of this domain's children
            ID=MY_DOMAIN_ID
!
            DO N=1,NUM_CHILDREN
              NUM_CHILD_TASKS=FTASKS_DOMAIN(CHILD_ID(N))
!
              DO NT=1,NUM_CHILD_TASKS
!
                IF(MYPE_LOCAL==0)THEN
                  CALL MPI_WAIT(HANDLE_CHILD_LIMITS(ID)%CHILDREN(N)%DATA(NT) &
                               ,JSTAT                                        &
                               ,IERR)
!
                  DO N2=2,FTASKS_DOMAIN(ID)
                    CALL MPI_SEND(CTASK_LIMITS(ID)%CHILDREN(N)%DATA(1,NT) &  !<-- Subdomain limits of child N's task NT
                                 ,4                                       &  !<-- Consists of 4 words
                                 ,MPI_INTEGER                             &  !<-- Data are integers
                                 ,N2-1                                    &  !<-- Send to parent task N2-1
                                 ,N2-1                                    &  !<-- Use target task rank as the MPI tag
                                 ,comms_domain(ID)%TO_FCST_TASKS          &  !<-- Intracomm between parent fcst tasks
                                 ,IERR )
                  ENDDO
!
                ELSE
!
                  IF(I_AM_A_FCST_TASK)THEN
                    CALL MPI_RECV(CTASK_LIMITS(ID)%CHILDREN(N)%DATA(1:4,NT) &  !<-- Subdomain limits of child N's task NT
                                 ,4                                         &  !<-- Consists of 4 words
                                 ,MPI_INTEGER                               &  !<-- Data are integers
                                 ,0                                          &  !<-- Parent task 0 is sending the data
                                 ,MYPE_LOCAL                                &  !<-- Current parent task's local rank
                                 ,comms_domain(ID)%TO_FCST_TASKS            &  !<-- Intracomm between parent fcst tasks
                                 ,JSTAT                                     &
                                 ,IERR )
                  ENDIF
!
                ENDIF
!
              ENDDO
!
            ENDDO
!
!-----------------------------------------------------------------------
!
          ENDIF parent_waits_limits
!
!-----------------------------------------------------------------------
!***  Now we can proceed in executing phase 2 of the Parent-Child
!***  coupler initialization.
!-----------------------------------------------------------------------
!
          PARENT_CHILD_COUPLER_COMP=>nmm_int_state%PC_CPL_COMP(MY_DOMAIN_ID) !<-- P-C coupler associated with this domain
          IMP_STATE_CPL_NEST=>nmm_int_state%IMP_STATE_PC_CPL(MY_DOMAIN_ID)   !<-- The P-C coupler's import state
          EXP_STATE_CPL_NEST=>nmm_int_state%EXP_STATE_PC_CPL(MY_DOMAIN_ID)   !<-- The P-C coupler's export state
!
          I_AM_A_FCST_TASK=nmm_int_state%I_AM_A_FCST_TASK(MY_DOMAIN_ID)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Phase 2 Initialization of the Parent-Child Coupler"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_CplCompInitialize(cplcomp    =PARENT_CHILD_COUPLER_COMP &  !<-- The parent-child coupler component
                                     ,importState=IMP_STATE_CPL_NEST        &  !<-- The parent-child coupler import state
                                     ,exportState=EXP_STATE_CPL_NEST        &  !<-- The parent-child coupler export state
                                     ,clock      =CLOCK_NMM(MY_DOMAIN_ID)   &  !<-- The DOMAIN Clock
                                     ,phase      =2                         &  !<-- The phase (see P-C Register routine)
                                     ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  The number of 2-way children is a required argument in the call
!***  to NMM_INTEGRATE.  Extract its value from the export state of
!***  the P-C coupler.
!-----------------------------------------------------------------------
!
          NUM_2WAY_CHILDREN=>nmm_int_state%NUM_2WAY_CHILDREN(MY_DOMAIN_ID)
!
          IF(I_AM_A_FCST_TASK.AND.NUM_CHILDREN>0)THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="NMM_Init: Extract # of 2-Way Children from P-C Export State"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_AttributeGet(state=EXP_STATE_CPL_NEST             &  !<-- The Parent-Child coupler export state
                                  ,name ='NUM_2WAY_CHILDREN'            &  !<-- Name of the attribute to extract
                                  ,value=NUM_2WAY_CHILDREN              &  !<-- How many 2-way children in the current domain?
                                  ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          ENDIF
!
!-----------------------------------------------------------------------
!
        ENDIF nesting_block_5
!
!-----------------------------------------------------------------------
!
        INTERVAL_CLOCKTIME=>nmm_int_state%INTERVAL_CLOCKTIME(MY_DOMAIN_ID)
        NPE_PRINT=>nmm_int_state%NPE_PRINT(MY_DOMAIN_ID)
!
!-----------------------------------------------------------------------
!***  Extract ID of the task that will print clocktimes on this domain.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Read MPI Task ID That Provides Clocktime Output"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)            &  !<-- The configure object
                                    ,value =NPE_PRINT                   &  !<-- Fill this variable (this task prints its clocktimes)
                                    ,label ='npe_print:'                &  !<-- Give the variable this label's value from the config file
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Get print_timing flag from config file.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Get print_timing flag from configure file"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)            &  !<-- The configure object
                                    ,value =PRINT_TIMING                &  !<-- Fill this variable (this task prints its clocktimes)
                                    ,label ='print_timing:'             &  !<-- Give the variable this label's value from the config file
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Set the forecast time interval (sec) between writes of the
!***  clocktime statistics by the task specified in the configure
!***  file.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Read Fcst Interval for Clocktime Output"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)            &  !<-- The configure object
                                    ,value =NHOURS_CLOCKTIME            &  !<-- Fill this variable (fcst hrs between clocktime prints)
                                    ,label ='nhours_clocktime:'         &  !<-- Give the variable this label's value from the config file
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Create ESMF Clocktime Output Interval"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_TimeIntervalSet(timeinterval=INTERVAL_CLOCKTIME       &  !<-- Time interval between clocktime writes (h) (ESMF)
                                 ,h           =NHOURS_CLOCKTIME         &  !<-- Hours between clocktime writes (INTEGER)
                                 ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  The tasks in each domain must synchronize before moving to
!***  a different generation.
!-----------------------------------------------------------------------
!
        COMM_MY_DOMAIN=nmm_int_state%COMM_MY_DOMAIN(MY_DOMAIN_ID)
!
        CALL MPI_BARRIER(COMM_MY_DOMAIN,IERR)
!
!-----------------------------------------------------------------------
!
      ENDDO gens_2
!
!-----------------------------------------------------------------------
!
      gens_3: DO NN=1,NUM_GENS
!
!-----------------------------------------------------------------------
!***  Before executing phase 3 of the initialization of the Parent-
!***  Child coupler we will clear the request handles associated with
!***  the nonblocking sends/recvs of the child topography in phase 1
!***  of the initialization.  Those sends/recvs consisted of cross-
!***  generational exchanges of data between parents and their
!***  children.  Because only one generation at a time is executed
!***  during the loop's iterations those data exchanges had to be
!***  non-blocking.  Now loop through the generations again, make
!***  certain the non-blocking sends/recvs have finished, then call
!***  phase 3 of the Parent-Child coupler's initialization that
!***  includes the use of the exchanged data.
!-----------------------------------------------------------------------
!
        MY_DOMAIN_ID=MY_DOMAINS_IN_GENS(NN)                                !<-- This task's (only) domain in generation NN
        IF(MY_DOMAIN_ID==0)CYCLE                                           !<-- This task not on a domain in generation NN
!
        I_AM_A_FCST_TASK=nmm_int_state%I_AM_A_FCST_TASK(MY_DOMAIN_ID)
        NUM_CHILDREN=>nmm_int_state%NUM_CHILDREN(MY_DOMAIN_ID)
!
!-----------------------------------------------------------------------
!
        parent_waits_topo: IF(NUM_CHILDREN>0                            &  !<-- If so this task is on a parent domain in generation NN
                                 .AND.                                  &
                              I_AM_A_FCST_TASK)THEN
!
!-----------------------------------------------------------------------
!
          DOMAIN_GRID_COMP=>nmm_int_state%DOMAIN_GRID_COMP(MY_DOMAIN_ID)
!
          CALL ESMF_GridCompGet(gridcomp=DOMAIN_GRID_COMP               &
                               ,vm      =VM                             &  !<-- Get the Virtual Machine for this domain
                               ,rc      =RC)
!
          CALL ESMF_VMGet(vm             =VM                            &  !<-- The virtual machine
                         ,localPet       =MYPE_LOCAL                    &  !<-- Rank of task in the domain's intracommunicator
                         ,rc             =RC)
!
!-----------------------------------------------------------------------
!***  Make sure the relevant parent tasks received boundary topography
!***  from their children.
!-----------------------------------------------------------------------
!
          ID=MY_DOMAIN_ID
!
          IF(ASSOCIATED(HANDLE_CHILD_TOPO_S(ID)%CHILDREN))THEN
!
            DO N=1,NUM_CHILDREN
              IF(ASSOCIATED(HANDLE_CHILD_TOPO_S(ID)%CHILDREN(N)%DATA))THEN
                INDX2=UBOUND(HANDLE_CHILD_TOPO_S(ID)%CHILDREN(N)%DATA,1)
                DO NT=1,INDX2
                  CALL MPI_WAIT(HANDLE_CHILD_TOPO_S(ID)%CHILDREN(N)%DATA(NT) &
                               ,JSTAT                                        &
                               ,IERR)
                ENDDO
              ENDIF
            ENDDO
          ENDIF
!
          IF(ASSOCIATED(HANDLE_CHILD_TOPO_N(ID)%CHILDREN))THEN
!
            DO N=1,NUM_CHILDREN
              IF(ASSOCIATED(HANDLE_CHILD_TOPO_N(ID)%CHILDREN(N)%DATA))THEN
                INDX2=UBOUND(HANDLE_CHILD_TOPO_N(ID)%CHILDREN(N)%DATA,1)
                DO NT=1,INDX2
                  CALL MPI_WAIT(HANDLE_CHILD_TOPO_N(ID)%CHILDREN(N)%DATA(NT) &
                               ,JSTAT                                        &
                               ,IERR)
                ENDDO
              ENDIF
            ENDDO
          ENDIF
!
          IF(ASSOCIATED(HANDLE_CHILD_TOPO_W(ID)%CHILDREN))THEN
!
            DO N=1,NUM_CHILDREN
              IF(ASSOCIATED(HANDLE_CHILD_TOPO_W(ID)%CHILDREN(N)%DATA))THEN
                INDX2=UBOUND(HANDLE_CHILD_TOPO_W(ID)%CHILDREN(N)%DATA,1)
                DO NT=1,INDX2
                  CALL MPI_WAIT(HANDLE_CHILD_TOPO_W(ID)%CHILDREN(N)%DATA(NT) &
                               ,JSTAT                                        &
                               ,IERR)
                ENDDO
              ENDIF
            ENDDO
          ENDIF
!
          IF(ASSOCIATED(HANDLE_CHILD_TOPO_E(ID)%CHILDREN))THEN
!
            DO N=1,NUM_CHILDREN
              IF(ASSOCIATED(HANDLE_CHILD_TOPO_E(ID)%CHILDREN(N)%DATA))THEN
                INDX2=UBOUND(HANDLE_CHILD_TOPO_E(ID)%CHILDREN(N)%DATA,1)
                DO NT=1,INDX2
                  CALL MPI_WAIT(HANDLE_CHILD_TOPO_E(ID)%CHILDREN(N)%DATA(NT) &
                               ,JSTAT                                        &
                               ,IERR)
                ENDDO
              ENDIF
            ENDDO
          ENDIF
!
!
          DO N=1,NUM_CHILDREN
!
            IF(ASSOCIATED(HANDLE_PARENT_DOM_LIMITS(ID)%DATA))THEN
              CALL MPI_WAIT(HANDLE_PARENT_DOM_LIMITS(ID)%DATA(N)        &
                           ,JSTAT                                       &
                           ,IERR)
            ENDIF
!
            IF(ASSOCIATED(HANDLE_PARENT_ITS(ID)%DATA))THEN
              CALL MPI_WAIT(HANDLE_PARENT_ITS(ID)%DATA(N)               &
                           ,JSTAT                                       &
                           ,IERR)
!
              CALL MPI_WAIT(HANDLE_PARENT_ITE(ID)%DATA(N)               &
                           ,JSTAT                                       &
                           ,IERR)
!
              CALL MPI_WAIT(HANDLE_PARENT_JTS(ID)%DATA(N)               &
                           ,JSTAT                                       &
                           ,IERR)
!
              CALL MPI_WAIT(HANDLE_PARENT_JTE(ID)%DATA(N)               &
                           ,JSTAT                                       &
                           ,IERR)
            ENDIF
!
          ENDDO
!
        ENDIF parent_waits_topo
!
!-----------------------------------------------------------------------
!***  Clear the request handles for the parents' ISends of the 
!***  boundary info packets in phase 1 of the Init step and
!***  deallocate memory we are finished with.
!-----------------------------------------------------------------------
!
        I_AM_A_FCST_TASK=nmm_int_state%I_AM_A_FCST_TASK(MY_DOMAIN_ID)
!
        parent_waits_bc_info: IF(NUM_CHILDREN>0.AND.                    &  !<-- Select fcst tasks on all the parents
                                 I_AM_A_FCST_TASK)THEN
!
          ID=MY_DOMAIN_ID
!
          IF(ASSOCIATED(HANDLE_PACKET_S_H(ID)%CHILDREN))THEN
!
            DO N=1,NUM_CHILDREN
              IF(ASSOCIATED(HANDLE_PACKET_S_H(ID)%CHILDREN(N)%DATA))THEN
                INDX2=UBOUND(HANDLE_PACKET_S_H(ID)%CHILDREN(N)%DATA,1)
                DO NT=1,INDX2
                  CALL MPI_WAIT(HANDLE_PACKET_S_H(ID)%CHILDREN(N)%DATA(NT) &
                               ,JSTAT                                      &
                               ,IERR)
                ENDDO
!
              ENDIF
            ENDDO
          ENDIF
!
          IF(ASSOCIATED(HANDLE_PACKET_S_V(ID)%CHILDREN))THEN
!
            DO N=1,NUM_CHILDREN
              IF(ASSOCIATED(HANDLE_PACKET_S_V(ID)%CHILDREN(N)%DATA))THEN
                INDX2=UBOUND(HANDLE_PACKET_S_V(ID)%CHILDREN(N)%DATA,1)
                DO NT=1,INDX2
                  CALL MPI_WAIT(HANDLE_PACKET_S_V(ID)%CHILDREN(N)%DATA(NT) &
                               ,JSTAT                                      &
                               ,IERR)
                ENDDO
!
              ENDIF
            ENDDO
          ENDIF
!
          IF(ASSOCIATED(HANDLE_PACKET_N_H(ID)%CHILDREN))THEN
!
            DO N=1,NUM_CHILDREN
              IF(ASSOCIATED(HANDLE_PACKET_N_H(ID)%CHILDREN(N)%DATA))THEN
                INDX2=UBOUND(HANDLE_PACKET_N_H(ID)%CHILDREN(N)%DATA,1)
                DO NT=1,INDX2
                  CALL MPI_WAIT(HANDLE_PACKET_N_H(ID)%CHILDREN(N)%DATA(NT) &
                               ,JSTAT                                      &
                               ,IERR)
                ENDDO
!
              ENDIF
            ENDDO
          ENDIF
!
          IF(ASSOCIATED(HANDLE_PACKET_N_V(ID)%CHILDREN))THEN
!
            DO N=1,NUM_CHILDREN
              IF(ASSOCIATED(HANDLE_PACKET_N_V(ID)%CHILDREN(N)%DATA))THEN
                INDX2=UBOUND(HANDLE_PACKET_N_V(ID)%CHILDREN(N)%DATA,1)
                DO NT=1,INDX2
                  CALL MPI_WAIT(HANDLE_PACKET_N_V(ID)%CHILDREN(N)%DATA(NT) &
                               ,JSTAT                                      &
                               ,IERR)
                ENDDO
!
              ENDIF
            ENDDO
          ENDIF
!
          IF(ASSOCIATED(HANDLE_PACKET_W_H(ID)%CHILDREN))THEN
!
            DO N=1,NUM_CHILDREN
              IF(ASSOCIATED(HANDLE_PACKET_W_H(ID)%CHILDREN(N)%DATA))THEN
                INDX2=UBOUND(HANDLE_PACKET_W_H(ID)%CHILDREN(N)%DATA,1)
                DO NT=1,INDX2
                  CALL MPI_WAIT(HANDLE_PACKET_W_H(ID)%CHILDREN(N)%DATA(NT) &
                               ,JSTAT                                      &
                               ,IERR)
                ENDDO
!
              ENDIF
            ENDDO
          ENDIF
!
          IF(ASSOCIATED(HANDLE_PACKET_W_V(ID)%CHILDREN))THEN
!
            DO N=1,NUM_CHILDREN
              IF(ASSOCIATED(HANDLE_PACKET_W_V(ID)%CHILDREN(N)%DATA))THEN
                INDX2=UBOUND(HANDLE_PACKET_W_V(ID)%CHILDREN(N)%DATA,1)
                DO NT=1,INDX2
                  CALL MPI_WAIT(HANDLE_PACKET_W_V(ID)%CHILDREN(N)%DATA(NT) &
                               ,JSTAT                                      &
                               ,IERR)
                ENDDO
!
              ENDIF
            ENDDO
          ENDIF
!
          IF(ASSOCIATED(HANDLE_PACKET_E_H(ID)%CHILDREN))THEN
!
            DO N=1,NUM_CHILDREN
              IF(ASSOCIATED(HANDLE_PACKET_E_H(ID)%CHILDREN(N)%DATA))THEN
                INDX2=UBOUND(HANDLE_PACKET_E_H(ID)%CHILDREN(N)%DATA,1)
                DO NT=1,INDX2
                  CALL MPI_WAIT(HANDLE_PACKET_E_H(ID)%CHILDREN(N)%DATA(NT) &
                               ,JSTAT                                      &
                               ,IERR)
                ENDDO
!
              ENDIF
            ENDDO
          ENDIF
!
          IF(ASSOCIATED(HANDLE_PACKET_E_V(ID)%CHILDREN))THEN
!
            DO N=1,NUM_CHILDREN
              IF(ASSOCIATED(HANDLE_PACKET_E_V(ID)%CHILDREN(N)%DATA))THEN
                INDX2=UBOUND(HANDLE_PACKET_E_V(ID)%CHILDREN(N)%DATA,1)
                DO NT=1,INDX2
                  CALL MPI_WAIT(HANDLE_PACKET_E_V(ID)%CHILDREN(N)%DATA(NT) &
                               ,JSTAT                                      &
                               ,IERR)
                ENDDO
!
              ENDIF
            ENDDO
          ENDIF
!
        ENDIF parent_waits_bc_info
!
!-----------------------------------------------------------------------
!***  The tasks in each domain must synchronize before moving to
!***  a different generation.
!-----------------------------------------------------------------------
!
        COMM_MY_DOMAIN=nmm_int_state%COMM_MY_DOMAIN(MY_DOMAIN_ID)
!
        CALL MPI_BARRIER(COMM_MY_DOMAIN,IERR)
!
!-----------------------------------------------------------------------
!***  The forecast tasks now execute phase 3 of the Parent-Child
!***  coupler initialization.
!-----------------------------------------------------------------------
!
        nesting_block_6: IF(NESTING_NMM)THEN                               !<-- All parents and children create the Coupler.
!
!-----------------------------------------------------------------------
!
          IMP_STATE_CPL_NEST=>nmm_int_state%IMP_STATE_PC_CPL(MY_DOMAIN_ID)
          EXP_STATE_CPL_NEST=>nmm_int_state%EXP_STATE_PC_CPL(MY_DOMAIN_ID)
          PARENT_CHILD_COUPLER_COMP=>nmm_int_state%PC_CPL_COMP(MY_DOMAIN_ID)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Phase 3 Initialization of Parent-Child Coupler"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_CplCompInitialize(cplcomp    =PARENT_CHILD_COUPLER_COMP &  !<-- The parent-child coupler component
                                     ,importState=IMP_STATE_CPL_NEST        &  !<-- The parent-child coupler import state
                                     ,exportState=EXP_STATE_CPL_NEST        &  !<-- The parent-child coupler export state
                                     ,clock      =CLOCK_NMM(MY_DOMAIN_ID)   &  !<-- The DOMAIN Clock
                                     ,phase      =3                         &  !<-- The phase (see P-C Register routine)
                                     ,rc         =RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        ENDIF nesting_block_6
!
!-----------------------------------------------------------------------
!
      ENDDO gens_3
!
!-----------------------------------------------------------------------
!***  The central lat/lon of the single domain / upper parent are
!***  needed.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!     MESSAGE_CHECK="NMM_INIT: Extract Total Domain Count Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF(1)                         &  !<-- The config object
                                  ,value =TPH0D                         &  !<-- Central geographic lat (deg) of uppermost domain.
                                  ,label ='tph0d:'                      &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
      CALL ESMF_ConfigGetAttribute(config=CF(1)                         &  !<-- The config object
                                  ,value =TLM0D                         &  !<-- Central geographic lon (deg) of uppermost domain.
                                  ,label ='tlm0d:'                      &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!     CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Save the values needed for NUOPC coupling to external models.
!-----------------------------------------------------------------------
!
      CALL STORE_DOMAIN_DESCRIPTORS
!
!-----------------------------------------------------------------------
!***  Create the ESMF grid.  See comments in FUNCTION NMMB_GridCreate.
!***  For multiple NMM domains each task must tell NMMB_GridCreate
!***  which domain it is on.  
!***  Announce the fields that will be part of the coupling.
!-----------------------------------------------------------------------
!
      DO N=1,NUM_DOMAINS_TOTAL
!
        ID_X=RANK_TO_DOMAIN_ID(N)                                      !<-- The domain ID for the Nth domain
        MY_DOMAIN_ID=ID_X
!
!--------------------------------------------------
!***  Extract the current domain's internal state.
!--------------------------------------------------
!
        IF(DOMAIN_DESCRIPTORS(MY_DOMAIN_ID)%TASK_ACTIVE)THEN
          CALL ESMF_GridCompGetInternalState(&
            nmm_int_state%DOMAIN_GRID_COMP(MY_DOMAIN_ID), &
            WRAP_DOMAIN, RC)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
          DOMAIN_INT_STATE=>wrap_domain%DOMAIN_INT_STATE
        ELSE
          DOMAIN_INT_STATE=>NULL()
        ENDIF
!
        pGrid_NMMB = NMMB_GridCreate(N, nmm_int_state%DOMAIN_GRID_COMP(MY_DOMAIN_ID), DOMAIN_DESCRIPTORS, TPH0D,TLM0D, RC=RC)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
!
        DOMAIN_DESCRIPTORS(N)%GRID = pGrid_NMMB
        DOMAIN_DESCRIPTORS(N)%PARENT_DOMAIN_ID = 1
!
        CALL NMMB_CreateRouteHandle(N, DOMAIN_DESCRIPTORS, RC=RC)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
!
        CALL NMMB_CreateDomainFields(N, DOMAIN_DESCRIPTORS, RC=RC)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
!
      ENDDO

      nmm_grid = DOMAIN_DESCRIPTORS(1)%GRID

!-----------------------------------------------------------------------
!***  Only some of the full list of potential export fields will
!***  be regridded from the nests to the upper parent.  Find those
!***  fields' indices in the full list of export fields and save them.
!-----------------------------------------------------------------------
!
      DO N1=1,nExportFields_NMMB
        N2=queryFieldList(exportFieldsList, EXPORT_FIELDS_BLEND(N1), rc=rc)
        EXPORT_FIELDS_INDX(N1)=N2
      ENDDO
!
!-----------------------------------------------------------------------
!
      IF(RC_INIT==ESMF_SUCCESS)THEN
!       WRITE(0,*)' NMM_INITIALIZE succeeded'
      ELSE
        WRITE(0,*)' NMM_INITIALIZE failed  RC_INIT=',RC_INIT
      ENDIF
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!
      SUBROUTINE STORE_DOMAIN_DESCRIPTORS
!
!-----------------------------------------------------------------------
!***  Store values into the derived-type object needed for coupling
!***  domains to external models.  This is an internal subroutine
!***  in NMM_INITIALIZE.
!-----------------------------------------------------------------------
!
      USE module_DOMAIN_INTERNAL_STATE,ONLY: DOMAIN_INTERNAL_STATE      &
                                            ,WRAP_DOMAIN_INTERNAL_STATE

      USE module_SOLVER_INTERNAL_STATE,ONLY: SOLVER_INTERNAL_STATE      &
                                            ,WRAP_SOLVER_INT_STATE
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!     
!-----------------------------------------------------------------------
!
!---------------------  
!***  Local variables
!---------------------  
!
      INTEGER(kind=KINT) :: I,ITE,ITS,J,JTE,JTS,NUM_PETS,PET_N,RC
      INTEGER(kind=KINT),DIMENSION(2) :: minIndex, maxIndex
!
      REAL(kind=KDBL) :: DEG2RAD,LAM0,PHI0,PI
!
      REAL(kind=KDBL),DIMENSION(2) :: CENTRAL_LATLON
!
      CHARACTER(4096)    :: TMPSTR
!
      TYPE(DOMAIN_INTERNAL_STATE),POINTER :: DOMAIN_INT_STATE
      TYPE(WRAP_DOMAIN_INTERNAL_STATE) :: WRAP_DOMAIN
      TYPE(SOLVER_INTERNAL_STATE),POINTER :: SOLVER_INT_STATE
      TYPE(WRAP_SOLVER_INT_STATE) :: WRAP_SOLVER
      TYPE(ESMF_VM)               :: VM
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Allocate the object that will store fundamental specifications
!***  of the domains and MPI task subdomains.  This object will be
!***  used for coupling to external models.
!-----------------------------------------------------------------------
!
      ALLOCATE(DOMAIN_DESCRIPTORS(1:NUM_DOMAINS_TOTAL),stat=ISTAT)
!
      IF(ISTAT/=0)THEN
        WRITE(0,*)' Failed to allocate NMMB_DOMAIN_DESCRIPTORS'
        WRITE(0,*)' stat=',ISTAT
        WRITE(0,*)' ABORTING!'
        CALL ESMF_Finalize(rc=RC,endflag=ESMF_END_ABORT)
      ENDIF

      call ESMF_VMGetCurrent(vm, rc=RC)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
!
!-----------------------------------------------------------------------
!
      PI=ACOS(-1._KDBL)
      DEG2RAD=PI/180._KDBL
!
!-----------------------------------------------------------------------
!
      domains: DO N=1,NUM_DOMAINS_TOTAL
!
!-----------------------------------------------------------------------
!***  First save the domains' compute task layouts.
!-----------------------------------------------------------------------
!
        ID_X=RANK_TO_DOMAIN_ID(N)                                          !<-- The domain ID for the Nth domain
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="NMM_INIT: Extract generation from Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=CF(ID_X)                    &  !<-- The config object
                                    ,value =INPES                       &  !<-- # of tasks in domain N in I direction
                                    ,label ='inpes:'                    &  !<-- Give this label's value to the previous variable
                                    ,rc    =RC)
!
        CALL ESMF_ConfigGetAttribute(config=CF(ID_X)                    &  !<-- The config object
                                    ,value =JNPES                       &  !<-- # of tasks in domain N in J direction
                                    ,label ='jnpes:'                    &  !<-- Give this label's value to the previous variable
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        DOMAIN_DESCRIPTORS(ID_X)%INPES=INPES
        DOMAIN_DESCRIPTORS(ID_X)%JNPES=JNPES
!
!-----------------------------------------------------------------------
!***  Save the IDs of all compute tasks on the domain.  Also save
!***  a flag indicating whether or not the current task is on the
!***  given domain.
!-----------------------------------------------------------------------
!
        NUM_PETS=INPES*JNPES
        DOMAIN_DESCRIPTORS(ID_X)%NUM_PETS=NUM_PETS
!
        ALLOCATE(DOMAIN_DESCRIPTORS(ID_X)%PET_MAP(1:NUM_PETS),stat=ISTAT)
!      
        DOMAIN_DESCRIPTORS(ID_X)%TASK_ACTIVE=.FALSE.
!
        DO N1=1,NUM_PETS
          PET_N=PETLIST_DOMAIN(N1,ID_X)
          DOMAIN_DESCRIPTORS(ID_X)%PET_MAP(N1)=PET_N                       !<-- The task IDs on domain ID_X relative to all tasks.
          IF(MYPE==PET_N)THEN
            DOMAIN_DESCRIPTORS(ID_X)%TASK_ACTIVE=.TRUE.                    !<-- The current task lies on domain ID_X.
          ENDIF
        ENDDO
!
!-----------------------------------------------------------------------
!***  We need to extract the Solver component's internal state
!***  in order to obtain this domain's full coordinate limits.
!-----------------------------------------------------------------------
!
        if(DOMAIN_DESCRIPTORS(ID_X)%TASK_ACTIVE) then
          CALL ESMF_GridCompGetInternalState(nmm_int_state%DOMAIN_GRID_COMP(ID_X) &                                                 
                                            ,WRAP_DOMAIN                          & 
                                            ,RC)
!
          DOMAIN_INT_STATE=>wrap_domain%DOMAIN_INT_STATE
!
          CALL ESMF_GridCompGetInternalState(domain_int_state%SOLVER_GRID_COMP &
                                          ,WRAP_SOLVER                       &
                                          ,RC)
!
          SOLVER_INT_STATE=>wrap_solver%INT_STATE
!
          minIndex(1) = solver_int_state%IDS
          minIndex(2) = solver_int_state%JDS
          maxIndex(1) = solver_int_state%IDE
          maxIndex(2) = solver_int_state%JDE
        endif

        ! Broadcast min/max Index from first active PET to all PETs on NMMB_GRID_COMP
        call ESMF_VMBroadcast(vm, minIndex, 2, PETLIST_DOMAIN(1,ID_X), rc=RC)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        call ESMF_VMBroadcast(vm, maxIndex, 2, PETLIST_DOMAIN(1,ID_X), rc=RC)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out

        DOMAIN_DESCRIPTORS(ID_X)%INDX_MIN(1)=minIndex(1)
        DOMAIN_DESCRIPTORS(ID_X)%INDX_MIN(2)=minIndex(2)
        DOMAIN_DESCRIPTORS(ID_X)%INDX_MAX(1)=maxIndex(1)
        DOMAIN_DESCRIPTORS(ID_X)%INDX_MAX(2)=maxIndex(2)

        write(tmpstr, *) "DOMAIN ID_X = ", ID_X
        call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=rc)
        write(tmpstr, *) "INPES = ", DOMAIN_DESCRIPTORS(ID_X)%INPES
        call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=rc)
        write(tmpstr, *) "JNPES = ", DOMAIN_DESCRIPTORS(ID_X)%JNPES
        call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=rc)
        write(tmpstr, *) "MIN_INDEX = ", DOMAIN_DESCRIPTORS(ID_X)%INDX_MIN
        call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=rc)
        write(tmpstr, *) "MAX_INDEX = ", DOMAIN_DESCRIPTORS(ID_X)%INDX_MAX
        call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=rc)
        write(tmpstr, *) "TASK_ACTIVE = ", DOMAIN_DESCRIPTORS(ID_X)%TASK_ACTIVE
        call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=rc)
        write(tmpstr, *) "NUM_PETS = ", DOMAIN_DESCRIPTORS(ID_X)%NUM_PETS
        call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=rc)
        write(tmpstr, *) "PET_MAP = ", size(DOMAIN_DESCRIPTORS(ID_X)%PET_MAP), DOMAIN_DESCRIPTORS(ID_X)%PET_MAP
        call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=rc)

        call ESMF_LogFlush()
!
!-----------------------------------------------------------------------
!***  Get the central geographic lat/lon (radians) of the
!***  upper parent's grid.  Also compute the angular increments of
!***  the grid cells in rotated lat/lon.
!-----------------------------------------------------------------------
!
        IF(nmm_int_state%I_AM_A_FCST_TASK(ID_X))THEN
          PHI0=solver_int_state%TPH0D*DEG2RAD
          LAM0=solver_int_state%TLM0D*DEG2RAD
          CENTRAL_LATLON(1)=PHI0
          CENTRAL_LATLON(2)=LAM0
        ENDIF
!
        IF(NUM_DOMAINS_TOTAL>1.AND.ID_X==1)THEN
          CALL ESMF_VMBroadcast(VM, CENTRAL_LATLON, 2, PETLIST_DOMAIN(1,1), rc=RC)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out

          PHI0=CENTRAL_LATLON(1)
          LAM0=CENTRAL_LATLON(2)
        ENDIF
!
!-----------------------------------------------------------------------
!***  Generate the angles to rotate vectors on the NMM-B native grid
!***  to geographic latitude/longitude.  Also compute the grid cell
!***  areas and save the sea masks.  These quantities are only valid
!***  on compute tasks.
!-----------------------------------------------------------------------
!
        compute_tasks: IF(nmm_int_state%I_AM_A_FCST_TASK(ID_X))THEN
!
          ITS=solver_int_state%ITS
          ITE=solver_int_state%ITE
          JTS=solver_int_state%JTS
          JTE=solver_int_state%JTE
!
!-----------------------------------------------------------------------
!***  First allocate the arrays if not already done.
!-----------------------------------------------------------------------
!
          IF(.NOT.ALLOCATED(DOMAIN_DESCRIPTORS(ID_X)%ROT_ANGLE))THEN
!
            ALLOCATE(DOMAIN_DESCRIPTORS(ID_X)%ROT_ANGLE(ITS:ITE,JTS:JTE) &
                    ,stat=RC)
            IF(RC/=0)THEN
              WRITE(0,*)' Failed to allocate rotation angles in'         &
                       ,' STORE_DOMAIN_DESCRIPTORS!'
              WRITE(0,101)ID_X,RC
  101         FORMAT(' Domain #',I2,' RC=',I3)
              WRITE(0,*)' ABORT!'
              CALL ESMF_Finalize(endflag=ESMF_END_ABORT                  &
                                ,rc     =RC )
            ENDIF
!
            ALLOCATE(DOMAIN_DESCRIPTORS(ID_X)%CELL_AREA(ITS:ITE,JTS:JTE) &
                    ,stat=RC)
            IF(RC/=0)THEN
              WRITE(0,*)' Failed to allocate cell areas in'              &
                       ,' STORE_DOMAIN_DESCRIPTORS!'
              WRITE(0,102)ID_X,RC
  102         FORMAT(' Domain #',I2,' RC=',I3)
              WRITE(0,*)' ABORT!'
              CALL ESMF_Finalize(endflag=ESMF_END_ABORT                  &
                                ,rc     =RC )
            ENDIF
!
            ALLOCATE(DOMAIN_DESCRIPTORS(ID_X)%SEA_MASK(ITS:ITE,JTS:JTE)  &
                    ,stat=RC)
            IF(RC/=0)THEN
              WRITE(0,*)' Failed to allocate sea mask in'                &
                       ,' STORE_DOMAIN_DESCRIPTORS!'
              WRITE(0,103)ID_X,RC
  103         FORMAT(' Domain #',I2,' RC=',I3)
              WRITE(0,*)' ABORT!'
              CALL ESMF_Finalize(endflag=ESMF_END_ABORT                  &
                                ,rc     =RC )
            ENDIF
!
          ENDIF
!
!-----------------------------------------------------------------------
!
        ELSE
!
          ALLOCATE(DOMAIN_DESCRIPTORS(ID_X)%ROT_ANGLE(1,1)              &  !<-- Dummy allocation for write tasks
                  ,stat=RC)
          DOMAIN_DESCRIPTORS(ID_X)%ROT_ANGLE(1,1)=-999.
!
          ALLOCATE(DOMAIN_DESCRIPTORS(ID_X)%CELL_AREA(1,1)              &  !<-- Dummy allocation for write tasks
                  ,stat=RC)
          DOMAIN_DESCRIPTORS(ID_X)%CELL_AREA(1,1)=-999.
!
          ALLOCATE(DOMAIN_DESCRIPTORS(ID_X)%SEA_MASK(1,1)               &  !<-- Dummy allocation for write tasks
                  ,stat=RC)
          DOMAIN_DESCRIPTORS(ID_X)%SEA_MASK(1,1)=-999.
!
        ENDIF compute_tasks
!
!-----------------------------------------------------------------------
!
      ENDDO domains
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE STORE_DOMAIN_DESCRIPTORS
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE NMM_INITIALIZE
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE NMM_RUN(NMM_GRID_COMP                                  &
                        ,IMP_STATE                                      &
                        ,EXP_STATE                                      &
                        ,CLOCK_NEMS                                     &
                        ,RC_RUN)
!
!-----------------------------------------------------------------------
!***  This routine executes the integration timeloop for the NMM
!***  through a call to subroutine NMM_INTEGRATE.
!***  That is preceded by digital filtering if it is requested.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(ESMF_GridComp) :: NMM_GRID_COMP                                 !<-- The NMM component
!
      TYPE(ESMF_State) :: IMP_STATE                                     &  !<-- The NMM import state
                         ,EXP_STATE                                        !<-- The NMM export state
!
      TYPE(ESMF_Clock) :: CLOCK_NEMS                                       !<-- The NEMS ESMF Clock
!
      INTEGER,INTENT(OUT) :: RC_RUN                                        !<-- Error return code
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: HDIFF_ON,MYPE_LOCAL                         &
                           ,N,NTIMESTEP                                 &
                           ,YY,MM,DD,H,M,S,Sn,Sd
!
      INTEGER(kind=KINT) :: ID_DOM,Domain_RunstepCount
!
      INTEGER(kind=KINT) :: IERR,RC
!
      INTEGER(kind=ESMF_KIND_I8) :: NTIMESTEP_ESMF
!
      LOGICAL(kind=KLOG) :: ADVANCED,I_AM_ACTIVE,FREE_FORECAST          &
                           ,LAST_GENERATION
!
      TYPE(ESMF_Time) :: CURRTIME
!
      REAL(kind=KDBL) :: gentimer3
!
      REAL(kind=KDBL),DIMENSION(99) :: gentimer1,gentimer2
!
      REAL(ESMF_KIND_R8) :: NTIMESTEP_ESMF_REAL !nuopc
!
      CHARACTER(2) :: INT_TO_CHAR
      CHARACTER(6) :: FMT='(I2.2)'
!
      TYPE(ESMF_TimeInterval)         :: Master_TimeStep
      TYPE(ESMF_TimeInterval),POINTER :: Domain_TimeStep
      TYPE(ESMF_Time)                 :: Master_CurrTime                &
                                        ,Master_StartTime               &
                                        ,Master_StopTime
!
      type(esmf_timeinterval) :: timestep_esmf
      type(esmf_time) :: starttime,stoptime
      INTEGER(kind=KINT) :: iyear_fcst &
                           ,imonth_fcst &
                           ,iday_fcst &
                           ,ihour_fcst &
                           ,iminute_fcst &
                           ,isecond_fcst &
                           ,isecond_num &
                           ,isecond_den
      type(esmf_time) :: ringtime
      type(esmf_timeinterval) :: ringinterval
      integer :: month,day
!!!   integer(esmf_kind_i4) :: yy,mm,dd,h,m,s,sn,sd
      integer(kind=kint) :: coupling_interval_int
      real(kind=kfpt) :: coupling_interval_real
      character(10) :: coupling_interval_char
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      btim0=timef()
!
      gentimer1=0.
      gentimer2=0.
      gentimer3=0.
!
      RC    =ESMF_SUCCESS
      RC_RUN=ESMF_SUCCESS
!
      MESSAGE_CHECK='Print NEMS Clock at start of NMM_RUN'
      if(I_AM_ROOT(RC)) then
        CALL NMM_CLOCKPRINT(CLOCK_NEMS, 'Driver NMM Clock', rc=RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
      endif
!
!-----------------------------------------------------------------------
!***  CLOCK_NEMS is used to control coupling of the NMM with other
!***  models via NUOPC.  Here inside the NMM we need CLOCK_NEMS
!***  for its timestep to compute the number of NMM timesteps per
!***  coupling timestep.
!-----------------------------------------------------------------------
! 
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="NMM_RUN: Get the NEMS Clock"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ClockGet(clock    =CLOCK_NEMS                           &
                        ,currTime =Master_CurrTime                      &  !<-- The beginning time of the current coupling timestep (ESMF)
                        ,StartTime=Master_StartTime                     &
                        ,StopTime =Master_StopTime                      &  !<-- The end time of the current coupling timestep (ESMF)
                        ,timestep =master_timestep                      &
                        ,rc       =RC)
!
!     call esmf_timeintervalget(master_timestep, h=h, m=m, s=s, rc=rc)
!     write(0,*)' NMM_RUN CLOCK_NEMS timestep is h=',h,' m=',m,' s=',s
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      Master_TimeStep = Master_StopTime - Master_CurrTime                  !<-- The coupling timestep length (sec) (ESMF)

!     call esmf_timeintervalget(master_timestep, h=h, m=m, s=s, rc=rc)
!     write(0,*)' master_timestep h=',h,' m=',m,' s=',s
!     call esmf_timeget(master_stoptime, dd=dd, h=h, m=m, s=s, rc=rc)
!     write(0,*)' master_stoptime d=',dd,' h=',h,' m=',m,' s=',s
!     call esmf_timeget(master_currtime, dd=dd, h=h, m=m, s=s, rc=rc)
!     write(0,*)' master_currtime d=',dd,' h=',h,' m=',m,' s=',s
!
      DO N = 1, NUM_DOMAINS_TOTAL
        ID_DOM=RANK_TO_DOMAIN_ID_ptr(N)
        Domain_TimeStep=>nmm_int_state%TIMESTEP(ID_DOM)                    !<-- The timestep length (DT) of this Domain (sec) (ESMF)
        Domain_RunstepCount = nint(Master_TimeStep/Domain_TimeStep)        !<-- # of Domain timesteps per coupling timestep

!        if(I_AM_ROOT(RC)) print *, 'Domain_RunstepCount: ', Domain_RunstepCount
!       call esmf_timeintervalget(domain_timestep, h=h, m=m, s=s, rc=rc)
!       write(0,*)' domain_timestep h=',h,' m=',m,' s=',s
!       write(0,*)' domain_runstepcount=',domain_runstepcount
!
      END DO
!
!-----------------------------------------------------------------------
!***  Extract the digital filter specification from the configure file.
!***  If it is >0 then the user is asking that one of the filters be
!***  used prior to the free forecast.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Digital Filter: Extract Filter Method"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF(1)                         &  !<-- Use uppermost parent domain; all domains the same
                                  ,value =FILTER_METHOD                 &  !<-- The digital filter flag
                                  ,label ='filter_method:'              &  !<-- Give this label's value to preceding variable
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  If the user has requested digital filtering then proceed with the
!***  selected method before performing the normal forecast integration.
!-----------------------------------------------------------------------
!
      IF(FILTER_METHOD>0)THEN                                              !<-- If true then filtering was selected
!
!-----------------------------------------------------------------------
!
        DO N=1,NUM_GENS
!
          MY_DOMAIN_ID=MY_DOMAINS_IN_GENS(N)                                 !<-- Task's domain in generation N

          IF(MY_DOMAIN_ID>0)THEN                                             !<-- Domain ID is 0 for 2-way nesting if task not in generation N
!
            IMP_STATE_DOMAIN=>nmm_int_state%IMP_STATE_DOMAIN(MY_DOMAIN_ID)   !<-- This domain's import state
            IMP_STATE_CPL_NEST=>nmm_int_state%IMP_STATE_PC_CPL(MY_DOMAIN_ID)
            FREE_FORECAST=.FALSE.
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="NMM_Run: Set Free Forecast flag in the Domain import state"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_AttributeSet(state=IMP_STATE_DOMAIN             &  !<-- This DOMAIN component's import state
                                  ,name ='Free Forecast'              &  !<-- The forecast is in the digital filter.
                                  ,value=FREE_FORECAST                &  !<-- Value of filter method flag
                                  ,rc   =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="NMM_Run: Set Free Forecast flag in the P-C import state"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_AttributeSet(state=IMP_STATE_CPL_NEST           &  !<-- This DOMAIN component's import state
                                  ,name ='Free Forecast'              &  !<-- The forecast is in the digital filter.
                                  ,value=FREE_FORECAST                &  !<-- Value of filter method flag
                                  ,rc   =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          ENDIF
!
        ENDDO
!
!-----------------------------------------------------------------------
!
        CALL RUN_DIGITAL_FILTER_NMM                                        !<-- See internal subroutine below.
!
!-----------------------------------------------------------------------
!
      ENDIF
!
!-----------------------------------------------------------------------
!***  If there was digital filtering it is finished so set the filter
!***  method to 0 for the free forecast on all this task's domains.
!-----------------------------------------------------------------------
!
      FILTER_METHOD=0                                                      !<-- Filter is done or was not run so set method to 0
      I_AM_ACTIVE=.TRUE.                                                   !<-- All domains are active in the free forecast.
!
      DO N=1,NUM_GENS
!
        MY_DOMAIN_ID=MY_DOMAINS_IN_GENS(N)                                 !<-- Task's domain in generation N

        IF(MY_DOMAIN_ID>0)THEN                                             !<-- Domain ID is 0 for 2-way nesting if task not in generation N
!
          IMP_STATE_DOMAIN=>nmm_int_state%IMP_STATE_DOMAIN(MY_DOMAIN_ID)   !<-- This domain's import state
          EXP_STATE_DOMAIN=>nmm_int_state%EXP_STATE_DOMAIN(MY_DOMAIN_ID)   !<-- This domain's export state
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="NMM_Run: Set Filter Method to 0 in DOMAIN import state"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeSet(state=IMP_STATE_DOMAIN             &  !<-- This DOMAIN component's import state 
                                ,name ='Filter_Method'              &  !<-- Flag for filter method
                                ,value=FILTER_METHOD                &  !<-- Value of filter method flag
                                ,rc   =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="NMM_Run: Set domain active flag in DOMAIN export state"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeSet(state=EXP_STATE_DOMAIN             &  !<-- This DOMAIN component's export state
                                ,name ='I Am Active'                &  !<-- This domain is active in the forecast.
                                ,value=I_AM_ACTIVE                  &  !<-- Value of filter method flag
                                ,rc   =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          IMP_STATE_DOMAIN=>nmm_int_state%IMP_STATE_DOMAIN(MY_DOMAIN_ID)   !<-- This domain's import state
          IMP_STATE_CPL_NEST=>nmm_int_state%IMP_STATE_PC_CPL(MY_DOMAIN_ID)
          FREE_FORECAST=.TRUE.
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="NMM_Run: Set Free Forecast flag in the Domain import state"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeSet(state=IMP_STATE_DOMAIN                 &  !<-- This DOMAIN component's import state
                                ,name ='Free Forecast'                  &  !<-- The forecast is now free.
                                ,value=FREE_FORECAST                    &  !<-- Is this the free forecast?
                                ,rc   =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="NMM_Run: Set Free Forecast flag in the P-C import state"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeSet(state=IMP_STATE_CPL_NEST               &  !<-- This P-C coupler component's import state
                                ,name ='Free Forecast'                  &  !<-- The forecast is now free.
                                ,value=FREE_FORECAST                    &  !<-- Is this the free forecast?
                                ,rc   =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        ENDIF
!
      ENDDO
!
!-----------------------------------------------------------------------
!***  Prepare to run the free forecast.
!-----------------------------------------------------------------------
!
      ALL_FORECASTS_COMPLETE=.FALSE.
!
      DO N=1,NUM_GENS
        IF(MY_DOMAINS_IN_GENS(N)>0)THEN
          GENERATION_FINISHED(N)=.FALSE.
        ELSE
          GENERATION_FINISHED(N)=.TRUE.                                !<-- Task not in this generation; consider it finished.
        ENDIF
      ENDDO
!
!-----------------------------------------------------------------------
!
      main_block: DO WHILE(.NOT.ALL_FORECASTS_COMPLETE)
!
!-----------------------------------------------------------------------
!***  The execution of the timestepping differs fundamentally between
!***  1-way and 2-way nesting.  In 1-way nesting each task belongs to
!***  only one domain and all domains run concurrently from the start
!***  to the end of the forecast.  In 2-way nesting some or all tasks
!***  will lie in more than one domain but never more than one domain
!***  per generation therefore a loop over the generations must exist
!***  above partial timestep loops allowing tasks to return after the
!***  timestep is finished so they can participate in a different
!***  generation's timestep(s) before switching generations again.
!***  Thus NUM_GENS in generations_loop is a relevant integer >1 only 
!***  for 2-way nesting.
!-----------------------------------------------------------------------
!
        btim0=timef()
        generations_loop: DO N=1,NUM_GENS                                  !<-- A single iteration for 1-way nesting
!
!-----------------------------------------------------------------------
!
          IF(GENERATION_FINISHED(N))THEN
            CYCLE generations_loop
          ENDIF
!
          LAST_GENERATION=.FALSE.
          IF(N==NUM_GENS)LAST_GENERATION=.TRUE.
          MY_DOMAIN_ID=MY_DOMAINS_IN_GENS(N)                               !<-- Multiple generations only apply to 2-way nesting
!
!-----------------------------------------------------------------------
!
          domain: IF(MY_DOMAIN_ID>0)THEN                                   !<-- Domain ID is 0 for 2-way nesting if task not in generation N
!
!-----------------------------------------------------------------------
!
            btim=timef()
            DOMAIN_GRID_COMP=>nmm_int_state%DOMAIN_GRID_COMP(MY_DOMAIN_ID) !<-- This domain's ESMF component
            IMP_STATE_DOMAIN=>nmm_int_state%IMP_STATE_DOMAIN(MY_DOMAIN_ID) !<-- Its import state
            EXP_STATE_DOMAIN=>nmm_int_state%EXP_STATE_DOMAIN(MY_DOMAIN_ID) !<-- Its export state
!
            I_AM_A_FCST_TASK   =>nmm_int_state%I_AM_A_FCST_TASK(MY_DOMAIN_ID)    !<-- Is this task a fcst task on this domain?
            I_AM_LEAD_FCST_TASK=>nmm_int_state%I_AM_LEAD_FCST_TASK(MY_DOMAIN_ID) !<-- Is this the lead fcst task on this domain?
!
!-----------------------------------------------------------------------
!***  Again obtain current information from the Clock.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="NMM_Run: Get current time info from the Clock"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_ClockGet(clock       =CLOCK_NMM(MY_DOMAIN_ID)     &
                              ,starttime   =STARTTIME                   &
                              ,currtime    =CURRTIME                    &
                              ,advanceCount=NTIMESTEP_ESMF              &
                              ,runduration =RUNDURATION                 &
                              ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  We need the local MPI task ID on the given NMM domain.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="NMM_RUN: Retrieve VM from DOMAIN Gridded Component"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_GridCompGet(gridcomp=DOMAIN_GRID_COMP             &  !<-- The DOMAIN gridded component
                                 ,vm      =VM                           &  !<-- Get the Virtual Machine from the DOMAIN component
                                 ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="NMM_RUN: Obtain the Local Task ID"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_VMGet(vm      =VM                                 &  !<-- The virtual machine for this DOMAIN component
                           ,localpet=MYPE_LOCAL                         &  !<-- Each task's local rank on this domain
                           ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Set default value for horizontal diffusion flag (1-->ON).
!-----------------------------------------------------------------------
!
            HDIFF_ON=1
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="NMM_Run: Put Horizontal Diffusion Flag into DOMAIN import state"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_AttributeSet(state=IMP_STATE_DOMAIN               &  !<-- This DOMAIN component's import state
                                  ,name ='HDIFF'                        &  !<-- Flag for diffusion on/off
                                  ,value=HDIFF_ON                       &  !<-- Value of horizontal diffusion flag
                                  ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Execute the normal forecast integration after dereferencing
!***  argument variables for this particular domain.
!-----------------------------------------------------------------------
!
            COMM_TO_MY_PARENT=>comms_domain(MY_DOMAIN_ID)%TO_PARENT        !<-- This domain's intracommunicator to its parent
            NUM_CHILDREN     =>nmm_int_state%NUM_CHILDREN(MY_DOMAIN_ID)    !<-- How many children does this domain have?
!
            NTIMESTEP=NTIMESTEP_ESMF                                       !<-- This domains' current timestep count (integer)
            TIMESTEP=>nmm_int_state%TIMESTEP(MY_DOMAIN_ID)                 !<-- This domain's fundamental timestep (sec) (ESMF)
!
            CALL ESMF_TimeGet(CURRTIME, dd=DD, h=H, m=M, s=S, rc=RC)
!           IF (I_AM_LEAD_FCST_TASK) WRITE(0,*) 'CURRTIME going into normal NMM_INTEG: ', DD, H, M, S
!
            CALL ESMF_TimeGet(STARTTIME, dd=DD, h=H, m=M, s=S, rc=RC)
!           IF (I_AM_LEAD_FCST_TASK) WRITE(0,*) 'STARTTIME going into normal NMM_INTEG: ', DD, H, M, S
!
            INTERVAL_CLOCKTIME=>nmm_int_state%INTERVAL_CLOCKTIME(MY_DOMAIN_ID) !<-- Time interval for this domain's clocktime prints
            INTERVAL_HISTORY  =>nmm_int_state%INTERVAL_HISTORY(MY_DOMAIN_ID)   !<-- Time interval for this domain's history output
            INTERVAL_RESTART  =>nmm_int_state%INTERVAL_RESTART(MY_DOMAIN_ID)   !<-- Time interval for this domain's restart output
!
            NPE_PRINT=>nmm_int_state%NPE_PRINT(MY_DOMAIN_ID)               !<-- Print clocktimes from this task
!
            RESTARTED_RUN=>nmm_int_state%RESTARTED_RUN(MY_DOMAIN_ID)         
            RST_OUT_00   =>nmm_int_state%RST_OUT_00(MY_DOMAIN_ID)
!
            I_AM_A_NEST     =>nmm_int_state%I_AM_A_NEST(MY_DOMAIN_ID)      !<-- Is this domain a nest?
!
            PARENT_CHILD_COUPLER_COMP=>nmm_int_state%PC_CPL_COMP(MY_DOMAIN_ID) !<-- The P-C coupler associated with this domain
            IMP_STATE_CPL_NEST=>nmm_int_state%IMP_STATE_PC_CPL(MY_DOMAIN_ID)   !<-- The P-C coupler import state
            EXP_STATE_CPL_NEST=>nmm_int_state%EXP_STATE_PC_CPL(MY_DOMAIN_ID)   !<-- The P-C coupler export state
!
            PARENT_CHILD_TIME_RATIO=>nmm_int_state%P_C_TIME_RATIO(MY_DOMAIN_ID) !<-- Ratio of this domain's timestep to its parent's
            MY_DOMAIN_MOVES=nmm_int_state%MY_DOMAIN_MOVES(MY_DOMAIN_ID)         !<-- Does this domain move?
            NEST_MODE=nmm_int_state%NEST_MODE(MY_DOMAIN_ID)                     !<-- Is this domain involved in any 2-way nesting?
            NUM_2WAY_CHILDREN=>nmm_int_state%NUM_2WAY_CHILDREN(MY_DOMAIN_ID)    !<-- How many 2-way children on this domain?
            NTRACK=nmm_int_state%NTRACK(MY_DOMAIN_ID)                           !<-- Storm locator flag
            NPHS=nmm_int_state%NPHS(MY_DOMAIN_ID)                               !<-- Physics timestep
            ADVANCED=.FALSE.                                                    !<-- Does the integration advance?
!
            gentimer1(my_domain_id)=gentimer1(my_domain_id)+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Call the timestepping routine.
!-----------------------------------------------------------------------
!     if(mype_local==0)then
!       write(0,77501)n,my_domain_id,mype_local,ntimestep
77501   format(' NMM_RUN before NMM_INTEGRATE generation #',i2,' my_domain_id=',i2,' mype_local=',i3,' ntimestep=',i5)
!     endif
!
            btim=timef()
            CALL NMM_INTEGRATE(clock_direction    ='Forward '                &
                              ,domain_grid_comp   =DOMAIN_GRID_COMP          &
                              ,imp_state_domain   =IMP_STATE_DOMAIN          &
                              ,exp_state_domain   =EXP_STATE_DOMAIN          &
                              ,clock_integrate    =CLOCK_NMM(MY_DOMAIN_ID)   &
                              ,currtime           =CURRTIME                  &
                              ,starttime          =STARTTIME                 &
                              ,timestep           =TIMESTEP                  &
                              ,ntimestep_ext      =NTIMESTEP                 &
                              ,runstepcount       =Domain_RunstepCount       &
                              ,dt                 =DT(MY_DOMAIN_ID)          &
                              ,interval_clocktime =INTERVAL_CLOCKTIME        &
                              ,interval_history   =INTERVAL_HISTORY          &
                              ,interval_restart   =INTERVAL_RESTART          &
                              ,filter_method      =FILTER_METHOD             &
                              ,restarted_run      =RESTARTED_RUN             &
                              ,rst_out_00         =RST_OUT_00                &
                              ,i_am_a_fcst_task   =I_AM_A_FCST_TASK          &
                              ,i_am_lead_fcst_task=I_AM_LEAD_FCST_TASK       &
                              ,nesting            =NESTING_NMM               &
                              ,nest_mode          =NEST_MODE                 &
                              ,task_mode          =TASK_MODE                 &
                              ,i_am_a_nest        =I_AM_A_NEST               &
                              ,my_domain_id       =MY_DOMAIN_ID              &
                              ,num_children       =NUM_CHILDREN              &
                              ,num_2way_children  =NUM_2WAY_CHILDREN         &
                              ,parent_child_cpl   =PARENT_CHILD_COUPLER_COMP &
                              ,imp_state_cpl_nest =IMP_STATE_CPL_NEST        &
                              ,exp_state_cpl_nest =EXP_STATE_CPL_NEST        &
                              ,par_chi_time_ratio =PARENT_CHILD_TIME_RATIO   &
                              ,my_domain_moves    =MY_DOMAIN_MOVES           &
                              ,ntrack             =NTRACK                    &
                              ,nphs               =NPHS                      &
                              ,last_generation    =LAST_GENERATION           &
                              ,advanced           =ADVANCED                  &
                              ,mype               =MYPE_LOCAL                &
                              ,comm_global        =COMM_GLOBAL               &
                              ,timers_domain      =TIMERS(MY_DOMAIN_ID)      &
                              ,npe_print          =NPE_PRINT                 &
                              ,print_timing       =PRINT_TIMING )
!     write(0,40403)n,my_domain_id,mype_local,ntimestep,advanced
40403 format(' NMM_RUN after NMM_INTEGRATE generation #',i2,' my_domain_id=',i2,' mype_local=',i4,' ntimestep=',i5,' advanced=',l1)
!     if(mype_local==0)then
!     call esmf_clockget(clock=clock_nmm(my_domain_id) &
!                       ,timeStep=timestep_esmf &
!                       ,startTime=starttime &
!                       ,currTime=currtime &
!                       ,stopTime=stoptime &
!                       ,rc=rc)
!     call esmf_timeintervalget(timeinterval=timestep_esmf &
!                      ,s=isecond_fcst  &
!                      ,sn=isecond_num  &
!                      ,sd=isecond_den &
!                      ,rc=rc)
!     write(0,55020)isecond_fcst,isecond_num,isecond_den
55020 format(' timestep: sec=',i4,' sn=',i3,' sd=',i3)
!     call esmf_timeget(time=starttime &
!                        ,yy  =iyear_fcst &
!                        ,mm  =imonth_fcst &
!                        ,dd  =iday_fcst  &
!                        ,h   =ihour_fcst &
!                        ,m   =iminute_fcst  &
!                        ,s   =isecond_fcst  &
!                        ,sN  =isecond_num  &
!                        ,sD  =isecond_den &
!                      ,rc=rc)
!     write(0,55021)iyear_fcst,imonth_fcst,iday_fcst &
!                  ,ihour_fcst,iminute_fcst &
!                  ,isecond_fcst,isecond_num,isecond_den
55021 format(' starttime: year=',i4,' month=',i2,' day=',i2,' hour=',i2 &
            ,' min=',i2,' sec=',i2,' sn=',i3,' sd=',i3)
!     call esmf_timeget(time=currtime &
!                        ,yy  =iyear_fcst &
!                        ,mm  =imonth_fcst &
!                        ,dd  =iday_fcst  &
!                        ,h   =ihour_fcst &
!                        ,m   =iminute_fcst  &
!                        ,s   =isecond_fcst  &
!                        ,sN  =isecond_num  &
!                        ,sD  =isecond_den &
!                      ,rc=rc)
!     write(0,55022)iyear_fcst,imonth_fcst,iday_fcst &
!                  ,ihour_fcst,iminute_fcst &
!                  ,isecond_fcst,isecond_num,isecond_den
55022 format(' currtime: year=',i4,' month=',i2,' day=',i2,' hour=',i2 &
            ,' min=',i2,' sec=',i2,' sn=',i3,' sd=',i3)
!     call esmf_timeget(time=stoptime &
!                        ,yy  =iyear_fcst &
!                        ,mm  =imonth_fcst &
!                        ,dd  =iday_fcst  &
!                        ,h   =ihour_fcst &
!                        ,m   =iminute_fcst  &
!                        ,s   =isecond_fcst  &
!                        ,sN  =isecond_num  &
!                        ,sD  =isecond_den &
!                      ,rc=rc)
!     write(0,55023)iyear_fcst,imonth_fcst,iday_fcst &
!                  ,ihour_fcst,iminute_fcst &
!                  ,isecond_fcst,isecond_num,isecond_den
55023 format(' stoptime: year=',i4,' month=',i2,' day=',i2,' hour=',i2 &
            ,' min=',i2,' sec=',i2,' sn=',i3,' sd=',i3)
!     endif
!
            gentimer2(my_domain_id)=gentimer2(my_domain_id)+(timef()-btim)
!
!-----------------------------------------------------------------------
!
           IF(ESMF_AlarmIsRinging(alarm=ALARM_CPL(MY_DOMAIN_ID)         &  !<-- Is it time to couple atmosphere and ocean?
                                 ,rc   =RC)                             &
              .AND. ADVANCED ) THEN                                        !<-- Did the integration actually advance?
             GENERATION_FINISHED(N)=.TRUE.                                 !<-- Task's fcst in generation N has finished this cpling interval
!       write(0,23230)n
23230   format(' generation_finished true after NMM_INTEGRATE generation=',i2)
!     else
!       if(.not.generation_finished(n))then
!         write(0,23231)n
23231     format(' generation_finished false after NMM_INTEGRATE generation=',i2)
!       endif
           ENDIF
!     call esmf_alarmget(alarm=ALARM_CPL(MY_DOMAIN_ID) &
!                       ,ringtime=ringtime &
!                       ,ringinterval=ringinterval &
!                       ,rc=rc)
!     call esmf_timeget(time=ringtime &
!                      ,yy=yy &
!                      ,mm=month &
!                      ,dd=dd &
!                      ,h=h &
!                      ,m=m &
!                      ,s=s &
!                      ,sn=sn &
!                      ,sd=sd)
!     write(0,*)' ringtime: y=',yy,' mm=',month,' dd=',dd,' h=',h &
!              ,' m=',m,' s=',s,' sn=',sn,'s d=',sd
!     call esmf_timeintervalget(timeinterval=ringinterval &
!                      ,m=m &
!                      ,s=s &
!                      ,sn=sn &
!                      ,sd=sd)
!     write(0,*)' ringinterval: m=',m,' s=',s,' sn=',sn,'s d=',sd
!
!     write(0,55024)generation_finished(n),n
55024 format(' NMM_RUN end of domain loop generation_finished=',l1,' generation #',i2)
!-----------------------------------------------------------------------
!
          ENDIF domain
!
!     write(0,55025)generation_finished(n),n,num_gens
55025 format(' NMM_RUN after domain block generation_finished=',l1,' generation #',i2,' num_gens=',i2)
!-----------------------------------------------------------------------
!***  All tasks that are finished on all generations may leave.
!-----------------------------------------------------------------------
!
          IF(ALL(GENERATION_FINISHED,NUM_GENS))THEN                        !<-- If true, all of this task's domains are finished
            ALL_FORECASTS_COMPLETE=.TRUE.
!     write(0,55125)n,generation_finished
55125 format(' all_forecasts_complete is true  generation #',i2,' generation_finished=',3(1x,l1))
            EXIT generations_loop
!  else
!    write(0,55026)n,generation_finished
55026 format(' all_forecasts_complete is false  generation #',i2,' generation_finished=',3(1x,l1))
          ENDIF
!
!-----------------------------------------------------------------------
!
        ENDDO generations_loop
!
        gentimer3=gentimer3+timef()-btim0
!
!-----------------------------------------------------------------------
!
      ENDDO main_block
!
!-----------------------------------------------------------------------
!***  Currently in coupled runs only the upper parent is coupled to
!***  the ocean.  The nests receive SST updates from interpolation
!***  by the parent at the start of the upcoming coupling interval
!***  so now generate regrid interpolation weights between parent
!***  and nests for the nests' current positions.  The weights are
!***  stored in DOMAIN_DESCRIPTORS.
!-----------------------------------------------------------------------
!
      DO N=2,NUM_DOMAINS_TOTAL
!
        MY_DOMAIN_ID=RANK_TO_DOMAIN_ID_PTR(N)                              !<-- The domain ID for the Nth domain
!
        CALL NMMB_CreateRouteHandle(MY_DOMAIN_ID                        &
                                   ,DOMAIN_DESCRIPTORS                  &
                                   ,RC=RC)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,  &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
!
      ENDDO
!
!-----------------------------------------------------------------------
!***  Print clocktimes for:
!***   (1) Parent tasks learning if they can proceed with integration
!***       based on information from their parent and their children.
!***   (2) Child tasks sharing information as to whether their parent
!***       has sent an all clear signal that the parent has received
!***       2-way exchange data from all its children thus the children
!***       may proceed with their integration.
!-----------------------------------------------------------------------
!
      IF (NUM_GENS > 1) THEN
        DO N=1,NUM_GENS
!
           MY_DOMAIN_ID = MY_DOMAINS_IN_GENS(N)
           IF(MY_DOMAIN_ID>0) THEN
             IF (MY_DOMAIN_ID == 1) THEN
              WRITE(0,896)my_domain_id,            &
                          gentimer1(my_domain_id), &
                          gentimer2(my_domain_id), &
                          gentimer3,               & 
                          pbcst_tim(my_domain_id)
            ELSE
              WRITE(0,897)my_domain_id,            &
                          gentimer1(my_domain_id), &
                          gentimer2(my_domain_id), &
                          gentimer3,               &
                          cbcst_tim(my_domain_id)
            ENDIF
          ENDIF
!
        ENDDO
      ENDIF
!
  896 format (' For domain ',i2,' t1,t2,t3,pb ',4(g10.3))
  897 format (' For domain ',i2,' t1,t2,t3,cb ',4(g10.3))
!
!-----------------------------------------------------------------------
!
      MESSAGE_CHECK='Print domain 1 Clock and NEMS Clock'
      if(I_AM_ROOT(RC)) then
        CALL NMM_CLOCKPRINT(CLOCK_NMM(1), 'Native NMM Clock', rc=RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
        CALL NMM_CLOCKPRINT(CLOCK_NEMS, 'Driver NMM Clock', rc=RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
      endif
!
      ALL_FORECASTS_COMPLETE = .FALSE. ! Force the next time step to run
      DO N = 1, NUM_DOMAINS_TOTAL
        GENERATION_FINISHED(N) = .FALSE.
      ENDDO
!
!-----------------------------------------------------------------------
!
      IF(RC_RUN==ESMF_SUCCESS)THEN
!       WRITE(0,*)' NMM_RUN succeeded'
      ELSE
        WRITE(0,*)' NMM_RUN failed  RC_RUN=',RC_RUN
      ENDIF
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
      SUBROUTINE RUN_DIGITAL_FILTER_NMM
!
!-----------------------------------------------------------------------
!***  This routine executes the digital filters for the NMM
!***  if specified by the user.
!-----------------------------------------------------------------------
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: HDIFF_ON,MEAN_ON                            &
                           ,N,NTIMESTEP                                 &
                           ,YY,MM,DD,H,M,S
!
      INTEGER(kind=KINT) :: RC
!
      INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE,SAVE :: NDFISTEP
!
      INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE :: LOC_PAR_CHILD_TIME_RATIO
!
      TYPE(ESMF_Clock),DIMENSION(:),ALLOCATABLE :: CLOCK_FILTER
!
      TYPE(ESMF_Time) :: SDFITIME
!
      TYPE(ESMF_Time),DIMENSION(:),ALLOCATABLE,SAVE :: HALFDFITIME
!
      TYPE(ESMF_TimeInterval) :: TIMESTEP_FILTER
!
      TYPE(ESMF_TimeInterval),DIMENSION(:),ALLOCATABLE,SAVE :: HALFDFIINTVAL
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Allocate the clocks to control the filtering.
!-----------------------------------------------------------------------
!
      ALLOCATE(CLOCK_FILTER(1:NUM_DOMAINS_TOTAL),stat=RC)
!
      IF(RC/=0)THEN
        WRITE(0,*)' Error allocating filter clocks; rc=',RC
      ENDIF
!
!-----------------------------------------------------------------------
!
      ALLOCATE(HALFDFITIME(1:NUM_DOMAINS_TOTAL),stat=RC)
      ALLOCATE(HALFDFIINTVAL(1:NUM_DOMAINS_TOTAL),stat=RC)
!
      ALLOCATE(NDFISTEP(1:NUM_DOMAINS_TOTAL),stat=RC)
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      method_block: IF(FILTER_METHOD==1)THEN                               !<-- The DFL digital filter.
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  First loop through the generations to set preliminary variables
!***  specific to the domains in each generation.
!-----------------------------------------------------------------------
!
        gens_f1_1: DO N=1,NUM_GENS   
!
          MY_DOMAIN_ID=MY_DOMAINS_IN_GENS(N)                               !<-- This task's (only) domain in generation N
          IF(MY_DOMAIN_ID==0)CYCLE                                         !<-- This task is not on a domain in generation N
!
!-----------------------------------------------------------------------
! 
          IMP_STATE_DOMAIN=>nmm_int_state%IMP_STATE_DOMAIN(MY_DOMAIN_ID)
!
          FILT_TIMESTEP=>nmm_int_state%FILT_TIMESTEP(MY_DOMAIN_ID)         !<-- This domain's timestep for digital filtering
!
          I_AM_LEAD_FCST_TASK=>nmm_int_state%I_AM_LEAD_FCST_TASK(MY_DOMAIN_ID) !<-- Is this the lead fcst task on this domain?
!
          IF(I_AM_LEAD_FCST_TASK)WRITE(0,*)' Beginning DFL Filter'
!
!-----------------------------------------------------------------------
!***  Extract the length of the half forward filter window
!***  for the DFL filter.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Digital Filter: Extract DFIHR Value for DFL"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)          &  !<-- This domain's config object
                                      ,value =DFIHR                     &  !<-- Half foward filter window (s) 
                                      ,label ='nsecs_dfl:'              &  !<-- Give this label's value to preceding variable
                                      ,rc    =RC)
!
          CALL ESMF_TimeIntervalSet(timeinterval=HALFDFIINTVAL(MY_DOMAIN_ID) &
                                   ,s           =DFIHR                       &
                                   ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Is DFIHR divided evenly by the timestep?  We cannot proceed
!***  unless it is.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="For DFL Get Actual Timestep from ESMF Value"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_TimeIntervalGet(timeinterval=FILT_TIMESTEP          &  !<-- The filter timestep of this domain (sec) (ESMF)
                                   ,s           =S                      &  !<-- Integer part of timestep
                                   ,sn          =Sn                     &  !<-- Numerator of fractional part
                                   ,sd          =Sd                     &  !<-- Denominator of fractional part
                                   ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          NDFISTEP(MY_DOMAIN_ID)=INT( 0.1+DFIHR/(S+REAL(Sn)/REAL(Sd)))
          DFIHR_CHK=INT(0.1+NDFISTEP(MY_DOMAIN_ID)*(S+REAL(Sn)/REAL(Sd)))
!
          IF (DFIHR /= DFIHR_CHK) THEN
            WRITE(0,*)' DFIHR=',DFIHR,' DFIHR_CHK=',DFIHR_CHK,' for domain #',my_domain_id
            WRITE(0,*)' nsecs_dfl in configure MUST be integer multiple of the timestep'
            WRITE(0,*)' User must reset the value'
            WRITE(0,*)' ABORTING!!'
            CALL ESMF_Finalize(rc=RC,endflag=ESMF_END_ABORT)
          ENDIF
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="DFL Filter: Get current time info from NMM Clock"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_ClockGet(clock   =CLOCK_NMM(MY_DOMAIN_ID)           &
                            ,currtime=CURRTIME                          &
                            ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          STARTTIME=CURRTIME
          HALFDFITIME(MY_DOMAIN_ID)=STARTTIME+HALFDFIINTVAL(MY_DOMAIN_ID)
          SDFITIME=STARTTIME
          DFITIME=HALFDFITIME(MY_DOMAIN_ID)+HALFDFIINTVAL(MY_DOMAIN_ID)
!
          TIMESTEP_FILTER=FILT_TIMESTEP
!
!-----------------------------------------------------------------------
!***  In preparation for this filter's forward integration
!***  create a clock to control the filter's timestepping.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Create the Clock for the DFL Digital Filter."
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CLOCK_FILTER(MY_DOMAIN_ID)=ESMF_ClockCreate(name     ='CLOCK_DFL'     &  !<-- The clock for the DFI filter
                                                     ,timeStep =TIMESTEP_FILTER &  !<-- The filter timestep in this domain
                                                     ,startTime=STARTTIME       &  !<-- Start time of filter
                                                     ,stopTime =DFITIME         &  !<-- Stop time of the filter
                                                     ,rc       =RC)
!
	  CALL ESMF_ClockSet(clock    =CLOCK_FILTER(MY_DOMAIN_ID)       &
	                    ,currtime =CURRTIME                         &
	                    ,starttime=CURRTIME                         &
	                    ,rc       =RC)
!
          CALL ESMF_TimeGet(CURRTIME, dd=DD, h=H, m=M, s=S, rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          HDIFF_ON=1                                                       !<-- Forward integration so we want horiz diffusion.
          MEAN_ON =0
!
          CALL ESMF_AttributeSet(state=IMP_STATE_DOMAIN                 &  !<-- This DOMAIN component's import state for filter
                                ,name ='Clock_Direction'                &
                                ,value='Forward '                       &
                                ,rc   =RC)
!
          CALL ESMF_AttributeSet(state=IMP_STATE_DOMAIN                 &  !<-- This DOMAIN component's import state for filter
                                ,name ='HDIFF'                          &  !<-- Flag for horizontal diffusion on/off
                                ,value=HDIFF_ON                         &  !<-- Value of horizontal diffusion flag
                                ,rc   =RC)
!
          CALL ESMF_AttributeSet(state=IMP_STATE_DOMAIN                 &  !<-- This DOMAIN component's import state for filter
                                ,name ='MEAN_ON'                        &
                                ,value=MEAN_ON                          &
                                ,rc   =RC)
!
          CALL ESMF_AttributeSet(state=IMP_STATE_DOMAIN                 &  !<-- This DOMAIN component's import state for filter
                                ,name ='NDFISTEP'                       &
                                ,value=NDFISTEP(MY_DOMAIN_ID)           &
                                ,rc   =RC)
!
!-----------------------------------------------------------------------
!
        ENDDO gens_f1_1
!
!-----------------------------------------------------------------------
!***  Execute the DFL filter's integration for all domains after
!***  dereferencing argument variables for the given domain.
!***  See fuller explanation in subroutine NMM_RUN.
!-----------------------------------------------------------------------
!
        dfl_int: DO WHILE(.NOT.ALL_FORECASTS_COMPLETE)
!
!-----------------------------------------------------------------------
!
          gens_f1_2: DO N=1,NUM_GENS   
!
!-----------------------------------------------------------------------
!
            IF(GENERATION_FINISHED(N))THEN
              CYCLE gens_f1_2
            ENDIF
!
            LAST_GENERATION=.FALSE.
            IF(N==NUM_GENS)LAST_GENERATION=.TRUE.
!
            MY_DOMAIN_ID=MY_DOMAINS_IN_GENS(N)                             !<-- This task's (only) domain in generation N
            IF(MY_DOMAIN_ID==0)CYCLE                                       !<-- This task is not on a domain in generation N
!
            DOMAIN_GRID_COMP=>nmm_int_state%DOMAIN_GRID_COMP(MY_DOMAIN_ID) !<-- This domain's ESMF component
!
!-----------------------------------------------------------------------
!***  We need the task's rank on the current domain.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="NMM_RUN: Retrieve VM from DOMAIN Gridded Component"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_GridCompGet(gridcomp=DOMAIN_GRID_COMP             &  !<-- The DOMAIN gridded component
                                 ,vm      =VM                           &  !<-- Get the Virtual Machine from the DOMAIN component
                                 ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="NMM_RUN: Obtain the Local Task ID"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_VMGet(vm      =VM                                 &  !<-- The virtual machine for this DOMAIN component
                           ,localpet=MYPE_LOCAL                         &  !<-- Each task's local rank on this domain
                           ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
            IMP_STATE_DOMAIN=>nmm_int_state%IMP_STATE_DOMAIN(MY_DOMAIN_ID) !<-- Its import state
            EXP_STATE_DOMAIN=>nmm_int_state%EXP_STATE_DOMAIN(MY_DOMAIN_ID) !<-- Its export state
!
            I_AM_A_FCST_TASK   =>nmm_int_state%I_AM_A_FCST_TASK(MY_DOMAIN_ID)    !<-- Is this task a fcst task on this domain?
            I_AM_LEAD_FCST_TASK=>nmm_int_state%I_AM_LEAD_FCST_TASK(MY_DOMAIN_ID) !<-- Is this the lead fcst task on this domain?
!
            FILT_TIMESTEP=>nmm_int_state%FILT_TIMESTEP(MY_DOMAIN_ID)       !<-- This domain's timestep for digital filtering
!
            COMM_TO_MY_PARENT=>comms_domain(MY_DOMAIN_ID)%TO_PARENT        !<-- This domain's intracommunicator to its parent
            NUM_CHILDREN=>nmm_int_state%NUM_CHILDREN(MY_DOMAIN_ID)         !<-- How many children does this domain have?
!
            NPE_PRINT=>nmm_int_state%NPE_PRINT(MY_DOMAIN_ID)               !<-- Print clocktimes from this task
!
            RESTARTED_RUN=>nmm_int_state%RESTARTED_RUN(MY_DOMAIN_ID)
            RST_OUT_00=>nmm_int_state%RST_OUT_00(MY_DOMAIN_ID)
!
            I_AM_A_NEST=>nmm_int_state%I_AM_A_NEST(MY_DOMAIN_ID)           !<-- Is this domain a nest?
!
            PARENT_CHILD_COUPLER_COMP=>nmm_int_state%PC_CPL_COMP(MY_DOMAIN_ID) !<-- The P-C coupler associated with this domain
            IMP_STATE_CPL_NEST=>nmm_int_state%IMP_STATE_PC_CPL(MY_DOMAIN_ID)   !<-- The P-C coupler import state
            EXP_STATE_CPL_NEST=>nmm_int_state%EXP_STATE_PC_CPL(MY_DOMAIN_ID)   !<-- The P-C coupler export state
!
            PARENT_CHILD_TIME_RATIO=>nmm_int_state%P_C_TIME_RATIO(MY_DOMAIN_ID) !<-- Ratio of this domain's timestep to its parent's
            MY_DOMAIN_MOVES=nmm_int_state%MY_DOMAIN_MOVES(MY_DOMAIN_ID)         !<-- Does this domain move?
            NEST_MODE=nmm_int_state%NEST_MODE(MY_DOMAIN_ID)                     !<-- Is this domain involved in any 2-way nesting?
            NUM_2WAY_CHILDREN=>nmm_int_state%NUM_2WAY_CHILDREN(MY_DOMAIN_ID)    !<-- How many 2-way children on this domain?
            NTRACK=nmm_int_state%NTRACK(MY_DOMAIN_ID)                           !<-- Storm locator flag
            NPHS=nmm_int_state%NPHS(MY_DOMAIN_ID)                               !<-- Physics timestep
!
!-----------------------------------------------------------------------
!***  Obtain current information from the filter clock.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="DFL: Get time info from the Clock"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_ClockGet(clock       =CLOCK_FILTER(MY_DOMAIN_ID)  &
                              ,starttime   =STARTTIME                   &
                              ,currtime    =CURRTIME                    &
                              ,advanceCount=NTIMESTEP_ESMF              &
!                             ,runTimeStepCount=NTIMESTEP_ESMF_REAL     &
                              ,runduration =RUNDURATION                 &
                              ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            NTIMESTEP=NTIMESTEP_ESMF
!           NTIMESTEP=NINT(NTIMESTEP_ESMF_REAL)-1                          !<-- ESMF timestep count starts with 1
            TIMESTEP=>nmm_int_state%TIMESTEP(MY_DOMAIN_ID)
!
!-----------------------------------------------------------------------
!
            CALL NMM_INTEGRATE(clock_direction    ='Forward '                  &  !<-- This filter only integrates forward
                              ,domain_grid_comp   =DOMAIN_GRID_COMP            &
                              ,imp_state_domain   =IMP_STATE_DOMAIN            &
                              ,exp_state_domain   =EXP_STATE_DOMAIN            &
                              ,clock_integrate    =CLOCK_FILTER(MY_DOMAIN_ID)  &
                              ,currtime           =CURRTIME                    &
                              ,starttime          =STARTTIME                   &
                              ,timestep           =FILT_TIMESTEP               &
                              ,ntimestep_ext      =NTIMESTEP                   &
                              ,runstepcount       =Domain_RunstepCount         &
                              ,dt                 =FILT_DT(MY_DOMAIN_ID)       &
                              ,filter_method      =FILTER_METHOD               &
                              ,halfdfiintval      =HALFDFIINTVAL(MY_DOMAIN_ID) &
                              ,halfdfitime        =HALFDFITIME(MY_DOMAIN_ID)   &
                              ,restarted_run      =RESTARTED_RUN               &
                              ,rst_out_00         =RST_OUT_00                  &
                              ,i_am_a_fcst_task   =I_AM_A_FCST_TASK            &
                              ,i_am_lead_fcst_task=I_AM_LEAD_FCST_TASK         &
                              ,nesting            =NESTING_NMM                 &
                              ,nest_mode          =NEST_MODE                   &
                              ,task_mode          =TASK_MODE                   &
                              ,i_am_a_nest        =I_AM_A_NEST                 &
                              ,my_domain_id       =MY_DOMAIN_ID                &
                              ,num_children       =NUM_CHILDREN                &
                              ,num_2way_children  =NUM_2WAY_CHILDREN           &
                              ,parent_child_cpl   =PARENT_CHILD_COUPLER_COMP   &
                              ,imp_state_cpl_nest =IMP_STATE_CPL_NEST          &
                              ,exp_state_cpl_nest =EXP_STATE_CPL_NEST          &
                              ,par_chi_time_ratio =PARENT_CHILD_TIME_RATIO     &
                              ,my_domain_moves    =MY_DOMAIN_MOVES             &
                              ,ntrack             =NTRACK                      &
                              ,nphs               =NPHS                        &
                              ,last_generation    =LAST_GENERATION             &
                              ,advanced           =ADVANCED                    &
                              ,mype               =MYPE_LOCAL                  &
                              ,comm_global        =COMM_GLOBAL                 &
                              ,generation_finished=GENERATION_FINISHED(N)      &
                              ,timers_domain      =TIMERS(MY_DOMAIN_ID)        &
                              ,npe_print          =NPE_PRINT                   &
                              ,print_timing       =PRINT_TIMING )
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  All tasks that are finished on all generations may leave.
!-----------------------------------------------------------------------
!
            IF(ALL(GENERATION_FINISHED,NUM_GENS))THEN                      !<-- If true, all of this task's domains are finished
              ALL_FORECASTS_COMPLETE=.TRUE.
              EXIT gens_f1_2
            ENDIF
!
!-----------------------------------------------------------------------
!
          ENDDO gens_f1_2
!
!-----------------------------------------------------------------------
!
        ENDDO  dfl_int
!
!-----------------------------------------------------------------------
!***  The filter is now finished integrating.  Reset the actual
!***  integration clock.
!-----------------------------------------------------------------------
!
        gens_f1_3: DO N=1,NUM_GENS   
!
          MY_DOMAIN_ID=MY_DOMAINS_IN_GENS(N)                               !<-- This task's (only) domain in generation N
          IF(MY_DOMAIN_ID==0)CYCLE                                         !<-- This task is not on a domain in generation N
!
          I_AM_LEAD_FCST_TASK=>nmm_int_state%I_AM_LEAD_FCST_TASK(MY_DOMAIN_ID) !<-- Is this the lead fcst task on this domain?
          FILT_TIMESTEP=>nmm_int_state%FILT_TIMESTEP(MY_DOMAIN_ID)         !<-- This domain's timestep for digital filtering
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="DFL: Get CURRTIME from Filter Clock When Finished"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_ClockGet(clock           =CLOCK_FILTER(MY_DOMAIN_ID)    &
                            ,currtime        =CURRTIME                      &
!                           ,advanceCount    =NTIMESTEP_ESMF                &
                            ,runTimeStepCount=NTIMESTEP_ESMF_REAL       &
                            ,rc              =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_TimeGet(CURRTIME, dd=DD, h=H, m=M, s=S, rc=RC)
          CURRTIME=CURRTIME+FILT_TIMESTEP
!
          CALL ESMF_TimeGet(CURRTIME, dd=DD, h=H, m=M, s=S, rc=RC)
          STARTTIME=CURRTIME-HALFDFIINTVAL(MY_DOMAIN_ID)                   !<-- Start time set to halfway point of filter period
          CALL ESMF_TimeGet(STARTTIME, dd=DD, h=H, m=M, s=S, rc=RC)
!
!         NTIMESTEP=NTIMESTEP_ESMF
          NTIMESTEP=NINT(NTIMESTEP_ESMF_REAL)-1
          NTIMESTEP=NTIMESTEP+1
!         NTIMESTEP_ESMF=NTIMESTEP
          NTIMESTEP_ESMF_REAL=REAL(NTIMESTEP)+1
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Set Clock After DFL Filter"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_ClockSet(clock           =CLOCK_NMM(MY_DOMAIN_ID)   &  !<-- For DFL filter, the starttime of the free forecast
!                           ,starttime       =STARTTIME                 &  !    moves ahead to the halfway point of the filter
                            ,currtime        =CURRTIME                  &  !    interval.
!                           ,advancecount    =NTIMESTEP_ESMF            &
                            ,runTimeStepCount=NTIMESTEP+1               &  !<-- ESMF timestep count starts with 1, not 0
                            ,rc              =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
          IF(I_AM_LEAD_FCST_TASK) write(0,*) 'steps to increment for DFL: ', RESTVAL/DT(MY_DOMAIN_ID)
!
!nuopc    NTIMESTEP_ESMF=NTIMESTEP_ESMF*(FILT_DT(MY_DOMAIN_ID)/DT(MY_DOMAIN_ID))+0.1
!nuopc    NTIMESTEP_ESMF=NTIMESTEP_ESMF + (RESTVAL/DT(MY_DOMAIN_ID))
          NTIMESTEP=NTIMESTEP*(FILT_DT(MY_DOMAIN_ID)/DT(MY_DOMAIN_ID))+0.1
          NTIMESTEP=NTIMESTEP + (RESTVAL/DT(MY_DOMAIN_ID))
!
          CALL ESMF_ClockSet(clock           =CLOCK_NMM(MY_DOMAIN_ID)   &  !<-- The NEMS ESMF Clock
                            ,starttime       =STARTTIME                 &  !<-- The simulation start time (ESMF)
!                           ,advanceCount    =NTIMESTEP_ESMF            &
                            ,runTimeStepCount=NTIMESTEP+1               &
                            ,rc              =RC)
!
!-----------------------------------------------------------------------
!
          IF(I_AM_LEAD_FCST_TASK)THEN
            WRITE(0,*)' Completed filter method ',filter_method
            WRITE(0,*)' Now reset filter method to 0.'
          ENDIF
!
!-----------------------------------------------------------------------
!
        ENDDO gens_f1_3
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      ELSEIF(FILTER_METHOD==2)THEN  method_block                           !<-- The DDFI digital filter.
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  First loop through the generations to set preliminary variables
!***  specific to the domains in them.
!-----------------------------------------------------------------------
!
        gens_f2_1: DO N=1,NUM_GENS
!
          MY_DOMAIN_ID=MY_DOMAINS_IN_GENS(N)                               !<-- This task's (only) domain in generation N
          IF(MY_DOMAIN_ID==0)CYCLE                                         !<-- This task is not on a domain in generation N
!
!-----------------------------------------------------------------------
!
          FILT_TIMESTEP=>nmm_int_state%FILT_TIMESTEP(MY_DOMAIN_ID)         !<-- This domain's timestep for digital filtering
!
          I_AM_LEAD_FCST_TASK=>nmm_int_state%I_AM_LEAD_FCST_TASK(MY_DOMAIN_ID) !<-- Is this the lead fcst task on this domain?
!
          IF(I_AM_LEAD_FCST_TASK)WRITE(0,*)' Beginning DDFI Filter'
!
!-----------------------------------------------------------------------
!
!--------------------------------
!***  The initial backward step.
!--------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Digital Filter: Extract DFIHR Value for DDFI"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)          &  !<-- The config object
                                      ,value = DFIHR                    &  !<-- Half foward filter window (s)
                                      ,label ='nsecs_bckddfi:'          &  !<-- Time duration of this backward part of filter
                                      ,rc    =RC)
!
          CALL ESMF_TimeIntervalSet(timeinterval=HALFDFIINTVAL(MY_DOMAIN_ID) &
                                   ,s           =DFIHR                       &
                                   ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Is DFIHR divided evenly by the timestep?  We cannot proceed
!***  unless it is.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="For DDFI Get Actual Timestep from ESMF Value"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_TimeIntervalGet(timeinterval=FILT_TIMESTEP          &  !<-- The filter timestep on this domain (sec) (ESMF)
                                   ,s           =S                      &
                                   ,sn          =Sn                     &
                                   ,sd          =Sd                     &
                                   ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          NDFISTEP(MY_DOMAIN_ID)=INT( 0.1+DFIHR/(S+REAL(Sn)/REAL(Sd)))
          DFIHR_CHK=INT(0.1+NDFISTEP(MY_DOMAIN_ID)*(S+REAL(Sn)/REAL(Sd)))
!
          IF (DFIHR_CHK /= DFIHR) THEN
            WRITE(0,*)' DFIHR=',DFIHR,' DFIHR_CHK=',DFIHR_CHK,' on domain #',my_domain_id
            WRITE(0,*)'nsecs_bckddfi in configure MUST be integer multiple of the timestep'
            WRITE(0,*)' User must reset the value'
            WRITE(0,*)' *** ABORTING MODEL RUN *** '
            CALL ESMF_Finalize(RC=RC,endflag=ESMF_END_ABORT)
          ENDIF
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="DDFI Filter: Get current time info from NMM Clock"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_ClockGet(clock   =CLOCK_NMM(MY_DOMAIN_ID)           &
                            ,currtime=CURRTIME                          &
                            ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          STARTTIME=CURRTIME
          HALFDFITIME(MY_DOMAIN_ID)=STARTTIME-HALFDFIINTVAL(MY_DOMAIN_ID)
          DFITIME=HALFDFITIME(MY_DOMAIN_ID)
!
          TIMESTEP_FILTER=-FILT_TIMESTEP                                   !<-- Prepare for backward part of integration
!
!-----------------------------------------------------------------------
!***  In preparation for this filter's forward integration
!***  create a clock to control the filter's timestepping.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Create the Clock for the DDFI Digital Filter."
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CLOCK_FILTER(MY_DOMAIN_ID)=ESMF_ClockCreate(name     ='CLOCK_DDFI'      &  !<-- The Clock for the DDFI filter
                                                     ,timeStep =TIMESTEP_FILTER   &  !<-- The filter timestep in this component
                                                     ,startTime=STARTTIME         &  !<-- Start time of filter
                                                     ,stopTime =DFITIME           &  !<-- Stop time of the filter
                                                     ,rc       =RC)
!
          CALL ESMF_ClockSet(clock    =CLOCK_FILTER(MY_DOMAIN_ID)       &
                            ,currtime =CURRTIME                         &
                            ,starttime=CURRTIME                         &
                            ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          HDIFF_ON=0                                                       !<-- Turn off horiz diffusion for backward integration.
          MEAN_ON =0
!
          IMP_STATE_DOMAIN=>nmm_int_state%IMP_STATE_DOMAIN(MY_DOMAIN_ID)
!
          IF(I_AM_LEAD_FCST_TASK)WRITE(0,*)' Set Clock direction to backward for DDFI'
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="For DDFI Set Import State Attributes for Backward Integration"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeSet(state=IMP_STATE_DOMAIN                 &  !<-- This DOMAIN component's import state
                                ,name ='Clock_Direction'                &
                                ,value='Bckward '                       &
                                ,rc   =RC)
!
          CALL ESMF_AttributeSet(state=IMP_STATE_DOMAIN                 &  !<-- This DOMAIN component's import state for filter
                                ,name ='HDIFF'                          &  !<-- Flag for horizontal diffusion on/off
                                ,value=HDIFF_ON                         &  !<-- Value of horizontal diffusion flag
                                ,rc   =RC)
!
          CALL ESMF_AttributeSet(state=IMP_STATE_DOMAIN                 &  !<-- This DOMAIN component's import state for filter
                                ,name ='MEAN_ON'                        &
                                ,value=MEAN_ON                          &
                                ,rc   =RC)
!
          CALL ESMF_AttributeSet(state=IMP_STATE_DOMAIN                 &  !<-- This DOMAIN component's import state for filter
                                ,name ='NDFISTEP'                       &
                                ,value=NDFISTEP(MY_DOMAIN_ID)           &
                                ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
        ENDDO gens_f2_1
!
!-----------------------------------------------------------------------
!***  Execute the DDFI filter's backward integration for all domains 
!***  after dereferencing argument variables for the given domain.
!***  See fuller explanation in subroutine NMM_RUN.
!-----------------------------------------------------------------------
!
        ddfi_backward: DO WHILE(.NOT.ALL_FORECASTS_COMPLETE)
!
!-----------------------------------------------------------------------
!
          gens_f2_2: DO N=1,NUM_GENS
!
!-----------------------------------------------------------------------
!
            IF(GENERATION_FINISHED(N))THEN
              CYCLE gens_f2_2
            ENDIF
!
            LAST_GENERATION=.FALSE.
            IF(N==NUM_GENS)LAST_GENERATION=.TRUE.
!
            MY_DOMAIN_ID=MY_DOMAINS_IN_GENS(N)                             !<-- This task's (only) domain in generation N
            IF(MY_DOMAIN_ID==0)CYCLE                                       !<-- This task is not on a domain in generation N
!
            DOMAIN_GRID_COMP=>nmm_int_state%DOMAIN_GRID_COMP(MY_DOMAIN_ID) !<-- This domain's ESMF component
!
!-----------------------------------------------------------------------
!***  We need the task's rank on the current domain.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="NMM_RUN: Retrieve VM from DOMAIN Gridded Component"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_GridCompGet(gridcomp=DOMAIN_GRID_COMP             &  !<-- The DOMAIN gridded component
                                 ,vm      =VM                           &  !<-- Get the Virtual Machine from the DOMAIN component
                                 ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="NMM_RUN: Obtain the Local Task ID"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_VMGet(vm      =VM                                 &  !<-- The virtual machine for this DOMAIN component
                           ,localpet=MYPE_LOCAL                         &  !<-- Each task's local rank on this domain
                           ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
            IMP_STATE_DOMAIN=>nmm_int_state%IMP_STATE_DOMAIN(MY_DOMAIN_ID) !<-- Its import state
            EXP_STATE_DOMAIN=>nmm_int_state%EXP_STATE_DOMAIN(MY_DOMAIN_ID) !<-- Its export state
!
            I_AM_A_FCST_TASK   =>nmm_int_state%I_AM_A_FCST_TASK(MY_DOMAIN_ID)    !<-- Is this task a fcst task on this domain?
            I_AM_LEAD_FCST_TASK=>nmm_int_state%I_AM_LEAD_FCST_TASK(MY_DOMAIN_ID) !<-- Is this the lead fcst task on this domain?
!
            FILT_TIMESTEP=>nmm_int_state%FILT_TIMESTEP(MY_DOMAIN_ID)       !<-- This domain's timestep for digital filtering
!
            COMM_TO_MY_PARENT=>comms_domain(MY_DOMAIN_ID)%TO_PARENT        !<-- This domain's intracommunicator to its parent
            NUM_CHILDREN=>nmm_int_state%NUM_CHILDREN(MY_DOMAIN_ID)         !<-- How many children does this domain have?
!
            NPE_PRINT=>nmm_int_state%NPE_PRINT(MY_DOMAIN_ID)               !<-- Print clocktimes from this task
!
            RESTARTED_RUN=>nmm_int_state%RESTARTED_RUN(MY_DOMAIN_ID)
            RST_OUT_00=>nmm_int_state%RST_OUT_00(MY_DOMAIN_ID)
!
            I_AM_A_NEST=>nmm_int_state%I_AM_A_NEST(MY_DOMAIN_ID)           !<-- Is this domain a nest?
!
            PARENT_CHILD_COUPLER_COMP=>nmm_int_state%PC_CPL_COMP(MY_DOMAIN_ID) !<-- The P-C coupler associated with this domain
            IMP_STATE_CPL_NEST=>nmm_int_state%IMP_STATE_PC_CPL(MY_DOMAIN_ID)   !<-- The P-C coupler import state
            EXP_STATE_CPL_NEST=>nmm_int_state%EXP_STATE_PC_CPL(MY_DOMAIN_ID)   !<-- The P-C coupler export state
!
            PARENT_CHILD_TIME_RATIO=>nmm_int_state%P_C_TIME_RATIO(MY_DOMAIN_ID) !<-- Ratio of this domain's timestep to its parent's
            MY_DOMAIN_MOVES=nmm_int_state%MY_DOMAIN_MOVES(MY_DOMAIN_ID)         !<-- Does this domain move?
            NEST_MODE=nmm_int_state%NEST_MODE(MY_DOMAIN_ID)                     !<-- Is this domain involved in any 2-way nesting?
            NUM_2WAY_CHILDREN=>nmm_int_state%NUM_2WAY_CHILDREN(MY_DOMAIN_ID)    !<-- How many 2-way children on this domain?
            NTRACK=nmm_int_state%NTRACK(MY_DOMAIN_ID)                           !<-- Storm locator flag
            NPHS=nmm_int_state%NPHS(MY_DOMAIN_ID)                               !<-- Physics timestep
!
!-----------------------------------------------------------------------
!***  Obtain current information from the filter clock.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="DDFI Backward: Get time info from the Clock"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_ClockGet(clock           =CLOCK_FILTER(MY_DOMAIN_ID) &
                              ,starttime       =STARTTIME                  &
                              ,currtime        =CURRTIME                   &
!                             ,advanceCount    =NTIMESTEP_ESMF             &
                              ,runTimeStepCount=NTIMESTEP_ESMF_REAL        &
                              ,runduration     =RUNDURATION                &
                              ,rc              =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!           NTIMESTEP=NTIMESTEP_ESMF
            NTIMESTEP=NINT(NTIMESTEP_ESMF_REAL)-1
!
!-----------------------------------------------------------------------
!
            CALL NMM_INTEGRATE(clock_direction    ='Bckward '                 &  !<-- The initial backward piece of the filter
                              ,domain_grid_comp   =DOMAIN_GRID_COMP           &
                              ,imp_state_domain   =IMP_STATE_DOMAIN           &
                              ,exp_state_domain   =EXP_STATE_DOMAIN           &
                              ,clock_integrate    =CLOCK_FILTER(MY_DOMAIN_ID) &
                              ,currtime           =CURRTIME                   &
                              ,starttime          =STARTTIME                  &
                              ,timestep           =FILT_TIMESTEP              &
                              ,ntimestep_ext      =NTIMESTEP                  &
                              ,runstepcount       =Domain_RunstepCount        &
                              ,dt                 =FILT_DT(MY_DOMAIN_ID)      &
                              ,filter_method      =FILTER_METHOD              &
                              ,ndfistep           =NDFISTEP(MY_DOMAIN_ID)     &
                              ,restarted_run      =RESTARTED_RUN              &
                              ,rst_out_00         =RST_OUT_00                 &
                              ,i_am_a_fcst_task   =I_AM_A_FCST_TASK           &
                              ,i_am_lead_fcst_task=I_AM_LEAD_FCST_TASK        &
                              ,nesting            =NESTING_NMM                &
                              ,nest_mode          =NEST_MODE                  &
                              ,task_mode          =TASK_MODE                  &
                              ,i_am_a_nest        =I_AM_A_NEST                &
                              ,my_domain_id       =MY_DOMAIN_ID               &
                              ,num_children       =NUM_CHILDREN               &
                              ,num_2way_children  =NUM_2WAY_CHILDREN          &
                              ,parent_child_cpl   =PARENT_CHILD_COUPLER_COMP  &
                              ,imp_state_cpl_nest =IMP_STATE_CPL_NEST         &
                              ,exp_state_cpl_nest =EXP_STATE_CPL_NEST         &
                              ,par_chi_time_ratio =PARENT_CHILD_TIME_RATIO    &
                              ,my_domain_moves    =MY_DOMAIN_MOVES            &
                              ,ntrack             =NTRACK                     &
                              ,nphs               =NPHS                       &
                              ,last_generation    =LAST_GENERATION            &
                              ,advanced           =ADVANCED                   &
                              ,mype               =MYPE_LOCAL                 &
                              ,comm_global        =COMM_GLOBAL                &
                              ,generation_finished=GENERATION_FINISHED(N)     &
                              ,timers_domain      =TIMERS(MY_DOMAIN_ID)       &
                              ,npe_print          =NPE_PRINT                  &
                              ,print_timing       =PRINT_TIMING )
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  All tasks that are finished on all generations may leave.
!-----------------------------------------------------------------------
!
            IF(ALL(GENERATION_FINISHED,NUM_GENS))THEN                        !<-- If true, all of this task's domains are finished
              ALL_FORECASTS_COMPLETE=.TRUE.
              EXIT gens_f2_2
            ENDIF
!
!-----------------------------------------------------------------------
!
          ENDDO gens_f2_2
!
!-----------------------------------------------------------------------
!
        ENDDO ddfi_backward
!
!-----------------------------------------------------------------------
!***  Reset the completion flags for the forward integration.
!-----------------------------------------------------------------------
!
        IF(NUM_GENS==1)THEN
          GENERATION_FINISHED(1)=.FALSE.
!
        ELSE
          DO N=1,NUM_GENS
            IF(MY_DOMAINS_IN_GENS(N)>0)THEN
              GENERATION_FINISHED(N)=.FALSE.
            ELSE
              GENERATION_FINISHED(N)=.TRUE.                                  !<-- Task not in this generation; consider it finished.
            ENDIF
          ENDDO
!
        ENDIF
!
        ALL_FORECASTS_COMPLETE=.FALSE.
!
!-----------------------------------------------------------------------
!***  Prepare to do the final forward step of the DDFI filter.
!-----------------------------------------------------------------------
!
        gens_f2_3: DO N=1,NUM_GENS
!
          MY_DOMAIN_ID=MY_DOMAINS_IN_GENS(N)                               !<-- This task's (only) domain in generation N
          IF(MY_DOMAIN_ID==0)CYCLE                                         !<-- This task is not on a domain in generation N
!
!-----------------------------------------------------------------------
!
          NTIMESTEP=0
!nuopc    NTIMESTEP_ESMF=NTIMESTEP
          NTIMESTEP_ESMF_REAL=REAL(NTIMESTEP)+1                            !<-- ESMF timestep count starts with 1, not 0
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="For DDFI Get DFIHR for Forward Integration"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)          &  !<-- The config object for this domain
                                      ,value =DFIHR                     &  !<-- Time duration of this forward integration
                                      ,label ='nsecs_fwdddfi:'          &  !<-- The configure name
                                      ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="For Forward Part of DDFI Set HALFDFIINTVAL"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_TimeIntervalSet(timeinterval=HALFDFIINTVAL(MY_DOMAIN_ID) &
                                   ,s           =DFIHR                       &
                                   ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="For Forward Part of DDFI Get Starttime"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_ClockGet(clock    =CLOCK_FILTER(MY_DOMAIN_ID)       &  !<-- The Clock for the DFI filter
                            ,startTime=STARTTIME                        &  !<-- The simulation start time (ESMF)
                            ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          HALFDFITIME(MY_DOMAIN_ID)=CURRTIME+HALFDFIINTVAL(MY_DOMAIN_ID)
          DFITIME=HALFDFITIME(MY_DOMAIN_ID)+HALFDFIINTVAL(MY_DOMAIN_ID)
!
          FILT_TIMESTEP=>nmm_int_state%FILT_TIMESTEP(MY_DOMAIN_ID)         !<-- This domain's timestep for digital filtering
          TIMESTEP_FILTER=FILT_TIMESTEP                                    !<-- Prepare for forward part of integration
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Reset the Clock for Forward DDFI Digital Filter."
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_ClockSet(clock           =CLOCK_FILTER(MY_DOMAIN_ID) &  !<-- Reset the stoptime for the forward part of the filter
                            ,timeStep        =TIMESTEP_FILTER            &  !<-- The fundamental timestep in this component
                            ,starttime       =CURRTIME                   &  !<-- Start backward integration at current time
                            ,stoptime        =DFITIME                    &  !<-- End backward integration at DFITIME
!                           ,advancecount    =NTIMESTEP_ESMF             &
                            ,runTimeStepCount=NTIMESTEP+1                &
                            ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          HDIFF_ON=1                                                       !<-- Forward integration so we want horiz diffusion.
          MEAN_ON =0
!
          IMP_STATE_DOMAIN=>nmm_int_state%IMP_STATE_DOMAIN(MY_DOMAIN_ID)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="For DDFI Set Import State Attributes for Forward Integration"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeSet(state=IMP_STATE_DOMAIN                 &  !<-- This DOMAIN component's import state for filter
                                ,name ='Clock_Direction'                &
                                ,value='Forward '                       &
                                ,rc   =RC)
!
          CALL ESMF_AttributeSet(state=IMP_STATE_DOMAIN                 &  !<-- This DOMAIN component's import state for filter
                                ,name ='HDIFF'                          &  !<-- Flag for horizontal diffusion on/off
                                ,value=HDIFF_ON                         &  !<-- Value of horizontal diffusion flag
                                ,rc   =RC)
!
          CALL ESMF_AttributeSet(state=IMP_STATE_DOMAIN                 &  !<-- This DOMAIN component's import state for filter
                                ,name ='MEAN_ON'                        &
                                ,value=MEAN_ON                          &
                                ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
        ENDDO gens_f2_3
!
!-----------------------------------------------------------------------
!***  Now execute the final forward step of the DDFI filter.
!***  See fuller explanation in subroutine NMM_RUN.
!-----------------------------------------------------------------------
!
        ddfi_forward: DO WHILE(.NOT.ALL_FORECASTS_COMPLETE)
!
!-----------------------------------------------------------------------
!
          gens_f2_4: DO N=1,NUM_GENS
!
!-----------------------------------------------------------------------
!
            IF(GENERATION_FINISHED(N))THEN
              CYCLE gens_f2_4
            ENDIF
!
            LAST_GENERATION=.FALSE.
            IF(N==NUM_GENS)LAST_GENERATION=.TRUE.
!
            MY_DOMAIN_ID=MY_DOMAINS_IN_GENS(N)                             !<-- This task's (only) domain in generation N
            IF(MY_DOMAIN_ID==0)CYCLE                                       !<-- This task is not on a domain in generation N
!
            DOMAIN_GRID_COMP=>nmm_int_state%DOMAIN_GRID_COMP(MY_DOMAIN_ID) !<-- This domain's ESMF component
!
!-----------------------------------------------------------------------
!***  We need the task's rank on the current domain.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="NMM_RUN: Retrieve VM from DOMAIN Gridded Component"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_GridCompGet(gridcomp=DOMAIN_GRID_COMP             &  !<-- The DOMAIN gridded component
                                 ,vm      =VM                           &  !<-- Get the Virtual Machine from the DOMAIN component
                                 ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="NMM_RUN: Obtain the Local Task ID"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_VMGet(vm      =VM                                 &  !<-- The virtual machine for this DOMAIN component
                           ,localpet=MYPE_LOCAL                         &  !<-- Each task's local rank on this domain
                           ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Get Starttime and Currtime from Clock for Forward DDFI"
!           CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_ClockGet(clock         =CLOCK_FILTER(MY_DOMAIN_ID) &  !<-- The ESMF Clock for the digital filter
                            ,startTime       =STARTTIME                  &  !<-- The simulation start time (ESMF)
                            ,currTime        =CURRTIME                   &  !<-- The simulation current time (ESMF)
!!!                         ,runDuration     =RUNDURATION                &  !<-- The simulation run duration (ESMF)
!                           ,advanceCount    =NTIMESTEP_ESMF             &  !<-- Timestep count
                            ,runTimeStepCount=NTIMESTEP_ESMF_REAL        &
                            ,rc              =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!           NTIMESTEP=NTIMESTEP_ESMF
            NTIMESTEP=NINT(NTIMESTEP_ESMF_REAL)-1                          !<-- ESMF timestep count starts with 1, not 0
!
!-----------------------------------------------------------------------
!
            IMP_STATE_DOMAIN=>nmm_int_state%IMP_STATE_DOMAIN(MY_DOMAIN_ID)
            EXP_STATE_DOMAIN=>nmm_int_state%EXP_STATE_DOMAIN(MY_DOMAIN_ID) !<-- Its export state
!
            I_AM_A_FCST_TASK   =>nmm_int_state%I_AM_A_FCST_TASK(MY_DOMAIN_ID)    !<-- Is this task a fcst task on this domain?
            I_AM_LEAD_FCST_TASK=>nmm_int_state%I_AM_LEAD_FCST_TASK(MY_DOMAIN_ID) !<-- Is this the lead fcst task on this domain?
!
            FILT_TIMESTEP=>nmm_int_state%FILT_TIMESTEP(MY_DOMAIN_ID)       !<-- This domain's timestep for digital filtering
!
            COMM_TO_MY_PARENT=>comms_domain(MY_DOMAIN_ID)%TO_PARENT        !<-- This domain's intracommunicator to its parent
            NUM_CHILDREN=>nmm_int_state%NUM_CHILDREN(MY_DOMAIN_ID)         !<-- How many children does this domain have?
!
            NPE_PRINT=>nmm_int_state%NPE_PRINT(MY_DOMAIN_ID)               !<-- Print clocktimes from this task
!
            RESTARTED_RUN=>nmm_int_state%RESTARTED_RUN(MY_DOMAIN_ID)
            RST_OUT_00=>nmm_int_state%RST_OUT_00(MY_DOMAIN_ID)
!
            I_AM_A_NEST=>nmm_int_state%I_AM_A_NEST(MY_DOMAIN_ID)           !<-- Is this domain a nest?
!
            PARENT_CHILD_COUPLER_COMP=>nmm_int_state%PC_CPL_COMP(MY_DOMAIN_ID) !<-- The P-C coupler associated with this domain
            IMP_STATE_CPL_NEST=>nmm_int_state%IMP_STATE_PC_CPL(MY_DOMAIN_ID)   !<-- The P-C coupler import state
            EXP_STATE_CPL_NEST=>nmm_int_state%EXP_STATE_PC_CPL(MY_DOMAIN_ID)   !<-- The P-C coupler export state
!
            PARENT_CHILD_TIME_RATIO=>nmm_int_state%P_C_TIME_RATIO(MY_DOMAIN_ID) !<-- Ratio of this domain's timestep to its parent's
            MY_DOMAIN_MOVES=nmm_int_state%MY_DOMAIN_MOVES(MY_DOMAIN_ID)         !<-- Does this domain move?
            NEST_MODE=nmm_int_state%NEST_MODE(MY_DOMAIN_ID)                     !<-- Is this domain involved in any 2-way nesting?
            NUM_2WAY_CHILDREN=>nmm_int_state%NUM_2WAY_CHILDREN(MY_DOMAIN_ID)    !<-- How many 2-way children on this domain?
            NTRACK=nmm_int_state%NTRACK(MY_DOMAIN_ID)                           !<-- Storm locator flag
            NPHS=nmm_int_state%NPHS(MY_DOMAIN_ID)                               !<-- Physics timestep
!
!-----------------------------------------------------------------------
!
            CALL NMM_INTEGRATE(clock_direction    ='Forward '                  &  !<-- The final forward piece of the filter
                              ,domain_grid_comp   =DOMAIN_GRID_COMP            &
                              ,imp_state_domain   =IMP_STATE_DOMAIN            &
                              ,exp_state_domain   =EXP_STATE_DOMAIN            &
                              ,clock_integrate    =CLOCK_FILTER(MY_DOMAIN_ID)  &
                              ,currtime           =CURRTIME                    &
                              ,starttime          =CURRTIME                    &  !<-- CURRTIME was set or reset at end of backward piece
                              ,timestep           =FILT_TIMESTEP               &
                              ,ntimestep_ext      =NTIMESTEP                   &
                              ,ndfistep           =NDFISTEP(MY_DOMAIN_ID)      &
                              ,runstepcount       =Domain_RunstepCount         &
                              ,dt                 =FILT_DT(MY_DOMAIN_ID)       &
                              ,filter_method      =FILTER_METHOD               &
                              ,halfdfiintval      =HALFDFIINTVAL(MY_DOMAIN_ID) &
                              ,halfdfitime        =HALFDFITIME(MY_DOMAIN_ID)   &
                              ,restarted_run      =RESTARTED_RUN               &
                              ,rst_out_00         =RST_OUT_00                  &
                              ,i_am_a_fcst_task   =I_AM_A_FCST_TASK            &
                              ,i_am_lead_fcst_task=I_AM_LEAD_FCST_TASK         &
                              ,nesting            =NESTING_NMM                 &
                              ,nest_mode          =NEST_MODE                   &
                              ,task_mode          =TASK_MODE                   &
                              ,i_am_a_nest        =I_AM_A_NEST                 &
                              ,my_domain_id       =MY_DOMAIN_ID                &
                              ,num_children       =NUM_CHILDREN                &
                              ,num_2way_children  =NUM_2WAY_CHILDREN           &
                              ,parent_child_cpl   =PARENT_CHILD_COUPLER_COMP   &
                              ,imp_state_cpl_nest =IMP_STATE_CPL_NEST          &
                              ,exp_state_cpl_nest =EXP_STATE_CPL_NEST          &
                              ,par_chi_time_ratio =PARENT_CHILD_TIME_RATIO     &
                              ,my_domain_moves    =MY_DOMAIN_MOVES             &
                              ,ntrack             =NTRACK                      &
                              ,nphs               =NPHS                        &
                              ,last_generation    =LAST_GENERATION             &
                              ,advanced           =ADVANCED                    &
                              ,mype               =MYPE_LOCAL                  &
                              ,comm_global        =COMM_GLOBAL                 &
                              ,timers_domain      =TIMERS(MY_DOMAIN_ID)        &
                              ,npe_print          =NPE_PRINT                   &
                              ,print_timing       =PRINT_TIMING )
!
!-----------------------------------------------------------------------
!
            IF(ESMF_ClockIsStopTime(CLOCK_FILTER(MY_DOMAIN_ID),rc=RC))THEN
              GENERATION_FINISHED(N)=.TRUE.                                !<-- Generation N's filter has finished
!
              IF(I_AM_LEAD_FCST_TASK)THEN
                WRITE(0,*)' Completed filter method ',filter_method
                WRITE(0,*)' Now reset filter method to 0.'
              ENDIF
!
            ENDIF
!
!-----------------------------------------------------------------------
!***  All tasks that are finished on all generations may leave.
!-----------------------------------------------------------------------
!
            IF(ALL(GENERATION_FINISHED,NUM_GENS))THEN                        !<-- If true, all of this task's domains are finished
              ALL_FORECASTS_COMPLETE=.TRUE.
              EXIT gens_f2_4
            ENDIF
!
!-----------------------------------------------------------------------
!
          ENDDO gens_f2_4
!
!-----------------------------------------------------------------------
!
        ENDDO ddfi_forward
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      ELSEIF(FILTER_METHOD==3)THEN  method_block                           !<-- The TDFI digital filter.
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  First loop through the generations to set preliminary variables
!***  specific to the domains in them.
!-----------------------------------------------------------------------
!
        gens_f3_1: DO N=1,NUM_GENS
!
          MY_DOMAIN_ID=MY_DOMAINS_IN_GENS(N)                               !<-- This task's (only) domain in generation N
          IF(MY_DOMAIN_ID==0)CYCLE                                         !<-- This task is not on a domain in generation N
!
!-----------------------------------------------------------------------
!
          FILT_TIMESTEP=>nmm_int_state%FILT_TIMESTEP(MY_DOMAIN_ID)         !<-- This domain's timestep for digital filtering
!
          I_AM_LEAD_FCST_TASK=>nmm_int_state%I_AM_LEAD_FCST_TASK(MY_DOMAIN_ID) !<-- Is this the lead fcst task on this domain?
!
          IF(I_AM_LEAD_FCST_TASK)WRITE(0,*)' Beginning TDFI Filter'
!
!-----------------------------------------------------------------------
!
!--------------------------------
!***  The initial backward step.
!--------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Digital Filter: Extract DFIHR Value for TDFI"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)          &  !<-- The config object
                                      ,value = DFIHR                    &  !<-- The digital filter flag
                                      ,label ='nsecs_bcktdfi:'          &  !<-- Give this label's value to preceding variable
                                      ,rc    =RC)
!
          CALL ESMF_TimeIntervalSet(timeinterval=HALFDFIINTVAL(MY_DOMAIN_ID) &
                                   ,s           =DFIHR                       &
                                   ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="For TDFI Get Actual Timestep from ESMF Value"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_TimeIntervalGet(timeinterval=FILT_TIMESTEP          &  !<-- The fundamental timestep of this domain (sec) (ESMF)
                                   ,s           =S                      &
                                   ,sn          =Sn                     &
                                   ,sd          =Sd                     &
                                   ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          NDFISTEP(MY_DOMAIN_ID)=INT( 0.1+DFIHR/(S+REAL(Sn)/REAL(Sd)))
          DFIHR_CHK=INT(0.1+NDFISTEP(MY_DOMAIN_ID)*(REAL(S)+REAL(Sn)/REAL(Sd)))
!
          IF (DFIHR_CHK /= DFIHR) THEN
            WRITE(0,*)' DFIHR=',DFIHR,' DFIHR_CHK=',DFIHR_CHK
            WRITE(0,*)'nsecs_bcktdfi in configure MUST be integer multiple of the timestep'
            WRITE(0,*)' User must reset the value'
            WRITE(0,*)' *** ABORTING MODEL RUN *** '
            CALL ESMF_Finalize(RC=RC,endflag=ESMF_END_ABORT)
          ENDIF
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="TDFI Filter: Get current time info from NMM Clock"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_ClockGet(clock   =CLOCK_NMM(MY_DOMAIN_ID)           &
                            ,currtime=CURRTIME                          &
                            ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          STARTTIME=CURRTIME
          HALFDFITIME(MY_DOMAIN_ID)=STARTTIME-HALFDFIINTVAL(MY_DOMAIN_ID)
          DFITIME=HALFDFITIME(MY_DOMAIN_ID)-HALFDFIINTVAL(MY_DOMAIN_ID)
!
          TIMESTEP_FILTER=-FILT_TIMESTEP                                  !<-- Prepare for initial backward part of integration
!
          CALL ESMF_TimeGet(STARTTIME, dd=DD, h=H, m=M, s=S, rc=RC)
!
          IF(I_AM_LEAD_FCST_TASK)WRITE(0,*)' STARTTIME in TDFI DD H M S: ', DD, H, M, S
          IF(I_AM_LEAD_FCST_TASK)WRITE(0,*)' NDFISTEP=',NDFISTEP(MY_DOMAIN_ID),' DFIHR=',DFIHR
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Create the Clock for the TDFI Digital Filter."
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CLOCK_FILTER(MY_DOMAIN_ID)=ESMF_ClockCreate(name     ='CLOCK_TDFI'      &  !<-- The Clock for the DFI filter
                                                     ,timeStep =TIMESTEP_FILTER   &  !<-- The fundamental timestep in this component
                                                     ,startTime=STARTTIME         &  !<-- Start time of filter
!!!!!!!!                                             ,direction=ESMF_MODE_REVERSE &  !<-- Reverse the Clock for backward integration
                                                     ,stopTime =DFITIME           &  !<-- Stop time of the filter
                                                     ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          HDIFF_ON=0                                                       !<-- Turn off horiz diffusion for backward integration.
          MEAN_ON =0                                                       !<-- Turn off horiz diffusion for backward integration.
!
          IMP_STATE_DOMAIN=>nmm_int_state%IMP_STATE_DOMAIN(MY_DOMAIN_ID)   !<-- This domain's import state
!
          IF(I_AM_LEAD_FCST_TASK)WRITE(0,*)' Set Clock direction to backward for TDFI'
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="For TDFI Set Import State Attributes for Backward Integration"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeSet(state=IMP_STATE_DOMAIN                   &  !<-- This DOMAIN component's import state for filter
                                ,name ='Clock_Direction'                  &
                                ,value='Bckward '                         &
                                ,rc   =RC)
!
          CALL ESMF_AttributeSet(state=IMP_STATE_DOMAIN                   &  !<-- This DOMAIN component's import state for filter
                                ,name ='HDIFF'                            &  !<-- Flag for horizontal diffusion on/off
                                ,value=HDIFF_ON                           &  !<-- Value of horizontal diffusion flag
                                ,rc   =RC)
!
          CALL ESMF_AttributeSet(state=IMP_STATE_DOMAIN                   &  !<-- This DOMAIN component's import state for filter
                                ,name ='MEAN_ON'                          &
                                ,value=MEAN_ON                            &
                                ,rc   =RC)
!
          CALL ESMF_AttributeSet(state=IMP_STATE_DOMAIN                   &  !<-- This DOMAIN component's import state for filter
                                ,name ='NDFISTEP'                         &
                                ,value=NDFISTEP(MY_DOMAIN_ID)             &
                                ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
        ENDDO gens_f3_1
!
!-----------------------------------------------------------------------
!***  Execute the TDFI filter's backward integration for all domains
!***  after dereferencing argument variables for the given domain.
!***  See fuller explanation in subroutine NMM_RUN.
!-----------------------------------------------------------------------
!
        tdfi_backward: DO WHILE(.NOT.ALL_FORECASTS_COMPLETE)
!
!-----------------------------------------------------------------------
!
          gens_f3_2: DO N=1,NUM_GENS
!
!-----------------------------------------------------------------------
!
            IF(GENERATION_FINISHED(N))THEN
              CYCLE gens_f3_2
            ENDIF
!
            LAST_GENERATION=.FALSE.
            IF(N==NUM_GENS)LAST_GENERATION=.TRUE.
!
            MY_DOMAIN_ID=MY_DOMAINS_IN_GENS(N)                             !<-- This task's (only) domain in generation N
            IF(MY_DOMAIN_ID==0)CYCLE                                       !<-- This task is not on a domain in generation N
!
            DOMAIN_GRID_COMP=>nmm_int_state%DOMAIN_GRID_COMP(MY_DOMAIN_ID) !<-- This domain's ESMF component
!
!-----------------------------------------------------------------------
!***  We need the task's rank on the current domain.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="NMM_RUN: Retrieve VM from DOMAIN Gridded Component"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_GridCompGet(gridcomp=DOMAIN_GRID_COMP             &  !<-- The DOMAIN gridded component
                                 ,vm      =VM                           &  !<-- Get the Virtual Machine from the DOMAIN component
                                 ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="NMM_RUN: Obtain the Local Task ID"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_VMGet(vm      =VM                                 &  !<-- The virtual machine for this DOMAIN component
                           ,localpet=MYPE_LOCAL                         &  !<-- Each task's local rank on this domain
                           ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
            IMP_STATE_DOMAIN=>nmm_int_state%IMP_STATE_DOMAIN(MY_DOMAIN_ID) !<-- Its import state
            EXP_STATE_DOMAIN=>nmm_int_state%EXP_STATE_DOMAIN(MY_DOMAIN_ID) !<-- Its export state
!
            I_AM_A_FCST_TASK   =>nmm_int_state%I_AM_A_FCST_TASK(MY_DOMAIN_ID)    !<-- Is this task a fcst task on this domain?
            I_AM_LEAD_FCST_TASK=>nmm_int_state%I_AM_LEAD_FCST_TASK(MY_DOMAIN_ID) !<-- Is this the lead fcst task on this domain?
!
            FILT_TIMESTEP=>nmm_int_state%FILT_TIMESTEP(MY_DOMAIN_ID)       !<-- This domain's timestep for digital filtering
!
            COMM_TO_MY_PARENT=>comms_domain(MY_DOMAIN_ID)%TO_PARENT        !<-- This domain's intracommunicator to its parent
            NUM_CHILDREN=>nmm_int_state%NUM_CHILDREN(MY_DOMAIN_ID)         !<-- How many children does this domain have?
!
            NPE_PRINT=>nmm_int_state%NPE_PRINT(MY_DOMAIN_ID)               !<-- Print clocktimes from this task
!
            RESTARTED_RUN=>nmm_int_state%RESTARTED_RUN(MY_DOMAIN_ID)
            RST_OUT_00=>nmm_int_state%RST_OUT_00(MY_DOMAIN_ID)
!
            I_AM_A_NEST=>nmm_int_state%I_AM_A_NEST(MY_DOMAIN_ID)           !<-- Is this domain a nest?
!
            PARENT_CHILD_COUPLER_COMP=>nmm_int_state%PC_CPL_COMP(MY_DOMAIN_ID) !<-- The P-C coupler associated with this domain
            IMP_STATE_CPL_NEST=>nmm_int_state%IMP_STATE_PC_CPL(MY_DOMAIN_ID)   !<-- The P-C coupler import state
            EXP_STATE_CPL_NEST=>nmm_int_state%EXP_STATE_PC_CPL(MY_DOMAIN_ID)   !<-- The P-C coupler export state
!
            PARENT_CHILD_TIME_RATIO=>nmm_int_state%P_C_TIME_RATIO(MY_DOMAIN_ID) !<-- Ratio of this domain's timestep to its parent's
            MY_DOMAIN_MOVES=nmm_int_state%MY_DOMAIN_MOVES(MY_DOMAIN_ID)         !<-- Does this domain move?
            NEST_MODE=nmm_int_state%NEST_MODE(MY_DOMAIN_ID)                     !<-- Is this domain involved in any 2-way nesting?
            NUM_2WAY_CHILDREN=>nmm_int_state%NUM_2WAY_CHILDREN(MY_DOMAIN_ID)    !<-- How many 2-way children on this domain?
            NTRACK=nmm_int_state%NTRACK(MY_DOMAIN_ID)                       !<-- Storm locator flag
            NPHS=nmm_int_state%NPHS(MY_DOMAIN_ID)                           !<-- Physics timestep
!
!-----------------------------------------------------------------------
!***  Obtain current information from the filter clock.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="TDFI Backward: Get time info from the Clock"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_ClockGet(clock           =CLOCK_FILTER(MY_DOMAIN_ID) &
                              ,starttime       =STARTTIME                  &
                              ,currtime        =CURRTIME                   &
!                             ,advanceCount    =NTIMESTEP_ESMF             &
                              ,runTimeStepCount=NTIMESTEP_ESMF_REAL        &
                              ,rc              =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!           NTIMESTEP=NTIMESTEP_ESMF
            NTIMESTEP=NINT(NTIMESTEP_ESMF_REAL)-1
!
!-----------------------------------------------------------------------
!
            CALL NMM_INTEGRATE(clock_direction    ='Bckward '                  &  !<-- The initial backward piece of the filter
                              ,domain_grid_comp   =DOMAIN_GRID_COMP            &
                              ,imp_state_domain   =IMP_STATE_DOMAIN            &
                              ,exp_state_domain   =EXP_STATE_DOMAIN            &
                              ,clock_integrate    =CLOCK_FILTER(MY_DOMAIN_ID)  &
                              ,currtime           =CURRTIME                    &
                              ,starttime          =STARTTIME                   &
                              ,timestep           =FILT_TIMESTEP               &
                              ,ntimestep_ext      =NTIMESTEP                   &
                              ,runstepcount       =Domain_RunstepCount         &
                              ,dt                 =FILT_DT(MY_DOMAIN_ID)       &
                              ,filter_method      =FILTER_METHOD               &
                              ,halfdfiintval      =HALFDFIINTVAL(MY_DOMAIN_ID) &
                              ,ndfistep           =NDFISTEP(MY_DOMAIN_ID)      &
                              ,restarted_run      =RESTARTED_RUN               &
                              ,rst_out_00         =RST_OUT_00                  &
                              ,i_am_a_fcst_task   =I_AM_A_FCST_TASK            &
                              ,i_am_lead_fcst_task=I_AM_LEAD_FCST_TASK         &
                              ,nesting            =NESTING_NMM                 &
                              ,nest_mode          =NEST_MODE                   &
                              ,task_mode          =TASK_MODE                   &
                              ,i_am_a_nest        =I_AM_A_NEST                 &
                              ,my_domain_id       =MY_DOMAIN_ID                &
                              ,num_children       =NUM_CHILDREN                &
                              ,num_2way_children  =NUM_2WAY_CHILDREN           &
                              ,parent_child_cpl   =PARENT_CHILD_COUPLER_COMP   &
                              ,imp_state_cpl_nest =IMP_STATE_CPL_NEST          &
                              ,exp_state_cpl_nest =EXP_STATE_CPL_NEST          &
                              ,par_chi_time_ratio =PARENT_CHILD_TIME_RATIO     &
                              ,my_domain_moves    =MY_DOMAIN_MOVES             &
                              ,ntrack             =NTRACK                      &
                              ,nphs               =NPHS                        &
                              ,last_generation    =LAST_GENERATION             &
                              ,advanced           =ADVANCED                    &
                              ,mype               =MYPE_LOCAL                  &
                              ,comm_global        =COMM_GLOBAL                 &
                              ,generation_finished=GENERATION_FINISHED(N)      &
                              ,timers_domain      =TIMERS(MY_DOMAIN_ID)        &
                              ,npe_print          =NPE_PRINT                   &
                              ,print_timing       =PRINT_TIMING )
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  All tasks that are finished on all generations may leave.
!-----------------------------------------------------------------------
!
            IF(ALL(GENERATION_FINISHED,NUM_GENS))THEN                      !<-- If true, all of this task's domains are finished
              ALL_FORECASTS_COMPLETE=.TRUE.
              EXIT gens_f3_2
            ENDIF
!
!-----------------------------------------------------------------------
!
          ENDDO gens_f3_2
!
!-----------------------------------------------------------------------
!
        ENDDO tdfi_backward
!
!-----------------------------------------------------------------------
!***  Reset the completion flags for the forward integration.
!-----------------------------------------------------------------------
!
        IF(NUM_GENS==1)THEN
          GENERATION_FINISHED(1)=.FALSE.
!
        ELSE
          DO N=1,NUM_GENS
            IF(MY_DOMAINS_IN_GENS(N)>0)THEN
              GENERATION_FINISHED(N)=.FALSE.
            ELSE
              GENERATION_FINISHED(N)=.TRUE.                                  !<-- Task not in this generation; consider it finished.
            ENDIF
          ENDDO
!
        ENDIF
!
        ALL_FORECASTS_COMPLETE=.FALSE.
!
!-----------------------------------------------------------------------
!***  Prepare to do the final forward step of the TDFI filter.
!-----------------------------------------------------------------------
!
        gens_f3_3: DO N=1,NUM_GENS
!
          MY_DOMAIN_ID=MY_DOMAINS_IN_GENS(N)                               !<-- This task's (only) domain in generation N
          IF(MY_DOMAIN_ID==0)CYCLE                                         !<-- This task is not on a domain in generation N
!
          IMP_STATE_DOMAIN=>nmm_int_state%IMP_STATE_DOMAIN(MY_DOMAIN_ID)   !<-- Its import state
          FILT_TIMESTEP=>nmm_int_state%FILT_TIMESTEP(MY_DOMAIN_ID)         !<-- This domain's timestep for digital filtering
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="For TDFI Get DFIHR for Forward Integration"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)          &  !<-- The config object
                                      ,value =DFIHR                     &  !<-- The digital filter duration time
                                      ,label ='nsecs_fwdtdfi:'          &  !<-- Give this label's value to preceding variable
                                      ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="For Forward Part of TDFI Set HALFDFIINTVAL"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_TimeIntervalSet(timeinterval=HALFDFIINTVAL(MY_DOMAIN_ID) &
                                   ,s           =DFIHR                       &
                                   ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="For Forward Part of TDFI Set Clock Direction"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeSet(state=IMP_STATE_DOMAIN                 &  !<-- This DOMAIN component's import state for filter
                                ,name ='Clock_Direction'                &
                                ,value='Forward '                       &
                                ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          HALFDFITIME(MY_DOMAIN_ID)=CURRTIME+HALFDFIINTVAL(MY_DOMAIN_ID)
          SDFITIME=CURRTIME
          DFITIME=HALFDFITIME(MY_DOMAIN_ID)+HALFDFIINTVAL(MY_DOMAIN_ID)
!
          TIMESTEP_FILTER=FILT_TIMESTEP                                    !<-- Prepare for forward part of integration
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Reset the Clock for Forward TDFI Digital Filter."
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_TimeIntervalGet(timeinterval=TIMESTEP_FILTER        &
                                   ,s           =S                      &
                                   ,rc          =RC)
!
          NTIMESTEP=0
!         NTIMESTEP_ESMF=NTIMESTEP
          NTIMESTEP_ESMF_REAL=REAL(NTIMESTEP)+1                            !<-- ESMF timestep count starts with 1, not 0

!         NTIMESTEP=-NTIMESTEP/2
!         NTIMESTEP_ESMF=NTIMESTEP
!
          CALL ESMF_ClockSet(clock           =CLOCK_FILTER(MY_DOMAIN_ID) &  !<-- Reset the stoptime for the forward part of the filter 
                            ,timeStep        =TIMESTEP_FILTER            &  !<-- The fundamental timestep in this component
                            ,starttime       =CURRTIME                   &  !<-- Start forward integration at the current time
                            ,currtime        =CURRTIME                   &
!                           ,advancecount    =NTIMESTEP_ESMF             &
                            ,runTimeStepCount=NTIMESTEP+1                &  !<-- ESMF timestep count starts with 1, not 0
                            ,stoptime        =DFITIME                    &  !<-- Stop forward integration at DFITIME
                            ,rc              =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          HDIFF_ON=1                                                       !<-- Forward integration so we want horiz diffusion.
          MEAN_ON =0
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="For TDFI Set Import State Attributes for Forward Integration"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeSet(state=IMP_STATE_DOMAIN                 &  !<-- This DOMAIN component's import state for filter
                                ,name ='HDIFF'                          &  !<-- Flag for horizontal diffusion on/off
                                ,value=HDIFF_ON                         &  !<-- Value of horizontal diffusion flag
                                ,rc   =RC)
!
          CALL ESMF_AttributeSet(state=IMP_STATE_DOMAIN                 &  !<-- This DOMAIN component's import state for filter
                                ,name ='MEAN_ON'                        &
                                ,value=MEAN_ON                          &
                                ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
        ENDDO gens_f3_3
!
!-----------------------------------------------------------------------
!***  Now execute the final forward step of the TDFI filter.
!***  See fuller explanation in subroutine NMM_RUN.
!-----------------------------------------------------------------------
!
        tdfi_forward: DO WHILE(.NOT.ALL_FORECASTS_COMPLETE)
!
!-----------------------------------------------------------------------
!
          gens_f3_4: DO N=1,NUM_GENS
!
!-----------------------------------------------------------------------
!
            IF(GENERATION_FINISHED(N))THEN
              CYCLE gens_f3_4
            ENDIF
!
            LAST_GENERATION=.FALSE.
            IF(N==NUM_GENS)LAST_GENERATION=.TRUE.
!
            MY_DOMAIN_ID=MY_DOMAINS_IN_GENS(N)                             !<-- This task's (only) domain in generation N
            IF(MY_DOMAIN_ID==0)CYCLE                                       !<-- This task is not on a domain in generation N
!
            DOMAIN_GRID_COMP=>nmm_int_state%DOMAIN_GRID_COMP(MY_DOMAIN_ID) !<-- This domain's ESMF component
!
!-----------------------------------------------------------------------
!***  We need the task's rank on the current domain.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="NMM_RUN: Retrieve VM from DOMAIN Gridded Component"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_GridCompGet(gridcomp=DOMAIN_GRID_COMP             &  !<-- The DOMAIN gridded component
                                 ,vm      =VM                           &  !<-- Get the Virtual Machine from the DOMAIN component
                                 ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="NMM_RUN: Obtain the Local Task ID"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_VMGet(vm      =VM                                 &  !<-- The virtual machine for this DOMAIN component
                           ,localpet=MYPE_LOCAL                         &  !<-- Each task's local rank on this domain
                           ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Get Starttime and Currtime from Clock for Forward TDFI"
!           CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_ClockGet(clock         =CLOCK_FILTER(MY_DOMAIN_ID) &  !<-- The ESMF Clock for the digital filter
                            ,startTime       =STARTTIME                  &  !<-- The simulation start time (ESMF)
                            ,currTime        =CURRTIME                   &  !<-- The simulation current time (ESMF)
!!!                         ,runDuration     =RUNDURATION                &  !<-- The simulation run duration (ESMF)
!                           ,advanceCount    =NTIMESTEP_ESMF             &  !<-- Timestep count
                            ,runTimeStepCount=NTIMESTEP_ESMF_REAL        &
                            ,rc              =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!           NTIMESTEP=NTIMESTEP_ESMF
            NTIMESTEP=NINT(NTIMESTEP_ESMF_REAL)-1
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="TDFI Forward: Get Actual Timestep from ESMF Timestep"
!           CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_TimeIntervalGet(timeinterval=TIMESTEP             &
                                     ,s           =S                    &
                                     ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------   
!
            IMP_STATE_DOMAIN=>nmm_int_state%IMP_STATE_DOMAIN(MY_DOMAIN_ID)
            EXP_STATE_DOMAIN=>nmm_int_state%EXP_STATE_DOMAIN(MY_DOMAIN_ID) !<-- Its export state
!
            I_AM_A_FCST_TASK   =>nmm_int_state%I_AM_A_FCST_TASK(MY_DOMAIN_ID)    !<-- Is this task a fcst task on this domain?
            I_AM_LEAD_FCST_TASK=>nmm_int_state%I_AM_LEAD_FCST_TASK(MY_DOMAIN_ID) !<-- Is this the lead fcst task on this domain?
!
            FILT_TIMESTEP=>nmm_int_state%FILT_TIMESTEP(MY_DOMAIN_ID)       !<-- This domain's timestep for digital filtering
!
            COMM_TO_MY_PARENT=>comms_domain(MY_DOMAIN_ID)%TO_PARENT        !<-- This domain's intracommunicator to its parent
            NUM_CHILDREN=>nmm_int_state%NUM_CHILDREN(MY_DOMAIN_ID)         !<-- How many children does this domain have?
!
            NPE_PRINT=>nmm_int_state%NPE_PRINT(MY_DOMAIN_ID)               !<-- Print clocktimes from this task
!
            RESTARTED_RUN=>nmm_int_state%RESTARTED_RUN(MY_DOMAIN_ID)
            RST_OUT_00=>nmm_int_state%RST_OUT_00(MY_DOMAIN_ID)
!
            I_AM_A_NEST=>nmm_int_state%I_AM_A_NEST(MY_DOMAIN_ID)           !<-- Is this domain a nest?
!
            PARENT_CHILD_COUPLER_COMP=>nmm_int_state%PC_CPL_COMP(MY_DOMAIN_ID) !<-- The P-C coupler associated with this domain
            IMP_STATE_CPL_NEST=>nmm_int_state%IMP_STATE_PC_CPL(MY_DOMAIN_ID)   !<-- The P-C coupler import state
            EXP_STATE_CPL_NEST=>nmm_int_state%EXP_STATE_PC_CPL(MY_DOMAIN_ID)   !<-- The P-C coupler export state
!
            PARENT_CHILD_TIME_RATIO=>nmm_int_state%P_C_TIME_RATIO(MY_DOMAIN_ID) !<-- Ratio of this domain's timestep to its parent's
            MY_DOMAIN_MOVES=nmm_int_state%MY_DOMAIN_MOVES(MY_DOMAIN_ID)         !<-- Does this domain move?
            NEST_MODE=nmm_int_state%NEST_MODE(MY_DOMAIN_ID)                     !<-- Is this domain involved in any 2-way nesting?
            NUM_2WAY_CHILDREN=>nmm_int_state%NUM_2WAY_CHILDREN(MY_DOMAIN_ID)    !<-- How many 2-way children on this domain?
            NTRACK=nmm_int_state%NTRACK(MY_DOMAIN_ID)                           !<-- Storm locator flag
            NPHS=nmm_int_state%NPHS(MY_DOMAIN_ID)                               !<-- Physics timestep
!
!-----------------------------------------------------------------------
!
            CALL NMM_INTEGRATE(clock_direction    ='Forward '                  &  !<-- The final forward piece of the filter
                              ,domain_grid_comp   =DOMAIN_GRID_COMP            &
                              ,imp_state_domain   =IMP_STATE_DOMAIN            &
                              ,exp_state_domain   =EXP_STATE_DOMAIN            &
                              ,clock_integrate    =CLOCK_FILTER(MY_DOMAIN_ID)  &
                              ,currtime           =CURRTIME                    &
                              ,starttime          =CURRTIME                    &  !<-- CURRTIME was set or reset at end of backwward piece
                              ,timestep           =FILT_TIMESTEP               &
                              ,ntimestep_ext      =NTIMESTEP                   &
                              ,runstepcount       =Domain_RunstepCount         &
                              ,ndfistep           =NDFISTEP(MY_DOMAIN_ID)      &
                              ,dt                 =FILT_DT(MY_DOMAIN_ID)       &
                              ,filter_method      =FILTER_METHOD               &
                              ,halfdfiintval      =HALFDFIINTVAL(MY_DOMAIN_ID) &
                              ,halfdfitime        =HALFDFITIME(MY_DOMAIN_ID)   &
                              ,restarted_run      =RESTARTED_RUN               &
                              ,rst_out_00         =RST_OUT_00                  &
                              ,i_am_a_fcst_task   =I_AM_A_FCST_TASK            &
                              ,i_am_lead_fcst_task=I_AM_LEAD_FCST_TASK         &
                              ,nesting            =NESTING_NMM                 &
                              ,nest_mode          =NEST_MODE                   &
                              ,task_mode          =TASK_MODE                   &
                              ,i_am_a_nest        =I_AM_A_NEST                 &
                              ,my_domain_id       =MY_DOMAIN_ID                &
                              ,num_children       =NUM_CHILDREN                &
                              ,num_2way_children  =NUM_2WAY_CHILDREN           &
                              ,parent_child_cpl   =PARENT_CHILD_COUPLER_COMP   &
                              ,imp_state_cpl_nest =IMP_STATE_CPL_NEST          &
                              ,exp_state_cpl_nest =EXP_STATE_CPL_NEST          &
                              ,par_chi_time_ratio =PARENT_CHILD_TIME_RATIO     &
                              ,my_domain_moves    =MY_DOMAIN_MOVES             &
                              ,ntrack             =NTRACK                      &
                              ,nphs               =NPHS                        &
                              ,last_generation    =LAST_GENERATION             &
                              ,advanced           =ADVANCED                    &
                              ,mype               =MYPE_LOCAL                  &
                              ,comm_global        =COMM_GLOBAL                 &
                              ,timers_domain      =TIMERS(MY_DOMAIN_ID)        &
                              ,npe_print          =NPE_PRINT                   &
                              ,print_timing       =PRINT_TIMING )
!
!-----------------------------------------------------------------------
!
            IF(ESMF_ClockIsStopTime(CLOCK_FILTER(MY_DOMAIN_ID),rc=RC))THEN
              GENERATION_FINISHED(N)=.TRUE.                                !<-- Generation N's filter has finished
              IF(I_AM_LEAD_FCST_TASK)THEN
                WRITE(0,*)' Completed filter method ',filter_method
                WRITE(0,*)' Now reset filter method to 0.'
              ENDIF
!
            ENDIF
!
!-----------------------------------------------------------------------
!***  All tasks that are finished on all generations may leave.
!-----------------------------------------------------------------------
!
            IF(ALL(GENERATION_FINISHED,NUM_GENS))THEN                        !<-- If true, all of this task's domains are finished
              ALL_FORECASTS_COMPLETE=.TRUE.
              EXIT gens_f3_4
            ENDIF
!
!-----------------------------------------------------------------------
!
          ENDDO gens_f3_4
!
!-----------------------------------------------------------------------
!
        ENDDO tdfi_forward
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      ENDIF  method_block
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Remove working objects we are finished with.
!-----------------------------------------------------------------------
!
      DEALLOCATE(CLOCK_FILTER,stat=RC)
      DEALLOCATE(HALFDFITIME,stat=RC)
      DEALLOCATE(HALFDFIINTVAL,stat=RC)
      DEALLOCATE(NDFISTEP,stat=RC)
!
!-----------------------------------------------------------------------
!***  Reset the completion flags for integration of the free forecast.
!-----------------------------------------------------------------------
!
      IF(NUM_GENS==1)THEN
        GENERATION_FINISHED(1)=.FALSE.
!
      ELSE
        DO N=1,NUM_GENS
          IF(MY_DOMAINS_IN_GENS(N)>0)THEN
            GENERATION_FINISHED(N)=.FALSE.
          ELSE
            GENERATION_FINISHED(N)=.TRUE.                                  !<-- Task not in this generation; consider it finished.
          ENDIF
        ENDDO
!
      ENDIF
!
      ALL_FORECASTS_COMPLETE=.FALSE.
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE RUN_DIGITAL_FILTER_NMM
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE NMM_RUN
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE NMM_FINALIZE(NMM_GRID_COMP                             &
                             ,IMP_STATE                                 &
                             ,EXP_STATE                                 &
                             ,CLOCK_NMM                                 &
                             ,RC_FINALIZE)
!
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(ESMF_GridComp) :: NMM_GRID_COMP                                 !<-- The NMM component
!
      TYPE(ESMF_State) :: IMP_STATE                                     &  !<-- The NMM import state
                         ,EXP_STATE                                        !<-- The NMM export state
!
      TYPE(ESMF_Clock) :: CLOCK_NMM                                        !<-- The NMM component's ESMF Clock

      INTEGER,INTENT(OUT) :: RC_FINALIZE                                   !<-- Error return code
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER :: I,J,N
      INTEGER(kind=KINT) :: RC
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------

      RC         =ESMF_SUCCESS
      RC_FINALIZE=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!
      IF(RC_FINALIZE==ESMF_SUCCESS)THEN
!       WRITE(0,*)' NMM_FINALIZE succeeded'
      ELSE
        WRITE(0,*)' NMM_FINALIZE failed  RC_FINALIZE=',RC_FINALIZE
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE NMM_FINALIZE
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      RECURSIVE SUBROUTINE CALL_DOMAIN_INITIALIZE(ID_DOMAIN,CLOCK_NMM)
!
!-----------------------------------------------------------------------
!***  This routine calls DOMAIN_INITIALIZE for all DOMAIN components.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      INTEGER,INTENT(IN) :: ID_DOMAIN                                      !<-- ID of the DOMAIN Component to initialize
!
      TYPE(ESMF_Clock),DIMENSION(1:NUM_DOMAINS_TOTAL),INTENT(INOUT) :: CLOCK_NMM !<-- The NMM ESMF Clock
!
!-----------------------------------------------------------------------
!***  Local Variables
!-----------------------------------------------------------------------
!
      INTEGER :: ID_CHILD,IRTN,N,N_CHILDREN
      INTEGER :: RC,RC_CALL_INIT
!
      CHARACTER(2)  :: INT_TO_CHAR
      CHARACTER(6)  :: FMT
!
      integer :: i_par_sta,j_par_sta,next_move_timestep
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC          =ESMF_SUCCESS
      RC_CALL_INIT=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!
      FMT='(I2.2)'
      WRITE(INT_TO_CHAR,FMT)ID_DOMAIN
!
!-----------------------------------------------------------------------
!***  Initialize the DOMAIN component with the ID of ID_DOMAIN.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Initialize DOMAIN Component "//INT_TO_CHAR
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompInitialize(gridcomp   =nmm_int_state%DOMAIN_GRID_COMP(ID_DOMAIN)  &  !<-- The DOMAIN component
                                  ,importState=nmm_int_state%IMP_STATE_DOMAIN(ID_DOMAIN)  &  !<-- The DOMAIN import state
                                  ,exportState=nmm_int_state%EXP_STATE_DOMAIN(ID_DOMAIN)  &  !<-- The DOMAIN export state
                                  ,clock      =CLOCK_NMM(ID_DOMAIN)                       &  !<-- The DOMAIN clock
                                  ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!     CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CALL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  If the domain being initialized has children, force those 
!***  children to wait.  The parent might be generating input for
!***  the children therefore children should not be trying to read
!***  the input in their own initialize steps prematurely.
!-----------------------------------------------------------------------
!
      CALL MPI_BARRIER(COMM_GLOBAL,IRTN)
!
!-----------------------------------------------------------------------
!***  If there are children, initialize them.
!-----------------------------------------------------------------------
!
      N_CHILDREN=nmm_int_state%NUM_CHILDREN(ID_DOMAIN)
!
      IF(N_CHILDREN>0)THEN                                                 !<-- Does the current DOMAIN have any children?
        DO N=1,N_CHILDREN                                                  !<-- If so, loop through the children to Initialize them
          ID_CHILD=ID_CHILDREN(N,ID_DOMAIN)
          CALL CALL_DOMAIN_INITIALIZE(ID_CHILD,CLOCK_NMM)
        ENDDO
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE CALL_DOMAIN_INITIALIZE
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      SUBROUTINE BOUNDARY_DATA_TO_DOMAIN(EXP_STATE_CPL                  &
                                        ,IMP_STATE_DOMAIN )
!
!-----------------------------------------------------------------------
!***  This routine moves new boundary data for nested domains from the
!***  export state of the Parent-Child coupler to the import state of
!***  the NMM nests' DOMAIN components.
!-----------------------------------------------------------------------
!
      TYPE(ESMF_State),INTENT(INOUT) :: EXP_STATE_CPL                      !<-- The Parent-Child Coupler's export state

      TYPE(ESMF_State),INTENT(OUT)   :: IMP_STATE_DOMAIN                   !<-- The nests' DOMAIN import state
!
!-----------------------------------------------------------------------
!***  Local variables
!-----------------------------------------------------------------------
!
      TYPE SIDES_1D_REAL
        REAL,DIMENSION(:),ALLOCATABLE :: SOUTH
        REAL,DIMENSION(:),ALLOCATABLE :: NORTH
        REAL,DIMENSION(:),ALLOCATABLE :: WEST
        REAL,DIMENSION(:),ALLOCATABLE :: EAST
      END TYPE SIDES_1D_REAL
!
      INTEGER :: ISTAT,KOUNT,RC,RC_BND_MV
!
      TYPE(SIDES_1D_REAL),SAVE :: BOUNDARY_H                            &
                                 ,BOUNDARY_V
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Check each side of the child boundary.  If data is present from
!***  that side in the Parent-Child coupler export state then move it
!***  to the DOMAIN component's import state.
!-----------------------------------------------------------------------
!
!-------------
!***  South H
!-------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Check Parent-Child Cpl Export State for South H Data"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state    =EXP_STATE_CPL                    &   !<-- Look at the Parent-Child Coupler's export state
                            ,name     ='SOUTH_H'                        &   !<-- Is this name present?
                            ,itemCount=KOUNT                            &   !<-- How many items present?
                            ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      south_h: IF(KOUNT>0.AND.RC==ESMF_SUCCESS)THEN                         !<-- True => South boundary H point data is present
!
        IF(.NOT.ALLOCATED(BOUNDARY_H%SOUTH))THEN
          ALLOCATE(BOUNDARY_H%SOUTH(1:KOUNT),stat=ISTAT)
          IF(ISTAT/=0)WRITE(0,*)' Failed to allocate BOUNDARY_H%SOUTH stat=',ISTAT
        ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract South H Data from Parent-Child Cpl Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

        CALL ESMF_AttributeGet(state    =EXP_STATE_CPL                  &   !<-- Extract data from Parent-Child Coupler's export state
                              ,name     ='SOUTH_H'                      &   !<-- The name of the data
                              ,itemCount=KOUNT                          &   !<-- The data has this many items
                              ,valueList=BOUNDARY_H%SOUTH               &   !<-- The new combined boundary data
                              ,rc=RC )

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert South H Data into Nest DOMAIN Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

        CALL ESMF_AttributeSet(state    =IMP_STATE_DOMAIN               &   !<-- Insert data into nest's DOMAIN import state
                              ,name     ='SOUTH_H'                      &   !<-- The name of the data
                              ,itemCount=KOUNT                          &   !<-- The data has this many items
                              ,valueList=BOUNDARY_H%SOUTH               &   !<-- The new combined boundary data
                              ,rc=RC )

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF south_h
!
!-------------
!***  South V
!-------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Check Parent-Child Cpl Export State for South V Data"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

      CALL ESMF_AttributeGet(state    =EXP_STATE_CPL                    &   !<-- Look at the Parent-Child Coupler's export state
                            ,name     ='SOUTH_V'                        &   !<-- Is this name present?
                            ,itemCount=KOUNT                            &   !<-- How many items present?
                            ,rc       =RC )

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      south_v: IF(KOUNT>0.AND.RC==ESMF_SUCCESS)THEN                         !<-- True => South boundary V point data is present
!
        IF(.NOT.ALLOCATED(BOUNDARY_V%SOUTH))THEN
          ALLOCATE(BOUNDARY_V%SOUTH(1:KOUNT),stat=ISTAT)
          IF(ISTAT/=0)WRITE(0,*)' Failed to allocate BOUNDARY_V%SOUTH stat=',ISTAT
        ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract South V Data from Parent-Child Cpl Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

        CALL ESMF_AttributeGet(state    =EXP_STATE_CPL                  &   !<-- Extract data from Parent-Child Coupler's export state
                              ,name     ='SOUTH_V'                      &   !<-- The name of the data
                              ,itemCount=KOUNT                          &   !<-- The data has this many items
                              ,valueList=BOUNDARY_V%SOUTH               &   !<-- The new combined boundary data
                              ,rc=RC )

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert South V Data into Nest DOMAIN Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

        CALL ESMF_AttributeSet(state    =IMP_STATE_DOMAIN               &   !<-- Insert data into nest's DOMAIN import state
                              ,name     ='SOUTH_V'                      &   !<-- The name of the data
                              ,itemCount=KOUNT                          &   !<-- The data has this many items
                              ,valueList=BOUNDARY_V%SOUTH               &   !<-- The new combined boundary data
                              ,rc=RC )

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF south_v
!
!-------------
!***  North H
!-------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Check Parent-Child Cpl Export State for North H Data"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

      CALL ESMF_AttributeGet(state    =EXP_STATE_CPL                    &   !<-- Look at the Parent-Child Coupler's export state
                            ,name     ='NORTH_H'                        &   !<-- Is this name present?
                            ,itemCount=KOUNT                            &   !<-- How many items present?
                            ,rc       =RC )

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      north_h: IF(KOUNT>0.AND.RC==ESMF_SUCCESS)THEN                         !<-- True => North boundary H point data is present
!
        IF(.NOT.ALLOCATED(BOUNDARY_H%NORTH))THEN
          ALLOCATE(BOUNDARY_H%NORTH(1:KOUNT),stat=ISTAT)
          IF(ISTAT/=0)WRITE(0,*)' Failed to allocate BOUNDARY_H%NORTH stat=',ISTAT
        ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract North H Data from Parent-Child Cpl Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

        CALL ESMF_AttributeGet(state    =EXP_STATE_CPL                  &   !<-- Extract data from Parent-Child Coupler's export state
                              ,name     ='NORTH_H'                      &   !<-- The name of the data
                              ,itemCount=KOUNT                          &   !<-- The data has this many items
                              ,valueList=BOUNDARY_H%NORTH               &   !<-- The new combined boundary data
                              ,rc=RC )

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert North H Data into Nest DOMAIN Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

        CALL ESMF_AttributeSet(state    =IMP_STATE_DOMAIN               &   !<-- Insert data into nest's DOMAIN import state
                              ,name     ='NORTH_H'                      &   !<-- The name of the data
                              ,itemCount=KOUNT                          &   !<-- The data has this many items
                              ,valueList=BOUNDARY_H%NORTH               &   !<-- The new combined boundary data
                              ,rc=RC )

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF north_h
!
!-------------
!***  North V
!-------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Check Parent-Child Cpl Export State for North V Data"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

      CALL ESMF_AttributeGet(state    =EXP_STATE_CPL                    &   !<-- Look at the Parent-Child Coupler's export state
                            ,name     ='NORTH_V'                        &   !<-- Is this name present?
                            ,itemCount=KOUNT                            &   !<-- How many items present?
                            ,rc       =RC )

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      north_v: IF(KOUNT>0.AND.RC==ESMF_SUCCESS)THEN                         !<-- True => North boundary V point data is present
!
        IF(.NOT.ALLOCATED(BOUNDARY_V%NORTH))THEN
          ALLOCATE(BOUNDARY_V%NORTH(1:KOUNT),stat=ISTAT)
          IF(ISTAT/=0)WRITE(0,*)' Failed to allocate BOUNDARY_V%NORTH stat=',ISTAT
        ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract North V Data from Parent-Child Cpl Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

        CALL ESMF_AttributeGet(state    =EXP_STATE_CPL                  &   !<-- Extract data from Parent-Child Coupler's export state
                              ,name     ='NORTH_V'                      &   !<-- The name of the data
                              ,itemCount=KOUNT                          &   !<-- The data has this many items
                              ,valueList=BOUNDARY_V%NORTH               &   !<-- The new combined boundary data
                              ,rc=RC )

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert North V Data into Nest DOMAIN Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

        CALL ESMF_AttributeSet(state    =IMP_STATE_DOMAIN               &   !<-- Insert data into nest's DOMAIN import state
                              ,name     ='NORTH_V'                      &   !<-- The name of the data
                              ,itemCount=KOUNT                          &   !<-- The data has this many items
                              ,valueList=BOUNDARY_V%NORTH               &   !<-- The new combined boundary data
                              ,rc=RC )

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF north_v
!
!------------
!***  West H
!------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Check Parent-Child Cpl Export State for West H Data"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

      CALL ESMF_AttributeGet(state    =EXP_STATE_CPL                    &   !<-- Look at the Parent-Child Coupler's export state
                            ,name     ='WEST_H'                         &   !<-- Is this name present?
                            ,itemCount=KOUNT                            &   !<-- How many items present?
                            ,rc       =RC )

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      west_h: IF(KOUNT>0.AND.RC==ESMF_SUCCESS)THEN                          !<-- True => West boundary H point data is present
!
        IF(.NOT.ALLOCATED(BOUNDARY_H%WEST))THEN
          ALLOCATE(BOUNDARY_H%WEST(1:KOUNT),stat=ISTAT)
          IF(ISTAT/=0)WRITE(0,*)' Failed to allocate BOUNDARY_H%WEST stat=',ISTAT
        ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract West H Data from Parent-Child Cpl Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

        CALL ESMF_AttributeGet(state    =EXP_STATE_CPL                  &   !<-- Extract data from Parent-Child Coupler's export state
                              ,name     ='WEST_H'                       &   !<-- The name of the data
                              ,itemCount=KOUNT                          &   !<-- The data has this many items
                              ,valueList=BOUNDARY_H%WEST                &   !<-- The new combined boundary data
                              ,rc=RC )

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert West H Data into Nest DOMAIN Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

        CALL ESMF_AttributeSet(state    =IMP_STATE_DOMAIN               &   !<-- Insert data into nest's DOMAIN import state
                              ,name     ='WEST_H'                       &   !<-- The name of the data
                              ,itemCount=KOUNT                          &   !<-- The data has this many items
                              ,valueList=BOUNDARY_H%WEST                &   !<-- The new combined boundary data
                              ,rc=RC )

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF west_h
!
!------------
!***  West V
!------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Check Parent-Child Cpl Export State for West V Data"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

      CALL ESMF_AttributeGet(state    =EXP_STATE_CPL                    &   !<-- Look at the Parent-Child Coupler's export state
                            ,name     ='WEST_V'                         &   !<-- Is this name present?
                            ,itemCount=KOUNT                            &   !<-- How many items present?
                            ,rc       =RC )

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      west_v: IF(KOUNT>0.AND.RC==ESMF_SUCCESS)THEN                          !<-- True => West boundary V point data is present
!
        IF(.NOT.ALLOCATED(BOUNDARY_V%WEST))THEN
          ALLOCATE(BOUNDARY_V%WEST(1:KOUNT),stat=ISTAT)
          IF(ISTAT/=0)WRITE(0,*)' Failed to allocate BOUNDARY_V%WEST stat=',ISTAT
        ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract West V Data from Parent-Child Cpl Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

        CALL ESMF_AttributeGet(state    =EXP_STATE_CPL                  &   !<-- Extract data from Parent-Child Coupler's export state
                              ,name     ='WEST_V'                       &   !<-- The name of the data
                              ,itemCount=KOUNT                          &   !<-- The data has this many items
                              ,valueList=BOUNDARY_V%WEST                &   !<-- The new combined boundary data
                              ,rc=RC )

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert West V Data into Nest DOMAIN Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

        CALL ESMF_AttributeSet(state    =IMP_STATE_DOMAIN               &   !<-- Insert data into nest's DOMAIN import state
                              ,name     ='WEST_V'                       &   !<-- The name of the data
                              ,itemCount=KOUNT                          &   !<-- The data has this many items
                              ,valueList=BOUNDARY_V%NORTH               &   !<-- The new combined boundary data
                              ,rc=RC )

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF west_v
!
!------------
!***  East H
!------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Check Parent-Child Cpl Export State for East H Data"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

      CALL ESMF_AttributeGet(state    =EXP_STATE_CPL                    &   !<-- Look at the Parent-Child Coupler's export state
                            ,name     ='EAST_H'                         &   !<-- Is this name present?
                            ,itemCount=KOUNT                            &   !<-- How many items present?
                            ,rc       =RC )

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      east_h: IF(KOUNT>0.AND.RC==ESMF_SUCCESS)THEN                          !<-- True => East boundary H point data is present
!
        IF(.NOT.ALLOCATED(BOUNDARY_H%EAST))THEN
          ALLOCATE(BOUNDARY_H%EAST(1:KOUNT),stat=ISTAT)
          IF(ISTAT/=0)WRITE(0,*)' Failed to allocate BOUNDARY_H%EAST stat=',ISTAT
        ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract East H Data from Parent-Child Cpl Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

        CALL ESMF_AttributeGet(state    =EXP_STATE_CPL                  &   !<-- Extract data from Parent-Child Coupler's export state
                              ,name     ='EAST_H'                       &   !<-- The name of the data
                              ,itemCount=KOUNT                          &   !<-- The data has this many items
                              ,valueList=BOUNDARY_H%EAST                &   !<-- The new combined boundary data
                              ,rc=RC )

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert East H Data into Nest DOMAIN Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

        CALL ESMF_AttributeSet(state    =IMP_STATE_DOMAIN               &   !<-- Insert data into nest's DOMAIN import state
                              ,name     ='EAST_H'                       &   !<-- The name of the data
                              ,itemCount=KOUNT                          &   !<-- The data has this many items
                              ,valueList=BOUNDARY_H%EAST                &   !<-- The new combined boundary data
                              ,rc=RC )

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF east_h
!
!------------
!***  East V
!------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Check Parent-Child Cpl Export State for East V Data"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

      CALL ESMF_AttributeGet(state    =EXP_STATE_CPL                    &   !<-- Look at the Parent-Child Coupler's export state
                            ,name     ='EAST_V'                         &   !<-- Is this name present?
                            ,itemCount=KOUNT                            &   !<-- How many items present?
                            ,rc       =RC )

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      east_v: IF(KOUNT>0.AND.RC==ESMF_SUCCESS)THEN                          !<-- True => East boundary V point data is present
!
        IF(.NOT.ALLOCATED(BOUNDARY_V%EAST))THEN
          ALLOCATE(BOUNDARY_V%EAST(1:KOUNT),stat=ISTAT)
          IF(ISTAT/=0)WRITE(0,*)' Failed to allocate BOUNDARY_V%EAST stat=',ISTAT
        ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract East V Data from Parent-Child Cpl Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

        CALL ESMF_AttributeGet(state    =EXP_STATE_CPL                  &   !<-- Extract data from Parent-Child Coupler's export state
                              ,name     ='EAST_V'                       &   !<-- The name of the data
                              ,itemCount=KOUNT                          &   !<-- The data has this many items
                              ,valueList=BOUNDARY_V%EAST                &   !<-- The new combined boundary data
                              ,rc=RC )

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert East V Data into Nest DOMAIN Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

        CALL ESMF_AttributeSet(state    =IMP_STATE_DOMAIN               &   !<-- Insert data into nest's DOMAIN import state
                              ,name     ='EAST_V'                       &   !<-- The name of the data
                              ,itemCount=KOUNT                          &   !<-- The data has this many items
                              ,valueList=BOUNDARY_V%EAST                &   !<-- The new combined boundary data
                              ,rc=RC )

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF east_v
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE BOUNDARY_DATA_TO_DOMAIN

!-----------------------------------------------------------------------
      SUBROUTINE NMM_CLOCKPRINT(CLOCK, msg, RC)
!-----------------------------------------------------------------------
        type(ESMF_CLOCK), intent(inout) :: CLOCK
        character(len=*), intent(in)    :: msg
        integer, intent(out)            :: RC

        RC = ESMF_SUCCESS
        call ESMF_ClockPrint(CLOCK, &
          preString=trim(msg)// ' CurrentTime: ', unit=nuopcMsg)
        call ESMF_LogWrite(nuopcMsg, ESMF_LOGMSG_INFO)
        call ESMF_ClockPrint(CLOCK, &
          preString=trim(msg)// ' StartTime: ', unit=nuopcMsg)
        call ESMF_LogWrite(nuopcMsg, ESMF_LOGMSG_INFO)
        call ESMF_ClockPrint(CLOCK, &
          preString=trim(msg)// ' StopTime: ', unit=nuopcMsg)
        call ESMF_LogWrite(nuopcMsg, ESMF_LOGMSG_INFO)

!-----------------------------------------------------------------------
      END SUBROUTINE
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
  subroutine AdvertiseFields(state, nfields, field_defs, rc)
!
!-----------------------------------------------------------------------

    type(ESMF_State), intent(inout)             :: state
    integer,intent(in)                          :: nfields
    type(fld_list_type), intent(inout)          :: field_defs(:)
    integer, intent(inout)                      :: rc

    integer                                     :: i
    character(len=*),parameter  :: subname='(nmmb_cap:AdvertiseFields)'

!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    do i = 1, nfields

      call ESMF_LogWrite('NMMB Advertise: '//trim(field_defs(i)%stdname), ESMF_LOGMSG_INFO, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      call NUOPC_Advertise(state, &
        standardName=field_defs(i)%stdname, &
        name=field_defs(i)%shortname, &
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

    enddo

!-----------------------------------------------------------------------
!
  end subroutine AdvertiseFields
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
  subroutine RealizeFields(state, grid, nfields, field_defs, tag, rc)
!
!-----------------------------------------------------------------------

    type(ESMF_State), intent(inout)             :: state
    type(ESMF_Grid), intent(in)                 :: grid
    integer, intent(in)                         :: nfields
    type(fld_list_type), intent(inout)          :: field_defs(:)
    character(len=*), intent(in)                :: tag
    integer, intent(inout)                      :: rc

    integer                                     :: i
    type(ESMF_Field)                            :: field
    integer                                     :: npet, nx, ny, pet, &
        elb(2), eub(2), clb(2), cub(2), tlb(2), tub(2)
    type(ESMF_VM)                               :: vm
    character(len=*),parameter  :: subname='(nmmb_cap:RealizeFields)'
 
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
 
    rc = ESMF_SUCCESS

      !call ESMF_VMGetCurrent(vm, rc=rc)
      !if (rc /= ESMF_SUCCESS) call ESMF_Finalize()

      !call ESMF_VMGet(vm, petcount=npet, localPet=pet, rc=rc)
      !if (rc /= ESMF_SUCCESS) call ESMF_Finalize()

      !call ESMF_GridGet(grid, exclusiveLBound=elb, exclusiveUBound=eub, &
      !                        computationalLBound=clb, computationalUBound=cub, &
      !                        totalLBound=tlb, totalUBound=tub, rc=rc)
      !if (rc /= ESMF_SUCCESS) call ESMF_Finalize()

      !write(info, *) pet, 'exc', elb, eub, 'comp', clb, cub, 'total', tlb, tub
      !call ESMF_LogWrite(subname // tag // " Grid "// info, &
      !  ESMF_LOGMSG_INFO, &
      !  line=__LINE__, &
      !  file=__FILE__, &
      !  rc=rc)

    do i = 1, nfields

      if (field_defs(i)%assoc) then
        write(info, *) subname, tag, ' Field ', field_defs(i)%shortname, ':', &
          lbound(field_defs(i)%farrayPtr,1), ubound(field_defs(i)%farrayPtr,1), &
          lbound(field_defs(i)%farrayPtr,2), ubound(field_defs(i)%farrayPtr,2), &
          lbound(field_defs(i)%farrayPtr,3), ubound(field_defs(i)%farrayPtr,3)
        call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=rc)
        field = ESMF_FieldCreate(grid=grid, &
          farray=field_defs(i)%farrayPtr, indexflag=ESMF_INDEX_DELOCAL, &
          name=field_defs(i)%shortname, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      else
        field = ESMF_FieldCreate(grid, ESMF_TYPEKIND_R8, indexflag=ESMF_INDEX_DELOCAL, &
          name=field_defs(i)%shortname, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

      if (NUOPC_IsConnected(state, fieldName=field_defs(i)%shortname)) then
        call NUOPC_Realize(state, field=field, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        call ESMF_LogWrite(subname // tag // " Field "// field_defs(i)%stdname // " is connected.", &
          ESMF_LOGMSG_INFO, &
          line=__LINE__, &
          file=__FILE__, &
          rc=rc)
!        call ESMF_FieldPrint(field=field, rc=rc)
!        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!          line=__LINE__, &
!          file=__FILE__)) &
!          return  ! bail out
      else
        call ESMF_LogWrite(subname // tag // " Field "// field_defs(i)%stdname // " is not connected.", &
          ESMF_LOGMSG_INFO, &
          line=__LINE__, &
          file=__FILE__, &
          rc=rc)
        ! TODO: Initialize the value in the pointer to 0 after proper restart is setup
        !if(associated(field_defs(i)%farrayPtr) ) field_defs(i)%farrayPtr = 0.0
        ! remove a not connected Field from State
        call ESMF_StateRemove(state, (/field_defs(i)%shortname/), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

    enddo

!-----------------------------------------------------------------------
!
  end subroutine RealizeFields
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
  subroutine fld_list_add(num, fldlist, stdname, transferOffer, data, shortname)
!
!-----------------------------------------------------------------------
    ! ----------------------------------------------
    ! Set up a list of field information
    ! ----------------------------------------------
    integer,             intent(inout)  :: num
    type(fld_list_type), intent(inout)  :: fldlist(:)
    character(len=*),    intent(in)     :: stdname
    character(len=*),    intent(in)     :: transferOffer
    real(ESMF_KIND_R8), dimension(:,:,:), optional, target :: data
    character(len=*),    intent(in),optional :: shortname

    ! local variables
    integer :: rc
    character(len=*), parameter :: subname='(cice_cap:fld_list_add)'

!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
    ! fill in the new entry

    num = num + 1
    if (num > fldsMax) then
      call ESMF_LogWrite(trim(subname)//": ERROR num gt fldsMax "//trim(stdname), &
        ESMF_LOGMSG_ERROR, line=__LINE__, file=__FILE__, rc=rc)
      return
    endif

    fldlist(num)%stdname        = trim(stdname)
    if (present(shortname)) then
       fldlist(num)%shortname   = trim(shortname)
    else
       fldlist(num)%shortname   = trim(stdname)
    endif
    fldlist(num)%transferOffer  = trim(transferOffer)
    if (present(data)) then
      fldlist(num)%assoc        = .true.
      fldlist(num)%farrayPtr    => data
    else
      fldlist(num)%assoc        = .false.
    endif

  end subroutine fld_list_add
!  
!-----------------------------------------------------------------------
!
      END MODULE module_NMM_GRID_COMP
!
!-----------------------------------------------------------------------
