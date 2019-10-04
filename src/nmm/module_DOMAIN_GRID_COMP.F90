!-----------------------------------------------------------------------
!
      MODULE module_DOMAIN_GRID_COMP
!
!-----------------------------------------------------------------------
!
!***  This is the Domain gridded component module.
!***  It will set up Solver subcomponents
!***  and run their Initialize, Run, and Finalize routines.
!
!-----------------------------------------------------------------------
!
! PROGRAM HISTORY LOG:
!   2007-       Black - Modified from Wei-yu's version
!   2007-11-20  Black/Henderson - Created an ESMF Clock for the
!                                 ATM Component independent of
!                                 the Main Clock.
!   2007-12-11  Black - Generalized for easier use by any dynamics core.
!   2008-08     Colon - Added conditional checks multiple dynamics cores.
!   2008-10-14  Vasic - Added restart Alarm.
!   2009-08-03  Black - Merging with nesting.
!   2009-08-12  Black - Fixed logic for Physics export when direction of
!                       integration switches from backward to forward.
!   2009-10-05  Wang  - Added GFS ensemble member name and output data at
!                       every nsout timesteps.
!   2010-03-24  Black - Converted to Domain component for NMM-B only.
!   2010-11-03  Pyle  - Revised for digital filters.
!   2010-12-16  Pyle  - Change to nemsio library
!   2011-02     Yang  - Updated to use both the ESMF 4.0.0rp2 library,
!                       ESMF 5 series library and the the
!                       ESMF 3.1.0rp2 library.
!   2011-05-11  Yang  - Modified for using the ESMF 5.2.0r_beta_snapshot_07.
!   2011-07-16  Black - Moving nest capability.
!   2012-02-08  Yang  - Modified for using the ESMF 5.2.0rp1 library.
!   2012-07-20 Black  - Revised for generational task usage.
!   2016-07    Black  - Modifications for ocean coupling

!
! USAGE: Domain gridded component parts called from subroutines within
!        module_NMM_GRID_COMP.F90.
!
!-----------------------------------------------------------------------
!
      USE MPI
      USE ESMF
      USE netcdf
!
      USE MODULE_KINDS
!
      USE MODULE_CONSTANTS,ONLY : A,CP,G,P608,PI,TWOM 
!
      USE MODULE_DERIVED_TYPES,ONLY: BC_H_ALL,BC_V_ALL,MIXED_DATA
!
      USE MODULE_DOMAIN_INTERNAL_STATE,ONLY: DOMAIN_INTERNAL_STATE      &
                                            ,WRAP_DOMAIN_INTERNAL_STATE

      USE MODULE_SOLVER_GRID_COMP,ONLY: SOLVER_REGISTER
!
      USE MODULE_DM_PARALLEL,ONLY : IDS,IDE,JDS,JDE                     &
                                   ,IMS,IME,JMS,JME                     &
                                   ,ITS,ITE,JTS,JTE                     &
                                   ,IHALO,JHALO                         &
                                   ,MPI_COMM_COMP                       &
                                   ,MY_NEB
!
      USE MODULE_CONTROL,ONLY: BOUNDARY_INIT                            &
                              ,GRID_CONSTS                              &
                              ,NUM_DOMAINS_MAX                          &
                              ,TIMEF
!
      USE MODULE_DOMAIN_NUOPC_SET,ONLY: DOMAIN_DESCRIPTORS              &
                                       ,DUMP_DOMAIN_DESCRIPTOR          &
                                       ,NMMB_CreateDomainFields         &
                                       ,NMMB_CreateRouteHandle          &
                                       ,NMMB_GridCreate                 &
                                       ,NMMB_GridUpdate                 &
                                       ,ROTANGLE_CELLAREA_SEAMASK
!
      USE MODULE_DIAGNOSE,ONLY: FIELD_STATS
      USE NEMSIO_MODULE
! 
      USE MODULE_ERROR_MSG,ONLY: ERR_MSG,MESSAGE_CHECK
!
      USE MODULE_VARS,ONLY: FIND_VAR_INDX                               &
                           ,TKR_I0D,TKR_I1D                             &
                           ,TKR_R0D,TKR_R1D                             &
                           ,VAR
!
      USE MODULE_SOLVER_INTERNAL_STATE,ONLY: SOLVER_INTERNAL_STATE      &
                                            ,WRAP_SOLVER_INT_STATE
!
      USE MODULE_NESTING,ONLY: CHECK                                    &
                              ,LATLON_TO_IJ                             &
                              ,INTERNAL_DATA_TO_DOMAIN                  &
                              ,PARENT_TO_CHILD_INIT_NMM                 &
                              ,SUFFIX_MOVE                              &
                              ,SUFFIX_NESTBC                            &
                              ,SUFFIX_TWOWAY
!
      USE MODULE_OUTPUT,ONLY: POINT_OUTPUT
!
      USE MODULE_CLOCKTIMES,ONLY : TIMERS
!
!-----------------------------------------------------------------------
!***  List other modules with non-generic routines used by DOMAIN.
!-----------------------------------------------------------------------
!
      USE MODULE_WRITE_ROUTINES ,ONLY: WRITE_INIT,WRITE_ASYNC             !<-- These are routines used only when asynchronous
      USE MODULE_WRITE_GRID_COMP,ONLY: WRITE_SETUP                      & !    quilting is specified by the user in the
                                      ,WRITE_DESTROY                      !    configure file for history output.
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
! 
      PRIVATE
!
      PUBLIC :: DOMAIN_REGISTER,MY_DOMAIN_ID
!
!-----------------------------------------------------------------------
!
      INTEGER(kind=KINT),SAVE :: N8=8
!
      INTEGER(kind=KINT) :: INPES                                       &  !<-- Forecast tasks in I on this domain's grid
                           ,JNPES                                       &  !<-- Forecast tasks in J on this domain's grid
                           ,LM                                          &  !<-- # of model layers
                           ,MYPE                                        &  !<-- Each MPI task ID
                           ,NPE_PRINT                                   &
                           ,NLAYRS                                      &  !<-- Number of model layers
                           ,NTIMESTEP                                   &  !<-- Integration timestep
                           ,NUM_TRACERS_CHEM                            &  !<-- Number of chemistry tracer variables
                           ,WRITE_GROUP_READY_TO_GO                        !<-- The write group to use
!
      INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE,SAVE :: COMM_FCST_TASKS  !<-- Hold the intracommunicator for each domain's fcst tasks.
!
      REAL(kind=DOUBLE) :: D2R,D_ONE,D_180,PI_LOC
!
      REAL(kind=KFPT),SAVE :: ATM_OCN_CPL_INT=-999.                        !<-- The atmos-ocean coupling interval (sec)
!
      REAL(kind=KFPT),SAVE :: SBD_1,TPH0D_1,TLM0D_1,WBD_1                  !<-- SW corner & center (degrees) of upper parent
!
      LOGICAL(kind=KLOG) :: QUILTING                                    &  !<-- Is asynchronous quilting specified?
                           ,WRITE_LAST_RESTART                          &  !<-- Write last restart file?
                           ,RESTARTED_RUN                                  !<-- Restarted run logical flag
!
      TYPE(ESMF_VM),SAVE :: VM,VM_LOCAL                                    !<-- The ESMF virtual machine.
!
      TYPE(DOMAIN_INTERNAL_STATE),POINTER :: DOMAIN_INT_STATE              !<-- Any given Domain internal state
!
      TYPE(DOMAIN_INTERNAL_STATE),DIMENSION(:),POINTER,SAVE ::          &
                                                 DOMAIN_INT_STATE_ALL      !<-- The NMM Domain internal state pointer
!
      TYPE(WRAP_DOMAIN_INTERNAL_STATE),SAVE :: WRAP                        !<-- The F90 wrap of the NMM Domain internal state
!
      TYPE(SOLVER_INTERNAL_STATE),POINTER :: SOLVER_INT_STATE
!
      TYPE(ESMF_Time),SAVE :: DFITIME                                   &
                             ,HALFDFITIME 
!
      TYPE(ESMF_TimeInterval),SAVE :: HALFDFIINTVAL                        !<-- The ESMF time interval for filtering
!
      TYPE(ESMF_Config),SAVE :: CF_ATMOS                                   !<-- The atmos configure file
      TYPE(ESMF_Config),DIMENSION(99),SAVE :: CF                           !<-- The configure objects for all NMM domains
!
      LOGICAL(kind=KLOG) :: PHYSICS_ON                                     !<-- Is physics active?
!
!-----------------------------------------------------------------------
!
!---------------------
!***  For NMM Nesting
!---------------------
!
      INTEGER(kind=KINT),SAVE :: COMM_FULL_DOMAIN                       &  !<-- Communicator for ALL tasks on domain to be split
                                ,COMM_MY_DOMAIN                         &  !<-- Each domain's local intracommunicator
                                ,IM_1                                   &  !<-- I dimension of uppermost parent domain
                                ,I_SHIFT_CHILD                          &  !<-- Nest's shift in I in the nest's grid space
                                ,JM_1                                   &  !<-- J dimension of uppermost parent domain
                                ,J_SHIFT_CHILD                          &  !<-- Nest's shift in J in the nest's grid space
                                ,MY_DOMAIN_ID                           &  !<-- Domain IDs; begin with uppermost parent=1
                                ,NROWS_P_UPD_E                          &  !<-- # of footprint E bndry rows using parent updates
                                ,NROWS_P_UPD_N                          &  !<-- # of footprint N bndry rows using parent updates
                                ,NROWS_P_UPD_S                          &  !<-- # of footprint S bndry rows using parent updates
                                ,NROWS_P_UPD_W                          &  !<-- # of footprint W bndry rows using parent updates
                                ,NUM_CHILDREN                           &  !<-- Number of (1st generation) children within a domain
                                ,NUM_FIELDS_MOVE_2D_H_I                 &  !<-- # of 2-D integer H variables updated for moving nests
                                ,NUM_FIELDS_MOVE_2D_H_R                 &  !<-- # of 2-D real H variables updated for moving nests
                                ,NUM_FIELDS_MOVE_3D_H                   &  !<-- Number of 3-D H variables updated for moving nests
                                ,NUM_LEVELS_MOVE_3D_H                   &  !<-- Number of 2-D levels in all 3-D H update variables
                                ,NUM_FIELDS_MOVE_2D_V                   &  !<-- Number of 2-D V variables updated for moving nests
                                ,NUM_FIELDS_MOVE_3D_V                   &  !<-- Number of 3-D V variables updated for moving nests
                                ,NUM_LEVELS_MOVE_3D_V                   &  !<-- Number of 2-D levels in all 3-D V update variables
                                ,PARENT_CHILD_SPACE_RATIO               &  !<-- Ratio of parent's space increment to the nest's
                                ,PARENT_CHILD_TIME_RATIO                   !<-- # of child timesteps per parent timestep
!
      INTEGER(kind=KINT),PARAMETER :: N_PTS_SEARCH_WIDTH=50                !<-- Search this far east/west/south/north from problem
                                                                           !    point to fix moving nest sfc-type conflicts.
!
      INTEGER(kind=KINT),PARAMETER :: N_PTS_SEARCH=                     &  !<-- Search this many surrounding pts to fix moving nest 
                                       (2*N_PTS_SEARCH_WIDTH+1)         &  !    conflicts in sfc-type
                                      *(2*N_PTS_SEARCH_WIDTH+1) 
!
      INTEGER(kind=KINT),DIMENSION(1:N_PTS_SEARCH),SAVE :: I_SEARCH_INC &  !<-- I increment to search pt when fixing moving nest conflicts
                                                          ,J_SEARCH_INC    !<-- J increment to search pt when fixing moving nest conflicts
!
      INTEGER(kind=KINT),DIMENSION(:),POINTER,SAVE :: MY_CHILDREN_ID       !<-- A parent's children's domain IDs
!
      INTEGER(kind=KINT),DIMENSION(:),POINTER :: NTIMESTEP_CHILD_MOVES
!
      REAL(kind=KFPT),SAVE :: EPS=1.E-4
!
      REAL(kind=KFPT),SAVE :: ACDT                                      &  !<-- A divergence damping coefficient
                             ,CDDAMP                                    &  !<-- Divergence damping coefficient
                             ,DEG_TO_RAD                                &  !<-- To convert from degrees to radians
                             ,DLM                                       &  !<-- Nest grid increment in X (radians)
                             ,DPH                                       &  !<-- Nest grid increment in Y (radians)
                             ,DT_REAL                                   &  !<-- The dynamical timestep (s)
                             ,RECIP_DPH_1,RECIP_DLM_1                   &  !<-- Reciprocals of upper parent grid increments (radians)
                             ,SB_1,WB_1                                 &  !<-- Rotated S/W bndries of upper parent grid (radians, N/E)
                             ,TLM0                                      &  !<-- Central longitude of domain (radians)
                             ,TPH0                                      &  !<-- Central latitude of domain (radians)
                             ,TPH0_1,TLM0_1                             &  !<-- Central lat/lon of upper parent domain (radians, N/E)
                             ,WCOR
!
      LOGICAL(kind=KLOG),SAVE :: DOMAIN_MOVES                           &  !<-- Does my nested domain move?
                                ,GLOBAL_TOP_PARENT                         !<-- Is the uppermost parent a global domain?
!
      LOGICAL(kind=KLOG) :: I_AM_A_NEST                                 &  !<-- Is the domain a nest?
                           ,INPUT_READY                                 &  !<-- If a nest, does its input file already exist?
                           ,MY_DOMAIN_MOVES                                !<-- Does this domain move?
!
      CHARACTER(len=7) :: SFC_CONFLICT
!
      TYPE :: DIST
        REAL(kind=KFPT) :: VALUE
        INTEGER(kind=KINT) :: I_INC
        INTEGER(kind=KINT) :: J_INC
        TYPE(DIST),POINTER :: NEXT_VALUE
      END TYPE
!
      TYPE(DIST),DIMENSION(:),POINTER,SAVE :: LARGEX,SMALLX
!
!---------------------------------
!***  For determining clocktimes.
!---------------------------------
!
      REAL(kind=KDBL) :: btim,btim0
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE DOMAIN_REGISTER(DOMAIN_GRID_COMP,RC_REG)
! 
!-----------------------------------------------------------------------
!***  Register the Domain component's Initialize, Run, and Finalize
!***  routines.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(ESMF_GridComp) :: DOMAIN_GRID_COMP                              !<-- Domain gridded component
!
      INTEGER,INTENT(OUT) :: RC_REG                                        !<-- Return code for register
!     
!-----------------------------------------------------------------------
!***  Local Variables
!-----------------------------------------------------------------------
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
!***  Register the Domain Initialize subroutine.  Since it is just one 
!***  subroutine, use ESMF_SINGLEPHASE.  The second argument is
!***  a pre-defined subroutine type, such as ESMF_SETINIT, ESMF_SETRUN, 
!***  or ESMF_SETFINAL.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create/Load the Configure Object"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!
      CALL ESMF_GridCompSetEntryPoint(DOMAIN_GRID_COMP                  &  !<-- Domain gridded component
                                     ,ESMF_METHOD_INITIALIZE            &  !<-- Subroutine type (Initialize)
                                     ,DOMAIN_INITIALIZE                 &  !<-- User's subroutine name
                                     ,phase=1                           &
                                     ,rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
!-----------------------------------------------------------------------
!***  Register the Run step of the Domain component.
!***  The NMM needs three phases of Run.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
      MESSAGE_CHECK="Set 1st Entry Point for the Domain Run Step"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
      CALL ESMF_GridCompSetEntryPoint(DOMAIN_GRID_COMP                  &  !<-- The Domain component
                                     ,ESMF_METHOD_RUN                   &  !<-- Subroutine type (Run)
                                     ,DOMAIN_RUN                        &  !<-- The user's subroutine name for primary integration
                                     ,phase=1                           &
                                     ,rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
      MESSAGE_CHECK="Set 2nd Entry Point for the Domain Run Step"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
      CALL ESMF_GridCompSetEntryPoint(DOMAIN_GRID_COMP                  &  !<-- The Domain component
                                     ,ESMF_METHOD_RUN                   &  !<-- Subroutine type (Run)
                                     ,NMM_FILTERING                     &  !<-- Routine to govern digital filtering each timestep
                                     ,phase=2                           &
                                     ,rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
      MESSAGE_CHECK="Set 3rd Entry Point for the Domain Run Step"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
      CALL ESMF_GridCompSetEntryPoint(DOMAIN_GRID_COMP                  &  !<-- The Domain component
                                     ,ESMF_METHOD_RUN                   &  !<-- Subroutine type (Run)
                                     ,CALL_WRITE_ASYNC                  &  !<-- Routine to call asynchronous output
                                     ,phase=3                           &
                                     ,rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
!-----------------------------------------------------------------------
!***  Register the Domain Finalize subroutine.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
      MESSAGE_CHECK="Set Entry Point for Domain Finalize"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
      CALL ESMF_GridCompSetEntryPoint(DOMAIN_GRID_COMP                  &  !<-- The Domain component
                                     ,ESMF_METHOD_FINALIZE              &  !<-- Subroutine type (Finalize)
                                     ,DOMAIN_FINALIZE                   &  !<-- User's subroutine name
                                     ,phase=1                           &
                                     ,rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
!-----------------------------------------------------------------------
!***  Check the error signal variable and print out the result.
!-----------------------------------------------------------------------
!
      IF(RC_REG==ESMF_SUCCESS)THEN
!       WRITE(0,*)' DOMAIN_REGISTER succeeded'
      ELSE
        WRITE(0,*)' DOMAIN_REGISTER failed  RC_REG=',RC_REG
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE DOMAIN_REGISTER
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE DOMAIN_INITIALIZE(DOMAIN_GRID_COMP                     &
                                  ,IMP_STATE                            &
                                  ,EXP_STATE                            &
                                  ,CLOCK_DOMAIN                         &
                                  ,RC_INIT)
!
!-----------------------------------------------------------------------
!***  This routine sets up fundamental aspects of the model run.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(ESMF_GridComp) :: DOMAIN_GRID_COMP                              !<-- The Domain component
!
      TYPE(ESMF_State) :: IMP_STATE                                     &  !<-- The Domain component's import state
                         ,EXP_STATE                                        !<-- The Domain component's export state

      TYPE(ESMF_Clock) :: CLOCK_DOMAIN                                     !<-- The ESMF Clock from the NMM component.
!
      INTEGER,INTENT(OUT) :: RC_INIT                                       !<-- Return code for Initialize step
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: CONFIG_ID,FILT_METHOD_CHILD                 &
                           ,I,ISTAT,J,LB,LNSH,LNSV                      &
                           ,MAX_DOMAINS,N,NLEV_H,NLEV_V                 &
                           ,NFCST,NPHS,NTSD,NUM_DOMAINS,NUM_FIELDS      &
                           ,NUM_PES_FCST,NV                             &
                           ,NVARS_BC_2D_H,NVARS_BC_3D_H,NVARS_BC_4D_H   &
                           ,NVARS_BC_2D_V,NVARS_BC_3D_V                 &
                           ,NVARS_NESTBC_H,NVARS_NESTBC_V               &
                           ,SFC_FILE_RATIO,UB,UBOUND_VARS,NTRACK
!
      INTEGER(kind=KINT) :: IYEAR_FCST                                  &  !<-- Current year from restart file
                           ,IMONTH_FCST                                 &  !<-- Current month from restart file
                           ,IDAY_FCST                                   &  !<-- Current day from restart file
                           ,IHOUR_FCST                                  &  !<-- Current hour from restart file
                           ,IMINUTE_FCST                                &  !<-- Current minute from restart file
                           ,ISECOND_FCST                                   !<-- Current second from restart file
!
      INTEGER(kind=KINT) :: DT_INT,DT_DEN,DT_NUM                        &  !<-- Integer,fractional parts of timestep
                           ,NHOURS_CLOCKTIME                               !<-- Hours between clocktime prints

      INTEGER(kind=KINT) :: FILT_DT_INT,FILT_DT_DEN,FILT_DT_NUM         &  !<-- Integer,fractional parts of timestep used
                           ,FILT_METHOD                                    !    by the digital filter, plus the method
!
      INTEGER(kind=KINT) :: IHI,ILO,JHI,JLO
!
      INTEGER(kind=KINT) :: I_PAR_STA, J_PAR_STA                        &
                           ,LAST_STEP_MOVED,NEXT_MOVE_TIMESTEP          &
                           ,TRACKER_IFIX,TRACKER_JFIX
!
      INTEGER(kind=KINT) :: IERR,IRTN,RC          
!
      INTEGER(ESMF_KIND_I8) :: NTSD_START                                  !<-- Timestep count (>0 for restarted runs)
!
      INTEGER(kind=KINT),DIMENSION(2) :: STORM_CENTER
!
      INTEGER(kind=KINT),DIMENSION(7) :: FCSTDATE
!
      INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE :: DOMAIN_ID_TO_RANK     !<-- Associate configure file IDs with domains
!
      REAL(kind=KFPT) :: CODAMP                                         &  
                        ,DLMD,DPHD                                      &  !<-- Current second from restart file
                        ,DPH_1,DLM_1                                    &
                        ,SECOND_FCST                                    &  !<-- Current second from restart file
                        ,SMAG2                                          &  !<-- Smagorinsky constant
                        ,TLM0D                                          &  !<-- Central longitude of uppermost parent (degrees)
                        ,TPH0D                                             !<-- Central latitude of uppermost parent (degrees)
!
      REAL(kind=KFPT),DIMENSION(:,:),POINTER :: SEA_MASK=>NULL()
!
      LOGICAL(kind=KLOG),SAVE :: ALLOC_FLAG=.FALSE.
!
      LOGICAL(kind=KLOG) :: CALL_BUILD_MOVE_BUNDLE                      &
                           ,CFILE_EXIST                                 &
                           ,I_AM_ACTIVE                                 &
                           ,INPUT_READY_MY_CHILD                        &
                           ,NEMSIO_INPUT                                &
                           ,OPENED
!
      LOGICAL(kind=KLOG) :: I_AM_A_FCST_TASK                            &
                           ,I_AM_A_PARENT 
!
      LOGICAL(kind=KLOG),DIMENSION(:),ALLOCATABLE :: CHILD_ACTIVE
!
      CHARACTER(2)  :: INT_TO_CHAR
      CHARACTER(5)  :: NEST_MODE    
      CHARACTER(6)  :: FMT='(I2.2)'
      CHARACTER(64) :: RESTART_FILENAME
      CHARACTER(99) :: BUNDLE_NAME                                      &
                      ,CONFIG_FILE_NAME                                 &
                      ,FIELD_NAME
!
      TYPE(ESMF_Time) :: CURRTIME                                       &  !<-- The ESMF current time.
                        ,STARTTIME                                         !<-- The ESMF start time.
!
      TYPE(ESMF_Grid) :: GRID_DOMAIN                                       !<-- The ESMF GRID for the integration attached to
                                                                           !     the NMM Domain component.
      TYPE(ESMF_Grid) :: GRID_SOLVER                                       !<-- The ESMF GRID for the integration attached to
                                                                           !     the NMM Solver gridded component.
      TYPE(ESMF_Field) :: HOLD_FIELD
!
      TYPE(ESMF_State) :: IMP_STATE_WRITE
!
      TYPE(NEMSIO_GFILE) :: GFILE
!
      TYPE(WRAP_SOLVER_INT_STATE) :: WRAP_SOLVER
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC     =ESMF_SUCCESS
      RC_INIT=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  Initialize timing variables.
!-----------------------------------------------------------------------
!
      btim0=timef()
!
!-----------------------------------------------------------------------
!***  Take the domain count and this domain's ID from the
!***  import state.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Domain Count from Domain Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state       =IMP_STATE                     &  !<-- The Domain import state
                            ,name        ='NUM_DOMAINS'                 &  !<-- Name of the domain count
                            ,value       =NUM_DOMAINS                   &  !<-- The # of domains
                            ,defaultValue=1                             &  !<-- The default value
                            ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Domain_Initialize: Extract Domain ID from Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state       =IMP_STATE                     &  !<-- The Domain import state
                            ,name        ='DOMAIN_ID'                   &  !<-- Name of the attribute to extract
                            ,value       =MY_DOMAIN_ID                  &  !<-- The ID of this domain
                            ,defaultValue=1                             &  !<-- The default value
                            ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Initialize timers for this domain.
!-----------------------------------------------------------------------
!
      IF(.NOT.ALLOCATED(TIMERS) ) ALLOCATE(TIMERS(1:NUM_DOMAINS))
      timers(my_domain_id)%total_integ_tim=0.
      timers(my_domain_id)%update_interior_from_nest_tim  =0.
      timers(my_domain_id)%update_interior_from_parent_tim=0.
!
!-----------------------------------------------------------------------
!***  To allow a given MPI task to lie on more than one domain
!***  the domain's internal state will be an element of an array
!***  so that each internal state is unique.
!***  It might be more straight forward to allocate the domain 
!***  internal state array alongside the creation of the domain
!***  components themselves but that takes place in the Init step
!***  of the NMM component and we do not want the Domain internal 
!***  state module to be visible there.
!-----------------------------------------------------------------------
!
      IF(.NOT.ALLOC_FLAG)THEN
        ALLOCATE(DOMAIN_INT_STATE_ALL(1:NUM_DOMAINS),stat=RC)
        IF(RC/=0)THEN
          WRITE(0,*)' Failed to allocate DOMAIN_INT_STATE_ALL array  rc=',RC
        ENDIF
        ALLOC_FLAG=.TRUE.
      ENDIF
!
      wrap%DOMAIN_INT_STATE=>DOMAIN_INT_STATE_ALL(MY_DOMAIN_ID)
      DOMAIN_INT_STATE     =>DOMAIN_INT_STATE_ALL(MY_DOMAIN_ID)
!
!-----------------------------------------------------------------------
!***  Attach the Domain internal state to the Domain component.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Attach Domain Internal State to Gridded Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompSetInternalState(DOMAIN_GRID_COMP               &  !<-- The Domain gridded component
                                        ,WRAP                           &  !<-- Pointer to the Domain internal state
                                        ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Retrieve the VM (Virtual Machine) of the Domain component.
!***  Call ESMF_GridCompGet to retrieve the VM anywhere you need it.
!***  We need VM now to obtain the MPI task IDs and the local MPI
!***  communicator.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="DOMAIN_INIT: Retrieve VM from Domain Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompGet(gridcomp=DOMAIN_GRID_COMP                   &  !<-- The Domain component
                           ,vm      =VM                                 &  !<-- Get the Virtual Machine from the Domain component
                           ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="DOMAIN_INIT: Obtain Task IDs and Communicator"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_VMGet(vm             =VM                                &  !<-- The virtual machine
                     ,localpet       =MYPE                              &  !<-- Each MPI task ID
                     ,mpiCommunicator=COMM_MY_DOMAIN                    &  !<-- This domain's intracommunicator
                     ,rc             =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Extract the maximum number of domains from the import state.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract MAX_DOMAINS from Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state       =IMP_STATE                     &  !<-- The Domain import state
                            ,name        ='MAX_DOMAINS'                 &  !<-- Name of the attribute to extract
                            ,value       =MAX_DOMAINS                   &  !<-- Maximum # of domains
                            ,defaultValue=1                             &  !<-- The default value
                            ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Extract the configure file IDs associated with each domain.
!-----------------------------------------------------------------------
!
      ALLOCATE(DOMAIN_ID_TO_RANK(1:MAX_DOMAINS),stat=ISTAT)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Association of Configure Files with Domains"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state           =IMP_STATE                 &  !<-- The Domain import state
                            ,name            ='DOMAIN_ID_TO_RANK'       &  !<-- Name of the attribute to extract
                            ,itemCount       =MAX_DOMAINS               &  !<-- Name of the attribute to extract
                            ,valueList       =DOMAIN_ID_TO_RANK         &  !<-- The ID of this domain
                            ,defaultvalueList=[1]                       &  !<-- The default valueList
                            ,rc              =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Now we can load configure files for all domains into memory.
!***  The file name of the uppermost domain is 'configure_file_01'
!***  and is identical to the primary file called 'configure_file'
!***  which is needed in some early parts of the setup.
!-----------------------------------------------------------------------
!
      DO N=1,MAX_DOMAINS                                                   !<-- The number of config files cannot exceed 99
!
        CONFIG_ID=DOMAIN_ID_TO_RANK(N)
        WRITE(INT_TO_CHAR,FMT)CONFIG_ID 
        CONFIG_FILE_NAME='configure_file_'//INT_TO_CHAR                    !<-- Prepare the config file names
!
        CFILE_EXIST=.FALSE.
        INQUIRE(file=CONFIG_FILE_NAME,exist=CFILE_EXIST)
!
        IF(CFILE_EXIST)THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Create the Nest Configure Object"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CF(N)=ESMF_ConfigCreate(rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Load the Nest Configure Object"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_ConfigLoadFile(config  =CF(N)                       &
                                  ,filename=CONFIG_FILE_NAME            &
                                  ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        ELSE
!
          EXIT
!
        ENDIF
!
      ENDDO
!
      DEALLOCATE(DOMAIN_ID_TO_RANK)
!
!-----------------------------------------------------------------------
!***  If the atmosphere is coupled to an ocean then the 
!***  atm.configure file holds the coupling time interval.
!***  Extract it so we know when to update the atmosphere's
!***  SST in DOMAIN_RUN.
!-----------------------------------------------------------------------
!
      CONFIG_FILE_NAME='atmos.configure'                                   !<-- The config file name
      CFILE_EXIST=.FALSE.
      INQUIRE(file=CONFIG_FILE_NAME,exist=CFILE_EXIST)
!
      IF(CFILE_EXIST)THEN
!
        CF_ATMOS=ESMF_ConfigCreate(rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Load the Nest Configure Object"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigLoadFile(config  =CF_ATMOS                      &
                                ,filename=CONFIG_FILE_NAME              &
                                ,rc      =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract Atmos-Ocean Coupling Interval from Config File"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=CF_ATMOS                      &  !<-- The config object
                                    ,value =ATM_OCN_CPL_INT               &  !<-- The atm-ocean coupling interval (sec) 
                                    ,label ='atm_coupling_interval_sec:'  &  !<-- Give this label's value to the previous variable
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF
!
!-----------------------------------------------------------------------
!***  Shall we write last time step restart file?
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Write_last_restart Flag from Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)              &  !<-- The config object
                                  ,value =WRITE_LAST_RESTART            &  !<-- The quilting flag
                                  ,label ='write_last_restart:'         &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      domain_int_state%WRITE_LAST_RESTART=WRITE_LAST_RESTART               !<-- Save this for the write_async
!
!
!-----------------------------------------------------------------------
!***  Will the Write components with asynchronous quilting be used?
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Quilting Flag from Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)              &  !<-- The config object
                                  ,value =QUILTING                      &  !<-- The quilting flag
                                  ,label ='quilting:'                   &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      domain_int_state%QUILTING=QUILTING                                   !<-- Save this for the Run step
!
!-----------------------------------------------------------------------
!***  Initialize the flag indicating if the first history output has
!***  been written out.
!-----------------------------------------------------------------------
!
      domain_int_state%WROTE_1ST_HIST=.FALSE.
!
!-----------------------------------------------------------------------
!***  The task layout on this domain.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract MPI Task Layout in Domain Init"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)              &  !<-- The config object
                                  ,value =INPES                         &  !<-- The # of forecast tasks in I
                                  ,label ='inpes:'                      &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
      CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)              &  !<-- The config object
                                  ,value =JNPES                         &  !<-- The # of forecast tasks in J
                                  ,label ='jnpes:'                      &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Extract the physics call frequency from the configure file.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract NPHS in Domain Init"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)              &  !<-- The config object
                                  ,value =NPHS                          &  !<-- The physics call frequency
                                  ,label ='nphs:'                       &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Insert NPHS into the Domain export state"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(state=EXP_STATE                            &  !<-- The config object
                            ,name ='NPHS'                               &  !<-- The name in the export state
                            ,value=NPHS                                 &  !<-- The physics call frequency
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!-----------------------------------------------------------------------
!***  Extract the domain boundary blending width.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract LNSH,LNSV in Domain Init"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)              &  !<-- The config object
                                  ,value =LNSH                          &  !<-- Domain bndry blending width for H points
                                  ,label ='lnsh:'                       &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
      CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)              &  !<-- The config object
                                  ,value =LNSV                          &  !<-- Domain bndry blending width for V points
                                  ,label ='lnsv:'                       &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Extract this Domain component's Nest/Not-a-Nest flag
!***  from the import state.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Nest/Not-a-Nest Flag from Domain Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state       =IMP_STATE                     &  !<-- The Domain import state
                            ,name        ='I-Am-A-Nest Flag'            &  !<-- Name of the attribute to extract
                            ,value       =I_AM_A_NEST                   &  !<-- The flag indicating if this domain is a nest
                            ,defaultValue=.false.                       &  !<-- The default value
                            ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Save the nest flag in the Domain's internal state so it can be
!***  referred to in the Run step.
!-----------------------------------------------------------------------
!
      domain_int_state%I_AM_A_NEST=I_AM_A_NEST
!
!-----------------------------------------------------------------------
!***  Extract the ratio of the parent timestep to the child's if
!***  this domain is a nest.
!***  Also extract the flag indicating whether or not the nest's
!***  input file has already been generated by NPS.
!-----------------------------------------------------------------------
!
      IF(I_AM_A_NEST)THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract Parent-Child Time Ratio from Domain Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state=IMP_STATE                          &  !<-- The Domain import state
                              ,name ='Parent-Child Time Ratio'          &  !<-- Name of Attribute
                              ,value=PARENT_CHILD_TIME_RATIO            &  !<-- Ratio of this domain's parent's timestep to its own
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract Input Ready Flag from Configure File"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)            &  !<-- The config object
                                    ,value =INPUT_READY                 &  !<-- The variable filled (does nest input file exist?
                                    ,label ='input_ready:'              &  !<-- The input datafile for this domain does or does not exist
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract Parent-Child Space Ratio from Configure File"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)            &  !<-- The config object
                                    ,value =PARENT_CHILD_SPACE_RATIO    &  !<-- The variable filled (child grid increment / parent's)
                                    ,label ='parent_child_space_ratio:' &  !<-- Give this label's value to the previous variable
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  The nest must know whether or not it moves.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract Move Flag From Nest Configure Files"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)            &  !<-- This domain's configure object
                                    ,value =MY_DOMAIN_MOVES             &  !<-- Does this domain move?
                                    ,label ='my_domain_moves:'          &  !<-- The label in the configure file
                                    ,rc    =rc)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        domain_int_state%MY_DOMAIN_MOVES=MY_DOMAIN_MOVES
!
!-----------------------------------------------------------------------
!
      ENDIF
!
!-----------------------------------------------------------------------
!***  Does this forecast use 2-way exchange?
!-----------------------------------------------------------------------
!
      NEST_MODE=' '
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract 2-way Flag From Nest Configure Files"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)              &  !<-- This domain's configure object
                                  ,value =NEST_MODE                     &  !<-- Is there 2-way exchange from child to parent?
                                  ,label ='nest_mode:'                  &  !<-- The label in the configure file
                                  ,rc    =rc)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Extract the start time from the clock.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="DOMAIN_INIT: Start Time from Domain Clock"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ClockGet(clock    =CLOCK_DOMAIN                         &  !<-- The ESMF Clock of this domain
                        ,startTime=STARTTIME                            &  !<-- The simulation start time
                        ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CURRTIME=STARTTIME
      NTSD_START=0
!
!-----------------------------------------------------------------------
!***  Extract the NEMSIO_INPUT flag from the configure file.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract NEMSIO Flag from Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)            &  !<-- The config object
                                  ,value =NEMSIO_INPUT                &  !<-- The input datafile does or does not have NEMSIO metadata
                                  ,label ='nemsio_input:'             &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Extract the RESTART flag from the configure file.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Restart Flag from Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)              &  !<-- The config object
                                  ,value =RESTARTED_RUN                 &  !<-- True => restart; False => cold start
                                  ,label ='restart:'                    &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      domain_int_state%RESTARTED_RUN=RESTARTED_RUN
!
      domain_int_state%RESTARTED_RUN_FIRST=.TRUE.                          !<-- Prepare for the initial output for a restarted run.
!
!-----------------------------------------------------------------------
!***  If this is a restarted run then read:
!***    (1) The forecast time that the file was written.
!***    (2) The forecast timestep at which the file was written.
!-----------------------------------------------------------------------
!
      NTSD_START=0
!
      restart: IF(RESTARTED_RUN)THEN                                       !<-- If this is a restarted run, set the current time
!
        WRITE(INT_TO_CHAR,FMT)MY_DOMAIN_ID
!
!----------------------------------------------------------------------
!***  Read the restart data from either pure binary or NEMSIO file.
!-----------------------------------------------------------------------
!
        input: IF(NEMSIO_INPUT)THEN
!
          CALL NEMSIO_INIT()
!
          RESTART_FILENAME='restart_file_'//INT_TO_CHAR//'_nemsio'
          CALL NEMSIO_OPEN(GFILE,RESTART_FILENAME,'read',iret=IRTN)
          IF(IRTN/=0)THEN
            WRITE(0,*)' Unable to open nemsio file '                    &
                     ,TRIM(RESTART_FILENAME),' in DOMAIN_INITIALIZE'
            WRITE(0,*)' ABORTING!'
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT               &
                              ,rc             =RC)
          ENDIF
!
          CALL NEMSIO_GETHEADVAR(GFILE,'FCSTDATE',FCSTDATE,iret=irtn)
!
          IYEAR_FCST  =FCSTDATE(1)
          IMONTH_FCST =FCSTDATE(2)
          IDAY_FCST   =FCSTDATE(3)
          IHOUR_FCST  =FCSTDATE(4)
          IMINUTE_FCST=FCSTDATE(5)
          SECOND_FCST =0.
!
          IF(FCSTDATE(7)/=0)THEN
            SECOND_FCST=FCSTDATE(6)/(FCSTDATE(7)*1.)
          ENDIF
!
          CALL NEMSIO_GETHEADVAR(gfile,'NTIMESTEP',NTSD,iret=irtn)

          CALL NEMSIO_CLOSE(GFILE,iret=IERR)
!
        ELSE                                                               !<-- Pure binary input
!
          select_unit: DO N=51,59
            INQUIRE(N,OPENED=OPENED)
            IF(.NOT.OPENED)THEN
              NFCST=N
              EXIT select_unit
            ENDIF
          ENDDO select_unit
!
          RESTART_FILENAME='restart_file_'//INT_TO_CHAR
          OPEN(unit=NFCST,file=RESTART_FILENAME,status='old'            &
              ,form='unformatted',iostat=IRTN)
          IF(IRTN/=0)THEN
            WRITE(0,*)' Unable to open pure binary file '               &
                     ,TRIM(RESTART_FILENAME),' in DOMAIN_INITIALIZE'
            WRITE(0,*)' ABORTING!'
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT               &
                              ,rc             =RC)
          ENDIF
!
          READ(NFCST) IYEAR_FCST                                           !<-- Read time form restart file
          READ(NFCST) IMONTH_FCST                                          !
          READ(NFCST) IDAY_FCST                                            !
          READ(NFCST) IHOUR_FCST                                           !
          READ(NFCST) IMINUTE_FCST                                         !
          READ(NFCST) SECOND_FCST                                          !<--
!
          READ(NFCST) NTSD                                                 !<-- Read timestep from restart file
!
          CLOSE(NFCST)
!
        ENDIF  input
!
!-----------------------------------------------------------------------
!
        ISECOND_FCST=NINT(SECOND_FCST)                                     !<-- ESMF clock needs integer seconds
        NTSD_START=NTSD
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="RESTART: Set the Current Time of the Forecast"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_TimeSet(time=CURRTIME                                 &  !<-- Current time of the forecast (ESMF)
                         ,yy  =IYEAR_FCST                               &  !<-- Year from restart file
                         ,mm  =IMONTH_FCST                              &  !<-- Month from restart file
                         ,dd  =IDAY_FCST                                &  !<-- Day from restart file
                         ,h   =IHOUR_FCST                               &  !<-- Hour from restart file
                         ,m   =IMINUTE_FCST                             &  !<-- Minute from restart file
                         ,s   =ISECOND_FCST                             &  !<-- Second from restart file
                         ,rc  =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
      ENDIF restart
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  With data from above set the local ESMF Clock
!***  to its correct time and timestep.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set the Current Time on the Domain Clock"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ClockSet(clock       =CLOCK_DOMAIN                      &  !<-- The Domain Component's Clock
                        ,currtime    =CURRTIME                          &  !<-- Current time of simulation
                        ,advanceCount=NTSD_START                        &  !<-- Timestep at this current time
                        ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
! 
!-----------------------------------------------------------------------
!***  Create the time interval for printing clocktimes used by model 
!***  sections.  Read in forecast time interval for clocktime output 
!***  as well as the selected task ID that will provide the clocktimes.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Read Fcst Interval for Clocktime Output"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)              &  !<-- The configure object
                                  ,value =NHOURS_CLOCKTIME              &  !<-- Fill this variable
                                  ,label ='nhours_clocktime:'           &  !<-- Give the variable this label's value from the config file
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Read MPI Task ID That Provides Clocktime Output"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)              &  !<-- The configure object
                                  ,value =NPE_PRINT                     &  !<-- Fill this variable
                                  ,label ='npe_print:'                  &  !<-- Give the variable this label's value from the config file
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  How many tracer species are there?
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="DOMAIN_Init: Extract # of tracers from Config file"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)              &  !<-- The config object
                                  ,value =NUM_TRACERS_CHEM              &  !<-- The variable filled (number of chemical tracers)
                                  ,label ='num_tracers_chem:'           &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  Model-specific routines must be invoked in order to establish
!***  the ESMF Grid.  The different integration grids necessitate
!***  different ways of setting up both the parallelism for
!***  distributed memory runs and the ESMF Grid itself.
!***  When the parallelism is constructed, the local domain limits
!***  need to be inserted into the Domain component's internal state
!***  if quilting is to be used.  See 'IF(QUILTING)THEN' below.
!-----------------------------------------------------------------------
!
      WRITE(0,*)' '
      WRITE(0,11110)MY_DOMAIN_ID
11110 format(' DOMAIN_SETUP my_domain_id=',i2)
      CALL DOMAIN_SETUP(MYPE                                            &
                       ,COMM_MY_DOMAIN                                  &
                       ,QUILTING                                        &
                       ,CF(MY_DOMAIN_ID)                                &
                       ,DOMAIN_GRID_COMP                                &
                       ,DOMAIN_INT_STATE                                &
                       ,GRID_DOMAIN )
!
!-----------------------------------------------------------------------
!***  Save the intracommunicator for this domain's forecast tasks
!***  since it is used each timestep in DOMAIN_RUN.
!-----------------------------------------------------------------------
!
      IF(.NOT.ALLOCATED(COMM_FCST_TASKS))THEN
        ALLOCATE(COMM_FCST_TASKS(1:NUM_DOMAINS))
      ENDIF
!
      COMM_FCST_TASKS(MY_DOMAIN_ID)=MPI_COMM_COMP
!
!-----------------------------------------------------------------------
!***  Attach the NMM-specific ESMF Grid to the Domain component.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Attach the NMM ESMF Grid to the Domain Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompSet(gridcomp=DOMAIN_GRID_COMP                   & !<-- The Domain component
                           ,grid    =GRID_DOMAIN                        & !<-- Attach the ESMF grid to the Domain component
                           ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  Create the Solver gridded subcomponent.
!***  Register the Initialize, Run, and Finalize steps for it.
!***  Since there is only a single integration grid, give the
!***  Solver the Domain component's grid.
!***  Note that this subcomponent is part of the Domain component's
!***  internal state.  This will be convenient if we need to reach
!***  the Solver component via the Domain component such as happens
!***  when Write components are established.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!---------------------------------
!***  Create the Solver component
!---------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create the NMM Solver Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      domain_int_state%SOLVER_GRID_COMP=ESMF_GridCompCreate(            &
                                  name   ="Solver component"            &  !<-- Name of the new Solver gridded component
                                 ,config =CF(MY_DOMAIN_ID)              &  !<-- Attach this configure file to the component
                                 ,petList=domain_int_state%PETLIST_FCST &  !<-- The LOCAL forecast task IDs
                                 ,rc     =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!------------------------------------------------
!***  Register the Init, Run, and Finalize steps
!------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Register the NMM Solver Init, Run, Finalize"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompSetServices(domain_int_state%SOLVER_GRID_COMP      &  ! <-- The Solver gridded component
                                   ,SOLVER_REGISTER                        &  ! <-- The user's subroutineName for Register
                                   ,rc = RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!----------------------------------------------------
!***  Attach the ESMF Grid to the Solver component
!----------------------------------------------------
!
      GRID_SOLVER=GRID_DOMAIN                                              !<-- The Solver grid is the same as the Domain grid
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Attach the ESMF Grid to the Solver Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompSet(gridcomp=domain_int_state%SOLVER_GRID_COMP  &  !<-- The Solver component
                           ,grid    =GRID_SOLVER                        &  !<-- The Solver ESMF grid
                           ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Create empty import and export states for the Solver component.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create Empty Import/Export States for the Solver"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      domain_int_state%IMP_STATE_SOLVER=ESMF_StateCreate(                       &
                                            name  ="Solver Import"         &  !<-- The Solver import state name
                                           ,stateintent=ESMF_STATEINTENT_IMPORT &
                                           ,rc         =RC)
!
      domain_int_state%EXP_STATE_SOLVER=ESMF_StateCreate(                       &
                                            name  ="Solver Export"         &  !<-- The Solver export state name
                                           ,stateintent=ESMF_STATEINTENT_EXPORT &
                                           ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Extract the flag from the Domain import state indicating if the
!***  user wants physics to be active.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Physics Flag from Domain Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state       =IMP_STATE                     &
                            ,name        ='PHYSICS_ON'                  &
                            ,value       =PHYSICS_ON                    &
                            ,defaultValue=.true.                        &
                            ,rc          =rc)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Add fundamental domain characteristics to the Solver's
!***  import state that will be needed by that component.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Add Domain Dimensions to the Solver Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(state=domain_int_state%IMP_STATE_SOLVER    &  !<-- The Solver component import state
                            ,name ='ITS'                                &  !<-- Use this name inside the state
                            ,value=ITS                                  &  !<-- The scalar being inserted into the import state
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=domain_int_state%IMP_STATE_SOLVER    &  !<-- The Solver component import state
                            ,name ='ITE'                                &  !<-- Use this name inside the state
                            ,value=ITE                                  &  !<-- The scalar being inserted into the import state
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=domain_int_state%IMP_STATE_SOLVER    &  !<-- The Solver component import state
                            ,name ='JTS'                                &  !<-- Use this name inside the state
                            ,value=JTS                                  &  !<-- The scalar being inserted into the import state
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=domain_int_state%IMP_STATE_SOLVER    &  !<-- The Solver component import state
                            ,name ='JTE'                                &  !<-- Use this name inside the state
                            ,value=JTE                                  &  !<-- The scalar being inserted into the import state
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=domain_int_state%IMP_STATE_SOLVER    &  !<-- The Solver component import state
                            ,name ='IMS'                                &  !<-- Use this name inside the state
                            ,value=IMS                                  &  !<-- The scalar being inserted into the import state
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=domain_int_state%IMP_STATE_SOLVER    &  !<-- The Solver component import state
                            ,name ='IME'                                &  !<-- Use this name inside the state
                            ,value=IME                                  &  !<-- The scalar being inserted into the import state
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=domain_int_state%IMP_STATE_SOLVER    &  !<-- The Solver component import state
                            ,name ='JMS'                                &  !<-- Use this name inside the state
                            ,value=JMS                                  &  !<-- The scalar being inserted into the import state
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=domain_int_state%IMP_STATE_SOLVER    &  !<-- The Solver component import state
                            ,name ='JME'                                &  !<-- Use this name inside the state
                            ,value=JME                                  &  !<-- The scalar being inserted into the import state
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=domain_int_state%IMP_STATE_SOLVER    &  !<-- The Solver component import state
                            ,name ='IDS'                                &  !<-- Use this name inside the state
                            ,value=IDS                                  &  !<-- The scalar being inserted into the import state
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=domain_int_state%IMP_STATE_SOLVER    &  !<-- The Solver component import state
                            ,name ='IDE'                                &  !<-- Use this name inside the state
                            ,value=IDE                                  &  !<-- The scalar being inserted into the import state
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=domain_int_state%IMP_STATE_SOLVER    &  !<-- The Solver component import state
                            ,name ='JDS'                                &  !<-- Use this name inside the state
                            ,value=JDS                                  &  !<-- The scalar being inserted into the import state
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=domain_int_state%IMP_STATE_SOLVER    &  !<-- The Solver component import state
                            ,name ='JDE'                                &  !<-- Use this name inside the state
                            ,value=JDE                                  &  !<-- The scalar being inserted into the import state
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Add Halo Widths to Solver Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(state=domain_int_state%IMP_STATE_SOLVER    &  !<-- The Solver component import state
                            ,name ='IHALO'                              &  !<-- Use this name inside the state
                            ,value=IHALO                                &  !<-- The scalar being inserted into the import state
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=domain_int_state%IMP_STATE_SOLVER    &  !<-- The Solver component import state
                            ,name ='JHALO'                              &  !<-- Use this name inside the state
                            ,value=JHALO                                &  !<-- The scalar being inserted into the import state
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Add Fcst/Quilt Intracomms to the Solver Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(state=domain_int_state%IMP_STATE_SOLVER    &  !<-- The Solver component import state
                            ,name ='Fcst/Quilt Intracommunicators'      &  !<-- Use this name inside the state
                            ,value=MPI_COMM_COMP                        &  !<-- The scalar being inserted into the import state
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Add Domain ID to the Solver Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(state=domain_int_state%IMP_STATE_SOLVER    &  !<-- The Solver component import state
                            ,name ='DOMAIN_ID'                          &  !<-- Use this name inside the state
                            ,value=MY_DOMAIN_ID                         &  !<-- The scalar being inserted into the import state
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Add Task Neighbors to the Solver Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
! 
      CALL ESMF_AttributeSet(state    =domain_int_state%IMP_STATE_SOLVER &  !<-- The Solver component import state
                            ,name     ='MY_NEB'                          &  !<-- Use this name inside the state
                            ,itemCount=N8                                &  !<-- # of items in Attribute
                            ,valueList=MY_NEB                            &  !<-- The scalar being inserted into the import state
                            ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Add Task Neighbors to the Domain Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
! 
      CALL ESMF_AttributeSet(state    =EXP_STATE                         &  !<-- The Domain component export state
                            ,name     ='MY_NEB'                          &  !<-- Use this name inside the state
                            ,itemCount=N8                                &  !<-- # of items in Attribute
                            ,valueList=MY_NEB                            &  !<-- The scalar being inserted into the import state
                            ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!------------------------------------------------------------------------
!***  Insert the flag indicating if the Domain component is a nest.
!***  The Solver component needs to know this regarding BC's in
!***  order to properly compute fundamental aspects of the
!***  nested grids.
!------------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="DOMAIN_INIT: Add Nest Flag to the Solver Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(state=domain_int_state%IMP_STATE_SOLVER    &  !<-- The Solver component import state
                            ,name ='I-Am-A-Nest Flag'                   &  !<-- Use this name inside the state
                            ,value=I_AM_A_NEST                          &  !<-- The logical being inserted into the import state
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!------------------------------------------------------------------------
!***  Is the uppermost parent on a global domain?  We must know this
!***  for moving nests' reading the external surface files that span
!***  that domain.
!------------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Is Domain #1 Global?"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF(1)                         &  !<-- The config object of domain #1
                                  ,value =GLOBAL_TOP_PARENT             &  !<-- The variable filled
                                  ,label ='global:'                     &  !<-- True--> uppermost parent is on a global domain.
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!------------------------------------------------------------------------
!***  Add the transformed lat/lon (degrees) of the SW corner of domain #1
!***  domain #1 and the geographic lat/lon of its center to the Solver
!***  import state.  That information will be used if this is a restarted
!***  run containing moving nests in order to precisely determine the
!***  location of those nests.
!------------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Southern/Western Boundary of Uppermost Domain"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF(1)                         &  !<-- The config object of domain #1
                                  ,value =SBD_1                         &  !<-- The variable filled
                                  ,label ='sbd:'                        &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
      CALL ESMF_ConfigGetAttribute(config=CF(1)                         &  !<-- The config object of domain #1
                                  ,value =WBD_1                         &  !<-- The variable filled
                                  ,label ='wbd:'                        &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="DOMAIN_INIT: Add SW Corner of Domain #1 to Solver Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(state=domain_int_state%IMP_STATE_SOLVER    &  !<-- The Solver component import state
                            ,name ='SBD_1'                              &  !<-- Attribute's name
                            ,value=SBD_1                                &  !<-- Transformed lat (degrees) of domain #1's south bndry
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=domain_int_state%IMP_STATE_SOLVER    &  !<-- The Solver component import state
                            ,name ='WBD_1'                              &  !<-- Attribute's name
                            ,value=WBD_1                                &  !<-- Transformed lon (degrees) of domain #1's west bndry
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Central Lat/Lon of Uppermost Domain"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF(1)                         &  !<-- The config object of domain #1
                                  ,value =TPH0D_1                       &  !<-- Geographic lat (degrees) of center of domain #1
                                  ,label ='tph0d:'                      &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
      CALL ESMF_ConfigGetAttribute(config=CF(1)                         &  !<-- The config object of domain #1
                                  ,value =TLM0D_1                       &  !<-- Geographic lon (degrees) of center of domain #1
                                  ,label ='tlm0d:'                      &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="DOMAIN_INIT: Add Center of Domain #1 to Solver Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(state=domain_int_state%IMP_STATE_SOLVER    &  !<-- The Solver component import state
                            ,name ='TPH0D_1'                            &  !<-- Attribute's name
                            ,value=TPH0D_1                              &  !<-- Geographic lat (degrees) of domain #1's center
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=domain_int_state%IMP_STATE_SOLVER    &  !<-- The Solver component import state
                            ,name ='TLM0D_1'                            &  !<-- Attribute's name
                            ,value=TLM0D_1                              &  !<-- Geographic lon (degrees) of domain #1's center
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!------------------------------------------------------------------------
!***  Add the local domain index limits to the Solver import state
!***  on the compute tasks.
!------------------------------------------------------------------------
!
      NUM_PES_FCST=INPES*JNPES
!
      IF(MYPE<NUM_PES_FCST)THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert Local Domain Limits in Solver Imp State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

        CALL ESMF_AttributeSet(state    =domain_int_state%IMP_STATE_SOLVER &  !<-- Solver import state receives an attribute
                              ,name     ='LOCAL_ISTART'                    &  !<-- The attribute's name
                              ,itemCount=NUM_PES_FCST                      &  !<-- The attribute's length
                              ,valueList=domain_int_state%LOCAL_ISTART     &  !<-- Insert this quantity as an attribute
                              ,rc       =RC)
!
        CALL ESMF_AttributeSet(state    =domain_int_state%IMP_STATE_SOLVER &  !<-- Solver import state receives an attribute
                              ,name     ='LOCAL_IEND'                      &  !<-- The attribute's name
                              ,itemCount=NUM_PES_FCST                      &  !<-- The attribute's length
                              ,valueList=domain_int_state%LOCAL_IEND       &  !<-- Insert this quantity as an attribute
                              ,rc       =RC)
!
        CALL ESMF_AttributeSet(state    =domain_int_state%IMP_STATE_SOLVER &  !<-- Solver import state receives an attribute
                              ,name     ='LOCAL_JSTART'                    &  !<-- The attribute's name
                              ,itemCount=NUM_PES_FCST                      &  !<-- The attribute's length
                              ,valueList=domain_int_state%LOCAL_JSTART     &  !<-- Insert this quantity as an attribute
                              ,rc       =RC)
!
        CALL ESMF_AttributeSet(state    =domain_int_state%IMP_STATE_SOLVER &  !<-- Solver import state receives an attribute
                              ,name     ='LOCAL_JEND'                      &  !<-- The attribute's name
                              ,itemCount=NUM_PES_FCST                      &  !<-- The attribute's length
                              ,valueList=domain_int_state%LOCAL_JEND       &  !<-- Insert this quantity as an attribute
                              ,rc       =RC)

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)  
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF
!
!------------------------------------------------------------------------
!***  If this is a nest domain then insert the Parent-Child timestep
!***  ratio into the Solver import state since that will be needed
!***  to tell the Solver Run step how often to update the boundary
!***  tendencies.
!***  Also insert the flag indicating whether or not the nest domain
!***  already has an input file or if one needs to be generated by
!***  its parent.
!***  Finally insert the flag indicating if the nest moves.
!------------------------------------------------------------------------
!
      IF(I_AM_A_NEST)THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="DOMAIN_INIT: Add Parent-Child Time Ratio to the Solver Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=domain_int_state%IMP_STATE_SOLVER  &  !<-- The Solver component import state
                              ,name ='Parent-Child Time Ratio'          &  !<-- Use this name inside the state
                              ,value=PARENT_CHILD_TIME_RATIO            &  !<-- Put the Attribute here
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="DOMAIN_INIT: Add Input-Ready Flag to the Solver Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=domain_int_state%IMP_STATE_SOLVER  &  !<-- The Solver component import state
                              ,name ='Input Ready'                      &  !<-- Use this name inside the state
                              ,value=INPUT_READY                        &  !<-- Does this nest's input file already exist?
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="DOMAIN_INIT: Add Move Flag to the Solver Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=domain_int_state%IMP_STATE_SOLVER  &  !<-- The Solver component import state
                              ,name ='My Domain Moves'                  &  !<-- Use this name inside the state
                              ,value=MY_DOMAIN_MOVES                    &  !<-- Does this nest move?
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="DOMAIN_INIT: Add Move Flag to the Dom Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=IMP_STATE                          &  !<-- The Domain component import state
                              ,name ='My Domain Moves'                  &  !<-- Use this name inside the state
                              ,value=MY_DOMAIN_MOVES                    &  !<-- Does this nest move?
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  If quilting was selected for the generation of output,
!***  set up the Write component(s).
!***  This must be done prior to executing the Initialize steps
!***  of the Solver component because the Write components'
!***  import states are required by those steps when quilting
!***  is active and WRITE_SETUP is the routine in which the 
!***  Write components' import states are inserted into the
!***  Solver export state which is in turn part of the DOMAIN
!***  internal state.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="DOMAIN_INIT: Add Quilting Flag to the Solver Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(state=domain_int_state%IMP_STATE_SOLVER  &  !<-- The Solver component import state
                            ,name ='Quilting'                         &  !<-- Use this name inside the state
                            ,value=QUILTING                           &  !<-- Was quilting specified in the configure file?
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF(QUILTING)THEN
!
        CALL WRITE_SETUP(DOMAIN_GRID_COMP                               &
                        ,DOMAIN_INT_STATE                               &
                        ,CLOCK_DOMAIN )
!
      ENDIF
!
!-----------------------------------------------------------------------
!***  The user specifies the variables to be used for nest boundaries.
!***  The boundary variables for the uppermost parent (or for single
!***  domain runs) are still hard-wired.  All boundary variables must
!***  be part of the Solver component's internal state.  The internal
!***  state is set up early in the Initialize step of the Solver
!***  component.  Later in SOLVER_INITIALIZE the routine is called
!***  that reads in the input or restart files.  For restarted runs
!***  the read routine must allocate a working array to hold boundary
!***  values from the restart file.  All boundary variables are
!***  carried in an ESMF Bundle and in a generalized boundary object.
!***  The quantities needed by the read routine to allocate its 
!***  BC working array are determined when the Bundle and general BC
!***  object are created and filled.  Therefore the routine that
!***  produces the Bundle and general BC object must be called in
!***  SOLVER_INITIALIZE after the Solver internal state is set up
!***  and before the read routine.  All of the user-specified
!***  variables are carried in ESMF Bundles that lie within the
!***  Domain component's internal state.  Create the Bundle for
!***  BC variables now and send it into the Initialize step of the
!***  Solver component.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create the BC Variable Bundle"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      BUNDLE_NAME='Bundle_nestbc'
!
      domain_int_state%BUNDLE_NESTBC=ESMF_FieldBundleCreate(name=BUNDLE_NAME  &  !<-- The Bundle's name for nest BC variables
                                                           ,rc  =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="DOMAIN_INIT: Add Move Flag to the Dom Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateAddReplace(domain_int_state%IMP_STATE_SOLVER            &  !<-- The Solver component's import state
                        ,(/domain_int_state%BUNDLE_NESTBC/)  &  !<-- The Bundle of Solver int state pointers for nest BC vbls
                        ,rc =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  Execute the Initialize step for the Solver subcomponent.
!***  This is the Initialize subroutine specified in the
!***  Register routine called in ESMF_GridCompSetServices above.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Initialize the NMM Solver Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompInitialize(gridcomp   =domain_int_state%SOLVER_GRID_COMP  &  !<-- The Solver gridded component
                                  ,importState=domain_int_state%IMP_STATE_SOLVER  &  !<-- The Solver import state
                                  ,exportState=domain_int_state%EXP_STATE_SOLVER  &  !<-- The Solver export state
                                  ,clock      =CLOCK_DOMAIN                       &  !<-- The Domain clock
                                  ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Extract the Solver internal state so we can access it.
!-----------------------------------------------------------------------
!
      IF(MYPE<NUM_PES_FCST)THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
       MESSAGE_CHECK="Domain_Init: Extract Solver Internal State"
!      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_GridCompGetInternalState(domain_int_state%SOLVER_GRID_COMP &  !<-- The Solver component
                                          ,WRAP_SOLVER                  &
                                          ,RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        SOLVER_INT_STATE=>wrap_solver%INT_STATE
!
!-----------------------------------------------------------------------
!
        LM=solver_int_state%LM                                             !<-- We need LM later in the routine.
!
!-----------------------------------------------------------------------
!***  Tell the Solver whether quilting was selected.
!-----------------------------------------------------------------------
!
!
      ENDIF
!
!-----------------------------------------------------------------------
!***  Parents update the boundaries of their children every parent
!***  timestep.  The nest boundary variables in the Solver internal 
!***  state that are updated are selected by the user.
!-----------------------------------------------------------------------
!
      bc_variables: IF(MYPE<domain_int_state%NUM_PES_FCST)THEN             !<-- Select the compute tasks.
!
        NLEV_H=solver_int_state%NLEV_H
        NLEV_V=solver_int_state%NLEV_V
!
        NVARS_BC_2D_H=solver_int_state%NVARS_BC_2D_H
        NVARS_BC_3D_H=solver_int_state%NVARS_BC_3D_H
        NVARS_BC_4D_H=solver_int_state%NVARS_BC_4D_H
        NVARS_BC_2D_V=solver_int_state%NVARS_BC_2D_V
        NVARS_BC_3D_V=solver_int_state%NVARS_BC_3D_V
!
!-----------------------------------------------------------------------
!***  The parents send and the children receive boundary data in the
!***  Parent-Child coupler.  Insert the Bundle of nest boundary update
!***  variables into the Domain export state so it will be available
!***  to the coupler.  It will be transferred to the P-C coupler's
!***  import state in subroutine PARENT_CHILD_COUPLER_SETUP.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert Bundle of Nest BC Vbls into Domain Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_StateAddReplace(            EXP_STATE                        &  !<-- The Domain export state
                          ,(/domain_int_state%BUNDLE_NESTBC/)  &  !<-- Insert the nest BC Bundle into the state
                          ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Insert the total number of model layers for all H-pt and V-pt
!***  boundary variables into the Domain component export state.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert Total Lyr Count of BC Vbls into Domain Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Domain export state
                              ,name ='NLEV_H'                           &  !<-- Name of the attribute to add
                              ,value=NLEV_H                             &  !<-- # of model layers for all BC H-pt variables
                              ,rc   =RC)
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Domain export state
                              ,name ='NLEV_V'                           &  !<-- Name of the attribute to add
                              ,value=NLEV_V                             &  !<-- # of model layers for all BC V-pt variables
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Also insert the # of H-pt and V-pt nest boundary variables
!***  into the Domain export state so they can be transferred to
!***  the Parent-Child coupler where they are needed.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert # of H-pt Nest Bndry Vbls into the Domain Exp State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Domain export state
                              ,name ='NVARS_BC_2D_H'                    &  !<-- Name of the attribute to add
                              ,value=NVARS_BC_2D_H                      &  !<-- # of domain's 2-D H-pt boundary variables
                              ,rc   =RC)
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Domain export state
                              ,name ='NVARS_BC_3D_H'                    &  !<-- Name of the attribute to add
                              ,value=NVARS_BC_3D_H                      &  !<-- # of domain's 3-D H-pt boundary variables
                              ,rc   =RC)
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Domain export state
                              ,name ='NVARS_BC_4D_H'                    &  !<-- Name of the attribute to add
                              ,value=NVARS_BC_4D_H                      &  !<-- # of domain's 4-D H-pt boundary variables
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert # of V-pt Nest Bndry Vbls into the Domain Exp State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Domain export state
                              ,name ='NVARS_BC_2D_V'                    &  !<-- Name of the attribute to extract
                              ,value=NVARS_BC_2D_V                      &  !<-- # of domain's 2-D V-pt boundary variables
                              ,rc   =RC)
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Domain export state
                              ,name ='NVARS_BC_3D_V'                    &  !<-- Name of the attribute to extract
                              ,value=NVARS_BC_3D_V                      &  !<-- # of domain's 3-D V-pt boundary variables
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        IF(NVARS_BC_4D_H==0)THEN
          ALLOCATE(solver_int_state%LBND_4D(1:1))
          ALLOCATE(solver_int_state%UBND_4D(1:1))
!
        ELSEIF(NVARS_BC_4D_H>0)THEN
          ALLOCATE(solver_int_state%LBND_4D(1:NVARS_BC_4D_H))
          ALLOCATE(solver_int_state%UBND_4D(1:NVARS_BC_4D_H))
          DO NV=1,NVARS_BC_4D_H
            LB=LBOUND(solver_int_state%BND_VARS_H%VAR_4D(NV)%FULL_VAR,4)
            solver_int_state%LBND_4D(NV)=LB
            UB=UBOUND(solver_int_state%BND_VARS_H%VAR_4D(NV)%FULL_VAR,4)
            solver_int_state%UBND_4D(NV)=UB
          ENDDO
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Insert Lower Bnds of 4-D H-pt Nest Bndry Vbls into the Domain Exp State"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeSet(state    =EXP_STATE                    &  !<-- The Domain component export state
                                ,name     ='LBND_4D'                    &  !<-- Use this name inside the state
                                ,itemCount=NVARS_BC_4D_H                &  !<-- # of items in Attribute
                                ,valueList=solver_int_state%LBND_4D     &  !<-- Lower bnds of 4-D H-pt boundary variablesmport state
                                ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Insert Upper Bnds of 4-D H-pt Nest Bndry Vbls into the Domain Exp State"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeSet(state    =EXP_STATE                    &  !<-- The Domain component export state
                                ,name     ='UBND_4D'                    &  !<-- Use this name inside the state
                                ,itemCount=NVARS_BC_4D_H                &  !<-- # of items in Attribute
                                ,valueList=solver_int_state%UBND_4D     &  !<-- Upper bnds of 4-D H-pt boundary variables
                                ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        ENDIF
!
!-----------------------------------------------------------------------
!***  Fill the Solver's boundary objects with values from the full
!***  arrays in the Solver internal state.
!-----------------------------------------------------------------------
!
        IF(.NOT.solver_int_state%RESTART)THEN
!
          CALL BOUNDARY_INIT(ITS,ITE,JTS,JTE,LM                         &
                            ,IMS,IME,JMS,JME                            &
                            ,IDS,IDE,JDS,JDE                            &
                            ,solver_int_state%LNSH                      &
                            ,solver_int_state%LNSV                      &
!
                            ,solver_int_state%NVARS_BC_2D_H             &
                            ,solver_int_state%NVARS_BC_3D_H             &
                            ,solver_int_state%NVARS_BC_4D_H             &
                            ,solver_int_state%NVARS_BC_2D_V             &
                            ,solver_int_state%NVARS_BC_3D_V             &
!
                            ,solver_int_state%BND_VARS_H                &
                            ,solver_int_state%BND_VARS_V                &
!
                              )
        ENDIF
!
!-----------------------------------------------------------------------
!***  The restart output file must contain the boundary array data
!***  .... BUT:
!***   (1) They must be passed from the Solver to the Write
!***       component and since they are not on the ESMF Grid
!***       they must be passed as 1-D Attributes.
!***   (2) We do not want to waste clocktime inserting these
!***       BC winds into the 1-D arrays every timestep when
!***       they are only needed at restart output times
!***       so we must inform the Solver when to fill those
!***       arrays.
!
!***  The 1-D arrays are placed into the Write component's import
!***  state in SAVE_BC_DATA called during SOLVER_RUN.  They are
!***  unloaded in WRT_RUN and sent to the lead forecast task to
!***  assemble into a full-domain 1-D datastring that can be sent
!***  to the lead write task for insertion into the restart file.
!-----------------------------------------------------------------------
!
        solver_int_state%NSTEPS_BC_RESTART=NINT((solver_int_state%MINUTES_RESTART*60) &  !<-- Timestep frequency for BC data insertion into
                                               /solver_int_state%DT)                     !    1-D local datastrings
!
        LNSH=solver_int_state%LNSH
        LNSV=solver_int_state%LNSV
!
!       IF(JTS==JDS)THEN                                                      !<-- South boundary tasks
        solver_int_state%NUM_WORDS_BC_SOUTH=(solver_int_state%NLEV_H*LNSH  &
                                            +solver_int_state%NLEV_V*LNSV) &
                                            *2*(ITE-ITS+1)
        ALLOCATE(solver_int_state%RST_BC_DATA_SOUTH(1:solver_int_state%NUM_WORDS_BC_SOUTH))
        DO N=1,solver_int_state%NUM_WORDS_BC_SOUTH
          solver_int_state%RST_BC_DATA_SOUTH(N)=0.
        ENDDO
!       ENDIF
!
!       IF(JTE==JDE)THEN                                                      !<-- North boundary tasks
        solver_int_state%NUM_WORDS_BC_NORTH=(solver_int_state%NLEV_H*LNSH  &
                                            +solver_int_state%NLEV_V*LNSV) &
                                            *2*(ITE-ITS+1)
        ALLOCATE(solver_int_state%RST_BC_DATA_NORTH(1:solver_int_state%NUM_WORDS_BC_NORTH))
        DO N=1,solver_int_state%NUM_WORDS_BC_NORTH
          solver_int_state%RST_BC_DATA_NORTH(N)=0.
        ENDDO
!       ENDIF
!
!       IF(ITS==IDS)THEN                                                      !<-- West boundary tasks
        solver_int_state%NUM_WORDS_BC_WEST=(solver_int_state%NLEV_H*LNSH   &
                                           +solver_int_state%NLEV_V*LNSV)  &
                                           *2*(JTE-JTS+1)
        ALLOCATE(solver_int_state%RST_BC_DATA_WEST(1:solver_int_state%NUM_WORDS_BC_WEST))
        DO N=1,solver_int_state%NUM_WORDS_BC_WEST
          solver_int_state%RST_BC_DATA_WEST(N)=0.
        ENDDO
!       ENDIF
!
!       IF(ITE==IDE)THEN                                                      !<-- East boundary tasks
        solver_int_state%NUM_WORDS_BC_EAST=(solver_int_state%NLEV_H*LNSH   &
                                           +solver_int_state%NLEV_V*LNSV)  &
                                           *2*(JTE-JTS+1)
        ALLOCATE(solver_int_state%RST_BC_DATA_EAST(1:solver_int_state%NUM_WORDS_BC_EAST))
        DO N=1,solver_int_state%NUM_WORDS_BC_EAST
          solver_int_state%RST_BC_DATA_EAST(N)=0.
        ENDDO
!       ENDIF
!
!-----------------------------------------------------------------------
!
      ENDIF bc_variables
!
!-----------------------------------------------------------------------
!***  All compute tasks will now insert history and restart pointers
!***  from the Solver internal state into the Write component's
!***  import state.  This makes the output variables available to
!***  the Write component.
!-----------------------------------------------------------------------
!
      IF(MYPE<domain_int_state%NUM_PES_FCST.AND.QUILTING)THEN              !<-- Select the compute tasks.
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Get Write Import State from Solver Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_StateGet(state      =domain_int_state%EXP_STATE_SOLVER &  !<-- The Solver export state
                          ,itemName   ='Write Import State'              &  !<-- Name of the state to get from Solver export state
                          ,nestedState=IMP_STATE_WRITE                   &  !<-- Extract Write component import state from Solver export
                          ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL POINT_OUTPUT(GRID_DOMAIN,SOLVER_INT_STATE,IMP_STATE_WRITE)
!
      ENDIF
!
!-----------------------------------------------------------------------
!***  If quilting was selected for the generation of output,
!***  execute the Initialize step of the Write component(s).
!***  This must be done after the initialization of the
!***  Solver component because that component's internal
!***  state contains the history output variables.
!***  Pointers to those variables are set during the INITIALIZE
!***  steps and are then loaded into the Write components'
!***  import states which themselves reside in the Solver
!***  export states.
!-----------------------------------------------------------------------
!
      I_AM_A_FCST_TASK=.TRUE.
!
      IF(QUILTING)THEN
!
        CALL WRITE_INIT(DOMAIN_GRID_COMP                                &
                       ,DOMAIN_INT_STATE                                &
                       ,IMP_STATE                                       &
                       ,CLOCK_DOMAIN,MYPE)
!
        IF(MYPE>=domain_int_state%NUM_PES_FCST)THEN
!
          I_AM_A_FCST_TASK=.FALSE.
!
        ENDIF
!
      ENDIF
!
!-----------------------------------------------------------------------
!***  Add some key variables to the Domain export state.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Add Fcst-or-Write Task Flag to the Domain Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(state=EXP_STATE                            &  !<-- The Domain component export state
                            ,name ='Fcst-or-Write Flag'                 &  !<-- Use this name inside the state
                            ,value=I_AM_A_FCST_TASK                     &  !<-- The logical being inserted into the import state
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Add NUM_PES_FCST to the Domain Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(state=EXP_STATE                            &  !<-- The Domain component export state
                            ,name ='NUM_PES_FCST'                       &  !<-- Use this name inside the state
                            ,value=domain_int_state%NUM_PES_FCST        &  !<-- The value being set
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Add Fcst Task Intracomm to Domain Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(state=EXP_STATE                            &  !<-- The Domain component export state
                            ,name ='Comm Fcst Tasks'                    &  !<-- Use this name inside the state
                            ,value=MPI_COMM_COMP                        &  !<-- This domain's intracomm for fcst tasks
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Now extract number of children on this Domain component's domain.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
      MESSAGE_CHECK="Extract Number of Children from Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
      CALL ESMF_AttributeGet(state       =IMP_STATE                     &  !<-- The Domain import state
                            ,name        ='NUM_CHILDREN'                &  !<-- Name of the attribute to extract
                            ,value       =NUM_CHILDREN                  &  !<-- Put the Attribute here
                            ,defaultValue=0                             &  !<-- The default value
                            ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
!-----------------------------------------------------------------------
!***  If this domain has children then retrieve their domain IDs.
!***  They are needed for various aspects of initialization and
!***  integration.
!-----------------------------------------------------------------------
!
      IF(NUM_CHILDREN>0)THEN
!
        ALLOCATE(MY_CHILDREN_ID(1:NUM_CHILDREN),stat=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract Children's IDs from Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state    =IMP_STATE                      &  !<-- The Domain import state
                              ,name     ='CHILD_IDs'                    &  !<-- Name of the attribute to extract
                              ,itemCount=NUM_CHILDREN                   &  !<-- # of items in the Attribute
                              ,valueList=MY_CHILDREN_ID                 &  !<-- Put the Attribute here
                              ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF
!
!-----------------------------------------------------------------------
!***  With the Solver internal state available set the location of
!***  the SW corner of this domain if it is a nest.  This provides 
!***  the corner location values to the Write (output) components
!***  and therefore must precede the creation of those components
!***  which immediately follows.
!
!***  If the input file was generated by NPS or this is a restarted
!***  run then the value of the SW corner has already been placed into
!***  the Solver internal state.  If this is a free forecast without
!***  a pre-generated input file then read the SW corner location
!***  from the configure file.
!-----------------------------------------------------------------------
!
      IF(I_AM_A_FCST_TASK.AND.I_AM_A_NEST)THEN
!
        IF(.NOT.INPUT_READY.AND..NOT.RESTARTED_RUN)THEN                      !<-- If so then must get values from configure file.
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Domain Init: Child Gets SW Corner Point from Config File"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)            &  !<-- The nest domain's config object
                                      ,value =I_PAR_STA                   &  !<-- The variable filled (parent I of nest SW corner)
                                      ,label ='i_parent_start:'           &  !<-- Give this label's value to the previous variable
                                      ,rc    =RC)
!
          CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)            &  !<-- The nest domain's config object
                                      ,value =J_PAR_STA                   &  !<-- The variable filled (parent J of nest SW corner)
                                      ,label ='j_parent_start:'           &  !<-- Give this label's value to the previous variable
                                      ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          solver_int_state%I_PAR_STA = I_PAR_STA   
          solver_int_state%J_PAR_STA = J_PAR_STA  
!
          NEXT_MOVE_TIMESTEP=-999
          solver_int_state%NMTS=NEXT_MOVE_TIMESTEP 
          solver_int_state%LAST_STEP_MOVED=0
!
        ELSE                                                                 !<-- If so values were read from input or restart file.
!
          I_PAR_STA=solver_int_state%I_PAR_STA
          J_PAR_STA=solver_int_state%J_PAR_STA
          NEXT_MOVE_TIMESTEP=solver_int_state%NMTS
!
          TRACKER_IFIX=solver_int_state%TRACKER_IFIX
          TRACKER_JFIX=solver_int_state%TRACKER_JFIX
          STORM_CENTER(1)=solver_int_state%TRACKER_IFIX
          STORM_CENTER(2)=solver_int_state%TRACKER_JFIX
!
          IF(INPUT_READY.AND..NOT.RESTARTED_RUN)THEN
            solver_int_state%LAST_STEP_MOVED=0
!
          ENDIF
!
          LAST_STEP_MOVED=solver_int_state%LAST_STEP_MOVED
!
        ENDIF
!
!-----------------------------------------------------------------------
!***  Transfer the domain's SW corner and next move timestep to the
!***  Domain export state.  These values are dummies if not relevant.
!***  The values are obtained directly from the Solver internal state
!***  if this is a restarted run since they were read from the restart 
!***  file in that case.
!-----------------------------------------------------------------------

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert SW corner of Nest into Domain Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Domain export state
                              ,name ='I_PAR_STA'                        &  !<-- Name of the attribute to extract
                              ,value=I_PAR_STA                          &  !<-- Put the Attribute here
                              ,rc   =RC)
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Domain export state
                              ,name ='J_PAR_STA'                        &  !<-- Name of the attribute to extract
                              ,value=J_PAR_STA                          &  !<-- Put the Attribute here
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert NTRACK flag into the Domain Export State."
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        IF(ASSOCIATED(solver_int_state%NTRACK))THEN
          NTRACK=solver_int_state%NTRACK
        ELSE
          NTRACK=-99
        ENDIF
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Domain export state
                              ,name ='NTRACK'                           &  !<-- Name of the attribute to extract
                              ,value=NTRACK                             &  !<-- Put the Attribute here
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert Next Move Timestep into Domain Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Domain export state
                              ,name ='NEXT_MOVE_TIMESTEP'               &  !<-- Name of the attribute to extract
                              ,value=NEXT_MOVE_TIMESTEP                 &  !<-- Put the Attribute here
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert Last Move Timestep into Domain Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Domain export state
                              ,name ='LAST_STEP_MOVED'                  &  !<-- Name of the attribute to extract
                              ,value=LAST_STEP_MOVED                    &  !<-- Put the Attribute here
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert Storm Center into Domain Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Domain export state
                              ,name ='Storm Center'                     &  !<-- Name of the attribute to extract
                              ,itemCount=2                              &  !<-- Number of items in the array
                              ,valueList=STORM_CENTER                   &  !<-- Put the Attribute here
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
!***  If the current component/domain is the parent of nests then:
!
!***  (1) Extract the arrays from the Solver export state that
!***      is required for the children's boundaries and insert
!***      them into the Domain export state since ultimately 
!***      they must be available to the parent in the
!***      Parent-Child Coupler.
!
!***  (2) Check to see if the children have input data ready for them.
!***      If not, do simple nearest neighbor and bilinear interpolation
!***      from the parent's grid to the children's.  Write out that
!***      interpolated data into files that are waiting for the children 
!***      when they recursively execute DOMAIN_INITIALIZE themselves.
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      I_AM_A_PARENT=.FALSE.
!
!-----------------------------------------------------------------------
!
      fcst_tasks_init: IF(I_AM_A_FCST_TASK)THEN
!
!-----------------------------------------------------------------------
!***  Extract fundamental meteorological variables from the Solver
!***  export state.  Insert them into the Domain export state so that 
!***  NMM_INITIALIZE can take them and send them to the Parent-Child
!***  coupler.  Only the forecast tasks participate in doing this
!***  since the Write tasks never loaded data into the Solver export
!***  state.
!-----------------------------------------------------------------------
!
        CALL INTERNAL_DATA_TO_DOMAIN(domain_int_state%EXP_STATE_SOLVER  &  !<-- The Solver export state
                                    ,EXP_STATE                          &  !<-- The Domain export state
                                    ,NLAYRS )                              !<-- # of model layers
!
!-----------------------------------------------------------------------
!***  If there are moving nests then the Parent-Child coupler will
!***  need pointers to all the required Solver arrays that must be
!***  updated after a nest moves.  Create ESMF Bundles to hold those
!***  pointers then insert the designated pointers from the internal
!***  state into the Bundles (one for H-pt variables and one for
!***  V-pt variables).  These Bundles simply remain empty if there
!***  are no moving nests.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Create the empty Move Bundles.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Create the Move Bundles"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        BUNDLE_NAME='Move_Bundle H'
!
        domain_int_state%MOVE_BUNDLE_H=ESMF_FieldBundleCreate(name=BUNDLE_NAME  &  !<-- The H-pt Bundle's name
                                                             ,rc  =RC)
!
        BUNDLE_NAME='Move_Bundle V'
!
        domain_int_state%MOVE_BUNDLE_V=ESMF_FieldBundleCreate(name=BUNDLE_NAME  &  !<-- The V-pt Bundle's name
                                                             ,rc  =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        NUM_FIELDS_MOVE_2D_H_I=0
        NUM_FIELDS_MOVE_2D_H_R=0
        NUM_FIELDS_MOVE_3D_H=0
        NUM_FIELDS_MOVE_2D_V=0
        NUM_FIELDS_MOVE_3D_V=0
!
        NUM_LEVELS_MOVE_3D_H=0
        NUM_LEVELS_MOVE_3D_V=0
!
!-----------------------------------------------------------------------
!***  Fill the Bundles with the variables to be shifted after nests
!***  move.  All moving domains and parents of moving domains must
!***  fill the Move Bundles.
!-----------------------------------------------------------------------
!
        CALL_BUILD_MOVE_BUNDLE=.FALSE.
!
!-----------------------------------------------------------------------
!***  Does this domain have any moving children?
!-----------------------------------------------------------------------
!
        IF(NUM_CHILDREN>0)THEN
!
          child_loop: DO N=1,NUM_CHILDREN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Extract the Child's Flag Indicating Movability"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_ConfigGetAttribute(config=CF(MY_CHILDREN_ID(N))   &  !<-- The child's config object
                                        ,value =DOMAIN_MOVES            &  !<-- The variable filled (will the child move?)
                                        ,label ='my_domain_moves:'      &  !<-- Give this label's value to the previous variable
                                        ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            IF(DOMAIN_MOVES)THEN                                             !<-- If true then child N moves.
!
              CALL_BUILD_MOVE_BUNDLE=.TRUE.
              EXIT child_loop
!
            ENDIF
!
          ENDDO child_loop
!
        ENDIF
!
!-----------------------------------------------------------------------
!
        IF(MY_DOMAIN_MOVES)THEN
!
          CALL_BUILD_MOVE_BUNDLE=.TRUE.
!
        ENDIF
!
!-----------------------------------------------------------------------
!***  Now build the Move Bundles if any children move or if the parent
!***  itself moves.
!-----------------------------------------------------------------------
!
        IF(CALL_BUILD_MOVE_BUNDLE)THEN
!
          UBOUND_VARS=SIZE(solver_int_state%VARS)
!
          CALL BUILD_MOVE_BUNDLE(GRID_DOMAIN                            &  !<-- Add Solver variables to H and V Move Bundles
                                ,UBOUND_VARS                            &
                                ,solver_int_state%VARS                  &
                                ,domain_int_state%MOVE_BUNDLE_H         &
                                ,NUM_FIELDS_MOVE_2D_H_I                 &
                                ,NUM_FIELDS_MOVE_2D_H_R                 &
                                ,NUM_FIELDS_MOVE_3D_H                   &
                                ,NUM_LEVELS_MOVE_3D_H                   &
                                ,domain_int_state%MOVE_BUNDLE_V         &
                                ,NUM_FIELDS_MOVE_2D_V                   &
                                ,NUM_FIELDS_MOVE_3D_V                   &
                                ,NUM_LEVELS_MOVE_3D_V                   &
                                  ) 
!
!         CALL ESMF_FieldBundlePrint(domain_int_state%MOVE_BUNDLE_H)
!         CALL ESMF_FieldBundlePrint(domain_int_state%MOVE_BUNDLE_V)
!
        ENDIF
!
!-----------------------------------------------------------------------
!***  Since the parents will also update some of the moving nests'
!***  interior points, the Bundles will be moved into the Parent-Child
!***  coupler import state in subroutine PARENT_CHILD_COUPLER_SETUP.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert Move Bundles into Domain Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_StateAddReplace(            EXP_STATE                        &  !<-- The Domain export state
                          ,(/domain_int_state%MOVE_BUNDLE_H/)  &  !<-- Insert H-point MOVE_BUNDLE into the state
                          ,rc         =RC)
!
        CALL ESMF_StateAddReplace(            EXP_STATE                        &  !<-- The Domain export state
                          ,(/domain_int_state%MOVE_BUNDLE_V/)  &  !<-- Insert V-point MOVE_BUNDLE into the state
                          ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  What is the ratio of the uppermost parent's grid increment to
!***  this moving nest's?  That ratio is needed as part of the unique
!***  name of external files that contain nest-resolution data spanning
!***  the upper parent's domain that the nest reads directly.
!
!***  At the same time the moving nests must know the lateral
!***  dimensions of the uppermost parent domain so that they can
!***  properly read those external files which span that domain.
!-----------------------------------------------------------------------
!
        i_move: IF(MY_DOMAIN_MOVES)THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Extract Moving Child's Sfc File Ratio"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)          &  !<-- This domain's configure object
                                      ,value =SFC_FILE_RATIO            &  !<-- Ratio of upper parent's grid increment to this nest's
                                      ,label ='ratio_sfc_files:'        &  !<-- The variable read from the configure file
                                      ,rc    =RC)
!
          domain_int_state%SFC_FILE_RATIO=SFC_FILE_RATIO
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Extract Uppermost Parent Dimensions"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_ConfigGetAttribute(config=CF(1)                     &  !<-- The uppermost domain's configure object
                                      ,value =IM_1                      &  !<-- # of that domain's gridpoints in I direction 
                                      ,label ='im:'                     &  !<-- The variable read from the configure file
                                      ,rc    =RC)
!
          CALL ESMF_ConfigGetAttribute(config=CF(1)                     &  !<-- The uppermost domain's configure object
                                      ,value =JM_1                      &  !<-- # of that domain's gridpoints in J direction 
                                      ,label ='jm:'                     &  !<-- The variable read from the configure file
                                      ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Compute some values needed by moving nests for reading their
!***  full-resolution data that spans the uppermost parent.
!-----------------------------------------------------------------------
!
          D_ONE=1.
          D_180=180.
          PI_LOC=DACOS(-D_ONE)
          D2R=PI_LOC/D_180
!
          TPH0_1=TPH0D_1*D2R                                               !<-- Central geo lat of domain (radians, positive north)
          TLM0_1=TLM0D_1*D2R                                               !<-- Central geo lon of domain (radians, positive east)
          WB_1=WBD_1*D2R                                                   !<-- Rotated lon of west boundary (radians, positive east)
          SB_1=SBD_1*D2R                                                   !<-- Rotated lat of south boundary (radians, positive north)
!
          DPH_1=-2.*SB_1/(JM_1-1)                                          !<-- Uppermost parent's grid increment in J (radians)
          DLM_1=-2.*WB_1/(IM_1-1)                                          !<-- Uppermost parent's grid increment in I (radians)
!
          RECIP_DPH_1=1./DPH_1
          RECIP_DLM_1=1./DLM_1
!
!-----------------------------------------------------------------------
!***  When a parent sends interior update data to a moving child then
!***  at coastlines a parent may send data valid for water/land to a
!***  point on the child that is land/water on the child's sea mask.
!***  The user sets a configure flag to select one of two options to
!***  handle this situation.  In the general case the value 'nearest'
!***  is selected.  Then when a conflict point is encountered the 
!***  given child task searches on its subdomain for the nearest point
!***  to the conflict point that has the same sfc type (water or land)
!***  and uses that point's sfc values for the conflict point.  If no
!***  other point with the same sfc type can be found on the subdomain
!***  then a dummy value is assigned.  Note that this can lead to 
!***  different answers when different task layouts are used.  The
!***  other choice is 'dummy'.  When the user chooses that option then
!***  children automatically always set values at conflict points to
!***  dummy values.  Points on the earth will thus always have dummy
!***  values with 'dummy' whereas with 'nearest' the values at conflict
!***  points will likely have appropriate values during most of the
!***  time those locations lie within the moving child domain.  If
!***  identical answers are required for different task layouts then
!***  'dummy' must be used.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Extract SFC_CONFLICT from Config File"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)          &  !<-- This domain's configure object
                                      ,value =SFC_CONFLICT              &  !<-- Flag for handling parent-child sfc-type conflicts
                                      ,label ='sfc_conflict:'           &  !<-- The variable read from the configure file
                                      ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          domain_int_state%SFC_CONFLICT=SFC_CONFLICT
!
          IF(SFC_CONFLICT=='nearest')THEN
!
!-----------------------------------------------------------------------
!***  Generate the I,J increments needed to search for neighboring
!***  points to fix values at moving nest points where land points
!***  receive water point values from the parent and vice versa.
!***  Create empty objects for sorting distances between points on
!***  moving nests for patching mismatches between parent and child
!***  water and land points in 2-way exchange.
!-----------------------------------------------------------------------
!
            IF(.NOT.ASSOCIATED(SMALLX))THEN
              ALLOCATE(SMALLX(1:NUM_DOMAINS))
              ALLOCATE(LARGEX(1:NUM_DOMAINS))
            ENDIF
!
            CALL SEARCH_INIT
!
          ENDIF
!
!-----------------------------------------------------------------------
!***  Extract the nest grid increments for later use.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Extract DPHD and DLMD from Solver Export State"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeGet(state=domain_int_state%EXP_STATE_SOLVER &  !<-- The Solver export state
                                ,name ='DPHD'                            &  !<-- Name of the Attribute to extract
                                ,value=DPHD                              &  !<-- Angular grid increment in X (degrees)
                                ,rc   =RC)
!
          CALL ESMF_AttributeGet(state=domain_int_state%EXP_STATE_SOLVER &  !<-- The Solver export state
                                ,name ='DLMD'                            &  !<-- Name of the Attribute to extract
                                ,value=DLMD                              &  !<-- Angular grid increment in Y (degrees)
                                ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          DEG_TO_RAD=PI_LOC/180.
          DPH=DPHD*DEG_TO_RAD
          DLM=DLMD*DEG_TO_RAD
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Extract TPH0D and TLM0D from Solver Export State"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeGet(state=domain_int_state%EXP_STATE_SOLVER  &  !<-- The Solver export state
                                ,name ='TPH0D'                            &  !<-- Name of the Attribute to extract
                                ,value=TPH0D                              &  !<-- Central latitude (degrees) of rotated system
                                ,rc   =RC)
!
          CALL ESMF_AttributeGet(state=domain_int_state%EXP_STATE_SOLVER  &  !<-- The Solver export state
                                ,name ='TLM0D'                            &  !<-- Name of the Attribute to extract
                                ,value=TLM0D                              &  !<-- Central longitude (degrees) of rotated system
                                ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          TPH0=TPH0D*DEG_TO_RAD
          TLM0=TLM0D*DEG_TO_RAD
!
!-----------------------------------------------------------------------
!***  Moving nests explicitly compute the lat/lon in their parent
!***  update regions following each move as well as the HDAC variables
!***  that are directly dependent upon the lat/lon.  The Smagorinsky
!***  constant supplied by the user in the configure file is needed
!***  for the HDAC computation so extract it now.  The fundamental
!***  dynamical timestep length is also needed for this purpose as
!***  are grid constants.
!***  (Note:  DOMAIN_RUN extracts the timestep during the integration
!***          but that is after its sign may have changed if digital
!***          filtering is being used.)
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Extract Moving Child's Smagorinsky Constant"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)          &  !<-- This domain's configure object
                                      ,value =SMAG2                     &  !<-- Smagorinsky constant
                                      ,label ='smag2:'                  &  !<-- The variable read from the configure file
                                      ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Extract Moving Child's WCOR Constant"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)          &  !<-- This domain's configure object
                                      ,value =WCOR                      &
                                      ,label ='wcor:'                   &  !<-- The variable read from the configure file
                                      ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Extract CODAMP from Configure File"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)          &  !<-- This domain's configure object
                                      ,value =CODAMP                    &  !<-- Divergence damping coefficient
                                      ,label ='codamp:'                 &  !<-- The variable read from the configure file
                                      ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Extract Moving Child's Fundamental Timestep"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)          &  !<-- This domain's configure object
                                      ,value =DT_INT                    &  !<-- Integer part of time step.
                                      ,label ='dt_int:'                 &  !<-- The variable read from the configure file
                                      ,rc    =RC)
!
          CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)          &  !<-- This domain's configure object
                                      ,value =DT_NUM                    &  !<-- Numerator of fractional part of time step.
                                      ,label ='dt_num:'                 &  !<-- The variable read from the configure file
                                      ,rc    =RC)
!
          CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)          &  !<-- This domain's configure object
                                      ,value =DT_DEN                    &  !<-- Denominator of fractional part of time step.
                                      ,label ='dt_den:'                 &  !<-- The variable read from the configure file
                                      ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          DT_REAL=REAL(DT_INT)+REAL(DT_NUM)/REAL(DT_DEN)
!
          ACDT  =SMAG2*SMAG2*DT_REAL
          CDDAMP=CODAMP*DT_REAL
!
!-----------------------------------------------------------------------
!***  Due to the nature of the B-grid and the computations within 
!***  the NMM-B, locations corresponding to a minimum of the outer
!***  two rows of the pre-move footprint of the nest domain cannot
!***  use intra- or inter-task updates and instead must be updated
!***  by the parent.  Read in configure variables that specify the
!***  number of rows on each side of the nest's pre-move footprint 
!***  for which the parent will provide update data after the nest
!***  moves.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Extract # of Rows Parent Will Update"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)          &  !<-- This domain's configure object
                                      ,value =NROWS_P_UPD_W             &  !<-- # of rows parent will update on west bndry
                                      ,label ='nrows_p_upd_w:'          &  !<-- The variable read from the configure file
                                      ,rc    =RC)
!
          CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)          &  !<-- This domain's configure object
                                      ,value =NROWS_P_UPD_E             &  !<-- # of rows parent will update on east bndry
                                      ,label ='nrows_p_upd_e:'          &  !<-- The variable read from the configure file
                                      ,rc    =RC)
!
          CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)          &  !<-- This domain's configure object
                                      ,value =NROWS_P_UPD_S             &  !<-- # of rows parent will update on south bndry
                                      ,label ='nrows_p_upd_s:'          &  !<-- The variable read from the configure file
                                      ,rc    =RC)
!
          CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)          &  !<-- This domain's configure object
                                      ,value =NROWS_P_UPD_N             &  !<-- # of rows parent will update on north bndry
                                      ,label ='nrows_p_upd_n:'          &  !<-- The variable read from the configure file
                                      ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Due to the complexities of the parent's updating of the child
!***  for nest motion that is due east, west, south, or north 
!***  task subdomains cannot be too thin.  Stop the run if they are.
!-----------------------------------------------------------------------
!
          IF(IHALO>ITE-ITS+1-NROWS_P_UPD_W.OR.                          &
             IHALO>ITE-ITS+1-NROWS_P_UPD_E.OR.                          &
             JHALO>JTE-JTS+1-NROWS_P_UPD_S.OR.                          &
             JHALO>JTE-JTS+1-NROWS_P_UPD_N )THEN
            WRITE(0,*)' Task subdomains cannot be narrower than '
            WRITE(0,*)' the width of the halo plus the width of '
            WRITE(0,*)' the parent update region on the outer '
            WRITE(0,*)' edge of a nest pre-move footprint.'
            WRITE(0,11111)IHALO,JHALO
            WRITE(0,*)' The width of the parent update region on the '
            WRITE(0,*)' south, north, west, and east side of a moving '
            WRITE(0,*)' nest pre-move footprint are:'
            WRITE(0,11112)NROWS_P_UPD_S,NROWS_P_UPD_N                   &
                         ,NROWS_P_UPD_W,NROWS_P_UPD_E
            WRITE(0,11113)ITE-ITS+1,JTE-JTS+1
            WRITE(0,*)' The user must reset the domain decomposition.'
            WRITE(0,*)' Aborting!!'
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
11111       FORMAT(' The halo width in I is ',I2,' and in J is ',I2)
11112       FORMAT(4(1X,I2))
11113       FORMAT(' The subdomain widths in I and J are ',I4,1X,I4)
          ENDIF
!
!-----------------------------------------------------------------------
!***  Initialize the handle used in ISSends of intertask data when
!***  a moving domain shifts.
!-----------------------------------------------------------------------
!
          DO N=1,9
            domain_int_state%HANDLE_SEND_INTER_INT(N) =MPI_REQUEST_NULL
            domain_int_state%HANDLE_SEND_INTER_REAL(N)=MPI_REQUEST_NULL
          ENDDO
!
!-----------------------------------------------------------------------
!
        ENDIF i_move
!
!-----------------------------------------------------------------------
!***  If there is 2-way exchange from the children to the parents
!***  then the Parent-Child coupler will need pointers to all the
!***  required Solver arrays that are updated on the parents by the
!***  children each parent timestep.  Create an ESMF Bundle to hold
!***  those pointers then insert the desgnated pointers from the
!***  Solver component's internal state into the Bundle.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Does this forecast use 2-way exchange?
!-----------------------------------------------------------------------
!
        NEST_MODE=' '
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract 2-way Flag From Nest Configure Files"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)            &  !<-- This domain's configure object
                                    ,value =NEST_MODE                   &  !<-- Is there 2-way exchange from child to parent?
                                    ,label ='nest_mode:'                &  !<-- The label in the configure file
                                    ,rc    =rc)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Create the 2-way exchange Bundle.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Create the 2-way Exchange Bundle"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        BUNDLE_NAME='Bundle_2way'
!
        domain_int_state%BUNDLE_2WAY=ESMF_FieldBundleCreate(name=BUNDLE_NAME  &  !<-- The 2-way Bundle's name
                                                           ,rc  =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Now build the 2-way Bundle if 2-way exchange has been specified
!***  in the configure files.
!-----------------------------------------------------------------------
!
        IF(NEST_MODE=='2-way')THEN
!
          UBOUND_VARS=SIZE(solver_int_state%VARS)
!
          CALL BUILD_2WAY_BUNDLE(GRID_DOMAIN                            &  !<-- Add Solver int state variables to 2-way Bundle
                                ,LM                                     &
                                ,UBOUND_VARS                            &
                                ,solver_int_state%VARS                  &
                                ,domain_int_state%BUNDLE_2WAY           &
                                     ) 
!
        ENDIF
!
!-----------------------------------------------------------------------
!***  The 2-way exchange itself takes place in the Parent-Child coupler 
!***  so insert the 2-way Bundle into the Domain component's export
!***  state in order to transfer it to the P-C coupler import state
!***  in subroutine PARENT_CHILD_COUPLER_SETUP.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert 2-way Bundle into Domain Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_StateAddReplace(            EXP_STATE                        &  !<-- The Domain export state
                          ,(/domain_int_state%BUNDLE_2WAY/)    &  !<-- Insert BUNDLE of 2-way vbls into the state
                          ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
      ENDIF fcst_tasks_init
!
!-----------------------------------------------------------------------
!
      child_init_block: IF(NUM_CHILDREN>0)THEN                             !<-- Only parents participate
!
!-----------------------------------------------------------------------
!
        I_AM_A_PARENT=.TRUE.
!
!-----------------------------------------------------------------------
!***  Initialize the children's data directly from the parent if
!***  there are no pre-processed input files ready for them.
!***  Files will be written for the children to read in as usual.
!***  Only parent tasks participate.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        child_init_loop: DO N=1,NUM_CHILDREN                               !<-- Loop through the children
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Extract Children's Input Flag from Config File"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_ConfigGetAttribute(config=CF(MY_CHILDREN_ID(N))     &  !<-- The config object
                                      ,value =INPUT_READY_MY_CHILD      &  !<-- Child's flag for existence of its input file
                                      ,label ='input_ready:'            &  !<-- Give this label's value to the previous variable
                                      ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          IF(.NOT.INPUT_READY_MY_CHILD)THEN                                    !<-- INPUT_READY=false -> This child has no input file
                                                                               !      so parent will generate input.
            CALL PARENT_TO_CHILD_INIT_NMM(MYPE                              &  !<-- This task's rank (in)
                                         ,CF                                &  !<-- Array of configure files (in)
                                         ,MY_DOMAIN_ID                      &  !<-- Each domain's ID (in)
                                         ,MY_CHILDREN_ID(N)                 &  !<-- The child's domain ID
                                         ,domain_int_state%SOLVER_GRID_COMP &  !<-- The parent's Solver Component (inout)
                                         ,COMM_MY_DOMAIN )                     !<-- Each domain's intracommunicator
!
          ENDIF
!
!-----------------------------------------------------------------------
!
        ENDDO child_init_loop
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Send the next move timestep of the moving children to the
!***  Parent-Child coupler.  If this is a restarted run then the
!***  values are were read from the restart file in Solver Init.
!-----------------------------------------------------------------------
!
        IF(I_AM_A_FCST_TASK)THEN
!
          NTIMESTEP_CHILD_MOVES=>solver_int_state%NTSCM
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Initialize Next Timestep Children Move in Domain Import State"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeSet(state    =IMP_STATE                    &  !<-- The Domain import state
                                ,name     ='NEXT_TIMESTEP_CHILD_MOVES'  &  !<-- Name of the attribute to insert
                                ,itemCount=NUM_DOMAINS_MAX              &  !<-- Number of items in the array
                                ,valueList=NTIMESTEP_CHILD_MOVES        &  !<-- Put the Attribute here
                                ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Insert Next Timestep Children Move into Domain Export State"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeSet(state    =EXP_STATE                    &  !<-- The Domain export state
                                ,name     ='NEXT_TIMESTEP_CHILD_MOVES'  &  !<-- Name of the attribute to extract
                                ,itemCount=NUM_DOMAINS_MAX              &  !<-- Number of items in the array
                                ,valueList=NTIMESTEP_CHILD_MOVES        &  !<-- Put the Attribute here
                                ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        ENDIF
!
!-----------------------------------------------------------------------
!
      ENDIF child_init_block
!
!-----------------------------------------------------------------------
!***  For moving nests there are external files with nest-resolution
!***  sfc data spanning the uppermost parent.  If the parent generated
!***  the nest's initial conditions from its own then replace the
!***  values in those nest arrays with data from the external files.
!-----------------------------------------------------------------------
!   
      IF(MY_DOMAIN_MOVES                                                &
            .AND.                                                       &
         I_AM_A_FCST_TASK                                               &
            .AND.                                                       &
         .NOT.INPUT_READY)THEN
!
        CALL RESET_SFC_VARS(domain_int_state%SFC_FILE_RATIO             &
                           ,solver_int_state%GLAT                       &
                           ,solver_int_state%GLON                       &
                           ,domain_int_state%MOVE_BUNDLE_H)
!
        CALL RESET_SFC_VARS(domain_int_state%SFC_FILE_RATIO             &
                           ,solver_int_state%GLAT                       &
                           ,solver_int_state%GLON                       &
                           ,domain_int_state%MOVE_BUNDLE_V)
!
!-----------------------------------------------------------------------
!***  Now the nest's sea mask array contains the nest-resolution 
!***  data from the external file.  That means the nest's sea mask
!***  is at nest resolution while other land/sea variables were
!***  simply interpolated from the parent domain so near coastlines
!***  some points in those variables will not agree with the nest's
!***  sea mask.  Therefore call the same routine that must be called 
!***  after every move of the nest during the integration that will
!***  force various land/water variables to agree with the nest's 
!***  sea mask.
!-----------------------------------------------------------------------
!
        FIELD_NAME='SM-move'
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract Seamask Field from Move Bundle H"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_FieldBundleGet(FIELDBUNDLE=domain_int_state%MOVE_BUNDLE_H  &  !<-- Bundle holding the H arrays for move updates
                                ,fieldName  =FIELD_NAME                      &  !<-- Name of the seamask Field in the Bundle
                                ,field      =HOLD_FIELD                      &  !<-- Field containing the seamask
                                ,rc         =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract Seamask Array from Field"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_FieldGet(field    =HOLD_FIELD                         &  !<-- Field N_FIELD in the Bundle
                          ,localDe  =0                                  &
                          ,farrayPtr=SEA_MASK                           &  !<-- Dummy 2-D array with Field's Real data
                          ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        ILO=LBOUND(SEA_MASK,1)
        IHI=UBOUND(SEA_MASK,1)
        JLO=LBOUND(SEA_MASK,2)
        JHI=UBOUND(SEA_MASK,2)
!
        NUM_FIELDS=NUM_FIELDS_MOVE_2D_H_I                               &
                  +NUM_FIELDS_MOVE_2D_H_R                               &
                  +NUM_FIELDS_MOVE_3D_H
!
        CALL FIX_SFC(domain_int_state%MOVE_BUNDLE_H                     &
                    ,NUM_FIELDS                                         &
                    ,SEA_MASK                                           &
                    ,ILO,IHI,JLO,JHI                                    &
                    ,ILO,IHI,JLO,JHI)
!
      ENDIF
!
!-----------------------------------------------------------------------
!***  Insert into the Domain export state the flag indicating if the
!***  current domain is a parent.  The Domain Driver wants to know this
!***  since most Parent-Child work can be ignored by domains with
!***  no children.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Add Parent/Not-a-Parent Flag to the Domain Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(state=EXP_STATE                            &  !<-- The Domain component export state
                            ,name ='I-Am-A-Parent Flag'                 &  !<-- Use this name inside the state
                            ,value=I_AM_A_PARENT                        &  !<-- The logical being inserted into the import state
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      domain_int_state%I_AM_A_PARENT=I_AM_A_PARENT
!
!-----------------------------------------------------------------------
!***  Initialize some variables used in NMM_INTEGRATE.
!-----------------------------------------------------------------------
!
      domain_int_state%FIRST_PASS=.TRUE.                                   !<-- Note the first time NMM_INTEGRATE is entered.
      domain_int_state%TS_INITIALIZED=.FALSE.                              !<-- Note whether time series variables are initialized.
      domain_int_state%KOUNT_TIMESTEPS=0                                   !<-- Timestep counter
!
      domain_int_state%RECV_ALL_CHILD_DATA=.FALSE.                         !<-- Parent has recvd 2-way data from all children
      domain_int_state%ALLCLEAR_FROM_PARENT=.FALSE.                        !<-- Child told that parent has recvd all 2-way data
!
!-----------------------------------------------------------------------
!
      CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)              &  !<-- This domain's configure object
                                  ,value =FILT_METHOD                   &  !<-- The filter method
                                  ,label ='filter_method:'              &  !<-- The variable read from the configure file
                                  ,rc    =RC)
!
!-----------------------------------------------------------------------
!
      IF(I_AM_A_FCST_TASK)THEN
!
        IF(NUM_CHILDREN>0)THEN
          ALLOCATE(CHILD_ACTIVE(1:NUM_CHILDREN))
        ENDIF
!
!-----------------------------------------------------------------------
!***  Create and fill the Filter Bundles.  Since the current task
!***  might lie on more than one domain if the user selects 2-way
!***  nesting there needs to be a unique bundle for each of those
!***  domains.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract Filter Method beofre Creating Filter Bundles"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)            &  !<-- This domain's configure object
                                    ,value =FILT_METHOD                 &  !<-- The filter method
                                    ,label ='filter_method:'            &  !<-- The variable read from the configure file
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        IF (FILT_METHOD > 0) THEN
!
          domain_int_state%FILT_BUNDLE_FILTER=ESMF_FieldBundleCreate(name='Filt_Bundle Filter' &  !<-- The Bundle's name
                                                                    ,rc  =RC)
!
          domain_int_state%FILT_BUNDLE_RESTORE=ESMF_FieldBundleCreate(name='Filt_Bundle Restore' &  !<-- The Bundle's name
                                                                     ,rc  =RC)
!
          UBOUND_VARS=SIZE(solver_int_state%VARS)
!
          CALL BUILD_FILT_BUNDLE(GRID_DOMAIN                            &
                                ,UBOUND_VARS                            &
                                ,solver_int_state%VARS                  &
                                ,domain_int_state%FILT_BUNDLE_FILTER    &
                                ,domain_int_state%NUM_FIELDS_FILTER_2D  &
                                ,domain_int_state%NUM_FIELDS_FILTER_3D  &
                                ,domain_int_state%FILT_BUNDLE_RESTORE   &
                                ,domain_int_state%NUM_FIELDS_RESTORE_2D &
                                ,domain_int_state%NUM_FIELDS_RESTORE_3D & 
                                ,RESTARTED_RUN)

!
          domain_int_state%FIRST_FILTER=.TRUE.
!
          NULLIFY(domain_int_state%DOLPH_WGTS)
          NULLIFY(domain_int_state%SAVE_2D)
          NULLIFY(domain_int_state%SAVE_3D)
          NULLIFY(domain_int_state%SAVE_2D_PHYS)
          NULLIFY(domain_int_state%SAVE_3D_PHYS)
!
!-----------------------------------------------------------------------
!***  We want to be able to run the digital filter on a group of
!***  parents and children where some of the children are not 
!***  active in the filtering.  When the free forecast begins then
!***  all domains will be active.  Set flags in the Domain export
!***  state indicating if a domain will be inactive if the digital
!***  filter runs as well as if any of its children will not be
!***  active.
!-----------------------------------------------------------------------
!
          I_AM_ACTIVE=.TRUE.
!
          IF(NUM_CHILDREN>0)THEN
!
            DO N=1,NUM_CHILDREN                               !<-- Loop through the children
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Extract Children's DFI method from Config File"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
              CALL ESMF_ConfigGetAttribute(config=CF(MY_CHILDREN_ID(N))  &  !<-- The child's config object
                                          ,value =FILT_METHOD_CHILD      &  !<-- Child's digital filter methodigital filter method
                                          ,label ='filter_method:'       &  !<-- Give this label's value to the previous variable
                                          ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
              IF(FILT_METHOD_CHILD>0)THEN
                CHILD_ACTIVE(N)=.TRUE.
              ELSE
                CHILD_ACTIVE(N)=.FALSE.
              ENDIF
!
            ENDDO
!
          ENDIF
!
!-----------------------------------------------------------------------
!
        ELSEIF(FILT_METHOD==0)THEN
!
          I_AM_ACTIVE=.FALSE.        
!
          IF(NUM_CHILDREN>0)THEN
            DO N=1,NUM_CHILDREN
              CHILD_ACTIVE(N)=.FALSE.                                      !<-- Children can run DFI only if their parent does.
            ENDDO
          ENDIF
!
        ENDIF
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Add DFI flag for this domain into Domain export state"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Domain component export state
                              ,name ='I Am Active'                      &  !<-- Use this name inside the state
                              ,value=I_AM_ACTIVE                        &  !<-- The logical being inserted into the export state
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        IF(NUM_CHILDREN>0)THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Insert child DFI flags into Domain Export State"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeSet(state    =EXP_STATE                    &  !<-- The Domain export state
                                ,name     ='Child Active'               &  !<-- Name of the attribute to insert
                                ,itemCount=NUM_CHILDREN                 &  !<-- Number of items in the array
                                ,valueList=CHILD_ACTIVE                 &  !<-- Insert this attribute.
                                ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          DEALLOCATE(CHILD_ACTIVE)
!
        ENDIF
!
      ENDIF
!
      IF(ASSOCIATED(MY_CHILDREN_ID))THEN
        DEALLOCATE(MY_CHILDREN_ID)
      ENDIF
!
!-----------------------------------------------------------------------
!***  Coupled runs will use counters for accumulating export fields.
!***  Initialize those counters.
!-----------------------------------------------------------------------
!
      domain_int_state%KOUNT_NPRECIP=0
      domain_int_state%KOUNT_NPHS   =0
!
!-----------------------------------------------------------------------
!
      timers(my_domain_id)%total_integ_tim=(timef()-btim0)
!
!-----------------------------------------------------------------------
!
      IF(RC_INIT==ESMF_SUCCESS)THEN
!       WRITE(0,*)'DOMAIN INITIALIZE step succeeded'
      ELSE
        WRITE(0,*)'DOMAIN INITIALIZE step failed RC_INIT=',RC_INIT
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE DOMAIN_INITIALIZE
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE DOMAIN_RUN(DOMAIN_GRID_COMP                            &
                           ,IMP_STATE                                   &
                           ,EXP_STATE                                   &
                           ,CLOCK_DOMAIN                                &
                           ,RC_RUN)
!
!-----------------------------------------------------------------------
!***  The Run step of the Domain component for the NMM.
!***  The forecast tasks execute a single timestep in the Run step
!***  of the NMM-B Solver.  That is the Run subroutine specified
!***  in the Solver Register routine and is called SOLVER_RUN.
!-----------------------------------------------------------------------
!
      USE MODULE_NESTING,ONLY: BOUNDARY_DATA_STATE_TO_STATE
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(ESMF_GridComp) :: DOMAIN_GRID_COMP                              !<-- The Domain gridded component
!
      TYPE(ESMF_State) :: IMP_STATE                                     &  !<-- The Domain component's import state
                         ,EXP_STATE                                        !<-- The Domain component's export state
!
      TYPE(ESMF_Clock) :: CLOCK_DOMAIN                                     !<-- The Domain ESMF Clock
!
      INTEGER,INTENT(OUT) :: RC_RUN                                        !<-- Return code for the Run step
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(ESMF_KIND_I4) :: INTEGER_DT,NEXT_MOVE_TIMESTEP            &
                              ,NUMERATOR_DT,IDENOMINATOR_DT     
!
      INTEGER(kind=ESMF_KIND_I8) :: NTIMESTEP_ESMF                         !<-- The current forecast timestep
!
      INTEGER(kind=KINT) :: I,I_INC,ITE,ITS                             &
                           ,J,J_INC,JTE,JTS
!
      INTEGER(kind=KINT) :: I_SW_PARENT_NEW,J_SW_PARENT_NEW
!
      INTEGER(kind=KINT) :: FILTER_METHOD,HDIFF_ON,IERR,J_CENTER        &
                           ,LAST_STEP_MOVED,RC,NC,YY,MM,DD,H,M,S
!
      INTEGER(kind=KINT),DIMENSION(2) :: STORM_CENTER
!
      INTEGER(kind=KINT),DIMENSION(1:NUM_DOMAINS_MAX) :: NTIMESTEPCHILD_MOVES
!
      REAL(kind=KFPT) :: DLM,DPH,GLATX,GLONX,RAD2DEG,TLATX,TLONX,X,Y,Z
!
      REAL(kind=KFPT),DIMENSION(1:2) :: SW_X
!
      REAL(kind=KDBL),DIMENSION(:,:),POINTER :: GLAT_DBL,GLON_DBL       &
                                               ,VLAT_DBL,VLON_DBL
!
      LOGICAL(kind=KLOG) :: E_BDY,N_BDY,S_BDY,W_BDY                        !<-- Are tasks on a domain boundary?
      LOGICAL(kind=KLOG) :: DIG_FILTER,FREE_FORECAST                    &
                           ,I_AM_ACTIVE,MOVE_NOW                        &
                           ,MOVED_THIS_TIMESTEP
!
      REAL(kind=KFPT) :: DT,NPRECIP_STEP,NPHS_STEP,RECIP_KOUNT,RECIP_NPRECIP
!
      TYPE(ESMF_TimeInterval) :: DT_ESMF
      TYPE(ESMF_Time) :: CURRTIME
!
      TYPE(ESMF_Config) :: CF  
!
      TYPE(WRAP_SOLVER_INT_STATE) :: WRAP_SOLVER
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      btim0=timef()
!
      RC    =ESMF_SUCCESS
      RC_RUN=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  What is this domain's ID?
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Domain_Run: Extract Domain ID from Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state       =IMP_STATE                     &  !<-- The Domain import state
                            ,name        ='DOMAIN_ID'                   &  !<-- Name of the attribute to extract
                            ,value       =MY_DOMAIN_ID                  &  !<-- The ID of this domain
                            ,defaultValue=1                             &  !<-- The default value
                            ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Extract this domain's Virtual Machine so we can distinguish
!***  its MPI task specifications from those of other domains that
!***  the current task may also lie on.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="DOMAIN_RUN: Retrieve VM from Domain component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompGet(gridcomp=DOMAIN_GRID_COMP                   &  !<-- The Domain component
                           ,vm      =VM                                 &  !<-- Get the Virtual Machine from the Domain component
                           ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  What is this task's rank in this domain's set of tasks?
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="DOMAIN_RUN: Obtain Task IDs"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_VMGet(vm             =VM                                &  !<-- The virtual machine
                     ,localpet       =MYPE                              &  !<-- Each MPI task ID
                     ,mpiCommunicator=COMM_MY_DOMAIN                    &  !<-- This domain's intracommunicator
                     ,rc             =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Extract this domain's internal state so we can access
!***  its variables.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="DOMAIN_RUN: Extract the Domain's internal state"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompGetInternalState(DOMAIN_GRID_COMP               &  !<-- The Domain component
                                        ,WRAP                           &  !<-- The F90 wrap of the domain's internal state
                                        ,RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      DOMAIN_INT_STATE=>wrap%DOMAIN_INT_STATE                              !<-- The domain's internal state
!
!-----------------------------------------------------------------------
!***  Extract the timestep from the Clock so that we know the direction
!***  of the integration.  We skip all aspects of physics if the time
!***  step is negative.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
      MESSAGE_CHECK="DOMAIN_Run: Extract the ESMF Timestep"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
      CALL ESMF_ClockGet(clock       =CLOCK_DOMAIN                      &
                        ,timeStep    =DT_ESMF                           &
                        ,currTime    =CURRTIME                          &
                        ,advanceCount=NTIMESTEP_ESMF                    &  !<-- # of times the clock has advanced
                        ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
      CALL ESMF_TimeGet(time=CURRTIME,mm=MM,dd=DD,h=H,m=M,s=S,rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
      MESSAGE_CHECK="DOMAIN_Run: Extract Components of the Timestep" 
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
      CALL ESMF_TimeIntervalGet(timeinterval=DT_ESMF                    &  !<-- the ESMF timestep
                               ,s           =INTEGER_DT                 &  !<-- the integer part of the timestep in seconds
                               ,sN          =NUMERATOR_DT               &  !<-- the numerator of the fractional second
                               ,sD          =IDENOMINATOR_DT            &  !<-- the denominator of the fractional second
                               ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
!-----------------------------------------------------------------------
!
      NTIMESTEP=NTIMESTEP_ESMF
!
      MOVED_THIS_TIMESTEP=.FALSE.
!
!-----------------------------------------------------------------------
!
      fcst_pes: IF(MYPE<domain_int_state%NUM_PES_FCST)THEN                 !<-- Only the forecast tasks integrate
!
!-----------------------------------------------------------------------
!***  We must transfer the horizontal diffusion flag from the
!***  Domain import state to the Solver import state.  The Solver
!***  import state is not available outside of the integration time
!***  loop in NMM_INTEGRATE which lies in NMM_GRID_COMP therefore we
!***  need to perform this transfer each timestep.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract Horizontal Diffusion Flag from Domain Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state       =IMP_STATE                   &  !<-- The Domain import state
                              ,name        ='HDIFF'                     &  !<-- Name of the attribute to extract
                              ,value       =HDIFF_ON                    &  !<-- The ID of this domain
                              ,defaultValue=1                           &  !<-- The default value
                              ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Add Horizontal Diffusion Flag to the SOLVER Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=domain_int_state%IMP_STATE_SOLVER  &  !<-- The Solver component import state
                              ,name ='HDIFF'                            &  !<-- Use this name inside the state
                              ,value=HDIFF_ON                           &  !<-- The scalar being inserted into the import state
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  For digital filtering check to see if we are in the free forecast
!***  and if this domain is active during digital filtering.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract Free Forecast flag from Domain Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state       =IMP_STATE                   &  !<-- The Domain import state
                              ,name        ='Free Forecast'             &  !<-- Name of the attribute to extract
                              ,value       =FREE_FORECAST               &  !<-- Is this the free forecast?
                              ,defaultValue=.true.                      &  !<-- The default value
                              ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        DIG_FILTER=.FALSE.
        IF(.NOT.FREE_FORECAST)THEN
          DIG_FILTER=.TRUE.
        ENDIF       
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract DFI Active flag from Domain Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state=EXP_STATE                          &  !<-- The Domain import state
                              ,name ='I Am Active'                      &  !<-- Name of the attribute to extract
                              ,value=I_AM_ACTIVE                        &  !<-- Does this domain participate in the digital filter?
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Extract the Solver internal state.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract Solver Internal State for Bndry Info"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_GridCompGetInternalState(domain_int_state%SOLVER_GRID_COMP &  !<-- The Solver component
                                          ,WRAP_SOLVER                       &
                                          ,RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        SOLVER_INT_STATE=>wrap_solver%INT_STATE
!
!-----------------------------------------------------------------------
!***  Determine if forecast tasks are on the domain boundary.
!-----------------------------------------------------------------------
!
        S_BDY=(solver_int_state%JTS==solver_int_state%JDS)                 ! This task is on the southern boundary
        N_BDY=(solver_int_state%JTE==solver_int_state%JDE)                 ! This task is on the northern boundary
        W_BDY=(solver_int_state%ITS==solver_int_state%IDS)                 ! This task is on the western boundary
        E_BDY=(solver_int_state%ITE==solver_int_state%IDE)                 ! This task is on the eastern boundary       
!
!-----------------------------------------------------------------------
!***  If this is a nested run then we need to consider two things:
!***   (1) For all nests new boundary data must be moved from the
!***       Domain import state to the Solver import state every
!***       N timesteps where N is the number of the nest's timesteps
!***       within one timestep of its parent.  This is done before
!***       the Run step of the Solver is executed in order that
!***       the nests have correct boundary conditions for integrating
!***       through the next N timesteps.
!***   (2) If this is a moving nest and it has just moved then it
!***       updates those of its interior points that still lie
!***       within the footprint of its domain's pre-move location.
!***       These are the points that are NOT updated by its parent.
!***       Then it must also incorporate any interior update data
!***       that was sent to it by its parent that lie outside of
!***       the nest domain's pre-move footprint.
!-----------------------------------------------------------------------
!
        nests: IF(domain_int_state%I_AM_A_NEST)THEN
!
          IF(FREE_FORECAST.OR.(DIG_FILTER.AND.I_AM_ACTIVE))THEN
            CALL BOUNDARY_DATA_STATE_TO_STATE(s_bdy    =S_BDY                            &  !<-- This task lies on a south boundary?
                                             ,n_bdy    =N_BDY                            &  !<-- This task lies on a north boundary?
                                             ,w_bdy    =W_BDY                            &  !<-- This task lies on a west boundary?
                                             ,e_bdy    =E_BDY                            &  !<-- This task lies on an east boundary?
                                             ,clock    =CLOCK_DOMAIN                     &  !<-- The Domain Clock
                                             ,nest     =domain_int_state%I_AM_A_NEST     &  !<-- The nest flag (yes or no)
                                             ,ratio    =PARENT_CHILD_TIME_RATIO          &  !<-- # of child timesteps per parent timestep
                                             ,state_in =IMP_STATE                        &  !<-- Domain component's import state
                                             ,state_out=domain_int_state%IMP_STATE_SOLVER)  !<-- The Solver import state
          ENDIF
!
!-----------------------------------------------------------------------
!
          domain_moves: IF(domain_int_state%MY_DOMAIN_MOVES)THEN           !<-- Select the moving nests
!
!-----------------------------------------------------------------------
!***  If a nest is moving in its next timestep then update the
!***  SW corner location in the Solver internal state now at the
!***  end of the preceding timestep.  This is necessary for the
!***  correct location to be in place if a restart file is to be
!***  written that is valid for the beginning of the next timestep.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="DOMAIN_RUN: Get NEXT_MOVE_TIMESTEP from Import State"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_AttributeGet(state=IMP_STATE                      &  !<-- Domain component's import state
                                  ,name ='NEXT_MOVE_TIMESTEP'           &  !<-- Extract Attribute with this name
                                  ,value=NEXT_MOVE_TIMESTEP             &  !<-- When does this nest move again?
                                  ,rc   =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            solver_int_state%NMTS=NEXT_MOVE_TIMESTEP                       !<-- Save in the Solver internal state
!
            IF(NTIMESTEP==NEXT_MOVE_TIMESTEP-1)THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
              MESSAGE_CHECK="Domain_Run: Get I_SHIFT,J_SHIFT from Domain Import State"
!             CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
              CALL ESMF_AttributeGet(state=IMP_STATE                    &  !<-- The Domain import state
                                    ,name ='I_SHIFT'                    &  !<-- Get Attribute with this name
                                    ,value=I_SHIFT_CHILD                &  !<-- Motion of the nest in I on its grid
                                    ,rc   =RC )
              CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
!
              CALL ESMF_AttributeGet(state=IMP_STATE                    &  !<-- The Domain import state
                                    ,name ='J_SHIFT'                    &  !<-- Get Attribute with this name
                                    ,value=J_SHIFT_CHILD                &  !<-- Motion of the nest in J on its grid
                                    ,rc   =RC )
              CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
!
              CALL ESMF_AttributeGet(state=IMP_STATE                    &  !<-- The Domain import state
                                    ,name ='LAST_STEP_MOVED'            &  !<-- Get Attribute with this name
                                    ,value=LAST_STEP_MOVED              &  !<-- Motion of the nest in J on its grid
                                    ,rc   =RC )
              CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
              solver_int_state%LAST_STEP_MOVED=LAST_STEP_MOVED
!
            ENDIF
!
!-----------------------------------------------------------------------
!***  We need to update the Solver internal state's values for the
!***  'new' location of the nest's SW corner every timestep to 
!***  also include those cases where the parent shifts when the
!***  nest does not or else the restart file will be incorrect.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Domain_Run: Get I_SW_PARENT_NEW,J_SW_PARENT_NEW from Domain Import State"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_AttributeGet(state=IMP_STATE                      &  !<-- The Domain import state
                                  ,name ='I_SW_PARENT_NEW'              &  !<-- Get Attribute with this name
                                  ,value=I_SW_PARENT_NEW                &  !<-- Motion of the nest in I on its grid
                                  ,rc   =RC   )
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
!
            CALL ESMF_AttributeGet(state=IMP_STATE                      &  !<-- The Domain import state
                                  ,name ='J_SW_PARENT_NEW'              &  !<-- Get Attribute with this name
                                  ,value=J_SW_PARENT_NEW                &  !<-- Motion of the nest in J on its grid
                                  ,rc   =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            IF(I_SW_PARENT_NEW>-999999)THEN
!
              solver_int_state%I_PAR_STA=I_SW_PARENT_NEW
!
              solver_int_state%J_PAR_STA=J_SW_PARENT_NEW
!
            ENDIF
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="DOMAIN_RUN: Get the MOVE_NOW Flag from Import State"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_AttributeGet(state=IMP_STATE                      &  !<-- Domain component's import state
                                  ,name ='MOVE_NOW'                     &  !<-- Extract Attribute with this name
                                  ,value=MOVE_NOW                       &  !<-- Is the child moving right now?
                                  ,rc   =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="DOMAIN_RUN: Add MOVE_NOW Flag to the Solver Import State"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_AttributeSet(state=domain_int_state%IMP_STATE_SOLVER &  !<-- The Solver component import state
                                  ,name ='MOVE_NOW'                        &  !<-- Use this name inside the state
                                  ,value=MOVE_NOW                          &  !<-- Did this nest move this timestep?
                                  ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
            moving_now: IF(MOVE_NOW)THEN                                   !<-- Select moving nests that move this timestep
!
!-----------------------------------------------------------------------
!***  What are the nest's new transformed lat/lon of its south and
!***  west boundaries following the move?
!-----------------------------------------------------------------------
!
              MOVED_THIS_TIMESTEP=.TRUE.
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
              MESSAGE_CHECK="Domain_Run: Get I_SHIFT,J_SHIFT from Domain Import State"
!             CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
              CALL ESMF_AttributeGet(state=IMP_STATE                    &  !<-- The Domain import state
                                    ,name ='I_SHIFT'                    &  !<-- Get Attribute with this name
                                    ,value=I_SHIFT_CHILD                &  !<-- Motion of the nest in I on its grid
                                    ,rc   =RC )
              CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
!
              CALL ESMF_AttributeGet(state=IMP_STATE                    &  !<-- The Domain import state
                                    ,name ='J_SHIFT'                    &  !<-- Get Attribute with this name
                                    ,value=J_SHIFT_CHILD                &  !<-- Motion of the nest in J on its grid
                                    ,rc   =RC )
              CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  What are the transformed lat/lon of the nest's SW corner at
!***  the new nest location?
!-----------------------------------------------------------------------
!
              IF(MYPE==0)THEN
!
                D_ONE=1.
                D_180=180.
                PI_LOC=DACOS(-D_ONE)
                D2R=PI_LOC/D_180
!
                TPH0_1=TPH0D_1*D2R                                         !<-- The central lat/lon of domain #1 is the center
                TLM0_1=TLM0D_1*D2R                                         !    for all grid-associated nests
!
                ITS=solver_int_state%ITS
                JTS=solver_int_state%JTS
                GLATX=solver_int_state%GLAT(ITS,JTS)                       !<-- Geographic lat (radians) of nest's pre-move SW corner
                GLONX=solver_int_state%GLON(ITS,JTS)                       !<-- Geographic lon (radians) of nest's pre-move SW corner
!
                X=COS(TPH0_1)*COS(GLATX)*COS(GLONX-TLM0_1)+SIN(TPH0_1)*SIN(GLATX)
                Y=COS(GLATX)*SIN(GLONX-TLM0_1)
                Z=-SIN(TPH0_1)*COS(GLATX)*COS(GLONX-TLM0_1)+COS(TPH0_1)*SIN(GLATX)
!
                TLATX=ATAN(Z/SQRT(X*X+Y*Y))                                !<-- Transformed lat (radians) of nest's pre-move SW corner
                TLONX=ATAN(Y/X)                                            !<-- Transformed lon (radians) of nest's pre-move SW corner
                IF(X<0)TLONX=TLONX+PI
!
                SB_1=SBD_1*D2R                                             !<-- Transformed lat (radians) of domain #1's S bndry
                WB_1=WBD_1*D2R                                             !<-- Transformed lon (radians) of domain #1's W bndry
!
                DPH=solver_int_state%DPHD*D2R                              !<-- Nest's angular grid increment in J (radians)
                DLM=solver_int_state%DLMD*D2R                              !<-- Nest's angular grid increment in I (radians)
!
                TLATX=TLATX+J_SHIFT_CHILD*DPH                              !<-- Transformed lat (radians) of nest's post-move SW corner
                TLONX=TLONX+I_SHIFT_CHILD*DLM                              !<-- Transformed lon (radians) of nest's post-move SW corner
!
                I_INC=NINT((TLONX-WB_1)/DLM)                               !<-- Nest grid increments (integer) between west/south
                J_INC=NINT((TLATX-SB_1)/DPH)                               !    boundaries of the nest and domain #1.
!
                SW_X(1)=(SB_1+J_INC*DPH)/D2R                               !<-- Transformed lat (degrees) of nest domain's S bndry
                SW_X(2)=(WB_1+I_INC*DLM)/D2R                               !<-- Transformed lon (degrees) of nest domain's S bndry
!
              ENDIF
!
!-----------------------------------------------------------------------
!***  Local task 0 shares the transformed lat/lon of the nest domain's
!***  south and west boundaries with all other fcst tasks.
!-----------------------------------------------------------------------
!
              CALL MPI_BCAST(SW_X                                       &
                            ,2                                          &
                            ,MPI_REAL                                   &
                            ,0                                          &
                            ,COMM_FCST_TASKS(MY_DOMAIN_ID)              &
                            ,IERR )
!
              solver_int_state%SBD=SW_X(1)
              solver_int_state%WBD=SW_X(2)
!
!-----------------------------------------------------------------------
!***  Update quantities that are directly related to the nest's grid.
!-----------------------------------------------------------------------
!
              CALL GRID_CONSTS(solver_int_state%GLOBAL                           &
                              ,solver_int_state%DT                               &
                              ,solver_int_state%SMAG2                            &
                              ,solver_int_state%CODAMP,solver_int_state%WCOR     &
                              ,solver_int_state%TPH0D,solver_int_state%TLM0D     &
                              ,solver_int_state%SBD,solver_int_state%WBD         &
                              ,solver_int_state%DPHD,solver_int_state%DLMD       &
                              ,solver_int_state%DXH,solver_int_state%RDXH        &
                              ,solver_int_state%DXV,solver_int_state%RDXV        &
                              ,solver_int_state%DYH,solver_int_state%RDYH        &
                              ,solver_int_state%DYV,solver_int_state%RDYV        &
                              ,solver_int_state%DDV,solver_int_state%RDDV        &
                              ,solver_int_state%DDMPU,solver_int_state%DDMPV     &
                              ,solver_int_state%WPDAR                            &
                              ,solver_int_state%FCP,solver_int_state%FDIV        &
                              ,solver_int_state%CURV,solver_int_state%F          &
                              ,solver_int_state%FAD,solver_int_state%FAH         &
                              ,solver_int_state%DARE,solver_int_state%RARE       &
                              ,solver_int_state%GLAT,solver_int_state%GLON       &
                              ,solver_int_state%GLAT_SW,solver_int_state%GLON_SW &
                              ,solver_int_state%VLAT,solver_int_state%VLON       &
                              ,solver_int_state%HDACX,solver_int_state%HDACY     &
                              ,solver_int_state%HDACVX,solver_int_state%HDACVY   &
                              ,solver_int_state%E_BDY,solver_int_state%N_BDY     &
                              ,solver_int_state%S_BDY,solver_int_state%W_BDY     &
                              ,solver_int_state%ITS,solver_int_state%ITE         &
                              ,solver_int_state%JTS,solver_int_state%JTE         &
                              ,solver_int_state%IMS,solver_int_state%IME         &
                              ,solver_int_state%JMS,solver_int_state%JME         &
                              ,solver_int_state%IDS,solver_int_state%IDE         &
                              ,solver_int_state%JDS,solver_int_state%JDE )
!
!-----------------------------------------------------------------------
!***  Update all nest points that remain inside of its domain's 
!***  pre-move footprint.
!-----------------------------------------------------------------------
!
      btim=timef()
!
              CALL UPDATE_INTERIOR_FROM_NEST(IMP_STATE                       &  !<-- Domain import state (for nest's I_SHIFT and J_SHIFT)
                                            ,domain_int_state%MOVE_BUNDLE_H  &  !<-- The Bundle of pointers to update H variables
                                            ,NUM_FIELDS_MOVE_2D_H_I          &  !<-- Total # of 2-D integer H Fields in the Bundle
                                            ,NUM_FIELDS_MOVE_2D_H_R          &  !<-- Total # of 2-D real H Fields in the Bundle
                                            ,NUM_FIELDS_MOVE_3D_H            &  !<-- Total # of 3-D H Fields in the Bundle
                                            ,NUM_LEVELS_MOVE_3D_H            &  !<-- Total # of 2-D levels in all 3-D H update arrays
                                            ,domain_int_state%MOVE_BUNDLE_V  &  !<-- The Bundle of pointers to update V variables
                                            ,NUM_FIELDS_MOVE_2D_V            &  !<-- Total # of 2-D V Fields in the Bundle
                                            ,NUM_FIELDS_MOVE_3D_V            &  !<-- Total # of 3-D V Fields in the Bundle
                                            ,NUM_LEVELS_MOVE_3D_V            &  !<-- Total # of 2-D levels in all 3-D V update arrays
                                            ,solver_int_state%INPES          &  !<-- # of tasks in east-west on this domain
                                            ,solver_int_state%JNPES          &  !<-- # of tasks in north-south on this domain
                                            ,solver_int_state%ITS            &  !<-- Starting integration index in I
                                            ,solver_int_state%ITE            &  !<-- Ending integration index in I
                                            ,solver_int_state%JTS            &  !<-- Starting integration index in J
                                            ,solver_int_state%JTE            &  !<-- Ending integration index in J
                                            ,solver_int_state%IMS            &  !<-- Starting memory index in I
                                            ,solver_int_state%IME            &  !<-- Ending memory index in I
                                            ,solver_int_state%JMS            &  !<-- Starting memory index in J
                                            ,solver_int_state%JME            &  !<-- Ending memory index in J
                                            ,solver_int_state%IDS            &  !<-- Starting domain index in I
                                            ,solver_int_state%IDE            &  !<-- Ending domain index in I
                                            ,solver_int_state%JDS            &  !<-- Starting domain index in J
                                            ,solver_int_state%JDE            &  !<-- Ending domain index in J
                                              )
!                      
      timers(my_domain_id)%update_interior_from_nest_tim=               &
           timers(my_domain_id)%update_interior_from_nest_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Update all nest points that have moved outside of its domain's
!***  pre-move footprint.
!-----------------------------------------------------------------------
!
      btim=timef()
!
              CALL UPDATE_INTERIOR_FROM_PARENT(IMP_STATE                       &  !<-- The Domain import state
                                              ,domain_int_state%SFC_FILE_RATIO &  !<-- Ratio of upper parent grid increment to this domain's
                                              ,domain_int_state%MOVE_BUNDLE_H  &  !<-- The Bundle of pointers to update H variables
                                              ,NUM_FIELDS_MOVE_2D_H_I          &  !<-- Total # of 2-D integer H Fields in the Bundle
                                              ,NUM_FIELDS_MOVE_2D_H_R          &  !<-- Total # of 2-D real H Fields in the Bundle
                                              ,NUM_FIELDS_MOVE_3D_H            &  !<-- Total # of 3-D H Fields in the Bundle
                                              ,domain_int_state%MOVE_BUNDLE_V  &  !<-- The Bundle of pointers to update V variables
                                              ,NUM_FIELDS_MOVE_2D_V            &  !<-- Total # of 2-D V Fields in the Bundle
                                              ,NUM_FIELDS_MOVE_3D_V            &  !<-- Total # of 3-D V Fields in the Bundle
                                              ,solver_int_state%GLAT           &  !<-- This domain's geographic latitude (radians)
                                              ,solver_int_state%GLON           &  !<-- This domain's geographic longitude (radians)
                                              ,solver_int_state%ITS            &  !<-- Starting integration index in I
                                              ,solver_int_state%ITE            &  !<-- Ending integration index in I
                                              ,solver_int_state%JTS            &  !<-- Starting integration index in J
                                              ,solver_int_state%JTE            &  !<-- Ending integration index in J
                                              ,solver_int_state%IMS            &  !<-- Starting memory index in I
                                              ,solver_int_state%IME            &  !<-- Ending memory index in I
                                              ,solver_int_state%JMS            &  !<-- Starting memory index in J
                                              ,solver_int_state%JME )             !<-- Ending memory index in J
!
      timers(my_domain_id)%update_interior_from_parent_tim=             &
           timers(my_domain_id)%update_interior_from_parent_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="DOMAIN_RUN: Reset the MOVE_NOW Flag to False"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
              MOVE_NOW=.FALSE.
!
              CALL ESMF_AttributeSet(state=IMP_STATE                    &  !<-- Domain component's import state
                                    ,name ='MOVE_NOW'                   &  !<-- Set Attribute with this name
                                    ,value=MOVE_NOW                     &  !<-- Value is reset to false
                                    ,rc   =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
            ENDIF moving_now
!
!-----------------------------------------------------------------------
!
          ENDIF domain_moves
!
!-----------------------------------------------------------------------
!
        ENDIF nests
!
!-----------------------------------------------------------------------
!***  Parents update the next timesteps their children move for the
!***  Solver internal state in case it is needed for restart output.
!-----------------------------------------------------------------------
!
        IF(domain_int_state%I_AM_A_PARENT)THEN      
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="DOMAIN_RUN: Get NTIMESTEP_CHILD_MOVES from Import State"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeGet(state    =IMP_STATE                      &  !<-- Domain component's import state
                                ,name     ='NEXT_TIMESTEP_CHILD_MOVES'    &  !<-- Extract Attribute with this name
                                ,valueList=NTIMESTEP_CHILD_MOVES          &  !<-- When do the children move again?
                                ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          solver_int_state%NTSCM=NTIMESTEP_CHILD_MOVES                       !<-- Save in the Solver internal state
!
        ENDIF
!
!-----------------------------------------------------------------------
!***  Set the filter method in the Solver internal state with the
!***  value put into the Domain import state during NMM_RUN.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract Filter method from import state"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state       =IMP_STATE                   &  !<-- The Domain import state
                              ,name        ='Filter_Method'             &  !<-- Name of the attribute to extract
                              ,value       =FILTER_METHOD               &  !<-- The scalar being extracted from the import state
                              ,defaultValue=0                           &  !<-- The default value
                              ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        solver_int_state%FILTER_METHOD=FILTER_METHOD
!
!-----------------------------------------------------------------------
!***  If this is a coupling timestep then update the fields in the
!***  Solver internal state with the imported fields in the Domain
!***  internal state that are connected to the NMM-B cap.
!-----------------------------------------------------------------------
!
        DT=solver_int_state%DT
!
        IF(ATM_OCN_CPL_INT>0)THEN
          IF(MOD(NTIMESTEP*DT,ATM_OCN_CPL_INT)<DT)THEN                     !<-- Is this a coupling timestep?
!
            DO J=solver_int_state%JTS,solver_int_state%JTE
            DO I=solver_int_state%ITS,solver_int_state%ITE
              IF(solver_int_state%SM(I,J)>0.5)THEN
                IF(domain_int_state%SST_COUPLED(I,J)>265.               &
                                 .AND.                                  &
                   domain_int_state%SST_COUPLED(I,J)<325.)THEN
                  solver_int_state%SST(I,J)=domain_int_state%SST_COUPLED(I,J)  !<-- Insert imported SST (K) into the Solver internal state
                ENDIF
              ENDIf
            ENDDO
            ENDDO
!
          ENDIF
        ENDIF
!
!-----------------------------------------------------------------------
!***  We are now ready to execute the Run step of the Solver component.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
        MESSAGE_CHECK="Execute the Run Step for the Solver"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
!       call print_memory()
!
        CALL ESMF_GridCompRun(gridcomp   =domain_int_state%SOLVER_GRID_COMP  &  !<-- The Solver component
                             ,importState=domain_int_state%IMP_STATE_SOLVER  &  !<-- The Solver import state
                             ,exportState=domain_int_state%EXP_STATE_SOLVER  &  !<-- The Solver export state
                             ,clock      =CLOCK_DOMAIN                       &  !<-- The Domain Clock
                             ,rc         =RC)        
!
!       call print_memory()
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  If this is a coupled run then prepare the export fields.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
        export: IF(ATM_OCN_CPL_INT>0)THEN                                  !<-- If true then this is a coupled run.
!
!-----------------------------------------------------------------------
!***  Accumulate those fields to be exported from the atmosphere for
!***  which mean values through the coupling interval are needed.
!***  If a new accumulation period is starting then first zero out
!***  the accumulation arrays.
!-----------------------------------------------------------------------
!
!-----------------------------------------
!***  Accumulate at the physics interval.
!-----------------------------------------
!
          IF(MOD(NTIMESTEP,solver_int_state%NPHS)==0)THEN
            IF(domain_int_state%KOUNT_NPHS==0)THEN
              DO J=solver_int_state%JTS,solver_int_state%JTE
              DO I=solver_int_state%ITS,solver_int_state%ITE
                domain_int_state%MEAN_ZONAL_MOM_FLX_COUPLED(I,J)=0.
                domain_int_state%MEAN_MERID_MOM_FLX_COUPLED(I,J)=0.
              ENDDO
              ENDDO
            ENDIF
!
            domain_int_state%KOUNT_NPHS=domain_int_state%KOUNT_NPHS+1
!
            DO J=solver_int_state%JTS,solver_int_state%JTE
            DO I=solver_int_state%ITS,solver_int_state%ITE
!
              domain_int_state%MEAN_ZONAL_MOM_FLX_COUPLED(I,J)=                            &   !<-- Zonal momentum flx (Pa)
                domain_int_state%MEAN_ZONAL_MOM_FLX_COUPLED(I,J)+solver_int_state%TAUX(I,J)
!
              domain_int_state%MEAN_MERID_MOM_FLX_COUPLED(I,J)=                            &   !<-- Meridional momentum flx (Pa)
                domain_int_state%MEAN_MERID_MOM_FLX_COUPLED(I,J)+solver_int_state%TAUY(I,J)
!
            ENDDO
            ENDDO
!
          ENDIF
!
!-----------------------------------------------
!***  Accumulate at the precipitation interval.
!-----------------------------------------------
!
          IF(MOD(NTIMESTEP,solver_int_state%NPRECIP)==0)THEN
!
            IF(domain_int_state%KOUNT_NPRECIP==0)THEN
              DO J=solver_int_state%JTS,solver_int_state%JTE
              DO I=solver_int_state%ITS,solver_int_state%ITE
                domain_int_state%MEAN_PREC_RATE_COUPLED(I,J)=0.
              ENDDO
              ENDDO
            ENDIF
!
            domain_int_state%KOUNT_NPRECIP=domain_int_state%KOUNT_NPRECIP+1
            RECIP_NPRECIP=1./REAL(solver_int_state%NPRECIP*solver_int_state%DT)
!
            DO J=solver_int_state%JTS,solver_int_state%JTE
            DO I=solver_int_state%ITS,solver_int_state%ITE
              domain_int_state%MEAN_PREC_RATE_COUPLED(I,J)=                              &   !<-- Precipitation rate (kg m-2 s-1)
                domain_int_state%MEAN_PREC_RATE_COUPLED(I,J)+solver_int_state%PREC(I,J)  &
                                                            *1.E3*RECIP_NPRECIP
            ENDDO
            ENDDO
          ENDIF
!
!-----------------------------------------------------------------------
!***  Now prepare the instantaneous fields to be exported if this is
!***  at the end of the coupling interval.  NTIMESTEP+1 is used for
!***  the arithmetic since we at the END of the current timestep,
!***  i.e., at the end of timestep 0 we have traversed 1 timestep.
!-----------------------------------------------------------------------
!
          IF(MOD((NTIMESTEP+1)*solver_int_state%DT-EPS,ATM_OCN_CPL_INT)<DT)THEN  !<-- At or just before the end of a coupling interval.
!
            DO J=solver_int_state%JTS,solver_int_state%JTE
            DO I=solver_int_state%ITS,solver_int_state%ITE
!
              domain_int_state%INST_SFC_PRESSURE_COUPLED(I,J)=solver_int_state%PINT(I,J,LM) !<-- Export inst sfc pressure (Pa)
!
              domain_int_state%INST_NET_LW_FLX_COUPLED(I,J)=solver_int_state%RLWIN(I,J) &  !<-- Export inst net LW flux (W m-2)
                                                           -solver_int_state%RADOT(I,J)
!
              domain_int_state%INST_NET_SW_FLX_COUPLED(I,J)=solver_int_state%RSWIN(I,J) &  !<-- Export inst net SW flux (W m-2)
                                                           -solver_int_state%RSWOUT(I,J)
!
              domain_int_state%INST_SENS_HT_FLX_COUPLED(I,J)=solver_int_state%TWBS(I,J)   !<-- Export inst sensible heat flux (W m-2)
!
              domain_int_state%INST_LAT_HT_FLX_COUPLED(I,J)=solver_int_state%QWBS(I,J)    !<-- Export inst sensible heat flux (W m-2)
!
            ENDDO
            ENDDO
!
          ENDIF
!
!-----------------------------------------------------------------------
!***  If this is at or immediately before the end of the coupling
!***  interval then prepare the mean export fields accumulated at
!***  the physics interval.
!-----------------------------------------------------------------------
!
          NPHS_STEP=solver_int_state%NPHS*solver_int_state%DT
!
          IF(MOD(NTIMESTEP,solver_int_state%NPHS)==0                    &                              !    At or just before
                           .AND.                                        &                              !<-- a coupling
             ATM_OCN_CPL_INT-MOD((NTIMESTEP+1)*solver_int_state%DT-EPS,ATM_OCN_CPL_INT)<NPHS_STEP)THEN !    interval.
!
            RECIP_KOUNT=1./REAL(domain_int_state%KOUNT_NPHS)
!
            DO J=solver_int_state%JTS,solver_int_state%JTE
            DO I=solver_int_state%ITS,solver_int_state%ITE
!
              domain_int_state%MEAN_ZONAL_MOM_FLX_COUPLED(I,J)=               &   !<-- Export mean zonal momentum flx (N m-2)
                  RECIP_KOUNT*domain_int_state%MEAN_ZONAL_MOM_FLX_COUPLED(I,J)
!
              domain_int_state%MEAN_MERID_MOM_FLX_COUPLED(I,J)=               &   !<-- Export mean meridional momentum flx (N m-2)
                  RECIP_KOUNT*domain_int_state%MEAN_MERID_MOM_FLX_COUPLED(I,J)
!
            ENDDO
            ENDDO
!
            domain_int_state%KOUNT_NPHS=0
!
          ENDIF
!
!-----------------------------------------------------------------------
!***  If this is at the end of the coupling interval or immediately
!***  before it then prepare the mean export fields accumulated at
!***  the precipitation interval.
!-----------------------------------------------------------------------
!
          NPRECIP_STEP=solver_int_state%NPRECIP*solver_int_state%DT
!
          IF(MOD(NTIMESTEP,solver_int_state%NPRECIP)==0                 &                                 !    At or just before
                           .AND.                                        &                                 !<-- a coupling
             ATM_OCN_CPL_INT-MOD((NTIMESTEP+1)*solver_int_state%DT-EPS,ATM_OCN_CPL_INT)<NPRECIP_STEP)THEN !    interval.
!
            RECIP_KOUNT=1./REAL(domain_int_state%KOUNT_NPRECIP)
!
            DO J=solver_int_state%JTS,solver_int_state%JTE
            DO I=solver_int_state%ITS,solver_int_state%ITE
!
              domain_int_state%MEAN_PREC_RATE_COUPLED(I,J)=               &   !<-- Export mean precipitation rate (kg s m-2)
                  RECIP_KOUNT*domain_int_state%MEAN_PREC_RATE_COUPLED(I,J)
            ENDDO
            ENDDO
!
            domain_int_state%KOUNT_NPRECIP=0
!
          ENDIF
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
        ENDIF export
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  For moving nests the hurricane storm tracker determines the
!***  center of the storm on this domain's grid every NPHS*NTRACK 
!***  timesteps to see if the domain should shift.  If this is an
!***  appropriate timestep then load the storm center location into
!***  the Domain component's export state so it can be transferred
!***  to the Parent-Child coupler in NMM_INTEGRATE.
!-----------------------------------------------------------------------
!
        IF(domain_int_state%MY_DOMAIN_MOVES)THEN                           !<-- Select the moving nests
!
          IF(solver_int_state%NTRACK>0                                  &
                  .AND.                                                 &
             (NTIMESTEP==0                                              &
                  .OR.                                                  &
              MOD(NTIMESTEP+1,solver_int_state%NTRACK*solver_int_state%NPHS)==0))THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="DOMAIN_RUN: Set the storm center location."
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            STORM_CENTER(1)=solver_int_state%TRACKER_IFIX
            STORM_CENTER(2)=solver_int_state%TRACKER_JFIX
!
            CALL ESMF_AttributeSet(state    =EXP_STATE                  &  !<-- Domain component's export state
                                  ,name     ='Storm Center'             &  !<-- Set Attribute with this name
                                  ,itemCount=2                          &  !<-- Two words in the Attribute
                                  ,valueList=STORM_CENTER               &  !<-- Load these Attribute values
                                  ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          ENDIF
!
        ENDIF
!
!-----------------------------------------------------------------------
!***  Update ESMF Grid, Fields, and regrid interpolation weights
!***  following a nest's shift for coupled runs.  All the domain's
!***  compute tasks and write tasks need to participate.
!-----------------------------------------------------------------------
!
        IF(ATM_OCN_CPL_INT>0)THEN                                          !<-- If true then this is a coupled run.
!
!-----------------------------------------------------------------------
!***  Write tasks do not know anything about domain motion yet they 
!***  must participate in these ESMF updates when the nest shifts.
!***  The lead compute task therefore must inform the write tasks
!***  (via a broadcast) whether or not motion just occurred.
!-----------------------------------------------------------------------
!
          CALL MPI_BCAST(MOVED_THIS_TIMESTEP                            &  !<-- Did the nest just shift?
                        ,1                                              &  !<-- The signal is 1 word
                        ,MPI_LOGICAL                                    &  !<-- The signal is logical
                        ,0                                              &  !<-- Broadcast from the lead compute task
                        ,COMM_MY_DOMAIN                                 &  !<-- Intracommunicator for all tasks on this domain.
                        ,IERR )
!
          IF(MOVED_THIS_TIMESTEP)THEN
!
!-----------------------------------------------------------------------
!***  Update the moving nest's wind rotation angle, cell areas,
!***  and sea mask in the DOMAIN_DESCRIPTORS object.  Transfer
!***  post-shift values to the ESMF Grid.
!-----------------------------------------------------------------------
!
            ITS=solver_int_state%ITS
            ITE=solver_int_state%ITE
            JTS=solver_int_state%JTS
            JTE=solver_int_state%JTE
!
            CALL ROTANGLE_CELLAREA_SEAMASK(SOLVER_INT_STATE                           &
                                          ,ITS,ITE,JTS,JTE                            &
                                          ,DOMAIN_DESCRIPTORS(MY_DOMAIN_ID)%ROT_ANGLE &
                                          ,DOMAIN_DESCRIPTORS(MY_DOMAIN_ID)%CELL_AREA &
                                          ,DOMAIN_DESCRIPTORS(MY_DOMAIN_ID)%SEA_MASK )
!
            ALLOCATE(GLAT_DBL(ITS:ITE,JTS:JTE))
            ALLOCATE(GLON_DBL(ITS:ITE,JTS:JTE))
            ALLOCATE(VLAT_DBL(ITS:ITE,JTS:JTE))
            ALLOCATE(VLON_DBL(ITS:ITE,JTS:JTE))
!
            RAD2DEG = 180._kdbl/ACOS(-1._kdbl)
            DO J=JTS,JTE
            DO I=ITS,ITE
              GLAT_DBL(I,J)=solver_int_state%GLAT(I,J)*RAD2DEG
              GLON_DBL(I,J)=solver_int_state%GLON(I,J)*RAD2DEG
              VLAT_DBL(I,J)=solver_int_state%VLAT(I,J)*RAD2DEG
              VLON_DBL(I,J)=solver_int_state%VLON(I,J)*RAD2DEG
            ENDDO
            ENDDO
!
            CALL NMMB_GridUpdate(MY_DOMAIN_ID                               &
                                ,DOMAIN_DESCRIPTORS                         &
                                ,DOMAIN_DESCRIPTORS(MY_DOMAIN_ID)%CELL_AREA &
                                ,DOMAIN_DESCRIPTORS(MY_DOMAIN_ID)%SEA_MASK  &
                                ,GLON_DBL                                   &
                                ,GLAT_DBL                                   &
                                ,VLON_DBL                                   &
                                ,VLAT_DBL                                   &
                                ,RC )
!
            DEALLOCATE(GLAT_DBL)
            DEALLOCATE(GLON_DBL)
            DEALLOCATE(VLAT_DBL)
            DEALLOCATE(VLON_DBL)
!
            MOVED_THIS_TIMESTEP=.FALSE.
!
          ENDIF
!
!-----------------------------------------------------------------------
!
        ENDIF
!
!-----------------------------------------------------------------------
!
      ENDIF fcst_pes
!
!-----------------------------------------------------------------------
!
      timers(my_domain_id)%total_integ_tim=timers(my_domain_id)%total_integ_tim+(timef()-btim0)
!
!-----------------------------------------------------------------------
!***  The final error signal information.
!-----------------------------------------------------------------------
!
      IF(RC_RUN==ESMF_SUCCESS)THEN
!       WRITE(0,*)'DOMAIN RUN step succeeded'
      ELSE
        WRITE(0,*)'DOMAIN RUN step failed RC_RUN=',RC_RUN
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE DOMAIN_RUN
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE DOMAIN_FINALIZE(DOMAIN_GRID_COMP                       &
                                ,IMP_STATE                              &
                                ,EXP_STATE                              &
                                ,CLOCK_DOMAIN                           &
                                ,RC_FINALIZE)
!
!-----------------------------------------------------------------------
!***  This routine Finalizes the Domain gridded component.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(ESMF_GridComp) :: DOMAIN_GRID_COMP                              !<-- The Domain gridded component
!
      TYPE(ESMF_State) :: IMP_STATE                                     &  !<-- The Domain finalize step's import state
                         ,EXP_STATE                                        !<-- The Domain finalize step's export state
!
      TYPE(ESMF_Clock) :: CLOCK_DOMAIN                                     !<-- The Domain ESMF Clock
!
      INTEGER,INTENT(OUT) :: RC_FINALIZE                                   !<-- Return code for the Finalize step
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER :: I,J
      INTEGER :: RC                                                        ! The final error signal variables.
!
      CHARACTER(50):: MODE
!
      LOGICAL(kind=KLOG) :: PHYSICS_ON
!
      TYPE(ESMF_Config) :: CF                                              !<-- The config object
!
      TYPE(WRAP_DOMAIN_INTERNAL_STATE) :: WRAP                             !<-- The F90 wrap of the Domain internal state
!
      TYPE(DOMAIN_INTERNAL_STATE),POINTER :: DOMAIN_INT_STATE              !<-- The Domain internal state pointer
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC         =ESMF_SUCCESS
      RC_FINALIZE=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  Retrieve the config object CF from the Domain component.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Retrieve Config Object from Domain Component"
!      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompGet(gridcomp=DOMAIN_GRID_COMP                   &  !<-- The Domain component
                           ,config  =CF                                 &  !<-- The config object (~namelist)
                           ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINALIZE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Extract the diabatic/adiabatic flag from the configure file.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Retrieve Adiabatic Flag from Config Object"
!      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &
                                  ,value =MODE                          &
                                  ,label ='adiabatic:'                  &
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINALIZE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
      IF(TRIM(MODE)=='TRUE')THEN
        PHYSICS_ON=.FALSE.
        WRITE(0,*)' Finalize without physics coupling. '
      ELSE
        PHYSICS_ON=.TRUE.
        WRITE(0,*)' Finalize with physics coupling. '
      ENDIF
!
!-----------------------------------------------------------------------
!***  Retrieve the Domain component's internal state.
!-----------------------------------------------------------------------
!
      CALL ESMF_GridCompGetInternalState(DOMAIN_GRID_COMP               &  !<-- The Domain component
                                        ,WRAP                           &  !<-- The F90 wrap of the Domain internal state
                                        ,RC)
!
      DOMAIN_INT_STATE=>wrap%DOMAIN_INT_STATE
!
!-----------------------------------------------------------------------
!***  Finalize the Solver subcomponent.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Finalize the Solver Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompFinalize(gridcomp   =domain_int_state%SOLVER_GRID_COMP &
                                ,importState=domain_int_state%IMP_STATE_SOLVER &
                                ,exportState=domain_int_state%EXP_STATE_SOLVER &
                                ,clock      =CLOCK_DOMAIN                &
                                ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!ratko    CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINALIZE)
! - FIX later
 RC=ESMF_SUCCESS
 RC_FINALIZE=ESMF_SUCCESS
!ratko
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Destroy all States.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Destroy States"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateDestroy(state=domain_int_state%IMP_STATE_SOLVER       &
                            ,rc   =RC)
!
      CALL ESMF_StateDestroy(state=domain_int_state%EXP_STATE_SOLVER       &
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINALIZE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  If quilting was selected for the generation of output,
!***  finalize and destroy objects related to the Write components.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      IF(domain_int_state%QUILTING)THEN
        CALL WRITE_DESTROY(DOMAIN_GRID_COMP                             &
                          ,DOMAIN_INT_STATE                             &
                          ,CLOCK_DOMAIN)
      ENDIF
!
!-----------------------------------------------------------------------
!***  Destroy the Domain Clock.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Destroy Domain Clock"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ClockDestroy(clock=CLOCK_DOMAIN                         &
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINALIZE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Destroy all subcomponents.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Destroy the Solver Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompDestroy(gridcomp=domain_int_state%SOLVER_GRID_COMP & 
                               ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINALIZE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  The final error signal information.
!-----------------------------------------------------------------------
!
      IF(RC_FINALIZE==ESMF_SUCCESS)THEN
        WRITE(0,*)'DOMAIN FINALIZE step succeeded'
      ELSE
        WRITE(0,*)'DOMAIN FINALIZE step failed'
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE DOMAIN_FINALIZE
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE DOMAIN_SETUP(MYPE_IN                                   &
                             ,MPI_INTRA                                 &
                             ,QUILTING                                  &
                             ,CF                                        &
                             ,DOMAIN_GRID_COMP                          &
                             ,DOMAIN_INT_STATE                          &
                             ,GRID_DOMAIN)
! 
!-----------------------------------------------------------------------
!***  This routine contains NMM-specific code for the Domain component:
!***    (1) Setting up distributed memory parallelism in the NMM;
!***    (2) Creating the ESMF Grid for the Domain components;
!***    (3) Sharing local subdomain index limits among tasks.
!-----------------------------------------------------------------------
!
      USE module_DOMAIN_INTERNAL_STATE
!
      USE module_DM_PARALLEL,ONLY : DECOMP                              &
                                   ,LOCAL_ISTART,LOCAL_IEND             &
                                   ,LOCAL_JSTART,LOCAL_JEND             &
                                   ,SETUP_SERVERS
!
      USE module_KINDS
!
!------------------------
!***  Argument variables
!------------------------
!
      INTEGER(kind=KINT),INTENT(IN) :: MYPE_IN                          &  !<-- Each MPI task's rank
                                      ,MPI_INTRA                           !<-- The communicator with the domain's fcst and quilt tasks.
!
      LOGICAL(kind=KLOG),INTENT(IN) :: QUILTING                            !<-- Has output via quilt tasks been specified?
!
      TYPE(ESMF_Config),INTENT(INOUT) :: CF                                !<-- This domain's configure object
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: DOMAIN_GRID_COMP                !<-- The Domain component
!
      TYPE(DOMAIN_INTERNAL_STATE),INTENT(INOUT) :: DOMAIN_INT_STATE        !<-- The Domain Internal State
!
      TYPE(ESMF_Grid),INTENT(OUT) :: GRID_DOMAIN                           !<-- The ESMF Grid for the NMM integration grid
!
!---------------------
!***  Local variables
!---------------------
!
      INTEGER(kind=KINT) :: I,IERR,J,K,N,NUM_PES,RC,RC_DOMAIN
!
      INTEGER(kind=KINT) :: IM,JM                                       &  !<-- Horizontal dimensions of the full integration grid
                           ,INPES,JNPES                                 &  !<-- MPI tasks in I and J directions
                           ,LM                                          &  !<-- Number of atmospheric model layers
                           ,MPI_INTRA_B                                 &  !<-- The MPI intra-communicator
                           ,MYPE                                        &  !<-- My MPI task ID
                           ,NUM_PES_FCST                                &  !<-- Number of MPI tasks applied to the forecast
                           ,NUM_PES_TOT                                 &  !<-- Total # of MPI tasks in the job
                           ,WRITE_GROUPS                                &  !<-- Number of groups of write tasks
                           ,WRITE_TASKS_PER_GROUP                          !<-- #of tasks in each write group
!
      INTEGER(kind=KINT),DIMENSION(2) :: I1                             &  !<-- # of I and J points in each fcst task's subdomain
                                        ,MIN,MAX                        &  !<-- Set start/end of each Grid dimension
                                        ,NCOUNTS                           !<-- Array with I/J limits of MPI task subdomains
!
      CHARACTER(50) :: MODE                                                !<-- Flag for global or regional run
!
      LOGICAL(kind=KLOG) :: GLOBAL                                         !<-- .TRUE. => global ; .FALSE. => regional
!
      TYPE(ESMF_VM) :: VM                                                  !<-- The ESMF virtual machine.
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      MYPE=MYPE_IN
!
!-----------------------------------------------------------------------
!***  Set up parameters for MPI communications on this domain's grid.
!-----------------------------------------------------------------------
!
      CALL MPI_COMM_SIZE(MPI_INTRA,NUM_PES_TOT,IERR)
!
      NUM_PES=NUM_PES_TOT
!
!-----------------------------------------------------------------------
!***  Establish the task layout including the Write tasks.
!***  The MPI communicator was provided as input and
!***  the forecast tasks in the I and J directions are
!***  extracted from a configure file.
!***  Give those to SETUP_SERVERS which will split the
!***  communicator between Forecast and Quilt/Write tasks.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="NMM_SETUP: Get INPES/JNPES from Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL MPI_COMM_DUP(MPI_INTRA,MPI_INTRA_B,RC)                          !<-- Use a duplicate of the communicator for safety
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &
                                  ,value =INPES                         &  !<-- # of fcst tasks in I direction
                                  ,label ='inpes:'                      &
                                  ,rc    =RC)
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &
                                  ,value =JNPES                         &  !<-- # of fcst tasks in J direction
                                  ,label ='jnpes:'                      &
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_DOMAIN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Set up Quilt/Write task specifications.
!***  First retrieve the task and group counts from the config file.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="NMM_SETUP: Get Write Task/Group Info from Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(CF                                   &  !<-- The configure file
                                  ,WRITE_GROUPS                         &  !<-- Number of write groups from config file
                                  ,label ='write_groups:'               &
                                  ,rc    =RC)
!
      CALL ESMF_ConfigGetAttribute(CF                                   &  !<-- The configure file
                                  ,WRITE_TASKS_PER_GROUP                &  !<-- Number of write tasks per group from config file
                                  ,label ='write_tasks_per_group:'      &
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_DOMAIN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Segregate the Forecast tasks from the Quilt/Write tasks.
!-----------------------------------------------------------------------
!
      CALL SETUP_SERVERS(MYPE,INPES,JNPES,NUM_PES                       &
                        ,WRITE_GROUPS,WRITE_TASKS_PER_GROUP             &
                        ,MPI_INTRA_B,QUILTING)
!
!***
!***  NOTE: At this point, NUM_PES is the number of Forecast tasks only.
!***
!-----------------------------------------------------------------------
!
      NUM_PES_FCST=INPES*JNPES                                             !<-- Number of forecast tasks
      domain_int_state%NUM_PES_FCST=NUM_PES_FCST                           !<-- Save this for DOMAIN's Run step
!
!-----------------------------------------------------------------------
!***  Allocate and fill the task list that holds the IDs of
!***  the Forecast tasks.
!-----------------------------------------------------------------------
!
      ALLOCATE(domain_int_state%PETLIST_FCST(NUM_PES_FCST))                !<-- Task IDs of the forecast tasks
!
      DO N=0,NUM_PES_FCST-1
        domain_int_state%PETLIST_FCST(N+1)=N                               !<-- Collect just the forecast task IDs
      ENDDO
!
!-----------------------------------------------------------------------
!***  Retrieve the VM (Virtual Machine) of the Domain component.
!***  We need VM now to set up the DE layout.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="NMM_SETUP: Retrieve VM from Domain Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompGet(gridcomp=DOMAIN_GRID_COMP                   &  !<-- The Domain component
                           ,vm      =VM                                 &  !<-- The ESMF Virtual Machine
                           ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_DOMAIN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Create DE layout based on the I tasks by J tasks specified in
!***  the config file.
!***  This refers only to Forecast tasks.
!-----------------------------------------------------------------------
!
!d      IF(MYPE<NUM_PES_FCST)THEN                                            !<-- Select only the forecast tasks
!d        MY_DE_LAYOUT=ESMF_DELayoutCreate(            VM                 &  !<-- The ESMF virtual machine
!d                                        ,deCountList=(/INPES,JNPES/)    &  !<-- User-specified I-task by J-task layout
!d                                        ,rc         =RC)
!d      ENDIF
!
!-----------------------------------------------------------------------
!***  Create the ESMF Grid.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Extract the dimensions of the domain from the configure file.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="NMM_SETUP: Get IM,JM,LM from Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure file
                                  ,value =IM                            &  !<-- I dimension of full domain
                                  ,label ='im:'                         &  !<-- The label in the configure file
                                  ,rc    =RC)
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure file
                                  ,value =JM                            &  !<-- J dimension of full domain
                                  ,label ='jm:'                         &  !<-- The label in the configure file
                                  ,rc    =RC)
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure file
                                  ,value =LM                            &  !<-- Vertical dimension of full domain
                                  ,label ='lm:'                         &  !<-- The label in the configure file
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_DOMAIN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!------------------------------------------------------------
!***  Retrieve the forecast domain mode from the config file.
!------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="NMM_SETUP: Get GLOBAL/REGIONAL Mode from Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure file
                                  ,value =MODE                          &  !<-- Flag for global (true) or regional (false) run
                                  ,label ='global:'                     &  !<-- The label in the configure file
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_DOMAIN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF(TRIM(MODE)=='true')THEN
        GLOBAL=.TRUE.
      ELSE
        GLOBAL=.FALSE.
      ENDIF
!
!-----------------------------------------------------------------------
!***  If this is a global mode forecast, extend IM and JM.
!***  The first dimension of NCOUNTS is the I dimension for parallelization.
!***  The second dimension of NCOUNTS is the J dimension.
!-----------------------------------------------------------------------
!
      IF(GLOBAL)THEN      !<-- Global mode horizontal dimensions.
        NCOUNTS(1)=IM+2
        NCOUNTS(2)=JM+2
      ELSE                !<-- Regional mode horizontal dimensions.
        NCOUNTS(1)=IM
        NCOUNTS(2)=JM
      ENDIF
!
      MAX(1)=NCOUNTS(1)
      MAX(2)=NCOUNTS(2)
!
      MIN(1)=1
      MIN(2)=1
!
!-----------------------------------------------------------------------
!***  Now create the Domain component's ESMF Grid
!***  for the NMM's integration grid.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="NMM_SETUP: Create the ESMF Grid"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      GRID_DOMAIN=ESMF_GRIDCREATE      (regDecomp     =(/INPES,JNPES/)    &  !<-- I x J task layout
                                       ,minIndex      =(/MIN(1),MIN(2)/)  &  !<-- Min indices in I and J
                                       ,maxIndex      =(/MAX(1),MAX(2)/)  &  !<-- Max indices in I and J
                                       ,gridEdgeLWidth=(/0,0/)            &  !<-- Padding, lower edges for noncentered stagger
                                       ,gridEdgeUWidth=(/0,0/)            &  !<-- Padding, upper edges for noncentered stagger
                                       ,name          ="GRID"             &  !<-- Name of the Grid
                                       ,indexflag     =ESMF_INDEX_GLOBAL  &
                                       ,rc            =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_DOMAIN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Get the local array sizes for the Domain Grid.
!***  Only forecast tasks are relevant here.
!-----------------------------------------------------------------------
!
      IF(MYPE<NUM_PES_FCST)THEN                                              !<-- Select only fcst tasks
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="NMM_SETUP: Get EMSF Sizes of Local Subdomains"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_GridGet(grid              =GRID_DOMAIN                &
                         ,localDe           =0                          &
                         ,staggerloc        =ESMF_STAGGERLOC_CENTER     &
                         ,computationalCount=I1                         & !<-- # of local points in I and J on each task
                         ,rc                =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_DOMAIN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF
!
!-----------------------------------------------------------------------
!***  Using 'computationalCount' from array I1 obtained in the
!***  previous call, generate all of the local task index limits
!***  for all Forecast tasks.  
!***  The user, not ESMF, does this work.
!-----------------------------------------------------------------------
!
      IF(MYPE<NUM_PES_FCST)THEN                                            !<-- Select only the forecast tasks
        CALL DECOMP(MYPE,INPES,JNPES,NUM_PES_FCST,IM,JM,LM,GLOBAL,I1)
!
        ALLOCATE(domain_int_state%LOCAL_ISTART(0:NUM_PES_FCST-1))
        ALLOCATE(domain_int_state%LOCAL_IEND  (0:NUM_PES_FCST-1))
        ALLOCATE(domain_int_state%LOCAL_JSTART(0:NUM_PES_FCST-1))
        ALLOCATE(domain_int_state%LOCAL_JEND  (0:NUM_PES_FCST-1))
!
        DO N=0,NUM_PES_FCST-1
          domain_int_state%LOCAL_ISTART(N)=LOCAL_ISTART(N)                 !<-- Starting I for all forecast tasks' subdomains
          domain_int_state%LOCAL_IEND  (N)=LOCAL_IEND  (N)                 !<-- Ending I for all forecast tasks' subdomains
          domain_int_state%LOCAL_JSTART(N)=LOCAL_JSTART(N)                 !<-- Starting J for all forecast tasks' subdomains
          domain_int_state%LOCAL_JEND  (N)=LOCAL_JEND  (N)                 !<-- Ending J for all forecast tasks' subdomains
        ENDDO
!
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE DOMAIN_SETUP
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE NMM_FILTERING(DOMAIN_GRID_COMP                         &
                              ,IMP_STATE                                &
                              ,EXP_STATE                                &
                              ,CLOCK_DOMAIN                             &
                              ,RC_FILT)
!
!-----------------------------------------------------------------------
!***  Phase 2 of the Run step of the Domain component.
!***  This phase is only relevant when digital filtering is
!***  in effect and executes at the end of each timestep
!***  after the Solver (in phase 1). 
!
!***  Called from subroutine NMM_INTEGRATE.
!-----------------------------------------------------------------------
!
      USE module_DIGITAL_FILTER_NMM
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(ESMF_GridComp) :: DOMAIN_GRID_COMP                              !<-- The Domain gridded component
!
      TYPE(ESMF_State) :: IMP_STATE                                     &  !<-- The Domain import state
                         ,EXP_STATE                                        !<-- The Domain export state
!
      TYPE(ESMF_Clock) :: CLOCK_DOMAIN                                     !<-- The Domain ESMF Clock
!
      INTEGER,INTENT(OUT) :: RC_FILT                                       !<-- Return code for this step
!
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: DFIHR,FILTER_METHOD,MEAN_ON,NDFISTEP        &
                           ,NUM_TRACERS_CHEM
!
      INTEGER(kind=KINT) :: YY, MM, DD, H, M, S
!
      INTEGER(kind=KINT) :: RC
!
      CHARACTER(8) :: CLOCK_DIRECTION
!
      TYPE(ESMF_Time) :: CURRTIME                                       &  !<-- The current time of Clock_DOMAIN
                        ,STARTTIME                                      &  !<-- The start time of Clock_DOMAIN
                        ,TESTTIME
!
      TYPE(ESMF_TimeInterval) :: DT_ESMF
!
      INTEGER(ESMF_KIND_I4) :: INTEGER_DT,NUMERATOR_DT,IDENOMINATOR_DT
!
      TYPE(DOMAIN_INTERNAL_STATE),POINTER,SAVE :: DOMAIN_INT_STATE         !<-- The Domain internal state pointer
!
      TYPE(WRAP_DOMAIN_INTERNAL_STATE),SAVE :: WRAP                        !<-- The F90 wrap of the Domain internal state
!
      TYPE(SOLVER_INTERNAL_STATE),POINTER :: SOLVER_INT_STATE              !<-- The Solver internal state pointer
!
      TYPE(WRAP_SOLVER_INT_STATE) :: WRAP_SOLVER                           !<-- The F90 wrap of the Solver internal state
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC     =ESMF_SUCCESS
      RC_FILT=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  What is this domain's ID?
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="NMM_FILTERING: Extract Domain ID from Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=IMP_STATE                            &  !<-- The Domain import state
                            ,name ='DOMAIN_ID'                          &  !<-- Name of the attribute to extract
                            ,value=MY_DOMAIN_ID                         &  !<-- The ID of this domain
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FILT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Extract the Domain internal state.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
      MESSAGE_CHECK="NMM_Filtering: Extract the Domain Internal State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
      CALL ESMF_GridCompGetInternalState(DOMAIN_GRID_COMP               &  !<-- The Domain component
                                        ,WRAP                           &  !<-- The F90 wrap of the Domain internal state
                                        ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FILT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
      DOMAIN_INT_STATE=>wrap%DOMAIN_INT_STATE
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  What are the start time and the current time?
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="NMM_FILTERING: Extract StartTime,CurrentTime"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ClockGet(clock    =CLOCK_DOMAIN                        &
                        ,startTime=STARTTIME                           &
                        ,currTime =CURRTIME                            &
                        ,timeStep =DT_ESMF                             &
                        ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FILT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="NMM_FILTERING: Get Actual Timestep from ESMF Variable"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_TimeIntervalGet(timeinterval=DT_ESMF                   &  !<-- the ESMF timestep
                               ,s           =INTEGER_DT                &  !<-- the integer part of the timestep in seconds
                               ,sN          =NUMERATOR_DT              &  !<-- the numerator of the fractional second
                               ,sD          =IDENOMINATOR_DT           &  !<-- the denominator of the fractional second
                               ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FILT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  What is the Clock direction?
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="NMM_FILTERING: Extract Clock Direction."
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=IMP_STATE                            &  !<-- Extract the direction of the Clock from the import state
                            ,name ='Clock_Direction'                    &
                            ,value=CLOCK_DIRECTION                      &
                            ,rc   =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FILT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="NMM_FILTERING: Extract Mean_On Flag from Imp State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=IMP_STATE                            &  !<-- Extract MEAN_ON flag from import state
                            ,name ='MEAN_ON'                            &
                            ,value=MEAN_ON                              &
                            ,rc   =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FILT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="NMM_FILTERING: Extract Filter Method from Imp State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=IMP_STATE                            &  !<- Extract FILTER_METHOD flag from import state
                            ,name ='Filter_Method'                      &
                            ,value=FILTER_METHOD                        &
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FILT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
      fcst_tasks: IF(MYPE<domain_int_state%NUM_PES_FCST)THEN               !<-- Only forecast tasks deal with the integration
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Some subroutines called here will in turn call halo exchange
!***  routines.  The halo exchange routines need to know a set of
!***  15 domain-related variables.  Set those variables now.
!***  They all lie in the Solver internal state so first extract
!***  that object.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Extract the Solver internal state.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="NMM_FILT: Extract Solver Internal State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_GridCompGetInternalState(domain_int_state%SOLVER_GRID_COMP &  !<-- The Solver component
                                          ,WRAP_SOLVER                       &
                                          ,RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FILT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        SOLVER_INT_STATE=>wrap_solver%INT_STATE
!
!-----------------------
!***  The initial stage
!-----------------------
!
        IF(CURRTIME==STARTTIME .and. domain_int_state%FIRST_FILTER)THEN
          domain_int_state%FIRST_FILTER=.FALSE.
          CALL DIGITAL_FILTER_PHY_INIT_NMM(domain_int_state%FILT_BUNDLE_RESTORE   &
                                          ,solver_int_state%ITS                   &
                                          ,solver_int_state%ITE                   &
                                          ,solver_int_state%JTS                   &
                                          ,solver_int_state%JTE                   &
                                          ,solver_int_state%LM                    &
                                          ,domain_int_state%NUM_FIELDS_RESTORE_2D &
                                          ,domain_int_state%NUM_FIELDS_RESTORE_3D &
                                          ,domain_int_state%SAVE_2D_PHYS          &
                                          ,domain_int_state%SAVE_3D_PHYS )
        ENDIF

        IF(CURRTIME==STARTTIME)THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Extract NDFISTEP from Domain Import State"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeGet(state=IMP_STATE                        &  !<-- Extract the filter value NDFISTEP from the import state
                                ,name ='NDFISTEP'                       &
                                ,value=NDFISTEP                         &
                                ,rc   =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FILT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          INTEGER_DT=ABS(INTEGER_DT)
!
          CALL DIGITAL_FILTER_DYN_INIT_NMM(domain_int_state%FILT_BUNDLE_FILTER   &
                                          ,NDFISTEP                              &
                                          ,INTEGER_DT                            &
                                          ,NUMERATOR_DT                          &
                                          ,IDENOMINATOR_DT                       &
                                          ,solver_int_state%ITS                  &
                                          ,solver_int_state%ITE                  &
                                          ,solver_int_state%JTS                  &
                                          ,solver_int_state%JTE                  &
                                          ,solver_int_state%LM                   &
                                          ,domain_int_state%NUM_FIELDS_FILTER_2D &
                                          ,domain_int_state%NUM_FIELDS_FILTER_3D &
                                          ,domain_int_state%KSTEP                &
                                          ,domain_int_state%NSTEP                &
                                          ,domain_int_state%TOTALSUM             &
                                          ,domain_int_state%DOLPH_WGTS           &
                                          ,domain_int_state%SAVE_2D              &
                                          ,domain_int_state%SAVE_3D )
!
        ENDIF
!
!-----------------------------------------------------------------------
        direction: IF(CLOCK_DIRECTION(1:7)=='Forward')THEN
!-----------------------------------------------------------------------
!
!-------------------------
!***  The summation stage
!-------------------------
!
          startdef: IF(CURRTIME == STARTTIME)THEN
!
            DFIHR=NDFISTEP*(INTEGER_DT+(float(NUMERATOR_DT)/IDENOMINATOR_DT))
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Set HALFDFIINTVAL in Summation State"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_TimeIntervalSet(timeinterval=HALFDFIINTVAL        &
                                     ,s           =DFIHR                &
                                     ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FILT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            HALFDFITIME=CURRTIME+HALFDFIINTVAL
            DFITIME=HALFDFITIME+HALFDFIINTVAL
!
          ENDIF startdef
!
!---------------------
!
          IF(CURRTIME>=STARTTIME)THEN
            CALL DIGITAL_FILTER_DYN_SUM_NMM(domain_int_state%FILT_BUNDLE_FILTER   &
                                           ,MEAN_ON                               &
                                           ,solver_int_state%ITS                  &
                                           ,solver_int_state%ITE                  &
                                           ,solver_int_state%JTS                  &
                                           ,solver_int_state%JTE                  &
                                           ,solver_int_state%LM                   &
                                           ,domain_int_state%NUM_FIELDS_FILTER_2D &
                                           ,domain_int_state%NUM_FIELDS_FILTER_3D &
                                           ,domain_int_state%KSTEP                &
                                           ,domain_int_state%NSTEP                &
                                           ,domain_int_state%TOTALSUM             &
                                           ,domain_int_state%DOLPH_WGTS           &
                                           ,domain_int_state%SAVE_2D              &
                                           ,domain_int_state%SAVE_3D )
          ENDIF
!
!---------------------
!
          IF(CURRTIME==HALFDFITIME .AND. FILTER_METHOD == 1)THEN
            CALL DIGITAL_FILTER_PHY_SAVE_NMM(domain_int_state%FILT_BUNDLE_RESTORE   &
                                            ,solver_int_state%ITS                   &
                                            ,solver_int_state%ITE                   &
                                            ,solver_int_state%JTS                   &
                                            ,solver_int_state%JTE                   &
                                            ,domain_int_state%NUM_FIELDS_RESTORE_2D &
                                            ,domain_int_state%NUM_FIELDS_RESTORE_3D &
                                            ,domain_int_state%SAVE_2D_PHYS          &
                                            ,domain_int_state%SAVE_3D_PHYS )
          ENDIF
!
!---------------------
!***  The final stage
!---------------------
!
          TESTTIME=CURRTIME+DT_ESMF
!
          IF(TESTTIME==DFITIME)THEN
!
            CALL DIGITAL_FILTER_DYN_AVERAGE_NMM(domain_int_state%FILT_BUNDLE_FILTER   &
                                               ,solver_int_state%ITS                  &
                                               ,solver_int_state%ITE                  &
                                               ,solver_int_state%JTS                  &
                                               ,solver_int_state%JTE                  &
                                               ,solver_int_state%LM                   &
                                               ,domain_int_state%NUM_FIELDS_FILTER_2D &
                                               ,domain_int_state%NUM_FIELDS_FILTER_3D &
                                               ,domain_int_state%KSTEP                &
                                               ,domain_int_state%NSTEP                &
                                               ,domain_int_state%TOTALSUM             &
                                               ,domain_int_state%SAVE_2D              &
                                               ,domain_int_state%SAVE_3D )
!
            CALL DIGITAL_FILTER_PHY_RESTORE_NMM(domain_int_state%FILT_BUNDLE_RESTORE   &
                                               ,solver_int_state%ITS                   &
                                               ,solver_int_state%ITE                   &
                                               ,solver_int_state%JTS                   &
                                               ,solver_int_state%JTE                   &
                                               ,domain_int_state%NUM_FIELDS_RESTORE_2D &
                                               ,domain_int_state%NUM_FIELDS_RESTORE_3D &
                                               ,domain_int_state%SAVE_2D_PHYS          &
                                               ,domain_int_state%SAVE_3D_PHYS )
!
            CALL ESMF_ClockPrint(clock  =CLOCK_DOMAIN                   &
                                ,options="currtime string"              &
                                ,rc     =RC)
          ENDIF
!
!-----------------------------------------------------------------------
        ELSEIF(CLOCK_DIRECTION(1:7)=='Bckward')THEN
!-----------------------------------------------------------------------
!
          IF(CURRTIME == STARTTIME)THEN

            DFIHR=NDFISTEP*(ABS(INTEGER_DT)                             &
                           +ABS(REAL(NUMERATOR_DT)/IDENOMINATOR_DT) )
!
            CALL ESMF_TimeIntervalSet(timeinterval=HALFDFIINTVAL        &
                                    ,s            =DFIHR                &
                                    ,rc           =RC)
            CALL DIGITAL_FILTER_PHY_SAVE_NMM(domain_int_state%FILT_BUNDLE_RESTORE   &
                                            ,solver_int_state%ITS                   &
                                            ,solver_int_state%ITE                   &
                                            ,solver_int_state%JTS                   &
                                            ,solver_int_state%JTE                   &
                                            ,domain_int_state%NUM_FIELDS_RESTORE_2D &
                                            ,domain_int_state%NUM_FIELDS_RESTORE_3D &
                                            ,domain_int_state%SAVE_2D_PHYS          &
                                            ,domain_int_state%SAVE_3D_PHYS )
!
            HALFDFITIME=CURRTIME-HALFDFIINTVAL
            DFITIME=HALFDFITIME-HALFDFIINTVAL
!
          ENDIF
!
!-------------------------
!***  The summation stage
!-------------------------
!
          IF(CURRTIME<=STARTTIME)THEN
            CALL DIGITAL_FILTER_DYN_SUM_NMM(domain_int_state%FILT_BUNDLE_FILTER   &
                                           ,MEAN_ON                               &
                                           ,solver_int_state%ITS                  &
                                           ,solver_int_state%ITE                  &
                                           ,solver_int_state%JTS                  &
                                           ,solver_int_state%JTE                  &
                                           ,solver_int_state%LM                   &
                                           ,domain_int_state%NUM_FIELDS_FILTER_2D &
                                           ,domain_int_state%NUM_FIELDS_FILTER_3D &
                                           ,domain_int_state%KSTEP                &
                                           ,domain_int_state%NSTEP                &
                                           ,domain_int_state%TOTALSUM             &
                                           ,domain_int_state%DOLPH_WGTS           &
                                           ,domain_int_state%SAVE_2D              &
                                           ,domain_int_state%SAVE_3D )
          ENDIF
!
!---------------------
!***  The final stage
!---------------------
!
          TESTTIME=CURRTIME+DT_ESMF
!
          IF(TESTTIME==DFITIME)THEN
            IF (FILTER_METHOD == 3) THEN
            CALL DIGITAL_FILTER_DYN_AVERAGE_NMM(domain_int_state%FILT_BUNDLE_FILTER   &
                                               ,solver_int_state%ITS                  &
                                               ,solver_int_state%ITE                  &
                                               ,solver_int_state%JTS                  &
                                               ,solver_int_state%JTE                  &
                                               ,solver_int_state%LM                   &
                                               ,domain_int_state%NUM_FIELDS_FILTER_2D &
                                               ,domain_int_state%NUM_FIELDS_FILTER_3D &
                                               ,domain_int_state%KSTEP                &
                                               ,domain_int_state%NSTEP                &
                                               ,domain_int_state%TOTALSUM             &
                                               ,domain_int_state%SAVE_2D              &
                                               ,domain_int_state%SAVE_3D )
            ENDIF
!
            CALL ESMF_ClockPrint(clock  =CLOCK_DOMAIN                   &
                                ,options="currtime string"              &
                                ,rc     =RC)
          ENDIF
!
!-----------------------------------------------------------------------
!
        ENDIF direction
!
      ENDIF fcst_tasks
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE NMM_FILTERING
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE CALL_WRITE_ASYNC(DOMAIN_GRID_COMP                      &
                                 ,IMP_STATE                             &
                                 ,EXP_STATE                             &
                                 ,CLOCK_DOMAIN                          &
                                 ,RC_RUN2)
!
!-----------------------------------------------------------------------
!***  Phase 3 of the Run step of the NMM Domain component.
!***  It initiates the writing of history/restart files
!***  from each Domain component.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(ESMF_GridComp) :: DOMAIN_GRID_COMP                              !<-- The Domain gridded component
!
      TYPE(ESMF_State) :: IMP_STATE                                     &  !<-- The Domain component's import state
                         ,EXP_STATE                                        !<-- The Domain component's export state
!
      TYPE(ESMF_Clock) :: CLOCK_DOMAIN                                     !<-- The Domain ESMF Clock
!
      TYPE(ESMF_Time) :: CURRTIME                                       &  !<-- The ESMF current time.
                        ,STOPTIME                                          !<-- The ESMF start time.
!
      INTEGER,INTENT(OUT) :: RC_RUN2                                       !<-- Return code for the Run step 
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: LENGTH,N,NB,RC
!
      TYPE(DOMAIN_INTERNAL_STATE),POINTER :: DOMAIN_INT_STATE              !<-- The Domain internal state pointer
!
      TYPE(WRAP_DOMAIN_INTERNAL_STATE) :: WRAP                             !<-- The F90 wrap of the Domain internal state
!
      CHARACTER(ESMF_MAXSTR) :: CWRT
!
      LOGICAL(kind=KINT) :: LAST_TIME                                      !<-- Test time logical
!
      LOGICAL(kind=KLOG) :: I_AM_A_FCST_TASK
!
      TYPE(SOLVER_INTERNAL_STATE),POINTER :: SOLVER_INT_STATE
!
      TYPE(WRAP_SOLVER_INT_STATE) :: WRAP_SOLVER
!
      type(esmf_time) :: ringtime
      type(esmf_timeinterval) :: ringinterval
      integer :: month,day
      integer(esmf_kind_i4) :: yy,mm,dd,h,m,s,sn,sd
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC     =ESMF_SUCCESS
      RC_RUN2=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  Write a history file at the end of the appropriate timesteps.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  First retrieve the Domain component's internal state.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
      MESSAGE_CHECK="Run2: Retrieve Domain Component's Internal State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
      CALL ESMF_GridCompGetInternalState(DOMAIN_GRID_COMP               &  !<-- The Domain gridded component
                                        ,WRAP                           &  !<-- The F90 wrap of the Domain internal state
                                        ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN2)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
      DOMAIN_INT_STATE=>wrap%DOMAIN_INT_STATE
!
!-----------------------------------------------------------------------
!***  If quilting was not specified then exit.
!-----------------------------------------------------------------------
!
      IF(.NOT.domain_int_state%QUILTING)THEN
        RETURN
      ENDIF
!
      call esmf_gridcompget(gridcomp=domain_grid_comp  &
                           ,vm=vm &
                           ,rc=rc)
      call esmf_vmget(vm=vm  &
                     ,localpet=mype  &
                     ,rc=rc)
!-----------------------------------------------------------------------
!***  Is this a forecast task?
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Get Fcst-or-Write Task Flag from the Domain Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=EXP_STATE                            &  !<-- The Domain component export state
                            ,name ='Fcst-or-Write Flag'                 &  !<-- Use this name inside the state
                            ,value=I_AM_A_FCST_TASK                     &  !<-- The logical being inserted into the import state
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN2)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Check to see if the history alarm is ringing and if so then
!***  prepare to call the Write subroutine WRITE_ASYNC.
!-----------------------------------------------------------------------
!
!     write(0,80551)ESMF_AlarmIsRinging(alarm=domain_int_state%ALARM_HISTORY,rc=rc)
80551 format(' CALL_WRITE_ASYNC ALARM_HISTORY AlarmIsRinging=',l1)
!     call esmf_alarmget(alarm=domain_int_state%ALARM_HISTORY &
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
!     write(0,*)' ALARM_HISTORY ringtime: y=',yy,' mm=',month,' dd=',dd,' h=',h &
!              ,' m=',m,' s=',s,' sn=',sn,'s d=',sd
!     call esmf_timeintervalget(timeinterval=ringinterval &
!                      ,m=m &
!                      ,s=s &
!                      ,sn=sn &
!                      ,sd=sd)
!     write(0,*)' ringinterval: m=',m,' s=',s,' sn=',sn,'s d=',sd
      alarms: IF(ESMF_AlarmIsRinging(alarm=domain_int_state%ALARM_HISTORY  &  !<-- The history output alarm
                                    ,rc   =RC)                             &
                                .OR.                                       &
                 ESMF_AlarmIsRinging(alarm=domain_int_state%ALARM_RESTART  &  !<-- The restart output alarm
                                    ,rc   =RC))THEN
!
!     write(0,80552)
80552 format(' CALL_WRITE_ASYNC ALARM_HISTORY AlarmIsRinging= T ')
!-----------------------------------------------------------------------
!***  Extract the Solver internal state in order to access 
!***  output-related variables within it.
!-----------------------------------------------------------------------
!
        fcst_tasks: IF(I_AM_A_FCST_TASK)THEN
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="CALL_WRITE_ASYNC: Extract Solver Internal State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_GridCompGetInternalState(domain_int_state%SOLVER_GRID_COMP &  !<-- The Solver component
                                            ,WRAP_SOLVER                       &
                                            ,RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN2)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          SOLVER_INT_STATE=>wrap_solver%INT_STATE
!
!-----------------------------------------------------------------------
!***  Refresh the values of the ESMF Attributes in the history/restart
!***  Bundles because Attributes are not updated automatically as Fields
!***  are.
!-----------------------------------------------------------------------
!
          all_vars: DO N=1,solver_int_state%NUM_VARS                       !<-- Loop through all output variables.
!
!---------------------
!***  Integer scalars
!---------------------
!
            IF (solver_int_state%VARS(N)%TKR == TKR_I0D) THEN              !<-- Select integer scalars
!
              IF (solver_int_state%VARS(N)%HISTORY) THEN                   !<-- Is integer scalar specified for history output?
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="CALL_WRITE_ASYNC: Update Integer Scalar in History Bundle"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
                CALL ESMF_AttributeSet(FIELDBUNDLE=solver_int_state%BUNDLE_ARRAY(1)   &  !<-- The output Bundle
                                      ,name       =solver_int_state%VARS(N)%VBL_NAME  &  !<-- Name of the integer scalar
                                      ,value      =solver_int_state%VARS(N)%I0D       &  !<-- The scalar inserted into the Bundle
                                      ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN2)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
              ENDIF
!
              IF(solver_int_state%VARS(N)%RESTART)THEN                     !<-- Is integer scalar specified for restart output
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="CALL_WRITE_ASYNC: Update Integer Attribute in Restart Bundle"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
                CALL ESMF_AttributeSet(FIELDBUNDLE=solver_int_state%BUNDLE_ARRAY(2)  &  !<-- The restart Bundle
                                      ,name       =solver_int_state%VARS(N)%VBL_NAME &  !<-- Name of the integer scalar
                                      ,value      =solver_int_state%VARS(N)%I0D      &  !<-- The scalar being inserted into the import state
                                      ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
                CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN2)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
              ENDIF
!
            ENDIF
!
!------------------------
!***  1-D integer arrays
!------------------------
!
            IF (solver_int_state%VARS(N)%TKR == TKR_I1D) THEN              !<-- Select 1-D integer arrays
!
              LENGTH=SIZE(solver_int_state%VARS(N)%I1D)
!
              IF (solver_int_state%VARS(N)%HISTORY) THEN                   !<-- Is the array specified for history output?
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="CALL_WRITE_ASYNC: Update 1-D Integer Array in History Bundle"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
                CALL ESMF_AttributeSet(FIELDBUNDLE=solver_int_state%BUNDLE_ARRAY(1)  &  !<-- The history Bundle
                                      ,name       =solver_int_state%VARS(N)%VBL_NAME &  !<-- Name of the 1-D integer array
                                      ,itemCount  =LENGTH                            &  !<-- # of elements in this attribute
                                      ,valueList  =solver_int_state%VARS(N)%I1D      &  !<-- The array being inserted into the import state
                                      ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
                CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN2)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
              ENDIF
!
              IF(solver_int_state%VARS(N)%RESTART)THEN                     !<-- Is 1-D integer array specified for restart output?
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="CALL_WRITE_ASYNC: Update 1-D Integer Array in Restart Bundle"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
                CALL ESMF_AttributeSet(FIELDBUNDLE=solver_int_state%BUNDLE_ARRAY(2)  &  !<-- The restart Bundle
                                      ,name       =solver_int_state%VARS(N)%VBL_NAME &  !<-- Name of the integer scalar
                                      ,itemCount  =LENGTH                            &  !<-- # of elements in this attribute
                                      ,valueList  =solver_int_state%VARS(N)%I1D      &  !<-- The array being inserted into the import state
                                      ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
                CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN2)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
              ENDIF
!
            ENDIF
!
!------------------
!***  Real scalars
!------------------
!
            IF(solver_int_state%VARS(N)%TKR == TKR_R0D)THEN                !<-- Select real scalars
!
              IF(solver_int_state%VARS(N)%HISTORY)THEN                     !<-- Is real scalar specified for history output?
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="CALL_WRITE_ASYNC: Update Real Scalar in History Bundle"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
                CALL ESMF_AttributeSet(FIELDBUNDLE=solver_int_state%BUNDLE_ARRAY(1)  &  !<-- The history Bundle
                                      ,name       =solver_int_state%VARS(N)%VBL_NAME &  !<-- Name of the real scalar
                                      ,value      =solver_int_state%VARS(N)%R0D      &  !<-- The scalar being inserted into the import state
                                      ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
                CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN2)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
              ENDIF
!
              IF(solver_int_state%VARS(N)%RESTART)THEN                     !<-- Is real scalar specified for restart output
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="CALL_WRITE_ASYNC: Update Real Scalar in Restart Bundle"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
                CALL ESMF_AttributeSet(FIELDBUNDLE=solver_int_state%BUNDLE_ARRAY(2)  &  !<-- The restart Bundle
                                      ,name       =solver_int_state%VARS(N)%VBL_NAME &  !<-- Name of the real scalar
                                      ,value      =solver_int_state%VARS(N)%R0D      &  !<-- The scalar being inserted into the import state
                                      ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
                CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN2)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
              ENDIF
!
            ENDIF
!
!---------------------
!***  1-D real arrays
!---------------------
!
            IF(solver_int_state%VARS(N)%TKR == TKR_R1D)THEN                !<-- Select 1-D real arrays
!
              LENGTH=SIZE(solver_int_state%VARS(N)%R1D)
!
              IF(solver_int_state%VARS(N)%HISTORY)THEN                     !<-- Is the array specified for history output?
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="CALL_WRITE_ASYNC: Update 1-D Real Array in History Bundle"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
                CALL ESMF_AttributeSet(FIELDBUNDLE=solver_int_state%BUNDLE_ARRAY(1)  &  !<-- The history Bundle
                                      ,name       =solver_int_state%VARS(N)%VBL_NAME &  !<-- Name of the 1-D real array
                                      ,itemCount  =LENGTH                            &  !<-- # of elements in this attribute
                                      ,valueList  =solver_int_state%VARS(N)%R1D      &  !<-- The array being inserted into the import state
                                      ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
                CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN2)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
              ENDIF
!
              IF(solver_int_state%VARS(N)%RESTART)THEN                     !<-- Is the array specified for restart output?
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="CALL_WRITE_ASYNC: Update 1-D Real Array in Restart Bundle"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
                CALL ESMF_AttributeSet(FIELDBUNDLE=solver_int_state%BUNDLE_ARRAY(2)  &  !<-- The restart Bundle
                                      ,name       =solver_int_state%VARS(N)%VBL_NAME &  !<-- Name of the 1-D real array
                                      ,itemCount  =LENGTH                            &  !<-- # of elements in this attribute
                                      ,valueList  =solver_int_state%VARS(N)%R1D      &  !<-- The array being inserted into the import state
                                      ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
                CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN2)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
              ENDIF
!
            ENDIF
!
!-----------------------------------------------------------------------
!
          ENDDO all_vars
!
!-----------------------------------------------------------------------
!
        ENDIF fcst_tasks
!
!-----------------------------------------------------------------------
!
      ELSE
!
        RETURN                                                             !<-- No output alarm is ringing
!
      ENDIF alarms
!
!-----------------------------------------------------------------------
!***  Check to see if the history alarm is ringing and if so then
!***  call the Write subroutine WRITE_ASYNC to execute the writing
!***  of a history file.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="CALL_WRITE_ASYNC: Is ALARM_HISTORY ringing?"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF(ESMF_AlarmIsRinging(alarm=domain_int_state%ALARM_HISTORY  &  !<-- The history output alarm
                            ,rc   =RC))THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN2)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!     call print_memory()
        IF(domain_int_state%QUILTING)THEN
          CWRT='History'
          CALL WRITE_ASYNC(DOMAIN_GRID_COMP                             &
                          ,DOMAIN_INT_STATE                             &
                          ,CLOCK_DOMAIN                                 &
                          ,MYPE                                         &
                          ,CWRT)
        ENDIF
!     call print_memory()
!
!-----------------------------------------------------------------------
!
      ENDIF
!
!-----------------------------------------------------------------------
!***  Write a restart file at the end of the appropriate timesteps.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="CALL_WRITE_ASYNC: Is ALARM_RESTART ringing?"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF(ESMF_AlarmIsRinging(alarm=domain_int_state%ALARM_RESTART       &  !<-- The restart output alarm
                            ,rc   =RC))THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN2)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Retrieve the Domain component's internal state.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
        MESSAGE_CHECK="Run2: Retrieve Domain Component's Internal State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
        CALL ESMF_GridCompGetInternalState(DOMAIN_GRID_COMP             &  !<-- The Domain gridded component
                                          ,WRAP                         &  !<-- The F90 wrap of the Domain internal state
                                          ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN2)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
        DOMAIN_INT_STATE=>wrap%DOMAIN_INT_STATE
!
!-----------------------------------------------------------------------
!***  Execute the writing of a restart file.
!-----------------------------------------------------------------------
!
        CALL ESMF_ClockGet(clock    =CLOCK_DOMAIN                       &  !<-- The ESMF Clock of this domain
                          ,stopTime =STOPTIME                           &  !<-- The simulation stop time
                          ,currTime =CURRTIME                           &  !<-- Current time of simulation
                          ,rc       =RC)
!
        LAST_TIME = (STOPTIME==CURRTIME)                                   !<-- Is it last write step?
!
        IF(domain_int_state%QUILTING                                    &
              .AND.                                                     &
           (.NOT.LAST_TIME                                              &
              .OR.                                                      &
           domain_int_state%WRITE_LAST_RESTART) )THEN

          CWRT='Restart'
          CALL WRITE_ASYNC(DOMAIN_GRID_COMP                             &
                          ,DOMAIN_INT_STATE                             &
                          ,CLOCK_DOMAIN                                 &
                          ,MYPE                                         &
                          ,CWRT)
        ENDIF
!     call print_memory()
!
      ENDIF 
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE CALL_WRITE_ASYNC
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------

      SUBROUTINE BUILD_FILT_BUNDLE(GRID_DOMAIN                          &
                                  ,UBOUND_VARS                          &
                                  ,VARS                                 &
                                  ,FILT_BUNDLE_FILTER                   &
                                  ,NUM_VARS_2D_FILTER                   &
                                  ,NUM_VARS_3D_FILTER                   &
                                  ,FILT_BUNDLE_RESTORE                  &
                                  ,NUM_VARS_2D_RESTORE                  &
                                  ,NUM_VARS_3D_RESTORE                  & 
                                  ,RESTART)

!-----------------------------------------------------------------------
!***  For digital filtering purposes, the model needs to know both which
!***  variables should be filtered, and which variables need to be saved,
!***  before filtering begins and restored to their original state after
!***  filtering has occurred.  Insert the appropriate internal states into
!***  the appropriate bundles.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      INTEGER(kind=KINT),INTENT(IN) :: UBOUND_VARS                         !<-- Upper dimension of the VARS array
!
      LOGICAL(kind=KLOG),INTENT(IN):: RESTART
!
      TYPE(ESMF_Grid),INTENT(IN) :: GRID_DOMAIN                            !<-- The ESMF Grid for this domain
!
      TYPE(VAR),DIMENSION(1:UBOUND_VARS),INTENT(INOUT) :: VARS             !<-- Variables in the internal state
!
      TYPE(ESMF_FieldBundle),INTENT(INOUT) :: FILT_BUNDLE_FILTER        &  !<-- The Filter Bundle for variables to be filtered
                                             ,FILT_BUNDLE_RESTORE          !<-- The Filter Bundle for variables to be restored
!
!
      INTEGER(kind=KINT),INTENT(INOUT) :: NUM_VARS_2D_FILTER            &  !<-- # of 2-D variables to filter
                                         ,NUM_VARS_3D_FILTER            &  !<-- # of 3-D variables to filter
                                         ,NUM_VARS_2D_RESTORE           &  !<-- # of 2-D variables to restore
                                         ,NUM_VARS_3D_RESTORE              !<-- # of 3-D variables to restore

!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: IOS,N,RC,RC_CMB,FILT_TYP_INT
!
      CHARACTER(len=1) :: UPDATE_TYPE_CHAR
!
      CHARACTER(len=2) :: CH_FILTREST
!
      CHARACTER(len=25),SAVE :: FNAME_FILT='filt_vars.txt'   &
                               ,VBL_NAME
!
      CHARACTER(len=256) :: STRING
!
      TYPE(ESMF_Field) :: FIELD_X
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
        NUM_VARS_2D_FILTER=0
        NUM_VARS_3D_FILTER=0
        NUM_VARS_2D_RESTORE=0
        NUM_VARS_3D_RESTORE=0
!
        OPEN(unit=10,file=FNAME_FILT,status='OLD',action='READ'          &  !<-- Open the filtering text file with user specifications
            ,iostat=IOS)
!
        IF(IOS/=0)THEN
          WRITE(0,*)' Failed to open ',FNAME_FILT,' so ABORT!'
          CALL ESMF_FINALIZE(endflag=ESMF_END_ABORT                 &
                            ,rc             =RC)
        ENDIF

!-----------------------------------------------------------------------
      bundle_loop: DO
!-----------------------------------------------------------------------
!
        READ(UNIT=10,FMT="(A)",iostat=IOS)STRING                           !<-- Read in the next specification line
        IF(IOS/=0)EXIT                                                     !<-- Finished reading the specification lines
!
        IF(STRING(1:1)=='#'.OR.TRIM(STRING)=='')THEN
          CYCLE                                                            !<-- Read past comments and blanks.
        ENDIF
!
!-----------------------------------------------------------------------
!***  Read the text line containing the filtering requirements for
!***  variable N then find that variables' place within the VARS
!***  object.
!-----------------------------------------------------------------------
!
        READ(UNIT=STRING,FMT=*,iostat=IOS)VBL_NAME ,CH_FILTREST
!
        CALL FIND_VAR_INDX(VBL_NAME,VARS,UBOUND_VARS,N)
!

        IF (CH_FILTREST(1:1)=='F') THEN
          FILT_TYP_INT=1
        ELSEIF (CH_FILTREST(1:1)=='R') THEN
          FILT_TYP_INT=2
        ELSEIF (CH_FILTREST(1:1)=='B') THEN
          IF (RESTART) THEN
            FILT_TYP_INT=2
          ELSE
            FILT_TYP_INT=1
          ENDIF
        ELSE 
          FILT_TYP_INT=-999
        ENDIF

        build_bundle: IF (FILT_TYP_INT==1) THEN 

        IF (ASSOCIATED(VARS(N)%R2D))THEN                              !<-- 2-D real array 
!
            NUM_VARS_2D_FILTER=NUM_VARS_2D_FILTER+1

            FIELD_X=ESMF_FieldCreate(grid         =GRID_DOMAIN          &  !<-- The ESMF Grid for this domain
                                    ,farray       =VARS(N)%R2D          &  !<-- Nth variable in the VARS array
                                    ,totalUWidth=(/IHALO,JHALO/)        &  !<-- Upper bound of halo region
                                    ,totalLWidth=(/IHALO,JHALO/)        &  !<-- Lower bound of halo region
                                    ,name         =VARS(N)%VBL_NAME     &  !<-- The name of this variable
                                    ,indexFlag    =ESMF_INDEX_GLOBAL    &  !<-- The variable uses global indexing
                                    ,rc           =RC)
!
        ELSEIF (ASSOCIATED(VARS(N)%R3D))THEN
!
           NUM_VARS_3D_FILTER=NUM_VARS_3D_FILTER+1

           FIELD_X=ESMF_FieldCreate(grid            =GRID_DOMAIN                    &  !<-- The ESMF Grid for this domain
                                    ,farray         =VARS(N)%R3D                    &  !<-- Nth variable in the VARS array
                                    ,totalUWidth    =(/IHALO,JHALO/)                &  !<-- Upper bound of halo region
                                    ,totalLWidth    =(/IHALO,JHALO/)                &  !<-- Lower bound of halo region
                                    ,ungriddedLBound=(/lbound(VARS(N)%R3D,dim=3)/)  &
                                    ,ungriddedUBound=(/ubound(VARS(N)%R3D,dim=3)/)  &
                                    ,name           =VARS(N)%VBL_NAME               &  !<-- The name of this variable
                                    ,indexFlag      =ESMF_INDEX_GLOBAL              &  !<-- The variable uses global indexing
                                    ,rc             =RC)
!
!------------------------
!***  No others expected
!------------------------
!
          ELSE
            WRITE(0,*)' SELECTED FILTER VARIABLE IS NOT 2-D or 3-D REAL'
            WRITE(0,*)' Variable name is ',VARS(N)%VBL_NAME,' for variable #',N
            WRITE(0,*)' ABORT!!'
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
          ENDIF
!
!-----------------------------------------------------------------------
!***  Add this Field to the Filt Bundle that holds all of the
!***  Fields that must be processed through digital filtering.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Add Field to the Filtering Bundle"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_FieldBundleAdd(FILT_BUNDLE_FILTER            &  !<-- The Filt Bundle for Filtered variables
                                  ,(/FIELD_X/)          &  !<-- Add this Field to the Bundle
                                  ,rc    =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CMB)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!       IF (MYPE == 0) THEN
!         WRITE(0,*)' added variable ',trim(VARS(N)%VBL_NAME),' to Filt Bundle Filter'
!       ENDIF
!-----------------------------------------------------------------------


        ELSEIF (FILT_TYP_INT==2) THEN 

        IF (ASSOCIATED(VARS(N)%R2D))THEN                              !<-- 2-D real array 
!
            NUM_VARS_2D_RESTORE=NUM_VARS_2D_RESTORE+1

            FIELD_X=ESMF_FieldCreate(grid         =GRID_DOMAIN          &  !<-- The ESMF Grid for this domain
                                    ,farray       =VARS(N)%R2D          &  !<-- Nth variable in the VARS array
                                    ,totalUWidth=(/IHALO,JHALO/)        &  !<-- Upper bound of halo region
                                    ,totalLWidth=(/IHALO,JHALO/)        &  !<-- Lower bound of halo region
                                    ,name         =VARS(N)%VBL_NAME     &  !<-- The name of this variable
                                    ,indexFlag    =ESMF_INDEX_GLOBAL    &  !<-- The variable uses global indexing
                                    ,rc           =RC)
!
        ELSEIF (ASSOCIATED(VARS(N)%R3D))THEN
!
           NUM_VARS_3D_RESTORE=NUM_VARS_3D_RESTORE+1

           FIELD_X=ESMF_FieldCreate(grid           =GRID_DOMAIN                     &  !<-- The ESMF Grid for this domain
                                    ,farray         =VARS(N)%R3D                    &  !<-- Nth variable in the VARS array
                                    ,totalUWidth    =(/IHALO,JHALO/)                &  !<-- Upper bound of halo region
                                    ,totalLWidth    =(/IHALO,JHALO/)                &  !<-- Lower bound of halo region
                                    ,ungriddedLBound=(/lbound(VARS(N)%R3D,dim=3)/)  &
                                    ,ungriddedUBound=(/ubound(VARS(N)%R3D,dim=3)/)  &
                                    ,name           =VARS(N)%VBL_NAME               &  !<-- The name of this variable
                                    ,indexFlag      =ESMF_INDEX_GLOBAL              &  !<-- The variable uses global indexing
                                    ,rc             =RC)
!----------------
!***  No others expected
!----------------
!
          ELSE
            WRITE(0,*)' SELECTED FILTER VARIABLE IS NOT 2-D OR 3-D REAL'
            WRITE(0,*)' Variable name is ',VARS(N)%VBL_NAME,' for variable #',N
            WRITE(0,*)' ABORT!!'
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
          ENDIF

!-----------------------------------------------------------------------
!***  Add this Field to the Filt Bundle that holds all of the
!***  Fields that must be processed through digital filtering 
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Add Field to the Filtering Bundle"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_FieldBundleAdd(FILT_BUNDLE_RESTORE           &  !<-- The Filt Bundle for variables to be restored
                                  ,(/FIELD_X/)          &  !<-- Add this Field to the Bundle
                                  ,rc    =RC )
!         IF (MYPE == 0) THEN
!           WRITE(0,*)' added variable ',trim(VARS(N)%VBL_NAME),' to Filt Bundle Restore'
!         ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CMB)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------


        ELSE
!          if (MYPE==0) then
!          write(0,*) 'will ignore ' , VARS(N)%VBL_NAME , ' for Bundle'
!          endif

        ENDIF build_bundle



      ENDDO bundle_loop

      close(unit=10)


      END SUBROUTINE BUILD_FILT_BUNDLE
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE BUILD_MOVE_BUNDLE(GRID_DOMAIN                          &
                                  ,UBOUND_VARS                          &
                                  ,VARS                                 &
                                  ,MOVE_BUNDLE_H                        &
                                  ,NUM_VARS_2D_H_I                      &
                                  ,NUM_VARS_2D_H_R                      &
                                  ,NUM_VARS_3D_H                        &
                                  ,NUM_LEVELS_3D_H                      &
                                  ,MOVE_BUNDLE_V                        &
                                  ,NUM_VARS_2D_V                        &
                                  ,NUM_VARS_3D_V                        &
                                  ,NUM_LEVELS_3D_V                      &
                                    )
!
!-----------------------------------------------------------------------
!***  Following a nest's move its appropriate variables will be
!***  updated.  The Solver internal state variables lie within 
!***  their respective VARS composite arrays.  Insert the update
!***  variables from the internal state into the Bundles.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      INTEGER(kind=KINT),INTENT(IN) :: UBOUND_VARS                         !<-- Upper dimension of the VARS array
!
      TYPE(ESMF_Grid),INTENT(IN) :: GRID_DOMAIN                            !<-- The ESMF Grid for this domain
!
      TYPE(VAR),DIMENSION(1:UBOUND_VARS),INTENT(INOUT) :: VARS             !<-- Variables in the internal state
!
      TYPE(ESMF_FieldBundle),INTENT(INOUT) :: MOVE_BUNDLE_H             &  !<-- The Move Bundle for H-point update variables
                                             ,MOVE_BUNDLE_V                !<-- The Move Bundle for V-point update variables
!
      INTEGER(kind=KINT),INTENT(INOUT) :: NUM_LEVELS_3D_H               &  !<-- # of 2-D levels for all 3-D H-point variables
                                         ,NUM_LEVELS_3D_V                  !<-- # of 2-D levels for all 3-D V-point variables
!
      INTEGER(kind=KINT),INTENT(INOUT) :: NUM_VARS_2D_H_I               &  !<-- # of 2-D integer H variables updated for moving nests
                                         ,NUM_VARS_2D_H_R               &  !<-- # of 2-D real H variables updated for moving nests
                                         ,NUM_VARS_3D_H                 &  !<-- # of 3-D real H variables updated for moving nests
                                         ,NUM_VARS_2D_V                 &  !<-- # of 2-D V variables updated for moving nests
                                         ,NUM_VARS_3D_V                    !<-- # of 3-D V variables updated for moving nests
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: IOS,N,RC,RC_CMB,UPDATE_TYPE_INT
!
      CHARACTER(len=1) :: CH_B,CHECK_EXCH,UPDATE_TYPE_CHAR
!           
      CHARACTER(len=2) :: CH_M
!           
      CHARACTER(len=9),SAVE :: FNAME='nests.txt'
!
      CHARACTER(len=99) :: FIELD_NAME,VBL_NAME
!
      CHARACTER(len=256) :: STRING
!
      TYPE(ESMF_Field) :: FIELD_X
!
      LOGICAL(kind=KLOG) :: EXCH_NEEDED
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Loop through all internal state variables.
!-----------------------------------------------------------------------
!
      OPEN(unit=10,file=FNAME,status='OLD',action='READ'                &  !<-- Open the text file with user specifications
            ,iostat=IOS)
!
      IF(IOS/=0)THEN
        WRITE(0,*)' Failed to open ',FNAME,' so ABORT!'
        CALL ESMF_FINALIZE(endflag=ESMF_END_ABORT                   &
                          ,rc             =RC)
      ENDIF
!
!-----------------------------------------------------------------------
      bundle_loop: DO
!-----------------------------------------------------------------------
!
        READ(UNIT=10,FMT="(A)",iostat=IOS)STRING                           !<-- Read in the next specification line
        IF(IOS/=0)EXIT                                                     !<-- Finished reading the specification lines
!
        IF(STRING(1:1)=='#'.OR.TRIM(STRING)=='')THEN
          CYCLE                                                            !<-- Read past comments and blanks.
        ENDIF
!
!-----------------------------------------------------------------------
!***  Read the text line containing the shift specifications for 
!***  variable N then find that variables' place within the VARS 
!***  object.
!-----------------------------------------------------------------------
!
        READ(UNIT=STRING,FMT=*,iostat=IOS)VBL_NAME                      &  !<-- The variable's name in the text file.
                                         ,CH_B                          &  !<-- Not relevant here (specification for nest BC vbls)
                                         ,CH_M                             !<-- Specification for moving nests
!
!
        CALL FIND_VAR_INDX(VBL_NAME,VARS,UBOUND_VARS,N)
!
        FIELD_NAME=TRIM(VARS(N)%VBL_NAME)//SUFFIX_MOVE
!
!-----------------------------------------------------------------------
!***  Find the 2-D and 3-D arrays in the internal state that need
!***  updating in moving nests and add them to the Move Bundle.
!***  We will also specify whether the Field is a simple H-pt variable,
!***  an H-pt land surface variable that needs to use the sea mask,
!***  a variable that is read in from an external file, or a simple
!***  V-pt variable.
!                                NOTE
!***  Currently ESMF will not allow the use of Attributes that are
!***  characters therefore we must translate the character codes from
!***  the txt files into something that ESMF can use.  In this case
!***  we will use integers.
!-----------------------------------------------------------------------
!
        UPDATE_TYPE_CHAR=CH_M(1:1)                                         !<-- Specification flag for this Field
!
        IF(UPDATE_TYPE_CHAR=='H')THEN
          UPDATE_TYPE_INT=1                                                !<-- Ordinary H-pt variable
        ELSEIF(UPDATE_TYPE_CHAR=='L')THEN
          UPDATE_TYPE_INT=2                                                !<-- H-pt land surface variable
        ELSEIF(UPDATE_TYPE_CHAR=='W')THEN
          UPDATE_TYPE_INT=3                                                !<-- H-pt water surface variable
        ELSEIF(UPDATE_TYPE_CHAR=='F')THEN
          UPDATE_TYPE_INT=4                                                !<-- H-pt variable read from external file
        ELSEIF(UPDATE_TYPE_CHAR=='V')THEN
          UPDATE_TYPE_INT=5                                                !<-- Ordinary V-pt variable
        ELSE
          UPDATE_TYPE_INT=-999                                             !<-- Variable not specified for moving nest shifts
        ENDIF
!
!-----------------------------------------------------------------------
!***  Does the variable need to have its halos exchanged so parents
!***  are able to properly update nest points?
!-----------------------------------------------------------------------
!
        CHECK_EXCH=CH_M(2:2)
        IF(CHECK_EXCH=='x')THEN
          EXCH_NEEDED=.TRUE.
        ELSE
          EXCH_NEEDED=.FALSE.
        ENDIF
!
!-----------------------------------------------------------------------
!
        build_bundle: IF(UPDATE_TYPE_CHAR=='H'                          &
                                 .OR.                                   &
                         UPDATE_TYPE_CHAR=='L'                          &
                                 .OR.                                   &
                         UPDATE_TYPE_CHAR=='W'                          &
                                 .OR.                                   &
                         UPDATE_TYPE_CHAR=='F'                          &
                                              )THEN
!
!-----------------------------------------------------------------------
!
!---------------------
!***  2-D H Variables
!---------------------
!
!-------------
!***  Integer
!-------------
!
          IF(ASSOCIATED(VARS(N)%I2D))THEN                                  !<-- 2-D integer array on mass points
!
            NUM_VARS_2D_H_I=NUM_VARS_2D_H_I+1                              !<-- ALL 2-D integer variables updated on H points
!
            FIELD_X=ESMF_FieldCreate(grid       =GRID_DOMAIN            &  !<-- The ESMF Grid for this domain
                                    ,farray     =VARS(N)%I2D            &  !<-- Nth variable in the VARS array
                                    ,totalUWidth=(/IHALO,JHALO/)        &  !<-- Upper bound of halo region
                                    ,totalLWidth=(/IHALO,JHALO/)        &  !<-- Lower bound of halo region
                                    ,name       =FIELD_NAME             &  !<-- The name of this variable
                                    ,indexFlag  =ESMF_INDEX_GLOBAL      &  !<-- The variable uses global indexing
                                    ,rc         =RC)
!
!----------
!***  Real
!----------
!
          ELSEIF(ASSOCIATED(VARS(N)%R2D))THEN                              !<-- 2-D real array on mass points
!
            NUM_VARS_2D_H_R=NUM_VARS_2D_H_R+1                              !<-- ALL 2-D real variables updated on H points
!
            FIELD_X=ESMF_FieldCreate(grid       =GRID_DOMAIN            &  !<-- The ESMF Grid for this domain
                                    ,farray     =VARS(N)%R2D            &  !<-- Nth variable in the VARS array
                                    ,totalUWidth=(/IHALO,JHALO/)        &  !<-- Upper bound of halo region
                                    ,totalLWidth=(/IHALO,JHALO/)        &  !<-- Lower bound of halo region
                                    ,name       =FIELD_NAME             &  !<-- The name of this variable
                                    ,indexFlag  =ESMF_INDEX_GLOBAL      &  !<-- The variable uses global indexing
                                    ,rc         =RC)
!
!---------------------
!***  3-D H Variables
!---------------------
!
!----------
!***  Real
!----------
!
          ELSEIF(ASSOCIATED(VARS(N)%R3D))THEN                              !<-- 3-D real array on mass points
!
            NUM_VARS_3D_H=NUM_VARS_3D_H+1
!
            FIELD_X=ESMF_FieldCreate(grid           =GRID_DOMAIN                    &  !<-- The ESMF Grid for this domain
                                    ,farray         =VARS(N)%R3D                    &  !<-- Nth variable in the VARS array
                                    ,totalUWidth    =(/IHALO,JHALO/)                &  !<-- Upper bound of halo region
                                    ,totalLWidth    =(/IHALO,JHALO/)                &  !<-- Lower bound of halo region
                                    ,ungriddedLBound=(/lbound(VARS(N)%R3D,dim=3)/)  &
                                    ,ungriddedUBound=(/ubound(VARS(N)%R3D,dim=3)/)  &
                                    ,name           =FIELD_NAME                     &  !<-- The name of this variable
                                    ,indexFlag      =ESMF_INDEX_GLOBAL              &  !<-- The variable uses global indexing
                                    ,rc             =RC)
!
            NUM_LEVELS_3D_H=(UBOUND(VARS(N)%R3D,3)-LBOUND(VARS(N)%R3D,3)+1)         &
                            +NUM_LEVELS_3D_H
!
!----------------
!***  All Others
!----------------
!
          ELSE
            WRITE(0,*)' SELECTED UPDATE H VARIABLE IS NOT 2-D OR 3-D.'
            WRITE(0,*)' Variable name is ',VARS(N)%VBL_NAME,' for variable #',N
            WRITE(0,*)' UPDATE_TYPE=',UPDATE_TYPE_CHAR
            WRITE(0,*)' ABORT!!'
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
!
          ENDIF
!
!-----------------------------------------------------------------------
!***  Attach the specification flag to this Field that indicates
!***  how it must be handled in the parent-child update region.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Add Specification Flag to Move Bundle H Field"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeSet(field=FIELD_X                          &  !<-- The Field to be added to H-pt Move Bundle
                                ,name ='UPDATE_TYPE'                    &  !<-- The name of the Attribute to set
                                ,value=UPDATE_TYPE_INT                  &  !<-- The Attribute to be set
                                ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CMB)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Attach the halo exchange flag to this Field that indicates
!***  to the parent if it must perform exchanges prior to updating
!***  its moving nests.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Add Halo Exchange Flag to Move Bundle H Field"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeSet(field=FIELD_X                          &  !<-- The Field to be added to H-pt Move Bundle
                                ,name ='EXCH_NEEDED'                    &  !<-- The name of the Attribute to set
                                ,value=EXCH_NEEDED                      &  !<-- The Attribute to be set
                                ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CMB)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Add this Field to the Move Bundle that holds all the H-point
!***  Fields that must be shifted after a nest moves. 
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Add Field to the H-pt Move Bundle"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_FieldBundleAdd(            MOVE_BUNDLE_H            &  !<-- The Move Bundle for H point variables
                                  ,            (/FIELD_X/)     &  !<-- Add this Field to the Bundle
                                  ,rc         =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CMB)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
        ELSEIF(UPDATE_TYPE_CHAR=='V')THEN                                  !<-- If so, V variable is designated for moving nest updates
!
!---------------------
!***  2-D V Variables
!---------------------
!
!----------
!***  Real
!----------
!
          IF(ASSOCIATED(VARS(N)%R2D))THEN                                  !<-- 2-D reall array on velocity points
!
            NUM_VARS_2D_V=NUM_VARS_2D_V+1
!
            FIELD_X=ESMF_FieldCreate(grid       =GRID_DOMAIN            &  !<-- The ESMF Grid for this domain
                                    ,farray     =VARS(N)%R2D            &  !<-- Nth variable in the VARS array
                                    ,totalUWidth=(/IHALO,JHALO/)        &  !<-- Upper bound of halo region
                                    ,totalLWidth=(/IHALO,JHALO/)        &  !<-- Lower bound of halo region
                                    ,name       =FIELD_NAME             &  !<-- The name of this variable
                                    ,indexFlag  =ESMF_INDEX_GLOBAL      &  !<-- The variable uses global indexing
                                    ,rc         =RC)
!
!---------------------
!***  3-D V Variables
!---------------------
!
!----------
!***  Real
!----------
!
          ELSEIF(ASSOCIATED(VARS(N)%R3D))THEN                              !<-- 3-D real array on velocity points
!
            NUM_VARS_3D_V=NUM_VARS_3D_V+1
!
            FIELD_X=ESMF_FieldCreate(grid           =GRID_DOMAIN                    &  !<-- The ESMF Grid for this domain
                                    ,farray         =VARS(N)%R3D                    &  !<-- Nth variable in the VARS array
                                    ,totalUWidth    =(/IHALO,JHALO/)                &  !<-- Upper bound of halo region
                                    ,totalLWidth    =(/IHALO,JHALO/)                &  !<-- Lower bound of halo region
                                    ,ungriddedLBound=(/lbound(VARS(N)%R3D,dim=3)/)  &
                                    ,ungriddedUBound=(/ubound(VARS(N)%R3D,dim=3)/)  &
                                    ,name           =FIELD_NAME                     &  !<-- The name of this variable
                                    ,indexFlag      =ESMF_INDEX_GLOBAL              &  !<-- The variable uses global indexing
                                    ,rc             =RC)
!
            NUM_LEVELS_3D_V=(UBOUND(VARS(N)%R3D,3)-LBOUND(VARS(N)%R3D,3)+1)         &
                            +NUM_LEVELS_3D_V
!
!
!------------
!***  Others
!------------
!
          ELSE
            WRITE(0,*)' SELECTED UPDATE V VARIABLE IS NOT 2-D OR 3-D.  ABORT!!'
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
!
          ENDIF
!
!-----------------------------------------------------------------------
!***  Attach the specification flag to this Field that indicates
!***  how it must be handled in the parent-child update region.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Add Specification Flag to Move Bundle V Field"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeSet(field=FIELD_X                          &  !<-- The Field to be added to V-pt Move Bundle
                                ,name ='UPDATE_TYPE'                    &  !<-- The name of the Attribute to set
                                ,value=UPDATE_TYPE_INT                  &  !<-- The Attribute to be set
                                ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CMB)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Attach the halo exchange flag to this Field that indicates
!***  to the parent if it must perform exchanges prior to updating
!***  its moving nests.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Add Specification Flag to Move Bundle V Field"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeSet(field=FIELD_X                          &  !<-- The Field to be added to V-pt Move Bundle
                                ,name ='EXCH_NEEDED'                    &  !<-- The name of the Attribute to set
                                ,value=EXCH_NEEDED                      &  !<-- The Attribute to be set
                                ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CMB)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Add this Field to the Move Bundle that holds all the H-point
!***  Fields that must be shifted after a nest moves. 
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Add Field to the H-pt Move Bundle"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_FieldBundleAdd(            MOVE_BUNDLE_V            &  !<-- The Move Bundle for V-point variables
                                  ,            (/FIELD_X/)     &  !<-- Add this Field to the Bundle
                                  ,rc         =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CMB)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!
!-----------------------------------------------------------------------
!
        ENDIF build_bundle
!
!-----------------------------------------------------------------------
!
      ENDDO bundle_loop
!
!-----------------------------------------------------------------------
!
      CLOSE(10)
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE BUILD_MOVE_BUNDLE
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE BUILD_2WAY_BUNDLE(GRID_DOMAIN                          &
                                  ,LM                                   &
                                  ,UBOUND_VARS                          &
                                  ,VARS                                 &
                                  ,BUNDLE_2WAY                          &
                                    )
!
!-----------------------------------------------------------------------
!***  When 2-way exchange is invoked in the configure file then
!***  a child domain will interpolate specified variables from
!***  the Solver component's internal state on its grid to its
!***  parent's grid and send that data to its parent.  Parents
!***  receive that data and incorporate it.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      INTEGER(kind=KINT),INTENT(IN) :: LM                               &  !<-- # of model layers
                                      ,UBOUND_VARS                         !<-- Upper dimension of the VARS array
!
      TYPE(ESMF_Grid),INTENT(IN) :: GRID_DOMAIN                            !<-- The ESMF Grid for this domain
!
      TYPE(VAR),DIMENSION(1:UBOUND_VARS),INTENT(INOUT) :: VARS             !<-- Variables in the Solver internal state
!
      TYPE(ESMF_FieldBundle),INTENT(INOUT) :: BUNDLE_2WAY                  !<-- The Bundle of Solver internal state vbls for 2-way exchange
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: H_OR_V_INT,IOS,N,NLEV,RC,RC_CMB
!
      CHARACTER(len=1) :: CH_2,CH_B,CHECK_EXCH,H_OR_V
!
      CHARACTER(len=2) :: CH_M
!
      CHARACTER(len=9),SAVE :: FNAME='nests.txt'
!
      CHARACTER(len=99) :: FIELD_NAME,VBL_NAME
!
      CHARACTER(len=256) :: STRING
!
      TYPE(ESMF_Field) :: FIELD_X
!
      integer(kind=kint) :: lbnd1,lbnd2,lbnd3,ubnd1,ubnd2,ubnd3,nx,ny,nz
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Loop through all Solver internal state variables.
!-----------------------------------------------------------------------
!
      OPEN(unit=10,file=FNAME,status='OLD',action='READ'                &  !<-- Open the text file with user specifications
            ,iostat=IOS)
!
      IF(IOS/=0)THEN
        WRITE(0,*)' Failed to open ',FNAME,' so ABORT!'
        CALL ESMF_FINALIZE(endflag=ESMF_END_ABORT                   &
                          ,rc             =RC)
      ENDIF
!
      NLEV=0                                                               !<-- Counter for total # of levels in all 2-way vbls
!
!-----------------------------------------------------------------------
      bundle_loop: DO
!-----------------------------------------------------------------------
!
        READ(UNIT=10,FMT='(A)',iostat=IOS)STRING                           !<-- Read in the next specification line
        IF(IOS/=0)EXIT                                                     !<-- Finished reading the specification lines
!
        IF(STRING(1:1)=='#'.OR.TRIM(STRING)=='')THEN
          CYCLE                                                            !<-- Read past comments and blanks.
        ENDIF
!
!-----------------------------------------------------------------------
!***  Read the text line containing the nest specifications for
!***  variable N then find that variables' place within the VARS
!***  object.
!-----------------------------------------------------------------------
!
        READ(UNIT=STRING,FMT=*,iostat=IOS)VBL_NAME                      &  !<-- The variable's name in the text file.
                                         ,CH_B                          &  !<-- Not relevant here (specification for BC vbls)
                                         ,CH_M                          &  !<-- Not relevant here (specification for motion vbls)
                                         ,CH_2                             !<-- Specification for 2-way variables
!
        CALL FIND_VAR_INDX(VBL_NAME,VARS,UBOUND_VARS,N)
!
        FIELD_NAME=TRIM(VARS(N)%VBL_NAME)//SUFFIX_TWOWAY
!
!-----------------------------------------------------------------------
!***  Find the variables in the Solver internal state that will be
!***  used for 2-way exchange and add them to the 2-way Bundle.
!***  We will also specify whether the Field's variable lies on
!***  H points or V points.
!                                NOTE
!***  Currently ESMF will not allow the use of Attributes that are
!***  characters therefore we must translate the character codes from
!***  the txt files into something that ESMF can use.  In this case
!***  we will use integers:  H-->1 and V-->2 .
!-----------------------------------------------------------------------
!
        H_OR_V=CH_2                                                        !<-- H-V flag for this Field
!
        IF(H_OR_V=='H')THEN
          H_OR_V_INT=1                                                     !<-- H-pt variable
        ELSEIF(H_OR_V=='V')THEN
          H_OR_V_INT=2                                                     !<-- V-pt variable
        ELSE
          H_OR_V_INT=-999                                                  !<-- Variable not specified for 2-way exchange.
        ENDIF
!
!-----------------------------------------------------------------------
!
        build_bundle: IF(H_OR_V=='H'                                    &
                            .OR.                                        &
                         H_OR_V=='V'                                    &
                                     )THEN
!
!-----------------------------------------------------------------------
!
!-------------------
!***  2-D Variables
!-------------------
!
!-------------
!***  Integer
!-------------
!
          IF(ASSOCIATED(VARS(N)%I2D))THEN                                  !<-- 2-D integer array on mass points
!
!           FIELD_X=ESMF_FieldCreate(grid       =GRID_DOMAIN            &  !<-- The ESMF Grid for this domain
!                                   ,farray     =VARS(N)%I2D            &  !<-- Nth variable in the VARS array
!                                   ,totalUWidth=(/IHALO,JHALO/)        &  !<-- Upper bound of halo region
!                                   ,totalLWidth=(/IHALO,JHALO/)        &  !<-- Lower bound of halo region
!                                   ,name       =FIELD_NAME             &  !<-- The name of this variable
!                                   ,indexFlag  =ESMF_INDEX_GLOBAL      &  !<-- The variable uses global indexing
!                                   ,rc         =RC)
!
            WRITE(0,*)' The scheme does not consider integer variables.'
            WRITE(0,*)' Variable name is ',VARS(N)%VBL_NAME,' for variable #',N
            WRITE(0,*)' H_OR_V_INT=',H_OR_V_INT
            WRITE(0,*)' ABORT!!'
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!----------
!***  Real
!----------
!
          ELSEIF(ASSOCIATED(VARS(N)%R2D))THEN                              !<-- 2-D real array on mass points
!
            FIELD_X=ESMF_FieldCreate(grid       =GRID_DOMAIN            &  !<-- The ESMF Grid for this domain
                                    ,farray     =VARS(N)%R2D            &  !<-- Nth variable in the VARS array
                                    ,totalUWidth=(/IHALO,JHALO/)        &  !<-- Upper bound of halo region
                                    ,totalLWidth=(/IHALO,JHALO/)        &  !<-- Lower bound of halo region
                                    ,name       =FIELD_NAME             &  !<-- The name of this variable
                                    ,indexFlag  =ESMF_INDEX_GLOBAL      &  !<-- The variable uses global indexing
                                    ,rc         =RC)
!
            NLEV=NLEV+1                                                    !<-- Sum the levels for all real 2-way variables.
!
!---------------------
!***  3-D H Variables
!---------------------
!
!----------
!***  Real
!----------
!
          ELSEIF(ASSOCIATED(VARS(N)%R3D))THEN                              !<-- 3-D real array on mass points
!
            FIELD_X=ESMF_FieldCreate(grid           =GRID_DOMAIN                    &  !<-- The ESMF Grid for this domain
                                    ,farray         =VARS(N)%R3D                    &  !<-- Nth variable in the VARS array
                                    ,totalUWidth    =(/IHALO,JHALO/)                &  !<-- Upper bound of halo region
                                    ,totalLWidth    =(/IHALO,JHALO/)                &  !<-- Lower bound of halo region
                                    ,ungriddedLBound=(/lbound(VARS(N)%R3D,dim=3)/)  &
                                    ,ungriddedUBound=(/ubound(VARS(N)%R3D,dim=3)/)  &
                                    ,name           =FIELD_NAME                     &  !<-- The name of this variable
                                    ,indexFlag      =ESMF_INDEX_GLOBAL              &  !<-- The variable uses global indexing
                                    ,rc             =RC)
!
            NLEV=NLEV+LM                                                   !<-- Sum the levels for all real 2-way variables.
!
!----------------
!***  All Others
!----------------
!
          ELSE
            WRITE(0,*)' Selected update H variable is NOT 2-D OR 3-D Real.'
            WRITE(0,*)' Variable name is ',VARS(N)%VBL_NAME,' for variable #',N
            WRITE(0,*)' H_OR_V_INT=',H_OR_V_INT
            WRITE(0,*)' ABORT!!'
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
!
          ENDIF
!
!-----------------------------------------------------------------------
!***  Attach the specification flag to this Field that indicates
!***  how it must be handled in the parent-child update region.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Add Specification Flag to Move Bundle H Field"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeSet(field=FIELD_X                          &  !<-- The Field to be added to H-pt Move Bundle
                                ,name ='H_OR_V_INT'                     &  !<-- The name of the Attribute to set
                                ,value=H_OR_V_INT                       &  !<-- The Attribute to be set
                                ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CMB)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Add this Field to the 2-way Bundle that holds pointers to all
!***  variables in the Solver internal state that are used in 2-way
!***  exchange.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Add Field to the H-pt Move Bundle"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_FieldBundleAdd(            BUNDLE_2WAY              &  !<-- The Move Bundle for H point variables
                                  ,            (/FIELD_X/)     &  !<-- Add this Field to the Bundle
                                  ,rc         =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CMB)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
        ENDIF build_bundle
!
!-----------------------------------------------------------------------
!
      ENDDO bundle_loop
!
!-----------------------------------------------------------------------
!***  Attach the total number of levels in the 2-way variables,
!***  i.e., one level for each 2-D variable and LM levels for
!***  each 3-D variable.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Add total # of levels for all 2-way variables"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(FIELDBUNDLE=BUNDLE_2WAY                    &  !<-- The Bundle of 2-way variable pointers
                            ,name       ='NLEV 2-way'                   &  !<-- The name of the Attribute to set
                            ,value      =NLEV                           &  !<-- The # of 2-way Real BC variables
                            ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CMB)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CLOSE(10)
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE BUILD_2WAY_BUNDLE
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE UPDATE_GRID_ARRAYS(DOMAIN_IMP_STATE                    &
                                   ,SOLVER_GRID_COMP)
!
!-----------------------------------------------------------------------
!***  When a nest moves we must update the 1-D (in J) grid-dependent
!***  arrays which span the entire nest north-south dimension.
!***  In addition update the value of the nest's SW corner on its
!***  parent's grid.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(ESMF_State),INTENT(INOUT) :: DOMAIN_IMP_STATE                   !<-- The Domain component's import state
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: SOLVER_GRID_COMP                !<-- The Solver Component
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: I_SW_NEW,J,J_SW_NEW
!
      INTEGER(kind=KINT) :: RC,RC_FINAL
!
      REAL(kind=KFPT) :: ARG1,ARG2,DY,TLAT_SW,TLON_SW
!
      REAL(kind=KFPT),DIMENSION(JDS:JDE) :: TLAT_H,TLAT_V
!
      TYPE(WRAP_SOLVER_INT_STATE) :: WRAP_SOLVER
!
      TYPE(SOLVER_INTERNAL_STATE),POINTER :: SOLVER_INT_STATE
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Extract this domain's ID and the shifts in I and J that the
!***  nest is executing.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Domain ID from the Domain Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=DOMAIN_IMP_STATE                     &  !<-- The Domain import state
                            ,name ='DOMAIN_ID'                          &  !<-- Get Attribute with this name
                            ,value=MY_DOMAIN_ID                         &  !<-- This domain's ID
                            ,rc   =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Get I_SHIFT and J_SHIFT from Domain Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=DOMAIN_IMP_STATE                     &  !<-- The Domain import state
                            ,name ='I_SHIFT'                            &  !<-- Get Attribute with this name
                            ,value=I_SHIFT_CHILD                        &  !<-- Motion of the nest in I on its grid
                            ,rc   =RC )
!
      CALL ESMF_AttributeGet(state=DOMAIN_IMP_STATE                     &  !<-- The Domain import state
                            ,name ='J_SHIFT'                            &  !<-- Get Attribute with this name
                            ,value=J_SHIFT_CHILD                        &  !<-- Motion of the nest in J on its grid
                             ,rc   =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Extract the Solver internal state so we can access its contents.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Solver Internal State for Move Bundle"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompGetInternalState(SOLVER_GRID_COMP               &  !<-- The Solver component
                                        ,WRAP_SOLVER                    &
                                        ,RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
      SOLVER_INT_STATE=>wrap_solver%INT_STATE
!
!-----------------------------------------------------------------------
!***  What are the new coordinates on the parent's grid of the nest's
!***  SW corner after the shift?
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Get Parent-Child Space Ratio from Configure File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)              &  !<-- The config object
                                  ,value =PARENT_CHILD_SPACE_RATIO      &  !<-- The variable filled (child grid increment / parent's)
                                  ,label ='parent_child_space_ratio:'   &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      solver_int_state%I_PAR_STA=solver_int_state%I_PAR_STA             &
                                +I_SHIFT_CHILD/PARENT_CHILD_SPACE_RATIO
!
      solver_int_state%J_PAR_STA=solver_int_state%J_PAR_STA             &
                                +J_SHIFT_CHILD/PARENT_CHILD_SPACE_RATIO
!
!-----------------------------------------------------------------------
!***  The arrays are tied to the nest grid's transformed latitude.
!***  After determining the transformed latitude of the subdomain's
!***  SW corner following the move the rest can be filled in.
!-----------------------------------------------------------------------
!
      CALL GEO_TO_ROT(solver_int_state%GLAT_SW                          &  !<-- The pre-move geographic lat of nest's SW corner (radians)
                     ,solver_int_state%GLON_SW                          &  !<-- The pre-move geographic lon of nest's SW corner (radians)
                     ,TLAT_SW                                           &  !<-- The pre-move rotated lat of nest's SW corner (radians)
                     ,TLON_SW )                                            !<-- The pre-move rotated lon of nest's SW corner (radians)
!
      DPH=solver_int_state%DPHD*DEG_TO_RAD
      DLM=solver_int_state%DLMD*DEG_TO_RAD
!
      TLAT_H(JDS)=TLAT_SW+J_SHIFT_CHILD*DPH
      TLAT_V(JDS)=TLAT_H(JDS)+0.5*DPH
!
      DO J=JDS+1,JDE
        TLAT_H(J)=TLAT_H(JDS)+(J-JDS)*DPH
        TLAT_V(J)=TLAT_V(JDS)+(J-JDS)*DPH
      ENDDO
!
      DY=A*DPH
!
!-----------------------------------------------------------------------
!***  Update the relevant 1-D arrays.
!-----------------------------------------------------------------------
!
      DO J=JDS,JDE
        solver_int_state%DXH(J)=A*DLM*COS(TLAT_H(J))
        solver_int_state%RDXH(J)=1./solver_int_state%DXH(J)
        solver_int_state%DXV(J)=A*DLM*COS(TLAT_V(J))
        solver_int_state%RDXV(J)=1./solver_int_state%DXV(J)
        solver_int_state%DARE(J)=solver_int_state%DXH(J)*DY
        solver_int_state%RARE(J)=1./solver_int_state%DARE(J)
        solver_int_state%WPDAR(J)=-1.E-5*WCOR*DY*DY                        &
                               /(DT_REAL*solver_int_state%DXH(J)*DY)
        solver_int_state%CURV(J)=TAN(TLAT_V(J))/A
        solver_int_state%FAH(J)=-DT_REAL/(3.*solver_int_state%DXH(J)*DY)
        solver_int_state%FAD(J)=-0.25*DT_REAL/(3.*solver_int_state%DXV(J)*DY)
        solver_int_state%FCP(J)=DT_REAL/(3.*solver_int_state%DXH(J)*DY*CP)
        solver_int_state%FDIV(J)=2./(3.*solver_int_state%DXH(J)*DY)
        solver_int_state%DDV(J)=SQRT(solver_int_state%DXV(J)**2+DY*DY)
        solver_int_state%RDDV(J)=1./solver_int_state%DDV(J)
        solver_int_state%DDMPU(J)=0.5*CDDAMP*DY/solver_int_state%DXV(J)
      ENDDO
!
!-----------------------------------------------------------------------
!***  Compute the new geographic coordinates of the nest's SW corner
!***  after it shifts.
!-----------------------------------------------------------------------
!
      TPH0=solver_int_state%TPH0D*DEG_TO_RAD
      TLM0=solver_int_state%TLM0D*DEG_TO_RAD
!
      TLAT_SW=TLAT_H(JDS)                                                  !<-- Transformed lat (radians) of SW corner after shift
      TLON_SW=TLON_SW+I_SHIFT_CHILD*DLM                                    !<-- Transformed lon (radians) of SW corner after shift
!
      solver_int_state%GLAT_SW=ASIN(SIN(TLAT_SW)*COS(TPH0)              &  !<-- Geographic lat (radians) of SW corner after shift
                                   +COS(TLAT_SW)*SIN(TPH0)*COS(TLON_SW))
!
      ARG1=(COS(TLAT_SW)*COS(TLON_SW))/(COS(solver_int_state%GLAT_SW)   &
                                       *COS(TPH0))
      ARG2=TAN(solver_int_state%GLAT_SW)*TAN(TPH0)
      solver_int_state%GLON_SW=TLM0+SIGN(1.,TLON_SW)*ACOS(ARG1-ARG2)       !<-- Geographic lon (radians) of SW corner after shift
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE UPDATE_GRID_ARRAYS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE UPDATE_INTERIOR_FROM_NEST(IMP_STATE                    &
                                          ,MOVE_BUNDLE_H                &
                                          ,NUM_FIELDS_2D_H_I            &
                                          ,NUM_FIELDS_2D_H_R            &
                                          ,NUM_FIELDS_3D_H              &
                                          ,NUM_LEVELS_3D_H              &
                                          ,MOVE_BUNDLE_V                &
                                          ,NUM_FIELDS_2D_V              &
                                          ,NUM_FIELDS_3D_V              &
                                          ,NUM_LEVELS_3D_V              &
                                          ,INPES,JNPES                  &
                                          ,ITS,ITE,JTS,JTE              &
                                          ,IMS,IME,JMS,JME              &
                                          ,IDS,IDE,JDS,JDE              &
                                            )
!
!-----------------------------------------------------------------------
!***  After the nest has moved update all nest gridpoints in that
!***  domain's interior that remain within the footprint of its
!***  pre-move location.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      INTEGER(kind=KINT),INTENT(IN) :: NUM_FIELDS_2D_H_I                 &  !<-- # of 2-D integer H variables to update
                                      ,NUM_FIELDS_2D_H_R                 &  !<-- # of 2-D real H variables to update
                                      ,NUM_FIELDS_3D_H                   &  !<-- # of 3-D H variables to update
                                      ,NUM_LEVELS_3D_H                   &  !<-- # of 2-D levels in all 3-D H update variables
                                      ,NUM_FIELDS_2D_V                   &  !<-- # of 3-D V variables to update
                                      ,NUM_FIELDS_3D_V                   &  !<-- # of 3-D V variables to update
                                      ,NUM_LEVELS_3D_V                      !<-- # of 2-D levels in all 3-D V update variables
!
      INTEGER(kind=KINT),INTENT(IN) :: INPES,JNPES                          !<-- # of tasks west-east,north-south on this domain
!
      INTEGER(kind=KINT),INTENT(IN) :: ITS,ITE,JTS,JTE                   &  !<-- This domain's integration,
                                      ,IMS,IME,JMS,JME                   &  !    memory,
                                      ,IDS,IDE,JDS,JDE                      !<-- and domain limits.
!
      TYPE(ESMF_State),INTENT(INOUT) :: IMP_STATE                           !<-- The Domain import state
!
      TYPE(ESMF_FieldBundle),INTENT(INOUT) :: MOVE_BUNDLE_H              &  !<-- Bundle of internal state H arrays needing updates
                                             ,MOVE_BUNDLE_V                 !<-- Bundle of internal state V arrays needing updates
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: I,I_DIFF_MAX,I_END,I_START                   &
                           ,J,J_DIFF_MAX,J_END,J_START                   &
                           ,L,N
!
      INTEGER(kind=KINT) :: RC,RC_UPD
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC    =ESMF_SUCCESS
      RC_UPD=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  Extract the shifts in I and J that the nest is executing.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Get I_SHIFT and J_SHIFT from Domain Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=IMP_STATE                            &  !<-- The Domain import state
                            ,name ='I_SHIFT'                            &  !<-- Get Attribute with this name
                            ,value=I_SHIFT_CHILD                        &  !<-- Motion of the nest in I on its grid
                            ,rc   =RC )
!
      CALL ESMF_AttributeGet(state=IMP_STATE                            &  !<-- The Domain import state
                            ,name ='J_SHIFT'                            &  !<-- Get Attribute with this name
                            ,value=J_SHIFT_CHILD                        &  !<-- Motion of the nest in J on its grid
                             ,rc   =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_UPD)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  After a nested domain changes its position there are two ways
!***  in which its internal gridpoints can be updated that do not
!***  involve the parent.  The points updated are those that remain
!***  within the footprint of the domain's previous position.  The
!***  two types of updates are:
!
!***   (a) The new value of each variable will come from a different
!***       location on the same nest task's subdomain (intra-task).
!***   (b) The new values will be received from a different nest
!***       task (inter-task).
!
!***  These actions cannot be done in a simple sequence because if
!***  they were then data would be lost that was needed for either
!***  the intra- or inter-task shift.  Therefore they are done in
!***  the following order:
!
!***   (1) The data is first gathered into ISend buffers for the 
!***       inter-task shift and then sent.  
!***   (2) The intra-task update is done.
!***   (3) The inter-task data is Recvd and applied.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Gather data into buffers for the inter-task shift within the
!***  nest domain and send it.
!-----------------------------------------------------------------------
!
      CALL SEND_INTER_TASK_DATA(I_SHIFT_CHILD                           &
                               ,J_SHIFT_CHILD                           &
                               ,MYPE                                    &
                               ,INPES                                   &
                               ,JNPES                                   &
                               ,MOVE_BUNDLE_H                           &
                               ,NUM_FIELDS_2D_H_I                       &
                               ,NUM_FIELDS_2D_H_R                       &
                               ,NUM_FIELDS_3D_H                         &
                               ,NUM_LEVELS_3D_H                         &
                               ,MOVE_BUNDLE_V                           &
                               ,NUM_FIELDS_2D_V                         &
                               ,NUM_FIELDS_3D_V                         &
                               ,NUM_LEVELS_3D_V                         &
                               ,ITS,ITE,JTS,JTE                         &
                               ,IMS,IME,JMS,JME                         &
                               ,IDS,IDE,JDS,JDE                         &
                                )
!
!-----------------------------------------------------------------------
!***  Carry out the shift on all points that remains on the same task
!***  after the domain moves.
!-----------------------------------------------------------------------
!
      CALL SHIFT_INTRA_TASK_DATA(I_SHIFT_CHILD                          &
                                ,J_SHIFT_CHILD                          &
                                ,MOVE_BUNDLE_H                          &
                                ,NUM_FIELDS_2D_H_I                      &
                                ,NUM_FIELDS_2D_H_R                      &
                                ,NUM_FIELDS_3D_H                        &
                                ,MOVE_BUNDLE_V                          &
                                ,NUM_FIELDS_2D_V                        &
                                ,NUM_FIELDS_3D_V                        &
                                ,ITS,ITE,JTS,JTE                        &
                                ,IMS,IME,JMS,JME                        &
                                ,IDS,IDE,JDS,JDE                        &
                                 )
!
!-----------------------------------------------------------------------
!***  Receive data for the inter-task shift within the nest domain
!***  and apply it.
!-----------------------------------------------------------------------
!
      CALL RECV_INTER_TASK_DATA(I_SHIFT_CHILD                           &
                               ,J_SHIFT_CHILD                           &
                               ,MYPE                                    &
                               ,INPES                                   &
                               ,JNPES                                   &
                               ,MOVE_BUNDLE_H                           &
                               ,NUM_FIELDS_2D_H_I                       &
                               ,NUM_FIELDS_2D_H_R                       &
                               ,NUM_FIELDS_3D_H                         &
                               ,NUM_LEVELS_3D_H                         &
                               ,MOVE_BUNDLE_V                           &
                               ,NUM_FIELDS_2D_V                         &
                               ,NUM_FIELDS_3D_V                         &
                               ,NUM_LEVELS_3D_V                         &
                               ,ITS,ITE,JTS,JTE                         &
                               ,IMS,IME,JMS,JME                         &
                               ,IDS,IDE,JDS,JDE                         &
                                 )
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE UPDATE_INTERIOR_FROM_NEST
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE SHIFT_INTRA_TASK_DATA(I_SHIFT_CHILD                    &
                                      ,J_SHIFT_CHILD                    &
                                      ,MOVE_BUNDLE_H                    &
                                      ,NUM_FIELDS_2D_H_I                &
                                      ,NUM_FIELDS_2D_H_R                &
                                      ,NUM_FIELDS_3D_H                  &
                                      ,MOVE_BUNDLE_V                    &
                                      ,NUM_FIELDS_2D_V                  &
                                      ,NUM_FIELDS_3D_V                  &
                                      ,ITS,ITE,JTS,JTE                  &
                                      ,IMS,IME,JMS,JME                  &
                                      ,IDS,IDE,JDS,JDE                  &
                                        )
!
!-----------------------------------------------------------------------
!***  After the nest has moved update all nest gridpoints in the
!***  domain's interior whose new locations still lie within the
!***  same MPI task subdomain (same processor memory) as before 
!***  the move.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      INTEGER(kind=KINT),INTENT(IN) :: I_SHIFT_CHILD                    &  !<-- Nest domain moved this far in I in nest space
                                      ,J_SHIFT_CHILD                    &  !<-- Nest domain moved this far in J in nest space
                                      ,NUM_FIELDS_2D_H_I                &  !<-- # of 2-D integer H variables to update
                                      ,NUM_FIELDS_2D_H_R                &  !<-- # of 2-D real H variables to update
                                      ,NUM_FIELDS_3D_H                  &  !<-- # of 3-D H variables to update
                                      ,NUM_FIELDS_2D_V                  &  !<-- # of 2-D V variables to update
                                      ,NUM_FIELDS_3D_V                  &  !<-- # of 3-D V variables to update
!
                                      ,ITS,ITE,JTS,JTE                  &  !<-- Subdomain integration limits of this nest task
                                      ,IMS,IME,JMS,JME                  &  !<-- Subdomain memory limits of this nest task
                                      ,IDS,IDE,JDS,JDE                     !<-- Index limits of this nest's full domain
!
      TYPE(ESMF_FieldBundle),INTENT(INOUT) :: MOVE_BUNDLE_H             &  !<-- Bundle of internal state H arrays for updates
                                             ,MOVE_BUNDLE_V                !<-- Bundle of internal state V arrays for updates
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: I,I_DIFF_MAX,I_END,I_INC                    &
                           ,I_START,ITE_X,ITS_X                         &
                           ,J,J_DIFF_MAX,J_END,J_INC                    &
                           ,J_START,JTE_X,JTS_X                         &
                           ,L,N,N_FIELD,N_REMOVE,NF1,NF2                &
                           ,NUM_DIMS,NUM_FIELDS
!
      INTEGER(kind=KINT) :: RC
!
      INTEGER(kind=KINT),DIMENSION(1:3) :: LIMITS_HI                    &
                                          ,LIMITS_LO
!
      INTEGER(kind=KINT),DIMENSION(:,:),POINTER :: IARRAY_2D
!
      REAL(kind=KFPT),DIMENSION(:,:),POINTER :: ARRAY_2D
!
      REAL(kind=KFPT),DIMENSION(:,:,:),POINTER :: ARRAY_3D
!
      CHARACTER(len=99) :: FIELD_NAME
!
      TYPE(ESMF_Field) :: HOLD_FIELD
!
      TYPE(ESMF_TypeKind_Flag) :: DATATYPE
!
      integer(kind=kint),save :: kount_moves=0
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  It is critical to realize that neither H-pt nor V-pt variables
!***  on the nest domain's north and east limits can be used in
!***  the intra/intertask updates.  That is because the V-pt variables
!***  at IDE and JDE are not part of the integration.  Although the
!***  H-pt variables at the domain's IDE and JDE are part of the
!***  integration if a task is on the nest domain's east or
!***  north boundary we must exclude them in the Send too otherwise
!***  there would be occasions when index limits for shift updates
!***  of H-pt data would not be identical to the index limits for
!***  shift updates of V-pt data and we must avoid that situation
!***  since it would greatly complicate the already very complicated
!***  bookkeeping.  This fact must also be accounted for in the inter-
!***  task shifts and in the sending/recving of parent update data 
!***  after the nest moves.  At the same time there are a variety of
!***  variables that do not have valid integration values on the nest
!***  domain boundary so we do not want to shift those into the 
!***  interior via intra- or inter-task shifts.  Therefore the updates
!***  of the nest's south and west boundary points must also be done
!***  by the parent.  Moreover the dynamical tendencies of temperature
!***  and the wind components are not computed in the 2nd row in from
!***  the domain boundary which means nest points shifted to those
!***  pre-move locations cannot use intra- or inter-task updates.
!***  Therefore those locations must also be updated by the parent.
!***  In general the gridpoint locations on the outer two rows of the
!***  nest's pre-move footprint will be updated by the parent although
!***  the actual depth into the footprint that the parent will provide
!***  data will be specified by configure variables for generality.
!***  Use the following quantities in searches for points involved
!***  in the intra-task update.
!-----------------------------------------------------------------------
!
      ITS_X=MAX(ITS,IDS+NROWS_P_UPD_W)                                     !<-- These quantities indicate the 
      ITE_X=MIN(ITE,IDE-NROWS_P_UPD_E)                                     !    outermost locations in the nest
      JTS_X=MAX(JTS,JDS+NROWS_P_UPD_S)                                     !    domain subject to the intra-task
      JTE_X=MIN(JTE,JDE-NROWS_P_UPD_N)                                     !<-- updates.
!
!-----------------------------------------------------------------------
!***  What is the maximum shift in I and J that can occur for which
!***  there will be points that require an intra-task shift?
!-----------------------------------------------------------------------
!
      I_DIFF_MAX=ITE_X-ITS_X+IHALO                                         !<-- Maximum I shift for intra-task update
      J_DIFF_MAX=JTE_X-JTS_X+JHALO                                         !<-- Maximum J shift for intra-task update
!
!-----------------------------------------------------------------------
!
      IF(ABS(I_SHIFT_CHILD)>I_DIFF_MAX                                  &  !<-- If true, gridpoints cannot be updated from
                  .OR.                                                  &  !    the same child task following a move.
         ABS(J_SHIFT_CHILD)>J_DIFF_MAX)THEN                                !<--
!
        RETURN                                                             !<-- Therefore exit.
!
      ENDIF
!
      kount_moves=kount_moves+1
!-----------------------------------------------------------------------
!***  Update those interior nest gridpoints that receive their new
!***  values from within their tasks' subdomain (memory).  Update 
!***  points include the task subdomain haloes but source points
!***  do not.
!-----------------------------------------------------------------------
!***                                 NOTE:
!***  If J_SHIFT_CHILD > 0 then we can shift data within each nest task
!***  in the normal sense for J, i.e., from smaller to larger.
!***  However if J_SHIFT_CHILD < 0 then we must shift data going from
!***  larger to smaller J or else we would cover up original data that 
!***  needed to be shifted later in the loop.  If J_SHIFT_CHILD = 0
!***  then this same notion is needed for I, i.e., we must loop from
!***  larger to smaller I for a westward shift.
!
!***  First establish some default values then refine them for specific
!***  directions of nest motion.
!-----------------------------------------------------------------------
!
      I_START=MAX(IMS,IDS)
      I_END  =MIN(IME,IDE)
      I_INC  =1
!
      J_START=MAX(JMS,JDS)
      J_END  =MIN(JME,JDE)
      J_INC  =1
!
!-------------------------------------
!***  Motion has southward component.
!-------------------------------------
!
      IF(J_SHIFT_CHILD<0 )THEN
        J_START=MIN(J_END,JTE_X-J_SHIFT_CHILD)                             !<-- Starting J (after move) for updates on same task
        J_END  =JTS_X-J_SHIFT_CHILD                                        !<-- Ending J (after move) for updates on same task
        J_INC  =-1                                                         !<-- J loop increment
        IF(I_SHIFT_CHILD==0)THEN
          I_START=MAX(IMS,IDS+NROWS_P_UPD_W)                               !<-- See introductory note above.
          I_END  =MIN(IME,IDE-NROWS_P_UPD_E)                               !<-- See introductory note above.
        ENDIF
      ENDIF
!
!-------------------------------------
!***  Motion has northward component.
!-------------------------------------
!
      IF(J_SHIFT_CHILD>0)THEN
        J_START=MAX(J_START,JTS_X-J_SHIFT_CHILD)                           !<-- Starting J (after move) for updates on same task
        J_END  =JTE_X-J_SHIFT_CHILD                                        !<-- Ending J (after move) for updates on same task
        IF(I_SHIFT_CHILD==0)THEN
          I_START=MAX(IMS,IDS+NROWS_P_UPD_W)                               !<-- See introductory note above.
          I_END  =MIN(IME,IDE-NROWS_P_UPD_E)                               !<-- See introductory note above.
        ENDIF
      ENDIF
!
!------------------------------------
!***  Motion has westward component.
!------------------------------------
!
      IF(I_SHIFT_CHILD<0)THEN
        IF(J_SHIFT_CHILD/=0)THEN                                           !<-- There is a north or south component of motion
          I_START=ITS_X-I_SHIFT_CHILD                                      !<-- Starting I (after move) for updates on same task
          I_END  =MIN(I_END,ITE_X-I_SHIFT_CHILD)                           !<-- Ending I (after move) for updates on same task
!
        ELSE                                                               !<-- Motion is due west
          I_START=MIN(I_END,ITE_X-I_SHIFT_CHILD)                           !<-- Starting I (after move) for updates on same task
          I_END  =ITS_X-I_SHIFT_CHILD                                      !<-- Ending I (after move) for updates on same task
          I_INC  =-1                                                       !<-- I loop increment
          J_START=MAX(JMS,JDS+NROWS_P_UPD_S)                               !<-- See introductory note above.
          J_END  =MIN(JME,JDE-NROWS_P_UPD_N)                               !<-- See introductory note above.
        ENDIF
      ENDIF
!
!------------------------------------
!***  Motion has eastward component.
!------------------------------------
!
      IF(I_SHIFT_CHILD>0)THEN
        I_START=MAX(I_START,ITS_X-I_SHIFT_CHILD)                           !<-- Starting I (after move) for updates on same task
        I_END  =ITE_X-I_SHIFT_CHILD                                        !<-- Ending I (after move) for updates on same task
        IF(J_SHIFT_CHILD==0)THEN
          J_START=MAX(JMS,JDS+NROWS_P_UPD_S)                               !<-- See introductory note above.
          J_END  =MIN(JME,JDE-NROWS_P_UPD_N)                               !<-- See introductory note above.
        ENDIF
      ENDIF
!
!-----------------------------------------------------------------------
!***  Shift the internal points that stay within the same MPI task.
!***  Loop through all pertinent 2-D and 3-D internal state variables
!***  on the moving domain.  
!-----------------------------------------------------------------------
!
!--------------
!***  H points
!--------------
!
      NUM_FIELDS=NUM_FIELDS_2D_H_I+NUM_FIELDS_2D_H_R+NUM_FIELDS_3D_H
!
      DO N_FIELD=1,NUM_FIELDS   
!
        CALL ESMF_FieldBundleGet(FIELDBUNDLE=MOVE_BUNDLE_H              &  !<-- Bundle holding the arrays for move updates
                                ,fieldIndex =N_FIELD                    &  !<-- Index of the Field in the Bundle
                                ,field      =HOLD_FIELD                 &  !<-- Field N_FIELD in the Bundle
                                ,rc         =RC )
!
        CALL ESMF_FieldGet(field   =HOLD_FIELD                          &  !<-- Field N_FIELD in the Bundle
                          ,dimCount=NUM_DIMS                            &  !<-- Is this Field 2-D or 3-D?
                          ,typeKind=DATATYPE                            &  !<-- Does this Field contain an integer or real array?
                          ,name    =FIELD_NAME                          &  !<-- The name of the Field
                          ,rc       =RC )
!
        N_REMOVE=INDEX(FIELD_NAME,SUFFIX_MOVE)
        FIELD_NAME=FIELD_NAME(1:N_REMOVE-1)
!
!-----------------------------------------------------------------------
        IF(NUM_DIMS==2)THEN
!-----------------------------------------------------------------------
!
          IF(DATATYPE==ESMF_TYPEKIND_R4)THEN                               !<-- The Real 2-D H-point arrays
!
            CALL ESMF_FieldGet(field    =HOLD_FIELD                     &  !<-- Field N_FIELD in the Bundle
                              ,localDe  =0                              &
                              ,farrayPtr=ARRAY_2D                       &  !<-- Dummy 2-D real array with Field's data
                              ,rc       =RC )
!
            DO J=J_START,J_END,J_INC
            DO I=I_START,I_END,I_INC
              ARRAY_2D(I,J)=ARRAY_2D(I+I_SHIFT_CHILD,J+J_SHIFT_CHILD)
            ENDDO
            ENDDO
!
          ELSEIF(DATATYPE==ESMF_TYPEKIND_I4)THEN                           !<-- The Integer 2-D H-point arrays
!
            CALL ESMF_FieldGet(field    =HOLD_FIELD                     &  !<-- Field N_FIELD in the Bundle
                              ,localDe  =0                              &
                              ,farrayPtr=IARRAY_2D                      &  !<-- Dummy 2-D integer array with Field's data
                              ,rc       =RC )
!
            DO J=J_START,J_END,J_INC
            DO I=I_START,I_END,I_INC
              IARRAY_2D(I,J)=IARRAY_2D(I+I_SHIFT_CHILD,J+J_SHIFT_CHILD)
            ENDDO
            ENDDO
!
          ENDIF
!
!-----------------------------------------------------------------------
        ELSEIF(NUM_DIMS==3)THEN                                            !<-- The (Real) 3-D H-point arrays
!-----------------------------------------------------------------------
!
          CALL ESMF_FieldGet(field      =HOLD_FIELD                     &  !<-- Field N in the Bundle
                            ,localDe    =0                              &
                            ,farrayPtr  =ARRAY_3D                       &  !<-- Dummy 3-D array with Field's data
                            ,totalLBound=LIMITS_LO                      &  !<-- Starting index in each dimension
                            ,totalUBound=LIMITS_HI                      &  !<-- Ending index in each dimension
                            ,rc         =RC )
!
          DO N=LIMITS_LO(3),LIMITS_HI(3)
            DO J=J_START,J_END,J_INC
            DO I=I_START,I_END,I_INC
              ARRAY_3D(I,J,N)=ARRAY_3D(I+I_SHIFT_CHILD,J+J_SHIFT_CHILD,N)
            ENDDO
            ENDDO
          ENDDO
!
        ENDIF
!
      ENDDO 
!
!
!--------------
!***  V points
!--------------
!
      NUM_FIELDS=NUM_FIELDS_2D_V+NUM_FIELDS_3D_V
!
      DO N_FIELD=1,NUM_FIELDS   
!
        CALL ESMF_FieldBundleGet(FIELDBUNDLE=MOVE_BUNDLE_V              &  !<-- Bundle holding the arrays for move updates
                                ,fieldIndex =N_FIELD                    &  !<-- Index of the Field in the Bundle
                                ,field      =HOLD_FIELD                 &  !<-- Field N_FIELD in the Bundle
                                ,rc         =RC )
!
        CALL ESMF_FieldGet(field   =HOLD_FIELD                          &  !<-- Field N_FIELD in the Bundle
                          ,dimCount=NUM_DIMS                            &  !<-- Is this Field 2-D or 3-D?
                          ,name    =FIELD_NAME                          &  !<-- The name of the Field
                          ,rc      =RC )
!
        N_REMOVE=INDEX(FIELD_NAME,SUFFIX_MOVE)
        FIELD_NAME=FIELD_NAME(1:N_REMOVE-1)
!
!-----------------------------------------------------------------------
        IF(NUM_DIMS==2)THEN
!-----------------------------------------------------------------------
!
          CALL ESMF_FieldGet(field    =HOLD_FIELD                       &  !<-- Field N_FIELD in the Bundle
                            ,localDe  =0                                &
                            ,farrayPtr=ARRAY_2D                         &  !<-- Dummy 2-D array with Field's data
                            ,rc       =RC )
!
          DO J=J_START,J_END,J_INC
          DO I=I_START,I_END,I_INC
            ARRAY_2D(I,J)=ARRAY_2D(I+I_SHIFT_CHILD,J+J_SHIFT_CHILD)
          ENDDO
          ENDDO
!
!-----------------------------------------------------------------------
        ELSEIF(NUM_DIMS==3)THEN                                            !<-- The (Real) 3-D V-point arrays
!-----------------------------------------------------------------------
!
          CALL ESMF_FieldGet(field      =HOLD_FIELD                     &  !<-- Field N in the Bundle
                            ,localDe    =0                              &
                            ,farrayPtr  =ARRAY_3D                       &  !<-- Dummy 3-D array with Field's data
                            ,totalLBound=LIMITS_LO                      &  !<-- Starting index in each dimension
                            ,totalUBound=LIMITS_HI                      &  !<-- Ending index in each dimension
                            ,rc         =RC )
!
          DO N=LIMITS_LO(3),LIMITS_HI(3)
            DO J=J_START,J_END,J_INC
            DO I=I_START,I_END,I_INC
              ARRAY_3D(I,J,N)=ARRAY_3D(I+I_SHIFT_CHILD,J+J_SHIFT_CHILD,N)
            ENDDO
            ENDDO
          ENDDO
!
        ENDIF
!
      ENDDO 
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE SHIFT_INTRA_TASK_DATA
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE SEND_INTER_TASK_DATA(I_SHIFT                           &
                                     ,J_SHIFT                           &
                                     ,MYPE                              &
                                     ,INPES                             &
                                     ,JNPES                             &
                                     ,MOVE_BUNDLE_H                     &
                                     ,NUM_FIELDS_2D_H_I                 &
                                     ,NUM_FIELDS_2D_H_R                 &
                                     ,NUM_FIELDS_3D_H                   &
                                     ,NUM_LEVELS_3D_H                   &
                                     ,MOVE_BUNDLE_V                     &
                                     ,NUM_FIELDS_2D_V                   &
                                     ,NUM_FIELDS_3D_V                   &
                                     ,NUM_LEVELS_3D_V                   &
                                     ,ITS,ITE,JTS,JTE                   &
                                     ,IMS,IME,JMS,JME                   &
                                     ,IDS,IDE,JDS,JDE                   &
                                      )
!
!-----------------------------------------------------------------------
!***  After a nest has moved, update those of its interior points 
!***  that still lie inside the footprint of the nest domain prior
!***  to the move but which now lie at an earth location previously 
!***  occupied by a point in a different one of the nest's tasks.
!***  In this subroutine each of those nest tasks with subdomains
!***  following a move that overlap the location of the nest domain
!***  preceding the move send data to other nest tasks whose
!***  subdomains now overlap its pre-move location.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      INTEGER(kind=KINT),INTENT(IN) :: I_SHIFT                           &  !<-- Nest moved this far in I on its grid.
                                      ,INPES                             &  !<-- # of fcst tasks in I on child grid
                                      ,J_SHIFT                           &  !<-- Nest moved this far in J on its grid.
                                      ,JNPES                             &  !<-- # of fcst tasks in J on child grid
                                      ,MYPE                              &  !<-- This task's local rank
                                      ,NUM_FIELDS_2D_H_I                 &  !<-- # of 2-D integer H variables to update
                                      ,NUM_FIELDS_2D_H_R                 &  !<-- # of 2-D real H variables to update
                                      ,NUM_FIELDS_3D_H                   &  !<-- # of 3-D H variables to update
                                      ,NUM_LEVELS_3D_H                   &  !<-- # of 2-D levels in all 3-D H update variables
                                      ,NUM_FIELDS_2D_V                   &  !<-- # of 2-D V variables to update
                                      ,NUM_FIELDS_3D_V                   &  !<-- # of 3-D V variables to update
                                      ,NUM_LEVELS_3D_V                   &  !<-- # of 2-D levels in all 3-D V update variables
!
                                      ,ITS,ITE,JTS,JTE                   &  !<-- Subdomain integration limits of this nest task
                                      ,IMS,IME,JMS,JME                   &  !<-- Subdomain memory limits of this nest task
                                      ,IDS,IDE,JDS,JDE                      !<-- Index limits of this nest's full domain
!
      TYPE(ESMF_FieldBundle),INTENT(INOUT) :: MOVE_BUNDLE_H              &  !<-- Bundle of internal state H arrays for updates
                                             ,MOVE_BUNDLE_V                 !<-- Bundle of internal state V arrays for updates
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: I,I_END,I_ID_END_SEARCH,I_ID_INC_SEARCH     &
                           ,I_ID_STA_SEARCH                             &
                           ,I_INC,I_START,I_TASK                        &
                           ,I_TASK_EAST,I_TASK_WEST                     &
                           ,I1,I2                                       &
                           ,ISEND_END,ISEND_START,ITE_X,ITS_X           & 
                           ,J,J_END,J_ID_END_SEARCH,J_ID_INC_SEARCH     &
                           ,J_ID_STA_SEARCH                             &
                           ,J_INC,J_START,J_TASK                        &
                           ,J_TASK_NORTH,J_TASK_SOUTH                   &
                           ,J1,J2                                       &
                           ,JSEND_END,JSEND_START,JTE_X,JTS_X           & 
                           ,KOUNT,KOUNT_INTEGER,KOUNT_REAL              &
                           ,L,N,N_FIELD,N_REMOVE,NF1,NF2                &
                           ,NUM_DIMS,NUM_FIELDS                         &
                           ,NUM_WORDS_IJ                                &
                           ,NUM_WORDS_INTEGER,NUM_WORDS_REAL
!
      INTEGER(kind=KINT) :: IDS_BND,IDE_BND,JDS_BND,JDE_BND
!
      INTEGER(kind=KINT) :: IERR,ISTAT,RC
!
      INTEGER(kind=KINT),DIMENSION(1:3) :: LIMITS_LO                    &
                                          ,LIMITS_HI
!
      INTEGER(kind=KINT),DIMENSION(1:9) :: ID_RECV
!
      INTEGER(kind=KINT),DIMENSION(MPI_STATUS_SIZE) :: JSTAT
!
      INTEGER(kind=KINT),DIMENSION(:,:),POINTER :: IARRAY_2D
!
      REAL(kind=KFPT),DIMENSION(:,:),POINTER :: ARRAY_2D
!
      REAL(kind=KFPT),DIMENSION(:,:,:),POINTER :: ARRAY_3D
!
      CHARACTER(len=99) :: FIELD_NAME
!
      TYPE(ESMF_Field) :: HOLD_FIELD
!
      TYPE(ESMF_TypeKind_Flag) :: DATATYPE
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  If the footprint of a nest task prior to a move has no points
!***  in common with the nest domain following the move then that task
!***  will have no data to send to other nest tasks and thus this
!***  routine is not relevant.  Remember that sending tasks will send
!***  to the outer 2 rows of recving tasks that lie on the domain
!***  boundary.  Those outer two rows cannot provide update data
!***  following a shift but they do receive update data.
!-----------------------------------------------------------------------
!
      IF(ITE<=IDS-1+I_SHIFT                                             &  !<-- Task footprint lies west of domain after east shift.
              .OR.                                                      &
         ITS>=IDE+1+I_SHIFT                                             &  !<-- Task footprint lies east of domain after west shift.
              .OR.                                                      &
         JTE<=JDS-1+J_SHIFT                                             &  !<-- Task footprint lies south of domain after north shift.
              .OR.                                                      &
         JTS>=JDE+1+J_SHIFT )THEN                                          !<-- Task footprint lies north of domain after south shift.
!
        RETURN                                                             !<-- Therefore exit.
!
      ENDIF
!
!-----------------------------------------------------------------------
!***  Update those interior nest gridpoints that receive their new 
!***  values from a different task within the nest domain.  All nest 
!***  tasks must determine which of their points' data must be sent to
!***  which other nest tasks.  Under normal circumstances the number
!***  of grid increments the nest shifts on its grid will not exceed
!***  a forecast task's dimensions.  If that is the case and if the
!***  nest motion has both I and J components then each task except
!***  those on the trailing edge will send to three nest tasks.  All
!***  trailing edge tasks except the trailing corner will send to one
!***  task; the corner will send to none.  If the motion has only an 
!***  I or only a J component then each non-trailing edge task will
!***  send to just one task and all trailing edge tasks will send to
!***  none.
!
!***  However in the general sense if the distance of the nest's
!***  motion exceeds the dimensions of any of its forecast tasks
!***  and the halo points of the receivers are included in those
!***  points to be updated (to avoid doing repeated halo exchanges
!***  after the move) then there are nine tasks that can potentially
!***  receive data from a given task who is sending.  If the location
!***  of the footprint of the sending task's pre-move position is
!***  represented by 'X' then the nine tasks are:
!
!
!               7      8     9
!
!
!               4      5     6
!                    X            
!                                
!               1      2     3  
!                              
!
!***  After the move note that we include as target points the 
!***  receiving tasks' halo points that lie within the sending task's
!***  integration domain before the move.  We account for the
!***  possibility that the width of the halo is greater than the 
!***  magnitude of the shift.
!
!***********************************************************************
!**************************   NOTE   ***********************************
!***********************************************************************
!
!***  HOWEVER it is critical to realize that neither H-pt nor V-pt
!***  variables at the north and east domain limits can be used
!***  in the send.  That is because the V-pt variables at
!***  IDE and JDE are not part of the integration.  Although the
!***  H-pt variables at the domain's IDE and JDE are part of the
!***  integration if the sender is on the nest domain's east or
!***  north boundary we must exclude them in the send too otherwise
!***  there would be occasions when the receiving tasks for the
!***  sender's H-pt data would not be identical to the receiving
!***  tasks for the sender's V-pt data and we must avoid that
!***  situation since it would greatly complicate the already very
!***  complicated bookkeeping.  This fact must also be accounted for
!***  in the intra-task shifts and in the sending/recving of parent
!***  update data after the nest moves.  At the same time there are
!***  a variety of variables that do not have valid integration values
!***  on the domain boundary so we do not want to let the inter-task
!***  shift process move those points into the interior.  Moreover
!***  the dynamical tendencies of temperature and the wind components
!***  are not defined on the 2nd row of the domain from the boundary
!***  so nest points that shift onto those locations cannot use the
!***  the intra- or inter-task updating.  Therefore the parent will
!***  update the outer two boundary rows of the pre-move footprint
!***  location in general.  For flexibility the code will use
!***  variables from the configure file to specify the depth into
!***  the footprint that the parent will supply update data.
!***  Use the following quantities to search for points that will
!***  be updated via the inter-task sends/recvs.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  The 'sender' is the current task that is executing this routine.
!***  The range of points within its subdomain that are valid for
!***  sending to other tasks are the following:
!-----------------------------------------------------------------------
!
      IDS_BND=IDS+NROWS_P_UPD_W
      IDE_BND=IDE-NROWS_P_UPD_E
      JDS_BND=JDS+NROWS_P_UPD_S
      JDE_BND=JDE-NROWS_P_UPD_N
!
      ITS_X=MAX(ITS,IDS_BND)
      ITE_X=MIN(ITE,IDE_BND)
      JTS_X=MAX(JTS,JDS_BND)
      JTE_X=MIN(JTE,JDE_BND)
!
!-----------------------------------------------------------------------
!***  Initialize to nonsense the ranks of tasks who will receive.
!***  There can be no more than nine.
!-----------------------------------------------------------------------
!
      DO N=1,9
        ID_RECV(N)=-1      
      ENDDO
!
!-----------------------------------------------------------------------
!***  Search for the tasks on this nest domain that will receive 
!***  intertask update data from this current task.  First look
!***  west/east then north/south.
!
!***  NOTE:  The outer two rows of points on the nest grid DO RECEIVE
!***         intra- and inter-task updates after the nest moves.  Those
!***         two rows of points DO NOT PROVIDE intra- and inter-task
!***         update data.
!-----------------------------------------------------------------------
!
      I_TASK_EAST=(MYPE/INPES+1)*INPES-1                                   !<-- Task on east end of sender's row.
      I_TASK_WEST=(MYPE/INPES)*INPES                                       !<-- Task on west end of sender's row.
!
      I_INC=SIGN(1,I_SHIFT)                                                !<-- +1 for eastward motion; -1 for westward motion
      J_INC=SIGN(1,J_SHIFT)                                                !<-- +1 for northward motion; -1 for southward motion
!
      IF(I_SHIFT>0)THEN                                                    !<-- For eastward move, search to the west.
        I_ID_STA_SEARCH=MYPE                                               !<-- Begin search with sender's column
        I_ID_END_SEARCH=I_TASK_WEST                                        !<-- Task on west end of sender's row.
        I_ID_INC_SEARCH=-I_INC                                             !<-- Task rank search increment in I (westward).
!
      ELSEIF(I_SHIFT<0)THEN                                                !<-- For westward move, search to the east.
        I_ID_STA_SEARCH=MYPE                                               !<-- Begin search with sender's column
        I_ID_END_SEARCH=I_TASK_EAST                                        !<-- Task on east end of sender's row.
        I_ID_INC_SEARCH=-I_INC                                             !<-- Task rank search increment in I (eastward).
!
      ELSEIF(I_SHIFT==0)THEN                                               !<-- No west/east motion
        I_ID_STA_SEARCH=MAX(MYPE-1,I_TASK_WEST)                            !<-- Search inc is +1 so begin 1 task to the west
        I_ID_END_SEARCH=MIN(MYPE+1,I_TASK_EAST)                            !<-- Search 3 columns due to halos on west/east sides
        I_ID_INC_SEARCH=1                                                  !<-- Task rank search increment
      ENDIF
!
      KOUNT=0                                                              !<-- Initialize counter of tasks that receive 
!                                                                               intertask updates from the current sender
!-----------------------------------------------------------------------
      search: DO I_TASK=I_ID_STA_SEARCH,I_ID_END_SEARCH,I_ID_INC_SEARCH
!-----------------------------------------------------------------------
!
        I_START=MAX(domain_int_state%LOCAL_ISTART(I_TASK)-IHALO,IDS)    &  !<-- West limit of potential receiver task subdomain
                +I_SHIFT
        I_END  =MIN(domain_int_state%LOCAL_IEND(I_TASK)  +IHALO,IDE)    &  !<-- East limit of potential receiver task subdomain
                +I_SHIFT
!
        IF(I_END>=ITS_X.AND.I_START<=ITE_X)THEN                            !<-- If so, task I_TASK's subdomain has moved onto searcher's
!
          J_TASK_NORTH=(JNPES-1)*INPES+MOD(I_TASK,INPES)                   !<-- Task on north end of I_TASK's column.
          J_TASK_SOUTH=MOD(I_TASK,INPES)                                   !<-- Task on south end of I_TASK's column.
!
          IF(J_SHIFT>0)THEN                                                !<-- For northward move, search to the south.
            J_ID_STA_SEARCH=I_TASK                                         !<-- Begin search in I_TASK's column
            J_ID_END_SEARCH=J_TASK_SOUTH                                   !<-- Task on south end of I_TASK's column.
            J_ID_INC_SEARCH=-J_INC*INPES                                   !<-- Task rank search increment in J (southward).
!
          ELSEIF(J_SHIFT<0)THEN                                            !<-- For southward move, search to the north.
            J_ID_STA_SEARCH=I_TASK                                         !<-- Begin search in I_TASK's column
            J_ID_END_SEARCH=J_TASK_NORTH                                   !<-- Task on north end of I_TASK's column.
            J_ID_INC_SEARCH=-J_INC*INPES                                   !<-- Task rank search increment in J (northward).
!
          ELSEIF(J_SHIFT==0)THEN                                           !<-- No south/north motion
            J_ID_STA_SEARCH=MAX(I_TASK-INPES,J_TASK_SOUTH)                 !<-- Begin search 1 task to the south due to halos
            J_ID_END_SEARCH=MIN(I_TASK+INPES,J_TASK_NORTH)                 !<-- Search 3 rows due to halos on north/south sides
            J_ID_INC_SEARCH=INPES                                          !<-- Task rank search increment
          ENDIF
!
          DO J_TASK=J_ID_STA_SEARCH,J_ID_END_SEARCH,J_ID_INC_SEARCH        !<-- If so then search north/south.
!         
            J_START=MAX(domain_int_state%LOCAL_JSTART(J_TASK)-JHALO,JDS) & !<-- South limit of potential receiver task subdomain
                    +J_SHIFT
            J_END  =MIN(domain_int_state%LOCAL_JEND(J_TASK)  +JHALO,JDE) & !<-- North limit of potential receiver task subdomain
                    +J_SHIFT
!
            IF(J_END>=JTS_X.AND.J_START<=JTE_X)THEN                        !<-- If so, task J_TASK's subdomain has moved onto searcher's
!
              KOUNT=KOUNT+1
              ID_RECV(KOUNT)=J_TASK                                        !<-- Save this task ID as a definite receiver of intertask data
!
            ENDIF
! 
          ENDDO
!
        ENDIF
!
!-----------------------------------------------------------------------
      ENDDO search
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Loop through the receive tasks to determine precisely which
!***  of their points need to be updated by this sender.  The sender
!***  does not send to itself.  It will update its own internal points
!***  in subroutine SHIFT_INTRA_TASK_DATA.
!-----------------------------------------------------------------------
!
      send_loop: DO N=1,9
!
!-----------------------------------------------------------------------
!
        check: IF(ID_RECV(N)>=0.AND.ID_RECV(N)/=MYPE)THEN                  !<-- Select the genuine receive tasks.
!
!-----------------------------------------------------------------------
!
          I1=MAX(domain_int_state%LOCAL_ISTART(ID_RECV(N))-IHALO,IDS)+I_SHIFT  !<-- West side of potential receiver N relative to footprint
          I2=MIN(domain_int_state%LOCAL_IEND  (ID_RECV(N))+IHALO,IDE)+I_SHIFT  !<-- East side of potential receiver N relative to footprint
          J1=MAX(domain_int_state%LOCAL_JSTART(ID_RECV(N))-JHALO,JDS)+J_SHIFT  !<-- South side of potential receiver N relative to footprint
          J2=MIN(domain_int_state%LOCAL_JEND  (ID_RECV(N))+JHALO,JDE)+J_SHIFT  !<-- North side of potential receiver N relative to footprint
!
!-----------------------------------------------------------------------
!
          sending: IF(I1<=ITE_X.AND.I2>=ITS_X                           &  !<-- Do any points in potential receiver task ID_RECV(N)
                              .AND.                                     &  !    lie within the footprint of the sender's location
                      J1<=JTE_X.AND.J2>=JTS_X) THEN                        !    prior to the move?
!
!-----------------------------------------------------------------------
!
            ISEND_START=MAX(I1,ITS_X)                                      !<-- West limit of task N's overlap within sender's footprint
            ISEND_END  =MIN(I2,ITE_X)                                      !<-- East limit of task N's overlap within sender's footprint
            JSEND_START=MAX(J1,JTS_X)                                      !<-- South limit of task N's overlap within sender's footprint
            JSEND_END  =MIN(J2,JTE_X)                                      !<-- North limit of task N's overlap within sender's footprint
            NUM_WORDS_IJ=(ISEND_END-ISEND_START+1)*                     &  !<-- Number of points (in the horizontal) in the overlap region.
                         (JSEND_END-JSEND_START+1)
!
!-----------------------------------------------------------------------
!***  Make sure the buffers have been received from the previous move
!***  so we can deallocate then reallocate them for the current move.
!-----------------------------------------------------------------------
!
!-------------------------------
!***  Real intertask shift data
!-------------------------------
!
            CALL MPI_WAIT(domain_int_state%HANDLE_SEND_INTER_REAL(N)    &  !<-- Handle for the ISend of inter-task real data on nest
                         ,JSTAT                                         &  !<-- MPI status
                         ,IERR )
!
            IF(ASSOCIATED(domain_int_state%SHIFT_DATA(N)%DATA_REAL))THEN
              DEALLOCATE(domain_int_state%SHIFT_DATA(N)%DATA_REAL)
            ENDIF
!
            NUM_WORDS_REAL=NUM_WORDS_IJ*(NUM_FIELDS_2D_H_R              &  !<-- Total # of real words in receiving task N's overlap
                                        +NUM_FIELDS_2D_V                &  !    region with sender task's pre-move footprint.
                                        +NUM_LEVELS_3D_H                &  !
                                        +NUM_LEVELS_3D_V)                  !<--
!
            ALLOCATE(domain_int_state%SHIFT_DATA(N)%DATA_REAL(1:NUM_WORDS_REAL))
!
!----------------------------------
!***  Integer intertask shift data
!----------------------------------
!
            CALL MPI_WAIT(domain_int_state%HANDLE_SEND_INTER_INT(N)     &  !<-- Handle for the ISend of inter-task integer data on nest
                         ,JSTAT                                         &  !<-- MPI status
                         ,IERR )
!
            IF(ASSOCIATED(domain_int_state%SHIFT_DATA(N)%DATA_INTEGER))THEN
              DEALLOCATE(domain_int_state%SHIFT_DATA(N)%DATA_INTEGER)
            ENDIF
!
            NUM_WORDS_INTEGER=NUM_WORDS_IJ*NUM_FIELDS_2D_H_I               !<-- Total # of integer words in receiving task N's overlap
!
            ALLOCATE(domain_int_state%SHIFT_DATA(N)%DATA_INTEGER(1:NUM_WORDS_INTEGER))
!
!-----------------------------------------------------------------------
!***  Loop through the internal state variables lifting out the 
!***  data that lies in each receiving task's overlap region in 
!***  the sender's footprint.  The indices below are with respect
!***  to the sender's footprint.  Store the data in a 1-D array
!***  so that it can be given to MPI_ISEND for transfer to the
!***  receiver tasks.
!-----------------------------------------------------------------------
!
            KOUNT_REAL=0
            KOUNT_INTEGER=0
!
!--------------
!***  H points
!--------------
!
            NUM_FIELDS=NUM_FIELDS_2D_H_I+NUM_FIELDS_2D_H_R+NUM_FIELDS_3D_H
!
            DO N_FIELD=1,NUM_FIELDS
!
              CALL ESMF_FieldBundleGet(FIELDBUNDLE=MOVE_BUNDLE_H        &  !<-- Bundle holding the arrays for move updates
                                      ,fieldIndex =N_FIELD              &  !<-- Index of the Field in the Bundle
                                      ,field      =HOLD_FIELD           &  !<-- Field N_FIELD in the Bundle
                                      ,rc         =RC )
!
              CALL ESMF_FieldGet(field    =HOLD_FIELD                   &  !<-- Field N_FIELD in the Bundle
                                 ,dimCount=NUM_DIMS                     &  !<-- Is this Field 2-D or 3-D?
                                 ,typeKind=DATATYPE                     &  !<-- Does the Field contain an integer or real array?
                                 ,name    =FIELD_NAME                   &  !<-- This Field's name
                                 ,rc      =RC )
!
              N_REMOVE=INDEX(FIELD_NAME,SUFFIX_MOVE)
              FIELD_NAME=FIELD_NAME(1:N_REMOVE-1)
!
              IF(NUM_DIMS==2)THEN
!
                IF(DATATYPE==ESMF_TYPEKIND_R4)THEN                         !<-- Real 2-D H-point arrays
!
                  CALL ESMF_FieldGet(field    =HOLD_FIELD               &  !<-- Field N_FIELD in the Bundle
                                    ,localDe  =0                        &
                                    ,farrayPtr=ARRAY_2D                 &  !<-- Dummy 2-D real array with Field's data
                                    ,rc       =RC )
!
                  DO J=JSEND_START,JSEND_END
                  DO I=ISEND_START,ISEND_END
                    KOUNT_REAL=KOUNT_REAL+1
                    domain_int_state%SHIFT_DATA(N)%DATA_REAL(KOUNT_REAL)=ARRAY_2D(I,J) !<-- Sender collects its 2-D Real H data
                                                                                       !    in overlap region.
                  ENDDO
                  ENDDO
!
                ELSEIF(DATATYPE==ESMF_TYPEKIND_I4)THEN                     !<-- Integer 2-D H-point arrays
!
                  CALL ESMF_FieldGet(field    =HOLD_FIELD               &  !<-- Field N_FIELD in the Bundle
                                    ,localDe  =0                        &
                                    ,farrayPtr=IARRAY_2D                &  !<-- Dummy 2-D integer array with Field's data
                                    ,rc       =RC )
!
                  DO J=JSEND_START,JSEND_END
                  DO I=ISEND_START,ISEND_END
                    KOUNT_INTEGER=KOUNT_INTEGER+1
                    domain_int_state%SHIFT_DATA(N)%DATA_INTEGER(KOUNT_INTEGER)=IARRAY_2D(I,J)  !<-- Sender collects its 2-D Integer H data
                                                                                               !    in overlap region.
                  ENDDO
                  ENDDO
!
                ENDIF
!
              ELSEIF(NUM_DIMS==3)THEN                                      !<-- (Real) 3-D H-point arrays
!
                CALL ESMF_FieldGet(field      =HOLD_FIELD               &  !<-- Field N_FIELD in the Bundle
                                  ,localDe    =0                        &
                                  ,farrayPtr  =ARRAY_3D                 &  !<-- Dummy 3-D array with Field's data
                                  ,totalLBound=LIMITS_LO                &  !<-- Starting index in each dimension
                                  ,totalUBound=LIMITS_HI                &  !<-- Ending index in each dimension
                                  ,rc         =RC )
!
                DO L=LIMITS_LO(3),LIMITS_HI(3)
                  DO J=JSEND_START,JSEND_END
                  DO I=ISEND_START,ISEND_END
                    KOUNT_REAL=KOUNT_REAL+1
                    domain_int_state%SHIFT_DATA(N)%DATA_REAL(KOUNT_REAL)=ARRAY_3D(I,J,L)  !<-- Sender collects its 3-D (Real) H data
                                                                                          !    in overlap region.
                  ENDDO
                  ENDDO
                ENDDO
!
              ENDIF
!
            ENDDO
!
!
!--------------
!***  V points
!--------------
!
            NUM_FIELDS=NUM_FIELDS_2D_V+NUM_FIELDS_3D_V
!
            DO N_FIELD=1,NUM_FIELDS
!
              CALL ESMF_FieldBundleGet(FIELDBUNDLE=MOVE_BUNDLE_V        &  !<-- Bundle holding the arrays for move updates
                                      ,fieldIndex =N_FIELD              &  !<-- Index of the Field in the Bundle
                                      ,field      =HOLD_FIELD           &  !<-- Field N_FIELD in the Bundle
                                      ,rc         =RC )
!
              CALL ESMF_FieldGet(field    =HOLD_FIELD                   &  !<-- Field N_FIELD in the Bundle
                                 ,dimCount=NUM_DIMS                     &  !<-- Is this Field 2-D or 3-D?
                                 ,name      =FIELD_NAME                 &  !<-- This Field's name
                                 ,rc      =RC )
!
              N_REMOVE=INDEX(FIELD_NAME,SUFFIX_MOVE)
              FIELD_NAME=FIELD_NAME(1:N_REMOVE-1)
!
              IF(NUM_DIMS==2)THEN
!
                CALL ESMF_FieldGet(field    =HOLD_FIELD                 &  !<-- Field N_FIELD in the Bundle
                                  ,localDe  =0                          &
                                  ,farrayPtr=ARRAY_2D                   &  !<-- Dummy 2-D array with Field's data
                                  ,rc       =RC )
!
                DO J=JSEND_START,JSEND_END
                DO I=ISEND_START,ISEND_END
                  KOUNT_REAL=KOUNT_REAL+1
                  domain_int_state%SHIFT_DATA(N)%DATA_REAL(KOUNT_REAL)=ARRAY_2D(I,J)  !<-- Sender collects its 2-D (Real) V data
                                                                                      !    in overlap region.
                           
                ENDDO
                ENDDO
!
              ELSEIF(NUM_DIMS==3)THEN
!
                CALL ESMF_FieldGet(field      =HOLD_FIELD               &  !<-- Field N_FIELD in the Bundle
                                  ,localDe    =0                        &
                                  ,farrayPtr  =ARRAY_3D                 &  !<-- Dummy 3-D array with Field's data
                                  ,totalLBound=LIMITS_LO                &  !<-- Starting index in each dimension
                                  ,totalUBound=LIMITS_HI                &  !<-- Ending index in each dimension
                                  ,rc         =RC )
!
                DO L=LIMITS_LO(3),LIMITS_HI(3)
                  DO J=JSEND_START,JSEND_END
                  DO I=ISEND_START,ISEND_END
                    KOUNT_REAL=KOUNT_REAL+1
                    domain_int_state%SHIFT_DATA(N)%DATA_REAL(KOUNT_REAL)=ARRAY_3D(I,J,L)  !<-- Sender collects its 3-D (Real) V data
                                                                                          !    in overlap region.
                  ENDDO
                  ENDDO
                ENDDO
!
              ENDIF
!
            ENDDO
!
!-----------------------------------------------------------------------
!***  Send all the real data.
!-----------------------------------------------------------------------
!
            CALL MPI_ISSEND(domain_int_state%SHIFT_DATA(N)%DATA_REAL    &  !<-- All inter-task shift Real data for task N
                           ,KOUNT_REAL                                  &  !<-- # of words in the Real data
                           ,MPI_REAL                                    &  !<-- The words are real
                           ,ID_RECV(N)                                  &  !<-- The nest task to which the sender is sending
                           ,KOUNT_REAL                                  &  !<-- Use the word count as the tag
                           ,COMM_FCST_TASKS(MY_DOMAIN_ID)               &  !<-- The MPI intracommunicator for this domain's fcst tasks
                           ,domain_int_state%HANDLE_SEND_INTER_REAL(N)  &  !<-- Handle for this ISend
                           ,IERR )
!
!-----------------------------------------------------------------------
!***  Send all the integer data.
!-----------------------------------------------------------------------
!
            CALL MPI_ISSEND(domain_int_state%SHIFT_DATA(N)%DATA_INTEGER &  !<-- All inter-task shift Integer data for task N
                           ,KOUNT_INTEGER                               &  !<-- # of words in the Integer data
                           ,MPI_INTEGER                                 &  !<-- The words are integer
                           ,ID_RECV(N)                                  &  !<-- The nest task to which the sender is sending
                           ,KOUNT_INTEGER                               &  !<-- Use the word count as the tag
                           ,COMM_FCST_TASKS(MY_DOMAIN_ID)               &  !<-- The MPI intracommunicator for this domain's fcst tasks
                           ,domain_int_state%HANDLE_SEND_INTER_INT(N)   &  !<-- Handle for this ISend
                           ,IERR )
!
!-----------------------------------------------------------------------
!
          ENDIF sending
!
!-----------------------------------------------------------------------
!
        ENDIF check
!
!-----------------------------------------------------------------------
!
      ENDDO send_loop
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE SEND_INTER_TASK_DATA
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE RECV_INTER_TASK_DATA(I_SHIFT                           &
                                     ,J_SHIFT                           &
                                     ,MYPE                              &
                                     ,INPES                             &
                                     ,JNPES                             &
                                     ,MOVE_BUNDLE_H                     &
                                     ,NUM_FIELDS_2D_H_I                 &
                                     ,NUM_FIELDS_2D_H_R                 &
                                     ,NUM_FIELDS_3D_H                   &
                                     ,NUM_LEVELS_3D_H                   &
                                     ,MOVE_BUNDLE_V                     &
                                     ,NUM_FIELDS_2D_V                   &
                                     ,NUM_FIELDS_3D_V                   &
                                     ,NUM_LEVELS_3D_V                   &
                                     ,ITS,ITE,JTS,JTE                   &
                                     ,IMS,IME,JMS,JME                   &
                                     ,IDS,IDE,JDS,JDE                   &
                                      )
!
!-----------------------------------------------------------------------
!***  After a nest has moved, update those of its interior points 
!***  that still lie inside the footprint of the nest domain prior
!***  to the move but which now lie at an earth location previously 
!***  occupied by a point in a different one of the nest's tasks.
!***  In this subroutine those nest tasks (including their halo points)
!***  after a move that overlap the pre-move location of the integration
!***  subdomain of another nest task now receive their updata date from
!***  the other nest task.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      INTEGER(kind=KINT),INTENT(IN) :: I_SHIFT                           &  !<-- Nest moved this far in I on its grid.
                                      ,J_SHIFT                           &  !<-- Nest moved this far in J on its grid.
                                      ,INPES                             &  !<-- # of fcst tasks in I on child grid
                                      ,JNPES                             &  !<-- # of fcst tasks in J on child grid
                                      ,MYPE                              &  !<-- This task's local rank
                                      ,NUM_FIELDS_2D_H_I                 &  !<-- # of 2-D integer H variables to update
                                      ,NUM_FIELDS_2D_H_R                 &  !<-- # of 2-D real H variables to update
                                      ,NUM_FIELDS_3D_H                   &  !<-- # of 3-D internal state H variables to update
                                      ,NUM_LEVELS_3D_H                   &  !<-- # of 2-D levels in all 3-D H update variables
                                      ,NUM_FIELDS_2D_V                   &  !<-- # of 2-D internal state V variables to update
                                      ,NUM_FIELDS_3D_V                   &  !<-- # of 3-D internal state V variables to update
                                      ,NUM_LEVELS_3D_V                   &  !<-- # of 2-D levels in all 3-D V update variables
!
                                      ,ITS,ITE,JTS,JTE                   &  !<-- Subdomain integration limits of this nest task
                                      ,IMS,IME,JMS,JME                   &  !<-- Subdomain memory limits of this nest task
                                      ,IDS,IDE,JDS,JDE                      !<-- Index limits of this nest's full domain
!
      TYPE(ESMF_FieldBundle),INTENT(INOUT) :: MOVE_BUNDLE_H              &  !<-- Bundle of internal state H arrays to update
                                             ,MOVE_BUNDLE_V                 !<-- Bundle of internal state V arrays to update
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: I,I_END_X,I_ID_END_SEARCH,I_ID_INC_SEARCH    &
                           ,I_ID_STA_SEARCH                              &
                           ,I_INC,I_START,I_START_X,I_TASK               &
                           ,I_TASK_EAST,I_TASK_WEST                      &
                           ,I1,I2                                        &
                           ,IRECV_END,IRECV_START                        &
                           ,ITS_X,ITE_X                                  &
                           ,J,J_END_X,J_ID_END_SEARCH,J_ID_INC_SEARCH    &
                           ,J_ID_STA_SEARCH                              &
                           ,J_INC,J_START,J_START_X,J_TASK               &
                           ,J_TASK_NORTH,J_TASK_SOUTH                    &
                           ,J1,J2                                        &
                           ,JRECV_END,JRECV_START                        &
                           ,JTE_X,JTS_X                                  &
                           ,KOUNT,KOUNT_INTEGER,KOUNT_REAL               &
                           ,L,N,N_FIELD,N_REMOVE,NF1,NF2                 &
                           ,NUM_DIMS,NUM_FIELDS                          &
                           ,NUM_WORDS_IJ                                 &
                           ,NUM_WORDS_INTEGER,NUM_WORDS_REAL
!
      INTEGER(kind=KINT) :: IDE_BND,IDS_BND,JDE_BND,JDS_BND
!
      INTEGER(kind=KINT) :: IERR,RC
!
      INTEGER(kind=KINT),DIMENSION(1:3) :: LIMITS_LO                     &
                                          ,LIMITS_HI
!
      INTEGER(kind=KINT),DIMENSION(1:9) :: ID_SEND
!
      INTEGER(kind=KINT),DIMENSION(MPI_STATUS_SIZE) :: JSTAT
!
      INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE :: RECV_INTEGER_DATA
!
      INTEGER(kind=KINT),DIMENSION(:,:),POINTER :: IARRAY_2D
!
      REAL(kind=KFPT),DIMENSION(:),ALLOCATABLE :: RECV_REAL_DATA
!
      REAL(kind=KFPT),DIMENSION(:,:),POINTER :: ARRAY_2D
!
      REAL(kind=KFPT),DIMENSION(:,:,:),POINTER :: ARRAY_3D
!
      CHARACTER(len=99) :: FIELD_NAME
!
      TYPE(ESMF_Field) :: HOLD_FIELD
!
      TYPE(ESMF_TypeKind_Flag) :: DATATYPE
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  We will essentially use the inverse of the logic used in 
!***  SEND_INTER_TASK_DATA in order to find the set of nine 
!***  potential sending tasks from which to receive update data
!***  in each nest task that after a move at least partially remains
!***  within the outline of the nest domain prior to the move.
!***  See the comments in that subroutine.  'Footprint' will always
!***  refer to domain/subdomain positions prior to a move.
!
!***  The following four variables are the limits on the receiving
!***  task subdomain of the points that can receive inter-task
!***  updates.  Remember that while the outer two rows of a nest
!***  domain cannot provide intra- and inter-task update data
!***  the points in those two rows certainly do receive intra- and
!***  inter-task update data.
!-----------------------------------------------------------------------
!
      I_START_X=MAX(IMS,IDS)
      I_END_X  =MIN(IME,IDE)
      J_START_X=MAX(JMS,JDS)
      J_END_X  =MIN(JME,JDE)
!
!-----------------------------------------------------------------------
!***  If the subdomain of a nest task after a move has no points in
!***  common with the nest domain footprint prior to the move then 
!***  that task will have no data to receive from other nest tasks 
!***  and thus this routine is not relevant.
!-----------------------------------------------------------------------
!
      IF(I_END_X+I_SHIFT<=IDS+NROWS_P_UPD_W-1                           &  !<-- Task subdomain lies west of footprint after west shift.
                 .OR.                                                   &
         I_START_X+I_SHIFT>=IDE-NROWS_P_UPD_E+1                         &  !<-- Task subdomain lies east of footprint after east shift.
                 .OR.                                                   &
         J_END_X+J_SHIFT<=JDS+NROWS_P_UPD_S-1                           &  !<-- Task subdomain lies south of footprint after south shift.
                 .OR.                                                   &
         J_START_X+J_SHIFT>=JDE-NROWS_P_UPD_N+1 )THEN                      !<-- Task subdomain lies north of footprint after north shift.
!
        RETURN                                                             !<-- Therefore exit.
!
      ENDIF
!
!-----------------------------------------------------------------------
!***  Initialize to nonsense the ranks of tasks who might send data.
!***  There can be no more than nine.
!-----------------------------------------------------------------------
!
      DO N=1,9
        ID_SEND(N)=-1
      ENDDO
!
!-----------------------------------------------------------------------
!***  Search for the tasks that will send intertask update data to
!***  the current search task.  First look west/east then north/south.
!
!***  NOTE:  The search is done with respect to the grid indices on
!***         the subdomain of the receiver's position after the move.
!-----------------------------------------------------------------------
!
      I_INC=SIGN(1,I_SHIFT)                                                !<-- +1 for eastward motion; -1 for westward motion
      J_INC=SIGN(1,J_SHIFT)                                                !<-- +1 for northward motion; -1 for southward motion
!
      I_TASK_EAST=(MYPE/INPES+1)*INPES-1                                   !<-- Task on east end of receiver's row.
      I_TASK_WEST=(MYPE/INPES)*INPES                                       !<-- Task on west end of receiver's row.
!
      IF(I_SHIFT>0)THEN                                                    !<-- For eastward move, search to the east.
        I_ID_STA_SEARCH=MYPE                                               !<-- Begin search with current task's column.
        I_ID_END_SEARCH=I_TASK_EAST                                        !<-- Task on east end of receiver's row.
        I_ID_INC_SEARCH=I_INC                                              !<-- Task rank search increment in I (eastward).
!
      ELSEIF(I_SHIFT<0)THEN                                                !<-- For westward move, search to the west.
        I_ID_STA_SEARCH=MYPE                                               !<-- Begin search with current task's column.
        I_ID_END_SEARCH=I_TASK_WEST                                        !<-- Task on west end of receiver's row.
        I_ID_INC_SEARCH=I_INC                                              !<-- Task rank search increment in I (westward).
!
      ELSEIF(I_SHIFT==0)THEN                                               !<-- No west/east motion
        I_ID_STA_SEARCH=MAX(MYPE-1,I_TASK_WEST)                            !<-- Begin search 1 task to the west due to halos
        I_ID_END_SEARCH=MIN(MYPE+1,I_TASK_EAST)                            !<-- End search 1 task to the east due to halos
        I_ID_INC_SEARCH=1                                                  !<-- Task rank search increment
      ENDIF
!
      KOUNT=0                                                              !<-- Initialize counter of tasks that will send intertask data
!
      IDS_BND=IDS+NROWS_P_UPD_W
      IDE_BND=IDE-NROWS_P_UPD_E
      JDS_BND=JDS+NROWS_P_UPD_S
      JDE_BND=JDE-NROWS_P_UPD_N
!
!-----------------------------------------------------------------------
      search: DO I_TASK=I_ID_STA_SEARCH,I_ID_END_SEARCH,I_ID_INC_SEARCH
!-----------------------------------------------------------------------
!
        ITS_X=MAX(domain_int_state%LOCAL_ISTART(I_TASK),IDS_BND)           !<-- West limit of task I_TASK's integration region
        ITE_X=MIN(domain_int_state%LOCAL_IEND(I_TASK)  ,IDE_BND)           !<-- East limit of task I_TASK's integration region
!
        IF(I_END_X+I_SHIFT>=ITS_X.AND.I_START_X+I_SHIFT<=ITE_X)THEN        !<-- If so, some of current task's subdomain moved onto I_TASK's
!
          J_TASK_NORTH=(JNPES-1)*INPES+MOD(I_TASK,INPES)                   !<-- Task on north end of I_TASK's column.
          J_TASK_SOUTH=MOD(I_TASK,INPES)                                   !<-- Task on south end of I_TASK's column.
!
          IF(J_SHIFT>0)THEN                                                !<-- For northward move, search to the north.
            J_ID_STA_SEARCH=I_TASK                                         !<-- Begin search with I_TASK
            J_ID_END_SEARCH=J_TASK_NORTH                                   !<-- Task on north end of I_TASK's column.
            J_ID_INC_SEARCH=J_INC*INPES                                    !<-- Task rank search increment in J (northward).
!
          ELSEIF(J_SHIFT<0)THEN                                            !<-- For southward move, search to the south.
            J_ID_STA_SEARCH=I_TASK                                         !<-- Begin search with I_TASK
            J_ID_END_SEARCH=J_TASK_SOUTH                                   !<-- Task on south end of I_TASK's column.
            J_ID_INC_SEARCH=J_INC*INPES                                    !<-- Task rank search increment in J (southward).
!
          ELSEIF(J_SHIFT==0)THEN                                           !<-- No south/north motion
            J_ID_STA_SEARCH=MAX(I_TASK-INPES,J_TASK_SOUTH)                 !<-- Due to halos begin 1 task to the south
            J_ID_END_SEARCH=MIN(I_TASK+INPES,J_TASK_NORTH)                 !<-- And end search 1 task to the north     
            J_ID_INC_SEARCH=INPES                                          !<-- Task rank search increment
          ENDIF
!
          DO J_TASK=J_ID_STA_SEARCH,J_ID_END_SEARCH,J_ID_INC_SEARCH        !<-- If so then search north/south.
!
            JTS_X=MAX(domain_int_state%LOCAL_JSTART(J_TASK),JDS_BND)       !<-- South limit of task J_TASK integration region
            JTE_X=MIN(domain_int_state%LOCAL_JEND(J_TASK)  ,JDE_BND)       !<-- North limit of task J_TASK integration region
!
            IF(J_END_X+J_SHIFT>=JTS_X.AND.J_START_X+J_SHIFT<=JTE_X)THEN    !<-- If so, current task has moved onto J_TASK's subdomain
!
              KOUNT=KOUNT+1
              ID_SEND(KOUNT)=J_TASK                                        !<-- Save this task ID as a definite sender of intertask data
            ENDIF
!
          ENDDO
!
        ENDIF
!
!-----------------------------------------------------------------------
      ENDDO search
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Loop through the nine potential send tasks to determine which
!***  of their points are needed for updating points in this receiver.
!***  The current task executing this routine will not receive from
!***  itself.  It will update its own internal points in subroutine
!***  SHIFT_INTRA_TASK_DATA.
!-----------------------------------------------------------------------
!
      recv_loop: DO N=1,9
!
!-----------------------------------------------------------------------
!
        check: IF(ID_SEND(N)>=0.AND.ID_SEND(N)/=MYPE)THEN                  !<-- Potential send task has points inside nest boundary.
!
!-----------------------------------------------------------------------
!
          I1=MAX(domain_int_state%LOCAL_ISTART(ID_SEND(N))              &  !<-- West side of potential sender N relative to subdomain
                                           ,IDS+NROWS_P_UPD_W)-I_SHIFT
          I2=MIN(domain_int_state%LOCAL_IEND(ID_SEND(N))                &  !<-- East side of potential sender N relative to subdomain
                                           ,IDE-NROWS_P_UPD_E)-I_SHIFT
          J1=MAX(domain_int_state%LOCAL_JSTART(ID_SEND(N))              &  !<-- South side of potential sender N relative to subdomain
                                           ,JDS+NROWS_P_UPD_S)-J_SHIFT
          J2=MIN(domain_int_state%LOCAL_JEND(ID_SEND(N))                &  !<-- North side of potential sender N relative to subdomain
                                           ,JDE-NROWS_P_UPD_N)-J_SHIFT
!
!-----------------------------------------------------------------------
!
          recving: IF(I_START_X<=I2.AND.I_END_X>=I1                     &  !<-- Does any of the receiver's subdomain after the move
                                  .AND.                                 &  !    intersect potential sender N's subdomain location
                      J_START_X<=J2.AND.J_END_X>=J1) THEN                  !    at its pre-move location?
!
!-----------------------------------------------------------------------
!
            IRECV_START=MAX(I1,I_START_X)                                  !<-- West limit of task N's overlap within receiver's subdomain
            IRECV_END  =MIN(I2,I_END_X)                                    !<-- East limit of task N's overlap within receiver's subdomain
            JRECV_START=MAX(J1,J_START_X)                                  !<-- South limit of task N's overlap within receiver's subdomain
            JRECV_END  =MIN(J2,J_END_X)                                    !<-- North limit of task N's overlap within receiver's subdomain
!
            NUM_WORDS_IJ=(IRECV_END-IRECV_START+1)*                     &  !<-- Number of points (in the horizontal) in the overlap region.
                         (JRECV_END-JRECV_START+1)
            NUM_WORDS_REAL=NUM_WORDS_IJ*(NUM_FIELDS_2D_H_R              &  !<-- Total # of Real words in sending task N's overlap
                                        +NUM_FIELDS_2D_V                &  !    with receiver task's subdomain for all update
                                        +NUM_LEVELS_3D_H                &  !    variables.
                                        +NUM_LEVELS_3D_V)                  !<--   
!
            NUM_WORDS_INTEGER=NUM_WORDS_IJ*NUM_FIELDS_2D_H_I               !<-- Total # of Integer words in sending task N's overlap
!
            IF(ALLOCATED(RECV_REAL_DATA))DEALLOCATE(RECV_REAL_DATA)
            ALLOCATE(RECV_REAL_DATA(1:NUM_WORDS_REAL))                     !<-- Allocate the Recv buffer for Real data
!
            IF(ALLOCATED(RECV_INTEGER_DATA))DEALLOCATE(RECV_INTEGER_DATA)
            ALLOCATE(RECV_INTEGER_DATA(1:NUM_WORDS_INTEGER))               !<-- Allocate the Recv buffer for Integer data
!
!-----------------------------------------------------------------------
!***  Receive all Real update data from nest task N.
!-----------------------------------------------------------------------
!
            CALL MPI_RECV(RECV_REAL_DATA                                &  !<-- All Real inter-task shift data from task N
                         ,NUM_WORDS_REAL                                &  !<-- # of words in the Real data
                         ,MPI_REAL                                      &  !<-- The words are Real
                         ,ID_SEND(N)                                    &  !<-- The nest task who is sending
                         ,NUM_WORDS_REAL                                &  !<-- Use the word count as the tag
                         ,COMM_FCST_TASKS(MY_DOMAIN_ID)                 &  !<-- The MPI intracommunicator for this domain's fcst tasks
                         ,JSTAT                                         &  !<-- MPI status object
                         ,IERR )
!
!-----------------------------------------------------------------------
!***  Receive all Integer update data from nest task N.
!-----------------------------------------------------------------------
!
            CALL MPI_RECV(RECV_INTEGER_DATA                             &  !<-- All Integer inter-task shift data from task N
                         ,NUM_WORDS_INTEGER                             &  !<-- # of Integer words in the data
                         ,MPI_INTEGER                                   &  !<-- The words are Integer
                         ,ID_SEND(N)                                    &  !<-- The nest task who is sending
                         ,NUM_WORDS_INTEGER                             &  !<-- Use the word count as the tag
                         ,COMM_FCST_TASKS(MY_DOMAIN_ID)                 &  !<-- The MPI intracommunicator for this domain's fcst tasks
                         ,JSTAT                                         &  !<-- MPI status object
                         ,IERR )
!
!-----------------------------------------------------------------------
!***  Incorporate all update data.
!-----------------------------------------------------------------------
!
            KOUNT_REAL=0
            KOUNT_INTEGER=0
!
!--------------
!***  H points
!--------------
!
            NUM_FIELDS=NUM_FIELDS_2D_H_I+NUM_FIELDS_2D_H_R              &
                      +NUM_FIELDS_3D_H
!
            DO N_FIELD=1,NUM_FIELDS
!
              CALL ESMF_FieldBundleGet(FIELDBUNDLE=MOVE_BUNDLE_H        &  !<-- Bundle holding the H arrays for move updates
                                      ,fieldIndex =N_FIELD              &  !<-- Index of the Field in the Bundle
                                      ,field      =HOLD_FIELD           &  !<-- Field N_FIELD in the Bundle
                                      ,rc         =RC )
!
              CALL ESMF_FieldGet(field   =HOLD_FIELD                    &  !<-- Field N_FIELD in the Bundle
                                ,dimCount=NUM_DIMS                      &  !<-- Is this Field 2-D or 3-D?
                                ,typeKind=DATATYPE                      &  !<-- Does this Field contain an integer or real array?
                                ,name    =FIELD_NAME                    &  !<-- This Field's name
                                ,rc      =RC )
!
              N_REMOVE=INDEX(FIELD_NAME,SUFFIX_MOVE)
              FIELD_NAME=FIELD_NAME(1:N_REMOVE-1)
!
              IF(NUM_DIMS==2)THEN
!
                IF(DATATYPE==ESMF_TYPEKIND_I4)THEN
!
                  CALL ESMF_FieldGet(field    =HOLD_FIELD               &  !<-- Field N_FIELD in the Bundle
                                    ,localDe  =0                        &
                                    ,farrayPtr=IARRAY_2D                &  !<-- Dummy 2-D integer array with the Field's data
                                    ,rc       =RC )
!
                  DO J=JRECV_START,JRECV_END
                  DO I=IRECV_START,IRECV_END
                    KOUNT_INTEGER=KOUNT_INTEGER+1
                    IARRAY_2D(I,J)=RECV_INTEGER_DATA(KOUNT_INTEGER)        !<-- Task updates 2-D Integer H data in overlap region
                  ENDDO
                  ENDDO
!
                ELSEIF(DATATYPE==ESMF_TYPEKIND_R4)THEN
!
                  CALL ESMF_FieldGet(field    =HOLD_FIELD               &  !<-- Field N_FIELD in the Bundle
                                    ,localDe  =0                        &
                                    ,farrayPtr=ARRAY_2D                 &  !<-- Dummy 2-D real array with the Field's data
                                    ,rc       =RC )
!
                  DO J=JRECV_START,JRECV_END
                  DO I=IRECV_START,IRECV_END
                    KOUNT_REAL=KOUNT_REAL+1
                    ARRAY_2D(I,J)=RECV_REAL_DATA(KOUNT_REAL)               !<-- Task updates 2-D Real H data in overlap region
                  ENDDO
                  ENDDO
!
                ENDIF
!
              ELSEIF(NUM_DIMS==3)THEN
!
                CALL ESMF_FieldGet(field      =HOLD_FIELD               &  !<-- Field N_FIELD in the Bundle
                                  ,localDe    =0                        &
                                  ,farrayPtr  =ARRAY_3D                 &  !<-- Dummy 3-D array with the Field's data
                                  ,totalLBound=LIMITS_LO                &  !<-- Starting index in each dimension
                                  ,totalUBound=LIMITS_HI                &  !<-- Ending index in each dimension
                                  ,rc         =RC )
!
                DO L=LIMITS_LO(3),LIMITS_HI(3)
                  DO J=JRECV_START,JRECV_END
                  DO I=IRECV_START,IRECV_END
                    KOUNT_REAL=KOUNT_REAL+1
                    ARRAY_3D(I,J,L)=RECV_REAL_DATA(KOUNT_REAL)             !<-- Task updates 3-D Real H data in overlap region
                  ENDDO
                  ENDDO
                ENDDO
!
              ENDIF
!
            ENDDO
!
!--------------
!***  V points
!--------------
!
            NUM_FIELDS=NUM_FIELDS_2D_V+NUM_FIELDS_3D_V
!
            DO N_FIELD=1,NUM_FIELDS
!
              CALL ESMF_FieldBundleGet(FIELDBUNDLE=MOVE_BUNDLE_V        &  !<-- Bundle holding the V arrays for move updates
                                      ,fieldIndex =N_FIELD              &  !<-- Index of the Field in the Bundle
                                      ,field      =HOLD_FIELD           &  !<-- Field N_FIELD in the Bundle
                                      ,rc         =RC )
!
              CALL ESMF_FieldGet(field   =HOLD_FIELD                    &  !<-- Field N_FIELD in the Bundle
                                ,dimCount=NUM_DIMS                      &  !<-- Is this Field 2-D or 3-D?
                                ,name      =FIELD_NAME                  &  !<-- This Field's name
                                ,rc      =RC )
! 
              N_REMOVE=INDEX(FIELD_NAME,SUFFIX_MOVE)
              FIELD_NAME=FIELD_NAME(1:N_REMOVE-1)
!
              IF(NUM_DIMS==2)THEN
!
                CALL ESMF_FieldGet(field    =HOLD_FIELD                 &  !<-- Field N_FIELD in the Bundle
                                  ,localDe  =0                          &
                                  ,farrayPtr=ARRAY_2D                   &  !<-- Dummy 2-D array with the Field's data
                                  ,rc       =RC )
!
                DO J=JRECV_START,JRECV_END
                DO I=IRECV_START,IRECV_END
                  KOUNT_REAL=KOUNT_REAL+1
                  ARRAY_2D(I,J)=RECV_REAL_DATA(KOUNT_REAL)                 !<-- Task updates (REAL) 2-D V data in overlap region
                ENDDO
                ENDDO
!
              ELSEIF(NUM_DIMS==3)THEN
!
                CALL ESMF_FieldGet(field      =HOLD_FIELD               &  !<-- Field N_FIELD in the Bundle
                                  ,localDe    =0                        &
                                  ,farrayPtr  =ARRAY_3D                 &  !<-- Dummy 3-D array with the Field's data
                                  ,totalLBound=LIMITS_LO                &  !<-- Starting index in each dimension
                                  ,totalUBound=LIMITS_HI                &  !<-- Ending index in each dimension
                                  ,rc         =RC )
!
                DO L=LIMITS_LO(3),LIMITS_HI(3)
                  DO J=JRECV_START,JRECV_END
                  DO I=IRECV_START,IRECV_END
                    KOUNT_REAL=KOUNT_REAL+1
                    ARRAY_3D(I,J,L)=RECV_REAL_DATA(KOUNT_REAL)             !<-- Task updates (Real) 3-D V data in overlap region
                  ENDDO
                  ENDDO
                ENDDO
!
              ENDIF
!
            ENDDO
!
!-----------------------------------------------------------------------
!
          ENDIF recving
!
!-----------------------------------------------------------------------
!
        ENDIF check
!
!-----------------------------------------------------------------------
!
      ENDDO recv_loop
!
!-----------------------------------------------------------------------
!
      IF(ALLOCATED(RECV_REAL_DATA))DEALLOCATE(RECV_REAL_DATA)
      IF(ALLOCATED(RECV_INTEGER_DATA))DEALLOCATE(RECV_INTEGER_DATA)
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE RECV_INTER_TASK_DATA
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE UPDATE_INTERIOR_FROM_PARENT(IMP_STATE                  &
                                            ,SFC_FILE_RATIO             &
                                            ,MOVE_BUNDLE_H              &
                                            ,NUM_FIELDS_2D_H_I          &
                                            ,NUM_FIELDS_2D_H_R          &
                                            ,NUM_FIELDS_3D_H            &
                                            ,MOVE_BUNDLE_V              &
                                            ,NUM_FIELDS_2D_V            &
                                            ,NUM_FIELDS_3D_V            &
                                            ,GLAT_H                     &
                                            ,GLON_H                     &
                                            ,ITS,ITE,JTS,JTE            &
                                            ,IMS,IME,JMS,JME )
!
!-----------------------------------------------------------------------
!***  After the nest has moved update all nest gridpoints in that
!***  domain's interior that lie outside of the footprint of its
!***  pre-move location.  The update data comes from the parent.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      INTEGER(kind=KINT),INTENT(IN) :: ITS,ITE,JTS,JTE                  &  !<-- Integration limits 
                                      ,IMS,IME,JMS,JME                     !<-- Memory limits 
!
      INTEGER(kind=KINT),INTENT(IN) :: NUM_FIELDS_2D_H_I                &  !<-- # of 2-D integer H variables to update
                                      ,NUM_FIELDS_2D_H_R                &  !<-- # of 2-D real H variables to update
                                      ,NUM_FIELDS_3D_H                  &  !<-- # of 3-D H variables to update
                                      ,NUM_FIELDS_2D_V                  &  !<-- # of 2-D V variables to update
                                      ,NUM_FIELDS_3D_V                     !<-- # of 3-D V variables to update
!
      INTEGER(kind=KINT),INTENT(IN) :: SFC_FILE_RATIO                      !<-- Ratio of upper parent grid increment to this domain's
!
      REAL(kind=KFPT),DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: GLAT_H   &  !<-- This domain's geographic latitude (radians) on H pts
                                                              ,GLON_H      !<-- This domain's geographic longitude (radians) on H pts
!
      TYPE(ESMF_State),INTENT(INOUT) :: IMP_STATE                          !<-- The Domain import state
!
      TYPE(ESMF_FieldBundle),INTENT(INOUT) :: MOVE_BUNDLE_H             &  !<-- Bundle of internal state H arrays needing updates
                                             ,MOVE_BUNDLE_V                !<-- Bundle of internal state V arrays needing updates
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: I,I_END,I_OFFSET,I_START                    &
                           ,IHI,ILO,INPUT_NEST                          &
                           ,J,J_END,J_OFFSET,J_START                    &
                           ,JHI,JLO,KHI,KLO                             &
                           ,KOUNT_INTEGER,KOUNT_REAL                    &
                           ,N,N_FIELD,N_ITER,N_REMOVE,NI,NL,NN          &
                           ,NUM_DIMS,NUM_FIELDS                         &
                           ,NUM_INTEGER_WORDS                           &
                           ,NUM_PTASK_UPDATE,NUM_REAL_WORDS             &
                           ,UPDATE_TYPE_INT
!
      INTEGER(kind=KINT) :: I_COUNT_DATA,I_START_DATA                   &
                           ,J_COUNT_DATA,J_START_DATA                   &
                           ,NCID,NCTYPE,NDIMS,NX,NY,VAR_ID
!
      INTEGER(kind=KINT) :: N_FIELD_T,N_FIELD_TP                        &
                           ,N_FIELD_U,N_FIELD_UP                        &
                           ,N_FIELD_V,N_FIELD_VP
!
      INTEGER(kind=KINT) :: IERR,RC,RC_UPDATE
!
      INTEGER(kind=KINT),DIMENSION(1:2) :: DIM_IDS,LBND,UBND
!
      INTEGER(kind=KINT),DIMENSION(1:3) :: LIMITS_HI                    &
                                          ,LIMITS_LO
!
      INTEGER(kind=KINT),DIMENSION(1:8) :: INDICES_H,INDICES_V
!
      INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE :: UPDATE_INTEGER_DATA
!
      INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE :: I_INDX,J_INDX
!
      INTEGER(kind=KINT),DIMENSION(:,:),ALLOCATABLE :: SFC_IDATA
!
      INTEGER(kind=KINT),DIMENSION(:,:),POINTER :: IARRAY_2D=>NULL()
!
      REAL(kind=KFPT) :: GBL,REAL_I,REAL_J
!
      REAL(kind=KFPT),DIMENSION(:),ALLOCATABLE :: UPDATE_REAL_DATA
!
      REAL(kind=KFPT),DIMENSION(:,:),POINTER :: ARRAY_2D=>NULL()        &
!
                                               ,ALBASE=>NULL()          &
                                               ,SEA_MASK=>NULL()        &
                                               ,SSTX=>NULL()
!
      REAL(kind=KFPT),DIMENSION(:,:,:),POINTER :: ARRAY_3D=>NULL()      &
                                                 ,ARRAY_3D_X=>NULL()
!
      CHARACTER(len=1)  :: N_PTASK                                      &
                          ,UPDATE_TYPE_CHAR
      CHARACTER(len=2)  :: ID_SFC_FILE
      CHARACTER(len=12) :: NAME
      CHARACTER(len=15) :: VNAME
      CHARACTER(len=17) :: NAME_REAL
      CHARACTER(len=20) :: NAME_INTEGER
      CHARACTER(len=99) :: FIELD_NAME                                   &
                          ,FILENAME
!
      LOGICAL(kind=KLOG) :: OPENED
!
      TYPE(ESMF_Field) :: HOLD_FIELD                                    &
                         ,HOLD_FIELD_X
!
      TYPE(ESMF_TypeKind_Flag) :: DATATYPE
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Unload the number of parent tasks who provide update data
!***  for this nest task.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="How Many Parent Tasks Sent Interior Updates?"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=IMP_STATE                            &  !<-- The Domain import state
                            ,name ='Num Parent Tasks Update'            &  !<-- Name of the variable
                            ,value=NUM_PTASK_UPDATE                     &  !<-- # of parent tasks that update this nest task
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_UPDATE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  If no parent tasks provide data then there is nothing to do
!***  so RETURN.
!-----------------------------------------------------------------------
!
      IF(NUM_PTASK_UPDATE==0)RETURN
!
!-----------------------------------------------------------------------
!***  Unload each piece of data that was sent by each parent task
!***  and apply it to given locations within the arrays in the
!***  Move Bundles for H-pt and V-pt variables.
!-----------------------------------------------------------------------
!
      parent_loop: DO N=1,NUM_PTASK_UPDATE
!
!-----------------------------------------------------------------------
!
        KOUNT_INTEGER=0                                                    !<-- Count the integer update points from this parent task
        KOUNT_REAL=0                                                       !<-- Count the real update points from this parent task
!
        WRITE(N_PTASK,'(I1)')N
        NAME_INTEGER='PTASK_INTEGER_DATA_'//N_PTASK
        NAME_REAL   ='PTASK_REAL_DATA_'//N_PTASK
        NAME        ='PTASK_DATA_'//N_PTASK
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Unload # of Words in Integer Update Data from Parent"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state=IMP_STATE                          &  !<-- The Domain import state
                              ,name =NAME_INTEGER//' Words'             &  !<-- Name of the variable
                              ,value=NUM_INTEGER_WORDS                  &  !<-- # of words in integer update data from Nth parent task
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_UPDATE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
        unload_int: IF(NUM_INTEGER_WORDS>0)THEN
!
!-----------------------------------------------------------------------
!
          ALLOCATE(UPDATE_INTEGER_DATA(1:NUM_INTEGER_WORDS))
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Unload Interior Update Integer Data from Input State"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeGet(state    =IMP_STATE                    &  !<-- The Domain import state
                                ,name     =NAME_INTEGER                 &  !<-- Name of the variable
                                ,itemCount=NUM_INTEGER_WORDS            &  !<-- # of integer words in update data from Nth parent task
                                ,valueList=UPDATE_INTEGER_DATA          &  !<-- The integer update data from Nth parent task
                                ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_UPDATE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
        ENDIF unload_int
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Unload # of Words in Real Update Data from Parent"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state=IMP_STATE                          &  !<-- The Domain import state
                              ,name =NAME_REAL//' Words'                &  !<-- Name of the variable
                              ,value=NUM_REAL_WORDS                     &  !<-- # of words in real update data from Nth parent task
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_UPDATE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        ALLOCATE(UPDATE_REAL_DATA(1:NUM_REAL_WORDS))
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Unload Interior Update Real Data from Input State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state    =IMP_STATE                      &  !<-- The Domain import state
                              ,name     =NAME_REAL                      &  !<-- Name of the variable
                              ,itemCount=NUM_REAL_WORDS                 &  !<-- # of real words in update data from Nth parent task
                              ,valueList=UPDATE_REAL_DATA               &  !<-- The real update data from Nth parent task
                              ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_UPDATE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Index Limits for Update Data from Domain Input State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state    =IMP_STATE                      &  !<-- The Domain import state
                              ,name     =NAME//' Indices H'             &  !<-- Name of the variable
                              ,itemCount=N8                             &  !<-- # of words in index limits of update data
                              ,valueList=INDICES_H                      &  !<-- The update data index specifications for H
                              ,rc       =RC)
!
        CALL ESMF_AttributeGet(state    =IMP_STATE                      &  !<-- The Domain import state
                              ,name     =NAME//' Indices V'             &  !<-- Name of the variable
                              ,itemCount=N8                             &  !<-- # of words in index limits of update data
                              ,valueList=INDICES_V                      &  !<-- The update data index specifications for V
                              ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_UPDATE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Update the nest's H points from parent task N.
!-----------------------------------------------------------------------
!
        N_ITER=1
        IF(INDICES_H(2)>-99)THEN                                           !<-- If true, there are 2 update regions from parent task
          N_ITER=2
        ENDIF
!
        NUM_FIELDS=NUM_FIELDS_2D_H_I+NUM_FIELDS_2D_H_R+NUM_FIELDS_3D_H
!
!-----------------------------------------------------------------------
!***  Typically there is only one update region within the task
!***  unless the task lies on a corner of the pre-move footprint
!***  in which case there are two update regions.  N_ITER in the
!***  following loop refers to the number of update regions for 
!***  this nest task.
!-----------------------------------------------------------------------
!
        iterations_h: DO NI=1,N_ITER                                       !<-- Update each region with parent information
!
!-----------------------------------------------------------------------
!
          I_START=INDICES_H(1)
          I_END  =INDICES_H(3)
          J_START=INDICES_H(5)
          J_END  =INDICES_H(7)
!
          IF(NI==2)THEN
            I_START=INDICES_H(2)
            I_END  =INDICES_H(4)
            J_START=INDICES_H(6)
            J_END  =INDICES_H(8)
          ENDIF
!
!-----------------------------------------------------------------------
!***  For those 2-D surface variables that must be read from
!***  external files after nests shift we will need the regions
!***  within each relevant nest task that must be updated.
!-----------------------------------------------------------------------
!
          IF(GLOBAL_TOP_PARENT)THEN
            GBL=1.                                                         !<-- Account for the extra row that surrounds the global domain.
          ELSE
            GBL=0.
          ENDIF
!
          CALL LATLON_TO_IJ(GLAT_H(I_START,J_START)                     &  !<-- Geographic latitude of nest task's 1st update point
                           ,GLON_H(I_START,J_START)                     &  !<-- Geographic longitude of nest task's 1st update point
                           ,TPH0_1,TLM0_1                               &  !<-- Central lat/lon (radians, N/E) of uppermost parent
                           ,SB_1,WB_1                                   &  !<-- Rotated lat/lon of upper parent's S/W bndry (radians, N/E)
                           ,RECIP_DPH_1,RECIP_DLM_1                     &  !<-- Reciprocal of I/J grid increments (radians) on upper parent
                           ,GLOBAL_TOP_PARENT                           &  !<-- Is the uppermost parent on a global grid?
                           ,REAL_I                                      &  !<-- Corresponding I index on uppermost parent grid
                           ,REAL_J)                                        !<-- Corresponding J index on uppermost parent grid
!
          I_OFFSET=NINT((REAL_I-1.-GBL)*SFC_FILE_RATIO)                    !<-- Offset in I between sfc file index and nest index
          J_OFFSET=NINT((REAL_J-1.-GBL)*SFC_FILE_RATIO)                    !<-- Offset in J between sfc file index and nest index
!
          I_START_DATA=I_OFFSET+1                                          !<-- Start reading at this I in the external file array
          I_COUNT_DATA=I_END-I_START+1                                     !<-- Read this many points in I
          J_START_DATA=J_OFFSET+1                                          !<-- Start reading at this J in the external file array
          J_COUNT_DATA=J_END-J_START+1                                     !<-- Read this many points in J
!
!-----------------------------------------------------------------------
!***  Now proceed with the updating of all H-pt variables with data
!***  sent from the parent.  All Fields in the Move Bundles have names
!***  containing the suffix '-move' to distinguish them from the same
!***  Fields that occur in the same ESMF States for different reasons.
!***  The suffix is removed to obtain the actual name.
!-----------------------------------------------------------------------
!
          fields_h: DO N_FIELD=1,NUM_FIELDS                                !<-- Update all H-point arrays
!
!-----------------------------------------------------------------------
!
            CALL ESMF_FieldBundleGet(FIELDBUNDLE=MOVE_BUNDLE_H          &  !<-- Bundle holding the H arrays for move updates
                                    ,fieldIndex =N_FIELD                &  !<-- Index of the Field in the Bundle
                                    ,field      =HOLD_FIELD             &  !<-- Field N_FIELD in the Bundle
                                    ,rc         =RC )
!
            CALL ESMF_FieldGet(field   =HOLD_FIELD                      &  !<-- Field N_FIELD in the Bundle
                              ,dimCount=NUM_DIMS                        &  !<-- Is this Field 2-D or 3-D?
                              ,typekind=DATATYPE                        &  !<-- Is the data integer or real?
                              ,name    =FIELD_NAME                      &  !<-- Name of the Field
                              ,rc      =RC )
!
            N_REMOVE=INDEX(FIELD_NAME,SUFFIX_MOVE)
            FIELD_NAME=FIELD_NAME(1:N_REMOVE-1)
!
!-----------------------------------------------------------------------
!***                    ***** TEMPORARY *****
!***                    ****  T/TP, etc. ****
!-----------------------------------------------------------------------
!***  The variables PDO, TP, UP, and VP are valid one timestep before
!***  the current one.  Since the timestep of the parent is larger than
!***  that of its nests then if the parent simply interpolated those
!***  variables spatially to the nest update points then they would
!***  not be valid at one timestep before the nest's current time.
!***  Fixing this additional temporal interpolation will be done
!***  at a later date therefore for the moment the parent will send
!***  its spatially interpolated values of PD, T, U, and V for those
!***  four variables.  The parent already substituted its interpolated
!***  PD values for PDO in subroutine PARENT_UPDATES_MOVING.  Now save
!***  the locations of T and TP in the H Move Bundle so the analogous
!***  substitution can be made after the fields_h loop immediately
!***  following the updates of all H-point variables.  The same is
!***  done in the field_v loop for U,UP and V,VP.
!-----------------------------------------------------------------------
!
            IF(FIELD_NAME=='T')THEN
              N_FIELD_T=N_FIELD
            ELSEIF(FIELD_NAME=='TP')THEN
              N_FIELD_TP=N_FIELD
            ENDIF
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Extract UPDATE_TYPE from Field"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_AttributeGet(field=HOLD_FIELD                     &  !<-- Get Attribute from this Field
                                  ,name ='UPDATE_TYPE'                  &  !<-- Name of the attribute to extract
                                  ,value=UPDATE_TYPE_INT                &  !<-- Value of the Attribute
                                  ,rc   =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_UPDATE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            IF(UPDATE_TYPE_INT==1)THEN
              UPDATE_TYPE_CHAR='H'                                         !<-- Ordinary H-pt variable
            ELSEIF(UPDATE_TYPE_INT==2)THEN
              UPDATE_TYPE_CHAR='L'                                         !<-- H-pt land surface variable
            ELSEIF(UPDATE_TYPE_INT==3)THEN
              UPDATE_TYPE_CHAR='W'                                         !<-- H-pt water surface variable
            ELSEIF(UPDATE_TYPE_INT==4)THEN
              UPDATE_TYPE_CHAR='F'                                         !<-- H-pt variable obtained from an external file
            ELSEIF(UPDATE_TYPE_INT==5)THEN
              UPDATE_TYPE_CHAR='V'                                         !<-- Ordinary V-pt variable
            ENDIF
!
!-----------------------------------------------------------------------
!***  Updated 2-D H-point arrays include both Integer and Real.
!-----------------------------------------------------------------------
!
            IF(NUM_DIMS==2)THEN
!
!-----------------------------------------------------------------------
!***  Some surface-related variables that do not change in time are
!***  read directly by the relevant moving nest tasks from external
!***  files.  Those files contain the data at the nest's resolution
!***  but span the entire uppermost parent domain.  Such variables
!***  are all 2-D.
!-----------------------------------------------------------------------
!
              update_type: IF(UPDATE_TYPE_CHAR=='F')THEN                   !<-- If so, the variable is updated from an external file
!
!-----------------------------------------------------------------------
!
                IF(SFC_FILE_RATIO<=9)THEN
                  WRITE(ID_SFC_FILE,'(I1.1)')SFC_FILE_RATIO
                ELSEIF(SFC_FILE_RATIO>=10)THEN
                  WRITE(ID_SFC_FILE,'(I2.2)')SFC_FILE_RATIO
                ENDIF
!
                FILENAME=TRIM(FIELD_NAME)//'_'//TRIM(ID_SFC_FILE)//'.nc'
!
                CALL CHECK(NF90_OPEN(FILENAME,NF90_NOWRITE,NCID))          !<-- Open the current field's external netCDF file.
!
!-----------------------------------------------------------------------
!
!----------
!***  Real
!----------
!
                IF(DATATYPE==ESMF_TYPEKIND_R4)THEN                         !<-- The 2-D H-point external file data is Real
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
                  MESSAGE_CHECK="Extract 2-D Real Array for Type F"
!                 CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
                  CALL ESMF_FieldGet(field    =HOLD_FIELD               &  !<-- Field N_FIELD in the Bundle
                                    ,localDe  =0                        &
                                    ,farrayPtr=ARRAY_2D                 &  !<-- Dummy 2-D array with Field's Real data
                                    ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
                  CALL ERR_MSG(RC,MESSAGE_CHECK,RC_UPDATE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
                  CALL CHECK(NF90_INQUIRE_VARIABLE(NCID,3,VNAME,NCTYPE  &
                                                  ,NDIMS,DIM_IDS))
                  CALL CHECK(NF90_INQ_VARID(NCID,VNAME,VAR_ID))
!
                  CALL CHECK(NF90_GET_VAR(NCID,VAR_ID                           &  !<-- Extract the desired real values from the
                                         ,ARRAY_2D(I_START:I_END,J_START:J_END) &  !    current field's external file.
                                         ,start=(/I_START_DATA,J_START_DATA/)   &  !    Nest points that have moved beyond the
                                         ,count=(/I_COUNT_DATA,J_COUNT_DATA/)))    !    pre-move footprint are updated.
!
!-----------------------------------------------------------------------
!***  Save the base albedo to provide values for the dynamic albedo
!***  in parent update regions when needed.
!-----------------------------------------------------------------------
!
!                 IF(TRIM(FIELDNAME)=='ALBASE')THEN
!                   ALBASE=>ARRAY_2D                                       !<-- Save the base albedo for fixing conflict points 
!                 ENDIF
!
!-----------------------------------------------------------------------
!***  Save the sea mask for use in cleaning up surface variables 
!***  following this primary update.  This needs to be done when
!***  the parent uses a land(water) point to update a nest water(land)
!***  point.
!-----------------------------------------------------------------------
!
                  IF(TRIM(FIELD_NAME)=='SM')THEN
                    SEA_MASK=>ARRAY_2D 
                  ENDIF
!
!-------------
!***  Integer
!-------------
!
                ELSEIF(DATATYPE==ESMF_TYPEKIND_I4)THEN                     !<-- The 2-D H-point external file data is Integer
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
                  MESSAGE_CHECK="Extract 2-D Integer Array for Type F"
!                 CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
                  CALL ESMF_FieldGet(field    =HOLD_FIELD               &  !<-- Field N_FIELD in the Bundle
                                    ,localDe  =0                        &
                                    ,farrayPtr=IARRAY_2D                &  !<-- Dummy 2-D array with Field's Real data
                                    ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
                  CALL ERR_MSG(RC,MESSAGE_CHECK,RC_UPDATE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
                  CALL CHECK(NF90_INQUIRE_VARIABLE(NCID,3,VNAME,NCTYPE  &
                                                  ,NDIMS,DIM_IDS))
                  CALL CHECK(NF90_INQ_VARID(NCID,VNAME,VAR_ID))
!
                  CALL CHECK(NF90_GET_VAR(NCID,VAR_ID                            &  !<-- Extract the desired integer values from the
                                         ,IARRAY_2D(I_START:I_END,J_START:J_END) &  !    current field's external file.
                                         ,start=(/I_START_DATA,J_START_DATA/)    &  !    Nest points that have moved beyond the
                                         ,count=(/I_COUNT_DATA,J_COUNT_DATA/)))     !    pre-move footprint are updated.
!
                ENDIF
!
                CALL CHECK(NF90_CLOSE(NCID))                               !<-- Close the external netCDF file.
!
                CYCLE fields_h
!
!-----------------------------------------------------------------------
!***  All other variables that are not specified with UPDATE_TYPE='F'
!***  are updated from the data sent from the parent.
!-----------------------------------------------------------------------
!
              ELSE update_type
!
!-------------
!***  Integer
!-------------
!
                IF(DATATYPE==ESMF_TYPEKIND_I4                           &  !<-- The 2-D H-point array to be updated is Integer
                             .AND.                                      &
                   NUM_INTEGER_WORDS>0)THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
                  MESSAGE_CHECK="Extract General 2-D Integer Array"
!                 CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
                  CALL ESMF_FieldGet(field    =HOLD_FIELD               &  !<-- Field N_FIELD in the Bundle
                                    ,localDe  =0                        &
                                    ,farrayPtr=IARRAY_2D                &  !<-- Dummy 2-D array with Field's integer data
                                    ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
                  CALL ERR_MSG(RC,MESSAGE_CHECK,RC_UPDATE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  This nest task incorporates the standard 2-D Integer update data
!***  sent to it from the parent.
!-----------------------------------------------------------------------
!
                  DO J=J_START,J_END 
                  DO I=I_START,I_END
                    KOUNT_INTEGER=KOUNT_INTEGER+1
                    IARRAY_2D(I,J)=UPDATE_INTEGER_DATA(KOUNT_INTEGER)
                  ENDDO
                  ENDDO
!
!----------
!***  Real
!----------
!
                ELSEIF(DATATYPE==ESMF_TYPEKIND_R4)THEN                     !<-- The 2-D H-point array to be updated is Real
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
                  MESSAGE_CHECK="Extract General 2-D Real Array"
!                 CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
                  CALL ESMF_FieldGet(field    =HOLD_FIELD               &  !<-- Field N_FIELD in the Bundle
                                    ,localDe  =0                        &
                                    ,farrayPtr=ARRAY_2D                 &  !<-- Dummy 2-D array with Field's Real data
                                    ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
                  CALL ERR_MSG(RC,MESSAGE_CHECK,RC_UPDATE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  This nest task incorporates the standard 2-D Real update
!***  data sent to it from the parent.
!-----------------------------------------------------------------------
!
                  DO J=J_START,J_END 
                  DO I=I_START,I_END
                    KOUNT_REAL=KOUNT_REAL+1
                    ARRAY_2D(I,J)=UPDATE_REAL_DATA(KOUNT_REAL)
                  ENDDO
                  ENDDO
!
!-----------------------------------------------------------------------
!***  Save the sea surface temperature as a flag for discriminating
!***  between sea and land points when there is a conflict between
!***  the nest's sea mask and the type of data (sea or land) it is
!***  sent by its parent.  The SST should always equal 0.0 at land
!***  points.
!-----------------------------------------------------------------------
!
!                 IF(TRIM(FIELD_NAME)=='SST')THEN   
!                   SSTX=>ARRAY_2D
!                 ENDIF
!
                ENDIF
!
!-----------------------------------------------------------------------
!
              ENDIF update_type  
!
!-----------------------------------------------------------------------
!
            ELSEIF(NUM_DIMS==3)THEN                                        !<-- The parent's update of all 3-D variables
!
!-----------------------------------------------------------------------
!
              CALL ESMF_FieldGet(field      =HOLD_FIELD                 &  !<-- Field N in the Bundle
                                ,localDe    =0                          &
                                ,farrayPtr  =ARRAY_3D                   &  !<-- Dummy 3-D array with Field's data
                                ,totalLBound=LIMITS_LO                  &  !<-- Starting index in each dimension
                                ,totalUBound=LIMITS_HI                  &  !<-- Ending index in each dimension
                                ,rc         =RC )
!
              DO NL=LIMITS_LO(3),LIMITS_HI(3)
                DO J=J_START,J_END
                DO I=I_START,I_END
                  KOUNT_REAL=KOUNT_REAL+1
                  ARRAY_3D(I,J,NL)=UPDATE_REAL_DATA(KOUNT_REAL)
                ENDDO
                ENDDO
              ENDDO
!
            ENDIF
!
!-----------------------------------------------------------------------
!
          ENDDO fields_h
!
!-----------------------------------------------------------------------
!***  The temporary substitution of T into TP.  See note above ('T/TP').
!-----------------------------------------------------------------------
!
          CALL ESMF_FieldBundleGet(FIELDBUNDLE=MOVE_BUNDLE_H            &  !<-- Bundle holding the H arrays for move updates
                                  , fieldIndex=N_FIELD_T                &  !<-- Index of the T Field in the H Bundle
                                  ,field      =HOLD_FIELD               &  !<-- The T Field in the H Bundle
                                  ,rc         =RC )
!
          CALL ESMF_FieldBundleGet(FIELDBUNDLE=MOVE_BUNDLE_H            &  !<-- Bundle holding the H arrays for move updates
                                  ,fieldIndex =N_FIELD_TP               &  !<-- Index of the TP Field in the H Bundle
                                  ,field      =HOLD_FIELD_X             &  !<-- The TP Field in the H Bundle
                                  ,rc         =RC )
!
          CALL ESMF_FieldGet(field      =HOLD_FIELD                     &  !<-- The T Field in the H Bundle
                            ,localDe    =0                              &
                            ,farrayPtr  =ARRAY_3D                       &  !<-- Dummy 3-D array with Field's data
                            ,totalLBound=LIMITS_LO                      &  !<-- Starting index in each dimension
                            ,totalUBound=LIMITS_HI                      &  !<-- Ending index in each dimension
                            ,rc         =RC )
          CALL ESMF_FieldGet(field      =HOLD_FIELD_X                   &  !<-- The TP Field in the H Bundle
                            ,localDe    =0                              &
                            ,farrayPtr  =ARRAY_3D_X                     &  !<-- Dummy 3-D array with Field's data
                            ,totalLBound=LIMITS_LO                      &  !<-- Starting index in each dimension
                            ,totalUBound=LIMITS_HI                      &  !<-- Ending index in each dimension
                            ,rc         =RC )
!
          DO NL=LIMITS_LO(3),LIMITS_HI(3)
            DO J=J_START,J_END
            DO I=I_START,I_END
              ARRAY_3D_X(I,J,NL)=ARRAY_3D(I,J,NL)                          !<-- For now fill TP with T update values
            ENDDO
            ENDDO  
          ENDDO
!
!-----------------------------------------------------------------------
!***  Now go back and correct any mismatches between the nest's
!***  sea/land mask and the type of data (land or water) received
!***  from the parent.
!-----------------------------------------------------------------------
!
!
          ILO=LBOUND(SEA_MASK,1)
          IHI=UBOUND(SEA_MASK,1)
          JLO=LBOUND(SEA_MASK,2)
          JHI=UBOUND(SEA_MASK,2)
!
          CALL FIX_SFC(MOVE_BUNDLE_H,NUM_FIELDS                         &
                      ,SEA_MASK                                         &
                      ,ILO,IHI,JLO,JHI                                  &
                      ,I_START,I_END,J_START,J_END)
!
!-----------------------------------------------------------------------
!
        ENDDO iterations_h
!
!-----------------------------------------------------------------------
!***  Update the nest's V points from parent task N.
!-----------------------------------------------------------------------
!
        N_ITER=1
        IF(INDICES_V(2)>-99)THEN                                             !<-- If true, there are 2 update regions from parent task
          N_ITER=2
        ENDIF
!
        NUM_FIELDS=NUM_FIELDS_2D_V+NUM_FIELDS_3D_V
!
!-----------------------------------------------------------------------
!
        iterations_v: DO NI=1,N_ITER
!
!-----------------------------------------------------------------------
!
          I_START=INDICES_V(1)
          I_END  =INDICES_V(3)
          J_START=INDICES_V(5)
          J_END  =INDICES_V(7)
!
          IF(NI==2)THEN
            I_START=INDICES_V(2)
            I_END  =INDICES_V(4)
            J_START=INDICES_V(6)
            J_END  =INDICES_V(8)
          ENDIF
!
!-----------------------------------------------------------------------
!
          fields_v: DO N_FIELD=1,NUM_FIELDS
!
!-----------------------------------------------------------------------
!
            CALL ESMF_FieldBundleGet(FIELDBUNDLE=MOVE_BUNDLE_V          &  !<-- Bundle holding the V arrays for move updates
                                    ,fieldIndex =N_FIELD                &  !<-- Index of the Field in the Bundle
                                    ,field      =HOLD_FIELD             &  !<-- Field N_FIELD in the Bundle
                                    ,rc         =RC )
!
            CALL ESMF_FieldGet(field   =HOLD_FIELD                      &  !<-- Field N_FIELD in the Bundle
                              ,dimCount=NUM_DIMS                        &  !<-- Is this Field 2-D or 3-D?
                              ,name      =FIELD_NAME                    &  !<-- Name of the Field
                              ,rc      =RC )
!
            N_REMOVE=INDEX(FIELD_NAME,SUFFIX_MOVE)
            FIELD_NAME=FIELD_NAME(1:N_REMOVE-1)
!
!-----------------------------------------------------------------------
!***                     ***** TEMPORARY *****
!***  See above regarding T/TP, etc.
!-----------------------------------------------------------------------
!
            IF(FIELD_NAME=='U')THEN
              N_FIELD_U=N_FIELD
            ELSEIF(FIELD_NAME=='UP')THEN
              N_FIELD_UP=N_FIELD
            ELSEIF(FIELD_NAME=='V')THEN
              N_FIELD_V=N_FIELD
            ELSEIF(FIELD_NAME=='VP')THEN
              N_FIELD_VP=N_FIELD
            ENDIF
!
!-----------------------------------------------------------------------
!
            IF(NUM_DIMS==2)THEN
!
!-----------------------------------------------------------------------
!
              CALL ESMF_FieldGet(field    =HOLD_FIELD                   &  !<-- Field N_FIELD in the Bundle
                                ,localDe  =0                            &
                                ,farrayPtr=ARRAY_2D                     &  !<-- Dummy 2-D array with Field's data
                                ,rc       =RC )
!
!-----------------------------------------------------------------------
!***  Update this 2-D array with values from the parent.
!-----------------------------------------------------------------------
!
              DO J=J_START,J_END 
              DO I=I_START,I_END
                KOUNT_REAL=KOUNT_REAL+1
                ARRAY_2D(I,J)=UPDATE_REAL_DATA(KOUNT_REAL)
              ENDDO
              ENDDO
!
!-----------------------------------------------------------------------
!
            ELSEIF(NUM_DIMS==3)THEN
!
!-----------------------------------------------------------------------
!
              CALL ESMF_FieldGet(field      =HOLD_FIELD                 &  !<-- Field N in the Bundle
                                ,localDe    =0                          &
                                ,farrayPtr  =ARRAY_3D                   &  !<-- Dummy 3-D array with Field's data
                                ,totalLBound=LIMITS_LO                  &  !<-- Starting index in each dimension
                                ,totalUBound=LIMITS_HI                  &  !<-- Ending index in each dimension
                                ,rc         =RC )
!
              DO NL=LIMITS_LO(3),LIMITS_HI(3)
                DO J=J_START,J_END
                DO I=I_START,I_END
                  KOUNT_REAL=KOUNT_REAL+1
                  ARRAY_3D(I,J,NL)=UPDATE_REAL_DATA(KOUNT_REAL)
                ENDDO
                ENDDO
              ENDDO
!
            ENDIF
!
!-----------------------------------------------------------------------
!
          ENDDO fields_v
!
!-----------------------------------------------------------------------
!***  The temporary substitution of U/V into UP/VP.  See note above.
!-----------------------------------------------------------------------
!
!--------
!***  UP
!--------
!
          CALL ESMF_FieldBundleGet(FIELDBUNDLE=MOVE_BUNDLE_V            &  !<-- Bundle holding the V arrays for move updates
                                  ,fieldIndex =N_FIELD_U                &  !<-- Index of the U Field in the V Bundle
                                  ,field      =HOLD_FIELD               &  !<-- The U Field in the V Bundle
                                  ,rc         =RC )
!
          CALL ESMF_FieldBundleGet(FIELDBUNDLE=MOVE_BUNDLE_V            &  !<-- Bundle holding the V arrays for move updates
                                  ,fieldIndex =N_FIELD_UP               &  !<-- Index of the UP Field in the V Bundle
                                  ,field      =HOLD_FIELD_X             &  !<-- The UP Field in the V Bundle
                                  ,rc         =RC )
!
          CALL ESMF_FieldGet(field      =HOLD_FIELD                     &  !<-- The U Field in the V Bundle
                            ,localDe    =0                              &
                            ,farrayPtr  =ARRAY_3D                       &  !<-- Dummy 3-D array with Field's data
                            ,totalLBound=LIMITS_LO                      &  !<-- Starting index in each dimension
                            ,totalUBound=LIMITS_HI                      &  !<-- Ending index in each dimension
                            ,rc         =RC )
          CALL ESMF_FieldGet(field      =HOLD_FIELD_X                   &  !<-- The UP Field in the V Bundle
                            ,localDe    =0                              &
                            ,farrayPtr  =ARRAY_3D_X                     &  !<-- Dummy 3-D array with Field's data
                            ,totalLBound=LIMITS_LO                      &  !<-- Starting index in each dimension
                            ,totalUBound=LIMITS_HI                      &  !<-- Ending index in each dimension
                            ,rc         =RC )
!
          DO NL=LIMITS_LO(3),LIMITS_HI(3)
            DO J=J_START,J_END
            DO I=I_START,I_END
              ARRAY_3D_X(I,J,NL)=ARRAY_3D(I,J,NL)                          !<-- For now fill UP with U update values
            ENDDO
            ENDDO  
          ENDDO
!
!--------
!***  VP
!--------
!
          CALL ESMF_FieldBundleGet(FIELDBUNDLE=MOVE_BUNDLE_V            &  !<-- Bundle holding the V arrays for move updates
                                  ,fieldIndex =N_FIELD_V                &  !<-- Index of the V Field in the V Bundle
                                  ,field      =HOLD_FIELD               &  !<-- The U Field in the V Bundle
                                  ,rc         =RC )
!
          CALL ESMF_FieldBundleGet(FIELDBUNDLE=MOVE_BUNDLE_V            &  !<-- Bundle holding the V arrays for move updates
                                  ,fieldIndex =N_FIELD_VP               &  !<-- Index of the VP Field in the V Bundle
                                  ,field      =HOLD_FIELD_X             &  !<-- The UP Field in the V Bundle
                                  ,rc         =RC )
!
          CALL ESMF_FieldGet(field      =HOLD_FIELD                     &  !<-- The V Field in the V Bundle
                            ,localDe    =0                              &
                            ,farrayPtr  =ARRAY_3D                       &  !<-- Dummy 3-D array with Field's data
                            ,totalLBound=LIMITS_LO                      &  !<-- Starting index in each dimension
                            ,totalUBound=LIMITS_HI                      &  !<-- Ending index in each dimension
                            ,rc         =RC )
          CALL ESMF_FieldGet(field      =HOLD_FIELD_X                   &  !<-- The VP Field in the V Bundle
                            ,localDe    =0                              &
                            ,farrayPtr  =ARRAY_3D_X                     &  !<-- Dummy 3-D array with Field's data
                            ,totalLBound=LIMITS_LO                      &  !<-- Starting index in each dimension
                            ,totalUBound=LIMITS_HI                      &  !<-- Ending index in each dimension
                            ,rc         =RC )
!
          DO NL=LIMITS_LO(3),LIMITS_HI(3)
            DO J=J_START,J_END
            DO I=I_START,I_END
              ARRAY_3D_X(I,J,NL)=ARRAY_3D(I,J,NL)                          !<-- For now fill VP with V update values
            ENDDO
            ENDDO  
          ENDDO
!
!-----------------------------------------------------------------------
!
        ENDDO iterations_v
!
!-----------------------------------------------------------------------
!
        IF(ALLOCATED(UPDATE_INTEGER_DATA))THEN
          DEALLOCATE(UPDATE_INTEGER_DATA)
        ENDIF
!
        DEALLOCATE(UPDATE_REAL_DATA)
!
!-----------------------------------------------------------------------
!
      ENDDO parent_loop
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE UPDATE_INTERIOR_FROM_PARENT
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE UPDATE_LATLON(IMP_STATE                                &
                              ,EXP_STATE_SOLVER                         &
                              ,I_LBND,I_UBND                            &
                              ,J_LBND,J_UBND                            &
                              ,I_START,I_END                            &
                              ,J_START,J_END                            &
                              ,GLAT_X                                   &
                              ,GLON_X                                   &
                              ,HDACX                                    &
                              ,HDACY                                    &
                              ,F )
!
!-----------------------------------------------------------------------
!***  After the nest has moved recompute the geographic latitude and
!***  longitude on only those points that lie in the parent update
!***  region.  
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      INTEGER(kind=KINT),INTENT(IN) :: I_START,I_END                    &  !<-- Compute lat/lons over this range of I
                                      ,J_START,J_END                       !<-- Compute lat/lons over this range of J
!
      INTEGER(kind=KINT),INTENT(IN) :: I_LBND,I_UBND                    &  !<-- Lower/upper bounds of I in lat/lon arrays
                                      ,J_LBND,J_UBND                       !<-- Lower/upper bounds of J in lat/lon arrays
!
      TYPE(ESMF_State),INTENT(IN) :: IMP_STATE                             !<-- The Domain import state
!
      TYPE(ESMF_State),INTENT(IN) :: EXP_STATE_SOLVER                      !<-- The Solver export state
!
      REAL(kind=KFPT),DIMENSION(I_LBND:I_UBND,J_LBND:J_UBND)            &
                                              ,INTENT(INOUT) :: GLAT_X  &  !<-- Geographic latitude on nest (radians)
                                                               ,GLON_X     !<-- Geographic longitude on nest (radians)
!
      REAL(kind=KFPT),DIMENSION(I_LBND:I_UBND,J_LBND:J_UBND)            &
                                                ,INTENT(OUT) :: HDACX   &  !<-- Lateral diffusion coefficients
                                                               ,HDACY      !<-- 
!
      REAL(kind=KFPT),DIMENSION(I_LBND:I_UBND,J_LBND:J_UBND)            &
                                       ,INTENT(OUT),OPTIONAL :: F
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: I,I_SHIFT,INC_LAT,INC_LON                   &
                           ,J,J_SHIFT                                   &
                           ,KOUNT
!
      INTEGER(kind=KINT) :: RC,RC_RUN
!
      REAL(kind=KFPT) :: A_DLM,ADD,ARG1,ARG2,ARG3                       &
                        ,COS_TPH,COS_TPH0,DY                            &
                        ,SB_PARENT1                                     &
                        ,SIN_TPH,SIN_TPH0,TAN_TPH0                      &
                        ,WB_PARENT1
!
      REAL(kind=KFPT),DIMENSION(J_LBND:J_UBND) :: DX
!
      REAL(kind=KFPT),DIMENSION(I_LBND:I_UBND                           &
                               ,J_LBND:J_UBND) :: TLAT_X,TLON_X
!
      LOGICAL(kind=KLOG) :: VELOCITY                                       !<-- Are we computing for the V points?
!
      TYPE(ESMF_Field) :: FIELD_X
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      IF(.NOT.PRESENT(F))THEN
        VELOCITY=.FALSE.                                                   !<-- H points
        SB_PARENT1=SB_1                                                    !<-- Transformed latitude of upper parent S boundary
        WB_PARENT1=WB_1                                                    !<-- Transformed longitude of upper parent W boundary
      ELSEIF(PRESENT(F))THEN
        VELOCITY=.TRUE.                                                    !<-- V points
        SB_PARENT1=SB_1+0.5*DPH                                            !<-- Transformed latitude of upper parent S boundary
        WB_PARENT1=WB_1+0.5*DLM                                            !<-- Transformed longitude of upper parent W boundary
      ENDIF
!
      DY=A*DPH
      A_DLM=A*DLM
!
      COS_TPH0=COS(TPH0)
      SIN_TPH0=SIN(TPH0)
      TAN_TPH0=TAN(TPH0)
!
!-----------------------------------------------------------------------
!***  Compute the pre-move rotated latitude/longitude of the
!***  update region's SW corner H and V points.  Remember that
!***  GLAT and GLON still have their pre-move values.
!-----------------------------------------------------------------------
!
      CALL GEO_TO_ROT(GLAT_X(I_START,J_START),GLON_X(I_START,J_START)   &
                     ,TLAT_X(I_START,J_START),TLON_X(I_START,J_START))
!
!-----------------------------------------------------------------------
!***  What are the transformed coordinates of the SW corner 
!***  after the nest moved?  And then how many equivalent grid
!***  increments from the upper parent's grid's southern and
!***  western boundary to that SW corner?  By anchoring the
!***  update on the upper parent's domain then the nests will
!***  always generate bit identical results no matter what 
!***  their task layouts are.
!-----------------------------------------------------------------------
!
      TLAT_X(I_START,J_START)=TLAT_X(I_START,J_START)+J_SHIFT_CHILD*DPH
      TLON_X(I_START,J_START)=TLON_X(I_START,J_START)+I_SHIFT_CHILD*DLM
!
      INC_LAT=NINT((TLAT_X(I_START,J_START)-SB_PARENT1)/DPH)
      INC_LON=NINT((TLON_X(I_START,J_START)-WB_PARENT1)/DLM)
!
!-----------------------------------------------------------------------
!***  Now fill in the transformed coordinates on the task subdomain's
!***  update region following the nest's shift.
!-----------------------------------------------------------------------
!
      KOUNT=-1
      DO J=J_START,J_END
        KOUNT=KOUNT+1
        TLAT_X(I_START,J)=SB_PARENT1+(INC_LAT+KOUNT)*DPH
        TLON_X(I_START,J)=WB_PARENT1+INC_LON*DLM
      ENDDO
!
      DO J=J_START,J_END
        KOUNT=0
        DO I=I_START+1,I_END
          KOUNT=KOUNT+1
          TLAT_X(I,J)=TLAT_X(I_START,J)
!         TLON_X(I,J)=TLON_X(I_START,J)+KOUNT*DLM
          TLON_X(I,J)=WB_PARENT1+(INC_LON+KOUNT)*DLM
        ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
!***  Convert from transformed to geographic coordinates.
!-----------------------------------------------------------------------
!
      DO J=J_START,J_END
        DX(J)=A_DLM*COS(TLAT_X(I_START,J))
      ENDDO

      DO J=J_START,J_END
      DO I=I_START,I_END
!
        ARG1=SIN(TLAT_X(I,J))*COS_TPH0                                  &
            +COS(TLAT_X(I,J))*SIN_TPH0*COS(TLON_X(I,J))
        GLAT_X(I,J)=ASIN(ARG1)
!
        ARG1=COS(TLAT_X(I,J))*COS(TLON_X(I,J))/(COS(GLAT_X(I,J))*COS_TPH0)
        ARG2=TAN(GLAT_X(I,J))*TAN_TPH0
        ADD=SIGN(1.,TLON_X(I,J))
        ARG3=ARG1-ARG2
        ARG3=SIGN(1.,ARG3)*MIN(ABS(ARG3),1.)                               !<-- Bound the argument of ACOS
        GLON_X(I,J)=TLM0+ADD*ACOS(ARG3)
!
        HDACX(I,J)=ACDT*DY*MAX(DX(J),DY)/(4.*DX(J)*DY)
        HDACY(I,J)=HDACX(I,J)
!
      ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
!***  Update the Coriolis parameter for V points.
!-----------------------------------------------------------------------
!
      IF(VELOCITY)THEN
!
        DO J=J_START,J_END
        DO I=I_START,I_END
!
          SIN_TPH=SIN(TLAT_X(I,J))
          COS_TPH=COS(TLAT_X(I,J))
          F(I,J)=TWOM*(COS_TPH0*SIN_TPH+SIN_TPH0*COS_TPH*COS(TLON_X(I,J)))
!
        ENDDO
        ENDDO
!
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE UPDATE_LATLON
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE GEO_TO_ROT(GLATX,GLONX,RLATX,RLONX)
!
!-----------------------------------------------------------------------
!***  Convert from geographic to rotated coordinates.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      REAL(kind=KFPT),INTENT(IN) :: GLATX                               &  !<-- Geographic latitude (radians)
                                   ,GLONX                                  !<-- Geographic longitude (radians)
!
      REAL(kind=KFPT),INTENT(OUT) :: RLATX                              &  !<-- Rotated latitude (radians)
                                    ,RLONX                                 !<-- Rotated longitude (radians)
!
!---------------------
!***  Local Variables
!---------------------
!
      REAL(kind=KFPT) :: X,Y,Z
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      X=COS(TPH0)*COS(GLATX)*COS(GLONX-TLM0)+SIN(TPH0)*SIN(GLATX)
      Y=COS(GLATX)*SIN(GLONX-TLM0)
      Z=-SIN(TPH0)*COS(GLATX)*COS(GLONX-TLM0)+COS(TPH0)*SIN(GLATX)
!
      RLATX=ATAN(Z/SQRT(X*X+Y*Y))
      RLONX=ATAN(Y/X)
      IF(X<0)RLONX=RLONX+PI
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE GEO_TO_ROT
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE FIX_SFC(MOVE_BUNDLE_H,NUM_FIELDS                       &
                        ,SEA_MASK                                       &
                        ,ILO,IHI,JLO,JHI                                &
                        ,I_START,I_END,J_START,J_END)
!
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      INTEGER(kind=KINT),INTENT(IN) :: ILO,IHI,JLO,JHI                  &  !<-- I,J subdomain limits
                                      ,I_START,I_END                    &  !<-- I limits of parent update region
                                      ,J_START,J_END                    &  !<-- J limits of parent update region
                                      ,NUM_FIELDS                          !<-- # of H-pt update variables after shift
!
      REAL(kind=KFPT),DIMENSION(ILO:IHI,JLO:JHI),INTENT(IN) :: SEA_MASK    !<-- This nest's sea mask (1=>water)
!
      TYPE(ESMF_FieldBundle),INTENT(INOUT) :: MOVE_BUNDLE_H                !<-- Bundle of internal state H arrays that shift
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: I,J,K,N_FIELD,NUM_DIMS,UPDATE_TYPE_INT
!
      INTEGER(kind=KINT) :: RC,RC_FIX
!
      INTEGER(kind=KINT),DIMENSION(1:3) :: LIMITS_HI,LIMITS_LO
!
      INTEGER(kind=KINT),DIMENSION(:,:),POINTER :: IARRAY_2D
!
      REAL(kind=KFPT) :: CHECK
!
      REAL(kind=KFPT),DIMENSION(:,:),POINTER :: ARRAY_2D
!
      REAL(kind=KFPT),DIMENSION(:,:,:),POINTER :: ARRAY_3D
!
      CHARACTER(len=1) :: UPDATE_TYPE_CHAR
!
      CHARACTER(len=25) :: FNAME 
!
      CHARACTER(len=99) :: FIELD_NAME
!
      LOGICAL(kind=KLOG) :: FOUND
!
      TYPE(ESMF_Field) :: HOLD_FIELD
!
      TYPE(ESMF_TypeKind_Flag) :: DATATYPE
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  We must clean the surface-related variables.  If the parent
!***  sends gridpoint data relevant for land but the nest reads its
!***  sea mask and the gridpoint is a water point then the nest
!***  will search for its own nearest water point and use that point's
!***  values for the variable at the conflict point.  Conversely if the
!***  parent sends gridpoint data relevant for water but the nest reads
!***  its sea mask and the gridpoint is land then the nest will search 
!***  for its own nearest land point that has land surface values and
!***  use those at the conflict point.
!
!***  This work could not be done earlier during the execution of the
!***  fields_h loop in subroutine UPDATE_INTERIOR_FROM_PARENT because
!***  the H-pt variables specified for updates after a domain moves
!***  are not listed in any particular order and thus all updates
!***  from the parent must be complete before this clean up can begin.
!
!***  NOTE: The calls to SEARCH_NEAR are being excluded temporarily
!           in order to ensure identical answers in moving nests
!           in the simplest manner for different task layouts.
!-----------------------------------------------------------------------
!
      all_fields: DO N_FIELD=1,NUM_FIELDS                                  !<-- Loop through H-pt variables again
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Again Extract Field from MOVE_BUNDLE_H"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_FieldBundleGet(FIELDBUNDLE=MOVE_BUNDLE_H              &  !<-- Bundle holding the H arrays for move updates
                                ,fieldIndex =N_FIELD                    &  !<-- Index of the Field in the Bundle
                                ,field      =HOLD_FIELD                 &  !<-- Field N_FIELD in the Bundle
                                ,rc         =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FIX)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Again Extract Field Information"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_FieldGet(field   =HOLD_FIELD                          &  !<-- Field N_FIELD in the Bundle
                          ,dimCount=NUM_DIMS                            &  !<-- Is this Field 2-D or 3-D?
                          ,typekind=DATATYPE                            &  !<-- Is the data integer or real?
                          ,name    =FIELD_NAME                          &  !<-- Name of the Field
                          ,rc      =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FIX)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Again Extract UPDATE_TYPE from Field"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(field=HOLD_FIELD                         &  !<-- Get Attribute from this Field
                              ,name ='UPDATE_TYPE'                      &  !<-- Name of the attribute to extract
                              ,value=UPDATE_TYPE_INT                    &  !<-- Value of the Attribute
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FIX)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        FNAME=TRIM(FIELD_NAME)
!
        IF(UPDATE_TYPE_INT==1)THEN  
          UPDATE_TYPE_CHAR='H'                                             !<-- Ordinary H-pt variable
        ELSEIF(UPDATE_TYPE_INT==2)THEN  
          UPDATE_TYPE_CHAR='L'                                             !<-- H-pt land surface variable
        ELSEIF(UPDATE_TYPE_INT==3)THEN
          UPDATE_TYPE_CHAR='W'                                             !<-- H-pt water surface variable
        ELSEIF(UPDATE_TYPE_INT==4)THEN
          UPDATE_TYPE_CHAR='F'                                             !<-- H-pt variable obtained from an external file
        ELSEIF(UPDATE_TYPE_INT==5)THEN
          UPDATE_TYPE_CHAR='V'                                             !<-- Ordinary V-pt variable
        ENDIF
!
!-----------------------------------------------------------------------
!***  We are only interested in water/land sfc variables, albedo,
!***  and the deep ground temperature.
!-----------------------------------------------------------------------
!
        IF(UPDATE_TYPE_CHAR/='W'                                        &
                  .AND.                                                 &
           UPDATE_TYPE_CHAR/='L'                                        &
                  .AND.                                                 &
           INDEX(FNAME,'ALB')==0                                        &
                  .AND.                                                 &
           INDEX(FNAME,'TYP')==0                                        &
                  .AND.                                                 &
           FNAME/='MXSNAL-move'                                         &
                  .AND.                                                 &
           FNAME/='TG-move' )THEN
!
          CYCLE all_fields
!
        ENDIF
!
!-----------------------------------------------------------------------
!***  Near coastlines the parent can generate valid values for both
!***  SST and for land variables.  The nest now sorts things out
!***  for each relevant variable based on its own sea mask.  We
!***  consider each variable separately since their default values
!***  differ and we do not want to fill a single DO loop with many
!***  IF tests.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Consider the 2-D Integer surface variables.
!-----------------------------------------------------------------------
!
        IF(NUM_DIMS==2                                                  &
             .AND.                                                      &
           DATATYPE==ESMF_TYPEKIND_I4)THEN
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Extract 2-D Integer Sfc Array"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_FieldGet(field    =HOLD_FIELD                       &  !<-- Field N_FIELD in the Bundle
                            ,localDe  =0                                &
                            ,farrayPtr=IARRAY_2D                        &  !<-- Dummy 2-D array with Field's Integer data
                            ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FIX)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
!---------------------
!***  Vegetation type
!---------------------
!
          IF(FNAME=='IVGTYP-move')THEN
            DO J=J_START,J_END
            DO I=I_START,I_END
              IF(SEA_MASK(I,J)>0.5)THEN 
                IARRAY_2D(I,J)=17                                          !<-- Set value at nest water point.
              ENDIF
            ENDDO
            ENDDO
!
          ENDIF
!
!---------------
!***  Soil type
!---------------
!
          IF(FNAME=='ISLTYP-move')THEN
            DO J=J_START,J_END
            DO I=I_START,I_END
              IF(SEA_MASK(I,J)>0.5)THEN 
                IARRAY_2D(I,J)=14                                          !<-- Set value at nest water point.
              ENDIF
            ENDDO
            ENDDO
!
          ENDIF
!
        ENDIF
!
!-----------------------------------------------------------------------
!***  Consider the 2-D Real surface variables.
!-----------------------------------------------------------------------
!
        IF(NUM_DIMS==2                                                  &
             .AND.                                                      &
           DATATYPE==ESMF_TYPEKIND_R4)THEN
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Extract 2-D Real Sfc Array"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_FieldGet(field    =HOLD_FIELD                       &  !<-- Field N_FIELD in the Bundle
                            ,localDe  =0                                &
                            ,farrayPtr=ARRAY_2D                         &  !<-- Dummy 2-D array with Field's Real data
                            ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FIX)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
!---------
!***  SST
!---------
!
          IF(FNAME=='SST-move')THEN
            DO J=J_START,J_END
            DO I=I_START,I_END
!
              FOUND=.TRUE.
              IF(SEA_MASK(I,J)<0.5)THEN 
                ARRAY_2D(I,J)=0.                                           !<-- Set dummy value at nest land point.
!
              ELSEIF(SEA_MASK(I,J)>0.5.AND.ARRAY_2D(I,J)<1.)THEN           !<-- Parent sent land value to nest water point.
!
                FOUND=.FALSE.
!
                IF(domain_int_state%SFC_CONFLICT=='nearest')THEN
                  CALL SEARCH_NEAR(FNAME,SEA_MASK,I,J                   & 
                                  ,ILO,IHI,JLO,JHI                      &
                                  ,I_START,I_END,J_START,J_END          &
                                  ,LIMITS_LO(3),LIMITS_HI(3)            &
                                  ,FOUND                                &
                                  ,array_2d=ARRAY_2D )
                ENDIF
!
              ENDIF
!
              IF(.NOT.FOUND)THEN                           
                ARRAY_2D(I,J)=300.                                         !<-- Made-up sea sfc temperature
              ENDIF                               
!
            ENDDO
            ENDDO
!
          ENDIF
!
!-----------------
!***  Base Albedo
!-----------------
!
          IF(FNAME=='ALBASE-move')THEN
            DO J=J_START,J_END
            DO I=I_START,I_END
!
              FOUND=.TRUE.
              IF(SEA_MASK(I,J)>0.5)THEN 
                ARRAY_2D(I,J)=0.06                                         !<-- Set value at nest water point.
              ENDIF
!
            ENDDO
            ENDDO
!
          ENDIF
!
!--------------------
!***  Dynamic Albedo
!--------------------
!
          IF(FNAME=='ALBEDO-move')THEN
            DO J=J_START,J_END
            DO I=I_START,I_END
!
              FOUND=.TRUE.
              IF(SEA_MASK(I,J)>0.5)THEN 
                ARRAY_2D(I,J)=0.06                                         !<-- Set water value at nest water point.
!
              ELSEIF(SEA_MASK(I,J)<0.5)THEN
                CHECK=ABS(ARRAY_2D(I,J)-0.06)
                IF(CHECK<1.E-5)THEN                                        !<-- Parent sent water value to nest land point.
!
                  FOUND=.FALSE.
!
                  IF(domain_int_state%SFC_CONFLICT=='nearest')THEN
                    CALL SEARCH_NEAR(FNAME,SEA_MASK,I,J                 & 
                                    ,ILO,IHI,JLO,JHI                    &
                                    ,I_START,I_END,J_START,J_END        &
                                    ,LIMITS_LO(3),LIMITS_HI(3)          &
                                    ,FOUND                              &
                                    ,array_2d=ARRAY_2D )
                  ENDIF
!
                ENDIF
!
              ENDIF
!
              IF(.NOT.FOUND)THEN                           
                ARRAY_2D(I,J)=0.25                                         !<-- Made-up albedo over land
              ENDIF                               
!
            ENDDO
            ENDDO
!
          ENDIF
!
!-----------------------------
!***  Deep ground temperature
!-----------------------------
!
          IF(FNAME=='TG-move')THEN
            DO J=J_START,J_END
            DO I=I_START,I_END
!
              IF(SEA_MASK(I,J)>0.5)THEN 
                ARRAY_2D(I,J)=273.16                                       !<-- Set water value at nest water point.
              ENDIF
!
            ENDDO
            ENDDO
!
          ENDIF
!
!------------
!***  MXSNAL
!------------
!
          IF(FNAME=='MXSNAL-move')THEN
            DO J=J_START,J_END
            DO I=I_START,I_END
!
              IF(SEA_MASK(I,J)>0.5)THEN 
                ARRAY_2D(I,J)=0.08                                         !<-- Set water value at nest water point.
              ENDIF
!
            ENDDO
            ENDDO
!
          ENDIF
!
        ENDIF
!
!-----------------------------------------------------------------------
!***  Consider the 3-D Real surface variables.
!-----------------------------------------------------------------------
!
        IF(NUM_DIMS==3                                                  &
             .AND.                                                      &
           DATATYPE==ESMF_TYPEKIND_R4)THEN
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Extract 3-D Real Sfc Array"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_FieldGet(field      =HOLD_FIELD                     &  !<-- Field N_FIELD in the Bundle
                            ,localDe    =0                              &
                            ,farrayPtr  =ARRAY_3D                       &  !<-- Dummy 2-D array with Field's Real data
                            ,totalLBound=LIMITS_LO                      &  !<-- Starting index in each dimension
                            ,totalUBound=LIMITS_HI                      &  !<-- Ending index in each dimension
                            ,rc         =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FIX)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!---------
!***  STC
!---------
!
          IF(FNAME=='STC-move')THEN
            DO J=J_START,J_END
            DO I=I_START,I_END
!
              FOUND=.TRUE.
              IF(SEA_MASK(I,J)>0.5)THEN                
                DO K=LIMITS_LO(3),LIMITS_HI(3)
                  ARRAY_3D(I,J,K)=273.16                                   !<-- Set dummy value at nest water point.
                ENDDO
!
              ELSEIF(SEA_MASK(I,J)<0.5)THEN                     
                CHECK=ABS(ARRAY_3D(I,J,1)-273.16)
                IF(CHECK<1.E-2)THEN                                        !<-- Parent sent water value to nest land point.
!
                  FOUND=.FALSE.
!
                  IF(domain_int_state%SFC_CONFLICT=='nearest')THEN
                    CALL SEARCH_NEAR(FNAME,SEA_MASK,I,J                 & 
                                    ,ILO,IHI,JLO,JHI                    &
                                    ,I_START,I_END,J_START,J_END        &
                                    ,LIMITS_LO(3),LIMITS_HI(3)          &
                                    ,FOUND                              &
                                    ,array_3d=ARRAY_3D )
                  ENDIF
!
                ENDIF
!
              ENDIF
!
              IF(.NOT.FOUND)THEN                           
                DO K=LIMITS_LO(3),LIMITS_HI(3)
                  ARRAY_3D(I,J,K)=285.+K*2.                                !<-- Made-up soil temperature
                ENDDO
              ENDIF                               
!
            ENDDO
            ENDDO
!
          ENDIF
!
!--------------
!***  SMC/SH2O
!--------------
!
          IF(FNAME=='SMC-move'.OR.FNAME=='SH2O-move')THEN
            DO J=J_START,J_END
            DO I=I_START,I_END
!
              FOUND=.TRUE.
              IF(SEA_MASK(I,J)>0.5)THEN                
                DO K=LIMITS_LO(3),LIMITS_HI(3)
                  ARRAY_3D(I,J,K)=1.0                                      !<-- Set dummy value at nest water point.
                ENDDO
!
              ELSEIF(SEA_MASK(I,J)<0.5.AND.ARRAY_3D(I,J,1)>0.9)THEN        !<-- Parent sent water value to nest land point.
!
!*** Temporary exclusion of SEARCH_NEAR
                FOUND=.FALSE.
!
                IF(domain_int_state%SFC_CONFLICT=='nearest')THEN
                  CALL SEARCH_NEAR(FNAME,SEA_MASK,I,J                   & 
                                  ,ILO,IHI,JLO,JHI                      &
                                  ,I_START,I_END,J_START,J_END          &
                                  ,LIMITS_LO(3),LIMITS_HI(3)            &
                                  ,FOUND                                &
                                  ,array_3d=ARRAY_3D )
                ENDIF
!
              ENDIF
!
              IF(.NOT.FOUND)THEN                           
                DO K=LIMITS_LO(3),LIMITS_HI(3)
                  ARRAY_3D(I,J,K)=0.2                                      !<-- Made-up soil moisture
                ENDDO
              ENDIF                               
!
            ENDDO
            ENDDO
!
          ENDIF
!
        ENDIF
!
!-----------------------------------------------------------------------
!
        CYCLE all_fields
!
!-----------------------------------------------------------------------
!
      ENDDO all_fields
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE FIX_SFC
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE SEARCH_NEAR(FNAME,SEA_MASK,I_IN,J_IN                   &
                            ,ILO,IHI,JLO,JHI                            &
                            ,I_START,I_END,J_START,J_END                &
                            ,KLO,KHI                                    &
                            ,FOUND                                      &
                            ,ARRAY_2D                                   &
                            ,ARRAY_3D )
!
!-----------------------------------------------------------------------
!***  Search for nearest points to given conflict points on a nest
!***  after it has moved.  The search begins at the point nearest to
!***  the one in question then moves increasingly farther away.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      INTEGER(kind=KINT),INTENT(IN) :: I_IN,J_IN                        &  !<-- The conflict point 
                                      ,I_START,I_END,J_START,J_END      &  !<-- Limits of nest update region by parent
                                      ,ILO,IHI,JLO,JHI                  &  !<-- Nest subdomain dimensions
                                      ,KHI,KLO                             !<-- Vertical dimension limits of 3-D soil arrays
!
      REAL(kind=KFPT),DIMENSION(ILO:IHI,JLO:JHI),INTENT(IN) :: SEA_MASK    !<-- Nest's sea mask (1=>water)
!
      CHARACTER(len=*),INTENT(IN) :: FNAME                                 !<-- Name of the variable being considered
!
      LOGICAL(kind=KLOG),INTENT(OUT) :: FOUND                              !<-- Was a valid point found by the search?
!
      REAL(kind=KFPT),DIMENSION(ILO:IHI,JLO:JHI),INTENT(INOUT)          &
                                                  ,OPTIONAL :: ARRAY_2D    !<-- 2-D land/water array to repair
!
      REAL(kind=KFPT),DIMENSION(ILO:IHI,JLO:JHI,KLO:KHI),INTENT(INOUT)  &  !<-- 3-D soil array to repair
                                                  ,OPTIONAL :: ARRAY_3D
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: I1,I2,J1,J2
!
      INTEGER(kind=KINT) :: I_SEARCH,J_SEARCH,K,N_SEARCH
!
      REAL(kind=KFPT) :: CHECK
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      FOUND=.FALSE.
!
      I1=ILO+solver_int_state%NHALO                                        !<--  Local integration
      I2=IHI-solver_int_state%NHALO                                        !     limits of this
      J1=JLO+solver_int_state%NHALO                                        !     task's subdomain.
      J2=JHI-solver_int_state%NHALO                                        !<--
!
!-----------------------------------------------------------------------
!***  If the given nest point following the move is a water point
!***  based on the nest's reading its own sea mask but the value
!***  from the parent is a land value then the nest searches for
!***  its nearest legitimate water value and gives that to the point
!***  in question.
!-----------------------------------------------------------------------
!
!---------
!***  SST
!---------
!
      IF(FNAME=='SST-move')THEN       
!
        DO N_SEARCH=2,N_PTS_SEARCH
          I_SEARCH=I_IN+I_SEARCH_INC(N_SEARCH)
          J_SEARCH=J_IN+J_SEARCH_INC(N_SEARCH)
!
          IF(I_SEARCH<I1.OR.I_SEARCH>I2                                 &  !<-- Keep the search on the task subdomain
                         .OR.                                           &  !
             J_SEARCH<J1.OR.J_SEARCH>J2)CYCLE                              !<--
!
          IF(ARRAY_2D(I_SEARCH,J_SEARCH)>1.)THEN                           !<-- If true, the nest found its own nearest water point
            ARRAY_2D(I_IN,J_IN)=ARRAY_2D(I_SEARCH,J_SEARCH)
            FOUND=.TRUE.
            EXIT
          ENDIF
!
        ENDDO
!
      ENDIF
!
!-----------------------------------------------------------------------
!***  If the given nest point following the move is a land point
!***  based on the nest's reading its own sea mask but the value
!***  from the parent is a water value then the nest searches for
!***  its nearest legitimate land value and gives that to the point
!***  in question.  The varying dummy values for the different land
!***  variables forces individual searches.
!-----------------------------------------------------------------------
!
!---------
!***  STC
!---------
!
      IF(FNAME=='STC-move')THEN       
!
        DO N_SEARCH=2,N_PTS_SEARCH
          I_SEARCH=I_IN+I_SEARCH_INC(N_SEARCH)
          J_SEARCH=J_IN+J_SEARCH_INC(N_SEARCH)
!
          IF(I_SEARCH<I1.OR.I_SEARCH>I2                                 &  !<-- Keep the search on the task subdomain
                         .OR.                                           &  !
             J_SEARCH<J1.OR.J_SEARCH>J2)CYCLE                              !<--
!
          CHECK=ABS(ARRAY_3D(I_SEARCH,J_SEARCH,1)-273.16) 
!
          IF(CHECK>1.E-2)THEN                                              !<-- Make sure the search point has a valid land value
            DO K=KLO,KHI
              ARRAY_3D(I_IN,J_IN,K)=ARRAY_3D(I_SEARCH,J_SEARCH,K)
            ENDDO
            FOUND=.TRUE.
            EXIT
          ENDIF
!
        ENDDO
!
      ENDIF
!
!--------------
!***  SMC/SH2O
!--------------
!
      IF(FNAME=='STC-move'.OR.FNAME=='SH2O-move')THEN       
!
        DO N_SEARCH=2,N_PTS_SEARCH
          I_SEARCH=I_IN+I_SEARCH_INC(N_SEARCH)
          J_SEARCH=J_IN+J_SEARCH_INC(N_SEARCH)
!
          IF(I_SEARCH<I1.OR.I_SEARCH>I2                                 &  !<-- Keep the search on the task subdomain
                         .OR.                                           &  !
             J_SEARCH<J1.OR.J_SEARCH>J2)CYCLE                              !<--
!
          IF(ARRAY_3D(I_SEARCH,J_SEARCH,1)<0.9)THEN                        !<-- Make sure the search point has a valid land value
            DO K=KLO,KHI
              ARRAY_3D(I_IN,J_IN,K)=ARRAY_3D(I_SEARCH,J_SEARCH,K)
            ENDDO
            FOUND=.TRUE.
            EXIT
          ENDIF
!
        ENDDO
!
      ENDIF
!
!------------
!***  Albedo
!------------
!
      IF(FNAME=='ALBEDO-move')THEN       
!
        DO N_SEARCH=2,N_PTS_SEARCH
          I_SEARCH=I_IN+I_SEARCH_INC(N_SEARCH)
          J_SEARCH=J_IN+J_SEARCH_INC(N_SEARCH)
!
          IF(I_SEARCH<I1.OR.I_SEARCH>I2                                 &  !<-- Keep the search on the task subdomain
                         .OR.                                           &  !
             J_SEARCH<J1.OR.J_SEARCH>J2)CYCLE                              !<--
!
          CHECK=ABS(ARRAY_2D(I_SEARCH,J_SEARCH)-0.06)
!
          IF(CHECK>1.E-2)THEN                                              !<-- Make sure the search point has a valid land value
            ARRAY_2D(I_IN,J_IN)=ARRAY_2D(I_SEARCH,J_SEARCH)
            FOUND=.TRUE.
            EXIT
          ENDIF
!
        ENDDO
!
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE SEARCH_NEAR
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE SEARCH_INIT
!
!-----------------------------------------------------------------------
!***  Generate I and J increments from any given point to all other
!***  points surrounding it in a square of given width based on the
!***  distances from the central point ranging from smallest to
!***  largest distance.
!-----------------------------------------------------------------------
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: KOUNT=0
!
      INTEGER(kind=KINT) :: I,I_CENTER,ISTAT,J,J_CENTER                 &
                           ,N_WIDTH,RC
!
      TYPE(DIST),POINTER :: LARGE=>NULL()                               &
                           ,SMALL=>NULL()                               &
                           ,PTR  =>NULL()                               &
                           ,PTR1 =>NULL()                               &
                           ,PTR2 =>NULL()                               &
                           ,PTRX =>NULL()
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      N_WIDTH=2*N_PTS_SEARCH_WIDTH+1                                       !<-- Width of the search region
      I_CENTER=N_PTS_SEARCH_WIDTH+1                                        !<-- Relative I at the center of the search region
      J_CENTER=N_PTS_SEARCH_WIDTH+1                                        !<-- Relative J at the center of the search region
!
      SMALL=>SMALLX(MY_DOMAIN_ID)
      LARGE=>LARGEX(MY_DOMAIN_ID)
!
!-----------------------------------------------------------------------
!
      DO J=1,N_WIDTH
      DO I=1,N_WIDTH
!
!-----------------------------------------------------------------------
!***  We must allocate a pointer to each gridpoint in the search.
!-----------------------------------------------------------------------
!
        ALLOCATE(PTR,stat=ISTAT)
        IF(ISTAT/=0)THEN
          WRITE(0,*)' Failed to allocate search pointer  stat=',ISTAT
          CALL ESMF_Finalize(endflag=ESMF_END_ABORT                 &
                            ,rc             =RC)
        ENDIF
!
!-----------------------------------------------------------------------
!***  Compute the distance from the center point to all of the other
!***  points in the square.  We are not accounting for a map projection
!***  at this time.
!-----------------------------------------------------------------------
!
        PTR%I_INC=I-I_CENTER
        PTR%J_INC=J-J_CENTER
        PTR%VALUE=SQRT(REAL(PTR%I_INC*PTR%I_INC                         &
                           +PTR%J_INC*PTR%J_INC))
!
!-----------------------------------------------------------------------
!***  Sort the distances to each point in the square as they are
!***  computed going from smallest to largest.
!-----------------------------------------------------------------------
!
        new_val: IF(I==1.AND.J==1)THEN                                     !<-- 1st value is both smallest/largest
          SMALL=>PTR
          LARGE=>SMALL
          NULLIFY(PTR%NEXT_VALUE)
!
        ELSE                                                               !<-- All subsequent values
!
          IF(PTR%VALUE<SMALL%VALUE)THEN                                    !<-- New value smaller than previous smallest
            PTR%NEXT_VALUE=>SMALL
            SMALL=>PTR
!
          ELSEIF(PTR%VALUE>=LARGE%VALUE)THEN                               !<-- New value same or larger than previous largest
            LARGE%NEXT_VALUE=>PTR
            LARGE=>PTR
            NULLIFY(LARGE%NEXT_VALUE)
!
          ELSE                                                             !<-- New value between current smallest and largest
            PTR1=>SMALL
            PTR2=>PTR1%NEXT_VALUE
!
            search:DO                                                      !<-- Find new value's proper place in the list
!
              IF(PTR%VALUE>=PTR1%VALUE.AND.PTR%VALUE<PTR2%VALUE)THEN
                PTR%NEXT_VALUE=>PTR2
                PTR1%NEXT_VALUE=>PTR
                EXIT search
              ENDIF
!
              PTR1=>PTR2
              PTR2=>PTR2%NEXT_VALUE
!
            ENDDO search
!
          ENDIF
!
        ENDIF new_val
!
      ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
!***  We now have all the distances from the central point in order
!***  from smallest to largest.  Fill the 1-D index increment arrays
!***  for I and J that will allow the search to go from the smallest
!***  distance to the largest when those arrays are stepped through.
!***  As the index increments to each gridpoint are stored we can
!***  deallocate the pointer to that gridpoint since it is no
!***  longer needed.
!-----------------------------------------------------------------------
!
      PTR=>SMALL                                                           !<-- Begin with the nearest gridpoint
      KOUNT=0
!
      DO 
        KOUNT=KOUNT+1
        IF(.NOT.ASSOCIATED(PTR))EXIT
!       WRITE(0,23331)KOUNT,PTR%VALUE
23331   FORMAT(' Value #',I6,' is ',F10.6)    
        I_SEARCH_INC(KOUNT)=PTR%I_INC                                      !<-- Store the increments of I and J to the next
        J_SEARCH_INC(KOUNT)=PTR%J_INC                                      !    gridpoint in the distance list.
        PTRX=>PTR
        PTR=>PTR%NEXT_VALUE                                                !<-- Proceed to the next gridpoint further away
        DEALLOCATE(PTRX)                                                   !<-- Deallocate the previous gridpoint's pointer
      ENDDO
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE SEARCH_INIT
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE RESET_SFC_VARS(SFC_FILE_RATIO                          &
                               ,GLAT_H                                  &
                               ,GLON_H                                  &
                               ,MOVE_BUNDLE)
!
!-----------------------------------------------------------------------
!***  If the parent initialized this moving nest from the parent's
!***  own initial state then the nest will now reinitialize those
!***  2-D sfc fields that are constant in time and that the nest
!***  reads from external files for data replacement in parent 
!***  update regions during the integration.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument variables   
!------------------------
!
      INTEGER(kind=KINT),INTENT(IN) :: SFC_FILE_RATIO                      !<-- Ratio of upper parent grid increment to this domain's
!
      REAL(kind=KFPT),DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: GLAT_H   &  !<-- Geographic latitude (radians) at H pts on nest domain.
                                                              ,GLON_H      !<-- Geographic longitude (radians) at H pts on nest domain.
!
      TYPE(ESMF_FieldBundle),INTENT(INOUT) :: MOVE_BUNDLE
!
!---------------------
!***  Local variables
!---------------------
!
      INTEGER(kind=KINT),SAVE :: I_OFFSET,J_OFFSET
!
      INTEGER(kind=KINT) :: I,I_CORNER,IEND,ILOC                         &
                           ,INPUT_NEST,ISTART,ITE_X                      &
                           ,J,J_CORNER,JEND,JSTART,JTE_Y                 &
                           ,N_FIELD,N_REMOVE,NN,NUM_FIELDS               &
                           ,UPDATE_TYPE_INT
!
      INTEGER(kind=KINT) :: I_COUNT_DATA,I_START_DATA                   &
                           ,J_COUNT_DATA,J_START_DATA                   &
                           ,NCID,NCTYPE,NDIMS,VAR_ID
!
      INTEGER(kind=KINT) :: IERR,RC,RC_RES
!
      INTEGER(kind=KINT),DIMENSION(:,:),POINTER :: IARRAY_2D=>NULL()
!
      INTEGER(kind=KINT),DIMENSION(1:2) :: DIM_IDS
!
      REAL(kind=KFPT) :: GBL,REAL_I,REAL_J
!
      REAL(kind=KFPT),DIMENSION(:,:),POINTER :: ARRAY_2D=>NULL()
!
      CHARACTER(len=2)  :: ID_SFC_FILE
      CHARACTER(len=15) :: VNAME
      CHARACTER(len=99) :: FIELD_NAME,FILENAME
!
      LOGICAL(kind=KLOG),SAVE :: FIRST=.TRUE.
      LOGICAL(kind=KLOG) :: OPENED
!
      TYPE(ESMF_Field) :: HOLD_FIELD
      TYPE(ESMF_TypeKind_Flag) :: DATATYPE
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      IF(GLOBAL_TOP_PARENT)THEN
        GBL=1.
      ELSE
        GBL=0.
      ENDIF
!
!-----------------------------------------------------------------------
!***  The nest uses its GLAT and GLON to determine exactly where it
!***  lies on the uppermost parent grid and thus where its grid lies
!***  within the nest-resolution sfc data in the external file.
!***  Find the I,J on the uppermost parent grid on which the SW corner
!***  of this nest task lies.
!-----------------------------------------------------------------------
!
      I_CORNER=MAX(IMS,IDS)+GBL                                          !<-- Nest task halos are covered with data
      J_CORNER=MAX(JMS,JDS)+GBL                                          !
!
      CALL LATLON_TO_IJ(GLAT_H(I_CORNER,J_CORNER)                     &  !<-- Geographic latitude of nest task subdomain SW corner
                       ,GLON_H(I_CORNER,J_CORNER)                     &  !<-- Geographic longitude of nest task subdomain SW corner
                       ,TPH0_1,TLM0_1                                 &  !<-- Central lat/lon (radians, N/E) of uppermost parent
                       ,SB_1,WB_1                                     &  !<-- Rotated lat/lon of upper parent's S/W bndry (radians, N/E)
                       ,RECIP_DPH_1,RECIP_DLM_1                       &  !<-- Reciprocal of I/J grid increments (radians) on upper parent
                       ,GLOBAL_TOP_PARENT                             &  !<-- Is the uppermost daomin on a global grid?
                       ,REAL_I                                        &  !<-- Corresponding I index on uppermost parent grid
                       ,REAL_J)
!
      I_OFFSET=NINT((REAL_I-1.-GBL)*SFC_FILE_RATIO)                      !<-- Offset in I between sfc file index and nest index
      J_OFFSET=NINT((REAL_J-1.-GBL)*SFC_FILE_RATIO)                      !<-- Offset in J between sfc file index and nest index
!
      ITE_X =MIN(IME,IDE)                                                !<-- Last task point to update in I
      I_START_DATA=I_OFFSET+1                                            !<-- Start reading at this I in external data array
      I_COUNT_DATA=ITE_X-I_CORNER+1                                      !<-- # of points to read in I
!
      JTE_Y =MIN(JME,JDE)                                                !<-- Last task point to update in J
      J_START_DATA=J_OFFSET+1                                            !<-- Start reading at this J in external data array
      J_COUNT_DATA=JTE_Y-J_CORNER+1                                      !<-- # of points to read in J
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="How many Fields in the Move_Bundle?"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_FieldBundleGet(FIELDBUNDLE=MOVE_BUNDLE                  &  !<-- Bundle holding the arrays for move updates
                              ,fieldCount =NUM_FIELDS                   &  !<-- # of Fields in this Bundle
                              ,rc         =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RES)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Loop through this Bundles's Fields and replace the values of
!***  those that are associated with the external nest-resolution
!***  surface data files.
!-----------------------------------------------------------------------
!
      field_loop: DO N_FIELD=1,NUM_FIELDS
!
!-----------------------------------------------------------------------
!
        CALL ESMF_FieldBundleGet(FIELDBUNDLE=MOVE_BUNDLE                &  !<-- Bundle holding the arrays for move updates
                                ,fieldIndex =N_FIELD                    &  !<-- Index of the Field in the Bundle
                                ,field      =HOLD_FIELD                 &  !<-- Field N_FIELD in the Bundle
                                ,rc         =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RES)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_FieldGet(field   =HOLD_FIELD                          &  !<-- Field N_FIELD in the Bundle
                          ,typeKind=DATATYPE                            &  !<-- Does this Field contain an integer or real array?
                          ,name    =FIELD_NAME                          &  !<-- The name of the Field
                          ,rc      =RC )
!
        N_REMOVE=INDEX(FIELD_NAME,SUFFIX_MOVE)
        FIELD_NAME=FIELD_NAME(1:N_REMOVE-1)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="RESET_SFC_VARS: Extract UPDATE_TYPE from Field"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(field=HOLD_FIELD                         &  !<-- Get Attribute from this Field
                              ,name ='UPDATE_TYPE'                      &  !<-- Name of the attribute to extract
                              ,value=UPDATE_TYPE_INT                    &  !<-- Value of the Attribute
                              ,rc   =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RES)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  If this Field has data in an external sfc file then open the file.
!-----------------------------------------------------------------------
!
        filedata: IF(UPDATE_TYPE_INT==4)THEN                               !<-- This means the variable has an external file
!
!-----------------------------------------------------------------------
!
          IF(SFC_FILE_RATIO<=9)THEN
            WRITE(ID_SFC_FILE,'(I1.1)')SFC_FILE_RATIO
          ELSEIF(SFC_FILE_RATIO>=10)THEN
            WRITE(ID_SFC_FILE,'(I2.2)')SFC_FILE_RATIO
          ENDIF
!
          FILENAME=TRIM(FIELD_NAME)//'_'//TRIM(ID_SFC_FILE)//'.nc'
!
          CALL CHECK(NF90_OPEN(FILENAME,NF90_NOWRITE,NCID))                !<-- Open the current field's external netCDF file.
!
!-----------------------------------------------------------------------
!***  Extract the array from the Field.
!-----------------------------------------------------------------------
!
!----------
!***  Real
!----------
!
          IF(DATATYPE==ESMF_TYPEKIND_R4)THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="RESET_SFC_VARS: Extract 2-D Real Array for Type F"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_FieldGet(field    =HOLD_FIELD                     &  !<-- Field N_FIELD in the Bundle
                              ,localDe  =0                              &
                              ,farrayPtr=ARRAY_2D                       &  !<-- Dummy 2-D array with Field's Real data
                              ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RES)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL CHECK(NF90_INQUIRE_VARIABLE(NCID,3,VNAME,NCTYPE  &
                                            ,NDIMS,DIM_IDS))
            CALL CHECK(NF90_INQ_VARID(NCID,VNAME,VAR_ID))
!
            CALL CHECK(NF90_GET_VAR(NCID,VAR_ID                             &  !<-- Extract the desired real values from the
                                   ,ARRAY_2D(I_CORNER:ITE_X,J_CORNER:JTE_Y) &  !    current field's external file.
                                   ,start=(/I_START_DATA,J_START_DATA/)     &  !    Nest points that have moved beyond the
                                   ,count=(/I_COUNT_DATA,J_COUNT_DATA/)))      !    pre-move footprint are updated.
!
!-------------
!***  Integer
!-------------
!
          ELSEIF(DATATYPE==ESMF_TYPEKIND_I4)THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="RESET_SFC_VARS: Extract 2-D Integer Array for Type F"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_FieldGet(field    =HOLD_FIELD                     &  !<-- Field N_FIELD in the Bundle
                              ,localDe  =0                              &
                              ,farrayPtr=IARRAY_2D                      &  !<-- Dummy 2-D array with Field's Real data
                              ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RES)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL CHECK(NF90_INQUIRE_VARIABLE(NCID,3,VNAME,NCTYPE  &
                                            ,NDIMS,DIM_IDS))
            CALL CHECK(NF90_INQ_VARID(NCID,VNAME,VAR_ID))
!
            CALL CHECK(NF90_GET_VAR(NCID,VAR_ID                              &  !<-- Extract the desired integer values from the
                                   ,IARRAY_2D(I_CORNER:ITE_X,J_CORNER:JTE_Y) &  !    current field's external file.
                                   ,start=(/I_START_DATA,J_START_DATA/)      &  !    Nest points that have moved beyond the
                                   ,count=(/I_COUNT_DATA,J_COUNT_DATA/)))       !    pre-move footprint are updated.
!
          ENDIF
!
!-----------------------------------------------------------------------
!***  Read the appropriate section of the external data to reset the
!***  values of this sfc variable from those that were originally
!***  interpolated from the parent.
!-----------------------------------------------------------------------
!
          CALL CHECK(NF90_CLOSE(NCID))                                     !<-- Close the external netCDF file.
!
!-----------------------------------------------------------------------
!
        ENDIF  filedata
!
        CALL CHECK(NF90_CLOSE(NCID))                                       !<-- Close the external netCDF file.
!
!-----------------------------------------------------------------------
!
      ENDDO  field_loop
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE RESET_SFC_VARS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      END MODULE module_DOMAIN_GRID_COMP
!
!-----------------------------------------------------------------------
