!-----------------------------------------------------------------------
!
      MODULE MODULE_WRITE_GRID_COMP
!
!-----------------------------------------------------------------------
!***  The Write gridded component.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Data was put into this component's import state destined for
!***  history output.  This component extracts that information
!***  from the import state whose contents are seen only by the
!***  forecast tasks and transfers 2-D data to groups of write tasks
!***  where it is partially reassembled.  The write tasks then
!***  transfer their subsections to the lead write task which
!***  assembles the 2-D data onto the full domain and writes out
!***  all scalar/1-D/2-D data to a history file.
!-----------------------------------------------------------------------
!***
!***  HISTORY   
!***
!       xx Feb 2007:  W. Yang - Originator
!       13 Jun 2007:  T. Black - Name revisions in CPL_REGISTER  
!                                and CPL_INITIALIZE
!       14 Aug 2007:  T. Black - Revised CPL_RUN for general output
!                                selection and added documentation
!                                for users.
!       12 Sep 2007:  T. Black - Replaced the write component and the
!                                write gridded component with only
!                                a gridded component that contains
!                                quilting.
!          Mar 2008:  R. Vasic - Convert from ESMF 3.0.1 to 3.1.0
!       15 Aug 2008:  J. Wang  - Revised for addition of NEMS-IO
!       16 Sep 2008:  J. Wang  - Output array reverts from 3-D to 2-D
!       14 Oct 2008:  R. Vasic - Add restart capability
!       05 Jan 2009:  J. Wang  - Add 10-m wind factor into NMMB
!                                runhistory and restart files
!       02 Jul 2009:  J. Wang  - Added fcstdone/restartdone files
!       04 Sep 2009:  T. Black - Merged trunk and NMM-B nesting 
!                                versions.
!       07 May 2010:  T. Black - Change output frequency to minutes.
!       16 Dec 2010:  J. Wang  - Change to nemsio library
!          Feb 2011:  W. Yang  - Updated to use both the ESMF 4.0.0rp2 library,
!                                ESMF 5 library and the the ESMF 3.1.0rp2 library.
!       12 May 2011   W. Yang  - Modified for using the ESMF 5.2.0r_beta_snapshot_07.
!       23 May 2011   J. Wang  - add do post option
!       27 SEP 2011   W. Yang  - Modified for using the ESMF 5.2.0r library.
!       26 JAN 2016   J. Carley - Add error checking after calls to nemsio_writerec
!---------------------------------------------------------------------------------
!
      USE MPI
      USE ESMF
      USE MODULE_WRITE_INTERNAL_STATE
      USE MODULE_WRITE_ROUTINES,ONLY : OPEN_HST_FILE                    &
                                      ,OPEN_RST_FILE                    &
                                      ,WRITE_RUNHISTORY_OPEN            &
                                      ,SEND_UPDATED_ATTRIBUTES          &
                                      ,WRITE_NEMSIO_RUNHISTORY_OPEN     &
                                      ,WRITE_RUNRESTART_OPEN            &
                                      ,WRITE_NEMSIO_RUNRESTART_OPEN     &
                                      ,TIME_FOR_HISTORY                 &
                                      ,TIME_FOR_RESTART
!
      USE MODULE_DM_PARALLEL,ONLY : PARA_RANGE                          &
                                   ,MAX_GROUPS                          &
                                   ,MPI_COMM_COMP                       &
                                   ,MPI_INTERCOMM_ARRAY
!
      USE MODULE_CONTROL,ONLY : TIMEF
      USE MODULE_GET_CONFIG_WRITE
      USE MODULE_ERROR_MSG,ONLY : ERR_MSG,MESSAGE_CHECK
      USE MODULE_KINDS
      USE MODULE_CONSTANTS,ONLY : G
      USE NEMSIO_MODULE
      USE MODULE_BGRID_INTERP,ONLY: V_TO_H_BGRID
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------

      PRIVATE
!
      PUBLIC :: WRITE_REGISTER                                          &
               ,WRITE_SETUP                                             &
               ,WRITE_DESTROY
!
      PUBLIC :: write_init_tim                                          &
               ,write_run_tim                                           &
               ,write_first_tim                                         &
               ,write_recv_outp_tim                                     &
               ,write_send_outp_tim                                     &
               ,write_get_fields_tim
!
!-----------------------------------------------------------------------
!
      INTEGER(kind=KINT),PARAMETER :: MAX_LENGTH_I1D=5000               &  !<-- Max words in all 1-D integer history variables
                                     ,MAX_LENGTH_R1D=25000              &  !<-- Max words in all 1-D real history variables
                                     ,MAX_LENGTH_LOG=MAX_DATA_LOG          !<-- Max logical variables
!
      INTEGER(kind=KINT),SAVE :: ITS,ITE,JTS,JTE                        &  !<-- Integration grid limits on each task subdomain
                                ,IDS,IDE,JDS,JDE                        &  !<-- Full domain horizontal index limits
                                ,LAST_FCST_TASK                         &  !<-- Rank of the last Forecast task
                                ,LEAD_WRITE_TASK                        &  !<-- Rank of the lead (first) Write task in this Write group
                                ,LAST_WRITE_TASK                        &  !<-- Rank of the last Write task the Write group
                                ,LM                                     &  !<-- # of model layers
                                ,LNSH                                   &  !<-- H Rows in boundary region
                                ,LNSV                                   &  !<-- V Rows in boundary region
                                ,NLEV_H                                 &  !<-- Total # of levels in H-pt boundary variables
                                ,NLEV_V                                 &  !<-- Total # of levels in V-pt boundary variables
                                ,NTASKS                                 &  !<-- # of Write tasks in the current group + all Fcst tasks
                                ,NVARS_BC_2D_H                          &  !<-- # of 2-D H-pt boundary variables
                                ,NVARS_BC_3D_H                          &  !<-- # of 3-D H-pt boundary variables
                                ,NVARS_BC_4D_H                          &  !<-- # of 4-D H-pt boundary variables
                                ,NVARS_BC_2D_V                          &  !<-- # of 2-D V-pt boundary variables
                                ,NVARS_BC_3D_V                          &  !<-- # of 3-D V-pt boundary variables
                                ,NUM_DOMAINS_TOTAL                      &  !<-- Total # of domains
                                ,NUM_PES_FCST                           &  !<-- # of Fcst tasks in the current group + all Fcst tasks
                                ,NWTPG                                     !<-- # of Write tasks (servers) per group 
!
      INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE,SAVE :: NCURRENT_GROUP   !<-- The currently active write group
!
      CHARACTER(len=17) :: CONFIGFILE_01_NAME                              !<-- The name of domain #1's configure file
!
      TYPE(ESMF_Config),SAVE :: CF                                         !<-- The configure object for a domain
!
      TYPE(WRITE_INTERNAL_STATE),POINTER :: WRT_INT_STATE                  ! The internal state pointer.
!
      TYPE(ESMF_Config),SAVE :: CF_1                                       !<-- The uppermost domain's (#1's) configure file object
!
!-----------------------------------------------------------------------
!
      REAL(kind=KDBL) :: btim,btim0
!
      REAL(kind=KDBL),SAVE :: write_init_tim                            &
                             ,write_run_tim                             &
                             ,write_first_tim                           &
                             ,write_recv_outp_tim                       &
                             ,write_send_outp_tim                       &
                             ,write_get_fields_tim
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE WRITE_REGISTER(WRITE_COMP,RC_WRT)
! 
!-----------------------------------------------------------------------
!***  Register the Write component's Initialize, Run, and Finalize
!***  subroutine names.
!-----------------------------------------------------------------------
!
!***  HISTORY   
!       xx Feb 2007:  W. Yang  - Originator
!       30 Jun 2007:  T. Black - Modified to share same traits as
!                                rest of code.
!
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument variables
!------------------------
!
      TYPE(ESMF_GridComp) :: WRITE_COMP                                    ! The write component
!
      INTEGER,INTENT(OUT) :: RC_WRT                                        ! Final return code
!     
!---------------------
!***  Local variables
!---------------------
!
      INTEGER :: RC
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC_WRT=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set Entry Point for Initialize Step of Write Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompSetEntryPoint(WRITE_COMP                        &  !<-- The write component
                                     ,ESMF_METHOD_INITIALIZE            &  !<-- Predefined subroutine type (INIT)
                                     ,WRITE_INITIALIZE                  &  !<-- User's subroutineName
                                     ,phase=1                           &
                                     ,rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set Entry Point for Run Step of Write Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompSetEntryPoint(WRITE_COMP                        &  !<-- The write component
                                     ,ESMF_METHOD_RUN                   &  !<-- Predefined subroutine type (RUN)
                                     ,WRITE_RUN                         &  !<-- User's subroutineName
                                     ,phase=1                           &
                                     ,rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set Entry Point for Finalize Step of Write Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
     CALL ESMF_GridCompSetEntryPoint(WRITE_COMP                         &  !<-- The write component
                                    ,ESMF_METHOD_FINALIZE               &  !<-- Predefined subroutine type (FINALIZE)
                                    ,WRITE_FINALIZE                     &  !<-- User's subroutineName
                                     ,phase=1                           &
                                     ,rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!-----------------------------------------------------------------------
!
      IF(RC_WRT==ESMF_SUCCESS)THEN
!       WRITE(6,*)"PASS: Write_Register."
      ELSE
        WRITE(0,*)"FAIL: Write_Register."
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE WRITE_REGISTER
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE WRITE_INITIALIZE(WRITE_COMP                            &
                                 ,IMP_STATE_WRITE                       &
                                 ,EXP_STATE_WRITE                       &
                                 ,CLOCK                                 &
                                 ,RC_INIT)
!
!-----------------------------------------------------------------------
!***  Initialize the Write gridded component.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument variables
!------------------------
!
      TYPE(ESMF_State) :: IMP_STATE_WRITE                               &   !<-- The Write component import state
                         ,EXP_STATE_WRITE                                   !<-- The Write component export state
!
      TYPE(ESMF_GridComp) :: WRITE_COMP                                     !<-- The Write component
!
      TYPE(ESMF_Clock) :: CLOCK                                             !<-- The Write component Clock
!
      INTEGER,INTENT(OUT) :: RC_INIT
!
!---------------------
!***  Local variables
!---------------------
!
      INTEGER(kind=KINT) :: ID_DOMAIN,IEND,IM                           &
                           ,INTERCOMM_WRITE_GROUP,IONE                  &
                           ,JEND,JM,LB,LBND,MYPE                        &
                           ,N,NL,NUM_WORDS_TOT,NV,UB,UBND
!
      INTEGER(kind=KINT) :: IERR,ISTAT,RC
!
      INTEGER(kind=KINT),DIMENSION(MPI_STATUS_SIZE) :: JSTAT
!
      TYPE(WRITE_WRAP) :: WRAP
!
      TYPE(WRITE_INTERNAL_STATE),POINTER :: WRT_INT_STATE
!
      TYPE(ESMF_VM) :: VM
!
!----------------------------------------------------------------------- 
!*********************************************************************** 
!----------------------------------------------------------------------- 
!
      btim0=timef()
!
!-----------------------------------------------------------------------
!
      RC     =ESMF_SUCCESS
      RC_INIT=ESMF_SUCCESS
!
!----------------------------------------------------------------------- 
!***  Initialize the Write component timers.
!----------------------------------------------------------------------- 
!
      write_init_tim=0.
      write_run_tim=0.
      write_first_tim=0.
      write_recv_outp_tim=0.
      write_send_outp_tim=0.
      write_get_fields_tim=0.
!
!----------------------------------------------------------------------- 
!***  Allocate the Write component's internal state.
!----------------------------------------------------------------------- 
!
      ALLOCATE(WRT_INT_STATE,stat=RC)
!
!----------------------------------------------------------------------- 
!***  Attach the internal state to the Write component.
!----------------------------------------------------------------------- 
!
      wrap%WRITE_INT_STATE=>WRT_INT_STATE
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Attach the Write Component's Internal State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompSetInternalState(WRITE_COMP                     &  !<-- The write component
                                        ,WRAP                           &  !<-- Pointer to the internal state
                                        ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!----------------------------------------------------------------------- 
!***  Retrieve the local VM.
!----------------------------------------------------------------------- 
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Retrieve the Local VM"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_VMGetCurrent(vm=VM                                      &  !<-- The ESMF virtual machine for this group of tasks
                            ,rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!----------------------------------------------------------------------- 
!***  We must keep track of the currently active Write group.
!***  The value pertinent to the current domain is the one we
!***  are interested in now therefore we need to extract the
!***  current domain ID from this domain's configure file and
!***  save the ID in the Write component's internal state.
!----------------------------------------------------------------------- 
!
      IF(.NOT.ALLOCATED(NCURRENT_GROUP))THEN
!
        ALLOCATE(NCURRENT_GROUP(1:NUM_DOMAINS_TOTAL),stat=ISTAT)
        IF(ISTAT/=0)THEN
          WRITE(0,*)' Failed to allocate NCURRENT_GROUP in Write Initialize'
          WRITE(0,*)' Aborting!'
          CALL ESMF_Finalize(rc=RC,endflag=ESMF_END_ABORT)
        ENDIF
!
        DO N=1,NUM_DOMAINS_TOTAL
          NCURRENT_GROUP(N)=0
        ENDDO
      ENDIF
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The current domain's configure file object
                                  ,value =wrt_int_state%ID_DOMAIN       &  !<-- Extract the current domain's ID and save it
                                  ,label ='my_domain_id:'               &  !<-- The quantity's label in the configure file
                                  ,rc    =RC)
!
      ID_DOMAIN=wrt_int_state%ID_DOMAIN
!
!----------------------------------------------------------------------- 
!***  Extract the task IDs and the number of tasks present.
!----------------------------------------------------------------------- 
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Get MPI Task IDs and Count from VM"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_VMGet(vm      =VM                                       &  !<-- The local VM
                     ,localPet=wrt_int_state%MYPE                       &  !<-- My task ID
                     ,petCount=wrt_int_state%NTASKS                     &  !<-- Number of MPI tasks present in current group (fcst+write)
                     ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      MYPE=wrt_int_state%MYPE
      NTASKS=wrt_int_state%NTASKS
!
!-----------------------------------------------------------------------
!***  All tasks allocate buffer data arrays that will hold scalar/1-D
!***  history/restart data and will be used to Send/Recv that data
!***  between the forecast tasks that know it initially to the
!***  write tasks that obtain it from the forecast tasks for writing.
!***  Logical data buffers are also handled here.
!-----------------------------------------------------------------------
!
      IF(.NOT.ASSOCIATED(wrt_int_state%ALL_DATA_I1D))THEN
        ALLOCATE(wrt_int_state%ALL_DATA_I1D(MAX_LENGTH_I1D),stat=ISTAT)
      ENDIF
!
      IF(.NOT.ASSOCIATED(wrt_int_state%ALL_DATA_R1D))THEN
        ALLOCATE(wrt_int_state%ALL_DATA_R1D(MAX_LENGTH_R1D),stat=ISTAT)
      ENDIF
!
      IF(.NOT.ASSOCIATED(wrt_int_state%ALL_DATA_LOG))THEN
        ALLOCATE(wrt_int_state%ALL_DATA_LOG(MAX_LENGTH_LOG),stat=ISTAT)
      ENDIF
!
      IF(.NOT.ASSOCIATED(wrt_int_state%RST_ALL_DATA_I1D))THEN
        ALLOCATE(wrt_int_state%RST_ALL_DATA_I1D(MAX_LENGTH_I1D),stat=ISTAT)
      ENDIF
!
      IF(.NOT.ASSOCIATED(wrt_int_state%RST_ALL_DATA_R1D))THEN
        ALLOCATE(wrt_int_state%RST_ALL_DATA_R1D(MAX_LENGTH_R1D),stat=ISTAT)
      ENDIF
!
      IF(.NOT.ASSOCIATED(wrt_int_state%RST_ALL_DATA_LOG))THEN
        ALLOCATE(wrt_int_state%RST_ALL_DATA_LOG(MAX_LENGTH_LOG),stat=ISTAT)
      ENDIF
!
!-----------------------------------------------------------------------
!***  Allocate dimensions as 1-WORD arrays since ESMF needs
!***  contiguous data arrays for ESMF_Sends/ESMF_Recvs when
!***  those dimensions are transmitted to the write tasks.
!-----------------------------------------------------------------------
!
      IF(.NOT.ALLOCATED(wrt_int_state%IM))THEN
        ALLOCATE(wrt_int_state%IM(1)                                    &
                ,wrt_int_state%IDS(1)                                   &
                ,wrt_int_state%IDE(1)                                   &
                ,wrt_int_state%JM(1)                                    &
                ,wrt_int_state%JDS(1)                                   &
                ,wrt_int_state%JDE(1)                                   &
                ,wrt_int_state%LM(1))
      ENDIF
!
!-----------------------------------------------------------------------
!***  The number of Attributes (for scalars and 1-D arrays) and
!***  Fields (for gridded 2-D arrays) in the Write component's
!***  import state are not known a priori.
!
!***  Even though these counts are just scalar integers we must
!***  allocate their pointers to length 1 since they will be
!***  used in ESMF_Send/Recv which require them to be contiguous
!***  data arrays.
!-----------------------------------------------------------------------
!
      IF(.NOT.ASSOCIATED(wrt_int_state%NCOUNT_FIELDS))THEN
!
        ALLOCATE(wrt_int_state%NCOUNT_FIELDS(1),stat=ISTAT)
!
        ALLOCATE(wrt_int_state%KOUNT_I1D(1),stat=ISTAT)
        ALLOCATE(wrt_int_state%KOUNT_I2D(1),stat=ISTAT)
!
        ALLOCATE(wrt_int_state%KOUNT_R1D(1),stat=ISTAT)
        ALLOCATE(wrt_int_state%KOUNT_R2D(1),stat=ISTAT)
!
        ALLOCATE(wrt_int_state%KOUNT_LOG(1),stat=ISTAT)
!
        ALLOCATE(wrt_int_state%RST_NCOUNT_FIELDS(1),stat=ISTAT)
!
        ALLOCATE(wrt_int_state%RST_KOUNT_I1D(1),stat=ISTAT)
        ALLOCATE(wrt_int_state%RST_KOUNT_I2D(1),stat=ISTAT)
!
        ALLOCATE(wrt_int_state%RST_KOUNT_R1D(1),stat=ISTAT)
        ALLOCATE(wrt_int_state%RST_KOUNT_R2D(1),stat=ISTAT)
!
        ALLOCATE(wrt_int_state%RST_KOUNT_LOG(1),stat=ISTAT)
!
!-----------------------------------------------------------------------
!***  All integer quantities (as 1-D arrays) and 1-D and 2-D real
!***  quantities will be strung together in single arrays of
!***  each particular type.  We need to allocate the arrays that will
!***  hold the length of each of the quantities in these 'strings'
!***  as the 'strings' themselves.
!-----------------------------------------------------------------------
!
        ALLOCATE(wrt_int_state%LENGTH_DATA_I1D(100),stat=ISTAT)            !<-- Lengths of each individual 1-D integer array
        ALLOCATE(wrt_int_state%LENGTH_DATA_R1D(100),stat=ISTAT)            !<-- Lengths of each individual 1-D real array
        ALLOCATE(wrt_int_state%LENGTH_SUM_I1D(1),stat=ISTAT)               !<-- Length of string of data of ALL 1-D integer arrays
        ALLOCATE(wrt_int_state%LENGTH_SUM_R1D(1),stat=ISTAT)               !<-- Length of string of data of ALL 1-D real arrays
        ALLOCATE(wrt_int_state%LENGTH_SUM_LOG(1),stat=ISTAT)               !<-- Length of string of data of ALL logical variables
!
        ALLOCATE(wrt_int_state%RST_LENGTH_DATA_I1D(100),stat=ISTAT)        !<-- Lengths of each restart individual 1-D integer array
        ALLOCATE(wrt_int_state%RST_LENGTH_DATA_R1D(100),stat=ISTAT)        !<-- Lengths of each restart individual 1-D real array
        ALLOCATE(wrt_int_state%RST_LENGTH_SUM_I1D(1),stat=ISTAT)           !<-- Length of string of restart data of ALL 1-D integer arrays
        ALLOCATE(wrt_int_state%RST_LENGTH_SUM_R1D(1),stat=ISTAT)           !<-- Length of string of restart data of ALL 1-D real arrays
        ALLOCATE(wrt_int_state%RST_LENGTH_SUM_LOG(1),stat=ISTAT)           !<-- Length of string of restart data of ALL logical variables
!
      ENDIF
!
!-----------------------------------------------------------------------
!***  Retrieve information regarding output from the configuration file.
!***  Values include write tasks per group, the number of write groups,
!***  output channels, and names of disk files.
!-----------------------------------------------------------------------
!
      CALL GET_CONFIG_WRITE(WRITE_COMP,WRT_INT_STATE,RC)                   !<-- User's routine to extract configfile data
!
!-----------------------------------------------------------------------
!***  The MPI intracommunicator for the forecast and quilt tasks
!***  associated with this Write component is naturally valid only
!***  for the domain (i.e., the Domain component) in whose internal
!***  state the Write component lies.  A task can lie on more than
!***  one domain in 2-way nesting therefore save the current intracomm
!***  in this Write component's internal state for use whenever it is
!***  active.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Write Initialize: Get the Current Intercomm"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!!!   CALL ESMF_VMGet(vm             =VM                                &  !<-- The local VM
!!!                  ,mpiCommunicator=INTERCOMM_WRITE_GROUP             &  !<-- This Write groups's (component's) intercommunicator
!!!                  ,rc             =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Increment the value of the current Write group.
!-----------------------------------------------------------------------
!
      NCURRENT_GROUP(ID_DOMAIN)=NCURRENT_GROUP(ID_DOMAIN)+1
!
!-----------------------------------------------------------------------
!***  The forecast tasks enter WRITE_INITIALIZE for every Write group
!***  but quilt tasks enter it only once.  The lead forecast task
!***  broadcasts the value of the group counter so all the quilt
!***  tasks have the correct value.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Write_Initialize: Broadcast Current Write Group Number"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_VMBroadcast(VM                                          &  !<-- The local VM
                           ,NCURRENT_GROUP                              &  !<-- The array of current active write groups
                           ,NUM_DOMAINS_TOTAL                           &  !<-- # of elements in the array
                           ,0                                           &  !<-- Root sender is fcst task 0
                           ,rc=RC) 
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Now we can save this Write group's intercommunicator and the
!***  local intracommunicator between fcst tasks among themselves
!***  and quilt tasks among themselves.
!***  NOTE:  We must use the value of the intercommunicator that was
!***         just generated in SETUP_SERVERS called from DOMAIN_SETUP
!***         called from DOMAIN_INITIALIZE prior to WRITE_INITIALIZE
!***         being called from WRITE_INIT called by DOMAIN_INITIALIZE.
!-----------------------------------------------------------------------
!
      IF(.NOT.ALLOCATED(wrt_int_state%MPI_INTERCOMM_ARRAY))THEN
        ALLOCATE(wrt_int_state%MPI_INTERCOMM_ARRAY(1:MAX_GROUPS))
      ENDIF
!
      INTERCOMM_WRITE_GROUP=MPI_INTERCOMM_ARRAY(NCURRENT_GROUP(ID_DOMAIN)) !<-- The intercommunicator from SETUP_SERVERS
!
      wrt_int_state%MPI_INTERCOMM_ARRAY(NCURRENT_GROUP(ID_DOMAIN))=     &
                                                 INTERCOMM_WRITE_GROUP
!
      wrt_int_state%MPI_COMM_COMP=MPI_COMM_COMP
!
!-----------------------------------------------------------------------
!
      NWTPG=wrt_int_state%WRITE_TASKS_PER_GROUP
      wrt_int_state%LAST_FCST_TASK=NTASKS-NWTPG-1
      wrt_int_state%LEAD_WRITE_TASK=wrt_int_state%LAST_FCST_TASK+1
      wrt_int_state%LAST_WRITE_TASK=NTASKS-1
!
      LAST_FCST_TASK=wrt_int_state%LAST_FCST_TASK
      LEAD_WRITE_TASK=wrt_int_state%LEAD_WRITE_TASK
      LAST_WRITE_TASK=wrt_int_state%LAST_WRITE_TASK
!
!-----------------------------------------------------------------------
!***  The forecast tasks now reset the current group to zero since
!***  that value will cycle through each group during the Run step.
!***  The quilt tasks only know about their own Write group so
!***  their value for the current group can be saved since it will
!***  never change.
!-----------------------------------------------------------------------
!
      IF(MYPE<=LAST_FCST_TASK)THEN
!
        IF(NCURRENT_GROUP(ID_DOMAIN)==wrt_int_state%WRITE_GROUPS)THEN  
          NCURRENT_GROUP(ID_DOMAIN)=0                                      !<-- Reset this counter for the Run step
        ENDIF
!
      ENDIF
!
!-----------------------------------------------------------------------
!***  Allocate the pointers that hold the local limits 
!***  of all the forecast tasks' subdomains then fill
!***  them with the values inserted into the Write
!***  import state in POINT_DYNAMICS_OUTPUT.
!-----------------------------------------------------------------------
!
      IF(.NOT.ALLOCATED(wrt_int_state%LOCAL_ISTART))THEN
        ALLOCATE(wrt_int_state%LOCAL_ISTART(0:LAST_FCST_TASK),stat=ISTAT)  !<-- Local starting I for each fcst task's subdomain
        ALLOCATE(wrt_int_state%LOCAL_IEND  (0:LAST_FCST_TASK),stat=ISTAT)  !<-- Local ending I for each fcst task's subdomain
        ALLOCATE(wrt_int_state%LOCAL_JSTART(0:LAST_FCST_TASK),stat=ISTAT)  !<-- Local starting J for each fcst task's subdomain
        ALLOCATE(wrt_int_state%LOCAL_JEND  (0:LAST_FCST_TASK),stat=ISTAT)  !<-- Local ending J for each fcst task's subdomain
      ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="WRITE_INITIALIZE: Extract Local Domain Limits"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF(MYPE<=LAST_FCST_TASK)THEN
!
        CALL ESMF_AttributeGet(state    =IMP_STATE_WRITE                &  !<-- The Write component's import state
                              ,name     ='LOCAL_ISTART'                 &  !<-- Name of the Attribute to extract
                              ,valueList=wrt_int_state%LOCAL_ISTART     &  !<-- Extract local subdomain starting I's
                              ,rc=RC)
!
        CALL ESMF_AttributeGet(state    =IMP_STATE_WRITE                &  !<-- The Write component's import state
                              ,name     ='LOCAL_IEND'                   &  !<-- Name of the Attribute to extract
                              ,valueList=wrt_int_state%LOCAL_IEND       &  !<-- Extract local subdomain ending I's
                              ,rc=RC)
!
        CALL ESMF_AttributeGet(state    =IMP_STATE_WRITE                &  !<-- The Write component's import state
                              ,name     ='LOCAL_JSTART'                 &  !<-- Name of the Attribute to extract
                              ,valueList=wrt_int_state%LOCAL_JSTART     &  !<-- Extract local subdomain starting J's
                              ,rc=RC)
!
        CALL ESMF_AttributeGet(state    =IMP_STATE_WRITE                &  !<-- The Write component's import state
                              ,name     ='LOCAL_JEND'                   &  !<-- Name of the Attribute to extract
                              ,valueList=wrt_int_state%LOCAL_JEND       &  !<-- Extract local subdomain ending J's
                              ,rc=RC)
!
      ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="WRITE_INITIALIZE: Extract Full Domain Limits"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IONE = 1
!
      IF(MYPE<=LAST_FCST_TASK)THEN
!
        CALL ESMF_AttributeGet(state    =IMP_STATE_WRITE                &  !<-- The Write component's import state
                              ,name     ='IDS'                          &  !<-- Name of the Attribute to extract
                              ,valueList=wrt_int_state%IDS              &  !<-- Extract full subdomain starting I
                              ,rc=RC)
!
        CALL ESMF_AttributeGet(state    =IMP_STATE_WRITE                &  !<-- The Write component's import state
                              ,name     ='IDE'                          &  !<-- Name of the Attribute to extract
                              ,valueList=wrt_int_state%IDE              &  !<-- Extract full subdomain ending I
                              ,rc=RC)
!
        CALL ESMF_AttributeGet(state    =IMP_STATE_WRITE                &  !<-- The Write component's import state
                              ,name     ='JDS'                          &  !<-- Name of the Attribute to extract
                              ,valueList=wrt_int_state%JDS              &  !<-- Extract full subdomain starting J
                              ,rc=RC)
!
        CALL ESMF_AttributeGet(state    =IMP_STATE_WRITE                &  !<-- The Write component's import state
                              ,name     ='JDE'                          &  !<-- Name of the Attribute to extract
                              ,valueList=wrt_int_state%JDE              &  !<-- Extract full subdomain ending J
                              ,rc=RC)
!
        IDS=wrt_int_state%IDS(1)
        IDE=wrt_int_state%IDE(1)
        JDS=wrt_int_state%JDS(1)
        JDE=wrt_int_state%JDE(1)
!
      ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Boundary data is required in the restart file.  Because it must
!***  come from the Solver component via export/import states and
!***  because it is not on the ESMF computational grid, it must be
!***  handled as ESMF Attributes which will be 1-D data strings.
!***  They must be 1-D data strings because ESMF Attributes can have
!***  no more than 1 dimension.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  First, what are the domain dimensions?
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="IM,JM,LM from Config Object for Write Initialize"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The ESMF configure object
                                  ,value =wrt_int_state%IM(1)           &  !<-- I dimension of full domain
                                  ,label ='im:'                         &  !<-- Label in configure file
                                  ,rc    =RC)
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The ESMF configure object
                                  ,value =wrt_int_state%JM(1)           &  !<-- J dimension of full domain
                                  ,label ='jm:'                         &  !<-- Label in configure file
                                  ,rc    =RC)
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The ESMF configure object
                                  ,value =wrt_int_state%LM(1)           &  !<-- # of model layers
                                  ,label ='lm:'                         &  !<-- Label in configure file
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IM=wrt_int_state%IM(1)
      JM=wrt_int_state%JM(1)
      LM=wrt_int_state%LM(1)
!
!-----------------------------------------------------------------------
!***  # of words in BC data on each side of the domain on this
!***  forecast task.
!-----------------------------------------------------------------------
!
      fcst_tasks: IF(MYPE<=LAST_FCST_TASK)THEN
!
!-----------------------------------------------------------------------
!
        ALLOCATE(wrt_int_state%NUM_WORDS_BC_SOUTH(0:NUM_PES_FCST-1))
        ALLOCATE(wrt_int_state%NUM_WORDS_BC_NORTH(0:NUM_PES_FCST-1))
        ALLOCATE(wrt_int_state%NUM_WORDS_BC_WEST(0:NUM_PES_FCST-1))
        ALLOCATE(wrt_int_state%NUM_WORDS_BC_EAST(0:NUM_PES_FCST-1))
!
        DO N=0,NUM_PES_FCST-1
          wrt_int_state%NUM_WORDS_BC_SOUTH(N)=-1
          wrt_int_state%NUM_WORDS_BC_NORTH(N)=-1
          wrt_int_state%NUM_WORDS_BC_WEST(N) =-1
          wrt_int_state%NUM_WORDS_BC_EAST(N) =-1
        ENDDO
!
        ITS=wrt_int_state%LOCAL_ISTART(MYPE)
        ITE=wrt_int_state%LOCAL_IEND(MYPE)
        JTS=wrt_int_state%LOCAL_JSTART(MYPE)
        JTE=wrt_int_state%LOCAL_JEND(MYPE)
!
        IEND=ITE
        JEND=JTE
!
!-----------------------------------------------------------------------
!***  And how many boundary rows for H and V points?
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract LNSH from Write Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state=IMP_STATE_WRITE                    &  !<-- The Write component import state
                              ,name ='LNSH'                             &  !<-- Name of the Attribute to extract
                              ,value=wrt_int_state%LNSH                 &  !<-- Extract this Attribute from import state
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        LNSH=wrt_int_state%LNSH
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract LNSV from Write Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state=IMP_STATE_WRITE                    &  !<-- The Write component import state
                              ,name ='LNSV'                             &  !<-- Name of the Attribute to extract
                              ,value=wrt_int_state%LNSV                 &  !<-- Extract this Attribute from import state
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        LNSV=wrt_int_state%LNSV
!
!-----------------------------------------------------------------------
!***  How many 2-D,3-D,4-D H-pt and 2-D,3-D V-pt boundary variables?
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract BC Vbl Dim Counts from Write Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state=IMP_STATE_WRITE                    &  !<-- The Write component import state
                              ,name ='NVARS_BC_2D_H'                    &  !<-- Name of the Attribute to extract
                              ,value=wrt_int_state%NVARS_BC_2D_H        &  !<-- Extract this Attribute from import state
                              ,rc   =RC)
!
        NVARS_BC_2D_H=wrt_int_state%NVARS_BC_2D_H
!
        CALL ESMF_AttributeGet(state=IMP_STATE_WRITE                    &  !<-- The Write component import state
                              ,name ='NVARS_BC_3D_H'                    &  !<-- Name of the Attribute to extract
                              ,value=wrt_int_state%NVARS_BC_3D_H        &  !<-- Extract this Attribute from import state
                              ,rc   =RC)
!
        NVARS_BC_3D_H=wrt_int_state%NVARS_BC_3D_H
!
        CALL ESMF_AttributeGet(state=IMP_STATE_WRITE                    &  !<-- The Write component import state
                              ,name ='NVARS_BC_4D_H'                    &  !<-- Name of the Attribute to extract
                              ,value=wrt_int_state%NVARS_BC_4D_H        &  !<-- Extract this Attribute from import state
                              ,rc   =RC)
!
        NVARS_BC_4D_H=wrt_int_state%NVARS_BC_4D_H
!
        CALL ESMF_AttributeGet(state=IMP_STATE_WRITE                    &  !<-- The Write component import state
                              ,name ='NVARS_BC_2D_V'                    &  !<-- Name of the Attribute to extract
                              ,value=wrt_int_state%NVARS_BC_2D_V        &  !<-- Extract this Attribute from import state
                              ,rc   =RC)
!
        NVARS_BC_2D_V=wrt_int_state%NVARS_BC_2D_V
!
        CALL ESMF_AttributeGet(state=IMP_STATE_WRITE                    &  !<-- The Write component import state
                              ,name ='NVARS_BC_3D_V'                    &  !<-- Name of the Attribute to extract
                              ,value=wrt_int_state%NVARS_BC_3D_V        &  !<-- Extract this Attribute from import state
                              ,rc   =RC)
!
        NVARS_BC_3D_V=wrt_int_state%NVARS_BC_3D_V
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Lower/upper bounds of the count of the # of 3-D arrays within
!***  each 4-D variable.
!-----------------------------------------------------------------------
!
        IF(NVARS_BC_4D_H==0)THEN
          ALLOCATE(wrt_int_state%LBND_4D(1:1))
          ALLOCATE(wrt_int_state%UBND_4D(1:1))
!
!
        ELSEIF(NVARS_BC_4D_H>0)THEN
!
          ALLOCATE(wrt_int_state%LBND_4D(1:NVARS_BC_4D_H))
          ALLOCATE(wrt_int_state%UBND_4D(1:NVARS_BC_4D_H))
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Extract Lower Bnds of 3-D Arrays in 4-D Vbls"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeGet(state    =IMP_STATE_WRITE              &  !<-- The Write component import state
                                ,name     ='LBND_4D'                    &  !<-- Name of the Attribute to extract
                                ,valueList=wrt_int_state%LBND_4D        &  !<-- Lower bounds of 3-D array count in each 4-D vbl
                                ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Extract Upper Bnds of 3-D Arrays in 4-D Vbls"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeGet(state    =IMP_STATE_WRITE              &  !<-- The Write component import state
                                ,name     ='UBND_4D'                    &  !<-- Name of the Attribute to extract
                                ,valueList=wrt_int_state%UBND_4D        &  !<-- Upper bounds of 3-D array count in each 4-D vbl
                                ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        ENDIF
!
!-----------------------------------------------------------------------
!***  How many total levels in H-pt and V-pt boundary variables?
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract NLEV_H from Write Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state=IMP_STATE_WRITE                    &  !<-- The Write component import state
                              ,name ='NLEV_H'                           &  !<-- Name of the Attribute to extract
                              ,value=wrt_int_state%NLEV_H               &  !<-- Extract this Attribute from import state
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        NLEV_H=wrt_int_state%NLEV_H
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract NLEV_V from Write Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state=IMP_STATE_WRITE                    &  !<-- The Write component import state
                              ,name ='NLEV_V'                           &  !<-- Name of the Attribute to extract
                              ,value=wrt_int_state%NLEV_V               &  !<-- Extract this Attribute from import state
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        NLEV_V=wrt_int_state%NLEV_V
!
!-----------------------------------------------------------------------
!***  Compute the number of boundary words on boundary tasks then
!***  send the information to forecast task 0 who will ultimately
!***  assemble the full-domain boundary data.
!-----------------------------------------------------------------------
!
        IF(JTS==JDS)THEN                                                   !<-- Fcst tasks on south boundary
          wrt_int_state%NUM_WORDS_BC_SOUTH(MYPE)=(NLEV_H*LNSH           &
                                                 +NLEV_V*LNSV)          &
                                                 *2*(ITE-ITS+1)
          ALLOCATE(wrt_int_state%RST_BC_DATA_SOUTH(1:wrt_int_state%NUM_WORDS_BC_SOUTH(MYPE)),stat=ISTAT)
          IF(ISTAT/=0)WRITE(0,*)' Failed to allocate RST_BC_DATA_SOUTH'
!
          IF(MYPE/=0)THEN
            CALL MPI_SEND(wrt_int_state%NUM_WORDS_BC_SOUTH(MYPE)        &  !<-- Send this data
                         ,1                                             &  !<-- Number of words sent
                         ,MPI_INTEGER                                   &  !<-- Datatype
                         ,0                                             &  !<-- Send to fcst task 0
                         ,MYPE                                          &  !<-- An MPI tag
                         ,MPI_COMM_COMP                                 &  !<-- MPI communicator
                         ,IERR)
          ENDIF
!
        ENDIF
!
!-----------------------------------------------------------------------
!
        IF(JTE==JDE)THEN                                                   !<-- Fcst tasks on north boundary
          wrt_int_state%NUM_WORDS_BC_NORTH(MYPE)=(NLEV_H*LNSH           &
                                                 +NLEV_V*LNSV)          &
                                                 *2*(ITE-ITS+1)
          ALLOCATE(wrt_int_state%RST_BC_DATA_NORTH(1:wrt_int_state%NUM_WORDS_BC_NORTH(MYPE)),stat=ISTAT)
          IF(ISTAT/=0)WRITE(0,*)' Failed to allocate RST_BC_DATA_NORTH'
!
          IF(MYPE/=0)THEN
            CALL MPI_SEND(wrt_int_state%NUM_WORDS_BC_NORTH(MYPE)        &  !<-- Send this data
                         ,1                                             &  !<-- Number of words sent
                         ,MPI_INTEGER                                   &  !<-- Datatype
                         ,0                                             &  !<-- Send to fcst task 0
                         ,MYPE                                          &  !<-- An MPI tag
                         ,MPI_COMM_COMP                                 &  !<-- MPI communicator
                         ,IERR)
          ENDIF
!
        ENDIF
!
!-----------------------------------------------------------------------
!
        IF(ITS==IDS)THEN                                                   !<-- Fcst tasks on west boundary
          wrt_int_state%NUM_WORDS_BC_WEST(MYPE)=(NLEV_H*LNSH            &
                                                +NLEV_V*LNSV)           &
                                                *2*(JTE-JTS+1)
          ALLOCATE(wrt_int_state%RST_BC_DATA_WEST(1:wrt_int_state%NUM_WORDS_BC_WEST(MYPE)),stat=ISTAT)
          IF(ISTAT/=0)WRITE(0,*)' Failed to allocate RST_BC_DATA_WEST'
!
          IF(MYPE/=0)THEN
            CALL MPI_SEND(wrt_int_state%NUM_WORDS_BC_WEST(MYPE)         &  !<-- Send this data
                         ,1                                             &  !<-- Number of words sent
                         ,MPI_INTEGER                                   &  !<-- Datatype
                         ,0                                             &  !<-- Send to fcst task 0
                         ,MYPE                                          &  !<-- An MPI tag
                         ,MPI_COMM_COMP                                 &  !<-- MPI communicator
                         ,IERR)
          ENDIF
!
        ENDIF
!
!-----------------------------------------------------------------------
!
        IF(ITE==IDE)THEN                                                   !<-- Fcst tasks on east boundary
          wrt_int_state%NUM_WORDS_BC_EAST(MYPE)=(NLEV_H*LNSH            &
                                                +NLEV_V*LNSV)           &
                                                *2*(JTE-JTS+1)
          ALLOCATE(wrt_int_state%RST_BC_DATA_EAST(1:wrt_int_state%NUM_WORDS_BC_EAST(MYPE)),stat=ISTAT)
          IF(ISTAT/=0)WRITE(0,*)' Failed to allocate RST_BC_DATA_EAST'
!
          IF(MYPE/=0)THEN
            CALL MPI_SEND(wrt_int_state%NUM_WORDS_BC_EAST(MYPE)         &  !<-- Send this data
                         ,1                                             &  !<-- Number of words sent
                         ,MPI_INTEGER                                   &  !<-- Datatype
                         ,0                                             &  !<-- Send to fcst task 0
                         ,MYPE                                          &  !<-- An MPI tag
                         ,MPI_COMM_COMP                                 &  !<-- MPI communicator
                         ,IERR)
          ENDIF
!
        ENDIF
!
!-----------------------------------------------------------------------
!***    Forecast task 0 receives local BC word counts.
!-----------------------------------------------------------------------
!
        IF(MYPE==0.AND.NUM_PES_FCST>1)THEN
!
          DO N=1,LAST_FCST_TASK
!
            IF(wrt_int_state%LOCAL_JSTART(N)==1)THEN                       !<-- Recv from south bndry tasks
              CALL MPI_RECV(wrt_int_state%NUM_WORDS_BC_SOUTH(N)         &  !<-- Recv this data
                           ,1                                           &  !<-- Words received
                           ,MPI_INTEGER                                 &  !<-- Datatype
                           ,N                                           &  !<-- Recv from fcst task N
                           ,N                                           &  !<-- An MPI tag
                           ,MPI_COMM_COMP                               &  !<-- MPI communicator
                           ,JSTAT                                       &  !<-- MPI status object
                           ,IERR)
            ENDIF
!
            IF(wrt_int_state%LOCAL_JEND(N)==JDE)THEN                       !<-- Recv from north bndry tasks
              CALL MPI_RECV(wrt_int_state%NUM_WORDS_BC_NORTH(N)         &  !<-- Recv this data
                           ,1                                           &  !<-- Words received
                           ,MPI_INTEGER                                 &  !<-- Datatype
                           ,N                                           &  !<-- Recv from fcst task N
                           ,N                                           &  !<-- An MPI tag
                           ,MPI_COMM_COMP                               &  !<-- MPI communicator
                           ,JSTAT                                       &  !<-- MPI status object
                           ,IERR)
            ENDIF
!
            IF(wrt_int_state%LOCAL_ISTART(N)==1)THEN                       !<-- Recv from west bndry tasks
              CALL MPI_RECV(wrt_int_state%NUM_WORDS_BC_WEST(N)          &  !<-- Recv this data
                           ,1                                           &  !<-- Words received
                           ,MPI_INTEGER                                 &  !<-- Datatype
                           ,N                                           &  !<-- Recv from fcst task N
                           ,N                                           &  !<-- An MPI tag
                           ,MPI_COMM_COMP                               &  !<-- MPI communicator
                           ,JSTAT                                       &  !<-- MPI status object
                           ,IERR)
            ENDIF
!
            IF(wrt_int_state%LOCAL_IEND(N)==IDE)THEN                       !<-- Recv from east bndry tasks
              CALL MPI_RECV(wrt_int_state%NUM_WORDS_BC_EAST(N)          &  !<-- Recv this data
                           ,1                                           &  !<-- Words received
                           ,MPI_INTEGER                                 &  !<-- Datatype
                           ,N                                           &  !<-- Recv from fcst task N
                           ,N                                           &  !<-- An MPI tag
                           ,MPI_COMM_COMP                               &  !<-- MPI communicator
                           ,JSTAT                                       &  !<-- MPI status object
                           ,IERR)
            ENDIF
!
          ENDDO
!
        ENDIF
!
!-----------------------------------------------------------------------
!
      ENDIF fcst_tasks
!
!-----------------------------------------------------------------------
!***  Full-domain boundary data.
!-----------------------------------------------------------------------
!
      IF(MYPE==0)THEN                                                      !<-- Fcst task 0 needs all the full arrays
!
!-------------------------
!***  Each side of domain
!-------------------------
!
        IF(NVARS_BC_2D_H>0)THEN
!
          ALLOCATE(wrt_int_state%BND_VARS_H%VAR_2D(1:NVARS_BC_2D_H))
!
          DO NV=1,NVARS_BC_2D_H
            ALLOCATE(wrt_int_state%BND_VARS_H%VAR_2D(NV)%SOUTH(1:IDE,1:LNSH,1:2))
            ALLOCATE(wrt_int_state%BND_VARS_H%VAR_2D(NV)%NORTH(1:IDE,1:LNSH,1:2))
            ALLOCATE(wrt_int_state%BND_VARS_H%VAR_2D(NV)%WEST(1:LNSH,1:JDE,1:2))
            ALLOCATE(wrt_int_state%BND_VARS_H%VAR_2D(NV)%EAST(1:LNSH,1:JDE,1:2))
          ENDDO
        ENDIF
!
        IF(NVARS_BC_3D_H>0)THEN
!
          ALLOCATE(wrt_int_state%BND_VARS_H%VAR_3D(1:NVARS_BC_3D_H))
!
          DO NV=1,NVARS_BC_3D_H
            ALLOCATE(wrt_int_state%BND_VARS_H%VAR_3D(NV)%SOUTH(1:IDE,1:LNSH,1:LM,1:2))
            ALLOCATE(wrt_int_state%BND_VARS_H%VAR_3D(NV)%NORTH(1:IDE,1:LNSH,1:LM,1:2))
            ALLOCATE(wrt_int_state%BND_VARS_H%VAR_3D(NV)%WEST(1:LNSH,1:JDE,1:LM,1:2))
            ALLOCATE(wrt_int_state%BND_VARS_H%VAR_3D(NV)%EAST(1:LNSH,1:JDE,1:LM,1:2))
          ENDDO
        ENDIF
!
        IF(NVARS_BC_4D_H>0)THEN
!
          ALLOCATE(wrt_int_state%BND_VARS_H%VAR_4D(1:NVARS_BC_4D_H))
!
          DO NV=1,NVARS_BC_4D_H
            LB=wrt_int_state%LBND_4D(NV)
            UB=wrt_int_state%UBND_4D(NV)
            DO NL=LB,UB
              ALLOCATE(wrt_int_state%BND_VARS_H%VAR_4D(NV)%SOUTH(1:IDE,1:LNSH,1:LM,1:2,NL))
              ALLOCATE(wrt_int_state%BND_VARS_H%VAR_4D(NV)%NORTH(1:IDE,1:LNSH,1:LM,1:2,NL))
              ALLOCATE(wrt_int_state%BND_VARS_H%VAR_4D(NV)%WEST(1:LNSH,1:JDE,1:LM,1:2,NL))
              ALLOCATE(wrt_int_state%BND_VARS_H%VAR_4D(NV)%EAST(1:LNSH,1:JDE,1:LM,1:2,NL))
            ENDDO
          ENDDO
        ENDIF
!
        IF(NVARS_BC_2D_V>0)THEN
!
          ALLOCATE(wrt_int_state%BND_VARS_V%VAR_2D(1:NVARS_BC_2D_V))
!
          DO NV=1,NVARS_BC_2D_V
            ALLOCATE(wrt_int_state%BND_VARS_V%VAR_2D(NV)%SOUTH(1:IDE,1:LNSV,1:2))
            ALLOCATE(wrt_int_state%BND_VARS_V%VAR_2D(NV)%NORTH(1:IDE,1:LNSV,1:2))
            ALLOCATE(wrt_int_state%BND_VARS_V%VAR_2D(NV)%WEST(1:LNSV,1:JDE,1:2))
            ALLOCATE(wrt_int_state%BND_VARS_V%VAR_2D(NV)%EAST(1:LNSV,1:JDE,1:2))
          ENDDO
        ENDIF
!
        IF(NVARS_BC_3D_V>0)THEN
!
          ALLOCATE(wrt_int_state%BND_VARS_V%VAR_3D(1:NVARS_BC_3D_V))
!
          DO NV=1,NVARS_BC_3D_V
            ALLOCATE(wrt_int_state%BND_VARS_V%VAR_3D(NV)%SOUTH(1:IDE,1:LNSV,1:LM,1:2))
            ALLOCATE(wrt_int_state%BND_VARS_V%VAR_3D(NV)%NORTH(1:IDE,1:LNSV,1:LM,1:2))
            ALLOCATE(wrt_int_state%BND_VARS_V%VAR_3D(NV)%WEST(1:LNSV,1:JDE,1:LM,1:2))
            ALLOCATE(wrt_int_state%BND_VARS_V%VAR_3D(NV)%EAST(1:LNSV,1:JDE,1:LM,1:2))
          ENDDO
        ENDIF
!
!-----------------
!***  Full domain
!-----------------
!
        NUM_WORDS_TOT=(NLEV_H*LNSH                                      &  !<-- Total # of words
                      +NLEV_V*LNSV)                                     &  !    in full-domain
                      *2*2*(IDE-IDS+JDE-JDS+2)                             !<-- boundary arrays
!
        ALLOCATE(wrt_int_state%NUM_WORDS_SEND_BC(1))
        wrt_int_state%NUM_WORDS_SEND_BC(1)=NUM_WORDS_TOT
!
        ALLOCATE(wrt_int_state%RST_ALL_BC_DATA(1:NUM_WORDS_TOT),stat=ISTAT)
        IF(ISTAT/=0)THEN
          WRITE(0,*)' Failed to allocate RST_ALL_BC_DATA in Write Initialize'
          WRITE(0,*)' Aborting!'
          CALL ESMF_Finalize(rc=RC,endflag=ESMF_END_ABORT)
        ENDIF
!
      ENDIF
!
!-----------------------------------------------------------------------
!***  Set the IO_BaseTime to the initial Clock time.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set the Output Base Time to the Initial Clock Time"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ClockGet(clock    =CLOCK                                &  !<-- The ESMF Clock
                        ,startTime=wrt_int_state%IO_BASETIME            &  !<-- The Clock's starting time
                        ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Set the first history file's time index.
!-----------------------------------------------------------------------
!
      wrt_int_state%NFHOURS=0
!
!-----------------------------------------------------------------------
!
      IF(RC_INIT==ESMF_SUCCESS)THEN
!!!     WRITE(0,*)"PASS: Write_Initialize."
      ELSE
        WRITE(0,*)"FAIL: Write_Initialize."
      ENDIF
!
      write_init_tim=(timef()-btim0)
!
!----------------------------------------------------------------------- 
!
      END SUBROUTINE WRITE_INITIALIZE
!
!----------------------------------------------------------------------- 
!####################################################################### 
!----------------------------------------------------------------------- 
!
      SUBROUTINE WRITE_RUN(WRITE_COMP                                   &
                          ,IMP_STATE_WRITE                              &
                          ,EXP_STATE_WRITE                              &
                          ,CLOCK                                        &
                          ,RC_RUN)
!
!----------------------------------------------------------------------- 
!***  The Run step for the Write gridded component.  
!***  Move data intended for history output from the import state
!***  to the Write tasks.
!----------------------------------------------------------------------- 
!
      USE ESMF_FieldGetMOD
!
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument variables
!------------------------
!
      TYPE(ESMF_GridComp) :: WRITE_COMP                                    !<-- The Write component
!
      TYPE(ESMF_Clock) :: CLOCK                                            !<-- The Write component Clock
! 
      TYPE(ESMF_State) :: IMP_STATE_WRITE                               &  !<-- The Write component import state
                         ,EXP_STATE_WRITE                                  !<-- The Write component export state.
!
      INTEGER(kind=KINT),INTENT(OUT) :: RC_RUN 
!
!---------------------
!***  Local variables
!---------------------
!
      INTEGER(kind=KINT),SAVE :: IH_INT =MPI_REQUEST_NULL               &
                                ,IH_REAL=MPI_REQUEST_NULL
!
      INTEGER(kind=KINT),SAVE :: RST_IH_BC  =MPI_REQUEST_NULL           &
                                ,RST_IH_INT =MPI_REQUEST_NULL           &
                                ,RST_IH_REAL=MPI_REQUEST_NULL
!
      INTEGER(kind=KINT),SAVE :: IMS,IME,JMS,JME                        &
                                ,IHALO,JHALO                            &
                                ,INTERCOMM_WRITE_GROUP                  &
                                ,MPI_COMM_COMP                          &
                                ,N_START,N_END                          &
                                ,NPOSN_1,NPOSN_2                        &
                                ,NSUM_WORDS    
!
      INTEGER(kind=KINT),SAVE :: NSUM_WORDS_NEW=0
!
      INTEGER(kind=KINT) :: I,I1,IJ,IJG,IM,IMJM                         &
                           ,J,JM,K,KOUNT                                &
                           ,L,LB,LEAD_TASK_DOMAIN                       &
                           ,MY_LOCAL_ID,MY_RANK                         &
                           ,MYPE,MYPE_ROW                               &
                           ,N,N1,N2,NF,NN,NX                            &
                           ,NN_INTEGER,NN_REAL                          &
                           ,N_POSITION                                  &
                           ,NA,NB,NC,ND,NL                              &
                           ,NT,NTASK                                    &
                           ,NUM_ATTRIB,NV                               &
                           ,UB,WRITE_GROUP
!
      INTEGER(kind=KINT) :: DIM1                                        &
                           ,DIM2                                        &
                           ,FIELDSIZE                                   &
                           ,NBDR
!
      INTEGER(kind=KINT) :: IYEAR_FCST                                  &
                           ,IMONTH_FCST                                 &
                           ,IDAY_FCST                                   &
                           ,IHOUR_FCST                                  &
                           ,IMINUTE_FCST                                &
                           ,ISECOND_FCST                                &
                           ,ISECOND_NUM                                 &
                           ,ISECOND_DEN
!
      INTEGER(kind=KINT) :: NTIMESTEP
!
      INTEGER(kind=KINT) :: FRAC_SEC                                    &
                           ,INT_SEC                                     &
                           ,NF_HOURS                                    &
                           ,NF_MINUTES                                  &
                           ,NSECONDS                                    &
                           ,NSECONDS_NUM                                &
                           ,NSECONDS_DEN
!
      INTEGER(kind=KINT) :: ID_DUMMY                                    &
                           ,ID_RECV                                     &
                           ,ID_SEND                                     &
                           ,NFCST_TASKS                                 &
                           ,NFIELD                                      &
                           ,NPE_WRITE                                   &
                           ,NUM_FIELD_NAMES
!
      INTEGER(kind=KINT) :: ID_START,ID_END                             &
                           ,ISTART,IEND                                 &
                           ,JSTART,JEND
!
      INTEGER(kind=KINT) :: JROW_FIRST,JROW_LAST                        &
                           ,JSTA_WRITE,JEND_WRITE
!
      INTEGER(kind=KINT) :: KOUNT_I2D                                   &
                           ,KOUNT_I2D_DATA                              &
                           ,KOUNT_R2D                                   &
                           ,KOUNT_R2D_DATA
!
      INTEGER(kind=KINT) :: RST_KOUNT_I2D                               &
                           ,RST_KOUNT_I2D_DATA                          &
                           ,RST_KOUNT_R2D                               &
                           ,RST_KOUNT_R2D_DATA
!
      INTEGER(kind=KINT) :: ID_DOMAIN,LENGTH                            &
                           ,MPI_COMM,MPI_COMM2
!
      INTEGER(kind=KINT) :: IO_HST_UNIT,IO_RST_UNIT,FFSYNC
!
      INTEGER(kind=KINT) :: IERR,ISTAT,RC
!
!--posts
      INTEGER(kind=KINT) :: IEOF,NSOIL,POST_MAPTYPE
      CHARACTER(1)       :: POST_GRIDTYPE
      LOGICAL(kind=KLOG) :: LOG_PESET
      INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE :: JSTAGRP,JENDGRP
      INTEGER(kind=KINT) :: KPO,KTH,KPV
      real(kind=KFPT),dimension(70) :: PO,TH,PV
!--posts
!
      INTEGER(kind=ESMF_KIND_I8) :: NTIMESTEP_ESMF
!
      INTEGER(kind=KINT),DIMENSION(1) :: NCURRENT_GROUP_BCAST
!
      INTEGER(kind=KINT),DIMENSION(MPI_STATUS_SIZE) :: JSTAT
!
      INTEGER(kind=KINT),DIMENSION(:),POINTER :: WORK_ARRAY_I1D
!
      INTEGER(kind=KINT),DIMENSION(:,:),POINTER :: WORK_ARRAY_I2D
!
      REAL(kind=KFPT) :: NF_SECONDS                                     &
                        ,SECOND_FCST
!
      REAL(kind=KFPT),DIMENSION(:),POINTER :: GLAT1D,GLON1D,TMP         &
                                             ,WORK_ARRAY_R1D
!
      REAL(kind=KFPT),DIMENSION(:,:),POINTER  :: WORK_ARRAY_R2D
!
      REAL(KIND=KFPT),DIMENSION(:),ALLOCATABLE:: BUFF_NTASK             &
                                                ,HOLD_RST_DATA_R1D
!
      REAL(kind=KFPT),DIMENSION(:,:),ALLOCATABLE:: FACT10               &
                                                  ,FACT10TMPU           &
                                                  ,FACT10TMPV           &
                                                  ,HGT
!
      real(esmf_kind_i8) :: RUN_TIMESTEP_COUNT
!
      LOGICAL(kind=KLOG) :: GLOBAL                                      &
                           ,OPENED                                      &
                           ,WRITE_LOGICAL
!
      LOGICAL(kind=KLOG),SAVE :: FIRST=.TRUE.                           &
                                ,HST_FIRST=.TRUE.                       &
                                ,RST_FIRST=.TRUE.
!
      CHARACTER(3)           :: MODEL_LEVEL
      CHARACTER(ESMF_MAXSTR) :: FILENAME,GFNAME,NAME
!
      TYPE(WRITE_WRAP)                   :: WRAP
      TYPE(WRITE_INTERNAL_STATE),POINTER :: WRT_INT_STATE
!
      TYPE(ESMF_VM)        :: VM
      TYPE(ESMF_Grid),SAVE :: GRID1
      TYPE(ESMF_DELayout)  :: MY_DE_LAYOUT
      TYPE(ESMF_Field)     :: FIELD_WORK1
!
      TYPE(ESMF_Time)     :: CURRTIME
!
      TYPE(ESMF_TypeKind_Flag) :: DATATYPE
!
      TYPE(ESMF_LOGICAL),DIMENSION(:),POINTER :: FIRST_IO_PE
!
      TYPE(ESMF_FieldBundle) :: HISTORY_BUNDLE                          &  !<-- The history output data Bundle
                               ,RESTART_BUNDLE                             !<-- The restart output data Bundle
!
      TYPE(NEMSIO_GFILE)  :: NEMSIOFILE
!
!-----------------------------------------------------------------------
!
      real(kind=kfpt) :: wait_time
!
      integer(kind=kint),dimension(8) :: values
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
!-----------------------------------------------------------------------
!***  It is important to note that while the tasks executing this
!***  step include all forecast tasks plus those write tasks in
!***  this write group, the import state was filled in the Solver
!***  component only by the forecast tasks.  Therefore any information
!***  extracted from the import state that is needed by the write 
!***  tasks must be sent to them by the forecast tasks.
!
!***  Also note that history data consisting of scalars or 1D arrays
!***  are present in the Write component's import state as Attributes.
!***  All 2D (gridded) history data are present as Fields.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Retrieve the Write component's ESMF internal state.
!-----------------------------------------------------------------------
!
      btim=timef()
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Retrieve Write Component's Internal State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompGetInternalState(WRITE_COMP                     &  !<-- The write component
                                        ,WRAP                           &  !<-- Pointer to internal state
                                        ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      WRT_INT_STATE=>wrap%WRITE_INT_STATE                                  !<-- Local working pointer to internal state
!
!-----------------------------------------------------------------------
!***  Get the current local VM.
!***  This comes from the PetList used to create
!***  the Write components in DOMAIN_INITIALIZE.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Retrieve the Current VM for WRITE_RUN"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_VMGetCurrent(VM,rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      MYPE=wrt_int_state%MYPE
      ID_DOMAIN=wrt_int_state%ID_DOMAIN
!
!-----------------------------------------------------------------------
!
      IM=wrt_int_state%IM(1)
      JM=wrt_int_state%JM(1)
!
      KOUNT_I2D=wrt_int_state%KOUNT_I2D(1)
      KOUNT_R2D=wrt_int_state%KOUNT_R2D(1)
!
      LAST_FCST_TASK=wrt_int_state%LAST_FCST_TASK
      LEAD_WRITE_TASK=wrt_int_state%LEAD_WRITE_TASK
      LAST_WRITE_TASK=wrt_int_state%LAST_WRITE_TASK
!
      NWTPG=wrt_int_state%WRITE_TASKS_PER_GROUP
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Lead Fcst Task from Write Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=IMP_STATE_WRITE                      &
                            ,name ='Lead Task Domain'                   &
                            ,value=LEAD_TASK_DOMAIN                     &
                            ,rc   =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  The forecast tasks increment the write group so they know 
!***  which one is active.  Only the forecast tasks enter the 
!***  Run step of the Write component every output time therefore 
!***  only they need to increment the number of the current write
!***  group.  The quilt tasks belong to only a single write group
!***  so if they present here then it must be their write group 
!***  that is active.
!-----------------------------------------------------------------------
!
      IF(MYPE<=LAST_FCST_TASK)THEN
!
        NCURRENT_GROUP(ID_DOMAIN)=NCURRENT_GROUP(ID_DOMAIN)+1
        IF(NCURRENT_GROUP(ID_DOMAIN)>wrt_int_state%WRITE_GROUPS)THEN
          NCURRENT_GROUP(ID_DOMAIN)=1
        ENDIF
!
        NCURRENT_GROUP_BCAST(1)=NCURRENT_GROUP(ID_DOMAIN)
        NCURRENT_GROUP(ID_DOMAIN)=NCURRENT_GROUP_BCAST(1)
!
      ENDIF
!
!-----------------------------------------------------------------------
!***  Take the relevant intercommunicator for this Write group.
!-----------------------------------------------------------------------
!
      INTERCOMM_WRITE_GROUP=wrt_int_state%MPI_INTERCOMM_ARRAY(NCURRENT_GROUP(ID_DOMAIN))
!
!-----------------------------------------------------------------------
!***  Take the relevant intracommunicator for these fcst/quilt tasks.
!-----------------------------------------------------------------------
!
      MPI_COMM_COMP=wrt_int_state%MPI_COMM_COMP
!
!-----------------------------------------------------------------------
!
      write_first_tim=write_first_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  The elapsed forecast time (hours) will be appended to the name
!***  of each history output file.  Extract that value now.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Current Time for Output"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ClockGet(clock   =CLOCK                                 &  !<-- The ESMF Clock
                        ,currTime=CURRTIME                              &  !<-- The current time (ESMF) on the clock
                        ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  The ESMF time difference between start time and current time.
!-----------------------------------------------------------------------
!
      wrt_int_state%IO_CURRTIMEDIFF=CURRTIME-wrt_int_state%IO_BASETIME
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Elapsed Forecast Time for Output"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_TimeIntervalGet(timeinterval=wrt_int_state%IO_CURRTIMEDIFF &
                               ,h           =wrt_int_state%NFHOURS         &  !<-- The elapsed time in hours (REAL)
                               ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Now pull the 2d history data from the import state.
!***  This includes all individual 2D history quantities as well as
!***  all model levels of the 3D Real history arrays.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
      hst_fcst_tasks: IF(TIME_FOR_HISTORY.AND.MYPE<=LAST_FCST_TASK)THEN    !<-- Only the forecast tasks can see this data so far
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Extract the history data Bundle from the import state.
!***  The Bundle was created during the Init step of the SOLVER
!***  since subroutine POINT_OUTPUT must have it available for
!***  inserting data pointers into it.  Only the forecast tasks 
!***  can extract it properly since it was they who inserted it.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract History Bundle from Write Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_StateGet(state      =IMP_STATE_WRITE                  &  !<-- The write component's import state
                          ,itemName   ='History Bundle'                 &  !<-- The name of the history data Bundle
                          ,fieldbundle=HISTORY_BUNDLE                   &  !<-- The history data Bundle inside the import state
                          ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        NN_INTEGER=0
        NN_REAL   =0
!
        ITS=wrt_int_state%LOCAL_ISTART(MYPE)                               !<-- Starting I of this task's integration region
        ITE=wrt_int_state%LOCAL_IEND(MYPE)                                 !<-- Ending I of this task's integration region
        JTS=wrt_int_state%LOCAL_JSTART(MYPE)                               !<-- Starting J of this task's integration region
        JTE=wrt_int_state%LOCAL_JEND(MYPE)                                 !<-- Ending J of this task's integration region
!
        IHALO=wrt_int_state%IHALO                                          !<-- Halo depth in I
        JHALO=wrt_int_state%JHALO                                          !<-- Halo depth in J
!
        IDE=wrt_int_state%IDE(1)
        JDE=wrt_int_state%JDE(1)
!
!-----------------------------------------------------------------------
!***  Collect and send the updated Attributes (scalars and 1-D arrays)
!***  to the lead Write task for history output.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Get Attribute Count from History Bundle"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(FIELDBUNDLE=HISTORY_BUNDLE               &  !<-- The history data Bundle
                              ,count      =NUM_ATTRIB                   &  !<-- # of Attributes in the history Bundle
                              ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        IF(MYPE==0)THEN
          WRITE_GROUP=NCURRENT_GROUP(ID_DOMAIN)
          CALL SEND_UPDATED_ATTRIBUTES(HISTORY_BUNDLE                     &
                                      ,wrt_int_state%ALL_DATA_I1D         &
                                      ,wrt_int_state%ALL_DATA_R1D         &
                                      ,wrt_int_state%ALL_DATA_LOG         &
                                      ,MAX_LENGTH_I1D                     &
                                      ,MAX_LENGTH_R1D                     &
                                      ,MAX_LENGTH_LOG                     &
                                      ,MAX_GROUPS                         &
                                      ,WRITE_GROUP                        &
                                      ,INTERCOMM_WRITE_GROUP )
        ENDIF
!
!-----------------------------------------------------------------------
!***  Be sure the Integer and Real buffers are available for ISends.
!-----------------------------------------------------------------------
!
        btim=timef()
        CALL MPI_WAIT(IH_INT,JSTAT,IERR) 
        wait_time=(timef()-btim)
        if(wait_time>1.e3)write(0,*)' Long integer buffer WAIT =',wait_time*1.e-3
!     call date_and_time(values=values)
!     write(0,555)values(5),values(6),values(7),values(8)
  555 format(' Write Run after Wait IH_INT at ',i2.2,':',i2.2,':',i2.2,'.',i3.3)
!
        btim=timef()
        CALL MPI_WAIT(IH_REAL,JSTAT,IERR) 
        wait_time=(timef()-btim)
        if(wait_time>1.e3)write(0,*)' Long real buffer WAIT =',wait_time*1.e-3
!     call date_and_time(values=values)
!     write(0,556)values(5),values(6),values(7),values(8)
  556 format(' Write Run after Wait IH_REAL at ',i2.2,':',i2.2,':',i2.2,'.',i3.3)
!
!-----------------------------------------------------------------------
!
        btim=timef()
!
        field_block: DO N=1,wrt_int_state%NCOUNT_FIELDS(1)                 !<-- Loop through all Fields in the import state
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Extract 2-D Fields from History Bundle"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_FieldBundleGet(FIELDBUNDLE=HISTORY_BUNDLE                &  !<-- The write component's history data Bundle
                                  ,fieldName  =wrt_int_state%FIELD_NAME(N)   &  !<-- The ESMF Field's name
                                  ,field      =FIELD_WORK1                   &  !<-- The ESMF Field data pointer
                                  ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Does this extracted Field hold Integer or Real data?
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Check Datatype of Field from History Bundle"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_FieldGet(field   =FIELD_WORK1                       &  !<-- The ESMF Field
                            ,typekind=DATATYPE                          &  !<-- ESMF specifier of variable type and kind
                            ,rc      =RC)
 
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!--------------------------------------------------------------------
!                      -- INTEGER FIELDS --
!--------------------------------------------------------------------
!
          IF(DATATYPE==ESMF_TYPEKIND_I4)THEN                               !<-- Extract integer gridded data from each ESMF Field
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Extract Pointer from 2-D Integer Field"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

            CALL ESMF_FieldGet(field     =FIELD_WORK1                   &  !<-- The ESMF Field
                              ,localDe   =0                             &  !<-- # of DEs in this grid
                              ,farrayPtr =WORK_ARRAY_I2D                &  !<-- Put the 2D integer data from the Field here
                              ,rc        =RC)

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            ISTART=LBOUND(WORK_ARRAY_I2D,1)
            IEND  =UBOUND(WORK_ARRAY_I2D,1)
            JSTART=LBOUND(WORK_ARRAY_I2D,2)
            JEND  =UBOUND(WORK_ARRAY_I2D,2)
!
            IF(NN_INTEGER+(IEND-ISTART+1)*(JEND-JSTART+1)>wrt_int_state%NUM_WORDS_SEND_I2D_HST)THEN
              WRITE(0,*)' WARNING:  THE NUMBER OF INTEGER WORDS YOU'    &
                       ,' ARE SENDING FROM FCST TO WRITE TASKS HAS'     &
                       ,' EXCEEDED THE ORIGINAL COUNT WHICH SHOULD'     &
                       ,' NOT CHANGE.  CHECK YOUR WORK'
            ENDIF
!
            DO J=JSTART,JEND
            DO I=ISTART,IEND
              NN_INTEGER=NN_INTEGER+1
              wrt_int_state%ALL_DATA_I2D(NN_INTEGER)=WORK_ARRAY_I2D(I,J)   !<-- String together this task's 2D integer data
            ENDDO
            ENDDO
!
!--------------------------------------------------------------------
!                        -- REAL FIELDS --
!--------------------------------------------------------------------
!
          ELSEIF(DATATYPE==ESMF_TYPEKIND_R4)THEN                           !<-- Extract real gridded data from each ESMF Field
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Extract Pointer from 2-D Real Field"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

            CALL ESMF_FieldGet(field     =FIELD_WORK1                   &  !<-- The ESMF Field
                              ,localDe   =0                             &  !<-- # of DEs in this grid
                              ,farrayPtr =WORK_ARRAY_R2D                &  !<-- Put the 2D real data from the Field here
                              ,rc        =RC)

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            ISTART=LBOUND(WORK_ARRAY_R2D,1)
            IEND  =UBOUND(WORK_ARRAY_R2D,1)
            JSTART=LBOUND(WORK_ARRAY_R2D,2)
            JEND  =UBOUND(WORK_ARRAY_R2D,2)
!
            IF(NN_REAL+(IEND-ISTART+1)*(JEND-JSTART+1)>wrt_int_state%NUM_WORDS_SEND_R2D_HST)THEN
              WRITE(0,*)' WARNING:  THE NUMBER OF REAL WORDS YOU'       &
                       ,' ARE SENDING FROM FCST TO WRITE TASKS HAS'     &
                       ,' EXCEEDED THE ORIGINAL COUNT WHICH SHOULD'     &
                       ,' NOT CHANGE.  CHECK YOUR WORK'
            ENDIF
!
            DO J=JSTART,JEND
            DO I=ISTART,IEND
              NN_REAL=NN_REAL+1
              wrt_int_state%ALL_DATA_R2D(NN_REAL)=WORK_ARRAY_R2D(I,J)      !<-- String together this task's 2D real data
            ENDDO
            ENDDO
!
          ENDIF
!
!-----------------------------------------------------------------------
!
        ENDDO field_block
!
        write_get_fields_tim=write_get_fields_tim+(timef()-btim)
!
        btim=timef()
!
!-----------------------------------------------------------------------
!***  All forecast tasks now send their strings of 2D history data
!***  to the appropriate write tasks.
!-----------------------------------------------------------------------
!
        KOUNT_I2D_DATA=wrt_int_state%NUM_WORDS_SEND_I2D_HST                !<-- Total #of words in 2D integer history data on this fcst task
        KOUNT_R2D_DATA=wrt_int_state%NUM_WORDS_SEND_R2D_HST                !<-- Total #of words in 2D real history data on this fcst task
!
        MYPE_ROW=MYPE/wrt_int_state%INPES+1                                !<-- Each fcst task's row among all rows of fcst tasks
!
        DO N=1,NWTPG                                                       !<-- Loop through the write tasks in this group
          CALL PARA_RANGE(wrt_int_state%JNPES,NWTPG,N                   &  !<-- Find each write task's first and last rows of
                         ,JROW_FIRST,JROW_LAST)                            !<--   fcst tasks from which it will recv
!
          NPE_WRITE=N-1                                                    !<-- Consider the write task with this local ID
                                                                           !    beginning with 0
!
          IF(MYPE_ROW>=JROW_FIRST.AND.MYPE_ROW<=JROW_LAST)THEN             !<-- This fcst task associated with this write task
!
!-----------------------------------------------------------------------
!***  First the 2-D Integer data.
!-----------------------------------------------------------------------
!
            IF(KOUNT_I2D>0)THEN
              CALL MPI_ISSEND(wrt_int_state%ALL_DATA_I2D                &  !<-- Fcst tasks' string of 2D integer history data
                             ,KOUNT_I2D_DATA                            &  !<-- #of words in the data string
                             ,MPI_INTEGER                               &  !<-- The datatype
                             ,NPE_WRITE                                 &  !<-- The target write task
                             ,wrt_int_state%NFHOURS                     &  !<-- An MPI tag
                             ,INTERCOMM_WRITE_GROUP                     &  !<-- The MPI intercommunicator between fcst and quilt tasks
                             ,IH_INT                                    &  !<-- MPI communication request handle
                             ,IERR )
!
              IF(IERR/=0)WRITE(0,*)' ISend of integer data by fcst task 0 has failed.  IERR=',IERR
            ENDIF
!
!-----------------------------------------------------------------------
!***  Then the 2-D Real data.
!-----------------------------------------------------------------------
!
            IF(KOUNT_R2D>0)THEN
              CALL MPI_ISSEND(wrt_int_state%ALL_DATA_R2D                  &  !<-- Fcst tasks' string of 2D real history data
                             ,KOUNT_R2D_DATA                              &  !<-- #of words in the data string
                             ,MPI_REAL                                    &  !<-- The datatype
                             ,NPE_WRITE                                   &  !<-- The target write task
                             ,wrt_int_state%NFHOURS                       &  !<-- An MPI tag
                             ,INTERCOMM_WRITE_GROUP                       &  !<-- The MPI intercommunicator between fcst and quilt tasks
                             ,IH_REAL                                     &  !<-- MPI communication request handle
                             ,IERR )
!
              IF(IERR/=0)WRITE(0,*)' ISend of real data by fcst task 0 has failed.  IERR=',IERR
            ENDIF
!
          ENDIF
!
        ENDDO
!
        write_send_outp_tim=write_send_outp_tim+timef()-btim
!
!-----------------------------------------------------------------------
!
      ENDIF hst_fcst_tasks
!
!-----------------------------------------------------------------------
!***  Now pull the 2D restart data from the import state.
!***  This includes all individual 2D restart quantities as well as
!***  all model levels of the 3D Real restart arrays.
!-----------------------------------------------------------------------
!
      RST_KOUNT_I2D=wrt_int_state%RST_KOUNT_I2D(1)
      RST_KOUNT_R2D=wrt_int_state%RST_KOUNT_R2D(1)
!
!-----------------------------------------------------------------------
      rst_fcst_tasks: IF(TIME_FOR_RESTART.AND.MYPE<=LAST_FCST_TASK)THEN    !<-- Only the forecast tasks can see this data so far
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Extract the restart data Bundle from the import state.
!***  The Bundle was created during the Init step of the Solver 
!***  since subroutine POINT_OUTPUT must have it available for
!***  inserting data pointers into it.  Only the forecast tasks
!***  can extract it properly since it was they who inserted it.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract Restart Bundle from Write Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_StateGet(state      =IMP_STATE_WRITE                  &  !<-- The write component's import state
                          ,itemName   ='Restart Bundle'                 &  !<-- The name of the restart data Bundle
                          ,fieldbundle=RESTART_BUNDLE                   &  !<-- The restart data Bundle inside the import state
                          ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        NN_INTEGER=0
        NN_REAL   =0
!
        ITS=wrt_int_state%LOCAL_ISTART(MYPE)                               !<-- Starting I of this task's integration region
        ITE=wrt_int_state%LOCAL_IEND(MYPE)                                 !<-- Ending I of this task's integration region
        JTS=wrt_int_state%LOCAL_JSTART(MYPE)                               !<-- Starting J of this task's integration region
        JTE=wrt_int_state%LOCAL_JEND(MYPE)                                 !<-- Ending J of this task's integration region
!
        IHALO=wrt_int_state%IHALO                                          !<-- Halo depth in I
        JHALO=wrt_int_state%JHALO                                          !<-- Halo depth in J
!
        IDE=wrt_int_state%IDE(1)
        JDE=wrt_int_state%JDE(1)
!
!-----------------------------------------------------------------------
!***  Collect and send the updated Attributes (scalars and 1-D arrays)
!***  to the lead Write task for restart output.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Get Attribute Count from Restart Bundle"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(FIELDBUNDLE=RESTART_BUNDLE               &  !<-- The restart data Bundle
                              ,count      =NUM_ATTRIB                   &  !<-- # of Attributes in the restart Bundle
                              ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        IF(MYPE==0)THEN
          WRITE_GROUP=NCURRENT_GROUP(ID_DOMAIN)
          CALL SEND_UPDATED_ATTRIBUTES(RESTART_BUNDLE                     &
                                      ,wrt_int_state%RST_ALL_DATA_I1D     &
                                      ,wrt_int_state%RST_ALL_DATA_R1D     &
                                      ,wrt_int_state%RST_ALL_DATA_LOG     &
                                      ,MAX_LENGTH_I1D                     &
                                      ,MAX_LENGTH_R1D                     &
                                      ,MAX_LENGTH_LOG                     &
                                      ,MAX_GROUPS                         &
                                      ,WRITE_GROUP                        &
                                      ,INTERCOMM_WRITE_GROUP )
        ENDIF
!
!-----------------------------------------------------------------------
!***  Be sure the Integer and Real buffers are available for ISends.
!-----------------------------------------------------------------------
!
        btim=timef()
        CALL MPI_WAIT(RST_IH_INT,JSTAT,IERR) 
        wait_time=(timef()-btim)
        if(wait_time>1.e3)write(0,*)' Long integer buffer WAIT =',wait_time*1.e-3
!
        btim=timef()
        CALL MPI_WAIT(RST_IH_REAL,JSTAT,IERR) 
        wait_time=(timef()-btim)
        if(wait_time>1.e3)write(0,*)' Long real buffer WAIT =',wait_time*1.e-3
!
!-----------------------------------------------------------------------
!
        btim=timef()
!
        rst_field_block: DO N=1,wrt_int_state%RST_NCOUNT_FIELDS(1)         !<-- Loop through all Fields in the import state
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Extract 2-D Fields from Restart Bundle"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_FieldBundleGet(FIELDBUNDLE=RESTART_BUNDLE                  &  !<-- The write component's restart data Bundle
                                  ,fieldName  =wrt_int_state%RST_FIELD_NAME(N) &  !<-- The ESMF Field's name
                                  ,field      =FIELD_WORK1                     &  !<-- The ESMF Field data pointer
                                  ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Does this extracted Field hold Integer or Real restart data?
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Check Datatype of Field from Restart Bundle"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_FieldGet(field   =FIELD_WORK1                       &  !<-- The ESMF Field
                            ,typekind=DATATYPE                          &  !<-- ESMF specifier of variable type and kind
                            ,rc      =RC)
 
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!--------------------------------------------------------------------
!                      -- INTEGER FIELDS --
!--------------------------------------------------------------------
!
          IF(DATATYPE==ESMF_TYPEKIND_I4)THEN                               !<-- Extract integer gridded data from each ESMF Field
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Extract Pointer from 2-D Integer Field"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

            CALL ESMF_FieldGet(field     =FIELD_WORK1                   &  !<-- The ESMF Field
                              ,localDe   =0                             &  !<-- # of DEs in this grid
                              ,farrayPtr =WORK_ARRAY_I2D                &  !<-- Put the 2D integer data from the Field here
                              ,rc        =RC)

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            ISTART=LBOUND(WORK_ARRAY_I2D,1)
            IEND  =UBOUND(WORK_ARRAY_I2D,1)
            JSTART=LBOUND(WORK_ARRAY_I2D,2)
            JEND  =UBOUND(WORK_ARRAY_I2D,2)
!
            IF(NN_INTEGER+(IEND-ISTART+1)*(JEND-JSTART+1)>wrt_int_state%NUM_WORDS_SEND_I2D_RST)THEN
              WRITE(0,*)' WARNING:  THE NUMBER OF INTEGER WORDS YOU'    &
                       ,' ARE SENDING FROM FCST TO WRITE TASKS HAS'     &
                       ,' EXCEEDED THE ORIGINAL COUNT WHICH SHOULD'     &
                       ,' NOT CHANGE.  CHECK YOUR WORK'
            ENDIF
!
            DO J=JSTART,JEND
            DO I=ISTART,IEND
              NN_INTEGER=NN_INTEGER+1
              wrt_int_state%RST_ALL_DATA_I2D(NN_INTEGER)=WORK_ARRAY_I2D(I,J)   !<-- String together this task's 2D integer data
            ENDDO
            ENDDO
!
!--------------------------------------------------------------------
!                        -- REAL FIELDS --
!--------------------------------------------------------------------
!
          ELSEIF(DATATYPE==ESMF_TYPEKIND_R4)THEN                           !<-- Extract real gridded data from each ESMF Field
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Extract Pointer from 2-D Real Field"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

            CALL ESMF_FieldGet(field     =FIELD_WORK1                   &  !<-- The ESMF Field
                              ,localDe   =0                             &  !<-- # of DEs in this grid
                              ,farrayPtr =WORK_ARRAY_R2D                &  !<-- Put the 2D real data from the Field here
                              ,rc        =RC)

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            RST_KOUNT_R2D=RST_KOUNT_R2D+1                                  !<-- Count # of 2D real Fields (as opposed to 2D integer Fields)
!
            ISTART=LBOUND(WORK_ARRAY_R2D,1)
            IEND  =UBOUND(WORK_ARRAY_R2D,1)
            JSTART=LBOUND(WORK_ARRAY_R2D,2)
            JEND  =UBOUND(WORK_ARRAY_R2D,2)
!
            IF(NN_REAL+(IEND-ISTART+1)*(JEND-JSTART+1)>wrt_int_state%NUM_WORDS_SEND_R2D_RST)THEN
              WRITE(0,*)' WARNING:  THE NUMBER OF REAL WORDS YOU'       &
                       ,' ARE SENDING FROM FCST TO WRITE TASKS HAS'     &
                       ,' EXCEEDED THE ORIGINAL COUNT WHICH SHOULD'     &
                       ,' NOT CHANGE.  CHECK YOUR WORK'
            ENDIF
!
            DO J=JSTART,JEND
            DO I=ISTART,IEND
              NN_REAL=NN_REAL+1
              wrt_int_state%RST_ALL_DATA_R2D(NN_REAL)=WORK_ARRAY_R2D(I,J)  !<-- String together this task's 2D real data
            ENDDO
            ENDDO
!
          ENDIF
!
!-----------------------------------------------------------------------
!
        ENDDO rst_field_block
!
        write_get_fields_tim=write_get_fields_tim+(timef()-btim)
!
        btim=timef()
!
!-----------------------------------------------------------------------
!***  All forecast tasks now send their strings of 2D restart data
!***  to the appropriate write tasks.
!-----------------------------------------------------------------------
!
        RST_KOUNT_I2D_DATA=wrt_int_state%NUM_WORDS_SEND_I2D_RST            !<-- # of words in 2D integer restart data on this fcst task
        RST_KOUNT_R2D_DATA=wrt_int_state%NUM_WORDS_SEND_R2D_RST            !<-- # of words in 2D real restart data on this fcst task
!
        MYPE_ROW=MYPE/wrt_int_state%INPES+1                                !<-- Each fcst task's row among all rows of fcst tasks
!
        DO N=1,NWTPG                                                       !<-- Loop through the write tasks in this group
          CALL PARA_RANGE(wrt_int_state%JNPES,NWTPG,N                   &  !<-- Find each write task's first and last rows of
                         ,JROW_FIRST,JROW_LAST)                            !<--   fcst tasks from which it will recv
!
          NPE_WRITE=N-1                                                    !<-- Consider the write task with this local ID
                                                                           !    beginning with 0
!
          IF(MYPE_ROW>=JROW_FIRST.AND.MYPE_ROW<=JROW_LAST)THEN             !<-- This fcst task associated with this write task
!
!-----------------------------------------------------------------------
!***  First the 2-D Integer restart data.
!-----------------------------------------------------------------------
!
            IF(RST_KOUNT_I2D>0)THEN
              CALL MPI_ISSEND(wrt_int_state%RST_ALL_DATA_I2D            &  !<-- Fcst tasks' string of 2D integer restart data
                             ,RST_KOUNT_I2D_DATA                        &  !<-- # of words in the data string
                             ,MPI_INTEGER                               &  !<-- The datatype
                             ,NPE_WRITE                                 &  !<-- The target write task
                             ,wrt_int_state%NFHOURS                     &  !<-- An MPI tag
                             ,INTERCOMM_WRITE_GROUP                     &  !<-- The MPI intercommunicator between fcst and quilt tasks
                             ,RST_IH_INT                                &  !<-- MPI communication request handle
                             ,IERR )
!
              IF(IERR/=0)WRITE(0,*)' ISend of integer data by fcst task 0 has failed.  IERR=',IERR
            ENDIF
!
!-----------------------------------------------------------------------
!***  Then the 2-D Real restart data.
!-----------------------------------------------------------------------
!
            IF(RST_KOUNT_R2D>0)THEN
              CALL MPI_ISSEND(wrt_int_state%RST_ALL_DATA_R2D              &  !<-- Fcst tasks' string of 2D real restart data
                             ,RST_KOUNT_R2D_DATA                          &  !<-- # of words in the data string
                             ,MPI_REAL                                    &  !<-- The datatype
                             ,NPE_WRITE                                   &  !<-- The target write task
                             ,wrt_int_state%NFHOURS                       &  !<-- An MPI tag
                             ,INTERCOMM_WRITE_GROUP                       &  !<-- The MPI intercommunicator between fcst and quilt tasks
                             ,RST_IH_REAL                                 &  !<-- MPI communication request handle
                             ,IERR )
!
              IF(IERR/=0)WRITE(0,*)' ISend of real data by fcst task 0 has failed.  IERR=',IERR
            ENDIF
!
          ENDIF
!
        ENDDO
!
        write_send_outp_tim=write_send_outp_tim+timef()-btim
!
!-----------------------------------------------------------------------
!***  The restart files need to contain the full-domain BC data
!***  in order for nests to produce bitwise identical results when
!***  restarting as compared to their free forecasts.
!***  Thus all boundary forecast tasks need to unload their pieces 
!***  of the data from the Write import state and send them to task 0 
!***  for assembly.
!-----------------------------------------------------------------------
!
!--------------------------
!***  South Boundary Tasks
!--------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract Local Pieces of South BC Wind Data"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        IF(wrt_int_state%NUM_WORDS_BC_SOUTH(MYPE)>0)THEN                   !<-- Fcst tasks along south boundary
!
          CALL ESMF_AttributeGet(state    =IMP_STATE_WRITE                  &  !<-- The Write component import state
                                ,name     ='RST_BC_DATA_SOUTH'              &  !<-- Name of south BC data on this task
                                ,valueList=wrt_int_state%RST_BC_DATA_SOUTH  &  !<-- Place the data here
                                ,rc       =RC)
!
          IF(MYPE/=0)THEN
            CALL MPI_SEND(wrt_int_state%RST_BC_DATA_SOUTH               &  !<-- Send this string of subdomain data 
                         ,wrt_int_state%NUM_WORDS_BC_SOUTH(MYPE)        &  !<-- Number of words sent
                         ,MPI_REAL                                      &  !<-- Datatype
                         ,0                                             &  !<-- Send the data to fcst task 0
                         ,MYPE                                          &  !<-- An MPI tag
                         ,MPI_COMM_COMP                                 &  !<-- MPI communicator
                         ,IERR )
          ENDIF
!
        ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!--------------------------
!***  North Boundary Tasks
!--------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Unload Local Pieces of North BC Wind Data"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        IF(wrt_int_state%NUM_WORDS_BC_NORTH(MYPE)>0)THEN                   !<-- Fcst tasks along north boundary
!
          CALL ESMF_AttributeGet(state    =IMP_STATE_WRITE                  &  !<-- The Write component import state
                                ,name     ='RST_BC_DATA_NORTH'              &  !<-- Name of north BC data on this task
                                ,valueList=wrt_int_state%RST_BC_DATA_NORTH  &  !<-- Place the data here
                                ,rc       =RC)
!
          IF(MYPE/=0)THEN
            CALL MPI_SEND(wrt_int_state%RST_BC_DATA_NORTH               &  !<-- Send this string of subdomain data 
                         ,wrt_int_state%NUM_WORDS_BC_NORTH(MYPE)        &  !<-- Number of words sent
                         ,MPI_REAL                                      &  !<-- Datatype
                         ,0                                             &  !<-- Send the data to the fcst task 0
                         ,MYPE                                          &  !<-- An MPI tag
                         ,MPI_COMM_COMP                                 &  !<-- MPI communicator
                         ,IERR )
          ENDIF
!
        ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-------------------------
!***  West Boundary Tasks
!-------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract Local Pieces of West BC Wind Data"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        IF(wrt_int_state%NUM_WORDS_BC_WEST(MYPE)>0)THEN                    !<-- Fcst tasks along west boundary
!
          CALL ESMF_AttributeGet(state    =IMP_STATE_WRITE                 &  !<-- The Write component import state
                                ,name     ='RST_BC_DATA_WEST'              &  !<-- Name of west BC data on this task
                                ,valueList=wrt_int_state%RST_BC_DATA_WEST  &  !<-- Place the data here
                                ,rc       =RC)
!
          IF(MYPE/=0)THEN
            CALL MPI_SEND(wrt_int_state%RST_BC_DATA_WEST                &  !<-- Send this string of subdomain data 
                         ,wrt_int_state%NUM_WORDS_BC_WEST(MYPE)         &  !<-- Number of words sent
                         ,MPI_REAL                                      &  !<-- Datatype
                         ,0                                             &  !<-- Send the data to the fcst task 0
                         ,MYPE                                          &  !<-- An MPI tag
                         ,MPI_COMM_COMP                                 &  !<-- MPI communicator
                         ,IERR )
          ENDIF
!
        ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-------------------------
!***  East Boundary Tasks
!-------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract Local Pieces of East BC Wind Data"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        IF(wrt_int_state%NUM_WORDS_BC_EAST(MYPE)>0)THEN                    !<-- Fcst tasks along east boundary
!
          CALL ESMF_AttributeGet(state    =IMP_STATE_WRITE                 &  !<-- The Write component import state
                                ,name     ='RST_BC_DATA_EAST'              &  !<-- Name of west BC data on this task
                                ,valueList=wrt_int_state%RST_BC_DATA_EAST  &  !<-- Place the data here
                                ,rc       =RC)
!
          IF(MYPE/=0)THEN
            CALL MPI_SEND(wrt_int_state%RST_BC_DATA_EAST                &  !<-- Send this string of subdomain data 
                         ,wrt_int_state%NUM_WORDS_BC_EAST(MYPE)         &  !<-- Number of words sent
                         ,MPI_REAL                                      &  !<-- Datatype
                         ,0                                             &  !<-- Send the data to the fcst task 0
                         ,MYPE                                          &  !<-- An MPI tag
                         ,MPI_COMM_COMP                                 &  !<-- MPI communicator
                         ,IERR )
          ENDIF
!
        ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Forecast task 0 receives the boundary data from the
!***  boundary tasks and assembles it into full-domain arrays.
!-----------------------------------------------------------------------
!
        fcst_task0: IF(MYPE==0)THEN
!
!-----------------------------------------------------------------------
!
!------------------------------------
!***  Recv from South boundary tasks
!------------------------------------
!
          recv_south: IF(wrt_int_state%INPES>1)THEN
!
            N1=1
            N2=wrt_int_state%INPES-1
!
            DO NTASK=N1,N2                                                 !<-- Task IDs of south boundary tasks    
              ALLOCATE(BUFF_NTASK(1:wrt_int_state%NUM_WORDS_BC_SOUTH(NTASK)))
!
              CALL MPI_RECV(BUFF_NTASK                                  &  !<-- Fcst tasks' string of local BC data
                           ,wrt_int_state%NUM_WORDS_BC_SOUTH(NTASK)     &  !<-- # of integer words in the local BC data string
                           ,MPI_REAL                                    &  !<-- Datatype
                           ,NTASK                                       &  !<-- Recv from this fcst task
                           ,NTASK                                       &  !<-- An MPI tag
                           ,MPI_COMM_COMP                               &  !<-- The MPI intercommunicator
                           ,JSTAT                                       &  !<-- MPI status object
                           ,IERR )
!
              NX=0
!
!-----------------------------------------------------------------------
!***  Extract south boundary points for 2-D H-pt variables.
!-----------------------------------------------------------------------
!
              IF(NVARS_BC_2D_H>0)THEN
                DO NV=1,NVARS_BC_2D_H
                  DO NT=1,2
                  DO NB=1,LNSH
                  DO NA=wrt_int_state%LOCAL_ISTART(NTASK),wrt_int_state%LOCAL_IEND(NTASK)
                    NX=NX+1
                    wrt_int_state%BND_VARS_H%VAR_2D(NV)%SOUTH(NA,NB,NT)=BUFF_NTASK(NX)
                  ENDDO
                  ENDDO
                  ENDDO
                ENDDO 
              ENDIF
!
!-----------------------------------------------------------------------
!***  Extract south boundary points for 3-D H-pt variables.
!-----------------------------------------------------------------------
!
              IF(NVARS_BC_3D_H>0)THEN
                DO NV=1,NVARS_BC_3D_H
                  DO NT=1,2
                  DO NC=1,LM
                  DO NB=1,LNSH
                  DO NA=wrt_int_state%LOCAL_ISTART(NTASK),wrt_int_state%LOCAL_IEND(NTASK)
                    NX=NX+1
                    wrt_int_state%BND_VARS_H%VAR_3D(NV)%SOUTH(NA,NB,NC,NT)=BUFF_NTASK(NX)
                  ENDDO 
                  ENDDO 
                  ENDDO 
                  ENDDO 
                ENDDO
              ENDIF
!
!-----------------------------------------------------------------------
!***  Extract south boundary points for 4-D H-pt variables.
!-----------------------------------------------------------------------
!
              IF(NVARS_BC_4D_H>0)THEN
                DO NV=1,NVARS_BC_4D_H
                  LB=wrt_int_state%LBND_4D(NV)
                  UB=wrt_int_state%UBND_4D(NV)
                  DO NL=LB,UB
                  DO NT=1,2
                  DO NC=1,LM
                  DO NB=1,LNSH
                  DO NA=wrt_int_state%LOCAL_ISTART(NTASK),wrt_int_state%LOCAL_IEND(NTASK)
                    NX=NX+1
                    wrt_int_state%BND_VARS_H%VAR_4D(NV)%SOUTH(NA,NB,NC,NT,NL)=BUFF_NTASK(NX)
                  ENDDO 
                  ENDDO 
                  ENDDO 
                  ENDDO 
                  ENDDO 
                ENDDO
              ENDIF
!
!-----------------------------------------------------------------------
!***  Extract south boundary points for 2-D V-pt variables.
!-----------------------------------------------------------------------
!
              IF(NVARS_BC_2D_V>0)THEN
                DO NV=1,NVARS_BC_2D_V
                  DO NT=1,2
                  DO NB=1,LNSV
                  DO NA=wrt_int_state%LOCAL_ISTART(NTASK),wrt_int_state%LOCAL_IEND(NTASK)
                    NX=NX+1
                    wrt_int_state%BND_VARS_V%VAR_2D(NV)%SOUTH(NA,NB,NT)=BUFF_NTASK(NX)
                  ENDDO
                  ENDDO
                  ENDDO
                ENDDO 
              ENDIF
!
!-----------------------------------------------------------------------
!***  Extract south boundary points for 3-D V-pt variables.
!-----------------------------------------------------------------------
!
              IF(NVARS_BC_3D_V>0)THEN
                DO NV=1,NVARS_BC_3D_V
                  DO NT=1,2
                  DO NC=1,LM
                  DO NB=1,LNSV
                  DO NA=wrt_int_state%LOCAL_ISTART(NTASK),wrt_int_state%LOCAL_IEND(NTASK)
                    NX=NX+1
                    wrt_int_state%BND_VARS_V%VAR_3D(NV)%SOUTH(NA,NB,NC,NT)=BUFF_NTASK(NX)
                  ENDDO 
                  ENDDO 
                  ENDDO 
                  ENDDO 
                ENDDO
              ENDIF
!
              DEALLOCATE(BUFF_NTASK)
            ENDDO
!
          ENDIF recv_south
!
!-----------------------------------------------------------------------
!
!------------------------------------
!***  Recv from North boundary tasks
!------------------------------------
!
          recv_north: IF(LAST_FCST_TASK>0)THEN

            IF(wrt_int_state%JNPES>1)THEN
              N1=(wrt_int_state%JNPES-1)*wrt_int_state%INPES
              N2=N1+wrt_int_state%INPES-1
            ELSEIF(wrt_int_state%JNPES==1)THEN
              N1=1
              N2=LAST_FCST_TASK
            ENDIF
!
            DO NTASK=N1,N2                                                 !<-- Task IDs of north boundary tasks
              ALLOCATE(BUFF_NTASK(1:wrt_int_state%NUM_WORDS_BC_NORTH(NTASK)))
!
              CALL MPI_RECV(BUFF_NTASK                                  &  !<-- Fcst tasks' string of local BC data
                           ,wrt_int_state%NUM_WORDS_BC_NORTH(NTASK)     &  !<-- # of integer words in the local BC data string
                           ,MPI_REAL                                    &  !<-- Datatype
                           ,NTASK                                       &  !<-- Recv from this fcst task
                           ,NTASK                                       &  !<-- An MPI tag
                           ,MPI_COMM_COMP                               &  !<-- The MPI intercommunicator
                           ,JSTAT                                       &  !<-- MPI status object
                           ,IERR )
!
              NX=0
!
!-----------------------------------------------------------------------
!***  Extract north boundary points for 2-D H-pt variables.
!-----------------------------------------------------------------------
!
              IF(NVARS_BC_2D_H>0)THEN
                DO NV=1,NVARS_BC_2D_H
                  DO NT=1,2
                  DO NB=1,LNSH
                  DO NA=wrt_int_state%LOCAL_ISTART(NTASK),wrt_int_state%LOCAL_IEND(NTASK)
                    NX=NX+1
                    wrt_int_state%BND_VARS_H%VAR_2D(NV)%NORTH(NA,NB,NT)=BUFF_NTASK(NX)
                  ENDDO
                  ENDDO
                  ENDDO
                ENDDO
              ENDIF
!
!-----------------------------------------------------------------------
!***  Extract north boundary points for 3-D H-pt variables.
!-----------------------------------------------------------------------
!
              IF(NVARS_BC_3D_H>0)THEN
                DO NV=1,NVARS_BC_3D_H
                  DO NT=1,2
                  DO NC=1,LM
                  DO NB=1,LNSH
                  DO NA=wrt_int_state%LOCAL_ISTART(NTASK),wrt_int_state%LOCAL_IEND(NTASK)
                    NX=NX+1
                    wrt_int_state%BND_VARS_H%VAR_3D(NV)%NORTH(NA,NB,NC,NT)=BUFF_NTASK(NX)
                  ENDDO
                  ENDDO
                  ENDDO
                  ENDDO
                ENDDO
              ENDIF
!
!-----------------------------------------------------------------------
!***  Extract north boundary points for 4-D H-pt variables.
!-----------------------------------------------------------------------
!
              IF(NVARS_BC_4D_H>0)THEN
                DO NV=1,NVARS_BC_4D_H
                  LB=wrt_int_state%LBND_4D(NV)
                  UB=wrt_int_state%UBND_4D(NV)
                  DO NL=LB,UB
                  DO NT=1,2
                  DO NC=1,LM
                  DO NB=1,LNSH
                  DO NA=wrt_int_state%LOCAL_ISTART(NTASK),wrt_int_state%LOCAL_IEND(NTASK)
                    NX=NX+1
                    wrt_int_state%BND_VARS_H%VAR_4D(NV)%NORTH(NA,NB,NC,NT,NL)=BUFF_NTASK(NX)
                  ENDDO
                  ENDDO
                  ENDDO
                  ENDDO
                  ENDDO
                ENDDO
              ENDIF
!
!-----------------------------------------------------------------------
!***  Extract north boundary points for 2-D V-pt variables.
!-----------------------------------------------------------------------
!
              IF(NVARS_BC_2D_V>0)THEN
                DO NV=1,NVARS_BC_2D_V
                  DO NT=1,2
                  DO NB=1,LNSV
                  DO NA=wrt_int_state%LOCAL_ISTART(NTASK),wrt_int_state%LOCAL_IEND(NTASK)
                    NX=NX+1
                    wrt_int_state%BND_VARS_V%VAR_2D(NV)%NORTH(NA,NB,NT)=BUFF_NTASK(NX)
                  ENDDO
                  ENDDO
                  ENDDO
                ENDDO
              ENDIF
!
!-----------------------------------------------------------------------
!***  Extract north boundary points for 3-D V-pt variables.
!-----------------------------------------------------------------------
!
              IF(NVARS_BC_3D_V>0)THEN
                DO NV=1,NVARS_BC_3D_V
                  DO NT=1,2
                  DO NC=1,LM
                  DO NB=1,LNSV
                  DO NA=wrt_int_state%LOCAL_ISTART(NTASK),wrt_int_state%LOCAL_IEND(NTASK)
                    NX=NX+1
                    wrt_int_state%BND_VARS_V%VAR_3D(NV)%NORTH(NA,NB,NC,NT)=BUFF_NTASK(NX)
                  ENDDO
                  ENDDO
                  ENDDO
                  ENDDO
                ENDDO
              ENDIF
!
              DEALLOCATE(BUFF_NTASK)
            ENDDO
!
          ENDIF recv_north
!
!-----------------------------------------------------------------------
!
!-----------------------------------
!***  Recv from West boundary tasks
!-----------------------------------
!
          recv_west: IF(wrt_int_state%JNPES>1)THEN
!
            N1=wrt_int_state%INPES
            N2=(wrt_int_state%JNPES-1)*wrt_int_state%INPES
!
            DO NTASK=N1,N2,wrt_int_state%INPES                             !<-- Task IDs of west boundary tasks
              ALLOCATE(BUFF_NTASK(1:wrt_int_state%NUM_WORDS_BC_WEST(NTASK)))
!
              CALL MPI_RECV(BUFF_NTASK                                  &  !<-- Fcst tasks' string of local BC data
                           ,wrt_int_state%NUM_WORDS_BC_WEST(NTASK)      &  !<-- # of integer words in the local BC data string
                           ,MPI_REAL                                    &  !<-- Datatype
                           ,NTASK                                       &  !<-- Recv from this fcst task
                           ,NTASK                                       &  !<-- An MPI tag
                           ,MPI_COMM_COMP                               &  !<-- The MPI intercommunicator
                           ,JSTAT                                       &  !<-- MPI status object
                           ,IERR )
!
              NX=0
!
!-----------------------------------------------------------------------
!***  Extract west boundary points for 2-D H-pt variables.
!-----------------------------------------------------------------------
!
              IF(NVARS_BC_2D_H>0)THEN
                DO NV=1,NVARS_BC_2D_H
                  DO NT=1,2
                  DO NB=wrt_int_state%LOCAL_JSTART(NTASK),wrt_int_state%LOCAL_JEND(NTASK)
                  DO NA=1,LNSH
                    NX=NX+1
                    wrt_int_state%BND_VARS_H%VAR_2D(NV)%WEST(NA,NB,NT)=BUFF_NTASK(NX)
                  ENDDO
                  ENDDO
                  ENDDO
                ENDDO
              ENDIF
!
!-----------------------------------------------------------------------
!***  Extract west boundary points for 3-D H-pt variables.
!-----------------------------------------------------------------------
!
              IF(NVARS_BC_3D_H>0)THEN
                DO NV=1,NVARS_BC_3D_H
                  DO NT=1,2
                  DO NC=1,LM
                  DO NB=wrt_int_state%LOCAL_JSTART(NTASK),wrt_int_state%LOCAL_JEND(NTASK)
                  DO NA=1,LNSH
                    NX=NX+1
                    wrt_int_state%BND_VARS_H%VAR_3D(NV)%WEST(NA,NB,NC,NT)=BUFF_NTASK(NX)
                  ENDDO
                  ENDDO
                  ENDDO
                  ENDDO
                ENDDO
              ENDIF
!
!-----------------------------------------------------------------------
!***  Extract west boundary points for 4-D H-pt variables.
!-----------------------------------------------------------------------
!
              IF(NVARS_BC_4D_H>0)THEN
                DO NV=1,NVARS_BC_4D_H
                  LB=wrt_int_state%LBND_4D(NV)
                  UB=wrt_int_state%UBND_4D(NV)
                  DO NL=LB,UB
                  DO NT=1,2
                  DO NC=1,LM
                  DO NB=wrt_int_state%LOCAL_JSTART(NTASK),wrt_int_state%LOCAL_JEND(NTASK)
                  DO NA=1,LNSH
                    NX=NX+1
                    wrt_int_state%BND_VARS_H%VAR_4D(NV)%WEST(NA,NB,NC,NT,NL)=BUFF_NTASK(NX)
                  ENDDO
                  ENDDO
                  ENDDO
                  ENDDO
                  ENDDO
                ENDDO
              ENDIF
!
!-----------------------------------------------------------------------
!***  Extract west boundary points for 2-D V-pt variables.
!-----------------------------------------------------------------------
!
              IF(NVARS_BC_2D_V>0)THEN
                DO NV=1,NVARS_BC_2D_V
                  DO NT=1,2
                  DO NB=wrt_int_state%LOCAL_JSTART(NTASK),wrt_int_state%LOCAL_JEND(NTASK)
                  DO NA=1,LNSV
                    NX=NX+1
                    wrt_int_state%BND_VARS_V%VAR_2D(NV)%WEST(NA,NB,NT)=BUFF_NTASK(NX)
                  ENDDO
                  ENDDO
                  ENDDO
                ENDDO
              ENDIF
!
!-----------------------------------------------------------------------
!***  Extract west boundary points for 3-D V-pt variables.
!-----------------------------------------------------------------------
!
              IF(NVARS_BC_3D_V>0)THEN
                DO NV=1,NVARS_BC_3D_V
                  DO NT=1,2
                  DO NC=1,LM
                  DO NB=wrt_int_state%LOCAL_JSTART(NTASK),wrt_int_state%LOCAL_JEND(NTASK)
                  DO NA=1,LNSV
                    NX=NX+1
                    wrt_int_state%BND_VARS_V%VAR_3D(NV)%WEST(NA,NB,NC,NT)=BUFF_NTASK(NX)
                  ENDDO
                  ENDDO
                  ENDDO
                  ENDDO
                ENDDO
              ENDIF
!
              DEALLOCATE(BUFF_NTASK)
            ENDDO
!
          ENDIF recv_west
!
!-----------------------------------------------------------------------
!
!-----------------------------------
!***  Recv from East boundary tasks
!-----------------------------------
!
          recv_east: IF(LAST_FCST_TASK>0)THEN     
!
            IF(wrt_int_state%INPES>1)THEN
              N1=wrt_int_state%INPES-1
              N2=wrt_int_state%JNPES*wrt_int_state%INPES-1
            ELSEIF(wrt_int_state%INPES==1)THEN
              N1=1
              N2=LAST_FCST_TASK
            ENDIF
!
            DO NTASK=N1,N2,wrt_int_state%INPES                             !<-- Task IDs of east boundary tasks
              ALLOCATE(BUFF_NTASK(1:wrt_int_state%NUM_WORDS_BC_EAST(NTASK)))
!
              CALL MPI_RECV(BUFF_NTASK                                  &  !<-- Fcst tasks' string of local BC data
                           ,wrt_int_state%NUM_WORDS_BC_EAST(NTASK)      &  !<-- # of integer words in the local BC data string
                           ,MPI_REAL                                    &  !<-- Datatype
                           ,NTASK                                       &  !<-- Recv from this fcst task
                           ,NTASK                                       &  !<-- An MPI tag
                           ,MPI_COMM_COMP                               &  !<-- The MPI intercommunicator
                           ,JSTAT                                       &  !<-- MPI status object
                           ,IERR )
!
              NX=0
!
!-----------------------------------------------------------------------
!***  Extract east boundary points for 2-D H-pt variables.
!-----------------------------------------------------------------------
!
              IF(NVARS_BC_2D_H>0)THEN
                DO NV=1,NVARS_BC_2D_H
                  DO NT=1,2
                  DO NB=wrt_int_state%LOCAL_JSTART(NTASK),wrt_int_state%LOCAL_JEND(NTASK)
                  DO NA=1,LNSH
                    NX=NX+1
                    wrt_int_state%BND_VARS_H%VAR_2D(NV)%EAST(NA,NB,NT)=BUFF_NTASK(NX)
                  ENDDO
                  ENDDO
                  ENDDO
                ENDDO
              ENDIF
!
!-----------------------------------------------------------------------
!***  Extract east boundary points for 3-D H-pt variables.
!-----------------------------------------------------------------------
!
              IF(NVARS_BC_3D_H>0)THEN
                DO NV=1,NVARS_BC_3D_H
                  DO NT=1,2
                  DO NC=1,LM
                  DO NB=wrt_int_state%LOCAL_JSTART(NTASK),wrt_int_state%LOCAL_JEND(NTASK)
                  DO NA=1,LNSH
                    NX=NX+1
                    wrt_int_state%BND_VARS_H%VAR_3D(NV)%EAST(NA,NB,NC,NT)=BUFF_NTASK(NX)
                  ENDDO
                  ENDDO
                  ENDDO
                  ENDDO
                ENDDO
              ENDIF
!
!-----------------------------------------------------------------------
!***  Extract east boundary points for 4-D H-pt variables.
!-----------------------------------------------------------------------
!
              IF(NVARS_BC_4D_H>0)THEN
                DO NV=1,NVARS_BC_4D_H
                  LB=wrt_int_state%LBND_4D(NV)
                  UB=wrt_int_state%UBND_4D(NV)
                  DO NL=LB,UB
                  DO NT=1,2
                  DO NC=1,LM
                  DO NB=wrt_int_state%LOCAL_JSTART(NTASK),wrt_int_state%LOCAL_JEND(NTASK)
                  DO NA=1,LNSH
                    NX=NX+1
                    wrt_int_state%BND_VARS_H%VAR_4D(NV)%EAST(NA,NB,NC,NT,NL)=BUFF_NTASK(NX)
                  ENDDO
                  ENDDO
                  ENDDO
                  ENDDO
                  ENDDO
                ENDDO
              ENDIF
!
!-----------------------------------------------------------------------
!***  Extract east boundary points for 2-D V-pt variables.
!-----------------------------------------------------------------------
!
              IF(NVARS_BC_2D_V>0)THEN
                DO NV=1,NVARS_BC_2D_V
                  DO NT=1,2
                  DO NB=wrt_int_state%LOCAL_JSTART(NTASK),wrt_int_state%LOCAL_JEND(NTASK)
                  DO NA=1,LNSV
                    NX=NX+1
                    wrt_int_state%BND_VARS_V%VAR_2D(NV)%EAST(NA,NB,NT)=BUFF_NTASK(NX)
                  ENDDO
                  ENDDO
                  ENDDO
                ENDDO
              ENDIF
!
!-----------------------------------------------------------------------
!***  Extract east boundary points for 3-D V-pt variables.
!-----------------------------------------------------------------------
!
              IF(NVARS_BC_3D_V>0)THEN
                DO NV=1,NVARS_BC_3D_V
                  DO NT=1,2
                  DO NC=1,LM
                  DO NB=wrt_int_state%LOCAL_JSTART(NTASK),wrt_int_state%LOCAL_JEND(NTASK)
                  DO NA=1,LNSV
                    NX=NX+1
                    wrt_int_state%BND_VARS_V%VAR_3D(NV)%EAST(NA,NB,NC,NT)=BUFF_NTASK(NX)
                  ENDDO
                  ENDDO
                  ENDDO
                  ENDDO
                ENDDO
              ENDIF
!
              DEALLOCATE(BUFF_NTASK)
            ENDDO
!
          ENDIF recv_east
!
!-----------------------------------------------------------------------
!***  Forecast task 0's own BC data.  First the south boundary.
!-----------------------------------------------------------------------
!
          IF(JTS==JDS)THEN
!
            NX=0
!
!-----------------------------------------------------------------------
!***  Extract south boundary points for 2-D H-pt variables.
!-----------------------------------------------------------------------
!
            IF(NVARS_BC_2D_H>0)THEN
              DO NV=1,NVARS_BC_2D_H
                DO NT=1,2
                DO NB=1,LNSH
                DO NA=ITS,ITE
                  NX=NX+1
                  wrt_int_state%BND_VARS_H%VAR_2D(NV)%SOUTH(NA,NB,NT)=  &
                      wrt_int_state%RST_BC_DATA_SOUTH(NX)
                ENDDO
                ENDDO
                ENDDO
              ENDDO
            ENDIF
!
!-----------------------------------------------------------------------
!***  Extract south boundary points for 3-D H-pt variables.
!-----------------------------------------------------------------------
!
            IF(NVARS_BC_3D_H>0)THEN
              DO NV=1,NVARS_BC_3D_H
                DO NT=1,2
                DO NC=1,LM
                DO NB=1,LNSH
                DO NA=ITS,ITE
                  NX=NX+1
                  wrt_int_state%BND_VARS_H%VAR_3D(NV)%SOUTH(NA,NB,NC,NT)= &
                     wrt_int_state%RST_BC_DATA_SOUTH(NX)
                ENDDO
                ENDDO
                ENDDO
                ENDDO
              ENDDO
            ENDIF
!
!-----------------------------------------------------------------------
!***  Extract south boundary points for 4-D H-pt variables.
!-----------------------------------------------------------------------
!
            IF(NVARS_BC_4D_H>0)THEN
              DO NV=1,NVARS_BC_4D_H
                LB=wrt_int_state%LBND_4D(NV)
                UB=wrt_int_state%UBND_4D(NV)
                DO NL=LB,UB
                DO NT=1,2
                DO NC=1,LM
                DO NB=1,LNSH
                DO NA=ITS,ITE
                  NX=NX+1
                  wrt_int_state%BND_VARS_H%VAR_4D(NV)%SOUTH(NA,NB,NC,NT,NL)= &
                     wrt_int_state%RST_BC_DATA_SOUTH(NX)
                ENDDO
                ENDDO
                ENDDO
                ENDDO
                ENDDO
              ENDDO
            ENDIF
!
!-----------------------------------------------------------------------
!***  Extract south boundary points for 2-D V-pt variables.
!-----------------------------------------------------------------------
!
            IF(NVARS_BC_2D_V>0)THEN
              DO NV=1,NVARS_BC_2D_V
                DO NT=1,2
                DO NB=1,LNSV
                DO NA=ITS,ITE
                  NX=NX+1
                  wrt_int_state%BND_VARS_V%VAR_2D(NV)%SOUTH(NA,NB,NT)= &
                     wrt_int_state%RST_BC_DATA_SOUTH(NX)
                ENDDO
                ENDDO
                ENDDO
              ENDDO
            ENDIF
!
!-----------------------------------------------------------------------
!***  Extract south boundary points for 3-D V-pt variables.
!-----------------------------------------------------------------------
!
            IF(NVARS_BC_3D_V>0)THEN
              DO NV=1,NVARS_BC_3D_V
                DO NT=1,2
                DO NC=1,LM
                DO NB=1,LNSV
                DO NA=ITS,ITE
                  NX=NX+1
                  wrt_int_state%BND_VARS_V%VAR_3D(NV)%SOUTH(NA,NB,NC,NT)= &
                     wrt_int_state%RST_BC_DATA_SOUTH(NX)
                ENDDO
                ENDDO
                ENDDO
                ENDDO
              ENDDO
            ENDIF
!
          ENDIF
!
!-----------------------------------------------------------------------
!***  Forecast task 0's north boundary.
!-----------------------------------------------------------------------
!
          IF(JTE==JDE)THEN     
!
            NX=0
!
!-----------------------------------------------------------------------
!***  Extract north boundary points for 2-D H-pt variables.
!-----------------------------------------------------------------------
!
            IF(NVARS_BC_2D_H>0)THEN
              DO NV=1,NVARS_BC_2D_H
                DO NT=1,2
                DO NB=1,LNSH
                DO NA=ITS,ITE
                  NX=NX+1
                  wrt_int_state%BND_VARS_H%VAR_2D(NV)%NORTH(NA,NB,NT)=  &
                      wrt_int_state%RST_BC_DATA_NORTH(NX)
                ENDDO
                ENDDO
                ENDDO
              ENDDO
            ENDIF
!
!-----------------------------------------------------------------------
!***  Extract north boundary points for 3-D H-pt variables.
!-----------------------------------------------------------------------
!
            IF(NVARS_BC_3D_H>0)THEN
              DO NV=1,NVARS_BC_3D_H
                DO NT=1,2
                DO NC=1,LM
                DO NB=1,LNSH
                DO NA=ITS,ITE
                  NX=NX+1
                  wrt_int_state%BND_VARS_H%VAR_3D(NV)%NORTH(NA,NB,NC,NT)= &
                     wrt_int_state%RST_BC_DATA_NORTH(NX)
                ENDDO
                ENDDO
                ENDDO
                ENDDO
              ENDDO
            ENDIF
!
!-----------------------------------------------------------------------
!***  Extract north boundary points for 4-D H-pt variables.
!-----------------------------------------------------------------------
!
            IF(NVARS_BC_4D_H>0)THEN
              DO NV=1,NVARS_BC_4D_H
                LB=wrt_int_state%LBND_4D(NV)
                UB=wrt_int_state%UBND_4D(NV)
                DO NL=LB,UB
                DO NT=1,2
                DO NC=1,LM
                DO NB=1,LNSH
                DO NA=ITS,ITE
                  NX=NX+1
                  wrt_int_state%BND_VARS_H%VAR_4D(NV)%NORTH(NA,NB,NC,NT,NL)= &
                     wrt_int_state%RST_BC_DATA_NORTH(NX)
                ENDDO
                ENDDO
                ENDDO
                ENDDO
                ENDDO
              ENDDO
            ENDIF
!
!-----------------------------------------------------------------------
!***  Extract north boundary points for 2-D V-pt variables.
!-----------------------------------------------------------------------
!
            IF(NVARS_BC_2D_V>0)THEN
              DO NV=1,NVARS_BC_2D_V
                DO NT=1,2
                DO NB=1,LNSV
                DO NA=ITS,ITE
                  NX=NX+1
                  wrt_int_state%BND_VARS_V%VAR_2D(NV)%NORTH(NA,NB,NT)= &
                     wrt_int_state%RST_BC_DATA_NORTH(NX)
                ENDDO
                ENDDO
                ENDDO
              ENDDO
            ENDIF
!
!-----------------------------------------------------------------------
!***  Extract north boundary points for 3-D V-pt variables.
!-----------------------------------------------------------------------
!
            IF(NVARS_BC_3D_V>0)THEN
              DO NV=1,NVARS_BC_3D_V
                DO NT=1,2
                DO NC=1,LM
                DO NB=1,LNSV
                DO NA=ITS,ITE
                  NX=NX+1
                  wrt_int_state%BND_VARS_V%VAR_3D(NV)%NORTH(NA,NB,NC,NT)= &
                     wrt_int_state%RST_BC_DATA_NORTH(NX)
                ENDDO
                ENDDO
                ENDDO
                ENDDO
              ENDDO
            ENDIF
!
          ENDIF
!
!-----------------------------------------------------------------------
!***  Forecast task 0's west boundary.
!-----------------------------------------------------------------------
!
          IF(ITS==IDS)THEN  
!
            NX=0
!
!-----------------------------------------------------------------------
!***  Extract west boundary points for 2-D H-pt variables.
!-----------------------------------------------------------------------
!
            IF(NVARS_BC_2D_H>0)THEN
              DO NV=1,NVARS_BC_2D_H
                DO NT=1,2
                DO NB=JTS,JTE
                DO NA=1,LNSH
                  NX=NX+1
                  wrt_int_state%BND_VARS_H%VAR_2D(NV)%WEST(NA,NB,NT)=   &
                      wrt_int_state%RST_BC_DATA_WEST(NX)
                ENDDO
                ENDDO
                ENDDO
              ENDDO
            ENDIF
!
!-----------------------------------------------------------------------
!***  Extract west boundary points for 3-D H-pt variables.
!-----------------------------------------------------------------------
!
            IF(NVARS_BC_3D_H>0)THEN
              DO NV=1,NVARS_BC_3D_H
                DO NT=1,2
                DO NC=1,LM
                DO NB=JTS,JTE
                DO NA=1,LNSH
                  NX=NX+1
                  wrt_int_state%BND_VARS_H%VAR_3D(NV)%WEST(NA,NB,NC,NT)=  &
                     wrt_int_state%RST_BC_DATA_WEST(NX)
                ENDDO
                ENDDO
                ENDDO
                ENDDO
              ENDDO
            ENDIF
!
!-----------------------------------------------------------------------
!***  Extract west boundary points for 4-D H-pt variables.
!-----------------------------------------------------------------------
!
            IF(NVARS_BC_4D_H>0)THEN
              DO NV=1,NVARS_BC_3D_H
                LB=wrt_int_state%LBND_4D(NV)
                UB=wrt_int_state%UBND_4D(NV)
                DO NL=LB,UB
                DO NT=1,2
                DO NC=1,LM
                DO NB=JTS,JTE
                DO NA=1,LNSH
                  NX=NX+1
                  wrt_int_state%BND_VARS_H%VAR_4D(NV)%WEST(NA,NB,NC,NT,NL)=  &
                     wrt_int_state%RST_BC_DATA_WEST(NX)
                ENDDO
                ENDDO
                ENDDO
                ENDDO
                ENDDO
              ENDDO
            ENDIF
!
!-----------------------------------------------------------------------
!***  Extract west boundary points for 2-D V-pt variables.
!-----------------------------------------------------------------------
!
            IF(NVARS_BC_2D_V>0)THEN
              DO NV=1,NVARS_BC_2D_V
                DO NT=1,2
                DO NB=JTS,JTE
                DO NA=1,LNSV
                  NX=NX+1
                  wrt_int_state%BND_VARS_V%VAR_2D(NV)%WEST(NA,NB,NT)=   &
                      wrt_int_state%RST_BC_DATA_WEST(NX)
                ENDDO
                ENDDO
                ENDDO
              ENDDO
            ENDIF
!
!-----------------------------------------------------------------------
!***  Extract west boundary points for 3-D V-pt variables.
!-----------------------------------------------------------------------
!
            IF(NVARS_BC_3D_V>0)THEN
              DO NV=1,NVARS_BC_3D_V
                DO NT=1,2
                DO NC=1,LM
                DO NB=JTS,JTE
                DO NA=1,LNSV
                  NX=NX+1
                  wrt_int_state%BND_VARS_V%VAR_3D(NV)%WEST(NA,NB,NC,NT)=  &
                     wrt_int_state%RST_BC_DATA_WEST(NX)
                ENDDO
                ENDDO
                ENDDO
                ENDDO
              ENDDO
            ENDIF
!
          ENDIF
!
!-----------------------------------------------------------------------
!***  Forecast task 0's east boundary.
!-----------------------------------------------------------------------
!
          IF(ITE==IDE)THEN
!
            NX=0
!
!-----------------------------------------------------------------------
!***  Extract east boundary points for 2-D H-pt variables.
!-----------------------------------------------------------------------
!
            IF(NVARS_BC_2D_H>0)THEN
              DO NV=1,NVARS_BC_2D_H
                DO NT=1,2
                DO NB=JTS,JTE
                DO NA=1,LNSH
                  NX=NX+1
                  wrt_int_state%BND_VARS_H%VAR_2D(NV)%EAST(NA,NB,NT)=   &
                      wrt_int_state%RST_BC_DATA_EAST(NX)
                ENDDO
                ENDDO
                ENDDO
              ENDDO
            ENDIF
!
!-----------------------------------------------------------------------
!***  Extract east boundary points for 3-D H-pt variables.
!-----------------------------------------------------------------------
!
            IF(NVARS_BC_3D_H>0)THEN
              DO NV=1,NVARS_BC_3D_H
                DO NT=1,2
                DO NC=1,LM
                DO NB=JTS,JTE
                DO NA=1,LNSH
                  NX=NX+1
                  wrt_int_state%BND_VARS_H%VAR_3D(NV)%EAST(NA,NB,NC,NT)=  &
                     wrt_int_state%RST_BC_DATA_EAST(NX)
                ENDDO
                ENDDO
                ENDDO
                ENDDO
              ENDDO
            ENDIF
!
!-----------------------------------------------------------------------
!***  Extract east boundary points for 4-D H-pt variables.
!-----------------------------------------------------------------------
!
            IF(NVARS_BC_4D_H>0)THEN
              DO NV=1,NVARS_BC_3D_H
                LB=wrt_int_state%LBND_4D(NV)
                UB=wrt_int_state%UBND_4D(NV)
                DO NL=LB,UB
                DO NT=1,2
                DO NC=1,LM
                DO NB=JTS,JTE
                DO NA=1,LNSH
                  NX=NX+1
                  wrt_int_state%BND_VARS_H%VAR_4D(NV)%EAST(NA,NB,NC,NT,NL)=  &
                     wrt_int_state%RST_BC_DATA_EAST(NX)
                ENDDO
                ENDDO
                ENDDO
                ENDDO
                ENDDO
              ENDDO
            ENDIF
!
!-----------------------------------------------------------------------
!***  Extract east boundary points for 2-D V-pt variables.
!-----------------------------------------------------------------------
!
            IF(NVARS_BC_2D_V>0)THEN
              DO NV=1,NVARS_BC_2D_V
                DO NT=1,2
                DO NB=JTS,JTE
                DO NA=1,LNSV
                  NX=NX+1
                  wrt_int_state%BND_VARS_V%VAR_2D(NV)%EAST(NA,NB,NT)=   &
                      wrt_int_state%RST_BC_DATA_EAST(NX)
                ENDDO
                ENDDO
                ENDDO
              ENDDO
            ENDIF
!
!-----------------------------------------------------------------------
!***  Extract east boundary points for 3-D V-pt variables.
!-----------------------------------------------------------------------
!
            IF(NVARS_BC_3D_V>0)THEN
              DO NV=1,NVARS_BC_3D_V
                DO NT=1,2
                DO NC=1,LM
                DO NB=JTS,JTE
                DO NA=1,LNSV
                  NX=NX+1
                  wrt_int_state%BND_VARS_V%VAR_3D(NV)%EAST(NA,NB,NC,NT)=  &
                     wrt_int_state%RST_BC_DATA_EAST(NX)
                ENDDO
                ENDDO
                ENDDO
                ENDDO
              ENDDO
            ENDIF
!
          ENDIF
!
!-----------------------------------------------------------------------
!***  Forecast task 0 renders all the boundary data into 
!***  a single 1-D string to send to the lead write task.
!***  Be sure the BC data buffer for the ISend to the lead
!***  quilt task is clear.
!-----------------------------------------------------------------------
!
          btim=timef()
          CALL MPI_WAIT(RST_IH_BC,JSTAT,IERR)
          wait_time=timef()-btim
          if(wait_time>1.e3)write(0,*)' Long BC buffer WAIT =',wait_time*1.e-3
!
          KOUNT=0
!
!--------------------
!***  South Boundary
!--------------------
!
          IF(NVARS_BC_2D_H>0)THEN
            DO NV=1,NVARS_BC_2D_H
              DO NT=1,2
              DO NB=1,LNSH
              DO NA=1,IDE
                KOUNT=KOUNT+1
                wrt_int_state%RST_ALL_BC_DATA(KOUNT)=                   &
                  wrt_int_state%BND_VARS_H%VAR_2D(NV)%SOUTH(NA,NB,NT)
              ENDDO
              ENDDO
              ENDDO
            ENDDO
          ENDIF
!
          IF(NVARS_BC_3D_H>0)THEN
            DO NV=1,NVARS_BC_3D_H
              DO NT=1,2
              DO NC=1,LM
              DO NB=1,LNSH
              DO NA=1,IDE
                KOUNT=KOUNT+1
                wrt_int_state%RST_ALL_BC_DATA(KOUNT)=                   &
                  wrt_int_state%BND_VARS_H%VAR_3D(NV)%SOUTH(NA,NB,NC,NT)
              ENDDO
              ENDDO
              ENDDO
              ENDDO
            ENDDO
          ENDIF
!
          IF(NVARS_BC_4D_H>0)THEN
            DO NV=1,NVARS_BC_4D_H
              LB=wrt_int_state%LBND_4D(NV)
              UB=wrt_int_state%UBND_4D(NV)
              DO NL=LB,UB
              DO NT=1,2
              DO NC=1,LM
              DO NB=1,LNSH
              DO NA=1,IDE
                KOUNT=KOUNT+1
                wrt_int_state%RST_ALL_BC_DATA(KOUNT)=                   &
                  wrt_int_state%BND_VARS_H%VAR_4D(NV)%SOUTH(NA,NB,NC,NT,NL)
              ENDDO
              ENDDO
              ENDDO
              ENDDO
              ENDDO
            ENDDO
          ENDIF
!
          IF(NVARS_BC_2D_V>0)THEN
            DO NV=1,NVARS_BC_2D_V
              DO NT=1,2
              DO NB=1,LNSV
              DO NA=1,IDE
                KOUNT=KOUNT+1
                wrt_int_state%RST_ALL_BC_DATA(KOUNT)=                   &
                  wrt_int_state%BND_VARS_V%VAR_2D(NV)%SOUTH(NA,NB,NT)
              ENDDO
              ENDDO
              ENDDO
            ENDDO
          ENDIF
!
          IF(NVARS_BC_3D_V>0)THEN
            DO NV=1,NVARS_BC_3D_V
              DO NT=1,2
              DO NC=1,LM
              DO NB=1,LNSV
              DO NA=1,IDE
                KOUNT=KOUNT+1
                wrt_int_state%RST_ALL_BC_DATA(KOUNT)=                   &
                  wrt_int_state%BND_VARS_V%VAR_3D(NV)%SOUTH(NA,NB,NC,NT)
              ENDDO
              ENDDO
              ENDDO
              ENDDO
            ENDDO
          ENDIF
!
!--------------------
!***  North Boundary
!--------------------
!
          IF(NVARS_BC_2D_H>0)THEN
            DO NV=1,NVARS_BC_2D_H
              DO NT=1,2
              DO NB=1,LNSH
              DO NA=1,IDE
                KOUNT=KOUNT+1
                wrt_int_state%RST_ALL_BC_DATA(KOUNT)=                   &
                  wrt_int_state%BND_VARS_H%VAR_2D(NV)%NORTH(NA,NB,NT)
              ENDDO
              ENDDO
              ENDDO
            ENDDO
          ENDIF
!
          IF(NVARS_BC_3D_H>0)THEN
            DO NV=1,NVARS_BC_3D_H
              DO NT=1,2
              DO NC=1,LM
              DO NB=1,LNSH
              DO NA=1,IDE
                KOUNT=KOUNT+1
                wrt_int_state%RST_ALL_BC_DATA(KOUNT)=                   &
                  wrt_int_state%BND_VARS_H%VAR_3D(NV)%NORTH(NA,NB,NC,NT)
              ENDDO
              ENDDO
              ENDDO
              ENDDO
            ENDDO
          ENDIF
!
          IF(NVARS_BC_4D_H>0)THEN
            DO NV=1,NVARS_BC_4D_H
              LB=wrt_int_state%LBND_4D(NV)
              UB=wrt_int_state%UBND_4D(NV)
              DO NL=LB,UB
              DO NT=1,2
              DO NC=1,LM
              DO NB=1,LNSH
              DO NA=1,IDE
                KOUNT=KOUNT+1
                wrt_int_state%RST_ALL_BC_DATA(KOUNT)=                   &
                  wrt_int_state%BND_VARS_H%VAR_4D(NV)%NORTH(NA,NB,NC,NT,NL)
              ENDDO
              ENDDO
              ENDDO
              ENDDO
              ENDDO
            ENDDO
          ENDIF
!
          IF(NVARS_BC_2D_V>0)THEN
            DO NV=1,NVARS_BC_2D_V
              DO NT=1,2
              DO NB=1,LNSV
              DO NA=1,IDE
                KOUNT=KOUNT+1
                wrt_int_state%RST_ALL_BC_DATA(KOUNT)=                   &
                  wrt_int_state%BND_VARS_V%VAR_2D(NV)%NORTH(NA,NB,NT)
              ENDDO
              ENDDO
              ENDDO
            ENDDO
          ENDIF
!
          IF(NVARS_BC_3D_V>0)THEN
            DO NV=1,NVARS_BC_3D_V
              DO NT=1,2
              DO NC=1,LM
              DO NB=1,LNSV
              DO NA=1,IDE
                KOUNT=KOUNT+1
                wrt_int_state%RST_ALL_BC_DATA(KOUNT)=                   &
                  wrt_int_state%BND_VARS_V%VAR_3D(NV)%NORTH(NA,NB,NC,NT)
              ENDDO
              ENDDO
              ENDDO
              ENDDO
            ENDDO
          ENDIF
!
!-------------------
!***  West Boundary
!-------------------
!
          IF(NVARS_BC_2D_H>0)THEN
            DO NV=1,NVARS_BC_2D_H
              DO NT=1,2
              DO NB=1,JDE
              DO NA=1,LNSH
                KOUNT=KOUNT+1
                wrt_int_state%RST_ALL_BC_DATA(KOUNT)=                   &
                  wrt_int_state%BND_VARS_H%VAR_2D(NV)%WEST(NA,NB,NT)
              ENDDO
              ENDDO
              ENDDO
            ENDDO
          ENDIF
!
          IF(NVARS_BC_3D_H>0)THEN
            DO NV=1,NVARS_BC_3D_H
              DO NT=1,2
              DO NC=1,LM
              DO NB=1,JDE
              DO NA=1,LNSH
                KOUNT=KOUNT+1
                wrt_int_state%RST_ALL_BC_DATA(KOUNT)=                   &
                  wrt_int_state%BND_VARS_H%VAR_3D(NV)%WEST(NA,NB,NC,NT)
              ENDDO
              ENDDO
              ENDDO
              ENDDO
            ENDDO
          ENDIF
!
          IF(NVARS_BC_4D_H>0)THEN
            DO NV=1,NVARS_BC_4D_H
              LB=wrt_int_state%LBND_4D(NV)
              UB=wrt_int_state%UBND_4D(NV)
              DO NL=LB,UB
              DO NT=1,2
              DO NC=1,LM
              DO NB=1,JDE
              DO NA=1,LNSH
                KOUNT=KOUNT+1
                wrt_int_state%RST_ALL_BC_DATA(KOUNT)=                   &
                  wrt_int_state%BND_VARS_H%VAR_4D(NV)%WEST(NA,NB,NC,NT,NL)
              ENDDO
              ENDDO
              ENDDO
              ENDDO
              ENDDO
            ENDDO
          ENDIF
!
          IF(NVARS_BC_2D_V>0)THEN
            DO NV=1,NVARS_BC_2D_V
              DO NT=1,2
              DO NB=1,JDE
              DO NA=1,LNSV
                KOUNT=KOUNT+1
                wrt_int_state%RST_ALL_BC_DATA(KOUNT)=                   &
                  wrt_int_state%BND_VARS_V%VAR_2D(NV)%WEST(NA,NB,NT)
              ENDDO
              ENDDO
              ENDDO
            ENDDO
          ENDIF
!
          IF(NVARS_BC_3D_V>0)THEN
            DO NV=1,NVARS_BC_3D_V
              DO NT=1,2
              DO NC=1,LM
              DO NB=1,JDE
              DO NA=1,LNSV
                KOUNT=KOUNT+1
                wrt_int_state%RST_ALL_BC_DATA(KOUNT)=                   &
                  wrt_int_state%BND_VARS_V%VAR_3D(NV)%WEST(NA,NB,NC,NT)
              ENDDO
              ENDDO
              ENDDO
              ENDDO
            ENDDO
          ENDIF
!
!-------------------
!***  East Boundary
!-------------------
!
          IF(NVARS_BC_2D_H>0)THEN
            DO NV=1,NVARS_BC_2D_H
              DO NT=1,2
              DO NB=1,JDE
              DO NA=1,LNSH
                KOUNT=KOUNT+1
                wrt_int_state%RST_ALL_BC_DATA(KOUNT)=                   &
                  wrt_int_state%BND_VARS_H%VAR_2D(NV)%EAST(NA,NB,NT)
              ENDDO
              ENDDO
              ENDDO
            ENDDO
          ENDIF
!
          IF(NVARS_BC_3D_H>0)THEN
            DO NV=1,NVARS_BC_3D_H
              DO NT=1,2
              DO NC=1,LM
              DO NB=1,JDE
              DO NA=1,LNSH
                KOUNT=KOUNT+1
                wrt_int_state%RST_ALL_BC_DATA(KOUNT)=                   &
                  wrt_int_state%BND_VARS_H%VAR_3D(NV)%EAST(NA,NB,NC,NT)
              ENDDO
              ENDDO
              ENDDO
              ENDDO
            ENDDO
          ENDIF
!
          IF(NVARS_BC_4D_H>0)THEN
            DO NV=1,NVARS_BC_4D_H
              LB=wrt_int_state%LBND_4D(NV)
              UB=wrt_int_state%UBND_4D(NV)
              DO NL=LB,UB
              DO NT=1,2
              DO NC=1,LM
              DO NB=1,JDE
              DO NA=1,LNSH
                KOUNT=KOUNT+1
                wrt_int_state%RST_ALL_BC_DATA(KOUNT)=                   &
                  wrt_int_state%BND_VARS_H%VAR_4D(NV)%EAST(NA,NB,NC,NT,NL)
              ENDDO
              ENDDO
              ENDDO
              ENDDO
              ENDDO
            ENDDO
          ENDIF
!
          IF(NVARS_BC_2D_V>0)THEN
            DO NV=1,NVARS_BC_2D_V
              DO NT=1,2
              DO NB=1,JDE
              DO NA=1,LNSV
                KOUNT=KOUNT+1
                wrt_int_state%RST_ALL_BC_DATA(KOUNT)=                   &
                  wrt_int_state%BND_VARS_V%VAR_2D(NV)%EAST(NA,NB,NT)
              ENDDO
              ENDDO
              ENDDO
            ENDDO
          ENDIF
!
          IF(NVARS_BC_3D_V>0)THEN
            DO NV=1,NVARS_BC_3D_V
              DO NT=1,2
              DO NC=1,LM
              DO NB=1,JDE
              DO NA=1,LNSV
                KOUNT=KOUNT+1
                wrt_int_state%RST_ALL_BC_DATA(KOUNT)=                   &
                  wrt_int_state%BND_VARS_V%VAR_3D(NV)%EAST(NA,NB,NC,NT)
              ENDDO
              ENDDO
              ENDDO
              ENDDO
            ENDDO
          ENDIF
!
!-----------------------------------------------------------------------
!***  Now forecast task 0 must send the full-domain boundary data
!***  to the lead write task for inserting it into the restart file.
!-----------------------------------------------------------------------
!
          CALL MPI_ISSEND(wrt_int_state%RST_ALL_BC_DATA                 &  !<-- 1-D String of full domain BC data
                         ,wrt_int_state%NUM_WORDS_SEND_BC(1)            &  !<-- # of words in the BC data string
                         ,MPI_REAL                                      &  !<-- The datatype
                         ,0                                             &  !<-- Local ID of lead write task
                         ,wrt_int_state%NFHOURS                         &  !<-- An MPI tag
                         ,INTERCOMM_WRITE_GROUP                         &  !<-- The MPI intercommunicator between fcst and quilt tasks
                         ,RST_IH_BC                                     &  !<-- MPI communication request handle
                         ,IERR )
!
          IF(IERR/=0)WRITE(0,*)' ISend of BC data by fcst task 0 has failed.  IERR=',IERR
!
        ENDIF fcst_task0
!
!-----------------------------------------------------------------------
!
      ENDIF rst_fcst_tasks
!
      write_run_tim=write_run_tim+(timef()-btim0)
!
!-----------------------------------------------------------------------
!***  The forecast tasks are completely finished with history and
!***  restart output now so they will exit the routine and resume
!***  the integration.
!-----------------------------------------------------------------------
!
      IF(MYPE<=LAST_FCST_TASK)RETURN
!
!-----------------------------------------------------------------------
!
      btim0=timef()
!
      history_time: IF(TIME_FOR_HISTORY) THEN
!
!-----------------------------------------------------------------------
!***  The lead Write task receives the latest Attributes from the
!***  lead fcst task.
!-----------------------------------------------------------------------
!
      IF(MYPE==LEAD_WRITE_TASK)THEN
!
        WRITE_GROUP=NCURRENT_GROUP(ID_DOMAIN)
!
!--------------
!***  Integers
!--------------
!
        IF(wrt_int_state%LENGTH_SUM_I1D(1)>0)THEN
!
          CALL MPI_RECV(wrt_int_state%ALL_DATA_I1D                    &  !<-- Recv string of integer history Attributes
                       ,wrt_int_state%LENGTH_SUM_I1D(1)               &  !<-- Words received
                       ,MPI_INTEGER                                   &  !<-- Data is integer
                       ,0                                             &  !<-- Sending task (lead fcst task)
                       ,WRITE_GROUP                                   &  !<-- MPI tag
                       ,INTERCOMM_WRITE_GROUP                         &  !<-- MPI domain commumicator
                       ,JSTAT                                         &  !<-- MPI status object
                       ,IERR )
!
        ENDIF
!
!-----------
!***  Reals
!-----------
!
        IF(wrt_int_state%LENGTH_SUM_R1D(1)>0)THEN
!
          CALL MPI_RECV(wrt_int_state%ALL_DATA_R1D                    &  !<-- Recv string of real history Attributes
                       ,wrt_int_state%LENGTH_SUM_R1D(1)               &  !<-- Words received
                       ,MPI_REAL                                      &  !<-- Data is real
                       ,0                                             &  !<-- Sending task (lead fcst task)
                       ,WRITE_GROUP                                   &  !<-- MPI tag
                       ,INTERCOMM_WRITE_GROUP                         &  !<-- MPI domain commumicator
                       ,JSTAT                                         &  !<-- MPI status object
                       ,IERR )
!
        ENDIF
!
!--------------
!***    Logicals
!--------------
!
        IF(wrt_int_state%LENGTH_SUM_LOG(1)>0)THEN
!
          CALL MPI_RECV(wrt_int_state%ALL_DATA_LOG                    &  !<-- Recv string of logical history Attributes
                       ,wrt_int_state%LENGTH_SUM_LOG(1)               &  !<-- Words received
                       ,MPI_LOGICAL                                   &  !<-- Data is logical
                       ,0                                             &  !<-- Sending task (lead fcst task)
                       ,WRITE_GROUP                                   &  !<-- MPI tag
                       ,INTERCOMM_WRITE_GROUP                         &  !<-- MPI domain commumicator
                       ,JSTAT                                         &  !<-- MPI status object
                       ,IERR )
!
        ENDIF
!
      ENDIF
!
!-----------------------------------------------------------------------
!***  Each write task in the active write group receives the
!***  strings of 2D history data from the appropriate fcst tasks.
!-----------------------------------------------------------------------
!
      ID_START=wrt_int_state%ID_FTASK_RECV_STA(MYPE)                       !<-- First fcst task that sends to this write task
      ID_END  =wrt_int_state%ID_FTASK_RECV_END(MYPE)                       !<-- Last fcst task that sends to this write task
      NFCST_TASKS=ID_END-ID_START+1                                        !<-- Number of fcst tasks sending to this write task
!
!-----------------------------------------------------------------------
      hst_from_fcst_tasks: DO N=1,NFCST_TASKS                              !<-- Loop through fcst tasks sending to this write task
!-----------------------------------------------------------------------
!
        ID_RECV=ID_START+N-1
!
!-----------------------------------------------------------------------
!***  Receive 2-D Integer history data if there is any.
!-----------------------------------------------------------------------
!
        IF(KOUNT_I2D>0)THEN
          CALL MPI_RECV(wrt_int_state%ALL_DATA_I2D                      &  !<-- Fcst tasks' string of 2D integer history data
                       ,wrt_int_state%NUM_WORDS_RECV_I2D_HST(ID_RECV)   &  !<-- # of integer words in the data string
                       ,MPI_INTEGER                                     &  !<-- The datatype
                       ,ID_RECV                                         &  !<-- Recv from this fcst task
                       ,wrt_int_state%NFHOURS                           &  !<-- An MPI tag
                       ,INTERCOMM_WRITE_GROUP                           &  !<-- The MPI intercommunicator between fcst and quilt tasks
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR )
!
          IF(IERR/=0)WRITE(0,*)' Recv by write task from fcst task has failed.  IERR=',IERR
        ENDIF
!
!-----------------------------------------------------------------------
!***  Receive 2-D Real history data if there is any.
!-----------------------------------------------------------------------
!
        IF(KOUNT_R2D>0)THEN
          CALL MPI_RECV(wrt_int_state%ALL_DATA_R2D                      &  !<-- Fcst tasks' string of 2D real history data
                       ,wrt_int_state%NUM_WORDS_RECV_R2D_HST(ID_RECV)   &  !<-- # of real words in the data string
                       ,MPI_REAL                                        &  !<-- The datatype
                       ,ID_RECV                                         &  !<-- Recv from this fcst task
                       ,wrt_int_state%NFHOURS                           &  !<-- An MPI tag
                       ,INTERCOMM_WRITE_GROUP                           &  !<-- The MPI intercommunicator between fcst and quilt tasks
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR )
!
          IF(IERR/=0)WRITE(0,*)' Recv by write task from fcst task has failed.  IERR=',IERR
        ENDIF
        write_recv_outp_tim=write_recv_outp_tim+timef()-btim
!
!-----------------------------------------------------------------------
!***  Each write task needs to insert the pieces of the various
!***  2D history arrays received from the individual fcst tasks
!***  into arrays that span the write tasks' own subsection of
!***  the full 2D domain.  That subsection always spans the 
!***  entire east-west dimension of the full domain (since full
!***  rows of fcst tasks always send to write tasks, never 
!***  partial rows) and as much of the north-south dimension of
!***  the full domain as covered by those fcst tasks sending to
!***  a given write task.
!-----------------------------------------------------------------------
!
        ITS=wrt_int_state%LOCAL_ISTART(ID_RECV)                            !<-- Local domain integration limits of sending fcst task
        ITE=wrt_int_state%LOCAL_IEND(ID_RECV)                              !<--
        JTS=wrt_int_state%LOCAL_JSTART(ID_RECV)                            !<--
        JTE=wrt_int_state%LOCAL_JEND(ID_RECV)                              !<--
!
        IHALO=wrt_int_state%IHALO                                          !<-- Subdomain halo depth in I
        JHALO=wrt_int_state%JHALO                                          !<-- Subdomain halo depth in J
!
        IMS=ITS-IHALO
        IME=ITE+IHALO
        JMS=JTS-JHALO
        JME=JTE+JHALO
!
        NN=0
!
        DO NF=1,KOUNT_I2D                                                 !<-- Loop through all the 2D integer fields
!
          DO J=JMS,JME
          DO I=IMS,IME
            NN=NN+1
            IF(I<ITS.OR.I>ITE.OR.J<JTS.OR.J>JTE)CYCLE                           !<-- Exclude halo points
            wrt_int_state%WRITE_SUBSET_I(I,J,NF)=wrt_int_state%ALL_DATA_I2D(NN) !<-- Put data into write task's domain subsection
          ENDDO
          ENDDO
        ENDDO
!
        NN=0
!
        DO NF=1,KOUNT_R2D                                                 !<-- Loop through all the 2D real fields
!
          DO J=JMS,JME
          DO I=IMS,IME
            NN=NN+1
            IF(I<ITS.OR.I>ITE.OR.J<JTS.OR.J>JTE)CYCLE                           !<-- Exclude halo points
            wrt_int_state%WRITE_SUBSET_R(I,J,NF)=wrt_int_state%ALL_DATA_R2D(NN) !<-- Put data into write task's domain subsection
          ENDDO
          ENDDO
        ENDDO
!
!-----------------------------------------------------------------------
!
      ENDDO hst_from_fcst_tasks
!            
!-----------------------------------------------------------------------
!***  At this point, all write tasks have received all of the history
!***  data from their associated fcst tasks and assembled it onto 
!***  their own subsections of the full 2D domain.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  The lead write task obtains the current forecast time and
!***  the elapsed forecast time for its writing of output.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  For dopost, all write tasks need to know the time.
!-----------------------------------------------------------------------
!

       IF(wrt_int_state%WRITE_DOPOST) THEN
         LOG_PESET=MYPE>=LEAD_WRITE_TASK
       ELSE
         LOG_PESET=MYPE==LEAD_WRITE_TASK
       ENDIF
!
!-----------------------------------------------------------------------
      hst_time_get: IF(LOG_PESET)THEN                                   !<-- The lead write task
!-----------------------------------------------------------------------
!
        IF(wrt_int_state%WRITE_HST_BIN.OR.                              &
           wrt_int_state%WRITE_HST_NEMSIO)THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Lead Write Task Gets Current ESMF Time from Clock"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_ClockGet(clock           =CLOCK                     &  !<-- The ESMF Clock
                            ,currTime        =CURRTIME                  &  !<-- The current time (ESMF) on the clock
                            ,rc              =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  The current forecast time.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Lead Write Task Gets Actual Current Time from Clock"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_TimeGet(time=CURRTIME                               &  !<-- The cuurent forecast time (ESMF)
                           ,yy  =IYEAR_FCST                             &  !<-- The current forecast year (integer)
                           ,mm  =IMONTH_FCST                            &  !<-- The current forecast month (integer)
                           ,dd  =IDAY_FCST                              &  !<-- The current forecast day (integer)
                           ,h   =IHOUR_FCST                             &  !<-- The current forecast hour (integer)
                           ,m   =IMINUTE_FCST                           &  !<-- The current forecast minute (integer)
                           ,s   =ISECOND_FCST                           &  !<-- The current forecast second (integer)
                           ,sN  =ISECOND_NUM                            &  !<-- Numerator of current fractional second (integer)
                           ,sD  =ISECOND_DEN                            &  !<-- Denominator of current fractional second (integer)
                           ,rc  =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          SECOND_FCST=ISECOND_FCST+REAL(ISECOND_NUM)/REAL(ISECOND_DEN)     !<-- Current forecast seconds (real)
!
!-----------------------------------------------------------------------
!***  The elapsed forecast time.
!-----------------------------------------------------------------------
!
          wrt_int_state%IO_CURRTIMEDIFF=CURRTIME-wrt_int_state%IO_BASETIME
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Lead Write Task Gets Actual Elapsed Fcst Time"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_TimeIntervalGet(timeinterval=wrt_int_state%IO_CURRTIMEDIFF &
                                   ,h           =NF_HOURS               &  !<-- Hours of elapsed time
                                   ,m           =NF_MINUTES             &  !<-- Minutes of elapsed time
                                   ,s           =NSECONDS               &  !<-- Seconds of elapsed time
                                   ,sN          =NSECONDS_NUM           &  !<-- Numerator of fractional elapsed seconds
                                   ,sD          =NSECONDS_DEN           &  !<-- Denominator of fractional elapsed seconds
                                   ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          NF_SECONDS=NSECONDS+REAL(NSECONDS_NUM)/REAL(NSECONDS_DEN)
!
          wrt_int_state%NFHOURS  =NF_HOURS
          wrt_int_state%NFMINUTES=NF_MINUTES
          wrt_int_state%NFSECONDS=NF_SECONDS
!
        ENDIF
!
      ENDIF hst_time_get
!
!-----------------------------------------------------------------------
!***  DO POST:
!***  Call post processors to compute post variables
!----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
      hst_dopost: IF(wrt_int_state%WRITE_DOPOST.and.                    &
                     wrt_int_state%WRITE_HST_NEMSIO.and.NF_HOURS>0 )THEN               !<-- do post
!-----------------------------------------------------------------------
!

        IF(LOG_PESET)THEN                                               !<-- The write tasks
!
          POST_gridtype='B'
          POST_MAPTYPE=205
          NSOIL=4
!
          CALL POST_RUN_NMM(wrt_int_state,MYPE,MPI_COMM_COMP,               &
                        LEAD_WRITE_TASK,post_gridtype,   &
                        post_maptype,NSOIL,NF_HOURS,NF_MINUTES)
            print *,'af post_run_nmm,NF_HOURS=',NF_HOURS
!
        ENDIF

      ENDIF hst_dopost
!
!-----------------------------------------------------------------------
!***  The lead Write task now opens the history file(s) and writes
!***  the scalar/1D quantities.
!-----------------------------------------------------------------------
!
      IF(MYPE==LEAD_WRITE_TASK)THEN
!
        IF(wrt_int_state%WRITE_HST_BIN)THEN
!
          CALL WRITE_RUNHISTORY_OPEN(WRT_INT_STATE                      &
                                    ,IYEAR_FCST                         &
                                    ,IMONTH_FCST                        &
                                    ,IDAY_FCST                          &
                                    ,IHOUR_FCST                         &
                                    ,IMINUTE_FCST                       &
                                    ,SECOND_FCST                        &
                                    ,NF_HOURS                           &
                                    ,NF_MINUTES                         &
                                    ,NF_SECONDS                         &
                                    ,HST_FIRST                          &
                                    ,LEAD_WRITE_TASK)
        ENDIF
!
        IF(wrt_int_state%WRITE_HST_NEMSIO)THEN
!
          CALL WRITE_NEMSIO_RUNHISTORY_OPEN(WRT_INT_STATE               &
                                           ,NEMSIOFILE                  &
                                           ,IYEAR_FCST                  &
                                           ,IMONTH_FCST                 &
                                           ,IDAY_FCST                   &
                                           ,IHOUR_FCST                  &
                                           ,IMINUTE_FCST                &
                                           ,SECOND_FCST                 &
                                           ,NF_HOURS                    &
                                           ,NF_MINUTES                  &
                                           ,NF_SECONDS                  &
                                           ,DIM1,DIM2,NBDR,GLOBAL       &
                                           ,LEAD_WRITE_TASK,ID_DOMAIN)
!
          FIELDSIZE=(DIM1+2*NBDR)*(DIM2+2*NBDR)
          ALLOCATE(TMP(FIELDSIZE))
        ENDIF
!
      ENDIF
!
!-----------------------------------------------------------------------
!***  We will now assemble the full domain 2-D history data onto
!***  the lead Write task from the subsections on all Write tasks
!***  then the lead task will write each 2D field to the history
!***  file(s). 
!
!***  NOTE:  The lead Write task assembles and writes to history only
!***         one 2-D Field at a time.  
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  First loop through all of the integer Fields.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
      field_loop_int: DO NFIELD=1,KOUNT_I2D                                !<-- Loop through all 2D integer gridded history data
!-----------------------------------------------------------------------
!
        IF(MYPE>LEAD_WRITE_TASK)THEN                                       !<-- All write tasks except the lead one
          JSTA_WRITE=wrt_int_state%LOCAL_JSTART(wrt_int_state%ID_FTASK_RECV_STA(MYPE))  !<-- Starting J of this write task's subsection
          JEND_WRITE=wrt_int_state%LOCAL_JEND  (wrt_int_state%ID_FTASK_RECV_END(MYPE))  !<-- Ending J of this write task's subsection
!
          NN=0
!
          DO J=JSTA_WRITE,JEND_WRITE
          DO I=1,IM
            NN=NN+1
            wrt_int_state%BUFF_INT(NN)=wrt_int_state%WRITE_SUBSET_I(I,J,NFIELD)
          ENDDO
          ENDDO
!
          CALL MPI_RECV(ID_DUMMY                                        &  !<-- Blocking Recv keeps the following sends in line
                       ,1                                               &  !<-- Length of ID_DUMMY
                       ,MPI_INTEGER                                     &  !<-- Datatype
                       ,0                                               &  !<-- The lead write task sent this 
                       ,0                                               &  !<-- An MPI tag
                       ,MPI_COMM_COMP                                   &  !<-- The communicator
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR )
!
          MY_LOCAL_ID=MYPE-LAST_FCST_TASK-1                                !<-- This write task's local ID (between 0 and NWTPG-1)
!
          CALL MPI_SEND(wrt_int_state%BUFF_INT                          &  !<-- Send this string of subsection data 
                       ,NN                                              &  !<-- Number of words sent
                       ,MPI_INTEGER                                     &  !<-- Datatype
                       ,0                                               &  !<-- Send the data to the lead write task with local ID of 0
                       ,MY_LOCAL_ID                                     &  !<-- An MPI tag
                       ,MPI_COMM_COMP                                   &  !<-- MPI communicator
                       ,IERR )
!
!-----------------------------------------------------------------------
!
        ELSEIF(MYPE==LEAD_WRITE_TASK)THEN                                  !<-- The lead write task
!
          JSTA_WRITE=wrt_int_state%LOCAL_JSTART(wrt_int_state%ID_FTASK_RECV_STA(MYPE))  !<-- Starting J of lead write task's subsection
          JEND_WRITE=wrt_int_state%LOCAL_JEND  (wrt_int_state%ID_FTASK_RECV_END(MYPE))  !<-- Ending J of lead write task's subsection
!
          DO J=JSTA_WRITE,JEND_WRITE
          DO I=1,IM
            wrt_int_state%OUTPUT_ARRAY_I2D(I,J)=wrt_int_state%WRITE_SUBSET_I(I,J,NFIELD) !<-- Lead write task fills its part of full domain
          ENDDO
          ENDDO
!
          IF(LAST_WRITE_TASK>LEAD_WRITE_TASK)THEN                          !<-- Recv output subsections if more than 1 write task
            DO N=1,NWTPG-1                                                 !<-- Loop through local IDs of all other write tasks
                                                                           !    that send to the lead task
!
              ID_SEND=N+LEAD_WRITE_TASK
              CALL MPI_SEND(N                                           &  !<-- Send to other write tasks to keep their sends in line
                           ,1                                           &  !<-- Number of words sent
                           ,MPI_INTEGER                                 &  !<-- Datatype
                           ,N                                           &  !<-- Send to each of the other write tasks
                           ,0                                           &  !<-- An MPI tag
                           ,MPI_COMM_COMP                               &  !<-- MPI communicator
                           ,IERR )
!
              ID_RECV=ID_SEND
              CALL MPI_RECV(wrt_int_state%BUFF_INT                      &  !<-- Recv string of subsection data from other write tasks
                           ,IM*JM                                       &  !<-- Maximum number of words sent
                           ,MPI_INTEGER                                 &  !<-- Datatype
                           ,N                                           &  !<-- Recv from this write task
                           ,N                                           &  !<-- An MPI tag
                           ,MPI_COMM_COMP                               &  !<-- MPI communicator
                           ,JSTAT                                       &  !<-- MPI status object
                           ,IERR )
!
              NN=0
              JSTA_WRITE=wrt_int_state%LOCAL_JSTART(wrt_int_state%ID_FTASK_RECV_STA(N+LEAD_WRITE_TASK)) !<-- Starting J of sending write task
              JEND_WRITE=wrt_int_state%LOCAL_JEND  (wrt_int_state%ID_FTASK_RECV_END(N+LEAD_WRITE_TASK)) !<-- Ending J of sending write task
!
              DO J=JSTA_WRITE,JEND_WRITE
              DO I=1,IM
                NN=NN+1
                wrt_int_state%OUTPUT_ARRAY_I2D(I,J)=wrt_int_state%BUFF_INT(NN)  !<-- Insert other write tasks' subsections into full domain
              ENDDO
              ENDDO
!
            ENDDO
          ENDIF
!
          NPOSN_1=(NFIELD-1)*ESMF_MAXSTR+1
          NPOSN_2=NFIELD*ESMF_MAXSTR
          NAME=wrt_int_state%NAMES_I2D_STRING(NPOSN_1:NPOSN_2)                        !<-- The name of this 2D integer history quantity
! 
          IF(wrt_int_state%WRITE_HST_BIN)THEN
!
            WRITE(wrt_int_state%IO_HST_UNIT,iostat=RC)wrt_int_state%OUTPUT_ARRAY_I2D  !<-- Lead write task writes out the 2D real data
!
            IF(HST_FIRST .AND. (wrt_int_state%PRINT_OUTPUT .OR. wrt_int_state%PRINT_ALL) )THEN
              WRITE(0,*)'Wrote ',TRIM(NAME),' to history file unit ',wrt_int_state%IO_HST_UNIT
            ENDIF
          ENDIF
!
!-----------------------------------------------------------------------
!***  For the NEMSIO history file.
!-----------------------------------------------------------------------
!
          IF(wrt_int_state%WRITE_HST_NEMSIO)THEN
!
            IF(FIELDSIZE/=IM*JM)THEN
              WRITE(0,*)'WRONG: input data dimension ',IM*JM,           &
               ' does not match data size in NEMSIO file ',FIELDSIZE
            ENDIF
!
            TMP=RESHAPE(wrt_int_state%OUTPUT_ARRAY_I2D(1:IM,1:JM)       &
                       ,(/FIELDSIZE/))
!
            CALL NEMSIO_WRITEREC(NEMSIOFILE,NFIELD,TMP,IRET=IERR)           !<-- Lead write task writes out the 2D int data!
            IF(IERR/=0)THEN
              WRITE(0,*)' Failed to write output to file! Aborting!'
              CALL ESMF_Finalize(rc=RC,endflag=ESMF_END_ABORT)
            ENDIF
!
         ENDIF
!-----------------------------------------------------------------------
!
        ENDIF
!
!-----------------------------------------------------------------------
!
      ENDDO field_loop_int
!
!-----------------------------------------------------------------------
!***  Now loop through all the real Fields.
!-----------------------------------------------------------------------
!
      WRITE(MODEL_LEVEL,'(I3.3)')wrt_int_state%LM(1)
!
!-----------------------------------------------------------------------
!
      field_loop_real: DO NFIELD=1,KOUNT_R2D                               !<-- Loop through all 2D real gridded history data
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
        IF(MYPE>LEAD_WRITE_TASK)THEN                                       !<-- All write tasks except the lead one
!-----------------------------------------------------------------------
!
          JSTA_WRITE=wrt_int_state%LOCAL_JSTART(wrt_int_state%ID_FTASK_RECV_STA(MYPE))  !<-- Starting J of this write task's subsection
          JEND_WRITE=wrt_int_state%LOCAL_JEND  (wrt_int_state%ID_FTASK_RECV_END(MYPE))  !<-- Ending J of this write task's subsection
!
          NN=0
!
          DO J=JSTA_WRITE,JEND_WRITE
          DO I=1,IM
            NN=NN+1
            wrt_int_state%BUFF_REAL(NN)=wrt_int_state%WRITE_SUBSET_R(I,J,NFIELD)
          ENDDO
          ENDDO
!
          CALL MPI_RECV(ID_DUMMY                                        &  !<-- Blocking Recv keeps the following sends in line
                       ,1                                               &  !<-- Length of ID_DUMMY
                       ,MPI_INTEGER                                     &  !<-- Datatype
                       ,0                                               &  !<-- The lead write task sent this 
                       ,0                                               &  !<-- An MPI tag
                       ,MPI_COMM_COMP                                   &  !<-- The communicator
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR )
!
          MY_LOCAL_ID=MYPE-LAST_FCST_TASK-1                                !<-- This write task's local ID (between 0 and NWTPG-1)
!
          CALL MPI_SEND(wrt_int_state%BUFF_REAL                         &  !<-- Send this string of subsection data 
                       ,NN                                              &  !<-- Number of words sent
                       ,MPI_REAL                                        &  !<-- Datatype
                       ,0                                               &  !<-- Send the data to the lead write task with local ID of 0
                       ,MY_LOCAL_ID                                     &  !<-- An MPI tag
                       ,MPI_COMM_COMP                                   &  !<-- MPI communicator
                       ,IERR )
!
!-----------------------------------------------------------------------
        ELSEIF(MYPE==LEAD_WRITE_TASK)THEN                                  !<-- The lead write task
!-----------------------------------------------------------------------
!
          JSTA_WRITE=wrt_int_state%LOCAL_JSTART(wrt_int_state%ID_FTASK_RECV_STA(MYPE))  !<-- Starting J of lead write task's subsection
          JEND_WRITE=wrt_int_state%LOCAL_JEND  (wrt_int_state%ID_FTASK_RECV_END(MYPE))  !<-- Ending J of lead write task's subsection
!
          DO J=JSTA_WRITE,JEND_WRITE
          DO I=1,IM
            wrt_int_state%OUTPUT_ARRAY_R2D(I,J)=wrt_int_state%WRITE_SUBSET_R(I,J,NFIELD) !<-- Lead write task fills its part of full domain
          ENDDO
          ENDDO
!
          IF(LAST_WRITE_TASK>LEAD_WRITE_TASK)THEN                          !<-- Recv output subsections if more than 1 write task
            DO N=1,NWTPG-1                                                 !<-- Loop through local IDs of all other write tasks
                                                                           !    that send to the lead task
!
              ID_SEND=N+LEAD_WRITE_TASK
              CALL MPI_SEND(N                                           &  !<-- Send to other write tasks to keep their sends in line
                           ,1                                           &  !<-- Number of words sent
                           ,MPI_INTEGER                                 &  !<-- Datatype
                           ,N                                           &  !<-- Send to each of the other write tasks
                           ,0                                           &  !<-- An MPI tag
                           ,MPI_COMM_COMP                               &  !<-- MPI communicator
                           ,IERR )
!
              ID_RECV=ID_SEND
              CALL MPI_RECV(wrt_int_state%BUFF_REAL                     &  !<-- Recv string of subsection data from other write tasks
                           ,IM*JM                                       &  !<-- Maximum number of words sent
                           ,MPI_REAL                                    &  !<-- Datatype
                           ,N                                           &  !<-- Recv from this write task
                           ,N                                           &  !<-- An MPI tag
                           ,MPI_COMM_COMP                               &  !<-- MPI communicator
                           ,JSTAT                                       &  !<-- MPI status object
                           ,IERR )
!
              NN=0
              JSTA_WRITE=wrt_int_state%LOCAL_JSTART(wrt_int_state%ID_FTASK_RECV_STA(N+LEAD_WRITE_TASK)) !<-- Starting J of sending write task
              JEND_WRITE=wrt_int_state%LOCAL_JEND  (wrt_int_state%ID_FTASK_RECV_END(N+LEAD_WRITE_TASK)) !<-- Ending J of sending write task
!
              DO J=JSTA_WRITE,JEND_WRITE
              DO I=1,IM
                NN=NN+1
                wrt_int_state%OUTPUT_ARRAY_R2D(I,J)=wrt_int_state%BUFF_REAL(NN)   !<-- Insert other write tasks' subsections into full domain
              ENDDO
              ENDDO
!
            ENDDO
          ENDIF
!
          NPOSN_1=(NFIELD-1)*ESMF_MAXSTR+1
          NPOSN_2=NFIELD*ESMF_MAXSTR
          NAME=wrt_int_state%NAMES_R2D_STRING(NPOSN_1:NPOSN_2)                       !<-- The name of this 2D real history quantity
!
!-----------------------------------------------------------------------
!***  Begin computation of the 10-m wind factor for GSI.
!-----------------------------------------------------------------------
!
          IF(TRIM(NAME)=='U10') THEN
            IF(.NOT.ALLOCATED(FACT10)) THEN
              ALLOCATE(FACT10(1:IM,1:JM))
!
              DO J=1,JM
              DO I=1,IM
                FACT10(I,J)=0.
              ENDDO
              ENDDO
            ENDIF
!
            DO J=1,JM
            DO I=1,IM
              FACT10(I,J)=FACT10(I,J)+                                  &
                          wrt_int_state%OUTPUT_ARRAY_R2D(I,J)*          &
                          wrt_int_state%OUTPUT_ARRAY_R2D(I,J)
            ENDDO
            ENDDO
          ENDIF
!
          IF(TRIM(NAME)=='V10') THEN
            IF(.NOT.ALLOCATED(FACT10)) THEN
              ALLOCATE(FACT10(1:IM,1:JM))
              FACT10=0.
            ENDIF
!
            DO J=1,JM
            DO I=1,IM
              FACT10(I,J)=FACT10(I,J)+                                 &
                          wrt_int_state%OUTPUT_ARRAY_R2D(I,J)*         &
                          wrt_int_state%OUTPUT_ARRAY_R2D(I,J)
            ENDDO
            ENDDO
          ENDIF
!
          IF(TRIM(NAME)=='U_'//MODEL_LEVEL//'_2D') THEN
            ALLOCATE(FACT10TMPU(1:IM,1:JM))
            CALL V_TO_H_BGRID(wrt_int_state%OUTPUT_ARRAY_R2D(1:IM,1:JM) &
                             ,IM,JM,GLOBAL,FACT10TMPU)
          ENDIF
!
          IF(TRIM(NAME)=='V_'//MODEL_LEVEL//'_2D') THEN
            ALLOCATE(FACT10TMPV(1:IM,1:JM))
            CALL V_TO_H_BGRID(wrt_int_state%OUTPUT_ARRAY_R2D(1:IM,1:JM) &
                             ,IM,JM,GLOBAL,FACT10TMPV)
          ENDIF
!
          IF(wrt_int_state%WRITE_HST_BIN)THEN
!
            WRITE(wrt_int_state%IO_HST_UNIT,iostat=RC)wrt_int_state%OUTPUT_ARRAY_R2D   !<-- Lead write task writes out the 2D real data
!
            IF(HST_FIRST .AND. (wrt_int_state%PRINT_OUTPUT .OR. wrt_int_state%PRINT_ALL) )THEN
              WRITE(0,*)'Wrote ',TRIM(NAME)                                &
                       ,' to history file unit ',wrt_int_state%IO_HST_UNIT &
                       ,MAXVAL(wrt_int_state%OUTPUT_ARRAY_R2D)             &
                       ,MINVAL(wrt_int_state%OUTPUT_ARRAY_R2D)
            ENDIF
          ENDIF
!
!-----------------------------------------------------------------------
!***  The same for the NEMSIO history file.
!-----------------------------------------------------------------------
!
          IF(wrt_int_state%WRITE_HST_NEMSIO)THEN
!
            IF(FIELDSIZE/=IM*JM)THEN
              WRITE(0,*)'WRONG: data dimension ',IM*JM,                 &
               ' does not match data size in NEMSIO file,',FIELDSIZE
            ENDIF
!
            IF(TRIM(NAME)=='FIS')wrt_int_state%OUTPUT_ARRAY_R2D(1:IM,1:JM)=    &
                  wrt_int_state%OUTPUT_ARRAY_R2D(1:IM,1:JM)/G
!
            IF(TRIM(NAME)=='GLAT')THEN
               ALLOCATE(GLAT1D(FIELDSIZE))
               GLAT1D(1:FIELDSIZE)=RESHAPE(wrt_int_state%OUTPUT_ARRAY_R2D(1:IM,1:JM),(/FIELDSIZE/))
            ENDIF
!
            IF(TRIM(NAME)=='GLON')THEN
               ALLOCATE(GLON1D(FIELDSIZE))
               GLON1D(1:FIELDSIZE)=RESHAPE(wrt_int_state%OUTPUT_ARRAY_R2D(1:IM,1:JM),(/FIELDSIZE/))
            ENDIF
!
            N=NFIELD+wrt_int_state%KOUNT_I2D(1)
            TMP=RESHAPE(wrt_int_state%OUTPUT_ARRAY_R2D(1:IM,1:JM),(/FIELDSIZE/))
!
            CALL NEMSIO_WRITEREC(NEMSIOFILE,N,TMP,IRET=IERR)
            IF(IERR/=0)THEN
              WRITE(0,*)' Failed to write output to file! Aborting!'
              CALL ESMF_Finalize(rc=RC,endflag=ESMF_END_ABORT)
            ENDIF
!
!           IF(HST_FIRST)THEN
!             WRITE(0,*)'Wrote ',TRIM(NAME),' to nemsio history file iret=',ierr
!           ENDIF
!
          ENDIF
!
        ENDIF
!
!-----------------------------------------------------------------------
!
      ENDDO field_loop_real
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Complete computation of 10-m wind factor and write it out.
!-----------------------------------------------------------------------
!
      IF( MYPE==LEAD_WRITE_TASK )THEN
        IF(ALLOCATED(FACT10).AND.ALLOCATED(FACT10TMPU).AND.ALLOCATED(FACT10TMPV)) THEN
          DO J=1,JM
          DO I=1,IM
            FACT10TMPV(I,J)=SQRT(FACT10TMPU(I,J)*FACT10TMPU(I,J)+       &
                                 FACT10TMPV(I,J)*FACT10TMPV(I,J))
          ENDDO
          ENDDO
!
!         write(0,*)'wind mgn=',maxval(FACT10TMPV(1:IM,1:JM)),minval(FACT10TMPV(1:IM,1:JM))
!         write(0,*)'wind10 mgn=',maxval(sqrt(FACT10(1:IM,1:JM))),minval(sqrt(FACT10(1:IM,1:JM)))
!
          DO J=1,JM
          DO I=1,IM
            IF(FACT10TMPV(I,J)/=0) THEN
              FACT10(I,J)=SQRT(FACT10(I,J))/FACT10TMPV(I,J)
            ELSE
              FACT10(I,J)=1.
            ENDIF
          ENDDO
          ENDDO
!
          DEALLOCATE(FACT10TMPU)
          DEALLOCATE(FACT10TMPV)
!
          IF(wrt_int_state%WRITE_HST_BIN)THEN
            WRITE(wrt_int_state%IO_HST_UNIT,iostat=RC)FACT10                    !<-- Lead write task writes out the 2D real data
            IF(HST_FIRST .AND. (wrt_int_state%PRINT_OUTPUT .OR. wrt_int_state%PRINT_ALL) )THEN
              WRITE(0,*)'Wrote FACT10 to history file unit ',wrt_int_state%IO_HST_UNIT &
                        ,maxval(fact10),minval(fact10)
            ENDIF
          ENDIF
!
          IF(wrt_int_state%WRITE_HST_NEMSIO)THEN
            N=N+1
            TMP=RESHAPE(FACT10(1:IM,1:JM),(/FIELDSIZE/))
            CALL NEMSIO_WRITEREC(NEMSIOFILE,N,TMP,IRET=IERR)
            IF(IERR/=0)THEN
              WRITE(0,*)' Failed to write output to file! Aborting!'
              CALL ESMF_Finalize(rc=RC,endflag=ESMF_END_ABORT)
            ENDIF
!           write(0,*)'after nemsio_writerec,n=',n,'fact10=',maxval(tmp),minval(tmp),'iret=',ierr
          ENDIF
!
          DEALLOCATE(FACT10)
!
        ENDIF
!
      ENDIF
!
      HST_FIRST=.FALSE.
!
!-----------------------------------------------------------------------
!***  Close the history file if needed.
!-----------------------------------------------------------------------
!
      IF(wrt_int_state%WRITE_HST_BIN.and.MYPE==LEAD_WRITE_TASK)THEN
        CLOSE(wrt_int_state%IO_HST_UNIT)
!       write(0,*)' Closed history file with unit=',wrt_int_state%IO_HST_UNIT
      ENDIF
!
!-----------------------------------------------------------------------
!***  Close the NEMSIO history file.
!-----------------------------------------------------------------------
!
      IF(wrt_int_state%WRITE_HST_NEMSIO.AND.wrt_int_state%MYPE==LEAD_WRITE_TASK)THEN
!
        IF(ASSOCIATED(GLAT1D).AND.ASSOCIATED(GLON1D)) THEN
          CALL NEMSIO_SETFILEHEAD(NEMSIOFILE,IERR,GLAT1D,GLON1D)
          DEALLOCATE(GLAT1D,GLON1D)
        ENDIF
!
        CALL NEMSIO_GETFILEHEAD(NEMSIOFILE,IERR,gfname=GFNAME)
!
        DEALLOCATE(TMP)
!
        CALL NEMSIO_CLOSE(NEMSIOFILE)
!       WRITE(0,*)' Closed nemsio_history file, ', gfname
!
        CALL NEMSIO_FINALIZE()
!
! ffsync
        IF(WRT_INT_STATE%WRITE_FSYNCFLAG) THEN
          DO N=51,99
            INQUIRE(N,opened=OPENED)
            IF(.NOT.OPENED)THEN
              IO_HST_UNIT=N
              EXIT
            ENDIF
          ENDDO
!
          OPEN(unit=IO_HST_UNIT, file=trim(GFNAME) )
!
          RC=FFSYNC(IO_HST_UNIT)

! Handle possible error
          IF (RC.NE.0) THEN
            print *,"Error returned from ffsync, file=",                &
     &            trim(GFNAME),"rc=",RC
          ENDIF
!
          CLOSE(IO_HST_UNIT)
!
        ENDIF
!
!ffsync end
!
      ENDIF
!
!-----------------------------------------------------------------------
!***  Get this domain's configure object.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Get Config Object for Write Setup"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompGet(gridcomp=WRITE_COMP                         &  !<-- The Write component
                           ,config  =CF                                 &  !<-- The configure object on this domain
                           ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)   
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-------------------
!***  Fcstdone file
!-------------------
!
      IF(wrt_int_state%WRITE_DONEFILEFLAG                               &
          .AND.                                                         &
         wrt_int_state%MYPE==LEAD_WRITE_TASK                            &
          .AND.                                                         &
        (wrt_int_state%WRITE_HST_BIN                                    &
          .OR.                                                          &
         wrt_int_state%WRITE_HST_NEMSIO))THEN
!
        CALL ESMF_ConfigGetAttribute(config=CF                          &  !<-- The configure file object
                                    ,value =ID_DOMAIN                   &  !<-- Put extracted quantity here
                                    ,label ='my_domain_id:'             &  !<-- The quantity's label in the configure file
                                    ,rc    =RC)
!
        INT_SEC=INT(NF_SECONDS)
        FRAC_SEC=NINT((NF_SECONDS-INT_SEC)*100.)
!
        WRITE(FILENAME,'(A,I2.2,A,I4.4,A,I2.2,A,I2.2,A,I2.2,A)' )       &  
                                         'fcstdone.',ID_DOMAIN,'.'      &  !<-- Insert the domain ID for nests
                                         ,NF_HOURS,'h_'                 & 
                                         ,NF_MINUTES,'m_'               &
                                         ,INT_SEC,'.',FRAC_SEC,'s'
!
        DO N=51,99
          INQUIRE(N,opened=OPENED)
          IF(.NOT.OPENED)THEN
            IO_HST_UNIT=N
            EXIT
          ENDIF
        ENDDO
!
        OPEN(unit  =IO_HST_UNIT                                         &
            ,file  =trim(FILENAME)                                      &
            ,form  ='formatted'                                         &
            ,status='REPLACE')
!
        WRITE(IO_HST_UNIT,'(A4)')'DONE'
        CLOSE(IO_HST_UNIT)
!
      ENDIF
!
!-----------------------------------------------------------------------
!
      ENDIF history_time
!
!-----------------------------------------------------------------------
!
      restart_time: IF(TIME_FOR_RESTART) THEN
!
!-----------------------------------------------------------------------
!***  The lead Write task receives the latest Attributes from the
!***  lead fcst task for restart output.
!-----------------------------------------------------------------------
!
      IF(MYPE==LEAD_WRITE_TASK)THEN
!
        WRITE_GROUP=NCURRENT_GROUP(ID_DOMAIN)
!
!--------------
!***  Integers
!--------------
!
        IF(wrt_int_state%RST_LENGTH_SUM_I1D(1)>0)THEN
!
          CALL MPI_RECV(wrt_int_state%RST_ALL_DATA_I1D                &  !<-- Recv string of integer restart Attributes
                       ,wrt_int_state%RST_LENGTH_SUM_I1D(1)           &  !<-- Words received
                       ,MPI_INTEGER                                   &  !<-- Data is integer
                       ,0                                             &  !<-- Sending task (lead fcst task)
                       ,WRITE_GROUP                                   &  !<-- MPI tag
                       ,INTERCOMM_WRITE_GROUP                         &  !<-- MPI domain commumicator
                       ,JSTAT                                         &  !<-- MPI status object
                       ,IERR )
!     write(0,*)' Write Run lead write task recvd RST_ALL_DATA_I1D with length=',wrt_int_state%RST_LENGTH_SUM_I1D(1)
!     write(0,*)' NMTS=wrt_int_state%RST_ALL_DATA_I1D(10)=',wrt_int_state%RST_ALL_DATA_I1D(10)
!
        ENDIF
!
!-----------
!***    Reals
!-----------
!
        IF(wrt_int_state%RST_LENGTH_SUM_R1D(1)>0)THEN
!
          CALL MPI_RECV(wrt_int_state%RST_ALL_DATA_R1D                &  !<-- Recv string of real restart Attributes
                       ,wrt_int_state%RST_LENGTH_SUM_R1D(1)           &  !<-- Words received
                       ,MPI_REAL                                      &  !<-- Data is real
                       ,0                                             &  !<-- Sending task (lead fcst task)
                       ,WRITE_GROUP                                   &  !<-- MPI tag
                       ,INTERCOMM_WRITE_GROUP                         &  !<-- MPI domain commumicator
                       ,JSTAT                                         &  !<-- MPI status object
                       ,IERR )
!
        ENDIF
!
!--------------
!***    Logicals
!--------------
!
        IF(wrt_int_state%RST_LENGTH_SUM_LOG(1)>0)THEN
!
          CALL MPI_RECV(wrt_int_state%RST_ALL_DATA_LOG                &  !<-- Recv string of logical restart Attributes
                       ,wrt_int_state%RST_LENGTH_SUM_LOG(1)           &  !<-- Words received
                       ,MPI_LOGICAL                                   &  !<-- Data is logical
                       ,0                                             &  !<-- Sending task (lead fcst task)
                       ,WRITE_GROUP                                   &  !<-- MPI tag
                       ,INTERCOMM_WRITE_GROUP                         &  !<-- MPI domain commumicator
                       ,JSTAT                                         &  !<-- MPI status object
                       ,IERR )
!
        ENDIF
!
      ENDIF
!
!-----------------------------------------------------------------------
!***  Each Write task in the active Write group receives the
!***  strings of 2D restart data from the appropriate fcst tasks.
!-----------------------------------------------------------------------
!
      ID_START=wrt_int_state%ID_FTASK_RECV_STA(MYPE)                       !<-- First fcst task that sends to this write task
      ID_END  =wrt_int_state%ID_FTASK_RECV_END(MYPE)                       !<-- Last fcst task that sends to this write task
      NFCST_TASKS=ID_END-ID_START+1                                        !<-- Number of fcst tasks sending to this write task
!
!-----------------------------------------------------------------------
      rst_from_fcst_tasks: DO N=1,NFCST_TASKS                              !<-- Loop through fcst tasks sending to this write task
!-----------------------------------------------------------------------
!
        ID_RECV=ID_START+N-1
!
!-----------------------------------------------------------------------
!***  Receive 2-D integer data if there is any.
!-----------------------------------------------------------------------
!
        IF(RST_KOUNT_I2D>0)THEN
          CALL MPI_RECV(wrt_int_state%RST_ALL_DATA_I2D                  &  !<-- Fcst tasks' string of 2D integer restart data
                       ,wrt_int_state%NUM_WORDS_RECV_I2D_RST(ID_RECV)   &  !<-- # of words in the data string
                       ,MPI_INTEGER                                     &  !<-- The datatype
                       ,ID_RECV                                         &  !<-- Recv from this fcst task
                       ,wrt_int_state%NFHOURS                           &  !<-- An MPI tag
                       ,INTERCOMM_WRITE_GROUP                           &  !<-- The MPI intercommunicator between fcst and quilt tasks
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR )
!
          IF(IERR/=0)WRITE(0,*)' Recv by write task from fcst task has failed.  IERR=',IERR
        ENDIF
!
!-----------------------------------------------------------------------
!***  Receive 2-D real data if there is any.
!-----------------------------------------------------------------------
!
        IF(RST_KOUNT_R2D>0)THEN
          CALL MPI_RECV(wrt_int_state%RST_ALL_DATA_R2D                  &  !<-- Fcst tasks' string of 2D real restart data
                       ,wrt_int_state%NUM_WORDS_RECV_R2D_RST(ID_RECV)   &  !<-- # of words in the data string
                       ,MPI_REAL                                        &  !<-- The datatype
                       ,ID_RECV                                         &  !<-- Recv from this fcst task
                       ,wrt_int_state%NFHOURS                           &  !<-- An MPI tag
                       ,INTERCOMM_WRITE_GROUP                           &  !<-- The MPI intercommunicator between fcst and quilt tasks
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR )
!
          IF(IERR/=0)WRITE(0,*)' Recv by write task from fcst task has failed.  IERR=',IERR
        ENDIF
!
!-----------------------------------------------------------------------
!***  Each Write task needs to insert the pieces of the various
!***  2D restart arrays received from the individual fcst tasks
!***  into arrays that span the Write tasks' own subsection of
!***  the full 2D domain.  That subsection always spans the 
!***  entire East-West dimension of the full domain (since full
!***  rows of Fcst tasks always send to Write tasks, never 
!***  partial rows) and as much of the North-South dimension of
!***  the full domain as covered by those Fcst tasks sending to
!***  a given Write task.
!-----------------------------------------------------------------------
!
        ITS=wrt_int_state%LOCAL_ISTART(ID_RECV)                            !<-- Local domain integration limits of sending fcst task
        ITE=wrt_int_state%LOCAL_IEND(ID_RECV)                              !<--
        JTS=wrt_int_state%LOCAL_JSTART(ID_RECV)                            !<--
        JTE=wrt_int_state%LOCAL_JEND(ID_RECV)                              !<--
!
        IHALO=wrt_int_state%IHALO                                          !<-- Subdomain halo depth in I
        JHALO=wrt_int_state%JHALO                                          !<-- Subdomain halo depth in J
!
        IMS=ITS-IHALO
        IME=ITE+IHALO
        JMS=JTS-JHALO
        JME=JTE+JHALO
!
        NN=0
!
        DO NF=1,RST_KOUNT_I2D                                              !<-- Loop through all the 2D integer fields
!
          DO J=JMS,JME
          DO I=IMS,IME
            NN=NN+1
            IF(I<ITS.OR.I>ITE.OR.J<JTS.OR.J>JTE)CYCLE                      !<-- Exclude halo points
            wrt_int_state%RST_WRITE_SUBSET_I(I,J,NF)=wrt_int_state%RST_ALL_DATA_I2D(NN) !<-- Put data into write task's domain subsection
          ENDDO
          ENDDO
        ENDDO
!
        NN=0
!
        DO NF=1,RST_KOUNT_R2D                                              !<-- Loop through all the 2D real fields
!
          DO J=JMS,JME
          DO I=IMS,IME
            NN=NN+1
            IF(I<ITS.OR.I>ITE.OR.J<JTS.OR.J>JTE)CYCLE                      !<-- Exclude halo points
            wrt_int_state%RST_WRITE_SUBSET_R(I,J,NF)=wrt_int_state%RST_ALL_DATA_R2D(NN) !<-- Put data into write task's domain subsection
          ENDDO
          ENDDO
        ENDDO
!
!-----------------------------------------------------------------------
!
      ENDDO rst_from_fcst_tasks
!
!-----------------------------------------------------------------------
!***  The lead Write task must receive the full-domain boundary 
!***  data from the lead Forecast task.  The data is required in
!***  the restart files when nests are present.
!-----------------------------------------------------------------------
!
      insert_bc_data: IF(MYPE==LEAD_WRITE_TASK)THEN                        !<-- The lead write task
!
        CALL MPI_RECV(wrt_int_state%RST_ALL_BC_DATA                     &  !<-- 1-D string of BC data for restart
                     ,wrt_int_state%NUM_WORDS_SEND_BC(1)                &  !<-- # of words in the data string
                     ,MPI_REAL                                          &  !<-- The datatype
                     ,0                                                 &  !<-- Recv from this fcst task
                     ,wrt_int_state%NFHOURS                             &  !<-- An MPI tag
                     ,INTERCOMM_WRITE_GROUP                             &  !<-- The MPI intercommunicator between fcst and quilt tasks
                     ,JSTAT                                             &  !<-- MPI status object
                     ,IERR )
!
        IF(IERR/=0)WRITE(0,*)' Recv by write task of fcst task BC data has failed.  IERR=',IERR
!
!-----------------------------------------------------------------------
!***  Add the BC data, word count, and name into the appropriate arrays.
!-----------------------------------------------------------------------
!
        once: IF(NSUM_WORDS_NEW==0)THEN                                    !<-- Increase size of RST_ALL_DATA_R1D only once
          NSUM_WORDS=0
          DO NX=1,wrt_int_state%RST_KOUNT_R1D(1)                           !<-- Loop through all 1-D Real items in string
            NSUM_WORDS=NSUM_WORDS+wrt_int_state%RST_LENGTH_DATA_R1D(NX)    !<-- # of words so far in 1-D real data string
          ENDDO
!
          ALLOCATE(HOLD_RST_DATA_R1D(1:NSUM_WORDS))
          DO NX=1,NSUM_WORDS
            HOLD_RST_DATA_R1D(NX)=wrt_int_state%RST_ALL_DATA_R1D(NX)       !<-- Save those words temporarily
          ENDDO
!
          NSUM_WORDS_NEW=NSUM_WORDS+wrt_int_state%NUM_WORDS_SEND_BC(1)     !<-- # of 1-D real words including the BC data
          DEALLOCATE(wrt_int_state%RST_ALL_DATA_R1D)
          ALLOCATE(wrt_int_state%RST_ALL_DATA_R1D(1:NSUM_WORDS_NEW)     &  !<-- Reallocate the 1-D real storage to new length
                  ,stat=ISTAT)                                             !<--
!
          DO NX=1,NSUM_WORDS
            wrt_int_state%RST_ALL_DATA_R1D(NX)=HOLD_RST_DATA_R1D(NX)       !<-- Transfer the original data back
          ENDDO
!
          wrt_int_state%RST_KOUNT_R1D(1)=wrt_int_state%RST_KOUNT_R1D(1)+1  !<-- Increment the 1-D real item count by 1
!
          wrt_int_state%RST_LENGTH_DATA_R1D(wrt_int_state%RST_KOUNT_R1D)=   &  !<-- Insert the BC data's word count
                                         wrt_int_state%NUM_WORDS_SEND_BC(1)
!
          NPOSN_2=wrt_int_state%RST_KOUNT_R1D(1)*ESMF_MAXSTR
          NPOSN_1=NPOSN_2-ESMF_MAXSTR+1
          wrt_int_state%RST_NAMES_R1D_STRING(NPOSN_1:NPOSN_2)='ALL_BC_DATA'    !<-- Insert the data's name
!
        ENDIF once
!
        KOUNT=0
        DO NX=NSUM_WORDS+1,NSUM_WORDS_NEW
          KOUNT=KOUNT+1
          wrt_int_state%RST_ALL_DATA_R1D(NX)=wrt_int_state%RST_ALL_BC_DATA(KOUNT)  !<-- Insert the new BC data
        ENDDO
!
      ENDIF insert_bc_data
!
!-----------------------------------------------------------------------
!***  At this point, all Write tasks have received all of the restart
!***  data from their associated Fcst tasks and assembled it onto 
!***  their own subsections of the full 2D domain.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  The lead Write task obtains the current and elapsed forecast 
!***  times for its writing of the restart output.
!-----------------------------------------------------------------------
!
      rst_write_begin: IF(MYPE==LEAD_WRITE_TASK)THEN                       !<-- The lead write task
!
!-----------------------------------------------------------------------
!
        IF(wrt_int_state%WRITE_RST_BIN.OR.                              &
           wrt_int_state%WRITE_RST_NEMSIO)THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Lead Write Task Gets Current ESMF Time from Clock"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_ClockGet(clock       =CLOCK                         &  !<-- The ESMF Clock
                            ,currTime    =CURRTIME                      &  !<-- The current time (ESMF) on the clock
                            ,runTimeStepCount=RUN_TIMESTEP_COUNT        &  !<-- # of times the clock has advanced
                            ,rc          =RC)
!
          NTIMESTEP=NINT(RUN_TIMESTEP_COUNT)-1
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  The current forecast time.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Lead Write Task Gets Actual Current Time from Clock"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_TimeGet(time=CURRTIME                               &  !<-- The cuurent forecast time (ESMF)
                           ,yy  =IYEAR_FCST                             &  !<-- The current forecast year (integer)
                           ,mm  =IMONTH_FCST                            &  !<-- The current forecast month (integer)
                           ,dd  =IDAY_FCST                              &  !<-- The current forecast day (integer)
                           ,h   =IHOUR_FCST                             &  !<-- The current forecast hour (integer)
                           ,m   =IMINUTE_FCST                           &  !<-- The current forecast minute (integer)
                           ,s   =ISECOND_FCST                           &  !<-- The current forecast second (integer)
                           ,sN  =ISECOND_NUM                            &  !<-- Numerator of current fractional second (integer)
                           ,sD  =ISECOND_DEN                            &  !<-- Denominator of current fractional second (integer)
                           ,rc  =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          SECOND_FCST=ISECOND_FCST+REAL(ISECOND_NUM)/REAL(ISECOND_DEN)     !<-- Current forecast seconds (real)
!
!-----------------------------------------------------------------------
!***  Elapsed forecast time.
!-----------------------------------------------------------------------
!
          wrt_int_state%IO_CURRTIMEDIFF=CURRTIME-wrt_int_state%IO_BASETIME
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Lead Write Task Gets Actual Elapsed Fcst Time"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_TimeIntervalGet(timeinterval=wrt_int_state%IO_CURRTIMEDIFF &
                                   ,h           =NF_HOURS               &  !<-- Hours of elapsed time
                                   ,m           =NF_MINUTES             &  !<-- Minutes of elapsed time
                                   ,s           =NSECONDS               &  !<-- Seconds of elapsed time
                                   ,sN          =NSECONDS_NUM           &  !<-- Numerator of fractional elapsed seconds
                                   ,sD          =NSECONDS_DEN           &  !<-- denominator of fractional elapsed seconds
                                   ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          NF_SECONDS=NSECONDS+REAL(NSECONDS_NUM)/REAL(NSECONDS_DEN)
!
          wrt_int_state%NFHOURS  =NF_HOURS
          wrt_int_state%NFMINUTES=NF_MINUTES
          wrt_int_state%NFSECONDS=NF_SECONDS
!
        ENDIF
!
!-----------------------------------------------------------------------
!
      ENDIF rst_write_begin
!
!-----------------------------------------------------------------------
!***  The lead Write task opens the restart file(s) and writes out
!***  scalar/1D quantities.
!-----------------------------------------------------------------------
!
      IF(MYPE==LEAD_WRITE_TASK)THEN
!
        IF(wrt_int_state%WRITE_RST_BIN)THEN
!
          CALL WRITE_RUNRESTART_OPEN(WRT_INT_STATE                      &
                                    ,IYEAR_FCST                         &
                                    ,IMONTH_FCST                        &
                                    ,IDAY_FCST                          &
                                    ,IHOUR_FCST                         &
                                    ,IMINUTE_FCST                       &
                                    ,SECOND_FCST                        &
                                    ,NTIMESTEP                          &
                                    ,NF_HOURS                           &
                                    ,NF_MINUTES                         &
                                    ,NF_SECONDS                         &
                                    ,RST_FIRST                          &
                                    ,LEAD_WRITE_TASK)
        ENDIF
!
        IF(wrt_int_state%WRITE_RST_NEMSIO)THEN
!
          CALL WRITE_NEMSIO_RUNRESTART_OPEN(WRT_INT_STATE               &
                                           ,NEMSIOFILE                  &
                                           ,IYEAR_FCST                  &
                                           ,IMONTH_FCST                 &
                                           ,IDAY_FCST                   &
                                           ,IHOUR_FCST                  &
                                           ,IMINUTE_FCST                &
                                           ,SECOND_FCST                 &
                                           ,NTIMESTEP                   &
                                           ,NF_HOURS                    &
                                           ,NF_MINUTES                  &
                                           ,NF_SECONDS                  &
                                           ,DIM1,DIM2,NBDR,GLOBAL       &
                                           ,ID_DOMAIN                   &
                                           ,LEAD_WRITE_TASK)
!
          FIELDSIZE=(DIM1+2*NBDR)*(DIM2+2*NBDR)
          ALLOCATE(TMP(FIELDSIZE))
        ENDIF
!
      ENDIF
!
!-----------------------------------------------------------------------
!***  We will now assemble the full domain 2-D restart data onto
!***  the lead Write task from the subsections on all Write tasks
!***  then the lead task will write each 2D Field to the restart
!***  file(s). 
!
!***  NOTE:  The lead Write task assembles and writes to restart only
!***         one 2-D Field at a time.  
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  First loop through all of the integer Fields.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
      rst_field_loop_int: DO NFIELD=1,RST_KOUNT_I2D                        !<-- Loop through all 2D integer gridded restart data
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
        IF(MYPE>LEAD_WRITE_TASK)THEN                                       !<-- All write tasks except the lead one
!-----------------------------------------------------------------------
!
          JSTA_WRITE=wrt_int_state%LOCAL_JSTART(wrt_int_state%ID_FTASK_RECV_STA(MYPE))  !<-- Starting J of this write task's subsection
          JEND_WRITE=wrt_int_state%LOCAL_JEND  (wrt_int_state%ID_FTASK_RECV_END(MYPE))  !<-- Ending J of this write task's subsection
!
          NN=0
!
          DO J=JSTA_WRITE,JEND_WRITE
          DO I=1,IM
            NN=NN+1
            wrt_int_state%RST_BUFF_INT(NN)=wrt_int_state%RST_WRITE_SUBSET_I(I,J,NFIELD)
          ENDDO
          ENDDO
!
          CALL MPI_RECV(ID_DUMMY                                        &  !<-- Blocking Recv keeps the following sends in line
                       ,1                                               &  !<-- Length of ID_DUMMY
                       ,MPI_INTEGER                                     &  !<-- Datatype
                       ,0                                               &  !<-- The lead write task sent this 
                       ,0                                               &  !<-- An MPI tag
                       ,MPI_COMM_COMP                                   &  !<-- The communicator
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR )
!
          MY_LOCAL_ID=MYPE-LAST_FCST_TASK-1                                !<-- This write task's local ID (between 0 and NWTPG-1)
!
          CALL MPI_SEND(wrt_int_state%RST_BUFF_INT                      &  !<-- Send this string of subsection data 
                       ,NN                                              &  !<-- Number of words sent
                       ,MPI_INTEGER                                     &  !<-- Datatype
                       ,0                                               &  !<-- Send the data to the lead write task with local ID of 0
                       ,MY_LOCAL_ID                                     &  !<-- An MPI tag
                       ,MPI_COMM_COMP                                   &  !<-- MPI communicator
                       ,IERR )
!
!-----------------------------------------------------------------------
        ELSEIF(MYPE==LEAD_WRITE_TASK)THEN                                  !<-- The lead write task
!-----------------------------------------------------------------------
!
!         write(0,*)'RST lead_write_task, send and recv,LAST_WRITE_TASK=',  &
!          'LEAD_WRITE_TASK=',LEAD_WRITE_TASK
!
          JSTA_WRITE=wrt_int_state%LOCAL_JSTART(wrt_int_state%ID_FTASK_RECV_STA(MYPE))  !<-- Starting J of lead write task's subsection
          JEND_WRITE=wrt_int_state%LOCAL_JEND  (wrt_int_state%ID_FTASK_RECV_END(MYPE))  !<-- Ending J of lead write task's subsection
!
          DO J=JSTA_WRITE,JEND_WRITE
          DO I=1,IM
            wrt_int_state%RST_OUTPUT_ARRAY_I2D(I,J)=wrt_int_state%RST_WRITE_SUBSET_I(I,J,NFIELD) !<-- Lead write task fills its part
                                                                                                 !    of full domain
          ENDDO
          ENDDO
!
          IF(LAST_WRITE_TASK>LEAD_WRITE_TASK)THEN                          !<-- Recv output subsections if more than 1 write task
            DO N=1,NWTPG-1                                                 !<-- Loop through local IDs of all other write tasks
                                                                           !    that send to the lead task
!
              CALL MPI_SEND(N                                           &  !<-- Send to other write tasks to keep their sends in line
                           ,1                                           &  !<-- Number of words sent
                           ,MPI_INTEGER                                 &  !<-- Datatype
                           ,N                                           &  !<-- Send to each of the other write tasks
                           ,0                                           &  !<-- An MPI tag
                           ,MPI_COMM_COMP                               &  !<-- MPI communicator
                           ,IERR )
!
              CALL MPI_RECV(wrt_int_state%RST_BUFF_INT                  &  !<-- Recv string of subsection data from other write tasks
                           ,IM*JM                                       &  !<-- Maximum number of words sent
                           ,MPI_INTEGER                                 &  !<-- Datatype
                           ,N                                           &  !<-- Recv from this write task
                           ,N                                           &  !<-- An MPI tag
                           ,MPI_COMM_COMP                               &  !<-- MPI communicator
                           ,JSTAT                                       &  !<-- MPI status object
                           ,IERR )
!
              NN=0
              JSTA_WRITE=wrt_int_state%LOCAL_JSTART(wrt_int_state%ID_FTASK_RECV_STA(N+LEAD_WRITE_TASK)) !<-- Starting J of sending write task
              JEND_WRITE=wrt_int_state%LOCAL_JEND  (wrt_int_state%ID_FTASK_RECV_END(N+LEAD_WRITE_TASK)) !<-- Ending J of sending write task
!
              DO J=JSTA_WRITE,JEND_WRITE
              DO I=1,IM
                NN=NN+1
                wrt_int_state%RST_OUTPUT_ARRAY_I2D(I,J)=wrt_int_state%RST_BUFF_INT(NN)   !<-- Insert other write tasks' subsections
                                                                                         !    into full domain
              ENDDO
              ENDDO
!
            ENDDO
          ENDIF
!
         NPOSN_1=(NFIELD-1)*ESMF_MAXSTR+1
         NPOSN_2=NFIELD*ESMF_MAXSTR
         NAME=wrt_int_state%RST_NAMES_I2D_STRING(NPOSN_1:NPOSN_2)                       !<-- The name of this 2D integer restart quantity
!
         IF(wrt_int_state%WRITE_RST_BIN)THEN
!
          WRITE(wrt_int_state%IO_RST_UNIT,iostat=RC)wrt_int_state%RST_OUTPUT_ARRAY_I2D   !<-- Lead write task writes out the 2D integer data
!
          IF(RST_FIRST .AND. (wrt_int_state%PRINT_OUTPUT .OR. wrt_int_state%PRINT_ALL) )THEN
            WRITE(0,*)'Wrote ',TRIM(NAME),' to restart file unit ',wrt_int_state%IO_RST_UNIT, &
             maxval(wrt_int_state%RST_OUTPUT_ARRAY_I2D),minval(wrt_int_state%RST_OUTPUT_ARRAY_I2D)
          ENDIF
         ENDIF
!
!-----------------------------------------------------------------------
!***  For the NEMSIO restart file
!-----------------------------------------------------------------------
!
         IF(wrt_int_state%WRITE_RST_NEMSIO)THEN
!
           IF(FIELDSIZE/=IM*JM)THEN
             WRITE(0,*)'WRONG: input data dimension ',IM*JM,                &
              ' does not match data size in NEMSIO file ',FIELDSIZE
           ENDIF
           TMP=RESHAPE(wrt_int_state%RST_OUTPUT_ARRAY_I2D(1:IM,1:JM),(/FIELDSIZE/))
!
           CALL NEMSIO_WRITEREC(NEMSIOFILE,NFIELD,TMP,IRET=IERR)           !<-- Lead write task writes out the 2D int data!
           IF(IERR/=0)THEN
             WRITE(0,*)' Failed to write output to file! Aborting!'
             CALL ESMF_Finalize(rc=RC,endflag=ESMF_END_ABORT)
           ENDIF
!
         ENDIF
!-----------------------------------------------------------------------
!
        ENDIF
!
!-----------------------------------------------------------------------
!
      ENDDO rst_field_loop_int
!
!-----------------------------------------------------------------------
!***  Now loop through all the real Fields.
!-----------------------------------------------------------------------
!
      WRITE(MODEL_LEVEL,'(I3.3)')wrt_int_state%LM(1)
!
!-----------------------------------------------------------------------
!
      rst_field_loop_real: DO NFIELD=1,wrt_int_state%RST_KOUNT_R2D(1)      !<-- Loop through all 2D real gridded restart data
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
        IF(MYPE>LEAD_WRITE_TASK)THEN                                       !<-- All write tasks except the lead one
!-----------------------------------------------------------------------
!
          JSTA_WRITE=wrt_int_state%LOCAL_JSTART(wrt_int_state%ID_FTASK_RECV_STA(MYPE))  !<-- Starting J of this write task's subsection
          JEND_WRITE=wrt_int_state%LOCAL_JEND  (wrt_int_state%ID_FTASK_RECV_END(MYPE))  !<-- Ending J of this write task's subsection
!
          NN=0
!
          DO J=JSTA_WRITE,JEND_WRITE
          DO I=1,IM
            NN=NN+1
            wrt_int_state%RST_BUFF_REAL(NN)=wrt_int_state%RST_WRITE_SUBSET_R(I,J,NFIELD)
          ENDDO
          ENDDO
!
          CALL MPI_RECV(ID_DUMMY                                        &  !<-- Blocking Recv keeps the following sends in line
                       ,1                                               &  !<-- Length of ID_DUMMY
                       ,MPI_INTEGER                                     &  !<-- Datatype
                       ,0                                               &  !<-- The lead write task sent this 
                       ,0                                               &  !<-- An MPI tag
                       ,MPI_COMM_COMP                                   &  !<-- The communicator
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR )
!
          MY_LOCAL_ID=MYPE-LAST_FCST_TASK-1                                !<-- This write task's local ID (between 0 and NWTPG-1)
!
          CALL MPI_SEND(wrt_int_state%RST_BUFF_REAL                     &  !<-- Send this string of subsection data 
                       ,NN                                              &  !<-- Number of words sent
                       ,MPI_REAL                                        &  !<-- Datatype
                       ,0                                               &  !<-- Send the data to the lead write task with local ID of 0
                       ,MY_LOCAL_ID                                     &  !<-- An MPI tag
                       ,MPI_COMM_COMP                                   &  !<-- MPI communicator
                       ,IERR )
!
!-----------------------------------------------------------------------
        ELSEIF(MYPE==LEAD_WRITE_TASK)THEN                                  !<-- The lead write task
!-----------------------------------------------------------------------
!
          JSTA_WRITE=wrt_int_state%LOCAL_JSTART(wrt_int_state%ID_FTASK_RECV_STA(MYPE))  !<-- Starting J of lead write task's subsection
          JEND_WRITE=wrt_int_state%LOCAL_JEND  (wrt_int_state%ID_FTASK_RECV_END(MYPE))  !<-- Ending J of lead write task's subsection
!
          DO J=JSTA_WRITE,JEND_WRITE
          DO I=1,IM
            wrt_int_state%RST_OUTPUT_ARRAY_R2D(I,J)=wrt_int_state%RST_WRITE_SUBSET_R(I,J,NFIELD) !<-- Lead write task fills its part
                                                                                                 !    of full domain
          ENDDO
          ENDDO
!
          IF(LAST_WRITE_TASK>LEAD_WRITE_TASK)THEN                          !<-- Recv output subsections if more than 1 write task
            DO N=1,NWTPG-1                                                 !<-- Loop through local IDs of all other write tasks
                                                                           !    that send to the lead task
!
              CALL MPI_SEND(N                                           &  !<-- Send to other write tasks to keep their sends in line
                           ,1                                           &  !<-- Number of words sent
                           ,MPI_INTEGER                                 &  !<-- Datatype
                           ,N                                           &  !<-- Send to each of the other write tasks
                           ,0                                           &  !<-- An MPI tag
                           ,MPI_COMM_COMP                               &  !<-- MPI communicator
                           ,IERR )
!
              CALL MPI_RECV(wrt_int_state%RST_BUFF_REAL                 &  !<-- Recv string of subsection data from other write tasks
                           ,IM*JM                                       &  !<-- Maximum number of words sent
                           ,MPI_REAL                                    &  !<-- Datatype
                           ,N                                           &  !<-- Recv from this write task
                           ,N                                           &  !<-- An MPI tag
                           ,MPI_COMM_COMP                               &  !<-- MPI communicator
                           ,JSTAT                                       &  !<-- MPI status object
                           ,IERR )
!
              NN=0
              JSTA_WRITE=wrt_int_state%LOCAL_JSTART(wrt_int_state%ID_FTASK_RECV_STA(N+LEAD_WRITE_TASK)) !<-- Starting J of sending write task
              JEND_WRITE=wrt_int_state%LOCAL_JEND  (wrt_int_state%ID_FTASK_RECV_END(N+LEAD_WRITE_TASK)) !<-- Ending J of sending write task
!
              DO J=JSTA_WRITE,JEND_WRITE
              DO I=1,IM
                NN=NN+1
                wrt_int_state%RST_OUTPUT_ARRAY_R2D(I,J)=wrt_int_state%RST_BUFF_REAL(NN)   !<-- Insert other write tasks' subsections
                                                                                          !    into full domain
              ENDDO
              ENDDO
!
            ENDDO
          ENDIF
!
!
          NPOSN_1=(NFIELD-1)*ESMF_MAXSTR+1
          NPOSN_2=NFIELD*ESMF_MAXSTR
          NAME=wrt_int_state%RST_NAMES_R2D_STRING(NPOSN_1:NPOSN_2)                       !<-- The name of this 2D real restart quantity
!
!-----------------------------------------------------------------------
!***  Begin computation of the 10-m wind factor for GSI.
!-----------------------------------------------------------------------
!
          IF(TRIM(NAME)=='U10') THEN
            IF(.NOT.ALLOCATED(FACT10)) THEN
              ALLOCATE(FACT10(1:IM,1:JM))
!
              DO J=1,JM
              DO I=1,IM
                FACT10(I,J)=0.
              ENDDO
              ENDDO
            ENDIF
!
            DO J=1,JM
            DO I=1,IM
              FACT10(I,J)=FACT10(I,J)+                                 &
                          wrt_int_state%RST_OUTPUT_ARRAY_R2D(I,J)*     &
                          wrt_int_state%RST_OUTPUT_ARRAY_R2D(I,J)
            ENDDO
            ENDDO
          ENDIF
!
          IF(TRIM(NAME)=='V10') THEN
            IF(.NOT.ALLOCATED(FACT10)) THEN
              ALLOCATE(FACT10(1:IM,1:JM))
              FACT10=0.
            ENDIF
!
            DO J=1,JM
            DO I=1,IM
             FACT10(I,J)=FACT10(I,J)+                                  &
                         wrt_int_state%RST_OUTPUT_ARRAY_R2D(I,J)*      &
                         wrt_int_state%RST_OUTPUT_ARRAY_R2D(I,J)
            ENDDO
            ENDDO
          ENDIF
!
          IF(TRIM(NAME)=='U_'//MODEL_LEVEL//'_2D') THEN
            ALLOCATE(FACT10TMPU(1:IM,1:JM))
            CALL V_TO_H_BGRID(wrt_int_state%RST_OUTPUT_ARRAY_R2D(1:IM,1:JM) &
                             ,IM,JM,GLOBAL,FACT10TMPU)
          ENDIF
!
          IF(TRIM(NAME)=='V_'//MODEL_LEVEL//'_2D') THEN
            ALLOCATE(FACT10TMPV(1:IM,1:JM))
            CALL V_TO_H_BGRID(wrt_int_state%RST_OUTPUT_ARRAY_R2D(1:IM,1:JM) &
                             ,IM,JM,GLOBAL,FACT10TMPV)
          ENDIF

!
          IF(wrt_int_state%WRITE_RST_BIN)THEN
!
            WRITE(wrt_int_state%IO_RST_UNIT,iostat=RC)wrt_int_state%RST_OUTPUT_ARRAY_R2D   !<-- Lead write task writes out the 2D real data
!
            IF(RST_FIRST .AND. (wrt_int_state%PRINT_OUTPUT .OR. wrt_int_state%PRINT_ALL) )THEN
              WRITE(0,*)'Wrote ',TRIM(NAME)                                &
                       ,' to restart file unit ',wrt_int_state%IO_RST_UNIT &
                       ,MAXVAL(wrt_int_state%RST_OUTPUT_ARRAY_R2D)         &
                       ,MINVAL(wrt_int_state%RST_OUTPUT_ARRAY_R2D)
            ENDIF
!
          ENDIF
!
!-----------------------------------------------------------------------
!***  The same for the NEMSIO restart file.
!-----------------------------------------------------------------------
!
          IF(wrt_int_state%WRITE_RST_NEMSIO)THEN
!
            IF(FIELDSIZE/=IM*JM)THEN
              WRITE(0,*)'WRONG: data dimension ',IM*JM,                 &
               ' does not match data size in NEMSIO file,',FIELDSIZE
            ENDIF
            IF(TRIM(NAME)=='FIS') THEN
              IF(.NOT.ALLOCATED(HGT)) ALLOCATE(HGT(1:IM,1:JM))
              HGT(1:IM,1:JM)=wrt_int_state%RST_OUTPUT_ARRAY_R2D(1:IM,1:JM)/G
            ENDIF
!
            IF(TRIM(NAME)=='GLAT')THEN
               ALLOCATE(GLAT1D(FIELDSIZE))
               GLAT1D(1:FIELDSIZE)=RESHAPE(wrt_int_state%RST_OUTPUT_ARRAY_R2D(1:IM,1:JM),(/FIELDSIZE/))
            ENDIF
!
            IF(TRIM(NAME)=='GLON')THEN
               ALLOCATE(GLON1D(FIELDSIZE))
               GLON1D(1:FIELDSIZE)=RESHAPE(wrt_int_state%RST_OUTPUT_ARRAY_R2D(1:IM,1:JM),(/FIELDSIZE/))
            ENDIF
!
            N=NFIELD+wrt_int_state%RST_KOUNT_I2D(1)
            TMP=RESHAPE(wrt_int_state%RST_OUTPUT_ARRAY_R2D(1:IM,1:JM),(/FIELDSIZE/))
!
            CALL NEMSIO_WRITEREC(NEMSIOFILE,N,TMP,IRET=IERR)
            IF(IERR/=0)THEN
              WRITE(0,*)' Failed to write output to file! Aborting!'
              CALL ESMF_Finalize(rc=RC,endflag=ESMF_END_ABORT)
            ENDIF
!
!           IF(RST_FIRST)THEN
!             WRITE(0,*)'Wrote ',TRIM(NAME),' to nemsio restart file iret=',ierr
!           ENDIF

          ENDIF
!
!-----------------------------------------------------------------------
!
        ENDIF
!
!-----------------------------------------------------------------------
!
      ENDDO rst_field_loop_real
!
!-----------------------------------------------------------------------
!***  Complete the computation of 10-m wind factor and write out 
!***  FACT10 AND HGT.
!-----------------------------------------------------------------------
!
      IF(MYPE==LEAD_WRITE_TASK) THEN
        IF(ALLOCATED(FACT10).AND.ALLOCATED(FACT10TMPU).AND.ALLOCATED(FACT10TMPV)) THEN
!
          DO J=1,JM
          DO I=1,IM
            FACT10TMPV(I,J)=SQRT(FACT10TMPU(I,J)*FACT10TMPU(I,J)+       &
                                 FACT10TMPV(I,J)*FACT10TMPV(I,J) )
!
            IF(FACT10TMPV(I,J)/=0) THEN
              FACT10(I,J)=SQRT(FACT10(I,J))/FACT10TMPV(I,J)
            ELSE
              FACT10(I,J)=1.
            ENDIF
          ENDDO
          ENDDO
!
          DEALLOCATE(FACT10TMPU)
          DEALLOCATE(FACT10TMPV)
!
          IF(wrt_int_state%WRITE_RST_BIN)THEN
            WRITE(wrt_int_state%IO_RST_UNIT,iostat=RC)FACT10                !<-- Lead write task writes out the 2D real data
            IF(RST_FIRST .AND. (wrt_int_state%PRINT_OUTPUT .OR. wrt_int_state%PRINT_ALL) )THEN
              WRITE(0,*)'Wrote FACT10 to restart file unit ',wrt_int_state%IO_RST_UNIT &
                        ,maxval(fact10),minval(fact10)
            ENDIF
          ENDIF

          IF(wrt_int_state%WRITE_RST_NEMSIO)THEN
            N=N+1
            TMP=RESHAPE(FACT10(1:IM,1:JM),(/FIELDSIZE/))
            CALL NEMSIO_WRITEREC(NEMSIOFILE,N,TMP,IRET=IERR)
            IF(IERR/=0)THEN
              WRITE(0,*)' Failed to write output to file! Aborting!'
              CALL ESMF_Finalize(rc=RC,endflag=ESMF_END_ABORT)
            ENDIF
!           write(0,*)'after nemsio_writerec,n=',n,'fact10=',maxval(tmp),minval(tmp),'iret=',ierr
!           WRITE(0,*)'Wrote FACT10 to nemsio restart file iret=',ierr
          ENDIF
!
          DEALLOCATE(FACT10)

        ENDIF
!
        IF( wrt_int_state%WRITE_RST_NEMSIO)THEN
          N=N+1
          TMP=RESHAPE(HGT(1:IM,1:JM),(/FIELDSIZE/))
          CALL NEMSIO_WRITEREC(NEMSIOFILE,N,TMP,IRET=IERR)
          IF(IERR/=0)THEN
            WRITE(0,*)' Failed to write output to file! Aborting!'
            CALL ESMF_Finalize(rc=RC,endflag=ESMF_END_ABORT)
          ENDIF
          DEALLOCATE(HGT)
        ENDIF
!
      ENDIF
!
      RST_FIRST=.FALSE.
!
!-----------------------------------------------------------------------
!***  Close the restart file.
!-----------------------------------------------------------------------
!
      IF(wrt_int_state%WRITE_RST_BIN.and.MYPE==LEAD_WRITE_TASK)THEN
        CLOSE(wrt_int_state%IO_RST_UNIT)
!       write(0,*)' Closed restart file with unit=',wrt_int_state%IO_RST_UNIT
      ENDIF
!
!-----------------------------------------------------------------------
!***  Close the NEMSIO restart file.
!-----------------------------------------------------------------------
!
      IF(wrt_int_state%WRITE_RST_NEMSIO.AND.wrt_int_state%MYPE==LEAD_WRITE_TASK)THEN
!
        IF(ASSOCIATED(GLAT1D).AND.ASSOCIATED(GLON1D)) THEN
          CALL NEMSIO_SETFILEHEAD(NEMSIOFILE,IERR,GLAT1D,GLON1D)
          DEALLOCATE(GLAT1D,GLON1D)
        ENDIF
!  
        CALL NEMSIO_GETFILEHEAD(NEMSIOFILE,IERR,gfname=GFNAME)
!
        DEALLOCATE(TMP)
!
        CALL NEMSIO_CLOSE(NEMSIOFILE)
!       WRITE(0,*)' Closed nemsio_restart file, ', gfname
!
        CALL NEMSIO_FINALIZE()
!
      ENDIF
!
!----------------------
!***  Restartdone file
!----------------------
!
!-----------------------------------------------------------------------
      IF(wrt_int_state%WRITE_DONEFILEFLAG                               &
           .and.                                                        &
         wrt_int_state%MYPE==LEAD_WRITE_TASK                            &
           .and.                                                        &
         (wrt_int_state%WRITE_RST_BIN                                   &
           .or.                                                         &
          wrt_int_state%WRITE_RST_NEMSIO))THEN
!
        CALL ESMF_ConfigGetAttribute(config=CF                          &  !<-- The configure file object
                                    ,value =ID_DOMAIN                   &  !<-- Put extracted quantity here
                                    ,label ='my_domain_id:'             &  !<-- The quantity's label in the configure file
                                    ,rc    =RC)
!
        INT_SEC=INT(NF_SECONDS)
        FRAC_SEC=NINT((NF_SECONDS-INT_SEC)*100.)
!
        WRITE(FILENAME,'(A,I2.2,A,I4.4,A,I2.2,A,I2.2,A,I2.2,A)' )       &
                                         'restartdone.',ID_DOMAIN,'.'   &  !<-- Insert domain ID for nests
                                         ,NF_HOURS,'h_'                 &  
                                         ,NF_MINUTES,'m_'               &
                                         ,INT_SEC,'.',FRAC_SEC,'s'
!
        DO N=51,99
          INQUIRE(N,opened=OPENED)
          IF(.NOT.OPENED)THEN
            IO_RST_UNIT=N
            EXIT
          ENDIF
        ENDDO
!
        OPEN(unit  =IO_RST_UNIT                                         &
            ,file  =trim(FILENAME)                                      &
            ,form  ='formatted'                                         &
            ,status='REPLACE')
!
        WRITE(IO_RST_UNIT,'(A4)')'DONE'
        CLOSE(IO_RST_UNIT)
!
      ENDIF
!
!-----------------------------------------------------------------------
!
      ENDIF restart_time
!
!-----------------------------------------------------------------------
!
      IF(RC_RUN==ESMF_SUCCESS)THEN
!       WRITE(0,*)"PASS: WRITE_RUN"
      ELSE
        WRITE(0,*)"FAIL: WRITE_RUN"
      ENDIF
!
      write_run_tim=write_run_tim+(timef()-btim0)
!
      IF(MYPE==LEAD_WRITE_TASK)THEN
        IF(wrt_int_state%PRINT_OUTPUT .OR. wrt_int_state%PRINT_ALL)   &
        WRITE(0,*)' Write Time is ',write_run_tim*1.e-3               &
                 ,' at Fcst ',NF_HOURS,':',NF_MINUTES,':',NF_SECONDS
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE WRITE_RUN
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE WRITE_FINALIZE(WRITE_COMP                              &
                               ,IMP_STATE_WRITE                         &
                               ,EXP_STATE_WRITE                         &
                               ,CLOCK                                   &
                               ,RCFINAL)
!
!-----------------------------------------------------------------------
!***  Finalize the Write gridded component.
!-----------------------------------------------------------------------
!
!***  HISTORY
!       xx Feb 2007:  W. Yang - Originator
!       13 Jun 2007:  T. Black - Name revisions
!
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(ESMF_GridComp) :: WRITE_COMP                                    !<-- The Write component
!
      TYPE(ESMF_State) :: IMP_STATE_WRITE                               &  !<-- The Write component import state 
                         ,EXP_STATE_WRITE                                  !<-- The Write component export state
!
      TYPE(ESMF_Clock) :: CLOCK                                            !<-- The Write component Clock
!
      INTEGER,INTENT(OUT) :: RCFINAL
!
!---------------------
!***  Local Variables
!---------------------
!
      TYPE(WRITE_WRAP) :: WRAP
!
      INTEGER :: RC
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC     =ESMF_SUCCESS
      RCFINAL=ESMF_SUCCESS
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Retrieve Write Component's Internal State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompGetInternalState(WRITE_COMP                     &  !<-- The write component
                                        ,WRAP                           &  !<-- Pointer to internal state
                                        ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RCFINAL)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      DEALLOCATE(WRAP%WRITE_INT_STATE)
!
!-----------------------------------------------------------------------
!
      IF(RCFINAL==ESMF_SUCCESS)THEN
!       WRITE(0,*)'PASS: Write_Finalize.'
      ELSE
        WRITE(0,*)'FAIL: Write_Finalize.'
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE WRITE_FINALIZE
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE WRITE_SETUP(DOMAIN_GRID_COMP                           &
                            ,DOMAIN_INT_STATE                           &
                            ,CLOCK_DOMAIN)
! 
!-----------------------------------------------------------------------
!***  Set up the Write components with the forecast tasks and
!***  the groups of write tasks needed for quilting the output
!***  and writing it to history/restart files.
!-----------------------------------------------------------------------
!
      USE MODULE_DOMAIN_INTERNAL_STATE,ONLY: DOMAIN_INTERNAL_STATE      
!
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: DOMAIN_GRID_COMP                !<-- The DOMAIN component
!
      TYPE(DOMAIN_INTERNAL_STATE),INTENT(INOUT) :: DOMAIN_INT_STATE        !<-- The DOMAIN component's Internal State
!
      TYPE(ESMF_Clock),INTENT(INOUT) :: CLOCK_DOMAIN                       !<-- The DOMAIN component's ESMF Clock
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: INPES,JNPES                                 &  !<-- Number of fcst tasks in I and J directions
                           ,MYPE                                        &  !<-- My task ID
                           ,WRITE_GROUPS                                &  !<-- Number of groups of write tasks
                           ,WRITE_TASKS_PER_GROUP                          !<-- #of tasks in each write group
!
      INTEGER(kind=KINT) :: I,ISTAT,J,K,RC,RC_SETUP
!
      CHARACTER(2) :: MY_WRITE_GROUP
      CHARACTER(6) :: FMT='(I2.2)'
      CHARACTER(ESMF_MAXSTR) :: WRITE_NAME
!
      TYPE(ESMF_VM)          :: VM                                         !<-- The ESMF virtual machine.
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Save the total number of domains.  This is always given in
!***  domain #1's configure file.
!-----------------------------------------------------------------------
!
      CF_1=ESMF_ConfigCreate(rc=RC)                                        !<-- Create the confure object for domain #1
      CONFIGFILE_01_NAME='configure_file_01'
!
      CALL ESMF_ConfigLoadFile(config  =CF_1                            &  !<-- Load the configure file object for domain #1
                              ,filename=CONFIGFILE_01_NAME              &  !<-- The configure file's name
                              ,rc      =RC)
!
      CALL ESMF_ConfigGetAttribute(config=CF_1                          &  !<-- The configure file object for domain #1
                                  ,value =NUM_DOMAINS_TOTAL             &  !<-- Extract the total # of domains
                                  ,label ='num_domains_total:'          &  !<-- The quantity's label in the configure file
                                  ,rc    =RC)
!
!-----------------------------------------------------------------------
!***  Retrieve the config object CF from the DOMAIN component for
!***  this domain.  Save it in the Write component's internal state.
!***  This will allow a given task to reference it for any of the
!***  multiple domains a given task may lie on.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Get Config Object for Write Setup"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompGet(gridcomp=DOMAIN_GRID_COMP                   &  !<-- The DOMAIN component
                           ,config  =CF                                 &  !<-- The configure object on this domain
                           ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_SETUP)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Retrieve task and group counts from the config file.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Write Tasks and Groups from Config File"
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
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_SETUP)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      domain_int_state%WRITE_GROUPS=WRITE_GROUPS                           !<-- Save for the DOMAIN's Finalize step
      domain_int_state%WRITE_TASKS_PER_GROUP=WRITE_TASKS_PER_GROUP         !<-- Save for the DOMAIN's Finalize step
!
!-----------------------------------------------------------------------
!***  How many forecast tasks do we have?
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="INPES/JNPES from Config Object for Write Setup"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The ESMF configure object
                                  ,value =INPES                         &  !<-- # of fcst tasks in I direction
                                  ,label ='inpes:'                      &  !<-- Give the value of this label to INPES
                                  ,rc    =RC)
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The ESMF configure object
                                  ,value =JNPES                         &  !<-- # of fcst tasks in J direction
                                  ,label ='jnpes:'                      &  !<-- Give the value of this label to JNPES
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_SETUP)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      NUM_PES_FCST=INPES*JNPES                                             !<-- Total number of forecast tasks
!
!-----------------------------------------------------------------------
!***  Retrieve the current VM.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Retrieve the Local VM"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_VMGetCurrent(vm=VM                                      &  !<-- The ESMF virtual machine
                            ,rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_SETUP)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  What is my MPI task ID?
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Get MPI Task IDs for Write Setup"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_VMGet(vm      =VM                                       &  !<-- The virtual machine
                     ,localpet=MYPE                                     &  !<-- Local PE rank
                     ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_SETUP)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Associate all of the forecast tasks with the write tasks
!***  in each write group.
!-----------------------------------------------------------------------
!
      ALLOCATE(domain_int_state%PETLIST_WRITE(NUM_PES_FCST+WRITE_TASKS_PER_GROUP,WRITE_GROUPS)) !<-- Task IDs of all fcst tasks
                                                                                                !    plus the write tasks
                                                                                                !    by write group
!
!-----------------------------------------------------------------------
!***  Collect the task IDs for the write tasks and the associated
!***  forecast-write tasks.
!-----------------------------------------------------------------------
!
      DO I=0,NUM_PES_FCST-1
        DO J=1,WRITE_GROUPS
          domain_int_state%PETLIST_WRITE(I+1,J)=I                          !<-- Collect forecast task IDs to be associated with
                                                                           !    write tasks by write group
        ENDDO
!
      ENDDO
!
      K=NUM_PES_FCST
!
      DO J=1,WRITE_GROUPS
        DO I=1,WRITE_TASKS_PER_GROUP
          domain_int_state%PETLIST_WRITE(NUM_PES_FCST+I,J)=K               !<-- Append Write task IDs to associated forecast task IDs by group
          K=K+1
        ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
!***  Create the Write gridded component(s).
!***  There are as many Write components as there are groups of
!***  write tasks specified in the configure file.
!***  Register their Init, Run, and Finalize steps.
!-----------------------------------------------------------------------
!
      IF(WRITE_GROUPS>0)THEN
        ALLOCATE(domain_int_state%WRITE_COMPS(WRITE_GROUPS))               !<-- The Write gridded components
      ENDIF
!
!---------------------------------
!***  Create the Write components
!---------------------------------
!
      DO I=1,WRITE_GROUPS
        WRITE(MY_WRITE_GROUP,FMT)I
        WRITE_NAME='write_GridComp_'//MY_WRITE_GROUP
!
        domain_int_state%WRITE_COMPS(I)=ESMF_GridCompCreate(                       &
                                       name   =WRITE_NAME                          &  !<-- Name of this group's Write gridded component
                                      ,config =CF                                  &  !<-- The configure file for writes
                                      ,petList=domain_int_state%PETLIST_WRITE(:,I) &  !<-- The task IDs of the write tasks in this group
                                                                                      !    provide the local VM information per component.
                                      ,rc     =RC)
!
      ENDDO
!
!-----------------------------------
!***  Register the Write components
!-----------------------------------
!
      DO I=1,WRITE_GROUPS
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Register Write Components"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_GridCompSetServices(domain_int_state%WRITE_COMPS(I)    &  !<-- The Write gridded components
                                     ,WRITE_REGISTER                     &  !<-- The user's subroutine name
                                     ,rc = RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_SETUP)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDDO
!
!------------------------------------------------------------------------
!***  Create empty Import and Export states for the Write subcomponent(s)
!------------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create Empty Import/Export States for Write Components"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      domain_int_state%IMP_STATE_WRITE=ESMF_StateCreate(name  ='Write Import State' &  !<-- Import state name for writes
                                                       ,stateintent= ESMF_STATEINTENT_IMPORT &
                                                       ,rc         = RC)
!
      domain_int_state%EXP_STATE_WRITE=ESMF_StateCreate(name  ='Write Export State' &  !<-- Export state names for writes
                                                       ,stateintent= ESMF_STATEINTENT_EXPORT &
                                                       ,rc         = RC) 
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_SETUP)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Insert the Write components' import state into the
!***  Solver export state since history data itself must 
!***  come from the Solver.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Insert Write Import State into Solver Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateAddReplace(domain_int_state%EXP_STATE_SOLVER               & !<-- Solver export state receives a state
                        ,(/domain_int_state%IMP_STATE_WRITE/)   & !<-- Add the write components' import state
                        ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_SETUP)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE WRITE_SETUP
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE WRITE_DESTROY(DOMAIN_GRID_COMP                         &
                              ,DOMAIN_INT_STATE                         &
                              ,CLOCK_DOMAIN)
! 
!-----------------------------------------------------------------------
!***  Destroy all objects related to the Write components.
!-----------------------------------------------------------------------
!
      USE MODULE_DOMAIN_INTERNAL_STATE,ONLY: DOMAIN_INTERNAL_STATE      
!
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: DOMAIN_GRID_COMP                !<-- The DOMAIN component
!
      TYPE(DOMAIN_INTERNAL_STATE),INTENT(INOUT) :: DOMAIN_INT_STATE        !<-- The Domain component's Internal State
!
      TYPE(ESMF_Clock),INTENT(INOUT) :: CLOCK_DOMAIN                       !<-- The DOMAIN component's ESMF Clock
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER :: I,J,MYPE,N,RC,RC_DES
!
      TYPE(ESMF_VM) :: VM
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC    =ESMF_SUCCESS
      RC_DES=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  Retrieve the current VM.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Retrieve the Local VM in Write Destroy"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_VMGetCurrent(vm=VM                                      &  !<-- The ESMF virtual machine
                            ,rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_DES)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  What is my MPI task ID?
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Get MPI Task IDs for Write Destroy"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_VMGet(vm      =VM                                       &  !<-- The virtual machine
                     ,localpet=MYPE                                     &  !<-- Local PE rank
                     ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_DES)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Finalize the Write gridded components in each write group.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Finalize Write Components"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      DO N=1,domain_int_state%WRITE_GROUPS              
        IF(MYPE>=domain_int_state%PETLIST_WRITE(1,N).AND.                          &
           MYPE<=domain_int_state%PETLIST_WRITE(domain_int_state%WRITE_TASKS_PER_GROUP,N))THEN
!
           CALL ESMF_GridCompFinalize(gridcomp   =domain_int_state%WRITE_COMPS(N)  &
                                     ,importstate=domain_int_state%EXP_STATE_WRITE &
                                     ,exportstate=domain_int_state%IMP_STATE_WRITE &
                                     ,clock      =CLOCK_DOMAIN                     &
                                     ,rc         =RC)
        ENDIF
      ENDDO
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_DES)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Destroy the Write components' import/export states.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Destroy Write Component Import/Export States"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateDestroy(domain_int_state%IMP_STATE_WRITE,rc=RC)
      CALL ESMF_StateDestroy(domain_int_state%EXP_STATE_WRITE,rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_DES)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Destroy the Write components.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Destroy Write Components"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      DO J=1,domain_int_state%WRITE_GROUPS
!
        CALL ESMF_VMBarrier(vm=VM,rc=RC)
!
        DO I=1,domain_int_state%NUM_PES_FCST+domain_int_state%WRITE_TASKS_PER_GROUP
!
          IF(MYPE==domain_int_state%PETLIST_WRITE(I,J))THEN
            CALL ESMF_GridCompDestroy(gridcomp=domain_int_state%WRITE_COMPS(J) &
                                     ,rc      =RC)
          ENDIF
!
        ENDDO
!
      ENDDO
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_DES)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  The final error signal information.
!-----------------------------------------------------------------------
!
      IF(RC_DES==ESMF_SUCCESS)THEN
!       WRITE(0,*)' WRITE Finalize step succeeded'
      ELSE
        WRITE(0,*)' WRITE Finalize step failed  RC_DES=',RC_DES
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE WRITE_DESTROY
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      END MODULE MODULE_WRITE_GRID_COMP
!
!-----------------------------------------------------------------------
