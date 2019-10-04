!-----------------------------------------------------------------------
!
      MODULE MODULE_WRITE_ROUTINES
!
!-----------------------------------------------------------------------
!***  This module contains routines needed by the Run step of the
!***  Write gridded component in which history output data from
!***  the forecast tasks are assembled and written to history files
!***  by the write tasks.
!-----------------------------------------------------------------------
!***
!***  HISTORY   
!***
!       06 Sep 2007:  T. Black - Created module.
!                                Moved the FIRST block from the old
!                                Write component here for clarity.
!       15 Aug 2008:  J. Wang  - Revised for NEMS-IO
!       20 Aug 2008:  J. Wang  - Output start date first instead of
!                                forecast date.
!       16 Sep 2008:  J. Wang  - WRITE_NEMSIO_RUNHISTORY_OPEN only
!                                opens file and writes metadata.
!       30 Sep 2008:  E. Colon - Generalize counts for nemsio
!       14 Oct 2008:  R. Vasic - Added restart capability
!       05 Jan 2009:  J. Wang  - Added 10-m wind factor to NMMB
!                                runhistory and restart files.
!       06 Jan 2009:  T. Black - Replace Max # of words recv'd by
!                                Write tasks with actual # of words.
!       03 Sep 2009:  T. Black - Merged with NMM-B nesting code.
!       24 Mar 2010:  T. Black - Revised for NEMS restructing.
!       07 May 2010:  T. Black - Change output frequency to minutes.
!       16 Dec 2010:  J. Wang  - Change to nemsio library
!          Feb 2011:  W. Yang  - Updated to use both the ESMF 4.0.0rp2 library,
!                                ESMF 5 library and the the ESMF 3.1.0rp2 library.
!          May 2011:  J. Wang  - add run post on quilt option
!       27 SEP 2011:  W. Yang  - Modified for using the ESMF 5.2.0r library.
!
!---------------------------------------------------------------------------------
!
      USE MPI
      USE ESMF
!
      USE MODULE_WRITE_INTERNAL_STATE,ONLY: WRITE_INTERNAL_STATE        &
                                           ,WRITE_WRAP                  &
                                           ,MAX_DATA_I1D                &
                                           ,MAX_DATA_I2D                &
                                           ,MAX_DATA_R1D                &
                                           ,MAX_DATA_R2D                &
                                           ,MAX_DATA_LOG
!
      USE MODULE_DOMAIN_INTERNAL_STATE,ONLY: DOMAIN_INTERNAL_STATE
!
      USE MODULE_DM_PARALLEL,ONLY : PARA_RANGE                          &
                                   ,MPI_COMM_COMP
!
      USE MODULE_CONTROL,ONLY : TIMEF
!
      USE MODULE_KINDS
!
      USE MODULE_ERROR_MSG,ONLY: ERR_MSG,MESSAGE_CHECK
!
      USE NEMSIO_MODULE
!
      USE MODULE_CONSTANTS,ONLY : A,PI,G
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------

      PRIVATE
!
      PUBLIC :: SEND_UPDATED_ATTRIBUTES                                 &
               ,WRITE_ASYNC                                             &
               ,WRITE_INIT                                              &
               ,WRITE_NEMSIO_RUNHISTORY_OPEN                            &
               ,WRITE_NEMSIO_RUNRESTART_OPEN                            &
               ,WRITE_NEMSIOCTL                                         &
               ,WRITE_RUNHISTORY_OPEN                                   &
               ,WRITE_RUNRESTART_OPEN                                   &
               ,OPEN_HST_FILE                                           &
               ,OPEN_RST_FILE
!
!-----------------------------------------------------------------------
!
      LOGICAL,PUBLIC,SAVE :: TIME_FOR_HISTORY = .FALSE.                 &
                            ,TIME_FOR_RESTART = .FALSE.
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE PRELIM_INFO_FOR_OUTPUT(OUTPUT_FLAG                     &
                                       ,NUM_PES_FCST                    &
                                       ,NUM_WRITE_GROUPS                &
                                       ,WRITE_COMPS                     &
                                       ,IMP_STATE_WRITE )
! 
!-----------------------------------------------------------------------
!***  This routine will perform certain computations and operations
!***  that only need to be done once for each Write group.  Data that
!***  does not change with forecast time is unloaded from the Write
!***  components' import state.  The data that is unloaded consists
!***  of everything except the 2D/3D gridded forecast arrays (whose
!***  contents do change with time).  Also basic information is 
!***  provided to forecast and write tasks that will be necessary in
!***  the quilting of local 2-D gridded history data into full domain
!***  arrays.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      CHARACTER(len=7),INTENT(IN) :: OUTPUT_FLAG
!
      INTEGER(kind=KINT),INTENT(INOUT) :: NUM_PES_FCST                  &  !<-- # of forecast tasks
                                         ,NUM_WRITE_GROUPS                 !<-- # of write groups
!
      TYPE(ESMF_GridComp),DIMENSION(NUM_WRITE_GROUPS),INTENT(IN) ::     &
                                                          WRITE_COMPS      !<-- The array of Write components
!
      TYPE(ESMF_State),INTENT(INOUT) :: IMP_STATE_WRITE                    !<-- The import state of the Write components
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KDIN) :: NUM_WORDS_TOT
!
      INTEGER(kind=KINT) :: COMM_MY_DOMAIN                              &
                           ,MYPE_DOMAIN
!
      INTEGER(kind=KINT),SAVE :: ITS,ITE,JTS,JTE                        &
                                ,NCHAR_I1D                              &
                                ,NCHAR_R1D                              &
                                ,NCHAR_LOG
!
      INTEGER(kind=KINT) :: I,IERR,IM,ISTAT,J,JM,L,M1,M2                &
                           ,N,N1,N2,NN,NUM_ATTRIB,NWTPG,NX,NY           &
                           ,RC,RC_WRT
!
      INTEGER(kind=KINT) :: JROW_FIRST,JROW_LAST,JROWS                  &
                           ,LAST_FCST_TASK                              &
                           ,LOCAL_ID                                    &
                           ,N_END                                       &
                           ,N_STA
!
      INTEGER(kind=KINT) :: JEND_WRITE,JSTA_WRITE
!
      INTEGER(kind=KINT) :: KOUNT_I1D_X                                 &
                           ,KOUNT_I2D_X                                 &
                           ,KOUNT_R1D_X                                 &
                           ,KOUNT_R2D_X                                 &
                           ,KOUNT_LOG_X
!
      INTEGER(kind=KINT) :: LENGTH                                      &
                           ,LENGTH_SUM_I1D_X                            &
                           ,LENGTH_SUM_R1D_X                            &
                           ,LENGTH_SUM_LOG_X
!
      INTEGER(kind=KINT) :: MAX_WORDS                                   &
                           ,NPOSN_START,NPOSN_END                       &
                           ,NUM_FIELD_NAMES                             &
                           ,NUM_WORDS
!
      INTEGER(kind=KINT),DIMENSION(NUM_WRITE_GROUPS) :: LEAD_WRITE_TASK &
                                                       ,LAST_WRITE_TASK
!
      INTEGER(kind=KINT),DIMENSION(MPI_STATUS_SIZE) :: JSTAT
!
      INTEGER(kind=KINT),DIMENSION(:),POINTER :: IHALO,INPES            &
                                                ,JHALO,JNPES            &
                                                ,LOCAL_ISTART           &
                                                ,LOCAL_IEND             &
                                                ,LOCAL_JSTART           &
                                                ,LOCAL_JEND
!
      INTEGER(kind=KINT),DIMENSION(:),POINTER :: NCHAR_I2D              &
                                                ,NCHAR_R2D              &
                                                ,WORK_ARRAY_I1D
!
      REAL(kind=KFPT),DIMENSION(:),POINTER :: WORK_ARRAY_R1D
!
      CHARACTER(len=14) :: BUNDLE_NAME
!
      CHARACTER(ESMF_MAXSTR) :: ATTRIB_NAME
! 
      LOGICAL(kind=KLOG) :: NO_FIELDS
!
      TYPE(WRITE_WRAP) :: WRAP
!
      TYPE(WRITE_INTERNAL_STATE),POINTER :: WRT_INT_STATE
!
      TYPE(ESMF_TypeKind_Flag) :: DATATYPE
!
      TYPE(ESMF_Field) :: FIELD_WORK1
!
      LOGICAL(kind=KLOG) :: WORK_LOGICAL
!
      TYPE(ESMF_VM) :: VM_DOMAIN
!
      TYPE(ESMF_FieldBundle) :: OUTPUT_BUNDLE 
!
      TYPE :: TEMP_INT_STATE
        TYPE(WRITE_INTERNAL_STATE),POINTER :: LOC_INT_STATE
      END TYPE
!
      TYPE(TEMP_INT_STATE),DIMENSION(:),ALLOCATABLE :: WRT_INT_STATE_X
!
!----------------------------------------------------
!***  Local pointers to History or Restart variables
!***  in the Write component's internal state.
!----------------------------------------------------
!
      INTEGER(kind=KINT),POINTER :: NUM_WORDS_SEND_I2D                  &
                                   ,NUM_WORDS_SEND_R2D
!
      INTEGER(kind=KINT),DIMENSION(:),POINTER :: KOUNT_I1D              &
                                                ,KOUNT_I2D              &
                                                ,KOUNT_R1D              &
                                                ,KOUNT_R2D              &
                                                ,KOUNT_LOG
!
      INTEGER(kind=KINT),DIMENSION(:),POINTER :: LENGTH_DATA_I1D        &
                                                ,LENGTH_DATA_R1D        &
                                                ,LENGTH_DATA_R2D
!
      INTEGER(kind=KINT),DIMENSION(:),POINTER :: LENGTH_SUM_I1D         &
                                                ,LENGTH_SUM_R1D         &
                                                ,LENGTH_SUM_R2D         &
                                                ,LENGTH_SUM_LOG
!
      INTEGER(kind=KINT),DIMENSION(:),POINTER :: NUM_WORDS_RECV_I2D     &
                                                ,NUM_WORDS_RECV_R2D
!
      INTEGER(kind=KINT),DIMENSION(:),POINTER :: NCOUNT_FIELDS
!
      INTEGER(kind=KINT),DIMENSION(:),POINTER :: ALL_DATA_I1D           &
                                                ,ALL_DATA_I2D
!
      INTEGER(kind=KINT),DIMENSION(:),POINTER :: BUFF_INT
!
      INTEGER(kind=KINT),DIMENSION(:,:),POINTER :: OUTPUT_ARRAY_I2D
!
      INTEGER(kind=KINT),DIMENSION(:,:,:),POINTER :: WRITE_SUBSET_I
!
      REAL(kind=KFPT),DIMENSION(:),POINTER :: ALL_DATA_R1D              &
                                             ,ALL_DATA_R2D
!
      REAL(kind=KFPT),DIMENSION(:),POINTER :: BUFF_REAL
!
      REAL(kind=KFPT),DIMENSION(:,:),POINTER :: OUTPUT_ARRAY_R2D
!
      REAL(kind=KFPT),DIMENSION(:,:,:),POINTER :: WRITE_SUBSET_R
!
      CHARACTER(ESMF_MAXSTR),DIMENSION(:),POINTER :: FIELD_NAME
!
      CHARACTER(ESMF_MAXSTR*MAX_DATA_I1D),POINTER :: NAMES_I1D_STRING   
      CHARACTER(ESMF_MAXSTR*MAX_DATA_I2D),POINTER :: NAMES_I2D_STRING
      CHARACTER(ESMF_MAXSTR*MAX_DATA_R1D),POINTER :: NAMES_R1D_STRING
      CHARACTER(ESMF_MAXSTR*MAX_DATA_R2D),POINTER :: NAMES_R2D_STRING
      CHARACTER(ESMF_MAXSTR*MAX_DATA_LOG),POINTER :: NAMES_LOG_STRING
!
      LOGICAL,           DIMENSION(:),POINTER :: ALL_DATA_LOG

!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Get this domain's VM (ESMF Virtual Machine) and the domain's
!***  MPI communicator so we can address all the tasks as needed.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="PRELIM_INFO_FOR_HISTORY: Get the Global VM"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!
      CALL ESMF_VMGetCurrent(vm=VM_DOMAIN                               &  !<-- The VM for all tasks on this domain
                            ,rc=RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract the Communicator from the Global VM"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_VMGet(vm             =VM_DOMAIN                         &
                     ,mpiCommunicator=COMM_MY_DOMAIN                    &  !<-- The communicator for all tasks on this domain
                     ,localPet       =MYPE_DOMAIN                       &  !<-- Each task's full-domain rank
                     ,rc             =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Let the forecast tasks do all the work they need to do first.
!***  While the number of Write components is equal to WRITE_GROUPS
!***  the contents of the internal state of each of those components
!***  is identical for the forecast tasks.  The forecast tasks begin
!***  by extracting the internal state from the Write component for
!***  write group #1.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  NOTE:  All the forecast tasks are part of ALL write groups
!***         while each write task belongs to only one of the write
!***         groups.  In the work below the forecast tasks will
!***         fill all their relevant variables within the internal
!***         state of the Write component for write group #1 while
!***         providing the write tasks the quantities they need for
!***         the their variables within the internal state of the
!***         Write component associated for their particular write
!***         group.  At the end of this subroutine the forecast
!***         tasks then need to set the same internal state variables
!***         for each of the other Write components with which they
!***         are associated.
!-----------------------------------------------------------------------
!
      LAST_FCST_TASK=NUM_PES_FCST-1                                        !<-- The rank of the final forecast task
!
!-----------------------------------------------------------------------
!
      fcst_tasks_1: IF(MYPE_DOMAIN<=LAST_FCST_TASK)THEN                    !<-- Select only the forecast tasks
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="PRELIM_INFO: Fcst Tasks Get the Write Internal State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_GridCompGetInternalState(WRITE_COMPS(1)               &
                                          ,WRAP                         &
                                          ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        WRT_INT_STATE=>wrap%WRITE_INT_STATE                                !<-- Local working pointer to internal state
!
!-----------------------------------------------------------------------
!***  Extract the number of Write tasks in each group.
!-----------------------------------------------------------------------
!
        NWTPG=wrt_int_state%WRITE_TASKS_PER_GROUP
!
      ENDIF fcst_tasks_1
!
!-----------------------------------------------------------------------
!
      CALL MPI_BCAST(NWTPG                                              &  !<-- Broadcast the # of write tasks per group
                    ,1                                                  &  !<-- It is a scalar
                    ,MPI_INTEGER                                        &  !<-- It is an integer
                    ,0                                                  &  !<-- The lead forecast task broadcasts
                    ,COMM_MY_DOMAIN                                     &  !<-- The domain communicator
                    ,IERR)
!
!-----------------------------------------------------------------------
!***  The ranks of key forecast and write tasks. 
!-----------------------------------------------------------------------
!
      LEAD_WRITE_TASK(1)=LAST_FCST_TASK+1                                  !<-- The rank of the lead write task in group 1
      LAST_WRITE_TASK(1)=LAST_FCST_TASK+NWTPG                              !<-- The rank of the last write task in group 1
!
      IF(NUM_WRITE_GROUPS>=2)THEN
        DO N=2,NUM_WRITE_GROUPS
          LEAD_WRITE_TASK(N)=LEAD_WRITE_TASK(N-1)+NWTPG                    !<-- The rank of the lead write task in group N
          LAST_WRITE_TASK(N)=LEAD_WRITE_TASK(N)+NWTPG-1                    !<-- The rank of the final write task in group N
        ENDDO
      ENDIF
!
!-----------------------------------------------------------------------
!***  Now all the tasks in each Write group can extract the
!***  internal state from their respective Write component.
!-----------------------------------------------------------------------
!
      DO N=1,NUM_WRITE_GROUPS
!
        IF(MYPE_DOMAIN>=LEAD_WRITE_TASK(N).AND.                         &
           MYPE_DOMAIN<=LAST_WRITE_TASK(N))THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="PRELIM_INFO: Write Tasks Get the Write Internal State"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_GridCompGetInternalState(WRITE_COMPS(N)             &
                                            ,WRAP                       &
                                            ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          WRT_INT_STATE=>wrap%WRITE_INT_STATE                              !<-- Local working pointer to internal state
!
        ENDIF
!
      ENDDO
!
!-----------------------------------------------------------------------
!***  Now that all forecast and write tasks can see the internal state
!***  of their Write component, point the local working pointers at the
!***  appropriate components within the internal state depending upon
!***  whether we are preparing work for History or Restart output.
!***  In that way the following code can be used for both cases
!***  without a huge number of IF tests.
!-----------------------------------------------------------------------
!
      CALL POINT_LOCAL
!
!-----------------------------------------------------------------------
!
      ALLOCATE(INPES(1))
      ALLOCATE(JNPES(1))
      ALLOCATE(IHALO(1))
      ALLOCATE(JHALO(1))
      ALLOCATE(NCHAR_I2D(1))
      ALLOCATE(NCHAR_R2D(1))
!
!-----------------------------------------------------------------------
!***  Extract the appropriate output Bundle from the Write component's
!***  import state and then from it extract the full domain limits.
!***  These are needed for allocating the working arrays that will
!***  move the 2-D and 3-D Fields from the import to the export state.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
      domain_limits: IF(MYPE_DOMAIN<=LAST_FCST_TASK)THEN                   !<-- This selects only forecast tasks to do extractions
                                                                           !    since only they know what is in the import state.
!-----------------------------------------------------------------------
!
        IF(OUTPUT_FLAG=='History')THEN
          BUNDLE_NAME='History Bundle'
        ELSEIF(OUTPUT_FLAG=='Restart')THEN
          BUNDLE_NAME='Restart Bundle'
        ENDIF
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract the History Bundle from the Write Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_StateGet(state      =IMP_STATE_WRITE                  &  !<-- The Write component's import state
                          ,itemName   =BUNDLE_NAME                      &  !<-- The name of the data Bundle
                          ,fieldbundle=OUTPUT_BUNDLE                    &  !<-- The data Bundle inside the import state
                          ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Allocate the arrays that hold all start and end points in I,J
!***  for all forecast tasks.
!-----------------------------------------------------------------------
!
        ALLOCATE(LOCAL_ISTART(0:LAST_FCST_TASK),stat=RC)
        ALLOCATE(LOCAL_IEND  (0:LAST_FCST_TASK),stat=RC)
        ALLOCATE(LOCAL_JSTART(0:LAST_FCST_TASK),stat=RC)
        ALLOCATE(LOCAL_JEND  (0:LAST_FCST_TASK),stat=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract Global Parameters from Output Bundle" 
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
! 
        CALL ESMF_AttributeGet(FIELDBUNDLE=OUTPUT_BUNDLE                &  !<-- The Bundle of output data
                              ,name       ='IM'                         &  !<-- Name of the Attribute to extract
                              ,valueList  =wrt_int_state%IM             &  !<-- Extract this Attribute from History Bundle
                              ,rc         =RC)
!
        CALL ESMF_AttributeGet(FIELDBUNDLE=OUTPUT_BUNDLE                &  !<-- The Bundle of output data
                              ,name       ='JM'                         &  !<-- Name of the Attribute to extract
                              ,valueList  =wrt_int_state%JM             &  !<-- Extract this Attribute from History Bundle
                              ,rc         =RC)
!
        CALL ESMF_AttributeGet(FIELDBUNDLE=OUTPUT_BUNDLE                &  !<-- The Bundle of output data
                              ,name       ='LM'                         &  !<-- Name of the Attribute to extract
                              ,valueList  =wrt_int_state%LM             &  !<-- Extract this Attribute from History Bundle
                              ,rc         =RC)
!
        CALL ESMF_AttributeGet(FIELDBUNDLE=OUTPUT_BUNDLE                &  !<-- The Bundle of output data
                              ,name       ='GLOBAL'                     &  !<-- Name of the Attribute to extract
                              ,value      =wrt_int_state%GLOBAL         &  !<-- Extract this Attribute from History Bundle
                              ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        IF(wrt_int_state%GLOBAL) THEN                                      !<-- Increase lateral dimensions by 2 for global runs
          wrt_int_state%IM(1)=wrt_int_state%IM(1)+2
          wrt_int_state%JM(1)=wrt_int_state%JM(1)+2
        ENDIF
!
!-----------------------------------------------------------------------
!***  Now extract local subdomain limits.
!***  These will be used to allocate the working array to hold fields
!***  on each subdomain prior to quilting them together.
!***  We first need the number of forecast tasks since that
!***  determines the size of the arrays holding the local
!***  subdomain limits.
!
!***  Also extract the halo depths since they are needed for
!***  excluding halo points from the final output data.
!
!***  These values are not to be written to the output files so
!***  they were not inserted into the data Bundles inside the
!***  Write component's import state.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract Local Quilting Info from Write Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!

        CALL ESMF_AttributeGet(state=IMP_STATE_WRITE                    &  !<-- The Write component's import state
                              ,name ='INPES'                            &  !<-- Name of the Attribute to extract
                              ,value=INPES(1)                           &  !<-- Extract this Attribute from import state
                              ,rc   =RC)
!
        CALL ESMF_AttributeGet(state=IMP_STATE_WRITE                    &  !<-- The Write component's import state
                              ,name ='JNPES'                            &  !<-- Name of the Attribute to extract
                              ,value=JNPES(1)                           &  !<-- Extract this Attribute from import state
                              ,rc   =RC)
!
        wrt_int_state%INPES=INPES(1)                                       !<-- Place in internal state for later use
        wrt_int_state%JNPES=JNPES(1)                                       !<-- Place in internal state for later use
!
        CALL ESMF_AttributeGet(state=IMP_STATE_WRITE                    &  !<-- The Write component's import state
                              ,name ='IHALO'                            &  !<-- Name of the Attribute to extract
                              ,value=IHALO(1)                           &  !<-- Extract this Attribute from import state
                              ,rc   =RC)
!
        CALL ESMF_AttributeGet(state=IMP_STATE_WRITE                    &  !<-- The Write component's import state
                              ,name ='JHALO'                            &  !<-- Name of the Attribute to extract
                              ,value=JHALO(1)                           &  !<-- Extract this Attribute from import state
                              ,rc   =RC)
!
        wrt_int_state%IHALO=IHALO(1)                                       !<-- Place in internal state for later use
        wrt_int_state%JHALO=JHALO(1)                                       !<-- Place in internal state for later use
!
        CALL ESMF_AttributeGet(state    =IMP_STATE_WRITE                &  !<-- The Write component's import state
                              ,name     ='LOCAL_ISTART'                 &  !<-- Name of the Attribute to extract
                              ,valueList=LOCAL_ISTART                   &  !<-- Extract local subdomain starting I's
                              ,itemCount=NUM_PES_FCST                   &  !<-- Length of Attribute
                              ,rc       =RC)
!
        CALL ESMF_AttributeGet(state    =IMP_STATE_WRITE                &  !<-- The Write component's import state
                              ,name     ='LOCAL_IEND'                   &  !<-- Name of the Attribute to extract
                              ,valueList= LOCAL_IEND                    &  !<-- Extract local subdomain ending I's
                              ,itemCount=NUM_PES_FCST                   &  !<-- Length of Attribute
                              ,rc       =RC)
!
        CALL ESMF_AttributeGet(state    =IMP_STATE_WRITE                &  !<-- The Write component's import state
                              ,name     ='LOCAL_JSTART'                 &  !<-- Name of the Attribute to extract
                              ,valueList=LOCAL_JSTART                   &  !<-- Extract local subdomain starting J's
                              ,itemCount=NUM_PES_FCST                   &  !<-- Length of Attribute
                              ,rc       =RC)
!
        CALL ESMF_AttributeGet(state    =IMP_STATE_WRITE                &  !<-- The Write component's import state
                              ,name     ='LOCAL_JEND'                   &  !<-- Name of the Attribute to extract
                              ,valueList=LOCAL_JEND                     &  !<-- Extract local subdomain ending J's
                              ,itemCount=NUM_PES_FCST                   &  !<-- Length of Attribute
                              ,rc       =RC)

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        DO N=0,LAST_FCST_TASK
          wrt_int_state%LOCAL_ISTART(N)=LOCAL_ISTART(N)
          wrt_int_state%LOCAL_IEND  (N)=LOCAL_IEND(N)
          wrt_int_state%LOCAL_JSTART(N)=LOCAL_JSTART(N)
          wrt_int_state%LOCAL_JEND  (N)=LOCAL_JEND(N)
        ENDDO
!
        ITS=LOCAL_ISTART(MYPE_DOMAIN)
        ITE=LOCAL_IEND(MYPE_DOMAIN)
        JTS=LOCAL_JSTART(MYPE_DOMAIN)
        JTE=LOCAL_JEND(MYPE_DOMAIN)
!
        DEALLOCATE(LOCAL_ISTART)
        DEALLOCATE(LOCAL_IEND  )
        DEALLOCATE(LOCAL_JSTART)
        DEALLOCATE(LOCAL_JEND  )
!
!-----------------------------------------------------------------------
!
      ENDIF domain_limits
!
!-----------------------------------------------------------------------
!***  Forecast task 0 sends the domain size information
!***  to the first Write task in each Write group because 
!***  the Write tasks need to know this to assemble the
!***  final gridded data.
!-----------------------------------------------------------------------
!
      n_groups_1: DO N=1,NUM_WRITE_GROUPS
!
!-----------------------------------------------------------------------
!
        N1=LEAD_WRITE_TASK(N)
        N2=LAST_WRITE_TASK(N)
!
!-----------------------------------------------------------------------
!                            -- IM --
!-----------------------------------------------------------------------
!
        IF(MYPE_DOMAIN==0)THEN                                             !<-- Forecast task 0 sends
          DO NN=N1,N2     
            CALL MPI_SEND(wrt_int_state%IM                              &  !<-- Send this data
                         ,1                                             &  !<-- Number of words sent
                         ,MPI_INTEGER                                   &  !<-- Datatype
                         ,NN                                            &  !<-- Send to each of the write tasks (domain IDs)
                         ,0                                             &  !<-- An MPI tag
                         ,COMM_MY_DOMAIN                                &  !<-- MPI domain communicator
                         ,IERR)
!
            IF(IERR/=0)WRITE(0,*)' Failed to send IM from fcst task 0 to write tasks'
!
          ENDDO
        ENDIF 
!
        IF(MYPE_DOMAIN>=N1.AND.MYPE_DOMAIN<=N2)THEN                        !<-- Write tasks in this group receive
          CALL MPI_RECV(wrt_int_state%IM                                &  !<-- Recv this data
                       ,1                                               &  !<-- Words received
                       ,MPI_INTEGER                                     &  !<-- Datatype
                       ,0                                               &  !<-- Recv from fcst task 0
                       ,0                                               &  !<-- An MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- MPI domain communicator
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR)
!
          IF(IERR/=0)WRITE(0,*)' Write tasks failed to receive IM from fcst task0'
!
        ENDIF
!
!-----------------------------------------------------------------------
!                          -- JM --
!-----------------------------------------------------------------------
!
        IF(MYPE_DOMAIN==0)THEN                                             !<-- Forecast task 0 sends
          DO NN=N1,N2                                               
            CALL MPI_SEND(wrt_int_state%JM                              &  !<-- Send this data
                         ,1                                             &  !<-- Number of words sent
                         ,MPI_INTEGER                                   &  !<-- Datatype
                         ,NN                                            &  !<-- Send to each of the write tasks (local IDs)
                         ,0                                             &  !<-- An MPI tag
                         ,COMM_MY_DOMAIN                                &  !<-- MPI domain communicator
                         ,IERR)
!
            IF(IERR/=0)WRITE(0,*)' Failed to send JM from fcst task0 to write tasks'
!
          ENDDO
        ENDIF 
!
        IF(MYPE_DOMAIN>=N1.AND.MYPE_DOMAIN<=N2)THEN                        !<-- Write tasks in this group receive
          CALL MPI_RECV(wrt_int_state%JM                                &  !<-- Recv this data
                       ,1                                               &  !<-- Words received
                       ,MPI_INTEGER                                     &  !<-- Datatype
                       ,0                                               &  !<-- Recv from fcst task 0
                       ,0                                               &  !<-- An MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- MPI domain communicator
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR)
!
          IF(IERR/=0)WRITE(0,*)' Write tasks failed to receive JM from fcst task0'
!
        ENDIF
!
!-----------------------------------------------------------------------
!                          -- LM --
!-----------------------------------------------------------------------
!
        IF(MYPE_DOMAIN==0)THEN                                             !<-- Forecast task 0 sends
          DO NN=N1,N2                                               
            CALL MPI_SEND(wrt_int_state%LM                              &  !<-- Send this data
                         ,1                                             &  !<-- Number of words sent
                         ,MPI_INTEGER                                   &  !<-- Datatype
                         ,NN                                            &  !<-- Send to each of the write tasks (local IDs)
                         ,0                                             &  !<-- An MPI tag
                         ,COMM_MY_DOMAIN                                &  !<-- MPI domain communicator
                         ,IERR)
!
            IF(IERR/=0)WRITE(0,*)' Failed to send LM from fcst task0 to write tasks'
!
          ENDDO
        ENDIF 
!
        IF(MYPE_DOMAIN>=N1.AND.MYPE_DOMAIN<=N2)THEN                        !<-- Write tasks in this group receive
          CALL MPI_RECV(wrt_int_state%LM                                &  !<-- Recv this data
                       ,1                                               &  !<-- Words received
                       ,MPI_INTEGER                                     &  !<-- Datatype
                       ,0                                               &  !<-- Recv from fcst task 0
                       ,0                                               &  !<-- An MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- MPI domain communicator
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR)
!
          IF(IERR/=0)WRITE(0,*)' Write tasks failed to receive LM from fcst task0'
!
        ENDIF
!
!-----------------------------------------------------------------------
!
      ENDDO n_groups_1
!
!-----------------------------------------------------------------------
!
      IM=wrt_int_state%IM(1)
      JM=wrt_int_state%JM(1)
!
!-----------------------------------------------------------------------
!***  The number of Attributes (for scalars and 1D arrays) and
!***  Fields (for gridded 2D arrays) in the Write component's
!***  import state are not known a priori.  In order to transfer
!***  them to the Write tasks, extract the number of each of
!***  them along with their names.  The scalars can be lumped in
!***  with the 1D arrays at this point.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  All integer quantities (as 1D arrays) and 1D and 2D real
!***  quantities will be strung together in single arrays of 
!***  each particular type.  Arrays that will hold the length of  
!***  each of the quantities in these 'strings' were allocated
!***  in WRT_INIT.
!-----------------------------------------------------------------------
!
      KOUNT_I1D_X=0
      KOUNT_R1D_X=0
      KOUNT_LOG_X=0
!
      LENGTH_SUM_I1D_X=0
      LENGTH_SUM_R1D_X=0
      LENGTH_SUM_LOG_X=0
!
      NCHAR_I1D=0
      NCHAR_R1D=0
      NCHAR_LOG=0
!
!-----------------------------------------------------------------------
!
      fcst_tasks_2: IF(MYPE_DOMAIN<=LAST_FCST_TASK)THEN                    !<-- Only forecast tasks will extract output information
                                                                           !    from the import state because only they participated
                                                                           !    in filling the import state in the Solver component.
!
!-----------------------------------------------------------------------
!***  First find the number of Attributes in the output data Bundle
!***  in the import state and then find their names, lengths, and
!***  datatypes.
!***  Extract the integer and real data and pack it into integer
!***  and real buffers.  Later the buffers will be sent from the
!***  Forecast tasks (the only ones who can see the original
!***  data) to the write tasks.
!
!***  The fact that the Attribute output data is being collected
!***  here in a block that executes only once per Write group
!***  implies the assumption that only the 2D/3D data
!***  associated with the forecast grid can change with time.
!***  If any scalar/1D Attribute data change with time then
!***  this must be moved out of this routine and into the Run step.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Get Attribute Count from "//BUNDLE_NAME
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(FIELDBUNDLE=OUTPUT_BUNDLE                &  !<-- The write component's history data Bundle
                              ,count      =NUM_ATTRIB                   &  !<-- # of Attributes in the history data Bundle
                              ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
        attribute_loop: DO N=1,NUM_ATTRIB                                  !<-- Loop through all the Attributes
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Get Attribute Names, Datatypes, Lengths"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeGet(FIELDBUNDLE   =OUTPUT_BUNDLE           &  !<-- The write component's history data Bundle
                                ,attributeIndex=N                       &  !<-- Index of each Attribute
                                ,name          =ATTRIB_NAME             &  !<-- Each Attribute's name
                                ,typekind      =DATATYPE                &  !<-- Each Attribute's ESMF Datatype
                                ,itemCount     =LENGTH                  &  !<-- Each Attribute's length
                                ,rc            =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!                 -- Scalar and 1-D Integer Output Data --
!-----------------------------------------------------------------------
!
          IF(DATATYPE==ESMF_TYPEKIND_I4)THEN                               !<-- Extract integer data with rank <2
!
            ALLOCATE(WORK_ARRAY_I1D(LENGTH),stat=RC)                       !<-- This length is from the preceding call 
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Get Scalar/1-D Integer Data from "//BUNDLE_NAME
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_AttributeGet(FIELDBUNDLE=OUTPUT_BUNDLE            &  !<-- The write component's history data Bundle
                                  ,name       =ATTRIB_NAME              &  !<-- Name of the Attribute to extract
                                  ,itemCount  =LENGTH                   &  !<-- Length of Attribute
                                  ,valueList  =WORK_ARRAY_I1D           &  !<-- Place the Attribute here
                                  ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            KOUNT_I1D_X=KOUNT_I1D_X+1                                      !<-- Count # of integer Attributes
!
            NPOSN_END=KOUNT_I1D_X*ESMF_MAXSTR
            NPOSN_START=NPOSN_END-ESMF_MAXSTR+1     
            NCHAR_I1D=NCHAR_I1D+ESMF_MAXSTR                                !<-- Save #of characters in all scalar/1D integer names.
                                                                           !    Note that each name is being given
                                                                           !    EMSF_MAXSTR total spaces.
!
            DO L=1,LENGTH
              ALL_DATA_I1D(LENGTH_SUM_I1D_X+L)=WORK_ARRAY_I1D(L)           !<-- String together the integer data
            ENDDO
!
            LENGTH_SUM_I1D_X=LENGTH_SUM_I1D_X+LENGTH                       !<-- Total word sum of scalar/1D integer data
!
            NAMES_I1D_STRING(NPOSN_START:NPOSN_END)=ATTRIB_NAME            !<-- Save the 1D integer names
            LENGTH_DATA_I1D(KOUNT_I1D_X)=LENGTH                            !<-- Store length of each individual integer variable
!
            DEALLOCATE(WORK_ARRAY_I1D)
!
!-----------------------------------------------------------------------
!                  -- Scalar and 1-D Real Output Data --
!-----------------------------------------------------------------------
!
          ELSEIF(DATATYPE==ESMF_TYPEKIND_R4)THEN                           ! <-- Extract real data with rank <2 as Attributes
!
            ALLOCATE(WORK_ARRAY_R1D(LENGTH),stat=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Get Scalar/1-D Real Data from "//BUNDLE_NAME
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_AttributeGet(FIELDBUNDLE=OUTPUT_BUNDLE            &  !<-- The write component's history data Bundle
                                  ,name       =ATTRIB_NAME              &  !<-- Name of the Attribute to extract
                                  ,itemCount  =LENGTH                   &  !<-- Length of Attribute
                                  ,valueList  =WORK_ARRAY_R1D           &  !<-- Place the Attribute here
                                  ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            KOUNT_R1D_X=KOUNT_R1D_X+1                                      !<-- Count # of real Attributes
!
            NPOSN_END=KOUNT_R1D_X*ESMF_MAXSTR
            NPOSN_START=NPOSN_END-ESMF_MAXSTR+1     
            NCHAR_R1D=NCHAR_R1D+ESMF_MAXSTR                                !<-- Save #of characters in all scalar/1D real names
                                                                           !    Note that each name is being given
                                                                           !    EMSF_MAXSTR total spaces
!
            DO L=1,LENGTH
              ALL_DATA_R1D(LENGTH_SUM_R1D_X+L)=WORK_ARRAY_R1D(L)           !<-- String together the real data
            ENDDO
!
            LENGTH_SUM_R1D_X=LENGTH_SUM_R1D_X+LENGTH                       !<-- Total word sum of scalar/1D real data
!
            NAMES_R1D_STRING(NPOSN_START:NPOSN_END)=ATTRIB_NAME            !<-- Save the scalar/1D real names
            LENGTH_DATA_R1D(KOUNT_R1D_X)=LENGTH                            !<-- Store length of each individual real variable
!
            DEALLOCATE(WORK_ARRAY_R1D)
!
!-----------------------------------------------------------------------
!                          -- Logical Data --                       
!-----------------------------------------------------------------------
!
!ratko    ELSEIF(DATATYPE==ESMF_DATA_LOGICAL)THEN                          ! <-- Extract logical data
! --- nothing else than I4 and R4; should FIX later
          ELSE
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Get Logical Data from "//BUNDLE_NAME
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_AttributeGet(FIELDBUNDLE=OUTPUT_BUNDLE            &  !<-- The write component's history data Bundle
                                  ,name  =ATTRIB_NAME                   &  !<-- Name of the Attribute to extract
                                  ,value =WORK_LOGICAL                  &  !<-- Place the Attribute here
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            KOUNT_LOG_X=KOUNT_LOG_X+1                                      !<-- Count # of logical Attributes
!
            NPOSN_END=KOUNT_LOG_X*ESMF_MAXSTR
            NPOSN_START=NPOSN_END-ESMF_MAXSTR+1     
            NCHAR_LOG=NCHAR_LOG+ESMF_MAXSTR                                !<-- Save #of characters in all logical names
                                                                           !    Note that each name is being given
                                                                           !    EMSF_MAXSTR total spaces
!
            LENGTH_SUM_LOG_X=LENGTH_SUM_LOG_X+1                            !<-- Total length of all logical data variables
!
            NAMES_LOG_STRING(NPOSN_START:NPOSN_END)=ATTRIB_NAME            !<-- Save the logical names
!
            ALL_DATA_LOG(KOUNT_LOG_X)=WORK_LOGICAL                         !<-- String together the logical data
!
          ENDIF
!
!-----------------------------------------------------------------------
!
        ENDDO attribute_loop
!
!-----------------------------------------------------------------------
!***  Insert number and lengths of scalar/1D integer and real quantities
!***  and logicals into the Write component's internal state.
!-----------------------------------------------------------------------
!
        KOUNT_I1D(1)=KOUNT_I1D_X
        KOUNT_R1D(1)=KOUNT_R1D_X
        KOUNT_LOG(1)=KOUNT_LOG_X
!
        LENGTH_SUM_I1D(1)=LENGTH_SUM_I1D_X
        LENGTH_SUM_R1D(1)=LENGTH_SUM_R1D_X
        LENGTH_SUM_LOG(1)=LENGTH_SUM_LOG_X
!
!-----------------------------------------------------------------------
!***  Now extract the number of ESMF Fields in the output data Bundle
!***  along with their names.  Save the Field information of the 2-D 
!***  data since it will be needed for extraction from the import
!***  state and the transfer to the write tasks.
!***  Then extract the names of all the Fields in the bundle.  
!***  Also, the number of Field names returned should equal 
!***  the Field count in the preceding call.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Get Field Count from "//BUNDLE_NAME
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_FieldBundleGet(           OUTPUT_BUNDLE               &  !<-- The write component's data Bundle
                                ,fieldCount=NCOUNT_FIELDS(1)            &  !<-- Get total # of Fields in the data Bundle
                                ,rc        =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract Field Names from "//BUNDLE_NAME
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_FieldBundleGet(              OUTPUT_BUNDLE            &  !<-- The write component's data Bundle
                                ,fieldNameList=FIELD_NAME               &  !<-- Array of ESMF Field names in the Bundle
                                ,fieldCount   =NUM_FIELD_NAMES          &  !<-- Number of Field names in the Bundle
                                ,rc           =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
        IF(NUM_FIELD_NAMES/=NCOUNT_FIELDS(1))THEN
          WRITE(0,*)' WARNING: Number of Fields in '//BUNDLE_NAME       &
                   ,' output does not equal the number of Field names'
          WRITE(0,*)' They are ',NUM_FIELD_NAMES,' and '            &
                   ,NCOUNT_FIELDS(1),', respectively'
        ENDIF
!
!-----------------------------------------------------------------------
!***  Do a preliminary extraction of the Fields themselves in order to
!***  count the number of real and integer 2D arrays.  
!-----------------------------------------------------------------------
!
        KOUNT_R2D_X=0
        KOUNT_I2D_X=0
        NCHAR_I2D(1)=0
        NCHAR_R2D(1)=0
!
        DO N=1,NCOUNT_FIELDS(1)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Extract Fields from Output Bundle for Counting"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_FieldBundleGet(FIELDBUNDLE=OUTPUT_BUNDLE               &  !<-- The write component's data Bundle
                                  ,fieldName  =FIELD_NAME(N)               &  !<-- The ESMF Field's name
                                  ,field      =FIELD_WORK1                 &  !<-- The ESMF Field taken from the Bundle
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Get Datatype of Fields for Counting Real/Integer"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_FieldGet(field   =FIELD_WORK1                       &  !<-- The ESMF 2D Field
                            ,typekind=DATATYPE                          &  !<-- The Field's ESMF Datatype
                            ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          IF(DATATYPE==ESMF_TYPEKIND_I4)THEN
            KOUNT_I2D_X=KOUNT_I2D_X+1                                      !<-- Add up the total number of integer 2D Fields
!
            IF(KOUNT_I2D_X>MAX_DATA_I2D)THEN
              WRITE(0,*)' FATAL: YOU HAVE EXCEEDED MAX NUMBER OF INTEGER 2D FIELDS FOR OUTPUT'
              WRITE(0,*)' YOU MUST INCREASE VALUE OF MAX_DATA_I2D WHICH NOW EQUALS ',MAX_DATA_I2D
            ENDIF
!
            NPOSN_END  =KOUNT_I2D_X*ESMF_MAXSTR
            NPOSN_START=NPOSN_END-ESMF_MAXSTR+1
            NAMES_I2D_STRING(NPOSN_START:NPOSN_END)=FIELD_NAME(N)          !<-- Save the 2D integer Field names 
                                                                           !    in one long string.
            NCHAR_I2D(1)=NCHAR_I2D(1)+ESMF_MAXSTR
!
          ELSEIF(DATATYPE==ESMF_TYPEKIND_R4)THEN
            KOUNT_R2D_X=KOUNT_R2D_X+1                                      !<-- Add up the total number of real 2D Fields
!
            IF(KOUNT_R2D_X>MAX_DATA_R2D)THEN
              WRITE(0,*)' FATAL: YOU HAVE EXCEEDED MAX NUMBER OF REAL 2D FIELDS FOR OUTPUT'
              WRITE(0,*)' YOU MUST INCREASE VALUE OF MAX_DATA_R2D WHICH NOW EQUALS ',MAX_DATA_R2D
              CALL ESMF_Finalize(endflag=ESMF_END_ABORT             &
                                ,rc             =RC )
            ENDIF
!
            NPOSN_END  =KOUNT_R2D_X*ESMF_MAXSTR
            NPOSN_START=NPOSN_END-ESMF_MAXSTR+1
            NAMES_R2D_STRING(NPOSN_START:NPOSN_END)=FIELD_NAME(N)          !<-- Save the 2D real Field names 
                                                                           !    in one long string.
            NCHAR_R2D(1)=NCHAR_R2D(1)+ESMF_MAXSTR
!
          ENDIF
!
        ENDDO
!
        KOUNT_R2D(1)=KOUNT_R2D_X
        KOUNT_I2D(1)=KOUNT_I2D_X
!
!-----------------------------------------------------------------------
!***  Compute the total number of words for all 2D and 3D real data
!***  and allocate a datastring to that length.  It will transfer
!***  the 2D/3D real data from forecast to Write tasks.
!-----------------------------------------------------------------------
!
        NUM_WORDS_TOT=(ITE-ITS+1+2*IHALO(1))*(JTE-JTS+1+2*JHALO(1))     &
                      *KOUNT_R2D_X
!
        IF(NUM_WORDS_TOT>2147483647)THEN
          WRITE(0,*)' You have TOO MANY words in your datastring.'
          WRITE(0,*)' You must increase the number of tasks.'
          CALL ESMF_Finalize(endflag=ESMF_END_ABORT                 &
                            ,rc             =RC  )
        ELSE
          NUM_WORDS_SEND_R2D=NUM_WORDS_TOT
          NUM_WORDS=NUM_WORDS_SEND_R2D
!
          IF(OUTPUT_FLAG=='History')THEN
            IF(.NOT.ASSOCIATED(wrt_int_state%ALL_DATA_R2D))THEN
              ALLOCATE(wrt_int_state%ALL_DATA_R2D(NUM_WORDS),stat=ISTAT)
              IF(ISTAT/=0)THEN
                WRITE(0,*)' Fcst task FAILED to allocate ALL_DATA_R2D'
                CALL ESMF_Finalize(endflag=ESMF_END_ABORT           &
                                  ,rc             =RC  )
              ENDIF
              ALL_DATA_R2D=>wrt_int_state%ALL_DATA_R2D
            ENDIF
!
          ELSEIF(OUTPUT_FLAG=='Restart')THEN
            IF(.NOT.ASSOCIATED(wrt_int_state%RST_ALL_DATA_R2D))THEN
              ALLOCATE(wrt_int_state%RST_ALL_DATA_R2D(NUM_WORDS),stat=ISTAT)
              IF(ISTAT/=0)THEN
                WRITE(0,*)' Fcst task FAILED to allocate RST_ALL_DATA_R2D'
                CALL ESMF_Finalize(endflag=ESMF_END_ABORT           &
                                  ,rc             =RC  )
              ENDIF
              ALL_DATA_R2D=>wrt_int_state%RST_ALL_DATA_R2D
            ENDIF
!
          ENDIF
!
        ENDIF
!
!-----------------------------------------------------------------------
!***  Likewise for 2D integer data.
!-----------------------------------------------------------------------
!
        NUM_WORDS_TOT=(ITE-ITS+1+2*IHALO(1))*(JTE-JTS+1+2*JHALO(1))     &
                      *KOUNT_I2D_X
!
        IF(NUM_WORDS_TOT>2147483647)THEN
          WRITE(0,*)' You have TOO MANY words in your datastring.'
          WRITE(0,*)' You must increase the number of tasks.'
          CALL ESMF_Finalize(endflag=ESMF_END_ABORT                 &
                            ,rc             =RC  )
        ELSE
          NUM_WORDS_SEND_I2D=NUM_WORDS_TOT
          NUM_WORDS=NUM_WORDS_SEND_I2D
!
          IF(OUTPUT_FLAG=='History')THEN
            IF(.NOT.ASSOCIATED(wrt_int_state%ALL_DATA_I2D))THEN
              ALLOCATE(wrt_int_state%ALL_DATA_I2D(NUM_WORDS),stat=ISTAT)
              IF(ISTAT/=0)THEN
                WRITE(0,*)' Fcst task FAILED to allocate ALL_DATA_I2D'
                CALL ESMF_Finalize(endflag=ESMF_END_ABORT           &
                                  ,rc             =RC  )
              ENDIF
              ALL_DATA_I2D=>wrt_int_state%ALL_DATA_I2D
            ENDIF
!
          ELSEIF(OUTPUT_FLAG=='Restart')THEN
            IF(.NOT.ASSOCIATED(wrt_int_state%RST_ALL_DATA_I2D))THEN
              ALLOCATE(wrt_int_state%RST_ALL_DATA_I2D(NUM_WORDS),stat=ISTAT)
              IF(ISTAT/=0)THEN
                WRITE(0,*)' Fcst task FAILED to allocate RST_ALL_DATA_I2D'
                CALL ESMF_Finalize(endflag=ESMF_END_ABORT           &
                                  ,rc             =RC  )
              ENDIF
              ALL_DATA_I2D=>wrt_int_state%RST_ALL_DATA_I2D
            ENDIF
!
          ENDIF
!
        ENDIF
!
!-----------------------------------------------------------------------
!
      ENDIF fcst_tasks_2
!
!-----------------------------------------------------------------------
!***  If there are no quantities specified for history output,
!***  forecast task 0 will inform the write tasks and then
!***  everyone will return.
!-----------------------------------------------------------------------
!
      NO_FIELDS=.FALSE.
!
!-----------------------------------------------------------------------
      n_groups_2: DO N=1,NUM_WRITE_GROUPS
!-----------------------------------------------------------------------
!
        N1=LEAD_WRITE_TASK(N)                                              !<-- The lead write task in group N
        N2=LAST_WRITE_TASK(N)                                              !<-- The last write task in group N
!
        IF(MYPE_DOMAIN==0)THEN
          IF(NCOUNT_FIELDS(1)==0)NO_FIELDS=.TRUE.                          !<-- Reset flag saying there are no history quantities
!
          DO NN=N1,N2                                                      !<-- Loop through all the write tasks in the write group
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Fcst Task0 Informs All That There Are No Fields"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL MPI_SEND(NO_FIELDS                                     &  !<-- Send flag that there is no history data
                         ,1                                             &  !<-- The flag is a scalar
                         ,MPI_LOGICAL                                   &  !<-- The flag is logical
                         ,NN                                            &  !<-- Receiving task in write group N
                         ,0                                             &  !<-- MPI tag
                         ,COMM_MY_DOMAIN                                &  !<-- MPI domain communicator
                         ,IERR )          
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          ENDDO
!
        ENDIF
!
        IF(MYPE_DOMAIN>=N1.AND.MYPE_DOMAIN<=N2)THEN                        !<-- All write tasks in this group receive
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Write Tasks Told By Fcst Task0 There Are No Fields"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL MPI_RECV(NO_FIELDS                                       &  !<-- Recv flag that there is no history data
                       ,1                                               &  !<-- The flag is a scalar
                       ,MPI_LOGICAL                                     &  !<-- The flag is logicallar
                       ,0                                               &  !<-- Sending task (forecast task 0)
                       ,0                                               &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- MPI domain communicator
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR )         
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        ENDIF
!
      ENDDO n_groups_2
!
!-----------------------------------------------------------------------
!
      IF(NO_FIELDS)THEN
        IF(MYPE_DOMAIN==0)THEN
          WRITE(6,*)'WARNING: No Import ESMF quantities for the Write Component'
          WRITE(0,*)'WARNING: No Import ESMF quantities for the Write Component'
        ENDIF
!
        RETURN                                                             !<-- All tasks return if there is no history output
!
      ENDIF
!
!-----------------------------------------------------------------------
!***  Forecast task 0 sends all the write tasks the number of
!***  real and integer 2D gridded quantities plus all of the
!***  local horizontal domain limits in preparation for the
!***  Write tasks' receiving and assembling the local history
!***  data they receive from the Forecast tasks.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
      n_groups_3: DO N=1,NUM_WRITE_GROUPS
!-----------------------------------------------------------------------
!
        N1=LEAD_WRITE_TASK(N)
        N2=LAST_WRITE_TASK(N)
!
!-----------------------------------------------------------------------
!
        task_0: IF(MYPE_DOMAIN==0)THEN                                     !<-- Forecast task 0 sends
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Fcst Task0 Sends Write Tasks Info for Quilting"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          DO NN=N1,N2                                                      !<-- Loop through all the write tasks in the write group
!
            CALL MPI_SEND(INPES                                         &  !<-- Send # of fcst tasks in I
                         ,1                                             &  !<-- Number is one word
                         ,MPI_INTEGER                                   &  !<-- Data is integer
                         ,NN                                            &  !<-- Receiving task in active write group
                         ,0                                             &  !<-- MPI tag
                         ,COMM_MY_DOMAIN                                &  !<-- MPI domain communicator
                         ,IERR )
!
            CALL MPI_SEND(JNPES                                         &  !<-- Send # of fcst tasks in J
                         ,1                                             &  !<-- Number is one word
                         ,MPI_INTEGER                                   &  !<-- Data is integer
                         ,NN                                            &  !<-- Receiving task in active write group
                         ,0                                             &  !<-- MPI tag
                         ,COMM_MY_DOMAIN                                &  !<-- MPI domain communicator
                         ,IERR )
!
            CALL MPI_SEND(IHALO                                         &  !<-- Send # of points in the east-west halo
                         ,1                                             &  !<-- Number is one word
                         ,MPI_INTEGER                                   &  !<-- Data is integer
                         ,NN                                            &  !<-- Receiving task in active write group
                         ,0                                             &  !<-- MPI tag
                         ,COMM_MY_DOMAIN                                &  !<-- MPI domain communicator
                         ,IERR )
!
            CALL MPI_SEND(JHALO                                         &  !<-- Send # of points in the south-north halo
                         ,1                                             &  !<-- Number is one word
                         ,MPI_INTEGER                                   &  !<-- Data is integer
                         ,NN                                            &  !<-- Receiving task in active write group
                         ,0                                             &  !<-- MPI tag
                         ,COMM_MY_DOMAIN                                &  !<-- MPI domain communicator
                         ,IERR )
!
            CALL MPI_SEND(NCOUNT_FIELDS                                 &  !<-- Send # of Fields in history data
                         ,1                                             &  !<-- Number is one word
                         ,MPI_INTEGER                                   &  !<-- Data is integer
                         ,NN                                            &  !<-- Receiving task in active write group
                         ,0                                             &  !<-- MPI tag
                         ,COMM_MY_DOMAIN                                &  !<-- MPI domain communicator
                         ,IERR )
!
            CALL MPI_SEND(KOUNT_R2D                                     &  !<-- Send # of real Fields in history data
                         ,1                                             &  !<-- Number is one word
                         ,MPI_INTEGER                                   &  !<-- Data is integer
                         ,NN                                            &  !<-- Receiving task in active write group
                         ,0                                             &  !<-- MPI tag
                         ,COMM_MY_DOMAIN                                &  !<-- MPI domain communicator
                         ,IERR )
!
            CALL MPI_SEND(KOUNT_I2D                                     &  !<-- Send # of integer Fields in history data
                         ,1                                             &  !<-- Number is one word
                         ,MPI_INTEGER                                   &  !<-- Data is integer
                         ,NN                                            &  !<-- Receiving task in active write group
                         ,0                                             &  !<-- MPI tag
                         ,COMM_MY_DOMAIN                                &  !<-- MPI domain communicator
                         ,IERR )
!
            CALL MPI_SEND(wrt_int_state%LOCAL_ISTART                    &  !<-- Send starting I for all fcst task subdomains
                         ,NUM_PES_FCST                                  &  !<-- Number of words
                         ,MPI_INTEGER                                   &  !<-- Data is integer
                         ,NN                                            &  !<-- Receiving task in active write group
                         ,0                                             &  !<-- MPI tag
                         ,COMM_MY_DOMAIN                                &  !<-- MPI domain communicator
                         ,IERR )
!
            CALL MPI_SEND(wrt_int_state%LOCAL_IEND                      &  !<-- Send ending I for all fcst task subdomains
                         ,NUM_PES_FCST                                  &  !<-- Number of words
                         ,MPI_INTEGER                                   &  !<-- Data is integer
                         ,NN                                            &  !<-- Receiving task in active write group
                         ,0                                             &  !<-- MPI tag
                         ,COMM_MY_DOMAIN                                &  !<-- MPI domain communicator
                         ,IERR )
!
            CALL MPI_SEND(wrt_int_state%LOCAL_JSTART                    &  !<-- Send starting J for all fcst task subdomains
                         ,NUM_PES_FCST                                  &  !<-- Number of words
                         ,MPI_INTEGER                                   &  !<-- Data is integer
                         ,NN                                            &  !<-- Receiving task in active write group
                         ,0                                             &  !<-- MPI tag
                         ,COMM_MY_DOMAIN                                &  !<-- MPI domain communicator
                         ,IERR )
!
            CALL MPI_SEND(wrt_int_state%LOCAL_JEND                      &  !<-- Send starting J for all fcst task subdomains
                         ,NUM_PES_FCST                                  &  !<-- Number of words
                         ,MPI_INTEGER                                   &  !<-- Data is integer
                         ,NN                                            &  !<-- Receiving task in active write group
                         ,0                                             &  !<-- MPI tag
                         ,COMM_MY_DOMAIN                                &  !<-- MPI domain communicator
                         ,IERR )
!
          ENDDO
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Fcst Task0 Sends Lead Write Task # of BC Words"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL MPI_SEND(wrt_int_state%NUM_WORDS_SEND_BC                 &  !<-- Send # of words in BC data
                       ,1                                               &  !<-- Send one word
                       ,MPI_INTEGER                                     &  !<-- Data is integer
                       ,N1                                              &  !<-- Receiving task in active write group
                       ,0                                               &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- MPI domain communicator
                       ,IERR )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!
        ENDIF task_0
!
!-----------------------------------------------------------------------
!
!
        IF(MYPE_DOMAIN>=N1.AND.MYPE_DOMAIN<=N2)THEN                        !<-- All write tasks in this group receive
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Write Tasks Recv Quilting Info From Fcst Task0"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL MPI_RECV(INPES                                           &  !<-- Recv # of fcst tasks in I
                       ,1                                               &  !<-- Data is a scalar
                       ,MPI_INTEGER                                     &  !<-- Data is integer
                       ,0                                               &  !<-- Sending task (forecast task 0)
                       ,0                                               &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- MPI domain communicator
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR )
!
          CALL MPI_RECV(JNPES                                           &  !<-- Recv # of fcst tasks in J
                       ,1                                               &  !<-- Data is a scalar
                       ,MPI_INTEGER                                     &  !<-- Data is integer
                       ,0                                               &  !<-- Sending task (forecast task 0)
                       ,0                                               &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- MPI domain communicator
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR )
!
          CALL MPI_RECV(IHALO                                           &  !<-- Recv width of halo in I
                       ,1                                               &  !<-- Data is a scalar
                       ,MPI_INTEGER                                     &  !<-- Data is integer
                       ,0                                               &  !<-- Sending task (forecast task 0)
                       ,0                                               &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- MPI domain communicator
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR )
!
          CALL MPI_RECV(JHALO                                           &  !<-- Recv width of halo in J
                       ,1                                               &  !<-- Data is a scalar
                       ,MPI_INTEGER                                     &  !<-- Data is integer
                       ,0                                               &  !<-- Sending task (forecast task 0)
                       ,0                                               &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- MPI domain communicator
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR )
!
          wrt_int_state%IHALO=IHALO(1)
          wrt_int_state%JHALO=JHALO(1)
!
          CALL MPI_RECV(NCOUNT_FIELDS                                   &  !<-- Recv # of Fields in history data
                       ,1                                               &  !<-- Data is a scalar
                       ,MPI_INTEGER                                     &  !<-- Data is integer
                       ,0                                               &  !<-- Sending task (forecast task 0)
                       ,0                                               &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- MPI domain communicator
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR )
!
          CALL MPI_RECV(KOUNT_R2D                                       &  !<-- Recv # of real Fields in history data
                       ,1                                               &  !<-- Data is a scalar
                       ,MPI_INTEGER                                     &  !<-- Data is integer
                       ,0                                               &  !<-- Sending task (forecast task 0)
                       ,0                                               &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- MPI domain communicator
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR )
!
          CALL MPI_RECV(KOUNT_I2D                                       &  !<-- Recv # of integer Fields in history data
                       ,1                                               &  !<-- Data is a scalar
                       ,MPI_INTEGER                                     &  !<-- Data is integer
                       ,0                                               &  !<-- Sending task (forecast task 0)
                       ,0                                               &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- MPI domain communicator
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR )
!
          CALL MPI_RECV(wrt_int_state%LOCAL_ISTART                      &  !<-- Recv starting I of each fcst task subdomain
                       ,NUM_PES_FCST                                    &  !<-- # of words in data
                       ,MPI_INTEGER                                     &  !<-- Data is integer
                       ,0                                               &  !<-- Sending task (forecast task 0)
                       ,0                                               &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- MPI domain communicator
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR )
!
          CALL MPI_RECV(wrt_int_state%LOCAL_IEND                        &  !<-- Recv ending I of each fcst task subdomain
                       ,NUM_PES_FCST                                    &  !<-- # of words in data
                       ,MPI_INTEGER                                     &  !<-- Data is integer
                       ,0                                               &  !<-- Sending task (forecast task 0)
                       ,0                                               &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- MPI domain communicator
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR )
!
          CALL MPI_RECV(wrt_int_state%LOCAL_JSTART                      &  !<-- Recv starting J of each fcst task subdomain
                       ,NUM_PES_FCST                                    &  !<-- # of words in data
                       ,MPI_INTEGER                                     &  !<-- Data is integer
                       ,0                                               &  !<-- Sending task (forecast task 0)
                       ,0                                               &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- MPI domain communicator
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR )
!
          CALL MPI_RECV(wrt_int_state%LOCAL_JEND                        &  !<-- Recv ending J of each fcst task subdomain
                       ,NUM_PES_FCST                                    &  !<-- # of words in data
                       ,MPI_INTEGER                                     &  !<-- Data is integer
                       ,0                                               &  !<-- Sending task (forecast task 0)
                       ,0                                               &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- MPI domain communicator
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          IF(MYPE_DOMAIN==N1)THEN                                          !<-- Lead write task in each group must recv BC word count
!
            IF(ALLOCATED(wrt_int_state%NUM_WORDS_SEND_BC))THEN
!             write(0,*)' PRELIM_INFO wrt_int_state%NUM_WORDS_SEND_BC already allocated for write group ',n
            ELSE
!             WRITE(0,*)' PRELIM_INFO must allocate wrt_int_state%NUM_WORDS_SEND_BC'
              ALLOCATE(wrt_int_state%NUM_WORDS_SEND_BC(1),stat=ISTAT)
              IF(ISTAT/=0)THEN
                WRITE(0,*)' Failed to allocate wrt_int_state%NUM_WORDS_SEND_BC in PRELIM_INFO'
              ENDIF
            ENDIF
!
            CALL MPI_RECV(wrt_int_state%NUM_WORDS_SEND_BC               &  !<-- Recv # of BC words
                         ,1                                             &  !<-- Data is a scalar
                         ,MPI_INTEGER                                   &  !<-- Data is integer
                         ,0                                             &  !<-- Sending task (forecast task 0)
                         ,0                                             &  !<-- MPI tag
                         ,COMM_MY_DOMAIN                                &  !<-- MPI domain communicator
                         ,JSTAT                                         &  !<-- MPI status object
                         ,IERR )
!
            IF(ALLOCATED(wrt_int_state%RST_ALL_BC_DATA))THEN
!             WRITE(0,*)' PRELIM_INFO wrt_int_state%RST_ALL_BC_DATA already allocated for write group ',n
            ELSE
!             WRITE(0,*)' PRELIM_INFO must allocate wrt_int_state%RST_ALL_BC_DATA'
              ALLOCATE(wrt_int_state%RST_ALL_BC_DATA(1:wrt_int_state%NUM_WORDS_SEND_BC(1)),stat=ISTAT)
              IF(ISTAT/=0)THEN
                WRITE(0,*)' Failed to allocate wrt_int_state%wrt_int_state%RST_ALL_BC_DATA in PRELIM_INFO'
              ENDIF
            ENDIF
!
          ENDIF
!
        ENDIF
!
!-----------------------------------------------------------------------
!***  Forecast task 0 sends the 2D data names to the lead Write task.
!-----------------------------------------------------------------------
!
        IF(MYPE_DOMAIN==0)THEN                                             !<-- Fcst task 0 sends data names
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Fcst Task0 Sends Lead Write Task 2D Data Names"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL MPI_SEND(NCHAR_I2D                                       &  !<-- Send total length of the names of 2D integer data
                       ,1                                               &  !<-- Words sent
                       ,MPI_INTEGER                                     &  !<-- Data is integer
                       ,N1                                              &  !<-- Receiving task (lead write task in group)
                       ,0                                               &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- MPI domain communicator
                       ,IERR )
!
          CALL MPI_SEND(NAMES_I2D_STRING                                &  !<-- Send names of 2D integer history variables
                       ,NCHAR_I2D(1)                                    &  !<-- Words sent
                       ,MPI_CHARACTER                                   &  !<-- Data is character
                       ,N1                                              &  !<-- Receiving task (lead write task in group)
                       ,0                                               &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- MPI domain communicator
                       ,IERR )
!
          CALL MPI_SEND(NCHAR_R2D                                       &  !<-- Send total length of the names of 2D real data
                       ,1                                               &  !<-- Words sent
                       ,MPI_INTEGER                                     &  !<-- Data is integer
                       ,N1                                              &  !<-- Receiving task (lead write task in group)
                       ,0                                               &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- MPI domain communicator
                       ,IERR )
!
          CALL MPI_SEND(NAMES_R2D_STRING                                &  !<-- Send names of 2D real history variables
                       ,NCHAR_R2D(1)                                    &  !<-- Words sent
                       ,MPI_CHARACTER                                   &  !<-- Data is character
                       ,N1                                              &  !<-- Receiving task (lead write task in group)
                       ,0                                               &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- MPI domain communicator
                       ,IERR )
!
        ELSEIF(MYPE_DOMAIN==N1)THEN                                        !<-- Lead write task receives 2D preliminary info
!
          CALL MPI_RECV(NCHAR_I2D                                       &  !<-- Recv total length of the names of 2D integer data
                       ,1                                               &  !<-- Words sent
                       ,MPI_INTEGER                                     &  !<-- Data is integer
                       ,0                                               &  !<-- Sending task (lead fcst task)
                       ,0                                               &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- MPI domain communicator
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR )
!
          CALL MPI_RECV(NAMES_I2D_STRING                                &  !<-- Recv names of 2D integer history variables
                       ,NCHAR_I2D(1)                                    &  !<-- Words sent
                       ,MPI_CHARACTER                                   &  !<-- Data is character
                       ,0                                               &  !<-- Sending task (lead fcst task)
                       ,0                                               &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- MPI domain communicator
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR )
!
          CALL MPI_RECV(NCHAR_R2D                                       &  !<-- Recv total length of the names of 2D real data
                       ,1                                               &  !<-- Words sent
                       ,MPI_INTEGER                                     &  !<-- Data is integer
                       ,0                                               &  !<-- Sending task (lead fcst task)
                       ,0                                               &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- MPI domain communicator
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR )
!
          CALL MPI_RECV(NAMES_R2D_STRING                                &  !<-- Recv names of 2D real history variables
                       ,NCHAR_R2D(1)                                    &  !<-- Words sent
                       ,MPI_CHARACTER                                   &  !<-- Data is character
                       ,0                                               &  !<-- Sending task (lead fcst task)
                       ,0                                               &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- MPI domain communicator
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        ENDIF
!
!-----------------------------------------------------------------------
!***  Each Write task must know the IDs of the Forecast tasks
!***  from which it will receive 2D gridded history data.
!-----------------------------------------------------------------------
!
        IF(MYPE_DOMAIN>=N1.AND.MYPE_DOMAIN<=N2)THEN                        !<-- The write tasks
!
          M1=LEAD_WRITE_TASK(1)                                            !<-- Local rank of lead write task
          M2=LAST_WRITE_TASK(1)                                            !<-- Local rank of last write task
          LOCAL_ID=M1+MYPE_DOMAIN-LEAD_WRITE_TASK(N)                       !<-- Local rank of every write task
!
          IF(.NOT.ALLOCATED(wrt_int_state%ID_FTASK_RECV_STA))THEN
            ALLOCATE(wrt_int_state%ID_FTASK_RECV_STA(M1:M2))
            ALLOCATE(wrt_int_state%ID_FTASK_RECV_END(M1:M2))
          ENDIF
!
          NY=0
          DO NX=M1,M2
            NY=NY+1
            CALL PARA_RANGE(JNPES(1),NWTPG,NY                           &  !<-- Find each write task's first and last rows of
                           ,JROW_FIRST,JROW_LAST)                          !    fcst tasks from which it will recv.
!
            wrt_int_state%ID_FTASK_RECV_STA(NX)=(JROW_FIRST-1)*INPES(1)    !<-- First fcst task that sends to this write task
            wrt_int_state%ID_FTASK_RECV_END(NX)=JROW_LAST*INPES(1)-1       !<-- Last fcst task that sends to this write task
          ENDDO
!
!-----------------------------------------------------------------------
!***  Each write task computes the number of words in the datastring
!***  of 2D/3D real history data it will receive from each Forecast
!***  task it is associated with.  Then allocate that datastring.
!-----------------------------------------------------------------------
!
          IF(OUTPUT_FLAG=='History')THEN
            IF(.NOT.ASSOCIATED(wrt_int_state%NUM_WORDS_RECV_R2D_HST))THEN
              N_STA=wrt_int_state%ID_FTASK_RECV_STA(LOCAL_ID)
              N_END=wrt_int_state%ID_FTASK_RECV_END(LOCAL_ID)
              ALLOCATE(wrt_int_state%NUM_WORDS_RECV_R2D_HST(N_STA:N_END))
              NUM_WORDS_RECV_R2D=>wrt_int_state%NUM_WORDS_RECV_R2D_HST
            ENDIF
!
          ELSEIF(OUTPUT_FLAG=='Restart')THEN
            IF(.NOT.ASSOCIATED(wrt_int_state%NUM_WORDS_RECV_R2D_RST))THEN
              N_STA=wrt_int_state%ID_FTASK_RECV_STA(LOCAL_ID)
              N_END=wrt_int_state%ID_FTASK_RECV_END(LOCAL_ID)
              ALLOCATE(wrt_int_state%NUM_WORDS_RECV_R2D_RST(N_STA:N_END))
              NUM_WORDS_RECV_R2D=>wrt_int_state%NUM_WORDS_RECV_R2D_RST
            ENDIF
!
          ENDIF
!
          MAX_WORDS=0
!
          DO NX=wrt_int_state%ID_FTASK_RECV_STA(LOCAL_ID)               &  !<-- The fcst tasks sending to this write task
               ,wrt_int_state%ID_FTASK_RECV_END(LOCAL_ID)                  !<--
!
            ITS=wrt_int_state%LOCAL_ISTART(NX)
            ITE=wrt_int_state%LOCAL_IEND  (NX)
            JTS=wrt_int_state%LOCAL_JSTART(NX)
            JTE=wrt_int_state%LOCAL_JEND  (NX)
!
            NUM_WORDS_TOT=(ITE-ITS+1+2*IHALO(1))*(JTE-JTS+1+2*JHALO(1)) &  !<-- # of words of 2D/3D real history data from fcst task NX
                          *KOUNT_R2D(1)
!
            NUM_WORDS_RECV_R2D(NX)=NUM_WORDS_TOT
            MAX_WORDS=MAX(MAX_WORDS                                     &  !<-- Max # of integer words from any fcst tasks
                         ,NUM_WORDS_RECV_R2D(NX))                          !
          ENDDO
!
          IF(OUTPUT_FLAG=='History')THEN
            IF(.NOT.ASSOCIATED(wrt_int_state%ALL_DATA_R2D))THEN
              ALLOCATE(wrt_int_state%ALL_DATA_R2D(MAX_WORDS),stat=ISTAT)
              IF(ISTAT/=0)THEN
                WRITE(0,*)' Write task FAILED to allocate ALL_DATA_R2D'
                CALL ESMF_Finalize(endflag=ESMF_END_ABORT             &
                                  ,rc             =RC  )
              ENDIF
              ALL_DATA_R2D=>wrt_int_state%ALL_DATA_R2D
            ENDIF
!
          ELSEIF(OUTPUT_FLAG=='Restart')THEN
            IF(.NOT.ASSOCIATED(wrt_int_state%RST_ALL_DATA_R2D))THEN
              ALLOCATE(wrt_int_state%RST_ALL_DATA_R2D(MAX_WORDS),stat=ISTAT)
              IF(ISTAT/=0)THEN
                WRITE(0,*)' Write task FAILED to allocate RST_ALL_DATA_R2D'
                CALL ESMF_Finalize(endflag=ESMF_END_ABORT             &
                                  ,rc             =RC  )
              ENDIF
              ALL_DATA_R2D=>wrt_int_state%RST_ALL_DATA_R2D
            ENDIF
!
          ENDIF
!
!-----------------------------------------------------------------------
!***  Likewise for the 2-D integer data.
!-----------------------------------------------------------------------
!
          N_STA=wrt_int_state%ID_FTASK_RECV_STA(LOCAL_ID)
          N_END=wrt_int_state%ID_FTASK_RECV_END(LOCAL_ID)
!
          IF(OUTPUT_FLAG=='History')THEN
            IF(.NOT.ASSOCIATED(wrt_int_state%NUM_WORDS_RECV_I2D_HST))THEN
              ALLOCATE(wrt_int_state%NUM_WORDS_RECV_I2D_HST(N_STA:N_END))
              NUM_WORDS_RECV_I2D=>wrt_int_state%NUM_WORDS_RECV_I2D_HST
            ENDIF
          ELSEIF(OUTPUT_FLAG=='Restart')THEN
            IF(.NOT.ASSOCIATED(wrt_int_state%NUM_WORDS_RECV_I2D_RST))THEN
              ALLOCATE(wrt_int_state%NUM_WORDS_RECV_I2D_RST(N_STA:N_END))
              NUM_WORDS_RECV_I2D=>wrt_int_state%NUM_WORDS_RECV_I2D_RST
            ENDIF
          ENDIF
!
          MAX_WORDS=0
!
          DO NX=wrt_int_state%ID_FTASK_RECV_STA(LOCAL_ID)               &  !<-- The fcst tasks sending to this write task
              ,wrt_int_state%ID_FTASK_RECV_END(LOCAL_ID)                   !
!
            ITS=wrt_int_state%LOCAL_ISTART(NX)
            ITE=wrt_int_state%LOCAL_IEND  (NX)
            JTS=wrt_int_state%LOCAL_JSTART(NX)
            JTE=wrt_int_state%LOCAL_JEND  (NX)
!
            NUM_WORDS_TOT=(ITE-ITS+1+2*IHALO(1))*(JTE-JTS+1+2*JHALO(1)) &  !<-- # of words of 2D integer history data from fcst task NX
                          *KOUNT_I2D(1)
!
            NUM_WORDS_RECV_I2D(NX)=NUM_WORDS_TOT
            MAX_WORDS=MAX(MAX_WORDS                                     &  !<-- Max # of real words from all fcst tasks
                         ,NUM_WORDS_RECV_I2D(NX))                          !
          ENDDO
!
          IF(OUTPUT_FLAG=='History')THEN
            IF(.NOT.ASSOCIATED(wrt_int_state%ALL_DATA_I2D))THEN
              ALLOCATE(wrt_int_state%ALL_DATA_I2D(MAX_WORDS),stat=ISTAT)
              IF(ISTAT/=0)THEN
                WRITE(0,*)' Write task FAILED to allocate ALL_DATA_I2D'
                CALL ESMF_Finalize(endflag=ESMF_END_ABORT             &
                                  ,rc             =RC  )
              ENDIF
              ALL_DATA_I2D=>wrt_int_state%ALL_DATA_I2D
            ENDIF
!
          ELSEIF(OUTPUT_FLAG=='Restart')THEN
            IF(.NOT.ASSOCIATED(wrt_int_state%RST_ALL_DATA_I2D))THEN
              ALLOCATE(wrt_int_state%RST_ALL_DATA_I2D(MAX_WORDS),stat=ISTAT)
              IF(ISTAT/=0)THEN
                WRITE(0,*)' Write task FAILED to allocate RST_ALL_DATA_I2D'
                CALL ESMF_Finalize(endflag=ESMF_END_ABORT             &
                                  ,rc             =RC  )
              ENDIF
              ALL_DATA_I2D=>wrt_int_state%RST_ALL_DATA_I2D
            ENDIF
!
          ENDIF
!
!-----------------------------------------------------------------------
!***  Each Write task also must know the North-South extent of the
!***  full 2D domain that it will handle.  This is determined by
!***  the coverage of the Fcst tasks that send to it.
!-----------------------------------------------------------------------
!
          JSTA_WRITE=wrt_int_state%LOCAL_JSTART(wrt_int_state%ID_FTASK_RECV_STA(LOCAL_ID))  !<-- JTS/JTE of the first fcst task
          JEND_WRITE=wrt_int_state%LOCAL_JEND  (wrt_int_state%ID_FTASK_RECV_END(LOCAL_ID))  !    that sends to this write task.
!
!-----------------------------------------------------------------------
!***  Now each Write task allocates its own section of the 2D domain
!***  for all the 2D variables it will receive and its 1D equivalent
!***  used to transfer the data to the lead Write task.
!-----------------------------------------------------------------------
!
          LENGTH=IM*(JEND_WRITE-JSTA_WRITE+1)
!
          IF(OUTPUT_FLAG=='History')THEN
!
            ALLOCATE(wrt_int_state%WRITE_SUBSET_I(1:IM,JSTA_WRITE:JEND_WRITE  &
                                                 ,KOUNT_I2D(1)))
            ALLOCATE(wrt_int_state%BUFF_INT(LENGTH))
            WRITE_SUBSET_I=>wrt_int_state%WRITE_SUBSET_I
            BUFF_INT=>wrt_int_state%BUFF_INT
!
            ALLOCATE(wrt_int_state%WRITE_SUBSET_R(1:IM,JSTA_WRITE:JEND_WRITE  &
                                                 ,KOUNT_R2D(1)))
            ALLOCATE(wrt_int_state%BUFF_REAL(LENGTH))
            WRITE_SUBSET_R=>wrt_int_state%WRITE_SUBSET_R
            BUFF_REAL=>wrt_int_state%BUFF_REAL
!
          ELSEIF(OUTPUT_FLAG=='Restart')THEN
!
            ALLOCATE(wrt_int_state%RST_WRITE_SUBSET_I(1:IM,JSTA_WRITE:JEND_WRITE  &
                                                     ,KOUNT_I2D(1)))
            ALLOCATE(wrt_int_state%RST_BUFF_INT(LENGTH))
            WRITE_SUBSET_I=>wrt_int_state%RST_WRITE_SUBSET_I
            BUFF_INT=>wrt_int_state%RST_BUFF_INT
!
            ALLOCATE(wrt_int_state%RST_WRITE_SUBSET_R(1:IM,JSTA_WRITE:JEND_WRITE  &
                                                     ,KOUNT_R2D(1)))
            ALLOCATE(wrt_int_state%RST_BUFF_REAL(LENGTH))
            WRITE_SUBSET_R=>wrt_int_state%RST_WRITE_SUBSET_R
            BUFF_REAL=>wrt_int_state%RST_BUFF_REAL
!
          ENDIF
!       
        ENDIF
!
!-----------------------------------------------------------------------
!***  The lead write task allocates its working arrays into which
!***  it will assemble each individual 2D field that will be
!***  written to the history files.
!-----------------------------------------------------------------------
!
        IF(MYPE_DOMAIN==LEAD_WRITE_TASK(N))THEN
!
          IF(OUTPUT_FLAG=='History')THEN
            ALLOCATE(wrt_int_state%OUTPUT_ARRAY_I2D(1:IM,1:JM))
            ALLOCATE(wrt_int_state%OUTPUT_ARRAY_R2D(1:IM,1:JM))
            OUTPUT_ARRAY_I2D=>wrt_int_state%OUTPUT_ARRAY_I2D
            OUTPUT_ARRAY_R2D=>wrt_int_state%OUTPUT_ARRAY_R2D
!
          ELSEIF(OUTPUT_FLAG=='Restart')THEN
            ALLOCATE(wrt_int_state%RST_OUTPUT_ARRAY_I2D(1:IM,1:JM))
            ALLOCATE(wrt_int_state%RST_OUTPUT_ARRAY_R2D(1:IM,1:JM))
            OUTPUT_ARRAY_I2D=>wrt_int_state%RST_OUTPUT_ARRAY_I2D
            OUTPUT_ARRAY_R2D=>wrt_int_state%RST_OUTPUT_ARRAY_R2D
!
          ENDIF
!
        ENDIF
!
!-----------------------------------------------------------------------
!***  Since all scalar/1D data is identical on all forecast tasks,
!***  task 0 alone can send the information to the lead write task
!***  that will later write it to the history file.
!--------------------------------------------------------------------
!
!--------------------------------------------------------------------
        task_0_sends: IF(MYPE_DOMAIN==0)THEN                            !<-- Forecast task 0 sends
!--------------------------------------------------------------------
!
!------------------------------------------------
!***  Send scalar/1D integer history information.
!------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Fcst Task0 Sends Scalar/1D Integer History Data"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL MPI_SEND(KOUNT_I1D                                       &  !<-- Send # of scalar/1D integer history variables
                       ,1                                               &  !<-- Words sent
                       ,MPI_INTEGER                                     &  !<-- Data is integer
                       ,LEAD_WRITE_TASK(N)                              &  !<-- Receiving task (1st write task in group)
                       ,0                                               &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- MPI domain communicator
                       ,IERR )
!
          CALL MPI_SEND(LENGTH_SUM_I1D                                  &  !<-- Send length of string of all such integer history variables
                       ,1                                               &  !<-- Words sent
                       ,MPI_INTEGER                                     &  !<-- Data is integer
                       ,LEAD_WRITE_TASK(N)                              &  !<-- Receiving task (1st write task in group)
!                      ,0                                               &  !<-- MPI tag
                       ,678                                             &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- MPI domain commumicator
                       ,IERR )
!
          CALL MPI_SEND(LENGTH_DATA_I1D                                 &  !<-- Send lengths of each scalar/1D integer history variable
                       ,KOUNT_I1D(1)                                    &  !<-- Words sent
                       ,MPI_INTEGER                                     &  !<-- Data is integer
                       ,LEAD_WRITE_TASK(N)                              &  !<-- Receiving task (1st write task in group)
                       ,0                                               &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- MPI domain commumicator
                       ,IERR )
!
          CALL MPI_SEND(NAMES_I1D_STRING                                &  !<-- Send names of each scalar/1D integer history variable
                       ,NCHAR_I1D                                       &  !<-- Words sent
                       ,MPI_CHARACTER                                   &  !<-- Words are character
                       ,LEAD_WRITE_TASK(N)                              &  !<-- Receiving task (1st write task in group)
                       ,0                                               &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- MPI domain commumicator
                       ,IERR )
!
          CALL MPI_SEND(ALL_DATA_I1D                                    &  !<-- Send the full string of all scalar/1D integer history data
                       ,LENGTH_SUM_I1D(1)                               &  !<-- Words sent
                       ,MPI_INTEGER                                     &  !<-- Data is integer
                       ,LEAD_WRITE_TASK(N)                              &  !<-- Receiving task (1st write task in group)
                       ,0                                               &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- MPI domain commumicator
                       ,IERR )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!---------------------------------------------
!***  Send scalar/1D real history information.
!---------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Fcst Task 0 Sends Scalar/1D Real History Data"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL MPI_SEND(KOUNT_R1D                                       &  !<-- Send # of scalar/1D real history variables
                       ,1                                               &  !<-- Words sent
                       ,MPI_INTEGER                                     &  !<-- Data is integer
                       ,LEAD_WRITE_TASK(N)                              &  !<-- Receiving task (1st write task in group)
                       ,0                                               &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- MPI domain commumicator
                       ,IERR )
!
          CALL MPI_SEND(LENGTH_SUM_R1D                                  &  !<-- Send length of string of all such real history variables
                       ,1                                               &  !<-- Words sent
                       ,MPI_INTEGER                                     &  !<-- Data is integer     
                       ,LEAD_WRITE_TASK(N)                              &  !<-- Receiving task (1st write task in group)
                       ,0                                               &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- MPI domain commumicator
                       ,IERR )
!
          CALL MPI_SEND(LENGTH_DATA_R1D                                 &  !<-- Send lengths of each scalar/1D real history variable
                       ,KOUNT_R1D(1)                                    &  !<-- Words sent
                       ,MPI_INTEGER                                     &  !<-- Data is integer
                       ,LEAD_WRITE_TASK(N)                              &  !<-- Receiving task (1st write task in group)
                       ,0                                               &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- MPI domain commumicator
                       ,IERR )
!
          CALL MPI_SEND(NAMES_R1D_STRING                                &  !<-- Send names of each scalar/1D real history variable
                       ,NCHAR_R1D                                       &  !<-- Words sent
                       ,MPI_CHARACTER                                   &  !<-- Words are character
                       ,LEAD_WRITE_TASK(N)                              &  !<-- Receiving task (1st write task in group)
                       ,0                                               &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- MPI domain commumicator
                       ,IERR )
!
          CALL MPI_SEND(ALL_DATA_R1D                                    &  !<-- Send the full string of all scalar/1D real history data
                       ,LENGTH_SUM_R1D(1)                               &  !<-- Words sent
                       ,MPI_REAL                                        &  !<-- Words are real
                       ,LEAD_WRITE_TASK(N)                              &  !<-- Receiving task (1st write task in group)
                       ,0                                               &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- MPI domain commumicator
                       ,IERR )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!--------------------------------------
!***  Send logical history information.
!--------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Fcst Task0 Sends Logical History Data"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL MPI_SEND(KOUNT_LOG                                       &  !<-- Send # of logical history variables
                       ,1                                               &  !<-- Words sent
                       ,MPI_INTEGER                                     &  !<-- Data is integer
                       ,LEAD_WRITE_TASK(N)                              &  !<-- Receiving task (1st write task in group)
                       ,0                                               &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- MPI domain commumicator
                       ,IERR )
!
          CALL MPI_SEND(LENGTH_SUM_LOG                                  &  !<-- Send length of string of all logical variables
                       ,1                                               &  !<-- Words sent
                       ,MPI_INTEGER                                     &  !<-- Data is integer
                       ,LEAD_WRITE_TASK(N)                              &  !<-- Receiving task (1st write task in group)
                       ,0                                               &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- MPI domain commumicator
                       ,IERR )
!
          CALL MPI_SEND(NAMES_LOG_STRING                                &  !<-- Send names of each logical history variable
                       ,NCHAR_LOG                                       &  !<-- Words sent
                       ,MPI_CHARACTER                                   &  !<-- Data is character
                       ,LEAD_WRITE_TASK(N)                              &  !<-- Receiving task (1st write task in group)
                       ,0                                               &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- MPI domain commumicator
                       ,IERR )
!
          CALL MPI_SEND(ALL_DATA_LOG                                    &  !<-- Send the full string of all logical history data
                       ,LENGTH_SUM_LOG(1)                               &  !<-- Words sent
                       ,MPI_LOGICAL                                     &  !<-- Data is logical 
                       ,LEAD_WRITE_TASK(N)                              &  !<-- Receiving task (1st write task in group)
                       ,0                                               &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- MPI domain commumicator
                       ,IERR )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
        ENDIF task_0_sends
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
        write_task_recvs: IF(MYPE_DOMAIN==LEAD_WRITE_TASK(N))THEN          !<-- 1st write task in this group receives
                                                                           !    all of the data just sent to it by
                                                                           !    fcst task 0.
!-----------------------------------------------------------------------
!***  Receive scalar/1D integer history information
!***  from Forecast task 0.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Write Tasks Recv Scalar/1D Integer History Data"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL MPI_RECV(KOUNT_I1D                                       &  !<-- Recv # of scalar/1D integer history variables
                       ,1                                               &  !<-- Words received
                       ,MPI_INTEGER                                     &  !<-- Data is integer
                       ,0                                               &  !<-- Sending task (lead fcst task)
                       ,0                                               &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- MPI domain commumicator
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR )
!
          CALL MPI_RECV(LENGTH_SUM_I1D                                  &  !<-- Recv length of string of all such integer history variables
                       ,1                                               &  !<-- Words received
                       ,MPI_INTEGER                                     &  !<-- Data is integer
                       ,0                                               &  !<-- Sending task (lead fcst task)
                       ,678                                             &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- MPI domain commumicator
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR )
!
          CALL MPI_RECV(LENGTH_DATA_I1D                                 &  !<-- Recv lengths of each scalar/1D integer history variable
                       ,KOUNT_I1D(1)                                    &  !<-- Words received
                       ,MPI_INTEGER                                     &  !<-- Data is integer
                       ,0                                               &  !<-- Sending task (lead fcst task)
                       ,0                                               &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- MPI domain commumicator
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR )
!
          NCHAR_I1D=KOUNT_I1D(1)*ESMF_MAXSTR
!
          CALL MPI_RECV(NAMES_I1D_STRING                                &  !<-- Recv names of each scalar/1D integer history variable
                       ,NCHAR_I1D                                       &  !<-- Words received
                       ,MPI_CHARACTER                                   &  !<-- Words are character
                       ,0                                               &  !<-- Sending task (lead fcst task)
                       ,0                                               &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- MPI domain commumicator
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR )
!
          CALL MPI_RECV(ALL_DATA_I1D                                    &  !<-- Recv the full string of all scalar/1D integer history data
                       ,LENGTH_SUM_I1D(1)                               &  !<-- Words received
                       ,MPI_INTEGER                                     &  !<-- Data is integer
                       ,0                                               &  !<-- Sending task (lead fcst task)
                       ,0                                               &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- MPI domain commumicator
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Receive scalar/1D real history information.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Write Tasks Recv Scalar/1D Real History Data"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL MPI_RECV(KOUNT_R1D                                       &  !<-- Recv # of scalar/1D real history variables
                       ,1                                               &  !<-- Words received
                       ,MPI_INTEGER                                     &  !<-- Data is integer
                       ,0                                               &  !<-- Sending task (lead fcst task)
                       ,0                                               &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- MPI domain commumicator
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR )
!
          CALL MPI_RECV(LENGTH_SUM_R1D                                  &  !<-- Recv length of string of all such real history variables
                       ,1                                               &  !<-- Words received
                       ,MPI_INTEGER                                     &  !<-- Data is integer
                       ,0                                               &  !<-- Sending task (lead fcst task)
                       ,0                                               &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- MPI domain commumicator
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR )
!
          CALL MPI_RECV(LENGTH_DATA_R1D                                 &  !<-- Recv lengths of each scalar/1D real history variable
                       ,KOUNT_R1D(1)                                    &  !<-- Words received
                       ,MPI_INTEGER                                     &  !<-- Data is integer
                       ,0                                               &  !<-- Sending task (lead fcst task)
                       ,0                                               &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- MPI domain commumicator
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR )
!
          NCHAR_R1D=KOUNT_R1D(1)*ESMF_MAXSTR
          CALL MPI_RECV(NAMES_R1D_STRING                                &  !<-- Recv names of scalar/1D real history variables
                       ,NCHAR_R1D                                       &  !<-- Words received
                       ,MPI_CHARACTER                                   &  !<-- Data is character
                       ,0                                               &  !<-- Sending task (lead fcst task)
                       ,0                                               &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- MPI domain commumicator
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR )
!
          CALL MPI_RECV(ALL_DATA_R1D                                    &  !<-- Recv the string of all scalar/1D real history data
                       ,LENGTH_SUM_R1D(1)                               &  !<-- Words received
                       ,MPI_REAL                                        &  !<-- Data is real
                       ,0                                               &  !<-- Sending task (lead fcst task)
                       ,0                                               &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- MPI domain commumicator
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Receive logical history information.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Write Tasks Recv Logical Real History Data"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL MPI_RECV(KOUNT_LOG                                       &  !<-- Recv # of logical history variables
                       ,1                                               &  !<-- Words received
                       ,MPI_INTEGER                                     &  !<-- Data is integer  
                       ,0                                               &  !<-- Sending task (lead fcst task)
                       ,0                                               &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- MPI domain commumicator
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR )
!
          CALL MPI_RECV(LENGTH_SUM_LOG                                  &  !<-- Recv length of string of all logical history variables
                       ,1                                               &  !<-- Words received
                       ,MPI_INTEGER                                     &  !<-- Data is integer  
                       ,0                                               &  !<-- Sending task (lead fcst task)
                       ,0                                               &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- MPI domain commumicator
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR )
!
          NCHAR_LOG=KOUNT_LOG(1)*ESMF_MAXSTR
          CALL MPI_RECV(NAMES_LOG_STRING                                &  !<-- Recv names of logical history variables
                       ,NCHAR_LOG                                       &  !<-- Words received
                       ,MPI_CHARACTER                                   &  !<-- Data is character
                       ,0                                               &  !<-- Sending task (lead fcst task)
                       ,0                                               &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- MPI domain commumicator
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR )
!
          CALL MPI_RECV(ALL_DATA_LOG                                    &  !<-- Recv the string of all logical history data
                       ,LENGTH_SUM_LOG(1)                               &  !<-- Words received
                       ,MPI_LOGICAL                                     &  !<-- Data is logical  
                       ,0                                               &  !<-- Sending task (lead fcst task)
                       ,0                                               &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- MPI domain commumicator
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
        ENDIF write_task_recvs
!
!-----------------------------------------------------------------------
!
      ENDDO n_groups_3
!
!-----------------------------------------------------------------------
!***  When run post on quilt, all the write tasks need to know the
!***  variables in history bundle that first forecast task sent to
!***  lead write tasks.
!-----------------------------------------------------------------------
!
      IF(WRT_INT_STATE%WRITE_DOPOST.and.OUTPUT_FLAG=='History')THEN
!
      n_groups_4: DO N=1,NUM_WRITE_GROUPS
!-----------------------------------------------------------------------
!
        N1=LEAD_WRITE_TASK(N)
        N2=LAST_WRITE_TASK(N)
!
        IF(MYPE_DOMAIN>=N1.AND.MYPE_DOMAIN<=N2)THEN

!          write(0,*)'bf broadcast NCHAR_I2D'
          CALL MPI_BCAST(NCHAR_I2D,1,MPI_INTEGER,0,MPI_COMM_COMP,IERR)
!          write(0,*)'af broadcast NCHAR_I2D=',NCHAR_I2D
          CALL MPI_BCAST(NAMES_I2D_STRING,NCHAR_I2D(1),MPI_CHARACTER,0,   &
                         MPI_COMM_COMP,IERR)
!          write(0,*)'af broadcast NAMES_I2D_STRING=',NAMES_I2D_STRING(1:30)
          CALL MPI_BCAST(NCHAR_R2D,1,MPI_INTEGER,0,MPI_COMM_COMP,IERR)
!          write(0,*)'af broadcast NCHAR_R2D=',NCHAR_R2D
          CALL MPI_BCAST(NAMES_R2D_STRING,NCHAR_R2D(1),MPI_CHARACTER,0,   &
                         MPI_COMM_COMP,IERR)
!integer
          CALL MPI_BCAST(KOUNT_I1D,1,MPI_INTEGER,0,MPI_COMM_COMP,IERR)
          CALL MPI_BCAST(LENGTH_SUM_I1D,1,MPI_INTEGER,0,MPI_COMM_COMP,   &
                          IERR)
!          write(0,*)'af broadcast LENGTH_SUM_I1D=',LENGTH_SUM_I1D(1)
          CALL MPI_BCAST(LENGTH_DATA_I1D,KOUNT_I1D(1),MPI_INTEGER,0,     &
                          MPI_COMM_COMP,IERR)
!          write(0,*)'af broadcast LENGTH_DATA_I1D=',LENGTH_DATA_I1D(1)
          NCHAR_I1D=KOUNT_I1D(1)*ESMF_MAXSTR
          CALL MPI_BCAST(NAMES_I1D_STRING,NCHAR_I1D,MPI_CHARACTER,0,     &
                          MPI_COMM_COMP,IERR)
!          write(0,*)'af broadcast NAMES_I1D_STRING=',NAMES_I1D_STRING(1:30)
          CALL MPI_BCAST(ALL_DATA_I1D,LENGTH_SUM_I1D(1),MPI_LOGICAL,0,   &
                          MPI_COMM_COMP,IERR)
!          write(0,*)'af broadcast ALL_DATA_I1D=',ALL_DATA_I1D(1:5)
!real
          CALL MPI_BCAST(KOUNT_R1D,1,MPI_INTEGER,0,MPI_COMM_COMP,IERR)
!          write(0,*)'af broadcast NCHAR_R1D=',KOUNT_R1D
          CALL MPI_BCAST(LENGTH_SUM_R1D,1,MPI_INTEGER,0,MPI_COMM_COMP,   &
                          IERR)
          CALL MPI_BCAST(LENGTH_DATA_R1D,KOUNT_R1D(1),MPI_INTEGER,0,     &
                          MPI_COMM_COMP,IERR)
          NCHAR_R1D=KOUNT_R1D(1)*ESMF_MAXSTR
          CALL MPI_BCAST(NAMES_R1D_STRING,NCHAR_R1D,MPI_CHARACTER,0,     &
                          MPI_COMM_COMP,IERR)
          CALL MPI_BCAST(ALL_DATA_R1D,LENGTH_SUM_R1D(1),MPI_LOGICAL,0,   &
                 MPI_COMM_COMP,IERR)
!logical
          CALL MPI_BCAST(KOUNT_LOG,1,MPI_INTEGER,0,MPI_COMM_COMP,IERR)
!          write(0,*)'af broadcast NCHAR_LOG=',KOUNT_LOG
          CALL MPI_BCAST(LENGTH_SUM_LOG,1,MPI_INTEGER,0,MPI_COMM_COMP,   &
                          IERR)
          NCHAR_LOG=KOUNT_LOG(1)*ESMF_MAXSTR
          CALL MPI_BCAST(NAMES_LOG_STRING,NCHAR_LOG,MPI_CHARACTER,0,     &
                          MPI_COMM_COMP,IERR)
          CALL MPI_BCAST(ALL_DATA_LOG,LENGTH_SUM_LOG(1),MPI_LOGICAL,0,   &
                          MPI_COMM_COMP,IERR)
!
          ENDIF
!
        ENDDO n_groups_4
!
      ENDIF
!
!-----------------------------------------------------------------------
!
      DEALLOCATE(INPES)
      DEALLOCATE(JNPES)
      DEALLOCATE(IHALO)
      DEALLOCATE(JHALO)
      DEALLOCATE(NCHAR_I2D)
      DEALLOCATE(NCHAR_R2D)
!
!-----------------------------------------------------------------------
!***  All the write tasks now have everything they need with respect to
!***  the nature of the data they will be writing out whenever their
!***  particular write group is invoked.  However for the forecast
!***  tasks only the variables in the internal state of the Write 
!***  component associated with write group #1 have been filled.
!***  So now the write tasks fill those same variables in each of 
!***  the other internal states of the Write components associated 
!***  with the remaining write groups.
!-----------------------------------------------------------------------
!
      IF(MYPE_DOMAIN<=LAST_FCST_TASK                                    &  !<-- Select only the forecast tasks
                 .AND.                                                  &
         NUM_WRITE_GROUPS>1)THEN                                           !<-- Variables in write group #1 already done
!
!-----------------------------------------------------------------------
!
        ALLOCATE(WRT_INT_STATE_X(2:NUM_WRITE_GROUPS),stat=RC)              !<-- Allocate working array of internal states
!
        IF(RC/=0)THEN
          WRITE(0,*)' Failed to allocate working array of'              &
                   ,' internal state pointers!'
          WRITE(0,*)' ABORT!'
          CALL ESMF_Finalize(endflag=ESMF_END_ABORT                 &
                            ,rc             =RC )
        ENDIF
!
!-----------------------------------------------------------------------
!
        DO N=2,NUM_WRITE_GROUPS                                            !<-- Loop through remaining write groups
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="PRELIM_INFO: Fcst Tasks Get Nth Write Internal State"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_GridCompGetInternalState(WRITE_COMPS(N)             &
                                            ,WRAP                       &
                                            ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          WRT_INT_STATE_X(N)%LOC_INT_STATE=>wrap%WRITE_INT_STATE           !<-- Pointer to internal state of Nth Write component
!
!-----------------------------------------------------------------------
!
          CALL FILL_GROUP_STATES                                           !<-- Fill variables in internal state N from state #1
!
!-----------------------------------------------------------------------
!
        ENDDO
!
        DEALLOCATE(WRT_INT_STATE_X)
!
      ENDIF
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      SUBROUTINE POINT_LOCAL
!
!-----------------------------------------------------------------------
!
      IF(OUTPUT_FLAG=='History')THEN
!
        ALLOCATE(wrt_int_state%NUM_WORDS_SEND_I2D_HST)
        ALLOCATE(wrt_int_state%NUM_WORDS_SEND_R2D_HST)
        ALLOCATE(wrt_int_state%FIELD_NAME(1:5000))
        ALLOCATE(wrt_int_state%NAMES_I1D_STRING)
        ALLOCATE(wrt_int_state%NAMES_I2D_STRING)
        ALLOCATE(wrt_int_state%NAMES_R1D_STRING)
        ALLOCATE(wrt_int_state%NAMES_R2D_STRING)
        ALLOCATE(wrt_int_state%NAMES_LOG_STRING)
!
        NUM_WORDS_SEND_I2D=>wrt_int_state%NUM_WORDS_SEND_I2D_HST
        NUM_WORDS_SEND_R2D=>wrt_int_state%NUM_WORDS_SEND_R2D_HST
        KOUNT_I1D=>wrt_int_state%KOUNT_I1D
        KOUNT_I2D=>wrt_int_state%KOUNT_I2D
        KOUNT_R1D=>wrt_int_state%KOUNT_R1D
        KOUNT_R2D=>wrt_int_state%KOUNT_R2D
        KOUNT_LOG=>wrt_int_state%KOUNT_LOG
        LENGTH_DATA_I1D=>wrt_int_state%LENGTH_DATA_I1D
        LENGTH_DATA_R1D=>wrt_int_state%LENGTH_DATA_R1D
        LENGTH_DATA_R2D=>wrt_int_state%LENGTH_DATA_R2D
        LENGTH_SUM_I1D=>wrt_int_state%LENGTH_SUM_I1D
        LENGTH_SUM_R1D=>wrt_int_state%LENGTH_SUM_R1D
        LENGTH_SUM_R2D=>wrt_int_state%LENGTH_SUM_R2D
        LENGTH_SUM_LOG=>wrt_int_state%LENGTH_SUM_LOG
        NCOUNT_FIELDS=>wrt_int_state%NCOUNT_FIELDS
        ALL_DATA_I1D=>wrt_int_state%ALL_DATA_I1D
        ALL_DATA_R1D=>wrt_int_state%ALL_DATA_R1D
        FIELD_NAME=>wrt_int_state%FIELD_NAME
        NAMES_I1D_STRING=>wrt_int_state%NAMES_I1D_STRING
        NAMES_I2D_STRING=>wrt_int_state%NAMES_I2D_STRING
        NAMES_R1D_STRING=>wrt_int_state%NAMES_R1D_STRING
        NAMES_R2D_STRING=>wrt_int_state%NAMES_R2D_STRING
        NAMES_LOG_STRING=>wrt_int_state%NAMES_LOG_STRING
        ALL_DATA_LOG=>wrt_int_state%ALL_DATA_LOG
!
      ELSEIF(OUTPUT_FLAG=='Restart')THEN
!
        ALLOCATE(wrt_int_state%NUM_WORDS_SEND_I2D_RST)
        ALLOCATE(wrt_int_state%NUM_WORDS_SEND_R2D_RST)
        ALLOCATE(wrt_int_state%RST_FIELD_NAME(1:5000))
        ALLOCATE(wrt_int_state%RST_NAMES_I1D_STRING)
        ALLOCATE(wrt_int_state%RST_NAMES_I2D_STRING)
        ALLOCATE(wrt_int_state%RST_NAMES_R1D_STRING)
        ALLOCATE(wrt_int_state%RST_NAMES_R2D_STRING)
        ALLOCATE(wrt_int_state%RST_NAMES_LOG_STRING)
!
        NUM_WORDS_SEND_I2D=>wrt_int_state%NUM_WORDS_SEND_I2D_RST
        NUM_WORDS_SEND_R2D=>wrt_int_state%NUM_WORDS_SEND_R2D_RST
        KOUNT_I1D=>wrt_int_state%RST_KOUNT_I1D
        KOUNT_I2D=>wrt_int_state%RST_KOUNT_I2D
        KOUNT_R1D=>wrt_int_state%RST_KOUNT_R1D
        KOUNT_R2D=>wrt_int_state%RST_KOUNT_R2D
        KOUNT_LOG=>wrt_int_state%RST_KOUNT_LOG
        LENGTH_DATA_I1D=>wrt_int_state%RST_LENGTH_DATA_I1D
        LENGTH_DATA_R1D=>wrt_int_state%RST_LENGTH_DATA_R1D
        LENGTH_DATA_R2D=>wrt_int_state%RST_LENGTH_DATA_R2D
        LENGTH_SUM_I1D=>wrt_int_state%RST_LENGTH_SUM_I1D
        LENGTH_SUM_R1D=>wrt_int_state%RST_LENGTH_SUM_R1D
        LENGTH_SUM_R2D=>wrt_int_state%RST_LENGTH_SUM_R2D
        LENGTH_SUM_LOG=>wrt_int_state%RST_LENGTH_SUM_LOG
        NCOUNT_FIELDS=>wrt_int_state%RST_NCOUNT_FIELDS
        ALL_DATA_I1D=>wrt_int_state%RST_ALL_DATA_I1D
        ALL_DATA_R1D=>wrt_int_state%RST_ALL_DATA_R1D
        FIELD_NAME=>wrt_int_state%RST_FIELD_NAME
        NAMES_I1D_STRING=>wrt_int_state%RST_NAMES_I1D_STRING
        NAMES_I2D_STRING=>wrt_int_state%RST_NAMES_I2D_STRING
        NAMES_R1D_STRING=>wrt_int_state%RST_NAMES_R1D_STRING
        NAMES_R2D_STRING=>wrt_int_state%RST_NAMES_R2D_STRING
        NAMES_LOG_STRING=>wrt_int_state%RST_NAMES_LOG_STRING
        ALL_DATA_LOG=>wrt_int_state%RST_ALL_DATA_LOG
!
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE POINT_LOCAL
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
      SUBROUTINE FILL_GROUP_STATES
!
!-----------------------------------------------------------------------
!
      IF(OUTPUT_FLAG=='History')THEN
!
        ALLOCATE(wrt_int_state_x(N)%loc_int_state%NUM_WORDS_SEND_I2D_HST)
        ALLOCATE(wrt_int_state_x(N)%loc_int_state%NUM_WORDS_SEND_R2D_HST)
        ALLOCATE(wrt_int_state_x(N)%loc_int_state%FIELD_NAME(1:5000))
        ALLOCATE(wrt_int_state_x(N)%loc_int_state%NAMES_I1D_STRING)
        ALLOCATE(wrt_int_state_x(N)%loc_int_state%NAMES_I2D_STRING)
        ALLOCATE(wrt_int_state_x(N)%loc_int_state%NAMES_R1D_STRING)
        ALLOCATE(wrt_int_state_x(N)%loc_int_state%NAMES_R2D_STRING)
        ALLOCATE(wrt_int_state_x(N)%loc_int_state%NAMES_LOG_STRING)
        ALLOCATE(wrt_int_state_x(N)%loc_int_state%ALL_DATA_I2D(wrt_int_state%NUM_WORDS_SEND_I2D_HST),stat=ISTAT)
        ALLOCATE(wrt_int_state_x(N)%loc_int_state%ALL_DATA_R2D(wrt_int_state%NUM_WORDS_SEND_R2D_HST),stat=ISTAT)
!
        wrt_int_state_x(N)%loc_int_state%INPES=wrt_int_state%INPES
        wrt_int_state_x(N)%loc_int_state%JNPES=wrt_int_state%JNPES
        wrt_int_state_x(N)%loc_int_state%IHALO=wrt_int_state%IHALO
        wrt_int_state_x(N)%loc_int_state%JHALO=wrt_int_state%JHALO
        wrt_int_state_x(N)%loc_int_state%NUM_WORDS_SEND_I2D_HST=wrt_int_state%NUM_WORDS_SEND_I2D_HST
        wrt_int_state_x(N)%loc_int_state%NUM_WORDS_SEND_R2D_HST=wrt_int_state%NUM_WORDS_SEND_R2D_HST
        wrt_int_state_x(N)%loc_int_state%KOUNT_I1D=wrt_int_state%KOUNT_I1D
        wrt_int_state_x(N)%loc_int_state%KOUNT_I2D=wrt_int_state%KOUNT_I2D
        wrt_int_state_x(N)%loc_int_state%KOUNT_R1D=wrt_int_state%KOUNT_R1D
        wrt_int_state_x(N)%loc_int_state%KOUNT_R2D=wrt_int_state%KOUNT_R2D
        wrt_int_state_x(N)%loc_int_state%KOUNT_LOG=wrt_int_state%KOUNT_LOG
        wrt_int_state_x(N)%loc_int_state%LENGTH_DATA_I1D=wrt_int_state%LENGTH_DATA_I1D
        wrt_int_state_x(N)%loc_int_state%LENGTH_DATA_R1D=wrt_int_state%LENGTH_DATA_R1D
!       wrt_int_state_x(N)%loc_int_state%LENGTH_DATA_R2D=wrt_int_state%LENGTH_DATA_R2D
        wrt_int_state_x(N)%loc_int_state%LENGTH_SUM_I1D=wrt_int_state%LENGTH_SUM_I1D
        wrt_int_state_x(N)%loc_int_state%LENGTH_SUM_R1D=wrt_int_state%LENGTH_SUM_R1D
!       wrt_int_state_x(N)%loc_int_state%LENGTH_SUM_R2D=wrt_int_state%LENGTH_SUM_R2D
        wrt_int_state_x(N)%loc_int_state%LENGTH_SUM_LOG=wrt_int_state%LENGTH_SUM_LOG
        wrt_int_state_x(N)%loc_int_state%NCOUNT_FIELDS=wrt_int_state%NCOUNT_FIELDS
        wrt_int_state_x(N)%loc_int_state%ALL_DATA_I1D=wrt_int_state%ALL_DATA_I1D
        wrt_int_state_x(N)%loc_int_state%ALL_DATA_R1D=wrt_int_state%ALL_DATA_R1D
        wrt_int_state_x(N)%loc_int_state%FIELD_NAME=wrt_int_state%FIELD_NAME
        wrt_int_state_x(N)%loc_int_state%NAMES_I1D_STRING=wrt_int_state%NAMES_I1D_STRING
        wrt_int_state_x(N)%loc_int_state%NAMES_I2D_STRING=wrt_int_state%NAMES_I2D_STRING
        wrt_int_state_x(N)%loc_int_state%NAMES_R1D_STRING=wrt_int_state%NAMES_R1D_STRING
        wrt_int_state_x(N)%loc_int_state%NAMES_R2D_STRING=wrt_int_state%NAMES_R2D_STRING
        wrt_int_state_x(N)%loc_int_state%NAMES_LOG_STRING=wrt_int_state%NAMES_LOG_STRING
        wrt_int_state_x(N)%loc_int_state%ALL_DATA_LOG=wrt_int_state%ALL_DATA_LOG
!
        DO NN=0,LAST_FCST_TASK
          wrt_int_state_x(N)%loc_int_state%LOCAL_ISTART(NN)=wrt_int_state%LOCAL_ISTART(NN)
          wrt_int_state_x(N)%loc_int_state%LOCAL_IEND(NN)=wrt_int_state%LOCAL_IEND(NN)
          wrt_int_state_x(N)%loc_int_state%LOCAL_JSTART(NN)=wrt_int_state%LOCAL_JSTART(NN)
          wrt_int_state_x(N)%loc_int_state%LOCAL_JEND(NN)=wrt_int_state%LOCAL_JEND(NN)
        ENDDO
!
      ELSEIF(OUTPUT_FLAG=='Restart')THEN
!
        wrt_int_state_x(N)%loc_int_state%INPES=wrt_int_state%INPES
        wrt_int_state_x(N)%loc_int_state%JNPES=wrt_int_state%JNPES
        wrt_int_state_x(N)%loc_int_state%IHALO=wrt_int_state%IHALO
        wrt_int_state_x(N)%loc_int_state%JHALO=wrt_int_state%JHALO
        ALLOCATE(wrt_int_state_x(N)%loc_int_state%NUM_WORDS_SEND_I2D_RST)
        ALLOCATE(wrt_int_state_x(N)%loc_int_state%NUM_WORDS_SEND_R2D_RST)
        ALLOCATE(wrt_int_state_x(N)%loc_int_state%RST_FIELD_NAME(1:5000))
        ALLOCATE(wrt_int_state_x(N)%loc_int_state%RST_NAMES_I1D_STRING)
        ALLOCATE(wrt_int_state_x(N)%loc_int_state%RST_NAMES_I2D_STRING)
        ALLOCATE(wrt_int_state_x(N)%loc_int_state%RST_NAMES_R1D_STRING)
        ALLOCATE(wrt_int_state_x(N)%loc_int_state%RST_NAMES_R2D_STRING)
        ALLOCATE(wrt_int_state_x(N)%loc_int_state%RST_NAMES_LOG_STRING)
        ALLOCATE(wrt_int_state_x(N)%loc_int_state%RST_ALL_DATA_I2D(wrt_int_state%NUM_WORDS_SEND_I2D_RST),stat=ISTAT)
        ALLOCATE(wrt_int_state_x(N)%loc_int_state%RST_ALL_DATA_R2D(wrt_int_state%NUM_WORDS_SEND_R2D_RST),stat=ISTAT)
!
        wrt_int_state_x(N)%loc_int_state%NUM_WORDS_SEND_I2D_RST=wrt_int_state%NUM_WORDS_SEND_I2D_RST
        wrt_int_state_x(N)%loc_int_state%NUM_WORDS_SEND_R2D_RST=wrt_int_state%NUM_WORDS_SEND_R2D_RST
        wrt_int_state_x(N)%loc_int_state%RST_KOUNT_I1D=wrt_int_state%RST_KOUNT_I1D
        wrt_int_state_x(N)%loc_int_state%RST_KOUNT_I2D=wrt_int_state%RST_KOUNT_I2D
        wrt_int_state_x(N)%loc_int_state%RST_KOUNT_R1D=wrt_int_state%RST_KOUNT_R1D
        wrt_int_state_x(N)%loc_int_state%RST_KOUNT_R2D=wrt_int_state%RST_KOUNT_R2D
        wrt_int_state_x(N)%loc_int_state%RST_KOUNT_LOG=wrt_int_state%RST_KOUNT_LOG
        wrt_int_state_x(N)%loc_int_state%RST_LENGTH_DATA_I1D=wrt_int_state%RST_LENGTH_DATA_I1D
        wrt_int_state_x(N)%loc_int_state%RST_LENGTH_DATA_R1D=wrt_int_state%RST_LENGTH_DATA_R1D
!       wrt_int_state_x(N)%loc_int_state%RST_LENGTH_DATA_R2D=wrt_int_state%RST_LENGTH_DATA_R2D
        wrt_int_state_x(N)%loc_int_state%RST_LENGTH_SUM_I1D=wrt_int_state%RST_LENGTH_SUM_I1D
        wrt_int_state_x(N)%loc_int_state%RST_LENGTH_SUM_R1D=wrt_int_state%RST_LENGTH_SUM_R1D
!       wrt_int_state_x(N)%loc_int_state%RST_LENGTH_SUM_R2D=wrt_int_state%RST_LENGTH_SUM_R2D
        wrt_int_state_x(N)%loc_int_state%RST_LENGTH_SUM_LOG=wrt_int_state%RST_LENGTH_SUM_LOG
        wrt_int_state_x(N)%loc_int_state%RST_NCOUNT_FIELDS=wrt_int_state%RST_NCOUNT_FIELDS
        wrt_int_state_x(N)%loc_int_state%RST_ALL_DATA_I1D=wrt_int_state%RST_ALL_DATA_I1D
        wrt_int_state_x(N)%loc_int_state%RST_ALL_DATA_R1D=wrt_int_state%RST_ALL_DATA_R1D
        wrt_int_state_x(N)%loc_int_state%RST_FIELD_NAME=wrt_int_state%RST_FIELD_NAME
        wrt_int_state_x(N)%loc_int_state%RST_NAMES_I1D_STRING=wrt_int_state%RST_NAMES_I1D_STRING
        wrt_int_state_x(N)%loc_int_state%RST_NAMES_I2D_STRING=wrt_int_state%RST_NAMES_I2D_STRING
        wrt_int_state_x(N)%loc_int_state%RST_NAMES_R1D_STRING=wrt_int_state%RST_NAMES_R1D_STRING
        wrt_int_state_x(N)%loc_int_state%RST_NAMES_R2D_STRING=wrt_int_state%RST_NAMES_R2D_STRING
        wrt_int_state_x(N)%loc_int_state%RST_NAMES_LOG_STRING=wrt_int_state%RST_NAMES_LOG_STRING
        wrt_int_state_x(N)%loc_int_state%RST_ALL_DATA_LOG=wrt_int_state%RST_ALL_DATA_LOG
!
        DO NN=0,LAST_FCST_TASK
          wrt_int_state_x(N)%loc_int_state%LOCAL_ISTART(NN)=wrt_int_state%LOCAL_ISTART(NN)
          wrt_int_state_x(N)%loc_int_state%LOCAL_IEND(NN)=wrt_int_state%LOCAL_IEND(NN)
          wrt_int_state_x(N)%loc_int_state%LOCAL_JSTART(NN)=wrt_int_state%LOCAL_JSTART(NN)
          wrt_int_state_x(N)%loc_int_state%LOCAL_JEND(NN)=wrt_int_state%LOCAL_JEND(NN)
        ENDDO
!
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE FILL_GROUP_STATES
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE PRELIM_INFO_FOR_OUTPUT
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE OPEN_HST_FILE(WRT_INT_STATE)
!
!-----------------------------------------------------------------------
!***  Open a history disk file.
!-----------------------------------------------------------------------
!
!-----------------------
!*** Argument variables
!-----------------------
!
      TYPE(WRITE_INTERNAL_STATE),INTENT(INOUT) :: WRT_INT_STATE            !<-- The I/O component's internal state
!
!---------------------
!***  Local variables
!---------------------
!
      INTEGER :: FRAC_SEC,INT_SEC,N,RC
!
      LOGICAL :: OPENED
!
      CHARACTER(ESMF_MAXSTR) :: FILENAME
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Specifying 'DEFERRED' as the filename in the configure file
!***  means that we want to construct the output filename from
!***  HST_NAME_BASE (from the config file) appended with the
!***  forecast hour.
!-----------------------------------------------------------------------
!
      INT_SEC=INT(wrt_int_state%NFSECONDS)
      FRAC_SEC=NINT((wrt_int_state%NFSECONDS-INT_SEC)*100.)
      WRITE(FILENAME,100)TRIM(wrt_int_state%HST_NAME_BASE)//'_bin_'     &
                        ,wrt_int_state%NFHOURS,'h_'                     &
                        ,wrt_int_state%NFMINUTES,'m_'                   &
                        ,INT_SEC,'.',FRAC_SEC,'s'
  100 FORMAT(A,I4.4,A,I2.2,A,I2.2,A,I2.2,A)
!
!-----------------------------------------------------------------------
!***  Find an unopened unit number if one was not designated in
!***  the configure file.
!-----------------------------------------------------------------------
!
      DO N=51,99
        INQUIRE(N,opened=OPENED)
        IF(.NOT.OPENED)THEN
          wrt_int_state%IO_HST_UNIT=N
          EXIT
        ENDIF
      ENDDO
!
!-----------------------------------------------------------------------
!***  Open the file now.
!-----------------------------------------------------------------------
!
      OPEN(unit  =wrt_int_state%IO_HST_UNIT                             &
          ,file  =FILENAME                                              &
          ,status='REPLACE'                                             &
          ,access='SEQUENTIAL'                                          &
          ,form  ='UNFORMATTED'                                         &
          ,iostat=RC)
!
      IF(RC==0)THEN
        IF(wrt_int_state%PRINT_OUTPUT .OR. wrt_int_state%PRINT_ALL) THEN
          WRITE(0,*)' Opened IO_HST_UNIT=',wrt_int_state%IO_HST_UNIT,' for history'
          write(0,*)' iostat=',rc,' file=',trim(filename)
        ENDIF
      ELSE
        WRITE(0,*)' Failed to OPEN IO_HST_UNIT=',wrt_int_state%IO_HST_UNIT,' for history'
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE OPEN_HST_FILE
!
!-----------------------------------------------------------------------
!#######################################################################
!
!-----------------------------------------------------------------------
!
      SUBROUTINE OPEN_RST_FILE(WRT_INT_STATE)
!
!-----------------------------------------------------------------------
!***  Open a restart disk file.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument variables
!------------------------
!
      TYPE(WRITE_INTERNAL_STATE),INTENT(INOUT) :: WRT_INT_STATE            !<-- The I/O component's internal state
!
!---------------------
!***  Local variables
!---------------------
!
      INTEGER :: FRAC_SEC,INT_SEC,N,RC
!
      LOGICAL :: OPENED
!
      CHARACTER(ESMF_MAXSTR) :: FILENAME
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Specifying 'DEFERRED' as the filename in the configure file
!***  means that we want to construct the output filename from
!***  RST_NAME_BASE (from the config file) appended with the
!***  forecast hour.
!-----------------------------------------------------------------------
!
      INT_SEC=INT(wrt_int_state%NFSECONDS)
      FRAC_SEC=NINT((wrt_int_state%NFSECONDS-INT_SEC)*100.)
      WRITE(FILENAME,100)TRIM(wrt_int_state%RST_NAME_BASE)//'_bin_'     &
                        ,wrt_int_state%NFHOURS,'h_'                     &
                        ,wrt_int_state%NFMINUTES,'m_'                   &
                        ,INT_SEC,'.',FRAC_SEC,'s'
  100 FORMAT(A,I4.4,A,I2.2,A,I2.2,A,I2.2,A)
!
!-----------------------------------------------------------------------
!***  Find an unopened unit number if one was not designated in
!***  the configure file.
!-----------------------------------------------------------------------
!
      DO N=51,99
        INQUIRE(N,opened=OPENED)
        IF(.NOT.OPENED)THEN
          wrt_int_state%IO_RST_UNIT=N
          EXIT
        ENDIF
      ENDDO
!
!-----------------------------------------------------------------------
!***  Open the file now.
!-----------------------------------------------------------------------
!
      OPEN(unit  =wrt_int_state%IO_RST_UNIT                             &
          ,file  =FILENAME                                              &
          ,status='REPLACE'                                             &
          ,access='SEQUENTIAL'                                          &
          ,form  ='UNFORMATTED'                                         &
          ,iostat=RC)
!
      IF(RC==0) THEN
        IF(wrt_int_state%PRINT_OUTPUT .OR. wrt_int_state%PRINT_ALL) THEN
          WRITE(0,*)' Opened IO_RST_UNIT=',wrt_int_state%IO_RST_UNIT,' for restart'
          write(0,*)' iostat=',rc,' file=',trim(filename)
        ENDIF
      ELSE
        WRITE(0,*)' Failed to OPEN IO_RST_UNIT=',wrt_int_state%IO_RST_UNIT,' for restart'
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE OPEN_RST_FILE
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE WRITE_INIT(DOMAIN_GRID_COMP                            &
                           ,DOMAIN_INT_STATE                            &
                           ,DOMAIN_IMP_STATE                            &
                           ,CLOCK_DOMAIN                                &
                           ,MYPE)
! 
!-----------------------------------------------------------------------
!***  Execute the Initialize step of the Write components.
!-----------------------------------------------------------------------
!
!-----------------------
!*** Argument Variables
!-----------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT)         :: DOMAIN_GRID_COMP        !<-- The DOMAIN component
!
      TYPE(DOMAIN_INTERNAL_STATE),INTENT(INOUT) :: DOMAIN_INT_STATE        !<-- The DOMAIN Internal State
!
      TYPE(ESMF_State),INTENT(INOUT)            :: DOMAIN_IMP_STATE        !<-- The DOMAIN import state
!
      TYPE(ESMF_Clock),INTENT(INOUT)            :: CLOCK_DOMAIN            !<-- The DOMAIN Component's ESMF Clock
!
      INTEGER,INTENT(IN)                        :: MYPE
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: I,INPES,J,JNPES,LEAD_TASK,NUM_PES_FCST
!
      INTEGER(kind=KINT) :: RC,RC_INIT
!
      TYPE(ESMF_Config) :: CF                                              !<-- The config object
!
      TYPE(ESMF_VM) :: VM                                                  !<-- The ESMF virtual machine.
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Retrieve the config object CF from the DOMAIN component.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Get Config Object from DOMAIN Component in Write Init"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompGet(gridcomp=DOMAIN_GRID_COMP                   &  !<-- The DOMAIN component
                           ,config  =CF                                 &  !<-- The config object (~namelist)
                           ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  How many forecast tasks do we have?
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Write Init: Get INPES from Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The ESMF configure object
                                  ,value =INPES                         &  !<-- # of fcst tasks in I direction
                                  ,label ='inpes:'                      &  !<-- Give the value of this label to INPES
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Write Init: Get JNPES from Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The ESMF configure object
                                  ,value =JNPES                         &  !<-- # of fcst tasks in J direction
                                  ,label ='jnpes:'                      &  !<-- Give the value of this label to JNPES
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      NUM_PES_FCST=INPES*JNPES                                             !<-- Total number of forecast tasks
!
!-----------------------------------------------------------------------
!***  Transfer the rank of the lead task on each domain into the
!***  Write component import states since that component needs
!***  this information.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Write Init: Extract Lead Task on Domain from Domain Import"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state=DOMAIN_IMP_STATE                   &  !<-- The Domain component's import state
                              ,name ='Lead Task Domain'                 &  !<-- Name of the Attribute to extract
                              ,value=LEAD_TASK                          &  !<-- Global rank of lead task on this domain
                              ,defaultValue=0                           &  !<-- The default value
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Write Init: Insert Lead Task into Write Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=domain_int_state%IMP_STATE_WRITE   &  !<-- The Write component's import state
                              ,name ='Lead Task Domain'                 &  !<-- Name of the Attribute to extract
                              ,value=LEAD_TASK                          &  !<-- Global rank of lead task on this domain
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
       CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Execute the Initialize step for the Write components.
!***  These are the Initialize subroutines specified in the
!***  register routines called in ESMF_GridCompSetServices.
!-----------------------------------------------------------------------
!
      DO J=1,domain_int_state%WRITE_GROUPS
!
!!!!    CALL ESMF_VMBarrier(VM,rc=RC)    ! Insert barrier since fcst tasks are involved in each iteration of write groups
!
        DO I=1,NUM_PES_FCST+domain_int_state%WRITE_TASKS_PER_GROUP
          IF(MYPE==domain_int_state%PETLIST_WRITE(I,J))THEN                !<--  Forecast tasks plus the Write tasks in each write group
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Execute Initialize Step of Write Component"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_GridCompInitialize(gridcomp   =domain_int_state%WRITE_COMPS(J)   &  !<-- The Write gridded components
                                        ,importState=domain_int_state%IMP_STATE_WRITE  &  !<-- The Write import state
                                        ,exportState=domain_int_state%EXP_STATE_WRITE  &  !<-- The Write export state
                                        ,clock      =CLOCK_DOMAIN                      &  !<-- The ESMF clock of the DOMAIN component
                                        ,rc         =RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          ENDIF
        ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
!***  Some aspects of the output data do not change with forecast time
!***  for both the forecast and the quilt tasks.  Determine all such
!***  information now and save it.  Do this work once for history
!***  and once for restart output.
!-----------------------------------------------------------------------
!
      CALL PRELIM_INFO_FOR_OUTPUT('History'                             &
                                 ,NUM_PES_FCST                          &
                                 ,domain_int_state%WRITE_GROUPS         &
                                 ,domain_int_state%WRITE_COMPS          &
                                 ,domain_int_state%IMP_STATE_WRITE )
!
      CALL PRELIM_INFO_FOR_OUTPUT('Restart'                             &
                                 ,NUM_PES_FCST                          &
                                 ,domain_int_state%WRITE_GROUPS         &
                                 ,domain_int_state%WRITE_COMPS          &
                                 ,domain_int_state%IMP_STATE_WRITE )
!
!-----------------------------------------------------------------------
!***  Set the first Write group as the first one to act.
!-----------------------------------------------------------------------
!
      domain_int_state%WRITE_GROUP_READY_TO_GO=1
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE WRITE_INIT
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE SEND_UPDATED_ATTRIBUTES(OUTPUT_BUNDLE                  &
                                        ,ALL_DATA_INT_ATT               &
                                        ,ALL_DATA_REAL_ATT              &
                                        ,ALL_DATA_LOG_ATT               &
                                        ,LENGTH_INT_DATA                &
                                        ,LENGTH_REAL_DATA               &
                                        ,LENGTH_LOG_DATA                &
                                        ,MAX_GROUPS                     &
                                        ,WRITE_GROUP                    &
                                        ,INTERCOMM )
!
!-----------------------------------------------------------------------
!***  The lead forecast task sends the updated ESMF Attributes in
!***  the input Bundle to the lead quilt task.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      INTEGER(kind=KINT),INTENT(IN) :: INTERCOMM                        &  !<-- Intercommunicator between fcst/quilt tasks
                                      ,LENGTH_INT_DATA                  &  !<-- Length of integer Attribute datastring
                                      ,LENGTH_LOG_DATA                  &  !<-- Length of logical Attribute datastring
                                      ,LENGTH_REAL_DATA                 &  !<-- Length of real Attribute datastring
                                      ,MAX_GROUPS                       &  !<-- Max # of Write groups
                                      ,WRITE_GROUP                         !<-- The current Write group
!
      INTEGER(kind=KINT),DIMENSION(1:LENGTH_INT_DATA),INTENT(OUT) ::    &
                                                       ALL_DATA_INT_ATT    !<-- The integer Attributes in the output Bundle
!
      REAL(kind=KFPT),DIMENSION(1:LENGTH_REAL_DATA),INTENT(OUT) ::      &
                                                       ALL_DATA_REAL_ATT   !<-- The real Attributes in the output Bundle
!
      LOGICAL(kind=KLOG),DIMENSION(1:LENGTH_LOG_DATA),INTENT(OUT) ::    &
                                                       ALL_DATA_LOG_ATT    !<-- The logical Attributes in the output Bundle
!
      TYPE(ESMF_FieldBundle),INTENT(INOUT) :: OUTPUT_BUNDLE                !<-- The ESMF output data Bundle (history/restart)
!
!---------------------
!***  Local variables
!---------------------
!
      INTEGER(kind=KINT) :: KOUNT_INT_ATT,KOUNT_LOG_ATT,KOUNT_REAL_ATT  &
                           ,L,LENGTH,LENGTH_SUM_INT_ATT                 &
                           ,LENGTH_SUM_LOG_ATT,LENGTH_SUM_REAL_ATT      &
                           ,N,NUM_ATTRIB
!
      INTEGER(kind=KINT) :: IERR,RC,RC_ATT
!
      INTEGER(kind=KINT),DIMENSION(MPI_STATUS_SIZE) :: JSTAT
!
      INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE :: WORK_ARRAY_INT_ATT
!
      INTEGER(kind=KINT),DIMENSION(:),POINTER :: HANDLE_INT_ATT         &
                                                ,HANDLE_LOG_ATT         &
                                                ,HANDLE_REAL_ATT
!
      INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE,TARGET,SAVE ::        &
                                                    HANDLE_INT_ATT_HST  &
                                                   ,HANDLE_INT_ATT_RST  &
                                                   ,HANDLE_LOG_ATT_HST  &
                                                   ,HANDLE_LOG_ATT_RST  &
                                                   ,HANDLE_REAL_ATT_HST &
                                                   ,HANDLE_REAL_ATT_RST
!
      REAL(kind=KFPT),DIMENSION(:),ALLOCATABLE :: WORK_ARRAY_REAL_ATT
!
      LOGICAL(kind=KLOG) :: WORK_LOGICAL
!
      CHARACTER(len=14) :: BUNDLE_NAME
!
      CHARACTER(ESMF_MAXSTR) :: ATTRIB_NAME
!
      TYPE(ESMF_TypeKind_Flag) :: DATATYPE
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Get Name of the Output Bundle"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_FieldBundleGet(FIELDBUNDLE   =OUTPUT_BUNDLE             &
                              ,name          =BUNDLE_NAME               &
                              ,rc            =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_ATT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Get Attribute Count from Output Bundle"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(FIELDBUNDLE=OUTPUT_BUNDLE                  &  !<-- The write component's history data Bundle
                            ,count      =NUM_ATTRIB                     &  !<-- # of Attributes in the history data Bundle
                            ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_ATT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      KOUNT_INT_ATT=0                                                      !<-- Initialize counter of integer Attributes
      KOUNT_REAL_ATT=0                                                     !<-- Initialize counter of real Attributes
      KOUNT_LOG_ATT=0                                                      !<-- Initialize counter of logical Attributes
!
      LENGTH_SUM_INT_ATT=0                                                 !<-- Initialize length of the integer Attribute datastring
      LENGTH_SUM_REAL_ATT=0                                                !<-- Initialize length of the real Attribute datastring
      LENGTH_SUM_LOG_ATT=0                                                 !<-- Initialize length of the logical Attribute datastring
!
!-----------------------------------------------------------------------
!***  Allocate the ISend handles if not done already.
!-----------------------------------------------------------------------
!
      IF(BUNDLE_NAME=='History Bundle')THEN
!
        IF(.NOT.ALLOCATED(HANDLE_INT_ATT_HST))THEN
!
          ALLOCATE(HANDLE_INT_ATT_HST(1:MAX_GROUPS))
          ALLOCATE(HANDLE_REAL_ATT_HST(1:MAX_GROUPS))
          ALLOCATE(HANDLE_LOG_ATT_HST(1:MAX_GROUPS))
!
          DO N=1,MAX_GROUPS
            HANDLE_INT_ATT_HST(N)=MPI_REQUEST_NULL
            HANDLE_REAL_ATT_HST(N)=MPI_REQUEST_NULL
            HANDLE_LOG_ATT_HST(N)=MPI_REQUEST_NULL
          ENDDO
!
        ENDIF
!
        HANDLE_INT_ATT=>HANDLE_INT_ATT_HST
        HANDLE_REAL_ATT=>HANDLE_REAL_ATT_HST
        HANDLE_LOG_ATT=>HANDLE_LOG_ATT_HST
!
      ELSEIF(BUNDLE_NAME=='Restart Bundle')THEN
!
        IF(.NOT.ALLOCATED(HANDLE_INT_ATT_RST))THEN
!
          ALLOCATE(HANDLE_INT_ATT_RST(1:MAX_GROUPS))
          ALLOCATE(HANDLE_REAL_ATT_RST(1:MAX_GROUPS))
          ALLOCATE(HANDLE_LOG_ATT_RST(1:MAX_GROUPS))
!
          DO N=1,MAX_GROUPS
            HANDLE_INT_ATT_RST(N)=MPI_REQUEST_NULL
            HANDLE_REAL_ATT_RST(N)=MPI_REQUEST_NULL
            HANDLE_LOG_ATT_RST(N)=MPI_REQUEST_NULL
          ENDDO
!
        ENDIF
!
        HANDLE_INT_ATT=>HANDLE_INT_ATT_RST
        HANDLE_REAL_ATT=>HANDLE_REAL_ATT_RST
        HANDLE_LOG_ATT=>HANDLE_LOG_ATT_RST
!
      ENDIF
!
!-----------------------------------------------------------------------
!
      attrib_loop: DO N=1,NUM_ATTRIB
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Get Attribute Names, Datatypes, and Lengths"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(FIELDBUNDLE   =OUTPUT_BUNDLE             &  !<-- The history data Bundle
                              ,attributeIndex=N                         &  !<-- Index of each Attribute
                              ,name          =ATTRIB_NAME               &  !<-- Each Attribute's name
                              ,typekind      =DATATYPE                  &  !<-- Each Attribute's ESMF Datatype
                              ,itemCount     =LENGTH                    &  !<-- Each Attribute's length
                              ,rc            =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_ATT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!                 -- Scalar and 1-D Integer Output Data --
!-----------------------------------------------------------------------
!
        IF(DATATYPE==ESMF_TYPEKIND_I4)THEN                                 !<-- Extract integer data with rank <2
!
          ALLOCATE(WORK_ARRAY_INT_ATT(LENGTH),stat=RC)                     !<-- Allocate array to hold integer Attribute N
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Get Scalar/1-D Integer Data from Output Bundle"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeGet(FIELDBUNDLE=OUTPUT_BUNDLE              &  !<-- The history data Bundle
                                ,name       =ATTRIB_NAME                &  !<-- Name of the Attribute to extract
                                ,itemCount  =LENGTH                     &  !<-- Length of Attribute
                                ,valueList  =WORK_ARRAY_INT_ATT         &  !<-- Place the Attribute here
                                ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_ATT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          KOUNT_INT_ATT=KOUNT_INT_ATT+1                                    !<-- Count # of integer Attributes
!
          DO L=1,LENGTH
            ALL_DATA_INT_ATT(LENGTH_SUM_INT_ATT+L)=WORK_ARRAY_INT_ATT(L)   !<-- String together the integer Attributes
          ENDDO
!
          LENGTH_SUM_INT_ATT=LENGTH_SUM_INT_ATT+LENGTH                     !<-- Total word sum of integer Attributesd
          DEALLOCATE(WORK_ARRAY_INT_ATT)
!
!-----------------------------------------------------------------------
!                 -- Scalar and 1-D Real Output Data --
!-----------------------------------------------------------------------
!
        ELSEIF(DATATYPE==ESMF_TYPEKIND_R4)THEN                             !<-- Extract real data with rank <2
!
          ALLOCATE(WORK_ARRAY_REAL_ATT(LENGTH),stat=RC)                   !<-- Allocate array to hold real Attribute N
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Get Scalar/1-D Real Data from Output Bundle"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeGet(FIELDBUNDLE=OUTPUT_BUNDLE            &  !<-- The history data Bundle
                                ,name       =ATTRIB_NAME              &  !<-- Name of the Attribute to extract
                                ,itemCount  =LENGTH                   &  !<-- Length of Attribute
                                ,valueList  =WORK_ARRAY_REAL_ATT      &  !<-- Place the Attribute here
                                ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_ATT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          KOUNT_REAL_ATT=KOUNT_REAL_ATT+1                                 !<-- Count # of real Attributes
!
          DO L=1,LENGTH
            ALL_DATA_REAL_ATT(LENGTH_SUM_REAL_ATT+L)=WORK_ARRAY_REAL_ATT(L) !<-- String together the real Attributes
          ENDDO
!
          LENGTH_SUM_REAL_ATT=LENGTH_SUM_REAL_ATT+LENGTH                  !<-- Total word sum of real Attributes
          DEALLOCATE(WORK_ARRAY_REAL_ATT)
!
!-----------------------------------------------------------------------
!                       -- Logical Output Data --
!-----------------------------------------------------------------------
!
!!!     ELSEIF(DATATYPE==ESMF_TYPEKIND_LOGICAL)THEN                       !<-- Extract logical data
        ELSE                                                              !<-- Extract logical data
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Get Scalar/1-D Logical Data from Output Bundle"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeGet(FIELDBUNDLE=OUTPUT_BUNDLE            &  !<-- The history data Bundle
                                ,name       =ATTRIB_NAME              &  !<-- Name of the Attribute to extract
                                ,value      =WORK_LOGICAL             &  !<-- Place the Attribute here
                                ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_ATT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          KOUNT_LOG_ATT=KOUNT_LOG_ATT+1                                   !<-- Count # of logical Attributes
!
          ALL_DATA_LOG_ATT(LENGTH_SUM_LOG_ATT+1)=WORK_LOGICAL             !<-- String together the logical Attributes
!
          LENGTH_SUM_LOG_ATT=LENGTH_SUM_LOG_ATT+1                         !<-- Total word sum of logical Attributes
!
!-----------------------------------------------------------------------
!
        ENDIF
!
!-----------------------------------------------------------------------
!
      ENDDO attrib_loop
!
!-----------------------------------------------------------------------
!***  Lead fcst task sends the integer Attributes to the lead
!***  quilt task.
!-----------------------------------------------------------------------
!
      IF(LENGTH_SUM_INT_ATT>0)THEN
!
        CALL MPI_WAIT(HANDLE_INT_ATT(WRITE_GROUP),JSTAT,IERR)
!
        CALL MPI_ISSEND(ALL_DATA_INT_ATT                                &  !<-- String of integer Attribute output data
                       ,LENGTH_SUM_INT_ATT                              &  !<-- # of words in the data string
                       ,MPI_INTEGER                                     &  !<-- The datatype
                       ,0                                               &  !<-- Send to the lead Write task in the group
                       ,WRITE_GROUP                                     &  !<-- Use the Write group as the tag
                       ,INTERCOMM                                       &  !<-- The MPI intercommunicator between fcst and quilt tasks
                       ,HANDLE_INT_ATT(WRITE_GROUP)                     &  !<-- MPI communication request handle
                       ,IERR )
!
      ENDIF
!
!-----------------------------------------------------------------------
!***  Lead fcst task sends the real Attributes to the lead quilt task.
!-----------------------------------------------------------------------
!
      IF(LENGTH_SUM_REAL_ATT>0)THEN
!
        CALL MPI_WAIT(HANDLE_REAL_ATT(WRITE_GROUP),JSTAT,IERR)
!
        CALL MPI_ISSEND(ALL_DATA_REAL_ATT                               &  !<-- String of real Attribute output data
                       ,LENGTH_SUM_REAL_ATT                             &  !<-- # of words in the data string
                       ,MPI_REAL                                        &  !<-- The datatype
                       ,0                                               &  !<-- Send to the lead Write task in the group
                       ,WRITE_GROUP                                     &  !<-- Use the Write group as the tag
                       ,INTERCOMM                                       &  !<-- The MPI intercommunicator between fcst and quilt tasks
                       ,HANDLE_REAL_ATT(WRITE_GROUP)                    &  !<-- MPI communication request handle
                       ,IERR )
!
      ENDIF
!
!-----------------------------------------------------------------------
!***  Lead fcst task sends the logical Attributes to the 
!***  lead quilt task.
!-----------------------------------------------------------------------
!
      IF(LENGTH_SUM_LOG_ATT>0)THEN
!
        CALL MPI_WAIT(HANDLE_LOG_ATT(WRITE_GROUP),JSTAT,IERR)
!
        CALL MPI_ISSEND(ALL_DATA_LOG_ATT                                &  !<-- String of logical Attribute output data
                       ,LENGTH_SUM_LOG_ATT                              &  !<-- # of words in the data string
                       ,MPI_LOGICAL                                     &  !<-- The datatype
                       ,0                                               &  !<-- Send to the lead Write task in the group
                       ,WRITE_GROUP                                     &  !<-- Use the Write group as the tag
                       ,INTERCOMM                                       &  !<-- The MPI intercommunicator between fcst and quilt tasks
                       ,HANDLE_LOG_ATT(WRITE_GROUP)                     &  !<-- MPI communication request handle
                       ,IERR )
!
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE SEND_UPDATED_ATTRIBUTES
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE WRITE_ASYNC(DOMAIN_GRID_COMP                           &
                            ,DOMAIN_INT_STATE                           &
                            ,CLOCK_DOMAIN                               &
                            ,MYPE                                       &
                            ,CWRT)
!
!-----------------------------------------------------------------------
!***  Write out a history file using the asynchronous quilting.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: DOMAIN_GRID_COMP               !<-- The DOMAIN component
!
      TYPE(DOMAIN_INTERNAL_STATE),INTENT(INOUT) :: DOMAIN_INT_STATE       !<-- The DOMAIN Internal State
!
      TYPE(ESMF_Clock),INTENT(INOUT) :: CLOCK_DOMAIN                      !<-- The DOMAIN Component's ESMF Clock
!
      CHARACTER(ESMF_MAXSTR),INTENT(IN) :: CWRT                           !<-- Restart/History label
!
!---------------------
!***  Local Variables
!---------------------
!
      TYPE(ESMF_Config) :: CF                                             !<-- The configure object (~namelist)
!
      TYPE(ESMF_Time)   :: CURRTIME                                       !<-- The current forecast time (ESMF)
!
      INTEGER,INTENT(IN) :: MYPE
!
      INTEGER :: YY,MM,DD,H,M,S                                           !<-- Year, Month, Day, Hour, Minute, Second (integer)
!
      INTEGER :: I,INPES,JNPES,N_GROUP,NUM_PES_FCST                     &
                ,WRITE_GROUPS,WRITE_TASKS_PER_GROUP                     &
                ,RC,RC_ASYNC
!
      LOGICAL :: PRINT_OUTPUT, PRINT_ALL                                  !<-- Prints to err file flags
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Check whether this is history or restart call.
!-----------------------------------------------------------------------
!
      IF(CWRT=='History') TIME_FOR_HISTORY = .TRUE.
      IF(CWRT=='Restart') TIME_FOR_RESTART = .TRUE.
!     write(0,*)' enter WRITE_ASYNC time_for_history=',time_for_history,' time_for_restart=',time_for_restart
!
!-----------------------------------------------------------------------
!***  Extract the configure object in order to know the number
!***  of forecast tasks, Write groups, and write tasks per group.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="WRITE_ASYNC: Extract Config Object from DOMAIN Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompGet(gridcomp=DOMAIN_GRID_COMP                   &  !<-- The DOMAIN component
                           ,config  =CF                                 &  !<-- The config object (~namelist)
                           ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_ASYNC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="WRITE_ASYNC: Get General Info from Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
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
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure file
                                  ,value =WRITE_GROUPS                  &  !<-- Number of write groups
                                  ,label ='write_groups:'               &
                                  ,rc    =RC)
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure file
                                  ,value =WRITE_TASKS_PER_GROUP         &  !<-- Number of write tasks per group
                                  ,label ='write_tasks_per_group:'      &
                                  ,rc    =RC)
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure file
                                  ,value =PRINT_OUTPUT                  &  !<-- Print output flag
                                  ,label ='print_output:'               &
                                  ,rc    =RC)
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure file
                                  ,value =PRINT_ALL                     &  !<-- Print output flag
                                  ,label ='print_all:'                  &
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_ASYNC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      NUM_PES_FCST=INPES*JNPES                                             !<-- Number of forecast tasks
!
!-----------------------------------------------------------------------
!***  What is the current forecast time?
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="WRITE_ASYNC: Get Current Time from DOMAIN Clock"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ClockGet(clock   =CLOCK_DOMAIN                          &  !<-- The DOMAIN component's ESMF Clock
                        ,currTime=CURRTIME                              &  !<-- The current forecast time (ESMF)
                        ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_ASYNC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="WRITE_ASYNC: Convert ESMF Time to Real Time"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_TimeGet (time=CURRTIME                                  &  !<-- The current forecast time (ESMF)
                        ,yy  =YY                                        &  !<-- The current year (integer)
                        ,mm  =MM                                        &  !<-- The current month (integer)
                        ,dd  =DD                                        &  !<-- The current day (integer)
                        ,h   =H                                         &  !<-- The current hour (integer)
                        ,m   =M                                         &  !<-- The current minute (integer)
                        ,s   =S                                         &  !<-- The current second (integer)
                        ,rc  =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_ASYNC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  The export state of the SOLVER component lies within the
!***  internal state of the DOMAIN component and holds the
!***  import state of the Write component.
!***  Extract that Write component's import state since we are 
!***  about to execute the Run step of the Write component.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="WRITE_ASYNC: Extract Write Import State from Dyn Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateGet(state      =domain_int_state%EXP_STATE_SOLVER  &  !<-- The Solver component's export state
                        ,itemName   ="Write Import State"               &  !<-- Name of state to be extracted
                        ,nestedState=domain_int_state%IMP_STATE_WRITE   &  !<-- The extracted state
                        ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_ASYNC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!!!   CALL ESMF_VMBarrier(vm=VM,rc=RC)
!
!-----------------------------------------------------------------------
!***  All forecast tasks plus those write tasks in the appropriate
!***  Write group execute the Run step of a Write component.
!-----------------------------------------------------------------------
!
      N_GROUP=domain_int_state%WRITE_GROUP_READY_TO_GO                       !<-- The active write group
!     write(0,*)' WRITE_ASYNC calling WRITE_RUN for write group ',n_group &
!              ,' num_pes_fcst=',num_pes_fcst,' write_tasks_per_group=',write_tasks_per_group
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="WRITE_ASYNC: Execute Run Step of Write Components" 
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      DO I=1,NUM_PES_FCST+WRITE_TASKS_PER_GROUP
        IF(MYPE==domain_int_state%PETLIST_WRITE(I,N_GROUP))THEN
          CALL ESMF_GridCompRun(gridcomp=domain_int_state%WRITE_COMPS(N_GROUP) &  !<-- The write gridded component
                               ,importState=domain_int_state%IMP_STATE_WRITE   &  !<-- Its import state
                               ,exportState=domain_int_state%EXP_STATE_WRITE   &  !<-- Its export state
                               ,clock      =CLOCK_DOMAIN                       &  !<-- The DOMAIN Clock
                               ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_ASYNC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          IF(I==NUM_PES_FCST+1)THEN                                          !<-- The first write task tells us the history output time
            IF(PRINT_OUTPUT .OR. PRINT_ALL) &
            WRITE(0,101)TRIM(CWRT),YY,MM,DD,H,M,S
  101       FORMAT(' Wrote ',A7,' File at ',I4.4,'_',I2.2,'_',I2.2,'_',I2.2,':',I2.2,':',I2.2)
          ENDIF
!
        ENDIF
!
      ENDDO
!
!-----------------------------------------------------------------------
!***  Prepare to use the next write group at the next output time.
!***  Return to the 1st group if we have cycled through all of them.
!-----------------------------------------------------------------------
!
      IF(domain_int_state%WRITE_GROUP_READY_TO_GO==WRITE_GROUPS)THEN
        domain_int_state%WRITE_GROUP_READY_TO_GO=1
      ELSE
        domain_int_state%WRITE_GROUP_READY_TO_GO=domain_int_state%WRITE_GROUP_READY_TO_GO+1
      ENDIF
!
!-----------------------------------------------------------------------
!***  Restore TIME_FOR_HISTORY and TIME_FOR_RESTART to false.
!-----------------------------------------------------------------------
!
      TIME_FOR_HISTORY = .FALSE.
      TIME_FOR_RESTART = .FALSE.
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE WRITE_ASYNC
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE WRITE_RUNHISTORY_OPEN(WRT_INT_STATE                    &
                                      ,IYEAR_FCST                       &
                                      ,IMONTH_FCST                      &
                                      ,IDAY_FCST                        &
                                      ,IHOUR_FCST                       &
                                      ,IMINUTE_FCST                     &
                                      ,SECOND_FCST                      &
                                      ,NF_HOURS                         &
                                      ,NF_MINUTES                       &
                                      ,NF_SECONDS                       &
                                      ,HST_FIRST                        &
                                      ,LEAD_WRITE_TASK )
!
!-----------------------------------------------------------------------
!***  Write out a binary run history file.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument variables
!------------------------
!
      TYPE(WRITE_INTERNAL_STATE),INTENT(INOUT) :: WRT_INT_STATE            !<-- The Write component's internal state
!
      INTEGER,INTENT(IN) :: IYEAR_FCST                                  &
                           ,IMONTH_FCST                                 &
                           ,IDAY_FCST                                   &
                           ,IHOUR_FCST                                  &
                           ,IMINUTE_FCST                                &
                           ,NF_HOURS                                    &
                           ,NF_MINUTES                                  &
                           ,LEAD_WRITE_TASK
!
      REAL,INTENT(IN) :: NF_SECONDS,SECOND_FCST
!
      LOGICAL,INTENT(IN) :: HST_FIRST
!
!---------------------
!***  Local variables
!---------------------
!
      INTEGER :: N,N1,N2,NPOSN_1,NPOSN_2,LENGTH
      INTEGER :: NFIELD,RC
      CHARACTER(ESMF_MAXSTR)                :: NAME
      INTEGER,DIMENSION(:),POINTER          :: WORK_ARRAY_I1D
      REAL(4),DIMENSION(:),POINTER          :: WORK_ARRAY_R1D
      LOGICAL                               :: WRITE_LOGICAL
      LOGICAL                               :: WORK_LOGICAL

!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Open the history file and write the current forecast time
!***  and elapsed time.
!-----------------------------------------------------------------------
!
      CALL OPEN_HST_FILE(WRT_INT_STATE)
!
      WRITE(wrt_int_state%IO_HST_UNIT,iostat=RC)IYEAR_FCST
      WRITE(wrt_int_state%IO_HST_UNIT,iostat=RC)IMONTH_FCST
      WRITE(wrt_int_state%IO_HST_UNIT,iostat=RC)IDAY_FCST
      WRITE(wrt_int_state%IO_HST_UNIT,iostat=RC)IHOUR_FCST
      WRITE(wrt_int_state%IO_HST_UNIT,iostat=RC)IMINUTE_FCST
      WRITE(wrt_int_state%IO_HST_UNIT,iostat=RC)SECOND_FCST
      WRITE(wrt_int_state%IO_HST_UNIT,iostat=RC)NF_HOURS
      WRITE(wrt_int_state%IO_HST_UNIT,iostat=RC)NF_MINUTES
      WRITE(wrt_int_state%IO_HST_UNIT,iostat=RC)NF_SECONDS
!
      IF(HST_FIRST .AND. (wrt_int_state%PRINT_OUTPUT .OR. wrt_int_state%PRINT_ALL) )THEN
        WRITE(0,*)' Wrote IYEAR_FCST to history file unit ',wrt_int_state%IO_HST_UNIT
        WRITE(0,*)' Wrote IMONTH_FCST to history file unit ',wrt_int_state%IO_HST_UNIT
        WRITE(0,*)' Wrote IDAY_FCST to history file unit ',wrt_int_state%IO_HST_UNIT
        WRITE(0,*)' Wrote IHOUR_FCST to history file unit ',wrt_int_state%IO_HST_UNIT
        WRITE(0,*)' Wrote IMINUTE_FCST to history file unit ',wrt_int_state%IO_HST_UNIT
        WRITE(0,*)' Wrote SECOND_FCST to history file unit ',wrt_int_state%IO_HST_UNIT
        WRITE(0,*)' Wrote NF_HOURS to history file unit ',wrt_int_state%IO_HST_UNIT
        WRITE(0,*)' Wrote NF_MINUTES to history file unit ',wrt_int_state%IO_HST_UNIT
        WRITE(0,*)' Wrote NF_SECONDS to history file unit ',wrt_int_state%IO_HST_UNIT
      ENDIF
!
!-----------------------------------------------------------------------
!***  Integer scalar/1-D history variables
!-----------------------------------------------------------------------
!
      N2=0                                                                !<-- Word counter for full string of integer scalar/1D data
!
      DO N=1,wrt_int_state%KOUNT_I1D(1)                                   !<-- Loop through all scalar/1D integer data
!
        NPOSN_1=(N-1)*ESMF_MAXSTR+1
        NPOSN_2=N*ESMF_MAXSTR
        NAME=wrt_int_state%NAMES_I1D_STRING(NPOSN_1:NPOSN_2)              !<-- The variable's name
        LENGTH=wrt_int_state%LENGTH_DATA_I1D(N)                           !<-- The variable's length in words
        ALLOCATE(WORK_ARRAY_I1D(LENGTH),stat=RC)
!
        DO N1=1,LENGTH
          N2=N2+1
          WORK_ARRAY_I1D(N1)=wrt_int_state%ALL_DATA_I1D(N2)               !<-- Extract the individual data from the data string
        ENDDO
!
        WRITE(wrt_int_state%IO_HST_UNIT,iostat=RC)WORK_ARRAY_I1D          !<-- Write out the data
!
        IF(HST_FIRST .AND. (wrt_int_state%PRINT_OUTPUT .OR. wrt_int_state%PRINT_ALL) )THEN
          WRITE(0,*)'Wrote ',TRIM(NAME),' to history file unit ',wrt_int_state%IO_HST_UNIT
        ENDIF
!
        DEALLOCATE(WORK_ARRAY_I1D)
!
      ENDDO
!
!-----------------------------------------------------------------------
!***  Real scalar/1-D history variables
!-----------------------------------------------------------------------
!
      N2=0                                                                !<-- Word counter for full string of real scalar/1D data
!
      DO N=1,wrt_int_state%KOUNT_R1D(1)                                   !<-- Loop through all scalar/1D real data
!
        NPOSN_1=(N-1)*ESMF_MAXSTR+1
        NPOSN_2=N*ESMF_MAXSTR
        NAME=wrt_int_state%NAMES_R1D_STRING(NPOSN_1:NPOSN_2)              !<-- The variable's name
        LENGTH=wrt_int_state%LENGTH_DATA_R1D(N)                           !<-- The variable's length
        ALLOCATE(WORK_ARRAY_R1D(LENGTH),stat=RC)
!
        DO N1=1,LENGTH
          N2=N2+1
          WORK_ARRAY_R1D(N1)=wrt_int_state%ALL_DATA_R1D(N2)               !<-- Extract the individual data from the data string
        ENDDO
!
        WRITE(wrt_int_state%IO_HST_UNIT,iostat=RC)WORK_ARRAY_R1D          !<-- Write out the data
!
        IF(HST_FIRST .AND. (wrt_int_state%PRINT_OUTPUT .OR. wrt_int_state%PRINT_ALL) )THEN
          WRITE(0,*)'Wrote ',TRIM(NAME),' to history file unit ',wrt_int_state%IO_HST_UNIT
        ENDIF
!
        DEALLOCATE(WORK_ARRAY_R1D)
!
      ENDDO
!
!-----------------------------------------------------------------------
!***  Logical history variables
!-----------------------------------------------------------------------
!
      N2=0                                                                !<-- Counter for full string of logical data
!
      DO N=1,wrt_int_state%KOUNT_LOG(1)                                   !<-- Loop through all logical data
!
        NPOSN_1=(N-1)*ESMF_MAXSTR+1
        NPOSN_2=N*ESMF_MAXSTR
        NAME=wrt_int_state%NAMES_LOG_STRING(NPOSN_1:NPOSN_2)              !<-- The variable's name
!
        N2=N2+1

        WORK_LOGICAL = wrt_int_state%ALL_DATA_LOG(N2)
        WRITE_LOGICAL=WORK_LOGICAL                                        !<-- Convert from ESMF_Logical to F90 logical
!
        WRITE(wrt_int_state%IO_HST_UNIT,iostat=RC)WRITE_LOGICAL           !<-- Write out the data
!
        IF(HST_FIRST .AND. (wrt_int_state%PRINT_OUTPUT .OR. wrt_int_state%PRINT_ALL) )THEN
          WRITE(0,*)'Wrote ',TRIM(NAME),' to history file unit ',wrt_int_state%IO_HST_UNIT
        ENDIF
!
      ENDDO
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE WRITE_RUNHISTORY_OPEN
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE WRITE_NEMSIO_RUNHISTORY_OPEN(WRT_INT_STATE             &
                                             ,NEMSIOFILE                &
                                             ,IYEAR_FCST                &
                                             ,IMONTH_FCST               &
                                             ,IDAY_FCST                 &
                                             ,IHOUR_FCST                &
                                             ,IMINUTE_FCST              &
                                             ,SECOND_FCST               &
                                             ,NF_HOURS                  &
                                             ,NF_MINUTES                &
                                             ,NF_SECONDS                &
                                             ,DIM1,DIM2,NFRAME,GLOBAL   &
                                             ,LEAD_WRITE_TASK           &
                                             ,ID_DOMAIN)
!
!-----------------------------------------------------------------------
!***  Write out a NEMSIO binary run history file.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument variables
!------------------------
!
      TYPE(WRITE_INTERNAL_STATE),INTENT(INOUT) :: WRT_INT_STATE             !<-- The Write component's internal state
!
      TYPE(NEMSIO_GFILE),INTENT(INOUT)         :: NEMSIOFILE                !<-- The nemsio file handler
!
      INTEGER,INTENT(IN)  :: IYEAR_FCST                                 &
                            ,IMONTH_FCST                                &
                            ,IDAY_FCST                                  &
                            ,IHOUR_FCST                                 &
                            ,IMINUTE_FCST                               &
                            ,ID_DOMAIN                                  &
                            ,NF_HOURS                                   &
                            ,NF_MINUTES                                 &
                            ,LEAD_WRITE_TASK

      INTEGER,INTENT(OUT) :: DIM1,DIM2,NFRAME
      LOGICAL,INTENT(OUT) :: GLOBAL
!
      REAL,INTENT(IN)     :: NF_SECONDS                                 &
                            ,SECOND_FCST
!
!---------------------
!***  Local variables
!---------------------
!
      INTEGER :: I,J,N,N1,N2,NPOSN_1,NPOSN_2,LENGTH,MAXLENGTH
!
      INTEGER :: FIELDSIZE,IM,JM,LM,IDATE(7),FCSTDATE(7)                &
                ,INDX_2D,IRET,IND1,IND2,IND3,IND4,CNT                   &
 		,INI1,INI2,INI3                                         &
                ,N2ISCALAR,N2IARY,N2RSCALAR,N2RARY,N2LSCALAR            &
                ,NMETA,NSOIL,TLMETA,VLEV
!
      INTEGER :: FRAC_SEC,INT_SEC,NFIELD,RC
!
      INTEGER,DIMENSION(:),POINTER :: ARYILEN                           &
                                     ,ARYRLEN                           &
                                     ,RECLEV                            &
                                     ,VARIVAL
!
      INTEGER,DIMENSION(:,:),POINTER :: ARYIVAL
!
      REAL(4) :: DEGRAD,DXCTL,DYCTL,TPH0D,TLM0D
!
      REAL(4),DIMENSION(:),POINTER :: DX,DY,DXH                         &
                                     ,GLAT1D,GLON1D
!
      REAL(4),DIMENSION(:,:,:),POINTER :: VCOORD
!
      REAL(KIND=KFPT),DIMENSION(:)  ,POINTER :: VARRVAL
      REAL(KIND=KFPT),DIMENSION(:,:),POINTER :: ARYRVAL
!
      LOGICAL,DIMENSION(:),POINTER :: VARLVAL
!
      CHARACTER(6)  :: MODEL_LEVEL
      CHARACTER(16) :: VLEVTYP,FILE_ENDIAN
!
      CHARACTER(16),DIMENSION(:),POINTER :: ARYINAME                    &
                                           ,ARYRNAME                    &
                                           ,RECNAME                     &
                                           ,VARINAME                    &
                                           ,VARRNAME                    &
                                           ,VARLNAME
!
      CHARACTER(16),DIMENSION(:),POINTER :: RECLEVTYP
!
      CHARACTER(ESMF_MAXSTR) :: NAME,FILENAME
!
      LOGICAL            :: WORK_LOGICAL
!
      INTEGER :: NDYH=0,NDXH=0,NPT=0,NPDTOP=0,NREC=0                    &
                ,NSG1=0,NSG2=0,NSGML1=0,NSGML2=0
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      FCSTDATE(1)=IYEAR_FCST
      FCSTDATE(2)=IMONTH_FCST
      FCSTDATE(3)=IDAY_FCST
      FCSTDATE(4)=IHOUR_FCST
      FCSTDATE(5)=IMINUTE_FCST
      FCSTDATE(6)=NINT(SECOND_FCST*100.)
      FCSTDATE(7)=100
!
!-----------------------------------------------------------------------
!***  Integer scalar/1-D history variables
!-----------------------------------------------------------------------
!
!-------------------------------------------------------------
!*** Find out the total number of int scalars and int arrays.
!-------------------------------------------------------------
!
      N2ISCALAR=0
      N2IARY=0
      MAXLENGTH=1
!
      DO N=1,wrt_int_state%KOUNT_I1D(1)                                    !<-- Loop through all scalar/1D integer data
        LENGTH=wrt_int_state%LENGTH_DATA_I1D(N)
!
        IF(LENGTH==1)THEN
          N2ISCALAR=N2ISCALAR+1
        ELSE
          N2IARY=N2IARY+1
          MAXLENGTH=MAX(LENGTH,MAXLENGTH)
        ENDIF
!
      ENDDO
!
      N2IARY=N2IARY+1
      MAXLENGTH=MAX(MAXLENGTH,7)
      ALLOCATE(VARINAME(N2ISCALAR),VARIVAL(N2ISCALAR))
      ALLOCATE(ARYINAME(N2IARY),ARYILEN(N2IARY),ARYIVAL(MAXLENGTH,N2IARY))
!
!---------------------------------------
!***  Set value to AVRIVAL and ARYIVAL.
!---------------------------------------
!
      N2=0                                                                 !<-- Word counter for full string of integer scalar/1D data
      N2ISCALAR=0
      N2IARY=0
      IDATE=0
!
      DO N=1,wrt_int_state%KOUNT_I1D(1)                                    !<-- Loop through all scalar/1D integer data
!
        NPOSN_1=(N-1)*ESMF_MAXSTR+1
        NPOSN_2=N*ESMF_MAXSTR
        NAME=wrt_int_state%NAMES_I1D_STRING(NPOSN_1:NPOSN_2)               !<-- The variable's name
        LENGTH=wrt_int_state%LENGTH_DATA_I1D(N)                            !<-- The variable's length in words
!
        IF(LENGTH==1)THEN
          N2=N2+1
          N2ISCALAR=N2ISCALAR+1
          VARINAME(N2ISCALAR)=TRIM(NAME)
          VARIVAL(N2ISCALAR)=wrt_int_state%ALL_DATA_I1D(N2)
          IF(VARINAME(N2ISCALAR)=='IHRST') then
            IDATE(4)=VARIVAL(N2ISCALAR)
          ENDIF
        ELSE
          N2IARY=N2IARY+1
          ARYINAME(N2IARY)=TRIM(NAME)
          ARYILEN(N2IARY)=LENGTH
!            write(0,*)'in I1D array,aryiname=',aryiname(N2IARY),'len=',aryilen(N2IARY),  &
!              wrt_int_state%ALL_DATA_I1D(N2+1:N2+length)
!
          DO N1=1,LENGTH
            N2=N2+1
            ARYIVAL(N1,N2IARY)=wrt_int_state%ALL_DATA_I1D(N2)              !<-- Extract the individual data from the data string
          ENDDO

          IF(ARYINAME(N2IARY)=='IDAT') THEN
            IDATE(1)=ARYIVAL(3,N2IARY)
            IDATE(2)=ARYIVAL(2,N2IARY)
            IDATE(3)=ARYIVAL(1,N2IARY)
            IDATE(7)=100.
          ENDIF
!
!           write(0,*)'in I1D array,aryival=',aryival(:,N2IARY)
        ENDIF
!
      ENDDO
!
!-----------------------------
!***  Add fcst_date into ARYI
!-----------------------------
!
      N2IARY=N2IARY+1
      ARYINAME(N2IARY)='FCSTDATE'
      ARYILEN(N2IARY)=7
      ARYIVAL(1:7,N2IARY)=FCSTDATE(1:7)
!
!-----------------------------------------------------------------------
!***  Real scalar/1-D history variables
!-----------------------------------------------------------------------
!
!------------------------------------------------------------
!***  Find the total number of real scalars and real arrays.
!------------------------------------------------------------
!
      N2RSCALAR=0
      N2RARY=0
      MAXLENGTH=1
!
      DO N=1,wrt_int_state%KOUNT_R1D(1)                                    !<-- Loop through all scalar/1D real data
        LENGTH=wrt_int_state%LENGTH_DATA_R1D(N)                            !<-- The variable's length
        IF(LENGTH==1)THEN
           N2RSCALAR=N2RSCALAR+1
        ELSE
          N2RARY=N2RARY+1
          MAXLENGTH=MAX(LENGTH,MAXLENGTH)
        ENDIF
      ENDDO
!
      ALLOCATE(VARRNAME(N2RSCALAR),VARRVAL(N2RSCALAR))
      ALLOCATE(ARYRNAME(N2RARY),ARYRLEN(N2RARY),ARYRVAL(MAXLENGTH,N2RARY))
!
!------------------------------------------------------
!***  Set values for the real scalars and real arrays.
!------------------------------------------------------
      N2=0                                                                 !<-- Word counter for full string of real scalar/1D data
      N2RSCALAR=0
      N2RARY=0
!
      DO N=1,wrt_int_state%KOUNT_R1D(1)                                    !<-- Loop through all scalar/1D real data
!
        NPOSN_1=(N-1)*ESMF_MAXSTR+1
        NPOSN_2=N*ESMF_MAXSTR
        NAME=wrt_int_state%NAMES_R1D_STRING(NPOSN_1:NPOSN_2)               !<-- The variable's name
        LENGTH=wrt_int_state%LENGTH_DATA_R1D(N)                            !<-- The variable's length
!
        IF(LENGTH==1)THEN
          N2=N2+1
          N2RSCALAR=N2RSCALAR+1
          VARRNAME(N2RSCALAR)=TRIM(NAME)
          VARRVAL(N2RSCALAR)=wrt_int_state%ALL_DATA_R1D(N2)
!
          IF( TRIM(NAME)=='PT') THEN
            NPT=N2RSCALAR
          ELSEIF ( TRIM(NAME)=='PDTOP') THEN
            NPDTOP=N2RSCALAR
          ELSEIF ( TRIM(NAME)=='DYH') THEN
            NDYH=N2RSCALAR
          ELSEIF ( TRIM(NAME)=='TPH0D') THEN
            TPH0D=VARRVAL(N2RSCALAR)
          ELSEIF ( trim(NAME)=='TLM0D') THEN
            TLM0D=VARRVAL(N2RSCALAR)
          ENDIF
!
        ELSE
          N2RARY=N2RARY+1
          ARYRNAME(N2RARY)=TRIM(NAME)
          ARYRLEN(N2RARY)=LENGTH
!
          DO N1=1,LENGTH
            N2=N2+1
            ARYRVAL(N1,N2RARY)=wrt_int_state%ALL_DATA_R1D(N2)              !<-- Extract the individual data from the data string
          ENDDO
!
          IF( TRIM(NAME)=='SG1') THEN
            NSG1=N2RARY
          ELSEIF ( TRIM(NAME)=='SG2') THEN
            NSG2=N2RARY
          ELSEIF ( TRIM(NAME)=='SGML1' ) THEN
            NSGML1=N2RARY
          ELSEIF (TRIM(NAME)=='SGML2') THEN
            NSGML2=N2RARY
          ELSEIF (TRIM(NAME)=='DXH') THEN
            NDXH=N2RARY
          ENDIF

        ENDIF
!
      ENDDO
!
!-----------------------------------------------------------------------
!***  Logical history variables
!-----------------------------------------------------------------------
!
      N2LSCALAR=wrt_int_state%KOUNT_LOG(1)                                 !<-- Counter for full string of logical data
!
      ALLOCATE(VARLNAME(N2LSCALAR),VARLVAL(N2LSCALAR))
      N2LSCALAR=0
!
      DO N=1,wrt_int_state%KOUNT_LOG(1)                                    !<-- Loop through all logical data
!
        NPOSN_1=(N-1)*ESMF_MAXSTR+1
        NPOSN_2=N*ESMF_MAXSTR
        NAME=wrt_int_state%NAMES_LOG_STRING(NPOSN_1:NPOSN_2)               !<-- The variable's name
!
        N2LSCALAR=N2LSCALAR+1
        WORK_LOGICAL = wrt_int_state%ALL_DATA_LOG(N2LSCALAR)

        VARLNAME(N2LSCALAR)=NAME

        VARLVAL(N2LSCALAR)=wrt_int_state%ALL_DATA_LOG(N2LSCALAR)

        IF(TRIM(NAME)=='GLOBAL') GLOBAL=WORK_LOGICAL
!
      ENDDO
!
!-----------------------------------------------------------------------
!***  Now open NEMSIO file.
!-----------------------------------------------------------------------
!
      N=LEN_TRIM(wrt_int_state%HST_NAME_BASE)
      INT_SEC=INT(wrt_int_state%NFSECONDS)
      FRAC_SEC=NINT((wrt_int_state%NFSECONDS-INT_SEC)*100.)
      WRITE(FILENAME,100)wrt_int_state%HST_NAME_BASE(1:N)//'_nio_'      &
                        ,wrt_int_state%NFHOURS,'h_'                     &
                        ,wrt_int_state%NFMINUTES,'m_'                   &
                        ,INT_SEC,'.',FRAC_SEC,'s'
      IF(wrt_int_state%PRINT_OUTPUT .OR. wrt_int_state%PRINT_ALL)       &
      write(0,*)'FILENAME=',trim(FILENAME),'n=',n
  100 FORMAT(A,I4.4,A,I2.2,A,I2.2,A,I2.2,A)
!
!----------------------------------------------------
!***  Prepare variables needed by the nemsip header:
!----------------------------------------------------
!
!dimension
      IF(GLOBAL) THEN
!for global im/jm for data field
        NFRAME=1
      ELSE
!for regional
        NFRAME=0
      ENDIF
      IM=wrt_int_state%im(1)
      JM=wrt_int_state%jm(1)
      DIM1=wrt_int_state%im(1)-2*NFRAME
      DIM2=wrt_int_state%jm(1)-2*NFRAME
!
      LM=wrt_int_state%LM(1)
!
!for nmmb whole domain
      FIELDSIZE=IM*JM
      NREC=wrt_int_state%kount_I2D(1)+wrt_int_state%kount_R2D(1)+1        !add fact10 for GSI
!
!vcoord
      ALLOCATE(VCOORD(LM+1,3,2))
      VCOORD=0.
      IF(NSG1>0.and.NSG2>0.and.NPDTOP>0.and.NPT>0) then
        VCOORD(1:LM+1,1,1)=0.1*(ARYRVAL(1:LM+1,NSG1)*VARRVAL(NPDTOP)      &
          -ARYRVAL(1:LM+1,NSG2)*(VARRVAL(NPDTOP)+VARRVAL(NPT))            &
          +VARRVAL(NPT) )
        VCOORD(1:LM+1,2,1)=ARYRVAL(1:LM+1,NSG2)
        VCOORD(1:LM+1,3,1)=0
      ENDIF
      IF(NSGML1>0.and.NSGML2>0.and.NPDTOP>0.and.NPT>0) then
        VCOORD(1:LM,1,2)=0.1*(ARYRVAL(1:LM,NSGML1)*VARRVAL(NPDTOP)        &
          -ARYRVAL(1:LM,NSGML2)*(VARRVAL(NPDTOP)+VARRVAL(NPT))            &
          +VARRVAL(NPT) )
        VCOORD(1:LM,2,2)=ARYRVAL(1:LM,NSGML2)
        VCOORD(1:LM,3,2)=0
      ENDIF
!!!   write(0,*)'after vcoord,count_I2d=',wrt_int_state%kount_I2D(1),'nrec=',nrec
!
!-----------------------------------------------------------------------
!***  Cut the output I2D array.
!-----------------------------------------------------------------------
!
      ALLOCATE(RECNAME(NREC),RECLEVTYP(NREC),RECLEV(NREC))
      NREC=0
      INI1=0                                                               !<-- # of 1 layer vars 
      INI2=0                                                               !<-- # of total layers of vars with lm layer
      INI3=0                                                               !<-- # of total layers of vars with lm+1 layer
!
      DO NFIELD=1,wrt_int_state%KOUNT_I2D(1)
!
        NREC=NREC+1
        NPOSN_1=(NFIELD-1)*ESMF_MAXSTR+1
        NPOSN_2=NFIELD*ESMF_MAXSTR
        NAME=wrt_int_state%NAMES_I2D_STRING(NPOSN_1:NPOSN_2)               !<-- The name of this 2D integer history quantity
        INDX_2D=index(NAME,"_2D")
!
        IF (INDX_2D > 0) THEN
          MODEL_LEVEL=NAME(INDX_2D-3:INDX_2D-1)
          RECLEV(NREC)=(ICHAR(MODEL_LEVEL(1:1))-48)*100+(ICHAR(MODEL_LEVEL(2:2))-48)*10+ICHAR(MODEL_LEVEL(3:3))-48
          RECNAME(NREC)=NAME(1:INDX_2D-5)
          RECLEVTYP(NREC)='mid_layer'
          IF (RECLEV(NREC)==LM+1) THEN 
            RECLEVTYP(NREC-LM:NREC)='layer'
            INI3=INI3+LM+1
            INI2=INI2-(LM+1)
          ENDIF
          INI2=INI2+1
        ELSE
          RECNAME(NREC)=TRIM(NAME)
          RECLEV(NREC)=1
          RECLEVTYP(NREC)='sfc'
          INI1=INI1+1
        ENDIF
!
        IF (RECNAME(NREC)=='ISLTYP') RECNAME(NREC)='sltyp'
        IF (RECNAME(NREC)=='IVGTYP') RECNAME(NREC)='vgtyp'
        IF (RECNAME(NREC)=='NCFRCV') RECNAME(NREC)='cfrcv'
        IF (RECNAME(NREC)=='NCFRST') RECNAME(NREC)='cfrst'
        CALL LOWERCASE(RECNAME(NREC))
      ENDDO
!
!!!   write(0,*)'after I2D,nrec=',nrec
!
!-----------------------------------------------------------------------
!*** Cut the output R2D array.
!-----------------------------------------------------------------------
!
      NSOIL=0
      IND1=0                                                               !<-- # of 1-layer vars
      IND2=0                                                               !<-- # of total layers for vars with lm layers 
      IND3=0                                                               !<-- # of total layers for vars with lm+1 layers
      IND4=0                                                               !<-- # of total layers for vars with nsoil layers
!
      DO NFIELD=1,wrt_int_state%KOUNT_R2D(1)
!
        NREC=NREC+1
        NPOSN_1=(NFIELD-1)*ESMF_MAXSTR+1
        NPOSN_2=NFIELD*ESMF_MAXSTR
        NAME=wrt_int_state%NAMES_R2D_STRING(NPOSN_1:NPOSN_2)  !<-- The name of this 2D integer history quantity
        INDX_2D=INDEX(NAME,"_2D")
!
        IF (INDX_2D > 0) THEN
          MODEL_LEVEL=NAME(INDX_2D-3:INDX_2D-1)
          RECLEV(NREC)=(ICHAR(MODEL_LEVEL(1:1))-48)*100+(ICHAR(MODEL_LEVEL(2:2))-48)*10+ICHAR(MODEL_LEVEL(3:3))-48
          RECNAME(NREC)=NAME(1:INDX_2D-5)
          RECLEVTYP(NREC)='mid layer'
          IF (RECNAME(NREC)=='SMC') NSOIL=NSOIL+1
          IF (RECNAME(NREC)=='W_TOT') RECNAME(NREC)='vvel'
          IF (RECNAME(NREC)=='CW') RECNAME(NREC)='clwmr'
          IF (RECNAME(NREC)=='U') RECNAME(NREC)='ugrd'
          IF (RECNAME(NREC)=='V') RECNAME(NREC)='vgrd'
          IF (RECNAME(NREC)=='T') RECNAME(NREC)='tmp'
          IF (RECNAME(NREC)=='Q') RECNAME(NREC)='spfh'
          IF (RECNAME(NREC)=='O3') RECNAME(NREC)='o3mr'
          IF (RECLEV(NREC)==LM+1) THEN
            RECLEVTYP(NREC-LM:NREC)='layer'
            IND3=IND3+LM+1
            IND2=IND2-(LM+1)
          ENDIF
          IF (RECNAME(NREC)=='PINT') THEN
          RECNAME(NREC)='pres'
          ELSE IF (RECNAME(NREC)=='SMC'.OR.RECNAME(NREC)=='SH2O'.or.RECNAME(NREC)=='STC') THEN 
             RECLEVTYP(NREC)='soil layer' 
             IND4=IND4+1
             IND2=IND2-1
          ENDIF
          IND2=IND2+1
        ELSE
          RECLEV(NREC)=1
          RECNAME(NREC)=TRIM(NAME)
          RECLEVTYP(NREC)='sfc'
!
          IF (INDEX(RECNAME(NREC),"10")>0) RECLEVTYP(NREC)='10 m above gnd'
          IF (RECNAME(NREC)=='PD') THEN
            RECNAME(NREC)='dpres'
            RECLEVTYP(NREC)='hybrid sig lev'
          ENDIF
!
          IF (RECNAME(NREC)=='SST') RECNAME(NREC)='tsea'
          IF (RECNAME(NREC)=='FIS') RECNAME(NREC)='hgt'
          IF (RECNAME(NREC)=='USTAR') RECNAME(NREC)='uustar'
          IF (RECNAME(NREC)=='Z0') RECNAME(NREC)='zorl'
          IND1=IND1+1
        ENDIF
!
        CALL LOWERCASE(RECNAME(NREC))
      ENDDO
!
!for fact10
      NREC=NREC+1
      RECNAME(NREC)='fact10'
      RECLEVTYP(NREC)='10 m above gnd'
      RECLEV(NREC)=1

!glat1d and glon1d
      ALLOCATE(GLAT1D(FIELDSIZE),GLON1D(FIELDSIZE))
      DEGRAD=90./ASIN(1.)
      glon1d=0.
      glat1d=0.
      NMETA=12
!
!dx and dy
      ALLOCATE(DX(FIELDSIZE),DY(FIELDSIZE))
!
      if(NDXH>0) then
       DO J=1,JM
       DO I=1,IM
         DX(I+(J-1)*IM)=ARYRVAL(J,NDXH)
       ENDDO
       ENDDO
!       write(0,*)'after dx=',maxval(dx),minval(dx)
      else
       NMETA=7
      endif
!
      if(NDYH>0) then
       DO I=1,FIELDSIZE
        DY(I)=VARRVAL(NDYH)
       ENDDO
!       write(0,*)'after dy=',maxval(dy),minval(dy)
      endif
!
!-----------------------------------------------------------------------
!                      SET UP NEMSIO WRITE
!-----------------------------------------------------------------------
!
      CALL NEMSIO_INIT(IRET=IRET)
!
!-----------------------------------------------------------------------
!***  Open NEMSIO run history file.
!-----------------------------------------------------------------------
!
      CALL NEMSIO_OPEN(NEMSIOFILE,trim(FILENAME),'write',iret,           &
        modelname="NMMB", gdatatype="bin4", idate=IDATE,nfhour=NF_HOURS, &
        nfminute=NF_MINUTES,nfsecondn=nint(NF_SECONDS*100),              &
        nfsecondd=100,dimx=DIM1,dimy=DIM2,dimz=LM,nframe=NFRAME,         &
        nmeta=NMETA,                                                     &
        nsoil=NSOIL,ntrac=3,nrec=nrec, ncldt=1,rlon_min=minval(glon1d),  &
        rlon_max=maxval(glon1d), rlat_max=maxval(glat1d),                &
        rlat_min=minval(glat1d),vcoord=vcoord,lon=glon1d,lat=glat1d,     &
        dx=dx,dy=dy,extrameta=.true.,nmetavari=N2ISCALAR,                &
        nmetavarr=N2RSCALAR,nmetavarl=N2LSCALAR,nmetaaryi=N2IARY,        &
        nmetaaryr=N2RARY,variname=VARINAME,varival=VARIVAL,              &
        varrname=VARRNAME,varrval=VARRVAL,varlname=VARLNAME,             &
        varlval=VARLVAL,aryiname=ARYINAME,aryilen=ARYILEN,               &
        aryival=ARYIVAL,aryrname=ARYRNAME,aryrlen=ARYRLEN,               &
        aryrval=ARYRVAL,recname=RECNAME,reclevtyp=RECLEVTYP,reclev=RECLEV)
!
!-----------------------------------------------------------------------
!***  Get variables needed by the .ctl file.
!-----------------------------------------------------------------------
!
      IF(wrt_int_state%WRITE_NEMSIOCTL.AND.wrt_int_state%MYPE==LEAD_WRITE_TASK)THEN
        CALL NEMSIO_GETFILEHEAD(NEMSIOFILE,TLMETA=TLMETA,FILE_ENDIAN=FILE_ENDIAN)
        DXCTL=MAXVAL(DX)*180./(A*PI)
        DYCTL=MAXVAL(DY)*180./(A*PI)
        CNT=INI1           & ! # of integer 1-layer fields
           +(INI2/LM)      & ! # of integer lm-layer fields
           +(INI3/(LM+1))  & ! # of integer lm+1-layer fields
           +IND1           & ! # of real 1-layer fields
           +(IND2/LM)      & ! # of real lm-layer fields
           +(IND3/(LM+1))  & ! # of real lm+1-layer fields
           +(IND4/NSOIL)   & ! # of real nsoil-layer fields
           +1                ! fact10
!
!-----------------------------------------------------------------------
!***  Write out NEMSIO ctl file.
!-----------------------------------------------------------------------
!
        CALL WRITE_NEMSIOCTL(GLOBAL,IHOUR_FCST,IDAY_FCST,IMONTH_FCST,   &
          IYEAR_FCST,FILENAME,TLMETA,IM,JM,LM,NSOIL,TLM0D,TPH0D,DXCTL,  &
          DYCTL,NF_HOURS,NREC,RECNAME,RECLEVTYP,CNT,FILE_ENDIAN,        &
          ID_DOMAIN)
      ENDIF
!
!-----------------------------------------------------------------------
!***  Clean up
!-----------------------------------------------------------------------
!
      DEALLOCATE(VCOORD,DX,DY)
      DEALLOCATE(VARINAME,VARIVAL,ARYINAME,ARYILEN,ARYIVAL)
      DEALLOCATE(VARRNAME,VARRVAL,ARYRNAME,ARYRLEN,ARYRVAL)
      DEALLOCATE(VARLNAME,VARLVAL)
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE WRITE_NEMSIO_RUNHISTORY_OPEN
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE WRITE_RUNRESTART_OPEN(WRT_INT_STATE                    &
                                      ,IYEAR_FCST                       &
                                      ,IMONTH_FCST                      &
                                      ,IDAY_FCST                        &
                                      ,IHOUR_FCST                       &
                                      ,IMINUTE_FCST                     &
                                      ,SECOND_FCST                      &
                                      ,NTIMESTEP                        &
                                      ,NF_HOURS                         &
                                      ,NF_MINUTES                       &
                                      ,NF_SECONDS                       &
                                      ,RST_FIRST                        &
                                      ,LEAD_WRITE_TASK )
!
!-----------------------------------------------------------------------
!***  Write out a binary run restart file.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(WRITE_INTERNAL_STATE),INTENT(INOUT) :: WRT_INT_STATE            !<-- The Write component's internal state
!
      INTEGER,INTENT(IN) :: IYEAR_FCST                                  &
                           ,IMONTH_FCST                                 &
                           ,IDAY_FCST                                   &
                           ,IHOUR_FCST                                  &
                           ,IMINUTE_FCST                                &
                           ,LEAD_WRITE_TASK                             &
                           ,NF_HOURS                                    &
                           ,NF_MINUTES                                  &
                           ,NTIMESTEP 
!
      LOGICAL,INTENT(IN) :: RST_FIRST
!
      REAL,INTENT(IN) :: NF_SECONDS,SECOND_FCST
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER :: N,N1,N2,NPOSN_1,NPOSN_2,LENGTH
      INTEGER :: NFIELD,RC
      CHARACTER(ESMF_MAXSTR)       :: NAME
      INTEGER,DIMENSION(:),POINTER :: WORK_ARRAY_I1D
      REAL(4),DIMENSION(:),POINTER :: WORK_ARRAY_R1D
      LOGICAL                      :: WRITE_LOGICAL
      LOGICAL                      :: WORK_LOGICAL
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Open the restart file and write the current forecast time
!***  and elapsed time.
!-----------------------------------------------------------------------
!
      CALL OPEN_RST_FILE(WRT_INT_STATE)
!
      WRITE(wrt_int_state%IO_RST_UNIT,iostat=RC)IYEAR_FCST
      WRITE(wrt_int_state%IO_RST_UNIT,iostat=RC)IMONTH_FCST
      WRITE(wrt_int_state%IO_RST_UNIT,iostat=RC)IDAY_FCST
      WRITE(wrt_int_state%IO_RST_UNIT,iostat=RC)IHOUR_FCST
      WRITE(wrt_int_state%IO_RST_UNIT,iostat=RC)IMINUTE_FCST
      WRITE(wrt_int_state%IO_RST_UNIT,iostat=RC)SECOND_FCST
      WRITE(wrt_int_state%IO_RST_UNIT,iostat=RC)NTIMESTEP
!
      IF(RST_FIRST .AND. (wrt_int_state%PRINT_OUTPUT .OR. wrt_int_state%PRINT_ALL) )THEN
        WRITE(0,*)' Wrote IYEAR_FCST to restart file unit ',wrt_int_state%IO_RST_UNIT
        WRITE(0,*)' Wrote IMONTH_FCST to restart file unit ',wrt_int_state%IO_RST_UNIT
        WRITE(0,*)' Wrote IDAY_FCST to restart file unit ',wrt_int_state%IO_RST_UNIT
        WRITE(0,*)' Wrote IHOUR_FCST to restart file unit ',wrt_int_state%IO_RST_UNIT
        WRITE(0,*)' Wrote IMINUTE_FCST to restart file unit ',wrt_int_state%IO_RST_UNIT
        WRITE(0,*)' Wrote SECOND_FCST to restart file unit ',wrt_int_state%IO_RST_UNIT
        WRITE(0,*)' Wrote NTIMESTEP to restart file unit ',wrt_int_state%IO_RST_UNIT
      ENDIF
!
!-----------------------------------------------------------------------
!***  Integer scalar/1-D restart variables
!-----------------------------------------------------------------------
!
        N2=0                                                               !<-- Word counter for full string of integer scalar/1D data
!
        DO N=1,wrt_int_state%RST_KOUNT_I1D(1)                              !<-- Loop through all scalar/1D integer data
!
          NPOSN_1=(N-1)*ESMF_MAXSTR+1
          NPOSN_2=N*ESMF_MAXSTR
          NAME=wrt_int_state%RST_NAMES_I1D_STRING(NPOSN_1:NPOSN_2)         !<-- The variable's name
          LENGTH=wrt_int_state%RST_LENGTH_DATA_I1D(N)                      !<-- The variable's length in words
          ALLOCATE(WORK_ARRAY_I1D(LENGTH),stat=RC)
!
          DO N1=1,LENGTH
            N2=N2+1
            WORK_ARRAY_I1D(N1)=wrt_int_state%RST_ALL_DATA_I1D(N2)          !<-- Extract the individual data from the data string
          ENDDO
!
          WRITE(wrt_int_state%IO_RST_UNIT,iostat=RC)WORK_ARRAY_I1D         !<-- Write out the data
!
          IF(RST_FIRST .AND. (wrt_int_state%PRINT_OUTPUT .OR. wrt_int_state%PRINT_ALL) )THEN
            WRITE(0,*)'Wrote ',TRIM(NAME),' to restart file unit ',wrt_int_state%IO_RST_UNIT
          ENDIF
!
          DEALLOCATE(WORK_ARRAY_I1D)
!
        ENDDO
!
!-----------------------------------------------------------------------
!***  Real scalar/1-D restart variables
!-----------------------------------------------------------------------
!
        N2=0                                                               !<-- Word counter for full string of real scalar/1D data
!
        DO N=1,wrt_int_state%RST_KOUNT_R1D(1)                              !<-- Loop through all scalar/1D real data
!
          NPOSN_1=(N-1)*ESMF_MAXSTR+1
          NPOSN_2=N*ESMF_MAXSTR
          NAME=wrt_int_state%RST_NAMES_R1D_STRING(NPOSN_1:NPOSN_2)         !<-- The variable's name
          LENGTH=wrt_int_state%RST_LENGTH_DATA_R1D(N)                      !<-- The variable's length
          ALLOCATE(WORK_ARRAY_R1D(LENGTH),stat=RC)
!
          DO N1=1,LENGTH
            N2=N2+1
            WORK_ARRAY_R1D(N1)=wrt_int_state%RST_ALL_DATA_R1D(N2)          !<-- Extract the individual data from the data string
          ENDDO
!
          WRITE(wrt_int_state%IO_RST_UNIT,iostat=RC)WORK_ARRAY_R1D         !<-- Write out the data
!
          IF(RST_FIRST .AND. (wrt_int_state%PRINT_OUTPUT .OR. wrt_int_state%PRINT_ALL) )THEN
            WRITE(0,*)'Wrote ',TRIM(NAME),' to restart file unit ',wrt_int_state%IO_RST_UNIT
          ENDIF
!
          DEALLOCATE(WORK_ARRAY_R1D)
!
        ENDDO
!
!-----------------------------------------------------------------------
!***  Logical restart variables
!-----------------------------------------------------------------------
!
        N2=0                                                               !<-- Counter for full string of logical data
!
        DO N=1,wrt_int_state%RST_KOUNT_LOG(1)                              !<-- Loop through all logical data
!
          NPOSN_1=(N-1)*ESMF_MAXSTR+1
          NPOSN_2=N*ESMF_MAXSTR
          NAME=wrt_int_state%RST_NAMES_LOG_STRING(NPOSN_1:NPOSN_2)         !<-- The variable's name
!
          N2=N2+1
          WORK_LOGICAL = wrt_int_state%RST_ALL_DATA_LOG(N2)
          WRITE_LOGICAL=WORK_LOGICAL                                       !<-- Convert from ESMF_Logical to F90 logical
!
          WRITE(wrt_int_state%IO_RST_UNIT,iostat=RC)WRITE_LOGICAL          !<-- Write out the data
!
          IF(RST_FIRST .AND. (wrt_int_state%PRINT_OUTPUT .OR. wrt_int_state%PRINT_ALL) )THEN
            WRITE(0,*)'Wrote ',TRIM(NAME),' to restart file unit ',wrt_int_state%IO_RST_UNIT
          ENDIF
!
        ENDDO
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE WRITE_RUNRESTART_OPEN
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE WRITE_NEMSIO_RUNRESTART_OPEN(WRT_INT_STATE             &
                                        ,NEMSIOFILE                     &
                                        ,IYEAR_FCST                     &
                                        ,IMONTH_FCST                    &
                                        ,IDAY_FCST                      &
                                        ,IHOUR_FCST                     &
                                        ,IMINUTE_FCST                   &
                                        ,SECOND_FCST                    &
                                        ,NTIMESTEP                      &
                                        ,NF_HOURS                       &
                                        ,NF_MINUTES                     &
                                        ,NF_SECONDS                     &
                                        ,DIM1,DIM2,NFRAME,GLOBAL        &
                                        ,ID_DOMAIN                      &
                                        ,LEAD_WRITE_TASK)
!
!-----------------------------------------------------------------------
!***  Write out a NEMSIO binary run history file.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(WRITE_INTERNAL_STATE),INTENT(INOUT) :: WRT_INT_STATE             !<-- The Write component's internal state
!
      TYPE(NEMSIO_GFILE),INTENT(INOUT)         :: NEMSIOFILE                !<-- The nemsio file handler
!
      INTEGER,INTENT(IN)  :: IYEAR_FCST                                 &
                            ,IMONTH_FCST                                &
                            ,IDAY_FCST                                  &
                            ,IHOUR_FCST                                 &
                            ,IMINUTE_FCST                               &
                            ,NF_HOURS                                   &
                            ,NF_MINUTES                                 &
                            ,ID_DOMAIN                                  &
                            ,LEAD_WRITE_TASK                            &
			    ,NTIMESTEP

      INTEGER,INTENT(OUT) :: DIM1,DIM2,NFRAME         
      LOGICAL,INTENT(OUT) :: GLOBAL
!
      REAL,INTENT(IN)     :: NF_SECONDS                                 &
                            ,SECOND_FCST
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER :: I,J,N,N1,N2,NPOSN_1,NPOSN_2,LENGTH,MAXLENGTH
!
      INTEGER :: FIELDSIZE,IM,JM,LM,IDATE(7),FCSTDATE(7)                &
                ,INDX_2D,IRET,IND1,IND2,IND3,IND4,IND5,CNT              &
 		,INI1,INI2,INI3                                         &
                ,N2ISCALAR,N2IARY,N2RSCALAR,N2RARY,N2LSCALAR            &
                ,NMETA,NSOIL,TLMETA,VLEV    
!
      INTEGER :: FRAC_SEC,INT_SEC,NFIELD,RC
!
      INTEGER,DIMENSION(:),POINTER :: ARYILEN                           &
                                     ,ARYRLEN                           &
                                     ,RECLEV                            &
                                     ,VARIVAL
!
      INTEGER,DIMENSION(:,:),POINTER :: ARYIVAL
!
      REAL(4) :: DEGRAD,DXCTL,DYCTL,TPH0D,TLM0D
!
      REAL(4),DIMENSION(:),POINTER :: DX,DY,DXH                         &
                                     ,GLAT1D,GLON1D
!
      REAL(4),DIMENSION(:,:,:),POINTER :: VCOORD
!
      REAL(KIND=KFPT),DIMENSION(:)  ,POINTER :: VARRVAL
      REAL(KIND=KFPT),DIMENSION(:,:),POINTER :: ARYRVAL
!
      LOGICAL,DIMENSION(:),POINTER :: VARLVAL
!
      CHARACTER(6)  :: MODEL_LEVEL
      CHARACTER(16) :: VLEVTYP,FILE_ENDIAN
!
      CHARACTER(16),DIMENSION(:),POINTER :: ARYINAME                    &
                                           ,ARYRNAME                    &
                                           ,RECNAME                     &
                                           ,VARINAME                    &
                                           ,VARRNAME                    &
                                           ,VARLNAME
!
      CHARACTER(16),DIMENSION(:),POINTER :: RECLEVTYP
!
      CHARACTER(ESMF_MAXSTR) :: NAME,FILENAME

      LOGICAL            :: WORK_LOGICAL

      INTEGER :: NDYH=0,NDXH=0,NPT=0,NPDTOP=0,NREC=0                    &
                ,NSG1=0,NSG2=0,NSGML1=0,NSGML2=0
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      FCSTDATE(1)=IYEAR_FCST
      FCSTDATE(2)=IMONTH_FCST
      FCSTDATE(3)=IDAY_FCST
      FCSTDATE(4)=IHOUR_FCST
      FCSTDATE(5)=IMINUTE_FCST
      FCSTDATE(6)=NINT(SECOND_FCST*100.)
      FCSTDATE(7)=100
!
!-----------------------------------------------------------------------
!***  Integer scalar/1-D history variables
!-----------------------------------------------------------------------
!
!-------------------------------------------------------------
!*** Find out the total number of int scalars and int arrays.
!-------------------------------------------------------------
!
      N2ISCALAR=0
      N2IARY=0
      MAXLENGTH=1
!
      DO N=1,wrt_int_state%RST_KOUNT_I1D(1)                                    !<-- Loop through all scalar/1D integer data
        LENGTH=wrt_int_state%RST_LENGTH_DATA_I1D(N)
!
        IF(LENGTH==1)THEN
          N2ISCALAR=N2ISCALAR+1
        ELSE
          N2IARY=N2IARY+1
          MAXLENGTH=MAX(LENGTH,MAXLENGTH)
        ENDIF
!
      ENDDO
!
      N2ISCALAR=N2ISCALAR+1
      N2IARY=N2IARY+1
      MAXLENGTH=MAX(MAXLENGTH,7)
      ALLOCATE(VARINAME(N2ISCALAR),VARIVAL(N2ISCALAR))
      ALLOCATE(ARYINAME(N2IARY),ARYILEN(N2IARY),ARYIVAL(MAXLENGTH,N2IARY))
!
!---------------------------------------
!***  Set value to AVRIVAL and ARYIVAL.
!---------------------------------------
!
      N2=0                                                                 !<-- Word counter for full string of integer scalar/1D data
      N2ISCALAR=0
      N2IARY=0
      IDATE=0
!
      DO N=1,wrt_int_state%RST_KOUNT_I1D(1)                                    !<-- Loop through all scalar/1D integer data
!
        NPOSN_1=(N-1)*ESMF_MAXSTR+1
        NPOSN_2=N*ESMF_MAXSTR
        NAME=wrt_int_state%RST_NAMES_I1D_STRING(NPOSN_1:NPOSN_2)               !<-- The variable's name
        LENGTH=wrt_int_state%RST_LENGTH_DATA_I1D(N)                            !<-- The variable's length in words
!
        IF(LENGTH==1)THEN
          N2=N2+1
          N2ISCALAR=N2ISCALAR+1
          VARINAME(N2ISCALAR)=TRIM(NAME)
          VARIVAL(N2ISCALAR)=wrt_int_state%RST_ALL_DATA_I1D(N2)
!     if(n2==10)then
!       write(0,*)' WRITE_NEMSIO_RUNRESTART_OPEN n2iscalar=',n2iscalar &
!                ,' VARIVAL(N2ISCALAR)=',VARIVAL(N2ISCALAR),' variname=',variname(n2iscalar)
!     endif
          IF(VARINAME(N2ISCALAR)=='IHRST') then
            IDATE(4)=VARIVAL(N2ISCALAR)
          ENDIF
        ELSE
          N2IARY=N2IARY+1
          ARYINAME(N2IARY)=TRIM(NAME)
          ARYILEN(N2IARY)=LENGTH
!            write(0,*)'in I1D array,aryiname=',aryiname(N2IARY),'len=',aryilen(N2IARY),  &
!              wrt_int_state%RST_ALL_DATA_I1D(N2+1:N2+length)
!
          DO N1=1,LENGTH
            N2=N2+1
            ARYIVAL(N1,N2IARY)=wrt_int_state%RST_ALL_DATA_I1D(N2)              !<-- Extract the individual data from the data string
          ENDDO

          IF(ARYINAME(N2IARY)=='IDAT') THEN
            IDATE(1)=ARYIVAL(3,N2IARY)
            IDATE(2)=ARYIVAL(2,N2IARY)
            IDATE(3)=ARYIVAL(1,N2IARY)
            IDATE(7)=100.
          ENDIF
!
!           write(0,*)'in I1D array,aryival=',aryival(:,N2IARY)
        ENDIF
!
      ENDDO
!------------------------------
!***  Add ntimestep into VARI.
!------------------------------
!
     N2ISCALAR= N2ISCALAR+1
     VARINAME(N2ISCALAR)='NTIMESTEP'
     VARIVAL(N2ISCALAR)=NTIMESTEP
!     write(0,*)'in I1D scalar,varival=',varival,'varname=',variname
!
!
!------------------------------
!***  Add fcst_date into ARYI.
!------------------------------
!
      N2IARY=N2IARY+1
      ARYINAME(N2IARY)='FCSTDATE'
      ARYILEN(N2IARY)=7
      ARYIVAL(1:7,N2IARY)=FCSTDATE(1:7)
!
!-----------------------------------------------------------------------
!***  Real scalar/1-D history variables
!-----------------------------------------------------------------------
!
!------------------------------------------------------------
!***  Find the total number of real scalars and real arrays.
!------------------------------------------------------------
!
      N2RSCALAR=0
      N2RARY=0
      MAXLENGTH=1
!
      DO N=1,wrt_int_state%RST_KOUNT_R1D(1)                                    !<-- Loop through all scalar/1D real data
        LENGTH=wrt_int_state%RST_LENGTH_DATA_R1D(N)                            !<-- The variable's length
        IF(LENGTH==1)THEN
           N2RSCALAR=N2RSCALAR+1
        ELSE
          N2RARY=N2RARY+1
          MAXLENGTH=MAX(LENGTH,MAXLENGTH)
        ENDIF
      ENDDO
!
      ALLOCATE(VARRNAME(N2RSCALAR),VARRVAL(N2RSCALAR))
      ALLOCATE(ARYRNAME(N2RARY),ARYRLEN(N2RARY),ARYRVAL(MAXLENGTH,N2RARY))
!
!------------------------------------------------------
!***  Set values for the real scalars and real arrays.
!------------------------------------------------------
      N2=0                                                                 !<-- Word counter for full string of real scalar/1D data
      N2RSCALAR=0
      N2RARY=0
!
      DO N=1,wrt_int_state%RST_KOUNT_R1D(1)                                    !<-- Loop through all scalar/1D real data
!
        NPOSN_1=(N-1)*ESMF_MAXSTR+1
        NPOSN_2=N*ESMF_MAXSTR
        NAME=wrt_int_state%RST_NAMES_R1D_STRING(NPOSN_1:NPOSN_2)               !<-- The variable's name
        LENGTH=wrt_int_state%RST_LENGTH_DATA_R1D(N)                            !<-- The variable's length
!
        IF(LENGTH==1)THEN
          N2=N2+1
          N2RSCALAR=N2RSCALAR+1
          VARRNAME(N2RSCALAR)=TRIM(NAME)
          VARRVAL(N2RSCALAR)=wrt_int_state%RST_ALL_DATA_R1D(N2)
!
          IF( TRIM(NAME)=='PT') THEN
            NPT=N2RSCALAR
          ELSEIF ( TRIM(NAME)=='PDTOP') THEN
            NPDTOP=N2RSCALAR
          ELSEIF ( TRIM(NAME)=='DYH') THEN
            NDYH=N2RSCALAR
          ELSEIF ( TRIM(NAME)=='TPH0D') THEN
            TPH0D=VARRVAL(N2RSCALAR)
          ELSEIF ( trim(NAME)=='TLM0D') THEN
            TLM0D=VARRVAL(N2RSCALAR)
          ENDIF
!
        ELSE
          N2RARY=N2RARY+1
          ARYRNAME(N2RARY)=TRIM(NAME)
          ARYRLEN(N2RARY)=LENGTH
!
          DO N1=1,LENGTH
            N2=N2+1
            ARYRVAL(N1,N2RARY)=wrt_int_state%RST_ALL_DATA_R1D(N2)              !<-- Extract the individual data from the data string
          ENDDO
!
          IF( TRIM(NAME)=='SG1') THEN
            NSG1=N2RARY
          ELSEIF ( TRIM(NAME)=='SG2') THEN
            NSG2=N2RARY
          ELSEIF ( TRIM(NAME)=='SGML1' ) THEN
            NSGML1=N2RARY
          ELSEIF (TRIM(NAME)=='SGML2') THEN
            NSGML2=N2RARY
          ELSEIF (TRIM(NAME)=='DXH') THEN
            NDXH=N2RARY
          ENDIF

        ENDIF
!
      ENDDO
!
!-----------------------------------------------------------------------
!***  Logical history variables
!-----------------------------------------------------------------------
!
      N2LSCALAR=wrt_int_state%RST_KOUNT_LOG(1)                                 !<-- Counter for full string of logical data
!
      ALLOCATE(VARLNAME(N2LSCALAR),VARLVAL(N2LSCALAR))
      N2LSCALAR=0
!
      DO N=1,wrt_int_state%RST_KOUNT_LOG(1)                                    !<-- Loop through all logical data
!
        NPOSN_1=(N-1)*ESMF_MAXSTR+1
        NPOSN_2=N*ESMF_MAXSTR
        NAME=wrt_int_state%RST_NAMES_LOG_STRING(NPOSN_1:NPOSN_2)               !<-- The variable's name
!
        N2LSCALAR=N2LSCALAR+1

        WORK_LOGICAL = wrt_int_state%RST_ALL_DATA_LOG(N2LSCALAR)

        VARLNAME(N2LSCALAR)=NAME

        VARLVAL(N2LSCALAR)=wrt_int_state%RST_ALL_DATA_LOG(N2LSCALAR)

        IF(TRIM(NAME)=='GLOBAL') GLOBAL=WORK_LOGICAL
!
      ENDDO
!
!-----------------------------------------------------------------------
!***  Now open NEMSIO file.
!-----------------------------------------------------------------------
!
      N=LEN_TRIM(wrt_int_state%RST_NAME_BASE)
      INT_SEC=INT(wrt_int_state%NFSECONDS)
      FRAC_SEC=NINT((wrt_int_state%NFSECONDS-INT_SEC)*100.)
      WRITE(FILENAME,100)wrt_int_state%RST_NAME_BASE(1:N)//'_nio_'    &
                        ,wrt_int_state%NFHOURS,'h_'                   &
                        ,wrt_int_state%NFMINUTES,'m_'                 &
                        ,INT_SEC,'.',FRAC_SEC,'s'
      IF(wrt_int_state%PRINT_OUTPUT .OR. wrt_int_state%PRINT_ALL)     &
      write(0,*)'FILENAME=',trim(FILENAME),'n=',n
  100 FORMAT(A,I4.4,A,I2.2,A,I2.2,A,I2.2,A)
!
!----------------------------------------------------
!***  Prepare variables needed by the nemsip header:
!----------------------------------------------------
!
!dimension
      IF(GLOBAL) THEN
!for global im/jm for data field
        NFRAME=1
      ELSE
!for regional
        NFRAME=0
      ENDIF
      IM=wrt_int_state%im(1)
      JM=wrt_int_state%jm(1)
      DIM1=wrt_int_state%im(1)-2*NFRAME
      DIM2=wrt_int_state%jm(1)-2*NFRAME
!
      LM=wrt_int_state%LM(1)
!
!for nmmb trimmed domain
      FIELDSIZE=IM*JM
      NREC=wrt_int_state%RST_KOUNT_I2D(1)+wrt_int_state%RST_KOUNT_R2D(1)+2 !add fact10 for GSI
                                                                           !add hgt for unified code
!
!vcoord
      ALLOCATE(VCOORD(LM+1,3,2))
      VCOORD=0.
      if(NSG1>0.and.NSG2>0.and.NPT>0.and.NPDTOP>0 ) THEN
        VCOORD(1:LM+1,1,1)=0.1*(ARYRVAL(1:LM+1,NSG1)*VARRVAL(NPDTOP)      &
          -ARYRVAL(1:LM+1,NSG2)*(VARRVAL(NPDTOP)+VARRVAL(NPT))            &
          +VARRVAL(NPT) )
        VCOORD(1:LM+1,2,1)=ARYRVAL(1:LM+1,NSG2)
        VCOORD(1:LM+1,3,1)=0
      ENDIF
      if(NSGML1>0.and.NSGML2>0.and.NPT>0.and.NPDTOP>0 ) THEN
        VCOORD(1:LM,1,2)=0.1*(ARYRVAL(1:LM,NSGML1)*VARRVAL(NPDTOP)        &
          -ARYRVAL(1:LM,NSGML2)*(VARRVAL(NPDTOP)+VARRVAL(NPT))            &
          +VARRVAL(NPT) )
        VCOORD(1:LM,2,2)=ARYRVAL(1:LM,NSGML2)
        VCOORD(1:LM,3,2)=0
      ENDIF
!!!     write(0,*)'after vcoord,count_I2d=',wrt_int_state%rst_kount_I2D(1),'nrec=',nrec
!
!-----------------------------------------------------------------------
!***  Cut the output I2D array.
!-----------------------------------------------------------------------
!
      ALLOCATE(RECNAME(NREC),RECLEVTYP(NREC),RECLEV(NREC))
      NREC=0
      INI1=0                                                               !<-- # of 1-layer vars
      INI2=0                                                               !<-- # of total layes of vars with lm layers
      INI3=0                                                               !<-- # of total layes of vars with lm+1 layers
!
      DO NFIELD=1,wrt_int_state%RST_KOUNT_I2D(1)
!
        NREC=NREC+1
        NPOSN_1=(NFIELD-1)*ESMF_MAXSTR+1
        NPOSN_2=NFIELD*ESMF_MAXSTR
        NAME=wrt_int_state%RST_NAMES_I2D_STRING(NPOSN_1:NPOSN_2)           !<-- The name of this 2D integer history quantity
        INDX_2D=index(NAME,"_2D")
!
        IF (INDX_2D > 0) THEN
          MODEL_LEVEL=NAME(INDX_2D-3:INDX_2D-1)
          RECLEV(NREC)=(ICHAR(MODEL_LEVEL(1:1))-48)*100+(ICHAR(MODEL_LEVEL(2:2))-48)*10+ICHAR(MODEL_LEVEL(3:3))-48
          RECNAME(NREC)=NAME(1:INDX_2D-5)
          RECLEVTYP(NREC)='mid_layer'
          IF (RECLEV(NREC)==LM+1) THEN
            RECLEVTYP(NREC-LM:NREC)='layer'
            INI3=INI3+LM+1
            INI2=INI2-(LM+1)
          ENDIF
          INI2=INI2+1
        ELSE
          RECNAME(NREC)=TRIM(NAME)
          RECLEV(NREC)=1
          RECLEVTYP(NREC)='sfc'
          INI1=INI1+1
        ENDIF
!
        IF (RECNAME(NREC)=='ISLTYP') RECNAME(NREC)='sltyp'
        IF (RECNAME(NREC)=='IVGTYP') RECNAME(NREC)='vgtyp'
        IF (RECNAME(NREC)=='NCFRCV') RECNAME(NREC)='cfrcv'
        IF (RECNAME(NREC)=='NCFRST') RECNAME(NREC)='cfrst'
        CALL LOWERCASE(RECNAME(NREC))
      ENDDO
!
!     write(0,*)'after I2D,nrec=',nrec
!
!-----------------------------------------------------------------------
!*** Cut the output R2D array.
!-----------------------------------------------------------------------
!
      NSOIL=0
      IND1=0                                                               !<-- # of 1-layer vars
      IND2=0                                                               !<-- # of total layers for vars with lm layers 
      IND3=0                                                               !<-- # of total layers for vars with lm+1 layers
      IND4=0                                                               !<-- # of total layers for vars with nsoil layers
      IND5=0                                                               !<-- # of total layers for vars with lm-1 layers
!
      DO NFIELD=1,wrt_int_state%RST_KOUNT_R2D(1)
!
        NREC=NREC+1
        NPOSN_1=(NFIELD-1)*ESMF_MAXSTR+1
        NPOSN_2=NFIELD*ESMF_MAXSTR
        NAME=wrt_int_state%RST_NAMES_R2D_STRING(NPOSN_1:NPOSN_2)  !<-- The name of this 2D integer history quantity
        INDX_2D=INDEX(NAME,"_2D")
!
        IF (INDX_2D > 0) THEN
          MODEL_LEVEL=NAME(INDX_2D-3:INDX_2D-1)
          RECLEV(NREC)=(ICHAR(MODEL_LEVEL(1:1))-48)*100+(ICHAR(MODEL_LEVEL(2:2))-48)*10+ICHAR(MODEL_LEVEL(3:3))-48
          RECNAME(NREC)=NAME(1:INDX_2D-5)
          RECLEVTYP(NREC)='mid layer'
          IF (RECNAME(NREC)=='SMC') NSOIL=NSOIL+1
          IF (RECNAME(NREC)=='W') RECNAME(NREC)='vvel'
          IF (RECNAME(NREC)=='CW') RECNAME(NREC)='clwmr'
          IF (RECNAME(NREC)=='U') RECNAME(NREC)='ugrd'
          IF (RECNAME(NREC)=='V') RECNAME(NREC)='vgrd'
          IF (RECNAME(NREC)=='T') RECNAME(NREC)='tmp'
          IF (RECNAME(NREC)=='Q') RECNAME(NREC)='spfh'
          IF (RECNAME(NREC)=='O3') RECNAME(NREC)='o3mr'
          IF (RECLEV(NREC)==LM+1) THEN
            RECLEVTYP(NREC-LM:NREC)='layer'
            IND3=IND3+LM+1
            IND2=IND2-(LM+1)
          ENDIF
          IF (RECNAME(NREC)=='PSGDT') THEN
            RECLEVTYP(NREC)='layerm1'
            IF (RECLEV(NREC)==LM-1) THEN
              IND5=IND5+LM-1
              IND2=IND2-(LM-1)
            ENDIF
          ENDIF
          IF (RECNAME(NREC)=='PINT') THEN
            RECNAME(NREC)='pres'
          ELSE IF (RECNAME(NREC)=='SMC'.OR.RECNAME(NREC)=='SH2O'.or.RECNAME(NREC)=='STC') THEN
            RECLEVTYP(NREC)='soil layer'
            IND4=IND4+1
            IND2=IND2-1
          ENDIF
          IND2=IND2+1
        ELSE
          RECLEV(NREC)=1
          RECNAME(NREC)=TRIM(NAME)
          RECLEVTYP(NREC)='sfc'
!
          IF (INDEX(RECNAME(NREC),"10")>0) RECLEVTYP(NREC)='10 m above gnd'
          IF (RECNAME(NREC)=='PD') THEN
            RECNAME(NREC)='dpres'
            RECLEVTYP(NREC)='hybrid sig lev'
          ENDIF
!
          IF (RECNAME(NREC)=='SST') RECNAME(NREC)='tsea'
          IF (RECNAME(NREC)=='USTAR') RECNAME(NREC)='uustar'
          IF (RECNAME(NREC)=='Z0') RECNAME(NREC)='zorl'
          IND1=IND1+1
        ENDIF
!
        CALL LOWERCASE(RECNAME(NREC))
      ENDDO
!     write(0,*)'after R2D,nrec=',nrec,'kount_r2d=',wrt_int_state%RST_KOUNT_R2D(1)
!
!for fact10
      NREC=NREC+1
      RECNAME(NREC)='fact10'
      RECLEVTYP(NREC)='10 m above gnd'
      RECLEV(NREC)=1
!for hgt
      NREC=NREC+1
      RECNAME(NREC)='hgt'
      RECLEVTYP(NREC)='sfc'
      RECLEV(NREC)=1
!
!glat1d and glon1d
      ALLOCATE(GLAT1D(FIELDSIZE),GLON1D(FIELDSIZE))
      DEGRAD=90./ASIN(1.)
      glon1d=0.
      glat1d=0.
      NMETA=12
!     write(0,*)'after glat1d,NDYH=',ndyh
!
!dx and dy
      ALLOCATE(DX(FIELDSIZE),DY(FIELDSIZE))
!
      if(NDXH>0) then
       DO J=1,JM
       DO I=1,IM
         DX(I+(J-1)*IM)=ARYRVAL(J,NDXH)
       ENDDO
       ENDDO
!      write(0,*)'after dx=',maxval(dx),minval(dx),'dy=',maxval(dy),minval(dy)
      else
       NMETA=7
      endif
!
      if(NDYH>0) then
       DO I=1,FIELDSIZE
        DY(I)=VARRVAL(NDYH)
       ENDDO
      endif
!     write(0,*)'after DY,nrec=',nrec

!
!-----------------------------------------------------------------------
!                      SET UP NEMSIO WRITE
!-----------------------------------------------------------------------
!
      CALL NEMSIO_INIT(IRET=IRET)
!
!-----------------------------------------------------------------------
!***  Open NEMSIO file
!-----------------------------------------------------------------------
!
      CALL NEMSIO_OPEN(NEMSIOFILE,trim(FILENAME),'write',iret,           &
        modelname="NMMB", gdatatype="bin4", idate=IDATE,nfhour=NF_HOURS, &
        nfminute=NF_MINUTES,nfsecondn=nint(NF_SECONDS*100),              &
        nfsecondd=100,dimx=DIM1,dimy=DIM2,dimz=LM,nframe=NFRAME,         &
        nmeta=NMETA,                                                     &
        nsoil=NSOIL,ntrac=3,nrec=nrec, ncldt=1,rlon_min=minval(glon1d),  &
        rlon_max=maxval(glon1d), rlat_max=maxval(glat1d),                &
        rlat_min=minval(glat1d),vcoord=vcoord,lon=glon1d,lat=glat1d,     &
        dx=dx,dy=dy,extrameta=.true.,nmetavari=N2ISCALAR,                &
        nmetavarr=N2RSCALAR,nmetavarl=N2LSCALAR,nmetaaryi=N2IARY,        &
        nmetaaryr=N2RARY,variname=VARINAME,varival=VARIVAL,              &
        varrname=VARRNAME,varrval=VARRVAL,varlname=VARLNAME,             &
        varlval=VARLVAL,aryiname=ARYINAME,aryilen=ARYILEN,               &
        aryival=ARYIVAL,aryrname=ARYRNAME,aryrlen=ARYRLEN,               &
        aryrval=ARYRVAL,recname=RECNAME,reclevtyp=RECLEVTYP,reclev=RECLEV)
!
!     write(0,*)' WRITE_NEMSIO_RUNRESTART_OPEN after NEMSIO_OPEN variname(10)=',variname(10) &
!              ,' varival(10)=',varival(10)
!-----------------------------------------------------------------------
!***  Get variables needed by the .ctl file.
!-----------------------------------------------------------------------
!
      IF(wrt_int_state%WRITE_NEMSIOCTL.AND.wrt_int_state%MYPE==LEAD_WRITE_TASK)THEN
        CALL NEMSIO_GETFILEHEAD(NEMSIOFILE,TLMETA=TLMETA,FILE_ENDIAN=FILE_ENDIAN)
        DXCTL=MAXVAL(DX)*180./(A*PI)
        DYCTL=MAXVAL(DY)*180./(A*PI)
        CNT=INI1           & ! # of integer 1-layer fields
           +(INI2/LM)      & ! # of integer lm-layer fields
           +(INI3/(LM+1))  & ! # of integer lm+1-layer fields
           +IND1           & ! # of real 1-layer fields
           +(IND2/LM)      & ! # of real lm-layer fields
           +(IND3/(LM+1))  & ! # of real lm+1-layer fields
           +(IND4/NSOIL)   & ! # of real nsoil-layer fields
           +(IND5/(LM-1))  & ! # of real lm-1-layer fields
           +2                ! fact10 and hgt
!
!-----------------------------------------------------------------------
!***  Write out NEMSIO ctl file.
!-----------------------------------------------------------------------
!
        CALL WRITE_NEMSIOCTL(GLOBAL,IHOUR_FCST,IDAY_FCST,IMONTH_FCST,    &
          IYEAR_FCST,FILENAME,TLMETA,IM,JM,LM,NSOIL,TLM0D,TPH0D,DXCTL,   &
          DYCTL,NF_HOURS,NREC,RECNAME,RECLEVTYP,CNT,FILE_ENDIAN,         &
          ID_DOMAIN)
      ENDIF
!
!-----------------------------------------------------------------------
!***  Clean up
!-----------------------------------------------------------------------
!
      DEALLOCATE(VCOORD,DX,DY)
      DEALLOCATE(VARINAME,VARIVAL,ARYINAME,ARYILEN,ARYIVAL)
      DEALLOCATE(VARRNAME,VARRVAL,ARYRNAME,ARYRLEN,ARYRVAL)
      DEALLOCATE(VARLNAME,VARLVAL)
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE WRITE_NEMSIO_RUNRESTART_OPEN
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE WRITE_NEMSIOCTL(GLOBAL,IHOUR_FCST,IDAY_FCST,IMONTH_FCST, &
        IYEAR_FCST,FILENAME,TLMETA,DIM1,DIM2,LM,NSOIL,TLM0D,TPH0D,DXCTL,  &
        DYCTL,NF_HOURS,NREC,RECNAME,RECLEVTYP,KOUNT_R2D,FILE_ENDIAN,      &
        ID_DOMAIN)
!
!-----------------------------------------------------------------------
!***  Write out ctl file.
!-----------------------------------------------------------------------
!
!-----------------------
!*** Argument Variables
!-----------------------
!
      INTEGER,INTENT(IN) :: DIM1,DIM2                                   &
                           ,IHOUR_FCST,IDAY_FCST,IMONTH_FCST,IYEAR_FCST &
                           ,LM,NF_HOURS,NREC,NSOIL,TLMETA,KOUNT_R2D
!
      REAL,INTENT(IN) :: DXCTL,DYCTL,TLM0D,TPH0D
!
      LOGICAL,INTENT(IN) :: GLOBAL
!
      CHARACTER(*) ,INTENT(IN) :: FILENAME
      CHARACTER(16),INTENT(IN) :: RECNAME(:)
      CHARACTER(16),INTENT(IN) :: RECLEVTYP(:)
      CHARACTER(16),INTENT(IN) :: FILE_ENDIAN
!
!---------------------
!***  Local Variables
!---------------------
!
!
      INTEGER(KIND=KINT) :: IERR,RC 
      INTEGER ID_DOMAIN,N,NEWDIM1,NEWDIM2,IO_UNIT
!
      REAL CRSDXCTL,CRSDYCTL,RATIO,XBARLON,YBARLAT,MAXLAT,MAXLON,MINLAT,MINLON
!
      CHARACTER(3)  CMON
      CHARACTER(32) DATE
      CHARACTER(64) INFILE
!
      LOGICAL OPENED
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!---------------------
!***  Get unit number
!---------------------
!
      DO N=51,99
        INQUIRE(N,opened=OPENED)
          IF(.NOT.OPENED)THEN
            IO_UNIT=N
            EXIT
        ENDIF
      ENDDO
!
      CALL CMONTH(IMONTH_FCST,CMON)
      WRITE(DATE,'(I2.2,A,I2.2,A3,I4.4)')IHOUR_FCST,'Z',IDAY_FCST       &
           ,CMON,IYEAR_FCST
      OPEN(IO_UNIT,file=TRIM(FILENAME)//'.ctl',form='formatted')
!
      WRITE(IO_UNIT,105)TRIM(FILENAME)
      WRITE(IO_UNIT,106)
      WRITE(IO_UNIT,107)FILE_ENDIAN
      WRITE(IO_UNIT,108)TLMETA
      WRITE(IO_UNIT,109)
!
      IF (GLOBAL) THEN
        WRITE(IO_UNIT,121)DIM1,DXCTL
        WRITE(IO_UNIT,122)DIM2,DYCTl
      ELSE
        WRITE(IO_UNIT,110)DIM1,DIM2,TLM0D,TPH0D,DXCTL,DYCTL
!
!** Read min and max lat and lon values for each grid from a file
!
        WRITE(INFILE,'(A,I2.2)')'lat_lon_bnds_',ID_DOMAIN
        OPEN(UNIT=178,FILE=INFILE,STATUS='OLD',FORM='UNFORMATTED')
          READ (178)MINLAT,MAXLAT,MINLON,MAXLON 
        CLOSE(178)
!
!** Create proper nx(newdim1) and ny(newdim2) values for nests
!
        IF(ID_DOMAIN.GT.1)THEN
          XBARLON=(MAXLON*(180./pi)-MINLON*(180./pi))/2.
          YBARLAT=(MAXLAT*(180./pi)-MINLAT*(180./pi))/2.
          NEWDIM1=(MAXLON*(180./pi)-MINLON*(180./pi))/DXCTL
          NEWDIM2=(MAXLAT*(180./pi)-MINLAT*(180./pi))/DYCTL
          WRITE(IO_UNIT,131)NEWDIM1,TLM0D-XBARLON,DXCTL
          WRITE(IO_UNIT,132)NEWDIM2,TPH0D-YBARLAT,DYCTL
        ELSE
!
!** Create proper nx(newdim1) and ny(newdim2) values for parent domain
!
          NEWDIM1=(MAXLON*(180./pi)-MINLON*(180./pi))/DXCTL
          NEWDIM2=(MAXLAT*(180./pi)-MINLAT*(180./pi))/DYCTL
!
!** If parent domain nx and ny values are too big, reduce their values,
!** and force either nx or ny to be 1500 depending on which dimension is
!** the largest.  This also requires a coarser resolution.  Keep the aspect 
!** ratio of the domain the same (nx/ny).
!
!          IF(NEWDIM1.GT.1500.OR.NEWDIM2.GT.1500)THEN
!            IF(NEWDIM1.GT.NEWDIM2)THEN
!              RATIO=REAL(NEWDIM1)/REAL(NEWDIM2)
!              CRSDXCTL=(NEWDIM1*DXCTL)/1500.
!              CRSDYCTL=(NEWDIM2*DYCTL)/(1500./RATIO)
!              NEWDIM1=1500
!              NEWDIM2=REAL(NEWDIM1)/RATIO
!            ELSEIF(NEWDIM2.GE.NEWDIM1)THEN
!              RATIO=REAL(NEWDIM2)/REAL(NEWDIM1)
!              CRSDYCTL=(REAL(NEWDIM2)*DYCTL)/1500.
!              CRSDXCTL=(REAL(NEWDIM1)*DXCTL)/(1500./RATIO)
!              NEWDIM2=1500
!              NEWDIM1=REAL(NEWDIM2)/RATIO
!            ENDIF
!            WRITE(IO_UNIT,131)NEWDIM1,MINLON*(180./pi),CRSDXCTL
!            WRITE(IO_UNIT,132)NEWDIM2,MINLAT*(180./pi),CRSDYCTL
!          ELSE
            WRITE(IO_UNIT,131)NEWDIM1,MINLON*(180./pi),DXCTL
            WRITE(IO_UNIT,132)NEWDIM2,MINLAT*(180./pi),DYCTL
!          ENDIF
        ENDIF
      ENDIF ! global/regional
!
      WRITE(IO_UNIT,113)LM
      WRITE(IO_UNIT,114)1,TRIM(DATE)
!
 105  FORMAT('dset ^',A)
 106  FORMAT('undef -9.E+20')
 107  FORMAT('options ',A16,' sequential')
 108  FORMAT('fileheader',I12.0)
 109  FORMAT('title EXP1')

 110  FORMAT('pdef ',I6,I6,' eta.u ',f8.1,f8.1,f12.6,f12.6)
 121  FORMAT('xdef ',I6,' linear  -180.000 ',f12.6)
 122  FORMAT('ydef ',I6,' linear   -90.000 ',f12.6)
 131  FORMAT('xdef ',I6,' linear  ',f8.3,' ',f12.6)
 132  FORMAT('ydef ',I6,' linear  ',f8.3,' ',f12.6)
 113  FORMAT('zdef ',I6,' linear 1 1 ')
 114  FORMAT('tdef ',I6,' linear ',A12,' 6hr')
!
      WRITE(IO_UNIT,'(A,I6)')'VARS ',KOUNT_R2D

      N=1

      DO WHILE (N<=NREC)
        IF(RECLEVTYP(N)=='mid layer') THEN
          WRITE(IO_UNIT,'(A16,I3,A)')RECNAME(N),LM,' 99 mid layer'
          N=N+LM
        ELSEIF(RECLEVTYP(N)=='layerm1') THEN
          WRITE(IO_UNIT,'(A16,I3,A)')RECNAME(N),LM-1,' 99 layer'
          N=N+LM-1
        ELSEIF(RECLEVTYP(N)=='layer') THEN
          WRITE(IO_UNIT,'(A16,I3,A)')RECNAME(N),LM+1,' 99 layer'
          N=N+LM+1
        ELSEIF(RECLEVTYP(N)=='soil layer') THEN
          WRITE(IO_UNIT,'(A16,I3,A)')RECNAME(N),NSOIL,' 99 soil layer'
          N=N+NSOIL
        ELSE
          WRITE(IO_UNIT,'(A16,A)')RECNAME(N),'  0 99 sfc'
          N=N+1
        ENDIF
      ENDDO

      WRITE(IO_UNIT,'(A8)')'endvars'
      CLOSE(IO_UNIT)
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE WRITE_NEMSIOCTL
!
!-----------------------------------------------------------------------
!
      elemental subroutine lowercase(word)
!
!-----------------------------------------------------------------------
!***  convert a word to lower case
!-----------------------------------------------------------------------
!
      character (len=*) , intent(inout) :: word
      integer :: i,ic,nlen
      nlen = len(word)
!
      do i=1,nlen
        ic = ichar(word(i:i))
        if (ic >= 65 .and. ic < 91) word(i:i) = char(ic+32)
      end do
!
!
!-----------------------------------------------------------------------
!
      end subroutine lowercase
!
!-----------------------------------------------------------------------
!
      SUBROUTINE CMONTH(IMON,CMON)
!
!-----------------------------------------------------------------------
!***  Convert month
!-----------------------------------------------------------------------
!
      INTEGER,INTENT(IN) :: IMON
      CHARACTER(LEN=3)   :: CMON
!
!-----------------------------------------------------------------------
!
      SELECT CASE (IMON)
        CASE(1)
            CMON='Jan'
        CASE(2)
            CMON='Feb'
        CASE(3)
            CMON='Mar'
        CASE(4)
            CMON='Apr'
        CASE(5)
            CMON='May'
        CASE(6)
            CMON='Jun'
        CASE(7)
            CMON='Jul'
        CASE(8)
            CMON='Aug'
        CASE(9)
            CMON='Sep'
        CASE(10)
            CMON='Oct'
        CASE(11)
            CMON='Nov'
        CASE(12)
            CMON='Dec'
      END SELECT
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE CMONTH
!
!-----------------------------------------------------------------------
!
      END MODULE MODULE_WRITE_ROUTINES
!
!-----------------------------------------------------------------------
