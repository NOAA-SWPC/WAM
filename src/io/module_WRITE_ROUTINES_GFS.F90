!#include "../../ESMFVersionDefine.h"

!-----------------------------------------------------------------------
!
      MODULE MODULE_WRITE_ROUTINES_GFS
!
!-----------------------------------------------------------------------
!***  THIS MODULE CONTAINS ROUTINES NEEDED BY THE RUN STEP OF THE
!***  WRITE GRIDDED COMPONENT IN WHICH HISTORY OUTPUT DATA FROM
!***  THE FORECAST TASKS ARE ASSEMBLED AND WRITTEN TO HISTORY FILES
!***  BY THE WRITE TASKS.
!-----------------------------------------------------------------------
!***
!***  HISTORY   
!***
!       07 May 2009:  J. Wang - adopt write _routine from NMMB_io
!                                modified for GFS
!       03 Sep 2009:  W. Yang - Ensemble GEFS.
!       29 Sep 2010:  J. Wang - set up data mapping between fcst/write ps only once
!       16 Dec 2010:  J. Wang - change to nemsio library
!          Feb 2011:  W. Yang - Updated to use both the ESMF 4.0.0rp2 library,
!                               ESMF 5 library and the the ESMF 3.1.0rp2 library.
!       07 Mar 2011:  S. Lu - idvm is determined from sfcpress and thermodyn
!       27 Mar 2011:  J. Wang - set idsl for hybrid and sigms coord
!       13 May 2011:  W. yang - Modified for using the ESMF 5.2.0r_beta_snapshot_07.
!       05 May 2011:  J. Wang - add run post option on write quilt
!       28 Sep 2011:  W. yang - Modified for using the ESMF 5.2.0r library.
!       07 Nov 2012:  J. Wang  - generalize io for atmosphere
!--------------------------------------------------------------------------------
!
      USE ESMF
!
      USE MODULE_WRITE_INTERNAL_STATE_GFS
!
      USE MODULE_IO_MPI_DEF, ONLY :   MPI_COMM_INTER_ARRAY,MPI_COMM_COMP   &
                                     ,N_GROUP,ENSMEM_NAME
!
      USE MODULE_INCLUDE_IO 
!
      USE MODULE_ERR_MSG,ONLY: ERR_MSG,MESSAGE_CHECK
!
      USE NEMSIO_MODULE
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------

      PRIVATE
!
      PUBLIC :: FIRST_PASS_GFS, WRITE_NEMSIO_OPEN
!
!-----------------------------------------------------------------------
!
      INTEGER :: MYPE                                                     !<-- My MPI task ID
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE FIRST_PASS_GFS(IMP_STATE_WRITE, WRT_INT_STATE          &
                               ,NTASKS, MYPE, NBDL, NCURRENT_GROUP)
! 
!-----------------------------------------------------------------------
!***  EACH TIME A NEW GROUP OF WRITE TASKS IS INVOKED FOR THE FIRST
!***  TIME THIS ROUTINE WILL PERFORM CERTAIN COMPUTATIONS AND TASKS
!***  THAT ONLY NEED TO BE DONE ONCE FOR EACH WRITE GROUP.
!***  THE ROUTINE WILL UNLOAD SOME VALUES FROM THE WRITE COMPONENT'S
!***  IMPORT STATE.  BECAUSE THESE QUANTITIES DO NOT CHANGE WITH
!***  FORECAST TIME, THE ROUTINE IS NOT NEEDED FOR SUBSEQUENT USE
!***  OF EACH WRITE GROUP.  THE DATA BEING UNLOADED CONSISTS OF
!***  EVERYTHING EXCEPT THE 2D/3D GRIDDED FORECAST ARRAYS.
!***  ALSO BASIC INFORMATION IS PROVIDED TO FORECAST AND WRITE TASKS
!***  THAT WILL BE NECESSARY IN THE QUILTING OF LOCAL 2-D GRIDDED
!***  HISTORY DATA INTO FULL DOMAIN 2-D FIELDS.
!-----------------------------------------------------------------------
!
      INTEGER,INTENT(IN) :: NTASKS, MYPE, NCURRENT_GROUP, NBDL
!
      TYPE(ESMF_State)              ,INTENT(INOUT) :: IMP_STATE_WRITE
      TYPE(WRITE_INTERNAL_STATE_GFS),INTENT(INOUT) :: WRT_INT_STATE
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      TYPE(ESMF_FieldBundle)       :: FILE_BUNDLE
      INTEGER                      :: I,IERR,IM,J,JM,L,K,N,NN,NUM_ATTRIB   &
                                     ,NWTPG,RC,RC_WRT
!
      INTEGER,DIMENSION(:),POINTER :: ITMP
      REAL,   DIMENSION(:),POINTER :: RTMP
      TYPE(ESMF_Logical),DIMENSION(:),POINTER :: LTMP
!
      INTEGER                      :: JROW_FIRST,JROW_LAST,JROWS,          &
                                      LAST_FCST_TASK, LEAD_WRITE_TASK,     &
                                      LAST_WRITE_TASK
!
      INTEGER,SAVE                 :: NCHAR_I1D, NCHAR_R1D, NCHAR_LOG
!
      INTEGER                      :: JEND_WRITE,JSTA_WRITE, NUM_PES_FCST 
!
      INTEGER                      :: KOUNT_I1D, KOUNT_I2D, KOUNT_R1D      &
                                     ,KOUNT_R2D, KOUNT_LOG
!
      INTEGER                      :: LENGTH, LENGTH_SUM_I1D, LENGTH_SUM_R1D &
                                     ,LENGTH_SUM_LOG
!
      INTEGER                      :: NPOSN_START,NPOSN_END
!
      INTEGER                      :: NUM_FIELD_NAMES, MYPE_LOCAL   
!
      INTEGER,DIMENSION(MPI_STATUS_SIZE) :: JSTAT
!
      INTEGER,DIMENSION(:),POINTER :: LOCAL_ISTART, LOCAL_IEND          &
                                     ,LOCAL_JSTART, LOCAL_JEND
!jw:gfs
      INTEGER                      :: nbelt,nremain,istart,lat,nlat,    &
                                      TARGET_WRT, MPI_COMMUN
      INTEGER,DIMENSION(:),allocatable :: fcst_lat_for_write_task
      CHARACTER(NAME_MAXSTR),DIMENSION(:),allocatable :: field_name
      CHARACTER(NAME_MAXSTR*MAX_DATA_R2D)             :: NAMETMP
!
      INTEGER,DIMENSION(:),POINTER :: NCHAR_I2D, NCHAR_R2D
!
      INTEGER,DIMENSION(:),POINTER :: WORK_ARRAY_I1D
!
      REAL(4),DIMENSION(:),POINTER :: WORK_ARRAY_R1D
      REAL(8),DIMENSION(:),POINTER :: WORK_ARRAY_R1D8
!
      CHARACTER(ESMF_MAXSTR)       :: ATTRIB_NAME
!
      TYPE(ESMF_TypeKind_Flag)     :: DATATYPE
!
      TYPE(ESMF_Field)             :: FIELD_WORK1
 
      LOGICAL                         :: WORK_LOGICAL

      TYPE(ESMF_Logical),DIMENSION(1) :: NO_FIELDS
!
      TYPE(ESMF_VM)                   :: VM
!
!-----------------------------------------------------------------------
!
      REAL(KIND=kind_evod)              :: btim,btim0,timef,time1,time2
!
      btim0 = timef()
!-----------------------------------------------------------------------
!***** part 1: necessary info for fcst tasks to send data to write tasks
!***********************************************************************
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!***  FIRST WE NEED THE NUMBER OF WRITE TASKS IN EACH GROUP.
!-----------------------------------------------------------------------
!
      NUM_PES_FCST = wrt_int_state%NUM_PES_FCST
      NWTPG        = wrt_int_state%WRITE_TASKS_PER_GROUP
!
      IF(wrt_int_state%quilting) then
        LAST_FCST_TASK  = NTASKS-NWTPG-1
        LEAD_WRITE_TASK = LAST_FCST_TASK+1
        LAST_WRITE_TASK = NTASKS-1
        MYPE_LOCAL      = MOD(MYPE-NUM_PES_FCST,NWTPG)
      else                                                               !<-- for QUILTING=.false.,last pe will be the io pe
        LAST_FCST_TASK  = NTASKS-1
        LEAD_WRITE_TASK = LAST_FCST_TASK
        LAST_WRITE_TASK = LAST_FCST_TASK
        MYPE_LOCAL      = 0
      endif
!
!       write(0,*)'in first pass, LAST_FCST_TASK=',LAST_FCST_TASK,  &
!         'LEAD_WRITE_TASK=',LEAD_WRITE_TASK,'LAST_WRITE_TASK=',  &
!         LAST_WRITE_TASK,'NWTPG=',NWTPG,'NUM_PES_FCST=',NUM_PES_FCST, &
!        'NBDL=',NBDL,'NTASKS=',NTASKS

      ALLOCATE(ITMP(10000), RTMP(50000), LTMP(5000))
      ALLOCATE(NCHAR_I2D(1), NCHAR_R2D(1))
!
!-----------------------------------------------------------------------
!***  EXTRACT THE FULL DOMAIN LIMITS FROM THE COMPONENT'S IMPORT 
!***  STATE.  THESE ARE NEEDED FOR ALLOCATING THE WORKING ARRAYS 
!***  THAT WILL MOVE THE 2-D AND 3-D FIELDS FROM THE IMPORT TO THE
!***  EXPORT STATE.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
      domain_limits: IF(MYPE <= LAST_FCST_TASK) THEN                       !<-- This selects only forecast tasks to do extractions
                                                                           !    since only they know what is in the import state
!-----------------------------------------------------------------------
!*** get bundle
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK = "Extract History Bundle from Write Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_StateGet(state      =IMP_STATE_WRITE                  &  !<-- The write component's import state
                          ,itemName   =wrt_int_state%filename_base(NBDL)&  !<-- The name of the history data Bundle
                          ,fieldbundle=FILE_BUNDLE                      &  !<-- The history data Bundle inside the import state
                          ,rc         =RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!--------------------------------------------------------------------
        IF(NBDL == 1) THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK = "Extract Global Parameters from History Bundle" 
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!        CALL ESMF_AttributeGet(FIELDBUNDLE =FILE_BUNDLE                   &  !<-- The Bundle of history data
!                              ,name        ='lonf'                        &  !<-- Name of the Attribute to extract
!                              ,valueList   =wrt_int_state%IM              &  !<-- Extract this Attribute from History Bundle
!                              ,rc          =RC)
!
        CALL ESMF_AttributeGet(FIELDBUNDLE =FILE_BUNDLE                   &  !<-- The Bundle of history data
                              ,name        ='levs'                        &  !<-- Name of the Attribute to extract
                              ,valueList   =wrt_int_state%LM              &  !<-- Extract this Attribute from History Bundle
                              ,rc          =RC)

!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  NOW EXTRACT LOCAL SUBDOMAIN LIMITS.
!***  THESE WILL BE USED TO ALLOCATE THE WORKING ARRAY TO HOLD FIELDS
!***  ON EACH SUBDOMAIN PRIOR TO QUILTING THEM TOGETHER.
!***  WE FIRST NEED THE NUMBER OF FORECAST TASKS SINCE THAT
!***  DETERMINES THE SIZE OF THE ARRAYS HOLDING THE LOCAL
!***  SUBDOMAIN LIMITS.
!
!***  ALSO EXTRACT THE HALO DEPTHS SINCE THEY ARE NEEDED FOR
!***  EXCLUDING HALO POINTS FROM THE FINAL HISTORY DATA.
!
!***  THESE VALUES ARE NOT TO BE WRITTEN TO THE HISTORY FILES SO
!***  THEY WERE NOT INSERTED INTO THE HISTORY DATA Bundle INSIDE
!***  THE WRITE COMPONENT'S IMPORT STATE.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK = "Extract Local Quilting Info from Write Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
        if1tasks:  if(NTASKS >= 1 ) then
!
!-----------------------------------------------------------------------
!*** set up an array of the write_tasks for each fcst task to send data so
!*** the data on write PEs will be continuously distributed in the belted 
!*** subdomain.
!-----------------------------------------------------------------------
        if(.not.allocated(wrt_int_state%nlat_to_write_task)) then
          allocate(wrt_int_state%nlat_to_write_task(NWTPG))
          allocate(wrt_int_state%nstart_to_write_task(NWTPG))
          allocate(wrt_int_state%fcst_lat_to_write_task(wrt_int_state%lats_node_a))
          allocate(wrt_int_state%nwrttask_on_fcst(wrt_int_state%lats_node_a))
        endif

        allocate(fcst_lat_for_write_task(wrt_int_state%lats_node_a))
        NBELT   = wrt_int_state%jm(1)/NWTPG
        Nremain = mod(wrt_int_state%jm(1),NWTPG)
!
        do i=1,wrt_int_state%lats_node_a
          lat = wrt_int_state%global_lats_a(wrt_int_state%ipt_lats_node_a-1+i)
          if(lat <= nremain*(nbelt+1) ) then
            fcst_lat_for_write_task(i) = (lat-1)/(nbelt+1)+1
          else
            fcst_lat_for_write_task(i) = nremain+ (lat-nremain*(nbelt+1)-1)/Nbelt+1
          endif
        enddo
!       write(0,*)'nbelt=',nbelt,'nremain=',nremain,'lats_nodes=', &
!          wrt_int_state%lats_node_a,'ipt_lats_node_a=',wrt_int_state%ipt_lats_node_a, &
!          'global_lats_a=', &
!          wrt_int_state%global_lats_a(wrt_int_state%ipt_lats_node_a: &
!          wrt_int_state%ipt_lats_node_a-1+wrt_int_state%lats_node_a),'fcst_lat_to_write_task=',  &
!          wrt_int_state%fcst_lat_to_write_task
!
!*** on each forecast task, specify the starting point and the number of lats in global_lat_a to
!*** to each write task, which will decide the position of lat on write tasks
!*** 
        n = 0
        wrt_int_state%nstart_to_write_task(1) = 1
        wrt_int_state%fcst_lat_to_write_task(:) = -1.0E6
        do i=1,NWTPG
          wrt_int_state%nlat_to_write_task(i) = 0
          do j=1,wrt_int_state%lats_node_a
            if( fcst_lat_for_write_task(j) == i) then
              n = n+1
              wrt_int_state%nlat_to_write_task(i) = wrt_int_state%nlat_to_write_task(i)+1
              lat = j-1+wrt_int_state%ipt_lats_node_a
              wrt_int_state%fcst_lat_to_write_task(n) = wrt_int_state%global_lats_a(lat)
              wrt_int_state%nwrttask_on_fcst(n) = j
            endif
          enddo 
          if(i > 1) then
            wrt_int_state%nstart_to_write_task(i) = wrt_int_state%nstart_to_write_task(i-1)+ &
                        wrt_int_state%nlat_to_write_task(i-1)
          endif
       enddo
!        write(0,*)'nstart_write_task=',wrt_int_state%nstart_to_write_task(1:NWTPG),'nlat_write_task=', &
!           wrt_int_state%nlat_to_write_task(1:nwtpg),'fcst_lat_for_write_task=', &
!           wrt_int_state%fcst_lat_to_write_task(1:wrt_int_state%lats_node_a), &
!           'nwrttask_on_fcst=',wrt_int_state%nwrttask_on_fcst(1:wrt_int_state%lats_node_a)
!
        deallocate(fcst_lat_for_write_task)
!
!-----------------------------------------------------------------------
!
        endif if1tasks
!
      ENDIF !nbdl==1
!-----------------------------------------------------------------------
!
      ENDIF domain_limits
!
      time1 = timef()
!
!-----------------------------------------------------------------------
!----------------------------------------------------------------------
!***  FORECAST TASK 0 SENDS THE DOMAIN SIZE INFORMATION
!***  TO THE WRITE TASKS IN EACH WRITE GROUP BECAUSE
!***  THE WRITE TASKS NEED TO KNOW THIS TO ASSEMBLE THE
!***  FINAL GRIDDED DATA.
!***  FIRST WE NEED THE VM.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK = "Get the Current VM"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_VMGetCurrent(vm=VM                                      &  !<-- The ESMF virtual machine for this group of tasks
                            ,rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      if1task1:  if(wrt_int_state%quilting ) then
!
!-----------------------------------------------------------------------
!                            -- IM --
!-----------------------------------------------------------------------
!
      IF(MYPE == 0) THEN                                                   !<-- Forecast task 0 sends
        DO N=0,NWTPG-1
          CALL MPI_SEND(wrt_int_state%IM                                &  !<-- Send this data
                       ,1                                               &  !<-- Number of words sent
                       ,MPI_INTEGER                                     &  !<-- Datatype
                       ,N                                               &  !<-- Send to each of the write tasks (local IDs)
                       ,0                                               &  !<-- An MPI tag
                       ,MPI_COMM_INTER_ARRAY(NCURRENT_GROUP)            &  !<-- MPI communicator
                       ,IERR)
!
          IF(IERR /= 0) WRITE(0,*)' Failed to send IM from fcst task0 to write tasks'
!
        ENDDO
      ENDIF 
!
      IF(MYPE >= LEAD_WRITE_TASK) THEN                                        !<-- Write tasks in this group receive
!        write(0,*)'get dimension,NCURRENT_GROUP=',NCURRENT_GROUP,'MPI_COMM=',  &
!          MPI_COMM_INTER_ARRAY(NCURRENT_GROUP)
        CALL MPI_RECV(wrt_int_state%IM                                  &  !<-- Recv this data
                     ,1                                                 &  !<-- Words received
                     ,MPI_INTEGER                                       &  !<-- Datatype
                     ,0                                                 &  !<-- Recv from fcst 0
                     ,0                                                 &  !<-- An MPI tag
                     ,MPI_COMM_INTER_ARRAY(NCURRENT_GROUP)              &  !<-- MPI communicator
                     ,JSTAT                                             &  !<-- MPI status object
                     ,IERR)
!
        IF(IERR /= 0) WRITE(0,*)' Write tasks failed to receive IM from fcst task0'
!
      ENDIF
!
!-----------------------------------------------------------------------
!                          -- JM --
!-----------------------------------------------------------------------
!
      IF(MYPE == 0) THEN                                                      !<-- Forecast task 0 sends
        DO N=0,NWTPG-1                                            
          CALL MPI_SEND(wrt_int_state%JM                                &  !<-- Send this data
                       ,1                                               &  !<-- Number of words sent
                       ,MPI_INTEGER                                     &  !<-- Datatype
                       ,N                                               &  !<-- Send to each of the write tasks (local IDs)
                       ,0                                               &  !<-- An MPI tag
                       ,MPI_COMM_INTER_ARRAY(NCURRENT_GROUP)            &  !<-- MPI communicator
                       ,IERR)
!
          IF(IERR /= 0) WRITE(0,*)' Failed to send JM from fcst task0 to write tasks'
!
        ENDDO
      ENDIF 
!
      IF(MYPE >= LEAD_WRITE_TASK) THEN                                        !<-- Write tasks in this group receive
        CALL MPI_RECV(wrt_int_state%JM                                  &  !<-- Recv this data
                     ,1                                                 &  !<-- Words received
                     ,MPI_INTEGER                                       &  !<-- Datatype
                     ,0                                                 &  !<-- Recv from fcst 0
                     ,0                                                 &  !<-- An MPI tag
                     ,MPI_COMM_INTER_ARRAY(NCURRENT_GROUP)              &  !<-- MPI communicator
                     ,JSTAT                                             &  !<-- MPI status object
                     ,IERR)
!
        IF(IERR /= 0) WRITE(0,*)' Write tasks failed to receive JM from fcst task0'
!
      ENDIF
!
!-----------------------------------------------------------------------
!                          -- LM --
!-----------------------------------------------------------------------
!
      IF(MYPE == 0) THEN                                                   !<-- Forecast task 0 sends
        DO N=0,NWTPG-1                                            
          CALL MPI_SEND(wrt_int_state%LM                                &  !<-- Send this data
                       ,1                                               &  !<-- Number of words sent
                       ,MPI_INTEGER                                     &  !<-- Datatype
                       ,N                                               &  !<-- Send to each of the write tasks (local IDs)
                       ,0                                               &  !<-- An MPI tag
                       ,MPI_COMM_INTER_ARRAY(NCURRENT_GROUP)            &  !<-- MPI communicator
                       ,IERR)
!
          IF(IERR/=0)WRITE(0,*)' Failed to send LM from fcst task0 to write tasks'
!
        ENDDO
      ENDIF 
!
      IF(MYPE >= LEAD_WRITE_TASK)THEN                                        !<-- Write tasks in this group receive
        CALL MPI_RECV(wrt_int_state%LM                                  &  !<-- Recv this data
                     ,1                                                 &  !<-- Words received
                     ,MPI_INTEGER                                       &  !<-- Datatype
                     ,0                                                 &  !<-- Recv from fcst 0
                     ,0                                                 &  !<-- An MPI tag
                     ,MPI_COMM_INTER_ARRAY(NCURRENT_GROUP)              &  !<-- MPI communicator
                     ,JSTAT                                             &  !<-- MPI status object
                     ,IERR)
!
        IF(IERR /= 0) WRITE(0,*)' Write tasks failed to receive LM from fcst task0'
!
      ENDIF
!
      endif if1task1

      IM = wrt_int_state%IM(1)
      JM = wrt_int_state%JM(1)
!
!-----------------------------------------------------------------------
!*** set up the dimension for write pes
!-----------------------------------------------------------------------
!
      IF(MYPE >= LEAD_WRITE_TASK)THEN                                    !<-- Write tasks in this group receive

       if(.not. allocated(wrt_int_state%nlat_from_fcst_task) )then
         allocate(wrt_int_state%jstart_write(NWTPG))
         allocate(wrt_int_state%jend_write(NWTPG))
         allocate(wrt_int_state%nlat_from_fcst_task(NUM_PES_FCST))
         allocate(wrt_int_state%nstart_from_fcst_task(NUM_PES_FCST))
       endif
!
       nbelt   = jm/NWTPG
       nremain = mod(jm,NWTPG)

       DO I=1,NWTPG
         if(mod(i-1,NWTPG) < nremain) then
             wrt_int_state%JSTART_WRITE(i) = (i-1)*(nbelt+1) +1
             wrt_int_state%JEND_WRITE(i)   = wrt_int_state%JSTART_WRITE(i) + nbelt
          else
             wrt_int_state%JSTART_WRITE(i) = nremain*(nbelt+1)+(i-1-nremain)*nbelt+1
             wrt_int_state%JEND_WRITE(i)   = wrt_int_state%JSTART_WRITE(i)+nbelt-1
          endif
        ENDDO
!
       JSTA_WRITE = wrt_int_state%JSTART_WRITE(mype-lead_write_task+1)
       JEND_WRITE = wrt_int_state%JEND_WRITE(mype-lead_write_task+1)
       if(.not. allocated(wrt_int_state%fcst_lat_on_write_task) ) then
         allocate(wrt_int_state%fcst_lat_on_write_task(JSTA_WRITE:JEND_WRITE))
       endif
!
      ENDIF
!
!-----------------------------------------------------------------------
!***  get the subdomain information for write tasks
!-----------------------------------------------------------------------
!
      IF(wrt_int_state%quilting) THEN
       MPI_COMMUN = MPI_COMM_INTER_ARRAY(NCURRENT_GROUP)
      ELSE
       MPI_COMMUN = MPI_COMM_COMP
      ENDIF

      IF(MYPE < NUM_PES_FCST .and. MYPE /= LEAD_WRITE_TASK) THEN
        DO N=0,NWTPG-1
          TARGET_WRT = N
          IF(.not.wrt_int_state%quilting) TARGET_WRT = LEAD_WRITE_TASK
          CALL MPI_SEND(wrt_int_state%nlat_to_write_task(N+1)           &  !<-- Send this data
                       ,1                                               &  !<-- Number of words sent
                       ,MPI_INTEGER                                     &  !<-- Datatype
                       ,TARGET_WRT                                      &  !<-- Send to each of the write tasks (local IDs)
                       ,mype+1001                                       &  !<-- An MPI tag
                       ,MPI_COMMUN                                      &  !<-- MPI communicator
                       ,IERR)
!
          if (wrt_int_state%nlat_to_write_task(N+1) > 0) then
            ISTART = wrt_int_state%nstart_to_write_task(N+1)
            CALL MPI_SEND(wrt_int_state%fcst_lat_to_write_task(ISTART)  &  !<-- Send this data
                       ,wrt_int_state%nlat_to_write_task(N+1)           &  !<-- Number of words sent
                       ,MPI_INTEGER                                     &  !<-- Datatype
                       ,TARGET_WRT                                      &  !<-- Send to each of the write tasks (local IDs)
                       ,mype+1001                                       &  !<-- An MPI tag
                       ,MPI_COMMUN                                      &  !<-- MPI communicator
                       ,IERR)
!            write(0,*)'nlat on fcst=',wrt_int_state%nlat_write_task(N+1), &
!             'nstart=',wrt_int_state%nstart_write_task(N+1),              &
!             'fcst_lat=',wrt_int_state%fcst_lat_for_write_task(           &
!             wrt_int_state%nstart_write_task(N+1):wrt_int_state%nstart_write_task(N+1)+ &
!             wrt_int_state%nlat_write_task(N+1)-1)
!
            IF(IERR /= 0) WRITE(0,*)' Failed to send IM from fcst task0 to write tasks'
          endif
!
        ENDDO
     ENDIF
!
     IF(MYPE >= LEAD_WRITE_TASK) THEN
       DO I=1,NUM_PES_FCST
         IF(I /= LEAD_WRITE_TASK+1) THEN
           CALL MPI_RECV(wrt_int_state%nlat_from_fcst_task(I)           &  !<-- Recv this data
                     ,1                                                 &  !<-- Words received
                     ,MPI_INTEGER                                       &  !<-- Datatype
                     ,I-1                                               &  !<-- Recv from fcst 0
                     ,I+1000                                            &  !<-- An MPI tag
                     ,MPI_COMMUN                                        &  !<-- MPI communicator
                     ,JSTAT                                             &  !<-- MPI status object
                     ,IERR)
!
           if(I == 1) THEN
             wrt_int_state%nstart_from_fcst_task(I) = JSTA_WRITE
           else
             wrt_int_state%nstart_from_fcst_task(I) =                   &
                wrt_int_state%nstart_from_fcst_task(I-1) +              &
                wrt_int_state%nlat_from_fcst_task(I-1)
           endif
           ISTART = wrt_int_state%nstart_from_fcst_task(I)
           if(wrt_int_state%nlat_from_fcst_task(I) > 0) then
            CALL MPI_RECV(wrt_int_state%fcst_lat_on_write_task(ISTART)  &  !<-- Recv this data
                     ,wrt_int_state%nlat_from_fcst_task(I)              &  !<-- Words received
                     ,MPI_INTEGER                                       &  !<-- Datatype
                     ,I-1                                               &  !<-- Recv from fcst 0
                     ,I+1000                                            &  !<-- An MPI tag
                     ,MPI_COMMUN                                        &  !<-- MPI communicator
                     ,JSTAT                                             &  !<-- MPI status object
                     ,IERR) 
!
            IF(IERR /= 0)WRITE(0,*)' Failed to send IM from fcst task0 to write tasks'

          endif
         ELSE
           wrt_int_state%nlat_from_fcst_task(I) = wrt_int_state%nlat_to_write_task(1)
           if(ntasks > 1) then
              wrt_int_state%nstart_from_fcst_task(I) =                  &
                wrt_int_state%nstart_from_fcst_task(I-1) +              &
                wrt_int_state%nlat_from_fcst_task(I-1)
           else
              wrt_int_state%nstart_from_fcst_task(I) = 1
           endif
           ISTART = wrt_int_state%nstart_from_fcst_task(I)
           NLAT   = wrt_int_state%nlat_from_fcst_task(I)
           wrt_int_state%fcst_lat_on_write_task(ISTART:ISTART+NLAT-1) = &
             wrt_int_state%fcst_lat_to_write_task(1:NLAT)
         ENDIF
!
       ENDDO
!       
      ENDIF
!
      time2 = timef()
!
!-----------------------------------------------------------------------
!***************** part 2: bundles *************************************
!***********************************************************************
!*** each bundle is corresponding to one file, and 
!*** num_file is the total number of files that need to be opened
!-----------------------------------------------------------------------
!*** get bundle
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  THE NUMBER OF Attributes (FOR SCALARS AND 1D ARRAYS) AND
!***  Fields (FOR GRIDDED 2D ARRAYS) IN THE WRITE COMPONENT'S
!***  IMPORT STATE ARE NOT KNOWN A PRIORI.  IN ORDER TO TRANSFER
!***  THEM TO THE WRITE TASKS, EXTRACT THE NUMBER OF EACH OF
!***  THEM ALONG WITH THEIR NAMES.  THE SCALARS CAN BE LUMPED IN
!***  WITH THE 1D ARRAYS AT THIS POINT.
!
!***  EVEN THOUGH THESE COUNTS ARE JUST SCALAR INTEGERS THEIR
!***  POINTERS WERE ALLOCATED IN WRT_INIT TO LENGTH 1 SINCE THEY
!***  WILL BE USED IN ESMF_Send/Recv WHICH REQUIRE THEM TO BE 
!***  CONTIGUOUS DATA ARRAYS.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  ALL INTEGER QUANTITIES (AS 1D ARRAYS) AND 1D AND 2D REAL
!***  QUANTITIES WILL BE STRUNG TOGETHER IN SINGLE ARRAYS OF 
!***  EACH PARTICULAR TYPE.  ARRAYS THAT WILL HOLD THE LENGTH OF  
!***  EACH OF THE QUANTITIES IN THESE 'STRINGS' WERE ALLOCATED
!***  IN WRT_INIT.
!-----------------------------------------------------------------------
!
      KOUNT_I1D = 0
      KOUNT_R1D = 0
      KOUNT_LOG = 0
!
      LENGTH_SUM_I1D = 0
      LENGTH_SUM_R1D = 0
      LENGTH_SUM_LOG = 0
!
      NCHAR_I1D = 0
      NCHAR_R1D = 0
      NCHAR_LOG = 0
!
!-----------------------------------------------------------------------
!
      fcst_tasks: IF(MYPE <= LAST_FCST_TASK) THEN                          !<-- Only forecast tasks will extract output information
                                                                           !    from the import state because only they participated
                                                                           !    in filling the import state in the Dynamics/Physics
                                                                           !    components.
!
!-----------------------------------------------------------------------
!***  FIRST FIND THE NUMBER OF Attributes IN THE HISTORY DATA Bundle
!***  IN THE IMPORT STATE AND THEN FIND THEIR NAMES, LENGTHS, AND
!***  DATATYPES.
!***  EXTRACT THE INTEGER AND REAL DATA AND PACK IT INTO INTEGER
!***  AND REAL BUFFERS.  LATER THE BUFFERS WILL BE SENT FROM THE
!***  FORECAST TASKS (THE ONLY ONES WHO CAN SEE THE ORIGINAL
!***  DATA) TO THE WRITE TASKS.
!
!***  THE FACT THAT THE Attribute HISTORY DATA IS BEING COLLECTED
!***  HERE IN A BLOCK THAT EXECUTES ONLY ONCE PER WRITE GROUP
!***  IMPLIES THE ASSUMPTION THAT ONLY THE 2D/3D DATA
!***  ASSOCIATED WITH THE FORECAST GRID CAN CHANGE WITH TIME.
!***  IF ANY SCALAR/1D Attribute DATA CHANGE WITH TIME THEN
!***  THIS MUST BE MOVED OUT OF THIS 'FIRST' BLOCK.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK = "Get Attribute Count from History Bundle"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(FIELDBUNDLE =FILE_BUNDLE                      &  !<-- The write component's history data Bundle
                              ,count       =NUM_ATTRIB                       &  !<-- # of Attributes in the history data Bundle
                              ,rc          =RC)

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
          MESSAGE_CHECK = "Get Attribute Names, Datatypes, Lengths"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeGet(FIELDBUNDLE    =FILE_BUNDLE            &  !<-- The write component's history data Bundle
                                ,attributeIndex =N                      &  !<-- Index of each Attribute
                                ,name           =ATTRIB_NAME            &  !<-- Each Attribute's name
                                ,typekind       =DATATYPE               &  !<-- Each Attribute's ESMF Datatype
                                ,itemCount      =LENGTH                 &  !<-- Each Attribute's length
                                ,rc             =RC)

!    write(0,*)' length=',length,' ATTRIB_NAME=',&
!          ATTRIB_NAME,' n=',n,' mype=',mype
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!-----------------------------------------------------------------------
!                 -- SCALAR AND 1D INTEGER HISTORY DATA --
!-----------------------------------------------------------------------
!
          IF(DATATYPE == ESMF_TYPEKIND_I4) THEN                            !<-- Extract integer data with rank <2
!
            ALLOCATE(WORK_ARRAY_I1D(LENGTH),stat=RC)                       !<-- This length is from the preceding call 
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK = "Get Scalar/1-D Integer Data from History Bundle"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_AttributeGet(FIELDBUNDLE =FILE_BUNDLE               &  !<-- The write component's history data Bundle
                                  ,name        =ATTRIB_NAME               &  !<-- Name of the Attribute to extract
                                  ,itemCount   =LENGTH                    &  !<-- Length of Attribute
                                  ,valueList   =WORK_ARRAY_I1D            &  !<-- Place the Attribute here
                                  ,rc          =RC)

!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            KOUNT_I1D = KOUNT_I1D + 1                                      !<-- Count # of integer Attributes
!
            NPOSN_END   = KOUNT_I1D*NAME_MAXSTR
            NPOSN_START = NPOSN_END-NAME_MAXSTR + 1     
            wrt_int_state%NAMES_I1D_STRING(NBDL)(NPOSN_START:NPOSN_END) =       &
              ATTRIB_NAME(1:NAME_MAXSTR)                                   !<-- Save the 1D integer names
            NCHAR_I1D = NCHAR_I1D+NAME_MAXSTR                              !<-- Save #of characters in all scalar/1D integer names
                                                                           !    Note that each name is being given
                                                                           !    EMSF_MAXSTR total spaces
!
            DO L=1,LENGTH
              wrt_int_state%ALL_DATA_I1D(LENGTH_SUM_I1D+L,NBDL) = WORK_ARRAY_I1D(L)!<-- String together the integer data
            ENDDO
!
            LENGTH_SUM_I1D = LENGTH_SUM_I1D+LENGTH                         !<-- Total word sum of scalar/1D integer data
            wrt_int_state%LENGTH_DATA_I1D(KOUNT_I1D,NBDL) = LENGTH         !<-- Store length of each individual integer variable
!
            DEALLOCATE(WORK_ARRAY_I1D)
!
!-----------------------------------------------------------------------
!                  -- SCALAR AND 1D REAL HISTORY DATA --
!-----------------------------------------------------------------------
!
          ELSEIF(DATATYPE == ESMF_TYPEKIND_R4) THEN                           ! <-- Extract real data with rank <2
!
            ALLOCATE(WORK_ARRAY_R1D(LENGTH),stat=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK = "Get Scalar/1-D Real Data from History Bundle"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_AttributeGet(FIELDBUNDLE =FILE_BUNDLE               &  !<-- The write component's history data Bundle
                                  ,name        =ATTRIB_NAME               &  !<-- Name of the Attribute to extract
                                  ,itemCount   =LENGTH                    &  !<-- Length of AttributeME
                                  ,valueList   =WORK_ARRAY_R1D            &  !<-- Place the Attribute here
                                  ,rc          =RC)

!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            KOUNT_R1D = KOUNT_R1D + 1                                          !<-- Count # of real Attributes
!
            NPOSN_END   = KOUNT_R1D*NAME_MAXSTR
            NPOSN_START = NPOSN_END-NAME_MAXSTR+1     
            wrt_int_state%NAMES_R1D_STRING(NBDL)(NPOSN_START:NPOSN_END) =  &
               ATTRIB_NAME(1:NAME_MAXSTR)                                  !<-- sclar/1D real names
            NCHAR_R1D = NCHAR_R1D + NAME_MAXSTR                            !<-- Save #of characters in all scalar/1D real names
                                                                           !    Note that each name is being given
                                                                           !    NAME_MAXSTR total spaces
!
            DO L=1,LENGTH
              wrt_int_state%ALL_DATA_R1D(LENGTH_SUM_R1D+L,NBDL) = WORK_ARRAY_R1D(L)  !<-- String together the real data
            ENDDO
!
            LENGTH_SUM_R1D = LENGTH_SUM_R1D + LENGTH                       !<-- Total word sum of scalar/1D real data
            wrt_int_state%LENGTH_DATA_R1D(KOUNT_R1D,NBDL) = LENGTH         !<-- Store length of each individual real variable
!
            DEALLOCATE(WORK_ARRAY_R1D)
!
!
!-----------------------------------------------------------------------
!                  -- SCALAR AND 1D REAL8 HISTORY DATA --
!-----------------------------------------------------------------------
!
          ELSEIF(DATATYPE == ESMF_TYPEKIND_R8)THEN                         ! <-- Extract real data with rank <2
!
            ALLOCATE(WORK_ARRAY_R1D8(LENGTH),stat=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK = "Get Scalar/1-D Real8 Data from History Bundle"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_AttributeGet(FIELDBUNDLE =FILE_BUNDLE               &  !<-- The write component's history data Bundle
                                  ,name        =ATTRIB_NAME               &  !<-- Name of the Attribute to extract
                                  ,itemCount   =LENGTH                    &  !<-- Length of AttributeME
                                  ,valueList   =WORK_ARRAY_R1D8           &  !<-- Place the Attribute here
                                  ,rc          =RC)

!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            KOUNT_R1D   = KOUNT_R1D+1                                      !<-- Count # of real Attributes
!
            NPOSN_END   = KOUNT_R1D*NAME_MAXSTR
            NPOSN_START = NPOSN_END - NAME_MAXSTR + 1
            wrt_int_state%NAMES_R1D_STRING(NBDL)(NPOSN_START:NPOSN_END) =  &
               ATTRIB_NAME(1:NAME_MAXSTR)                                  !<-- sclar/1D real names
            NCHAR_R1D=NCHAR_R1D+NAME_MAXSTR                                !<-- Save #of characters in all scalar/1D real names
                                                                           !    Note that each name is being given
                                                                           !    NAME_MAXSTR total spaces
!change r8-->r4
            DO L=1,LENGTH

!     write(0,*)' mype=',mype,' l=',l,' length=',length,' work=',&
!               WORK_ARRAY_R1D8(L)

              wrt_int_state%ALL_DATA_R1D(LENGTH_SUM_R1D+L,NBDL) = WORK_ARRAY_R1D8(L)  !<-- String together the real data
            ENDDO
!
            LENGTH_SUM_R1D = LENGTH_SUM_R1D + LENGTH                       !<-- Total word sum of scalar/1D real data
            wrt_int_state%LENGTH_DATA_R1D(KOUNT_R1D,NBDL) = LENGTH         !<-- Store length of each individual real variable
!
            DEALLOCATE(WORK_ARRAY_R1D8)
!
!-----------------------------------------------------------------------
!                          -- LOGICAL DATA --                       
!-----------------------------------------------------------------------
!
!jw datatype== 9 for logical
!          ELSEIF(DATATYPE == ESMF_TYPEKIND_I1) THEN                       ! <-- Extract logical data
           ELSE
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Get Logical Data from History Bundle"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_AttributeGet(FIELDBUNDLE =FILE_BUNDLE               &!<-- The write component's history data Bundle
                                  ,name        =ATTRIB_NAME               &!<-- Name of the Attribute to extract
                                  ,value       =WORK_LOGICAL              &!<-- Place the Attribute here
                                  ,rc          =RC)

!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            KOUNT_LOG   = KOUNT_LOG + 1                                          !<-- Count # of logical Attributes
!
            NPOSN_END   = KOUNT_LOG*NAME_MAXSTR
            NPOSN_START = NPOSN_END - NAME_MAXSTR + 1     
            wrt_int_state%NAMES_LOG_STRING(NBDL)(NPOSN_START:NPOSN_END) = &
              ATTRIB_NAME(1:NAME_MAXSTR)                                   !<-- Save the logical names
            NCHAR_LOG = NCHAR_LOG + NAME_MAXSTR                            !<-- Save #of characters in all logical names
                                                                           !    Note that each name is being given
                                                                           !    EMSF_MAXSTR total spaces
!
            IF(WORK_LOGICAL) THEN
                wrt_int_state%ALL_DATA_LOG(KOUNT_LOG,NBDL) = ESMF_TRUE
            ELSE
                wrt_int_state%ALL_DATA_LOG(KOUNT_LOG,NBDL) = ESMF_FALSE
            END IF

!
            LENGTH_SUM_LOG = LENGTH_SUM_LOG + 1                            !<-- Total length of all logical data variables
!
          ENDIF
!
!-----------------------------------------------------------------------
!
        ENDDO attribute_loop
!
!-----------------------------------------------------------------------
!***  INSERT NUMBER AND LENGTHS OF SCALAR/1D INTEGER AND REAL QUANTITIES
!***  AND LOGICALS INTO THE WRITE COMPONENT'S INTERNAL STATE.
!-----------------------------------------------------------------------
!
        wrt_int_state%KOUNT_I1D(NBDL) = KOUNT_I1D
        wrt_int_state%KOUNT_R1D(NBDL) = KOUNT_R1D
        wrt_int_state%KOUNT_LOG(NBDL) = KOUNT_LOG
!
        wrt_int_state%LENGTH_SUM_I1D(NBDL) = LENGTH_SUM_I1D 
        wrt_int_state%LENGTH_SUM_R1D(NBDL) = LENGTH_SUM_R1D
        wrt_int_state%LENGTH_SUM_LOG(NBDL) = LENGTH_SUM_LOG
!
!-----------------------------------------------------------------------
!***  NOW EXTRACT THE NUMBER OF ESMF Fields IN THE HISTORY DATA Bundle
!***  WRITE COMPONENT'S IMPORT STATE ALONG WITH THEIR NAMES.
!***  SAVE THE Field INFORMATION OF THE 2D HISTORY DATA SINCE 
!***  IT WILL BE NEEDED FOR DATA EXTRACTION FROM THE IMPORT STATE
!***  AND THE TRANSFER TO THE WRITE TASKS.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  FIND OUT THE NUMBER OF ESMF Fields IN THE HISTORY DATA Bundle.
!***  IT WAS Fields INTO WHICH THE 2D GRIDDED HISTORY DATA WAS PLACED.
!
!***  THIS INFORMATION WILL BE SAVED IN THE INTERNAL STATE
!***  FOR ALL OUTPUT TIMES AND NOT BE RETRIEVED OVER AND OVER AGAIN.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK = "Get Field Count from History Bundle"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_FieldBundleGet(FILE_BUNDLE                   &  !<-- The write component's history data Bundle
                                ,fieldCount=wrt_int_state%NCOUNT_FIELDS(NBDL)   &  !<-- Get total # of Fields in the history data Bundle
                                ,rc        =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  NOW EXTRACT THE NAMES OF ALL THE Fields IN THE BUNDLE.  
!***  ALSO, THE NUMBER OF Field NAMES RETURNED SHOULD EQUAL 
!***  THE FIELD COUNT IN THE PRECEDING CALL.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK = "Extract Field Names from History Bundle"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        allocate(field_name(5000))
        CALL ESMF_FieldBundleGet(FILE_BUNDLE                         &  !<-- The write component's history data Bundle
                                ,FIELDNAMELIST =FIELD_NAME           &  !<-- Array of ESMF Field names in the Bundle
                                ,FIELDCOUNT    =NUM_FIELD_NAMES      &  !<-- Number of Field names in the Bundle
                                ,rc            =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        IF(NUM_FIELD_NAMES /= wrt_int_state%NCOUNT_FIELDS(nbdl)) THEN
          WRITE(0,*)' WARNING: Number of Fields in Bundle of history'   &
                   ,' output does not equal the number of Field names'
          WRITE(0,*)' They are ',NUM_FIELD_NAMES,' and '                &
                   ,wrt_int_state%NCOUNT_FIELDS(1),', respectively'
        ENDIF
!
!-----------------------------------------------------------------------
!***  DO A PRELIMINARY EXTRACTION OF THE Fields THEMSELVES IN ORDER TO
!***  COUNT THE NUMBER OF REAL AND INTEGER 2D ARRAYS.  
!-----------------------------------------------------------------------
!
        KOUNT_R2D    = 0
        KOUNT_I2D    = 0
        NCHAR_I2D(1) = 0
        NCHAR_R2D(1) = 0
!
        DO N=1,wrt_int_state%NCOUNT_FIELDS(NBDL)

!       write(0,*)' in the Nloop N=',N,' NCOUNT_FIELDS=',wrt_int_state%NCOUNT_FIELDS(NBDL)&
!    ,' fieldname=',FIELD_NAME(N),' mype=',mype
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Extract Fields from History Bundle for Counting"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_FieldBundleGet(FILE_BUNDLE               &  !<-- The write component's history data Bundle
                                  ,FIELDNAME =FIELD_NAME(N)  &  !<-- The ESMF Field's name
                                  ,field     =FIELD_WORK1    &  !<-- The ESMF Field taken from the Bundle
                                  ,rc        =RC)
         wrt_int_state%FIELD_NAME(N,NBDL) = field_name(N)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK = "Get Datatype of Fields for Counting Real/Integer"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_FieldGet(field=FIELD_WORK1                          &  !<-- The ESMF 2D Field
                            ,typekind=DATATYPE                          &  !<-- The Field's ESMF Datatype
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          IF(DATATYPE == ESMF_TYPEKIND_I4)THEN
            KOUNT_I2D = KOUNT_I2D + 1                                      !<-- Add up the total number of integer 2D Fields
!
            IF(KOUNT_I2D > MAX_DATA_I2D)THEN
              WRITE(0,*)' FATAL: YOU HAVE EXCEEDED MAX NUMBER OF INTEGER 2D FIELDS FOR OUTPUT'
              WRITE(0,*)' YOU MUST INCREASE VALUE OF MAX_DATA_I2D WHICH NOW EQUALS ',MAX_DATA_I2D
            ENDIF
!
            NPOSN_END   = KOUNT_I2D*NAME_MAXSTR
            NPOSN_START = NPOSN_END-NAME_MAXSTR + 1
            wrt_int_state%NAMES_I2D_STRING(NBDL)(NPOSN_START:NPOSN_END) =  &
              wrt_int_state%FIELD_NAME(N,NBDL)(1:NAME_MAXSTR) !<-- Save the 2D integer Field names 
                                                                                              !<-- in one long string
            NCHAR_I2D(1) = NCHAR_I2D(1) + NAME_MAXSTR
!
          ELSEIF(DATATYPE == ESMF_TYPEKIND_R4) THEN
            KOUNT_R2D = KOUNT_R2D + 1                                      !<-- Add up the total number of real 2D Fields
!
            IF(KOUNT_R2D > MAX_DATA_R2D) THEN
              WRITE(0,*)' FATAL: YOU HAVE EXCEEDED MAX NUMBER OF REAL 2D FIELDS FOR OUTPUT'
              WRITE(0,*)' KOUNT_R2D=',KOUNT_R2D
              WRITE(0,*)' YOU MUST INCREASE VALUE OF MAX_DATA_R2D WHICH NOW EQUALS ',MAX_DATA_R2D
            ENDIF
!
            NPOSN_END   = KOUNT_R2D*NAME_MAXSTR
            NPOSN_START = NPOSN_END - NAME_MAXSTR + 1
            wrt_int_state%NAMES_R2D_STRING(NBDL)(NPOSN_START:NPOSN_END) =  &
              wrt_int_state%FIELD_NAME(N,NBDL)(1:NAME_MAXSTR) !<-- Save the 2D real Field names 
                                                                                              !<-- in one long string
            NCHAR_R2D(1) = NCHAR_R2D(1) + NAME_MAXSTR
!
          ENDIF
!
        ENDDO
        deallocate(field_name)
!
        wrt_int_state%KOUNT_R2D(NBDL) = KOUNT_R2D
        wrt_int_state%KOUNT_I2D(NBDL) = KOUNT_I2D
!
      ENDIF fcst_tasks
!-----------------------------------------------------------------------
!***  IF THERE ARE NO QUANTITIES SPECIFIED FOR HISTORY OUTPUT,
!***  FORECAST TASK 0 WILL INFORM THE WRITE TASKS AND THEN
!***  EVERYONE WILL RETURN.
!-----------------------------------------------------------------------
!
      NO_FIELDS(1) = ESMF_FALSE
!
!
!-----------------------------------------------------------------------
!
        if1task2:  if(wrt_int_state%quilting) then
!

      IF(MYPE==0)THEN
        IF(wrt_int_state%NCOUNT_FIELDS(1) == 0)NO_FIELDS(1) = ESMF_TRUE    !<-- Reset flag saying there are no history quantities
        LAST_WRITE_TASK = NTASKS-1                                         !<-- The last write task in this group
!
        DO N=LEAD_WRITE_TASK,LAST_WRITE_TASK                               !<-- Loop through all the write tasks in the write group
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK = "Fcst Task0 Informs All That There Are No Fields"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_VMSend(VM                                  &  !<-- ESMF Virtual Machine
                          ,NO_FIELDS                           &  !<-- Send this data
                          ,1                                   &  !<-- Words sent
                          ,N                                   &  !<-- Receiving task in active write group
                          ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        ENDDO
      ENDIF
!
      IF(MYPE >= LEAD_WRITE_TASK)THEN                                      !<-- All write tasks in this group receive
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK = "Write Tasks Told By Fcst Task0 There Are No Fields"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_VMRecv(VM                                    &  !<-- ESMF Virtual Machine
                        ,NO_FIELDS                             &  !<-- Recv this data
                        ,1                                     &  !<-- Words received
                        ,0                                     &  !<-- Sending task (forecast task 0)
                        ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF
!
!
      else   !for 1 task
!
        IF(wrt_int_state%NCOUNT_FIELDS(1)==0)NO_FIELDS(1) = ESMF_TRUE      !<-- Reset flag saying there are no history quantities

      endif  if1task2
!-----------------------------------------------------------------------
!
      IF(NO_FIELDS(1) == ESMF_TRUE) THEN
        IF(MYPE == 0) THEN
          WRITE(6,*)'WARNING: No Import ESMF quantities for the Write Component'
          WRITE(0,*)'WARNING: No Import ESMF quantities for the Write Component'
        ENDIF
!
        RETURN                                                             !<-- All tasks return if there is no history output
!
      ENDIF
!
!
!-----------------------------------------------------------------------
!
        if1task3:   if(wrt_int_state%quilting) then
!-----------------------------------------------------------------------
!***  FORECAST TASK 0 SENDS ALL THE WRITE TASKS THE NUMBER OF
!***  REAL AND INTEGER 2D GRIDDED QUANTITIES PLUS ALL OF THE
!***  LOCAL HORIZONTAL DOMAIN LIMITS IN PREPARATION FOR THE
!***  WRITE TASKS' RECEIVING AND ASSEMBLING THE LOCAL HISTORY
!***  DATA THEY RECEIVE FROM THE FORECAST TASKS.
!-----------------------------------------------------------------------
!
       IF(MYPE == 0)THEN                                                   !<-- Forecast task 0 sends
!
        LAST_WRITE_TASK = NTASKS - 1                                       !<-- The last write task in this group
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK = "Fcst Task0 Sends Write Tasks Info for Quilting"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        DO N=LEAD_WRITE_TASK,LAST_WRITE_TASK                               !<-- Loop through all the write tasks in the write group
!
          ITMP(1)=wrt_int_state%NCOUNT_FIELDS(NBDL)
          CALL ESMF_VMSend(VM                                  &  !<-- ESMF Virtual Machine
                          ,ITMP                                &  !<-- Send this data
                          ,1                                   &  !<-- Words sent
                          ,N                                   &  !<-- Receiving task in active write group
                          ,rc      =RC)
!
          ITMP(1)=wrt_int_state%KOUNT_R2D(NBDL)    
          CALL ESMF_VMSend(VM                                  &  !<-- ESMF Virtual Machine
                          ,ITMP                                &  !<-- Send this data
                          ,1                                   &  !<-- Words sent
                          ,N                                   &  !<-- Receiving task in active write group
                          ,rc      =RC)
!
          ITMP(1)=wrt_int_state%KOUNT_I2D(NBDL)    
          CALL ESMF_VMSend(VM                                  &  !<-- ESMF Virtual Machine
                          ,ITMP                                &  !<-- Send this data
                          ,1                                   &  !<-- Words sent
                          ,N                                   &  !<-- Receiving task in active write group
                          ,rc      =RC)
!
        ENDDO
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
       ENDIF
!
!-----------------------------------------------------------------------
!
       IF(MYPE >= LEAD_WRITE_TASK) THEN                                        !<-- All write tasks in this group receive
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK = "Write Tasks Recv Quilting Info From Fcst Task0"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_VMRecv(VM                                    &  !<-- ESMF Virtual Machine
                        ,ITMP                                  &  !<-- Recv this data
                        ,1                                     &  !<-- Words received
                        ,0                                     &  !<-- Sending task (forecast task 0)
                        ,rc      =RC)
        wrt_int_state%NCOUNT_FIELDS(NBDL) = ITMP(1)
!
        CALL ESMF_VMRecv(VM                                    &  !<-- ESMF Virtual Machine
                        ,ITMP                                  &  !<-- Recv this data
                        ,1                                     &  !<-- Words received
                        ,0                                     &  !<-- Sending task (forecast task 0)
                        ,rc      =RC)
        wrt_int_state%KOUNT_R2D(NBDL) = ITMP(1)
!
        CALL ESMF_VMRecv(VM                                    &  !<-- ESMF Virtual Machine
                        ,ITMP                                  &  !<-- Recv this data
                        ,1                                     &  !<-- Words received
                        ,0                                     &  !<-- Sending task (forecast task 0)
                        ,rc      =RC)
        wrt_int_state%KOUNT_I2D(NBDL) = ITMP(1)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF
!
!-----------------------------------------------------------------------
!***  FORECAST TASK 0 SENDS THE 2D DATA NAMES TO ALL LEAD WRITE TASKS.
!-----------------------------------------------------------------------
!
      IF(MYPE == 0)THEN                                                    !<-- Fcst task0 alone can send write task preliminary info
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK = "Fcst Task0 Sends Write Tasks 2D Data Names"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
       DO N=LEAD_WRITE_TASK,LAST_WRITE_TASK
!
         CALL ESMF_VMSend(VM                                    & !<-- ESMF Virtual Machine
                         ,NCHAR_I2D                             & !<-- Send total length of the names of 2D integer data
                         ,1                                     & !<-- Words sent
                         ,N                                     & !<-- Receiving task (1st write task in group)
                         ,rc      =RC)
!
        if(NCHAR_I2D(1) > 0) then
          CALL ESMF_VMSend(VM                                   & !<-- ESMF Virtual Machine
                          ,wrt_int_state%NAMES_I2D_STRING(nbdl) & !<-- Send names of 2D integer history variables
                          ,NCHAR_I2D(1)                         & !<-- Words sent
                          ,N                                    & !<-- Receiving task (1st write task in group)
                          ,rc      =RC)
        endif
!
         CALL ESMF_VMSend(VM                                   &  !<-- ESMF Virtual Machine
                         ,NCHAR_R2D                            &  !<-- Send total length of the names of 2D real data
                         ,1                                    &  !<-- Words sent
                         ,N                                    &  !<-- Receiving task (1st write task in group)
                         ,rc      =RC)
!
        if(NCHAR_R2D(1) > 0) then
          CALL ESMF_VMSend(VM                                  &  !<-- ESMF Virtual Machine
                          ,wrt_int_state%NAMES_R2D_STRING(NBDL)&  !<-- Send names of 2D real history variables
                          ,NCHAR_R2D(1)                        &  !<-- Words sent
                          ,N                                   &  !<-- Receiving task (1st write task in group)
                          ,rc      =RC)
         endif
!
       ENDDO
!
      ELSEIF(MYPE >= LEAD_WRITE_TASK) then                        !<-- 1st write task receives 2D preliminary info
!
        CALL ESMF_VMRecv(VM                                    &  !<-- ESMF Virtual Machine
                        ,NCHAR_I2D                             &  !<-- Recv total length of the names of 2D integer data
                        ,1                                     &  !<-- Words sent
                        ,0                                     &  !<-- Sending task (fcst task 0)
                        ,rc      =RC)
!
        if(NCHAR_I2D(1)>0) then
        CALL ESMF_VMRecv(VM                                    &  !<-- ESMF Virtual Machine
                        ,wrt_int_state%NAMES_I2D_STRING(NBDL)        &  !<-- Recv names of 2D integer history variables
                        ,NCHAR_I2D(1)                          &  !<-- Words sent
                        ,0                                     &  !<-- Sending task (fcst task 0)
                        ,rc      =RC)
         endif
!
        CALL ESMF_VMRecv(VM                                    &  !<-- ESMF Virtual Machine
                        ,NCHAR_R2D                             &  !<-- Recv total length of the names of 2D gridded data
                        ,1                                     &  !<-- Words sent
                        ,0                                     &  !<-- Sending task (fcst task 0)
                        ,rc      =RC)
!
        if(NCHAR_R2D(1) > 0) then
        CALL ESMF_VMRecv(VM                                    &  !<-- ESMF Virtual Machine
                        ,wrt_int_state%NAMES_R2D_STRING(NBDL)        &  !<-- Recv names of 2D real history variables
                        ,NCHAR_R2D(1)                          &  !<-- Words sent
                        ,0                                     &  !<-- Sending task (fcst task 0)
                        ,rc      =RC)
         endif
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF  !write task
!
!-----------------------------------------------------------------------
!***  EACH WRITE TASK MUST KNOW THE IDs OF THE FORECAST TASKS
!***  FROM WHICH IT WILL RECEIVE 2D GRIDDED HISTORY DATA.
!########################################################################
!
      IF(MYPE >= LEAD_WRITE_TASK)THEN                                      !<-- The write tasks
!
!-----------------------------------------------------------------------
!***  EACH WRITE TASK ALSO MUST KNOW THE NORTH-SOUTH EXTENT OF THE
!***  FULL 2D DOMAIN THAT IT WILL HANDLE.  THIS IS DETERMINED BY
!***  THE COVERAGE OF THE FCST TASKS THAT SEND TO IT.
!-----------------------------------------------------------------------
!
        JSTA_WRITE = wrt_int_state%JSTART_WRITE(mype_local+1)  !<-- JTS of 1st fcst task that sends to this write task
        JEND_WRITE = wrt_int_state%JEND_WRITE(mype_local+1)    !<-- JTE of last fcst task that sends to this write task
!
!-----------------------------------------------------------------------
!***  NOW EACH WRITE TASK ALLOCATES ITS OWN SECTION OF THE 2D DOMAIN
!***  FOR ALL THE 2D VARIABLES IT WILL RECEIVE AND ITS 1D EQUIVALENT
!***  USED TO TRANSFER THE DATA TO THE LEAD WRITE TASK.
!-----------------------------------------------------------------------
!
       if(.not. allocated(wrt_int_state%BUFF_INT)) then
        LENGTH = IM*JM
        ALLOCATE(wrt_int_state%BUFF_INT(LENGTH))
        ALLOCATE(wrt_int_state%BUFF_INT_TMP(LENGTH))
        ALLOCATE(wrt_int_state%BUFF_REAL(LENGTH))
        ALLOCATE(wrt_int_state%BUFF_REAL_TMP(LENGTH))
       endif
!
       ENDIF  !MYPE >= LEAD_WRITE_TASK
!
!-----------------------------------------------------------------------
!
      endif if1task3
!
!-----------------------------------------------------------------------
!***  THE LEAD WRITE TASK ALLOCATES ITS WORKING ARRAYS INTO WHICH
!***  IT WILL ASSEMBLE EACH INDIVIDUAL 2D FIELD THAT WILL BE
!***  WRITTEN TO THE HISTORY FILES.
!-----------------------------------------------------------------------
!
     IF(MYPE >= LEAD_WRITE_TASK .and. MYPE_LOCAL == 0 )THEN
       if(.not. allocated(wrt_int_state%OUTPUT_ARRAY_I2D)) then
        ALLOCATE(wrt_int_state%OUTPUT_ARRAY_I2D(1:IM,1:JM))
        ALLOCATE(wrt_int_state%OUTPUT_ARRAY_R2D(1:IM,1:JM))
       endif
     ENDIF
! 
!-----------------------------------------------------------------------
!
     if1task4:  if(wrt_int_state%quilting) then
!
!-----------------------------------------------------------------------
!***  SINCE ALL SCALAR/1D DATA IS IDENTICAL ON ALL FORECAST TASKS,
!***  TASK 0 ALONE CAN SEND THE INFORMATION TO THE LEAD WRITE TASK
!***  THAT WILL LATER WRITE IT TO THE HISTORY FILE.
!--------------------------------------------------------------------
!
!--------------------------------------------------------------------
      task_0_sends: IF(MYPE == 0) THEN                                    !<-- Forecast task 0 sends
!--------------------------------------------------------------------
!
!------------------------------------------------
!***  SEND SCALAR/1D INTEGER HISTORY INFORMATION. to all the write tasks
!------------------------------------------------
!
       DO N=LEAD_WRITE_TASK,LAST_WRITE_TASK

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK = "Fcst Task0 Sends Scalar/1D Integer History Data"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        ITMP(1) = wrt_int_state%KOUNT_I1D(NBDL)
        CALL ESMF_VMSend(VM                                    &  !<-- ESMF Virtual Machine
                        ,ITMP                                  &  !<-- Send # of scalar/1D integer history variables
                        ,1                                     &  !<-- Words sent
                        ,N                                     &  !<-- Receiving task (1st write task in group)
                        ,rc      =RC)
!
!if there are any I 1D data
!---------------------------------------------------------------------
        if(wrt_int_state%KOUNT_I1D(NBDL) > 0 )then
!
        ITMP(1) = wrt_int_state%LENGTH_SUM_I1D(NBDL)
        CALL ESMF_VMSend(VM                                    &  !<-- ESMF Virtual Machine
                        ,ITMP                                  &  !<-- Send length of string of all such integer history variables
                        ,1                                     &  !<-- Words sent
                        ,N                                     &  !<-- Receiving task (1st write task in group)
                        ,rc      =RC)
!
        ITMP(1:wrt_int_state%KOUNT_I1D(NBDL)) =                &
         wrt_int_state%LENGTH_DATA_I1D(1:wrt_int_state%KOUNT_I1D(NBDL),NBDL)
        CALL ESMF_VMSend(VM                                    &  !<-- ESMF Virtual Machine
                        ,ITMP                                  &  !<-- Send lengths of each scalar/1D integer history variable
                        ,wrt_int_state%KOUNT_I1D(NBDL)         &  !<-- Words sent
                        ,N                                     &  !<-- Receiving task (1st write task in group)
                        ,rc      =RC)
!
        NAMETMP(1:NCHAR_I1D)=wrt_int_state%NAMES_I1D_STRING(nbdl)(1:NCHAR_I1D)
        CALL ESMF_VMSend(VM                                    &  !<-- ESMF Virtual Machine
                        ,NAMETMP                               &  !<-- Send names of each scalar/1D integer history variable
                        ,NCHAR_I1D                             &  !<-- Words sent
                        ,N                                     &  !<-- Receiving task (1st write task in group)
                        ,rc      =RC)
!
        ITMP(1:wrt_int_state%LENGTH_SUM_I1D(NBDL)) =           &
           wrt_int_state%ALL_DATA_I1D(1:wrt_int_state%LENGTH_SUM_I1D(NBDL),NBDL)
        CALL ESMF_VMSend(VM                                    &  !<-- ESMF Virtual Machine
                        ,ITMP                                  &  !<-- Send the full string of all scalar/1D integer history data
                        ,wrt_int_state%LENGTH_SUM_I1D(NBDL)    &  !<-- Words sent
                        ,N                                     &  !<-- Receiving task (1st write task in group)
                        ,rc      =RC)
!
         endif
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!---------------------------------------------
!***  SEND SCALAR/1D REAL HISTORY INFORMATION.
!---------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK = "Fcst Task0 Sends Scalar/1D Real History Data"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        ITMP(1) = wrt_int_state%KOUNT_R1D(NBDL)
        CALL ESMF_VMSend(VM                                    &  !<-- ESMF Virtual Machine
                        ,ITMP                                  &  !<-- Send # of scalar/1D real history variables
                        ,1                                     &  !<-- Words sent
                        ,N                                     &  !<-- Receiving task (1st write task in group)
                        ,rc      =RC)
!
!
!if there are any R 1D data
!---------------------------------------------------------------------
       if(wrt_int_state%KOUNT_R1D(NBDL) > 0 )then
!
        ITMP(1) = wrt_int_state%LENGTH_SUM_R1D(NBDL)
        CALL ESMF_VMSend(VM                                    &  !<-- ESMF Virtual Machine
                        ,ITMP                                  &  !<-- Send length of string of all such real history variables
                        ,1                                     &  !<-- Words sent
                        ,N                                     &  !<-- Receiving task (1st write task in group)
                        ,rc      =RC)
!
        ITMP(1:wrt_int_state%KOUNT_R1D(NBDL)) = wrt_int_state%LENGTH_DATA_R1D &
            (1:wrt_int_state%KOUNT_R1D(NBDL),NBDL)
        CALL ESMF_VMSend(VM                                    &  !<-- ESMF Virtual Machine
                        ,ITMP                                  &!<-- Send lengths of each scalar/1D real history variable
                        ,wrt_int_state%KOUNT_R1D(NBDL)         &  !<-- Words sent
                        ,N                                     &  !<-- Receiving task (1st write task in group)
                        ,rc      =RC)
!
        NAMETMP(1:NCHAR_R1D) = wrt_int_state%NAMES_R1D_STRING(NBDL)(1:NCHAR_R1D)
        CALL ESMF_VMSend(VM                                    &  !<-- ESMF Virtual Machine
                        ,NAMETMP                               &  !<-- Send names of each scalar/1D real history variable
                        ,NCHAR_R1D                             &  !<-- Words sent
                        ,N                                     &  !<-- Receiving task (1st write task in group)
                        ,rc      =RC)
!
        RTMP(1:wrt_int_state%LENGTH_SUM_R1D(NBDL)) = wrt_int_state%ALL_DATA_R1D &
            (1:wrt_int_state%LENGTH_SUM_R1D(NBDL),NBDL)
        CALL ESMF_VMSend(VM                                    &  !<-- ESMF Virtual Machine
                        ,RTMP                                  &  !<-- Send the full string of all scalar/1D real history data
                        ,wrt_int_state%LENGTH_SUM_R1D(NBDL)    &  !<-- Words sent
                        ,N                                     &  !<-- Receiving task (1st write task in group)
                        ,rc      =RC)
!
         ENDIF
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!--------------------------------------
!***  SEND LOGICAL HISTORY INFORMATION.
!--------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK = "Fcst Task0 Sends Logical History Data"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        ITMP(1)=wrt_int_state%KOUNT_LOG(NBDL)
        CALL ESMF_VMSend(VM                                    &  !<-- ESMF Virtual Machine
                        ,ITMP                                  &  !<-- Send # of logical history variables
                        ,1                                     &  !<-- Words sent
                        ,N                                     &  !<-- Receiving task (1st write task in group)
                        ,rc      =RC)
!
!
!if there are any R 1D data
!---------------------------------------------------------------------
        if(wrt_int_state%KOUNT_LOG(NBDL) > 0 )then
!
        ITMP(1) = wrt_int_state%LENGTH_SUM_LOG(NBDL)
        CALL ESMF_VMSend(VM                                    &  !<-- ESMF Virtual Machine
                        ,ITMP                                  &  !<-- Send length of string of all logical variables
                        ,1                                     &  !<-- Words sent
                        ,N                                     &  !<-- Receiving task (1st write task in group)
                        ,rc      =RC)
!
        NAMETMP(1:NCHAR_R1D) = wrt_int_state%NAMES_LOG_STRING(NBDL)(1:NCHAR_LOG)
        CALL ESMF_VMSend(VM                                    &  !<-- ESMF Virtual Machine
                        ,NAMETMP                               & !<-- Send names of each logical history variable
                        ,NCHAR_LOG                             &  !<-- Words sent
                        ,N                                     &  !<-- Receiving task (1st write task in group)
                        ,rc      =RC)
!
        DO I=1,wrt_int_state%LENGTH_SUM_LOG(NBDL)
          if(wrt_int_state%ALL_DATA_LOG(i,NBDL) == ESMF_TRUE) then
            LTMP(I) = ESMF_TRUE
          else
            LTMP(I) = ESMF_FALSE
          endif
        ENDDO

        CALL ESMF_VMSend(VM                                    &  !<-- ESMF Virtual Machine
                        ,LTMP                                  &  !<-- Send the full string of all logical history data
                        ,wrt_int_state%LENGTH_SUM_LOG(NBDL)    &  !<-- Words sent
                        ,N                                     &  !<-- Receiving task (1st write task in group)
                        ,rc      =RC)
!
        ENDIF
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDDO
!
!-----------------------------------------------------------------------
      ENDIF task_0_sends
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
      write_task_recvs: IF(MYPE >= LEAD_WRITE_TASK) THEN                   !<-- write tasks in this group receives
                                                                           !    all of the data just sent to it by
                                                                           !    fcst task 0
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  RECEIVE SCALAR/1D INTEGER HISTORY INFORMATION
!***  FROM FORECAST TASK 0.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Write Tasks Recv Scalar/1D Integer History Data"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_VMRecv(VM                                    &  !<-- ESMF Virtual Machine
                        ,ITMP                                  &  !<-- Recv # of integer history variables
                        ,1                                     &  !<-- Words received
                        ,0                                     &  !<-- Sending task (forecast task 0)
                        ,rc      =RC)
        wrt_int_state%KOUNT_I1D(NBDL)=ITMP(1)
!
!if there are any R 1D data
!---------------------------------------------------------------------
        if(wrt_int_state%KOUNT_I1D(NBDL) > 0 )then
!
         CALL ESMF_VMRecv(VM                                   &  !<-- ESMF Virtual Machine
                        ,itmp                                  &  !<-- Recv length of string of all integer history variables
                        ,1                                     &  !<-- Words received
                        ,0                                     &  !<-- Sending task (forecast task 0)
                        ,rc      =RC)
         wrt_int_state%LENGTH_SUM_I1D(NBDL) = itmp(1)
!
         CALL ESMF_VMRecv(VM                                   &  !<-- ESMF Virtual Machine
                        ,itmp                                  &  !<-- Recv lengths of each integer history variable
                        ,wrt_int_state%KOUNT_I1D(NBDL)         &  !<-- Words received
                        ,0                                     &  !<-- Sending task (forecast task 0)
                        ,rc      =RC)
        wrt_int_state%LENGTH_DATA_I1D(1:wrt_int_state%KOUNT_I1D(NBDL),NBDL) = &
                                 itmp(1:wrt_int_state%KOUNT_I1D(NBDL))
!
        NCHAR_I1D = wrt_int_state%KOUNT_I1D(NBDL)*NAME_MAXSTR
        CALL ESMF_VMRecv(VM                                    &  !<-- ESMF Virtual Machine
                        ,NAMETMP                               &  !<-- Recv names of integer history variables
                        ,NCHAR_I1D                             &  !<-- Words received
                        ,0                                     &  !<-- Sending task (forecast task 0)
                        ,rc      =RC)
        wrt_int_state%NAMES_I1D_STRING(NBDL) = trim(NAMETMP)
!
        CALL ESMF_VMRecv(VM                                    &  !<-- ESMF Virtual Machine
                        ,ITMP                                  &  !<-- Recv the string of integer history data
                        ,wrt_int_state%LENGTH_SUM_I1D(NBDL)    &  !<-- Words received
                        ,0                                     &  !<-- Sending task (forecast task 0)
                        ,rc      =RC)
        wrt_int_state%ALL_DATA_I1D(1:wrt_int_state%LENGTH_SUM_I1D(NBDL),NBDL)= &
          ITMP(1:wrt_int_state%LENGTH_SUM_I1D(NBDL))
!
        ENDIF
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  RECEIVE SCALAR/1D REAL HISTORY INFORMATION.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK = "Write Tasks Recv Scalar/1D Real History Data"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_VMRecv(VM                                    &  !<-- ESMF Virtual Machine
                        ,itmp                                  &  !<-- Recv # of scalar/1D real history variables
                        ,1                                     &  !<-- Words received
                        ,0                                     &  !<-- Sending task (forecast task 0)
                        ,rc      =RC)
        wrt_int_state%KOUNT_R1D(NBDL)=itmp(1)
!
!if there are any R 1D data
!---------------------------------------------------------------------
        if(wrt_int_state%KOUNT_R1D(NBDL) > 0 )then
!
        CALL ESMF_VMRecv(VM                                    &  !<-- ESMF Virtual Machine
                        ,itmp                                  &  !<-- Recv length of string of all such real history variables
                        ,1                                     &  !<-- Words received
                        ,0                                     &  !<-- Sending task (forecast task 0)
                        ,rc      =RC)
        wrt_int_state%LENGTH_SUM_R1D(NBDL) = itmp(1)
!
        CALL ESMF_VMRecv(VM                                    &  !<-- ESMF Virtual Machine
                        ,itmp                                  &  !<-- Recv lengths of each scalar/1D real history variable
                        ,wrt_int_state%KOUNT_R1D(NBDL)         &  !<-- Words received
                        ,0                                     &  !<-- Sending task (forecast task 0)
                        ,rc      =RC)
        wrt_int_state%LENGTH_DATA_R1D(1:wrt_int_state%KOUNT_R1D(NBDL),NBDL) = &
                                 itmp(1:wrt_int_state%KOUNT_R1D(NBDL))
!
        NCHAR_R1D = wrt_int_state%KOUNT_R1D(NBDL)*ESMF_MAXSTR
        CALL ESMF_VMRecv(VM                                    &  !<-- ESMF Virtual Machine
                        ,nametmp                               &  !<-- Recv names of scalar/1D real history variables
                        ,NCHAR_R1D                             &  !<-- Words received
                        ,0                                     &  !<-- Sending task (forecast task 0)
                        ,rc      =RC)
        wrt_int_state%NAMES_R1D_STRING(NBDL) = trim(nametmp)
!
        CALL ESMF_VMRecv(VM                                    &  !<-- ESMF Virtual Machine
                        ,rtmp                                  &  !<-- Recv the string of all scalar/1D real history data
                        ,wrt_int_state%LENGTH_SUM_R1D(NBDL)    &  !<-- Words received
                        ,0                                     &  !<-- Sending task (forecast task 0)
                        ,rc      =RC)
        wrt_int_state%ALL_DATA_R1D(1:wrt_int_state%LENGTH_SUM_R1D(NBDL),NBDL) = &
                              rtmp(1:wrt_int_state%LENGTH_SUM_R1D(NBDL))
!
        ENDIF
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  RECEIVE LOGICAL HISTORY INFORMATION.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK = "Write Tasks Recv Logical Real History Data"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_VMRecv(VM                                    &  !<-- ESMF Virtual Machine
                        ,itmp                                  &  !<-- Recv # of logical history variables
                        ,1                                     &  !<-- Words received
                        ,0                                     &  !<-- Sending task (forecast task 0)
                        ,rc      =RC)
        wrt_int_state%KOUNT_LOG(NBDL) = itmp(1)
!
!if there are any R 1D data
!---------------------------------------------------------------------
        if(wrt_int_state%KOUNT_LOG(NBDL)>0 )then

!
        CALL ESMF_VMRecv(VM                                    &  !<-- ESMF Virtual Machine
                        ,itmp                                  &  !<-- Recv length of string of all logical history variables
                        ,1                                     &  !<-- Words received
                        ,0                                     &  !<-- Sending task (forecast task 0)
                        ,rc      =RC)
        wrt_int_state%LENGTH_SUM_LOG(NBDL) = itmp(1)
!
        NCHAR_LOG=wrt_int_state%KOUNT_LOG(NBDL)*NAME_MAXSTR
        CALL ESMF_VMRecv(VM                                    &  !<-- ESMF Virtual Machine
                        ,NAMETMP                               &  !<-- Recv names of logical history variables
                        ,NCHAR_LOG                             &  !<-- Words received
                        ,0                                     &  !<-- Sending task (forecast task 0)
                        ,rc      =RC)
        wrt_int_state%NAMES_LOG_STRING(NBDL) = trim(NAMETMP)
!
        CALL ESMF_VMRecv(VM                                    &  !<-- ESMF Virtual Machine
                        ,LTMP                                  &  !<-- Recv the string of all logical history data
                        ,wrt_int_state%LENGTH_SUM_LOG(NBDL)    &  !<-- Words received
                        ,0                                     &  !<-- Sending task (forecast task 0)
                        ,rc      =RC)
        wrt_int_state%ALL_DATA_LOG(1:wrt_int_state%LENGTH_SUM_LOG(NBDL),NBDL) = &
                              LTMP(1:wrt_int_state%LENGTH_SUM_LOG(NBDL))
        ENDIF

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
      ENDIF write_task_recvs
!      write(0,*)'before first_pass_gfs end,KOUNT_I1D=',wrt_int_state%KOUNT_I1D(NBDL), &
!       'I1D length=',wrt_int_state%LENGTH_DATA_I1D(1:wrt_int_state%KOUNT_I1D(NBDL),NBDL), &
!       'I1D value=',wrt_int_state%ALL_DATA_I1D(1:wrt_int_state%KOUNT_I1D(NBDL),NBDL), &
!       'KOUNT_R1D=',wrt_int_state%KOUNT_R1D(NBDL),'KOUNT_LOG=', &
!       wrt_int_state%KOUNT_LOG(NBDL)
!
!-----------------------------------------------------------------------
!
      endif  if1task4
!
!-----------------------------------------------------------------------
!
      DEALLOCATE(ITMP)
      DEALLOCATE(RTMP)
      DEALLOCATE(LTMP)
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE FIRST_PASS_GFS
!
!-----------------------------------------------------------------------
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!-----------------------------------------------------------------------
!
      SUBROUTINE WRITE_NEMSIO_OPEN(WRT_INT_STATE                  &
                                  ,NBDL                           &
                                  ,NEMSIOFILE                     &
                                  ,IYEAR_FCST                     &
                                  ,IMONTH_FCST                    &
                                  ,IDAY_FCST                      &
                                  ,IHOUR_FCST                     &
                                  ,IMINUTE_FCST                   &
                                  ,SECOND_FCST                    &
                                  ,NF_HOURS                       &
                                  ,NF_MINUTES                     &
                                  ,NF_SECONDS                     &
                                  ,NF_HOURS_IAU                   &
                                  ,DIM1,DIM2,NFRAME               &
                                  ,LEAD_WRITE_TASK)

!-----------------------------------------------------------------------
!***  WRITE OUT A NEMSIO BINARY RUN HISTORY FILE.
!-----------------------------------------------------------------------
!
      TYPE(WRITE_INTERNAL_STATE_GFS),INTENT(INOUT) :: WRT_INT_STATE             !<-- The Write component's internal state
!
      TYPE(NEMSIO_GFILE),INTENT(INOUT)             :: NEMSIOFILE                !<-- The nemsio file handler
!
      INTEGER,INTENT(IN)  :: IYEAR_FCST, IMONTH_FCST, IDAY_FCST, IHOUR_FCST  &
                            ,IMINUTE_FCST, NF_HOURS, NF_MINUTES,NF_HOURS_IAU &
                            ,LEAD_WRITE_TASK, NBDL

      INTEGER,INTENT(OUT) :: DIM1,DIM2,NFRAME
!
      REAL,INTENT(IN)     :: NF_SECONDS, SECOND_FCST
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      INTEGER,DIMENSION(:),    POINTER :: ARYILEN, ARYRLEN, RECLEV, VARIVAL
!
      INTEGER,DIMENSION(:,:),  POINTER :: ARYIVAL
!
      REAL(4),DIMENSION(:,:,:),POINTER :: VCOORD
!
!jw
      REAL(KIND=kind_io4),DIMENSION(:)  ,POINTER :: VARRVAL,CPI,RI
      REAL(KIND=kind_io4),DIMENSION(:,:),POINTER :: ARYRVAL
!
      LOGICAL                      :: LFNHR, GLOBAL, HYBRID, GEN_COORD_HYBRID
      LOGICAL,DIMENSION(:),POINTER :: VARLVAL
!
      CHARACTER(40)                       :: CFHOUR,CFORM
      CHARACTER(16)                       :: VLEVTYP
!
      CHARACTER(16),DIMENSION(:) ,POINTER :: ARYINAME => null()          &
                                            ,ARYRNAME => null()          &
                                            ,RECNAME  => null()          &
                                            ,VARINAME => null()          &
                                            ,VARRNAME => null()          &
                                            ,VARLNAME => null()
!
      CHARACTER(16),DIMENSION(:),POINTER :: RECLEVTYP
!
      CHARACTER(ESMF_MAXSTR) :: NAME,FILENAME
!
      TYPE(ESMF_Logical)     :: WORK_LOGICAL
!
      INTEGER                :: NREC=0    
!
      INTEGER :: NAK5,NBK5,NCK5,NSI,NIDVC                               &
                ,NIDVM,NIDSL,NIDRT,Nthermodyn_id                        &
                ,Nsfcpress_id,Nvertcoord_id,NTRAC,NCLD                  &
                ,NGEN_COORD_HYBRID,NHYBRID                              &
                ,FIELDSIZE,IM,JM,LM,IDATE(7),FCSTDATE(7)                &
                ,INDX_2D,INDX_2D2,INDX_2D3,INDX_2DA                     &
                ,IRET,IND1,IND2,IND3,IND4,CNT,INI1,INI2                 &
                ,N2ISCALAR,N2IARY,N2RSCALAR,N2RARY,N2LSCALAR            &
                ,NMETA,TLMETA,VLEV,LSOIL,NFIELD,NDIG,RC                 &
                ,I,J,N,N1,N2,NPOSN_1,NPOSN_2,LENGTH,MAXLENGTH           &
                ,JCAP,IDVC,IDVM,IDSL,IDRT,NVCOORD,NSOIL
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!initialization
!
      NAK5  = 0  ; NBK5  = 0  ; NCK5  = 0  ; NSI = 0  ; NIDVC = 0
      NIDVM = 0  ; NIDSL = 0  ; NIDRT = 0  ; Nthermodyn_id = 0
      Nsfcpress_id = 0  ; Nvertcoord_id = 0  ; NTRAC = 0  ; NCLD = 0
      NGEN_COORD_HYBRID = 0 ; NHYBRID = 0
!
      JCAP  = -9999 ; IDVC  = -9999  ; IDVM = -9999 ; IDSL =-9999
      IDRT = -9999  ; NSOIL = -9999  ; NVCOORD = -9999
      GEN_COORD_HYBRID = .false. ; HYBRID = .false.
!
      FCSTDATE(1) = IYEAR_FCST
      FCSTDATE(2) = IMONTH_FCST
      FCSTDATE(3) = IDAY_FCST
      FCSTDATE(4) = IHOUR_FCST
      FCSTDATE(5) = IMINUTE_FCST
      FCSTDATE(6) = nint(SECOND_FCST*100.)
      FCSTDATE(7) = 100

!      write(0,*)'fcstdate=',fcstdate(1:6),'kount_I1d=',wrt_int_state%KOUNT_I1D(NBDL), &
!       'KOUNT_R1D=',wrt_int_state%KOUNT_R1D(NBDL),'LOG=',wrt_int_state%KOUNT_LOG(NBDL),'NAK5=',NAK5
!
!-----------------------------------------------------------------------
!***  INTEGER SCALAR/1D HISTORY VARIABLES
!-----------------------------------------------------------------------
!
!-------------------------------------------------------------
!*** Find out the total number of int scalars and int arrays.
!-------------------------------------------------------------
!
      N2ISCALAR = 0
      N2IARY    = 0
      MAXLENGTH = 1
!
      DO N=1,wrt_int_state%KOUNT_I1D(NBDL)                                    !<-- Loop through all scalar/1D integer data
        LENGTH = wrt_int_state%LENGTH_DATA_I1D(N,NBDL)
!
        IF(LENGTH == 1) THEN
          N2ISCALAR = N2ISCALAR + 1
        ELSE
          N2IARY    = N2IARY + 1
          MAXLENGTH = MAX(LENGTH,MAXLENGTH)
        ENDIF
!
      ENDDO
!
      N2IARY    = N2IARY + 1
      MAXLENGTH = MAX(MAXLENGTH,7)
      ALLOCATE(VARINAME(N2ISCALAR),VARIVAL(N2ISCALAR))
      ALLOCATE(ARYINAME(N2IARY),ARYILEN(N2IARY),ARYIVAL(MAXLENGTH,N2IARY))
!
!---------------------------------------
!***  SET VALUE TO AVRIVAL and ARYIVAL.
!---------------------------------------
!
      N2 = 0                                                             !<-- Word counter for full string of integer scalar/1D data
      N2ISCALAR = 0
      N2IARY    = 0
      IDATE     = 0
      ARYIVAL   = 0
!
      DO N=1,wrt_int_state%KOUNT_I1D(NBDL)                                  !<-- Loop through all scalar/1D integer data
!
        NPOSN_1 = (N-1)*NAME_MAXSTR + 1
        NPOSN_2 = N*NAME_MAXSTR
        NAME    = wrt_int_state%NAMES_I1D_STRING(NBDL)(NPOSN_1:NPOSN_2)      !<-- The variable's name
        LENGTH  = wrt_int_state%LENGTH_DATA_I1D(N,NBDL)                      !<-- The variable's length in words
!
        IF(LENGTH == 1) THEN
          N2 = N2 + 1
          N2ISCALAR = N2ISCALAR + 1
          VARINAME(N2ISCALAR) = TRIM(NAME)
          VARIVAL(N2ISCALAR)  = wrt_int_state%ALL_DATA_I1D(N2,NBDL)
          IF(VARINAME(N2ISCALAR) == 'IHRST') then
            IDATE(4) = VARIVAL(N2ISCALAR)
          ELSEIF(VARINAME(N2ISCALAR) == 'jcap') then
            JCAP     = VARIVAL(N2ISCALAR)
          ELSEIF(VARINAME(N2ISCALAR) == 'idrt') then
            IDRT     = VARIVAL(N2ISCALAR)
          ELSEIF(VARINAME(N2ISCALAR) == 'idsl') then
            IDSL     = VARIVAL(N2ISCALAR)
          ELSEIF(VARINAME(N2ISCALAR) == 'idvc') then
            IDVC     = VARIVAL(N2ISCALAR)
          ELSEIF(VARINAME(N2ISCALAR) == 'sfcpress_id') then
            Nsfcpress_id  = VARIVAL(N2ISCALAR)
          ELSEIF(VARINAME(N2ISCALAR) == 'vertcoord_id') then
            Nvertcoord_id = VARIVAL(N2ISCALAR)
          ELSEIF(VARINAME(N2ISCALAR) == 'thermodyn_id') then
            Nthermodyn_id = VARIVAL(N2ISCALAR)
          ELSEIF(VARINAME(N2ISCALAR) == 'ntrac') then
            NTRAC = VARIVAL(N2ISCALAR)
          ELSEIF(VARINAME(N2ISCALAR) == 'ncld') then
            NCLD = VARIVAL(N2ISCALAR)
          ELSEIF(VARINAME(N2ISCALAR) == 'nsoil') then
            NSOIL = VARIVAL(N2ISCALAR)

          ENDIF
        ELSE
          N2IARY           = N2IARY + 1
          ARYINAME(N2IARY) = TRIM(NAME)
          ARYILEN(N2IARY)  = LENGTH
!            write(0,*)'in I1D array,aryiname=',aryiname(N2IARY),'len=',aryilen(N2IARY),  &
!              wrt_int_state%ALL_DATA_I1D(N2+1:N2+length)
!
          DO N1=1,LENGTH
            N2 = N2 + 1
            ARYIVAL(N1,N2IARY) = wrt_int_state%ALL_DATA_I1D(N2,NBDL)            !<-- Extract the individual data from the data string
          ENDDO

          IF(ARYINAME(N2IARY) == 'IDAT') THEN
            IDATE(4) = ARYIVAL(1,N2IARY)
            IDATE(2) = ARYIVAL(2,N2IARY)
            IDATE(3) = ARYIVAL(3,N2IARY)
            IDATE(1) = ARYIVAL(4,N2IARY)
          ENDIF
          IDATE(7) = 100.
!
        ENDIF
!
      ENDDO
!adjust idrt to default value 4
      if(IDRT == -9999) IDRT = 4
!
!-----------------------------
!***  Add fcst_date into ARYI
!-----------------------------
!
      N2IARY = N2IARY + 1
      ARYINAME(N2IARY) = 'FCSTDATE'
      ARYILEN(N2IARY)  = 7
      ARYIVAL(1:7,N2IARY) = FCSTDATE(1:7)
!
!*** prepare cpi,ri
      allocate(cpi(ntrac+1),ri(ntrac+1))
      cpi = 0.
      ri  = 0.
!
!-----------------------------------------------------------------------
!***  REAL SCALAR/1D HISTORY VARIABLES
!-----------------------------------------------------------------------
!
!------------------------------------------------------------
!***  Find the total number of real scalars and real arrays.
!------------------------------------------------------------
!
      N2RSCALAR = 0
      N2RARY    = 0
      MAXLENGTH = 1
!
      DO N=1,wrt_int_state%KOUNT_R1D(NBDL)                                   !<-- Loop through all scalar/1D real data
        LENGTH = wrt_int_state%LENGTH_DATA_R1D(N,NBDL)                       !<-- The variable's length
        IF(LENGTH == 1)THEN
           N2RSCALAR = N2RSCALAR + 1
        ELSE
          N2RARY     = N2RARY + 1
          MAXLENGTH  = MAX(LENGTH,MAXLENGTH)
        ENDIF
      ENDDO
!      write(0,*)'NBDL=',NBDL,'N2RSCALAR=',N2RSCALAR,'N2RARY=',N2RARY, &
!        'pdryini=',wrt_int_state%pdryini,'R1D=',wrt_int_state%KOUNT_R1D(NBDL)
!
!add pdryini to real scalar

      IF(N2RSCALAR > 0) then
        ALLOCATE(VARRNAME(N2RSCALAR),VARRVAL(N2RSCALAR))
      else
        ALLOCATE(VARRNAME(1),VARRVAL(1))
      endif
      IF(N2RARY    > 0) then
        ALLOCATE(ARYRNAME(N2RARY),ARYRLEN(N2RARY),ARYRVAL(MAXLENGTH,N2RARY))
      else
        ALLOCATE(ARYRNAME(1),ARYRLEN(1),ARYRVAL(1,1))
      endif
!
!------------------------------------------------------
!***  Set values for the real scalars and real arrays.
!------------------------------------------------------
      N2        = 0                                                           !<-- Word counter for full string of real scalar/1D data
      N2RSCALAR = 0
      N2RARY    = 0
!
      DO N=1,wrt_int_state%KOUNT_R1D(NBDL)                                    !<-- Loop through all scalar/1D real data
!
        NPOSN_1 = (N-1)*NAME_MAXSTR+1
        NPOSN_2 = N*NAME_MAXSTR
        NAME    = wrt_int_state%NAMES_R1D_STRING(NBDL)(NPOSN_1:NPOSN_2)       !<-- The variable's name
        LENGTH  = wrt_int_state%LENGTH_DATA_R1D(N,NBDL)                       !<-- The variable's length
!
        IF(LENGTH==1)THEN
          N2 = N2 + 1
          N2RSCALAR = N2RSCALAR + 1
          VARRNAME(N2RSCALAR) = TRIM(NAME)
          VARRVAL(N2RSCALAR)  = wrt_int_state%ALL_DATA_R1D(N2,NBDL)
          if(trim(VARRNAME(N2RSCALAR)) == 'pdryini') then
            VARRVAL(N2RSCALAR) = wrt_int_state%pdryini
          else if(trim(VARRNAME(N2RSCALAR)) == 'zhour') then
            VARRVAL(N2RSCALAR) = wrt_int_state%zhour
          endif
!
        ELSE
!
          N2RARY = N2RARY + 1
          ARYRNAME(N2RARY) = TRIM(NAME)
          ARYRLEN(N2RARY)  = LENGTH
!
          DO N1=1,LENGTH
            N2 = N2 + 1
            ARYRVAL(N1,N2RARY) = wrt_int_state%ALL_DATA_R1D(N2,NBDL)          !<-- Extract the individual data from the data string
          ENDDO
!
          IF( TRIM(NAME) == 'AK5') THEN
            NAK5 = N2RARY
          ELSEIF (TRIM(NAME) == 'BK5') THEN
            NBK5 = N2RARY
          ELSEIF (TRIM(NAME) == 'CK5') THEN
            NCK5 = N2RARY
          ELSEIF (TRIM(NAME) == 'SI') THEN
            NSI = N2RARY
          ELSEIF (TRIM(NAME) == 'CPI') THEN
            cpi(1:length) = aryrval(1:length,n2rary)
          ELSEIF (TRIM(NAME) == 'RI') THEN
            ri(1:length) = aryrval(1:length,n2rary)
          ENDIF

        ENDIF
!
      ENDDO
!
!-----------------------------------------------------------------------
!***  LOGICAL HISTORY VARIABLES
!-----------------------------------------------------------------------
!
      N2LSCALAR = wrt_int_state%KOUNT_LOG(NBDL)                               !<-- Counter for full string of logical data
!
      IF(N2LSCALAR > 0) then
        ALLOCATE(VARLNAME(N2LSCALAR),VARLVAL(N2LSCALAR))
      else
        ALLOCATE(VARLNAME(1),VARLVAL(1))
      ENDIF
      VARLVAL = .false.
      N2LSCALAR = 0
!
      DO N=1,wrt_int_state%KOUNT_LOG(NBDL)                                    !<-- Loop through all logical data
!
        NPOSN_1 = (N-1)*NAME_MAXSTR+1
        NPOSN_2 = N*NAME_MAXSTR
        NAME    = wrt_int_state%NAMES_LOG_STRING(NBDL)(NPOSN_1:NPOSN_2)       !<-- The variable's name
!
        N2LSCALAR = N2LSCALAR+1
        WORK_LOGICAL = wrt_int_state%ALL_DATA_LOG(N2LSCALAR,NBDL)                 !<-- Extract the individual data from the data string
        VARLNAME(N2LSCALAR) = NAME
        if(WORK_LOGICAL == ESMF_TRUE) VARLVAL(N2LSCALAR) = .true.
        IF(TRIM(NAME) == 'GLOBAL') GLOBAL = VARLVAL(N2LSCALAR)
        IF(TRIM(NAME) == 'GEN_COORD_HYBRID') THEN
          GEN_COORD_HYBRID  = VARLVAL(N2LSCALAR)
          NGEN_COORD_HYBRID = N2LSCALAR
        endif
        IF(TRIM(NAME) == 'HYBRID') THEN
          HYBRID  = VARLVAL(N2LSCALAR)
          NHYBRID = N2LSCALAR
        ENDIF
!
      ENDDO
!      write(0,*)'in nemsio, GEN_COORD_HYBRID=',GEN_COORD_HYBRID, &
!          'HYgrid=',hybrid
!
!-----------------------------------------------------------------------
!***  NOW OPEN NEMSIO FILE
!-----------------------------------------------------------------------
!
!      write(0,*)' OPEN_NEMSIO_FILE wrt_int_state%IO_NEMSIOFILE=',        &
!          trim(wrt_int_state%FILENAME_BASE(NBDL)),'idate=',idate,        &
!          'NF_HOURS=',NF_HOURS,'NF_MINUTES=',NF_MINUTES,'NF_SECONDS=',   &
!          NF_SECONDS
!
      IF(wrt_int_state%IO_FILE(NBDL) == 'DEFERRED') THEN
        LFNHR = .true.    ! no output
!       lfnhr=3600*abs(fhour-nint(fhour)).le.1.or.phour.eq.0
        LFNHR = (NF_MINUTES == 0 .and. NINT(NF_SECONDS) == 0)
        IF(LFNHR) THEN
          NDIG = MAX(LOG10(NF_HOURS+0.5)+1.,2.)
          WRITE(CFORM,'("(I",I1,".",I1,")")') NDIG,NDIG
          WRITE(CFHOUR,CFORM) NF_HOURS
        ELSE
          NDIG = MAX(LOG10(NF_HOURS+0.5)+1.,2.)
          WRITE(CFORM,'("(I",I1,".",I1,",A1,I2.2,A1,I2.2)")') NDIG,NDIG
          WRITE(CFHOUR,CFORM) NF_HOURS,':',NF_MINUTES,':',NINT(NF_SECONDS)
        ENDIF
        CFHOUR = trim(CFHOUR) // trim(ensmem_name)
!
        N = LEN_TRIM(wrt_int_state%FILENAME_BASE(NBDL))
        FILENAME = wrt_int_state%FILENAME_BASE(NBDL)(1:n)//trim(CFHOUR)
      ELSE
        FILENAME = wrt_int_state%FILENAME_BASE(NBDL)//'_nemsio'
      ENDIF
!
!----------------------------------------------------
!***  Prepare variables needed by the nemsip header:
!----------------------------------------------------
!
!dimension
      IF(wrt_int_state%core == 'NMMB' .and. GLOBAL) THEN
        NFRAME = 1           !for global im/jm for data field
      ELSE
        NFRAME = 0           !for regional
      ENDIF
      IM   = wrt_int_state%im(1)
      JM   = wrt_int_state%jm(1)
      DIM1 = wrt_int_state%im(1)-2*NFRAME
      DIM2 = wrt_int_state%jm(1)-2*NFRAME
!
      LM   = wrt_int_state%LM(1)
!
!for nmmb whole domain
      FIELDSIZE = IM*JM
      NREC      = wrt_int_state%kount_I2D(NBDL)+wrt_int_state%kount_R2D(NBDL)

!      write(0,*)'before vcoord,nbdl=',nbdl,'Kount_i2d=',wrt_int_state%kount_I2D(NBDL) &
!       ,'Kount_r2d=',wrt_int_state%kount_R2D(NBDL)
!
!vcoord
      ALLOCATE(VCOORD(LM+1,3,2))
      VCOORD = 0.
      if(wrt_int_state%core == 'gfs') then
        idvm    = Nthermodyn_id*10 + Nsfcpress_id    ! user specified
!for output:
!       idvm    = 22                                 ! 1:  ln(ps) 2:ps   ! hmhj
                                                     ! 1: Tv, 2: T, 3:Th
        IF(gen_coord_hybrid) then
          idvc    = 3
          idsl    = 2    ! idsl=2 for middle of layer                   ! hmhj
          nvcoord = 3
          if(NAK5 > 0 .and. NBK5 > 0 .and. NCK5 > 0 ) then
            VCOORD(1:LM+1,1,1) = ARYRVAL(1:LM+1,NAK5)*1000.
            VCOORD(1:LM+1,2,1) = ARYRVAL(1:LM+1,NBK5)
            VCOORD(1:LM+1,3,1) = ARYRVAL(1:LM+1,NCK5)*1000.
          endif
        ELSEIF(hybrid) then
          idvc    = 2                        ! for hybrid vertical coord.
          idsl    = 1
          nvcoord = 2
          if(NAK5>0 .and. NBK5>0 ) then
            do i=1,LM+1
              vcoord(i,1,1) = ARYRVAL(LM+1+1-i,NAK5)*1000.
              vcoord(i,2,1) = ARYRVAL(LM+1+1-i,NBK5)
            enddo
            VCOORD(1:LM+1,3,1) = 0
          endif
        ELSE if(NGEN_COORD_HYBRID /= 0 .and. NHYBRID /= 0) THEN
          idvc    = 1    ! for sigma vertical coord. (default)
          idsl    = 1
          nvcoord = 1
          if(NSI > 0 ) then
            vcoord(:,1,1) = ARYRVAL(:,NSI)
          endif
        ENDIF 
!
      endif 
!
!-----------------------------------------------------------------------
!***  Cut the output I2D array.
!-----------------------------------------------------------------------
!
      ALLOCATE(RECNAME(NREC),RECLEVTYP(NREC),RECLEV(NREC))
      NREC = 0
      INI1 = 0
      INI2 = 0
!
      DO NFIELD=1,wrt_int_state%KOUNT_I2D(NBDL)
!
        NREC    = NREC + 1
        NPOSN_1 = (NFIELD-1)*NAME_MAXSTR+1
        NPOSN_2 = NFIELD*NAME_MAXSTR
        NAME    = wrt_int_state%NAMES_I2D_STRING(NBDL)(NPOSN_1:NPOSN_2)          !<-- The name of this 2D integer history quantity
        INDX_2D = index(NAME,"_2D")
!
        IF (INDX_2D > 0) THEN
          INDX_2DA     = INDEX(NAME(1:INDX_2D-1),"_",back=.true.)
          RECLEV(NREC) = 0
          DO I=1, INDX_2D-INDX_2DA-1
            RECLEV(NREC) = (ICHAR(NAME(INDX_2D-i:INDX_2D-i))-48)*10**(I-1)+RECLEV(NREC)
          ENDDO
          INDX_2D2     = INDEX(NAME,"_")
          if(INDX_2D2 > 0) RECNAME(NREC) = NAME(1:INDX_2D2-1)
          CALL LOWERCASE(RECNAME(NREC))
          if(INDX_2D-4 > INDX_2D2)  then
             RECLEVTYP(NREC) = NAME(INDX_2D2+1:INDX_2D-4)
          else
             RECLEVTYP(NREC) = 'mid layer'
          endif
          IF (RECLEV(NREC) == LM+1) RECLEVTYP(NREC-LM:NREC) = 'layer'

          INI1 = INI1 + 1
        ELSE
          RECLEV(NREC) = 1
          INDX_2D2 = INDEX(NAME,"_")
          INDX_2D3 = INDEX(NAME,"_",BACK=.true.)
          if(INDX_2D2 > 0) then
            RECNAME(NREC) = trim(NAME(1:INDX_2D2-1))
          else
            RECNAME(NREC) = TRIM(NAME)
          endif
          if(INDX_2D3 > 0) then
            RECLEVTYP(NREC) = NAME(INDX_2D3+1:)
          else
            RECLEVTYP(NREC) = 'sfc'
          endif
          CALL LOWERCASE(RECNAME(NREC))

          RECLEV(NREC) = 1
          INI2         = INI2 + 1
        ENDIF
!
        IF (RECNAME(NREC) == 'ISLTYP') RECNAME(NREC) = 'sltyp'
        IF (RECNAME(NREC) == 'IVGTYP') RECNAME(NREC) = 'vgtyp'
        IF (RECNAME(NREC) == 'NCFRCV') RECNAME(NREC) = 'cfrcv'
        IF (RECNAME(NREC) == 'NCFRST') RECNAME(NREC) = 'cfrst'
        CALL LOWERCASE(RECNAME(NREC))
      ENDDO
!
!-----------------------------------------------------------------------
!*** Cut the output R2D array.
!-----------------------------------------------------------------------
!
      LSOIL = 0
      IND4  = 0
      IND3  = 0
      IND2  = 0
      IND1  = 0
!
      DO NFIELD=1,wrt_int_state%KOUNT_R2D(NBDL)
!
        NREC    = NREC + 1
        NPOSN_1 = (NFIELD-1)*NAME_MAXSTR+1
        NPOSN_2 = NFIELD*NAME_MAXSTR
        NAME    = wrt_int_state%NAMES_R2D_STRING(NBDL)(NPOSN_1:NPOSN_2)  !<-- The name of this 2D integer history quantity
        INDX_2D = INDEX(NAME,"_2D")
!
        IF (INDX_2D > 0) THEN
          INDX_2DA     = INDEX(NAME(1:INDX_2D-1),"_",back=.true.)
          RECLEV(NREC) = 0
          DO I=1, INDX_2D-INDX_2DA-1
            RECLEV(NREC) = (ICHAR(NAME(INDX_2D-i:INDX_2D-i))-48)*10**(I-1)+RECLEV(NREC)
          ENDDO
          INDX_2D2     = INDEX(NAME,"_")
          if(INDX_2D2 > 0) RECNAME(NREC) = NAME(1:INDX_2D2-1)
          CALL LOWERCASE(RECNAME(NREC))
          if(INDX_2D-4 > INDX_2D2)  then
             RECLEVTYP(NREC) = NAME(INDX_2D2+1:INDX_2D-4)
          else
             RECLEVTYP(NREC) = 'mid layer'
          endif
          IF (RECLEV(NREC) == LM+1) RECLEVTYP(NREC-LM:NREC) = 'layer'
!
          IF (RECNAME(NREC) == 'smc'.or. RECNAME(NREC) == 'soilw') LSOIL = LSOIL + 1
          IF (RECNAME(NREC) == 'w')  RECNAME(NREC) = 'vvel'
          IF (RECNAME(NREC) == 'cw') RECNAME(NREC) = 'clwmr'
          IF (RECNAME(NREC) == 'u')  RECNAME(NREC) = 'ugrd'
          IF (RECNAME(NREC) == 'v')  RECNAME(NREC) = 'vgrd'
          IF (RECNAME(NREC) == 't')  RECNAME(NREC) = 'tmp'
          IF (RECNAME(NREC) == 'q')  RECNAME(NREC) = 'spfh'
          IF (RECNAME(NREC) == 'pint') THEN
          RECNAME(NREC) = 'pres'
          IND1 = IND1 + 1
          ELSE IF (RECNAME(NREC) == 'smc' .OR. RECNAME(NREC) == 'sh2o' .or. RECNAME(NREC) == 'stc' &
              .or. RECNAME(NREC) == 'slc' ) THEN 
            RECLEVTYP(NREC) = 'soil layer' 
            IND2 = IND2 + 1
          ELSE
            IND3 = IND3 + 1
          ENDIF
        ELSE
          RECLEV(NREC) = 1
          INDX_2D2 = INDEX(NAME,"_")
          INDX_2D3 = INDEX(NAME,"_",BACK=.true.)
          if(INDX_2D2 > 0) then
            RECNAME(NREC) = trim(NAME(1:INDX_2D2-1))
          else
            RECNAME(NREC) = TRIM(NAME)
          endif
          if(INDX_2D3 > 0) then
            RECLEVTYP(NREC) = NAME(INDX_2D3+1:)
          else
            RECLEVTYP(NREC) = 'sfc'
          endif
          CALL LOWERCASE(RECNAME(NREC))
!
          IF (RECNAME(NREC) == 'pd') THEN
            RECNAME(NREC)   = 'dpres'
            RECLEVTYP(NREC) = 'hybrid sig lev'
          ENDIF
!
          IF (RECNAME(NREC) == 'pressfc') RECNAME(NREC) = 'pres'
          IF (RECNAME(NREC) == 'sst')     RECNAME(NREC) = 'tsea'
          IF (RECNAME(NREC) == 'fis')     RECNAME(NREC) = 'hgt'
          IF (RECNAME(NREC) == 'ustar')   RECNAME(NREC) = 'uustar'
          IF (RECNAME(NREC) == 'z0')      RECNAME(NREC) = 'zorl'
          IND4 = IND4 + 1
        ENDIF
!
      ENDDO
      if(lsoil /= 0) NSOIL = LSOIL
!
!set meta data record
      NMETA = 12
      IF (NAK5 == 0 .and. NBK5 == 0 .and. NCK5 == 0) NMETA = 8
!
!-----------------------------------------------------------------------
!                      SET UP NEMSIO WRITE
!-----------------------------------------------------------------------
!
      CALL NEMSIO_INIT(IRET=IRET)
!      write(0,*)'after nemsio_init, iret=',iret,'dim1=',dim1,'dim2=',dim2, &
!       'nsoil=',nsoil,'ntrac=',ntrac,'ncldt=',ncld,'nrec=',nrec,'cpi=',cpi, &
!       'ri=',ri,'nmetavarl=',n2lscalar,'nmetaaryr=',n2rary,'idsl=',idsl,    &
!       'idvc=',idvc,'idvm=',idvm
!      write(0,*)'before nemsio_open,variname=',variname,varival
!
!-----------------------------------------------------------------------
!***  OPEN NEMSIO FILE
!-----------------------------------------------------------------------
!
      CALL NEMSIO_OPEN(NEMSIOFILE,trim(FILENAME),'write',iret,           &
        modelname="GFS", gdatatype=wrt_int_state%io_form(NBDL),          &
        idate=idate,nfhour=NF_HOURS_IAU,                                 &
        nfminute=NF_MINUTES,nfsecondn=nint(NF_SECONDS*100),              &
        nfsecondd=100,dimx=DIM1,dimy=DIM2,dimz=LM,nframe=NFRAME,         &
        nmeta=NMETA,jcap=JCAP,idsl=IDSL,idvm=IDVM,idvc=IDVC,idrt=IDRT,   &
        nsoil=NSOIL,ntrac=ntrac,nrec=nrec, ncldt=ncld,                   &
        vcoord=vcoord,cpi=cpi,ri=ri,                                     &
        extrameta=.true.,nmetavari=N2ISCALAR,                            &
        nmetavarr=N2RSCALAR,nmetavarl=N2LSCALAR,nmetaaryi=N2IARY,        &
        nmetaaryr=N2RARY,variname=VARINAME,varival=VARIVAL,              &
        varrname=VARRNAME,varrval=VARRVAL,varlname=VARLNAME,             &
        varlval=VARLVAL,aryiname=ARYINAME,aryilen=ARYILEN,               &
        aryival=ARYIVAL,aryrname=ARYRNAME,aryrlen=ARYRLEN,               &
        aryrval=ARYRVAL,recname=RECNAME,reclevtyp=RECLEVTYP,reclev=RECLEV)

       if(iret /= 0) print *,'nemsio_open, file=',trim(filename),' iret=',iret
!
!-----------------------------------------------------------------------
!***  CLEAN UP
!-----------------------------------------------------------------------
!
      DEALLOCATE(VCOORD,CPI,RI)
      if(associated(VARINAME))  DEALLOCATE(VARINAME,VARIVAL)
      if(associated(ARYINAME))  DEALLOCATE(ARYINAME,ARYILEN,ARYIVAL)
      if(associated(VARRNAME))  DEALLOCATE(VARRNAME,VARRVAL)
      if(associated(ARYRNAME))  DEALLOCATE(ARYRNAME,ARYRLEN,ARYRVAL)
      if(associated(varlname))  DEALLOCATE(VARLNAME,VARLVAL)
!
      NULLIFY(VARINAME)
      NULLIFY(VARIVAL)
      NULLIFY(VARRNAME)
      NULLIFY(VARRVAL)
      NULLIFY(VARLNAME)
      NULLIFY(VARLVAL)
      NULLIFY(ARYINAME)
      NULLIFY(ARYILEN)
      NULLIFY(ARYIVAL)
      NULLIFY(ARYRNAME)
      NULLIFY(ARYRLEN)
      NULLIFY(ARYRVAL)
      write(0,*)'end of write_nemsio_open'
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE WRITE_NEMSIO_OPEN
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
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
!***  CONVERT MONTH
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
      END MODULE MODULE_WRITE_ROUTINES_GFS
!
!-----------------------------------------------------------------------
