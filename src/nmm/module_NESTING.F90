!-----------------------------------------------------------------------
!
      MODULE MODULE_NESTING
!
!-----------------------------------------------------------------------
!
!***  This module contains routines that perform various interactions
!***  between parent domains and their children.
!
!-----------------------------------------------------------------------
!
! PROGRAM HISTORY LOG:
!
!   2008-02-07  Black - PARENT_TO_CHILD_FILL
!   2008-03-05  Black - PARENT_CHILD_SPLIT
!   2008-03-25  Black - PARENT_TO_CHILD_INIT_NMM
!   2008-04-22  Black - Replace PARENT_CHILD_SPLIT with _COMMS
!   2008-06-18  Black - PARENT_TO_CHILD_COMPUTE
!   2008-06-18  Black - PREPARE_PARENT_TO_CHILD_INTERP
!   2008-08-14  Black - Added BOUNDARY_DATA_STATE_TO_STATE
!   2009-03-12  Black - Added Z0BASE and STDH now needed for NPS.
!   2009-10-12  Black - Fix for generalized of parent-child space ratios.
!   2010-03-31  Black - Add parent computation of child boundary topo.
!   2011-05-17  Yang  - Modified for using the ESMF 5.2.0r_beta_snapshot_07.
!   2011-07-16  Black - Moving nest capability.
!   2012-07-20  Black - Generational use of MPI tasks.
!
!-----------------------------------------------------------------------
!
! USAGE: 
!
!-----------------------------------------------------------------------
!
      USE MPI
      USE ESMF
      USE netcdf
!
      USE module_KINDS
!
      USE module_DERIVED_TYPES,ONLY: BNDS_2D                            &
                                    ,CHILD_UPDATE_LINK                  &
                                    ,COMMS_FAMILY                       &
                                    ,DOMAIN_DATA                        &
                                    ,INTEGER_DATA                       &
                                    ,INTERIOR_DATA_FROM_PARENT          &
                                    ,MIXED_DATA_TASKS                   &
                                    ,REAL_DATA_2D
!
      USE module_VARS,ONLY: VAR
!
      USE module_LS_NOAHLSM,ONLY: NUM_SOIL_LAYERS
!
      USE module_CONSTANTS,ONLY: P608,R_D
!
      USE module_CONTROL,ONLY: NUM_DOMAINS_MAX,TIMEF
!
      USE module_EXCHANGE,ONLY: HALO_EXCH
!
      USE module_ERR_MSG,ONLY: ERR_MSG,MESSAGE_CHECK
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      PRIVATE
!
      PUBLIC :: BOUNDARY_DATA_STATE_TO_STATE                            &
               ,CHECK                                                   &
               ,CHECK_REAL                                              &
               ,CHILD_2WAY_BOOKKEEPING                                  &
               ,CHILD_RANKS                                             &
               ,GENERATE_2WAY_DATA                                      &
               ,HYPERBOLA                                               &
               ,INTERNAL_DATA_TO_DOMAIN                                 &
               ,INTERIOR_DATA_STATE_TO_STATE                            &
               ,LAG_STEPS                                               &
               ,LATLON_TO_IJ                                            &
               ,MOVING_NEST_BOOKKEEPING                                 &
               ,MOVING_NEST_RECV_DATA                                   &
               ,PARENT_2WAY_BOOKKEEPING                                 &
               ,PARENT_BOOKKEEPING_MOVING                               &
               ,PARENT_CHILD_COMMS                                      &
               ,PARENT_READS_MOVING_CHILD_TOPO                          &
               ,PARENT_TO_CHILD_INIT_NMM                                &
               ,PARENT_UPDATES_HALOS                                    &
               ,PARENT_UPDATES_MOVING                                   &
               ,REAL_IJ_TO_LATLON                                       &
               ,SET_NEST_GRIDS                                          &
               ,STENCIL_H_EVEN,STENCIL_SFC_H_EVEN                       &
               ,STENCIL_V_EVEN,STENCIL_SFC_V_EVEN                       &
               ,STENCIL_H_ODD,STENCIL_SFC_H_ODD                         &
               ,STENCIL_V_ODD,STENCIL_SFC_V_ODD                         &
               ,SUFFIX_MOVE                                             &
               ,SUFFIX_NESTBC                                           &
               ,SUFFIX_TWOWAY
!
!-----------------------------------------------------------------------
!
      INTEGER(kind=KINT),PARAMETER :: NEAREST=0                         &  !<-- Flag for nearest neighbor interpolation (parent to child)
                                     ,BILINEAR=1                           !<-- Flag for bilinear interpolation (parent to child)
!
      INTEGER(kind=KINT),SAVE :: LM,N8=8
!
      INTEGER(kind=KINT),SAVE :: LAG_STEPS=4                               !<-- Nest moves this many parent timesteps after deciding
!
      INTEGER(kind=KINT),SAVE :: STENCIL_H_EVEN=3                       &
                                ,STENCIL_V_EVEN=2                       &
                                ,STENCIL_SFC_H_EVEN=3                   &
                                ,STENCIL_SFC_V_EVEN=3                   &
                                ,STENCIL_H_ODD=3                        &
                                ,STENCIL_V_ODD=3                        &
                                ,STENCIL_SFC_H_ODD=3                    &
                                ,STENCIL_SFC_V_ODD=2
!
      REAL(kind=KFPT),SAVE :: CHILD_PARENT_SPACE_RATIO                  &
                             ,EPS=1.E-4 
!
      CHARACTER(len=5) :: SUFFIX_MOVE='-move'
      CHARACTER(len=5) :: SUFFIX_TWOWAY='-2way'
      CHARACTER(len=7) :: SUFFIX_NESTBC='-nestbc'
!
!-----------------------------------------------------------------------
!
      REAL(kind=KDBL) :: btim,btim0
!
      TYPE(CHILD_UPDATE_LINK),POINTER,SAVE :: TAIL
!
      TYPE(DOMAIN_DATA),DIMENSION(:),POINTER,SAVE :: CHILD_RANKS            !<-- Lists of child tasks' local ranks in p-c intracomms
!
      integer(kind=kint) :: iprt=01,jprt=61
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE PARENT_CHILD_COMMS(MYPE                                &
                                   ,NUM_DOMAINS_TOTAL                   &
                                   ,NUM_TASKS_TOTAL                     &
                                   ,COMM_WORLD                          &
                                   ,RANK_TO_DOMAIN_ID                   &
                                   ,CF                                  &
                                   ,TASK_MODE                           &
                                   ,QUILTING                            &
                                   ,DOMAIN_GEN                          &
                                   ,FULL_GEN                            &
                                   ,MY_DOMAIN_ID_N                      &
                                   ,ID_DOMAINS                          &
                                   ,ID_PARENTS                          &
                                   ,NUM_CHILDREN                        &
                                   ,ID_CHILDREN                         &
                                   ,COMMS_DOMAIN                        &
                                   ,FTASKS_DOMAIN                       &
                                   ,NTASKS_DOMAIN                       &
                                   ,PETLIST_DOM                         &
                                   ,NUM_GENS                            &
                                         )
!
!-----------------------------------------------------------------------
!***  Create MPI intracommunicators between the tasks of a parent domain
!***  and those of all its 1st generation nests (children).
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      INTEGER(kind=KINT),INTENT(IN) :: COMM_WORLD                       &  !<-- MPI intracommunicator for ALL tasks
                                      ,FULL_GEN                         &  !<-- For 2-way nesting, the generation using all tasks
                                      ,MYPE                             &  !<-- My task ID (global)
                                      ,NUM_DOMAINS_TOTAL                &  !<-- Total number of domains
                                      ,NUM_TASKS_TOTAL                     !<-- Total number of tasks in the run
!
      INTEGER(kind=KINT),DIMENSION(*),INTENT(IN) :: DOMAIN_GEN          &  !<-- For 2-way nesting, each domain's generation
                                                   ,RANK_TO_DOMAIN_ID      !<-- Domain ID for each configure file
!
      CHARACTER(len=12),INTENT(IN) :: TASK_MODE                            !<-- Unique or generational task assignment
!
      LOGICAL(kind=KLOG),INTENT(IN) :: QUILTING                            !<-- Was quilting specified in the configure files?
!
      TYPE(ESMF_Config),DIMENSION(*),INTENT(INOUT) :: CF                   !<-- The config objects (one per domain)
!
      INTEGER(kind=KINT),INTENT(OUT) :: NUM_GENS                           !<-- The # of generations of domains
!
      INTEGER(kind=KINT),DIMENSION(:,:),POINTER,INTENT(OUT) :: ID_CHILDREN &  !<-- Domain IDs of all domains' children
                                                              ,PETLIST_DOM    !<-- List of task IDs on each domain
!
      INTEGER(kind=KINT),DIMENSION(:),POINTER,INTENT(OUT) :: ID_DOMAINS      &  !<-- Array of the domain IDs
                                                            ,ID_PARENTS      &  !<-- Array of the domains' parent IDs
                                                            ,FTASKS_DOMAIN   &  !<-- # of forecast tasks on each domain
                                                            ,MY_DOMAIN_ID_N  &  !<-- IDs of the domains on which current task resides
                                                            ,NTASKS_DOMAIN   &  !<-- # of tasks on each domain excluding descendents
                                                            ,NUM_CHILDREN       !<-- # of children on each domain
!
      TYPE(COMMS_FAMILY),DIMENSION(:),POINTER,INTENT(OUT) :: COMMS_DOMAIN  !<-- Intracommunicators between parent and child domains
                                                                           !    and between each domains' forecast tasks.
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: I,IERR,ISTAT                                &
                           ,N,N1,N2,N3,NN,RC
!
      INTEGER(kind=KINT) :: COMM_INTRA                                  &
                           ,GROUP_UNION                                 &
                           ,GROUP_WORLD                                 &
                           ,ID_CHILD                                    &
                           ,ID_DOM                                      &
                           ,ID_FULL                                     &
                           ,ID_PARENT                                   &
                           ,INPES                                       &
                           ,JNPES                                       &
                           ,KOUNT                                       &
                           ,KOUNT_DOMS                                  &
                           ,KOUNT_TASKS                                 &
                           ,LAST_FCST_TASK_X                            &
                           ,LAST_WRITE_TASK_X                           &
                           ,LEAD_REMOTE                                 &
                           ,N_CHILDREN                                  &
                           ,N_GEN                                       &
                           ,NDOMS_FULL                                  &
                           ,NMAX                                        &
                           ,NSAVE                                       &
                           ,NTASKS_CONTRIB                              &
                           ,NTASKS_PARENT                               &
                           ,NTASKS_X                                    &
                           ,NUM_FCST_TASKS                              &
                           ,NUM_TASKS_FULL                              &
                           ,NUM_WRITE_TASKS                             &
                           ,RC_COMMS                                    &
                           ,TASK_X                                      &
                           ,WRITE_GROUPS                                &
                           ,WRITE_TASKS_PER_GROUP
!
      INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE :: DOMS_PER_GEN       &  !<-- Domain count per generation
                                                    ,DOMS_FULL          &  !<-- IDs of domains in the full generation
                                                    ,GLOBAL_UNION       &  !<-- Union of parent and child tasks in intracomms
                                                    ,GROUP              &  !<-- MPI group for each domain
                                                    ,KOUNT_FULL         &  !<-- 1st task on each domain of the full generation
                                                    ,LAST_FCST_TASK     &  !<-- ID of last forecast task on each domain
                                                    ,LAST_WRITE_TASK    &  !<-- ID of last write task on each domain
                                                    ,LAST_TASK          &  !<-- ID of last task on each domain
                                                    ,LEAD_FCST_TASK     &  !<-- ID of first task on each domain
                                                    ,LEAD_WRITE_TASK    &  !<-- ID of first write on each domain
                                                    ,LEAD_TASK          &  !<-- ID of first task on each domain
                                                    ,WTASKS_DOMAIN         !<-- # of write/quilt tasks on each domain
!
      REAL(kind=KFPT) :: RECIP_NUM_TASKS_FULL
!
      REAL(kind=KFPT),DIMENSION(:),ALLOCATABLE :: FRAC_FULL                !<-- Fraction of tasks on each domain in full generation
!
      CHARACTER(2)      :: NUM_DOMAIN
      CHARACTER(6),SAVE :: FMT='(I2.2)'
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC_COMMS=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!
      ALLOCATE(ID_DOMAINS   (1:NUM_DOMAINS_TOTAL))       
      ALLOCATE(ID_PARENTS   (1:NUM_DOMAINS_TOTAL))       
      ALLOCATE(LEAD_TASK    (1:NUM_DOMAINS_TOTAL))       
      ALLOCATE(LAST_TASK    (1:NUM_DOMAINS_TOTAL))       
      ALLOCATE(FTASKS_DOMAIN(1:NUM_DOMAINS_TOTAL))       
      ALLOCATE(WTASKS_DOMAIN(1:NUM_DOMAINS_TOTAL))       
      ALLOCATE(NTASKS_DOMAIN(1:NUM_DOMAINS_TOTAL))       
      ALLOCATE(NUM_CHILDREN (1:NUM_DOMAINS_TOTAL))                          
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  Incoming tasks extract relevant information from all config files.
!
!***  This loop is general thus the domain IDs do not need to correspond
!***  to the number in the configure file name.  The user may assign
!***  IDs monotonically to the domains starting with 1 and in any order
!***  desired except that the uppermost parent must have an ID of 1.
!***  However the rank/element of each domain in the CF array is equal
!***  to the given domain's ID.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      read_configs: DO N=1,NUM_DOMAINS_TOTAL                               !<-- Loop through all configure objects
!
!-----------------------------------------------------------------------
!***  Save the domain IDs.
!***  These are simply integers each domain will use to keep track
!***  of itself with respect to others.
!-----------------------------------------------------------------------
!
        ID_DOMAINS(N)=RANK_TO_DOMAIN_ID(N)
        ID_DOM=ID_DOMAINS(N)
!
!-----------------------------------------------------------------------
!***  Who is the parent of each domain?
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract ID of Parent of this Domain"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=CF(ID_DOM)                  &  !<-- The config object
                                    ,value =ID_PARENTS(ID_DOM)          &  !<-- The ID of the parent of this domain
                                    ,label ='my_parent_id:'             &  !<-- Take values from this config label 
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_COMMS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!-----------------------------------------------------------------------
!***  How many children does each domain have?
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract # of Children of this Domain"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=CF(ID_DOM)                  &  !<-- The config object
                                    ,value =NUM_CHILDREN(ID_DOM)        &  !<-- # of children on this domain
                                    ,label ='n_children:'               &  !<-- Take value from this config label 
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_COMMS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  How many Forecast/Write tasks will be active on each domain?
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Parent_Child_Comms: Extract INPES From Config File"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=CF(ID_DOM)                  &  !<-- The config object
                                    ,value =INPES                       &  !<-- The domain's fcst tasks in I
                                    ,label ='inpes:'                    &  !<-- Give this label's value to the previous variable
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_COMMS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Parent_Child_Comms: Extract JNPES From Config File"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=CF(ID_DOM)                  &  !<-- The config object
                                    ,value =JNPES                       &  !<-- The domain's fcst tasks in J
                                    ,label ='jnpes:'                    &  !<-- Give this label's value to the previous variable
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_COMMS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Parent_Child_Comms: Extract Write_Groups From Config File"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=CF(ID_DOM)                  &  !<-- The config object
                                    ,value =WRITE_GROUPS                &  !<-- The number of Write groups on this domain
                                    ,label ='write_groups:'             &  !<-- Give this label's value to the previous variable
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_COMMS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Parent_Child_Comms: Extract Write_Task_Per_Group From Config File"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=CF(ID_DOM)                  &  !<-- The config object
                                    ,value =WRITE_TASKS_PER_GROUP       &  !<-- The number of tasks per Write group 
                                    ,label ='write_tasks_per_group:'    &  !<-- Give this label's value to the previous variable
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_COMMS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        IF(.NOT.QUILTING)THEN
          WRITE_GROUPS=0
          WRITE_TASKS_PER_GROUP=0
        ENDIF
!
!-----------------------------------------------------------------------
!
        FTASKS_DOMAIN(ID_DOM)=INPES*JNPES                                  !<-- # of compute/forecast tasks on domain ID_DOM
!
        WTASKS_DOMAIN(ID_DOM)=WRITE_GROUPS*WRITE_TASKS_PER_GROUP           !<-- # of write/quilt tasks on domain ID_DOM
!
        NTASKS_DOMAIN(ID_DOM)=FTASKS_DOMAIN(ID_DOM)                     &  !<-- Total # of tasks on each domain
                             +WTASKS_DOMAIN(ID_DOM)              
!
!-----------------------------------------------------------------------
!
      ENDDO read_configs
!
      ALLOCATE(PETLIST_DOM(1:NUM_TASKS_TOTAL,1:NUM_DOMAINS_TOTAL))
!
      DO N2=1,NUM_DOMAINS_TOTAL
      DO N1=1,NUM_TASKS_TOTAL
        PETLIST_DOM(N1,N2)=-999
      ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
!***  Assign tasks to all domains.  
!***  For 1-way nesting each task is uniquely assigned to a single
!***  domain.  For 2-way nesting each forecast task can be assigned
!***  to more than one domain but cannot lie on more than one domain
!***  in each generation.  The write/quilt tasks in 2-way nesting
!***  must be assigned uniquely to a single domain and they cannot 
!***  also be forecast tasks or else the asynchronous writing of
!***  the history/restart files would not always be independent of
!***  the forecast integration.
!-----------------------------------------------------------------------
!
      ALLOCATE(MY_DOMAIN_ID_N(1:NUM_DOMAINS_TOTAL),stat=ISTAT)
!
      DO N=1,NUM_DOMAINS_TOTAL
        MY_DOMAIN_ID_N(N)=0
      ENDDO
!
      IF(ISTAT/=0)THEN
        WRITE(0,*)' Failed to allocate MY_DOMAIN_ID_N  rc=',ISTAT
      ENDIF
!
!-----------------------------------------------------------------------
!
      task_assign: IF(TASK_MODE=='unique')THEN
!
!-----------------------------------------------------------------------
!
        DO N=1,NUM_DOMAINS_TOTAL
!
          ID_DOM=RANK_TO_DOMAIN_ID(N)
!
!-----------------------------------------------------------------------
!***  Determine the global IDs of the lead and last tasks on each
!***  domain which in turn lets us fill the PETLIST for each domain.
!***  These include the I/O tasks.
!-----------------------------------------------------------------------
!
          IF(N==1)THEN
            LEAD_TASK(N)=0                                                   !<-- Task 0 is first in line
          ELSE
            LEAD_TASK(ID_DOM)=LAST_TASK(ID_DOMAINS(N-1))+1                   !<-- Lead task on domain follows last task on previous domain
          ENDIF
!
          LAST_TASK(ID_DOM)=LEAD_TASK(ID_DOM)+NTASKS_DOMAIN(ID_DOM)-1        !<-- The last task on each domain
!
          IF(MYPE>=LEAD_TASK(ID_DOM).AND.MYPE<=LAST_TASK(ID_DOM))THEN        !<-- Associate tasks with each domain
            MY_DOMAIN_ID_N(1)=ID_DOM                                         !<-- Tell this task the ID of the single domain it is on
          ENDIF
!
          KOUNT_TASKS=0
          DO N2=LEAD_TASK(ID_DOM),LAST_TASK(ID_DOM)
            KOUNT_TASKS=KOUNT_TASKS+1
            PETLIST_DOM(KOUNT_TASKS,ID_DOM)=N2
          ENDDO
!
        ENDDO
!
        NUM_GENS=1                                                           !<-- This is a dummy value; only relevant for 2-way
!
!-----------------------------------------------------------------------
!
      ELSEIF(TASK_MODE=='generational')THEN
!
!-----------------------------------------------------------------------
!***  First determine how many domains are in each generation.
!-----------------------------------------------------------------------
!
        NUM_GENS=0
        NUM_WRITE_TASKS=0
!
        ALLOCATE(DOMS_PER_GEN(1:NUM_DOMAINS_TOTAL))
!
        DO N=1,NUM_DOMAINS_TOTAL
          DOMS_PER_GEN(N)=0
        ENDDO
!
        DO N=1,NUM_DOMAINS_TOTAL
          ID_DOM=RANK_TO_DOMAIN_ID(N)
          N_GEN=DOMAIN_GEN(ID_DOM)                                         !<-- The generation that domain ID_DOM is in
          IF(N_GEN>NUM_GENS)NUM_GENS=N_GEN                                 !<-- Determining the # of generations
          DOMS_PER_GEN(N_GEN)=DOMS_PER_GEN(N_GEN)+1                        !<-- Determining the # of domains per generation
          NUM_WRITE_TASKS=NUM_WRITE_TASKS+WTASKS_DOMAIN(ID_DOM)            !<-- Sum write tasks for all domains
        ENDDO
!
!-----------------------------------------------------------------------
!***  Assign all the run's forecast tasks across the first generation
!***  that uses all of them.  This is the first 'full' generation.
!-----------------------------------------------------------------------
!
        ALLOCATE(LEAD_FCST_TASK(1:NUM_DOMAINS_TOTAL))
        ALLOCATE(LAST_FCST_TASK(1:NUM_DOMAINS_TOTAL))
        ALLOCATE(LEAD_WRITE_TASK(1:NUM_DOMAINS_TOTAL))
        ALLOCATE(LAST_WRITE_TASK(1:NUM_DOMAINS_TOTAL))
!
        KOUNT_DOMS=0
!
        NUM_FCST_TASKS=NUM_TASKS_TOTAL-NUM_WRITE_TASKS                     !<-- Total # of forecast tasks available
!
        DO N=1,NUM_DOMAINS_TOTAL
!
          ID_DOM=RANK_TO_DOMAIN_ID(N)
!
          IF(DOMAIN_GEN(ID_DOM)/=FULL_GEN)CYCLE                            !<-- Only interested in the domains of the 'full' generation
          KOUNT_DOMS=KOUNT_DOMS+1                                          !<-- Counting domains in the 'full' generation
!
          IF(KOUNT_DOMS==1)THEN
            LEAD_FCST_TASK(ID_DOM)=0                                       !<-- Task 0 is first in line
            LEAD_WRITE_TASK(ID_DOM)=NUM_FCST_TASKS                         !<-- 1st write task follows the last forecast task
          ELSE
            LEAD_FCST_TASK(ID_DOM)=LAST_FCST_TASK_X+1                      !<-- Lead fcst task on domain follows last on previous domain
            LEAD_WRITE_TASK(ID_DOM)=LAST_WRITE_TASK_X+1                    !<-- Lead write task on domain follows last on previous domain
          ENDIF
!
          LAST_FCST_TASK_X=LEAD_FCST_TASK(ID_DOM)+FTASKS_DOMAIN(ID_DOM)-1 
          LAST_FCST_TASK(ID_DOM)=LAST_FCST_TASK_X                          !<-- The last forecast task on this domain
!
          LAST_WRITE_TASK_X=LEAD_WRITE_TASK(ID_DOM)+WTASKS_DOMAIN(ID_DOM)-1 
          LAST_WRITE_TASK(ID_DOM)=LAST_WRITE_TASK_X                        !<-- The last write task on this domain
!
          IF(MYPE>=LEAD_FCST_TASK(ID_DOM)                               &  !<-- 
                       .AND.                                            &  !
             MYPE<=LAST_FCST_TASK(ID_DOM)                               &  !  Associate tasks with each domain.
                       .OR.                                             &  !  Write tasks can be tied to only one domain.
             MYPE>=LEAD_WRITE_TASK(ID_DOM)                              &  !
                       .AND.                                            &  !
             MYPE<=LAST_WRITE_TASK(ID_DOM))THEN                            !<---
!
            MY_DOMAIN_ID_N(ID_DOM)=ID_DOM                                  !<-- This task collects its domain ID in this generation.
          ENDIF
!
          KOUNT_TASKS=0
          DO N2=LEAD_FCST_TASK(ID_DOM),LAST_FCST_TASK(ID_DOM)              !<-- Loop through this domain's fcst tasks.
            KOUNT_TASKS=KOUNT_TASKS+1
            PETLIST_DOM(KOUNT_TASKS,ID_DOM)=N2                             !<-- Insert this fcst task into the domain's task list.
          ENDDO
!
          DO N2=LEAD_WRITE_TASK(ID_DOM),LAST_WRITE_TASK(ID_DOM)            !<-- Loop through this domain's quilt/write tasks.
            KOUNT_TASKS=KOUNT_TASKS+1
            PETLIST_DOM(KOUNT_TASKS,ID_DOM)=N2                             !<-- Insert this write task into the domain's task list.
          ENDDO
!
        ENDDO
!
!-----------------------------------------------------------------------
!***  Now assign the tasks on all the domains in the remaining
!***  generations.  In order to balance the work load as evenly 
!***  as possible the domains in each of these generations will
!***  take their tasks equally from each of the domains in the
!***  full generation.
!-----------------------------------------------------------------------
!
        NDOMS_FULL=DOMS_PER_GEN(FULL_GEN)                                  !<-- # of domains in 1st full generation
!
        ALLOCATE(DOMS_FULL(1:NDOMS_FULL))
        ALLOCATE(FRAC_FULL(1:NDOMS_FULL))
        ALLOCATE(KOUNT_FULL(1:NDOMS_FULL))
!
        DO N=1,NDOMS_FULL
          FRAC_FULL(N)=0.
          KOUNT_FULL(N)=-1
        ENDDO
!
!-----------------------------------------------------------------------
!
        NUM_TASKS_FULL=0
        DO N=1,NUM_DOMAINS_TOTAL
          ID_DOM=RANK_TO_DOMAIN_ID(N)                                      !<-- Domain #N's domain ID (selected by the user)
          IF(DOMAIN_GEN(ID_DOM)/=FULL_GEN)CYCLE
          NUM_TASKS_FULL=NUM_TASKS_FULL+FTASKS_DOMAIN(ID_DOM)              !<-- Add up the # of compute tasks in the full generation.
        ENDDO
!
        RECIP_NUM_TASKS_FULL=1./REAL(NUM_TASKS_FULL)
        KOUNT=0
!
        DO N=1,NUM_DOMAINS_TOTAL
          ID_DOM=RANK_TO_DOMAIN_ID(N)
          IF(DOMAIN_GEN(ID_DOM)/=FULL_GEN)CYCLE                            !<-- Consider only the first full generation
!
          KOUNT=KOUNT+1
          DOMS_FULL(KOUNT)=ID_DOM                                          !<-- IDs of domains in the full generation
          FRAC_FULL(KOUNT)=FTASKS_DOMAIN(ID_DOM)*RECIP_NUM_TASKS_FULL      !<-- Fraction of all tasks on each domain in full generation
          KOUNT_FULL(KOUNT)=PETLIST_DOM(1,ID_DOM)                          !<-- 1st task on each domain of the full generation
        ENDDO
!
!-----------------------------------------------------------------------
!***  In each of the remaining generations fill the domain's compute
!***  tasks proportionately with tasks from each domain in the full
!***  generation in order to spread the work load evenly.
!-----------------------------------------------------------------------
!
        gens_loop: DO N=1,NUM_GENS                                         !<-- Loop through the generations
!
          dom_loop: DO N1=1,NUM_DOMAINS_TOTAL
!
            ID_DOM=RANK_TO_DOMAIN_ID(N1)                                   !<-- Domain ID of domain #N1
            IF(DOMAIN_GEN(ID_DOM)/=N.OR.DOMAIN_GEN(ID_DOM)==FULL_GEN)THEN  !<-- Consider domains in gen #N and not in the full generation
              CYCLE dom_loop
            ENDIF   
!
!-----------------------------------------------------------------------
!***  First assign the write/quilt tasks since they are totally
!***  separate from the fcst/compute tasks and are always assigned
!***  monotonically.
!-----------------------------------------------------------------------
!
            LEAD_WRITE_TASK(ID_DOM)=LAST_WRITE_TASK_X+1                    !<-- Lead write task on domain follows last on previous domain
            LAST_WRITE_TASK(ID_DOM)=LEAD_WRITE_TASK(ID_DOM)             &  !<-- Last write task on this domain
                                   +WTASKS_DOMAIN(ID_DOM)-1
            KOUNT_TASKS=FTASKS_DOMAIN(ID_DOM)                              !<-- Write/quilt tasks follow the compute tasks in the PETList
!
            DO N2=LEAD_WRITE_TASK(ID_DOM),LAST_WRITE_TASK(ID_DOM)
              KOUNT_TASKS=KOUNT_TASKS+1
              PETLIST_DOM(KOUNT_TASKS,ID_DOM)=N2                           !<-- Insert write task into domain's task list
              IF(MYPE==PETLIST_DOM(KOUNT_TASKS,ID_DOM))THEN
                MY_DOMAIN_ID_N(ID_DOM)=ID_DOM                              !<-- This write task collects its domain ID in this generation.
              ENDIF
            ENDDO
            LAST_WRITE_TASK_X=LAST_WRITE_TASK(ID_DOM)
!
!-----------------------------------------------------------------------
!***  Now proceed with the assignment of forecast/compute tasks.
!-----------------------------------------------------------------------
!
            KOUNT_TASKS=0
            DO N2=1,NDOMS_FULL                                             !<-- Loop through the domains in the full generation.
              ID_FULL=DOMS_FULL(N2)                                        !<-- ID of domain #N2 in full generation
              NTASKS_CONTRIB=NINT(FRAC_FULL(N2)*FTASKS_DOMAIN(ID_DOM))     !<-- # of tasks contributed by domain #N2 in full generation
!
              DO N3=1,NTASKS_CONTRIB                                       !<-- Apply the contributed tasks to domain ID_DOM.
                KOUNT_TASKS=KOUNT_TASKS+1
                PETLIST_DOM(KOUNT_TASKS,ID_DOM)=KOUNT_FULL(N2)             !<-- Add this fcst task to this domain's task list.
                IF(MYPE==PETLIST_DOM(KOUNT_TASKS,ID_DOM))THEN
                  MY_DOMAIN_ID_N(ID_DOM)=ID_DOM                            !<-- This fcst task collects its domain ID in this generation.
                ENDIF
!
                KOUNT_FULL(N2)=KOUNT_FULL(N2)+1
                IF(KOUNT_FULL(N2)>PETLIST_DOM(FTASKS_DOMAIN(ID_FULL),ID_FULL))THEN
                  KOUNT_FULL(N2)=PETLIST_DOM(1,ID_FULL)                    !<-- Cycle around contributed tasks from domain ID_FULL.
                ENDIF
!
                IF(KOUNT_TASKS==FTASKS_DOMAIN(ID_DOM))THEN                 !<-- If so, domain ID_DOM has filled its compute tasks.
                  LEAD_FCST_TASK(ID_DOM)=PETLIST_DOM(1,ID_DOM)             !<-- Save identity of this domain's lead compute task.
                  LAST_FCST_TASK(ID_DOM)=PETLIST_DOM(KOUNT_TASKS,ID_DOM)   !<-- Save identity of this domain's last compute task.
                  CYCLE dom_loop
                ENDIF
!
              ENDDO
            ENDDO
!
!-----------------------------------------------------------------------
!***  If we reach this point then all of domain ID_DOM's compute tasks
!***  still have not been assigned.  This is simply due to fractional
!***  roundoff in computing the number of tasks contributed by each
!***  of the domains in the full generation.  Finish assigning this
!***  domain's tasks by taking them from the first domain in the full
!***  generation.
!-----------------------------------------------------------------------
!
            ID_FULL=DOMS_FULL(1)                                           !<-- Take the tasks from the 1st domain in the full generation.
!
            DO N2=KOUNT_TASKS+1,FTASKS_DOMAIN(ID_DOM)
              PETLIST_DOM(N2,ID_DOM)=KOUNT_FULL(1)                         !<-- Add remaining fcst tasks to this domain's task list.
              IF(MYPE==PETLIST_DOM(N2,ID_DOM))THEN
                MY_DOMAIN_ID_N(ID_DOM)=ID_DOM                              !<-- This fcst task collects its domain ID in this generation.
              ENDIF
!
              KOUNT_FULL(1)=KOUNT_FULL(1)+1
              IF(KOUNT_FULL(1)>PETLIST_DOM(FTASKS_DOMAIN(ID_FULL),ID_FULL))THEN
                KOUNT_FULL(1)=PETLIST_DOM(1,ID_FULL)                       !<-- Cycle around contributed tasks from domain ID_FULL.
              ENDIF
            ENDDO
!
            LEAD_FCST_TASK(ID_DOM)=PETLIST_DOM(1,ID_DOM)                     !<-- Save identity of this domain's lead compute task.
            LAST_FCST_TASK(ID_DOM)=PETLIST_DOM(FTASKS_DOMAIN(ID_DOM),ID_DOM) !<-- Save identity of this domain's last compute task.
!
!-----------------------------------------------------------------------
!
          ENDDO dom_loop
!
        ENDDO gens_loop
!
        DEALLOCATE(DOMS_FULL)
        DEALLOCATE(FRAC_FULL)
        DEALLOCATE(KOUNT_FULL)
!
!-----------------------------------------------------------------------
!
      ENDIF  task_assign
!
!-----------------------------------------------------------------------
!***  All tasks know the task counts and IDs of all domains as well as
!***  the parents of each domain.
!
!***  Loop through all domains in order to associate all parents
!***  with their children through intracommunicators.  We cannot use
!***  intercommunicators in general since parent and child domains 
!***  can contain some of the same forecast tasks in 2-way nesting 
!***  and MPI dictates that intercommunicators can only link disjoint
!***  sets of tasks.
!-----------------------------------------------------------------------
!
      ALLOCATE(ID_CHILDREN(1:NUM_DOMAINS_TOTAL,1:NUM_DOMAINS_TOTAL))       !<-- Array to hold all domains' children's IDs
!
      DO N1=1,NUM_DOMAINS_TOTAL
      DO N2=1,NUM_DOMAINS_TOTAL
        ID_CHILDREN(N1,N2)=0                                               !<-- All valid Domain IDs are >0
      ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
!***  Allocate intracommunicators between parents and children for
!***  all of the domains since some forecast tasks may lie on more 
!***  than one parent and/or child domain.  The same is true for 
!***  the lists of ranks of children's local task ranks in the 
!***  intracommunicators.
!-----------------------------------------------------------------------
!
      ALLOCATE(COMMS_DOMAIN(1:NUM_DOMAINS_TOTAL),stat=ISTAT)
      IF(ISTAT/=0)THEN
        WRITE(0,*)' Failed to allocate COMMS_DOMAIN!'
        CALL ESMF_Finalize(rc=RC,endflag=ESMF_END_ABORT)
      ENDIF
!
      DO N=1,NUM_DOMAINS_TOTAL
        comms_domain(N)%TO_PARENT=-999                                     !<-- Initialize to nonsense the intracommunicator to parent
      ENDDO
!
      ALLOCATE(CHILD_RANKS(1:NUM_DOMAINS_TOTAL),stat=ISTAT)
      IF(ISTAT/=0)THEN
        WRITE(0,*)' Failed to allocate CHILD_RANKS!'
        CALL ESMF_Finalize(rc=RC,endflag=ESMF_END_ABORT)
      ENDIF
!
      DO N=1,NUM_DOMAINS_TOTAL
        child_ranks(N)%CHILDREN=>NULL()
      ENDDO
!
!-----------------------------------------------------------------------
!***  Next we need to create MPI groups for the task sets on each
!***  of the domains.
!-----------------------------------------------------------------------
!
      ALLOCATE(GROUP(1:NUM_DOMAINS_TOTAL))                          
!
      CALL MPI_COMM_GROUP(COMM_WORLD                                    &  !<-- Intracommunicator between all tasks in the run
                         ,GROUP_WORLD                                   &  !<-- The MPI group of all tasks in the run
                         ,IERR )
!
      DO N=1,NUM_DOMAINS_TOTAL
        ID_DOM=RANK_TO_DOMAIN_ID(N)
        NTASKS_X=NTASKS_DOMAIN(ID_DOM)                                     !<-- Total # of tasks on domain ID_DOM
!
        CALL MPI_GROUP_INCL(GROUP_WORLD                                 &  !<-- MPI group with all tasks in the run
                           ,NTASKS_X                                    &  !<-- # of fcst tasks on domain ID_DOM
                           ,PETLIST_DOM(1:NTASKS_X,ID_DOM)              &  !<-- The global ranks of tasks that lie on ID_DOM
                           ,GROUP(ID_DOM)                               &  !<-- The new group containing the tasks on ID_DOM
                           ,IERR )
      ENDDO
!
!-----------------------------------------------------------------------
!***  Loop through all domains.  Parent domains will create
!***  intracommunicators with each of their children and vice versa.
!-----------------------------------------------------------------------
!
      main_loop: DO N=1,NUM_DOMAINS_TOTAL
!
!-----------------------------------------------------------------------
!
        ID_DOM=RANK_TO_DOMAIN_ID(N)
!
!-----------------------------------------------------------------------
!
        ID_PARENT=-999                                                     !<-- Initialize to nonsense the parent's domain ID
!
        N_CHILDREN=NUM_CHILDREN(ID_DOM)                                    !<-- The # of children on this domain
!
!-----------------------------------------------------------------------
!
        has_children: IF(N_CHILDREN>0)THEN              
!
!-----------------------------------------------------------------------
!
          ID_PARENT=ID_DOM                                                 !<-- ID_DOM is a parent domain
          NTASKS_PARENT=NTASKS_DOMAIN(ID_PARENT)                           !<-- Total # of fcst and write tasks on this parent domain
!
!-----------------------------------------------------------------------
!***  All domain IDs will be searched to find matches between the
!***  current domain's ID and the parent IDs of the other domains. 
!***  Matches will identify Parent-Child couplets.
!-----------------------------------------------------------------------
!
          NN=0
!
          DO N2=1,NUM_DOMAINS_TOTAL                                        !<-- Search for children who have parent ID_PARENT
            ID_CHILD=ID_DOMAINS(N2)                                        !<-- Check if this domain ID is that of a child
!
            IF(ID_PARENTS(ID_CHILD)==ID_PARENT.AND.ID_PARENT/=-999)THEN    !<-- If yes then we found a nest that is this domain's child
              NN=NN+1                                                      !<-- Increment index of children of the parent domain
              ID_CHILDREN(NN,ID_PARENT)=ID_CHILD                           !<-- IDs of this parent's (ID_PARENT's) children's domains
            ENDIF
!
            IF(NN==N_CHILDREN)THEN                                         !<-- We have found all of this domain's children
              ALLOCATE(comms_domain(ID_PARENT)%TO_CHILDREN(1:N_CHILDREN)  &  !<-- Parent allocates intracommunicators with each child
                      ,stat=ISTAT)
!
              DO N3=1,N_CHILDREN
                comms_domain(ID_PARENT)%TO_CHILDREN(N3)=-999               !<-- Parent initializes intracommunicators with each child
              ENDDO
!
              EXIT
!
            ENDIF
!
          ENDDO
!
!-----------------------------------------------------------------------
!
        ENDIF has_children
!
!-----------------------------------------------------------------------
!***  Now create groups that are unions of each parent with each
!***  of their children.  From those unions create the final 
!***  parent-child intracommunicators.
!-----------------------------------------------------------------------
!
        IF(N_CHILDREN>0)THEN
!
          ALLOCATE(child_ranks(ID_PARENT)%CHILDREN(1:N_CHILDREN)        &
                  ,stat=ISTAT)
!
          IF(ISTAT/=0)THEN
            WRITE(0,*)' Failed to allocate child_ranks%CHILDREN!'
            CALL ESMF_Finalize(rc=RC,endflag=ESMF_END_ABORT)
          ENDIF
!
!-----------------------------------------------------------------------
!
          intra_comm: DO N2=1,N_CHILDREN                                   !<-- Loop through the given parent's children
!
!-----------------------------------------------------------------------
!
            ID_CHILD=ID_CHILDREN(N2,ID_PARENT)                             !<-- Domain ID of child N2 of domain ID_PARENT
!
            CALL MPI_GROUP_UNION(GROUP(ID_PARENT)                       &  !<-- The group containing the parent tasks (in)
                                ,GROUP(ID_CHILD)                        &  !<-- The group containing the child tasks (in)
                                ,GROUP_UNION                            &  !<-- The union of the parent and child groups (out)
                                ,IERR )
!
            CALL MPI_COMM_CREATE(COMM_WORLD                             &  !<-- Intracommunicator between all tasks in the run (in)
                                ,GROUP_UNION                            &  !<-- The union of the parent and child groups (in)
                                ,COMM_INTRA                             &  !<-- Intracommunicator between tasks in the union (out)
                                ,IERR )
!
            comms_domain(ID_PARENT)%TO_CHILDREN(N2)=COMM_INTRA             !<-- Parent: The intracommunicator with its child N2
            comms_domain(ID_CHILD)%TO_PARENT=COMM_INTRA                    !<-- Child: The intracommunicator with its parent
!
!-----------------------------------------------------------------------
!***  The parent's tasks were listed first in the creation of the union
!***  with child tasks so the parent task ranks in the union go from
!***  0 to FTASKS_DOMAIN(ID_PARENT)-1.  However the child task ranks
!***  in the union can be rather jumbled depending on how they overlie
!***  the parent tasks.  Therefore parents must store the union ranks
!***  of their children's forecast tasks in order to use the
!***  intracommunicators.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  First we need to produce the list of global parent and child
!***  tasks equivalent to that which MPI produces but is not seen
!***  when the union of the parent and child groups is created.
!-----------------------------------------------------------------------
!
            NMAX=NTASKS_DOMAIN(ID_PARENT)+NTASKS_DOMAIN(ID_CHILD)          !<-- Max # of tasks that can be in parent-child union
            ALLOCATE(GLOBAL_UNION(1:NMAX))                                 !<-- For holding global ranks in the union
!
            DO N3=1,NTASKS_DOMAIN(ID_PARENT)
              GLOBAL_UNION(N3)=PETLIST_DOM(N3,ID_PARENT)                   !<-- Insert parent's global task ranks into list first
            ENDDO
!
            KOUNT=NTASKS_DOMAIN(ID_PARENT)                                 !<-- We just inserted this many values into GLOBAL_UNION
!
            child_loop1: DO N3=1,NTASKS_DOMAIN(ID_CHILD)                   !<-- Now loop through all child tasks
!
              DO NN=1,NTASKS_DOMAIN(ID_PARENT)                             !<-- Compare against the parent task ranks
                IF(PETLIST_DOM(N3,ID_CHILD)==GLOBAL_UNION(NN))THEN
                  CYCLE child_loop1                                        !<-- No task rank can appear twice in the union
                ENDIF
              ENDDO
!
              KOUNT=KOUNT+1                                                !<-- Accumulating # of unique global task ranks in the union
              GLOBAL_UNION(KOUNT)=PETLIST_DOM(N3,ID_CHILD)                 !<-- Add this child global rank to the union list
!
            ENDDO child_loop1
!
!-----------------------------------------------------------------------
!***  The GLOBAL_UNION array now holds the union of the parent and
!***  child global task ranks with no ranks appearing more than once.
!***  Now the parent creates a list of the child's tasks in the union
!***  but using ranks as they exist in the intracommunicator which
!***  start with 0 for the parent's lead task and simply increase
!***  one by one in a monotonic sequence.
!-----------------------------------------------------------------------
!
            ALLOCATE(child_ranks(ID_PARENT)%CHILDREN(N2)%DATA(0:NTASKS_DOMAIN(ID_CHILD)-1))  !<-- Local ranks of child N2's tasks
!                                                                                                 in parent-child intracomm
            child_loop2: DO N3=0,NTASKS_DOMAIN(ID_CHILD)-1
!
              DO NN=1,KOUNT                                                !<-- Loop through all global task ranks in the union list
                IF(PETLIST_DOM(N3+1,ID_CHILD)==GLOBAL_UNION(NN))THEN       !<-- Search for child task N3's global rank in the union list
                  child_ranks(ID_PARENT)%CHILDREN(N2)%DATA(N3)=NN-1        !<-- Save its local rank in the parent-child intracommunicator
                  CYCLE child_loop2
                ENDIF
                IF(NN<KOUNT)CYCLE
!
                WRITE(0,*)' ERROR: Child global task rank not found'
                WRITE(0,*)' ABORTING!'
                CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
              ENDDO
!
            ENDDO child_loop2
!
            DEALLOCATE(GLOBAL_UNION)
!
!-----------------------------------------------------------------------
!
            CALL MPI_BARRIER(COMM_WORLD,IERR)
!
!-----------------------------------------------------------------------
!
          ENDDO intra_comm
!
        ENDIF
!
!-----------------------------------------------------------------------
!
      ENDDO main_loop
!
!-----------------------------------------------------------------------
!
      DEALLOCATE(LEAD_TASK)
      DEALLOCATE(LAST_TASK)
      DEALLOCATE(GROUP)
      IF(ALLOCATED(DOMS_PER_GEN))DEALLOCATE(DOMS_PER_GEN)
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE PARENT_CHILD_COMMS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE PARENT_TO_CHILD_INIT_NMM(MYPE                          &
                                         ,CF                            &
                                         ,MY_DOMAIN_ID                  &
                                         ,THIS_CHILD_ID                 &
                                         ,SOLVER_GRID_COMP              &
                                         ,COMM_MY_DOMAIN )
!
!-----------------------------------------------------------------------
!***  Initialize data in a child's domain with data from its parent.
!***  Rows and columns on the child's grid lie parallel to rows and
!***  colums on the parent's grid.
!***  Only parent tasks are needed for this.
!-----------------------------------------------------------------------
!
      USE MODULE_SOLVER_INTERNAL_STATE,ONLY: SOLVER_INTERNAL_STATE  &
                                            ,WRAP_SOLVER_INT_STATE
!
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      INTEGER(kind=KINT),INTENT(IN) :: MYPE                             &  !<-- My MPI task rank
                                      ,MY_DOMAIN_ID                     &  !<-- IDs of each domain
                                      ,THIS_CHILD_ID                    &  !<-- ID of the current child
                                      ,COMM_MY_DOMAIN                      !<-- MPI intracommunicator for individual domains
!
      TYPE(ESMF_Config),DIMENSION(*),INTENT(INOUT) :: CF                   !<-- The config objects (one per domain)
!
      TYPE(ESMF_GridComp),INTENT(INOUT)          :: SOLVER_GRID_COMP       !<-- The parent's Solver gridded component
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: I,J,L,M,N,NCHILD,NFCST,RC,RC_CHILD,LNSH=1    ! Tom, change lnsh value later
      INTEGER(kind=KINT) :: INPES,JNPES
      INTEGER(kind=KINT) :: IDS,IDE,JDS,JDE
      INTEGER(kind=KINT) :: IMS,IME,JMS,JME
      INTEGER(kind=KINT) :: ITS,ITE,JTS,JTE
      INTEGER(kind=KINT) :: NUM_PES_PARENT
      INTEGER(kind=KINT) :: IM_CHILD,JM_CHILD
      INTEGER(kind=KINT) :: IM_PARENT,JM_PARENT
      INTEGER(kind=KINT) :: I_PARENT_START,J_PARENT_START
      INTEGER(kind=KINT) :: IHREND,NLEV,NTSD
      INTEGER(kind=KINT) :: PARENT_CHILD_SPACE_RATIO
!
      INTEGER(kind=KINT),DIMENSION(:),POINTER :: LOCAL_ISTART           &
                                                ,LOCAL_IEND             &
                                                ,LOCAL_JSTART           &
                                                ,LOCAL_JEND
!
      INTEGER(kind=KINT),ALLOCATABLE,DIMENSION(:,:) :: IDUMMY_2D
!
      REAL(kind=KFPT) :: CHILD_PARENT_SPACE_RATIO                       &
                        ,COL_0                                          &
                        ,DLMD                                           &
                        ,DLMD_CHILD                                     &
                        ,DLMD_PARENT                                    &
                        ,DPHD                                           &
                        ,DPHD_CHILD                                     &
                        ,DPHD_PARENT                                    &
                        ,ROW_0                                          &
                        ,SBD                                            &
                        ,SBD_CHILD                                      &
                        ,SBD_PARENT                                     &
                        ,SW_LATD_CHILD                                  &
                        ,SW_LOND_CHILD                                  &
                        ,TLM0D_CHILD                                    &
                        ,TLM0D_PARENT                                   &
                        ,TPH0D_CHILD                                    &
                        ,TPH0D_PARENT                                   &
                        ,WBD                                            &
                        ,WBD_CHILD                                      &
                        ,WBD_PARENT
!
      REAL(kind=KFPT),ALLOCATABLE,DIMENSION(:) :: DUMMY_SOIL
!
      REAL(kind=KFPT),ALLOCATABLE,DIMENSION(:,:) :: SEA_MASK,SEA_ICE
!
      REAL(kind=KFPT),ALLOCATABLE,DIMENSION(:,:,:) :: DUMMY_2D_IN       &
                                                     ,DUMMY_2D_OUT      &
                                                     ,DUMMY_3D          &
                                                     ,DUMMY_3DS         &
                                                     ,PD_BILINEAR       &
                                                     ,PD_NEAREST        &
                                                     ,TEMPSOIL
!
      CHARACTER(len=2)  :: INT_TO_CHAR
      CHARACTER(len=6)  :: FMT
      CHARACTER(len=50) :: GLOBAL_FLAG                                  &
                          ,OUTFILE    
!
      LOGICAL(kind=KLOG) :: GLOBAL,OPENED
      LOGICAL(kind=KLOG),ALLOCATABLE,DIMENSION(:,:) :: LOWER_TOPO
!
      TYPE(WRAP_SOLVER_INT_STATE)  :: WRAP_SOLVER
!
      TYPE(SOLVER_INTERNAL_STATE),POINTER :: SOLVER_INT_STATE
!
!-----------------------------------------------------------------------
!***  This routine provides data to a child domain when no normal
!***  pre-processed input has been prepared.  This is done by simple
!***  bilinear interpolation from the parent domain's data.
!
!***  The record order of the parent's input is duplicated and those 
!***  records are written to a file so that the child can read in and 
!***  distribute the data just as if a normal pre-processed file were 
!***  being used.
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC      =ESMF_SUCCESS
      RC_CHILD=ESMF_SUCCESS
!
      NULLIFY(solver_int_state)
!
!-----------------------------------------------------------------------
!***  Only forecast tasks are relevant and have correct data
!***  for this work.  Find the number of forecast tasks from
!***  configure file (all forecast tasks already have this
!***  in their internal state but the write tasks do not).
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Child_Init: Extract INPES From Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)              &  !<-- The parent's config object
                                  ,value =INPES                         &  !<-- The variable filled (parent fcst tasks in I)
                                  ,label ='inpes:'                      &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CHILD)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Child_Init: Extract JNPES From Config File"
!!!   CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)              &  !<-- The parent's config object
                                  ,value =JNPES                         &  !<-- The variable filled (parent fcst tasks in J)
                                  ,label ='jnpes:'                      &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CHILD)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
      NUM_PES_PARENT=INPES*JNPES                                           !<-- # of parent fcst tasks
!
      IF(MYPE>=NUM_PES_PARENT)RETURN                                       !<-- Parent's quilt/write tasks may leave
!
!-----------------------------------------------------------------------
!***  We need the spatial resolution of the parent grid so extract
!***  its dimensions and its bounds.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Child_Init: Extract IM From Config File"
!!!   CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)              &  !<-- The parent's config object
                                  ,value =IM_PARENT                     &  !<-- The variable filled (I dimension of parent grid)
                                  ,label ='im:'                         &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CHILD)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Child_Init: Extract JM From Config File"
!!!   CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)              &  !<-- The parent's config object
                                  ,value =JM_PARENT                     &  !<-- The variable filled (J dimension of parent grid)
                                  ,label ='jm:'                         &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CHILD)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Child_Init: Extract SBD From Config File"
!!!   CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!!!!  CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)              &  !<-- The parent's config object
!!!!                              ,value =SBD_PARENT                    &  !<-- The variable filled (South boundary of parent grid)
!!!!                              ,label ='sbd:'                        &  !<-- Give this label's value to the previous variable
!!!!                              ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CHILD)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Child_Init: Extract WBD From Config File"
!!!   CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!!!!  CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)              &  !<-- The parent's config object
!!!!                              ,value =WBD_PARENT                    &  !<-- The variable filled (West boundary of parent grid)
!!!!                              ,label ='wbd:'                        &  !<-- Give this label's value to the previous variable
!!!!                              ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CHILD)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Child_Init: Extract TPH0D From Config File"
!!!   CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!!!!  CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)              &  !<-- The parent's config object
!!!!                              ,value =TPH0D_PARENT                  &  !<-- The variable filled (Central lat of parent grid)
!!!!                              ,label ='tph0d:'                      &  !<-- Give this label's value to the previous variable
!!!!                              ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CHILD)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Child_Init: Extract TLM0D From Config File"
!!!   CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!!!!  CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)              &  !<-- The parent's config object
!!!!                              ,value =TLM0D_PARENT                  &  !<-- The variable filled (Central lon of parent grid)
!!!!                              ,label ='tlm0d:'                      &  !<-- Give this label's value to the previous variable
!!!!                              ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CHILD)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Global Flag for Parent Domain"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)             &  !<-- The configure object of my parent
                                  ,value =GLOBAL_FLAG                  &  !<-- The variable filled 
                                  ,label ='global:'                    &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CHILD) 
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  The parent grid resolution.
!-----------------------------------------------------------------------
!
      IF(TRIM(GLOBAL_FLAG)=='true')THEN                                    !<-- Parent is global 
        GLOBAL=.TRUE.
        IDE=IM_PARENT+2
        JDE=JM_PARENT+2
!!!!    DPHD_PARENT=-SBD_PARENT*2./REAL(JDE-3)
!!!!    DLMD_PARENT=-WBD_PARENT*2./REAL(IDE-3)
      ELSE                                                                 !<-- Parent is regional
        GLOBAL=.FALSE.
        IDE=IM_PARENT
        JDE=JM_PARENT
!!!!    DPHD_PARENT=-SBD_PARENT*2./REAL(JDE-1)
!!!!    DLMD_PARENT=-WBD_PARENT*2./REAL(IDE-1)
      ENDIF
!
      ROW_0=0.5*(JDE+1)
      COL_0=0.5*(IDE+1)
!
!-----------------------------------------------------------------------
!***  Extract the Solver internal state of the parent
!***  so we can use their data for the nests.
!-----------------------------------------------------------------------
!
      CALL ESMF_GridCompGetInternalState(SOLVER_GRID_COMP               &
                                        ,WRAP_SOLVER                    &
                                        ,RC )
!
!-----------------------------------------------------------------------
!
      SOLVER_INT_STATE=>wrap_solver%INT_STATE
!
      IMS=solver_int_state%IMS                                                !<-- Horizontal memory limits on parent tasks
      IME=solver_int_state%IME                                                !
      JMS=solver_int_state%JMS                                                !
      JME=solver_int_state%JME                                                !<--
!
      ITS=solver_int_state%ITS                                                !<-- Horizontal integration limits on parent tasks
      ITE=solver_int_state%ITE                                                !
      JTS=solver_int_state%JTS                                                !
      JTE=solver_int_state%JTE                                                !<--
!
      LM=solver_int_state%LM                                                  !<-- Number of atmospheric layers
!
      LOCAL_ISTART=>solver_int_state%LOCAL_ISTART                             !<-- Local integration limits for all parent tasks
      LOCAL_IEND  =>solver_int_state%LOCAL_IEND                               !
      LOCAL_JSTART=>solver_int_state%LOCAL_JSTART                             !
      LOCAL_JEND  =>solver_int_state%LOCAL_JEND                               !<--
!
!-----------------------------------------------------------------------
!***  DPHD/DLMD and SBD/WBD are used only for stand-alone, independent
!***  rotated parent/nest grids (i.e., not grid-associated nests).
!-----------------------------------------------------------------------
!
      DPHD_PARENT=solver_int_state%DPHD
      DLMD_PARENT=solver_int_state%DLMD
      SBD_PARENT=solver_int_state%SBD
      WBD_PARENT=solver_int_state%WBD
      TPH0D_PARENT=solver_int_state%TPH0D
      TLM0D_PARENT=solver_int_state%TLM0D
!
!-----------------------------------------------------------------------
!***  Extract relevant information from this child's configure file.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Child_Init: Extract I of SW Point on Parent"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF(THIS_CHILD_ID)             &  !<-- The child's config object
                                  ,value =I_PARENT_START                &  !<-- The variable filled (parent I of child's SW corner)
                                  ,label ='i_parent_start:'             &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CHILD)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Child_Init: Extract J of SW Point on Parent"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF(THIS_CHILD_ID)             &  !<-- The child's config object
                                  ,value =J_PARENT_START                &  !<-- The variable filled (parent J of child's SW corner)
                                  ,label ='j_parent_start:'             &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CHILD)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Child_Init: Extract Child/Parent Grid Ratio"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF(THIS_CHILD_ID)             &  !<-- The child's config object
                                  ,value =PARENT_CHILD_SPACE_RATIO      &  !<-- The variable filled (child grid increment / parent's)
                                  ,label ='parent_child_space_ratio:'   &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CHILD)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CHILD_PARENT_SPACE_RATIO=1./REAL(PARENT_CHILD_SPACE_RATIO)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Child_Init: Extract Global IM of Child"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF(THIS_CHILD_ID)             &  !<-- The child's config object
                                  ,value =IM_CHILD                      &  !<-- The variable filled (IM of child domain)
                                  ,label ='im:'                         &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CHILD)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Child_Init: Extract Global JM of Child"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF(THIS_CHILD_ID)             &  !<-- The child's config object
                                  ,value =JM_CHILD                      &  !<-- The variable filled (JM of child domain)
                                  ,label ='jm:'                         &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CHILD)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Only for free-standing nests:
!
!***  What is the parent lat/lon of the SW corner H point of the
!***  child grid?
!***  Find the resolution, bounds, and center of the child grid.
!-----------------------------------------------------------------------
!
!     CALL CONVERT_IJ_TO_LATLON  (I_PARENT_START                        &
!                                ,J_PARENT_START                        &
!                                ,IM_PARENT                             &
!                                ,JM_PARENT                             &
!                                ,TPH0D_PARENT                          &
!                                ,TLM0D_PARENT                          &
!                                ,DPHD_PARENT                           &
!                                ,DLMD_PARENT                           &
!                                ,SW_LATD_CHILD                         &
!                                ,SW_LOND_CHILD )
!
!!!   DPHD_CHILD=DPHD_PARENT*CHILD_PARENT_SPACE_RATIO
!!!   DLMD_CHILD=DLMD_PARENT*CHILD_PARENT_SPACE_RATIO
!
!!!   SBD_CHILD=-0.5*(JM_CHILD-1)*DPHD_CHILD
!!!   WBD_CHILD=-0.5*(IM_CHILD-1)*DLMD_CHILD
!
!!!   CALL CENTER_NEST(SBD_CHILD                                        &
!!!                   ,WBD_CHILD                                        &
!!!                   ,SW_LATD_CHILD                                    &
!!!                   ,SW_LOND_CHILD                                    &
!!!                   ,TPH0D_CHILD                                      &
!!!                   ,TLM0D_CHILD )
!
!-----------------------------------------------------------------------
!***  Allocate 2-D and 3-D dummy arrays for child quantities.
!-----------------------------------------------------------------------
!
      ALLOCATE(SEA_MASK(1:IM_CHILD,1:JM_CHILD))
      ALLOCATE(SEA_ICE(1:IM_CHILD,1:JM_CHILD))
!
      ALLOCATE(IDUMMY_2D(1:IM_CHILD,1:JM_CHILD))
      ALLOCATE(DUMMY_2D_IN (IMS:IME,JMS:JME,1:1))
      ALLOCATE(DUMMY_2D_OUT(1:IM_CHILD,1:JM_CHILD,1:1))
      ALLOCATE(DUMMY_3D (1:IM_CHILD,1:JM_CHILD,1:LM))
      ALLOCATE(DUMMY_3DS(1:IM_CHILD,1:JM_CHILD,1:NUM_SOIL_LAYERS))
      ALLOCATE(DUMMY_SOIL(1:NUM_SOIL_LAYERS))
      ALLOCATE(TEMPSOIL (1:NUM_SOIL_LAYERS,1:IM_CHILD,1:JM_CHILD))
!
      ALLOCATE(PD_NEAREST (1:IM_CHILD,1:JM_CHILD,1:1))
      ALLOCATE(PD_BILINEAR(1:IM_CHILD,1:JM_CHILD,1:1))
!
      ALLOCATE(LOWER_TOPO(1:IM_CHILD,1:JM_CHILD))
      DO J=1,JM_CHILD
      DO I=1,IM_CHILD
        LOWER_TOPO(I,J)=.FALSE.
      ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
!***  Parent task 0 opens a file for writing out the child's input.
!-----------------------------------------------------------------------
!
      IF(MYPE==0)THEN
!
        select_unit: DO N=51,59
          INQUIRE(N,OPENED=OPENED)
          IF(.NOT.OPENED)THEN
            NFCST=N
            EXIT select_unit
          ENDIF
        ENDDO select_unit
!
        FMT='(I2.2)'
        WRITE(INT_TO_CHAR,FMT)THIS_CHILD_ID
        OUTFILE='input_domain_'//INT_TO_CHAR
!
        OPEN(unit=NFCST,file=OUTFILE,status='replace',form='unformatted')
!
!-----------------------------------------------------------------------
!***  The following variables are for the vertical grid structure
!***  and are shared by the parent and its children.
!-----------------------------------------------------------------------
!
        IHREND=0                                                           !<-- Not used 
        NTSD  =0                                                           !<-- Not used
!
        WRITE(NFCST)solver_int_state%RUN                                   &
                   ,solver_int_state%IDAT                                  &
                   ,solver_int_state%IHRST                                 &
!                  ,solver_int_state%IHREND                                &
!                  ,solver_int_state%NTSD
                   ,IHREND                                              &
                   ,NTSD
!
        WRITE(NFCST)solver_int_state%PT                                    &
                   ,solver_int_state%PDTOP                                 &
                   ,solver_int_state%LPT2                                  &
                   ,solver_int_state%SGM                                   &
                   ,solver_int_state%SG1                                   &
                   ,solver_int_state%DSG1                                  &
                   ,solver_int_state%SGML1                                 &
                   ,solver_int_state%SG2                                   &
                   ,solver_int_state%DSG2                                  &
                   ,solver_int_state%SGML2
!
        WRITE(NFCST)I_PARENT_START,J_PARENT_START
!
        DLMD=DLMD_PARENT*CHILD_PARENT_SPACE_RATIO
        DPHD=DPHD_PARENT*CHILD_PARENT_SPACE_RATIO
!
        IF(GLOBAL)THEN
          SBD=SBD_PARENT+(J_PARENT_START-2)*DPHD_PARENT
          WBD=WBD_PARENT+(I_PARENT_START-2)*DLMD_PARENT
        ELSE 
          SBD=SBD_PARENT+(J_PARENT_START-1)*DPHD_PARENT
          WBD=WBD_PARENT+(I_PARENT_START-1)*DLMD_PARENT
        ENDIF
!
        WRITE(NFCST)DLMD,DPHD                                           &
                   ,WBD,SBD                                             &
                   ,TLM0D_PARENT,TPH0D_PARENT
!
        WRITE(NFCST)IM_CHILD,JM_CHILD,LM,LNSH
!
!-----------------------------------------------------------------------
!
      ENDIF
!
      NLEV=1
!
!-----------------------------------------------------------------------
!***  Sea Mask
!-----------------------------------------------------------------------
!***  The Sea Mask is needed for the Sfc Geopotential so compute it now.
!***  If there are adjacent water points with different elevations
!***  after Sfc Geopotential is computed then the WATERFALL routine
!***  will level them by changing the sfc elevations.  At such points
!***  the atmospheric column will need adjusting so save the locations
!***  of those points along with the preliminary values of the nest's
!***  PD, T, Q, CW, U, and V which will then be modified.  
!***  The Sea Mask will be written out in its proper place following
!***  the Stnd Deviation of Sfc Height.  
!-----------------------------------------------------------------------
!
      DO J=JMS,JME
      DO I=IMS,IME
        DUMMY_2D_IN(I,J,1)=solver_int_state%SM(I,J)
      ENDDO
      ENDDO
!
      CALL PARENT_TO_CHILD_FILL(DUMMY_2D_IN                             &
                               ,NLEV                                    &
                               ,'SeaMask'                               &
                               ,DUMMY_2D_OUT                            &
                               ,NEAREST)
!
      IF(MYPE==0)THEN
        DO J=1,JM_CHILD
        DO I=1,IM_CHILD
          SEA_MASK(I,J)=REAL(NINT(DUMMY_2D_OUT(I,J,1)))
        ENDDO
        ENDDO
      ENDIF
!
!-----------------------------------------------------------------------
!***  Sfc Geopotential
!-----------------------------------------------------------------------
!
      DO J=JMS,JME
      DO I=IMS,IME
        DUMMY_2D_IN(I,J,1)=solver_int_state%FIS(I,J)
      ENDDO
      ENDDO
!
      CALL PARENT_TO_CHILD_FILL(DUMMY_2D_IN                             &
                               ,NLEV                                    &
                               ,'FIS'                                   &
                               ,DUMMY_2D_OUT                            &
                               ,BILINEAR)
!
      IF(MYPE==0)THEN
        CALL WATERFALLS(DUMMY_2D_OUT                                    &  !<-- Level adjacent water points with different elevations
                       ,SEA_MASK                                        &
                       ,LOWER_TOPO                                      &
                       ,1,IM_CHILD,1,JM_CHILD)
!
        WRITE(NFCST)DUMMY_2D_OUT
      ENDIF
!
!     write(0,*)' after Sfc Geo'
!
!-----------------------------------------------------------------------
!***  Stnd Deviation of Sfc Height
!-----------------------------------------------------------------------

      DO J=JMS,JME
      DO I=IMS,IME
        DUMMY_2D_IN(I,J,1)=solver_int_state%STDH(I,J)
      ENDDO
      ENDDO
!
      CALL PARENT_TO_CHILD_FILL(DUMMY_2D_IN                             &
                               ,NLEV                                    &
                               ,'STDH'                                  &
                               ,DUMMY_2D_OUT                            &
                               ,BILINEAR)
!
      IF(MYPE==0)WRITE(NFCST)DUMMY_2D_OUT
!     write(0,*)' after STDH'
!
!-----------------------------------------------------------------------
!***  Sea Mask
!-----------------------------------------------------------------------
!
      IF(MYPE==0)THEN
        WRITE(NFCST)SEA_MASK
      ENDIF
!     write(0,*)' after Sea Mask'
!
!-----------------------------------------------------------------------
!***  PD
!-----------------------------------------------------------------------
!
      DO J=JMS,JME
      DO I=IMS,IME
        DUMMY_2D_IN(I,J,1)=solver_int_state%PD(I,J)
      ENDDO
      ENDDO
!
      CALL PARENT_TO_CHILD_FILL(DUMMY_2D_IN                             &
                               ,NLEV                                    &
                               ,'PD'                                    &
                               ,PD_BILINEAR                             &
                               ,BILINEAR)
!
      CALL PARENT_TO_CHILD_FILL(DUMMY_2D_IN                             &  !<-- Save nearest neighbors for topo adjustment
                               ,NLEV                                    &
                               ,'PD'                                    &
                               ,PD_NEAREST                              &
                               ,NEAREST) 
!
      IF(MYPE==0)THEN
        DO J=1,JM_CHILD
        DO I=1,IM_CHILD
          IF(LOWER_TOPO(I,J))THEN
            DUMMY_2D_OUT(I,J,1)=PD_NEAREST(I,J,1)
          ELSE
            DUMMY_2D_OUT(I,J,1)=PD_BILINEAR(I,J,1)
          ENDIF
        ENDDO
        ENDDO
      ENDIF
!
!     write(0,*)' after PD'
!
      IF(MYPE==0)WRITE(NFCST)DUMMY_2D_OUT
!
!-----------------------------------------------------------------------
!***  U
!-----------------------------------------------------------------------
!
      CALL PARENT_TO_CHILD_FILL(solver_int_state%U, LM                     &
                               ,'Uwind'                                 &
                               ,DUMMY_3D                                &
                               ,BILINEAR)
!
      IF(MYPE==0)THEN
        CALL ADJUST_COLUMNS(PD_NEAREST                                  &
                           ,PD_BILINEAR                                 &
                           ,LOWER_TOPO                                  &
                           ,DUMMY_3D                                    &
                           ,solver_int_state%PT                            &
                           ,solver_int_state%PDTOP                         &
                           ,solver_int_state%SG1                           &
                           ,solver_int_state%SG2                           &
                           ,IM_CHILD,JM_CHILD)
      ENDIF
!
      DO L=1,LM
        IF(MYPE==0)WRITE(NFCST)DUMMY_3D(:,:,L)
      ENDDO
!     write(0,*)' after U'
!
!-----------------------------------------------------------------------
!***  V
!-----------------------------------------------------------------------
!
      CALL PARENT_TO_CHILD_FILL(solver_int_state%V, LM                     &
                               ,'Vwind'                                 &
                               ,DUMMY_3D                                &
                               ,BILINEAR)
!
      IF(MYPE==0)THEN
        CALL ADJUST_COLUMNS(PD_NEAREST                                  &
                           ,PD_BILINEAR                                 &
                           ,LOWER_TOPO                                  &
                           ,DUMMY_3D                                    &
                           ,solver_int_state%PT                            &
                           ,solver_int_state%PDTOP                         &
                           ,solver_int_state%SG1                           &
                           ,solver_int_state%SG2                           &
                           ,IM_CHILD,JM_CHILD)
      ENDIF
!
      DO L=1,LM
        IF(MYPE==0)WRITE(NFCST)DUMMY_3D(:,:,L)
      ENDDO
!     write(0,*)' after V'
!
!-----------------------------------------------------------------------
!***  T
!-----------------------------------------------------------------------
!
      CALL PARENT_TO_CHILD_FILL(solver_int_state%T, LM                     &
                               ,'Temperature'                           &
                               ,DUMMY_3D                                &
                               ,BILINEAR)
!
      IF(MYPE==0)THEN
        CALL ADJUST_COLUMNS(PD_NEAREST                                  &
                           ,PD_BILINEAR                                 &
                           ,LOWER_TOPO                                  &
                           ,DUMMY_3D                                    &
                           ,solver_int_state%PT                            &
                           ,solver_int_state%PDTOP                         &
                           ,solver_int_state%SG1                           &
                           ,solver_int_state%SG2                           &
                           ,IM_CHILD,JM_CHILD)
      ENDIF
!
      DO L=1,LM
        IF(MYPE==0)WRITE(NFCST)DUMMY_3D(:,:,L)
      ENDDO
!     write(0,*)' after T'
!
!-----------------------------------------------------------------------
!***  Q
!-----------------------------------------------------------------------
!
      CALL PARENT_TO_CHILD_FILL(solver_int_state%Q, LM                     &
                               ,'SpecHum'                               &
                               ,DUMMY_3D                                &
                               ,BILINEAR)
!
      IF(MYPE==0)THEN
        CALL ADJUST_COLUMNS(PD_NEAREST                                  &
                           ,PD_BILINEAR                                 &
                           ,LOWER_TOPO                                  &
                           ,DUMMY_3D                                    &
                           ,solver_int_state%PT                            &
                           ,solver_int_state%PDTOP                         &
                           ,solver_int_state%SG1                           &
                           ,solver_int_state%SG2                           &
                           ,IM_CHILD,JM_CHILD)
      ENDIF
!
      DO L=1,LM
        IF(MYPE==0)WRITE(NFCST)DUMMY_3D(:,:,L)
      ENDDO
!     write(0,*)' after Q'
!
!-----------------------------------------------------------------------
!***  CW
!-----------------------------------------------------------------------
!
      CALL PARENT_TO_CHILD_FILL(solver_int_state%CW, LM                    &
                               ,'CW'                                    &
                               ,DUMMY_3D                                &
                               ,BILINEAR)
!
      IF(MYPE==0)THEN
        CALL ADJUST_COLUMNS(PD_NEAREST                                  &
                           ,PD_BILINEAR                                 &
                           ,LOWER_TOPO                                  &
                           ,DUMMY_3D                                    &
                           ,solver_int_state%PT                            &
                           ,solver_int_state%PDTOP                         &
                           ,solver_int_state%SG1                           &
                           ,solver_int_state%SG2                           &
                           ,IM_CHILD,JM_CHILD)
      ENDIF
!
      DO L=1,LM
        IF(MYPE==0)WRITE(NFCST)DUMMY_3D(:,:,L)
      ENDDO
!     write(0,*)' after CW'
!
!-----------------------------------------------------------------------
!***  O3
!-----------------------------------------------------------------------
!
!     CALL PARENT_TO_CHILD_FILL(solver_int_state%O3, LM                    &
!                              ,'O3'                                    &
!                              ,DUMMY_3D                                &
!                              ,BILINEAR)
!
!     IF(MYPE==0)THEN
!       CALL ADJUST_COLUMNS(PD_NEAREST                                  &
!                          ,PD_BILINEAR                                 &
!                          ,LOWER_TOPO                                  &
!                          ,DUMMY_3D                                    &
!                          ,solver_int_state%PT                            &
!                          ,solver_int_state%PDTOP                         &
!                          ,solver_int_state%SG1                           &
!                          ,solver_int_state%SG2                           &
!                          ,IM_CHILD,JM_CHILD)
!     ENDIF
!
      DO L=1,LM
        DO J=1,JM_CHILD
        DO I=1,IM_CHILD
          DUMMY_3D(I,J,L)=0.    ! for now keep O3 = 0.
        ENDDO
        ENDDO
        IF(MYPE==0)WRITE(NFCST)DUMMY_3D(:,:,L)
      ENDDO
!     write(0,*)' after O3'
!
!-----------------------------------------------------------------------
!***  ALBEDO
!-----------------------------------------------------------------------
!
      DO J=JMS,JME
      DO I=IMS,IME
        DUMMY_2D_IN(I,J,1)=solver_int_state%ALBEDO(I,J)
      ENDDO
      ENDDO
!
      CALL PARENT_TO_CHILD_FILL(DUMMY_2D_IN                             &
                               ,NLEV                                    &
                               ,'ALBEDO'                                &
                               ,DUMMY_2D_OUT                            &
                               ,BILINEAR)
!
      IF(MYPE==0)WRITE(NFCST)DUMMY_2D_OUT
!     write(0,*)' after Albedo'
!
!-----------------------------------------------------------------------
!***  ALBASE
!-----------------------------------------------------------------------
!
      DO J=JMS,JME
      DO I=IMS,IME
        DUMMY_2D_IN(I,J,1)=solver_int_state%ALBASE(I,J)
      ENDDO
      ENDDO
!
      CALL PARENT_TO_CHILD_FILL(DUMMY_2D_IN, 1                          &
                               ,'ALBASE'                                &
                               ,DUMMY_2D_OUT                            &
                               ,BILINEAR)
!
      IF(MYPE==0)WRITE(NFCST)DUMMY_2D_OUT
!
!-----------------------------------------------------------------------
!***  EPSR
!-----------------------------------------------------------------------
!
      DO J=JMS,JME
      DO I=IMS,IME
        DUMMY_2D_IN(I,J,1)=solver_int_state%EPSR(I,J)
      ENDDO
      ENDDO
!
      CALL PARENT_TO_CHILD_FILL(DUMMY_2D_IN                             &
                               ,NLEV                                    &
                               ,'EPSR'                                  &
                               ,DUMMY_2D_OUT                            &
                               ,BILINEAR)
!
      IF(MYPE==0)WRITE(NFCST)DUMMY_2D_OUT
!     write(0,*)' after EPSR'
!
!-----------------------------------------------------------------------
!***  MXSNAL
!-----------------------------------------------------------------------
!
      DO J=JMS,JME
      DO I=IMS,IME
        DUMMY_2D_IN(I,J,1)=solver_int_state%MXSNAL(I,J)
      ENDDO
      ENDDO
!
      CALL PARENT_TO_CHILD_FILL(DUMMY_2D_IN                             &
                               ,NLEV                                    &
                               ,'MXSNAL'                                &
                               ,DUMMY_2D_OUT                            &
                               ,BILINEAR)
!
      IF(MYPE==0)WRITE(NFCST)DUMMY_2D_OUT
!     write(0,*)' after MXSNAL'
!
!-----------------------------------------------------------------------
!***  TSKIN  
!-----------------------------------------------------------------------
!
!     write(0,*)' before TSKIN'
      DO J=JMS,JME
      DO I=IMS,IME
        DUMMY_2D_IN(I,J,1)=solver_int_state%TSKIN(I,J)
      ENDDO
      ENDDO
!
      CALL PARENT_TO_CHILD_FILL(DUMMY_2D_IN                             &
                               ,NLEV                                    &
                               ,'TSKIN'                                 &
                               ,DUMMY_2D_OUT                            &
                               ,BILINEAR)
!
      IF(MYPE==0)THEN
        DO J=1,JM_CHILD
        DO I=1,IM_CHILD
          IF(DUMMY_2D_OUT(I,J,1)<150.)THEN
            SEA_MASK(I,J)=1.0
            DUMMY_2D_OUT(I,J,1)=0.
          ENDIF
          if(dummy_2d_out(i,j,1)<173..and.sea_mask(i,j)<0.5)then
            write(0,*)' Very cold TSKIN=',dummy_2d_out(i,j,1) &
                     ,' at (',i,',',j,')'
          endif
        ENDDO
        ENDDO
      ENDIF
!
      IF(MYPE==0)WRITE(NFCST)DUMMY_2D_OUT
!     write(0,*)' after TSKIN'
!
!-----------------------------------------------------------------------
!***  SST
!-----------------------------------------------------------------------
!
      DO J=JMS,JME
      DO I=IMS,IME
        DUMMY_2D_IN(I,J,1)=solver_int_state%SST(I,J)
      ENDDO
      ENDDO
!
      CALL PARENT_TO_CHILD_FILL(DUMMY_2D_IN                             &
                               ,NLEV                                    &
                               ,'SST'                                   &
                               ,DUMMY_2D_OUT                            &
                               ,BILINEAR)
!
      IF(MYPE==0)WRITE(NFCST)DUMMY_2D_OUT
!     write(0,*)' after SST'
!
!-----------------------------------------------------------------------
!***  SNO
!-----------------------------------------------------------------------
!
      DO J=JMS,JME
      DO I=IMS,IME
        DUMMY_2D_IN(I,J,1)=solver_int_state%SNO(I,J)
      ENDDO
      ENDDO
!
      CALL PARENT_TO_CHILD_FILL(DUMMY_2D_IN                             &
                               ,NLEV                                    &
                               ,'SNO'                                   &
                               ,DUMMY_2D_OUT                            &
                               ,BILINEAR)
!
      IF(MYPE==0)WRITE(NFCST)DUMMY_2D_OUT
!     write(0,*)' after SNO'
!
!-----------------------------------------------------------------------
!***  SI
!-----------------------------------------------------------------------
!
      DO J=JMS,JME
      DO I=IMS,IME
        DUMMY_2D_IN(I,J,1)=solver_int_state%SI(I,J)
      ENDDO
      ENDDO
!
      CALL PARENT_TO_CHILD_FILL(DUMMY_2D_IN                             &
                               ,NLEV                                    &
                               ,'SI'                                    &
                               ,DUMMY_2D_OUT                            &
                               ,BILINEAR)
!
      IF(MYPE==0)WRITE(NFCST)DUMMY_2D_OUT
!     write(0,*)' after SI'
!
!-----------------------------------------------------------------------
!***  SICE
!-----------------------------------------------------------------------
!
      DO J=JMS,JME
      DO I=IMS,IME
        DUMMY_2D_IN(I,J,1)=solver_int_state%SICE(I,J)
      ENDDO
      ENDDO
!
!     write(0,*)' PARENT_TO_CHILD_INIT SICE max=',maxval(DUMMY_2D_IN(IMS:IME,JMS:JME,1)) &
!              ,' min=',minval(DUMMY_2D_IN(IMS:IME,JMS:JME,1))
      CALL PARENT_TO_CHILD_FILL(DUMMY_2D_IN                             &
                               ,NLEV                                    &
                               ,'SICE'                                  &
                               ,DUMMY_2D_OUT                            &
                               ,BILINEAR)
!
      IF(MYPE==0)THEN
        DO J=1,JM_CHILD
        DO I=1,IM_CHILD
          SEA_ICE(I,J)=DUMMY_2D_OUT(I,J,1)
        ENDDO
        ENDDO
      ENDIF
!
      IF(MYPE==0)WRITE(NFCST)DUMMY_2D_OUT
!     write(0,*)' after SICE'
!
!-----------------------------------------------------------------------
!***  TG  
!-----------------------------------------------------------------------
!
      DO J=JMS,JME
      DO I=IMS,IME
        DUMMY_2D_IN(I,J,1)=solver_int_state%TG(I,J)
      ENDDO
      ENDDO
!
      CALL PARENT_TO_CHILD_FILL(DUMMY_2D_IN                             &
                               ,NLEV                                    &
                               ,'TG'                                    &
                               ,DUMMY_2D_OUT                            &
                               ,BILINEAR)
!
      IF(MYPE==0)WRITE(NFCST)DUMMY_2D_OUT
!     write(0,*)' after TG'
!
!-----------------------------------------------------------------------
!***  CMC
!-----------------------------------------------------------------------
!
      DO J=JMS,JME
      DO I=IMS,IME
        DUMMY_2D_IN(I,J,1)=solver_int_state%CMC(I,J)
      ENDDO
      ENDDO
!
      CALL PARENT_TO_CHILD_FILL(DUMMY_2D_IN                             &
                               ,NLEV                                    &
                               ,'CMC'                                   &
                               ,DUMMY_2D_OUT                            &
                               ,BILINEAR)
!
      IF(MYPE==0)WRITE(NFCST)DUMMY_2D_OUT
!     write(0,*)' after CMC'
!
!-----------------------------------------------------------------------
!***  SR
!-----------------------------------------------------------------------
!
      DO J=JMS,JME
      DO I=IMS,IME
        DUMMY_2D_IN(I,J,1)=solver_int_state%SR(I,J)
      ENDDO
      ENDDO
!
      CALL PARENT_TO_CHILD_FILL(DUMMY_2D_IN                             &
                               ,NLEV                                    &
                               ,'SR'                                    &
                               ,DUMMY_2D_OUT                            &
                               ,BILINEAR)
!
      IF(MYPE==0)WRITE(NFCST)DUMMY_2D_OUT
!     write(0,*)' after SR'
!
!-----------------------------------------------------------------------
!***  USTAR
!-----------------------------------------------------------------------
!
      DO J=JMS,JME
      DO I=IMS,IME
        DUMMY_2D_IN(I,J,1)=solver_int_state%USTAR(I,J)
      ENDDO
      ENDDO
!
      CALL PARENT_TO_CHILD_FILL(DUMMY_2D_IN                             &
                               ,NLEV                                    &
                               ,'USTAR'                                 &
                               ,DUMMY_2D_OUT                            &
                               ,BILINEAR)
!
      IF(MYPE==0)WRITE(NFCST)DUMMY_2D_OUT
!     write(0,*)' after USTAR'
!
!-----------------------------------------------------------------------
!***  Z0
!-----------------------------------------------------------------------
!
      DO J=JMS,JME
      DO I=IMS,IME
        DUMMY_2D_IN(I,J,1)=solver_int_state%Z0(I,J)
      ENDDO
      ENDDO
!
      CALL PARENT_TO_CHILD_FILL(DUMMY_2D_IN                             &
                               ,NLEV                                    &
                               ,'Z0'                                    &
                               ,DUMMY_2D_OUT                            &
                               ,BILINEAR)
!
      IF(MYPE==0)WRITE(NFCST)DUMMY_2D_OUT
!     write(0,*)' after Z0'
!
!-----------------------------------------------------------------------
!***  Z0BASE
!-----------------------------------------------------------------------
!
      DO J=JMS,JME
      DO I=IMS,IME
        DUMMY_2D_IN(I,J,1)=solver_int_state%Z0BASE(I,J)
      ENDDO
      ENDDO
!
      CALL PARENT_TO_CHILD_FILL(DUMMY_2D_IN                             &
                               ,NLEV                                    &
                               ,'Z0BASE'                                &
                               ,DUMMY_2D_OUT                            &
                               ,BILINEAR)
!
      IF(MYPE==0)WRITE(NFCST)DUMMY_2D_OUT
!     write(0,*)' after Z0BASE'
!
!-----------------------------------------------------------------------
!***  STC
!-----------------------------------------------------------------------
!
      CALL PARENT_TO_CHILD_FILL(solver_int_state%STC, NUM_SOIL_LAYERS      &
                               ,'STC'                                   &
                               ,DUMMY_3DS                               &
                               ,BILINEAR)
!
      IF(MYPE==0)THEN
        DO L=1,NUM_SOIL_LAYERS
        DO J=1,JM_CHILD
        DO I=1,IM_CHILD
          TEMPSOIL(L,I,J)=DUMMY_3DS(I,J,L)
        ENDDO
        ENDDO
        ENDDO
!
        WRITE(NFCST)TEMPSOIL
      ENDIF
!     write(0,*)' after STC'
!
!-----------------------------------------------------------------------
!***  SMC
!-----------------------------------------------------------------------
!
      CALL PARENT_TO_CHILD_FILL(solver_int_state%SMC, NUM_SOIL_LAYERS      &
                               ,'SMC'                                   &
                               ,DUMMY_3DS                               &
                               ,BILINEAR)
!
      IF(MYPE==0)THEN
        DO L=1,NUM_SOIL_LAYERS
        DO J=1,JM_CHILD
        DO I=1,IM_CHILD
          TEMPSOIL(L,I,J)=DUMMY_3DS(I,J,L)
        ENDDO
        ENDDO
        ENDDO
!
        WRITE(NFCST)TEMPSOIL
      ENDIF
!     write(0,*)' after SMC'
!
!-----------------------------------------------------------------------
!***  SH2O
!-----------------------------------------------------------------------
!
      CALL PARENT_TO_CHILD_FILL(solver_int_state%SH2O, NUM_SOIL_LAYERS     &
                               ,'SH2O'                                  &
                               ,DUMMY_3DS                               &
                               ,BILINEAR)
!
      IF(MYPE==0)THEN
        DO L=1,NUM_SOIL_LAYERS
        DO J=1,JM_CHILD
        DO I=1,IM_CHILD
          TEMPSOIL(L,I,J)=DUMMY_3DS(I,J,L)
        ENDDO
        ENDDO
        ENDDO
!
        WRITE(NFCST)TEMPSOIL
      ENDIF
!     write(0,*)' after SH2O'
!
!-----------------------------------------------------------------------
!***  ISLTYP
!-----------------------------------------------------------------------
!
      CALL PARENT_TO_CHILD_IFILL(solver_int_state%ISLTYP                   &
                                ,'ISLTYP'                               &
                                ,IDUMMY_2D )
!
      IF(MYPE==0)THEN
        DO J=1,JM_CHILD
        DO I=1,IM_CHILD
          IF(IDUMMY_2D(I,J)<1.AND.SEA_MASK(I,J)<0.5)THEN
            IDUMMY_2D(I,J)=1     !<--------- Bandaid for interpolated soil value=0 while interpolated seamask=0 (i.e, a land point)
!           if(abs(IDUMMY_2D(I,J))>50)write(0,*)' write ISLTYP i=',i,' j=',j,' sea_mask=',SEA_MASK(I,J),' isltyp=',IDUMMY_2D(I,J)
          ENDIF
        ENDDO
        ENDDO
!
        WRITE(NFCST)IDUMMY_2D
      ENDIF
!     write(0,*)' after ISLTYP'
!
!-----------------------------------------------------------------------
!***  IVGTYP
!-----------------------------------------------------------------------
!
!     write(0,*)' PARENT_TO_CHILD_INIT IVGTYP max=',maxval(solver_int_state%IVGTYP) &
!              ,' min=',minval(solver_int_state%IVGTYP),' maxloc=',maxloc(solver_int_state%IVGTYP) &
!              ,' minloc=',minloc(solver_int_state%IVGTYP)
      CALL PARENT_TO_CHILD_IFILL(solver_int_state%IVGTYP                   &
                                ,'IVGTYP'                               &
                                ,IDUMMY_2D )
!
      IF(MYPE==0)THEN
        DO J=1,JM_CHILD
        DO I=1,IM_CHILD
          IF(IDUMMY_2D(I,J)<1.AND.SEA_MASK(I,J)<0.5)THEN
            IDUMMY_2D(I,J)=1     !<--------- Bandaid for interpolated vegetation value=0 while interpolated seamask=0 (i.e, a land point)
          ENDIF
        ENDDO
        ENDDO
!
        WRITE(NFCST)IDUMMY_2D
      ENDIF
!     write(0,*)' after IVGTYP'
!
!-----------------------------------------------------------------------
!***  VEGFRC
!-----------------------------------------------------------------------
!
      DO J=JMS,JME
      DO I=IMS,IME
        DUMMY_2D_IN(I,J,1)=solver_int_state%VEGFRC(I,J)
      ENDDO
      ENDDO
!
      CALL PARENT_TO_CHILD_FILL(DUMMY_2D_IN                             &
                               ,NLEV                                    &
                               ,'VEGFRC'                                &
                               ,DUMMY_2D_OUT                            &
                               ,BILINEAR)
!
      IF(MYPE==0)THEN
        DO J=1,JM_CHILD
        DO I=1,IM_CHILD
          IF(DUMMY_2D_OUT(I,J,1)>0..AND.(SEA_MASK(I,J)>0.5.OR.SEA_ICE(I,J)>0.))THEN
            DUMMY_2D_OUT(I,J,1)=0.    !<--------- Bandaid for interpolated veg frac value >0 while interpolated seamask or sice >0
          ENDIF
        ENDDO
        ENDDO
      ENDIF
!
      IF(MYPE==0)WRITE(NFCST)DUMMY_2D_OUT
!     write(0,*)' after VEGFRC'
!
!-----------------------------------------------------------------------
!
      IF(MYPE==0)WRITE(NFCST)DUMMY_SOIL
      IF(MYPE==0)WRITE(NFCST)DUMMY_SOIL
!
!-----------------------------------------------------------------------
!
      IF(MYPE==0)CLOSE(NFCST)
!
!-----------------------------------------------------------------------
!
      DEALLOCATE(IDUMMY_2D)
      DEALLOCATE(DUMMY_2D_IN)
      DEALLOCATE(DUMMY_2D_OUT)
      DEALLOCATE(DUMMY_3D)
      DEALLOCATE(DUMMY_3DS)
      DEALLOCATE(DUMMY_SOIL)
      DEALLOCATE(TEMPSOIL)
      DEALLOCATE(SEA_MASK)
      DEALLOCATE(PD_BILINEAR)
      DEALLOCATE(PD_NEAREST)
      DEALLOCATE(LOWER_TOPO)
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!!!   SUBROUTINE PARENT_TO_CHILD_FILL_ASSOC(PARENT_ARRAY                &
      SUBROUTINE PARENT_TO_CHILD_FILL      (PARENT_ARRAY                &
                                           ,NLEV                        &
                                           ,VBL_NAME                    &
                                           ,CHILD_ARRAY                 &
                                           ,METHOD)
!
!-----------------------------------------------------------------------
!***  Rows and columns of the child's grid lie directly on top of
!***  rows and colums of the parent (thus 'ASSOCIATED').
!
!***  Fill a child's domain with data from the parent.  Only the parent
!***  tasks are needed in this routine.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      INTEGER(kind=KINT),INTENT(IN) :: METHOD                            &  !<-- Interpolaton method (bilinear or nearest neighbor)
                                      ,NLEV                                 !<-- Vertical dimension of the data array 
!
!!!   REAL(kind=KFPT),DIMENSION(IMS:IME,JMS:JME,1:NLEV),INTENT(INOUT) :: &
!!!                                                           DATA_ARRAY    !<-- The parent array that will initialize the child array
!
      REAL(kind=KFPT),DIMENSION(IMS:IME,JMS:JME,1:NLEV),INTENT(IN) ::    &
                                                              PARENT_ARRAY  !<-- The parent array that will initialize the child array
!
      CHARACTER(*),INTENT(IN) :: VBL_NAME
!
      REAL(kind=KFPT),DIMENSION(1:IM_CHILD,1:JM_CHILD,1:NLEV),INTENT(OUT) :: &  !<-- Data from parent tasks interpolated to child grid
                                                              CHILD_ARRAY       !    but still on parent task 0
!
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: I,IERR,II,IPE,IPE_LOCAL,ISTAT,J,JJ,L,N,NN
      INTEGER(kind=KINT) :: I_COPY,I_END,I_END_COPY,I_EXTENT,I_PARENT_END  &
                ,I_START_COPY
      INTEGER(kind=KINT) :: J_COPY,J_END,J_END_COPY,J_EXTENT,J_PARENT_END  &
                ,J_START_COPY
      INTEGER(kind=KINT) :: INDX_EAST,INDX_NORTH,INDX_SOUTH,INDX_WEST
      INTEGER(kind=KINT) :: NWORDS_RECV,NWORDS_SEND
!
      INTEGER(kind=KINT),DIMENSION(MPI_STATUS_SIZE) :: JSTAT
!
      REAL(kind=KFPT) :: DELTA_I_EAST,DELTA_I_WEST                       &
                        ,DELTA_J_NORTH,DELTA_J_SOUTH
      REAL(kind=KFPT) :: RATIO,REAL_INDX_I_PARENT,REAL_INDX_J_PARENT
      REAL(kind=KFPT) :: WEIGHT_EAST,WEIGHT_NORTH                        &
                        ,WEIGHT_SOUTH,WEIGHT_WEST
      REAL(kind=KFPT) :: WEIGHT_NE,WEIGHT_NW,WEIGHT_SE,WEIGHT_SW
      REAL(kind=KFPT) :: WEIGHT_MAX,WEIGHT_SUM
!
      REAL(kind=KFPT),DIMENSION(:)    ,ALLOCATABLE :: DATA_BUFFER
      REAL(kind=KFPT),DIMENSION(:,:,:),ALLOCATABLE :: ARRAY_STAGE_PARENT    
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  To simplify matters somewhat, isolate the minimum subset of
!***  points on the parent domain that underlie the child's domain.
!
!***  The southwest corner of the child always lies directly on a
!***  point in the parent domain.  We already know the I,J of that
!***  parent point since it was specified in the configure file.
!***  The number of parent points that are covered by the child is
!***  determined by the child-to-parent grid ratio and the lateral
!***  dimensions of the child's domain.
!-----------------------------------------------------------------------
!
      I_PARENT_END=I_PARENT_START                                       &  !<-- Easternmost I on parent domain surrounding child domain
                     +INT((IM_CHILD-1)*CHILD_PARENT_SPACE_RATIO)+1
!
      I_EXTENT=I_PARENT_END-I_PARENT_START+1
!
      J_PARENT_END=J_PARENT_START                                       &  !<-- Northernmost J on parent domain surrounding child domain
                     +INT((JM_CHILD-1)*CHILD_PARENT_SPACE_RATIO)+1
!
      J_EXTENT=J_PARENT_END-J_PARENT_START+1
!
!-----------------------------------------------------------------------
!***  Create a staging array on parent task 0 that will hold the entire
!***  subset of the parent domain underlying the child.
!***  Then all parent tasks with points in the intersecting region
!***  send their data to parent task 0.
!-----------------------------------------------------------------------
!
      parent_stage: IF(MYPE==0)THEN                                        !<-- Parent task 0
! 
!-----------------------------------------------------------------------
!
        ALLOCATE(ARRAY_STAGE_PARENT(1:I_EXTENT,1:J_EXTENT,1:NLEV))         !<-- Array holding all parent points in staging region
                                                                           !    Note that this array begins at (1,1,1), i.e.,
                                                                           !      its indices are relative to the nest.
!
!-----------------------------------------------------------------------
!***  If parent task 0 holds some of the staging region, copy it to
!***  the staging array.
!-----------------------------------------------------------------------
!
        IF(I_PARENT_START<=ITE.AND.J_PARENT_START<=JTE)THEN
          I_END=MIN(ITE,I_PARENT_END)
          J_END=MIN(JTE,J_PARENT_END)
!
          DO L=1,NLEV
            JJ=0
            DO J=J_PARENT_START,J_END
              JJ=JJ+1
!
              II=0
              DO I=I_PARENT_START,I_END
                II=II+1     
                ARRAY_STAGE_PARENT(II,JJ,L)=PARENT_ARRAY(I,J,L)
              ENDDO
            ENDDO
          ENDDO
!
        ENDIF
!
!-----------------------------------------------------------------------
!***  If there are points in the staging region outside of parent task 0
!***  then task 0 receives those points from the other parent tasks that
!***  contain those points.
!-----------------------------------------------------------------------
!
        parent_search: DO IPE=1,NUM_PES_PARENT-1                           !<-- Parent task 0 checks other parent fcst tasks for points
!
          remote_stage: IF(I_PARENT_START<=LOCAL_IEND  (IPE).AND.       &  !<-- Does remote parent task IPE contain any staging region?
                           I_PARENT_END  >=LOCAL_ISTART(IPE)            &  !
                            .AND.                                       &  !
                           J_PARENT_START<=LOCAL_JEND  (IPE).AND.       &  !
                           J_PARENT_END  >=LOCAL_JSTART(IPE))THEN          !<--
! 
            I_START_COPY=MAX(I_PARENT_START,LOCAL_ISTART(IPE))             !<-- I index of first point in staging region on remote parent task
            I_END_COPY  =MIN(I_PARENT_END  ,LOCAL_IEND  (IPE))             !<-- I index of last point in staging region on remote parent task
            I_COPY      =I_END_COPY-I_START_COPY+1                         !<-- I range of points to receive
!
            J_START_COPY=MAX(J_PARENT_START,LOCAL_JSTART(IPE))             !<-- J index of first point in staging region on remote parent task
            J_END_COPY  =MIN(J_PARENT_END  ,LOCAL_JEND  (IPE))             !<-- J index of last point in staging region on remote parent task
            J_COPY      =J_END_COPY-J_START_COPY+1                         !<-- J range of points to receive
!
            NWORDS_RECV=I_COPY*J_COPY*NLEV                                 !<-- Total # of words from remote parent task in staging region
!
            ALLOCATE(DATA_BUFFER(1:NWORDS_RECV))                           !<-- Allocate buffer array to hold remote task's staging data
            CALL MPI_RECV(DATA_BUFFER                                   &  !<-- The staging region data from remote parent task IPE
                         ,NWORDS_RECV                                   &  !<-- Total words received
                         ,MPI_REAL                                      &  !<-- Datatype
                         ,IPE                                           &  !<-- Receive from this parent task
                         ,IPE                                           &  !<-- MPI tag
                         ,COMM_MY_DOMAIN                                &  !<-- The MPI communicator
                         ,JSTAT                                         &  !<-- MPI status object
                         ,IERR )
!
            NN=0                                                           !<-- Counter for received words
!
            DO L=1,NLEV
!
              JJ=J_START_COPY-J_PARENT_START
              DO J=1,J_COPY
                JJ=JJ+1
!
                II=I_START_COPY-I_PARENT_START
                DO I=1,I_COPY
                  II=II+1
                  NN=NN+1
                  ARRAY_STAGE_PARENT(II,JJ,L)=DATA_BUFFER(NN)              !<-- Fill in array with staging region data from parent task IPE
                ENDDO
              ENDDO
            ENDDO
!
            DEALLOCATE(DATA_BUFFER)
!
          ENDIF remote_stage
!
        ENDDO parent_search
!
!-----------------------------------------------------------------------
!***  Now the remaining parent tasks check to see if they contain
!***  any points in the staging region.  If they do, gather them
!***  and send them to parent task 0.
!-----------------------------------------------------------------------
!
      ELSEIF(MYPE>0.AND.MYPE<=NUM_PES_PARENT-1)THEN  parent_stage          !<-- All parent forecast tasks other than 0
!
!-----------------------------------------------------------------------
        IF(I_PARENT_START<=ITE.AND.I_PARENT_END>=ITS                    &  !<-- Does this parent task contain any staging region?
           .AND.                                                        &  !
           J_PARENT_START<=JTE.AND.J_PARENT_END>=JTS)THEN                  !<--
! 
          I_START_COPY=MAX(I_PARENT_START,ITS)                             !<-- I index of first point in staging region on this parent task
          I_END_COPY  =MIN(I_PARENT_END  ,ITE)                             !<-- I index of last point in staging region on this parent task
          I_COPY=I_END_COPY-I_START_COPY+1                                 !<-- I range of points to send to parent task 0
!
          J_START_COPY=MAX(J_PARENT_START,JTS)                             !<-- J index of first point in staging region on remote parent task
          J_END_COPY  =MIN(J_PARENT_END  ,JTE)                             !<-- J index of last point in staging region on remote parent task
          J_COPY=J_END_COPY-J_START_COPY+1                                 !<-- J range of copied points
!
          NWORDS_SEND=I_COPY*J_COPY*NLEV                                   !<-- Total number of words from this parent task in staging region
          ALLOCATE(DATA_BUFFER(1:NWORDS_SEND),stat=ISTAT)                  !<-- Allocate the buffer array to hold this task's staging data
!
          NN=0
!
          DO L=1,NLEV
          DO J=J_START_COPY,J_END_COPY
          DO I=I_START_COPY,I_END_COPY
            NN=NN+1
            DATA_BUFFER(NN)=PARENT_ARRAY(I,J,L)
          ENDDO
          ENDDO
          ENDDO
!
          CALL MPI_SEND(DATA_BUFFER                                     &  !<-- The staging region data from this parent task to parent task 0
                       ,NWORDS_SEND                                     &  !<-- Total words sent
                       ,MPI_REAL                                        &  !<-- Datatype
                       ,0                                               &  !<-- Send to parent task 0
                       ,MYPE                                            &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- The MPI communicator
                       ,IERR )
!
          DEALLOCATE(DATA_BUFFER)
!
        ENDIF
!
!-----------------------------------------------------------------------
!
      ENDIF parent_stage
!
!-----------------------------------------------------------------------
!***  The subset of the input array on the parent domain that lies 
!***  under the child's domain has been mirrored onto parent task 0.
!***  Parent task 0 will fill out the array to match the child
!***  domain's horizontal grid increments and then parcel out the
!***  appropriate pieces to the corresponding tasks of the child.
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!***  First fill in the southern and western sides of the array.
!***  If bilinear interpolation is specified then only linear 
!***  interpolation needs to be used.
!-----------------------------------------------------------------------
!
      parent_task_0: IF(MYPE==0)THEN                                      !<-- Parent task 0
!
!-----------------------------------------------------------------------
!
        RATIO=CHILD_PARENT_SPACE_RATIO
!
        DO L=1,NLEV
!
          CHILD_ARRAY(1,1,L)=ARRAY_STAGE_PARENT(1,1,L)                    !<-- SW corner of child's array coincides with a parent point
!
          DO I=2,IM_CHILD                                                 !<-- Move along southern boundary of child's domain
            REAL_INDX_I_PARENT=1+(I-1)*RATIO                              !<-- Exact I index of child point on parent 
            INDX_WEST=INT(REAL_INDX_I_PARENT)                             !<-- The parent point's I index west of the child's point
            INDX_EAST=INDX_WEST+1                                         !<-- The parent point's I index east of the child's point
            WEIGHT_WEST=INDX_EAST-REAL_INDX_I_PARENT                      !<-- Interpolation weight given parent's point to the west
            WEIGHT_EAST=1.-WEIGHT_WEST                                    !<-- Interpolation weight given parent's point to the east
!
            IF(METHOD==NEAREST)THEN                                       !<-- Assign points using nearest neighbors
              WEIGHT_MAX=MAX(WEIGHT_WEST,WEIGHT_EAST)
              IF(WEIGHT_WEST==WEIGHT_MAX)THEN
                CHILD_ARRAY(I,1,L)=ARRAY_STAGE_PARENT(INDX_WEST,1,L)
              ELSEIF(WEIGHT_EAST==WEIGHT_MAX)THEN
                CHILD_ARRAY(I,1,L)=ARRAY_STAGE_PARENT(INDX_EAST,1,L)
              ENDIF
!
            ELSEIF(METHOD==BILINEAR)THEN                                         !<-- Assign points using (bi)linear interpolation
              CHILD_ARRAY(I,1,L)=WEIGHT_WEST*ARRAY_STAGE_PARENT(INDX_WEST,1,L) & !<-- Value at points along child's southern boundary
                                +WEIGHT_EAST*ARRAY_STAGE_PARENT(INDX_EAST,1,L)
            ELSE
              WRITE(0,*)" Attempting to use unknown interpolation method: ",METHOD
              CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
!
            ENDIF
!
          ENDDO
!
          DO J=2,JM_CHILD                                                 !<-- Move along western boundary of child's domain
            REAL_INDX_J_PARENT=1+(J-1)*RATIO                              !<-- Exact J index of child point on parent 
            INDX_SOUTH=INT(REAL_INDX_J_PARENT)                            !<-- The parent point's J index south of the child's point
            INDX_NORTH=INDX_SOUTH+1                                       !<-- The parent point's J index north of the child's point
            WEIGHT_SOUTH=INDX_NORTH-REAL_INDX_J_PARENT                    !<-- Interpolation weight of parent's point to the south
            WEIGHT_NORTH=1.-WEIGHT_SOUTH                                  !<-- Interpolation weight of parent's point to the north
!
            IF(METHOD==NEAREST)THEN                                       !<-- Assign points using nearest neighbors
              WEIGHT_MAX=MAX(WEIGHT_SOUTH,WEIGHT_NORTH)
              IF(WEIGHT_SOUTH==WEIGHT_MAX)THEN
                CHILD_ARRAY(1,J,L)=ARRAY_STAGE_PARENT(1,INDX_SOUTH,L)
              ELSE
                CHILD_ARRAY(1,J,L)=ARRAY_STAGE_PARENT(1,INDX_NORTH,L)
              ENDIF
!
            ELSEIF(METHOD==BILINEAR)THEN                                           !<-- Assign points using (bi)linear interpolation
              CHILD_ARRAY(1,J,L)=WEIGHT_SOUTH*ARRAY_STAGE_PARENT(1,INDX_SOUTH,L) & !<-- Value at points along child's western boundary
                              +WEIGHT_NORTH*ARRAY_STAGE_PARENT(1,INDX_NORTH,L)
            ELSE
              WRITE(0,*)" Attempting to use unknown interpolation method: ",METHOD
              CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
!
            ENDIF
!
          ENDDO
!
!-----------------------------------------------------------------------
!***  Fill in the interior of the staging array.
!-----------------------------------------------------------------------
!
          DO J=2,JM_CHILD
            REAL_INDX_J_PARENT=1+(J-1)*RATIO                              !<-- Exact J index of child point in parent staging region 
            INDX_SOUTH=INT(REAL_INDX_J_PARENT)                            !<-- The parent point's J index south of the child's point
            INDX_NORTH=INDX_SOUTH+1                                       !<-- The parent point's J index north of the child's point
!
            DELTA_J_NORTH=INDX_NORTH-REAL_INDX_J_PARENT                   !<-- Parent grid increment from child point to parent point north
            DELTA_J_SOUTH=REAL_INDX_J_PARENT-INDX_SOUTH                   !<-- Parent grid increment from child point to parent point south
!
            DO I=2,IM_CHILD
              REAL_INDX_I_PARENT=1+(I-1)*RATIO                            !<-- Exact I index of child point in parent staging region
              INDX_WEST=INT(REAL_INDX_I_PARENT)                           !<-- The parent point's I index west of the child's point
              INDX_EAST=INDX_WEST+1                                       !<-- The parent point's I index east of the child's point
!
              DELTA_I_EAST=INDX_EAST-REAL_INDX_I_PARENT
              DELTA_I_WEST=REAL_INDX_I_PARENT-INDX_WEST
!
              WEIGHT_SW=DELTA_I_EAST*DELTA_J_NORTH                        !<-- Interpolation weight of parent's point to SW 
              WEIGHT_SE=DELTA_I_WEST*DELTA_J_NORTH                        !<-- Interpolation weight of parent's point to SE
              WEIGHT_NW=DELTA_I_EAST*DELTA_J_SOUTH                        !<-- Interpolation weight of parent's point to NW
              WEIGHT_NE=DELTA_I_WEST*DELTA_J_SOUTH                        !<-- Interpolation weight of parent's point to NE
!
!-----------------------------------------------------------------------
!
              assign: IF(METHOD==NEAREST)THEN                             !<-- Assign points using nearest neighbors      
                WEIGHT_MAX=MAX(WEIGHT_SW,WEIGHT_SE                     &
                              ,WEIGHT_NW,WEIGHT_NE)
                IF(WEIGHT_SW==WEIGHT_MAX)THEN
                  CHILD_ARRAY(I,J,L)=ARRAY_STAGE_PARENT(INDX_WEST,INDX_SOUTH,L)
                ELSEIF(WEIGHT_SE==WEIGHT_MAX)THEN
                  CHILD_ARRAY(I,J,L)=ARRAY_STAGE_PARENT(INDX_EAST,INDX_SOUTH,L)
                ELSEIF(WEIGHT_NW==WEIGHT_MAX)THEN
                  CHILD_ARRAY(I,J,L)=ARRAY_STAGE_PARENT(INDX_WEST,INDX_NORTH,L)
                ELSEIF(WEIGHT_NE==WEIGHT_MAX)THEN
                  CHILD_ARRAY(I,J,L)=ARRAY_STAGE_PARENT(INDX_EAST,INDX_NORTH,L)
                ENDIF
!
              ELSEIF(METHOD==BILINEAR)THEN                                                 !<-- Assign points using bilinear interpolation
                IF(VBL_NAME/='FIS')THEN
                  IF(ABS(ARRAY_STAGE_PARENT(INDX_WEST,INDX_SOUTH,L))<1.E-12)WEIGHT_SW=0.
                  IF(ABS(ARRAY_STAGE_PARENT(INDX_EAST,INDX_SOUTH,L))<1.E-12)WEIGHT_SE=0.
                  IF(ABS(ARRAY_STAGE_PARENT(INDX_WEST,INDX_NORTH,L))<1.E-12)WEIGHT_NW=0.
                  IF(ABS(ARRAY_STAGE_PARENT(INDX_EAST,INDX_NORTH,L))<1.E-12)WEIGHT_NE=0.
                ENDIF
!
                CHILD_ARRAY(I,J,L)=WEIGHT_SW*ARRAY_STAGE_PARENT(INDX_WEST,INDX_SOUTH,L)  & !<-- Value at points in child's interior
                                  +WEIGHT_SE*ARRAY_STAGE_PARENT(INDX_EAST,INDX_SOUTH,L)  &
                                  +WEIGHT_NW*ARRAY_STAGE_PARENT(INDX_WEST,INDX_NORTH,L)  &
                                  +WEIGHT_NE*ARRAY_STAGE_PARENT(INDX_EAST,INDX_NORTH,L)
!
                WEIGHT_SUM=WEIGHT_SW+WEIGHT_SE+WEIGHT_NW+WEIGHT_NE
                IF(WEIGHT_SUM<0.99.AND.WEIGHT_SUM>0.01)THEN
                  CHILD_ARRAY(I,J,L)=CHILD_ARRAY(I,J,L)/WEIGHT_SUM        !<-- Normalize if some weights are zero (e.g., coastal land Temp)
                ENDIF
!
                IF(VBL_NAME=='SST')THEN                                   !<-- Include only realistic SST temperatures
                  WEIGHT_SUM=0.
                  CHILD_ARRAY(I,J,L)=0.
                  IF(ARRAY_STAGE_PARENT(INDX_WEST,INDX_SOUTH,L)>200.)THEN
                    WEIGHT_SUM=WEIGHT_SUM+WEIGHT_SW
                    CHILD_ARRAY(I,J,L)=WEIGHT_SW*ARRAY_STAGE_PARENT(INDX_WEST,INDX_SOUTH,L) &
                                       +CHILD_ARRAY(I,J,L)
                  ENDIF
                  IF(ARRAY_STAGE_PARENT(INDX_EAST,INDX_SOUTH,L)>200.)THEN
                    WEIGHT_SUM=WEIGHT_SUM+WEIGHT_SE
                    CHILD_ARRAY(I,J,L)=WEIGHT_SE*ARRAY_STAGE_PARENT(INDX_EAST,INDX_SOUTH,L) &
                                       +CHILD_ARRAY(I,J,L)
                  ENDIF
                  IF(ARRAY_STAGE_PARENT(INDX_WEST,INDX_NORTH,L)>200.)THEN
                    WEIGHT_SUM=WEIGHT_SUM+WEIGHT_NW
                    CHILD_ARRAY(I,J,L)=WEIGHT_NW*ARRAY_STAGE_PARENT(INDX_WEST,INDX_NORTH,L) &
                                       +CHILD_ARRAY(I,J,L)
                  ENDIF
                  IF(ARRAY_STAGE_PARENT(INDX_EAST,INDX_NORTH,L)>200.)THEN
                    WEIGHT_SUM=WEIGHT_SUM+WEIGHT_NE
                    CHILD_ARRAY(I,J,L)=WEIGHT_NE*ARRAY_STAGE_PARENT(INDX_EAST,INDX_NORTH,L) &
                                       +CHILD_ARRAY(I,J,L)
                  ENDIF
                  IF(WEIGHT_SUM<0.99.AND.WEIGHT_SUM>0.01)THEN
                    CHILD_ARRAY(I,J,L)=CHILD_ARRAY(I,J,L)/WEIGHT_SUM
                  ENDIF
                ENDIF
!
              ELSE
                WRITE(0,*)" Attempting to use unknown interpolation method: ",METHOD
                CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
!
              ENDIF assign
!
!-----------------------------------------------------------------------
!
            ENDDO
          ENDDO
!
        ENDDO
!
        DEALLOCATE(ARRAY_STAGE_PARENT)
!
!-----------------------------------------------------------------------
!
      ENDIF parent_task_0
!
!-----------------------------------------------------------------------
!
!!!   END SUBROUTINE PARENT_TO_CHILD_FILL_ASSOC
      END SUBROUTINE PARENT_TO_CHILD_FILL
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!!!   SUBROUTINE PARENT_TO_CHILD_IFILL_ASSOC(PARENT_ARRAY               &
      SUBROUTINE PARENT_TO_CHILD_IFILL      (PARENT_ARRAY               &
                                            ,VBL_NAME                   &
                                            ,CHILD_ARRAY)
!
!-----------------------------------------------------------------------
!***  Rows and columns of the child's grid lie directly on top of
!***  rows and colums of the parent (thus 'ASSOCIATED').
!***  Fill a child's domain with data from the parent.  Only the parent
!***  tasks are needed in this routine.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      INTEGER(kind=KINT),DIMENSION(IMS:IME,JMS:JME),INTENT(IN) ::       &
                                                          PARENT_ARRAY      !<-- The parent array that will initialize the child array
!
      INTEGER(kind=KINT),DIMENSION(1:IM_CHILD,1:JM_CHILD),INTENT(OUT) :: &  !<-- Data from parent tasks interpolated to child grid
                                                           CHILD_ARRAY      !      but still on parent task 0
!
      CHARACTER(*),INTENT(IN) :: VBL_NAME                                   !<-- The variable's name
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: I,IERR,II,IPE,IPE_LOCAL,J,JJ,N,NN
      INTEGER(kind=KINT) :: I_COPY,I_END,I_END_COPY,I_EXTENT,I_PARENT_END  &
                           ,I_START_COPY
      INTEGER(kind=KINT) :: J_COPY,J_END,J_END_COPY,J_EXTENT,J_PARENT_END  &
                           ,J_START_COPY
      INTEGER(kind=KINT) :: INDX_EAST,INDX_NORTH,INDX_SOUTH,INDX_WEST
      INTEGER(kind=KINT) :: NWORDS_RECV,NWORDS_SEND
!
      INTEGER(kind=KINT),DIMENSION(MPI_STATUS_SIZE) :: JSTAT
!
      INTEGER,DIMENSION(:)  ,ALLOCATABLE :: DATA_BUFFER
      INTEGER,DIMENSION(:,:),ALLOCATABLE :: ARRAY_STAGE_PARENT
!
      REAL(kind=KFPT) :: DELTA_I_EAST,DELTA_I_WEST,DELTA_J_NORTH,DELTA_J_SOUTH
      REAL(kind=KFPT) :: RATIO,REAL_INDX_I_PARENT,REAL_INDX_J_PARENT
      REAL(kind=KFPT) :: WEIGHT_EAST,WEIGHT_NORTH,WEIGHT_SOUTH,WEIGHT_WEST
      REAL(kind=KFPT) :: WEIGHT_NE,WEIGHT_NW,WEIGHT_SE,WEIGHT_SW
      REAL(kind=KFPT) :: WEIGHT_MAX
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  To simplify matters somewhat, isolate the minimum subset of
!***  points on the parent domain that underlie the child's domain.
!
!***  The southwest corner of the child always lies directly on a
!***  point in the parent domain.  We already know the I,J of that
!***  parent point since it was specified in the configure file.
!***  The number of parent points that are covered by the child is
!***  determined by the child-to-parent grid ratio and the lateral
!***  dimensions of the child's domain.
!-----------------------------------------------------------------------
!
      I_PARENT_END=I_PARENT_START                                       &  !<-- Easternmost I on parent domain surrounding child domain
                     +INT((IM_CHILD-1)*CHILD_PARENT_SPACE_RATIO)+1
!
      I_EXTENT=I_PARENT_END-I_PARENT_START+1
!
      J_PARENT_END=J_PARENT_START                                       &  !<-- Northernmost J on parent domain surrounding child domain
                     +INT((JM_CHILD-1)*CHILD_PARENT_SPACE_RATIO)+1
!
      J_EXTENT=J_PARENT_END-J_PARENT_START+1
!
!-----------------------------------------------------------------------
!***  Create a staging array on parent task 0 that will hold the entire
!***  subset of the parent domain underlying the child.
!***  Then all parent tasks with points in the intersecting region
!***  send their data to parent task 0.
!-----------------------------------------------------------------------
!
      parent_stage: IF(MYPE==0)THEN                                        !<-- Parent task 0
! 
!-----------------------------------------------------------------------
!
        ALLOCATE(ARRAY_STAGE_PARENT(1:I_EXTENT,1:J_EXTENT))                !<-- Array holding all parent points in staging region
                                                                           !    Note that this array begins at (1,1,1), i.e.,
                                                                           !      its indices are relative to the nest.
!
!-----------------------------------------------------------------------
!***  If parent task 0 holds some of the staging region, copy it to
!***  the staging array.
!-----------------------------------------------------------------------
!
        IF(I_PARENT_START<=ITE.AND.J_PARENT_START<=JTE)THEN
          I_END=MIN(ITE,I_PARENT_END)
          J_END=MIN(JTE,J_PARENT_END)
!
          JJ=0
          DO J=J_PARENT_START,J_END
            JJ=JJ+1
!
            II=0
            DO I=I_PARENT_START,I_END
              II=II+1     
              ARRAY_STAGE_PARENT(II,JJ)=PARENT_ARRAY(I,J)
            ENDDO
          ENDDO
!
        ENDIF
!
!-----------------------------------------------------------------------
!***  If there are points in the staging region outside of parent task 0
!***  then task 0 receives those points from the other parent tasks that
!***  contain those points.
!-----------------------------------------------------------------------
!
        parent_search: DO IPE=1,NUM_PES_PARENT-1                           !<-- Parent task 0 checks other parent tasks for points
!
          remote_stage: IF(I_PARENT_START<=LOCAL_IEND  (IPE).AND.       &  !<-- Does remote parent task IPE contain any staging region?
                           I_PARENT_END  >=LOCAL_ISTART(IPE)            &  !
                            .AND.                                       &  !
                           J_PARENT_START<=LOCAL_JEND  (IPE).AND.       &  !
                           J_PARENT_END  >=LOCAL_JSTART(IPE))THEN          !<--
! 
            I_START_COPY=MAX(I_PARENT_START,LOCAL_ISTART(IPE))             !<-- I index of first point in staging region on remote parent task
            I_END_COPY  =MIN(I_PARENT_END  ,LOCAL_IEND  (IPE))             !<-- I index of last point in staging region on remote parent task
            I_COPY      =I_END_COPY-I_START_COPY+1                         !<-- I range of points to receive
!
            J_START_COPY=MAX(J_PARENT_START,LOCAL_JSTART(IPE))             !<-- J index of first point in staging region on remote parent task
            J_END_COPY  =MIN(J_PARENT_END  ,LOCAL_JEND  (IPE))             !<-- J index of last point in staging region on remote parent task
            J_COPY      =J_END_COPY-J_START_COPY+1                         !<-- J range of points to receive
!
            NWORDS_RECV=I_COPY*J_COPY                                      !<-- Total # of words from remote parent task in staging region
!
            ALLOCATE(DATA_BUFFER(1:NWORDS_RECV))                           !<-- Allocate buffer array to hold remote task's staging data
            CALL MPI_RECV(DATA_BUFFER                                   &  !<-- The staging region data from remote parent task IPE
                         ,NWORDS_RECV                                   &  !<-- Total words received
                         ,MPI_INTEGER                                   &  !<-- Datatype
                         ,IPE                                           &  !<-- Receive from this parent task
                         ,IPE                                           &  !<-- MPI tag
                         ,COMM_MY_DOMAIN                                &  !<-- The MPI communicator
                         ,JSTAT                                         &  !<-- MPI status object
                         ,IERR )
!
            NN=0                                                           !<-- Counter for received words
!
              JJ=J_START_COPY-J_PARENT_START
              DO J=1,J_COPY
                JJ=JJ+1
!
                II=I_START_COPY-I_PARENT_START
                DO I=1,I_COPY
                  II=II+1
                  NN=NN+1
                  ARRAY_STAGE_PARENT(II,JJ)=DATA_BUFFER(NN)                !<-- Fill in array with staging region data from parent task IPE
                ENDDO
              ENDDO
!
            DEALLOCATE(DATA_BUFFER)
!
          ENDIF remote_stage
!
        ENDDO parent_search
!
!-----------------------------------------------------------------------
!***  Now the remaining parent tasks check to see if they contain
!***  any points in the staging region.  If they do, gather them
!***  and send them to parent task 0.
!-----------------------------------------------------------------------
!
      ELSEIF(MYPE>0.AND.MYPE<=NUM_PES_PARENT-1)THEN  parent_stage          !<-- All parent tasks other than 0
!
!-----------------------------------------------------------------------
        IF(I_PARENT_START<=ITE.AND.I_PARENT_END>=ITS                    &  !<-- Does this parent task contain any staging region?
           .AND.                                                        &  !
           J_PARENT_START<=JTE.AND.J_PARENT_END>=JTS)THEN                  !<--
! 
          I_START_COPY=MAX(I_PARENT_START,ITS)                             !<-- I index of first point in staging region on this parent task
          I_END_COPY  =MIN(I_PARENT_END  ,ITE)                             !<-- I index of last point in staging region on this parent task
          I_COPY=I_END_COPY-I_START_COPY+1                                 !<-- I range of points to send to parent task 0
!
          J_START_COPY=MAX(J_PARENT_START,JTS)                             !<-- J index of first point in staging region on remote parent task
          J_END_COPY  =MIN(J_PARENT_END  ,JTE)                             !<-- J index of last point in staging region on remote parent task
          J_COPY=J_END_COPY-J_START_COPY+1                                 !<-- J range of copied points
!
          NWORDS_SEND=I_COPY*J_COPY                                        !<-- Total number of words from this parent task in staging region
          ALLOCATE(DATA_BUFFER(1:NWORDS_SEND))                             !<-- Allocate the buffer array to hold this task's staging data
!
          NN=0
!
          DO J=J_START_COPY,J_END_COPY
          DO I=I_START_COPY,I_END_COPY
            NN=NN+1
            DATA_BUFFER(NN)=PARENT_ARRAY(I,J)
          ENDDO
          ENDDO
!
          CALL MPI_SEND(DATA_BUFFER                                     &  !<-- The staging region data from this parent task to parent task 0
                       ,NWORDS_SEND                                     &  !<-- Total words sent
                       ,MPI_INTEGER                                     &  !<-- Datatype
                       ,0                                               &  !<-- Send to parent task 0
                       ,MYPE                                            &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- The MPI communicator
                       ,IERR )
!
          DEALLOCATE(DATA_BUFFER)
!
        ENDIF
!
!-----------------------------------------------------------------------
!
      ENDIF parent_stage
!
!-----------------------------------------------------------------------
!***  The subset of the input array on the parent domain that lies 
!***  under the child's domain has been mirrored onto parent task 0.
!***  Parent task 0 will fill out the array to match the child
!***  domain's horizontal grid increments and then parcel out the
!***  appropriate pieces to the corresponding tasks of the child.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  First fill in the southern and western sides of the array
!***  choosing the nearest parent points.
!-----------------------------------------------------------------------
!
      parent_task_0: IF(MYPE==0)THEN                                      !<-- Parent task 0
!
!-----------------------------------------------------------------------
!
        RATIO=CHILD_PARENT_SPACE_RATIO
!
        CHILD_ARRAY(1,1)=ARRAY_STAGE_PARENT(1,1)                          !<-- SW corner of child's array coincides with a parent point
!
!***  Choose nearest parent point along child's southern boundary  
!
        DO I=2,IM_CHILD                                                   !<-- Move along southern boundary of child's domain
          REAL_INDX_I_PARENT=1+(I-1)*RATIO                                !<-- Exact I index of child point on parent 
          INDX_WEST=INT(REAL_INDX_I_PARENT)                               !<-- The parent point's I index west of the child's point
          INDX_EAST=INDX_WEST+1                                           !<-- The parent point's I index east of the child's point
          WEIGHT_WEST=INDX_EAST-REAL_INDX_I_PARENT                        !<-- Interpolation weight given parent's point to the west
          WEIGHT_EAST=1.-WEIGHT_WEST                                      !<-- Interpolation weight given parent's point to the east
!
          WEIGHT_MAX=MAX(WEIGHT_WEST,WEIGHT_EAST)
          IF(WEIGHT_WEST==WEIGHT_MAX)THEN
            CHILD_ARRAY(I,1)=ARRAY_STAGE_PARENT(INDX_WEST,1)
          ELSEIF(WEIGHT_EAST==WEIGHT_MAX)THEN
            CHILD_ARRAY(I,1)=ARRAY_STAGE_PARENT(INDX_EAST,1)
          ENDIF
        ENDDO
!
!***  Choose nearest parent point along child's western boundary
!
        DO J=2,JM_CHILD                                                   !<-- Move along western boundary of child's domain
          REAL_INDX_J_PARENT=1+(J-1)*RATIO                                !<-- Exact J index of child point on parent 
          INDX_SOUTH=INT(REAL_INDX_J_PARENT)                              !<-- The parent point's J index south of the child's point
          INDX_NORTH=INDX_SOUTH+1                                         !<-- The parent point's J index north of the child's point
          WEIGHT_SOUTH=INDX_NORTH-REAL_INDX_J_PARENT                      !<-- Interpolation weight of parent's point to the south
          WEIGHT_NORTH=1.-WEIGHT_SOUTH                                       !<-- Interpolation weight of parent's point to the north
          WEIGHT_MAX=MAX(WEIGHT_SOUTH,WEIGHT_NORTH)
!
          IF(WEIGHT_SOUTH==WEIGHT_MAX)THEN
            CHILD_ARRAY(1,J)=ARRAY_STAGE_PARENT(1,INDX_SOUTH)
          ELSE
            CHILD_ARRAY(1,J)=ARRAY_STAGE_PARENT(1,INDX_NORTH)
          ENDIF
!
        ENDDO
!
!-----------------------------------------------------------------------
!***  Fill in the interior of the staging array choosing the
!***  nearest parent point.
!-----------------------------------------------------------------------
!
        DO J=2,JM_CHILD
          REAL_INDX_J_PARENT=1+(J-1)*RATIO                                !<-- Exact J index of child point in parent staging region 
          INDX_SOUTH=INT(REAL_INDX_J_PARENT)                              !<-- The parent point's J index south of the child's point
          INDX_NORTH=INDX_SOUTH+1                                         !<-- The parent point's J index north of the child's point
!
          DELTA_J_NORTH=INDX_NORTH-REAL_INDX_J_PARENT
          DELTA_J_SOUTH=REAL_INDX_J_PARENT-INDX_SOUTH
!
          DO I=2,IM_CHILD
            REAL_INDX_I_PARENT=1+(I-1)*RATIO                              !<-- Exact I index of child point in parent staging region
            INDX_WEST=INT(REAL_INDX_I_PARENT)                             !<-- The parent point's I index west of the child's point
            INDX_EAST=INDX_WEST+1                                         !<-- The parent point's I index east of the child's point
!
            DELTA_I_EAST=INDX_EAST-REAL_INDX_I_PARENT
            DELTA_I_WEST=REAL_INDX_I_PARENT-INDX_WEST
!
            WEIGHT_SW=DELTA_I_EAST*DELTA_J_NORTH                          !<-- Interpolation weight of parent's point to SW 
            WEIGHT_SE=DELTA_I_WEST*DELTA_J_NORTH                          !<-- Interpolation weight of parent's point to SE
            WEIGHT_NW=DELTA_I_EAST*DELTA_J_SOUTH                          !<-- Interpolation weight of parent's point to NW
            WEIGHT_NE=DELTA_I_WEST*DELTA_J_SOUTH                          !<-- Interpolation weight of parent's point to NE
!
            WEIGHT_MAX=MAX(WEIGHT_SW,WEIGHT_SE                         & 
                          ,WEIGHT_NW,WEIGHT_NE)
!
            IF(WEIGHT_SW==WEIGHT_MAX)THEN
              CHILD_ARRAY(I,J)=ARRAY_STAGE_PARENT(INDX_WEST,INDX_SOUTH)
            ELSEIF(WEIGHT_SE==WEIGHT_MAX)THEN
              CHILD_ARRAY(I,J)=ARRAY_STAGE_PARENT(INDX_EAST,INDX_SOUTH)
            ELSEIF(WEIGHT_NW==WEIGHT_MAX)THEN
              CHILD_ARRAY(I,J)=ARRAY_STAGE_PARENT(INDX_WEST,INDX_NORTH)
            ELSEIF(WEIGHT_NE==WEIGHT_MAX)THEN
              CHILD_ARRAY(I,J)=ARRAY_STAGE_PARENT(INDX_EAST,INDX_NORTH)
            ENDIF
          ENDDO
        ENDDO
!
        DEALLOCATE(ARRAY_STAGE_PARENT)
!
!-----------------------------------------------------------------------
!
      ENDIF parent_task_0
!
!-----------------------------------------------------------------------
!
!!!   END SUBROUTINE PARENT_TO_CHILD_IFILL_ASSOC
      END SUBROUTINE PARENT_TO_CHILD_IFILL
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE PARENT_TO_CHILD_FILL_GENERAL(PARENT_ARRAY              &
!!!   SUBROUTINE PARENT_TO_CHILD_FILL        (PARENT_ARRAY              &
                                             ,NLEV                      &
                                             ,VBL_NAME                  &
                                             ,CHILD_ARRAY)
!
!-----------------------------------------------------------------------
!***  Parent tasks interpolate their data to the locations of their
!***  children's gridpoints.  The child grids are unique rotated
!***  lat/lon grids with their own centers.  The southwest H point
!***  of the child grid lies directly on an H point of the parent. 
!
!***  Only parent tasks participate in this work.
!-----------------------------------------------------------------------
!
      USE module_CONSTANTS,ONLY: PI
!
!------------------------
!***  Argument Variables
!------------------------
!
      INTEGER(kind=KINT),INTENT(IN) :: NLEV                                !<-- Vertical dimension of the data array 
!
      REAL(kind=KFPT),DIMENSION(IMS:IME,JMS:JME,1:NLEV),INTENT(IN) ::   &
                                                            PARENT_ARRAY   !<-- The parent array that will initialize the child array
!
      CHARACTER(*),INTENT(IN) :: VBL_NAME                                  !<-- The variable's name 
!
      REAL(kind=KFPT),DIMENSION(1:IM_CHILD,1:JM_CHILD,1:NLEV),INTENT(OUT) :: &
                                                            CHILD_ARRAY    !<-- Data from parent tasks interpolated to child grid
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: I,I_END,ISTART,J,J_END,JSTART               &
                           ,KOUNT,L,NIJ,NN,NTOT                         &
                           ,NUM_DATA                                    &
                           ,NUM_CHILD_POINTS                            &
                           ,NUM_IJ                                      &
                           ,NUM_POINTS_REMOTE
!
      INTEGER(kind=KINT) :: I_PARENT_SW,I_PARENT_SE                     &
                           ,I_PARENT_NW,I_PARENT_NE                     &
                           ,J_PARENT_SW,J_PARENT_SE                     &
                           ,J_PARENT_NW,J_PARENT_NE
!
      INTEGER(kind=KINT) :: IERR,IPE,ISTAT
!
      INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE :: CHILD_POINT_INDICES &
                                                    ,IJ_REMOTE
!
      INTEGER,DIMENSION(MPI_STATUS_SIZE) :: JSTAT
!
      REAL(kind=KFPT) :: CHILD_LATD_ON_PARENT                           &
                        ,CHILD_LOND_ON_PARENT                           &
                        ,DEG_TO_RAD                                     &
                        ,DIST                                           &
                        ,R_DLMD,R_DPHD                                  &
                        ,REAL_I_PARENT                                  &
                        ,REAL_J_PARENT                                  &
                        ,RLATD_SW,RLOND_SW                              &
                        ,RLATD_SE,RLOND_SE                              &
                        ,RLATD_NW,RLOND_NW                              &
                        ,RLATD_NE,RLOND_NE                              &
                        ,SUM,SUM_RECIP                                  &
                        ,WEIGHT_SW,WEIGHT_SE                            &
                        ,WEIGHT_NW,WEIGHT_NE                            &
                        ,WEIGHT_SUM,WEIGHT_SUM_RECIP
!
      REAL(kind=KFPT),DIMENSION(4) :: RLATD,RLOND,WGT
!
      REAL(kind=KFPT),DIMENSION(:),ALLOCATABLE :: CHILD_STRING          &
                                                 ,DATA_REMOTE
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      R_DPHD=1./DPHD_PARENT
      R_DLMD=1./DLMD_PARENT
      DEG_TO_RAD=PI/180.
!
      NUM_CHILD_POINTS=0
!
      ISTART=1
      JSTART=1
      IF(GLOBAL)THEN
        ISTART=2
        JSTART=2
      ENDIF
!
!-----------------------------------------------------------------------
!***  Each parent task is responsible for searching the parent domain 
!***  extending from ITS and JTS to ITE+1 and JTE+1.  Those latter +1's
!***  are needed in order to reach the next gridpoint in each direction.
!***  We cannot go outside the full domain of course plus the wind 
!***  points have no values at IDE and JDE.
!-----------------------------------------------------------------------
!
      IF(VBL_NAME=='Uwind'.OR.VBL_NAME=='Vwind')THEN
        I_END=MIN(ITE+1,IDE-1)
        J_END=MIN(JTE+1,JDE-1)
      ELSE
        I_END=MIN(ITE+1,IDE)
        J_END=MIN(JTE+1,JDE)
      ENDIF
!
!-----------------------------------------------------------------------
!
      NTOT=2*IM_CHILD*JM_CHILD
      ALLOCATE(CHILD_POINT_INDICES(1:NTOT),stat=ISTAT)
      IF(ISTAT/=0)WRITE(0,*)' Failed to allocate CHILD_POINT_INDICES in PARENT_TO_CHILD_FILL_GENERAL'
!
      NTOT=IM_CHILD*JM_CHILD*NLEV
      ALLOCATE(CHILD_STRING(1:NTOT),stat=ISTAT)
      IF(ISTAT/=0)WRITE(0,*)' Failed to allocate CHILD_STRING in PARENT_TO_CHILD_FILL_GENERAL'
!
!-----------------------------------------------------------------------
!***  Compute the parent's lat/lon of each child gridpoint in order to
!***  determine if that gridpoint lies on a given parent task. 
!***  Save the I,J of each gridpoint found since ultimately parent
!***  task 0 will need that to properly place all interpolated child 
!***  data onto the full child grid for writing out.
!-----------------------------------------------------------------------
!
      NN=0
!
      DO J=1,JM_CHILD
      DO I=1,IM_CHILD
!
        CALL CONVERT_IJ_TO_LATLON(I,J                                   &  !<-- A point on the child grid
                                 ,IM_CHILD,JM_CHILD                     &  !<-- Dimensions of child grid
                                 ,TPH0D_CHILD,TLM0D_CHILD               &  !<-- Parent lat/lon (deg) of child grid central point
                                 ,DPHD_CHILD,DLMD_CHILD                 &  !<-- Angular grid increments (deg) on child grid
                                 ,CHILD_LATD_ON_PARENT                  &  !<-- Parent latitude of child point
                                 ,CHILD_LOND_ON_PARENT )                   !<-- Parent longitude of child point
!
        REAL_I_PARENT=(CHILD_LOND_ON_PARENT-WBD_PARENT)*R_DLMD+ISTART      !<-- REAL I index of child point on parent grid
        REAL_J_PARENT=(CHILD_LATD_ON_PARENT-SBD_PARENT)*R_DPHD+JSTART      !<-- REAL J index of child point on parent grid
!
!-----------------------------------------------------------------------
!
        IF(REAL(ITS)<=REAL_I_PARENT.AND.REAL(I_END)>REAL_I_PARENT.AND.  &  !<-- Is child gridpoint on this parent task?
           REAL(JTS)<=REAL_J_PARENT.AND.REAL(J_END)>REAL_J_PARENT)THEN     !<--
!
          NUM_CHILD_POINTS=NUM_CHILD_POINTS+1                              !<-- Add up number of child points on this parent task
!
          CHILD_POINT_INDICES(2*NUM_CHILD_POINTS-1)=I                      !<-- Save I index of this child
          CHILD_POINT_INDICES(2*NUM_CHILD_POINTS  )=J                      !<-- Save J index of this child point
!
!-----------------------------------------------------------------------
!***  Compute the distance from the child point location to each of
!***  the four surrounding parent points and generate the bilinear
!***  interpolation weights.
!***  The indices 1-->4 indicate the parent points to the SW, SE,
!***  NW, and NE in that order.
!-----------------------------------------------------------------------
!
          I_PARENT_SW=INT(REAL_I_PARENT)
          J_PARENT_SW=INT(REAL_J_PARENT)
          RLATD(1)=(J_PARENT_SW-ROW_0)*DPHD_PARENT                         !<-- Parent latitude (deg) of parent point SW of child point
          RLOND(1)=(I_PARENT_SW-COL_0)*DLMD_PARENT                         !<-- Parent longitude (deg) of parent point SW of child point
!
          I_PARENT_SE=I_PARENT_SW+1         
          J_PARENT_SE=J_PARENT_SW         
          RLATD(2)=RLATD(1)                                                !<-- SE and SW on same line of parent latitude
          RLOND(2)=RLOND(1)+DLMD_PARENT                                    !<-- SE is one gridpoint east of SW parent point
!
          I_PARENT_NW=I_PARENT_SW           
          J_PARENT_NW=J_PARENT_SW+1       
          RLATD(3)=RLATD(1)+DPHD_PARENT                                    !<-- NW is one gridpoint north of SW parent point
          RLOND(3)=RLOND(1)                                                !<-- NW and SW on same line of parent longitude
!
          I_PARENT_NE=I_PARENT_SE           
          J_PARENT_NE=J_PARENT_NW         
          RLATD(4)=RLATD(3)                                                !<-- NE and NW on same line of parent latitude
          RLOND(4)=RLOND(2)                                                !<-- NE and SE on same line of parent longitude
!
          SUM=0.
!
          DO N=1,4                                                         !<-- Loop over SW, SE, NW, and NE parent points
!
            CALL DISTANCE_ON_SPHERE(CHILD_LATD_ON_PARENT*DEG_TO_RAD    &   !<-- Parent latitiude (deg) of child gridpoint
                                   ,CHILD_LOND_ON_PARENT*DEG_TO_RAD    &   !<-- Parent latitiude (deg) of child gridpoint
                                   ,RLATD(N)*DEG_TO_RAD                &   !<-- Latitude (deg) of surrounding parent point N
                                   ,RLOND(N)*DEG_TO_RAD                &   !<-- Longitude (deg) of surrounding parent point N
                                   ,DIST )                                 !<-- Distance (radians) from child point to parent point N
!
            WGT(N)=1./DIST
            SUM=SUM+WGT(N)
!
          ENDDO
!
          SUM_RECIP=1./SUM
!
!-----------------------------------------------------------------------
!***  The bilinear interpolation weights of the four parent points
!***  surrounding the child point.
!-----------------------------------------------------------------------
!
          WEIGHT_SW=WGT(1)*SUM_RECIP
          WEIGHT_SE=WGT(2)*SUM_RECIP
          WEIGHT_NW=WGT(3)*SUM_RECIP
          WEIGHT_NE=WGT(4)*SUM_RECIP
!
          IF(ABS(PARENT_ARRAY(I_PARENT_SW,J_PARENT_SW,1))<1.E-12)WEIGHT_SW=0.
          IF(ABS(PARENT_ARRAY(I_PARENT_SE,J_PARENT_SE,1))<1.E-12)WEIGHT_SE=0.
          IF(ABS(PARENT_ARRAY(I_PARENT_NW,J_PARENT_NW,1))<1.E-12)WEIGHT_NW=0.
          IF(ABS(PARENT_ARRAY(I_PARENT_NE,J_PARENT_NE,1))<1.E-12)WEIGHT_NE=0.
          WEIGHT_SUM=WEIGHT_SW+WEIGHT_SE+WEIGHT_NW+WEIGHT_NE
          WEIGHT_SUM_RECIP=1./WEIGHT_SUM
!
          DO L=1,NLEV
            NN=NN+1
!
            CHILD_STRING(NN)=WEIGHT_SW*PARENT_ARRAY(I_PARENT_SW,J_PARENT_SW,L)  & !<-- Value at points on child's grid
                            +WEIGHT_SE*PARENT_ARRAY(I_PARENT_SE,J_PARENT_SE,L)  &
                            +WEIGHT_NW*PARENT_ARRAY(I_PARENT_NW,J_PARENT_NW,L)  &
                            +WEIGHT_NE*PARENT_ARRAY(I_PARENT_NE,J_PARENT_NE,L)
!
            IF(WEIGHT_SUM<0.99.AND.WEIGHT_SUM>0.01)THEN
              CHILD_STRING(NN)=CHILD_STRING(NN)*WEIGHT_SUM_RECIP           !<-- Normalize if some weights are zero (e.g., coastal land Temp)
            ENDIF
!
            IF(VBL_NAME=='SeaMask')THEN
              IF(CHILD_STRING(NN)>=0.5)CHILD_STRING(NN)=1.0
              IF(CHILD_STRING(NN)< 0.5)CHILD_STRING(NN)=0.0
            ENDIF
!
          ENDDO
!
!-----------------------------------------------------------------------
!
        ENDIF
!
!-----------------------------------------------------------------------
!
      ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
!***  Each parent task that contains a child grid point has now done
!***  the horizontal interpolation of the parent variable to the child.
!***  Now parent task 0 receives all the interpolated data from the
!***  other parent tasks.
!-----------------------------------------------------------------------
!
      data_fill: IF(MYPE==0)THEN
!
!-----------------------------------------------------------------------
!
        remote_tasks: DO IPE=1,NUM_PES_PARENT-1
!        
          CALL MPI_RECV(NUM_POINTS_REMOTE                               &  !<-- # of child points on parent task IPE
                       ,1                                               &  !<-- Total words received
                       ,MPI_INTEGER                                     &  !<-- Datatype
                       ,IPE                                             &  !<-- Receive from this parent task
                       ,IPE                                             &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- The MPI communicator
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR )
!
          IF(NUM_POINTS_REMOTE==0)CYCLE remote_tasks
!
          NUM_IJ  =2*NUM_POINTS_REMOTE
          NUM_DATA=NUM_POINTS_REMOTE*NLEV
!
          ALLOCATE(DATA_REMOTE(1:NUM_DATA))
          ALLOCATE(IJ_REMOTE  (1:NUM_IJ  ))
!
          CALL MPI_RECV(DATA_REMOTE                                     &  !<-- Interpolated data on child grid from parent task IPE
                       ,NUM_DATA                                        &  !<-- Total words received
                       ,MPI_REAL                                        &  !<-- Datatype
                       ,IPE                                             &  !<-- Receive from this parent task
                       ,IPE                                             &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- The MPI communicator
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR )
!
          CALL MPI_RECV(IJ_REMOTE                                       &  !<-- Interpolated data on child grid from parent task IPE
                       ,NUM_IJ                                          &  !<-- Total words received
                       ,MPI_INTEGER                                     &  !<-- Datatype
                       ,IPE                                             &  !<-- Receive from this parent task
                       ,IPE                                             &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- The MPI communicator
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR )
!
!-----------------------------------------------------------------------
!***  Parent task 0 fills in section of child data from remote 
!***  parent task IPE.
!-----------------------------------------------------------------------
!
          KOUNT=0
!
          DO L=1,NLEV
            DO NIJ=1,NUM_POINTS_REMOTE
              KOUNT=KOUNT+1
              I=IJ_REMOTE(2*NIJ-1)
              J=IJ_REMOTE(2*NIJ  )
              CHILD_ARRAY(I,J,L)=DATA_REMOTE(KOUNT)
            ENDDO
          ENDDO
!
!-----------------------------------------------------------------------
!
          DEALLOCATE(DATA_REMOTE)
          DEALLOCATE(IJ_REMOTE)
!
!-----------------------------------------------------------------------
!
        ENDDO remote_tasks
!
!-----------------------------------------------------------------------
!***  Finally parent task 0 fills in its own section of the child array.
!-----------------------------------------------------------------------
!
        IF(NUM_CHILD_POINTS>0)THEN
!
          KOUNT=0          
          DO L=1,NLEV
            DO N=1,NUM_CHILD_POINTS
              KOUNT=KOUNT+1
              I=CHILD_POINT_INDICES(2*N-1)
              J=CHILD_POINT_INDICES(2*N  )
              CHILD_ARRAY(I,J,L)=CHILD_STRING(KOUNT)
            ENDDO
          ENDDO
!
        ENDIF
!
!-----------------------------------------------------------------------
!
      ELSE data_fill
!
!-----------------------------------------------------------------------
!***  Remote parent tasks send their sections of interpolated child
!***  data to parent task 0.
!-----------------------------------------------------------------------
!
        CALL MPI_SEND(NUM_CHILD_POINTS                                  &  !<-- # of child points on this parent task
                     ,1                                                 &  
                     ,MPI_INTEGER                                       &  !<-- Datatype
                     ,0                                                 &  !<-- Send to parent task 0
                     ,MYPE                                              &  !<-- MPI tag
                     ,COMM_MY_DOMAIN                                    &  !<-- The MPI communicator
                     ,IERR )
!
        IF(NUM_CHILD_POINTS>0)THEN
          NUM_IJ  =2*NUM_CHILD_POINTS
          NUM_DATA=NUM_CHILD_POINTS*NLEV
!
          CALL MPI_SEND(CHILD_STRING                                    &  !<-- Interpolated data on child grid for this parent task
                       ,NUM_DATA                                        &  !<-- Total words sent
                       ,MPI_REAL                                        &  !<-- Datatype
                       ,0                                               &  !<-- Send to parent task 0
                       ,MYPE                                            &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- The MPI communicator
                       ,IERR )
!
          CALL MPI_SEND(CHILD_POINT_INDICES                             &  !<-- Indices of child points for this parent task
                       ,NUM_IJ                                          &  !<-- Total words sent
                       ,MPI_INTEGER                                     &  !<-- Datatype
                       ,0                                               &  !<-- Send to parent task 0
                       ,MYPE                                            &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- The MPI communicator
                       ,IERR )
!
        ENDIF
!
!-----------------------------------------------------------------------
!
      ENDIF data_fill
!
!-----------------------------------------------------------------------
!
      DEALLOCATE(CHILD_POINT_INDICES)
      DEALLOCATE(CHILD_STRING)
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE PARENT_TO_CHILD_FILL_GENERAL
!!!   END SUBROUTINE PARENT_TO_CHILD_FILL
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE PARENT_TO_CHILD_IFILL_GENERAL(PARENT_ARRAY             &
!!!   SUBROUTINE PARENT_TO_CHILD_IFILL        (PARENT_ARRAY             &
                                              ,VBL_NAME                 &
                                              ,CHILD_ARRAY)
!
!-----------------------------------------------------------------------
!***  Parent tasks interpolate their data to the locations of their
!***  children's gridpoints.  The child grids are unique rotated
!***  lat/lon grids with their own centers.  The southwest H point
!***  of the child grid lies directly on an H point of the parent. 
!
!***  Only parent tasks participate in this work.
!-----------------------------------------------------------------------
!
      USE module_CONSTANTS,ONLY: PI
!
!------------------------
!***  Argument Variables
!------------------------
!
      INTEGER(kind=KINT),DIMENSION(IMS:IME,JMS:JME),INTENT(IN) ::        &
                                                           PARENT_ARRAY    !<-- The parent array that will initialize the child array
!
      CHARACTER(*),INTENT(IN) :: VBL_NAME                                  !<-- The variable's name 
!
      INTEGER(kind=KINT),DIMENSION(1:IM_CHILD,1:JM_CHILD),INTENT(OUT) :: &
                                                           CHILD_ARRAY     !<-- Data from parent tasks interpolated to child grid
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: I,I_END,ISTART,J,J_END,JSTART               &
                           ,KOUNT,NIJ,NN,NTOT                           &
                           ,NUM_DATA                                    &
                           ,NUM_CHILD_POINTS                            &
                           ,NUM_IJ                                      &
                           ,NUM_POINTS_REMOTE
!
      INTEGER(kind=KINT) :: I_PARENT_SW,I_PARENT_SE                     &
                           ,I_PARENT_NW,I_PARENT_NE                     &
                           ,J_PARENT_SW,J_PARENT_SE                     &
                           ,J_PARENT_NW,J_PARENT_NE
!
      INTEGER(kind=KINT) :: IERR,IPE,ISTAT
!
      INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE :: CHILD_POINT_INDICES &
                                                    ,CHILD_STRING        &
                                                    ,DATA_REMOTE         &
                                                    ,IJ_REMOTE  
!
      INTEGER,DIMENSION(MPI_STATUS_SIZE) :: JSTAT
!
      REAL(kind=KFPT) :: CHILD_LATD_ON_PARENT                           &
                        ,CHILD_LOND_ON_PARENT                           &
                        ,DEG_TO_RAD                                     &
                        ,DIST                                           &
                        ,R_DLMD,R_DPHD                                  &
                        ,REAL_I_PARENT                                  &
                        ,REAL_J_PARENT                                  &
                        ,RLATD_SW,RLOND_SW                              &
                        ,RLATD_SE,RLOND_SE                              &
                        ,RLATD_NW,RLOND_NW                              &
                        ,RLATD_NE,RLOND_NE                              &
                        ,SUM,SUM_RECIP                                  &
                        ,WEIGHT_SW,WEIGHT_SE                            &
                        ,WEIGHT_NW,WEIGHT_NE                            &
                        ,WEIGHT_MAX
!
      REAL(kind=KFPT),DIMENSION(4) :: RLATD,RLOND,WGT
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      R_DPHD=1./DPHD_PARENT
      R_DLMD=1./DLMD_PARENT
      DEG_TO_RAD=PI/180.
!
      NUM_CHILD_POINTS=0
!
      ISTART=1
      JSTART=1
      IF(GLOBAL)THEN
        ISTART=2
        JSTART=2
      ENDIF
!
!-----------------------------------------------------------------------
!***  Each parent task is responsible for searching the parent domain 
!***  extending from ITS and JTS to ITE+1 and JTE+1.  Those latter +1's
!***  are needed in order to reach the next gridpoint in each direction.
!***  We cannot go outside the full domain of course.
!-----------------------------------------------------------------------
!
      I_END=MIN(ITE+1,IDE)
      J_END=MIN(JTE+1,JDE)
!
!-----------------------------------------------------------------------
!
      NTOT=2*IM_CHILD*JM_CHILD
      ALLOCATE(CHILD_POINT_INDICES(1:NTOT),stat=ISTAT)
      IF(ISTAT/=0)WRITE(0,*)' Failed to allocate CHILD_POINT_INDICES in PARENT_TO_CHILD_IFILL_GENERAL'
!
      NTOT=IM_CHILD*JM_CHILD
      ALLOCATE(CHILD_STRING(1:NTOT),stat=ISTAT)
      IF(ISTAT/=0)WRITE(0,*)' Failed to allocate CHILD_STRING in PARENT_TO_CHILD_IFILL_GENERAL'
!
!-----------------------------------------------------------------------
!***  Compute the parent's lat/lon of each child gridpoint in order to
!***  determine if that gridpoint lies on a given parent task. 
!***  Save the I,J of each gridpoint found since ultimately parent
!***  task 0 will need that to properly place all interpolated child 
!***  data onto the full child grid for writing out.
!-----------------------------------------------------------------------
!
      NN=0
!
      DO J=1,JM_CHILD
      DO I=1,IM_CHILD
!
        CALL CONVERT_IJ_TO_LATLON(I,J                                   &  !<-- A point on the child grid
                                 ,IM_CHILD,JM_CHILD                     &  !<-- Dimensions of child grid
                                 ,TPH0D_CHILD,TLM0D_CHILD               &  !<-- Parent lat/lon (deg) of child grid central point
                                 ,DPHD_CHILD,DLMD_CHILD                 &  !<-- Angular grid increments (deg) on child grid
                                 ,CHILD_LATD_ON_PARENT                  &  !<-- Parent latitude of child point
                                 ,CHILD_LOND_ON_PARENT )                   !<-- Parent longitude of child point
!
        REAL_I_PARENT=(CHILD_LOND_ON_PARENT-WBD_PARENT)*R_DLMD+ISTART      !<-- REAL I index of child point on parent grid
        REAL_J_PARENT=(CHILD_LATD_ON_PARENT-SBD_PARENT)*R_DPHD+JSTART      !<-- REAL J index of child point on parent grid
!
!-----------------------------------------------------------------------
!
        IF(REAL(ITS)<=REAL_I_PARENT.AND.REAL(I_END)>REAL_I_PARENT.AND.  &  !<-- Is child gridpoint on this parent task?
           REAL(JTS)<=REAL_J_PARENT.AND.REAL(J_END)>REAL_J_PARENT)THEN     !<--
!
          NUM_CHILD_POINTS=NUM_CHILD_POINTS+1                              !<-- Add up number of child points on this parent task
!
          CHILD_POINT_INDICES(2*NUM_CHILD_POINTS-1)=I                      !<-- Save I index of this child
          CHILD_POINT_INDICES(2*NUM_CHILD_POINTS  )=J                      !<-- Save J index of this child point
!
!-----------------------------------------------------------------------
!***  Compute the distance from the child point location to each of
!***  the four surrounding parent points and generate the bilinear
!***  interpolation weights.
!***  The indices 1-->4 indicate the parent points to the SW, SE,
!***  NW, and NE in that order.
!-----------------------------------------------------------------------
!
          I_PARENT_SW=INT(REAL_I_PARENT)
          J_PARENT_SW=INT(REAL_J_PARENT)
          RLATD(1)=(J_PARENT_SW-ROW_0)*DPHD_PARENT                         !<-- Parent latitude (deg) of parent point SW of child point
          RLOND(1)=(I_PARENT_SW-COL_0)*DLMD_PARENT                         !<-- Parent longitude (deg) of parent point SW of child point
!
          I_PARENT_SE=I_PARENT_SW+1         
          J_PARENT_SE=J_PARENT_SW         
          RLATD(2)=RLATD(1)                                                !<-- SE and SW on same line of parent latitude
          RLOND(2)=RLOND(1)+DLMD_PARENT                                    !<-- SE is one gridpoint east of SW parent point
!
          I_PARENT_NW=I_PARENT_SW           
          J_PARENT_NW=J_PARENT_SW+1       
          RLATD(3)=RLATD(1)+DPHD_PARENT                                    !<-- NW is one gridpoint north of SW parent point
          RLOND(3)=RLOND(1)                                                !<-- NW and SW on same line of parent longitude
!
          I_PARENT_NE=I_PARENT_SE           
          J_PARENT_NE=J_PARENT_NW         
          RLATD(4)=RLATD(3)                                                !<-- NE and NW on same line of parent latitude
          RLOND(4)=RLOND(2)                                                !<-- NE and SE on same line of parent longitude
!
          SUM=0.
!
          DO N=1,4                                                         !<-- Loop over SW, SE, NW, and NE parent points
!
            CALL DISTANCE_ON_SPHERE(CHILD_LATD_ON_PARENT*DEG_TO_RAD    &   !<-- Parent latitiude (deg) of child gridpoint
                                   ,CHILD_LOND_ON_PARENT*DEG_TO_RAD    &   !<-- Parent latitiude (deg) of child gridpoint
                                   ,RLATD(N)*DEG_TO_RAD                &   !<-- Latitude (deg) of surrounding parent point N
                                   ,RLOND(N)*DEG_TO_RAD                &   !<-- Longitude (deg) of surrounding parent point N
                                   ,DIST )                                 !<-- Distance (radians) from child point to parent point N
!
            WGT(N)=1./DIST
            SUM=SUM+WGT(N)
!
          ENDDO
!
          SUM_RECIP=1./SUM
!
!-----------------------------------------------------------------------
!***  The bilinear interpolation weights of the four parent points
!***  surrounding the child point.
!-----------------------------------------------------------------------
!
          WEIGHT_SW=WGT(1)*SUM_RECIP
          WEIGHT_SE=WGT(2)*SUM_RECIP
          WEIGHT_NW=WGT(3)*SUM_RECIP
          WEIGHT_NE=WGT(4)*SUM_RECIP
          WEIGHT_MAX=MAX(WEIGHT_SW,WEIGHT_SE,WEIGHT_NW,WEIGHT_NE)
!
!-----------------------------------------------------------------------
!***  Using the bilinear interpolation weights, assign the value of
!***  the nearest parent point to the child point.
!-----------------------------------------------------------------------
!
          NN=NN+1
!
          IF(WEIGHT_SW==WEIGHT_MAX)THEN
            CHILD_STRING(NN)=PARENT_ARRAY(I_PARENT_SW,J_PARENT_SW)
!     write(0,*)' SW parent=',PARENT_ARRAY(I_PARENT_SW,J_PARENT_SW),' nn=',nn,' i=',i,' j=',j
          ELSEIF(WEIGHT_SE==WEIGHT_MAX)THEN
            CHILD_STRING(NN)=PARENT_ARRAY(I_PARENT_SE,J_PARENT_SE)
!     write(0,*)' SE parent=',PARENT_ARRAY(I_PARENT_SE,J_PARENT_SE),' nn=',nn,' i=',i,' j=',j
          ELSEIF(WEIGHT_NW==WEIGHT_MAX)THEN
            CHILD_STRING(NN)=PARENT_ARRAY(I_PARENT_NW,J_PARENT_NW)
!     write(0,*)' NW parent=',PARENT_ARRAY(I_PARENT_NW,J_PARENT_NW),' nn=',nn,' i=',i,' j=',j
          ELSEIF(WEIGHT_NE==WEIGHT_MAX)THEN
            CHILD_STRING(NN)=PARENT_ARRAY(I_PARENT_NE,J_PARENT_NE)
!     write(0,*)' NE parent=',PARENT_ARRAY(I_PARENT_NE,J_PARENT_NE),' nn=',nn,' i=',i,' j=',j
          ENDIF
!     if(i==01.and.j==01.and.vbl_name=='ISLTYP')then
!       write(0,*)' parent interp value to ISLTYP is ',CHILD_STRING(NN),' nn=',nn
!     endif
!     if(vbl_name=='ISLTYP'.and.child_string(nn)<1.and.sea_mask(i,j)<0.5)then
!       write(0,*)' Parent creating bad value of ISLTYP=',CHILD_STRING(NN),' nn=',nn,' i=',i,' j=',j &
!                ,' SeaMask=',sea_mask(i,j)
!     endif
        
!
!-----------------------------------------------------------------------
!
        ENDIF
!
!-----------------------------------------------------------------------
!
      ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
!***  Each parent task that contains a child grid point has now done
!***  the horizontal interpolation of the parent variable to the child.
!***  Now parent task 0 receives all the interpolated data from the
!***  other parent tasks.
!-----------------------------------------------------------------------
!
      data_fill: IF(MYPE==0)THEN
!
!-----------------------------------------------------------------------
!
        remote_tasks: DO IPE=1,NUM_PES_PARENT-1
!        
          CALL MPI_RECV(NUM_POINTS_REMOTE                               &  !<-- # of child points on parent task IPE
                       ,1                                               &  !<-- Total words received
                       ,MPI_INTEGER                                     &  !<-- Datatype
                       ,IPE                                             &  !<-- Receive from this parent task
                       ,IPE                                             &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- The MPI communicator
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR )
!
          IF(NUM_POINTS_REMOTE==0)CYCLE remote_tasks
!
          NUM_IJ  =2*NUM_POINTS_REMOTE
          NUM_DATA=NUM_POINTS_REMOTE
!
          ALLOCATE(DATA_REMOTE(1:NUM_DATA))
          ALLOCATE(IJ_REMOTE  (1:NUM_IJ  ))
!
          CALL MPI_RECV(DATA_REMOTE                                     &  !<-- Interpolated data on child grid from parent task IPE
                       ,NUM_DATA                                        &  !<-- Total words received
                       ,MPI_INTEGER                                     &  !<-- Datatype
                       ,IPE                                             &  !<-- Receive from this parent task
                       ,IPE                                             &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- The MPI communicator
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR )
!
          CALL MPI_RECV(IJ_REMOTE                                       &  !<-- Interpolated data on child grid from parent task IPE
                       ,NUM_IJ                                          &  !<-- Total words received
                       ,MPI_INTEGER                                     &  !<-- Datatype
                       ,IPE                                             &  !<-- Receive from this parent task
                       ,IPE                                             &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- The MPI communicator
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR )
!
!-----------------------------------------------------------------------
!***  Parent task 0 fills in section of child data from remote 
!***  parent task IPE.
!-----------------------------------------------------------------------
!
          KOUNT=0
!
          DO NIJ=1,NUM_POINTS_REMOTE
            KOUNT=KOUNT+1
            I=IJ_REMOTE(2*NIJ-1)
            J=IJ_REMOTE(2*NIJ  )
            CHILD_ARRAY(I,J)=DATA_REMOTE(KOUNT)
!     write(0,*)' new child i=',i,' j=',j,' kount=',kount,' data=',data_remote(kount)
!     if(vbl_name=='ISLTYP'.and.i==01.and.j==01)then
!       write(0,*)' parent fills ISLTYP with ',DATA_REMOTE(KOUNT),' kount=',kount
!     endif
          ENDDO
!
!-----------------------------------------------------------------------
!
          DEALLOCATE(DATA_REMOTE)
          DEALLOCATE(IJ_REMOTE)
!
!-----------------------------------------------------------------------
!
        ENDDO remote_tasks
!
!-----------------------------------------------------------------------
!***  Finally parent task 0 fills in its own section of the child array.
!-----------------------------------------------------------------------
!
        IF(NUM_CHILD_POINTS>0)THEN
!
          KOUNT=0          
          DO N=1,NUM_CHILD_POINTS
            KOUNT=KOUNT+1
            I=CHILD_POINT_INDICES(2*N-1)
            J=CHILD_POINT_INDICES(2*N  )
            CHILD_ARRAY(I,J)=CHILD_STRING(KOUNT)
          ENDDO
!
        ENDIF
!
!-----------------------------------------------------------------------
!
      ELSE data_fill
!
!-----------------------------------------------------------------------
!***  Remote parent tasks send their sections of interpolated child
!***  data to parent task 0.
!-----------------------------------------------------------------------
!
        CALL MPI_SEND(NUM_CHILD_POINTS                                  &  !<-- # of child points on this parent task
                     ,1                                                 &  
                     ,MPI_INTEGER                                       &  !<-- Datatype
                     ,0                                                 &  !<-- Send to parent task 0
                     ,MYPE                                              &  !<-- MPI tag
                     ,COMM_MY_DOMAIN                                    &  !<-- The MPI communicator
                     ,IERR )
!
        IF(NUM_CHILD_POINTS>0)THEN
          NUM_IJ  =2*NUM_CHILD_POINTS
          NUM_DATA=NUM_CHILD_POINTS
!
          CALL MPI_SEND(CHILD_STRING                                    &  !<-- Interpolated data on child grid for this parent task
                       ,NUM_DATA                                        &  !<-- Total words sent
                       ,MPI_INTEGER                                     &  !<-- Datatype
                       ,0                                               &  !<-- Send to parent task 0
                       ,MYPE                                            &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- The MPI communicator
                       ,IERR )
!
          CALL MPI_SEND(CHILD_POINT_INDICES                             &  !<-- Indices of child points for this parent task
                       ,NUM_IJ                                          &  !<-- Total words sent
                       ,MPI_INTEGER                                     &  !<-- Datatype
                       ,0                                               &  !<-- Send to parent task 0
                       ,MYPE                                            &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- The MPI communicator
                       ,IERR )
!
        ENDIF
!
!-----------------------------------------------------------------------
!
      ENDIF data_fill
!
!-----------------------------------------------------------------------
!
      DEALLOCATE(CHILD_POINT_INDICES)
      DEALLOCATE(CHILD_STRING)
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE PARENT_TO_CHILD_IFILL_GENERAL
!!!   END SUBROUTINE PARENT_TO_CHILD_IFILL
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE PARENT_TO_CHILD_INIT_NMM
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE CONVERT_IJ_TO_LATLON  (I_INDEX                         &
                                       ,J_INDEX                         &
                                       ,IM                              &
                                       ,JM                              &
                                       ,TPH0D                           &
                                       ,TLM0D                           &
                                       ,DPHD                            &
                                       ,DLMD                            &
                                       ,RLATD                           &
                                       ,RLOND)
!
!-----------------------------------------------------------------------
!***  Given the (I,J) of mass points on an Arakawa B-Grid,
!***  compute the latitudes and longitudes before rotation.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      INTEGER(kind=KINT),INTENT(IN) :: I_INDEX                          &  !<-- I value on the grid
                                      ,J_INDEX                          &  !<-- J value on the grid
                                      ,IM                               &  !<-- Full I dimension
                                      ,JM                                  !<-- Full J dimension
!
      REAL(kind=KFPT),INTENT(IN) :: DPHD                                &  !<-- Latitude grid increment (degrees)
                                   ,DLMD                                &  !<-- Longitude grid increment (degrees)
                                   ,TPH0D                               &  ! Central latitude (deg, positive north), unrotated system
                                   ,TLM0D                                  ! Central longitude (deg, positive east), unrotated system
!
      REAL(kind=KFPT),INTENT(OUT) :: RLATD                              &  !<-- Latitude (deg, positive north) of point, unrotated system
                                    ,RLOND                                 !<-- Longitude (deg, positive east) of point, unrotated system
!
!-----------------------------------------------------------------------
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: I,IEND,ISTART,J,JEND,JSTART
!
      REAL(kind=KDBL) :: ARG1,ARG2,COL_MID,D2R,FCTR,GLATR,GLATD,GLOND   &
                        ,HALF,ONE,PI,R2D,ROW_MID,TLATD,TLOND            &
                        ,TLATR,TLONR,TLM0,TPH0
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!***  Convert from transformed grid location (I,J) 
!***  to geographic coordinates (degrees).
!-----------------------------------------------------------------------
!
      ONE=1.0
      HALF=1./2.
      PI=DACOS(-ONE)
      D2R=PI/180.
      R2D=1./D2R
      TPH0=TPH0D*D2R
      TLM0=TLM0D*D2R
!
      ROW_MID=(JM+ONE)*HALF
      COL_MID=(IM+ONE)*HALF
!
!-----------------------------------------------------------------------
!
      J=J_INDEX
      I=I_INDEX
!
!-----------------------------------------------------------------------
!***  Find the rotated latitude (positive north) and 
!***  longitude (positive east).
!-----------------------------------------------------------------------
!
      TLATD=(J-ROW_MID)*DPHD
      TLOND=(I-COL_MID)*DLMD
!
!     WRITE(0,50)I,J,TLATD,TLOND
   50 FORMAT(' I=',I4,' J=',I4,' ROTATED LATITUDE IS',F8.3              &
                              ,4X,'LONGITUDE IS',F8.3)
!
!-----------------------------------------------------------------------
!***  Now convert to geographic latitude (positive north) and
!***  longitude (positive west) in degrees.
!-----------------------------------------------------------------------
!
      TLATR=TLATD*D2R
      TLONR=TLOND*D2R
      ARG1=SIN(TLATR)*COS(TPH0)+COS(TLATR)*SIN(TPH0)*COS(TLONR)
      GLATR=ASIN(ARG1)
!
      GLATD=GLATR*R2D
!
      ARG2=DCOS(TLATR)*DCOS(TLONR)/(DCOS(GLATR)*DCOS(TPH0))-            &
           DTAN(GLATR)*DTAN(TPH0)
      IF(ABS(ARG2)>1.)ARG2=ABS(ARG2)/ARG2
      FCTR=1.
      IF(TLOND>0.)FCTR=1.
      IF(TLOND>180.)FCTR=-1.
!
      GLOND=-TLM0D+FCTR*DACOS(ARG2)*R2D
!
!     WRITE(6,100)I,J,GLATD,GLOND
  100 FORMAT(' I=',I3,' J=',I3                                          &
            ,'  PARENT LATITUDE=',F9.5,'  LONGITUDE=',F10.5)
!-----------------------------------------------------------------------
!
      RLATD=GLATD
      RLOND=-GLOND
      IF(RLOND<-180.)RLOND=RLOND+360.
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE CONVERT_IJ_TO_LATLON
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE REAL_IJ_TO_LATLON (I_INDEX                             &
                                   ,J_INDEX                             &
                                   ,IM                                  &
                                   ,JM                                  &
                                   ,TPH0                                &
                                   ,TLM0                                &
                                   ,DPH                                 &
                                   ,DLM                                 &
                                   ,RLAT                                &
                                   ,RLON )
!
!-----------------------------------------------------------------------
!***  Given the (I,J) of mass points on an Arakawa B-Grid, compute
!***  the latitudes and longitudes on the given projection.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      INTEGER(kind=KINT),INTENT(IN) :: IM                               &  !<-- Full I dimension
                                      ,JM                                  !<-- Full J dimension
!
      REAL(kind=KFPT),INTENT(IN) :: I_INDEX                             &  !<-- Real I value on the grid
                                   ,J_INDEX                             &  !<-- Real J value on the grid
                                   ,DPH                                 &  !<-- Latitude grid increment (radians)
                                   ,DLM                                 &  !<-- Longitude grid increment (radians)
                                   ,TPH0                                &  !<-- Central latitude (rad, positive north) of projection
                                   ,TLM0                                   !<-- Central longitude (rad, positive east) of projection
!
      REAL(kind=KFPT),INTENT(OUT) :: RLAT                               &  !<-- Latitude (rad, positive north) of point on projection
                                    ,RLON                                  !<-- Longitude (rad, positive east) of point on projection
!
!-----------------------------------------------------------------------
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: I,IEND,ISTART,J,JEND,JSTART
!
      REAL(kind=KDBL) :: ARG1,ARG2,COL_MID,FCTR,GLATR,GLATD,GLOND       &
                        ,HALF,ONE,PI,R2D,ROW_MID,TLAT,TLON
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!***  Convert from transformed grid location (I,J) 
!***  to geographic coordinates (degrees).
!-----------------------------------------------------------------------
!
      ONE=1.0
      HALF=1./2.
      PI=DACOS(-ONE)
      R2D=180./PI
!
      ROW_MID=(JM+ONE)*HALF
      COL_MID=(IM+ONE)*HALF
!
!-----------------------------------------------------------------------
!
      J=J_INDEX
      I=I_INDEX
!
!-----------------------------------------------------------------------
!***  Find the rotated latitude (positive north) and 
!***  longitude (positive east).
!-----------------------------------------------------------------------
!
      TLAT=(J-ROW_MID)*DPH
      TLON=(I-COL_MID)*DLM
!
!     WRITE(0,50)I,J,TLAT*R2D,TLOND*R2D
   50 FORMAT(' I=',I4,' J=',I4,'  Projection latitude=',F8.3            &
                                       ,4X,'longitude=',F8.3)
!
!-----------------------------------------------------------------------
!***  Now convert to geographic latitude (positive north) and
!***  longitude (positive west) in degrees.
!-----------------------------------------------------------------------
!
      ARG1=DSIN(TLAT)*COS(TPH0)+DCOS(TLAT)*SIN(TPH0)*DCOS(TLON)
      RLAT=ASIN(ARG1)
!
      ARG2=DCOS(TLAT)*DCOS(TLON)/(DCOS(TLAT)*COS(TPH0))-                &
           DTAN(TLAT)*TAN(TPH0)
      IF(ABS(ARG2)>1.)ARG2=ABS(ARG2)/ARG2
      FCTR=1.
      IF(TLON>0.)FCTR=1.
      IF(TLON>PI)FCTR=-1.
!
      RLON=-TLM0+FCTR*DACOS(ARG2)
      RLON=-RLON
      IF(RLON<-PI)RLON=RLON+PI*2.
!
!     WRITE(6,100)I,J,RLAT*R2D,RLON*R2D
  100 FORMAT(' I=',I4,' J=',I4                                          &
            ,'  Geographic latitude=',F9.5,'  longitude=',F10.5)
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE REAL_IJ_TO_LATLON
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE GEO_TO_ROT_LATLON(GLAT,GLON,TPH0,TLM0                  &
                                  ,RLAT,RLON )
!
!-----------------------------------------------------------------------
!***  Convert from geographic coordinates to latitude/longitude on
!***  a rotated projection.
!-----------------------------------------------------------------------
!
      USE module_CONSTANTS,ONLY: PI
!
!------------------------
!***  Argument Variables
!------------------------
!
      REAL(kind=KFPT),INTENT(IN) :: GLAT,GLON                           &  !<-- Geographic lat/lon (rad, +east) of point
                                   ,TPH0,TLM0                              !<-- Geographic lat/lon (rad, +east) of projection center
!
      REAL(kind=KFPT),INTENT(OUT) :: RLAT,RLON                             !<-- Lat/lon (rad) of point on the projection
!
!-----------------------------------------------------------------------
!
!--------------------
!*** Local Variables
!--------------------
!
      REAL(kind=KFPT) :: X,Y,Z
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      X=COS(TPH0)*COS(GLAT)*COS(GLON-TLM0)+SIN(TPH0)*SIN(GLAT)
      Y=COS(GLAT)*SIN(GLON-TLM0)
      Z=COS(TPH0)*SIN(GLAT)-SIN(TPH0)*COS(GLAT)*COS(GLON-TLM0)
      RLAT=ATAN(Z/SQRT(X*X+Y*Y))
      RLON=ATAN(Y/X)
      IF(X<0.)THEN
        RLON=RLON+PI
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE GEO_TO_ROT_LATLON
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE DISTANCE_ON_SPHERE(RLAT_1,RLON_1                       &
                                   ,RLAT_2,RLON_2                       &
                                   ,DISTANCE )                  
!
!-----------------------------------------------------------------------
!***  Compute the great circle distance between two points on the
!***  spherical earth.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      REAL(kind=KFPT),INTENT(IN) :: RLAT_1,RLON_1                       &  !<-- Lat/lon (rad, +east) of point 1
                                   ,RLAT_2,RLON_2                          !<-- Lat/lon (rad, +east) of point 2
!
      REAL(kind=KFPT),INTENT(OUT) :: DISTANCE                              !<-- Distance (radians) between points 1 and 2
!
!-----------------------------------------------------------------------
!
!--------------------
!*** Local Variables
!--------------------
!
      REAL(kind=KDBL) :: ALPHA,ARG,BETA,CROSS,DLON,PI_H
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      PI_H=ACOS(0.)
!
!-----------------------------------------------------------------------
!
      DLON=RLON_2-RLON_1
!
      CROSS=ACOS(COS(DLON)*COS(RLAT_2))
      ARG=TAN(RLAT_2)/SIN(DLON)
      ALPHA=ATAN(ARG)
      IF(DLON<0.)ALPHA=-ALPHA
      BETA=PI_H-ALPHA
!
      DISTANCE=ACOS(COS(RLAT_1)*COS(RLAT_2)*COS(DLON)                   &
                   +SIN(RLAT_1)*SIN(CROSS)*COS(BETA))
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE DISTANCE_ON_SPHERE
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
      SUBROUTINE CENTER_NEST(SBD_DOMAIN                                 &
                            ,WBD_DOMAIN                                 &
                            ,SW_CORNER_LATD                             &
                            ,SW_CORNER_LOND                             &
                            ,TPH0D_DOMAIN                               &
                            ,TLM0D_DOMAIN )
!-----------------------------------------------------------------------
!***  Given the southern and western boundaries of a rotated lat/lon
!***  grid as well as the coordinates of the southwest corner point,
!***  find the coordinates of the grid's central point with respect
!***  to the grid upon which the rotated grid lies.
!-----------------------------------------------------------------------
!
!---------------
!***  Arguments
!---------------
!
      REAL(kind=KFPT),INTENT(IN) :: SBD_DOMAIN                          &  !<-- Latitude (deg) of domain's southern boundary
                                   ,WBD_DOMAIN                          &  !<-- Longitude (deg, +east) of domain's western boundary
                                   ,SW_CORNER_LATD                      &  !<-- Latitude (deg) of domain's southwest corner point
                                   ,SW_CORNER_LOND                         !<-- Longitude (deg, +east) of domain's southwest corner point
!
      REAL(kind=KFPT),INTENT(OUT) :: TPH0D_DOMAIN                       &  !<-- Latitude (deg) of domain's center
                                    ,TLM0D_DOMAIN                          !<-- Longitude (deg) of domain's center
!
!-----------------------------------------------------------------------
!
!---------------------
!***  Local Variables
!---------------------
!
      REAL(kind=KFPT) :: ALPHA,BETA,CENTRAL_LAT,CENTRAL_LON             &
                        ,DEG_RAD,DELTA,GAMMA                            &
                        ,PI_2,SB_R,SIDE1,SIDE2,SIDE3,SIDE4,SIDE5        &
                        ,SW_LAT,SW_LON,WB_R
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      PI_2=ACOS(0.)
      DEG_RAD=PI_2/90.
!
!-----------------------------------------------------------------------
!***  Southern and western boundaries of the rotated domain in radians.
!-----------------------------------------------------------------------
!
      SB_R=-SBD_DOMAIN*DEG_RAD
      WB_R=-WBD_DOMAIN*DEG_RAD
!
!-----------------------------------------------------------------------
!***  Southwest corner of the domain in radians.
!-----------------------------------------------------------------------
!
      SW_LAT=SW_CORNER_LATD*DEG_RAD
      SW_LON=SW_CORNER_LOND*DEG_RAD
!
!-----------------------------------------------------------------------
!***  SIDE1 is the arc from the southwest corner to the center 
!***  of the domain.
!-----------------------------------------------------------------------
!
      SIDE1=ACOS(COS(SB_R)*COS(WB_R))
!
!-----------------------------------------------------------------------
!***  ALPHA is the angle between SIDE1 and the domain's equator west of
!***  the central point.
!-----------------------------------------------------------------------
!
      ALPHA=ATAN(TAN(SB_R)/SIN(WB_R))
!
!-----------------------------------------------------------------------
!***  BETA is the angle between SIDE1 and the domain's prime meridian
!***  south of the central point.
!-----------------------------------------------------------------------
!
      BETA=PI_2-ALPHA
!
!-----------------------------------------------------------------------
!***  SIDE2 is the arc from the central point southward along the 
!***  domain's prime meridian to the great circle that intersects 
!***  both the SW and SE corners of the domain.
!-----------------------------------------------------------------------
!
      SIDE2=ATAN(COS(BETA)*TAN(SIDE1))
!
!-----------------------------------------------------------------------
!***  SIDE3 is the arc between the domain's prime meridian and the SW
!***  corner along the great circle that connects the domain's SW and
!***  SE corners.
!-----------------------------------------------------------------------
!
      SIDE3=ASIN(SIN(BETA)*SIN(SIDE1))
!
!-----------------------------------------------------------------------
!***  SIDE4 is the arc along the outer grid's equator that lies between 
!***  its western intersection with the above mentioned great circle 
!***  and the outer grid's meridian that passes through the domain's 
!***  SW corner.
!-----------------------------------------------------------------------
!
      SIDE4=ACOS(SIN(SIDE3)/COS(SW_LAT))
!
!-----------------------------------------------------------------------
!***  GAMMA is the angle between the outer grid's equator and the arc 
!***  that connects the domain's SW corner with the point where the 
!***  domain's central meridian crosses the outer grid's equator.
!-----------------------------------------------------------------------
!
      GAMMA=ATAN(TAN(SW_LAT)/COS(SIDE4))
!
!-----------------------------------------------------------------------
!***  DELTA is the angle between the arc that connects the domain's SW
!***  corner with the point where the domain's central meridian crosses
!***  the outer grid's equator and the domain's central meridian itself.
!-----------------------------------------------------------------------
!
      DELTA=PI_2-GAMMA
!
!-----------------------------------------------------------------------
!***  SIDE5 is the arc along the domain's central meridian that lies
!***  between the outer grid's equator and the great circle that passes
!***  through the SW and SE corners of the domain.
!-----------------------------------------------------------------------
!
      SIDE5=ASIN(TAN(SIDE3)/TAN(DELTA))
!
!-----------------------------------------------------------------------
!***  The central latitude and longitude of the domain in terms of
!***  the coordinates of the outer grid.
!-----------------------------------------------------------------------
!
      CENTRAL_LAT=SIDE2+SIDE5
      CENTRAL_LON=SW_LON+PI_2-SIDE4
!
      TPH0D_DOMAIN=CENTRAL_LAT/DEG_RAD
      TLM0D_DOMAIN=CENTRAL_LON/DEG_RAD
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE CENTER_NEST
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE SET_NEST_GRIDS(DOMAIN_ID_MINE                          &
                               ,TPH0D,TLM0D                             &
!!!                            ,SBD_MINE,WBD_MINE                       &
                               ,DPHD_MINE,DLMD_MINE)
! 
!-----------------------------------------------------------------------
!***  Basic grid characteristics for nests are based upon those of
!***  the uppermost parent grid.  Use those parent values to compute
!***  appropriate analogs for the nests.
!***  This subroutine is relevant only to grid-associated nests.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      INTEGER(kind=KINT),INTENT(IN) :: DOMAIN_ID_MINE                      !<-- Domain ID for this nested domain
!
      REAL(kind=KFPT),INTENT(OUT) :: DPHD_MINE                          &  !<-- Delta phi of this nested domain (degrees)
                                    ,DLMD_MINE                          &  !<-- Delta lambda of this nested domain (degrees)
                                    ,TLM0D                              &  !<-- Central rotated longitude of all domains (degrees)
                                    ,TPH0D                                 !<-- Central rotated latitude of all domains (degrees)
!!!                                 ,SBD_MINE                           &  !<-- Southern boundary this nested domain (degrees)
!!!                                 ,WBD_MINE                              !<-- Western boundary this nested domain (degrees)
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT),PARAMETER :: MAX_DOMAINS=99
!
      INTEGER(kind=KINT) :: IM_1,JM_1                                   &
                           ,ID_ANCESTOR,ID_DOMAIN                       &
                           ,IDE_1,JDE_1                                 &
                           ,I_BOUND,J_BOUND                             &
                           ,I_PARENT_SW,J_PARENT_SW                     &
                           ,I_START_SW,J_START_SW                       &
                           ,N,NUM_ANCESTORS
!
      INTEGER(kind=KINT) :: RC,RC_SET
!
      INTEGER(kind=KINT),DIMENSION(MAX_DOMAINS) :: ID_ANCESTORS=0
!
      INTEGER(kind=KINT),DIMENSION(MAX_DOMAINS) :: PARENT_CHILD_SPACE_RATIO
!
      INTEGER(kind=KINT),DIMENSION(2,MAX_DOMAINS) :: SW_CORNER
!
      REAL(kind=KFPT) :: DPHD_1,DLMD_1,TLM_BASE_1,TPH_BASE_1,SBD_1,WBD_1
      REAL(kind=KFPT) :: DPHD_X,DLMD_X,TLM_BASE,TPH_BASE,SBD_X,WBD_X
!
      CHARACTER(2)  :: INT_TO_CHAR
      CHARACTER(6)  :: FMT='(I2.2)'
      CHARACTER(50) :: GLOBAL
      CHARACTER(99) :: CONFIG_FILE_NAME
!
      TYPE(ESMF_Config),DIMENSION(MAX_DOMAINS) :: CF
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC    =ESMF_SUCCESS
      RC_SET=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  First load all of the domains' configure files.
!-----------------------------------------------------------------------
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*******  See NMM_ATM_INIT where
!*******  CF(N) is made to be
!*******  CF(ID of domain).
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      DO N=1,MAX_DOMAINS
        CF(N)=ESMF_ConfigCreate(rc=RC)
!
        WRITE(INT_TO_CHAR,FMT)N
        CONFIG_FILE_NAME='configure_file_'//INT_TO_CHAR                    !<-- Prepare the config file names
!
        CALL ESMF_ConfigLoadFile(config  =CF(N)                         &
                                ,filename=CONFIG_FILE_NAME              &
                                ,rc      =RC)
        IF(RC/=0)EXIT                                                      !<-- Exit loop after running out of config files
      ENDDO
!
!-----------------------------------------------------------------------
!***  We must loop through the configure files of all of the current
!***  domain's ancestors to collect information needed to properly
!***  describe the current grid.  This is necessary because all
!***  grids' rows and columns lie parallel to those of the uppermost
!***  grid.
!-----------------------------------------------------------------------
!
      ID_DOMAIN=DOMAIN_ID_MINE
!
      N=0
!
!-----------------------------------------------------------------------
      main_loop: DO
!-----------------------------------------------------------------------
!
        N=N+1
!
!-----------------------------
!***  Domain IDs of Ancestors
!-----------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract Domain ID of Ancestor"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=CF(ID_DOMAIN)               &  !<-- The config object
                                    ,value =ID_ANCESTOR                 &  !<-- The variable filled 
                                    ,label ='my_parent_id:'             &  !<-- Give this label's value to the previous variable
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_SET)   
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!--------------------------------------
!***  SW Corner Locations on Ancestors
!--------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Get SW Corner I and J on Ancestor Grid"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=CF(ID_DOMAIN)               &  !<-- The config object
                                    ,value =I_START_SW                  &  !<-- The variable filled 
                                    ,label ='i_parent_start:'           &  !<-- Give this label's value to the previous variable
                                    ,rc    =RC)
!
        CALL ESMF_ConfigGetAttribute(config=CF(ID_DOMAIN)               &  !<-- The config object
                                    ,value =J_START_SW                  &  !<-- The variable filled 
                                    ,label ='j_parent_start:'           &  !<-- Give this label's value to the previous variable
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_SET)   
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!----------------------------
!***  Parent-to-Child Ratios
!----------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract Child-to-Parent Ratio of Ancestor"  
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=CF(ID_DOMAIN)               &  !<-- The config object
                                    ,value =PARENT_CHILD_SPACE_RATIO(N) &  !<-- The variable filled 
                                    ,label ='parent_child_space_ratio:' &  !<-- Give this label's value to the previous variable
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_SET)   
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
        ID_ANCESTORS(N)=ID_ANCESTOR                                        !<-- Store domain IDs of all ancestors
        SW_CORNER(1,N)=I_START_SW                                          !<-- Store parent I of SW corner of its child
        SW_CORNER(2,N)=J_START_SW                                          !<-- Store parent J of SW corner of its child
!
        IF(ID_ANCESTOR==1)EXIT                                             !<-- We have reached the uppermost domain
!
        ID_DOMAIN=ID_ANCESTOR
!
!-----------------------------------------------------------------------
!
      ENDDO main_loop
!
!-----------------------------------------------------------------------
!
      NUM_ANCESTORS=N                                                        !<-- How many ancestors are there?
!
!-----------------------------------------------------------------------
!***  Rows and columns of all nests' grids lie parallel to those of 
!***  uppermost parent grid.  Thus the central rotated latitude and
!***  longitude of all nests must be those of the uppermost domain.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Central Lat/Lon of Uppermost Domain"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF(ID_ANCESTORS(N))           &  !<-- The config object
                                  ,value =TPH0D                         &  !<-- The variable filled 
                                  ,label ='tph0d:'                      &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
      CALL ESMF_ConfigGetAttribute(config=CF(ID_ANCESTORS(N))           &  !<-- The config object
                                  ,value =TLM0D                         &  !<-- The variable filled 
                                  ,label ='tlm0d:'                      &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_SET)   
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Get dimensions of uppermost domain as the baseline.
!***  We must also know southern and western boundary locations
!***  as well as whether it is global or not.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Baseline Dimensions of Uppermost Domain"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF(ID_ANCESTORS(N))           &  !<-- The config object
                                  ,value =IM_1                          &  !<-- The variable filled 
                                  ,label ='im:'                         &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
      CALL ESMF_ConfigGetAttribute(config=CF(ID_ANCESTORS(N))           &  !<-- The config object
                                  ,value =JM_1                          &  !<-- The variable filled 
                                  ,label ='jm:'                         &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_SET)   
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Southern/Western Boundary of Uppermost Domain"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF(ID_ANCESTORS(N))           &  !<-- The config object
                                  ,value =SBD_1                         &  !<-- The variable filled 
                                  ,label ='sbd:'                        &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
      CALL ESMF_ConfigGetAttribute(config=CF(ID_ANCESTORS(N))           &  !<-- The config object
                                  ,value =WBD_1                         &  !<-- The variable filled 
                                  ,label ='wbd:'                        &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_SET)   
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Global Flag for Uppermost Domain"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF(ID_ANCESTORS(N))           &  !<-- The config object
                                  ,value =GLOBAL                        &  !<-- The variable filled 
                                  ,label ='global:'                     &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_SET)   
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Full grid dimensions; delta phi and delta lambda 
!***  for uppermost domain.
!-----------------------------------------------------------------------
!
      IF(TRIM(GLOBAL)=='true')THEN                                         !<-- Uppermost domain is global
        IDE_1=IM_1+2
        JDE_1=JM_1+2
        DPHD_1=-SBD_1*2./REAL(JDE_1-3)
        DLMD_1=-WBD_1*2./REAL(IDE_1-3)
        TPH_BASE_1=SBD_1-2.*DPHD_1
        TLM_BASE_1=WBD_1-2.*DLMD_1
      ELSE                                                                 !<-- Uppermost domain is regional
        IDE_1=IM_1  
        JDE_1=JM_1  
        DPHD_1=-SBD_1*2./REAL(JDE_1-1)
        DLMD_1=-WBD_1*2./REAL(IDE_1-1)
        TPH_BASE_1=SBD_1-DPHD_1
        TLM_BASE_1=WBD_1-DLMD_1
      ENDIF
!
!-----------------------------------------------------------------------
!***  Loop through this nest's ancestors in order to obtain its:
!***  (1) delta phi and delta lambda
!***  (2) southern/western boundary locations 
!
!***  We must work downward through the ancestors because
!***  the uppermost domain is the foundation.
!-----------------------------------------------------------------------
!
      DPHD_X=DPHD_1
      DLMD_X=DLMD_1
      TPH_BASE=TPH_BASE_1
      TLM_BASE=TLM_BASE_1
!
!-----------------------------------------------------------------------
!
      work_loop: DO N=NUM_ANCESTORS,1,-1
!
        I_START_SW=SW_CORNER(1,N)
        J_START_SW=SW_CORNER(2,N)
!
        SBD_X=TPH_BASE+J_START_SW*DPHD_X                                   !<-- Southern boundary of ancestor N
        WBD_X=TLM_BASE+I_START_SW*DLMD_X                                   !<-- Western boundary of ancestor N
!
        DPHD_X=DPHD_X/REAL(PARENT_CHILD_SPACE_RATIO(N))                    !<-- Delta phi for child of ancestor N
        DLMD_X=DLMD_X/REAL(PARENT_CHILD_SPACE_RATIO(N))                    !<-- Delta lambda for child of ancestor N
!
        TPH_BASE=SBD_X-DPHD_X
        TLM_BASE=WBD_X-DLMD_X
!
      ENDDO work_loop
!
      DPHD_MINE=DPHD_X
      DLMD_MINE=DLMD_X
!!!   SBD_MINE=SBD_X
!!!   WBD_MINE=WBD_X
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE SET_NEST_GRIDS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
! 
      SUBROUTINE WATERFALLS(FIS                                         &
                           ,SEA_MASK                                    &
                           ,LOWER_TOPO                                  &
                           ,IDS,IDE,JDS,JDE)
!
!-----------------------------------------------------------------------
!***  When a parent initializes its child, the sea mask had to be done
!***  with nearest neighbor logic while FIS should be done bilinearly.
!***  This can lead to adjacent water points having different values
!***  of FIS.  when that is the case, make the elevation of all
!***  adjacent water points equal to the lowest of their values.
!***  Save the I,J of all lowered points so the atmospheric column
!***  can ultimately be adjusted.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      INTEGER(kind=KINT),INTENT(IN) :: IDS,IDE,JDS,JDE                     !<-- Lateral dimensions of nest grid
!
      REAL(kind=KFPT),DIMENSION(IDS:IDE,JDS:JDE),INTENT(IN) :: SEA_MASK    !<-- Sea mask of nest grid points
!
      REAL(kind=KFPT),DIMENSION(IDS:IDE,JDS:JDE,1),INTENT(INOUT) :: FIS    !<-- Sfc geopotential on nest grid points
!
      LOGICAL(kind=KLOG),DIMENSION(IDS:IDE,JDS:JDE),INTENT(OUT) ::      &
                                                              LOWER_TOPO   !<-- Flag points where topography is lowered
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: I,ITER,J,KOUNT_CHANGE
!
      REAL(kind=KFPT) :: FIS_0                                          &
                        ,FIS_E,FIS_N,FIS_W,FIS_S                        &
                        ,FIS_NE,FIS_NW,FIS_SW,FIS_SE                    &
                        ,FIS_NEW
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      iter_loop: DO ITER=1,500
!
        KOUNT_CHANGE=0
!
!-----------------------------------------------------------------------
!
        DO J=JDS,JDE
        DO I=IDS,IDE
!
          IF(SEA_MASK(I,J)<0.01)CYCLE                                      !<-- We are adjusting only water points 
!
!-----------------------------------------------------------------------
!
          FIS_0=FIS(I,J,1)
!
!----------
!***  East
!----------
!
          FIS_E=FIS_0
!
          IF(I+1<=IDE)THEN
            IF(SEA_MASK(I+1,J)>0.99)FIS_E=FIS(I+1,J,1)
          ENDIF
!
!---------------
!***  Northeast
!---------------
!
          FIS_NE=FIS_0
!
          IF(I+1<=IDE.AND.J+1<=JDE)THEN
            IF(SEA_MASK(I+1,J+1)>0.99)FIS_NE=FIS(I+1,J+1,1)
          ENDIF
!
!-----------
!***  North
!-----------
!
          FIS_N=FIS_0
!
          IF(J+1<=JDE)THEN
            IF(SEA_MASK(I,J+1)>0.99)FIS_N=FIS(I,J+1,1)
          ENDIF
!
!---------------
!***  Northwest
!---------------
!
          FIS_NW=FIS_0
!
          IF(I-1>=IDS.AND.J+1<=JDE)THEN
            IF(SEA_MASK(I-1,J+1)>0.99)FIS_NW=FIS(I-1,J+1,1)
          ENDIF
!
!----------
!***  West
!----------
!
          FIS_W=FIS_0
!
          IF(I-1>=IDS)THEN
            IF(SEA_MASK(I-1,J)>0.99)FIS_W=FIS(I-1,J,1)
          ENDIF
!
!---------------
!***  Southwest
!---------------
!
          FIS_SW=FIS_0
!
          IF(I-1>=IDS.AND.J-1>=JDS)THEN
            IF(SEA_MASK(I-1,J-1)>0.99)FIS_SW=FIS(I-1,J-1,1)
          ENDIF
!
!-----------
!***  South
!-----------
!
          FIS_S=FIS_0
!
          IF(J-1>=JDS)THEN
            IF(SEA_MASK(I,J-1)>0.99)FIS_S=FIS(I,J-1,1)
          ENDIF
!
!---------------
!***  Southeast
!---------------
!
          FIS_SE=FIS_0
!
          IF(I+1<=IDE.AND.J-1>=JDS)THEN
            IF(SEA_MASK(I+1,J-1)>0.99)FIS_SE=FIS(I+1,J-1,1)
          ENDIF
!
!-----------------------------------------------------------------------
!***  Lower the point in question to the lowest value of itself and 
!***  its neighbors if it is a water point.
!***  Also save all I,J locations where FIS is changed so that we
!***  can adjust the atmospheric column appropriately later.
!-----------------------------------------------------------------------
!
          FIS_NEW=MIN(FIS_0                                             &
                     ,FIS_E,FIS_N,FIS_W,FIS_E                           &
                     ,FIS_NE,FIS_NW,FIS_SW,FIS_SE)
!
          IF(FIS_NEW+0.1<FIS_0)THEN
            KOUNT_CHANGE=KOUNT_CHANGE+1
            FIS(I,J,1)=FIS_NEW
            LOWER_TOPO(I,J)=.TRUE.
          ENDIF
!
!-----------------------------------------------------------------------
!
        ENDDO
        ENDDO
!
!
        IF(KOUNT_CHANGE==0)EXIT iter_loop
        IF(ITER==100)THEN
          WRITE(0,*)' Reached 100 iterations and KOUNT_CHANGE='         &
                   ,KOUNT_CHANGE
        ENDIF
!
      ENDDO iter_loop
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE WATERFALLS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
! 
      SUBROUTINE ADJUST_COLUMNS(PD_NEAREST                              &
                               ,PD_BILINEAR                             &
                               ,LOWER_TOPO                              &
                               ,DUMMY_3D                                &
                               ,PT                                      &
                               ,PDTOP                                   &
                               ,SG1                                     &
                               ,SG2                                     &
                               ,IM_CHILD                                &
                               ,JM_CHILD                                &
                                         )
!
!-----------------------------------------------------------------------
!***  When the surface elevation of a nested domain is changed due to
!***  leveling of adjacent water points, adjust the atmospheric column
!***  at each of those points.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      INTEGER(kind=KINT),INTENT(IN) :: IM_CHILD                         &
                                      ,JM_CHILD
!
      REAL(kind=KFPT),INTENT(IN) :: PDTOP,PT
!
      REAL(kind=KFPT),DIMENSION(1:LM+1),INTENT(IN) :: SG1,SG2
!
      REAL(kind=KFPT),DIMENSION(1:IM_CHILD,1:JM_CHILD,1),INTENT(IN) ::  &
                                                            PD_NEAREST  &
                                                           ,PD_BILINEAR
!
      LOGICAL(kind=KLOG),DIMENSION(1:IM_CHILD,1:JM_CHILD),INTENT(IN) :: &
                                                            LOWER_TOPO
!
      REAL(kind=KFPT),DIMENSION(1:IM_CHILD,1:JM_CHILD,LM),INTENT(INOUT) :: &
                                                                  DUMMY_3D
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: I,J,L,NUM_LEVS_SEC,NUM_LEVS_SPLINE
!
      REAL(kind=KFPT) :: COEFF_1,DELP_EXTRAP,PBOT_IN,PBOT_TARGET        &
                        ,PDTOP_PT,PTOP_IN,PTOP_TARGET,R_DELP
!
      REAL(kind=KFPT),DIMENSION(1:LM) :: PMID_TARGET                    &
                                        ,VBL_COLUMN
!
      REAL(kind=KFPT),DIMENSION(1:LM+1) :: PMID_IN                      &
                                          ,SEC_DERIV                    &
                                          ,VBL_INPUT
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      NUM_LEVS_SPLINE=LM+1
      NUM_LEVS_SEC   =LM+1
!
!-----------------------------------------------------------------------
!***  Compute the input and target mid-layer pressures for the
!***  spline interpolation.
!-----------------------------------------------------------------------
!
      DO J=1,JM_CHILD
      DO I=1,IM_CHILD
!
!-----------------------------------------------------------------------
!
        adjust: IF(LOWER_TOPO(I,J))THEN
!
          PTOP_IN=SG2(1)*PD_BILINEAR(I,J,1)+SG1(1)*PDTOP+PT
          PTOP_TARGET=SG2(1)*PD_NEAREST(I,J,1)+SG1(1)*PDTOP+PT
!
          DO L=1,LM+1
            SEC_DERIV(L)=0.
          ENDDO
!
          DO L=1,LM
!
            PDTOP_PT=SG1(L+1)*PDTOP+PT
!
            PBOT_IN   =SG2(L+1)*PD_BILINEAR(I,J,1)+PDTOP_PT
            PMID_IN(L)=0.5*(PTOP_IN+PBOT_IN)
            PTOP_IN   =PBOT_IN
            VBL_INPUT(L)=DUMMY_3D(I,J,L)
!
            PBOT_TARGET   =SG2(L+1)*PD_NEAREST(I,J,1)+PDTOP_PT
            PMID_TARGET(L)=0.5*(PTOP_TARGET+PBOT_TARGET)
            PTOP_TARGET   =PBOT_TARGET
!
          ENDDO
!
!***  We know the target mid-layer pressure is greater than that in
!***  the original column since the sfc elevation has been lowered.
!***  Add a new input level by extrapolating linearly in pressure
!***  to obtain a value at the lowest output mid-layer then fill
!***  in all output levels with the spline.
!
          PMID_IN(LM+1)=PMID_TARGET(LM)
          R_DELP=1./(PMID_IN(LM)-PMID_IN(LM-1))
          DELP_EXTRAP=PMID_TARGET(LM)-PMID_IN(LM)
          COEFF_1=(VBL_INPUT(LM)-VBL_INPUT(LM-1))*R_DELP
          VBL_INPUT(LM+1)=VBL_INPUT(LM)+COEFF_1*DELP_EXTRAP                !<-- Create extrapolated value at nest's lowest mid-layer
!                                                                               in input array.
!-----------------------------------------------------------------------
!
          IF(ABS(PMID_IN(LM+1)-PMID_IN(LM))<10.)EXIT
!
          CALL SPLINE(NUM_LEVS_SPLINE                                   &  !<-- # of input levels
                     ,PMID_IN                                           &  !<-- Input mid-layer pressures
                     ,VBL_INPUT                                         &  !<-- Input mid-layer mass variable value
                     ,SEC_DERIV                                         &  !<-- Specified 2nd derivatives (=0) at input levels
                     ,NUM_LEVS_SEC                                      &  !<-- Vertical dimension of SEC_DERIV
                     ,LM                                                &  !<-- # of mid-layers to which to interpolate
                     ,PMID_TARGET                                       &  !<-- Mid-layer pressures to which to interpolate 
                     ,VBL_COLUMN )                                         !<-- Mid-layer variable value returned
!
          DO L=1,LM
            DUMMY_3D(I,J,L)=VBL_COLUMN(L)
          ENDDO
!
        ENDIF adjust
!
!-----------------------------------------------------------------------
!
      ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE ADJUST_COLUMNS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE INTERNAL_DATA_TO_DOMAIN(EXP_STATE_SOLVER               &
                                        ,EXP_STATE_DOMAIN               &
                                        ,LM )
!
!-----------------------------------------------------------------------
!***  Transfer from Solver export state to the DOMAIN export state
!***  the data needed for parent generation of child boundary data
!***  as well as internal updates and shift decisions in moving
!***  nests.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      INTEGER(kind=KINT),INTENT(OUT) :: LM                                 !<-- # of model layers
!
      TYPE(ESMF_State),INTENT(INOUT) :: EXP_STATE_SOLVER                   !<-- Solver export state
!
      TYPE(ESMF_State),INTENT(INOUT) :: EXP_STATE_DOMAIN                   !<-- DOMAIN export state into which fcst Arrays are transferred
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: ITS,ITE,JTS,JTE                             &
                           ,IDS,IDE,JDS,JDE                             &
                           ,J,JM                                        &
                           ,INDX_CW,INDX_Q                              &
                           ,LMP1,LNSH,LNSV                              &
                           ,N,NHALO,NKOUNT
!
      INTEGER(kind=KINT) :: RC,RC_TRANS
!
      REAL(kind=KFPT) :: DLMD,DPHD,DYH,PDTOP,PT
!
      REAL(kind=KFPT),DIMENSION(:),ALLOCATABLE :: ARRAY_1D
!
      TYPE(ESMF_StateItem_Flag) :: ITEMTYPE
!
      TYPE(ESMF_Field) :: HOLD_FIELD
!
      LOGICAL(kind=KLOG) :: RESTART
!
      INTEGER(kind=KINT) :: itemCount
!
      CHARACTER(LEN=14),DIMENSION(:),ALLOCATABLE :: itemnamelist
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC_TRANS=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  First find out the width of the haloes which is needed in
!***  renaming the ESMF Fields associated with nest boundary updates.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract NHALO from the Solver export state"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=EXP_STATE_SOLVER                     &  !<-- The Solver export state
                            ,name ='NHALO'                              &  !<-- Name of Attribute to extract
                            ,value=NHALO                                &  !<-- Put the extracted Attribute here
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_TRANS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateGet(state       =EXP_STATE_SOLVER                  &  !<-- The Solver export state
                        ,itemCount   =itemCount                         &  !<-- # of items in the state
                        ,rc          =RC)
!
      ALLOCATE(itemNameList(itemCount))
!
      CALL ESMF_StateGet(state       =EXP_STATE_SOLVER                  &  !<-- The Solver export state
                        ,itemNameList=itemNameList                      &  !<-- List of item names
                        ,rc          =RC)
!
!-----------------------------------------------------------------------
!
      item_loop: DO N=1,SIZE(itemNameList)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Test the presence of "//TRIM(itemNameList(N))//" in Solver Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_StateGet(state       =EXP_STATE_SOLVER                &  !<-- The Solver export state
                          ,itemName    =TRIM(itemNameList(N))           &  !<-- Check presence of this Field
                          ,itemType    =ITEMTYPE                        &  !<-- ESMF Type of the Field
                          ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_TRANS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        IF (ITEMTYPE == ESMF_STATEITEM_FIELD) THEN

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Extract "//TRIM(itemNameList(N))//" field from Solver Export State"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_StateGet(state   =EXP_STATE_SOLVER                  &  !<-- The Solver export state
                            ,itemName=TRIM(itemNameList(N))             &  !<-- Extract this Field
                            ,field   =HOLD_FIELD                        &  !<-- Put the extracted Field here
                            ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_TRANS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Insert "//TRIM(itemNameList(N))//" field into DOMAIN Export State"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_StateAddReplace(EXP_STATE_DOMAIN                  &  !<-- Insert field into DOMAIN export state
                                   ,(/HOLD_FIELD/)                    &  !<-- The Field to be inserted
                                   ,rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_TRANS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        END IF
!
      END DO  item_loop
!
!-----------------------------------------------------------------------
!
!---------------------
!***  Transfer INDX_Q
!---------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract INDX_Q from Solver Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=EXP_STATE_SOLVER                     &  !<-- The Solver export state
                            ,name ='INDX_Q'                             &  !<-- Name of Attribute to extract
                            ,value=INDX_Q                               &  !<-- Put the extracted Attribute here
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_TRANS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Insert INDX_Q into DOMAIN Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(state=EXP_STATE_DOMAIN                     &  !<-- The DOMAIN export state
                            ,name ='INDX_Q'                             &  !<-- The name of the Attribute to insert
                            ,value=INDX_Q                               &  !<-- The Attribute to be inserted
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_TRANS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!----------------------
!***  Transfer INDX_CW
!----------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract INDX_CW from Solver Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=EXP_STATE_SOLVER                     &  !<-- The Solver export state
                            ,name ='INDX_CW'                            &  !<-- Name of Attribute to extract
                            ,value=INDX_CW                              &  !<-- Put the extracted Attribute here
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_TRANS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Insert INDX_CW into DOMAIN Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(state=EXP_STATE_DOMAIN                     &  !<-- The DOMAIN export state
                            ,name ='INDX_CW'                            &  !<-- The name of the Attribute to insert
                            ,value=INDX_CW                              &  !<-- The Attribute to be inserted
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_TRANS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!--------------------------------------
!***  LM is needed but not transferred
!--------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract LM from DYN Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=EXP_STATE_SOLVER                     &  !<-- The Solver export state
                            ,name ='LM'                                 &  !<-- Name of Attribute to extract
                            ,value=LM                                   &  !<-- Put the extracted Attribute here
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_TRANS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-------------------------------------------------------
!***  Transfer PT,PDTOP,PSGML1,SG1,SG2,SGML2,DSG2,PDSG1
!-------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Insert PT into DOMAIN Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=EXP_STATE_SOLVER                     &  !<-- The Solver export state
                            ,name ='PT'                                 &  !<-- Extract PT
                            ,value=PT                                   &  !<-- Put the extracted Attribute here
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=EXP_STATE_DOMAIN                     &  !<-- The DOMAIN export state
                            ,name ='PT'                                 &  !<-- Set PT
                            ,value=PT                                   &  !<-- The Attribute to be inserted
                            ,rc   =RC)
!
      CALL ESMF_AttributeGet(state=EXP_STATE_SOLVER                     &  !<-- The Solver export state
                            ,name ='PDTOP'                              &  !<-- Extract PDTOP
                            ,value=PDTOP                                &  !<-- Put the extracted Attribute here
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=EXP_STATE_DOMAIN                     &  !<-- The Parent's DOMAIN export state
                            ,name ='PDTOP'                              &  !<-- Set PDTOP
                            ,value=PDTOP                                &  !<-- Put the extracted Attribute here
                            ,rc   =RC)
!
      ALLOCATE(ARRAY_1D(1:LM))
!
      CALL ESMF_AttributeGet(state    =EXP_STATE_SOLVER                 &  !<-- The Solver export state
                            ,name     ='PSGML1'                         &  !<-- Extract PGMSL1
                            ,itemCount=LM                               &  !<-- # of words in data list
                            ,valueList=ARRAY_1D                         &  !<-- Put extracted values here
                            ,rc       =RC)
!
      CALL ESMF_AttributeSet(state    =EXP_STATE_DOMAIN                 &  !<-- The Parent's DOMAIN export state
                            ,name     ='PSGML1'                         &  !<-- Extract PGMSL1
                            ,itemCount=LM                               &  !<-- # of words in data list
                            ,valueList=ARRAY_1D                         &  !<-- Put extracted values here
                            ,rc       =RC)
!
      CALL ESMF_AttributeGet(state    =EXP_STATE_SOLVER                 &  !<-- The Solver export state
                            ,name     ='SGML2'                          &  !<-- Extract SGML2
                            ,itemCount=LM                               &  !<-- # of words in data list
                            ,valueList=ARRAY_1D                         &  !<-- Put extracted values here
                            ,rc       =RC)
!
      CALL ESMF_AttributeSet(state    =EXP_STATE_DOMAIN                 &  !<-- The Parent's DOMAIN export state
                            ,name     ='SGML2'                          &  !<-- Set SGML2
                            ,itemCount=LM                               &  !<-- # of words in data list
                            ,valueList=ARRAY_1D                         &  !<-- Put extracted values here
                            ,rc       =RC)
!
      CALL ESMF_AttributeGet(state    =EXP_STATE_SOLVER                 &  !<-- The Solver export state
                            ,name     ='DSG2'                           &  !<-- Extract DSG2
                            ,itemCount=LM                               &  !<-- # of words in data list
                            ,valueList=ARRAY_1D                         &  !<-- Put extracted values here
                            ,rc       =RC)
!
      CALL ESMF_AttributeSet(state    =EXP_STATE_DOMAIN                 &  !<-- The Parent's DOMAIN export state
                            ,name     ='DSG2'                           &  !<-- Set DSG2
                            ,itemCount=LM                               &  !<-- # of words in data list
                            ,valueList=ARRAY_1D                         &  !<-- Put extracted values here
                            ,rc       =RC)
!
      CALL ESMF_AttributeGet(state    =EXP_STATE_SOLVER                 &  !<-- The Solver export state
                            ,name     ='PDSG1'                          &  !<-- Extract PDSG1
                            ,itemCount=LM                               &  !<-- # of words in data list
                            ,valueList=ARRAY_1D                         &  !<-- Put extracted values here
                            ,rc       =RC)
!
      CALL ESMF_AttributeSet(state    =EXP_STATE_DOMAIN                 &  !<-- The Parent's DOMAIN export state
                            ,name     ='PDSG1'                          &  !<-- Set PDSG1
                            ,itemCount=LM                               &  !<-- # of words in data list
                            ,valueList=ARRAY_1D                         &  !<-- Put extracted values here
                            ,rc       =RC)
!
      LMP1=LM+1
      DEALLOCATE(ARRAY_1D)
      ALLOCATE(ARRAY_1D(1:LMP1))
!
      CALL ESMF_AttributeGet(state    =EXP_STATE_SOLVER                 &  !<-- The Solver export state
                            ,name     ='SG1'                            &  !<-- Extract PGMSL1
                            ,itemCount=LMP1                             &  !<-- # of words in data list
                            ,valueList=ARRAY_1D                         &  !<-- Put extracted values here
                            ,rc       =RC)
!
      CALL ESMF_AttributeSet(state    =EXP_STATE_DOMAIN                 &  !<-- The Parent's DOMAIN export state
                            ,name     ='SG1'                            &  !<-- Extract PGMSL1
                            ,itemCount=LMP1                             &  !<-- # of words in data list
                            ,valueList=ARRAY_1D                         &  !<-- Put extracted values here
                            ,rc       =RC)
!
      CALL ESMF_AttributeGet(state    =EXP_STATE_SOLVER                 &  !<-- The Solver export state
                            ,name     ='SG2'                            &  !<-- Extract PGMSL1
                            ,itemCount=LMP1                             &  !<-- # of words in data list
                            ,valueList=ARRAY_1D                         &  !<-- Put extracted values here
                            ,rc       =RC)
!
      CALL ESMF_AttributeSet(state    =EXP_STATE_DOMAIN                 &  !<-- The Parent's DOMAIN export state
                            ,name     ='SG2'                            &  !<-- Extract PGMSL1
                            ,itemCount=LMP1                             &  !<-- # of words in data list
                            ,valueList=ARRAY_1D                         &  !<-- Put extracted values here
                            ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_TRANS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      DEALLOCATE(ARRAY_1D)
!
!-------------------------------------------
!***  Transfer Subdomain Integration Limits
!***  and Full Domain Limits
!-------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Integration Limits from Solver Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=EXP_STATE_SOLVER                     &  !<-- The Solver export state
                            ,name ='ITS'                                &  !<-- Name of Attribute to extract
                            ,value=ITS                                  &  !<-- Put the extracted Attribute here
                            ,rc   =RC)
!
      CALL ESMF_AttributeGet(state=EXP_STATE_SOLVER                     &  !<-- The Solver export state
                            ,name ='ITE'                                &  !<-- Name of Attribute to extract
                            ,value=ITE                                  &  !<-- Put the extracted Attribute here
                            ,rc   =RC)
!
      CALL ESMF_AttributeGet(state=EXP_STATE_SOLVER                     &  !<-- The Solver export state
                            ,name ='JTS'                                &  !<-- Name of Attribute to extract
                            ,value=JTS                                  &  !<-- Put the extracted Attribute here
                            ,rc   =RC)
!
      CALL ESMF_AttributeGet(state=EXP_STATE_SOLVER                     &  !<-- The Solver export state
                            ,name ='JTE'                                &  !<-- Name of Attribute to extract
                            ,value=JTE                                  &  !<-- Put the extracted Attribute here
                            ,rc   =RC)
!
      CALL ESMF_AttributeGet(state=EXP_STATE_SOLVER                     &  !<-- The Solver export state
                            ,name ='LM'                                 &  !<-- Name of Attribute to extract
                            ,value=LM                                   &  !<-- Put the extracted Attribute here
                            ,rc   =RC)
!
      CALL ESMF_AttributeGet(state=EXP_STATE_SOLVER                     &  !<-- The Solver export state
                            ,name ='NHALO'                              &  !<-- Name of Attribute to extract
                            ,value=NHALO                                &  !<-- Put the extracted Attribute here
                            ,rc   =RC)
!
      CALL ESMF_AttributeGet(state=EXP_STATE_SOLVER                     &  !<-- The Solver export state
                            ,name ='IDS'                                &  !<-- Name of Attribute to extract
                            ,value=IDS                                  &  !<-- Put the extracted Attribute here
                            ,rc   =RC)
!
      CALL ESMF_AttributeGet(state=EXP_STATE_SOLVER                     &  !<-- The Solver export state
                            ,name ='IDE'                                &  !<-- Name of Attribute to extract
                            ,value=IDE                                  &  !<-- Put the extracted Attribute here
                            ,rc   =RC)
!
      CALL ESMF_AttributeGet(state=EXP_STATE_SOLVER                     &  !<-- The Solver export state
                            ,name ='JDS'                                &  !<-- Name of Attribute to extract
                            ,value=JDS                                  &  !<-- Put the extracted Attribute here
                            ,rc   =RC)
!
      CALL ESMF_AttributeGet(state=EXP_STATE_SOLVER                     &  !<-- The Solver export state
                            ,name ='JDE'                                &  !<-- Name of Attribute to extract
                            ,value=JDE                                  &  !<-- Put the extracted Attribute here
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_TRANS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Insert Integration Limits into DOMAIN Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(state=EXP_STATE_DOMAIN                     &  !<-- The DOMAIN export state
                            ,name ='ITS'                                &  !<-- The name of the Attribute to insert
                            ,value=ITS                                  &  !<-- The Attribute to be inserted
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=EXP_STATE_DOMAIN                     &  !<-- The DOMAIN export state
                            ,name ='ITE'                                &  !<-- The name of the Attribute to insert
                            ,value=ITE                                  &  !<-- The Attribute to be inserted
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=EXP_STATE_DOMAIN                     &  !<-- The DOMAIN export state
                            ,name ='JTS'                                &  !<-- The name of the Attribute to insert
                            ,value=JTS                                  &  !<-- The Attribute to be inserted
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=EXP_STATE_DOMAIN                     &  !<-- The DOMAIN export state
                            ,name ='JTE'                                &  !<-- The name of the Attribute to insert
                            ,value=JTE                                  &  !<-- The Attribute to be inserted
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=EXP_STATE_DOMAIN                     &  !<-- The DOMAIN export state
                            ,name ='LM'                                 &  !<-- The name of the Attribute to insert
                            ,value=LM                                   &  !<-- The Attribute to be inserted
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=EXP_STATE_DOMAIN                     &  !<-- The DOMAIN export state
                            ,name ='NHALO'                              &  !<-- The name of the Attribute to insert
                            ,value=NHALO                                &  !<-- The Attribute to be inserted
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=EXP_STATE_DOMAIN                     &  !<-- The DOMAIN export state
                            ,name ='IDS'                                &  !<-- The name of the Attribute to insert
                            ,value=IDS                                  &  !<-- The Attribute to be inserted
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=EXP_STATE_DOMAIN                     &  !<-- The DOMAIN export state
                            ,name ='IDE'                                &  !<-- The name of the Attribute to insert
                            ,value=IDE                                  &  !<-- The Attribute to be inserted
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=EXP_STATE_DOMAIN                     &  !<-- The DOMAIN export state
                            ,name ='JDS'                                &  !<-- The name of the Attribute to insert
                            ,value=JDS                                  &  !<-- The Attribute to be inserted
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=EXP_STATE_DOMAIN                     &  !<-- The DOMAIN export state
                            ,name ='JDE'                                &  !<-- The name of the Attribute to insert
                            ,value=JDE                                  &  !<-- The Attribute to be inserted
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_TRANS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!--------------------------------------------------------
!***  Transfer Width of Blending Region Along Boundaries
!--------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract LNSH,LNSV from Solver Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=EXP_STATE_SOLVER                     &  !<-- The Solver export state
                            ,name ='LNSH'                               &  !<-- Name of Attribute to extract
                            ,value=LNSH                                 &  !<-- Put the extracted Attribute here
                            ,rc   =RC)
!
      CALL ESMF_AttributeGet(state=EXP_STATE_SOLVER                     &  !<-- The Solver export state
                            ,name ='LNSV'                               &  !<-- Name of Attribute to extract
                            ,value=LNSV                                 &  !<-- Put the extracted Attribute here
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_TRANS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Insert LNSH,LNSV into DOMAIN Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(state=EXP_STATE_DOMAIN                     &  !<-- The DOMAIN export state
                            ,name ='LNSH'                               &  !<-- The name of the Attribute to insert
                            ,value=LNSH                                 &  !<-- The Attribute to be inserted
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=EXP_STATE_DOMAIN                     &  !<-- The DOMAIN export state
                            ,name ='LNSV'                               &  !<-- The name of the Attribute to insert
                            ,value=LNSV                                 &  !<-- The Attribute to be inserted
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_TRANS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-------------------------------
!***  Transfer the restart flag
!-------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Restart Flag from Solver Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=EXP_STATE_SOLVER                     &  !<-- The Solver export state
                            ,name ='RESTART'                            &  !<-- Name of Attribute to extract
                            ,value=RESTART                              &  !<-- Put the extracted Attribute here
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=EXP_STATE_DOMAIN                     &  !<-- The DOMAIN export state
                            ,name ='RESTART'                            &  !<-- The name of the Attribute to insert
                            ,value=RESTART                              &  !<-- The Attribute to be inserted
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_TRANS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!------------------------
!***  Transfer DX and DY
!------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract DYH from Solver Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=EXP_STATE_SOLVER                     &  !<-- The Solver export state
                            ,name ='DYH'                                &  !<-- Name of DYH scalar
                            ,value=DYH                                  &  !<-- Put the extracted Attribute here
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=EXP_STATE_DOMAIN                     &  !<-- The DOMAIN export state
                            ,name ='DYH'                                &  !<-- Name of DYH scalar
                            ,value=DYH                                  &  !<-- Put the extracted Attribute here
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_TRANS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      NKOUNT=JDE-JDS+1
      ALLOCATE(ARRAY_1D(JDS:JDE))
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract DXH from Solver Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state    =EXP_STATE_SOLVER                 &  !<-- The Solver export state
                            ,name     ='DXH'                            &  !<-- Name of DXH array
                            ,itemCount=NKOUNT                           &  !<-- # of words in data list
                            ,valueList=ARRAY_1D                         &  !<-- Put extracted values here
                            ,rc       =RC)
!
      CALL ESMF_AttributeSet(state    =EXP_STATE_DOMAIN                 &  !<-- The DOMAIN export state
                            ,name     ='DXH'                            &  !<-- Name of DXH array
                            ,itemCount=NKOUNT                           &  !<-- # of words in data list
                            ,valueList=ARRAY_1D                         &  !<-- Put extracted values here
                            ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_TRANS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      DEALLOCATE(ARRAY_1D)
!
!----------------------------
!***  Transfer DPHD,DLMD,,JM
!----------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract DPHD,DLMD,JM from Solver Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=EXP_STATE_SOLVER                     &  !<-- The Solver export state
                            ,name ='DPHD'                               &  !<-- Latitude grid increment (deg)
                            ,value=DPHD                                 &  !<-- Put the extracted Attribute here
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=EXP_STATE_DOMAIN                     &  !<-- The DOMAIN export state
                            ,name ='DPHD'                               &  !<-- Latitude grid increment (deg)
                            ,value=DPHD                                 &  !<-- Put the extracted Attribute here
                            ,rc   =RC)
!
      CALL ESMF_AttributeGet(state=EXP_STATE_SOLVER                     &  !<-- The Solver export state
                            ,name ='DLMD'                               &  !<-- Longitude grid increment (deg)
                            ,value=DLMD                                 &  !<-- Put the extracted Attribute here
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=EXP_STATE_DOMAIN                     &  !<-- The DOMAIN export state
                            ,name ='DLMD'                               &  !<-- Longitude grid increment (deg)
                            ,value=DLMD                                 &  !<-- Put the extracted Attribute here
                            ,rc   =RC)
!
      CALL ESMF_AttributeGet(state=EXP_STATE_SOLVER                     &  !<-- The Solver export state
                            ,name ='JM'                                 &  !<-- J index extent of domain
                            ,value=JM                                   &  !<-- Put the extracted Attribute here
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=EXP_STATE_DOMAIN                     &  !<-- The DOMAIN export state
                            ,name ='JM'                                 &  !<-- J index extent of domain
                            ,value=JM                                   &  !<-- Put the extracted Attribute here
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_TRANS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE INTERNAL_DATA_TO_DOMAIN
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE BOUNDARY_DATA_STATE_TO_STATE(S_BDY,N_BDY,W_BDY,E_BDY   &
                                             ,CLOCK                     &
                                             ,PARENT                    &
                                             ,NEST                      &
                                             ,RATIO                     &
                                             ,STATE_IN                  &
                                             ,STATE_OUT )
!
!-----------------------------------------------------------------------
!***  This routine moves new boundary data for nested domains from
!***  one import/export state to another.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      LOGICAL(kind=KLOG),INTENT(IN),OPTIONAL :: S_BDY,N_BDY,W_BDY,E_BDY    !<-- Is this task on the domain's boundary?
!
      TYPE(ESMF_Clock),INTENT(IN),OPTIONAL :: CLOCK                        !<-- ESMF Clock
!
      INTEGER(kind=KINT),INTENT(IN),OPTIONAL :: RATIO                      !<-- # of child timesteps per parent timestep          
!
      LOGICAL(kind=KLOG),INTENT(IN),OPTIONAL :: NEST                    &  !<-- Is task on a nest domain?
                                               ,PARENT                     !<-- Is task on a parent domain?
      TYPE(ESMF_State),INTENT(INOUT) :: STATE_IN                        &  !<-- Input ESMF State
                                       ,STATE_OUT                          !<-- Output ESMF State
!
!---------------------
!***  Local Variables
!---------------------
!
      TYPE SIDES_1D_REAL
        REAL(kind=KFPT),DIMENSION(:),ALLOCATABLE :: SOUTH
        REAL(kind=KFPT),DIMENSION(:),ALLOCATABLE :: NORTH
        REAL(kind=KFPT),DIMENSION(:),ALLOCATABLE :: WEST
        REAL(kind=KFPT),DIMENSION(:),ALLOCATABLE :: EAST
      END TYPE SIDES_1D_REAL
!
      INTEGER(kind=KINT) :: I_SHIFT,I_SW_PARENT_NEW,ISTAT               &
                           ,J_SHIFT,J_SW_PARENT_NEW,LAST_STEP_MOVED     &
                           ,KOUNT,LIMIT                                 &
                           ,NEXT_MOVE_TIMESTEP,NTIMESTEP,NTYPE          &
                           ,RC,RC_BND_MV
!
      INTEGER(kind=ESMF_KIND_I8) :: NTIMESTEP_ESMF
!
      INTEGER(kind=KINT),DIMENSION(1:NUM_DOMAINS_MAX) :: NTIMESTEP_CHILD_MOVES
!
      CHARACTER(len=7) :: TIME
!
      LOGICAL(kind=KLOG),SAVE :: EXTRACTED_FLAGS=.FALSE.
!
!!!   LOGICAL(kind=KLOG),SAVE :: I_AM_A_NEST                            &
      LOGICAL(kind=KLOG) :: I_AM_A_NEST                                 &
                           ,MY_DOMAIN_MOVES
!
      LOGICAL(kind=KLOG) :: MOVE_NOW
!
      TYPE(SIDES_1D_REAL),SAVE :: BOUNDARY_H                            &
                                 ,BOUNDARY_V
!
      integer :: nnnn
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  If the Clock was sent in then this transfer of data is 
!***  dependent upon the timestep.  
!-----------------------------------------------------------------------
!
      IF(PRESENT(CLOCK))THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="State_to_State Gets DOMAIN Clock"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ClockGet(clock       =CLOCK                           &  !<-- The ESMF clock
                          ,advanceCount=NTIMESTEP_ESMF                  &  !<-- The number of times the clock has been advanced
                          ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        NTIMESTEP=NTIMESTEP_ESMF
!
        IF(MOD(NTIMESTEP,RATIO)/=0)THEN
          RETURN                                                           !<-- There is new bndry data only at parent timesteps
        ENDIF
!
!     else
!       write(0,*)' enter BOUNDARY_DATA_STATE_TO_STATE clock is NOT present'
      ENDIF
!
!-----------------------------------------------------------------------
!***  Extract only once the Nest flag and the Moving Nest flag since
!***  they never change.
!-----------------------------------------------------------------------
!
!!!   IF(.NOT.EXTRACTED_FLAGS)THEN
!
        EXTRACTED_FLAGS=.TRUE.
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Check the Nest Flag in BOUNDARY_DATA_STATE_TO_STATE"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state=STATE_IN                           &  !<-- Look at the input state
                              ,name ='I-Am-A-Nest Flag'                 &  !<-- Extract Attribute with this name
                              ,value=I_AM_A_NEST                        &  !<-- Is this domain a nest?
                              ,rc   =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert the Nest Flag in BOUNDARY_DATA_STATE_TO_STATE"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=STATE_OUT                          &  !<-- Look at the output state
                              ,name ='I-Am-A-Nest Flag'                 &  !<-- Insert Attribute with this name
                              ,value=I_AM_A_NEST                        &  !<-- Is this domain a nest?
                              ,rc   =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        IF(I_AM_A_NEST)THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Extract the Moving Nest Flag in BOUNDARY_DATA_STATE_TO_STATE"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeGet(state=STATE_IN                         &  !<-- Look at the input state
                                ,name ='MY_DOMAIN_MOVES'                &  !<-- Extract Attribute with this name
                                ,value=MY_DOMAIN_MOVES                  &  !<-- Does the nest move?
                                ,rc   =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Insert the Moving Nest Flag in BOUNDARY_DATA_STATE_TO_STATE"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeSet(state=STATE_OUT                        &  !<-- Look at the output state
                                ,name ='MY_DOMAIN_MOVES'                &  !<-- Insert Attribute with this name
                                ,value=MY_DOMAIN_MOVES                  &  !<-- Does the nest move?
                                ,rc   =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        ENDIF
!
!!!   ENDIF
!
!-----------------------------------------------------------------------
!***  Parents transfer the next move timesteps of their moving
!***  children.
!-----------------------------------------------------------------------
!
      IF(PRESENT(PARENT))THEN
!
        IF(PARENT)THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract NTIMESTEP_CHILD_MOVES in BOUNDARY_DATA_STATE_TO_STATE"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeGet(state    =STATE_IN                     &  !<-- Look at the input state
                                ,name     ='NEXT_TIMESTEP_CHILD_MOVES'  &  !<-- Get this Attribute
                                ,valueList=NTIMESTEP_CHILD_MOVES        &  !<-- What are the children's next move timesteps?
                                ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Insert NTIMESTEP_CHILD_MOVES into the Output State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeSet(state    =STATE_OUT                    &  !<-- Insert data into output state
                                ,name     ='NEXT_TIMESTEP_CHILD_MOVES'  &  !<-- Get this Attribute
                                ,itemCount=NUM_DOMAINS_MAX              &  !<-- How many items?
                                ,valueList=NTIMESTEP_CHILD_MOVES        &  !<-- What are the children's next move timesteps?
                                ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        ENDIF
!
      ENDIF
!
!-----------------------------------------------------------------------
!***  There are two 'types' of nest boundary data that must be
!***  considered.  The first is the standard data coming from the
!***  parent at the end of each parent timestep that is sent back
!***  to all of the children at the beginning of that timestep so
!***  the children can form time tendencies of their boundary variables.
!***  The second is new boundary data from the parent generated at the
!***  end of a timestep at which it has received a move signal from
!***  a child.  Simply put, the first data is X(N+1) and the second is
!***  X(N) used in the children's computation of boundary tendencies
!***  where dX/dt=[X(N+1)-X(N)]/DT(parent).  If the MOVE_NOW flag from
!***  STATE_IN is false then only the first type of boundary data is
!***  present and needs to be moved into STATE_OUT but if MOVE_NOW is 
!***  true then both types are present and must be transfered.
!-----------------------------------------------------------------------
!
      MOVE_NOW=.FALSE.
!
      IF(I_AM_A_NEST.AND.MY_DOMAIN_MOVES.AND..NOT.PRESENT(PARENT))THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract NEXT_MOVE_TIMESTEP in BOUNDARY_DATA_STATE_TO_STATE"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state=STATE_IN                           &  !<-- Look at the input state
                              ,name ='NEXT_MOVE_TIMESTEP'               &  !<-- Is this name present?
                              ,value=NEXT_MOVE_TIMESTEP                 &  !<-- At which timestep does the nest move?
                              ,rc   =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert NEXT_MOVE_TIMESTEP into the Output State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=STATE_OUT                          &  !<-- Insert data into output state
                              ,name ='NEXT_MOVE_TIMESTEP'               &  !<-- The name of the data 
                              ,value=NEXT_MOVE_TIMESTEP                 &  !<-- At which timestep does the nest move?
                              ,rc   =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract the MOVE_NOW Flag in BOUNDARY_DATA_STATE_TO_STATE"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state=STATE_IN                           &  !<-- Look at the input state
                              ,name ='MOVE_NOW'                         &  !<-- Is this name present?
                              ,value=MOVE_NOW                           &  !<-- Is the child moving right now?
                              ,rc   =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Also move the nest's shift in I and J on the parent grid
!***  from the input state to the output state.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract Nest Shift from Input State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state=STATE_IN                           &  !<-- Extract data from input state
                              ,name ='I_SHIFT'                          &  !<-- The name of the data 
                              ,value=I_SHIFT                            &  !<-- Nest's shift in I on its grid
                              ,rc   =RC )
!
        CALL ESMF_AttributeGet(state=STATE_IN                           &  !<-- Extract data from input state
                              ,name ='J_SHIFT'                          &  !<-- The name of the data 
                              ,value=J_SHIFT                            &  !<-- Nest's shift in J on its grid
                              ,rc   =RC )
!
        CALL ESMF_AttributeGet(state=STATE_IN                           &  !<-- Extract data from input state
                              ,name ='I_SW_PARENT_NEW'                  &  !<-- The name of the data 
                              ,value=I_SW_PARENT_NEW                    &  !<-- Nest's shift in I on its grid
                              ,rc   =RC )
!
        CALL ESMF_AttributeGet(state=STATE_IN                           &  !<-- Extract data from input state
                              ,name ='J_SW_PARENT_NEW'                  &  !<-- The name of the data 
                              ,value=J_SW_PARENT_NEW                    &  !<-- Nest's shift in J on its grid
                              ,rc   =RC )
!
        CALL ESMF_AttributeGet(state=STATE_IN                           &  !<-- Extract data from input state
                              ,name ='LAST_STEP_MOVED'                  &  !<-- The name of the data 
                              ,value=LAST_STEP_MOVED                    &  !<-- Nest's shift in J on its grid
                              ,rc   =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert Nest Shift into Output State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=STATE_OUT                          &  !<-- Insert data into output state
                              ,name ='I_SHIFT'                          &  !<-- The name of the data 
                              ,value=I_SHIFT                            &  !<-- Nest's shift in I on its grid
                              ,rc   =RC )
!
        CALL ESMF_AttributeSet(state=STATE_OUT                          &  !<-- Insert data into output state
                              ,name ='J_SHIFT'                          &  !<-- The name of the data 
                              ,value=J_SHIFT                            &  !<-- Nest's shift in J on its grid
                              ,rc   =RC )
!
        CALL ESMF_AttributeSet(state=STATE_OUT                          &  !<-- Insert data into output state
                              ,name ='I_SW_PARENT_NEW'                  &  !<-- The name of the data 
                              ,value=I_SW_PARENT_NEW                    &  !<-- Nest's shift in I on its grid
                              ,rc   =RC )
!
        CALL ESMF_AttributeSet(state=STATE_OUT                          &  !<-- Insert data into output state
                              ,name ='J_SW_PARENT_NEW'                  &  !<-- The name of the data 
                              ,value=J_SW_PARENT_NEW                    &  !<-- Nest's shift in J on its grid
                              ,rc   =RC )
!
        CALL ESMF_AttributeSet(state=STATE_OUT                          &  !<-- Insert data into output state
                              ,name ='LAST_STEP_MOVED'                  &  !<-- The name of the data 
                              ,value=LAST_STEP_MOVED                    &  !<-- Nest's shift in J on its grid
                              ,rc   =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Is the nest ready to move now?
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert MOVE_NOW Flag into the Output State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=STATE_OUT                          &  !<-- Insert data into output state
                              ,name ='MOVE_NOW'                         &  !<-- The name of the data 
                              ,value=MOVE_NOW                           &  !<-- Is the child moving right now?
                              ,rc   =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
      ENDIF
!
!-----------------------------------------------------------------------
!***  Return if this is not a (child) boundary task.
!-----------------------------------------------------------------------
!
      IF(.NOT.PRESENT(S_BDY))THEN
        RETURN
      ENDIF
!
!-----------------------------------------------------------------------
!
      IF(.NOT.MOVE_NOW)THEN
        LIMIT=1                                                            !<-- Only normal boundary data present at time N+1
      ELSE
        LIMIT=2                                                            !<-- Boundary data also present after move for time N
      ENDIF
!
!-----------------------------------------------------------------------
!
      type_loop: DO NTYPE=1,LIMIT
!
!-----------------------------------------------------------------------
!
        IF(NTYPE==1)THEN
          TIME='Future'
        ELSE
          TIME='Current'
        ENDIF
!
!-----------------------------------------------------------------------
!***  Move boundary data on the nest boundary tasks from STATE_IN
!***  to STATE_OUT.
!-----------------------------------------------------------------------
!
        south_boundary: IF(S_BDY)THEN
!
!-------------
!***  South H
!-------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Find # of Words in South H Data"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeGet(state    =STATE_IN                     &  !<-- Look at the input state
                                ,name     ='SOUTH_H_'//TIME             &  !<-- Is this name present?
                                ,itemCount=KOUNT                        &  !<-- How many words present?
                                ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          IF(.NOT.ALLOCATED(BOUNDARY_H%SOUTH))THEN
            ALLOCATE(BOUNDARY_H%SOUTH(1:KOUNT),stat=ISTAT)
            IF(ISTAT/=0)WRITE(0,*)' Failed to allocate BOUNDARY_H%SOUTH stat=',ISTAT
          ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Extract South H Data from Input State"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeGet(state    =STATE_IN                     &  !<-- Extract data from input state
                                ,name     ='SOUTH_H_'//TIME             &  !<-- The name of the data
                                ,itemCount=KOUNT                        &  !<-- The data has this many items
                                ,valueList=BOUNDARY_H%SOUTH             &  !<-- The new combined boundary data
                                ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Insert South H Data into Output State"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeSet(state    =STATE_OUT                    &  !<-- Insert data into output state
                                ,name     ='SOUTH_H_'//TIME             &  !<-- The name of the data 
                                ,itemCount=KOUNT                        &  !<-- The data has this many words
                                ,valueList=BOUNDARY_H%SOUTH             &  !<-- The new combined boundary data
                                ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-------------
!***  South V
!-------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Find # of Words in South V Data"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeGet(state    =STATE_IN                     &  !<-- Look at the input state
                                ,name     ='SOUTH_V_'//TIME             &  !<-- Is this name present?
                                ,itemCount=KOUNT                        &  !<-- How many words present?
                                ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          IF(.NOT.ALLOCATED(BOUNDARY_V%SOUTH))THEN
            ALLOCATE(BOUNDARY_V%SOUTH(1:KOUNT),stat=ISTAT)
            IF(ISTAT/=0)WRITE(0,*)' Failed to allocate BOUNDARY_V%SOUTH stat=',ISTAT
          ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Extract South V Data from Input State"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeGet(state    =STATE_IN                     &  !<-- Extract data from input state
                                ,name     ='SOUTH_V_'//TIME             &  !<-- The name of the data
                                ,itemCount=KOUNT                        &  !<-- The data has this many words
                                ,valueList=BOUNDARY_V%SOUTH             &  !<-- The new combined boundary data
                                ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Insert South V Data into Output State"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeSet(state    =STATE_OUT                    &  !<-- Insert data into output state
                                ,name     ='SOUTH_V_'//TIME             &  !<-- The name of the data 
                                ,itemCount=KOUNT                        &  !<-- The data has this many words
                                ,valueList=BOUNDARY_V%SOUTH             &  !<-- The new combined boundary data
                                ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          IF(ALLOCATED(BOUNDARY_H%SOUTH))THEN
            DEALLOCATE(BOUNDARY_H%SOUTH,stat=ISTAT)
            DEALLOCATE(BOUNDARY_V%SOUTH,stat=ISTAT)
            IF(ISTAT/=0)THEN
              WRITE(0,*)' BOUNDARY_DATA_STATE_TO_STATE failed to'       &
                       ,' deallocate BOUNDARY_V%SOUTH stat=',ISTAT
            ENDIF
          ENDIF
!
        ENDIF south_boundary
!
!-----------------------------------------------------------------------
!
        north_boundary: IF(N_BDY)THEN
!
!-------------
!***  North H
!-------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Find # of Words in North H Data"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeGet(state    =STATE_IN                     &  !<-- Look at input state
                                ,name     ='NORTH_H_'//TIME             &  !<-- Is this name present?
                                ,itemCount=KOUNT                        &  !<-- How many words present?
                                ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          IF(.NOT.ALLOCATED(BOUNDARY_H%NORTH))THEN
            ALLOCATE(BOUNDARY_H%NORTH(1:KOUNT),stat=ISTAT)
            IF(ISTAT/=0)WRITE(0,*)' Failed to allocate BOUNDARY_H%NORTH stat=',ISTAT
          ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Extract North H Data from Input State"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeGet(state    =STATE_IN                     &  !<-- Extract data from input state
                                ,name     ='NORTH_H_'//TIME             &  !<-- The name of the data
                                ,itemCount=KOUNT                        &  !<-- The data has this many words
                                ,valueList=BOUNDARY_H%NORTH             &  !<-- The new combined boundary data
                                ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Insert North H Data into Output State"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeSet(state    =STATE_OUT                    &  !<-- Insert data into output state
                                ,name     ='NORTH_H_'//TIME             &  !<-- The name of the data 
                                ,itemCount=KOUNT                        &  !<-- The data has this many words
                                ,valueList=BOUNDARY_H%NORTH             &  !<-- The new combined boundary data
                                ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-------------
!***  North V
!-------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Find # of Words in North V Data"
!          CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeGet(state    =STATE_IN                     &  !<-- Look at the input state
                                ,name     ='NORTH_V_'//TIME             &  !<-- Is this name present?
                                ,itemCount=KOUNT                        &  !<-- How many words present?
                                ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          IF(.NOT.ALLOCATED(BOUNDARY_V%NORTH))THEN
            ALLOCATE(BOUNDARY_V%NORTH(1:KOUNT),stat=ISTAT)
            IF(ISTAT/=0)WRITE(0,*)' Failed to allocate BOUNDARY_V%NORTH stat=',ISTAT
          ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Extract North V Data from Input State"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeGet(state    =STATE_IN                     &  !<-- Extract data from input state
                                ,name     ='NORTH_V_'//TIME             &  !<-- The name of the data
                                ,itemCount=KOUNT                        &  !<-- The data has this many words
                                ,valueList=BOUNDARY_V%NORTH             &  !<-- The new combined boundary data
                                ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Insert North V Data into Output State"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeSet(state    =STATE_OUT                    &  !<-- Insert data into output state
                                ,name     ='NORTH_V_'//TIME             &  !<-- The name of the data 
                                ,itemCount=KOUNT                        &  !<-- The data has this many words
                                ,valueList=BOUNDARY_V%NORTH             &  !<-- The new combined boundary data
                                ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          IF(ALLOCATED(BOUNDARY_H%NORTH))THEN
            DEALLOCATE(BOUNDARY_H%NORTH,stat=ISTAT)
            DEALLOCATE(BOUNDARY_V%NORTH,stat=ISTAT)
            IF(ISTAT/=0)THEN
              WRITE(0,*)' BOUNDARY_DATA_STATE_TO_STATE failed to'       &
                       ,' deallocate BOUNDARY_V%NORTH stat=',ISTAT
            ENDIF
          ENDIF
!
        ENDIF north_boundary
!
!-----------------------------------------------------------------------
!
        west_boundary: IF(W_BDY)THEN
!
!------------
!***  West H
!------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Find # of Words in West H Data"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeGet(state    =STATE_IN                     &  !<-- Look at the input state
                                ,name     ='WEST_H_'//TIME              &  !<-- Is this name present?
                                ,itemCount=KOUNT                        &  !<-- How many words present?
                                ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          IF(.NOT.ALLOCATED(BOUNDARY_H%WEST))THEN
            ALLOCATE(BOUNDARY_H%WEST(1:KOUNT),stat=ISTAT)
            IF(ISTAT/=0)WRITE(0,*)' Failed to allocate BOUNDARY_H%WEST stat=',ISTAT
          ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Extract West H Data from Input State"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeGet(state    =STATE_IN                     &  !<-- Extract data from input state
                                ,name     ='WEST_H_'//TIME              &  !<-- The name of the data
                                ,itemCount=KOUNT                        &  !<-- The data has this many words
                                ,valueList=BOUNDARY_H%WEST              &  !<-- The new combined boundary data
                                ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Insert West H Data into Output State"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeSet(state    =STATE_OUT                    &  !<-- Insert data into output state
                                ,name     ='WEST_H_'//TIME              &  !<-- The name of the data 
                                ,itemCount=KOUNT                        &  !<-- The data has this many words
                                ,valueList=BOUNDARY_H%WEST              &  !<-- The new combined boundary data
                                ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!------------
!***    West V
!------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Find # of Words in West V Data"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeGet(state    =STATE_IN                     &  !<-- Look at the input state
                                ,name     ='WEST_V_'//TIME              &  !<-- Is this name present?
                                ,itemCount=KOUNT                        &  !<-- How many words present?
                                ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          IF(.NOT.ALLOCATED(BOUNDARY_V%WEST))THEN
            ALLOCATE(BOUNDARY_V%WEST(1:KOUNT),stat=ISTAT)
            IF(ISTAT/=0)WRITE(0,*)' Failed to allocate BOUNDARY_V%WEST stat=',ISTAT
          ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Extract West V Data from Input State"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeGet(state    =STATE_IN                     &  !<-- Extract data input state
                                ,name     ='WEST_V_'//TIME              &  !<-- The name of the data
                                ,itemCount=KOUNT                        &  !<-- The data has this many words
                                ,valueList=BOUNDARY_V%WEST              &  !<-- The new combined boundary data
                                ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Insert West V Data into Output State"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeSet(state    =STATE_OUT                    &  !<-- Insert data into output state
                                ,name     ='WEST_V_'//TIME              &  !<-- The name of the data 
                                ,itemCount=KOUNT                        &  !<-- The data has this many words
                                ,valueList=BOUNDARY_V%WEST              &  !<-- The new combined boundary data
                                ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          IF(ALLOCATED(BOUNDARY_H%WEST))THEN
            DEALLOCATE(BOUNDARY_H%WEST,stat=ISTAT)
            DEALLOCATE(BOUNDARY_V%WEST,stat=ISTAT)
            IF(ISTAT/=0)THEN
              WRITE(0,*)' BOUNDARY_DATA_STATE_TO_STATE failed to'       &
                       ,' deallocate BOUNDARY_V%WEST stat=',ISTAT
            ENDIF
          ENDIF
!
        ENDIF west_boundary
!
!-----------------------------------------------------------------------
!
        east_boundary: IF(E_BDY)THEN
!
!------------
!***  East H
!------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Find # of Words in East H Data"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeGet(state    =STATE_IN                     &  !<-- Look at the input state
                                ,name     ='EAST_H_'//TIME              &  !<-- Is this name present?
                                ,itemCount=KOUNT                        &  !<-- How many words present?
                                ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          IF(.NOT.ALLOCATED(BOUNDARY_H%EAST))THEN
            ALLOCATE(BOUNDARY_H%EAST(1:KOUNT),stat=ISTAT)
            IF(ISTAT/=0)WRITE(0,*)' Failed to allocate BOUNDARY_H%EAST stat=',ISTAT
          ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Extract East H Data from Input State"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeGet(state    =STATE_IN                     &  !<-- Extract data from input state
                                ,name     ='EAST_H_'//TIME              &  !<-- The name of the data
                                ,itemCount=KOUNT                        &  !<-- The data has this many words
                                ,valueList=BOUNDARY_H%EAST              &  !<-- The new combined boundary data
                                ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Insert East H Data into Output State"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeSet(state    =STATE_OUT                    &  !<-- Insert data into output state
                                ,name     ='EAST_H_'//TIME              &  !<-- The name of the data 
                                ,itemCount=KOUNT                        &  !<-- The data has this many words
                                ,valueList=BOUNDARY_H%EAST              &  !<-- The new combined boundary data
                                ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!------------
!***    East V
!------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Find # of Words in East V Data"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeGet(state    =STATE_IN                     &  !<-- Look at the input state
                                ,name     ='EAST_V_'//TIME              &  !<-- Is this name present?
                                ,itemCount=KOUNT                        &  !<-- How many words present?
                                ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          IF(.NOT.ALLOCATED(BOUNDARY_V%EAST))THEN
            ALLOCATE(BOUNDARY_V%EAST(1:KOUNT),stat=ISTAT)
            IF(ISTAT/=0)WRITE(0,*)' Failed to allocate BOUNDARY_V%EAST stat=',ISTAT
          ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Extract East V Data from Input State"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeGet(state    =STATE_IN                     &  !<-- Extract data from input state
                                ,name     ='EAST_V_'//TIME              &  !<-- The name of the data
                                ,itemCount=KOUNT                        &  !<-- The data has this many words
                                ,valueList=BOUNDARY_V%EAST              &  !<-- The new combined boundary data
                                ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Insert East V Data into Output State"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeSet(state    =STATE_OUT                    &  !<-- Insert data into output state
                                ,name     ='EAST_V_'//TIME              &  !<-- The name of the data 
                                ,itemCount=KOUNT                        &  !<-- The data has this many words
                                ,valueList=BOUNDARY_V%EAST              &  !<-- The new combined boundary data
                                ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          IF(ALLOCATED(BOUNDARY_H%EAST))THEN
            DEALLOCATE(BOUNDARY_H%EAST,stat=ISTAT)
            DEALLOCATE(BOUNDARY_V%EAST,stat=ISTAT)
            IF(ISTAT/=0)THEN
              WRITE(0,*)' BOUNDARY_DATA_STATE_TO_STATE failed to'       &
                       ,' deallocate BOUNDARY_V%EAST stat=',ISTAT
            ENDIF
          ENDIF
!
        ENDIF east_boundary
!
!-----------------------------------------------------------------------
!
      ENDDO type_loop
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE BOUNDARY_DATA_STATE_TO_STATE
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE INTERIOR_DATA_STATE_TO_STATE(STATE_IN                  &
                                             ,STATE_OUT )
!
!-----------------------------------------------------------------------
!***  This routine moves from the input to the output state the new
!***  interior update data that was sent from the parent to this
!***  moving nest.  As of now this routine is used only for moving 
!***  the data from the Parent-Child coupler export state (STATE_IN)
!***  to the DOMAIN import state (STATE_OUT).
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(ESMF_State),INTENT(INOUT) :: STATE_IN                           !<-- Input ESMF State
!
      TYPE(ESMF_State),INTENT(INOUT) :: STATE_OUT                          !<-- Output ESMF State
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: I_SHIFT,I_SW_PARENT_NEW                     &
                           ,J_SHIFT,J_SW_PARENT_NEW                     &
                           ,LAST_STEP_MOVED                             &
                           ,N,NEXT_MOVE_TIMESTEP                        &
                           ,NUM_PTASK_UPDATE                            &
                           ,NUM_INTEGER_WORDS                           &
                           ,NUM_REAL_WORDS
!
      INTEGER(kind=KINT) :: RC,RC_S2S
!
      INTEGER(kind=KINT),DIMENSION(1:8) :: INDICES_H,INDICES_V
!
      INTEGER(kind=KINT),DIMENSION(:),POINTER :: UPDATE_INTEGER_DATA
!
      REAL(kind=KFPT),DIMENSION(:),POINTER :: UPDATE_REAL_DATA
!
      CHARACTER(len=1)  :: N_PTASK
      CHARACTER(len=12) :: NAME
      CHARACTER(len=17) :: NAME_REAL
      CHARACTER(len=20) :: NAME_INTEGER
!
      LOGICAL(kind=KLOG) :: MOVE_NOW
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC    =ESMF_SUCCESS
      RC_S2S=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  The Run step of the DOMAIN component needs to know each timestep
!***  whether or not a moving domain wants to move at that time.
!***  So we must transfer the current value of the MOVE_NOW flag to the
!***  DOMAIN import state every timestep.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Unload MOVE_NOW Flag from P-C-Cpl Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=STATE_IN                             &  !<-- The input State
                            ,name ='MOVE_NOW'                           &  !<-- Extract MOVE_NOW flag
                            ,value=MOVE_NOW                             &  !<-- Put the flag here
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_S2S)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Load MOVE_NOW Flag into DOMAIN Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(state=STATE_OUT                            &  !<-- The output State
                            ,name ='MOVE_NOW'                           &  !<-- Name of MOVE_NOW flag
                            ,value=MOVE_NOW                             &  !<-- Inserting MOVE_NOW flag into STATE_OUT
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_S2S)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  If this moving nest does want to move now and this task is
!***  to receive interior update data from its parent then proceed
!***  with transferring additional required data from the Parent-Child
!***  coupler export state to the DOMAIN import state.
!-----------------------------------------------------------------------
!
      move_check: IF(MOVE_NOW)THEN
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Unload I_SHIFT from P-C Cpl Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state=STATE_IN                           &  !<-- The input State
                              ,name ='I_SHIFT'                          &  !<-- Name of the variable
                              ,value=I_SHIFT                            &  !<-- Nest moves this far in I in its space
                              ,rc   =RC)
!
        CALL ESMF_AttributeGet(state=STATE_IN                           &  !<-- The input State
                              ,name ='I_SW_PARENT_NEW'                  &  !<-- Name of the variable
                              ,value=I_SW_PARENT_NEW                    &  !<-- Nest moves this far in I in its space
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_S2S)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Load I_SHIFT into DOMAIN Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=STATE_OUT                          &  !<-- The output State
                              ,name ='I_SHIFT'                          &  !<-- Name of the variable
                              ,value=I_SHIFT                            &  !<-- Load this into STATE_OUT
                              ,rc   =RC)
!
        CALL ESMF_AttributeSet(state=STATE_OUT                          &  !<-- The output State
                              ,name ='I_SW_PARENT_NEW'                  &  !<-- Name of the variable
                              ,value=I_SW_PARENT_NEW                    &  !<-- Load this into STATE_OUT
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_S2S)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Unload J_SHIFT from P-C Cpl Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state=STATE_IN                           &  !<-- The input State
                              ,name ='J_SHIFT'                          &  !<-- Name of the variable
                              ,value=J_SHIFT                            &  !<-- Nest moves this far in J in its space
                              ,rc   =RC)
!
        CALL ESMF_AttributeGet(state=STATE_IN                           &  !<-- The input State
                              ,name ='J_SW_PARENT_NEW'                  &  !<-- Name of the variable
                              ,value=J_SW_PARENT_NEW                    &  !<-- Nest moves this far in J in its space
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_S2S)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Load J_SHIFT into DOMAIN Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=STATE_OUT                          &  !<-- The output State
                              ,name ='J_SHIFT'                          &  !<-- Name of the variable
                              ,value=J_SHIFT                            &  !<-- Load this into STATE_OUT
                              ,rc   =RC)
!
        CALL ESMF_AttributeSet(state=STATE_OUT                          &  !<-- The output State
                              ,name ='J_SW_PARENT_NEW'                  &  !<-- Name of the variable
                              ,value=J_SW_PARENT_NEW                    &  !<-- Load this into STATE_OUT
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_S2S)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Unload LAST_STEP_MOVED from P-C Cpl Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state=STATE_IN                           &  !<-- The input State
                              ,name ='LAST_STEP_MOVED'                  &  !<-- Name of the variable
                              ,value=LAST_STEP_MOVED                    &  !<-- Nest moves this far in J in its space
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_S2S)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Load LAST_STEP_MOVED into DOMAIN Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=STATE_OUT                          &  !<-- The output State
                              ,name ='LAST_STEP_MOVED'                          &  !<-- Name of the variable
                              ,value=LAST_STEP_MOVED                            &  !<-- Load this into STATE_OUT
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_S2S)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Unload # of Parent Tasks Sending Interior Updates"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state=STATE_IN                           &  !<-- The input State
                              ,name ='Num Parent Tasks Update'          &  !<-- Name of the variable
                              ,value=NUM_PTASK_UPDATE                   &  !<-- # of parent tasks that update this nest task
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_S2S)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="INTERIOR_DATA_STATE_TO_STATE: Load # of Parent Tasks Sending Interior Updates"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=STATE_OUT                          &  !<-- The output State
                              ,name ='Num Parent Tasks Update'          &  !<-- Name of the variable
                              ,value=NUM_PTASK_UPDATE                   &  !<-- # of parent tasks that update this nest task
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_S2S)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
        transfer: IF(NUM_PTASK_UPDATE>0)THEN
!
!-----------------------------------------------------------------------
!
          parent_tasks: DO N=1,NUM_PTASK_UPDATE
!
!-----------------------------------------------------------------------
!
            WRITE(N_PTASK,'(I1)')N
!
            NAME_INTEGER='PTASK_INTEGER_DATA_'//N_PTASK
            NAME_REAL   ='PTASK_REAL_DATA_'//N_PTASK
            NAME        ='PTASK_DATA_'//N_PTASK
!
!-----------------------------------------------------------------------
!
!------------
!*** Integer
!------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Unload # of Words in Integer Update Data from Parent"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_AttributeGet(state=STATE_IN                       &  !<-- The input State
                                  ,name =NAME_INTEGER//' Words'         &  !<-- Name of the variable
                                  ,value=NUM_INTEGER_WORDS              &  !<-- # of words in integer update data from Nth parent task
                                  ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_S2S)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Load # of Words in Integer Update Data from Parent"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_AttributeSet(state=STATE_OUT                      &  !<-- The output State
                                  ,name =NAME_INTEGER//' Words'         &  !<-- Name of the variable
                                  ,value=NUM_INTEGER_WORDS              &  !<-- # of words in integer update data from Nth parent task
                                  ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_S2S)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
            transfer_int: IF(NUM_INTEGER_WORDS>0)THEN
!
!-----------------------------------------------------------------------
!
              ALLOCATE(UPDATE_INTEGER_DATA(1:NUM_INTEGER_WORDS))
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
              MESSAGE_CHECK="Unload Interior Integer Update Data from Input State"
!             CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
              CALL ESMF_AttributeGet(state    =STATE_IN                 &  !<-- The input State
                                    ,name     =NAME_INTEGER             &  !<-- Name of the variable
                                    ,itemCount=NUM_INTEGER_WORDS        &  !<-- # of words in integer update data from Nth parent task
                                    ,valueList=UPDATE_INTEGER_DATA      &  !<-- The integer update data from Nth parent task
                                    ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
              CALL ERR_MSG(RC,MESSAGE_CHECK,RC_S2S)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
              MESSAGE_CHECK="Load Interior Integer Update Data into Output State"
!             CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
              CALL ESMF_AttributeSet(state    =STATE_OUT                &  !<-- The output State
                                    ,name     =NAME_INTEGER             &  !<-- Name of the variable
                                    ,itemCount=NUM_INTEGER_WORDS        &  !<-- # of words in integer update data from Nth parent task
                                    ,valueList=UPDATE_INTEGER_DATA      &  !<-- The integer update data from Nth parent task
                                    ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
              CALL ERR_MSG(RC,MESSAGE_CHECK,RC_S2S)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
              DEALLOCATE(UPDATE_INTEGER_DATA)
!
!-----------------------------------------------------------------------
!
            ENDIF transfer_int
!
!-----------------------------------------------------------------------
!
!----------
!***  Real
!----------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Unload # of Words in Real Update Data from Parent"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_AttributeGet(state=STATE_IN                       &  !<-- The input State
                                  ,name =NAME_REAL//' Words'            &  !<-- Name of the variable
                                  ,value=NUM_REAL_WORDS                 &  !<-- # of words in real update data from Nth parent task
                                  ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_S2S)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Load # of Words in Real Update Data from Parent"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_AttributeSet(state=STATE_OUT                      &  !<-- The output State
                                  ,name =NAME_REAL//' Words'            &  !<-- Name of the variable
                                  ,value=NUM_REAL_WORDS                 &  !<-- # of words in real update data from Nth parent task
                                  ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_S2S)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            ALLOCATE(UPDATE_REAL_DATA(1:NUM_REAL_WORDS))
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Unload Interior Real Update Data from Input State"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_AttributeGet(state    =STATE_IN                   &  !<-- The input State
                                  ,name     =NAME_REAL                  &  !<-- Name of the variable
                                  ,itemCount=NUM_REAL_WORDS             &  !<-- # of words in real update data from Nth parent task
                                  ,valueList=UPDATE_REAL_DATA           &  !<-- The real update data from Nth parent task
                                  ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_S2S)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Load Interior Real Update Data into Output State"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_AttributeSet(state    =STATE_OUT                  &  !<-- The output State
                                  ,name     =NAME_REAL                  &  !<-- Name of the variable
                                  ,itemCount=NUM_REAL_WORDS             &  !<-- # of words in real update data from Nth parent task
                                  ,valueList=UPDATE_REAL_DATA           &  !<-- The real update data from Nth parent task
                                  ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_S2S)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            DEALLOCATE(UPDATE_REAL_DATA)
!
!-----------------------------------------------------------------------
!***  Transfer the H and V loop limits for the nest update regions.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Unload Index Limits for H Update Data from Parent"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_AttributeGet(state    =STATE_IN                   &  !<-- The input State
                                  ,name     =NAME//' Indices H'         &  !<-- Name of the variable
                                  ,itemCount=N8                         &  !<-- # of words in index limits of update data
                                  ,valueList=INDICES_H                  &  !<-- The update data index specifications for H
                                  ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_S2S)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Load Index Limits for H Update Data from Parent"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_AttributeSet(state    =STATE_OUT                  &  !<-- The output State
                                  ,name     =NAME//' Indices H'         &  !<-- Name of the variable
                                  ,itemCount=N8                         &  !<-- # of words in index limits of update data
                                  ,valueList=INDICES_H                  &  !<-- The update data index specifications for H
                                  ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_S2S)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Unload Index Limits for V Update Data from Parent"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_AttributeGet(state    =STATE_IN                   &  !<-- The input State
                                  ,name     =NAME//' Indices V'         &  !<-- Name of the variable
                                  ,itemCount=N8                         &  !<-- # of words in index limits of update data
                                  ,valueList=INDICES_V                  &  !<-- The update data index specifications for V
                                  ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_S2S)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Load Index Limits for V Update Data from Parent"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_AttributeSet(state    =STATE_OUT                    &  !<-- The output State
                                  ,name     =NAME//' Indices V'           &  !<-- Name of the variable
                                  ,itemCount=N8                           &  !<-- # of words in index limits of update data
                                  ,valueList=INDICES_V                    &  !<-- The update data index specifications for V
                                  ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_S2S)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
          ENDDO parent_tasks
!
!-----------------------------------------------------------------------
!
        ENDIF transfer
!
!-----------------------------------------------------------------------
!
      ENDIF move_check
!
!-----------------------------------------------------------------------
!***  Finally transfer the value of the domain's next move timestep.
!***  This variable is part of the Solver internal state and is thus
!***  defined for all domains.  Its value is a dummy if not relevant.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="INTERIOR_DATA_STATE_TO_STATE: Unload the Next Move Timestep"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=STATE_IN                             &  !<-- The input State
                            ,name ='NEXT_MOVE_TIMESTEP'                 &  !<-- Name of the variable
                            ,value=NEXT_MOVE_TIMESTEP                   &  !<-- Timestep of domain's next shift.
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_S2S)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="INTERIOR_DATA_STATE_TO_STATE: Load the Next Move Timestep"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(state=STATE_OUT                            &  !<-- The output State
                            ,name ='NEXT_MOVE_TIMESTEP'                 &  !<-- Name of the variable
                            ,value=NEXT_MOVE_TIMESTEP                   &  !<-- Timestep of domain's next shift.
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_S2S)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE INTERIOR_DATA_STATE_TO_STATE
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE MOVING_NEST_BOOKKEEPING(I_SHIFT_CHILD                  &
                                        ,J_SHIFT_CHILD                  &
                                        ,I_SW_PARENT_NEW                &
                                        ,J_SW_PARENT_NEW                &
                                        ,NUM_TASKS_PARENT               &
                                        ,INPES_PARENT                   &
                                        ,ITS_PARENT                     &
                                        ,ITE_PARENT                     &
                                        ,JTS_PARENT                     &
                                        ,JTE_PARENT                     &
                                        ,SPACE_RATIO_MY_PARENT          &
                                        ,NROWS_P_UPD_W                  &
                                        ,NROWS_P_UPD_E                  &
                                        ,NROWS_P_UPD_S                  &
                                        ,NROWS_P_UPD_N                  &
                                        ,SEND_TASK                      &
                                        ,ITS,ITE,JTS,JTE                &
                                        ,IMS,IME,JMS,JME                &
                                        ,IDS,IDE,JDS,JDE                &
                                          )
!
!-----------------------------------------------------------------------
!***  Nest tasks determine which parent tasks will send them update
!***  data and on which points following a move by the nest domain.
!***  The data is for all nest task subdomain points including haloes 
!***  that lie outside of the footprint of the nest domain's location 
!***  preceding the move.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      INTEGER(kind=KINT),INTENT(IN) :: I_SHIFT_CHILD                    &  !<-- Nest domain moved this far in I in nest space    
                                      ,J_SHIFT_CHILD                    &  !<-- Nest domain moved this far in J in nest space    
                                      ,I_SW_PARENT_NEW                  &  !<-- SW corner of nest on this parent I after the move
                                      ,J_SW_PARENT_NEW                  &  !<-- SW corner of nest on this parent J after the move
                                      ,INPES_PARENT                     &  !<-- # of tasks in E/W direction on parent domain
                                      ,NROWS_P_UPD_W                    &  !<-- Moving nest footprint W bndry rows updated by parent
                                      ,NROWS_P_UPD_E                    &  !<-- Moving nest footprint E bndry rows updated by parent
                                      ,NROWS_P_UPD_S                    &  !<-- Moving nest footprint S bndry rows updated by parent
                                      ,NROWS_P_UPD_N                    &  !<-- Moving nest footprint N bndry rows updated by parent
                                      ,NUM_TASKS_PARENT                 &  !<-- Number of fcst tasks on this nest's parent domain
                                      ,SPACE_RATIO_MY_PARENT            &  !<-- Ratio of parent grid increment to this child's
!
                                      ,ITS,ITE,JTS,JTE                  &  !<-- Integration limits of nest task
                                      ,IMS,IME,JMS,JME                  &  !<-- Memory limits of nest task
                                      ,IDS,IDE,JDS,JDE                     !<-- Nest's domain limits
!
      INTEGER(kind=KINT),DIMENSION(0:NUM_TASKS_PARENT-1),INTENT(IN) ::  &
                                                           ITS_PARENT   &  !<-- Starting I of all parent integration subdomains 
                                                          ,ITE_PARENT   &  !<-- Ending I of all parent integration subdomains 
                                                          ,JTS_PARENT   &  !<-- Starting J of all parent integration subdomains 
                                                          ,JTE_PARENT      !<-- Ending J of all parent integration subdomains 
!
      TYPE(INTERIOR_DATA_FROM_PARENT),DIMENSION(1:4),INTENT(OUT) ::     &
                                                            SEND_TASK      !<-- Specifics about interior data from sending parent tasks
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: I_END_X,I_SHIFT,I_START_X                   &
                           ,ID_1,ID_E,ID_N,ID_NE                        &
                           ,J_END_X,J_SHIFT,J_START_X                   &
                           ,KOUNT_PARENT_TASKS,KP,N,NI,NJ
!
      INTEGER(kind=KINT),DIMENSION(1:4) :: I_UPDATE                     &
                                          ,J_UPDATE
!
      REAL(kind=KFPT),DIMENSION(:),ALLOCATABLE :: ITS_PARENT_ON_CHILD   &
                                                 ,ITE_PARENT_ON_CHILD   &
                                                 ,JTS_PARENT_ON_CHILD   &
                                                 ,JTE_PARENT_ON_CHILD 
!
      CHARACTER(2) :: CORNER
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Initialize working variables.
!-----------------------------------------------------------------------
!
      DO N=1,4                                                             !<-- 4 is maximum # of parent tasks that can send data
        SEND_TASK(N)%ID=-9999
        SEND_TASK(N)%ISTART(1)=-9999
        SEND_TASK(N)%IEND  (1)=-9999
        SEND_TASK(N)%JSTART(1)=-9999
        SEND_TASK(N)%JEND  (1)=-9999
        SEND_TASK(N)%ISTART(2)=-9999
        SEND_TASK(N)%IEND  (2)=-9999
        SEND_TASK(N)%JSTART(2)=-9999
        SEND_TASK(N)%JEND  (2)=-9999
      ENDDO
!
!-----------------------------------------------------------------------
!***  Each nest task needs to determine if it has moved outside of
!***  the footprint of the nest domain's previous position.  If it
!***  has then that task next finds out from which parent tasks it
!***  must receive data.  Finally it receives and incorporates that
!***  data.
!
!***  If no part of a nest task's subdomain has moved outside of the 
!***  footprint of the nest domain's previous location then that task
!***  may RETURN now from this routine since none of its points will
!***  be updated by the parent.
!***  
!***  Note that the north and east domain limits are not considered
!***  to be part of the nest's pre-move footprint because those points
!***  cannot be updated by intra-task or inter-task shifts.  The reason
!***  is that V-pt variables at those points are not part of the nest's
!***  integration thus their values are not valid.  Although the H-pt
!***  variables are valid at those points, we cannot use them for
!***  intra- or inter-task updates or else the nest tasks being updated
!***  for H points would sometimes differ from the nest tasks being
!***  updated for V points.  We do not allow that to happen or else
!***  the bookkeeping would be even more complex.  Therefore the
!***  parent updates nest points that would otherwise have been updated
!***  from the north and east domain limits of the nest's pre-move
!***  footprint.  But since a variety of variables do not have valid
!***  integration values on the domain boundary then we also must not
!***  allow the intra- and inter-task shift to handle the updating
!***  of the southern and western boundaries of the nest but instead
!***  must let the parent handle those points as well.  Moreover some
!***  key dynamical tendenies are not computed one row inside the
!***  domain boundary which thus means that the parent must provide
!***  updates for all nest points that not only move beyond the
!***  nest's pre-move footprint but also for those nest points that
!***  move onto IDE and IDE-1 and JDE and JDE-1.  Variables read
!***  from the configure file now specify how deeply the parent will
!***  update nest points with respect to the pre-move footprint.
!-----------------------------------------------------------------------
!
      I_START_X=MAX(IMS,IDS)
      I_END_X  =MIN(IME,IDE)
      J_START_X=MAX(JMS,JDS)
      J_END_X  =MIN(JME,JDE)
!
      IF(I_START_X>=IDS+NROWS_P_UPD_W-I_SHIFT_CHILD                     &  !<-- If the entire nest task subdomain including its
                   .AND.                                                &  !    halo is inside the footprint of the nest domain
         I_END_X<=IDE-NROWS_P_UPD_E-I_SHIFT_CHILD                       &  !    (the domain's position prior to the move) then
                   .AND.                                                &  !    no update from the parent is done.
         J_START_X>=JDS+NROWS_P_UPD_S-J_SHIFT_CHILD                     &  !   
                   .AND.                                                &  !
         J_END_X<=JDE-NROWS_P_UPD_N-J_SHIFT_CHILD )THEN                    !<---
!
        RETURN                                                             !<-- Therefore exit.
!     
      ENDIF
!
!-----------------------------------------------------------------------
!
      I_SHIFT=I_SHIFT_CHILD
      J_SHIFT=J_SHIFT_CHILD
      CORNER='  '
      KOUNT_PARENT_TASKS=0
!
!-----------------------------------------------------------------------
!***  The parent is going to send data interpolated to the nest grid
!***  since that is simpler than handling sparse parent grid data on
!***  the nest domain.  Which points on this nest task now lie outside
!***  the footprint of the pre-move nest domain?  Do a search relative
!***  to the indices on the post-move nest task position since it is
!***  that location at which parent data are received.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Following the move most nest tasks that lie along the boundary
!***  of the footprint of the nest domain will have simple
!***  rectangular update regions in which they will receive
!***  update data from parent tasks.  However if a nest task's new 
!***  position is over a corner of the nest domain's pre-move footprint
!***  then there will be an update region in the nest task that has a 
!***  wedge missing in the intersection with the footprint corner.
!***  This greatly complicates the situation.
!
!***  The diagram below shows a nest task after the nest has moved.
!***  That task now lies over the northeast corner of the footprint
!***  of the domain's pre-move position.  In this case the task
!***  receives update data from four parent tasks (the maximum).
!***  Note that the southwest update region in that nest task 
!***  subdomain has the missing wedge.  To handle this situation
!***  when it arises, the I and J limits on the nest task subdomain
!***  update region will be dimensioned (1:4).  See how the missing
!***  wedge in the update region goes from I_UPDATE(1) to I_UPDATE(2)
!***  and from J_UPDATE(1) to J_UPDATE(2).  The 3rd and 4th elements
!***  of these arrays are filled only for tasks that are on the
!***  footprint's corners.
!-----------------------------------------------------------------------
!
!                                         '                
!                                         '
!                                         '                
!                                         '<-- parent task boundary
!                                         '                       
!                                         '                      
!                                   + + + + + + + + + + + + + + + + + + + + + --- J_UPDATE(4)   
!                                   +     '                                 + 
!   ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' '
!             ^                     +     '                                 +
!             |                     +     '      nest task position         +
!        parent task                +     '        after the move           +
!         boundary                  +     '                                 + --- J_UPDATE(3)
!  ------------------------------------------                               + --- J_UPDATE(2)
!                                  /        /|                              +
!                       I_UPDATE(1)        / |                              +
!                               I_UPDATE(2)  |+ + + + + + + + + + + + + + + + --- J_UPDATE(1)
!                                            | \                             \
!         footprint of the nest domain       |  \                             \
!           in its pre-move location         | I_UPDATE(3)                 I_UPDATE(4)
!                                            |             
!                                            |          
!                                            |             
!
!-----------------------------------------------------------------------
!***  Compute I_UPDATE(1-2) and J_UPDATE(1-2) as well as I_UPDATE(3-4)
!***  and J_UPDATE(3-4) if they exist.
!-----------------------------------------------------------------------
!
      update_limits: IF(                                                  &
                        I_START_X<IDS+NROWS_P_UPD_W-I_SHIFT               &  !<-- If this IF test is true then this nest task lies
                                     .AND.                                &  !    entirely outside of the pre-move footprint and
                        I_END_X  <IDS+NROWS_P_UPD_W-I_SHIFT               &  !    thus all its points will be updated by parent
                                       .OR.                               &  !    tasks.  Since footprint corners are not involved
                        I_START_X>IDE-NROWS_P_UPD_E-I_SHIFT               &  !    then we only have I and J indices 1 and 2 which
                                     .AND.                                &  !    are the starting and ending indices for the
                        I_END_X  >IDE-NROWS_P_UPD_E-I_SHIFT               &  !    nest task's entire subdomain including the halo.
                                       .OR.                               &  !
                        J_START_X<JDS+NROWS_P_UPD_S-J_SHIFT               &  !
                                     .AND.                                &  ! 
                        J_END_X  <JDS+NROWS_P_UPD_S-J_SHIFT               &  !
                                       .OR.                               &  !
                        J_START_X>JDE-NROWS_P_UPD_N-J_SHIFT               &  !
                                     .AND.                                &  !
                        J_END_X  >JDE-NROWS_P_UPD_N-J_SHIFT               &  !<--
                                                                )THEN 
!
        I_UPDATE(1)=I_START_X
        I_UPDATE(2)=I_END_X
        I_UPDATE(3)=-9999
        I_UPDATE(4)=-9999
!
        J_UPDATE(1)=J_START_X
        J_UPDATE(2)=J_END_X
        J_UPDATE(3)=-9999
        J_UPDATE(4)=-9999
!
!-----------------------------------------------------------------------
!
      ELSE update_limits                                                   !<-- Nest task lies on footprint boundary after move
! 
!-----------------------------------------------------------------------
!
        i_block: IF(I_SHIFT>0)THEN                                         !<-- Shift has eastward component
! 
!-----------------------------------------------------------------------
!
!---------------------
!***  Northeast shift
!---------------------
!
          IF(J_SHIFT>0)THEN    
!
            IF(I_END_X>IDE-NROWS_P_UPD_E-I_SHIFT)THEN                      !<-- NE shift, nest task on east side of footprint
              I_UPDATE(1)=IDE-NROWS_P_UPD_E+1-I_SHIFT                      !<-- Begin on east edge of footprint
              I_UPDATE(2)=I_END_X
              J_UPDATE(1)=J_START_X
              J_UPDATE(2)=J_END_X
!
              IF(J_END_X>JDE-NROWS_P_UPD_N-J_SHIFT)THEN                    !<-- NE shift, nest task on NE corner of footprint
                CORNER='NE'
                I_UPDATE(1)=I_START_X
                I_UPDATE(2)=IDE-NROWS_P_UPD_E-I_SHIFT
                I_UPDATE(3)=I_UPDATE(2)+1
                I_UPDATE(4)=I_END_X
                J_UPDATE(2)=JDE-NROWS_P_UPD_N-J_SHIFT
                J_UPDATE(3)=J_UPDATE(2)+1
                J_UPDATE(4)=J_END_X
              ENDIF
!
            ELSEIF(I_END_X<=IDE-NROWS_P_UPD_E-I_SHIFT)THEN                 !<-- NE shift, nest task on north side of footprint; not corner
              I_UPDATE(1)=I_START_X
              I_UPDATE(2)=I_END_X
              J_UPDATE(1)=JDE-NROWS_P_UPD_N+1-J_SHIFT                      !<-- Begin on north edge of footprint
              J_UPDATE(2)=J_END_X
!
            ENDIF
!
!---------------------
!***  Southeast shift
!---------------------
!
          ELSEIF(J_SHIFT<0)THEN 
!
            IF(I_END_X>IDE-NROWS_P_UPD_E-I_SHIFT)THEN                      !<-- SE shift, nest task on east side of footprint; not corner
              I_UPDATE(1)=IDE-NROWS_P_UPD_E+1-I_SHIFT                      !<-- Begin on east edge of footprint
              I_UPDATE(2)=I_END_X
              J_UPDATE(1)=J_START_X
              J_UPDATE(2)=J_END_X
!
              IF(J_START_X<JDS+NROWS_P_UPD_S-J_SHIFT)THEN                  !<-- SE shift, nest task on SE corner of footprint
                CORNER='SE'
                I_UPDATE(1)=I_START_X
                I_UPDATE(2)=IDE-NROWS_P_UPD_E-I_SHIFT
                I_UPDATE(3)=I_UPDATE(2)+1
                I_UPDATE(4)=I_END_X
                J_UPDATE(2)=JDS+NROWS_P_UPD_S-1-J_SHIFT
                J_UPDATE(3)=J_UPDATE(2)+1
                J_UPDATE(4)=J_END_X
              ENDIF
!
            ELSEIF(I_END_X<=IDE-NROWS_P_UPD_E-I_SHIFT)THEN                 !<-- SE shift, nest tasks on south side of footprint; not corner
              I_UPDATE(1)=I_START_X
              I_UPDATE(2)=I_END_X
              J_UPDATE(1)=J_START_X
              J_UPDATE(2)=JDS+NROWS_P_UPD_S-1-J_SHIFT 
!
            ENDIF
!
!-----------------------
!***  Shift is due east
!-----------------------
!
          ELSEIF(J_SHIFT==0)THEN  
            IF(JTE==JDE)THEN
!-> general IF(JME>=JDE-NROWS_P_UPD_N+1)THEN
              IF(I_END_X<IDE-NROWS_P_UPD_E+1-I_SHIFT)THEN                  !<-- Nest task on N bndry of footprint; no part east of it
                I_UPDATE(1)=I_START_X
                I_UPDATE(2)=I_END_X  
                J_UPDATE(1)=JDE-NROWS_P_UPD_N+1
                J_UPDATE(2)=J_END_X  
              ELSE                                                         !<-- Nest task on N bndry of footprint; extends east of it
                CORNER='NE'
                I_UPDATE(1)=I_START_X
                I_UPDATE(2)=IDE-NROWS_P_UPD_E-I_SHIFT
                I_UPDATE(3)=I_UPDATE(2)+1
                I_UPDATE(4)=I_END_X
                J_UPDATE(1)=J_START_X
                J_UPDATE(2)=J_END_X-NROWS_P_UPD_N
                J_UPDATE(3)=J_UPDATE(2)+1
                J_UPDATE(4)=J_END_X
              ENDIF
!
            ELSEIF(JTS>JDS)THEN                                            !<-- Nest task only on east edge of footprint
!-> general ELSEIF(JMS>JDS+NROWS_P_UPD_S-1)THEN                            !<-- Nest task only on east edge of footprint
              I_UPDATE(1)=IDE-NROWS_P_UPD_E+1-I_SHIFT                      !<-- Begin on east edge of footprint
              I_UPDATE(2)=I_END_X
              J_UPDATE(1)=J_START_X
              J_UPDATE(2)=J_END_X
!
            ELSEIF(JTS==JDS)THEN
!-> general ELSEIF(JMS<=JDS+NROWS_P_UPD_S-1)THEN
              IF(I_END_X<IDE-NROWS_P_UPD_E+1-I_SHIFT)THEN                  !<-- Nest task on S bndry of footprint; no part east of it
                I_UPDATE(1)=I_START_X
                I_UPDATE(2)=I_END_X  
                J_UPDATE(1)=JDS
!-> general     J_UPDATE(1)=J_START_X
                J_UPDATE(2)=JDS+NROWS_P_UPD_S-1
              ELSE                                                         !<-- Nest task on S bndry of footprint; extends east of it
                CORNER='SE'
                I_UPDATE(1)=I_START_X
                I_UPDATE(2)=IDE-NROWS_P_UPD_E-I_SHIFT
                I_UPDATE(3)=I_UPDATE(2)+1
                I_UPDATE(4)=I_END_X
                J_UPDATE(1)=J_START_X
                J_UPDATE(2)=JDS+NROWS_P_UPD_S-1
                J_UPDATE(3)=J_UPDATE(2)+1
                J_UPDATE(4)=J_END_X
              ENDIF
!
            ENDIF
!
          ENDIF
!
!-----------------------------------------------------------------------
!
        ELSEIF(I_SHIFT<0)THEN i_block                                      !<-- Shift has westard component
!
!-----------------------------------------------------------------------
!
!---------------------
!***  Northwest shift
!---------------------
!
          IF(J_SHIFT>0)THEN
!
            IF(I_START_X<IDS+NROWS_P_UPD_W-I_SHIFT)THEN                    !<-- NW shift, nest tasks on west side of footprint; not corner
              I_UPDATE(1)=I_START_X        
              I_UPDATE(2)=IDS+NROWS_P_UPD_W-1-I_SHIFT 
              J_UPDATE(1)=J_START_X
              J_UPDATE(2)=J_END_X
!
              IF(J_END_X>JDE-NROWS_P_UPD_N-J_SHIFT)THEN                    !<-- NW shift, nest task on NW corner of footprint
                CORNER='NW'
                I_UPDATE(3)=I_UPDATE(2)+1
                I_UPDATE(4)=I_END_X
                J_UPDATE(2)=JDE-NROWS_P_UPD_N-J_SHIFT
                J_UPDATE(3)=J_UPDATE(2)+1
                J_UPDATE(4)=J_END_X
              ENDIF
!
            ELSEIF(I_START_X>=IDS+NROWS_P_UPD_W-I_SHIFT)THEN               !<-- NW shift, nest tasks on north side of footprint; not corner
              I_UPDATE(1)=I_START_X
              I_UPDATE(2)=I_END_X
              J_UPDATE(1)=JDE-NROWS_P_UPD_N+1-J_SHIFT                      !<-- Begin on north edge of footprint
              J_UPDATE(2)=J_END_X
!
            ENDIF
!
!---------------------
!***  Southwest shift
!---------------------
!
          ELSEIF(J_SHIFT<0)THEN 
!
            IF(I_START_X<IDS+NROWS_P_UPD_W-I_SHIFT)THEN                    !<-- SW shift, nest tasks on west side of footprint
              I_UPDATE(1)=I_START_X     
              I_UPDATE(2)=IDS+NROWS_P_UPD_W-1-I_SHIFT
              J_UPDATE(1)=J_START_X
              J_UPDATE(2)=J_END_X
!
              IF(J_START_X<JDS+NROWS_P_UPD_S-J_SHIFT)THEN                  !<-- SW shift, nest task(s) on SW corner of footprint
                CORNER='SW'
                I_UPDATE(3)=I_UPDATE(2)+1
                I_UPDATE(4)=I_END_X
                J_UPDATE(2)=JDS+NROWS_P_UPD_S-1-J_SHIFT
                J_UPDATE(3)=J_UPDATE(2)+1
                J_UPDATE(4)=J_END_X
              ENDIF
!
            ELSEIF(I_START_X>=IDS+NROWS_P_UPD_W-I_SHIFT)THEN               !<-- SW shift, nest tasks on south side; not corner
              I_UPDATE(1)=I_START_X
              I_UPDATE(2)=I_END_X
              J_UPDATE(1)=J_START_X
              J_UPDATE(2)=JDS+NROWS_P_UPD_S-1-J_SHIFT
!
            ENDIF
!
!-----------------------
!***  Shift is due west
!-----------------------
!
          ELSEIF(J_SHIFT==0)THEN
            IF(JTE==JDE)THEN
!-> general IF(JME>=JDE-NROWS_P_UPD_N+1)THEN
              IF(I_START_X>=IDS+NROWS_P_UPD_W-I_SHIFT)THEN                 !<-- Nest task on N bndry of footprint; no part west of it
                I_UPDATE(1)=I_START_X
                I_UPDATE(2)=I_END_X  
                J_UPDATE(1)=JDE-NROWS_P_UPD_N+1
                J_UPDATE(2)=J_END_X  
              ELSE                                                         !<-- Nest task on N bndry of footprint; extends west of it
                CORNER='NW'
                I_UPDATE(1)=I_START_X
                I_UPDATE(2)=IDS+NROWS_P_UPD_W-1-I_SHIFT
                I_UPDATE(3)=I_UPDATE(2)+1
                I_UPDATE(4)=I_END_X
                J_UPDATE(1)=J_START_X
                J_UPDATE(2)=J_END_X-1
                J_UPDATE(2)=JDE-NROWS_P_UPD_N
                J_UPDATE(3)=J_UPDATE(2)+1
                J_UPDATE(4)=J_END_X
              ENDIF
!
            ELSEIF(JTS>JDS)THEN                                            !<-- Nest task only on west edge of footprint
!-> general ELSEIF(JMS>JDS+NROWS_P_UPD_S-1)THEN                            !<-- Nest task only on west edge of footprint
              I_UPDATE(1)=I_START_X
              I_UPDATE(2)=IDS+NROWS_P_UPD_W-1-I_SHIFT                      !<-- End on the west edge of footprint
              J_UPDATE(1)=J_START_X
              J_UPDATE(2)=J_END_X
!
            ELSEIF(JTS==JDS)THEN
!-> general ELSEIF(JMS<=JDS+NROWS_P_UPD_S-1)THEN
              IF(I_START_X>=IDS+NROWS_P_UPD_W-I_SHIFT)THEN                 !<-- Nest task on S bndry of footprint; no part west of it
                I_UPDATE(1)=I_START_X
                I_UPDATE(2)=I_END_X  
                J_UPDATE(1)=JDS
!-> general     J_UPDATE(1)=J_START_X
                J_UPDATE(2)=JDS+NROWS_P_UPD_S-1
              ELSE                                                         !<-- Nest task on S bndry of footprint; extends west of it
                CORNER='SW'
                I_UPDATE(1)=I_START_X
                I_UPDATE(2)=IDS+NROWS_P_UPD_W-1-I_SHIFT
                I_UPDATE(3)=I_UPDATE(2)+1
                I_UPDATE(4)=I_END_X
                J_UPDATE(1)=JDS
!-> general     J_UPDATE(1)=J_START_X
                J_UPDATE(2)=JDS+NROWS_P_UPD_S-1
                J_UPDATE(3)=J_UPDATE(2)+1
                J_UPDATE(4)=J_END_X
              ENDIF
!
            ENDIF
!
          ENDIF
!
!-----------------------------------------------------------------------
!
        ELSEIF(I_SHIFT==0)THEN                                             !<-- Shift has no eastward or westward component
!
!------------------------
!***  Shift is due north
!------------------------
!
          IF(J_SHIFT>0)THEN 
            IF(ITE==IDE)THEN   
!-> general IF(IME>=IDE-NROWS_P_UPD_E+1)THEN   
              IF(J_END_X<JDE-NROWS_P_UPD_N+1-J_SHIFT)THEN                  !<-- Nest task on E bndry of footprint; no part north of it
                I_UPDATE(1)=IDE-NROWS_P_UPD_E+1
                I_UPDATE(2)=I_END_X
                J_UPDATE(1)=J_START_X
                J_UPDATE(2)=J_END_X
              ELSE                                                         !<-- Nest task on E bndry of footprint; extends north of it
                CORNER='NE'
                I_UPDATE(1)=I_START_X
                I_UPDATE(2)=IDE-NROWS_P_UPD_E
                I_UPDATE(3)=I_UPDATE(2)+1
                I_UPDATE(4)=I_END_X
                J_UPDATE(1)=J_START_X
                J_UPDATE(2)=JDE-NROWS_P_UPD_N-J_SHIFT
                J_UPDATE(3)=J_UPDATE(2)+1
                J_UPDATE(4)=J_END_X
              ENDIF
!
            ELSEIF(ITS>IDS)THEN                                            !<-- Nest task only on north edge of footprint
!-> general ELSEIF(IMS>IDS+NROWS_P_UPD_W-1)THEN                            !<-- Nest task only on north edge of footprint
              I_UPDATE(1)=I_START_X
              I_UPDATE(2)=I_END_X
              J_UPDATE(1)=JDE-NROWS_P_UPD_N+1-J_SHIFT                      !<-- Begin on north edge of footprint
              J_UPDATE(2)=J_END_X
!
            ELSEIF(ITS==IDS)THEN   
!-> general ELSEIF(IMS<=IDS+NROWS_P_UPD_W-1)THEN   
              IF(J_END_X<JDE-NROWS_P_UPD_N+1-J_SHIFT)THEN                  !<-- Nest task on W bndry of footprint; no part north of it
                I_UPDATE(1)=IDS     
!-> general     I_UPDATE(1)=I_START_X
                I_UPDATE(2)=IDS+NROWS_P_UPD_W-1
                J_UPDATE(1)=J_START_X
                J_UPDATE(2)=J_END_X
              ELSE                                                         !<-- Nest task on W bndry of footprint; extends north of it
                CORNER='NW'
                I_UPDATE(1)=IDS
!-> general     I_UPDATE(1)=I_START_X
                I_UPDATE(2)=IDS+NROWS_P_UPD_W-1
                I_UPDATE(3)=I_UPDATE(2)+1
                I_UPDATE(4)=I_END_X
                J_UPDATE(1)=J_START_X
                J_UPDATE(2)=JDE-NROWS_P_UPD_N-J_SHIFT
                J_UPDATE(3)=J_UPDATE(2)+1
                J_UPDATE(4)=J_END_X
              ENDIF
!
            ENDIF
!
!------------------------
!***  Shift is due south
!------------------------
!
          ELSEIF(J_SHIFT<0)THEN 
            IF(ITE==IDE)THEN   
!-> general IF(IME>=IDE-NROWS_P_UPD_E+1)THEN   
              IF(J_START_X>=JDS+NROWS_P_UPD_S-J_SHIFT)THEN                 !<-- Nest task on E bndry of footprint; no part south of it
                I_UPDATE(1)=IDE-NROWS_P_UPD_E+1
                I_UPDATE(2)=I_END_X
                J_UPDATE(1)=J_START_X
                J_UPDATE(2)=J_END_X
              ELSE                                                         !<-- Nest task on E bndry of footprint; extends south of it
                CORNER='SE'
                I_UPDATE(1)=I_START_X
                I_UPDATE(2)=IDE-NROWS_P_UPD_E
                I_UPDATE(3)=I_UPDATE(2)+1
                I_UPDATE(4)=I_END_X
                J_UPDATE(1)=J_START_X
                J_UPDATE(2)=JDS+NROWS_P_UPD_S-1-J_SHIFT
                J_UPDATE(2)=JDS+NROWS_P_UPD_S-1-J_SHIFT
                J_UPDATE(3)=J_UPDATE(2)+1
                J_UPDATE(4)=J_END_X
              ENDIF
!
            ELSEIF(ITS>IDS)THEN                                            !<-- Nest task only on south edge of footprint
!-> general ELSEIF(IMS>IDS+NROWS_P_UPD_W-1)THEN                            !<-- Nest task only on south edge of footprint
              I_UPDATE(1)=I_START_X
              I_UPDATE(2)=I_END_X   
              J_UPDATE(1)=J_START_X
              J_UPDATE(2)=JDS+NROWS_P_UPD_S-1-J_SHIFT                      !<-- End on south edge of footprint
!
            ELSEIF(ITS==IDS)THEN   
!-> general ELSEIF(IMS<=IDS+NROWS_P_UPD_W-1)THEN   
              IF(J_START_X>=JDS+NROWS_P_UPD_S-J_SHIFT)THEN                 !<-- Nest task on W bndry of footprint; no part south of it
                I_UPDATE(1)=IDS    
!-> general     I_UPDATE(1)=I_START_X
                I_UPDATE(2)=IDS+NROWS_P_UPD_W-1
                J_UPDATE(1)=J_START_X
                J_UPDATE(2)=J_END_X
              ELSE                                                         !<-- Nest task on W bndry of footprint; extends south of it
                CORNER='SW'
                I_UPDATE(1)=I_START_X
                I_UPDATE(2)=IDS+NROWS_P_UPD_W-1
                I_UPDATE(3)=I_UPDATE(2)+1
                I_UPDATE(4)=I_END_X
                J_UPDATE(1)=J_START_X
                J_UPDATE(2)=JDS+NROWS_P_UPD_S-1-J_SHIFT
                J_UPDATE(3)=J_UPDATE(2)+1
                J_UPDATE(4)=J_END_X
              ENDIF
!
            ENDIF
!
          ENDIF 
!
!-----------------------------------------------------------------------
!
        ENDIF i_block
!
!-----------------------------------------------------------------------
!
      ENDIF update_limits
!
!-----------------------------------------------------------------------
!***  Now we know which portion of each task's subdomain on the
!***  moving nest lies outside of the nest domain's pre-move
!***  footprint location and it is that portion that must be
!***  updated from the parent tasks.
!
!***  To receive data from its parent, each moving nest task must
!***  know how many parent tasks it will receive from.  Nest tasks
!***  could do this blindly by receiving a message from all parent
!***  tasks which would inform them which parent tasks had actual
!***  data, or they could receive all that information from parent
!***  task 0 if that task had first been sent all that information
!***  from the other parent tasks, or each nest task could compute
!***  which parent tasks will send it data and how much.  The third
!***  option involves the least overall communication and serves to
!***  double check the parent's bookkeeping therefore that option 
!***  is the one used here.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  The nest tasks determine all of their parent tasks' integration
!***  limits in terms of their own (the nest tasks') indices so they
!***  will know where to put the data they receive from the parent.
!***  This must be done for all parent tasks since the nest tasks
!***  can make no assumptions about which parent tasks will have
!***  update data for them and that is because there is no limit
!***  imposed on the distance the nest can traverse in any single
!***  move.
!
!***  See the explanation and accompanying diagrams in subroutine
!***  PARENT_BOOKKEEPING_MOVING for more details.  The results must
!***  be the same for both H and V points.
!-----------------------------------------------------------------------
!
      ALLOCATE(ITS_PARENT_ON_CHILD(0:NUM_TASKS_PARENT-1))
      ALLOCATE(ITE_PARENT_ON_CHILD(0:NUM_TASKS_PARENT-1))
      ALLOCATE(JTS_PARENT_ON_CHILD(0:NUM_TASKS_PARENT-1))
      ALLOCATE(JTE_PARENT_ON_CHILD(0:NUM_TASKS_PARENT-1))
!
      DO N=0,NUM_TASKS_PARENT-1
!
        ITS_PARENT_ON_CHILD(N)=REAL(IDS                                 &  !<-- ITS of parent task N in child's coordinate space
                                   -(I_SW_PARENT_NEW-ITS_PARENT(N))     &
                                     *SPACE_RATIO_MY_PARENT)
!
        ITE_PARENT_ON_CHILD(N)=REAL(IDS                                 &  !<-- ITE of parent task N in child's coordinate space
                                   -(I_SW_PARENT_NEW-ITE_PARENT(N))     & 
                                     *SPACE_RATIO_MY_PARENT             &
                                     +SPACE_RATIO_MY_PARENT-1)             !<-- Filling in gap beyond last parent point on nest grid
!
        JTS_PARENT_ON_CHILD(N)=REAL(JDS                                 &  !<-- JTS of parent task N in child's coordinate space
                                   -(J_SW_PARENT_NEW-JTS_PARENT(N))     &
                                     *SPACE_RATIO_MY_PARENT)
!
        JTE_PARENT_ON_CHILD(N)=REAL(JDS                                 &  !<-- JTE of parent task N in child's coordinate space
                                   -(J_SW_PARENT_NEW-JTE_PARENT(N))     &
                                     *SPACE_RATIO_MY_PARENT             &
                                     +SPACE_RATIO_MY_PARENT-1)             !<-- Filling in gap beyond last parent point on nest grid
!
      ENDDO
!
!-----------------------------------------------------------------------
!***  Find the parent task whose subdomain contains the nest point
!***  [I_UPDATE(1),J_UPDATE(1)].  This parent task will be referred
!***  to as parent task #1 since up to four parent tasks might
!***  provide update data to the nest task.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
      search_i:  DO NI=0,INPES_PARENT-1                                    !<-- Look eastward for parent task.
!-----------------------------------------------------------------------
!
         update_i: IF(REAL(I_UPDATE(1))>=ITS_PARENT_ON_CHILD(NI)-EPS    &  !<-- Search for 1st parent task that covers
                                  .AND.                                 &  !    nest index I_UPDATE(1).
                      REAL(I_UPDATE(1))<=ITE_PARENT_ON_CHILD(NI)+EPS)   &
                                                                   THEN    !<--
!
!-----------------------------------------------------------------------
          search_j: DO NJ=NI,NUM_TASKS_PARENT-1,INPES_PARENT               !<-- Look northward for parent task.
!-----------------------------------------------------------------------
!
            update_j: IF(REAL(J_UPDATE(1))>=JTS_PARENT_ON_CHILD(NJ)-EPS  &  !<-- Search for parent task that covers nest point
                                     .AND.                               &  !    I_UPDATE(1),J_UPDATE(1).
                         REAL(J_UPDATE(1))<=JTE_PARENT_ON_CHILD(NJ)+EPS) &  !<--
                                                                    THEN    !<--
!
              SEND_TASK(1)%ID=NJ                                           !<-- Store the ID of this parent task who will send data.
              ID_1=NJ                                                      !<-- Local task ID of the identified parent task.
              KOUNT_PARENT_TASKS=1                                         !<-- Count how many parent tasks send to this nest task.
!
!-----------------------------------------------------------------------
!***  First consider all nest tasks that either lie totally outside 
!***  of the footprint or lie on the footprint's edge but not on a 
!***  corner.  Corners can be very complicated and are each treated
!***  separately.
!-----------------------------------------------------------------------
!
              not_a_corner: IF(CORNER=='  ')THEN
!
!-----------------------------------------------------------------------
!***  I and J limits on the nest task of data received from the
!***  parent task #1 that covers [I_UPDATE(1),J_UPDATE(1)].
!-----------------------------------------------------------------------
!
                SEND_TASK(1)%ISTART(1)=I_UPDATE(1)                         !<-- Nest index limits updated by parent task 1.
                SEND_TASK(1)%IEND  (1)=MIN(I_UPDATE(2)                  &  !
                                      ,NINT(ITE_PARENT_ON_CHILD(ID_1)))    !
                SEND_TASK(1)%JSTART(1)=J_UPDATE(1)                         !
                SEND_TASK(1)%JEND  (1)=MIN(J_UPDATE(2)                  &  !
                                      ,NINT(JTE_PARENT_ON_CHILD(ID_1)))    !<--
!
!-----------------------------------------------------------------------
!***  Is there a parent task to the the north of the first that covers
!***  points on this nest task's subdomain?
!-----------------------------------------------------------------------
!
                IF(JTE_PARENT_ON_CHILD(ID_1)+EPS<J_UPDATE(2))THEN          !<-- If so, there is a parent task to the north with data.
                  KOUNT_PARENT_TASKS=KOUNT_PARENT_TASKS+1                  !<-- Increment parent task counter.
                  KP=KOUNT_PARENT_TASKS
                  SEND_TASK(KP)%ID=ID_1+INPES_PARENT                       !<-- Store the ID of this parent task to the north.
                  ID_N=SEND_TASK(KP)%ID
                  SEND_TASK(KP)%ISTART(1)=SEND_TASK(1)%ISTART(1)           !<-- Nest index limits updated by parent task ID_N.
                  SEND_TASK(KP)%IEND  (1)=SEND_TASK(1)%IEND  (1)           !
                  SEND_TASK(KP)%JSTART(1)=SEND_TASK(1)%JEND  (1)+1         !
                  SEND_TASK(KP)%JEND  (1)=J_UPDATE(2)                      !<--
!
!-----------------------------------------------------------------------
!***  Is there a parent task to the the northeast of the first that 
!***  covers points on this nest task's subdomain?  If there is then
!***  its southern and northern update limits are the same as the 
!***  parent task's to the north.
!-----------------------------------------------------------------------
!
                  IF(ITE_PARENT_ON_CHILD(ID_N)+EPS<I_UPDATE(2))THEN        !<-- If so, there is a parent task to the NE with data. 
                    KOUNT_PARENT_TASKS=KOUNT_PARENT_TASKS+1                !<-- Increment parent task counter.
                    KP=KOUNT_PARENT_TASKS
                    SEND_TASK(KP)%ID=ID_N+1                                !<-- Store the ID of this parent task to the northeast.
                    ID_NE=SEND_TASK(KP)%ID
                    SEND_TASK(KP)%ISTART(1)=SEND_TASK(KP-1)%IEND(1)+1      !<-- Nest index limits updated by parent task ID_NE.
                    SEND_TASK(KP)%IEND  (1)=I_UPDATE(2)                    !
                    SEND_TASK(KP)%JSTART(1)=SEND_TASK(KP-1)%JSTART(1)      !
                    SEND_TASK(KP)%JEND  (1)=SEND_TASK(KP-1)%JEND  (1)      !<--
                  ENDIF
!
                ENDIF
!
!-----------------------------------------------------------------------
!***  Finally is there a parent task to the east of the first parent
!***  task?  Its southern and northern update limits are the same as
!***  the first parent task's.
!-----------------------------------------------------------------------
!
                IF(ITE_PARENT_ON_CHILD(ID_1)+EPS<I_UPDATE(2))THEN          !<-- If so, there is a parent task to the E with data. 
                  KOUNT_PARENT_TASKS=KOUNT_PARENT_TASKS+1                  !<-- Increment parent task counter.
                  KP=KOUNT_PARENT_TASKS
                  SEND_TASK(KP)%ID=ID_1+1                                  !<-- Store the ID of this parent task to the east.
                  ID_E=SEND_TASK(KP)%ID
                  SEND_TASK(KP)%ISTART(1)=SEND_TASK(1)%IEND(1)+1           !<-- Nest index limits updated by parent task ID_E.
                  SEND_TASK(KP)%IEND  (1)=I_UPDATE(2)                      !
                  SEND_TASK(KP)%JSTART(1)=SEND_TASK(1)%JSTART(1)           !
                  SEND_TASK(KP)%JEND  (1)=SEND_TASK(1)%JEND  (1)           !<--
!
                ENDIF
!
                EXIT search_i
!
              ENDIF not_a_corner
!
!------------------------------------------------------
!------------------------------------------------------
!***  The nest task on the SW corner of the footprint.
!------------------------------------------------------
!------------------------------------------------------
!
              sw: IF(CORNER=='SW')THEN
!
                SEND_TASK(1)%ISTART(1)=I_START_X                           !<-- Nest I where parent task 1 begins updating nest task.
                SEND_TASK(1)%IEND(1)=MIN(I_UPDATE(4)                    &  !<-- Nest I where parent task 1 ends updating nest task.
                                     ,NINT(ITE_PARENT_ON_CHILD(ID_1))) 
!
                SEND_TASK(1)%JSTART(1)=J_START_X                           !<-- Nest J where parent task 1 begins updating nest task.
!
                IF(JTE_PARENT_ON_CHILD(ID_1)-EPS<=J_UPDATE(2))THEN
                  SEND_TASK(1)%JEND(1)=NINT(JTE_PARENT_ON_CHILD(ID_1))     !<-- Nest J where parent task 1 ends updating nest task.
!
                ELSEIF(JTE_PARENT_ON_CHILD(ID_1)+EPS>=J_UPDATE(3))THEN
!
                  IF(ITE_PARENT_ON_CHILD(ID_1)-EPS<=I_UPDATE(2))THEN
                    SEND_TASK(1)%JEND(1)=MIN(J_END_X                    &  !<-- Nest J where parent task 1 ends updating nest task.
                                         ,NINT(JTE_PARENT_ON_CHILD(ID_1))) !<-- 
!
                  ELSEIF(ITE_PARENT_ON_CHILD(ID_1)+EPS>=I_UPDATE(3))THEN   !<-- Parent task 1 covers SW corner of footprint too.
!
                    SEND_TASK(1)%JEND(1)=J_UPDATE(2)                       !<-- Nest J where parent task 1 ends updating nest task's
                                                                           !    first region.
                    SEND_TASK(1)%ISTART(2)=I_START_X                       !<-- 2nd region on nest task updated by parent task 1
                    SEND_TASK(1)%IEND(2)=I_UPDATE(2)                       !
                    SEND_TASK(1)%JSTART(2)=J_UPDATE(3)                     !
                    SEND_TASK(1)%JEND(2)=MIN(J_END_X                    &  !
                                         ,NINT(JTE_PARENT_ON_CHILD(ID_1))) !<--
!
                  ENDIF
!
                ENDIF
!
!-----------------------------------------------------------------------
!***  The points in this nest task subdomain being updated by the first
!***  identified parent task have been demarcated.  Now identify any
!***  other parent tasks updating this nest task lying on the SW corner
!***  of the footprint.
!-----------------------------------------------------------------------
!                  
!------------------------------------------------------
!***  Is there a parent task to the north of the first
!***  that provides update data?
!------------------------------------------------------
!
                sw_north: IF(JTE_PARENT_ON_CHILD(ID_1)+EPS<J_END_X)THEN     !<-- If true there is another parent task north of the first.
!
                  KOUNT_PARENT_TASKS=KOUNT_PARENT_TASKS+1                   !<-- Increment parent task counter.
                  KP=KOUNT_PARENT_TASKS
                  SEND_TASK(KP)%ID=ID_1+INPES_PARENT                        !<-- Store the ID of this parent task to the north.
                  ID_N=SEND_TASK(KP)%ID  
!
                  SEND_TASK(KP)%ISTART(1)=I_START_X                         !<-- Starting I on nest task where parent task ID_N updates.
!
                  IF(ITE_PARENT_ON_CHILD(ID_N)-EPS<=I_UPDATE(2).OR.     &   !<-- Parent task ID_N does not cover SW corner of footprint
                     JTS_PARENT_ON_CHILD(ID_N)+EPS>=J_UPDATE(3))THEN        !    if either of these statements is true.
!
                    SEND_TASK(KP)%IEND(1)=MIN(I_UPDATE(2)               &   !<-- Ending I on nest task where parent task ID_N updates.
                                         ,NINT(ITE_PARENT_ON_CHILD(ID_N)))  !<--
                    SEND_TASK(KP)%JSTART(1)=NINT(JTS_PARENT_ON_CHILD(ID_N)) !<-- Starting J on nest task where parent task ID_N updates.
                    SEND_TASK(KP)%JEND(1)=J_END_X                           !<-- Ending J on nest task where parent task ID_N updates.
!
                  ELSEIF(ITE_PARENT_ON_CHILD(ID_N)+EPS>=I_UPDATE(3).AND. &  !<-- Parent task ID_N covers SW corner of
                         JTS_PARENT_ON_CHILD(ID_N)-EPS<=J_UPDATE(2))THEN    !    footprint too.
!
!                                         |
!                                         |
!                                         |
!                                         |
!                 + + + + + + + + + + + +.|     footprint of nest domain
!                 +                      .|     in its pre-move location
!                 +   parent task 2's    .|
!                 +   2nd update region  .|
!                 +                      .|
!                 +.......................----------------------------------
!                 +..............................       +
!                 +     parent task 2's         '    <------- parent task 3's update region
!                 +     1st update region       '       +
!   ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' <--- parent task boundary
!                 +                             '       +
!                 +      parent task 1's        '    <------- parent task 4's update region
!                 +      update region          '       +
!                 +                             '       +
!                 + + + + + + + + + + + + + + + ' + + + +
!                                 ^             '
!                                 |             '
!                             nest task         '<--- parent task boundary
!                             boundary          '
!                             after move        '
!                                               '
!
                    SEND_TASK(KP)%IEND(1)=MIN(I_END_X                     &  !<-- Ending I on nest task where parent task ID_N updates
                                         ,NINT(ITE_PARENT_ON_CHILD(ID_N)))   !<-- in nest task's 1st region.
                    SEND_TASK(KP)%JSTART(1)=NINT(JTS_PARENT_ON_CHILD(ID_N))  !<-- Starting J of parent task ID_N in nest task 1st region.
                    SEND_TASK(KP)%JEND  (1)=J_UPDATE(2)                      !<-- Ending J of parent task ID_N in nest task 1st region.
                    SEND_TASK(KP)%ISTART(2)=I_START_X                        !<-- Starting I of parent task ID_N in nest task 2nd region.
                    SEND_TASK(KP)%IEND  (2)=MIN(I_UPDATE(2)               &  !<-- Ending I of parent task ID_N in nest task 2nd region.
                                           ,NINT(ITE_PARENT_ON_CHILD(ID_N))) !<--
                    SEND_TASK(KP)%JSTART(2)=J_UPDATE(3)                      !<-- Starting J of parent task ID_N in nest task 2nd region.
                    SEND_TASK(KP)%JEND  (2)=J_END_X                          !<-- Ending J of parent task ID_N in nest task 2nd region.
!
                  ENDIF                                                      !<-- End contingencies of parent task north of first one.
!
!-----------------------------------------------
!***  Does a parent task northeast of the first 
!***  provide any update data?  This can only
!***  happen if there was already a parent task
!***  to the north of the first one providing
!***  update data.
!-----------------------------------------------
!
!
                  IF(ITE_PARENT_ON_CHILD(ID_N)+EPS<I_UPDATE(2)            &  !<-- 1st scenario of parent update task to northeast of 
                                  .AND.                                   &  !    the first parent update task.  The NE parent does
                     JTS_PARENT_ON_CHILD(ID_N)+EPS>=J_UPDATE(3))THEN         !<-- not cover the SW corner of the footprint.
!
                    KOUNT_PARENT_TASKS=KOUNT_PARENT_TASKS+1                  !<-- Increment parent task counter.
                    KP=KOUNT_PARENT_TASKS
                    SEND_TASK(KP)%ID=ID_N+1                                  !<-- Store the ID of this parent task northeast of the first
                    ID_NE=SEND_TASK(KP)%ID                                   !    (i.e. east of the parent task to the north of the first).
                    SEND_TASK(KP)%ISTART(1)=NINT(ITS_PARENT_ON_CHILD(ID_NE)) !<-- Starting I where parent task ID_NE updates nest task.
                    SEND_TASK(KP)%IEND  (1)=I_UPDATE(2)                      !<-- Ending I where parent task ID_NE updates nest task.
                    SEND_TASK(KP)%JSTART(1)=NINT(JTS_PARENT_ON_CHILD(ID_NE)) !<-- Starting J where parent task ID_NE updates nest task.
                    SEND_TASK(KP)%JEND  (1)=J_END_X                          !<-- Ending J where parent task ID_NE updates nest task.
!
                  ENDIF
!
                  IF(ITE_PARENT_ON_CHILD(ID_N)+EPS>=I_UPDATE(2)           &  !<-- 2nd scenario of parent update task to northeast of
                                  .AND.                                   &  !    the first update parent task.  The NE parent task
                     ITE_PARENT_ON_CHILD(ID_N)+EPS< I_END_X               &  !    does not cover the SW corner of the footprint.
                                  .AND.                                   &  !   
                     JTS_PARENT_ON_CHILD(ID_N)-EPS<=J_UPDATE(2))THEN         !<--
!
                    KOUNT_PARENT_TASKS=KOUNT_PARENT_TASKS+1                  !<-- Increment parent task counter.
                    KP=KOUNT_PARENT_TASKS
                    SEND_TASK(KP)%ID=ID_N+1                                  !<-- Store the ID of this parent task northeast of the first.
                    ID_NE=SEND_TASK(KP)%ID
                    SEND_TASK(KP)%ISTART(1)=NINT(ITS_PARENT_ON_CHILD(ID_NE)) !<-- Starting I where parent task ID_NE updates nest task.
                    SEND_TASK(KP)%IEND  (1)=I_END_X                          !<-- Ending I where parent task ID_NE updates nest task.
                    SEND_TASK(KP)%JSTART(1)=NINT(JTS_PARENT_ON_CHILD(ID_NE)) !<-- Starting J where parent task ID_NE updates nest task.
                    SEND_TASK(KP)%JEND  (1)=J_UPDATE(2)                      !<-- Ending J where parent task ID_NE updates nest task.
!
                  ENDIF
!
                  IF(ITE_PARENT_ON_CHILD(ID_N)+EPS<I_UPDATE(2)            &  !<-- 3rd scenario of parent update task to northeast of
                                  .AND.                                   &  !    the first update parent task.  In this case the
                     JTS_PARENT_ON_CHILD(ID_N)-EPS<=J_UPDATE(2))THEN         !    parent task to the NE is on the footprint corner.
!
!
!                                                           |
!                                   '                       |
!                                   '                       |
!                                   '                       |
!                         + + + + + ' + + + + + + + + + + +.|     footprint of nest domain
!                         +         '                      .|     in its pre-move location
!                         +         '   parent task 3's    .|
!   parent task 2's       +         '   2nd update region  .|
!    update region  ---------->     '                      .|
!                         +         '.......................----------------------------------
!                         +         '.................................+
!                         +         '          parent task 3's        +
!                         +         '         1st update region       +
!           ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' '  ' ' ' ' <--- parent task boundary
!                         +         '                                 +
!   parent task 1's       +         '          parent task 4's        +
!    update region  ---------->     '           update region         +
!                         +         '                                 + <--- nest task boundary after move
!                         + + + + + ' + + + + + + + + + + + + + + + + +
!                                   '
!                                   '
!                                   '<----- parent task boundary
!                                   '
!                                   '
!
!
                    KOUNT_PARENT_TASKS=KOUNT_PARENT_TASKS+1                  !<-- Increment the parent task counter.
                    KP=KOUNT_PARENT_TASKS
                    SEND_TASK(KP)%ID=ID_N+1                                  !<-- Store the ID of this parent task northeast of the first.
                    ID_NE=SEND_TASK(KP)%ID
                    SEND_TASK(KP)%ISTART(1)=NINT(ITS_PARENT_ON_CHILD(ID_NE)) !<-- Starting I in nest task's 1st region updated by parent.
                    SEND_TASK(KP)%IEND  (1)=I_END_X                          !<-- Ending I in nest task's 1st region updated by parent.
                    SEND_TASK(KP)%JSTART(1)=NINT(JTS_PARENT_ON_CHILD(ID_NE)) !<-- Starting J in nest task's 1st region updated by parent.
                    SEND_TASK(KP)%JEND  (1)=J_UPDATE(2)                      !<-- Ending J in nest task's 1st region updated by parent.
                    SEND_TASK(KP)%ISTART(2)=NINT(ITS_PARENT_ON_CHILD(ID_NE)) !<-- Starting I in nest task's 2nd region updated by parent.
                    SEND_TASK(KP)%IEND  (2)=I_UPDATE(2)                      !<-- Ending I in nest task's 2nd region updated by parent.
                    SEND_TASK(KP)%JSTART(2)=J_UPDATE(3)                      !<-- Starting J in nest task's 2nd region updated by parent.
                    SEND_TASK(KP)%JEND  (2)=J_END_X                          !<-- Ending J in nest task's 2nd region updated by parent.
!
                  ENDIF
!
                ENDIF sw_north
!
!----------------------------------------------
!***  Is there a parent task east of the first
!***  that provides update data?
!----------------------------------------------
!
                sw_east: IF(ITE_PARENT_ON_CHILD(ID_1)+EPS<I_END_X)THEN      !<-- If true there is another parent task east of the first.
!
                  KOUNT_PARENT_TASKS=KOUNT_PARENT_TASKS+1                   !<-- Increment the parent task counter.
                  KP=KOUNT_PARENT_TASKS
                  SEND_TASK(KP)%ID=ID_1+1                                   !<-- Store the ID of this parent task to the east.
                  ID_E=SEND_TASK(KP)%ID
                  SEND_TASK(KP)%ISTART(1)=NINT(ITS_PARENT_ON_CHILD(ID_E))   !<-- Starting I where parent task ID_E updates nest task.
!
                  IF(JTE_PARENT_ON_CHILD(ID_E)-EPS<=J_UPDATE(2))THEN        !<-- 1st scenario of parent update task to east of first.
!
                    SEND_TASK(KP)%IEND(1)=I_END_X
                    SEND_TASK(KP)%JSTART(1)=J_START_X
                    SEND_TASK(KP)%JEND(1)=NINT(JTE_PARENT_ON_CHILD(ID_E))
!
                  ELSEIF(JTE_PARENT_ON_CHILD(ID_E)+EPS>=J_UPDATE(3))THEN 
!
                    IF(ITS_PARENT_ON_CHILD(ID_E)+EPS>=I_UPDATE(3))THEN      !<-- 2nd scenario of parent update task to east of first.
!
!                                         |
!                                         |
!                                         |
!                                         |
!                 + + + + + + + + + + + +.|     footprint of nest domain
!                 +                      .|     in its pre-move location
!                 +   parent task 1's    .|
!                 +   2nd update region  .|
!                 +                      .|
!                 +.......................----------------------------------
!                 +..............................       +
!                 +                             '       +
!                 +                             '       +
!                 +     parent task 1's         '       +        
!                 +     1st upate region        '    <------- parent task 2's update region
!                 +                             '       +
!                 +                             '       +
!                 +                             '       +
!                 + + + + + + + + + + + + + + + ' + + + +
!                                 ^             '
!                                 |             '
!                             nest task         '<--- parent task boundary
!                             boundary          '
!                             after move        '
!                                               '
!
                      SEND_TASK(KP)%IEND(1)=I_END_X
                      SEND_TASK(KP)%JSTART(1)=J_START_X
                      SEND_TASK(KP)%JEND(1)=J_UPDATE(2)                       
!
                    ELSEIF(ITS_PARENT_ON_CHILD(ID_E)-EPS<=I_UPDATE(2))THEN  !<-- 3rd scenario of parent update task to east of first.
!
!                                                           |
!                                   '                       |
!                                   '                       |
!                                   '                       |
!                         + + + + + ' + + + + + + + + + + +.|
!    parent task 2's      +         '                      .|
!    update region  ---------->     '   parent task 3's    .|     footprint of nest domain
!                         +         '    update region     .|     in its pre-move location
!                         +         '                      .|
!           ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' '
!                  ^      +         '                      .|
!      parent task |      +         '   parent task 4's    .|
!         boundary        +         '  2nd update region   .|
!                         +         '                      .|
!                         +         '.......................----------------------------------
!    parent task 1's      +         '.................................+
!    update region  ---------->     '                                 +
!                         +         '                                 +
!                         +         '         parent task 4's         +
!                         +         '        1st update region        +
!                         +         '                                 + <--- nest task boundary after move
!                         +         '                                 +
!                         +         '                                 +
!                         + + + + + ' + + + + + + + + + + + + + + + + +
!                                   '
!                                   '
!                                   '<----- parent task boundary
!                                   '
!                                   '
!
!
                      SEND_TASK(KP)%IEND  (1)=I_END_X                          !<-- Ending I of 1st update region in nest task by parent.
                      SEND_TASK(KP)%JSTART(1)=J_START_X                        !<-- Starting J of 1st update region in nest task by parent.
                      SEND_TASK(KP)%JEND  (1)=J_UPDATE(2)                      !<-- Ending J of 1st update region in nest task by parent.
                      SEND_TASK(KP)%ISTART(2)=NINT(ITS_PARENT_ON_CHILD(ID_E))  !<-- Starting I of 2nd update region in nest task by parent.
                      SEND_TASK(KP)%IEND  (2)=I_UPDATE(2)                      !<-- Ending I of 2nd update region in nest task by parent.
                      SEND_TASK(KP)%JSTART(2)=J_UPDATE(3)                      !<-- Starting J of 2nd update region in nest task by parent.
                      SEND_TASK(KP)%JEND  (2)=MIN(J_END_X                    & !<-- Ending J of 2nd update region in nest task by parent.
                                             ,NINT(JTE_PARENT_ON_CHILD(ID_E)))
                    ENDIF
!
                  ENDIF 
!
                ENDIF sw_east
!
                EXIT search_i
!
              ENDIF sw
!
!------------------------------------------------------
!------------------------------------------------------
!***  The nest task on the SE corner of the footprint.
!------------------------------------------------------
!------------------------------------------------------
!
              se: IF(CORNER=='SE')THEN                                      !<-- This nest task lies on the SE corner of the footprint.
!
                SEND_TASK(1)%ISTART(1)=I_START_X                            !<-- Nest I where parent task 1 begins updating nest task.
                SEND_TASK(1)%IEND(1)=MIN(I_UPDATE(4)                     &  !<-- Nest I where parent task 1 ends updating nest task.
                                    ,NINT(ITE_PARENT_ON_CHILD(ID_1))) 
!
                SEND_TASK(1)%JSTART(1)=J_START_X                            !<-- Nest J where parent task 1 begins updating nest task.
!
                IF(JTE_PARENT_ON_CHILD(ID_1)-EPS<=J_UPDATE(2))THEN
                  SEND_TASK(1)%JEND(1)=NINT(JTE_PARENT_ON_CHILD(ID_1))      !<-- Nest J where parent task 1 ends updating nest task.
!
                ELSEIF(JTE_PARENT_ON_CHILD(ID_1)+EPS>=J_UPDATE(3))THEN
                  SEND_TASK(1)%JEND(1)=J_UPDATE(2)                          !<-- Nest J where parent task 1 ends updating nest task.
!
                  IF(ITE_PARENT_ON_CHILD(ID_1)+EPS>=I_UPDATE(3))THEN        !<-- Parent task 1 covers SE corner of footprint too.
!
                    SEND_TASK(1)%ISTART(2)=I_UPDATE(3)                      !<-- 2nd region on nest task updated by parent task 1
                    SEND_TASK(1)%IEND(2)=MIN(I_END_X                     &  !
                                        ,NINT(ITE_PARENT_ON_CHILD(ID_1)))   !
                    SEND_TASK(1)%JSTART(2)=J_UPDATE(3)                      !
                    SEND_TASK(1)%JEND(2)=MIN(J_END_X                     &  !
                                        ,NINT(JTE_PARENT_ON_CHILD(ID_1)))   !<--
!
                  ENDIF
                ENDIF
!
!-----------------------------------------------------------------------
!***  The points in this nest task subdomain being updated by the first
!***  identified parent task have been demarcated.  Now identify any
!***  other parent tasks updating this nest task lying on the SE corner
!***  of the footprint.
!-----------------------------------------------------------------------
!                  
!------------------------------------------------------
!***  Is there a parent task to the north of the first
!***  that provides update data?
!------------------------------------------------------
!
                se_north: IF(JTE_PARENT_ON_CHILD(ID_1)+EPS<J_END_X          &  !<-- If true there is another parent task north of the first
                                      .AND.                                 &  !    on the SE corner of the footprint or only on its
                             ITE_PARENT_ON_CHILD(ID_1)+EPS>=I_UPDATE(3))THEN   !<-- east side.
!
                  KOUNT_PARENT_TASKS=KOUNT_PARENT_TASKS+1                    !<-- Increment parent task counter.
                  KP=KOUNT_PARENT_TASKS
                  SEND_TASK(KP)%ID=ID_1+INPES_PARENT                         !<-- Store the ID of this parent task to the north.
                  ID_N=SEND_TASK(KP)%ID 
!
                  IF(JTS_PARENT_ON_CHILD(ID_N)+EPS>=J_UPDATE(3))THEN         !<-- Parent task ID_N does not cover SE corner of footprint.
                    SEND_TASK(KP)%ISTART(1)=I_UPDATE(3)                      !<-- Starting I on nest task where parent task ID_N updates.
                    SEND_TASK(KP)%IEND(1)=MIN(I_UPDATE(4)                 &  !<-- Ending I on nest task where parent task ID_N updates.
                                         ,NINT(ITE_PARENT_ON_CHILD(ID_N)))   !<--
                    SEND_TASK(KP)%JSTART(1)=NINT(JTS_PARENT_ON_CHILD(ID_N))  !<-- Starting J on nest task where parent task ID_N updates.
                    SEND_TASK(KP)%JEND(1)=J_END_X                            !<-- Ending J on nest task where parent task ID_N updates.
!
                  ELSEIF(JTS_PARENT_ON_CHILD(ID_N)-EPS<=J_UPDATE(2))THEN     !<-- Parent task ID_N covers SE corner of footprint too.
!
                    SEND_TASK(KP)%ISTART(1)=I_START_X                        !<-- Starting I on nest task where parent task ID_N updates
                    SEND_TASK(KP)%IEND(1)=MIN(I_END_X                     &  !<-- Ending I on nest task where parent task ID_N updates
                                         ,NINT(ITE_PARENT_ON_CHILD(ID_N)))   !<-- in nest task's 1st region.
                    SEND_TASK(KP)%JSTART(1)=NINT(JTS_PARENT_ON_CHILD(ID_N))  !<-- Starting J of parent task ID_N in nest task 1st region.
                    SEND_TASK(KP)%JEND  (1)=J_UPDATE(2)                      !<-- Ending J of parent task ID_N in nest task 1st region.
                    SEND_TASK(KP)%ISTART(2)=I_UPDATE(3)                      !<-- Starting I of parent task ID_N in nest task 2nd region.
                    SEND_TASK(KP)%IEND  (2)=MIN(I_UPDATE(4)               &  !<-- Ending I of parent task ID_N in nest task 2nd region.
                                           ,NINT(ITE_PARENT_ON_CHILD(ID_N))) !<--
                    SEND_TASK(KP)%JSTART(2)=J_UPDATE(3)                      !<-- Starting J of parent task ID_N in nest task 2nd region.
                    SEND_TASK(KP)%JEND  (2)=J_END_X                          !<-- Ending J of parent task ID_N in nest task 2nd region.
!
                  ENDIF
!
                ELSEIF(JTE_PARENT_ON_CHILD(ID_1)+EPS<J_UPDATE(2))THEN        !<-- Parent task ID_N only on south side of footprint.
!
                  KOUNT_PARENT_TASKS=KOUNT_PARENT_TASKS+1                    !<-- Increment parent task counter.
                  KP=KOUNT_PARENT_TASKS
                  SEND_TASK(KP)%ID=ID_1+INPES_PARENT                         !<-- Store the ID of this parent task to the north.
                  ID_N=SEND_TASK(KP)%ID
!
                  SEND_TASK(KP)%ISTART(1)=I_START_X                          !<-- Starting I on nest task where parent task ID_N updates
                  SEND_TASK(KP)%IEND  (1)=MIN(I_UPDATE(2)                 &  !<-- Ending I on nest task where parent task ID_N updates
                                         ,NINT(ITE_PARENT_ON_CHILD(ID_N)))
                  SEND_TASK(KP)%JSTART(1)=NINT(JTS_PARENT_ON_CHILD(ID_N))    !<-- Starting J on nest task where parent task ID_N updates.
                  SEND_TASK(KP)%JEND  (1)=J_UPDATE(2)                        !<-- Ending J on nest task where parent task ID_N updates.
!
                ENDIF se_north
!
!-----------------------------------------------
!***  Does a parent task northeast of the first 
!***  provide any update data?  The presence of
!***  a parent task to the north of the first 
!***  that is providing update data is not 
!***  required for a parent task to exist to
!***  the northeast that sends update data
!***  to this nest task.  That is because a
!***  parent task to the north might be totally
!***  covered by the footprint's SE corner.
!-----------------------------------------------
!
                IF(ITE_PARENT_ON_CHILD(ID_1)+EPS<I_END_X                &  !<-- 1st scenario of parent update task to northeast of 
                                .AND.                                   &  !    the first parent update task.  The NE parent task
                   ITE_PARENT_ON_CHILD(ID_1)+EPS>=I_UPDATE(2)           &  !    does not cover the SE corner of the footprint.
                                .AND.                                   &  !
                   JTE_PARENT_ON_CHILD(ID_1)+EPS<J_END_X)THEN              !<--
!
                  KOUNT_PARENT_TASKS=KOUNT_PARENT_TASKS+1                  !<-- Increment parent task counter.
                  KP=KOUNT_PARENT_TASKS
                  SEND_TASK(KP)%ID=ID_1+INPES_PARENT+1                     !<-- Store the ID of this parent task northeast of the first
                  ID_NE=SEND_TASK(KP)%ID
                  SEND_TASK(KP)%ISTART(1)=NINT(ITS_PARENT_ON_CHILD(ID_NE)) !<-- Starting I where parent task ID_NE updates nest task.
                  SEND_TASK(KP)%IEND  (1)=I_END_X                          !<-- Ending I where parent task ID_NE updates nest task.
                  SEND_TASK(KP)%JSTART(1)=NINT(JTS_PARENT_ON_CHILD(ID_NE)) !<-- Starting J where parent task ID_NE updates nest task.
                  SEND_TASK(KP)%JEND  (1)=J_END_X                          !<-- Ending J where parent task ID_NE updates nest task.
!
                ENDIF
!
                IF(ITE_PARENT_ON_CHILD(ID_1)+EPS<I_UPDATE(2)            &  !<-- 2nd scenario of parent update task to northeast of
                                .AND.                                   &  !    the first update parent task.  The NE parent task
                   JTE_PARENT_ON_CHILD(ID_1)+EPS<J_UPDATE(2))THEN          !<-- covers the SE corner of the footprint.
!
                  KOUNT_PARENT_TASKS=KOUNT_PARENT_TASKS+1                  !<-- Increment parent task counter.
                  KP=KOUNT_PARENT_TASKS
                  SEND_TASK(KP)%ID=ID_1+INPES_PARENT+1                     !<-- Store the ID of this parent task northeast of the first.
                  ID_NE=SEND_TASK(KP)%ID
                  SEND_TASK(KP)%ISTART(1)=NINT(ITS_PARENT_ON_CHILD(ID_NE)) !<-- Starting I in nest task's 1st region updated by parent.
                  SEND_TASK(KP)%IEND  (1)=I_END_X                          !<-- Ending I in nest task's 1st region updated by parent.
                  SEND_TASK(KP)%JSTART(1)=NINT(JTS_PARENT_ON_CHILD(ID_NE)) !<-- Starting J in nest task's 1st region updated by parent.
                  SEND_TASK(KP)%JEND  (1)=J_UPDATE(2)                      !<-- Ending J in nest task's 1st region updated by parent.
                  SEND_TASK(KP)%ISTART(2)=I_UPDATE(3)                      !<-- Starting I in nest task's 2nd region updated by parent.
                  SEND_TASK(KP)%IEND  (2)=I_END_X                          !<-- Ending I in nest task's 2nd region updated by parent.
                  SEND_TASK(KP)%JSTART(2)=J_UPDATE(3)                      !<-- Starting J in nest task's 2nd region updated by parent.
                  SEND_TASK(KP)%JEND  (2)=J_END_X                          !<-- Ending J in nest task's 2nd region updated by parent.
!
                ENDIF
!
                IF(ITE_PARENT_ON_CHILD(ID_1)+EPS<I_UPDATE(2)            &  !<-- 3rd scenario of parent update task to northeast of
                                .AND.                                   &  !    the first update parent task.  The NE parent task
                   JTE_PARENT_ON_CHILD(ID_1)+EPS>=J_UPDATE(2)           &  !<-- is north of the SE corner of the footprint.
                                .AND.                                   &  !
                   JTE_PARENT_ON_CHILD(ID_1)+EPS<J_END_X)THEN    
!
                  KOUNT_PARENT_TASKS=KOUNT_PARENT_TASKS+1                  !<-- Increment parent task counter.
                  KP=KOUNT_PARENT_TASKS
                  SEND_TASK(KP)%ID=ID_1+INPES_PARENT+1                     !<-- Store the ID of this parent task northeast of the first.
                  ID_NE=SEND_TASK(KP)%ID
                  SEND_TASK(KP)%ISTART(1)=I_UPDATE(3)                      !<-- Starting I in nest task's 1st region updated by parent.
                  SEND_TASK(KP)%IEND  (1)=I_END_X                          !<-- Ending I in nest task's 1st region updated by parent.
                  SEND_TASK(KP)%JSTART(1)=NINT(JTS_PARENT_ON_CHILD(ID_NE)) !<-- Starting J in nest task's 1st region updated by parent.
                  SEND_TASK(KP)%JEND  (1)=J_END_X                          !<-- Ending J in nest task's 1st region updated by parent.
!
                ENDIF
!
!----------------------------------------------
!***  Is there a parent task east of the first
!***  that provides update data?
!----------------------------------------------
!
                se_east: IF(ITE_PARENT_ON_CHILD(ID_1)+EPS<I_END_X)THEN     !<-- If true there is another parent task east of the first.
!
                  KOUNT_PARENT_TASKS=KOUNT_PARENT_TASKS+1                  !<-- Increment the parent task counter.
                  KP=KOUNT_PARENT_TASKS
                  SEND_TASK(KP)%ID=ID_1+1                                  !<-- Store the ID of this parent task to the east.
                  ID_E=SEND_TASK(KP)%ID
                  SEND_TASK(KP)%ISTART(1)=NINT(ITS_PARENT_ON_CHILD(ID_E))  !<-- Starting I where parent task ID_E updates nest task.
                  SEND_TASK(KP)%IEND(1)=I_END_X                            !<-- Ending I where parent task ID_E updates nest task.
                  SEND_TASK(KP)%JSTART(1)=J_START_X                        !<-- Starting J where parent task ID_E updates nest task.
!
                  IF(JTE_PARENT_ON_CHILD(ID_E)-EPS<=J_UPDATE(2))THEN       !<-- 1st scenario of parent update task E of first.  No corner.
!
                    SEND_TASK(KP)%JEND(1)=NINT(JTE_PARENT_ON_CHILD(ID_E))
!
                  ELSEIF(ITS_PARENT_ON_CHILD(ID_E)+EPS>=I_UPDATE(3))THEN     !<-- 2nd scenario of parent update task E of first.  No corner.
!
                    SEND_TASK(KP)%JEND(1)=MIN(J_END_X                     &  !<-- Ending J of 2nd update region in nest task by parent.
                                         ,NINT(JTE_PARENT_ON_CHILD(ID_E)))   !<--
!
                  ELSEIF(ITS_PARENT_ON_CHILD(ID_E)-EPS<=I_UPDATE(2)       &  !<-- 2nd scenario of parent update task to E of first.
                                     .AND.                                &  !    The east parent task covers the SE corner of the
                         JTE_PARENT_ON_CHILD(ID_E)+EPS>=J_UPDATE(3))THEN     !<-- footprint.
!
                    SEND_TASK(KP)%JEND  (1)=J_UPDATE(2)                      !<-- Ending J of 1st update region in nest task by parent.
                    SEND_TASK(KP)%ISTART(2)=I_UPDATE(3)                      !<-- Starting I of 2nd update region in nest task by parent.
                    SEND_TASK(KP)%IEND  (2)=I_END_X                          !<-- Ending I of 2nd update region in nest task by parent.
                    SEND_TASK(KP)%JSTART(2)=J_UPDATE(3)                      !<-- Starting J of 2nd update region in nest task by parent.
                    SEND_TASK(KP)%JEND  (2)=MIN(J_END_X                   &  !<-- Ending J of 2nd update region in nest task by parent.
                                           ,NINT(JTE_PARENT_ON_CHILD(ID_E))) !<--
                  ENDIF                                                      !<-- Finished with parent task east of the first one.
!
                ENDIF se_east
!
                EXIT search_i
!
              ENDIF se
!
!-----------------------------------------------------------------------
!
!------------------------------------------------------
!------------------------------------------------------
!***  The nest task on the NW corner of the footprint.
!------------------------------------------------------
!------------------------------------------------------
!
              nw: IF(CORNER=='NW')THEN                                     !<-- This nest task lies on the NW corner of the footprint.
!
                SEND_TASK(1)%ISTART(1)=I_START_X                           !<-- Nest I where parent task 1 begins updating nest task.
                SEND_TASK(1)%JSTART(1)=J_START_X                           !<-- Nest J where parent task 1 begins updating nest task.
!
                IF(ITE_PARENT_ON_CHILD(ID_1)-EPS<=I_UPDATE(2))THEN
                  SEND_TASK(1)%IEND(1)=NINT(ITE_PARENT_ON_CHILD(ID_1))     !<-- Nest I where parent task 1 ends updating nest task.
                  SEND_TASK(1)%JEND(1)=MIN(J_END_X                      &  !<-- Nest J where parent task 1 ends updating nest task.
                                      ,NINT(JTE_PARENT_ON_CHILD(ID_1)))    !<--
!
                ELSEIF(ITE_PARENT_ON_CHILD(ID_1)+EPS>=I_UPDATE(3))THEN
                  SEND_TASK(1)%IEND(1)=I_UPDATE(2)                         !<-- Nest J where parent task 1 ends updating nest task.
                  SEND_TASK(1)%JEND(1)=MIN(J_UPDATE(2)                  &  !<-- Nest J where parent task 1 ends updating nest task.
                                      ,NINT(JTE_PARENT_ON_CHILD(ID_1)))    !<--
!
                  IF(JTE_PARENT_ON_CHILD(ID_1)+EPS>=J_UPDATE(3))THEN       !<-- Parent task 1 covers NW corner of footprint too.
                    SEND_TASK(1)%ISTART(2)=I_START_X                       !<-- 2nd region on nest task updated by parent task 1
                    SEND_TASK(1)%IEND(2)=MIN(I_END_X                    &  !
                                        ,NINT(ITE_PARENT_ON_CHILD(ID_1)))  !
                    SEND_TASK(1)%JSTART(2)=J_UPDATE(3)                     !
                    SEND_TASK(1)%JEND(2)=MIN(J_END_X                    &  !
                                        ,NINT(JTE_PARENT_ON_CHILD(ID_1)))  !<--
                  ENDIF
                ENDIF
!
!-----------------------------------------------------------------------
!***  The points in this nest task subdomain being updated by the first
!***  identified parent task have been demarcated.  Now identify any
!***  other parent tasks updating this nest task lying on the NW corner
!***  of the footprint.
!-----------------------------------------------------------------------
!                  
!------------------------------------------------------
!***  Is there a parent task to the north of the first
!***  that provides update data?
!------------------------------------------------------
!
                nw_north: IF(JTE_PARENT_ON_CHILD(ID_1)+EPS<J_END_X)THEN      !<-- If true there is another parent task north of the first.
!
                  KOUNT_PARENT_TASKS=KOUNT_PARENT_TASKS+1                    !<-- Increment parent task counter.
                  KP=KOUNT_PARENT_TASKS
                  SEND_TASK(KP)%ID=ID_1+INPES_PARENT                         !<-- Store the ID of this parent task to the north.
                  ID_N=SEND_TASK(KP)%ID 
                  SEND_TASK(KP)%ISTART(1)=I_START_X                          !<-- Starting I on nest task where parent task ID_N updates.
                  SEND_TASK(KP)%JSTART(1)=NINT(JTS_PARENT_ON_CHILD(ID_N))    !<-- Starting J on nest task where parent task ID_N updates.
!
                  IF(JTS_PARENT_ON_CHILD(ID_N)+EPS>=J_UPDATE(3)           &  !<-- Parent task ID_N does not cover NW corner of footprint.
                                     .OR.                                 &  !
                     ITE_PARENT_ON_CHILD(ID_N)-EPS<=I_UPDATE(2))THEN         !<--
!
                    SEND_TASK(KP)%IEND(1)=MIN(I_UPDATE(4)                 &  !<-- Ending I on nest task where parent task ID_N updates.
                                         ,NINT(ITE_PARENT_ON_CHILD(ID_N)))   !<--
                    SEND_TASK(KP)%JEND(1)=J_END_X                            !<-- Ending J on nest task where parent task ID_N updates.
!
                  ELSEIF(JTS_PARENT_ON_CHILD(ID_N)-EPS<=J_UPDATE(2)       &  !<-- Parent task ID_N covers NW corner of footprint too.
                                     .AND.                                &
                         ITE_PARENT_ON_CHILD(ID_N)+EPS>=I_UPDATE(3))THEN
!
                    SEND_TASK(KP)%IEND(1)=MIN(I_UPDATE(2)                 &  !<-- Ending I on nest task where parent task ID_N updates
                                         ,NINT(ITE_PARENT_ON_CHILD(ID_N)))   !<-- in nest task's 1st region.
                    SEND_TASK(KP)%JEND  (1)=J_UPDATE(2)                      !<-- Ending J of parent task ID_N in nest task 1st region.
                    SEND_TASK(KP)%ISTART(2)=I_START_X                        !<-- Starting I of parent task ID_N in nest task 2nd region.
                    SEND_TASK(KP)%IEND  (2)=MIN(I_UPDATE(4)               &  !<-- Ending I of parent task ID_N in nest task 2nd region.
                                           ,NINT(ITE_PARENT_ON_CHILD(ID_N))) !<--
                    SEND_TASK(KP)%JSTART(2)=J_UPDATE(3)                      !<-- Starting J of parent task ID_N in nest task 2nd region.
                    SEND_TASK(KP)%JEND  (2)=J_END_X                          !<-- Ending J of parent task ID_N in nest task 2nd region.
                  ENDIF                                                      !<-- End contingencies of parent task north of first one.
!
!-----------------------------------------------
!***  Does a parent task northeast of the first 
!***  provide any update data?  For this to be
!***  true there must be a parent task to the
!***  north of the first so we remain in the
!***  nw_north IF block.
!-----------------------------------------------
!
                  IF(ITE_PARENT_ON_CHILD(ID_N)+EPS<I_UPDATE(2)             &  !<-- 1st scenario of parent update task to northeast of
                                  .AND.                                    &  !    the first parent update task. It is on NW corner.
                     JTS_PARENT_ON_CHILD(ID_N)-EPS<=J_UPDATE(2))THEN          !<-- It is a 'flag-shape' scenario.
!
                    KOUNT_PARENT_TASKS=KOUNT_PARENT_TASKS+1                   !<-- Increment parent task counter.
                    KP=KOUNT_PARENT_TASKS
                    SEND_TASK(KP)%ID=ID_N+1                                   !<-- Store the ID of this parent task northeast of the first
                    ID_NE=SEND_TASK(KP)%ID
                    SEND_TASK(KP)%ISTART(1)=NINT(ITS_PARENT_ON_CHILD(ID_NE))  !<-- Starting I where parent task ID_NE updates nest task.
                    SEND_TASK(KP)%IEND  (1)=I_UPDATE(2)                       !<-- Ending I where parent task ID_NE updates nest task.
                    SEND_TASK(KP)%JSTART(1)=NINT(JTS_PARENT_ON_CHILD(ID_NE))  !<-- Starting J where parent task ID_NE updates nest task.
                    SEND_TASK(KP)%JEND  (1)=J_UPDATE(2)                       !<-- Ending J where parent task ID_NE updates nest task.
                    SEND_TASK(KP)%ISTART(2)=SEND_TASK(KP)%ISTART(1)           !<-- Starting I where parent task ID_NE updates 2nd region.
                    SEND_TASK(KP)%IEND  (2)=I_END_X                           !<-- Ending I where parent task ID_NE updates 2nd region.
                    SEND_TASK(KP)%JSTART(2)=J_UPDATE(3)                       !<-- Starting J where parent task ID_NE updates 2nd region.
                    SEND_TASK(KP)%JEND  (2)=J_END_X                           !<-- Ending J where parent task ID_NE updates 2nd region.
!
                  ENDIF
!
                  IF(ITE_PARENT_ON_CHILD(ID_N)+EPS< I_END_X                &  !<-- 2nd scenario of parent update task to northeast of
                                  .AND.                                    &  !    the first parent update task.  A simple rectangle.
                    (JTS_PARENT_ON_CHILD(ID_N)+EPS>=J_UPDATE(3)            &  !
                                  .OR.                                     &  !
                     ITE_PARENT_ON_CHILD(ID_N)+EPS>=I_UPDATE(2)            &  !
                                  .AND.                                    &  !
                     JTS_PARENT_ON_CHILD(ID_N)-EPS<=J_UPDATE(2)))THEN         !<--
!
                    KOUNT_PARENT_TASKS=KOUNT_PARENT_TASKS+1                   !<-- Increment parent task counter.
                    KP=KOUNT_PARENT_TASKS
                    SEND_TASK(KP)%ID=ID_N+1                                   !<-- Store the ID of this parent task northeast of the first.
                    ID_NE=SEND_TASK(KP)%ID
                    SEND_TASK(KP)%ISTART(1)=NINT(ITS_PARENT_ON_CHILD(ID_NE))  !<-- Starting I in nest task's update region.
                    SEND_TASK(KP)%IEND  (1)=I_END_X                           !<-- Ending I in nest task's update region.
                    SEND_TASK(KP)%JSTART(1)=MAX(J_UPDATE(3)                &  !<-- Starting J where parent task ID_NE updates nest.
                                           ,NINT(JTS_PARENT_ON_CHILD(ID_NE)))
                    SEND_TASK(KP)%JEND  (1)=J_END_X                           !<-- Ending J in nest task's update region.
!
                  ENDIF
!
                ENDIF nw_north
!
!----------------------------------------------
!***  Is there a parent task east of the first
!***  that provides update data?
!----------------------------------------------
!
                nw_east: IF(ITE_PARENT_ON_CHILD(ID_1)+EPS<I_END_X)THEN         !<-- Necessary, not sufficient for parent task east of first.
!
                  IF(ITE_PARENT_ON_CHILD(ID_1)+EPS<I_UPDATE(2))THEN            !<-- 1st scenario of parent task east of the first.
!
                    KOUNT_PARENT_TASKS=KOUNT_PARENT_TASKS+1                    !<-- Increment the parent task counter.
                    KP=KOUNT_PARENT_TASKS
                    SEND_TASK(KP)%ID=ID_1+1                                    !<-- Store the ID of this parent task to the east.
                    ID_E=SEND_TASK(KP)%ID
                    SEND_TASK(KP)%ISTART(1)=NINT(ITS_PARENT_ON_CHILD(ID_E))    !<-- Starting I where parent task ID_E updates nest task.
                    SEND_TASK(KP)%IEND  (1)=I_UPDATE(2)                        !<-- Ending I where parent task ID_E updates nest task.
                    SEND_TASK(KP)%JSTART(1)=J_START_X                          !<-- Starting J where parent task ID_E updates nest task.
                    SEND_TASK(KP)%JEND  (1)=MIN(J_UPDATE(2)                 &  !<-- Ending J where parent task ID_E updates nest task.
                                           ,NINT(JTE_PARENT_ON_CHILD(ID_E)))   !
!
                    IF(JTE_PARENT_ON_CHILD(ID_E)+EPS>=J_UPDATE(3))THEN
                      SEND_TASK(KP)%ISTART(2)=NINT(ITS_PARENT_ON_CHILD(ID_E))  !<-- Starting I for east parent in 2nd update region.
                      SEND_TASK(KP)%IEND  (2)=I_END_X                          !<-- Ending I for east parent in 2nd update region.
                      SEND_TASK(KP)%JSTART(2)=J_UPDATE(3)                      !<-- Starting J for east parent in 2nd update region.
                      SEND_TASK(KP)%JEND  (2)=MIN(J_END_X                   &  !<-- Ending J for east parent in 2nd update region.
                                             ,NINT(JTE_PARENT_ON_CHILD(ID_E)))
                    ENDIF
!
                  ENDIF
!
                  IF(ITE_PARENT_ON_CHILD(ID_1)+EPS>=I_UPDATE(2)             &  !<-- 2nd scenario of a parent task to the east of
                                    .AND.                                   &  !    the first parent task.
                     JTE_PARENT_ON_CHILD(ID_1)+EPS>=J_UPDATE(3))THEN           !<--
!
                    KOUNT_PARENT_TASKS=KOUNT_PARENT_TASKS+1                    !<-- Increment the parent task counter.
                    KP=KOUNT_PARENT_TASKS
                    SEND_TASK(KP)%ID=ID_1+1                                    !<-- Store the ID of this parent task to the east.
                    ID_E=SEND_TASK(KP)%ID
                    SEND_TASK(KP)%ISTART(1)=NINT(ITS_PARENT_ON_CHILD(ID_E))    !<-- Starting I where parent task ID_E updates nest task.
                    SEND_TASK(KP)%IEND  (1)=I_END_X                            !<-- Ending I where parent task ID_E updates nest task.
                    SEND_TASK(KP)%JSTART(1)=J_UPDATE(3)                        !<-- Starting J where parent task ID_E updates nest task.
                    SEND_TASK(KP)%JEND  (1)=MIN(J_END_X                     &  !<-- Ending J where parent task ID_E updates nest task.
                                           ,NINT(JTE_PARENT_ON_CHILD(ID_E))) 
                  ENDIF
!
                ENDIF nw_east
!
                EXIT search_i
!
              ENDIF nw
!
!------------------------------------------------------
!------------------------------------------------------
!***  The nest task on the NE corner of the footprint.
!------------------------------------------------------
!------------------------------------------------------
!
              ne: IF(CORNER=='NE')THEN                                     !<-- This nest task lies on the NE corner of the footprint.
!
!-----------------------------------------------------------------------
!***  The northeast corner of the footprint is even more involved
!***  than the other three because [I_UPDATE(1),J_UPDATE(1)] on this
!***  nest task is within the footprint of the nest's previous location
!***  and thus is not updated by parent task 1.  In fact parent task 1
!***  might not update any points on this nest task if that region of
!***  the nest task covered by parent task 1 lies entirely within the
!***  footprint.  
!-----------------------------------------------------------------------
!
                KOUNT_PARENT_TASKS=0                                       !<-- Parent task ID_1 might not send any data.
!
                IF(ITE_PARENT_ON_CHILD(ID_1)+EPS>=I_UPDATE(3))THEN               
!
                  KOUNT_PARENT_TASKS=1                                     !<-- Parent task ID_1 does send data.
                  SEND_TASK(1)%ISTART(1)=I_UPDATE(3)                       !<-- Nest I where parent task 1 begins updating nest task.
                  SEND_TASK(1)%IEND(1)=MIN(I_END_X                      &  !<-- Nest I where parent task 1 ends updating nest task.
                                      ,NINT(ITE_PARENT_ON_CHILD(ID_1))) 
                  SEND_TASK(1)%JSTART(1)=J_START_X                         !<-- Nest J where parent task 1 begins updating nest task.
                  SEND_TASK(1)%JEND  (1)=MIN(J_UPDATE(2)                &  !<-- Nest J where parent task 1 ends updating nest task.
                                        ,NINT(JTE_PARENT_ON_CHILD(ID_1)))
!
                  IF(JTE_PARENT_ON_CHILD(ID_1)+EPS>=J_UPDATE(3))THEN
                    SEND_TASK(1)%ISTART(2)=I_START_X                       !<-- Nest I where parent task 1 begins updating 2nd region.
                    SEND_TASK(1)%IEND  (2)=SEND_TASK(1)%IEND(1)            !<-- Nest I where parent task 1 ends updating 2nd region.
                    SEND_TASK(1)%JSTART(2)=J_UPDATE(3)                     !<-- Nest J where parent task 1 begins updating 2nd region.
                    SEND_TASK(1)%JEND  (2)=MIN(J_END_X                  &  !<-- Nest J where parent task 1 ends updating 2nd region.
                                          ,NINT(JTE_PARENT_ON_CHILD(ID_1)))
                  ENDIF
!
                ELSEIF(JTE_PARENT_ON_CHILD(ID_1)+EPS>=J_UPDATE(3))THEN
!
                  KOUNT_PARENT_TASKS=1                                     !<-- Parent task ID_1 does send data.
                  SEND_TASK(1)%ISTART(1)=I_UPDATE(1)                       !<-- Nest I where parent task 1 begins updating nest task.
                  SEND_TASK(1)%IEND(1)=NINT(ITE_PARENT_ON_CHILD(ID_1))     !<-- Nest I where parent task 1 ends updating nest task.
                  SEND_TASK(1)%JSTART(1)=J_UPDATE(3)                       !<-- Nest J where parent task 1 begins updating nest task.
                  SEND_TASK(1)%JEND(1)=MIN(J_END_X                      &  !<-- Nest J where parent task 1 ends updating 2nd region.
                                      ,NINT(JTE_PARENT_ON_CHILD(ID_1)))
!
                ENDIF
!
!-----------------------------------------------------------------------
!***  The points in this nest task subdomain being updated by the first
!***  identified parent task have been demarcated.  Now identify any
!***  other parent tasks updating this nest task lying on the NE corner
!***  of the footprint.
!-----------------------------------------------------------------------
!                  
!------------------------------------------------------
!***  Is there a parent task to the north of the first
!***  that provides update data?  
!------------------------------------------------------
!
                ne_north: IF(JTE_PARENT_ON_CHILD(ID_1)+EPS<J_END_X)THEN     !<-- If true there is another parent task north of the first.
!
                  KOUNT_PARENT_TASKS=KOUNT_PARENT_TASKS+1                   !<-- Increment parent task counter.
                  KP=KOUNT_PARENT_TASKS
                  SEND_TASK(KP)%ID=ID_1+INPES_PARENT                        !<-- Store the ID of this parent task to the north.
                  ID_N=SEND_TASK(KP)%ID
!
                  IF(ITE_PARENT_ON_CHILD(ID_N)<=I_UPDATE(2))THEN            !<-- 1st scenario of parent task north of the first.
                    SEND_TASK(KP)%ISTART(1)=I_START_X                       !<-- Nest I where parent task ID_N begins updating this region.
                    SEND_TASK(KP)%IEND  (1)=NINT(ITE_PARENT_ON_CHILD(ID_N)) !<-- Nest I where parent task ID_N ends updating this region.
                    SEND_TASK(KP)%JSTART(1)=MAX(J_UPDATE(3)              &  !<-- Nest J where parent task ID_N begins updating this region.
                                           ,NINT(JTS_PARENT_ON_CHILD(ID_N)))
                    SEND_TASK(KP)%JEND  (1)=J_END_X                         !<-- Nest J where parent task ID_N ends updating this region.
                  ENDIF
!
                  IF(ITE_PARENT_ON_CHILD(ID_N)>=I_UPDATE(3))THEN 
                    SEND_TASK(KP)%IEND  (1)=MIN(I_END_X                    &  !<-- Nest I where parent task ID_N ends updating this region.
                                           ,NINT(ITE_PARENT_ON_CHILD(ID_N)))
                    IF(JTS_PARENT_ON_CHILD(ID_N)<=J_UPDATE(2))THEN            !<-- 2nd scenario of parent task north of the first.
                      SEND_TASK(KP)%ISTART(1)=I_UPDATE(3)                     !<-- Nest I where parent task ID_N begins updating this region.
                      SEND_TASK(KP)%JSTART(1)=NINT(JTS_PARENT_ON_CHILD(ID_N)) !<-- Nest J where parent task ID_N begins updating this region.
                      SEND_TASK(KP)%JEND  (1)=J_UPDATE(2)                     !<-- Nest J where parent task ID_N ends updating this region.
                      SEND_TASK(KP)%ISTART(2)=I_START_X                       !<-- Nest I where parent task ID_N begins updating 2nd region.
                      SEND_TASK(KP)%IEND  (2)=SEND_TASK(KP)%IEND(1)           !<-- Nest I where parent task ID_N ends updating 2nd region.
                      SEND_TASK(KP)%JSTART(2)=J_UPDATE(3)                     !<-- Nest J where parent task ID_N begins updating 2nd region.
                      SEND_TASK(KP)%JEND  (2)=J_END_X                         !<-- Nest J where parent task ID_N ends updating 2nd region.
                    ELSEIF(JTS_PARENT_ON_CHILD(ID_N)>=J_UPDATE(3))THEN        !<-- 3rd scenario of parent task north of the first.
                      SEND_TASK(KP)%ISTART(1)=I_START_X                       !<-- Nest I where parent task ID_N begins updating this region.
                      SEND_TASK(KP)%JSTART(1)=MAX(J_UPDATE(3)              &  !<-- Nest J where parent task ID_N begins updating this region.
                                             ,NINT(JTS_PARENT_ON_CHILD(ID_N)))
                      SEND_TASK(KP)%JEND  (1)=J_END_X                         !<-- Nest J where parent task ID_N begins updating this region.
                    ENDIF
                  ENDIF
!
!-----------------------------------------------
!***  Does a parent task northeast of the first
!***  provide any update data?  This can only
!***  happen if there was a parent task to
!***  the north of the first therefore we
!***  remain in the ne_north IF block.
!-----------------------------------------------
!
                  ne_ne: IF(ITE_PARENT_ON_CHILD(ID_N)+EPS<I_END_X)THEN         !<-- If so, a parent task NE of the first provides data.
!
                    KOUNT_PARENT_TASKS=KOUNT_PARENT_TASKS+1                    !<-- Increment parent task counter.
                    KP=KOUNT_PARENT_TASKS
                    SEND_TASK(KP)%ID=ID_N+1                                    !<-- Store the ID of this parent task to the northeast.
                    ID_NE=SEND_TASK(KP)%ID
!
                    IF(JTS_PARENT_ON_CHILD(ID_NE)+EPS>=J_UPDATE(3)          &  !<-- 1st scenario of parent task to the NE of the first.
                                      .OR.                                  &  !
                       ITS_PARENT_ON_CHILD(ID_NE)+EPS>=I_UPDATE(3)          &  !
                                      .AND.                                 &  !
                       JTS_PARENT_ON_CHILD(ID_NE)-EPS<=J_UPDATE(2))THEN        !<--
!
                      SEND_TASK(KP)%ISTART(1)=NINT(ITS_PARENT_ON_CHILD(ID_NE)) !<-- Nest I where parent task ID_NE begins updating this region.
                      SEND_TASK(KP)%IEND  (1)=I_END_X                          !<-- Nest I where parent task ID_NE ends updating this region.
                      SEND_TASK(KP)%JSTART(1)=NINT(JTS_PARENT_ON_CHILD(ID_NE)) !<-- Nest J where parent task ID_NE begins updating this region.
                      SEND_TASK(KP)%JEND  (1)=J_END_X                          !<-- Nest J where parent task ID_NE ends updating this region.
                    ENDIF
!
                    IF(ITS_PARENT_ON_CHILD(ID_NE)-EPS<=I_UPDATE(2)          &  !<-- 2nd scenario of parent task to the NE of the first.
                                   .AND.                                    &  !
                       JTS_PARENT_ON_CHILD(ID_NE)-EPS<=J_UPDATE(2))THEN        !<--
!                    
                      SEND_TASK(KP)%ISTART(1)=I_UPDATE(3)                      !<-- Nest I where parent task ID_NE begins updating this region
                      SEND_TASK(KP)%IEND  (1)=I_END_X                          !<-- Nest I where parent task ID_NE ends updating this region.
                      SEND_TASK(KP)%JSTART(1)=NINT(JTS_PARENT_ON_CHILD(ID_NE)) !<-- Nest J where parent task ID_NE begins updating this region
                      SEND_TASK(KP)%JEND  (1)=J_UPDATE(2)                      !<-- Nest J where parent task ID_NE ends updating this region.
                      SEND_TASK(KP)%ISTART(2)=NINT(ITS_PARENT_ON_CHILD(ID_NE)) !<-- Nest I where parent task ID_NE begins updating 2nd region.
                      SEND_TASK(KP)%IEND  (2)=I_END_X                          !<-- Nest I where parent task ID_NE ends updating 2nd region.
                      SEND_TASK(KP)%JSTART(2)=J_UPDATE(3)                      !<-- Nest J where parent task ID_NE begins updating 2nd region.
                      SEND_TASK(KP)%JEND  (2)=J_END_X                          !<-- Nest J where parent task ID_NE ends updating 2nd region.
                    ENDIF
!
                  ENDIF ne_ne
!
                ENDIF ne_north
!
!----------------------------------------------
!***  Is there a parent task east of the first
!***  that provides update data?
!----------------------------------------------
!
                ne_east: IF(ITE_PARENT_ON_CHILD(ID_1)+EPS<I_END_X)THEN      !<-- If so, a parent task east of the first provides data.
!
                  KOUNT_PARENT_TASKS=KOUNT_PARENT_TASKS+1                   !<-- Increment parent task counter.
                  KP=KOUNT_PARENT_TASKS
                  SEND_TASK(KP)%ID=ID_1+1                                   !<-- Store the ID of this parent task to the north.
                  ID_E=SEND_TASK(KP)%ID
!
                  IF(ITS_PARENT_ON_CHILD(ID_E)+EPS>=I_UPDATE(3))THEN        !<-- 1st scenario of parent task east of the first.
                    SEND_TASK(KP)%ISTART(1)=NINT(ITS_PARENT_ON_CHILD(ID_E)) !<-- Starting I where parent task ID_E updates nest task.
                    SEND_TASK(KP)%IEND  (1)=I_END_X                         !<-- Ending I where parent task ID_E updates nest task.
                    SEND_TASK(KP)%JSTART(1)=J_START_X                       !<-- Starting J where parent task ID_E updates nest task.
                    SEND_TASK(KP)%JEND  (1)=MIN(J_END_X                  &  !<-- Ending J where parent task ID_E updates nest task.
                                           ,NINT(JTE_PARENT_ON_CHILD(ID_E))) 
                  ENDIF
!
                  IF(ITS_PARENT_ON_CHILD(ID_E)-EPS<=I_UPDATE(2))THEN 
                    SEND_TASK(KP)%ISTART(1)=I_UPDATE(3)                       !<-- Nest I where parent task ID_E begins updating 1st region.
                    SEND_TASK(KP)%IEND  (1)=I_END_X                           !<-- Nest I where parent task ID_E ends updating 1st region.
                    SEND_TASK(KP)%JSTART(1)=J_START_X                         !<-- Nest J where parent task ID_E begins updating 1st region.
                    SEND_TASK(KP)%JEND  (1)=MIN(J_UPDATE(2)                &  !<-- Nest J where parent task ID_E ends updating 1st region.
                                           ,NINT(JTE_PARENT_ON_CHILD(ID_E)))
!
                    IF(JTE_PARENT_ON_CHILD(ID_E)+EPS>=J_UPDATE(3))THEN
                      SEND_TASK(KP)%ISTART(2)=NINT(ITS_PARENT_ON_CHILD(ID_E)) !<-- Nest I where parent task ID_E begins updating 2nd region.
                      SEND_TASK(KP)%IEND  (2)=I_END_X                         !<-- Nest I where parent task ID_E ends updating 2nd region.
                      SEND_TASK(KP)%JSTART(2)=J_UPDATE(3)                     !<-- Nest J where parent task ID_E begins updating 2nd region.
                      SEND_TASK(KP)%JEND  (2)=MIN(J_END_X                  &  !<-- Nest I where parent task ID_E ends updating 2nd region.
                                             ,NINT(JTE_PARENT_ON_CHILD(ID_E))) 
                    ENDIF
                  ENDIF
!
                ENDIF ne_east
!
                EXIT search_i
!
              ENDIF ne
!
!-----------------------------------------------------------------------
!
            ENDIF update_j
!
          ENDDO search_j
!
        ENDIF update_i
!
      ENDDO search_i
!
!-----------------------------------------------------------------------
!***  Add up the number of points being updated by each parent task.
!-----------------------------------------------------------------------
!
      DO KP=1,KOUNT_PARENT_TASKS
!
        SEND_TASK(KP)%NPTS=(SEND_TASK(KP)%IEND(1)                       & 
                           -SEND_TASK(KP)%ISTART(1)+1)*                 & 
                           (SEND_TASK(KP)%JEND(1)                       &
                           -SEND_TASK(KP)%JSTART(1)+1)
!
        IF(SEND_TASK(KP)%ISTART(2)>=IMS)THEN                               !<-- Add points for 2nd regions on corners if present.
          SEND_TASK(KP)%NPTS=SEND_TASK(KP)%NPTS                         &
                           +(SEND_TASK(KP)%IEND(2)                      &
                            -SEND_TASK(KP)%ISTART(2)+1)*                &
                            (SEND_TASK(KP)%JEND(2)                      &
                            -SEND_TASK(KP)%JSTART(2)+1)
        ENDIF
!
      ENDDO
!
!-----------------------------------------------------------------------
!
      DEALLOCATE(ITS_PARENT_ON_CHILD)
      DEALLOCATE(ITE_PARENT_ON_CHILD)
      DEALLOCATE(JTS_PARENT_ON_CHILD)
      DEALLOCATE(JTE_PARENT_ON_CHILD)
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE MOVING_NEST_BOOKKEEPING
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE MOVING_NEST_RECV_DATA(COMM_TO_MY_PARENT                &
                                      ,NTIMESTEP                        &
                                      ,NUM_FIELDS_MOVE_2D_H_I           &
                                      ,NUM_FIELDS_MOVE_2D_X_I           &
                                      ,NUM_FIELDS_MOVE_2D_H_R           &
                                      ,NUM_FIELDS_MOVE_2D_X_R           &
                                      ,NUM_LEVELS_MOVE_3D_H             &
                                      ,NUM_FIELDS_MOVE_2D_V             &
                                      ,NUM_LEVELS_MOVE_3D_V             &
                                      ,SEND_TASK                        &
                                      ,EXPORT_STATE                     &
                                          )
!
!-----------------------------------------------------------------------
!***  After having determined which of their internal gridpoints
!***  need to be updated by which parent tasks following a nest's
!***  move, the nest's forecast tasks now receive the update data
!***  from the parent.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      INTEGER(kind=KINT),INTENT(IN) :: COMM_TO_MY_PARENT                &  !<-- MPI communicator from this nest to its parent
                                      ,NTIMESTEP                        &  !<-- Nest's current timestep
                                      ,NUM_FIELDS_MOVE_2D_H_I           &  !<-- # of 2-D internal state integer H variables to be updated
                                      ,NUM_FIELDS_MOVE_2D_X_I           &  !<-- # of 2-D integer H variables updated from external files
                                      ,NUM_FIELDS_MOVE_2D_H_R           &  !<-- # of 2-D internal state real H variables to be updated
                                      ,NUM_FIELDS_MOVE_2D_X_R           &  !<-- # of 2-D real H variables updated from external files
                                      ,NUM_LEVELS_MOVE_3D_H             &  !<-- # of 2-D levels in all 3-D H update variables
                                      ,NUM_FIELDS_MOVE_2D_V             &  !<-- # of 2-D internal state V variables to be updated
                                      ,NUM_LEVELS_MOVE_3D_V                !<-- # of 2-D levels in all 3-D V update variables
!
      TYPE(INTERIOR_DATA_FROM_PARENT),DIMENSION(1:4),INTENT(IN) ::      &
                                                              SEND_TASK    !<-- Specifics about interior data from sending parent tasks
!
      TYPE(ESMF_State),INTENT(INOUT) :: EXPORT_STATE                       !<-- The Parent-Child coupler export state
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: ITAG,N,NUM_PTASK_UPDATE,NUM_WORDS
!
      INTEGER(kind=KINT) :: IERR,RC,RC_RECV
!
      INTEGER(kind=KINT),DIMENSION(1:8) :: INDICES_H,INDICES_V
!
      INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE :: UPDATE_INTEGER_DATA
!
      INTEGER,DIMENSION(MPI_STATUS_SIZE) :: JSTAT
!
      REAL(kind=KFPT),DIMENSION(:),ALLOCATABLE :: UPDATE_REAL_DATA
!
      CHARACTER(len=1)  :: N_PTASK
      CHARACTER(len=12) :: NAME
      CHARACTER(len=17) :: NAME_REAL
      CHARACTER(len=20) :: NAME_INTEGER
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC     =ESMF_SUCCESS
      RC_RECV=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  First load into the Parent-Child coupler export state the
!***  number of parent tasks that send update data to this nest task.
!***  We insist that parent tasks will update the same H and V points
!***  with respect to their I,J indices.
!-----------------------------------------------------------------------
!
      NUM_PTASK_UPDATE=0
!
      DO N=1,4                                                             !<-- No more than 4 parent tasks will send data.
        IF(SEND_TASK(N)%ID<0)THEN
          EXIT
        ELSE
          NUM_PTASK_UPDATE=NUM_PTASK_UPDATE+1
        ENDIF
      ENDDO
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="MOVING_NEST_RECV_DATA: Load # of Parent Tasks Sending Interior Updates"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(state=EXPORT_STATE                         &  !<-- The Parent-Child coupler export state
                            ,name ='Num Parent Tasks Update'            &  !<-- Name of the variable
                            ,value=NUM_PTASK_UPDATE                     &  !<-- # of parent tasks that update this nest task
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RECV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  If no parent tasks are sending update data to this nest task
!***  then there is nothing more to do so RETURN.
!-----------------------------------------------------------------------
!
      IF(NUM_PTASK_UPDATE==0)RETURN                    
!
!-----------------------------------------------------------------------
!
      parent_tasks: DO N=1,NUM_PTASK_UPDATE
!
!-----------------------------------------------------------------------
!
        NUM_WORDS=(NUM_FIELDS_MOVE_2D_H_R-NUM_FIELDS_MOVE_2D_X_R        &  !<-- Total # of real words coming from Nth parent task
                  +NUM_LEVELS_MOVE_3D_H)*SEND_TASK(N)%NPTS              &
                 +(NUM_FIELDS_MOVE_2D_V+NUM_LEVELS_MOVE_3D_V)           &                                 
                  *SEND_TASK(N)%NPTS
!
        ALLOCATE(UPDATE_REAL_DATA(1:NUM_WORDS))                            !<-- Allocate the Recv buffer
!
        ITAG=NUM_WORDS+NTIMESTEP                                           !<-- Tag that changes for both data size and time
!
!-----------------------------------------------------------------------
!***  Receive the interior H and V real update data sent by
!***  parent task N.
!-----------------------------------------------------------------------
!
        CALL MPI_RECV(UPDATE_REAL_DATA                                  &  !<-- Real update data from Nth parent task
                     ,NUM_WORDS                                         &  !<-- # of real words received
                     ,MPI_REAL                                          &  !<-- The data is Real
                     ,SEND_TASK(N)%ID                                   &  !<-- Receive from parent task with this rank
                     ,ITAG                                              &  !<-- Unique MPI tag
                     ,COMM_TO_MY_PARENT                                 &  !<-- MPI communicator from this nest to its parent
                     ,JSTAT                                             &  !<-- MPI status object
                     ,IERR )
!
!-----------------------------------------------------------------------
!***  Load the update data and associated index limits into the
!***  Parent-Child coupler export state so it can be sent back into
!***  the DOMAIN component for incorporation.
!-----------------------------------------------------------------------
!
        WRITE(N_PTASK,'(I1)')N
        NAME_REAL='PTASK_REAL_DATA_'//N_PTASK
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Load # of Words in Real Update Data from Parent"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXPORT_STATE                       &  !<-- The Parent-Child coupler export state
                              ,name =NAME_REAL//' Words'                &  !<-- Name of the variable
                              ,value=NUM_WORDS                          &  !<-- Put # of real words here
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RECV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Load Real Update Data from Parent into P-C Cpl Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state    =EXPORT_STATE                   &  !<-- The Parent-Child coupler export state
                              ,name     =NAME_REAL                      &  !<-- Name of the variable
                              ,itemCount=NUM_WORDS                      &  !<-- # of words in real update data from parent task N
                              ,valueList=UPDATE_REAL_DATA               &  !<-- The real update data from parent task N
                              ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RECV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        DEALLOCATE(UPDATE_REAL_DATA)
!
!-----------------------------------------------------------------------
!***  There may or may not be integer variable updates at this time.
!-----------------------------------------------------------------------
!
        NUM_WORDS=(NUM_FIELDS_MOVE_2D_H_I-NUM_FIELDS_MOVE_2D_X_I)       &   !<-- Total # of integer words coming from
                  *SEND_TASK(N)%NPTS                                        !    the Nth parent task
!
!-----------------------------------------------------------------------
!***  Load into the Parent-Child coupler export state the number
!***  of integer words to be updated so the value can be sent to
!***  the DOMAIN component for incorporation of the integer data.
!-----------------------------------------------------------------------
!
        WRITE(N_PTASK,'(I1)')N
        NAME_INTEGER='PTASK_INTEGER_DATA_'//N_PTASK
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Load # of Words in Integer Update Data from Parent"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXPORT_STATE                       &  !<-- The Parent-Child coupler export state
                              ,name =NAME_INTEGER//' Words'             &  !<-- Name of the variable
                              ,value=NUM_WORDS                          &  !<-- Put # of integer words here
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RECV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
        recv_int: IF(NUM_WORDS>0)THEN
!
!-----------------------------------------------------------------------
!
          ALLOCATE(UPDATE_INTEGER_DATA(1:NUM_WORDS))                       !<-- Allocate the Recv buffer
!
          ITAG=NUM_WORDS+NTIMESTEP                                         !<-- Tag that changes for both data size and time
!
!-----------------------------------------------------------------------
!***  Receive the interior integer update data for H and V points
!***  sent by parent task N.
!-----------------------------------------------------------------------
!
          CALL MPI_RECV(UPDATE_INTEGER_DATA                             &  !<-- Integer update data from Nth parent task
                       ,NUM_WORDS                                       &  !<-- # of integer words received
                       ,MPI_INTEGER                                     &  !<-- The data is Integer
                       ,SEND_TASK(N)%ID                                 &  !<-- Receive from parent task with this rank
                       ,ITAG                                            &  !<-- Unique MPI tag
                       ,COMM_TO_MY_PARENT                               &  !<-- MPI communicator from this nest to its parent
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR )
!
!-----------------------------------------------------------------------
!***  Load the update data and associated index limits into the
!***  Parent-Child coupler export state so it can be sent back into
!***  the DOMAIN component for incorporation.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Load Integer Update Data from Parent into P-C Cpl Export State"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeSet(state    =EXPORT_STATE                 &  !<-- The Parent-Child coupler export state
                                ,name     =NAME_INTEGER                 &  !<-- Name of the variable
                                ,itemCount=NUM_WORDS                    &  !<-- # of words in integer update data from parent task N
                                ,valueList=UPDATE_INTEGER_DATA          &  !<-- The integer update data from parent task N
                                ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RECV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          DEALLOCATE(UPDATE_INTEGER_DATA)
!
!-----------------------------------------------------------------------
!
        ENDIF recv_int
!
!-----------------------------------------------------------------------
!
        INDICES_H(1)=SEND_TASK(N)%ISTART(1)
        INDICES_H(2)=SEND_TASK(N)%ISTART(2)
        INDICES_H(3)=SEND_TASK(N)%IEND(1)
        INDICES_H(4)=SEND_TASK(N)%IEND(2)
        INDICES_H(5)=SEND_TASK(N)%JSTART(1)
        INDICES_H(6)=SEND_TASK(N)%JSTART(2)
        INDICES_H(7)=SEND_TASK(N)%JEND(1)
        INDICES_H(8)=SEND_TASK(N)%JEND(2)
!
        INDICES_V(1)=SEND_TASK(N)%ISTART(1)
        INDICES_V(2)=SEND_TASK(N)%ISTART(2)
        INDICES_V(3)=SEND_TASK(N)%IEND(1)
        INDICES_V(4)=SEND_TASK(N)%IEND(2)
        INDICES_V(5)=SEND_TASK(N)%JSTART(1)
        INDICES_V(6)=SEND_TASK(N)%JSTART(2)
        INDICES_V(7)=SEND_TASK(N)%JEND(1)
        INDICES_V(8)=SEND_TASK(N)%JEND(2)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Load Index Limits for Update Data from Parent"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        WRITE(N_PTASK,'(I1)')N
        NAME='PTASK_DATA_'//N_PTASK
!
        CALL ESMF_AttributeSet(state    =EXPORT_STATE                   &  !<-- The Parent-Child coupler export state
                              ,name     =NAME//' Indices H'             &  !<-- Name of the variable
                              ,itemCount=N8                             &  !<-- # of words in index limits of update data
                              ,valueList=INDICES_H                      &  !<-- The update data index specifications for H
                              ,rc       =RC)
!
        CALL ESMF_AttributeSet(state    =EXPORT_STATE                   &  !<-- The Parent-Child coupler export state
                              ,name     =NAME//' Indices V'             &  !<-- Name of the variable
                              ,itemCount=N8                             &  !<-- # of words in index limits of update data
                              ,valueList=INDICES_V                      &  !<-- The update data index specifications for V
                              ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RECV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
      ENDDO parent_tasks
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE MOVING_NEST_RECV_DATA
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE PARENT_UPDATES_HALOS(FLAG_H_OR_V                       &
                                     ,MOVE_BUNDLE                       &
                                     ,NFLDS_3DR                         &
                                     ,NFLDS_2DR                         &
                                     ,NFLDS_2DI                         &
                                       )
!
!-----------------------------------------------------------------------
!***  Before a parent can update locations on its moving nests' domains
!***  it must perform halo exchanges for all those variables specified
!***  for use in updates but which do not have their halos exchanged
!***  during the normal integration.  Use of the parent tasks halo
!***  regions cannot be avoided during the nest point updates.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      CHARACTER(len=1),INTENT(IN) :: FLAG_H_OR_V
!
      TYPE(ESMF_FieldBundle),INTENT(INOUT) :: MOVE_BUNDLE                  !<-- ESMF Bundle of 2-D and 3-D arrays specified for updating
!
      INTEGER(kind=KINT),INTENT(IN) :: NFLDS_2DR                        &  !<-- # of 2-D real arrays specified for updating
                                      ,NFLDS_3DR                           !<-- # of 3-D real arrays specified for updating
!
      INTEGER(kind=KINT),INTENT(IN),OPTIONAL :: NFLDS_2DI                  !<-- # of 2-D integer arrays specified for updating
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: N_FIELD,N_REMOVE,NUM_DIMS                   &
                           ,NUM_FIELDS_MOVE,NUM_LEVELS                  &
                           ,RC,RC_FINAL
!
      INTEGER(kind=KINT),DIMENSION(1:3) :: LIMITS_HI                    &
                                          ,LIMITS_LO
!
      REAL(kind=KFPT),DIMENSION(:,:),POINTER :: ARRAY_2D
!
      REAL(kind=KFPT),DIMENSION(:,:,:),POINTER :: ARRAY_3D
!
      CHARACTER(len=30) :: FIELD_NAME
!
      TYPE(ESMF_Field) :: HOLD_FIELD
!
      TYPE(ESMF_TypeKind_Flag) :: DATATYPE
!
      LOGICAL(kind=KLOG) :: EXCH_NEEDED
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  What is the total number of Fields in the update data BUNDLE?
!-----------------------------------------------------------------------
!
      IF(FLAG_H_OR_V=='H')THEN
        NUM_FIELDS_MOVE=NFLDS_2DI+NFLDS_2DR+NFLDS_3DR
!
      ELSEIF(FLAG_H_OR_V=='V')THEN
        NUM_FIELDS_MOVE=NFLDS_2DR+NFLDS_3DR
      ENDIF
!
!-----------------------------------------------------------------------
!***  Check each Field to see if its array has its halo exchanged
!***  during the integration.
!-----------------------------------------------------------------------
!
      field_loop: DO N_FIELD=1,NUM_FIELDS_MOVE
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract Each Field From Move_Bundle"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_FieldBundleGet(FIELDBUNDLE=MOVE_BUNDLE                &  !<-- Bundle holding the arrays for move updates
                                ,fieldIndex =N_FIELD                    &  !<-- Index of the Field in the Bundle
                                ,field      =HOLD_FIELD                 &  !<-- Field N_FIELD in the Bundle
                                ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Type, Dimensions, Name of the Field"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_FieldGet(field   =HOLD_FIELD                          &  !<-- Field N_FIELD in the Bundle
                          ,dimCount=NUM_DIMS                            &  !<-- Is this Field 2-D or 3-D?
                          ,typeKind=DATATYPE                            &  !<-- Is this Field integer or real?
                          ,name    =FIELD_NAME                          &  !<-- This Field's name
                          ,rc      =RC )
!
        N_REMOVE=INDEX(FIELD_NAME,SUFFIX_MOVE)
        FIELD_NAME=FIELD_NAME(1:N_REMOVE-1)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  None of the variables needed by the parent for updating its
!***  moving nests are type integer so we can skip those outright.
!-----------------------------------------------------------------------
!
        IF(DATATYPE==ESMF_TYPEKIND_I4)THEN
          CYCLE field_loop
        ENDIF
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract the Halo Exchange Flag"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(field=HOLD_FIELD                         &  !<-- Take Attribute from this Field
                              ,name ='EXCH_NEEDED'                      &  !<-- The Attribute's name
                              ,value=EXCH_NEEDED                        &  !<-- The Attribute's value
                              ,rc   =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Move to the next Field if a halo exchange is not needed.
!-----------------------------------------------------------------------
!
        IF(.NOT.EXCH_NEEDED)THEN
          CYCLE field_loop
        ENDIF
!
!-----------------------------------------------------------------------
!***  2-D Fields 
!-----------------------------------------------------------------------
!
        dims_2_or_3: IF(NUM_DIMS==2)THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Extract the 2-D Array from Field"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_FieldGet(field    =HOLD_FIELD                       &  !<-- Field N_FIELD in the Bundle
                            ,localDe  =0                                &
                            ,farrayPtr=ARRAY_2D                         &  !<-- Dummy 2-D real array with Field's data
                            ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL HALO_EXCH(ARRAY_2D,1,1,1)
!
!-----------------------------------------------------------------------
!***  3-D Fields 
!-----------------------------------------------------------------------
!
        ELSEIF(NUM_DIMS==3)THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Extract the 3-D Array from Field"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_FieldGet(field      =HOLD_FIELD                     &  !<-- Field N_FIELD in the Bundle
                            ,localDe    =0                              &
                            ,farrayPtr  =ARRAY_3D                       &  !<-- Dummy 3-D real array with Field's data
                            ,totalLBound=LIMITS_LO                      &  !<-- Starting index in each dimension
                            ,totalUBound=LIMITS_HI                      &  !<-- Ending index in each dimension
                            ,rc     =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          NUM_LEVELS=LIMITS_HI(3)-LIMITS_LO(3)+1
!
          CALL HALO_EXCH(ARRAY_3D,NUM_LEVELS,1,1)
!
!-----------------------------------------------------------------------
!
        ENDIF dims_2_or_3
!
!-----------------------------------------------------------------------
!
      ENDDO field_loop
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE PARENT_UPDATES_HALOS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE PARENT_BOOKKEEPING_MOVING(I_PARENT_SW_NEW              &
                                          ,J_PARENT_SW_NEW              &
                                          ,I_PARENT_SW_OLD              &
                                          ,J_PARENT_SW_OLD              &
                                          ,ITS,ITE,JTS,JTE              &
                                          ,NUM_CHILD_TASKS              &
                                          ,CHILD_TASK_LIMITS            &
                                          ,PARENT_CHILD_SPACE_RATIO     &
                                          ,NHALO                        &
                                          ,NROWS_P_UPD_W                &
                                          ,NROWS_P_UPD_E                &
                                          ,NROWS_P_UPD_S                &
                                          ,NROWS_P_UPD_N                &
                                          ,N_UPDATE_CHILD_TASKS         &
                                          ,TASK_UPDATE_SPECS            &
                                          ,HANDLE_UPDATE                &
                                          ,CHILD_UPDATE_DATA            &
                                            )
!
!-----------------------------------------------------------------------
!***  This parent has learned that one of its children wants to move
!***  to a new location therefore the parent must determine which 
!***  points on which child tasks need to be updated by which of its
!***  own tasks.  Update points on nests are those that lie outside
!***  of the nest's pre-move footprint following the move.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      INTEGER(kind=KINT),INTENT(IN) :: I_PARENT_SW_NEW                  &  !<-- SW corner of nest on this parent I after move
                                      ,I_PARENT_SW_OLD                  &  !<-- SW corner of nest on this parent I before move
                                      ,J_PARENT_SW_NEW                  &  !<-- SW corner of nest on this parent J after move
                                      ,J_PARENT_SW_OLD                  &  !<-- SW corner of nest on this parent J before move
!
                                      ,ITS,ITE,JTS,JTE                  &  !<-- Subdomain integration limits of parent task
!
                                      ,NHALO                            &  !<-- # of halo points
                                      ,NROWS_P_UPD_W                    &  !<-- Moving nest footprint W bndry rows updated by parent
                                      ,NROWS_P_UPD_E                    &  !<-- Moving nest footprint E bndry rows updated by parent
                                      ,NROWS_P_UPD_S                    &  !<-- Moving nest footprint S bndry rows updated by parent
                                      ,NROWS_P_UPD_N                    &  !<-- Moving nest footprint N bndry rows updated by parent
                                      ,NUM_CHILD_TASKS                  &  !<-- # of child forecast tasks
                                      ,PARENT_CHILD_SPACE_RATIO            !<-- # of child grid increments in one of parent's
!
      INTEGER(kind=KINT),DIMENSION(1:4,NUM_CHILD_TASKS),INTENT(IN) ::   &
                                                     CHILD_TASK_LIMITS     !<-- ITS,ITE,JTS,JTE for each child forecast task
!
      INTEGER(kind=KINT),INTENT(INOUT) :: N_UPDATE_CHILD_TASKS             !<-- # of moving nest tasks to be updated by this parent task
!
      INTEGER(kind=KINT),DIMENSION(1:NUM_CHILD_TASKS),INTENT(IN) ::     &
                                                      HANDLE_UPDATE        !<-- MPI Handles for ISends to the child tasks
!
      TYPE(CHILD_UPDATE_LINK),TARGET,INTENT(INOUT) :: TASK_UPDATE_SPECS    !<-- Linked list with nest task update region specifications
!
      TYPE(MIXED_DATA_TASKS),INTENT(INOUT) :: CHILD_UPDATE_DATA            !<-- Composite of all update data from parent for nest tasks
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: I_SHIFT,I1,I2                               &
                           ,IDE_CHILD,IDS_CHILD                         &
                           ,IMS_CHILD,IME_CHILD                         &
                           ,IDE_FOOTPRINT,IDS_FOOTPRINT                 &
                           ,ITE_PARENT_ON_CHILD,ITS_PARENT_ON_CHILD     &
                           ,J_SHIFT,J1,J2                               &
                           ,JDE_CHILD,JDS_CHILD                         &
                           ,JMS_CHILD,JME_CHILD                         &
                           ,JDE_FOOTPRINT,JDS_FOOTPRINT                 &
                           ,JTE_PARENT_ON_CHILD,JTS_PARENT_ON_CHILD     &
                           ,KOUNT_TASKS,N,NN
!
      INTEGER(kind=KINT) :: IERR,ISTAT
!
      INTEGER,DIMENSION(MPI_STATUS_SIZE) :: JSTAT
!
      REAL(kind=KFPT) :: I1R,I2R                                        &
                        ,ITE_PARENT_ON_CHILD_R,ITS_PARENT_ON_CHILD_R    &
                        ,J1R,J2R                                        &
                        ,JTE_PARENT_ON_CHILD_R,JTS_PARENT_ON_CHILD_R
!
      TYPE(CHILD_UPDATE_LINK),POINTER :: PTR,PTR_X
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Prior to doing anything related to updating nest tasks 
!***  following the latest move, be sure that all update data
!***  that this parent task might have sent to any nest tasks
!***  following the preceding move has indeed been received
!***  by all of those tasks whether or not this parent task
!***  will send to any of the same nest tasks this time.
!-----------------------------------------------------------------------
!
      PTR=>TASK_UPDATE_SPECS                                               !<-- Start at the top of the list of updated nest tasks
!
      DO WHILE(ASSOCIATED(PTR%TASK_ID))                                    !<-- A link exists if TASK_ID is associated.
        CALL MPI_WAIT(HANDLE_UPDATE(PTR%TASK_ID)                        &  !<-- Handle for ISend from parent task to child task
                     ,JSTAT                                             &  !<-- MPI status
                     ,IERR )
        IF(ASSOCIATED(PTR%NEXT_LINK))THEN
          PTR=>PTR%NEXT_LINK
        ELSE
          EXIT
        ENDIF
      ENDDO
!
!-----------------------------------------------------------------------
!***  All nest tasks have received data from the previous move
!***  so proceed with deleting those data objects.
!-----------------------------------------------------------------------
!
      PTR_X=>TASK_UPDATE_SPECS                                            !<-- Go back to the top of the list of updated nest tasks
      KOUNT_TASKS=0
!
      DO WHILE(ASSOCIATED(PTR_X%TASK_ID))                                 !<-- An old link exists if TASK_ID is associated.
        KOUNT_TASKS=KOUNT_TASKS+1
        DEALLOCATE(PTR_X%TASK_ID,stat=ISTAT)
        IF(ISTAT/=0)THEN
          WRITE(0,*)' Failed to deallocate TASK_UPDATE_SPECS%TASK_ID for nest task #',KOUNT_TASKS,' stat=',istat
        ENDIF
        DEALLOCATE(PTR_X%NUM_PTS_UPDATE_HZ,stat=ISTAT)
        DEALLOCATE(PTR_X%IL,stat=ISTAT)
        DEALLOCATE(PTR_X%JL,stat=ISTAT)
!
        TAIL=>NULL()
        IF(ASSOCIATED(PTR_X%NEXT_LINK))THEN                               !<-- If another links exists, point to it.
          TAIL=>PTR_X%NEXT_LINK
        ENDIF
!
        IF(KOUNT_TASKS>1)THEN                                      !<-- The top of TASK_UPDATE_SPECS is allocatable array element N 
          DEALLOCATE(PTR_X,stat=ISTAT)                             !    (for the Nth moving child) and is not a pointer.
          IF(ISTAT/=0)THEN
            WRITE(0,*)' Failed to deallocate TASK_UPDATE_SPECS for nest task #',KOUNT_TASKS,' stat=',istat
          ENDIF
        ENDIF
!
!---------------------------------------------------------------
!***  Precisely the same nest tasks are updated for both
!***  H and V points therefore the deallocation of working
!***  pointers for nest tasks in the following block is 
!***  removing all data for both types of points and not
!***  leaving some behind of one type or the other.
!---------------------------------------------------------------
!
        IF(ASSOCIATED(CHILD_UPDATE_DATA%TASKS(KOUNT_TASKS)%DATA_INTEGER))THEN
          DEALLOCATE(CHILD_UPDATE_DATA%TASKS(KOUNT_TASKS)%DATA_INTEGER,stat=ISTAT)
          IF(ISTAT/=0)then
            WRITE(0,*)' Failed to deallocate CHILD_UPDATE_DATA%TASKS(KOUNT_TASKS)%DATA_INTEGER' &
                     ,' for KOUNT_TASKS=',kount_tasks,' stat=',istat
          ENDIF
        ENDIF
!
        IF(ASSOCIATED(CHILD_UPDATE_DATA%TASKS(KOUNT_TASKS)%DATA_REAL))THEN
          DEALLOCATE(CHILD_UPDATE_DATA%TASKS(KOUNT_TASKS)%DATA_REAL,stat=ISTAT)
          IF(ISTAT/=0)then
            WRITE(0,*)' Failed to deallocate CHILD_UPDATE_DATA%TASKS(KOUNT_TASKS)%DATA_REAL' &
                     ,' for KOUNT_TASKS=',kount_tasks,' stat=',istat
          ENDIF
        ENDIF
!
        IF(ASSOCIATED(TAIL))THEN
          PTR_X=>TAIL                                                      !<-- There is still another old link
        ELSE
          EXIT                                                             !<-- The last link in this list has been deallocated
        ENDIF
!
      ENDDO
!
      IF(ASSOCIATED(CHILD_UPDATE_DATA%TASKS))THEN
        DEALLOCATE(CHILD_UPDATE_DATA%TASKS)
      ENDIF
!
!-----------------------------------------------------------------------
!***  How far did the nest move on the parent grid?
!-----------------------------------------------------------------------
!
      I_SHIFT=I_PARENT_SW_NEW-I_PARENT_SW_OLD
      J_SHIFT=J_PARENT_SW_NEW-J_PARENT_SW_OLD
!
!-----------------------------------------------------------------------
!***  What are this parent task's integration limits
!***  in terms of the moving nest's grid indices?
!***  To figure that out begin with the values of the
!***  index limits of the entire moving nest domain.
!-----------------------------------------------------------------------
!
      IDS_CHILD=CHILD_TASK_LIMITS(1,1)                                     !<-- Index limits of the moving nest on
      IDE_CHILD=CHILD_TASK_LIMITS(2,NUM_CHILD_TASKS)                       !    its own grid.
      JDS_CHILD=CHILD_TASK_LIMITS(3,1)                                     !
      JDE_CHILD=CHILD_TASK_LIMITS(4,NUM_CHILD_TASKS)                       !<--
!
!-----------------------------------------------------------------------
!***  In the following diagram 'H' represents mass points on the
!***  parent grid while 'h' represents mass points on the nest grid.
!***  Gridpoint values on the top are with respect to the nest.
!***  Gridpoint values on the bottom are with respect to the parent.
!***  The Parent-Child space ratio is 3.  The given parent task must
!***  cover the gap between its ITE and the next parent task's ITS.
!***  'Hh' indicates that parent and nest points coincide.
!-----------------------------------------------------------------------
!
!
!  ITS_PARENT_ON_CHILD=-5     I=1                     ITE_PARENT_ON_CHILD=9
!    |                         |                                 |    
!    |                         |                                 |          
!   Hh   h   h   Hh   h   h   Hh   h   h   Hh   h   h   Hh   h   h   Hh
!   |                         |                         |
!   |                         |                         |<--gap-->
!  ITS_PARENT=1           I_PARENT_SW=3           ITE_PARENT=5
!
!
!-----------------------------------------------------------------------
!
      ITS_PARENT_ON_CHILD=IDS_CHILD-(I_PARENT_SW_NEW-ITS)               &  !<-- ITS of parent task in child's coordinate space
                                    *PARENT_CHILD_SPACE_RATIO              !    for H points
!
      ITE_PARENT_ON_CHILD=IDS_CHILD-(I_PARENT_SW_NEW-ITE)               &  !<-- ITE of parent task in child's coordinate space
                                    *PARENT_CHILD_SPACE_RATIO           &  !    for H points
                                    +PARENT_CHILD_SPACE_RATIO-1            !<-- Filling in gap beyond last parent point on nest grid
!
      JTS_PARENT_ON_CHILD=JDS_CHILD-(J_PARENT_SW_NEW-JTS)               &  !<-- JTS of parent task in child's coordinate space
                                    *PARENT_CHILD_SPACE_RATIO              !    for H points
!
      JTE_PARENT_ON_CHILD=JDS_CHILD-(J_PARENT_SW_NEW-JTE)               &  !<-- JTE of parent task in child's coordinate space
                                    *PARENT_CHILD_SPACE_RATIO           &  !    for H points
                                    +PARENT_CHILD_SPACE_RATIO-1            !<-- Filling in gap beyond last parent point on nest grid
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  The situation for V points is necessarily more complex.
!***  In the following diagram 'H' represents mass points on the
!***  parent grid and 'h' represents mass points on the nest grid
!***  while 'V' and 'v' represent the velocity points on the respective
!***  grids.  Gridpoint values on the top are with respect to the
!***  nest's v points.  The Parent-Child space ratio is 3.
!***  'Hh' and 'Vv' indicate that parent and nest points coincide.
!***  Note the correspondence of the V diagram below with the H diagram
!***  above.  The nest's v points for which a parent task is responsible
!***  have exactly the same indices as the nest h points for which that
!***  parent task is responsible.  Although doing this means that 
!***  ITS_PARENT_ON_CHILD is not at the same location as ITS_PARENT,
!***  it is required for exactly the same nest tasks to be updated by
!***  a parent task for both the h and v points.  Likewise for
!***  ITE_PARENT_ON_CHILD, etc.
!***  The reason the  relationships on velocity points are much more
!***  complicated than on mass points is that the SW corner point
!***  which serves as the anchor of the nest is always an H/h point.
!-----------------------------------------------------------------------
!
!
!  ITS_PARENT_ON_CHILD_R=-5.    I=1                     ITE_PARENT_ON_CHILD_R=9.
!      |                         |                                 |        
!      |                         |                                 |        
!   Hh | h   h   Hh   h   h   Hh | h   h   Hh   h   h   Hh   h   h | Hh     
!   |  |                         v   v   v               |<-gap->| |        
!   |  |                                                           |        
!   |  v   Vv  v    v   Vv  v    v   Vv  v   v   Vv   v    v   Vv  v   v   Vv
!   |      |                         |                         |
!   |      |                     v   |   v                     |            
!   Hh   h | h   Hh   h   h   Hh   h | h   Hh   h   h   Hh   h | h   Hh
!   |      |                  |      |                  |      |
!   |      |                  |      |                  |      |
!   |  ITS_PARENT=1           | I_PARENT_SW=3           | ITE_PARENT=5
!   |     on V                |     on V                |    on V
!   |                         |                         |
!   |                         |                         |
!  ITS_PARENT=1           I_PARENT_SW=3           ITE_PARENT=5
!     on H                    on H                   on H
!
!
!-----------------------------------------------------------------------
!***  However the logic has been constructed such that the index limits
!***  on each moving nest task subdomain for which each parent task
!***  must provide update data are identical for H and V points so we
!***  need only use the simpler perspective of H points to find those
!***  limits.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Boundary of the nest's pre-move footprint in terms of the
!***  nest's new position.
!-----------------------------------------------------------------------
!
      IDS_FOOTPRINT=IDS_CHILD-I_SHIFT*PARENT_CHILD_SPACE_RATIO
      IDE_FOOTPRINT=IDE_CHILD-I_SHIFT*PARENT_CHILD_SPACE_RATIO
      JDS_FOOTPRINT=JDS_CHILD-J_SHIFT*PARENT_CHILD_SPACE_RATIO
      JDE_FOOTPRINT=JDE_CHILD-J_SHIFT*PARENT_CHILD_SPACE_RATIO
!
!-----------------------------------------------------------------------
!***  Loop through the nest's task subdomains.
!-----------------------------------------------------------------------
!
      child_tasks: DO N=1,NUM_CHILD_TASKS
!
!-----------------------------------------------------------------------
!***  What are child task N's memory limits?  We use those limits
!***  since the parent task updates both integration and halo points
!***  on the child's subdomains in order to avoid all of the 
!***  communication involved in doing halo exchanges following the
!***  updates.  The parent uses only its integration points (no halo
!***  points) to do the updating.
!-----------------------------------------------------------------------
!
        IMS_CHILD=MAX(CHILD_TASK_LIMITS(1,N)-NHALO,IDS_CHILD)
        IME_CHILD=MIN(CHILD_TASK_LIMITS(2,N)+NHALO,IDE_CHILD)
        JMS_CHILD=MAX(CHILD_TASK_LIMITS(3,N)-NHALO,JDS_CHILD)
        JME_CHILD=MIN(CHILD_TASK_LIMITS(4,N)+NHALO,JDE_CHILD)
!
!-----------------------------------------------------------------------
!***  Do any of child task N's H points lie within this parent task's
!***  subdomain for the new nest position?
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
        limits: IF((IMS_CHILD>=ITS_PARENT_ON_CHILD                      &
                                .AND.                                   &
                    IMS_CHILD<=ITE_PARENT_ON_CHILD                      &
                                .OR.                                    &
                    IME_CHILD>=ITS_PARENT_ON_CHILD                      &
                                .AND.                                   &
                    IME_CHILD<=ITE_PARENT_ON_CHILD)                     &
                                .AND.                                   &
                   (JMS_CHILD>=JTS_PARENT_ON_CHILD                      &
                                .AND.                                   &
                    JMS_CHILD<=JTE_PARENT_ON_CHILD                      &
                                .OR.                                    &
                    JME_CHILD>=JTS_PARENT_ON_CHILD                      &
                                .AND.                                   &
                    JME_CHILD<=JTE_PARENT_ON_CHILD))THEN                   !<-- If so, some of child task N's points are within
!                                                                          !    this parent task's region of responsibility for
                                                                           !    updating post-move nest points.        
!-----------------------------------------------------------------------
!***  The intersection of child task N's subdomain with this parent
!***  task's region.
!-----------------------------------------------------------------------
!
          I1=MAX(IMS_CHILD,ITS_PARENT_ON_CHILD)                            !<-- I limits of child task N's subdomain that lies
          I2=MIN(IME_CHILD,ITE_PARENT_ON_CHILD)                            !    within this parent task's subdomain.
!
          J1=MAX(JMS_CHILD,JTS_PARENT_ON_CHILD)                            !<-- J limits of child task N's subdomain that lies
          J2=MIN(JME_CHILD,JTE_PARENT_ON_CHILD)                            !    within this parent task's subdomain.
!
!-----------------------------------------------------------------------
!***  The parent task will update only those nest H points that lie
!***  outside of the footprint of the nest domain's pre-move position.
!***  If all the nest points in child task N's subdomain lie within
!***  the footprint then the parent task has nothing to do so move on
!***  to the next child task.
!
!***  NOTE:  The north and east limits of the nest domain's pre-move
!***         footprint cannot be used as a source for post-move updates
!***         in the intra-task and inter-task shifts of data.  That is
!***         because the V-pt variables there are not part of the nest
!***         integration therefore their values are not valid.  So we 
!***         also must not use the H-pt variables at those same limits
!***         or else occasions would arise when nest tasks receiving
!***         H-pt updates would not be exactly the same as the nest
!***         tasks receiving V-pt updates.  That situation is avoided
!***         or else the bookkeeping would be even more complicated.
!***         The parent will update H-pt and V-pt variables along the
!***         nest domain's pre-move footprint's north and east limits.
!***         But the intra- and inter-task shifts also cannot do the
!***         updating of the nest domain's southern and western boundary
!***         because many of the nest variables do not have valid
!***         integration values there so the parent must also update
!***         those nest boundaries following a shift.  Moreover the
!***         dynamical tendencies for T, U, and V are not computed in
!***         the next to the outermost row of the domain which means the
!***         parent will have to update all nest points that move to
!***         IDE and IDE-1 and JDE and JDE-1 on the pre-move footprint.
!***         Use variables for the depth to which the parent will 
!***         provide update data to nest points within the footprint
!***         in case that depth needs to change in the future.
!-----------------------------------------------------------------------
!
          IF(I1>=IDS_FOOTPRINT+NROWS_P_UPD_W                            &
                         .AND.                                          &
             I2<=IDE_FOOTPRINT-NROWS_P_UPD_E                            &
                         .AND.                                          &
             J1>=JDS_FOOTPRINT+NROWS_P_UPD_S                            &
                         .AND.                                          &
             J2<=JDE_FOOTPRINT-NROWS_P_UPD_N )THEN                         !<-- If so, these nest points lie entirely within the footprint.
!
              CYCLE child_tasks                                            !<-- So this child task receives no updating from this
!                                                                          !    parent task.
          ENDIF
!
!-----------------------------------------------------------------------
!***  Now we know this parent task is updating at least some H points
!***  within child task N's subdomain so allocate a link in the
!***  linked list that holds update information about task N.
!***  We use a linked list because we do not know a priori how many
!***  child tasks need updates from each parent task and that number
!***  will change with each shift of the nest.
!-----------------------------------------------------------------------
!
          N_UPDATE_CHILD_TASKS=N_UPDATE_CHILD_TASKS+1
!
          CALL PARENT_FINDS_UPDATE_LIMITS
!
!-----------------------------------------------------------------------
!
        ENDIF limits
!
!-----------------------------------------------------------------------
!***  The parent task now knows which H points on child task N's
!***  subdomain that it must update following the nest's move.
!***  Those same index limits will apply to V point updates even
!***  though the physical locations differ.
!-----------------------------------------------------------------------
!
      ENDDO child_tasks
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!
      SUBROUTINE PARENT_FINDS_UPDATE_LIMITS
!
!-----------------------------------------------------------------------
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: ISTAT,NLOC
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Add a link to this parent task's linked list of moving nest
!***  specifications.  Each new link is associated with another
!***  nest task that needs updating by this parent task on the
!***  current moving nest.
!-----------------------------------------------------------------------
!
      IF(N_UPDATE_CHILD_TASKS==1)THEN
        TAIL=>TASK_UPDATE_SPECS                                            !<-- For the 1st link, point at the top of the list.
        NULLIFY(TAIL%NEXT_LINK)
      ELSE
        ALLOCATE(TAIL%NEXT_LINK,stat=ISTAT)                                !<-- Add a new link for each additional child task
        TAIL=>TAIL%NEXT_LINK                                               !<-- Point at the new link so it is ready to use
        NULLIFY(TAIL%NEXT_LINK)
      ENDIF
!
      ALLOCATE(TAIL%TASK_ID)                                               !<-- Allocate the pieces of data in this link
      ALLOCATE(TAIL%NUM_PTS_UPDATE_HZ)                                     !
      ALLOCATE(TAIL%IL(1:4))                                               !
      ALLOCATE(TAIL%JL(1:4))                                               !<--
!
      TAIL%TASK_ID=N                                                       !<-- Task is Nth among all tasks on this child
!
      DO NLOC=1,4
        TAIL%IL(NLOC)=-999
        TAIL%JL(NLOC)=-999
      ENDDO
!
!-----------------------------------------------------------------------
!***  The simplest case occurs when all of child task N's subdomain
!***  that intersects this parent task's subdomain lies outside of
!***  the pre-move footprint.  The parent task then just updates 
!***  all those nest points.
!-----------------------------------------------------------------------
!
      parent_updates: IF(I2<=IDS_FOOTPRINT+NROWS_P_UPD_W-1              &
                                     .OR.                               &
                         I1>=IDE_FOOTPRINT-NROWS_P_UPD_E+1              &
                                     .OR.                               &
                         J2<=JDS_FOOTPRINT+NROWS_P_UPD_S-1              &
                               .OR.                                     &
                         J1>=JDE_FOOTPRINT-NROWS_P_UPD_N+1 )THEN
!
!-----------------------------------------------------------------------
!
        TAIL%IL(1)=I1                                                      !<-- I limits of nest task N's update region by parent task
        TAIL%IL(2)=I2                                                      !    in terms of the nest's grid.
!
        TAIL%JL(1)=J1                                                      !<-- J limits of nest task N's update region by parent task
        TAIL%JL(2)=J2                                                      !    in terms of the nest's grid.
!
!-----------------------------------------------------------------------
!***  What remains are intersections between child task N's subdomain
!***  and this parent task's subdomain that lie along the edge of the
!***  pre-move footprint.  Usually these regions will be a rectangle.
!***  However if both child task N and this parent task cover a corner 
!***  of the footprint then the update region of the child task's
!***  subdomain is not a simple rectangle; essentially it is two 
!***  rectangles.
!***  See diagrams in subroutine RECV_INTERIOR_DATA_FROM_PARENT
!***  in this module.
!-----------------------------------------------------------------------
!
      ELSE parent_updates
!
        IF(I1>=IDS_FOOTPRINT+NROWS_P_UPD_W                              &
                  .AND.                                                 &  
           I2<=IDE_FOOTPRINT-NROWS_P_UPD_E )THEN                           !<-- Rectangular update region on S/N edge of footprint.
!
          TAIL%IL(1)=I1
          TAIL%IL(2)=I2
!
          IF(J1<=JDS_FOOTPRINT+NROWS_P_UPD_S-1)THEN                        !<-- Rectangular update region on south edge of footprint.
            TAIL%JL(1)=J1
            TAIL%JL(2)=JDS_FOOTPRINT+NROWS_P_UPD_S-1
!
          ELSEIF(J2>=JDE_FOOTPRINT-NROWS_P_UPD_N+1)THEN                    !<-- Rectangular update region on north edge of footprint.
            TAIL%JL(1)=JDE_FOOTPRINT-NROWS_P_UPD_N+1
            TAIL%JL(2)=J2
          ENDIF
!
        ELSEIF(J1>=JDS_FOOTPRINT+NROWS_P_UPD_S                          &  
                     .AND.                                              & 
               J2<=JDE_FOOTPRINT-NROWS_P_UPD_N )THEN                       !<-- Rectangular update region on W/E edge of footprint.
!
          TAIL%JL(1)=J1
          TAIL%JL(2)=J2
!
          IF(I1<=IDS_FOOTPRINT+NROWS_P_UPD_W-1)THEN                        !<-- Rectangular update region on west edge of footprint.
            TAIL%IL(1)=I1
            TAIL%IL(2)=IDS_FOOTPRINT+NROWS_P_UPD_W-1
!
          ELSEIF(I2>=IDE_FOOTPRINT-NROWS_P_UPD_N+1)THEN                    !<-- Rectangular update region on east edge of footprint.
            TAIL%IL(1)=IDE_FOOTPRINT-NROWS_P_UPD_N+1
            TAIL%IL(2)=I2
          ENDIF
!
        ELSEIF(I1<=IDS_FOOTPRINT+NROWS_P_UPD_W-1                        &
                            .AND.                                       &
               I2>=IDS_FOOTPRINT+NROWS_P_UPD_W )THEN                       !<-- Child task update region on SW/NW corner of footprint.
!
          IF(J1<=JDS_FOOTPRINT+NROWS_P_UPD_S-1)THEN                        !<-- Child task update region on SW corner of footprint.
            TAIL%IL(1)=I1
            TAIL%IL(2)=I2
            TAIL%IL(3)=I1
            TAIL%IL(4)=IDS_FOOTPRINT+NROWS_P_UPD_W-1
            TAIL%JL(1)=J1
            TAIL%JL(2)=JDS_FOOTPRINT+NROWS_P_UPD_S-1
            TAIL%JL(3)=TAIL%JL(2)+1
            TAIL%JL(4)=J2
!
          ELSEIF(J2>=JDE_FOOTPRINT-NROWS_P_UPD_N-1)THEN                    !<-- Child task update region on NW corner of footprint.
            TAIL%IL(1)=I1
            TAIL%IL(2)=IDS_FOOTPRINT+NROWS_P_UPD_W-1
            TAIL%IL(3)=I1
            TAIL%IL(4)=I2
            TAIL%JL(1)=J1
            TAIL%JL(2)=JDE_FOOTPRINT-NROWS_P_UPD_N
            TAIL%JL(3)=TAIL%JL(2)+1
            TAIL%JL(4)=J2
          ENDIF
!
        ELSEIF(I1<=IDE_FOOTPRINT-NROWS_P_UPD_E                          &
                       .AND.                                            &
               I2>=IDE_FOOTPRINT-NROWS_P_UPD_E+1 )THEN                     !<-- Child task update region on SE/NE corner of footprint
!
          IF(J1<=JDS_FOOTPRINT+NROWS_P_UPD_S-1)THEN                        !<-- Child task update region on SE corner of footprint.
            TAIL%IL(1)=I1
            TAIL%IL(2)=I2
            TAIL%IL(3)=IDE_FOOTPRINT-NROWS_P_UPD_E+1
            TAIL%IL(4)=I2
            TAIL%JL(1)=J1
            TAIL%JL(2)=JDS_FOOTPRINT+NROWS_P_UPD_S-1
            TAIL%JL(3)=TAIL%JL(2)+1
            TAIL%JL(4)=J2
!
          ELSEIF(J2>=JDE_FOOTPRINT-NROWS_P_UPD_N+1)THEN                    !<-- Child task update region on NE corner of footprint.
            TAIL%IL(1)=IDE_FOOTPRINT-NROWS_P_UPD_E+1
            TAIL%IL(2)=I2
            TAIL%IL(3)=I1
            TAIL%IL(4)=I2
            TAIL%JL(1)=J1
            TAIL%JL(2)=JDE_FOOTPRINT-NROWS_P_UPD_N
            TAIL%JL(3)=TAIL%JL(2)+1
            TAIL%JL(4)=J2
          ENDIF
!
        ENDIF
!
!-----------------------------------------------------------------------
!
      ENDIF parent_updates
!
!-----------------------------------------------------------------------
!
      TAIL%NUM_PTS_UPDATE_HZ=(TAIL%IL(2)-TAIL%IL(1)+1)                  &
                            *(TAIL%JL(2)-TAIL%JL(1)+1)
!
      IF(TAIL%IL(3)>0)THEN
        TAIL%NUM_PTS_UPDATE_HZ=(TAIL%IL(4)-TAIL%IL(3)+1)                &
                              *(TAIL%JL(4)-TAIL%JL(3)+1)                &
                              +TAIL%NUM_PTS_UPDATE_HZ
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE PARENT_FINDS_UPDATE_LIMITS
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE PARENT_BOOKKEEPING_MOVING
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE PARENT_UPDATES_MOVING(FLAG_H_OR_V                      &
                                      ,N_UPDATE_CHILD_TASKS             &
                                      ,PARENT_CHILD_SPACE_RATIO         &
                                      ,PARENT_CHILD_TIME_RATIO          &
                                      ,NTIMESTEP_CHILD                  &
                                      ,I_PARENT_SW                      &
                                      ,J_PARENT_SW                      &
                                      ,PT,PDTOP,PSGML1,SGML2,SG1,SG2    &
                                      ,DSG2,PDSG1                       &
                                      ,FIS,PD                           &
                                      ,T,Q,CW                           &
                                      ,NUM_PARENT_TASKS                 &
                                      ,NUM_CHILD_TASKS                  &
                                      ,CHILD_TASK_RANKS                 &
                                      ,CHILD_TASK_LIMITS                &
                                      ,HYPER_A                          &
                                      ,IMS,IME,JMS,JME                  &
                                      ,IDS,IDE,JDS,JDE                  &
                                      ,NUM_LYRS                         &
                                      ,LBND1,UBND1,LBND2,UBND2          &
                                      ,FIS_CHILD                        &
                                      ,COMM_TO_MY_CHILD                 &
                                      ,HANDLE_UPDATE                    &
                                      ,MOVE_BUNDLE                      &
                                      ,NUM_FIELDS_MOVE_2D_H_I           &
                                      ,NUM_FIELDS_MOVE_2D_X_I           &
                                      ,NUM_FIELDS_MOVE_2D_H_R           &
                                      ,NUM_FIELDS_MOVE_2D_X_R           &
                                      ,NUM_FIELDS_MOVE_3D_H             &
                                      ,NUM_LEVELS_MOVE_3D_H             &
                                      ,NUM_FIELDS_MOVE_2D_V             &
                                      ,NUM_FIELDS_MOVE_3D_V             &
                                      ,NUM_LEVELS_MOVE_3D_V             &
                                      ,TASK_UPDATE_SPECS                &
                                      ,CHILD_UPDATE_DATA                &
                                        )
!
!-----------------------------------------------------------------------
!***  Each parent task knows which moving nest tasks if any that it
!***  must update and which points on those tasks.  Now the bilinear
!***  interpolation weights can be computed and then all specified
!***  2-D and 3-D variables are interpolated from the parent grid
!***  to the nest's.  Finally the parent tasks send the data to the
!***  appropriate nest tasks.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      INTEGER(kind=KINT),INTENT(IN) :: COMM_TO_MY_CHILD                 &  !<-- MPI communicator to the current nest/child
                                      ,I_PARENT_SW,J_PARENT_SW          &  !<-- SW corner of nest on this parent I,J after move
                                      ,IDS,IDE,JDS,JDE                  &  !<-- Parent domain index limits
                                      ,IMS,IME,JMS,JME                  &  !<-- Parent task memory index limits
                                      ,N_UPDATE_CHILD_TASKS             &  !<-- # of moving nest tasks updated by this parent task
                                      ,NTIMESTEP_CHILD                  &  !<-- Child's timestep at which it recvs parent data
                                      ,NUM_LYRS                         &  !<-- # of model layers
                                      ,NUM_CHILD_TASKS                  &  !<-- # of forecast tasks on all of this parent's children
                                      ,NUM_FIELDS_MOVE_2D_H_I           &  !<-- # of 2-D integer H arrays specified for updating
                                      ,NUM_FIELDS_MOVE_2D_X_I           &  !<-- # of 2-D integer H arrays updated from external files
                                      ,NUM_FIELDS_MOVE_2D_H_R           &  !<-- # of 2-D real H arrays specified for updating
                                      ,NUM_FIELDS_MOVE_2D_X_R           &  !<-- # of 2-D real H arrays updated from external files
                                      ,NUM_FIELDS_MOVE_3D_H             &  !<-- # of 3-D H arrays specified for updating
                                      ,NUM_LEVELS_MOVE_3D_H             &  !<-- # of 2-D levels in all 3-D H update arrays
                                      ,NUM_FIELDS_MOVE_2D_V             &  !<-- # of 2-D V arrays specified for updating
                                      ,NUM_FIELDS_MOVE_3D_V             &  !<-- # of 3-D V arrays specified for updating
                                      ,NUM_LEVELS_MOVE_3D_V             &  !<-- # of 2-D levels in all 3-D V update arrays
                                      ,NUM_PARENT_TASKS                 &  !<-- # of forecast tasks on this parent
                                      ,PARENT_CHILD_SPACE_RATIO         &  !<-- Ratio of parent's grid increment to its child's
                                      ,PARENT_CHILD_TIME_RATIO          &  !<-- Ratio of parent's time step to its child's
                                      ,LBND1,UBND1,LBND2,UBND2             !<-- Array bounds of nest-resolution FIS on parent
!
      INTEGER(kind=KINT),DIMENSION(:),POINTER,INTENT(IN) ::             &
                                                  HANDLE_UPDATE            !<-- MPI Handles for ISends to the child tasks
!
      INTEGER(kind=KINT),DIMENSION(1:NUM_CHILD_TASKS),INTENT(IN) ::     &
                                                  CHILD_TASK_RANKS         !<-- Child task local ranks in p-c intracomm
!
      INTEGER(kind=KINT),DIMENSION(1:4,NUM_CHILD_TASKS),INTENT(IN) ::   &
                                                     CHILD_TASK_LIMITS     !<-- ITS,ITE,JTS,JTE for each child forecast task
!
      REAL(kind=KFPT),INTENT(IN) :: PDTOP                               &  !<-- Pressure at top of sigma domain (Pa)
                                   ,PT                                     !<-- Top pressure of model domain (Pa)
!
      REAL(kind=KDBL),INTENT(IN) :: HYPER_A                                !<-- Underground extrapolation quantity
!
      REAL(kind=KFPT),DIMENSION(1:NUM_LYRS),INTENT(IN) :: DSG2          &  !<-- Vertical structure coefficients for midlayers
                                                         ,PDSG1         &  !
                                                         ,PSGML1        &  !
                                                         ,SGML2            !<--
!
      REAL(kind=KFPT),DIMENSION(1:NUM_LYRS+1),INTENT(IN) :: SG1,SG2        !<-- Vertical structure coefficients for interfaces
!
      REAL(kind=KFPT),DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: FIS      &  !<-- Sfc geopotential on parent mass points
                                                              ,PD          !<-- Parent PD
!
      REAL(kind=KFPT),DIMENSION(LBND1:UBND1,LBND2:UBND2),INTENT(IN) ::  &
                                                             FIS_CHILD     !<-- Moving nest's full res FIS distributed on the parent
!
      REAL(kind=KFPT),DIMENSION(IMS:IME,JMS:JME,1:NUM_LYRS)             &
                                                      ,INTENT(IN) :: T  &  !<-- Parent sensible temperature (K)
                                                                    ,Q  &  !<-- Parent specific humidity (kg/kg)
                                                                    ,CW    !<-- Parent cloud condensate (kg/kg)
!
      CHARACTER(len=1),INTENT(IN) :: FLAG_H_OR_V                           !<-- Are we updating H or V points?
!
      TYPE(MIXED_DATA_TASKS),INTENT(INOUT) :: CHILD_UPDATE_DATA            !<-- Composite of all update data from parent for each nest task
!
      TYPE(CHILD_UPDATE_LINK),TARGET,INTENT(INOUT) ::                   &
                                                  TASK_UPDATE_SPECS        !<-- Linked list with nest task update specifications
!
      TYPE(ESMF_FieldBundle),INTENT(INOUT) :: MOVE_BUNDLE                  !<-- ESMF Bundle of 2-D and 3-D arrays specified for updating
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT),SAVE :: I,I1,I2                                &
                                ,I_EAST,I_OFFSET,I_WEST                 &
                                ,IDS_CHILD,ISTART,ITAG,ITER             &
                                ,J,J1,J2                                &
                                ,J_NORTH,J_OFFSET,J_SOUTH               &
                                ,JDS_CHILD,JSTART                       &
                                ,KHI,KLO                                &
                                ,L,LOC_1,LOC_2                          &
                                ,N,N_ADD,N_FIELD,N_REMOVE,N_STRIDE      &
                                ,NPOINTS_HORIZ                          &
                                ,NPOINTS_HORIZ_H                        &
                                ,NPOINTS_HORIZ_V                        &
                                ,NUM_DIMS                               &
                                ,NUM_FIELDS_MOVE                        &
                                ,NUM_LEVS_IN                            &
                                ,NUM_LEVS_SEC                           &
                                ,NUM_LEVELS                             &
                                ,NUM_INTEGER_WORDS_SEND                 &
                                ,NUM_REAL_WORDS_SEND                    &
                                ,UPDATE_TYPE_INT
!
      INTEGER(kind=KINT) :: CHILDTASK,I_TRANS,IERR,ISTAT,IVAL,J_TRANS   &
                           ,KNT_DUMMY,RC,RC_UPDATE
!
      INTEGER(kind=KINT),DIMENSION(1:3) :: LIMITS_HI                    &
                                          ,LIMITS_LO
!
      INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE,TARGET ::             &
                                                         I_PARENT_EAST  &
                                                        ,I_PARENT_WEST  &
                                                        ,J_PARENT_NORTH &
                                                        ,J_PARENT_SOUTH  
!
      INTEGER(kind=KINT),DIMENSION(:),POINTER :: I_PARENT_EAST_H        &
                                                ,I_PARENT_WEST_H        &
                                                ,J_PARENT_NORTH_H       &
                                                ,J_PARENT_SOUTH_H
!
      INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE,SAVE ::               &
                                                       KNT_INTEGER_PTS  &
                                                      ,KNT_REAL_PTS     &
                                                      ,NUM_ITER
!
      INTEGER(kind=KINT),DIMENSION(:,:),POINTER :: IARRAY_2D
!
      REAL(kind=KFPT) :: CHILD_PARENT_SPACE_RATIO                       &
                        ,COEFF_1,COEFF_2,CW_INTERP                      &
                        ,D_LNP_DFI,DELP_EXTRAP,DP,FACTOR                &
                        ,IDIFF_EAST,IDIFF_WEST                          &
                        ,JDIFF_NORTH,JDIFF_SOUTH                        &
                        ,LOG_P1_PARENT                                  &
                        ,MAX_WGHT                                       &
                        ,PDTOP_PT,PHI_DIFF                              &
                        ,PSFC_CHILD                                     &
                        ,PSFC_PARENT_NE,PSFC_PARENT_NW                  &
                        ,PSFC_PARENT_SE,PSFC_PARENT_SW                  &
                        ,PX_NE,PX_NW,PX_SE,PX_SW                        &
                        ,Q_INTERP,R_DELP,R_INC                          &
                        ,RECIP_SUM_WGT,SUM_PROD,SUM_WGT                 &
                        ,T_INTERP,TMP                                   &
                        ,X_NE,X_NW,X_SE,X_SW
!
      REAL(kind=KFPT),DIMENSION(:),ALLOCATABLE :: P_OUTPUT
!
      REAL(kind=KFPT),DIMENSION(1:NUM_LYRS+2) :: P_INPUT                &
                                                ,VBL_INPUT
!
      REAL(kind=KFPT),DIMENSION(:),ALLOCATABLE :: I_PARENT              &
                                                 ,J_PARENT              &
                                                 ,SEC_DERIV
!
      REAL(kind=KFPT),DIMENSION(:),POINTER :: I_CHILD_ON_PARENT_H       &
                                             ,J_CHILD_ON_PARENT_H
!
      REAL(kind=KFPT),DIMENSION(:),POINTER :: VBL_COL_CHILD             &
                                             ,VBL_COL_X
!
      REAL(kind=KFPT),DIMENSION(IMS:IME,JMS:JME) :: LMASK
!
      REAL(kind=KFPT),DIMENSION(1:NUM_LYRS+2,1:4) :: C_TMP                 !<-- Working array for ESSL spline call
!
      REAL(kind=KFPT),DIMENSION(:,:),ALLOCATABLE :: LOG_PBOT            &
                                                   ,LOG_PTOP            &
                                                   ,PD_CHILD            &
                                                   ,PD_INTERP           &
                                                   ,PROD_LWGT_NE        &
                                                   ,PROD_LWGT_NW        &
                                                   ,PROD_LWGT_SE        &
                                                   ,PROD_LWGT_SW        &
                                                   ,PROD_SWGT_NE        &
                                                   ,PROD_SWGT_NW        &
                                                   ,PROD_SWGT_SE        &
                                                   ,PROD_SWGT_SW
!
      REAL(kind=KFPT),DIMENSION(:,:),ALLOCATABLE,TARGET :: WGHT_NE      &
                                                          ,WGHT_NW      &
                                                          ,WGHT_SE      &
                                                          ,WGHT_SW
!
      REAL(kind=KFPT),DIMENSION(:,:),POINTER,SAVE :: PDO                &
                                                    ,SMASK
!
      REAL(kind=KFPT),DIMENSION(:,:),POINTER :: ARRAY_2D                &
                                               ,WGHT_NE_H               &
                                               ,WGHT_NW_H               &
                                               ,WGHT_SE_H               &
                                               ,WGHT_SW_H
!
      REAL(kind=KFPT),DIMENSION(:,:,:),ALLOCATABLE,TARGET :: PINT_CHILD  &
                                                            ,PINT_INTERP &
                                                            ,PMID_CHILD  &
                                                            ,PMID_INTERP
!
      REAL(kind=KFPT),DIMENSION(:,:,:),ALLOCATABLE :: PHI_INTERP        &
                                                     ,VBL_INTERP
!
      REAL(kind=KFPT),DIMENSION(:,:,:),POINTER :: ARRAY_3D              &
                                                 ,P3D_INPUT             &
                                                 ,P3D_OUTPUT
!
      LOGICAL(kind=KLOG) :: INTERFACES                                  &
                           ,MIDLAYERS
!
      CHARACTER(len=1) :: UPDATE_TYPE_CHAR
!
      CHARACTER(len=4) :: FNAME
!
      CHARACTER(len=30) :: FIELD_NAME
!
      TYPE(CHILD_UPDATE_LINK),POINTER,SAVE :: PTR_H,PTR_V
!
      TYPE(CHILD_UPDATE_LINK),POINTER :: PTR_X
!
      TYPE(ESMF_Field) :: HOLD_FIELD
!
      TYPE(ESMF_TypeKind_Flag) :: DATATYPE
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      btim=timef()
!
      RC       =ESMF_SUCCESS
      RC_UPDATE=ESMF_SUCCESS
!
      CHILD_PARENT_SPACE_RATIO=1./PARENT_CHILD_SPACE_RATIO
!
!-----------------------------------------------------------------------
!***  This update routine is called first for H points and then 
!***  a second time for V points.  To save time on communication
!***  all H and V point data will be sent together at the end of
!***  the 2nd (V-point) call.  First do some prep work that only 
!***  needs to be done once for both H and V at this move.
!-----------------------------------------------------------------------
!
      prep_block: IF(FLAG_H_OR_V=='H')THEN  
!
!-----------------------------------------------------------------------
!
        ALLOCATE(KNT_REAL_PTS(1:N_UPDATE_CHILD_TASKS)                   &
                                                     ,stat=ISTAT)
        ALLOCATE(KNT_INTEGER_PTS(1:N_UPDATE_CHILD_TASKS)                &
                                                     ,stat=ISTAT)
        ALLOCATE(CHILD_UPDATE_DATA%TASKS(1:N_UPDATE_CHILD_TASKS)        &
                                                     ,stat=ISTAT)
        ALLOCATE(NUM_ITER(1:N_UPDATE_CHILD_TASKS)                       &
                                                     ,stat=ISTAT)
!
        LM=NUM_LYRS
!
!-----------------------------------------------------------------------
!***  Start at the top of the linked lists that hold the task ID
!***  and index limits for all update H and V points on each nest
!***  task for the current nest.  Remember that each link in the
!***  lists corresponds to a nest task that this parent task must
!***  update.
!-----------------------------------------------------------------------
!
        PTR_H=>TASK_UPDATE_SPECS
        PTR_V=>TASK_UPDATE_SPECS
!
!-----------------------------------------------------------------------
!***  Find the total number of words to be updated on each nest task
!***  for both H and V points.
!-----------------------------------------------------------------------
!
        prep_loop: DO N=1,N_UPDATE_CHILD_TASKS
!
!-----------------------------------------------------------------------
!
          IF(N>1)THEN                                                      !<-- Point to the next link (the next task to be updated).
            PTR_H=>PTR_H%NEXT_LINK
            PTR_V=>PTR_V%NEXT_LINK
          ENDIF
!
          NPOINTS_HORIZ_H=(PTR_H%IL(2)-PTR_H%IL(1)+1)                   &
                         *(PTR_H%JL(2)-PTR_H%JL(1)+1)
!
          NPOINTS_HORIZ_V=(PTR_V%IL(2)-PTR_V%IL(1)+1)                   &
                         *(PTR_V%JL(2)-PTR_V%JL(1)+1)
!
          NUM_INTEGER_WORDS_SEND=(NUM_FIELDS_MOVE_2D_H_I                &
                                 -NUM_FIELDS_MOVE_2D_X_I)               &
                                 *NPOINTS_HORIZ_H
!
          NUM_REAL_WORDS_SEND=(NUM_FIELDS_MOVE_2D_H_R                   &
                              -NUM_FIELDS_MOVE_2D_X_R                   &
                              +NUM_LEVELS_MOVE_3D_H)                    &
                              *NPOINTS_HORIZ_H                          &
                             +(NUM_FIELDS_MOVE_2D_V                     &
                              +NUM_LEVELS_MOVE_3D_V)                    &
                              *NPOINTS_HORIZ_V
!
!-----------------------------------------------------------------------
!***  If there is a 2nd region on nest task N updated by the current
!***  parent task then we need to iterate twice through the updating
!***  process.  These 2nd regions exist only when the parent task 
!***  and the nest task it is updating both lie on the corner of the
!***  nest's pre-move footprint.
!-----------------------------------------------------------------------
!
          NUM_ITER(N)=1
!
          IF(PTR_H%IL(3)>0)THEN                                            !<-- If true then there must be a 2nd update region.               
            IF(PTR_V%IL(3)<0)THEN
              WRITE(0,*)' A 2nd update region exists for H points but not V!!  ABORT!!'
              CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
            ENDIF
            NUM_ITER(N)=2                           
            NPOINTS_HORIZ_H=(PTR_H%IL(4)-PTR_H%IL(3)+1)                 &
                           *(PTR_H%JL(4)-PTR_H%JL(3)+1)
!
            NPOINTS_HORIZ_V=(PTR_V%IL(4)-PTR_V%IL(3)+1)                 &
                           *(PTR_V%JL(4)-PTR_V%JL(3)+1)
!
            NUM_INTEGER_WORDS_SEND=(NUM_FIELDS_MOVE_2D_H_I              &  !<-- Total # of integer words in parent's update 
                                   -NUM_FIELDS_MOVE_2D_X_I)             &  !    of nest task N.
                                   *NPOINTS_HORIZ_H                     &
                                   +NUM_INTEGER_WORDS_SEND
!
            NUM_REAL_WORDS_SEND=(NUM_FIELDS_MOVE_2D_H_R                 &  !<-- Total # of real words in parent's update
                                -NUM_FIELDS_MOVE_2D_X_R                 &  !    of nest task N.
                                +NUM_LEVELS_MOVE_3D_H)                  &
                                *NPOINTS_HORIZ_H                        &
                               +(NUM_FIELDS_MOVE_2D_V                   &
                                +NUM_LEVELS_MOVE_3D_V)                  &
                                *NPOINTS_HORIZ_V                        &
                                +NUM_REAL_WORDS_SEND
          ENDIF                                 
!                                            
!-----------------------------------------------------------------------
!***  Now we know how many words will be sent from the current parent 
!***  task to nest task N so allocate the objects that will hold this
!***  data.  There may or may not be any integer variables updated at
!***  this time.
!-----------------------------------------------------------------------
!
          ALLOCATE(CHILD_UPDATE_DATA%TASKS(N)%DATA_REAL(1:NUM_REAL_WORDS_SEND)) 
!
          IF(NUM_INTEGER_WORDS_SEND>0)THEN
            ALLOCATE(CHILD_UPDATE_DATA%TASKS(N)%DATA_INTEGER(1:NUM_INTEGER_WORDS_SEND))
          ELSE
            CHILD_UPDATE_DATA%TASKS(N)%DATA_INTEGER=>NULL()
          ENDIF
!
!-----------------------------------------------------------------------
!
          KNT_REAL_PTS(N)=0                                                !<-- Initialize the counter of real update data words.
          KNT_INTEGER_PTS(N)=0                                             !<-- Initialize the counter of integer update data words.
!
          ISTART=MAX(IMS,IDS)
          JSTART=MAX(JMS,JDS)
!
          I_OFFSET=(I_PARENT_SW-ISTART)*PARENT_CHILD_SPACE_RATIO        &  !<-- I offset of child SW corner in full topo array on parent
                   +LBND1-1
          J_OFFSET=(J_PARENT_SW-JSTART)*PARENT_CHILD_SPACE_RATIO        &  !<-- J offset of child SW corner in full topo array on parent
                   +LBND2-1
!
!-----------------------------------------------------------------------
!
        ENDDO prep_loop
!
!-----------------------------------------------------------------------
!***  We need PD and PDO on the parent for building the
!***  pressure structure from which to interpolate to the nest
!***  update points.  PD was already sent into this routine via
!***  the argument list since it was needed earlier in the coupler.
!***  Unload PDO now.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract PDO Field From H Move_Bundle"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_FieldBundleGet(FIELDBUNDLE=MOVE_BUNDLE                &  !<-- Bundle holding the arrays for move updates
                                ,fieldName  ='PDO'//SUFFIX_MOVE         &  !<-- Get the Field with this name
                                ,field      =HOLD_FIELD                 &  !<-- Put the Field here
                                ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_UPDATE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract PDO Array from Field"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_FieldGet(field    =HOLD_FIELD                         &  !<-- Field holding PDO
                          ,localDe  =0                                  &
                          ,farrayPtr=PDO                                &  !<-- Put array here
                          ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_UPDATE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  We need the parent's Sea Mask for generating surface variable
!***  updates in order to exclude either sea or land point values 
!***  in the bilinear interpolation.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract the Sea Mask from the H Move_Bundle"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_FieldBundleGet(FIELDBUNDLE=MOVE_BUNDLE                &  !<-- Bundle holding the arrays for move updates
                                ,fieldName  ='SM'//SUFFIX_MOVE          &  !<-- The parent's sea mask
                                ,field      =HOLD_FIELD                 &  !<-- Put the Field here
                                ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_UPDATE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract Sea Mask Array from Field"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_FieldGet(field    =HOLD_FIELD                           &  !<-- Field holding PDO
                          ,localDe  =0                                    &
                          ,farrayPtr=SMASK                                &  !<-- Put the sea mask array here
                          ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_UPDATE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  'Flip' the seamask values (1=>sea) for use as a landmask (0=>sea).
!-----------------------------------------------------------------------
!
        DO J=JMS,JME
        DO I=IMS,IME
          IF(SMASK(I,J)>0.5)THEN
            LMASK(I,J)=0.
          ELSE
            LMASK(I,J)=1.
          ENDIF
        ENDDO
        ENDDO
!
!-----------------------------------------------------------------------
!
      ENDIF prep_block
!
!-----------------------------------------------------------------------
!***  As we prepare to interpolate from the parent to the child,
!***  we need to be aware of another distinction between the H
!***  and the V point locations on those domains' grids regarding
!***  I=1 on the nest relative to I_PARENT_SW.  For H points
!***  I=1 on the nest coincides with I_PARENT_SW on the parent.
!***  However for V points I=1 is to the west of I_PARENT_SW.
!***  See the diagram at the beginning of the FLAG_H_OR_V==V
!***  section of PARENT_BOOKKEEPING_MOVING.  Specifically
!***  v(1) on the nest is
!***  (0.5*PARENT_CHILD_SPACE_RATIO-0.5)/PARENT_CHILD_SPACE_RATIO
!***  to the west of V(I_PARENT_SW) on the parent grid.  Compute
!***  that increment here then use it below when we need the
!***  Real values for parent I's and J's on the parent grid
!***  that coincide with the update locations on the nest grid.
!-----------------------------------------------------------------------
!
      PTR_X=>TASK_UPDATE_SPECS
!
      IF(FLAG_H_OR_V=='H')THEN
        NUM_FIELDS_MOVE=NUM_FIELDS_MOVE_2D_H_I                          &
                       +NUM_FIELDS_MOVE_2D_H_R                          &
                       +NUM_FIELDS_MOVE_3D_H
        R_INC=0.
!
      ELSEIF(FLAG_H_OR_V=='V')THEN
        NUM_FIELDS_MOVE=NUM_FIELDS_MOVE_2D_V                            &
                       +NUM_FIELDS_MOVE_3D_V
        R_INC=-(0.5*PARENT_CHILD_SPACE_RATIO-0.5)                       &
               *CHILD_PARENT_SPACE_RATIO
      ENDIF
!
!-----------------------------------------------------------------------
!
      DO L=1,NUM_LYRS+2
        P_INPUT(L)=0.
        VBL_INPUT(L)=0.
      ENDDO
!
!-----------------------------------------------------------------------
!***  Loop through each of the moving nest tasks whose subdomains 
!***  contain points that must be updated by this parent task
!***  after the nest moved.
!-----------------------------------------------------------------------
!
      ctask_loop: DO N=1,N_UPDATE_CHILD_TASKS
!
!-----------------------------------------------------------------------
!
        iter_loop: DO ITER=1,NUM_ITER(N)                                   !<-- Either one or two regions on the nest task must be updated.
!
!-----------------------------------------------------------------------
!
          IF(ITER==1)THEN
            I1=PTR_X%IL(1)                                                 !<-- I limits of nest task's update region by parent task
            I2=PTR_X%IL(2)                                                 !     in terms of the nest's grid.
            J1=PTR_X%JL(1)                                                 !<-- J limits of nest task's update region by parent task
            J2=PTR_X%JL(2)                                                 !     in terms of the nest's grid.
          ELSE
            I1=PTR_X%IL(3)                                                 !<-- I limits of nest task's update region by parent task
            I2=PTR_X%IL(4)                                                 !     in terms of the nest's grid for 2nd update region.
            J1=PTR_X%JL(3)                                                 !<-- J limits of nest task's update region by parent task
            J2=PTR_X%JL(4)                                                 !     in terms of the nest's grid for 2nd update region.
          ENDIF
!
          ALLOCATE(I_PARENT(I1:I2))
          ALLOCATE(I_PARENT_EAST(I1:I2))
          ALLOCATE(I_PARENT_WEST(I1:I2))
!
          ALLOCATE(J_PARENT(J1:J2))
          ALLOCATE(J_PARENT_NORTH(J1:J2))
          ALLOCATE(J_PARENT_SOUTH(J1:J2))
!
          ALLOCATE(WGHT_SW(I1:I2,J1:J2))
          ALLOCATE(WGHT_NW(I1:I2,J1:J2))
          ALLOCATE(WGHT_NE(I1:I2,J1:J2))
          ALLOCATE(WGHT_SE(I1:I2,J1:J2))
!
          ALLOCATE(PINT_INTERP(I1:I2,J1:J2,1:NUM_LYRS+1))
          ALLOCATE( PHI_INTERP(I1:I2,J1:J2,1:NUM_LYRS+1))
          ALLOCATE(   PD_CHILD(I1:I2,J1:J2))
          ALLOCATE(  PD_INTERP(I1:I2,J1:J2))
          ALLOCATE(   LOG_PBOT(I1:I2,J1:J2))
          ALLOCATE(   LOG_PTOP(I1:I2,J1:J2))
!
          NPOINTS_HORIZ=(I2-I1+1)*(J2-J1+1)
          ALLOCATE(PMID_INTERP(I1:I2,J1:J2,1:NUM_LYRS))
          ALLOCATE( PMID_CHILD(I1:I2,J1:J2,1:NUM_LYRS))
          ALLOCATE( PINT_CHILD(I1:I2,J1:J2,1:NUM_LYRS+1))
!
          IDS_CHILD=CHILD_TASK_LIMITS(1,1)                                 !<-- Child task's starting I on grid of moving nest
          JDS_CHILD=CHILD_TASK_LIMITS(3,1)                                 !<-- Child task's starting J on grid of moving nest
!
          DO I=I1,I2
            I_PARENT(I)=I_PARENT_SW+R_INC                               &  !<-- Real Parent I's on parent grid for these nest I's
                        +(I-IDS_CHILD)*CHILD_PARENT_SPACE_RATIO            !    in the nest grid's Update region.
          ENDDO
!
          DO J=J1,J2
            J_PARENT(J)=J_PARENT_SW+R_INC                               &  !<-- Real Parent J's on parent grid for these nest J's
                        +(J-JDS_CHILD)*CHILD_PARENT_SPACE_RATIO            !    in the nest grid's update region.
          ENDDO
!
!-----------------------------------------------------------------------
!***  Loop through this nest's update points and determine the four
!***  parent points that surround each nest point as well and the
!***  bilinear interpolation weight associated with each of those
!***  four parent points for the given nest point.
!-----------------------------------------------------------------------
!
          DO J=J1,J2
            J_PARENT_SOUTH(J)=INT(J_PARENT(J)+EPS)                         !<-- Parent J at or immediately south of nest point
            J_PARENT_NORTH(J)=J_PARENT_SOUTH(J)+1                          !<-- Parent J immediately north of nest point
!
            DO I=I1,I2
              I_PARENT_WEST(I)=INT(I_PARENT(I)+EPS)                        !<-- Parent I at or immediately west of nest point
              I_PARENT_EAST(I)=I_PARENT_WEST(I)+1                          !<-- Parent I immediately east of nest point
!
              IDIFF_EAST=I_PARENT_EAST(I)-I_PARENT(I)
              IDIFF_WEST=I_PARENT(I)-I_PARENT_WEST(I)
              JDIFF_NORTH=J_PARENT_NORTH(J)-J_PARENT(J)
              JDIFF_SOUTH=J_PARENT(J)-J_PARENT_SOUTH(J)
!
              WGHT_SW(I,J)=IDIFF_EAST*JDIFF_NORTH                          !<-- Bilinear weight for parent's point SW of child's point
              WGHT_NW(I,J)=IDIFF_EAST*JDIFF_SOUTH                          !<-- Bilinear weight for parent's point NW of child's point
              WGHT_NE(I,J)=IDIFF_WEST*JDIFF_SOUTH                          !<-- Bilinear weight for parent's point NE of child's point
              WGHT_SE(I,J)=IDIFF_WEST*JDIFF_NORTH                          !<-- Bilinear weight for parent's point SE of child's point
!
            ENDDO
          ENDDO
!
!-----------------------------------------------------------------------
!***  The parent computes its layer interface pressures at the 
!***  locations of the moving nest update points.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  If we are updating mass point variables then those variables
!***  obviously coincide with the pressure information.  In other
!***  words we are interpolating from parent H points to nest h points.
!-----------------------------------------------------------------------
!
          h_v_block: IF(FLAG_H_OR_V=='H')THEN
!
            I_PARENT_EAST_H=>I_PARENT_EAST
            I_PARENT_WEST_H=>I_PARENT_WEST
            J_PARENT_NORTH_H=>J_PARENT_NORTH
            J_PARENT_SOUTH_H=>J_PARENT_SOUTH
!
            WGHT_NE_H=>WGHT_NE
            WGHT_NW_H=>WGHT_NW
            WGHT_SE_H=>WGHT_SE
            WGHT_SW_H=>WGHT_SW
!
!-----------------------------------------------------------------------
!***  If we are updating wind components then we need to know
!***  T, Q, and FIS at V points in order to compute sfc pressure
!***  and ultimately midlayer pressure at the V points.
!***  Base the bilinear interpolation to nest v points on the 
!***  values at parent H points (where T, Q, and FIS are defined)
!***  in order to minimize horizontal interpolation.
!-----------------------------------------------------------------------
!
          ELSEIF(FLAG_H_OR_V=='V')THEN
!
            ALLOCATE(I_PARENT_EAST_H(I1:I2))
            ALLOCATE(I_PARENT_WEST_H(I1:I2))
            ALLOCATE(J_PARENT_NORTH_H(J1:J2))
            ALLOCATE(J_PARENT_SOUTH_H(J1:J2))
!
            ALLOCATE(I_CHILD_ON_PARENT_H(I1:I2))
            ALLOCATE(J_CHILD_ON_PARENT_H(J1:J2))
!
            ALLOCATE(WGHT_NE_H(I1:I2,J1:J2))
            ALLOCATE(WGHT_NW_H(I1:I2,J1:J2))
            ALLOCATE(WGHT_SE_H(I1:I2,J1:J2))
            ALLOCATE(WGHT_SW_H(I1:I2,J1:J2))
!
            DO I=I1,I2     
!
              I_CHILD_ON_PARENT_H(I)=I_PARENT_SW                        &
                                   +(I-IDS_CHILD+0.5)                   &
                                    *CHILD_PARENT_SPACE_RATIO
!
              I_PARENT_WEST_H(I)=INT(I_CHILD_ON_PARENT_H(I))               !<-- Parent I on H immediately west of nest V point
              I_PARENT_EAST_H(I)=I_PARENT_WEST_H(I)+1                      !<-- Parent I on H immediately east of nest V point
!
            ENDDO    
!
            DO J=J1,J2     
!
              J_CHILD_ON_PARENT_H(J)=J_PARENT_SW                        &
                                   +(J-JDS_CHILD+0.5)                   &
                                    *CHILD_PARENT_SPACE_RATIO
!
              J_PARENT_SOUTH_H(J)=INT(J_CHILD_ON_PARENT_H(J))              !<-- Parent J on H immediately south of nest V point
              J_PARENT_NORTH_H(J)=J_PARENT_SOUTH_H(J)+1                    !<-- Parent J on H immediately north of nest V point
!
            ENDDO
!
            DO J=J1,J2
            DO I=I1,I2
              WGHT_SW_H(I,J)=(I_PARENT_EAST_H(I)-I_CHILD_ON_PARENT_H(I)) &
                            *(J_PARENT_NORTH_H(J)-J_CHILD_ON_PARENT_H(J))
              WGHT_SE_H(I,J)=(I_CHILD_ON_PARENT_H(I)-I_PARENT_WEST_H(I)) &
                            *(J_PARENT_NORTH_H(J)-J_CHILD_ON_PARENT_H(J))
              WGHT_NW_H(I,J)=(I_PARENT_EAST_H(I)-I_CHILD_ON_PARENT_H(I)) &
                            *(J_CHILD_ON_PARENT_H(J)-J_PARENT_SOUTH_H(J))
              WGHT_NE_H(I,J)=(I_CHILD_ON_PARENT_H(I)-I_PARENT_WEST_H(I)) &
                            *(J_CHILD_ON_PARENT_H(J)-J_PARENT_SOUTH_H(J))
            ENDDO
            ENDDO
!
          ENDIF h_v_block
!
!-----------------------------------------------------------------------
!***  When the parent generates Real soil variable updates for its
!***  moving nests it uses bilinear interpolation but also must use
!***  the sea/land mask in order to avoid including sea values in
!***  land variables and vice versa.  This means the bilinear
!***  interpolation weighting needs to be adjusted to account for
!***  the exclusion of sea or land points in the 4-pt summation.
!-----------------------------------------------------------------------
!
          soil_wgts: IF(FLAG_H_OR_V=='H')THEN
!
            ALLOCATE(PROD_LWGT_SW(I1:I2,J1:J2))
            ALLOCATE(PROD_LWGT_SE(I1:I2,J1:J2))
            ALLOCATE(PROD_LWGT_NW(I1:I2,J1:J2))
            ALLOCATE(PROD_LWGT_NE(I1:I2,J1:J2))
            ALLOCATE(PROD_SWGT_SW(I1:I2,J1:J2))
            ALLOCATE(PROD_SWGT_SE(I1:I2,J1:J2))
            ALLOCATE(PROD_SWGT_NW(I1:I2,J1:J2))
            ALLOCATE(PROD_SWGT_NE(I1:I2,J1:J2))
!
            DO J=J1,J2
              J_SOUTH=J_PARENT_SOUTH(J)
              J_NORTH=J_PARENT_NORTH(J)
!
              DO I=I1,I2
                I_WEST=I_PARENT_WEST(I)
                I_EAST=I_PARENT_EAST(I)
!
                X_SW=WGHT_SW(I,J)*LMASK(I_WEST,J_SOUTH)
                X_SE=WGHT_SE(I,J)*LMASK(I_EAST,J_SOUTH)
                X_NW=WGHT_NW(I,J)*LMASK(I_WEST,J_NORTH)
                X_NE=WGHT_NE(I,J)*LMASK(I_EAST,J_NORTH)
!
                SUM_WGT=X_SW+X_SE+X_NW+X_NE
!
                IF(ABS(SUM_WGT)>1.E-6)THEN
                  RECIP_SUM_WGT=1./(X_SW+X_SE+X_NW+X_NE)
                ELSE
                  RECIP_SUM_WGT=0.
                ENDIF
!
                PROD_LWGT_SW(I,J)=X_SW*RECIP_SUM_WGT                       !<-- These are the adjusted bilinear interpolation
                PROD_LWGT_SE(I,J)=X_SE*RECIP_SUM_WGT                       !    weights that take into account the presence
                PROD_LWGT_NW(I,J)=X_NW*RECIP_SUM_WGT                       !    of sea points that must be excluded in the
                PROD_LWGT_NE(I,J)=X_NE*RECIP_SUM_WGT                       !    summation.
!
                X_SW=WGHT_SW(I,J)*SMASK(I_WEST,J_SOUTH)
                X_SE=WGHT_SE(I,J)*SMASK(I_EAST,J_SOUTH)
                X_NW=WGHT_NW(I,J)*SMASK(I_WEST,J_NORTH)
                X_NE=WGHT_NE(I,J)*SMASK(I_EAST,J_NORTH)
!
                SUM_WGT=X_SW+X_SE+X_NW+X_NE
!
                IF(ABS(SUM_WGT)>1.E-6)THEN
                  RECIP_SUM_WGT=1./(X_SW+X_SE+X_NW+X_NE)
                ELSE
                  RECIP_SUM_WGT=0.
                ENDIF
!
                PROD_SWGT_SW(I,J)=X_SW*RECIP_SUM_WGT                       !<-- These are the adjusted bilinear interpolation
                PROD_SWGT_SE(I,J)=X_SE*RECIP_SUM_WGT                       !    weights that take into account the presence
                PROD_SWGT_NW(I,J)=X_NW*RECIP_SUM_WGT                       !    of land points that must be excluded in the
                PROD_SWGT_NE(I,J)=X_NE*RECIP_SUM_WGT                       !    summation.
              ENDDO
            ENDDO
!
          ENDIF soil_wgts
!
!-----------------------------------------------------------------------
!***  Some of the primary dynamics integration variables are valid
!***  at the previous time step (PDO,TP,UP,VP).  Update those at
!***  the appropriate nest gridpoints as well.  Since the parent's
!***  time step is larger than the child's, approximate the child's
!***  previous time step value as
!***  (PARENT_CHILD_TIME_RATIO-1.)/PARENT_CHILD_TIME_RATIO 
!***  between the old and current parent values.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Compute the parent's Psfc at the nest H or V update points.
!-----------------------------------------------------------------------
!
          DO J=J1,J2
          DO I=I1,I2
!
            PSFC_PARENT_SW=PD(I_PARENT_WEST_H(I),J_PARENT_SOUTH_H(J))+PT
            PSFC_PARENT_SE=PD(I_PARENT_EAST_H(I),J_PARENT_SOUTH_H(J))+PT
            PSFC_PARENT_NW=PD(I_PARENT_WEST_H(I),J_PARENT_NORTH_H(J))+PT
            PSFC_PARENT_NE=PD(I_PARENT_EAST_H(I),J_PARENT_NORTH_H(J))+PT
!
            PINT_INTERP(I,J,LM+1)=WGHT_SW_H(I,J)*PSFC_PARENT_SW          &  !<-- Parent's Psfc at nest point at parent's sfc elevation
                                 +WGHT_SE_H(I,J)*PSFC_PARENT_SE          &  !
                                 +WGHT_NW_H(I,J)*PSFC_PARENT_NW          &  !
                                 +WGHT_NE_H(I,J)*PSFC_PARENT_NE             !<--
!
            LOG_PBOT(I,J)=LOG(PINT_INTERP(I,J,LM+1))
!
            PHI_INTERP(I,J,LM+1)=WGHT_SW_H(I,J)*FIS(I_PARENT_WEST_H(I)   &  !<-- Parent's sfc geopotential at nest point
                                                   ,J_PARENT_SOUTH_H(J)) &  !
                                +WGHT_SE_H(I,J)*FIS(I_PARENT_EAST_H(I)   &  !
                                                   ,J_PARENT_SOUTH_H(J)) &  !
                                +WGHT_NW_H(I,J)*FIS(I_PARENT_WEST_H(I)   &  !
                                                   ,J_PARENT_NORTH_H(J)) &  !
                                +WGHT_NE_H(I,J)*FIS(I_PARENT_EAST_H(I)   &  !
                                                   ,J_PARENT_NORTH_H(J))    !<--
!
          ENDDO
          ENDDO
!
!-----------------------------------------------------------------------
!***  Parent computes its layer interface pressures and geopotentials
!***  at the locations of the moving nest update points.  The input
!***  and target values of pressure locations for the vertical
!***  interpolations from parent to nest are the hydrostatic midlayer
!***  and interface pressures.
!-----------------------------------------------------------------------
!
          DO J=J1,J2                                                       !<-- J limits of child task update region on parent task
            J_SOUTH=J_PARENT_SOUTH_H(J)
            J_NORTH=J_PARENT_NORTH_H(J)
!
            DO I=I1,I2                                                     !<-- I limits of child task update region on parent task
              I_WEST=I_PARENT_WEST_H(I)
              I_EAST=I_PARENT_EAST_H(I)
!
              PD_INTERP(I,J)=WGHT_SW_H(I,J)*PD(I_WEST,J_SOUTH)          &  !<-- Parent's PD interp'd to child task update points
                            +WGHT_SE_H(I,J)*PD(I_EAST,J_SOUTH)          &
                            +WGHT_NW_H(I,J)*PD(I_WEST,J_NORTH)          &
                            +WGHT_NE_H(I,J)*PD(I_EAST,J_NORTH)
!
            ENDDO
          ENDDO
!
          DO L=NUM_LYRS,1,-1
!
            PDTOP_PT=SG1(L+1)*PDTOP+PT
!
            DO J=J1,J2                                                     !<-- J limits of child task update region on parent task
              J_SOUTH=J_PARENT_SOUTH_H(J)
              J_NORTH=J_PARENT_NORTH_H(J)
!
              DO I=I1,I2                                                   !<-- I limits of child task update region on parent task
                I_WEST=I_PARENT_WEST_H(I)
                I_EAST=I_PARENT_EAST_H(I)
!
                PX_SW=SG2(L)*PD(I_WEST,J_SOUTH)+PDTOP_PT                   !<-- Pressure, top of layer L, parent point SW of nest point
                PX_SE=SG2(L)*PD(I_EAST,J_SOUTH)+PDTOP_PT                   !<-- Pressure, top of layer L, parent point SE of nest point
                PX_NW=SG2(L)*PD(I_WEST,J_NORTH)+PDTOP_PT                   !<-- Pressure, top of layer L, parent point NW of nest point
                PX_NE=SG2(L)*PD(I_EAST,J_NORTH)+PDTOP_PT                   !<-- Pressure, top of layer L, parent point NE of nest point
!
                PINT_INTERP(I,J,L)=WGHT_SW_H(I,J)*PX_SW                 &  !<-- Top interface hydrostatic pressure interpolated to
                                  +WGHT_SE_H(I,J)*PX_SE                 &  !    update point for child task N.  These are the source
                                  +WGHT_NW_H(I,J)*PX_NW                 &  !    pressures for interface variables.
                                  +WGHT_NE_H(I,J)*PX_NE                    !<--
!
                PMID_INTERP(I,J,L)=0.5*(PINT_INTERP(I,J,L)              &  !<-- Parent's midlayer hydrostatic pressure at nest update
                                       +PINT_INTERP(I,J,L+1))              !    points.  Source pressures for midlayer variables. 
!
                T_INTERP=WGHT_SW_H(I,J)*T(I_WEST,J_SOUTH,L)             &  !<-- T interp'd to update point for child task N
                        +WGHT_SE_H(I,J)*T(I_EAST,J_SOUTH,L)             &  !
                        +WGHT_NW_H(I,J)*T(I_WEST,J_NORTH,L)             &  !
                        +WGHT_NE_H(I,J)*T(I_EAST,J_NORTH,L)                !<--
!
                Q_INTERP=WGHT_SW_H(I,J)*Q(I_WEST,J_SOUTH,L)             &  !<-- Q interp'd to update point for child task N
                        +WGHT_SE_H(I,J)*Q(I_EAST,J_SOUTH,L)             &  !
                        +WGHT_NW_H(I,J)*Q(I_WEST,J_NORTH,L)             &  !
                        +WGHT_NE_H(I,J)*Q(I_EAST,J_NORTH,L)                !<--
!
                CW_INTERP=WGHT_SW_H(I,J)*CW(I_WEST,J_SOUTH,L)           &  !<-- CW interp'd to update point for child task N
                         +WGHT_SE_H(I,J)*CW(I_EAST,J_SOUTH,L)           &  !
                         +WGHT_NW_H(I,J)*CW(I_WEST,J_NORTH,L)           &  !
                         +WGHT_NE_H(I,J)*CW(I_EAST,J_NORTH,L)              !<--
!
                DP=DSG2(L)*PD_INTERP(I,J)+PDSG1(L)
!
                TMP=R_D*T_INTERP*((1.-CW_INTERP)+P608*Q_INTERP)
                LOG_PTOP(I,J)=LOG(PINT_INTERP(I,J,L))
!
                PHI_INTERP(I,J,L)=PHI_INTERP(I,J,L+1)                   &  !<-- Top interface geopotl of parent at child update point I,J
                                 +TMP*(LOG_PBOT(I,J)-LOG_PTOP(I,J))
!
                LOG_PBOT(I,J)=LOG_PTOP(I,J)
!
              ENDDO
            ENDDO
!
          ENDDO
!
!-----------------------------------------------------------------------
!***  Use the sfc geopotential at the nest points to derive the 
!***  value of PD at the nest points based on the parent's heights
!***  and pressures on the parent's layer interfaces over the
!***  child's points.
!
!***  If the child's terrain is lower than the value of the parent's
!***  terrain interpolated to the child point then extrapolate the
!***  parent's interpolated sfc pressure down to the child's terrain
!***  quadratically.
!-----------------------------------------------------------------------
!
          DO J=J1,J2
            J_TRANS=J+J_OFFSET                                             !<-- J on full nest resolution of parent at given nest J
!
            DO I=I1,I2
              I_TRANS=I+I_OFFSET                                           !<-- I on full nest resolution of parent at given nest I
!
              IF(FIS_CHILD(I_TRANS,J_TRANS)<PHI_INTERP(I,J,LM+1))THEN      !<-- If so, child's terrain lies below parent's 
!
                COEFF_1=(PINT_INTERP(I,J,LM+1)-PINT_INTERP(I,J,LM))     &  !<-- Coefficient for linear term
                       /(PHI_INTERP(I,J,LM+1)-PHI_INTERP(I,J,LM))
!
                COEFF_2=(COEFF_1                                        &  !<-- Coefficient for quadratic term
                        -(PINT_INTERP(I,J,LM)-PINT_INTERP(I,J,LM-1))    &  !
                        /(PHI_INTERP(I,J,LM)-PHI_INTERP(I,J,LM-1)))     &  !
                        /(PHI_INTERP(I,J,LM+1)-PHI_INTERP(I,J,LM-1))       !
!
                PHI_DIFF=FIS_CHILD(I_TRANS,J_TRANS)-PHI_INTERP(I,J,LM+1)   !<-- Diff between nest sfc geopotential and parent's interp'd
!
                PSFC_CHILD=COEFF_2*PHI_DIFF*PHI_DIFF                    &  !<-- Parent pressure at child's surface elevation
                          +COEFF_1*PHI_DIFF                             &
                          +PINT_INTERP(I,J,LM+1)
!
              ELSE                                                         !<-- Child's terrain is at or above parent's
                DO L=NUM_LYRS,1,-1
                  IF(FIS_CHILD(I_TRANS,J_TRANS)<PHI_INTERP(I,J,L))THEN
                    LOG_P1_PARENT=LOG(PINT_INTERP(I,J,L+1))                !<-- Log pressure on bottom of parent model layer
                    D_LNP_DFI=(LOG_P1_PARENT-LOG(PINT_INTERP(I,J,L)))   &  !<-- d[ln(p)]/d[fi] in parent layer L
                             /(PHI_INTERP(I,J,L+1)-PHI_INTERP(I,J,L))
                    PSFC_CHILD=EXP(LOG_P1_PARENT                        &  !<-- Parent pressure at child's surface elevation
                                  +D_LNP_DFI                            &
                                   *(FIS_CHILD(I_TRANS,J_TRANS)         &
                                    -PHI_INTERP(I,J,L+1)))
                    EXIT
                  ENDIF
                ENDDO
              ENDIF
!
              PD_CHILD(I,J)=PSFC_CHILD-PT                                  !<-- Parent's approx. of child PD on child's update point
!
              PINT_CHILD(I,J,NUM_LYRS+1)=PSFC_CHILD                        !<-- Parent's approx. of child pressure at child's sfc
!
            ENDDO
          ENDDO
!
!-----------------------------------------------------------------------
!***  The nest's interface and midlayer hydrostatic pressures
!***  at the update locations.  These are the target pressures
!***  for the interface and midlayer variables, respectively.
!-----------------------------------------------------------------------
!
          DO L=NUM_LYRS,1,-1
            PDTOP_PT=SG1(L+1)*PDTOP+PT
!
            DO J=J1,J2
            DO I=I1,I2
              PINT_CHILD(I,J,L)=SG2(L)*PD_CHILD(I,J)+PDTOP_PT
              PMID_CHILD(I,J,L)=0.5*(PINT_CHILD(I,J,L)                  &
                                    +PINT_CHILD(I,J,L+1))
            ENDDO
            ENDDO
!
          ENDDO
!
!-----------------------------------------------------------------------
!***  The parent extracts each update variable from the Move Bundle
!***  that will be sent to those nest points that have moved beyond
!***  the footprint of the nest domain's previous location and onto
!***  a new region of the parent domain then updates these variables 
!***  by: (1) selecting the nearest neighbor (integers); (2) simple
!***  bilinear interpolation; (3) bilinear interpolation that accounts
!***  for the presence of the sea mask; (4) bilinear interpolation that
!***  accounts for the presence of the land mask.  A fifth method of
!***  updates in these regions is done by reading the values directly
!***  from external files but that is done by the moving nest tasks
!***  themselves in UPDATE_INTERIOR_FROM_PARENT.
!-----------------------------------------------------------------------
!
          field_loop: DO N_FIELD=1,NUM_FIELDS_MOVE
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Extract N_Fieldth Field From Move_Bundle"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_FieldBundleGet(FIELDBUNDLE=MOVE_BUNDLE            &  !<-- Bundle holding the arrays for move updates
                                    ,fieldIndex =N_FIELD                &  !<-- Index of the Field in the Bundle
                                    ,field      =HOLD_FIELD             &  !<-- Field N_FIELD in the Bundle
                                    ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_UPDATE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Type and Dimensions of the Field"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_FieldGet(field   =HOLD_FIELD                      &  !<-- Field N_FIELD in the Bundle
                              ,dimCount=NUM_DIMS                        &  !<-- Is this Field 2-D or 3-D?
                              ,typeKind=DATATYPE                        &  !<-- Is this Field integer or real?
                              ,name    =FIELD_NAME                      &  !<-- This Field's name
                              ,rc      =RC )
!
            N_REMOVE=INDEX(FIELD_NAME,SUFFIX_MOVE)
            FIELD_NAME=FIELD_NAME(1:N_REMOVE-1)                            !<-- Remove Move Bundle Fieldname's suffix '-move'
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_UPDATE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Extract Type of Update by the Parent"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_AttributeGet(field=HOLD_FIELD                     &  !<-- Take Attribute from this Field
                                  ,name ='UPDATE_TYPE'                  &  !<-- The Attribute's name
                                  ,value=UPDATE_TYPE_INT                &  !<-- The Attribute's value
                                  ,rc   =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_UPDATE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            IF(UPDATE_TYPE_INT==1)THEN
              UPDATE_TYPE_CHAR='H'                                         !<-- Ordinary H-pt variable 
            ELSEIF(UPDATE_TYPE_INT==2)THEN
              UPDATE_TYPE_CHAR='L'                                         !<-- H-pt land sfc variable 
            ELSEIF(UPDATE_TYPE_INT==3)THEN
              UPDATE_TYPE_CHAR='W'                                         !<-- H-pt water sfc variable 
            ELSEIF(UPDATE_TYPE_INT==4)THEN
              UPDATE_TYPE_CHAR='F'                                         !<-- H-pt variable updated from external file
            ELSEIF(UPDATE_TYPE_INT==5)THEN
              UPDATE_TYPE_CHAR='V'                                         !<-- Ordinary V-pt variable
            ENDIF
!
!-----------------------------------------------------------------------
!***  Variables that nests will read from external files are 
!***  simply skipped here.
!-----------------------------------------------------------------------
!
            IF(UPDATE_TYPE_CHAR=='F')THEN
!
              CYCLE field_loop
!
            ENDIF
!
!-----------------------------------------------------------------------
!***  Interpolate the parent's update variables to the nest update
!***  locations.  Of course this is simple bilinear interpolation
!***  for many 2-D variables.  Since we now have the parent's and the
!***  nest's vertical pressure distributions at each of the nest's
!***  horizontal update locations the 3-D update variables can be
!***  interpolated to the nest's locations using bilinear in the
!***  horizontal and cubic spline in the vertical.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  2-D Fields 
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
            dims_2_or_3: IF(NUM_DIMS==2)THEN
!
!--------------------------------
!***  Integer - Nearest Neighbor
!--------------------------------
!
              type_2d : IF(DATATYPE==ESMF_TYPEKIND_I4)THEN                 !<-- Field N holds an integer array
!
                CALL ESMF_FieldGet(field    =HOLD_FIELD                 &  !<-- Field N_FIELD in the Bundle
                                  ,localDe  =0                          &
                                  ,farrayPtr=IARRAY_2D                  &  !<-- Dummy 2-D integer array with Field's data
                                  ,rc       =RC )
!
                DO J=J1,J2
                  J_SOUTH=J_PARENT_SOUTH(J)
                  J_NORTH=J_PARENT_NORTH(J)
!
                  DO I=I1,I2
                    I_WEST =I_PARENT_WEST(I)
                    I_EAST =I_PARENT_EAST(I)

                    KNT_INTEGER_PTS(N)=KNT_INTEGER_PTS(N)+1                !<-- Total integer points updated in composite output.
!
                    MAX_WGHT=MIN(WGHT_SW(I,J),WGHT_SE(I,J)              &  !<-- What is the greatest weight of parent points?
                                ,WGHT_NW(I,J),WGHT_NE(I,J))
!
                    IF(MAX_WGHT==WGHT_SW(I,J))THEN                         !<-- Assign integer parent value with greatest weight.
                      IVAL=IARRAY_2D(I_WEST,J_SOUTH)                       !
                    ELSEIF(MAX_WGHT==WGHT_SE(I,J))THEN                     !
                      IVAL=IARRAY_2D(I_EAST,J_SOUTH)                       !
                    ELSEIF(MAX_WGHT==WGHT_NW(I,J))THEN                     !
                      IVAL=IARRAY_2D(I_WEST,J_NORTH)                       !
                    ELSEIF(MAX_WGHT==WGHT_NE(I,J))THEN                     !
                      IVAL=IARRAY_2D(I_EAST,J_NORTH)                       !
                    ENDIF                                                  !<--
!
                    CHILD_UPDATE_DATA%TASKS(N)%DATA_INTEGER(KNT_INTEGER_PTS(N))=IVAL
                  ENDDO
                ENDDO
!
!---------------------
!***  Real - Bilinear
!---------------------
!
              ELSEIF(DATATYPE==ESMF_TYPEKIND_R4)THEN  type_2d              !<-- Field N holds a real array
!
                CALL ESMF_FieldGet(field    =HOLD_FIELD                 &  !<-- Field N_FIELD in the Bundle
                                  ,localDe  =0                          &
                                  ,farrayPtr=ARRAY_2D                   &  !<-- Dummy 2-D real array with Field's data
                                  ,rc       =RC )
!
!-----------------------------------------------------------------------
!***  Latitude and longitude and variables directly related to them
!***  will be computed by the moving nests in the parent update regions
!***  so the parent does not update those variables.  However they must
!***  be designated as shift variables in order for the intra- and
!***  inter-task shifts to be executed by the nest therefore we must
!***  insert dummy values and increment the total word count of the
!***  parent update data object.
!-----------------------------------------------------------------------
!
                IF(FIELD_NAME=='GLAT'.OR.FIELD_NAME=='GLON'             &
                         .OR.                                           &
                   FIELD_NAME=='VLAT'.OR.FIELD_NAME=='VLON'             &
                         .OR.                                           &
                   FIELD_NAME=='HDACX'.OR.FIELD_NAME=='HDACY'           &
                         .OR.                                           &
                   FIELD_NAME=='HDACVX'.OR.FIELD_NAME=='HDACVY'         &
                         .OR.                                           &
                   FIELD_NAME=='F') THEN
!
                  KNT_DUMMY=0
!
                  DO J=J1,J2
                  DO I=I1,I2
                    KNT_DUMMY=KNT_DUMMY+1
                    CHILD_UPDATE_DATA%TASKS(N)%DATA_REAL(KNT_REAL_PTS(N)+KNT_DUMMY)=  &  !<-- Dummy values needed for correct word count
                                                                               -999.0
                  ENDDO
                  ENDDO
!
                  KNT_REAL_PTS(N)=KNT_REAL_PTS(N)+(I2-I1+1)*(J2-J1+1)        !<-- # of dummy values added to composite output.
!
                  CYCLE field_loop
!
                ENDIF
!
!-----------------------------------------------------------------------
!***  The parent has already computed its approximation of the
!***  child's PD so just insert it into the composite data object.
!***  NOTE:  For the time being the parent will use its generated
!***         values of PD on the nest in filling the nest's PDO.
!-----------------------------------------------------------------------
!
!--------
!***  PD
!--------
!
                IF(FIELD_NAME=='PD'.OR.FIELD_NAME=='PDO')THEN
!
                  DO J=J1,J2
                    J_SOUTH=J_PARENT_SOUTH(J)
                    J_NORTH=J_PARENT_NORTH(J)
!
                    DO I=I1,I2
                      I_WEST=I_PARENT_WEST(I)
                      I_EAST=I_PARENT_EAST(I)
!
                      KNT_REAL_PTS(N)=KNT_REAL_PTS(N)+1                          !<-- Total real points updated in composite output.
!
                      CHILD_UPDATE_DATA%TASKS(N)%DATA_REAL(KNT_REAL_PTS(N))=  &  !<-- Save the parent's computation of the nest's
                                                                 PD_CHILD(I,J)   !    new PD into the data object.
!
                    ENDDO
                  ENDDO
!
                  CYCLE field_loop
!
                ENDIF
!
!-----------------------------------------------------------------------
!***  Nest-resolution values of FIS and SM were read by the parent
!***  tasks and sent into this routine since they are required by
!***  the parent for generating update values for its moving children.
!***  Insert them directly into the data object of update values.
!-----------------------------------------------------------------------
!
!!!             ELSEIF(TRIM(FIELD_NAME)=='FIS')THEN
!
!!!               DO J=J1,J2
!!!                 J_SOUTH=J_PARENT_SOUTH(J)
!!!                 J_NORTH=J_PARENT_NORTH(J)
!!!                 J_TRANS=J+(J_PARENT_SW-JDS_CHILD)*PARENT_CHILD_SPACE_RATIO   !<-- J on full nest resolution of parent at given nest J
!
!!!                 DO I=I1,I2
!!!                   I_WEST=I_PARENT_WEST(I)
!!!                   I_EAST=I_PARENT_EAST(I)
!!!                   I_TRANS=I+(I_PARENT_SW-IDS_CHILD)*PARENT_CHILD_SPACE_RATIO !<-- I on full nest resolution of parent at given nest I
!
!!!                   KNT_REAL_PTS(N)=KNT_REAL_PTS(N)+1                          !<-- Total real points updated in composite output.
!
!!!                   CHILD_UPDATE_DATA%TASKS(N)%DATA_REAL(KNT_REAL_PTS(N))=  &  !<-- Insert the nest-resolution values of FIS
!!!                                               FIS_CHILD(I_TRANS,J_TRANS)     !    that the parent had read in earlier.
!
!!!                 ENDDO
!!!               ENDDO
!
!               ELSEIF(TRIM(FIELD_NAME)=='SM')THEN
!
!                 DO J=J1,J2
!                   J_SOUTH=J_PARENT_SOUTH(J)
!                   J_NORTH=J_PARENT_NORTH(J)
!
!                   DO I=I1,I2
!                     I_WEST=I_PARENT_WEST(I)
!                     I_EAST=I_PARENT_EAST(I)
!
!                     KNT_REAL_PTS(N)=KNT_REAL_PTS(N)+1                          !<-- Total real points updated in composite output.
!
!                     CHILD_UPDATE_DATA%TASKS(N)%DATA_REAL(KNT_REAL_PTS(N))=  &  !<-- Insert the nest-resolution values of SM
!                                                                SM_CHILD(I,J)   !    that the parent had read in earlier.
!
!                   ENDDO
!                 ENDDO
!
!----------------------------------
!***  All other 2-D Real variables
!----------------------------------
!
!-----------------------------------------------------------------------
!***  If the 2-D Real variable is not related to the land/sea surface
!***  then bilinearly interpolate it to the nest points' locations
!***  and insert it into the data object.
!-----------------------------------------------------------------------
!
                IF(UPDATE_TYPE_CHAR/='L'.AND.UPDATE_TYPE_CHAR/='W')THEN
!
!-----------------------------------------------------------------------
!
                  DO J=J1,J2
                    J_SOUTH=J_PARENT_SOUTH(J)
                    J_NORTH=J_PARENT_NORTH(J)
!
                    DO I=I1,I2
                      I_WEST=I_PARENT_WEST(I)
                      I_EAST=I_PARENT_EAST(I)
!
                      KNT_REAL_PTS(N)=KNT_REAL_PTS(N)+1                          !<-- Total real points updated in composite output.
!
                      CHILD_UPDATE_DATA%TASKS(N)%DATA_REAL(KNT_REAL_PTS(N))=  &  !<-- Parent's 2-D real variable interpolated 
                                       WGHT_SW(I,J)*ARRAY_2D(I_WEST,J_SOUTH)  &  !    horizontally to the moving nest's
                                      +WGHT_SE(I,J)*ARRAY_2D(I_EAST,J_SOUTH)  &  !    update location on nest task N's
                                      +WGHT_NW(I,J)*ARRAY_2D(I_WEST,J_NORTH)  &  !    subdomain.
                                      +WGHT_NE(I,J)*ARRAY_2D(I_EAST,J_NORTH)     !<--
                    ENDDO
                  ENDDO
!
                  CYCLE field_loop
!
                ENDIF
!
!-----------------------------------------------------------------------
!***  If the 2-D Real variable is related to the land/sea surface
!***  then the bilinear interpolation must also take into account
!***  the exclusion of either sea or land points in the final
!***  4-pt summation.
!-----------------------------------------------------------------------
!
                IF(UPDATE_TYPE_CHAR=='L')THEN                              !<-- If true, a 2-D land surface real variable
!
!-----------------------------------------------------------------------
!
                  DO J=J1,J2
                    J_SOUTH=J_PARENT_SOUTH(J)
                    J_NORTH=J_PARENT_NORTH(J)
!
                    DO I=I1,I2
                      I_WEST=I_PARENT_WEST(I)
                      I_EAST=I_PARENT_EAST(I)
!
                      KNT_REAL_PTS(N)=KNT_REAL_PTS(N)+1                          !<-- Total real points updated in composite output.
!
                      CHILD_UPDATE_DATA%TASKS(N)%DATA_REAL(KNT_REAL_PTS(N))=  &  !<-- Parent's 2-D real land variable interpolated 
                                 PROD_LWGT_SW(I,J)*ARRAY_2D(I_WEST,J_SOUTH)   &  !    horizontally to the moving nest's
                                +PROD_LWGT_SE(I,J)*ARRAY_2D(I_EAST,J_SOUTH)   &  !    update location using the land mask.
                                +PROD_LWGT_NW(I,J)*ARRAY_2D(I_WEST,J_NORTH)   &  !
                                +PROD_LWGT_NE(I,J)*ARRAY_2D(I_EAST,J_NORTH)      !<--
                    ENDDO
                  ENDDO
!
!-----------------------------------------------------------------------
!
                ELSEIF(UPDATE_TYPE_CHAR=='W')THEN                          !<-- If true, a 2-D water surface real variable
!
!-----------------------------------------------------------------------
!
                  DO J=J1,J2
                    J_SOUTH=J_PARENT_SOUTH(J)
                    J_NORTH=J_PARENT_NORTH(J)
!
                    DO I=I1,I2
                      I_WEST=I_PARENT_WEST(I)
                      I_EAST=I_PARENT_EAST(I)
!
                      KNT_REAL_PTS(N)=KNT_REAL_PTS(N)+1                          !<-- Total real points updated in composite output.
!
                      CHILD_UPDATE_DATA%TASKS(N)%DATA_REAL(KNT_REAL_PTS(N))=  &  !<-- Parent's 2-D real sea variable interpolated 
                                 PROD_SWGT_SW(I,J)*ARRAY_2D(I_WEST,J_SOUTH)   &  !    horizontally to the moving nest's
                                +PROD_SWGT_SE(I,J)*ARRAY_2D(I_EAST,J_SOUTH)   &  !    update location using the land mask.
                                +PROD_SWGT_NW(I,J)*ARRAY_2D(I_WEST,J_NORTH)   &  !
                                +PROD_SWGT_NE(I,J)*ARRAY_2D(I_EAST,J_NORTH)      !<--
                    ENDDO
                  ENDDO
!
!-----------------------------------------------------------------------
!
                ENDIF
!
!-----------------------------------------------------------------------
!
              ENDIF  type_2d
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  3-D Fields --> Bilinear interpolation
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
            ELSEIF(NUM_DIMS==3)THEN  dims_2_or_3
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
              KLO=LIMITS_LO(3)
              KHI=LIMITS_HI(3)
              NUM_LEVELS=KHI-KLO+1                                         !<-- # of levels in this 3-D Real variable
!
!-----------------------------------------------------------------------
!***  The nature of the unique Q2 array complicates how the parent
!***  must interpolate its values to the nest gridpoints.  Q2 is a
!***  3-D array that lie on the layer interfaces BUT its level K=1
!***  is the BOTTOM of the uppermost model layer and not the top of
!***  the uppermost layer.  Thus while there are NUM_LYRS+1 layer
!***  interfaces there are only NUM_LYRS levels in Q2 that correspond
!***  with interfaces 2->NUM_LYRS+1.  Rather than insert an assortment
!***  of confusing IF tests to make a single set of code be generic,
!***  separate Q2 from the rest of the variables and deal with it
!***  inside its own block.
!-----------------------------------------------------------------------
!
              q2: IF(FIELD_NAME=='Q2')THEN
!
!-----------------------------------------------------------------------
!***  We must add a new level to the top of the Q2 or E2 data in case
!***  their pressures at the bottom of the uppermost layer are less
!***  than on the bottom of the uppermost layer in the parent.
!-----------------------------------------------------------------------
!
                IF(ALLOCATED(P_OUTPUT))DEALLOCATE(P_OUTPUT)
                ALLOCATE(P_OUTPUT(KLO:KHI+1))
!
                IF(ALLOCATED(VBL_INTERP))DEALLOCATE(VBL_INTERP)
                ALLOCATE(VBL_INTERP(I1:I2,J1:J2,KLO:KHI+1+1))
                ALLOCATE(VBL_COL_X(KLO:KHI+1)) 
!
                DO L=KLO,KHI+2
                  P_INPUT(L)=0.         
                  VBL_INPUT(L)=0.      
                ENDDO
!
                DO J=J1,J2                        
                DO I=I1,I2                       
                  VBL_INTERP(I,J,KHI+1+1)=0.    
                ENDDO
                ENDDO                          
!
                DO L=KLO,KHI
                  DO J=J1,J2
                    J_SOUTH=J_PARENT_SOUTH(J)
                    J_NORTH=J_PARENT_NORTH(J)
!
                    DO I=I1,I2
                      I_WEST=I_PARENT_WEST(I)
                      I_EAST=I_PARENT_EAST(I)
!
                      VBL_INTERP(I,J,L+1)=                              &  !<-- Parent's 3-D variable interpolated
                               WGHT_SW(I,J)*ARRAY_3D(I_WEST,J_SOUTH,L)  &  !    horizontally to the moving nest's
                              +WGHT_SE(I,J)*ARRAY_3D(I_EAST,J_SOUTH,L)  &  !    update location.
                              +WGHT_NW(I,J)*ARRAY_3D(I_WEST,J_NORTH,L)  &  !
                              +WGHT_NE(I,J)*ARRAY_3D(I_EAST,J_NORTH,L)     !<--
                    ENDDO
                  ENDDO
                ENDDO
!
                DO J=J1,J2
                DO I=I1,I2
                  VBL_INTERP(I,J,1)=VBL_INTERP(I,J,2)                      !<-- Fill in the artificial top level.
                ENDDO
                ENDDO
!
!-----------------------------------------------------------------------
!***  Use cubic spline interpolation to move variables to child update
!***  point levels from their original vertical locations in the column
!***  following horizontal interpolation from the surrounding parent 
!***  points.  The target locations are the new interface pressures 
!***  in the nest update point columns based on the new surface
!***  pressure for the nest's terrain. 
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  If the target location lies below the lowest parent pressure 
!***  level in the newly created child column then extrapolate linearly
!***  in pressure to obtain a value at the lowest child level and
!***  fill in the remaining 'underground' levels using the call to
!***  'SPLINE' just as with all the other levels above it.
!-----------------------------------------------------------------------
!
                N_STRIDE=NPOINTS_HORIZ
                N_ADD   =NPOINTS_HORIZ*(NUM_LEVELS-1)
!
                NUM_LEVS_SEC=NUM_LEVELS+1+1                                !<-- Use this many levels in the 2nd derivative array
                ALLOCATE(SEC_DERIV(1:NUM_LEVS_SEC))                        !<-- Allocate 1 longer in case we increase the
!                                                                          !    # of input levels below.
                LOC_1=KNT_REAL_PTS(N)
!
                P3D_INPUT=>PINT_INTERP
                P3D_OUTPUT=>PINT_CHILD
!
                DO J=J1,J2
                DO I=I1,I2
!
                  DO L=1,NUM_LEVELS+1                                      !<-- We are adding a temporary top level to Q2 and E2
                    P_INPUT  (L)=P3D_INPUT(I,J,L)                          !<-- Parent input pressures over nest update point
                    P_OUTPUT (L)=P3D_OUTPUT(I,J,L)                         !<-- Nest target pressures over nest update point
                    VBL_INPUT(L)=VBL_INTERP(I,J,L)                         !<-- Values of parent variable values over nest update point
                  ENDDO
!
                  NUM_LEVS_IN=NUM_LEVELS+1
                  LOC_1=LOC_1+1
                  LOC_2=LOC_1+N_ADD
                  VBL_COL_CHILD=>CHILD_UPDATE_DATA%TASKS(N)%DATA_REAL(LOC_1:LOC_2:N_STRIDE) !<-- Point working column pointer into
                                                                                            !    the 1-D rendering of this 3-D real 
                                                                                            !    update variable in the composite output.
!
                  IF(P_OUTPUT(NUM_LEVELS+1)>P_INPUT(NUM_LEVELS+1))THEN     !<-- The nest's bottom level is below the parent's
                    NUM_LEVS_IN=NUM_LEVELS+1+1                             !    so add another input level that is the same
                    P_INPUT(NUM_LEVELS+1+1)=P_OUTPUT(NUM_LEVELS+1)         !    as the nest's lowest level.
                    R_DELP=1./(P_INPUT(NUM_LEVELS+1)-P_INPUT(NUM_LEVELS))
                    DELP_EXTRAP=P_OUTPUT(NUM_LEVELS+1)                  &
                               -P_INPUT(NUM_LEVELS+1)
!
                    COEFF_1=(VBL_INPUT(NUM_LEVELS+1)                    &
                            -VBL_INPUT(NUM_LEVELS))*R_DELP
                    FACTOR=HYPER_A/(DELP_EXTRAP+HYPER_A)
                    VBL_INPUT(NUM_LEVELS+1+1)=VBL_INPUT(NUM_LEVELS+1)   &  !<-- Create extrapolated value at parent's new lowest
                                             +COEFF_1*DELP_EXTRAP*FACTOR   !    level for input to the spline.
                  ENDIF
!
                  DO L=1,NUM_LEVS_SEC  
                    SEC_DERIV(L)=0.                                        !<-- Initialize 2nd derivatives of the spline to zero.
                  ENDDO
!
                  CALL SPLINE(NUM_LEVS_IN                               &  !<-- # of input levels
                             ,P_INPUT                                   &  !<-- Input variable is at these input pressure values
                             ,VBL_INPUT                                 &  !<-- The column of input variable values
                             ,SEC_DERIV                                 &  !<-- Specified 2nd derivatives (=0) at parent levels
                             ,NUM_LEVS_SEC                              &  !<-- Vertical dimension of SEC_DERIV
                             ,NUM_LEVELS+1                              &  !<-- # of child target levels to interpolate to
                             ,P_OUTPUT                                  &  !<-- Child target pressure values to interpolate to
                             ,VBL_COL_X)                                   !<-- Child values of variable returned on P_OUTPUT levels
!
                  DO L=KLO,KHI
                    VBL_COL_CHILD(L)=VBL_COL_X(L+1)                        !<-- Eliminate the artificial level on top of layer 1.
                  ENDDO
!
                ENDDO
                ENDDO
!
                KNT_REAL_PTS(N)=KNT_REAL_PTS(N)+NPOINTS_HORIZ*NUM_LEVELS   !<-- Total points updated in composite output after 
!                                                                          !    this 3-D real variable was done.
                DEALLOCATE(VBL_COL_X)
                DEALLOCATE(SEC_DERIV)
!
!-----------------------------------------------------------------------
!
              ELSE q2                                                      !<-- All 3-D variables that are not Q2
!
!-----------------------------------------------------------------------
!
                MIDLAYERS=.FALSE.
                INTERFACES=.FALSE.
!
                IF(NUM_LEVELS==NUM_LYRS)THEN
                  MIDLAYERS=.TRUE.
                ELSEIF(NUM_LEVELS==NUM_LYRS+1)THEN
                  INTERFACES=.TRUE.
                ENDIF
!
                IF(ALLOCATED(P_OUTPUT))DEALLOCATE(P_OUTPUT)
                ALLOCATE(P_OUTPUT(KLO:KHI))
!
                IF(ALLOCATED(VBL_INTERP))DEALLOCATE(VBL_INTERP)
                ALLOCATE(VBL_INTERP(I1:I2,J1:J2,KLO:KHI+1))
!
!-----------------------------------------------------------------------
!***  Use cubic spline interpolation to move variables to child update
!***  point levels from their original vertical locations in the column
!***  following horizontal interpolation from the surrounding parent 
!***  points.  The target locations are the new pressures values
!***  in the nest update point columns based on the new surface
!***  pressure for the nest's terrain.  However this is obviously
!***  done only for atmospheric variables.  Of course it is not
!***  done for 3-D land surface variables.
!-----------------------------------------------------------------------
!
                soil_or_not: IF(UPDATE_TYPE_CHAR/='L')THEN                 !<-- 3-D H-pt variable that is not soil.
!
!-----------------------------------------------------------------------
!
                  DO L=1,NUM_LYRS+2                                        !<-- Maximum # of levels to be used.
                    P_INPUT(L)=0.  
                    VBL_INPUT(L)=0.
                  ENDDO
!
                  DO J=J1,J2                           
                  DO I=I1,I2                          
                    VBL_INTERP(I,J,KHI+1)=0.         
                  ENDDO
                  ENDDO                             
!
                  DO L=KLO,KHI
                    DO J=J1,J2
                      J_SOUTH=J_PARENT_SOUTH(J)
                      J_NORTH=J_PARENT_NORTH(J)
!
                      DO I=I1,I2
                        I_WEST=I_PARENT_WEST(I)
                        I_EAST=I_PARENT_EAST(I)
!
                        VBL_INTERP(I,J,L)=                              &  !<-- Parent's 3-D variable interpolated
                               WGHT_SW(I,J)*ARRAY_3D(I_WEST,J_SOUTH,L)  &  !    horizontally to the moving nest's
                              +WGHT_SE(I,J)*ARRAY_3D(I_EAST,J_SOUTH,L)  &  !    update location.
                              +WGHT_NW(I,J)*ARRAY_3D(I_WEST,J_NORTH,L)  &  !
                              +WGHT_NE(I,J)*ARRAY_3D(I_EAST,J_NORTH,L)     !<--
                      ENDDO
                    ENDDO
                  ENDDO
!
!-----------------------------------------------------------------------
!***  If the target location lies below the lowest parent pressure
!***  level in the newly created child column then extrapolate linearly
!***  in pressure to obtain a value at the lowest child level and
!***  fill in the remaining 'underground' levels using the call to
!***  'SPLINE' just as with all the other levels above it.
!-----------------------------------------------------------------------
!
                  N_STRIDE=NPOINTS_HORIZ
                  N_ADD   =NPOINTS_HORIZ*(NUM_LEVELS-1)
!
                  NUM_LEVS_SEC=NUM_LEVELS+1                                !<-- Use this many levels in the 2nd derivative array
                  ALLOCATE(SEC_DERIV(1:NUM_LEVS_SEC))                      !<-- Allocate 1 longer in case we increase the
!                                                                          !    # of input levels below.
                  LOC_1=KNT_REAL_PTS(N)
!
                  IF(MIDLAYERS)THEN                                        !<-- Input/output pressures are at midlayers
                    P3D_INPUT=>PMID_INTERP
                    P3D_OUTPUT=>PMID_CHILD
!
                  ELSEIF(INTERFACES)THEN                                   !<-- Input/output pressures are at interfaces
                    P3D_INPUT=>PINT_INTERP
                    P3D_OUTPUT=>PINT_CHILD
!
                  ELSE
                    WRITE(0,*)' # of levels in 3-D variable is ',NUM_LEVELS
                    WRITE(0,*)' That is not midlayer, interface, or soil.'
                    WRITE(0,*)' ABORT!!'
                    CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
!
                  ENDIF
!
                  DO J=J1,J2
                  DO I=I1,I2
!
                    DO L=1,NUM_LEVELS                                      !<-- Variable has NUM_LEVELS levels in parent and nest
                      P_INPUT  (L)=P3D_INPUT(I,J,L)                        !<-- Parent input pressures over nest update point
                      P_OUTPUT (L)=P3D_OUTPUT(I,J,L)                       !<-- Nest target pressures over nest update point
                      VBL_INPUT(L)=VBL_INTERP(I,J,L)                       !<-- Values of parent variable values over nest update point
                    ENDDO
!
                    NUM_LEVS_IN=NUM_LEVELS
                    LOC_1=LOC_1+1
                    LOC_2=LOC_1+N_ADD
                    VBL_COL_CHILD=>CHILD_UPDATE_DATA%TASKS(N)%DATA_REAL(LOC_1:LOC_2:N_STRIDE) !<-- Point working column pointer into
                                                                                              !    the 1-D rendering of this 3-D real 
                                                                                              !    update variable in the composite output.
!
                    IF(P_OUTPUT(NUM_LEVELS)>P_INPUT(NUM_LEVELS))THEN       !<-- The nest's bottom level is below the parent's
                      NUM_LEVS_IN=NUM_LEVELS+1                             !    so add another input level that is the same
                      P_INPUT(NUM_LEVELS+1)=P_OUTPUT(NUM_LEVELS)           !    as the nest's lowest level.
                      R_DELP=1./(P_INPUT(NUM_LEVELS)-P_INPUT(NUM_LEVELS-1))
                      DELP_EXTRAP=P_OUTPUT(NUM_LEVELS)                  &
                                 -P_INPUT(NUM_LEVELS)
!
                      COEFF_1=(VBL_INPUT(NUM_LEVELS)                    &
                              -VBL_INPUT(NUM_LEVELS-1))*R_DELP
                      FACTOR=HYPER_A/(DELP_EXTRAP+HYPER_A)
                      VBL_INPUT(NUM_LEVELS+1)=VBL_INPUT(NUM_LEVELS)     &  !<-- Create extrapolated value at parent's new lowest
                                             +COEFF_1*DELP_EXTRAP*FACTOR   !    level for input to the spline.
                    ENDIF
!
                    DO L=1,NUM_LEVS_SEC  
                      SEC_DERIV(L)=0.                                      !<-- Initialize 2nd derivatives of the spline to zero.
                    ENDDO
!
                    CALL SPLINE(NUM_LEVS_IN                             &  !<-- # of input levels
                               ,P_INPUT                                 &  !<-- Input variable is at these input pressure values
                               ,VBL_INPUT                               &  !<-- The column of input variable values
                               ,SEC_DERIV                               &  !<-- Specified 2nd derivatives (=0) at parent levels
                               ,NUM_LEVS_SEC                            &  !<-- Vertical dimension of SEC_DERIV
                               ,NUM_LEVELS                              &  !<-- # of child target levels to interpolate to
                               ,P_OUTPUT                                &  !<-- Child target pressure values to interpolate to
                               ,VBL_COL_CHILD)                             !<-- Child values of variable returned on P_OUTPUT levels
!
                  ENDDO
                  ENDDO
!
                  KNT_REAL_PTS(N)=KNT_REAL_PTS(N)+NPOINTS_HORIZ*NUM_LEVELS !<-- Total points updated in composite output after 
!                                                                          !    this 3-D real variable was done.
                  DEALLOCATE(SEC_DERIV)
!
!-----------------------------------------------------------------------
!***  For 3-D soil variables the parent uses bilinear interpolation
!***  but must also use the land mask.  The bilinear interpolation
!***  weighting needs to be adjusted to account for water points that
!***  are excluded from the summation.  The code assumes there are
!***  no 3-D water point variables to update (UPDATE_TYPE_CHAR=='S').
!-----------------------------------------------------------------------
!
                ELSEIF(UPDATE_TYPE_CHAR=='L')THEN                          !<-- 3-D H-pt variable that is soil
!
!-----------------------------------------------------------------------
!
!!!               FNAME=TRIM(FIELD_NAME)
                  FNAME=FIELD_NAME
!
                  DO L=KLO,KHI                                             !<-- Loop through the soil layers
!
                    DO J=J1,J2
                      J_SOUTH=J_PARENT_SOUTH(J)
                      J_NORTH=J_PARENT_NORTH(J)
!
                      DO I=I1,I2
                        I_WEST=I_PARENT_WEST(I)
                        I_EAST=I_PARENT_EAST(I)
!
                        KNT_REAL_PTS(N)=KNT_REAL_PTS(N)+1                  !<-- Total real points updated in composite output.
!
!!!                     CHILD_UPDATE_DATA%TASKS(N)%DATA_REAL(KNT_REAL_PTS(N))=  &  !<-- Parent's 3-D soil variable interpolated 
                        SUM_PROD=PROD_LWGT_SW(I,J)+PROD_LWGT_SE(I,J)            &
                                +PROD_LWGT_NW(I,J)+PROD_LWGT_NE(I,J)
!
                        IF(ABS(SUM_PROD)<1.E-5)THEN
                          IF(FNAME=='SMC'.OR.FNAME=='SH2O')THEN
                            CHILD_UPDATE_DATA%TASKS(N)%DATA_REAL(KNT_REAL_PTS(N))=1.0
                          ELSEIF(FNAME=='STC')THEN
                            CHILD_UPDATE_DATA%TASKS(N)%DATA_REAL(KNT_REAL_PTS(N))=273.16
                          ENDIF
                        ELSE
                          CHILD_UPDATE_DATA%TASKS(N)%DATA_REAL(KNT_REAL_PTS(N))= &  !<---
                                    PROD_LWGT_SW(I,J)*ARRAY_3D(I_WEST,J_SOUTH,L) &  !    Parent's 3-D soil variable interpolated 
                                   +PROD_LWGT_SE(I,J)*ARRAY_3D(I_EAST,J_SOUTH,L) &  !    horizontally to the moving nest's
                                   +PROD_LWGT_NW(I,J)*ARRAY_3D(I_WEST,J_NORTH,L) &  !    update location using the land mask.
                                   +PROD_LWGT_NE(I,J)*ARRAY_3D(I_EAST,J_NORTH,L)    !<--
                        ENDIF
!
                    ENDDO
                  ENDDO
                ENDDO
!
!-----------------------------------------------------------------------
!
                ENDIF  soil_or_not
!
!-----------------------------------------------------------------------
!
              ENDIF q2
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
            ENDIF dims_2_or_3
!
!-----------------------------------------------------------------------
!
          ENDDO field_loop
!
!-----------------------------------------------------------------------
!
          DEALLOCATE(I_PARENT)
          DEALLOCATE(I_PARENT_EAST)
          DEALLOCATE(I_PARENT_WEST)
!
          DEALLOCATE(J_PARENT)
          DEALLOCATE(J_PARENT_NORTH)
          DEALLOCATE(J_PARENT_SOUTH)
!
          DEALLOCATE(WGHT_SW)
          DEALLOCATE(WGHT_NW)
          DEALLOCATE(WGHT_NE)
          DEALLOCATE(WGHT_SE)
!
          DEALLOCATE(LOG_PBOT   )
          DEALLOCATE(LOG_PTOP   )
          DEALLOCATE(PINT_INTERP)
          DEALLOCATE( PHI_INTERP)
          DEALLOCATE(PMID_INTERP)
          DEALLOCATE(PMID_CHILD )
          DEALLOCATE(PINT_CHILD )
          DEALLOCATE(PD_CHILD   )
          DEALLOCATE(PD_INTERP  )
          DEALLOCATE(VBL_INTERP )
!
          NULLIFY(P3D_INPUT)
          NULLIFY(P3D_OUTPUT)
!
          IF(FLAG_H_OR_V=='H')THEN
            DEALLOCATE(PROD_LWGT_NE)
            DEALLOCATE(PROD_LWGT_NW)
            DEALLOCATE(PROD_LWGT_SE)
            DEALLOCATE(PROD_LWGT_SW)
            DEALLOCATE(PROD_SWGT_NE)
            DEALLOCATE(PROD_SWGT_NW)
            DEALLOCATE(PROD_SWGT_SE)
            DEALLOCATE(PROD_SWGT_SW)
          ENDIF
!
          IF(FLAG_H_OR_V=='V')THEN
!
            DEALLOCATE(I_PARENT_EAST_H)
            DEALLOCATE(I_PARENT_WEST_H)
            DEALLOCATE(J_PARENT_NORTH_H)
            DEALLOCATE(J_PARENT_SOUTH_H)
!
            DEALLOCATE(I_CHILD_ON_PARENT_H)
            DEALLOCATE(J_CHILD_ON_PARENT_H)
!
            DEALLOCATE(WGHT_NE_H)
            DEALLOCATE(WGHT_NW_H)
            DEALLOCATE(WGHT_SE_H)
            DEALLOCATE(WGHT_SW_H)
!
          ENDIF
!
        ENDDO iter_loop
!
!-----------------------------------------------------------------------
!***  The parent task sends its update data to this moving nest task.
!***  The parent only sends to a moving nest after updating both H
!***  and V points so that all data can be sent to each nest task
!***  in a single message.
!-----------------------------------------------------------------------
!
        IF(FLAG_H_OR_V=='V')THEN
!
          CHILDTASK=CHILD_TASK_RANKS(PTR_X%TASK_ID)
          ITAG=KNT_REAL_PTS(N)+NTIMESTEP_CHILD                             !<-- Tag that changes for data size and time
!
          CALL MPI_ISSEND(CHILD_UPDATE_DATA%TASKS(N)%DATA_REAL          &  !<-- Internal real update data for moving nest task N
                         ,KNT_REAL_PTS(N)                               &  !<-- # of real words in the data string
                         ,MPI_REAL                                      &  !<-- Datatype
                         ,CHILDTASK                                     &  !<-- Local intracom rank of nest task to recv data
                         ,ITAG                                          &  !<-- Unique MPI tag
                         ,COMM_TO_MY_CHILD                              &  !<-- MPI intracommunicator
                         ,HANDLE_UPDATE(PTR_X%TASK_ID)                  &  !<-- Handle for ISend to child task
                         ,IERR )
!
          IF(KNT_INTEGER_PTS(N)>0)THEN
            ITAG=KNT_INTEGER_PTS(N)+NTIMESTEP_CHILD                        !<-- Tag that changes for data size and time
!
            CALL MPI_ISSEND(CHILD_UPDATE_DATA%TASKS(N)%DATA_INTEGER     &  !<-- Internal integer update data for moving nest task N
                           ,KNT_INTEGER_PTS(N)                          &  !<-- # of integer words in the data string
                           ,MPI_INTEGER                                 &  !<-- Datatype
                           ,CHILDTASK                                   &  !<-- Local intracom rank of nest task to recv data
                           ,ITAG                                        &  !<-- Unique MPI tag
                           ,COMM_TO_MY_CHILD                            &  !<-- MPI intracommunicator
                           ,HANDLE_UPDATE(PTR_X%TASK_ID)                &  !<-- Handle for ISend to child task
                           ,IERR )
          ENDIF
!
        ENDIF
!
!-----------------------------------------------------------------------
!***  Point at the next link of the linked list holding the 
!***  update task ID and index limits on the next task.
!-----------------------------------------------------------------------
!
        PTR_X=>PTR_X%NEXT_LINK
!
!-----------------------------------------------------------------------
!
      ENDDO ctask_loop
!
!-----------------------------------------------------------------------
!***  All of the combined H and V update data has been sent by this
!***  parent task to each appropriate task on this nest so deallocate
!***  the array holding the number of update points on each nest task.
!-----------------------------------------------------------------------
!
      IF(FLAG_H_OR_V=='V')THEN
        DEALLOCATE(KNT_REAL_PTS)
        DEALLOCATE(KNT_INTEGER_PTS)
        DEALLOCATE(NUM_ITER)
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE PARENT_UPDATES_MOVING
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE PARENT_READS_MOVING_CHILD_TOPO(MY_DOMAIN_ID            &
                                               ,NUM_MOVING_CHILDREN     &
                                               ,LINK_MRANK_RATIO        &
                                               ,LIST_OF_RATIOS          &
                                               ,M_NEST_RATIO            &
                                               ,KOUNT_RATIOS_MN         &
                                               ,GLOBAL_TOP_PARENT       &
                                               ,IM_1,JM_1               &
                                               ,TPH0_1,TLM0_1           &
                                               ,SB_1,WB_1               &
                                               ,RECIP_DPH_1,RECIP_DLM_1 &
                                               ,GLAT,GLON               &
                                               ,NEST_FIS_ON_PARENT_BNDS &
                                               ,NEST_FIS_ON_PARENT      &
                                               ,NEST_FIS_V_ON_PARENT    &
                                               ,IDS,IDE,IMS,IME,ITS,ITE &
                                               ,JDS,JDE,JMS,JME,JTS,JTE)
!
!-----------------------------------------------------------------------
!***  Parents of moving nests must fill their own domains with the
!***  full resolution topography of those children.  That data spans
!***  the entire domain of the uppermost parent.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      INTEGER(kind=KINT),INTENT(IN) :: IM_1,JM_1                        &  !<-- Dimensions of the uppermost parent domain
                                      ,IDS,IDE,JDS,JDE                  &  !<-- This parent domain's index limits
                                      ,IMS,IME,JMS,JME                  &  !<-- This parent tasks's memory limits
                                      ,ITS,ITE,JTS,JTE                  &  !<-- This parent tasks's integration limits
                                      ,KOUNT_RATIOS_MN                  &  !<-- # of space ratios of children to upper parent
                                      ,MY_DOMAIN_ID                     &  !<-- This parent domain's ID
                                      ,NUM_MOVING_CHILDREN                 !<-- # of moving children on this parent domain
!
      INTEGER(kind=KINT),DIMENSION(1:NUM_MOVING_CHILDREN),INTENT(IN) :: &
                                                       LINK_MRANK_RATIO &  !<-- Each child asociated with rank of space ratio in list
                                                      ,LIST_OF_RATIOS   &  !<-- The list of different space ratios
                                                      ,M_NEST_RATIO        !<-- Associate each child with its upper parent space ratio
!
      REAL(kind=KFPT),INTENT(IN) :: RECIP_DPH_1,RECIP_DLM_1             &  !<-- Reciprocal of uppermost domain grid increments (radians)
                                   ,TLM0_1,TPH0_1                       &  !<-- Central geo lat/lon of uppermost domain (radians; east/north)
                                   ,SB_1,WB_1                              !<-- Rotated lat/lon of south/west boundary (radians; north/east)
!
      REAL(kind=KFPT),DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: GLAT,GLON   !<-- Geographic lat/lon (radians) on parent grid
!
      LOGICAL(kind=KLOG),INTENT(IN) :: GLOBAL_TOP_PARENT                   !<-- Is the uppermost parent domain global?
!
      TYPE(BNDS_2D),DIMENSION(1:KOUNT_RATIOS_MN),INTENT(OUT) ::         &
                                               NEST_FIS_ON_PARENT_BNDS     !<-- Parent subdomain index limits of nest-res topo data
!
      TYPE(REAL_DATA_2D),DIMENSION(1:KOUNT_RATIOS_MN),INTENT(INOUT) ::  &
                                                   NEST_FIS_ON_PARENT   &  !<-- Nest-res topo data on the parent task subdomain
                                                  ,NEST_FIS_V_ON_PARENT    !<-- Nest-res topo data at V on the parent task subdomain
!
!--------------------
!***  Local Variables
!--------------------
!
      INTEGER(kind=KINT) :: I,ICORNER,IDIM,IEND,ISTART,IUNIT_FIS_NEST   &
                           ,J,JCORNER,JDIM,JEND,JSTART,JSTOP,LOR,N,NN
!
      INTEGER(kind=KINT) :: I_COUNT_DATA,J_COUNT_DATA                   &
                           ,I_EXTRA_DATA,J_EXTRA_DATA                   &
                           ,NCID,NCTYPE,NDIMS,VAR_ID
!
      INTEGER(kind=KINT) :: IERR,ISTAT
!
      INTEGER(kind=KINT),DIMENSION(1:2) :: DIM_IDS
!
      REAL(kind=KFPT) :: GBL,REAL_I_NE,REAL_I_SW,REAL_J_NE,REAL_J_SW    &
                        ,VAL_NE
!
      REAL(kind=KFPT),DIMENSION(:),ALLOCATABLE :: COL,ROW
!
      CHARACTER(len=2) :: ID_TOPO_FILE
      CHARACTER(len=9) :: FILENAME
      CHARACTER(len=15) :: VNAME
!
      LOGICAL(kind=KLOG) :: OPENED
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Parents with moving nests need to know those nests' topography
!***  at the nests' own resolutions for the hydrostatic adjustment
!***  that must take place when the parents interpolate their data
!***  to moving nest grid points.  For the sake of generality all
!***  of those nest-resolution datasets must span the domain of the
!***  uppermost parent which is the true maximum range of any nest's
!***  motion.
!
!***  So each parent with moving nests must:
!***   (1) Know how many different space resolutions its moving
!***       children use;
!***   (2) Associate each resolution with the appropriate moving
!***       child using the nest-to-uppermost parent space ratio
!***       that the user specified in each moving nest's configure
!***       file;
!***   (3) Have each of its forecast tasks read in its own piece of
!***       each different resolution of topography data needed by
!***       all of its moving children.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Each parent task obtains the real I,J on the uppermost parent
!***  of the SW and NE corners of each of their subdomains.  Given
!***  the known resolution of the nest topography data the parent
!***  tasks can then extract and save only that data which covers
!***  their own subdomain.
!-----------------------------------------------------------------------
!
!----------------------------------------
!***  SW corner of parent task subdomain
!----------------------------------------
!
      ICORNER=MAX(IMS,IDS)                                                 !<-- Parent task covers its halos with data too since
      JCORNER=MAX(JMS,JDS)                                                 !    the moving nest boundaries can extend into them.
!
      IF(MY_DOMAIN_ID==1.AND.GLOBAL_TOP_PARENT)THEN                        !<-- The current parent domain is global.
        IF(ITS==IDS)THEN
          ICORNER=ICORNER+1.                                               !<-- Past buffer row to Intl Dateline
        ENDIF
        IF(JTS==JDS)THEN
          JCORNER=JCORNER+1.
        ENDIF
      ENDIF
!
      CALL LATLON_TO_IJ(GLAT(ICORNER,JCORNER)                           &  !<-- Geographic lat (radians) of parent task's SW corner
                       ,GLON(ICORNER,JCORNER)                           &  !<-- Geographic lon (radians) of parent task's SW corner
                       ,TPH0_1,TLM0_1                                   &  !<-- Geographic lat,lon of upper parent's central point
                       ,SB_1,WB_1                                       &  !<-- Rotated lat/lon of upper parent's SW corner
                       ,RECIP_DPH_1,RECIP_DLM_1                         &
                       ,GLOBAL_TOP_PARENT                               &  !<-- Is the uppermost parent domain global?
                       ,REAL_I_SW                                       &  !<-- Uppermost parent I of this task's SW corner
                       ,REAL_J_SW)                                         !<-- Uppermost parent J of this task's SW corner
!
!----------------------------------------
!***  NE corner of parent task subdomain
!----------------------------------------
!
      ICORNER=MIN(IME,IDE)                                                 !<-- Parent task covers its halos with data too since
      JCORNER=MIN(JME,JDE)                                                 !    the moving nest boundaries can extend into them.
!
      IF(MY_DOMAIN_ID==1.AND.GLOBAL_TOP_PARENT)THEN                        !<-- The current parent domain is global.
        IF(ITE==IDE)THEN
          ICORNER=ICORNER-1.                                               !<-- Past buffer row to Intl Dateline
        ENDIF
        IF(JTE==JDE)THEN
          JCORNER=JCORNER-1.
        ENDIF
      ENDIF
!
      CALL LATLON_TO_IJ(GLAT(ICORNER,JCORNER)                           &
                       ,GLON(ICORNER,JCORNER)                           &
                       ,TPH0_1,TLM0_1                                   &
                       ,SB_1,WB_1                                       &
                       ,RECIP_DPH_1,RECIP_DLM_1                         &
                       ,GLOBAL_TOP_PARENT                               &
                       ,REAL_I_NE                                       &
                       ,REAL_J_NE)
!
!-----------------------------------------------------------------------
      nr_loop: DO N=1,KOUNT_RATIOS_MN                                      !<-- Loop through the different parent-child space ratios
!-----------------------------------------------------------------------
!
        LOR=LIST_OF_RATIOS(N)
!
        IF(GLOBAL_TOP_PARENT)THEN
          GBL=1.                                                           !<-- Account for the extra row that surrounds global domains.
        ELSE
          GBL=0.
        ENDIF
!
        ISTART=NINT((REAL_I_SW-1.-GBL)*LOR+1.)                             !<-- I index in sfc data at W bndry of this parent task
        JSTART=NINT((REAL_J_SW-1.-GBL)*LOR+1.)                             !<-- J index in sfc data at S bndry of this parent task
!
        IEND=NINT((REAL_I_NE-1.-GBL)*LOR+1.)                               !<-- I index in nest sfc data at E bndry (H) of this parent task
        JEND=NINT((REAL_J_NE-1.-GBL)*LOR+1.)                               !<-- J index in nest sfc data at N bndry (H) of this parent task
!
        I_COUNT_DATA=IEND-ISTART+1
        J_COUNT_DATA=JEND-JSTART+1
!
!-----------------------------------------------------------------------
!
        NEST_FIS_ON_PARENT_BNDS(N)%LBND1=ISTART                            !<-- Array limits in nest-resolution topography data
        NEST_FIS_ON_PARENT_BNDS(N)%UBND1=IEND                              !    for region covering this parent task's subdomain.
        NEST_FIS_ON_PARENT_BNDS(N)%LBND2=JSTART                            !
        NEST_FIS_ON_PARENT_BNDS(N)%UBND2=JEND                              !<--
!
!-----------------------------------------------------------------------
!***  Each parent task opens and reads the topography file.
!-----------------------------------------------------------------------
!
        IF(N<=9)THEN
          NN=LOR
          IF(NN<=9)THEN
            WRITE(ID_TOPO_FILE,'(I1.1)')NN
          ELSEIF(NN>=10)THEN
            WRITE(ID_TOPO_FILE,'(I2.2)')NN
          ENDIF
        ELSE
          WRITE(0,*)' User specified more than 9 different'            &
                   ,' moving nest resolutions!!!'
          WRITE(0,*)' ABORTING'
          CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
        ENDIF
!
        FILENAME='FIS_'//TRIM(ID_TOPO_FILE)//'.nc'
!
        CALL CHECK(NF90_OPEN(FILENAME,NF90_NOWRITE,NCID))                  !<-- Open the FIS external netCDF file for Nth space ratio.
!
!-----------------------------------------------------------------------
!***  Each task allocates its space for holding its moving children's
!***  topography at their resolution.
!-----------------------------------------------------------------------
!
        IF(ASSOCIATED(NEST_FIS_ON_PARENT(N)%DATA))THEN
          DEALLOCATE(NEST_FIS_ON_PARENT(N)%DATA,stat=ISTAT)
        ENDIF
        IF(ASSOCIATED(NEST_FIS_V_ON_PARENT(N)%DATA))THEN
          DEALLOCATE(NEST_FIS_V_ON_PARENT(N)%DATA,stat=ISTAT)
        ENDIF
!
        ALLOCATE(NEST_FIS_ON_PARENT(N)%DATA(ISTART:IEND,JSTART:JEND))
        ALLOCATE(NEST_FIS_V_ON_PARENT(N)%DATA(ISTART:IEND,JSTART:JEND))
!
!-----------------------------------------------------------------------
!***  Save only those points in the topography data for resolution N
!***  that cover this parent task's subdomain.
!***  Begin with the nest-resolution topography at H points.
!-----------------------------------------------------------------------
!
        CALL CHECK(NF90_INQUIRE_VARIABLE(NCID,3,VNAME,NCTYPE  &            !<-- Topography is the 3rd variable in the file.
                                        ,NDIMS,DIM_IDS))
        CALL CHECK(NF90_INQ_VARID(NCID,VNAME,VAR_ID))
!
        CALL CHECK(NF90_GET_VAR(NCID,VAR_ID                                         &  !<-- Extract the values
                               ,NEST_FIS_ON_PARENT(N)%DATA(ISTART:IEND,JSTART:JEND) &  !    of nest-resolution
                               ,start=(/ISTART,JSTART/)                             &  !    topography from the
                               ,count=(/I_COUNT_DATA,J_COUNT_DATA/)))                  !    external file.
!
!-----------------------------------------------------------------------
!***  For the nest topography values at V points we can begin by
!***  averaging the values at H points.
!-----------------------------------------------------------------------
!
        DO J=JSTART,JEND-1
        DO I=ISTART,IEND-1
          NEST_FIS_V_ON_PARENT(N)%DATA(I,J)=(NEST_FIS_ON_PARENT(N)%DATA(I,J)     &
                                            +NEST_FIS_ON_PARENT(N)%DATA(I+1,J)   &
                                            +NEST_FIS_ON_PARENT(N)%DATA(I,J+1)   &
                                            +NEST_FIS_ON_PARENT(N)%DATA(I+1,J+1) &
                                            )*0.25
        ENDDO
        ENDDO
!                
!-----------------------------------------------------------------------
!***  The V row at J=J_END is north of the H row at J=J_END.  
!***  The V column at I=I_END is east of the H column at I=I_END.
!***  This means we need to read in extra values to get those 
!***  V points on the north and east edges of the parent tasks.
!-----------------------------------------------------------------------
!
        ALLOCATE(ROW(ISTART:IEND))
        ALLOCATE(COL(JSTART:JEND))
!
        I_EXTRA_DATA=ISTART+I_COUNT_DATA                                   !<-- 1 column east of task's saved H data
        J_EXTRA_DATA=JSTART+J_COUNT_DATA                                   !<-- 1 row north of task's saved H data
!
!-----------------------------------------------------------------------
!***  Fill in values of nest topography on V points one row north
!***  of the northern limit of the H-point topography saved on
!***  this parent task.
!-----------------------------------------------------------------------
!
        IF(JTE<JDE)THEN
          CALL CHECK(NF90_GET_VAR(NCID,VAR_ID                           &  !<-- Extract the values
                                 ,ROW(ISTART:IEND)                      &  !    1 row north of the
                                 ,start=(/ISTART,J_EXTRA_DATA/)         &  !    northernmost H pts
                                 ,count=(/I_COUNT_DATA,1/)))               !    on this parent task.
!
          DO I=ISTART,IEND-1
            NEST_FIS_V_ON_PARENT(N)%DATA(I,JEND)=(NEST_FIS_ON_PARENT(N)%DATA(I,JEND)    &
                                                  +NEST_FIS_ON_PARENT(N)%DATA(I+1,JEND) &
                                                  +ROW(I)+ROW(I+1) )*0.25
          ENDDO
!
        ELSE
!
          DO I=ISTART,IEND-1
            NEST_FIS_V_ON_PARENT(N)%DATA(I,JEND)=NEST_FIS_V_ON_PARENT(N)%DATA(I,JEND-1)
          ENDDO
!
        ENDIF
!
!-----------------------------------------------------------------------
!***  Fill in values of nest topography on V points one column east
!***  of the eastern limit of the H-point topography saved on this
!***  parent task.
!-----------------------------------------------------------------------
!
        IF(ITE<IDE)THEN
          CALL CHECK(NF90_GET_VAR(NCID,VAR_ID                           &  !<-- Extract the values
                                 ,COL(JSTART:JEND)                      &  !    1 column east of the
                                 ,start=(/I_EXTRA_DATA,JSTART/)         &  !    easternmost H pts
                                 ,count=(/1,J_COUNT_DATA/)))               !    on this parent task.
!
          DO J=JSTART,JEND-1
            NEST_FIS_V_ON_PARENT(N)%DATA(IEND,J)=(NEST_FIS_ON_PARENT(N)%DATA(IEND,J)   &
                                                 +NEST_FIS_ON_PARENT(N)%DATA(IEND,J+1) &
                                                 +COL(J)+COL(J+1) )*0.25
          ENDDO
!
        ELSE
!          
          DO J=JSTART,JEND-1
            NEST_FIS_V_ON_PARENT(N)%DATA(IEND,J)=NEST_FIS_V_ON_PARENT(N)%DATA(IEND-1,J)
          ENDDO
!
        ENDIF
!
!-----------------------------------------------------------------------
!***  The extreme northeast V point on all parent tasks that do not
!***  lie along the northern and eastern side of the parent domain.
!-----------------------------------------------------------------------
!
        IF(ITE<IDE.AND.JTE<JDE)THEN
          CALL CHECK(NF90_GET_VAR(NCID,VAR_ID                            & 
                                 ,VAL_NE                                 &
                                 ,start=(/I_EXTRA_DATA,J_EXTRA_DATA/))) 
!
          NEST_FIS_V_ON_PARENT(N)%DATA(IEND,JEND)=(NEST_FIS_ON_PARENT(N)%DATA(IEND,JEND) &
                                                    +ROW(IEND)+COL(JEND)+VAL_NE)*0.25
!
!-----------------------------------------------------------------------
!***  The extreme northeast V point on all parent tasks that DO
!***  lie along the northern and eastern side of the parent domain.
!-----------------------------------------------------------------------
!
        ELSE
!
          NEST_FIS_V_ON_PARENT(N)%DATA(IEND,JEND)=0.
!
        ENDIF
!
!-----------------------------------------------------------------------
!***  Set all tiny values identically to zero.
!-----------------------------------------------------------------------
!
        DO J=JSTART,JEND
        DO I=ISTART,IEND
          IF(ABS(NEST_FIS_ON_PARENT(N)%DATA(I,J))<1.E-2)THEN
            NEST_FIS_ON_PARENT(N)%DATA(I,J)=0.
          ENDIF
!
          IF(ABS(NEST_FIS_V_ON_PARENT(N)%DATA(I,J))<1.E-2)THEN
            NEST_FIS_V_ON_PARENT(N)%DATA(I,J)=0.
          ENDIF
        ENDDO
        ENDDO
!
        DEALLOCATE(ROW,COL)
!
        CALL CHECK(NF90_CLOSE(NCID))                                       !<-- Close the external netCDF file.
!
!-----------------------------------------------------------------------
!
      ENDDO nr_loop
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE PARENT_READS_MOVING_CHILD_TOPO
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE LATLON_TO_IJ(GLAT,GLON                                 &
                             ,TPH0_1,TLM0_1                             &
                             ,SB_1,WB_1                                 &
                             ,RECIP_DPH_1,RECIP_DLM_1                   &
                             ,GLOBAL_TOP_PARENT                         &
                             ,REAL_I,REAL_J)
! 
!-----------------------------------------------------------------------
!***  Determine the Real values of I and J on a rotated B Grid given
!***  the geographic latitude/longitude.  In this routine the domain
!***  is specifically that of the uppermost parent.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      REAL(kind=KFPT),INTENT(IN) :: GLAT                                &  !<-- Geographic latitude (radians, positive north)
                                   ,GLON                                &  !<-- Geographic longitude (radians, positive east)
                                   ,TPH0_1                              &  !<-- Central latitude of uppermost parent (radians, north)
                                   ,TLM0_1                              &  !<-- Central longitude of uppermost parent (radians, east)
                                   ,SB_1                                &  !<-- Rotated lat of upper parent's south boundary (radians, north)
                                   ,WB_1                                &  !<-- Rotated lon of upper parent's west boundary (radians, east)
                                   ,RECIP_DPH_1                         &  !<-- Reciprocal of upper parent's grid increment (radians) in J
                                   ,RECIP_DLM_1                            !<-- Reciprocal of upper parent's grid increment (radians) in I
!
      LOGICAL(kind=KLOG),INTENT(IN) :: GLOBAL_TOP_PARENT                   !<-- Is the uppermost parent domain global?
!
      REAL(kind=KFPT),INTENT(OUT) :: REAL_I                             &  !<-- Real I on uppermost parent grid for GLAT,GLON
                                    ,REAL_J                                !<-- Real J on uppermost parent grid for GLAT,GLON
!
!---------------------
!***  Local Variables
!---------------------
!
      REAL(kind=KDBL),SAVE :: PI= 3.14159265359                         &
                             ,PI2=6.28318530718
!
      REAL(kind=KDBL) :: X,Y,Z
!
      REAL(kind=KFPT) :: SIN_DIFF,TLAT,TLON 
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      IF(GLOBAL_TOP_PARENT)THEN
        TLAT=GLAT
        TLON=GLON
!
      ELSE
        X=COS(TPH0_1)*COS(GLAT)*COS(GLON-TLM0_1)+SIN(TPH0_1)*SIN(GLAT)
        SIN_DIFF=SIN(GLON-TLM0_1)
        IF(ABS(SIN_DIFF)<1.E-6)SIN_DIFF=0.
        Y=COS(GLAT)*SIN_DIFF
        Z=COS(TPH0_1)*SIN(GLAT)-SIN(TPH0_1)*COS(GLAT)*COS(GLON-TLM0_1)
        TLAT=ATAN(Z/SQRT(X*X+Y*Y))
        TLON=ATAN(Y/X)
        IF(X<0.)THEN
          TLON=TLON+PI
          IF(Y<=0..AND.GLON<0.)THEN
            TLON=TLON-PI2
          ENDIF
        ENDIF
!
      ENDIF
!
      REAL_I=(TLON-WB_1)*RECIP_DLM_1+1
      REAL_J=(TLAT-SB_1)*RECIP_DPH_1+1
!
      IF(GLOBAL_TOP_PARENT)THEN
        REAL_I=REAL_I+1.                                                   !<--  Account for the extra row that
        REAL_J=REAL_J+1.                                                   !     surrounds global domains.
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE LATLON_TO_IJ
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE PARENT_2WAY_BOOKKEEPING(PARENT_CHILD_SPACE_RATIO       &
                                        ,NUM_CHILD_TASKS                &
                                        ,CHILD_TASK_LIMITS              &
                                        ,IM_CHILD                       &
                                        ,JM_CHILD                       &
                                        ,I_PARENT_SW                    &
                                        ,J_PARENT_SW                    &
                                        ,N_BLEND_H_CHILD                &
                                        ,N_BLEND_V_CHILD                &
                                        ,N_STENCIL_H_CHILD              & ! input
                                        ,N_STENCIL_V_CHILD              & !   ^
                                        ,N_STENCIL_SFC_H_CHILD          & !   |
                                        ,N_STENCIL_SFC_V_CHILD          & !   |
                                        ,ITS,ITE,JTS,JTE                & !   |
!                                                                         ! -----
                                        ,NTASKS_UPDATE_CHILD            & !   |
                                        ,CHILD_TASKS_2WAY_UPDATE        & ! output
                                                             )
!
!-----------------------------------------------------------------------
!***  The parent domain determines which of its task subdomains and
!***  which points on those subdomains are updated by its children.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument variables
!------------------------
!
      INTEGER(kind=KINT),INTENT(IN) :: I_PARENT_SW                      &  !<-- Parent I of nest domain's SW corner
                                      ,J_PARENT_SW                      &  !<-- Parent J of nest domain's SW corner
                                      ,N_STENCIL_H_CHILD                &  !<-- Stencil width of child h to parent H averaging
                                      ,N_STENCIL_V_CHILD                &  !<-- Stencil width of child v to parent V averaging
                                      ,N_STENCIL_SFC_H_CHILD            &  !<-- Stencil width of child pd to parent H averaging
                                      ,N_STENCIL_SFC_V_CHILD            &  !<-- Stencil width of child pd to parent V averaging
                                      ,NUM_CHILD_TASKS                  &  !<-- # of forecast tasks on this child domain
                                      ,PARENT_CHILD_SPACE_RATIO            !<-- Ratio of parent grid increment to child's
!
      INTEGER(kind=KINT),INTENT(IN) :: ITS,ITE,JTS,JTE                     !<-- Integration limits of this parent task's subdomain     
!
      INTEGER(kind=KINT),INTENT(IN) :: IM_CHILD                         &  !<-- I dimension of this child domain
                                      ,JM_CHILD                         &  !<-- J dimension of this child domain
                                      ,N_BLEND_H_CHILD                  &  !<-- Blending region width for this child's domain H pts
                                      ,N_BLEND_V_CHILD                     !<-- Blending region width for this child's domain V pts
!
      INTEGER(kind=KINT),DIMENSION(1:4,NUM_CHILD_TASKS),INTENT(IN) ::   &
                                                     CHILD_TASK_LIMITS     !<-- ITS,ITE,JTS,JTE for each child forecast task
!
      INTEGER(kind=KINT),INTENT(INOUT) :: NTASKS_UPDATE_CHILD              !<-- Parent task recvs 2-way update from this many child tasks
!
      TYPE(CHILD_UPDATE_LINK),TARGET,INTENT(INOUT) ::                   &
                                             CHILD_TASKS_2WAY_UPDATE       !<-- Ranks of tasks on 2-way child that update this parent task
!                                                                               and the parent points updated.
!---------------------
!***  Local variables
!---------------------
!
      INTEGER(kind=KINT) :: ISTAT,KOUNT,N,N_BLEND_CHILD,N_STENCIL_0     &
                           ,NPTS_PARENT_UPDATE,NT
!
      INTEGER(kind=KINT) :: I1,I2,J1,J2
!
      INTEGER(kind=KINT),DIMENSION(1:4) :: N_STENCIL_X
!
      REAL(kind=KFPT) :: RECIP_RATIO
!
      REAL(kind=KFPT) :: CHILD_ISTART_ON_PARENT,CHILD_IEND_ON_PARENT    &
                        ,CHILD_JSTART_ON_PARENT,CHILD_JEND_ON_PARENT    &
!
                        ,IDE_CHILD_ON_PARENT,IDS_CHILD_ON_PARENT        &
                        ,JDE_CHILD_ON_PARENT,JDS_CHILD_ON_PARENT        &
!
                        ,ITE_CHILD_ON_PARENT,ITS_CHILD_ON_PARENT        &
                        ,JTE_CHILD_ON_PARENT,JTS_CHILD_ON_PARENT        &
!
                        ,LIMIT_EAST_H,LIMIT_WEST_H                      &
                        ,LIMIT_NORTH_H,LIMIT_SOUTH_H
!
      TYPE(CHILD_UPDATE_LINK),POINTER :: HEAD,TAIL
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RECIP_RATIO=1./PARENT_CHILD_SPACE_RATIO
!
      N_BLEND_CHILD=N_BLEND_H_CHILD                                        !<-- Domain blending regions must be the same for H,V.
!
!-----------------------------------------------------------------------
!***  What are this child's domain limits on the parent grid?
!-----------------------------------------------------------------------
!
      IDS_CHILD_ON_PARENT=REAL(I_PARENT_SW)                                !<-- 2-way child's W boundary in parent's grid space       
      IDE_CHILD_ON_PARENT=REAL(I_PARENT_SW)                             &  !<-- 2-way child's E boundary in parent's grid space
                          +(IM_CHILD-1)*RECIP_RATIO
      JDS_CHILD_ON_PARENT=REAL(J_PARENT_SW)                                !<-- 2-way child's S boundary in parent's grid space
      JDE_CHILD_ON_PARENT=REAL(J_PARENT_SW)                             &  !<-- 2-way child's N boundary in parent's grid space
                          +(JM_CHILD-1)*RECIP_RATIO
!
!-----------------------------------------------------------------------
!***  We want to limit the child points that can be used for
!***  computing the 2-way data.  For now do not use any child
!***  points in the averaging stencil that lie in the child's
!***  boundary blending rows.  Stencils can vary for h-->H,
!***  v-->V, pd-->H, and pd-->V.  Determine the set of parent
!***  target I's and J's common to all the stencils and use
!***  that to ensure that the same parent I,J indices are used
!***  for both H and V variables.
!-----------------------------------------------------------------------
!
      N_STENCIL_X(1)=N_STENCIL_H_CHILD
      N_STENCIL_X(2)=N_STENCIL_V_CHILD
      N_STENCIL_X(3)=N_STENCIL_SFC_H_CHILD
      N_STENCIL_X(4)=N_STENCIL_SFC_V_CHILD
!
!-----------------------------------------------------------------------
!***  Deallocate the linked list of child update tasks if it already
!***  exists.  This is relevant only for moving nests.  This routine
!***  is called only once for static nests when the linked list does
!***  not yet exist.
!-----------------------------------------------------------------------
!
      KOUNT=0
      HEAD=>CHILD_TASKS_2WAY_UPDATE                                        !<-- Point at the top of the linked list
!
      dealloc: DO
!
        KOUNT=KOUNT+1
        TAIL=>NULL()
        IF(ASSOCIATED(HEAD%NEXT_LINK))THEN
          TAIL=>HEAD%NEXT_LINK                                             !<-- If another link exists, point at it.
        ENDIF
!
        IF(KOUNT>1)THEN                                                    !<-- Do not deallocate the topmost object's memory
          DEALLOCATE(HEAD%TASK_ID)
          DEALLOCATE(HEAD%IL)     
          DEALLOCATE(HEAD%JL)     
          DEALLOCATE(HEAD%NUM_PTS_UPDATE_HZ)
          DEALLOCATE(HEAD,stat=ISTAT)                                      !<-- Deallocate the current link.
          IF(ISTAT/=0)THEN
            WRITE(0,*)' Failed to deallocate link #',KOUNT,' in 2-way linked list!'
            WRITE(0,*)' Aborting!!'
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
          ENDIF
        ENDIF
!
        IF(ASSOCIATED(TAIL))THEN                                           !<-- If true, another link exists.
          HEAD=>TAIL                                                       !<-- Reset so that the head is at the new link.
        ELSE
          EXIT dealloc                                                     !<-- No further links exist.
        ENDIF
!
      ENDDO dealloc
!
!-----------------------------------------------------------------------
!***  Only the top of the list remains.  Point at it.
!-----------------------------------------------------------------------
!
      HEAD=>CHILD_TASKS_2WAY_UPDATE
      HEAD%NEXT_LINK=>NULL()                                               !<-- There is no 'next link' in the list yet.
!
!-----------------------------------------------------------------------
!***  Which if any child tasks will be updating this parent task? 
!***  Begin by finding the subdomain limits of each child task on
!***  the parent domain.
!-----------------------------------------------------------------------
!
      child_tasks: DO NT=1,NUM_CHILD_TASKS
!
!-----------------------------------------------------------------------
!
        ITS_CHILD_ON_PARENT=(CHILD_TASK_LIMITS(1,NT)-CHILD_TASK_LIMITS(1,1))*RECIP_RATIO  &  !<-- Child task NT's starting I on parent grid
                            +REAL(I_PARENT_SW)
        ITE_CHILD_ON_PARENT=(CHILD_TASK_LIMITS(2,NT)-CHILD_TASK_LIMITS(1,1))*RECIP_RATIO  &  !<-- Child task NT's ending I on parent grid
                            +REAL(I_PARENT_SW)
!
        JTS_CHILD_ON_PARENT=(CHILD_TASK_LIMITS(3,NT)-CHILD_TASK_LIMITS(3,1))*RECIP_RATIO  &  !<-- Child task NT's starting J on parent grid
                            +REAL(J_PARENT_SW)
        JTE_CHILD_ON_PARENT=(CHILD_TASK_LIMITS(4,NT)-CHILD_TASK_LIMITS(3,1))*RECIP_RATIO  &  !<-- Child task NT's ending J on parent grid
                            +REAL(J_PARENT_SW)
!
        CHILD_ISTART_ON_PARENT=ITS_CHILD_ON_PARENT
        CHILD_IEND_ON_PARENT  =ITE_CHILD_ON_PARENT
        CHILD_JSTART_ON_PARENT=JTS_CHILD_ON_PARENT
        CHILD_JEND_ON_PARENT  =JTE_CHILD_ON_PARENT
!
!-----------------------------------------------------------------------
!***  Find the common parent target points for all averaging stencils.
!-----------------------------------------------------------------------
!
        DO N=1,4                                                           !<-- There are 4 averaging stencils; see above.
!                          
          N_STENCIL_0=N_STENCIL_X(N)/2                                     !<-- Child's delta I,J from parent update pt to edge of
!                                                                          !    stencil region that will update the parent point.
          LIMIT_WEST_H=REAL(I_PARENT_SW)                                &
                       +(N_BLEND_CHILD+N_STENCIL_0)*RECIP_RATIO
!
          CHILD_ISTART_ON_PARENT=MAX(CHILD_ISTART_ON_PARENT             &  !<--
                                    ,ITS_CHILD_ON_PARENT                &
                                    ,LIMIT_WEST_H )
!
          LIMIT_EAST_H=REAL(I_PARENT_SW)                                &
                   +(IM_CHILD-1-N_BLEND_CHILD-N_STENCIL_0)*RECIP_RATIO
!
          CHILD_IEND_ON_PARENT=MIN(CHILD_IEND_ON_PARENT                 &  !<--
                                  ,ITE_CHILD_ON_PARENT                  &
                                  ,LIMIT_EAST_H ) 
!
          LIMIT_SOUTH_H=REAL(J_PARENT_SW)                               &
                       +(N_BLEND_CHILD+N_STENCIL_0)*RECIP_RATIO
!
          CHILD_JSTART_ON_PARENT=MAX(CHILD_JSTART_ON_PARENT             &  !<--
                                    ,JTS_CHILD_ON_PARENT                &
                                    ,LIMIT_SOUTH_H )
!
          LIMIT_NORTH_H=REAL(J_PARENT_SW)                               &
                   +(JM_CHILD-1-N_BLEND_CHILD-N_STENCIL_0)*RECIP_RATIO
!
          CHILD_JEND_ON_PARENT=MIN(CHILD_JEND_ON_PARENT                 &  !<--
                                  ,JTE_CHILD_ON_PARENT                  &
                                  ,LIMIT_NORTH_H ) 
!
        ENDDO
!
!-----------------------------------------------------------------------
!***  Which if any of this parent task's points are updated by
!***  this child's task NT?
!-----------------------------------------------------------------------
!
        IF(REAL(ITS)<CHILD_IEND_ON_PARENT+EPS                           &
                          .AND.                                         &
           REAL(ITE)>CHILD_ISTART_ON_PARENT-EPS                         &
                          .AND.                                         &
           REAL(JTS)<CHILD_JEND_ON_PARENT+EPS                           &
                          .AND.                                         &
           REAL(JTE)>CHILD_JSTART_ON_PARENT-EPS )THEN
!
!-----------------------------------------------------------------------
!***  Which points on this parent task will be updated by this
!***  child task?  See examples of the logic used here in 
!***  subroutine CHILD_2WAY_BOOKKEEPING.
!-----------------------------------------------------------------------
!
          I1=MAX(ITS,INT(CHILD_ISTART_ON_PARENT+1.-EPS))                   !<-- Lower I limit on parent update region by child task NT.
          I2=MIN(ITE,INT(CHILD_IEND_ON_PARENT+EPS))                        !<-- Upper I limit on parent update region by child task NT.
          J1=MAX(JTS,INT(CHILD_JSTART_ON_PARENT+1.-EPS))                   !<-- Lower J limit on parent update region by child task NT.
          J2=MIN(JTE,INT(CHILD_JEND_ON_PARENT+EPS))                        !<-- Upper J limit on parent update region by child task NT. 
!
          NPTS_PARENT_UPDATE=(I2-I1+1)*(J2-J1+1)
!
          IF(NPTS_PARENT_UPDATE<=0)THEN
            CYCLE child_tasks                                              !<-- No usable 2-way exchange region on this child task.
          ENDIF
!
          NTASKS_UPDATE_CHILD=NTASKS_UPDATE_CHILD+1                        !<-- Save # of child tasks that send 2-way update
!
          IF(NTASKS_UPDATE_CHILD>1)THEN                                    !<-- We need another link in the list.
            ALLOCATE(HEAD%NEXT_LINK)                                       !<-- Create the new link
            HEAD=>HEAD%NEXT_LINK                                           !<-- Point at the new link.
            HEAD%NEXT_LINK=>NULL()                                         !<-- Nullify the link that would follow the new link.
!
            ALLOCATE(HEAD%TASK_ID)                                         !<--
            ALLOCATE(HEAD%IL(1:2))                                         !    Create the components
            ALLOCATE(HEAD%JL(1:2))                                         !    of the new link.
            ALLOCATE(HEAD%NUM_PTS_UPDATE_HZ)                               !<--
          ENDIF
!
!-----------------------------------------------------------------------
!***  In this link of the list save the updating child task's local
!***  rank on its domain as well as the index limits on this parent
!***  task that this child task will update along with the total
!***  number of updated parent task points.
!-----------------------------------------------------------------------
!
          HEAD%TASK_ID=NT-1                                                !<-- Local rank of child task sending 2-way update
!
          HEAD%IL(1)=I1
          HEAD%IL(2)=I2
          HEAD%JL(1)=J1
          HEAD%JL(2)=J2
!
          HEAD%NUM_PTS_UPDATE_HZ=NPTS_PARENT_UPDATE
!
!-----------------------------------------------------------------------
!
        ENDIF
!
!-----------------------------------------------------------------------
!
      ENDDO child_tasks
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE PARENT_2WAY_BOOKKEEPING
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE CHILD_2WAY_BOOKKEEPING(I_SW_PARENT_CURRENT             &
                                       ,J_SW_PARENT_CURRENT             &
                                       ,SPACE_RATIO_MY_PARENT           &
                                       ,NUM_FCST_TASKS_PARENT           &
                                       ,ITS_PARENT_TASKS                &
                                       ,ITE_PARENT_TASKS                &
                                       ,JTS_PARENT_TASKS                &
                                       ,JTE_PARENT_TASKS                &
                                       ,N_BLEND_H                       &
                                       ,N_BLEND_V                       &
                                       ,N_STENCIL_H                     &
                                       ,N_STENCIL_V                     &
                                       ,N_STENCIL_SFC_H                 &
                                       ,N_STENCIL_SFC_V                 &
                                       ,ITS,ITE,JTS,JTE                 &
                                       ,IDS,IDE,JDS,JDE                 &
!
                                       ,NTASKS_UPDATE_PARENT            &
                                       ,ID_PARENT_UPDATE_TASKS          &
                                       ,NPTS_UPDATE_PARENT              &
                                       ,I_2WAY_UPDATE                   &
                                       ,J_2WAY_UPDATE                   &
                                                          )
!
!-----------------------------------------------------------------------
!***  In 2-way mode each child domain must determine to which parent
!***  tasks and to which points on those tasks update data must be
!***  sent.  The method used here is taking the mean of the points
!***  on a stencil of child points that surround a given parent
!***  point.
!***  This routine is called from CHILDREN_SEND_PARENTS_2WAY_DATA.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument variables
!------------------------
!
      INTEGER(kind=KINT),INTENT(IN) :: I_SW_PARENT_CURRENT              &  !<-- Child domain SW corner on this parent I 
                                      ,J_SW_PARENT_CURRENT              &  !<-- Child domain SW corner on this parent J
                                      ,N_BLEND_H                        &  !<-- # of nest blending rows for H pts
                                      ,N_BLEND_V                        &  !<-- # of nest blending rows for V pts
                                      ,N_STENCIL_H                      &  !<-- Width of stencil for averaging h to parent H
                                      ,N_STENCIL_V                      &  !<-- Width of stencil for averaging v to parent V
                                      ,N_STENCIL_SFC_H                  &  !<-- Width of stencil for averaging fis,pd to parent H
                                      ,N_STENCIL_SFC_V                  &  !<-- Width of stencil for averaging fis,pd to parent V
                                      ,NUM_FCST_TASKS_PARENT            &  !<-- # of fcst tasks on this nest's parent
                                      ,SPACE_RATIO_MY_PARENT               !<-- Parent-to-child gridspace ratio
!
      INTEGER(kind=KINT),DIMENSION(0:NUM_FCST_TASKS_PARENT-1),INTENT(IN) :: &
                                                        ITS_PARENT_TASKS    &  !<-- Starting I on all parent tasks in parent space
                                                       ,ITE_PARENT_TASKS    &  !<-- Ending I on all parent tasks in parent space
                                                       ,JTS_PARENT_TASKS    &  !<-- Starting J on all parent tasks in parent space
                                                       ,JTE_PARENT_TASKS       !<-- Ending J on all parent tasks in parent space
!
      INTEGER(kind=KINT),INTENT(IN) :: IDE,IDS,ITE,ITS                  &
                                      ,JDE,JDS,JTE,JTS
!
      INTEGER(kind=KINT),INTENT(INOUT) :: NTASKS_UPDATE_PARENT             !<-- How many parent tasks does this child task update?
!
      INTEGER(kind=KINT),DIMENSION(1:4),INTENT(OUT) :: NPTS_UPDATE_PARENT  !<-- # of points to update on each parent task subdomain
!
      INTEGER(kind=KINT),DIMENSION(1:4),INTENT(OUT) :: &
                                                  ID_PARENT_UPDATE_TASKS   !<-- Local ID in P-C intracom of parent tasks to update
!
      TYPE(INTEGER_DATA),DIMENSION(1:4),INTENT(OUT) :: I_2WAY_UPDATE    &  !<-- I indices of parent points to update 
                                                      ,J_2WAY_UPDATE       !<-- J indices of parent points to update
!
!---------------------
!***  Local variables
!---------------------
!
      INTEGER(kind=KINT) :: I,ISTAT,J,KOUNT,N,N_BLEND,N_STENCIL_0       &
                           ,NPTS_PARENT_UPDATE
!
      INTEGER(kind=KINT) :: I1,I2,J1,J2
!
      INTEGER(kind=KINT),DIMENSION(1:4) :: N_STENCIL_X
!
      REAL(kind=KFPT) :: LIMIT_EAST,LIMIT_NORTH                         &
                        ,LIMIT_SOUTH,LIMIT_WEST
!
      REAL(kind=KFPT) :: MY_IDE_ON_PARENT,MY_IDS_ON_PARENT              &
                        ,MY_JDE_ON_PARENT,MY_JDS_ON_PARENT              &
                        ,MY_ITE_ON_PARENT,MY_ITS_ON_PARENT              &
                        ,MY_JTE_ON_PARENT,MY_JTS_ON_PARENT
!
      REAL(kind=KFPT) :: MY_ISTART_ON_PARENT,MY_IEND_ON_PARENT          &
                        ,MY_JSTART_ON_PARENT,MY_JEND_ON_PARENT
!
      REAL(kind=KFPT) :: RECIP_RATIO
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RECIP_RATIO=1./REAL(SPACE_RATIO_MY_PARENT)                           !<-- Reciprocal of parent-to-child gridspace ratio
!
      NTASKS_UPDATE_PARENT=0                                               !<-- Initialize the # of parent tasks this child task
!                                                                          !    will update.
!-----------------------------------------------------------------------
!***  The domain blending region must be the same for H and V points
!***  at the current time.
!-----------------------------------------------------------------------
!
      N_BLEND=N_BLEND_H
!
!-----------------------------------------------------------------------
!***  What are this child's domain limits in terms of its parent's
!***  grid?
!-----------------------------------------------------------------------
!
      MY_IDS_ON_PARENT=REAL(I_SW_PARENT_CURRENT)
      MY_IDE_ON_PARENT=REAL(I_SW_PARENT_CURRENT)+(IDE-1)*RECIP_RATIO
      MY_JDS_ON_PARENT=REAL(J_SW_PARENT_CURRENT)
      MY_JDE_ON_PARENT=REAL(J_SW_PARENT_CURRENT)+(JDE-1)*RECIP_RATIO
!
!-----------------------------------------------------------------------
!***  What are this child task's subdomain integration limits in 
!***  terms of its parent's grid?
!-----------------------------------------------------------------------
!
      MY_ITS_ON_PARENT=REAL(I_SW_PARENT_CURRENT)+(ITS-IDS)*RECIP_RATIO     !<-- Child task starting I in parent grid space
      MY_ITE_ON_PARENT=REAL(I_SW_PARENT_CURRENT)+(ITE-IDS)*RECIP_RATIO     !<-- Child task ending I in parent grid space
      MY_JTS_ON_PARENT=REAL(J_SW_PARENT_CURRENT)+(JTS-JDS)*RECIP_RATIO     !<-- Child task starting J in parent grid space
      MY_JTE_ON_PARENT=REAL(J_SW_PARENT_CURRENT)+(JTE-JDS)*RECIP_RATIO     !<-- Child task ending J in parent grid space
!
!-----------------------------------------------------------------------
!***  We want to limit the child points that can be used for 
!***  computing the 2-way data.  For now do not use any child
!***  points in the averaging stencil that lie in the child's
!***  boundary blending region.  Stencils can vary for h-->H,
!***  v-->V, fis,pd-->H, and fis,pd-->V.  Determine the set of
!***  parent target I's and J's common to all the stencils and 
!***  use that to ensure that the same parent I,J indices are 
!***  used for both H and V variables.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Loop through the four stencils.  They will be considered in
!***  this order:  h-->H, v-->V, fis,pd-->H, fis,pd-->V where small
!***  letters refer to the child and capitals refer to the parent.
!-----------------------------------------------------------------------
!
      N_STENCIL_X(1)=N_STENCIL_H
      N_STENCIL_X(2)=N_STENCIL_V
      N_STENCIL_X(3)=N_STENCIL_SFC_H
      N_STENCIL_X(4)=N_STENCIL_SFC_V
!
      MY_ISTART_ON_PARENT=MY_IDS_ON_PARENT                                 !<--
      MY_IEND_ON_PARENT  =MY_IDE_ON_PARENT                                 !   | Child domain limits in terms of
      MY_JSTART_ON_PARENT=MY_JDS_ON_PARENT                                 !   | the parent I,J.
      MY_JEND_ON_PARENT  =MY_JDE_ON_PARENT                                 !<--
!
!-----------------------------------------------------------------------
!
      DO N=1,4                                                             !<-- Loop through the four stencils.
!
!-----------------------------------------------------------------------
!
        N_STENCIL_0=N_STENCIL_X(N)/2                                       !<-- Child's delta I,J from parent update pt to 
!                                                                               west/south edge of stencil.
        LIMIT_WEST=REAL(MY_IDS_ON_PARENT)                               &
                   +(N_BLEND+N_STENCIL_0)*RECIP_RATIO
        MY_ISTART_ON_PARENT=MAX(MY_ISTART_ON_PARENT                     &  !<-- Westernmost parent I that this child task 
                               ,MY_ITS_ON_PARENT                        &  !    will update on the parent domain.
                               ,LIMIT_WEST)
!
        LIMIT_EAST=REAL(MY_IDE_ON_PARENT)                               &
                   -(N_BLEND+N_STENCIL_0)*RECIP_RATIO
        MY_IEND_ON_PARENT=MIN(MY_IEND_ON_PARENT                         &  !<-- Easternmost parent I that this child task
                             ,MY_ITE_ON_PARENT                          &  !    will update on the parent domain.
                             ,LIMIT_EAST)
!
        LIMIT_SOUTH=REAL(MY_JDS_ON_PARENT)                              &
                   +(N_BLEND+N_STENCIL_0)*RECIP_RATIO
        MY_JSTART_ON_PARENT=MAX(MY_JSTART_ON_PARENT                     &  !<-- Southernmost parent J that this child task
                               ,MY_JTS_ON_PARENT                        &  !    will update on the parent domain.
                               ,LIMIT_SOUTH)
!
        LIMIT_NORTH=REAL(MY_JDE_ON_PARENT)                              &
                   -(N_BLEND+N_STENCIL_0)*RECIP_RATIO
        MY_JEND_ON_PARENT=MIN(MY_JEND_ON_PARENT                         &  !<-- Northernmost parent J that this child task
                             ,MY_JTE_ON_PARENT                          &  !    will update on the parent domain.
                             ,LIMIT_NORTH)
!
      ENDDO
!
!-----------------------------------------------------------------------
!***  Find how many parent tasks will be updated by this child task
!***  and save their local IDs from the P-C intracommunicator.
!-----------------------------------------------------------------------
!
      find: DO N=0,NUM_FCST_TASKS_PARENT-1
!
!-----------------------------------------------------------------------
!
        IF(REAL(ITS_PARENT_TASKS(N))<MY_IEND_ON_PARENT+EPS              &
                          .AND.                                         &
           REAL(ITE_PARENT_TASKS(N))>MY_ISTART_ON_PARENT-EPS            &
                          .AND.                                         &
           REAL(JTS_PARENT_TASKS(N))<MY_JEND_ON_PARENT+EPS              &
                          .AND.                                         &
           REAL(JTE_PARENT_TASKS(N))>MY_JSTART_ON_PARENT-EPS )THEN
!
!-----------------------------------------------------------------------
!***  Now determine which points on each parent task will be updated.
!
!     Example 1: The child task's MY_ISTART_ON_PARENT is 10.666667 and
!                the parent task's ITS is 10.  Then the first parent I
!                to be updated by the child task is
!                INT(10.66667+1.-EPS)=INT(11.66667-EPS)=11.
!     Example 2: The child task's MY_ISTART_ON_PARENT is 10.999999 and
!                the parent task's ITS is 10.  Then the first parent I
!                to be updated by the child task is
!                INT(10.999999+1.-EPS)=INT(11.999999-EPS)=11.
!     Example 3: The child task's MY_ISTART_ON_PARENT is 11.000001 and
!                the parent task's ITS is 10.  Then the first parent I
!                to be updated by the child task is
!                INT(11.000001+1.-EPS)=INT(12.000001-EPS)=11.
!     Example 4: The child task's MY_IEND_ON_PARENT is 18.999999 and
!                the parent task's ITE is 20.  Then the last parent I
!                to be updated by the child task is
!                INT(18.999999+EPS)=19.
!-----------------------------------------------------------------------
!
          I1=MAX(ITS_PARENT_TASKS(N),INT(MY_ISTART_ON_PARENT+1.-EPS))      !<-- Starting parent I to update on parent task N
          I2=MIN(ITE_PARENT_TASKS(N),INT(MY_IEND_ON_PARENT+EPS))           !<-- Ending parent I to update on parent task N
          J1=MAX(JTS_PARENT_TASKS(N),INT(MY_JSTART_ON_PARENT+1.-EPS))      !<-- Starting parent J to update on parent task N
          J2=MIN(JTE_PARENT_TASKS(N),INT(MY_JEND_ON_PARENT+EPS))           !<-- Ending parent J to update on parent task N
!
          NPTS_PARENT_UPDATE=(I2-I1+1)*(J2-J1+1)                           !<-- # of points to update on parent task N
!
          IF(NPTS_PARENT_UPDATE<=0)THEN
            CYCLE find                                                     !<-- No usable 2-way exchange region on this child task.
          ENDIF
!
          NTASKS_UPDATE_PARENT=NTASKS_UPDATE_PARENT+1                      !<-- Count the # of parent tasks to update.
!
          IF(NTASKS_UPDATE_PARENT>4)THEN
            WRITE(0,11101)NTASKS_UPDATE_PARENT
11101       FORMAT(' Child task is updating ',I3,' parent tasks which is > 4')
            WRITE(0,*)' Aborting!!'
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
          ENDIF
!
          ID_PARENT_UPDATE_TASKS(NTASKS_UPDATE_PARENT)=N                   !<-- Local rank of the parent task.
!
          NPTS_UPDATE_PARENT(NTASKS_UPDATE_PARENT)=NPTS_PARENT_UPDATE      !<-- # of points to update on parent task N
!
          IF(ASSOCIATED(I_2WAY_UPDATE(NTASKS_UPDATE_PARENT)%DATA))THEN
            DEALLOCATE(I_2WAY_UPDATE(NTASKS_UPDATE_PARENT)%DATA,stat=ISTAT)
            IF(ISTAT/=0)THEN
              WRITE(0,11102)NTASKS_UPDATE_PARENT
11102         FORMAT(' Failed to deallocate I_2WAY_UPDATE(',I1,')%DATA')
              WRITE(0,*)' Aborting!!'
              CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
            ENDIF
          ENDIF
!
          ALLOCATE(I_2WAY_UPDATE(NTASKS_UPDATE_PARENT)%DATA(1:NPTS_UPDATE_PARENT(NTASKS_UPDATE_PARENT)) &
                  ,stat=ISTAT)
          IF(ISTAT/=0)THEN
            WRITE(0,11103)NTASKS_UPDATE_PARENT
11103       FORMAT(' Failed to allocate I_2WAY_UPDATE(',I1,')%DATA')
            WRITE(0,*)' Aborting!!'
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
          ENDIF
!
          IF(ASSOCIATED(J_2WAY_UPDATE(NTASKS_UPDATE_PARENT)%DATA))THEN
            DEALLOCATE(J_2WAY_UPDATE(NTASKS_UPDATE_PARENT)%DATA,stat=ISTAT)
            IF(ISTAT/=0)THEN
              WRITE(0,11104)NTASKS_UPDATE_PARENT
11104         FORMAT(' Failed to deallocate J_2WAY_UPDATE(',I1,')%DATA')
              WRITE(0,*)' Aborting!!'
              CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
            ENDIF
          ENDIF
!
          ALLOCATE(J_2WAY_UPDATE(NTASKS_UPDATE_PARENT)%DATA(1:NPTS_UPDATE_PARENT(NTASKS_UPDATE_PARENT)) &
                  ,stat=ISTAT)
          IF(ISTAT/=0)THEN
            WRITE(0,11105)NTASKS_UPDATE_PARENT
11105       FORMAT(' Failed to allocate J_2WAY_UPDATE(',I1,')%DATA')
            WRITE(0,*)' Aborting!!'
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
          ENDIF
!
!-----------------------------------------------------------------------
!***  This child task saves the parent I's and J's it will update
!***  on parent task N which is update task #NTASKS_UPDATE_PARENT).
!***  Recall that NTASKS_UPDATE_PARENT ranges from 1 to 4.
!-----------------------------------------------------------------------
!
          KOUNT=0
          DO J=J1,J2
          DO I=I1,I2
            KOUNT=KOUNT+1
            I_2WAY_UPDATE(NTASKS_UPDATE_PARENT)%DATA(KOUNT)=I
            J_2WAY_UPDATE(NTASKS_UPDATE_PARENT)%DATA(KOUNT)=J
          ENDDO
          ENDDO
!
        ENDIF
!
!-----------------------------------------------------------------------
!
      ENDDO find
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE CHILD_2WAY_BOOKKEEPING
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE GENERATE_2WAY_DATA(VAR_CHILD                           &
                                   ,PD_CHILD                            &
                                   ,FIS_CHILD                           &
                                   ,IMS,IME,JMS,JME,NVERT               &
                                   ,I_2WAY                              &
                                   ,J_2WAY                              &
                                   ,N_STENCIL                           &
                                   ,N_STENCIL_SFC                       &
                                   ,NPTS_UPDATE_PARENT                  &
                                   ,VAR_2WAY                            &
                                   ,INTERPOLATE_SFC                     &
                                   ,CHILD_SFC_ON_PARENT                 &
                                                         )
!
!-----------------------------------------------------------------------
!***  When there is 2-way nesting the children interpolate data in
!***  their domains to gridpoints in their parents' domains.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument variables
!------------------------
!
      INTEGER(kind=KINT),INTENT(IN) :: IMS,IME,JMS,JME                     !<-- Child task subdomain horizontal memory dimensions
!
      INTEGER(kind=KINT),INTENT(IN) :: N_STENCIL                        &  !<-- Use N_STENCILxN_STENCIL child pts for each parent point
                                      ,N_STENCIL_SFC                    &  !<-- Stencil width for interpolating child FIS,PD to parent
                                      ,NPTS_UPDATE_PARENT               &  !<-- # of parent points (I,J) updated on given parent task
                                      ,NVERT                               !<-- Vertical dimension of VAR_CHILD
!
      INTEGER(kind=KINT),DIMENSION(1:NPTS_UPDATE_PARENT),INTENT(IN) ::  &
                                                                I_2WAY  &  !<-- Child I on each parent update point (H or V)
                                                               ,J_2WAY     !<-- Child J on each parent update point (H or V)
!
      REAL(kind=KFPT),DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: FIS_CHILD &  !<-- The child's sfc geopotential
                                                              ,PD_CHILD     !<-- The child's PD array
!
      REAL(kind=KFPT),DIMENSION(IMS:IME,JMS:JME,1:NVERT),INTENT(IN) ::   &
                                                             VAR_CHILD     !<-- The child array of the 3-D update variable
!
      REAL(kind=KFPT),DIMENSION(1:NPTS_UPDATE_PARENT,1:2),INTENT(OUT) :: &
                                                    CHILD_SFC_ON_PARENT    !<-- Child's FIS,PD interpolated to parent update points
!
      REAL(kind=KFPT),DIMENSION(1:NPTS_UPDATE_PARENT*NVERT),INTENT(OUT) :: & 
                                                                VAR_2WAY   !<-- 2-way variable interp'd from child grid to parent's
!
      LOGICAL(kind=KLOG),INTENT(IN) :: INTERPOLATE_SFC                     !<-- Should FIS,PD be interpolated this call?
!
!---------------------
!***  Local variables
!---------------------
!
      INTEGER(kind=KINT),SAVE :: KNT_PTS
!
      INTEGER(kind=KINT) :: I,IC,J,JC,KNT_PTS_HORZ,L                    &
                           ,N_STENCIL_0,N_STENCIL_TOT,NP
!
      INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE :: I_START,I_END      &
                                                    ,J_START,J_END
!
      REAL(kind=KFPT) :: FIS_SUM,PD_SUM,RECIP_N_STENCIL_TOT,VSUM
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      N_STENCIL_0=N_STENCIL/2                                              !<-- 2->1; 3->1; 4->2; 5->2, etc.
      N_STENCIL_TOT=N_STENCIL*N_STENCIL                                    !<-- # of points in the stencil
      RECIP_N_STENCIL_TOT=1./REAL(N_STENCIL_TOT)                           !<-- Reciprocal of # of points in the stencil
!
!-----------------------------------------------------------------------
!***  Parent-child gridspace ratios can be any positive integer (>1 of
!***  course).  On the B-grid a child H point will lie on a parent H 
!***  point no matter what the ratio is.  A child V point will lie on
!***  a parent V point only for odd values of the parent-child gridspace
!***  ratio.  If the ratio is even then child H points will lie on
!***  parent V points.  The I,J of the child H point on the parent
!***  V point in that case will be the same as the child V point's
!***  immediately to the NE on the B grid.  This implies that the
!***  stencil will always be even (2x2, 4x4, etc) for interpolating
!***  to parent V points when the gridspace ratio is even while it 
!***  will be odd (3x3, 5x5, etc.) for all other cases.  The previous 
!***  statement is true for stencils that are oriented north-south.
!***  New code will need to be added if stencils rotated 45 degrees 
!***  are desired.
!
!***  The diagram below exemplifies odd and even ratios.  The capital
!***  letters are parent gridpoints and the small letters are child
!***  gridpoints.
!-----------------------------------------------------------------------
!
!           Parent-child                          Parent-child
!          gridspace ratio                       gridspace ratio
!             is odd (3)                           is even (2)
!
!
!     Hh      h       h        Hh            Hh          h         Hh
!
!          v      v       v
!                                                  v          v
!      h      h       h        h
!
!          v      Vv      v                   h          Vh        h
!
!      h      h       h        h
!                                                  v          v
!          v      v       v                       
!
!     Hh      h       h       Hh             Hh          h         Hh
!
!
!
!    Child h points lie on parent H         Child h points lie on parent H
!    points and child v points lie          points but child h points also
!    on parent V points.                    lie on parent V points. 
!
!
!-----------------------------------------------------------------------
!***  Recall that the I,J of a V point on the B grid is the same as
!***  that of the neighboring H point to the southwest.  Therefore
!***  from the diagrams above one can see that if a child point I,J
!***  coincides with a parent point to be interpolated to then
!***  the SW corner of the interpolation stencil will always be
!***  at I-N_STENCIL_0, J-N_STENCIL_0 where N_STENCIL_0 is equal to
!***  N_STENCIL/2 (integer division).
!-----------------------------------------------------------------------
!
      ALLOCATE(I_START(1:NPTS_UPDATE_PARENT))
      ALLOCATE(I_END  (1:NPTS_UPDATE_PARENT))
      ALLOCATE(J_START(1:NPTS_UPDATE_PARENT))
      ALLOCATE(J_END  (1:NPTS_UPDATE_PARENT))
!
      DO NP=1,NPTS_UPDATE_PARENT                                           !<-- Loop through this parent task subdomain's update points
!
        IC=I_2WAY(NP)                                                      !<-- Child I at parent's NP'th update point
        I_START(NP)=IC-N_STENCIL_0                                         !<-- Child I on west side of averaging stencil
        I_END(NP)  =I_START(NP)+N_STENCIL-1                                !<-- Child I on east side of averaging stencil
!
        JC=J_2WAY(NP)                                                      !<-- Child J at parent's NP'th update point
        J_START(NP)=JC-N_STENCIL_0                                         !<-- Child J on south side of averaging stencil
        J_END(NP)  =J_START(NP)+N_STENCIL-1                                !<-- Child J on north side of averaging stencil
!
      ENDDO
!
!-----------------------------------------------------------------------
!***  This child task loops through the parent points for which it is
!***  responsible on the given parent task.
!-----------------------------------------------------------------------
!
      KNT_PTS=0
      DO L=1,NVERT
!
        DO NP=1,NPTS_UPDATE_PARENT                                         !<-- Loop over update points on the given parent task
!
          VSUM=0.
!
          DO J=J_START(NP),J_END(NP)
          DO I=I_START(NP),I_END(NP)
            VSUM=VSUM+VAR_CHILD(I,J,L)                                     !<-- Sum the variable over the averaging stencil for
          ENDDO                                                            !    parent point NP.
          ENDDO
!
          KNT_PTS=KNT_PTS+1
          VAR_2WAY(KNT_PTS)=VSUM*RECIP_N_STENCIL_TOT                       !<-- Child's update value at parent point stored as 1-D
!
        ENDDO
!
      ENDDO
!
!-----------------------------------------------------------------------
!***  The child interpolates its sfc geopotential and sfc pressure
!***  to the parent points to be updated as it did for the primary
!***  prognostic variables.  If either the parent's sfc geopotential
!***  or the child's interpolated sfc geopotential is above sea level
!***  then the parent will interpolate vertically the update values
!***  received from the child to account for differences in the 
!***  domains' topographies.
!***  Note that the value of N_STENCIL_0 (the distance in I or J
!***  from the child I,J lying on the target parent H or V point to
!***  the west/south edge of the stencil) is different than above
!***  since now child H-pt values (FIS,PD) are always being averaged
!***  onto both H and V parent points.
!-----------------------------------------------------------------------
!
      IF(INTERPOLATE_SFC)THEN
!
        N_STENCIL_0=(N_STENCIL_SFC+1)/2-1                                  !<-- 2-->0; 3-->1; 4-->1; 5-->2, etc.
        N_STENCIL_TOT=N_STENCIL_SFC*N_STENCIL_SFC                          !<-- # of points in the sfc stencil
        RECIP_N_STENCIL_TOT=1./REAL(N_STENCIL_TOT)                         !<-- Reciprocal of # of points in the sfc stencil
!
        KNT_PTS_HORZ=0
!
        DO NP=1,NPTS_UPDATE_PARENT                                         !<-- Loop over update points on the given parent task
!
          PD_SUM=0.
          FIS_SUM=0.
!
          IC=I_2WAY(NP)                                                    !<-- Child I at parent's NP'th update point
          I_START(NP)=IC-N_STENCIL_0                                       !<-- Child I on west side of sfc averaging stencil
          I_END(NP)  =I_START(NP)+N_STENCIL_SFC-1                          !<-- Child I on east side of sfc averaging stencil
!
          JC=J_2WAY(NP)                                                    !<-- Child J at parent's NP'th update point
          J_START(NP)=JC-N_STENCIL_0                                       !<-- Child J on south side of sfc averaging stencil
          J_END(NP)  =J_START(NP)+N_STENCIL_SFC-1                          !<-- Child J on north side of sfc averaging stencil
!
          DO J=J_START(NP),J_END(NP)
          DO I=I_START(NP),I_END(NP)
            PD_SUM=PD_SUM+PD_CHILD(I,J)
            FIS_SUM=FIS_SUM+FIS_CHILD(I,J)
          ENDDO
          ENDDO
!
          KNT_PTS_HORZ=KNT_PTS_HORZ+1
          CHILD_SFC_ON_PARENT(KNT_PTS_HORZ,1)=FIS_SUM*RECIP_N_STENCIL_TOT  !<-- Child's mean sfc geopotential within stencil
          CHILD_SFC_ON_PARENT(KNT_PTS_HORZ,2)=PD_SUM*RECIP_N_STENCIL_TOT   !<-- Child's mean PD within stencil
!
        ENDDO
!
      ENDIF
!
!-----------------------------------------------------------------------
!
      DEALLOCATE(I_START)
      DEALLOCATE(I_END)
      DEALLOCATE(J_START)
      DEALLOCATE(J_END)
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE GENERATE_2WAY_DATA
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
      SUBROUTINE SPLINE(NOLD,XOLD,YOLD,Y2,Y2_K,NNEW,XNEW,YNEW)
!-----------------------------------------------------------------------
!
!     ******************************************************************
!     *                                                                *
!     *  This is a one-dimensional cubic spline fitting routine        *
!     *  programmed for a small scalar machine.                        *
!     *                                                                *
!     *  Programmer: Z. Janjic, Yugoslav Fed. Hydromet. Inst., Beograd *
!     *                                                                *
!     *  NOLD - Number of given values of the function.  Must be >= 3. *
!     *  XOLD - Locations of the points at which the values of the     *
!     *         function are given.  Must be in ascending order.       *
!     *  YOLD - The given values of the function at the points XOLD.   *
!     *  Y2   - The second derivatives at the points XOLD.  If natural *
!     *         spline is fitted Y2(1)=0 and Y2(nold)=0. Must be       *
!     *         specified.                                             *
!     *  Y2_K - Vertical dimension of Y2 array.                        *
!     *  NNEW - Number of values of the function to be calculated.     *
!     *  XNEW - Locations of the points at which the values of the     *
!     *         function are calculated.  XNEW(K) must be >= XOLD(1)   *
!     *         and <= XOLD(NOLD).                                     *
!     *  YNEW - The values of the function to be calculated.           *
!     *                                                                *
!     ******************************************************************
!
!-----------------------------------------------------------------------
!***  Arguments
!-----------------------------------------------------------------------
!
      INTEGER,INTENT(IN) :: NNEW,NOLD,Y2_K
!
      REAL,DIMENSION(1:NOLD),INTENT(IN) :: XOLD,YOLD
      REAL,DIMENSION(1:NNEW),INTENT(IN) :: XNEW
!
      REAL,DIMENSION(1:Y2_K),INTENT(INOUT) :: Y2
!
      REAL,DIMENSION(1:NNEW),INTENT(OUT) :: YNEW
!
!-----------------------------------------------------------------------
!***  Local Variables
!-----------------------------------------------------------------------
!
      INTEGER :: K,K1,K2,KOLD,NOLDM1
!
      REAL :: AK,BK,CK,DEN,DX,DXC,DXL,DXR,DYDXL,DYDXR,RDX,RTDXC         &
             ,X,XK,XSQ,Y2K,Y2KP1
!
      REAL,DIMENSION(1:NOLD-2) :: P,Q
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      NOLDM1=NOLD-1
!
      DXL=XOLD(2)-XOLD(1)
      DXR=XOLD(3)-XOLD(2)
      DYDXL=(YOLD(2)-YOLD(1))/DXL
      DYDXR=(YOLD(3)-YOLD(2))/DXR
      RTDXC=0.5/(DXL+DXR)
!
      P(1)= RTDXC*(6.*(DYDXR-DYDXL)-DXL*Y2(1))
      Q(1)=-RTDXC*DXR
!
      IF(NOLD==3) GO TO 700
!
!-----------------------------------------------------------------------
      K=3
!
  100 CONTINUE
      DXL=DXR
      DYDXL=DYDXR
      DXR=XOLD(K+1)-XOLD(K)
      DYDXR=(YOLD(K+1)-YOLD(K))/DXR
      DXC=DXL+DXR
      DEN=1./(DXL*Q(K-2)+DXC+DXC)
!
      P(K-1)= DEN*(6.*(DYDXR-DYDXL)-DXL*P(K-2))
      Q(K-1)=-DEN*DXR
!
      K=K+1
      IF(K<NOLD) GO TO 100
!
!-----------------------------------------------------------------------
!
  700 CONTINUE
      K=NOLDM1
!
  200 CONTINUE
      Y2(K)=P(K-1)+Q(K-1)*Y2(K+1)
!
      K=K-1
      IF(K>1) GO TO 200
!
!-----------------------------------------------------------------------
!
      K1=1
!
  300 CONTINUE
      XK=XNEW(K1)
!
      DO 400 K2=2,NOLD
        IF(XOLD(K2)<=XK) GO TO 400
        KOLD=K2-1
        GO TO 450
  400 CONTINUE
!
      YNEW(K1)=YOLD(NOLD)
      GO TO 600
!
  450 CONTINUE
      IF(K1==1)   GO TO 500
      IF(K==KOLD) GO TO 550
!
  500 CONTINUE
      K=KOLD
!
      Y2K=Y2(K)
      Y2KP1=Y2(K+1)
      DX=XOLD(K+1)-XOLD(K)
      RDX=1./DX
!
      AK=0.1666667*RDX*(Y2KP1-Y2K)
      BK=0.5*Y2K
      CK=RDX*(YOLD(K+1)-YOLD(K))-0.1666667*DX*(Y2KP1+Y2K+Y2K)
!
  550 CONTINUE
      X=XK-XOLD(K)
      XSQ=X*X
!
      YNEW(K1)=AK*XSQ*X+BK*XSQ+CK*X+YOLD(K)
!
  600 CONTINUE
      K1=K1+1
!
      IF(K1<=NNEW) GO TO 300
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE SPLINE
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE HYPERBOLA(A)
!
!-----------------------------------------------------------------------
!***  Generate a hyperbola that will reduce the magnitude of the
!***  source domain's underground extrapolation in those instances
!***  when the target domain's ground surface lies below the source
!***  domain's.  The hyperbola has the formula:
!
!     Y=A/(X+A)
!
!***  The value of Y is the fraction between 1 and 0 that provides
!***  the reduction in the amount added to the source domain's lowest
!***  layer value to account for the extrapolation underground.
!***  The value of X is the difference in pressure (Pa) between the
!***  source domain's lowest pressure level and the target pressure
!***  of the extrapolation.  When the pressure difference is zero then
!***  there is no reduction in the source domain's extrapolation and
!***  so the value of Y is 1.0.  For very large extrapolations then
!***  the amount added to the source domain's lowest layer value to
!***  account for the extrapolation is reduced by a factor approaching
!***  zero.
!***  The formula gives the user 1 degree of freedom.  Specify one
!***  extrapolated underground pressure depth and the amount desired
!***  for the reduction in the linear extrapolation of the source
!***  domain's lowest layer value through that depth.
!***  For example, if X1=10000.0 and Y1=0.05 then when the lowest
!***  layer value in the source domain is linearly extrapolated
!***  through an underground depth of 10000 Pa then the amount added
!***  to that lowest layer value to account for the extrapolation is
!***  first multiplied by 0.05.
!-----------------------------------------------------------------------
!
      REAL(kind=KDBL),PARAMETER :: X1=10000.0, Y1=0.05
!
!------------------------
!***  Argument Variables
!------------------------
!
      REAL(kind=KDBL),INTENT(OUT) :: A                                     !<-- Constant in the hyperbola Y=A/(X+A)
!
!---------------------
!***  Local Variables
!---------------------
!
      REAL(kind=KDBL) :: F,G,H,DISCRIM,PROD1,PROD2
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      A=(X1*Y1)/(1.-Y1)
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE HYPERBOLA
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE CHECK_REAL(P_IN,NAME)
!
!-----------------------------------------------------------------------
!***  Check the status of pointer P_IN and deallocate or nullify.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      REAL(kind=KFPT),DIMENSION(:),POINTER,INTENT(INOUT) :: P_IN
!
      CHARACTER(len=*),INTENT(IN) :: NAME
!
!--------------------
!*** Local Variables
!--------------------
!
      INTEGER(kind=KINT) :: ISTAT
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      IF(ASSOCIATED(P_IN))THEN
        DEALLOCATE(P_IN,stat=ISTAT)
        IF(ISTAT/=0)THEN
          NULLIFY(P_IN)
          WRITE(0,*)NAME,' was associated but not allocated. '          &
                        ,' It has now been nullified.'
        ELSE
          WRITE(0,*)' Forced to deallocate ',NAME
        ENDIF
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE CHECK_REAL
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE CHECK(RC)
!
      IMPLICIT NONE
!
      INTEGER,INTENT(IN) :: RC
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      IF(RC/=NF90_NOERR)THEN
        WRITE(*,*)TRIM(ADJUSTL(NF90_STRERROR(RC)))
!       WRITE(0,11101)RC
11101   FORMAT(' ERROR: RC=',I5)
      ENDIF
!
      END SUBROUTINE CHECK
!
!-----------------------------------------------------------------------
!
      END MODULE MODULE_NESTING
!
!-----------------------------------------------------------------------
