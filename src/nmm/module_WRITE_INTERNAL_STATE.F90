!-----------------------------------------------------------------------
!
      MODULE MODULE_WRITE_INTERNAL_STATE
!
!-----------------------------------------------------------------------
!***  The internal state of the Write component.
!-----------------------------------------------------------------------
!***
!***  HISTORY
!***
!       xx Feb 2007:  W. Yang - Originator
!       14 Jun 2007:  T. Black - Name revisions
!       14 Aug 2007:  T. Black - Some pointers changed to arrays
!       11 Sep 2007:  T. Black - Updates for quilting
!       15 Aug 2008:  J. Wang  - Add NEMSIO variables
!       16 Sep 2008:  J. Wang  - 3-D output arrays revert to 2-D
!       04 Sep 2009:  T. Black - Add the 1-D boundary restart arrays
!       22 Apr 2010:  T. Black - Add minutes and seconds to elapsed
!                                forecast time.
!       15 Sep 2010:  T. Black - Changed many components to pointers
!          Feb 2011:  W. Yang  - Updated to use both the ESMF 4.0.0rp2 library,
!                                ESMF 5 library and the the ESMF 3.1.0rp2 library.
!          May 2011:  J. Wang  - add dopost and post_gribversion option to run post
!       27 SEP 2011   W. Yang, - Modified for using the ESMF 5.2.0r library.
!
!---------------------------------------------------------------------------------
!
      USE ESMF
      USE MODULE_KINDS
      USE MODULE_DERIVED_TYPES,ONLY: BC_H_ALL,BC_V_ALL
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      PRIVATE
!
      PUBLIC :: WRITE_INTERNAL_STATE,WRITE_WRAP                         &
               ,MAX_DATA_I1D,MAX_DATA_I2D                               &
               ,MAX_DATA_R1D,MAX_DATA_R2D                               &
               ,MAX_DATA_LOG
!
!-----------------------------------------------------------------------
!
      INTEGER(kind=KINT),PARAMETER :: MAX_DATA_I1D=50                      !<-- Max # of 1D integer arrays
      INTEGER(kind=KINT),PARAMETER :: MAX_DATA_I2D=50                      !<-- Max # of 2D integer arrays
      INTEGER(kind=KINT),PARAMETER :: MAX_DATA_R1D=50                      !<-- Max # of 1D real arrays
      INTEGER(kind=KINT),PARAMETER :: MAX_DATA_R2D=10000                   !<-- Max # of 2D real arrays and layers
                                                                           !    of all real 3D arrays combined
      INTEGER(kind=KINT),PARAMETER :: MAX_DATA_LOG=10
!
!-----------------------------------------------------------------------

      TYPE WRITE_INTERNAL_STATE

!------------------------------------
!***  PE information and task layout
!------------------------------------
!
      INTEGER(kind=KINT) :: MYPE
      INTEGER(kind=KINT) :: INPES,JNPES
      INTEGER(kind=KINT) :: IHALO,JHALO
      INTEGER(kind=KINT) :: LAST_FCST_TASK,LAST_WRITE_TASK
      INTEGER(kind=KINT) :: LEAD_WRITE_TASK
      INTEGER(kind=KINT) :: NTASKS
      INTEGER(kind=KINT) :: WRITE_GROUPS,WRITE_TASKS_PER_GROUP
!
!-----------------------------
!***  Full domain information
!-----------------------------
!
      INTEGER(kind=KINT) :: ID_DOMAIN
!
      INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE :: IM
      INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE :: IDS
      INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE :: IDE
      INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE :: JM
      INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE :: JDS
      INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE :: JDE
      INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE :: LM
!
      LOGICAL(kind=KLOG) :: GLOBAL
!
!----------------------------------------------------
!***  The forecast or quilt tasks' intracommunicator
!----------------------------------------------------
!
      INTEGER(kind=KINT) :: MPI_COMM_COMP
!
!---------------------------------------------
!***  Array of Write group intercommunicators
!---------------------------------------------
!
      INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE :: MPI_INTERCOMM_ARRAY
!
!--------------------
!***  Subdomain size
!--------------------
!
      INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE :: LOCAL_ISTART       &
                                                    ,LOCAL_IEND         &
                                                    ,LOCAL_JSTART       &
                                                    ,LOCAL_JEND
!
!----------------------------------------------------
!***  IDs of fcst tasks that send to each write task
!----------------------------------------------------
!
      INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE :: ID_FTASK_RECV_STA  &
                                                    ,ID_FTASK_RECV_END
!
!----------------------------------------------------
!***  # of words sent by each forecast task 
!***  to its designated write task.
!----------------------------------------------------
!
      INTEGER(kind=KINT),POINTER :: NUM_WORDS_SEND_I2D_HST              &
                                   ,NUM_WORDS_SEND_R2D_HST              &
                                   ,NUM_WORDS_SEND_I2D_RST              & 
                                   ,NUM_WORDS_SEND_R2D_RST
!
!----------------------------------------------------
!***  # of words received by each write task from
!***  all of its designated forecast tasks.
!----------------------------------------------------
!
      INTEGER(kind=KINT),DIMENSION(:),POINTER :: NUM_WORDS_RECV_I2D_HST =>NULL() &
                                                ,NUM_WORDS_RECV_R2D_HST =>NULL() &
                                                ,NUM_WORDS_RECV_I2D_RST =>NULL() &
                                                ,NUM_WORDS_RECV_R2D_RST =>NULL()
!
!--------------------------------------
!***  History/restart data information
!--------------------------------------
!
      INTEGER(kind=KINT),DIMENSION(:),POINTER :: KOUNT_I1D =>NULL()              &
                                                ,KOUNT_I2D =>NULL()              &
                                                ,KOUNT_R1D =>NULL()              &
                                                ,KOUNT_R2D =>NULL()              &
                                                ,KOUNT_LOG =>NULL()
!
      INTEGER(kind=KINT),DIMENSION(:),POINTER :: RST_KOUNT_I1D =>NULL()          &
                                                ,RST_KOUNT_I2D =>NULL()          &
                                                ,RST_KOUNT_R1D =>NULL()          &
                                                ,RST_KOUNT_R2D =>NULL()          &
                                                ,RST_KOUNT_LOG =>NULL()
!
      INTEGER(kind=KINT),DIMENSION(:),POINTER :: LENGTH_DATA_I1D =>NULL()        &
                                                ,LENGTH_DATA_R1D =>NULL()        &
                                                ,LENGTH_DATA_R2D =>NULL()
!
      INTEGER(kind=KINT),DIMENSION(:),POINTER :: RST_LENGTH_DATA_I1D =>NULL()    &
                                                ,RST_LENGTH_DATA_R1D =>NULL()    &
                                                ,RST_LENGTH_DATA_R2D =>NULL()
!
      INTEGER(kind=KINT),DIMENSION(:),POINTER :: LENGTH_SUM_I1D =>NULL()         &
                                                ,LENGTH_SUM_R1D =>NULL()         &
                                                ,LENGTH_SUM_R2D =>NULL()         &
                                                ,LENGTH_SUM_LOG =>NULL()
!
      INTEGER(kind=KINT),DIMENSION(:),POINTER :: RST_LENGTH_SUM_I1D =>NULL()     &
                                                ,RST_LENGTH_SUM_R1D =>NULL()     &
                                                ,RST_LENGTH_SUM_R2D =>NULL()     &
                                                ,RST_LENGTH_SUM_LOG =>NULL()
!
      INTEGER(kind=KINT),DIMENSION(:),POINTER :: NCOUNT_FIELDS =>NULL()          &
                                                ,RST_NCOUNT_FIELDS =>NULL()
!
      INTEGER(kind=KINT),DIMENSION(:),POINTER :: ALL_DATA_I1D =>NULL()           &
                                                ,ALL_DATA_I2D =>NULL()
!
      INTEGER(kind=KINT),DIMENSION(:,:),POINTER :: OUTPUT_ARRAY_I2D =>NULL()
!
      INTEGER(kind=KINT),DIMENSION(:),POINTER :: RST_ALL_DATA_I1D =>NULL()       &
                                                ,RST_ALL_DATA_I2D =>NULL()
!
      INTEGER(kind=KINT),DIMENSION(:,:),POINTER :: RST_OUTPUT_ARRAY_I2D =>NULL()
!
      REAL(kind=KFPT),DIMENSION(:),POINTER :: ALL_DATA_R1D =>NULL()              &
                                             ,ALL_DATA_R2D =>NULL()
!
      REAL(kind=KFPT),DIMENSION(:,:),POINTER :: OUTPUT_ARRAY_R2D =>NULL()
!
      REAL(kind=KFPT),DIMENSION(:),POINTER :: RST_ALL_DATA_R1D =>NULL()          &
                                             ,RST_ALL_DATA_R2D =>NULL()
!
      REAL(kind=KFPT),DIMENSION(:,:),POINTER :: RST_OUTPUT_ARRAY_R2D =>NULL()
!
!----------------------
!***  Boundary restart 
!----------------------
!
      INTEGER(kind=KINT) :: LNSH,LNSV                                      !<-- # of H,V bndry rows, respectively
      INTEGER(kind=KINT) :: NLEV_H,NLEV_V                                  !<-- Total # of 2-D levels in all H-pt,V-pt bndry vbls
!
!-----------
!***  Local
!-----------
!
      INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE :: NUM_WORDS_BC_SOUTH &  !<-- Word counts of 1-D boundary data strings
                                                    ,NUM_WORDS_BC_NORTH &  !    for each side of the domain.
                                                    ,NUM_WORDS_BC_WEST  &  !
                                                    ,NUM_WORDS_BC_EAST                !<--
!
      REAL(kind=KFPT),DIMENSION(:),ALLOCATABLE :: RST_BC_DATA_SOUTH     &  !<-- 1-D strings of boundary data 
                                                 ,RST_BC_DATA_NORTH     &  !    for each side of the domain.
                                                 ,RST_BC_DATA_WEST      &  !
                                                 ,RST_BC_DATA_EAST         !<--
!
!-----------------
!***  Full-domain
!-----------------
!
      INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE :: NUM_WORDS_SEND_BC     !<-- Word count of full-domain 1-D boundary data string
!
      INTEGER(kind=KINT) :: NVARS_BC_2D_H                               &  !<-- # of 2-D H-pt boundary variables
                           ,NVARS_BC_3D_H                               &  !<-- # of 3-D H-pt boundary variables
                           ,NVARS_BC_4D_H                               &  !<-- # of 4-D H-pt boundary variables
                           ,NVARS_BC_2D_V                               &  !<-- # of 2-D V-pt boundary variables
                           ,NVARS_BC_3D_V                                     !<-- # of 3-D V-pt boundary variables
!
      INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE :: LBND_4D            &  !<-- Lower,upper bounds of the counts of the # of 
                                                    ,UBND_4D               !    3-D arrays in each 4-D boundary variable
!
      REAL(kind=KFPT),DIMENSION(:),ALLOCATABLE :: RST_ALL_BC_DATA          !<-- 1-D string of full-domain boundary data
!
      TYPE(BC_H_ALL) :: BND_VARS_H                                         !<-- All H-pt boundary data/tendencies
!
      TYPE(BC_V_ALL) :: BND_VARS_V                                         !<-- All V-pt boundary data/tendencies
!
!--------------------
!***  Storage arrays
!--------------------
!
      INTEGER(kind=KINT),DIMENSION(:),POINTER :: BUFF_INT =>NULL()              &
                                                ,RST_BUFF_INT =>NULL()
!
      INTEGER(kind=KINT),DIMENSION(:,:,:),POINTER :: WRITE_SUBSET_I =>NULL()    &
                                                    ,RST_WRITE_SUBSET_I=>NULL()
!
      REAL(kind=KFPT),DIMENSION(:),POINTER :: BUFF_REAL =>NULL()                &
                                             ,RST_BUFF_REAL=>NULL()
!
      REAL(kind=KFPT),DIMENSION(:,:,:),POINTER :: WRITE_SUBSET_R =>NULL()       &
                                                 ,RST_WRITE_SUBSET_R=>NULL()
!
      LOGICAL,           DIMENSION(:),POINTER :: ALL_DATA_LOG =>NULL()     
      LOGICAL,           DIMENSION(:),POINTER :: RST_ALL_DATA_LOG=>NULL()
!
      CHARACTER(ESMF_MAXSTR),DIMENSION(:),POINTER :: FIELD_NAME=>NULL()         &
                                                    ,RST_FIELD_NAME=>NULL()
!
      CHARACTER(ESMF_MAXSTR*MAX_DATA_I1D),POINTER :: NAMES_I1D_STRING
      CHARACTER(ESMF_MAXSTR*MAX_DATA_I2D),POINTER :: NAMES_I2D_STRING
      CHARACTER(ESMF_MAXSTR*MAX_DATA_R1D),POINTER :: NAMES_R1D_STRING
      CHARACTER(ESMF_MAXSTR*MAX_DATA_R2D),POINTER :: NAMES_R2D_STRING
      CHARACTER(ESMF_MAXSTR*MAX_DATA_LOG),POINTER :: NAMES_LOG_STRING
!
      CHARACTER(ESMF_MAXSTR*MAX_DATA_I1D),POINTER :: RST_NAMES_I1D_STRING
      CHARACTER(ESMF_MAXSTR*MAX_DATA_I2D),POINTER :: RST_NAMES_I2D_STRING
      CHARACTER(ESMF_MAXSTR*MAX_DATA_R1D),POINTER :: RST_NAMES_R1D_STRING
      CHARACTER(ESMF_MAXSTR*MAX_DATA_R2D),POINTER :: RST_NAMES_R2D_STRING
      CHARACTER(ESMF_MAXSTR*MAX_DATA_LOG),POINTER :: RST_NAMES_LOG_STRING
!
!---------------------
!***  The output file
!---------------------
!
      INTEGER(kind=KINT) :: IO_HST_UNIT,IO_RST_UNIT
      INTEGER(kind=KINT) :: IO_RECL
      INTEGER(kind=KINT) :: NFHOURS
      INTEGER(kind=KINT) :: NFMINUTES
!
      REAL(kind=KFPT)    :: NFSECONDS
!
      CHARACTER(ESMF_MAXSTR) :: HST_NAME_BASE,RST_NAME_BASE
!
      CHARACTER(ESMF_MAXSTR) :: POST_GRIBVERSION
!
!-------------------------------------
!***  Times used in history filenames
!-------------------------------------
!
      TYPE(ESMF_Time)         :: IO_BASETIME
      TYPE(ESMF_TimeInterval) :: IO_CURRTIMEDIFF
!
!-----------------------------------------
!***  I/O direction flags (Read or Write)
!-----------------------------------------
!
      LOGICAL(kind=KLOG) :: WRITE_HST_BIN,WRITE_HST_NEMSIO
      LOGICAL(kind=KLOG) :: WRITE_RST_BIN,WRITE_RST_NEMSIO
      LOGICAL(kind=KLOG) :: WRITE_NEMSIOCTL
      LOGICAL(kind=KLOG) :: WRITE_FSYNCFLAG
      LOGICAL(kind=KLOG) :: WRITE_DONEFILEFLAG
      LOGICAL(kind=KLOG) :: WRITE_DOPOST
      LOGICAL(kind=KLOG) :: PRINT_ALL
      LOGICAL(kind=KLOG) :: PRINT_DIAG
      LOGICAL(kind=KLOG) :: PRINT_OUTPUT
      LOGICAL(kind=KLOG) :: PRINT_ESMF

      integer            :: nlunit        ! post namelist unit number - Moorthi
      character(80)      :: post_namelist
 
!-----------------------------------------------------------------------
!
      END TYPE WRITE_INTERNAL_STATE
!
!-----------------------------------------------------------------------
!***  This state is supported by C pointers but not F90 pointers
!***  therefore we need this wrap.
!-----------------------------------------------------------------------
!
      TYPE WRITE_WRAP
        TYPE(WRITE_INTERNAL_STATE),POINTER :: WRITE_INT_STATE
      END TYPE WRITE_WRAP

!-----------------------------------------------------------
!
      END MODULE MODULE_WRITE_INTERNAL_STATE
!
!-----------------------------------------------------------
