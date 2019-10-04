!#include "../../ESMFVersionDefine.h"

!-----------------------------------------------------------------------
      MODULE MODULE_WRITE_INTERNAL_STATE_GFS
!-----------------------------------------------------------------------
!***  THE INTERNAL STATE OF THE WRITE COMPONENT.
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
!       03 Sep 2009:  W. Yang  - Ensemble GEFS.
!       29 Sep 2010:  J. Wang  - change all_data_I/R2d from 1D to 2D 
!       28 Sep 2011:  W. Yang  - Modified for using the ESMF 5.2.0r library.
!
!-----------------------------------------------------------------------
!
      USE ESMF
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      INTEGER,PARAMETER :: MAX_DATA_I1D=50             !<-- Max # of 1D integer arrays
      INTEGER,PARAMETER :: MAX_DATA_I2D=50             !<-- Max # of 2D integer arrays
      INTEGER,PARAMETER :: MAX_DATA_R1D=50             !<-- Max # of 1D real arrays
      INTEGER,PARAMETER :: MAX_DATA_R2D=2100           !<-- Max # of 2D real arrays and layers of all real 3D arrays combined
      INTEGER,PARAMETER :: MAX_DATA_LOG=10
      INTEGER,PARAMETER :: NAME_MAXSTR=32
!
!-----------------------------------------------------------------------

      TYPE WRITE_INTERNAL_STATE_GFS

!--------------------------------
! PE INFORMATION AND TASK LAYOUT
!--------------------------------
!
      INTEGER :: MYPE
      INTEGER :: NUM_PES_FCST
      INTEGER :: IHALO,JHALO
      INTEGER :: NTASKS
      INTEGER :: WRITE_GROUPS,WRITE_TASKS_PER_GROUP
!
!-----------------------------
!***  FULL DOMAIN INFORMATION
!-----------------------------
!
      INTEGER,DIMENSION(:),ALLOCATABLE :: IM
      INTEGER,DIMENSION(:),ALLOCATABLE :: JM
      INTEGER,DIMENSION(:),ALLOCATABLE :: LM
      TYPE(ESMF_Logical)           :: GLOBAL
!
!--------------------
!***  SUBDOMAIN SIZE
!--------------------
!--- GFS
      INTEGER                          :: ipt_lats_node_a               &
                                         ,lats_node_a
      INTEGER,DIMENSION(:),ALLOCATABLE :: global_lats_a
      INTEGER,DIMENSION(:),ALLOCATABLE :: fcst_lat_to_write_task        &
                                         ,nwrttask_on_fcst              &
                                         ,nlat_to_write_task            &
                                         ,nstart_to_write_task          &
                                         ,nlat_from_fcst_task           &
                                         ,nstart_from_fcst_task         &
                                         ,fcst_lat_on_write_task
      INTEGER,DIMENSION(:),ALLOCATABLE :: JSTART_WRITE                  &
                                         ,JEND_WRITE
!
!----------------------------------------------------
!***  IDs OF FCST TASKS THAT SEND TO EACH WRITE TASK
!----------------------------------------------------
!
      INTEGER,DIMENSION(:),ALLOCATABLE :: ID_FTASK_RECV_STA             &
                                         ,ID_FTASK_RECV_END
!
!------------------------------
!***  BUNDLE DATA INFORMATION
!------------------------------
!
      INTEGER,DIMENSION(:),ALLOCATABLE :: KOUNT_I1D                     &
                                         ,KOUNT_I2D                     &
                                         ,KOUNT_R1D                     &
                                         ,KOUNT_R2D                     &
                                         ,KOUNT_LOG
!
      INTEGER,DIMENSION(:,:),ALLOCATABLE :: LENGTH_DATA_I1D             &
                                           ,LENGTH_DATA_R1D             &
                                           ,LENGTH_DATA_R2D
!
      INTEGER,DIMENSION(:),ALLOCATABLE :: LENGTH_SUM_I1D                &
                                         ,LENGTH_SUM_R1D                &
                                         ,LENGTH_SUM_I2D                &
                                         ,LENGTH_SUM_R2D                &
                                         ,LENGTH_SUM_LOG
!
      INTEGER,DIMENSION(:),ALLOCATABLE :: NCOUNT_FIELDS
!
      INTEGER,DIMENSION(:,:)  ,POINTER :: ALL_DATA_I1D => null()
      INTEGER,DIMENSION(:,:)  ,POINTER :: ALL_DATA_I2D => null()
      INTEGER,DIMENSION(:,:),ALLOCATABLE :: OUTPUT_ARRAY_I2D
!
      REAL   ,DIMENSION(:,:)  ,POINTER :: ALL_DATA_R1D => null()
      REAL(kind=4)   ,DIMENSION(:,:)  ,POINTER :: ALL_DATA_R2D => null()
      REAL(kind=4)   ,DIMENSION(:,:),ALLOCATABLE :: OUTPUT_ARRAY_R2D
!
!-----------------------------------------------------------------------
!*** STORAGE ARRAYS
!-----------------------------------------------------------------------
!
      INTEGER,DIMENSION(:)    ,ALLOCATABLE :: BUFF_INT,BUFF_INT_TMP
      INTEGER,DIMENSION(:,:,:),ALLOCATABLE :: WRITE_SUBSET_I
      REAL(KIND=4),DIMENSION(:)    ,ALLOCATABLE :: BUFF_REAL,BUFF_REAL_TMP
      REAL(KIND=4),DIMENSION(:,:,:),ALLOCATABLE :: WRITE_SUBSET_R
!
      TYPE(ESMF_Logical),DIMENSION(:,:),POINTER :: ALL_DATA_LOG  => null()
!
      CHARACTER(NAME_MAXSTR),DIMENSION(:,:),ALLOCATABLE :: FIELD_NAME
!
      CHARACTER(ESMF_MAXSTR*MAX_DATA_I1D),DIMENSION(:),ALLOCATABLE :: NAMES_I1D_STRING
      CHARACTER(ESMF_MAXSTR*MAX_DATA_I2D),DIMENSION(:),ALLOCATABLE :: NAMES_I2D_STRING
      CHARACTER(ESMF_MAXSTR*MAX_DATA_R1D),DIMENSION(:),ALLOCATABLE :: NAMES_R1D_STRING
      CHARACTER(ESMF_MAXSTR*MAX_DATA_R2D),DIMENSION(:),ALLOCATABLE :: NAMES_R2D_STRING
      CHARACTER(ESMF_MAXSTR*MAX_DATA_LOG),DIMENSION(:),ALLOCATABLE :: NAMES_LOG_STRING
!
!-----------------------------------------------------------------------
!***  THE OUTPUT FILE
!-----------------------------------------------------------------------
!
      CHARACTER(3)           :: CORE
!
      INTEGER                :: IO_UNIT
      INTEGER                :: IO_RECL
      INTEGER                :: NFHOUR
      REAL(KIND=4)           :: ZHOUR
      REAL(KIND=4)           :: PDRYINI
!
!jws
      INTEGER                :: NUM_FILE
      CHARACTER(ESMF_MAXSTR),DIMENSION(:),ALLOCATABLE :: FILENAME_BASE
      CHARACTER(ESMF_MAXSTR),DIMENSION(:),ALLOCATABLE :: IO_FORM
      CHARACTER(ESMF_MAXSTR),DIMENSION(:),ALLOCATABLE :: IO_FILE
!jwe
      CHARACTER(ESMF_MAXSTR) :: IO_STATUS
      CHARACTER(ESMF_MAXSTR) :: IO_ACCESS
      CHARACTER(ESMF_MAXSTR) :: IO_POSITION
      CHARACTER(ESMF_MAXSTR) :: IO_ACTION
      CHARACTER(ESMF_MAXSTR) :: IO_DELIM
      CHARACTER(ESMF_MAXSTR) :: IO_PAD
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
      LOGICAL :: QUILTING
      LOGICAL :: WRITE_NEMSIOFLAG
      LOGICAL :: WRITE_NEMSIOCTL
!
!-----------------------------------------
!***  POST flags 
!-----------------------------------------
!
      LOGICAL                :: WRITE_DOPOST
      LOGICAL                :: IAU
      CHARACTER(ESMF_MAXSTR) :: POST_GRIBVERSION
      LOGICAL                :: GOCART_AER2POST
      integer                :: nlunit             ! post namelist unit number - Moorthi
      character(80)          :: post_namelist
 
!-----------------------------------------------------------------------
!
      END TYPE WRITE_INTERNAL_STATE_GFS
!
!-----------------------------------------------------------------------
!***  THIS STATE IS SUPPORTED BY C POINTERS BUT NOT F90 POINTERS
!***  THEREFORE WE NEED THIS WRAP.
!-----------------------------------------------------------
!
      TYPE WRITE_WRAP_GFS
        TYPE(WRITE_INTERNAL_STATE_GFS),POINTER :: WRITE_INT_STATE
      END TYPE WRITE_WRAP_GFS

!-----------------------------------------------------------
!
      END MODULE MODULE_WRITE_INTERNAL_STATE_GFS
!
!-----------------------------------------------------------
