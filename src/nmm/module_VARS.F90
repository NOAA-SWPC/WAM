!-----------------------------------------------------------------------

      MODULE MODULE_VARS

!-----------------------------------------------------------------------
!***  This module contains the routine that reads in the text files
!***  in which the user has specified internal state variables from
!***  the Solver component to be:
!***    (1) In the history output
!***    (2) In the restart output
!***    (3) 'Owned' (allocated) by the component
!***    (4) In the component's import state
!***    (5) In the component's export state
!
!***  and contains the routines that allocate memory within the
!***  general VARS composite variable and then point the 'Owned'
!***  internal state variables into that memory.
!-----------------------------------------------------------------------

      USE MODULE_KINDS

!-----------------------------------------------------------------------
      IMPLICIT NONE

      PRIVATE
      PUBLIC :: VAR
      PUBLIC :: TKR_I0D,TKR_I1D,TKR_I2D,TKR_R0D,TKR_R1D,TKR_R2D,TKR_R3D,TKR_R4D
      PUBLIC :: READ_CONFIG
      PUBLIC :: SET_VAR_PTR
      PUBLIC :: DEALLOC_VARS
      PUBLIC :: FIND_VAR_INDX

      TYPE VAR

        CHARACTER(LEN=32)                  :: VBL_NAME = ''
        LOGICAL                            :: HISTORY  = .FALSE.
        LOGICAL                            :: HISTORY_P= .FALSE.
        LOGICAL                            :: RESTART  = .FALSE.
        LOGICAL                            :: RESTART_P= .FALSE.
        LOGICAL                            :: OWNED    = .FALSE.
        LOGICAL                            :: IMPORT   = .FALSE.
        LOGICAL                            :: EXPORT   = .FALSE.
        LOGICAL                            :: TSERIES  = .FALSE.
        LOGICAL                            :: TSERIES_P= .FALSE.
        CHARACTER(LEN=128)                 :: DESCRIPTION = ''

        INTEGER                            :: TKR      = 0

        INTEGER                   ,POINTER :: I0D => NULL()
        INTEGER,DIMENSION(:)      ,POINTER :: I1D => NULL()
        INTEGER,DIMENSION(:,:)    ,POINTER :: I2D => NULL()
        REAL                      ,POINTER :: R0D => NULL()
        REAL   ,DIMENSION(:)      ,POINTER :: R1D => NULL()
        REAL   ,DIMENSION(:,:)    ,POINTER :: R2D => NULL()
        REAL   ,DIMENSION(:,:,:)  ,POINTER :: R3D => NULL()
        REAL   ,DIMENSION(:,:,:,:),POINTER :: R4D => NULL()

      END TYPE VAR

      INTEGER,PARAMETER :: TKR_I0D =1000       !<-- These are simply codes that designate
      INTEGER,PARAMETER :: TKR_I1D =1001       !    whether an internal state variable
      INTEGER,PARAMETER :: TKR_I2D =1002       !    inside of the VARS array is integer
      INTEGER,PARAMETER :: TKR_R0D =1003       !    scalar, 1-d integer array, etc.
      INTEGER,PARAMETER :: TKR_R1D =1004       !    so the type can be referred to easily.
      INTEGER,PARAMETER :: TKR_R2D =1005       !
      INTEGER,PARAMETER :: TKR_R3D =1006       !
      INTEGER,PARAMETER :: TKR_R4D =1007       !<--

      INTERFACE SET_VAR_PTR
        MODULE PROCEDURE SET_VAR_PTR_I0D
        MODULE PROCEDURE SET_VAR_PTR_I1D
        MODULE PROCEDURE SET_VAR_PTR_I2D
        MODULE PROCEDURE SET_VAR_PTR_R1D
        MODULE PROCEDURE SET_VAR_PTR_R0D
        MODULE PROCEDURE SET_VAR_PTR_R2D
        MODULE PROCEDURE SET_VAR_PTR_R3D
        MODULE PROCEDURE SET_VAR_PTR_R4D
      END INTERFACE SET_VAR_PTR

      CONTAINS

!#######################################################################

      SUBROUTINE READ_CONFIG(FNAME,MYPE,VARS,NUM_VARS,RC)

!-----------------------------------------------------------------------
!***  Read the text file for the Solver component that specifies
!***  internal state variables for History, Restart, Import, or
!***  eXport.
!-----------------------------------------------------------------------

        IMPLICIT NONE
        INTEGER, INTENT(IN) :: MYPE
        CHARACTER(LEN=*), INTENT(IN) :: FNAME
        TYPE(VAR), DIMENSION(:), INTENT(INOUT) :: VARS
        INTEGER, INTENT(OUT) :: NUM_VARS
        INTEGER, INTENT(OUT) :: RC

        INTEGER :: IERR,N,IOS,NVARS
        CHARACTER(LEN=256) :: STRING
        CHARACTER(LEN=1) :: CH_H,CH_R,CH_O,CH_I,CH_X,CH_T

!-----------------------------------------------------------------------

        RC = 0

        NVARS = SIZE(VARS)                                                 !<-- Max # of variables MAX_VARS set in internal state modules
        N = 0
        OPEN(UNIT=10,FILE=FNAME,STATUS='OLD',ACTION='READ',IOSTAT=IERR)    !<-- Open the text file with user specifications
        IF(IERR/=0)THEN
          WRITE(0,*)' Unable to open file ',TRIM(FNAME)                 &
                   ,' in READ_CONFIG'
          RC = -1
          RETURN
        ENDIF

        read_specs: DO WHILE(.TRUE.)

          READ(UNIT=10,FMT="(A)",IOSTAT=ios) STRING                        !<-- Insert each variable's specification line into STRING
          IF (IOS/=0) EXIT                                                 !<-- We have now read specification lines for all variables
          IF (STRING(1:1)=='#') CYCLE
          IF (TRIM(STRING)=='') CYCLE
          N = N + 1
          IF (N > NVARS) THEN                                              !<-- # of variables exceeds MAX_VARS
            WRITE(0,*)' increase the size of VARS array ', NVARS
            STOP
          END IF
!
!---------------------------------------------------------------------
!***  Read the text line containing the specifications for variable N
!---------------------------------------------------------------------
!
          READ(UNIT=STRING,FMT=*,IOSTAT=ios) VARS(N)%VBL_NAME,CH_H,CH_R,CH_O,CH_I,CH_X,CH_T,VARS(N)%DESCRIPTION

          IF (IOS/=0) THEN
            IF (IOS>0) THEN
              WRITE(0,*)' error while reading ',FNAME,' on line : ',TRIM(STRING),' iostat = ',ios
              STOP
            ELSE
              N = N - 1
              EXIT
            END IF
          END IF
!
!------------------------------------------------
!***  'Turn on' and store the those qualities
!***  that the user specified in the text file
!***  for variable N.
!------------------------------------------------

          IF (CH_H=='H') VARS(N)%HISTORY=.TRUE.
          IF (CH_H=='P') VARS(N)%HISTORY_P=.TRUE.
          IF (CH_R=='R') VARS(N)%RESTART=.TRUE.
          IF (CH_R=='S') VARS(N)%RESTART_P=.TRUE.
          IF (CH_O=='O') VARS(N)%OWNED  =.TRUE.
          IF (CH_I=='I') VARS(N)%IMPORT =.TRUE.
          IF (CH_X=='X') VARS(N)%EXPORT =.TRUE.
          IF (CH_T=='T') VARS(N)%TSERIES=.TRUE.
          IF (CH_T=='t') VARS(N)%TSERIES_P=.TRUE.

        END DO  read_specs

        NUM_VARS = N
        IF (MYPE==0) THEN
          WRITE(0,*)' NUM_VARS in ',TRIM(FNAME),' ',NUM_VARS
        ENDIF

        CLOSE(UNIT=10)

      END SUBROUTINE READ_CONFIG

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      SUBROUTINE SET_VAR_PTR_I0D (VARS,NUM_VARS,VBL_NAME,I0D)

!--------------------------------------------------------------------
!***  Allocate memory for 'Owned' integer scalars in the Solver
!***  internal state if so directed by ALLOC_FLAG and point those
!***  variables into that memory.
!--------------------------------------------------------------------

        IMPLICIT NONE
        TYPE(VAR), DIMENSION(:), INTENT(INOUT) :: VARS
        INTEGER, INTENT(IN)                    :: NUM_VARS
        CHARACTER(LEN=*), INTENT(IN)           :: VBL_NAME
        INTEGER    ,POINTER                    :: I0D

        INTEGER :: INDX

!-----------------------------------------------------------------------

        CALL FIND_VAR_INDX(VBL_NAME,VARS,NUM_VARS,INDX)
        IF (VARS(INDX)%OWNED) THEN
          ALLOCATE(VARS(INDX)%I0D)                                         !<-- Allocate an integer scalar pointer for input I0D
          VARS(INDX)%I0D=I4_IN
        END IF
        I0D => VARS(INDX)%I0D                                              !<-- Internal state variable I0D uses the newly allocated space
        VARS(INDX)%TKR = TKR_I0D

      END SUBROUTINE SET_VAR_PTR_I0D

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      SUBROUTINE SET_VAR_PTR_R0D (VARS,NUM_VARS,VBL_NAME,R0D)

!--------------------------------------------------------------------
!***  Allocate memory for 'Owned' real scalars in the Solver
!***  internal state if so directed by ALLOC_FLAG and point
!***  those variables into that memory.
!--------------------------------------------------------------------

        IMPLICIT NONE
        TYPE(VAR), DIMENSION(:), INTENT(INOUT) :: VARS
        INTEGER, INTENT(IN)                    :: NUM_VARS
        CHARACTER(LEN=*), INTENT(IN)           :: VBL_NAME
        REAL       ,POINTER                    :: R0D

        INTEGER :: INDX

!-----------------------------------------------------------------------

        CALL FIND_VAR_INDX(VBL_NAME,VARS,NUM_VARS,INDX)
        IF (VARS(INDX)%OWNED) THEN
          ALLOCATE(VARS(INDX)%R0D)
          VARS(INDX)%R0D=R4_IN
        END IF
        R0D => VARS(INDX)%R0D
        VARS(INDX)%TKR = TKR_R0D

      END SUBROUTINE SET_VAR_PTR_R0D

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      SUBROUTINE SET_VAR_PTR_I1D (VARS,NUM_VARS,VBL_NAME,I1D,lowbound,upbound)

!-----------------------------------------------------------------------
!***  Allocate memory for 'Owned' integer 1-D arrays in the Solver
!***  internal state if so directed by ALLOC_FLAG and point those
!***  variables into that memory.
!-----------------------------------------------------------------------

        IMPLICIT NONE
        TYPE(VAR), DIMENSION(:), INTENT(INOUT) :: VARS
        INTEGER, INTENT(IN)                    :: NUM_VARS
        CHARACTER(LEN=*), INTENT(IN)           :: VBL_NAME
        INTEGER, DIMENSION(:), POINTER         :: I1D
        INTEGER, INTENT(IN)                    :: lowbound,upbound

        INTEGER :: INDX

!-----------------------------------------------------------------------

        CALL FIND_VAR_INDX(VBL_NAME,VARS,NUM_VARS,INDX)
        IF (VARS(INDX)%OWNED) THEN
          ALLOCATE (VARS(INDX)%I1D(lowbound:upbound))
          VARS(INDX)%I1D=I4_IN
        END IF
        I1D => VARS(INDX)%I1D
        VARS(INDX)%TKR = TKR_I1D

      END SUBROUTINE SET_VAR_PTR_I1D

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      SUBROUTINE SET_VAR_PTR_I2D (VARS,NUM_VARS,VBL_NAME,I2D,lowbound,upbound)

!-----------------------------------------------------------------------
!***  Allocate memory for 'Owned' integer 2-D arrays in the Solver
!***  internal state if so directed by ALLOC_FLAG and point those
!***  variables into that memory.
!-----------------------------------------------------------------------

        IMPLICIT NONE
        TYPE(VAR), DIMENSION(:), INTENT(INOUT) :: VARS
        INTEGER, INTENT(IN)                    :: NUM_VARS
        CHARACTER(LEN=*), INTENT(IN)           :: VBL_NAME
        INTEGER,DIMENSION(:,:)  ,POINTER       :: I2D
        INTEGER, DIMENSION(2), INTENT(IN)      :: lowbound,upbound

        INTEGER :: INDX

!-----------------------------------------------------------------------

        CALL FIND_VAR_INDX(VBL_NAME,VARS,NUM_VARS,INDX)
        IF (VARS(INDX)%OWNED) THEN
          ALLOCATE (VARS(INDX)%I2D(lowbound(1):upbound(1), &
                                   lowbound(2):upbound(2) ))
          VARS(INDX)%I2D=I4_IN
        END IF
        I2D => VARS(INDX)%I2D
        VARS(INDX)%TKR = TKR_I2D

      END SUBROUTINE SET_VAR_PTR_I2D

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      SUBROUTINE SET_VAR_PTR_R1D (VARS,NUM_VARS,VBL_NAME,R1D,lowbound,upbound)

!-----------------------------------------------------------------------
!***  Allocate memory for 'Owned' real 1-D arrays in the Solver
!***  internal state if so directed by ALLOC_FLAG and point those
!***  variables into that memory.
!-----------------------------------------------------------------------

        IMPLICIT NONE
        TYPE(VAR), DIMENSION(:), INTENT(INOUT) :: VARS
        INTEGER, INTENT(IN)                    :: NUM_VARS
        CHARACTER(LEN=*), INTENT(IN)           :: VBL_NAME
        REAL, DIMENSION(:), POINTER            :: R1D
        INTEGER, INTENT(IN)                    :: lowbound,upbound

        INTEGER :: INDX

!-----------------------------------------------------------------------

        CALL FIND_VAR_INDX(VBL_NAME,VARS,NUM_VARS,INDX)
        IF (VARS(INDX)%OWNED) THEN
          ALLOCATE (VARS(INDX)%R1D(lowbound:upbound))
          VARS(INDX)%R1D=R4_IN
        END IF
        R1D => VARS(INDX)%R1D
        VARS(INDX)%TKR = TKR_R1D

      END SUBROUTINE SET_VAR_PTR_R1D

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      SUBROUTINE SET_VAR_PTR_R2D (VARS,NUM_VARS,VBL_NAME,R2D,lowbound,upbound)

!-----------------------------------------------------------------------
!***  Allocate memory for 'Owned' real 2-D arrays in the Solver
!***  internal state if so directed by ALLOC_FLAG and point those
!***  variables into that memory.
!-----------------------------------------------------------------------

        IMPLICIT NONE
        TYPE(VAR), DIMENSION(:), INTENT(INOUT) :: VARS
        INTEGER, INTENT(IN)                    :: NUM_VARS
        CHARACTER(LEN=*), INTENT(IN)           :: VBL_NAME
        REAL   ,DIMENSION(:,:)  ,POINTER       :: R2D
        INTEGER, DIMENSION(2), INTENT(IN)      :: lowbound,upbound

        INTEGER :: INDX

!-----------------------------------------------------------------------

        CALL FIND_VAR_INDX(VBL_NAME,VARS,NUM_VARS,INDX)
        IF (VARS(INDX)%OWNED) THEN
          ALLOCATE (VARS(INDX)%R2D(lowbound(1):upbound(1), &
                                   lowbound(2):upbound(2) ))
          VARS(INDX)%R2D=R4_IN
        END IF
        R2D => VARS(INDX)%R2D
        VARS(INDX)%TKR = TKR_R2D

      END SUBROUTINE SET_VAR_PTR_R2D

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      SUBROUTINE SET_VAR_PTR_R3D (VARS,NUM_VARS,VBL_NAME,R3D,lowbound,upbound)

!-----------------------------------------------------------------------
!***  Allocate memory for 'Owned' real 3-D arrays in the Solver
!***  internal state if so directed by ALLOC_FLAG and point those
!***  variables into that memory.
!-----------------------------------------------------------------------

        IMPLICIT NONE
        TYPE(VAR), DIMENSION(:), INTENT(INOUT) :: VARS
        INTEGER, INTENT(IN)                    :: NUM_VARS
        CHARACTER(LEN=*), INTENT(IN)           :: VBL_NAME
        REAL   ,DIMENSION(:,:,:),POINTER       :: R3D
        INTEGER, DIMENSION(3), INTENT(IN)      :: lowbound,upbound

        INTEGER :: INDX

!-----------------------------------------------------------------------

        CALL FIND_VAR_INDX(VBL_NAME,VARS,NUM_VARS,INDX)
        IF (VARS(INDX)%OWNED) THEN
          ALLOCATE (VARS(INDX)%R3D(lowbound(1):upbound(1), &
                                   lowbound(2):upbound(2), &
                                   lowbound(3):upbound(3) ))
          VARS(INDX)%R3D=R4_IN
        END IF
        R3D => VARS(INDX)%R3D
        VARS(INDX)%TKR = TKR_R3D

      END SUBROUTINE SET_VAR_PTR_R3D

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      SUBROUTINE SET_VAR_PTR_R4D (VARS,NUM_VARS,VBL_NAME,R4D,lowbound,upbound)

!-----------------------------------------------------------------------
!***  Allocate memory for 'Owned' real 4-D arrays in the Solver
!***  internal state if so directed by ALLOC_FLAG and point those
!***  variables into that memory.
!-----------------------------------------------------------------------

        IMPLICIT NONE
        TYPE(VAR), DIMENSION(:), INTENT(INOUT) :: VARS
        INTEGER, INTENT(IN)                    :: NUM_VARS
        CHARACTER(LEN=*), INTENT(IN)           :: VBL_NAME
        REAL   ,DIMENSION(:,:,:,:),POINTER     :: R4D
        INTEGER, DIMENSION(4), INTENT(IN)      :: lowbound,upbound

        INTEGER :: INDX

!-----------------------------------------------------------------------

        CALL FIND_VAR_INDX(VBL_NAME,VARS,NUM_VARS,INDX)
        IF (VARS(INDX)%OWNED) THEN
          ALLOCATE (VARS(INDX)%R4D(lowbound(1):upbound(1), &
                                   lowbound(2):upbound(2), &
                                   lowbound(3):upbound(3), &
                                   lowbound(4):upbound(4) ))
          VARS(INDX)%R4D=R4_IN
        END IF
        R4D => VARS(INDX)%R4D
        VARS(INDX)%TKR = TKR_R4D

      END SUBROUTINE SET_VAR_PTR_R4D

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      SUBROUTINE DEALLOC_VARS(VARS,NUM_VARS)

!----------------------------------------------------------------------
!***  Deallocate the memory that had been allocated within the VARS
!***  composite array into which Solver internal state variables are
!***  pointing.
!----------------------------------------------------------------------

      IMPLICIT NONE

      TYPE(VAR), DIMENSION(:), INTENT(INOUT) :: VARS
      INTEGER, INTENT(IN) :: NUM_VARS

      INTEGER :: N
      INTEGER :: ISTAT

!-----------------------------------------------------------------------

      DO N=1,NUM_VARS
        IF (VARS(N)%OWNED) THEN
          SELECT CASE(VARS(N)%TKR)
            CASE(TKR_I0D)
              DEALLOCATE(VARS(N)%I0D,STAT=ISTAT)
            CASE(TKR_I1D)
              DEALLOCATE(VARS(N)%I1D,STAT=ISTAT)
            CASE(TKR_I2D)
              DEALLOCATE(VARS(N)%I2D,STAT=ISTAT)
            CASE(TKR_R0D)
              DEALLOCATE(VARS(N)%R0D,STAT=ISTAT)
            CASE(TKR_R1D)
              DEALLOCATE(VARS(N)%R1D,STAT=ISTAT)
            CASE(TKR_R2D)
              DEALLOCATE(VARS(N)%R2D,STAT=ISTAT)
            CASE(TKR_R3D)
              DEALLOCATE(VARS(N)%R3D,STAT=ISTAT)
            CASE(TKR_R4D)
              DEALLOCATE(VARS(N)%R4D,STAT=ISTAT)
            CASE DEFAULT
              write(0,*)' Unknown TKR in DEALLOC_VARS ', VARS(N)%VBL_NAME, VARS(N)%TKR
              stop 9
          END SELECT
        END IF
      END DO

      END SUBROUTINE DEALLOC_VARS

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      SUBROUTINE FIND_VAR_INDX(VBL_NAME,VARS,NUM_VARS,INDX)

!----------------------------------------------------------------
!***  Find the location (index) within the VARS array with which
!***  the internal state variable called VBL_NAME is associated.
!----------------------------------------------------------------
      USE MPI
      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN)        :: VBL_NAME
      TYPE(VAR), DIMENSION(:), INTENT(IN) :: VARS
      INTEGER, INTENT(IN)                 :: NUM_VARS
      INTEGER, INTENT(OUT)                :: INDX

      INTEGER :: I,IERR

      INDX = 0
      DO I=1,NUM_VARS
        IF (TRIM(VBL_NAME) == TRIM(VARS(I)%VBL_NAME) ) THEN
          INDX = I
          EXIT
        END IF
      END DO
      IF (INDX == 0) THEN
138    format(' can not find |',A,'| in solver state text file.')
       write(0,138) trim(VBL_NAME)
       call MPI_Abort(MPI_COMM_WORLD,2,ierr)
       stop 2
      END IF

      END SUBROUTINE FIND_VAR_INDX

!#######################################################################

      END MODULE MODULE_VARS
