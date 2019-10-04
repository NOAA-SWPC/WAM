!-------------------------------------------------------------------------------

      MODULE MODULE_VARS_STATE

!-------------------------------------------------------------------------------

      USE ESMF
      USE MODULE_VARS
      USE module_CONTROL,ONLY: TIMEF

      IMPLICIT NONE

      PRIVATE
      PUBLIC :: PUT_VARS_IN_STATE
      PUBLIC :: GET_VARS_FROM_STATE
      PUBLIC :: DELETE_FIELDS_FROM_STATE
      PUBLIC :: PUT_VARS_IN_BUNDLES

      integer,parameter :: double=selected_real_kind(p=13,r=200)
      integer,parameter:: kdbl=double
      real(kind=kdbl) :: btim,btim0

      CONTAINS

!###############################################################################

      SUBROUTINE PUT_VARS_IN_STATE(VARS, NUM_VARS, STATE_TYPE, GRID, STATE)

!-------------------------------------------------------------------------------
!***  Put pointers of variables in the composite VARS array into a given
!***  ESMF import/export state.
!-------------------------------------------------------------------------------

      IMPLICIT NONE

      TYPE(VAR), DIMENSION(:), INTENT(IN)    :: VARS
      INTEGER,                 INTENT(IN)    :: NUM_VARS
      CHARACTER(LEN=1),        INTENT(IN)    :: STATE_TYPE
      TYPE(ESMF_Grid),         INTENT(IN)    :: GRID
      TYPE(ESMF_State ),       INTENT(INOUT) :: STATE

      TYPE(ESMF_StateItem_Flag):: stateItemType
      TYPE(ESMF_Field)         :: FIELD
      INTEGER                  :: KOUNT,N, RC
      INTEGER                  :: IHALO,JHALO

!-------------------------------------------------------------------------------

!!!!!!!!! FIX this later
!!!!!!!!! FIX this later
                IHALO = 3
                JHALO = 3
!!!!!!!!! FIX this later
!!!!!!!!! FIX this later

      DO N=1,NUM_VARS
        IF (( (STATE_TYPE == 'I') .AND. VARS(N)%IMPORT ) .OR. &
            ( (STATE_TYPE == 'X') .AND. VARS(N)%EXPORT ) ) THEN

          SELECT CASE(VARS(N)%TKR)
            CASE(TKR_I0D)
              CALL ESMF_AttributeSet(state=STATE ,name=VARS(N)%VBL_NAME, value=VARS(N)%I0D, rc=RC)
            CASE(TKR_I1D)
            CASE(TKR_I2D)
            CASE(TKR_R0D)
              CALL ESMF_AttributeSet(state=STATE ,name=VARS(N)%VBL_NAME, value=VARS(N)%R0D, rc=RC)
            CASE(TKR_R1D)
              KOUNT=SIZE(VARS(N)%R1D)
              CALL ESMF_AttributeSet(state=STATE ,name=VARS(N)%VBL_NAME, itemCount=KOUNT, valueList=VARS(N)%R1D, rc=RC)
            CASE(TKR_R2D)
              CALL ESMF_StateGet(STATE ,VARS(N)%VBL_NAME , stateItemType, rc=RC)
              IF (stateItemType==ESMF_STATEITEM_NOTFOUND) THEN
                FIELD = ESMF_FieldCreate(grid       =GRID                       &
                                        ,farray     =VARS(N)%R2D                &
                                        ,totalUWidth=(/IHALO,JHALO/)            &  !<-- Upper bound of halo region
                                        ,totalLWidth=(/IHALO,JHALO/)            &  !<-- Lower bound of halo region
                                        ,name       =VARS(N)%VBL_NAME           &
                                        ,indexFlag  =ESMF_INDEX_GLOBAL          &
                                        ,rc         =RC)
                CALL ESMF_StateAddReplace(STATE ,(/FIELD/) ,rc=RC)
              ENDIF
            CASE(TKR_R3D)
              CALL ESMF_StateGet(STATE ,VARS(N)%VBL_NAME ,stateItemType, rc=RC)
              IF (stateItemType==ESMF_STATEITEM_NOTFOUND .and. ASSOCIATED(VARS(N)%R3D)) THEN
                FIELD = ESMF_FieldCreate(grid            =GRID                              &
                                        ,farray          =VARS(N)%R3D                       &
                                        ,totalUWidth     =(/IHALO,JHALO/)                   &  !<-- Upper bound of halo region
                                        ,totalLWidth     =(/IHALO,JHALO/)                   &  !<-- Lower bound of halo region
                                        ,ungriddedLBound =(/ lbound(VARS(N)%R3D,dim=3) /)   &
                                        ,ungriddedUBound =(/ ubound(VARS(N)%R3D,dim=3) /)   &
                                        ,name            =VARS(N)%VBL_NAME                  &
                                        ,indexFlag       =ESMF_INDEX_GLOBAL                 &
                                        ,rc              =RC)
                CALL ESMF_StateAddReplace(STATE ,(/FIELD/) ,rc=RC)
              ENDIF
            CASE(TKR_R4D)
              CALL ESMF_StateGet(STATE ,VARS(N)%VBL_NAME ,stateItemType, rc=RC)
              IF (stateItemType==ESMF_STATEITEM_NOTFOUND) THEN
                FIELD = ESMF_FieldCreate(grid            =GRID                              &
                                        ,farray          =VARS(N)%R4D                       &
                                        ,totalUWidth     =(/IHALO,JHALO/)                   &  !<-- Upper bound of halo region
                                        ,totalLWidth     =(/IHALO,JHALO/)                   &  !<-- Lower bound of halo region
                                        ,ungriddedLBound =(/ lbound(VARS(N)%R4D,dim=3),lbound(VARS(N)%R4D,dim=4) /)   &
                                        ,ungriddedUBound =(/ ubound(VARS(N)%R4D,dim=3),ubound(VARS(N)%R4D,dim=4) /)   &
                                        ,name            =VARS(N)%VBL_NAME                  &
                                        ,indexFlag       =ESMF_INDEX_GLOBAL                 &
                                        ,rc              =RC)
                CALL ESMF_StateAddReplace(STATE ,(/FIELD/),rc=RC)
              ENDIF
            CASE DEFAULT
              write(0,*)' TKR = ', VARS(N)%TKR, TRIM(VARS(N)%VBL_NAME)
              write(0,*)' This TKR Case not available in PUT_VARS_IN_STATE'
              write(0,*)' ABORTING!'
              CALL ESMF_FINALIZE(endflag=ESMF_END_ABORT                     &
                                ,rc             =RC)
          END SELECT

        ENDIF
      END DO

      END SUBROUTINE PUT_VARS_IN_STATE

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      SUBROUTINE GET_VARS_FROM_STATE(VARS, NUM_VARS, STATE)

!-------------------------------------------------------------------------------
!***  Take allocated pointers from a given ESMF state and point VARS 
!***  locations at them if the VARS variable is unowned/unallocated
!***  or move the pointer data into the VARS location if owned by 
!***  multiple components.
!-------------------------------------------------------------------------------

      IMPLICIT NONE

      TYPE(VAR), DIMENSION(:), INTENT(INOUT) :: VARS
      INTEGER,                 INTENT(IN)    :: NUM_VARS

      TYPE(ESMF_State ),       INTENT(INOUT) :: STATE

      TYPE(ESMF_Field) :: FIELD
      INTEGER :: KOUNT,N, RC
      INTEGER                            :: HOLD_I0D
      INTEGER,DIMENSION(:)      ,POINTER :: HOLD_I1D
      INTEGER,DIMENSION(:,:)    ,POINTER :: HOLD_I2D
      REAL                               :: HOLD_R0D
      REAL   ,DIMENSION(:)      ,POINTER :: HOLD_R1D
      REAL   ,DIMENSION(:,:)    ,POINTER :: HOLD_R2D
      REAL   ,DIMENSION(:,:,:)  ,POINTER :: HOLD_R3D
      REAL   ,DIMENSION(:,:,:,:),POINTER :: HOLD_R4D
      CHARACTER(LEN=32) :: VBL_NAME

!-------------------------------------------------------------------------------

      DO N=1,NUM_VARS
        IF ( VARS(N)%IMPORT ) THEN
          VBL_NAME = TRIM(VARS(N)%VBL_NAME)
          SELECT CASE(VARS(N)%TKR)
            CASE(TKR_I0D)
              CALL ESMF_AttributeGet(state=STATE ,name=VARS(N)%VBL_NAME, value=VARS(N)%I0D, rc=RC)
              if (rc/=ESMF_SUCCESS) then
                write(0,*)' Unable to get VBL_NAME for CASE(TKR_I0D) in GET_VARS_FROM_STATE'
                write(0,*)' ABORTING!'
                CALL ESMF_FINALIZE(endflag=ESMF_END_ABORT                   &
                                  ,rc             =RC)
              end if
            CASE(TKR_I1D)
              write(0,*)' not implemented TKR_I1D in GET_VARS_FROM_STATE '
              write(0,*)' ABORTING!'
              CALL ESMF_FINALIZE(endflag=ESMF_END_ABORT                     &
                                ,rc             =RC)
            CASE(TKR_I2D)
              CALL ESMF_StateGet(STATE ,VBL_NAME ,FIELD ,rc=RC)
              if (rc/=ESMF_SUCCESS) then
                write(0,*)' Unable to get VBL_NAME for CASE(TKR_I2D) in GET_VARS_FROM_STATE'
                write(0,*)' ABORTING!'
                CALL ESMF_FINALIZE(endflag=ESMF_END_ABORT                   &
                                  ,rc             =RC)
              end if

              CALL ESMF_FieldGet(field=FIELD ,localDe=0 ,farrayPtr=HOLD_I2D ,rc=RC)

              if (rc/=ESMF_SUCCESS) then
                write(0,*)' Unable to get 2D integer array from Field in GET_VARS_FROM_STATE'
                write(0,*)' ABORTING!'
                CALL ESMF_FINALIZE(endflag=ESMF_END_ABORT                   &
                                  ,rc             =RC)
              end if
              IF (VARS(N)%OWNED ) THEN
                if (size(VARS(N)%I2D) /= size(HOLD_I2D) ) then
                  write(0,*)' size(VARS(N)%I2D) /= size(HOLD_I2D) in GET_VARS_FROM_STATE'
                  write(0,*)' ABORTING!'
                  CALL ESMF_FINALIZE(endflag=ESMF_END_ABORT                 &
                                    ,rc             =RC)
                end if
                VARS(N)%I2D =  HOLD_I2D        !<-- Transfer data since multiply owned
              ELSE
                VARS(N)%I2D => HOLD_I2D        !<-- Point the appropriate unallocated VARS location at allocated pointer
              END IF
            CASE(TKR_R0D)
              CALL ESMF_AttributeGet(state=STATE ,name=VARS(N)%VBL_NAME, value=VARS(N)%R0D, rc=RC)
              if (rc/=ESMF_SUCCESS) then
                write(0,*)' Unable to get VBL_NAME for CASE(TKR_R0D) in GET_VARS_FROM_STATE'
                write(0,*)' ABORTING!'
                CALL ESMF_FINALIZE(endflag=ESMF_END_ABORT                   &
                                  ,rc             =RC)
              end if
            CASE(TKR_R1D)
              KOUNT=SIZE(VARS(N)%R1D)
              CALL ESMF_AttributeGet(state=STATE ,name=VARS(N)%VBL_NAME, itemCount=KOUNT,  valueList=VARS(N)%R1D, rc=RC)
!!!           write(0,*)' not implemented TKR_R1D in UPDATE_VARS '
!!!           stop
            CASE(TKR_R2D)
              CALL ESMF_StateGet(state=STATE ,itemName=VBL_NAME ,field=FIELD ,rc=RC)
              if (rc/=ESMF_SUCCESS) then
                write(0,*)' Unable to get ',trim(VBL_NAME),' for CASE(TKR_R2D) in GET_VARS_FROM_STATE'
                write(0,*)' ABORTING!'
                CALL ESMF_FINALIZE(endflag=ESMF_END_ABORT                   &
                                  ,rc             =RC)
              end if

              CALL ESMF_FieldGet(field=FIELD ,localDe=0 ,farrayPtr=HOLD_R2D ,rc=RC)

              if (rc/=ESMF_SUCCESS) then
                write(0,*)' Unable to get 2D real array from Field in GET_VARS_FROM_STATE'
                write(0,*)' ABORTING!'
                CALL ESMF_FINALIZE(endflag=ESMF_END_ABORT                   &
                                  ,rc             =RC)
              end if
              IF (VARS(N)%OWNED ) THEN
                if (size(VARS(N)%R2D) /= size(HOLD_R2D) ) then
                  write(0,*)' size(VARS(N)%R2D) /= size(HOLD_R2D) in GET_VARS_FROM_STATE'
                  write(0,*)' ABORTING!'
                  CALL ESMF_FINALIZE(endflag=ESMF_END_ABORT                 &
                                    ,rc             =RC)
                end if
                VARS(N)%R2D =  HOLD_R2D         !<-- Transfer data since multiply owned
              ELSE
                VARS(N)%R2D => HOLD_R2D         !<-- Point the appropriate unallocated VARS location at allocated pointer
              END IF
            CASE(TKR_R3D)
              CALL ESMF_StateGet(STATE ,VBL_NAME ,FIELD ,rc=RC)
              if (rc/=ESMF_SUCCESS) then
                write(0,*)' Unable to get VBL_NAME for CASE(TKR_R3D) in GET_VARS_FROM_STATE'
                write(0,*)' ABORTING!'
                CALL ESMF_FINALIZE(endflag=ESMF_END_ABORT                   &
                                  ,rc             =RC)
              end if

              CALL ESMF_FieldGet(field=FIELD ,localDe=0 ,farrayPtr=HOLD_R3D ,rc=RC)

              if (rc/=ESMF_SUCCESS) then
                write(0,*)' Unable to get 3D real array from Field in GET_VARS_FROM_STATE'
                write(0,*)' ABORTING!'
                CALL ESMF_FINALIZE(endflag=ESMF_END_ABORT                   &
                                  ,rc             =RC)
              end if
              IF (VARS(N)%OWNED ) THEN
                if (size(VARS(N)%R3D) /= size(HOLD_R3D) ) then
                  write(0,*)TRIM(VARS(N)%VBL_NAME), size(VARS(N)%R3D), size(HOLD_R3D)
                  write(0,*)TRIM(VARS(N)%VBL_NAME), 'lbound ',lbound(VARS(N)%R3D), lbound(HOLD_R3D)
                  write(0,*)TRIM(VARS(N)%VBL_NAME), 'ubound ',ubound(VARS(N)%R3D), ubound(HOLD_R3D)
                  write(0,*)' VARS(N)%R3D) /= size(HOLD_R3D in GET_VARS_FROM_STATE'
                  write(0,*)' ABORTING!'
                  CALL ESMF_FINALIZE(endflag=ESMF_END_ABORT                 &
                                    ,rc             =RC)
                end if
                VARS(N)%R3D =  HOLD_R3D         !<-- Transfer data since multiply owned
              ELSE
                VARS(N)%R3D => HOLD_R3D         !<-- Point the appropriate unallocated VARS location at allocated pointer
              END IF
            CASE(TKR_R4D)
              CALL ESMF_StateGet(STATE ,VBL_NAME ,FIELD ,rc=RC)
              if (rc/=ESMF_SUCCESS) then
                write(0,*)' Unable to get VBL_NAME for CASE(TKR_R4D) in GET_VARS_FROM_STATE'
                write(0,*)' ABORTING!'
                CALL ESMF_FINALIZE(endflag=ESMF_END_ABORT                   &
                                  ,rc             =RC)
              end if

              CALL ESMF_FieldGet(field=FIELD ,localDe=0 ,farrayPtr=HOLD_R4D ,rc=RC)

              if (rc/=ESMF_SUCCESS) then
                write(0,*)' Unable to get 4D real array from Field in GET_VARS_FROM_STATE'
                write(0,*)' ABORTING!'
                CALL ESMF_FINALIZE(endflag=ESMF_END_ABORT                   &
                                  ,rc             =RC)
              end if
              IF (VARS(N)%OWNED ) THEN
                if (size(VARS(N)%R4D) /= size(HOLD_R4D) ) then
                  write(0,*)TRIM(VARS(N)%VBL_NAME), size(VARS(N)%R4D), size(HOLD_R4D)
                  write(0,*)TRIM(VARS(N)%VBL_NAME), 'lbound ',lbound(VARS(N)%R4D), lbound(HOLD_R4D)
                  write(0,*)TRIM(VARS(N)%VBL_NAME), 'ubound ',ubound(VARS(N)%R4D), ubound(HOLD_R4D)
                  write(0,*)' size(VARS(N)%R4D) /= size(HOLD_R4D) in GET_VARS_FROM_STATE'
                  write(0,*)' ABORTING!'
                  CALL ESMF_FINALIZE(endflag=ESMF_END_ABORT                 &
                                    ,rc             =RC)
                end if
                VARS(N)%R4D =  HOLD_R4D         !<-- Transfer data since multiply owned
              ELSE
                VARS(N)%R4D => HOLD_R4D         !<-- Point the appropriate unallocated VARS location at allocated pointer
              END IF
            CASE DEFAULT
              write(0,*)' TKR = ', VARS(N)%TKR, TRIM(VARS(N)%VBL_NAME)
              write(0,*)' This TKR Case is not available in GET_VARS_FROM_STATE'
              write(0,*)' ABORTING!'
              CALL ESMF_FINALIZE(endflag=ESMF_END_ABORT                     &
                                ,rc             =RC)
          END SELECT

        END IF
      END DO

      END SUBROUTINE GET_VARS_FROM_STATE

!-------------------------------------------------------------------------------
!###############################################################################
!-------------------------------------------------------------------------------

      SUBROUTINE DELETE_FIELDS_FROM_STATE(STATE)

      IMPLICIT NONE

      TYPE(ESMF_State ), INTENT(INOUT) :: STATE

      INTEGER                                :: RC
      INTEGER                                :: i, itemcount
      TYPE(ESMF_Field)                       :: FIELD

      CHARACTER(LEN=10), DIMENSION(:), ALLOCATABLE :: itemNameList

      CALL ESMF_StateGet(STATE,itemcount=itemcount, rc=RC)
      ALLOCATE(itemNameList(itemcount))
      CALL ESMF_StateGet(STATE,itemcount=itemcount,itemNameList=itemNameList, rc=RC)
      DO i=1,itemcount
        CALL ESMF_StateGet(STATE ,itemNameList(i) ,FIELD ,rc=RC)
        CALL ESMF_FieldDestroy(field=FIELD, rc=RC)
      END DO
      DEALLOCATE(itemNameList)

      END SUBROUTINE DELETE_FIELDS_FROM_STATE

!-------------------------------------------------------------------------------
!###############################################################################
!-------------------------------------------------------------------------------

      SUBROUTINE PUT_VARS_IN_BUNDLES(VARS                                       &
                                    ,NUM_VARS                                   &
                                    ,GRID                                       &
                                    ,HISTORY_BUNDLE                             &
                                    ,RESTART_BUNDLE)

      USE MODULE_ERROR_MSG,ONLY: ERR_MSG,MESSAGE_CHECK

      IMPLICIT NONE

!------------------------
!***  Argument Variables
!------------------------

      TYPE(VAR), DIMENSION(:), INTENT(IN)    :: VARS
      INTEGER,                 INTENT(IN)    :: NUM_VARS
      TYPE(ESMF_Grid),         INTENT(IN)    :: GRID
      TYPE(ESMF_FieldBundle),  INTENT(INOUT) :: HISTORY_BUNDLE
      TYPE(ESMF_FieldBundle),  INTENT(INOUT) :: RESTART_BUNDLE

!---------------------
!***  Local variables
!---------------------

      INTEGER                      :: K,LENGTH                          &
                                     ,N,M,NDIM3,NFIND                   &
                                     ,RC,RC_OUT

      INTEGER                  :: IHALO,JHALO

      INTEGER :: LDIM1,LDIM2,LDIM3,LDIM4,UDIM1,UDIM2,UDIM3,UDIM4

      CHARACTER(3)           :: MODEL_LEVEL
      CHARACTER(3)           :: TRACERS_KIND
      CHARACTER(6)           :: FMT3='(I3.3)'
      CHARACTER(6)           :: FMT2='(I2.2)'
      CHARACTER(ESMF_MAXSTR) :: VBL_NAME,VBL_NAME_X

      TYPE(ESMF_Field)       :: FIELD

      TYPE(ESMF_DataCopy_Flag)    :: COPYFLAG=ESMF_DATACOPY_REFERENCE
!     TYPE(ESMF_DataCopy_Flag)    :: COPYFLAG=ESMF_DATA_COPY

      INTEGER :: INDX_Q2  !! FIXME

!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Begin with the integer scalars.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Insert Integer Scalars into History Bundle"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      DO N=1,NUM_VARS
        IF (VARS(N)%TKR == TKR_I0D) THEN

          IF (VARS(N)%HISTORY) THEN                                        !<-- Take integer scalar data specified for history output
            CALL ESMF_AttributeSet(FIELDBUNDLE=HISTORY_BUNDLE           &  !<-- The Write component output history Bundle
                                  ,name       =VARS(N)%VBL_NAME         &  !<-- Name of the integer scalar
                                  ,value      =VARS(N)%I0D              &  !<-- The scalar being inserted into the import state
                                  ,rc         =RC)
          END IF
          IF (VARS(N)%RESTART) THEN                                        !<-- Take integer scalar data specified for restart output
            CALL ESMF_AttributeSet(FIELDBUNDLE=RESTART_BUNDLE           &  !<-- The Write component output restart Bundle
                                  ,name       =VARS(N)%VBL_NAME         &  !<-- Name of the integer scalar
                                  ,value      =VARS(N)%I0D              &  !<-- The scalar being inserted into the import state
                                  ,rc         =RC)
          END IF

        END IF
      END DO
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  The real scalars.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Insert Real Scalars into History Bundle"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      DO N=1,NUM_VARS
        IF (VARS(N)%TKR == TKR_R0D) THEN

          IF (VARS(N)%HISTORY) THEN                                        !<-- Take real scalar data specified for history output
            CALL ESMF_AttributeSet(FIELDBUNDLE=HISTORY_BUNDLE           &  !<-- The Write component output history Bundle
                                  ,name       =VARS(N)%VBL_NAME         &  !<-- Name of the real scalar
                                  ,value      =VARS(N)%R0D              &  !<-- The scalar being inserted into the history Bundle
                                  ,rc         =RC)
          END IF
          IF (VARS(N)%RESTART) THEN                                        !<-- Take real scalar data specified for restart output
            CALL ESMF_AttributeSet(FIELDBUNDLE=RESTART_BUNDLE           &  !<-- The Write component output restart Bundle
                                  ,name       =VARS(N)%VBL_NAME         &  !<-- Name of the real scalar
                                  ,value      =VARS(N)%R0D              &  !<-- The scalar being inserted into the restart Bundle
                                  ,rc         =RC)
          END IF

        END IF
      END DO
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  The 1-D integer arrays.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Insert 1-D Integer Arrays into History Bundle"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      DO N=1,NUM_VARS
        IF (VARS(N)%TKR == TKR_I1D) THEN
          LENGTH=SIZE(VARS(N)%I1D)

          IF (VARS(N)%HISTORY) THEN                                        !<-- Take 1D integer array data specified for history output
            CALL ESMF_AttributeSet(FIELDBUNDLE   =HISTORY_BUNDLE        &  !<-- The Write component output history Bundle
                                  ,name          =VARS(N)%VBL_NAME      &  !<-- Name of the integer array
                                  ,itemCount     =LENGTH                &  !<-- # of elements in this attribute
                                  ,valueList     =VARS(N)%I1D           &  !<-- The 1D integer being inserted into the history Bundle
                                  ,rc            =RC)
          END IF
          IF (VARS(N)%RESTART) THEN                                        !<-- Take 1D integer array data specified for restart output
            CALL ESMF_AttributeSet(FIELDBUNDLE   =RESTART_BUNDLE        &  !<-- The Write component output restart Bundle
                                  ,name          =VARS(N)%VBL_NAME      &  !<-- Name of the integer array
                                  ,itemCount     =LENGTH                &  !<-- # of elements in this attribute
                                  ,valueList     =VARS(N)%I1D           &  !<-- The 1D integer being inserted into the restart Bundle
                                  ,rc            =RC)
          END IF

        END IF
      END DO
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  The 1-D real arrays.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Insert 1-D Real Arrays into History Bundle"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      DO N=1,NUM_VARS
        IF (VARS(N)%TKR == TKR_R1D) THEN
          LENGTH=SIZE(VARS(N)%R1D)

          IF (VARS(N)%HISTORY) THEN                                        !<-- Take 1D real array data specified for history output
            CALL ESMF_AttributeSet(FIELDBUNDLE   =HISTORY_BUNDLE        &  !<-- The Write component output history Bundle
                                  ,name          =VARS(N)%VBL_NAME      &  !<-- Name of the real array
                                  ,itemCount     =LENGTH                &  !<-- # of elements in this attribute
                                  ,valueList     =VARS(N)%R1D           &  !<-- The 1D real being inserted into the history Bundle
                                  ,rc            =RC)
          END IF
          IF (VARS(N)%RESTART) THEN                                        !<-- Take 1D real array data specified for restart output
            CALL ESMF_AttributeSet(FIELDBUNDLE   =RESTART_BUNDLE        &  !<-- The Write component output restart Bundle
                                  ,name          =VARS(N)%VBL_NAME      &  !<-- Name of the real array
                                  ,itemCount     =LENGTH                &  !<-- # of elements in this attribute
                                  ,valueList     =VARS(N)%R1D           &  !<-- The 1D real being inserted into the restart Bundle
                                  ,rc            =RC)
          END IF

        END IF
      END DO
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  The 2-D integer arrays.
!-----------------------------------------------------------------------
!
      IHALO=3
      JHALO=3
!
      DO N=1,NUM_VARS
        IF (VARS(N)%TKR == TKR_I2D) THEN
          IF (VARS(N)%HISTORY) THEN                                        !<-- Take 2D integer array data specified for history output
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Insert 2-D Integer Data into Field"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          FIELD=ESMF_FieldCreate(grid         =GRID                     &  !<-- The ESMF grid
                                ,farray       =VARS(N)%I2D              &  !<-- The 2D integer array being inserted into history Bundle
                                ,datacopyflag =COPYFLAG                 &
                                ,totalUWidth  =(/IHALO,JHALO/)          &
                                ,totalLWidth  =(/IHALO,JHALO/)          &
                                ,name         =VARS(N)%VBL_NAME         &  !<-- Name of the 2D real array
                                ,indexFlag    =ESMF_INDEX_GLOBAL        &
                                ,rc           =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Insert Solver 2-D Integer Field into History Bundles"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_FieldBundleAdd(            HISTORY_BUNDLE           &  !<-- The Write component output history Bundle
                                  ,            (/FIELD/)       &  !<-- ESMF Field holding the 2D real array
                                  ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          END IF
          IF (VARS(N)%RESTART) THEN                                        !<-- Take 2D integer array data specified for restart output
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Insert 2-D Integer Data into Field"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          FIELD=ESMF_FieldCreate(grid         =GRID                     &  !<-- The ESMF grid
                                ,farray       =VARS(N)%I2D              &  !<-- The 2D integer array being inserted into restart Bundle
                                ,datacopyflag =COPYFLAG                 &
                                ,totalUWidth  =(/IHALO,JHALO/)          &
                                ,totalLWidth  =(/IHALO,JHALO/)          &
                                ,name         =VARS(N)%VBL_NAME         &  !<-- Name of the 2D real array
                                ,indexFlag    =ESMF_INDEX_GLOBAL        &
                                ,rc           =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Insert 2-D Integer Field into Restart Bundles"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_FieldBundleAdd(            RESTART_BUNDLE           &  !<-- The Write component output restart Bundle
                                  ,            (/FIELD/)       &  !<-- ESMF Field holding the 2D real array
                                  ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          END IF
        END IF
      END DO

!
!-----------------------------------------------------------------------
!***  The 2-D real arrays.
!-----------------------------------------------------------------------
!
      DO N=1,NUM_VARS
        IF (VARS(N)%TKR == TKR_R2D) THEN
          IF (VARS(N)%HISTORY) THEN                              !<-- Take 2D real array data specified for history output
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Insert Solver 2-D Real Data into Field"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          FIELD=ESMF_FieldCreate(grid         =GRID                     &  !<-- The ESMF grid
                                ,farray       =VARS(N)%R2D              &  !<-- The 2D real array being inserted into history Bundle
                                ,datacopyflag =COPYFLAG                 &
                                ,totalUWidth  =(/IHALO,JHALO/)          &
                                ,totalLWidth  =(/IHALO,JHALO/)          &
                                ,name         =VARS(N)%VBL_NAME         &  !<-- Name of the 2D real array
                                ,indexFlag    =ESMF_INDEX_GLOBAL        &
                                ,rc           =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Insert Solver 2-D Real Field into History Bundles"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_FieldBundleAdd(            HISTORY_BUNDLE           &  !<-- The Write component output history Bundle
                                  ,            (/FIELD/)       &  !<-- ESMF Field holding the 2D real array
                                  ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          END IF
          IF (VARS(N)%RESTART) THEN                                        !<-- Take 2D real array data specified for restart output
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Insert Solver 2-D Real Data into Field"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          FIELD=ESMF_FieldCreate(grid         =GRID                     &  !<-- The ESMF grid
                                ,farray       =VARS(N)%R2D              &  !<-- The 2D real array being inserted into restart Bundle
                                ,datacopyflag =COPYFLAG                 &
                                ,totalUWidth  =(/IHALO,JHALO/)          &
                                ,totalLWidth  =(/IHALO,JHALO/)          &
                                ,name         =VARS(N)%VBL_NAME         &  !<-- Name of the 2D real array
                                ,indexFlag    =ESMF_INDEX_GLOBAL        &
                                ,rc           =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Insert Solver 2-D Real Field into Restart Bundles"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_FieldBundleAdd(            RESTART_BUNDLE           &  !<-- The Write component output restart Bundle
                                  ,            (/FIELD/)       &  !<-- ESMF Field holding the 2D real array
                                  ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          END IF
        END IF
      END DO
!
!-----------------------------------------------------------------------
!***  The 3-D real arrays.
!***  We are working with 3-D arrays but they are loaded layer by layer
!***  into 2-D Fields.
!-----------------------------------------------------------------------
!
      DO N=1,NUM_VARS
        IF (VARS(N)%TKR == TKR_R3D .and. ASSOCIATED(VARS(N)%R3D)) THEN
          IF (VARS(N)%HISTORY) THEN                                        !<-- Take 3D real array data specified for history output
          NDIM3=UBOUND(VARS(N)%R3D,3)
          LDIM1=LBOUND(VARS(N)%R3D,1)
          UDIM1=UBOUND(VARS(N)%R3D,1)
          LDIM2=LBOUND(VARS(N)%R3D,2)
          UDIM2=UBOUND(VARS(N)%R3D,2)
!
          DO K=1,NDIM3
            WRITE(MODEL_LEVEL,FMT3)K
            VBL_NAME=TRIM(VARS(N)%VBL_NAME)//'_'//MODEL_LEVEL//'_2D'
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Fill 2-D Fields with Each Level of 3-D Data"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            FIELD=ESMF_FieldCreate(grid         =GRID                   &  !<-- The ESMF grid
                                  ,farray       =VARS(N)%R3D(:,:,K)     &  !<-- Level K of 3D real array being inserted into history Bundle
                                  ,datacopyflag =COPYFLAG               &
                                  ,totalUWidth  =(/IHALO,JHALO/)        &
                                  ,totalLWidth  =(/IHALO,JHALO/)        &
                                  ,name         =VBL_NAME               &  !<-- Name of this level of the 3D real array
                                  ,indexFlag    =ESMF_INDEX_GLOBAL      &
                                  ,rc           =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Insert 3-D Data into History Bundle"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_FieldBundleAdd(            HISTORY_BUNDLE         &  !<-- The Write component output history Bundle
                                    ,            (/FIELD/)     &  !<-- ESMF Field holding the 3D real array
                                    ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          ENDDO
!
          END IF
          IF (VARS(N)%RESTART) THEN                                        !<-- Take 3D real array data specified for restart output
          NDIM3=UBOUND(VARS(N)%R3D,3)
          LDIM1=LBOUND(VARS(N)%R3D,1)
          UDIM1=UBOUND(VARS(N)%R3D,1)
          LDIM2=LBOUND(VARS(N)%R3D,2)
          UDIM2=UBOUND(VARS(N)%R3D,2)
!
          DO K=1,NDIM3
            WRITE(MODEL_LEVEL,FMT3)K
            VBL_NAME=TRIM(VARS(N)%VBL_NAME)//'_'//MODEL_LEVEL//'_2D'
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Fill 2-D Fields with Each Level of 3-D Data"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            FIELD=ESMF_FieldCreate(grid         =GRID                   &  !<-- The ESMF grid
                                  ,farray       =VARS(N)%R3D(:,:,K)               &  !<-- Level K of 3D real array being inserted into restart Bundle
                                  ,datacopyflag =COPYFLAG               &
                                  ,totalUWidth  =(/IHALO,JHALO/)        &
                                  ,totalLWidth  =(/IHALO,JHALO/)        &
                                  ,name         =VBL_NAME               &  !<-- Name of this level of the 3D real array
                                  ,indexFlag    =ESMF_INDEX_GLOBAL      &
                                  ,rc           =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Insert 3-D Data into Restart Bundle"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_FieldBundleAdd(            RESTART_BUNDLE         &  !<-- The Write component output restart Bundle
                                    ,            (/FIELD/)     &  !<-- ESMF Field holding the 3D real array
                                    ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          ENDDO
!
          END IF
        END IF
      END DO
!
!-----------------------------------------------------------------------
!***  The 4-D real arrays.
!***  We are working with 4-D arrays but they are loaded layer by layer
!***  into 2-D Fields.
!-----------------------------------------------------------------------
!
!!!!!!!!! FIX this later
!!!!!!!!! FIX this later
      INDX_Q2 = 3
!!!!!!!!! FIX this later
!!!!!!!!! FIX this later

      DO N=1,NUM_VARS
        IF (VARS(N)%TKR == TKR_R4D) THEN
          IF (VARS(N)%HISTORY) THEN                                        !<-- Take 4D real array data specified for history output
          LDIM1=LBOUND(VARS(N)%R4D,1)
          UDIM1=UBOUND(VARS(N)%R4D,1)
          LDIM2=LBOUND(VARS(N)%R4D,2)
          UDIM2=UBOUND(VARS(N)%R4D,2)
          LDIM3=LBOUND(VARS(N)%R4D,3)
          UDIM3=UBOUND(VARS(N)%R4D,3)
          LDIM4=LBOUND(VARS(N)%R4D,4)
          UDIM4=UBOUND(VARS(N)%R4D,4)
!
          IF(TRIM(VARS(N)%VBL_NAME)=='TRACERS_PREV' .OR. &
             TRIM(VARS(N)%VBL_NAME)=='TRACERS') THEN
            LDIM4=INDX_Q2+1                                                !<-- TRACERS bounds:      INDX_Q2+1 - UDIM4
          ENDIF

          DO M=LDIM4,UDIM4                                                 !<-- Loop through the tracers (skip unallocated pointers)
          DO K=LDIM3,UDIM3                                                 !<-- Loop through the levels of the array
            WRITE(TRACERS_KIND,FMT3)M
            WRITE(MODEL_LEVEL,FMT3)K
!
            VBL_NAME=TRIM(VARS(N)%VBL_NAME)//'_'//TRACERS_KIND//'_'//MODEL_LEVEL//'_2D'
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Fill 2-D Fields with Each Level of 4-D Data"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            FIELD=ESMF_FieldCreate(grid         =GRID                   &  !<-- The ESMF grid
                                  ,farray       =VARS(N)%R4D(:,:,K,M)   &  !<-- Level K of 4D real array being inserted into history Bundle
                                  ,datacopyflag =COPYFLAG               &
                                  ,totalUWidth  =(/IHALO,JHALO/)        &
                                  ,totalLWidth  =(/IHALO,JHALO/)        &
                                  ,name         =VBL_NAME               &  !<-- Name of this level of the 4D real array
                                  ,indexFlag    =ESMF_INDEX_GLOBAL      &
                                  ,rc           =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Insert 4-D Data into History Bundle"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_FieldBundleAdd(            HISTORY_BUNDLE         &  !<-- The Write component output history Bundle
                                    ,            (/FIELD/)     &  !<-- ESMF Field holding the 4D real array
                                    ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          ENDDO
          ENDDO
!
          END IF
          IF (VARS(N)%RESTART) THEN                                        !<-- Take 4D real array data specified for restart output
          LDIM1=LBOUND(VARS(N)%R4D,1)
          UDIM1=UBOUND(VARS(N)%R4D,1)
          LDIM2=LBOUND(VARS(N)%R4D,2)
          UDIM2=UBOUND(VARS(N)%R4D,2)
          LDIM3=LBOUND(VARS(N)%R4D,3)
          UDIM3=UBOUND(VARS(N)%R4D,3)
          LDIM4=LBOUND(VARS(N)%R4D,4)
          UDIM4=UBOUND(VARS(N)%R4D,4)
!
          IF( TRIM(VARS(N)%VBL_NAME)=='TRACERS') THEN
            LDIM4=INDX_Q2+1                                                !<-- TRACERS bounds:      INDX_Q2+1 - UDIM4
          ENDIF
!
          DO M=LDIM4,UDIM4                                                 !<-- Loop through the tracers (skip unallocated pointers)
          DO K=LDIM3,UDIM3                                                 !<-- Loop through the levels of the array
            WRITE(TRACERS_KIND,FMT3)M
            WRITE(MODEL_LEVEL,FMT3)K
!
            VBL_NAME=TRIM(VARS(N)%VBL_NAME)//'_'//TRACERS_KIND//'_'//MODEL_LEVEL//'_2D'
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Fill 2-D Fields with Each Level of 4-D Data"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            FIELD=ESMF_FieldCreate(grid         =GRID                   &  !<-- The ESMF grid
                                  ,farray       =VARS(N)%R4D(:,:,K,M)               &  !<-- Level K of 4D real array being inserted into restart Bundle
                                  ,datacopyflag =COPYFLAG               &
                                  ,totalUWidth  =(/IHALO,JHALO/)        &
                                  ,totalLWidth  =(/IHALO,JHALO/)        &
                                  ,name         =VBL_NAME               &  !<-- Name of this level of the 4D real array
                                  ,indexFlag    =ESMF_INDEX_GLOBAL      &
                                  ,rc           =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Insert 4-D Data into Restart Bundle"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_FieldBundleAdd(            RESTART_BUNDLE         &  !<-- The Write component output restart Bundle
                                    ,            (/FIELD/)     &  !<-- ESMF Field holding the 4D real array
                                    ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          ENDDO
          ENDDO
!
          END IF
        END IF
      END DO
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE PUT_VARS_IN_BUNDLES

!-------------------------------------------------------------------------------
!###############################################################################
!-------------------------------------------------------------------------------

      END MODULE MODULE_VARS_STATE
