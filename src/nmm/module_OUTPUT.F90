!-----------------------------------------------------------------------
!
      MODULE MODULE_OUTPUT
!
!-----------------------------------------------------------------------
!***  Insert quantities from the Solver internal state into the
!***  Write import state for output.
!-----------------------------------------------------------------------
!
      USE ESMF
      USE MODULE_KINDS
      USE MODULE_SOLVER_INTERNAL_STATE,ONLY: SOLVER_INTERNAL_STATE
      USE MODULE_ERROR_MSG,ONLY: ERR_MSG,MESSAGE_CHECK
      USE MODULE_VARS
      USE MODULE_VARS_STATE
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      PRIVATE
!
      PUBLIC :: POINT_OUTPUT
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      SUBROUTINE POINT_OUTPUT(GRID,INT_STATE,IMP_STATE_WRITE)
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  This routine takes the user's selections for output quantities,
!***  points at them, and inserts those pointers into the import state
!***  of the Write components.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument variables
!------------------------
!
      TYPE(ESMF_Grid) ,INTENT(IN) :: GRID                                  !<-- The ESMF Grid
!
      TYPE(SOLVER_INTERNAL_STATE),POINTER,INTENT(INOUT) :: INT_STATE       !<-- The Solver internal state
!
      TYPE(ESMF_State),INTENT(INOUT) :: IMP_STATE_WRITE                    !<-- Import state for the Write components
!
!---------------------
!***  Local variables
!---------------------
!
      INTEGER :: IHALO,JHALO
!
      INTEGER :: K,LENGTH,MYPE                                          &
                ,N,NDIM3,NFIND,NUM_2D_FIELDS,NV                         &
                ,RC,RC_DYN_OUT
!
      INTEGER :: LDIM1,LDIM2                                            &
                ,UDIM1,UDIM2
!
      INTEGER :: ITWO=2
!
      REAL(KIND=KFPT),DIMENSION(:,:),POINTER :: TEMP_R2D
!
      CHARACTER(2)                :: MODEL_LEVEL
      CHARACTER(6)                :: FMT='(I2.2)'
      CHARACTER(ESMF_MAXSTR)      :: VBL_NAME
!
      TYPE(ESMF_FieldBundle) :: HISTORY_BUNDLE                          &
                               ,RESTART_BUNDLE
!
      TYPE(ESMF_Field)            :: FIELD
!
      TYPE(ESMF_DataCopy_Flag)    :: COPYFLAG=ESMF_DATACOPY_REFERENCE
!     TYPE(ESMF_DataCopy_Flag)    :: COPYFLAG=ESMF_DATA_COPY
!
!-----------------------------------------------------------------------
!
      MYPE=int_state%MYPE
!
!-----------------------------------------------------------------------
!***  Create an ESMF Bundle that will hold history output data
!***  and nothing else.  This will serve to isolate the output
!***  data from everything else inside the Write component's
!***  import state.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create History Data Bundle"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      HISTORY_BUNDLE=ESMF_FieldBundleCreate(name='History Bundle'       &  !<-- The Bundle's name
                                           ,rc  =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_DYN_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Create an ESMF Bundle that will hold restart data
!***  and nothing else.  This will serve to isolate the restart
!***  data from everything else inside the Write component's
!***  import state.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create Restart Data Bundle"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      RESTART_BUNDLE=ESMF_FieldBundleCreate(name='Restart Bundle'       &  !<-- The Bundle's name
                                           ,rc  =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_DYN_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  First add the local subdomain limits to the Write component's
!***  import state as Attributes along with the global/regional mode.
!***  This information is needed for quilting the local domain data
!***  into full domain fields.
!***  The local domain limits go directly into the Write component's
!***  import state to keep them separate from the history data that
!***  will be inserted into a Bundle.
!
!***  Do the same with the number of fcst tasks (INPESxJNPES) since
!***  the Write component also needs that information as well as
!***  the halo depths.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Add Local Subdomain Limits to the Write Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

      CALL ESMF_AttributeSet(state    =IMP_STATE_WRITE                  &  !<-- The Write component import state
                            ,name     ='LOCAL_ISTART'                   &  !<-- Name of the integer array
                            ,itemCount=int_state%NUM_PES                &  !<-- Length of array being inserted into the import state
                            ,valueList=int_state%LOCAL_ISTART           &  !<-- The array being inserted into the import state
                            ,rc       =RC)
!
      CALL ESMF_AttributeSet(state    =IMP_STATE_WRITE                  &  !<-- The Write component import state
                            ,name     ='LOCAL_IEND'                     &  !<-- Name of the integer array
                            ,itemCount=int_state%NUM_PES                &  !<-- Length of array being inserted into the import state
                            ,valueList=int_state%LOCAL_IEND             &  !<-- The array being inserted into the import state
                            ,rc       =RC)
!
      CALL ESMF_AttributeSet(state    =IMP_STATE_WRITE                  &  !<-- The Write component import state
                            ,name     ='LOCAL_JSTART'                   &  !<-- Name of the integer array
                            ,itemCount=int_state%NUM_PES                &  !<-- Length of array being inserted into the import state
                            ,valueList=int_state%LOCAL_JSTART           &  !<-- The array being inserted into the import state
                            ,rc       =RC)
!
      CALL ESMF_AttributeSet(state    =IMP_STATE_WRITE                  &  !<-- The Write component import state
                            ,name     ='LOCAL_JEND'                     &  !<-- Name of the integer array
                            ,itemCount=int_state%NUM_PES                &  !<-- Length of array being inserted into the import state
                            ,valueList=int_state%LOCAL_JEND             &  !<-- The array being inserted into the import state
                            ,rc       =RC)
!
      CALL ESMF_AttributeSet(state=IMP_STATE_WRITE                      &  !<-- The Write component import state
                            ,name ='INPES'                              &  !<-- Name of the integer scalar
                            ,value=int_state%INPES                      &  !<-- The value being inserted into the import state
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=IMP_STATE_WRITE                      &  !<-- The Write component import state
                            ,name ='JNPES'                              &  !<-- Name of the integer scalar
                            ,value=int_state%JNPES                      &  !<-- The value being inserted into the import state
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=IMP_STATE_WRITE                      &  !<-- The Write component import state
                            ,name ='IHALO'                              &  !<-- Name of the integer scalar
                            ,value=int_state%IHALO                      &  !<-- The value being inserted into the import state
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=IMP_STATE_WRITE                      &  !<-- The Write component import state
                            ,name ='JHALO'                              &  !<-- Name of the integer scalar
                            ,value=int_state%JHALO                      &  !<-- The value being inserted into the import state
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=IMP_STATE_WRITE                      &  !<-- The Write component import state
                            ,name ='WRITE_TASKS_PER_GROUP'              &  !<-- Name of the integer scalar
                            ,value=int_state%WRITE_TASKS_PER_GROUP      &  !<-- The value being inserted into the import state
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=IMP_STATE_WRITE                      &  !<-- The Write component import state
                            ,name ='WRITE_GROUPS'                       &  !<-- Name of the integer scalar
                            ,value=int_state%WRITE_GROUPS               &  !<-- The value being inserted into the import state
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=IMP_STATE_WRITE                      &  !<-- The Write component import state
                            ,name ='LNSH'                               &  !<-- Name of the integer scalar
                            ,value=int_state%LNSH                       &  !<-- The value being inserted into the import state
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=IMP_STATE_WRITE                      &  !<-- The Write component import state
                            ,name ='LNSV'                               &  !<-- Name of the integer scalar
                            ,value=int_state%LNSV                       &  !<-- The value being inserted into the import state
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=IMP_STATE_WRITE                      &  !<-- The Write component import state
                            ,name ='NVARS_BC_2D_H'                      &  !<-- Name of the integer scalar
                            ,value=int_state%NVARS_BC_2D_H              &  !<-- The value being inserted into the import state
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=IMP_STATE_WRITE                      &  !<-- The Write component import state
                            ,name ='NVARS_BC_3D_H'                      &  !<-- Name of the integer scalar
                            ,value=int_state%NVARS_BC_3D_H              &  !<-- The value being inserted into the import state
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=IMP_STATE_WRITE                      &  !<-- The Write component import state
                            ,name ='NVARS_BC_4D_H'                      &  !<-- Name of the integer scalar
                            ,value=int_state%NVARS_BC_4D_H              &  !<-- The value being inserted into the import state
                            ,rc   =RC)
!
      IF(int_state%NVARS_BC_4D_H>0)THEN
        CALL ESMF_AttributeSet(state    =IMP_STATE_WRITE                &  !<-- The Write component import state
                              ,name     ='LBND_4D'                      &  !<-- Name of the integer scalar
                              ,itemCount=int_state%NVARS_BC_4D_H        &  !<-- Length of array being inserted into the import state
                              ,valuelist=int_state%LBND_4D              &  !<-- The array being inserted into the import state
                              ,rc   =RC)
        CALL ESMF_AttributeSet(state    =IMP_STATE_WRITE                &  !<-- The Write component import state
                              ,name     ='UBND_4D'                      &  !<-- Name of the integer scalar
                              ,itemCount=int_state%NVARS_BC_4D_H        &  !<-- Length of array being inserted into the import state
                              ,valuelist=int_state%UBND_4D              &  !<-- The array being inserted into the import state
                              ,rc   =RC)
      ENDIF
!
      CALL ESMF_AttributeSet(state=IMP_STATE_WRITE                      &  !<-- The Write component import state
                            ,name ='NVARS_BC_2D_V'                      &  !<-- Name of the integer scalar
                            ,value=int_state%NVARS_BC_2D_V              &  !<-- The value being inserted into the import state
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=IMP_STATE_WRITE                      &  !<-- The Write component import state
                            ,name ='NVARS_BC_3D_V'                      &  !<-- Name of the integer scalar
                            ,value=int_state%NVARS_BC_3D_V              &  !<-- The value being inserted into the import state
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=IMP_STATE_WRITE                      &  !<-- The Write component import state
                            ,name ='NLEV_H'                             &  !<-- Name of the integer scalar
                            ,value=int_state%NLEV_H                     &  !<-- The value being inserted into the import state
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=IMP_STATE_WRITE                      &  !<-- The Write component import state
                            ,name ='NLEV_V'                             &  !<-- Name of the integer scalar
                            ,value=int_state%NLEV_V                     &  !<-- The value being inserted into the import state
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state    =IMP_STATE_WRITE                  &  !<-- The Write component import state
                            ,name     ='LOCAL_JEND'                     &  !<-- Name of the integer array
                            ,itemCount=int_state%NUM_PES                &  !<-- Length of array being inserted into the import state
                            ,valueList=int_state%LOCAL_JEND             &  !<-- The array being inserted into the import state
                            ,rc       =RC)
!
      CALL ESMF_AttributeSet(state=IMP_STATE_WRITE                      &  !<-- The Write component import state
                            ,name ='IDS'                                &  !<-- Name of the integer scalar
                            ,value=int_state%IDS                        &  !<-- The value being inserted into the import state
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=IMP_STATE_WRITE                      &  !<-- The Write component import state
                            ,name ='IDE'                                &  !<-- Name of the integer scalar
                            ,value=int_state%IDE                        &  !<-- The value being inserted into the import state
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=IMP_STATE_WRITE                      &  !<-- The Write component import state
                            ,name ='JDS'                                &  !<-- Name of the integer scalar
                            ,value=int_state%JDS                        &  !<-- The value being inserted into the import state
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=IMP_STATE_WRITE                      &  !<-- The Write component import state
                            ,name ='JDE'                                &  !<-- Name of the integer scalar
                            ,value=int_state%JDE                        &  !<-- The value being inserted into the import state
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_DYN_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  The following logical variables are to be part of the
!***  history output therefore place them into the history Bundle.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Insert Global and Run Logicals into History Bundle"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

      CALL ESMF_AttributeSet(FIELDBUNDLE=HISTORY_BUNDLE                 &  !<-- The Write component output history Bundle
                            ,name       ='GLOBAL'                       &  !<-- Name of the logical
                            ,value      =int_state%GLOBAL               &  !<-- The logical being inserted into the Bundle
                            ,rc         =RC)
!
      CALL ESMF_AttributeSet(FIELDBUNDLE=HISTORY_BUNDLE                 &  !<-- The Write component output history Bundle
                            ,name       ='RUN'                          &  !<-- Name of the logical
                            ,value      =int_state%RUN                  &  !<-- The logical being inserted into the Bundle
                            ,rc         =RC)
!
      CALL ESMF_AttributeSet(FIELDBUNDLE=HISTORY_BUNDLE                 &  !<-- The Write component output history Bundle
                            ,name       ='ADIABATIC'                    &  !<-- Name of the logical
                            ,value      =int_state%ADIABATIC            &  !<-- The logical being inserted into the Bundle
                            ,rc         =RC)

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_DYN_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  The following logical variables are to be part of the
!***  restart output therefore place them into the restart Bundle.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Insert Global and Run Logicals into Restart Bundle"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

      CALL ESMF_AttributeSet(FIELDBUNDLE=RESTART_BUNDLE                 &  !<-- The Write component restart Bundle
                            ,name       ='GLOBAL'                       &  !<-- Name of the logical
                            ,value      =int_state%GLOBAL               &  !<-- The logical being inserted into the Bundle
                            ,rc         =RC)
!
      CALL ESMF_AttributeSet(FIELDBUNDLE=RESTART_BUNDLE                 &  !<-- The Write component restart Bundle
                            ,name       ='RUN'                          &  !<-- Name of the logical
                            ,value      =int_state%RUN                  &  !<-- The logical being inserted into the Bundle
                            ,rc         =RC)
!
      CALL ESMF_AttributeSet(FIELDBUNDLE=RESTART_BUNDLE                 &  !<-- The Write component restart Bundle
                            ,name       ='ADIABATIC'                    &  !<-- Name of the logical
                            ,value      =int_state%ADIABATIC            &  !<-- The logical being inserted into the Bundle
                            ,rc         =RC)

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_DYN_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Now insert into the Write components' import state the pointers
!***  of only those quantities that are specified by the user for
!***  history and restart output.  The data is placed into an ESMF
!***  Bundle which itself will be placed into the import state at
!***  the end of the routine.
!-----------------------------------------------------------------------
!
      CALL PUT_VARS_IN_BUNDLES(int_state%VARS                           &
                              ,int_state%NUM_VARS                       &
                              ,GRID                                     &
                              ,HISTORY_BUNDLE                           &
                              ,RESTART_BUNDLE)
!
!-----------------------------------------------------------------------
!***  Load the two output Bundles into the Solver's internal state
!***  array needed to add them to the Write component's import state.
!-----------------------------------------------------------------------
!
      int_state%BUNDLE_ARRAY(1)=HISTORY_BUNDLE
      int_state%BUNDLE_ARRAY(2)=RESTART_BUNDLE
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Solver: Insert Bundle Array into the Write Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateAddReplace(state          =IMP_STATE_WRITE                     &  !<-- The Write component's import state
                        ,fieldbundlelist=(/int_state%BUNDLE_ARRAY/) &  !<-- Array holding the History/Restart Bundles
                        ,rc             =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_DYN_OUT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE POINT_OUTPUT
!
!-----------------------------------------------------------------------
!
      END MODULE MODULE_OUTPUT
!
!-----------------------------------------------------------------------
