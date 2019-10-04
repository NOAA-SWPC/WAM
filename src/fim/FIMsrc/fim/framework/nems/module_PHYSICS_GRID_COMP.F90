!JR copied from fimlatest
!TODO:  DRY out all of this code.  Initial NEMS was not DRY.  We can be.  
!-----------------------------------------------------------------------
!
      MODULE MODULE_PHYSICS_GRID_COMP
!
!-----------------------------------------------------------------------
!
!***  THIS MODULE HOLDS THE PHYSICS REGISTER, INIT, RUN, AND FINALIZE 
!***  ROUTINES.  THEY ARE CALLED FROM THE MAIN GRIDDED COMPONENT
!***  (ATM INITIALIZE CALLS PHYSICS INITIALIZE, ETC.) 
!***  IN MODULE_ATM_GRID_COMP.F.
!
!-----------------------------------------------------------------------
!
      USE ESMF_MOD
!
      USE MODULE_ERR_MSG,ONLY: ERR_MSG,MESSAGE_CHECK
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      PRIVATE
!
      PUBLIC :: PHY_REGISTER
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE PHY_REGISTER(GRID_COMP,RC_REG)
!
!-----------------------------------------------------------------------
!***  REGISTER THE PHYSICS COMPONENT'S INITIALIZE, RUN, AND FINALIZE
!***  ROUTINES.
!-----------------------------------------------------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: GRID_COMP                      !<-- The Physics Gridded Component
!
      INTEGER,INTENT(OUT) :: RC_REG                                       !<-- Return code for Register
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      INTEGER :: RC=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC_REG=ESMF_SUCCESS                                                 !<-- Initialize error signal variable
                                                                                                                                              
!-----------------------------------------------------------------------
!***  REGISTER THE PHYSICS INITIALIZE SUBROUTINE.  SINCE IT IS JUST ONE
!***  SUBROUTINE, USE ESMF_SINGLEPHASE.  THE SECOND ARGUMENT IS
!***  A PRE-DEFINED SUBROUTINE TYPE, SUCH AS ESMF_SETINIT, ESMF_SETRUN,
!***  OR ESMF_SETFINAL.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set Entry Point for Physics Initialize"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompSetEntryPoint(GRID_COMP                         &  !<-- Physics gridcomp
                                     ,ESMF_SETINIT                      &  !<-- Subroutine type
                                     ,PHY_INITIALIZE                    &  !<-- User's subroutine name
                                     ,ESMF_SINGLEPHASE                  &  !<-- Phase
                                     ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  REGISTER THE PHYSICS RUN SUBROUTINE.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set Entry Point for Physics Run"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompSetEntryPoint(GRID_COMP                         &  !<-- Physics gridcomp
                                     ,ESMF_SETRUN                       &  !<-- Subroutine type
                                     ,PHY_RUN                           &  !<-- User's subroutine name
                                     ,ESMF_SINGLEPHASE                  &  !<-- Phase
                                     ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  REGISTER THE PHYSICS FINALIZE SUBROUTINE.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set Entry Point for Physics Finalize"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompSetEntryPoint(GRID_COMP                         &  !<-- Physics gridcomp
                                     ,ESMF_SETFINAL                     &  !<-- Subroutine type
                                     ,PHY_FINALIZE                      &  !<-- User's subroutine name
                                     ,ESMF_SINGLEPHASE                  &  !<-- Phase
                                     ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  CHECK THE ERROR SIGNAL VARIABLE.
!-----------------------------------------------------------------------
!
      IF(RC_REG==ESMF_SUCCESS)THEN
!       WRITE(0,*)' PHY_REGISTER SUCCEEDED'
      ELSE
        WRITE(0,*)' PHY_REGISTER FAILED RC_REG=',RC_REG
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE PHY_REGISTER
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE PHY_INITIALIZE(GRID_COMP                               &
                               ,IMP_STATE,EXP_STATE                     &
                               ,CLOCK_ATM                               &
                               ,RC_INIT)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  SET UP THE MODEL PHYSICS.
!-----------------------------------------------------------------------
!
      USE module_fim_phy_init ,only: PHY_INITIALIZE_FIM => phy_init
      USE gfs_physics_internal_state_mod, only:  &
        gfs_physics_internal_state,              &
        gis_phy,                                 &
        WRAP_INTERNAL_STATE => gfs_phy_wrap
!
!-----------------------------------------------------------------------
!***  ARGUMENT VARIABLES.
!-----------------------------------------------------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: GRID_COMP  !<-- The Physics gridded component
      TYPE(ESMF_State),   INTENT(INOUT) :: IMP_STATE  !<-- The Physics Initialize step's import state
      TYPE(ESMF_State),   INTENT(INOUT) :: EXP_STATE  !<-- The Physics Initialize step's export state
      TYPE(ESMF_Clock),   INTENT(IN)    :: CLOCK_ATM  !<-- The ATM's ESMF Clock
!
      INTEGER,            INTENT(OUT)   :: RC_INIT
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!***  WRAP_INTERNAL_STATE IS DEFINED IN THE INTERNAL STATE MODULE.
!-----------------------------------------------------------------------
!
      TYPE(WRAP_INTERNAL_STATE)    :: WRAP                               !<-- This wrap is a derived type which contains
                                                                         !    only a pointer to the internal state.  It is needed
                                                                         !    for using different architectures or compilers.
!
!      TYPE(ESMF_Field) :: TMP_FIELD
      TYPE(ESMF_Grid) :: GRID
      TYPE(ESMF_DistGrid) :: DISTGRID
      TYPE(ESMF_Array) :: TMP_ARRAY
      INTEGER :: RC
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  INITIALIZE THE ERROR SIGNAL VARIABLES.
!-----------------------------------------------------------------------
!
      RC     =ESMF_SUCCESS
      RC_INIT=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  INITIALIZE THE PHYSICS SCHEMES. 
!-----------------------------------------------------------------------
!
! Allocate internal state and set up initial values for some fields.  
!
      CALL PHY_INITIALIZE_FIM
!
!-----------------------------------------------------------------------
!***  ATTACH THE INTERNAL STATE TO THE PHYSICS GRIDDED COMPONENT.
!-----------------------------------------------------------------------
!
      WRAP%INT_STATE=>gis_phy
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK= &
        "Attach Physics Internal State to the Gridded Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompSetInternalState(GRID_COMP                      &  !<-- Physics gridcomp
                                        ,WRAP                           &  !<-- Data pointer to internal state
                                        ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!

      MESSAGE_CHECK="PHY:  Extract GRID from GRID_COMP"
      call esmf_gridcompget(GRID_COMP, grid = GRID, rc = RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)

      call esmf_gridvalidate(grid=GRID, rc=rc)
      CALL ERR_MSG(RC,'PHY:  validate GRID',RC_INIT)

      MESSAGE_CHECK="PHY:  Extract DISTGRID from GRID"
      CALL ESMF_GridGet(grid=GRID, distgrid=DISTGRID, rc=RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)

      call esmf_distgridvalidate(distgrid=DISTGRID, rc=rc)
      CALL ERR_MSG(RC,'PHY:  validate DISTGRID',RC_INIT)

! Set flags to enable import/export of each field, hard-coded for the moment
!TODO:  Read these from a config file as in GFS, adapting 
!TODO:  gfs_physics_getcf_mod.f
      gis_phy%esmf_sta_list%idate1_import = 0
      gis_phy%esmf_sta_list%idate1_export = 0
      gis_phy%esmf_sta_list%z_import      = 0
      gis_phy%esmf_sta_list%z_export      = 0
      gis_phy%esmf_sta_list%ps_import     = 1
      gis_phy%esmf_sta_list%ps_export     = 1
      gis_phy%esmf_sta_list%temp_import   = 1
      gis_phy%esmf_sta_list%temp_export   = 1
      gis_phy%esmf_sta_list%u_import      = 1
      gis_phy%esmf_sta_list%u_export      = 1
      gis_phy%esmf_sta_list%v_import      = 1
      gis_phy%esmf_sta_list%v_export      = 1
      gis_phy%esmf_sta_list%q_import      = 1
      gis_phy%esmf_sta_list%q_export      = 1
      gis_phy%esmf_sta_list%oz_import     = 1
      gis_phy%esmf_sta_list%oz_export     = 1
      gis_phy%esmf_sta_list%cld_import    = 1
      gis_phy%esmf_sta_list%cld_export    = 1
      gis_phy%esmf_sta_list%p_import      = 1
      gis_phy%esmf_sta_list%p_export      = 1
      gis_phy%esmf_sta_list%dp_import     = 1
      gis_phy%esmf_sta_list%dp_export     = 1
      gis_phy%esmf_sta_list%dpdt_import   = 1
      gis_phy%esmf_sta_list%dpdt_export   = 1

!-----------------------------------------------------------------------
!***  Attach gfs fields in the internal state 
!***  to the esmf import and export states.  
!TBH:  I use GFS naming conventions, *not* NMMB conventions.  As 
!TBH:  of NEMS r3038 they do indeed differ, by case at least!  
!TBH:  Creation of unique ESMF_Field objects for import and export 
!TBH:  states should require little additional memory since the pointers 
!TBH:  to Fortran arrays are shared.  This approach makes object deletion 
!TBH:  easier.  
!-------------------------------------------------------

      MESSAGE_CHECK= &
        "initial internal state to esmf import and export states"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)

      IF (gis_phy%esmf_sta_list%idate1_import == 1) THEN
        WRITE(0,*)' PHY_INITIALIZE import of idate1 not supported'
        RC_INIT = ESMF_FAILURE
        RETURN
      ENDIF
      IF (gis_phy%esmf_sta_list%idate1_export == 1) THEN
        WRITE(0,*)' PHY_INITIALIZE export of idate1 not supported'
        RC_INIT = ESMF_FAILURE
        RETURN
      ENDIF
      IF (gis_phy%esmf_sta_list%z_import == 1) THEN
        WRITE(0,*)' PHY_INITIALIZE import of z not supported'
        RC_INIT = ESMF_FAILURE
        RETURN
      ENDIF
      IF (gis_phy%esmf_sta_list%z_export == 1) THEN
        WRITE(0,*)' PHY_INITIALIZE export of z not supported'
        RC_INIT = ESMF_FAILURE
        RETURN
      ENDIF

!TODO:  Need to add gridToFieldMap to ESMF_FieldCreate() to address 
!TODO:  differences between 2D and 3D arrays.  At present this is 
!TODO:  irrelevant since we do not use ESMF to do any re-grid or 
!TODO:  re-dist operations.  This must be fixed before we use these 
!TODO:  ESMF features.  

      IF (gis_phy%esmf_sta_list%ps_import == 1) THEN
        MESSAGE_CHECK="Create ps array for import state"
!TBH:  Note that the following call to ESMF_FieldCreate() yields the 
!TBH:  stunningly informative error code 540 which maps to string 
!TBH:  "Not valid" in ESMC_ErrMsgs.C.  Backed off to ESMF_ArrayCreate().  
!TODO:  Switch back to ESMF_FieldCreate() since future NEMS will use 
!TODO:  ESMF_Fields.  
!        TMP_FIELD=ESMF_FieldCreate(grid       =GRID                     &
!                                  ,farray     =gis_phy%ps          &
!                                  ,distgridToArrayMap=(/1/)             &
!                                  ,name       ='ps'                     &
!                                  ,rc         =RC)
        ! create the ESMF_Array
        TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                  ,farray     =gis_phy%ps          &
                                  ,distgridToArrayMap=(/1/)             &
                                  ,indexflag=ESMF_INDEX_GLOBAL          &
                                  ,name       ='ps'                     &
                                  ,rc         =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
        ! attach array to state
        MESSAGE_CHECK="Add ps array to import state"
        CALL ESMF_StateAdd(state=IMP_STATE                              &
                          ,array=TMP_ARRAY                              &
                          ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ENDIF
      IF (gis_phy%esmf_sta_list%ps_export == 1) THEN
        MESSAGE_CHECK="Create ps array for export state"
        ! create the ESMF_Array
        TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                  ,farray     =gis_phy%ps          &
                                  ,distgridToArrayMap=(/1/)             &
                                  ,indexflag=ESMF_INDEX_GLOBAL          &
                                  ,name       ='ps'                     &
                                  ,rc         =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
        ! attach array to state
        MESSAGE_CHECK="Add ps array to export state"
        CALL ESMF_StateAdd(state=EXP_STATE                              &
                          ,array=TMP_ARRAY                              &
                          ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ENDIF

      IF (gis_phy%esmf_sta_list%temp_import == 1) THEN
        MESSAGE_CHECK="Create t array for import state"
        ! create the ESMF_Array
        TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                  ,farray     =gis_phy%t           &
                                  ,distgridToArrayMap=(/1/)             &
                                  ,indexflag=ESMF_INDEX_GLOBAL          &
                                  ,name       ='t'                      &
                                  ,rc         =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
        ! attach array to state
        MESSAGE_CHECK="Add t array to import state"
        CALL ESMF_StateAdd(state=IMP_STATE                              &
                          ,array=TMP_ARRAY                              &
                          ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ENDIF
      IF (gis_phy%esmf_sta_list%temp_export == 1) THEN
        MESSAGE_CHECK="Create t array for export state"
        ! create the ESMF_Array
        TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                  ,farray     =gis_phy%t           &
                                  ,distgridToArrayMap=(/1/)             &
                                  ,indexflag=ESMF_INDEX_GLOBAL          &
                                  ,name       ='t'                      &
                                  ,rc         =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
        ! attach array to state
        MESSAGE_CHECK="Add t array to export state"
        CALL ESMF_StateAdd(state=EXP_STATE                              &
                          ,array=TMP_ARRAY                              &
                          ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ENDIF

      IF (gis_phy%esmf_sta_list%u_import == 1) THEN
        MESSAGE_CHECK="Create u array for import state"
        ! create the ESMF_Array
        TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                  ,farray     =gis_phy%u           &
                                  ,distgridToArrayMap=(/1/)             &
                                  ,indexflag=ESMF_INDEX_GLOBAL          &
                                  ,name       ='u'                      &
                                  ,rc         =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
        ! attach array to state
        MESSAGE_CHECK="Add u array to import state"
        CALL ESMF_StateAdd(state=IMP_STATE                              &
                          ,array=TMP_ARRAY                              &
                          ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ENDIF
      IF (gis_phy%esmf_sta_list%u_export == 1) THEN
        MESSAGE_CHECK="Create u array for export state"
        ! create the ESMF_Array
        TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                  ,farray     =gis_phy%u           &
                                  ,distgridToArrayMap=(/1/)             &
                                  ,indexflag=ESMF_INDEX_GLOBAL          &
                                  ,name       ='u'                      &
                                  ,rc         =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
        ! attach array to state
        MESSAGE_CHECK="Add u array to export state"
        CALL ESMF_StateAdd(state=EXP_STATE                              &
                          ,array=TMP_ARRAY                              &
                          ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ENDIF

      IF (gis_phy%esmf_sta_list%v_import == 1) THEN
        MESSAGE_CHECK="Create v array for import state"
        ! create the ESMF_Array
        TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                  ,farray     =gis_phy%v           &
                                  ,distgridToArrayMap=(/1/)             &
                                  ,indexflag=ESMF_INDEX_GLOBAL          &
                                  ,name       ='v'                      &
                                  ,rc         =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
        ! attach array to state
        MESSAGE_CHECK="Add v array to import state"
        CALL ESMF_StateAdd(state=IMP_STATE                              &
                          ,array=TMP_ARRAY                              &
                          ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ENDIF
      IF (gis_phy%esmf_sta_list%v_export == 1) THEN
        MESSAGE_CHECK="Create v array for export state"
        ! create the ESMF_Array
        TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                  ,farray     =gis_phy%v           &
                                  ,distgridToArrayMap=(/1/)             &
                                  ,indexflag=ESMF_INDEX_GLOBAL          &
                                  ,name       ='v'                      &
                                  ,rc         =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
        ! attach array to state
        MESSAGE_CHECK="Add v array to export state"
        CALL ESMF_StateAdd(state=EXP_STATE                              &
                          ,array=TMP_ARRAY                              &
                          ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ENDIF

      IF (gis_phy%esmf_sta_list%q_import == 1) THEN
        MESSAGE_CHECK="Create q array for import state"
        ! create the ESMF_Array
        TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                  ,farray     =gis_phy%q           &
                                  ,distgridToArrayMap=(/1/)             &
                                  ,indexflag=ESMF_INDEX_GLOBAL          &
                                  ,name       ='shum'                   &
                                  ,rc         =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
        ! attach array to state
        MESSAGE_CHECK="Add q array to import state"
        CALL ESMF_StateAdd(state=IMP_STATE                              &
                          ,array=TMP_ARRAY                              &
                          ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ENDIF
      IF (gis_phy%esmf_sta_list%q_export == 1) THEN
        MESSAGE_CHECK="Create q array for export state"
        ! create the ESMF_Array
        TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                  ,farray     =gis_phy%q           &
                                  ,distgridToArrayMap=(/1/)             &
                                  ,indexflag=ESMF_INDEX_GLOBAL          &
                                  ,name       ='shum'                   &
                                  ,rc         =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
        ! attach array to state
        MESSAGE_CHECK="Add q array to export state"
        CALL ESMF_StateAdd(state=EXP_STATE                              &
                          ,array=TMP_ARRAY                              &
                          ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ENDIF

      IF (gis_phy%esmf_sta_list%oz_import == 1) THEN
        MESSAGE_CHECK="Create oz array for import state"
        ! create the ESMF_Array
        TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                  ,farray     =gis_phy%oz          &
                                  ,distgridToArrayMap=(/1/)             &
                                  ,indexflag=ESMF_INDEX_GLOBAL          &
                                  ,name       ='oz'                     &
                                  ,rc         =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
        ! attach array to state
        MESSAGE_CHECK="Add oz array to import state"
        CALL ESMF_StateAdd(state=IMP_STATE                              &
                          ,array=TMP_ARRAY                              &
                          ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ENDIF
      IF (gis_phy%esmf_sta_list%oz_export == 1) THEN
        MESSAGE_CHECK="Create oz array for export state"
        ! create the ESMF_Array
        TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                  ,farray     =gis_phy%oz          &
                                  ,distgridToArrayMap=(/1/)             &
                                  ,indexflag=ESMF_INDEX_GLOBAL          &
                                  ,name       ='oz'                     &
                                  ,rc         =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
        ! attach array to state
        MESSAGE_CHECK="Add oz array to export state"
        CALL ESMF_StateAdd(state=EXP_STATE                              &
                          ,array=TMP_ARRAY                              &
                          ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ENDIF

      IF (gis_phy%esmf_sta_list%cld_import == 1) THEN
        MESSAGE_CHECK="Create cld array for import state"
        ! create the ESMF_Array
        TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                  ,farray     =gis_phy%cld         &
                                  ,distgridToArrayMap=(/1/)             &
                                  ,indexflag=ESMF_INDEX_GLOBAL          &
                                  ,name       ='cld'                    &
                                  ,rc         =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
        ! attach array to state
        MESSAGE_CHECK="Add cld array to import state"
        CALL ESMF_StateAdd(state=IMP_STATE                              &
                          ,array=TMP_ARRAY                              &
                          ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ENDIF
      IF (gis_phy%esmf_sta_list%cld_export == 1) THEN
        MESSAGE_CHECK="Create cld array for export state"
        ! create the ESMF_Array
        TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                  ,farray     =gis_phy%cld         &
                                  ,distgridToArrayMap=(/1/)             &
                                  ,indexflag=ESMF_INDEX_GLOBAL          &
                                  ,name       ='cld'                    &
                                  ,rc         =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
        ! attach array to state
        MESSAGE_CHECK="Add cld array to export state"
        CALL ESMF_StateAdd(state=EXP_STATE                              &
                          ,array=TMP_ARRAY                              &
                          ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ENDIF

      IF (gis_phy%esmf_sta_list%p_import == 1) THEN
        MESSAGE_CHECK="Create p array for import state"
        ! create the ESMF_Array
        TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                  ,farray     =gis_phy%p           &
                                  ,distgridToArrayMap=(/1/)             &
                                  ,indexflag=ESMF_INDEX_GLOBAL          &
                                  ,name       ='p'                      &
                                  ,rc         =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
        ! attach array to state
        MESSAGE_CHECK="Add p array to import state"
        CALL ESMF_StateAdd(state=IMP_STATE                              &
                          ,array=TMP_ARRAY                              &
                          ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ENDIF
      IF (gis_phy%esmf_sta_list%p_export == 1) THEN
        MESSAGE_CHECK="Create p array for export state"
        ! create the ESMF_Array
        TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                  ,farray     =gis_phy%p           &
                                  ,distgridToArrayMap=(/1/)             &
                                  ,indexflag=ESMF_INDEX_GLOBAL          &
                                  ,name       ='p'                      &
                                  ,rc         =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
        ! attach array to state
        MESSAGE_CHECK="Add p array to export state"
        CALL ESMF_StateAdd(state=EXP_STATE                              &
                          ,array=TMP_ARRAY                              &
                          ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ENDIF

      IF (gis_phy%esmf_sta_list%dp_import == 1) THEN
        MESSAGE_CHECK="Create dp array for import state"
        ! create the ESMF_Array
        TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                  ,farray     =gis_phy%dp          &
                                  ,distgridToArrayMap=(/1/)             &
                                  ,indexflag=ESMF_INDEX_GLOBAL          &
                                  ,name       ='dp'                     &
                                  ,rc         =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
        ! attach array to state
        MESSAGE_CHECK="Add dp array to import state"
        CALL ESMF_StateAdd(state=IMP_STATE                              &
                          ,array=TMP_ARRAY                              &
                          ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ENDIF
      IF (gis_phy%esmf_sta_list%dp_export == 1) THEN
        MESSAGE_CHECK="Create dp array for export state"
        ! create the ESMF_Array
        TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                  ,farray     =gis_phy%dp          &
                                  ,distgridToArrayMap=(/1/)             &
                                  ,indexflag=ESMF_INDEX_GLOBAL          &
                                  ,name       ='dp'                     &
                                  ,rc         =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
        ! attach array to state
        MESSAGE_CHECK="Add dp array to export state"
        CALL ESMF_StateAdd(state=EXP_STATE                              &
                          ,array=TMP_ARRAY                              &
                          ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ENDIF

      IF (gis_phy%esmf_sta_list%dpdt_import == 1) THEN
        MESSAGE_CHECK="Create dpdt array for import state"
        ! create the ESMF_Array
        TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                  ,farray     =gis_phy%dpdt        &
                                  ,distgridToArrayMap=(/1/)             &
                                  ,indexflag=ESMF_INDEX_GLOBAL          &
                                  ,name       ='dpdt'                   &
                                  ,rc         =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
        ! attach array to state
        MESSAGE_CHECK="Add dpdt array to import state"
        CALL ESMF_StateAdd(state=IMP_STATE                              &
                          ,array=TMP_ARRAY                              &
                          ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ENDIF
      IF (gis_phy%esmf_sta_list%dpdt_export == 1) THEN
        MESSAGE_CHECK="Create dpdt array for export state"
        ! create the ESMF_Array
        TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                  ,farray     =gis_phy%dpdt        &
                                  ,distgridToArrayMap=(/1/)             &
                                  ,indexflag=ESMF_INDEX_GLOBAL          &
                                  ,name       ='dpdt'                   &
                                  ,rc         =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
        ! attach array to state
        MESSAGE_CHECK="Add dpdt array to export state"
        CALL ESMF_StateAdd(state=EXP_STATE                              &
                          ,array=TMP_ARRAY                              &
                          ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ENDIF

!TBH:  New arrays needed by FIM DYN component.  
!TODO: Need to reach agreement with NCEP about exporting these arrays from 
!TODO: GFS PHY component for use by FIM diagnostics.  
!     IF (gis_phy%esmf_sta_list%geshem_import == 1) THEN
      IF (.FALSE.) THEN
        MESSAGE_CHECK="Create geshem array for import state"
        ! create the ESMF_Array
        TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                  ,farray     =gis_phy%flx_fld%GESHEM   &
                                  ,distgridToArrayMap=(/1/)             &
                                  ,indexflag=ESMF_INDEX_GLOBAL          &
                                  ,name       ='geshem'                 &
                                  ,rc         =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
        ! attach array to state
        MESSAGE_CHECK="Add geshem array to import state"
        CALL ESMF_StateAdd(state=IMP_STATE                              &
                          ,array=TMP_ARRAY                              &
                          ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ENDIF
!     IF (gis_phy%esmf_sta_list%geshem_export == 1) THEN
      IF (.TRUE.) THEN
        MESSAGE_CHECK="Create geshem array for export state"
        ! create the ESMF_Array
        TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                  ,farray     =gis_phy%flx_fld%GESHEM   &
                                  ,distgridToArrayMap=(/1/)             &
                                  ,indexflag=ESMF_INDEX_GLOBAL          &
                                  ,name       ='geshem'                 &
                                  ,rc         =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
        ! attach array to state
        MESSAGE_CHECK="Add geshem array to export state"
        CALL ESMF_StateAdd(state=EXP_STATE                              &
                          ,array=TMP_ARRAY                              &
                          ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ENDIF

!     IF (gis_phy%esmf_sta_list%rainc_import == 1) THEN
      IF (.FALSE.) THEN
        MESSAGE_CHECK="Create rainc array for import state"
        ! create the ESMF_Array
        TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                  ,farray     =gis_phy%flx_fld%RAINC    &
                                  ,distgridToArrayMap=(/1/)             &
                                  ,indexflag=ESMF_INDEX_GLOBAL          &
                                  ,name       ='rainc'                  &
                                  ,rc         =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
        ! attach array to state
        MESSAGE_CHECK="Add rainc array to import state"
        CALL ESMF_StateAdd(state=IMP_STATE                              &
                          ,array=TMP_ARRAY                              &
                          ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ENDIF
!     IF (gis_phy%esmf_sta_list%rainc_export == 1) THEN
      IF (.TRUE.) THEN
        MESSAGE_CHECK="Create rainc array for export state"
        ! create the ESMF_Array
        TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                  ,farray     =gis_phy%flx_fld%RAINC    &
                                  ,distgridToArrayMap=(/1/)             &
                                  ,indexflag=ESMF_INDEX_GLOBAL          &
                                  ,name       ='rainc'                  &
                                  ,rc         =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
        ! attach array to state
        MESSAGE_CHECK="Add rainc array to export state"
        CALL ESMF_StateAdd(state=EXP_STATE                              &
                          ,array=TMP_ARRAY                              &
                          ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ENDIF

!     IF (gis_phy%esmf_sta_list%tsea_import == 1) THEN
      IF (.FALSE.) THEN
        MESSAGE_CHECK="Create tsea array for import state"
        ! create the ESMF_Array
        TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                  ,farray     =gis_phy%sfc_fld%TSEA     &
                                  ,distgridToArrayMap=(/1/)             &
                                  ,indexflag=ESMF_INDEX_GLOBAL          &
                                  ,name       ='tsea'                   &
                                  ,rc         =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
        ! attach array to state
        MESSAGE_CHECK="Add tsea array to import state"
        CALL ESMF_StateAdd(state=IMP_STATE                              &
                          ,array=TMP_ARRAY                              &
                          ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ENDIF
!     IF (gis_phy%esmf_sta_list%tsea_export == 1) THEN
      IF (.TRUE.) THEN
        MESSAGE_CHECK="Create tsea array for export state"
        ! create the ESMF_Array
        TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                  ,farray     =gis_phy%sfc_fld%TSEA     &
                                  ,distgridToArrayMap=(/1/)             &
                                  ,indexflag=ESMF_INDEX_GLOBAL          &
                                  ,name       ='tsea'                   &
                                  ,rc         =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
        ! attach array to state
        MESSAGE_CHECK="Add tsea array to export state"
        CALL ESMF_StateAdd(state=EXP_STATE                              &
                          ,array=TMP_ARRAY                              &
                          ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ENDIF

!     IF (gis_phy%esmf_sta_list%uustar_import == 1) THEN
      IF (.FALSE.) THEN
        MESSAGE_CHECK="Create uustar array for import state"
        ! create the ESMF_Array
        TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                  ,farray     =gis_phy%sfc_fld%UUSTAR   &
                                  ,distgridToArrayMap=(/1/)             &
                                  ,indexflag=ESMF_INDEX_GLOBAL          &
                                  ,name       ='uustar'                 &
                                  ,rc         =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
        ! attach array to state
        MESSAGE_CHECK="Add uustar array to import state"
        CALL ESMF_StateAdd(state=IMP_STATE                              &
                          ,array=TMP_ARRAY                              &
                          ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ENDIF
!     IF (gis_phy%esmf_sta_list%uustar_export == 1) THEN
      IF (.TRUE.) THEN
        MESSAGE_CHECK="Create uustar array for export state"
        ! create the ESMF_Array
        TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                  ,farray     =gis_phy%sfc_fld%UUSTAR   &
                                  ,distgridToArrayMap=(/1/)             &
                                  ,indexflag=ESMF_INDEX_GLOBAL          &
                                  ,name       ='uustar'                 &
                                  ,rc         =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
        ! attach array to state
        MESSAGE_CHECK="Add uustar array to export state"
        CALL ESMF_StateAdd(state=EXP_STATE                              &
                          ,array=TMP_ARRAY                              &
                          ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ENDIF

!     IF (gis_phy%esmf_sta_list%hflx_import == 1) THEN
      IF (.FALSE.) THEN
        MESSAGE_CHECK="Create hflx array for import state"
        ! create the ESMF_Array
        TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                  ,farray     =gis_phy%flx_fld%HFLX     &
                                  ,distgridToArrayMap=(/1/)             &
                                  ,indexflag=ESMF_INDEX_GLOBAL          &
                                  ,name       ='hflx'                   &
                                  ,rc         =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
        ! attach array to state
        MESSAGE_CHECK="Add hflx array to import state"
        CALL ESMF_StateAdd(state=IMP_STATE                              &
                          ,array=TMP_ARRAY                              &
                          ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ENDIF
!     IF (gis_phy%esmf_sta_list%hflx_export == 1) THEN
      IF (.TRUE.) THEN
        MESSAGE_CHECK="Create hflx array for export state"
        ! create the ESMF_Array
        TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                  ,farray     =gis_phy%flx_fld%HFLX     &
                                  ,distgridToArrayMap=(/1/)             &
                                  ,indexflag=ESMF_INDEX_GLOBAL          &
                                  ,name       ='hflx'                   &
                                  ,rc         =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
        ! attach array to state
        MESSAGE_CHECK="Add hflx array to export state"
        CALL ESMF_StateAdd(state=EXP_STATE                              &
                          ,array=TMP_ARRAY                              &
                          ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ENDIF

!     IF (gis_phy%esmf_sta_list%evap_import == 1) THEN
      IF (.FALSE.) THEN
        MESSAGE_CHECK="Create evap array for import state"
        ! create the ESMF_Array
        TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                  ,farray     =gis_phy%flx_fld%EVAP     &
                                  ,distgridToArrayMap=(/1/)             &
                                  ,indexflag=ESMF_INDEX_GLOBAL          &
                                  ,name       ='evap'                   &
                                  ,rc         =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
        ! attach array to state
        MESSAGE_CHECK="Add evap array to import state"
        CALL ESMF_StateAdd(state=IMP_STATE                              &
                          ,array=TMP_ARRAY                              &
                          ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ENDIF
!     IF (gis_phy%esmf_sta_list%evap_export == 1) THEN
      IF (.TRUE.) THEN
        MESSAGE_CHECK="Create evap array for export state"
        ! create the ESMF_Array
        TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                  ,farray     =gis_phy%flx_fld%EVAP     &
                                  ,distgridToArrayMap=(/1/)             &
                                  ,indexflag=ESMF_INDEX_GLOBAL          &
                                  ,name       ='evap'                   &
                                  ,rc         =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
        ! attach array to state
        MESSAGE_CHECK="Add evap array to export state"
        CALL ESMF_StateAdd(state=EXP_STATE                              &
                          ,array=TMP_ARRAY                              &
                          ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ENDIF

!     IF (gis_phy%esmf_sta_list%sheleg_import == 1) THEN
      IF (.FALSE.) THEN
        MESSAGE_CHECK="Create sheleg array for import state"
        ! create the ESMF_Array
        TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                  ,farray     =gis_phy%sfc_fld%SHELEG   &
                                  ,distgridToArrayMap=(/1/)             &
                                  ,indexflag=ESMF_INDEX_GLOBAL          &
                                  ,name       ='sheleg'                 &
                                  ,rc         =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
        ! attach array to state
        MESSAGE_CHECK="Add sheleg array to import state"
        CALL ESMF_StateAdd(state=IMP_STATE                              &
                          ,array=TMP_ARRAY                              &
                          ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ENDIF
!     IF (gis_phy%esmf_sta_list%sheleg_export == 1) THEN
      IF (.TRUE.) THEN
        MESSAGE_CHECK="Create sheleg array for export state"
        ! create the ESMF_Array
        TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                  ,farray     =gis_phy%sfc_fld%SHELEG   &
                                  ,distgridToArrayMap=(/1/)             &
                                  ,indexflag=ESMF_INDEX_GLOBAL          &
                                  ,name       ='sheleg'                 &
                                  ,rc         =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
        ! attach array to state
        MESSAGE_CHECK="Add sheleg array to export state"
        CALL ESMF_StateAdd(state=EXP_STATE                              &
                          ,array=TMP_ARRAY                              &
                          ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ENDIF

!     IF (gis_phy%esmf_sta_list%canopy_import == 1) THEN
      IF (.FALSE.) THEN
        MESSAGE_CHECK="Create canopy array for import state"
        ! create the ESMF_Array
        TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                  ,farray     =gis_phy%sfc_fld%CANOPY   &
                                  ,distgridToArrayMap=(/1/)             &
                                  ,indexflag=ESMF_INDEX_GLOBAL          &
                                  ,name       ='canopy'                 &
                                  ,rc         =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
        ! attach array to state
        MESSAGE_CHECK="Add canopy array to import state"
        CALL ESMF_StateAdd(state=IMP_STATE                              &
                          ,array=TMP_ARRAY                              &
                          ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ENDIF
!     IF (gis_phy%esmf_sta_list%canopy_export == 1) THEN
      IF (.TRUE.) THEN
        MESSAGE_CHECK="Create canopy array for export state"
        ! create the ESMF_Array
        TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                  ,farray     =gis_phy%sfc_fld%CANOPY   &
                                  ,distgridToArrayMap=(/1/)             &
                                  ,indexflag=ESMF_INDEX_GLOBAL          &
                                  ,name       ='canopy'                 &
                                  ,rc         =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
        ! attach array to state
        MESSAGE_CHECK="Add canopy array to export state"
        CALL ESMF_StateAdd(state=EXP_STATE                              &
                          ,array=TMP_ARRAY                              &
                          ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ENDIF

!     IF (gis_phy%esmf_sta_list%hice_import == 1) THEN
      IF (.FALSE.) THEN
        MESSAGE_CHECK="Create hice array for import state"
        ! create the ESMF_Array
        TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                  ,farray     =gis_phy%sfc_fld%HICE     &
                                  ,distgridToArrayMap=(/1/)             &
                                  ,indexflag=ESMF_INDEX_GLOBAL          &
                                  ,name       ='hice'                   &
                                  ,rc         =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
        ! attach array to state
        MESSAGE_CHECK="Add hice array to import state"
        CALL ESMF_StateAdd(state=IMP_STATE                              &
                          ,array=TMP_ARRAY                              &
                          ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ENDIF
!     IF (gis_phy%esmf_sta_list%hice_export == 1) THEN
      IF (.TRUE.) THEN
        MESSAGE_CHECK="Create hice array for export state"
        ! create the ESMF_Array
        TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                  ,farray     =gis_phy%sfc_fld%HICE     &
                                  ,distgridToArrayMap=(/1/)             &
                                  ,indexflag=ESMF_INDEX_GLOBAL          &
                                  ,name       ='hice'                   &
                                  ,rc         =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
        ! attach array to state
        MESSAGE_CHECK="Add hice array to export state"
        CALL ESMF_StateAdd(state=EXP_STATE                              &
                          ,array=TMP_ARRAY                              &
                          ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ENDIF

!     IF (gis_phy%esmf_sta_list%fice_import == 1) THEN
      IF (.FALSE.) THEN
        MESSAGE_CHECK="Create fice array for import state"
        ! create the ESMF_Array
        TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                  ,farray     =gis_phy%sfc_fld%FICE     &
                                  ,distgridToArrayMap=(/1/)             &
                                  ,indexflag=ESMF_INDEX_GLOBAL          &
                                  ,name       ='fice'                   &
                                  ,rc         =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
        ! attach array to state
        MESSAGE_CHECK="Add fice array to import state"
        CALL ESMF_StateAdd(state=IMP_STATE                              &
                          ,array=TMP_ARRAY                              &
                          ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ENDIF
!     IF (gis_phy%esmf_sta_list%fice_export == 1) THEN
      IF (.TRUE.) THEN
        MESSAGE_CHECK="Create fice array for export state"
        ! create the ESMF_Array
        TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                  ,farray     =gis_phy%sfc_fld%FICE     &
                                  ,distgridToArrayMap=(/1/)             &
                                  ,indexflag=ESMF_INDEX_GLOBAL          &
                                  ,name       ='fice'                   &
                                  ,rc         =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
        ! attach array to state
        MESSAGE_CHECK="Add fice array to export state"
        CALL ESMF_StateAdd(state=EXP_STATE                              &
                          ,array=TMP_ARRAY                              &
                          ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ENDIF

!     IF (gis_phy%esmf_sta_list%stc_import == 1) THEN
      IF (.FALSE.) THEN
        MESSAGE_CHECK="Create stc array for import state"
        ! create the ESMF_Array
        TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                  ,farray     =gis_phy%sfc_fld%STC      &
                                  ,distgridToArrayMap=(/2/)             &
                                  ,indexflag=ESMF_INDEX_GLOBAL          &
                                  ,name       ='stc'                    &
                                  ,rc         =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
        ! attach array to state
        MESSAGE_CHECK="Add stc array to import state"
        CALL ESMF_StateAdd(state=IMP_STATE                              &
                          ,array=TMP_ARRAY                              &
                          ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ENDIF
!     IF (gis_phy%esmf_sta_list%stc_export == 1) THEN
      IF (.TRUE.) THEN
        MESSAGE_CHECK="Create stc array for export state"
        ! create the ESMF_Array
        TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                  ,farray     =gis_phy%sfc_fld%STC      &
                                  ,distgridToArrayMap=(/2/)             &
                                  ,indexflag=ESMF_INDEX_GLOBAL          &
                                  ,name       ='stc'                    &
                                  ,rc         =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
        ! attach array to state
        MESSAGE_CHECK="Add stc array to export state"
        CALL ESMF_StateAdd(state=EXP_STATE                              &
                          ,array=TMP_ARRAY                              &
                          ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ENDIF

!     IF (gis_phy%esmf_sta_list%smc_import == 1) THEN
      IF (.FALSE.) THEN
        MESSAGE_CHECK="Create smc array for import state"
        ! create the ESMF_Array
        TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                  ,farray     =gis_phy%sfc_fld%SMC      &
                                  ,distgridToArrayMap=(/2/)             &
                                  ,indexflag=ESMF_INDEX_GLOBAL          &
                                  ,name       ='smc'                    &
                                  ,rc         =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
        ! attach array to state
        MESSAGE_CHECK="Add smc array to import state"
        CALL ESMF_StateAdd(state=IMP_STATE                              &
                          ,array=TMP_ARRAY                              &
                          ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ENDIF
!     IF (gis_phy%esmf_sta_list%smc_export == 1) THEN
      IF (.TRUE.) THEN
        MESSAGE_CHECK="Create smc array for export state"
        ! create the ESMF_Array
        TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                  ,farray     =gis_phy%sfc_fld%SMC      &
                                  ,distgridToArrayMap=(/2/)             &
                                  ,indexflag=ESMF_INDEX_GLOBAL          &
                                  ,name       ='smc'                    &
                                  ,rc         =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
        ! attach array to state
        MESSAGE_CHECK="Add smc array to export state"
        CALL ESMF_StateAdd(state=EXP_STATE                              &
                          ,array=TMP_ARRAY                              &
                          ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ENDIF

!     IF (gis_phy%esmf_sta_list%sfcdsw_import == 1) THEN
      IF (.FALSE.) THEN
        MESSAGE_CHECK="Create sfcdsw array for import state"
        ! create the ESMF_Array
        TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                  ,farray     =gis_phy%flx_fld%SFCDSW   &
                                  ,distgridToArrayMap=(/1/)             &
                                  ,indexflag=ESMF_INDEX_GLOBAL          &
                                  ,name       ='sfcdsw'                 &
                                  ,rc         =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
        ! attach array to state
        MESSAGE_CHECK="Add sfcdsw array to import state"
        CALL ESMF_StateAdd(state=IMP_STATE                              &
                          ,array=TMP_ARRAY                              &
                          ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ENDIF
!     IF (gis_phy%esmf_sta_list%sfcdsw_export == 1) THEN
      IF (.TRUE.) THEN
        MESSAGE_CHECK="Create sfcdsw array for export state"
        ! create the ESMF_Array
        TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                  ,farray     =gis_phy%flx_fld%SFCDSW   &
                                  ,distgridToArrayMap=(/1/)             &
                                  ,indexflag=ESMF_INDEX_GLOBAL          &
                                  ,name       ='sfcdsw'                 &
                                  ,rc         =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
        ! attach array to state
        MESSAGE_CHECK="Add sfcdsw array to export state"
        CALL ESMF_StateAdd(state=EXP_STATE                              &
                          ,array=TMP_ARRAY                              &
                          ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ENDIF

!     IF (gis_phy%esmf_sta_list%sfcdlw_import == 1) THEN
      IF (.FALSE.) THEN
        MESSAGE_CHECK="Create sfcdlw array for import state"
        ! create the ESMF_Array
        TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                  ,farray     =gis_phy%flx_fld%SFCDLW   &
                                  ,distgridToArrayMap=(/1/)             &
                                  ,indexflag=ESMF_INDEX_GLOBAL          &
                                  ,name       ='sfcdlw'                 &
                                  ,rc         =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
        ! attach array to state
        MESSAGE_CHECK="Add sfcdlw array to import state"
        CALL ESMF_StateAdd(state=IMP_STATE                              &
                          ,array=TMP_ARRAY                              &
                          ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ENDIF
!     IF (gis_phy%esmf_sta_list%sfcdlw_export == 1) THEN
      IF (.TRUE.) THEN
        MESSAGE_CHECK="Create sfcdlw array for export state"
        ! create the ESMF_Array
        TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                  ,farray     =gis_phy%flx_fld%SFCDLW   &
                                  ,distgridToArrayMap=(/1/)             &
                                  ,indexflag=ESMF_INDEX_GLOBAL          &
                                  ,name       ='sfcdlw'                 &
                                  ,rc         =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
        ! attach array to state
        MESSAGE_CHECK="Add sfcdlw array to export state"
        CALL ESMF_StateAdd(state=EXP_STATE                              &
                          ,array=TMP_ARRAY                              &
                          ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ENDIF

!     IF (gis_phy%esmf_sta_list%t2m_import == 1) THEN
      IF (.FALSE.) THEN
        MESSAGE_CHECK="Create t2m array for import state"
        ! create the ESMF_Array
        TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                  ,farray     =gis_phy%sfc_fld%T2M      &
                                  ,distgridToArrayMap=(/1/)             &
                                  ,indexflag=ESMF_INDEX_GLOBAL          &
                                  ,name       ='t2m'                    &
                                  ,rc         =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
        ! attach array to state
        MESSAGE_CHECK="Add t2m array to import state"
        CALL ESMF_StateAdd(state=IMP_STATE                              &
                          ,array=TMP_ARRAY                              &
                          ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ENDIF
!     IF (gis_phy%esmf_sta_list%t2m_export == 1) THEN
      IF (.TRUE.) THEN
        MESSAGE_CHECK="Create t2m array for export state"
        ! create the ESMF_Array
        TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                  ,farray     =gis_phy%sfc_fld%T2M      &
                                  ,distgridToArrayMap=(/1/)             &
                                  ,indexflag=ESMF_INDEX_GLOBAL          &
                                  ,name       ='t2m'                    &
                                  ,rc         =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
        ! attach array to state
        MESSAGE_CHECK="Add t2m array to export state"
        CALL ESMF_StateAdd(state=EXP_STATE                              &
                          ,array=TMP_ARRAY                              &
                          ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ENDIF

!     IF (gis_phy%esmf_sta_list%q2m_import == 1) THEN
      IF (.FALSE.) THEN
        MESSAGE_CHECK="Create q2m array for import state"
        ! create the ESMF_Array
        TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                  ,farray     =gis_phy%sfc_fld%Q2M      &
                                  ,distgridToArrayMap=(/1/)             &
                                  ,indexflag=ESMF_INDEX_GLOBAL          &
                                  ,name       ='q2m'                    &
                                  ,rc         =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
        ! attach array to state
        MESSAGE_CHECK="Add q2m array to import state"
        CALL ESMF_StateAdd(state=IMP_STATE                              &
                          ,array=TMP_ARRAY                              &
                          ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ENDIF
!     IF (gis_phy%esmf_sta_list%q2m_export == 1) THEN
      IF (.TRUE.) THEN
        MESSAGE_CHECK="Create q2m array for export state"
        ! create the ESMF_Array
        TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                  ,farray     =gis_phy%sfc_fld%Q2M      &
                                  ,distgridToArrayMap=(/1/)             &
                                  ,indexflag=ESMF_INDEX_GLOBAL          &
                                  ,name       ='q2m'                    &
                                  ,rc         =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
        ! attach array to state
        MESSAGE_CHECK="Add q2m array to export state"
        CALL ESMF_StateAdd(state=EXP_STATE                              &
                          ,array=TMP_ARRAY                              &
                          ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ENDIF

!     IF (gis_phy%esmf_sta_list%slmsk_import == 1) THEN
      IF (.FALSE.) THEN
        MESSAGE_CHECK="Create slmsk array for import state"
        ! create the ESMF_Array
        TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                  ,farray     =gis_phy%sfc_fld%SLMSK    &
                                  ,distgridToArrayMap=(/1/)             &
                                  ,indexflag=ESMF_INDEX_GLOBAL          &
                                  ,name       ='slmsk'                  &
                                  ,rc         =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
        ! attach array to state
        MESSAGE_CHECK="Add slmsk array to import state"
        CALL ESMF_StateAdd(state=IMP_STATE                              &
                          ,array=TMP_ARRAY                              &
                          ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ENDIF
!     IF (gis_phy%esmf_sta_list%slmsk_export == 1) THEN
      IF (.TRUE.) THEN
        MESSAGE_CHECK="Create slmsk array for export state"
        ! create the ESMF_Array
        TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                  ,farray     =gis_phy%sfc_fld%SLMSK    &
                                  ,distgridToArrayMap=(/1/)             &
                                  ,indexflag=ESMF_INDEX_GLOBAL          &
                                  ,name       ='slmsk'                  &
                                  ,rc         =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
        ! attach array to state
        MESSAGE_CHECK="Add slmsk array to export state"
        CALL ESMF_StateAdd(state=EXP_STATE                              &
                          ,array=TMP_ARRAY                              &
                          ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ENDIF

!     IF (gis_phy%esmf_sta_list%hprime_import == 1) THEN
      IF (.FALSE.) THEN
        MESSAGE_CHECK="Create hprime array for import state"
        ! create the ESMF_Array
        TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                  ,farray     =gis_phy%hprime           &
                                  ,distgridToArrayMap=(/2/)             &
                                  ,indexflag=ESMF_INDEX_GLOBAL          &
                                  ,name       ='hprime'                 &
                                  ,rc         =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
        ! attach array to state
        MESSAGE_CHECK="Add hprime array to import state"
        CALL ESMF_StateAdd(state=IMP_STATE                              &
                          ,array=TMP_ARRAY                              &
                          ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ENDIF
!     IF (gis_phy%esmf_sta_list%hprime_export == 1) THEN
      IF (.TRUE.) THEN
        MESSAGE_CHECK="Create hprime array for export state"
        ! create the ESMF_Array
        TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                  ,farray     =gis_phy%hprime           &
                                  ,distgridToArrayMap=(/2/)             &
                                  ,indexflag=ESMF_INDEX_GLOBAL          &
                                  ,name       ='hprime'                 &
                                  ,rc         =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
        ! attach array to state
        MESSAGE_CHECK="Add hprime array to export state"
        CALL ESMF_StateAdd(state=EXP_STATE                              &
                          ,array=TMP_ARRAY                              &
                          ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ENDIF

!     IF (gis_phy%esmf_sta_list%fluxr_import == 1) THEN
      IF (.FALSE.) THEN
        MESSAGE_CHECK="Create fluxr array for import state"
        ! create the ESMF_Array
        TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                  ,farray     =gis_phy%fluxr            &
                                  ,distgridToArrayMap=(/2/)             &
                                  ,indexflag=ESMF_INDEX_GLOBAL          &
                                  ,name       ='fluxr'                  &
                                  ,rc         =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
        ! attach array to state
        MESSAGE_CHECK="Add fluxr array to import state"
        CALL ESMF_StateAdd(state=IMP_STATE                              &
                          ,array=TMP_ARRAY                              &
                          ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ENDIF
!     IF (gis_phy%esmf_sta_list%fluxr_export == 1) THEN
      IF (.TRUE.) THEN
        MESSAGE_CHECK="Create fluxr array for export state"
        ! create the ESMF_Array
        TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                  ,farray     =gis_phy%fluxr            &
                                  ,distgridToArrayMap=(/2/)             &
                                  ,indexflag=ESMF_INDEX_GLOBAL          &
                                  ,name       ='fluxr'                  &
                                  ,rc         =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
        ! attach array to state
        MESSAGE_CHECK="Add fluxr array to export state"
        CALL ESMF_StateAdd(state=EXP_STATE                              &
                          ,array=TMP_ARRAY                              &
                          ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ENDIF

!TBH:  validate states
      MESSAGE_CHECK="PHY_INITIALIZE:  Validate import state"
      call ESMF_StateValidate(state=IMP_STATE,rc=rc)
      IF(RC==ESMF_SUCCESS)THEN
!JR        WRITE(0,*)'PHY INITIALIZE import state valid'
      ENDIF
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      MESSAGE_CHECK="PHY_INITIALIZE:  Validate export state"
      call ESMF_StateValidate(state=EXP_STATE,rc=rc)
      IF(RC==ESMF_SUCCESS)THEN
!JR        WRITE(0,*)'PHY INITIALIZE export state valid'
      ENDIF
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)

!
!-----------------------------------------------------------------------
!
      RC_INIT = RC
      IF(RC_INIT==ESMF_SUCCESS)THEN
!       WRITE(0,*)'PHY INITIALIZE STEP SUCCEEDED'
      ELSE
        WRITE(0,*)'PHY INITIALIZE STEP FAILED RC_INIT=',RC_INIT
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE PHY_INITIALIZE
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE PHY_RUN(GRID_COMP                                      &
                        ,IMP_STATE,EXP_STATE                            &
                        ,CLOCK_ATM                                      &
                        ,RC_RUN)
!
!-----------------------------------------------------------------------
!***  THE INTEGRATION OF THE MODEL PHYSICS IS DONE
!***  THROUGH THIS ROUTINE.
!-----------------------------------------------------------------------
!
      USE module_fim_phy_run ,only: PHY_RUN_FIM => phy_run
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!***  ARGUMENT VARIABLES
!-----------------------------------------------------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: GRID_COMP  !<-- The Physics gridded component
      TYPE(ESMF_State)   ,INTENT(INOUT) :: IMP_STATE  !<-- The Physics import state
      TYPE(ESMF_State)   ,INTENT(INOUT) :: EXP_STATE  !<-- The Physics export state
      TYPE(ESMF_Clock)   ,INTENT(IN)    :: CLOCK_ATM  !<-- The ATM's ESMF Clock
!
      INTEGER            ,INTENT(OUT)   :: RC_RUN
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      INTEGER :: NTIMESTEP,RC
!
      INTEGER(KIND=ESMF_KIND_I8) :: NTIMESTEP_ESMF
!
      INTEGER :: its
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC_RUN=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  EXTRACT THE TIMESTEP COUNT FROM THE CLOCK.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Retrieve Timestep from ATM Clock in Physics Run"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ClockGet(clock       =CLOCK_ATM                         &  !<-- The ESMF clock
                        ,advanceCount=NTIMESTEP_ESMF                    &  !<-- The number of times the clock has been advanced
                        ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      NTIMESTEP=NTIMESTEP_ESMF
!
! NOTE:  Pointers in import and export states point to internal state as 
! NOTE:  set up in the init phase, consistent with future plans for NEMS.  
! NOTE:  So wrap%int_state is not needed here at present, nor are explicit 
! NOTE:  transfers between internal and import/export states.  
!TODO:  adjust as plans evolve

!-----------------------------------------------------------------------
!***  CALL THE INDIVIDUAL PHYSICAL PROCESSES
!-----------------------------------------------------------------------
!
      its = NTIMESTEP + 1
      CALL PHY_RUN_FIM (its)
!
      RC_RUN=RC
!
!-----------------------------------------------------------------------
!
      IF(RC_RUN==ESMF_SUCCESS)THEN
!       WRITE(0,*)'PHY RUN STEP SUCCEEDED'
      ELSE
        WRITE(0,*)'PHY RUN STEP FAILED RC_RUN=',RC_RUN
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE PHY_RUN
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE PHY_FINALIZE(GRID_COMP                                      &
                             ,IMP_STATE,EXP_STATE                            &
                             ,CLOCK_ATM                                      &
                             ,RC_FINAL)
!
!-----------------------------------------------------------------------
!***  FINALIZE THE PHYSICS COMPONENT.
!-----------------------------------------------------------------------
!
      USE module_fim_phy_finalize ,only: PHY_FINALIZE_FIM => phy_finalize
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!***  ARGUMENT VARIABLES
!-----------------------------------------------------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: GRID_COMP  !<-- The Physics gridded component
      TYPE(ESMF_State)   ,INTENT(INOUT) :: IMP_STATE  !<-- The Physics import state
      TYPE(ESMF_State)   ,INTENT(INOUT) :: EXP_STATE  !<-- The Physics export state
      TYPE(ESMF_Clock)   ,INTENT(IN)    :: CLOCK_ATM  !<-- The ATM's ESMF Clock
!
      INTEGER            ,INTENT(OUT)   :: RC_FINAL
!      
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC_FINAL=ESMF_SUCCESS
!
      CALL PHY_FINALIZE_FIM
!
!      WRITE(0,*)' Physics Completed Normally.'
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE PHY_FINALIZE
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      END MODULE MODULE_PHYSICS_GRID_COMP
!
!-----------------------------------------------------------------------
