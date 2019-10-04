!TODO:  DRY out all of this code.  Follow NMM "array object" approach 
!TODO:  developed by Tom Black and I (TBH - we called it the "ownership" 
!TODO:  proposal) and implemented within NEMS-NMMB by Dusan Jovic.  
!-----------------------------------------------------------------------
!
      MODULE MODULE_DYNAMICS_GRID_COMP
!
!-----------------------------------------------------------------------
!
!***  THIS MODULE HOLDS THE DYNAMICS REGISTER, INIT, RUN, AND FINALIZE 
!***  ROUTINES.  THEY ARE CALLED FROM THE FIM GRIDDED COMPONENT
!***  (FIM INITIALIZE CALLS DYNAMICS INITIALIZE, ETC.).  
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
#ifdef MANUALGPTL
#include <gptl.inc>
      integer :: ret
#endif

! TODO:  put this in "dyn_internal_state"
      TYPE(ESMF_Grid), SAVE :: GRID_FIM_DYN  !<-- The ESMF GRID for FIM "nip" dimension

!
!-----------------------------------------------------------------------
!
      PRIVATE
!
      PUBLIC :: DYN_REGISTER
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE DYN_REGISTER(DYN_GRID_COMP,RC_REG)
!
!-----------------------------------------------------------------------
!***  REGISTER THE DYNAMICS COMPONENT'S INITIALIZE, RUN, AND FINALIZE
!***  SUBROUTINE NAMES.
!-----------------------------------------------------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: DYN_GRID_COMP                  !<-- The Dynamics gridded component
!
      INTEGER,INTENT(OUT) :: RC_REG                                       !<-- Return code for Dyn register
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      INTEGER :: RC
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC_REG=ESMF_SUCCESS                                                 !<-- Initialize error signal variable
                                                                                                                                              
!-----------------------------------------------------------------------
!***  REGISTER THE DYNAMICS INITIALIZE SUBROUTINE.  SINCE IT IS JUST ONE
!***  SUBROUTINE, USE ESMF_SINGLEPHASE.  THE SECOND ARGUMENT IS
!***  A PRE-DEFINED SUBROUTINE TYPE, SUCH AS ESMF_SETINIT, ESMF_SETRUN,
!***  OR ESMF_SETFINAL.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set Entry Point for Dynamics Initialize"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifdef MANUALGPTL
      ret = gptlstart ('ESMF_GridCompSetEntryPoint:dyn_initialize')
#endif
      CALL ESMF_GridCompSetEntryPoint(DYN_GRID_COMP                     &  !<-- The gridded component
                                     ,ESMF_SETINIT                      &  !<-- Predefined subroutine type
                                     ,DYN_INITIALIZE                    &  !<-- User's subroutineName
                                     ,ESMF_SINGLEPHASE                  &  !<-- phase
                                     ,RC)
#ifdef MANUALGPTL
      ret = gptlstop ('ESMF_GridCompSetEntryPoint:dyn_initialize')
#endif
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  REGISTER THE DYNAMICS RUN SUBROUTINE.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set Entry Point for Dynamics Run"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifdef MANUALGPTL
      ret = gptlstart ('ESMF_GridCompSetEntryPoint:dyn_run')
#endif
      CALL ESMF_GridCompSetEntryPoint(DYN_GRID_COMP                     &  !<-- gridcomp
                                     ,ESMF_SETRUN                       &  !<-- subroutineType
                                     ,DYN_RUN                           &  !<-- user's subroutineName
                                     ,ESMF_SINGLEPHASE                  &  !<-- phase
                                     ,RC)
#ifdef MANUALGPTL
      ret = gptlstop ('ESMF_GridCompSetEntryPoint:dyn_run')
#endif
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  REGISTER THE DYNAMICS FINALIZE SUBROUTINE.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set Entry Point for Dynamics Finalize"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompSetEntryPoint(DYN_GRID_COMP                     &  !<-- gridcomp
                                     ,ESMF_SETFINAL                     &  !<-- subroutineType
                                     ,DYN_FINALIZE                      &  !<-- user's subroutineName
                                     ,ESMF_SINGLEPHASE                  &  !<-- phase
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
!       WRITE(0,*)" DYN_REGISTER SUCCEEDED"
      ELSE
        WRITE(0,*)" DYN_REGISTER FAILED"
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE DYN_REGISTER
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE DYN_INITIALIZE(DYN_GRID_COMP                           &
                               ,IMP_STATE,EXP_STATE                     &
                               ,CLOCK_FIM                               &
                               ,RC_INIT)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  CARRY OUT ALL NECESSARY SETUPS FOR THE MODEL DYNAMICS.
! This includes creating the FIM ESMF_Grid and attaching it to dyn_grid_comp.  
!-----------------------------------------------------------------------
!
      use module_control,only: nip,nvl
      USE module_fim_dyn_init ,only: DYN_INITIALIZE_FIM => dyn_init
      use module_variables,only: us3d,vs3d,pr3d,tr3d,ws3d
      use module_sfc_variables,only: rn2d,rc2d,ts2d,us2d,hf2d,qf2d,        &
                                     sheleg2d, canopy2d, hice2d, fice2d,   &
                                     st3d, sm3d, sw2d, lw2d, t2m2d, q2m2d, &
                                     slmsk2d, hprm2d, flxlwtoa2d
!
!-----------------------------------------------------------------------
!***  ARGUMENT VARIABLES.
!-----------------------------------------------------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: DYN_GRID_COMP                   !<-- The Dynamics gridded component
      TYPE(ESMF_State),   INTENT(INOUT) :: IMP_STATE                       !<-- The Dynamics Initialize step's import state
      TYPE(ESMF_State),   INTENT(INOUT) :: EXP_STATE                       !<-- The Dynamics Initialize step's export state
      TYPE(ESMF_Clock),   INTENT(IN)    :: CLOCK_FIM                       !<-- The FIM's ESMF Clock
!
      INTEGER,            INTENT(OUT)   :: RC_INIT
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      ! temporary field object
!      TYPE(ESMF_Field) :: TMP_FIELD
      TYPE(ESMF_DistGrid) :: DISTGRID
      TYPE(ESMF_Array) :: TMP_ARRAY
      type(esmf_vm),save :: vm_local   ! TODO:  is SAVE needed?  
      INTEGER :: NUM_PES_FCST,mype
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
#ifdef MANUALGPTL
      ret = gptlstart ('dyn_initialize')
#endif

      RC     =ESMF_SUCCESS
      RC_INIT=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  PRIMARY INITIALIZATION OF SCALARS/ARRAYS.
!***  Also sets up SMS decomposition.  
!-----------------------------------------------------------------------
!
#ifdef MANUALGPTL
      ret = gptlstart ('dyn_init')
#endif
      CALL DYN_INITIALIZE_FIM(.false.)
#ifdef MANUALGPTL
      ret = gptlstop ('dyn_init')
#endif
!
!-----------------------------------------------------------------------
!***  CREATE THE ESMF GRID.
!-----------------------------------------------------------------------
! Create ESMF_Grid and attach it to dyn_grid_comp.  
!TBH:  For FIM grid creation the ESMF_Grid is currently a dummy.  This 
!TBH:  routine creates a bogus grid and distributes it across the 
!TBH:  correct number of FIM compute tasks.  It then attaches it to 
!TBH:  dyn_grid_comp.  
! TODO:  Match decomposition to SMS.  
! TODO:  Load lat-lon into this grid for later extraction by PHY, 
! TODO:  may have to defer this to newer ESMF version!  
!-----------------------------------------------------------------------
!
#ifdef MANUALGPTL
      ret = gptlstart ('ESMF_gridcreate')
#endif
!TODO:  Match decomposition to SMS via new args passed out of 
!TODO:  DYN_INITIALIZE_FIM.  
!TODO:  Replace with a real icos constructor once we switch to a version 
!TODO:  of ESMF that supports it.  
!TODO:  Load lat-lon into this grid for later extraction by PHY, 
!TODO:  may have to defer this to newer ESMF version!  

      call esmf_gridcompget(gridcomp=dyn_grid_comp, vm=vm_local, rc=rc)
      CALL err_msg(RC,"DYN_INITIALIZE: get local vm",rc_init)
      call esmf_vmget(vm=vm_local, pecount=NUM_PES_FCST, localpet=mype,&
                      rc=rc)
      CALL err_msg(RC,"DYN_INITIALIZE: get NUM_PES_FCST",rc_init)

!TODO:  Add deBlockList to specify start and end indicies for each MPI task, 
!TODO:  need NUM_PES_FCST for deBlockList(0:NUM_PES_FCST-1) ...  
      ! Create "1D" ESMF_DistGrid per ESMF Reference Manual section 26.2.1
      DISTGRID=ESMF_DistGridCreate(minIndex=(/1/),maxIndex=(/nip/), &
                                   rc=rc_init)
      CALL err_msg(RC,"DYN_INITIALIZE: create DISTGRID",rc_init)

      ! Create "1D" ESMF_Grid from ESMF_DistGrid
      GRID_FIM_DYN=ESMF_GridCreate(name="FIM_GRID",distgrid=DISTGRID, &
                                   rc=rc_init)
      CALL err_msg(RC,"DYN_INITIALIZE: create GRID from DISTGRID", &
                   rc_init)
#ifdef MANUALGPTL
      ret = gptlstop ('ESMF_gridcreate')
#endif

      call esmf_gridvalidate(grid=GRID_FIM_DYN, rc=RC)
      CALL ERR_MSG(RC,"DYN_INITIALIZE:  Validate new GRID",RC_INIT)

!TODO:  Add guard to prevent memory leak if parent has already attached 
!TODO:  an ESMF_Grid to dyn_grid_comp.  
      ! attach ESMF_Grid to ESMF_GridComp
      call esmf_gridcompset(dyn_grid_comp, grid=grid_fim_dyn, &
                            rc=rc)
      CALL err_msg(rc, "attach FIM grid to DYN component", &
                   rc_init)
!
!-----------------------------------------------------------------------
!***  Attach FIM fields in the internal state 
!***  to the esmf import and export states.  
!TBH:  I use GFS naming conventions, *not* NMMB conventions.  As 
!TBH:  of NEMS r3038 they do indeed differ, by case at least!  
!TBH:  Creation of unique ESMF_Field objects for import and export 
!TBH:  states should require little additional memory since the pointers 
!TBH:  to Fortran arrays are shared.  This approach makes object deletion 
!TBH:  easier.  
!-------------------------------------------------------

!TODO:  Implement esmf_sta_list to allow config control of coupling once 
!TODO:  this settles down between NMMB and GFS in NEMS.  

!TODO:  Need to add gridToFieldMap to ESMF_FieldCreate() to address 
!TODO:  differences between 2D and 3D arrays.  At present this is 
!TODO:  irrelevant since we do not use ESMF to do any re-grid or 
!TODO:  re-dist operations.  This must be fixed before we use these 
!TODO:  ESMF features.  

!
! pr3d
!
      MESSAGE_CHECK="Create pr3d array for import state"
      ! create the ESMF_Field
!TBH:  Note that the following call to ESMF_FieldCreate() yields the 
!TBH:  stunningly informative error code 540 which maps to string 
!TBH:  "Not valid" in ESMC_ErrMsgs.C.  Backed off to ESMF_ArrayCreate().  
!TODO:  Switch back to ESMF_FieldCreate() since future NEMS will use 
!TODO:  ESMF_Fields.  
!      TMP_FIELD=ESMF_FieldCreate(grid       =grid_fim_dyn             &
!                                ,farray     =pr3d                     &
!                                ,distgridToArrayMap=(/2/)             &
!                                ,name       ='pr3d'                   &
!                                ,rc         =RC)
#ifdef MANUALGPTL
      ret = gptlstart ('esmf_arraycreate')
#endif
      TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                ,farray     =pr3d                     &
                                ,distgridToArrayMap=(/2/)             &
                                ,indexflag=ESMF_INDEX_GLOBAL          &
                                ,name       ='pr3d'                   &
                                ,rc         =RC)
#ifdef MANUALGPTL
      ret = gptlstop ('esmf_arraycreate')
#endif
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ! attach array to state
      MESSAGE_CHECK="Add pr3d array to import state"
#ifdef MANUALGPTL
      ret = gptlstart ('esmf_stateadd')
#endif
      CALL ESMF_StateAdd(state=IMP_STATE                              &
                        ,array=TMP_ARRAY                              &
                        ,rc   =RC)
#ifdef MANUALGPTL
      ret = gptlstop ('esmf_stateadd')
#endif
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      MESSAGE_CHECK="Create pr3d array for export state"
      ! create the ESMF_Field
      TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                ,farray     =pr3d                     &
                                ,distgridToArrayMap=(/2/)             &
                                ,indexflag=ESMF_INDEX_GLOBAL          &
                                ,name       ='pr3d'                   &
                                ,rc         =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ! attach array to state
      MESSAGE_CHECK="Add pr3d array to export state"
      CALL ESMF_StateAdd(state=EXP_STATE                              &
                        ,array=TMP_ARRAY                              &
                        ,rc   =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
!
! us3d
!
      MESSAGE_CHECK="Create us3d array for import state"
      ! create the ESMF_Field
      TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                ,farray     =us3d                     &
                                ,distgridToArrayMap=(/2/)             &
                                ,indexflag=ESMF_INDEX_GLOBAL          &
                                ,name       ='us3d'                   &
                                ,rc         =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ! attach array to state
      MESSAGE_CHECK="Add us3d array to import state"
      CALL ESMF_StateAdd(state=IMP_STATE                              &
                        ,array=TMP_ARRAY                              &
                        ,rc   =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      MESSAGE_CHECK="Create us3d array for export state"
      ! create the ESMF_Field
      TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                ,farray     =us3d                     &
                                ,distgridToArrayMap=(/2/)             &
                                ,indexflag=ESMF_INDEX_GLOBAL          &
                                ,name       ='us3d'                   &
                                ,rc         =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ! attach array to state
      MESSAGE_CHECK="Add us3d array to export state"
      CALL ESMF_StateAdd(state=EXP_STATE                              &
                        ,array=TMP_ARRAY                              &
                        ,rc   =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
!
! vs3d
!
      MESSAGE_CHECK="Create vs3d array for import state"
      ! create the ESMF_Field
      TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                ,farray     =vs3d                     &
                                ,distgridToArrayMap=(/2/)             &
                                ,indexflag=ESMF_INDEX_GLOBAL          &
                                ,name       ='vs3d'                   &
                                ,rc         =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ! attach array to state
      MESSAGE_CHECK="Add vs3d array to import state"
      CALL ESMF_StateAdd(state=IMP_STATE                              &
                        ,array=TMP_ARRAY                              &
                        ,rc   =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      MESSAGE_CHECK="Create vs3d array for export state"
      ! create the ESMF_Field
      TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                ,farray     =vs3d                     &
                                ,distgridToArrayMap=(/2/)             &
                                ,indexflag=ESMF_INDEX_GLOBAL          &
                                ,name       ='vs3d'                   &
                                ,rc         =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ! attach array to state
      MESSAGE_CHECK="Add vs3d array to export state"
      CALL ESMF_StateAdd(state=EXP_STATE                              &
                        ,array=TMP_ARRAY                              &
                        ,rc   =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
!
! ws3d
!
      MESSAGE_CHECK="Create ws3d array for import state"
      ! create the ESMF_Field
      TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                ,farray     =ws3d                     &
                                ,distgridToArrayMap=(/2/)             &
                                ,indexflag=ESMF_INDEX_GLOBAL          &
                                ,name       ='ws3d'                   &
                                ,rc         =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ! attach array to state
      MESSAGE_CHECK="Add ws3d array to import state"
      CALL ESMF_StateAdd(state=IMP_STATE                              &
                        ,array=TMP_ARRAY                              &
                        ,rc   =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      MESSAGE_CHECK="Create ws3d array for export state"
      ! create the ESMF_Field
      TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                ,farray     =ws3d                     &
                                ,distgridToArrayMap=(/2/)             &
                                ,indexflag=ESMF_INDEX_GLOBAL          &
                                ,name       ='ws3d'                   &
                                ,rc         =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ! attach array to state
      MESSAGE_CHECK="Add ws3d array to export state"
      CALL ESMF_StateAdd(state=EXP_STATE                              &
                        ,array=TMP_ARRAY                              &
                        ,rc   =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
!
! tr3d
!
      MESSAGE_CHECK="Create tr3d array for import state"
      ! create the ESMF_Field
      TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                ,farray     =tr3d                     &
                                ,distgridToArrayMap=(/2/)             &
                                ,indexflag=ESMF_INDEX_GLOBAL          &
                                ,name       ='tr3d'                   &
                                ,rc         =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ! attach array to state
      MESSAGE_CHECK="Add tr3d array to import state"
      CALL ESMF_StateAdd(state=IMP_STATE                              &
                        ,array=TMP_ARRAY                              &
                        ,rc   =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      MESSAGE_CHECK="Create tr3d array for export state"
      ! create the ESMF_Field
      TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                ,farray     =tr3d                     &
                                ,distgridToArrayMap=(/2/)             &
                                ,indexflag=ESMF_INDEX_GLOBAL          &
                                ,name       ='tr3d'                   &
                                ,rc         =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ! attach array to state
      MESSAGE_CHECK="Add tr3d array to export state"
      CALL ESMF_StateAdd(state=EXP_STATE                              &
                        ,array=TMP_ARRAY                              &
                        ,rc   =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
!
! rn2d
!
      MESSAGE_CHECK="Create rn2d array for import state"
      ! create the ESMF_Field
      TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                ,farray     =rn2d                     &
                                ,distgridToArrayMap=(/1/)             &
                                ,indexflag=ESMF_INDEX_GLOBAL          &
                                ,name       ='rn2d'                   &
                                ,rc         =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ! attach array to state
      MESSAGE_CHECK="Add rn2d array to import state"
      CALL ESMF_StateAdd(state=IMP_STATE                              &
                        ,array=TMP_ARRAY                              &
                        ,rc   =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      MESSAGE_CHECK="Create rn2d array for export state"
      ! create the ESMF_Field
      TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                ,farray     =rn2d                     &
                                ,distgridToArrayMap=(/1/)             &
                                ,indexflag=ESMF_INDEX_GLOBAL          &
                                ,name       ='rn2d'                   &
                                ,rc         =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ! attach array to state
      MESSAGE_CHECK="Add rn2d array to export state"
      CALL ESMF_StateAdd(state=EXP_STATE                              &
                        ,array=TMP_ARRAY                              &
                        ,rc   =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
!
! rc2d
!
      MESSAGE_CHECK="Create rc2d array for import state"
      ! create the ESMF_Field
      TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                ,farray     =rc2d                     &
                                ,distgridToArrayMap=(/1/)             &
                                ,indexflag=ESMF_INDEX_GLOBAL          &
                                ,name       ='rc2d'                   &
                                ,rc         =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ! attach array to state
      MESSAGE_CHECK="Add rc2d array to import state"
      CALL ESMF_StateAdd(state=IMP_STATE                              &
                        ,array=TMP_ARRAY                              &
                        ,rc   =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      MESSAGE_CHECK="Create rc2d array for export state"
      ! create the ESMF_Field
      TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                ,farray     =rc2d                     &
                                ,distgridToArrayMap=(/1/)             &
                                ,indexflag=ESMF_INDEX_GLOBAL          &
                                ,name       ='rc2d'                   &
                                ,rc         =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ! attach array to state
      MESSAGE_CHECK="Add rc2d array to export state"
      CALL ESMF_StateAdd(state=EXP_STATE                              &
                        ,array=TMP_ARRAY                              &
                        ,rc   =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
!
! ts2d
!
      MESSAGE_CHECK="Create ts2d array for import state"
      ! create the ESMF_Field
      TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                ,farray     =ts2d                     &
                                ,distgridToArrayMap=(/1/)             &
                                ,indexflag=ESMF_INDEX_GLOBAL          &
                                ,name       ='ts2d'                   &
                                ,rc         =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ! attach array to state
      MESSAGE_CHECK="Add ts2d array to import state"
      CALL ESMF_StateAdd(state=IMP_STATE                              &
                        ,array=TMP_ARRAY                              &
                        ,rc   =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      MESSAGE_CHECK="Create ts2d array for export state"
      ! create the ESMF_Field
      TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                ,farray     =ts2d                     &
                                ,distgridToArrayMap=(/1/)             &
                                ,indexflag=ESMF_INDEX_GLOBAL          &
                                ,name       ='ts2d'                   &
                                ,rc         =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ! attach array to state
      MESSAGE_CHECK="Add ts2d array to export state"
      CALL ESMF_StateAdd(state=EXP_STATE                              &
                        ,array=TMP_ARRAY                              &
                        ,rc   =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
!
! us2d
!
      MESSAGE_CHECK="Create us2d array for import state"
      ! create the ESMF_Field
      TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                ,farray     =us2d                     &
                                ,distgridToArrayMap=(/1/)             &
                                ,indexflag=ESMF_INDEX_GLOBAL          &
                                ,name       ='us2d'                   &
                                ,rc         =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ! attach array to state
      MESSAGE_CHECK="Add us2d array to import state"
      CALL ESMF_StateAdd(state=IMP_STATE                              &
                        ,array=TMP_ARRAY                              &
                        ,rc   =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      MESSAGE_CHECK="Create us2d array for export state"
      ! create the ESMF_Field
      TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                ,farray     =us2d                     &
                                ,distgridToArrayMap=(/1/)             &
                                ,indexflag=ESMF_INDEX_GLOBAL          &
                                ,name       ='us2d'                   &
                                ,rc         =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ! attach array to state
      MESSAGE_CHECK="Add us2d array to export state"
      CALL ESMF_StateAdd(state=EXP_STATE                              &
                        ,array=TMP_ARRAY                              &
                        ,rc   =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
!
! hf2d
!
      MESSAGE_CHECK="Create hf2d array for import state"
      ! create the ESMF_Field
      TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                ,farray     =hf2d                     &
                                ,distgridToArrayMap=(/1/)             &
                                ,indexflag=ESMF_INDEX_GLOBAL          &
                                ,name       ='hf2d'                   &
                                ,rc         =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ! attach array to state
      MESSAGE_CHECK="Add hf2d array to import state"
      CALL ESMF_StateAdd(state=IMP_STATE                              &
                        ,array=TMP_ARRAY                              &
                        ,rc   =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      MESSAGE_CHECK="Create hf2d array for export state"
      ! create the ESMF_Field
      TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                ,farray     =hf2d                     &
                                ,distgridToArrayMap=(/1/)             &
                                ,indexflag=ESMF_INDEX_GLOBAL          &
                                ,name       ='hf2d'                   &
                                ,rc         =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ! attach array to state
      MESSAGE_CHECK="Add hf2d array to export state"
      CALL ESMF_StateAdd(state=EXP_STATE                              &
                        ,array=TMP_ARRAY                              &
                        ,rc   =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
!
! qf2d
!
      MESSAGE_CHECK="Create qf2d array for import state"
      ! create the ESMF_Field
      TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                ,farray     =qf2d                     &
                                ,distgridToArrayMap=(/1/)             &
                                ,indexflag=ESMF_INDEX_GLOBAL          &
                                ,name       ='qf2d'                   &
                                ,rc         =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ! attach array to state
      MESSAGE_CHECK="Add qf2d array to import state"
      CALL ESMF_StateAdd(state=IMP_STATE                              &
                        ,array=TMP_ARRAY                              &
                        ,rc   =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      MESSAGE_CHECK="Create qf2d array for export state"
      ! create the ESMF_Field
      TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                ,farray     =qf2d                     &
                                ,distgridToArrayMap=(/1/)             &
                                ,indexflag=ESMF_INDEX_GLOBAL          &
                                ,name       ='qf2d'                   &
                                ,rc         =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ! attach array to state
      MESSAGE_CHECK="Add qf2d array to export state"
      CALL ESMF_StateAdd(state=EXP_STATE                              &
                        ,array=TMP_ARRAY                              &
                        ,rc   =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
!
! sheleg2d
!
      MESSAGE_CHECK="Create sheleg2d array for import state"
      ! create the ESMF_Field
      TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                ,farray     =sheleg2d                 &
                                ,distgridToArrayMap=(/1/)             &
                                ,indexflag=ESMF_INDEX_GLOBAL          &
                                ,name       ='sheleg2d'               &
                                ,rc         =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ! attach array to state
      MESSAGE_CHECK="Add sheleg2d array to import state"
      CALL ESMF_StateAdd(state=IMP_STATE                              &
                        ,array=TMP_ARRAY                              &
                        ,rc   =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      MESSAGE_CHECK="Create sheleg2d array for export state"
      ! create the ESMF_Field
      TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                ,farray     =sheleg2d                 &
                                ,distgridToArrayMap=(/1/)             &
                                ,indexflag=ESMF_INDEX_GLOBAL          &
                                ,name       ='sheleg2d'               &
                                ,rc         =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ! attach array to state
      MESSAGE_CHECK="Add sheleg2d array to export state"
      CALL ESMF_StateAdd(state=EXP_STATE                              &
                        ,array=TMP_ARRAY                              &
                        ,rc   =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
!
! canopy2d
!
      MESSAGE_CHECK="Create canopy2d array for import state"
      ! create the ESMF_Field
      TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                ,farray     =canopy2d                 &
                                ,distgridToArrayMap=(/1/)             &
                                ,indexflag=ESMF_INDEX_GLOBAL          &
                                ,name       ='canopy2d'               &
                                ,rc         =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ! attach array to state
      MESSAGE_CHECK="Add canopy2d array to import state"
      CALL ESMF_StateAdd(state=IMP_STATE                              &
                        ,array=TMP_ARRAY                              &
                        ,rc   =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      MESSAGE_CHECK="Create canopy2d array for export state"
      ! create the ESMF_Field
      TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                ,farray     =canopy2d                 &
                                ,distgridToArrayMap=(/1/)             &
                                ,indexflag=ESMF_INDEX_GLOBAL          &
                                ,name       ='canopy2d'               &
                                ,rc         =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ! attach array to state
      MESSAGE_CHECK="Add canopy2d array to export state"
      CALL ESMF_StateAdd(state=EXP_STATE                              &
                        ,array=TMP_ARRAY                              &
                        ,rc   =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
!
! hice2d
!
      MESSAGE_CHECK="Create hice2d array for import state"
      ! create the ESMF_Field
      TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                ,farray     =hice2d                   &
                                ,distgridToArrayMap=(/1/)             &
                                ,indexflag=ESMF_INDEX_GLOBAL          &
                                ,name       ='hice2d'                 &
                                ,rc         =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ! attach array to state
      MESSAGE_CHECK="Add hice2d array to import state"
      CALL ESMF_StateAdd(state=IMP_STATE                              &
                        ,array=TMP_ARRAY                              &
                        ,rc   =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      MESSAGE_CHECK="Create hice2d array for export state"
      ! create the ESMF_Field
      TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                ,farray     =hice2d                   &
                                ,distgridToArrayMap=(/1/)             &
                                ,indexflag=ESMF_INDEX_GLOBAL          &
                                ,name       ='hice2d'                 &
                                ,rc         =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ! attach array to state
      MESSAGE_CHECK="Add hice2d array to export state"
      CALL ESMF_StateAdd(state=EXP_STATE                              &
                        ,array=TMP_ARRAY                              &
                        ,rc   =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
!
! fice2d
!
      MESSAGE_CHECK="Create fice2d array for import state"
      ! create the ESMF_Field
      TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                ,farray     =fice2d                   &
                                ,distgridToArrayMap=(/1/)             &
                                ,indexflag=ESMF_INDEX_GLOBAL          &
                                ,name       ='fice2d'                 &
                                ,rc         =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ! attach array to state
      MESSAGE_CHECK="Add fice2d array to import state"
      CALL ESMF_StateAdd(state=IMP_STATE                              &
                        ,array=TMP_ARRAY                              &
                        ,rc   =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      MESSAGE_CHECK="Create fice2d array for export state"
      ! create the ESMF_Field
      TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                ,farray     =fice2d                   &
                                ,distgridToArrayMap=(/1/)             &
                                ,indexflag=ESMF_INDEX_GLOBAL          &
                                ,name       ='fice2d'                 &
                                ,rc         =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ! attach array to state
      MESSAGE_CHECK="Add fice2d array to export state"
      CALL ESMF_StateAdd(state=EXP_STATE                              &
                        ,array=TMP_ARRAY                              &
                        ,rc   =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
!
! st3d
!
      MESSAGE_CHECK="Create st3d array for import state"
      ! create the ESMF_Field
      TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                ,farray     =st3d                     &
                                ,distgridToArrayMap=(/2/)             &
                                ,indexflag=ESMF_INDEX_GLOBAL          &
                                ,name       ='st3d'                   &
                                ,rc         =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ! attach array to state
      MESSAGE_CHECK="Add st3d array to import state"
      CALL ESMF_StateAdd(state=IMP_STATE                              &
                        ,array=TMP_ARRAY                              &
                        ,rc   =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      MESSAGE_CHECK="Create st3d array for export state"
      ! create the ESMF_Field
      TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                ,farray     =st3d                     &
                                ,distgridToArrayMap=(/2/)             &
                                ,indexflag=ESMF_INDEX_GLOBAL          &
                                ,name       ='st3d'                   &
                                ,rc         =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ! attach array to state
      MESSAGE_CHECK="Add st3d array to export state"
      CALL ESMF_StateAdd(state=EXP_STATE                              &
                        ,array=TMP_ARRAY                              &
                        ,rc   =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
!
! sm3d
!
      MESSAGE_CHECK="Create sm3d array for import state"
      ! create the ESMF_Field
      TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                ,farray     =sm3d                     &
                                ,distgridToArrayMap=(/2/)             &
                                ,indexflag=ESMF_INDEX_GLOBAL          &
                                ,name       ='sm3d'                   &
                                ,rc         =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ! attach array to state
      MESSAGE_CHECK="Add sm3d array to import state"
      CALL ESMF_StateAdd(state=IMP_STATE                              &
                        ,array=TMP_ARRAY                              &
                        ,rc   =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      MESSAGE_CHECK="Create sm3d array for export state"
      ! create the ESMF_Field
      TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                ,farray     =sm3d                     &
                                ,distgridToArrayMap=(/2/)             &
                                ,indexflag=ESMF_INDEX_GLOBAL          &
                                ,name       ='sm3d'                   &
                                ,rc         =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ! attach array to state
      MESSAGE_CHECK="Add sm3d array to export state"
      CALL ESMF_StateAdd(state=EXP_STATE                              &
                        ,array=TMP_ARRAY                              &
                        ,rc   =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
!
! sw2d
!
      MESSAGE_CHECK="Create sw2d array for import state"
      ! create the ESMF_Field
      TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                ,farray     =sw2d                     &
                                ,distgridToArrayMap=(/1/)             &
                                ,indexflag=ESMF_INDEX_GLOBAL          &
                                ,name       ='sw2d'                   &
                                ,rc         =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ! attach array to state
      MESSAGE_CHECK="Add sw2d array to import state"
      CALL ESMF_StateAdd(state=IMP_STATE                              &
                        ,array=TMP_ARRAY                              &
                        ,rc   =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      MESSAGE_CHECK="Create sw2d array for export state"
      ! create the ESMF_Field
      TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                ,farray     =sw2d                     &
                                ,distgridToArrayMap=(/1/)             &
                                ,indexflag=ESMF_INDEX_GLOBAL          &
                                ,name       ='sw2d'                   &
                                ,rc         =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ! attach array to state
      MESSAGE_CHECK="Add sw2d array to export state"
      CALL ESMF_StateAdd(state=EXP_STATE                              &
                        ,array=TMP_ARRAY                              &
                        ,rc   =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
!
! lw2d
!
      MESSAGE_CHECK="Create lw2d array for import state"
      ! create the ESMF_Field
      TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                ,farray     =lw2d                     &
                                ,distgridToArrayMap=(/1/)             &
                                ,indexflag=ESMF_INDEX_GLOBAL          &
                                ,name       ='lw2d'                   &
                                ,rc         =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ! attach array to state
      MESSAGE_CHECK="Add lw2d array to import state"
      CALL ESMF_StateAdd(state=IMP_STATE                              &
                        ,array=TMP_ARRAY                              &
                        ,rc   =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      MESSAGE_CHECK="Create lw2d array for export state"
      ! create the ESMF_Field
      TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                ,farray     =lw2d                     &
                                ,distgridToArrayMap=(/1/)             &
                                ,indexflag=ESMF_INDEX_GLOBAL          &
                                ,name       ='lw2d'                   &
                                ,rc         =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ! attach array to state
      MESSAGE_CHECK="Add lw2d array to export state"
      CALL ESMF_StateAdd(state=EXP_STATE                              &
                        ,array=TMP_ARRAY                              &
                        ,rc   =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
!
! t2m2d
!
      MESSAGE_CHECK="Create t2m2d array for import state"
      ! create the ESMF_Field
      TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                ,farray     =t2m2d                    &
                                ,distgridToArrayMap=(/1/)             &
                                ,indexflag=ESMF_INDEX_GLOBAL          &
                                ,name       ='t2m2d'                  &
                                ,rc         =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ! attach array to state
      MESSAGE_CHECK="Add t2m2d array to import state"
      CALL ESMF_StateAdd(state=IMP_STATE                              &
                        ,array=TMP_ARRAY                              &
                        ,rc   =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      MESSAGE_CHECK="Create t2m2d array for export state"
      ! create the ESMF_Field
      TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                ,farray     =t2m2d                    &
                                ,distgridToArrayMap=(/1/)             &
                                ,indexflag=ESMF_INDEX_GLOBAL          &
                                ,name       ='t2m2d'                  &
                                ,rc         =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ! attach array to state
      MESSAGE_CHECK="Add t2m2d array to export state"
      CALL ESMF_StateAdd(state=EXP_STATE                              &
                        ,array=TMP_ARRAY                              &
                        ,rc   =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
!
! q2m2d
!
      MESSAGE_CHECK="Create q2m2d array for import state"
      ! create the ESMF_Field
      TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                ,farray     =q2m2d                    &
                                ,distgridToArrayMap=(/1/)             &
                                ,indexflag=ESMF_INDEX_GLOBAL          &
                                ,name       ='q2m2d'                  &
                                ,rc         =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ! attach array to state
      MESSAGE_CHECK="Add q2m2d array to import state"
      CALL ESMF_StateAdd(state=IMP_STATE                              &
                        ,array=TMP_ARRAY                              &
                        ,rc   =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      MESSAGE_CHECK="Create q2m2d array for export state"
      ! create the ESMF_Field
      TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                ,farray     =q2m2d                    &
                                ,distgridToArrayMap=(/1/)             &
                                ,indexflag=ESMF_INDEX_GLOBAL          &
                                ,name       ='q2m2d'                  &
                                ,rc         =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ! attach array to state
      MESSAGE_CHECK="Add q2m2d array to export state"
      CALL ESMF_StateAdd(state=EXP_STATE                              &
                        ,array=TMP_ARRAY                              &
                        ,rc   =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
!
! slmsk2d
!
      MESSAGE_CHECK="Create slmsk2d array for import state"
      ! create the ESMF_Field
      TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                ,farray     =slmsk2d                  &
                                ,distgridToArrayMap=(/1/)             &
                                ,indexflag=ESMF_INDEX_GLOBAL          &
                                ,name       ='slmsk2d'                &
                                ,rc         =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ! attach array to state
      MESSAGE_CHECK="Add slmsk2d array to import state"
      CALL ESMF_StateAdd(state=IMP_STATE                              &
                        ,array=TMP_ARRAY                              &
                        ,rc   =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      MESSAGE_CHECK="Create slmsk2d array for export state"
      ! create the ESMF_Field
      TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                ,farray     =slmsk2d                  &
                                ,distgridToArrayMap=(/1/)             &
                                ,indexflag=ESMF_INDEX_GLOBAL          &
                                ,name       ='slmsk2d'                &
                                ,rc         =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ! attach array to state
      MESSAGE_CHECK="Add slmsk2d array to export state"
      CALL ESMF_StateAdd(state=EXP_STATE                              &
                        ,array=TMP_ARRAY                              &
                        ,rc   =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
!
! hprm2d
!
      MESSAGE_CHECK="Create hprm2d array for import state"
      ! create the ESMF_Field
      TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                ,farray     =hprm2d                  &
                                ,distgridToArrayMap=(/2/)             &
                                ,indexflag=ESMF_INDEX_GLOBAL          &
                                ,name       ='hprm2d'                &
                                ,rc         =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ! attach array to state
      MESSAGE_CHECK="Add hprm2d array to import state"
      CALL ESMF_StateAdd(state=IMP_STATE                              &
                        ,array=TMP_ARRAY                              &
                        ,rc   =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      MESSAGE_CHECK="Create hprm2d array for export state"
      ! create the ESMF_Field
      TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                ,farray     =hprm2d                  &
                                ,distgridToArrayMap=(/2/)             &
                                ,indexflag=ESMF_INDEX_GLOBAL          &
                                ,name       ='hprm2d'                &
                                ,rc         =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ! attach array to state
      MESSAGE_CHECK="Add hprm2d array to export state"
      CALL ESMF_StateAdd(state=EXP_STATE                              &
                        ,array=TMP_ARRAY                              &
                        ,rc   =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)

!
! flxlwtoa2d
!
      MESSAGE_CHECK="Create flxlwtoa2d array for import state"
      ! create the ESMF_Field
      TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                ,farray     =flxlwtoa2d               &
                                ,distgridToArrayMap=(/1/)             &
                                ,indexflag=ESMF_INDEX_GLOBAL          &
                                ,name       ='flxlwtoa2d'             &
                                ,rc         =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ! attach array to state
      MESSAGE_CHECK="Add flxlwtoa2d array to import state"
      CALL ESMF_StateAdd(state=IMP_STATE                              &
                        ,array=TMP_ARRAY                              &
                        ,rc   =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      MESSAGE_CHECK="Create flxlwtoa2d array for export state"
      ! create the ESMF_Field
      TMP_ARRAY=ESMF_ArrayCreate(distgrid   =DISTGRID                 &
                                ,farray     =flxlwtoa2d               &
                                ,distgridToArrayMap=(/1/)             &
                                ,indexflag=ESMF_INDEX_GLOBAL          &
                                ,name       ='flxlwtoa2d'             &
                                ,rc         =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      ! attach array to state
      MESSAGE_CHECK="Add flxlwtoa2d array to export state"
      CALL ESMF_StateAdd(state=EXP_STATE                              &
                        ,array=TMP_ARRAY                              &
                        ,rc   =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)

!TBH:  validate states
      MESSAGE_CHECK="DYN_INITIALIZE:  Validate import state"
      call ESMF_StateValidate(state=IMP_STATE,rc=rc)
      IF(RC==ESMF_SUCCESS)THEN
!JR        WRITE(0,*)'DYN INITIALIZE import state valid'
      ENDIF
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
      MESSAGE_CHECK="DYN_INITIALIZE:  Validate export state"
      call ESMF_StateValidate(state=EXP_STATE,rc=rc)
      IF(RC==ESMF_SUCCESS)THEN
!JR        WRITE(0,*)'DYN INITIALIZE export state valid'
      ENDIF
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)

!
      IF(RC_INIT==ESMF_SUCCESS)THEN
!       WRITE(0,*)'DYN INITIALIZE STEP SUCCEEDED'
      ELSE
        WRITE(0,*)'DYN INITIALIZE STEP FAILED RC_INIT=',RC_INIT
      ENDIF
!
!-----------------------------------------------------------------------
!
#ifdef MANUALGPTL
      ret = gptlstop ('dyn_initialize')
#endif
      END SUBROUTINE DYN_INITIALIZE
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE DYN_RUN(DYN_GRID_COMP                                  &
                        ,IMP_STATE,EXP_STATE                            &
                        ,CLOCK_FIM                                      &
                        ,RC_RUN)
!
!-----------------------------------------------------------------------
!***  THE INTEGRATION OF THE MODEL DYNAMICS IS DONE
!***  THROUGH THIS ROUTINE.
!-----------------------------------------------------------------------
!
      USE module_fim_dyn_run ,only: DYN_RUN_FIM => dyn_run
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!***  ARGUMENT VARIABLES
!-----------------------------------------------------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: DYN_GRID_COMP                 !<-- The Dynamics gridded component
      TYPE(ESMF_State)   ,INTENT(INOUT) :: IMP_STATE                     !<-- The Dynamics import state
      TYPE(ESMF_State)   ,INTENT(INOUT) :: EXP_STATE                     !<-- The Dynamics export state
      TYPE(ESMF_Clock)   ,INTENT(IN)    :: CLOCK_FIM                     !<-- The FIM's ESMF Clock
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
#ifdef MANUALGPTL
      ret = gptlstart ('dyn_run')
#endif
      RC_RUN=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  EXTRACT THE TIMESTEP COUNT FROM THE CLOCK.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Retrieve Timestep from FIM Clock in Dynamics Run"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ClockGet(clock       =CLOCK_FIM                         &  !<-- The ESMF clock
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
!-----------------------------------------------------------------------
!***  THE MAIN DYNAMICS INTEGRATION LOOP.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      its = NTIMESTEP + 1
      CALL DYN_RUN_FIM (its)
!
!-----------------------------------------------------------------------
!
      RC=0
!
      IF(RC_RUN==ESMF_SUCCESS)THEN
!       WRITE(0,*)'DYN RUN STEP SUCCEEDED'
      ELSE
        WRITE(0,*)'DYN RUN STEP FAILED RC_RUN=',RC_RUN
      ENDIF
!
!-----------------------------------------------------------------------
!
#ifdef MANUALGPTL
      ret = gptlstop ('dyn_run')
#endif
      END SUBROUTINE DYN_RUN
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE DYN_FINALIZE(DYN_GRID_COMP                             &
                             ,IMP_STATE_WRITE                           &
                             ,EXP_STATE_WRITE                           &
                             ,CLOCK_FIM                                 &
                             ,RCFINAL)
!
!-----------------------------------------------------------------------
!***  FINALIZE THE DYNAMICS COMPONENT.
!-----------------------------------------------------------------------
!
      USE module_fim_dyn_finalize ,only: DYN_FINALIZE_FIM => dyn_finalize
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!***  ARGUMENT VARIABLES
!-----------------------------------------------------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: DYN_GRID_COMP                   !<-- The Dynamics gridded component
      TYPE(ESMF_State)   ,INTENT(INOUT) :: IMP_STATE_WRITE                 !<-- The Dynamics import state
      TYPE(ESMF_State),   INTENT(INOUT) :: EXP_STATE_WRITE                 !<-- The Dynamics export state
      TYPE(ESMF_Clock)   ,INTENT(INOUT) :: CLOCK_FIM                       !<-- The FIM component's ESMF Clock.
!
      INTEGER            ,INTENT(OUT)   :: RCFINAL
!      
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      INTEGER :: RC
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RCFINAL=ESMF_SUCCESS
!
      CALL DYN_FINALIZE_FIM
!
      ! destroy the grid created during DYN_INITIALIZE
      call ESMF_GridDestroy(GRID_FIM_DYN,rc=rc)
      CALL ERR_MSG(RC,"DYN_FINALIZE:  destroy GRID_FIM_DYN",RCFINAL)
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE DYN_FINALIZE


      END MODULE MODULE_DYNAMICS_GRID_COMP

