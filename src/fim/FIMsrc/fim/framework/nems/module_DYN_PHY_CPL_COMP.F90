!JR Copied from fimlatest
!TODO:  DRY out all of this code.  Initial NEMS was not DRY.  We can be.  
!-----------------------------------------------------------------------
!
      MODULE MODULE_DYN_PHY_CPL_COMP
!
!-----------------------------------------------------------------------
!
!***  THIS MODULE HOLDS THE COUPLER'S REGISTER, INIT, RUN, AND FINALIZE 
!***  ROUTINES.  THEY ARE CALLED FROM THE FIM GRIDDED COMPONENT
!***  IN module_FIM_GRID_COMP.F90.
!
!***  THE COUPLER PROVIDES 2-WAY COUPLING BETWEEN THE DYNAMICS AND
!***  PHYSICS GRIDDED COMPONENTS BY TRANSFERING THEIR EXPORT AND
!***  IMPORT STATES BETWEEN THE TWO.
!
!-----------------------------------------------------------------------
!
      USE ESMF_MOD
!
      USE MODULE_ERR_MSG,ONLY: ERR_MSG,MESSAGE_CHECK
!JR module_nems_share doesn't seem to exist anymore. Comment out since we may go to GPTL anyway
!JR      USE MODULE_NEMS_SHARE,ONLY : TIMEF
      USE MACHINE         ,only: kind_evod,kind_phys,kind_rad
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      PRIVATE
!
      PUBLIC :: DYN_PHY_CPL_REGISTER
      public :: cpl_initialize   !JR Make this public so can call directly
!
!-----------------------------------------------------------------------
!
      REAL*8 :: btim0
      REAL*8, PUBLIC :: cpl_dyn_phy_tim
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------

! Pointers to arrays to be coupled.  These pointers are set to point to 
! fields in corresponding states prior to use.  

!-----------------------------------------------------------------------
! FIM import state
!-----------------------------------------------------------------------
real,pointer :: fim_imp_us3d(:,:)   ! FIM zonal wind (m/s), layer
real,pointer :: fim_imp_vs3d(:,:)   ! FIM meridional wind (m/s), layer
real,pointer :: fim_imp_tr3d(:,:,:) ! FIM tracers
real,pointer :: fim_imp_rn2d(:)     ! FIM accumulated total precipitation/rainfall
real,pointer :: fim_imp_rc2d(:)     ! FIM accumulated convective precipitation/rainfall
real,pointer :: fim_imp_ts2d(:)     ! FIM skin temperature
real,pointer :: fim_imp_us2d(:)     ! FIM friction velocity/equivalent momentum flux
real,pointer :: fim_imp_hf2d(:)     ! FIM sensible heat flux
real,pointer :: fim_imp_qf2d(:)     ! FIM water vapor/equivalent latent heat flux
real,pointer :: fim_imp_sheleg2d(:)
real,pointer :: fim_imp_canopy2d(:)
real,pointer :: fim_imp_hice2d(:)
real,pointer :: fim_imp_fice2d(:)
real,pointer :: fim_imp_st3d(:,:)   ! FIM soil temperature
real,pointer :: fim_imp_sm3d(:,:)   ! FIM soil moisture
real,pointer :: fim_imp_sw2d(:)     ! FIM downward short-wave radiation flux
real,pointer :: fim_imp_lw2d(:)     ! FIM downward long-wave radiation flux
real,pointer :: fim_imp_t2m2d(:)    ! FIM 2-meter temp.
real,pointer :: fim_imp_q2m2d(:)    ! FIM 2-meter spfh
real,pointer :: fim_imp_slmsk2d(:)
real,pointer :: fim_imp_hprm2d(:,:) ! FIM soil temperature
real,pointer :: fim_imp_flxlwtoa2d(:)
!-----------------------------------------------------------------------
! FIM export state
!-----------------------------------------------------------------------
real,pointer :: fim_exp_pr3d(:,:)   ! FIM pressure (pascal)
real,pointer :: fim_exp_us3d(:,:)   ! FIM zonal wind (m/s), layer
real,pointer :: fim_exp_vs3d(:,:)   ! FIM meridional wind (m/s), layer
real,pointer :: fim_exp_ws3d(:,:)   ! FIM vertical   wind (m/s), layer
real,pointer :: fim_exp_tr3d(:,:,:) ! FIM tracers
!-----------------------------------------------------------------------
! GFS import state
!-----------------------------------------------------------------------
real(kind=kind_evod) ,pointer :: gfs_imp_ps(:)
real(kind=kind_evod) ,pointer :: gfs_imp_dp(:,:)
real(kind=kind_evod) ,pointer :: gfs_imp_p(:,:)
real(kind=kind_evod) ,pointer :: gfs_imp_u(:,:)
real(kind=kind_evod) ,pointer :: gfs_imp_v(:,:)
real(kind=kind_evod) ,pointer :: gfs_imp_dpdt(:,:)
real(kind=kind_evod) ,pointer :: gfs_imp_q(:,:)
real(kind=kind_evod) ,pointer :: gfs_imp_oz(:,:)
real(kind=kind_evod) ,pointer :: gfs_imp_cld(:,:)
real(kind=kind_evod) ,pointer :: gfs_imp_t(:,:)
!-----------------------------------------------------------------------
! GFS export state
!-----------------------------------------------------------------------
real(kind=kind_evod) ,pointer :: gfs_exp_p(:,:)
real(kind=kind_evod) ,pointer :: gfs_exp_u(:,:)
real(kind=kind_evod) ,pointer :: gfs_exp_v(:,:)
real(kind=kind_evod) ,pointer :: gfs_exp_q(:,:)
real(kind=kind_evod) ,pointer :: gfs_exp_cld(:,:)
real(kind=kind_evod) ,pointer :: gfs_exp_t(:,:)
real(kind=kind_phys) ,pointer :: gfs_exp_geshem(:,:)
real(kind=kind_phys) ,pointer :: gfs_exp_rainc(:,:)
real(kind=kind_phys) ,pointer :: gfs_exp_tsea(:,:)
real(kind=kind_phys) ,pointer :: gfs_exp_uustar(:,:)
real(kind=kind_phys) ,pointer :: gfs_exp_hflx(:,:)
real(kind=kind_phys) ,pointer :: gfs_exp_evap(:,:)
real(kind=kind_phys) ,pointer :: gfs_exp_sheleg(:,:)
real(kind=kind_phys) ,pointer :: gfs_exp_canopy(:,:)
real(kind=kind_phys) ,pointer :: gfs_exp_hice(:,:)
real(kind=kind_phys) ,pointer :: gfs_exp_fice(:,:)
real(kind=kind_phys) ,pointer :: gfs_exp_stc(:,:,:)
real(kind=kind_phys) ,pointer :: gfs_exp_smc(:,:,:)
real(kind=kind_phys) ,pointer :: gfs_exp_sfcdsw(:,:)
real(kind=kind_phys) ,pointer :: gfs_exp_sfcdlw(:,:)
real(kind=kind_phys) ,pointer :: gfs_exp_t2m(:,:)
real(kind=kind_phys) ,pointer :: gfs_exp_q2m(:,:)
real(kind=kind_phys) ,pointer :: gfs_exp_slmsk(:,:)
real(kind=kind_rad)  ,pointer :: gfs_exp_hprime(:,:,:)
real(kind=kind_rad)  ,pointer :: gfs_exp_fluxr(:,:,:)

!-----------------------------------------------------------------------

      CONTAINS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE DYN_PHY_CPL_REGISTER(CPL_COMP,IRC_REG)
!
!-----------------------------------------------------------------------
!***  REGISTER THE COUPLER COMPONENT'S INITIALIZE, RUN, AND FINALIZE
!***  ROUTINES.
!-----------------------------------------------------------------------
!
      TYPE(ESMF_CplComp),INTENT(INOUT) :: CPL_COMP ! Coupler component
!
      INTEGER,INTENT(OUT) :: IRC_REG               ! Return code for register
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      INTEGER :: IRC=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      IRC_REG=ESMF_SUCCESS  ! The error signal variable
                                                                                                                                              
!-----------------------------------------------------------------------
!***  REGISTER THE COUPLER INITIALIZE SUBROUTINE.  SINCE IT IS JUST ONE
!***  SUBROUTINE, USE ESMF_SINGLEPHASE.  THE SECOND ARGUMENT IS
!***  A PRE-DEFINED SUBROUTINE TYPE, SUCH AS ESMF_SETINIT, ESMF_SETRUN,
!***  OR ESMF_SETFINAL.
!-----------------------------------------------------------------------
!
      CALL ESMF_LogWrite("Set Entry Point for Coupler Initialize"       &
                        ,ESMF_LOG_INFO,RC=IRC)
!
      CALL ESMF_CplCompSetEntryPoint(CPL_COMP                           &  !<-- cplcomp
                                    ,ESMF_SETINIT                       &  !<-- subroutineType
                                    ,CPL_INITIALIZE                     &  !<-- user's subroutineName
                                    ,ESMF_SINGLEPHASE                   &  !<-- phase
                                    ,IRC)
!
      IF(ESMF_LogMsgFoundError(IRC,"Set Entry Point for Coupler Initialize"))THEN
        IRC_REG=ESMF_FAILURE
        WRITE(0,*)'Error Setting the Entry Point for Coupler Initialize, RC =',IRC
        IRC=ESMF_SUCCESS
      ENDIF
!
!-----------------------------------------------------------------------
!***  REGISTER THE COUPLER RUN SUBROUTINE.
!-----------------------------------------------------------------------
!
      CALL ESMF_LogWrite("Set Entry Point for Coupler Run"              &
                        ,ESMF_LOG_INFO,RC=IRC)
!
      CALL ESMF_CplCompSetEntryPoint(CPL_COMP                           &  !<-- cplcomp
                                    ,ESMF_SETRUN                        &  !<-- subroutineType
                                    ,CPL_RUN                            &  !<-- user's subroutineName
                                    ,ESMF_SINGLEPHASE                   &  !<-- phase
                                    ,IRC)
!
      IF(ESMF_LogMsgFoundError(IRC,"Set Entry Point for Coupler Run"))THEN
        IRC_REG=ESMF_FAILURE
        WRITE(0,*)'Error Setting the Entry Point for Coupler Run, RC =',IRC
        IRC=ESMF_SUCCESS
      ENDIF
!
!-----------------------------------------------------------------------
!***  REGISTER THE COUPLER FINALIZE SUBROUTINE.
!-----------------------------------------------------------------------
!
      CALL ESMF_LogWrite("Set Entry Point for Coupler Finalize"         &
                        ,ESMF_LOG_INFO,RC=IRC)
!
      CALL ESMF_CplCompSetEntryPoint(CPL_COMP                           &  !<-- cplcomp
                                    ,ESMF_SETFINAL                      &  !<-- subroutineType
                                    ,CPL_FINALIZE                       &  !<-- user's subroutineName
                                    ,ESMF_SINGLEPHASE                   &  !<-- phase
                                    ,IRC)
!
      IF(ESMF_LogMsgFoundError(IRC,"Set Entry Point for Coupler Finalize"))THEN
        IRC_REG=ESMF_FAILURE
        WRITE(0,*)'Error Setting the Entry Point for Coupler Finalize, RC =',IRC
      ENDIF
!
!-----------------------------------------------------------------------
!***  CHECK THE ERROR SIGNAL VARIABLE.
!-----------------------------------------------------------------------
!
      IF(IRC_REG==ESMF_SUCCESS)THEN
!       WRITE(0,*)" COUPLER_REGISTER SUCCEEDED"
      ELSE
        WRITE(0,*)" COUPLER_REGISTER FAILED"
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE DYN_PHY_CPL_REGISTER
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE CPL_INITIALIZE(CPL_COMP,IMP_STATE,EXP_STATE,CLOCK      &
                               ,IRC_CPL)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  SET UP THE COUPLER.
!-----------------------------------------------------------------------
!
      ! FIM dynamics + GFS physics
      USE module_fim_cpl_init ,only: CPL_INITIALIZE_FIM_GFS => cpl_init
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!***  ARGUMENT VARIABLES.
!-----------------------------------------------------------------------
!
      TYPE(ESMF_CplComp),INTENT(INOUT) :: CPL_COMP
      TYPE(ESMF_State),  INTENT(INOUT) :: IMP_STATE
      TYPE(ESMF_State),  INTENT(INOUT) :: EXP_STATE
      TYPE(ESMF_Clock),  INTENT(IN)    :: CLOCK
!
      INTEGER,           INTENT(OUT)   :: IRC_CPL
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      INTEGER :: IRCFINAL
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      cpl_dyn_phy_tim=0.
!-----------------------------------------------------------------------
!***  INITIALIZE THE ERROR SIGNAL VARIABLES.
!-----------------------------------------------------------------------
!
      IRCFINAL=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  INITIALIZE THE PHYSICS SCHEMES. 
!-----------------------------------------------------------------------
!
!TODO:  Note that the latest FIM does delegate creation of the ESMF_Grid to 
!TODO:  the DYN component.  This grid is then used by the FIM and PHY 
!TODO:  components.  However, coordinates are not set in the ESMF_Grid so 
!TODO:  we still need a better approach for handling lat/lon.  
!TODO:
!TODO:  We appear to have three cases of interest for handling lat,lon during 
!TODO:  component initialization:
!TODO:    A)  PHY runs on the same lat,lon grid as DYN and FIM.  ATM also uses 
!TODO:        this grid (i.e. for coupling with OCN, etc.).
!TODO:    B)  PHY and DYN run on independent grids.  ATM/FIM and DYN share
!TODO:        the same grid.
!TODO:    C)  PHY and DYN run on independent grids.  ATM/FIM and PHY share
!TODO:        the same grid.  
!TODO:  NMMB and FIM use case "A".  GFS also uses "A" for the physical-space 
!TODO:  grid.  Spectral grids appear to be handled internally to GFS DYN.  
!TODO:  "B" and "C" are impractical due to the excessive cost of re-gridding 
!TODO:  full 3D arrays every time step.  
!TODO:  
!TODO:  For FIM, we'd like to have the option of using case "A" with the NCEP 
!TODO:  GFS physics component.  That way we could grab the latest GFS PHY 
!TODO:  component and hook it up without going into 
!TODO:  gfs_phy_initialize()->gfs_physics_initialize()->fix_fields()->
!TODO:  LONLAT_PARA() and inserting our own lat,lon computation (or performing 
!TODO:  some other ugly hackery.)  It would make it easier to swap GFS PHY 
!TODO:  between GFS DYN and FIM DYN.
!TODO:  
!TODO:  Logic of ESMF_Grid creation in the FIM component now proceeds as 
!TODO:  follows:  
!TODO:    1)  The ATM component does not attach an ESMF_Grid to the FIM 
!TODO:        component.  
!TODO:    2)  The FIM component does not attach an ESMF_Grid to the DYN 
!TODO:        component.  
!TODO:    3)  The DYN component creates an ESMF_Grid and attaches it to 
!TODO:        itself.  DYN does *not* yet fill in lat,lon coordinates.  
!TODO:    4)  The FIM extracts the ESMF_Grid from the DYN component after 
!TODO:        dyn_initialize() and attaches it to itself and to the PHY 
!TODO:        component.  
!TODO:    5)  The ATM component does not yet use the ESMF_Grid available 
!TODO:        via the FIM component.  
!TODO:
!TODO:  Since this pattern is not yet implemented in NEMS, we pass lat,lon 
!TODO:  from FIM DYN to PHY in cpl_init.F90 via cpl_init_dyn_to_phy().  We 
!TODO:  need to extend step #3 to include initialization of ESMF_Grid 
!TODO:  coordinates to pass lat/lon via the ESMF_Grid instead of via 
!TODO:  ESMF_State objects.  Also, step #5 needs to be addressed in the 
!TODO:  NCEP code to allow ATM to use the ESMF_Grid attached to the 
!TODO:  FIM/NMMB/GFS component.  
!TODO:  
!TODO:  Also there is a desire to write some of the PHY variables to the FIM 
!TODO:  history output stream during the first (0h) write.  This is currently 
!TODO:  handled by passing these fields from PHY to DYN via a call to 
!TODO:  cpl_init_phy_to_dyn().  Future NEMS I/O component(s) may allow this 
!TODO:  to be handled in another way.  
!TODO:
      CALL CPL_INITIALIZE_FIM_GFS
!
!-----------------------------------------------------------------------
!
      IF(IRCFINAL==ESMF_SUCCESS)THEN
!       WRITE(0,*)"CPL INITIALIZE STEP SUCCEEDED"
      ELSE
        WRITE(0,*)"CPL INITIALIZE STEP FAILED"
      ENDIF
!
      IRC_CPL=IRCFINAL
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE CPL_INITIALIZE
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE CPL_RUN(CPL_COMP,IMP_STATE,EXP_STATE,CLOCK,IRC_CPL)
!
!-----------------------------------------------------------------------
!***  RUN THE COUPLER TO TRANSFER DATA BETWEEN THE GRIDDED COMPONENTS.
!-----------------------------------------------------------------------
!
      ! FIM dynamics + GFS physics
      USE module_fim_cpl_run ,only: CPL_DYN_TO_PHY, CPL_PHY_TO_DYN

!TODO:  REMOVE THIS use-association!  Refactor or replace with ESMF_Alarms
      use module_control  ,only: nts,CallPhysics

!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!***  ARGUMENT VARIABLES
!-----------------------------------------------------------------------
!
      TYPE(ESMF_CplComp),INTENT(INOUT) :: CPL_COMP
      TYPE(ESMF_State),  INTENT(INOUT) :: IMP_STATE
      TYPE(ESMF_State),  INTENT(INOUT) :: EXP_STATE
      TYPE(ESMF_Clock),  INTENT(IN)    :: CLOCK
!
      INTEGER,           INTENT(OUT)   :: IRC_CPL
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      INTEGER :: IRCFINAL
      INTEGER :: NTIMESTEP,RC
      INTEGER(KIND=ESMF_KIND_I8) :: NTIMESTEP_ESMF
      INTEGER :: its
      character(esmf_maxstr) :: import_statename
      character(esmf_maxstr) :: export_statename
      ! temporary field object
!      TYPE(ESMF_Field) :: TMP_FIELD
      TYPE(ESMF_Array) :: TMP_ARRAY
!
!     TYPE(ESMF_RouteHandle) :: ROUTEHANDLE
!
!-----------------------------------------------------------------------
!***********************************************************************
!
!JR      btim0=timef()
!-----------------------------------------------------------------------
!***  INITIALIZE THE ERROR SIGNAL VARIABLES.
!-----------------------------------------------------------------------
!
      IRCFINAL=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  EXTRACT THE TIMESTEP COUNT FROM THE CLOCK.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Retrieve Timestep from FIM Clock in Physics Run"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ClockGet(clock       =CLOCK                             &  !<-- The ESMF clock
                        ,advanceCount=NTIMESTEP_ESMF                    &  !<-- The number of times the clock has been advanced
                        ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      NTIMESTEP=NTIMESTEP_ESMF
!
!-----------------------------------------------------------------------
!***  COUPLE from DYN->PHY or PHY->DYN
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  determine the direction of the transfer by extracting
!***  the statename from the import state.
!-----------------------------------------------------------------------
!
      call esmf_logwrite("retrieve state name in coupler"               &
                         ,esmf_log_info,rc=rc)
!
      call esmf_stateget(imp_state                                      &
                        ,name =import_statename                         &
                        ,rc   =rc)
      call esmf_stateget(exp_state                                      &
                        ,name =export_statename                         &
                        ,rc   =rc)
!
!  print *,'CPL_RUN: move data from(',trim(import_statename), ') to (', &
!                                     trim(export_statename), ')'
!
      its = NTIMESTEP + 1
!TODO:  Replace all of this "its" stuff with ESMF_Alarms.  
!TODO:  Then eliminate dependence on nts and CallPhysics.  
!TODO:  Note that GFS has no concept of not calling physics every time step!  
!TODO:  OR, shove "if" statements down into CPL_DYN_TO_PHY and CPL_PHY_TO_DYN 
      if (its <= nts ) then
        !TODO:  Eliminate duplication by encapsulating this logic
        if(mod(its,CallPhysics)==0.or.its==1) then ! Do physics
          ! Note that state names are set in fim_initialize()
          if      ( (trim(import_statename).eq.'FIM dynamics export') .and. &
                    (trim(export_statename).eq.'FIM physics import') ) then
!
! extract pointers from ESMF_States and stuff in fim_* and gfs_* pointers
!
            MESSAGE_CHECK="Get pr3d array from import state"
            CALL ESMF_StateGet(state    =IMP_STATE                      &
                              ,itemName ='pr3d'                         &
                              ,array    =TMP_ARRAY                      &
                              ,rc       =RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)
            MESSAGE_CHECK="Get pr3d pointer from imported array"
!TBH:  Note that the following calls with named parameters causes ifort to 
!TBH:  complain whine that "There is no matching specific subroutine for 
!TBH:  this generic subroutine call."  Removing the first three dummy 
!TBH:  argument named makes ifort happy.  
!            CALL ESMF_FieldGet(field    =TMP_FIELD                      &
!                              ,localDe=0                                &
!                              ,farrayPtr=fim_exp_pr3d                   &
!                              ,rc       =RC)
!TBH:  Futher note that ESMF_Field did not work for reasons not yet known.  
!TBH:  (See comments in module_DYNAMICS_GRID_COMP.F90 and 
!TBH:  module_PHYSICS_GRID_COMP.F90).  Switched to ESMF_Array until this 
!TBH:  problem is resolved.  
!            CALL ESMF_FieldGet(TMP_FIELD,0,fim_exp_pr3d,rc=RC)
!TODO:  Switch back to ESMF_Field since future NEMS will use ESMF_Fields.  
            CALL ESMF_ArrayGet(TMP_ARRAY,0,fim_exp_pr3d,rc=RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)

            MESSAGE_CHECK="Get us3d array from import state"
            CALL ESMF_StateGet(state    =IMP_STATE                      &
                              ,itemName ='us3d'                         &
                              ,array    =TMP_ARRAY                      &
                              ,rc       =RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)
            MESSAGE_CHECK="Get us3d pointer from imported array"
            CALL ESMF_ArrayGet(TMP_ARRAY,0,fim_exp_us3d,rc=RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)

            MESSAGE_CHECK="Get vs3d array from import state"
            CALL ESMF_StateGet(state    =IMP_STATE                      &
                              ,itemName ='vs3d'                         &
                              ,array    =TMP_ARRAY                      &
                              ,rc       =RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)
            MESSAGE_CHECK="Get vs3d pointer from imported array"
            CALL ESMF_ArrayGet(TMP_ARRAY,0,fim_exp_vs3d,rc=RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)

            MESSAGE_CHECK="Get ws3d array from import state"
            CALL ESMF_StateGet(state    =IMP_STATE                      &
                              ,itemName ='ws3d'                         &
                              ,array    =TMP_ARRAY                      &
                              ,rc       =RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)
            MESSAGE_CHECK="Get ws3d pointer from imported array"
            CALL ESMF_ArrayGet(TMP_ARRAY,0,fim_exp_ws3d,rc=RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)

            MESSAGE_CHECK="Get tr3d array from import state"
            CALL ESMF_StateGet(state    =IMP_STATE                      &
                              ,itemName ='tr3d'                         &
                              ,array    =TMP_ARRAY                      &
                              ,rc       =RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)
            MESSAGE_CHECK="Get tr3d pointer from imported array"
            CALL ESMF_ArrayGet(TMP_ARRAY,0,fim_exp_tr3d,rc=RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)

            MESSAGE_CHECK="Get ps array from export state"
            CALL ESMF_StateGet(state    =EXP_STATE                      &
                              ,itemName ='ps'                           &
                              ,array    =TMP_ARRAY                      &
                              ,rc       =RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)
            MESSAGE_CHECK="Get ps pointer from exported array"
            CALL ESMF_ArrayGet(TMP_ARRAY,0,gfs_imp_ps,rc=RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)

            MESSAGE_CHECK="Get t array from export state"
            CALL ESMF_StateGet(state    =EXP_STATE                      &
                              ,itemName ='t'                            &
                              ,array    =TMP_ARRAY                      &
                              ,rc       =RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)
            MESSAGE_CHECK="Get t pointer from exported array"
            CALL ESMF_ArrayGet(TMP_ARRAY,0,gfs_imp_t,rc=RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)

            MESSAGE_CHECK="Get u array from export state"
            CALL ESMF_StateGet(state    =EXP_STATE                      &
                              ,itemName ='u'                            &
                              ,array    =TMP_ARRAY                      &
                              ,rc       =RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)
            MESSAGE_CHECK="Get u pointer from exported array"
            CALL ESMF_ArrayGet(TMP_ARRAY,0,gfs_imp_u,rc=RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)

            MESSAGE_CHECK="Get v array from export state"
            CALL ESMF_StateGet(state    =EXP_STATE                      &
                              ,itemName ='v'                            &
                              ,array    =TMP_ARRAY                      &
                              ,rc       =RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)
            MESSAGE_CHECK="Get v pointer from exported array"
            CALL ESMF_ArrayGet(TMP_ARRAY,0,gfs_imp_v,rc=RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)

            MESSAGE_CHECK="Get q array from export state"
            CALL ESMF_StateGet(state    =EXP_STATE                      &
                              ,itemName ='shum'                         &
                              ,array    =TMP_ARRAY                      &
                              ,rc       =RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)
            MESSAGE_CHECK="Get q pointer from exported array"
            CALL ESMF_ArrayGet(TMP_ARRAY,0,gfs_imp_q,rc=RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)

            MESSAGE_CHECK="Get oz array from export state"
            CALL ESMF_StateGet(state    =EXP_STATE                      &
                              ,itemName ='oz'                           &
                              ,array    =TMP_ARRAY                      &
                              ,rc       =RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)
            MESSAGE_CHECK="Get oz pointer from exported array"
            CALL ESMF_ArrayGet(TMP_ARRAY,0,gfs_imp_oz,rc=RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)

            MESSAGE_CHECK="Get cld array from export state"
            CALL ESMF_StateGet(state    =EXP_STATE                      &
                              ,itemName ='cld'                          &
                              ,array    =TMP_ARRAY                      &
                              ,rc       =RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)
            MESSAGE_CHECK="Get cld pointer from exported array"
            CALL ESMF_ArrayGet(TMP_ARRAY,0,gfs_imp_cld,rc=RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)

            MESSAGE_CHECK="Get p array from export state"
            CALL ESMF_StateGet(state    =EXP_STATE                      &
                              ,itemName ='p'                            &
                              ,array    =TMP_ARRAY                      &
                              ,rc       =RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)
            MESSAGE_CHECK="Get p pointer from exported array"
            CALL ESMF_ArrayGet(TMP_ARRAY,0,gfs_imp_p,rc=RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)

            MESSAGE_CHECK="Get dp array from export state"
            CALL ESMF_StateGet(state    =EXP_STATE                      &
                              ,itemName ='dp'                           &
                              ,array    =TMP_ARRAY                      &
                              ,rc       =RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)
            MESSAGE_CHECK="Get dp pointer from exported array"
            CALL ESMF_ArrayGet(TMP_ARRAY,0,gfs_imp_dp,rc=RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)

            MESSAGE_CHECK="Get dpdt array from export state"
            CALL ESMF_StateGet(state    =EXP_STATE                      &
                              ,itemName ='dpdt'                         &
                              ,array    =TMP_ARRAY                      &
                              ,rc       =RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)
            MESSAGE_CHECK="Get dpdt pointer from exported array"
            CALL ESMF_ArrayGet(TMP_ARRAY,0,gfs_imp_dpdt,rc=RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)

            ! couple from dynamics to physics
            CALL CPL_DYN_TO_PHY(its,                                        &
                   ! IN args
                   fim_exp_pr3d, fim_exp_us3d, fim_exp_vs3d, fim_exp_ws3d,  &
                   fim_exp_tr3d,                                            &
                   ! OUT args
                   gfs_imp_ps, gfs_imp_dp, gfs_imp_p, gfs_imp_u, gfs_imp_v, &
                   gfs_imp_dpdt, gfs_imp_q,  gfs_imp_oz, gfs_imp_cld,       &
                   gfs_imp_t)
          else if ( (trim(import_statename).eq.'FIM physics export')  .and. &
                    (trim(export_statename).eq.'FIM dynamics import') ) then
!
! extract pointers from ESMF_States and stuff in fim_* and gfs_* pointers
!
            MESSAGE_CHECK="Get p array from import state"
            CALL ESMF_StateGet(state    =IMP_STATE                      &
                              ,itemName ='p'                            &
                              ,array    =TMP_ARRAY                      &
                              ,rc       =RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)
            MESSAGE_CHECK="Get p pointer from imported array"
            CALL ESMF_ArrayGet(TMP_ARRAY,0,gfs_exp_p,rc=RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)

            MESSAGE_CHECK="Get u array from import state"
            CALL ESMF_StateGet(state    =IMP_STATE                      &
                              ,itemName ='u'                            &
                              ,array    =TMP_ARRAY                      &
                              ,rc       =RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)
            MESSAGE_CHECK="Get u pointer from imported array"
            CALL ESMF_ArrayGet(TMP_ARRAY,0,gfs_exp_u,rc=RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)

            MESSAGE_CHECK="Get v array from import state"
            CALL ESMF_StateGet(state    =IMP_STATE                      &
                              ,itemName ='v'                            &
                              ,array    =TMP_ARRAY                      &
                              ,rc       =RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)
            MESSAGE_CHECK="Get v pointer from imported array"
            CALL ESMF_ArrayGet(TMP_ARRAY,0,gfs_exp_v,rc=RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)

            MESSAGE_CHECK="Get q array from import state"
            CALL ESMF_StateGet(state    =IMP_STATE                      &
                              ,itemName ='shum'                         &
                              ,array    =TMP_ARRAY                      &
                              ,rc       =RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)
            MESSAGE_CHECK="Get q pointer from imported array"
            CALL ESMF_ArrayGet(TMP_ARRAY,0,gfs_exp_q,rc=RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)

            MESSAGE_CHECK="Get cld array from import state"
            CALL ESMF_StateGet(state    =IMP_STATE                      &
                              ,itemName ='cld'                          &
                              ,array    =TMP_ARRAY                      &
                              ,rc       =RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)
            MESSAGE_CHECK="Get cld pointer from imported array"
            CALL ESMF_ArrayGet(TMP_ARRAY,0,gfs_exp_cld,rc=RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)

            MESSAGE_CHECK="Get t array from import state"
            CALL ESMF_StateGet(state    =IMP_STATE                      &
                              ,itemName ='t'                            &
                              ,array    =TMP_ARRAY                      &
                              ,rc       =RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)
            MESSAGE_CHECK="Get t pointer from imported array"
            CALL ESMF_ArrayGet(TMP_ARRAY,0,gfs_exp_t,rc=RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)

            MESSAGE_CHECK="Get geshem array from import state"
            CALL ESMF_StateGet(state    =IMP_STATE                      &
                              ,itemName ='geshem'                       &
                              ,array    =TMP_ARRAY                      &
                              ,rc       =RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)
            MESSAGE_CHECK="Get geshem pointer from imported array"
            CALL ESMF_ArrayGet(TMP_ARRAY,0,gfs_exp_geshem,rc=RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)

            MESSAGE_CHECK="Get rainc array from import state"
            CALL ESMF_StateGet(state    =IMP_STATE                      &
                              ,itemName ='rainc'                        &
                              ,array    =TMP_ARRAY                      &
                              ,rc       =RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)
            MESSAGE_CHECK="Get rainc pointer from imported array"
            CALL ESMF_ArrayGet(TMP_ARRAY,0,gfs_exp_rainc,rc=RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)

            MESSAGE_CHECK="Get tsea array from import state"
            CALL ESMF_StateGet(state    =IMP_STATE                      &
                              ,itemName ='tsea'                         &
                              ,array    =TMP_ARRAY                      &
                              ,rc       =RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)
            MESSAGE_CHECK="Get tsea pointer from imported array"
            CALL ESMF_ArrayGet(TMP_ARRAY,0,gfs_exp_tsea,rc=RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)

            MESSAGE_CHECK="Get uustar array from import state"
            CALL ESMF_StateGet(state    =IMP_STATE                      &
                              ,itemName ='uustar'                       &
                              ,array    =TMP_ARRAY                      &
                              ,rc       =RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)
            MESSAGE_CHECK="Get uustar pointer from imported array"
            CALL ESMF_ArrayGet(TMP_ARRAY,0,gfs_exp_uustar,rc=RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)

            MESSAGE_CHECK="Get hflx array from import state"
            CALL ESMF_StateGet(state    =IMP_STATE                      &
                              ,itemName ='hflx'                         &
                              ,array    =TMP_ARRAY                      &
                              ,rc       =RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)
            MESSAGE_CHECK="Get hflx pointer from imported array"
            CALL ESMF_ArrayGet(TMP_ARRAY,0,gfs_exp_hflx,rc=RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)

            MESSAGE_CHECK="Get evap array from import state"
            CALL ESMF_StateGet(state    =IMP_STATE                      &
                              ,itemName ='evap'                         &
                              ,array    =TMP_ARRAY                      &
                              ,rc       =RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)
            MESSAGE_CHECK="Get evap pointer from imported array"
            CALL ESMF_ArrayGet(TMP_ARRAY,0,gfs_exp_evap,rc=RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)

            MESSAGE_CHECK="Get sheleg array from import state"
            CALL ESMF_StateGet(state    =IMP_STATE                      &
                              ,itemName ='sheleg'                       &
                              ,array    =TMP_ARRAY                      &
                              ,rc       =RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)
            MESSAGE_CHECK="Get sheleg pointer from imported array"
            CALL ESMF_ArrayGet(TMP_ARRAY,0,gfs_exp_sheleg,rc=RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)

            MESSAGE_CHECK="Get canopy array from import state"
            CALL ESMF_StateGet(state    =IMP_STATE                      &
                              ,itemName ='canopy'                       &
                              ,array    =TMP_ARRAY                      &
                              ,rc       =RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)
            MESSAGE_CHECK="Get canopy pointer from imported array"
            CALL ESMF_ArrayGet(TMP_ARRAY,0,gfs_exp_canopy,rc=RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)

            MESSAGE_CHECK="Get hice array from import state"
            CALL ESMF_StateGet(state    =IMP_STATE                      &
                              ,itemName ='hice'                         &
                              ,array    =TMP_ARRAY                      &
                              ,rc       =RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)
            MESSAGE_CHECK="Get hice pointer from imported array"
            CALL ESMF_ArrayGet(TMP_ARRAY,0,gfs_exp_hice,rc=RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)

            MESSAGE_CHECK="Get fice array from import state"
            CALL ESMF_StateGet(state    =IMP_STATE                      &
                              ,itemName ='fice'                         &
                              ,array    =TMP_ARRAY                      &
                              ,rc       =RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)
            MESSAGE_CHECK="Get fice pointer from imported array"
            CALL ESMF_ArrayGet(TMP_ARRAY,0,gfs_exp_fice,rc=RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)

            MESSAGE_CHECK="Get stc array from import state"
            CALL ESMF_StateGet(state    =IMP_STATE                      &
                              ,itemName ='stc'                          &
                              ,array    =TMP_ARRAY                      &
                              ,rc       =RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)
            MESSAGE_CHECK="Get stc pointer from imported array"
            CALL ESMF_ArrayGet(TMP_ARRAY,0,gfs_exp_stc,rc=RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)

            MESSAGE_CHECK="Get smc array from import state"
            CALL ESMF_StateGet(state    =IMP_STATE                      &
                              ,itemName ='smc'                          &
                              ,array    =TMP_ARRAY                      &
                              ,rc       =RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)
            MESSAGE_CHECK="Get smc pointer from imported array"
            CALL ESMF_ArrayGet(TMP_ARRAY,0,gfs_exp_smc,rc=RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)

            MESSAGE_CHECK="Get sfcdsw array from import state"
            CALL ESMF_StateGet(state    =IMP_STATE                      &
                              ,itemName ='sfcdsw'                       &
                              ,array    =TMP_ARRAY                      &
                              ,rc       =RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)
            MESSAGE_CHECK="Get sfcdsw pointer from imported array"
            CALL ESMF_ArrayGet(TMP_ARRAY,0,gfs_exp_sfcdsw,rc=RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)

            MESSAGE_CHECK="Get sfcdlw array from import state"
            CALL ESMF_StateGet(state    =IMP_STATE                      &
                              ,itemName ='sfcdlw'                       &
                              ,array    =TMP_ARRAY                      &
                              ,rc       =RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)
            MESSAGE_CHECK="Get sfcdlw pointer from imported array"
            CALL ESMF_ArrayGet(TMP_ARRAY,0,gfs_exp_sfcdlw,rc=RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)

            MESSAGE_CHECK="Get t2m array from import state"
            CALL ESMF_StateGet(state    =IMP_STATE                      &
                              ,itemName ='t2m'                          &
                              ,array    =TMP_ARRAY                      &
                              ,rc       =RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)
            MESSAGE_CHECK="Get t2m pointer from imported array"
            CALL ESMF_ArrayGet(TMP_ARRAY,0,gfs_exp_t2m,rc=RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)

            MESSAGE_CHECK="Get q2m array from import state"
            CALL ESMF_StateGet(state    =IMP_STATE                      &
                              ,itemName ='q2m'                          &
                              ,array    =TMP_ARRAY                      &
                              ,rc       =RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)
            MESSAGE_CHECK="Get q2m pointer from imported array"
            CALL ESMF_ArrayGet(TMP_ARRAY,0,gfs_exp_q2m,rc=RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)

            MESSAGE_CHECK="Get slmsk array from import state"
            CALL ESMF_StateGet(state    =IMP_STATE                      &
                              ,itemName ='slmsk'                        &
                              ,array    =TMP_ARRAY                      &
                              ,rc       =RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)
            MESSAGE_CHECK="Get slmsk pointer from imported array"
            CALL ESMF_ArrayGet(TMP_ARRAY,0,gfs_exp_slmsk,rc=RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)

            MESSAGE_CHECK="Get hprime array from import state"
            CALL ESMF_StateGet(state    =IMP_STATE                      &
                              ,itemName ='hprime'                          &
                              ,array    =TMP_ARRAY                      &
                              ,rc       =RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)
            MESSAGE_CHECK="Get hprime pointer from imported array"
            CALL ESMF_ArrayGet(TMP_ARRAY,0,gfs_exp_hprime,rc=RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)

            MESSAGE_CHECK="Get fluxr array from import state"
            CALL ESMF_StateGet(state    =IMP_STATE                      &
                              ,itemName ='fluxr'                          &
                              ,array    =TMP_ARRAY                      &
                              ,rc       =RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)
            MESSAGE_CHECK="Get fluxr pointer from imported array"
            CALL ESMF_ArrayGet(TMP_ARRAY,0,gfs_exp_fluxr,rc=RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)


            MESSAGE_CHECK="Get us3d array from export state"
            CALL ESMF_StateGet(state    =EXP_STATE                      &
                              ,itemName ='us3d'                         &
                              ,array    =TMP_ARRAY                      &
                              ,rc       =RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)
            MESSAGE_CHECK="Get us3d pointer from exported array"
            CALL ESMF_ArrayGet(TMP_ARRAY,0,fim_imp_us3d,rc=RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)

            MESSAGE_CHECK="Get vs3d array from export state"
            CALL ESMF_StateGet(state    =EXP_STATE                      &
                              ,itemName ='vs3d'                         &
                              ,array    =TMP_ARRAY                      &
                              ,rc       =RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)
            MESSAGE_CHECK="Get vs3d pointer from exported array"
            CALL ESMF_ArrayGet(TMP_ARRAY,0,fim_imp_vs3d,rc=RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)

            MESSAGE_CHECK="Get tr3d array from export state"
            CALL ESMF_StateGet(state    =EXP_STATE                      &
                              ,itemName ='tr3d'                         &
                              ,array    =TMP_ARRAY                      &
                              ,rc       =RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)
            MESSAGE_CHECK="Get tr3d pointer from exported array"
            CALL ESMF_ArrayGet(TMP_ARRAY,0,fim_imp_tr3d,rc=RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)

            MESSAGE_CHECK="Get rn2d array from export state"
            CALL ESMF_StateGet(state    =EXP_STATE                      &
                              ,itemName ='rn2d'                         &
                              ,array    =TMP_ARRAY                      &
                              ,rc       =RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)
            MESSAGE_CHECK="Get rn2d pointer from exported array"
            CALL ESMF_ArrayGet(TMP_ARRAY,0,fim_imp_rn2d,rc=RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)

            MESSAGE_CHECK="Get rc2d array from export state"
            CALL ESMF_StateGet(state    =EXP_STATE                      &
                              ,itemName ='rc2d'                         &
                              ,array    =TMP_ARRAY                      &
                              ,rc       =RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)
            MESSAGE_CHECK="Get rc2d pointer from exported array"
            CALL ESMF_ArrayGet(TMP_ARRAY,0,fim_imp_rc2d,rc=RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)

            MESSAGE_CHECK="Get ts2d array from export state"
            CALL ESMF_StateGet(state    =EXP_STATE                      &
                              ,itemName ='ts2d'                         &
                              ,array    =TMP_ARRAY                      &
                              ,rc       =RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)
            MESSAGE_CHECK="Get ts2d pointer from exported array"
            CALL ESMF_ArrayGet(TMP_ARRAY,0,fim_imp_ts2d,rc=RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)

            MESSAGE_CHECK="Get us2d array from export state"
            CALL ESMF_StateGet(state    =EXP_STATE                      &
                              ,itemName ='us2d'                         &
                              ,array    =TMP_ARRAY                      &
                              ,rc       =RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)
            MESSAGE_CHECK="Get us2d pointer from exported array"
            CALL ESMF_ArrayGet(TMP_ARRAY,0,fim_imp_us2d,rc=RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)

            MESSAGE_CHECK="Get hf2d array from export state"
            CALL ESMF_StateGet(state    =EXP_STATE                      &
                              ,itemName ='hf2d'                         &
                              ,array    =TMP_ARRAY                      &
                              ,rc       =RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)
            MESSAGE_CHECK="Get hf2d pointer from exported array"
            CALL ESMF_ArrayGet(TMP_ARRAY,0,fim_imp_hf2d,rc=RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)

            MESSAGE_CHECK="Get qf2d array from export state"
            CALL ESMF_StateGet(state    =EXP_STATE                      &
                              ,itemName ='qf2d'                         &
                              ,array    =TMP_ARRAY                      &
                              ,rc       =RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)
            MESSAGE_CHECK="Get qf2d pointer from exported array"
            CALL ESMF_ArrayGet(TMP_ARRAY,0,fim_imp_qf2d,rc=RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)

            MESSAGE_CHECK="Get sheleg2d array from export state"
            CALL ESMF_StateGet(state    =EXP_STATE                      &
                              ,itemName ='sheleg2d'                     &
                              ,array    =TMP_ARRAY                      &
                              ,rc       =RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)
            MESSAGE_CHECK="Get sheleg2d pointer from exported array"
            CALL ESMF_ArrayGet(TMP_ARRAY,0,fim_imp_sheleg2d,rc=RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)

            MESSAGE_CHECK="Get canopy2d array from export state"
            CALL ESMF_StateGet(state    =EXP_STATE                      &
                              ,itemName ='canopy2d'                     &
                              ,array    =TMP_ARRAY                      &
                              ,rc       =RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)
            MESSAGE_CHECK="Get canopy2d pointer from exported array"
            CALL ESMF_ArrayGet(TMP_ARRAY,0,fim_imp_canopy2d,rc=RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)

            MESSAGE_CHECK="Get hice2d array from export state"
            CALL ESMF_StateGet(state    =EXP_STATE                      &
                              ,itemName ='hice2d'                       &
                              ,array    =TMP_ARRAY                      &
                              ,rc       =RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)
            MESSAGE_CHECK="Get hice2d pointer from exported array"
            CALL ESMF_ArrayGet(TMP_ARRAY,0,fim_imp_hice2d,rc=RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)

            MESSAGE_CHECK="Get fice2d array from export state"
            CALL ESMF_StateGet(state    =EXP_STATE                      &
                              ,itemName ='fice2d'                       &
                              ,array    =TMP_ARRAY                      &
                              ,rc       =RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)
            MESSAGE_CHECK="Get fice2d pointer from exported array"
            CALL ESMF_ArrayGet(TMP_ARRAY,0,fim_imp_fice2d,rc=RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)

            MESSAGE_CHECK="Get st3d array from export state"
            CALL ESMF_StateGet(state    =EXP_STATE                      &
                              ,itemName ='st3d'                         &
                              ,array    =TMP_ARRAY                      &
                              ,rc       =RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)
            MESSAGE_CHECK="Get st3d pointer from exported array"
            CALL ESMF_ArrayGet(TMP_ARRAY,0,fim_imp_st3d,rc=RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)

            MESSAGE_CHECK="Get sm3d array from export state"
            CALL ESMF_StateGet(state    =EXP_STATE                      &
                              ,itemName ='sm3d'                         &
                              ,array    =TMP_ARRAY                      &
                              ,rc       =RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)
            MESSAGE_CHECK="Get sm3d pointer from exported array"
            CALL ESMF_ArrayGet(TMP_ARRAY,0,fim_imp_sm3d,rc=RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)

            MESSAGE_CHECK="Get sw2d array from export state"
            CALL ESMF_StateGet(state    =EXP_STATE                      &
                              ,itemName ='sw2d'                         &
                              ,array    =TMP_ARRAY                      &
                              ,rc       =RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)
            MESSAGE_CHECK="Get sw2d pointer from exported array"
            CALL ESMF_ArrayGet(TMP_ARRAY,0,fim_imp_sw2d,rc=RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)

            MESSAGE_CHECK="Get lw2d array from export state"
            CALL ESMF_StateGet(state    =EXP_STATE                      &
                              ,itemName ='lw2d'                         &
                              ,array    =TMP_ARRAY                      &
                              ,rc       =RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)
            MESSAGE_CHECK="Get lw2d pointer from exported array"
            CALL ESMF_ArrayGet(TMP_ARRAY,0,fim_imp_lw2d,rc=RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)

            MESSAGE_CHECK="Get t2m2d array from export state"
            CALL ESMF_StateGet(state    =EXP_STATE                      &
                              ,itemName ='t2m2d'                        &
                              ,array    =TMP_ARRAY                      &
                              ,rc       =RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)
            MESSAGE_CHECK="Get t2m2d pointer from exported array"
            CALL ESMF_ArrayGet(TMP_ARRAY,0,fim_imp_t2m2d,rc=RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)

            MESSAGE_CHECK="Get q2m2d array from export state"
            CALL ESMF_StateGet(state    =EXP_STATE                      &
                              ,itemName ='q2m2d'                        &
                              ,array    =TMP_ARRAY                      &
                              ,rc       =RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)
            MESSAGE_CHECK="Get q2m2d pointer from exported array"
            CALL ESMF_ArrayGet(TMP_ARRAY,0,fim_imp_q2m2d,rc=RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)

            MESSAGE_CHECK="Get slmsk2d array from export state"
            CALL ESMF_StateGet(state    =EXP_STATE                      &
                              ,itemName ='slmsk2d'                      &
                              ,array    =TMP_ARRAY                      &
                              ,rc       =RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)
            MESSAGE_CHECK="Get slmsk2d pointer from exported array"
            CALL ESMF_ArrayGet(TMP_ARRAY,0,fim_imp_slmsk2d,rc=RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)

            MESSAGE_CHECK="Get hprm2d array from export state"
            CALL ESMF_StateGet(state    =EXP_STATE                      &
                              ,itemName ='hprm2d'                       &
                              ,array    =TMP_ARRAY                      &
                              ,rc       =RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)
            MESSAGE_CHECK="Get hprm2d pointer from exported array"
            CALL ESMF_ArrayGet(TMP_ARRAY,0,fim_imp_hprm2d,rc=RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)

            MESSAGE_CHECK="Get flxlwtoa2d array from export state"
            CALL ESMF_StateGet(state    =EXP_STATE                      &
                              ,itemName ='flxlwtoa2d'                   &
                              ,array    =TMP_ARRAY                      &
                              ,rc       =RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)
            MESSAGE_CHECK="Get flxlwtoa2d pointer from exported array"
            CALL ESMF_ArrayGet(TMP_ARRAY,0,fim_imp_flxlwtoa2d,rc=RC)
            CALL ERR_MSG(RC,MESSAGE_CHECK,IRCFINAL)

            ! couple from physics to dynamics 
            CALL CPL_PHY_TO_DYN(its,                      &
                   ! IN args
                   gfs_exp_p, gfs_exp_u  , gfs_exp_v,     &
                   gfs_exp_q, gfs_exp_cld, gfs_exp_t,     &
                   ! these GFS PHY fields are passed to FIM DYN for 
                   ! output and diagnostics only.  
                   gfs_exp_geshem, gfs_exp_rainc,         &
                   gfs_exp_tsea,   gfs_exp_uustar,        &
                   gfs_exp_hflx,   gfs_exp_evap,          &
                   gfs_exp_sheleg, gfs_exp_canopy,        &
                   gfs_exp_hice,   gfs_exp_fice,          &
                   gfs_exp_stc,    gfs_exp_smc,           &
                   gfs_exp_sfcdsw, gfs_exp_sfcdlw,        &
                   gfs_exp_t2m,    gfs_exp_q2m,           &
                   gfs_exp_slmsk,  gfs_exp_hprime,        &
                   gfs_exp_fluxr,                         &
                   ! OUT args
                   fim_imp_us3d, fim_imp_vs3d,            &
                   fim_imp_tr3d,                          &
                   fim_imp_rn2d, fim_imp_rc2d,            &
                   fim_imp_ts2d, fim_imp_us2d,            &
                   fim_imp_hf2d, fim_imp_qf2d,            &
                   fim_imp_sheleg2d, fim_imp_canopy2d,    &
                   fim_imp_hice2d, fim_imp_fice2d,        &
                   fim_imp_st3d, fim_imp_sm3d,            &
                   fim_imp_sw2d, fim_imp_lw2d,            &
                   fim_imp_t2m2d, fim_imp_q2m2d,          &
                   fim_imp_slmsk2d, fim_imp_hprm2d,       &
                   fim_imp_flxlwtoa2d )
          else
            WRITE(0,*)"ERROR:  UNEXPECTED STATE NAME IN CPL_RUN"
            IRCFINAL = esmf_failure
          endif

        endif ! CallPhysics

      endif
!
      IF(IRCFINAL==ESMF_SUCCESS)THEN
!       WRITE(0,*)"CPL RUN STEP SUCCEEDED"
      ELSE
        WRITE(0,*)"CPL RUN STEP FAILED"
      ENDIF
!
      IRC_CPL=IRCFINAL
!
!JR      cpl_dyn_phy_tim=cpl_dyn_phy_tim+timef()-btim0
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE CPL_RUN
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE CPL_FINALIZE(CPL_COMP,IMP_STATE,EXP_STATE,CLOCK,IRC_CPL)
!
!-----------------------------------------------------------------------
!***  FINALIZE THE COUPLER.
!-----------------------------------------------------------------------
!
      ! FIM dynamics + GFS physics
      USE module_fim_cpl_finalize ,only: CPL_FINALIZE_FIM_GFS => cpl_finalize
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!***  ARGUMENT VARIABLES.
!-----------------------------------------------------------------------
!
      TYPE(ESMF_CplComp),INTENT(INOUT) :: CPL_COMP
      TYPE(ESMF_State),  INTENT(INOUT) :: IMP_STATE
      TYPE(ESMF_State),  INTENT(INOUT) :: EXP_STATE
      TYPE(ESMF_Clock),  INTENT(IN)    :: CLOCK
!
      INTEGER,           INTENT(OUT)   :: IRC_CPL
!      
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      INTEGER :: IRCFINAL
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      IRCFINAL=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  FINALIZE THE PHYSICS COMPONENT.
!-----------------------------------------------------------------------
!
      CALL CPL_FINALIZE_FIM_GFS
!
!-----------------------------------------------------------------------
!
      IF(IRCFINAL==ESMF_SUCCESS)THEN
!       WRITE(0,*)"CPL FINALIZE STEP SUCCEEDED"
      ELSE
        WRITE(0,*)"CPL FINALIZE STEP FAILED"
      ENDIF
!
      IRC_CPL=IRCFINAL
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE CPL_FINALIZE
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
      END MODULE MODULE_DYN_PHY_CPL_COMP
!
!-----------------------------------------------------------------------
