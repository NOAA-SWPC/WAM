      MODULE module_FIM_GRID_COMP

      USE ESMF_MOD
      USE MODULE_FIM_INTERNAL_STATE ,ONLY: FIM_INTERNAL_STATE            &
                                          ,WRAP_FIM_INTERNAL_STATE
      USE MODULE_FIM_INTEGRATE      ,ONLY: FIM_INTEGRATE
!TBH:      USE MODULE_ERR_MSG            ,ONLY: ERR_MSG,SET_IPRINT
      USE MODULE_ERR_MSG            ,ONLY: ERR_MSG
      USE MODULE_DYNAMICS_GRID_COMP ,ONLY: DYN_REGISTER
      USE MODULE_PHYSICS_GRID_COMP  ,ONLY: PHY_REGISTER
      USE MODULE_DYN_PHY_CPL_COMP   ,ONLY: DYN_PHY_CPL_REGISTER

      IMPLICIT NONE

      PRIVATE
      PUBLIC :: FIM_REGISTER

      TYPE(FIM_INTERNAL_STATE),POINTER,SAVE :: FIM_INT_STATE
      TYPE(WRAP_FIM_INTERNAL_STATE)   ,SAVE :: WRAP


      CONTAINS

      SUBROUTINE FIM_REGISTER(FIM_GRID_COMP,RC_REG)
      TYPE(ESMF_GridComp),INTENT(INOUT) :: FIM_GRID_COMP
      INTEGER            ,INTENT(OUT)   :: RC_REG

      INTEGER :: RC

!      write(0,*) "    FIM_REGISTER"

!-----------------------------------------------------------------------
!***  Register the fim initialize subroutine.  Since it is just one 
!***  subroutine, use esmf_singlephase.  The second argument is
!***  a pre-defined subroutine type, such as ESMF_SETINIT, ESMF_SETRUN, 
!***  or ESMF_SETFINAL.
!-----------------------------------------------------------------------

      CALL ESMF_GridCompSetEntryPoint(FIM_GRID_COMP ,ESMF_SETINIT ,FIM_INITIALIZE ,ESMF_SINGLEPHASE ,RC)
      CALL ERR_MSG (RC, 'set fim initialize entry point', RC_REG)

!-----------------------------------------------------------------------
!***  Register the Run step of the FIM component.
!-----------------------------------------------------------------------

      CALL ESMF_GridCompSetEntryPoint(FIM_GRID_COMP ,ESMF_SETRUN ,FIM_RUN ,ESMF_SINGLEPHASE ,RC)
      CALL ERR_MSG (RC, 'set fim run entry point', RC_REG)

!-----------------------------------------------------------------------
!***  Register the FIM FINALIZE subroutine.
!-----------------------------------------------------------------------

      CALL ESMF_GridCompSetEntryPoint(FIM_GRID_COMP ,ESMF_SETFINAL ,FIM_FINALIZE ,ESMF_SINGLEPHASE ,RC)
      CALL ERR_MSG (RC, 'set fim finalize entry point', RC_REG)

!      write(0,*) "    END OF FIM_REGISTER"

      END SUBROUTINE FIM_REGISTER


      SUBROUTINE FIM_INITIALIZE(FIM_GRID_COMP ,IMP_STATE ,EXP_STATE ,CLOCK_ATM ,RC_INIT)

      ! TODO:  move logical flags to internal state
      use module_core_setup,only: core_setup_fim,iam_fim_task,iam_write_task
      use module_fim_dyn_init,only:dyn_init

      TYPE(ESMF_GridComp),INTENT(INOUT) :: FIM_GRID_COMP
      TYPE(ESMF_State)   ,INTENT(INOUT) :: IMP_STATE
      TYPE(ESMF_State)   ,INTENT(INOUT) :: EXP_STATE
      TYPE(ESMF_Clock)   ,INTENT(INOUT) :: CLOCK_ATM
      INTEGER            ,INTENT(OUT)   :: RC_INIT

      INTEGER :: RC
      TYPE(ESMF_Config) :: CF

      type(esmf_vm),save :: vm_global, vm_local            ! the esmf virtual machine.
      type(esmf_time) :: currtime                         ! the esmf current time.
      type(esmf_time) :: starttime                        ! the esmf start time.
      type(esmf_timeinterval) :: timestep

      ! note that this is just a reference so no "save" is needed
      type(esmf_grid) :: grid_fim  ! FIM grid, created by DYN, used by PHY and 
                                   ! FIM components

      integer :: total_tasks
      integer :: mype_global
      integer :: timestep_sec_whole
      integer :: timestep_sec_numerator
      integer :: timestep_sec_denominator

      integer :: nfhout, nfmout, nfsout, nsout
      real :: deltim
      integer :: MPI_COMM_FIM   ! MPI communicator for this FIM component
      integer :: ppp__status
      logical :: iprint_lcl
      character(128) :: comp_name


!      write(0,*) "        FIM_INITIALIZE"
      RC_INIT = ESMF_SUCCESS

! Start SMS.  Every build of NEMS.x with FIM must use SMS.  
! Normally this code is created via SMS directive, but PPP is 
! not run on this file.  
      CALL sms_start(ppp__status)
      IF (ppp__status.NE.0) THEN
        ! Follow NMM template from Tom Black for non-ESMF error-abort
        write(0,*) "ERROR IN FIM_INITIALIZE:  sms_start FAILED"
        CALL ESMF_Finalize(RC=RC,terminationflag=ESMF_ABORT)
      ENDIF
!
!-----------------------------------------------------------------------
!***  Allocate the FIM component's internal state, point at it,
!***  and attach it to the FIM component.
!-----------------------------------------------------------------------
!
      ALLOCATE(FIM_INT_STATE,stat=RC)
      WRAP%FIM_INT_STATE=>FIM_INT_STATE

!JR This call stuffs "wrap" into "fim_grid_comp" for later retrieval
      CALL ESMF_GridCompSetInternalState(FIM_GRID_COMP ,WRAP ,RC)
      CALL ERR_MSG (RC, 'ESMF_GridCompSetInternalState', RC_INIT)

!
!-----------------------------------------------------------------------
!***  For the moment, use a direct copy of the ATM Clock.
!-----------------------------------------------------------------------
!
!JR This should probably use ESMF's copy constructor method instead
      fim_int_state%clock_fim = clock_atm

!TBH:  This bit is useful for debugging
      call esmf_gridcompget (gridcomp=fim_grid_comp, name=comp_name, rc=rc)
      call err_msg (rc, 'get name of fim_grid_comp', rc_init)
!JR print *,'DEBUG:  name of fim_grid_comp is [',TRIM(comp_name),']'

!-----------------------------------------------------------------------
!***  Attach the configure file to the FIM component.
!-----------------------------------------------------------------------
      CF=ESMF_ConfigCreate(rc=RC)
!TODO:  hard-coded 'fim.configure' is not DRY, work with NCEP to fix this
      CALL ESMF_ConfigLoadFile(config=CF ,filename='fim.configure' ,rc=RC)
      CALL ERR_MSG (rc, 'load configure file fim.configure into configure object', rc_init)
      CALL ESMF_GridCompSet(gridcomp=FIM_GRID_COMP, config=CF, rc=RC)
      CALL ERR_MSG (rc, 'attach configure object to fim component', rc_init)

!
!-----------------------------------------------------------------------
!***  Set verbostiy of err_msg prints
!TODO:  move this setting up into MAIN_NEMS.F90, work with NCEP
!TBH:  Temporarily removed this until we can merge set_iprint() into 
!TBH:  nems repository.  
!-----------------------------------------------------------------------
!
!TBH    call esmf_configgetattribute (config=cf, value=iprint_lcl, label='iprint:', rc=rc)
!TBH    call err_msg (rc, "extract iprint information from fim config file", rc_init)
!TBH    call set_iprint(iprint_lcl)
!
!-----------------------------------------------------------------------
!***  Retrieve the VM from the FIM component.
!-----------------------------------------------------------------------
!
    call esmf_gridcompget (gridcomp=fim_grid_comp, vm=vm_local, rc=rc)
    call err_msg (rc, "retrieve the cf and vm from fim component", rc_init)
!-----------------------------------------------------------------------
!***  Retrieve global VM then the total number of tasks for
!***  then entire system.
!-----------------------------------------------------------------------
    call esmf_vmgetglobal (vm=vm_global, rc=rc)
    call err_msg (rc, "retrieve global vm_global for fim", rc_init)

    call esmf_vmget (vm=vm_global, pecount=total_tasks, localpet=mype_global, rc=rc)
    call err_msg (rc, "fim_initialize: obtain global mpi task id from vm_global", rc_init)
!JR print *,'DEBUG:  mype_global = ',mype_global
    call esmf_vmget (vm=vm_local, localpet=fim_int_state%mype, rc=rc)
    call err_msg (rc, "fim_initialize: obtain local mpi task id from vm_local", rc_init)
!JR print *,'DEBUG:  fim_int_state%mype = ',fim_int_state%mype
!-----------------------------------------------------------------------
!***  Extract fundamental timestep information from the config file.
!-----------------------------------------------------------------------
    call esmf_configgetattribute (config=cf, value=timestep_sec_whole, label ='dt_int:', rc=rc)
    call esmf_configgetattribute (config=cf, value=timestep_sec_numerator, label ='dt_num:', rc=rc)
    call esmf_configgetattribute (config=cf, value =timestep_sec_denominator, label ='dt_den:', rc=rc)
    call err_msg (rc, "extract timestep information from fim config file", rc_init)
!-----------------------------------------------------------------------
!***  Establish the timestep for the FIM Clock.
!-----------------------------------------------------------------------
    call esmf_timeintervalset (timeinterval=timestep, s=timestep_sec_whole, sn=timestep_sec_numerator, &
                               sd=timestep_sec_denominator, rc=rc)
    call esmf_clockset (clock=fim_int_state%clock_fim, timestep = timestep, rc=rc)
    call err_msg (rc, "set time step interval in fim clock", rc_init)

!TBH:    Note that DYN, PHY, and CPL must currently live on the same 
!TBH:    MPI tasks via the use of "petlist=fim_int_state%petlist_fcst" 
!TBH:    during ESMF_*CompCreate() calls.   

!
!-----------------------------------------------------------------------
!***  SEGREGATE THE FORECAST TASKS FROM THE QUILT/WRITE TASKS.
!***  VIA CALL TO core_setup_fim WHICH RETURNS PETLISTS.  
!-----------------------------------------------------------------------
!
    CALL ESMF_VMGet(vm=vm_local,mpiCommunicator=MPI_COMM_FIM,rc=RC)
    call err_msg (rc, "extract mpiCommunicator from vm_local", rc_init)

#ifdef MANUALGPTL
    ret = gptlstart ('core_setup_fim')
#endif

! Split VM between compute and write tasks via petlists.  
! core_setup_fim() allocates and initializes the petlists.  
! NOTE:  Executable SMS directives must not be placed before 
! NOTE:  this call!  
    CALL core_setup_fim(MPI_COMM_FIM,fim_int_state%petlist_fcst, &
                                     fim_int_state%petlist_write)
#ifdef MANUALGPTL
    ret = gptlstop ('core_setup_fim')
#endif

!-----------------------------------------------------------------------
!***  Will the Write components with asynchronous quilting be used?
!-----------------------------------------------------------------------
!TBH    call esmf_configgetattribute (config=cf,                     &  ! the fim config object
!TBH                                  value =fim_int_state%quilting, &  ! the quilting flag
!TBH                                  label ='quilting:',            &  ! give label's value to the previous variable
!TBH                                  rc    =rc)
!TBH    call err_msg (rc, "extract quilting flag from fim config file", rc_init)

! Note:  At the moment, the 'quilting' field in the config file is ignored in 
! Note:  favor of FIMnamelist settings read by core_setup_fim.  
!TODO:  Reconcile 'quilting' field and FIM write task settings.  
    if (size(fim_int_state%petlist_write) > 0) then
      fim_int_state%quilting = .true.
    else
      fim_int_state%quilting = .false.
    endif

!TODO:  Connect fim_int_state%quilting to creation of NEMS write components.  
!TODO:  FIM write tasks are not yet ESMF components.  Per UMIG discussion 
!TODO:  with Mark Iredell during UMIG meeting on 1/14/2011, it is more 
!TODO:  important to get write tasks working with NEMSIO before encapsulating 
!TODO:  them as components.  

!-----------------------------------------------------------------------
!***  Establish the frequency of forecast output.
!TODO:  Reconcile FIMnamelist and config file here
!-----------------------------------------------------------------------
    call esmf_configgetattribute (config=cf, value=nfhout, label ='nfhout:', rc=rc)
    call err_msg (rc, "extract nfhout from fim config file", rc_init)
    call esmf_configgetattribute (config=cf, value=nsout, label ='nsout:', rc=rc)
    call err_msg (rc, "extract nsout from fim config file", rc_init)
    call esmf_configgetattribute (config=cf, value =deltim, label ='deltim:', rc=rc)
    call err_msg (rc, "extract history output interval from fim config file", rc_init)
    call esmf_timeintervalset (fim_int_state%timeinterval_fim_output, h=nfhout, m=nfmout, &
                               s=nfsout, rc=rc)
    call err_msg (rc, "set fim history output interval", rc_init)
!-----------------------------------------------------------------------
!***  Extract the start time from the clock.
!-----------------------------------------------------------------------
    call esmf_clockget (clock=fim_int_state%clock_fim, starttime=starttime, rc=rc)
    call err_msg (rc, "fim_atm_init: start time from fim clock", rc_init)
    currtime = starttime
!-----------------------------------------------------------------------
!***  No need to extract the RESTART flag from the configure file: FIM currently has no restart capability
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  Create the Dynamics gridded subcomponent.
! all tasks must participate in component creation even if they do 
! not participate in a component
!-----------------------------------------------------------------------
    fim_int_state%gc_fim_dyn = esmf_gridcompcreate (name="dynamics component", &
                                                    configfile='fim.configure', &
                                                    petlist=fim_int_state%petlist_fcst, &
                                                    rc=rc)
    call err_msg (rc, "create the fim dynamics component", rc_init)

!-----------------------------------------------------------------------
!***  Always create the Physics gridded subcomponent.
!-----------------------------------------------------------------------
    fim_int_state%gc_fim_phy = esmf_gridcompcreate (name      ="physics component", &
                                                    configfile='fim.configure', &
                                                    petlist   =fim_int_state%petlist_fcst, &
                                                    rc        =rc)
    call err_msg (rc, "create the fim physics component", rc_init)
!-----------------------------------------------------------------------
!***  Create the Dynamics-Physics coupler subcomponent.
!-----------------------------------------------------------------------
    fim_int_state%gc_fim_cpl = esmf_cplcompcreate (name="coupler component",           &
                                                   petlist=fim_int_state%petlist_fcst, &
                                                   rc     =rc)
    call err_msg (rc, "create the fim dynamics-physics coupler component", rc_init)

! all tasks must participate in setservices calls even if they do 
! not participate in a component (ESMF refdoc 15.3.8)
!***  Register the Initialize, Run, and Finalize steps for created components
    call esmf_gridcompsetservices (fim_int_state%gc_fim_dyn, dyn_register, rc)
    call err_msg (rc, "register fim dynamics init, run, finalize", rc_init)
    call esmf_gridcompsetservices (fim_int_state%gc_fim_phy, phy_register, rc)
    call err_msg (rc, "register physics init, run, finalize", rc_init)
    call esmf_cplcompsetservices (fim_int_state%gc_fim_cpl, dyn_phy_cpl_register, rc)
    call err_msg (rc, "register the dyn-phy coupler's init, run, finalize", rc_init)

!-----------------------------------------------------------------------
!***  Create empty Import and Export states for the Dynamics component
!-----------------------------------------------------------------------
    fim_int_state%IMP_fim_DYN = ESMF_StateCreate (statename="FIM dynamics import",      &
                                                  statetype=esmf_state_import,          &
                                                  rc       =RC)
    call err_msg (rc, "create empty import state for fim dynamics", rc_init)
    fim_int_state%exp_fim_DYN = ESMF_StateCreate (statename="FIM dynamics export",      &
                                                  statetype=esmf_state_export,          &
                                                  rc       =RC)
    call err_msg (rc, "create empty export state for fim dynamics", rc_init)
!----------------------------------------------------------------------------------
! Add the FIM dynamics states as the nested states into the FIM parent states.
!----------------------------------------------------------------------------------
    call esmf_stateadd (imp_state, fim_int_state%imp_fim_dyn, rc)
    call esmf_stateadd (exp_state, fim_int_state%exp_fim_dyn, rc)
!JR I don't know what cpl_flag does, so copy from GFS for now
    fim_int_state%cpl_flag = esmf_false
    call esmf_attributeset (fim_int_state%imp_fim_dyn, 'cpl_flag', fim_int_state%cpl_flag, rc)
    call err_msg (rc, "fim set cpl_flag", rc_init)
!------------------------------------------------------------------------
!***  Create empty Import and Export states for the Physics subcomponent.
! Note that statenames are used by CPL_RUN() to determine direction of 
! coupling (DYN->PHY or PHY->DYN).  Any changes must be matched by 
! equivalent changes in CPL_RUN().  
!------------------------------------------------------------------------
    fim_int_state%imp_fim_phy = esmf_statecreate (statename="FIM physics import", &
                                                  statetype=esmf_state_import,    &
                                                  rc       =rc)
    call err_msg (rc, "create empty import state for fim physics", rc_init)
    fim_int_state%exp_fim_phy = esmf_statecreate (statename="FIM physics export", &
                                                  statetype=esmf_state_export,    &
                                                  rc       =rc)
    call err_msg (rc, "create empty export state for fim physics", rc_init)

!-----------------------------------------------------------------------
!***  Setup the Write component(s) (which may run without quilting).
!JR Keep this: FIM will require mods vs GFS
!-----------------------------------------------------------------------
!
!JR    write(0,*)'before write_setup_fim, allocate,write_groups=', fim_int_state%write_groups
!    allocate (fim_int_state%wrt_comps(fim_int_state%write_groups))
!JR Comment out for now
!JR    call write_setup_fim (fim_grid_comp, fim_int_state%wrt_comps, fim_int_state%exp_fim_dyn, &
!JR                          fim_int_state%exp_fim_phy, fim_int_state%imp_fim_wrt, fim_int_state%exp_fim_wrt)

!-----------------------------------------------------------------------
!***  Execute the Initialize steps for the gridded subcomponents.
!***  These are the Initialize subroutines specified in the
!***  Register routines called in ESMF_GridCompSetServices above.
! Note:  Only tasks listed in the petlist during component creation 
! Note:  participate in the esmf_gridcompinitialize (etc.) calls.  
! Note:  For example, the write tasks return immediately from the calls 
! Note:  to esmf_gridcompinitialize for fim_int_state%gc_fim_??? below.  
!-----------------------------------------------------------------------

! TODO:  Make sure that parent does *not* set an ESMF_Grid in the FIM component!  
! TODO:  Per Gerhard we will need to call esmf_gridcompget() to grab the attached 
! TODO:  ESMF_Grid and then call ESMF_GridValidate() on it.  Unfortunately at the 
! TODO:  moment, ESMF regards calling of ESMF_GridValidate() on an un-initialized 
! TODO:  ESMF_Grid as an error!  Straighten this out later.  

! TODO:  Make DYN set lat & lon values in the ESMF_Grid and make PHY extract 
! TODO:  them, then eliminate passing of lat/lon arrays during CPL INIT.  
! TODO:  Ultimately we'll want to modify the NEMS GFS PHY component to accept 
! TODO:  and ESMF_Grid from its parent, at least optionally, and extract lat/lon 
! TODO:  from it.  Tom Black agrees with this approach.  
! TODO:  See module_DYN_PHY_CPL_COMP.F90::CPL_INITIALIZE() for a long-winded 
! TODO:  explanation.  

!--------------
!***  DYNAMICS
!--------------
    ! DYN initialize creates the FIM ESMF_Grid and attaches it to 
    ! fim_int_state%gc_fim_dyn.  This delegation is appropriate since 
    ! the design and implementation of DYN is intimately related to 
    ! the grid.  The PHY and FIM components just use this grid.  Note 
    ! that because of this dependence, DYN initialize must be called 
    ! before PHY initialize.  
    ! Note also that DYN is responsible for destroying the grid.  
!TBH:  This is what the ESMF_GridCompInitialize() call does, but do *not* 
!TBH:  call it directly because ESMF_GridCompInitialize() has needed 
!TBH:  side-effects.  
!TBH    call dyn_initialize (fim_int_state%gc_fim_dyn, &
!TBH                         fim_int_state%imp_fim_dyn, &
!TBH                         fim_int_state%exp_fim_dyn, &
!TBH                         fim_int_state%clock_fim, &
!TBH                         rc)
    call esmf_gridcompinitialize (fim_int_state%gc_fim_dyn,               &
                                  importstate=fim_int_state%imp_fim_dyn,  &
                                  exportstate=fim_int_state%exp_fim_dyn,  &
                                  clock      =fim_int_state%clock_fim,    &
                                  phase      =esmf_singlephase,           &
                                  rc         =rc)
    call err_msg (rc, "initialize fim dynamics component", rc_init)

    ! extract ESMF_Grid from fim_int_state%gc_fim_dyn and set in 
    ! fim_int_state%gc_fim_phy and fim_grid_comp
    ! it is only safe to do this on the FIM compute tasks
    if (iam_fim_task) then
      ! grab FIM grid from DYN
      call esmf_gridcompget(fim_int_state%gc_fim_dyn, grid=grid_fim, rc=rc)
      CALL err_msg(rc, "extract FIM grid from DYN component", rc_init)
      ! is this a good grid?  
      call esmf_gridvalidate(grid=grid_fim, rc=rc)
      CALL err_msg(RC,"validate FIM grid",rc_init)
      ! attach grid to FIM component
      call esmf_gridcompset(fim_grid_comp, grid=grid_fim, rc=rc)
      call err_msg (rc, "attach FIM grid to FIM component", rc_init)
      ! attach grid to PHY component
      call esmf_gridcompset(fim_int_state%gc_fim_phy, grid=grid_fim, rc=rc)
      call err_msg (rc, "attach grid to fim physics component", rc_init)
    endif

!-------------
!***  PHYSICS
!-------------
!TBH:  This is what the ESMF_GridCompInitialize() call does, but do *not* 
!TBH:  call it directly because ESMF_GridCompInitialize() has needed 
!TBH:  side-effects.  
!TBH      call phy_initialize (fim_int_state%gc_fim_phy, fim_int_state%imp_fim_phy, &
!TBH                           fim_int_state%exp_fim_phy, fim_int_state%clock_fim, rc)
    call esmf_gridcompinitialize (fim_int_state%gc_fim_phy,               &
                                  importstate=fim_int_state%imp_fim_phy,  &
                                  exportstate=fim_int_state%exp_fim_phy,  &
                                  clock      =fim_int_state%clock_fim,    &
                                  phase      =esmf_singlephase,           &
                                  rc         =rc)
    call err_msg (rc, "initialize fim physics component", rc_init)

!--------------
!***  DYN-PHY COUPLER COMPONENT
!--------------
!TBH:  This is what the ESMF_CplCompInitialize() call does, but do *not* 
!TBH:  call it directly because ESMF_CplCompInitialize() has needed 
!TBH:  side-effects.  
!TBH    call cpl_initialize (fim_int_state%gc_fim_cpl, fim_int_state%exp_fim_dyn, &
!TBH                         fim_int_state%imp_fim_phy, fim_int_state%clock_fim, rc)
    call esmf_cplcompinitialize (cplcomp    =fim_int_state%gc_fim_cpl,   &
                                 importstate=fim_int_state%exp_fim_dyn,  &
                                 exportstate=fim_int_state%imp_fim_phy,  &
                                 clock      =fim_int_state%clock_fim,    &
                                 rc         =rc)
    call err_msg (rc, "initialize dyn-phy coupler", rc_init)

    !TBH:  Recall that the above call to 
    !TBH:  esmf_gridcompinitialize (fim_int_state%gc_fim_dyn ... ) is a NOOP 
    !TBH:  on write tasks because gc_fim_dyn was created with a "petlist" 
    !TBH:  optional argument that excluded the write tasks from participation 
    !TBH:  in this DYN components method calls.  See for example the above 
    !TBH:  call to:  
    !TBH:  fim_int_state%gc_fim_dyn = esmf_gridcompcreate (...)
    ! TODO:  move this into FIM write component
    if (iam_write_task) then
      call dyn_init(.false.)
    endif

!      write(0,*) "        END OF FIM_INITIALIZE"
      END SUBROUTINE FIM_INITIALIZE


      SUBROUTINE FIM_RUN(FIM_GRID_COMP ,IMP_STATE ,EXP_STATE ,CLOCK_ATM ,RC_RUN)

      ! TODO:  move these to internal state
      use module_core_setup,only: iam_fim_task,iam_write_task
      ! TODO:  move this into FIM write component
      use icosio,only:icosio_run

      TYPE(ESMF_GridComp),INTENT(INOUT) :: FIM_GRID_COMP
      TYPE(ESMF_State)   ,INTENT(INOUT) :: IMP_STATE
      TYPE(ESMF_State)   ,INTENT(INOUT) :: EXP_STATE
      TYPE(ESMF_Clock)   ,INTENT(INOUT) :: CLOCK_ATM
      INTEGER            ,INTENT(OUT)   :: RC_RUN

!---------------------
!***  Local variables
!---------------------
!
      type(esmf_config) :: cf
      integer :: rc                             ! error signal variables.
      type(esmf_timeinterval) :: runduration    ! the forecast length
      type(esmf_timeinterval) :: timestep       ! the fundamental timestep
      type(esmf_time) :: currtime               ! the esmf current time.
      type(esmf_time) :: starttime              ! the esmf start time.

!      write(0,*) "        FIM_RUN"
      RC_RUN=ESMF_SUCCESS

      if (iam_fim_task) then
        ! compute tasks execute this branch

!-----------------------------------------------------------------------
!***  For the moment, use a direct copy of the ATM Clock.
!JR This isn't good
!-----------------------------------------------------------------------
!
      fim_int_state%clock_fim = clock_atm
!-----------------------------------------------------------------------
!***  Extract the fundamental time information from the Clock.
!-----------------------------------------------------------------------
      call esmf_clockget (clock       =fim_int_state%clock_fim,     &  ! the esmf clock
                          timestep    =timestep,                     &  ! the model's timestep length
                          starttime   =starttime,                    &  ! the forecast start time
                          currtime    =currtime,                     &  ! the clock's current time
                          runduration =runduration,                  &  ! the length of the forecast
                          rc          =rc)
      call err_msg (rc, "retrieve fim timestep from the atm clock", rc_run)
!-----------------------------------------------------------------------
!***  Extract the configure file from the ATM component.
!***  GFS needed DFI info here but FIM doesn't
!-----------------------------------------------------------------------
      call esmf_gridcompget (gridcomp=fim_grid_comp, config=cf, rc=rc)
!-----------------------------------------------------------------------
!***  Execute the FIM forecast runstream.
!-----------------------------------------------------------------------
      call fim_integrate (fim_int_state%gc_fim_dyn,                  &
                          fim_int_state%gc_fim_phy,                  &
                          fim_int_state%gc_fim_cpl,                  &
                          fim_int_state%imp_fim_dyn,                 &
                          fim_int_state%exp_fim_dyn,                 &
                          fim_int_state%imp_fim_phy,                 &
                          fim_int_state%exp_fim_phy,                 &
                          fim_int_state%clock_fim,                   &
                          rc_run)

      else if (iam_write_task) then
        ! TODO:  move this into FIM write component
        ! write tasks execute this branch, if present
!
!-----------------------------------------------------------------------
!***  Call the run method for the optional write tasks.
!-----------------------------------------------------------------------
!
        call icosio_run

      endif

      ! do-nothing tasks just print along with the others
!      write(0,*) "        END OF FIM_RUN"
      END SUBROUTINE FIM_RUN


      SUBROUTINE FIM_FINALIZE(FIM_GRID_COMP ,IMP_STATE ,EXP_STATE ,CLOCK_ATM ,RC_FINALIZE)

      ! TODO:  move to internal state
      use module_core_setup,only: iam_fim_task

      TYPE(ESMF_GridComp),INTENT(INOUT) :: FIM_GRID_COMP
      TYPE(ESMF_State)   ,INTENT(INOUT) :: IMP_STATE
      TYPE(ESMF_State)   ,INTENT(INOUT) :: EXP_STATE
      TYPE(ESMF_Clock)   ,INTENT(INOUT) :: CLOCK_ATM
      INTEGER            ,INTENT(OUT)   :: RC_FINALIZE

!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      type(esmf_config)                 :: cf                              ! config object
      integer :: i,j
      integer :: rc,rc_final                                               ! the final error signal variables.
      TYPE(ESMF_VM) :: VM

!      write(0,*) "        FIM_FINALIZE"
      RC_FINALIZE=ESMF_SUCCESS

    rc          = esmf_success
    rc_final    = esmf_success

!-----------------------------------------------------------------------
!***  For the moment, use a direct copy of the ATM Clock.
!JR Fix this
!-----------------------------------------------------------------------
!
    fim_int_state%clock_fim = clock_atm

    call esmf_gridcompget (gridcomp=fim_grid_comp, config=cf, rc=rc)
    call err_msg (rc, "Retrieve Config Object from FIM Component", rc_final)

!-----------------------------------------------------------------------
!***  Finalize each of the subcomponents.
!-----------------------------------------------------------------------
!
if (iam_fim_task) then
!-----------------------------
!***  DYNAMICS-PHYSICS COUPLER
!-----------------------------
    call esmf_cplcompfinalize (fim_int_state%gc_fim_cpl,               &
                               importstate=fim_int_state%exp_fim_dyn,  &
                               exportstate=fim_int_state%imp_fim_phy,  &
                               clock      =fim_int_state%clock_fim,    &
                               rc         =rc)
    call err_msg (rc, "finalize dynamics-physics coupler", rc_final)
!--------------
!***  Dynamics
!--------------
    call esmf_gridcompfinalize (fim_int_state%gc_fim_dyn,             &
                                importstate=fim_int_state%imp_fim_dyn, &
                                exportstate=fim_int_state%exp_fim_dyn, &
                                clock      =fim_int_state%clock_fim,   &
                                rc         =rc)
    call err_msg (rc, "finalize dynamics component", rc_final)
!--------------
!***  Physics
!--------------
    call esmf_gridcompfinalize (fim_int_state%gc_fim_phy,              &
                                importstate=fim_int_state%imp_fim_phy, &
                                exportstate=fim_int_state%exp_fim_phy, &
                                clock      =fim_int_state%clock_fim,   &
                                rc         =rc)
    call err_msg (rc, "finalize physics component", rc_final)
endif

    ! ensure that "write" and "do-nothing" tasks do not race ahead and 
    ! destroy components before other tasks are ready
    CALL ESMF_VMGetCurrent(vm=VM,rc=RC)
    call err_msg (rc, "FIM_FINALIZE:  get current vm", rc_final)
    CALL ESMF_VMBarrier(vm=VM,rc=RC)
    call err_msg (rc, "FIM_FINALIZE:  barrier on current vm", rc_final)

!-----------------------------------------------------------------------
!***  DESTROY ALL STATES.
!-----------------------------------------------------------------------
    call esmf_statedestroy (fim_int_state%imp_fim_dyn, rc=rc)
    call esmf_statedestroy (fim_int_state%exp_fim_dyn, rc=rc)
    call esmf_statedestroy (state=fim_int_state%imp_fim_phy, rc=rc)
    call esmf_statedestroy (state=fim_int_state%exp_fim_phy, rc=rc)
!-----------------------------------------------------------------------
!***  IF QUILTING WAS SELECTED FOR THE GENERATION OF OUTPUT,
!***  FINALIZE AND DESTROY OBJECTS RELATED TO THE WRITE COMPONENTS.
!-----------------------------------------------------------------------
    if (fim_int_state%quilting) then
!JR turn this off until quilting is implemented
!JR      call write_destroy_fim (fim_grid_comp, fim_int_state%wrt_comps, fim_int_state%imp_fim_wrt, &
!JR                              fim_int_state%exp_fim_wrt, fim_int_state%clock_fim) 
    end if

!-----------------------------------------------------------------------
!***  DESTROY ALL SUBCOMPONENTS.
!-----------------------------------------------------------------------
    call esmf_gridcompdestroy (fim_int_state%gc_fim_dyn, rc=rc)
    call err_msg (rc, "Destroy Dynamics Component", rc_final)
!-------------
!***  PHYSICS
!-------------
    call esmf_gridcompdestroy (fim_int_state%gc_fim_phy, rc=rc)
    call err_msg (rc, "destroy physics component", rc_final)
!------------------------------
!***  DYNAMICS-PHYSICS COUPLER
!------------------------------
    call esmf_cplcompdestroy (fim_int_state%gc_fim_cpl, rc=rc)
    call err_msg (rc, "destroy dynamics-physics coupler", rc_final)

!TODO:  make sure all allocated bits of fim_int_state are deallocated
! deallocate other components of fim_int_state
    deallocate(fim_int_state%petlist_fcst)
    deallocate(fim_int_state%petlist_write)

      rc_finalize = rc_final
!      write(0,*) "        END OF FIM_FINALIZE"
      END SUBROUTINE FIM_FINALIZE

      END MODULE module_FIM_GRID_COMP
