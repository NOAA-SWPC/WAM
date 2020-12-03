!#include "../../../ESMFVersionDefine.h"
#include "wam_defs.h"

!! !module: gfs_physics_grid_comp_mod --- 
!                       esmf gridded component of gfs physics
!
! !description: gfs physics gridded component main module.
!
! !revision history:
!
!  july     2007     shrinivas moorthi
!  november 2007     hann-ming henry juang 
!  may      2009     jun wang, add write grid component
!  october  2009     jun wang, output every time step, add no quilting option 
!  oct 09   2009     sarah lu, 3D Gaussian grid (DistGrid5) added
!  oct 12   2009     sarah lu, set the association between imp/exp states
!                    and internal state grid_fld; reset start_step
!  oct 17 2009      Sarah Lu, add debug print to check imp/exp state
!  dec 10 2009      Sarah Lu, add debug print to chekc fcld
!  dec 15 2009      Sarah Lu, add debug print to chekc 3d diag fld (fcld, dqdt)
!  Feb 05 2010      Jun Wang, set init time for restart
!  Mar 02 2010      Sarah Lu, associate export state with internal state in init
!  Apr 11 2010      Sarah Lu, debug print removed
!  Aug 25 2010      Jun Wang, output half dfi filted fields
!  Oct 16 2010      Sarah Lu, retrieve fscav from exp state
!  Dec 23 2010      Sarah Lu, modify fscav initialization 
!  Mar 30 2011      Weiyu Yang, modified code to avoid ESMF log error.
!  May 05 2011      Weiyu Yang, modified for using the ESMF 5.2.0r_beta_snapshot_07.
!  Nov 27 2011      Sarah Lu, add kdt to physics export state
!  Feb 06 2012      Weiyu Yang, modified for using the ESMF 5.2.0r library.
!  Mar 31 2014      S Moorthi Add code to extract ocean surface temperature exported
!                             by the mediator - code provided by Gerhard Theurich
!                             Code under NUOPC ifdef
!  Sep 29 2014      Sarah Lu, remove the code to retrieve the scavenging coefficients 
!                             from physics export state
!  Oct 20 2015      Weiyu Yang, add inputted f10.7 and kp data for the WAM model.
!  Jan    2016      P. Tripp, Coupling use importFields for NUOPC/GSM merge
!                           
!
! !interface:
!
      module gfs_physics_grid_comp_mod
 
!!uses:
!------
      use ESMF
      use NUOPC
      use module_CPLFIELDS, only: importFields, NImportFields, &
                                global_lats_ptr, lonsperlat_ptr
! modules for internal coupling
      use resol_def,                      only: lonr, latr, kdt_start
      use layout1,                        only: me, lats_node_r, lats_node_r_max
      USE ESMF

      use gfs_physics_err_msg_mod,        ONLY: gfs_physics_err_msg,        &
                                                gfs_physics_err_msg_final
      use gfs_physics_initialize_mod,     ONLY: gfs_physics_initialize
      use gfs_physics_run_mod,            ONLY: gfs_physics_run
      use gfs_physics_finalize_mod,       ONLY: gfs_physics_finalize
      USE gfs_physics_getcf_mod,          ONLY: gfs_physics_getcf
      USE gfs_physics_internal_state_mod, ONLY: gfs_physics_internal_state, &
                                                gfs_phy_wrap
      USE mpi_def,                        ONLY: mpi_comm_all,quilting
      USE layout1,                        ONLY: me
      USE date_def,                       ONLY: idate, fhour
      USE namelist_physics_def,           ONLY: fhini, fhmax, lssav,ndfi,ldfi
!jw
      USE gfs_physics_output,             ONLY: point_physics_output_gfs
!
      use GFS_Phy_States_Mod,             ONLY: gfs_physics_import2internal, &  
                                                gfs_physics_internal2export  
      use gfs_phy_tracer_config,          ONLY: gfs_phy_tracer
      use wam_ifp_mod,                    ONLY: farr, params, read_ifp
!
      implicit none

      integer                            :: timestep_sec

      private   ! by default, data is private to this module

      public gfs_phy_setservices      ! only set service is public

!eop
!-------------------------------------------------------------------------


      contains


!----------------------------------------------------------------------
!bop
!
! !routine: gfs_phy_setservices --- 
!           set services for gfs physics gridded component.
! 
! !interface:
!
      subroutine gfs_phy_setservices (gc_gfs_phy, rc)
 
! !arguments:
!------------

      type(esmf_gridcomp)              :: gc_gfs_phy 	! gridded component
      integer,             intent(out) :: rc    	! return code
     
! !description: set services (register) for the gfs physics grid component.
!         
!eop         
!----------------------------------------------------------------------
  
      integer                            :: rc1     = esmf_success

! initializing the error signal variable rc.
!-------------------------------------------
      rc = esmf_success

! register services for this component
! ------------------------------------

! register the initialize subroutine.  since it is just one subroutine
! for the initialize, no need to specify phase.  the second argument is
! a pre-defined subroutine type, such as esmf_setinit, esmf_setrun, 
! esmf_setfinal.
!---------------------------------------------------------------------
      call esmf_logwrite("set entry point for initialize",              &
                          ESMF_LOGMSG_INFO, rc = rc1)
      call esmf_gridcompsetentrypoint (gc_gfs_phy,                      &
                                      ESMF_METHOD_INITIALIZE,           &
                                      gfs_phy_initialize,               &
                                      rc=rc1)
      call gfs_physics_err_msg(rc1,'set entry point for initialize',rc)

! register the run subroutine.
!-----------------------------
      call esmf_logwrite("set entry point for run", ESMF_LOGMSG_INFO, rc = rc1)
      call esmf_gridcompsetentrypoint (gc_gfs_phy,                      &
                                       ESMF_METHOD_RUN,                 &
                                       gfs_phy_run,                     &
                                       rc=rc1)
      call gfs_physics_err_msg(rc1,'set entry point for run',rc)


! register the finalize subroutine.
!----------------------------------
      call esmf_logwrite("set entry point for finalize",                &
                          ESMF_LOGMSG_INFO, rc = rc1)
      call esmf_gridcompsetentrypoint (gc_gfs_phy,                      &
                                       ESMF_METHOD_FINALIZE,            &
                                       gfs_phy_finalize,                &
                                       rc=rc1)
      call gfs_physics_err_msg(rc1,'set entry point for finalize',rc)

! check the error signal variable and print out the result.
!----------------------------------------------------------
      call gfs_physics_err_msg_final(rc1,                               &
                        'setservice for gfs physics grid comp.',rc)

      end subroutine gfs_phy_setservices





!----------------------------------------------------------------------
!bop
! !routine:  gfs_phy_initialize --- initialize routine to initialize 
!                                   and set up the gfs running job.
!
! !description: this subroutine initializes the gfs running before
!               the main running loop.
!
!
! !revision history:
!
!  november 2004     weiyu yang initial code.
!  may      2005     weiyu yang for the updated gfs version.
!  february 2006     moorthi
!  february 2007     h.-m. h. juang
!  may      2009     j. wang
!
! !interface:
!

! this argument list is a standard list for all the initialize,
! the run and finalize routines for an esmf system.
!--------------------------------------------------------------
      subroutine gfs_phy_initialize(gc_gfs_phy,                         &
                                    imp_gfs_phy, exp_gfs_phy, clock, rc)

! user code, for computations related to the esmf interface states.
!------------------------------------------------------------------
!*    use gfs_physics_states_mod
      use gfs_physics_grid_create_mod
!
! !input/output variables and parameters:
!----------------------------------------

      type(esmf_gridcomp)                :: gc_gfs_phy 
      type(esmf_state)                   :: imp_gfs_phy, exp_gfs_phy
      type(esmf_clock)                   :: clock
!
! !output variables and parameters:
!----------------------------------

      integer, intent(out) :: rc

! !eop
!------------------------------------------------------------------------- 
 
! !working arrays and local parameters.  
!--------------------------------------
      type(gfs_phy_wrap)                :: wrap
! this wrap is a derived type which contains
! only a pointer to the internal state.  it is needed
! for using different architectures or compliers.
      type(gfs_physics_internal_state), pointer  :: int_state    
      type(esmf_vm)                      :: vm_local
!jw
      type(esmf_state)                   :: imp_wrt_state
      type(esmf_timeinterval)            :: timestep, runduration
      type(esmf_time)                    :: starttime, stoptime, currtime
      type(esmf_timeinterval)            :: reftimeinterval
      type(esmf_delayout)                :: mydelayout
      integer(kind=esmf_kind_i4)         :: yy, mm, dd   ! time variables for date
      integer(kind=esmf_kind_i4)         :: hh, mns, sec ! time variables for time
      integer                            :: advancecount4
      integer                            :: atm_timestep_s, phy_timestep_s
      integer(esmf_kind_i8)              :: advancecount

      TYPE(ESMF_DistGrid)                :: DistGrid5    ! the ESMF DistGrid.

      integer                            :: rc1, rcfinal, grib_inp, ifhmax, &
                                            runduration_hour

      type(ESMF_FieldBundle)             :: Bundle     ! debug check

! initialize the error signal variables.
!---------------------------------------
      rc1     = esmf_success
      rcfinal = esmf_success

! allocate the internal state pointer.
!-------------------------------------
      call esmf_logwrite("allocate the internal state",                 &
                          ESMF_LOGMSG_INFO, rc = rc1)

      allocate(int_state, stat = rc1)

      call gfs_physics_err_msg(rc1,' - allocate the internal state',rc)

      wrap%int_state => int_state

! attach internal state to the gfs physics grid component.
!-------------------------------------------------
      call esmf_logwrite("set up the internal state",                   &
                          ESMF_LOGMSG_INFO, rc = rc1)

      call esmf_gridcompsetinternalstate(gc_gfs_phy, wrap, rc1)

      call gfs_physics_err_msg(rc1,'set up the internal state',rc)

! use esmf utilities to get information from the configuration file.
! the function is similar to reading the namelist in the original gfs.
!---------------------------------------------------------------------
      call esmf_logwrite("getting information from the configure file", &
                          ESMF_LOGMSG_INFO, rc = rc1)

      call gfs_physics_getcf(gc_gfs_phy, int_state,  rc = rc1)

      call gfs_physics_err_msg(rc1,'get configure file information',rc)
!jws
!-----------------------------------------------------------------------
!***  retrieve the import state of the write gridded component
!***  from the physics export state.
!-----------------------------------------------------------------------
      call esmf_logwrite("Retrieve Write Import State from Physics Export State", &
                          ESMF_LOGMSG_INFO, rc = rc1)
 
      CALL ESMF_StateGet(state      =exp_gfs_phy                        &  !<-- The Physics export state
                        ,itemName   ='Write Import State'               &  !<-- Name of the state to get from Physics export state
                        ,nestedState=imp_wrt_state                      &  !<-- Extract Write component import state from Physics export
                        ,rc         =RC1)
!
      CALL gfs_physics_err_msg(rc1,"Retrieve Write Import State from Physics Export State",RC)
!jwe
!
! initialize time interval to the parameter from the configure file.
!-------------------------------------------------------------------
!     call esmf_logwrite("set up time step interval",                   &
!                          ESMF_LOGMSG_INFO, rc = rc1)

!     call esmf_clockget(clock,            				&
!                        timestep    = timestep,                	&
!                        rc          = rc1)

!     call esmf_timeintervalget(timestep,                               &
!                               s  = atm_timestep_s,                    &
!                               rc = rc1)

!     phy_timestep_s = nint(int_state%nam_gfs_phy%deltim)

!     int_state%nam_gfs_phy%deltim = atm_timestep_s /			&
!           min( 1, atm_timestep_s / phy_timestep_s  )

!     call gfs_physics_err_msg(rc1,'set up time step interval',rc)

! get the start time from reading the surface file.
!----------------------------------------------------------
      call esmf_logwrite("getting the start time",                      &
                          ESMF_LOGMSG_INFO, rc = rc1)

      call gfs_physics_start_time_get(                                  &
                        yy, mm, dd, hh, mns, sec, int_state%kfhour,     &
!--get init time
                        fhini,                                          &
                        int_state%n3, int_state%nam_gfs_phy%sfc_ini,rc1)
 
      call gfs_physics_err_msg(rc1,'getting the start time',rc)
 
      advancecount4    = nint(real(int_state%kfhour) * 3600.0 /         &
                              int_state%nam_gfs_phy%deltim)
      int_state%phour  = advancecount4 *                                &
                         int_state%nam_gfs_phy%deltim / 3600.0
      fhour            = int_state%phour
!     write(0,*)' in phys gridcomp fhour=',fhour,' phour=',int_state%phour
      int_state%kfhour = nint(int_state%phour)
      int_state%restart_step = .false.
      kdt_start              = advancecount4
      params%kdt_start = kdt_start
      if (fhini /= 0 .and. int_state%restart_run) int_state%restart_step = .true.

! initialize the clock with the start time based on the information
! from calling starttimeget.
!------------------------------------------
      call esmf_logwrite("set up the esmf time",                        &
                          ESMF_LOGMSG_INFO, rc = rc1)

! in dynamics, we had already timeset, this is redo, update later.

      call esmf_timeset(starttime, yy = yy, mm = mm,  dd = dd,          &
                              h  = hh, m  = mns, s  = sec, rc = rc1)

      call gfs_physics_err_msg(rc1,'set up the esmf time',rc)

      call esmf_logwrite("set up the reference time interval",          &
                          ESMF_LOGMSG_INFO, rc = rc1)

      call esmf_timeintervalset(reftimeinterval, h = int_state%kfhour,  &
                           m = 0, rc = rc1)

! re-set up the start time based on the kfhour value in the sigma file.
!----------------------------------------------------------------------
!starttime = starttime + reftimeinterval

      call gfs_physics_err_msg(rc1,                                     &
                         'set up the reference time interval',rc)

! set the esmf clock which will control the gfs run do loop.
!--------------------------------------------------------------

! in dynamics, clock is set, this is redo, update later

      currtime = starttime + reftimeinterval
      call esmf_clockset(clock, currtime = currtime, rc = rc1)
!
! get the grid component vm.
! this esmf_gridcompget vm can be used at any where you need it.
!---------------------------------------------------------------
      call esmf_logwrite("get the local vm", ESMF_LOGMSG_INFO, rc = rc1)

      call esmf_vmgetcurrent(vm_local, rc = rc1)

      call gfs_physics_err_msg(rc1,'get the vm',rc)


! set up parameters of mpi communications.
! use esmf utility to get pe identification and total number of pes.
!-------------------------------------------------------------------
      call esmf_logwrite("get me and nodes from vm",                    &
                          ESMF_LOGMSG_INFO, rc = rc1)

      call esmf_vmget(vm_local, localpet = int_state%me,                &
                      mpicommunicator = mpi_comm_all,                   &
                      petcount = int_state%nodes,                       &
                      rc       = rc1)
      me = int_state%me

      call gfs_physics_err_msg(rc1,'get me and nodes from vm',rc)

! initialize the gfs, including set up the internal state
! variables and some local parameter short names, aloocate
! internal state arrays.
!---------------------------------------------------------
      call esmf_logwrite("run the gfs_physics_initialize",              &
                          ESMF_LOGMSG_INFO, rc = rc1)

! ----------------- gfs physics related initialize --------------------
! ----------------------------------------------------------------------
!      write(0,*)'in after init, size fhour_idate=',size(int_state%fhour_idate,1), &
!        size(int_state%fhour_idate,2),'idate size=',size(idate)                   &
!       ,' fhour=',fhour,' phour=',int_state%phour

      call gfs_physics_initialize(int_state, rc1)

! ----------------------------------------------------------------------
!
      call gfs_physics_err_msg(rc1,'run the gfs_physics_initialize',rc)

      call esmf_clockget(clock, timestep    = timestep,                 &
                         runduration = runduration,                     &
                         starttime   = starttime,                       &
                         currtime    = currtime,                        &
                         rc          = rc1)
!
!
      call esmf_timeintervalget(runduration,                            &
                                h = runduration_hour, rc = rc1)
!
!
!moor ifhmax = nint(int_state%nam_gfs_phy%fhmax)
      ifhmax = nint(fhmax)
      if(runduration_hour <= 0    .or.                                  &
          ifhmax /= 0             .and.                                 &
          ifhmax <= int_state%kfhour + runduration_hour) then
          ifhmax            = nint(fhmax)
          runduration_hour  = nint(fhmax) - nint(fhini)
          call esmf_timeintervalset(runduration,                        &
                                    h = runduration_hour, rc = rc1)
      end if
      if (runduration_hour < 0) then
        print *,' fhini=',fhini, ' > fhmax=',fhmax,' job aborted'
        call mpi_quit(444)
      endif
      stoptime = currtime  + runduration
                           
      call esmf_clockset(clock, stoptime = stoptime, rc = rc1)
!
      call esmf_timeintervalget(timestep, s = timestep_sec, rc = rc1)
!!
      if (me==0) then
        print *,' timestep_sec=',timestep_sec,' rc1=',rc1
        call out_para(real(timestep_sec))
        print *,' gsm physics will forecast ',runduration_hour,' hours',  &
                ' from hour ',int_state%kfhour,' to hour ',               &
                  runduration_hour+int_state%kfhour
      endif
!
!
      call synchro

!
! Create 3D Gaussian grid                                   
!-----------------------
!
      call gfs_physics_grid_create_Gauss3D(vm_local,int_state,DistGrid5,rc1) 

      call gfs_physics_err_msg(rc1,'gfs_physics_grid_create_gauss',rc)     

!
! Define Physics Import and Export states   
!
      if ( .not. int_state%grid_aldata ) then
        call gfs_physics_import2internal( imp_gfs_phy,    &         
                                          int_state, rc = rc1)      
        call gfs_physics_err_msg(rc1,'PHY INIT-call import2internal',rc)
      endif

      call gfs_physics_internal2export ( int_state,        &
                                         exp_gfs_phy, rc = rc1)
      call gfs_physics_err_msg(rc1,'PHY INIT-call internal2export',rc)

!! debug check
!     call ESMF_StatePrint( exp_gfs_phy, rc=rc)

!

! set pointer the gfs export fields in the internal state 
! to the esmf exprot state which is the public interface
! for other esmf grid components.
!-------------------------------------------------------
!     call esmf_logwrite("internal state link to esmf export state", 	&
!                          ESMF_LOGMSG_INFO, rc = rc1)

      int_state%fhour_idate(1,1)   = int_state%kfhour
      int_state%fhour_idate(1,2:5) = idate(1:4)

!      call gfs_physics_internal2export(gc_gfs_phy, int_state,  	&
!                                          exp_gfs_phy, rc = rc1)

!     call gfs_physics_err_msg(rc1,'internal state to esmf export state',rc)

!-------------------------------------------------------
!jw send all the head info to write tasks
!-------------------------------------------------------
!
        call point_physics_output_gfs(int_state,imp_wrt_state)
!
!*******************************************************************
! print out the final error signal variable and put it to rc.
!------------------------------------------------------------
      call gfs_physics_err_msg_final(rcfinal,                           &
                        'initialize from gfs physics grid comp.',rc)

      end subroutine gfs_phy_initialize





!----------------------------------------------------------------------
!bop
!
! !routine: gfs_phy_run --- 
!           main grid component routine to run the gfs physics.
!
! !description: this subroutine will run the most part computations 
!               of the gfs physics.
!
! !revision history:
!
!  november 2004     weiyu yang initial code.
!  may      2005     weiyu yang for the updated gfs version.
!  february 2006     moorthi
!  december 2007     juang
!  oct 12 2009       Sarah Lu, call gfs_physics_import2internal_mgrid and
!                    gfs_physics_internal2export_mgrid to associate imp/exp
!                    states with internal state grid_fld
!  oct 17 2009       Sarah Lu, debug print added to track imp/exp states
!  mar 05 2010       Sarah Lu, internal2export_mgrid is called in init step
!
! !interface:
!

      subroutine gfs_phy_run(gc_gfs_phy, imp_gfs_phy, exp_gfs_phy, clock, rc)

!*     use gfs_physics_states_mod
!
! !input variables and parameters:
!---------------------------------
      type(esmf_gridcomp)                :: gc_gfs_phy   
      type(esmf_state)                   :: imp_gfs_phy 
 
! !output variables and parameters:
!----------------------------------
      type(esmf_clock)                   :: clock
      type(esmf_timeinterval)            :: timestep, donetime    
      type(esmf_time)                    :: starttime, currtime, stoptime
      type(esmf_state)                   :: exp_gfs_phy
      integer,             intent(out)   :: rc   
!
!eop
!-------------------------------------------------------------------------

! local variables
      real(kind=ESMF_KIND_R8)           :: buffer(lonr,latr)
      real(kind=ESMF_KIND_R8)           :: bufferSplit(lonr,lats_node_r)
      integer                           :: kmsk(lonr,lats_node_r)

!    validation variables
     type(ESMF_Grid)                   :: grid
     type(ESMF_Field), save            :: validationField
     integer, save                     :: slice=1


!
! !working arrays and local parameters.
!--------------------------------------
      type(gfs_phy_wrap)                :: wrap         
! this wrap is a derived type which contains
! only a pointer to the internal state.  it is needed
! for using different architectures or compliers.
      type(gfs_physics_internal_state), pointer   :: int_state   
      integer                                     :: rc1, rcfinal
      real(8) :: zhour1
      real(4) :: wgt
      integer :: kint
!
!jw
      type(esmf_state)                   :: imp_wrt_state
!lu
      type(ESMF_Field)                   :: Field
      type(ESMF_FieldBundle)             :: Bundle
      integer                            :: i
      character*10                       :: vname
!
! initialize the error signal variables.
!---------------------------------------
      rc1     = esmf_success
      rcfinal = esmf_success

! retrieve the esmf internal state.
!---------------------------------- 
      call esmf_logwrite("get the internal state in the run routine",   &
                          ESMF_LOGMSG_INFO, rc = rc1)

      call esmf_gridcompgetinternalstate(gc_gfs_phy, wrap, rc1)

      call gfs_physics_err_msg(rc1,                                     &
                  'get the internal state in the run routine',rc)

! pointing the local internal state pointer to the esmf internal state pointer.
!------------------------------------------------------------------------------
      int_state => wrap%int_state


!     write(0,*)' int_state%lonsperlar=',int_state%lonsperlar
! get the esmf import state and over-write the gfs internal state.
! update the initial condition arrays in the internal state based on
! the information of the esmf import state. 
!------------------------------------------------------------------
      call esmf_logwrite("esmf import state to internal state",         &
                          ESMF_LOGMSG_INFO, rc = rc1)
!
! the pointer/copy option (Sarah Lu)
!  get the esmf import state and over-write the gfs internal state         
!  for one-copy option, import2internal is called in the init step
!  for two-copy option, import2internal is called every time step        
!
!*    call gfs_physics_import2internal(gc_gfs_phy, imp_gfs_phy, 	&
!*                                        int_state, rc = rc1)
!
!       print *,'in gfs phys grid comp, run'

      if ( int_state%grid_aldata ) then
        call gfs_physics_import2internal( imp_gfs_phy, int_state, rc = rc1)
        call gfs_physics_err_msg(rc1,'import2internal',rc)     
      endif                                                         

      idate(1:4)=int_state%fhour_idate(1,2:5)
      fhour     =int_state%fhour_idate(1,1)

      call gfs_physics_err_msg(rc1,'esmf import state to internal state',rc)

!
! get clock times
! ------------------
      call esmf_clockget(clock,                                         &
                         timestep    = timestep,                        &
                         starttime   = starttime,                       &
                         currtime    = currtime,                        &
                         stoptime    = stoptime,                        &
                         rc          = rc1)

      call gfs_physics_err_msg(rc1,'esmf clockget',rc)

      donetime      = currtime - starttime
      int_state%kdt = nint(donetime/timeStep) + 1

      if (params % ifp_realtime_interval .gt. 0 .and. &
          mod(int_state%kdt * timestep_sec, params % ifp_realtime_interval) .eq. 0 ) then
          call read_ifp
      end if

      kint = ((int_state%kdt - 1 + params % skip - params % kdt_start) * timestep_sec / params % ifp_interval) + 1
      write(6,*)'wic', int_state%kdt, timestep_sec, params % skip, params % kdt_start, kint
      if ( kint + 1 .le. size(farr % f107)) then
        wgt = 1 - real(mod((int_state%kdt-1)*timestep_sec, params % ifp_interval))/params % ifp_interval

        int_state % forcing % f107  = farr % f107 (kint) * wgt + farr % f107 (kint+1) * (1 - wgt)
        int_state % forcing % f107d = farr % f107d(kint) * wgt + farr % f107d(kint+1) * (1 - wgt)
        int_state % forcing % kp    = farr % kp   (kint) * wgt + farr % kp   (kint+1) * (1 - wgt)
        int_state % forcing % kpa   = farr % kpa  (kint) * wgt + farr % kpa  (kint+1) * (1 - wgt)
        int_state % forcing % nhp   = farr % nhp  (kint) * wgt + farr % nhp  (kint+1) * (1 - wgt)
        int_state % forcing % nhpi  = farr % nhpi (kint) * wgt + farr % nhpi (kint+1) * (1 - wgt)
        int_state % forcing % shp   = farr % shp  (kint) * wgt + farr % shp  (kint+1) * (1 - wgt)
        int_state % forcing % shpi  = farr % shpi (kint) * wgt + farr % shpi (kint+1) * (1 - wgt)
        int_state % forcing % swang = farr % swang(kint) * wgt + farr % swang(kint+1) * (1 - wgt)
        int_state % forcing % swden = farr % swden(kint) * wgt + farr % swden(kint+1) * (1 - wgt)
        int_state % forcing % swvel = farr % swvel(kint) * wgt + farr % swvel(kint+1) * (1 - wgt)
        int_state % forcing % swbz  = farr % swbz (kint) * wgt + farr % swbz (kint+1) * (1 - wgt)
        int_state % forcing % swbt  = farr % swbt (kint) * wgt + farr % swbt (kint+1) * (1 - wgt)
      else
        kint = size(farr % f107)

        int_state % forcing % f107  = farr % f107 (kint)
        int_state % forcing % f107d = farr % f107d(kint)
        int_state % forcing % kp    = farr % kp   (kint)
        int_state % forcing % kpa   = farr % kpa  (kint)
        int_state % forcing % nhp   = farr % nhp  (kint)
        int_state % forcing % nhpi  = farr % nhpi (kint)
        int_state % forcing % shp   = farr % shp  (kint)
        int_state % forcing % shpi  = farr % shpi (kint)
        int_state % forcing % swang = farr % swang(kint)
        int_state % forcing % swden = farr % swden(kint)
        int_state % forcing % swvel = farr % swvel(kint)
        int_state % forcing % swbz  = farr % swbz (kint)
        int_state % forcing % swbt  = farr % swbt (kint)
      end if
!      PRINT*, 'In phys grid comp, kdt, kdt_3h=', int_state%kdt, kdt_3h

!      write(0,*)' in physics kdt=',int_state%kdt,' kdt_start=',kdt_start &
!                ,' ndfi=',ndfi,' zhour=',int_state%zhour,                &
!                 ' zhour_dfi=',int_state%zhour_dfi
!      write(0,*)' in physics kdt=',int_state%kdt,' donetime=',&
!      donetime,' currtime=',currtime,' starttime=',starttime, stoptime

!     if( currtime .eq. stoptime ) then
!         print *,' currtime equals to stoptime '
!         int_state%end_step=.true.
!     endif
!
!-----------------------------------------------------------------------
!***  retrieve the import state of the write gridded component
!***  from the physics export state.
!-----------------------------------------------------------------------
      call esmf_logwrite("Retrieve Write Import State from Physics Export State", &
                          ESMF_LOGMSG_INFO, rc = rc1)


      CALL ESMF_StateGet(state      =exp_gfs_phy                        &  !<-- The Physics export state
                        ,itemName   ='Write Import State'               &  !<-- Name of the state to get from Physics export state
                        ,nestedState=imp_wrt_state                      &  !<-- Extract Write component import state from Physics export
                        ,rc         =RC1)
!
      CALL gfs_physics_err_msg(rc1,"Retrieve Write Import State from Physics Export State",RC)
!-----------------------------------------------------------------------
!      CALL ESMF_AttributeSet(state    =imp_wrt_state                    &  !<-- The Write component import state
!                           ,name     ='zhour'                          &  !<-- Name of the var
!                           ,value    =int_state%zhour                  &  !<-- The var being inserted into the impo
!                            ,rc       =RC)
      if(ldfi .and. int_state%kdt-kdt_start-1 == ndfi) then
         CALL ESMF_AttributeSet(state    =imp_wrt_state                    &  !<-- The Write component import state
                               ,name     ='zhour'                          &  !<-- Name of the var
                               ,value    =int_state%zhour_dfi              &  !<-- The var being inserted into the import state
                               ,rc       =RC)
      else
         CALL ESMF_AttributeSet(state    =imp_wrt_state                    &  !<-- The Write component import state
                               ,name     ='zhour'                          &  !<-- Name of the var
                               ,value    =int_state%zhour                  &  !<-- The var being inserted into the impo
                               ,rc       =RC)

      endif
!       write(0,*)'in physgrid comp,kdt=',int_state%kdt,'zhour=',int_state%zhour,'zhour_dfi=',  &
!           int_state%zhour_dfi,'ldfi=',ldfi,'ndfi=',ndfi,'rc=',rc
!
      CALL ESMF_AttributeSet(state=exp_gfs_phy        &  !<-- The physics export state
                            ,name ='kdt'              &  !<-- Name of the attribute to insert
                            ,value= int_state%kdt     & !<-- Value of the attribute
                            ,rc   =RC)

!

! Whenever the GSM run method is called, the "importFields" array holds the
! latest imported field data

  do i=1, NImportFields

    ! Each import field is only available if it was connected in the 
    ! import state.
    if (ESMF_FieldIsCreated(importFields(i))) then

      ! The Field is defined on a non-reduced Gaussian grid
      ! and unshuffled latitudes. Take 3 steps to change the SST data into the
      ! native GSM representation (reduced Gaussian grid with shuffled
      ! latitudes).

      ! Step 1: Gather the entire SST data onto local PET 0.

      buffer = -999000.*real(i) - real(me)  ! identifiable init value
      call ESMF_FieldGather(importFields(i), buffer, rootPet=0, rc=rc)
      ESMF_ERR_RETURN(rc,rc)
    
!print *, "importData",i," after gather:", me, minval(buffer), maxval(buffer)

      ! Step 2: Split the buffer across all PETs, taking into account the GSM
      ! distribution as well as latitude shuffle.

      bufferSplit = -888000.*real(i) - real(me) ! identifiable init value
      call split2d_phys_r8(buffer,bufferSplit,global_lats_ptr,0)

!print *, "importData",i," after split:", me, minval(bufferSplit), maxval(bufferSplit)

      ! Step 3: Reduce local blocks along the longitude dimension

      int_state%importData(:,:,i) = -777000.*real(i) - real(me) ! identifiable init value
      kmsk = 0  ! no masking
      call interpred_phys(1, kmsk, bufferSplit, int_state%importData(:,:,i), &
                          global_lats_ptr, lonsperlat_ptr)
                        
!print *, "importData",i," after reduce:", me, minval(int_state%importData(:,:,i)), maxval(int_state%importData(:,:,i))

      ! -> DONE: now "importData(:,:,i) holds the incoming data
    endif
  enddo
!
! run the gfs.
!--------------------------
      call esmf_logwrite("run the gfs_physics_run", ESMF_LOGMSG_INFO, rc = rc1)

      call gfs_physics_run(int_state, rc = rc1)
!
      call gfs_physics_err_msg(rc1,'run the gfs_physics_run',rc)
!
!
! transfer the gfs export fields in the internal state 
! to the esmf exprot state which is the public interface
! for other esmf grid components.
!-------------------------------------------------------
     call esmf_logwrite("internal state to esmf export state", &
                          ESMF_LOGMSG_INFO, rc = rc1)

! the pointer/copy option (Sarah Lu)
!  point export state to internal state grid_fld in init step

!   Moorthi testing on Nov 24 2015
!     if ( int_state%start_step ) then                             
!       call gfs_physics_internal2export( int_state,             & 
!                                          exp_gfs_phy, rc = rc1)  
!       call gfs_physics_err_msg(rc1,'internal2export',rc)  
!
!       int_state%start_step = .false. 
!     endif                                                     

!     call gfs_physics_err_msg(rc1,'internal state to esmf export state',rc)

!
!-----------------------------------------------------------------------
!***  retrieve the import state of the write gridded component
!***  from the physics export state.
!-----------------------------------------------------------------------
!      call esmf_logwrite("Retrieve Write Import State from Physics Export State", &
!                          ESMF_LOGMSG_INFO, rc = rc1)
!
!      CALL ESMF_StateGet(state      =exp_gfs_phy                        &  !<-- The Physics export state
!                        ,itemName   ='Write Import State'               &  !<-- Name of the state to get from Physics export state
!                        ,nestedState=imp_wrt_state                      &  !<-- Extract Write component import state from Physics export
!                        ,rc         =RC1)
!!
!      CALL gfs_physics_err_msg(rc1,"Retrieve Write Import State from Physics Export State",RC)
!-----------------------------------------------------------------------
!      CALL ESMF_AttributeSet(state    =imp_wrt_state                    &  !<-- The Write component import state
!                            ,name     ='zhour'                          &  !<-- Name of the var
!                            ,value    =int_state%zhour                  &  !<-- The var being inserted into the import state
!                            ,rc       =RC)
!
!*******************************************************************
!
! print out the final error signal information and put it to rc.
!---------------------------------------------------------------
      call gfs_physics_err_msg_final(rcfinal,                          &
                        'run from gfs physics grid comp.',rc)

      end subroutine gfs_phy_run


!----------------------------------------------------------------------
!bop
!
! !routine: finalize --- finalizing routine to finish the 
!                        gfs running job.
!
! !description: this subroutine will finish the gfs computations,
! !             and will release the memory space.
!
! !revision history:
!
!  november 2004     weiyu yang initial code.
!  may      2005     weiyu yang for the updated gfs version.
!  february 2006     moorthi
!  february 2007     juang for dynamics only
!  july     2007     juang for physics only
!
! !interface:

      subroutine gfs_phy_finalize(gc_gfs_phy,                           &
                                 imp_gfs_phy, exp_gfs_phy, clock, rc)

!
! !input variables and parameters:
!---------------------------------
      type(esmf_gridcomp)                 :: gc_gfs_phy
      type(esmf_state)                    :: imp_gfs_phy, exp_gfs_phy
      type(esmf_clock)                    :: clock

! !output variables and parameters:
!----------------------------------
      integer,             intent(out)    :: rc

! !working arrays and local parameters.
!--------------------------------------
      type(gfs_phy_wrap)                            :: wrap   
      type(gfs_physics_internal_state), pointer     :: int_state  
      integer                                       :: rc1, rcfinal

!eop
!-------------------------------------------------------------------------

! initialize the error signal variables.
!---------------------------------------
      rc1     = esmf_success
      rcfinal = esmf_success

! retrieve the esmf internal state.
!----------------------------------
     call esmf_logwrite(                                                &
                      "get the internal state in the finalize routine", &
                          ESMF_LOGMSG_INFO, rc = rc1)

     call esmf_gridcompgetinternalstate(gc_gfs_phy, wrap, rc1)

     call gfs_physics_err_msg(rc1,                                      &
              'get the internal state in the finalize routine',rc)

! point the local internal state pointer to the esmf internal state pointer.
!------------------------------------------------------------------------------
      int_state => wrap%int_state

! run the gfs finalize routine to release the memory space, etc. 
!----------------------------------------------------------------------------
      call esmf_logwrite("run the gfs_physics_finalize",                &
                          ESMF_LOGMSG_INFO, rc = rc1)

      call gfs_physics_finalize(int_state, rc = rc1)

      call gfs_physics_err_msg(rc1,'run the gfs_physics_finalize',rc)

! print out the final error signal information and put it to rc.
!---------------------------------------------------------------
      call gfs_physics_err_msg_final(rcfinal,                           &
                        'finalize from gfs physics grid comp.',rc)

      end subroutine gfs_phy_finalize

! end of the gfs esmf grid component module.
!-------------------------------------------
      end module gfs_physics_grid_comp_mod
