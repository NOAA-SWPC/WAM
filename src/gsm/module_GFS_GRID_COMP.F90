!#include "../../ESMFVersionDefine.h"

!----------------------------------------------------------------------
!
      MODULE MODULE_GFS_GRID_COMP
!
!-----------------------------------------------------------------------
!
!***  THIS IS THE GFS GRIDDED COMPONENT MODULE.
!***  IT WILL SET UP DYNAMICS, PHYSICS, AND COUPLER SUBCOMPONENTS
!***  AND RUN THEIR INITIALIZE, RUN, AND FINALIZE ROUTINES.
!
!-----------------------------------------------------------------------
!
! PROGRAM HISTORY LOG:
!   2007-       Black - Modified from Wei-yu's version
!   2007-11-20  Black/Henderson - Created an ESMF Clock for the
!                                 ATM Component independent of
!                                 the Main Clock.
!   2007-12-11  Black - Generalized for easier use by any dynamics core.
!   2008        Juang - create GFS grid component
!   2008-08     Colon - Added conditional checks multiple dynamics cores.
!   2008-10-14  Vasic - Added restart Alarm.
!   2009-05-29  Wang  - Added GFS write grid component
!   2009-07-23  Lu    - Added GOCART grid component
!   2009-08-03  Black - Merging with nesting.
!   2009-08-12  Black - Fixed logic for Physics export when direction of
!                       integration switches from backward to forward.
!   2009-10-05  Wang  - Added GFS ensemble member name and output data at
!                       every nsout timesteps.
!   2009-10-08  W Yang - Ensemble GEFS.
!   2009-11-03  Lu    - Add GOCART and ChemRegistry modules
!   2009-12-17  Lu    - Modify GFS_ATM_INIT routine to create, register,
!                       and initialize GOCART grid component
!   2009-12-22  Lu    - Modify GFS_ATM_INIT routine to create and register
!                       dyn2chem and phy2chem coupler components
!   2009-12-23  Lu    - Modify GFS_INTEGRATE routine to loop thru dyn, phy,
!   2010-02-01  Lu    - Remove dyn2chem coupler component
!   2010-02-05  Wang  - Added restart file for GFS
!   2010-03-04  Lu    - Modify GFS_ATM_INIT (initialization is changed from
!                       DYN-PHY-CPL to DYN-CPL-PHY)
!   2010-03-05  Lu    - Add GOCART_SETUP (to create and register GOCART) and
!                       GOCART_INIT (to initialize GOCART)
!   2010-03-23  Lu    - Add passive_tracer option
!   2010-08-17  Lu    - Make debug print optional
!   2010-08-25  Wang  - Add 3hr dfi filtered fields output option
!   2010-12-16  Wang  - Change to nemsio library
!   2011-02   W Yang  - Updated to use both the ESMF 4.0.0rp2 library,
!                       ESMF 5 series library and the the
!                       ESMF 3.1.0rp2 library.
!   2011-05-11  Yang  - Modified for using the ESMF 5.2.0r_beta_snapshot_07.
!   2011-10-01  Wang/Lu  - MYPE added to GOCART_INIT argument
!   2012-02-08  Yang  - Modified for using the ESMF 5.2.0rp1.
!   2012-11-06  Wang  - separate io from gfs for generalization
!   2013-01-12  Moorthi - Initialize total_integ_tim
!   2015-10-22  S Moorthi - Modify digital filter to work for restart step
!                           e.g. second segment of forecast
!   2016-01-12  Wang    - pass config variable "adiabatic" to dyn/phy imp/exp state
!   2016-01-13  P Tripp - NUOPC/GSM merge - NUOPC_ClockPrintStartTime to ESMF_ClockPrint
!   2016-01-30  Wang    - add high frequency output
!
! USAGE: GFS Gridded component parts called from subroutines within
!        module_ATM_GRID_COMP.F90.
!
!-----------------------------------------------------------------------
!
      USE ESMF
      use NUOPC
      USE MODULE_INCLUDE
!
      USE MODULE_GFS_INTERNAL_STATE ,ONLY: GFS_INTERNAL_STATE            &
                                          ,WRAP_GFS_INTERNAL_STATE
!
      USE NEMSIO_MODULE
! 
      USE GFS_DYNAMICS_GRID_COMP_MOD,ONLY: GFS_DYN_SETSERVICES
      USE GFS_PHYSICS_GRID_COMP_MOD ,ONLY: GFS_PHY_SETSERVICES
      USE MODULE_GFS_INTEGRATE      ,ONLY: GFS_INTEGRATE
      USE MODULE_GFS_WRITE          ,ONLY: WRITE_ASYNC_GFS
!
      USE gfs_dyn_phy_cpl_comp_mod,  ONLY: gfs_cpl_setservices
!
      USE MODULE_GFS_MPI_DEF        ,ONLY: PETLIST_FCST
      USE MODULE_IO_MPI_DEF         ,ONLY: ENSMEM_NAME, WRITE_GROUPS
!
      USE MODULE_ERR_MSG            ,ONLY: ERR_MSG,MESSAGE_CHECK
!
!-----------------------------------------------------------------------
!***  LIST MODULES FOR GSFC CHEMISTRY PACKAGE
!-----------------------------------------------------------------------
!
      USE MODULE_GOCART_ROUTINES    ,ONLY: GOCART_SETUP, GOCART_INIT
!
!-----------------------------------------------------------------------
!***  LIST OTHER MODULES WITH NON-GENERIC ROUTINES USED BY GFS.
!-----------------------------------------------------------------------
!
      USE MODULE_GFS_WRITE          ,ONLY: WRITE_INIT_GFS             & !<-- These are routines used only when asynchronous
                                          ,WRITE_SETUP_GFS            & !    quilting is specified by the user in the
                                          ,WRITE_DESTROY_GFS            !    configure file for history output.
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
      PRIVATE
!
      PUBLIC :: GFS_REGISTER 
!
!-----------------------------------------------------------------------
!
!---------------------------------
!***  For determining clocktimes.
!---------------------------------
!
      TYPE(GFS_INTERNAL_STATE), POINTER, SAVE :: gfs_int_state
      TYPE(WRAP_GFS_INTERNAL_STATE),     SAVE :: WRAP
      TYPE(ESMF_TimeInterval),           SAVE :: TIMESTEP
      TYPE(ESMF_VM),                     SAVE :: VM, VM_LOCAL           !<-- The ESMF virtual machine.
!WY bug fix.
!-----------
      real(kind=8)                            :: total_integ_tim, btim0
      INTEGER      :: NFHMAX_HF

      logical      :: LWRTGRDCMP               !<-- run wrt grid comp
!
!-----------------------------------------------------------------------
!
      character(len=160) :: nuopcMsg
      CONTAINS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE GFS_REGISTER(GFS_GRID_COMP,RC_REG)
! 
!-----------------------------------------------------------------------
!***  Register the GFS gridded component's initialize, run, and finalize
!***  routines.
!-----------------------------------------------------------------------
!
      TYPE(ESMF_GridComp) :: GFS_GRID_COMP                !<-- GFS gridded component
!
      INTEGER,INTENT(OUT) :: RC_REG                       !<-- Return code for register
!     
      INTEGER             :: RC                           ! local variable
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC     = ESMF_SUCCESS   ! Error signal variable
      RC_REG = ESMF_SUCCESS   ! Error signal variable
!
!-----------------------------------------------------------------------
!***  Register the GFS INITIALIZE subroutine.  Since it is just one 
!***  subroutine, no need to specify phase.  The second argument is
!***  a pre-defined subroutine type, such as ESMF_SETINIT, ESMF_SETRUN, 
!***  or ESMF_SETFINAL.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK = "Create/Load the Configure Object"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!
      CALL ESMF_GridCompSetEntryPoint(GFS_GRID_COMP                  &  !<-- GFS gridded component
                                     ,ESMF_METHOD_INITIALIZE         &  !<-- Subroutine type
                                     ,GFS_INITIALIZE                 &  !<-- User's subroutine name
                                     ,rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
!-----------------------------------------------------------------------
!***  Register the Run step of the GFS component.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
      MESSAGE_CHECK = "Set 1st Entry Point for GFS Run"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
      CALL ESMF_GridCompSetEntryPoint(GFS_GRID_COMP         &  !<-- GFS gridded component
                                     ,ESMF_METHOD_RUN       &  !<-- Subroutine type
                                     ,GFS_RUN               &  !<-- The primary Dynamics / Physics /Coupler sequence
                                     ,rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
!-----------------------------------------------------------------------
!***  Register the GFS FINALIZE subroutine.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
      MESSAGE_CHECK = "Set Entry Point for GFS Finalize"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
      CALL ESMF_GridCompSetEntryPoint(GFS_GRID_COMP         &  !<-- GFS gridded component
                                     ,ESMF_METHOD_FINALIZE  &  !<-- Subroutine type
                                     ,GFS_FINALIZE          &  !<-- User's subroutine name
                                     ,rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
!-----------------------------------------------------------------------
!***  Check the error signal variable and print out the result.
!-----------------------------------------------------------------------
!
      IF(RC_REG == ESMF_SUCCESS) THEN
!       WRITE(0,*)' GFS_SET_SERVICES SUCCEEDED'
      ELSE
        WRITE(0,*)' GFS_SET_SERVICES FAILED  RC_REG=',RC_REG
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE GFS_REGISTER
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE GFS_INITIALIZE(GFS_GRID_COMP, IMP_STATE, EXP_STATE,    &
                                CLOCK_ATM,     RC_INIT)
!
!-----------------------------------------------------------------------
!***  THIS ROUTINE SETS UP FUNDAMENTAL ASPECTS OF THE MODEL RUN.
!-----------------------------------------------------------------------
!
      USE MODULE_GFS_CORE_SETUP,ONLY: GFS_SETUP

!------------------------
!***  Argument variables
!------------------------
!
      TYPE(ESMF_GridComp)               :: GFS_GRID_COMP                !<-- The GFS gridded component
      TYPE(ESMF_State)                  :: IMP_STATE                 &  !<-- The GFS component's import state
                                          ,EXP_STATE                    !<-- The GFS component's export state
      TYPE(ESMF_Clock)                  :: CLOCK_ATM                    !<-- The ESMF Clock from the ATM component.
      INTEGER            ,INTENT(OUT)   :: RC_INIT                      !<-- Return code for Initialize step
!
!---------------------
!***  Local variables
!---------------------

      INTEGER            :: MEMBER_ID,MODENS,NENS                    &
                           ,NFHOUT,NFMOUT,NFSOUT,NSOUT,NTSD,IRTN     &
                           ,NFHOUT_HF,NFMOUT_HF,NFSOUT_HF            &
                           ,TOTAL_MEMBER,TOTAL_TASKS,IERR,NFHOUR
!
      REAL               :: DELTIM
!
      CHARACTER(50)      :: MODE
!
      TYPE(ESMF_Config)  :: CF                                          !<-- The configure object for all GFS members
!
      TYPE(ESMF_Grid)    :: GRID_GFS_DYN                             &  !<-- The ESMF grid for the integration attached to
                                                                        !     the GFS dynamics gridded component.
                          , GRID_GFS_PHY                             &  !<-- The ESMF grid for the integration attached to
                                                                        !     the GFS physics gridded component.
                          , GRID_GFS_ATM                                !<-- The ESMF grid for the integration attached to
                                                                        !     the GFS ATM gridded component
!* restart file
!
      TYPE(ESMF_Time)    :: CURRTIME                                 &  !<-- The ESMF current time.
                           ,STARTTIME                                   !<-- The ESMF start time.
      LOGICAL            :: RESTARTED_RUN, ADIABATIC                    !<-- Original/restarted run logical flag
!
      INTEGER            :: IYEAR_FCST                               &  !<-- Current year from restart file
                           ,IMONTH_FCST                              &  !<-- Current month from restart file
                           ,IDAY_FCST                                &  !<-- Current day from restart file
                           ,IHOUR_FCST                               &  !<-- Current hour from restart file
                           ,IMINUTE_FCST                             &  !<-- Current minute from restart file
                           ,ISECOND_FCST                                !<-- Current second from restart file
!
      INTEGER,DIMENSION(7)     :: FCSTDATE
!
      INTEGER(ESMF_KIND_I8)    :: NTSD_START                            !<-- Timestep count (>0 for restarted runs)
      INTEGER                  :: MYPE_GLOBAL                           !<-- Each MPI task ID
!
      REAL(kind=KFPT)          :: SECOND_FCST                           !<-- Current second from restart file
      CHARACTER(64)            :: RESTART_FILENAME
!
      TYPE(NEMSIO_GFILE)       :: GFILE
                                                                        !     the GFS ATM gridded component.
      INTEGER, DIMENSION(100)  :: pe_member
      CHARACTER(12)            :: PELAB
      INTEGER(kind=KINT)       :: TIMESTEP_SEC_WHOLE                 &
                                , TIMESTEP_SEC_NUMERATOR             &
                                , TIMESTEP_SEC_DENOMINATOR           &
                                , fhrot                              & !<-- restart time
                                , i, j, i1, rc
!
      real(kind=8)             :: rtim0, rtim1, rtc
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC      = ESMF_SUCCESS
      RC_INIT = ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  Initialize timing variables.
!-----------------------------------------------------------------------
!
      rtim0           = rtc()
      btim0           = rtc()
      total_integ_tim = 0.
!
!-----------------------------------------------------------------------
!***  Allocate the GFS component's internal state, point at it,
!***  and attach it to the GFS component.
!-----------------------------------------------------------------------
!
      ALLOCATE(gfs_int_state, stat = RC)
      wrap%GFS_INT_STATE => gfs_int_state
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK = "Set the GFS Internal State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompSetInternalState(GFS_GRID_COMP, WRAP, RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC, MESSAGE_CHECK, RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  For the moment, use a direct alias to the ATM Clock.
!-----------------------------------------------------------------------
!
      gfs_int_state%CLOCK_GFS = CLOCK_ATM
!
!-----------------------------------------------------------------------
!***  Attach the configure file to the GFS component.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK = "Create Configure Object for GFS Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CF = ESMF_ConfigCreate(rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC, MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK = "Attach Configure File to the GFS Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK = "Load Configure File inot Configure Object"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigLoadFile(config=CF, filename='atm_namelist.rc', rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC, MESSAGE_CHECK, RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK = "Attach Configure File to the GFS Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompSet(gridcomp = GFS_GRID_COMP                 &  !<-- The GFS gridded component
                           ,config   = CF                            &  !<-- The configure object (~namelist)
                           ,rc       = RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC, MESSAGE_CHECK, RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Retrieve the VM from the GFS component.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK = "Retrieve the CF and VM from GFS Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompGet(gridcomp = GFS_GRID_COMP                 &  !<-- The GFS gridded component
                           ,VM       = VM_LOCAL                      &  !<-- The Virtual Machine
                           ,rc       = RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC, MESSAGE_CHECK, RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Retrieve global VM then the total number of tasks for
!***  then entire system.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK = "Retrieve global VM for GFS"
!     CALL ESMF_LogWrite(MESSAGE_CHECK, ESMF_LOGMSG_INFO, rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

      CALL ESMF_VMGetGlobal(vm=VM, rc=RC)                             !<-- The virtual machine
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC, MESSAGE_CHECK, RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK = "GFS_INITIALIZE: Obtain MPI Task IDs from VM"
!     CALL ESMF_LogWrite(MESSAGE_CHECK, ESMF_LOGMSG_INFO, rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_VMGet(vm       = VM                                  &  !<-- The virtual machine
                     ,pecount  = TOTAL_TASKS                         &  !<-- # of MPI tasks for entire GFS system
                     ,localpet = MYPE_GLOBAL                         &  !<-- Each MPI task ID
                     ,rc       = RC)

      CALL ESMF_VMGet(vm       = VM_LOCAL                            &  !<-- The virtual machine
                     ,localpet = gfs_int_state%MYPE                  &  !<-- Each MPI task ID
                     ,rc       = RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC, MESSAGE_CHECK, RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Obtain the total number of GFS ensemble members from the 
!***  configure file then create the members' IDs and names.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK = "GFS_INITIALIZE: Extract # of members from Configure File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK, ESMF_LOGMSG_INFO, rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                         &  !<-- The configure object
                                  ,value =TOTAL_MEMBER               &  !<-- Fill this variable 
                                  ,label ='total_member:'            &  !<-- Give the variable this label's value from the config file
                                  ,rc    =RC)

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      pe_member = 0
      DO i = 1, TOTAL_MEMBER
          WRITE(PELAB, '("PE_MEMBER", I2.2, ":")') i
          CALL ESMF_ConfigGetAttribute(Cf, pe_member(i), label = PELAB, rc = RC)
          IF(pe_member(i) == 0) pe_member(i) = TOTAL_TASKS / TOTAL_MEMBER
      END DO

      i1 = 0
      DO j = 1, TOTAL_MEMBER
          DO i = 1, pe_member(j)
              IF(MYPE_GLOBAL == i1) THEN
                  member_id = j
              END IF
              i1 = i1 + 1
          END DO
      END DO

      IF(total_member == 1) THEN
          ENSMEM_NAME = ' '
      ELSE
          WRITE(ENSMEM_NAME, '("_",i2.2)') MEMBER_ID
      END IF

!-----------------------------------------------------------------------
!***  Extract fundamental timestep information from the config file.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK = "Extract Timestep Information from GFS Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK, ESMF_LOGMSG_INFO, rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                         &  !<-- The config object
                                  ,value =TIMESTEP_SEC_WHOLE         &  !<-- The variable filled (integer part of timestep (sec))
                                  ,label ='dt_int:'                  &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
      CALL ESMF_ConfigGetAttribute(config=CF                         &  !<-- The config object
                                  ,value =TIMESTEP_SEC_NUMERATOR     &  !<-- The variable filled (numerator of timestep fraction)
                                  ,label ='dt_num:'                  &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
      CALL ESMF_ConfigGetAttribute(config=CF                         &  !<-- The config object
                                  ,value =TIMESTEP_SEC_DENOMINATOR   &  !<-- The variable filled (denominator of timestep fraction)
                                  ,label ='dt_den:'                  &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC, MESSAGE_CHECK, RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Establish the timestep for the GFS Clock.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK = "Set Time Step Interval in GFS Clock"
!     CALL ESMF_LogWrite(MESSAGE_CHECK, ESMF_LOGMSG_INFO, rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_TimeIntervalSet(timeinterval=TIMESTEP                 &  !<-- GFS clock's fundamental timestep (sec) (ESMF)
                               ,s           =TIMESTEP_SEC_WHOLE       &  !<-- Whole part of timestep
                               ,sn          =TIMESTEP_SEC_NUMERATOR   &  !<-- Numerator of fractional part
                               ,sd          =TIMESTEP_SEC_DENOMINATOR &  !<-- Denominator of fractional part
                               ,rc          =RC)

      CALL ESMF_ClockSet(clock=gfs_int_state%CLOCK_GFS, timeStep=TIMESTEP, rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC, MESSAGE_CHECK, RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Model-specific routines must be invoked in order to establish
!***  the ESMF Grid.  The different integration grids necessitate
!***  different ways of setting up both the parallelism for
!***  distributed memory runs and the ESMF Grid itself.
!***  When the parallelism is constructed, the local domain limits
!***  need to be inserted into the ATM component's internal state
!***  if quilting is to be used.
!-----------------------------------------------------------------------
!
      CALL GFS_SETUP(GFS_GRID_COMP, GRID_GFS_ATM)
!
!-----------------------------------------------------------------------
!***  Now establish the frequency of forecast output.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK = "Extract History Output Interval from GFS Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure object
                                  ,value =NFHOUT                        &  !<-- Hours between GFS history output
                                  ,label ='nfhout:'                     &  !<-- Give the variable this label's value
                                  ,rc    =RC)
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure object
                                  ,value =NFHOUT_HF                     &  !<-- Hours between GFS high frequency history output
                                  ,label ='nfhout_hf:'                  &  !<-- Give the variable this label's value
                                  ,rc    =RC)
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure object
                                  ,value =NFHMAX_HF                     &  !<-- max Hours for GFS high frequency history output
                                  ,label ='nfhmax_hf:'                  &  !<-- Give the variable this label's value
                                  ,rc    =RC)

      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure object
                                  ,value =NSOUT                         &  !<-- Fill this variable bel's value from the config file
                                  ,label ='nsout:'                      &  !<-- Give the variable this label's value from the config file
                                  ,rc    =RC)
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure object
                                  ,value =DELTIM                        &  !<-- Fill this variable
                                  ,label ='deltim:'                     &  !<-- Give the variable this label's value from the config file
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC, MESSAGE_CHECK, RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF(NSOUT>0) THEN
        NFHOUT = INT(NSOUT*DELTIM/3600.)
        NFMOUT = INT((NSOUT*DELTIM-NFHOUT*3600.)/60.)
        NFSOUT = INT(NSOUT*DELTIM-NFHOUT*3600.-NFMOUT*60)
      ELSE
        NFMOUT = 0
        NFSOUT = 0
      ENDIF
      IF (gfs_int_state%MYPE == 0 )                                          &
      write(0,*)'nfhout=',nfhout,'nfmout=',nfmout,'nfsout=',nfsout
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK = "Set GFS History Output Interval"
!     CALL ESMF_LogWrite(MESSAGE_CHECK, ESMF_LOGMSG_INFO, rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_TimeIntervalSet(gfs_int_state%TIMEINTERVAL_GFS_OUTPUT         &  !<-- ESMF time interval between GFS history output
                               ,h           =NFHOUT                           &  !<-- Hours between GFS history output
                               ,m           =NFMOUT                           &  !<-- Minutes between GFS history output
                               ,s           =NFSOUT                           &  !<-- Seconds between GFS history output
                               ,rc          =RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set GFS History High Frequency Output Interval"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF(NFHMAX_HF > 0) THEN
         NFMOUT_HF = 0
         NFSOUT_HF = 0
         CALL ESMF_TimeIntervalSet(gfs_int_state%TIMEINTERVAL_GFS_OUTPUT_HF      &  !<-- ESMF time interval between GFS history output
             ,h           =NFHOUT_HF                                             &  !<-- Hours between GFS history output
             ,m           =NFMOUT_HF                                             &  !<-- Minutes between GFS history output
             ,s           =NFSOUT_HF                                             &  !<-- Seconds between GFS history output
             ,rc          =RC)
      END IF
!
!-----------------------------------------------------------------------
      CALL ERR_MSG(RC, MESSAGE_CHECK, RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Extract the start time from the clock.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK = "NMM_ATM_INIT: Start Time from GFS Clock"
!     CALL ESMF_LogWrite(MESSAGE_CHECK, ESMF_LOGMSG_INFO, rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ClockGet(clock    =gfs_int_state%CLOCK_GFS              &  !<-- The ESMF Clock of this domain
                        ,startTime=STARTTIME                            &  !<-- The simulation start time
                        ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC, MESSAGE_CHECK, RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CURRTIME   = STARTTIME
      NTSD_START = 0
!
!-----------------------------------------------------------------------
!***  Extract the RESTART flag from the configure file.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK = "Extract Restart Flag from Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK, ESMF_LOGMSG_INFO, rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The config object
                                  ,value =RESTARTED_RUN                 &  !<-- True => restart; False => cold start
                                  ,label ='restart:'                    &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC, MESSAGE_CHECK, RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!     if(restarted_run)print *,'restart_run is true'
!
!-----------------------------------------------------------------------
!***  If this is a restarted run then read:
!***    (1) The forecast time that the file was written.
!***    (2) The forecast timestep at which the file was written.
!-----------------------------------------------------------------------
!
      MESSAGE_CHECK = "GFS get FHROT from CF"
      CALL ESMF_ConfigGetAttribute(config=CF, value=FHROT, label='fhrot:',rc=RC)
      CALL ERR_MSG(RC, MESSAGE_CHECK, RC_INIT)


      restart: IF(RESTARTED_RUN .or. fhrot > 0) THEN                       !<-- If this is a restarted run, set the current time
!
        RESTART_FILENAME = 'grid_ini'
        CALL NEMSIO_INIT()
        CALL NEMSIO_OPEN(GFILE,trim(RESTART_FILENAME),'read',iret=IRTN)
        if ( irtn /= 0 )                                       &
         print *,'after nesmio open restartfile, irtn=',irtn
!
        CALL NEMSIO_GETHEADVAR(GFILE,'FCSTDATE',FCSTDATE,iret=irtn)
!
        IYEAR_FCST   = FCSTDATE(1)
        IMONTH_FCST  = FCSTDATE(2)
        IDAY_FCST    = FCSTDATE(3)
        IHOUR_FCST   = FCSTDATE(4)
        IMINUTE_FCST = FCSTDATE(5)
        SECOND_FCST  = 0.
!
        IF(FCSTDATE(7) /= 0 )THEN
          SECOND_FCST = FCSTDATE(6)/(FCSTDATE(7)*1.)
        ENDIF
!
        NTSD = 0
        CALL NEMSIO_GETHEADVAR(gfile,'NTIMESTEP',NTSD,iret=irtn)
        if ( irtn /= 0 ) then    
           CALL NEMSIO_GETHEADVAR(gfile,'NFHOUR',NFHOUR,iret=irtn)
           if ( irtn == 0 .and. NFHOUR > 0 ) then
             NTSD = NFHOUR*3600./deltim
           endif
!           print *,'from history file,ntsd=',ntsd,'NFHOUR=',NFHOUR,'deltim=',deltim
        endif
!
        CALL NEMSIO_CLOSE(GFILE,iret=IERR)
        CALL NEMSIO_finalize()
!
        ISECOND_FCST = NINT(SECOND_FCST)                                   !<-- ESMF clock needs integer seconds
        NTSD_START   = NTSD
      IF (gfs_int_state%MYPE == 0 )                                          &
        write(0,*)' NTSD_START=',NTSD_START,' NTSD=',NTSD
!
!~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK = "RESTART: Set the Current Time of the Forecast"
!       CALL ESMF_LogWrite(MESSAGE_CHECK, ESMF_LOGMSG_INFO, rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_TimeSet(time=CURRTIME                                 &  !<-- Current time of the forecast (ESMF)
                         ,yy  =IYEAR_FCST                               &  !<-- Year from restart file
                         ,mm  =IMONTH_FCST                              &  !<-- Month from restart file
                         ,dd  =IDAY_FCST                                &  !<-- Day from restart file
                         ,h   =IHOUR_FCST                               &  !<-- Hour from restart file
                         ,m   =IMINUTE_FCST                             &  !<-- Minute from restart file
                         ,s   =ISECOND_FCST                             &  !<-- Second from restart file
                         ,rc  =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC, MESSAGE_CHECK, RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
      ENDIF restart
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK = "Set the Current Time on the GFS Clock"
!     CALL ESMF_LogWrite(MESSAGE_CHECK, ESMF_LOGMSG_INFO, rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ClockSet(clock       =gfs_int_state%CLOCK_GFS           &  !<-- The ATM Component's Clock
                        ,currtime    =CURRTIME                          &  !<-- Current time of simulation
                        ,advanceCount=NTSD_START                        &  !<-- Timestep at this current time
                        ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC, MESSAGE_CHECK, RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF( ESMF_ClockIsStopTime(gfs_int_state%CLOCK_GFS,rc=RC)) then
         print *,'in GFS_INITIALIZE, It is stop time'
      ENDIF
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  Create the Dynamics gridded subcomponent.
!***  Register the Initialize, Run, and Finalize steps for it.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-------------------------------
!***  Create Dynamics component
!-------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK = "Create the GFS Dynamics Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      gfs_int_state%GC_GFS_DYN = ESMF_GridCompCreate                 &
                                 (name      ="dynamics component"    &
                                 ,configFile='dyn_namelist.rc'       &
                                 ,petList   =PETLIST_FCST            &
                                 ,rc        =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC, MESSAGE_CHECK, RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!------------------------------------------------
!***  Register the Init, Run, and Finalize steps
!------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK = "Register GFS Dynamics Init, Run, Finalize"
!     CALL ESMF_LogWrite(MESSAGE_CHECK, ESMF_LOGMSG_INFO, rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompSetServices(gfs_int_state%GC_GFS_DYN            &  !<-- The GFS Dynamics gridded component
                                   ,GFS_DYN_SETSERVICES                 &  !<-- The user's subroutineName for Register
                                   ,rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC, MESSAGE_CHECK, RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Create empty Import and Export states for the Dynamics component
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK = "Create Empty Import/Export States for Dynamics"
!     CALL ESMF_LogWrite(MESSAGE_CHECK, ESMF_LOGMSG_INFO, rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      gfs_int_state%IMP_GFS_DYN = ESMF_StateCreate                     &
                                 (name="GFS dynamics import"      &
                                 ,stateintent = ESMF_STATEINTENT_IMPORT&
                                 ,rc       =RC)
!
      gfs_int_state%EXP_GFS_DYN = ESMF_StateCreate                     &
                                 (name="GFS dynamics export"      &
                                 ,stateintent = ESMF_STATEINTENT_EXPORT&
                                 ,rc       =RC)
! Add the GFS dynamics ESMF states as the nested states into the ATM parent states.
!----------------------------------------------------------------------------------
      CALL ESMF_StateAddReplace(IMP_STATE, (/gfs_int_state%IMP_GFS_DYN/), rc = RC)
      CALL ESMF_StateAddReplace(EXP_STATE, (/gfs_int_state%EXP_GFS_DYN/), rc = RC)

      MESSAGE_CHECK = "GFS set Cpl_flag"
!     CALL ESMF_LogWrite(MESSAGE_CHECK, ESMF_LOGMSG_INFO, rc=RC)

      gfs_int_state%Cpl_flag = .false.

      CALL ESMF_AttributeSet(gfs_int_state%IMP_GFS_DYN,          &
                             'Cpl_flag',                         &
                             gfs_int_state%Cpl_flag,             &
                             rc = RC)

      CALL ERR_MSG(RC, MESSAGE_CHECK, RC_INIT)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC, MESSAGE_CHECK, RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Check if the run is adiabatic (i.e. no physics)
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK = "Extract Physics On/Off Switch from GFS Config Object"
!     CALL ESMF_LogWrite(MESSAGE_CHECK, ESMF_LOGMSG_INFO, rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF, value=adiabatic, label='adiabatic:',rc=rc)

      IF (gfs_int_state%MYPE == 0 )                                          &
        write(0,*)' from configure file adiabatic=',adiabatic
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC, MESSAGE_CHECK, RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      call ESMF_AttributeSet( gfs_int_state%IMP_GFS_DYN,"adiabatic",adiabatic,rc=rc)
      call ESMF_AttributeSet( gfs_int_state%EXP_GFS_DYN,"adiabatic",adiabatic,rc=rc)

      IF(adiabatic)THEN                                         !<-- Adiabatic => Physics off
        gfs_int_state%PHYSICS_ON = ESMF_False
        IF (gfs_int_state%MYPE == 0 )                                  &
        write(0,*)' Adiabatic run - Initialize without physics coupling'
      ELSE                                                                 !<-- Not adiabatic => Physics on
        gfs_int_state%PHYSICS_ON = ESMF_True
        IF (gfs_int_state%MYPE == 0 )                                  &
        write(0,*)' Initialize with physics coupling '
      ENDIF
!

!-----------------------------------------------------------------------
!***  Is this a passive_tracer (no chemistry) run?
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK = "Extract Chemistry On/Off Switch from GFS Config Object"
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF, value=MODE               &
                                  ,label='passive_tracer:', rc=rc)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC, MESSAGE_CHECK, RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF(TRIM(MODE) == '.true.')THEN                                         !<-- passive tracer => chemistry off
        gfs_int_state%CHEMISTRY_ON=ESMF_False
        IF (gfs_int_state%MYPE == 0 )                                        &
          write(0,*)' Initialize without chemistry coupling '
        ELSE                                                                 !<-- non-passive tracer => chemistry on
        gfs_int_state%CHEMISTRY_ON=ESMF_True
        IF (gfs_int_state%MYPE == 0 )                                        &
          write(0,*)' Initialize with chemistry coupling '
      ENDIF

                                         
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  Create the Physics gridded subcomponent if physics is turned on.
!***  Register the Initialize, Run, and Finalize steps for it.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      IF(gfs_int_state%PHYSICS_ON == ESMF_True) THEN
!
!-------------------------------
!***  Create Physics component
!-------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK = "Create the GFS Physics Component"
!       CALL ESMF_LogWrite(MESSAGE_CHECK, ESMF_LOGMSG_INFO, rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        gfs_int_state%GC_GFS_PHY = ESMF_GridCompCreate              &
                                  (name      ="physics component"   &
                                  ,configfile='phy_namelist.rc'     &
                                  ,petList   =PETLIST_FCST          &
                                  ,rc        =RC)
!       write(0,*)'in GFS_INITIALIZE after phys comp created, petlist_fcst=',petlist_fcst
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC, MESSAGE_CHECK, RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-------------------------------------------------
!***  Register the Init, Run, and Finalize steps.
!-------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK = "Register Physics Init, Run, Finalize"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_GridCompSetServices(gfs_int_state%GC_GFS_PHY     &
                                     ,GFS_PHY_SETSERVICES          &
                                     ,rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC, MESSAGE_CHECK, RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF
!
!------------------------------------------------------------------------
!***  Create empty Import and Export states for the Physics subcomponent.
!------------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK = "Create Empty Import/Export States for GFS Physics"
!     CALL ESMF_LogWrite(MESSAGE_CHECK, ESMF_LOGMSG_INFO, rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      gfs_int_state%IMP_GFS_PHY = ESMF_StateCreate                     &
                                 (name="physics import"           &
                                 ,stateintent = ESMF_STATEINTENT_IMPORT&
                                 ,rc       =RC)

      call ESMF_AttributeSet( gfs_int_state%IMP_GFS_PHY,"adiabatic",MODE,rc=rc)
!
      gfs_int_state%EXP_GFS_PHY = ESMF_StateCreate                     &
                                 (name="physics export"           &
                                 ,stateintent = ESMF_STATEINTENT_EXPORT&
                                 ,rc       =RC)
      call ESMF_AttributeSet( gfs_int_state%EXP_GFS_PHY,"adiabatic",MODE,rc=rc)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC, MESSAGE_CHECK, RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  Create the Dynamics-Physics coupler subcomponent.
!***  Register the Initialize, Run, and Finalize steps for it.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!----------------------------
!***  Create Dyn-Phy Coupler
!----------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK = "Create the GFS Dynamics-Physics Coupler Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      gfs_int_state%GC_GFS_CPL = ESMF_CplCompCreate          &
                                (name="coupler component"    &
                                ,petList=PETLIST_FCST, rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC, MESSAGE_CHECK, RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-------------------------------------------------
!***  Register the Init, Run, and Finalize steps.
!-------------------------------------------------
!
! ~ ~   ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK = "Register the Dyn-Phy Coupler's Init, Run, Finalize"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_CplCompSetServices(gfs_int_state%GC_GFS_CPL        &  !<-- The GFS Dynamics/Physics coupler component
                                  ,GFS_CPL_SETSERVICES             &  !<-- The user's subroutine name for Register
                                  ,rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC, MESSAGE_CHECK, RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  Setup GOCART grid component and PHY-CHEM/CHEM-PHY coupler components
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!

      IF(gfs_int_state%CHEMISTRY_ON == ESMF_True) THEN

        MESSAGE_CHECK = "Setup GOCART and PHY-CHEM/CHEM-PHY coupler"

        CALL GOCART_SETUP (gfs_int_state%GC_GFS_CHEM               &
                          ,gfs_int_state%IMP_GFS_CHEM              &
                          ,gfs_int_state%EXP_GFS_CHEM              &
                          ,gfs_int_state%GC_PHY2CHEM_CPL           &
                          ,gfs_int_state%GC_CHEM2PHY_CPL           &
                          ,gfs_int_state%CHEMISTRY_ON              &
                          ,gfs_int_state%MYPE                      &
                          ,RC                                      &
                          )

        CALL ERR_MSG(RC, MESSAGE_CHECK, RC_INIT)
      ENDIF
!
!-----------------------------------------------------------------------
!***  Will the Write components with asynchronous quilting be used?
!-----------------------------------------------------------------------

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK = "Extract Quilting Flag from GFS Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK, ESMF_LOGMSG_INFO, rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                       &  !<-- The GFS config object
                                  ,value =gfs_int_state%QUILTING   &  !<-- The quilting flag
                                  ,label ='quilting:'              &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC, MESSAGE_CHECK, RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  Setup the Write component(s) (which may run without quilting).
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      IF (gfs_int_state%MYPE == 0 )                                &
      write(0,*)'before write_setup_gfs, allocate,write_groups=',write_groups

      ALLOCATE(gfs_int_state%WRT_COMPS(WRITE_GROUPS))
      CALL WRITE_SETUP_GFS(GFS_GRID_COMP                           &
                          ,gfs_int_state%WRT_COMPS                 &
                          ,gfs_int_state%EXP_GFS_DYN               &
                          ,gfs_int_state%EXP_GFS_PHY               &
                          ,gfs_int_state%IMP_GFS_WRT               &
                          ,gfs_int_state%EXP_GFS_WRT)
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  Execute the Initialize steps for the gridded subcomponents.
!***  These are the Initialize subroutines specified in the
!***  Register routines called in ESMF_GridCompSetServices above.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!--------------
!***  DYNAMICS
!--------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK = "Initialize GFS Dynamics Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK, ESMF_LOGMSG_INFO, rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompInitialize(gfs_int_state%GC_GFS_DYN               &
                                  ,importstate=gfs_int_state%IMP_GFS_DYN  &
                                  ,exportstate=gfs_int_state%EXP_GFS_DYN  &
                                  ,clock      =gfs_int_state%CLOCK_GFS    &  
                                  ,rc         =RC)
!
      GRID_GFS_DYN = GRID_GFS_ATM                                            !<-- Use the ATM Grid for the Dynamics
!
      CALL ESMF_GridCompSet(gfs_int_state%GC_GFS_DYN, grid=GRID_GFS_DYN, rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC, MESSAGE_CHECK, RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!--------------
!***  DYN-PHY COUPLER COMPONENT
!--------------

!     IF(gfs_int_state%PHYSICS_ON == ESMF_True) THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK = "Initialize Dyn-Phy Coupler"
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!
        CALL ESMF_CplCompInitialize(cplcomp    =gfs_int_state%GC_GFS_CPL   &
                                   ,importstate=gfs_int_state%EXP_GFS_DYN  &
                                   ,exportstate=gfs_int_state%IMP_GFS_PHY  &
                                   ,clock      =gfs_int_state%CLOCK_GFS    &
                                   ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC, MESSAGE_CHECK, RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!-------------
!***  PHYSICS
!-------------
!
      IF(gfs_int_state%PHYSICS_ON == ESMF_True) THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK = "Initialize GFS Physics Component"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_GridCompInitialize(gfs_int_state%GC_GFS_PHY               &
                                    ,importstate=gfs_int_state%IMP_GFS_PHY  &
                                    ,exportstate=gfs_int_state%EXP_GFS_PHY  &
                                    ,clock      =gfs_int_state%CLOCK_GFS    &
                                    ,rc         =RC)
!
        GRID_GFS_PHY = GRID_GFS_ATM                                        !<-- Use the ATM Grid for the Physics
!
        CALL ESMF_GridCompSet(gfs_int_state%GC_GFS_PHY             &
                             ,grid    =GRID_GFS_PHY                &
                             ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!     ENDIF
!
!-------------
!***  CHEMISTRY
!-------------
!
        IF(gfs_int_state%CHEMISTRY_ON == ESMF_True) THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK = "Initialize GFS Chemistry Component"
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL GOCART_INIT ( gfs_int_state%GC_GFS_CHEM               &
                            ,gfs_int_state%EXP_GFS_PHY               &
                            ,gfs_int_state%IMP_GFS_CHEM              &
                            ,gfs_int_state%EXP_GFS_CHEM              &
                            ,gfs_int_state%GC_PHY2CHEM_CPL           &
                            ,gfs_int_state%GC_CHEM2PHY_CPL           &
                            ,gfs_int_state%CLOCK_GFS                 &
                            ,gfs_int_state%MYPE                      &
                            ,RC                                      &
                           )

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        ENDIF                   ! CHEMISTRY_ON
      ENDIF                     ! PHYSICS_ON
!
!-----------------------------------------------------------------------
!***  Execute the Initialize step of the Write component(s).
!-----------------------------------------------------------------------
!
      CALL WRITE_INIT_GFS(GFS_GRID_COMP                            &
                         ,gfs_int_state%WRT_COMPS                  &
                         ,gfs_int_state%IMP_GFS_WRT                &
                         ,gfs_int_state%EXP_GFS_WRT                &
                         ,gfs_int_state%CLOCK_GFS                  &
                         ,WRITE_GROUP_READY_TO_GO=                 &
                          gfs_int_state%WRITE_GROUP_READY_TO_GO)
!
!-----------------------------------------------------------------------
!***  WRITE THE FINAL ERROR SIGNAL.
!-----------------------------------------------------------------------
!
      IF(RC_INIT == ESMF_SUCCESS) THEN
!       WRITE(0,*)'GFS INITIALIZE STEP SUCCEEDED'
      ELSE
        WRITE(0,*)'GFS INITIALIZE STEP FAILED RC_INIT=',RC_INIT
      ENDIF
!
!-----------------------------------------------------------------------
!
      rtim1           = rtc()
      total_integ_tim = (rtc()-btim0)
!
      IF(gfs_int_state%MYPE == 0 )  &
          WRITE(0, *) ' GFS_init_tim=', total_integ_tim, (rtim1-rtim0)
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE GFS_INITIALIZE
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE GFS_RUN(GFS_GRID_COMP, IMP_STATE, EXP_STATE            &
                        ,CLOCK_ATM,     RC_RUN)
!
!-----------------------------------------------------------------------
!***  RUN THE GFS GRIDDED COMPONENT.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument variables
!------------------------
!
      TYPE(ESMF_GridComp)               :: GFS_GRID_COMP       !<-- GFS gridded component
      TYPE(ESMF_State)                  :: IMP_STATE         & !<-- GFS Run step's import
                                         , EXP_STATE           !<-- GFS Run step's export
      TYPE(ESMF_Clock)                  :: CLOCK_ATM           !<-- ATM ESMF Clock
      INTEGER,            INTENT(OUT)   :: RC_RUN              !<-- Return code for the Run step
!
!---------------------
!***  Local variables
!---------------------
!
      TYPE(ESMF_Config)          :: CF  
!
      real(kind=8)               :: rtim0, rtim1, rtc
!
      INTEGER(kind=KINT)         :: DFIHR                    &
                                   ,NTIMESTEP                & !<-- current forecast timestep (integer)
                                   ,DFILEVS                  & !<-- vertical lvl above which no change from digital filter
                                                               !    default: levs, wam:80
                                   ,fhrot                    & !<-- restart time
                                   ,RC                         !<-- Error signal variables.
!
      INTEGER(kind=ESMF_KIND_I8) :: NTIMESTEP_ESMF             !<-- current forecast timestep (ESMF) (integer)
!
      TYPE(ESMF_TimeInterval)    :: RUNDURATION              & !<-- forecast length (ESMF)
                                   ,TIMESTEP                 & !<-- fundamental timestep, seconds (ESMF)
                                   ,FHROTINTVAL
!
      TYPE(ESMF_Time)            :: CURRTIME                 & !<-- ESMF current time.
                                   ,STARTTIME                & !<-- ESMF start time.
                                   ,RESTRTTIME                 !<-- ESMF restart time
!
      LOGICAL                    :: LDFIFLTO                 & !<-- output filtered fields.
                                   ,LDFI_GRD                   !<-- grid point digital filter
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      rtim0           = rtc()
      btim0           = rtc()
      total_integ_tim = 0.0
!
!-----------------------------------------------------------------------
!***  For the moment, use a direct alias to the ATM Clock.
!-----------------------------------------------------------------------
!
      gfs_int_state%CLOCK_GFS = CLOCK_ATM
!
      call ESMF_ClockPrint(CLOCK_ATM, options="currTime", &
        preString="entering GFS_RUN with CLOCK current: ", &
        unit=nuopcMsg)
      call ESMF_LogWrite(nuopcMsg, ESMF_LOGMSG_INFO)
      call ESMF_ClockPrint(CLOCK_ATM, options="startTime", &
        preString="entering GFS_RUN with CLOCK start:   ", &
        unit=nuopcMsg)
      call ESMF_LogWrite(nuopcMsg, ESMF_LOGMSG_INFO)
      call ESMF_ClockPrint(CLOCK_ATM, options="stopTime", &
        preString="entering GFS_RUN with CLOCK stop:    ", &
        unit=nuopcMsg)
      call ESMF_LogWrite(nuopcMsg, ESMF_LOGMSG_INFO)
!
!-----------------------------------------------------------------------
!***  Extract the fundamental time information from the Clock.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK = "Retrieve GFS Timestep from the ATM Clock"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ClockGet(clock       = gfs_int_state%CLOCK_GFS     &  !<-- ESMF Clock
                        ,advanceCount=NTIMESTEP_ESMF               &  !<-- # of times the clock has advanced
                        ,timestep    =TIMESTEP                     &  !<-- model's timestep length
                        ,starttime   =STARTTIME                    &  !<-- forecast start time
                        ,currtime    =CURRTIME                     &  !<-- Clock's current time
                        ,runduration =RUNDURATION                  &  !<-- length of the forecast
                        ,rc          =RC)

!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC, MESSAGE_CHECK, RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      NTIMESTEP = NTIMESTEP_ESMF                                        !<-- Convert timestep from ESMF to integer
!
!-----------------------------------------------------------------------
!***  We need the DFI filter duration.  Extract the configure file 
!***  from the ATM component and then obtain the filter duration.
!-----------------------------------------------------------------------
!
      MESSAGE_CHECK = "GFS get DFIHR from CF"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)

      CALL ESMF_GridCompGet(gridcomp=GFS_GRID_COMP                 &  !<-- ATM component
                           ,config  =CF                            &  !<-- configure object
                           ,rc      =RC)
!
      CALL ESMF_ConfigGetAttribute(config=CF                       &  !<-- config object
                                  ,value =DFIHR                    &  !<-- DFI filter duration
                                  ,label ='nhours_dfini:'          &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
!
      RESTRTTIME = STARTTIME
!
! Get FHROT from the config file
!
      MESSAGE_CHECK = "GFS get FHROT from CF"
      CALL ESMF_ConfigGetAttribute(config=CF, value=FHROT, label='fhrot:',rc=RC)
      CALL ERR_MSG(RC, MESSAGE_CHECK, RC_RUN)

!     write(0,*)' in gridcomp fhrot=',fhrot
      if (fhrot > 0) then
        CALL ESMF_TimeIntervalSet(timeinterval=FHROTINTVAL, h=FHROT, rc=RC)
        RESTRTTIME = STARTTIME + FHROTINTVAL
      endif
!
!------------
!
      MESSAGE_CHECK = "GFS get DFILEVS from CF"
      CALL ESMF_ConfigGetAttribute(config=CF                       &  !<-- config object
                                  ,value =DFILEVS                  &  !<-- vertical lvl above which no change from digital filter, default: levs, wam:80
                                  ,label ='dfilevs:'               &
                                  ,rc    =RC)
      CALL ERR_MSG(RC, MESSAGE_CHECK, RC_RUN)
!
!------------
!
      MESSAGE_CHECK = "GFS get LDFIFLTO from CF"
      CALL ESMF_ConfigGetAttribute(config=CF                       &  !<-- config object
                                  ,value =LDFIFLTO                 &  !<-- OUPUT filtered fields for fh>=dfihr
                                  ,label ='ldfiflto:'              &  !<-- Give this label's value to the previous var iable
                                  ,rc    =RC)
      CALL ERR_MSG(RC, MESSAGE_CHECK, RC_RUN)
!
!------------
!
      MESSAGE_CHECK = "GFS get LDFI_GRD from CF"
      CALL ESMF_ConfigGetAttribute(config=CF                       &  !<-- The config object
                                  ,value =LDFI_GRD                 &  !<-- grid point digital filter
                                  ,label ='ldfi_grd:'              &  !<-- Give this label's value to the previous var iable
                                  ,rc    =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
!
!------------
!
      MESSAGE_CHECK = "GFS get LWRTGRDCMP from CF"
      CALL ESMF_ConfigGetAttribute(config=CF                       &  !<-- config object
                                  ,value =LWRTGRDCMP               &  !<-- run wrt grd comp
                                  ,label ='lwrtgrdcmp:'            &  !<-- Give this label's value
                                  ,rc    =RC)
      CALL ERR_MSG(RC, MESSAGE_CHECK, RC_RUN)
!
!-----------------------------------------------------------------------
!***  Execute the GFS forecast runstream.
!-----------------------------------------------------------------------
!
      CALL GFS_INTEGRATE(gfs_int_state%GC_GFS_DYN                  &
                        ,gfs_int_state%GC_GFS_PHY                  &
                        ,gfs_int_state%GC_GFS_CHEM                 &
                        ,gfs_int_state%GC_GFS_CPL                  &
                        ,gfs_int_state%GC_PHY2CHEM_CPL             &
                        ,gfs_int_state%GC_CHEM2PHY_CPL             &
                        ,gfs_int_state%WRT_COMPS                   &
                        ,gfs_int_state%IMP_GFS_DYN                 &
                        ,gfs_int_state%EXP_GFS_DYN                 &
                        ,gfs_int_state%IMP_GFS_PHY                 &
                        ,gfs_int_state%EXP_GFS_PHY                 &
                        ,gfs_int_state%IMP_GFS_CHEM                &
                        ,gfs_int_state%EXP_GFS_CHEM                &
                        ,gfs_int_state%IMP_GFS_WRT                 &
                        ,gfs_int_state%EXP_GFS_WRT                 &
                        ,gfs_int_state%CLOCK_GFS                   &
                        ,gfs_int_state%TIMEINTERVAL_GFS_OUTPUT     &
                        ,gfs_int_state%TIMEINTERVAL_GFS_OUTPUT_HF  &
                        ,gfs_int_state%QUILTING                    &
                        ,gfs_int_state%WRITE_GROUP_READY_TO_GO     &
                        ,CURRTIME                                  &
                        ,STARTTIME                                 &
                        ,RESTRTTIME                                &
                        ,NTIMESTEP                                 &
                        ,TIMESTEP                                  &
                        ,NFHMAX_HF                                 &
                        ,DFIHR                                     &
                        ,DFILEVS                                   &
                        ,LDFI_GRD                                  &
                        ,LDFIFLTO                                  &
                        ,LWRTGRDCMP                                &
                        ,gfs_int_state%MYPE                        &
                        ,gfs_int_state%PHYSICS_ON                  &
                        ,gfs_int_state%CHEMISTRY_ON, VM_local)
!
!
!-----------------------------------------------------------------------
!***  WRITE THE FINAL ERROR SIGNAL.
!-----------------------------------------------------------------------
!
      IF(RC_RUN == ESMF_SUCCESS) THEN
          WRITE(0, *) 'GFS ATM RUN STEP SUCCEEDED'
      ELSE
          WRITE(0, *) 'GFS ATM RUN STEP FAILED RC_RUN=', RC_RUN
      ENDIF
!
!-----------------------------------------------------------------------
!
      rtim1           = rtc()
      total_integ_tim = total_integ_tim+(rtc()-btim0)
!     write(0,*)'exit ATM_RUN integration time ',total_integ_tim,(rtim1-rtim0)
!
!-----------------------------------------------------------------------
!
      call ESMF_ClockPrint(CLOCK_ATM, options="currTime", &
        preString="exiting GFS_RUN with CLOCK current: ", &
        unit=nuopcMsg)
      call ESMF_LogWrite(nuopcMsg, ESMF_LOGMSG_INFO)
      call ESMF_ClockPrint(CLOCK_ATM, options="startTime", &
        preString="exiting GFS_RUN with CLOCK start:   ", &
        unit=nuopcMsg)
      call ESMF_LogWrite(nuopcMsg, ESMF_LOGMSG_INFO)
      call ESMF_ClockPrint(CLOCK_ATM, options="stopTime", &
        preString="exiting GFS_RUN with CLOCK stop:    ", &
        unit=nuopcMsg)
      call ESMF_LogWrite(nuopcMsg, ESMF_LOGMSG_INFO)
      END SUBROUTINE GFS_RUN
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE GFS_FINALIZE(GFS_GRID_COMP, IMP_STATE, EXP_STATE      &
                             ,CLOCK_ATM, RC_FINAL)
!
!-----------------------------------------------------------------------
!***  THIS ROUTINE FINALIZES THE ATM GRIDDED COMPONENT.
!-----------------------------------------------------------------------
!
      TYPE(ESMF_GridComp)               :: GFS_GRID_COMP                   !<-- ATM gridded component
      TYPE(ESMF_State)                  :: IMP_STATE                   &   !<-- ATM finalize step's import state
                                         , EXP_STATE                       !<-- ATM finalize step's export state
      TYPE(ESMF_Clock)                  :: CLOCK_ATM                       !<-- main ESMF Clock
      TYPE(ESMF_Config)                 :: CF                              !<-- config object
      INTEGER,            INTENT(OUT)   :: RC_FINAL                        !<-- Return code for the Finalize step
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      INTEGER :: I,J,rc
!
      CHARACTER(50):: MODE
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC         = ESMF_SUCCESS
      RC_FINAL   = ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  For the moment, use a direct alias to the ATM Clock.
!-----------------------------------------------------------------------
!
      gfs_int_state%CLOCK_GFS = CLOCK_ATM

!-----------------------------------------------------------------------
!***  Call into GFS Run() method one last time to handle end_step output
!-----------------------------------------------------------------------

      CALL ESMF_GridCompRun(gfs_int_state%GC_GFS_DYN              &
                           ,importstate=gfs_int_state%imp_gfs_dyn &
                           ,exportstate=gfs_int_state%exp_gfs_dyn &
                           ,clock      =gfs_int_state%CLOCK_GFS   &
                           ,rc         =RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC, MESSAGE_CHECK, RC_FINAL)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      if (LWRTGRDCMP) then
        CALL WRITE_ASYNC_GFS(gfs_int_state%WRT_COMPs              &
                            ,gfs_int_state%exp_gfs_dyn            &
                            ,gfs_int_state%imp_gfs_wrt            &
                            ,gfs_int_state%exp_gfs_wrt            &
                            ,gfs_int_state%CLOCK_GFS              &
                            ,gfs_int_state%MYPE                   &
                            ,WRITE_GROUP_READY_TO_GO=             &
                             gfs_int_state%WRITE_GROUP_READY_TO_GO)
      endif

!-----------------------------------------------------------------------
!***  RETRIEVE THE CONFIG OBJECT CF FROM THE GFS GRIDDED COMPONENT.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK = "Retrieve Config Object from GFS Component"
!      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompGet(gridcomp=GFS_GRID_COMP                      &  !<-- GFS gridded component
                           ,config  =CF                                 &  !<-- config object (~namelist)
                           ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC, MESSAGE_CHECK, RC_FINAL)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!---------------------------------------------------------
!***  condition to run only adiabatic (dynamics only)
!---------------------------------------------------------
!
      CALL ESMF_ConfigGetAttribute(CF, value=MODE, label='adiabatic:', rc=RC)
!
      IF(TRIM(MODE) == '.true.') THEN
        gfs_int_state%PHYSICS_ON = ESMF_False
        write(0,*)' Finalize without physics coupling. '
       ELSE
        gfs_int_state%PHYSICS_ON = ESMF_True
        write(0,*)' Finalize with physics coupling. '
      ENDIF

!-----------------------------------------------------------------------

!
!-----------------------------------------------------------------------
!***  FINALIZE EACH OF THE SUBCOMPONENTS.
!-----------------------------------------------------------------------
!
!--------------
!***  DYNAMICS
!--------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK = "Finalize Dynamics Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompFinalize(gfs_int_state%GC_GFS_DYN             &
                               ,importstate=gfs_int_state%imp_gfs_dyn &
                               ,exportstate=gfs_int_state%exp_gfs_dyn &
                               ,clock      =gfs_int_state%CLOCK_GFS   &
                               ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!ratko    CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL)
! - FIX later
      RC=ESMF_SUCCESS
      RC_FINAL=ESMF_SUCCESS
!ratko
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF(gfs_int_state%PHYSICS_ON == ESMF_True) THEN
!
!-------------
!***  PHYSICS
!-------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK = "Finalize Physics Component"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_GridCompFinalize(gfs_int_state%GC_GFS_PHY              &
                                  ,importstate=gfs_int_state%IMP_GFS_PHY &
                                  ,exportstate=gfs_int_state%EXP_GFS_PHY &
                                  ,clock      =gfs_int_state%CLOCK_GFS   &
                                  ,rc         =RC)

      ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!ratko    CALL ERR_MSG(RC, MESSAGE_CHECK, RC_FINAL)
! - FIX later
        RC       = ESMF_SUCCESS
        RC_FINAL = ESMF_SUCCESS
!ratko
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------
!***  DYNAMICS-PHYSICS COUPLER
!-----------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK = "Finalize Dynamics-Physics Coupler"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_CplCompFinalize(gfs_int_state%GC_GFS_CPL               &
                                 ,importstate=gfs_int_state%EXP_GFS_DYN  &
                                 ,exportstate=gfs_int_state%IMP_GFS_PHY  &
                                 ,clock      =gfs_int_state%CLOCK_GFS    &
                                 ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC, MESSAGE_CHECK, RC_FINAL)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!     ENDIF 

!
!-----------------------------------------------------------------------
!***  Remove nested STATES
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK = "Remove nested States"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!

      CALL ESMF_StateRemove(IMP_STATE, itemNameList=(/"GFS dynamics import"/), &
                            rc=RC)
      CALL ERR_MSG(RC, MESSAGE_CHECK, RC_FINAL)
      CALL ESMF_StateRemove(EXP_STATE, itemNameList=(/"GFS dynamics export"/), &
                            rc=RC)
      CALL ERR_MSG(RC, MESSAGE_CHECK, RC_FINAL)
!
!-----------------------------------------------------------------------
!***  DESTROY ALL STATES.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK = "Destroy States"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateDestroy(gfs_int_state%IMP_GFS_DYN, rc=RC)
      CALL ESMF_StateDestroy(gfs_int_state%EXP_GFS_DYN, rc=RC)
      IF(gfs_int_state%PHYSICS_ON == ESMF_True) THEN
        CALL ESMF_StateDestroy(state=gfs_int_state%IMP_GFS_PHY, rc=RC)
        CALL ESMF_StateDestroy(state=gfs_int_state%EXP_GFS_PHY, rc=RC)
      ENDIF
!
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC, MESSAGE_CHECK, RC_FINAL)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  IF QUILTING WAS SELECTED FOR THE GENERATION OF OUTPUT,
!***  FINALIZE AND DESTROY OBJECTS RELATED TO THE WRITE COMPONENTS.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      IF(gfs_int_state%QUILTING) THEN
        CALL WRITE_DESTROY_GFS(GFS_GRID_COMP,                 &
                               gfs_int_state%WRT_COMPS,       &
                               gfs_int_state%IMP_GFS_WRT,     &
                               gfs_int_state%EXP_GFS_WRT,     &
                               gfs_int_state%CLOCK_GFS) 
      ENDIF
!
!-----------------------------------------------------------------------
!***  DESTROY ALL SUBCOMPONENTS.
!-----------------------------------------------------------------------
!
!--------------
!***  DYNAMICS
!--------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK = "Destroy Dynamics Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompDestroy(gfs_int_state%GC_GFS_DYN, rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC, MESSAGE_CHECK, RC_FINAL)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF (gfs_int_state%PHYSICS_ON == ESMF_True) THEN
!
!-------------
!***  PHYSICS
!-------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK = "Destroy Physics Component"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO, rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_GridCompDestroy(gfs_int_state%GC_GFS_PHY, rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC, MESSAGE_CHECK, RC_FINAL)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!------------------------------
!***  DYNAMICS-PHYSICS COUPLER
!------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK = "Destroy Dynamics-Physics Coupler"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_CplCompDestroy(gfs_int_state%GC_GFS_CPL, rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC, MESSAGE_CHECK, RC_FINAL)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      ENDIF
!
!-----------------------------------------------------------------------
!***  THE FINAL ERROR SIGNAL INFORMATION.
!-----------------------------------------------------------------------
!
      IF(RC_FINAL == ESMF_SUCCESS) THEN
        WRITE(0,*)'GFS FINALIZE STEP SUCCEEDED'
      ELSE
        WRITE(0,*)'GFS FINALIZE STEP FAILED'
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE GFS_FINALIZE
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      END MODULE MODULE_GFS_GRID_COMP
!
!-----------------------------------------------------------------------
