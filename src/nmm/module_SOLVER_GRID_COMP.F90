!-----------------------------------------------------------------------
!
      MODULE module_SOLVER_GRID_COMP
!
!-----------------------------------------------------------------------
!
!***  This module holds the Solver component's  Register, Init, Run, 
!***  and Finalize routines.  They are called from the DOMAIN component
!***  (DOMAIN_INITIALIZE calls SOLVER_INITIALIZE, etc.) 
!***  in MODULE_DOMAIN_GRID_COMP.F90.
!
!-----------------------------------------------------------------------
! HISTORY LOG:
!
!   2008-07-30  Janjic - Add CONVECTION='none' to OPERATIONAL_PHYSICS.
!               Janjic - Fix lower J limit in FFTFHN(WATER).
!   2008-08-23  Janjic - General pressure-sigma hybrid
!               Janjic - Consistent nonhydrostatic correction in the
!                        first term of the pressure gradient force
!   2008-09-03  Black  - Added initialization of boundary arrays
!                        for nests.
!   2009-03-12  Black  - Changes for general hybrid coordinate.
!   2009-11     Jovic  - Modified for ownership/import/export specification
!   2010-11-03  Pyle   - Modifications/corrections for digital filter.
!   2011-02     Yang   - Updated to use both the ESMF 4.0.0rp2 library,
!                        ESMF 5 series library and the the
!                        ESMF 3.1.0rp2 library.
!   2011-05-12  Yang   - Modified for using the ESMF 5.2.0r_beta_snapshot_07.
!   2011-12-22  Jovic  - Combined Dyn and Phy into single component.
!
!   2012-02-08  Yang   - Modified for using the ESMF 5.2.0rp1 library.
!   2012-04-06  Juang  - add passing argument for gbphys for idea
!   2012-07-20  Black  - Modified for generational usage.
!   2013-09-09 Moorthi - Adding SR, DTDT, and TRIGGERPERTS for GBPHYS call
!   2013-11-09 Xingren Wu - Adding DUSFCI/DVSFCI for GBPHYS call
!   2014-03-28 Xingren Wu - Add "_CPL" field for GBPHYS call
!   2014-05-14 J. Wang - Adding cgwf,prslrd0 and levr to gbphys call
!   2014-06-26 Weiguo Wang -- Add HURRICANE PBL and SFCLAY calls
!   2016-02-16 J. Wang - change newsas/sashal from logical to integer
!   2016-05-09 Ferrier/Janjic - Constants epsq2,epsl function of level,
!                       added subroutine TQadjust
!   2016-08-29 Weiguo wang -- add scale-aware convection schemes
!-----------------------------------------------------------------------
!
      USE MPI
      USE ESMF
      USE MODULE_KINDS
      USE MODULE_VARS,ONLY : FIND_VAR_INDX
      USE MODULE_VARS_STATE
      USE MODULE_SOLVER_INTERNAL_STATE                                     !<-- Horizontal loop limits obtained here
!
      USE MODULE_MY_DOMAIN_SPECS, IDS_share=>IDS,IDE_share=>IDE         &
                                 ,IMS_share=>IMS,IME_share=>IME         &
                                 ,ITS_share=>ITS,ITE_share=>ITE         &
                                 ,JDS_share=>JDS,JDE_share=>JDE         &
                                 ,JMS_share=>JMS,JME_share=>JME         &
                                 ,JTS_share=>JTS,JTE_share=>JTE 
!
      USE MODULE_EXCHANGE,ONLY: HALO_EXCH
!
      USE MODULE_GET_CONFIG
!
      USE MODULE_DERIVED_TYPES,ONLY : BC_H_ALL,BC_V_ALL
!
      USE MODULE_CONTROL,ONLY : NUM_DOMAINS_MAX,TIMEF,NMMB_FINALIZE
!
      USE MODULE_CONSTANTS,ONLY : A2,A3,A4,CAPPA,CP,ELIV,ELWV,EPSQ,G &
                                 ,P608,PQ0,R_D,TIW,DBZmin
!
      USE MODULE_DIAGNOSE,ONLY : EXIT,FIELD_STATS                       &
                                ,MAX_FIELDS,MAX_FIELDS_HR,MAX_FIELDS_W6 &
                                ,MAX_FIELDS_THO                         &
                                ,HMAXMIN,TWR,VMAXMIN,VWR,WRT_PCP        &
                                ,LAT_LON_BNDS
!
      USE MODULE_CLOCKTIMES,ONLY : INTEGRATION_TIMERS,TIMERS
!
      USE MODULE_ERROR_MSG,ONLY: ERR_MSG,MESSAGE_CHECK
!
      USE MODULE_FLTBNDS,ONLY : POLEHN,POLEWN,SWAPHN,SWAPWN
!
      USE MODULE_VARS,ONLY : VAR
!
      USE MODULE_NESTING,ONLY : SUFFIX_NESTBC
!
      USE MODULE_RADIATION  ,ONLY : RADIATION
      USE MODULE_RA_GFDL    ,ONLY : GFDL_INIT,RDTEMP,TIME_MEASURE
      USE MODULE_RA_RRTM    ,ONLY : RRTM_INIT
      USE MODULE_TURBULENCE
      USE MODULE_SF_JSFC    ,ONLY : JSFC_INIT
      USE MODULE_SF_GFDL    ,ONLY : JSFC_INIT4GFDL
      USE MODULE_BL_MYJPBL  ,ONLY : MYJPBL_INIT
      USE MODULE_LS_NOAHLSM ,ONLY : DZSOIL,NOAH_LSM_INIT                &
                                   ,NUM_SOIL_LAYERS,SLDPTH
      USE MODULE_CU_BMJ     ,ONLY : BMJ_INIT
      USE MODULE_CU_SAS     ,ONLY : SAS_INIT
      USE MODULE_CU_SASHUR  ,ONLY : SASHUR_INIT
      USE MODULE_CU_SCALE   ,ONLY : SCALECU_INIT
      USE MODULE_CONVECTION

      USE MODULE_MICROPHYSICS_NMM ,ONLY : GSMDRIVE                      &
                                         ,MICRO_RESTART
      USE MODULE_MP_ETANEW    ,ONLY : FERRIER_INIT
      USE MODULE_MP_FER_HIRES ,ONLY : FERRIER_INIT_HR
      USE MODULE_MP_WSM6      ,ONLY : WSM6INIT
      USE MODULE_MP_THOMPSON, ONLY  : thompson_init
      USE MODULE_MP_GFS       ,ONLY : GFSMP_INIT

      USE MODULE_H_TO_V ,ONLY : H_TO_V,H_TO_V_TEND
      USE MODULE_GWD    ,ONLY : GWD_INIT
      USE MODULE_PRECIP_ADJUST
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      PRIVATE
!
      PUBLIC :: SOLVER_REGISTER                                            
!
      INTEGER(kind=KINT),PUBLIC :: IM,JM,LM,RESTVAL
!
      INTEGER(kind=KINT) :: START_YEAR,START_MONTH,START_DAY            &
                           ,START_HOUR,START_MINUTE,START_SECOND
!
      INTEGER(kind=KINT),SAVE :: JC
!
      INTEGER(kind=KINT) :: NUM_PES
!
      LOGICAL(kind=KLOG),SAVE :: QUILTING                                  !<-- Was quilting specified by the user?
!
      LOGICAL(kind=KLOG) :: I_AM_A_NEST                                    !<-- Flag indicating if DOMAIN Component is a nest
!
      LOGICAL(kind=KLOG),SAVE :: MOVE_NOW                                  !<-- Flag indicating if nested moves this timestep
!     LOGICAL(kind=KLOG) :: MOVE_NOW                                    &  !<-- Flag indicating if nested moves this timestep
!                          ,MY_DOMAIN_MOVES                                !<-- Flag indicating if nested domain moves
 
      REAL(kind=KFPT),SAVE :: PT
!
      TYPE(SOLVER_INTERNAL_STATE),POINTER :: INT_STATE                     !<-- The Solver component internal state pointer.
!
!-----------------------------------------------------------------------
!***  For determining clocktimes of various pieces of the Solver.
!-----------------------------------------------------------------------
!
      REAL(kind=KDBL) :: btim,btim0
!
      TYPE(INTEGRATION_TIMERS),POINTER :: TD
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE SOLVER_REGISTER(GRID_COMP,RC_REG)
!
!-----------------------------------------------------------------------
!***  Register the Solver component's Initialize, Run, and Finalize
!***  subroutine names.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument variables
!------------------------
!
      TYPE(ESMF_GridComp) :: GRID_COMP                                    !<-- The Solver Gridded Component
!
      INTEGER(kind=KINT),INTENT(OUT) :: RC_REG                            !<-- Return code for Solver register
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: RC
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC    =ESMF_SUCCESS
      RC_REG=ESMF_SUCCESS                                                 !<-- Initialize error signal variable
!
!-----------------------------------------------------------------------
!***  Register the Solver initialize subroutine.  Since it is just one
!***  subroutine, use ESMF_SINGLEPHASE.  The second argument is
!***  a pre-defined subroutine type, such as ESMF_SETINIT, ESMF_SETRUN,
!***  or ESMF_SETFINAL.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set Entry Point for Solver Initialize"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompSetEntryPoint(GRID_COMP                         &  !<-- The gridded component
                                     ,ESMF_METHOD_INITIALIZE            &  !<-- Predefined subroutine type
                                     ,SOLVER_INITIALIZE                 &  !<-- User's subroutineName
                                     ,phase=1                           &
                                     ,rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Register the Solver Run subroutine.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set Entry Point for Solver Run"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompSetEntryPoint(GRID_COMP                         &  !<-- gridcomp
                                     ,ESMF_METHOD_RUN                   &  !<-- subroutineType
                                     ,SOLVER_RUN                        &  !<-- user's subroutineName
                                     ,phase=1                           &
                                     ,rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Register the Solver Finalize subroutine.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set Entry Point for Solver Finalize"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompSetEntryPoint(GRID_COMP                         &  !<-- gridcomp
                                     ,ESMF_METHOD_FINALIZE              &  !<-- subroutineType
                                     ,SOLVER_FINALIZE                   &  !<-- user's subroutineName
                                     ,phase=1                           &
                                     ,rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Check the error signal variable.
!-----------------------------------------------------------------------
!
      IF(RC_REG==ESMF_SUCCESS)THEN
!       WRITE(0,*)" SOLVER_REGISTER SUCCEEDED"
      ELSE
        WRITE(0,*)" SOLVER_REGISTER FAILED"
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE SOLVER_REGISTER
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE SOLVER_INITIALIZE (GRID_COMP                           &
                                   ,IMP_STATE                           &
                                   ,EXP_STATE                           &
                                   ,CLOCK_ATM                           &
                                   ,RC_INIT)
!
!-----------------------------------------------------------------------
!***  Carry out all necessary setups for the model Solver.
!-----------------------------------------------------------------------
!
      USE MODULE_CONTROL,ONLY : CONSTS
!
      USE MODULE_INIT_READ_BIN,ONLY : READ_BINARY
      USE MODULE_INIT_READ_NEMSIO,ONLY : READ_NEMSIO
!
      USE MODULE_FLTBNDS,ONLY : PREFFT, PRESMUD
      USE MODULE_TRACKER

!------------------------
!***  Argument variables
!------------------------
!
      TYPE(ESMF_GridComp) :: GRID_COMP                                     !<-- The Solver gridded component
!
      TYPE(ESMF_State) :: IMP_STATE                                     &  !<-- The Solver Initialize step's import state
                         ,EXP_STATE                                        !<-- The Solver Initialize step's export state
!
      TYPE(ESMF_Clock) :: CLOCK_ATM                                        !<-- The ATM's ESMF Clock
!
      INTEGER,INTENT(OUT) :: RC_INIT
!
!---------------------
!***  Local variables
!---------------------
!
      INTEGER(kind=KINT),SAVE :: N8=8
!
      INTEGER(kind=KINT) :: IDE,IDS,IME,IMS,ITE,ITS                     &
                           ,JDE,JDS,JME,JMS,JTE,JTS
!
      INTEGER(kind=KINT) :: IHALO,JHALO,MPI_COMM_COMP,MY_DOMAIN_ID      &
                           ,MY_DOMAIN_ID_LOC,MYPE,NUM_PES,UBOUND_VARS
!
      INTEGER(kind=KINT) :: I,I_INC,IDENOMINATOR_DT                     &
                           ,IEND,IERR,INTEGER_DT                        &
                           ,J,J_INC,JEND,KK,KOUNT,KSE,KSS               &
                           ,L,LL,LMP1,LNSH,LNSV                         &
                           ,N,NUMERATOR_DT,NV,RC
!
      INTEGER(kind=KINT) :: ITE_H2,ITS_H2,JTE_H2,JTS_H2
!
      INTEGER(kind=KINT),DIMENSION(1:8) :: MY_NEB
!
      REAL(kind=KFPT) :: DPH,DLM,DT,GLATX,GLONX,SB_1,SBD_1,TLATX,TLONX  &
                        ,TPH0_1,TPH0D_1,TLM0_1,TLM0D_1,WB_1,WBD_1       &
                        ,X,Y,Z
!
      REAL(kind=KFPT),DIMENSION(1:2) :: SW_X
!
      REAL(kind=DOUBLE) :: D2R,D_ONE,D_180,PI
!
      LOGICAL(kind=KLOG) :: RUN_LOCAL
!
      CHARACTER(20) :: FIELD_NAME
!
      TYPE(WRAP_SOLVER_INT_STATE) :: WRAP                                  ! <-- This wrap is a derived type which contains
                                                                           !     only a pointer to the internal state.  It is needed
                                                                           !     for using different architectures or compilers.
!
      TYPE(ESMF_Grid) :: GRID                                              !<-- The ESMF Grid
!
      TYPE(ESMF_VM) :: VM                                                  !<-- The ESMF Virtual Machine
!
      TYPE(ESMF_Field) :: FIELD
!
      TYPE(ESMF_FieldBundle) :: BUNDLE_NESTBC
!
      TYPE(ESMF_TimeInterval) :: DT_ESMF                                   !<-- The ESMF fundamental timestep (s)
!
      TYPE(ESMF_Config) :: CF                                              !<-- ESMF configure object
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      btim0=timef()
!
!-----------------------------------------------------------------------
!***  Initialize the error signal variables.
!-----------------------------------------------------------------------
!
      RC     =ESMF_SUCCESS
      RC_INIT=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  Allocate the Solver internal state pointer.
!-----------------------------------------------------------------------
!
      ALLOCATE(INT_STATE,STAT=RC)
!
!-----------------------------------------------------------------------
!***  Attach the internal state to the Solver gridded component.
!-----------------------------------------------------------------------
!
      WRAP%INT_STATE=>INT_STATE
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Attach Solver Internal State to the Gridded Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompSetInternalState(GRID_COMP                      &  !<-- The Solver gridded component
                                        ,WRAP                           &  !<-- Pointer to the Solver internal state
                                        ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Retrieve fundamental domain characteristics from the Solver   
!***  import state and set them in the internal state so they will
!***  always be available to this component.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Get Domain Dimensions from Solver Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=IMP_STATE                            &  !<-- The Solver import state
                            ,name ='ITS'                                &  !<-- Name of variable to get from Solver import state
                            ,value=int_state%ITS                        &  !<-- Put extracted value here
                            ,rc   =RC)
!
      CALL ESMF_AttributeGet(state=IMP_STATE                            &  !<-- The Solver import state
                            ,name ='ITE'                                &  !<-- Name of variable to get from Solver import state
                            ,value=int_state%ITE                        &  !<-- Put extracted value here
                            ,rc   =RC)
!
      CALL ESMF_AttributeGet(state=IMP_STATE                            &  !<-- The Solver import state
                            ,name ='JTS'                                &  !<-- Name of variable to get from Solver import state
                            ,value=int_state%JTS                        &  !<-- Put extracted value here
                            ,rc   =RC)
!
      CALL ESMF_AttributeGet(state=IMP_STATE                            &  !<-- The Solver import state
                            ,name ='JTE'                                &  !<-- Name of variable to get from Solver import state
                            ,value=int_state%JTE                        &  !<-- Put extracted value here
                            ,rc   =RC)
!
      CALL ESMF_AttributeGet(state=IMP_STATE                            &  !<-- The Solver import state
                            ,name ='IMS'                                &  !<-- Name of variable to get from Solver import state
                            ,value=int_state%IMS                        &  !<-- Put extracted value here
                            ,rc   =RC)
!
      CALL ESMF_AttributeGet(state=IMP_STATE                            &  !<-- The Solver import state
                            ,name ='IME'                                &  !<-- Name of variable to get from Solver import state
                            ,value=int_state%IME                        &  !<-- Put extracted value here
                            ,rc   =RC)
!
      CALL ESMF_AttributeGet(state=IMP_STATE                            &  !<-- The Solver import state
                            ,name ='JMS'                                &  !<-- Name of variable to get from Solver import state
                            ,value=int_state%JMS                        &  !<-- Put extracted value here
                            ,rc   =RC)
!
      CALL ESMF_AttributeGet(state=IMP_STATE                            &  !<-- The Solver import state
                            ,name ='JME'                                &  !<-- Name of variable to get from Solver import state
                            ,value=int_state%JME                        &  !<-- Put extracted value here
                            ,rc   =RC)
!
      CALL ESMF_AttributeGet(state=IMP_STATE                            &  !<-- The Solver import state
                            ,name ='IDS'                                &  !<-- Name of variable to get from Solver import state
                            ,value=int_state%IDS                        &  !<-- Put extracted value here
                            ,rc   =RC)
!
      CALL ESMF_AttributeGet(state=IMP_STATE                            &  !<-- The Solver import state
                            ,name ='IDE'                                &  !<-- Name of variable to get from Solver import state
                            ,value=int_state%IDE                        &  !<-- Put extracted value here
                            ,rc   =RC)
!
      CALL ESMF_AttributeGet(state=IMP_STATE                            &  !<-- The Solver import state
                            ,name ='JDS'                                &  !<-- Name of variable to get from Solver import state
                            ,value=int_state%JDS                        &  !<-- Put extracted value here
                            ,rc   =RC)
!
      CALL ESMF_AttributeGet(state=IMP_STATE                            &  !<-- The Solver import state
                            ,name ='JDE'                                &  !<-- Name of variable to get from Solver import state
                            ,value=int_state%JDE                        &  !<-- Put extracted value here
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Get Halo Widths from Solver Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=IMP_STATE                            &  !<-- The Solver import state
                            ,name ='IHALO'                              &  !<-- Name of variable to get from Solver import state
                            ,value=int_state%IHALO                      &  !<-- Put extracted value here
                            ,rc   =RC)
!
      CALL ESMF_AttributeGet(state=IMP_STATE                            &  !<-- The Solver import state
                            ,name ='JHALO'                              &  !<-- Name of variable to get from Solver import state
                            ,value=int_state%JHALO                      &  !<-- Put extracted value here
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Get Fcst/Quilt Task Intracomm from Solver Imp State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=IMP_STATE                            &  !<-- The Solver import state
                            ,name ='Fcst/Quilt Intracommunicators'      &  !<-- Name of variable to get from Solver import state
                            ,value=int_state%MPI_COMM_COMP              &  !<-- Put extracted value here
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Task Neighbors from Solver Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state    =IMP_STATE                      &  !<-- The Solver import state
                            ,name     ='MY_NEB'                       &  !<-- Name of the attribute to extract
                            ,valueList=int_state%MY_NEB               &  !<-- Insert Attribute into Solver internal state
                            ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Insert the local domain starting limits and the halo width into
!***  the Solver internal state.
!-----------------------------------------------------------------------
!
      ITS=int_state%ITS
      ITE=int_state%ITE
      IMS=int_state%IMS
      IME=int_state%IME
      IDS=int_state%IDS
      IDE=int_state%IDE
!
      JTS=int_state%JTS
      JTE=int_state%JTE
      JMS=int_state%JMS
      JME=int_state%JME
      JDS=int_state%JDS
      JDE=int_state%JDE
!
      int_state%ITS_B1=MAX(ITS,IDS+1)
      int_state%ITE_B1=MIN(ITE,IDE-1)
      int_state%ITS_B2=MAX(ITS,IDS+2)
      int_state%ITE_B2=MIN(ITE,IDE-2)
      int_state%ITS_B1_H1=MAX(ITS-1,IDS+1)
      int_state%ITE_B1_H1=MIN(ITE+1,IDE-1)
      int_state%ITE_B1_H2=MIN(ITE+2,IDE-1)
      int_state%ITS_H1=MAX(ITS-1,IDS)
      int_state%ITE_H1=MIN(ITE+1,IDE)
      int_state%ITS_H2=MAX(ITS-2,IDS)
      int_state%ITE_H2=MIN(ITE+2,IDE)
      int_state%JTS_B1=MAX(JTS,JDS+1)
      int_state%JTE_B1=MIN(JTE,JDE-1)
      int_state%JTS_B2=MAX(JTS,JDS+2)
      int_state%JTE_B2=MIN(JTE,JDE-2)
      int_state%JTS_B1_H1=MAX(JTS-1,JDS+1)
      int_state%JTE_B1_H1=MIN(JTE+1,JDE-1)
      int_state%JTE_B1_H2=MIN(JTE+2,JDE-1)
      int_state%JTS_H1=MAX(JTS-1,JDS)
      int_state%JTE_H1=MIN(JTE+1,JDE)
      int_state%JTS_H2=MAX(JTS-2,JDS)
      int_state%JTE_H2=MIN(JTE+2,JDE)
!
      IHALO=int_state%IHALO
      JHALO=int_state%JHALO
!
      ! Disable the tracker by default.  This may be overridden below
      ! when reading the configure file.
      int_state%NTRACK_trigger=0
!
      IF(IHALO==JHALO)THEN
        int_state%NHALO=IHALO
      ELSE
        RC_INIT=ESMF_FAILURE
        WRITE(0,*)'Error due to ihalo /= jhalo'
      ENDIF
!
!-----------------------------------------------------------------------
!***  Use ESMF utilities to get information from the configuration file.
!***  The function is similar to reading a namelist.  The GET_CONFIG
!***  routine is the user's.  It extracts values from the config file
!***  and places them in the namelist components of the internal state.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Get Configure File Parameters for Solver"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL GET_CONFIG_DIMS (GRID_COMP                                   &
                           ,int_state%INPES,int_state%JNPES             &
                           ,LM                                          &
                           ,int_state%NUM_TRACERS_CHEM                  &
                           ,int_state%PCPHR                             &
                           ,int_state%GFS                               &
                           ,int_state%MICROPHYSICS                      &
                           ,int_state%SHORTWAVE                         &
                           ,int_state%LONGWAVE                          &
                           ,int_state%LMPRATE                           &
                           ,int_state%LNSH, int_state%LNSV              &
                           ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      LNSH=int_state%LNSH
      LNSV=int_state%LNSV
!
!-----------------------------------------------------------------------
!***  We must know whether or not this is a global domain.  Get the
!***  configure object from the Solver component and extract the
!***  value of 'global'.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Solver_Init: Retrieve Config Object from Solver Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompGet(gridcomp=GRID_COMP                          &   !<--- The Solver component
                           ,config  =CF                                 &   !<--- The configure (namelist) object
                           ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Solver_Init: Extract GLOBAL from Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure file object
                                  ,value =int_state%GLOBAL              &  !<-- Put extracted quantity here
                                  ,label ='global:'                     &  !<-- The quantity's label in the configure file
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
! Is the tracker enabled?  Triggers allocations later on.  (Or, rather,
! it would do that if such a thing was supported by the current
! framework.)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure file object
                                  ,value =int_state%NTRACK_TRIGGER      &  !<-- Put extracted quantity here
                                  ,label ='ntrack:'                     &  !<-- The quantity's label in the configure file
                                  ,rc    =RC)
      if(RC/=0) then
         print '(A)','Disabling tracker and tracker vars for domain.'
         int_state%NTRACK_TRIGGER=0
      endif

      int_state%HIFREQ_file=' '
      int_state%PATCF_file=' '
      if(int_state%NTRACK_TRIGGER /= 0) then ! Check for additional tracker options.
         ! Per-timestep output.
         call ESMF_ConfigGetAttribute(config=CF,value=int_state%HIFREQ_file,label='hifreq:',rc=RC)
         if(RC/=0) then
            print '(A)','Disabling per-timestep output because "hifreq:" was not specified.'
            int_state%HIFREQ_file=' '
         endif
         ! Per-tracker-step output.
         call ESMF_ConfigGetAttribute(config=CF,value=int_state%PATCF_file,label='patcf:',rc=RC)
         if(RC/=0) then
            print '(A)','Disabling tracker output because "patcf:" was not specified.'
            int_state%PATCF_file=' '
         endif
      endif
!
!-----------------------------------------------------------------------
!***  Retrieve the VM to obtain the task ID and total number of tasks
!***  for the internal state.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Get VM from the Solver Gridded Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompGet(gridcomp=GRID_COMP                          &  !<-- The Solver gridded component
                           ,vm      =VM                                 &  !<-- The ESMF Virtual Machine
                           ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Get Task IDs and Number of MPI Tasks from VM"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_VMGet(vm      =VM                                       &  !<-- The ESMF virtual machine
                     ,localpet=int_state%MYPE                           &  !<-- My task's local rank on this domain
                     ,petcount=int_state%NUM_PES                        &  !<-- Total number of MPI tasks
                     ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  int_state%NUM_PES taken from VM is the total number of tasks 
!***  on this domain including Write/Quilt tasks.  We want only the
!***  number of forecast tasks.
!-----------------------------------------------------------------------
!
      int_state%NUM_PES=int_state%INPES*int_state%JNPES
!
      NUM_PES=int_state%NUM_PES                                            !<-- The number of forecast tasks
      MYPE=int_state%MYPE                                                  !<-- The local task ID
!
!-----------------------------------------------------------------------
!***  Only forecast tasks are needed for the remaining
!***  initialization process.
!-----------------------------------------------------------------------
!
      fcst_tasks: IF(MYPE<NUM_PES)THEN                                     !<-- Select only forecast tasks
!
!-----------------------------------------------------------------------
!***  Allocate all necessary internal state variables.  Those that
!***  are owned/exported are pointed into allocated memory within
!***  the Solver's composite VARS array.  
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Solver_Init: Allocate internal state variables"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL SET_INTERNAL_STATE_SOLVER(INT_STATE                        &
                                      ,LM                               &
                                      ,ITS,ITE,JTS,JTE                  &
                                      ,IMS,IME,JMS,JME                  &
                                      ,IDS,IDE,JDS,JDE                  &
                                      ,IHALO,JHALO                      &
                                      ,MYPE                             &
                                      ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract the ESMF Grid from the Solver Component"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_GridCompGet(gridcomp=GRID_COMP                        &  !<-- The Solver gridded component
                             ,grid    =GRID                             &  !<-- The ESMF Grid
                             ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Put the allocated pointers of all export/import variables
!***  into the Solver export/import states.  
!-----------------------------------------------------------------------
!
        CALL PUT_VARS_IN_STATE(int_state%VARS,int_state%NUM_VARS,'X',GRID,EXP_STATE)
!
        CALL PUT_VARS_IN_STATE(int_state%VARS,int_state%NUM_VARS,'I',GRID,IMP_STATE)
!
!-----------------------------------------------------------------------
!
      ENDIF fcst_tasks
!
!-----------------------------------------------------------------------
!***  Use ESMF utilities to get information from the configuration file.
!***  The function is similar to reading a namelist.  The GET_CONFIG
!***  routine is the user's.  It extracts values from the config file
!***  and places them in the namelist components of the internal state.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Get Configure File Parameters for Solver"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL GET_CONFIG(GRID_COMP,INT_STATE,RC)                             !<-- User's routine to extract config file information
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Only forecast tasks are needed for the remaining
!***  initialization process.
!-----------------------------------------------------------------------
!
      fcst_tasks2: IF(int_state%MYPE<int_state%NUM_PES)THEN                !<-- Select only forecast tasks
!
!-----------------------------------------------------------------------
!***  Assign the fundamental timestep retrieved from the clock.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract Fundamental Timestep from ATM's Clock"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ClockGet(clock   =CLOCK_ATM                           &  !<-- The ATM Clock
                          ,timeStep=DT_ESMF                             &  !<-- Fundamental timestep (s) (ESMF)
                          ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Get Real Timestep from ESMF Timestep"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!

        CALL ESMF_TimeIntervalGet(timeinterval=DT_ESMF                  &  !<-- the ESMF timestep
                                 ,s           =INTEGER_DT               &  !<-- the integer part of the timestep in seconds
                                 ,sN          =NUMERATOR_DT             &  !<-- the numerator of the fractional second
                                 ,sD          =IDENOMINATOR_DT          &  !<-- the denominator of the fractional second
                                 ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        int_state%DT=REAL(INTEGER_DT)+REAL(NUMERATOR_DT)                &  !<-- Fundamental tiemstep (s) (REAL)
                                     /REAL(IDENOMINATOR_DT)
        DT=int_state%DT
!
        int_state%NSTEPS_PER_HOUR=NINT(3600./DT)
        int_state%NSTEPS_PER_RESET=NINT(int_state%AVGMAXLEN/DT)
        int_state%NSTEPS_PER_CHECK=MAX(2,NINT(40/DT))
!
!-----------------------------------------------------------------------
!***  Save fundamental timestep to distinguish from filter timestep
!***  which may be shorter
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Set Dyn Timestep to Distinguish from Filter DT"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=IMP_STATE                          &  !<-- The Solver import state
                              ,name ='FUND_DT'                          &  !<-- Name of variable to get from Solver import state
                              ,value=DT                                 &  !<-- Put extracted value here
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        int_state%FIRST_NMM=.TRUE.
!
        int_state%DT_LAST=0.                                               !<-- For use in digital filtering in SOLVE_RUN
        int_state%DT_TEST_RATIO=0.                                         !<-- For use in digital filtering in SOLVE_RUN
!
!-----------------------------------------------------------------------
!***  Retrieve the domain ID from the Solver import state.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Get Domain ID from Solver Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state=IMP_STATE                          &  !<-- The Solver import state
                              ,name ='DOMAIN_ID'                        &  !<-- Name of variable to get from Solver import state
                              ,value=MY_DOMAIN_ID_LOC                   &  !<-- Put extracted value here
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        int_state%MY_DOMAIN_ID=MY_DOMAIN_ID_LOC
!
        int_state%MY_DOMAIN_MOVES=.FALSE.
!
!-----------------------------------------------------------------------
!***  Was quilting specified by the user?
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Solver_Init: Was Quilting Specified?"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state=IMP_STATE                          &  !<-- The Solver import state
                              ,name ='Quilting'                         &  !<-- Name of variable to get from Solver import state
                              ,value=QUILTING                           &  !<-- Put extracted value here
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Initialize allocated arrays.
!-----------------------------------------------------------------------
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%PD(I,J)=0.
          int_state%PDO(I,J)=0.
        ENDDO
        ENDDO
!
        DO L=1,LM-1
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%PSGDT(I,J,L)=0.
        ENDDO
        ENDDO
        ENDDO
!
!-- Initialize 4D microphysics rates (diagnostic arrays)
!
        DO KK=1,int_state%d_ss
          DO L=1,LM
            DO J=JMS,JME
            DO I=IMS,IME
              int_state%MPRATES(I,J,L,KK)=0.
            ENDDO
            ENDDO
          ENDDO
        ENDDO
!
!-- Initialize all tracer-related arrays to zero. Water vapor mixing 
!   ratio array (int_state%QV) is no longer in solver_run but is
!   calculated when needed in the physics drivers.
!   
        DO N=1,int_state%NUM_TRACERS_TOTAL
          DO L=1,LM
            DO J=JMS,JME
            DO I=IMS,IME
              int_state%TRACERS     (I,J,L,N)=1.E-20
              int_state%TRACERS_SQRT(I,J,L,N)=1.E-20
              int_state%TRACERS_PREV(I,J,L,N)=1.E-20
              int_state%TRACERS_TEND(I,J,L,N)=1.E-20
            ENDDO
            ENDDO
          ENDDO
        ENDDO
!
!-- Initialize all "normal" 3D arrays
!
        DO L=1,LM
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%Told(I,J,L)=0.
          int_state%Tadj(I,J,L)=0.
          int_state%F_ICE(I,J,L)=0.
          int_state%F_RAIN(I,J,L)=0.
          int_state%F_RIMEF(I,J,L)=0.
          int_state%refl_10cm(I,J,L)=DBZmin
          int_state%Q2(I,J,L)=0.02       !=> int_state%TRACERS(:,:,:,INDX_Q2)
          int_state%OMGALF(I,J,L)=0.
          int_state%T(I,J,L)=-1.E6
          int_state%U(I,J,L)=-1.E6
          int_state%V(I,J,L)=-1.E6
          int_state%RLWTT(I,J,L)=0.
          int_state%RSWTT(I,J,L)=0.
          int_state%EXCH_H(I,J,L)=0.
          int_state%XLEN_MIX(I,J,L)=0.
          int_state%CLDFRA(I,J,L)=0.
          int_state%TRAIN(I,J,L) =0.
          int_state%TCUCN(I,J,L) =0.
          int_state%TCT(I,J,L) =-1.E6
          int_state%TCU(I,J,L) =-1.E6
          int_state%TCV(I,J,L) =-1.E6
          int_state%W_TOT(I,J,L)=0.
        ENDDO
        ENDDO
        ENDDO
!
        int_state%I_PAR_STA=0
        int_state%J_PAR_STA=0
        int_state%NMTS=-999
!
        DO L=1,NUM_SOIL_LAYERS
          int_state%SLDPTH(L)=SLDPTH(L)
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%SMC(I,J,L)=-1.E6
          int_state%STC(I,J,L)=-1.E6
          int_state%SH2O(I,J,L)=-1.E6
        ENDDO
        ENDDO
        ENDDO
!
        DO L=1,MICRO_RESTART
          int_state%MP_RESTART_STATE(L)=0.
          int_state%TBPVS_STATE(L)=0.
          int_state%TBPVS0_STATE(L)=0.
        ENDDO
        DO L=1, int_state%MDRMAXout-int_state%MDRMINout+1
           int_state%MASSRout(L)=0.
        ENDDO
        DO L=1, int_state%MDIMAXout-int_state%MDIMINout+1
           int_state%MASSIout(L)=0.
        ENDDO
!
        int_state%NSOIL=NUM_SOIL_LAYERS
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%LPBL(I,J)    =-999
          int_state%NCFRCV(I,J)  =-999
          int_state%NCFRST(I,J)  =-999
          int_state%ACFRCV(I,J)  =-1.E6
          int_state%ACFRST(I,J)  =-1.E6
          int_state%AKHS(I,J)    = 0.
          int_state%AKHS_OUT(I,J)= 0.
          int_state%AKMS(I,J)    = 0.
          int_state%AKMS_OUT(I,J)= 0.
          int_state%ALBASE(I,J)  =-1.E6
          int_state%ALBEDO(I,J)  =-1.E6
          int_state%ALWIN(I,J)   =-1.E6
          int_state%ALWOUT(I,J)  =-1.E6
          int_state%ALWTOA(I,J)  =-1.E6
          int_state%ASWIN(I,J)   =-1.E6
          int_state%ASWOUT(I,J)  =-1.E6
          int_state%ASWTOA(I,J)  =-1.E6
          int_state%BGROFF(I,J)  =-1.E6
          int_state%CFRACH(I,J)  =-1.E6
          int_state%CFRACM(I,J)  =-1.E6
          int_state%CFRACL(I,J)  =-1.E6
          int_state%CNVBOT(I,J)  =0.0
          int_state%CNVTOP(I,J)  =0.0
          int_state%CMC(I,J)     =-1.E6
          int_state%CPRATE(I,J)  =0.0
          int_state%CUPPT(I,J)   =-1.E6
          int_state%CZMEAN(I,J)  =-1.E6
          int_state%CZEN(I,J)    =-1.E6
          int_state%LSPA(I,J)    =-1.E6
          int_state%EPSR(I,J)    =-1.E6
          int_state%FIS(I,J)     =-1.E6
          int_state%HBOT(I,J)    =-1.E6
          int_state%HBOTD(I,J)   =-1.E6
          int_state%HBOTS(I,J)   =-1.E6
          int_state%HTOP(I,J)    =-1.E6
          int_state%HTOPD(I,J)   =-1.E6
          int_state%HTOPS(I,J)   =-1.E6
          int_state%GRNFLX(I,J)  = 0.
          int_state%MAVAIL(I,J)  = 1.
          int_state%MXSNAL(I,J)  =-1.E6
          int_state%PBLH(I,J)    =-1.E6
          int_state%MIXHT(I,J)   =0.
          int_state%PD(I,J)      =-1.E6
          int_state%POTEVP(I,J)  = 0.
          int_state%POTFLX(I,J)  =-1.E6
          int_state%QSH(I,J)     = 0.
          int_state%QWBS(I,J)    =-1.E6
          int_state%QZ0(I,J)     = 0.
          int_state%RADOT(I,J)   = 0.
          int_state%RLWIN(I,J)   = 0.
          int_state%RMOL(I,J)    =-1.E6
          int_state%RSWIN(I,J)   = 0.
          int_state%RSWINC(I,J)  = 0.
          int_state%RSWOUT(I,J)  = 0.
          int_state%RLWTOA(I,J)  = 0.
          int_state%RSWTOA(I,J)  = 0.
          int_state%SFCEVP(I,J)  = 0.
          int_state%SFCEXC(I,J)  = 0.
          int_state%SFCLHX(I,J)  =-1.E6
          int_state%SFCSHX(I,J)  =-1.E6
          int_state%SICE(I,J)    =-1.E6
          int_state%SIGT4(I,J)   =-1.E6
          int_state%SM(I,J)      =-1.E6
          int_state%SMSTAV(I,J)  = 0.
          int_state%SMSTOT(I,J)  = 0.
          int_state%SNO(I,J)     = 0.
          int_state%SNOWC(I,J)   = 0.
          int_state%SNOPCX(I,J)  =-1.E6
          int_state%SOILTB(I,J)  = 273.
          int_state%SR(I,J)      =-1.E6
          int_state%SSROFF(I,J)  = 0.
          int_state%SST(I,J)     = 273.
          int_state%SUBSHX(I,J)  =-1.E6
          int_state%THS(I,J)     =-1.E6
          int_state%THZ0(I,J)    = 273.
          int_state%TSKIN(I,J)   =-1.E6
          int_state%TWBS(I,J)    =-1.E6
          int_state%USTAR(I,J)   = 0.1
          int_state%UZ0(I,J)     = 0.
          int_state%VEGFRC(I,J)  =-1.E6
          int_state%VZ0(I,J)     = 0.
          int_state%Z0(I,J)      =-1.E6
          int_state%Z0BASE(I,J)  =-1.E6
          int_state%STDH(I,J)    =-1.E6
          int_state%CROT(I,J)    = 0.
          int_state%SROT(I,J)    = 0.
          int_state%HSTDV(I,J)   = 0.
          int_state%HCNVX(I,J)   = 0.
          int_state%HASYW(I,J)   = 0.
          int_state%HASYS(I,J)   = 0.
          int_state%HASYSW(I,J)  = 0.
          int_state%HASYNW(I,J)  = 0.
          int_state%HLENW(I,J)   = 0.
          int_state%HLENS(I,J)   = 0.
          int_state%HLENSW(I,J)  = 0.
          int_state%HLENNW(I,J)  = 0.
          int_state%HANGL(I,J)   = 0.
          int_state%HANIS(I,J)   = 0.
          int_state%HSLOP(I,J)   = 0.
          int_state%HZMAX(I,J)   = 0.
        ENDDO
        ENDDO
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%ACSNOM(I,J)= 0.
          int_state%ACSNOW(I,J)= 0.
          int_state%ACPREC(I,J)= 0.
          int_state%ACPREC_TOT(I,J)= 0.
          int_state%acpcp_ra(I,J)= 0.
          int_state%acpcp_sn(I,J)= 0.
          int_state%acpcp_gr(I,J)= 0.
          int_state%CUPREC(I,J)= 0.
          int_state%PREC(I,J)  = 0.
          int_state%CLDEFI(I,J)= 0.
          int_state%PSHLTR(I,J)= 1.E5
          int_state%P10(I,J)   = 1.E5
          int_state%PSFC(I,J)  = 1.E5
          int_state%Q02(I,J)   = 0.
          int_state%Q10(I,J)   = 0.
          int_state%QSHLTR(I,J)= 0.
          int_state%T2(I,J)    = 273.
          int_state%TH02(I,J)  = 0.
          int_state%TH10(I,J)  = 273.
          int_state%TSHLTR(I,J)= 273.
          int_state%U10(I,J)   = 0.
          int_state%V10(I,J)   = 0.
          int_state%TLMIN(I,J) = 0.
          int_state%TLMAX(I,J) = 0.

          int_state%ACUTIM(I,J)= 0.
          int_state%APHTIM(I,J)= 0.
          int_state%ARDLW(I,J) = 0.
          int_state%ARDSW(I,J) = 0.
          int_state%ASRFC(I,J) = 0.
          int_state%AVRAIN(I,J)= 0.
          int_state%AVCNVC(I,J)= 0.
        ENDDO
        ENDDO
!
        IF (int_state%has_reqc.eq.1 .and. int_state%has_reqi.eq.1 .and. int_state%has_reqs.eq.1) THEN
        DO L=1,LM
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%re_cloud(I,J,L)=2.51E-6
          int_state%re_ice(I,J,L)=10.1E-6
          int_state%re_snow(I,J,L)=20.1E-6
        ENDDO
        ENDDO
        ENDDO
        ENDIF
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%TLMAX(I,J)=-999.
          int_state%TLMIN(I,J)=999.
          int_state%T02MAX(I,J)=-999.
          int_state%T02MIN(I,J)=999.
          int_state%RH02MAX(I,J)=-999.
          int_state%RH02MIN(I,J)=999.
          int_state%SPD10MAX(I,J)=-999.
          int_state%UPHLMAX(I,J)=0.
          int_state%U10MAX(I,J)=-999.
          int_state%V10MAX(I,J)=-999.
          int_state%UPVVELMAX(I,J)=-999.
          int_state%DNVVELMAX(I,J)=999.
          int_state%T10AVG(I,J)=0.
          int_state%T10(I,J)=0.
          int_state%PSFCAVG(I,J)=0.
          int_state%AKHSAVG(I,J)=0.
          int_state%AKMSAVG(I,J)=0.
          int_state%SNOAVG(I,J)=0.
          int_state%REFDMAX(I,J)=DBZmin
          int_state%UPHLMAX(I,J)=-999.
        ENDDO
        ENDDO
        int_state%NCOUNT=0
!
        DO N=1,NUM_DOMAINS_MAX
          int_state%NTSCM(N)=-999
        ENDDO
!
        int_state%BDY_WAS_READ=.FALSE.
!
!! End of tracker variables
!###    Tracker scalar integer
!rv     int_state%NTRACK=0
        int_state%TRACK_HAVE_GUESS=0
        int_state%TRACK_N_OLD=0
        int_state%TRACKER_HAVEFIX=0
        int_state%TRACKER_GAVE_UP=0
!###    Tracker scalar real
        int_state%TRACK_LAST_HOUR=0.
        int_state%TRACK_GUESS_LAT=0.
        int_state%TRACK_GUESS_LON=0.
        int_state%TRACK_EDGE_DIST=0.
        int_state%TRACK_STDERR_M1=0.
        int_state%TRACK_STDERR_M2=0.
        int_state%TRACK_STDERR_M3=0.
        int_state%TRACKER_FIXLAT=0.
        int_state%TRACKER_FIXLON=0.
        int_state%TRACKER_IFIX=0.
        int_state%TRACKER_JFIX=0.
        int_state%TRACKER_RMW=0.
        int_state%TRACKER_PMIN=0.
        int_state%TRACKER_VMAX=0.
!###    Tracker 1D integer
        DO I=1,TRACK_MAX_OLD
          int_state%TRACK_OLD_NTSD(I)=0
        ENDDO
!###    Tracker 1D real
        DO I=1,TRACK_MAX_OLD
          int_state%TRACK_OLD_LAT(I)=0.
          int_state%TRACK_OLD_LON(I)=0.
        ENDDO
!###    Tracker 2D integer
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%TRACKER_FIXES(I,J)=-999.
        ENDDO
        ENDDO
!###    Tracker 2D real
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%MEMBRANE_MSLP(I,J)=0.
          int_state%P850RV(I,J)=0.
          int_state%P700RV(I,J)=0.
          int_state%P850WIND(I,J)=0.
          int_state%P700WIND(I,J)=0.
          int_state%P500U(I,J)=0.
          int_state%P500V(I,J)=0.
          int_state%P700U(I,J)=0.
          int_state%P700V(I,J)=0.
          int_state%P850U(I,J)=0.
          int_state%P850V(I,J)=0.
          int_state%P850Z(I,J)=0.
          int_state%P700Z(I,J)=0.
          int_state%M10WIND(I,J)=0.
          int_state%M10RV(I,J)=0.
          int_state%SP850RV(I,J)=0.
          int_state%SP700RV(I,J)=0.
          int_state%SP850WIND(I,J)=0.
          int_state%SP700WIND(I,J)=0.
          int_state%SP850Z(I,J)=0.
          int_state%SP700Z(I,J)=0.
          int_state%SM10WIND(I,J)=0.
          int_state%SM10RV(I,J)=0.
          int_state%SMSLP(I,J)=0.
          int_state%TRACKER_ANGLE(I,J)=0.
          int_state%TRACKER_DISTSQ(I,J)=0.
        ENDDO
        ENDDO
!! End of tracker variables
!
!-----------------------------------------------------------------------
!***  Initialize the timer variables now.
!-----------------------------------------------------------------------
!
        TD=>TIMERS(MY_DOMAIN_ID_LOC)                                       !<-- Abbreviate the name of this domain's timers
!
        td%adv1_tim=0.
        td%adv2_tim=0.
        td%bocoh_tim=0.
        td%bocov_tim=0.
        td%cdwdt_tim=0.
        td%cdzdt_tim=0.
        td%consts_tim=0.
        td%ddamp_tim=0.
        td%dht_tim=0.
        td%exch_dyn=0.
        td%exch_phy=0.
        td%exch_tim=0.
        td%fftfhn_tim=0.
        td%fftfwn_tim=0.
        td%hdiff_tim=0.
        td%mono_tim=0.
        td%pdtsdt_tim=0.
        td%pgforce_tim=0.
        td%poavhn_tim=0.
        td%polehn_tim=0.
        td%pole_swap_tim=0.
        td%polewn_tim=0.
        td%prefft_tim=0.
        td%presmud_tim=0.
        td%solver_init_tim=0.
        td%solver_dyn_tim=0.
        td%solver_phy_tim=0.
        td%swaphn_tim=0.
        td%swapwn_tim=0.
        td%updatet_tim=0.
        td%updateuv_tim=0.
        td%updates_tim=0.
        td%vsound_tim=0.
        td%vtoa_tim=0.
!
        td%cucnvc_tim=0.
        td%gsmdrive_tim=0.
        td%cltend_tim=0.
        td%rfupdate_tim=0.
        td%tqadjust_tim=0.
        td%h_to_v_tim=0.
        td%radiation_tim=0.
        td%rdtemp_tim=0.
        td%turbl_tim=0.
        td%adjppt_tim=0.
        td%gfs_phy_tim=0.
!
!-----------------------------------------------------------------------
!
        ITS=int_state%ITS
        ITE=int_state%ITE
        JTS=int_state%JTS
        JTE=int_state%JTE
        IMS=int_state%IMS
        IME=int_state%IME
        JMS=int_state%JMS
        JME=int_state%JME
        IDS=int_state%IDS
        IDE=int_state%IDE
        JDS=int_state%JDS
        JDE=int_state%JDE
!
        IHALO=int_state%IHALO    
        JHALO=int_state%JHALO    
!
        MYPE=int_state%MYPE
        MY_DOMAIN_ID=int_state%MY_DOMAIN_ID
        MPI_COMM_COMP=int_state%MPI_COMM_COMP
        NUM_PES=int_state%NUM_PES
!
        DO N=1,8
          MY_NEB(N)=int_state%MY_NEB(N)
        ENDDO
!
!-----------------------------------------------------------------------
!***  Extract all forecast tasks' horizontal subdomain limits
!***  from the Solver import state and give them to the
!***  Solver internal state.
!***  This is necessary if quilting is selected because these
!***  limits will be taken from the Solver internal state,
!***  placed into the Write components' import states and
!***  used for the combining of local domain data onto the
!***  global domain.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Local Domain Limits to Solver Internal State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state    =IMP_STATE                      &  !<-- The Solver import state
                              ,name     ='LOCAL_ISTART'                 &  !<-- Name of the attribute to extract
                              ,valueList=int_state%LOCAL_ISTART         &  !<-- Insert Attribute into Solver internal state
                              ,rc       =RC)
!
        CALL ESMF_AttributeGet(state    =IMP_STATE                      &  !<-- The Solver import state
                              ,name     ='LOCAL_IEND'                   &  !<-- Name of the attribute to extract
                              ,valueList=int_state%LOCAL_IEND           &  !<-- Insert Attribute into Solver internal state
                              ,rc       =RC)
!
        CALL ESMF_AttributeGet(state    =IMP_STATE                      &  !<-- The Solver import state
                              ,name     ='LOCAL_JSTART'                 &  !<-- Name of the attribute to extract
                              ,valueList=int_state%LOCAL_JSTART         &  !<-- Insert Attribute into Solver internal state
                              ,rc       =RC)
!
        CALL ESMF_AttributeGet(state    =IMP_STATE                      &  !<-- The Solver import state
                              ,name     ='LOCAL_JEND'                   &  !<-- Name of the attribute to extract
                              ,valueList=int_state%LOCAL_JEND           &  !<-- Insert Attribute into Solver internal state
                              ,rc       =RC)

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Fill the ESMF Bundle with the user-selected boundary variables
!***  and also the generalized boundary object that Solver must use
!***  when handling those boundary variables.  This must take place 
!***  here because it must follow the creation of the Solver's
!***  internal state but precede the call to the read routine.  For
!***  restarted runs the read routine must allocate an object to
!***  hold special boundary data from the restart files and the size
!***  of tha object depends on values determined when the ESMF Bundle
!***  is filled.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Solver Init Extracts BC Bundle"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_StateGet(state      =IMP_STATE                        &  !<-- The Solver component import state
                          ,itemname   ='Bundle_nestbc'                  &  !<-- Name of Bundle of selected BC variables
                          ,fieldbundle=BUNDLE_NESTBC                    &  !<-- The Bundle 
                          ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        UBOUND_VARS=SIZE(int_state%VARS)
!
        CALL BUILD_BC_BUNDLE(GRID                                       &  !<-- Add Solver int state variables to the nest BC Bundle
                            ,LNSH,LNSV                                  &
                            ,IHALO,JHALO                                &
                            ,UBOUND_VARS                                &
                            ,int_state%VARS                             &
                            ,MY_DOMAIN_ID                               &
                            ,BUNDLE_NESTBC                              &
                            ,int_state%BND_VARS_H                       &
                            ,int_state%BND_VARS_V                       &
                            ,int_state%NVARS_BC_2D_H                    &
                            ,int_state%NVARS_BC_3D_H                    &
                            ,int_state%NVARS_BC_4D_H                    &
                            ,int_state%NVARS_BC_2D_V                    &
                            ,int_state%NVARS_BC_3D_V                    &
                            ,int_state%NLEV_H                           &
                            ,int_state%NLEV_V                           &
                            ,int_state%N_BC_3D_H                        &
                               )
!
!-----------------------------------------------------------------------
!***  The input file is about to be read and halo exchanges will be
!***  done in conjunction with that process.  The halo exchange
!***  routines require 15 domain-related variables so set them now.
!-----------------------------------------------------------------------
!
        CALL SET_DOMAIN_SPECS(int_state%ITS,int_state%ITE               &
                             ,int_state%JTS,int_state%JTE               &
                             ,int_state%IMS,int_state%IME               &
                             ,int_state%JMS,int_state%JME               &
                             ,int_state%IDS,int_state%IDE               &
                             ,int_state%JDS,int_state%JDE               &
                             ,int_state%IHALO,int_state%JHALO           &
                             ,int_state%MY_DOMAIN_ID                    &
                             ,int_state%MYPE                            &
                             ,int_state%MY_NEB                          &
                             ,int_state%MPI_COMM_COMP                   &
                             ,int_state%NUM_PES                         &
                             ,LOCAL_ISTART_IN=int_state%LOCAL_ISTART    &
                             ,LOCAL_IEND_IN=int_state%LOCAL_IEND        &
                             ,LOCAL_JSTART_IN=int_state%LOCAL_JSTART    &
                             ,LOCAL_JEND_IN=int_state%LOCAL_JEND        &
                              )
!
!-----------------------------------------------------------------------
!***  Read the input file.
!-----------------------------------------------------------------------
!
        KSS=1        
        KSE=int_state%NUM_TRACERS_MET
!
        ITS_H2=MAX(ITS-2,int_state%IDS)
        ITE_H2=MIN(ITE+2,int_state%IDE)
        JTS_H2=MAX(JTS-2,int_state%JDS)
        JTE_H2=MIN(JTE+2,int_state%JDE)
!
        btim=timef()
!

!       write(0,*)'int_state%NEMSIO_INPUT=',int_state%NEMSIO_INPUT  !wang
        IF(.NOT.int_state%NEMSIO_INPUT)THEN
!
          CALL READ_BINARY(INT_STATE                                    &
                          ,MY_DOMAIN_ID                                 &
                          ,MPI_COMM_COMP                                &
                          ,int_state%MYPE                               &
                          ,int_state%ITS,int_state%ITE                  &
                          ,int_state%JTS,int_state%JTE                  &
                          ,int_state%IMS,int_state%IME                  &
                          ,int_state%JMS,int_state%JME                  &
                          ,int_state%IDS,int_state%IDE                  &
                          ,int_state%JDS,int_state%JDE                  &
                          ,ITS_H2,ITE_H2,JTS_H2,JTE_H2                  &
                          ,LM                                           &
                          ,RC)
!
          IF (RC /= 0) THEN
            RC_INIT = RC
            RETURN
          END IF
!
        ELSE
!
!        write(0,*) 'mype=',mype,'call read_nemsio'
          CALL READ_NEMSIO(int_state,MY_DOMAIN_ID,RC)
!
          IF (RC /= 0) THEN
            RC_INIT = RC
            RETURN
          END IF
!
        ENDIF
!rv
!  Use this (OPER) for operational run, for having vertical velocity
!  in history file (00hr) when starting from restart file
!rv
        IF(int_state%OPER) THEN
          DO L=1,LM
            DO J=JMS,JME
              DO I=IMS,IME
                int_state%W_TOT(I,J,L)=int_state%W(I,J,L)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
!rv
!
        if (mype==-9999) then
          write(0,*)'solver'
          write(0,*)'ihr,ihrst,lpt2,ntstm=',int_state%ihr,int_state%ihrst,int_state%lpt2,int_state%ntstm
          write(0,*)'idat=',int_state%idat(1),int_state%idat(2),int_state%idat(3)
          write(0,*)'dsg1=',minval(int_state%dsg1),maxval(int_state%dsg1)
          write(0,*)'pdsg1=',minval(int_state%pdsg1),maxval(int_state%pdsg1)
          write(0,*)'psgml1=',minval(int_state%psgml1),maxval(int_state%psgml1)
          write(0,*)'sgml1=',minval(int_state%sgml1),maxval(int_state%sgml1)
          write(0,*)'sgml2=',minval(int_state%sgml2),maxval(int_state%sgml2)
          write(0,*)'psg1=',minval(int_state%psg1),maxval(int_state%psg1)
          write(0,*)'sg1=',minval(int_state%sg1),maxval(int_state%sg1)
          write(0,*)'sg2=',minval(int_state%sg2),maxval(int_state%sg2)
          write(0,*)'fis=',minval(int_state%fis),maxval(int_state%fis)
          write(0,*)'pd=',minval(int_state%pd),maxval(int_state%pd)
          write(0,*)'pdo=',minval(int_state%pdo),maxval(int_state%pdo)
          write(0,*)'sice=',minval(int_state%sice),maxval(int_state%sice)
          write(0,*)'sm=',minval(int_state%sm),maxval(int_state%sm)
          write(0,*)'cw=',minval(int_state%cw),maxval(int_state%cw)
          write(0,*)'dwdt=',minval(int_state%dwdt),maxval(int_state%dwdt)
          write(0,*)'q=',minval(int_state%q),maxval(int_state%q)
          write(0,*)'q2=',minval(int_state%q2),maxval(int_state%q2)
          write(0,*)'o3=',minval(int_state%o3),maxval(int_state%o3)
          write(0,*)'omgalf=',minval(int_state%omgalf),maxval(int_state%omgalf)
          write(0,*)'div=',minval(int_state%div),maxval(int_state%div)
          write(0,*)'z=',minval(int_state%z),maxval(int_state%z)
          write(0,*)'rtop=',minval(int_state%rtop),maxval(int_state%rtop)
          write(0,*)'tcu=',minval(int_state%tcu),maxval(int_state%tcu)
          write(0,*)'tcv=',minval(int_state%tcv),maxval(int_state%tcv)
          write(0,*)'tct=',minval(int_state%tct),maxval(int_state%tct)
          write(0,*)'t=',minval(int_state%t),maxval(int_state%t)
          write(0,*)'tp=',minval(int_state%tp),maxval(int_state%tp)
          write(0,*)'u=',minval(int_state%u),maxval(int_state%u)
          write(0,*)'up=',minval(int_state%up),maxval(int_state%up)
          write(0,*)'v=',minval(int_state%v),maxval(int_state%v)
          write(0,*)'vp=',minval(int_state%vp),maxval(int_state%vp)
          write(0,*)'w=',minval(int_state%w),maxval(int_state%w)
          write(0,*)'w_tot=',minval(int_state%w_tot),maxval(int_state%w_tot)
          write(0,*)'pint=',minval(int_state%pint),maxval(int_state%pint)
          write(0,*)'tracers=',minval(int_state%tracers),maxval(int_state%tracers)
!         write(0,*)'sp=',minval(int_state%sp),maxval(int_state%sp)
          write(0,*)'run=',int_state%run 
        endif
!
!-----------------------------------------------------------------------
!***  Check if starting Date/Time in input data file agrees with
!***  the configure file.
!-----------------------------------------------------------------------
!
        IF(.NOT.int_state%RESTART.AND.MYPE==0)THEN
          IF(int_state%START_HOUR /=int_state%IHRST.OR.                 &
             int_state%START_DAY  /=int_state%IDAT(1).OR.               &
             int_state%START_MONTH/=int_state%IDAT(2).OR.               &
             int_state%START_YEAR /=int_state%IDAT(3))THEN
            WRITE(0,*)' *** WARNING *** WARNING *** WARNING *** '
            WRITE(0,*)' *** WARNING *** WARNING *** WARNING *** '
            WRITE(0,*)' DATES IN INPUT AND CONFIGURE FILES DISAGREE!!'
            WRITE(0,*)' INPUT: HOUR=',int_state%IHRST                   &
                      ,       ' DAY=',int_state%IDAT(1)                 &
                      ,     ' MONTH=',int_state%IDAT(2)                 &
                      ,      ' YEAR=',int_state%IDAT(3)
            WRITE(0,*)' CONFIG: HOUR=',int_state%START_HOUR             &
                      ,        ' DAY=',int_state%START_DAY              &
                      ,      ' MONTH=',int_state%START_MONTH            &
                      ,       ' YEAR=',int_state%START_YEAR
            WRITE(0,*)' *** WARNING *** WARNING *** WARNING *** '
            WRITE(0,*)' *** WARNING *** WARNING *** WARNING *** '
          ENDIF
        ENDIF
!
!-----------------------------------------------------------------------
!
        td%solver_init_tim=td%solver_init_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Nested domains do not have boundary condition files since the
!***  boundary values come from their parents.  However the boundary
!***  variable arrays need to contain initial values before tendencies
!***  from the parent can be added.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Retrieve the Nest/Not_A_Nest flag from the Solver import state.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Get Nest/Not-a-Nest Flag from Solver Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state=IMP_STATE                          &  !<-- The Solver import state
                              ,name ='I-Am-A-Nest Flag'                 &  !<-- Name of variable to get from Solver import state
                              ,value=I_AM_A_NEST                        &  !<-- Put extracted value here
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        int_state%I_AM_A_NEST=I_AM_A_NEST
!
        IF(I_AM_A_NEST)THEN
!
!-----------------------------------------------------------------------
!***  Also we need to retrieve the Parent-Child timestep ratio in order
!***  to know how often to update the boundary tendencies.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Get Parent-Child Time Ratio from Solver Import State"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeGet(state=IMP_STATE                         &  !<-- The Solver import state
                                ,name ='Parent-Child Time Ratio'         &  !<-- Name of variable to get from Solver import state
                                ,value=int_state%PARENT_CHILD_TIME_RATIO &  !<-- Put extracted value here
                                ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Does this nested domain move?
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Get Nest Move Flag from Solver Import State"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeGet(state=IMP_STATE                        &  !<-- The Solver import state
                                ,name ='My Domain Moves'                &  !<-- Name of variable to get from Solver import state
                                ,value=int_state%MY_DOMAIN_MOVES        &  !<-- Put extracted value here
                                ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Currently moving nests are not allowed to use gravity wave drag.
!***  One quantity used in that parameterization is the mountains'
!***  angle with respect to east.  From the moving nest's perspective
!***  the mountains are moving and thus that angle would need to be
!***  updated with each shift of the domain.  That is not handled
!***  yet in the code.
!-----------------------------------------------------------------------
!
          IF(int_state%MY_DOMAIN_MOVES)THEN
!
            int_state%GWDFLG=.FALSE.
!
          ENDIF
!
!-----------------------------------------------------------------------
!
        ENDIF
!
!-----------------------------------------------------------------------
!***  Assign grid-related constants after dereferencing needed variables.
!-----------------------------------------------------------------------
!
        btim=timef()
!
        CALL CONSTS(int_state%GLOBAL                                    &
                   ,int_state%DT                                        &
                   ,int_state%SMAG2                                     &
                   ,int_state%CODAMP,int_state%WCOR                     &
                   ,int_state%PT                                        &
                   ,int_state%TPH0D,int_state%TLM0D                     &
                   ,int_state%SBD,int_state%WBD                         &
                   ,int_state%DPHD,int_state%DLMD                       &
                   ,int_state%DXH,int_state%RDXH                        &
                   ,int_state%DXV,int_state%RDXV                        &
                   ,int_state%DYH,int_state%RDYH                        &
                   ,int_state%DYV,int_state%RDYV                        &
                   ,int_state%DDV,int_state%RDDV                        &
                   ,int_state%DDMPU,int_state%DDMPV                     &
                   ,int_state%EF4T,int_state%WPDAR                      &
                   ,int_state%FCP,int_state%FDIV                        &
                   ,int_state%CURV,int_state%F                          &
                   ,int_state%FAD,int_state%FAH                         &
                   ,int_state%DARE,int_state%RARE                       &
                   ,int_state%GLAT,int_state%GLON                       &
                   ,int_state%GLAT_SW,int_state%GLON_SW                 &
                   ,int_state%VLAT,int_state%VLON                       &
                   ,int_state%HDACX,int_state%HDACY                     &
                   ,int_state%HDACVX,int_state%HDACVY                   &
                   ,int_state%LNSH,int_state%LNSAD                      &
                   ,int_state%ADV_STANDARD,int_state%ADV_UPSTREAM       &
                   ,int_state%E_BDY,int_state%N_BDY                     &
                   ,int_state%S_BDY,int_state%W_BDY                     &
                   ,int_state%NBOCO,int_state%TBOCO                     &
                   ,MY_DOMAIN_ID,MYPE                                   &
                   ,ITS,ITE,JTS,JTE                                     &
                   ,IMS,IME,JMS,JME                                     &
                   ,IDS,IDE,JDS,JDE )
!
        td%consts_tim=td%consts_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Exchange haloes for some grid-related arrays in case there are
!***  moving nests.
!-----------------------------------------------------------------------
!
        CALL HALO_EXCH                                                  &
             (int_state%GLAT,1                                          &
             ,int_state%GLON,1                                          &
             ,int_state%VLAT,1                                          &
             ,int_state%VLON,1                                          &
             ,3,3)
!
        CALL HALO_EXCH                                                  &
             (int_state%HDACX,1                                         &
             ,int_state%HDACY,1                                         &
             ,int_state%HDACVX,1                                        &
             ,int_state%HDACVY,1                                        &
             ,3,3)
!
        CALL HALO_EXCH                                                  &
             (int_state%F,1                                             &
             ,3,3)
!
!-----------------------------------------------------------------------
!*** Search for lat/lon min/max values and store it in file for
!*** later use in creating GrADS ctl file
!-----------------------------------------------------------------------
!
       CALL LAT_LON_BNDS(int_state%GLAT,int_state%GLON                  &
                       ,mype,num_pes,mpi_comm_comp                      &
                       ,ids,ide,jds,jde                                 &
                       ,ims,ime,jms,jme                                 &
                       ,its,ite,jts,jte                                 &
                       ,my_domain_id )
!
!
!-----------------------------------------------------------------------
!***  Initialize the FFT filters.
!-----------------------------------------------------------------------
!
        IF(int_state%GLOBAL)THEN
          btim=timef()
!
          CALL PREFFT(int_state%DLMD,int_state%DPHD,int_state%SBD,LM      &
                     ,int_state%KHFILT,int_state%KVFILT                   &
                     ,int_state%HFILT,int_state%VFILT                     &
                     ,int_state%WFFTRH,int_state%NFFTRH                   &
                     ,int_state%WFFTRW,int_state%NFFTRW                   &
                     ,int_state%INPES,int_state%JNPES,int_state%MYPE)
!
          td%prefft_tim=td%prefft_tim+(timef()-btim)
!
        ENDIF
!
!-----------------------------------------------------------------------
!***  Initialize the physics schemes.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Initialize the Physics Schemes"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL PHYSICS_INITIALIZE(int_state%GFS                           &
                               ,int_state%SHORTWAVE                     &
                               ,int_state%LONGWAVE                      &
                               ,int_state%CONVECTION                    &
                               ,int_state%MICROPHYSICS                  &
                               ,int_state%SFC_LAYER                     &
                               ,int_state%TURBULENCE                    &
                               ,int_state%LAND_SURFACE                  &
                               ,int_state%CO2TF                         &
                               ,int_state%NP3D                          &
                               ,int_state%SBD                           &
                               ,int_state%WBD                           &
                               ,int_state%DPHD                          &
                               ,int_state%DLMD                          &
                               ,int_state%TPH0D                         &
                               ,int_state%TLM0D                         &
                               ,MY_DOMAIN_ID                            &
                               ,MYPE                                    &
                               ,MPI_COMM_COMP                           &
                               ,IDS,IDE,JDS,JDE,LM                      &
                               ,IMS,IME,JMS,JME                         &
                               ,ITS,ITE,JTS,JTE                         &
                               ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
! Initialize the storm tracker if needed.
!-----------------------------------------------------------------------
!
        IF(int_state%MYPE<int_state%NUM_PES                             &
                  .AND.                                                 &
           .NOT.int_state%RESTART)THEN
!
           CALL TRACKER_INIT(int_state)
!
        ENDIF
!
!-----------------------------------------------------------------------
!***  Retrieve the ESMF Grid then create the ESMF Fields on that Grid
!***  for the Solver import/export states.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Retrieve ESMF Grid in Solver Initialize"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Insert the value of NUM_TRACERS_TOTAL into the export state.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert NUM_TRACERS into Solver Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Solver export state
                              ,name ='NUM_TRACERS_TOTAL'                &  !<-- The inserted quantity will have this name
                              ,value=int_state%NUM_TRACERS_TOTAL        &  !<-- The value of this is associated with the preceding name
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Also insert the index values of the 4-D Tracers array where
!***  Q and CW reside.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert INDX_Q into Physics Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Physics export state
                              ,name ='INDX_Q'                           &  !<-- The inserted quantity will have this name
                              ,value=int_state%INDX_Q                   &  !<-- The location of Q in TRACERS
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert INDX_CW into Physics Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Physics export state
                              ,name ='INDX_CW'                          &  !<-- The inserted quantity will have this name
                              ,value=int_state%INDX_CW                  &  !<-- The value of this is associated with the preceding name
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert INDX_QC into Physics Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Physics export state
                              ,name ='INDX_QC'                          &  !<-- The inserted quantity will have this name
                              ,value=int_state%INDX_QC                  &  !<-- The value of this is associated with the preceding name
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert INDX_QI into Physics Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Physics export state
                              ,name ='INDX_QI'                          &  !<-- The inserted quantity will have this name
                              ,value=int_state%INDX_QI                  &  !<-- The value of this is associated with the preceding name
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert INDX_QR into Physics Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Physics export state
                              ,name ='INDX_QR'                          &  !<-- The inserted quantity will have this name
                              ,value=int_state%INDX_QR                  &  !<-- The value of this is associated with the preceding name
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert INDX_QS into Physics Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Physics export state
                              ,name ='INDX_QS'                          &  !<-- The inserted quantity will have this name
                              ,value=int_state%INDX_QS                  &  !<-- The value of this is associated with the preceding name
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert INDX_QG into Physics Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Physics export state
                              ,name ='INDX_QG'                          &  !<-- The inserted quantity will have this name
                              ,value=int_state%INDX_QG                  &  !<-- The value of this is associated with the preceding name
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert INDX_NI into Physics Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Physics export state
                              ,name ='INDX_NI'                          &  !<-- The inserted quantity will have this name
                              ,value=int_state%INDX_NI                  &  !<-- The value of this is associated with the preceding name
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert INDX_NR into Physics Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Physics export state
                              ,name ='INDX_NR'                          &  !<-- The inserted quantity will have this name
                              ,value=int_state%INDX_NR                  &  !<-- The value of this is associated with the preceding name
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Insert this task's integration index limits into the
!***  export state along with the full domain limits.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Add Task Integration Limits to Solver Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Solver export state
                              ,name ='ITS'                              &  !<-- The inserted quantity will have this name
                              ,value=int_state%ITS                      &  !<-- The value of this is associated with the preceding name
                              ,rc   =RC)
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Solver export state
                              ,name ='ITE'                              &  !<-- The inserted quantity will have this name
                              ,value=int_state%ITE                      &  !<-- The value of this is associated with the preceding name
                              ,rc   =RC)
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Solver export state
                              ,name ='JTS'                              &  !<-- The inserted quantity will have this name
                              ,value=int_state%JTS                      &  !<-- The value of this is associated with the preceding name
                              ,rc   =RC)
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Solver export state
                              ,name ='JTE'                              &  !<-- The inserted quantity will have this name
                              ,value=int_state%JTE                      &  !<-- The value of this is associated with the preceding name
                              ,rc   =RC)
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Solver export state
                              ,name ='LM'                               &  !<-- The inserted quantity will have this name
                              ,value=int_state%LM                       &  !<-- The value of this is associated with the preceding name
                              ,rc   =RC)
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Solver export state
                              ,name ='NHALO'                            &  !<-- The inserted quantity will have this name
                              ,value=int_state%NHALO                    &  !<-- The value of this is associated with the preceding name
                              ,rc   =RC)
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Solver export state
                              ,name ='IDS'                              &  !<-- The inserted quantity will have this name
                              ,value=int_state%IDS                      &  !<-- The value of this is associated with the preceding name
                              ,rc   =RC)
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Solver export state
                              ,name ='IDE'                              &  !<-- The inserted quantity will have this name
                              ,value=int_state%IDE                      &  !<-- The value of this is associated with the preceding name
                              ,rc   =RC)
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Solver export state
                              ,name ='JDS'                              &  !<-- The inserted quantity will have this name
                              ,value=int_state%JDS                      &  !<-- The value of this is associated with the preceding name
                              ,rc   =RC)
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Solver export state
                              ,name ='JDE'                              &  !<-- The inserted quantity will have this name
                              ,value=int_state%JDE                      &  !<-- The value of this is associated with the preceding name
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Insert the domain's top pressure, the pressure thickness of the
!***  pressure domain, the mid-layer pressures in the pressure domain
!***  and the mid-layer sigmas in the sigma domain.
!-----------------------------------------------------------------------
!
        LMP1=LM+1
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert PT into Solver Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Solver export state
                              ,name ='PT'                               &  !<-- The inserted quantity will have this name
                              ,value=int_state%PT                       &  !<-- The value of this is associated with the preceding name
                              ,rc   =RC)
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Solver export state
                              ,name ='PDTOP'                            &  !<-- The inserted quantity will have this name
                              ,value=int_state%PDTOP                    &  !<-- The value of this is associated with the preceding name
                              ,rc   =RC)

        CALL ESMF_AttributeSet(state    =EXP_STATE                      &  !<-- The Solver export state
                              ,name     ='PSGML1'                       &  !<-- The inserted quantity will have this name
                              ,itemCount=LM                             &  !<-- The data has this many items
                              ,valueList=int_state%PSGML1               &  !<-- The value of this is associated with the preceding name
                              ,rc       =RC)
!
        CALL ESMF_AttributeSet(state    =EXP_STATE                      &  !<-- The Solver export state
                              ,name     ='SGML2'                        &  !<-- The inserted quantity will have this name
                              ,itemCount=LM                             &  !<-- The data has this many items
                              ,valueList=int_state%SGML2                &  !<-- The value of this is associated with the preceding name
                              ,rc       =RC)
!
        CALL ESMF_AttributeSet(state    =EXP_STATE                      &  !<-- The Solver export state
                              ,name     ='SG1'                          &  !<-- The inserted quantity will have this name
                              ,itemCount=LMP1                           &  !<-- The data has this many items
                              ,valueList=int_state%SG1                  &  !<-- The value of this is associated with the preceding name
                              ,rc       =RC)
!
        CALL ESMF_AttributeSet(state    =EXP_STATE                      &  !<-- The Solver export state
                              ,name     ='SG2'                          &  !<-- The inserted quantity will have this name
                              ,itemCount=LMP1                           &  !<-- The data has this many items
                              ,valueList=int_state%SG2                  &  !<-- The value of this is associated with the preceding name
                              ,rc       =RC)
!
        CALL ESMF_AttributeSet(state    =EXP_STATE                      &  !<-- The Solver export state
                              ,name     ='DSG2'                         &  !<-- The inserted quantity will have this name
                              ,itemCount=LM                             &  !<-- The data has this many items
                              ,valueList=int_state%DSG2                 &  !<-- The value of this is associated with the preceding name
                              ,rc       =RC)
!
        CALL ESMF_AttributeSet(state    =EXP_STATE                      &  !<-- The Solver export state
                              ,name     ='PDSG1'                        &  !<-- The inserted quantity will have this name
                              ,itemCount=LM                             &  !<-- The data has this many items
                              ,valueList=int_state%PDSG1                &  !<-- The value of this is associated with the preceding name
                              ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Insert DXH and DYH into the Solver export state.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert DYH into the Solver Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Solver export state
                              ,name ='DYH'                              &  !<-- The inserted quantity will have this name
                              ,value=int_state%DYH                      &  !<-- The value of this is associated with the preceding name
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        KOUNT=SIZE(int_state%DXH)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert DXH into the Solver Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state    =EXP_STATE                      &  !<-- The Solver export state
                              ,name     ='DXH'                          &  !<-- The inserted quantity will have this name
                              ,itemCount=KOUNT                          &  !<-- The data has this many items
                              ,valueList=int_state%DXH                  &  !<-- The value of this is associated with the preceding name
                              ,rc       =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Insert DPHD and JM into the export state.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert DPHD,DLMD,JM into the Solver Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Solver export state
                              ,name ='DLMD'                             &  !<-- The inserted quantity will have this name
                              ,value=int_state%DLMD                     &  !<-- The value of this is associated with the preceding name
                              ,rc   =RC)
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Solver export state
                              ,name ='DPHD'                             &  !<-- The inserted quantity will have this name
                              ,value=int_state%DPHD                     &  !<-- The value of this is associated with the preceding name
                              ,rc   =RC)
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Solver export state
                              ,name ='JM'                               &  !<-- The inserted quantity will have this name
                              ,value=int_state%JM                       &  !<-- The value of this is associated with the preceding name
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Insert the value of LNSH and LNSV (the width of the
!***  blending region along the boundaries for H and V points).
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert LNSH, LNSV into Solver Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Solver export state
                              ,name ='LNSH'                             &  !<-- The inserted quantity will have this name
                              ,value=int_state%LNSH                     &  !<-- The value of this is associated with the preceding name
                              ,rc   =RC)
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Solver export state
                              ,name ='LNSV'                             &  !<-- The inserted quantity will have this name
                              ,value=int_state%LNSV                     &  !<-- The value of this is associated with the preceding name
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Insert the geographic latitude and longitude of the grid points
!***  into the export state.  From there they will be updated in 
!***  DOMAIN_RUN when a moving nest moves.  The central lat/lon
!***  of the nest's rotated system and the angular grid increments
!***  are also needed.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Create Field from H-pt Geographic Latitude"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        FIELD=ESMF_FieldCreate(            GRID                         &  !<-- The ESMF Grid
                              ,            int_state%GLAT               &  !<-- The geographic latitude on H points
                              ,totalUWidth=(/IHALO,JHALO/)              &  !<-- Upper bound of halo region
                              ,totalLWidth=(/IHALO,JHALO/)              &  !<-- Lower bound of halo region
                              ,name       ='GLAT'                       &  !<-- Name of Field
                              ,indexFlag  =ESMF_INDEX_GLOBAL            &
                              ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Add GLAT to the Solver Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_StateAddReplace(EXP_STATE                              &  !<-- The Solver export state
                          ,(/FIELD/)                     &  !<-- Field with H-pt geographic lat
                          ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Create Field from H-pt Geographic Longitude"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        FIELD=ESMF_FieldCreate(            GRID                         &  !<-- The ESMF Grid
                              ,            int_state%GLON               &  !<-- The geographic longitude on H points
                              ,totalUWidth=(/IHALO,JHALO/)              &  !<-- Upper bound of halo region
                              ,totalLWidth=(/IHALO,JHALO/)              &  !<-- Lower bound of halo region
                              ,name       ='GLON'                       &  !<-- Name of Field
                              ,indexFlag  =ESMF_INDEX_GLOBAL            &
                              ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Add GLON to the Solver Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_StateAddReplace(EXP_STATE                              &  !<-- The Solver export state
                          ,(/FIELD/)                     &  !<-- Field with H-pt geographic lon
                          ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Create Field from V-pt Geographic Latitude"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        FIELD=ESMF_FieldCreate(            GRID                         &  !<-- The ESMF Grid
                              ,            int_state%VLAT               &  !<-- The geographic latitude on V points
                              ,totalUWidth=(/IHALO,JHALO/)              &  !<-- Upper bound of halo region
                              ,totalLWidth=(/IHALO,JHALO/)              &  !<-- Lower bound of halo region
                              ,name       ='VLAT'                       &  !<-- Name of Field
                              ,indexFlag  =ESMF_INDEX_GLOBAL            &
                              ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Add VLAT to the Solver Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_StateAddReplace(EXP_STATE                              &  !<-- The Solver export state
                          ,(/FIELD/)                     &  !<-- Field with V-pt geographic lat
                          ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Create Field from V-pt Geographic Longitude"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        FIELD=ESMF_FieldCreate(            GRID                         &  !<-- The ESMF Grid
                              ,            int_state%VLON               &  !<-- The geographic longitude on V points
                              ,totalUWidth=(/IHALO,JHALO/)              &  !<-- Upper bound of halo region
                              ,totalLWidth=(/IHALO,JHALO/)              &  !<-- Lower bound of halo region
                              ,name       ='VLON'                       &  !<-- Name of Field
                              ,indexFlag  =ESMF_INDEX_GLOBAL            &
                              ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Add VLON to the Solver Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_StateAddReplace(EXP_STATE                              &  !<-- The Solver export state
                          ,(/FIELD/)                     &  !<-- Field with V-pt geographic lon
                          ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert TPH0D, TLM0D into the Solver Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Solver export state
                              ,name ='TPH0D'                            &  !<-- Name of the Attribute
                              ,value=int_state%TPH0D                    &  !<-- The central geo lat of the rotated system
                              ,rc   =RC)
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Solver export state
                              ,name ='TLM0D'                            &  !<-- Name of the Attribute
                              ,value=int_state%TLM0D                    &  !<-- The central geo lon of the rotated system
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert Restart Flag into the Solver Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Solver export state
                              ,name ='RESTART'                          &  !<-- Name of the Attribute
                              ,value=int_state%RESTART                  &  !<-- Is this a restarted run?
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  If this is a nested domain being restarted then it will have
!***  read in the latest values for its SW corner on its parent grid.
!***  Load those into the export state to transfer to the Parent-
!***  Child coupler.  They are only relevant for nests in restarted
!***  runs.  If this is not a nest the values will be dummies and are
!***  never used.  Likewise a moving nest's next move timestep will
!***  have been read from the restart file for a restarted run.
!***  If this is a parent being restarted then it will have read in
!***  the latest value of the next timestep that its moving children
!***  will move.  Add those to the export state to transfer to the
!***  parent-Child coupler.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert SW Corner of Nest into the Solver Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Solver export state
                              ,name ='I_PAR_STA'                        &  !<-- Name of the Attribute
                              ,value=int_state%I_PAR_STA                &  !<-- Parent I of SW corner of this nest
                              ,rc   =RC)
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Solver export state
                              ,name ='J_PAR_STA'                        &  !<-- Name of the Attribute
                              ,value=int_state%J_PAR_STA                &  !<-- Parent J of SW corner of this nest
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert Next Move Timestep into the Solver Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Solver export state
                              ,name ='NEXT_MOVE_TIMESTEP'               &  !<-- Name of the Attribute
                              ,value=int_state%NMTS                     &  !<-- Timestep of the nest's next move
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Let SOLVER_RUN know that the first timestep is special as well
!***  as the first time SOLVER_RUN is executed (which might not be the 
!***  first timestep).
!-----------------------------------------------------------------------
!
        int_state%FIRST_STEP=.TRUE.
        int_state%FIRST_PASS=.TRUE.
!
!-----------------------------------------------------------------------
!***  Set flag for the operational physics suite.
!***  This will be used to save clocktime by skipping
!***  frequent updates of the moist array and instead
!***  update it only when it is needed for physics.
!-----------------------------------------------------------------------
!
        int_state%OPERATIONAL_PHYSICS=.FALSE.
!
        IF((int_state%SHORTWAVE   =='gfdl' .OR.                         &
            int_state%SHORTWAVE   =='rrtm').AND.                        &
           (int_state%LONGWAVE    =='gfdl' .OR.                         &
            int_state%LONGWAVE    =='rrtm').AND.                        &
            int_state%SFC_LAYER   =='myj'  .AND.                        &
            int_state%TURBULENCE  =='myj'  .AND.                        &
           (int_state%CONVECTION  =='bmj'  .OR.                         &
            int_state%CONVECTION  =='none').AND.                        &
           (int_state%MICROPHYSICS=='fer'  .OR.                         &
            int_state%MICROPHYSICS=='fer_hires') ) THEN
!
          int_state%OPERATIONAL_PHYSICS=.TRUE.
!
        ENDIF
!
!-----------------------------------------------------------------------
!
      ENDIF fcst_tasks2
!
!-----------------------------------------------------------------------
!
      td%solver_init_tim=td%solver_init_tim+(timef()-btim0)
!
!-----------------------------------------------------------------------
!
      IF(RC_INIT==ESMF_SUCCESS)THEN
!       WRITE(0,*)'SOLVER INITIALIZE STEP SUCCEEDED'
      ELSE
        WRITE(0,*)'SOLVER INITIALIZE STEP FAILED RC_INIT=',RC_INIT
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE SOLVER_INITIALIZE
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE SOLVER_RUN (GRID_COMP                                  &
                            ,IMP_STATE                                  &
                            ,EXP_STATE                                  &
                            ,CLOCK_ATM                                  &
                            ,RC_RUN)
!
!-----------------------------------------------------------------------
!***  The integration of each timestep of the model Solver is done
!***  through this routine.
!-----------------------------------------------------------------------
!
      USE MODULE_CONSTANTS,ONLY : CP,G,R,RHOWATER,STBOLT,XLV,R_D,R_V,PI
!
      USE MODULE_DYNAMICS_ROUTINES,ONLY: ADV1,ADV2                      &
                                        ,CDWDT,CDZDT,DDAMP,DHT          &
                                        ,HDIFF                          &
                                        ,MONO,PDTSDT,PGFORCE            &
                                        ,UPDATES,UPDATET,UPDATEUV       &
                                        ,VSOUND,VTOA
!
      USE MODULE_FLTBNDS,ONLY: BOCOH,BOCOV,FFTFHN,FFTFUVN               &
                              ,POAVHN,READ_BC                           &
                              ,WRITE_BC
!
!-----------------------------------------------------------------------
!***  The following USEs are needed only for GFS physics:
!-----------------------------------------------------------------------
!
      USE N_NAMELIST_PHYSICS_DEF,      ONLY: FHSWR, FDAER               &
                                            ,IAER,IALB,ICO2,IEMS,ICTM   &
                                            ,IOVR_LW,IOVR_SW,ISOL       &
                                            ,LDIAG3D,LSCCA,LGGFS3D      &
                                            ,LSLWR,LSM,LSSAV,LSSWR      &
                                            ,PRE_RAD,RAS,SASHAL         &
                                            ,SHAL_CNV                   &
                                            ,GEN_COORD_HYBRID           &
                                            ,CDMBGWD,DLQF,CTEI_RM       &
                                            ,BKGD_VDIF_M                &
                                            ,BKGD_VDIF_H,BKGD_VDIF_S    &
                                            ,PSAUTCO,PRAUTCO,EVPCO      &
                                            ,CAL_PRE,MOM4ICE,MSTRAT     &
                                            ,TRANS_TRAC,NST_FCST        &
                                            ,MOIST_ADJ,WMINCO

      USE N_LAYOUT1,                  ONLY : IPT_LATS_NODE_R            &
                                            ,LATS_NODE_R
!
#ifdef USE_GFS_PHYS
      USE DATE_DEF,                   ONLY : FHOUR
!jm      USE MODULE_RADIATION_DRIVER_gfs,    ONLY : GRRAD_gfs,RADINIT_gfs
!jm      USE MODULE_RADIATION_ASTRONOMY_gfs, ONLY : ASTRONOMY
      USE MERSENNE_TWISTER
      USE N_RESOL_DEF,                ONLY : LATR,LONR,LEVR             &
                                            ,NCLD,NFXR,NMTVR            &
                                            ,NTCW,NTOZ                  &
                                            ,THERMODYN_ID, SFCPRESS_ID  &
                                            ,NUM_P2D,NUM_P3D

      USE OZNE_DEF,                   ONLY : LEVOZP,PL_COEFF,PL_PRES
      USE MODULE_RADSW_PARAMETERS_nmmb,    ONLY : TOPFSW_TYPE, SFCFSW_TYPE
      USE MODULE_RADLW_PARAMETERS_nmmb,    ONLY : TOPFLW_TYPE, SFCFLW_TYPE
#endif
      USE MODULE_TRACKER
      USE MODULE_QUASIPOST
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument variables
!------------------------
!
      TYPE(ESMF_GridComp) :: GRID_COMP                                     !<-- The Solver gridded component
!
      TYPE(ESMF_State) :: IMP_STATE                                     &  !<-- The Solver import state
                         ,EXP_STATE                                        !<-- The Solver export state
!
      TYPE(ESMF_Clock) :: CLOCK_ATM                                        !<-- The ATM's ESMF Clock
!
      INTEGER,INTENT(OUT) :: RC_RUN
!
!---------------------
!***  Local variables
!---------------------
!
      INTEGER(kind=KINT) :: IDE,IDS,IME,IMS,ITE,ITS                     &
                           ,JDE,JDS,JME,JMS,JTE,JTS
!
      INTEGER(kind=KINT) :: ITE_B1,ITE_B2,ITE_B1_H1,ITE_B1_H2           &
                           ,ITE_H1,ITE_H2                               &
                           ,ITS_B1,ITS_B2,ITS_B1_H1,ITS_B1_H2           &
                           ,ITS_H1,ITS_H2                               &
                           ,JTE_B1,JTE_B2,JTE_B1_H1,JTE_B1_H2           &
                           ,JTE_H1,JTE_H2                               &
                           ,JTS_B1,JTS_B2,JTS_B1_H1,JTS_B1_H2           &
                           ,JTS_H1,JTS_H2
!
      INTEGER(kind=KINT) :: IHALO,JHALO,MPI_COMM_COMP,MY_DOMAIN_ID      &
                           ,MYPE,NUM_PES
!
      INTEGER(kind=KINT) :: DFIHR,I,IER,INPES,IRTN,ISTAT,J,JNPES        &
                           ,K,KFLIP,KS,KSE1,L,N,NSTEPS_HISTORY          &
                           ,NTIMESTEP,NTIMESTEP_BC,NTIMESTEP_RAD        &
                           ,RC,ICLTEND                                  &
                           ,WRITE_BC_FLAG,WRITE_BC_FLAG_NEST
!
      INTEGER(kind=KINT) :: FILTER_METHOD,FILTER_METHOD_LAST            &
                           ,JULDAY,JULYR                                &
                           ,NPRECIP,NSTEPS_PER_CHECK,NSTEPS_PER_HOUR    &
                           ,NSTEPS_PER_RESET !rh_hold ,USE_RADAR_FIRST
!
      LOGICAL(kind=KLOG) :: USE_RADAR
!
      INTEGER(kind=KINT),SAVE :: HDIFF_ON                               &
                                ,PARENT_CHILD_TIME_RATIO
!
      INTEGER(kind=ESMF_KIND_I8) :: NTIMESTEP_ESMF
!
      LOGICAL(kind=KLOG) :: READBC                                      &
                           ,E_BDY,N_BDY,S_BDY,W_BDY
!
      TYPE(ESMF_TimeInterval) :: DT_ESMF                                   !<-- The ESMF fundamental timestep (s)
!
      TYPE(SOLVER_INTERNAL_STATE),POINTER :: INT_STATE                     !<-- The Solver internal state pointer 
!
      TYPE(WRAP_SOLVER_INT_STATE) :: WRAP                                  !<-- The F90 'wrap' for the Solver internal state
!
!-----------------------------------------------------------------------
!***  SAVEs are for dereferenced constant variables.
!-----------------------------------------------------------------------
!
      INTEGER(kind=KINT),SAVE :: IDTADT,IDTADTQ,IFACT,IHRSTBC                   &
                                ,INTEGER_DT                             &
                                ,KSE,KSS                                &
                                ,LNSAD,LNSH,LNSV,LPT2,NBOCO             &
                                ,N_PRINT_STATS                          &  !<--- Timesteps between statistics prints
                                ,NUMERATOR_DT                           &
                                ,IDENOMINATOR_DT,IFLAG
!
      INTEGER(kind=KINT),DIMENSION(3),SAVE :: IDATBC
!
      INTEGER(kind=KINT),DIMENSION(8)  :: IDAT,JDAT
      INTEGER(kind=KINT),DIMENSION(13) :: DAYS
!
      REAL(kind=KFPT) :: FICE,FRAIN,QI,QR,QW,SECONDS_TOTAL,WC
!
      REAL(kind=KFPT) :: DT,DT_TEST,DT_TEST_RATIO,DTPHY
!
      REAL(kind=KFPT),SAVE :: DDMPV                                     &
                             ,DYH,DYV,EF4T,PDTOP,PT                     &
                             ,RDYH,RDYV,TBOCO
!
!rh_hold      REAL(kind=KFPT),DIMENSION(:,:,:),ALLOCATABLE,SAVE :: RH_HOLD
!
      LOGICAL(kind=KLOG),SAVE :: GLOBAL,HYDRO,RUNBC,SECADV
!
      REAL(kind=KFPT),DIMENSION(:,:),ALLOCATABLE,SAVE :: HDACX_SV       &
                                                       , HDACY_SV       &                       
                                                       , HDACVX_SV      &
                                                       , HDACVY_SV
!
      REAL(kind=KFPT),DIMENSION(:),ALLOCATABLE,SAVE :: DDMPU_SV,FAD_SV, &
                                                       FAH_SV,FCP_SV
!
      REAL(kind=KFPT),SAVE :: DDMPV_SV,EF4T_SV
!

      LOGICAL(kind=KLOG) :: COMPUTE_BC,FIRST_PASS
!
      REAL(kind=KFPT) :: JULIAN,XTIME, FILT_DT, FUND_DT, DTRATIO
!
      INTEGER :: KK
!
      LOGICAL(kind=KLOG) :: CALL_LONGWAVE                               &
                           ,CALL_SHORTWAVE                              &
                           ,CALL_TURBULENCE                             &
                           ,CALL_PRECIP                                 &
                           ,CALL_GFS_PHY                                &
!aligo
                           ,RIME_FACTOR_ADVECT                          &
                           ,RIME_FACTOR_INPUT                           &
!aligo                           
                           ,LOC_PCPFLG
!
      TYPE(ESMF_Time) :: STARTTIME,CURRTIME,SIMULATION_START_TIME
!
      TYPE(ESMF_TimeInterval),SAVE:: REST_OFFSET
!
      TYPE(ESMF_Field) :: HOLD_FIELD
!
      DATA DAYS / 31,28,31,30,31,30,31,31,30,31,30,31,30 /
!
#ifdef USE_GFS_PHYS
!---------------------------------
!***  GFS physics local variables
!---------------------------------
!
      LOGICAL,SAVE                                 :: FIRST=.true.
      LOGICAL                                      :: LPRNT=.false.,NORAD_PRECIP=.false.,CRICK_PROOF=.false., CCNORM=.false.
      LOGICAL                                      :: LSFWD,OPENED,FLIPV,CHANGE,LSSAV_CC,LGOCART=.FALSE.
      LOGICAL                                      :: LSSAV_CPL
      INTEGER,PARAMETER                            :: IFLIP=1,NTRAC=3            !!!!!! later ntrac read form namelist
      INTEGER                                      :: ICWP,IMJM, IDATE(4)
      INTEGER                                      :: ISEED,IDE_GR
      INTEGER ,SAVE                                :: ID,IDAY,IMON,MIDMON,MIDM,MIDP,K1OZ,K2OZ,SEED0
      INTEGER ,DIMENSION(1)                        :: ICSDSW,ICSDLW
      INTEGER ,DIMENSION(:),ALLOCATABLE            :: LONSPERLAR, GLOBAL_LATS_R,NLNSP
!
      REAL (kind=KDBL)                             :: T850,FACOZ,DTLW,DTSW,DTLWI,DTSWI,RTvR,CLSTP,DTP,DTF,SOLHR,RADDT
      REAL (kind=KDBL)                             :: XLVRW,XLVRWI,DTPHS,DTPHSI,RoCP,MINDT
      REAL (kind=KDBL) ,DIMENSION(1)               :: RCS2_V,XLAT,FLGMIN_L,CV,CVB,CVT  ! (cv, cvb, cvt not in use when ntcw-1 > 0)
      REAL (kind=KDBL) ,DIMENSION(1)               :: TSEA,TISFC,ZORL,SLMSK,SNWDPH,WEASD,SNCOVR,SNOALB
      REAL (kind=KDBL) ,DIMENSION(1)               :: XSIHFCS,XSICFCS,XSLPFCS,XTG3FCS,XVEGFCS,XVETFCS,XSOTFCS
      REAL (kind=KDBL) ,DIMENSION(1)               :: ALVSF,ALNSF,ALVWF,ALNWF,FACSF,FACWF
      REAL (kind=KDBL) ,DIMENSION(1)               :: WRK, DPSHC, GQ, RANNUM_V
      REAL (kind=KDBL) ,DIMENSION(1)               :: ORO, EVAP, HFLX, CDQ, QSS
      REAL (kind=KDBL) ,DIMENSION(:),ALLOCATABLE   :: CLDCOV_V,PRSL,PRSLK,GU,GV,GT,GR,VVEL,F_ICE,F_RAIN,R_RIME
      REAL (kind=KDBL) ,DIMENSION(:),ALLOCATABLE   :: ADT,ADU,ADV,PHIL
      REAL (kind=KDBL) ,DIMENSION(:,:),ALLOCATABLE :: GR3,ADR
      REAL (kind=KDBL) ,DIMENSION(:),ALLOCATABLE   :: PRSI,PRSIK,RSGM,PHII
      REAL (kind=KDBL) ,DIMENSION(:),ALLOCATABLE   :: SINLAT_R,COSLAT_R
      REAL (kind=KDBL) ,DIMENSION(:,:),ALLOCATABLE :: XLON,COSZEN,COSZDG,RANN,SINLAT_V,COSLAT_V
      REAL (kind=KDBL) ,DIMENSION(:),ALLOCATABLE   :: RANNUM
!
      REAL (kind=KDBL) ,DIMENSION(39)              :: FLUXR_V
      REAL (kind=KDBL) ,DIMENSION(:,:,:),ALLOCATABLE :: GR1
!
      REAL (kind=KDBL) ,DIMENSION(1)               :: SFALB,TSFLW,SEMIS,SFCDLW,SFCDSW,SFCNSW
      REAL (kind=KDBL) ,DIMENSION(1)               :: NIRBMD_CPL,NIRDFD_CPL,VISBMD_CPL,VISDFD_CPL
      REAL (kind=KDBL) ,DIMENSION(1)               :: NIRBMU_CPL,NIRDFU_CPL,VISBMU_CPL,VISDFU_CPL


      type (topfsw_type), dimension(1) :: topfsw
      type (sfcfsw_type), dimension(1) :: sfcfsw
      type (topflw_type), dimension(1) :: topflw
      type (sfcflw_type), dimension(1) :: sfcflw

      REAL (kind=KDBL) ,DIMENSION(:),ALLOCATABLE   :: SWH,HLW,DKH,RNP
!--- gbphys ---
      LOGICAL                                      :: OLD_MONIN, CNVGWD
      INTEGER                                      :: NEWSAS
      INTEGER ,DIMENSION(2)                        :: NCW
      REAL (kind=KDBL)                             :: CGWF(2),PRSLRD0
      REAL (kind=KDBL)                             :: CCWF,FAC
      REAL (kind=KDBL) ,DIMENSION(1)               :: CNVPRCP, TOTPRCP, TPRCP, SRFLAG, SHDMIN, SHDMAX, CANOPY
      REAL (kind=KDBL) ,DIMENSION(1)               :: RAIN, RAINC
      REAL (kind=KDBL) ,DIMENSION(1)               :: ACV, ACVB, ACVT
      REAL (kind=KDBL) ,DIMENSION(2)               :: FLGMIN
      REAL (kind=KDBL) ,DIMENSION(3)               :: CRTRH
      REAL (kind=KDBL) ,DIMENSION(NUM_SOIL_LAYERS) :: SMC_V, STC_V, SLC_V
      REAL (kind=KDBL) ,DIMENSION(14)              :: HPRIME
      REAL (kind=KDBL) ,DIMENSION(:),ALLOCATABLE   :: UPD_MF, DWN_MF, DET_MF   !!!!!!!!!!! not in use
      REAL (kind=KDBL) ,DIMENSION(:),ALLOCATABLE   :: DQDT                     !!!!!!!!!!! not in use
      REAL (kind=KDBL) ,DIMENSION(:,:),ALLOCATABLE :: DQ3DT                    !!!!!!!!!!!  (9=5+pl_coeff)
      REAL (kind=KDBL) ,DIMENSION(:,:),ALLOCATABLE :: DT3DT                    !!!!!!!!!!! while
      REAL (kind=KDBL) ,DIMENSION(:,:),ALLOCATABLE :: DU3DT, DV3DT             !!!!!!!!!!! LDIAG3D =.FALSE.

      REAL (kind=KDBL) ,DIMENSION(:,:)  ,ALLOCATABLE :: OZPLOUT_V
      REAL (kind=KDBL) ,DIMENSION(:,:,:),ALLOCATABLE :: OZPLOUT

      REAL (kind=KDBL) ,DIMENSION(3)               :: PHY_F2DV   ! NUM_P2D for Zhao =3, Ferr=1 (fix later)
      REAL (kind=KDBL) ,DIMENSION(:,:),ALLOCATABLE :: PHY_F3DV   ! NUM_P3D for Zhao =4, Ferr=3 (fix later)
!--- gbphys output
      REAL (kind=KDBL) ,DIMENSION(1)               :: EVBSA, EVCWA, TRANSA, SBSNOA, SNOWCA, CLDWRK, PSMEAN
      REAL (kind=KDBL) ,DIMENSION(1)               :: CHH, CMM, EP, EPI, DLWSFCI, ULWSFCI, USWSFCI, DSWSFCI
      REAL (kind=KDBL) ,DIMENSION(1)               :: DLWSFC, ULWSFC, DTSFC, DQSFC, DUSFC, DVSFC, GFLUX
      REAL (kind=KDBL) ,DIMENSION(1)               :: DUSFCI, DVSFCI
      REAL (kind=KDBL) ,DIMENSION(1)               :: DTSFCI, DQSFCI, GFLUXI, T1, Q1, U1, V1
      REAL (kind=KDBL) ,DIMENSION(1)               :: ZLVL, SOILM, RUNOFF, SRUNOFF, SUNTIM
      REAL (kind=KDBL) ,DIMENSION(1)               :: F10M, UUSTAR, FFMM, FFHH, SPFHMIN, SPFHMAX
      REAL (kind=KDBL) ,DIMENSION(1)               :: PSURF, U10M, V10M, T2M, Q2M, HPBL, PWAT, SNOHFA
      REAL (kind=KDBL) ,DIMENSION(1)               :: DLWSFC_CC, ULWSFC_CC, DTSFC_CC, SWSFC_CC
      REAL (kind=KDBL) ,DIMENSION(1)               :: DUSFC_CC, DVSFC_CC, DQSFC_CC, PRECR_CC
      REAL (kind=KDBL) ,DIMENSION(1)               :: DUSFC_CPL,DVSFC_CPL,DTSFC_CPL,DQSFC_CPL
      REAL (kind=KDBL) ,DIMENSION(1)               :: DLWSFC_CPL,DSWSFC_CPL,DNIRBM_CPL,DNIRDF_CPL
      REAL (kind=KDBL) ,DIMENSION(1)               :: DVISBM_CPL,DVISDF_CPL,RAIN_CPL
      REAL (kind=KDBL) ,DIMENSION(1)               :: NLWSFC_CPL,NSWSFC_CPL,NNIRBM_CPL,NNIRDF_CPL
      REAL (kind=KDBL) ,DIMENSION(1)               :: NVISBM_CPL,NVISDF_CPL
      REAL (kind=KDBL) ,DIMENSION(1)               :: XT, XS, XU, XV, XZ, ZM, XTTS
      REAL (kind=KDBL) ,DIMENSION(1)               :: XZTS, D_CONV, IFD, DT_COOL, QRAIN
      REAL (kind=KDBL) ,DIMENSION(1)               :: SMCWLT2, SMCREF2, GSOIL, GTMP2M, GUSTAR, GPBLH, WET1
      REAL (kind=KDBL) ,DIMENSION(1)               :: GU10M, GV10M, GZORL, GORO, SR
      REAL (kind=KDBL) ,DIMENSION(1)               :: XMU_CC, DLW_CC, DSW_CC, SNW_CC, LPREC_CC, TREF
      REAL (kind=KDBL) ,DIMENSION(1)               :: DUSFCI_CPL,DVSFCI_CPL,DTSFCI_CPL,DQSFCI_CPL
      REAL (kind=KDBL) ,DIMENSION(1)               :: DLWSFCI_CPL,DSWSFCI_CPL,DNIRBMI_CPL,DNIRDFI_CPL
      REAL (kind=KDBL) ,DIMENSION(1)               :: DVISBMI_CPL,DVISDFI_CPL
      REAL (kind=KDBL) ,DIMENSION(1)               :: NLWSFCI_CPL,NSWSFCI_CPL,NNIRBMI_CPL,NNIRDFI_CPL
      REAL (kind=KDBL) ,DIMENSION(1)               :: NVISBMI_CPL,NVISDFI_CPL,T2MI_CPL,Q2MI_CPL
      REAL (kind=KDBL) ,DIMENSION(1)               :: U10MI_CPL,V10MI_CPL,TSEAI_CPL,PSURFI_CPL,ORO_CPL
      REAL (kind=KDBL) ,DIMENSION(1)               :: SLMSK_CPL
      REAL (kind=KDBL) ,DIMENSION(1)               :: Z_C, C_0, C_D, W_0, W_D, RQTK
      REAL (kind=KDBL) ,DIMENSION(1)               :: HLWD
      REAL (kind=KDBL) ,DIMENSION(LM)              :: DTDT
      REAL (kind=KDBL) ,DIMENSION(1)               :: TRIGGERPERTS
      LOGICAL, PARAMETER                           :: LSIDEA  = .FALSE.
!
#endif
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      btim0=timef()
!
!-----------------------------------------------------------------------
!
      RC_RUN=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  Extract the Solver internal state.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="SOLVER_RUN: Extract Solver Internal State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompGetInternalState(GRID_COMP                      &  !<-- The Solver component
                                        ,WRAP                           &
                                        ,RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      INT_STATE=>wrap%INT_STATE
!
!-----------------------------------------------------------------------
!***  The total number of forecast tasks.
!-----------------------------------------------------------------------
!
      INPES=int_state%INPES                                                !<-- I fcst tasks
      JNPES=int_state%JNPES                                                !<-- J fcst tasks
      NUM_PES=INPES*JNPES                                                  !<-- # of fcst tasks
!
!-----------------------------------------------------------------------
!***  Is this task on a domain boundary?
!-----------------------------------------------------------------------
!
      S_BDY=int_state%S_BDY
      N_BDY=int_state%N_BDY
      W_BDY=int_state%W_BDY
      E_BDY=int_state%E_BDY
!
!-----------------------------------------------------------------------
!***  Dereference fundamental variables for the dynamics routines.
!-----------------------------------------------------------------------
!
      ITS=int_state%ITS
      ITE=int_state%ITE
      JTS=int_state%JTS
      JTE=int_state%JTE
      IMS=int_state%IMS
      IME=int_state%IME
      JMS=int_state%JMS
      JME=int_state%JME
      IDS=int_state%IDS
      IDE=int_state%IDE
      JDS=int_state%JDS
      JDE=int_state%JDE
!
      ITS_B1=int_state%ITS_B1
      ITE_B1=int_state%ITE_B1
      ITS_B2=int_state%ITS_B2
      ITE_B2=int_state%ITE_B2
      ITS_B1_H1=int_state%ITS_B1_H1
      ITE_B1_H1=int_state%ITE_B1_H1
      ITE_B1_H2=int_state%ITE_B1_H2
      ITS_H1=int_state%ITS_H1
      ITE_H1=int_state%ITE_H1
      ITS_H2=int_state%ITS_H2
      ITE_H2=int_state%ITE_H2
      JTS_B1=int_state%JTS_B1
      JTE_B1=int_state%JTE_B1
      JTS_B2=int_state%JTS_B2
      JTE_B2=int_state%JTE_B2
      JTS_B1_H1=int_state%JTS_B1_H1
      JTE_B1_H1=int_state%JTE_B1_H1
      JTE_B1_H2=int_state%JTE_B1_H2
      JTS_H1=int_state%JTS_H1
      JTE_H1=int_state%JTE_H1
      JTS_H2=int_state%JTS_H2
      JTE_H2=int_state%JTE_H2
!
      LM=int_state%LM
!
      IHALO=int_state%IHALO    
      JHALO=int_state%JHALO    
!
      MYPE=int_state%MYPE                                                  !<-- The local task rank on this domain
      MY_DOMAIN_ID=int_state%MY_DOMAIN_ID
      MPI_COMM_COMP=int_state%MPI_COMM_COMP
!
!-----------------------------------------------------------------------
!***  Nested domains
!-----------------------------------------------------------------------
!
      I_AM_A_NEST=int_state%I_AM_A_NEST
!
!-----------------------------------------------------------------------
!***  Dereference more variables for shorter names.
!-----------------------------------------------------------------------
!
!     firstpass: IF(FIRST_PASS)THEN
!
      DDMPV=int_state%DDMPV
      DT=int_state%DT
      DYH=int_state%DYH
      DYV=int_state%DYV
      EF4T=int_state%EF4T
      GLOBAL=int_state%GLOBAL
!      HYDRO=int_state%HYDRO
      IDTADT=int_state%IDTADT
      IF(GLOBAL) THEN
        IDTADTQ=IDTADT  !global
      ELSE
        IDTADTQ=1       !regional
      ENDIF
      IHRSTBC=int_state%IHRSTBC
      KSE=int_state%NUM_TRACERS_MET
      KSS=1
      LNSAD=int_state%LNSAD
      LNSH=int_state%LNSH
      LNSV=int_state%LNSV
      LPT2=int_state%LPT2
      NBOCO=int_state%NBOCO
      NSTEPS_PER_CHECK=int_state%NSTEPS_PER_CHECK
      NSTEPS_PER_HOUR=int_state%NSTEPS_PER_HOUR
      NSTEPS_PER_RESET=int_state%NSTEPS_PER_RESET
      PDTOP=int_state%PDTOP
      PT=int_state%PT
      RDYH=int_state%RDYH
      RDYV=int_state%RDYV
      RUNBC=int_state%RUNBC
      SECADV=int_state%SECADV
      TBOCO=int_state%TBOCO
      FILTER_METHOD=int_state%FILTER_METHOD      
      FILTER_METHOD_LAST=int_state%FILTER_METHOD_LAST
!
      RIME_FACTOR_ADVECT=.FALSE.
      RIME_FACTOR_INPUT=.FALSE.
      IF (TRIM(int_state%MICROPHYSICS) == 'fer_hires' .AND.         &
          int_state%F_QG .AND. int_state%SPEC_ADV) THEN
         RIME_FACTOR_ADVECT=.TRUE.
      ENDIF
!
      PARENT_CHILD_TIME_RATIO=int_state%PARENT_CHILD_TIME_RATIO
!
      DO N=1,3
        IDATBC(N)=int_state%IDATBC(N)
      ENDDO
!
      CALL SET_DOMAIN_SPECS(int_state%ITS,int_state%ITE                 &          
                           ,int_state%JTS,int_state%JTE                 &
                           ,int_state%IMS,int_state%IME                 &
                           ,int_state%JMS,int_state%JME                 &
                           ,int_state%IDS,int_state%IDE                 &
                           ,int_state%JDS,int_state%JDE                 &
                           ,int_state%IHALO,int_state%JHALO             &
                           ,int_state%MY_DOMAIN_ID                      &
                           ,int_state%MYPE                              &
                           ,int_state%MY_NEB                            &
                           ,int_state%MPI_COMM_COMP                     &
                           ,int_state%NUM_PES                           &
!
                           ,LOCAL_ISTART_IN=int_state%LOCAL_ISTART      &
                           ,LOCAL_IEND_IN=int_state%LOCAL_IEND          &
                           ,LOCAL_JSTART_IN=int_state%LOCAL_JSTART      &
                           ,LOCAL_JEND_IN=int_state%LOCAL_JEND          &
                           ,ADV_STANDARD_IN=int_state%ADV_STANDARD      &
                           ,ADV_UPSTREAM_IN=int_state%ADV_UPSTREAM      &
                           ,S_BDY_IN=int_state%S_BDY                    &
                           ,N_BDY_IN=int_state%N_BDY                    &
                           ,W_BDY_IN=int_state%W_BDY                    &
                           ,E_BDY_IN=int_state%E_BDY                    &
                             )
!
!-----------------------------------------------------------------------
!***  Extract the timestep count from the Clock.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Solver Run Gets Timestep from the ATM Clock"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ClockGet(clock       =CLOCK_ATM                         &  !<-- The ESMF Clock
                        ,timeStep    =DT_ESMF                           &  !<-- Fundamental timestep (s) (ESMF)
                        ,currtime    =CURRTIME                          &  !<-- current time
                        ,advanceCount=NTIMESTEP_ESMF                    &  !<-- The number of times the clock has been advanced
                        ,rc          =RC)

!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_TimeIntervalGet(timeinterval=DT_ESMF                    &  !<-- the ESMF timestep
                               ,s           =INTEGER_DT                 &  !<-- the integer part of the timestep in seconds
                               ,sN          =NUMERATOR_DT               &  !<-- the numerator of the fractional second
                               ,sD          =IDENOMINATOR_DT            &  !<-- the denominator of the fractional second
                               ,rc          =RC)
!
      int_state%DT=REAL(INTEGER_DT)+REAL(NUMERATOR_DT)                  &  !<-- Fundamental timestep (s) (REAL)
                                   /REAL(IDENOMINATOR_DT)
      DT=int_state%DT
!
      CALL ESMF_AttributeGet(state=IMP_STATE                            &
                            ,name ='FUND_DT'                            &
                            ,value=FUND_DT                              &
                            ,rc   =RC)
!
      DTRATIO=ABS(DT/FUND_DT)
!
      NTIMESTEP=NTIMESTEP_ESMF
      int_state%NTSD=NTIMESTEP

!----------------------
      if (int_state%RADAR_INIT==0) then
          USE_RADAR=.false.
      else
         if(filter_method==0) then
            USE_RADAR=.false.
         else 
            USE_RADAR=.true.
         end if
      endif
!----------------------
!     write(6,*)'filter method::',NTIMESTEP,int_state%RADAR_INIT,FILTER_METHOD,FILTER_METHOD_LAST,USE_RADAR

        if (DT .lt. 0 .and. FILTER_METHOD .ge. 2 .and. &
           (int(NTIMESTEP*DT) .le. int_state%DFIHR_BOCO/2.)) then
      HYDRO=.true.
        else
      HYDRO=int_state%HYDRO
        endif



!     
      FIRST_PASS=int_state%FIRST_PASS
!
      NSTEPS_PER_HOUR=NINT(3600./DT)
!
      N_PRINT_STATS=NINT(3600./DT)                                         !<-- Print layer statistics once per forecast hour
!
!-----------------------------------------------------------------------
!***  Extract the horizontal diffusion flag from the import state.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Solver Run Extracts Horizontal Diffusion Flag "
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=IMP_STATE                            &
                            ,name ='HDIFF'                              &
                            ,value=HDIFF_ON                             &
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Extract the digital filter method from the import state.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!     MESSAGE_CHECK="Solver Run Extracts Horizontal Diffusion Flag "
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!     CALL ESMF_AttributeGet(state=IMP_STATE                            &  !<-- The Solver import state
!                           ,name ='Filter_Method'                      &  !<-- Name of the attribute to extract
!                           ,value=int_state%FILTER_METHOD              &  !<-- The scalar being extracted from the import state
!                           ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!     FILTER_METHOD=int_state%FILTER_METHOD      
!     FILTER_METHOD_LAST=int_state%FILTER_METHOD_LAST
!
!-----------------------------------------------------------------------
!rh_hold        IF (USE_RADAR_FIRST == 1 .and. FIRST_PASS ) THEN
!rh_hold           ALLOCATE(RH_HOLD(IMS:IME,JMS:JME,1:LM))
!rh_hold           IFLAG=1 ! <----   IFLAG=1 takes T,Q,P and returns RH_HOLD
!rh_hold           CALL CALC_RH_RADAR_DFI(int_state%T,int_state%Q,int_state%PD  &
!rh_hold                                  ,int_state%PSGML1,int_state%SGML2     &
!rh_hold                                  ,R_D,R_V,RH_HOLD                      & 
!rh_hold                                  ,IMS,IME,JMS,JME,LM                   & 
!rh_hold                                  ,IFLAG)
!rh_hold        ENDIF

!
!     ENDIF firstpass
!
!-----------------------------------------------------------------------
!***  The following set of internal state arrays never changes unless
!***  the domain moves in which case they must be dereferenced again.
!-----------------------------------------------------------------------
!
      MOVE_NOW=.FALSE.
      IF(int_state%MY_DOMAIN_MOVES)THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract the MOVE_NOW flag in SOLVER_RUN"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state=IMP_STATE                          &  !<-- Solver import state
                              ,name ='MOVE_NOW'                         &  !<-- Name of the flag for current domain motion
                              ,value=MOVE_NOW                           &  !<-- Did the nest move this timestep?
                              ,rc   =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF

!
!-----------------------------------------------------------------------
!***  If this is a moving nest and it moved this timestep then we
!***  need to update the haloes of the geographic lat/lon and the
!***  HDAC variables because like all variables they are updated
!***  only in the integration region when a nest shifts.
!-----------------------------------------------------------------------
!
      IF(MOVE_NOW)THEN
!
        btim=timef()
        CALL HALO_EXCH                                                  &
           (int_state%GLAT,1                                            &
           ,int_state%GLON,1                                            &
           ,int_state%VLAT,1                                            &
           ,int_state%VLAT,1                                            &
           ,2,2)
!
        CALL HALO_EXCH                                                  &
           (int_state%HDACX,1                                           &
           ,int_state%HDACY,1                                           &
           ,int_state%HDACVX,1                                          &
           ,int_state%HDACVY,1                                          &
           ,2,2)
        td%exch_dyn=td%exch_dyn+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Also the geography information for the gravity wave drag
!***  must be updated to account for the domain's new position.
!
!***  NOTE:  Currently the gravity wave drag is turned off in
!***         moving nests.  A quantity used by the parameterization
!***         is mountains' angle with respect to east.  From the
!***         moving nest's perspective the mountains are moving
!***         and thus those angles would need to be updated.
!***         Such updating is not yet included.
!-----------------------------------------------------------------------
!
        IF(int_state%GWDFLG)THEN
!
          DTPHY=int_state%DT*int_state%NPHS
!
          CALL GWD_init(DTPHY,int_state%RESTART                         &
                       ,int_state%CLEFFAMP,int_state%DPHD               &
                       ,int_state%CLEFF                                 &
                       ,int_state%TPH0D,int_state%TLM0D                 &
                       ,int_state%GLAT,int_state%GLON                   &
                       ,int_state%CROT,int_state%SROT,int_state%HANGL   &
                       ,IDS,IDE,JDS,JDE                                 &
                       ,IMS,IME,JMS,JME                                 &
                       ,ITS,ITE,JTS,JTE,LM)
        ENDIF
!
        IF(int_state%NTRACK>0 .AND. int_state%MYPE<int_state%NUM_PES) THEN
           CALL UPDATE_TRACKER_POST_MOVE(INT_STATE)
        ENDIF
!
!-----------------------------------------------------------------------
!
      ENDIF
!
!-----------------------------------------------------------------------
!
      IF (INTEGER_DT >= 0) IFACT=1
      IF (INTEGER_DT <  0) IFACT=-1
!
      IF(FIRST_PASS)THEN
!
        IF (FILTER_METHOD /= 0) THEN
!
!*** Save copies of the internal state variables scaled by
!*** DTRATIO below so we can restore them precisely after the filter 
!*** (*_SV variables).  Needed for bit-wise identical restarts
!
          ALLOCATE(HDACX_SV(ITS:ITE,JTS:JTE),HDACY_SV(ITS:ITE,JTS:JTE),   &
                   HDACVX_SV(ITS:ITE,JTS:JTE),HDACVY_SV(ITS:ITE,JTS:JTE))
          ALLOCATE(DDMPU_SV(JDS:JDE),FAD_SV(JDS:JDE),FAH_SV(JDS:JDE),     &
                   FCP_SV(JDS:JDE))
 
          DDMPV_SV=int_state%DDMPV
          EF4T_SV=int_state%EF4T
        END IF

        int_state%DDMPV=IFACT*DTRATIO*int_state%DDMPV
        int_state%EF4T=IFACT*DTRATIO*int_state%EF4T
!
        DO J=JDS,JDE

          IF (FILTER_METHOD /= 0) THEN
             DDMPU_SV(J) = int_state%DDMPU(J)
             FAD_SV(J)   = int_state%FAD(J)
             FAH_SV(J)   = int_state%FAH(J)
             FCP_SV(J)   = int_state%FCP(J)
          END IF

          int_state%DDMPU(J)=IFACT*int_state%DDMPU(J)
          int_state%FAD(J)=IFACT*DTRATIO*int_state%FAD(J)
          int_state%FAH(J)=IFACT*DTRATIO*int_state%FAH(J)
          int_state%FCP(J)=IFACT*DTRATIO*int_state%FCP(J)
          int_state%WPDAR(J)=IFACT*int_state%WPDAR(J)
        ENDDO
!
        DO J=JTS,JTE
        DO I=ITS,ITE

         IF (FILTER_METHOD /= 0) THEN
            HDACX_SV(I,J)=int_state%HDACX(I,J)
            HDACY_SV(I,J)=int_state%HDACY(I,J)
            HDACVX_SV(I,J)=int_state%HDACVX(I,J)
            HDACVY_SV(I,J)=int_state%HDACVY(I,J)
         END IF

          int_state%HDACX(I,J)=IFACT*DTRATIO*int_state%HDACX(I,J)
          int_state%HDACY(I,J)=IFACT*DTRATIO*int_state%HDACY(I,J)
          int_state%HDACVX(I,J)=IFACT*DTRATIO*int_state%HDACVX(I,J)
          int_state%HDACVY(I,J)=IFACT*DTRATIO*int_state%HDACVY(I,J)
        ENDDO
        ENDDO
      ENDIF
!
      DDMPV=int_state%DDMPV
      EF4T=int_state%EF4T
!
      NBOCO=int(0.5+NBOCO/DTRATIO)
!     IF (MYPE == 0) WRITE(0,*) 'NBOCO reset to : ', NBOCO
!
!-----------------------------------------------------------------------
!***  Now we need to do some things related to digital filtering
!***  that are only relevant after the first pass through the
!***  Run step.
!-----------------------------------------------------------------------
!
      DT_TEST=INTEGER_DT
      DT_TEST_RATIO=int_state%DT_TEST_RATIO
!
!-----------------------------------------------------------------------
!
      not_firstpass: IF (.NOT. FIRST_PASS) THEN
!
!-----------------------------------------------------------------------
! 
        changedir: IF (int_state%DT_LAST /= DT_TEST                     &
                                 .AND.                                  &
                       ABS(int_state%DT_LAST) == ABS(DT_TEST) ) THEN
!
!-----------------------------------------------------------------------
!
          IF(MYPE == 0)WRITE(0,*)' Change in integration direction...'  &
                                ,' dt_last=',int_state%dt_last          &
                                ,' dt_test=',dt_test
!
!-----------------------------------------------------------------------
!***  Setting previous time level variables (Adams-Bashforth scheme)
!***  to the current time level.  Seems safer than potentially leaving them
!***  defined as values at a very different point in the time integration.
!-----------------------------------------------------------------------
!
          int_state%FIRST_STEP=.TRUE.
!
          int_state%TP=int_state%T
          int_state%UP=int_state%U
          int_state%VP=int_state%V
!
          IFACT=-1
!
          int_state%DDMPV=IFACT*int_state%DDMPV
          int_state%EF4T=IFACT*int_state%EF4T
          DDMPV=int_state%DDMPV
          EF4T=int_state%EF4T
!
          DO J=JDS,JDE
            int_state%DDMPU(J)=IFACT*int_state%DDMPU(J)
            int_state%FAD(J)=IFACT*int_state%FAD(J)
            int_state%FAH(J)=IFACT*int_state%FAH(J)
            int_state%FCP(J)=IFACT*int_state%FCP(J)
            int_state%WPDAR(J)=IFACT*int_state%WPDAR(J)
          ENDDO
!
          DO J=JTS,JTE
          DO I=ITS,ITE
            int_state%HDACX(I,J)=IFACT*int_state%HDACX(I,J)
            int_state%HDACY(I,J)=IFACT*int_state%HDACY(I,J)
            int_state%HDACVX(I,J)=IFACT*int_state%HDACVX(I,J)
            int_state%HDACVY(I,J)=IFACT*int_state%HDACVY(I,J)
          ENDDO
          ENDDO
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Solver Run Gets HDIFF from Import State"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeGet(state=IMP_STATE                        &  !<-- The Solver import state
                                ,name ='HDIFF'                          &  !<-- Name of the Attribute to extract
                                ,value=HDIFF_ON                         &  !<-- Put the Attribute here
                                ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          btim=timef()
          CALL HALO_EXCH                                                &
             (int_state%T,LM                                            &
             ,int_state%Q,LM                                            &
             ,int_state%CW,LM                                           &
             ,2,2)
!..What about other items in the TRACER array?
!
          CALL HALO_EXCH                                                &
             (int_state%U,LM                                            &
             ,int_state%V,LM                                            &
             ,2,2)
!
          CALL HALO_EXCH                                                &
             (int_state%PD,1                                            &
             ,2,2)
          td%exch_dyn=td%exch_dyn+(timef()-btim)
!
          IF(.NOT.int_state%GLOBAL)THEN
            CALL WRITE_BC(LM,LNSH,LNSV,NTIMESTEP,DT                     &
                         ,RUNBC                                         &
                         ,TBOCO+int_state%DFIHR_BOCO/2.                 &
                         ,int_state%NVARS_BC_2D_H                       &
                         ,int_state%NVARS_BC_3D_H                       &
                         ,int_state%NVARS_BC_4D_H                       &
                         ,int_state%NVARS_BC_2D_V                       &
                         ,int_state%NVARS_BC_3D_V                       &
                         ,int_state%BND_VARS_H                          &
                         ,int_state%BND_VARS_V                          &
                         ,.TRUE.)                                          !<-- Recompute tendencies at this stage?
          ENDIF
!
!-----------------------------------------------------------------------
!
        ENDIF changedir
!
!-----------------------------------------------------------------------
!
        end_filt: IF (FILTER_METHOD /= FILTER_METHOD_LAST) THEN
!
!-----------------------------------------------------------------------
!
          DTRATIO=ABS(FUND_DT/DT_TEST_RATIO)
          IF(MYPE == 0) WRITE(0,*) ' RESTORING PRE-FILTER CONSTANTS with DTRATIO: ', DTRATIO
!
!-----------------------------------------------------------------------
!***  Setting previous time level variables (Adams-Bashforth scheme)
!***  to the current time level.  Seems safer than potentially leaving them
!***  defined as values at a very different point in the time integration.
!-----------------------------------------------------------------------
!
          int_state%TP=int_state%T
          int_state%UP=int_state%U
          int_state%VP=int_state%V
!
          IFACT=1
!
          int_state%DDMPV=DDMPV_SV
          int_state%EF4T=EF4T_SV

          DDMPV=int_state%DDMPV
          EF4T=int_state%EF4T

          DDMPV=int_state%DDMPV
          EF4T=int_state%EF4T
          NBOCO=int(0.5+NBOCO/DTRATIO)
!
!         IF (MYPE == 0) WRITE(0,*) 'NBOCO reset to : ', NBOCO
!
          DO J=JDS,JDE
            int_state%DDMPU(J)=DDMPU_SV(J)
            int_state%FAD(J)=FAD_SV(J)
            int_state%FAH(J)=FAH_SV(J)
            int_state%FCP(J)=FCP_SV(J)
            int_state%WPDAR(J)=IFACT*int_state%WPDAR(J)
          ENDDO
!
          DO J=JTS,JTE
          DO I=ITS,ITE
            int_state%HDACX(I,J)=HDACX_SV(I,J)
            int_state%HDACY(I,J)=HDACY_SV(I,J)
            int_state%HDACVX(I,J)=HDACVX_SV(I,J)
            int_state%HDACVY(I,J)=HDACVY_SV(I,J)
          ENDDO
          ENDDO
!
!-----------------------------------------------------------------------
!
          DEALLOCATE(HDACX_SV,HDACY_SV,HDACVX_SV,HDACVY_SV)
          DEALLOCATE(DDMPU_SV,FAD_SV,FAH_SV,FCP_SV)
!
          int_state%FIRST_STEP=.TRUE.
!
        ENDIF end_filt
!
!-----------------------------------------------------------------------
!
      ENDIF not_firstpass
!
!-----------------------------------------------------------------------
!
      IF(FIRST_PASS)THEN
        int_state%FIRST_PASS=.FALSE.
        FIRST_PASS=int_state%FIRST_PASS
      ENDIF
!
!-----------------------------------------------------------------------
!
      TD=>TIMERS(MY_DOMAIN_ID)                                             !<-- Abbreviate the name of this domain's timers.
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  Begin the Solver calling sequence.
!***  Note that the first timestep begins differently
!***  than all subsequent timesteps.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      firststep: IF(int_state%FIRST_STEP.AND.                           &  !<--  The following block is used only for
                    .NOT.int_state%RESTART)THEN                            !     the first timestep and cold start
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
        IF(GLOBAL)THEN
!
          btim=timef()
          CALL SWAPHN                                                   &
           (int_state%T,IMS,IME,JMS,JME,LM                              &
           ,INPES)
          td%swaphn_tim=td%swaphn_tim+(timef()-btim)
!
          btim=timef()
          CALL POLEHN                                                   &
           (int_state%T                                                 &
           ,IMS,IME,JMS,JME,LM                                          &
           ,INPES,JNPES)
          td%polehn_tim=td%polehn_tim+(timef()-btim)
!
        ENDIF
!
        btim=timef()
        CALL HALO_EXCH(int_state%T,LM                                   &
                      ,2,2)
        td%exch_dyn=td%exch_dyn+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  The pressure gradient routine.
!-----------------------------------------------------------------------
!
        btim=timef()
!
        CALL PGFORCE                                                    &
          (int_state%FIRST_STEP,int_state%GLOBAL,int_state%RESTART      &
          ,LM,DT,NTIMESTEP                                              &
          ,RDYV,int_state%DSG2,int_state%PDSG1,int_state%RDXV           &
          ,int_state%WPDAR,int_state%FIS                                &
          ,int_state%PD                                                 &
          ,int_state%T,int_state%Q,int_state%CW                         & ! And how about other TRACER elements?
          ,int_state%PINT                                               &
          ,int_state%RTOP                                               &
          ,int_state%DIV                                                &
          ,int_state%PCNE,int_state%PCNW                                &
          ,int_state%PCX,int_state%PCY                                  &
          ,int_state%TCU,int_state%TCV )
!
        td%pgforce_tim=td%pgforce_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!
        btim=timef()
        CALL HALO_EXCH(int_state%DIV,LM                                 &
                      ,2,2)
        CALL HALO_EXCH(int_state%U,LM                                   &
                      ,int_state%V,LM                                   &
                      ,2,2)
        td%exch_dyn=td%exch_dyn+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Divergence and horizontal pressure advection in thermo eqn
!-----------------------------------------------------------------------
!
        btim=timef()
!
        CALL DHT                                                        &
          (GLOBAL,LM,DYV,int_state%DSG2,int_state%PDSG1,int_state%DXV   &
          ,int_state%FCP,int_state%FDIV                                 &
          ,int_state%PD,int_state%PDO                                   &
          ,int_state%U,int_state%V                                      &
          ,int_state%OMGALF                                             &
          ,int_state%PCNE,int_state%PCNW,int_state%PCX,int_state%PCY    &
          ,int_state%PFNE,int_state%PFNW,int_state%PFX,int_state%PFY    &
          ,int_state%DIV,int_state%TDIV)
!
        td%dht_tim=td%dht_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Filtering and boundary conditions for the global forecast.
!-----------------------------------------------------------------------
!
        IF(GLOBAL)THEN
!
          btim=timef()
          CALL FFTFHN                                                   &
           (LM                                                          &
           ,int_state%KHFILT                                            &
           ,int_state%HFILT                                             &
           ,int_state%DIV                                               &
           ,int_state%WFFTRH,int_state%NFFTRH                           &
           ,NUM_PES,MYPE,MPI_COMM_COMP)
          td%fftfhn_tim=td%fftfhn_tim+(timef()-btim)
!
          btim=timef()
          CALL SWAPHN                                                   &
           (int_state%DIV                                               &
           ,IMS,IME,JMS,JME,LM                                          &
           ,INPES)
!
          CALL SWAPHN                                                   &
           (int_state%OMGALF                                            &
           ,IMS,IME,JMS,JME,LM                                          &
           ,INPES)
          td%swaphn_tim=td%swaphn_tim+(timef()-btim)
!
          btim=timef()
          CALL POLEHN                                                   &
           (int_state%DIV                                               &
           ,IMS,IME,JMS,JME,LM                                          &
           ,INPES,JNPES)
!
          CALL POLEHN                                                   &
           (int_state%OMGALF                                            &
           ,IMS,IME,JMS,JME,LM                                          &
           ,INPES,JNPES)
          td%polehn_tim=td%polehn_tim+(timef()-btim)
!
          btim=timef()
          CALL SWAPWN                                                   &
            (int_state%U                                                &
            ,IMS,IME,JMS,JME,LM                                         &
            ,INPES)
!
          CALL SWAPWN                                                   &
            (int_state%V                                                &
            ,IMS,IME,JMS,JME,LM                                         &
            ,INPES)
          td%swapwn_tim=td%swapwn_tim+(timef()-btim)
!
          btim=timef()
          CALL POLEWN                                                   &
            (int_state%U,int_state%V                                    &
            ,IMS,IME,JMS,JME,LM                                         &
            ,INPES,JNPES)
          td%polewn_tim=td%polewn_tim+(timef()-btim)
!
        ENDIF
!
!-----------------------------------------------------------------------
!
        btim=timef()
        CALL HALO_EXCH                                                  &
         (int_state%T,LM                                                &
         ,int_state%U,LM                                                &
         ,int_state%V,LM                                                &
         ,2,2)
        td%exch_dyn=td%exch_dyn+(timef()-btim)
!
!-----------------------------------------------------------------------
!
      ENDIF firststep
!
!-----------------------------------------------------------------------
!
      not_firststep: IF(.NOT.int_state%FIRST_STEP                       &  !<-- The following block is for all timesteps after
                        .OR.int_state%RESTART)THEN                         !    the first or all steps in restart case
!

!rh_hold        IF(FILTER_METHOD==0.and.USE_RADAR_FIRST==1.and.USE_RADAR==0)THEN
!rh_hold            IFLAG=-1  !  <----   IFLAG=-1 takes RH_HOLD, and filtered
!rh_hold                      !          T,P, to restore Q to be consistent with
!rh_hold                      !          prefiltered humidity level
!rh_hold
!rh_hold!!! NOTE:  restoring down here means they are restored AFTER the 00 h
!rh_hold!!!        output is written.  Any way to restore them before the output
!rh_hold!!!        is written?
!rh_hold
!rh_hold            CALL CALC_RH_RADAR_DFI(int_state%T,int_state%Q,int_state%PD &
!rh_hold                                  ,int_state%PSGML1,int_state%SGML2     &
!rh_hold                                  ,R_D,R_V,RH_HOLD                      &
!rh_hold                                  ,IMS,IME,JMS,JME,LM                   &
!rh_hold                                  ,IFLAG)
!rh_hold
!rh_hold        USE_RADAR_FIRST=0
!rh_hold
!rh_hold        ENDIF

!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Horizontal diffusion (internal halo exchange for 4th order)
!-----------------------------------------------------------------------
!
        btim=timef()
!
        IF(HDIFF_ON>0)THEN
          CALL HDIFF                                                    &
            (GLOBAL,HYDRO                                               &
            ,INPES,JNPES,LM,LPT2                                        &
            ,DYH,RDYH                                                   &
            ,int_state%EPSQ2                                            &
            ,int_state%DXV,int_state%RARE,int_state%RDXH                &
            ,int_state%SICE,int_state%SM                                &
            ,int_state%HDACX,int_state%HDACY                            &
            ,int_state%HDACVX,int_state%HDACVY                          &
            ,int_state%W,int_state%Z                                    &
            ,int_state%CW,int_state%Q,int_state%Q2                      & ! And how about other TRACER elements?
            ,int_state%T,int_state%U,int_state%V,int_state%DEF)
        ENDIF
!
        td%hdiff_tim=td%hdiff_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Filtering and boundary conditions for the global forecast.
!-----------------------------------------------------------------------
!
        IF(GLOBAL)THEN
!
          btim=timef()
!
          CALL POAVHN                                                   &
            (IMS,IME,JMS,JME,LM                                         &
            ,int_state%T                                                &
            ,INPES,JNPES                                                &
            ,int_state%USE_ALLREDUCE)
!
          CALL POAVHN                                                   &
            (IMS,IME,JMS,JME,LM                                         &
            ,int_state%Q                                                &
            ,INPES,JNPES                                                &
            ,int_state%USE_ALLREDUCE)
!
          CALL POAVHN                                                   &
            (IMS,IME,JMS,JME,LM                                         &
            ,int_state%CW                                               &
            ,INPES,JNPES                                                &
            ,int_state%USE_ALLREDUCE)
!
          CALL POAVHN                                                   &
            (IMS,IME,JMS,JME,LM                                         &
            ,int_state%Q2                                               &
            ,INPES,JNPES                                                &
            ,int_state%USE_ALLREDUCE)
!
          td%poavhn_tim=td%poavhn_tim+(timef()-btim)
!
          btim=timef()
          CALL SWAPHN(int_state%T,IMS,IME,JMS,JME,LM,INPES)
          CALL SWAPHN(int_state%Q,IMS,IME,JMS,JME,LM,INPES)
          CALL SWAPHN(int_state%CW,IMS,IME,JMS,JME,LM,INPES)
          CALL SWAPHN(int_state%Q2,IMS,IME,JMS,JME,LM,INPES)
          td%swaphn_tim=td%swaphn_tim+(timef()-btim)
!
          btim=timef()
          CALL POLEHN(int_state%T,IMS,IME,JMS,JME,LM,INPES,JNPES)
          CALL POLEHN(int_state%Q,IMS,IME,JMS,JME,LM,INPES,JNPES)
          CALL POLEHN(int_state%CW,IMS,IME,JMS,JME,LM,INPES,JNPES)
          CALL POLEHN(int_state%Q2,IMS,IME,JMS,JME,LM,INPES,JNPES)
          td%polehn_tim=td%polehn_tim+(timef()-btim)
!
          btim=timef()
          CALL SWAPWN(int_state%U,IMS,IME,JMS,JME,LM,INPES)
          CALL SWAPWN(int_state%V,IMS,IME,JMS,JME,LM,INPES)
          td%swapwn_tim=td%swapwn_tim+(timef()-btim)
!
          btim=timef()
          CALL POLEWN(int_state%U,int_state%V                           &
                     ,IMS,IME,JMS,JME,LM,INPES,JNPES)
          td%polewn_tim=td%polewn_tim+(timef()-btim)
!
        ENDIF
!
!-----------------------------------------------------------------------
!
        btim=timef()
        CALL HALO_EXCH(int_state%T,LM                                   &
                      ,int_state%Q,LM                                   &
                      ,int_state%CW,LM                                  &
                      ,int_state%Q2,LM                                  & ! And how about other TRACER elements?
                      ,2,2)
        CALL HALO_EXCH(int_state%U,LM                                   &
                      ,int_state%V,LM                                   &
                      ,1,1)
        td%exch_dyn=td%exch_dyn+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Regional domains that have no children or are uppermost parents
!***  need to set a digital filter flag and exchange haloes.
!-----------------------------------------------------------------------
!
        IF(.NOT.I_AM_A_NEST.AND..NOT.GLOBAL)THEN                           !<-- For single domains or uppermost parents
!
          READBC=(NTIMESTEP==1.OR.MOD(NTIMESTEP,NBOCO)==0)
!
          bc_check: IF(READBC)THEN                                         !<-- Is it time to read BCs?
!
            IF(MYPE==0)THEN
              WRITE_BC_FLAG=0
!
              IF(FILTER_METHOD>0)THEN
                IF(NTIMESTEP<=1.AND.int_state%BDY_WAS_READ)THEN
                  WRITE_BC_FLAG=1
                ELSE
                  WRITE_BC_FLAG=0
                ENDIF
              ENDIF
!
            ENDIF
!
            CALL MPI_BCAST(WRITE_BC_FLAG,1,MPI_INTEGER,0                &
                          ,MPI_COMM_COMP,IRTN)
!
            IF(WRITE_BC_FLAG==1)THEN
              CALL HALO_EXCH                                            &
               (int_state%T,LM                                          &
               ,int_state%Q,LM                                          &
               ,int_state%CW,LM                                         & ! And how about other TRACER elements?
               ,2,2)
!
              CALL HALO_EXCH                                            &
               (int_state%U,LM                                          &
               ,int_state%V,LM                                          &
               ,2,2)
!
             CALL HALO_EXCH                                             &
              (int_state%PD,1                                           &
              ,2,2)
!
            ENDIF
!
          ENDIF  bc_check
!
        ENDIF
!
!-----------------------------------------------------------------------
!***  Update the boundary mass points.
!
!***  For non-nested regional domains, read new boundary tendencies
!***  at the appropriate times.
!
!***  If this is a nested domain then unload the new boundary data
!***  from the Solver import state and compute the time tendencies.
!-----------------------------------------------------------------------
!
        bc_update: IF(.NOT.GLOBAL)THEN
!
!-----------------------------------------------------------------------
!***  The following block is for digital filtering.
!-----------------------------------------------------------------------
!
          IF(I_AM_A_NEST)THEN
!
            IF(MYPE==0)THEN
              WRITE_BC_FLAG_NEST=0
!
              IF(FILTER_METHOD>0)THEN
                IF (S_BDY.AND.W_BDY                                     &
                         .AND.                                          &
                    NTIMESTEP <= 1                                      &
                         .AND.                                          &
                    int_state%BDY_WAS_READ) THEN
!
                  WRITE_BC_FLAG_NEST=1
                ENDIF
              ENDIF
!
            ENDIF
!
            CALL MPI_BCAST(WRITE_BC_FLAG_NEST,1,MPI_INTEGER             &
                          ,0,MPI_COMM_COMP,IRTN)
!
            IF (WRITE_BC_FLAG_NEST == 1) THEN
              CALL HALO_EXCH                                            &
               (int_state%T,LM                                          &
               ,int_state%Q,LM                                          &
               ,int_state%CW,LM                                         & ! And how about other TRACER elements?
               ,2,2)
!
              CALL HALO_EXCH                                            &
               (int_state%U,LM                                          &
               ,int_state%V,LM                                          &
               ,2,2)
!
              CALL HALO_EXCH                                            &
               (int_state%PD,1                                          &
               ,2,2)
            ENDIF
          ENDIF
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Set SIMULATION_START_TIME for Filter in Solver Run"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_TimeSet(time=SIMULATION_START_TIME                  &
                           ,yy  =START_YEAR                             &
                           ,mm  =START_MONTH                            &
                           ,dd  =START_DAY                              &
                           ,h   =START_HOUR                             &
                           ,calkindflag=ESMF_CALKIND_GREGORIAN          &
                           ,rc  =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          IF (FILTER_METHOD == 1 .and. NTIMESTEP == 0) THEN
!
            REST_OFFSET=CURRTIME-SIMULATION_START_TIME
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Get Time Offset for Filter in Solver Run"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_TimeIntervalGet(timeinterval=REST_OFFSET          &
                                     ,s           =JDAT(7)              &
                                     ,rc          =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            RESTVAL=JDAT(7)
            IF (MYPE == 0) WRITE(0,*) 'set RESTVAL to: ', RESTVAL
!
          ENDIF
!
!-----------------------------------------------------------------------
!
          boundary_tendencies: IF(S_BDY.OR.N_BDY.OR.W_BDY.OR.E_BDY)THEN
!
!-----------------------------------------------------------------------
!***  Nests update boundary tendencies based on data from parent.
!-----------------------------------------------------------------------
!
            nest_or_parent: IF(I_AM_A_NEST)THEN
!
!-----------------------------------------------------------------------
!***  The following block is for digital filtering.
!-----------------------------------------------------------------------
!
              IF(NTIMESTEP<=1.AND.WRITE_BC_FLAG_NEST==1)THEN
!
                TBOCO=PARENT_CHILD_TIME_RATIO*DT
                CALL WRITE_BC(LM,LNSH,LNSV,NTIMESTEP,DT                 &
                             ,RUNBC,TBOCO                               &
                             ,int_state%NVARS_BC_2D_H                   &
                             ,int_state%NVARS_BC_3D_H                   &
                             ,int_state%NVARS_BC_4D_H                   &
                             ,int_state%NVARS_BC_2D_V                   &
                             ,int_state%NVARS_BC_3D_V                   &
                             ,int_state%BND_VARS_H                      &
                             ,int_state%BND_VARS_V                      &
                             ,.FALSE.)                                     !<-- Are tendencies recomputed?
!
              ENDIF
!
!-----------------------------------------------------------------------
!
              COMPUTE_BC=(NTIMESTEP==1.OR.                              &
                          MOD(NTIMESTEP,PARENT_CHILD_TIME_RATIO)==0)
!
              IF(COMPUTE_BC)THEN
!   
                CALL UPDATE_BC_TENDS(IMP_STATE                          &
                                    ,LM,LNSH,LNSV                       &
                                    ,PARENT_CHILD_TIME_RATIO,DT         &
                                    ,S_BDY,N_BDY,W_BDY,E_BDY            &
                                    ,int_state%NLEV_H                   &
                                    ,int_state%NLEV_V                   &
                                    ,int_state%NVARS_BC_2D_H            &
                                    ,int_state%NVARS_BC_3D_H            &
                                    ,int_state%NVARS_BC_4D_H            &
                                    ,int_state%NVARS_BC_2D_V            &
                                    ,int_state%NVARS_BC_3D_V            &
                                    ,int_state%BND_VARS_H               &
                                    ,int_state%BND_VARS_V               &
                                    ,int_state%ITS,int_state%ITE        &
                                    ,int_state%JTS,int_state%JTE        &
                                    ,int_state%IMS,int_state%IME        &
                                    ,int_state%JMS,int_state%JME        &
                                    ,int_state%IDS,int_state%IDE        &
                                    ,int_state%JDS,int_state%JDE        &
                                                                 )
!
              ENDIF
!
!-----------------------------------------------------------------------
!***  Single/uppermost domain reads its own boundary input data
!-----------------------------------------------------------------------
!
            ELSE nest_or_parent
!
              CALL ESMF_TimeSet(time=SIMULATION_START_TIME              &
                               ,yy  =START_YEAR                         &
                               ,mm  =START_MONTH                        &
                               ,dd  =START_DAY                          &
                               ,h   =START_HOUR)
!
              IF (FILTER_METHOD > 0 .and. NTIMESTEP == 0) THEN
                REST_OFFSET=CURRTIME-SIMULATION_START_TIME
                CALL ESMF_TimeIntervalGet(timeinterval=REST_OFFSET, s=JDAT(7))
                NTIMESTEP_BC=(NTIMESTEP)+NINT(JDAT(7)/abs(DT))
              ELSE
                NTIMESTEP_BC=NTIMESTEP
              ENDIF
!
!-----------------------------------------------------------------------
!***  Set logical flag to read the BCs
!-----------------------------------------------------------------------
!
              READBC=( (NTIMESTEP==0 .AND. MOD(NTIMESTEP_BC,NBOCO)==0)     &  !<-- Filter related?
!
                                          .OR.                             &
!
                       NTIMESTEP_BC==1                                     &  !<-- First timestep
!
                                          .OR.                             &
!
                     ((MOD(NTIMESTEP_BC,NBOCO)==0) .AND. FILTER_METHOD==0) )  !<-- Non-filter, NBOCO coincident time
!
!-----------------------------------------------------------------------
!
              bc_read: IF(READBC)THEN
!
                bc_flag: IF(WRITE_BC_FLAG==0)THEN
!
                  CALL READ_BC(LM,LNSH,LNSV,NTIMESTEP_BC,DT             &
                              ,RUNBC,IDATBC,IHRSTBC,TBOCO               &
!
                              ,int_state%NVARS_BC_2D_H                  &
                              ,int_state%NVARS_BC_3D_H                  &
                              ,int_state%NVARS_BC_4D_H                  &
                              ,int_state%NVARS_BC_2D_V                  &
                              ,int_state%NVARS_BC_3D_V                  &
!
                              ,int_state%BND_VARS_H                     &
                              ,int_state%BND_VARS_V                     &
                              ,int_state%N_BC_3D_H                      &
                                )
!
                ELSE
!
                  IF (NTIMESTEP==0) THEN
                    CALL WRITE_BC(LM,LNSH,LNSV,NTIMESTEP,DT             &
                            ,RUNBC,TBOCO                                &
                            ,int_state%NVARS_BC_2D_H                    &
                            ,int_state%NVARS_BC_3D_H                    &
                            ,int_state%NVARS_BC_4D_H                    &
                            ,int_state%NVARS_BC_2D_V                    &
                            ,int_state%NVARS_BC_3D_V                    &
                            ,int_state%BND_VARS_H                       &
                            ,int_state%BND_VARS_V                       &
                            ,.TRUE.)                                       !<-- Are tendencies recomputed?
                 ENDIF
!
                ENDIF  bc_flag
!
              ENDIF  bc_read
!
            ENDIF  nest_or_parent
!
!-----------------------------------------------------------------------
!
          ENDIF boundary_tendencies
!
!-----------------------------------------------------------------------
!
          IF(.NOT.int_state%BDY_WAS_READ) THEN
            int_state%BDY_WAS_READ=.TRUE.
          ENDIF
!
!-----------------------------------------------------------------------
!
          btim=timef()
!
          CALL BOCOH                                                    &
            (LM,LNSH,DT,PT                                              &
             ,int_state%PD,int_state%DSG2,int_state%PDSG1               &
             ,int_state%NVARS_BC_2D_H                                   &
             ,int_state%NVARS_BC_3D_H                                   &
             ,int_state%NVARS_BC_4D_H                                   &
             ,int_state%LBND_4D                                         &
             ,int_state%UBND_4D                                         &
             ,int_state%BND_VARS_H                                      &
             ,int_state%PINT)
!
          td%bocoh_tim=td%bocoh_tim+(timef()-btim)
!
        ENDIF bc_update
!
!-----------------------------------------------------------------------
!***  The pressure gradient routine.
!-----------------------------------------------------------------------
!
        btim=timef()
!
        CALL PGFORCE                                                    &
          (int_state%FIRST_STEP,int_state%GLOBAL,int_state%RESTART      &
          ,LM,DT,NTIMESTEP                                              &
          ,RDYV,int_state%DSG2,int_state%PDSG1,int_state%RDXV           &
          ,int_state%WPDAR,int_state%FIS                                &
          ,int_state%PD                                                 &
          ,int_state%T,int_state%Q,int_state%CW                         & ! And how about other TRACER elements?
          ,int_state%PINT                                               &
          ,int_state%RTOP                                               &
          ,int_state%DIV                                                &
          ,int_state%PCNE,int_state%PCNW                                &
          ,int_state%PCX,int_state%PCY                                  &
          ,int_state%TCU,int_state%TCV )
!
        td%pgforce_tim=td%pgforce_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Filtering and boundary conditions for the global forecast.
!-----------------------------------------------------------------------
!
        IF(GLOBAL)THEN
!
          btim=timef()
          CALL FFTFUVN                                                  &
            (LM                                                         &
            ,int_state%KVFILT,int_state%VFILT                           &
            ,int_state%TCU,int_state%TCV                                &
            ,int_state%WFFTRW,int_state%NFFTRW                          &
            ,NUM_PES,MYPE,MPI_COMM_COMP)
          td%fftfwn_tim=td%fftfwn_tim+(timef()-btim)
!
        ENDIF
!
!-----------------------------------------------------------------------
!***  Update the wind field.
!-----------------------------------------------------------------------
!
        btim=timef()
        CALL UPDATEUV                                                   &
         (LM                                                            &
         ,int_state%U,int_state%V                                       &
         ,int_state%TCU,int_state%TCV )
!
        td%updateuv_tim=td%updateuv_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Filtering and boundary conditions for the global forecast.
!-----------------------------------------------------------------------
!
        IF(GLOBAL)THEN
!
          btim=timef()
          CALL SWAPWN(int_state%U,IMS,IME,JMS,JME,LM,INPES)
          CALL SWAPWN(int_state%V,IMS,IME,JMS,JME,LM,INPES)
          td%swapwn_tim=td%swapwn_tim+(timef()-btim)
!
          btim=timef()
          CALL POLEWN(int_state%U,int_state%V                           &
                     ,IMS,IME,JMS,JME,LM,INPES,JNPES)
          td%polewn_tim=td%polewn_tim+(timef()-btim)
!
        ENDIF
!
!-----------------------------------------------------------------------
!
        btim=timef()
        CALL HALO_EXCH(int_state%DIV,LM                                 &
                      ,2,2)
        CALL HALO_EXCH(int_state%U,LM                                   &
                      ,int_state%V,LM                                   &
                      ,2,2)
        td%exch_dyn=td%exch_dyn+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Update the boundary velocity points for the regional forecast.
!-----------------------------------------------------------------------
!
        IF(.NOT.GLOBAL)THEN
!
          btim=timef()
          CALL BOCOV                                                    &
            (LM,LNSV,DT                                                 &
            ,int_state%NVARS_BC_2D_V                                    &
            ,int_state%NVARS_BC_3D_V                                    &
            ,int_state%BND_VARS_V )
          td%bocov_tim=td%bocov_tim+(timef()-btim)
!
        ENDIF
!
!-----------------------------------------------------------------------
!***  The boundary winds have just been updated.  In order to replicate
!***  the integration of a restarted run compared to its free-forecast
!***  counterpart then we must save the wind data in the boundary
!***  arrays for the restart files at this place in the runstream.
!-----------------------------------------------------------------------
!
        IF(QUILTING)THEN
          IF(MOD(NTIMESTEP+1,int_state%NSTEPS_BC_RESTART)==0)THEN          !<-- Look ahead to the end of this timestep
            CALL SAVE_BC_DATA                                           &
              (LM,LNSH,LNSV                                             &
              ,int_state%NVARS_BC_2D_H                                  &
              ,int_state%NVARS_BC_3D_H                                  &
              ,int_state%NVARS_BC_4D_H                                  &
              ,int_state%NVARS_BC_2D_V                                  &
              ,int_state%NVARS_BC_3D_V                                  &
              ,int_state%BND_VARS_H                                     &
              ,int_state%BND_VARS_V                                     &
              ,int_state%NUM_WORDS_BC_SOUTH,int_state%RST_BC_DATA_SOUTH &
              ,int_state%NUM_WORDS_BC_NORTH,int_state%RST_BC_DATA_NORTH &
              ,int_state%NUM_WORDS_BC_WEST ,int_state%RST_BC_DATA_WEST  &
              ,int_state%NUM_WORDS_BC_EAST ,int_state%RST_BC_DATA_EAST  &
              ,EXP_STATE                                                &
              ,int_state%ITS,int_state%ITE,int_state%JTS,int_state%JTE  &
              ,int_state%IMS,int_state%IME,int_state%JMS,int_state%JME  &
              ,int_state%IDS,int_state%IDE,int_state%JDS,int_state%JDE  &
                )
!
          ENDIF
        ENDIF
!
!-----------------------------------------------------------------------
!***  Divergence and horizontal pressure advection in thermo eqn
!-----------------------------------------------------------------------
!
        btim=timef()
!
        CALL DHT                                                        &
          (GLOBAL,LM,DYV,int_state%DSG2,int_state%PDSG1,int_state%DXV   &
          ,int_state%FCP,int_state%FDIV                                 &
          ,int_state%PD,int_state%PDO                                   &
          ,int_state%U,int_state%V                                      &
          ,int_state%OMGALF                                             &
          ,int_state%PCNE,int_state%PCNW,int_state%PCX,int_state%PCY    &
          ,int_state%PFNE,int_state%PFNW,int_state%PFX,int_state%PFY    &
          ,int_state%DIV,int_state%TDIV)
!
        td%dht_tim=td%dht_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Filtering and boundary conditions for the global forecast.
!-----------------------------------------------------------------------
!
        IF(GLOBAL)THEN
!
          btim=timef()
          CALL FFTFHN                                                   &
           (LM                                                          &
           ,int_state%KHFILT                                            &
           ,int_state%HFILT                                             &
           ,int_state%DIV                                               &
           ,int_state%WFFTRH,int_state%NFFTRH                           &
           ,NUM_PES,MYPE,MPI_COMM_COMP)
          td%fftfhn_tim=td%fftfhn_tim+(timef()-btim)
!
          btim=timef()
          CALL SWAPHN                                                   &
           (int_state%DIV                                               &
           ,IMS,IME,JMS,JME,LM                                          &
           ,INPES)
!
          CALL SWAPHN                                                   &
           (int_state%OMGALF                                            &
           ,IMS,IME,JMS,JME,LM                                          &
           ,INPES)
          td%swaphn_tim=td%swaphn_tim+(timef()-btim)
!
          btim=timef()
          CALL POLEHN                                                   &
           (int_state%DIV                                               &
           ,IMS,IME,JMS,JME,LM                                          &
           ,INPES,JNPES)
!
          CALL POLEHN                                                   &
           (int_state%OMGALF                                            &
           ,IMS,IME,JMS,JME,LM                                          &
           ,INPES,JNPES)
          td%polehn_tim=td%polehn_tim+(timef()-btim)
!
        ENDIF
!
!-----------------------------------------------------------------------
!
        btim=timef()
        CALL HALO_EXCH(int_state%DIV,LM                                 &
                      ,int_state%OMGALF,LM                              &
                      ,2,2)
        td%exch_dyn=td%exch_dyn+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Divergence damping
!-----------------------------------------------------------------------
!
        btim=timef()
!
        IF(HDIFF_ON>0)THEN
          CALL DDAMP                                                    &
            (LM                                                         &
            ,DDMPV,PDTOP                                                &
            ,int_state%DSG2,int_state%PDSG1                             &
            ,int_state%SG1,int_state%SG2                                &
            ,int_state%DDMPU                                            &
            ,int_state%FREERUN                                          &
            ,int_state%PD,int_state%PDO                                 &
            ,int_state%U,int_state%V                                    &
            ,int_state%DIV)
        ENDIF
!
        td%ddamp_tim=td%ddamp_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Filtering and boundary conditions for the global forecast.
!-----------------------------------------------------------------------
!
        IF(GLOBAL)THEN
!
          btim=timef()
          CALL SWAPWN                                                   &
            (int_state%U                                                &
            ,IMS,IME,JMS,JME,LM                                         &
            ,INPES)
!
          CALL SWAPWN                                                   &
            (int_state%V                                                &
            ,IMS,IME,JMS,JME,LM                                         &
            ,INPES)
          td%swapwn_tim=td%swapwn_tim+(timef()-btim)
!
          btim=timef()
          CALL POLEWN                                                   &
            (int_state%U,int_state%V                                    &
            ,IMS,IME,JMS,JME,LM                                         &
            ,INPES,JNPES)
          td%polewn_tim=td%polewn_tim+(timef()-btim)
!
        ENDIF
!
!-----------------------------------------------------------------------
!
        btim=timef()
        CALL HALO_EXCH(int_state%U,int_state%LM                         &
                      ,int_state%V,int_state%LM                         &
                      ,2,2)
        td%exch_dyn=td%exch_dyn+(timef()-btim)
!
!-----------------------------------------------------------------------
!
      ENDIF not_firststep
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  The remainder of the Solver integration call sequence
!***  is the same for all timesteps.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      int_state%FIRST_STEP=.FALSE.
!
!-----------------------------------------------------------------------
!***  Update the surface pressure.
!-----------------------------------------------------------------------
!
      btim=timef()
!
      CALL PDTSDT                                                       &
        (LM,DT,int_state%SG2                                            &
        ,int_state%PD                                                   &
        ,int_state%PDO,int_state%PSDT                                   &
        ,int_state%PSGDT                                                &
!
!***  Temporary argument
!
       ,int_state%DIV,int_state%TDIV)
!
      td%pdtsdt_tim=td%pdtsdt_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Filtering and boundary conditions
!-----------------------------------------------------------------------
!
      IF(GLOBAL)THEN
        btim=timef()
        CALL SWAPHN(int_state%PD,IMS,IME,JMS,JME,1,INPES)
        CALL SWAPHN(int_state%PSDT,IMS,IME,JMS,JME,1,INPES)
        td%swaphn_tim=td%swaphn_tim+(timef()-btim)
!
        btim=timef()
        CALL POLEHN(int_state%PD,IMS,IME,JMS,JME,1,INPES,JNPES)
        CALL POLEHN(int_state%PSDT,IMS,IME,JMS,JME,1,INPES,JNPES)
        td%polehn_tim=td%polehn_tim+(timef()-btim)
!
        btim=timef()
        CALL SWAPHN(int_state%PSGDT,IMS,IME,JMS,JME,LM-1,INPES)
        td%swaphn_tim=td%swaphn_tim+(timef()-btim)
!
        btim=timef()
        CALL POLEHN(int_state%PSGDT,IMS,IME,JMS,JME,LM-1,INPES,JNPES)
        td%polehn_tim=td%polehn_tim+(timef()-btim)
      ENDIF
!
!-----------------------------------------------------------------------
!
      btim=timef()
      CALL HALO_EXCH(int_state%PD,1                                     &
                    ,int_state%PSDT,1                                   &
                    ,int_state%PSGDT,LM-1                               &
                    ,2,2)
      td%exch_dyn=td%exch_dyn+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Advection of T, U, and V
!-----------------------------------------------------------------------
!
      btim=timef()
!
      CALL ADV1                                                         &
        (GLOBAL,SECADV                                                  &
        ,LM,LNSAD,INPES,JNPES                                           &
        ,DT,DYV,RDYH,RDYV                                               &
        ,int_state%DSG2,int_state%PDSG1                                 &
        ,int_state%CURV,int_state%DXV,int_state%FAD,int_state%FAH       &
        ,int_state%RDXH,int_state%RDXV,int_state%F                      &
        ,int_state%PD,int_state%PDO                                     &
        ,int_state%OMGALF,int_state%PSGDT                               &
        ,int_state%T,int_state%U,int_state%V                            &
        ,int_state%TP,int_state%UP,int_state%VP                         &
!
!***  Temporary arguments
!
        ,int_state%PFNE,int_state%PFNW                                  &
        ,int_state%PFX,int_state%PFY                                    &
        ,int_state%TCT,int_state%TCU,int_state%TCV)
!
      td%adv1_tim=td%adv1_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!*** Advect specific humidity IDTADTQ time steps
!-----------------------------------------------------------------------
!
q_tracer: IF(MOD(ABS(NTIMESTEP),IDTADTQ)==0)THEN
        KSS=int_state%INDX_Q
        KSE1=KSS
  !
        btim=timef()
        CALL ADV2                                                         &
            (GLOBAL                                                       &
            ,IDTADTQ,KSS,KSE1,LM,LNSAD                                    &
            ,DT,RDYH                                                      &
            ,int_state%DSG2,int_state%PDSG1                               &
            ,int_state%EPSQ2                                              &
            ,int_state%FAH,int_state%RDXH                                 &
            ,int_state%PD,int_state%PDO                                   &
            ,int_state%PSGDT                                              &
            ,int_state%UP,int_state%VP                                    &
            ,int_state%INDX_Q2                                            &
            ,int_state%TRACERS                                            &
            ,int_state%TRACERS_PREV                                       &
!
!***  Temporary arguments
!
            ,int_state%PFNE,int_state%PFNW                                &
            ,int_state%PFX,int_state%PFY                                  &
            ,int_state%TRACERS_SQRT                                       &
            ,int_state%TRACERS_TEND)
        td%adv2_tim=td%adv2_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Filtering and boundary conditions for global forecasts of specific humidity
!-----------------------------------------------------------------------
!
        IF(GLOBAL)THEN
          btim=timef()
          DO KS=KSS,KSE1
             CALL FFTFHN                                                  &
                    (LM                                                   &
                    ,int_state%KHFILT                                     &
                    ,int_state%HFILT                                      &
                    ,int_state%TRACERS_TEND(IMS:IME,JMS:JME,1:LM,KS)      &
                    ,int_state%WFFTRH,int_state%NFFTRH                    &
                    ,NUM_PES,MYPE,MPI_COMM_COMP)
          ENDDO
          td%fftfhn_tim=td%fftfhn_tim+(timef()-btim)
        ENDIF
!
!-----------------------------------------------------------------------
!***  Tracer monotonization for specific humidity
!-----------------------------------------------------------------------
!
        btim=timef()
        CALL MONO                                                         &
            (IDTADTQ,KSS,KSE1,LM                                          &
            ,int_state%DSG2,int_state%PDSG1                               &
            ,int_state%EPSQ2                                              &
            ,int_state%DARE                                               &
            ,int_state%PD                                                 &
            ,int_state%INDX_Q2                                            &
            ,int_state%TRACERS                                            &
            ,INPES,JNPES                                                  &
            ,int_state%USE_ALLREDUCE                                      &
!
!***  Temporary arguments
!
            ,int_state%TRACERS_SQRT                                       &
            ,int_state%TRACERS_TEND)
        td%mono_tim=td%mono_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Update specific humidity
!-----------------------------------------------------------------------
!
          btim=timef()
          CALL UPDATES                                                     &
            (LM,int_state%NUM_TRACERS_TOTAL,KSS,KSE1                       &
            ,int_state%TRACERS,int_state%TRACERS_TEND)
          td%updates_tim=td%updates_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!
          IF(GLOBAL)THEN
            btim=timef()
            CALL SWAPHN(int_state%Q,IMS,IME,JMS,JME,LM,INPES)
            td%swaphn_tim=td%swaphn_tim+(timef()-btim)
            btim=timef()
            CALL POLEHN(int_state%Q,IMS,IME,JMS,JME,LM,INPES,JNPES)
            td%polehn_tim=td%polehn_tim+(timef()-btim)
          ENDIF
!
          btim=timef()
          CALL HALO_EXCH(int_state%Q,LM,2,2)
          td%exch_dyn=td%exch_dyn+(timef()-btim)
      ENDIF  q_tracer
!
!-----------------------------------------------------------------------
!***  Advection of tracers *other* than specific humidity
!-----------------------------------------------------------------------
! 
not_q_tracers: IF(MOD(ABS(NTIMESTEP),IDTADT)==0)THEN
!
!-----------------------------------------------------------------------
!
!-- Water vapor *mixing ratio* is not considered below, it will be 
!   calculated from specific humidity before entering the physics block
!
        IF(int_state%SPEC_ADV)THEN
!
!-- Separate species advection (SPEC_ADV=T): advect Q2 (TKE) and 
!   individual *condensate* species (QC,QI,QR,QS,QG,etc)
!
!-- At the initial time step (NTIMESTEP=0), subroutine UPDATE_WATER has not
!   been called yet, so the initial individual hydrometeor species (QC,QI,etc.) 
!   have not been calculated from the int_state%CW/F_ice/F_rain arrays when
!   initialized from NPS-generated input files.  In the trunk code, the total
!   condensate is advected along with the individual species, so this leads to
!   an initial discrepancy in the values for int_state%CW, which are passed in
!   as input to the radiation code. 
!
!          IF(NTIMESTEP<=0) THEN
!            KSS=int_state%INDX_CW
!          ELSE
!            KSS=int_state%INDX_Q2
!          ENDIF
          KSS=int_state%INDX_CW
          KSE1=int_state%NUM_TRACERS_TOTAL
          IF (RIME_FACTOR_ADVECT) THEN
             btim=timef()
!----------- QG(:,:,:)=F_RIMEF(:,:,:)*QS(:,:,:) for advection
             RIME_FACTOR_INPUT=.TRUE.
             CALL RIME_FACTOR_UPDATE (RIME_FACTOR_INPUT                 &
                                     ,int_state%QS,int_state%QG         &
                                     ,int_state%F_RIMEF                 &
                                     ,IDS,IDE,JDS,JDE,LM                &
                                     ,IMS,IME,JMS,JME                   &
                                     ,ITS,ITE,JTS,JTE)
             td%rfupdate_tim=td%rfupdate_tim+(timef()-btim)
          ENDIF
        ELSE
!
!-- Total condensate advection (SPEC_ADV=F): advect Q2 (TKE) & total condensate (CW)
!
          KSS=int_state%INDX_CW
          KSE1=int_state%INDX_Q2
        ENDIF

!       write(6,*) 'DEBUG-GT: 2nd call to ADV2, kss,kse=',kss,kse1
!
        btim=timef()
!
        CALL ADV2                                                       &
          (GLOBAL                                                       &
          ,IDTADT,KSS,KSE1,LM,LNSAD                                     &
          ,DT,RDYH                                                      &
          ,int_state%DSG2,int_state%PDSG1                               &
          ,int_state%EPSQ2                                              &
          ,int_state%FAH,int_state%RDXH                                 &
          ,int_state%PD,int_state%PDO                                   &
          ,int_state%PSGDT                                              &
          ,int_state%UP,int_state%VP                                    &
          ,int_state%INDX_Q2                                            &
          ,int_state%TRACERS                                            &
          ,int_state%TRACERS_PREV                                       &
!
!***  Temporary arguments
!
          ,int_state%PFNE,int_state%PFNW                                &
          ,int_state%PFX,int_state%PFY                                  &
          ,int_state%TRACERS_SQRT                                       &
          ,int_state%TRACERS_TEND)
!
        td%adv2_tim=td%adv2_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Filtering and boundary conditions for global forecasts
!-----------------------------------------------------------------------
!
        IF(GLOBAL)THEN
!
          btim=timef()
!
          DO KS=KSS,KSE1
                CALL FFTFHN                                             &
                  (LM                                                   &
                  ,int_state%KHFILT                                     &
                  ,int_state%HFILT                                      &
                  ,int_state%TRACERS_TEND(IMS:IME,JMS:JME,1:LM,KS)      &
                  ,int_state%WFFTRH,int_state%NFFTRH                    &
                  ,NUM_PES,MYPE,MPI_COMM_COMP)
          ENDDO
! 
          td%fftfhn_tim=td%fftfhn_tim+(timef()-btim)
!
        ENDIF
!
!-----------------------------------------------------------------------
!***  Tracer monotonization
!-----------------------------------------------------------------------
!
        btim=timef()
!
        CALL MONO                                                       &
          (IDTADT,KSS,KSE1,LM                                           &
          ,int_state%DSG2,int_state%PDSG1                               &
          ,int_state%EPSQ2                                              &
          ,int_state%DARE                                               &
          ,int_state%PD                                                 &
          ,int_state%INDX_Q2                                            &
          ,int_state%TRACERS                                            &
          ,INPES,JNPES                                                  &
          ,int_state%USE_ALLREDUCE                                      &
!
!***  Temporary arguments
!
          ,int_state%TRACERS_SQRT                                       &
          ,int_state%TRACERS_TEND)
!
        td%mono_tim=td%mono_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Update tracers
!-----------------------------------------------------------------------
!
        btim=timef()
!
        CALL UPDATES                                                     &
          (LM,int_state%NUM_TRACERS_TOTAL,KSS,KSE1                       &
          ,int_state%TRACERS,int_state%TRACERS_TEND)
!
        td%updates_tim=td%updates_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!
if_global: IF(GLOBAL)THEN    !-- Global NMMB
!
          btim=timef()
          CALL SWAPHN(int_state%CW,IMS,IME,JMS,JME,LM,INPES)
          CALL SWAPHN(int_state%O3,IMS,IME,JMS,JME,LM,INPES)
          CALL SWAPHN(int_state%Q2,IMS,IME,JMS,JME,LM,INPES)
!..Need a similar set of lines for the TRACERS array at some point.
!
          td%swaphn_tim=td%swaphn_tim+(timef()-btim)
!
          btim=timef()
          CALL POLEHN(int_state%CW,IMS,IME,JMS,JME,LM,INPES,JNPES)
          CALL POLEHN(int_state%O3,IMS,IME,JMS,JME,LM,INPES,JNPES)
          CALL POLEHN(int_state%Q2,IMS,IME,JMS,JME,LM,INPES,JNPES)
!..Need a similar set of lines for the TRACERS array at some point.
!
          td%polehn_tim=td%polehn_tim+(timef()-btim)
!
          btim=timef()
          CALL HALO_EXCH(int_state%CW,LM                                &
                      ,int_state%O3,LM                                  &
                      ,int_state%Q2,LM                                  &
                      ,2,2)
          td%exch_dyn=td%exch_dyn+(timef()-btim)
!
        ELSE  if_global      !-- Regional NMMB
!
          btim=timef()
          DO KS=KSS,KSE1
            CALL HALO_EXCH(                                             &
                int_state%TRACERS(IMS:IME,JMS:JME,1:LM,KS),LM           &
               ,2,2)
          ENDDO
          td%exch_dyn=td%exch_dyn+(timef()-btim)
!
          IF (RIME_FACTOR_ADVECT) THEN
            btim=timef()
!--------------- F_RIMEF(:,:,:)=QG(:,:,:)/QS(:,:,:) for physics (after advection)
            RIME_FACTOR_INPUT=.FALSE.
            CALL RIME_FACTOR_UPDATE (RIME_FACTOR_INPUT                  &
                                     ,int_state%QS,int_state%QG         &
                                     ,int_state%F_RIMEF                 &
                                     ,IDS,IDE,JDS,JDE,LM                &
                                     ,IMS,IME,JMS,JME                   &
                                     ,ITS,ITE,JTS,JTE)
             td%rfupdate_tim=td%rfupdate_tim+(timef()-btim)
           ENDIF
!
        ENDIF  if_global
!
!-----------------------------------------------------------------------
!
      ENDIF not_q_tracers
!
!-----------------------------------------------------------------------
!***  Interface pressures and horizontal part of Omega-Alpha term
!-----------------------------------------------------------------------
!
      btim=timef()
!
      CALL VTOA                                                         &
        (LM,DT,EF4T,PT,int_state%SG2                                    &
        ,int_state%PSDT                                                 &
        ,int_state%DWDT,int_state%RTOP                                  &
        ,int_state%OMGALF                                               &
        ,int_state%PINT                                                 &
!
!***  Temporary arguments
!
        ,int_state%TDIV,int_state%TCT)
!
      td%vtoa_tim=td%vtoa_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Filtering and boundary conditions for global forecasts
!-----------------------------------------------------------------------
!
      IF(GLOBAL)THEN
!
        btim=timef()
        CALL FFTFHN                                                     &
          (LM                                                           &
          ,int_state%KHFILT                                             &
          ,int_state%HFILT                                              &
          ,int_state%TCT                                                &
          ,int_state%WFFTRH,int_state%NFFTRH                            &
          ,NUM_PES,MYPE,MPI_COMM_COMP)
        td%fftfhn_tim=td%fftfhn_tim+(timef()-btim)
!
      ENDIF
!
!-----------------------------------------------------------------------
!***  Update the temperature field.
!-----------------------------------------------------------------------
!
      btim=timef()
!
      CALL UPDATET                                                      &
        (LM                                                             &
        ,int_state%T                                                    &
!
!***  Temporary argument
!
        ,int_state%TCT)
!
      td%updatet_tim=td%updatet_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Filtering and boundary conditions for global forecasts
!-----------------------------------------------------------------------
!
      IF(GLOBAL)THEN
!
        btim=timef()
        CALL SWAPHN(int_state%OMGALF,IMS,IME,JMS,JME,LM,INPES)
        CALL SWAPHN(int_state%PINT,IMS,IME,JMS,JME,LM+1,INPES)
        CALL SWAPHN(int_state%T,IMS,IME,JMS,JME,LM,INPES)
        td%swaphn_tim=td%swaphn_tim+(timef()-btim)
!
        btim=timef()
        CALL POLEHN(int_state%OMGALF,IMS,IME,JMS,JME,LM,INPES,JNPES)
        CALL POLEHN(int_state%PINT,IMS,IME,JMS,JME,LM+1,INPES,JNPES)
        CALL POLEHN(int_state%T,IMS,IME,JMS,JME,LM,INPES,JNPES)
        td%polehn_tim=td%polehn_tim+(timef()-btim)
!
      ENDIF
!
!-----------------------------------------------------------------------
!
      btim=timef()
      CALL HALO_EXCH(int_state%OMGALF,LM                                &
                    ,int_state%PINT,LM+1                                &
                    ,2,2)
      CALL HALO_EXCH(int_state%T,LM                                     &
                    ,2,2)
      td%exch_dyn=td%exch_dyn+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Nonhydrostatic advection of height
!-----------------------------------------------------------------------
!
      btim=timef()
!
      CALL CDZDT                                                        &
        (GLOBAL,HYDRO                                                   &
        ,LM,DT,int_state%DSG2,int_state%PDSG1                           &
        ,int_state%FAH,int_state%FIS                                    &
        ,int_state%PD,int_state%PDO                                     &
        ,int_state%PSGDT                                                &
        ,int_state%CW,int_state%Q,int_state%RTOP,int_state%T            & 
        ,int_state%PINT                                                 &
        ,int_state%DWDT,int_state%PDWDT,int_state%W,int_state%BARO      &
        ,int_state%Z                                                    &
!
!***  temporary arguments
!
        ,int_state%PFNE,int_state%PFNW,int_state%PFX,int_state%PFY)
!
      td%cdzdt_tim=td%cdzdt_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Filtering and boundary conditions for global forecasts
!-----------------------------------------------------------------------
!
      IF(GLOBAL)THEN
!
        btim=timef()
        CALL FFTFHN                                                     &
          (LM                                                           &
          ,int_state%KHFILT                                             &
          ,int_state%HFILT                                              &
          ,int_state%W                                                  &
          ,int_state%WFFTRH,int_state%NFFTRH                            &
          ,NUM_PES,MYPE,MPI_COMM_COMP)
        td%fftfhn_tim=td%fftfhn_tim+(timef()-btim)
!
        btim=timef()
        CALL SWAPHN(int_state%W,IMS,IME,JMS,JME,LM,INPES)
        td%swaphn_tim=td%swaphn_tim+(timef()-btim)
!
        btim=timef()
        CALL POLEHN(int_state%W,IMS,IME,JMS,JME,LM,INPES,JNPES)
        td%polehn_tim=td%polehn_tim+(timef()-btim)
!
      ENDIF
!
!-----------------------------------------------------------------------
!
      btim=timef()
      CALL HALO_EXCH(int_state%W,LM                                     &
                    ,3,3)
      td%exch_dyn=td%exch_dyn+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Advection of W (with internal halo exchange)
!-----------------------------------------------------------------------
!
      btim=timef()
!
      CALL CDWDT                                                        &
        (GLOBAL,HYDRO,int_state%RESTART                                 &
        ,INPES,JNPES,LM,ABS(NTIMESTEP)                                  &
        ,DT,G,int_state%DSG2,int_state%PDSG1,int_state%PSGML1           &
        ,int_state%FAH                                                  &
        ,int_state%HDACX,int_state%HDACY                                &
        ,int_state%PD,int_state%PDO                                     &
        ,int_state%PSGDT                                                &
        ,int_state%DWDT,int_state%PDWDT,int_state%W                     &
        ,int_state%PINT                                                 &
!
!***  External scratch areas
!
        ,int_state%DEF,int_state%PFX,int_state%PFY                      &
        ,int_state%PFNE,int_state%PFNW)
!
      td%cdwdt_tim=td%cdwdt_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Filtering and boundary conditions for global forecasts
!-----------------------------------------------------------------------
!
      IF(GLOBAL)THEN
!
        btim=timef()
        CALL FFTFHN                                                     &
          (LM                                                           &
          ,int_state%KHFILT                                             &
          ,int_state%HFILT                                              &
          ,int_state%DWDT                                               &
          ,int_state%WFFTRH,int_state%NFFTRH                            &
          ,NUM_PES,MYPE,MPI_COMM_COMP)
        td%fftfhn_tim=td%fftfhn_tim+(timef()-btim)
!
        btim=timef()
        CALL SWAPHN(int_state%DWDT,IMS,IME,JMS,JME,LM,INPES)
        td%swaphn_tim=td%swaphn_tim+(timef()-btim)
!
        btim=timef()
        CALL POLEHN(int_state%DWDT,IMS,IME,JMS,JME,LM,INPES,JNPES)
        td%polehn_tim=td%polehn_tim+(timef()-btim)
!
      ENDIF
!
!-----------------------------------------------------------------------
!
      btim=timef()
      CALL HALO_EXCH(int_state%DWDT,LM                                  &
                    ,2,2)
      td%exch_dyn=td%exch_dyn+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Vertically propagating fast waves
!-----------------------------------------------------------------------
!
      btim=timef()
!
      CALL VSOUND                                                       &
        (GLOBAL,HYDRO,int_state%RESTART                                 &
        ,LM,ABS(NTIMESTEP)                                              &
        ,CP,DT,PT,int_state%DSG2,int_state%PDSG1                        &
        ,int_state%PD                                                   &
        ,int_state%CW,int_state%Q,int_state%RTOP                        &
        ,int_state%DWDT,int_state%T,int_state%W,int_state%W_TOT         &
        ,int_state%BARO                                                 &
        ,int_state%PINT)
!
      td%vsound_tim=td%vsound_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Filtering and boundary conditions for global forecasts
!-----------------------------------------------------------------------
!
      IF(GLOBAL)THEN
!
        btim=timef()
        CALL POAVHN                                                     &
          (IMS,IME,JMS,JME,LM                                           &
          ,int_state%DWDT                                               &
          ,INPES,JNPES                                                  &
          ,int_state%USE_ALLREDUCE)
        CALL POAVHN                                                     &
          (IMS,IME,JMS,JME,LM                                           &
          ,int_state%W                                                  &
          ,INPES,JNPES                                                  &
          ,int_state%USE_ALLREDUCE)
        CALL POAVHN                                                     &
          (IMS,IME,JMS,JME,LM                                           &
          ,int_state%PINT                                               &
          ,INPES,JNPES                                                  &
          ,int_state%USE_ALLREDUCE)
        td%poavhn_tim=td%poavhn_tim+(timef()-btim)
!
        btim=timef()
        CALL SWAPHN(int_state%DWDT,IMS,IME,JMS,JME,LM,INPES)
        CALL SWAPHN(int_state%T,IMS,IME,JMS,JME,LM,INPES)
        CALL SWAPHN(int_state%W,IMS,IME,JMS,JME,LM,INPES)
        CALL SWAPHN(int_state%PINT,IMS,IME,JMS,JME,LM+1,INPES)
        td%swaphn_tim=td%swaphn_tim+(timef()-btim)
!
        btim=timef()
        CALL POLEHN(int_state%DWDT,IMS,IME,JMS,JME,LM,INPES,JNPES)
        CALL POLEHN(int_state%T,IMS,IME,JMS,JME,LM,INPES,JNPES)
        CALL POLEHN(int_state%W,IMS,IME,JMS,JME,LM,INPES,JNPES)
        CALL POLEHN(int_state%PINT,IMS,IME,JMS,JME,LM+1,INPES,JNPES)
        td%polehn_tim=td%polehn_tim+(timef()-btim)
!
      ENDIF
!
!-----------------------------------------------------------------------
!
      btim=timef()
      CALL HALO_EXCH(int_state%DWDT,LM                                  &
                    ,int_state%T,LM                                     &
                    ,2,2)
      CALL HALO_EXCH(int_state%W,LM                                     &
                    ,int_state%PINT,LM+1                                &
                    ,2,2)
      td%exch_dyn=td%exch_dyn+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Save DT to compare and see if sign has changed for filtering.
!-----------------------------------------------------------------------
!
      int_state%DT_LAST=DT_TEST
      int_state%DT_TEST_RATIO=REAL(INTEGER_DT)+REAL(NUMERATOR_DT)       &
                                              /REAL(IDENOMINATOR_DT)
      int_state%FILTER_METHOD_LAST=FILTER_METHOD
!
!-----------------------------------------------------------------------
!***  NOTE:  The Solver export state is fully updated now
!***         because subroutine SOLVER_INITIALIZE inserted the 
!***         appropriate ESMF Fields into it.  Those Fields 
!***         contain pointers to the actual data and those
!***         pointers are never re-directed, i.e., no explicit
!***         action is needed to update the Solver export state.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Write the layer statistics for temperature.
!-----------------------------------------------------------------------
!
      IF(MOD(ABS(NTIMESTEP)+1,N_PRINT_STATS)==0)THEN
!
        IF(int_state%PRINT_DIAG .OR. int_state%PRINT_ALL) &
        CALL FIELD_STATS(INT_STATE%T,MYPE,MPI_COMM_COMP,LM              &
                        ,ITS,ITE,JTS,JTE                                &
                        ,IMS,IME,JMS,JME                                &
                        ,IDS,IDE,JDS,JDE)
      ENDIF
!
      td%solver_dyn_tim=td%solver_dyn_tim+(timef()-btim0)
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!----PHY_RUN START -----------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!rv - please do not remove this template call:
!     call exit('dyn',int_state%pint,int_state%t,int_state%q            &
!                    ,int_state%u,int_state%v,int_state%q2,int_state%w  &
!                    ,ntimestep,mype,my_domain_id,mpi_comm_comp         &
!                    ,ids,ide,jds,jde,lm                                &
!                    ,ims,ime,jms,jme                                   &
!                    ,its,ite,jts,jte)
!     if(mod(nint(dt*ntimestep),3600)==0)then
!       call twr(int_state%t,lm,'tphy',ntimestep,mype,num_pes,mpi_comm_comp &
!               ,ids,ide,jds,jde &
!               ,ims,ime,jms,jme &
!               ,its,ite,jts,jte &
!               ,my_domain_id )
!       call vwr(int_state%u,lm,'uphy',ntimestep,mype,num_pes,mpi_comm_comp &
!               ,ids,ide,jds,jde &
!               ,ims,ime,jms,jme &
!               ,its,ite,jts,jte &
!               ,my_domain_id )
!     endif
!rv
!
      physics: IF(INTEGER_DT>0)THEN                                     !<-- Physics is active
!
      btim0=timef()
!
!-----------------------------------------------------------------------
!***  Call radiation so that updated fields are written to the
!***  history files after 0 hours.
!-----------------------------------------------------------------------
!
      IF(NTIMESTEP==0)THEN
         NTIMESTEP_RAD=NTIMESTEP
      ELSE
         NTIMESTEP_RAD=NTIMESTEP+1
      ENDIF
!
!-----------------------------------------------------------------------
!***  Dereference some internal state components for convenience.
!-----------------------------------------------------------------------
!
      NPRECIP=int_state%NPRECIP
      PDTOP=int_state%PDTOP
      PT=int_state%PT
!
!-----------------------------------------------------------------------
!
      gfs_phys_test: IF(.NOT.int_state%GFS)THEN                            !<-- NMM-B physics is NOT the GFS package
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  At the appropriate times, reset the various min/max/average
!***  diagnostic fields to begin accumulating for the next period
!-----------------------------------------------------------------------
!
      IF(NTIMESTEP == 0 .or. MOD(NTIMESTEP,NSTEPS_PER_RESET)==0) THEN
        DO J=JTS,JTE
        DO I=ITS,ITE
          int_state%TLMAX(I,J)=-999.
          int_state%TLMIN(I,J)=999.
          int_state%T02MAX(I,J)=-999.
          int_state%T02MIN(I,J)=999.
          int_state%RH02MAX(I,J)=-999.
          int_state%RH02MIN(I,J)=999.
          int_state%SPD10MAX(I,J)=-999.
          int_state%UPHLMAX(I,J)=0.
          int_state%U10MAX(I,J)=-999.
          int_state%V10MAX(I,J)=-999.
          int_state%UPVVELMAX(I,J)=-999.
          int_state%DNVVELMAX(I,J)=999.
          int_state%T10AVG(I,J)=0.
          int_state%T10(I,J)=0.
          int_state%PSFCAVG(I,J)=0.
          int_state%AKHSAVG(I,J)=0.
          int_state%AKMSAVG(I,J)=0.
          int_state%SNOAVG(I,J)=0.
          int_state%REFDMAX(I,J)=DBZmin
          int_state%PRATEMAX(I,J)=0
          int_state%FPRATEMAX(I,J)=0
          int_state%UPHLMAX(I,J)=-999.
        ENDDO
        ENDDO
!
        int_state%NCOUNT=0
      ENDIF
!
!     IF (mod(int_state%NTSD,NSTEPS_PER_CHECK) == 0) THEN
      IF (mod(int_state%NTSD,NSTEPS_PER_CHECK) == 0 .and. FILTER_METHOD==0 ) THEN
!
max_hrly: IF (TRIM(int_state%MICROPHYSICS) == 'fer') THEN
!
          CALL MAX_FIELDS(int_state%T,int_state%Q,int_state%U            &
                         ,int_state%V,int_state%CW                       &
                         ,int_state%F_RAIN,int_state%F_ICE               &
                         ,int_state%F_RIMEF,int_state%Z                  &
                         ,int_state%W_TOT,int_state%PINT                 &
                         ,int_state%PD,int_state%PREC                    &
                         ,int_state%CPRATE,int_state%HTOP                &
                         ,int_state%T2,int_state%U10,int_state%V10       &
                         ,int_state%PSHLTR,int_state%TSHLTR              &
                         ,int_state%QSHLTR                               &
                         ,int_state%SGML2,int_state%PSGML1               &
                         ,int_state%REFDMAX,int_state%PRATEMAX           &
                         ,int_state%FPRATEMAX,int_state%SR               &
                         ,int_state%UPVVELMAX,int_state%DNVVELMAX        &
                         ,int_state%TLMAX,int_state%TLMIN                &
                         ,int_state%T02MAX,int_state%T02MIN              &
                         ,int_state%RH02MAX,int_state%RH02MIN            &
                         ,int_state%U10MAX,int_state%V10MAX              &
                         ,int_state%TH10,int_state%T10                   &
                         ,int_state%SPD10MAX,int_state%T10AVG            &
                         ,int_state%PSFCAVG                              &
                         ,int_state%AKHS,int_state%AKMS                  &
                         ,int_state%AKHSAVG,int_state%AKMSAVG            &
                         ,int_state%SNO,int_state%SNOAVG                 &
                         ,int_state%UPHLMAX                              &
                         ,int_state%DT,int_state%NPHS,int_state%NTSD     &
                         ,int_state%DXH,int_state%DYH                    &
                         ,int_state%FIS                                  &
                         ,ITS,ITE,JTS,JTE                                &
                         ,IMS,IME,JMS,JME                                &
                         ,IDE,JDE                                        &
                         ,ITS_B1,ITE_B1,JTS_B1,JTE_B1                    &
                         ,LM,int_state%NCOUNT,int_state%FIRST_NMM        &
                         ,MY_DOMAIN_ID                                   &
                                        )
!
        ELSEIF (TRIM(int_state%MICROPHYSICS) == 'fer_hires') THEN  max_hrly
!
          CALL MAX_FIELDS_HR(int_state%T,int_state%Q,int_state%U         &
                            ,int_state%V,int_state%CW                    &
                            ,int_state%F_RAIN,int_state%F_ICE            &
                            ,int_state%F_RIMEF,int_state%Z               &
                            ,int_state%W_TOT,int_state%refl_10cm         &
                            ,int_state%PINT,int_state%PD,int_state%PREC  &
                            ,int_state%CPRATE,int_state%HTOP             &
                            ,int_state%T2,int_state%U10,int_state%V10    &
                            ,int_state%PSHLTR,int_state%TSHLTR           &
                            ,int_state%QSHLTR                            &
                            ,int_state%SGML2,int_state%PSGML1            &
                            ,int_state%REFDMAX,int_state%PRATEMAX        &
                            ,int_state%FPRATEMAX,int_state%SR            &
                            ,int_state%UPVVELMAX,int_state%DNVVELMAX     &
                            ,int_state%TLMAX,int_state%TLMIN             &
                            ,int_state%T02MAX,int_state%T02MIN           &
                            ,int_state%RH02MAX,int_state%RH02MIN         &
                            ,int_state%U10MAX,int_state%V10MAX           &
                            ,int_state%TH10,int_state%T10                &
                            ,int_state%SPD10MAX,int_state%T10AVG         &
                            ,int_state%PSFCAVG                           &
                            ,int_state%AKHS,int_state%AKMS               &
                            ,int_state%AKHSAVG,int_state%AKMSAVG         &
                            ,int_state%SNO,int_state%SNOAVG              &
                            ,int_state%UPHLMAX                           &
                            ,int_state%DT,int_state%NPHS,int_state%NTSD  &
                            ,int_state%DXH,int_state%DYH                 &
                            ,int_state%FIS                               &
                            ,ITS,ITE,JTS,JTE                             &
                            ,IMS,IME,JMS,JME                             &
                            ,IDE,JDE                                     &
                            ,ITS_B1,ITE_B1,JTS_B1,JTE_B1                 &
                            ,LM,int_state%NCOUNT,int_state%FIRST_NMM     &
                            ,MY_DOMAIN_ID                                &
                                           )
!
       ELSEIF (TRIM(int_state%MICROPHYSICS) == 'wsm6') THEN  max_hrly
!
         CALL MAX_FIELDS_W6(int_state%T,int_state%Q,int_state%U         &
                           ,int_state%V,int_state%Z,int_state%W_TOT     &
                           ,int_state%QR,int_state%QS,int_state%QG      &
                           ,int_state%PINT,int_state%PD,int_state%PREC  &
                           ,int_state%CPRATE,int_state%HTOP             &
                           ,int_state%T2,int_state%U10,int_state%V10    &
                           ,int_state%PSHLTR,int_state%TSHLTR           &
                           ,int_state%QSHLTR                            &
                           ,int_state%SGML2,int_state%PSGML1            &
                           ,int_state%REFDMAX,int_state%PRATEMAX        &
                           ,int_state%FPRATEMAX,int_state%SR            &
                           ,int_state%UPVVELMAX,int_state%DNVVELMAX     &
                           ,int_state%TLMAX,int_state%TLMIN             &
                           ,int_state%T02MAX,int_state%T02MIN           &
                           ,int_state%RH02MAX,int_state%RH02MIN         &
                           ,int_state%U10MAX,int_state%V10MAX           &
                           ,int_state%TH10,int_state%T10                &
                           ,int_state%SPD10MAX,int_state%T10AVG         &
                           ,int_state%PSFCAVG                           &
                           ,int_state%AKHS,int_state%AKMS               &
                           ,int_state%AKHSAVG,int_state%AKMSAVG         &
                           ,int_state%SNO,int_state%SNOAVG              &
                           ,int_state%UPHLMAX                           &
                           ,int_state%DT,int_state%NPHS,int_state%NTSD  &
                           ,int_state%DXH,int_state%DYH                 &
                           ,int_state%FIS                               &
                           ,ITS,ITE,JTS,JTE                             &
                           ,IMS,IME,JMS,JME                             &
                           ,IDE,JDE                                     &
                           ,ITS_B1,ITE_B1,JTS_B1,JTE_B1                 &
                           ,LM                                          &
                           ,int_state%NCOUNT,int_state%FIRST_NMM        &
                           ,MY_DOMAIN_ID                                &
                                           )
!
       ELSEIF (TRIM(int_state%MICROPHYSICS) == 'thompson') THEN  max_hrly
!
         CALL MAX_FIELDS_THO(int_state%T,int_state%Q,int_state%U        &
                           ,int_state%V,int_state%Z,int_state%W_TOT     &
                           ,int_state%refl_10cm                         &
                           ,int_state%PINT,int_state%PD,int_state%PREC  &
                           ,int_state%CPRATE,int_state%HTOP             &
                           ,int_state%T2,int_state%U10,int_state%V10    &
                           ,int_state%PSHLTR,int_state%TSHLTR           &
                           ,int_state%QSHLTR                            &
                           ,int_state%SGML2,int_state%PSGML1            &
                           ,int_state%REFDMAX,int_state%PRATEMAX        &
                           ,int_state%FPRATEMAX,int_state%SR            &
                           ,int_state%UPVVELMAX,int_state%DNVVELMAX     &
                           ,int_state%TLMAX,int_state%TLMIN             &
                           ,int_state%T02MAX,int_state%T02MIN           &
                           ,int_state%RH02MAX,int_state%RH02MIN         &
                           ,int_state%U10MAX,int_state%V10MAX           &
                           ,int_state%TH10,int_state%T10                &
                           ,int_state%SPD10MAX,int_state%T10AVG         &
                           ,int_state%PSFCAVG                           &
                           ,int_state%AKHS,int_state%AKMS               &
                           ,int_state%AKHSAVG,int_state%AKMSAVG         &
                           ,int_state%SNO,int_state%SNOAVG              &
                           ,int_state%UPHLMAX                           &
                           ,int_state%DT,int_state%NPHS,int_state%NTSD  &
                           ,int_state%DXH,int_state%DYH                 &
                           ,int_state%FIS                               &
                           ,ITS,ITE,JTS,JTE                             &
                           ,IMS,IME,JMS,JME                             &
                           ,IDE,JDE                                     &
                           ,ITS_B1,ITE_B1,JTS_B1,JTE_B1                 &
                           ,LM                                          &
                           ,int_state%NCOUNT,int_state%FIRST_NMM        &
                           ,MY_DOMAIN_ID                                &
                                           )
!
        ENDIF  max_hrly
!
      ENDIF
!
!-----------------------------------------------------------------------
!***  Set logical switches for calling each of the Physics schemes.
!-----------------------------------------------------------------------
!
        CALL_SHORTWAVE=MOD(NTIMESTEP_RAD,int_state%NRADS)==0
        CALL_LONGWAVE=MOD(NTIMESTEP_RAD,int_state%NRADL)==0
        CALL_TURBULENCE=MOD(NTIMESTEP,int_state%NPHS)==0
        CALL_PRECIP=MOD(NTIMESTEP,NPRECIP)==0
!
!-----------------------------------------------------------------------
!***  Update WATER array from CWM, F_ICE, F_RAIN for Ferrier 
!***  microphysics but only if any of the Physics subroutines 
!***  are called (subroutine UPDATE_WATER is after subroutine
!***  PHYSICS_INITIALIZE in this module).
!
!***  Expanded to also update CWM, F_ICE, F_RAIN, F_RIMEF for non-Ferrier
!***  microphysics.
!-----------------------------------------------------------------------
!
        update_wtr: IF((int_state%MICROPHYSICS=='fer'                   &
                                   .OR.                                 &
                        int_state%MICROPHYSICS=='fer_hires'             &
                                   .OR.                                 &
                        int_state%MICROPHYSICS=='gfs'                   &
                                   .OR.                                 &
                        int_state%MICROPHYSICS=='wsm6'                  &
                                   .OR.                                 &
                        int_state%MICROPHYSICS=='thompson')             &
                                   .AND.                                &
                       (CALL_SHORTWAVE .OR. CALL_LONGWAVE .OR.          &
                        CALL_TURBULENCE .OR. CALL_PRECIP) ) THEN
!
!          write(*,*) 'DEBUG-GT, now calling UPDATE_WATER'
           CALL UPDATE_WATER(int_state%CW                               &
                            ,int_state%F_ICE                            &
                            ,int_state%F_RAIN                           &
                            ,int_state%F_RIMEF                          &
                            ,int_state%T                                &
                            ,int_state%QC                               &
                            ,int_state%QR                               &
                            ,int_state%QS                               &
                            ,int_state%QI                               &
                            ,int_state%QG                               &
                            ,int_state%MICROPHYSICS                     &
                            ,int_state%SPEC_ADV                         &
                            ,NTIMESTEP                                  &
                            ,IDS,IDE,JDS,JDE,LM                         &
                            ,IMS,IME,JMS,JME                            &
                            ,ITS,ITE,JTS,JTE)
!
        ENDIF update_wtr
!
!---------------------------------------------------------------------
!***  Precipitation Adjustment
!-----------------------------------------------------------------------
!
!***
!***      Call READPCP to
!***            1) READ IN PRECIPITATION FOR HOURS 1, 2 and 3;
!***            2) Initialize DDATA to 999. (this is the amount
!***               of input precip allocated to each physics time step
!***               in ADJPPT; TURBL/SURFCE, which uses DDATA, is called
!***               before ADJPPT)
!***            3) Initialize LSPA to zero
!***
!-----------------------------------------------------------------------
!
        IF(int_state%NTSD==0)THEN
          IF(int_state%PCPFLG .and. FILTER_METHOD == 0)THEN
            CALL READPCP(MYPE,MPI_COMM_COMP                             &
                        ,int_state%PPTDAT                               &
                        ,int_state%DDATA                                &
                        ,int_state%LSPA                                 &
                        ,int_state%PCPHR                                &
                        ,MY_DOMAIN_ID                                   &
                        ,IDS,IDE,JDS,JDE,LM                             &
                        ,IMS,IME,JMS,JME                                &
                        ,ITS,ITE,JTS,JTE                                &
                        ,ITS_B1,ITE_B1,JTS_B2,JTE_B2)
          ENDIF
        ENDIF
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  Call the individual physical processes.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  Radiation
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!***  Radiation needs some specific time quantities.  Use NTIMESTEP_rad 
!***  for the next time step ahead of the current time so that the
!***  radiation fields can be updated prior to being written to
!***  output (BSF 10/6/2010).
!
        CALL TIME_MEASURE(START_YEAR,START_MONTH,START_DAY,START_HOUR   &
                         ,START_MINUTE,START_SECOND                     &
                         ,NTIMESTEP_rad,int_state%DT                    &
                         ,JULDAY,JULYR,JULIAN,XTIME)
!
!-----------------------------------------------------------------------
        radiatn: IF(CALL_SHORTWAVE.OR.CALL_LONGWAVE)THEN
!-----------------------------------------------------------------------
!
          btim=timef()
!         write(*,*) 'DEBUG-GT, now calling RADIATION ', btim
!
!-----------------------------------------------------------------------
!***  Temporary switch between radiation schemes placed in SOLVER_RUN
!***  rather than inside RADIATION_DRIVER (will be done later)
!-----------------------------------------------------------------------
!
          CALL ESMF_ClockGet(clock       =CLOCK_ATM                     &  !<-- The ESMF Clock
                            ,startTime   =STARTTIME                     &  !<-- The start time (ESMF) on the clock
                            ,currTime    =CURRTIME                      &  !<-- The current time (ESMF) on the clock
                            ,rc          =RC)
!
          CALL ESMF_TimeGet(time=STARTTIME                              &  !<-- The start forecast time (ESMF)
                           ,yy  =IDAT(1)                                &  !<-- The start forecast year (integer)
                           ,mm  =IDAT(2)                                &  !<-- The start forecast month (integer)
                           ,dd  =IDAT(3)                                &  !<-- The start forecast day (integer)
                           ,h   =IDAT(5)                                &  !<-- The start forecast hour (integer)
                           ,m   =IDAT(6)                                &  !<-- The start forecast minute (integer)
                           ,s   =IDAT(7)                                &  !<-- The start forecast second (integer)
                           ,rc  =RC)
          IDAT(4)=0
          IDAT(8)=0
!
          CALL ESMF_TimeGet(time=CURRTIME                               &  !<-- The cuurent forecast time (ESMF)
                           ,yy  =JDAT(1)                                &  !<-- The current forecast year (integer)
                           ,mm  =JDAT(2)                                &  !<-- The current forecast month (integer)
                           ,dd  =JDAT(3)                                &  !<-- The current forecast day (integer)
                           ,h   =JDAT(5)                                &  !<-- The current forecast hour (integer)
                           ,m   =JDAT(6)                                &  !<-- The current forecast minute (integer)
                           ,s   =JDAT(7)                                &  !<-- The current forecast second (integer)
                           ,rc  =RC)
          JDAT(4)=0
          JDAT(8)=0
!
          CALL RADIATION(NTIMESTEP_RAD                                  &
                        ,int_state%DT,JULDAY,JULYR,XTIME,JULIAN         &
                        ,START_HOUR,int_state%NPHS                      &
                        ,int_state%GLAT,int_state%GLON                  &
                        ,int_state%NRADS,int_state%NRADL                &
                        ,int_state%DSG2,int_state%SGML2,int_state%SG2   &
                        ,int_state%PDSG1,int_state%PSGML1               &
                        ,int_state%PSG1                                 &
                        ,int_state%PT,int_state%PD                      &
                        ,int_state%T,int_state%Q                        &
                        ,int_state%THS,int_state%ALBEDO                 &
                        ,int_state%QC,int_state%QR                      &
                        ,int_state%QI,int_state%QS,int_state%QG         &
                        ,int_state%NI                                   &
                        ,int_state%F_QC,int_state%F_QR                  &
                        ,int_state%F_QI,int_state%F_QS,int_state%F_QG   &
                        ,int_state%F_NI                                 &
                        ,int_state%NUM_WATER                            &
                        ,int_state%SM,int_state%CLDFRA                  &
                        ,int_state%RLWTT,int_state%RSWTT                &
                        ,int_state%RLWIN,int_state%RSWIN                &
                        ,int_state%RSWINC,int_state%RSWOUT              &
                        ,int_state%RLWTOA,int_state%RSWTOA              &
                        ,int_state%CZMEAN,int_state%SIGT4               &
                        ,int_state%CFRACL,int_state%CFRACM              &
                        ,int_state%CFRACH                               &
                        ,int_state%ACFRST,int_state%NCFRST              &
                        ,int_state%ACFRCV,int_state%NCFRCV              &
                        ,int_state%CUPPT,int_state%SNO                  &
                        ,int_state%HTOP,int_state%HBOT                  &
                        ,int_state%SHORTWAVE,int_state%LONGWAVE         &
                        ,int_state%CLDFRACTION                          &
                        ,int_state%DYH                                  &
!---- RRTM part ---------------------------------------------------------
                        ,int_state%DT_INT,JDAT                          &
                        ,int_state%CW,int_state%O3                      &
                        ,int_state%F_ICE,int_state%F_RAIN               &
                        ,int_state%F_RIMEF                              &
                        ,int_state%SI,int_state%TSKIN                   &
                        ,int_state%Z0,int_state%SICE                    &
                        ,int_state%MXSNAL,int_state%SGM                 &
                        ,int_state%STDH,int_state%OMGALF                &
                        ,int_state%SNOWC                                &
!------------------------------------------------------------------------
                        ,LM)
!
          td%radiation_tim=td%radiation_tim+(timef()-btim)
!
        ENDIF radiatn
!
!-----------------------------------------------------------------------
!***  Empty the ACFRST and ACFRCV accumulation arrays if it is time
!***  to do so prior to their being updated by the radiation.
!-----------------------------------------------------------------------
!
        IF(MOD(NTIMESTEP,int_state%NCLOD)==0)THEN
          DO J=JTS,JTE
          DO I=ITS,ITE
            int_state%ACFRST(I,J)=0.
            int_state%ACFRCV(I,J)=0.
            int_state%NCFRST(I,J)=0
            int_state%NCFRCV(I,J)=0
          ENDDO
          ENDDO
        ENDIF
!
!-----------------------------------------------------------------------
!***  Update the temperature with the radiative tendency.
!-----------------------------------------------------------------------
!
        btim=timef()
!
        CALL RDTEMP(NTIMESTEP,int_state%DT,JULDAY,JULYR,START_HOUR      &
                   ,int_state%GLAT,int_state%GLON                       &
                   ,int_state%CZEN,int_state%CZMEAN,int_state%T         &
                   ,int_state%RSWTT,int_state%RLWTT                     &
                   ,IDS,IDE,JDS,JDE,LM                                  &
                   ,IMS,IME,JMS,JME                                     &
                   ,ITS,ITE,JTS,JTE                                     &
                   ,ITS_B1,ITE_B1,JTS_B1,JTE_B1)
!
        td%rdtemp_tim=td%rdtemp_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Poles and East-West boundary.
!-----------------------------------------------------------------------
!
        IF(int_state%GLOBAL)THEN
          btim=timef()
!
          CALL SWAPHN(int_state%RSWIN,IMS,IME,JMS,JME,1,int_state%INPES)
          CALL POLEHN(int_state%RSWIN,IMS,IME,JMS,JME,1                  &
                     ,int_state%INPES,int_state%JNPES)
!
          CALL SWAPHN(int_state%T,IMS,IME,JMS,JME,LM,int_state%INPES)
          CALL POLEHN(int_state%T,IMS,IME,JMS,JME,LM                     &
                     ,int_state%INPES,int_state%JNPES)
!
          td%pole_swap_tim=td%pole_swap_tim+(timef()-btim)
        ENDIF
!
!-----------------------------------------------------------------------
!***  Empty the accumulators of sfc energy flux and sfc hydrology if
!***  it is time to do so prior to their being updated by turbulence.
!-----------------------------------------------------------------------
!
        IF(MOD(NTIMESTEP,int_state%NRDLW)==0)THEN
          DO J=JTS,JTE
          DO I=ITS,ITE
            int_state%ALWIN(I,J) =0.
            int_state%ALWOUT(I,J)=0.
            int_state%ALWTOA(I,J)=0.
            int_state%ARDLW(I,J) =0.                                       !<-- An artificial 2-D array
                                                                           !    (ESMF cannot have an evolving scalar Attribute)
          ENDDO
          ENDDO
        ENDIF
!
        IF(MOD(NTIMESTEP,int_state%NRDSW)==0)THEN
          DO J=JTS,JTE
          DO I=ITS,ITE
            int_state%ASWIN(I,J)=0.
            int_state%ASWOUT(I,J)=0.
            int_state%ASWTOA(I,J)=0.
            int_state%ARDSW(I,J) =0.                                       !<-- An artificial 2-D array 
                                                                           !    (ESMF cannot have an evolving scalar Attribute)
          ENDDO
          ENDDO
        ENDIF
!
        IF(MOD(NTIMESTEP,int_state%NSRFC)==0)THEN
          DO J=JTS,JTE
          DO I=ITS,ITE
            int_state%SFCSHX(I,J)=0.
            int_state%SFCLHX(I,J)=0.
            int_state%SUBSHX(I,J)=0.
            int_state%SNOPCX(I,J)=0.
            int_state%POTFLX(I,J)=0.
            int_state%ASRFC(I,J) =0.                                       !<-- An artificial 2-D array
                                                                           !    (ESMF cannot have an evolving scalar Attribute)
          ENDDO
          ENDDO
        ENDIF
!
        IF(MOD(NTIMESTEP,int_state%NPREC)==0)THEN
          DO J=JTS,JTE
          DO I=ITS,ITE
            int_state%ACSNOW(I,J)=0.
              int_state%ACSNOM(I,J)=0.
            int_state%SSROFF(I,J)=0.
            int_state%BGROFF(I,J)=0.
            int_state%SFCEVP(I,J)=0.
            int_state%POTEVP(I,J)=0.
          ENDDO
          ENDDO
        ENDIF
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  Turbulence, Sfc Layer, and Land Surface
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
        turbulence: IF(CALL_TURBULENCE)THEN
!
          btim=timef()
!         write(*,*) 'DEBUG-GT, now calling TURBL ', btim
!
          DO L=1,NUM_SOIL_LAYERS
            DZSOIL(L)=SLDPTH(L)
          ENDDO
!
          IF(int_state%PCPFLG .and. FILTER_METHOD == 0)THEN
            LOC_PCPFLG=int_state%PCPFLG
          ELSE
            LOC_PCPFLG=.FALSE.
          ENDIF
!
          CALL TURBL(NTIMESTEP,int_state%DT,int_state%NPHS              &
                    ,NUM_SOIL_LAYERS,SLDPTH,DZSOIL                      &
                    ,int_state%DSG2,int_state%SGML2,int_state%SG2       &
                    ,int_state%PDSG1,int_state%PSGML1,int_state%PSG1,PT &
                    ,int_state%EPSL,int_state%EPSQ2                     &
                    ,int_state%SM,int_state%CZEN,int_state%CZMEAN       &
                    ,int_state%SIGT4,int_state%RLWIN,int_state%RSWIN    &
                    ,int_state%RADOT                                    &
                    ,int_state%RLWTT,int_state%RSWTT                    &
                    ,int_state%PD,int_state%T                           &
                    ,int_state%Q,int_state%CW                           &
                    ,int_state%F_ICE,int_state%F_RAIN,int_state%F_RIMEF &
                    ,int_state%SR,int_state%Q2,int_state%U,int_state%V  &
                    ,int_state%DUDT,int_state%DVDT                      &
                    ,int_state%THS,int_state%TSKIN,int_state%SST        &
                    ,int_state%PREC,int_state%SNO                       &
                    ,int_state%SNOWC                                    &
                    ,int_state%QC,int_state%QR                          &
                    ,int_state%QI,int_state%QS,int_state%QG             &
                    ,int_state%F_QC,int_state%F_QR                      &
                    ,int_state%F_QI,int_state%F_QS,int_state%F_QG       &
                    ,int_state%FIS,int_state%Z0,int_state%Z0BASE        &
                    ,int_state%USTAR,int_state%PBLH,int_state%LPBL      &
                    ,int_state%XLEN_MIX,int_state%RMOL                  &
                    ,int_state%EXCH_H,int_state%AKHS,int_state%AKMS     &
                    ,int_state%AKHS_OUT,int_state%AKMS_OUT              &
                    ,int_state%THZ0,int_state%QZ0                       &
                    ,int_state%UZ0,int_state%VZ0                        &
                    ,int_state%QSH,int_state%MAVAIL                     &
                    ,int_state%STC,int_state%SMC,int_state%CMC          &
                    ,int_state%SMSTAV,int_state%SMSTOT                  &
                    ,int_state%SSROFF,int_state%BGROFF                  &
                    ,int_state%IVGTYP,int_state%ISLTYP,int_state%VEGFRC &
                    ,int_state%GRNFLX                                   &
                    ,int_state%SFCEXC,int_state%ACSNOW,int_state%ACSNOM &
                    ,int_state%SNOPCX,int_state%SICE                    &
                    ,int_state%TG,int_state%SOILTB                      &
                    ,int_state%ALBASE,int_state%MXSNAL,int_state%ALBEDO &
                    ,int_state%SH2O,int_state%SI,int_state%EPSR         &
                    ,int_state%U10,int_state%V10                        &
                    ,int_state%TH10,int_state%Q10                       &
                    ,int_state%TSHLTR,int_state%QSHLTR,int_state%PSHLTR &
                    ,int_state%PSFC,int_state%T2                        &
                    ,int_state%TWBS,int_state%QWBS                      &
                    ,int_state%SFCSHX,int_state%SFCLHX,int_state%SFCEVP &
                    ,int_state%TAUX,int_state%TAUY                      &
                    ,int_state%POTEVP,int_state%POTFLX,int_state%SUBSHX &
                    ,int_state%APHTIM                                   &
                    ,int_state%ARDSW,int_state%ARDLW                    &
                    ,int_state%ASRFC                                    &
                    ,int_state%CROT,int_state%SROT,int_state%MIXHT      &
                    ,int_state%HSTDV,int_state%HCNVX,int_state%HASYW    &
                    ,int_state%HASYS,int_state%HASYSW,int_state%HASYNW  &
                    ,int_state%HLENW,int_state%HLENS,int_state%HLENSW   &
                    ,int_state%HLENNW,int_state%HANGL,int_state%HANIS   &
                    ,int_state%HSLOP,int_state%HZMAX                    &
                    ,int_state%CDMB,int_state%CLEFF,int_state%SIGFAC    &
                    ,int_state%FACTOP,int_state%RLOLEV                  &
                    ,int_state%DPMIN                                    &
                    ,int_state%RSWOUT,int_state%RSWTOA,int_state%RLWTOA &
                    ,int_state%ASWIN,int_state%ASWOUT,int_state%ASWTOA  &
                    ,int_state%ALWIN,int_state%ALWOUT,int_state%ALWTOA  &
                    ,int_state%GWDFLG,LOC_PCPFLG                        &
                    ,int_state%DDATA,int_state%UCMCALL,int_state%IVEGSRC&
                    ,int_state%TURBULENCE,int_state%SFC_LAYER           &
                    ,int_state%LAND_SURFACE                             &
                    ,int_state%MICROPHYSICS                             &
                    ,int_state%LISS_RESTART                             &
                    ,int_state%GLOBAL                                   &
 !!! HURRICANE PBL/SFCLAY
                    ,int_state%VAR_RIC,int_state%COEF_RIC_L             &
                    ,int_state%COEF_RIC_S,int_state%DISHEAT             &
                    ,int_state%ALPHA,int_state%SFENTH                   &
!!! HURRICANE

                    ,IDS,IDE,JDS,JDE,LM                                 &
                    ,IMS,IME,JMS,JME                                    &
                    ,ITS,ITE,JTS,JTE)
!
          td%turbl_tim=td%turbl_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Exchange wind tendencies.
!-----------------------------------------------------------------------
!
          btim=timef()
!
          CALL HALO_EXCH(int_state%DUDT,LM,int_state%DVDT,LM,1,1)
!
          td%exch_phy=td%exch_phy+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Now interpolate wind tendencies from H to V points.
!-----------------------------------------------------------------------
!
          btim=timef()
!
          CALL H_TO_V_TEND(int_state%DUDT,int_state%DT,int_state%NPHS,LM &
                          ,int_state%U)
          CALL H_TO_V_TEND(int_state%DVDT,int_state%DT,int_state%NPHS,LM &
                          ,int_state%V)
!
          td%h_to_v_tim=td%h_to_v_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Poles and East-West boundary.
!-----------------------------------------------------------------------
!
          IF(int_state%GLOBAL)THEN
            btim=timef()
!
            CALL SWAPHN(int_state%T,IMS,IME,JMS,JME,LM,int_state%INPES)
            CALL POLEHN(int_state%T,IMS,IME,JMS,JME,LM                  &
                       ,int_state%INPES,int_state%JNPES)
!
            CALL SWAPHN(int_state%Q,IMS,IME,JMS,JME,LM,int_state%INPES)
            CALL POLEHN(int_state%Q,IMS,IME,JMS,JME,LM                  &
                       ,int_state%INPES,int_state%JNPES)
!
            CALL SWAPHN(int_state%CW,IMS,IME,JMS,JME,LM,int_state%INPES)
            CALL POLEHN(int_state%CW,IMS,IME,JMS,JME,LM                 &
                       ,int_state%INPES,int_state%JNPES)
!
            CALL SWAPHN(int_state%Q2,IMS,IME,JMS,JME,LM,int_state%INPES)
            CALL POLEHN(int_state%Q2,IMS,IME,JMS,JME,LM                 &
                       ,int_state%INPES,int_state%JNPES)
!
            CALL SWAPWN(int_state%U,IMS,IME,JMS,JME,LM,int_state%INPES)
            CALL SWAPWN(int_state%V,IMS,IME,JMS,JME,LM,int_state%INPES)
            CALL POLEWN(int_state%U,int_state%V,IMS,IME,JMS,JME,LM      &
                       ,int_state%INPES,int_state%JNPES)
!
            td%pole_swap_tim=td%pole_swap_tim+(timef()-btim)
          ENDIF
!
!-----------------------------------------------------------------------
!***  Exchange wind components and TKE.
!-----------------------------------------------------------------------
!
          btim=timef()
!
          CALL HALO_EXCH(int_state%U,LM,int_state%V,LM                  &
                        ,2,2)
!
          CALL HALO_EXCH(int_state%UZ0,1,int_state%VZ0,1                &
                        ,int_state%Q2,LM                                &
                        ,1,1)
!
!-----------------------------------------------------------------------
!***  Exchange other variables that are needed for parents' 
!***  interpolations to interior points of moving nests.
!-----------------------------------------------------------------------
!
          CALL HALO_EXCH(int_state%ALBEDO,1                             &
                        ,int_state%EPSR,1                               &
                        ,int_state%QSH,1                                &
                        ,int_state%QWBS,1,1,1)
          CALL HALO_EXCH(int_state%QZ0,1                                &
                        ,int_state%SOILTB,1                             &
                        ,int_state%THS,1                                &
                        ,int_state%THZ0,1,1,1)
          CALL HALO_EXCH(int_state%USTAR,1                              &
                        ,int_state%UZ0,1                                &
                        ,int_state%VZ0,1                                &
                        ,int_state%Z0,1,1,1)
          CALL HALO_EXCH(int_state%TSKIN,1                              &
                        ,int_state%CMC,1,1,1)
          CALL HALO_EXCH(int_state%SMC,NUM_SOIL_LAYERS                  &
                        ,int_state%SH2O,NUM_SOIL_LAYERS                 &
                        ,int_state%STC,NUM_SOIL_LAYERS,1,1)
!
          td%exch_phy=td%exch_phy+(timef()-btim)
!
!-----------------------------------------------------------------------
!
        ENDIF turbulence
!
!----------------------------------------------------------------------- 
!***  Empty the accumulators of precipitation and latent heating if is
!***  is time prior to their being updated by convection/microphysics.
!-----------------------------------------------------------------------
!
        IF(MOD(NTIMESTEP,int_state%NPREC)==0)THEN
          DO J=JTS,JTE
          DO I=ITS,ITE
            int_state%ACPREC(I,J)=0.
            int_state%CUPREC(I,J)=0.
          ENDDO
          ENDDO
        ENDIF
!
        IF(MOD(NTIMESTEP,int_state%NHEAT)==0)THEN
          DO J=JTS,JTE
          DO I=ITS,ITE
            int_state%AVCNVC(I,J)=0.   !- was a scalar, now 2D for ESMF
            int_state%AVRAIN(I,J)=0.   !- was a scalar, now 2D for ESMF
          ENDDO
          ENDDO
!
          DO L=1,LM
          DO J=JTS,JTE
          DO I=ITS,ITE
            int_state%TRAIN(I,J,L)=0.
            int_state%TCUCN(I,J,L)=0.
            do KK=1,int_state%d_ss
              int_state%MPRATES(I,J,L,KK)=0.
            enddo
          ENDDO
          ENDDO
          ENDDO
        ENDIF    !-- IF(MOD(NTSD_BUCKET,NHEAT)==0)THEN
!
!-----------------------------------------------------------------------
!***  1 of 3 calls to CLTEND, save Told array before convection & microphysics
!-----------------------------------------------------------------------
!
        cld_tend1: IF(CALL_PRECIP .AND. int_state%NPRECIP>1) THEN
            btim=timef()
            ICLTEND=-1
            CALL CLTEND(ICLTEND,int_state%NPRECIP,int_state%T           &
                       ,int_state%Told,int_state%Tadj                   &
                       ,IDS,IDE,JDS,JDE,LM                              &
                       ,IMS,IME,JMS,JME                                 &
                       ,ITS,ITE,JTS,JTE)
            td%cltend_tim=td%cltend_tim+(timef()-btim)
         ENDIF cld_tend1
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  Convection
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
        convection: IF(CALL_PRECIP.AND.int_state%CONVECTION/='none')THEN
!
          btim=timef()
!         write(*,*) 'DEBUG-GT, now calling CUCNVC ', btim
!
!-----------------------------------------------------------------------
          IF(int_state%CONVECTION=='bmj' .OR. &
             int_state%CONVECTION=='sas' .OR. & 
             int_state%CONVECTION=='scalecu' .OR. &
             int_state%CONVECTION=='sashur') THEN
!
            CALL CUCNVC(NTIMESTEP,int_state%DT,int_state%NPRECIP          &
                       ,int_state%NRADS,int_state%NRADL                   &
                       ,int_state%MINUTES_HISTORY                         &
                       ,int_state%ENTRAIN,int_state%NEWALL                &
                       ,int_state%NEWSWAP,int_state%NEWUPUP               &
                       ,int_state%NODEEP                                  &
                       ,int_state%FRES,int_state%FR                       &
                       ,int_state%FSL,int_state%FSS                       &
                       ,int_state%DYH,int_state%RESTART,int_state%HYDRO   &
                       ,int_state%CLDEFI                                  &
                       ,int_state%F_ICE,int_state%F_RAIN                  &
                       ,int_state%QC,int_state%QR                         &
                       ,int_state%QI,int_state%QS,int_state%QG            &
                       ,int_state%F_QC,int_state%F_QR                     &
                       ,int_state%F_QI,int_state%F_QS,int_state%F_QG      &
                       ,int_state%DSG2,int_state%SGML2,int_state%SG2      &
                       ,int_state%PDSG1,int_state%PSGML1,int_state%PSG1   &
                       ,int_state%DXH                                     &
                       ,int_state%PT,int_state%PD                         &
                       ,int_state%T,int_state%Q                           &
                       ,int_state%CW,int_state%TCUCN                      &
                       ,int_state%OMGALF                                  &
                       ,int_state%U,int_state%V                           &
                       ,int_state%FIS,int_state%W0AVG                     &
                       ,int_state%PREC,int_state%ACPREC                   &
                       ,int_state%CUPREC,int_state%ACPREC_TOT             &
                       ,int_state%CUPPT,int_state%CPRATE                  &
                       ,int_state%CNVBOT,int_state%CNVTOP                 &
                       ,int_state%SM,int_state%LPBL                       &
                       ,int_state%HTOP,int_state%HTOPD,int_state%HTOPS    &
                       ,int_state%HBOT,int_state%HBOTD,int_state%HBOTS    &
                       ,int_state%AVCNVC,int_state%ACUTIM                 &
                       ,int_state%RSWIN,int_state%RSWOUT                  &
                       ,int_state%CONVECTION,int_state%CU_PHYSICS         &
                       ,int_state%MICROPHYSICS                            &
                       ,int_state%SICE,int_state%QWBS,int_state%TWBS      &
                       ,int_state%PBLH,int_state%DUDT,int_state%DVDT      &
!!!  added for SAS-hurricane
                       ,int_state%SAS_MOMMIX,int_state%SAS_PGCON          &   !hwrf,namelist
                       ,int_state%SAS_MASS_FLUX                           &   !hwrf,namelist
                       ,int_state%SAS_SHALCONV,int_state%SAS_SHAL_PGCON   &   !hwrf,namelist
                       ,int_state%W_TOT,int_state%PSGDT                    &
!!!  SAS-huricane
                       ,A2,A3,A4,CAPPA,CP,ELIV,ELWV,EPSQ,G                &
                       ,P608,PQ0,R_D,TIW                                  &
                       ,IDS,IDE,JDS,JDE,LM                                &
                       ,IMS,IME,JMS,JME                                   &
                       ,ITS,ITE,JTS,JTE                                   &
                       ,ITS_B1,ITE_B1,JTS_B1,JTE_B1)
!
          ELSE
!
!           write(0,*)' Invalid selection for convection scheme'
          STOP
!
          ENDIF
!
          td%cucnvc_tim=td%cucnvc_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***    Poles and East-West boundary.
!-----------------------------------------------------------------------
!
          IF(int_state%GLOBAL)THEN
            btim=timef()
!
            CALL SWAPHN(int_state%T,IMS,IME,JMS,JME,LM,int_state%INPES)
            CALL POLEHN(int_state%T,IMS,IME,JMS,JME,LM                  &
                       ,int_state%INPES,int_state%JNPES)
!
            CALL SWAPHN(int_state%Q,IMS,IME,JMS,JME,LM,int_state%INPES)
            CALL POLEHN(int_state%Q,IMS,IME,JMS,JME,LM                  &
                       ,int_state%INPES,int_state%JNPES)
!
            td%pole_swap_tim=td%pole_swap_tim+(timef()-btim)
          ENDIF
!
!-----------------------------------------------------------------------
!***  Exchange wind tendencies for SAS and bmj schemes.
!-----------------------------------------------------------------------
!
          wind: IF (int_state%CONVECTION=='sas' .or. &
                    int_state%CONVECTION=='bmj') THEN !zj
!
!-----------------------------------------------------------------------
!
            btim=timef()
            CALL HALO_EXCH(int_state%DUDT,LM,int_state%DVDT,LM,1,1)
            td%exch_phy=td%exch_phy+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Now interpolate wind tendencies from H to V points.
!-----------------------------------------------------------------------
!
            btim=timef()
            CALL H_TO_V_TEND(int_state%DUDT,int_state%DT                &
                            ,int_state%NPRECIP,LM                       &
                            ,int_state%U)
            CALL H_TO_V_TEND(int_state%DVDT,int_state%DT                &
                            ,int_state%NPRECIP,LM                       &
                            ,int_state%V)
            td%h_to_v_tim=td%h_to_v_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Poles and East-West boundary.
!-----------------------------------------------------------------------
!
            IF(int_state%GLOBAL)THEN
              btim=timef()
!
              CALL SWAPWN(int_state%U,IMS,IME,JMS,JME,LM                &
                         ,int_state%INPES)
              CALL SWAPWN(int_state%V,IMS,IME,JMS,JME,LM                &
                         ,int_state%INPES)
              CALL POLEWN(int_state%U,int_state%V,IMS,IME,JMS,JME,LM    &
                         ,int_state%INPES,int_state%JNPES)
!
              td%pole_swap_tim=td%pole_swap_tim+(timef()-btim)
            ENDIF
!
!-----------------------------------------------------------------------
!***  Exchange wind components.
!-----------------------------------------------------------------------
!
            btim=timef()
            CALL HALO_EXCH(int_state%U,LM,int_state%V,LM                &
                          ,2,2)
            td%exch_phy=td%exch_phy+(timef()-btim)
!
!-----------------------------------------------------------------------
!
          ENDIF wind
!
!-----------------------------------------------------------------------
!
        ENDIF convection
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  Microphysics
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
        microphysics: IF(CALL_PRECIP)THEN
!
          btim=timef()
!         write(*,*) 'DEBUG-GT, now calling GSMDRIVE ', btim
!
          CALL GSMDRIVE(NTIMESTEP,int_state%DT                             &
                       ,NPRECIP                                            &
                       ,int_state%DXH(JC),int_state%DYH                    &
                       ,int_state%SM,int_state%FIS                         &
                       ,int_state%DSG2,int_state%SGML2                     &
                       ,int_state%PDSG1,int_state%PSGML1                   &
                       ,int_state%PT,int_state%PD                          &
                       ,int_state%T,int_state%Q                            &
                       ,int_state%CW,int_state%OMGALF                      &
                       ,int_state%TRAIN,int_state%SR                       &
                       ,int_state%F_ICE,int_state%F_RAIN,int_state%F_RIMEF &
                       ,int_state%QC,int_state%QR                          &
                       ,int_state%QI,int_state%QS,int_state%QG             &
                       ,int_state%NI,int_state%NR                          & ! G. Thompson
                       ,int_state%F_QC,int_state%F_QR                      &
                       ,int_state%F_QI,int_state%F_QS,int_state%F_QG       &
                       ,int_state%F_NI,int_state%F_NR                      & ! G. Thompson
                       ,int_state%PREC,int_state%ACPREC                    &
                       ,int_state%AVRAIN,int_state%ACPREC_TOT              &
                       ,int_state%acpcp_ra,int_state%acpcp_sn,int_state%acpcp_gr &  ! G. Thompson
                       ,int_state%refl_10cm                                &  !  G. Thompson
                       ,int_state%re_cloud,int_state%re_ice,int_state%re_snow  &  !  G. Thompson
                       ,int_state%has_reqc,int_state%has_reqi,int_state%has_reqs  &  !  G. Thompson
                       ,int_state%MP_RESTART_STATE                         &
                       ,int_state%TBPVS_STATE,int_state%TBPVS0_STATE       &
                       ,int_state%SPECIFIED,int_state%NESTED               &
                       ,int_state%MICROPHYSICS                             &
                       ,int_state%RHGRD                                    &  ! fer_hires only
                       ,int_state%TP1                                      &  !gfs mod-brad
                       ,int_state%QP1                                      &  !gfs mod-brad
                       ,int_state%PSP1                                     &  !gfs mod-brad
                       ,USE_RADAR                                          &
                       ,int_state%DFI_TTEN                                 &
                       ,IDS,IDE,JDS,JDE,LM                                 &
                       ,IMS,IME,JMS,JME                                    &
                       ,ITS,ITE,JTS,JTE                                    &
                       ,ITS_B1,ITE_B1,JTS_B1,JTE_B1,int_state%MPRATES      &
                       ,int_state%D_SS)
!
          td%gsmdrive_tim=td%gsmdrive_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  2 of 3 calls to CLTEND, calculate Tadj and replace T with Told
!-----------------------------------------------------------------------
!
        cld_tend2: IF(int_state%NPRECIP>1) THEN
          btim=timef()
          ICLTEND=0
          CALL CLTEND(ICLTEND,int_state%NPRECIP,int_state%T             &
                     ,int_state%Told,int_state%Tadj                     &
                     ,IDS,IDE,JDS,JDE,LM                                &
                     ,IMS,IME,JMS,JME                                   &
                     ,ITS,ITE,JTS,JTE)
          td%cltend_tim=td%cltend_tim+(timef()-btim)
        ENDIF  cld_tend2
!
!-----------------------------------------------------------------------
!***  Precipitation Assimilation
!-----------------------------------------------------------------------
!
          IF (int_state%PCPFLG .and. FILTER_METHOD == 0) THEN
!
            btim=timef()
            CALL CHKSNOW(MYPE                                           &
                        ,int_state%NTSD                                 &
                        ,int_state%DT                                   &
                        ,int_state%NPHS                                 &
                        ,int_state%SR                                   &
                        ,int_state%PPTDAT                               &
                        ,int_state%PCPHR                                &
                        ,IDS,IDE,JDS,JDE,LM                             &
                        ,IMS,IME,JMS,JME                                &
                        ,ITS,ITE,JTS,JTE                                &
                        ,ITS_B1,ITE_B1,JTS_B2,JTE_B2)
!
            CALL ADJPPT(MYPE                                            &
                       ,int_state%NTSD                                  &
                       ,int_state%DT                                    &
                       ,int_state%NPHS                                  &
                       ,int_state%PREC                                  &
                       ,int_state%LSPA                                  &
                       ,int_state%PPTDAT                                &
                       ,int_state%DDATA                                 &
                       ,int_state%PCPHR                                 &
                       ,IDS,IDE,JDS,JDE,LM                              &
                       ,IMS,IME,JMS,JME                                 &
                       ,ITS,ITE,JTS,JTE                                 &
                       ,ITS_B1,ITE_B1,JTS_B2,JTE_B2)
!
            td%adjppt_tim=td%adjppt_tim+(timef()-btim)
!
          ENDIF
!
!-----------------------------------------------------------------------
!***  Poles and East-West boundary.
!-----------------------------------------------------------------------
!
          IF(int_state%GLOBAL)THEN
            btim=timef()
!
!bsf: Apply these after last (3rd) call to CLTEND below
!
!            CALL SWAPHN(int_state%T,IMS,IME,JMS,JME,LM,int_state%INPES)
!            CALL POLEHN(int_state%T,IMS,IME,JMS,JME,LM                  &
!                       ,int_state%INPES,int_state%JNPES)
!
            CALL SWAPHN(int_state%Q,IMS,IME,JMS,JME,LM,int_state%INPES)
            CALL POLEHN(int_state%Q,IMS,IME,JMS,JME,LM                  &
                       ,int_state%INPES,int_state%JNPES)
!
            CALL SWAPHN(int_state%CW,IMS,IME,JMS,JME,LM,int_state%INPES)
            CALL POLEHN(int_state%CW,IMS,IME,JMS,JME,LM                 &
                       ,int_state%INPES,int_state%JNPES)
!
            td%pole_swap_tim=td%pole_swap_tim+(timef()-btim)
          ENDIF
!
!-----------------------------------------------------------------------
!
        ENDIF microphysics
!
!-----------------------------------------------------------------------
!***  3 of 3 calls to CLTEND, incremental updates of T using Told & Tadj
!-----------------------------------------------------------------------
!
        cld_tend3: IF(int_state%NPRECIP>1) THEN
          btim=timef()
          ICLTEND=1
          CALL CLTEND(ICLTEND,int_state%NPRECIP,int_state%T             &
                     ,int_state%Told,int_state%Tadj                     &
                     ,IDS,IDE,JDS,JDE,LM                                &
                     ,IMS,IME,JMS,JME                                   &
                     ,ITS,ITE,JTS,JTE)
          td%cltend_tim=td%cltend_tim+(timef()-btim)
        ENDIF  cld_tend3
!
!-----------------------------------------------------------------------
!***  Prevent supersaturation w/r/t water and smooth temperature profiles
!     if lapse rates are steeper than dry adiabatic above lowest levels.
!-----------------------------------------------------------------------
!
        btim=timef()
        CALL TQADJUST(int_state%T,int_state%Q,int_state%QC              &
                     ,int_state%CW,int_state%F_ICE,int_state%F_RAIN     &
                     ,int_state%PD,int_state%DSG2,int_state%PDSG1       &
                     ,int_state%PSGML1,int_state%SGML2                  &
                     ,int_state%SPEC_ADV,int_state%RHGRD                &
                     ,IDS,IDE,JDS,JDE,LM                                &
                     ,IMS,IME,JMS,JME                                   &
                     ,ITS,ITE,JTS,JTE)
        td%tqadjust_tim=td%tqadjust_tim+(timef()-btim)
!
!bsf: Call SWAPHN & POLEHN for temperature here after temperature update
!
        IF(int_state%GLOBAL)THEN
           btim=timef()
!
           CALL SWAPHN(int_state%T,IMS,IME,JMS,JME,LM,int_state%INPES)
           CALL POLEHN(int_state%T,IMS,IME,JMS,JME,LM                  &
                      ,int_state%INPES,int_state%JNPES)
           td%pole_swap_tim=td%pole_swap_tim+(timef()-btim)
        ENDIF
!
!-----------------------------------------------------------------------
!***  Exchange T, Q, CW now every timestep; also QC for species advection
!-----------------------------------------------------------------------
!
        btim=timef()
        CALL HALO_EXCH(int_state%T,LM,2,2)
        CALL HALO_EXCH(int_state%Q,LM,int_state%CW,LM,2,2)
        IF(int_state%SPEC_ADV) THEN
          CALL HALO_EXCH(int_state%QC,LM,2,2)
!
!-----------------------------------------------------------------------
!***    Exchange various cloud species for separate species advection
!-----------------------------------------------------------------------
!
          IF(CALL_PRECIP .OR. CALL_TURBULENCE) THEN
            IF(int_state%F_QR) CALL HALO_EXCH(int_state%QR,LM,2,2)
            IF(int_state%F_QS) CALL HALO_EXCH(int_state%QS,LM,2,2)
            IF(int_state%F_QI) CALL HALO_EXCH(int_state%QI,LM,2,2)
            IF(int_state%F_QG) CALL HALO_EXCH(int_state%QG,LM,2,2)
            IF(int_state%F_NI) CALL HALO_EXCH(int_state%NI,LM,2,2)
            IF(int_state%F_NR) CALL HALO_EXCH(int_state%NR,LM,2,2)
          ENDIF
        ENDIF
        td%exch_phy=td%exch_phy+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  NOTE:  The Physics export state is fully updated now
!***         because subroutine PHY_INITIALIZE inserted the
!***         appropriate ESMF Fields into it.  Those Fields
!***         contain pointers to the actual data and those
!***         pointers are never re-directed.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      ELSE gfs_phys_test                                                   !<-- Use GFS physics package
#if 1
        WRITE(0,*)'Init of GFS phys in NMMB disabled, 20140812, jm'
        CALL ESMF_FINALIZE(endflag=ESMF_END_ABORT                   &
                          ,rc             =RC)
#else
!
!-----------------------------------------------------------------------
!
!#######################################################################
!#######################################################################
!############ G F S   P H Y S I C S   D R I V E R ######################
!#######################################################################
!#######################################################################
!
        btim=timef()
!
        ALLOCATE(LONSPERLAR(JTS:JTE))
        ALLOCATE(NLNSP(JTS:JTE))
        ALLOCATE(GLOBAL_LATS_R(JTS:JTE))
        ALLOCATE(CLDCOV_V(LM))
        ALLOCATE(PRSL(LM))
        ALLOCATE(PRSLK(LM))
        ALLOCATE(GU(LM))
        ALLOCATE(GV(LM))
        ALLOCATE(GT(LM))
        ALLOCATE(GR(LM))
        ALLOCATE(VVEL(LM))
        ALLOCATE(F_ICE(LM))
        ALLOCATE(F_RAIN(LM))
        ALLOCATE(R_RIME(LM))
        ALLOCATE(ADT(LM))
        ALLOCATE(ADU(LM))
        ALLOCATE(ADV(LM))
        ALLOCATE(PHIL(LM))
        ALLOCATE(GR3(LM,NTRAC))
        ALLOCATE(ADR(LM,NTRAC))
        ALLOCATE(PRSI(LM+1))
        ALLOCATE(PRSIK(LM+1))
        ALLOCATE(RSGM(LM+1))
        ALLOCATE(PHII(LM+1))
        ALLOCATE(SINLAT_R(JTS:JTE))
        ALLOCATE(COSLAT_R(JTS:JTE))
        ALLOCATE(XLON(ITS:ITE,JTS:JTE))
        ALLOCATE(COSZEN(ITS:ITE,JTS:JTE))
        ALLOCATE(COSZDG(ITS:ITE,JTS:JTE))
        ALLOCATE(SINLAT_V(ITS:ITE,JTS:JTE))
        ALLOCATE(COSLAT_V(ITS:ITE,JTS:JTE))
        ALLOCATE(RANN(ITS:ITE,JTS:JTE))
        ALLOCATE(RANNUM((ITE-ITS+1)*(JTE-JTS+1)))
        ALLOCATE(GR1(1,LM,NTRAC-1))
        ALLOCATE(SWH(LM))
        ALLOCATE(HLW(LM))
        ALLOCATE(DKH(LM))
        ALLOCATE(RNP(LM))
        ALLOCATE(UPD_MF(LM))
        ALLOCATE(DWN_MF(LM))
        ALLOCATE(DET_MF(LM))
        ALLOCATE(DQDT(LM))
        ALLOCATE(DQ3DT(LM,9))
        ALLOCATE(DT3DT(LM,6))
        ALLOCATE(DU3DT(LM,4))
        ALLOCATE(DV3DT(LM,4))
        ALLOCATE(PHY_F3DV(LM,4))
!
        CALL ESMF_ClockGet(clock       =CLOCK_ATM                       &  !<-- The ESMF Clock
                          ,currTime    =CURRTIME                        &  !<-- The current time (ESMF) on the clock
                          ,rc          =RC)
!
        CALL ESMF_TimeGet(time=CURRTIME                                 &  !<-- The cuurent forecast time (ESMF)
                         ,yy  =JDAT(1)                                  &  !<-- The current forecast year (integer)
                         ,mm  =JDAT(2)                                  &  !<-- The current forecast month (integer)
                         ,dd  =JDAT(3)                                  &  !<-- The current forecast day (integer)
                         ,h   =JDAT(5)                                  &  !<-- The current forecast hour (integer)
                         ,m   =JDAT(6)                                  &  !<-- The current forecast minute (integer)
                         ,s   =JDAT(7)                                  &  !<-- The current forecast second (integer)
                         ,rc  =RC)
        JDAT(4)=0
        JDAT(8)=0
!
        DO J=JTS,JTE
          GLOBAL_LATS_R(J) = J-JTS+1
          LONSPERLAR(J)    = ITE-ITS+1
          SINLAT_R(J)      = SIN(int_state%GLAT( (ITS+ITE)/2 ,J))
          COSLAT_R(J)      = SQRT( 1.d0 - SINLAT_R(J)*SINLAT_R(J) )
          DO I=ITS,ITE
            XLON(I,J)        = int_state%GLON(I,J)
            IF(int_state%GLON(I,J)<0) &
             XLON(I,J)        = 2.0d0*3.14159d0+XLON(I,J)
            COSZEN(I,J)      = int_state%CZEN(I,J)
            COSZDG(I,J)      = int_state%CZMEAN(I,J)
          ENDDO
        ENDDO
!
!-----------------------------------------------------------------------
!***  GFS Radiation
!-----------------------------------------------------------------------
!
        CALL_GFS_PHY = MOD(NTIMESTEP,int_state%NPHS)==0

        FHSWR        = FLOAT(int_state%NRADS)*int_state%DT/3600.   ! [h]
        LSCCA        = MOD(NTIMESTEP+1,int_state%NRADS)==0         ! logical true during a step for which convective clouds
                                                                   ! are calculated from convective precipitation rates
        LSSWR        = MOD(NTIMESTEP,int_state%NRADS)==0
        LSLWR        = MOD(NTIMESTEP,int_state%NRADL)==0
!
!-----------------------------------------------------------------------
        lw_or_sw: IF (LSSWR .OR. LSLWR ) THEN
!-----------------------------------------------------------------------
!
          DO L=1,LM+1
            KFLIP=LM-L+2
            RSGM(KFLIP)=int_state%SGM(L)
          ENDDO
! 
          ICWP=0                  ! control flag for cloud generation schemes
          IF (NTCW > 0) ICWP = 1  ! 0: use diagnostic cloud scheme
                                  ! 1: use prognostic cloud scheme (default)
!
! ----
!rv - find IOVR_SW,IOVR_LW,isubc_sw, isubc_lw
          CALL RADINIT_gfs ( RSGM, LM, IFLIP, IDAT, JDAT, ICTM, ISOL, ICO2, &
                         IAER, IALB, IEMS, ICWP, NUM_P3D, 0, 0,             &
                         0, 0, MYPE, RADDT, FDAER )
! ----

          IF (NTOZ .LE. 0) THEN                ! Climatological Ozone
!
            IDAY   = JDAT(3)
            IMON   = JDAT(2)
            MIDMON = DAYS(IMON)/2 + 1
            CHANGE = FIRST .OR. ( (IDAY .EQ. MIDMON) .AND. (JDAT(5).EQ.0) )
!
            IF (CHANGE) THEN
              IF (IDAY .LT. MIDMON) THEN
                 K1OZ = MOD(IMON+10,12) + 1
                 MIDM = DAYS(K1OZ)/2 + 1
                 K2OZ = IMON
                 MIDP = DAYS(K1OZ) + MIDMON
              ELSE
                 K1OZ = IMON
                 MIDM = MIDMON
                 K2OZ = MOD(IMON,12) + 1
                 MIDP = DAYS(K2OZ)/2 + 1 + DAYS(K1OZ)
              ENDIF
            ENDIF
!
            IF (IDAY .LT. MIDMON) THEN
              ID = IDAY + DAYS(K1OZ)
            ELSE
              ID = IDAY
            ENDIF
!
            FACOZ = REAL (ID-MIDM) / REAL (MIDP-MIDM)
!
          ELSE
!
            K1OZ = 0
            K2OZ = 0
            FACOZ = 1.0D0
!
          ENDIF
!
          FLGMIN_L(1)     = 0.2D0      ! --- for ferrier (for now, any number)

          DO J=JTS,JTE
            DO I=ITS,ITE
              SINLAT_V(I,J) = SINLAT_R(J)
              COSLAT_V(I,J) = COSLAT_R(J)
            ENDDO
          ENDDO
          DO J=JTS,JTE
            NLNSP(J) = LONR
          ENDDO
         

! ----
          CALL ASTRONOMY                                                &
!  ---  inputs:
             ( SINLAT_V, COSLAT_V, XLON, FHSWR, JDAT,                   &
               LONR, LATS_NODE_R, NLNSP, LSSWR, MYPE,                   &
!  ---  outputs:
               int_state%SOLCON, int_state%SLAG, int_state%SDEC,        &
               int_state%CDEC, COSZEN, COSZDG )
!
!-----------------------------------------------------------------------
!
        ENDIF  lw_or_sw
!
!-----------------------------------------------------------------------
!
!---
        IF (FIRST) THEN
!
          SEED0 = JDAT(4) + JDAT(3) + JDAT(2) + JDAT(1)
          CALL RANDOM_SETSEED(SEED0)
          CALL RANDOM_NUMBER(WRK)
          SEED0 = SEED0 + NINT(WRK(1)*1000.0)
          FIRST = .FALSE.
!
        ENDIF
!---
        FHOUR=NTIMESTEP*int_state%DT/3600.d0
        ISEED = MOD(100.0*SQRT(FHOUR*3600),1.0d9) + 1 + SEED0
        CALL RANDOM_SETSEED(ISEED)
        CALL RANDOM_NUMBER(RANNUM)
        N=0
!
        DO J=JTS,JTE
        DO I=ITS,ITE
          N=N+1
          RANN(I,J) = RANNUM(N)
        ENDDO
        ENDDO
!---
        DTF=int_state%NPHS*int_state%DT
        DTP=int_state%NPHS*int_state%DT
!---
        SOLHR=MOD(FHOUR+START_HOUR,24.d0)
!---
!...  set switch for saving convective clouds
        IF(LSCCA.AND.LSSWR) THEN
          CLSTP=1100+MIN(FHSWR,FHOUR,99.d0)  !initialize,accumulate,convert
        ELSEIF(LSCCA) THEN
          CLSTP=0100+MIN(FHSWR,FHOUR,99.d0)  !accumulate,convert
        ELSEIF(LSSWR) THEN
          CLSTP=1100                         !initialize,accumulate
        ELSE
          CLSTP=0100                         !accumulate
        ENDIF
!---
!---- OZONE ------------------------------------------------------------
!
        IF(.NOT.ALLOCATED(OZPLOUT_V)) &
                           ALLOCATE (OZPLOUT_V(LEVOZP,        PL_COEFF))
        IF(.NOT.ALLOCATED(OZPLOUT  )) &
                           ALLOCATE (OZPLOUT  (LEVOZP,JTS:JTE,PL_COEFF))
!
        IDATE(1)=JDAT(5)
        IDATE(2)=JDAT(2)
        IDATE(3)=JDAT(3)
        IDATE(4)=JDAT(1)
!
        IF (NTOZ .GT. 0) THEN
          CALL OZINTERPOL(MYPE,LATS_NODE_R,LATS_NODE_R,IDATE,FHOUR,     &
                          int_state%JINDX1,int_state%JINDX2,            &
                          int_state%OZPLIN,OZPLOUT,int_state%DDY)
        ENDIF
!
!---- OZONE ------------------------------------------------------------
!-----------------------------------------------------------------------
!***  Set diagnostics to 0.
!-----------------------------------------------------------------------
!
        DT3DT=0.0d0
        DU3DT=0.0d0
        DV3DT=0.0d0
        DQ3DT=0.0d0
!
        CV (1) = 0.d0       !!!!! not in use if ntcw-1 > 0
        CVB(1) = 0.d0       !!!!! not in use if ntcw-1 > 0
        CVT(1) = 0.d0       !!!!! not in use if ntcw-1 > 0
!
!-----------------------------------------------------------------------
!***  Empty the radiation flux and precipitation arrays if it is time.
!-----------------------------------------------------------------------
!
        IF(MOD(NTIMESTEP,int_state%NRDLW)==0)THEN
!
          DO J=JTS,JTE
          DO I=ITS,ITE
            int_state%ALWIN(I,J) =0.
            int_state%ALWOUT(I,J)=0.
            int_state%ALWTOA(I,J)=0.
            int_state%ARDLW (I,J)=0.   !<-- An artificial 2-D array (ESMF
                                       !<-- cannot have evolving scalar Attributes)
          ENDDO
          ENDDO
!
        ENDIF
!
        IF(MOD(NTIMESTEP,int_state%NRDSW)==0)THEN
!
          DO J=JTS,JTE
          DO I=ITS,ITE
            int_state%ASWIN(I,J)=0.
            int_state%ASWOUT(I,J)=0.
            int_state%ASWTOA(I,J)=0.
            int_state%ARDSW (I,J)=0.   !<-- An artificial 2-D array (ESMF
                                       !<-- cannot have evolving scalar Attributes)
          ENDDO
          ENDDO
!
        ENDIF
!
        IF(MOD(NTIMESTEP,int_state%NSRFC)==0)THEN
!
          DO J=JTS,JTE
          DO I=ITS,ITE
            int_state%ACSNOW(I,J)=0.
            int_state%POTEVP(I,J)=0.
            int_state%SFCEVP(I,J)=0.
            int_state%SFCLHX(I,J)=0.
            int_state%SFCSHX(I,J)=0.
            int_state%SUBSHX(I,J)=0.
            int_state%BGROFF(I,J)=0.
            int_state%SSROFF(I,J)=0.
            int_state%ASRFC (I,J)=0.   !<-- An artificial 2-D array (ESMF
                                       !<-- cannot have evolving scalar Attributes)
          ENDDO
          ENDDO
!
        ENDIF
!
        IF(MOD(NTIMESTEP,int_state%NPREC)==0)THEN
          DO J=JTS,JTE
          DO I=ITS,ITE
            int_state%ACPREC(I,J)=0.
            int_state%CUPREC(I,J)=0.
          ENDDO
          ENDDO
        ENDIF
!
!-----------------------------------------------------------------------
        gfs_physics: IF(CALL_GFS_PHY)THEN
#if 1
        WRITE(0,*)' GFS physics option in NMMB disabled, 20140812, jm'
        CALL ESMF_FINALIZE(endflag=ESMF_END_ABORT                   &
                          ,rc             =RC)
#else
!-----------------------------------------------------------------------
!
          DTLW     = FLOAT(int_state%NRADL)*int_state%DT   ! [s]
          DTSW     = FLOAT(int_state%NRADS)*int_state%DT   ! [s]
          RADDT    = MIN(DTSW,DTLW)
          DTLWI    = 1./DTLW
          DTSWI    = 1./DTSW
          MINDT    = 1./MIN(DTLW,DTSW)
          DTPHS    = int_state%NPHS*int_state%DT
          DTPHSI   = 1./DTPHS
          XLVRW    = XLV*RHOWATER
          XLVRWI   = 1./XLVRW
          RoCP     = R/CP
          LSSAV_CC = LSSAV
          LSSAV_CPL= LSSAV
          IDE_GR   = IDE-1
          IF(int_state%GLOBAL) IDE_GR = IDE-3
!
!-----------------------------------------------------------------------
! ***  MAIN GFS-PHYS DOMAIN LOOP
!-----------------------------------------------------------------------
!
          j_loop: DO J=JTS,JTE
!
            i_loop: DO I=ITS,ITE
!
              int_state%ACUTIM(I,J) = int_state%ACUTIM(I,J) + 1.     ! advance counters
              int_state%APHTIM(I,J) = int_state%APHTIM(I,J) + 1.
              int_state%ARDLW(I,J)  = int_state%ARDLW(I,J)  + 1.
              int_state%ARDSW(I,J)  = int_state%ARDSW(I,J)  + 1.
              int_state%ASRFC(I,J)  = int_state%ASRFC(I,J)  + 1.
              int_state%AVRAIN(I,J) = int_state%AVRAIN(I,J) + 1.
              int_state%AVCNVC(I,J) = int_state%AVCNVC(I,J) + 1.
!
              T1(1)           = 0.0D0        ! initialize all local variables
              Q1(1)           = 0.0D0        ! used in gfs_physics
              U1(1)           = 0.0D0
              V1(1)           = 0.0D0
              QSS(1)          = 0.0D0
              CDQ(1)          = 0.0D0
              HFLX(1)         = 0.0D0
              EVAP(1)         = 0.0D0
              DTSFC(1)        = 0.0D0
              DQSFC(1)        = 0.0D0
              DUSFC(1)        = 0.0D0
              DVSFC(1)        = 0.0D0
              PSMEAN(1)       = 0.0D0
              EPI(1)          = 0.0D0
              EVBSA(1)        = 0.0D0
              EVCWA(1)        = 0.0D0
              TRANSA(1)       = 0.0D0
              SBSNOA(1)       = 0.0D0
              SOILM(1)        = 0.0D0
              SNOWCA(1)       = 0.0D0
              DLWSFC_CC(1)    = 0.0D0
              ULWSFC_CC(1)    = 0.0D0
              DTSFC_CC(1)     = 0.0D0
              SWSFC_CC(1)     = 0.0D0
              DUSFC_CC(1)     = 0.0D0
              DVSFC_CC(1)     = 0.0D0
              DQSFC_CC(1)     = 0.0D0
              PRECR_CC(1)     = 0.0D0
              DUSFC_CPL(1)    = 0.0D0
              DVSFC_CPL(1)    = 0.0D0
              DTSFC_CPL(1)    = 0.0D0
              DQSFC_CPL(1)    = 0.0D0
              DLWSFC_CPL(1)   = 0.0D0
              DSWSFC_CPL(1)   = 0.0D0
              DNIRBM_CPL(1)   = 0.0D0
              DNIRDF_CPL(1)   = 0.0D0
              DVISBM_CPL(1)   = 0.0D0
              DVISDF_CPL(1)   = 0.0D0
              NLWSFC_CPL(1)   = 0.0D0
              NSWSFC_CPL(1)   = 0.0D0
              NNIRBM_CPL(1)   = 0.0D0
              NNIRDF_CPL(1)   = 0.0D0
              NVISBM_CPL(1)   = 0.0D0
              NVISDF_CPL(1)   = 0.0D0
              RAIN_CPL(1)     = 0.0D0
              XT(1)           = 0.0D0
              XS(1)           = 0.0D0
              XU(1)           = 0.0D0
              XV(1)           = 0.0D0
              XZ(1)           = 0.0D0
              ZM(1)           = 0.0D0
              XTTS(1)         = 0.0D0
              XZTS(1)         = 0.0D0
              D_CONV(1)       = 0.0D0
              IFD(1)          = 0.0D0
              DT_COOL(1)      = 0.0D0
              QRAIN(1)        = 0.0D0
              XMU_CC(1)       = 0.0D0
              DLW_CC(1)       = 0.0D0
              DSW_CC(1)       = 0.0D0
              SNW_CC(1)       = 0.0D0
              LPREC_CC(1)     = 0.0D0
              DUSFCI_CPL(1)   = 0.0D0
              DVSFCI_CPL(1)   = 0.0D0
              DTSFCI_CPL(1)   = 0.0D0
              DQSFCI_CPL(1)   = 0.0D0
              DLWSFCI_CPL(1)  = 0.0D0
              DSWSFCI_CPL(1)  = 0.0D0
              DNIRBMI_CPL(1)  = 0.0D0
              DNIRDFI_CPL(1)  = 0.0D0
              DVISBMI_CPL(1)  = 0.0D0
              DVISDFI_CPL(1)  = 0.0D0
              NLWSFCI_CPL(1)  = 0.0D0
              NSWSFCI_CPL(1)  = 0.0D0
              NNIRBMI_CPL(1)  = 0.0D0
              NNIRDFI_CPL(1)  = 0.0D0
              NVISBMI_CPL(1)  = 0.0D0
              NVISDFI_CPL(1)  = 0.0D0
              T2MI_CPL(1)     = 0.0D0
              Q2MI_CPL(1)     = 0.0D0
              U10MI_CPL(1)    = 0.0D0
              V10MI_CPL(1)    = 0.0D0
              TSEAI_CPL(1)    = 0.0D0
              PSURFI_CPL(1)   = 0.0D0
              ORO_CPL(1)      = 0.0D0
              SLMSK_CPL(1)    = 0.0D0
              TREF(1)         = 0.0D0
              Z_C(1)          = 0.0D0
              C_0(1)          = 0.0D0
              C_D(1)          = 0.0D0
              W_0(1)          = 0.0D0
              W_D(1)          = 0.0D0
              RQTK(1)         = 0.0D0
              SNOHFA(1)       = 0.0D0
              SMCWLT2(1)      = 0.0D0
              SMCREF2(1)      = 0.0D0
              WET1(1)         = 0.0D0
              GSOIL(1)        = 0.0D0
              GTMP2M(1)       = 0.0D0
              GUSTAR(1)       = 0.0D0
              GPBLH(1)        = 0.0D0
              GU10M(1)        = 0.0D0
              GV10M(1)        = 0.0D0
              GZORL(1)        = 0.0D0
              GORO(1)         = 0.0D0
              SR(1)           = 0.0D0
              SPFHMIN(1)      = 0.0D0
              SPFHMAX(1)      = 0.0D0
              CLDWRK(1)       = 0.0D0
              ZLVL(1)         = 0.0D0
              PHII            = 0.0D0
              PHIL            = 0.0D0
              CHH(1)          = 0.0D0
              HPBL(1)         = 0.0D0
              PSURF(1)        = 100000.0D0
              T2M(1)          = 273.0D0
              Q2M(1)          = 0.0D0
              U10M(1)         = 0.0D0
              V10M(1)         = 0.0D0
              ADR             = 0.0D0
              ADT             = 0.0D0
              ADU             = 0.0D0
              ADV             = 0.0D0
              SUNTIM          = 0.0D0
              ICSDSW(1)       = 0
              ICSDLW(1)       = 0

         IF(int_state%TSKIN(I,J) .LT. 50. ) THEN
             TSEA(1)         = int_state%SST(I,J)
             TISFC(1)        = int_state%SST(I,J)
         ELSE
             TSEA(1)         = int_state%TSKIN(I,J)
             TISFC(1)        = int_state%TSKIN(I,J)
         ENDIF

         IF(int_state%SICE(I,J) > 0.5 ) THEN                                ! slmsk - ocean  - 0
             SLMSK(1)        = 2.0D0                                        !         land   - 1
         ELSE                                                               !         seaice - 2
             SLMSK(1)        = 1.0D0-int_state%SM(I,J)                      !
         ENDIF

         DO L=1,LM
            KFLIP=LM+1-L
!            CLDCOV_V(KFLIP) = 0.0D0                      ! GRRAD now returns instant cloud cover (Sarah Lu)
             F_ICE(KFLIP)    = int_state%F_ICE(I,J,L)                       ! for ferrier phy, do init first
             F_RAIN(KFLIP)   = int_state%F_RAIN(I,J,L)
             R_RIME(KFLIP)   = int_state%F_RIMEF(I,J,L)
         ENDDO

             XLAT(1)         = int_state%GLAT(I,J)
             ZORL(1)         = int_state%ZORFCS(I,J)
             SNCOVR(1)       = int_state%SNO(I,J)/(int_state%SNO(I,J)+70.)  ! FORMULATION OF MARSHALL ET AL. 1994
                                                                            ! change this later only initially, add new int_state
             SNWDPH(1)       = int_state%SI(I,J)                            ! snwdph[mm]
             WEASD(1)        = int_state%SNO(I,J)                           ! snow water eq.[mm]
             SNOALB(1)       = int_state%MXSNAL(I,J)
             ALVSF(1)        = int_state%ALBFC1(I,J,1)                      ! VIS, direct
             ALVWF(1)        = int_state%ALBFC1(I,J,2)                      ! VIS, diffuse
             ALNSF(1)        = int_state%ALBFC1(I,J,3)                      ! NIR, direct
             ALNWF(1)        = int_state%ALBFC1(I,J,4)                      ! NIR, diffuse
             FACSF(1)        = int_state%ALFFC1(I,J,1)                      ! direct
             FACWF(1)        = int_state%ALFFC1(I,J,2)                      ! diffuse
!
             PRSI (LM+1)     = int_state%PT                                 ! [ Pa]
             PRSIK(LM+1)     = (PRSI(LM+1)*0.00001d0)**RoCP
         DO L=1,LM
            KFLIP=LM+1-L
             PRSI (KFLIP)    = PRSI(KFLIP+1) + &
                                (int_state%DSG2(L)*int_state%PD(I,J)+int_state%PDSG1(L))  ! (pressure on interface) [ Pa]
             PRSIK(KFLIP)    = (PRSI(KFLIP)*0.00001d0)**RoCP

             PRSL (KFLIP)    = (PRSI(KFLIP)+PRSI(KFLIP+1))*0.5d0             ! (pressure on mid-layer) [kPa]
             PRSLK(KFLIP)    = (PRSL(KFLIP)*0.00001d0)**RoCP
!
             RTvR = 1. / ( R * (int_state%Q(I,J,L)*0.608+1.-int_state%CW(I,J,L) ) * int_state%T(I,J,L) )
             VVEL(KFLIP)     = int_state%OMGALF(I,J,L) * PRSL(KFLIP) * RTvR
!
             GU(KFLIP)       = (int_state%U(I,J  ,L) + int_state%U(I-1,J  ,L) +                    &
                                int_state%U(I,J-1,L) + int_state%U(I-1,J-1,L))*0.25d0
             GV(KFLIP)       = (int_state%V(I,J  ,L) + int_state%V(I-1,J  ,L) +                    &
                                int_state%V(I,J-1,L) + int_state%V(I-1,J-1,L))*0.25d0
             GT(KFLIP)       = int_state%T(I,J,L)
             GR(KFLIP)       = int_state%Q(I,J,L)
             GR3(KFLIP,1)    = int_state%Q(I,J,L)
           IF (NTIMESTEP == 0 ) THEN
             GR3(KFLIP,2)    = 0.0d0
             GR3(KFLIP,3)    = 0.0d0
           ELSE
             GR3(KFLIP,2)    = int_state%O3(I,J,L)
             GR3(KFLIP,3)    = int_state%CW(I,J,L)
           ENDIF
             GR1(1,KFLIP,1)  = GR3(KFLIP,2)
             GR1(1,KFLIP,2)  = int_state%CW(I,J,L)
         ENDDO
!---
             DLWSFC(1)       = int_state%ALWIN(I,J)
             ULWSFC(1)       = int_state%ALWOUT(I,J)
             DLWSFCI(1)      = int_state%RLWIN(I,J)
             ULWSFCI(1)      = int_state%RADOT(I,J)
             DSWSFCI(1)      = int_state%RSWIN(I,J)
             USWSFCI(1)      = int_state%RSWOUT(I,J)
!---
             GFLUX(1)        = 0.0D0
             DQSFCI(1)       = 0.0D0
             DTSFCI(1)       = 0.0D0
             GFLUXI(1)       = 0.0D0
             EP(1)           = int_state%POTEVP(I,J)*XLVRW
!---
             XSIHFCS(1)      = int_state%SIHFCS(I,J)
             XSICFCS(1)      = int_state%SICFCS(I,J)
             XSLPFCS(1)      = int_state%SLPFCS(I,J)
             XTG3FCS(1)      = int_state%TG3FCS(I,J)
             XVEGFCS(1)      = int_state%VEGFCS(I,J)
             XVETFCS(1)      = int_state%VETFCS(I,J)
             XSOTFCS(1)      = int_state%SOTFCS(I,J)
!---
             FLUXR_V         = 0.0D0
             IF(.NOT.LSLWR) FLUXR_V(1) = int_state%RLWTOA(I,J)*DTLW
             IF(.NOT.LSSWR) FLUXR_V(2) = int_state%RSWTOA(I,J)*DTSW
!---
             HPRIME (1)      = int_state%HSTDV(I,J)
             HPRIME (2)      = int_state%HCNVX(I,J)
             HPRIME (3)      = int_state%HASYW(I,J)
             HPRIME (4)      = int_state%HASYS(I,J)
             HPRIME (5)      = int_state%HASYSW(I,J)
             HPRIME (6)      = int_state%HASYNW(I,J)
             HPRIME (7)      = int_state%HLENW(I,J)
             HPRIME (8)      = int_state%HLENS(I,J)
             HPRIME (9)      = int_state%HLENSW(I,J)
             HPRIME(10)      = int_state%HLENNW(I,J)
             HPRIME(11)      = int_state%HANGL(I,J)*180.D0/3.14159D0
             HPRIME(12)      = int_state%HANIS(I,J)
             HPRIME(13)      = int_state%HSLOP(I,J)
             HPRIME(14)      = int_state%HZMAX(I,J)
!---
             RUNOFF(1)       = int_state%BGROFF(I,J)*0.001D0
             SRUNOFF(1)      = int_state%SSROFF(I,J)*0.001D0
!---
           DO L=1,LM
            KFLIP=LM+1-L
             DKH(L)          = 0.0D0
             RNP(L)          = 0.0D0
             SWH(KFLIP)      = int_state%RSWTT(I,J,L)
             HLW(KFLIP)      = int_state%RLWTT(I,J,L)
           ENDDO
           DO N=1,3                                    ! for Zhao =3, Ferr=1
             PHY_F2DV(N)     = int_state%PHY_F2DV (I,J,N)
           ENDDO
           DO N=1,4                                    ! for Zhao =4, Ferr=3
           DO L=1,LM
            PHY_F3DV(L,N)    = int_state%PHY_F3DV (I,J,L,N)
           ENDDO
           ENDDO
!
!-----------------------------------------------------------------------
          CALL GRRAD_gfs                                             &
!-----------------------------------------------------------------------
!  ---  inputs:
          (PRSI,PRSL,PRSLK,GT,GR,GR1,VVEL,SLMSK,                     &
           XLON(I,J),XLAT,TSEA,SNWDPH,SNCOVR,SNOALB,ZORL,HPRIME(1),  &
           ALVSF,ALNSF,ALVWF,ALNWF,FACSF,FACWF,XSICFCS,TISFC,        &
           int_state%SOLCON,COSZEN(I,J),COSZDG(I,J),K1OZ,K2OZ,FACOZ, &
           CV,CVT,CVB,IOVR_SW,IOVR_LW,F_ICE,F_RAIN,R_RIME,FLGMIN_L,  &
           ICSDSW,ICSDLW,NUM_P3D,NTCW-1,NCLD,NTOZ-1,NTRAC-1,NFXR,    &
           DTLW,DTSW,LSSWR,LSLWR,LSSAV,SASHAL,NORAD_PRECIP,          &
           CRICK_PROOF,CCNORM,                                       &
           1,1,LM,IFLIP,MYPE,LPRNT,1,NTIMESTEP,                      &
!  ---  outputs:
           SWH,TOPFSW,SFCFSW,SFALB,                                  &
           HLW,TOPFLW,SFCFLW,TSFLW,SEMIS,CLDCOV_V,                   &
!  ---  input/output:
           FLUXR_V                                                   &
          )
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  GBPHYS   !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---
        IF (LSSWR .OR. LSLWR ) THEN
          SFCDLW(1)             = SFCFLW(1)%DNFXC
          SFCDSW(1)             = SFCFSW(1)%DNFXC
          SFCNSW(1)             = SFCFSW(1)%DNFXC - SFCFSW(1)%UPFXC
          int_state%SFCDLW(I,J) = SFCDLW(1)
          int_state%SFCDSW(I,J) = SFCDSW(1)
          int_state%SFCNSW(I,J) = SFCNSW(1)
          int_state%SFALB(I,J)  = SFALB(1)
          int_state%TSFLW(I,J)  = TSFLW(1)
          int_state%SEMIS(I,J)  = SEMIS(1)
        ELSE
          SFCDLW(1)   = int_state%SFCDLW(I,J)
          SFCDSW(1)   = int_state%SFCDSW(I,J)
          SFCNSW(1)   = int_state%SFCNSW(I,J)
          SFALB(1)    = int_state%SFALB(I,J)
          TSFLW(1)    = int_state%TSFLW(I,J)
          SEMIS(1)    = int_state%SEMIS(I,J)
        ENDIF
!---
          NIRBMD_CPL(1) = 0.25 * SFCDSW(1)
          NIRDFD_CPL(1) = 0.26 * SFCDSW(1)
          VISBMD_CPL(1) = 0.16 * SFCDSW(1)
          VISDFD_CPL(1) = 0.33 * SFCDSW(1)
          NIRBMU_CPL(1) = NIRBMD_CPL(1) - 0.25 * SFCNSW(1)
          NIRDFU_CPL(1) = NIRDFD_CPL(1) - 0.26 * SFCNSW(1)
          VISBMU_CPL(1) = VISBMD_CPL(1) - 0.16 * SFCNSW(1)
          VISDFU_CPL(1) = VISDFD_CPL(1) - 0.33 * SFCNSW(1)
!---
          DPSHC(1)    = 0.3 * PRSI(1)
          GQ(1)       = PRSI(1)
!---
          RANNUM_V(1) = RANN(I,J)
!---
          RCS2_V(1)   = 1.0d0/(0.0001d0+COS(XLAT(1))*COS(XLAT(1))) ! fixed XLAT=+/-90
!---
          TPRCP(1)    = int_state%PREC(I,J)
          CNVPRCP(1)  = int_state%CUPREC(I,J)
          TOTPRCP(1)  = int_state%ACPREC(I,J)

         DO L=1,NUM_SOIL_LAYERS
          SMC_V(L)    = int_state%SMC(I,J,L)
          STC_V(L)    = int_state%STC(I,J,L)
          SLC_V(L)    = int_state%SH2O(I,J,L)
         ENDDO

         SHDMIN(1)   = int_state%SHDMIN(I,J)
         SHDMAX(1)   = int_state%SHDMAX(I,J)

         UUSTAR(1)   = int_state%USTAR(I,J)
         CANOPY(1)   = int_state%CMC(I,J)*1000.
         FLGMIN      = 0.2d0 !!! put this in ferrier init
         CRTRH       = 0.85d0
         FLIPV       = .FALSE.
         NCW(1)      = 50
         NCW(2)      = 150
         OLD_MONIN   = .FALSE.
         CNVGWD      = .FALSE.
         NEWSAS      = 0
         CCWF        = 0.5d0  ! only for RAS scheme
         CGWF        = 0.1    ! cloud top fraction for convective gwd scheme
         PRSLRD0     = 0.     ! pressure(pa) above which Raleigh damping applied
!
! ---- ESTIMATE T850 FOR RAIN-SNOW DECISION ----------------------------
!
          T850 = GT(1)
      DO L = 1, LM - 1
        IF(PRSL(L) .GT. 85000.d0 .AND. PRSL(L+1) .LE. 85000.d0) THEN
          T850 = GT(L) - (PRSL(L)-85000.d0) / (PRSL(L)-PRSL(L+1)) * (GT(L)-GT(L+1))
        ENDIF
      ENDDO
!
      SRFLAG(1) = 0.0d0
      IF(T850 .LE. 273.16d0) SRFLAG(1) = 1.0d0
!
!---- OZONE ------------------------------------------------------------
      IF (NTOZ .GT. 0) THEN
        DO N=1,PL_COEFF
          DO L=1,LEVOZP
              OZPLOUT_V(L,N) = OZPLOUT(L,j,N)
          ENDDO
        ENDDO
      ELSE
              OZPLOUT_V      = 0.0d0
      ENDIF
      DTDT = 0.0D0
      TRIGGERPERTS(1) = 0.0d0
!-----------------------------------------------------------------------
      CALL GBPHYS(1, 1, LM, NUM_SOIL_LAYERS, LSM, NTRAC, NCLD, NTOZ, NTCW,  &
           NMTVR, 1, LEVOZP, IDE_GR, LATR, 62, NUM_P3D, NUM_P2D,            &
           NTIMESTEP, J-JTS+1, MYPE, PL_COEFF, LONR, NCW, FLGMIN, CRTRH,    &
             CDMBGWD, &
           CCWF, DLQF, CTEI_RM, CLSTP, CGWF, PRSLRD0, DTP, DTF, FHOUR,      &
           SOLHR, int_state%SLAG, int_state%SDEC, int_state%CDEC,           &
           SINLAT_R(J), COSLAT_R(J), GQ, GU, GV,                            &
           GT, GR3, VVEL, PRSI, PRSL, PRSLK, PRSIK, PHII, PHIL,             &
           RANN, OZPLOUT_V, PL_PRES, DPSHC, HPRIME, XLON(I,J), XLAT,        &
           XSLPFCS, SHDMIN, SHDMAX, SNOALB, XTG3FCS, SLMSK, XVEGFCS,        &
           XVETFCS, XSOTFCS, UUSTAR, ORO, oro, COSZEN(I,J), SFCDSW, SFCNSW, &
           NIRBMD_CPL,       NIRDFD_CPL,  VISBMD_CPL,       VISDFD_CPL,     &
           NIRBMU_CPL,       NIRDFU_CPL,  VISBMU_CPL,       VISDFU_CPL,     &
           SFCDLW, TSFLW, SEMIS, SFALB, SWH, HLW, RAS, PRE_RAD,             &
           LDIAG3D, LGGFS3D, LGOCART, LSSAV, LSSAV_CC, LSSAV_CPL,           &
           BKGD_VDIF_M, BKGD_VDIF_H, BKGD_VDIF_S, PSAUTCO, PRAUTCO, EVPCO,  &
           WMINCO,                                                          &
           FLIPV, OLD_MONIN, CNVGWD, SHAL_CNV, SASHAL, NEWSAS, CAL_PRE,     &
           MOM4ICE, MSTRAT, TRANS_TRAC, NST_FCST, MOIST_ADJ,                &
           THERMODYN_ID, SFCPRESS_ID, GEN_COORD_HYBRID, LEVR,               &
           XSIHFCS, XSICFCS, TISFC, TSEA, TPRCP, CV, CVB, CVT,              &
           SRFLAG, SNWDPH, WEASD, SNCOVR, ZORL, CANOPY,                     &
           FFMM, FFHH, F10M, SRUNOFF, EVBSA, EVCWA, SNOHFA,                 &
           TRANSA, SBSNOA, SNOWCA, SOILM, int_state%TMPMIN(I,J),            &
           int_state%TMPMAX(I,J),                                           &
           DUSFC, DVSFC, DTSFC, DQSFC, TOTPRCP, GFLUX,                      &
           DLWSFC, ULWSFC, SUNTIM, RUNOFF, EP, CLDWRK,                      &
           int_state%DUGWD(I,J), int_state%DVGWD(I,J), PSMEAN, CNVPRCP,     &
           SPFHMIN, SPFHMAX, RAIN, RAINC,                                   &
           DT3DT, DQ3DT, DU3DT, DV3DT, DQDT, ACV, ACVB, ACVT,               &
           SLC_V, SMC_V, STC_V, UPD_MF, DWN_MF, DET_MF, DKH, RNP, PHY_F3DV, &
             PHY_F2DV,                                                      &
           DLWSFC_CC, ULWSFC_CC, DTSFC_CC, SWSFC_CC,                        &
           DUSFC_CC, DVSFC_CC, DQSFC_CC, PRECR_CC,                          &
           DUSFC_CPL, DVSFC_CPL, DTSFC_CPL, DQSFC_CPL,                      &
           DLWSFC_CPL,DSWSFC_CPL,DNIRBM_CPL,DNIRDF_CPL,                     &
           DVISBM_CPL,DVISDF_CPL,RAIN_CPL,                                  &
           NLWSFC_CPL,NSWSFC_CPL,NNIRBM_CPL,NNIRDF_CPL,                     &
           NVISBM_CPL,NVISDF_CPL,                                           &
           XT, XS, XU, XV, XZ, ZM, XTTS, XZTS, D_CONV, IFD, DT_COOL, QRAIN, &
           ADT, ADR, ADU, ADV, T2M, Q2M, U10M, V10M,                        &
           ZLVL, PSURF, HPBL, PWAT, T1, Q1, U1, V1,                         &
           CHH, CMM, DLWSFCI, ULWSFCI, DSWSFCI, USWSFCI, DUSFCI, DVSFCI,    &
           DTSFCI, DQSFCI, GFLUXI, EPI, SMCWLT2, SMCREF2, WET1,             &
           GSOIL, GTMP2M, GUSTAR, GPBLH, GU10M, GV10M, GZORL, GORO, SR,     &
           XMU_CC, DLW_CC, DSW_CC, SNW_CC, LPREC_CC,                        &
           DUSFCI_CPL, DVSFCI_CPL, DTSFCI_CPL, DQSFCI_CPL,                  &
           DLWSFCI_CPL,DSWSFCI_CPL,DNIRBMI_CPL,DNIRDFI_CPL,                 &
           DVISBMI_CPL,DVISDFI_CPL,                                         &
           NLWSFCI_CPL,NSWSFCI_CPL,NNIRBMI_CPL,NNIRDFI_CPL,                 &
           NVISBMI_CPL,NVISDFI_CPL,T2MI_CPL,Q2MI_CPL,                       &
           U10MI_CPL,V10MI_CPL,TSEAI_CPL,PSURFI_CPL,ORO_CPL,SLMSK_CPL,      &
           TREF, Z_C, C_0, C_D, W_0, W_D, RQTK, HLWD, LSIDEA,               &
           DTDT, TRIGGERPERTS)
!-----------------------------------------------------------------------
! ***     UPDATE AFTER PHYSICS
!-----------------------------------------------------------------------
             int_state%SIHFCS(I,J)        = XSIHFCS(1)
             int_state%SICFCS(I,J)        = XSICFCS(1)
             int_state%ZORFCS(I,J)        = ZORL(1)

             int_state%CZEN(I,J)          = COSZEN(I,J)
             int_state%CZMEAN(I,J)        = COSZDG(I,J)

             int_state%SI(I,J)            = SNWDPH(1)
             int_state%SNO(I,J)           = WEASD(1)
             int_state%MXSNAL(I,J)        = SNOALB(1)

         DO L=1,NUM_SOIL_LAYERS
             int_state%SMC(I,J,L)         = SMC_V(L)
             int_state%STC(I,J,L)         = STC_V(L)
             int_state%SH2O(I,J,L)        = SLC_V(L)
         ENDDO

             int_state%CMC(I,J)           = CANOPY(1)*0.001
             int_state%USTAR(I,J)         = UUSTAR(1)
             int_state%SMSTOT(I,J)        = SOILM(1)*1000.

         DO L=1,LM
            KFLIP=LM+1-L
             int_state%RSWTT(I,J,L)       = SWH(KFLIP)
             int_state%RLWTT(I,J,L)       = HLW(KFLIP)
         ENDDO

         DO N=1,3                                    ! for Zhao =3, Ferr=1
             int_state%PHY_F2DV (I,J,N)   = PHY_F2DV(N)
         ENDDO

         DO N=1,4                                    ! for Zhao =4, Ferr=3
         DO L=1,LM
             int_state%PHY_F3DV (I,J,L,N) = PHY_F3DV(L,N)
         ENDDO
         ENDDO

             int_state%ALWIN(I,J)         = int_state%ALWIN(I,J)  + int_state%RLWIN(I,J)
             int_state%ALWOUT(I,J)        = int_state%ALWOUT(I,J) - int_state%RADOT(I,J)
             int_state%ASWIN(I,J)         = int_state%ASWIN(I,J)  + int_state%RSWIN(I,J)
             int_state%ASWOUT(I,J)        = int_state%ASWOUT(I,J) - int_state%RSWOUT(I,J)
             int_state%RLWIN(I,J)         = DLWSFCI(1)
             int_state%RADOT(I,J)         = ULWSFCI(1)
             int_state%RSWIN(I,J)         = DSWSFCI(1)
             int_state%RSWOUT(I,J)        = USWSFCI(1)
             int_state%RSWINC(I,J)        = int_state%RSWIN(I,J)/(1.-int_state%ALBEDO(I,J))

             int_state%RLWTOA(I,J)        = FLUXR_V(1)*DTLWI
             int_state%RSWTOA(I,J)        = FLUXR_V(2)*DTSWI
             int_state%ALWTOA(I,J)        = int_state%ALWTOA(I,J) + FLUXR_V(1)*DTLWI
             int_state%ASWTOA(I,J)        = int_state%ASWTOA(I,J) + FLUXR_V(2)*DTSWI

             int_state%TWBS(I,J)          = -DQSFCI(1)
             int_state%QWBS(I,J)          = -DTSFCI(1)
             int_state%SFCSHX(I,J)        = int_state%SFCSHX(I,J) + int_state%QWBS(I,J)
             int_state%SFCLHX(I,J)        = int_state%SFCLHX(I,J) + int_state%TWBS(I,J)
             int_state%SUBSHX(I,J)        = int_state%SUBSHX(I,J) + GFLUXI(1)
             int_state%GRNFLX(I,J)        = GFLUXI(1)
             int_state%POTEVP(I,J)        = EP(1)*XLVRWI
             int_state%POTFLX(I,J)        = -EP(1)*DTPHSI
             int_state%SFCEVP(I,J)        = int_state%SFCEVP(I,J) + DQSFCI(1)*DTPHS*XLVRWI

             int_state%SFCEXC(I,J)        = CDQ(1)                                       !need CDQ  from GFS
             int_state%PBLH(I,J)          = HPBL(1)
             int_state%PSFC(I,J)          = PSURF(1)
             int_state%PREC(I,J)          = TPRCP(1)
             int_state%CUPPT(I,J)         = CNVPRCP(1)-int_state%CUPREC(I,J)
             int_state%CPRATE(I,J)        = CNVPRCP(1)-int_state%CUPREC(I,J)
             int_state%CUPREC(I,J)        = CNVPRCP(1)
             int_state%ACPREC(I,J)        = TOTPRCP(1)

             int_state%BGROFF(I,J)        = RUNOFF(1)*1000.
             int_state%SSROFF(I,J)        = SRUNOFF(1)*1000.

             int_state%TSKIN(I,J)         = TISFC(1)
             int_state%SST(I,J)           = TSEA(1)
             int_state%SOILTB(I,J)        = XTG3FCS(1)
             IF( SRFLAG(1) >= 0.5 .AND. SLMSK(1) >= 0.5 ) &
             int_state%ACSNOW(I,J)        = int_state%ACSNOW(I,J) + int_state%ACPREC(I,J)*100.

             int_state%PSHLTR(I,J)        = PSURF(1)*EXP(0.06823/T2M(1))
             int_state%TSHLTR(I,J)        = T2M(1)
             int_state%QSHLTR(I,J)        = Q2M(1)
             int_state%QSH(I,J)           = QSS(1)                                    !need QSS  from GFS
             int_state%T2(I,J)            = T2M(1)
             int_state%TH02(I,J)          = T2M(1)*(100000./PSURF(1))**RoCP
             int_state%Q02(I,J)           = Q2M(1)
             int_state%U10(I,J)           = U10M(1)
             int_state%V10(I,J)           = V10M(1)
             int_state%THS(I,J)           = TSFLW(1)*(100000./PSURF(1))**RoCP
             int_state%SIGT4(I,J)         = int_state%T(I,J,LM)*int_state%T(I,J,LM) * &
                                            int_state%T(I,J,LM)*int_state%T(I,J,LM) * STBOLT
         DO L=1,LM
            KFLIP=LM+1-L
             int_state%T(I,J,L)           = ADT(KFLIP)
             int_state%DUDT(I,J,L)        = (ADU(KFLIP) - GU(KFLIP)) / DTP
             int_state%DVDT(I,J,L)        = (ADV(KFLIP) - GV(KFLIP)) / DTP
             int_state%CLDFRA(I,J,L)      = CLDCOV_V(KFLIP) 
             int_state%Q (I,J,L)          = ADR(KFLIP,1)
             int_state%O3(I,J,L)          = ADR(KFLIP,2)
             int_state%CW(I,J,L)          = ADR(KFLIP,3)
         ENDDO
!
!-----------------------------------------------------------------------
! ***     End update after Physics
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
! ***  End GFS-PHYS domain loop
!-----------------------------------------------------------------------
!
            ENDDO  i_loop
!
          ENDDO    j_loop
!
!-----------------------------------------------------------------------
!
          td%gfs_phy_tim=td%gfs_phy_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Exchange wind tendencies
!-----------------------------------------------------------------------
!
          btim=timef()
!
          CALL HALO_EXCH(int_state%DUDT,LM,int_state%DVDT,LM,3,3)
!
          td%exch_phy=td%exch_phy+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Now interpolate wind tendencies from H to V points.
!-----------------------------------------------------------------------
!
          btim=timef()
!
          CALL H_TO_V_TEND(int_state%DUDT,int_state%DT,int_state%NPHS,LM &
                          ,int_state%U)
          CALL H_TO_V_TEND(int_state%DVDT,int_state%DT,int_state%NPHS,LM &
                          ,int_state%V)
!
          td%h_to_v_tim=td%h_to_v_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Poles and East-West boundary.
!-----------------------------------------------------------------------
!
          IF(int_state%GLOBAL)THEN
            btim=timef()
!
            CALL SWAPHN(int_state%T,IMS,IME,JMS,JME,LM,int_state%INPES)
            CALL POLEHN(int_state%T,IMS,IME,JMS,JME,LM                  &
                       ,int_state%INPES,int_state%JNPES)
!
            CALL SWAPHN(int_state%Q,IMS,IME,JMS,JME,LM,int_state%INPES)
            CALL POLEHN(int_state%Q,IMS,IME,JMS,JME,LM                  &
                       ,int_state%INPES,int_state%JNPES)
!
            CALL SWAPHN(int_state%CW,IMS,IME,JMS,JME,LM,int_state%INPES)
            CALL POLEHN(int_state%CW,IMS,IME,JMS,JME,LM                 &
                       ,int_state%INPES,int_state%JNPES)
!
            CALL SWAPHN(int_state%O3,IMS,IME,JMS,JME,LM,int_state%INPES)
            CALL POLEHN(int_state%O3,IMS,IME,JMS,JME,LM                 &
                       ,int_state%INPES,int_state%JNPES)
!
            CALL SWAPWN(int_state%U,IMS,IME,JMS,JME,LM,int_state%INPES)
            CALL SWAPWN(int_state%V,IMS,IME,JMS,JME,LM,int_state%INPES)
            CALL POLEWN(int_state%U,int_state%V,IMS,IME,JMS,JME,LM      &
                       ,int_state%INPES,int_state%JNPES)
!
            td%pole_swap_tim=td%pole_swap_tim+(timef()-btim)
          ENDIF
!
!-----------------------------------------------------------------------
!***  Exchange U, V, T, Q and CW
!-----------------------------------------------------------------------
!
          btim=timef()
!
          CALL HALO_EXCH(int_state%T,LM                                 &
                        ,3,3)
!
          CALL HALO_EXCH(int_state%Q,LM,int_state%CW,LM                 &
                        ,3,3)
!
          CALL HALO_EXCH(int_state%O3,LM                                &
                        ,3,3)
!
          CALL HALO_EXCH(int_state%U,LM,int_state%V,LM                  &
                        ,3,3)
!
          td%exch_phy=td%exch_phy+(timef()-btim)
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
#endif
        ENDIF gfs_physics
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
        DEALLOCATE(LONSPERLAR)
        DEALLOCATE(GLOBAL_LATS_R)
        DEALLOCATE(CLDCOV_V)
        DEALLOCATE(PRSL)
        DEALLOCATE(PRSLK)
        DEALLOCATE(GU)
        DEALLOCATE(GV)
        DEALLOCATE(GT)
        DEALLOCATE(GR)
        DEALLOCATE(VVEL)
        DEALLOCATE(F_ICE)
        DEALLOCATE(F_RAIN)
        DEALLOCATE(R_RIME)
        DEALLOCATE(ADT)
        DEALLOCATE(ADU)
        DEALLOCATE(ADV)
        DEALLOCATE(PHIL)
        DEALLOCATE(GR3)
        DEALLOCATE(ADR)
        DEALLOCATE(PRSI)
        DEALLOCATE(PRSIK)
        DEALLOCATE(RSGM)
        DEALLOCATE(PHII)
        DEALLOCATE(SINLAT_R)
        DEALLOCATE(COSLAT_R)
        DEALLOCATE(XLON)
        DEALLOCATE(COSZEN)
        DEALLOCATE(COSZDG)
        DEALLOCATE(RANN)
        DEALLOCATE(RANNUM)
        DEALLOCATE(GR1)
        DEALLOCATE(SWH)
        DEALLOCATE(HLW)
        DEALLOCATE(DKH)
        DEALLOCATE(RNP)
        DEALLOCATE(UPD_MF)
        DEALLOCATE(DWN_MF)
        DEALLOCATE(DET_MF)
        DEALLOCATE(DQDT)
        DEALLOCATE(DQ3DT)
        DEALLOCATE(DT3DT)
        DEALLOCATE(DU3DT)
        DEALLOCATE(DV3DT)
        DEALLOCATE(PHY_F3DV)
!
!#######################################################################
!#######################################################################
!######### E N D   O F   G F S   P H Y S I C S   D R I V E R ###########
!#######################################################################
!#######################################################################
!-----------------------------------------------------------------------
!
#endif
      ENDIF  gfs_phys_test 
!
!-----------------------------------------------------------------------
!***  Write precipitation files for ADJPPT regression test
!-----------------------------------------------------------------------
!
      IF( int_state%WRITE_PREC_ADJ   .AND.                              &
          MOD(XTIME,60.) <= 0.001    .AND.                              &
          INT(XTIME/60.) <= int_state%PCPHR ) THEN
        CALL WRT_PCP(int_state%PREC                                     &
                ,MYPE,NUM_PES,MPI_COMM_COMP,MY_DOMAIN_ID                &
                ,INT(XTIME/60.)+1                                       &
                ,IDS,IDE,JDS,JDE                                        &
                ,IMS,IME,JMS,JME                                        &
                ,ITS,ITE,JTS,JTE)
      ENDIF
!
      ENDIF  physics
!
!-----------------------------------------------------------------------
!***  Run the tracker
!-----------------------------------------------------------------------
!
      IF(int_state%NTRACK>0 .AND. int_state%MYPE<int_state%NUM_PES .and. &
           (int_state%NTSD==0 .or. &
           MOD(int_state%NTSD+1,int_state%NTRACK*int_state%NPHS)==0)) THEN
         CALL QUASIPOST(INT_STATE)
         CALL TRACKER_CENTER(INT_STATE)
      ENDIF
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!---- PHY_RUN END ------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      RC=0
!
      IF(RC_RUN==ESMF_SUCCESS)THEN
!       WRITE(0,*)'SOLVER RUN STEP SUCCEEDED'
      ELSE
        WRITE(0,*)'SOLVER RUN STEP FAILED RC_RUN=',RC_RUN
      ENDIF
!
!-----------------------------------------------------------------------
!
      td%solver_phy_tim=td%solver_phy_tim+(timef()-btim0)

!     write(*,*) 'DEBUG-GT,  ending SOLVER_RUN'
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE SOLVER_RUN
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE SOLVER_FINALIZE (GRID_COMP                             &
                                 ,IMP_STATE                             &
                                 ,EXP_STATE                             &
                                 ,CLOCK_ATM                             &
                                 ,RC_FINALIZE)
!
!-----------------------------------------------------------------------
!***  Finalize the Solver component.
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(ESMF_GridComp) :: GRID_COMP                                     !<-- The Solver gridded component
!
      TYPE(ESMF_State) :: IMP_STATE                                     &  !<-- The Solver import state
                         ,EXP_STATE                                        !<-- The Solver export state
!
      TYPE(ESMF_Clock) :: CLOCK_ATM                                        !<-- The ATM component's ESMF Clock.
!
      INTEGER,INTENT(OUT) :: RC_FINALIZE
!      
!---------------------
!***  Local Variables
!---------------------
!
      TYPE(SOLVER_INTERNAL_STATE),POINTER :: INT_STATE                     !<-- The Solver internal state pointer 
!
      TYPE(WRAP_SOLVER_INT_STATE) :: WRAP                                  !<-- The F90 'wrap' for the Solver internal state
!
      INTEGER(kind=KINT) :: MYPE,RC,RC_FINAL
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC      =ESMF_SUCCESS
      RC_FINAL=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  Extract the Solver internal state.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="SOLVER_FINALIZE: Extract Solver Internal State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompGetInternalState(GRID_COMP                      &  !<-- The Solver component
                                        ,WRAP                           &
                                        ,RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      INT_STATE=>wrap%INT_STATE
!
      MYPE=int_state%MYPE                                                  !<-- The local task rank
!
      IF(MYPE==0)THEN
        WRITE(0,*)' Solver Completed Normally.'
      ENDIF
!
!-----------------------------------------------------------------------
!***  DO NOT DEALLOCATE THE SOLVER INTERNAL STATE POINTER 
!***  WITHOUT DEALLOCATING ITS CONTENTS.
!-----------------------------------------------------------------------
!
!!!   DEALLOCATE(INT_STATE,stat=RC)
!
!-----------------------------------------------------------------------
!
      IF(RC_FINAL==ESMF_SUCCESS)THEN
        WRITE(0,*)'SOLVER FINALIZE STEP SUCCEEDED'
      ELSE
        WRITE(0,*)'SOLVER FINALIZE STEP FAILED'
      ENDIF
!
!     IF(PRESENT(RC_FINALIZE))THEN
        RC_FINALIZE=RC_FINAL
!     ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE SOLVER_FINALIZE
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE BUILD_BC_BUNDLE(GRID                                   &
                                ,LNSH,LNSV                              &
                                ,IHALO,JHALO                            &
                                ,UBOUND_VARS                            &
                                ,VARS                                   &
                                ,MY_DOMAIN_ID                           &
                                ,BUNDLE_NESTBC                          &
                                ,BND_VARS_H                             &
                                ,BND_VARS_V                             &
                                ,NVARS_BC_2D_H                          &
                                ,NVARS_BC_3D_H                          &
                                ,NVARS_BC_4D_H                          &
                                ,NVARS_BC_2D_V                          &
                                ,NVARS_BC_3D_V                          &
                                ,NLEV_H                                 &
                                ,NLEV_V                                 &
                                ,N_BC_3D_H                              &
                                   )
!
!-----------------------------------------------------------------------
!***  This routine builds an ESMF Bundle for holding groups of pointers
!***  to Solver internal state variables that are updated on the 
!***  domain boundaries during the integration.
!***  In addition the object that holds primary boundary information
!***  is partially allocated and pointed at the relevant variables.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      INTEGER(kind=KINT),INTENT(IN) :: IHALO,JHALO                      &  !<-- Subdomain halo widths
                                      ,LNSH,LNSV                           !<-- Domain boundary blending width
      INTEGER(kind=KINT),INTENT(IN) :: MY_DOMAIN_ID                     &  !<-- This domain's ID
                                      ,UBOUND_VARS                         !<-- Upper dimension of the VARS array
!
      INTEGER(kind=KINT),DIMENSION(1:3),INTENT(OUT) :: N_BC_3D_H           !<-- Hold order of domain #1's BC vbls from boco files
!
      TYPE(ESMF_Grid),INTENT(IN) :: GRID                                   !<-- The ESMF Grid for this domain
!
      TYPE(VAR),DIMENSION(1:UBOUND_VARS),INTENT(INOUT) :: VARS             !<-- Variables in the Solver internal state
!
      TYPE(ESMF_FieldBundle),INTENT(INOUT) :: BUNDLE_NESTBC                !<-- The Bundle of Solver internal state vbls to be
!                                                                          !    updated on the nest boundaries
      INTEGER(kind=KINT),INTENT(OUT) :: NVARS_BC_2D_H                   &  !<-- # of 2-D,3-D,4-D H-pt variables
                                       ,NVARS_BC_3D_H                   &  !    that are inserted
                                       ,NVARS_BC_4D_H                      !    into the Bundle.
!
      INTEGER(kind=KINT),INTENT(OUT) :: NVARS_BC_2D_V                   &  !<-- # of 2-D,3-D V-pt variables
                                       ,NVARS_BC_3D_V                      !    that are inserted into the Bundle.
!
      INTEGER(kind=KINT),INTENT(OUT) :: NLEV_H,NLEV_V                      !<-- # of model levels in all H-pt,V-pt variables used
!
      TYPE(BC_H_ALL),INTENT(OUT) :: BND_VARS_H                             !<-- Object holding H-pt variable info on domain boundaries
      TYPE(BC_V_ALL),INTENT(OUT) :: BND_VARS_V                             !<-- Object holding V-pt variable info on domain boundaries
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: H_OR_V_INT,IOS,LB3,LB4                      &
                           ,N,NSIZE,NUM_DIMS,NUM_FIELDS                 &
                           ,UB3,UB4
!
      INTEGER(kind=KINT) :: IMS,IME,JMS,JME
!
      INTEGER(kind=KINT) :: KNT_2D_H,KNT_3D_H,KNT_4D_H                  &
                           ,KNT_2D_V,KNT_3D_V
!
      INTEGER(kind=KINT) :: KNT_3D_DOM_01
!
      INTEGER(kind=KINT) :: ISTAT,RC,RC_CMB
!
      REAL(kind=KFPT),DIMENSION(:,:),POINTER :: ARRAY_2D=>NULL()
      REAL(kind=KFPT),DIMENSION(:,:,:),POINTER :: ARRAY_3D=>NULL()
      REAL(kind=KFPT),DIMENSION(:,:,:,:),POINTER :: ARRAY_4D=>NULL()
!
      CHARACTER(len=1) :: CH_2,CH_B,H_OR_V
!           
      CHARACTER(len=2) :: CH_M
!           
      CHARACTER(len=9),SAVE :: FNAME='nests.txt'
!
      CHARACTER(len=99) :: BUNDLE_NAME,FIELD_NAME,VBL_NAME
!
      CHARACTER(len=256) :: STRING
!
      LOGICAL(kind=KLOG) :: CASE_2WAY,CASE_NESTBC
!
      TYPE(ESMF_Field) :: FIELD_X
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      IMS=int_state%IMS
      IME=int_state%IME
      JMS=int_state%JMS
      JME=int_state%JME
!
      NVARS_BC_2D_H=0     
      NVARS_BC_3D_H=0     
      NVARS_BC_4D_H=0     
!
      NVARS_BC_2D_V=0     
      NVARS_BC_3D_V=0     
!
      KNT_3D_DOM_01=0
!
      DO N=1,3
        N_BC_3D_H(N)=-1
      ENDDO
!
!-----------------------------------------------------------------------
!***  Loop through all Solver internal state variables.
!-----------------------------------------------------------------------
!
      OPEN(unit=10,file=FNAME,status='OLD',action='READ'                &  !<-- Open the text file with user specifications
            ,iostat=IOS)
!
      IF(IOS/=0)THEN
        WRITE(0,*)' Failed to open ',FNAME,' so ABORT!'
        CALL ESMF_FINALIZE(endflag=ESMF_END_ABORT                   &
                          ,rc             =RC)
      ENDIF
!
      NLEV_H=0                                                               !<-- Counter for total # of levels in all 2-way vbls
      NLEV_V=0                                                               !<-- Counter for total # of levels in all 2-way vbls
!
!-----------------------------------------------------------------------
      bundle_loop: DO
!-----------------------------------------------------------------------
!
        READ(UNIT=10,FMT='(A)',iostat=IOS)STRING                           !<-- Read in the next specification line
        IF(IOS/=0)THEN                                                     !<-- Finished reading the specification lines
          CLOSE(10)
          EXIT
        ENDIF
!
        IF(STRING(1:1)=='#'.OR.TRIM(STRING)=='')THEN
          CYCLE                                                            !<-- Read past comments and blanks.
        ENDIF
!
!-----------------------------------------------------------------------
!***  Read the text line containing the H or V specification for 
!***  variable N then find that variable's place within the VARS 
!***  object.
!-----------------------------------------------------------------------
!
        READ(UNIT=STRING,FMT=*,iostat=IOS)VBL_NAME                      &  !<-- The variable's name in the text file.
                                         ,CH_B                          &  !<-- The flag for nest BC vbls in the text file.
                                         ,CH_M                          &  !<-- Not relevant here (flag for moving nests)
                                         ,CH_2                             !<-- The flag for 2-way vbls in the text file.
!
        CALL FIND_VAR_INDX(VBL_NAME,VARS,UBOUND_VARS,N)
!
        FIELD_NAME=TRIM(VARS(N)%VBL_NAME)//TRIM(SUFFIX_NESTBC)             !<-- Append the BC suffix to the Field
!                                                                 
!-----------------------------------------------------------------------
!***  Check the Bundle's name to determine which column of user
!***  specifications to read from the text file.
!-----------------------------------------------------------------------
!
        H_OR_V=CH_B                                                        !<-- H-V flag for this nest BC variable
!
!-----------------------------------------------------------------------
!***  Find the variables in the Solver internal state that have been
!***  selected to be placed into the Bundle.  The user has specified
!***  whether the variable lies on H points or V points.
!***  Currently ESMF will not allow the use of Attributes that are
!***  characters therefore we must translate the character codes from
!***  the txt file into something that ESMF can use.  In this case
!***  we will use integers:  H-->1 and V-->2 .
!-----------------------------------------------------------------------
!
        IF(H_OR_V=='H')THEN
          H_OR_V_INT=1                                                     !<-- H-pt variable
        ELSEIF(H_OR_V=='V')THEN
          H_OR_V_INT=2                                                     !<-- V-pt variable
        ELSE
          H_OR_V_INT=-999                                                  !<-- Variable not specified for use.
        ENDIF
!
!-----------------------------------------------------------------------
!
        build_bundle: IF(H_OR_V=='H'                                    &
                            .OR.                                        &
                         H_OR_V=='V'                                    &
                                     )THEN
!
!-----------------------------------------------------------------------
!
!-------------------
!***  2-D Variables
!-------------------
!
!-------------
!***  Integer
!-------------
!
          IF(ASSOCIATED(VARS(N)%I2D))THEN                                  !<-- 2-D integer array on mass points
!
!           FIELD_X=ESMF_FieldCreate(grid       =GRID                   &  !<-- The ESMF Grid for this domain
!                                   ,farray     =VARS(N)%I2D            &  !<-- Nth variable in the VARS array
!                                   ,totalUWidth=(/IHALO,JHALO/)        &  !<-- Upper bound of halo region
!                                   ,totalLWidth=(/IHALO,JHALO/)        &  !<-- Lower bound of halo region
!                                   ,name       =FIELD_NAME             &  !<-- The name of this variable
!                                   ,indexFlag  =ESMF_INDEX_GLOBAL      &  !<-- The variable uses global indexing
!                                   ,rc         =RC)
            WRITE(0,*)' MUST ADD THE CAPABILITY TO USE 2-D INTEGERS IN 2WAY/BC UPDATES!!'
            WRITE(0,*)' Variable name is ',VARS(N)%VBL_NAME,' for variable #',N
            WRITE(0,*)' H_OR_V_INT=',H_OR_V_INT
            WRITE(0,*)' ABORT!!'
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
!
!
!----------
!***  Real
!----------
!
          ELSEIF(ASSOCIATED(VARS(N)%R2D))THEN                              !<-- 2-D real array on mass points
!
            FIELD_X=ESMF_FieldCreate(grid       =GRID                   &  !<-- The ESMF Grid for this domain
                                    ,farray     =VARS(N)%R2D            &  !<-- Nth variable in the VARS array
                                    ,totalUWidth=(/IHALO,JHALO/)        &  !<-- Upper bound of halo region
                                    ,totalLWidth=(/IHALO,JHALO/)        &  !<-- Lower bound of halo region
                                    ,name       =FIELD_NAME             &  !<-- The name of this variable
                                    ,indexFlag  =ESMF_INDEX_GLOBAL      &  !<-- The variable uses global indexing
                                    ,rc         =RC)
!
            IF(H_OR_V=='H')THEN
              NVARS_BC_2D_H=NVARS_BC_2D_H+1                                !<-- Count # of 2-D H-pt variables
              NLEV_H=NLEV_H+1                                              !<-- Sum all levels for H-pt variables
            ELSEIF(H_OR_V=='V')THEN
              NVARS_BC_2D_V=NVARS_BC_2D_V+1                                !<-- Count # of 2-D V-pt variables
              NLEV_V=NLEV_V+1                                              !<-- Sum all levels for V-pt variables
            ENDIF
!
!-------------------
!***  3-D Variables
!-------------------
!
!----------
!***  Real
!----------
!
          ELSEIF(ASSOCIATED(VARS(N)%R3D))THEN                              !<-- 3-D real array on mass points
!
            FIELD_X=ESMF_FieldCreate(grid           =GRID                           &  !<-- The ESMF Grid for this domain
                                    ,farray         =VARS(N)%R3D                    &  !<-- Nth variable in the VARS array
                                    ,totalUWidth    =(/IHALO,JHALO/)                &  !<-- Upper bound of halo region
                                    ,totalLWidth    =(/IHALO,JHALO/)                &  !<-- Lower bound of halo region
                                    ,ungriddedLBound=(/lbound(VARS(N)%R3D,dim=3)/)  &
                                    ,ungriddedUBound=(/ubound(VARS(N)%R3D,dim=3)/)  &
                                    ,name           =FIELD_NAME                     &  !<-- The name of this variable
                                    ,indexFlag      =ESMF_INDEX_GLOBAL              &  !<-- The variable uses global indexing
                                    ,rc             =RC)
!
            LB3=LBOUND(VARS(N)%R3D,3)
            UB3=UBOUND(VARS(N)%R3D,3)
!
            IF(H_OR_V=='H')THEN
              NVARS_BC_3D_H=NVARS_BC_3D_H+1                                !<-- Count # of 3-D H-pt variables
              NLEV_H=NLEV_H+(UB3-LB3+1)                                    !<-- Sum all levels for H-pt variables
            ELSEIF(H_OR_V=='V')THEN
              NVARS_BC_3D_V=NVARS_BC_3D_V+1                                !<-- Count # of 3-D V-pt variables
              NLEV_V=NLEV_V+(UB3-LB3+1)                                    !<-- Sum all levels for V-pt variables
            ENDIF
!
!-------------------
!***  4-D Variables
!-------------------
!
!----------
!***  Real
!----------
!
          ELSEIF(ASSOCIATED(VARS(N)%R4D))THEN                              !<-- 4-D real array on mass points
!
            LB4=LBOUND(VARS(N)%R4D,dim=4)
            UB4=UBOUND(VARS(N)%R4D,dim=4)
!
            FIELD_X=ESMF_FieldCreate(grid           =GRID                           &  !<-- The ESMF Grid for this domain
                                    ,farray         =VARS(N)%R4D                    &  !<-- Nth variable in the VARS array
                                    ,totalUWidth    =(/IHALO,JHALO/)                &  !<-- Upper bound of halo region
                                    ,totalLWidth    =(/IHALO,JHALO/)                &  !<-- Lower bound of halo region
                                    ,ungriddedLBound =(/ LBOUND(VARS(N)%R4D,dim=3),LB4 /) &
                                    ,ungriddedUBound =(/ UBOUND(VARS(N)%R4D,dim=3),UB4 /) &
                                    ,name           =FIELD_NAME                     &  !<-- The name of this variable
                                    ,indexFlag      =ESMF_INDEX_GLOBAL              &  !<-- The variable uses global indexing
                                    ,rc             =RC)
!
            LB3=LBOUND(VARS(N)%R4D,3)
            UB3=UBOUND(VARS(N)%R4D,3)
!
            IF(H_OR_V=='H')THEN
              NVARS_BC_4D_H=NVARS_BC_4D_H+1                                !<-- Count # of 4-D H-pt variables
              NLEV_H=NLEV_H+(UB3-LB3+1)*(UB4-LB4+1)                        !<-- Sum all levels for H-pt variables
            ENDIF
!
!----------------
!***  All Others
!----------------
!
          ELSE
            WRITE(0,*)' SELECTED UPDATE H VARIABLE IS NOT 2,3,4-D REAL.'
            WRITE(0,*)' Variable name is ',VARS(N)%VBL_NAME,' for variable #',N
            WRITE(0,*)' H_OR_V_INT=',H_OR_V_INT
            WRITE(0,*)' ABORT!!'
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
!
          ENDIF
!
!-----------------------------------------------------------------------
!***  Attach the index of this variable within the Solver internal
!***  state so it can be referenced w/r to the boundary objects.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Add Solver Int State Indx to Bundle Field"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeSet(field=FIELD_X                          &  !<-- The Field to be added to the Bundle
                                ,name ='Solver Int State Indx'          &  !<-- The name of the Attribute to set
                                ,value=N                                &  !<-- The index of the Solver internal state vbl
                                ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CMB)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Attach the specification flag to this Field that indicates
!***  whether it is an H-pt or a V-pt variable.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Add H-or-V Specification Flag to Bundle Field"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeSet(field=FIELD_X                          &  !<-- The Field to be added to the Bundle
                                ,name ='H_OR_V_INT'                     &  !<-- The name of the Attribute to set
                                ,value=H_OR_V_INT                       &  !<-- H-pt or V-pt flag 
                                ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CMB)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Add this Field to the Bundle that holds pointers to all
!***  variables in the Solver internal state that have been
!***  selected.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Add Desired Field to the Bundle"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_FieldBundleAdd(            BUNDLE_NESTBC            &  !<-- The Bundle of Solver internal state BC variables
                                  ,            (/FIELD_X/)     &  !<-- Add this Field to the Bundle
                                  ,rc         =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CMB)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
        ENDIF build_bundle
!
!-----------------------------------------------------------------------
!
      ENDDO bundle_loop
!
!-----------------------------------------------------------------------
!***  Allocate the appropriate pieces of the boundary variable
!***  objects.  All nests use the same set of boundary variables
!***  that are specified by the user in the external text file.
!
!***  The upper parent domain uses its own set of boundary variables
!***  updated from the BC files generated during preprocessing.
!***  They are currently hardwired to PD,T,Q,CW,U,V.
!-----------------------------------------------------------------------
!
      IF(MY_DOMAIN_ID==1)THEN                                              !<-- The uppermost parent will hardwire its BC vbls
!
        NVARS_BC_2D_H=1                                                    !<-- PD
        NVARS_BC_3D_H=3                                                    !<-- T,Q,CW
        NVARS_BC_4D_H=0
        NVARS_BC_2D_V=0
        NVARS_BC_3D_V=2                                                    !<-- U,V
!
      ENDIF
!
      IF(NVARS_BC_2D_H>0)THEN
        ALLOCATE(BND_VARS_H%VAR_2D(1:NVARS_BC_2D_H))                       !<-- All 2-D H-pt nest boundary variables
      ENDIF
!
      IF(NVARS_BC_3D_H>0)THEN
        ALLOCATE(BND_VARS_H%VAR_3D(1:NVARS_BC_3D_H))                       !<-- All 3-D H-pt nest boundary variables
      ENDIF
!
      IF(NVARS_BC_4D_H>0)THEN
        ALLOCATE(BND_VARS_H%VAR_4D(1:NVARS_BC_4D_H))                       !<-- All 4-D H-pt nest boundary variables
      ENDIF
!
      IF(NVARS_BC_2D_V>0)THEN
        ALLOCATE(BND_VARS_V%VAR_2D(1:NVARS_BC_2D_V))                       !<-- All 2-D V-pt nest boundary variables
      ENDIF
!
      IF(NVARS_BC_3D_V>0)THEN
        ALLOCATE(BND_VARS_V%VAR_3D(1:NVARS_BC_3D_V))                       !<-- All 3-D V-pt nest boundary variables
      ENDIF
!
!-----------------------------------------------------------------------
!***  Now go through the boundary Bundle's variables and point the
!***  full variable pointer of the appropriate boundary object that
!***  was allocated immediately above at that variable in the Bundle
!***  in order to associate each piece of the boundary object with
!***  the actual boundary variable.
!-----------------------------------------------------------------------
!
      CALL ESMF_FieldBundleGet(FIELDBUNDLE=BUNDLE_NESTBC                &  !<-- Bundle holding the arrays for BC updates
                              ,fieldCount =NUM_FIELDS                   &  !<-- Number of Fields in the Bundle
                              ,rc         =RC )
!
!-----------------------------------------------------------------------
!
      KNT_2D_H=0
      KNT_3D_H=0
      KNT_4D_H=0
      KNT_2D_V=0
      KNT_3D_V=0
!
!-----------------------------------------------------------------------
!
      bc_fields: DO N=1,NUM_FIELDS
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract Field N from the Bundle"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_FieldBundleGet(FIELDBUNDLE=BUNDLE_NESTBC              &  !<-- Bundle holding the arrays for BC updates
                                ,fieldIndex =N                          &  !<-- Index of the Field in the Bundle
                                ,field      =FIELD_X                    &  !<-- Field N in the Bundle
                                ,rc         =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CMB)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract H_OR_V Flag from the Field"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(field=FIELD_X                            &  !<-- The Domain import state
                              ,name ='H_OR_V_INT'                       &  !<-- Name of the Attribute
                              ,value=H_OR_V_INT                         &  !<-- Is the Field on H or V points? (1 is H; 2 is V)
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CMB)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------
!***  H-pt boundary variables
!-----------------------------
!
        h_v: IF(H_OR_V_INT==1)THEN                                         !<-- If true, it is an H-pt variable
!
          CALL ESMF_FieldGet(field   =FIELD_X                           &  !<-- Field N in the Bundle
                            ,dimCount=NUM_DIMS                          &  !<-- How many dimensions?
                            ,rc      =RC )
!
!--------------
!***  2-D Real
!--------------
!
          IF(NUM_DIMS==2)THEN
!
            KNT_2D_H=KNT_2D_H+1
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Extract the 2-D H-pt Array from the Field"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_FieldGet(field    =FIELD_X                        &  !<-- Field N in the Bundle
                              ,localDe  =0                              &
                              ,farrayPtr=ARRAY_2D                       &  !<-- Dummy 2-D array with Field's Real data
                              ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CMB)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            BND_VARS_H%VAR_2D(KNT_2D_H)%FULL_VAR=>ARRAY_2D                 !<-- This variable becomes a boundary variable
!
            ALLOCATE(BND_VARS_H%VAR_2D(KNT_2D_H)%SOUTH(IMS:IME,1:LNSH,1:2) & !<-- 2-D H-pt boundary variable N on domain's south side
                    ,stat=ISTAT)
            IF(ISTAT/=0)THEN
              WRITE(0,45001)ISTAT
45001         FORMAT(' Failed to allocate BND_VARS_H%VAR_2D(N)%SOUTH  istat=',i5)
              WRITE(0,*)' Aborting!'
              CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
            ENDIF
!
            ALLOCATE(BND_VARS_H%VAR_2D(KNT_2D_H)%NORTH(IMS:IME,1:LNSH,1:2) & !<-- 2-D H-pt boundary variable N on domain's north side
                    ,stat=ISTAT)
            IF(ISTAT/=0)THEN
              WRITE(0,45002)ISTAT
45002         FORMAT(' Failed to allocate BND_VARS_H%VAR_2D(N)%NORTH  istat=',i5)
              WRITE(0,*)' Aborting!'
              CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
            ENDIF
!
            ALLOCATE(BND_VARS_H%VAR_2D(KNT_2D_H)%WEST(1:LNSH,JMS:JME,1:2)  & !<-- 2-D H-pt boundary variable N on domain's west side
                    ,stat=ISTAT)
            IF(ISTAT/=0)THEN
              WRITE(0,45003)ISTAT
45003         FORMAT(' Failed to allocate BND_VARS_H%VAR_2D(N)%WEST  istat=',i5)
              WRITE(0,*)' Aborting!'
              CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
            ENDIF
!
            ALLOCATE(BND_VARS_H%VAR_2D(KNT_2D_H)%EAST(1:LNSH,JMS:JME,1:2)  & !<-- 2-D H-pt boundary variable N on domain's east side
                    ,stat=ISTAT)
            IF(ISTAT/=0)THEN
              WRITE(0,45004)ISTAT
45004         FORMAT(' Failed to allocate BND_VARS_H%VAR_2D(N)%EasT  istat=',i5)
              WRITE(0,*)' Aborting!'
              CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
            ENDIF
!
            BND_VARS_H%VAR_2D(KNT_2D_H)%SOUTH=R4_IN
            BND_VARS_H%VAR_2D(KNT_2D_H)%NORTH=R4_IN
            BND_VARS_H%VAR_2D(KNT_2D_H)%WEST=R4_IN
            BND_VARS_H%VAR_2D(KNT_2D_H)%EAST=R4_IN
!
!--------------
!***  3-D Real
!--------------
!
          ELSEIF(NUM_DIMS==3)THEN
!
            KNT_3D_H=KNT_3D_H+1
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Extract the 3-D H-pt Array from the Field"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_FieldGet(field    =FIELD_X                        &  !<-- Field N in the Bundle
                              ,localDe  =0                              &
                              ,farrayPtr=ARRAY_3D                       &  !<-- Dummy 3-D array with Field's Real data
                              ,rc       =RC )
!
            CALL ESMF_FieldGet(field=FIELD_X                            &  !<-- Field N in the Bundle
                              ,name =FIELD_NAME                         &  !<-- This Field's name
                              ,rc   =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CMB)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            BND_VARS_H%VAR_3D(KNT_3D_H)%FULL_VAR=>ARRAY_3D                 !<-- This variable becomes a boundary variable
!
            LB3=LBOUND(ARRAY_3D,3)
            UB3=UBOUND(ARRAY_3D,3)
!
            ALLOCATE(BND_VARS_H%VAR_3D(KNT_3D_H)%SOUTH(IMS:IME,1:LNSH,LB3:UB3,1:2) & !<-- 3-D H-pt bndry vbl N on domain's south side
                    ,stat=ISTAT)
            IF(ISTAT/=0)THEN
              WRITE(0,45011)ISTAT
45011         FORMAT(' Failed to allocate BND_VARS_H%VAR_3D(N)%SOUTH  istat=',i5)
              WRITE(0,*)' Aborting!'
              CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
            ENDIF
!
            ALLOCATE(BND_VARS_H%VAR_3D(KNT_3D_H)%NORTH(IMS:IME,1:LNSH,LB3:UB3,1:2) & !<-- 3-D H-pt bndry vbl N on domain's north side
                    ,stat=ISTAT)
            IF(ISTAT/=0)THEN
              WRITE(0,45012)ISTAT
45012         FORMAT(' Failed to allocate BND_VARS_H%VAR_3D(N)%NORTH  istat=',i5)
              WRITE(0,*)' Aborting!'
              CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
            ENDIF
!
            ALLOCATE(BND_VARS_H%VAR_3D(KNT_3D_H)%WEST(1:LNSH,JMS:JME,LB3:UB3,1:2)  & !<-- 3-D H-pt bndry vbl N on domain's west side
                    ,stat=ISTAT)
            IF(ISTAT/=0)THEN
              WRITE(0,45013)ISTAT
45013         FORMAT(' Failed to allocate BND_VARS_H%VAR_3D(N)%WEST  istat=',i5)
              WRITE(0,*)' Aborting!'
              CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
            ENDIF
!
            ALLOCATE(BND_VARS_H%VAR_3D(KNT_3D_H)%EAST(1:LNSH,JMS:JME,LB3:UB3,1:2)  & !<-- 3-D H-pt bndry vbl N on domain's east side
                    ,stat=ISTAT)
            IF(ISTAT/=0)THEN
              WRITE(0,45014)ISTAT
45014         FORMAT(' Failed to allocate BND_VARS_H%VAR_3D(N)%EAST  istat=',i5)
              WRITE(0,*)' Aborting!'
              CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
            ENDIF
!
            BND_VARS_H%VAR_3D(KNT_3D_H)%SOUTH=R4_IN
            BND_VARS_H%VAR_3D(KNT_3D_H)%NORTH=R4_IN
            BND_VARS_H%VAR_3D(KNT_3D_H)%WEST=R4_IN
            BND_VARS_H%VAR_3D(KNT_3D_H)%EAST=R4_IN
!
!-----------------------------------------------------------------------
!***  Now some hardwiring is required.  The same boundary objects
!***  are of course used by all domains but the arrays for the
!***  upper domain will be read in from the external boco files
!***  when that domain is not global.  The boundary objects store
!***  the arrays in the order they are encountered in the nests.txt
!***  file and in general that order will be different than the 
!***  order they are read from the boco files.  Therefore we now
!***  save the order of the three 3-D H-pt boundary arrays used
!***  by the upper parent so they can be saved in the proper order
!***  when they are read in subroutine READ_BC.  The order in which
!***  READ_BC reads them from the boco files is T,Q,CW.
!***  REGRETTABLY THIS IS DIRTY but is needed since all domains 
!***  must use the same boundary objects but the dataread in
!***  READ_BC is fixed in its order in NPS.
!-----------------------------------------------------------------------
!
            IF(FIELD_NAME(1:1)=='T'                                     &
                         .OR.                                           &
               FIELD_NAME(1:1)=='Q'                                     &
                         .OR.                                           &
               FIELD_NAME(1:2)=='CW')THEN
!
              KNT_3D_DOM_01=KNT_3D_DOM_01+1
!
              IF(FIELD_NAME(1:1)=='T')THEN
                N_BC_3D_H(1)=KNT_3D_DOM_01
              ELSEIF(FIELD_NAME(1:1)=='Q')THEN
                N_BC_3D_H(2)=KNT_3D_DOM_01
              ELSEIF(FIELD_NAME(1:2)=='CW')THEN
                N_BC_3D_H(3)=KNT_3D_DOM_01
              ENDIF
!
            ENDIF
!
!--------------
!***  4-D Real
!--------------
!
          ELSEIF(NUM_DIMS==4)THEN
!
            KNT_4D_H=KNT_4D_H+1
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Extract the 4-D H-pt Array from the Field"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_FieldGet(field    =FIELD_X                        &  !<-- Field N in the Bundle
                              ,localDe  =0                              &
                              ,farrayPtr=ARRAY_4D                       &  !<-- Dummy 4-D array with Field's Real data
                              ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CMB)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            BND_VARS_H%VAR_4D(KNT_4D_H)%FULL_VAR=>ARRAY_4D                 !<-- This variable becomes a boundary variable
!
            LB3=LBOUND(ARRAY_4D,3)
            UB3=UBOUND(ARRAY_4D,3)
            LB4=LBOUND(ARRAY_4D,4)
            UB4=UBOUND(ARRAY_4D,4)
!
            ALLOCATE(BND_VARS_H%VAR_4D(KNT_4D_H)%SOUTH(IMS:IME,1:LNSH,LB4:UB4,1:2,LB4:UB4) & !<-- 4-D H-pt bndry vbl N on domain's south side
                    ,stat=ISTAT)
            IF(ISTAT/=0)THEN
              WRITE(0,45021)ISTAT
45021         FORMAT(' Failed to allocate BND_VARS_H%VAR_4D(N)%SOUTH  istat=',i5)
              WRITE(0,*)' Aborting!'
              CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
            ENDIF
!
            ALLOCATE(BND_VARS_H%VAR_4D(KNT_4D_H)%NORTH(IMS:IME,1:LNSH,LB3:UB3,1:2,LB4:UB4) & !<-- 4-D H-pt bndry vbl N on domain's north side
                    ,stat=ISTAT)
            IF(ISTAT/=0)THEN
              WRITE(0,45022)ISTAT
45022         FORMAT(' Failed to allocate BND_VARS_H%VAR_4D(N)%NORTH  istat=',i5)
              WRITE(0,*)' Aborting!'
              CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
            ENDIF
!
            ALLOCATE(BND_VARS_H%VAR_4D(KNT_4D_H)%WEST(1:LNSH,JMS:JME,LB3:UB3,1:2,LB4:UB4)  & !<-- 4-D H-pt bndry vbl N on domain's west side
                    ,stat=ISTAT)
            IF(ISTAT/=0)THEN
              WRITE(0,45023)ISTAT
45023         FORMAT(' Failed to allocate BND_VARS_H%VAR_4D(N)%WEST  istat=',i5)
              WRITE(0,*)' Aborting!'
              CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
            ENDIF
!
            ALLOCATE(BND_VARS_H%VAR_4D(KNT_4D_H)%EAST(1:LNSH,JMS:JME,LB3:UB3,1:2,LB4:UB4)  & !<-- 4-D H-pt bndry vbl N on domain's east side
                    ,stat=ISTAT)
            IF(ISTAT/=0)THEN
              WRITE(0,45024)ISTAT
45024         FORMAT(' Failed to allocate BND_VARS_H%VAR_4D(N)%EAST  istat=',i5)
              WRITE(0,*)' Aborting!'
              CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
            ENDIF
!
            BND_VARS_H%VAR_4D(KNT_4D_H)%SOUTH=R4_IN
            BND_VARS_H%VAR_4D(KNT_4D_H)%NORTH=R4_IN
            BND_VARS_H%VAR_4D(KNT_4D_H)%WEST=R4_IN
            BND_VARS_H%VAR_4D(KNT_4D_H)%EAST=R4_IN
!
          ENDIF
!
!-----------------------------
!***  V-pt boundary variables
!-----------------------------
!
        ELSEIF(H_OR_V_INT==2)THEN
!
          CALL ESMF_FieldGet(field   =FIELD_X                           &  !<-- Field N in the Bundle
                            ,dimCount=NUM_DIMS                          &  !<-- How many dimensions?
                            ,rc      =RC )
!
!--------------
!***  2-D Real
!--------------
!
          IF(NUM_DIMS==2)THEN
!
            KNT_2D_V=KNT_2D_V+1
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Extract the 2-D V-pt Array from the Field"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_FieldGet(field    =FIELD_X                        &  !<-- Field N in the Bundle
                              ,localDe  =0                              &
                              ,farrayPtr=ARRAY_2D                       &  !<-- Dummy 2-D array with Field's Real data
                              ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CMB)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            BND_VARS_V%VAR_2D(KNT_2D_V)%FULL_VAR=>ARRAY_2D                   !<-- This variable becomes a boundary variable
!
            ALLOCATE(BND_VARS_V%VAR_2D(KNT_2D_V)%SOUTH(IMS:IME,1:LNSV,1:2) & !<-- 2-D V-pt bndry vbl N on domain's south side
                    ,stat=ISTAT)
            IF(ISTAT/=0)THEN
              WRITE(0,45041)ISTAT
45041         FORMAT(' Failed to allocate BND_VARS_V%VAR_2D(N)%SOUTH  istat=',i5)
              WRITE(0,*)' Aborting!'
              CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
            ENDIF
!
            ALLOCATE(BND_VARS_V%VAR_2D(KNT_2D_V)%NORTH(IMS:IME,1:LNSV,1:2) & !<-- 2-D V-pt bndry vbl N on domain's north side
                    ,stat=ISTAT)
            IF(ISTAT/=0)THEN
              WRITE(0,45042)ISTAT
45042         FORMAT(' Failed to allocate BND_VARS_V%VAR_2D(N)%NORTH  istat=',i5)
              WRITE(0,*)' Aborting!'
              CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
            ENDIF
!
            ALLOCATE(BND_VARS_V%VAR_2D(KNT_2D_V)%WEST(1:LNSV,JMS:JME,1:2)  & !<-- 2-D V-pt bndry vbl N on domain's west side
                    ,stat=ISTAT)
            IF(ISTAT/=0)THEN
              WRITE(0,45043)ISTAT
45043         FORMAT(' Failed to allocate BND_VARS_V%VAR_2D(N)%WEST  istat=',i5)
              WRITE(0,*)' Aborting!'
              CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
            ENDIF
!
            ALLOCATE(BND_VARS_V%VAR_2D(KNT_2D_V)%EAST(1:LNSV,JMS:JME,1:2)  & !<-- 2-D V-pt bndry vbl N on domain's east side
                    ,stat=ISTAT)
            IF(ISTAT/=0)THEN
              WRITE(0,45044)ISTAT
45044         FORMAT(' Failed to allocate BND_VARS_V%VAR_2D(N)%EAST  istat=',i5)
              WRITE(0,*)' Aborting!'
              CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
            ENDIF
!
            BND_VARS_V%VAR_2D(KNT_2D_V)%SOUTH=R4_IN
            BND_VARS_V%VAR_2D(KNT_2D_V)%NORTH=R4_IN
            BND_VARS_V%VAR_2D(KNT_2D_V)%WEST=R4_IN
            BND_VARS_V%VAR_2D(KNT_2D_V)%EAST=R4_IN
!
!--------------
!***  3-D Real
!--------------
!
          ELSEIF(NUM_DIMS==3)THEN
!
            KNT_3D_V=KNT_3D_V+1
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Extract the 3-D V-pt Array from the Field"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_FieldGet(field    =FIELD_X                        &  !<-- Field N in the Bundle
                              ,localDe  =0                              &
                              ,farrayPtr=ARRAY_3D                       &  !<-- Dummy 3-D array with Field's Real data
                              ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CMB)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            BND_VARS_V%VAR_3D(KNT_3D_V)%FULL_VAR=>ARRAY_3D                 !<-- This variable becomes a boundary variable
!
            LB3=LBOUND(ARRAY_3D,3)
            UB3=UBOUND(ARRAY_3D,3)
!
            ALLOCATE(BND_VARS_V%VAR_3D(KNT_3D_V)%SOUTH(IMS:IME,1:LNSV,LB3:UB3,1:2) & !<-- 3-D V-pt bndry vbl N on domain's south side
                    ,stat=ISTAT)
            IF(ISTAT/=0)THEN
              WRITE(0,45051)ISTAT
45051         FORMAT(' Failed to allocate BND_VARS_V%VAR_3D(N)%SOUTH  istat=',i5)
              WRITE(0,*)' Aborting!'
              CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
            ENDIF
!
            ALLOCATE(BND_VARS_V%VAR_3D(KNT_3D_V)%NORTH(IMS:IME,1:LNSV,LB3:UB3,1:2) & !<-- 3-D V-pt bndry vbl N on domain's north side
                    ,stat=ISTAT)
            IF(ISTAT/=0)THEN
              WRITE(0,45052)ISTAT
45052         FORMAT(' Failed to allocate BND_VARS_V%VAR_3D(N)%NORTH  istat=',i5)
              WRITE(0,*)' Aborting!'
              CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
            ENDIF
!
            ALLOCATE(BND_VARS_V%VAR_3D(KNT_3D_V)%WEST(1:LNSV,JMS:JME,LB3:UB3,1:2)  & !<-- 3-D V-pt bndry vbl N on domain's west side
                    ,stat=ISTAT)
            IF(ISTAT/=0)THEN
              WRITE(0,45053)ISTAT
45053         FORMAT(' Failed to allocate BND_VARS_V%VAR_3D(N)%WEST  istat=',i5)
              WRITE(0,*)' Aborting!'
              CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
            ENDIF
!
            ALLOCATE(BND_VARS_V%VAR_3D(KNT_3D_V)%EAST(1:LNSV,JMS:JME,LB3:UB3,1:2)  & !<-- 3-D V-pt bndry vbl N on domain's east side
                    ,stat=ISTAT)
            IF(ISTAT/=0)THEN
              WRITE(0,45054)ISTAT
45054         FORMAT(' Failed to allocate BND_VARS_V%VAR_3D(N)%EAST  istat=',i5)
              WRITE(0,*)' Aborting!'
              CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
            ENDIF
!
            BND_VARS_V%VAR_3D(KNT_3D_V)%SOUTH=R4_IN
            BND_VARS_V%VAR_3D(KNT_3D_V)%NORTH=R4_IN
            BND_VARS_V%VAR_3D(KNT_3D_V)%WEST=R4_IN
            BND_VARS_V%VAR_3D(KNT_3D_V)%EAST=R4_IN
!
          ENDIF
!
!-----------------------------------------------------------------------
!
        ENDIF h_v
!
!-----------------------------------------------------------------------
!
      ENDDO bc_fields
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE BUILD_BC_BUNDLE
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE UPDATE_BC_TENDS(IMP_STATE                              &
                                ,LM,LNSH,LNSV                           &
                                ,PARENT_CHILD_TIME_RATIO,DT             &
                                ,S_BDY,N_BDY,W_BDY,E_BDY                &
                                ,NLEV_H,NLEV_V                          &
                                ,NVARS_BC_2D_H                          &
                                ,NVARS_BC_3D_H                          &
                                ,NVARS_BC_4D_H                          &
                                ,NVARS_BC_2D_V                          &
                                ,NVARS_BC_3D_V                          &
                                ,BND_VARS_H                             &
                                ,BND_VARS_V                             &
                                ,ITS,ITE,JTS,JTE                        &
                                ,IMS,IME,JMS,JME                        &
                                ,IDS,IDE,JDS,JDE )
! 
!-----------------------------------------------------------------------
!***  This routine extracts boundary data from the Solver import
!***  state of nested domains that was received from their parents.
!***  This data is then used to update the time tendencies of the
!***  boundary variables.  Those tendencies are valid through each
!***  timestep of the nested domain's parent.
!***  Note that this data was first loaded into the export state of
!***  the Parent-Child coupler in subroutine EXPORT_CHILD_BOUNDARY.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      INTEGER,INTENT(IN) :: LNSH                                        &  !<-- # of boundary blending rows for H points
                           ,LNSV                                        &  !<-- # of boundary blending rows for V points
                           ,NLEV_H,NLEV_V                               &  !<-- Total # of levels in H-pt,V-pt BC vbls
                           ,NVARS_BC_2D_H,NVARS_BC_3D_H,NVARS_BC_4D_H   &  !<-- # of multi-dim H-pt boundary variables
                           ,NVARS_BC_2D_V,NVARS_BC_3D_V                 &  !<-- # of multi-dim V-pt boundary variables
                           ,PARENT_CHILD_TIME_RATIO                        !<-- # of child timesteps per parent timestep
!
      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE                             &  !
                           ,IMS,IME,JMS,JME                             &  !<-- Array dimensions
                           ,ITS,ITE,JTS,JTE                             &  !
                           ,LM                                             !
!
      REAL,INTENT(IN) :: DT                                                !<-- This domain's fundamental timestep
!
      LOGICAL(kind=KLOG),INTENT(IN) :: E_BDY,N_BDY,S_BDY,W_BDY             !<-- Is this task on any side of its domain boundary?
!
      TYPE(ESMF_State),INTENT(INOUT) :: IMP_STATE                          !<-- Solver import state
!
      TYPE(BC_H_ALL),INTENT(INOUT) :: BND_VARS_H                           !<-- All H-pt boundary data/tendencies
!
      TYPE(BC_V_ALL),INTENT(INOUT) :: BND_VARS_V                           !<-- All V-pt boundary data/tendencies
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: I1,I2_H,I2_V,J1,J2_H,J2_V
!
      INTEGER(kind=KINT) :: KOUNT_S_H,KOUNT_S_V,KOUNT_N_H,KOUNT_N_V     &
                           ,KOUNT_W_H,KOUNT_W_V,KOUNT_E_H,KOUNT_E_V
!
      INTEGER(kind=KINT) :: I,J,K,KOUNT,LBND,NL,NV,UBND
      INTEGER(kind=KINT) :: ISTAT,RC,RC_BCT
!
      REAL,SAVE :: RECIP
!
      REAL,DIMENSION(:),ALLOCATABLE :: BND_DATA_S_H                     &
                                      ,BND_DATA_S_V                     & 
                                      ,BND_DATA_N_H                     & 
                                      ,BND_DATA_N_V                     & 
                                      ,BND_DATA_W_H                     & 
                                      ,BND_DATA_W_V                     &
                                      ,BND_DATA_E_H                     & 
                                      ,BND_DATA_E_V
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC    =ESMF_SUCCESS
      RC_BCT=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  Gridpoint index limits along the South/North and West/East
!***  boundaries for mass (H) and velocity (V) points.  Note that
!***  the boundary data goes two points into the halo.
!-----------------------------------------------------------------------
!
      I1  =MAX(ITS-2,IDS)
      I2_H=MIN(ITE+2,IDE)
      I2_V=MIN(ITE+2,IDE-1)
      J1  =MAX(JTS-2,JDS)
      J2_H=MIN(JTE+2,JDE)
      J2_V=MIN(JTE+2,JDE-1)
!
!-----------------------------------------------------------------------
!***  The following 'KOUNT' variables are the number of gridpoints
!***  on the given task subdomain's South/North/West/East boundaries
!***  for all 2-D,3-D,4-D quantities on mass and velocity points.
!-----------------------------------------------------------------------
!
      KOUNT_S_H=NLEV_H*(I2_H-I1+1)*LNSH
      KOUNT_N_H=NLEV_H*(I2_H-I1+1)*LNSH
      KOUNT_S_V=NLEV_V*(I2_V-I1+1)*LNSV
      KOUNT_N_V=NLEV_V*(I2_V-I1+1)*LNSV
      KOUNT_W_H=NLEV_H*(J2_H-J1+1)*LNSH
      KOUNT_E_H=NLEV_H*(J2_H-J1+1)*LNSH
      KOUNT_W_V=NLEV_V*(J2_V-J1+1)*LNSV
      KOUNT_E_V=NLEV_V*(J2_V-J1+1)*LNSV
!
!-----------------------------------------------------------------------
!***  Compute RECIP every time in case the sign of DT has changed 
!***  due to digital filtering.
!-----------------------------------------------------------------------
!
      RECIP=1./(DT*PARENT_CHILD_TIME_RATIO)
!
!-----------------------------------------------------------------------
!***  Unload the boundary data from the import state and compute
!***  the time tendencies for the time period spanning the number
!***  of this nest's timesteps needed to reach the end of its 
!***  parent's timestep (from which the data was sent).
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  If this is a moving nest SOLVER_RUN already knows if it moved at 
!***  the beginning of this timestep.  If it has then the import state 
!***  not only contains the usual boundary data from one parent timestep 
!***  in the future but it also contains boundary data for the current
!***  timestep for the domain's new location.  We would then need to
!***  fill the current time level of the boundary variable arrays 
!***  before differencing with the values from the future to obtain
!***  the tendencies.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
      south: IF(S_BDY)THEN
!
!-----------------------------------------------------------------------
!
!-------------
!***  South H
!-------------
!
        ALLOCATE(BND_DATA_S_H(1:KOUNT_S_H))                                !<-- For south boundary H-pt data from Solver import state 
!
        move_now_south_h: IF(MOVE_NOW)THEN
!
!-----------------------------------------------------------------------
!***  Time level 1 (current) south boundary H values for new location
!***  of this nest.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Extract South Boundary H Data in UPDATE_BC_TENDS for Time N"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeGet(state    =IMP_STATE                    &  !<-- Solver import state
                                ,name     ='SOUTH_H_Current'            &  !<-- Name of south boundary H data at time N
                                ,valueList=BND_DATA_S_H                 &  !<-- The south boundary H data at time N
                                ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BCT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          KOUNT=0
!
          IF(NVARS_BC_2D_H>0)THEN
            DO NV=1,NVARS_BC_2D_H
              DO J=1,LNSH
              DO I=I1,I2_H
                KOUNT=KOUNT+1
                BND_VARS_H%VAR_2D(NV)%SOUTH(I,J,1)=BND_DATA_S_H(KOUNT)
              ENDDO
              ENDDO
            ENDDO
          ENDIF
!
          IF(NVARS_BC_3D_H>0)THEN
            DO NV=1,NVARS_BC_3D_H
              DO K=1,LM
              DO J=1,LNSH
              DO I=I1,I2_H
                KOUNT=KOUNT+1
                BND_VARS_H%VAR_3D(NV)%SOUTH(I,J,K,1)=BND_DATA_S_H(KOUNT)
              ENDDO
              ENDDO
              ENDDO
            ENDDO
          ENDIF
!
          IF(NVARS_BC_4D_H>0)THEN
            DO NV=1,NVARS_BC_4D_H
              LBND=LBOUND(BND_VARS_H%VAR_4D(NV)%FULL_VAR,4)
              UBND=UBOUND(BND_VARS_H%VAR_4D(NV)%FULL_VAR,4)
              DO NL=LBND,UBND
              DO K=1,LM
              DO J=1,LNSH
              DO I=I1,I2_H
                KOUNT=KOUNT+1
                BND_VARS_H%VAR_4D(NV)%SOUTH(I,J,K,1,NL)=BND_DATA_S_H(KOUNT)
              ENDDO
              ENDDO
              ENDDO
              ENDDO
            ENDDO
          ENDIF
!
        ENDIF move_now_south_h
!
!-----------------------------------------------------------------------
!***  Use time level 2 (future) south boundary H values to compute
!***  new tendencies.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract South Boundary H Data in UPDATE_BC_TENDS for Time N+1"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

        CALL ESMF_AttributeGet(state    =IMP_STATE                      &  !<-- Solver import state
                              ,name     ='SOUTH_H_Future'               &  !<-- Name of south boundary H data at time N+1
                              ,valueList=BND_DATA_S_H                   &  !<-- The boundary data
                              ,rc       =RC )

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BCT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        KOUNT=0
!
        IF(NVARS_BC_2D_H>0)THEN
          DO NV=1,NVARS_BC_2D_H
            DO J=1,LNSH
            DO I=I1,I2_H
              KOUNT=KOUNT+1
              BND_VARS_H%VAR_2D(NV)%SOUTH(I,J,2)=                       &
                 (BND_DATA_S_H(KOUNT)-BND_VARS_H%VAR_2D(NV)%SOUTH(I,J,1))*RECIP
            ENDDO
            ENDDO
          ENDDO
        ENDIF
!
        IF(NVARS_BC_3D_H>0)THEN
          DO NV=1,NVARS_BC_3D_H
            DO K=1,LM
            DO J=1,LNSH
            DO I=I1,I2_H
              KOUNT=KOUNT+1
              BND_VARS_H%VAR_3D(NV)%SOUTH(I,J,K,2)=                    &
                (BND_DATA_S_H(KOUNT)-BND_VARS_H%VAR_3D(NV)%SOUTH(I,J,K,1))*RECIP
            ENDDO
            ENDDO
            ENDDO
          ENDDO
        ENDIF
!
        IF(NVARS_BC_4D_H>0)THEN
          DO NV=1,NVARS_BC_4D_H
            LBND=LBOUND(BND_VARS_H%VAR_4D(NV)%FULL_VAR,4)
            UBND=UBOUND(BND_VARS_H%VAR_4D(NV)%FULL_VAR,4)
            DO NL=LBND,UBND
            DO K=1,LM
            DO J=1,LNSH
            DO I=I1,I2_H
              KOUNT=KOUNT+1
              BND_VARS_H%VAR_4D(NV)%SOUTH(I,J,K,2,NL)=                  &
                (BND_DATA_S_H(KOUNT)-BND_VARS_H%VAR_4D(NV)%SOUTH(I,J,K,1,NL))*RECIP
            ENDDO
            ENDDO
            ENDDO
            ENDDO
          ENDDO
        ENDIF
!
        DEALLOCATE(BND_DATA_S_H)
!
!-------------
!***  South V
!-------------
!
        ALLOCATE(BND_DATA_S_V(1:KOUNT_S_V))                                !<-- For south boundary V-pt data from Solver import state
!
        move_now_south_v: IF(MOVE_NOW)THEN
!
!-----------------------------------------------------------------------
!***  Time level 1 (current) south boundary V values for new location
!***  of this nest.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Extract South Boundary V Data in UPDATE_BC_TENDS for Time N"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeGet(state    =IMP_STATE                    &  !<-- Solver import state
                                ,name     ='SOUTH_V_Current'            &  !<-- Name of south boundary V data at time N
                                ,valueList=BND_DATA_S_V                 &  !<-- The south boundary V data at time N
                                ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BCT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          KOUNT=0
!
          IF(NVARS_BC_2D_V>0)THEN
            DO NV=1,NVARS_BC_2D_V
              DO J=1,LNSV
              DO I=I1,I2_V
                KOUNT=KOUNT+1
                BND_VARS_V%VAR_2D(NV)%SOUTH(I,J,1)=BND_DATA_S_V(KOUNT)
              ENDDO
              ENDDO
            ENDDO
          ENDIF
!
          IF(NVARS_BC_3D_V>0)THEN
            DO NV=1,NVARS_BC_3D_V
              DO K=1,LM
              DO J=1,LNSV
              DO I=I1,I2_V
                KOUNT=KOUNT+1
                BND_VARS_V%VAR_3D(NV)%SOUTH(I,J,K,1)=BND_DATA_S_V(KOUNT)
              ENDDO
              ENDDO
              ENDDO
            ENDDO
          ENDIF
!
        ENDIF move_now_south_v
!
!-----------------------------------------------------------------------
!***  Use time level 2 (future) south boundary V values to compute
!***  new tendencies.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract South Boundary V Data in UPDATE_BC_TENDS"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

        CALL ESMF_AttributeGet(state    =IMP_STATE                      &  !<-- Solver import state
                              ,name     ='SOUTH_V_Future'               &  !<-- Name of south boundary V data at time N+1
                              ,valueList=BND_DATA_S_V                   &  !<-- The boundary data
                              ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BCT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        KOUNT=0
!
        IF(NVARS_BC_2D_V>0)THEN
          DO NV=1,NVARS_BC_2D_V
            DO J=1,LNSV
            DO I=I1,I2_V
              KOUNT=KOUNT+1
              BND_VARS_V%VAR_2D(NV)%SOUTH(I,J,2)=                       &
                 (BND_DATA_S_V(KOUNT)-BND_VARS_V%VAR_2D(NV)%SOUTH(I,J,1))*RECIP
            ENDDO
            ENDDO
          ENDDO
        ENDIF
!
        IF(NVARS_BC_3D_V>0)THEN
          DO NV=1,NVARS_BC_3D_V
            DO K=1,LM
            DO J=1,LNSV
            DO I=I1,I2_V
              KOUNT=KOUNT+1
              BND_VARS_V%VAR_3D(NV)%SOUTH(I,J,K,2)=                    &
                (BND_DATA_S_V(KOUNT)-BND_VARS_V%VAR_3D(NV)%SOUTH(I,J,K,1))*RECIP
            ENDDO
            ENDDO
            ENDDO
          ENDDO
        ENDIF
!
        DEALLOCATE(BND_DATA_S_V)
!
      ENDIF south
!
!-----------------------------------------------------------------------
!
      north: IF(N_BDY)THEN
!
!-----------------------------------------------------------------------
!
!-------------
!***  North H
!-------------
!
        ALLOCATE(BND_DATA_N_H(1:KOUNT_N_H),stat=ISTAT)                     !<-- For north boundary H-pt data from Solver import state
!
        move_now_north_h: IF(MOVE_NOW)THEN
!
!-----------------------------------------------------------------------
!***  Time level 1 (current) north boundary H values for new location
!***  of this nest.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Extract North Boundary H Data in UPDATE_BC_TENDS for Time N"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeGet(state    =IMP_STATE                    &  !<-- Solver import state
                                ,name     ='NORTH_H_Current'            &  !<-- Name of north boundary H data at time N
                                ,valueList=BND_DATA_N_H                 &  !<-- The north boundary H data at time N
                                ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BCT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          KOUNT=0
!
          IF(NVARS_BC_2D_H>0)THEN
            DO NV=1,NVARS_BC_2D_H
              DO J=1,LNSH
              DO I=I1,I2_H
                KOUNT=KOUNT+1
                BND_VARS_H%VAR_2D(NV)%NORTH(I,J,1)=BND_DATA_N_H(KOUNT)
              ENDDO
              ENDDO
            ENDDO
          ENDIF
!
          IF(NVARS_BC_3D_H>0)THEN
            DO NV=1,NVARS_BC_3D_H
              DO K=1,LM
              DO J=1,LNSH
              DO I=I1,I2_H
                KOUNT=KOUNT+1
                BND_VARS_H%VAR_3D(NV)%NORTH(I,J,K,1)=BND_DATA_N_H(KOUNT)
              ENDDO
              ENDDO
              ENDDO
            ENDDO
          ENDIF
!
          IF(NVARS_BC_4D_H>0)THEN
            DO NV=1,NVARS_BC_4D_H
              LBND=LBOUND(BND_VARS_H%VAR_4D(NV)%FULL_VAR,4)
              UBND=UBOUND(BND_VARS_H%VAR_4D(NV)%FULL_VAR,4)
              DO NL=LBND,UBND
              DO K=1,LM
              DO J=1,LNSH
              DO I=I1,I2_H
                KOUNT=KOUNT+1
                BND_VARS_H%VAR_4D(NV)%NORTH(I,J,K,1,NL)=BND_DATA_N_H(KOUNT)
              ENDDO
              ENDDO
              ENDDO
              ENDDO
            ENDDO
          ENDIF
!
        ENDIF move_now_north_h
!
!-----------------------------------------------------------------------
!***  Use time level 2 (future) north boundary H values to compute
!***  new tendencies.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract North Boundary H Data in UPDATE_BC_TENDS for time N+1"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state    =IMP_STATE                      &  !<-- Solver import state
                              ,name     ='NORTH_H_Future'               &  !<-- Name of north boundary H data for time N+1
                              ,valueList=BND_DATA_N_H                   &  !<-- The boundary data
                              ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BCT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        KOUNT=0
!
        IF(NVARS_BC_2D_H>0)THEN
          DO NV=1,NVARS_BC_2D_H
            DO J=1,LNSH
            DO I=I1,I2_H
              KOUNT=KOUNT+1
              BND_VARS_H%VAR_2D(NV)%NORTH(I,J,2)=                       &
                 (BND_DATA_N_H(KOUNT)-BND_VARS_H%VAR_2D(NV)%NORTH(I,J,1))*RECIP
            ENDDO
            ENDDO
          ENDDO
        ENDIF
!
        IF(NVARS_BC_3D_H>0)THEN
          DO NV=1,NVARS_BC_3D_H
            DO K=1,LM
            DO J=1,LNSH
            DO I=I1,I2_H
              KOUNT=KOUNT+1
              BND_VARS_H%VAR_3D(NV)%NORTH(I,J,K,2)=                    &
                (BND_DATA_N_H(KOUNT)-BND_VARS_H%VAR_3D(NV)%NORTH(I,J,K,1))*RECIP
            ENDDO
            ENDDO
            ENDDO
          ENDDO
        ENDIF
!
        IF(NVARS_BC_4D_H>0)THEN
          DO NV=1,NVARS_BC_4D_H
            LBND=LBOUND(BND_VARS_H%VAR_4D(NV)%FULL_VAR,4)
            UBND=UBOUND(BND_VARS_H%VAR_4D(NV)%FULL_VAR,4)
            DO NL=LBND,UBND
            DO K=1,LM
            DO J=1,LNSH
            DO I=I1,I2_H
              KOUNT=KOUNT+1
              BND_VARS_H%VAR_4D(NV)%NORTH(I,J,K,2,NL)=                  &
                (BND_DATA_N_H(KOUNT)-BND_VARS_H%VAR_4D(NV)%NORTH(I,J,K,1,NL))*RECIP
            ENDDO
            ENDDO
            ENDDO
            ENDDO
          ENDDO
        ENDIF
!
        DEALLOCATE(BND_DATA_N_H)
!
!-------------
!***  North V
!-------------
!
        ALLOCATE(BND_DATA_N_V(1:KOUNT_N_V))                                !<-- For north boundary V-pt data from Solver import state
!
        move_now_north_v: IF(MOVE_NOW)THEN
!
!-----------------------------------------------------------------------
!***  Time level 1 (current) north boundary V values for new location
!***  of this nest.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Extract North Boundary V Data in UPDATE_BC_TENDS for Time N"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeGet(state    =IMP_STATE                    &  !<-- Solver import state
                                ,name     ='NORTH_V_Current'            &  !<-- Name of north boundary V data at time N
                                ,valueList=BND_DATA_N_V                 &  !<-- The north boundary V data at time N
                                ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BCT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          KOUNT=0
!
          IF(NVARS_BC_2D_V>0)THEN
            DO NV=1,NVARS_BC_2D_V
              DO J=1,LNSV
              DO I=I1,I2_V
                KOUNT=KOUNT+1
                BND_VARS_V%VAR_2D(NV)%NORTH(I,J,1)=BND_DATA_S_V(KOUNT)
              ENDDO
              ENDDO
            ENDDO
          ENDIF
!
          IF(NVARS_BC_3D_V>0)THEN
            DO NV=1,NVARS_BC_3D_V
              DO K=1,LM
              DO J=1,LNSV
              DO I=I1,I2_V
                KOUNT=KOUNT+1
                BND_VARS_V%VAR_3D(NV)%NORTH(I,J,K,1)=BND_DATA_N_V(KOUNT)
              ENDDO
              ENDDO
              ENDDO
            ENDDO
          ENDIF
!
        ENDIF move_now_north_v
!
!-----------------------------------------------------------------------
!***  Use time level 2 (future) north boundary H values to compute
!***  new tendencies.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract North Boundary V Data in UPDATE_BC_TENDS for Time N+1"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

        CALL ESMF_AttributeGet(state    =IMP_STATE                      &  !<-- Solver import state
                              ,name     ='NORTH_V_Future'               &  !<-- Name of north boundary V data at time N+1
                              ,valueList=BND_DATA_N_V                   &  !<-- The boundary data
                              ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BCT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        KOUNT=0
!
        IF(NVARS_BC_2D_V>0)THEN
          DO NV=1,NVARS_BC_2D_V
            DO J=1,LNSV
            DO I=I1,I2_V
              KOUNT=KOUNT+1
              BND_VARS_V%VAR_2D(NV)%NORTH(I,J,2)=                       &
                 (BND_DATA_N_V(KOUNT)-BND_VARS_V%VAR_2D(NV)%NORTH(I,J,1))*RECIP
            ENDDO
            ENDDO
          ENDDO
        ENDIF
!
        IF(NVARS_BC_3D_V>0)THEN
          DO NV=1,NVARS_BC_3D_V
            DO K=1,LM
            DO J=1,LNSV
            DO I=I1,I2_V
              KOUNT=KOUNT+1
              BND_VARS_V%VAR_3D(NV)%NORTH(I,J,K,2)=                    &
                (BND_DATA_N_V(KOUNT)-BND_VARS_V%VAR_3D(NV)%NORTH(I,J,K,1))*RECIP
            ENDDO
            ENDDO
            ENDDO
          ENDDO
        ENDIF
!
        DEALLOCATE(BND_DATA_N_V)
!
      ENDIF north
!
!-----------------------------------------------------------------------
!
      west: IF(W_BDY)THEN
!
!-----------------------------------------------------------------------
!
!------------
!***  West H
!------------
!
        ALLOCATE(BND_DATA_W_H(1:KOUNT_W_H))                                !<-- For west boundary H-pt data from Solver import state
!
        move_now_west_h: IF(MOVE_NOW)THEN
!
!-----------------------------------------------------------------------
!***  Time level 1 (current) west boundary H values for new location
!***  of this nest.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Extract West Boundary H Data in UPDATE_BC_TENDS for Time N"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeGet(state    =IMP_STATE                    &  !<-- Solver import state
                                ,name     ='WEST_H_Current'             &  !<-- Name of west boundary H data at time N
                                ,valueList=BND_DATA_W_H                 &  !<-- The west boundary H data at time N
                                ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BCT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          KOUNT=0
!
          IF(NVARS_BC_2D_H>0)THEN
            DO NV=1,NVARS_BC_2D_H
              DO J=J1,J2_H
              DO I=1,LNSH
                KOUNT=KOUNT+1
                BND_VARS_H%VAR_2D(NV)%WEST(I,J,1)=BND_DATA_W_H(KOUNT)
              ENDDO
              ENDDO
            ENDDO
          ENDIF
!
          IF(NVARS_BC_3D_H>0)THEN
            DO NV=1,NVARS_BC_3D_H
              DO K=1,LM
              DO J=J1,J2_H
              DO I=1,LNSH
                KOUNT=KOUNT+1
                BND_VARS_H%VAR_3D(NV)%WEST(I,J,K,1)=BND_DATA_W_H(KOUNT)
              ENDDO
              ENDDO
              ENDDO
            ENDDO
          ENDIF
!
          IF(NVARS_BC_4D_H>0)THEN
            DO NV=1,NVARS_BC_4D_H
              LBND=LBOUND(BND_VARS_H%VAR_4D(NV)%FULL_VAR,4)
              UBND=UBOUND(BND_VARS_H%VAR_4D(NV)%FULL_VAR,4)
              DO NL=LBND,UBND
              DO K=1,LM
              DO J=J1,J2_H
              DO I=1,LNSH
                KOUNT=KOUNT+1
                BND_VARS_H%VAR_4D(NV)%WEST(I,J,K,1,NL)=BND_DATA_W_H(KOUNT)
              ENDDO
              ENDDO
              ENDDO
              ENDDO
            ENDDO
          ENDIF
!
        ENDIF move_now_west_h
!
!-----------------------------------------------------------------------
!***  Use time level 2 (future) west boundary H values to compute
!***  new tendencies.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract West Boundary H Data in UPDATE_BC_TENDS at Time N+1"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state    =IMP_STATE                      &  !<-- Solver import state
                              ,name     ='WEST_H_Future'                &  !<-- Name of west boundary H data at time N+1
                              ,valueList=BND_DATA_W_H                   &  !<-- The boundary data
                              ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BCT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        KOUNT=0
!
        IF(NVARS_BC_2D_H>0)THEN
          DO NV=1,NVARS_BC_2D_H
            DO J=J1,J2_H
            DO I=1,LNSH
              KOUNT=KOUNT+1
              BND_VARS_H%VAR_2D(NV)%WEST(I,J,2)=                        &
                 (BND_DATA_W_H(KOUNT)-BND_VARS_H%VAR_2D(NV)%WEST(I,J,1))*RECIP
            ENDDO
            ENDDO
          ENDDO
        ENDIF
!
        IF(NVARS_BC_3D_H>0)THEN
          DO NV=1,NVARS_BC_3D_H
            DO K=1,LM
            DO J=J1,J2_H
            DO I=1,LNSH
              KOUNT=KOUNT+1
              BND_VARS_H%VAR_3D(NV)%WEST(I,J,K,2)=                     &
                (BND_DATA_W_H(KOUNT)-BND_VARS_H%VAR_3D(NV)%WEST(I,J,K,1))*RECIP
            ENDDO
            ENDDO
            ENDDO
          ENDDO
        ENDIF
!
        IF(NVARS_BC_4D_H>0)THEN
          DO NV=1,NVARS_BC_4D_H
            LBND=LBOUND(BND_VARS_H%VAR_4D(NV)%FULL_VAR,4)
            UBND=UBOUND(BND_VARS_H%VAR_4D(NV)%FULL_VAR,4)
            DO NL=LBND,UBND
            DO K=1,LM
            DO J=J1,J2_H
            DO I=1,LNSH
              KOUNT=KOUNT+1
              BND_VARS_H%VAR_4D(NV)%WEST(I,J,K,2,NL)=                   &
                (BND_DATA_W_H(KOUNT)-BND_VARS_H%VAR_4D(NV)%WEST(I,J,K,1,NL))*RECIP
            ENDDO
            ENDDO
            ENDDO
            ENDDO
          ENDDO
        ENDIF
!
        DEALLOCATE(BND_DATA_W_H)
!
!------------
!***  West V
!------------
!
        ALLOCATE(BND_DATA_W_V(1:KOUNT_W_V))                                !<-- For west boundary V-pt data from Solver import state
!
        move_now_west_v: IF(MOVE_NOW)THEN
!
!-----------------------------------------------------------------------
!***  Time level 1 (current) west boundary V values for new location
!***  of this nest.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Extract West Boundary V Data in UPDATE_BC_TENDS for Time N"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeGet(state    =IMP_STATE                    &  !<-- Solver import state
                                ,name     ='WEST_V_Current'             &  !<-- Name of west boundary V data at time N
                                ,valueList=BND_DATA_W_V                 &  !<-- The west boundary V data at time N
                                ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BCT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          KOUNT=0
!
          IF(NVARS_BC_2D_V>0)THEN
            DO NV=1,NVARS_BC_2D_V
              DO J=J1,J2_V
              DO I=1,LNSV
                KOUNT=KOUNT+1
                BND_VARS_V%VAR_2D(NV)%WEST(I,J,1)=BND_DATA_W_V(KOUNT)
              ENDDO
              ENDDO
            ENDDO
          ENDIF
!
          IF(NVARS_BC_3D_V>0)THEN
            DO NV=1,NVARS_BC_3D_V
              DO K=1,LM
              DO J=J1,J2_V
              DO I=1,LNSV
                KOUNT=KOUNT+1
                BND_VARS_V%VAR_3D(NV)%WEST(I,J,K,1)=BND_DATA_W_V(KOUNT)
              ENDDO
              ENDDO
              ENDDO
            ENDDO
          ENDIF
!
        ENDIF move_now_west_v
!
!-----------------------------------------------------------------------
!***  Use time level 2 (future) west boundary V values to compute
!***  new tendencies.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract West Boundary V Data in UPDATE_BC_TENDS at Time N+1"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state    =IMP_STATE                      &  !<-- Solver import state
                              ,name     ='WEST_V_Future'                &  !<-- Name of west boundary V data at time N+1
                              ,valueList=BND_DATA_W_V                   &  !<-- The boundary data
                              ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BCT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        KOUNT=0
!
        IF(NVARS_BC_2D_V>0)THEN
          DO NV=1,NVARS_BC_2D_V
            DO J=J1,J2_V
            DO I=1,LNSV
              KOUNT=KOUNT+1
              BND_VARS_V%VAR_2D(NV)%WEST(I,J,2)=                        &
                 (BND_DATA_W_V(KOUNT)-BND_VARS_V%VAR_2D(NV)%WEST(I,J,1))*RECIP
            ENDDO
            ENDDO
          ENDDO
        ENDIF
!
        IF(NVARS_BC_3D_V>0)THEN
          DO NV=1,NVARS_BC_3D_V
            DO K=1,LM
            DO J=J1,J2_V
            DO I=1,LNSV
              KOUNT=KOUNT+1
              BND_VARS_V%VAR_3D(NV)%WEST(I,J,K,2)=                     &
                (BND_DATA_W_V(KOUNT)-BND_VARS_V%VAR_3D(NV)%WEST(I,J,K,1))*RECIP
            ENDDO
            ENDDO
            ENDDO
          ENDDO
        ENDIF
!
        DEALLOCATE(BND_DATA_W_V)
!
      ENDIF west
!
!-----------------------------------------------------------------------
!
      east: IF(E_BDY)THEN
!
!-----------------------------------------------------------------------
!
!------------
!***  East H
!------------
!
        ALLOCATE(BND_DATA_E_H(1:KOUNT_E_H))                                !<-- For east boundary H-pt data from Solver import state
!
        move_now_east_h: IF(MOVE_NOW)THEN
!
!-----------------------------------------------------------------------
!***  Time level 1 (current) east boundary H values for new location
!***  of this nest.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Extract East Boundary H Data in UPDATE_BC_TENDS for Time N"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeGet(state    =IMP_STATE                    &  !<-- Solver import state
                                ,name     ='EAST_H_Current'             &  !<-- Name of east boundary H data at time N
                                ,valueList=BND_DATA_E_H                 &  !<-- The east boundary H data at time N
                                ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BCT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          KOUNT=0
!
          IF(NVARS_BC_2D_H>0)THEN
            DO NV=1,NVARS_BC_2D_H
              DO J=J1,J2_H
              DO I=1,LNSH
                KOUNT=KOUNT+1
                BND_VARS_H%VAR_2D(NV)%EAST(I,J,1)=BND_DATA_E_H(KOUNT)
              ENDDO
              ENDDO
            ENDDO
          ENDIF
!
          IF(NVARS_BC_3D_H>0)THEN
            DO NV=1,NVARS_BC_3D_H
              DO K=1,LM
              DO J=J1,J2_H
              DO I=1,LNSH
                KOUNT=KOUNT+1
                BND_VARS_H%VAR_3D(NV)%EAST(I,J,K,1)=BND_DATA_E_H(KOUNT)
              ENDDO
              ENDDO
              ENDDO
            ENDDO
          ENDIF
!
          IF(NVARS_BC_4D_H>0)THEN
            DO NV=1,NVARS_BC_4D_H
              LBND=LBOUND(BND_VARS_H%VAR_4D(NV)%FULL_VAR,4)
              UBND=UBOUND(BND_VARS_H%VAR_4D(NV)%FULL_VAR,4)
              DO NL=LBND,UBND
              DO K=1,LM
              DO J=J1,J2_H
              DO I=1,LNSH
                KOUNT=KOUNT+1
                BND_VARS_H%VAR_4D(NV)%EAST(I,J,K,1,NL)=BND_DATA_E_H(KOUNT)
              ENDDO
              ENDDO
              ENDDO
              ENDDO
            ENDDO
          ENDIF
!
        ENDIF move_now_east_h
!
!-----------------------------------------------------------------------
!***  Use time level 2 (future) east boundary H values to compute
!***  new tendencies.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract East Boundary H Data in UPDATE_BC_TENDS at Time N+1"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state    =IMP_STATE                      &  !<-- Solver import state
                              ,name     ='EAST_H_Future'                &  !<-- Name of east boundary H data at time N+1
                              ,valueList=BND_DATA_E_H                   &  !<-- The boundary data
                              ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BCT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        KOUNT=0
!
        IF(NVARS_BC_2D_H>0)THEN
          DO NV=1,NVARS_BC_2D_H
            DO J=J1,J2_H
            DO I=1,LNSH
              KOUNT=KOUNT+1
              BND_VARS_H%VAR_2D(NV)%EAST(I,J,2)=                        &
                 (BND_DATA_E_H(KOUNT)-BND_VARS_H%VAR_2D(NV)%EAST(I,J,1))*RECIP
            ENDDO
            ENDDO
          ENDDO
        ENDIF
!
        IF(NVARS_BC_3D_H>0)THEN
          DO NV=1,NVARS_BC_3D_H
            DO K=1,LM
            DO J=J1,J2_H
            DO I=1,LNSH
              KOUNT=KOUNT+1
              BND_VARS_H%VAR_3D(NV)%EAST(I,J,K,2)=                     &
                (BND_DATA_E_H(KOUNT)-BND_VARS_H%VAR_3D(NV)%EAST(I,J,K,1))*RECIP
            ENDDO
            ENDDO
            ENDDO
          ENDDO
        ENDIF
!
        IF(NVARS_BC_4D_H>0)THEN
          DO NV=1,NVARS_BC_4D_H
            LBND=LBOUND(BND_VARS_H%VAR_4D(NV)%FULL_VAR,4)
            UBND=UBOUND(BND_VARS_H%VAR_4D(NV)%FULL_VAR,4)
            DO NL=LBND,UBND
            DO K=1,LM
            DO J=J1,J2_H
            DO I=1,LNSH
              KOUNT=KOUNT+1
              BND_VARS_H%VAR_4D(NV)%EAST(I,J,K,2,NL)=                   &
                (BND_DATA_E_H(KOUNT)-BND_VARS_H%VAR_4D(NV)%EAST(I,J,K,1,NL))*RECIP
            ENDDO
            ENDDO
            ENDDO
            ENDDO
          ENDDO
        ENDIF
!
        DEALLOCATE(BND_DATA_E_H)
!
!------------
!***  East V
!------------
!
        ALLOCATE(BND_DATA_E_V(1:KOUNT_E_V))                                !<-- For east boundary V-pt data from Solver import state
!
        move_now_east_v: IF(MOVE_NOW)THEN
!
!-----------------------------------------------------------------------
!***  Time level 1 (current) east boundary V values for new location
!***  of this nest.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Extract East Boundary V Data in UPDATE_BC_TENDS for Time N"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeGet(state    =IMP_STATE                    &  !<-- Solver import state
                                ,name     ='EAST_V_Current'             &  !<-- Name of esat boundary V data at time N
                                ,valueList=BND_DATA_E_V                 &  !<-- The east boundary V data at time N
                                ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BCT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          KOUNT=0
!
          IF(NVARS_BC_2D_V>0)THEN
            DO NV=1,NVARS_BC_2D_V
              DO J=J1,J2_V
              DO I=1,LNSV
                KOUNT=KOUNT+1
                BND_VARS_V%VAR_2D(NV)%EAST(I,J,1)=BND_DATA_E_V(KOUNT)
              ENDDO
              ENDDO
            ENDDO
          ENDIF
!
          IF(NVARS_BC_3D_V>0)THEN
            DO NV=1,NVARS_BC_3D_V
              DO K=1,LM
              DO J=J1,J2_V
              DO I=1,LNSV
                KOUNT=KOUNT+1
                BND_VARS_V%VAR_3D(NV)%EAST(I,J,K,1)=BND_DATA_E_V(KOUNT)
              ENDDO
              ENDDO
              ENDDO
            ENDDO
          ENDIF
!
        ENDIF move_now_east_v
!
!-----------------------------------------------------------------------
!***  Use time level 2 (future) east boundary V values to compute
!***  new tendencies.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract East Boundary V Data in UPDATE_BC_TENDS for Time N+1"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state    =IMP_STATE                      &  !<-- Solver import state
                              ,name     ='EAST_V_Future'                &  !<-- Name of east boundary V data at time N+1
                              ,valueList=BND_DATA_E_V                   &  !<-- The boundary data
                              ,rc       =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BCT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        KOUNT=0
!
        IF(NVARS_BC_2D_V>0)THEN
          DO NV=1,NVARS_BC_2D_V
            DO J=J1,J2_V
            DO I=1,LNSV
              KOUNT=KOUNT+1
              BND_VARS_V%VAR_2D(NV)%EAST(I,J,2)=                        &
                 (BND_DATA_E_V(KOUNT)-BND_VARS_V%VAR_2D(NV)%EAST(I,J,1))*RECIP
            ENDDO
            ENDDO
          ENDDO
        ENDIF
!
        IF(NVARS_BC_3D_V>0)THEN
          DO NV=1,NVARS_BC_3D_V
            DO K=1,LM
            DO J=J1,J2_V
            DO I=1,LNSV
              KOUNT=KOUNT+1
              BND_VARS_V%VAR_3D(NV)%EAST(I,J,K,2)=                     &
                (BND_DATA_E_V(KOUNT)-BND_VARS_V%VAR_3D(NV)%EAST(I,J,K,1))*RECIP
            ENDDO
            ENDDO
            ENDDO
          ENDDO
        ENDIF
!
        DEALLOCATE(BND_DATA_E_V)
!
      ENDIF east
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE UPDATE_BC_TENDS
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE SAVE_BC_DATA(LM,LNSH,LNSV                              &
                             ,NVARS_BC_2D_H                             &
                             ,NVARS_BC_3D_H                             &
                             ,NVARS_BC_4D_H                             &
                             ,NVARS_BC_2D_V                             &
                             ,NVARS_BC_3D_V                             &
                             ,BND_VARS_H                                &
                             ,BND_VARS_V                                &
                             ,NUM_WORDS_BC_SOUTH,RST_BC_DATA_SOUTH      &
                             ,NUM_WORDS_BC_NORTH,RST_BC_DATA_NORTH      &
                             ,NUM_WORDS_BC_WEST ,RST_BC_DATA_WEST       &
                             ,NUM_WORDS_BC_EAST ,RST_BC_DATA_EAST       &
                             ,EXP_STATE_SOLVER                          &
                             ,ITS,ITE,JTS,JTE                           &
                             ,IMS,IME,JMS,JME                           &
                             ,IDS,IDE,JDS,JDE                           &
                               )
! 
!-----------------------------------------------------------------------
!***  Boundary array winds are needed in the restart file in order to
!***  achieve bit identical answers between restarted runs and their
!***  free-forecast analogs.  The boundary arrays do not span the
!***  integration grid thus they can only be transmitted through
!***  ESMF States as Attributes.  Non-scalar Attributes can only
!***  contain one dimension therefore the boundary data is moved
!***  into 1-D arrays in this routine then inserted into the
!***  Write component's import state.
!-----------------------------------------------------------------------
!
!---------------------
!***  Input Arguments
!---------------------
!
      INTEGER(kind=KINT),INTENT(IN) :: LNSH                             &  !<-- # of boundary blending rows for H points
                                      ,LNSV                             &  !<-- # of boundary blending rows for V points
                                      ,NUM_WORDS_BC_SOUTH               &  !<-- Total # of words in south bndry winds, this fcst task
                                      ,NUM_WORDS_BC_NORTH               &  !<-- Total # of words in north bndry winds, this fcst task
                                      ,NUM_WORDS_BC_WEST                &  !<-- Total # of words in west bndry winds, this fcst task
                                      ,NUM_WORDS_BC_EAST                   !<-- Total # of words in east bndry winds, this fcst task
!
      INTEGER(kind=KINT),INTENT(IN) :: NVARS_BC_2D_H                    &
                                      ,NVARS_BC_3D_H                    &
                                      ,NVARS_BC_4D_H                    &
                                      ,NVARS_BC_2D_V                    &
                                      ,NVARS_BC_3D_V
!
      INTEGER(kind=KINT),INTENT(IN) :: IDS,IDE,JDS,JDE                  &  !<-- 
                                      ,IMS,IME,JMS,JME                  &  !<-- Array dimensions
                                      ,ITS,ITE,JTS,JTE                  &  !<-- 
                                      ,LM                                  !<--
!
      TYPE(BC_H_ALL),INTENT(IN) :: BND_VARS_H                              !<-- All H-pt boundary data/tendencies
!
      TYPE(BC_V_ALL),INTENT(IN) :: BND_VARS_V                              !<-- All V-pt boundary data/tendencies
!
!---------------------
!***  Inout Arguments
!---------------------
!
      TYPE(ESMF_State),INTENT(INOUT) :: EXP_STATE_SOLVER                   !<-- The Solver export state
!
!----------------------
!***  Output Arguments
!----------------------
!
      REAL(kind=KFPT),DIMENSION(1:NUM_WORDS_BC_SOUTH),INTENT(OUT) ::    &
                                                     RST_BC_DATA_SOUTH     !<-- All south bndry wind data on this fcst task
      REAL(kind=KFPT),DIMENSION(1:NUM_WORDS_BC_NORTH),INTENT(OUT) ::    &
                                                     RST_BC_DATA_NORTH     !<-- All north bndry wind data on this fcst task
      REAL(kind=KFPT),DIMENSION(1:NUM_WORDS_BC_WEST ),INTENT(OUT) ::    &
                                                     RST_BC_DATA_WEST      !<-- All west bndry wind data on this fcst task
      REAL(kind=KFPT),DIMENSION(1:NUM_WORDS_BC_EAST ),INTENT(OUT) ::    &
                                                     RST_BC_DATA_EAST      !<-- All east bndry wind data on this fcst task
!
!-----------------------------------------------------------------------
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: IB,JB,KOUNT,L,LBND,NL,NT,NV,UBND
!
      INTEGER(kind=KINT) :: RC,RC_SAVE
!
      TYPE(ESMF_State) :: IMP_STATE_WRITE
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Southern boundary data to 1-D
!-----------------------------------------------------------------------
!
      IF(JTS==JDS)THEN                                                     !<-- Tasks on south boundary
!
        KOUNT=0
!
        IF(NVARS_BC_2D_H>0)THEN
          DO NV=1,NVARS_BC_2D_H
            DO NT=1,2
            DO JB=1,LNSH
            DO IB=ITS,ITE
              KOUNT=KOUNT+1
              RST_BC_DATA_SOUTH(KOUNT)=BND_VARS_H%VAR_2D(NV)%SOUTH(IB,JB,NT)
            ENDDO
            ENDDO
            ENDDO
          ENDDO
        ENDIF
!
        IF(NVARS_BC_3D_H>0)THEN
          DO NV=1,NVARS_BC_3D_H
            DO NT=1,2
            DO L=1,LM
              DO JB=1,LNSH
              DO IB=ITS,ITE
                KOUNT=KOUNT+1
                RST_BC_DATA_SOUTH(KOUNT)=BND_VARS_H%VAR_3D(NV)%SOUTH(IB,JB,L,NT)
              ENDDO
              ENDDO
            ENDDO
            ENDDO
          ENDDO
        ENDIF
!
        IF(NVARS_BC_4D_H>0)THEN
          DO NV=1,NVARS_BC_4D_H
            LBND=LBOUND(BND_VARS_H%VAR_4D(NV)%FULL_VAR,4)
            UBND=UBOUND(BND_VARS_H%VAR_4D(NV)%FULL_VAR,4)
            DO NL=LBND,UBND
              DO NT=1,2
              DO L=1,LM
                DO JB=1,LNSH
                DO IB=ITS,ITE
                  KOUNT=KOUNT+1
                  RST_BC_DATA_SOUTH(KOUNT)=BND_VARS_H%VAR_4D(NV)%SOUTH(IB,JB,L,NT,NL)
                ENDDO
                ENDDO
              ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF
!
        IF(NVARS_BC_2D_V>0)THEN
          DO NV=1,NVARS_BC_2D_V
            DO NT=1,2
            DO JB=1,LNSV
            DO IB=ITS,ITE
              KOUNT=KOUNT+1
              RST_BC_DATA_SOUTH(KOUNT)=BND_VARS_V%VAR_2D(NV)%SOUTH(IB,JB,NT)
            ENDDO
            ENDDO
            ENDDO
          ENDDO
        ENDIF
!
        IF(NVARS_BC_3D_V>0)THEN
          DO NV=1,NVARS_BC_3D_V
            DO NT=1,2
            DO L=1,LM
              DO JB=1,LNSV
              DO IB=ITS,ITE
                KOUNT=KOUNT+1
                RST_BC_DATA_SOUTH(KOUNT)=BND_VARS_V%VAR_3D(NV)%SOUTH(IB,JB,L,NT)
              ENDDO
              ENDDO
            ENDDO
            ENDDO
          ENDDO
        ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract Write Import State in SAVE_BC_DATA"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_StateGet(state      =EXP_STATE_SOLVER                 &  !<-- The Solver export state
                          ,itemName   ='Write Import State'             &  !<-- Name of the state to get from Solver export state
                          ,nestedState=IMP_STATE_WRITE                  &  !<-- Extract Write Component import state from Solver export
                          ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_SAVE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Set BC South Data Attribute in SAVE_BC_DATA"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

        CALL ESMF_AttributeSet(state    =IMP_STATE_WRITE                &  !<-- The Write component import state
                              ,name     ='RST_BC_DATA_SOUTH'            &  !<-- Name of 1-D string of south boundary values
                              ,itemCount=NUM_WORDS_BC_SOUTH             &  !<-- # of south boundary words on this fcst task
                              ,valueList=RST_BC_DATA_SOUTH              &  !<-- The 1-D data being inserted into the Write import state
                              ,rc       =RC)

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_SAVE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF
!
!-----------------------------------------------------------------------
!***  Northern boundary data to 1-D
!-----------------------------------------------------------------------
!
      IF(JTE==JDE)THEN                                                     !<-- Tasks on north boundary
!
        KOUNT=0
!
        IF(NVARS_BC_2D_H>0)THEN
          DO NV=1,NVARS_BC_2D_H
            DO JB=1,LNSH
            DO IB=ITS,ITE
              KOUNT=KOUNT+1
              RST_BC_DATA_NORTH(KOUNT)=BND_VARS_H%VAR_2D(NV)%NORTH(IB,JB,1)
            ENDDO
            ENDDO
            DO JB=1,LNSH
            DO IB=ITS,ITE
              KOUNT=KOUNT+1
              RST_BC_DATA_NORTH(KOUNT)=BND_VARS_H%VAR_2D(NV)%NORTH(IB,JB,2)
            ENDDO
            ENDDO
          ENDDO
        ENDIF
!
        IF(NVARS_BC_3D_H>0)THEN
          DO NV=1,NVARS_BC_3D_H
            DO L=1,LM
              DO JB=1,LNSH
              DO IB=ITS,ITE
                KOUNT=KOUNT+1
                RST_BC_DATA_NORTH(KOUNT)=BND_VARS_H%VAR_3D(NV)%NORTH(IB,JB,L,1)
              ENDDO
              ENDDO
            ENDDO
            DO L=1,LM
              DO JB=1,LNSH
              DO IB=ITS,ITE
                KOUNT=KOUNT+1
                RST_BC_DATA_NORTH(KOUNT)=BND_VARS_H%VAR_3D(NV)%NORTH(IB,JB,L,2)
              ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF
!
        IF(NVARS_BC_4D_H>0)THEN
          DO NV=1,NVARS_BC_4D_H
            LBND=LBOUND(BND_VARS_H%VAR_4D(NV)%FULL_VAR,4)
            UBND=UBOUND(BND_VARS_H%VAR_4D(NV)%FULL_VAR,4)
            DO NL=LBND,UBND
              DO L=1,LM
                DO JB=1,LNSH
                DO IB=ITS,ITE
                  KOUNT=KOUNT+1
                  RST_BC_DATA_NORTH(KOUNT)=BND_VARS_H%VAR_4D(NV)%NORTH(IB,JB,L,1,NL)
                ENDDO
                ENDDO
              ENDDO
              DO L=1,LM
                DO JB=1,LNSH
                DO IB=ITS,ITE
                  KOUNT=KOUNT+1
                  RST_BC_DATA_NORTH(KOUNT)=BND_VARS_H%VAR_4D(NV)%NORTH(IB,JB,L,2,NL)
                ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF
!
        IF(NVARS_BC_2D_V>0)THEN
          DO NV=1,NVARS_BC_2D_V
            DO JB=1,LNSV
            DO IB=ITS,ITE
              KOUNT=KOUNT+1
              RST_BC_DATA_NORTH(KOUNT)=BND_VARS_V%VAR_2D(NV)%NORTH(IB,JB,1)
            ENDDO
            ENDDO
            DO JB=1,LNSV
            DO IB=ITS,ITE
              KOUNT=KOUNT+1
              RST_BC_DATA_NORTH(KOUNT)=BND_VARS_V%VAR_2D(NV)%NORTH(IB,JB,2)
            ENDDO
            ENDDO
          ENDDO
        ENDIF
!
        IF(NVARS_BC_3D_V>0)THEN
          DO NV=1,NVARS_BC_3D_V
            DO L=1,LM
              DO JB=1,LNSV
              DO IB=ITS,ITE
                KOUNT=KOUNT+1
                RST_BC_DATA_NORTH(KOUNT)=BND_VARS_V%VAR_3D(NV)%NORTH(IB,JB,L,1)
              ENDDO
              ENDDO
            ENDDO
            DO L=1,LM
              DO JB=1,LNSV
              DO IB=ITS,ITE
                KOUNT=KOUNT+1
                RST_BC_DATA_NORTH(KOUNT)=BND_VARS_V%VAR_3D(NV)%NORTH(IB,JB,L,2)
              ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract Write Import State in SAVE_BC_DATA"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_StateGet(state      =EXP_STATE_SOLVER                 &  !<-- The Solver export state
                          ,itemName   ='Write Import State'             &  !<-- Name of the state to get from Solver export state
                          ,nestedState=IMP_STATE_WRITE                  &  !<-- Extract Write Component import state from Solver export
                          ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_SAVE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Set BC North Data Attribute in SAVE_BC_DATA"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

        CALL ESMF_AttributeSet(state    =IMP_STATE_WRITE                &  !<-- The Write component import state
                              ,name     ='RST_BC_DATA_NORTH'            &  !<-- Name of 1-D string of north boundary values
                              ,itemCount=NUM_WORDS_BC_NORTH             &  !<-- # of north boundary words on this fcst task
                              ,valueList=RST_BC_DATA_NORTH              &  !<-- The 1-D data being inserted into the Write import state
                              ,rc       =RC)

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_SAVE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF
!
!-----------------------------------------------------------------------
!***  Western boundary data to 1-D
!-----------------------------------------------------------------------
!
      IF(ITS==IDS)THEN                                                     !<-- Tasks on west boundary
!
        KOUNT=0
!
        IF(NVARS_BC_2D_H>0)THEN
          DO NV=1,NVARS_BC_2D_H
            DO JB=JTS,JTE
            DO IB=1,LNSH
              KOUNT=KOUNT+1
              RST_BC_DATA_WEST(KOUNT)=BND_VARS_H%VAR_2D(NV)%WEST(IB,JB,1)
            ENDDO
            ENDDO
            DO JB=JTS,JTE
            DO IB=1,LNSH
              KOUNT=KOUNT+1
              RST_BC_DATA_WEST(KOUNT)=BND_VARS_H%VAR_2D(NV)%WEST(IB,JB,2)
            ENDDO
            ENDDO
          ENDDO
        ENDIF
!
        IF(NVARS_BC_3D_H>0)THEN
          DO NV=1,NVARS_BC_3D_H
            DO L=1,LM
              DO JB=JTS,JTE
              DO IB=1,LNSH
                KOUNT=KOUNT+1
                RST_BC_DATA_WEST(KOUNT)=BND_VARS_H%VAR_3D(NV)%WEST(IB,JB,L,1)
              ENDDO
              ENDDO
            ENDDO
            DO L=1,LM
              DO JB=JTS,JTE
              DO IB=1,LNSH
                KOUNT=KOUNT+1
                RST_BC_DATA_WEST(KOUNT)=BND_VARS_H%VAR_3D(NV)%WEST(IB,JB,L,2)
              ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF
!
        IF(NVARS_BC_4D_H>0)THEN
          DO NV=1,NVARS_BC_4D_H
            LBND=LBOUND(BND_VARS_H%VAR_4D(NV)%FULL_VAR,4)
            UBND=UBOUND(BND_VARS_H%VAR_4D(NV)%FULL_VAR,4)
            DO NL=LBND,UBND
              DO L=1,LM
                DO JB=JTS,JTE
                DO IB=1,LNSH
                  KOUNT=KOUNT+1
                  RST_BC_DATA_WEST(KOUNT)=BND_VARS_H%VAR_4D(NV)%WEST(IB,JB,L,1,NL)
                ENDDO
                ENDDO
              ENDDO
              DO L=1,LM
                DO JB=JTS,JTE
                DO IB=1,LNSH
                  KOUNT=KOUNT+1
                  RST_BC_DATA_WEST(KOUNT)=BND_VARS_H%VAR_4D(NV)%WEST(IB,JB,L,2,NL)
                ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF
!
        IF(NVARS_BC_2D_V>0)THEN
          DO NV=1,NVARS_BC_2D_V
            DO JB=JTS,JTE
            DO IB=1,LNSV
              KOUNT=KOUNT+1
              RST_BC_DATA_WEST(KOUNT)=BND_VARS_V%VAR_2D(NV)%WEST(IB,JB,1)
            ENDDO
            ENDDO
            DO JB=JTS,JTE
            DO IB=1,LNSV
              KOUNT=KOUNT+1
              RST_BC_DATA_WEST(KOUNT)=BND_VARS_V%VAR_2D(NV)%WEST(IB,JB,2)
            ENDDO
            ENDDO
          ENDDO
        ENDIF
!
        IF(NVARS_BC_3D_V>0)THEN
          DO NV=1,NVARS_BC_3D_V
            DO L=1,LM
              DO JB=JTS,JTE
              DO IB=1,LNSV
                KOUNT=KOUNT+1
                RST_BC_DATA_WEST(KOUNT)=BND_VARS_V%VAR_3D(NV)%WEST(IB,JB,L,1)
              ENDDO
              ENDDO
            ENDDO
            DO L=1,LM
              DO JB=JTS,JTE
              DO IB=1,LNSV
                KOUNT=KOUNT+1
                RST_BC_DATA_WEST(KOUNT)=BND_VARS_V%VAR_3D(NV)%WEST(IB,JB,L,2)
              ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract Write Import State in SAVE_BC_DATA"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_StateGet(state      =EXP_STATE_SOLVER                 &  !<-- The Solver export state
                          ,itemName   ='Write Import State'             &  !<-- Name of the state to get from Solver export state
                          ,nestedState=IMP_STATE_WRITE                  &  !<-- Extract Write Component import state from Solver export
                          ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_SAVE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Set BC West Data Attribute in SAVE_BC_DATA"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

        CALL ESMF_AttributeSet(state    =IMP_STATE_WRITE                &  !<-- The Write component import state
                              ,name     ='RST_BC_DATA_WEST'             &  !<-- Name of 1-D string of west boundary values
                              ,itemCount=NUM_WORDS_BC_WEST              &  !<-- # of west boundary words on this fcst task
                              ,valueList=RST_BC_DATA_WEST               &  !<-- The 1-D data being inserted into the Write import state
                              ,rc       =RC)

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_SAVE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF
!
!-----------------------------------------------------------------------
!***  Eastern boundary data to 1-D
!-----------------------------------------------------------------------
!
      IF(ITE==IDE)THEN                                                     !<-- Tasks on east boundary
!
        KOUNT=0
!
        IF(NVARS_BC_2D_H>0)THEN
          DO NV=1,NVARS_BC_2D_H
            DO JB=JTS,JTE
            DO IB=1,LNSH
              KOUNT=KOUNT+1
              RST_BC_DATA_EAST(KOUNT)=BND_VARS_H%VAR_2D(NV)%EAST(IB,JB,1)
            ENDDO
            ENDDO
            DO JB=JTS,JTE
            DO IB=1,LNSH
              KOUNT=KOUNT+1
              RST_BC_DATA_EAST(KOUNT)=BND_VARS_H%VAR_2D(NV)%EAST(IB,JB,2)
            ENDDO
            ENDDO
          ENDDO
        ENDIF
!
        IF(NVARS_BC_3D_H>0)THEN
          DO NV=1,NVARS_BC_3D_H
            DO L=1,LM
              DO JB=JTS,JTE
              DO IB=1,LNSH
                KOUNT=KOUNT+1
                RST_BC_DATA_EAST(KOUNT)=BND_VARS_H%VAR_3D(NV)%EAST(IB,JB,L,1)
              ENDDO
              ENDDO
            ENDDO
            DO L=1,LM
              DO JB=JTS,JTE
              DO IB=1,LNSH
                KOUNT=KOUNT+1
                RST_BC_DATA_EAST(KOUNT)=BND_VARS_H%VAR_3D(NV)%EAST(IB,JB,L,2)
              ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF
!
        IF(NVARS_BC_4D_H>0)THEN
          DO NV=1,NVARS_BC_4D_H
            LBND=LBOUND(BND_VARS_H%VAR_4D(NV)%FULL_VAR,4)
            UBND=UBOUND(BND_VARS_H%VAR_4D(NV)%FULL_VAR,4)
            DO NL=LBND,UBND
              DO L=1,LM
                DO JB=JTS,JTE
                DO IB=1,LNSH
                  KOUNT=KOUNT+1
                  RST_BC_DATA_EAST(KOUNT)=BND_VARS_H%VAR_4D(NV)%EAST(IB,JB,L,1,NL)
                ENDDO
                ENDDO
              ENDDO
              DO L=1,LM
                DO JB=JTS,JTE
                DO IB=1,LNSH
                  KOUNT=KOUNT+1
                  RST_BC_DATA_EAST(KOUNT)=BND_VARS_H%VAR_4D(NV)%EAST(IB,JB,L,2,NL)
                ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF
!
        IF(NVARS_BC_2D_V>0)THEN
          DO NV=1,NVARS_BC_2D_V
            DO JB=JTS,JTE
            DO IB=1,LNSV
              KOUNT=KOUNT+1
              RST_BC_DATA_EAST(KOUNT)=BND_VARS_V%VAR_2D(NV)%EAST(IB,JB,1)
            ENDDO
            ENDDO
            DO JB=JTS,JTE
            DO IB=1,LNSV
              KOUNT=KOUNT+1
              RST_BC_DATA_EAST(KOUNT)=BND_VARS_V%VAR_2D(NV)%EAST(IB,JB,2)
            ENDDO
            ENDDO
          ENDDO
        ENDIF
!
        IF(NVARS_BC_3D_V>0)THEN
          DO NV=1,NVARS_BC_3D_V
            DO L=1,LM
              DO JB=JTS,JTE
              DO IB=1,LNSV
                KOUNT=KOUNT+1
                RST_BC_DATA_EAST(KOUNT)=BND_VARS_V%VAR_3D(NV)%EAST(IB,JB,L,1)
              ENDDO
              ENDDO
            ENDDO
            DO L=1,LM
              DO JB=JTS,JTE
              DO IB=1,LNSV
                KOUNT=KOUNT+1
                RST_BC_DATA_EAST(KOUNT)=BND_VARS_V%VAR_3D(NV)%EAST(IB,JB,L,2)
              ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract Write Import State in SAVE_BC_DATA"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_StateGet(state      =EXP_STATE_SOLVER                 &  !<-- The Solver export state
                          ,itemName   ='Write Import State'             &  !<-- Name of the state to get from Solver export state
                          ,nestedState=IMP_STATE_WRITE                  &  !<-- Extract Write Component import state from Solver export
                          ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_SAVE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Set BC East Data Attribute in SAVE_BC_DATA"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

        CALL ESMF_AttributeSet(state    =IMP_STATE_WRITE                &  !<-- The Write component import state
                              ,name     ='RST_BC_DATA_EAST'             &  !<-- Name of 1-D string of east boundary values
                              ,itemCount=NUM_WORDS_BC_EAST              &  !<-- # of east boundary words on this fcst task
                              ,valueList=RST_BC_DATA_EAST               &  !<-- The 1-D data being inserted into the Write import state
                              ,rc       =RC)

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_SAVE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE SAVE_BC_DATA
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE PHYSICS_INITIALIZE(GFS                                 &
                                   ,SHORTWAVE                           &
                                   ,LONGWAVE                            &
                                   ,CONVECTION                          &
                                   ,MICROPHYSICS                        &
                                   ,SFC_LAYER                           &
                                   ,TURBULENCE                          &
                                   ,LAND_SURFACE                        &
                                   ,CO2TF                               &
                                   ,NP3D                                &
                                   ,SBD,WBD                             &
                                   ,DPHD,DLMD                           &
                                   ,TPH0D,TLM0D                         &
                                   ,MY_DOMAIN_ID                        &
                                   ,MYPE                                &
                                   ,MPI_COMM_COMP                       &
                                   ,IDS,IDE,JDS,JDE,LM                  &
                                   ,IMS,IME,JMS,JME                     &
                                   ,ITS,ITE,JTS,JTE                     &
                                   ,RC)
!
!-----------------------------------------------------------------------
!
      USE MODULE_CONSTANTS,ONLY : A,CLIQ,CV,DTR,PI                      &
                                 ,RHOAIR0,RHOWATER,RHOSNOW
!
      USE MODULE_INIT_READ_BIN,ONLY : physics_read_gwd
!
!-----------------------------------------------------------------------
!***  Only for GFS physics
!-----------------------------------------------------------------------
!
      USE FUNCPHYS
      USE MODULE_MP_FER_HIRES, ONLY : GPVS_HR

      USE MERSENNE_TWISTER
      USE N_LAYOUT1,        ONLY : LATS_NODE_R,IPT_LATS_NODE_R
!      USE TRACER_CONST,     ONLY : SET_TRACER_CONST
!      USE DATE_DEF,         ONLY : FHOUR
      USE N_RESOL_DEF,      ONLY : LSOIL,LEVR,NXPT,JCAP,LEVS,NYPT       &
                                  ,JINTMX,THERMODYN_ID,SFCPRESS_ID      &
                                  ,NUM_P3D,NUM_P2D,NTOZ,NTCW,NCLD       &
                                  ,NMTVR,NFXR,LONR,LATR
!
      USE OZNE_DEF,         ONLY: LEVOZC,LATSOZP,BLATC,TIMEOZC,TIMEOZ   &
                                 ,KOZPL,LEVOZP,PL_TIME,PL_LAT,PL_PRES   &
                                 ,KOZC,DPHIOZC,LATSOZC,PL_COEFF
 
      USE N_NAMELIST_PHYSICS_DEF, ONLY: ISOL,ICO2,IALB,IEMS,IAER,ICTM   &
                                       ,IOVR_SW,IOVR_LW,LSSAV,LDIAG3D   &
                                       ,FHCYC,SASHAL,PRE_RAD,RAS,LSM    &
                                       ,CDMBGWD,DLQF,CTEI_RM,LGGFS3D    &
                                       ,BKGD_VDIF_M, SHAL_CNV           &
                                       ,BKGD_VDIF_H,BKGD_VDIF_S         &
                                       ,PSAUTCO,PRAUTCO,EVPCO           &
                                       ,CAL_PRE,MOM4ICE,MSTRAT          &
                                       ,TRANS_TRAC,NST_FCST             &
                                       ,MOIST_ADJ

      USE MODULE_CONTROL,ONLY : NMMB_FINALIZE

!
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument variables
!------------------------
!
      INTEGER(kind=KINT),INTENT(IN) :: CO2TF                            &
                                      ,MPI_COMM_COMP                    &
                                      ,MY_DOMAIN_ID                     &
                                      ,MYPE, NP3D
!
      INTEGER(kind=KINT),INTENT(IN) :: IDS,IDE,JDS,JDE,LM               &
                                      ,IMS,IME,JMS,JME                  &
                                      ,ITS,ITE,JTS,JTE
!
      REAL(kind=KFPT),INTENT(INOUT) :: DLMD,DPHD                        &
                                      ,TPH0D,TLM0D                      &
                                      ,SBD,WBD
!
      LOGICAL,INTENT(IN) :: GFS
!
      CHARACTER(99),INTENT(IN) :: CONVECTION,LONGWAVE,MICROPHYSICS      &
                                 ,SFC_LAYER,SHORTWAVE,TURBULENCE        &
                                 ,LAND_SURFACE
!
      INTEGER(kind=KINT),INTENT(OUT) :: RC
!
!---------------------
!***  Local variables
!---------------------
!
      INTEGER :: I,I_HI,I_LO,IHRST,IRTN,J,J_HI,J_LO,JULDAY,JULYR        &
                ,K,KFLIP,L,LPT2,N,NFCST,NRECS_SKIP_FOR_PT               &
                ,NSOIL,NSTEPS_PER_HOUR,NTIMESTEP
!
      INTEGER :: LDIM1,LDIM2,UDIM1,UDIM2
!
      INTEGER :: IAER_MDL, ISUBCSW, ISUBCLW, IFLIP,                     &
                 ICLIQ_SW, ICICE_SW, ICLIQ_LW, ICICE_LW
!
      INTEGER,DIMENSION(3) :: IDAT
!
      INTEGER,DIMENSION(:,:),ALLOCATABLE :: ITEMP,LOWLYR
!
      REAL :: SECOND_FCST
!
      REAL :: SWRAD_SCAT=1.
!
      REAL :: DELX,DELY,DPH,DT,DT_MICRO,DTPHS                           &
             ,GMT,JULIAN,PDBOT,PDTOP,PDTOT,PT_CB,RELM,RPDTOT            &
             ,SB,THETA_HALF,TPV,XTIME
!
      REAL,DIMENSION(LM) :: DSG1,PDSG1,PSGML1,SGML1,SGML2
      REAL,DIMENSION(LM+1) :: PSG1,SG1,SG2,SGM                          &
                             ,SFULL,SFULL_FLIP,SMID,SMID_FLIP
      REAL(KIND=KDBL),DIMENSION(LM+1) :: SFULLD
!
      REAL,DIMENSION(IMS:IME,JMS:JME) :: EMISS
      REAL,DIMENSION(:,:),ALLOCATABLE :: TEMP1,TEMP_GWD
      REAL,DIMENSION(:,:,:),ALLOCATABLE :: TEMPSOIL
      REAL,DIMENSION(NUM_SOIL_LAYERS)   :: SOIL1DIN
!
      CHARACTER(LEN=256) :: INFILE
!
      LOGICAL,SAVE :: ALLOWED_TO_READ=.TRUE.
      LOGICAL :: OPENED
      LOGICAL :: LSASHAL

      LOGICAL :: CRICK_PROOF, CCNORM, NORAD_PRECIP
!
#ifdef USE_GFS_PHYS
!---------------------------------
!***  GFS physics local variables
!---------------------------------
!
      CHARACTER(80)   :: GFS_PHY_NAMELIST
      INTEGER         :: JDAT(8),NLUNIT,NTRAC,IRET,IMJM,NIJ
      REAL(kind=KDBL) :: DELTIM,GAUL

      REAL *4              :: BLATC4
      REAL *4, ALLOCATABLE :: PL_LAT4(:), PL_PRES4(:), PL_TIME4(:), TEMPIN(:)

      REAL(kind=KDBL),DIMENSION(:),ALLOCATABLE ::                             &
                     SIG1T, RLA, RLO, SLMASK, OROG, AISFCS,                   &
                     SIHFCS, SICFCS, SITFCS, SWDFCS, VMNFCS, VMXFCS, SLPFCS,  &
                     ABSFCS, TSFFCS, SNOFCS, ZORFCS, TG3FCS, CNPFCS, SLIFCS,  &
                     F10MFCS, VEGFCS, VETFCS, SOTFCS, CVFCS, CVBFCS, CVTFCS
!
      REAL(kind=KDBL),DIMENSION(:,:),ALLOCATABLE :: ALFFC1              &
                                                   ,SMCFC1              &
                                                   ,STCFC1              &
                                                   ,SLCFC1              &
                                                   ,ALBFC1
!
#endif
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Initialize allocated arrays
!-----------------------------------------------------------------------
!
      NSOIL=NUM_SOIL_LAYERS                                              !<-- From Landsurface module
!
#ifdef USE_GFS_PHYS
      IF(GFS)THEN

        int_state%SOLCON=0.0D0
        int_state%SLAG  =0.0D0
        int_state%SDEC  =0.0D0
        int_state%CDEC  =0.0D0
!
        DO J=JTS,JTE
          int_state%DDY   (J)=0.0D0
          int_state%JINDX1(J)=0
          int_state%JINDX2(J)=0
        ENDDO
!
        DO J=JMS,JME
        DO I=IMS,IME
!
          int_state%DUGWD  (I,J)=0.0D0
          int_state%DVGWD  (I,J)=0.0D0
!
          int_state%TMPMIN (I,J)=373.0D0
          int_state%TMPMAX (I,J)=173.0D0
!
          int_state%SHDMIN (I,J)=0.0
          int_state%SHDMAX (I,J)=0.0
!
          int_state%SFALB  (I,J)=0.0D0
          int_state%TSFLW  (I,J)=0.0D0
          int_state%SEMIS  (I,J)=0.0D0
          int_state%SFCDLW (I,J)=0.0D0
          int_state%SFCDSW (I,J)=0.0D0
          int_state%SFCNSW (I,J)=0.0D0
          int_state%ZORFCS (I,J)=-1.D6
          int_state%SIHFCS (I,J)=-1.D6
          int_state%SICFCS (I,J)=-1.D6
          int_state%SLPFCS (I,J)=-1.D6
          int_state%TG3FCS (I,J)=-1.D6
          int_state%VEGFCS (I,J)=-1.D6
          int_state%VETFCS (I,J)=-1.D6
          int_state%SOTFCS (I,J)=-1.D6
!
          DO N=1,4
            int_state%ALBFC1(I,J,N)=-1.D6
          ENDDO
!
          DO N=1,2
            int_state%ALFFC1(I,J,N)=-1.D6
          ENDDO
!
          DO N=1,3                                 ! for Zhao =3, Ferr=1
            int_state%PHY_F2DV (I,J,N)=0.0D0
          ENDDO
!
          DO N=1,4                                 ! for Zhao =4, Ferr=3
          DO L=1,LM
            int_state%PHY_F3DV (I,J,L,N)=0.0D0
          ENDDO
          ENDDO
!
        ENDDO
        ENDDO

        rewind (kozpl)
        read   (kozpl) pl_coeff, latsozp, levozp, timeoz

        DO N=1,TIMEOZ
        DO L=1,PL_COEFF
        DO J=1,LEVOZP
        DO I=1,LATSOZP
          int_state%OZPLIN(I,J,L,N)=-1.D6
        ENDDO
        ENDDO
        ENDDO
        ENDDO
!
      ENDIF
!
#endif
!-----------------------------------------------------------------------
!***  Dereference the start time.
!-----------------------------------------------------------------------
!
      START_YEAR=int_state%START_YEAR
      START_MONTH=int_state%START_MONTH
      START_DAY=int_state%START_DAY
      START_HOUR=int_state%START_HOUR
      START_MINUTE=int_state%START_MINUTE
      START_SECOND=int_state%START_SECOND
      DT=int_state%DT
!
!-----------------------------------------------------------------------
!***  Radiation needs some specific time quantities.
!-----------------------------------------------------------------------
!
      CALL TIME_MEASURE(START_YEAR,START_MONTH,START_DAY,START_HOUR     &
                       ,START_MINUTE,START_SECOND                       &
                       ,NTIMESTEP,DT                                    &
                       ,JULDAY,JULYR,JULIAN,XTIME)
!
!-----------------------------------------------------------------------
!***  Open and read GWD data file (14 orography fields)
!-----------------------------------------------------------------------
!
      gwd_read: IF(int_state%GWDFLG) THEN
!
        select_GWD_unit: DO N=51,59
          INQUIRE(N,OPENED=OPENED)
          IF(.NOT.OPENED)THEN
            NFCST=N
            EXIT select_GWD_unit
          ENDIF
        ENDDO select_GWD_unit
!
        WRITE(INFILE,'(A,I2.2)')'GWD_bin_',MY_DOMAIN_ID
!
!-----------------------------------------------------------------------
!
        CALL PHYSICS_READ_GWD(INFILE,NFCST,INT_STATE                    &
                             ,MYPE,MPI_COMM_COMP                        &
                             ,IDS,IDE,JDS,JDE,RC)
!
        IF (RC /= 0) THEN
          RETURN
        ENDIF
!
!-----------------------------------------------------------------------
!
      ENDIF gwd_read
!
!-----------------------------------------------------------------------
!
      PT_CB=int_state%PT*1.0E-3   !<-- Convert pascals to centibars for GFDL initialization
!
!-----------------------------------------------------------------------
!***  Make up a potential skin temperature.
!-----------------------------------------------------------------------
!
      IF(.NOT.int_state%RESTART) THEN
!
        DO J=JTS,JTE
        DO I=ITS,ITE
          int_state%THS(I,J)=int_state%TSKIN(I,J)                       &
                       *(100000./(int_state%SG2(LM+1)*int_state%PD(I,J) &
                                 +int_state%PSG1(LM+1)))**CAPPA
        ENDDO
        ENDDO
!
      ENDIF
!
!-----------------------------------------------------------------------
!*** Initializing TLMAX, TLMIN
!-----------------------------------------------------------------------
!
      DO J=JTS,JTE
        DO I=ITS,ITE
          int_state%TLMAX(I,J)=int_state%T(I,J,1)
          int_state%TLMIN(I,J)=int_state%T(I,J,1)
       ENDDO
     ENDDO
!
!-----------------------------------------------------------------------
!***  Recreate sigma values at layer interfaces for the full vertical
!***  domain. 
!-----------------------------------------------------------------------
!
      DO L=1,LM+1
        SFULL(L)=int_state%SGM(L)
      ENDDO
!
      DO L=1,LM
        SMID(L)=(SFULL(L)+SFULL(L+1))*0.5
      ENDDO
!
      SMID(LM+1)=-9999999.
!
!-----------------------------------------------------------------------
!***  The radiative emissivity
!-----------------------------------------------------------------------
!
      DO J=JMS,JME
      DO I=IMS,IME
        EMISS(I,J)=1.
      ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
!***  Choose a J index for an "average" DX.
!***  Select the J that divides the domain's area in half.
!-----------------------------------------------------------------------
!
      SB=int_state%SBD*DTR
      DPH=int_state%DPHD*DTR
!!!   THETA_HALF=ASIN(0.5*SIN(-SB))
      THETA_HALF=0.
      JC=NINT(0.5*(JDE-JDS+1)+THETA_HALF/DPH)
!
!-----------------------------------------------------------------------
!***  Set time variables needed for history output.
!-----------------------------------------------------------------------
!
      NSTEPS_PER_HOUR=NINT(3600./int_state%DT)
      int_state%NPREC=NSTEPS_PER_HOUR*int_state%NHRS_PREC
      int_state%NCLOD=NSTEPS_PER_HOUR*int_state%NHRS_CLOD
      int_state%NHEAT=NSTEPS_PER_HOUR*int_state%NHRS_HEAT
      int_state%NRDLW=NSTEPS_PER_HOUR*int_state%NHRS_RDLW
      int_state%NRDSW=NSTEPS_PER_HOUR*int_state%NHRS_RDSW
      int_state%NSRFC=NSTEPS_PER_HOUR*int_state%NHRS_SRFC
!
!-----------------------------------------------------------------------
!***  If this is a restarted run from timestep 0 then zero out
!***  the accumulated precip since they pass through the analysis
!***  with nonzero values from the first guess.
!-----------------------------------------------------------------------
!
      IF(int_state%RST_OUT_00)THEN
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%ACPREC(I,J)=0.
          int_state%ACPREC_TOT(I,J)=0.
          int_state%CUPREC(I,J)=0.
        ENDDO
        ENDDO
      ENDIF
!
!-----------------------------------------------------------------------
!***  Finally initialize individual schemes.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  The GFS physics suite is considered a single package here.
!----------------------------------------------------------------------
!
      package: IF(GFS)THEN
!
#ifdef USE_GFS_PHYS
!-----------------------------------------------------------------------
!***  The GFS physics suite is considered a single package here.
!----------------------------------------------------------------------
!
      namelist_unit: DO N=101,151
        INQUIRE(N,OPENED=OPENED)
        IF(.NOT.OPENED)THEN
          NLUNIT=N
          EXIT namelist_unit
        ENDIF
      ENDDO namelist_unit

      GFS_PHY_NAMELIST = 'atm_namelist'
      DELTIM           = int_state%DT
      LEVS             = LM
      CALL N_COMPNS_PHYSICS(DELTIM,    IRET,                 &
                   NTRAC,   NXPT,    NYPT,  JINTMX,          &
                   JCAP,    LEVS,    LEVR,    LONR,   LATR,  &
                   NTOZ,    NTCW,    NCLD,                   &
                   LSOIL,   NMTVR,   NUM_P3D, NUM_P2D,       &
                   THERMODYN_ID,     SFCPRESS_ID,            &
                   NLUNIT,  MYPE,    GFS_PHY_NAMELIST)
!--------------
        IALB            = 0
        DELTIM          = int_state%DT
        LONR            = ITE-ITS+1  ! this is changed in compns_physics from
        LATR            = JTE-JTS+1  ! atm_namelist (restore it back)
        LATS_NODE_R     = JTE-JTS+1
        IPT_LATS_NODE_R = 1
        NFXR            = 39
        LSSAV           = .TRUE.   ! logical flag for store 3-d cloud field
        LDIAG3D         = .FALSE.  ! logical flag for store 3-d diagnostic fields
        LGGFS3D         = .FALSE.
!--------------
      CALL SET_SOILVEG(MYPE,NLUNIT)
      CALL SET_TRACER_CONST(NTRAC,MYPE,NLUNIT)
      CALL GFUNCPHYS
!
!-----------------------------------------------------------------------
!***  Initialize ozone
!-----------------------------------------------------------------------
     IF( NTOZ .LE. 0 ) THEN        ! Diagnostic ozone
!
!!        rewind (kozc)            ! there is no header in global_o3clim.txt file
!!        read   (kozc,end=101) latsozc, levozc, timeozc, blatc4
!!101     if (levozc .lt. 10 .or. levozc .gt. 100) then
!!          rewind (kozc)
            levozc  = 17
            latsozc = 18
            blatc   = -85.0
!!        else
!!          blatc   = blatc4
!!        endif
          latsozp   = 2
          levozp    = 1
          timeoz    = 1
          pl_coeff  = 1  !!!  0 (MUST set > 0, used in GBPHYS for allocation)
          timeozc   = 12 !!!  this is not in header
!
     ELSE                          ! Prognostic Ozone
!
         rewind (kozpl)
         read   (kozpl) pl_coeff, latsozp, levozp, timeoz
       IF(.NOT.ALLOCATED(pl_lat))THEN
         allocate (pl_lat (latsozp), pl_pres (levozp),pl_time (timeoz+1))
       ENDIF
       IF(.NOT.ALLOCATED(pl_lat4))THEN
         allocate (pl_lat4(latsozp), pl_pres4(levozp),pl_time4(timeoz+1))
       ENDIF
       IF(.NOT.ALLOCATED(tempin)) allocate (tempin(latsozp))
         rewind (kozpl)
         read (kozpl) pl_coeff, latsozp, levozp, timeoz, pl_lat4, pl_pres4, pl_time4
         pl_pres(:) = pl_pres4(:)
         pl_lat(:)  = pl_lat4(:)
         pl_time(:) = pl_time4(:)
!
         DO J=JTS,JTE

           gaul=int_state%GLAT( (ITS+ITE)/2 ,J)*180.0d0/3.14159d0

           int_state%jindx2(j) = latsozp + 1
           do i=1,latsozp
             if (gaul.lt. pl_lat(i)) then
               int_state%jindx2(j) = i
               exit
             endif
           enddo
           int_state%jindx1(j) = max(int_state%jindx2(j)-1,1)
           int_state%jindx2(j) = min(int_state%jindx2(j),latsozp)

           if (int_state%jindx2(j) .ne. int_state%jindx1(j)) then
             int_state%ddy(j) = (gaul                        - pl_lat(int_state%jindx1(j)))  &
                              / (pl_lat(int_state%jindx2(j)) - pl_lat(int_state%jindx1(j)))
           else
             int_state%ddy(j) = 1.0
           endif

         ENDDO
!
         DO I=1,TIMEOZ
           DO N=1,PL_COEFF
             DO K=1,LEVOZP
               READ(KOZPL) TEMPIN
               int_state%OZPLIN(:,K,N,I) = TEMPIN(:)
             ENDDO
           ENDDO
         ENDDO
!
     ENDIF                          ! Diagnostic/Prognostic Ozone
!
        dphiozc = -(blatc+blatc)/(latsozc-1)
!-----------------------------------------------------------------------
!***  End initialization  of ozone
!-----------------------------------------------------------------------
!
       IMJM=(ITE-ITS+1)*(JTE-JTS+1)
       ALLOCATE(SIG1T(IMJM),RLA(IMJM),RLO(IMJM),SLMASK(IMJM),OROG(IMJM)   &
               ,AISFCS(IMJM),SIHFCS(IMJM),SICFCS(IMJM),SITFCS(IMJM)       &
               ,SWDFCS(IMJM),VMNFCS(IMJM),VMXFCS(IMJM),SLPFCS(IMJM)       &
               ,ABSFCS(IMJM),TSFFCS(IMJM),SNOFCS(IMJM),ZORFCS(IMJM)       &
               ,TG3FCS(IMJM),CNPFCS(IMJM),SLIFCS(IMJM),F10MFCS(IMJM)      &
               ,VEGFCS(IMJM),VETFCS(IMJM),SOTFCS(IMJM),CVFCS(IMJM)        &
               ,CVBFCS(IMJM),CVTFCS(IMJM),ALFFC1(IMJM,2),ALBFC1(IMJM,4)   &
               ,SMCFC1(IMJM,NSOIL),STCFC1(IMJM,NSOIL),SLCFC1(IMJM,NSOIL) )

       SIHFCS  = 0.0d0
       SICFCS  = 0.0d0
       SITFCS  = 0.0d0
       SWDFCS  = 0.0d0
       VMNFCS  = 0.0d0
       VMXFCS  = 0.0d0
       SLPFCS  = 0.0d0
       ABSFCS  = 0.0d0
       TSFFCS  = 0.0d0
       SNOFCS  = 0.0d0
       ZORFCS  = 0.0d0
       TG3FCS  = 0.0d0
       CNPFCS  = 0.0d0
       SLIFCS  = 0.0d0
       F10MFCS = 0.0d0
       VEGFCS  = 0.0d0
       VETFCS  = 0.0d0
       SOTFCS  = 0.0d0
       CVFCS   = 0.0d0
       CVBFCS  = 0.0d0
       CVTFCS  = 0.0d0
       ALFFC1  = 0.0d0
       ALBFC1  = 0.0d0
       SMCFC1  = 0.0d0
       STCFC1  = 0.0d0
       SLCFC1  = 0.0d0


       JDAT(1)=int_state%IDAT(3)
       JDAT(2)=int_state%IDAT(2)
       JDAT(3)=int_state%IDAT(1)
       JDAT(4)=0
       JDAT(5)=int_state%IHRST
       JDAT(6)=0
       JDAT(7)=0
       JDAT(8)=0
       FHOUR=FLOAT(int_state%IHRST)
!
    NIJ=0
    DO J=JTS,JTE
      DO I=ITS,ITE

       NIJ=NIJ+1

       SIG1T(NIJ)    = 0.d0
       RLA(NIJ)      = int_state%GLAT(I,J)*180.0d0/3.14159d0
       RLO(NIJ)      = int_state%GLON(I,J)*180.0d0/3.14159d0
       SLMASK(NIJ)   = 1.0d0-int_state%SM(I,J)
       OROG(NIJ)     = int_state%FIS(I,J)/9.81d0
       AISFCS(NIJ)   = int_state%SICE(I,J)

      ENDDO
    ENDDO
!
      CALL SFCCYCLE(204,(ITE-ITS+1)*(JTE-JTS+1),4,SIG1T,FHCYC              &
     &,             JDAT(1), JDAT(2), JDAT(3), JDAT(5), FHOUR              &
!    &,             RLA, RLO, SLMASK, OROG                                 &
     &,             RLA, RLO, SLMASK, OROG, OROG, .FALSE.                  &
     &,             SIHFCS,   SICFCS, SITFCS                               &
     &,             SWDFCS,   SLCFC1                                       &
     &,             VMNFCS,   VMXFCS, SLPFCS, ABSFCS                       &
     &,             TSFFCS,   SNOFCS, ZORFCS, ALBFC1, TG3FCS               &
     &,             CNPFCS,   SMCFC1, STCFC1, SLIFCS, AISFCS, F10MFCS      &
     &,             VEGFCS,   VETFCS, SOTFCS, ALFFC1                       &
     &,             CVFCS,    CVBFCS, CVTFCS, MYPE, NLUNIT, IALB)
!
     NIJ=0
     DO J=JTS,JTE
       DO I=ITS,ITE
        NIJ=NIJ+1

        SIHFCS(NIJ) = int_state%SICE(I,J) * 1.0d0  ! initialize like this
        SICFCS(NIJ) = int_state%SICE(I,J) * 0.9d0  ! initialize like this

        if(int_state%SICE(I,J) > 0.5 ) then
          SLPFCS(NIJ)=9.0d0
          VEGFCS(NIJ)=0.01d0
          SOTFCS(NIJ)=9.0d0
          VETFCS(NIJ)=13.0d0
        endif

        int_state%ZORFCS(I,J)   = ZORFCS(NIJ)
        int_state%SIHFCS(I,J)   = SIHFCS(NIJ)
        int_state%SICFCS(I,J)   = SICFCS(NIJ)
        int_state%SLPFCS(I,J)   = SLPFCS(NIJ)
        int_state%TG3FCS(I,J)   = TG3FCS(NIJ)
        int_state%VEGFCS(I,J)   = VEGFCS(NIJ)
        int_state%VETFCS(I,J)   = VETFCS(NIJ)
        int_state%SOTFCS(I,J)   = SOTFCS(NIJ)
!!!
        int_state%isltyp(I,J)   = nint(SOTFCS(NIJ))
        int_state%ivgtyp(I,J)   = nint(VETFCS(NIJ))
!!!

        int_state%ALBFC1(I,J,1) = ALBFC1(NIJ,1)
        int_state%ALBFC1(I,J,2) = ALBFC1(NIJ,2)
        int_state%ALBFC1(I,J,3) = ALBFC1(NIJ,3)
        int_state%ALBFC1(I,J,4) = ALBFC1(NIJ,4)

        int_state%ALFFC1(I,J,1) = ALFFC1(NIJ,1)
        int_state%ALFFC1(I,J,2) = ALFFC1(NIJ,2)

       ENDDO
     ENDDO

     DEALLOCATE(SIG1T,RLA,RLO,SLMASK,OROG     &
               ,AISFCS,SIHFCS,SICFCS,SITFCS   &
               ,SWDFCS,VMNFCS,VMXFCS,SLPFCS   &
               ,ABSFCS,TSFFCS,SNOFCS,ZORFCS   &
               ,TG3FCS,CNPFCS,SLIFCS,F10MFCS  &
               ,VEGFCS,VETFCS,SOTFCS,CVFCS    &
               ,CVBFCS,CVTFCS,ALFFC1,ALBFC1   &
               ,SMCFC1,STCFC1,SLCFC1)
!
!----------------------------------------------------------------------
!***  Set fluxes to zero
!----------------------------------------------------------------------
!
             int_state%ALWIN(I,J)         = 0.
             int_state%ALWOUT(I,J)        = 0.
             int_state%ASWIN(I,J)         = 0.
             int_state%ASWOUT(I,J)        = 0.
             int_state%RLWIN(I,J)         = 0.
             int_state%RADOT(I,J)         = 0.
             int_state%RSWIN(I,J)         = 0.
             int_state%RSWOUT(I,J)        = 0.

             int_state%ALWTOA(I,J)        = 0.
             int_state%ASWTOA(I,J)        = 0.
             int_state%RLWTOA(I,J)        = 0.
             int_state%RSWTOA(I,J)        = 0.

             int_state%SFCSHX(I,J)        = 0.
             int_state%SFCLHX(I,J)        = 0.
             int_state%TWBS(I,J)          = 0.
             int_state%QWBS(I,J)          = 0.

             int_state%BGROFF(I,J)        = 0.
             int_state%SSROFF(I,J)        = 0.
             int_state%ACSNOW(I,J)        = 0.

             int_state%CUPPT(I,J)         = 0.
!
!----------------------------------------------------------------------
!***  End GFS package init
!----------------------------------------------------------------------
#endif
!
      ELSE
!
!----------------------------------------------------------------------
!***  If not selecting the GFS suite, each of the physics groups is
!***  treated individually.
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
!***  Longwave radiation
!----------------------------------------------------------------------
!
        SELECT CASE (longwave)  
          CASE ('gfdl')
!
!***  We are calling a WRF routine thus flip the vertical.
!
            DO K=1,LM
              KFLIP=LM+1-K
              SFULL_FLIP(KFLIP)=SFULL(K+1)
              SMID_FLIP(KFLIP)=SMID(K)
            ENDDO
            SFULL_FLIP(LM+1)=SFULL(1)
!
            GMT=REAL(START_HOUR)
            CALL GFDL_INIT(EMISS,SFULL_FLIP,SMID_FLIP,PT_CB            &
                          ,JULYR,START_MONTH,START_DAY,GMT             &
                          ,CO2TF                                       &
                          ,IDS,IDE,JDS,JDE,1,LM+1                      &
                          ,IMS,IME,JMS,JME,1,LM+1                      &
                          ,ITS,ITE,JTS,JTE,1,LM)
          CASE ('rrtm')

            CALL GPKAP    ! for ozone by using the unified RRTM from GFS
            CALL GPVS     ! for aerosol by using the unified RRTM from GFS

            CALL GPVS_HR  !- Initialize regional version of FPVS, FPVS0 functions
!
!-----------------------------------------------------------------------
!***  For threading safe  (rad_initialize). Default value
!-----------------------------------------------------------------------
!
            ICTM=1     !  0: use data at initial cond time, if not available, use latest, no extrapolation.
                       !  1: use data at the forecast time, if not available, use latest and extrapolation.
                       ! -1: use user provided external data for the fcst time, no extrapolation.
                       ! -2: same as ictm=0, but add seasonal cycle from climatology. no extrapolation.
                       ! yyyy0: use yyyy data for the forecast time, no further data extrapolation.
                       ! yyyy1: use yyyy data for the fcst. if needed, do extrapolation to match the fcst time.
            ISOL=0     ! 0: use a fixed solar constant value 1.3660e+3 (default)
                       !10: use a fixed solar constant value 1.3608e+3
                       ! 1: use 11-year cycle solar constant table
            ICO2=1     ! 0: use prescribed global mean co2   (default)
                       ! 1: use observed co2 annual mean value only
                       ! 2: use obs co2 monthly data with 2-d variation
            IAER=11    ! flag for aerosols scheme selection (all options work for NMMB)
                       ! - 3-digit aerosol flag (volc,lw,sw)
                       !   0: turn all aeros effects off (sw,lw,volc)
                       !   1: use clim tropspheric aerosol for sw only
                       !  10: use clim tropspheric aerosol for lw only
                       !  11: use clim tropspheric aerosol for both sw and lw
                       ! 100: volc aerosol only for both sw and lw
                       ! 101: volc and clim trops aerosol for sw only
                       ! 110: volc and clim trops aerosol for lw only
                       ! 111: volc and clim trops aerosol for both sw and lw
                       !   2: gocart/BSC-Dust tropspheric aerosol for sw only
                       !  20: gocart/BSC-Dust tropspheric aerosol for lw only
                       !  22: gocart/BSC-Dust tropspheric aerosol for both sw and lw
                       ! 102: volc and gocart trops aerosol for sw only
                       ! 120: volc and gocart trops aerosol for lw only
                       ! 122: volc and gocart trops aerosol for both sw and lw
            IAER_MDL=0 !  default aerosol model is opac-climatology
                       !  > 0,  future gocart-clim/prog scheme (not ready)
            IALB=2     ! control flag for surface albedo schemes
                       ! 0: climatology, based on surface veg types  ! ONLY THIS ONE WORKS (GFS)
                       ! 1: modis retrieval based surface albedo scheme
                       ! 2: use externally provided albedoes directly. ! ONLY THIS ONE WORKS for regional
                       !    (CALCULATES ALBEDO FROM NMMB MONTHLY CLIMATOLOGY AS IN GFDL RADIATION)
            IEMS=0     ! control flag for surface emissivity schemes
                       ! 0: fixed value of 1.0   (default)
                       ! 1: varying value based on surface veg types
            NTCW=3     !  0: no cloud condensate calculated
                       ! >0: array index location for cloud condensate
          ! NP3D=3     ! 3: ferrier's microphysics cloud scheme (only stratiform cloud)
                       !    (set iflagliq>0 in radsw_param.f and radlw_param.f)
                       ! 4: zhao/carr/sundqvist microphysics cloud (now available in the NMMB)
                       ! 5: NAM stratiform + convective cloud optical depth and fraction
                       !    (set iflagliq=0 in radsw_param.f and radlw_param.f)
            NTOZ=0     !  0: climatological ozone profile
                       ! >0: interactive ozone profile
            IOVR_SW=1  !  0 sw: random overlap clouds
                       !  1 sw: max-random overlap clouds
            IOVR_LW=1  !  0 lw: random overlap clouds
                       !  1 lw: max-random overlap clouds
            ISUBCSW=0  !  isubcsw/isubclw
                       !  sub-column cloud approx control flag (sw/lw rad)
                       !  0: with out sub-column cloud approximation
                       !  1: mcica sub-col approx. prescribed random seed
                       !  2: mcica sub-col approx. provided random seed
            ISUBCLW=0

            !----------------------------------------------------------
            ! --- check physparam for detail of the following ---------

            ICLIQ_SW=1 ! sw optical property for liquid clouds
            ICICE_SW=3 ! sw optical property for ice clouds (only iswcliq>0)
            ICLIQ_LW=1 ! lw optical property for liquid clouds
            ICICE_LW=1 ! lw optical property for ice clouds (only ilwcliq>0)

            !----------------------------------------------------------

            IFLIP=0    !   0: input data from toa to sfc
                       !   1: input data from sfc to toa

            SASHAL=0              ! New Massflux based shallow convection  (Not in use for NMMB)
            LSASHAL=.false.
            if (SASHAL>0 .and. .not.RAS) LSASHAL=.true.
            CRICK_PROOF=.false.   ! flag for eliminating CRICK (smooths profiles)
            CCNORM=.true.         ! flag for incloud condensate mixing ratio
            NORAD_PRECIP=.false.  ! flag for precip in radiation
                                  ! .true. snow/rain has no impact on radiation

!-----------------------------------------------------------------------
!***  Initialize ozone
!-----------------------------------------------------------------------

!OZONE CLIMATOLOGY
!
! there is no header in global_o3clim.txt file

            IF (NTOZ .LE. 0) THEN     ! DIAGNOSTIC OZONE, ONLY THIS ONE WORKS
               LEVOZC  = 17
               LATSOZC = 18
               BLATC   = -85.0
               TIMEOZC = 12            !!!  this is not in header
               LATSOZP   = 2
               LEVOZP    = 1
               TIMEOZ    = 1
               PL_COEFF  = 0
            ENDIF

            DPHIOZC = -(BLATC+BLATC)/(LATSOZC-1)

!-----------------------------------------------------------------------
!***  End initialization  of ozone
!-----------------------------------------------------------------------

            DO L=1,LM+1
              SFULLD(L)=SFULL(L)    !-- double precision
            ENDDO

!==========================================================================
!  Similar to GFS "GFS_Initialize_ESMFMod.f" line #1103
!==========================================================================

!..Special case for altering microphysics coupling with RRTM radiation
!.. based on namelist settings.  The NP3Dx variable is incredibly convoluted
!.. and renamed many times, including icmphys, np3d, and num_p3d.  Extremely
!.. confusing and hard-wired and needs help to adapt to new physics couplings
!.. and choices for full flexibility.   G. Thompson 06Feb2013

!..SPECIAL TEST FOR THOMPSON MICROPHYSICS AND RRTM RADIATION.  It is strongly
!.. advised against using GFDL or other radiation in combination with Thompson
!.. microphysics because other schemes are not properly using the cloud data.

            IF (TRIM(int_state%SHORTWAVE)=='rrtm' .AND.                 &
     &          TRIM(int_state%MICROPHYSICS)=='thompson' ) THEN

              IF (NP3D /=8) THEN
                 WRITE(0,*)' User selected np3d=',NP3D
                 WRITE(0,*)' NP3D=8 for RRTM & THOMPSON MICROPHYSICS'
                 CALL NMMB_FINALIZE
              ENDIF

              ICICE_SW=4
              ICICE_LW=4

            ENDIF

!==========================================================================
!..For GFDL type diagnostic
!==========================================================================

            IF (NP3D == 5) THEN
              ICLIQ_SW=0
              ICLIQ_LW=0
            ENDIF

            IF(MYPE==0)THEN
              WRITE(0,*)' Model Proces np3d=',NP3D
            ENDIF

!==========================================================================

            call rad_initialize_nmmb                                   &
!        ---  inputs:
     &       ( SFULLD,LM,ICTM,ISOL,ICO2,IAER,IAER_MDL,IALB,IEMS,NTCW,  &
     &         NP3D,NTOZ,IOVR_SW,IOVR_LW,ISUBCSW,ISUBCLW,              &
     &         ICLIQ_SW,ICICE_SW,ICLIQ_LW,ICICE_LW,                    &
     &         LSASHAL,CRICK_PROOF,CCNORM,NORAD_PRECIP,IFLIP,MYPE )
!  ---        outputs:
!                ( none )

!==========================================================================
!==========================================================================

            DO K=1,LM
              KFLIP=LM+1-K
              SFULL_FLIP(KFLIP)=SFULL(K+1)
              SMID_FLIP(KFLIP)=SMID(K)
            ENDDO
            SFULL_FLIP(LM+1)=SFULL(1)
!
            GMT=REAL(START_HOUR)

!==========================================================================
! This following "RRTM_INIT" is only a L,M,H  DIAGNOSTIC cloud.
! It is not a real RRTM initialization
!==========================================================================


            CALL RRTM_INIT(EMISS,SFULL_FLIP,SMID_FLIP,PT_CB            &
                          ,JULYR,START_MONTH,START_DAY,GMT             &
                          ,CO2TF                                       &
                          ,IDS,IDE,JDS,JDE,1,LM+1                      &
                          ,IMS,IME,JMS,JME,1,LM+1                      &
                          ,ITS,ITE,JTS,JTE,1,LM)
!
          CASE DEFAULT
            WRITE(0,*)' BAD SELECTION OF LONGWAVE SCHEME: INIT '
        END SELECT
!
!----------------------------------------------------------------------
!***  Shortwave radiation
!----------------------------------------------------------------------
!
        SELECT CASE (shortwave)
          CASE ('gfdl')
!           WRITE(0,*)' Already called GFDL_INIT from LONGWAVE'
          CASE ('rrtm')
!           WRITE(0,*)' Already called RRTM_INIT from LONGWAVE'
!!!       CASE ('gsfc')
!!!         CALL GSFC_INIT
          CASE ('dudh')
!!!         CALL SWINIT(SWRAD_SCAT,int_state%RESTART                   &
!!!                    ,ALLOWED_TO_READ                                &
!!!                    ,IDS,IDE,JDS,JDE,1,LM+1                         &
!!!                    ,IMS,IME,JMS,JME,1,LM+1                         &
!!!                    ,ITS,ITE,JTS,JTE,1,LM)
          CASE DEFAULT
            WRITE(0,*)' BAD SELECTION OF SHORTWAVE SCHEME: INIT'
        END SELECT
!
!----------------------------------------------------------------------
!***  Surface layer
!----------------------------------------------------------------------
!
        ALLOCATE(LOWLYR(IMS:IME,JMS:JME),STAT=I)
!
        SELECT CASE (sfc_layer)
          CASE ('myj')
            CALL JSFC_INIT(LOWLYR                                      &  !<-- Placeholder (computed in TURBULENCE)
                          ,int_state%USTAR,int_state%Z0                &
                          ,int_state%SM,int_state%SICE                 &
                          ,int_state%IVGTYP,int_state%RESTART          &            
                          ,ALLOWED_TO_READ                             &
                          ,IDS,IDE,JDS,JDE,1,LM+1                      &
                          ,IMS,IME,JMS,JME,1,LM+1                      &
                          ,ITS,ITE,JTS,JTE,1,LM                        &
                          ,MPI_COMM_COMP )
          CASE ('gfdl')
            CALL JSFC_INIT4GFDL(LOWLYR                                      &  !<-- Placeholder (computed in TURBULENCE)
                          ,int_state%USTAR,int_state%Z0                &
                          ,int_state%SM,int_state%SICE                 &
                          ,int_state%IVGTYP,int_state%RESTART          &            
                          ,ALLOWED_TO_READ                             &
                          ,IDS,IDE,JDS,JDE,1,LM+1                      &
                          ,IMS,IME,JMS,JME,1,LM+1                      &
                          ,ITS,ITE,JTS,JTE,1,LM                        &
                          ,MPI_COMM_COMP )
!!!       CASE ('mm5')
!!!         CALL SFCLYR_INIT
          CASE DEFAULT
            WRITE(0,*)' BAD SELECTION OF SURFACE LAYER SCHEME: INIT'
        END SELECT
!
!----------------------------------------------------------------------
!***  Turbulence
!----------------------------------------------------------------------
!
        SELECT CASE (turbulence)
          CASE ('myj')
            CALL MYJPBL_INIT(int_state%EXCH_H,int_state%RESTART        &
                            ,IDS,IDE,JDS,JDE,LM                        &
                            ,IMS,IME,JMS,JME                           &
                            ,ITS,ITE,JTS,JTE)
          CASE ('gfs')
!!!       CASE ('ysu')
!!!         CALL YSU_INIT
          CASE ('gfshur')
          CASE ('gfsedmfhur')
          CASE DEFAULT
            WRITE(0,*)' BAD SELECTION OF TURBULENCE SCHEME: INIT'
        END SELECT
!
!----------------------------------------------------------------------
!***  Land surface
!----------------------------------------------------------------------
!
        SELECT CASE (land_surface)
          CASE ('noah')
          int_state%LSM_PHYSICS=LSMSCHEME
          CALL NOAH_LSM_INIT(int_state%CMC,     int_state%ISLTYP       &
                            ,int_state%STC,     int_state%SMC          &
                            ,int_state%IVEGSRC                         &
                            ,int_state%SH2O,    NUM_SOIL_LAYERS        &
                            ,int_state%RESTART, ALLOWED_TO_READ        &
                            ,IDS,IDE, JDS,JDE                          &
                            ,IMS,IME, JMS,JME                          &
                            ,ITS,ITE, JTS,JTE                          &
                            ,MYPE,MPI_COMM_COMP )
          CASE ('liss')
          int_state%LSM_PHYSICS=LISSSCHEME

          CASE ('gfdlslab')
          int_state%LSM_PHYSICS=GFDLSLABSCHEME
!            WRITE(0,*)'See GFDL Surface Layer SF_GFDL'

          CASE DEFAULT
            WRITE(0,*)' BAD SELECTION OF LAND SURFACE SCHEME: INIT'
        END SELECT
!
!----------------------------------------------------------------------
!****  Convection
!----------------------------------------------------------------------
!
        SELECT CASE (convection)
          CASE ('bmj')
            int_state%CU_PHYSICS=BMJSCHEME
            CALL BMJ_INIT(int_state%CLDEFI,int_state%RESTART &
                         ,a2,a3,a4,cappa,cp &
                         ,pq0,r_d &
                         ,IDS,IDE,JDS,JDE &
                         ,IMS,IME,JMS,JME &
                         ,ITS,ITE,JTS,JTE,LM)

          CASE ('sas')
            int_state%CU_PHYSICS=SASSCHEME
            CALL SAS_INIT
!
          CASE ('sashur')
            int_state%CU_PHYSICS=SASHURSCHEME
            CALL SASHUR_INIT
!
          CASE ('scalecu')
            int_state%CU_PHYSICS=SCALECUSCHEME
            CALL SCALECU_INIT( IMS,IME,JMS,JME &
                              ,ITS,ITE,JTS,JTE,lm &
                              ,int_state%DUDT,int_state%DVDT &
                               )
!
          CASE ('none')
!           WRITE(0,*)' User has chosen to run with no parameterized convection.'
          CASE DEFAULT
            WRITE(0,*)' BAD SELECTION OF CONVECTION SCHEME: INIT'
            WRITE(0,*)' User selected CONVECTION = ',TRIM(CONVECTION)
            CALL NMMB_FINALIZE
        END SELECT
!
!----------------------------------------------------------------------
!***  Microphysics
!----------------------------------------------------------------------
!
        SELECT CASE (microphysics)
!
          CASE ('fer')
            int_state%MP_PHYSICS=95
            DT_MICRO=int_state%NPRECIP*DT
            DELX=-2.*int_state%WBD*111.3/REAL(int_state%IM) !DX at rotated equator (km)
            DELY=-2.*int_state%SBD*111.3/REAL(int_state%JM) !DY at rotated equator (km)
!
            CALL FERRIER_INIT(DT_MICRO,DT,DELX,DELY,int_state%RESTART  &
                             ,int_state%F_ICE                          &
                             ,int_state%F_RAIN                         &
                             ,int_state%F_RIMEF                        &
                             ,int_state%MP_RESTART_STATE               &
                             ,int_state%TBPVS_STATE                    &
                             ,int_state%TBPVS0_STATE                   &
                             ,ALLOWED_TO_READ                          &
                             ,IDS,IDE,JDS,JDE,1,LM+1                   &
                             ,IMS,IME,JMS,JME,1,LM                     &
                             ,ITS,ITE,JTS,JTE,1,LM                     &
                             ,MPI_COMM_COMP,MYPE,int_state%MASSRout    &
                             ,int_state%MASSIout)
!
          CASE ('fer_hires')
            int_state%MP_PHYSICS=5
            DT_MICRO=int_state%NPRECIP*DT
            DELX=-2.*int_state%WBD*111.3/REAL(int_state%IM) !DX at rotated equator (km)
            DELY=-2.*int_state%SBD*111.3/REAL(int_state%JM) !DY at rotated equator (km)
!
            CALL FERRIER_INIT_HR(DT_MICRO,DT,DELX,DELY,int_state%RESTART  &
                                ,int_state%F_ICE                          &
                                ,int_state%F_RAIN                         &
                                ,int_state%F_RIMEF                        &
                                ,int_state%MP_RESTART_STATE               &
                                ,int_state%TBPVS_STATE                    &
                                ,int_state%TBPVS0_STATE                   &
                                ,ALLOWED_TO_READ                          &
                                ,IDS,IDE,JDS,JDE,1,LM+1                   &
                                ,IMS,IME,JMS,JME,1,LM                     &
                                ,ITS,ITE,JTS,JTE,1,LM                     &
                                ,MPI_COMM_COMP,MYPE,int_state%MASSRout    &
                                ,int_state%MASSIout)
!
          CASE ('gfs')
            int_state%MP_PHYSICS=99
            CALL GFSMP_INIT
!
          CASE ('wsm6')
            int_state%MP_PHYSICS=6
            CALL WSM6INIT(RHOAIR0,RHOWATER,RHOSNOW,CLIQ,CV                &
                         ,ALLOWED_TO_READ )
!
          CASE ('thompson')
            int_state%MP_PHYSICS=8
            CALL thompson_init()
!
          CASE DEFAULT
            WRITE(0,*)' BAD SELECTION OF MICROPHYSICS SCHEME: INIT'
            WRITE(0,*)' User selected MICROPHYSICS = ',TRIM(MICROPHYSICS)
            CALL NMMB_FINALIZE

        END SELECT
!
!----------------------------------------------------------------------
!****  Gravity wave drag (GWD) & mountain blocking (MB) initialization
!----------------------------------------------------------------------
!
        DTPHS=int_state%DT*int_state%NPHS
!
        CALL GWD_init(DTPHS,int_state%RESTART                           &
                      ,int_state%CLEFFAMP,int_state%DPHD                &
                       ,int_state%CLEFF                                 &
                      ,int_state%TPH0D,int_state%TLM0D                  &
                      ,int_state%GLAT,int_state%GLON                    &
                      ,int_state%CROT,int_state%SROT,int_state%HANGL    &
                      ,IDS,IDE,JDS,JDE                                  &
                      ,IMS,IME,JMS,JME                                  &
                      ,ITS,ITE,JTS,JTE,LM)
!
! uncomment this for output in future
!
!       IF(.NOT.int_state%RESTART)THEN
!         DO J=JMS,JME
!         DO I=IMS,IME
!           UGWDsfc(I,J)=0.
!           VGWDsfc(I,J)=0.
!         ENDDO
!         ENDDO
!       ENDIF
!
!
!----------------------------------------------------------------------
!
      ENDIF package
!
!----------------------------------------------------------------------
!
      END SUBROUTINE PHYSICS_INITIALIZE
!

!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE UPDATE_WATER(CWM,F_ICE,F_RAIN,F_RIMEF                  &
                             ,T,QC,QR,QS,QI,QG                          &
                             ,MICROPHYSICS,SPEC_ADV,NTIMESTEP           &
                             ,IDS,IDE,JDS,JDE,LM                        &
                             ,IMS,IME,JMS,JME                           &
                             ,ITS,ITE,JTS,JTE)
!***********************************************************************
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .
! SUBPROGRAM:    UPDATE_WATER          UPDATE WATER ARRAY
!   PRGRMMR: FERRIER         ORG: NP22     DATE: 3 AUG 2009
!
! ABSTRACT:
!     UPDATE WATER ARRAY FOR FERRIER MICROPHYSICS
!
! PROGRAM HISTORY LOG (with changes to called routines) :
!   2009-08     FERRIER     - Synchronize WATER array with CWM, F_rain, F_ice arrays
!
! USAGE: CALL UPDATE_WATER FROM PHY_RUN
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!-----------------------------------------------------------------------
      USE MODULE_CONSTANTS,ONLY : EPSQ,TIW
!-----------------------------------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------------------------------
!
!----------------------
!-- Argument Variables
!----------------------
!
      INTEGER,INTENT(IN) :: NTIMESTEP                                   &
                           ,IDS,IDE,JDS,JDE,LM                          &
                           ,IMS,IME,JMS,JME                             &
                           ,ITS,ITE,JTS,JTE
!
      CHARACTER(99),INTENT(IN) :: MICROPHYSICS
!
      LOGICAL,INTENT(IN) :: SPEC_ADV
!
      REAL,DIMENSION(IMS:IME,JMS:JME,1:LM),INTENT(INOUT) :: CWM         &
                                                           ,F_ICE       &
                                                           ,F_RAIN      &
                                                           ,F_RIMEF     &
                                                           ,T
!
      REAL,DIMENSION(:,:,:),POINTER,INTENT(INOUT) :: QC,QR,QS,QI,QG
!
!--------------------
!--  Local Variables
!--------------------
!
      INTEGER :: I,J,K, NW
      REAL :: FRACTION, LIQW, OLDCWM
      LOGICAL :: CLD_INIT
      LOGICAL :: deep_ice
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      IF(NTIMESTEP<=1)THEN
        CLD_INIT=.TRUE.
      ELSE
        CLD_INIT=.FALSE.
      ENDIF
!
!----------------------------------------------------------------------
!-- Couple 2 sets of condensed water arrays for different microphysics: 
!   QC,QR,QS, etc. arrays <=> CWM,F_ice,F_rain,F_RimeF 3D arrays
!----------------------------------------------------------------------
!
      SELECT CASE ( TRIM(MICROPHYSICS) )
!
!----------------------------------------------------------------------
        CASE ('fer','fer_hires')  !-- Update fields for Ferrier microphysics
!----------------------------------------------------------------------
!
          spec_adv_fer: IF (.NOT.SPEC_ADV .OR. CLD_INIT) THEN
!-- Update WATER arrays when advecting only total condensate (spec_adv=F)
!   or at the initial time step
            DO K=1,LM
             DO J=JMS,JME
              DO I=IMS,IME
                IF (CWM(I,J,K)>EPSQ) THEN
                  LIQW=(1.-F_ice(I,J,K))*CWM(I,J,K)
                  QC(I,J,K)=(1.-F_rain(I,J,K))*LIQW
                  QR(I,J,K)=F_rain(I,J,K)*LIQW
                  QS(I,J,K)=F_ice(I,J,K)*CWM(I,J,K)
                ELSE
                  QC(I,J,K)=0.
                  QR(I,J,K)=0.
                  QS(I,J,K)=0.
                ENDIF
              ENDDO
             ENDDO
            ENDDO
!
          ELSE spec_adv_fer
!-- Update CWM,F_ICE,F_RAIN arrays from separate species advection (spec_adv=T)
            DO K=1,LM
             DO J=JMS,JME
              DO I=IMS,IME
                CWM(I,J,K)=QC(I,J,K)+QR(I,J,K)+QS(I,J,K)
                IF (QS(I,J,K)>EPSQ) THEN
                  F_ICE(I,J,K)=QS(I,J,K)/CWM(I,J,K)
                ELSE
                  F_ICE(I,J,K)=0.0
                ENDIF
                IF (QR(I,J,K)>EPSQ) THEN
                  F_RAIN(I,J,K)=QR(I,J,K)/(QC(I,J,K)+QR(I,J,K))
                ELSE
                  F_RAIN(I,J,K)=0.
                ENDIF
              ENDDO
             ENDDO
            ENDDO
          ENDIF spec_adv_fer
!
!----------------------------------------------------------------------
        CASE ('gfs')       !-- Update fields for GFS microphysics
!----------------------------------------------------------------------
!
          spec_adv_gfs: IF (.NOT.SPEC_ADV .OR. CLD_INIT) THEN
            cld_init_gfs: IF (CLD_INIT) THEN
!-- Initialize F_ICE, F_RAIN, & F_RIMEF arrays
              IF (SPEC_ADV) THEN
                WRITE(0,*) 'Never ran GFS microphysics with SPEC_ADV=T.'   &
                          ,'  Use at your own risk.'
              ENDIF
              DO K=1,LM
               DO J=JMS,JME
                DO I=IMS,IME
                  F_RAIN(I,J,K)=0.
                  F_RIMEF(I,J,K)=1.
                  IF (CWM(I,J,K)>EPSQ .AND. T(I,J,K)<233.15) THEN
                    F_ICE(I,J,K)=1.
                  ELSE
                    F_ICE(I,J,K)=0.
                  ENDIF
                ENDDO
               ENDDO
              ENDDO
            ENDIF  cld_init_gfs
!-- Update WATER arrays (QC,QI) when advecting only total condensate (spec_adv=F)
!   or initialize them at the start of the forecast (CLD_INIT=T).
            DO K=1,LM
             DO J=JMS,JME
              DO I=IMS,IME
                IF (CWM(I,J,K)>EPSQ) THEN
                  QC(I,J,K)=(1.-F_ice(I,J,K))*CWM(I,J,K)
                  QI(I,J,K)=F_ice(I,J,K)*CWM(I,J,K)
                ELSE
                  QC(I,J,K)=0.
                  QI(I,J,K)=0.
                ENDIF
              ENDDO
             ENDDO
            ENDDO
          ELSE spec_adv_gfs
!-- Update CWM, F_ICE arrays from separate species advection (spec_adv=T)
            DO K=1,LM
             DO J=JMS,JME
              DO I=IMS,IME
                CWM(I,J,K)=QC(I,J,K)+QI(I,J,K)
                IF (CWM(I,J,K)>EPSQ) THEN
                  F_ICE(I,J,K)=QI(I,J,K)/CWM(I,J,K)
                ELSE
                  F_ICE(I,J,K)=0.
                ENDIF
              ENDDO
             ENDDO
            ENDDO
          ENDIF  spec_adv_gfs
!
!----------------------------------------------------------------------
        CASE ('wsm6')      !-- Update fields for WSM6 microphysics
!----------------------------------------------------------------------
!
          init_adv_wsm6: IF (CLD_INIT) THEN
!-- Assume only cloud ice is present at initial time
            DO K=1,LM
             DO J=JMS,JME
              DO I=IMS,IME
                QS(I,J,K)=0.0
                QG(I,J,K)=0.0
                IF (CWM(I,J,K)>EPSQ) THEN
                  LIQW=(1.-F_ice(I,J,K))*CWM(I,J,K)
                  QC(I,J,K)=(1.-F_rain(I,J,K))*LIQW
                  QR(I,J,K)=F_rain(I,J,K)*LIQW
                  QI(I,J,K)=F_ice(I,J,K)*CWM(I,J,K)
                ELSE
                  QC(I,J,K)=0.
                  QR(I,J,K)=0.
                  QI(I,J,K)=0.
                ENDIF
              ENDDO
             ENDDO
            ENDDO
          ELSE init_adv_wsm6
            notspec_adv_wsm6: IF (.NOT.SPEC_ADV) THEN
!-- Update WATER arrays (QC,QR,...) when advecting only total condensate (spec_adv=F).
!-- Assume fraction of each water category is unchanged by advection. 
              DO K=1,LM
               DO J=JMS,JME
                DO I=IMS,IME
                  OLDCWM=QC(I,J,K)+QR(I,J,K)   &
                        +QI(I,J,K)+QS(I,J,K)   &
                        +QG(I,J,K)
                  IF (OLDCWM>EPSQ) THEN
                    FRACTION=CWM(I,J,K)/OLDCWM
                    QC(I,J,K)=FRACTION*QC(I,J,K)
                    QR(I,J,K)=FRACTION*QR(I,J,K)
                    QI(I,J,K)=FRACTION*QI(I,J,K)
                    QS(I,J,K)=FRACTION*QS(I,J,K)
                    QG(I,J,K)=FRACTION*QG(I,J,K)
                  ELSE
                    QC(I,J,K)=0.0
                    QR(I,J,K)=0.0
                    QI(I,J,K)=0.0
                    QS(I,J,K)=0.0
                    QG(I,J,K)=0.0
                    IF (T(I,J,K)<233.15) THEN
                      QI(I,J,K)=CWM(I,J,K)
                    ELSE
                      QC(I,J,K)=CWM(I,J,K)
                    ENDIF
                  ENDIF
                ENDDO
               ENDDO
              ENDDO
            ENDIF  notspec_adv_wsm6
!
!-- Couple QC,QR,... <=> CWM,F_ice,F_rain,F_RimeF arrays
!-- Update CWM,F_XXX arrays from separate species advection (spec_adv=T)
!
            DO K=1,LM
             DO J=JMS,JME
              DO I=IMS,IME
                CWM(I,J,K)=QC(I,J,K)+QR(I,J,K)      &
                          +QI(I,J,K)+QS(I,J,K)      &
                          +QG(I,J,K)
                IF (CWM(I,J,K)>EPSQ) THEN
                  LIQW=QI(I,J,K)+QS(I,J,K)+QG(I,J,K)
                  F_ICE(I,J,K)=LIQW/CWM(I,J,K)
                ELSE
                  F_ICE(I,J,K)=0.
                ENDIF
                IF (QR(I,J,K)>EPSQ) THEN
                  F_RAIN(I,J,K)=QR(I,J,K)/(QC(I,J,K)+QR(I,J,K))
                ELSE
                  F_RAIN(I,J,K)=0.
                ENDIF
                IF (QG(I,J,K)>EPSQ) THEN
!-- Update F_RIMEF: assume 5x higher graupel density (500 kg/m**3) vs snow (100 kg/m**3)
                  LIQW=5.*QG(I,J,K)+QS(I,J,K)
                  F_RIMEF(I,J,K)=LIQW/(QS(I,J,K)+QG(I,J,K))
                ELSE
                  F_RIMEF(I,J,K)=1.
                ENDIF
              ENDDO
             ENDDO
            ENDDO
!
          ENDIF init_adv_wsm6
!
!----------------------------------------------------------------------
        CASE ('thompson')   !-- Update fields for Thompson microphysics
!----------------------------------------------------------------------
!
!+---+-----------------------------------------------------------------+
!..The CLD_INIT test provides a way to translate initial values of CWM
!.. into coomponent species of cloud water, rain, and ice, but not snow
!.. or graupel. Thompson MP will pretty rapidly make snow from the
!.. cloud ice field.  Next IF-test is whether individual species
!.. advection is enabled, which almost certainly should be the case when
!.. picking this scheme.  In this case, the separate species are summed
!.. into the CWM and ice, rain, and rime variables are computed only for
!.. consistency with other schemes.  But, if single species advection is
!.. not enabled, then each t-step the CWM array needs to be split into
!.. component species to prepare MP routine to have some semblance of
!.. proper individual species.  Again, this is strongly discouraged.
!+---+-----------------------------------------------------------------+
          spec_adv_thompson: IF (CLD_INIT) THEN
             DO K=1,LM
                DO J=JMS,JME
                DO I=IMS,IME
                   QS(I,J,K)=0.0
                   QG(I,J,K)=0.0
                   IF (CWM(I,J,K) .gt. EPSQ) THEN
                      LIQW=(1.-F_ice(I,J,K))*CWM(I,J,K)
                      QC(I,J,K)=(1.-F_rain(I,J,K))*LIQW
                      QR(I,J,K)=F_rain(I,J,K)*LIQW
                      QI(I,J,K)=F_ice(I,J,K)*CWM(I,J,K)
                   ELSE
                      QC(I,J,K)=0.
                      QR(I,J,K)=0.
                      QI(I,J,K)=0.
                   ENDIF
                ENDDO
                ENDDO
             ENDDO
          ELSE IF(SPEC_ADV) THEN  spec_adv_thompson
             DO K=1,LM
                DO J=JMS,JME
                DO I=IMS,IME
                   CWM(I,J,K) = QC(I,J,K)+QR(I,J,K)     &
                              + QI(I,J,K)                       &
                              + QS(I,J,K)+QG(I,J,K)
                   IF (CWM(I,J,K) .gt. EPSQ) THEN
                      LIQW = MAX(0., CWM(I,J,K) - QI(I,J,K)     &
                                                - QS(I,J,K)     &
                                                - QG(I,J,K))
                      F_ICE(I,J,K) = MAX(0., 1.0 - LIQW/CWM(I,J,K))
                      IF (QR(I,J,K) .gt. EPSQ) THEN
                         F_RAIN(I,J,K) = QR(I,J,K)              &
                                 / (QC(I,J,K)+QR(I,J,K))
                      ELSE
                         F_RAIN(I,J,K)=0.
                      ENDIF
                      IF (QG(I,J,K) .gt. EPSQ) THEN
                         F_RIMEF(I,J,K) = (5.*QG(I,J,K)         &
                                        +     QS(I,J,K))        &
                                        / (QS(I,J,K)            &
                                        +  QG(I,J,K))
                      ELSE
                         F_RIMEF(I,J,K)=1.
                      ENDIF
                   ELSE
                      F_ICE(I,J,K) = 0.
                      F_RAIN(I,J,K)=0.
                      F_RIMEF(I,J,K)=1.
                      CWM(I,J,K) = 0.
                   ENDIF
                ENDDO
                ENDDO
             ENDDO
          ELSE  spec_adv_thompson
            ! write(0,*) 'WARNING: This option is STRONGLY DISCOURAGED'
            ! write(0,*) '  please consider using full advection of all'
            ! write(0,*) '  species when picking Thompson microphysics.'
             DO J=JMS,JME
             DO I=IMS,IME
                DO K=LM,1,-1
                   deep_ice = .false.
                   IF (CWM(I,J,K) .gt. EPSQ) THEN
                      OLDCWM  = QC(I,J,K)+QR(I,J,K)     &
                              + QI(I,J,K)                       &
                              + QS(I,J,K)+QG(I,J,K)
                      IF (OLDCWM .gt. EPSQ) THEN
                         LIQW = MAX(0., OLDCWM - QI(I,J,K)      &
                                               - QS(I,J,K)      &
                                               - QG(I,J,K))
                         F_ICE(I,J,K) = MAX(0., 1.0 - LIQW/OLDCWM)
                         IF (QR(I,J,K) .gt. EPSQ) THEN
                            F_RAIN(I,J,K) = QR(I,J,K)           &
                                 / (QC(I,J,K)+QR(I,J,K))
                         ELSE
                            F_RAIN(I,J,K)=0.
                         ENDIF
                         IF (QG(I,J,K) .gt. EPSQ) THEN
                            F_RIMEF(I,J,K) = (5.*QG(I,J,K)      &
                                           +     QS(I,J,K))     &
                                           / (QS(I,J,K)         &
                                           +  QG(I,J,K))
                         ELSE
                            F_RIMEF(I,J,K)=1.
                         ENDIF
                         LIQW = MAX(0., (1.-F_ICE(I,J,K))*CWM(I,J,K))
                         QR(I,J,K) = LIQW*F_RAIN(I,J,K)*CWM(I,J,K)
                         QC(I,J,K) = LIQW*(1.-F_RAIN(I,J,K))*CWM(I,J,K)
                         IF (QG(I,J,K) .gt. EPSQ) THEN
                            FRACTION = MAX(0., MIN(QG(I,J,K)            &
                                       / (QG(I,J,K)+QS(I,J,K)), 1.) )
                         ELSE
                            FRACTION = 0.
                         ENDIF
                         QG(I,J,K) = FRACTION*F_ICE(I,J,K)*CWM(I,J,K)
                         QI(I,J,K) = 0.1*(1.-FRACTION)*F_ICE(I,J,K)*CWM(I,J,K)
                         QS(I,J,K) = 0.9*(1.-FRACTION)*F_ICE(I,J,K)*CWM(I,J,K)

                      ELSE       ! Below, the condensate is all new here
                         QC(I,J,K) = 0.0
                         QI(I,J,K) = 0.0
                         QR(I,J,K) = 0.0
                         QS(I,J,K) = 0.0
                         QG(I,J,K) = 0.0
                         IF (T(I,J,K) .le. 235.15) THEN
                            QI(I,J,K) = 0.5*CWM(I,J,K)
                            QS(I,J,K) = 0.5*CWM(I,J,K)
                         ELSEIF (T(I,J,K) .le. 258.15) THEN
                            QI(I,J,K) = 0.1*CWM(I,J,K)
                            QS(I,J,K) = 0.9*CWM(I,J,K)
                            deep_ice = .true.
                         ELSEIF (T(I,J,K) .le. 275.15) THEN
                            if (deep_ice .and. T(I,J,K).lt.273.15) then
                               QS(I,J,K) = CWM(I,J,K)
                            elseif (deep_ice .and. T(I,J,K).lt.274.15) then
                               QS(I,J,K) = 0.333*CWM(I,J,K)
                               QR(I,J,K) = 0.667*CWM(I,J,K)
                            elseif (deep_ice) then
                               QS(I,J,K) = 0.1*CWM(I,J,K)
                               QR(I,J,K) = 0.9*CWM(I,J,K)
                            else
                               QC(I,J,K) = CWM(I,J,K)
                            endif
                         ELSE
                            QC(I,J,K) = CWM(I,J,K)
                         ENDIF
                         LIQW = MAX(0., CWM(I,J,K) - QI(I,J,K)  &
                                                   - QS(I,J,K)  &
                                                   - QG(I,J,K))
                         IF (CWM(I,J,K) .gt. EPSQ) THEN
                            F_ICE(I,J,K) = (1.0-LIQW)/CWM(I,J,K)
                         ELSE
                            F_ICE(I,J,K) = 0.
                         ENDIF
                         IF (QR(I,J,K) .gt. EPSQ) THEN
                            F_RAIN(I,J,K) = QR(I,J,K)           &
                                    / (QC(I,J,K)+QR(I,J,K))
                         ELSE
                            F_RAIN(I,J,K)=0.
                         ENDIF
                         IF (QG(I,J,K) .gt. EPSQ) THEN
                            F_RIMEF(I,J,K) = (5.*QG(I,J,K)      &
                                           +     QS(I,J,K))     &
                                           / (QS(I,J,K)         &
                                           +  QG(I,J,K))
                         ELSE
                            F_RIMEF(I,J,K)=1.
                         ENDIF
                      ENDIF
                   ELSE
                      QC(I,J,K) = 0.0
                      QR(I,J,K) = 0.0
                      QI(I,J,K) = 0.0
                      QS(I,J,K) = 0.0
                      QG(I,J,K) = 0.0
                      F_ICE(I,J,K) = 0.0
                      F_RAIN(I,J,K) = 0.0
                      F_RIMEF(I,J,K) = 1.0
                   ENDIF
                ENDDO
             ENDDO
             ENDDO
          ENDIF  spec_adv_thompson

!
!----------------------------------------------------------------------
        CASE DEFAULT
!----------------------------------------------------------------------
!
          IF (CLD_INIT) THEN
            WRITE(0,*) 'Do nothing for default option'
          ENDIF
!
      END SELECT   ! MICROPHYSICS
!
!----------------------------------------------------------------------
!
      END SUBROUTINE UPDATE_WATER

!----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------

!
      SUBROUTINE CALC_RH_RADAR_DFI(T,Q,PD,PSGML1,SGML2      &
                                  ,R_D,R_V,RH_HOLD          &
                                  ,IMS,IME,JMS,JME,LM       &
                                  ,IFLAG)

       USE MODULE_MP_ETANEW, ONLY : FERRIER_INIT,GPVS,FPVS  &
                                   ,FPVS0,NX,RQR_DRmin      &
                                   ,RQR_DRmax,MASSI,CN0R0   &
                                   ,CN0r_DMRmin,CN0r_DMRmax

       IMPLICIT NONE

!tst
       INTEGER, INTENT(IN):: IMS,IME,JMS,JME, LM

       REAL :: T(IMS:IME,JMS:JME,1:LM)
       REAL :: Q(IMS:IME,JMS:JME,1:LM)
       REAL :: RH_HOLD(IMS:IME,JMS:JME,1:LM)
       REAL :: PD(IMS:IME,JMS:JME)
       REAL :: PSGML1(LM),SGML2(LM)
       REAL :: R_D,R_V,PMID,VPRES,SATVPRES, EPS, DEN
       INTEGER :: IFLAG,I,J,L

        EPS=R_D/R_V

        IF (IFLAG == 1) THEN
        DO L=1,LM
          DO J=JMS,JME
            DO I=IMS,IME
            PMID=SGML2(L)*PD(I,J)+PSGML1(L)
            DEN=EPS+Q(I,J,L)*(1.-EPS)
            VPRES=PMID*Q(I,J,L)/DEN
            SATVPRES=1.E3*FPVS0(T(I,J,L))
            RH_HOLD(I,J,L)=VPRES/SATVPRES
            ENDDO
          ENDDO
        ENDDO
        ENDIF

        IF (IFLAG == -1) THEN
        DO L=1,LM
          DO J=JMS,JME
            DO I=IMS,IME
            SATVPRES=1.E3*FPVS0(T(I,J,L))
            VPRES=SATVPRES*RH_HOLD(I,J,L)
            PMID=SGML2(L)*PD(I,J)+PSGML1(L)
            DEN=PMID-VPRES*(1.-EPS)
            Q(I,J,L)=VPRES*EPS/DEN
            ENDDO
          ENDDO
        ENDDO
        ENDIF

      END SUBROUTINE CALC_RH_RADAR_DFI

 
!----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
 
      SUBROUTINE CLTEND (ICLTEND,NPRECIP, T,Told,Tadj                    &
                        ,IDS,IDE,JDS,JDE,LM                              &
                        ,IMS,IME,JMS,JME                                 &
                        ,ITS,ITE,JTS,JTE)
!----------------------------------------------------------------------
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .     
! SUBPROGRAM:    CLTEND      TEMPERATURE CHANGE BY CLOUD PROCESSES
!   PRGRMMR: FERRIER         ORG: W/NP22     DATE: 01-09-26
!     
! ABSTRACT:
!     CLTEND GRADUALLY UPDATES TEMPERATURE TENDENCIES FROM CONVECTION 
!     AND GRID-SCALE MICROPHYSICS.
!     
! USAGE: CALL CLTEND FROM SOLVER_RUN
!   INPUT ARGUMENT LIST:
!     ICLTEND - FLAG SET TO -1 PRIOR TO PHYSICS CALLS, 0 AFTER PHYSICS
!               CALLS, AND 1 FOR UPDATING TEMPERATURES EVERY TIME STEP
!  
!   OUTPUT ARGUMENT LIST:  NONE
!     
!   OUTPUT FILES:  NONE
!     
!   SUBPROGRAMS CALLED:  NONE
!  
!   UNIQUE: NONE
!  
!   LIBRARY: NONE
!  
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!$$$  
!----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!----------------------------------------------------------------------
!
      INTEGER,INTENT(IN) :: ICLTEND,NPRECIP                            &
                           ,IDS,IDE,JDS,JDE,LM                         &
                           ,IMS,IME,JMS,JME                            &
                           ,ITS,ITE,JTS,JTE
!
      REAL,DIMENSION(IMS:IME,JMS:JME,1:LM),INTENT(INOUT) :: T          &
                                                           ,Tadj       &
                                                           ,Told
!
!***  LOCAL VARIABLES 
!
      INTEGER :: I,J,K
      REAL :: RDTPH
!
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!
      IF(ICLTEND<0)THEN
         DO K=1,LM
         DO J=JTS,JTE
         DO I=ITS,ITE
            Told(I,J,K)=T(I,J,K)
         ENDDO
         ENDDO
         ENDDO
      ELSE IF(ICLTEND==0)THEN
         RDTPH=1./REAL(NPRECIP)
         DO K=1,LM
         DO J=JTS,JTE
         DO I=ITS,ITE
            Tadj(I,J,K)=RDTPH*(T(I,J,K)-Told(I,J,K))
            T(I,J,K)=Told(I,J,K)
         ENDDO
         ENDDO
         ENDDO
      ELSE
         DO K=1,LM
         DO J=JTS,JTE
         DO I=ITS,ITE
            T(I,J,K)=T(I,J,K)+Tadj(I,J,K)
         ENDDO
         ENDDO
         ENDDO
      ENDIF
!----------------------------------------------------------------------
!
      END SUBROUTINE CLTEND
!
!-----------------------------------------------------------------------
!----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------

      SUBROUTINE RIME_FACTOR_UPDATE (RIME_FACTOR_INPUT                  &
                                    ,QS,QG,F_RIMEF                      &
                                    ,IDS,IDE,JDS,JDE,LM                 &
                                    ,IMS,IME,JMS,JME                    &
                                    ,ITS,ITE,JTS,JTE)
!----------------------------------------------------------------------
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .
! SUBPROGRAM:  RIME_FACTOR_UPDATE
!   PRGRMMR: FERRIER         ORG: W/NP22     DATE: 2013-06-14
!
! ABSTRACT:
!
!     UPDATES THE RIME FACTOR ARRAY AFTER 3D ADVECTION
!
! USAGE: CALL CLTEND FROM SOLVER_RUN
!   INPUT ARGUMENT LIST:
!     RIME_FACTOR_INPUT= TRUE BEFORE ADVECTION, RIME_FACTOR IS INPUT
!     RIME_FACTOR_INPUT=FALSE BEFORE ADVECTION, RIME FACTOR IS OUTPUT
!
!   OUTPUT ARGUMENT LIST:  NONE
!
!   OUTPUT FILES:  NONE
!
!   SUBPROGRAMS CALLED:  NONE
!
!   UNIQUE: NONE
!
!   LIBRARY: NONE
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE : IBM SP
!$$$
!----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!----------------------------------------------------------------------
!
      LOGICAL,INTENT(IN) :: RIME_FACTOR_INPUT
!
      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,LM                         &
                           ,IMS,IME,JMS,JME                            &
                           ,ITS,ITE,JTS,JTE
!
      REAL,DIMENSION(IMS:IME,JMS:JME,1:LM),INTENT(INOUT) :: F_RIMEF
!
      REAL,DIMENSION(IMS:IME,JMS:JME,1:LM),INTENT(INOUT) :: QS,QG
!
!***  LOCAL VARIABLES
!
      INTEGER :: I,J,K
      REAL :: RIMEF
!
!----------------------------------------------------------------------
      IF (RIME_FACTOR_INPUT) THEN     !-- Before advection
         DO K=1,LM
           DO J=JTS,JTE
             DO I=ITS,ITE
                QG(I,J,K)=QS(I,J,K)*F_RIMEF(I,J,K)
             ENDDO
           ENDDO
         ENDDO
!
         CALL HALO_EXCH(QG,LM,2,2)
!

      ELSE                            !-- After advection
         DO K=1,LM
           DO J=JMS,JME
             DO I=IMS,IME
                IF (QG(I,J,K)>EPSQ .AND.                        &
                    QS(I,J,K)>EPSQ) THEN
                   RIMEF=QG(I,J,K)/QS(I,J,K)
                   F_RIMEF(I,J,K)=MIN(50., MAX(1.,RIMEF) )
                ELSE
                   F_RIMEF(I,J,K)=1.
                ENDIF
             ENDDO
           ENDDO
         ENDDO
      ENDIF
      END SUBROUTINE RIME_FACTOR_UPDATE
!
!----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
 
      SUBROUTINE TQADJUST(T,Q,QC,CWM,F_ICE,F_RAIN                       &
                         ,PD,DSG2,PDSG1,PSGML1,SGML2                    &
                         ,SPEC_ADV,RHgrd                                &
                         ,IDS,IDE,JDS,JDE,LM                            &
                         ,IMS,IME,JMS,JME                               &
                         ,ITS,ITE,JTS,JTE)
!***********************************************************************
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .
! SUBPROGRAM:    TQADJUST             TQADJUST
!   PRGRMMR: FERRIER         ORG: NP22     DATE: 5 APR 2016
!
! ABSTRACT:
!     Smooth temperature profiles when lapse rates exceed dry adiabatic
!     above PBL, prevent supersaturation with respect to water.
!
! PROGRAM HISTORY LOG (with changes to called routines) :
!   2016-04     FERRIER, JANJIC  - Smooth T profiles, prevent supersaturation
!
! USAGE: CALL TQADJUST FROM PHY_RUN
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!-----------------------------------------------------------------------
!
      USE MODULE_CONSTANTS,ONLY : CAPPA,CP,EP_2,EPSQ,R_d,R_v,CPV,CLIQ,  &
                                  A2,A4,PSAT,XLV,TIW
      USE MODULE_MP_ETANEW, ONLY : FPVS0
!
!-----------------------------------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------------------------------
!
!----------------------
!-- Input argument variables
!----------------------
!
      REAL,DIMENSION(IMS:IME,JMS:JME,1:LM),INTENT(INOUT) ::             &
                                                T,Q,QC,CWM,F_ICE,F_RAIN
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: PD
      REAL,DIMENSION(1:LM),INTENT(IN) :: DSG2,PDSG1,PSGML1,SGML2
      REAL,INTENT(IN) :: RHgrd
      LOGICAL,INTENT(IN) :: SPEC_ADV
      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,LM                          &
                           ,IMS,IME,JMS,JME                             &
                           ,ITS,ITE,JTS,JTE
!
!--  Local Variables
!
      INTEGER :: I,J,K,LM2,KMIX,KTHmin,KBOT,KTOP,ITmax,ITER,ITRmax
      REAL :: RCP,RRV,TK,PP,QV,QCW,TREF,ESW,QSW,DQsat,SSat,DTHmin,COND, &
              Qrain,Qice,Qliq
      REAL,DIMENSION(1:LM) :: Tcol,Pcol,QVcol,QCcol,EXNcol,THcol,DPcol, &
                              DTHcol,Fcol
      LOGICAL :: LRFilt,SSFilt
!
      REAL,PARAMETER :: SupSat=0.001, SubSat=-SupSat, DTHthresh=-0.01,  &
        TTP=TIW+0.01, XA=(CLIQ-CPV)/R_V, XB=XA+XLV/(R_V*TTP),           &
        XLV1=XLV/CP, XLV2=XLV1*XLV/R_V
!
!-----------------------------------------------------------------------
!
      ITmax=LM/5
      LM2=LM-2
!
!-----------------------------------------------------------------------
!--  Main loop through I, J,  ------------------------------------------
!-----------------------------------------------------------------------
!
      DO J=JTS,JTE
        DO I=ITS,ITE
!
          LRFilt=.FALSE.       ! Lapse rate flag (full column)
          SSFilt=.FALSE.       ! Supersaturation flag
          IF(SPEC_ADV) THEN
            DO K=1,LM
              QCcol(K)=QC(I,J,K)
            ENDDO
          ELSE
            DO K=1,LM
              QCcol(K)=CWM(I,J,K)*(1.-F_ICE(I,J,K))*(1.-F_RAIN(I,J,K))
            ENDDO
          ENDIF
          DO K=1,LM
            Tcol(K)=T(I,J,K)
            Pcol(K)=SGML2(K)*PD(I,J)+PSGML1(K)
            QVcol(K)=Q(I,J,K)/(1.-Q(I,J,K))        ! Water vapor mixing ratio
            EXNcol(K)=(1.E5/Pcol(K))**CAPPA
          ENDDO
!
!-----------------------------------------------------------------------
!-- Ferrier-Aligo condensation/evaporation algorithm - 1st of 2 times
!-----------------------------------------------------------------------
!
SSadj1:   DO K=1,LM
            TK=Tcol(K)                                     ! Temperature (deg K)
            PP=Pcol(K)                                     ! Pressure (Pa)
            QV=QVcol(K)                                    ! Water vapor mixing ratio
            QCW=QCcol(K)                                   ! Cloud water mixing ratio
!            TREF=TTP/TK                                    ! WSM6
!            ESW=PSAT*EXP(LOG(TREF)*(XA))*EXP(XB*(1.-TREF)) ! WSM6
!            ESW=1000.*FPVS0(TK)                            ! Old global tables
            ESW=PSAT*EXP(A2*(TK-TTP)/(TK-A4))              ! Magnus Tetens
!            TREF=TK-TIW                                    ! Bolton (1980)
!           ESW=611.2*EXP(17.67*TREF/(TREF+243.5))          ! Bolton (1980)
            ESW=MIN(ESW,0.99*PP)                           ! Saturation vapor pressure (water)
            QSW=RHgrd*EP_2*ESW/(PP-ESW)                    ! Saturation mixing ratio (water)
            DQsat=QV-QSW                                   ! Excess QV above saturation
            SSat=DQsat/QSW                                 ! Grid-scale supersaturation ratio
SSrem1:     IF(SSat>SupSat .OR.                     &      ! Remove supersaturation if SSat>0.1%
               (QCW>EPSQ .AND. SSat<SubSat) ) THEN         ! Adjust to saturation if SSat<0.1% w/ cloud water
              SSFilt=.TRUE.                                ! Supersaturation flag
cond_iter1:   DO ITER=1,10                                 ! Usually converges in <=3 iterations
                COND=DQsat/(1.+XLV2*QSW/(TK*TK))           ! Asai (1965, J. Japan)
                COND=MAX(COND, -QCW)                       ! Limit cloud water evaporation
                TK=TK+XLV1*COND                            ! Update temperature
                QV=QV-COND                                 ! Update water vapor mixing ratio
                QCW=QCW+COND                               ! Update cloud water mixing ratio
!                TREF=TTP/TK                                ! WSM6
!                ESW=PSAT*EXP(LOG(TREF)*(XA))*EXP(XB*(1.-TREF)) ! WSM6
!                ESW=1000.*FPVS0(TK)                        ! Old global tables
                ESW=PSAT*EXP(A2*(TK-TTP)/(TK-A4))          ! Magnus Tetens
!                TREF=TK-TIW                                ! Bolton (1980)
!                ESW=611.2*EXP(17.67*TREF/(TREF+243.5))     ! Bolton (1980)
                ESW=MIN(ESW,0.99*PP)                       ! Saturation vapor pressure (water)
                QSW=RHgrd*EP_2*ESW/(PP-ESW)                ! Water saturation mixing ratio
                DQsat=QV-QSW                               ! Excess QV above saturation
                SSat=DQsat/QSW                             ! Grid-scale supersaturation ratio
                IF (SSat>=SubSat .AND. SSat<=SupSat) EXIT  ! Exit if -0.1%<SSat<0.1%
                IF (SSat<SubSat .AND. QCW<=EPSQ)     EXIT  ! Exit if SSat<-0.1% & no cloud water
              ENDDO  cond_iter1                            ! 1st *cond*ensation *iter*ation
              IF (ITER<=10) THEN
                Tcol(K)=TK
                QVcol(K)=QV
                QCcol(K)=QCW
              ENDIF
            ENDIF  SSrem1
          ENDDO  SSadj1
!
          DO K=1,LM
            THcol(K)=Tcol(K)*EXNcol(K)
          ENDDO
!
          DTHcol(1)=1.
          DO K=2,LM
            DTHcol(K)=THcol(K-1)-THcol(K)
          ENDDO
!
          KMIX=0
          DO K=LM2,2,-1
            IF(DTHcol(K)>0.) THEN
!-- Start above the well-mixed layer immediately above the
!   surface where theta may decrease with height 
              KMIX=K
              EXIT
            ENDIF
          ENDDO
!
!*************************
LRadjust: IF (KMIX>2) THEN
!*************************
!
            KTOP=0
            DO K=3,KMIX
              IF(DTHcol(K)<DTHthresh) THEN
                KTOP=K-1            !- Level at top of highest unstable layer
                EXIT
              ENDIF
            ENDDO
!
!-------------------------
Maybe_mix:  IF(KTOP>0) THEN
!-------------------------
!
              KBOT=0
              DO K=KMIX,2,-1
                IF(DTHcol(K)<DTHthresh) THEN
                  KBOT=K            !- Lowest unstable layer
                  EXIT
                ENDIF
              ENDDO
              IF(KBOT>0) THEN
                LRFilt=.TRUE.       !- For the full column (any layer)
                ITRmax=ITmax
              ELSE
                ITRmax=0            !- Do not mix
              ENDIF
!
              DO K=1,LM
                DPcol(K)=DSG2(K)*PD(I,J)+PDSG1(K)  ! Hydrostatic pressure thickness
                Fcol(K)=THcol(K)    !- Fcol, modified theta
              ENDDO
!
!- - - - - - - - - - - - -
Mix_lyrs:     DO ITER=1,ITRmax
!- - - - - - - - - - - - -
                DO K=KTOP,KBOT
                  IF(DTHcol(K)<DTHthresh) THEN
                    IF(DTHcol(K+1)<DTHthresh) THEN
!-- Mix 3 layers if current layer and layers above and below are unstable
                      Fcol(K)=(THcol(K-1)*DPcol(K-1)+THcol(K)*DPcol(K)    &
                               +THcol(K+1)*DPcol(K+1))/                   &
                              (DPcol(K-1)+DPcol(K)+DPcol(K+1))
                    ELSE
!-- Mix with higher layer if current layer is unstable
                      Fcol(K)=(THcol(K-1)*DPcol(K-1)+THcol(K)*DPcol(K))/  &
                              (DPcol(K-1)+DPcol(K))
                    ENDIF
                  ELSE IF(DTHcol(K+1)<DTHthresh) THEN
!-- Mix with lower layer if it is unstable
                    Fcol(K)=(THcol(K)*DPcol(K)+THcol(K+1)*DPcol(K+1))/    &
                            (DPcol(K)+DPcol(K+1))
                  ENDIF
!-- Do nothing if the current layer or the layer below is not unstable
                ENDDO
!
                DO K=KTOP,KBOT
                  THcol(K)=Fcol(K)
                ENDDO
                DO K=KTOP,KBOT
                  DTHcol(K)=THcol(K-1)-THcol(K)
                ENDDO
!
                KTOP=0
                DO K=3,KMIX
                  IF(DTHcol(K)<DTHthresh) THEN
                    KTOP=K-1         !- Level at top of highest unstable layer
                    EXIT
                  ENDIF
                ENDDO
!
                IF(KTOP<=0) EXIT Mix_lyrs   !- Exit with no unstable layer in column
!
                KBOT=0
                DO K=KMIX,2,-1
                  IF(DTHcol(K)<DTHthresh) THEN
                    KBOT=K           !- Lowest unstable layer
                    EXIT
                  ENDIF
                ENDDO
!- - - - - - - - - - - - -
              ENDDO  Mix_lyrs   !DO ITER
!- - - - - - - - - - - - -
              DO K=1,LM
                Tcol(K)=THcol(K)/EXNcol(K)     !-- Update T
              ENDDO
!-------------------------
            ENDIF  Maybe_mix    !IF(KTOP>0)
!-------------------------
!*************************
          ENDIF  LRadjust
!*************************
!
!-----------------------------------------------------------------------
!-- Ferrier-Aligo condensation/evaporation algorithm - 2nd of 2 times
!-----------------------------------------------------------------------
!
SSadj2:   DO K=1,LM
            TK=Tcol(K)                                     ! Temperature (deg K)
            PP=Pcol(K)                                     ! Pressure (Pa)
            QV=QVcol(K)                                    ! Water vapor mixing ratio
            QCW=QCcol(K)                                   ! Cloud water mixing ratio
!            TREF=TTP/TK                                    ! WSM6
!            ESW=PSAT*EXP(LOG(TREF)*(XA))*EXP(XB*(1.-TREF)) ! WSM6
!            ESW=1000.*FPVS0(TK)                            ! Old global tables
            ESW=PSAT*EXP(A2*(TK-TTP)/(TK-A4))              ! Magnus Tetens
!            TREF=TK-TIW                                    ! Bolton (1980)
!           ESW=611.2*EXP(17.67*TREF/(TREF+243.5))          ! Bolton (1980)
            ESW=MIN(ESW,0.99*PP)                           ! Saturation vapor pressure (water)
            QSW=RHgrd*EP_2*ESW/(PP-ESW)                    ! Saturation mixing ratio (water)
            DQsat=QV-QSW                                   ! Excess QV above saturation
            SSat=DQsat/QSW                                 ! Grid-scale supersaturation ratio
SSrem2:     IF(SSat>SupSat .OR.                     &      ! Remove supersaturation if SSat>0.1%
               (QCW>EPSQ .AND. SSat<SubSat) ) THEN         ! Adjust to saturation if SSat<0.1% w/ cloud water
              SSFilt=.TRUE.                                ! Supersaturation flag
cond_iter2:   DO ITER=1,10                                 ! Usually converges in <=3 iterations
                COND=DQsat/(1.+XLV2*QSW/(TK*TK))           ! Asai (1965, J. Japan)
                COND=MAX(COND, -QCW)                       ! Limit cloud water evaporation
                TK=TK+XLV1*COND                            ! Update temperature
                QV=QV-COND                                 ! Update water vapor mixing ratio
                QCW=QCW+COND                               ! Update cloud water mixing ratio
!                TREF=TTP/TK                                ! WSM6
!                ESW=PSAT*EXP(LOG(TREF)*(XA))*EXP(XB*(1.-TREF)) ! WSM6
!                ESW=1000.*FPVS0(TK)                        ! Old global tables
                ESW=PSAT*EXP(A2*(TK-TTP)/(TK-A4))          ! Magnus Tetens
!                TREF=TK-TIW                                ! Bolton (1980)
!                ESW=611.2*EXP(17.67*TREF/(TREF+243.5))     ! Bolton (1980)
                ESW=MIN(ESW,0.99*PP)                       ! Saturation vapor pressure (water)
                QSW=RHgrd*EP_2*ESW/(PP-ESW)                ! Water saturation mixing ratio
                DQsat=QV-QSW                               ! Excess QV above saturation
                SSat=DQsat/QSW                             ! Grid-scale supersaturation ratio
                IF (SSat>=SubSat .AND. SSat<=SupSat) EXIT  ! Exit if -0.1%<SSat<0.1%
                IF (SSat<SubSat .AND. QCW<=EPSQ)     EXIT  ! Exit if SSat<-0.1% & no cloud water
              ENDDO  cond_iter2                            ! 2nd *cond*ensation *iter*ation
              IF (ITER<=10) THEN
                Tcol(K)=TK
                QVcol(K)=QV
                QCcol(K)=QCW
              ENDIF
            ENDIF  SSrem2
          ENDDO  SSadj2
!
!#######################################################################
!-- Update 3D arrays
!#######################################################################
!
adjust1:  IF (LRFilt .OR. SSFilt) THEN
            DO K=1,LM
              T(I,J,K)=Tcol(K)    !- Update T
            ENDDO
          ENDIF  adjust1
!
adjust2:  IF (SSFilt) THEN
            DO K=1,LM             !- Update Q
              Q(I,J,K)=QVcol(K)/(1.+QVcol(K))
            ENDDO
            IF(SPEC_ADV) THEN
              DO K=1,LM           !- Update QC
                QC(I,J,K)=QCcol(K)
              ENDDO
            ELSE
              DO K=1,LM           !- Update CWM, F_ICE, F_RAIN
                Qrain=CWM(I,J,K)*(1.-F_ICE(I,J,K))*F_RAIN(I,J,K)
                Qliq=QCcol(K)+Qrain
                Qice=CWM(I,J,K)*F_ICE(I,J,K)
                CWM(I,J,K)=Qliq+Qice
                IF(CWM(I,J,K)>EPSQ) F_ICE(I,J,K)=Qice/CWM(I,J,K)
                IF(Qliq>EPSQ) F_RAIN(I,J,K)=Qrain/Qliq
              ENDDO
            ENDIF
          ENDIF  adjust2
!
        ENDDO   !- I
      ENDDO     !- J
!-----------------------------------------------------------------------
!
      END SUBROUTINE TQADJUST
 
!----------------------------------------------------------------------
!######################################################################
!-----------------------------------------------------------------------
 
      END MODULE MODULE_SOLVER_GRID_COMP
!
!-----------------------------------------------------------------------
