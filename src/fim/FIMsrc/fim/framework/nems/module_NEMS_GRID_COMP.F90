!TBH:  Matches r14147 of https://svnemc.ncep.noaa.gov/projects/nems/trunk/ 
!TBH:  except for pe_member issues noted below.  
!JR commenting out the setting of pe_member so ESMF can stop on failure,

#include "./ESMFVersionDefine.h"

#if (ESMF_MAJOR_VERSION < 5 || ESMF_MINOR_VERSION < 2)
#undef ESMF_520rbs
#else
#define ESMF_520rbs
#endif

!-----------------------------------------------------------------------
!
      MODULE module_NEMS_GRID_COMP
!
!-----------------------------------------------------------------------
!***  This module contains codes directly related to the NEMS component.
!-----------------------------------------------------------------------
!
!***  The NEMS component lies in the heirarchy seen here:
!
!          Main program
!               |
!               |
!          NEMS component
!               |     |________________________.
!               |                              |
!          EARTH component        Ensemble Coupler component
!               |
!               |
!          ATM/OCEAN/ICE components
!               |
!               |
!          CORE component (GFS, NMM, FIM, GEN, etc.)
!
!-----------------------------------------------------------------------
!  2011-05-11  Theurich & Yang  - Modified for using the ESMF 5.2.0r_beta_snapshot_07.
!-----------------------------------------------------------------------
!
      USE ESMF_MOD
!
      USE module_NEMS_INTERNAL_STATE,ONLY: NEMS_INTERNAL_STATE          &
                                          ,WRAP_NEMS_INTERNAL_STATE
!
      USE ENS_CplComp_ESMFMod,ONLY: ENS_CplCompSetServices
!
      USE module_EARTH_GRID_COMP
!
      USE module_ERR_MSG,ONLY: ERR_MSG,MESSAGE_CHECK
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      PRIVATE
!
      PUBLIC :: NEMS_REGISTER
!
!-----------------------------------------------------------------------
!
      INTEGER :: MEMBER_ID                                              &
                ,TOTAL_MEMBER
!
      INTEGER :: HH_INCREASE                                            &
                ,HH_START                                               &
                ,HH_FINAL
!
      INTEGER :: NUMBER_START                                           &
                ,NUMBER_FINAL
!
      INTEGER,DIMENSION(:),   ALLOCATABLE :: PE_MEMBER                     !<-- Tasks for each member
      INTEGER,DIMENSION(:, :),ALLOCATABLE :: PETLIST                       !<-- Task list for each member
!
      LOGICAL :: ENS_SPS                                                   !<-- Control of Stochastic Perturbation Scheme (SPS)
!
      CHARACTER(ESMF_MAXSTR),DIMENSION(:),ALLOCATABLE :: IMP_EARTH_NAME    !<-- Import state name of the EARTH components
      CHARACTER(ESMF_MAXSTR),DIMENSION(:),ALLOCATABLE :: EXP_EARTH_NAME    !<-- Export state name of the EARTH components
      CHARACTER(ESMF_MAXSTR),DIMENSION(:),ALLOCATABLE :: GC_EARTH_NAME     !<-- Name of the EARTH component
!
      TYPE(NEMS_INTERNAL_STATE),POINTER,SAVE :: NEMS_INT_STATE
      TYPE(WRAP_NEMS_INTERNAL_STATE)   ,SAVE :: WRAP
!
      TYPE(ESMF_Clock), SAVE :: CLOCK_NEMS                                 !<-- The ESMF Clock of the NEMS component
      TYPE(ESMF_Config),SAVE :: CF_NEMS                                    !<-- The configure object of the NEMS component
      TYPE(ESMF_VM),    SAVE :: VM_GLOBAL
      TYPE(ESMF_TIME),  SAVE :: STARTTIME
!
      TYPE(ESMF_CplComp),SAVE :: ENS_CPL_COMP                              !<-- Ensemble Coupler components
      TYPE(ESMF_State),  SAVE :: ENS_CPL_IMP_STATE                         !<-- Import state of the Ensemble Coupler component
      TYPE(ESMF_State),  SAVE :: ENS_CPL_EXP_STATE                         !<-- Export state of the Ensemble Coupler component
!
      TYPE(ESMF_GridComp),DIMENSION(:),ALLOCATABLE,SAVE :: EARTH_GRID_COMP !<-- EARTH components for each member
      TYPE(ESMF_State),   DIMENSION(:),ALLOCATABLE,SAVE :: EARTH_IMP_STATE !<-- Import state of the EARTH component
      TYPE(ESMF_State),   DIMENSION(:),ALLOCATABLE,SAVE :: EARTH_EXP_STATE !<-- Export state of the EARTH component
!
!-----------------------------------------------------------------------
!
      CONTAINS

!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE NEMS_REGISTER(NEMS_GRID_COMP,RC_REG)
!
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(ESMF_GridComp)               :: NEMS_GRID_COMP                  !<-- The NEMS gridded component
      INTEGER            ,INTENT(OUT)   :: RC_REG                          !<-- Error return code
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER :: RC
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Register the NEMS Initialize, Run, and Finalize routines.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set Entry Point for NEMS Initialize"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifdef ESMF_3
      CALL ESMF_GridCompSetEntryPoint(NEMS_GRID_COMP                    &  !<-- The NEMS component
                                     ,ESMF_SETINIT                      &  !<-- Subroutine type (Initialize)
                                     ,NEMS_INITIALIZE                   &  !<-- User's subroutine name
                                     ,ESMF_SINGLEPHASE                  &
                                     ,RC)
#else
      CALL ESMF_GridCompSetEntryPoint(NEMS_GRID_COMP                    &  !<-- The NEMS component
                                     ,ESMF_SETINIT                      &  !<-- Subroutine type (Initialize)
                                     ,NEMS_INITIALIZE                   &  !<-- User's subroutine name
                                     ,phase=ESMF_SINGLEPHASE            &
                                     ,rc=RC)
#endif
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set Entry Point for NEMS Run"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifdef ESMF_3
      CALL ESMF_GridCompSetEntryPoint(NEMS_GRID_COMP                    &  !<-- The NEMS component
                                     ,ESMF_SETRUN                       &  !<-- Subroutine type (Run)
                                     ,NEMS_RUN                          &  !<-- User's subroutine name
                                     ,ESMF_SINGLEPHASE                  &
                                     ,RC)
#else
      CALL ESMF_GridCompSetEntryPoint(NEMS_GRID_COMP                    &  !<-- The NEMS component
                                     ,ESMF_SETRUN                       &  !<-- Subroutine type (Run)
                                     ,NEMS_RUN                          &  !<-- User's subroutine name
                                     ,phase=ESMF_SINGLEPHASE            &
                                     ,rc=RC)
#endif
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set Entry Point for NEMS Finalize"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifdef ESMF_3
      CALL ESMF_GridCompSetEntryPoint(NEMS_GRID_COMP                    &  !<-- The NEMS component
                                     ,ESMF_SETFINAL                     &  !<-- Subroutine type (Finalize)
                                     ,NEMS_FINALIZE                     &  !<-- User's subroutine name
                                     ,ESMF_SINGLEPHASE                  &
                                     ,RC)
#else
      CALL ESMF_GridCompSetEntryPoint(NEMS_GRID_COMP                    &  !<-- The NEMS component
                                     ,ESMF_SETFINAL                     &  !<-- Subroutine type (Finalize)
                                     ,NEMS_FINALIZE                     &  !<-- User's subroutine name
                                     ,phase=ESMF_SINGLEPHASE            &
                                     ,rc=RC)
#endif
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
      IF(RC_REG==ESMF_SUCCESS)THEN
!       WRITE(0,*)' NEMS_REGISTER succeeded'
      ELSE
        WRITE(0,*)' NEMS_REGISTER failed  RC_REG=',RC_REG
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE NEMS_REGISTER
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE NEMS_INITIALIZE(NEMS_GRID_COMP                         &
                                ,IMP_STATE                              &
                                ,EXP_STATE                              &
                                ,CLOCK_MAIN                             &
                                ,RC_INIT)
!
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(ESMF_GridComp)               :: NEMS_GRID_COMP                  !<-- The NEMS component
      TYPE(ESMF_State)                  :: IMP_STATE                       !<-- The NEMS import state
      TYPE(ESMF_State)                  :: EXP_STATE                       !<-- The NEMS export state
      TYPE(ESMF_Clock)                  :: CLOCK_MAIN                      !<-- The main Clock
      INTEGER            ,INTENT(OUT)   :: RC_INIT                         !<-- Error return code
!
!-----------------------------------------------------------------------
!
!---------------------
!***  Local Variables
!---------------------
!
      TYPE(ESMF_TimeInterval) :: RUNDURATION
!
      CHARACTER(20) :: PELAB
!
      CHARACTER(ESMF_MAXSTR),DIMENSION(:),ALLOCATABLE :: EARTH_COMP_NAME &  !<-- Names of each member's EARTH component
                                                        ,IMP_EARTH_NAME  &  !<-- Import state name of the EARTH components
                                                        ,EXP_EARTH_NAME     !<-- Export state name of the EARTH components
!
      INTEGER :: I,IJ,J,RC
!
      INTEGER :: MYPE_GLOBAL                                            &
                ,NHOURS_FCST                                            &
                ,NSECONDS_FCST                                          &
                ,PE_MAX                                                 &
                ,TASKS
!
      INTEGER,DIMENSION(:,:),ALLOCATABLE :: PETLIST                        !<-- Task list for each ensemble member
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC     =ESMF_SUCCESS
      RC_INIT=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  For the moment, use a direct copy of the main Clock within
!***  the NEMS component.
!-----------------------------------------------------------------------
!
      CLOCK_NEMS=CLOCK_MAIN
!
!-----------------------------------------------------------------------
!***  What is the start time on the NEMS clock?
!-----------------------------------------------------------------------

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK = "Extract the start time of the NEMS clock"
!     CALL ESMF_LogWrite(MESSAGE_CHECK, ESMF_LOG_INFO, rc = RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

      CALL ESMF_ClockGet(clock     = CLOCK_NEMS                         &
                        ,startTime = STARTTIME                          &
                        ,rc = RC)

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Allocate the NEMS component's internal state, point at it,
!***  and attach it to the NEMS component.
!-----------------------------------------------------------------------
!
      ALLOCATE(NEMS_INT_STATE,stat=RC)
      wrap%NEMS_INT_STATE=>NEMS_INT_STATE
!
      CALL ESMF_GridCompSetInternalState(NEMS_GRID_COMP                 &  !<--The NEMS component
                                        ,WRAP                           &  !<-- Pointer to the NEMS internal state
                                        ,RC)     
!
!-----------------------------------------------------------------------
!***  Get the global VM (Virtual Machine).
!***  Obtain the total task count and the local task ID.
!-----------------------------------------------------------------------
!
      CALL ESMF_VMGetGlobal(vm = VM_GLOBAL                              &  !<-- The ESMF global Virtual Machine
                           ,rc = RC)
!
      CALL ESMF_VmGet(vm       = VM_GLOBAL                              &  !<-- The ESMF global Virtual Machine
                     ,pecount  = TASKS                                  &  !<-- Total # of MPI tasks
                     ,localpet = MYPE_GLOBAL                            &  !<-- This task's global rank
                     ,rc       = RC)
!
!-----------------------------------------------------------------------
!***  Create and load the Configure object which will hold the contents
!***  of the NEMS configure file.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create/Load the NEMS Configure Object"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CF_NEMS=ESMF_ConfigCreate(rc=RC)
!
      CALL ESMF_ConfigLoadFile(config   = CF_NEMS                       &  !<-- The configure object
                              ,filename = 'model_configure'             &  !<-- The name of the configure file
                              ,rc       = RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!-----------------------------------------------------------------------
!***  Get the ensemble stochastic coupling flag from the config file.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract the Ensemble Stochastic Coupling Flag from Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config = CF_NEMS                     &  !<-- The NEMS configure object
                                  ,value  = ENS_SPS                     &  !<-- Value of control flag for 
                                                                           !    stochastic perturbation scheme
                                  ,label  = 'ENS_SPS:'                  &  !<-- Flag's label in configure file
                                  ,rc     = RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Reset the run duration in the NEMS Clock for the ensemble 
!***  stochastic perturbation case.
!-----------------------------------------------------------------------
!
      IF(ENS_SPS) THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="NEMS: Extract Forecast Length from Config File"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=CF_NEMS                     &
                                    ,value =NHOURS_FCST                 &
                                    ,label ='nhours_fcst1:'             &
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        NSECONDS_FCST=NHOURS_FCST*3600                                     !<-- The forecast length (sec) (Integer)
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="NEMS: Set the Forecast Length"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_TimeIntervalSet(timeinterval=RUNDURATION              &  !<-- The forecast length (s) (ESMF)
                                 ,s           =NSECONDS_FCST            &  !<-- The forecast length (s) (Integer)
                                 ,rc          =RC)
!
        CALL ESMF_ClockSet(clock       = CLOCK_NEMS                     &  !<-- The NEMS Clock
                          ,runDuration = RUNDURATION                    &  !<-- The forecast length (s) (ESMF)
                          ,rc          = RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      END IF
!
!-----------------------------------------------------------------------
!***  Extract the total number of EARTH ensemble members.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract the Total Number of the EARTH Members from Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config = CF_NEMS                     &  !<-- The NEMS configure object
                                  ,value  = TOTAL_MEMBER                &  !<-- Total # of ensemble members
                                  ,label  = 'total_member:'             &  !<-- Flag in configure file 
                                  ,rc     = RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Allocate a standard set of arrays for each ensemble member.
!-----------------------------------------------------------------------
!
      ALLOCATE(EARTH_GRID_COMP (TOTAL_MEMBER))
      ALLOCATE(EARTH_IMP_STATE (TOTAL_MEMBER))
      ALLOCATE(EARTH_EXP_STATE (TOTAL_MEMBER))
      ALLOCATE(EARTH_COMP_NAME (TOTAL_MEMBER))
      ALLOCATE(IMP_EARTH_NAME  (TOTAL_MEMBER))
      ALLOCATE(EXP_EARTH_NAME  (TOTAL_MEMBER))
      ALLOCATE(PE_MEMBER       (TOTAL_MEMBER))
!
!-----------------------------------------------------------------------
!***  For each member create the names of the EARTH components and 
!***  ESMF states then fill in the task information.
!-----------------------------------------------------------------------
!
      PE_MEMBER = 0
!
!JR Talked with Tom Black about this: all the models have no entries for PE_MEMBER in their
!JR config file. And, there's no check on the return code.
      DO I = 1, TOTAL_MEMBER
!
        WRITE(EARTH_COMP_NAME(I), '("EARTH grid component", I2.2)') I
        WRITE(IMP_EARTH_NAME (I), '("EARTH import state",   I2.2)') I
        WRITE(EXP_EARTH_NAME (I), '("EARTH export state",   I2.2)') I
!
        WRITE(PELAB, '("PE_MEMBER", I2.2, ":")') I
! 
!JR Comment out the call so that we can get esmf to halt on the first error. Apparently
!JR all the NEMS models fail on this call.
!JR        CALL ESMF_ConfigGetAttribute(config = CF_NEMS                   &
!JR                                    ,value  = PE_MEMBER(I)              &
!JR                                    ,label  = PELAB                     &
!JR                                    ,rc     = RC)
!
!JR        IF(PE_MEMBER(I) == 0) PE_MEMBER(i) = TASKS / TOTAL_MEMBER
!
      END DO
!JR Just set PE_MEMBER(1) manually
      pe_member(1) = tasks / total_member
!
      PE_MAX = 1
!
      DO I = 1, TOTAL_MEMBER
        PE_MAX = MAX(PE_MAX, PE_MEMBER(I))
      END DO
!
!-----------------------------------------------------------------------
!***  Set up the PE list.
!-----------------------------------------------------------------------
!
      ALLOCATE(PETLIST(1:PE_MAX, 1:TOTAL_MEMBER))
!
      IJ = 0
!
      DO J = 1, TOTAL_MEMBER
!
        DO I = 1, PE_MEMBER(J)
          PETLIST(I, J) = IJ
!
          IF(MYPE_GLOBAL == IJ) THEN
            MEMBER_ID = J
          END IF
!
          IJ = IJ + 1
        END DO
!
      END DO
!
!-----------------------------------------------------------------------
!***  Get the clock parameters for the ensemble stochastic perturbation 
!***  cycles from the config file.
!-----------------------------------------------------------------------
!
      IF(ENS_SPS) THEN 
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract the Ensemble Clock Parameters from Config File"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config = CF_NEMS                   &
                                    ,value  = hh_increase               &
                                    ,label  = 'HH_INCREASE:'            &
                                    ,rc     = RC)
!
        CALL ESMF_ConfigGetAttribute(config = CF_NEMS                   &
                                    ,value  = hh_start                  &
                                    ,label  = 'HH_START:'               &
                                    ,rc     = RC)
!
        CALL ESMF_ConfigGetAttribute(config = CF_NEMS                   &
                                    ,value  = hh_final                  &
                                    ,label  = 'HH_FINAL:'               &
                                    ,rc     = RC)
!
        NUMBER_START = HH_START / HH_INCREASE + 1
        NUMBER_FINAL = HH_FINAL / HH_INCREASE - 1
!
        PRINT *, 'ENS Coup ', NUMBER_START, NUMBER_FINAL                &
               , HH_START, HH_FINAL, HH_INCREASE
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      END IF
!
!-----------------------------------------------------------------------
!***  Create the EARTH grid components.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create EARTH grid Components"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      DO I = 1, TOTAL_MEMBER
        EARTH_GRID_COMP(i) = ESMF_GridCompCreate (                      &
                             name         = EARTH_COMP_NAME(I)          &  !<-- Name of element I of the EARTH component array
                            ,petlist      = PETLIST(1:PE_MEMBER(I), I)  &  !<-- Element I's PE list
                            ,config       = CF_NEMS                     &  !<-- Associate the NEMS config object with this element
                            ,rc           = RC)
      END DO
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Create the Ensemble Coupler component inside of
!***  the NEMS internal state.
!-----------------------------------------------------------------------
!
      IF(ENS_SPS) THEN
        ENS_CPL_COMP=ESMF_CplCompCreate(name = "ENS Cpl component"      &  
                                       ,rc   = RC)
      END IF
!
!-----------------------------------------------------------------------
!***  Register the Initialize, Run, and Finalize routines of
!***  each element in the EARTH component array.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Register EARTH Init, Run, Finalize"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      DO I = 1, TOTAL_MEMBER
#ifdef ESMF_3
        CALL ESMF_GridCompSetServices(EARTH_GRID_COMP(I)                &  !<-- The EARTH gridded components
                                     ,EARTH_REGISTER                    &  !<-- User's name for the Register routine
                                     ,RC)
#else
        CALL ESMF_GridCompSetServices(EARTH_GRID_COMP(I)                &  !<-- The EARTH gridded components
                                     ,EARTH_REGISTER                    &  !<-- User's name for the Register routine
                                     ,rc=RC)
#endif
      END DO
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Register the Initialize, Run, and Finalize routines of
!***  the Ensemble Coupler component.
!-----------------------------------------------------------------------
!
      IF(ENS_SPS) THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Register Ensemble Coupler Init, Run, Finalize"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifdef ESMF_3
        CALL ESMF_CplCompSetServices(ENS_CPL_COMP                       &
                                    ,ENS_CplCompSetServices             &  !<-- The user's name for the Register routine
                                    ,RC)
#else
        CALL ESMF_CplCompSetServices(ENS_CPL_COMP                       &
                                    ,ENS_CplCompSetServices             &  !<-- The user's name for the Register routine
                                    ,rc=RC)
#endif
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      END IF
!
!-----------------------------------------------------------------------
!***  Create the EARTH import and export states.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create the EARTH import states"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      DO I = 1,TOTAL_MEMBER
#ifdef ESMF_520rbs
        EARTH_IMP_STATE(I) = ESMF_StateCreate(                          &
                                              NAME = IMP_EARTH_NAME(I)  &
                                        ,statetype = ESMF_STATE_IMPORT  &
                                        ,rc        = RC)
#else
        EARTH_IMP_STATE(I) = ESMF_StateCreate(                          &
                                         STATENAME = IMP_EARTH_NAME(I)  &
                                        ,statetype = ESMF_STATE_IMPORT  &
                                        ,rc        = RC)
#endif
      END DO
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create the EARTH export states"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      DO I = 1,TOTAL_MEMBER
#ifdef ESMF_520rbs
        EARTH_EXP_STATE(I) = ESMF_StateCreate(                          &
                                              NAME = EXP_EARTH_NAME(I)  &
                                        ,statetype = ESMF_STATE_EXPORT  &
                                        ,rc        = RC)
#else
        EARTH_EXP_STATE(I) = ESMF_StateCreate(                          &
                                         STATENAME = EXP_EARTH_NAME(I)  &
                                        ,statetype = ESMF_STATE_EXPORT  &
                                        ,rc        = RC)
#endif
      END DO
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Create the Ensemble Coupler import and export states.
!-----------------------------------------------------------------------
!
      IF(ENS_SPS) THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Create the Ensemble Coupler import state"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifdef ESMF_520rbs
        ENS_CPL_IMP_STATE=ESMF_StateCreate(     NAME = "ENS_CPL_Import"  &
                                          ,statetype = ESMF_STATE_IMPORT &
                                          ,rc        = RC)
#else
        ENS_CPL_IMP_STATE=ESMF_StateCreate(STATENAME = "ENS_CPL_Import"  &
                                          ,statetype = ESMF_STATE_IMPORT &
                                          ,rc        = RC)
#endif
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Create the Ensemble Coupler export state"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifdef ESMF_520rbs
        ENS_CPL_EXP_STATE=ESMF_StateCreate(     NAME = "ENS_CPL_Export"  &
                                          ,statetype = ESMF_STATE_EXPORT &
                                          ,rc        = RC)
#else
        ENS_CPL_EXP_STATE=ESMF_StateCreate(STATENAME = "ENS_CPL_Export"  &
                                          ,statetype = ESMF_STATE_EXPORT &
                                          ,rc        = RC)
#endif
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Nest the EARTH export/import states into the import/export states
!***  of the ensemble coupler.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK= "Add the EARTH states into the ENS_CPL states"
!       CALL ESMF_LogWrite(MESSAGE_CHECK, ESMF_LOG_INFO, rc = RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
! 
        DO I = 1, TOTAL_MEMBER
          IF(MEMBER_ID == I) THEN
            CALL ESMF_StateAdd(state       = ENS_CPL_IMP_STATE          &
                              ,nestedState = EARTH_EXP_STATE(I)         &
                              ,rc          = RC)
!
            CALL ESMF_StateAdd(state       = ENS_CPL_EXP_STATE          &
                              ,nestedState = EARTH_IMP_STATE(I)         &
                              ,rc          = RC)
          END IF
        END DO
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      END IF
!
!-----------------------------------------------------------------------
!***  Execute the Initialize step of each element of the EARTH
!***  component array.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Execute the Initialize step of the EARTH component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      DO I = 1, TOTAL_MEMBER
!
        IF(MEMBER_ID == I) THEN
          CALL ESMF_GridCompInitialize(gridcomp    = EARTH_GRID_COMP(I)  &
                                      ,importState = EARTH_IMP_STATE(I)  &
                                      ,exportState = EARTH_EXP_STATE(I)  &
                                      ,clock       = CLOCK_NEMS          &
                                      ,phase       = ESMF_SINGLEPHASE    &
                                      ,rc          = RC)
        END IF
!
      END DO
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Execute the Initialize step of the Ensemble Coupler component.
!-----------------------------------------------------------------------
!
      IF(ENS_SPS) THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Execute the Initialize step of the Ensemble Coupler component"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_CplCompInitialize(cplcomp    =ENS_CPL_COMP            &
                                   ,importState=ENS_CPL_IMP_STATE       &
                                   ,exportState=ENS_CPL_EXP_STATE       &
                                   ,clock      =CLOCK_NEMS              &
                                   ,phase      =ESMF_SINGLEPHASE        &
                                   ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      END IF
!
!-----------------------------------------------------------------------
!
      IF(RC_INIT==ESMF_SUCCESS)THEN
!       WRITE(0,*)' NEMS_INITIALIZE succeeded'
      ELSE
        WRITE(0,*)' NEMS_INITIALIZE failed  RC_INIT=',RC_INIT
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE NEMS_INITIALIZE
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE NEMS_RUN(NEMS_GRID_COMP                                &
                         ,IMP_STATE                                     &
                         ,EXP_STATE                                     &
                         ,CLOCK_MAIN                                    &
                         ,RC_RUN)
!
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(ESMF_GridComp)               :: NEMS_GRID_COMP                  !<-- The NEMS component
      TYPE(ESMF_State)                  :: IMP_STATE                       !<-- The NEMS import state
      TYPE(ESMF_State)                  :: EXP_STATE                       !<-- The NEMS export state
      TYPE(ESMF_Clock)                  :: CLOCK_MAIN                      !<-- The main Clock
      INTEGER            ,INTENT(OUT)   :: RC_RUN                          !<-- Error return code
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER :: HH,I,J,RC
!
      TYPE(ESMF_Time) :: CURRTIME
!
      TYPE(ESMF_TimeInterval) :: RUNDURATION
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  For the moment, use a direct copy of the main Clock within
!***  the NEMS component.
!-----------------------------------------------------------------------
!
      CLOCK_NEMS=CLOCK_MAIN
!
!-----------------------------------------------------------------------
!***  Execute the Run step of each element in the EARTH component
!***  array.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Execute the Run step of the EARTH component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      DO I = 1, TOTAL_MEMBER
!
        IF(MEMBER_ID == I) THEN
          CALL ESMF_GridCompRun(gridcomp    = EARTH_GRID_COMP(I)        &
                               ,importState = EARTH_IMP_STATE(I)        &
                               ,exportState = EARTH_EXP_STATE(I)        &
                               ,clock       = CLOCK_NEMS                &
                               ,phase       = ESMF_SINGLEPHASE          &
                               ,rc          = RC)
        END IF
!
      END DO
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Execute the Run step of the Ensemble Coupler component.
!-----------------------------------------------------------------------
!
      sps: IF(ENS_SPS) THEN
!
!-----------------------------------------------------------------------
!
        DO I = NUMBER_START, NUMBER_FINAL
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Execute the Run step of the Ensemble Coupler component"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_CplCompRun(cplcomp    =ENS_CPL_COMP                 &
                              ,importState=ENS_CPL_IMP_STATE            &
                              ,exportState=ENS_CPL_EXP_STATE            &
                              ,clock      =CLOCK_NEMS                   &
                              ,phase      =ESMF_SINGLEPHASE             &
                              ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_VMBarrier(vm = VM_GLOBAL, rc = RC)
!
!-----------------------------------------------------------------------
!***  Adjust the ESMF clock for the next run cycle.
!-----------------------------------------------------------------------

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK = "Update the current time of the NEMS clock"
!         CALL ESMF_LogWrite(MESSAGE_CHECK, ESMF_LOG_INFO, rc = RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_ClockGet(clock       = CLOCK_NEMS                   &
                            ,runDuration = RUNDURATION                  &
                            ,rc          = RC)
!
          CURRTIME = STARTTIME + RUNDURATION
!
          CALL ESMF_ClockSet(clock    = CLOCK_NEMS                      &
                            ,currTime = CURRTIME                        &
                            ,rc       = RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC, MESSAGE_CHECK, RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK = "Adjust clock - add one more cycle run duration"
!         CALL ESMF_LogWrite(MESSAGE_CHECK, ESMF_LOG_INFO, rc = RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_TimeIntervalGet(timeInterval = RUNDURATION          &
                                   ,h            = HH                   &
                                   ,rc           = RC)
!
          HH = HH + HH_INCREASE
!
          CALL ESMF_TimeIntervalSet(timeInterval = RUNDURATION          &
                                   ,h            = hh                   &
                                   ,rc           = RC)
!
          CALL ESMF_ClockSet(clock       = CLOCK_NEMS                   &
                            ,runDuration = RUNDURATION                  &
                            ,rc = RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC, MESSAGE_CHECK, RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Execute the Run step of each element in the EARTH component
!***  array.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Execute the Run step of the EARTH components"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          DO J = 1, TOTAL_MEMBER
!
            IF(MEMBER_ID == J) THEN
                CALL ESMF_GridCompRun(gridcomp    = EARTH_GRID_COMP(J)  &
                                     ,importState = EARTH_IMP_STATE(J)  &
                                     ,exportState = EARTH_EXP_STATE(J)  &
                                     ,clock       = CLOCK_NEMS          &
                                     ,phase       = ESMF_SINGLEPHASE    &
                                     ,rc          = RC)
            END IF
!
          END DO
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          PRINT*, 'Complete EARTH Run Cycle ', I + 1
        END DO
!
!-----------------------------------------------------------------------
!
      ELSE sps
!
         CALL ESMF_ClockGet(clock       = CLOCK_NEMS                    &
                           ,runDuration = RUNDURATION                   &
                           ,rc          = RC)
!
!-----------------------------------------------------------------------
!
      END IF sps
!
!-----------------------------------------------------------------------
!
      IF(RC_RUN==ESMF_SUCCESS)THEN
!       WRITE(0,*)' NEMS_RUN succeeded'
      ELSE
        WRITE(0,*)' NEMS_RUN failed  RC_RUN=',RC_RUN
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE NEMS_RUN
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE NEMS_FINALIZE(NEMS_GRID_COMP                           &
                              ,IMP_STATE                                &
                              ,EXP_STATE                                &
                              ,CLOCK_MAIN                               &
                              ,RC_FINALIZE)
!
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(ESMF_GridComp)               :: NEMS_GRID_COMP                  !<-- The NEMS component
      TYPE(ESMF_State)                  :: IMP_STATE                       !<-- The NEMS import state
      TYPE(ESMF_State)                  :: EXP_STATE                       !<-- The NEMS export state
      TYPE(ESMF_Clock)                  :: CLOCK_MAIN                      !<-- The main Clock
      INTEGER            ,INTENT(OUT)   :: RC_FINALIZE                     !<-- Error return code
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER :: I,RC
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Execute the Finalize step of each element of the 
!***  EARTH component array.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Execute the Finalize step of the EARTH component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      DO I = 1, TOTAL_MEMBER
!
        IF(MEMBER_ID == I) THEN
          CALL ESMF_GridCompFinalize(gridcomp    = EARTH_GRID_COMP(I)   &
                                    ,importState = EARTH_IMP_STATE(I)   &
                                    ,exportState = EARTH_EXP_STATE(I)   &
                                    ,clock       = CLOCK_NEMS           &
                                    ,phase       = ESMF_SINGLEPHASE     &
                                    ,rc          = RC)
        END IF
!
      END DO
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINALIZE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
      IF(ENS_SPS) THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Execute the Finalize step of the Ensemble Coupler component"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_CplCompFinalize(cplcomp    =ENS_CPL_COMP              &
                                 ,importState=ENS_CPL_IMP_STATE         &
                                 ,exportState=ENS_CPL_EXP_STATE         &
                                 ,clock      =CLOCK_NEMS                &
                                 ,phase      =ESMF_SINGLEPHASE          &
                                 ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINALIZE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      END IF
!
!-----------------------------------------------------------------------
!
      IF(RC_FINALIZE==ESMF_SUCCESS)THEN
!       WRITE(0,*)' NEMS_FINALIZE succeeded'
      ELSE
        WRITE(0,*)' NEMS_FINALIZE failed  RC_FINALIZE=',RC_FINALIZE
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE NEMS_FINALIZE
!
!-----------------------------------------------------------------------
!
      END MODULE module_NEMS_GRID_COMP
!
!-----------------------------------------------------------------------
