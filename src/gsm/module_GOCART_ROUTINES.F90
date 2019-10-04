#include "../../ESMFVersionDefine.h"

!-----------------------------------------------------------------------
!
      MODULE MODULE_GOCART_ROUTINES
!
!-----------------------------------------------------------------------
!
!*** THIS MODULE CONTAINS THE ROUTINES TO SETUP, INITIALIZE, AND RUN
!*** THE AEROSOL MODULE (GOCART)
!***
!*** THE SETUP AND INIT ROUTINES ARE CALLED FROM GFS_ATM_INIT 
!*** THE INTEGRATE ROUTINE IS CALLED FROM GFS_INTEGRATE
!*** 
!
!-----------------------------------------------------------------------
!
! PROGRAM HISTORY LOG:
!   2010-03-05  Lu    - Create the module
!   2010-03-05  Lu    - Change chemistry_on from out to inout
!   2010-03-06  Lu    - 2-phased GOCART initialization
!   2010-08-17  Lu    - Call Chem_RegistryPrint only for master PE
!   2011-05-11  Theurich & Yang - Modified for using the ESMF 5.2.0r_beta_snapshot_07.
!   2011-10-01  Wang/Lu - exit GOCART_INIT and GOCART_INTEGRATE for IO PE
!   2014-11-20  Wang  - Update to esmf6.3.0 
!-----------------------------------------------------------------------
!
      USE ESMF
!
      USE MODULE_ERR_MSG,               ONLY: ERR_MSG,MESSAGE_CHECK

      USE MODULE_GFS_MPI_DEF,           ONLY: PETLIST_FCST 

      USE gfs_physics_grid_create_mod,  ONLY: mgrid

!
!-----------------------------------------------------------------------
!***  LIST MODULES FOR GSFC CHEMISTRY PACKAGE
!-----------------------------------------------------------------------
!
!jw      USE GOCART_GridCompMod    , ONLY: GOCART_SETSERVICES
      USE GOCART_GridCompMod    , ONLY: GOCART_SETSERVICES => SETSERVICES 
!
      USE ATMOS_PHY_CHEM_CPL_COMP_MOD, ONLY: PHY2CHEM_SETSERVICES => SETSERVICES
      USE ATMOS_CHEM_PHY_CPL_COMP_MOD, ONLY: CHEM2PHY_SETSERVICES => SETSERVICES
!
      USE Chem_RegistryMod
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
      PRIVATE
!
      PUBLIC :: GOCART_SETUP, GOCART_INIT, GOCART_INTEGRATE
!
      CONTAINS


!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------

      SUBROUTINE GOCART_SETUP ( GC_GFS_CHEM                            &
                               ,IMP_GFS_CHEM                           &
                               ,EXP_GFS_CHEM                           &
                               ,GC_PHY2CHEM_CPL                        &
                               ,GC_CHEM2PHY_CPL                        &
                               ,CHEMISTRY_ON                           &
                               ,MYPE                                   &
                               ,RC_SETUP                               &
                                )
!
!-----------------------------------------------------------------------
!***  THIS SUBROUTINE PERFORMS THE SETUP STEP OF THE GOCART GRID
!***  COMPONENT AND THE ASSOCIATED COUPLER COMPONENTS
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(ESMF_GridComp)               :: GC_GFS_CHEM                     !<-- The gocart gridded component
      TYPE(ESMF_State)                  :: IMP_GFS_CHEM                    !<-- The gocart import state
      TYPE(ESMF_State)                  :: EXP_GFS_CHEM                    !<-- The gocart export state
      TYPE(ESMF_CplComp)                :: GC_PHY2CHEM_CPL                 !<-- Phy to Chem coupler component
      TYPE(ESMF_CplComp)                :: GC_CHEM2PHY_CPL                 !<-- Chem to Phy coupler component
      TYPE(ESMF_Logical)                :: CHEMISTRY_ON                    !<-- The option to activate gocart
      INTEGER,            INTENT(IN)    :: MYPE                            !<-- MPI task ID
      INTEGER,OPTIONAL,   INTENT(OUT)   :: RC_SETUP                        !<-- Return code for the SETUP step
      real gocart_setservices1
!
!---------------------
!***  Local Variables
!---------------------
!
      TYPE(Chem_Registry),SAVE :: REG                                      !<-- The GOCART Chem_Registry

      INTEGER                  :: RC, IERR
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC      =ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  Determine whether GOCART is active
!-----------------------------------------------------------------------
!
      REG = Chem_RegistryCreate ( IERR )                    !<-- read Chem_Registry

      IF (MYPE==0)   CALL Chem_RegistryPrint ( REG )

      IF(REG%doing_gocart)THEN                              !<-- GOCART => Chemistry on

         CHEMISTRY_ON=ESMF_True                                         
         write(0,*)' Initialize with gocart coupling '

      ELSE                                                  !<-- no  GOCART => Chemistry off

         CHEMISTRY_ON=ESMF_False                                     
         write(0,*)' Initialize without gocart coupling '    
      ENDIF

      CALL Chem_RegistryDestroy ( REG, IERR ) 
!
!
!------------------------------------------------------------------------
!***  Create empty Import and Export states 
!------------------------------------------------------------------------
!
      MESSAGE_CHECK="Create Empty Import/Export States for GFS Chemistry"
!
      IMP_GFS_CHEM =ESMF_StateCreate(name="chemistry import"         &
                                  ,stateintent=ESMF_STATEINTENT_IMPORT    &
                                  ,rc         =RC)
!
      EXP_GFS_CHEM =ESMF_StateCreate(name="chemistry export"         &
                                  ,stateintent=ESMF_STATEINTENT_EXPORT    &
                                  ,rc         =RC)

      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_SETUP)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  Create the chemistry gridded component if chemistry is turned on.
!***  Register the Initialize, Run, and Finalize steps for it.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

      IF (CHEMISTRY_ON==ESMF_True)THEN
!
!-------------------------------
!***  Create Chemistry component
!-------------------------------
!
        MESSAGE_CHECK="Create the GFS Chemistry Component"
!
        GC_GFS_CHEM =ESMF_GridCompCreate(name    ="chemistry component"   &
                                     ,ConfigFile ='MAPL.rc'               &
                                     ,petList   =PETLIST_FCST             &
                                     ,rc        =RC)
        write(0,*)'in GOCART_SETUP after chem comp created, petlist_fcst=',petlist_fcst
!
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_SETUP)
!
!-------------------------------------------------
!***  Register the Init, Run, and Finalize steps.
!-------------------------------------------------
!
        MESSAGE_CHECK="Register Chemistry Init, Run, Finalize"
!
        CALL ESMF_GridCompSetServices(GC_GFS_CHEM                        &
                                     ,GOCART_SETSERVICES                &
                                     ,RC=RC)
!
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_SETUP)
!
!
!----------------------------
!***  Create Phy-Chem Coupler
!----------------------------
!
        MESSAGE_CHECK="Create the GFS Phy2Chem Coupler Component"
!
        GC_PHY2CHEM_CPL=ESMF_CplCompCreate(name   ="phy2chem component" &
                                   ,petList=PETLIST_FCST                &
                                   ,rc     =RC)
!
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_SETUP)
!
!-------------------------------------------------
!***  Register the phy-to-chem coupler Init and Run steps.
!-------------------------------------------------
!
        MESSAGE_CHECK="Register the Phy2Chem Coupler's Init and Run"
!
        CALL ESMF_CplCompSetServices(GC_PHY2CHEM_CPL             &  !<-- The GFS Phys-to-Chem coupler component
                                  ,PHY2CHEM_SETSERVICES          &  !<-- The user's subroutine name for Register
                                  ,rc=RC)
!
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_SETUP)
!
!----------------------------
!***  Create Chem-Phy Coupler
!----------------------------
!
        MESSAGE_CHECK="Create the GFS Chem2Phy Coupler Component"
!
        GC_CHEM2PHY_CPL=ESMF_CplCompCreate(name   ="chem2phy component" &
                                   ,petList=PETLIST_FCST                &
                                   ,rc     =RC)
!
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_SETUP)
!
!-------------------------------------------------
!***  Register the chem-to-phy coupler Run steps.
!-------------------------------------------------
!
        MESSAGE_CHECK="Register the Chem2Phy Coupler's Run"
!
        CALL ESMF_CplCompSetServices(GC_CHEM2PHY_CPL             &  !<-- The GFS Chem-to-Phys coupler component
                                  ,CHEM2PHY_SETSERVICES          &  !<-- The user's subroutine name for Register
                                  ,rc=RC)
!
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_SETUP)

      ENDIF 

!-----------------------------------------------------------------------
!***  THE FINAL ERROR SIGNAL INFORMATION.
!-----------------------------------------------------------------------
!
      IF(RC_SETUP==ESMF_SUCCESS)THEN
!       WRITE(0,*)'GOCART_SETUP STEP SUCCEEDED'
      ELSE
        WRITE(0,*)'GOCART_SETUP STEP FAILED RC_SETUP=',RC_SETUP
      ENDIF

      END SUBROUTINE GOCART_SETUP


!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE GOCART_INIT ( GC_GFS_CHEM                              &
                              ,EXP_GFS_PHY                              &
                              ,IMP_GFS_CHEM                             &
                              ,EXP_GFS_CHEM                             &
                              ,GC_PHY2CHEM_CPL                          &
                              ,GC_CHEM2PHY_CPL                          &
                              ,CLOCK_ATM                                &
                              ,MYPE                                     &
                              ,RC_INIT                                  &
                              )


!------------------------
!***  Argument variables
!------------------------
!
      TYPE(ESMF_GridComp)               :: GC_GFS_CHEM                     !<-- The gocart gridded component
      TYPE(ESMF_State)                  :: EXP_GFS_PHY                     !<-- The physics component's export
      TYPE(ESMF_State)                  :: IMP_GFS_CHEM                    !<-- The chemistry component' import
      TYPE(ESMF_State)                  :: EXP_GFS_CHEM                    !<-- The chemistry component's export
      TYPE(ESMF_CplComp)                :: GC_PHY2CHEM_CPL                 !<-- The Phys to Chem coupler component
      TYPE(ESMF_CplComp)                :: GC_CHEM2PHY_CPL                 !<-- The Chem to Phys coupler component
      TYPE(ESMF_Clock)                  :: CLOCK_ATM                       !<-- The ESMF Clock from the ATM Driver component
      INTEGER                           :: MYPE                            !<-- MPI task ID
      INTEGER,OPTIONAL,   INTENT(OUT)   :: RC_INIT                         !<-- Return code for the INIT step
!
!---------------------
!***  Local variables
!---------------------
!
      INTEGER :: RC, gridRank, I
      LOGICAL :: IOPE

      RC=ESMF_SUCCESS       ! Error signal variable

!-----------------------------------------------------------------------
!***  If io pe, return 
!-----------------------------------------------------------------------
!
      IOPE=.true.
      DO I=1,size(PETLIST_FCST)
         IF (PETLIST_FCST(I)== MYPE) IOPE=.false.
      ENDDO
      IF ( IOPE ) RETURN

!------------------------
!*** Attach the grid to GOCART Component
!------------------------
!
!     mGrid (3D Gaussian grid) is created in the PHY init step and will 
!     be used to set the grid in GOCART grid component
!
      MESSAGE_CHECK="Retrive gridRank from mGrid "

      call ESMF_GridGet ( grid = mGrid, dimCount=gridRank, rc=RC)
      print *,'GOCART_INIT mGrid grid info: gridRank = ',gridRank

      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
!!
      MESSAGE_CHECK="Set grid for GOCART Component"

      CALL ESMF_GridCompSet(gridcomp=GC_GFS_CHEM                     &
                           ,grid    =mGrid                           &
                           ,rc      =RC)

      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)

!------------------------
!*** Initialize GOCART Component (phase 1)
!------------------------
!
      MESSAGE_CHECK="Initialize GOCART Component (phase-1)"
!
      CALL ESMF_GridCompInitialize(gridcomp   =GC_GFS_CHEM           &
                                  ,importstate=IMP_GFS_CHEM          &
                                  ,exportstate=EXP_GFS_CHEM          &
                                  ,clock      =CLOCK_ATM             &
                                  ,rc         =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)

!
!------------------------
!*** Initialize Phy-Chem Coupler Component 
!------------------------
!
      MESSAGE_CHECK="Initialize Phys-Chem Coupler"

      CALL ESMF_CplCompInitialize(cplcomp    =GC_PHY2CHEM_CPL         &
                                 ,importstate=EXP_GFS_PHY             &
                                 ,exportstate=IMP_GFS_CHEM            &
                                 ,clock      =CLOCK_ATM               &
                                 ,rc         =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)


!------------------------
!*** Initialize GOCART Component (phase 2)
!------------------------
!
      MESSAGE_CHECK="Initialize GOCART Component (phase-2)"
!
      CALL ESMF_GridCompInitialize(gridcomp   =GC_GFS_CHEM           &
                                  ,importstate=IMP_GFS_CHEM          &
                                  ,exportstate=EXP_GFS_CHEM          &
                                  ,clock      =CLOCK_ATM             &
                                  ,phase      =2                     &
                                  ,rc         =RC)
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)

!-----------------------------------------------------------------------
!***  THE FINAL ERROR SIGNAL INFORMATION.
!-----------------------------------------------------------------------
!
      IF(RC_INIT==ESMF_SUCCESS)THEN
!       WRITE(0,*)'GOCART_INIT STEP SUCCEEDED'
      ELSE
        WRITE(0,*)'GOCART_INIT STEP FAILED RC_INIT=',RC_INIT
      ENDIF

!
      END SUBROUTINE GOCART_INIT

!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------

      SUBROUTINE GOCART_INTEGRATE(                                     &
                                  GC_GFS_CHEM,                         &
                                  GC_PHY2CHEM_CPL,                     &
                                  GC_CHEM2PHY_CPL,                     &
                                  EXP_GFS_PHY,                         &
                                  IMP_GFS_CHEM, EXP_GFS_CHEM,          &
                                  CLOCK_ATM, MYPE, RC_LOOP             )

!-----------------------------------------------------------------------
!
! PROGRAM HISTORY LOG:
!   2010-02-04  Lu    First crack.
!
!-----------------------------------------------------------------------

      INTEGER,INTENT(IN)       :: MYPE
      TYPE(ESMF_GridComp)      :: GC_GFS_CHEM                 !<-- The GOCART grid component
      TYPE(ESMF_CplComp)       :: GC_PHY2CHEM_CPL,        &   !<-- The Phy-to-Chem coupler component
                                  GC_CHEM2PHY_CPL             !<-- The Chem-to-Phy coupler component
!
      TYPE(ESMF_State)         :: EXP_GFS_PHY,             &  !<-- The export states for Physics component
                                  IMP_GFS_CHEM,EXP_GFS_CHEM   !<-- The imp/exp states for Chemistry component
      TYPE(ESMF_Clock)         :: CLOCK_ATM                   !<-- The ATM Component's ESMF Clock


      INTEGER,INTENT(OUT)      :: RC_LOOP                                  !<-- Return code

! Locals
      INTEGER                  :: I, RC=ESMF_SUCCESS  
      TYPE(ESMF_FieldBundle)   :: Bundle
      LOGICAL                  :: IOPE

!
!-----------------------------------------------------------------------
!***  If io pe, return
!-----------------------------------------------------------------------
!
       IOPE=.true.
       DO I=1,size(PETLIST_FCST)
         IF (PETLIST_FCST(I)== MYPE) IOPE=.false.
       ENDDO
       IF (IOPE ) RETURN
!-----------------------------------------------------------------------
!***  Couple Physics export state to Chemistry import State
!-----------------------------------------------------------------------

       MESSAGE_CHECK="GOCART_INTEGRATE: couple phy_exp to chem_imp"

       CALL ESMF_CplCompRun(cplcomp     = GC_PHY2CHEM_CPL            &
                            ,importstate= EXP_GFS_PHY                &
                            ,exportstate= IMP_GFS_CHEM               &
                            ,clock      = CLOCK_ATM                  &
                            ,rc         = RC)
!
       CALL ERR_MSG(RC,MESSAGE_CHECK,RC_LOOP)

!-----------------------------------------------------------------------
!***  Execute the Run step of the Chemistry Component
!-----------------------------------------------------------------------

       MESSAGE_CHECK="GOCART_INTEGRATE: execute chemistry"

       CALL ESMF_GridCompRun(gridcomp   =GC_GFS_CHEM                 &
                            ,importstate=IMP_GFS_CHEM                &
                            ,exportstate=EXP_GFS_CHEM                &
                            ,clock      =CLOCK_ATM                   &
                            ,rc         =RC)

       CALL ERR_MSG(RC,MESSAGE_CHECK,RC_LOOP)

!-----------------------------------------------------------------------
!***  Couple Chemistry export state to Physics export State
!-----------------------------------------------------------------------

       MESSAGE_CHECK="GOCART_INTEGRATE: couple chem_exp to phy_exp"

       CALL ESMF_CplCompRun(cplcomp     = GC_CHEM2PHY_CPL            &
                            ,importstate= EXP_GFS_CHEM               &
                            ,exportstate= EXP_GFS_PHY                &
                            ,clock      = CLOCK_ATM                  &
                            ,rc         = RC)

      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_LOOP)

!-----------------------------------------------------------------------
!***  THE FINAL ERROR SIGNAL INFORMATION.
!-----------------------------------------------------------------------
!
      IF(RC_LOOP==ESMF_SUCCESS)THEN
!       WRITE(0,*)'GOCART_LOOP STEP SUCCEEDED'
      ELSE
        WRITE(0,*)'GOCART_LOOP STEP FAILED RC_LOOP=',RC_LOOP
      ENDIF

!
!
      END SUBROUTINE GOCART_INTEGRATE
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!

      END MODULE MODULE_GOCART_ROUTINES

!-----------------------------------------------------------------------
