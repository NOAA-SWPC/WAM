!#include "../../ESMFVersionDefine.h"

!-----------------------------------------------------------------------
!
      MODULE MODULE_GOCART_ROUTINES
!
!-----------------------------------------------------------------------
!
!*** THIS MODULE CONTAINS THE STUB ROUTINES for 
!*** THE AEROSOL MODULE (GOCART)
!***
!*** THE SETUP AND INIT ROUTINES ARE CALLED FROM GFS_ATM_INIT 
!*** THE INTEGRATE ROUTINE IS CALLED FROM GFS_INTEGRATE
!*** 
!
!-----------------------------------------------------------------------
!
! PROGRAM HISTORY LOG:
!   2010-03-05  Wang  - set up stub for MODULE_GOCART_ROUTINES
!
!-----------------------------------------------------------------------
!
      USE ESMF
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
      TYPE(ESMF_GridComp),INTENT(INOUT) :: GC_GFS_CHEM                     !<-- The gocart gridded component
      TYPE(ESMF_State),   INTENT(INOUT) :: IMP_GFS_CHEM                    !<-- The gocart import state
      TYPE(ESMF_State),   INTENT(INOUT) :: EXP_GFS_CHEM                    !<-- The gocart export state
      TYPE(ESMF_CplComp), INTENT(INOUT) :: GC_PHY2CHEM_CPL                 !<-- Phy to Chem coupler component
      TYPE(ESMF_CplComp), INTENT(INOUT) :: GC_CHEM2PHY_CPL                 !<-- Chem to Phy coupler component
      TYPE(ESMF_Logical), INTENT(INOUT) :: CHEMISTRY_ON                    !<-- The option to activate gocart
      INTEGER,            INTENT(IN)    :: MYPE                            !<-- MPI task ID
      INTEGER,OPTIONAL,   INTENT(OUT)   :: RC_SETUP                        !<-- Return code for the SETUP step
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC_SETUP      =ESMF_SUCCESS
!     write(0,*)' Initialize with gocart coupling '
!
!
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
      INTEGER,            INTENT(IN)    :: MYPE
      TYPE(ESMF_GridComp),INTENT(INOUT) :: GC_GFS_CHEM                     !<-- The gocart gridded component
      TYPE(ESMF_State),   INTENT(INOUT) :: EXP_GFS_PHY                     !<-- The physics component's export
      TYPE(ESMF_State),   INTENT(INOUT) :: IMP_GFS_CHEM                    !<-- The chemistry component' import
      TYPE(ESMF_State),   INTENT(INOUT) :: EXP_GFS_CHEM                    !<-- The chemistry component's export
      TYPE(ESMF_CplComp), INTENT(INOUT) :: GC_PHY2CHEM_CPL                 !<-- The Phys to Chem coupler component
      TYPE(ESMF_CplComp), INTENT(INOUT) :: GC_CHEM2PHY_CPL                 !<-- The Chem to Phys coupler component
      TYPE(ESMF_Clock)   ,INTENT(INOUT) :: CLOCK_ATM                       !<-- The ESMF Clock from the ATM Driver component
      INTEGER,OPTIONAL,   INTENT(OUT)   :: RC_INIT                         !<-- Return code for the INIT step
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC_INIT = ESMF_SUCCESS       ! Error signal variable

!      WRITE(0,*)'GOCART_INIT STEP SUCCEEDED'
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
!   2011-10-30  Lu    MYPE added
!
!-----------------------------------------------------------------------

      INTEGER,INTENT(IN)               :: MYPE
      TYPE(ESMF_GridComp),INTENT(INOUT):: GC_GFS_CHEM                 !<-- The GOCART grid component
      TYPE(ESMF_CplComp),INTENT(INOUT) :: GC_PHY2CHEM_CPL             !<-- The Phy-to-Chem coupler component
      TYPE(ESMF_CplComp),INTENT(INOUT) :: GC_CHEM2PHY_CPL             !<-- The Chem-to-Phy coupler component
!
      TYPE(ESMF_State),INTENT(INOUT)   :: EXP_GFS_PHY,             &  !<-- The export states for Physics component
                                          IMP_GFS_CHEM,EXP_GFS_CHEM   !<-- The imp/exp states for Chemistry component
      TYPE(ESMF_Clock),INTENT(INOUT)   :: CLOCK_ATM                   !<-- The ATM Component's ESMF Clock

      INTEGER,INTENT(OUT) :: RC_LOOP                                  !<-- Return code

!-----------------------------------------------------------------------
!
      RC_LOOP = ESMF_SUCCESS     
!     WRITE(0,*)'GOCART_LOOP STEP SUCCEEDED'
!
      END SUBROUTINE GOCART_INTEGRATE
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!

      END MODULE MODULE_GOCART_ROUTINES

!-----------------------------------------------------------------------
