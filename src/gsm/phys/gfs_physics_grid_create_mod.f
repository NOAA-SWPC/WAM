!#include "../../../ESMFVersionDefine.h"

      module gfs_physics_grid_create_mod

!
!-------------------------------------------------------------------
! this code is used to create the esmf grids for the gfs esmf model.
! weiyu yang, 09/2005.
! updated by henry juang 04/2007
! updated by shrinivas moorthi for physics on 07/2007
! updated by henry juang 11/2007
! weiyu yang, 02/2008, updated to use the ESMF 3.1.0 library.
! Sarah Lu, 10/09/09, add gfs_physics_grid_create_Gauss3D routine 
!                     that creates mGrid (3D Gaussian grids)
! Weiyu yang, 05/16/2011, modified for using the ESMF 5.2.0r_beta_snapshot_07.
! Weiyu yang, 09/29/2011, modified for using the ESMF 5.2.0r library.
!-------------------------------------------------------------------
!
!!uses:
!
      use ESMF,                             ONLY: esmf_grid, esmf_vm,                &
                                                  ESMF_DistGrid, esmf_success,       &
                                                  ESMF_LogWrite, ESMF_LOGMSG_INFO,   &
                                                  ESMF_DistGridCreate,               &
                                                  ESMF_LogFoundError,                &
                                                  ESMF_FAILURE,                      &
                                                  ESMF_GridCreate,                   &
                                                  ESMF_GridCreate,                   &
                                                  ESMF_INDEX_DELOCAL,                &
                                                  ESMF_GridAddCoord,                 &
                                                  ESMF_GridGetCoord,                 &
                                                  ESMF_STAGGERLOC_CENTER,            &
                                                  ESMF_GridValidate,ESMF_MAXSTR,     &
                                                  ESMF_DELayout,ESMF_GridGet,        &
                                                  ESMF_DistGridGet,ESMF_DELayoutGet
      use gfs_physics_internal_state_mod,   ONLY: gfs_physics_internal_state
      use machine,                          ONLY: kind_rad       !added for mGrid

      implicit none

      type(esmf_grid), save :: grid0   ! the esmf grid type array. for the 
                                       ! gfs start date and time information.
      type(esmf_grid), save :: grid3   ! the esmf grid type array.
                                       ! for the single gaussian grid arrays.
      type(esmf_grid), save :: grid4   ! the esmf grid type array.
                                       ! for the multiple gaussian grid arrays.
      type(ESMF_Grid), save :: mGrid   ! Mid-layer 3D gaussian grid (im,jm,km)


      contains

!------------------------------------------------------------------------
      subroutine gfs_physics_grid_create_gauss(vm, int_state,  &
                                               DistGrid0, DistGrid3, DistGrid4, rc)
!
! this routine create Gaussian grid type for single and multiple levels
! grid 3 (single) and grid4(multiple)
!
      type(esmf_vm),                     intent(inout) :: vm   
      type(gfs_physics_internal_state),  intent(inout) :: int_state
      integer,                           intent(out)   :: rc

! The ESMF DistGrids
      TYPE(ESMF_DistGrid),               INTENT(inout) :: DistGrid0, DistGrid3, &
                                                          DistGrid4

      integer                           :: rc1, rcfinal
      integer,            dimension(2)  :: counts, arraystr, arrayend

      rc1     = esmf_success
      rcfinal = esmf_success

! create grid.
!=====================================================================
! set up parameter arrays for the esmf grid of the gaussian grid space.
! the first dimension is the multiple of latitudian and longitudian
! the second dimension is single level.
!----------------------------------------------------------------------
      counts(1)   = int_state%lonr*int_state%lats_node_r_max
      counts(1)   = counts(1) * int_state%nodes
      counts(2)   = 1
      arraystr(1) = 1
      arraystr(2) = 1
      arrayend(1) = counts(1)
      arrayend(2) = counts(2)

! Create the ESMF DistGrid3 using the 1-D default decomposition.
!---------------------------------------------------------------
      CALL ESMF_LogWrite("Create DistGrid3", ESMF_LOGMSG_INFO, rc = rc1)

      DistGrid3 = ESMF_DistGridCreate(arraystr, arrayend, rc = rc1)

      IF(ESMF_LogFoundError(rc1, msg="Create DistGrid3")) THEN
          rcfinal = ESMF_FAILURE
          PRINT*, 'Error Happened When Creating DistGrid3, rc = ', rc1
          rc1     = ESMF_SUCCESS
      END IF

! Create the ESMF grid3 based on the created ESMF DistGrid3 information.
! Grid3 is the grid for the Gaussian grid space.
!-----------------------------------------------------------------------
      CALL ESMF_LogWrite("create gfs_phy grid3", ESMF_LOGMSG_INFO, rc = rc1)

      grid3 = ESMF_GridCreate(name = "gfs_phy grid3", distgrid = DistGrid3, rc = rc1)

      IF(ESMF_LogFoundError(rc1, msg="Create Grid3")) THEN
          rcfinal = ESMF_FAILURE
          PRINT*, 'Error Happened When Creating Grid3, rc = ', rc1
          rc1     = ESMF_SUCCESS
      END IF

!=====================================================================
! set up parameter arrays for the esmf grid of the gaussian grid space.
! the first dimension is the multiple of latitudian and longitudian
! the second dimension is single level.
!--------------------------------------------------
      counts(1)   = int_state%lonr*int_state%lats_node_r_max
      counts(1)   = counts(1) * int_state%nodes
      counts(2)   = int_state%levs
      arraystr(1) = 1
      arraystr(2) = 1
      arrayend(1) = counts(1)
      arrayend(2) = counts(2)

! Create the ESMF DistGrid4 using the 1-D default decomposition.
!---------------------------------------------------------------
      CALL ESMF_LogWrite("Create DistGrid4", ESMF_LOGMSG_INFO, rc = rc1)

      DistGrid4 = ESMF_DistGridCreate(arraystr, arrayend, rc = rc1)

      IF(ESMF_LogFoundError(rc1, msg="Create DistGrid4")) THEN
          rcfinal = ESMF_FAILURE
          PRINT*, 'Error Happened When Creating DistGrid4, rc = ', rc1
          rc1     = ESMF_SUCCESS
      END IF

! Create the ESMF grid4 based on the created ESMF DistGrid4 information.
! Grid4 is the grid for the multiple level Gaussian grid space.
!-----------------------------------------------------------------------
      CALL ESMF_LogWrite("create gfs_phy grid4", ESMF_LOGMSG_INFO, rc = rc1)

      grid4 = ESMF_GridCreate(name = "gfs_phy grid4", distgrid = DistGrid4, rc = rc1)

      IF(ESMF_LogFoundError(rc1, msg="Create Grid4")) THEN
          rcfinal = ESMF_FAILURE
          PRINT*, 'Error Happened When Creating Grid4, rc = ', rc1
          rc1     = ESMF_SUCCESS
      END IF

! set up parameter arrays for the esmf grid used for the date and time
! information to run the gfs.  all processors contains the same five date
! and time valus.
!------------------------------------------------------------------------
      counts(1)   = int_state%nodes
      counts(2)   = 5
      arraystr(1) = 1
      arraystr(2) = 1
      arrayend(1) = counts(1)
      arrayend(2) = counts(2)

! Create the ESMF DistGrid0 using the 1-D default decomposition.
!---------------------------------------------------------------
      CALL ESMF_LogWrite("Create DistGrid0", ESMF_LOGMSG_INFO, rc = rc1)

      DistGrid0 = ESMF_DistGridCreate(arraystr, arrayend, rc = rc1)

      IF(ESMF_LogFoundError(rc1, msg="Create DistGrid0")) THEN
          rcfinal = ESMF_FAILURE
          PRINT*, 'Error Happened When Creating DistGrid0, rc = ', rc1
          rc1     = ESMF_SUCCESS
      END IF

! create the esmf grid for the date and time information.
!--------------------------------------------------------
      CALL ESMF_LogWrite("create create gfs_phy grid0", ESMF_LOGMSG_INFO, rc = rc1)

      grid0 = ESMF_GridCreate(name = "gfs_phy grid0", distgrid = DistGrid0, rc = rc1)

      IF(ESMF_LogFoundError(rc1, msg="Create Grid0")) THEN
          rcfinal = ESMF_FAILURE
          PRINT*, 'Error Happened When Creating Grid0, rc = ', rc1
          rc1     = ESMF_SUCCESS
      END IF

      if(rcfinal == esmf_success) then
          print*, "pass: gfs_physics_grid_create_gauss."
      else
          print*, "fail: gfs_physics_grid_create_gauss."
      end if

      rc = rcfinal

      end subroutine gfs_physics_grid_create_gauss

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP
!
!    !IROUTINE: gfs_physics_grid_create_Gauss3D
!
!    !INTERFACE:
!
      subroutine gfs_physics_grid_create_Gauss3D(vm, iState, distGrid, rc)

!     !ARGUMENTS:

      type(ESMF_VM),                     intent(in)  :: vm   
      type(gfs_physics_internal_state),  intent(in)  :: iState
      TYPE(ESMF_DistGrid),               intent(out) :: distGrid

      integer,                           intent(out)  :: rc

!     !DESCRIPTION:
!
!     This routine creates a 3D Gaussian grid on mid-layer vertical levels.
!     The actual distribution has already been determined and such information
!     is contained in the internal state *iState*.
!
!     Consistent with the current design, the actual grid is returned in 
!     module-scoped variable *mGrid*.
!
!     !REVISION HISTORY:
!      da Silva  05Feb2009  Based of similar MAPL routine
!      Sarah Lu  13Feb2009  Remove MAPL exception handling; specify coord info 
!      Sarah Lu  09Oct2009  Port from local branch to the trunk
!
!EOP
!                                      ---

      integer                 :: I1, IN, J1, JN, rc1, i, j

!     Local space for coordinate information
!     --------------------------------------
      real(kind=kind_rad), pointer :: centerX(:,:), centerY(:,:)

! initialize the error signal variables.                            
!---------------------------------------
      rc1 = esmf_success

      CALL ESMF_LogWrite('CreateGauss3D', ESMF_LOGMSG_INFO, rc = rc1)

!     Create grid with index-space information from internal state
!     ------------------------------------------------------------
      mGrid = ESMF_GridCreate( name='mGrid',                  &
                   countsPerDEDim1=(/iState%lonr/),           &
                   countsPerDEDim2=iState%lats_nodes_r_fix,   &
                   countsPerDEDim3=(/iState%levs/),           &
                   coordDep1 = (/1,2/),                       &
                   coordDep2 = (/1,2/),                       &
                   coordDep3 = (/3/),                         &
                   gridEdgeLWidth = (/0,0,0/),                &
                   gridEdgeUWidth = (/0,0,0/),                &
                   indexflag = ESMF_INDEX_DELOCAL,            &
                   rc = rc1 )

!  Add coordinate information
!  --------------------------
     call ESMF_GridAddCoord(mGrid, rc = rc1 )

!  Retrieve the coordinates so we can set them
!  -------------------------------------------
     call ESMF_GridGetCoord (mGrid, 1, localDE=0,               &
                           staggerloc=ESMF_STAGGERLOC_CENTER, &
                           farrayPtr=centerX, rc = rc1 )

     call ESMF_GridGetCoord (mGrid, 2, localDE=0,               &
                             staggerloc=ESMF_STAGGERLOC_CENTER, &
                             farrayPtr=centerY, rc = rc1 )
 
     call GridGetInterior_ (mGrid,i1,in,j1,jn)
  
      do i = 1, iState%lonr                               
        do j = 1, iState%lats_node_r      
          centerX(i,j) = 57.29578* iState%xlon(i,j)
          centerY(i,j) = 57.29578* iState%xlat(i,j)
        enddo
      enddo

!  Make sure we've got it right
!  ----------------------------
      call ESMF_GridValidate(mGrid, rc=rc1 )
      rc = 0

      return
      end subroutine gfs_physics_grid_create_Gauss3D

!..................................................................................

!  This routine came from MAPL...
      subroutine GridGetInterior_(GRID,I1,IN,J1,JN)
      type (ESMF_Grid), intent(IN) :: grid
      integer, intent(OUT)         :: I1, IN, J1, JN

! local vars
      integer                               :: status
      character(len=ESMF_MAXSTR)            :: IAm='MAPL_GridGetInterior'

      type (ESMF_DistGrid)                  :: distGrid
      type(ESMF_DELayout)                   :: LAYOUT
      integer,               allocatable    :: AL(:,:)
      integer,               allocatable    :: AU(:,:)
      integer                               :: nDEs
      integer                               :: deId
      integer                               :: gridRank
      integer                               :: deList(1)

      integer                               :: ierr, rc1


      CALL ESMF_LogWrite("GridGetInterior_", ESMF_LOGMSG_INFO, rc=rc1)

      call ESMF_GridGet    (GRID, dimCount=gridRank, distGrid=distGrid, rc=rc1) 
      call ESMF_DistGridGet(distGRID, delayout=layout, rc=STATUS)
!     call ESMF_DELayoutGet(layout, deCount =nDEs, localDeList=deList, rc=rc1)
      call ESMF_DELayoutGet(layout, deCount =nDEs, localDeToDeMap=deList, rc=rc1)

      deId = deList(1)

      allocate (AL(gridRank,0:nDEs-1), stat = ierr )
      allocate (AU(gridRank,0:nDEs-1), stat = ierr )


      call ESMF_DistGridGet(distgrid, minIndexPDe=AL, maxIndexPDe=AU, rc=rc1 )

      I1 = AL(1, deId)
      IN = AU(1, deId)
      J1 = AL(2, deId)
      JN = AU(2, deId)
      deallocate(AU, AL)

      end subroutine GridGetInterior_

      end module gfs_physics_grid_create_mod
