MODULE module_constants
!********************************************************************
!       This module specifies model constants for fim
!       A. E. MacDonald         October 11, 2004
!       J. LEE                  September,  2005
!       J. LEE: values for physical constants are taken from RUC (METCON)
!********************************************************************
implicit none

!..................................................................
!	Sec. 1  Math and Physics Constants
!..................................................................

real, parameter   :: pi       = 3.14159265
real, parameter   :: degrad   = pi/180.
real, parameter   :: raddeg   = 180./pi
real, parameter   :: ae	= 6371220.         ! earth radius (m)
real, parameter   :: omegx2	= 1.4584e-4! 2 x earth rotation rate (s^-1)
real, parameter   :: grvity   = 9.80665    ! acceleration from earth gravity (m/(s^2))
real, parameter   :: p1000    = 100000.    ! p at 1000mb (pascals)
real, parameter   :: cp	= 1004.6855        ! specific heat at const pres
real, parameter   :: rd	= 287.0586         ! spec gas constant for dry air
real, parameter   :: rv	= 461.50           ! gas constant for H2O
real, parameter   :: qvmin    = 1.e-10
real, parameter   :: qwmin    = 1.e-10

!..................................................................
!	Sec. 2.  Grid Descriptive Variables
!..................................................................
real, allocatable :: dpsig (:)          ! list of minimum layer thknss (Pa)
real, allocatable :: thetac(:)          ! target theta for hybgen
real, allocatable :: sigak(:)           ! sigma coordinate formula (used
real, allocatable :: sigbk(:)           !     ...only if pure_sig = true)

!SMS$DISTRIBUTE(dh,3) BEGIN
! velocity transform constants for projection from cell edges
real,allocatable :: cs(:,:,:),sn(:,:,:)

! Variables to describe the icos grid in xy (local stereographic)
real,allocatable ::  sidevec_c(:,:,:) ! side vectors projected from center
real,allocatable ::  sidevec_e(:,:,:) ! side vectors projected from edge
!SMS$DISTRIBUTE END

!SMS$DISTRIBUTE(dh,2) BEGIN
real,allocatable ::  sideln   (  :,:) ! the length of side vectors (m)
real,allocatable ::  rsideln  (  :,:) ! reciprocal of "sideln" (m**-1)
real,allocatable ::  rprox_ln (  :,:) ! reciprocal of distance cell cent to prox pts
!SMS$DISTRIBUTE END
!SMS$DISTRIBUTE(dh,1) BEGIN
real,allocatable ::  area     (    :) ! the area of cell polygon (m**2)
real,allocatable ::  rarea    (    :) ! reciprocal of the "area"
!SMS$DISTRIBUTE END

!.................................................................
!	Sec. 3.  Geographic Indices
!.................................................................

!SMS$DISTRIBUTE(dh,2) BEGIN
integer,allocatable :: prox       (:,:)  ! Holds index of proximity points
integer,allocatable :: proxs      (:,:)  ! Holds index of proximity sides
! permedge stores a look-up table for edge indexes.  
! For a serial case, permedge does nothing:  
!   permedge(:,ipn) = 1, 2, 3, ... nprox(ipn)
! For a parallel case, permedge does nothing on "interior" cells.  
! For a parallel case, permedge skips "missing" edges on "halo" cells.  
integer,allocatable :: permedge   (:,:)
!SMS$DISTRIBUTE END
!SMS$DISTRIBUTE(dh,1) BEGIN
integer,allocatable :: nprox      (  :)  ! Holds number of proximity points
integer,allocatable,target :: inv_perm   (  :)  ! inverse permutation of the grid
! nedge holds the number of edges valid at each grid cell on this task.  
! For a serial ! case, nedge == nprox.  
! For a parallel case, nedge == nprox on "interior" cells 
! and, nedge < nprox on "halo" cells.  
integer,allocatable :: nedge  (:)
integer,allocatable :: actual (:) ! actual ipn of halo points (others too)

!.....................................................................
!	Sec. 4.  Geo Variables:
!.....................................................................

real,allocatable ::  corio  (:)			! Coriolis acceleration 
real,allocatable ::      lat(:),    lon(:)	! lat and lon in radians
real,allocatable ::  deg_lat(:),deg_lon(:)	! lat and lon in degrees
character(20) :: StartDate = '2000:07:26::15:25:28'
!SMS$DISTRIBUTE END

!.....................................................................
!      Other variables
!.....................................................................
real,parameter    :: spval_p   = 999999.

END MODULE module_constants
