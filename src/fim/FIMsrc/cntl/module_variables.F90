module module_variables
!*********************************************************************
!	This module specifies the variables used in the fim
!  	A. E. MacDonald		          October 11  2004
!  	J. Lee        		          September   2004
!  	R. Bleck     Cleanup              April       2008
!       Middlecoff   Dynamic allocation   October     2008
!       Henderson    Move PHY vars. to    June        2009
!                    module_sfc_variables.  
!*********************************************************************

implicit none
save
!.....................................................................
!	Sec. 1.  3D Primary Variables
!.....................................................................
!  State variables at center point of cell for 3D grid:

!  Layer variables are defined in the middle of the layer
!SMS$DISTRIBUTE(dh,2) BEGIN
real,pointer     :: us3d  (:,:)	  ! zonal wind (m/s)
real,pointer     :: vs3d  (:,:)	  ! meridional wind (m/s)
real,pointer     :: ws3d  (:,:)	  ! vertical wind (Pa/s)
real,allocatable :: dp3d  (:,:)	  ! del p between coord levels (pascals)
real,allocatable :: dpinit(:,:)	  ! lyr thknss for class B tracer transport
real,allocatable :: mp3d  (:,:)	  ! Montgomery Potential (m^2/s^2)
real,allocatable :: tk3d  (:,:)	  ! temperature, kelvin
real,allocatable :: vor   (:,:)	  ! absolute vorticity (s^-1)
real,pointer     :: tr3d  (:,:,:) ! 1=pot.temp,2=water vapor,3=cloud water,4=ozone
real,allocatable :: trdp  (:,:,:) ! (tracer x thknss) for tracer transport eq.
real,allocatable :: rh3d  (:,:)   ! relative humidity from 0 to 1
real,allocatable :: qs3d  (:,:)   ! saturation specific humidity
!SMS$DISTRIBUTE END
!SMS$DISTRIBUTE(dh,1) BEGIN
real,allocatable :: pw2d(:)       ! precipitable water
!SMS$DISTRIBUTE END

!SMS$DISTRIBUTE(dh,2) BEGIN
!  Level variables defined at layer interfaces
real,pointer     :: pr3d(:,:)     ! pressure (pascal)
real,allocatable :: ex3d(:,:)     ! exner function
real,allocatable :: ph3d(:,:)     ! geopotential (=gz), m^2/s^2
real,allocatable :: sdot(:,:)     ! mass flux across interfaces, sdot*(dp/ds)

!.....................................................................
! 	Sec. 3. Forcing (tendency) Variables
!.....................................................................
integer          :: curr_write_time = 0     ! For time avg'd output vars: restart will overwrite
integer          :: nf,of,vof		    ! Adams Bashforth time slots
real             :: adbash1,adbash2,adbash3 ! Adams Bashforth weights
real,allocatable :: u_tdcy  (:,:,:)	    ! forcing of u
real,allocatable :: v_tdcy  (:,:,:)	    ! forcing of v
real,allocatable :: dp_tdcy (:,:,:)	    ! forcing of dp
real,allocatable :: dpl_tdcy(:,:,:)	    ! forcing dp, low order
real,allocatable :: trc_tdcy(:,:,:,:)       ! forcing of tracers
real,allocatable :: trl_tdcy(:,:,:,:)       ! forcing of tracers, low order
real,allocatable :: massflx(:,:,:) 
real,allocatable :: u_tdcy_phy  (:,:)         ! physics forcing of u
real,allocatable :: v_tdcy_phy  (:,:)         ! physics forcing of v
real,allocatable :: trc_tdcy_phy(:,:,:)       ! physics forcing of tracers
!SMS$DISTRIBUTE END

!SMS$DISTRIBUTE(dh,3) BEGIN
!..................................................................
!	Sec. 4. Edge Variables
!..................................................................
!  Variables carried at the midpoints of the 6(5) sides of each cell
real,allocatable :: u_edg   (:,:,:)   ! u on edge
real,allocatable :: v_edg   (:,:,:)   ! v on edge
real,allocatable :: dp_edg  (:,:,:)   ! dp on edge
real,allocatable :: trc_edg (:,:,:,:) ! tracers on edge
real,allocatable :: lp_edg  (:,:,:)   ! mid-layer pressure on edge
real,allocatable :: bnll_edg(:,:,:)   ! bernoulli fct (montg + kin energy) on edge
real,allocatable :: massfx (:,:,:,:)  ! mass fluxes on edge at 3 time levels
real,allocatable :: cumufx (:,:,:)    ! time-integrated mass fluxes on edges
!SMS$DISTRIBUTE END

!....................................................................
!       Sec. 5.  Misc. arrays
!....................................................................
!SMS$DISTRIBUTE(dh,1) BEGIN
real   ,allocatable :: work2d (:  )
integer,allocatable :: iwork2d(:  )
real   ,allocatable :: psrf   (:  ) ! surface pressure
real   ,allocatable :: ptdcy  (:,:) ! sfc.pres.tdcy at 2 consec.time levels
!SMS$DISTRIBUTE END

!SMS$DISTRIBUTE(dh,2) BEGIN
real   ,allocatable :: worka  (:,:) ! 3d work array
real   ,allocatable :: workb  (:,:) ! 3d work array
!SMS$DISTRIBUTE END

!SMS$DISTRIBUTE(dh,2) BEGIN
real   ,allocatable :: diaga  (:,:) ! 3d diagnosic array
real   ,allocatable :: diagb  (:,:) ! 3d diagnosic array
!SMS$DISTRIBUTE END
end module module_variables
