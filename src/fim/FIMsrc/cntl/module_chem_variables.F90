!*********************************************************************
module module_chem_variables
!	This module specifies chem variables.  
!*********************************************************************

save

!SMS$DISTRIBUTE(dh,2) BEGIN
real,allocatable :: tr1_tavg (:,:)! tracer time average
real,allocatable :: pm25  (:,:)   ! pm2.5
real,allocatable :: p10   (:,:)   ! pm10
!real,allocatable :: exch  (:,:)   ! exchange coeffs for chemistry transport
real,allocatable :: oh_backgd(:,:) ! OH background for GOCART
real,allocatable :: h2o2_backgd(:,:) ! H2O2 background for GOCART
real,allocatable :: no3_backgd(:,:) ! NO3 background for GOCART
real,allocatable :: sscal (:,:,:)   ! aerosol single scattering albedo
real,allocatable :: ext_cof (:,:,:)   ! aerosol extinction coefficients
real,allocatable :: extlw_cof (:,:,:)   ! aerosol extinction coefficients
real,allocatable :: asymp (:,:,:)   ! aerosol asymetry parameter

!SMS$DISTRIBUTE END
!SMS$DISTRIBUTE(dh,1) BEGIN
real,allocatable :: rcav(:)       ! accumulated convective precipitation since last chem call
real,allocatable :: rnav(:)       ! accumulated ls precipitation since last chem call
!real,allocatable :: pb2d(:)       ! Boundary layer height
real,allocatable :: ero1(:)       !  dust erosion factor
real,allocatable :: ero2(:)       !  dust erosion factor
real,allocatable :: ero3(:)       !  dust erosion factor
real,allocatable :: clayfrac(:)   !  clay fraction (AFWA dust scheme)
real,allocatable :: sandfrac(:)   !  sand fraction (AFWA dust scheme)
real,allocatable :: dm0(:)        !  dms reference emissions
real,allocatable :: emiss_ab(:,:) !  emissions for all available species
real,allocatable :: emiss_abu(:,:)!  emissions for all available species
real,allocatable :: trfall(:,:)  !  emissions for all available species
real,allocatable :: emiss_ash_mass(:) !  emissions for 
real,allocatable :: emiss_ash_height(:) !  emissions for 
real,allocatable :: emiss_ash_dt(:) !  emissions for 
real,allocatable :: emiss_tr_mass(:) !  emissions for 
real,allocatable :: emiss_tr_height(:) !  emissions for 
real,allocatable :: emiss_tr_dt(:) !  emissions for 
real,allocatable :: emiss_co2(:)   !  emissions for co2, 8 time slices
real,allocatable :: emiss_sf6(:)   !  emissions for sf6, 8 time slices ?
real,allocatable :: plumestuff(:,:) !  fire info
real,allocatable :: aod2d(:)   !  aerosol optical depth
real,allocatable :: dustfall(:)   !  dust fall (g/m2)
real,allocatable :: ashfall(: )   !  volcanic ash fall (g/m2)
real,allocatable :: emiss_oc(:)   !  emissions for organic carbon
real,allocatable :: emiss_bc(:)   !  emissions for black carbon
real,allocatable :: emiss_sulf(:) !  emissions for sulfate
!SMS$DISTRIBUTE END

end module module_chem_variables
