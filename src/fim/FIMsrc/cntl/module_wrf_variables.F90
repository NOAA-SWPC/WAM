!*********************************************************************
module module_wrf_variables
!	This module specifies WRF physics variables.  
!*********************************************************************

implicit none

save

!SMS$DISTRIBUTE(dh,2) BEGIN
real,allocatable :: phys3dwrf(:,:,:) ! Physics diagnostic variable
real,allocatable :: exch  (:,:)   ! exchange coeffs
!SMS$DISTRIBUTE END
!SMS$DISTRIBUTE(dh,1) BEGIN
real,allocatable :: phys2dwrf(:,:)   ! Physics diagnostic variable
real,allocatable :: pb2d(:)       ! Boundary layer height
!SMS$DISTRIBUTE END

end module module_wrf_variables
