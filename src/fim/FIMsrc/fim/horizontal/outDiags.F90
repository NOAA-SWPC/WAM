module module_outDiags
contains
subroutine outDiags(its,nvl,nvlp1,nip,ntra,pr3d,ph3d,tr3d,rn2d,rc2d,sdot,dp3d,us3d,vs3d,rh3d,lw2d,sw2d,hf2d,qf2d)
use module_constants       ,only: deg_lat,deg_lon
USE module_sfc_variables ,only : slmsk2d
use module_outqv           ,only: outqv
use module_outqv_wsp       ,only: outqv_wsp
use module_outqv_mn        ,only: outqv_mn
use module_outqv_mn_lat    ,only: outqv_mn_lat
use module_outqv_mn_lat_abs,only: outqv_mn_lat_abs
use module_outqv_mn_lat_land,only: outqv_mn_lat_land
use module_out4d_mn        ,only: out4d_mn
use module_control         ,only: dt,ArchvTimeUnit

implicit none

! External variable declarations:
integer,intent(IN) :: its,nvl,nvlp1,nip,ntra
!SMS$DISTRIBUTE (dh,nip) BEGIN
real   ,intent(IN) :: us3d(nvl  ,nip),vs3d(nvl,nip),dp3d(nvl,nip)
real   ,intent(IN) :: pr3d(nvlp1,nip)
real   ,intent(IN) :: ph3d(nvlp1,nip)
real   ,intent(IN) :: tr3d(nvl    ,nip,ntra)
real   ,intent(IN) :: rh3d(nvl  ,nip)
real   ,intent(IN) :: sw2d(nip),lw2d(nip)
real   ,intent(IN) :: rn2d(nip),rc2d(nip)
real   ,intent(IN) :: sdot(nvlp1,nip)	! mass flux across interfaces, sdot*(dp/ds)
real   ,intent(IN) :: hf2d(nip),qf2d(nip)
!SMS$DISTRIBUTE END

!factors for changing units for outqv range, max/min
!water vapor - g/g to g/kg
real :: fact_qv = 1.e3
real :: fact_qc = 1.e5
!pressure - Pa to hPa
real :: fact_pres = 1.e-2
real :: fact    = 1.
real :: time

integer :: LB
integer, external :: its2time

time=its2time(its)
write (6,'(a,f9.2,1x,2a,i8)') ' Pressure at ',time,ArchvTimeUnit,', time step=',its
call outqv_mn (pr3d,nvlp1,nip,its,fact_pres)
write (6,'(a,f9.2,1x,2a,i8)') ' Geopotential height at ',time,ArchvTimeUnit,', time step=',its
call outqv_mn (ph3d,nvlp1,nip,its,1.)
call outqv    (ph3d ,deg_lat,deg_lon,nip,nvlp1 ,fact)
write (6,'(a,f9.2,1x,2a,i8)') ' THETA at ',time,ArchvTimeUnit,', time step=',its
call out4d_mn (tr3d,nvl,nip,ntra,its,1.,1)
write (6,'(a,f9.2,1x,2a,i8)') ' Precip at ',time,ArchvTimeUnit,', time step=',its
call out4d_mn (rn2d,1,nip,1,its,1.,1)
write (6,'(a,f9.2,1x,2a,i8)') ' Precip at ',time,ArchvTimeUnit,', time step=',its
call outqv_mn_lat (rn2d,deg_lat,deg_lon,30.,1,nip,its,fact)
write (6,'(a,f9.2,1x,2a,i8)') ' Precip at ',time,ArchvTimeUnit,', time step=',its
call outqv    (rn2d ,deg_lat,deg_lon,nip,1 ,fact)
write (6,'(a,f9.2,1x,2a,i8)') ' Precip-conv at ',time,ArchvTimeUnit,', time step=',its
call out4d_mn (rc2d,1,nip,1,its,1.,1)
call outqv_mn_lat (rc2d,deg_lat,deg_lon,30.,1,nip,its,fact)
write (6,'(a,f9.2,1x,2a,i8)') ' Precip-conv at ',time,ArchvTimeUnit,', time step=',its
call outqv    (rc2d ,deg_lat,deg_lon,nip,1 ,fact)
write (6,'(a,f9.2,1x,2a,i8)') ' Vertical velocity at ',time,ArchvTimeUnit,', time step=',its
call out4d_mn (sdot,nvlp1,nip,1,its,1.,1)
write (6,'(a,f9.2,1x,2a,i8)') ' Vertical velocity at ',time,ArchvTimeUnit,', time step=',its
call outqv_mn_lat (sdot,deg_lat,deg_lon,30.,nvlp1,nip,its,fact)

write (6,'(a,f9.2,1x,2a,i8)') ' Vertical velocity at ',time,ArchvTimeUnit,', time step=',its
call outqv    ( sdot ,deg_lat,deg_lon,nip,nvlp1 ,fact)
write (6,'(a,f9.2,1x,2a,i8)') ' Abs Vertical velocity at ',time,ArchvTimeUnit,', time step=',its
call outqv_mn_lat_abs (sdot,deg_lat,deg_lon,30.,nvlp1,nip,its,fact)

LB = LBOUND(tr3d,2)
write (6,'(a,f8.1,a,f9.2,1x,2a,i8)') ' Water vapor - tr(2), fact=',fact_qv,' at ',time,ArchvTimeUnit,', time step=',its
call outqv    ( tr3d(1,LB,2) ,deg_lat,deg_lon,nip,nvl,fact_qv )
write (6,'(a,f8.1,a,f9.2,1x,2a,i8)') ' Cloud water - tr(3), fact=',fact_qc,' at ',time,ArchvTimeUnit,', time step=',its
call outqv    ( tr3d(1,LB,3) ,deg_lat,deg_lon,nip,nvl,fact_qc )
write (6,'(a,f9.2,1x,2a,i8)') ' DP3d at ',time,ArchvTimeUnit,', time step=',its
call outqv    ( dp3d ,deg_lat,deg_lon,nip,nvl ,fact_pres)
write (6,'(a,f9.2,1x,2a,i8)') ' pot temp at ',time,ArchvTimeUnit,', time step=',its
call outqv    ( tr3d(1,LB,1) ,deg_lat,deg_lon,nip,nvl ,fact)
write (6,'(a,f9.2,1x,2a,i8)') ' Water vapor - tr(2) at ',time,ArchvTimeUnit,', time step=',its
call outqv_mn (tr3d(1,LB,2),nvl,nip,its,fact_qv)
write (6,'(a,f8.1,a,f9.2,1x,2a,i8)') ' Cloud water - tr(3), fact=',fact_qc,' at ',time,ArchvTimeUnit,', time step=',its
call outqv_mn (tr3d(1,LB,3),nvl,nip,its,fact_qc)
write (6,'(a,f8.1,a,f9.2,1x,2a,i8)') ' Cloud water - tr(3), fact=',fact_qc,' at ',time,ArchvTimeUnit,', time step=',its
call outqv_mn_lat (tr3d(1,LB,3),deg_lat,deg_lon,30.,nvl,nip,its,fact_qc)
write (6,'(a,f9.2,1x,2a,i8)') ' Pot temp at ',time,ArchvTimeUnit,', time step=',its
call outqv_mn (tr3d(1,LB,1),nvl,nip,its,fact)
write (6,'(a,f9.2,1x,2a,i8)') ' Zonal wind at ',time,ArchvTimeUnit,', time step=',its
call outqv_mn (us3d,nvl,nip,its,1.)
write (6,'(a,f9.2,1x,2a,i8)') ' Zonal wind at ',time,ArchvTimeUnit,', time step=',its
call outqv    (us3d ,deg_lat,deg_lon,nip,nvl ,fact)
write (6,'(a,f9.2,1x,2a,i8)') ' Meridional wind at ',time,ArchvTimeUnit,', time step=',its
call outqv_mn (vs3d,nvl,nip,its,1.)
write (6,'(a,f9.2,1x,2a,i8)') ' Meridional wind at ',time,ArchvTimeUnit,', time step=',its
call outqv    (vs3d ,deg_lat,deg_lon,nip,nvl ,fact)

!  Call for wind speed
write (6,'(a,f9.2,1x,2a,i8)') ' Wind speed at ',time,ArchvTimeUnit,', time step=',its
call outqv_wsp    (us3d, vs3d ,deg_lat,deg_lon,nip,nvl ,fact)

write (6,'(a,f9.2,1x,2a,i8)') ' Relative humidity at ',time,ArchvTimeUnit,', time step=',its
call outqv_mn (rh3d,nvl,nip,its,1.)
write (6,'(a,f9.2,1x,2a,i8)') ' Longwave at ',time,ArchvTimeUnit,', time step=',its
call outqv_mn (lw2d,1,nip,its,1.)
write (6,'(a,f9.2,1x,2a,i8)') ' Shortwave at ',time,ArchvTimeUnit,', time step=',its
call outqv_mn (sw2d,1,nip,its,1.)
write (6,'(a,f9.2,1x,2a,i8)') ' Sensible heat flux at ',time,ArchvTimeUnit,', time step=',its
call outqv_mn (hf2d,1,nip,its,1.)
call outqv    (hf2d ,deg_lat,deg_lon,nip,1,1.)
call outqv_mn_lat_land (hf2d  ,deg_lat,deg_lon,30.,slmsk2d,1,1,nip,its,1.)
write (6,'(a,f9.2,1x,2a,i8)') ' Latent heat flux at ',time,ArchvTimeUnit,', time step=',its
call outqv_mn (qf2d,1,nip,its,1.)
call outqv    (qf2d ,deg_lat,deg_lon,nip,1,1.)
call outqv_mn_lat_land (qf2d  ,deg_lat,deg_lon,30.,slmsk2d,1,1,nip,its,1.)

end subroutine outDiags
end module module_outDiags
