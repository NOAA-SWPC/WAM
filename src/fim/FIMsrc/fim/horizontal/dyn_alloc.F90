!*********************************************************************
	module module_dyn_alloc
!	This module allocates variables used in dynamics
!  	Middlecoff         October       2008
!*********************************************************************
contains
subroutine dyn_alloc
use module_control,only: nvl,nvlp1,npp,nip,nabl,ntra,ntrb,nd
use module_wrf_control,only: nbands
use module_constants
use module_variables
use module_sfc_variables
use module_chem_variables,only: sscal,ext_cof,asymp,extlw_cof
use infnan, only: inf, negint

implicit none

! Allocate variables from module_constants

!..................................................................
!	Sec. 1  Math and Physics Constants
!..................................................................

!..................................................................
!	Sec. 2.  Grid Descriptive Variables
!..................................................................

allocate(dpsig (nvl))           ! list of minimum layer thickness (Pa)
allocate(thetac(nvl))           ! target theta for hybgen

! velocity transform constants for projection from cell edges
allocate(cs(4,npp,nip),sn(4,npp,nip))
#ifndef LAHEY
dpsig(:) = inf
thetac(:) = inf
cs(:,:,:) = inf
sn(:,:,:) = inf
#endif

! Variables to describe the icos grid in xy (local stereographic)
allocate(sidevec_c(nd,npp,nip)) ! side vectors projected from center
allocate(sidevec_e(nd,npp,nip)) ! side vectors projected from edge
allocate(sideln   (   npp,nip)) ! the length of side vectors (m)
allocate(rsideln  (   npp,nip)) ! reciprocal of "sideln" (m**-1)
allocate(rprox_ln (   npp,nip)) ! reciprocal of distance cell cent to prox pts
allocate(area     (       nip)) ! the area of cell polygon (m**2)
allocate(rarea    (       nip)) ! reciprocal of the "area"

#ifndef LAHEY
sidevec_c(:,:,:) = inf
sidevec_e(:,:,:) = inf
sideln(:,:) = inf
rsideln(:,:) = inf
rprox_ln(:,:) = inf
area(:) = inf
rarea(:) = inf
#endif

!.................................................................
!	Sec. 3.  Neighbor Lookup Tables etc.
!.................................................................

allocate(prox       (npp,nip))  ! Holds index of proximity points
allocate(proxs      (npp,nip))  ! Holds index of proximity sides
allocate(nprox      (    nip))  ! Holds number of proximity points
allocate(inv_perm   (    nip))  ! inverse permutation of grid indices

#ifndef LAHEY
prox(:,:) = negint
proxs(:,:) = negint
nprox(:) = negint
inv_perm(:) = negint
#endif

! nedge holds the number of edges valid at each grid cell on this task.  
! For a serial ! case, nedge == nprox.  
! For a parallel case, nedge == nprox on "interior" cells 
! and, nedge < nprox on "halo" cells.  
allocate(nedge      (    nip))

! permedge stores a look-up table for edge indexes.  
! For a serial case, permedge does nothing:  
!   permedge(:,ipn) = 1, 2, 3, ... nprox(ipn)
! For a parallel case, permedge does nothing on "interior" cells.  
! For a parallel case, permedge skips "missing" edges on "halo" cells.  
allocate(permedge   (npp,nip))
allocate(actual     (    nip))  ! actual ipn of halo points (others too)

#ifndef LAHEY
nedge(:) = negint
permedge(:,:) = negint
actual(:) = negint
#endif

!.....................................................................
!	Sec. 4.  Geo Variables:
!.....................................................................

allocate(corio(nip))			! Coriolis acceleration 
allocate(lat(nip),lon(nip))		! lat and lon in radians
allocate(deg_lat(nip),deg_lon(nip))	! lat and lon in degrees

#ifndef LAHEY
corio(:) = inf
lat(:) = inf
lon(:) = inf
deg_lat(:) = inf
deg_lon(:) = inf
#endif

! Allocate variables from module_variables

!.....................................................................
!	Sec. 1.  3D Primary Variables
!.....................................................................
!  State variables at center point of cell for 3D grid:

!  Layer variables are defined in the middle of the layer

allocate( us3d(nvl,nip))      ! zonal wind (m/s)
allocate( vs3d(nvl,nip))      ! meridional wind (m/s)
allocate( ws3d(nvl,nip))      ! vertical wind (Pa/s)
allocate( dp3d(nvl,nip))      ! press.diff. between coord levels (Pa)
allocate( dpinit(nvl,nip))    ! lyr thknss for class B tracer transport
allocate( mp3d(nvl,nip))      ! Montgomery Potential (m^2/s^2)
allocate( tk3d(nvl,nip))      ! temperature, kelvin
allocate( vor (nvl,nip))      ! absolute vorticity (s^-1)
allocate( tr3d(nvl,nip,ntra+ntrb)) ! 1=pot.temp, 2=water vapor, 3=cloud water, 4=ozone
allocate( trdp(nvl,nip,ntra+ntrb)) ! (tracer times dp3d ) for tracer transport eq.
allocate( rh3d(nvl,nip))      ! relative humidity from 0 to 1
allocate( qs3d(nvl,nip))      ! saturation specific humidity
allocate( pw2d(nip))          ! precipitable water

#ifndef LAHEY
us3d(:,:) = inf
vs3d(:,:) = inf
ws3d(:,:) = inf
dp3d(:,:) = inf
dpinit(:,:) = inf
mp3d(:,:) = inf
tk3d(:,:) = inf
vor(:,:) = inf
tr3d(:,:,:) = inf
trdp(:,:,:) = inf
rh3d(:,:) = inf
qs3d(:,:) = inf
pw2d(:) = inf
#endif

!  Level variables defined at layer interfaces
allocate( pr3d(nvlp1,nip))    ! pressure (pascal)
allocate( ex3d(nvlp1,nip))    ! exner function
allocate( ph3d(nvlp1,nip))    ! geopotential (=gz), m^2/s^2
allocate( sdot(nvlp1,nip))    ! mass flux across interfaces, sdot*(dp/ds)

#ifndef LAHEY
pr3d(:,:) = inf
ex3d(:,:) = inf
ph3d(:,:) = inf
sdot(:,:) = inf
#endif

!..................................................................
!	Sec. 2. Edge Variables
!..................................................................
!  Variables carried at the midpoints of the 6(5) sides of each cell
allocate( u_edg   (nvl,npp,nip))    ! u on edge
allocate( v_edg   (nvl,npp,nip))    ! v on edge
allocate( dp_edg  (nvl,npp,nip))    ! dp on edge
allocate( trc_edg (nvl,npp,nip,ntra+ntrb))! tracers on edge
allocate( lp_edg  (nvl,npp,nip))    ! mid-layer pressure on edge
allocate( bnll_edg(nvl,npp,nip))    ! bernoulli fct (montg + kin energy) on edge
allocate( massfx (nvl,npp,nip,3))   ! mass fluxes on edge
allocate( cumufx  (nvl,npp,nip))    ! time-integrated mass flx on edge

#ifndef LAHEY
u_edg(:,:,:) = inf
v_edg(:,:,:) = inf
dp_edg(:,:,:) = inf
trc_edg(:,:,:,:) = inf
lp_edg(:,:,:) = inf
bnll_edg(:,:,:) = inf
massfx(:,:,:,:) = inf
cumufx(:,:,:) = inf
#endif

!.....................................................................
! 	Sec. 3. Forcing (tendency) Variables
!.....................................................................
allocate( u_tdcy  (nvl,nip,nabl))     ! forcing of u
allocate( v_tdcy  (nvl,nip,nabl))     ! forcing of v
allocate( dp_tdcy (nvl,nip,nabl))     ! forcing of dp
allocate( dpl_tdcy(nvl,nip,nabl))     ! forcing dp, low order
allocate( trc_tdcy(nvl,nip,nabl,ntra+ntrb))! forcing of tracers
allocate( trl_tdcy(nvl,nip,nabl,ntra+ntrb))! forcing of tracers, low order
allocate( u_tdcy_phy  (nvl,nip))     ! physics forcing of u
allocate( v_tdcy_phy  (nvl,nip))     ! physics forcing of v
allocate( trc_tdcy_phy(nvl,nip,ntra+ntrb))! physics forcing of tracers

#ifndef LAHEY
u_tdcy(:,:,:) = inf
v_tdcy(:,:,:) = inf
dp_tdcy(:,:,:) = inf
dpl_tdcy(:,:,:) = inf
trc_tdcy(:,:,:,:) = inf
trl_tdcy(:,:,:,:) = inf
u_tdcy_phy(:,:) = inf
v_tdcy_phy(:,:) = inf
trc_tdcy_phy(:,:,:) = inf
#endif

!....................................................................
!       Sec. 4. Misc. arrays
!....................................................................
allocate(  work2d(nip  ))
allocate( iwork2d(nip  ))
allocate( psrf   (nip  ))   ! surface pressure
allocate( ptdcy  (nip,2))   ! sfc.pres.tdcy at 2 consec. time levels
allocate( worka  (nvl,nip)) ! 3d work array
allocate( workb  (nvl,nip)) ! 3d work array
allocate( diaga(nvl,nip))   ! diagnostic, for output, fill with anything
allocate( diagb(nvl,nip))   ! diagnostic, for output, fill with anything

!TODO Initialize diag* arrays to inf instead of zero. Init to zero currently required or
!TODO da and db fields will be wrong on history files.
diaga(:,:) = 0.
diagb(:,:) = 0.

#ifndef LAHEY
work2d(:) = inf
iwork2d(:) = negint
psrf(:) = inf
#endif

! Allocate variants of physics variables from module module_sfc_variables
! These are single-precision copies of physics variables passed from physics 
! to dynamics via the coupler.  
!JR These 5 things were moved from output.F90 so they can be written to the restart file.
allocate (rn2d0(nip))
allocate (rc2d0(nip))
allocate (rg2d0(nip))
allocate (flxswavg2d(nip))
allocate (flxlwavg2d(nip))

allocate(rn2d(nip))       ! accumulated total precipitation/rainfall
allocate(rc2d(nip))       ! accumulated convective precipitation/rainfall 
allocate(ts2d(nip))       ! skin temperature
allocate(sst_prev(nip))   ! skin temperature
allocate(sst_next(nip))   ! skin temperature
allocate(us2d(nip))       ! friction velocity/equivalent momentum flux
allocate(hf2d(nip))       ! sensible heat flux
allocate(qf2d(nip))       ! water vapor/equivalent latent heat flux
allocate(sheleg2d(nip))
allocate(tg32d(nip))
allocate(zorl2d(nip))
allocate(vfrac2d(nip))
allocate(vtype2d(nip))
allocate(stype2d(nip))
allocate(cv2d(nip))
allocate(cvb2d(nip))
allocate(cvt2d(nip))
allocate(alvsf2d(nip))
allocate(alvwf2d(nip))  
allocate(alnsf2d(nip))
allocate(alnwf2d(nip))
allocate(f10m2d(nip))
allocate(facsf2d(nip))
allocate(facwf2d(nip))
allocate(uustar2d(nip))
allocate(ffmm2d(nip))
allocate(ffhh2d(nip))
allocate(srflag2d(nip))
allocate(snwdph2d(nip))
allocate(shdmin2d(nip))
allocate(shdmax2d(nip))
allocate(slope2d(nip))
allocate(snoalb2d(nip))
allocate(canopy2d(nip))
allocate(hice2d(nip))
allocate(fice2d(nip))
allocate(fice2d_prev(nip))
allocate(fice2d_next(nip))
allocate(sw2d(nip))       ! downward short-wave radiation flux
allocate(lw2d(nip))       ! downward long-wave radiation flux
allocate(slc2d(nip))
allocate(t2m2d(nip))
allocate(q2m2d(nip))
allocate(slmsk2d(nip))
allocate(tprcp2d(nip))  ! precip rate (1000*kg/m**2)
allocate(st3d(4,nip))   ! soil temperature
allocate(sm3d(4,nip))   ! total soil moisture
allocate(slc3d(4,nip))  ! soil liquid content
allocate(hprm2d(14,nip))
allocate(flxlwtoa2d(nip))

#ifndef LAHEY
rn2d0(:) = inf
rc2d0(:) = inf
rg2d0(:) = inf
flxswavg2d(:) = inf
flxlwavg2d(:) = inf

rn2d(:) = inf
rc2d(:) = inf
ts2d(:) = inf
sst_prev(:) = inf
sst_next(:) = inf
us2d(:) = inf
hf2d(:) = inf
qf2d(:) = inf
sheleg2d(:) = inf
tg32d(:) = inf
zorl2d(:) = inf
vfrac2d(:) = inf
vtype2d(:) = inf
stype2d(:) = inf
cv2d(:) = inf
cvb2d(:) = inf
cvt2d(:) = inf
alvsf2d(:) = inf
alvwf2d(:) = inf
alnsf2d(:) = inf
alnwf2d(:) = inf
f10m2d(:) = inf
facsf2d(:) = inf
facwf2d(:) = inf
uustar2d(:) = inf
ffmm2d(:) = inf
ffhh2d(:) = inf
srflag2d(:) = inf
snwdph2d(:) = inf
shdmin2d(:) = inf
shdmax2d(:) = inf
slope2d(:) = inf
snoalb2d(:) = inf
canopy2d(:) = inf
hice2d(:) = inf
fice2d(:) = inf
fice2d_prev(:) = inf
fice2d_next(:) = inf
sw2d(:) = inf
lw2d(:) = inf
slc2d(:) = inf
t2m2d(:) = inf
q2m2d(:) = inf
slmsk2d(:) = inf
tprcp2d(:) = inf
st3d(:,:) = inf
sm3d(:,:) = inf
slc3d(:,:) = inf
hprm2d(:,:) = inf
#endif

! chem arrays

! TODO We'd rather not allocate these arrays if chem is off. But Lahey does not
! TODO allow passing unallocated arrays as subroutine parameters. So, work needs
! TODO to be done to (maybe) have these arrays as optional parameters and wrap
! TODO accesses in "if (present())" conditionals -- or some other approach. For
! TODO now, we have to allocate these.

!if(aer_ra_feedback == 1 ) then
  allocate(sscal(nvl,nip,nbands))
  sscal=1.
  allocate(ext_cof(nvl,nip,nbands))
  ext_cof=0.
  allocate(asymp(nvl,nip,nbands))
  asymp=0.
  allocate(extlw_cof(nvl,nip,16))
  extlw_cof=0.
!endif

return
end subroutine dyn_alloc

end module module_dyn_alloc
