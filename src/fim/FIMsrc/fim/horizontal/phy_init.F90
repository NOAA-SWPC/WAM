module module_fim_phy_init

private
public :: phy_init, sst_init, sstunit

integer, save :: sstunit = -1

contains

subroutine phy_init
!*********************************************************************
!       Loads the initial variables and constants for the physics 
!       component.  
!       Alexander E. MacDonald  11/27/04
!       J. Lee                  September, 2005
!*********************************************************************

use module_control  ,only: nvl,nip,dt,numphr,                                &
                           readrestart,alt_topo,                             &
                           glvl,curve,NumCacheBlocksPerPE,                   &
                           PhysicsInterval,RadiationInterval,                &
                           CallPhysics,CallRadiation,GravityWaveDrag,        &
                           ras,num_p3d,gfsltln_file,mtnvar_file,aerosol_file,&
                           co2_2008_file,co2_glb_file,SSTInterval,UpdateSST
use funcphys   ! GFS physics
use module_sfc_variables
!SMS$IGNORE BEGIN
USE gfs_physics_internal_state_mod, only: gfs_physics_internal_state, gis_phy
USE gfs_physics_sfc_flx_set_mod, only: sfcvar_aldata, flxvar_aldata, flx_init
use units, only: getunit, returnunit
!SMS$IGNORE END

implicit none

! Local variables

integer       :: mype,idx,ipn,ierr,mylb,myub
integer       :: unitno
real*8        :: t0,t1=0.0d0
logical       :: first_loop
!real :: cv2d(nip),cvt2d(nip),cvb2d(nip),slmsk2d(nip),ts2d(nip)
!real :: st3d(4,nip),sheleg2d(nip),snoalb2d(nip)
!real :: hprm2d(14,nip),fice2d(nip),tprcp2d(nip)
!real :: slc3d(4,nip),sm3d(4,nip),snwdph2d(nip),slope2d(nip),shdmin2d(nip)
!real :: shdmax2d(nip),tg32d(nip),canopy2d(nip)
!real :: alvsf2d(nip),alnsf2d(nip),alvwf2d(nip),alnwf2d(nip)
!real :: facsf2d(nip),facwf2d(nip),t2m2d(nip),q2m2d(nip),uustar2d(nip)
!real :: work2d(nip),slc2d(nip)
!real :: ffmm2d(nip),ffhh2d(nip),f10m2d(nip)
!SMS$DISTRIBUTE (dh,nip) BEGIN
real :: work2d(nip)
!SMS$DISTRIBUTE END

namelist /PREPnamelist/ curve,NumCacheBlocksPerPE,alt_topo,gfsltln_file,mtnvar_file  &
                       ,aerosol_file,co2_2008_file,co2_glb_file

namelist /PHYSICSnamelist/ PhysicsInterval,RadiationInterval,SSTInterval,GravityWaveDrag,ras,num_p3d

! TODO:  Create new decomp "dhp" and use it for all physics declarations 
! TODO:  and loops.  Split modules that contain variables used by both too!  
!!SMS$CREATE_DECOMP(dh,<nip>,<HALO_SIZE>)

!!allocate( hice2d(nip) )
!!allocate( srflag2d(nip) )
!!allocate( zorl2d(nip) )
!!allocate( vfrac2d(nip) )
!!allocate( vtype2d(nip) )
!!allocate( stype2d(nip) )

!SMS$insert call nnt_me(mype)
call StartTimer(t0)

! TBH:  BEGIN DUPLICATION with dyn_init().  REFACTOR TO REMOVE DUPLICATION
!TODO:  call control() here before reading PHYSICSnamelist...  
! Note:  REWIND required by IBM!  
! TODO:  Using open-read-close in place of REWIND until SMS is updated
unitno = getunit ()
if (unitno < 0) then
  print*,'phy_init: getunit failed: stopping'
  stop
end if

open  (unitno, file='./FIMnamelist', form='formatted', action='read', err=70)
read  (unitno, PREPnamelist, err=90)
close (unitno)
open  (unitno, file='./FIMnamelist', form='formatted', action='read', err=70)
read  (unitno, PHYSICSnamelist, err=90)
close (unitno)

call returnunit (unitno)
! TBH:  END DUPLICATION with dyn_init().  

call IncrementTimer(t0,t1)

CallPhysics   = max(1,numphr*PhysicsInterval/3600)
CallRadiation = max(1,(max(1,numphr*RadiationInterval/3600)/CallPhysics) &
                *CallPhysics)
print "(' Calculate gravity wave drag        ',L10)",GravityWaveDrag
print "(' Call physics every'  ,I27,' timesteps (',I0,' seconds)')", &
          CallPhysics,nint(CallPhysics*dt)
print "(' Call radiation every',I25,' timesteps (',I0,' seconds)')", &
          CallRadiation,nint(CallRadiation*dt)
print "(' ras =',L10)",ras
print "(' num_p3d =',I10)",num_p3d

write(6,PHYSICSnamelist)

! On a restart run, these variables will be read in from the restart file so this part can be skipped.

if (.not.readrestart) then
!SMS$PARALLEL(dh, ipn) BEGIN
  call StartTimer(t0)
!SMS$SERIAL BEGIN
  unitno = getunit ()
  if (unitno < 0) then
    print*,'phy_init: getunit failed: stopping'
    stop
  end if

  open                (unitno, file="gfsfc.dat", form="unformatted", action='read', err=70)
  call TestGlvlHeader (unitno,     "gfsfc.dat",'phy_init',glvl )
  call TestCurveHeader(unitno,     "gfsfc.dat",'phy_init',curve)
!SMS$SERIAL END
  do idx = 1,SIZE(st3d,1)
!SMS$SERIAL BEGIN
    read(unitno, err=90) work2d
!SMS$SERIAL END
    do ipn=1,nip
      st3d(idx,ipn) = work2d(ipn)
    enddo
  enddo
! soil moisture
  do idx = 1,SIZE(sm3d,1)
!SMS$SERIAL BEGIN
    read(unitno, err=90) work2d
!SMS$SERIAL END
    do ipn=1,nip
      sm3d(idx,ipn) = work2d(ipn)
    enddo
  enddo
  do idx = 1,SIZE(slc3d,1)
!SMS$SERIAL BEGIN
    read(unitno, err=90) work2d
!SMS$SERIAL END
    do ipn=1,nip
      slc3d(idx,ipn) = work2d(ipn)
    enddo
  enddo
!
! skin temperature
!
!SMS$SERIAL BEGIN
  read(unitno, err=90) ts2d
!SMS$SERIAL END
!
! snow water equivalent
!
!SMS$SERIAL BEGIN
  read(unitno, err=90) sheleg2d
!SMS$SERIAL END
!SMS$SERIAL BEGIN
  read(unitno, err=90) tg32d
!SMS$SERIAL END
!
! surface roughness
!
!SMS$SERIAL BEGIN
  read(unitno, err=90) zorl2d
!SMS$SERIAL END
!
! maybe conv cloud fraction?
!
!SMS$SERIAL BEGIN
  read(unitno, err=90) cv2d
!SMS$SERIAL END
!
! maybe conv cloud bottom pressure?
!
!SMS$SERIAL BEGIN
  read(unitno, err=90) cvb2d
!SMS$SERIAL END
!
! maybe conv cloud top pressure?
!
!SMS$SERIAL BEGIN
  read(unitno, err=90) cvt2d
!SMS$SERIAL END
!
! mean visible albedo with strong cosz dependence...???...
!
!SMS$SERIAL BEGIN
  read(unitno, err=90) alvsf2d
!SMS$SERIAL END
!
! mean vis albedo with weak cosz dependence...???...
!
!SMS$SERIAL BEGIN
  read(unitno, err=90) alvwf2d
!SMS$SERIAL END
!
! mean nir albedo with strong cosz dependence...???...
!
!SMS$SERIAL BEGIN
  read(unitno, err=90) alnsf2d
!SMS$SERIAL END
!
! mean nir albedo with weak cosz dependence...???...
!
!SMS$SERIAL BEGIN
  read(unitno, err=90) alnwf2d
!SMS$SERIAL END
!
! land/sea/ice mask (0:SEA.1:LAND,2:ICE)
!SMS$SERIAL BEGIN
  read(unitno, err=90) slmsk2d
!SMS$SERIAL END
!
! assuming veg fraction
!
!SMS$SERIAL BEGIN
  read(unitno, err=90) vfrac2d
!SMS$SERIAL END
!
! canopy moisture content in mm
!
!SMS$SERIAL BEGIN
  read(unitno, err=90) canopy2d
!SMS$SERIAL END
!SMS$SERIAL BEGIN
  read(unitno, err=90) f10m2d
!SMS$SERIAL END
!SMS$SERIAL BEGIN
  read(unitno, err=90) t2m2d
!SMS$SERIAL END
!SMS$SERIAL BEGIN
  read(unitno, err=90) q2m2d
!SMS$SERIAL END
!
!   vegtype or landuse?
!
!SMS$SERIAL BEGIN
  read(unitno, err=90) vtype2d
!SMS$SERIAL END
!
! soilcategory
!
!SMS$SERIAL BEGIN
  read(unitno, err=90) stype2d
!SMS$SERIAL END
!
! fractional coverage with strong (facsf) and weak (facwf) cosz dependence
!
!SMS$SERIAL BEGIN
  read(unitno, err=90) facsf2d
!SMS$SERIAL END
!SMS$SERIAL BEGIN
  read(unitno, err=90) facwf2d
!SMS$SERIAL END
!
! ustar
!
!SMS$SERIAL BEGIN
  read(unitno, err=90) uustar2d
!SMS$SERIAL END
!
! looks like surface exchange coeffs for m and h
!SMS$SERIAL BEGIN
  read(unitno, err=90) ffmm2d
!SMS$SERIAL END
!SMS$SERIAL BEGIN
  read(unitno, err=90) ffhh2d
!SMS$SERIAL END
!
! ice fractions! fice is also used as something different (cloud ice fraction!!)
!
!SMS$SERIAL BEGIN
  read(unitno, err=90) hice2d
!SMS$SERIAL END
!SMS$SERIAL BEGIN
  read(unitno, err=90) fice2d
!SMS$SERIAL END
!SMS$SERIAL BEGIN
  read(unitno, err=90) tprcp2d
!SMS$SERIAL END
!SMS$SERIAL BEGIN
  read(unitno, err=90) srflag2d
!SMS$SERIAL END
!SMS$SERIAL BEGIN
  read(unitno, err=90) snwdph2d
!SMS$SERIAL END
!SMS$SERIAL BEGIN
  read(unitno, err=90) slc2d
!SMS$SERIAL END
!SMS$SERIAL BEGIN
  read(unitno, err=90) shdmin2d
!SMS$SERIAL END
!SMS$SERIAL BEGIN
  read(unitno, err=90) shdmax2d
!SMS$SERIAL END
!SMS$SERIAL BEGIN
  read(unitno, err=90) slope2d
!SMS$SERIAL END
!SMS$SERIAL BEGIN
  read(unitno, err=90) snoalb2d
!SMS$SERIAL END
!
! this thing contains such things as origraphuc stand deviation, convexity, asymetry,...for grav wave drag calc
!
  do idx = 1,SIZE(hprm2d,1)
!SMS$SERIAL BEGIN
    read(unitno, err=90) work2d
!SMS$SERIAL END
    do ipn=1,nip
      hprm2d(idx,ipn) = work2d(ipn)
    enddo
  enddo
!SMS$SERIAL BEGIN
  close(unitno)
  call returnunit (unitno)
!SMS$SERIAL END
  call IncrementTimer(t0,t1)
!SMS$PARALLEL END
  print"(' PHYSICS INPUT time:',F10.0)",t1
end if     ! .not. readrestart

call gfuncphys ()   ! GFS physics

!TODO:  encapsulate this in a new subroutine
! Allocate GFS internal state
! (only on first call, after digital filter it is
!  already allocated)
if (.not. associated(gis_phy)) then
allocate( gis_phy )
! set memory bounds from SMS distribution (decomposition)
mylb = LBOUND(work2d,1)
myub = UBOUND(work2d,1)
gis_phy%ims = mylb
gis_phy%ime = myub
! set patch (distributed-memory loop) bounds from SMS distribution
!SMS$PARALLEL(dh, ipn) BEGIN
first_loop = .true.
do ipn=1,nip
  if (first_loop) then
    gis_phy%ips = ipn
    first_loop = .false.
  endif
  gis_phy%ipe = ipn
enddo
!SMS$PARALLEL END
gis_phy%lsoil = 4
!TODO:  fill this in...  
! This will declare surface and flux arrays as (mylb:myub,1), etc.
! Note:  GFS ignores ierr...  
call sfcvar_aldata(mylb, myub, 1, gis_phy%lsoil, gis_phy%sfc_fld, ierr)
call flxvar_aldata(mylb, myub, 1, gis_phy%flx_fld, ierr)
gis_phy%NBLCK = 1
gis_phy%LEVS = nvl
gis_phy%num_p2d = 3
gis_phy%num_p3d = num_p3d
gis_phy%ras = ras
gis_phy%NMTVR = 14
gis_phy%lats_node_r = 1
!TODO:  move nfxr to module resol_def
gis_phy%nfxr = 14
allocate( gis_phy%SLAG(mylb:myub,gis_phy%LATS_NODE_R),                   &
          gis_phy%SDEC(mylb:myub,gis_phy%LATS_NODE_R),                   &
          gis_phy%CDEC(mylb:myub,gis_phy%LATS_NODE_R),                   &
          gis_phy%SWH(mylb:myub,gis_phy%LEVS,gis_phy%NBLCK,              &
                      gis_phy%LATS_NODE_R),                              &
          gis_phy%HLW(mylb:myub,gis_phy%LEVS,gis_phy%NBLCK,              &
                      gis_phy%LATS_NODE_R),                              &
          gis_phy%SFALB(mylb:myub,gis_phy%LATS_NODE_R),                  &
          gis_phy%ACV(mylb:myub,gis_phy%LATS_NODE_R),                    &
          gis_phy%ACVT(mylb:myub,gis_phy%LATS_NODE_R),                   &
          gis_phy%ACVB(mylb:myub,gis_phy%LATS_NODE_R),                   &
          gis_phy%phy_f2d(mylb:myub,gis_phy%LATS_NODE_R,                 &
                      gis_phy%num_p2d),                                  &
          gis_phy%XLON(mylb:myub,gis_phy%LATS_NODE_R),                   &
          gis_phy%XLAT(mylb:myub,gis_phy%LATS_NODE_R),                   &
          gis_phy%HPRIME(gis_phy%NMTVR,mylb:myub,gis_phy%LATS_NODE_R),   &
          gis_phy%phy_f3d(mylb:myub,gis_phy%LEVS,gis_phy%NBLCK,          &
                          gis_phy%LATS_NODE_R,gis_phy%num_p3d),          &
          gis_phy%COSZDG(mylb:myub,gis_phy%LATS_NODE_R),                 &
          gis_phy%CLDCOV(gis_phy%LEVS,mylb:myub,gis_phy%LATS_NODE_R),    &
          gis_phy%FLUXR(gis_phy%nfxr,mylb:myub,gis_phy%LATS_NODE_R) )

allocate( gis_phy%ps(mylb:myub),               &
          gis_phy%dp(mylb:myub,gis_phy%LEVS),  &
          gis_phy%dpdt(mylb:myub,gis_phy%LEVS),&
          gis_phy%p(mylb:myub,gis_phy%LEVS),   &
          gis_phy%u(mylb:myub,gis_phy%LEVS),   &
          gis_phy%v(mylb:myub,gis_phy%LEVS),   &
          gis_phy%t(mylb:myub,gis_phy%LEVS),   &
          gis_phy%q(mylb:myub,gis_phy%LEVS),   &
          gis_phy%oz(mylb:myub,gis_phy%LEVS),  &
          gis_phy%cld(mylb:myub,gis_phy%LEVS) )

endif
!TODO:  Connect this properly, at present values set here are overwritten 
!TODO:  in physics().  
call flx_init(gis_phy%flx_fld, ierr)

!SMS$PARALLEL(dh, ipn) BEGIN
  do ipn=1,nip
    gis_phy%sfc_fld%CV    (ipn,1) = cv2d    (ipn)
    gis_phy%sfc_fld%CVT   (ipn,1) = cvt2d   (ipn)
    gis_phy%sfc_fld%CVB   (ipn,1) = cvb2d   (ipn)
    gis_phy%sfc_fld%SLMSK (ipn,1) = slmsk2d (ipn)
    gis_phy%flx_fld%SFCDSW(ipn,1) = 0.0
    gis_phy%flx_fld%SFCDLW(ipn,1) = 0.0
    gis_phy%sfc_fld%TSEA  (ipn,1) = ts2d    (ipn)
    gis_phy%sfc_fld%STC (:,ipn,1) = st3d  (:,ipn)
    gis_phy%sfc_fld%SHELEG(ipn,1) = sheleg2d(ipn)
    gis_phy%sfc_fld%ZORL  (ipn,1) = zorl2d  (ipn)
    gis_phy%sfc_fld%snoalb(ipn,1) = snoalb2d(ipn)
    gis_phy%HPRIME      (:,ipn,1) = hprm2d(:,ipn)
    gis_phy%sfc_fld%HICE  (ipn,1) = hice2d  (ipn)
    gis_phy%sfc_fld%FICE  (ipn,1) = fice2d  (ipn)
    gis_phy%sfc_fld%TPRCP (ipn,1) = tprcp2d (ipn)
    gis_phy%sfc_fld%SRFLAG(ipn,1) = srflag2d(ipn)
    gis_phy%sfc_fld%SLC (:,ipn,1) = slc3d (:,ipn)
    gis_phy%sfc_fld%SMC (:,ipn,1) = sm3d  (:,ipn)
    gis_phy%sfc_fld%SNWDPH(ipn,1) = snwdph2d(ipn)
    gis_phy%sfc_fld%SLOPE (ipn,1) = slope2d (ipn)
    gis_phy%sfc_fld%SHDMIN(ipn,1) = shdmin2d(ipn)
    gis_phy%sfc_fld%SHDMAX(ipn,1) = shdmax2d(ipn)
    gis_phy%sfc_fld%TG3   (ipn,1) = tg32d   (ipn)
    gis_phy%sfc_fld%VFRAC (ipn,1) = vfrac2d (ipn)
    gis_phy%sfc_fld%CANOPY(ipn,1) = canopy2d(ipn)
    gis_phy%sfc_fld%VTYPE (ipn,1) = vtype2d (ipn)
    gis_phy%sfc_fld%STYPE (ipn,1) = stype2d (ipn)
    gis_phy%sfc_fld%F10M  (ipn,1) = f10m2d  (ipn)
    gis_phy%sfc_fld%FFMM  (ipn,1) = ffmm2d  (ipn)
    gis_phy%sfc_fld%FFHH  (ipn,1) = ffhh2d  (ipn)
    gis_phy%sfc_fld%ALVSF (ipn,1) = alvsf2d (ipn)
    gis_phy%sfc_fld%ALNSF (ipn,1) = alnsf2d (ipn)
    gis_phy%sfc_fld%ALVWF (ipn,1) = alvwf2d (ipn)
    gis_phy%sfc_fld%ALNWF (ipn,1) = alnwf2d (ipn)
    gis_phy%sfc_fld%FACSF (ipn,1) = facsf2d (ipn)
    gis_phy%sfc_fld%FACWF (ipn,1) = facwf2d (ipn)
    ! Note that T2M and Q2M are overwritten before use in do_physics_one_step()
    gis_phy%sfc_fld%T2M   (ipn,1) = t2m2d   (ipn)
    gis_phy%sfc_fld%Q2M   (ipn,1) = q2m2d   (ipn)
    gis_phy%phy_f3d(ipn,:,1,1,:)  = 0.0
    gis_phy%phy_f2d(ipn,1,:)      = 0.0
    gis_phy%CLDCOV(:,ipn,1)       = 0.0
    gis_phy%flx_fld%HFLX(ipn,1)   = 0.0
    gis_phy%flx_fld%EVAP(ipn,1)   = 0.0
    ! Note that UUSTAR is overwritten before use here.  uustar2d is not used.  
    gis_phy%sfc_fld%UUSTAR(ipn,1) = 0.01
  enddo
!SMS$PARALLEL END
  if ( UpdateSST) then
      call sst_init
  endif
print *,'... exiting phy_init'

return
70 write(6,*)'phy_init: error opening a file'
  stop
90 write(6,*)'phy_init: error reading a file'
  stop
end subroutine phy_init

subroutine sst_init
!*********************************************************************
!       Loads the initial variables and constants for the ocean
!       component.  
!       Alexander E. MacDonald  11/27/04
!       J. Lee                  September, 2005
!*********************************************************************

use module_control  ,only: nip,dt,numphr,                 &                           
      glvl,curve,NumCacheBlocksPerPE,CallSST,PhysicsInterval, &
      RadiationInterval,GravityWaveDrag,SSTInterval
use units, only: getunit, returnunit

implicit none

! Local variables

integer :: mype,idx,ierr,mylb,myub
real*8  :: t0,t1=0.0d0

!namelist /PHYSICSnamelist/ PhysicsInterval,RadiationInterval,GravityWaveDrag,ras,num_p3d,SSTInterval

! TODO:  Create new decomp "dhp" and use it for all physics declarations 
! TODO:  and loops.  Split modules that contain variables used by both too!  

!SMS$insert call nnt_me(mype)
call StartTimer(t0)

! TBH:  BEGIN DUPLICATION with dyn_init().  REFACTOR TO REMOVE DUPLICATION
!open (11,file='./FIMnamelist')
!read (11,PHYSICSnamelist)
!close(11)


! TODO:  Create new decomp "dhp" and use it for all physics declarations 
! TODO:  and loops.  Split modules that contain variables used by both too!  

!SMS$insert call nnt_me(mype)
call StartTimer(t0)
call IncrementTimer(t0,t1)

CallSST = max(1,numphr*SSTInterval/3600)
print "(' Call SST every'  ,I27,' timesteps (',I0,' seconds)')", &
          CallSST,nint(CallSST*dt)
call InitSST
return
end subroutine sst_init

subroutine InitSST

use module_control  ,only: nip, yyyymmddhhmm,prev_date,next_date,glvl,curve,nvl,ptop
use module_constants, only : lat,lon
use module_sfc_variables, only : ts2d,sst_prev,sst_next,fice2d,fice2d_prev,fice2d_next,slmsk2d,hice2d
use gfs_physics_internal_state_mod, only: gfs_physics_internal_state, gis_phy
use units, only: getunit, returnunit
!SMS$ignore begin
  USE slint, ONLY: bilinear_init, bl_int
!SMS$ignore end

implicit none

! locals  
  real            , allocatable   :: sst_prev_ll(:,:)
  real            , allocatable   :: fice2d_prev_ll(:,:)
  real            , allocatable   :: sst_next_ll(:,:)
  real            , allocatable   :: fice2d_next_ll(:,:)
  integer :: YEAR,MONTH,DAY  
  integer :: IM_OC
  integer :: JM_OC,ipn

  integer :: MIDMON
  integer :: ID
  integer :: M1
  integer :: M2
  integer :: MIDM
  integer :: MIDP
  integer :: N
  integer :: I
  integer :: NRECS

  integer :: THIS_YEAR,MDATE
  integer :: THIS_MONTH
  integer :: time_header(4),SDATE1,SDATE2,DAYS(12)
  real    :: REAL_VAR
  real    :: DX
  real    :: DY,FAC
  data    DAYS /31,28,31,30,31,30,31,31,30,31,30,31/
  logical :: skip
  integer :: unitno
! BEGIN
 READ(UNIT=yyyymmddhhmm(1:4), FMT='(I4)') YEAR
 READ(UNIT=yyyymmddhhmm(5:6), FMT='(I2)') MONTH
 READ(UNIT=yyyymmddhhmm(7:8), FMT='(I2)') DAY

IF ( mod(YEAR,4).EQ.0) THEN
   DAYS(2)=29.0
ENDIF
IF ( YEAR.EQ.1900) THEN
   DAYS(2)=28.0
ENDIF
!SMS$SERIAL (<sst_prev,sst_next,fice2d_prev,fice2d_next,OUT> : default=ignore) BEGIN 
  IM_OC=360.0
  JM_OC=180.0
  allocate(sst_prev_ll(IM_OC,JM_OC))
  allocate(fice2d_prev_ll(IM_OC,JM_OC))
  allocate(sst_next_ll(IM_OC,JM_OC))
  allocate(fice2d_next_ll(IM_OC,JM_OC))
  print*,'calling bilinear_init for SST'

  unitno = getunit ()
  if (unitno < 0) then
    print*,'initsst: getunit failed: stopping'
    stop
  end if

  open (unitno, file='glvl.dat', form='unformatted', err=70)
  call TestGlvlHeader (unitno, 'glvl.dat', 'sst_init', glvl)
  call TestCurveHeader(unitno, 'glvl.dat', 'sst_init', curve)
  CALL bilinear_init('ocean_bcs_ltln', IM_OC*JM_OC, unitno, nip)
  close(unitno)
  call returnunit (unitno)
 
  MDATE=YEAR*100+MONTH
  MIDMON = DAYS(MONTH)/2 + 1

  sstunit = getunit ()
  if (sstunit < 0) then
    print*,'initsst: getunit failed for sstunit: stopping'
    stop
  else
    print*,'initsst: sst unit is ', sstunit
  end if

  open (sstunit, file='sst_dat', form='unformatted', err=70)
  skip=.TRUE.
  do while (skip)
     read (sstunit, err=90) time_header
     print*,MDATE,time_header
     read(sstunit, err=90) sst_prev_ll
     read(sstunit, err=90) fice2d_prev_ll
     SDATE1=time_header(1)*100+time_header(2)
     SDATE2=time_header(3)*100+time_header(4)
!     IF (DAY .LT.MIDMON.AND.MDATE.LE.SDATE2) THEN
     IF (MDATE.LE.SDATE2) THEN
!       SDATE2 needs to equal MDATE for sst_prev_ll
        skip=.FALSE.
     ENDIF
  end do
 IF (DAY .GE.MIDMON) THEN
     read(sstunit, err=90) time_header
     read(sstunit, err=90) sst_prev_ll
     read(sstunit, err=90) fice2d_prev_ll
 ENDIF
  prev_date=time_header
  read(sstunit, err=90) time_header 
  read(sstunit, err=90) sst_next_ll
  read(sstunit, err=90) fice2d_next_ll
  next_date=time_header
  
  print*,'FOUND SSTs, calling bl_int',prev_date,next_date

  CALL bl_int (sst_prev_ll(:,:), sst_prev)
  CALL bl_int (sst_next_ll(:,:), sst_next)
  CALL bl_int (fice2d_prev_ll(:,:), fice2d_prev)
  CALL bl_int (fice2d_next_ll(:,:), fice2d_next)

 IF (DAY < MIDMON) THEN

    M1   = MOD(MONTH+10,12) + 1
    M2   = MONTH
    MIDM = DAYS(M1)/2 + 1
    MIDP = DAYS(M1) + MIDMON
    ID = DAY + DAYS(M1)

  ELSE

    M2   = MOD(MONTH,12) + 1
    M1   = MONTH
    MIDM = MIDMON
    MIDP = DAYS(M2)/2 + 1 + DAYS(M1)
    ID = DAY

  ENDIF
!SMS$SERIAL END
!SMS$SERIAL (<FAC,out> : default=ignore) BEGIN
  FAC = (real(ID -   MIDM)*86400)/ &
        (real(MIDP - MIDM)*86400         )
!SMS$SERIAL END
! replace ts2d over ocean points
!SMS$PARALLEL(dh, ipn) BEGIN
   DO ipn=1,nip
      !  need logic to keep ice's temperature the same, only update sst and ice fraction
      ! update ice fraction 1st
      ! if there is new ice, set it to -1.8 and hice=0.0
      ! if ice melts, then set ts2d to sst and hice=0.0
      IF (slmsk2d(ipn).NE.1) THEN
         fice2d(ipn)=fice2d_next(ipn)*(FAC)+fice2d_prev(ipn)*(1.0 - FAC)
         if (fice2d(ipn) .GT. 1) fice2d(ipn)=1.0
	 if (fice2d(ipn) .LT. 0) fice2d(ipn)=0.0
      ENDIF
      IF (fice2d(ipn) .GE. 0.5 .AND. slmsk2d(ipn) .EQ. 0) THEN ! freeze open ocean
         slmsk2d(ipn)=2.0
         ts2d(ipn)=271.35
         hice2d(ipn)=0.0
      ENDIF
      IF (fice2d(ipn) .LT. 0.5 .AND. slmsk2d(ipn) .EQ. 2) THEN ! melt sea-ice
         slmsk2d(ipn)=0.0
         hice2d(ipn)=0.0
      ENDIF
      IF (slmsk2d(ipn).EQ.0) THEN
         ts2d(ipn)=sst_next(ipn)*(FAC)+sst_prev(ipn)*(1.0 - FAC)
      ENDIF

      gis_phy%sfc_fld%TSEA  (ipn,1) = ts2d    (ipn)
      gis_phy%sfc_fld%HICE  (ipn,1) = hice2d  (ipn)
      gis_phy%sfc_fld%FICE  (ipn,1) = fice2d  (ipn)
      gis_phy%sfc_fld%SLMSK (ipn,1) = slmsk2d (ipn)
   ENDDO
!SMS$PARALLEL END

  RETURN

70 write(6,*) 'initsst: error opening file for reading'
  stop
90 write(6,*) 'initsst: error reading file'
  stop
end subroutine InitSST
end module module_fim_phy_init
