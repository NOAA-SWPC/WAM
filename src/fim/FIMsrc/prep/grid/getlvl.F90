!
!*********************************************************************
   program getlvl
!  Loads the initial variables and constants to start sgm
!  Alexander E. MacDonald  11/27/05
!  J. Lee                  September, 2005
!  N.W.     Added namelist definition and reading statements. 
!           Changed to dynamic memory allocation for big arrays.
!  N.W.     Changed the computation of nip to work with the icos grid
!           that has mixed bi-section and tri-section subdivisions. 
!  N.W.     Changed to a more accurate computation of the middle points 
!           of the Voronoi-cell edges.
!*********************************************************************
use module_control,only: control,glvl,nip,curve
use read_queue_namelist, only: ReturnNIP
implicit none

integer, parameter :: npp=6
integer, parameter ::  nd=2
real, parameter :: pi = 3.1415926535897
real, parameter :: ae = 6371220.   !earth radius

!
real*4 map
real*4 conr_xy(npp,2,2)
real*4 prox_xy(npp,2)    ! holds x and y locs for prox pts (m)

!
real*4 eltp(4),elnp(4) ! 4 lat/lon surrounding a particular edge
real*4 conr_tmp(1:6,1:2) 

!
real*4, allocatable :: conr_ll(:,:,:,:), lle(:,:,:)
integer, allocatable :: nprox(:), prox(:, :), proxs(:, :), inv_perm(:)
real*4, allocatable :: lat(:),lon(:),area(:)
real*4, allocatable :: sideln(:, :),rprox_ln(:, :)
real*4, allocatable :: sidevec_c(:, :, :),sidevec_e(:, :, :)
real*4, allocatable :: cs(:,:,:),sn(:,:,:)

integer :: ipn, isn, ism, ixy, ipt, iprox, ip1, im1,  j, idx, i2, i3
real :: xlat, xlon, xxp, yyp, xxm, yym, xx, yy, xltc, xlnc, rf
real :: p1(2), p2(2), pm(2)

call control

! allocate memory for arrays
ALLOCATE(conr_ll(npp,2,2,nip))
ALLOCATE(lle(npp,2,nip))
ALLOCATE(nprox(nip))
ALLOCATE(prox(npp, nip))
ALLOCATE(proxs(npp, nip))
ALLOCATE(lat(nip))
ALLOCATE(lon(nip))
ALLOCATE(area(nip))
ALLOCATE(sideln(npp, nip))
ALLOCATE(rprox_ln(npp, nip))
ALLOCATE(sidevec_c(nd, npp, nip))
ALLOCATE(sidevec_e(nd, npp, nip))
ALLOCATE(cs(4,npp,nip))
ALLOCATE(sn(4,npp,nip))
ALLOCATE(inv_perm(nip))

print*, 'start getlvl ... '
!...................................................................

cs       (:,6,:) = 0.0
sn       (:,6,:) = 0.0
sidevec_c(:,6,:) = 0.0
sidevec_e(:,6,:) = 0.0
sideln   (  6,:) = 0.0
rprox_ln (  6,:) = 0.0

OPEN                (10,file='icos_grid_info_level.dat',form='unformatted')
call TestGlvlHeader (10,     'icos_grid_info_level.dat','getlvl',glvl)
call TestCurveHeader(10,     'icos_grid_info_level.dat','getlvl',curve)
READ(10) lat,lon
do isn = 1,size(prox,1)
  READ(10) prox(isn,:)
enddo
READ(10) nprox
do i3 = 1,size(conr_ll,3)
  do i2 = 1,size(conr_ll,2)
    do isn = 1,size(conr_ll,1)
      READ(10) conr_ll(isn,i2,i3,:)
    enddo
  enddo
enddo
READ(10) inv_perm
CLOSE(10)
!
do ipn=1,nip
if(lon(ipn).lt.0.) lon(ipn)=lon(ipn)+2.*pi
end do
!
do ipn=1,nip
!
do isn=1,nprox(ipn)
do ixy=1,2
conr_tmp(isn,ixy)=conr_ll(isn,1,ixy,ipn)
conr_tmp(isn,ixy)=conr_ll(isn,1,ixy,ipn)
end do
end do
!
do isn=1,nprox(ipn)
ism=isn-1
if(isn.eq.1) ism=nprox(ipn)
do ixy=1,2
conr_ll(isn,1,ixy,ipn)=conr_tmp(ism,ixy)
conr_ll(isn,2,ixy,ipn)=conr_tmp(isn,ixy)
end do
!
end do
end do
!
proxs=-99
do ipn=1,nip
do 20 isn=1,nprox(ipn)
iprox=prox(isn,ipn)
do j=1,nprox(iprox)
if(prox(j,iprox).eq.ipn) then
proxs(isn,ipn)=j
goto 20
end if
end do
20 continue
end do
!
!...................................................................
!
do ipn=1,nip
do isn=1,nprox(ipn)
!do ixy=1,nd
!lle(isn,ixy,ipn)=.5*(conr_ll(isn,1,ixy,ipn)+conr_ll(isn,2,ixy,ipn))
!end do
!if ( abs( conr_ll(isn,1,2,ipn)-conr_ll(isn,2,2,ipn)).gt.pi) &
!lle(isn,2,ipn)=lle(isn,2,ipn)-pi

p1(1:2) = conr_ll(isn,1,1:2,ipn)
p2(1:2) = conr_ll(isn,2,1:2,ipn)
call middle_r(p1, p2, pm)
lle(isn, 1:2, ipn) = pm(1:2)

end do
end do
!
!...................................................................
!Caculate sidevec and lat/lon at edges
!
do ipn=1,nip
do isn=1,nprox(ipn)
xlon=lon(prox(isn,ipn))
xlat=lat(prox(isn,ipn))
call ll2xy(lon(ipn),lat(ipn),xlon,xlat,prox_xy(isn,1),prox_xy(isn,2))
rprox_ln(isn,ipn)=1./(ae*sqrt(prox_xy(isn,1)**2+prox_xy(isn,2)**2))
do ipt=1,2
xlon=conr_ll(isn,ipt,2,ipn)
xlat=conr_ll(isn,ipt,1,ipn)
call ll2xy(lon(ipn),lat(ipn),xlon,xlat,conr_xy(isn,ipt,1),conr_xy(isn,ipt,2))
end do
map=2./(1.+sin(lle(isn,1,ipn))*sin(lat(ipn)) &
+cos(lle(isn,1,ipn))*cos(lat(ipn)) &
       *cos(lle(isn,2,ipn)-lon(ipn)))
do ixy=1,nd
sidevec_c(ixy,isn,ipn)=ae*( conr_xy(isn,2,ixy) &
                         -conr_xy(isn,1,ixy)) *map
end do
call ll2xy(lle(isn,2,ipn),lle(isn,1,ipn),conr_ll(isn,2,2,ipn) &
,conr_ll(isn,2,1,ipn),xxp,yyp)
call ll2xy(lle(isn,2,ipn),lle(isn,1,ipn),conr_ll(isn,1,2,ipn) &
,conr_ll(isn,1,1,ipn),xxm,yym)
sidevec_e(1,isn,ipn)= ae*(xxp-xxm)
sidevec_e(2,isn,ipn)= ae*(yyp-yym)
sideln(isn,ipn)=sqrt(sidevec_e(1,isn,ipn)**2+sidevec_e(2,isn,ipn)**2)
end do ! isn loop
area(ipn)=0.   
do isn=1,nprox(ipn)
xx=ae*.5*(conr_xy(isn,2,1)+conr_xy(isn,1,1))
yy=ae*.5*(conr_xy(isn,2,2)+conr_xy(isn,1,2))
area(ipn)=area(ipn)+.5*(xx*sidevec_c(2,isn,ipn)-yy*sidevec_c(1,isn,ipn))
end do
end do ! ipn loop
!
!...................................................................
!
do ipn=1,nip
  do isn=1,nprox(ipn) 
    xltc=lle(isn,1,ipn)
    xlnc=lle(isn,2,ipn)
    ip1=mod(isn,nprox(ipn))+1
    im1=isn-1
    if(im1.eq.0) im1=nprox(ipn)
    eltp(1)=lat(ipn)
    elnp(1)=lon(ipn)
    eltp(2)=lat(prox(isn,ipn))
    elnp(2)=lon(prox(isn,ipn))
    eltp(3)=lat(prox(im1,ipn))
    elnp(3)=lon(prox(im1,ipn))
    eltp(4)=lat(prox(ip1,ipn))
    elnp(4)=lon(prox(ip1,ipn))
    do ipt=1,4
      rf=1.0/(1.0+sin(xltc)*sin(eltp(ipt))+cos(xltc)*cos(eltp(ipt))*cos(elnp(ipt)-xlnc))
      cs(ipt,isn,ipn)=rf*( cos(xltc)*cos(eltp(ipt))+(1.0+sin(xltc)*sin(eltp(ipt)))*cos(elnp(ipt)-xlnc))
      sn(ipt,isn,ipn)=-rf*sin(elnp(ipt)-xlnc)*(sin(xltc)+sin(eltp(ipt)))
    end do
  enddo
enddo

open(unit=28,file="glvl.dat", form="unformatted")
call WriteGlvlHeader (28,glvl )
call WriteCurveHeader(28,curve)
write(28) lat
write(28) lon
write(28) nprox
do isn=1,size(proxs,1)
  write(28) proxs(isn,:)
enddo
do isn=1,size(prox,1)
  write(28) prox(isn,:)
enddo
write(28) area
do isn=1,size(cs,2)
  do idx=1,size(cs,1)
    write(28) cs(idx,isn,:)
  enddo
enddo
do isn=1,size(sn,2)
  do idx=1,size(sn,1)
    write(28) sn(idx,isn,:)
  enddo
enddo
do isn=1,size(sidevec_c,2)
  do idx=1,size(sidevec_c,1)
    write(28) sidevec_c(idx,isn,:)
  enddo
enddo
do isn=1,size(sidevec_e,2)
  do idx=1,size(sidevec_e,1)
    write(28) sidevec_e(idx,isn,:)
  enddo
enddo
do isn=1,size(sideln,1)
  write(28) sideln(isn,:)
enddo
do isn=1,size(rprox_ln,1)
  write(28) rprox_ln(isn,:)
enddo
write(28) inv_perm
close(28)
print*, 'done saving glvl.dat'
!...................................................................
!
stop  
end program getlvl
!
!!!
!
!#############################################################
!    ll2xy.f
!    Convert lat/lon to (x,y) on General Stereographic Coordinate (GSTC).
!    Original program:  J.Lee - 2004
!    Program testing:   J.Lee - 2004
!    Modified for Non-Structure Grid:  J.Lee - 2004
!############################################################

!    Purpose:  Given latitude and longitude on Spherical coordinate,
!              this subroutine computes X and Y coordinates on GSTC.
!    Reference: J.Lee, G. Browning, and Y. Xie: 
!               TELLUS (1995), p.892-910.
!               
! Input Variables : Angles are assumed in unit of "radian" 
! 
!     (latc,lonc) : the GSTC projected point.
!     ( lat, lon) : Input lat/lon in radians.
!
! OUTPUT Variables: 
!
!             xm  : X-Coordinate values on GSTC.  
!                   positive to East of central longitude
!             ym  : Y-Coordinate values on GSTC.  
!                   positive to North of central latitude.
! Note:  Output variables of xm and ym are 
!               nondimensionalized with "ae", the radius of earth.
!
subroutine ll2xy(lonc,latc,lon,lat,xm,ym)
!
implicit none

integer i
real*4 lonc,latc,lon,lat,mf
real*4 xm,ym
!
mf=2.0/(1.0+sin(lat)*sin(latc)+cos(lat)*cos(latc) &
       *cos(lon-lonc))
xm=mf*(cos(lat)*sin(lon-lonc))
ym=mf*((sin(lat)*cos(latc)-cos(lat) &
        *sin(latc)*cos(lon-lonc)) )
!
return
end

!======================================================
!  This subroutine computes the latitude and longitude 
!  of the middle point between two given ponits. The 
!  subroutine is similar to what is in the middle.F90,
!  except that its input and output variables have 
!  single precision and use radians for lat-lon values.
!
!   Ning Wang,   March, 2006
!======================================================

SUBROUTINE middle_r(p1,p2,p)

     IMPLICIT NONE

     REAL :: pi, d2r

     ! Two given points in lat/lon:
     REAL :: p1(2),p2(2),p(2)

     REAL :: xyz1(3),xyz2(3),xyz(3)

     pi = atan(1.0) * 4.0
     ! Convert them into Cardesian coor:
     xyz1(1) = cos(p1(1)) * cos(p1(2))
     xyz1(2) = cos(p1(1)) * sin(p1(2))
     xyz1(3) = sin(p1(1))

     xyz2(1) = cos(p2(1)) * cos(p2(2))
     xyz2(2) = cos(p2(1)) * sin(p2(2))
     xyz2(3) = sin(p2(1))

     ! middle point:

     xyz = 0.5 * (xyz1 + xyz2)
     xyz = xyz / sqrt(xyz(1)*xyz(1) + xyz(2)*xyz(2) + xyz(3)*xyz(3))

     ! Convert the middle point to lat/lon coor:
     p(1) = atan2(xyz(3), sqrt(xyz(1) * xyz(1) + xyz(2) * xyz(2)))
     p(2) = atan2(xyz(2), xyz(1))

     IF (p(2) < 0.0) THEN
       p(2) = p(2) + 2 * pi
     END IF

END SUBROUTINE middle_r


