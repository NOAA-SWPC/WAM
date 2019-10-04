module module_GetGrid
contains
!
subroutine GetGrid(npp,nd,glvl,curve,lat,lon,prox,nprox,conr_ll,area,cs,sn,sidevec_c,sidevec_e,sideln,rprox_ln)
!  Loads the initial variables and constants to start sgm
!  Alexander E. MacDonald  11/27/05
!*********************************************************************
use module_control,only: nip
implicit none
!real*8 t0
real, parameter :: pi = 3.1415926535897
real, parameter :: ae = 6371220.   !earth radius

INTEGER,intent (IN)    :: npp,nd,glvl,curve
!SMS$DISTRIBUTE(dh,NIP) BEGIN
real   ,intent (IN)    :: lat(nip)
real   ,intent (INOUT) :: lon(nip)
integer,intent (IN)    :: prox(npp,nip),nprox(nip)
real   ,intent (INOUT) :: conr_ll(npp,2,nip,2)
real   ,intent (OUT)   :: area(nip)
real   ,intent (OUT)   :: cs(4,npp,nip),sn(4,npp,nip)
real   ,intent (OUT)   :: sidevec_c(nd, npp, nip),sidevec_e(nd, npp, nip)
real   ,intent (OUT)   :: sideln(npp, nip),rprox_ln(npp, nip)

real*4                 :: lle(npp,2,nip)
!SMS$DISTRIBUTE END

real*4 map
real*4 conr_xy(npp,2,2)
real*4 prox_xy(npp,2)    ! holds x and y locs for prox pts (m)

real*4 eltp(4),elnp(4) ! 4 lat/lon surrounding a particular edge
real*4 conr_tmp(1:6,1:2) 


integer :: ipn, isn, ism, ixy, ipt, ip1, im1,  j
real :: xlat, xlon, xxp, yyp, xxm, yym, xx, yy, xltc, xlnc, rf

!call StartTimer(t0)
!SMS$PARALLEL(dh, ipn) BEGIN
!SMS$EXCHANGE(prox)
do ipn=1,nip
  if(lon(ipn).lt.0.) lon(ipn)=lon(ipn)+2.*pi
end do

do ipn=1,nip
  do isn=1,nprox(ipn)
    conr_tmp(isn,1)=conr_ll(isn,1,ipn,1)
    conr_tmp(isn,2)=conr_ll(isn,1,ipn,2)
  end do

  do isn=1,nprox(ipn)
    ism=isn-1
    if(isn.eq.1) ism=nprox(ipn)
    do ixy=1,2
      conr_ll(isn,1,ipn,ixy)=conr_tmp(ism,ixy)
      conr_ll(isn,2,ipn,ixy)=conr_tmp(isn,ixy)
    end do
  end do
end do

do ipn=1,nip
  do isn=1,nprox(ipn)
    do ixy=1,nd
      lle(isn,ixy,ipn)=.5*(conr_ll(isn,1,ipn,ixy)+conr_ll(isn,2,ipn,ixy))
    end do
    if ( abs( conr_ll(isn,1,ipn,2)-conr_ll(isn,2,ipn,2)).gt.pi) lle(isn,2,ipn)=lle(isn,2,ipn)-pi
 end do
end do

!Caculate sidevec and lat/lon at edges

!SMS$EXCHANGE(lat,lon)
do ipn=1,nip
  do isn=1,nprox(ipn)
    xlon=lon(prox(isn,ipn))
    xlat=lat(prox(isn,ipn))
    call ll2xy(lon(ipn),lat(ipn),xlon,xlat,prox_xy(isn,1),prox_xy(isn,2))
    rprox_ln(isn,ipn)=1./(ae*sqrt(prox_xy(isn,1)**2+prox_xy(isn,2)**2))
    do ipt=1,2
      xlon=conr_ll(isn,ipt,ipn,2)
      xlat=conr_ll(isn,ipt,ipn,1)
      call ll2xy(lon(ipn),lat(ipn),xlon,xlat,conr_xy(isn,ipt,1),conr_xy(isn,ipt,2))
    end do
    map=2./(1.+sin(lle(isn,1,ipn))*sin(lat(ipn))+cos(lle(isn,1,ipn))*cos(lat(ipn)) &
                                                *cos(lle(isn,2,ipn)-lon(ipn)))
    do ixy=1,nd
      sidevec_c(ixy,isn,ipn)=ae*( conr_xy(isn,2,ixy)-conr_xy(isn,1,ixy)) *map
    end do
    call ll2xy(lle(isn,2,ipn),lle(isn,1,ipn),conr_ll(isn,2,ipn,2),conr_ll(isn,2,ipn,1),xxp,yyp)
    call ll2xy(lle(isn,2,ipn),lle(isn,1,ipn),conr_ll(isn,1,ipn,2),conr_ll(isn,1,ipn,1),xxm,yym)
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
!SMS$PARALLEL END
!call IncrementTimer(t0,t1)
!print*,'GetGrid time      =',t1

return
end subroutine GetGrid

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

mf=2.0/(1.0+sin(lat)*sin(latc)+cos(lat)*cos(latc) &
       *cos(lon-lonc))
xm=mf*(cos(lat)*sin(lon-lonc))
ym=mf*((sin(lat)*cos(latc)-cos(lat) &
        *sin(latc)*cos(lon-lonc)) )
!
return
end subroutine ll2xy
end module module_GetGrid

