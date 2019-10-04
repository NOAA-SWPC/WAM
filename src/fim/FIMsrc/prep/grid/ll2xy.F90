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
