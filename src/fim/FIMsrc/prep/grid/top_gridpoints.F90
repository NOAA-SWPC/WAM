!====================================================
!  This routine sets up the initial grid points for a
!  Icosahedron grid.
!
!  Author: Yuanfu Xie on Feb. 2002
!====================================================

SUBROUTINE set_top_gridpoints(top)

  USE DataStru

  IMPLICIT NONE

!  REAL*8, PARAMETER :: pi = atan(1.0) * 4.0
!  REAL*8, PARAMETER :: Lat0 = 90.0 - 2.0 *acos(0.5/sin(pi/5.0))*180.0/pi, Lon0 = -72.0
!  REAL*8, PARAMETER :: Lat0 = 26.56505118, Lon0 = -72.0

  REAL*8, PARAMETER :: Lat0 = 26.5650428671016, Lon0 = -72.0

  TYPE(GridPointWnb) :: top(12)

  ! Local:
  INTEGER :: i,n

  print*, 'lat0=', lat0
  DO i = 1,12
    top(i)%Neighbor(1:2,6) = 0.0
  END DO

  ! North pole:
  top(1)%LatLon(1) = 90.0
  top(1)%LatLon(2) = 0.0
  ! 5 neighbors of this pole:
  DO i=1,5
     top(1)%Neighbor(1,i) = Lat0
     top(1)%Neighbor(2,i) = Lon0+72.0*i
  END DO

  ! 5 neighbors:
  DO i=2,6
     top(i)%LatLon(1) = top(1)%Neighbor(1,i-1)
     top(i)%LatLon(2) = top(1)%Neighbor(2,i-1)

     ! North pole is one of its neighbor:
     top(i)%Neighbor(1,1) = top(1)%LatLon(1)
     top(i)%Neighbor(2,1) = top(1)%LatLon(2)

     ! Left neighbor:
     n = i-1
     if (n == 1) n = 5
     top(i)%Neighbor(1,2) = top(1)%Neighbor(1,n)
     top(i)%Neighbor(2,2) = top(1)%Neighbor(2,n)
     ! Right neighbor:
     n = i+1
     if (n == 7) n = 2
     top(i)%Neighbor(1,3) = top(1)%Neighbor(1,n)
     top(i)%Neighbor(2,3) = top(1)%Neighbor(2,n)
  END DO

  ! South pole:
  top(7)%LatLon(1) = -90.0
  top(7)%LatLon(2) = 0.0
  ! 5 neighbors:
  DO i=1,5
     top(7)%Neighbor(1,i) = -Lat0
     top(7)%Neighbor(2,i) = Lon0+36.0+72.0*i
  END DO

  ! 5 neighbors:
  DO i=8,12
     top(i)%LatLon(1) = top(7)%Neighbor(1,i-7)
     top(i)%LatLon(2) = top(7)%Neighbor(2,i-7)

     ! South pole is one of its neighbor:
     top(i)%Neighbor(1,1) = top(7)%LatLon(1)
     top(i)%Neighbor(2,1) = top(7)%LatLon(2)

     ! Left neighbor:
     n = i-1
     if (n == 7) n = 12
     top(i)%Neighbor(1,2) = top(7)%Neighbor(1,n-7)
     top(i)%Neighbor(2,2) = top(7)%Neighbor(2,n-7)
     ! Right neighbor:
     n = i
     if (n == 12) n = 8
     top(i)%Neighbor(1,3) = top(7)%Neighbor(1,n-7)
     top(i)%Neighbor(2,3) = top(7)%Neighbor(2,n-7)
  END DO

  ! Neighbors cross the equater:
  DO i=2,6

     ! South starts 36 degree ahead:
     n = i+6
     top(i)%Neighbor(1,4) = top(n)%LatLon(1)
     top(i)%Neighbor(2,4) = top(n)%LatLon(2)
     top(n)%Neighbor(1,4) = top(i)%LatLon(1)
     top(n)%Neighbor(2,4) = top(i)%LatLon(2)
     n = i+5
     if (n == 7) n = 12
     top(i)%Neighbor(1,5) = top(n)%LatLon(1)
     top(i)%Neighbor(2,5) = top(n)%LatLon(2)
     top(n)%Neighbor(1,5) = top(i)%LatLon(1)
     top(n)%Neighbor(2,5) = top(i)%LatLon(2)

  END DO

END SUBROUTINE set_top_gridpoints
