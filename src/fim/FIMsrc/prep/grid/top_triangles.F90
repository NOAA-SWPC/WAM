!====================================================
!  This routine sets up the initial triangles for a
!  Icosahedron grid.
!
!  Author: Yuanfu Xie on Feb. 2002
!  
!  Aug 2006, made modifications to fit new data 
!  structure, and changed the orders of the veritces
!  to make them consistent. Changed name of the 
!  subroutine to  set_top_triangle(); 
!  Added some comments.
!====================================================

SUBROUTINE set_top_triangle(high,top)

  USE DataStru

  IMPLICIT NONE

  TYPE(GridPointWnb) :: top(12)
  TYPE(Triangle) :: high(20)

  INTEGER :: i,n

  ! Top 5:
  do i=1,5

     ! first vertex is always the north pole
     high(i)%vertex(1)%LatLon(1) = top(1)%LatLon(1)
     high(i)%vertex(1)%LatLon(2) = top(1)%LatLon(2)
  
     ! second and third is are north pole's first level neighbors
     high(i)%vertex(2)%LatLon(1) = top(1)%Neighbor(1,i)
     high(i)%vertex(2)%LatLon(2) = top(1)%Neighbor(2,i)
     n = i+1
     if (n == 6) n = 1
     high(i)%vertex(3)%LatLon(1) = top(1)%Neighbor(1,n)
     high(i)%vertex(3)%LatLon(2) = top(1)%Neighbor(2,n)
  END DO

  ! Bottom 5:
  do i=6,10

     ! last vertex is always the south pole
     high(i)%vertex(1)%LatLon(1) = top(7)%LatLon(1)
     high(i)%vertex(1)%LatLon(2) = top(7)%LatLon(2)

     ! second and third is are south pole's first level neighbors
     high(i)%vertex(3)%LatLon(1) = top(7)%Neighbor(1,i-5)
     high(i)%vertex(3)%LatLon(2) = top(7)%Neighbor(2,i-5)
     n = i+1
     if (n == 11) n = 6
     high(i)%vertex(2)%LatLon(1) = top(7)%Neighbor(1,n-5)
     high(i)%vertex(2)%LatLon(2) = top(7)%Neighbor(2,n-5)

  END DO

  ! Middle 10:
  do i=1,5
     
     high(i+10)%vertex(1)%LatLon(1) = top(i+1)%LatLon(1)
     high(i+10)%vertex(1)%LatLon(2) = top(i+1)%LatLon(2)

     high(i+10)%vertex(3)%LatLon(1) = top(i+1)%Neighbor(1,4)
     high(i+10)%vertex(3)%LatLon(2) = top(i+1)%Neighbor(2,4)

     high(i+10)%vertex(2)%LatLon(1) = top(i+1)%Neighbor(1,5)
     high(i+10)%vertex(2)%LatLon(2) = top(i+1)%Neighbor(2,5)

  end do
  do i=8,12
     
     high(i+8)%vertex(1)%LatLon(1) = top(i)%LatLon(1)
     high(i+8)%vertex(1)%LatLon(2) = top(i)%LatLon(2)

     high(i+8)%vertex(3)%LatLon(1) = top(i)%Neighbor(1,4)
     high(i+8)%vertex(3)%LatLon(2) = top(i)%Neighbor(2,4)

     high(i+8)%vertex(2)%LatLon(1) = top(i)%Neighbor(1,5)
     high(i+8)%vertex(2)%LatLon(2) = top(i)%Neighbor(2,5)

  end do

END SUBROUTINE set_top_triangle
