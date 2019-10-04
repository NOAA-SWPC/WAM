!==================================================================
! The two subroutines in this file rotate the top icosahedral
! grid for the given amount of angles, in theta and lambda directions
!
! Ning Wang,   Jan. 2008
!
!===================================================================
SUBROUTINE rotate(theta,lambda,topgrid, n)
     USE DataStru

     IMPLICIT NONE

     REAL*8 :: theta, lambda
     TYPE(GridPointWnb) :: topgrid(n)
     INTEGER :: n

     REAL*8 :: pi, d2r, cos_t, sin_t, cos_l, sin_l
     REAL*8 :: latlon(2)
     INTEGER :: i, j

     ! Degree-Radian:
     pi = 4.0 * ATAN(1.0) 
     d2r = pi / 180.0

     ! Convert to radian:
     theta = theta * d2r
     lambda = lambda * d2r

     sin_t = SIN(theta)
     cos_t = COS(theta)
     sin_l = SIN(lambda)
     cos_l = COS(lambda)

     DO i = 1, n
       latlon(1) = topgrid(i)%latlon(1) * d2r
       latlon(2) = topgrid(i)%latlon(2) * d2r
       CALL rotate_latlon(latlon, sin_t, cos_t, sin_l, cos_l)
       topgrid(i)%latlon(1) = latlon(1) / d2r
       topgrid(i)%latlon(2) = latlon(2) / d2r
       DO j = 1, 6
         latlon(1) = topgrid(i)%neighbor(1,j) * d2r
         latlon(2) = topgrid(i)%neighbor(2,j) * d2r
         CALL rotate_latlon(latlon, sin_t, cos_t, sin_l, cos_l)
         topgrid(i)%neighbor(1,j) = latlon(1) / d2r
         topgrid(i)%neighbor(2,j) = latlon(2) / d2r
       ENDDO
     ENDDO

END SUBROUTINE rotate

SUBROUTINE rotate_latlon(latlon, sin_t, cos_t, sin_l, cos_l)
     IMPLICIT NONE
     REAL*8 :: latlon(2), sin_t, cos_t, sin_l, cos_l
     REAL*8 :: xyz(3), xyz1(3), xyz2(3)
     REAL*8 :: pi, c11, c12, c13, c21, c22, c23, c31, c32, c33

     pi = 4.0 * ATAN(1.0) 
     ! Convert them into Cardesian coor:
     xyz(1) = cos(latlon(1)) * cos(latlon(2))
     xyz(2) = cos(latlon(1)) * sin(latlon(2))
     xyz(3) = sin(latlon(1))

     ! rotate about Z axis, theta degrees
     xyz1(1) = cos_t * xyz(1) - sin_t * xyz(2)
     xyz1(2) = sin_t * xyz(1) + cos_t * xyz(2)
     xyz1(3) = xyz(3)
     ! ratate about new axis (sin_t, -cos_t, 0), lambda degrees 
     c11 = cos_l + (1 - cos_l) * sin_t * sin_t
     c12 = -(1 - cos_l) * sin_t * cos_t
     c13 = -sin_l * cos_t
     c21 = c12
     c22 = cos_l + (1 - cos_l) * cos_t * cos_t
     c23 = -sin_l * sin_t
     c31 = -c13
     c32 = -c23
     c33 = cos_l

     xyz2(1) = c11 * xyz1(1) + c12 * xyz1(2) + c13 * xyz1(3)
     xyz2(2) = c21 * xyz1(1) + c22 * xyz1(2) + c23 * xyz1(3)
     xyz2(3) = c31 * xyz1(1) + c32 * xyz1(2) + c33 * xyz1(3)

     ! Convert the grid point back to lat/lon coor:
     latlon(1) = atan2(xyz2(3), sqrt(xyz2(1) * xyz2(1) + xyz2(2) * xyz2(2))) 
     latlon(2) = atan2(xyz2(2), xyz2(1)) 

     IF (latlon(2) < - pi / 2.0) THEN 
       latlon(2) = latlon(2) + 2 * pi 
     END IF

END SUBROUTINE rotate_latlon
    

