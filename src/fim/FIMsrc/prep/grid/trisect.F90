!======================================================================
!  This subroutine computes the latitude and longitude of the one-third 
!  and two-third points between two given ponits, on the unit sphere.
!
!  From a general m-sect formula:
!
!  xyz = sin((1-f)*theta) / sin(theta) * xyz1 +
!        sin(f*theta) /sin(theta) * xyz2 ;
!  where theta is the angle between unit vector xyz1 and xyz2.
!
!  xyz_1 = (xyz1 * sin(2/3*theta) + xyz2 * sin(1/3*theta)) / sin(theta)
!  xyz_2 = (xyz1 * sin(1/3*theta) + xyz2 * sin(2/3*theta)) / sin(theta)
!
!  Author: Ning Wang,   June, 2009
!======================================================================

SUBROUTINE trisect(p1, p2, p_1, p_2)
   !  USE module_constants 

     IMPLICIT NONE

     REAL*8 :: pi, d2r

     ! Two given points in lat/lon:
     REAL*8 :: p1(2),p2(2),p_1(2),p_2(2)

     REAL*8 :: xyz1(3),xyz2(3),xyz_1(3), xyz_2(3)
     REAL*8 :: theta, sin_theta, sin_1thetaov3, sin_2thetaov3

     pi = 4.0 * ATAN(1.0) 
     d2r = pi / 180.0

     ! Convert to radian:
     p1 = p1 * d2r
     p2 = p2 * d2r

     ! Convert them into Cardesian coor:
     xyz1(1) = cos(p1(1)) * cos(p1(2))
     xyz1(2) = cos(p1(1)) * sin(p1(2))
     xyz1(3) = sin(p1(1))

     xyz2(1) = cos(p2(1)) * cos(p2(2))
     xyz2(2) = cos(p2(1)) * sin(p2(2))
     xyz2(3) = sin(p2(1))

     theta = acos(cos(p1(1)) * cos(p2(1)) * cos(p1(2) - p2(2)) + sin(p1(1)) * sin(p2(1)))
     sin_theta = sin(theta) 
     sin_1thetaov3 = sin(theta / 3.0) 
     sin_2thetaov3 = sin(theta / 1.5) 
     ! the one-third point:
     xyz_1 = (xyz1 * sin_2thetaov3 + xyz2 * sin_1thetaov3) / sin_theta
     ! the two-third point:
     xyz_2 = (xyz1 * sin_1thetaov3 + xyz2 * sin_2thetaov3) / sin_theta

     ! Convert to lat/lon coor:
     p_1(1) = atan2(xyz_1(3), sqrt(xyz_1(1) * xyz_1(1) + xyz_1(2) * xyz_1(2))) 
     p_1(2) = atan2(xyz_1(2), xyz_1(1)) 
     IF (p_1(2) < -pi / 2.0) THEN 
       p_1(2) = p_1(2) + 2 * pi
     END IF

     p_2(1) = atan2(xyz_2(3), sqrt(xyz_2(1) * xyz_2(1) + xyz_2(2) * xyz_2(2))) 
     p_2(2) = atan2(xyz_2(2), xyz_2(1)) 
     IF (p_2(2) < -pi / 2.0) THEN 
       p_2(2) = p_2(2) + 2 * pi
     END IF

     p1 = p1 / d2r
     p2 = p2 / d2r
     p_1 = p_1 / d2r
     p_2 = p_2 / d2r
    
END SUBROUTINE trisect

