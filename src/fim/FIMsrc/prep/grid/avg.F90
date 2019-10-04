SUBROUTINE average(p1,p2,p3,p)
   !  USE module_constants 

     IMPLICIT NONE

     REAL*8 :: pi, d2r

     ! Two given points in lat/lon:
     REAL*8 :: p1(2),p2(2),p3(2),p(2)

     REAL*8 :: xyz1(3),xyz2(3),xyz3(3),xyz(3)

     ! Radian-Degree:
     pi = 4.0 * ATAN(1.0) 
     d2r = pi / 180.0

     ! Convert to radian:
     p1 = p1 * d2r
     p2 = p2 * d2r
     p3 = p3 * d2r

     ! Convert them into Cardesian coor:
     xyz1(1) = cos(p1(1)) * cos(p1(2))
     xyz1(2) = cos(p1(1)) * sin(p1(2))
     xyz1(3) = sin(p1(1))

     xyz2(1) = cos(p2(1)) * cos(p2(2))
     xyz2(2) = cos(p2(1)) * sin(p2(2))
     xyz2(3) = sin(p2(1))

     xyz3(1) = cos(p3(1)) * cos(p3(2))
     xyz3(2) = cos(p3(1)) * sin(p3(2))
     xyz3(3) = sin(p3(1))

     ! middle point:

     xyz = (xyz1 + xyz2 + xyz3) /3.0
     xyz = xyz / sqrt(DOT_PRODUCT(xyz,xyz))

     ! Convert the middle point to lat/lon coor:
     p(1) = atan2(xyz(3), sqrt(xyz(1) * xyz(1) + xyz(2) * xyz(2))) 
     p(2) = atan2(xyz(2), xyz(1)) 

     IF (p(2) < -pi / 2.0) THEN 
       p(2) = p(2) + 2 * pi
     END IF

     p1 = p1 / d2r
     p2 = p2 / d2r
     p3 = p3 / d2r
     p = p / d2r
    
END SUBROUTINE average

