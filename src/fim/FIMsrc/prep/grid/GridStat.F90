      PROGRAM GridStat

!==========================================================
!  This program generates the important statistics of the 
!  given unstructed grid.
!
!  Note: all lat/lon are in radians.
!
!  HISTORY: May. 2006 by Ning Wang.
!==========================================================

      IMPLICIT NONE
      INTEGER, PARAMETER :: n=163842
!      INTEGER, PARAMETER :: n=40962
      INTEGER, PARAMETER :: mnc = 6
      INTEGER :: icos_prox(mnc,n), icos_nprox(n)
      REAL    :: icos_grid(2,n) ! grid location (ll)
      REAL    :: icos_edge(mnc,2,2,n) ! end points of edges (ll)
      REAL    :: icos_map(mnc,n) ! mapping factor

      OPEN(unit=10,file='icos_grid_info.dat',form='unformatted')
      READ(10) icos_grid(1,1:n),icos_grid(2,1:n), &
       icos_prox(1:mnc,1:n),icos_nprox(1:n), &
       icos_edge(1:mnc,1:2,1:2,1:n), &
       icos_map(1:mnc,1:n)
      CLOSE(10)
      CALL minmax_values(n, mnc, icos_grid, icos_prox, icos_nprox, icos_edge)

      END PROGRAM GridStat

SUBROUTINE minmax_values(n, mnc, icos_grid, icos_prox, icos_nprox, icos_edge)

     IMPLICIT NONE

     REAL    :: dist
     INTEGER i, j, k
     INTEGER :: n, mnc
     INTEGER :: icos_prox(mnc,n),icos_nprox(n)
     REAL    :: icos_grid(2,n) ! grid location (ll)
     REAL    :: icos_edge(mnc,2,2,n) ! end points of edges (ll)

     REAL max_value, min_value, value

     max_value = 0.0
     min_value = 1000000.0
     DO i = 1, n
       DO j = 1, icos_nprox(i)
         value = dist(icos_grid(1:1, i), icos_grid(2:2, i), icos_grid(1:1, icos_prox(j, i)), icos_grid(2:2, icos_prox(j, i))) 
         min_value = min(value, min_value)
         max_value = max(value, max_value)
       END DO
     END DO
     PRINT*, 'min max distance between neighboring grid points = ', min_value, max_value  
     
     max_value = 0.0
     min_value = 1000000.0
     DO i = 1, n
       DO j = 1, icos_nprox(i)
         value = dist(icos_edge(j, 1, 1, i), icos_edge(j, 1, 2, i), &
           icos_edge(j, 2, 1, i), icos_edge(j, 2, 2, i)) 
         min_value = min(value, min_value)
         max_value = max(value, max_value)
!print*,  icos_edge(j, 1, 2, i), icos_edge(j, 2, 2, i)
       END DO
!print*, '---------------'
     END DO
     PRINT*, 'min max distance of the Voronoi cell edges = ', min_value, max_value  

END SUBROUTINE minmax_values


REAL FUNCTION dist(lat1, lon1, lat2, lon2)

     IMPLICIT NONE

     REAL, INTENT(IN) :: lat1, lon1, lat2, lon2
     REAL             :: x

     x = COS(lat1) * COS(lat2) * COS(lon1 - lon2) + SIN(lat1) * SIN(lat2)
     if(abs(x) >= 1.0) then
       dist = 0.0
     else
       dist = 6371.220 * ACOS(x)
     endif

END FUNCTION dist

