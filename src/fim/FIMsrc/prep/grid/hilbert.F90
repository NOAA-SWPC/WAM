!======================================================
!  This subroutine computes the hilbert 2d filling
!  curve. 
!
! orientation 0:
!       ________
!  |   |        |
!  |___|     ___|
!           |
!   ___     |___
!  |   |        |
!  |   |________|
! orientation 1, 2,and 3 rotates 90 deg each clockwise.
!
!  Author: Ning Wang,   Oct., 2006
!======================================================

SUBROUTINE hilbert_curve(n, ort, start_gc)
     USE DataStru
     IMPLICIT NONE

     INTEGER :: n, ort, start_gc
     INTEGER :: level

     xdim = SQRT(REAL(n)) + 0.5
     ydim = SQRT(REAL(n)) + 0.5
     level = LOG(REAL(xdim)) / LOG(2.0)  + 0.5

     gc = start_gc
     offset = start_gc - 1
     orient = ort
     CALL hilbert(0.0, 0.0, REAL(xdim), 0.0, 0.0, REAL(ydim),level)

END SUBROUTINE hilbert_curve



RECURSIVE SUBROUTINE hilbert(x0, y0, xis, xjs, yis, yjs, level)
     USE DataStru
     IMPLICIT NONE

     REAL :: x0, y0, xis, xjs, yis, yjs
     INTEGER :: level
     INTEGER :: i, j, tmp
     IF (level == 0) THEN
       i = x0 + (xis + yis) / 2
       j = y0 + (xjs + yjs) / 2
       IF (orient == 1) THEN
         tmp = i
         i = j
         j = tmp
       END IF
       IF (orient == 2) THEN
         i = xdim - i - 1 
       END IF
       IF (orient == 3) THEN
         tmp = i
         i = j
         j = ydim - tmp - 1
       END IF
       perm(gc) = offset + i + xdim * j + 1
       gc = gc + 1
     ELSE
       CALL hilbert(x0, y0, yis / 2, yjs / 2, xis / 2, xjs / 2, level - 1)
       CALL hilbert(x0 + xis / 2, y0 + xjs / 2, xis / 2, xjs / 2, yis /2, yjs / 2, level - 1)
       CALL hilbert(x0 + xis / 2 + yis / 2, y0 + xjs / 2 + yjs / 2, xis / 2, xjs / 2, yis / 2, yjs / 2, level - 1)
       CALL hilbert(x0 + xis / 2 + yis, y0 + xjs / 2 + yjs, - yis / 2, - yjs / 2, -xis / 2, -xjs / 2, level - 1)
     END IF
END SUBROUTINE hilbert
  

