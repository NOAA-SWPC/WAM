subroutine IncrementTimer(t0,t1)
INCLUDE "mpif.h"
real*8,intent(in   ) :: t0
real*8,intent(inout) :: t1
t1 = t1 + mpi_wtime()-t0
return
end subroutine IncrementTimer
