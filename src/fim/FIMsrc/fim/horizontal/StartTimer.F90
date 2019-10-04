subroutine StartTimer(t0)
INCLUDE "mpif.h"
real*8,intent(out) :: t0
t0 = mpi_wtime()
return
end subroutine StartTimer
