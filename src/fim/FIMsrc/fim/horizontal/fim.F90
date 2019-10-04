!*********************************************************************
   program fim
!  icosahedral flow-following model
!  Authors: Alexander E. MacDonald & Jin-Luen Lee  11/12/05
!  Lead Developer:  J. L. Lee
!  Design Team:  J.L. Lee, R. Bleck, A. E. MacDonald, S. Benjamin
!  Computational design: A. E. MacDonald, J. Middlecoff, D. Schaffer
!*********************************************************************

implicit none

real*8 :: t0,t1=0.0d0

!TODO:  Strictly speaking, the MPI timers called from StartTimer should not 
!TODO:  be called prior to MPI_INIT.  For a serial build, MPI_INIT is never 
!TODO:  called.  In practice this has not been a problem.  If it becomes 
!TODO:  a problem fix it by using a different timer inside StartTimer and 
!TODO:  IncrementTimer for a serial build.  

#ifdef MANUALGPTL
!JR GPTL timers enabled, without iargc/getarg support. Therefore need
!JR to invoke initialize/start/stop/print functions manually
#include <gptl.inc>
integer :: ret, mype
call gptlprocess_namelist ('FIMnamelist', 77, ret)
ret = gptlinitialize ()
ret = gptlstart ('main')
#endif

call StartTimer(t0)

! NOTE:  Executable SMS directives must not be placed before this call!  
! NOTE:  This includes writes or prints without !SMS$ignore because they
!        generate SMS code.

call init
call run

#ifdef MANUALGPTL
!JR GPTL timers enabled, without iargc/getarg support. Therefore need
!JR to invoke stop/print functions manually
!SMS$insert call nnt_me(mype)
ret = gptlstop ('main')
ret = gptlpr (mype)
#endif

call finalize
call IncrementTimer(t0,t1)

print*,'Total time =', t1

end program fim
