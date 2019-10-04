module module_abstart
contains
!#############################################################
!    abstart.f90
!    "Adams-Bashforth start"
!    Original program:  A. E. MacDonald - August 1991
!    J. Lee  - September, 2005
!    Purpose:  This subroutine taken from QNH model initiates 
!    the third order Adams-Bashforth time differencing.  
!    Note that  although it is the third order scheme of A-B,
!    its accuracy is fourth order in time.
!#############################################################

subroutine abstart (its,itsStart,   &
  nf,of,vof,                &	! Adams Bash time dif indices
  adbash1,adbash2,adbash3)			! constants for Adams Bashforth

use module_control, only: dt

implicit none

integer,intent (IN)    :: its			! model time step
integer,intent (IN)    :: itsStart              ! first model time step
integer,intent (INOUT) :: nf,of,vof
real   ,intent (INOUT) :: adbash1,adbash2,adbash3

!.............................................................
!    Sec 1.  Adams-Bashforth Load
!.............................................................

!     Adams-Bashforth indexes
nf = nf + 1
if (nf == 4) nf = 1
of = of + 1
if (of == 4) of = 1
vof = vof + 1
if (vof == 4) vof = 1

!JR Should be able to do it this way, but be CAREFUL that it also works for digifilt=.t. in run.F90
!  vof = of
!  of = nf
!  nf = nf + 1
!  if (nf == 4) nf = 1

! Third order Adams-Bashforth must use Forward differencing 
! on the first time step, and Second Order AB on the second

if (its == 1)then
  adbash1 = dt
  adbash2 = 0.
  adbash3 = 0.
endif

if (its == 2)then
  adbash1 = dt*1.5
  adbash2 = -dt*.5
  adbash3 = 0.
endif

if (its == 3)then
  adbash1 = dt*23./12.
  adbash2 = -dt*16./12.
  adbash3 = dt*5./12.
endif

return
end subroutine abstart
end module module_abstart
