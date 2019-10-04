module infnan
!SMS$IGNORE BEGIN
  implicit none

  private
  public :: inf, nan, negint

! TODO: add big endian ifdef for e.g. IBM
! TODO: Should NaN be signaling or non-signaling?

  real,    parameter :: inf = z'7f800000'   ! Infinity
  real,    parameter :: nan = z'ffc00000'   ! NaN
  integer, parameter :: negint = -999       ! Bad integer value
!SMS$IGNORE END
end module infnan
