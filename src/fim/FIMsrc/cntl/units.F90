!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! module units contains these public entities
!   getunit: obtain an available fortran unit number
!   returnunit: return a fortran unit number to the list of available units
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module units
!SMS$IGNORE BEGIN
  implicit none

  private
  public :: getunit    ! obtain a Fortran unit number for use
  public :: returnunit ! return a Fortran unit number to the pot

  integer, parameter :: maxunits = 99 ! Not all Fortran compilers allow more than 2 digits for units
  integer :: i                        ! index for "isinuse" array which follows
  logical, save :: isinuse(maxunits) = (/(.false.,i=1,maxunits)/) ! internal state of unit assignments
! The following units will be given only to callers that specifically ask for them, e.g. for situations
! where invoking scripts need to specify special handling via unit numbers for byte swapping
  integer, parameter :: mustask(4) = (/11,21,30,82/)

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! getunit: provide the caller a unit number they can use for I/O
!   Arguments:
!     unitno: Optional: If present, see if the requested unit number is available. If not available,
!             return an error code.
!   Return value: unit number to use, or -1 if a unit number cannot be provided
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function getunit (unitno)
    integer, intent(in), optional :: unitno  ! requested unit number (if present)

    integer :: i   ! index over unit numbers

    getunit = -1   ! initialize to bad return value

! If optional argument "unitno" is present, give the requestor that unit if it is available.

    if (present (unitno)) then
      if (unitno > maxunits .or. unitno < 1 .or. unitno == 5 .or. unitno == 6) then
        write(6,*) 'getunit: Unit ', unitno, ' is not valid'
        return
      end if
      
      if (isinuse (unitno)) then
        write(6,*) 'getunit: Unit ', unitno, ' is already in use'
        return
      end if
        
      isinuse (unitno) = .true.
      getunit = unitno
      return
    end if

! Don't allocate units 5 (stdin), 6 (stdout), or any of the special units normally reserved
! for byte-swapping.

    do i=1,maxunits
      if (.not. isinuse (i) .and. i /= mustask(1) .and. i /= mustask(2) .and. &
          i /= mustask(3) .and. i /= mustask(4) .and. i /= 5 .and. i /= 6) then
        isinuse(i) = .true.
        getunit = i
        return
      end if
    end do

    write(6,*) 'getunit: No more Fortran unit numbers available!'
    return
  end function getunit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! returnunit: return unitno to the pile of available units
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine returnunit (unitno)
    integer, intent(in) :: unitno

    if (unitno > maxunits .or. unitno < 1) then
      write(6,*) 'returnunit: Unit ', unitno, ' is not valid'
      return
    end if

    if (.not. isinuse(unitno)) then
      write(6,*) 'returnunit: WARNING--unit ', unitno, ' is not in use'
      return
    end if

    isinuse(unitno) = .false.
    return
  end subroutine returnunit
!SMS$IGNORE END
end module units
