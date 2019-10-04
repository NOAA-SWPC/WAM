!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Used by restart.F90 to create restart file name from current time
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

character(len=8) function its2string (its)
  use module_control, only: ArchvTimeUnit, dt, hrs_in_month

  implicit none

  integer, intent(in) :: its
  integer :: time

! Changes in this block may require changes in dyn_init.F90 as well
  if (ArchvTimeUnit == 'ts') then
    time = its
    write (its2string,'(i6.6,a2)') time, 'ts'
  else if (ArchvTimeUnit == 'mi') then
    time = its*dt/60.
    write (its2string,'(i6.6,a2)') time, 'mi'
  else if (ArchvTimeUnit == 'hr') then
    time = its*dt/3600.
    write (its2string,'(i6.6,a2)') time, 'hr'
  else if (ArchvTimeUnit == 'dy') then
    time = its*dt/86400.
    write (its2string,'(i6.6,a2)') time, 'dy'
  else if (ArchvTimeUnit == 'mo') then
    time = its*dt/(hrs_in_month*3600.)
    write (its2string,'(i6.6,a2)') time, 'mo'
  else
    write (*,'(a,a)') 'ERROR in its2string: unrecognized output time unit: ',ArchvTimeUnit
    stop
  endif

  return
end function its2string
