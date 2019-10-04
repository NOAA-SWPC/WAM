integer function its2time(its)
  use module_control,only:ArchvTimeUnit,dt,hrs_in_month
  implicit none
  integer,intent(in)::its
! Changes in this block may require changes in dyn_init.F90 as well
  if (ArchvTimeUnit.eq.'ts') then
    its2time=its
  else if (ArchvTimeUnit.eq.'mi') then
    its2time=its*dt/60.
  else if (ArchvTimeUnit.eq.'hr') then
    its2time=its*dt/3600.
  else if (ArchvTimeUnit.eq.'dy') then
    its2time=its*dt/86400.
  else if (ArchvTimeUnit.eq.'mo') then
    its2time=its*dt/(hrs_in_month*3600.)
  else
    write (*,'(a,a)') 'ERROR in its2time: unrecognized output time unit: ',ArchvTimeUnit
    stop
  endif
end function its2time

