program get_gribout
  use postdata, only: post_read_namelist, fimout, gribout

  implicit none

  integer :: ioerr
  integer :: ret

! After calling post_read_namelist, its public variables (fimout, gribout) are available for use

  call post_read_namelist (ret)

  if (ret < 0) then
    write(*,*) 'get_gribout: failure from post_read_namelist'
    stop 999
  end if

  if (gribout) then
    write(*,'(a)') 'gribout:TRUE'
  else
    write(*,'(a)') 'gribout:FALSE'
  end if

  if (fimout) then
    write(*,'(a)') 'fimout:TRUE'
  else
    write(*,'(a)') 'fimout:FALSE'
  end if
  stop
end program get_gribout
