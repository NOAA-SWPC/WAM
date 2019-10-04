subroutine read_mtnvar(mdrag3d,imax,jmax,mvar,filename)
  integer,intent(in)::imax,jmax,mvar
  real,intent(out)::mdrag3d(imax,jmax,mvar)
  character(len=*),intent(in)::filename
  integer::i
  open(21,file=trim(filename),form="unformatted",status='old',iostat=i)
  if (i.ne.0) then
    write (*,'(a,a,a)') 'ERROR in read_mtnvar: Could not open ',trim(filename),'.'
    stop
  endif
  read(21,iostat=i) mdrag3d
  if (i.ne.0) then
    write (*,'(a)') 'ERROR in read_mtnvar: Could not read mdrag3d.'
    stop
  endif    
  close(21)
end subroutine read_mtnvar
