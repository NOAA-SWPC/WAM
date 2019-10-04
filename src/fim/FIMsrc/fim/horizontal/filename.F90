character(len=*) function filename(tag,its)
!sms$ignore begin
  use module_control,only:ArchvTimeUnit
  implicit none
  character(len=*),intent(in)::tag
  integer,intent(in)::its
  character(len=6)::timestr
  integer::its2time
  write (timestr,'(i6.6)') its2time(its)
  filename='fim_out_'//trim(tag)//timestr//ArchvTimeUnit
!sms$ignore end
end function filename

character(len=*) function flexflnm(tag,its)
! --- same as subr.filename but without the hardwired 'fim_out_' part.
! --- file name returned by flexflnm starts with string 'tag'.
!sms$ignore begin
  use module_control,only:ArchvTimeUnit
  implicit none
  character(len=*),intent(in)::tag
  integer,intent(in)::its
  character(len=6)::timestr
  integer::its2time
  write (timestr,'(i6.6)') its2time(its)
  flexflnm=trim(tag)//timestr//ArchvTimeUnit
!sms$ignore end
end function flexflnm
