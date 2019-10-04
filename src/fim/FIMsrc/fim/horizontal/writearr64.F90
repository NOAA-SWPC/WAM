!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Write 64-bit array to the restart file after gathering from other MPI tasks.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine writearr64 (arr, dim2siz, unitno)
  use module_control, only: nip
  
  implicit none

  integer, intent(in) :: dim2siz  ! size of 2nd dimension of arr
  integer, intent(in) :: unitno   ! unit number to write to
    
!SMS$DISTRIBUTE(dh,1) BEGIN
  real*8, intent(in) :: arr(nip,dim2siz)
!SMS$DISTRIBUTE END

! The <arr,in> keeps SMS from generating a needless scatter after the write
!SMS$SERIAL (<ARR,IN> : DEFAULT=IGNORE) BEGIN
  write (unitno, err=90) arr
!SMS$SERIAL END
  return

90 write(6,*) 'writearr64: error writing to unit ', unitno, ' Stopping'
  stop
end subroutine writearr64
