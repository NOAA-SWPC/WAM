!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Write 32-bit array to the restart file after gathering from other MPI tasks.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine writearr32 (arr, dim1siz, unitno)
  use module_control, only: nip
  
  implicit none

  integer, intent(in) :: dim1siz  ! size of 1st dimension of arr
  integer, intent(in) :: unitno   ! unit number to write to
    
!SMS$DISTRIBUTE(dh,2) BEGIN
  real*4, intent(in) :: arr(dim1siz,nip)
!SMS$DISTRIBUTE END

! The <arr,in> keeps SMS from generating a needless scatter after the write
!SMS$SERIAL (<ARR,IN> : DEFAULT=IGNORE) BEGIN
  write (unitno, err=90) arr
!  write(6,*)'writearr32: arr(1,1)=',arr(1,1),' arr(dim1siz,nip)=',arr(dim1siz,nip)
!SMS$SERIAL END

  return

90 write(6,*) 'writearr32: error writing to unit ', unitno, ' Stopping'
  call flush(6)
  stop
end subroutine writearr32
