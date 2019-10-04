!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read 64-bit array from the restart file and distribute to other MPI tasks.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine readarr64 (arr, dim2siz, unitno)
  use module_control, only: nip
  
  implicit none

  integer, intent(in) :: dim2siz  ! size of 2nd dimension of arr
  integer, intent(in) :: unitno   ! unit number to read from

  integer :: n, ipn               ! loop indices
    
!SMS$DISTRIBUTE(dh,1) BEGIN
  real*8, intent(inout) :: arr(nip,dim2siz)  ! array to be read from restart file
!SMS$DISTRIBUTE END

! The <arr,out> keeps SMS from generating a needless gather before the read
!SMS$SERIAL (<ARR,OUT> : DEFAULT=IGNORE) BEGIN
  read (unitno, err=90) arr
!SMS$SERIAL END
  return

90 write(6,*) 'readarr64: error reading from unit ', unitno, ' Stopping'
  stop
end subroutine readarr64
