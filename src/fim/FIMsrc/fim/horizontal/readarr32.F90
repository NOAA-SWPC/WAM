!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read 32-bit array and distribute to other MPI tasks.
! Used by read_restart_dyn and phy_init.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine readarr32 (arr, dim1siz, unitno)
  use module_control, only: nip
  
  implicit none

  integer, intent(in) :: dim1siz  ! size of 1st dimension of arr
  integer, intent(in) :: unitno   ! unit number to read from
    
!SMS$DISTRIBUTE(dh,2) BEGIN
  real*4, intent(inout) :: arr(dim1siz,nip)
!SMS$DISTRIBUTE END

! The <arr,out> keeps SMS from generating a needless gather before the read
!SMS$SERIAL (<ARR,OUT> : DEFAULT=IGNORE) BEGIN
  read (unitno, err=90) arr
!SMS$SERIAL END

  return

90 write(6,*) 'readarr32: error reading from unit ', unitno, ' Stopping'
  stop
end subroutine readarr32
