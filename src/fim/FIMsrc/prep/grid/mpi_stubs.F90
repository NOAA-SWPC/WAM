
!
! Stubs to attmept work around a ugly build problem...  
! Probably won't work!  
!

subroutine mpi_initialized(initialized, istatus)
  logical, intent(inout) :: initialized
  integer, intent(inout) :: istatus
  initialized = .false.
  istatus = 0
end subroutine mpi_initialized

subroutine mpi_abort(comm,istatus,ignore)
  integer, intent(inout) :: comm,istatus,ignore
  istatus = 0
end subroutine mpi_abort

