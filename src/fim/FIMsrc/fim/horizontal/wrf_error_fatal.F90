
!TODO:  move this into a wrf-specific directory

   subroutine wrf_error_fatal(msg)
   implicit none
   character,intent(IN) :: msg*(*)
   write(0,*) 'Fatal error:  ',TRIM(msg)
   stop
   end subroutine wrf_error_fatal

