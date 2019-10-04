      module module_error_msg

      USE ESMF

      implicit none

      private
      public :: err_msg,message_check

      logical, parameter :: iprint = .false.
      character(esmf_maxstr) :: message_check

      contains

      subroutine err_msg(rc1,msg,rc)
      integer, intent(inout)        :: rc1
      integer, intent(out)          :: rc
      character (len=*), intent(in) :: msg
      if (ESMF_LogFoundError(rc1, msg=msg)) then
          rc  = esmf_failure
          print*, 'error happened for ',msg, ' rc = ', rc1
          write(0,*)' ERROR: ',trim(msg),' rc = ', rc1
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
      else
          rc  = esmf_success
          if(iprint) print*, 'pass ',msg
      end if
      return
      end subroutine err_msg

      end module module_error_msg
