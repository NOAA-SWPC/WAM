!#include "../../../ESMFVersionDefine.h"

!
! !description: gfs dynamics gridded component error messages
!
! !revision history:
!
!  january 2007     hann-ming henry juang
!  May     2011     Weiyu yang, modified for using the ESMF 5.2.0r_beta_snapshot_07.
!
!
! !interface:
!
      module gfs_dynamics_err_msg_mod

!
!!uses:
!
      USE ESMF

      implicit none
      logical,parameter::lprint=.false.
      contains

      subroutine gfs_dynamics_err_msg_var(rc1,msg,var,rcfinal)
!
      integer, intent(inout)        :: rc1
      integer, intent(out)          :: rcfinal
      character (len=*), intent(in) :: msg, var
      if(esmf_logfounderror(rc1, msg=msg)) then
          rcfinal = esmf_failure
          print*, 'error happened in dynamics for ',msg,' ',var,' rc = ', rc1
          rc1     = esmf_success
      else
          if(lprint) print*, 'pass in dynamics for ',msg,' ',var
      end if
      end subroutine gfs_dynamics_err_msg_var

      subroutine gfs_dynamics_err_msg(rc1,msg,rc)
      integer, intent(inout)        :: rc1
      integer, intent(out)          :: rc
      character (len=*), intent(in) :: msg
      if(esmf_logfounderror(rc1, msg=msg)) then
          rc  = esmf_failure
          print*, 'error happened in dynamics for ',msg, ' rc = ', rc1
          rc1 = esmf_success
      else
          if(lprint) print*, 'pass in dynamics for ',msg
      end if
      end subroutine gfs_dynamics_err_msg

      subroutine gfs_dynamics_err_msg_final(rcfinal,msg,rc)
      integer, intent(inout)        :: rcfinal, rc
      character (len=*), intent(in) :: msg
      if(rcfinal == esmf_success) then
         if(lprint) print*, "final pass in dynamics for ",msg
      else
         print*, "final fail in dynamics for ",msg
      end if
      rc = rcfinal
      end subroutine gfs_dynamics_err_msg_final

      end module gfs_dynamics_err_msg_mod
