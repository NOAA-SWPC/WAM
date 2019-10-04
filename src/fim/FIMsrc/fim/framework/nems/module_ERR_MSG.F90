!JR Rev. 11555 of NEMS repo.
!
! !description: error messages
!
! !revision history:
!
!  january 2007 	hann-ming henry juang
!
!
! !interface:
!
      module module_err_msg

!
!!uses:
!
      use esmf_mod 

      implicit none

      private
      public :: err_msg,message_check,set_iprint

!TODO:  Talk with NCEP about new "iprint:" field in config file and 
!TODO:  setting module_ERR_MSG::iprint at run-time from MAIN_NEMS.F90.  
!TODO:  Note that GFS has its own version of this module that differs 
!TODO:  slightly (as of nems r12470).  

      logical :: iprint = .false.
      character(esmf_maxstr) :: message_check

      contains

      ! allow verbosity to be adjusted at run-time
      subroutine set_iprint(val)
        logical, intent(in) :: val
        iprint = val
      end subroutine set_iprint

      subroutine err_msg_int(rc1,msg,val,rcfinal)
!
      integer, intent(inout)        :: rc1
      integer, intent(out)          :: rcfinal
      character (len=*), intent(in) :: msg
      integer,           intent(in) :: val
      if(esmf_logmsgfounderror(rc1, msg)) then
          rcfinal = esmf_failure
          print*, 'error happened for ',msg,' ',val,' rc = ', rc1
          write(0,*)' ERROR: ',msg,' ',val,' rc = ', rc1
          rc1     = esmf_success
      else
          if(iprint) print*, 'pass ',msg,' ',val
      end if
      return
      end subroutine err_msg_int

      subroutine err_msg_val(rc1,msg,val,rcfinal)
!
      integer, intent(inout)        :: rc1
      integer, intent(out)          :: rcfinal
      character (len=*), intent(in) :: msg
      real,              intent(in) :: val
      if(esmf_logmsgfounderror(rc1, msg)) then
          rcfinal = esmf_failure
          print*, 'error happened for ',msg,' ',val,' rc = ', rc1
          write(0,*)' ERROR: ',msg,' ',val,' rc = ', rc1
          rc1     = esmf_success
      else
          if(iprint) print*, 'pass ',msg,' ',val
      end if
      return
      end subroutine err_msg_val

      subroutine err_msg_var(rc1,msg,chr,rcfinal)
!
      integer, intent(inout)        :: rc1
      integer, intent(out)          :: rcfinal
      character (len=*), intent(in) :: msg
      character (len=*), intent(in) :: chr
      if(esmf_logmsgfounderror(rc1, msg)) then
          rcfinal = esmf_failure
          print*, 'error happened for ',msg,' ',chr,' rc = ', rc1
          write(0,*)' ERROR: ',msg,' ',chr,' rc = ', rc1
          rc1     = esmf_success
      else
          if(iprint) print*, 'pass ',msg,' ',chr
      end if
      return
      end subroutine err_msg_var

      subroutine err_msg(rc1,msg,rc)
      ! use-association to make prints more informative
      use module_core_setup,only:core_setup_done,iam_fim_task, &
                                 iam_write_task,world_rank
      integer, intent(inout)        :: rc1
      integer, intent(out)          :: rc
      character (len=*), intent(in) :: msg
      character (len=10) :: task_type
      if      (iam_fim_task)   then
        task_type="compute"
      else if (iam_write_task) then
        task_type="write"
      else
        task_type="do-nothing"
      endif
      if(esmf_logmsgfounderror(rc1, msg)) then
          rc  = esmf_failure
          if (core_setup_done) then
            write(6,'(3a,i0,3a,i0)')                         &
              'ERROR [',trim(task_type),' task:',world_rank, &
              '] ',trim(msg),' rc = ',rc1
            write(0,'(3a,i0,3a,i0)')                         &
              'ERROR [',trim(task_type),' task:',world_rank, &
              '] ',trim(msg),' rc = ',rc1
          else
            write(6,'(3a,i0)') 'ERROR ',trim(msg),' rc = ',rc1
            write(0,'(3a,i0)') 'ERROR ',trim(msg),' rc = ',rc1
          endif
          rc1 = esmf_success
      else
          rc  = esmf_success
          if (iprint) then
            if (core_setup_done) then
              write(6,'(3a,i0,2a)')                           &
                'pass [',trim(task_type),' task:',world_rank, &
                 '] ',trim(msg)
            else
              write(6,'(2a)') 'pass ',trim(msg)
            endif
          endif
      end if
      return
      end subroutine err_msg

      subroutine err_msg_final(rcfinal,msg,rc)
      integer, intent(inout)        :: rcfinal
      integer, intent(inout)        :: rc
      character (len=*), intent(in) :: msg
      if(rcfinal == esmf_success) then
          if(iprint) print*, "final pass: ",msg
      else
          print*, "final fail: ",msg
          write(0,*)' FINAL ERROR: ',msg
      end if
!     if(present(rc)) then
          rc = rcfinal
!     end if
      return
      end subroutine err_msg_final

      end module module_err_msg
