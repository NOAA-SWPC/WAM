!-----------------------------------------------------------------------
                        module module_alarms
!-----------------------------------------------------------------------
! 
      use ESMF
      implicit none
!
! instantiate Alarm lists
!     integer, parameter, save :: ALARM_NUM = 2
!     type(ESMF_Alarm), save   :: alarm(ALARM_NUM)
      type(ESMF_Alarm), save   :: alarm(2)
!
      end module module_alarms
