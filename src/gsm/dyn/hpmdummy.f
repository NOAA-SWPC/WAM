      subroutine f_hpminit(me,desc)
	implicit none
	integer me
c	^^^^^^^^^^^^^
      	character*20 desc
      return
      end
      subroutine f_hpmterminate(me)
        implicit none
        integer  me
c       ^^^^^^^^^^^^^
      return
      end
      subroutine f_hpmstart(n,desc)
        implicit none
        integer  n
c       ^^^^^^^^^^^^^
     	character*20 desc
      return
      end
      subroutine f_hpmstop(n)
!        implicit none
        integer  n
c       ^^^^^^^^^^^^^
     return
      end
