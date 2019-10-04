      subroutine init_countperf(lat)
c	vvv
	implicit none
	integer lat
      return
      end
      subroutine countperf(flag,ic,nop)
c	vvv
	implicit none
      integer flag,ic
      real nop
      return
      end
      subroutine end_countperf()
	implicit none
      return
      end
      subroutine write_countperf(npes,me,hours)
	implicit none
      integer npes,me
      real hours
      return
      end
      subroutine synchro
	implicit none
      return
      end
      FUNCTION timer()
	implicit none
      integer timer
      timer=0.
      return
      end
      subroutine synctime()
      	implicit none
	return
      end
