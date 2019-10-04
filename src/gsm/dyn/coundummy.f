      subroutine init_countperf(lat)
c    intro of implicit none > >
	implicit none      
	integer lat
      return
      end
      subroutine countperf(flag,ic,nop)
c    intro of implicit none > >
        implicit none
      integer flag,ic
      real nop
      return
      end
      subroutine end_countperf()
c    intro of implicit none > >
        implicit none
      return
      end
      subroutine write_countperf(npes,me,hours)
c    intro of implicit none > >
        implicit none
      integer npes,me
      real hours
      return
      end
      subroutine synchro
c    intro of implicit none > >
        implicit none
      return
      end
      FUNCTION timer()
c    intro of implicit none > >
        implicit none
      integer timer
      timer=0.
      return
      end
      subroutine synctime()
c    intro of implicit none > >
        implicit none
      return
      end
