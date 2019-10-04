      integer function num_parthds()
      use omp_lib
c    intro of implicit none > >
	implicit none
	
!!$OMP PARALLEL
!     num_parthds = omp_get_num_threads()
      num_parthds = 8
!     num_parthds = 6
!     num_parthds = 4
!!$OMP END PARALLEL
      return
      end

 
