      function num_parthds()
!#ifdef _OPENMP
!     use omp_lib
!!$OMP PARALLEL
!     num_parthds=omp_get_num_threads()
!!$OMP END PARALLEL
!#else
      num_parthds=8
!     num_parthds=1
!#endif
      return
      end 

