      module bfilt_def
      use gfs_dyn_MACHINE, ONLY: KIND_EVOD
      implicit none
      save
      
      REAL(KIND=KIND_EVOD) ,ALLOCATABLE :: bfilte(:),bfilto(:)
      end module bfilt_def
