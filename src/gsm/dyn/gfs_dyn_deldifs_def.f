      module gfs_dyn_deldifs_def
      use gfs_dyn_machine
      implicit none
      
      REAL(KIND=KIND_EVOD),ALLOCATABLE :: DNE(:),DNO(:),
     . SF(:),RTRD(:),RTHK(:),BKLY(:),CKLY(:)			! hmhj
!hmhj idea Apr 06 2012
!jw idea Jan 2013
      REAL(KIND=KIND_EVOD),ALLOCATABLE :: DNEidea(:),DNOidea(:)
      end module gfs_dyn_deldifs_def
