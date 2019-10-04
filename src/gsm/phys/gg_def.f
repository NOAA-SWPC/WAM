      module gg_def
      use machine, ONLY: KIND_EVOD

      implicit none
!-----------------------------------------------------------------
!  revision
!
!  Dec 2014  J. Wang  add slat_r dlat for local Gaussian latitude
!
      REAL(KIND=KIND_EVOD) ,ALLOCATABLE ::  colrad_r(:),wgt_r(:),
     & wgtcs_r(:),rcs2_r(:),sinlat_r(:),coslat_r(:),slat_r(:),dlat_r(:)
      end module gg_def
