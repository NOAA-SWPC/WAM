      module gfs_dyn_gg_def
      use gfs_dyn_machine

      implicit none
      
      real(kind=kind_dbl_prec), allocatable, dimension(:) ::  colrad_a,
     &                      wgt_a, wgtcs_a, rcs2_a, sinlat_a, coslat_a
!
      integer ,allocatable, dimension(:) :: lats_nodes_h,global_lats_h
      end module gfs_dyn_gg_def
