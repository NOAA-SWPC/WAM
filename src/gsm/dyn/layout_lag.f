      module layout_lag
      use gfs_dyn_MACHINE, only : kind_evod
      implicit none
      save
cc
      integer lats_dim_h,
     x        lats_node_h,
     x        lats_node_h_max,
     x        ipt_lats_node_h,
     x        lon_dim_h
cc
      real(kind=kind_evod), allocatable :: cosphi(:),sinphi(:)
cc
      INTEGER ,ALLOCATABLE :: lat1s_h(:)
!mjr  INTEGER ,ALLOCATABLE :: lon_dims_h_coef(:)
!mjr  INTEGER ,ALLOCATABLE :: lon_dims_h_grid(:)
      end module layout_lag
