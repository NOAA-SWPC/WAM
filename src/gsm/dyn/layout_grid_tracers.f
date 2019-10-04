      module layout_grid_tracers
      use gfs_dyn_MACHINE, only : kind_evod
      implicit none
      save
!
      integer xhalo, yhalo
!
      real(kind=kind_evod) ,allocatable :: rgt_h(:,:,:,:)
!
      real(kind=kind_evod) ,allocatable :: rgt_a(:,:,:,:)
!
      end module layout_grid_tracers
