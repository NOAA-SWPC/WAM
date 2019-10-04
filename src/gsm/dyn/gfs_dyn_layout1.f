      module gfs_dyn_layout1
      implicit none
!
! program log:
! 20110220     henry jaung  add more indexes for mass_dp and ndslfv options
!

      integer           nodes, nodes_comp,nodes_io,
     &                  me,lon_dim_a,
     &                  ls_dim,
     &                  ls_max_node,
     &                  lats_dim_a,
     &                  lats_node_a,
     &                  lats_node_a_max,
     &                  ipt_lats_node_a,
     &                  len_trie_ls,
     &                  len_trio_ls,
     &                  len_trie_ls_max,
     &                  len_trio_ls_max,
     &                  me_l_0,

     &                  lats_dim_ext,
     &                  lats_node_ext,
     &                  ipt_lats_node_ext
!
      INTEGER ,ALLOCATABLE :: lat1s_a(:), lon_dims_a(:),lon_dims_ext(:)

!hmhj ndslfv
      integer   lonfull,lonhalf,lonpart,lonlenmax,mylonlen
     &,         latfull,lathalf,latpart,latlenmax,mylatlen
     &,         ndslhvar,ndslvvar

      integer, allocatable :: lonstr(:),lonlen(:),latstr(:),latlen(:)
      real,    allocatable :: gglat(:),gglati(:),gslati(:),ggfact(:,:)
      real,    allocatable :: gglon(:),ggloni(:),cosglat(:)

      end module gfs_dyn_layout1
