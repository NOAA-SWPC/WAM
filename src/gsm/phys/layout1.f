      module layout1
      implicit none
      
!
      integer           nodes, nodes_comp,nodes_io,
     &                  me,
     &                  ls_dim,
     &                  ls_max_node,
     &                  lats_dim_r,
     &                  lats_dim_ext,
     &                  lats_node_r,
     &                  lats_node_r_max,
     &                  lats_node_ext,
     &                  ipt_lats_node_r,
     &                  ipt_lats_node_ext,
!    &			lonf, latg,
     &                  len_trio_ls, len_trie_ls,
     &                  me_l_0
!
      integer           idrt                                            !jw:for flx file outfile
      logical           comp_task                                       ! moorthi

      integer ,allocatable :: lon_dims_r(:),lon_dims_ext(:)
      end module layout1
