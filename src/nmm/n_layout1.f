      module n_layout1
      implicit none
      
cc
      integer           nodes, nodes_comp,nodes_io,
     x                  me,
     x                  ls_dim,
     x                  ls_max_node,
     x                  lats_dim_r,
     x                  lats_dim_ext,
     x                  lats_node_r,
     x                  lats_node_r_max,
     x                  lats_node_ext,
     x                  ipt_lats_node_r,
     x                  ipt_lats_node_ext,
!    x			lonf, latg,
     x                  len_trio_ls, len_trie_ls,
     x                  me_l_0
cc
      integer           idrt                                            !jw:for flx file outfile

      integer ,allocatable :: lon_dims_r(:),lon_dims_ext(:)
      end module n_layout1
