       subroutine spect_to_grid_iau(
     &     trie_ls,trio_ls,datag,
     &     ls_node,ls_nodes,max_ls_nodes,
     &     lats_nodes_a,global_lats_a,lonsperlar,
     &     plnev_a,plnod_a,nlevs,lon1,lon2)

      use gfs_dyn_resol_def
      use gfs_dyn_layout1
      use gfs_dyn_vert_def
      use gfs_dyn_coordinate_def
      use gfs_dyn_mpi_def
      implicit none
      real(kind=kind_evod), intent(in) :: trie_ls(len_trie_ls,2,nlevs)
      real(kind=kind_evod), intent(in) :: trio_ls(len_trio_ls,2,nlevs)
      real(kind=kind_phys), intent(out) :: datag(lon2,lats_node_a,nlevs)
      integer, intent(in) :: ls_node(ls_dim,3),ls_nodes(ls_dim,nodes),
     &  nlevs,max_ls_nodes(nodes),lats_nodes_a(nodes),
     &  global_lats_a(latg),lonsperlar(latg),lon1,lon2
      real(kind=kind_evod),intent(in) :: plnev_a(len_trie_ls,latg2),
     & plnod_a(len_trio_ls,latg2)
! local vars
      real(kind=kind_evod) for_gr_a_1(lon1,nlevs,lats_dim_a)
      real(kind=kind_evod) for_gr_a_2(lon2,nlevs,lats_dim_a)
      integer              i,j,k
      integer              l,lan,lat
      integer              lons_lat

      call sumfln_slg_gg(trie_ls,
     &             trio_ls,
     &             lat1s_a,
     &             plnev_a,plnod_a,
     &             nlevs,ls_node,latg2,
     &             lats_dim_a,nlevs,for_gr_a_1,
     &             ls_nodes,max_ls_nodes,
     &             lats_nodes_a,global_lats_a,
     &             lats_node_a,ipt_lats_node_a,lon1,
     &             lonsperlar,lon1,latg,0)

      do lan=1,lats_node_a
         lat = global_lats_a(ipt_lats_node_a-1+lan)
         lons_lat = lonsperlar(lat)
         CALL FOUR_TO_GRID(for_gr_a_1(1,1,lan),for_gr_a_2(1,1,lan),
     &                          lon1,lon2,lons_lat,nlevs)
      enddo  

      datag = 0.
      do lan=1,lats_node_a
        lat      = global_lats_a(ipt_lats_node_a-1+lan)
        lons_lat = lonsperlar(lat)
        do k=1,nlevs
          do i=1,lons_lat
            datag(i,lan,k) = for_gr_a_2(i,k,lan)
          enddo
        enddo
      enddo
      end subroutine spect_to_grid_iau

