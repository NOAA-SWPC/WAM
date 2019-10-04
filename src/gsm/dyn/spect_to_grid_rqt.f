      subroutine spect_to_grid_rqt
     &    (trie_rq,trio_rq,rq_gr_a_1,rq_gr_a_2,lotx,
     &     ls_node,ls_nodes,max_ls_nodes,
     &     lats_nodes_a,global_lats_a,lonsperlat,
     &     plnev_a,plnod_a)
!!
!  routine to do spectral transform initially for moisture only
!  initally coded by hann-ming henry juang
!
! program log:
! 2013 06 20 :    Henry Juang, do spectral to grid transform on moisture for NDSL
!!
      use gfs_dyn_resol_def
      use gfs_dyn_layout1
      use gfs_dyn_gg_def
      use gfs_dyn_vert_def
      use gfs_dyn_date_def
      use namelist_dynamics_def
      use gfs_dyn_coordinate_def 
      implicit none
!!
      real(kind=kind_evod) trie_rq(len_trie_ls,2,lotx)
      real(kind=kind_evod) trio_rq(len_trio_ls,2,lotx)
!
      real(kind=kind_evod) 
     &              rq_gr_a_1(lonfx*lotx,lats_dim_a),
     &              rq_gr_a_2(lonfx*lotx,lats_dim_a)
cc
      integer              ls_node, ls_nodes(ls_dim,nodes)
      integer              max_ls_nodes(nodes)
      integer              lats_nodes_a(nodes)
      integer              global_lats_a(latg)
      integer                 lonsperlat(latg)
      integer dimg
cc
      real(kind=kind_evod)   plnev_a(len_trie_ls,latg2)
      real(kind=kind_evod)   plnod_a(len_trio_ls,latg2)
cc
      integer              i,j,k
      integer              l,lan,lat,lotx
      integer              lon_dim,lons_lat
!
cc
cc
      real(kind=kind_evod) cons0,cons2     !constant
cc
      real(kind=kind_evod), parameter :: cons_0=0.0,   cons_24=24.0
     &,                                  cons_99=99.0, cons_1p0d9=1.0E9
!
cc
cc--------------------------------------------------------------------
cc
cc
      dimg=0
cc
      call sumflna(trie_rq,trio_rq,
     x            lat1s_a,
     x            plnev_a,plnod_a,
     x            lotx,ls_node,latg2,
     x            lats_dim_a,lotx, rq_gr_a_1,
     x            ls_nodes,max_ls_nodes,
     x            lats_nodes_a,global_lats_a,
     x            lats_node_a,ipt_lats_node_a,lon_dims_a,dimg,
     x            lonsperlat,lonfx,latg)
c
      do lan=1,lats_node_a
cc
         lat = global_lats_a(ipt_lats_node_a-1+lan)
         lon_dim = lon_dims_a(lan)
         lons_lat = lonsperlat(lat)
cc
         call four2grid_thread( rq_gr_a_1(1,lan), rq_gr_a_2(1,lan),
     &                  lon_dim,lons_lat,lonfx,lotx,lan,me)
!
      enddo   !lan
cc
!     print *,' exit spect_to_grid_rqt '
!!
      return
      end
