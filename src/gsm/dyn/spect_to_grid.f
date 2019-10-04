      subroutine spect_to_grid
     x    (trie_ls,trio_ls,
     &     syn_gr_a_1,syn_gr_a_2,
     x     ls_node,ls_nodes,max_ls_nodes,
     x     lats_nodes_a,global_lats_a,lonsperlat,
     x     epse,epso,epsedn,epsodn,
     x     snnp1ev,snnp1od,plnev_a,plnod_a)
!
! H.-M. H. Juang: get from gloopa to have only spect to grid
!                 syn_gr_a_1 and syn_gr_a_2 will be used in
!                 the next spect_to_gridxy for saving time
! program log
! 2011 02 20 : henry jaung, updated code for mass_dp and ndsl advection
!
      use gfs_dyn_resol_def
      use gfs_dyn_layout1
      use gfs_dyn_gg_def
      use gfs_dyn_vert_def
      use gfs_dyn_date_def
      use namelist_dynamics_def
      use gfs_dyn_mpi_def
      implicit none
!
      real(kind=kind_evod) trie_ls(len_trie_ls,2,lotls)
      real(kind=kind_evod) trio_ls(len_trio_ls,2,lotls)
!
      integer              ls_node(ls_dim,3)
!
      integer              ls_nodes(ls_dim,nodes)
!
      integer              max_ls_nodes(nodes)
      integer              lats_nodes_a(nodes)
!
      integer              global_lats_a(latg)
      integer                 lonsperlat(latg)
      integer dimg
!
      real(kind=kind_evod)    epse(len_trie_ls)
      real(kind=kind_evod)    epso(len_trio_ls)
      real(kind=kind_evod)  epsedn(len_trie_ls)
      real(kind=kind_evod)  epsodn(len_trio_ls)
!
      real(kind=kind_evod) snnp1ev(len_trie_ls)
      real(kind=kind_evod) snnp1od(len_trio_ls)
!
      real(kind=kind_evod)   plnev_a(len_trie_ls,latg2)
      real(kind=kind_evod)   plnod_a(len_trio_ls,latg2)
!
      real(kind=kind_evod) syn_gr_a_1(lonfx*lots,lats_dim_a)
      real(kind=kind_evod) syn_gr_a_2(lonfx*lots,lats_dim_a)
!
      integer              i,j,k,l
      integer              lan,lat
      integer              lon_dim,lons_lat,n,node
!
! xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
      real(kind=kind_evod) cons0,cons2     !constant
!
      cons0 = 0.d0     !constant
      cons2 = 2.d0     !constant
!
! ................................................................
!

      call f_hpmstart(1,"ga delnpe")
      call delnpe(trie_ls(1,1,P_q   ),
     x            trio_ls(1,1,P_dphi),
     x            trie_ls(1,1,P_dlam),
     x            epse,epso,ls_node)
      call delnpe(trie_ls(1,1,P_gz  ),
     x            trio_ls(1,1,P_zsphi),
     x            trie_ls(1,1,P_zslam),
     x            epse,epso,ls_node)
      call f_hpmstop(1)
!
      call f_hpmstart(2,"ga delnpo")
      call delnpo(trio_ls(1,1,P_q   ),
     x            trie_ls(1,1,P_dphi),
     x            trio_ls(1,1,P_dlam),
     x            epse,epso,ls_node)
      call delnpo(trio_ls(1,1,P_gz  ),
     x            trie_ls(1,1,P_zsphi),
     x            trio_ls(1,1,P_zslam),
     x            epse,epso,ls_node)
      call f_hpmstop(2)
!
      call f_hpmstart(3,"ga dezouv dozeuv")
!
!$omp parallel do shared(trie_ls,trio_ls)
!$omp+shared(epsedn,epsodn,snnp1ev,snnp1od,ls_node)
!$omp+private(k)
      do k=1,levs
         call dezouv(trie_ls(1,1,P_di +k-1), trio_ls(1,1,P_ze +k-1),
     x               trie_ls(1,1,P_uln+k-1), trio_ls(1,1,P_vln+k-1),
     x               epsedn,epsodn,snnp1ev,snnp1od,ls_node)
!
         call dozeuv(trio_ls(1,1,P_di +k-1), trie_ls(1,1,P_ze +k-1),
     x               trio_ls(1,1,P_uln+k-1), trie_ls(1,1,P_vln+k-1),
     x               epsedn,epsodn,snnp1ev,snnp1od,ls_node)
      enddo
!
      call f_hpmstop(3)
!
      CALL countperf(0,3,0.)
      CALL synctime()
      CALL countperf(1,3,0.)
!!
      dimg=0
      CALL countperf(0,1,0.)
      call f_hpmstart(8,"ga sumflna")
!
      call sumflna(trie_ls(1,1,P_ze),
     x            trio_ls(1,1,P_ze),
     x            lat1s_a,
     x            plnev_a,plnod_a,
     x            lots,ls_node,latg2,
     x            lats_dim_a,lots,
     x            syn_gr_a_1,
     x            ls_nodes,max_ls_nodes,
     x            lats_nodes_a,global_lats_a,
     x            lats_node_a,ipt_lats_node_a,lon_dims_a,dimg,
     x            lonsperlat,lonfx,latg)
!
      call f_hpmstop(8)
      CALL countperf(1,1,0.)
!
      call f_hpmstart(10,"ga lat_loop")
!11111111111111111111111111111111111111111111111111111111111111111111
      do lan=1,lats_node_a  
 
         lon_dim = lon_dims_a(lan)
         lat = global_lats_a(ipt_lats_node_a-1+lan)
         lons_lat = lonsperlat(lat)

         CALL countperf(0,6,0.)
         CALL FOUR2GRID_thread(syn_gr_a_1(1,lan),syn_gr_a_2(1,lan),
     &                  lon_dim,lons_lat,lonfx,lots,lan,me)
         CALL countperf(1,6,0.)

      enddo 
!11111111111111111111111111111111111111111111111111111111111111111111
!
      return
      end
