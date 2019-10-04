      subroutine spectpz_to_gridxy
     x    (trie_ls,trio_ls,
     &     pyn_gr_a_1,pyn_gr_a_2,
     x     ls_node,ls_nodes,max_ls_nodes,
     x     lats_nodes_a,global_lats_a,lonsperlat,
     x     epse,epso,epsedn,epsodn,
     x     snnp1ev,snnp1od,plnev_a,plnod_a)
!
! H.-M. H. Juang:  compute dpx dpy phix and phiy
! program log
! 2011 02 20 : henry juang, initially implemented into nems for options of
!              mass_dp and ndsl advection
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
      real(kind=kind_evod) pyn_gr_a_1(lonfx*lotp,lats_dim_a)
      real(kind=kind_evod) pyn_gr_a_2(lonfx*lotp,lats_dim_a)
!
      integer              i,j,k,kk,l
      integer              lan,lat,lotx
      integer              lon_dim,lons_lat,n,node
!
! xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
      real(kind=kind_evod) cons0,cons2     !constant
!
      cons0 = 0.d0     !constant
      cons2 = 2.d0     !constant
      lotx = levs * 4
!     print *,' p_dp p_rq p_q ',p_dp,p_rq,p_q
!     print *,' p_p_dplam p_dpphi ',p_dplam,p_dpphi
!     print *,' p_zz p_zzlam p_zzphi ',p_zz,p_zzlam,p_zzphi
!
! ................................................................
!
!$omp parallel do shared(trie_ls,trio_ls)
!$omp+shared(epse,epso,ls_node)
!$omp+private(k)
      do k=1,levs

      call delnpe(trie_ls(1,1,P_dp   +k-1),
     x            trio_ls(1,1,P_dpphi+k-1),
     x            trie_ls(1,1,P_dplam+k-1),
     x            epse,epso,ls_node)
      call delnpe(trie_ls(1,1,P_zz   +k-1),
     x            trio_ls(1,1,P_zzphi+k-1),
     x            trie_ls(1,1,P_zzlam+k-1),
     x            epse,epso,ls_node)
!
      call delnpo(trio_ls(1,1,P_dp   +k-1),
     x            trie_ls(1,1,P_dpphi+k-1),
     x            trio_ls(1,1,P_dplam+k-1),
     x            epse,epso,ls_node)
      call delnpo(trio_ls(1,1,P_zz   +k-1),
     x            trie_ls(1,1,P_zzphi+k-1),
     x            trio_ls(1,1,P_zzlam+k-1),
     x            epse,epso,ls_node)

      enddo

! hmhj debug print
!     do k=1,levs
!     print *,' k= ',k
!     print *,' trie dp  ',(trie_ls(i,1,p_dp+k-1),i=1,3)
!     print *,' trio dp  ',(trio_ls(i,1,p_dp+k-1),i=1,3)
!     print *,' trie dpx ',(trie_ls(i,1,p_dplam+k-1),i=1,3)
!     print *,' trio dpx ',(trio_ls(i,1,p_dplam+k-1),i=1,3)
!     print *,' trie dpy ',(trie_ls(i,1,p_dpphi+k-1),i=1,3)
!     print *,' trio dpy ',(trio_ls(i,1,p_dpphi+k-1),i=1,3)
!     enddo
!
      CALL countperf(0,3,0.)
      CALL synctime()
      CALL countperf(1,3,0.)
!!
      lotx = levs*4
      dimg=0
      CALL countperf(0,1,0.)
      call f_hpmstart(8,"ga sumflna")
!
      call sumflna(trie_ls(1,1,P_dpphi),
     x            trio_ls(1,1,P_dpphi),
     x            lat1s_a,
     x            plnev_a,plnod_a,
     x            lotx,ls_node,latg2,
     x            lats_dim_a,lotp,
     x            pyn_gr_a_1,
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
         CALL FOUR2GRID_thread(pyn_gr_a_1(1,lan),pyn_gr_a_2(1,lan),
     &                  lon_dim,lons_lat,lonfx,lotx,lan,me)
         CALL countperf(1,6,0.)

! hmhj debug print
!        do k=1,levs
!        kk=k+levs
!     print *,' k = ',k
!     print *,' pyndpx ',(pyn_gr_a_2(i+(k-1)*lon_dim,lan),i=1,3)
!     print *,' pyndpy ',(pyn_gr_a_2(i+(kk-1)*lon_dim,lan),i=1,3)
!        enddo
      
      enddo 
!
      return
      end
