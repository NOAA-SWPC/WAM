      subroutine grid_to_spect_rqt
     &    (rqt_gr_a_1,rqt_gr_a_2,
     &     trie_ls_rqt,trio_ls_rqt,lotx,
     &     ls_node,ls_nodes,max_ls_nodes,
     &     lats_nodes_a,global_lats_a,lonsperlat,
     &     epse,epso,plnew_a,plnow_a)
!!
!! hmhj - this routine do grid to spectral transform in model partial reduced grid
!!        for moisture only
! program log:
! 2013 06 20 : henry juang, moisture grid to spectral transform for ndsl version
!!
      use gfs_dyn_resol_def
      use gfs_dyn_layout1
      use gfs_dyn_gg_def
      use gfs_dyn_vert_def
      use gfs_dyn_date_def
      use namelist_dynamics_def
      use gfs_dyn_coordinate_def 
      use gfs_dyn_tracer_const
      use gfs_dyn_physcons,  grav  => con_g
      implicit none
!!
!
      real(kind=kind_evod) trie_ls_rqt(len_trie_ls,2,lotx)
      real(kind=kind_evod) trio_ls_rqt(len_trio_ls,2,lotx)
!!
      real(kind=kind_evod) rqt_gr_a_1(lonfx*lotx,lats_dim_a)
      real(kind=kind_evod) rqt_gr_a_2(lonfx*lotx,lats_dim_a)
!
      integer              ls_node(ls_dim,3)
      integer              ls_nodes(ls_dim,nodes)
      integer              max_ls_nodes(nodes)
      integer              lats_nodes_a(nodes)
      integer              global_lats_a(latg)
      integer                 lonsperlat(latg)
      integer dimg
!
      real(kind=kind_evod)  epse(len_trie_ls)
      real(kind=kind_evod)  epso(len_trio_ls)
!
      real(kind=kind_evod)  rcs2
!
      real(kind=kind_evod)  plnew_a(len_trie_ls,latg2)
      real(kind=kind_evod)  plnow_a(len_trio_ls,latg2)
!
      integer              i,j,k, nn, nnl
      integer              l,lan,lat,lotx
      integer              lon_dim,lons_lat
!
      integer              locl,n
      integer              indev
      integer              indod
      integer              indev1,indev2
      integer              indod1,indod2
      INTEGER              INDLSEV,JBASEV
      INTEGER              INDLSOD,JBASOD
!

!timers______________________________________________________---
      real*8 rtc ,timer1,timer2
!timers______________________________________________________---
!
      INCLUDE 'function2'
!
!--------------------------------------------------------------------
!
      do lan=1,lats_node_a
!
         lon_dim = lon_dims_a(lan)
!
         lat = global_lats_a(ipt_lats_node_a-1+lan)
         lons_lat = lonsperlat(lat)

         call grid2four_thread(rqt_gr_a_2(1,lan),rqt_gr_a_1(1,lan),
     &                  lon_dim,lons_lat,lonfx,lotx)
!
      enddo
!
      dimg=0
      call four2fln(lats_dim_a,lotx,lotx,rqt_gr_a_1,
     x              ls_nodes,max_ls_nodes,
     x              lats_nodes_a,global_lats_a,lon_dims_a,
     x              lats_node_a,ipt_lats_node_a,dimg,
     x              lat1s_a,lonfx,latg,latg2,
     x              trie_ls_rqt, trio_ls_rqt,
     x              plnew_a, plnow_a,
     x              ls_node,0)
!
!
!     print *,' exit grid_to_spect '
!!
      return
      end
