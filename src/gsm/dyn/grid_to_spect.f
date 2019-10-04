      subroutine grid_to_spect
     &    (anl_gr_a_1,anl_gr_a_2,
     &     trie_ls,trio_ls,
     &     ls_node,ls_nodes,max_ls_nodes,
     &     lats_nodes_a,global_lats_a,lonsperlat,
     &     epse,epso,plnew_a,plnow_a)
!!
!! hmhj - this routine do spectral to grid transform in model partial reduced grid
! program log:
! 2011 02 20 : henry juang, update code for ndsl and mass_dp
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
!!!!  integer, parameter :: lota = 3*levs+1*levh+1 
!
      real(kind=kind_evod) trie_ls(len_trie_ls,2,lotls)
      real(kind=kind_evod) trio_ls(len_trio_ls,2,lotls)
!!
      real(kind=kind_evod) anl_gr_a_1(lonfx*lota,lats_dim_a)
      real(kind=kind_evod) anl_gr_a_2(lonfx*lota,lats_dim_a)
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
!
!     real(kind=kind_evod), parameter :: cons_0=0.0,   cons_24=24.0
!    &,                                  cons_99=99.0, cons_1p0d9=1.0E9
!    &,                                  qmin=1.0e-10
!
!     real(kind=kind_evod) ga2, tem
!
      INCLUDE 'function2'
!
!--------------------------------------------------------------------
!
      if( ndslfv ) then
        lotx    = 4*levs+1
      else
        lotx    = 4*levs+levh+1
      endif
!
!--------------------------------------------------------------------
      do lan=1,lats_node_a
        lon_dim = lon_dims_a(lan)
        lat = global_lats_a(ipt_lats_node_a-1+lan)
        lons_lat = lonsperlat(lat)
        rcs2 = rcs2_a(min(lat,latg-lat+1))
        do k=1,levs
          do i=1,lons_lat
            anl_gr_a_2(i+(kau+k-2)*lon_dim,lan) = 
     &      anl_gr_a_2(i+(kau+k-2)*lon_dim,lan) * rcs2
            anl_gr_a_2(i+(kav+k-2)*lon_dim,lan) = 
     &      anl_gr_a_2(i+(kav+k-2)*lon_dim,lan) * rcs2
          enddo
        enddo
      enddo
!
      do lan=1,lats_node_a
!
         lon_dim = lon_dims_a(lan)
!
         lat = global_lats_a(ipt_lats_node_a-1+lan)
         lons_lat = lonsperlat(lat)

         call grid2four_thread(anl_gr_a_2(1,lan),anl_gr_a_1(1,lan),
     &                  lon_dim,lons_lat,lonfx,lotx)
!
      enddo
!
      dimg=0
      call four2fln(lats_dim_a,lota,lotx,anl_gr_a_1,
     x              ls_nodes,max_ls_nodes,
     x              lats_nodes_a,global_lats_a,lon_dims_a,
     x              lats_node_a,ipt_lats_node_a,dimg,
     x              lat1s_a,lonfx,latg,latg2,
     x              trie_ls(1,1,P_w), trio_ls(1,1,P_w),
     x              plnew_a, plnow_a,
     x              ls_node,2*levs)
!
!
!$OMP parallel do shared(trie_ls,trio_ls)
!$OMP+shared(kau,kav,kat,epse,epso,ls_node)
!$OMP+private(k)
      do k=1,levs
         call uveodz(trie_ls(1,1,P_w  +k-1), trio_ls(1,1,P_x  +k-1),
     x               trie_ls(1,1,P_uln+k-1), trio_ls(1,1,P_vln+k-1),
     x               epse,epso,ls_node)
!
         call uvoedz(trio_ls(1,1,P_w  +k-1), trie_ls(1,1,P_x  +k-1),
     x               trio_ls(1,1,P_uln+k-1), trie_ls(1,1,P_vln+k-1),
     x               epse,epso,ls_node)
      enddo
!
!   move uln back to x
!   move vln back to w
!
      do k=1,levs
         do i=1,len_trie_ls
            trie_ls(i,1,P_x +k-1)= trie_ls(i,1,P_uln +k-1)
            trie_ls(i,2,P_x +k-1)= trie_ls(i,2,P_uln +k-1)
            trie_ls(i,1,P_w +k-1)= trie_ls(i,1,P_vln +k-1)
            trie_ls(i,2,P_w +k-1)= trie_ls(i,2,P_vln +k-1)
         enddo
         do i=1,len_trio_ls
            trio_ls(i,1,P_x +k-1)= trio_ls(i,1,P_uln +k-1)
            trio_ls(i,2,P_x +k-1)= trio_ls(i,2,P_uln +k-1)
            trio_ls(i,1,P_w +k-1)= trio_ls(i,1,P_vln +k-1)
            trio_ls(i,2,P_w +k-1)= trio_ls(i,2,P_vln +k-1)
         enddo
      enddo
!
!     print *,' exit grid_to_spect '
!!
      return
      end
