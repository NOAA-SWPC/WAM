      subroutine grid_to_spect_slg(anl_gr_a_1,anl_gr_a_2
     &,                            trie_ls,trio_ls,lsout
     &,                            ls_node,ls_nodes,max_ls_nodes
     &,                            lats_nodes_a,global_lats_a,lonsperlat
     &,                            epse,epso,plnew_a,plnow_a)
!!
!! hmhj - this routine do spectral to grid transform in model partial reduced grid
! program log:
! 2011 02 20 : henry juang, update code for ndsl and mass_dp
! 2014 08 ?? : Shrinivas Moorthi - update for sela slg - including gridded tracers
!!
      use gfs_dyn_resol_def
      use gfs_dyn_layout1
      use gfs_dyn_gg_def
      use gfs_dyn_vert_def
      use gfs_dyn_date_def
      use namelist_dynamics_def
      use gfs_dyn_coordinate_def 
      use gfs_dyn_tracer_const
      use layout_grid_tracers , only : rgt_a,rgt_h,xhalo,yhalo
      use gfs_dyn_physcons,  grav  => con_g
      implicit none
!!
!
!!!!  integer, parameter :: lota = 3*levs+1*levh+1 
!
      logical lsout
      real(kind=kind_evod) trie_ls(len_trie_ls,2,lotls)
     &,                    trio_ls(len_trio_ls,2,lotls)
     &,                    anl_gr_a_1(lon_dim_a,lota,lats_dim_a)
     &,                    anl_gr_a_2(lonf,lota,lats_dim_a)
!!
!     real(kind=kind_evod) anl_gr_a_1(lonfx*lota,lats_dim_a)
!     real(kind=kind_evod) anl_gr_a_2(lonfx*lota,lats_dim_a)
!
      integer              ls_node(ls_dim,3), ls_nodes(ls_dim,nodes)
      integer, dimension(nodes) :: max_ls_nodes, lats_nodes_a
      integer, dimension(latg)  :: global_lats_a, lonsperlat
!
      real(kind=kind_evod) epse(len_trie_ls), epso(len_trio_ls)
     &,                    plnew_a(len_trie_ls,latg2)
     &,                    plnow_a(len_trio_ls,latg2)
!
      real(kind=kind_evod)  rcs2, rcs
!
      integer  i,j,k, nn, nnl,l,lan,lat,lotx,item,jtem,ktem
     &,        lon_dim,lons_lat,dimg
     &,        locl,n,indev,indod,indev1,indev2,indod1,indod2
     &,        INDLSEV,JBASEV,INDLSOD,JBASOD,ylan
!

!timers______________________________________________________---
!     real*8 rtc ,timer1,timer2
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
!     if( ndslfv ) then
!       lotx    = 4*levs+1
!     else
!       lotx    = 4*levs+levh+1
!     endif
!
!--------------------------------------------------------------------
!       if (gg_tracers) then
!         do i=1,lons_lat
!           bak_gr_r_2(i,kap,lan) = rqtk(i)
!         enddo
!       else
!         do i=1,lons_lat
!           bak_gr_r_2(i,kap,lan) = 0.0
!         enddo
!       endif
!     enddo
!
      kar = 3*levs + 2
      do lan=1,lats_node_a

!     write(150+me,*)' in grid_2_spec anl_2=',anl_gr_a_2(1,kar+levh,lan)
!    &,' lan=',lan,' me=',me

         lon_dim  = lon_dims_a(lan)
         lat      = global_lats_a(ipt_lats_node_a-1+lan)
         lons_lat = lonsperlat(lat)

         call grid_to_four(anl_gr_a_2(1,1,lan),anl_gr_a_1(1,1,lan),
     &                  lon_dim_a-2,lon_dim_a,lons_lat,3*levs+1)

         if (.not. gg_tracers .or. lsout) then

!     write(150+me,*)' grid_to_spec anl2 maxgr=',
!    &maxval(anl_gr_a_2(:,kar+levh-1,lan)),
!    &  ' mingr=',minval(anl_gr_a_2(:,kar+levh-1,lan)),
!    &' lan=',lan,' kar+levh-1=',kar+levh-1

            call grid_to_four(anl_gr_a_2(1,kar,lan),
     &                        anl_gr_a_1(1,kar,lan),
     &                        lon_dim_a-2,lon_dim_a,lons_lat,levh)

!     write(150+me,*)'in grid_tospec anl1=',anl_gr_a_1(1,kar+levh-1,lan)
!    &,' anl2=',anl_gr_a_2(1,kar+levh-1,lan),'lan=',lan

         endif
!
        if (gg_tracers) then
!         if (.not. shuff_lats_r) then
            item = lats_node_a + 1 - lan
            ylan = yhalo + lan
            do n=1,ntrac
!$omp parallel do private(k,jtem,ktem,i)
              do k=1,levs
                jtem = levs + 1 - k
                ktem = kar  - 1 + k + (n-1)*levs
                do i=1,min(lonf,lons_lat)
                  rgt_h(xhalo+i,jtem,ylan,n) = anl_gr_a_2(i,ktem,item)

!$$$            if (kdt .eq. 1) write(888,*) 'rg1_h, = ',
!$$$     &      i,k,lan, rg1_h(xhalo+i,levs+1-k,lats_node_a+1-lan+yhalo)
                enddo
              enddo
            enddo
!         endif ! .not.shuff_lats_r
        endif ! gg_tracers

!       if (gg_tracers .and. shuff_lats_r) then
!!        print*,' gloopb mpi_tracers_b_to_a shuff_lats_r',shuff_lats_r
!!!mr     call mpi_barrier (mc_comp,ierr)
!         call mpi_tracers_b_to_a(anl_gr_r_2(1,1,1),
!     &                           lats_nodes_r,global_lats_r,
!     &                           rgt_h,lats_nodes_a,global_lats_a,kar,0)
!       endif ! gg_tracers .and. shuff_lats_r

      enddo
!
      call four2fln_gg(lats_dim_a,lota,3*levs+1,anl_gr_a_1,
     &                 ls_nodes,max_ls_nodes,
     &                 lats_nodes_a,global_lats_a,lon_dim_a,
     &                 lats_node_a,ipt_lats_node_a,
     &                 lat1s_a,lonf+2,latg,latg2,
     &                 trie_ls(1,1,p_ze), trio_ls(1,1,p_ze),
     &                 plnew_a, plnow_a,
     &                 ls_node,0,
     &                 2*levs+1,3*levs+1)

      if (.not. gg_tracers .or. lsout ) then
         call four2fln_gg(lats_dim_a,lota,levh,anl_gr_a_1,
     &                    ls_nodes,max_ls_nodes,
     &                    lats_nodes_a,global_lats_a,lon_dim_a,
     &                    lats_node_a,ipt_lats_node_a,
     &                    lat1s_a,lonf+2,latg,latg2,
     &                    trie_ls(1,1,p_rq), trio_ls(1,1,p_rq),
     &                    plnew_a, plnow_a,
     &                    ls_node,3*levs+1,
     &                    1,levh)
!     write(150+me,*)' in grid_tospec trie=',trie_ls(1,1,p_rq+levh-1)
!    &,' p_rq+levh-1=',p_rq+levh-1
      endif
      dimg=0
!     call four2fln(lats_dim_a,lota,lotx,anl_gr_a_1,
!    x              ls_nodes,max_ls_nodes,
!    x              lats_nodes_a,global_lats_a,lon_dims_a,
!    x              lats_node_a,ipt_lats_node_a,dimg,
!    x              lat1s_a,lonfx,latg,latg2,
!    x              trie_ls(1,1,P_w), trio_ls(1,1,P_w),
!    x              plnew_a, plnow_a,
!    x              ls_node,2*levs)
!
!
!$OMP parallel do shared(trie_ls,trio_ls)
!$OMP+shared(epse,epso,ls_node)
!$OMP+private(k)
      do k=1,levs
         call uveodz(trie_ls(1,1,p_ze +k-1), trio_ls(1,1,p_di +k-1),
     &               trie_ls(1,1,p_uln+k-1), trio_ls(1,1,p_vln+k-1),
     &               epse,epso,ls_node)
!!
         call uvoedz(trio_ls(1,1,p_ze +k-1), trie_ls(1,1,p_di +k-1),
     &               trio_ls(1,1,p_uln+k-1), trie_ls(1,1,p_vln+k-1),
     &               epse,epso,ls_node)

!        call uveodz(trie_ls(1,1,P_w  +k-1), trio_ls(1,1,P_x  +k-1),
!    x               trie_ls(1,1,P_uln+k-1), trio_ls(1,1,P_vln+k-1),
!    x               epse,epso,ls_node)
!
!        call uvoedz(trio_ls(1,1,P_w  +k-1), trie_ls(1,1,P_x  +k-1),
!    x               trio_ls(1,1,P_uln+k-1), trie_ls(1,1,P_vln+k-1),
!    x               epse,epso,ls_node)
      enddo
!
!   move uln back to x
!   move vln back to w
!
!     do k=1,levs
!        do i=1,len_trie_ls
!           trie_ls(i,1,P_x +k-1) = trie_ls(i,1,P_uln +k-1)
!           trie_ls(i,2,P_x +k-1) = trie_ls(i,2,P_uln +k-1)
!           trie_ls(i,1,P_w +k-1) = trie_ls(i,1,P_vln +k-1)
!           trie_ls(i,2,P_w +k-1) = trie_ls(i,2,P_vln +k-1)
!        enddo
!        do i=1,len_trio_ls
!           trio_ls(i,1,P_x +k-1) = trio_ls(i,1,P_uln +k-1)
!           trio_ls(i,2,P_x +k-1) = trio_ls(i,2,P_uln +k-1)
!           trio_ls(i,1,P_w +k-1) = trio_ls(i,1,P_vln +k-1)
!           trio_ls(i,2,P_w +k-1) = trio_ls(i,2,P_vln +k-1)
!        enddo
!     enddo
 
!     print *,' exit grid_to_spect_slg '
!!
      return
      end
