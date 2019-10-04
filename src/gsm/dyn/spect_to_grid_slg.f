      subroutine spect_to_grid_slg
     &    (trie_ls,trio_ls,
     &     syn_gr_a_1,syn_gr_a_2,
     &     ls_node,ls_nodes,max_ls_nodes,
     &     lats_nodes_a,global_lats_a,lonsperlat,
     &     epsedn,epsodn,
     &     snnp1ev,snnp1od,plnev_a,plnod_a)
!
! H.-M. H. Juang: get from gloopa to have only spect to grid
!                 syn_gr_a_1 and syn_gr_a_2 will be used in
!                 the next spect_to_gridxy for saving time
! program log
! 2011 02 20 : henry jaung, updated code for mass_dp and ndsl advection
! 2014 08 -- : shrinivas moorthi - updated for sela-semi-lagrangian
!
      use gfs_dyn_resol_def
      use gfs_dyn_layout1
      use gfs_dyn_gg_def
      use gfs_dyn_vert_def
      use gfs_dyn_date_def
      use namelist_dynamics_def
      use gfs_dyn_mpi_def
      use layout_grid_tracers , only : rgt_a,rgt_h,xhalo,yhalo
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
      real(kind=kind_evod)  epsedn(len_trie_ls)
      real(kind=kind_evod)  epsodn(len_trio_ls)
!
      real(kind=kind_evod) snnp1ev(len_trie_ls)
      real(kind=kind_evod) snnp1od(len_trio_ls)
!
      real(kind=kind_evod)   plnev_a(len_trie_ls,latg2)
      real(kind=kind_evod)   plnod_a(len_trio_ls,latg2)
!
!     real(kind=kind_evod) syn_gr_a_1(lonfx*lots,lats_dim_a)
!     real(kind=kind_evod) syn_gr_a_2(lonfx*lots,lats_dim_a)
!
      real(kind=kind_evod) syn_gr_a_1(lon_dim_a,lota,lats_dim_a)
      real(kind=kind_evod) syn_gr_a_2(lonf,lota,lats_dim_a)
!
      integer              i,j,k,l,item,jtem
      integer              lan,lat
      integer              lon_dim,lons_lat,n,node,lots_l
!
! xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
      real(kind=kind_evod) cons0,cons2     !constant
!
      cons0 = 0.d0     !constant
      cons2 = 2.d0     !constant
!     lots_l = 5*levs+levh+3
!
!$omp parallel do shared(trie_ls,trio_ls)
!$omp+shared(epsedn,epsodn,snnp1ev,snnp1od,ls_node)
!$omp+private(k)
      do k=1,levs
!        call dezouv(trie_ls(1,1,P_di +k-1), trio_ls(1,1,P_ze +k-1),
         call dezouv(trie_ls(1,1,P_x +k-1), trio_ls(1,1,P_w +k-1),
     x               trie_ls(1,1,P_uln+k-1), trio_ls(1,1,P_vln+k-1),
     x               epsedn,epsodn,snnp1ev,snnp1od,ls_node)
!
!        call dozeuv(trio_ls(1,1,P_di +k-1), trie_ls(1,1,P_ze +k-1),
         call dozeuv(trio_ls(1,1,P_x +k-1), trie_ls(1,1,P_w +k-1),
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
      call f_hpmstart(8,"ga sumflna_slg")
!
!     call sumflna_slg(trie_ls(1,1,P_ze),
!     call sumfln_slg_gg(trie_ls(1,1,P_ze),
!    &                   trio_ls(1,1,P_ze),
!    &                   lat1s_a,
!    &                   plnev_a,plnod_a,
!    &                   lots,lots_slg,ls_node,latg2,
!    &                   lats_dim_a,lots,lots_slg,
!    &                   syn_gr_a_1,
!    &                   ls_nodes,max_ls_nodes,
!    &                   lats_nodes_a,global_lats_a,
!    &                   lats_node_a,ipt_lats_node_a,lon_dims_a,dimg,
!    &                   lonsperlat,lonfx,latg)
!
!     write(150+me,*)' p_q=',p_q,' ls_node=',ls_node,' lats_dim_a=',
!    &lats_dim_a,' lots_slg=',lots_slg,' lon_dim_a=',lon_dim_a,' latg=',
!    &latg,latg2,' lats_node_a=',lats_node_a
!    &,' trie_ls=',trie_ls(1,1,p_q),trie_ls(1,1,p_rt),p_q,p_rt
!    &,' me=',me,' in spect_to_grid'
!     write(150+me,*)'trie_ls_clw=',trie_ls(1,1,p_rt+levh-1)
      call sumfln_slg_gg(trie_ls(1,1,p_q),
     &                   trio_ls(1,1,p_q),
     &                   lat1s_a,
     &                   plnev_a,plnod_a,
     &                   5*levs+3,ls_node,latg2,
!    &                   lats_dim_a,lots_slg,syn_gr_a_1,
     &                   lats_dim_a,lota,syn_gr_a_1,
     &                   ls_nodes,max_ls_nodes,
     &                   lats_nodes_a,global_lats_a,
     &                   lats_node_a,ipt_lats_node_a,lon_dim_a,
     &                   lonsperlat,lon_dim_a,latg,0)

!
      if(.not.gg_tracers)then
        call sumfln_slg_gg(trie_ls(1,1,p_rt),
     &                     trio_ls(1,1,p_rt),
     &                     lat1s_a,
     &                     plnev_a,plnod_a,
     &                     levh,ls_node,latg2,
!    &                     lats_dim_a,lots_slg,syn_gr_a_1,
     &                     lats_dim_a,lota,syn_gr_a_1,
     &                     ls_nodes,max_ls_nodes,
     &                     lats_nodes_a,global_lats_a,
     &                     lats_node_a,ipt_lats_node_a,lon_dim_a,
     &                     lonsperlat,lon_dim_a,latg,5*levs+3)
      endif ! if(.not.gg_tracers)then  call f_hpmstop(8)
      CALL countperf(1,1,0.)
!
      call f_hpmstart(10,"ga lat_loop")
!
!11111111111111111111111111111111111111111111111111111111111111111111
      do lan=1,lats_node_a  
 
         lon_dim = lon_dims_a(lan)
         lat = global_lats_a(ipt_lats_node_a-1+lan)
         lons_lat = lonsperlat(lat)

!        CALL countperf(0,6,0.)
!        CALL FOUR2GRID_thread(syn_gr_a_1(1,lan),syn_gr_a_2(1,lan),
!    &                  lon_dim,lons_lat,lonfx,lots,lan,me)
!        CALL countperf(1,6,0.)

        call four_to_grid(syn_gr_a_1(1,1,lan),syn_gr_a_2(1,1,lan),
     &                    lon_dim_a,lon_dim_a-2,lons_lat,5*levs+3)
!        do i=1,min(lonf,lons_lat)
!!         syn_gr_a_2(i,ksq,lan) = exp(syn_gr_a_2(i,ksq,lan))
!          syn_gr_a_2(i,1,lan) = exp(syn_gr_a_2(i,1,lan))
!        enddo
!     write(0,*)' lon_dim_a=',lon_dim_a,' lons_lat=',lons_lat
!    &,' syn_gr_a_2=',syn_gr_a_2(1,1,lan),' lan=',lan,' me=',me

        ksr = 5*levs + 4
        if(.not.gg_tracers)then
          call four_to_grid(syn_gr_a_1(1,ksr,lan),
     &                      syn_gr_a_2(1,ksr,lan),
     &                      lon_dim_a,lon_dim_a-2,lons_lat,levh)
!     write(150+me,*)' lon_dim_a=',lon_dim_a,' lons_lat=',lons_lat
!    &,' syn_gr_a_2=',syn_gr_a_2(1,ksr+levh-1,lan),' lan=',lan,' me=',me
!    &,' syn_gr_a_1=',syn_gr_a_1(1,ksr+levh-1,lan)
        else  ! gg_tracers
!         if (.not.shuff_lats_r) then
!                  set for_gr_a_2 to rg1_a rg2_a rg3_a from gloopa
            do n=1,ntrac
!$omp parallel do private(k,i,item,jtem)
              do k=1,levs
                item = ksr - 1 + k +(n-1)*levs
                jtem = lats_node_a+1-lan
                do i=1,min(lonf,lons_lat)
                  syn_gr_a_2(i,item,lan) = rgt_a(i,k,jtem,n)
                enddo
              enddo
            enddo
!         endif              ! not shuff_lats_r
        endif                ! gg_tracers

      enddo 
!
!     if(gg_tracers .and. shuff_lats_r) then
!        print*,' gloopb mpi_tracers_a_to_b shuff_lats_r',shuff_lats_r
!        call mpi_tracers_a_to_b(rgt_a,lats_nodes_a,global_lats_a,
!    &                           syn_gr_a_2(1,1,1),
!    &                           lats_nodes_a,global_lats_a,ksr,0)
!     endif                       ! gg_tracers .and.  shuff_lats_r
!11111111111111111111111111111111111111111111111111111111111111111111
!
      return
      end subroutine spect_to_grid_slg
