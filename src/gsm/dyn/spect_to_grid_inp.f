      subroutine spect_to_grid_inp
     x    (trie_zs,trio_zs,trie_ps,trio_ps,
     x     trie_di,trio_di,trie_ze,trio_ze,
     x     trie_te,trio_te,trie_rq,trio_rq,
     &     zsg,psg,uug,vvg,teg,rqg,
     x     ls_node,ls_nodes,max_ls_nodes,
     x     lats_nodes_a,global_lats_a,lonsperlat,
     x     epsedn,epsodn,snnp1ev,snnp1od,plnev_a,plnod_a)
!!
!  routine to do spectral transform initially for grid-point input
!  initally coded by hann-ming henry juang
!
! program log:
! 2011 02 20 :    Henry Juang, add NDSL and MASS_DP options for semi-Lagrangian
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
      real(kind=kind_evod) trie_zs(len_trie_ls,2)
      real(kind=kind_evod) trio_zs(len_trio_ls,2)
      real(kind=kind_evod) trie_ps(len_trie_ls,2)
      real(kind=kind_evod) trio_ps(len_trio_ls,2)
      real(kind=kind_evod) trie_di(len_trie_ls,2,levs)
      real(kind=kind_evod) trio_di(len_trio_ls,2,levs)
      real(kind=kind_evod) trie_ze(len_trie_ls,2,levs)
      real(kind=kind_evod) trio_ze(len_trio_ls,2,levs)
      real(kind=kind_evod) trie_te(len_trie_ls,2,levs)
      real(kind=kind_evod) trio_te(len_trio_ls,2,levs)
      real(kind=kind_evod) trie_rq(len_trie_ls,2,levh)
      real(kind=kind_evod) trio_rq(len_trio_ls,2,levh)
!
!!!!  parameter            ( lota = 3*levs+1*levh+1 )

      real(kind=kind_evod) trie_ls(len_trie_ls,2,lota+1)
      real(kind=kind_evod) trio_ls(len_trio_ls,2,lota+1)
!!
      real(kind=kind_evod) for_gr_a_1(lonfx*(lota+1),lats_dim_a)
      real(kind=kind_evod) for_gr_a_2(lonfx*(lota+1),lats_dim_a)
cc
      real(kind=kind_grid) zsg(lonf,lats_node_a)
      real(kind=kind_grid) psg(lonf,lats_node_a)
      real(kind=kind_grid) uug(lonf,lats_node_a,levs)
      real(kind=kind_grid) vvg(lonf,lats_node_a,levs)
      real(kind=kind_grid) teg(lonf,lats_node_a,levs)
      real(kind=kind_grid) rqg(lonf,lats_node_a,levh)
cc
      integer              ls_node, ls_nodes(ls_dim,nodes)
      integer              max_ls_nodes(nodes)
      integer              lats_nodes_a(nodes)
      integer              global_lats_a(latg)
      integer                 lonsperlat(latg)
      integer dimg
cc
      real(kind=kind_evod)  epsedn(len_trie_ls)
      real(kind=kind_evod)  epsodn(len_trio_ls)
cc
      real(kind=kind_evod) snnp1ev(len_trie_ls)
      real(kind=kind_evod) snnp1od(len_trio_ls)
cc
      real(kind=kind_evod)   plnev_a(len_trie_ls,latg2)
      real(kind=kind_evod)   plnod_a(len_trio_ls,latg2)
cc
      integer              i,j,k
      integer              l,lan,lat,lotdim,lotx
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
      lotdim  = lota+1
      lotx    = 4*levs+levh+2

cc
cc--------------------------------------------------------------------
cc
!$OMP parallel do shared(trie_ls,trio_ls)
!$OMP+shared(trie_di,trie_ze,trie_te,kau,kav,kat)
!$OMP+shared(trio_di,trio_ze,trio_te)
!$OMP+shared(epsedn,epsodn,snnp1ev,snnp1od,ls_node)
!$OMP+private(k)
      do k=1,levs
        call dezouv(trie_di(1,1,k),       trio_ze(1,1,k),
     x              trie_ls(1,1,kau+k-1), trio_ls(1,1,kav+k-1),
     x              epsedn,epsodn,snnp1ev,snnp1od,ls_node)
cc
        call dozeuv(trio_di(1,1,k),       trie_ze(1,1,k),
     x              trio_ls(1,1,kau+k-1), trie_ls(1,1,kav+k-1),
     x              epsedn,epsodn,snnp1ev,snnp1od,ls_node)
        trie_ls(:,:,kat+k-1)=trie_te(:,:,k)
        trio_ls(:,:,kat+k-1)=trio_te(:,:,k)
      enddo
      do k=1,levh
        trie_ls(:,:,kar+k-1)=trie_rq(:,:,k)
        trio_ls(:,:,kar+k-1)=trio_rq(:,:,k)
      enddo
      trie_ls(:,:,kazs)=trie_zs(:,:)
      trio_ls(:,:,kazs)=trio_zs(:,:)
      trie_ls(:,:,kaps)=trie_ps(:,:)
      trio_ls(:,:,kaps)=trio_ps(:,:)
!!
      dimg=0
cc
      call sumflna(trie_ls,trio_ls,
     x            lat1s_a,
     x            plnev_a,plnod_a,
     x            lotx,ls_node,latg2,
     x            lats_dim_a,lotdim,for_gr_a_1,
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
         call four2grid_thread(for_gr_a_1(1,lan),for_gr_a_2(1,lan),
     &                  lon_dim,lons_lat,lonfx,lotx,lan,me)
!
      enddo   !lan
cc
      do lan=1,lats_node_a
        lon_dim = lon_dims_a(lan)
        lat = global_lats_a(ipt_lats_node_a-1+lan)
        lons_lat = lonsperlat(lat)
        do k=1,levs
          do i=1,lons_lat
            uug(i,lan,k)=for_gr_a_2(i+(kau+k-2)*lon_dim,lan)
            vvg(i,lan,k)=for_gr_a_2(i+(kav+k-2)*lon_dim,lan)
            teg(i,lan,k)=for_gr_a_2(i+(kat+k-2)*lon_dim,lan)
          enddo
        enddo
        do k=1,levh
          do i=1,lons_lat
            rqg(i,lan,k)=for_gr_a_2(i+(kar+k-2)*lon_dim,lan)
          enddo
        enddo
        do i=1,lons_lat
          zsg(i,lan)=for_gr_a_2(i+(kazs-1)*lon_dim,lan)
          psg(i,lan)=for_gr_a_2(i+(kaps-1)*lon_dim,lan)
        enddo
      enddo
cc
!     print *,' exit spect_to_grid_inp '
!!
      return
      end
