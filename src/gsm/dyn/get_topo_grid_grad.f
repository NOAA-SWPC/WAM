      subroutine get_topo_grid_grad(cosf1,
     &                              grad_gzlam,grad_gzphi,
     &                              gzie_ln,gzio_ln,
     &                              gz_grid,
     &                              phige,phigo,
     &                              for_gr_a_1,for_gr_a_2,
     &                              trie_ls,trio_ls,
     &                              ls_node,ls_nodes,max_ls_nodes,
     &                              lats_nodes_a,global_lats_a,
     &                              lonsperlat,
     &                              epse,epso,
     &                              plnev_a,plnod_a,
     &                              i_write,r_rt)

      use gfs_dyn_machine  , only : kind_evod, kind_io8
      use gfs_dyn_resol_def, only : jcap,latg,latg2,levh,levs,lnt2,lonf,
     &                              p_dlam,p_dphi,p_zq
      use gfs_dyn_layout1  , only : ipt_lats_node_a,lat1s_a,lats_dim_a,
     &                              lats_node_a,len_trie_ls,len_trio_ls,
     &                              lon_dim_a,ls_dim,ls_max_node,me,
     &                              nodes
      use gfs_dyn_io_header, only : z
      use gfs_dyn_physcons ,only :  grav => con_g, rd => con_rd
      use namelist_dynamics_def, only : ref_temp

      implicit none

      real(kind=kind_evod)      cosf1(          lats_dim_a)
      real(kind=kind_evod), dimension(lon_dim_a,lats_dim_a) ::
     &                           grad_gzlam, grad_gzphi, gz_grid

      real(kind=kind_evod), dimension(len_trie_ls,2) :: gzie_ln, phige
      real(kind=kind_evod), dimension(len_trio_ls,2) :: gzio_ln, phigo


      real(kind=kind_evod) for_gr_a_1(lon_dim_a,2,lats_dim_a)
!mjr  real(kind=kind_evod) for_gr_a_2(lon_dim_a,2,lats_dim_a)
      real(kind=kind_evod) for_gr_a_2(lonf     ,2,lats_dim_a)

      real(kind=kind_evod) trie_ls(len_trie_ls,2,11*levs+3*levh+6)
      real(kind=kind_evod) trio_ls(len_trio_ls,2,11*levs+3*levh+6)

      integer              ls_node(ls_dim,3)
              
!     ls_node(1,1) ... ls_node(ls_max_node,1) : values of L
!     ls_node(1,2) ... ls_node(ls_max_node,2) : values of jbasev
!     ls_node(1,3) ... ls_node(ls_max_node,3) : values of jbasod

      integer              ls_nodes(ls_dim,nodes)
      integer, dimension(nodes) :: max_ls_nodes, lats_nodes_a
      integer, dimension(latg)  ::  global_lats_a, lonsperlat

      real(kind=kind_evod)    epse(len_trie_ls)
      real(kind=kind_evod)    epso(len_trio_ls)

      real(kind=kind_evod)   plnev_a(len_trie_ls,latg2)
      real(kind=kind_evod)   plnod_a(len_trio_ls,latg2)

      real(kind=kind_evod) trisca(lnt2)
      real(kind=kind_evod) tref,r_rt, tem 

      integer i, l, lan, lat, locl, lons_lat, n, i_write
     &,       indlsev,jbasev,indev,indev1,indev2
     &,       indlsod,jbasod,indod,indod1,indod2

!     include 'function2'
      include 'function_indlsev'
      include 'function_indlsod'

! --------------------------------------------------------------------

!    begin  calculation of grad(gz) +++++++++++++++++++++++++

      tref = ref_temp
!sela r_rt = 1./( rd*tref)

      trisca = z*grav ! z is in module sig_io ; set in treadeo

!     print*,' grav in get_topo_grid_grad =',grav,' r_rt=',r_rt

      call triseori(trisca,trie_ls(1,1,p_zq),
     &                     trio_ls(1,1,p_zq),1,ls_node)
       do i=1,len_trie_ls
        gzie_ln(i,1) = trie_ls(i,1,p_zq)*r_rt
        gzie_ln(i,2) = trie_ls(i,2,p_zq)*r_rt
       enddo
       do i=1,len_trio_ls
        gzio_ln(i,1) = trio_ls(i,1,p_zq)*r_rt
        gzio_ln(i,2) = trio_ls(i,2,p_zq)*r_rt
       enddo

!---------------------------------------------------------------------------------
!     if ( i_write == 1 ) then
!       do locl=1,ls_max_node
!          write(9000+me,'(" ")')
!               L=ls_node(locl,1)
!          jbasev=ls_node(locl,2)
!          indev1 = indlsev(L,L)
!          if (mod(L,2).eq.mod(jcap+1,2)) then
!             indev2 = indlsev(jcap+1,L)
!          else
!            indev2 = indlsev(jcap  ,L)
!          endif
!          do indev = indev1 , indev2
!             write(9000+me,
!    &        '("L=",i6," indev=",i6," trie_ls=",2e25.15)')
!    &           L, indev, trie_ls(indev,1,p_zq),
!    &                     trie_ls(indev,2,p_zq)
!          end do
!       end do

!       do locl=1,ls_max_node
!          write(9000+me,'(" ")')
!               L=ls_node(locl,1)
!          jbasod=ls_node(locl,3)
!          indod1 = indlsod(L+1,L)
!          if (mod(L,2).eq.mod(jcap+1,2)) then
!             indod2 = indlsod(jcap  ,L)
!          else
!             indod2 = indlsod(jcap+1,L)
!          endif
!          do indod = indod1 , indod2
!             write(9000+me,
!    &        '("L=",i6," indod=",i6," trio_ls=",2e25.15)')
!    &           L, indod, trio_ls(indod,1,p_zq),
!    &                     trio_ls(indod,2,p_zq)
!          end do
!       end do
!       close(9000+me)
!     endif  !  if ( i_write == 1 )
!---------------------------------------------------------------------------------

      call delnpe(trie_ls(1,1,p_zq  ), trio_ls(1,1,p_dphi),
     &            trie_ls(1,1,p_dlam), epse,epso,ls_node)

      call delnpo(trio_ls(1,1,p_zq  ), trie_ls(1,1,p_dphi),
     &            trio_ls(1,1,p_dlam), epse,epso,ls_node)

      call sumfln_slg_gg(trie_ls(1,1,p_dlam),
     &                   trio_ls(1,1,p_dlam),
     &                   lat1s_a,
     &                   plnev_a,plnod_a,
     &                   2     ,ls_node,latg2,
     &                   lats_dim_a,2,
     &                   for_gr_a_1,
     &                   ls_nodes,max_ls_nodes,
     &                   lats_nodes_a,global_lats_a,
     &                   lats_node_a,ipt_lats_node_a,lon_dim_a,
     &                   lonsperlat,lon_dim_a,latg,0)

      do lan=1,lats_node_a
        lat = global_lats_a(ipt_lats_node_a-1+lan)
        lons_lat = lonsperlat(lat)
 
        CALL FOUR_TO_GRID(for_gr_a_1(1,1,lan),
     &                    for_gr_a_2(1,1,lan),
     &                    lon_dim_a,lon_dim_a-2,lons_lat,2)

        tem = 1.0 / cosf1(lan)
        do i=1,lons_lat
          grad_gzlam(i,lats_node_a+1-lan) = for_gr_a_2(i,1,lan) * tem
          grad_gzphi(i,lats_node_a+1-lan) = for_gr_a_2(i,2,lan) * tem
        enddo

      enddo ! do lan

!---------------------------------------------------------------------------------
!     if ( i_write == 1 ) then

!       do lan=1,lats_node_a
!         lat = global_lats_a(ipt_lats_node_a-1+lan)
!         lons_lat = lonsperlat(lat)

!         write(9400+me,'(" ")')
!         do i=1,lons_lat
!           write(9400+me,
!    &           '("lan=",i6," i=",i6," grad lam phi=",2e25.15)')
!    &              lan, i, grad_gzlam(i,lan),
!    &                      grad_gzphi(i,lan)
!         enddo

!       enddo ! do lan

!       close(9400+me)

!     endif  !  if ( i_write == 1 )
!---------------------------------------------------------------------------------

!    fin get_topo_grad_resonan grad(gz)  and gzie_ln, gzio_ln

!--------------------------------------------------------------------

!mjr  trisca=z*grav ! z is in module sig_io ; set in treadeo
!mjr  print*,' grav in get_topo_grid_grad =',grav

      call triseori(trisca,phige(1,1),phigo(1,1),1,ls_node)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     if ( i_write == 1 ) then
!       do locl=1,ls_max_node
!         write(8000+me,'(" ")')
!              l=ls_node(locl,1)
!         jbasev=ls_node(locl,2)
!         indev1 = indlsev(l,l)
!         if (mod(l,2).eq.mod(jcap+1,2)) then
!            indev2 = indlsev(jcap+1,l)
!         else
!            indev2 = indlsev(jcap  ,l)
!         endif
!         do indev = indev1 , indev2
!            write(8000+me,
!    &       '("l=",i6," indev=",i6," phige=",2e25.15)')
!    &          l, indev, phige(indev,1),
!    &                    phige(indev,2)
!         end do
!       end do
!       do locl=1,ls_max_node
!         write(8000+me,'(" ")')
!              l=ls_node(locl,1)
!         jbasod=ls_node(locl,3)
!         indod1 = indlsod(l+1,l)
!         if (mod(l,2).eq.mod(jcap+1,2)) then
!            indod2 = indlsod(jcap  ,l)
!         else
!          indod2 = indlsod(jcap+1,l)
!         endif
!         do indod = indod1 , indod2
!            write(8000+me,
!    &       '("l=",i6," indod=",i6," phigo=",2e25.15)')
!    &          l, indod, phigo(indod,1),
!    &                    phigo(indod,2)
!         end do
!       end do
!       close(8000+me)
!     endif  !  if ( i_write == 1 )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call sumfln_slg_gg(phige(1,1),
     &                   phigo(1,1),
     &                   lat1s_a,
     &                   plnev_a,plnod_a,
     &                   1     ,ls_node,latg2,
     &                   lats_dim_a,2,
     &                   for_gr_a_1,
     &                   ls_nodes,max_ls_nodes,
     &                   lats_nodes_a,global_lats_a,
!mjr &                   lats_node_a,ipt_lats_node_a,lon_dims_a,
     &                   lats_node_a,ipt_lats_node_a,lon_dim_a,
     &                   lonsperlat,lon_dim_a,latg,0)

      do lan=1,lats_node_a
        lat = global_lats_a(ipt_lats_node_a-1+lan)
!mjr    lon_dim = lon_dims_a(lan)
        lons_lat = lonsperlat(lat)
        call four_to_grid(for_gr_a_1(1,1,lan),
     &                    for_gr_a_2(1,1,lan),
!mjr &                    lon_dim,  lon_dim    ,lons_lat,1)
     &                    lon_dim_a,lon_dim_a-2,lons_lat,1)

        do i=1,lons_lat
!         gz_grid(i,lan) = for_gr_a_2(i,1,lan)

          gz_grid(i,lats_node_a+1-lan) = for_gr_a_2(i,1,lan)
        enddo
      enddo ! do lan

!     if ( i_write == 1 ) then
!       do lan=1,lats_node_a
!         lat = global_lats_a(ipt_lats_node_a-1+lan)
!         lons_lat = lonsperlat(lat)
!         write(8400+me,'(" ")')
!         do i=1,lons_lat
!           write(8400+me,
!    &       '("lan=",i6," i=",i6," gz_grid  =",2e25.15)')
!    &          lan, i, gz_grid(i,lan)
!         enddo
!       enddo ! do lan
!       close(8400+me)
!     endif  !  if ( i_write == 1 )

      return
      end
