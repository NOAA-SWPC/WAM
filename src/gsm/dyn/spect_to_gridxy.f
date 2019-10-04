      subroutine spect_to_gridxy
     x    (trie_ls,trio_ls,
     x     syn_gr_a_1,syn_gr_a_2,
     x     dyn_gr_a_1,dyn_gr_a_2,
     x     ls_node,ls_nodes,max_ls_nodes,
     x     lats_nodes_a,global_lats_a,
     x     lonsperlat,
     x     pddev_a,pddod_a)
!
! program log:
! 20110220    Henry Juang update the code to fit mass_dp and ndslfv
! 20130620    Henry Juang add length for sumder in case of ndsl
!
!
      use gfs_dyn_resol_def
      use gfs_dyn_layout1
      use gfs_dyn_gg_def
      use gfs_dyn_vert_def
      use gfs_dyn_date_def
      use namelist_dynamics_def
      use gfs_dyn_mpi_def
      use gfs_dyn_physcons, rerth => con_rerth
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
      real(kind=kind_evod)   pddev_a(len_trie_ls,latg2)
      real(kind=kind_evod)   pddod_a(len_trio_ls,latg2)
!
      real(kind=kind_evod) syn_gr_a_1(lonfx*lots,lats_dim_a)
      real(kind=kind_evod) dyn_gr_a_1(lonfx*lotd,lats_dim_a)
!
      real(kind=kind_evod) syn_gr_a_2(lonfx*lots,lats_dim_a)
      real(kind=kind_evod) dyn_gr_a_2(lonfx*lotd,lats_dim_a)
!
      integer              i,j,k,l,n,lotx,lotxx
      integer              lon,lan,lat,lmax
      integer              lon_dim,lons_lat,node
!
      integer              ipt_ls
!
      real(kind=kind_evod) reall
      real(kind=kind_evod) rlcs2(jcap1)
!
      real(kind=kind_evod) cons0
!
      cons0 = 0.d0 
      if( .not. ndslfv ) then
        lotx  =   levs +   levh
        lotxx = 4*levs + 2*levh
      else
        lotx  =   levs
        lotxx = 4*levs
      endif
!
! ................................................................
!
      dimg=0
!
      CALL countperf(0,1,0.)
      call f_hpmstart(9,"ga sumdera")
!
      call sumdera(trie_ls(1,1,P_te),
     x             trio_ls(1,1,P_te),
     x             lat1s_a,
     x             pddev_a,pddod_a,
     x             lotx,lotxx,ls_node,latg2,
     x             lats_dim_a,lotd,
     x             dyn_gr_a_1,
     x             ls_nodes,max_ls_nodes,
     x             lats_nodes_a,global_lats_a,
     x             lats_node_a,ipt_lats_node_a,lon_dims_a,dimg,
     x             lonsperlat,lonfx,latg)
!
      call f_hpmstop(9)
      CALL countperf(1,1,0.)
!
      call f_hpmstart(10,"ga lat_loop")
!11111111111111111111111111111111111111111111111111111111111111111111
      do lan=1,lats_node_a   !sela begin lan loop 1
!
         lon_dim = lon_dims_a(lan)
         lat = global_lats_a(ipt_lats_node_a-1+lan)
         lons_lat = lonsperlat(lat)
!!
!!       calculate t rq u v zonal derivs. by multiplication with i*l
!!       note rlcs2=rcs2*L/rerth
!
         lmax = min(jcap,lons_lat/2)
!
         ipt_ls=min(lat,latg-lat+1)
 
         do i=1,lmax+1
            if ( ipt_ls .ge. lat1s_a(i-1) ) then
               reall=i-1
               rlcs2(i)=reall*rcs2_a(ipt_ls)/rerth
            else
               rlcs2(i)=cons0     !constant
            endif
         enddo
!
!$omp parallel do private(k,i)
         do k=1,levs
            do i=1,lmax+1
!
!           d(t)/d(lam)
               dyn_gr_a_1(2*i-1+(kdtlam-2+k)*lon_dim,lan)=
     x        -syn_gr_a_1(2*i  +(kst   -2+k)*lon_dim,lan)*rlcs2(i)
               dyn_gr_a_1(2*i  +(kdtlam-2+k)*lon_dim,lan)=
     x         syn_gr_a_1(2*i-1+(kst   -2+k)*lon_dim,lan)*rlcs2(i)
!
!           d(u)/d(lam)
               dyn_gr_a_1(2*i-1+(kdulam-2+k)*lon_dim,lan)=
     x        -syn_gr_a_1(2*i  +(ksu   -2+k)*lon_dim,lan)*rlcs2(i)
               dyn_gr_a_1(2*i  +(kdulam-2+k)*lon_dim,lan)=
     x         syn_gr_a_1(2*i-1+(ksu   -2+k)*lon_dim,lan)*rlcs2(i)
!
!           d(v)/d(lam)
               dyn_gr_a_1(2*i-1+(kdvlam-2+k)*lon_dim,lan)=
     x        -syn_gr_a_1(2*i  +(ksv   -2+k)*lon_dim,lan)*rlcs2(i)
               dyn_gr_a_1(2*i  +(kdvlam-2+k)*lon_dim,lan)=
     x         syn_gr_a_1(2*i-1+(ksv   -2+k)*lon_dim,lan)*rlcs2(i)
!
            enddo
         end do
!
         if( .not. ndslfv ) then

!$omp parallel do private(k,i)
         do k=1,levh
            do i=1,lmax+1
!
!           d(rq)/d(lam)
               dyn_gr_a_1(2*i-1+(kdrlam-2+k)*lon_dim,lan)=
     x        -syn_gr_a_1(2*i  +(ksr   -2+k)*lon_dim,lan)*rlcs2(i)
               dyn_gr_a_1(2*i  +(kdrlam-2+k)*lon_dim,lan)=
     x         syn_gr_a_1(2*i-1+(ksr   -2+k)*lon_dim,lan)*rlcs2(i)
!
            enddo
         enddo

         endif
!
         CALL countperf(0,6,0.)
         CALL FOUR2GRID_thread(dyn_gr_a_1(1,lan),dyn_gr_a_2(1,lan),
     &                  lon_dim,lons_lat,lonfx,lotxx,lan,me)
         CALL countperf(1,6,0.)
!
      enddo !sela fin lan loop 1
!11111111111111111111111111111111111111111111111111111111111111111111
!
!22222222222222222222222222222222222222222222222222222222222
      do lan=1,lats_node_a   !sela begin lan loop 2
!
         lat = global_lats_a(ipt_lats_node_a-1+lan)
!
         lon_dim = lon_dims_a(lan)
         lons_lat = lonsperlat(lat)
!!
!!  calculate grid meridional derivatives of u and v.
!!
!!  cos*d(u)/d(theta)= d(v)/d(lam)-a*zeta*cos**2
!!  cos*d(v)/d(theta)=-d(u)/d(lam)+a*divr*cos**2
!!
!$omp parallel do private(k,j)
      do k=1,levs
         do j=1,lons_lat
!
            dyn_gr_a_2(j+(kduphi-2+k)*lon_dim,lan)=
     x      dyn_gr_a_2(j+(kdvlam-2+k)*lon_dim,lan)-
     x      syn_gr_a_2(j+(ksz   -2+k)*lon_dim,lan)
!
            dyn_gr_a_2(j+(kdvphi-2+k)*lon_dim,lan)=
     x     -dyn_gr_a_2(j+(kdulam-2+k)*lon_dim,lan)+
     x      syn_gr_a_2(j+(ksd   -2+k)*lon_dim,lan)
!
         enddo
      enddo
!
      enddo
!
      return
      end
