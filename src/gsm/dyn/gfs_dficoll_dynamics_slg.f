!
! ==========================================================================
!
      subroutine gfs_dficoll_dynamics_slg(for_grid_a,ug_h,vg_h,
     &                                    kst, ksq, ksrg,
     &                                    ak5, bk5,
     &                                    global_lats_a,lonsperlat,
     &                                    gz_grid,grid_gr_dfi)
!       
!    Nov 2015  S Moorthi Initial code from gfs_dficoll_dynamics (by Jun W)
!                        semi-Lagrangian version
!    Feb 2016  S Moorthi Add  code for computing p and dp
!
      use namelist_dynamics_def,  only : ngptc
      use gfs_dyn_machine      ,  only : kind_evod, kind_grid
      use gfs_dyn_resol_def,      ONLY : lonf, levs, latg, levh
      use layout_lag            , only : lats_dim_h,lon_dim_h
      use gfs_dyn_layout1,        ONLY : lats_node_a, lats_dim_a
     &,                                  lon_dims_a
     &,                                  lon_dim_a, ipt_lats_node_a, me
!    &,                                  nodes_comp, me, ipt_lats_node_a
      use layout_grid_tracers ,   ONLY : xhalo, yhalo

      use gfs_dyn_dfi_mod,        only : gfs_dfi_grid_gr
!     use gfs_dyn_layout1, only : lats_node_a_max
!     use gfs_dyn_resol_def
!
      implicit none
!
!     REAL(KIND=KIND_EVOD)  :: grid_gr(lonf,lats_node_a_max,lotgr)

      integer, intent(in)   ::  kst, ksq, ksrg
      INTEGER, dimension(latg) :: GLOBAL_lats_a, lonsperlat
      real(kind=kind_evod)  :: for_grid_a(lonf,2*levs+levh+3,lats_dim_a)
     &,                        ug_h(lon_dim_h,levs,lats_dim_h)
     &,                        vg_h(lon_dim_h,levs,lats_dim_h)
     &,                        gz_grid(lon_dim_a,lats_dim_a)
      type(gfs_dfi_grid_gr) :: grid_gr_dfi
      real(kind=kind_evod), dimension(levs+1) :: ak5, bk5

!local vars
      real(kind=kind_grid)  gq(ngptc)
      real(kind=kind_grid), dimension(ngptc,levs) :: prsl, dprs
      real(kind=kind_evod), dimension(lonf)   :: tx1,  wrk
      real    :: tem1
      integer :: i,j,k,kk,lons_lat,lat,jj,lon,njeff
!     integer :: i,j,k,kk,kstr,kend,krq
!
!1: gz
      if(grid_gr_dfi%z_imp == 1) then
!$omp parallel do private(i,j)
        do j=1,lats_node_a
          do i=1,lonf
            grid_gr_dfi%hs(i,j,1) = gz_grid(i,j)
          enddo
        enddo
!      print*,'in dficoll,hs=',maxval(grid_gr_dfi%hs),
!    &                         minval(grid_gr_dfi%hs) 
      endif
!
!2: ps
      if(grid_gr_dfi%ps_imp == 1) then
!!$omp parallel do private(i,j,latlon_lat)
!$omp parallel do private(i,j,jj)
        do j=1,lats_node_a
          jj = lats_node_a + 1 -j
!         lat      = global_lats_a(ipt_lats_node_a-1+j)
!         lon_lat = lonsperlat(lat)
!         do i=1,lon_lat
          do i=1,lonf
!           grid_gr_dfi%ps(i,j,1) = exp(for_grid_a(i,ksq,j))
            grid_gr_dfi%ps(i,j,1) = for_grid_a(i,ksq,jj)
          enddo
        enddo
!      print*,'in dficoll,ps=',maxval(grid_gr_dfi%ps),
!    &                         minval(grid_gr_dfi%ps)   
      endif

!3: t
      if(grid_gr_dfi%temp_imp == 1) then
!$omp parallel do private(i,j,k,kk,jj)
        do j=1,lats_node_a
          jj = lats_node_a + 1 -j
          do k=1, levs 
            kk = kst + k -1
            do i=1,lonf
              grid_gr_dfi%t(i,j,k) = for_grid_a(i,kk,jj)
            enddo
          enddo
        enddo
         
!      print*,'in dficoll,t='
!    &,        maxval(grid_gr_dfi%t(1:lonf,1:lats_node_a,:))
!    &,        minval(grid_gr_dfi%t(1:lonf,1:lats_node_a,:))
      endif

!4: u
      if(grid_gr_dfi%u_imp == 1) then

!$omp parallel do private(i,j,k,jj)
        do j=1,lats_node_a
          jj = lats_node_a + 1 -j
          do k=1,levs
            do i=1,lonf
              grid_gr_dfi%u(i,jj,k) = ug_h(i+xhalo,k,j+yhalo)
            enddo
          enddo
        enddo
!      print*,'in dficoll,u=',maxval(grid_gr_dfi%u)
!    &,                       minval(grid_gr_dfi%u)
      endif

!5: v
      if(grid_gr_dfi%v_imp == 1) then
!$omp parallel do private(i,j,k,jj)
        do j=1,lats_node_a
          jj = lats_node_a + 1 -j
          do k=1,levs
            do i=1,lonf
              grid_gr_dfi%v(i,jj,k) = vg_h(i+xhalo,k,j+yhalo)

            enddo
          enddo
        enddo
!      print*,'in dficoll,v=',maxval(grid_gr_dfi%v)
!    &,                       minval(grid_gr_dfi%v)
      endif

!6: tracer
!       print *,'in dficoll0,tracer=',size(grid_gr_dfi%tracer,3),
!     &   'ntrac=',ntrac,'levs=',levs,'g_rt=',g_rt

      if(grid_gr_dfi%tracer_imp == 1) then
!$omp parallel do private(i,j,k,kk,jj)
        do k=1,levh
          kk = ksrg + k - 1
          do j=1,lats_node_a
            jj = lats_node_a + 1 -j
            do i=1,lonf
              grid_gr_dfi%tracer(i,j,k) = for_grid_a(i,kk,jj)
            enddo
          enddo
!         write(0,*)' in dficoll k=',k,' kk==',kk,
!    &          maxval(grid_gr_dfi%tracer(1:lonf,1:lats_node_a,k)),
!    &          minval(grid_gr_dfi%tracer(1:lonf,1:lats_node_a,k))
        enddo
!        print*,'in dficoll,q=',maxval(grid_gr_dfi%tracer),
!    &                          minval(grid_gr_dfi%tracer)
      endif
!       print *,'in dficoll'

!7: p
      if(grid_gr_dfi%p_imp == 1) then
        do j=1,lats_node_a
          jj = lats_node_a + 1 -j
          lat      = global_lats_a(ipt_lats_node_a-1+j)
!         lat      = global_lats_a(ipt_lats_node_a+lats_node_a-j)
          lons_lat = lonsperlat(lat)
!!$omp parallel do private(i)
!         do i=1,lons_lat
!!         write(0,*)' for_grid_a=',for_grid_a(i,ksq,j),' i=',i,' j=',j
!           tx1(i) = exp(for_grid_a(i,ksq,jj))
!           wrk(i) = ak5(1) + bk5(1) * tx1(i)
!         enddo
!$omp parallel do schedule(dynamic,1) private(lon,i,k,njeff)
!$omp+private(gq,prsl,dprs)
          do lon=1,lons_lat,ngptc
            njeff = min(ngptc,lons_lat-lon+1)
            do i=1,njeff
              gq(i) = for_grid_a(lon+i-1,ksq,jj)
            enddo
            call hyb2press(njeff,ngptc,gq, prsl, dprs)
            do k=1,levs
              do i=1,njeff
                grid_gr_dfi%p(lon+i-1,j,k)  = prsl(i,k)
                grid_gr_dfi%dp(lon+i-1,j,k) = dprs(i,k)
              enddo
            enddo
!     write(3000+me,*)' indif_col prsl=',prsl(njeff,:),' j=',j
!     write(3000+me,*)' indif_col dprs=',prsl(njeff,:),' j=',j
          enddo
!           do k=1, levs
!           kk = levs + 1 -k
!!$omp parallel do private(i,tem1)
!           do i=1,lons_lat
!             tem1   = wrk(i)
!             wrk(i) = ak5(k+1) + bk5(k+1) * tx1(i)
!             grid_gr_dfi%p(i,j,kk)  = wrk(i)
!             grid_gr_dfi%dp(i,j,kk) = wrk(i) - tem1
!           enddo
!         enddo
        enddo
!       kstr=g_p
!       kend=kstr + levs - 1
!       grid_gr_dfi%p=grid_gr(:,:,kstr:kend)
!      print*,'in dficoll,p=',maxval(grid_gr_dfi%p),minval(grid_gr_dfi%p)
!      print*,'p_imp == 1 option is not coded for slg digital filter'
      endif
!8: dp
!     if(grid_gr_dfi%dp_imp == 1) then
!!$omp parallel do private(i,j,k)
!       do j=1,lats_node_a
!         do k=1, levs
!           kk = g_dp + k -1
!           do i=1,lonf
!             grid_gr_dfi%dp(i,j,k) = for_grid_a(i,kk,j)
!           enddo
!         enddo
!       kstr=g_dp
!       kend=kstr + levs - 1
!      grid_gr_dfi%dp=grid_gr(:,:,kstr:kend)
!      print*,'in dficoll,p=',maxval(grid_gr_dfi%p),minval(grid_gr_dfi%p)
!      print*,'dp_imp == 1 option is not coded for slg digital filter'
!     endif
!9: dpdt
!     if(grid_gr_dfi%dpdt_imp == 1) then
!!$omp parallel do private(i,j,k)
!       do j=1,lats_node_a
!         do k=1, levs
!!          kk = g_dpdt + k -1
!           do i=1,lonf
!             grid_gr_dfi%dpdt(i,j,k) = for_grid_a(i,kk,j)
!             grid_gr_dfi%dpdt(i,j,k) = 0.0
!           enddo
!         enddo
!       enddo
!       kstr=g_dpdt
!       kend=kstr + levs - 1
!       grid_gr_dfi%dpdt=grid_gr(:,:,kstr:kend)
!      print*,'in dficoll,p=',maxval(grid_gr_dfi%dpdt),
!     &    minval(grid_gr_dfi%dpdt)
!      print*,'dpdt_imp == 1 option is not coded for slg digital filter'
!     endif

!       print *,'end of dfi coll'
!
      end subroutine gfs_dficoll_dynamics_slg
