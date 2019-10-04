      subroutine grid_to_spect_inp_slg
     &                     (zsg,psg,uug,vvg,ttg,rqg,
!    &                     (zsg,psg,uug,vvg,ttg,rqg,dpg,
     &                      GZE,QE,TEE,DIE,ZEE,RQE,
     &                      GZO,QO,TEO,DIO,ZEO,RQO,
!    &                      trie_ls,trio_ls,
     &                      ls_node,ls_nodes,max_ls_nodes,
     &                      lats_nodes_a,global_lats_a,lonsperlat,
     &                      epse,epso,plnew_a,plnow_a)
!    &,                     plnev_a,plnod_a,pwat,ptot,ptrc)
!!
!! hmhj - this routine do grid to spectral transform 
!!        from nemsio read in field, to model fields
!! input zsg,psg,uug,vvg,ttg,rqg (mapping wind, temp)
!! output zsg,psg,uug,vvg,ttg,rqg in model values (mapping wind, enthalpy)
!! aug 2010      sarah lu, modified to compute tracer global sum
!! feb 2011      henry juang updated to fit mass_dp and ndslfv
!! feb 2015      Jun Wang  add option for slg_flag
!! feb 2015      S Moorthi Changed Jun's modifications and simplified to
!!                         do just for Joe semi-Lagrangian.  Added gg_tracer option
!!                         removed slg_flag and many other cleanup.
!!
      use gfs_dyn_resol_def, only : lonf, latg, levs, levh, ntrac, latg2
      use gfs_dyn_layout1
      use gfs_dyn_gg_def
      use gfs_dyn_vert_def
      use gfs_dyn_date_def
      use namelist_dynamics_def
      use gfs_dyn_coordinate_def 
      use gfs_dyn_tracer_const
      use layout_lag          , only : lat1s_h,lats_dim_h,lon_dim_h
      use layout_grid_tracers , only : rgt_a,rgt_h,xhalo,yhalo
!     use gfs_dyn_tracer_config, only: glbsum                     !glbsum
      use gfs_dyn_physcons, fv => con_fvirt, rerth => con_rerth,
     &              grav => con_g,  cp => con_cp , rd => con_rd
      implicit none
!!
      real(kind=kind_grid), dimension(lonf,lats_node_a) :: zsg,psg
      real(kind=kind_grid), dimension(lonf,lats_node_a,levs) :: uug
     &,                              vvg, ttg
!    &,                              vvg, ttg, dpg
      real(kind=kind_grid) rqg(lonf,lats_node_a,levh)
!
      real(kind=kind_evod), dimension(len_trie_ls,2) :: gze, qe
      real(kind=kind_evod), dimension(len_trio_ls,2) :: gzo, qo
      real(kind=kind_evod), dimension(len_trie_ls,2,levs) :: tee,die,zee
      real(kind=kind_evod), dimension(len_trio_ls,2,levs) :: teo,dio,zeo
      real(kind=kind_evod), dimension(len_trie_ls,2,levs,ntrac) :: rqe
      real(kind=kind_evod), dimension(len_trio_ls,2,levs,ntrac) :: rqo
!
!     REAL(KIND=KIND_GRID) pwat   (lonf,lats_node_a)
!     REAL(KIND=KIND_GRID) ptot   (lonf,lats_node_a)
!     REAL(KIND=KIND_GRID) ptrc   (lonf,lats_node_a,ntrac)        !glbsum
!     REAL(KIND=KIND_GRID) work   (lonf)
!     REAL(KIND=KIND_GRID) tki    (lonf,levp1)
!     REAL(KIND=KIND_GRID) prsi   (lonf,levp1)

!     real(kind=kind_evod)  tkrt0
!     real(kind=kind_evod), parameter :: rkappa = cp / rd
!
      real(kind=kind_evod), allocatable  :: trie_ls(:,:,:),
     &                                      trio_ls(:,:,:)
       real(kind=kind_evod), allocatable :: for_gr_a_1(:,:,:),
     &                                      for_gr_a_2(:,:,:)
!!
      integer              ls_node(ls_dim,3), ls_nodes(ls_dim,nodes)
      integer, dimension(nodes) :: max_ls_nodes, lats_nodes_a
      integer, dimension(latg)  :: global_lats_a, lonsperlat
!
      real(kind=kind_evod) epse(len_trie_ls), epso(len_trio_ls)
     &,                    plnew_a(len_trie_ls,latg2)
     &,                    plnow_a(len_trio_ls,latg2)
!
!     real(kind=kind_evod) plnev_a(len_trie_ls,latg2)
!     real(kind=kind_evod) plnod_a(len_trio_ls,latg2)
!
!     real(kind=kind_evod) tfac(lonf,levs), sumq(lonf,levs), rcs2

      real(kind=kind_evod) rcs2, rcs, ga2, tem
!
      integer  i,j,k, nn, nnl,l,lan,lat,lotx,ylan,jtem,ktem
     &,        lon_dim,lons_lat,dimg,kk,lot_loc
     &,        locl,n,indev,indod,indev1,indev2,indod1,indod2
     &,        INDLSEV,JBASEV,INDLSOD,JBASOD
     &,        kau, kav, kat, kap, kaz

!     integer              lotdim,lotx
!

      real(kind=kind_evod), parameter :: one=1.0, pa2cb=0.001
!
!timers______________________________________________________---
!     real*8 rtc ,timer1,timer2
!timers______________________________________________________---
!
!     integer ,allocatable, dimension(:) :: lats_nodes_h,global_lats_h
!
      real(kind=kind_evod), parameter :: cons_0=0.0,   cons_24=24.0
     &,                                  cons_99=99.0, cons_1p0d9=1.0E9
     &,                                  qmin=1.0e-10
!
      INCLUDE 'function2'
!
!--------------------------------------------------------------------
!
!      write(0,*)' in  grid_to_spect_inp_slg '

       kau = 1
       kav = 1 + levs
       kat = 1 + levs + levs
       kap = 1 + levs + levs + levs
       kaz = 1 + kap

!      if (gg_tracers) then
!       if (.not. allocated(lats_nodes_h))
!    &                    allocate (lats_nodes_h(nodes))
!       if (.not. allocated(lat1s_h))
!    &                    allocate (lat1s_h(0:jcap))
!       if (.not. allocated(global_lats_h))
!    &                    allocate (global_lats_h(latg+2*yhalo*nodes))

!       call getcon_lag(lats_nodes_a,global_lats_a,
!    &                  lats_nodes_h,global_lats_h,
!    &                  lonsperlat,xhalo,yhalo)
!     endif

      if (gg_tracers) then
        lot_loc = kaz
      else
        lot_loc = kaz + levh
      endif
      allocate(trie_ls(len_trie_ls,2,lot_loc))
      allocate(trio_ls(len_trio_ls,2,lot_loc))

!     write(0,*)' before lan loop kau,kav,kat,kap,kaz=',
!    &            kau,kav,kat,kap,kaz

      allocate(for_gr_a_1(lon_dim_a,lot_loc,lats_dim_a))
      allocate(for_gr_a_2(lonf,     lot_loc,lats_dim_a))

!     write(0,*)' size of rgt_h=',size(rgt_h,dim=3)

      do lan=1,lats_node_a
        lon_dim  = lon_dims_a(lan)
        lat      = global_lats_a(ipt_lats_node_a-1+lan)
!       lat      = global_lats_a(ipt_lats_node_a+lats_node_a-lan)
        lons_lat = lonsperlat(lat)
        rcs2     = rcs2_a(min(lat,latg-lat+1))
        tem      = rcs2 * coslat_a(lat)
!
!      write(0,*)' in grid_to_spect_inp_slg coslat=',coslat_a(lat)
!    &,' rcs2=',rcs2,sqrt(rcs2),' lan=',lan,' tem=',tem,' me=',me
!
        do k=1,levs
          do i=1,lons_lat
            for_gr_a_2(i,kau+k-1,lan) = uug(i,lan,k) * tem
            for_gr_a_2(i,kav+k-1,lan) = vvg(i,lan,k) * tem
            for_gr_a_2(i,kat+k-1,lan) = ttg(i,lan,k)
     &                            * (one + fv*max(rqg(i,lan,k),qmin))
          enddo
        enddo
!
!       write(0,*)' me=',me,' psg=',psg(lons_lat,lan),' zsg=',
!    &zsg(lons_lat,lan),' lan=',lan,' lat=',lat
        do i=1,lons_lat
          for_gr_a_2(i,kap,lan) = log(psg(i,lan)*pa2cb)
          for_gr_a_2(i,kaz,lan) = zsg(i,lan)
        enddo

!     if (lan < 5)
!    &write(0,*)' in grid_to_spect_inp_slg ttg=',ttg(1,lan,:)
!    &,' lan=',lan,' ktd=',ktd

!       do k=1,levh
!         do i=1,lons_lat
!           for_gr_a_2(i,ktd+k,lan) = rqg(i,lan,k)
!         enddo
!       enddo
!       do i=1,lons_lat
!         ptot(i,lan) = psg(i,lan) * pa2cb
!       enddo
!       if (gen_coord_hybrid) then   ! Ps is the prognostic variable
!         do i=1,lons_lat
!           psg(i,lan) = psg(i,lan) * pa2cb
!         enddo
!       else                         ! ln(Ps) is the prognostic variable
!         do i=1,lons_lat
!           psg(i,lan) = log(psg(i,lan)*pa2cb)
!         enddo
!       endif
!
! get pressure at interfaces for pwat 
!       if (gen_coord_hybrid) then  
!         tki = 0.0
!         do k=2,levs
!           do i=1,lons_lat
!             tkrt0 = (ttg(i,lan,k-1)+ttg(i,lan,k))
!    &                           /(thref(k-1)+thref(k))
!             tki (i,k)=ck5(k)*tkrt0**rkappa
!           enddo
!         enddo
!         do k=1,levp1
!           do i=1,lons_lat
!             prsi(i,k)  = ak5(k)+bk5(k)*psg(i,lan)+tki(i,k) 
!           enddo
!         enddo
!       else if( hybrid ) then
!         do k=1,levp1
!           kk=levp1+1-k
!           do i=1,lons_lat
!             prsi(i,k)  = ak5(kk)+bk5(kk)*ptot(i,lan)
!           enddo
!         enddo
!       else
!         do k=1,levp1
!           do i=1,lons_lat
!             prsi(i,k)  = si(k)*ptot(i,lan)
!           enddo
!         enddo
!       endif                      
!
! get pwat (total vertical integrated water)
!       do i=1,lons_lat
!         pwat(i,lan) = 0.0
!       enddo
!       do k=1,levs
!         do i=1,lons_lat
!           work(i) = 0.0
!         enddo
!         if( ncld.gt.0 ) then
!           do nn=ntcw,ntcw+ncld-1
!             nnl = (nn-1)*levs
!             do i=1,lons_lat
!               work(i) = work(i) + rqg(i,lan,nnl+k)
!             enddo
!           enddo
!         endif
!         do i=1,lons_lat
! use definition for dpg instead of read in to have more accurate
! definition by th coordinates
!          if (.not.slg_flag) then
!           dpg (i,lan,k) = prsi(i,k)-prsi(i,k+1)
!           pwat(i,lan) = pwat(i,lan) + dpg(i,lan,k)
!    &                                * (rqg(i,lan,k) + work(i))
!          endif
!         enddo
!       enddo
!       if (.not.slg_flag) then
!         do k=1,levs
!           do i=1,lons_lat
!             for_gr_a_2(i+(kadp+k-2)*lon_dim,lan) = dpg(i,lan,k)
!           enddo
!         enddo
!         if( me==0 ) then
!           print *,' dpg in grid_to_spect_inp ',(dpg(1,lan,k),k=1,levs)
!         endif
!       endif

!
! compute ptrc (tracer global sum)                               !glbsum
!
!       if( glbsum ) then                                        !glbsum
!         do nn = 1, ntrac                                       !glbsum
!           nnl = (nn-1)*levs                                    !glbsum
!           do i=1,lons_lat                                      !glbsum
!            ptrc(i,lan,nn) = 0.0                                !glbsum
!            do k=1,levs                                         !glbsum
!              ptrc(i,lan,nn) = ptrc(i,lan,nn) +                 !glbsum
!    &         (prsi(i,k)-prsi(i,k+1))*rqg(i,lan,nnl+k)          !glbsum
!            enddo                                               !glbsum
!           enddo                                                !glbsum
!         enddo                                                  !glbsum
!       endif                                                    !glbsum

        
!       write(0,*)' me=',me,' calling grid_to_four'
        call grid_to_four(for_gr_a_2(1,1,lan),for_gr_a_1(1,1,lan),
     &                    lon_dim_a-2,lon_dim_a,lons_lat,kaz)
!       write(0,*)' me=',me,' after calling grid_to_four'
!
        if (.not. gg_tracers) then
          do k=1,levh
            do i=1,lons_lat
              for_gr_a_2(i,kaz+k,lan) = rqg(i,lan,k)
            enddo
          enddo
          call grid_to_four(for_gr_a_2(1,kaz+1,lan),
     &                      for_gr_a_1(1,kaz+1,lan),
     &                      lon_dim_a-2,lon_dim_a,lons_lat,levh)
        else
          ylan = lats_node_a + 1 - lan + yhalo
!         ylan = yhalo + lan
          do n=1,ntrac
!$omp parallel do private(k,jtem,ktem,i)
            do k=1,levs
              jtem = levs + 1 - k
              ktem = n*levs - k + 1
              do i=1,min(lonf,lons_lat)
                rgt_h(xhalo+i,k,ylan,n) = rqg(i,lan,ktem)
              enddo
            enddo
          enddo
!         write(0,*)' in grid rqg=',rqg(1,45,1:levs)*1000
!         write(0,*)' in grid rgt_h=',
!    &rgt_h(xhalo+1,1:levs,lats_node_a+yhalo-44,1)*1000
        endif
      enddo                  ! end of lan loop
!     write(0,*)' in input rtg_hmaxmin=',
!    &      maxval(rgt_h(xhalo+1:xhalo+30,1:64,yhalo+1:yhalo+30,1))
!    &,     minval(rgt_h(xhalo+1:xhalo+30,1:64,yhalo+1:yhalo+30,1))
!
!     write(0,*)' p_ze=',p_ze,' p_di=',p_di,' p_uln=',p_uln
!    &,' p_vln=',p_vln
!     write(0,*)' before calling four2fln_gg'

!     call four2fln_gg(lats_dim_a,lota,3*levs+1,for_gr_a_1,
      call four2fln_gg(lats_dim_a,lot_loc,kaz,for_gr_a_1,
     &                 ls_nodes,max_ls_nodes,
     &                 lats_nodes_a,global_lats_a,lon_dim_a,
     &                 lats_node_a,ipt_lats_node_a,
     &                 lat1s_a,lon_dim_a,latg,latg2,
     &                 trie_ls(1,1,1), trio_ls(1,1,1),
     &                 plnew_a, plnow_a,
     &                 ls_node,0,
     &                 kat,kaz)

      if (.not. gg_tracers) then
        call four2fln_gg(lats_dim_a,lot_loc,levh,for_gr_a_1,
     &                   ls_nodes,max_ls_nodes,
     &                   lats_nodes_a,global_lats_a,lon_dim_a,
     &                   lats_node_a,ipt_lats_node_a,
     &                   lat1s_a,lon_dim_a,latg,latg2,
     &                   trie_ls(1,1,kaz+1), trio_ls(1,1,kaz+1),
     &                   plnew_a, plnow_a,
     &                   ls_node,kaz,1,levh)

        do n=1,ntrac
!$omp parallel do private(i,k,kk)
          do k=1,levs
            kk = kaz + (n-1)*levs + k
            do i=1,len_trie_ls
              rqe(i,1,k,n) = trie_ls(i,1,kk)
              rqe(i,2,k,n) = trie_ls(i,2,kk)
            enddo
            do i=1,len_trio_ls
              rqo(i,1,k,n) = trio_ls(i,1,kk)
              rqo(i,2,k,n) = trio_ls(i,2,kk)
            enddo
          enddo
        enddo
      endif
!
      do i=1,len_trie_ls
        qe(i,1)  = trie_ls(i,1,kap)
        qe(i,2)  = trie_ls(i,2,kap)
        gze(i,1) = trie_ls(i,1,kaz)
        gze(i,2) = trie_ls(i,2,kaz)
      enddo
      do i=1,len_trio_ls
        qo(i,1)  = trio_ls(i,1,kap)
        qo(i,2)  = trio_ls(i,2,kap)
        gzo(i,1) = trio_ls(i,1,kaz)
        gzo(i,2) = trio_ls(i,2,kaz)
      enddo
!     write(0,*)' qe1=',qe(1,1),' gze1=',gze(1,1)
!
!$omp parallel do private(i,k,kk)
      do k=1,levs
        kk = kat + k - 1
        do i=1,len_trie_ls
          tee(i,1,k) = trie_ls(i,1,kk)
          tee(i,2,k) = trie_ls(i,2,kk)
          die(i,1,k) = 0.0
          die(i,2,k) = 0.0
          zee(i,1,k) = 0.0
          zee(i,2,k) = 0.0
        enddo
        do i=1,len_trio_ls
          teo(i,1,k) = trio_ls(i,1,kk)
          teo(i,2,k) = trio_ls(i,2,kk)
          dio(i,1,k) = 0.0
          dio(i,2,k) = 0.0
          zeo(i,1,k) = 0.0
          zeo(i,2,k) = 0.0
        enddo
!     write(0,*)' before k loop for uveodz tee=',tee(1,1,k),' lan=',lan
!    &,' k=',k
      enddo
!
!$omp parallel do private(k)
      do k=1,levs
         call uveodz(trie_ls(1,1,kau+k-1), trio_ls(1,1,kav+k-1),
     &               die(1,1,k),           zeo(1,1,k),
     &               epse,epso,ls_node)

         call uvoedz(trio_ls(1,1,kau+k-1), trie_ls(1,1,kav+k-1),
     &               dio(1,1,k),           zee(1,1,k),
     &               epse,epso,ls_node)

      enddo
      deallocate (trie_ls, trio_ls)
      deallocate (for_gr_a_1,for_gr_a_2)

! =======================================================================
!     do lan=1,lats_node_a
!
!        lon_dim = lon_dims_a(lan)
!
!        lat = global_lats_a(ipt_lats_node_a-1+lan)
!        lons_lat = lonsperlat(lat)

!        call grid2four_thread(for_gr_a_2(1,lan),for_gr_a_1(1,lan),
!    &                  lon_dim,lons_lat,lonfx,lotx)
!
!     enddo
!
!     dimg=0
!     call four2fln(lats_dim_a,lotdim,lotx,for_gr_a_1,
!    x              ls_nodes,max_ls_nodes,
!    x              lats_nodes_a,global_lats_a,lon_dims_a,
!    x              lats_node_a,ipt_lats_node_a,dimg,
!    x              lat1s_a,lonfx,latg,latg2,
!    x              trie_ls(1,1,p_w), trio_ls(1,1,p_w),
!    x              plnew_a, plnow_a,
!    x              ls_node,2*levs)
!
!!$OMP parallel do shared(trie_ls,trio_ls)
!!$OMP+shared(p_w,p_x,p_uln,p_vln,epse,epso,ls_node)
!!$OMP+private(k)
!     do k=1,levs
!        call uveodz(trie_ls(1,1,P_w  +k-1), trio_ls(1,1,P_x  +k-1),
!    x               trie_ls(1,1,P_uln+k-1), trio_ls(1,1,P_vln+k-1),
!    x               epse,epso,ls_node)
!
!        call uvoedz(trio_ls(1,1,P_w  +k-1), trie_ls(1,1,P_x  +k-1),
!    x               trio_ls(1,1,P_uln+k-1), trie_ls(1,1,P_vln+k-1),
!    x               epse,epso,ls_node)
!     enddo
!
!   move uln back to x
!   move vln back to w
!
!     do k=1,levs
!        do i=1,len_trie_ls
!           trie_ls(i,1,P_x +k-1)= trie_ls(i,1,P_uln +k-1)
!           trie_ls(i,2,P_x +k-1)= trie_ls(i,2,P_uln +k-1)
!           trie_ls(i,1,P_w +k-1)= trie_ls(i,1,P_vln +k-1)
!           trie_ls(i,2,P_w +k-1)= trie_ls(i,2,P_vln +k-1)
!        enddo
!        do i=1,len_trio_ls
!           trio_ls(i,1,P_x +k-1)= trio_ls(i,1,P_uln +k-1)
!           trio_ls(i,2,P_x +k-1)= trio_ls(i,2,P_uln +k-1)
!           trio_ls(i,1,P_w +k-1)= trio_ls(i,1,P_vln +k-1)
!           trio_ls(i,2,P_w +k-1)= trio_ls(i,2,P_vln +k-1)
!        enddo
!     enddo
!
! -------------------------------------------------------------------
! model realted filter such as reduced grid spectral transform for zs
!
!     if( fhour.eq.0.0 ) then

!     dimg=0
!
!     call sumflna(trie_ls(1,1,p_gz),trio_ls(1,1,p_gz),
!            lat1s_a,
!            plnev_a,plnod_a,
!            1,ls_node,latg2,
!            lats_dim_a,lotdim,for_gr_a_1,
!            ls_nodes,max_ls_nodes,
!            lats_nodes_a,global_lats_a,
!            lats_node_a,ipt_lats_node_a,lon_dims_a,dimg,
!            lonsperlat,lonfx,latg)

!     do lan=1,lats_node_a
!
!        lat = global_lats_a(ipt_lats_node_a-1+lan)
!        lon_dim = lon_dims_a(lan)
!        lons_lat = lonsperlat(lat)
!c
!        call four2grid_thread(for_gr_a_1(1,lan),for_gr_a_2(1,lan),
!    &                  lon_dim,lons_lat,lonfx,1,lan,me)
 
!        do i=1,lons_lat
!          zsg(i,lan) = for_gr_a_2(i,lan)
!        enddo
!     enddo   !lan
!     endif	! fhour=0.0
! -------------------------------------------------------------------
!     write(0,*)' exit grid_to_spect_inp_slg '
!!
      return
      end
