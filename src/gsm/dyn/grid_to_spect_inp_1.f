      subroutine grid_to_spect_inp_1
     &    (psg,uug,vvg,ttg,rqg,
     &     trie_ps,trio_ps,
     &     trie_di,trio_di,trie_ze,trio_ze,
     &     trie_te,trio_te,trie_rq,trio_rq,
     &     ls_node,ls_nodes,max_ls_nodes,
     &     lats_nodes_a,global_lats_a,lonsperlat,
     &     epse,epso, plnew_a,plnow_a,pwat,ptot,ptrc)
!!
!! hmhj - this routine do spectral to grid transform 
!!        from gfsio read in field, to model fields
!! input psg,uug,vvg,ttg,rqg (mapping wind, temp)
!! output psg,uug,vvg,ttg,rqg in model values (mapping wind, enthalpy)
!! aug 2010      sarah lu, modified to compute tracer global sum
!!
      use gfs_dyn_resol_def
      use gfs_dyn_layout1
      use gfs_dyn_gg_def
      use gfs_dyn_vert_def
      use gfs_dyn_date_def
      use namelist_dynamics_def
      use gfs_dyn_coordinate_def 
      use gfs_dyn_tracer_const
      use gfs_dyn_tracer_config, only: glbsum                     !glbsum
      use gfs_dyn_physcons, fv => con_fvirt, rerth => con_rerth,
     &              grav => con_g,  cp => con_cp , rd => con_rd
      implicit none
!!
      real(kind=kind_evod) psg(lonf,lats_node_a)
      real(kind=kind_evod) uug(lonf,lats_node_a,levs)
      real(kind=kind_evod) vvg(lonf,lats_node_a,levs)
      real(kind=kind_evod) ttg(lonf,lats_node_a,levs)
      real(kind=kind_evod) rqg(lonf,lats_node_a,levh)
!
      REAL(KIND=KIND_GRID) pwat   (lonf,lats_node_a)
      REAL(KIND=KIND_GRID) ptot   (lonf,lats_node_a)
      REAL(KIND=KIND_GRID) ptrc   (lonf,lats_node_a,ntrac)        !glbsum
      REAL(KIND=KIND_GRID) work   (lonf)
      REAL(KIND=KIND_GRID) tki    (lonf,levp1)
      REAL(KIND=KIND_GRID) prsi   (lonf,levp1)

      real(kind=kind_evod)  tkrt0
      real(kind=kind_evod), parameter :: rkappa = cp / rd
!
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
!!!!  integer, parameter :: lota = 3*levs+1*levh+1 
!
      real(kind=kind_evod) trie_ls(len_trie_ls,2,lota)
      real(kind=kind_evod) trio_ls(len_trio_ls,2,lota)
!!
      real(kind=kind_evod) for_gr_a_1(lonfx*(lota),lats_dim_a)
      real(kind=kind_evod) for_gr_a_2(lonfx*(lota),lats_dim_a)
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
      real(kind=kind_evod)   plnew_a(len_trie_ls,latg2)
      real(kind=kind_evod)   plnow_a(len_trio_ls,latg2)
!
      real(kind=kind_evod)   tfac(lonf,levs), sumq(lonf,levs), rcs2
!
      integer              i,j,k,kk, nn, nnl
      integer              l,lan,lat
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
      logical 	lslag
      logical , parameter :: repro = .false.
!

      real(kind=kind_evod), parameter :: one=1.0D0, pa2cb=0.001D0
!
!timers______________________________________________________---
      real*8 rtc ,timer1,timer2
!timers______________________________________________________---
!
!
      real(kind=kind_evod), parameter :: cons_0=0.0D0, cons_24=24.0D0
     &,                                cons_99=99.0D0, cons_1p0d9=1.0D9
     &,                                  qmin=1.0D-10
!
      real(kind=kind_evod) ga2, tem
!
      INCLUDE 'function2'

!
!--------------------------------------------------------------------
!
      lslag   = .false.
!
      trie_ls = 0.0D0
      trio_ls = 0.0D0
!
!--------------------------------------------------------------------
      do lan=1,lats_node_a
        lon_dim = lon_dims_a(lan)
        lat = global_lats_a(ipt_lats_node_a-1+lan)
        lons_lat = lonsperlat(lat)
        rcs2     = rcs2_a(min(lat,latg-lat+1))
!
        if (thermodyn_id == 3) then
          do k=1,levs
            do i=1,lons_lat
              tfac(i,k) = 0.0D0
              sumq(i,k) = 0.0D0
            enddo
          enddo
          do nn=1,ntrac
            nnl = (nn-1)*levs
            if (cpi(nn) .ne. 0.0) then
              do k=1,levs
                do i=1,lons_lat
                  sumq(i,k) = sumq(i,k) + rqg(i,lan,nnl+k)
                  tfac(i,k) = tfac(i,k) + cpi(nn)*rqg(i,lan,nnl+k)
                enddo
              enddo
            endif
          enddo
          do k=1,levs
            do i=1,lons_lat
              tfac(i,k) = (one-sumq(i,k))*cpi(0) + tfac(i,k)
            enddo
          enddo
        else
          do k=1,levs
            do i=1,lons_lat
              tfac(i,k) = one + fv*max(rqg(i,lan,k),qmin) 
            enddo
          enddo
        endif

        do k=1,levs
          do i=1,lons_lat
            uug(i,lan,k) = uug(i,lan,k) * coslat_a(lat)
            vvg(i,lan,k) = vvg(i,lan,k) * coslat_a(lat)
            ttg(i,lan,k) = ttg(i,lan,k) * tfac(i,k)
            for_gr_a_2(i+(kat+k-2)*lon_dim,lan) = ttg(i,lan,k)
            for_gr_a_2(i+(kau+k-2)*lon_dim,lan) = uug(i,lan,k) * rcs2
            for_gr_a_2(i+(kav+k-2)*lon_dim,lan) = vvg(i,lan,k) * rcs2
          enddo
        enddo
        do k=1,levh
          do i=1,lons_lat
            for_gr_a_2(i+(kar+k-2)*lon_dim,lan)=rqg(i,lan,k)
          enddo
        enddo
        do i=1,lons_lat
          ptot(i,lan) = psg(i,lan) * pa2cb
        enddo
        if (gen_coord_hybrid) then   ! Ps is the prognostic variable
          do i=1,lons_lat
            psg(i,lan) = psg(i,lan) * pa2cb
          enddo
        else                         ! ln(Ps) is the prognostic variable
          do i=1,lons_lat
            psg(i,lan) = log(psg(i,lan)*pa2cb)
          enddo
        endif
        do i=1,lons_lat
          for_gr_a_2(i+(kaps-1)*lon_dim,lan) = psg(i,lan)
        enddo
!
! get pressure at interfaces for pwat 
        if (gen_coord_hybrid) then  
          tki = 0.0D0
          do k=2,levs
            do i=1,lons_lat
              tkrt0 = (ttg(i,lan,k-1)+ttg(i,lan,k))
     &                           /(thref(k-1)+thref(k))
              tki (i,k)=ck5(k)*tkrt0**rkappa
            enddo
          enddo
          do k=1,levp1
            do i=1,lons_lat
              prsi(i,k)  = ak5(k)+bk5(k)*psg(i,lan)+tki(i,k) 
            enddo
          enddo
        else if (hybrid) then
          do k=1,levp1
            kk=levp1+1-k
            do i=1,lons_lat
              prsi(i,k)  = ak5(kk)+bk5(kk)*ptot(i,lan)
            enddo
          enddo
        else
          do k=1,levp1
            do i=1,lons_lat
              prsi(i,k)  = si(k)*ptot(i,lan)
            enddo
          enddo
        endif                      
!
! get pwat (total vertical integrated water)
        do i=1,lons_lat
          pwat(i,lan) = 0.0D0
        enddo
        do k=1,levs
          do i=1,lons_lat
            work(i) = 0.0D0
          enddo
          if( ncld.gt.0 ) then
            do nn=ntcw,ntcw+ncld-1
              nnl = (nn-1)*levs
              do i=1,lons_lat
                work(i) = work(i) + rqg(i,lan,nnl+k)
              enddo
            enddo
          endif
          do i=1,lons_lat
            pwat(i,lan) = pwat(i,lan) + (prsi(i,k)-prsi(i,k+1))
     &                                * (rqg(i,lan,k) + work(i))
          enddo
        enddo
!
! get ptrc (tracer global sum)                                   !glbsum
!
        if( glbsum ) then                                        !glbsum
          do nn = 1, ntrac                                       !glbsum
            nnl = (nn-1)*levs                                    !glbsum
            do i=1,lons_lat                                      !glbsum
             ptrc(i,lan,nn) = 0.0D0                              !glbsum
             do k=1,levs                                         !glbsum
               ptrc(i,lan,nn) = ptrc(i,lan,nn) +                 !glbsum
     &         (prsi(i,k)-prsi(i,k+1))*rqg(i,lan,nnl+k)          !glbsum
             enddo                                               !glbsum
            enddo                                                !glbsum
          enddo                                                  !glbsum
        endif                                                    !glbsum

!
      enddo
!
! =======================================================================
      do lan=1,lats_node_a
!
         lon_dim = lon_dims_a(lan)
!
         lat = global_lats_a(ipt_lats_node_a-1+lan)
         lons_lat = lonsperlat(lat)

         call grid2four_thread(for_gr_a_2(1,lan),for_gr_a_1(1,lan),
     &                  lon_dim,lons_lat,lonfx,lota)
!
      enddo
!
      dimg=0
      call four2fln(lslag,lats_dim_a,lota,lota,for_gr_a_1,
     x              ls_nodes,max_ls_nodes,
     x              lats_nodes_a,global_lats_a,lon_dims_a,
     x              lats_node_a,ipt_lats_node_a,dimg,
     x              lat1s_a,lonfx,latg,latg2,
     x              trie_ls(1,1,1), trio_ls(1,1,1),
     x              plnew_a, plnow_a,
     x              ls_node)
!
!
      trie_di = 0.0D0
      trio_di = 0.0D0
      trie_ze = 0.0D0
      trio_ze = 0.0D0
!
!$omp parallel do shared(trie_ls,trio_ls)
!$omp+shared(trie_di,trio_di,trie_ze,trio_ze,trie_te,trio_te)
!$omp+shared(kau,kav,kat,epse,epso,ls_node)
!$omp+private(k)
      do k=1,levs
         call uveodz(trie_ls(1,1,kau+k-1), trio_ls(1,1,kav+k-1),
     x               trie_di(1,1,k),       trio_ze(1,1,k),
     x               epse,epso,ls_node)
!
         call uvoedz(trio_ls(1,1,kau+k-1), trie_ls(1,1,kav+k-1),
     x               trio_di(1,1,k),       trie_ze(1,1,k),
     x               epse,epso,ls_node)
        trie_te(:,:,k)=trie_ls(:,:,kat+k-1)
        trio_te(:,:,k)=trio_ls(:,:,kat+k-1)
      enddo
      do k=1,levh
        trie_rq(:,:,k)=trie_ls(:,:,kar+k-1)
        trio_rq(:,:,k)=trio_ls(:,:,kar+k-1)
      enddo
      trie_ps(:,:)=trie_ls(:,:,kaps)
      trio_ps(:,:)=trio_ls(:,:,kaps)

      END SUBROUTINE grid_to_spect_inp_1
