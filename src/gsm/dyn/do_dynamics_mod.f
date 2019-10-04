      module do_dynamics_mod
!
! group of subroutine to take care model grid to temporary grid etc 
! reassignment. It is coded initially by hann-ming henry juang
!
! program log:
! 2011 02 20 :    Henry Juang, add NDSL and MASS_DP options
!                 for semi-Lagrangian advection
! 2012 11 26 :    Shrinivas Moorthi - Rewrote the routines
!                 do_dynamics_gridomega and do_dynamics_gridupdate 
!                 to make threading work at higher horizontal resolution
! 2013 January:   Henry Juang, add digital filter with
!                 spectral coefficient routine, _spectdfini
! 2013 June   :   Henry Juang, add preparation of moisture for transform in ndsl
! 2013 November:  Weiyu Yang, add subroutines for the slg code.
! 2014 June    : Shrinivas Moorthi - add additional code for Semi-Lagrangian (Sela)
!
      use gfs_dyn_resol_def
      use gfs_dyn_layout1
      use namelist_dynamics_def
      use gfs_dyn_mpi_def
      implicit none
      private
!
      public do_dynamics_gridaddlapgz
      public do_dynamics_gridc2syq
      public do_dynamics_gridc2syn
      public do_dynamics_advhn2anl
      public do_dynamics_gridn2anl
      public do_dynamics_gridt2anl
      public do_dynamics_gridt2rqt
      public do_dynamics_gridm2rqt
      public do_dynamics_gridc2rqt
      public do_dynamics_spectaddlapgz
      public do_dynamics_spectupdatewrt
      public do_dynamics_spectupdatexydpnzq
      public do_dynamics_spectn2c
      public do_dynamics_spectn2m
      public do_dynamics_spectc2n
      public do_dynamics_spectdfini
      public do_dynamics_spectdfini_slg
      public do_dynamics_syn2gridn
      public do_dynamics_syn2gridn_slg
      public do_dynamics_gridomega_slg
      public do_dynamics_gridn2anl_slg
      public do_dynamics_gridn2anl_slg_dfi
      public do_dynamics_gridomega
      public do_dynamics_gridfilter
      public do_dynamics_gridc2n
      public do_dynamics_gridn2c
      public do_dynamics_gridn2m
      public do_dynamics_gridn2p
      public do_dynamics_gridm2p
      public do_dynamics_gridp2n
      public do_dynamics_gridap2n
      public do_dynamics_gridupdate
      public do_dynamics_gridpdpn
      public do_dynamics_gridpdp
      public do_dynamics_gridp
      public do_dynamics_gridzz
      public do_dynamics_gridzk
      public do_dynamics_gridspdmax
      public do_dynamics_gridmean
      public do_dynamics_gridcheck

      contains


! --------------------------------------------------------------
      subroutine  do_dynamics_gridaddlapgz(anl,syn,
     &                                     global_lats_a,lonsperlat)

      use gfs_dyn_physcons, grav  => con_g

      real(kind=kind_grid) anl(lonfx*lota,lats_dim_ext)
      real(kind=kind_evod) syn(lonfx*lots,lats_dim_ext)
      integer,intent(in):: global_lats_a(latg)
      integer,intent(in):: lonsperlat(latg)

      integer   lan,lat,lon_dim,lons_lat,k,i
      integer   jlonf,ilan

!$omp parallel do private(lan)
!$omp+private(lat,lon_dim,lons_lat,i,k,ilan)

        do lan=1,lats_node_a
          lat = global_lats_a(ipt_lats_node_a-1+lan)
          lon_dim = lon_dims_a(lan)
          lons_lat = lonsperlat(lat)
!
          do k=1,levr
            do i=1,lons_lat
             anl(i+(kau+k-2)*lon_dim,lan)=anl(i+(kau+k-2)*lon_dim,lan)-
     &                 grav*syn(i+(kzslam-1)*lon_dim,lan)
             anl(i+(kav+k-2)*lon_dim,lan)=anl(i+(kav+k-2)*lon_dim,lan)-
     &                 grav*syn(i+(kzsphi-1)*lon_dim,lan)
            enddo
          enddo
!
        enddo

      return
      end subroutine do_dynamics_gridaddlapgz

! --------------------------------------------------------------
      subroutine  do_dynamics_gridc2syq(grid_gr,syn_gr_syq,
     &                                  global_lats_a,lonsperlat)

      real(kind=kind_grid) grid_gr(lonf*lats_node_a_max,lotgr)
      real(kind=kind_evod) syn_gr_syq(lonfx*levh,lats_dim_ext)
      integer,intent(in):: global_lats_a(latg)
      integer,intent(in):: lonsperlat(latg)
 
      integer	lan,lat,lon_dim,lons_lat,k,i
      integer   jlonf,ilan

!$omp parallel do private(lan)
!$omp+private(lat,lon_dim,lons_lat,jlonf,i,k,ilan)

        do lan=1,lats_node_a
          lat = global_lats_a(ipt_lats_node_a-1+lan)
          lon_dim = lon_dims_a(lan)
          lons_lat = lonsperlat(lat)
          jlonf = (lan-1)*lonf
!
          do k=1,levh
            do i=1,lons_lat
             ilan=i+jlonf
             syn_gr_syq(i+(k-1)*lon_dim,lan)=grid_gr(ilan,g_rq+k-1)
            enddo
          enddo
!
        enddo

      return
      end subroutine do_dynamics_gridc2syq

! --------------------------------------------------------------
      subroutine  do_dynamics_gridc2syn(grid_gr,syn_gr_a_2,
     &                                  global_lats_a,lonsperlat)

      real(kind=kind_grid) grid_gr(lonf*lats_node_a_max,lotgr)
      real(kind=kind_evod) syn_gr_a_2(lonfx*lots,lats_dim_ext)
      integer,intent(in):: global_lats_a(latg)
      integer,intent(in):: lonsperlat(latg)
 
      integer	lan,lat,lon_dim,lons_lat,k,i
      integer   jlonf,ilan

!$omp parallel do private(lan)
!$omp+private(lat,lon_dim,lons_lat,jlonf,i,k,ilan)

        do lan=1,lats_node_a
          lat = global_lats_a(ipt_lats_node_a-1+lan)
          lon_dim = lon_dims_a(lan)
          lons_lat = lonsperlat(lat)
          jlonf = (lan-1)*lonf
          do k=1,levs
            do i=1,lons_lat
             ilan=i+jlonf
             syn_gr_a_2(i+(ksu -2+k)*lon_dim,lan)=grid_gr(ilan,g_uu+k-1)
             syn_gr_a_2(i+(ksv -2+k)*lon_dim,lan)=grid_gr(ilan,g_vv+k-1)
             syn_gr_a_2(i+(kst -2+k)*lon_dim,lan)=grid_gr(ilan,g_tt+k-1)
             syn_gr_a_2(i+(ksdp-2+k)*lon_dim,lan)=grid_gr(ilan,g_dp+k-1)
            enddo
          enddo
          do i=1,lons_lat
             ilan=i+jlonf
             syn_gr_a_2(i+(ksq-1)*lon_dim,lan)=grid_gr(ilan,g_q)
          enddo
!
          if( .not. ndslfv ) then
          do k=1,levh
            do i=1,lons_lat
             ilan=i+jlonf
             syn_gr_a_2(i+(ksr-2+k)*lon_dim,lan)=grid_gr(ilan,g_rq+k-1)
            enddo
          enddo
          endif
!
        enddo

      return
      end subroutine do_dynamics_gridc2syn
!
! --------------------------------------------------------------
      subroutine  do_dynamics_advhn2anl(grid_gr,anl_gr_a_2,
     &                                  global_lats_a,lonsperlat)

      real(kind=kind_grid) grid_gr(lonf*lats_node_a_max,lotgr)
      real(kind=kind_evod) anl_gr_a_2(lonfx*lota,lats_dim_ext)
      integer,intent(in):: global_lats_a(latg)
      integer,intent(in):: lonsperlat(latg)
 
      integer	lan,lat,lon_dim,lons_lat,k,i
      integer   jlonf,ilan

!$omp parallel do private(lan)
!$omp+private(lat,lon_dim,lons_lat,jlonf,i,k,ilan)

      do lan=1,lats_node_a
        lat = global_lats_a(ipt_lats_node_a-1+lan)
        lon_dim = lon_dims_a(lan)
        lons_lat = lonsperlat(lat)
        jlonf = (lan-1)*lonf
        do k=1,levs
          do i=1,lons_lat
            ilan=i+jlonf
            anl_gr_a_2(i+(kau -2+k)*lon_dim,lan)=grid_gr(ilan,g_u  +k-1)
            anl_gr_a_2(i+(kav -2+k)*lon_dim,lan)=grid_gr(ilan,g_v  +k-1)
            anl_gr_a_2(i+(kat -2+k)*lon_dim,lan)=grid_gr(ilan,g_t  +k-1)
            anl_gr_a_2(i+(kadp-2+k)*lon_dim,lan)=grid_gr(ilan,g_dpn+k-1)
            anl_gr_a_2(i+(kap2-2+k)*lon_dim,lan)=grid_gr(ilan,g_p  +k-1)
          enddo
        enddo
        do i=1,lons_lat
            ilan=i+jlonf
            anl_gr_a_2(i+(kaps-1)*lon_dim,lan)=grid_gr(ilan,g_zq)
        enddo
!
        if( .not. ndslfv ) then
        do k=1,levh
          do i=1,lons_lat
            ilan=i+jlonf
            anl_gr_a_2(i+(kar-2+k)*lon_dim,lan)=grid_gr(ilan,g_rt +k-1)
          enddo
        enddo
        endif
!
      enddo

      return
      end subroutine do_dynamics_advhn2anl
!
! --------------------------------------------------------------
      subroutine  do_dynamics_gridn2anl(grid_gr,anl_gr_a_2,
     &                                  global_lats_a,lonsperlat)

      real(kind=kind_grid) grid_gr(lonf*lats_node_a_max,lotgr)
      real(kind=kind_evod) anl_gr_a_2(lonfx*lota,lats_dim_ext)
      integer,intent(in):: global_lats_a(latg)
      integer,intent(in):: lonsperlat(latg)
 
      integer	lan,lat,lon_dim,lons_lat,k,i
      integer   jlonf,ilan
!$omp parallel do private(lan)
!$omp+private(lat,lon_dim,lons_lat,jlonf,i,k,ilan)

      do lan=1,lats_node_a
        lat = global_lats_a(ipt_lats_node_a-1+lan)
        lon_dim = lon_dims_a(lan)
        lons_lat = lonsperlat(lat)
        jlonf = (lan-1)*lonf
        do k=1,levs
          do i=1,lons_lat
            ilan=i+jlonf
            anl_gr_a_2(i+(kau -2+k)*lon_dim,lan)=grid_gr(ilan,g_u  +k-1)
            anl_gr_a_2(i+(kav -2+k)*lon_dim,lan)=grid_gr(ilan,g_v  +k-1)
            anl_gr_a_2(i+(kat -2+k)*lon_dim,lan)=grid_gr(ilan,g_t  +k-1)
            anl_gr_a_2(i+(kadp-2+k)*lon_dim,lan)=grid_gr(ilan,g_dpn+k-1)
          enddo
        enddo
        do i=1,lons_lat
            ilan=i+jlonf
            anl_gr_a_2(i+(kaps-1)*lon_dim,lan)=grid_gr(ilan,g_zq)
        enddo
!
        if( .not. ndslfv ) then
        do k=1,levh
          do i=1,lons_lat
            ilan=i+jlonf
            anl_gr_a_2(i+(kar-2+k)*lon_dim,lan)=grid_gr(ilan,g_rt +k-1)
          enddo
        enddo
        endif
!
      enddo
      return
      end subroutine do_dynamics_gridn2anl
!
! --------------------------------------------------------------
      subroutine do_dynamics_gridn2anl_slg(grid_gr,anl_gr_a_2,
     &                                     rcs2,
     &                                     global_lats_a,lonsperlat
     &                                     ,iniauinterval)

      real(kind=kind_grid) grid_gr(lonf*lats_node_a_max,lotgr)
!     real(kind=kind_evod) anl_gr_a_2(lonfx*lota,lats_dim_ext)
      real(kind=kind_evod) anl_gr_a_2(lonf,lota,lats_dim_a)
      real(kind=kind_grid) rcs2(latg2)
      integer,intent(in):: global_lats_a(latg)
      integer,intent(in):: lonsperlat(latg)
      logical,intent(in) :: iniauinterval
      integer   jlonf, ilan, kap, kau, kav, kat, kar
     &,         lan,lat,lon_dim,lons_lat,k,i,kk
      real      rcs_loc

      kau = 0*levs + 0*levh + 1
      kav = 1*levs + 0*levh + 1
      kat = 2*levs + 0*levh + 1
      kap = 3*levs + 0*levh + 1
      kar = 3*levs + 0*levh + 2
 
!$omp parallel do private(lan,lat,lon_dim,lons_lat)
!$omp+private(jlonf,i,k,kk,ilan,rcs_loc)

      do lan=1,lats_node_a
        lat      = global_lats_a(ipt_lats_node_a-1+lan)
        lon_dim  = lon_dims_a(lan)
        lons_lat = lonsperlat(lat)
        jlonf    = (lan-1)*lonf
        rcs_loc  = sqrt(rcs2(min(lat,latg-lat+1)))

        do k=1,levs
          do i=1,lons_lat
            ilan = i+jlonf
            anl_gr_a_2(i,kau-1+k,lan) = grid_gr(ilan,g_u+k-1) * rcs_loc
            anl_gr_a_2(i,kav-1+k,lan) = grid_gr(ilan,g_v+k-1) * rcs_loc
            anl_gr_a_2(i,kat-1+k,lan) = grid_gr(ilan,g_t+k-1)
!           anl_gr_a_2(i,kadp-1+k,lan) = grid_gr(ilan,g_dpn+k-1)
          enddo
        enddo
!
        if (gg_tracers .and. .not. iniauinterval) then
          do i=1,lons_lat
            ilan = i+jlonf
            anl_gr_a_2(i,kap,lan) = grid_gr(ilan,g_rqtk)
          enddo
        else
          do i=1,lons_lat
!            anl_gr_a_2(i,kap,lan) = 0.0
            ilan = i + jlonf
            anl_gr_a_2(i,kap,lan) = grid_gr(ilan,g_zqp)  ! surface pressure with IAU increment
          enddo
        endif
!
        if( .not. ndslfv) then
          do k=1,levh
            kk = kar-1+k
            do i=1,lons_lat
              ilan = i + jlonf
              anl_gr_a_2(i,kk,lan) = grid_gr(ilan,g_rt +k-1)
            enddo
          enddo
        endif
!
      enddo

      return

      end subroutine do_dynamics_gridn2anl_slg
!
      subroutine do_dynamics_gridn2anl_slg_dfi(grid_gr,anl_gr_a_2, rcs2,
     &                                         global_lats_a,lonsperlat)

      real(kind=kind_grid) grid_gr(lonf*lats_node_a_max,lotgr)
!     real(kind=kind_evod) anl_gr_a_2(lonfx*lota,lats_dim_ext)
      real(kind=kind_evod) anl_gr_a_2(lonf,lota,lats_dim_a)
      real(kind=kind_grid) rcs2(latg2)
      integer, intent(in), dimension(latg) :: global_lats_a, lonsperlat
      integer   jlonf, ilan, kap, kau, kav, kat, kar
     &,         lan,lat,lon_dim,lons_lat,k,i,kk
      real      rcs_loc

      kau = 1
      kav = kau + levs
      kat = kav + levs
      kap = kat + levs
      kar = kap + 1
 
!$omp parallel do private(lan,lat,lon_dim,lons_lat,jlonf,i,k,kk,ilan)
!$omp+private(rcs_loc)

      do lan=1,lats_node_a
        lat = global_lats_a(ipt_lats_node_a-1+lan)
        lon_dim  = lon_dims_a(lan)
        lons_lat = lonsperlat(lat)
        jlonf    = (lan-1)*lonf
        rcs_loc  = sqrt(rcs2(min(lat,latg-lat+1)))

        do k=1,levs
          do i=1,lons_lat
            ilan = i+jlonf
            anl_gr_a_2(i,kau-1+k,lan) = grid_gr(ilan,g_u+k-1) * rcs_loc
            anl_gr_a_2(i,kav-1+k,lan) = grid_gr(ilan,g_v+k-1) * rcs_loc
            anl_gr_a_2(i,kat-1+k,lan) = grid_gr(ilan,g_t+k-1)
!           anl_gr_a_2(i,kadp-1+k,lan) = grid_gr(ilan,g_dpn+k-1)
          enddo
        enddo
!
        do i=1,lons_lat
          ilan = i+jlonf
          anl_gr_a_2(i,kap,lan) = grid_gr(ilan,g_zq)
        enddo
!
        do k=1,levh
          kk = kar-1+k
          do i=1,lons_lat
            ilan = i + jlonf
            anl_gr_a_2(i,kk,lan) = grid_gr(ilan,g_rt+k-1)
          enddo
        enddo
!
      enddo

      return
      end subroutine do_dynamics_gridn2anl_slg_dfi
! --------------------------------------------------------------
      subroutine do_dynamics_gridt2anl(grid_gr,anl_gr_a_2,rdt2,
     &                                 global_lats_a,lonsperlat)

      real(kind=kind_grid) grid_gr(lonf*lats_node_a_max,lotgr)
      real(kind=kind_evod) anl_gr_a_2(lonfx*lota,lats_dim_ext)
      integer,intent(in):: global_lats_a(latg)
      integer,intent(in):: lonsperlat(latg)

      real      rdt2
      integer   lan,lat,lon_dim,lons_lat,k,i,jlonf,ilan

!$omp parallel do private(lan)
!$omp+private(lat,lon_dim,lons_lat,jlonf,i,k,ilan)

      do lan=1,lats_node_a
        lat      = global_lats_a(ipt_lats_node_a-1+lan)
        lon_dim  = lon_dims_a(lan)
        lons_lat = lonsperlat(lat)
        jlonf    = (lan-1)*lonf
!!$omp parallel do private(i,k,ilan)
        do k=1,levs
          do i=1,lons_lat
            ilan = i+jlonf
            anl_gr_a_2(i+(kau -2+k)*lon_dim,lan) =
     &     (grid_gr(ilan,G_u  +k-1)-grid_gr(ilan,G_uum+k-1))*rdt2
            anl_gr_a_2(i+(kav -2+k)*lon_dim,lan) =
     &     (grid_gr(ilan,G_v  +k-1)-grid_gr(ilan,G_vvm+k-1))*rdt2
            anl_gr_a_2(i+(kat -2+k)*lon_dim,lan) =
     &     (grid_gr(ilan,G_t  +k-1)-grid_gr(ilan,G_ttm+k-1))*rdt2
            anl_gr_a_2(i+(kadp-2+k)*lon_dim,lan) =
     &     (grid_gr(ilan,G_dpn+k-1)-grid_gr(ilan,G_dpm+k-1))*rdt2
          enddo
        enddo
!!$omp parallel do private(i,ilan)
        do i=1,lons_lat
            ilan = i + jlonf
            anl_gr_a_2(i+(kaps-1)*lon_dim,lan) =
     &     (grid_gr(ilan,G_zq)-grid_gr(ilan,G_qm))*rdt2
        enddo
!
        if( .not. ndslfv ) then
!!$omp parallel do private(i,k,ilan)
          do k=1,levh
            do i=1,lons_lat
              ilan = i + jlonf
              anl_gr_a_2(i+(kar-2+k)*lon_dim,lan) =
     &       (grid_gr(ilan,G_rt +k-1)-grid_gr(ilan,G_rm +k-1))*rdt2
            enddo
          enddo
        endif
!
      enddo

      return
      end subroutine do_dynamics_gridt2anl
!
! --------------------------------------------------------------
      subroutine do_dynamics_gridt2rqt(grid_gr,rqt_gr_a_2,
     &                                 global_lats_a,lonsperlat)

      real(kind=kind_grid) grid_gr(lonf*lats_node_a_max,lotgr)
      real(kind=kind_evod) rqt_gr_a_2(lonfx*levs,lats_dim_ext)
      integer,intent(in):: global_lats_a(latg)
      integer,intent(in):: lonsperlat(latg)

      integer   lan,lat,lon_dim,lons_lat,k,i
      integer   jlonf,ilan

!$omp parallel do private(lan)
!$omp+private(lat,lon_dim,lons_lat,jlonf,i,k,ilan)

      do lan=1,lats_node_a
        lat = global_lats_a(ipt_lats_node_a-1+lan)
        lon_dim = lon_dims_a(lan)
        lons_lat = lonsperlat(lat)
        jlonf = (lan-1)*lonf
        do k=1,levs
          do i=1,lons_lat
            ilan=i+jlonf
            rqt_gr_a_2(i+(k-1)*lon_dim,lan)=
     &      grid_gr(ilan,G_rt +k-1)-grid_gr(ilan,G_rq +k-1)
          enddo
        enddo
!
      enddo

      return
      end subroutine do_dynamics_gridt2rqt
!
! --------------------------------------------------------------
      subroutine do_dynamics_gridm2rqt(grid_gr,rqt_gr_a_2,
     &                                 global_lats_a,lonsperlat)

      real(kind=kind_grid) grid_gr(lonf*lats_node_a_max,lotgr)
      real(kind=kind_evod) rqt_gr_a_2(lonfx*levh,lats_dim_ext)
      integer,intent(in):: global_lats_a(latg)
      integer,intent(in):: lonsperlat(latg)

      integer   lan,lat,lon_dim,lons_lat,k,i
      integer   jlonf,ilan

!$omp parallel do private(lan)
!$omp+private(lat,lon_dim,lons_lat,jlonf,i,k,ilan)

      do lan=1,lats_node_a
        lat = global_lats_a(ipt_lats_node_a-1+lan)
        lon_dim = lon_dims_a(lan)
        lons_lat = lonsperlat(lat)
        jlonf = (lan-1)*lonf
        do k=1,levh
          do i=1,lons_lat
            ilan=i+jlonf
            rqt_gr_a_2(i+(k-1)*lon_dim,lan)=
     &      grid_gr(ilan,G_rm +k-1)
          enddo
        enddo
!
      enddo

      return
      end subroutine do_dynamics_gridm2rqt
!
! --------------------------------------------------------------
      subroutine do_dynamics_gridc2rqt(grid_gr,rqt_gr_a_2,
     &                                 global_lats_a,lonsperlat)

      real(kind=kind_grid) grid_gr(lonf*lats_node_a_max,lotgr)
      real(kind=kind_evod) rqt_gr_a_2(lonfx*levh,lats_dim_ext)
      integer,intent(in):: global_lats_a(latg)
      integer,intent(in):: lonsperlat(latg)

      integer   lan,lat,lon_dim,lons_lat,k,i
      integer   jlonf,ilan

!$omp parallel do private(lan)
!$omp+private(lat,lon_dim,lons_lat,jlonf,i,k,ilan)

      do lan=1,lats_node_a
        lat = global_lats_a(ipt_lats_node_a-1+lan)
        lon_dim = lon_dims_a(lan)
        lons_lat = lonsperlat(lat)
        jlonf = (lan-1)*lonf
        do k=1,levh
          do i=1,lons_lat
            ilan=i+jlonf
            rqt_gr_a_2(i+(k-1)*lon_dim,lan)=
     &      grid_gr(ilan,G_rq +k-1)
          enddo
        enddo
!
      enddo

      return
      end subroutine do_dynamics_gridc2rqt

!----------------------------------------------------------
      subroutine do_dynamics_spectaddlapgz(trie_ls,trio_ls)

      real(kind=kind_evod) trie_ls(len_trie_ls,2,lotls)
      real(kind=kind_evod) trio_ls(len_trio_ls,2,lotls)

      integer k,i

!$omp parallel do private(i,k)
      do k=1,levs
         do i=1,len_trie_ls
            trie_ls(i,1,P_x  +k-1)=
     &      trie_ls(i,1,P_x  +k-1)+trie_ls(i,1,P_lapgz)
            trie_ls(i,2,P_x  +k-1)=
     &      trie_ls(i,2,P_x  +k-1)+trie_ls(i,2,P_lapgz)
         enddo
         do i=1,len_trio_ls
            trio_ls(i,1,P_x  +k-1)=
     &      trio_ls(i,1,P_x  +k-1)+trio_ls(i,1,P_lapgz)
            trio_ls(i,2,P_x  +k-1)=
     &      trio_ls(i,2,P_x  +k-1)+trio_ls(i,2,P_lapgz)
         enddo
      enddo

      return
      end subroutine do_dynamics_spectaddlapgz
!----------------------------------------------------------
      subroutine do_dynamics_spectupdatewrt(trie_ls,trio_ls,dt2)

      real(kind=kind_evod) trie_ls(len_trie_ls,2,lotls)
      real(kind=kind_evod) trio_ls(len_trio_ls,2,lotls)
      real,   intent(in):: dt2

      integer	k,i

!$omp parallel do private(i,k)
      do k=1,levs
         do i=1,len_trie_ls
            trie_ls(i,1,P_w  +k-1)=
     &      trie_ls(i,1,P_zem+k-1)+dt2*trie_ls(i,1,P_w+k-1)
            trie_ls(i,2,P_w  +k-1)=
     &      trie_ls(i,2,P_zem+k-1)+dt2*trie_ls(i,2,P_w+k-1)
         enddo
         do i=1,len_trio_ls
            trio_ls(i,1,P_w  +k-1)=
     &      trio_ls(i,1,P_zem+k-1)+dt2*trio_ls(i,1,P_w+k-1)
            trio_ls(i,2,P_w  +k-1)=
     &      trio_ls(i,2,P_zem+k-1)+dt2*trio_ls(i,2,P_w+k-1)
         enddo
      enddo
!
      if( .not. ndslfv ) then

!$omp parallel do private(i,k)
      do k=1,levh
         do i=1,len_trie_ls
            trie_ls(i,1,P_rt+k-1)=
     &      trie_ls(i,1,P_rm+k-1)+dt2* trie_ls(i,1,P_rt+k-1)  
            trie_ls(i,2,P_rt+k-1)=
     &      trie_ls(i,2,P_rm+k-1)+dt2* trie_ls(i,2,P_rt+k-1) 
         enddo
         do i=1,len_trio_ls
            trio_ls(i,1,P_rt+k-1)=
     &      trio_ls(i,1,P_rm+k-1)+dt2* trio_ls(i,1,P_rt+k-1)
            trio_ls(i,2,P_rt+k-1)=
     &      trio_ls(i,2,P_rm+k-1)+dt2* trio_ls(i,2,P_rt+k-1)
         enddo
      enddo

      endif

      return
      end subroutine do_dynamics_spectupdatewrt
!
!----------------------------------------------------------
      subroutine do_dynamics_spectupdatexydpnzq(trie_ls,trio_ls,dt2)
 
      real(kind=kind_evod) trie_ls(len_trie_ls,2,lotls)
      real(kind=kind_evod) trio_ls(len_trio_ls,2,lotls)
      real,   intent(in):: dt2

      integer	k,i

!$omp parallel do private(i,k)

      do k=1,levs                                                       
         do i=1,len_trie_ls                                             
            trie_ls(i,1,P_x  +k-1)=                                     
     &      trie_ls(i,1,P_dim+k-1)+dt2*trie_ls(i,1,P_x  +k-1)    
            trie_ls(i,2,P_x  +k-1)=                                     
     &      trie_ls(i,2,P_dim+k-1)+dt2*trie_ls(i,2,P_x  +k-1)    
            trie_ls(i,1,P_y  +k-1)=                                     
     &      trie_ls(i,1,P_tem+k-1)+dt2*trie_ls(i,1,P_y  +k-1)    
            trie_ls(i,2,P_y  +k-1)=                                     
     &      trie_ls(i,2,P_tem+k-1)+dt2*trie_ls(i,2,P_y  +k-1)    
            trie_ls(i,1,P_dpn+k-1)=                                     
     &      trie_ls(i,1,P_dpm+k-1)+dt2*trie_ls(i,1,P_dpn+k-1)    
            trie_ls(i,2,P_dpn+k-1)=                                     
     &      trie_ls(i,2,P_dpm+k-1)+dt2*trie_ls(i,2,P_dpn+k-1)    
         enddo                                                          
         do i=1,len_trio_ls                                             
            trio_ls(i,1,P_x  +k-1)=                                     
     &      trio_ls(i,1,P_dim+k-1)+dt2*trio_ls(i,1,P_x  +k-1)    
            trio_ls(i,2,P_x  +k-1)=                                     
     &      trio_ls(i,2,P_dim+k-1)+dt2*trio_ls(i,2,P_x  +k-1)    
            trio_ls(i,1,P_y  +k-1)=                                     
     &      trio_ls(i,1,P_tem+k-1)+dt2*trio_ls(i,1,P_y  +k-1)    
            trio_ls(i,2,P_y  +k-1)=                                     
     &      trio_ls(i,2,P_tem+k-1)+dt2*trio_ls(i,2,P_y  +k-1)    
            trio_ls(i,1,P_dpn+k-1)=                                     
     &      trio_ls(i,1,P_dpm+k-1)+dt2*trio_ls(i,1,P_dpn+k-1)    
            trio_ls(i,2,P_dpn+k-1)=                                     
     &      trio_ls(i,2,P_dpm+k-1)+dt2*trio_ls(i,2,P_dpn+k-1)    
         enddo                                                          
      enddo                                                             
!
         do i=1,len_trie_ls                                             
            trie_ls(i,1,P_zq)=                                          
     &      trie_ls(i,1,P_qm)+dt2*trie_ls(i,1,P_zq)            
            trie_ls(i,2,P_zq)=                                          
     &      trie_ls(i,2,P_qm)+dt2*trie_ls(i,2,P_zq)            
         enddo                                                          
         do i=1,len_trio_ls                                             
            trio_ls(i,1,P_zq)=                                          
     &      trio_ls(i,1,P_qm)+dt2*trio_ls(i,1,P_zq)            
            trio_ls(i,2,P_zq)=                                          
     &      trio_ls(i,2,P_qm)+dt2*trio_ls(i,2,P_zq)            
         enddo                                                          

      return 
      end subroutine do_dynamics_spectupdatexydpnzq

!--------------------------------------------
      subroutine do_dynamics_spectn2c(trie_ls,trio_ls) 

      real(kind=kind_evod) trie_ls(len_trie_ls,2,lotls)
      real(kind=kind_evod) trio_ls(len_trio_ls,2,lotls)

      integer	k,j

!$omp parallel do private(j,k)
      DO K=1,LEVS
       DO J=1,LEN_TRIE_LS
         TRIE_LS(J,1,P_DI+K-1)=TRIE_LS(J,1,P_X  +K-1)
         TRIE_LS(J,2,P_DI+K-1)=TRIE_LS(J,2,P_X  +K-1)
         TRIE_LS(J,1,P_ZE+K-1)=TRIE_LS(J,1,P_W  +K-1)
         TRIE_LS(J,2,P_ZE+K-1)=TRIE_LS(J,2,P_W  +K-1)
         TRIE_LS(J,1,P_TE+K-1)=TRIE_LS(J,1,P_Y  +K-1)
         TRIE_LS(J,2,P_TE+K-1)=TRIE_LS(J,2,P_Y  +K-1)
         TRIE_LS(J,1,P_dp+K-1)=TRIE_LS(J,1,P_dpn+K-1)
         TRIE_LS(J,2,P_dp+K-1)=TRIE_LS(J,2,P_dpn+K-1)
       ENDDO
       DO J=1,LEN_TRIO_LS
         TRIO_LS(J,1,P_DI+K-1)=TRIO_LS(J,1,P_X  +K-1)
         TRIO_LS(J,2,P_DI+K-1)=TRIO_LS(J,2,P_X  +K-1)
         TRIO_LS(J,1,P_ZE+K-1)=TRIO_LS(J,1,P_W  +K-1)
         TRIO_LS(J,2,P_ZE+K-1)=TRIO_LS(J,2,P_W  +K-1)
         TRIO_LS(J,1,P_TE+K-1)=TRIO_LS(J,1,P_Y  +K-1)
         TRIO_LS(J,2,P_TE+K-1)=TRIO_LS(J,2,P_Y  +K-1)
         TRIO_LS(J,1,P_dp+K-1)=TRIO_LS(J,1,P_dpn+K-1)
         TRIO_LS(J,2,P_dp+K-1)=TRIO_LS(J,2,P_dpn+K-1)
       ENDDO
      ENDDO

      DO J=1,LEN_TRIE_LS
         TRIE_LS(J,1,P_Q)=TRIE_LS(J,1,P_ZQ)
         TRIE_LS(J,2,P_Q)=TRIE_LS(J,2,P_ZQ)
      ENDDO
      DO J=1,LEN_TRIO_LS
         TRIO_LS(J,1,P_Q)=TRIO_LS(J,1,P_ZQ)
         TRIO_LS(J,2,P_Q)=TRIO_LS(J,2,P_ZQ)
      ENDDO
!
      if( .not. ndslfv ) then

!$omp parallel do private(j,k)
      DO K=1,LEVH
       DO J=1,LEN_TRIE_LS
         TRIE_LS(J,1,P_RQ+K-1)=TRIE_LS(J,1,P_RT+K-1)
         TRIE_LS(J,2,P_RQ+K-1)=TRIE_LS(J,2,P_RT+K-1)
       ENDDO
       DO J=1,LEN_TRIO_LS
         TRIO_LS(J,1,P_RQ+K-1)=TRIO_LS(J,1,P_RT+K-1)
         TRIO_LS(J,2,P_RQ+K-1)=TRIO_LS(J,2,P_RT+K-1)
       ENDDO
      ENDDO

      endif

      return
      end subroutine do_dynamics_spectn2c

!--------------------------------------------
      subroutine do_dynamics_spectn2m(trie_ls,trio_ls) 

      real(kind=kind_evod) trie_ls(len_trie_ls,2,lotls)
      real(kind=kind_evod) trio_ls(len_trio_ls,2,lotls)

      integer	k,j

!$omp parallel do private(j,k)
      DO K=1,LEVS
        DO J=1,LEN_TRIE_LS
          TRIE_LS(J,1,P_DIM+K-1) = TRIE_LS(J,1,P_X  +K-1)
          TRIE_LS(J,2,P_DIM+K-1) = TRIE_LS(J,2,P_X  +K-1)
          TRIE_LS(J,1,P_ZEM+K-1) = TRIE_LS(J,1,P_W  +K-1)
          TRIE_LS(J,2,P_ZEM+K-1) = TRIE_LS(J,2,P_W  +K-1)
          TRIE_LS(J,1,P_TEM+K-1) = TRIE_LS(J,1,P_Y  +K-1)
          TRIE_LS(J,2,P_TEM+K-1) = TRIE_LS(J,2,P_Y  +K-1)
        ENDDO
        DO J=1,LEN_TRIO_LS
          TRIO_LS(J,1,P_DIM+K-1) = TRIO_LS(J,1,P_X  +K-1)
          TRIO_LS(J,2,P_DIM+K-1) = TRIO_LS(J,2,P_X  +K-1)
          TRIO_LS(J,1,P_ZEM+K-1) = TRIO_LS(J,1,P_W  +K-1)
          TRIO_LS(J,2,P_ZEM+K-1) = TRIO_LS(J,2,P_W  +K-1)
          TRIO_LS(J,1,P_TEM+K-1) = TRIO_LS(J,1,P_Y  +K-1)
          TRIO_LS(J,2,P_TEM+K-1) = TRIO_LS(J,2,P_Y  +K-1)
        ENDDO
        if( ndslfv ) then
          DO J=1,LEN_TRIE_LS
            TRIE_LS(J,1,P_dpm+K-1) = TRIE_LS(J,1,P_dpn+K-1)
            TRIE_LS(J,2,P_dpm+K-1) = TRIE_LS(J,2,P_dpn+K-1)
          ENDDO
          DO J=1,LEN_TRIO_LS
            TRIO_LS(J,2,P_dpm+K-1) = TRIO_LS(J,2,P_dpn+K-1)
            TRIO_LS(J,2,P_dpm+K-1) = TRIO_LS(J,2,P_dpn+K-1)
          ENDDO
        endif
      ENDDO

      DO J=1,LEN_TRIE_LS
        TRIE_LS(J,1,P_QM) = TRIE_LS(J,1,P_ZQ)
        TRIE_LS(J,2,P_QM) = TRIE_LS(J,2,P_ZQ)
      ENDDO
      DO J=1,LEN_TRIO_LS
        TRIO_LS(J,1,P_QM) = TRIO_LS(J,1,P_ZQ)
        TRIO_LS(J,2,P_QM) = TRIO_LS(J,2,P_ZQ)
      ENDDO
!
      if( .not. ndslfv ) then

!$omp parallel do private(j,k)
        DO K=1,LEVH
          DO J=1,LEN_TRIE_LS
            TRIE_LS(J,1,P_RM+K-1) = TRIE_LS(J,1,P_RT+K-1)
            TRIE_LS(J,2,P_RM+K-1) = TRIE_LS(J,2,P_RT+K-1)
          ENDDO
          DO J=1,LEN_TRIO_LS
            TRIO_LS(J,1,P_RM+K-1) = TRIO_LS(J,1,P_RT+K-1)
            TRIO_LS(J,2,P_RM+K-1) = TRIO_LS(J,2,P_RT+K-1)
          ENDDO
        ENDDO

      endif

      return
      end subroutine do_dynamics_spectn2m

!--------------------------------------------
      subroutine do_dynamics_spectc2n(trie_ls,trio_ls) 

      real(kind=kind_evod) trie_ls(len_trie_ls,2,lotls)
      real(kind=kind_evod) trio_ls(len_trio_ls,2,lotls)

      integer	k,j

!$omp parallel do private(j,k)
      DO K=1,LEVS
        DO J=1,LEN_TRIE_LS
          TRIE_LS(J,1,P_X  +K-1) = TRIE_LS(J,1,P_DI+K-1)
          TRIE_LS(J,2,P_X  +K-1) = TRIE_LS(J,2,P_DI+K-1)
          TRIE_LS(J,1,P_W  +K-1) = TRIE_LS(J,1,P_ZE+K-1)
          TRIE_LS(J,2,P_W  +K-1) = TRIE_LS(J,2,P_ZE+K-1)
          TRIE_LS(J,1,P_Y  +K-1) = TRIE_LS(J,1,P_TE+K-1)
          TRIE_LS(J,2,P_Y  +K-1) = TRIE_LS(J,2,P_TE+K-1)
        ENDDO
        DO J=1,LEN_TRIO_LS
          TRIO_LS(J,1,P_X  +K-1) = TRIO_LS(J,1,P_DI+K-1)
          TRIO_LS(J,2,P_X  +K-1) = TRIO_LS(J,2,P_DI+K-1)
          TRIO_LS(J,1,P_W  +K-1) = TRIO_LS(J,1,P_ZE+K-1)
          TRIO_LS(J,2,P_W  +K-1) = TRIO_LS(J,2,P_ZE+K-1)
          TRIO_LS(J,1,P_Y  +K-1) = TRIO_LS(J,1,P_TE+K-1)
          TRIO_LS(J,2,P_Y  +K-1) = TRIO_LS(J,2,P_TE+K-1)
        ENDDO
        if( ndslfv ) then
          DO J=1,LEN_TRIE_LS
            TRIE_LS(J,1,P_dpn+K-1) = TRIE_LS(J,1,P_dp+K-1)
            TRIE_LS(J,2,P_dpn+K-1) = TRIE_LS(J,2,P_dp+K-1)
          ENDDO
          DO J=1,LEN_TRIO_LS
            TRIO_LS(J,1,P_dpn+K-1) = TRIO_LS(J,1,P_dp+K-1)
            TRIO_LS(J,2,P_dpn+K-1) = TRIO_LS(J,2,P_dp+K-1)
          ENDDO
        endif
      ENDDO

      DO J=1,LEN_TRIE_LS
        TRIE_LS(J,1,P_ZQ) = TRIE_LS(J,1,P_Q)
        TRIE_LS(J,2,P_ZQ) = TRIE_LS(J,2,P_Q)
      ENDDO
      DO J=1,LEN_TRIO_LS
        TRIO_LS(J,1,P_ZQ) = TRIO_LS(J,1,P_Q)
        TRIO_LS(J,2,P_ZQ) = TRIO_LS(J,2,P_Q)
      ENDDO
!
      if( .not. ndslfv ) then

!$omp parallel do private(j,k)
        DO K=1,LEVH
          DO J=1,LEN_TRIE_LS
            TRIE_LS(J,1,P_RT+K-1) = TRIE_LS(J,1,P_RQ+K-1)
            TRIE_LS(J,2,P_RT+K-1) = TRIE_LS(J,2,P_RQ+K-1)
          ENDDO
          DO J=1,LEN_TRIO_LS
            TRIO_LS(J,1,P_RT+K-1) = TRIO_LS(J,1,P_RQ+K-1)
            TRIO_LS(J,2,P_RT+K-1) = TRIO_LS(J,2,P_RQ+K-1)
          ENDDO
        ENDDO

      endif

      return
      end subroutine do_dynamics_spectc2n

!--------------------------------------------
      subroutine do_dynamics_spectdfini(kst,nst,trie_ls,trio_ls,
     &           lvl_change)

      real(kind=kind_evod) trie_ls(len_trie_ls,2,lotls)
      real(kind=kind_evod) trio_ls(len_trio_ls,2,lotls)
      integer kst,nst,lvl_change

      real, allocatable, save :: trie_save(:,:,:)
      real, allocatable, save :: trio_save(:,:,:)
      real, save :: totsum
      integer, save ::  kze,kdi,kte,krq,kdp,kq,lotsave

      integer i,k,kl,ks,kk,knodfis,knodfie
      real totsumi,sx,wx,digfil
      logical ldodfi,lsavehdfi


      if( me==0 ) then
        print *,' Enter do_dynamics_spectdfini kst=',kst,' nst=',nst
      endif

      if(kst.lt.-nst) then
        totsum = 0.0
        kze=1
        kdi=kze+levs
        kte=kdi+levs
        if( ndslfv ) then
          kdp=kte+levs
        else
          krq=kte+levs
          kdp=krq+levh
        endif
        kq =kdp+levs
        lotsave = kq
        allocate ( trie_save(len_trie_ls,2,lotsave) )
        allocate ( trio_save(len_trio_ls,2,lotsave) )
        trie_save(:,:,:) = 0.0
        trio_save(:,:,:) = 0.0
!        if( me==0 ) then
!          print *,' Save initial values in do_dynamics_spectdfini.'
!          print *,' totsum kze kdi kte kq lotsave ',
!     &              totsum,kze,kdi,kte,kq,lotsave
!        endif
        return
      endif

      if(kst.le.nst) then
        sx = acos(-1.)*kst/nst
        wx = acos(-1.)*kst/(nst+1)
        if(kst.ne.0) then
          digfil = sin(wx)/wx*sin(sx)/sx
        else
          digfil = 1.0
        endif
        totsum = totsum + digfil
!        if(me.eq.0)then
!          print*,'in dfini sx=',sx,'wx=',wx,'digfil=',digfil,
!     &      'at kst=',kst
!          print *,' totsum kze kdi kte kq lotsave ',
!     &              totsum,kze,kdi,kte,kq,lotsave
!        endif

        lsavehdfi=.false.
        if (kst*2==nst) lsavehdfi=.true.
        do ks=1,lotsave
          ldodfi=.true.
          if (lvl_change/=levs) then
            kk=(ks-1)/levs
            knodfis=kk*levs+lvl_change
            knodfie=kk*levs+levs
            if (ks>knodfis .and. ks<knodfie) then
              ldodfi=.false.
            endif
          endif
          if (ldodfi) then
            kl = p_ze + ks - 1
            do i=1,len_trie_ls
              trie_save(i,1,ks)=trie_save(i,1,ks)+digfil*trie_ls(i,1,kl)
              trie_save(i,2,ks)=trie_save(i,2,ks)+digfil*trie_ls(i,2,kl)
            enddo
            do i=1,len_trio_ls
              trio_save(i,1,ks)=trio_save(i,1,ks)+digfil*trio_ls(i,1,kl)
              trio_save(i,2,ks)=trio_save(i,2,ks)+digfil*trio_ls(i,2,kl)
            enddo
          endif

          if (.not.ldodfi .and. lsavehdfi) then
            do i=1,len_trie_ls
              trie_save(i,1,ks)=trie_ls(i,1,kl)
              trie_save(i,2,ks)=trie_ls(i,2,kl)
            enddo
            do i=1,len_trio_ls
              trio_save(i,1,ks)=trio_ls(i,1,kl)
              trio_save(i,2,ks)=trio_ls(i,2,kl)
            enddo
          endif

        enddo


        return
      endif

      if(kst.gt.nst) then
        totsumi = 1. /totsum
        if(me.eq.0)then
          print*,'in dfini last step totsumi=',totsumi,
     &      'at kst=',kst
        endif
        if (kst*2==nst) lsavehdfi=.true.
        do ks=1,lotsave
          ldodfi=.true.
          if (lvl_change/=levs) then
            kk=(ks-1)/levs
            knodfis=kk*levs+lvl_change
            knodfie=kk*levs+levs
            if (ks>knodfis .and. ks<knodfie) then
              ldodfi=.false.
            endif
          endif
          if (ldodfi) then
            kl = p_ze + ks - 1
            do i=1,len_trie_ls
              trie_ls(i,1,kl) = trie_save(i,1,ks) * totsumi
              trie_ls(i,2,kl) = trie_save(i,2,ks) * totsumi
            enddo
            do i=1,len_trio_ls
              trio_ls(i,1,kl) = trio_save(i,1,ks) * totsumi
              trio_ls(i,2,kl) = trio_save(i,2,ks) * totsumi
            enddo
          endif
        enddo

        return
      endif

      return
      end subroutine do_dynamics_spectdfini

      subroutine do_dynamics_spectdfini_slg(kst,nst,trie_ls,trio_ls)

      use gfs_dyn_physcons, pi  => con_pi

      real(kind=kind_evod) trie_ls(len_trie_ls,2,lotls)
      real(kind=kind_evod) trio_ls(len_trio_ls,2,lotls)
      integer kst,nst

      real, allocatable, save :: trie_save(:,:,:), trio_save(:,:,:)
      real,    save :: totsum
      integer, save ::  kze,kdi,kte,krq,kdp,kq,lotsave

      integer i,k,kl,ks
      real    totsumi,sx,wx,digfil


!      if( me==0 ) then
!        print *,' Enter do_dynamics_spectdfini_slg kst=',kst,' nst=',nst
!      endif

      if(kst < -nst) then
        totsum = 0.0
        kze     = 1
        kdi     = kze + levs
        kte     = kdi + levs
        kq      = kte + levs
        krq     = kq  + 1
        lotsave = krq + levh - 1
        allocate ( trie_save(len_trie_ls,2,lotsave) )
        allocate ( trio_save(len_trio_ls,2,lotsave) )

!$omp parallel do private(i,ks)
        do ks=1,lotsave
          do i=1,len_trie_ls
            trie_save(i,1,ks) = 0.0
            trie_save(i,2,ks) = 0.0
          enddo
          do i=1,len_trio_ls
            trio_save(i,1,ks) = 0.0
            trio_save(i,2,ks) = 0.0
          enddo
        enddo

!        if( me==0 ) then
!          print *,' Save initial values in do_dynamics_spectdfini_slg'
!          print *,' totsum kze kdi kte kq krq lotsave ',
!     &              totsum,kze,kdi,kte,kq,krq,lotsave
!        endif
        return
      endif

      if(kst <= nst) then
        sx = pi*kst/nst
        wx = pi*kst/(nst+1)
        if(kst /= 0) then
          digfil = sin(wx)/wx*sin(sx)/sx
        else
          digfil = 1.0
        endif
        totsum = totsum + digfil
!        if (me == 0) then
!          print*,'in dfini sx=',sx,'wx=',wx,'digfil=',digfil,
!     &      'at kst=',kst
!          print *,' totsum kze kdi kte krq lotsave ',
!     &              totsum,kze,kdi,kte,krq,lotsave
!        endif

!$omp parallel do private(i,ks,kl)
        do ks=1,lotsave
          kl = p_ze + ks - 1
          if (ks >= krq) kl = p_rq + ks - krq
          do i=1,len_trie_ls
            trie_save(i,1,ks) = trie_save(i,1,ks)
     &                        + digfil*trie_ls(i,1,kl)
            trie_save(i,2,ks) = trie_save(i,2,ks)
     &                        + digfil*trie_ls(i,2,kl)
          enddo
          do i=1,len_trio_ls
            trio_save(i,1,ks) = trio_save(i,1,ks)
     &                        + digfil*trio_ls(i,1,kl)
            trio_save(i,2,ks) = trio_save(i,2,ks)
     &                        + digfil*trio_ls(i,2,kl)
          enddo
        enddo

        return
      endif

      if(kst > nst) then
        totsumi = 1. /totsum
        if(me == 0)then
          print*,'in dfini last step totsumi=',totsumi,
     &           'at kst=',kst
        endif

!$omp parallel do private(i,ks,kl)
        do ks=1,lotsave
          kl = p_ze + ks - 1
          IF(ks >= krq) kl = p_rq + ks - krq
          do i=1,len_trie_ls
            trie_ls(i,1,kl) = trie_save(i,1,ks) * totsumi
            trie_ls(i,2,kl) = trie_save(i,2,ks) * totsumi
          enddo
          do i=1,len_trio_ls
            trio_ls(i,1,kl) = trio_save(i,1,ks) * totsumi
            trio_ls(i,2,kl) = trio_save(i,2,ks) * totsumi
          enddo
        enddo

        deallocate(trie_save,trio_save)
        return
      endif

      return
      end subroutine do_dynamics_spectdfini_slg

!--------------------------------------------
      subroutine do_dynamics_syn2gridn(syn_gr_a_2,grid_gr,
     &                                 global_lats_a,lonsperlat)


      real(kind=kind_grid) grid_gr(lonf*lats_node_a_max,lotgr)
      real(kind=kind_evod) syn_gr_a_2(lonfx*lots,lats_dim_ext)
      integer,intent(in):: global_lats_a(latg)
      integer,intent(in):: lonsperlat(latg)

      integer   lan,lat,lon_dim,lons_lat,k,i,jlonf,ilan

!!$omp parallel do private(lan)
!!$omp+private(lat,lon_dim,lons_lat,jlonf,i,k,ilan)

      do lan=1,lats_node_a
        lat = global_lats_a(ipt_lats_node_a-1+lan)
        lon_dim = lon_dims_a(lan)
        lons_lat = lonsperlat(lat)
        jlonf = (lan-1)*lonf
!$omp parallel do private(i,k,ilan)
        do k=1,levs
          do i=1,lons_lat
            ilan = i+jlonf
            grid_gr(ilan,G_u  +k-1)=syn_gr_a_2(i+(ksu -2+k)*lon_dim,lan)
            grid_gr(ilan,G_v  +k-1)=syn_gr_a_2(i+(ksv -2+k)*lon_dim,lan)
            grid_gr(ilan,G_t  +k-1)=syn_gr_a_2(i+(kst -2+k)*lon_dim,lan)
            grid_gr(ilan,G_dpn+k-1)=syn_gr_a_2(i+(ksdp-2+k)*lon_dim,lan)
          enddo
        enddo
!$omp parallel do private(i,ilan)
        do i=1,lons_lat
            ilan=i+jlonf
            grid_gr(ilan,G_zq)= syn_gr_a_2(i+(ksq-1)*lon_dim,lan)
        enddo
!
        if( .not. ndslfv ) then
!$omp parallel do private(i,k,ilan)
        do k=1,levh
          do i=1,lons_lat
            ilan = i+jlonf
            grid_gr(ilan,G_rt+k-1)= syn_gr_a_2(i+(ksr-2+k)*lon_dim,lan)
          enddo
        enddo
        endif
!
      enddo

      return
      end subroutine do_dynamics_syn2gridn

!--------------------------------------------
!--------------------------------------------
      subroutine do_dynamics_syn2gridn_slg(syn_gr_a_2,grid_gr,
     &                                 ak5,bk5,kdt,
!    &                                 lonf,lots_sl,
!    &                                 lats_dim_a,
     &                                 global_lats_a,lonsperlat)


      integer kdt
      real(kind=kind_grid) grid_gr(lonf*lats_node_a_max,lotgr)
      real(kind=kind_evod) syn_gr_a_2(lonf,lota,lats_dim_a)
      real(kind=kind_evod), dimension(levp1) :: ak5, bk5
      integer,intent(in), dimension(latg) :: global_lats_a, lonsperlat

      real(kind=kind_evod) wrk(lonf), tx1(lonf), tem1
      integer   lan,lat,lon_dim,lons_lat,k,i
      integer   jlonf, ilan, ksq, ksu, ksv, kst

      ksq = 0*levs + 0*levh + 1
      ksu = 0*levs + 0*levh + 4
      ksv = 1*levs + 0*levh + 4
      kst = 4*levs + 0*levh + 4
      ksr = 5*levs + 0*levh + 4

!$omp parallel do private(lan)
!$omp+private(lat,lon_dim,lons_lat,jlonf,i,k,ilan,tem1,wrk,tx1)

      do lan=1,lats_node_a
        lat      = global_lats_a(ipt_lats_node_a-1+lan)
        lon_dim  = lon_dims_a(lan)
        lons_lat = lonsperlat(lat)
        jlonf    = (lan-1)*lonf

        do i=1,lons_lat
            ilan = i + jlonf
            grid_gr(ilan,G_zq) = syn_gr_a_2(i,ksq,lan)
            tx1(i)             = exp(syn_gr_a_2(i,ksq,lan))
            wrk(i)             = ak5(1) + bk5(1) * tx1(i)
        enddo
        do k=1,levs
          do i=1,lons_lat
            ilan = i + jlonf
            grid_gr(ilan,G_u  +k-1)    = syn_gr_a_2(i,ksu -1+k,lan)
            grid_gr(ilan,G_v  +k-1)    = syn_gr_a_2(i,ksv -1+k,lan)
            grid_gr(ilan,G_t  +k-1)    = syn_gr_a_2(i,kst -1+k,lan)
!           grid_gr(ilan,G_dpn+k-1)    = syn_gr_a_2(i,ksdp-1+k,lan)
            tem1                       = wrk(i)
            wrk(i)                     = ak5(k+1) + bk5(k+1) * tx1(i)
            grid_gr(ilan,G_dpn+levs-k) = wrk(i) - tem1
          enddo
        enddo
        do k=1,levh
          do i=1,lons_lat
            ilan = i + jlonf
            grid_gr(ilan,G_rt+k-1) = syn_gr_a_2(i,ksr-1+k,lan)
          enddo
        enddo
!
      enddo

      return
      end subroutine do_dynamics_syn2gridn_slg

      subroutine do_dynamics_gridomega(syn_gr_a_2,dyn_gr_a_2,
     &                                 pyn_gr_a_2,grid_gr,rcs2,
     &                                 global_lats_a,lonsperlat)

      use namelist_dynamics_def
      real(kind=kind_grid) grid_gr(lonf*lats_node_a_max,lotgr)
      real(kind=kind_evod) syn_gr_a_2(lonfx*lots,lats_dim_ext)
      real(kind=kind_evod) dyn_gr_a_2(lonfx*lotd,lats_dim_ext)
      real(kind=kind_evod) pyn_gr_a_2(lonfx*lotp,lats_dim_ext)
      real(kind=kind_grid) rcs2(latg2)
      integer,intent(in) :: global_lats_a(latg)
      integer,intent(in) :: lonsperlat(latg)

      real(kind=kind_grid)  ugr (ngptc,levs), vgr (ngptc,levs)
      real(kind=kind_grid)  gtv (ngptc,levs), gd  (ngptc,levs)
      real(kind=kind_grid)  gtvx(ngptc,levs), gtvy(ngptc,levs)
      real(kind=kind_grid)  gdpf(ngptc,levs), gdpl(ngptc,levs)
      real(kind=kind_grid)  gphi(ngptc)     , glam(ngptc)    , gq(ngptc)
      real(kind=kind_grid)  gdp (ngptc,levs)
      real(kind=kind_grid)  vvel(ngptc,levs), rcs2_lan

      integer   lan,lat,lon_dim,lons_lat,k,i,lon,ii,njeff
      integer   jlonf

!!$omp parallel do private(lan)
!!$omp+private(lat,lon_dim,lons_lat,jlonf,i,k)
!!$omp+private(ugr,vgr,gd,gtv,gtvx,gtvy,gdp,gdpf,gdpl)
!!$omp+private(gq,gphi,glam,vvel,rcs2_lan)

      do lan=1,lats_node_a
        lat      = global_lats_a(ipt_lats_node_a-1+lan)
        lon_dim  = lon_dims_a(lan)
        lons_lat = lonsperlat(lat)
        jlonf    = (lan-1)*lonf
        rcs2_lan = rcs2(min(lat,latg-lat+1))

!$omp parallel do schedule(dynamic,1) private(lon)
!$omp+private(i,k,ii,njeff)
!$omp+private(ugr,vgr,gd,gtv,gtvx,gtvy,gdp,gdpf,gdpl)
!$omp+private(gq,gphi,glam,vvel)
        do lon=1,lons_lat,ngptc
          njeff = min(ngptc,lons_lat-lon+1)
          do k=1,levs
            do i=1,njeff
              ii = lon - 1 + i
              ugr (i,k) = syn_gr_a_2(ii+(ksu-2+k)*lon_dim,lan)
              vgr (i,k) = syn_gr_a_2(ii+(ksv-2+k)*lon_dim,lan)
              gd  (i,k) = syn_gr_a_2(ii+(ksd-2+k)*lon_dim,lan)
              gtv (i,k) = syn_gr_a_2(ii+(kst-2+k)*lon_dim,lan)
              gtvx(i,k) = dyn_gr_a_2(ii+(kdtlam-2+k)*lon_dim,lan)
              gtvy(i,k) = dyn_gr_a_2(ii+(kdtphi-2+k)*lon_dim,lan)
              gdp (i,k) = syn_gr_a_2(ii+(ksdp-2+k)*lon_dim,lan)
              gdpf(i,k) = pyn_gr_a_2(ii+(kdpphi-2+k)*lon_dim,lan)
              gdpl(i,k) = pyn_gr_a_2(ii+(kdplam-2+k)*lon_dim,lan)
            enddo
          enddo
          do i=1,njeff
            ii = lon - 1 + i
            gq  (i) = syn_gr_a_2(ii+(ksq   -1)*lon_dim,lan)
            gphi(i) = syn_gr_a_2(ii+(kspphi-1)*lon_dim,lan)
            glam(i) = syn_gr_a_2(ii+(ksplam-1)*lon_dim,lan)
          enddo

          if( gen_coord_hybrid ) then 
            if( mass_dp ) then
              call omega_gcdp(njeff,ngptc,rcs2_lan,
     &                     gdp,gdpf,gdpl,gd,ugr,vgr,vvel)
            else
              call omega_gch(njeff,ngptc,rcs2_lan,
     &                     gq,gphi,glam,gtv,gtvx,gtvy,gd,ugr,vgr,vvel)
            endif
          else if( hybrid )then
            call omega_hyb(njeff,ngptc,rcs2_lan,
     &                    gq,gphi,glam,gd,ugr,vgr,vvel)
          else
            call omega_sig(njeff,ngptc,rcs2_lan,
     &                    gq,gphi,glam,gd,ugr,vgr,vvel)
          endif

          do k=1,levs
            do i=1,njeff
              grid_gr(jlonf+lon-1+i,g_dpdt+k-1) = vvel(i,k)
            enddo
          enddo
        enddo

      enddo

      return
      end subroutine do_dynamics_gridomega

      subroutine do_dynamics_gridomega_slg(syn_gr_a_2,grid_gr,rcs2,
!    &                                     lon_dim_a,lots_sl,
!    &                                     lats_dim_a,
     &                                     global_lats_a,lonsperlat)

      use namelist_dynamics_def

      real(kind=kind_grid) grid_gr(lonf*lats_node_a_max,lotgr)
!     real(kind=kind_evod) syn_gr_a_2(lonfx*lots_slg,lats_dim_ext)
      real(kind=kind_evod) syn_gr_a_2(lonf,lota,lats_dim_a)
      real(kind=kind_grid) rcs2(latg2)
      integer,intent(in), dimension(latg) :: global_lats_a, lonsperlat

      real(kind=kind_grid), dimension(ngptc,levs) :: ugr, vgr, gd, vvel
      real(kind=kind_grid), dimension(ngptc)      :: gphi, glam, gq
!     real(kind=kind_grid)  gtv (ngptc,levs), gd  (ngptc,levs)
!     real(kind=kind_grid)  gtvx(ngptc,levs), gtvy(ngptc,levs)
      real(kind=kind_grid)  rcs2_lan

      integer   lan,lat,lon_dim,lons_lat,k,i,lon,ii,njeff
      integer   jlonf, ksq, ksplam, kspphi, ksu, ksv, ksz, ksd, kst

      ksq     = 0*levs+0*levh+1
      ksplam  = 0*levs+0*levh+2
      kspphi  = 0*levs+0*levh+3
      ksu     = 0*levs+0*levh+4
      ksv     = 1*levs+0*levh+4
      ksz     = 2*levs+0*levh+4
      ksd     = 3*levs+0*levh+4
      kst     = 4*levs+0*levh+4


      do lan=1,lats_node_a
        lat      = global_lats_a(ipt_lats_node_a-1+lan)
        lon_dim  = lon_dims_a(lan)
        lons_lat = lonsperlat(lat)
        jlonf    = (lan-1)*lonf
        rcs2_lan = rcs2(min(lat,latg-lat+1))

!$omp parallel do schedule(dynamic,1) private(lon)
!$omp+private(i,k,ii,njeff)
!$omp+private(ugr,vgr,gd,gq,gphi,glam,vvel)
        do lon=1,lons_lat,ngptc
          njeff = min(ngptc,lons_lat-lon+1)
          do k=1,levs
            do i=1,njeff
              ii = lon - 1 + i
              ugr (i,k) = syn_gr_a_2(ii, ksu-1+k, lan)
              vgr (i,k) = syn_gr_a_2(ii, ksv-1+k, lan)
              gd  (i,k) = syn_gr_a_2(ii, ksd-1+k, lan)
            enddo
          enddo
          do i=1,njeff
            ii = lon - 1 + i
            gq  (i) = syn_gr_a_2(ii, ksq,    lan)
            gphi(i) = syn_gr_a_2(ii, kspphi, lan)
            glam(i) = syn_gr_a_2(ii, ksplam, lan)
          enddo

          call omega_hyb(njeff,ngptc,rcs2_lan,
     &                   gq,gphi,glam,gd,ugr,vgr,vvel)

          do k=1,levs
            do i=1,njeff
              grid_gr(jlonf+lon-1+i,g_dpdt+k-1) = vvel(i,k)
            enddo
          enddo
        enddo

      enddo

      return
      end subroutine do_dynamics_gridomega_slg
!-------------------------------------------------------------------
      subroutine do_dynamics_gridfilter(grid_gr,filta,filtb,
     &                                  global_lats_a,lonsperlat)

      real(kind=kind_grid) grid_gr(lonf*lats_node_a_max,lotgr)
      integer,intent(in):: global_lats_a(latg)
      integer,intent(in):: lonsperlat(latg)
      real,   intent(in):: filta, filtb

      integer   lan,lat,lons_lat,k,i
      integer   jlonf,ilan

!$omp parallel do private(lan,lat,lons_lat,jlonf,i,k,ilan)
      do lan=1,lats_node_a
        lat = global_lats_a(ipt_lats_node_a-1+lan)
        lons_lat = lonsperlat(lat)
        jlonf = (lan-1)*lonf
        do k=1,levs
          do i=1,lons_lat
            ilan=i+jlonf

            grid_gr(ilan,G_uum+k-1)=grid_gr(ilan,G_uu +k-1) *filta+
     &     (grid_gr(ilan,G_uum+k-1)+grid_gr(ilan,G_u  +k-1))*filtb
            grid_gr(ilan,G_vvm+k-1)=grid_gr(ilan,G_vv +k-1) *filta+
     &     (grid_gr(ilan,G_vvm+k-1)+grid_gr(ilan,G_v  +k-1))*filtb
            grid_gr(ilan,G_ttm+k-1)=grid_gr(ilan,G_tt +k-1) *filta+
     &     (grid_gr(ilan,G_ttm+k-1)+grid_gr(ilan,G_t  +k-1))*filtb
            grid_gr(ilan,G_dpm+k-1)=grid_gr(ilan,G_dp +k-1) *filta+
     &     (grid_gr(ilan,G_dpm+k-1)+grid_gr(ilan,G_dpn+k-1))*filtb

            grid_gr(ilan,G_uu +k-1)=grid_gr(ilan,G_u  +k-1)
            grid_gr(ilan,G_vv +k-1)=grid_gr(ilan,G_v  +k-1)
            grid_gr(ilan,G_tt +k-1)=grid_gr(ilan,G_t  +k-1)
            grid_gr(ilan,G_dp +k-1)=grid_gr(ilan,G_dpn+k-1)

          enddo
        enddo

        do i=1,lons_lat
          ilan=i+jlonf
          grid_gr(ilan,G_qm)=grid_gr(ilan,G_q )
          grid_gr(ilan,G_q )=grid_gr(ilan,G_zq)
        enddo

        do k=1,levh
          do i=1,lons_lat
            ilan=i+jlonf
            grid_gr(ilan,G_rm +k-1)=grid_gr(ilan,G_rq +k-1) *filta+
     &     (grid_gr(ilan,G_rm +k-1)+grid_gr(ilan,G_rt +k-1))*filtb
            grid_gr(ilan,G_rq +k-1)=grid_gr(ilan,G_rt +k-1)
          enddo
        enddo

      enddo

      return
      end subroutine do_dynamics_gridfilter

!-------------------------------------------------------------------
      subroutine do_dynamics_gridc2n(grid_gr,
     &                               global_lats_a,lonsperlat)

      real(kind=kind_grid) grid_gr(lonf*lats_node_a_max,lotgr)
      integer,intent(in):: global_lats_a(latg)
      integer,intent(in):: lonsperlat(latg)

      integer   lan,lat,lons_lat,k,i
      integer   jlonf,ilan

!$omp parallel do private(lan)
!$omp+private(lat,lons_lat,jlonf,i,k,ilan)

      do lan=1,lats_node_a
        lat = global_lats_a(ipt_lats_node_a-1+lan)
        lons_lat = lonsperlat(lat)
        jlonf = (lan-1)*lonf
        do k=1,levs
          do i=1,lons_lat
            ilan=i+jlonf
            grid_gr(ilan,G_u  +k-1)=grid_gr(ilan,G_uu +k-1)
            grid_gr(ilan,G_v  +k-1)=grid_gr(ilan,G_vv +k-1)
            grid_gr(ilan,G_t  +k-1)=grid_gr(ilan,G_tt +k-1)
            grid_gr(ilan,G_dpn+k-1)=grid_gr(ilan,G_dp +k-1)
          enddo
        enddo
        do k=1,levh
          do i=1,lons_lat
            ilan=i+jlonf
            grid_gr(ilan,G_rt +k-1)=grid_gr(ilan,G_rq +k-1)
          enddo
        enddo
        do i=1,lons_lat
          ilan=i+jlonf
          grid_gr(ilan,G_zq)=grid_gr(ilan,G_q )
        enddo
      enddo

      return
      end subroutine do_dynamics_gridc2n

!-------------------------------------------------------------------
      subroutine do_dynamics_gridn2c(grid_gr,
     &                               global_lats_a,lonsperlat)

      real(kind=kind_grid) grid_gr(lonf*lats_node_a_max,lotgr)
      integer,intent(in):: global_lats_a(latg)
      integer,intent(in):: lonsperlat(latg)

      integer   lan,lat,lons_lat,k,i
      integer   jlonf,ilan

!$omp parallel do private(lan)
!$omp+private(lat,lons_lat,jlonf,i,k,ilan)

      do lan=1,lats_node_a
        lat = global_lats_a(ipt_lats_node_a-1+lan)
        lons_lat = lonsperlat(lat)
        jlonf = (lan-1)*lonf
        do k=1,levs
          do i=1,lons_lat
            ilan=i+jlonf
            grid_gr(ilan,G_uu +k-1)=grid_gr(ilan,G_u  +k-1)
            grid_gr(ilan,G_vv +k-1)=grid_gr(ilan,G_v  +k-1)
            grid_gr(ilan,G_tt +k-1)=grid_gr(ilan,G_t  +k-1)
            grid_gr(ilan,G_dp +k-1)=grid_gr(ilan,G_dpn+k-1)
          enddo
        enddo
        do k=1,levh
          do i=1,lons_lat
            ilan=i+jlonf
            grid_gr(ilan,G_rq +k-1)=grid_gr(ilan,G_rt +k-1)
          enddo
        enddo
        do i=1,lons_lat
          ilan=i+jlonf
          grid_gr(ilan,G_q )=grid_gr(ilan,G_zq)
        enddo
      enddo

      return
      end subroutine do_dynamics_gridn2c

!-------------------------------------------------------------------
      subroutine do_dynamics_gridn2m(grid_gr,
     &                               global_lats_a,lonsperlat)

      real(kind=kind_grid) grid_gr(lonf*lats_node_a_max,lotgr)
      integer,intent(in):: global_lats_a(latg)
      integer,intent(in):: lonsperlat(latg)

      integer   lan,lat,lons_lat,k,i
      integer   jlonf,ilan

!$omp parallel do private(lan)
!$omp+private(lat,lons_lat,jlonf,i,k,ilan)

      do lan=1,lats_node_a
        lat = global_lats_a(ipt_lats_node_a-1+lan)
        lons_lat = lonsperlat(lat)
        jlonf = (lan-1)*lonf
        do k=1,levs
          do i=1,lons_lat
            ilan=i+jlonf
            grid_gr(ilan,G_uum+k-1)=grid_gr(ilan,G_u  +k-1)
            grid_gr(ilan,G_vvm+k-1)=grid_gr(ilan,G_v  +k-1)
            grid_gr(ilan,G_ttm+k-1)=grid_gr(ilan,G_t  +k-1)
            grid_gr(ilan,G_dpm+k-1)=grid_gr(ilan,G_dpn+k-1)
          enddo
        enddo
        do k=1,levh
          do i=1,lons_lat
            ilan=i+jlonf
            grid_gr(ilan,G_rm +k-1)=grid_gr(ilan,G_rt +k-1)
          enddo
        enddo
        do i=1,lons_lat
          ilan=i+jlonf
          grid_gr(ilan,G_qm)=grid_gr(ilan,G_zq)
        enddo
      enddo

      return
      end subroutine do_dynamics_gridn2m

!-------------------------------------------------------------------
      subroutine do_dynamics_gridap2n(grid_gr,
     &                               global_lats_a,lonsperlat)

      real(kind=kind_grid) grid_gr(lonf*lats_node_a_max,lotgr)
      integer,intent(in):: global_lats_a(latg)
      integer,intent(in):: lonsperlat(latg)

      integer   lan,lat,lons_lat,k,i
      integer   jlonf,ilan

!$omp parallel do private(lan)
!$omp+private(lat,lons_lat,jlonf,i,k,ilan)

      do lan=1,lats_node_a
        lat = global_lats_a(ipt_lats_node_a-1+lan)
        lons_lat = lonsperlat(lat)
        jlonf = (lan-1)*lonf
        do k=1,levs
          do i=1,lons_lat
            ilan=i+jlonf
            grid_gr(ilan,G_u  +k-1)=grid_gr(ilan,G_u  +k-1)+
     &      grid_gr(ilan,G_uup+k-1)-grid_gr(ilan,G_uum+k-1)
            grid_gr(ilan,G_v  +k-1)=grid_gr(ilan,G_v  +k-1)+
     &      grid_gr(ilan,G_vvp+k-1)-grid_gr(ilan,G_vvm+k-1)
            grid_gr(ilan,G_t  +k-1)=grid_gr(ilan,G_t  +k-1)+
     &      grid_gr(ilan,G_ttp+k-1)-grid_gr(ilan,G_ttm+k-1)
            grid_gr(ilan,G_dpn+k-1)=grid_gr(ilan,G_dpn+k-1)+
     &      grid_gr(ilan,G_dpp+k-1)-grid_gr(ilan,G_dpm+k-1)
          enddo
        enddo
        do k=1,levh
          do i=1,lons_lat
            ilan=i+jlonf
            grid_gr(ilan,G_rt +k-1)=grid_gr(ilan,G_rt +k-1)+
     &      grid_gr(ilan,G_rqp+k-1)-grid_gr(ilan,G_rm +k-1)
          enddo
        enddo
        do i=1,lons_lat
          ilan=i+jlonf
          grid_gr(ilan,G_zq )=grid_gr(ilan,G_zq )+
     &    grid_gr(ilan,G_zqp)-grid_gr(ilan,G_qm )
        enddo
      enddo

      return
      end subroutine do_dynamics_gridap2n

!-------------------------------------------------------------------
      subroutine do_dynamics_gridp2n(grid_gr,
     &                               global_lats_a,lonsperlat)

      real(kind=kind_grid) grid_gr(lonf*lats_node_a_max,lotgr)
      integer,intent(in), dimension(latg) :: global_lats_a,lonsperlat

      integer lan,lat,lons_lat,k,i,jlonf,ilan

!$omp parallel do private(lan,lat,lons_lat,jlonf,i,k,ilan)
      do lan=1,lats_node_a
        lat      = global_lats_a(ipt_lats_node_a-1+lan)
        lons_lat = lonsperlat(lat)
        jlonf    = (lan-1)*lonf
        do k=1,levs
          do i=1,lons_lat
            ilan = i + jlonf
            grid_gr(ilan,G_u  +k-1) = grid_gr(ilan,G_uup+k-1)
            grid_gr(ilan,G_v  +k-1) = grid_gr(ilan,G_vvp+k-1)
            grid_gr(ilan,G_t  +k-1) = grid_gr(ilan,G_ttp+k-1)
            grid_gr(ilan,G_dpn+k-1) = grid_gr(ilan,G_dpp+k-1)
          enddo
        enddo
        do k=1,levh
          do i=1,lons_lat
            ilan = i + jlonf
            grid_gr(ilan,G_rt +k-1) = grid_gr(ilan,G_rqp+k-1)
          enddo
        enddo
        do i=1,lons_lat
          ilan = i + jlonf
          grid_gr(ilan,G_zq ) = grid_gr(ilan,G_zqp)
        enddo
      enddo

      return
      end subroutine do_dynamics_gridp2n

!-------------------------------------------------------------------
      subroutine do_dynamics_gridn2p(grid_gr, global_lats_a, lonsperlat)

      real(kind=kind_grid) grid_gr(lonf*lats_node_a_max,lotgr)
      integer,intent(in), dimension(latg) :: global_lats_a, lonsperlat

      integer   lan,lat,lons_lat,k,i
      integer   jlonf,ilan

!$omp parallel do private(lan)
!$omp+private(lat,lons_lat,jlonf,i,k,ilan)

      do lan=1,lats_node_a
        lat      = global_lats_a(ipt_lats_node_a-1+lan)
        lons_lat = lonsperlat(lat)
        jlonf    = (lan-1)*lonf
        do k=1,levs
          do i=1,lons_lat
            ilan = i+jlonf
            grid_gr(ilan,G_uup+k-1) = grid_gr(ilan,G_u  +k-1)
            grid_gr(ilan,G_vvp+k-1) = grid_gr(ilan,G_v  +k-1)
            grid_gr(ilan,G_ttp+k-1) = grid_gr(ilan,G_t  +k-1)
            grid_gr(ilan,G_dpp+k-1) = grid_gr(ilan,G_dpn+k-1)
          enddo
        enddo
        do k=1,levh
          do i=1,lons_lat
            ilan = i+jlonf
            grid_gr(ilan,G_rqp+k-1) = grid_gr(ilan,G_rt +k-1)
          enddo
        enddo
        do i=1,lons_lat
          ilan = i+jlonf
          grid_gr(ilan,G_zqp) = grid_gr(ilan,G_zq)
        enddo
      enddo

      return
      end subroutine do_dynamics_gridn2p

!-------------------------------------------------------------------
      subroutine do_dynamics_gridm2p(grid_gr,
     &                               global_lats_a,lonsperlat)

      real(kind=kind_grid) grid_gr(lonf*lats_node_a_max,lotgr)
      integer,intent(in):: global_lats_a(latg)
      integer,intent(in):: lonsperlat(latg)

      integer   lan,lat,lons_lat,k,i
      integer   jlonf,ilan

!$omp parallel do private(lan)
!$omp+private(lat,lons_lat,jlonf,i,k,ilan)

      do lan=1,lats_node_a
        lat = global_lats_a(ipt_lats_node_a-1+lan)
        lons_lat = lonsperlat(lat)
        jlonf = (lan-1)*lonf
        do k=1,levs
          do i=1,lons_lat
            ilan=i+jlonf
            grid_gr(ilan,G_uup+k-1)=grid_gr(ilan,G_uum+k-1)
            grid_gr(ilan,G_vvp+k-1)=grid_gr(ilan,G_vvm+k-1)
            grid_gr(ilan,G_ttp+k-1)=grid_gr(ilan,G_ttm+k-1)
            grid_gr(ilan,G_dpp+k-1)=grid_gr(ilan,G_dpm+k-1)
          enddo
        enddo
        do k=1,levh
          do i=1,lons_lat
            ilan=i+jlonf
            grid_gr(ilan,G_rqp+k-1)=grid_gr(ilan,G_rm +k-1)
          enddo
        enddo
        do i=1,lons_lat
          ilan=i+jlonf
          grid_gr(ilan,G_zqp)=grid_gr(ilan,G_qm)
        enddo
      enddo

      return
      end subroutine do_dynamics_gridm2p

! -------------------------------------------------------------------
      subroutine do_dynamics_gridupdate(grid_gr,anl_gr_a_2,dt2,
     &                                  global_lats_a,lonsperlat)

      use namelist_dynamics_def

      real(kind=kind_grid) grid_gr(lonf*lats_node_a_max,lotgr)
      real(kind=kind_evod) anl_gr_a_2(lonfx*lota,lats_dim_ext)
      integer,intent(in):: global_lats_a(latg)
      integer,intent(in):: lonsperlat(latg)
      real,   intent(in):: dt2

      integer   lan,lat,lon_dim,lons_lat,k,i
      integer   jlonf,ilan

!$omp parallel do private(lan)
!$omp+private(lat,lon_dim,lons_lat,jlonf,i,k,ilan)


      do lan=1,lats_node_a
        lat = global_lats_a(ipt_lats_node_a-1+lan)
        lon_dim = lon_dims_a(lan)
        lons_lat = lonsperlat(lat)
        jlonf = (lan-1)*lonf
        do k=1,levs
          do i=1,lons_lat
            ilan=i+jlonf
            grid_gr(ilan,G_u  +k-1)=grid_gr(ilan,G_uum+k-1)+
     &      anl_gr_a_2(i+(kau -2+k)*lon_dim,lan)*dt2
            grid_gr(ilan,G_v  +k-1)=grid_gr(ilan,G_vvm+k-1)+
     &      anl_gr_a_2(i+(kav -2+k)*lon_dim,lan)*dt2
            grid_gr(ilan,G_t  +k-1)=grid_gr(ilan,G_ttm+k-1)+
     &      anl_gr_a_2(i+(kat -2+k)*lon_dim,lan)*dt2
            grid_gr(ilan,G_dpn+k-1)=grid_gr(ilan,G_dpm+k-1)+
     &      anl_gr_a_2(i+(kadp-2+k)*lon_dim,lan)*dt2
          enddo
        enddo
        do i=1,lons_lat
            ilan=i+jlonf
            grid_gr(ilan,G_zq)=grid_gr(ilan,G_qm)+
     &      anl_gr_a_2(i+(kaps-1)*lon_dim,lan)*dt2
        enddo
!

        if( .not. ndslfv ) then
        do k=1,levh
          do i=1,lons_lat
            ilan=i+jlonf
            grid_gr(ilan,G_rt +k-1)=grid_gr(ilan,G_rm +k-1)+
     &      anl_gr_a_2(i+(kar-2+k)*lon_dim,lan)*dt2
          enddo
        enddo
        endif
!

      enddo
 
      return
      end subroutine do_dynamics_gridupdate
!
! -------------------------------------------------------------------
!
      subroutine do_dynamics_gridpdpn(grid_gr,global_lats_a,lonsperlat)

      use namelist_dynamics_def

      real(kind=kind_grid) grid_gr(lonf*lats_node_a_max,lotgr)
      integer,intent(in) :: global_lats_a(latg)
      integer,intent(in) :: lonsperlat(latg)

      real(kind=kind_grid)  gtv (ngptc,levs)
      real(kind=kind_grid)  gq  (ngptc)
      real(kind=kind_grid)  prsl(ngptc,levs), dprs(ngptc,levs)

      integer   lan,lat,lons_lat,k,i,lon,jlonf,ilan,njeff

!!$omp parallel do private(lan)
!!$omp+private(lat,lons_lat,jlonf,i,k,ilan)
!!$omp+private(gtv,gq,prsl,dprs)

      do lan=1,lats_node_a
        lat      = global_lats_a(ipt_lats_node_a-1+lan)
        lons_lat = lonsperlat(lat)
        jlonf    = (lan-1)*lonf

!$omp parallel do schedule(dynamic,1) private(lon)
!$omp+private(i,k,ilan,njeff)
!$omp+private(gq,gtv,prsl,dprs)
        do lon=1,lons_lat,ngptc
          njeff = min(ngptc,lons_lat-lon+1)
          ilan  = lon - 1 + jlonf
          do k=1,levs
            do i=1,njeff
              gtv(i,k) = grid_gr(ilan+i, G_t+k-1)
            enddo
          enddo
          do i=1,njeff
            gq(i) = grid_gr(ilan+i,G_zq)
          enddo

          if( gen_coord_hybrid ) then 
            call gch2press(njeff,ngptc,gq, gtv, prsl, dprs)
          else if( hybrid ) then
            call hyb2press(njeff,ngptc,gq, prsl, dprs)
          else
            call sig2press(njeff,ngptc,gq, prsl, dprs)
          endif

          do k=1,levs
            do i=1,njeff
              grid_gr(ilan+i,g_p   +k-1) = prsl(i,k)
              grid_gr(ilan+i,g_dpn +k-1) = dprs(i,k)
            enddo
          enddo

        enddo
      enddo
 
      return
      end subroutine do_dynamics_gridpdpn
! -------------------------------------------------------------------
!
      subroutine do_dynamics_gridpdp(grid_gr,global_lats_a,lonsperlat,
     &                               stp)

      use namelist_dynamics_def

      real(kind=kind_grid) grid_gr(lonf*lats_node_a_max,lotgr)
      integer,intent(in) :: global_lats_a(latg)
      integer,intent(in) :: lonsperlat(latg),stp

      real(kind=kind_grid)  gtv (ngptc,levs)
      real(kind=kind_grid)  gq  (ngptc)
      real(kind=kind_grid)  prsl(ngptc,levs), dprs(ngptc,levs)

      integer lan,lat,lons_lat,k,i,lon,jlonf,ilan,njeff
      integer kt, kq, kp, kdp

      if( stp == 1 ) then
        kt  = g_t
        kq  = g_zq
        kp  = g_p
        kdp = g_dpn
      else
        kt  = g_ttm
        kq  = g_qm
        kp  = g_p
        kdp = g_dpm 
      endif

      do lan=1,lats_node_a
        lat      = global_lats_a(ipt_lats_node_a-1+lan)
        lons_lat = lonsperlat(lat)
        jlonf    = (lan-1)*lonf

!$omp parallel do schedule(dynamic,1) private(lon)
!$omp+private(i,k,ilan,njeff,gq,gtv,prsl,dprs)
        do lon=1,lons_lat,ngptc
          njeff = min(ngptc,lons_lat-lon+1)
          ilan  = lon - 1 + jlonf
          do k=1,levs
            do i=1,njeff
              gtv(i,k) = grid_gr(ilan+i, kt+k-1)
            enddo
          enddo
          do i=1,njeff
            gq(i) = grid_gr(ilan+i,kq)
          enddo

          if( gen_coord_hybrid ) then 
            call gch2press(njeff,ngptc,gq, gtv, prsl, dprs)
          else if( hybrid ) then
            call hyb2press(njeff,ngptc,gq, prsl, dprs)
          else
            call sig2press(njeff,ngptc,gq, prsl, dprs)
          endif

          do k=1,levs
            do i=1,njeff
              grid_gr(ilan+i,kp  +k-1) = prsl(i,k)
              grid_gr(ilan+i,kdp +k-1) = dprs(i,k)
            enddo
          enddo

        enddo
      enddo
 
      return
      end subroutine do_dynamics_gridpdp
! -------------------------------------------------------------------
!
      subroutine do_dynamics_gridp(grid_gr,global_lats_a,lonsperlat,stp)

      use namelist_dynamics_def
      use gfs_dyn_physcons,  rk => con_rocp

      real(kind=kind_grid) grid_gr(lonf*lats_node_a_max,lotgr)
      integer,intent(in):: global_lats_a(latg),stp
      integer,intent(in):: lonsperlat(latg)

      real (kind=kind_grid), parameter :: rk1 = rk + 1.0, rkr = 1.0/rk
     &,                                   R100=100.0, PT01=0.01

      real(kind=kind_grid)  prsi(lonf,levs+1),ppiu(lonf),ppid(lonf),tm

      integer   lan,lat,lons_lat,k,i
      integer   jlonf,ilan,kdp
  
      if( stp.eq.1 ) then
        kdp=g_dpn
      else
        kdp=g_dpm 
      endif

!$omp parallel do private(lan)
!$omp+private(lat,lons_lat,jlonf,i,k,ilan)
!$omp+private(prsi,ppid,ppiu,tm)

      do lan=1,lats_node_a
        lat = global_lats_a(ipt_lats_node_a-1+lan)
        lons_lat = lonsperlat(lat)
        jlonf = (lan-1)*lonf
        prsi(:,levs+1) = 0.0
        do k=levs,1,-1
          do i=1,lons_lat
            ilan=i+jlonf
            prsi(i,k) = prsi(i,k+1) + grid_gr(ilan,kdp+k-1)
          enddo
        enddo

        if( gen_coord_hybrid ) then 
          do k=1,levs
            do i=1,lons_lat
              ilan=i+jlonf
              grid_gr(ilan,g_p+k-1)=0.5*(prsi(i,k)+prsi(i,k+1))
            enddo
          enddo
        else 
          do i=1,lons_lat
            ppid(i) = (prsi(i,1)*PT01) ** rk
          enddo
          do k=1,levs
            do i=1,lons_lat
              ilan=i+jlonf
              ppiu(i)    = (prsi(i,k+1)*PT01) ** rk
              tm         = rk1 * (prsi(i,k) - prsi(i,k+1))
              ppid(i)    = (ppid(i)*prsi(i,k)-ppiu(i)*prsi(i,k+1))
     &                     / tm
              grid_gr(ilan,g_p+k-1)= R100 * ppid(i) ** rkr
              ppid(i)    = ppiu(i)
            enddo
          enddo
        endif

      enddo
 
      return
      end subroutine do_dynamics_gridp
!
! -------------------------------------------------------------------
!
      subroutine do_dynamics_gridzz(grid_gr,
     &                                  global_lats_a,lonsperlat)

      use namelist_dynamics_def
      use gfs_dyn_tracer_const
      use gfs_dyn_physcons, grav  => con_g, rd => con_rd

      real(kind=kind_grid) grid_gr(lonf*lats_node_a_max,lotgr)
      integer,intent(in):: global_lats_a(latg)
      integer,intent(in):: lonsperlat(latg)

      real(kind=kind_grid)  xr   (lonf,levs)
      real(kind=kind_grid)  xcp  (lonf,levs)
      real(kind=kind_grid)  sumq (lonf,levs)
      real(kind=kind_grid)  kh   (lonf,levs)
      real(kind=kind_grid)  phil (lonf,levs)
      real(kind=kind_grid)  prsi (lonf,levs+1)
      real(kind=kind_grid)  phii (lonf),  dphi(lonf)

      integer   lan,lat,lon_dim,lons_lat,k,i,n,kk
      integer   jlonf,ilan,lon,njeff

!!$omp parallel do private(lan)
!!$omp+private(lat,lons_lat,jlonf,i,k,kk,n,ilan)
!!$omp+private(sumq,xr,xcp,kh,prsi,phii,phil,dphi)

      do lan=1,lats_node_a
        lat = global_lats_a(ipt_lats_node_a-1+lan)
        lon_dim = lon_dims_a(lan)
        lons_lat = lonsperlat(lat)
        jlonf = (lan-1)*lonf

!$omp parallel do schedule(dynamic,1) private(lon)
!$omp+private(i,k,ilan,njeff,n,kk)
!$omp+private(xr,xcp,sumq,kh,prsi,phii,phil,dphi)
        do lon=1,lons_lat,ngptc
          njeff = min(ngptc,lons_lat-lon+1)

          if( thermodyn_id == 3 ) then
            sumq = 0.0
            xr   = 0.0
            xcp  = 0.0
            do n=1,ntrac
              if( ri(n) > 0.0 ) then
                do k=1,levs
                  kk = g_rq+k-1+(n-1)*levs
                  do i=1,njeff
                    ilan      = lon-1+i+jlonf
                    xr(i,k)   = xr(i,k)   + grid_gr(ilan,kk) * ri(n)
                    xcp(i,k)  = xcp(i,k)  + grid_gr(ilan,kk) * cpi(n)
                    sumq(i,k) = sumq(i,k) + grid_gr(ilan,kk)
                  enddo
                enddo
              endif
            enddo
            do k=1,levs
              do i=1,njeff
                ilan    = lon-1+i+jlonf
                xr(i,k) = (1.-sumq(i,k))*ri(0)  + xr(i,k)
                xcp(i,k)= (1.-sumq(i,k))*cpi(0) + xcp(i,k)
                kh(i,k) = xr(i,k)/xcp(i,k)*grid_gr(ilan,g_tt+k-1)
              enddo
            enddo
          else

            do k=1,levs
              do i=1,njeff
                ilan    = lon-1+i+jlonf
                kh(i,k) = rd*grid_gr(ilan,g_tt+k-1)
              enddo
            enddo

          endif

          prsi(:,levs+1) = 0.0
          do k=levs,1,-1
            do i=1,njeff
              ilan      = lon-1+i+jlonf
              prsi(i,k) = prsi(i,k+1)+grid_gr(ilan,g_dp+k-1)
            enddo
          enddo

          do i=1,njeff
            ilan    = lon-1+i+jlonf
            phii(i) = grav*grid_gr(ilan,g_gz)
!           phii(i) = 0.0
          enddo
          do k=1,levs
            do i=1,njeff
              ilan     = lon-1+i+jlonf
              dphi(i  ) = (prsi(i,k)-prsi(i,k+1))
     &                  / (prsi(i,k)+prsi(i,k+1)) * kh(i,k)
              phil(i,k) = phii(i  ) + dphi(i)
              phii(i  ) = phil(i,k) + dphi(i)

              grid_gr(ilan,g_zz+k-1)= phil(i,k)
            enddo
          enddo

        enddo

      enddo

      return
      end subroutine do_dynamics_gridzz

! -------------------------------------------------------------------
!
      subroutine do_dynamics_gridzk(grid_gr,
     &                              global_lats_a,lonsperlat,rcs2)

      use namelist_dynamics_def
      use gfs_dyn_tracer_const
      use gfs_dyn_physcons, grav  => con_g, rd => con_rd

      real(kind=kind_grid) grid_gr(lonf*lats_node_a_max,lotgr)
      integer,intent(in):: global_lats_a(latg)
      integer,intent(in):: lonsperlat(latg)
      real(kind=kind_grid) rcs2(latg2)

      real(kind=kind_grid)  xr   (lonf,levs)
      real(kind=kind_grid)  xcp  (lonf,levs)
      real(kind=kind_grid)  sumq (lonf,levs)
      real(kind=kind_grid)  kh   (lonf,levs)
      real(kind=kind_grid)  phil (lonf,levs)
      real(kind=kind_grid)  prsi (lonf,levs+1)
      real(kind=kind_grid)  phii (lonf),  dphi(lonf)

      integer   lan,lat,lon_dim,lons_lat,k,i,n,kk
      integer   jlonf,ilan,lon,njeff
      real(kind=kind_grid)  rcl

!!$omp parallel do private(lan)
!!$omp+private(lat,lons_lat,jlonf,rcl,i,k,kk,n,ilan)
!!$omp+private(sumq,xr,xcp,kh,prsi,phii,phil,dphi)

      do lan=1,lats_node_a
        lat      = global_lats_a(ipt_lats_node_a-1+lan)
        lon_dim  = lon_dims_a(lan)
        lons_lat = lonsperlat(lat)
        jlonf    = (lan-1)*lonf
        rcl      = rcs2(min(lat,latg-lat+1))

!$omp parallel do schedule(dynamic,1) private(lon)
!$omp+private(i,k,ilan,njeff,n,kk)
!$omp+private(xr,xcp,sumq,kh,prsi,phii,phil,dphi)
        do lon=1,lons_lat,ngptc
          njeff = min(ngptc,lons_lat-lon+1)

          if( thermodyn_id == 3 ) then
            sumq = 0.0
            xr   = 0.0
            xcp  = 0.0
            do n=1,ntrac
              if( ri(n) > 0.0 ) then
                do k=1,levs
                  kk = g_rq+k-1+(n-1)*levs
                  do i=1,njeff
                    ilan      = lon-1+i+jlonf
                    xr(i,k)   = xr(i,k)   + grid_gr(ilan,kk) * ri(n)
                    xcp(i,k)  = xcp(i,k)  + grid_gr(ilan,kk) * cpi(n)
                    sumq(i,k) = sumq(i,k) + grid_gr(ilan,kk)
                  enddo
                enddo
              endif
            enddo
            do k=1,levs
              do i=1,njeff
                ilan     = lon-1+i+jlonf
                xr(i,k)  = (1.-sumq(i,k))*ri(0)  + xr(i,k)
                xcp(i,k) = (1.-sumq(i,k))*cpi(0) + xcp(i,k)
                kh(i,k)  = xr(i,k)/xcp(i,k)*grid_gr(ilan,g_tt+k-1)
              enddo
            enddo

          else

            do k=1,levs
              do i=1,njeff
                ilan    = lon-1+i+jlonf
                kh(i,k) = rd*grid_gr(ilan,g_tt+k-1)
              enddo
            enddo

          endif

          prsi(:,levs+1) = 0.0
          do k=levs,1,-1
            do i=1,njeff
              ilan    = lon-1+i+jlonf
              prsi(i,k) = prsi(i,k+1)+grid_gr(ilan,g_dp+k-1)
            enddo
          enddo

          do i=1,njeff
            ilan    = lon-1+i+jlonf
            phii(i) = grav*grid_gr(ilan,g_gz)
!           phii(i) = 0.0
          enddo
          do k=1,levs
            do i=1,njeff
              ilan      = lon-1+i+jlonf
              dphi(i  ) = (prsi(i,k)-prsi(i,k+1))
     &                  / (prsi(i,k)+prsi(i,k+1)) * kh(i,k)
              phil(i,k) = phii(i  ) + dphi(i)
              phii(i  ) = phil(i,k) + dphi(i)
              grid_gr(ilan,g_zz+k-1)= phil(i,k)
     &                       + 0.5*(grid_gr(ilan,g_uu+k-1)**2 +
     &                              grid_gr(ilan,g_vv+k-1)**2)*rcl
            enddo
          enddo
        enddo

      enddo

      return
      end subroutine do_dynamics_gridzk

! -------------------------------------------------------------------
      subroutine do_dynamics_gridspdmax(grid_gr,spdmax,rcs2,
     &                                 global_lats_a,lonsperlat)

      real(kind=kind_grid) grid_gr(lonf*lats_node_a_max,lotgr)
      real(kind=kind_grid) rcs2(latg2)
      real(kind=kind_evod) spdmax(levs)
      integer,intent(in):: global_lats_a(latg)
      integer,intent(in):: lonsperlat(latg)

      real      spd,spdlat,rcl
      integer   lan,lat,lons_lat,i,k,ilan,jlonf

      do lan=1,lats_node_a
        lat = global_lats_a(ipt_lats_node_a-1+lan)
        lons_lat = lonsperlat(lat)
        jlonf = (lan-1)*lonf
        rcl  = rcs2(min(lat,latg-lat+1))
        do k=1,levs
          spdlat = 0.0
          do i=1,lons_lat
            ilan=i+jlonf
            spd = grid_gr(ilan,g_u +k-1)**2.+grid_gr(ilan,g_v +k-1)**2.
            spd = spd * rcl
            spdlat=max(spdlat,spd)
          enddo
          spdmax(k) = max( spdlat, spdmax(k) )
        enddo
      enddo
 
      return
      end subroutine do_dynamics_gridspdmax
!
! -------------------------------------------------------------------
      subroutine do_dynamics_gridmean(grid,glbmean,wgt_a,lats_nodes_a,
     &                                 global_lats_a,lonsperlat,km)

      real(kind=kind_grid) grid(lonf*lats_node_a_max,km)
      real(kind=kind_grid) glbmean(km)
      real(kind=kind_dbl_prec) wgt_a(latg)
      integer,intent(in):: lats_nodes_a(nodes),km
      integer,intent(in):: global_lats_a(latg)
      integer,intent(in):: lonsperlat(latg)

      real(kind=kind_grid) tmps(lats_node_a_max,km,nodes)
      real(kind=kind_grid) tmpr(lats_node_a_max,km,nodes)
      real(kind=kind_grid) gmeanj(latg,km)
      real sum
      integer   lan,lat,lons_lat,i,k,node
      integer   ilan,jlonf,ilat,nsend,ierr

      do lan=1,lats_node_a
        lat = global_lats_a(ipt_lats_node_a-1+lan)
        lons_lat = lonsperlat(lat)
        jlonf = (lan-1)*lonf
        do k=1,km
          sum = 0.0
          do i=1,lons_lat
            ilan=i+jlonf
            sum = sum + grid(ilan,k)
          enddo
          gmeanj(lan,k) = sum / (2*lons_lat)
        enddo
      enddo
      if( nodes.ne.1 ) then
        do node=1,nodes
          do lan=1,lats_node_a
            tmps(lan,1:km,node) = gmeanj(lan,1:km)
          enddo
        enddo
        nsend = lats_node_a_max * km
        call mpi_alltoall(tmps,nsend,mpi_a_def,
     &                    tmpr,nsend,mpi_a_def,MC_COMP,ierr)
        ilat=1
        do node=1,nodes
          do lan=1,lats_nodes_a(node)
            lat=global_lats_a(ilat)
            gmeanj(lat,1:km)=tmpr(lan,1:km,node)
            ilat=ilat+1
          enddo
        enddo
      endif
      do k=1,km
        glbmean(k) = 0
        do lat=1,latg
          glbmean(k) = glbmean(k) +
     &                 wgt_a(min(lat,latg-lat+1))*gmeanj(lat,k)
        enddo
      enddo

      return
      end subroutine do_dynamics_gridmean
!
! -------------------------------------------------------------------
      subroutine do_dynamics_gridcheck(grid_gr,
     &                                 global_lats_a,lonsperlat,chr)

      real(kind=kind_grid) grid_gr(lonf*lats_node_a_max,lotgr)
      integer,intent(in):: global_lats_a(latg)
      integer,intent(in):: lonsperlat(latg)
      character*(*) chr

      integer   lan,lat,lons_lat,k

      print *,' check: g_ttm g_tt g_t ',g_ttm,g_tt,g_t
      do lan=1,lats_node_a
        lat = global_lats_a(ipt_lats_node_a-1+lan)
        lons_lat = lonsperlat(lat)
        print *,' gridcheck: lan lat lons_lat ',lan,lat,lons_lat
        do k=1,levs
          print *,' check grid of ttm tt t at k=',k,' lat=',lat
          call mymaxmin(grid_gr(1,g_ttm+k-1),lons_lat,lonf,1,chr)
          call mymaxmin(grid_gr(1,g_tt +k-1),lons_lat,lonf,1,chr)
          call mymaxmin(grid_gr(1,g_t  +k-1),lons_lat,lonf,1,chr)
        enddo
      enddo
 
      return
      end subroutine do_dynamics_gridcheck
!
      end module do_dynamics_mod
