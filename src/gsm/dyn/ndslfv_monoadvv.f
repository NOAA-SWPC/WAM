      subroutine ndslfv_monoadvv (grid_gr,pdot,
     &                           global_lats_a,lonsperlat,deltim)
!
! a routine to do non-iteration semi-Lagrangain advection
! considering advection  with monotonicity in interpolation
! contact: hann-ming henry juang
! program log:
! 2011 02 20 : henry juang, initial implemented into nems as NDSL with mass_dp
! 2013 09 30 : henry juang, add option of theta advection, (used later)
!
!
      use gfs_dyn_machine , only : kind_grid
      use gfs_dyn_resol_def
      use gfs_dyn_layout1
      use gfs_dyn_vert_def
      use gfs_dyn_coordinate_def
      use gfs_dyn_physcons
      use gfs_dyn_tracer_const
      use gfs_dyn_mpi_def
      implicit none

      real(kind=kind_grid) grid_gr(lonf*lats_node_a_max,lotgr)
      real(kind=kind_grid) pdot(lonf,levs+1,lats_node_a_max)
      real(kind=kind_grid) plev(lonfull,levs+1)
      integer,intent(in):: global_lats_a(latg)
      integer,intent(in):: lonsperlat(latg)
      real,   intent(in):: deltim
   
      real      qqlon(lonfull,levs*ndslvvar,latpart)

!     real      xksav(lonfull,levs      ,latpart)
!     real      stsav(lonfull,levs      ,latpart)
!     real      ttsav(lonfull,levs      ,latpart)

!     real      xr    (lonfull,levs)
!     real      xcp   (lonfull,levs)
!     real      sumrq (lonfull,levs)
!     real      xkappa(lonfull,levs)
!     real      kappa, pi, ply, hh
      real      cons0, cons1

      logical 	lprint

      integer mono,mass,nvars
      integer ilan,i,n,k,kk,lon,lan,lat,lons_lat,lon_dim,jlonf
      integer kss, kqq, ktt, kuu, kvv, nqq
      integer ks, kp , kq , kt , ku , kv
      integer kpg, kqg, ktg, kug, kvg
!
!     lprint = .false.

!     if( lprint ) then
!       print *,' enter ndslfv_advect  with monotonicity '
!     endif
!
      mono  = 1
      mass  = 0
      cons0 = 0.0
      cons1 = 1.0
!
      kuu = 1
      kvv = kuu + levs
      ktt = kvv + levs
      kqq = ktt + levs
      
      nvars = ndslvvar
!
!$omp parallel do schedule(dynamic,1) private(lan)
!$omp+private(lat,lons_lat,jlonf,i,k,plev,ilan,mass)
!$omp+private(kug,kvg,ktg,kpg,kqg)
!$omp+private(ku ,kv ,kt ,kp ,kq )

      do lan=1,lats_node_a

        lat = global_lats_a(ipt_lats_node_a-1+lan)
        lons_lat = lonsperlat(lat)
        jlonf = (lan-1)*lonf
!
        plev(:,levs+1) = 0.0
        do k=levs,1,-1
          kpg=g_dp+k-1
          do i=1,lons_lat
            ilan=i+jlonf
            plev(i,k)=plev(i,k+1)+grid_gr(ilan,kpg)
          enddo
        enddo

!       if( lprint ) then
!       do k=1,levs
!       print *,' k= ',k
!       call mymaxmin(pdot(1,k,lan),lons_lat,lonfull,1,' pdot ')
!       call mymaxmin(plev(1,k),lons_lat,lonfull,1,' plev ')
!       enddo
!       endif
!
! u v h at n+1* 
!       kappa = con_rd / con_cp
        do k=1,levs
          ku=kuu+k-1
          kv=kvv+k-1
          kt=ktt+k-1
          kug=g_u+k-1
          kvg=g_v+k-1
          ktg=g_t+k-1
          do i=1,lons_lat
            ilan=i+jlonf
            qqlon(i,ku,lan) = grid_gr(ilan,kug) 
            qqlon(i,kv,lan) = grid_gr(ilan,kvg) 
            qqlon(i,kt,lan) = grid_gr(ilan,ktg) 
          enddo
        enddo
! rq at n+1* 
        do k=1,levh
          kq=kqq+k-1
          kqg=g_rt+k-1
          do i=1,lons_lat
            ilan=i+jlonf
            qqlon(i,kq,lan) = grid_gr(ilan,kqg)
          enddo
        enddo
!
! ----- prepare xr, xcp, xkapa 
!
!       xr    = cons0
!       xcp   = cons0
!       sumrq = cons0
!       do n=1,ntrac
!         nqq = kqq + (n-1)*levs
!         if( ri(n) .ne. cons0 .and. cpi(n) .ne. cons0 ) then
!           do k=1,levs
!             kq=nqq+k-1
!             do i=1,lons_lat
!               xr   (i,k) = xr   (i,k) + qqlon(i,kq,lan)*ri(n)
!               xcp  (i,k) = xcp  (i,k) + qqlon(i,kq,lan)*cpi(n)
!               sumrq(i,k) = sumrq(i,k) + qqlon(i,kq,lan)
!             enddo
!           enddo
!         endif
!       enddo
!       do k=1,levs
!         do i=1,lons_lat
!           xr (i,k)   = ( cons1 - sumrq(i,k) )*ri(0)  + xr (i,k)
!           xcp(i,k)   = ( cons1 - sumrq(i,k) )*cpi(0) + xcp(i,k)
!           xkappa(i,k) = xr(i,k) / xcp(i,k)
!           xksav(i,k,lan) = xkappa(i,k)
!         enddo
!       enddo
!
! change h to theta
!
!       do k=1,levs
!         kt=ktt+k-1
!         do i=1,lons_lat
!           pi = ((plev(i,k)+plev(i,k+1))*0.005)**xkappa(i,k)
!           ttsav(i,k ,lan) = qqlon(i,kt,lan)
!           qqlon(i,kt,lan) = qqlon(i,kt,lan)/pi
!           stsav(i,k ,lan) = qqlon(i,kt,lan)
!         enddo
!       enddo
!
!
        mass=0
        call vertical_cell_advect (lons_lat,lonfull,levs,nvars,
     &            deltim,plev,pdot(1,1,lan),qqlon(1,1,lan),mass)
!
! dp with mass conserving
!       do k=1,levs
!         kp =g_dp +k-1
!         kpg=g_dpn+k-1
!         do i=1,lon_dim
!           ilan=i+jlonf
!           rrlon(i,k,lan) = grid_gr(ilan,kpg)/grid_gr(ilan,kp )
!           rrlon(i,k,lan) = 1.0
!         enddo
!       enddo
!       mass=1
!       call vertical_cell_advect (lons_lat,lonfull,levs,1,
!    &            deltim,plev,pdot(1,1,lan),rrlon(1,1,lan),mass)
!       do k=1,levs
!         kp =g_dp +k-1
!         kpg=g_dpn+k-1
!         do i=1,lon_dim
!           ilan=i+jlonf
!           grid_gr(ilan,kpg) = rrlon(i,k ,lan)*grid_gr(ilan,kp )
!         enddo
!       enddo
!
! ----- prepare xr, xcp, xkapa 
!
!       xr    = cons0
!       xcp   = cons0
!       sumrq = cons0
!       do n=1,ntrac
!         nqq = kqq + (n-1)*levs
!         if( ri(n) .ne. cons0 .and. cpi(n) .ne. cons0 ) then
!           do k=1,levs
!             kq=nqq+k-1
!             do i=1,lons_lat
!               xr   (i,k) = xr   (i,k) + qqlon(i,kq,lan)*ri(n)
!               xcp  (i,k) = xcp  (i,k) + qqlon(i,kq,lan)*cpi(n)
!               sumrq(i,k) = sumrq(i,k) + qqlon(i,kq,lan)
!             enddo
!           enddo
!         endif
!       enddo
!       do k=1,levs
!         do i=1,lons_lat
!           xr (i,k)   = ( cons1 - sumrq(i,k) )*ri(0)  + xr (i,k)
!           xcp(i,k)   = ( cons1 - sumrq(i,k) )*cpi(0) + xcp(i,k)
!           xkappa(i,k) = xr(i,k) / xcp(i,k)
!         enddo
!       enddo
!
! add change of p, theta and kappa to h
!
!       do k=1,levs
!         kt=ktt+k-1
!         do i=1,lons_lat
!           ply = (plev(i,k)+plev(i,k+1))*0.5
!           pi = ((plev(i,k)+plev(i,k+1))*0.005)**xkappa(i,k)
!           hh = ttsav(i,k,lan)
!           ttsav(i,k,lan) = ttsav(i,k,lan)
!    &                      + xksav(i,k,lan)*hh/ply*
!    &                       (pdot(i,k,lan)+pdot(i,k+1,lan))*deltim
!           ttsav(i,k,lan) = ttsav(i,k,lan) 
!    &                      + pi*(qqlon(i,kt,lan)-stsav(i,k,lan))
!    &                      + hh*log(ply/100.)*
!    &                        (xkappa(i,k)-xksav(i,k,lan))
!           qqlon(i,kt,lan)= ttsav(i,k,lan)
!         enddo
!       enddo
!
! u v h at n+1
        do k=1,levs
          ku=kuu+k-1
          kv=kvv+k-1
          kt=ktt+k-1
          kug=g_u+k-1
          kvg=g_v+k-1
          ktg=g_t+k-1
          do i=1,lons_lat
            ilan=i+jlonf
            grid_gr(ilan,kug) = qqlon(i,ku,lan)
            grid_gr(ilan,kvg) = qqlon(i,kv,lan)
            grid_gr(ilan,ktg) = qqlon(i,kt,lan)
          enddo
        enddo
! rq at n+1
        do k=1,levh
          kq=kqq+k-1
          kqg=g_rt+k-1
          do i=1,lons_lat
            ilan=i+jlonf
            grid_gr(ilan,kqg) = qqlon(i,kq,lan)
          enddo
        enddo
!
!       if( lprint ) then
!       call mymaxmin(qqlon(1,kqq,lan),lons_lat,lonfull,1,' q vertadv')
!       print *,' ------------------------------------------- '
!       print *,' finish updating n+1* in grid_gr at lan=',lan
!       endif
!
      enddo

! 
! ===============================
!
      return
      end
