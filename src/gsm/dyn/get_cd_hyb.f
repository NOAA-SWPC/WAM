      subroutine get_cd_hyb(dt)
      use gfs_dyn_machine , only : kind_grid
      use gfs_dyn_resol_def
      use gfs_dyn_coordinate_def                   ! hmhj
      use gfs_dyn_layout1, only : me

      implicit none

      real(kind=kind_evod), intent(in) ::  dt
!
      real(kind=kind_evod), dimension(levs,levs) ::  ym, rim
!sela real(kind=kind_evod) dm205_hyb(jcap1,levs,levs)
      real(kind=kind_evod) ddd(jcap1),ppp(jcap1),rrr(jcap1)
      real(kind=kind_evod), parameter :: cons0=0.d0, cons1=1.d0
      real(kind=kind_evod) rnn1
      integer              lu(levs),mu(levs)
      integer              i,j,k,n,nn
 
      call am_bm_hyb
 
      do k=1,levs
        do j=1,levs
          rim(j,k) = cons0
        enddo
      enddo
      do k=1,levs
        rim(k,k) = cons1
      enddo
 
!***********************************************************************
!
!       initialisation of D_HYB_m.
!
!***********************************************************************
!
!     computations which do not depend on n
!     *************************************
!
!-------------------------------------------------------
      do i=1,levs
        do j=1,levs
          ym(i,j) = tor_hyb(i)*svhyb(j)
        enddo
        do k=1,levs
          do j=1,levs
            ym(i,j) = ym(i,j) + amhyb(i,k)*bmhyb(k,j)
          enddo
        enddo
      enddo
!-------------------------------------------------------
!
!     computations which on n
!     ***********************
!..................................................................
      do nn=1,jcap1
        n = nn-1
        rnn1 = n*(n+1)
        do i=1,levs
          do j=1,levs
            dm205_hyb(nn,i,j) = rim(i,j) + rnn1*dt*dt*ym(i,j)
          enddo
        enddo
      enddo
!..................................................................
      call matinv(dm205_hyb,jcap1,levs,ddd,ppp,rrr)
      do nn=1,jcap1
        do i=1,levs
          do j=1,levs
            D_HYB_m(i,j,nn) = dm205_hyb(nn,i,j)
          enddo
        enddo
      enddo
      if (me == 0) then
        write(0,*)' completed hyb sicdif preparation getcd_hyb dt=',dt
      endif
      return
      end
