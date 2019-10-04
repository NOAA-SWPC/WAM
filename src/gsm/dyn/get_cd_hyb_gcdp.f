      subroutine get_cd_hyb_gcdp(dti)
!
! program log:
! 
! 20110220    Henry Juang initiated and wrote the code for mass_dp and ndslfv
! 20130930    Henry Juang add option of uncenter for verticoord_id=3
!
      use gfs_dyn_machine , only : kind_grid
      use gfs_dyn_resol_def
      use gfs_dyn_coordinate_def
      use namelist_dynamics_def , only : lsidea
      implicit none
      integer              i,j,k,n,nn
      real(kind=kind_evod) dti,dt,rnn1
      real(kind=kind_evod) ym(levs,levs)
      real(kind=kind_evod) rim(levs,levs)
      real(kind=kind_evod) ddd(jcap1),ppp(jcap1),rrr(jcap1)
      integer              lu(levs),mu(levs)
      real(kind=kind_evod) cons0,cons1     !constant

!     print *,' enter get_cd_hyb_gcdp ',dti

      cons0 = 0.d0     !constant
      cons1 = 1.d0     !constant

! hmhj forward-weighted semi-implicit if eps_si >0
!      center-averaged semi-implicit if eps_si=0
!
      if( vertcoord_id.eq.3. .or. lsidea ) then
        eps_si=0.20
      else
        eps_si=0.00
      endif
!hmhj debug
      eps_si=0.5
      dt=(cons1+eps_si)*dti
 
      call am_bm_hyb_gcdp
 
      do 250 k=1,levs
      do 200 j=1,levs
      rim(j,k)=cons0     !constant
200   continue
250   continue
 
      do 1 k=1,levs
      rim(k,k) = cons1     !constant
1     continue
 
c***********************************************************************
c
c       initialisation of D_HYB_m.
c
c***********************************************************************
c
c     computations which do not depend on n
c     *************************************
c
!-------------------------------------------------------
      ym = 0.0
      do 10 j=1,levs
 
      do  k=1,levs
        do  i=1,levs
          ym(i,j) = ym(i,j) + hmhyb(i,k)*smhyb(k,j)
        enddo
      enddo
 
      do k=1,levs
        do i=1,levs
          ym(i,j) = ym(i,j) + amhyb(i,k)*bmhyb(k,j)
        enddo
      enddo
 
10    continue
!-------------------------------------------------------
c
c     computations which on n
c     ***********************
!..................................................................
      do 2000 nn=1,jcap1
 
       n = nn-1
       rnn1 =       n*(n+1)
 
       do 14 i=1,levs
       do 13 j=1,levs
        dm205_hyb(nn,i,j) = rim(i,j) + rnn1*dt*dt*ym(i,j)
13     continue
14     continue
 
2000  continue
!..................................................................
      call matinv(dm205_hyb,jcap1,levs,ddd,ppp,rrr)
      do 23 nn=1,jcap1
      do 22 i=1,levs
      do 21 j=1,levs
      D_HYB_m(i,j,nn)=dm205_hyb(nn,i,j)
21    continue
22    continue
23    continue
!hmhj print 100,dt
100   format(1h ,'completed hyb sicdif preparation getcd_hyb dt=',f7.1)

!     print *,' end of get_cd_hyb_gcdp '

      return
      end
