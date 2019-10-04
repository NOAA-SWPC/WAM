      subroutine gfidi_gchdp_noadv_noq (lon_dim,lons_lat,lat,
     &  dg,zg,hg,ug,vg,rqg,dpg, psg,zsphi,zslam, 
     &  dpphi,dplam,zzphi,zzlam,
     &  rcl,spdmax,deltim,
     &  dhdf,dhdl,dudl,dvdl,dudf,dvdf,
     &  dpsdt,ddpdt,dhdt,dudt,dvdt,dppdt,pdot)
 
!
! version of ppi=ak+bk*psfc+ck*(h/h0)^(1/kappa)
! pressure gradient by gradients of pressure and geopotential
! this is modified hybrid by finite difference from hann-ming henry juang 
! hmhj : use h=CpT instead of Tv for thermodynamical equation
! hmhj : use pressure thickness as prognostic variables
!        use geopotential gradient with all gas-tracer effect
!        coordinate definition used only for vertical flux
! 2011 02 20 :  Henry jaung, no tracer advection for NDSL
!
 
      use gfs_dyn_machine , only : kind_grid
 
      use gfs_dyn_resol_def
      use namelist_dynamics_def
      use gfs_dyn_coordinate_def
      use gfs_dyn_tracer_const
      use gfs_dyn_mpi_def
      use gfs_dyn_physcons, rerth => con_rerth
     &,             rd => con_rd, cpd => con_cp
     &,             omega => con_omega, kappa => con_rocp
     &,             grav  => con_g
      implicit none

      real (kind=kind_grid), parameter :: rkappa=1.0/kappa
      integer lon_dim,lons_lat
      integer i,k,kk,n,lat
      real coriol,rcl,sinra,deltim,sinlat,absvor
      real det,hkrt0,wmkm1,wmkp1,rkt
      real
     1    dg(lon_dim,levs), zg(lon_dim,levs), hg(lon_dim,levs),  
     2    ug(lon_dim,levs), vg(lon_dim,levs), dpg(lon_dim,levs),
     2   rqg(lon_dim,levs,ntrac)
      real  psg(lon_dim), zsphi(lon_dim), zslam(lon_dim)
      real
     1  dpphi(lon_dim,levs),      zzphi(lon_dim,levs),
     1  dplam(lon_dim,levs),      zzlam(lon_dim,levs)
      real
     1  dhdf(lon_dim,levs),       dhdl(lon_dim,levs),
     1  dudf(lon_dim,levs),       dudl(lon_dim,levs),
     1  dvdf(lon_dim,levs),       dvdl(lon_dim,levs)
      real
     1  dudt(lon_dim,levs),       dvdt(lon_dim,levs),
     1  dpsdt(lon_dim),          ddpdt(lon_dim,levs),
     1  dhdt(lon_dim,levs),      dppdt(lon_dim,levs),
     1  spdmax(levs)
      real pdot(lonf,levs+1)
!

! -----------------------------
c
      real rdelp2(lons_lat,levs)
      real xr(lons_lat,levs),xcp(lons_lat,levs)
      real xkappa(lons_lat,levs)
      real xkappai(lons_lat,levs+1)
      real cg(lons_lat,levs),ek(lons_lat,levs)
      real fb(lons_lat,levs+1),fg(lons_lat,levs)
      real dpxi (lons_lat,levs+1),dpyi (lons_lat,levs+1)
      real dpxii(lons_lat,levs+1),dpyii(lons_lat,levs+1)
      real ppi(lons_lat,levs+1), ppl(lons_lat,levs)
      real hki (lons_lat,levs+1),hkci(lons_lat,levs+1)
      real dpp(lons_lat,levs),rpp(lons_lat,levs)
      real dlnpx(lons_lat,levs),dlnpy(lons_lat,levs)
      real wf(lons_lat,levs+1),wf1(lons_lat,levs+1)
      real wflx(lons_lat,levs+1)
      real wml(lons_lat,levs),wmm(lons_lat,levs),wmu(lons_lat,levs)
      real work(lons_lat,levs)
      real dup(lons_lat,levs),dum(lons_lat,levs)
      real alpha(lons_lat,levs),betta(lons_lat,levs)
      real gamma(lons_lat,levs),delta(lons_lat,levs), alnpk
      real zadv(lons_lat,levs,3+ntrac)
      real sumdf(lons_lat,levs), sumdl(lons_lat,levs)
      real sumrq(lons_lat,levs)
!
      real, parameter ::  cons0=0.0, cons1=1.0, cons2=2.0, cons0p5=0.5
      real rdt2
      integer levsb
      integer irc
      logical lprint

      lprint = .false.
!
! -------- prepare coriolis and gaussian weighting 
!
      rdt2   = cons1/(cons2*deltim)

      sinra  = sqrt(cons1-cons1/rcl)
      coriol = cons2*omega*sinra
      sinra  = sinra/rerth
                                                                                
      if(lat.gt.latg2) then
        coriol = -coriol
        sinra  = -sinra
      endif
!
! max wind
!
      spdmax = cons0
      do k=1,levs
       do i=1,lons_lat
        ek(i,k)   = (  ug(i,k)*ug(i,k)+vg(i,k)*vg(i,k) ) * rcl
        spdmax(k) = max( ek(i,k),spdmax(k) )
       enddo
      enddo
!
! ----- prepare xr, xcp, xkapa 
!
      xr    = cons0
      xcp   = cons0
      sumrq = cons0
!
      do n=1,ntrac
        if( ri(n) .ne. cons0 .and. cpi(n) .ne. cons0 ) then
          do k=1,levs
            do i=1,lons_lat
              xr   (i,k) = xr   (i,k) + rqg  (i,k,n)*ri(n)
              xcp  (i,k) = xcp  (i,k) + rqg  (i,k,n)*cpi(n)
              sumrq(i,k) = sumrq(i,k) + rqg  (i,k,n)
            enddo
          enddo
        endif
      enddo
      do k=1,levs
        do i=1,lons_lat
          xr (i,k)   = ( cons1 - sumrq(i,k) )*ri(0)  + xr (i,k)
          xcp(i,k)   = ( cons1 - sumrq(i,k) )*cpi(0) + xcp(i,k)
          xkappa(i,k) = xr(i,k) / xcp(i,k)
        enddo
      enddo
!
! ----- compute ppi, ppl, rpp, and dpp from dpg
!
      ppi(:,levs+1) = 0.0
      do k=levs,1,-1
        do i=1,lons_lat
          ppi(i,k) = ppi(i,k+1) + dpg(i,k)
        enddo
      enddo
!
      do k=1,levs
        do i=1,lons_lat
          ppl(i,k)  = cons0p5 * ( ppi(i,k) + ppi(i,k+1) )
          rpp(i,k)=1./(ppi(i,k)+ppi(i,k+1))
          dpp(i,k)=ppi(i,k)-ppi(i,k+1)
          rdelp2(i,k) = cons0p5/dpp(i,k)
          if( dpp(i,k) .lt. cons0 ) then
            print *,' ----- dpp < 0 in gfidi_hyb_gchdp at i k ',i,k
          endif
        enddo
      enddo
!
! ----------------------------------------------------------------------
! ----- prepare dpxi, dpyi
!
      dpxi (:,levs+1) = 0.0
      dpyi (:,levs+1) = 0.0
      do k=levs,1,-1
        do i=1,lons_lat
          dpxi (i,k) = dpxi (i,k+1) + dplam(i,k)
          dpyi (i,k) = dpyi (i,k+1) + dpphi(i,k)
        enddo
      enddo
      do k=1,levs
        do i=1,lons_lat
          dpxi (i,k)=dpxi (i,k)*rcl
          dpyi (i,k)=dpyi (i,k)*rcl
        enddo
      enddo
c

      rkt = grav / ( rd * 300. )
      do k=1,levs
        do i=1,lons_lat
          fg(i,k)=-ddpdt(i,k)
     &           +dpp(i,k)*dg(i,k)
          cg(i,k)=-dppdt(i,k)
        enddo
      enddo
c
      do i=1,lons_lat
        fb(i,levs+1)=cons0
      enddo
      do k=levs,1,-1
        do i=1,lons_lat
          fb(i,k)=fb(i,k+1)+fg(i,k)
        enddo
      enddo
c
c local change of surface pressure  d ps dt
c
      do i=1,lons_lat
        dpsdt(i) = - fb(i,1)
      enddo

c
c get dlnpx dlnpy from dp
c
      do k=1,levs
        do i=1,lons_lat
          dlnpx(i,k)=rpp(i,k)*(dpxi(i,k)+dpxi(i,k+1))/rcl
          dlnpy(i,k)=rpp(i,k)*(dpyi(i,k)+dpyi(i,k+1))/rcl
        enddo
      enddo
! --------------- horizontal advection first --------
c
c horizontal advection for all
c
!     do k=1,levs
!       do i=1,lons_lat
!         dudt(i,k)=
!    &               -ug(i,k)*dudl(i,k)-vg(i,k)*dudf(i,k)
!         dvdt(i,k)=
!    &               -ug(i,k)*dvdl(i,k)-vg(i,k)*dvdf(i,k)
!         dhdt(i,k)=
!    &               -ug(i,k)*dhdl(i,k)-vg(i,k)*dhdf(i,k)
!       enddo
!     enddo

!     if( lprint ) then
!     print *,' input lat=',lat
!     k=levs-1
!     call mymaxmin(dudt(1,k),lons_lat,lons_lat,1,' advh dudt ')
!     call mymaxmin(dvdt(1,k),lons_lat,lons_lat,1,' advh dvdt ')
!     call mymaxmin(dhdt(1,k),lons_lat,lons_lat,1,' advh dhdt ')
!     endif
c
c
c total derivative of horizontal wind
c
!

      do k=1,levs
       do i=1,lons_lat
         dudt(i,k)=        dudt(i,k)
     &                   - xkappa(i,k) *hg(i,k)*dlnpx(i,k)
     &                   - zzlam(i,k) 	
     &                   + vg(i,k)*coriol
         dvdt(i,k)=        dvdt(i,k)
     &                   - xkappa(i,k) *hg(i,k)*dlnpy(i,k)
     &                   - zzphi(i,k) 	
     &                   - ug(i,k)*coriol
     &                   - ek(i,k) * sinra
       enddo
      enddo

c
c total derivative of virtual temperature
c

      do k=1,levs
       do i=1,lons_lat
         dhdt(i,k)=  dhdt(i,k)
     &              +xkappa(i,k)*hg(i,k)*rpp(i,k)*
     &                       (cg(i,k)-fb(i,k)-fb(i,k+1))
       enddo
      enddo
c
      if( lprint ) then
      print *,' lat=',lat
      k=levs-1
      call mymaxmin(dudt(1,k),lons_lat,lons_lat,1,' tall dudt ')
      call mymaxmin(dvdt(1,k),lons_lat,lons_lat,1,' tall dvdt ')
      call mymaxmin(dhdt(1,k),lons_lat,lons_lat,1,' tall dhdt ')
      endif
c
!
! ------ hybrid to solve vertical flux ----------
!
!     alpha(:,levs)=cons0
!     do k=2,levs
!       do i=1,lons_lat
!         alpha(i,k-1)=(ppl(i,k-1)/100.)**xkappa(i,k-1)
!         alpha(i,k-1)=(ppl(i,k  )/100.)**xkappa(i,k)/alpha(i,k-1)
!       enddo
!     enddo
!     betta(:,1   )=cons0
!     do k=1,levs-1
!       do i=1,lons_lat
!         betta(i,k+1)=(ppl(i,k+1)/100.)**xkappa(i,k+1)
!         betta(i,k+1)=(ppl(i,k  )/100.)**xkappa(i,k)/betta(i,k+1)
!       enddo
!     enddo
!     gamma(:,1   )=cons0
!     do k=2,levs
!       do i=1,lons_lat
!         alnpk = log( ppl(i,k) / 100. )        ! cb
!         gamma(i,k)=cons1-(xkappa(i,k-1)-xkappa(i,k))*alnpk
!    &                    -xkappa(i,k)*dpp(i,k)*rpp(i,k)*cons2
!       enddo
!     enddo
!     delta(:,levs)=cons0
!     do k=1,levs-1
!       do i=1,lons_lat
!         alnpk = log( ppl(i,k) / 100. )        ! cb
!         delta(i,k)=cons1+(xkappa(i,k)-xkappa(i,k+1))*alnpk
!    &                    +xkappa(i,k)*dpp(i,k)*rpp(i,k)*cons2
!       enddo
!     enddo

! Mishra correction
! xkappai
      do i=1,lons_lat
        xkappai(i,levs+1) = xkappa(i,levs)
      enddo
      do k=2,levs
        do i=1,lons_lat
          xkappai(i,k) = (xkappa(i,k-1) + xkappa(i,k)) / 2.
        enddo
      enddo
!
      alpha(:,levs)=cons0
      do k=2,levs
        do i=1,lons_lat
          alpha(i,k-1)=(ppl(i,k)/ppl(i,k-1))**xkappai(i,k)
        enddo
      enddo
      betta(:,1   )=cons0
      do k=1,levs-1
        do i=1,lons_lat
          betta(i,k+1)=(ppl(i,k)/ppl(i,k+1))**xkappai(i,k+1)
        enddo
      enddo
      gamma(:,1   )=cons0
      do k=2,levs
        do i=1,lons_lat
          gamma(i,k)=cons1-xkappa(i,k)*dpp(i,k)*rpp(i,k)*cons2
        enddo
      enddo
      delta(:,levs)=cons0
      do k=1,levs-1
        do i=1,lons_lat
          delta(i,k)=cons1+xkappa(i,k)*dpp(i,k)*rpp(i,k)*cons2
        enddo
      enddo
!
      hki(:,1)       = cons0
      hkci(:,1)      = cons0
      hki(:,levs+1)  = cons0
      hkci(:,levs+1) = cons0
      do k=2,levs
        do i=1,lons_lat
          hkrt0     = (hg(i,k-1)+hg(i,k))/(thref(k-1)+thref(k))
          hki (i,k) = ck5(k)*hkrt0**rkappa
          hkci(i,k) = hki(i,k)*rkappa/(hg(i,k-1)+hg(i,k))
        enddo
      enddo
!
      do i=1,lons_lat
        dup(i,levs)=cons0
        dum(i,1   )=cons0
      enddo
      do k=1,levs-1
        do i=1,lons_lat
          dup(i,k  )=delta(i,k)*hg(i,k)-betta(i,k+1)*hg(i,k+1)
          dum(i,k+1)=alpha(i,k)*hg(i,k)-gamma(i,k+1)*hg(i,k+1)
        enddo
      enddo
!
      k=2
        do i=1,lons_lat
          wmkm1=hkci(i,k)*rdelp2(i,k-1)
          wmkp1=hkci(i,k)*rdelp2(i,  k)
          wmm(i,k-1)=wmkm1*dup(i,k-1)+wmkp1*dum(i,k)-cons1
          wmu(i,k-1)=wmkp1*dup(i,k)
        enddo
      do k=3,levs-1
        do i=1,lons_lat
          wmkm1=hkci(i,k)*rdelp2(i,k-1)
          wmkp1=hkci(i,k)*rdelp2(i,  k)
          wml(i,k-2)=wmkm1*dum(i,k-1)
          wmm(i,k-1)=wmkm1*dup(i,k-1)+wmkp1*dum(i,k)-cons1
          wmu(i,k-1)=wmkp1*dup(i,k)
        enddo
      enddo
      k=levs
        do i=1,lons_lat
          wmkm1=hkci(i,k)*rdelp2(i,k-1)
          wmkp1=hkci(i,k)*rdelp2(i,  k)
          wml(i,k-2)=wmkm1*dum(i,k-1)
          wmm(i,k-1)=wmkm1*dup(i,k-1)+wmkp1*dum(i,k)-cons1
        enddo
!
      levsb=levs
      do k=levs,2,-1
        if( ak5(k).eq.cons0 .and. bk5(k).eq.cons0 
     &                      .and. ck5(k).ne.cons0 ) levsb=k-1
      enddo
!
      wf(:,levsb:levs+1)=cons0
      do k=2,levsb
        do i=1,lons_lat
          wf(i,k-1)=bk5(k)*dpsdt(i)+fb(i,k)
     &              +hkci(i,k)*(dhdt(i,k-1)+dhdt(i,k))
        enddo
      enddo

      call tridim_gchdp_noadv_noq (lons_lat,lons_lat,levsb-1,levs+1,1,
     &                wml,wmm,wmu,wf,work,wf1)

      wflx(:,1     )=cons0
      wflx(:,levsb+1:levs+1)=cons0
      do k=2,levsb
        do i=1,lons_lat
          wflx(i,k)=wf1(i,k-1)
        enddo
      enddo

! -----------------------------------------------------------------------
      do k=1,levs
      do i=1,lons_lat
        ddpdt(i,k)= - fg(i,k) - wflx(i,k) + wflx(i,k+1)
        pdot (i,k)= wflx(i,k)
      enddo
      enddo
      pdot (1:lons_lat,     1) = 0.0
      pdot (1:lons_lat,levs+1) = 0.0

!
!     print *,' end of gfidi_hyb_gchdp. '
!!

      return
      end

c-----------------------------------------------------------------------
      subroutine tridim_gchdp_noadv_noq (l,lx,n,nx,m,cl,cm,cu,r,au,a)
c                .      .    .                                       .
c subprogram:    tridim_gchdp_noadv_noq  solves tridiagonal matrix problems.
c   prgmmr: iredell          org: w/nmc23    date: 91-05-07
c
c abstract: this routine solves multiple tridiagonal matrix problems
c   with multiple right-hand-side and solution vectors for every matrix.
c   the solutions are found by eliminating off-diagonal coefficients,
c   marching first foreward then backward along the matrix diagonal.
c   the computations are vectorized around the number of matrices.
c   no checks are made for zeroes on the diagonal or singularity.
c
c program history log:
c   97-07-30  iredell
c
c usage:    call tridim_gchdp_noadv_noq(l,lx,n,nx,m,cl,cm,cu,r,au,a)
c
c   input argument list:
c     l        - integer number of tridiagonal matrices
c     lx       - integer first dimension (lx>=l)
c     n        - integer order of the matrices
c     nx       - integer second dimension (nx>=n)
c     m        - integer number of vectors for every matrix
c     cl       - real (lx,2:n) lower diagonal matrix elements
c     cm       - real (lx,n) main diagonal matrix elements
c     cu       - real (lx,n-1) upper diagonal matrix elements
c                (may be equivalent to au if no longer needed)
c     r        - real (lx,nx,m) right-hand-side vector elements
c                (may be equivalent to a if no longer needed)
c
c   output argument list:
c     au       - real (lx,n-1) work array
c     a        - real (lx,nx,m) solution vector elements
c
c attributes:
c   language: fortran 77.
c   machine:  cray.
c
c	
c------------ implicit none vvv	       
	implicit none
		real :: fk
		integer :: i
		integer :: j
		integer :: k
		integer :: l
		integer :: lx
		integer :: m
		integer :: n
		integer :: nx
c    intro of implicit none > >
	 real cl(lx,2:n),cm(lx,n),cu(lx,n-1),r(lx,nx,m),
     &                         au(lx,n-1),a(lx,nx,m)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  march up
      do i=1,l
        fk=1./cm(i,1)
        au(i,1)=fk*cu(i,1)
      enddo
      do j=1,m
        do i=1,l
          fk=1./cm(i,1)
          a(i,1,j)=fk*r(i,1,j)
        enddo
      enddo
      do k=2,n-1
        do i=1,l
          fk=1./(cm(i,k)-cl(i,k)*au(i,k-1))
          au(i,k)=fk*cu(i,k)
        enddo
        do j=1,m
          do i=1,l
            fk=1./(cm(i,k)-cl(i,k)*au(i,k-1))
            a(i,k,j)=fk*(r(i,k,j)-cl(i,k)*a(i,k-1,j))
          enddo
        enddo
      enddo
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  march down
      do j=1,m
        do i=1,l
          fk=1./(cm(i,n)-cl(i,n)*au(i,n-1))
          a(i,n,j)=fk*(r(i,n,j)-cl(i,n)*a(i,n-1,j))
        enddo
      enddo
      do k=n-1,1,-1
        do j=1,m
          do i=1,l
            a(i,k,j)=a(i,k,j)-au(i,k)*a(i,k+1,j)
          enddo
        enddo
      enddo
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      return
      end
