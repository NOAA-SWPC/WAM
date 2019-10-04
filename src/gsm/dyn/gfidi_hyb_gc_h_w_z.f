      subroutine gfidi_hyb_gc_h_w_z(lon_dim,lons_lat,lat,
     &  dg,hg,ug,vg,rqg,dphi,dlam,ps,zsphi,zslam,zs,
     &  rcl,
     &  dhdf,dhdl,drqdf,drqdl,wg,zg,me)
 
!
! version of ppi=ak+bk*psfc+ck*(h/h0)^(1/kappa)
!
! hmhj : this is modified hybrid by finite difference from hann-ming henry juang 
! hmhj : use h=CpT instead of Tv for thermodynamical equation
! hmhj : for hydrostatic w initial value
!
 
      use gfs_dyn_machine , only : kind_grid
 
      use gfs_dyn_resol_def
      use namelist_dynamics_def
      use gfs_dyn_coordinate_def
      use gfs_dyn_tracer_const
      use gfs_dyn_physcons, rerth => con_rerth
     &,             rd => con_rd, cpd => con_cp
     &,             omega => con_omega, kappa => con_rocp
     &,             grav  => con_g
      implicit none

      real (kind=kind_grid), parameter :: rkappa=1.0/kappa
      integer lon_dim,lons_lat,me
      integer i,k,kk,n,lat
      real rcl
      real det,hkrt0,wmkm1,wmkp1
      real
     1    dg(lon_dim,levs), hg(lon_dim,levs),  
     2    ug(lon_dim,levs), vg(lon_dim,levs),
     2    wg(lon_dim,levs), zg(lon_dim,levs),
     2   rqg(lon_dim,levs,ntrac),
     3  dphi(lon_dim), dlam(lon_dim), ps(lon_dim),
     3 zsphi(lon_dim),zslam(lon_dim), zs(lon_dim)
      real
     1  dhdf(lon_dim,levs),       dhdl(lon_dim,levs),
     1  dudf(lon_dim,levs),       dudl(lon_dim,levs),
     1  dvdf(lon_dim,levs),       dvdl(lon_dim,levs),
     1  drqdf(lon_dim,levs,ntrac), drqdl(lon_dim,levs,ntrac)
      real
     1  dudt(lon_dim,levs),       dvdt(lon_dim,levs),
     1  dpsdt(lon_dim),
     1  dhdt(lon_dim,levs),
     1  drqdt(lon_dim,levs,ntrac)
      real drdf (lon_dim,levs), drdl (lon_dim,levs), drdt (lon_dim,levs)
      real dcpdf(lon_dim,levs), dcpdl(lon_dim,levs), dcpdt(lon_dim,levs)
      real dkhdf(lon_dim,levs), dkhdl(lon_dim,levs), dkhdt(lon_dim,levs)

! -----------------------------
c
      real fs(lons_lat),rm2(lons_lat),xm2(lons_lat)
      real rdelp2(lons_lat,levs)
      real xr(lons_lat,levs),xcp(lons_lat,levs)
      real xkappa(lons_lat,levs),xkappai(lons_lat,levs+1)
      real cg(lons_lat,levs),ek(lons_lat,levs)
      real fb(lons_lat,levs+1),fg(lons_lat,levs)
      real dpxi(lons_lat,levs+1),dpyi(lons_lat,levs+1)
      real dpti(lons_lat,levs+1)
      real ppi(lons_lat,levs+1), ppl(lons_lat,levs)
      real hki (lons_lat,levs+1),hkci(lons_lat,levs+1)
      real dpp(lons_lat,levs),rpp(lons_lat,levs)
      real dlnpx(lons_lat,levs),dlnpy(lons_lat,levs)
      real dphix (lons_lat,levs),dphiy (lons_lat,levs)
      real dphixk(lons_lat,levs),dphiyk(lons_lat,levs)
      real dphit (lons_lat,levs),dphitk(lons_lat,levs)
      real  phii (lons_lat,levs), phiik(lons_lat,levs)
      real wf(lons_lat,levs+1),wf1(lons_lat,levs+1)
      real wflx(lons_lat,levs+1)
      real wml(lons_lat,levs),wmm(lons_lat,levs),wmu(lons_lat,levs)
      real work(lons_lat,levs)
      real dup(lons_lat,levs),dum(lons_lat,levs)
      real alpha(lons_lat,levs),betta(lons_lat,levs)
      real gamma(lons_lat,levs),delta(lons_lat,levs), alnpk
      real zadv(lons_lat,levs,1+ntrac)
      real sumdf(lons_lat,levs),sumdl(lons_lat,levs)
      real sumrq(lons_lat,levs),sumdt(lons_lat,levs)
!
      real, parameter ::  cons0=0.0, cons1=1.0, cons2=2.0, cons0p5=0.5
      integer levsb
!
! ----- prepare xr, xcp, xkapa and their derivatives
!
      xr    = cons0
      drdl  = cons0
      drdf  = cons0
      xcp   = cons0
      dcpdl = cons0
      dcpdf = cons0
      sumdl = cons0
      sumdf = cons0
      sumrq = cons0
!
      do n=1,ntrac
        if( ri(n) .ne. cons0 .and. cpi(n) .ne. cons0 ) then
          do k=1,levs
            do i=1,lons_lat
              xr   (i,k) = xr   (i,k) + rqg  (i,k,n)*ri(n)
              drdl (i,k) = drdl (i,k) + drqdl(i,k,n)*ri(n)
              drdf (i,k) = drdf (i,k) + drqdf(i,k,n)*ri(n)
!
              xcp  (i,k) = xcp  (i,k) + rqg  (i,k,n)*cpi(n)
              dcpdl(i,k) = dcpdl(i,k) + drqdl(i,k,n)*cpi(n)
              dcpdf(i,k) = dcpdf(i,k) + drqdf(i,k,n)*cpi(n)
!
              sumrq(i,k) = sumrq(i,k) + rqg  (i,k,n)
              sumdl(i,k) = sumdl(i,k) + drqdl(i,k,n)
              sumdf(i,k) = sumdf(i,k) + drqdf(i,k,n)
            enddo
          enddo
        endif
      enddo
      do k=1,levs
        do i=1,lons_lat
          drdl (i,k) = drdl (i,k) - ri(0) *sumdl(i,k)
          drdf (i,k) = drdf (i,k) - ri(0) *sumdf(i,k)
          xr (i,k)   = ( cons1 - sumrq(i,k) )*ri(0)  + xr (i,k)
!
          dcpdl(i,k) = dcpdl(i,k) - cpi(0)*sumdl(i,k)
          dcpdf(i,k) = dcpdf(i,k) - cpi(0)*sumdf(i,k)
          xcp(i,k)   = ( cons1 - sumrq(i,k) )*cpi(0) + xcp(i,k)
          xkappa(i,k) = xr(i,k) / xcp(i,k)
        enddo
      enddo
!
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
! ----- prepare dpxi, dpyi
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
      do k=1,levs+1
        do i=1,lons_lat
          ppi(i,k)  = ak5(k) + bk5(k)*ps(i) + hki(i,k)
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
            print *,' ----- dpp < 0 in gfidi at i k ',i,k
          endif
        enddo
      enddo
!
! ----------------------------------------------------------------------

      do k=1,levs+1
        do i=1,lons_lat
          dpxi(i,k)=bk5(k)*dlam(i)*rcl
          dpyi(i,k)=bk5(k)*dphi(i)*rcl
        enddo
      enddo
      do k=2,levs
        do i=1,lons_lat
          dpxi(i,k)=dpxi(i,k)+hkci(i,k)*(dhdl(i,k-1)+dhdl(i,k))
          dpyi(i,k)=dpyi(i,k)+hkci(i,k)*(dhdf(i,k-1)+dhdf(i,k))
        enddo
      enddo
!
      alpha(:,levs)=cons0
      betta(:,1   )=cons0
      do k=2,levs
        do i=1,lons_lat

!change to Misha's revision
          alpha(i,k-1)=(ppl(i,k)/ppl(i,k-1))**xkappai(i,k)
        enddo
      enddo
      do k=1,levs-1
        do i=1,lons_lat

!change to Misha's revision
           betta(i,k+1)=(ppl(i,k)/ppl(i,k+1))**xkappai(i,k+1)

        enddo
      enddo
!
!-----------
      do k=2,levs
        do i=1,lons_lat
          gamma(i,k)=cons1-xkappa(i,k)*dpp(i,k)*rpp(i,k)*cons2
        enddo
      enddo
      do k=1,levs - 1
        do i=1,lons_lat
          delta(i,k)=cons1+xkappa(i,k)*dpp(i,k)*rpp(i,k)*cons2
        enddo
      enddo
!
! ----- prepare cg and fb
 
      do k=1,levs
        do i=1,lons_lat
          fg(i,k)=ug(i,k)*(dpxi(i,k)-dpxi(i,k+1))
     &           +vg(i,k)*(dpyi(i,k)-dpyi(i,k+1))
     &           +dpp(i,k)*dg(i,k)
          cg(i,k)=ug(i,k)*(dpxi(i,k)+dpxi(i,k+1))
     &           +vg(i,k)*(dpyi(i,k)+dpyi(i,k+1))
        enddo
      enddo
!
      do i=1,lons_lat
        fb(i,levs+1)=cons0
      enddo
      do k=levs,1,-1
        do i=1,lons_lat
          fb(i,k)=fb(i,k+1)+fg(i,k)
        enddo
      enddo
!
! local change of surface pressure  d ps dt
!
      do i=1,lons_lat
        dpsdt(i) = - fb(i,1)
      enddo
!
! total derivative of virtual temperature
!
      do k=1,levs
       do i=1,lons_lat
         dhdt(i,k)=
     &              +xkappa(i,k)*hg(i,k)*rpp(i,k)*
     &                       (cg(i,k)-fb(i,k)-fb(i,k+1))
       enddo
      enddo
!
!
! horizontal advection for all
!
      do k=1,levs
        do i=1,lons_lat
          dhdt(i,k)=dhdt(i,k)
     &               -ug(i,k)*dhdl(i,k)-vg(i,k)*dhdf(i,k)
        enddo
      enddo

        do n=1,ntrac
        do k=1,levs
          do i=1,lons_lat
            drqdt(i,k,n)=
     &               -ug(i,k)*drqdl(i,k,n)-vg(i,k)*drqdf(i,k,n)
          enddo
        enddo
        enddo
!
! ------ hybrid to solve vertical flux ----------

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
      levsb=levs
      do k=levs,2,-1
        if( ak5(k).eq.cons0 .and. bk5(k).eq.cons0 
     &                      .and. ck5(k).ne.cons0 ) levsb=k-1
      enddo
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
!hmhj wf(:,levs:levs+1)=cons0
!hmhj do k=2,levs
      wf(:,levsb:levs+1)=cons0
      do k=2,levsb
        do i=1,lons_lat
          wf(i,k-1)=bk5(k)*dpsdt(i)+fb(i,k)
     &              +hkci(i,k)*(dhdt(i,k-1)+dhdt(i,k))
        enddo
      enddo

      call tridim_hyb_gc_h(lons_lat,lons_lat,levsb-1,levs+1,1,
     &                wml,wmm,wmu,wf,work,wf1)

      wflx(:,1     )=cons0
      wflx(:,levs+1)=cons0
!hmhj do k=2,levs
      wflx(:,levsb+1:levs+1)=cons0
      do k=2,levsb
        do i=1,lons_lat
          wflx(i,k)=wf1(i,k-1)
        enddo
      enddo

! -----------------------------------------------------------------------
! do vertical advection of hg first, since dup and dum are obtained
!
      do k=1,levs
        do i=1,lons_lat
          zadv(i,k,1)=-rdelp2(i,k)*
     &               (wflx(i,k)*dum(i,k)+wflx(i,k+1)*dup(i,k))
        enddo
      enddo
!
      do n=1,ntrac
        do k=1,levs-1
          do i=1,lons_lat
            dup(i,k  )=rqg(i,k,n)-rqg(i,k+1,n)
            dum(i,k+1)=rqg(i,k,n)-rqg(i,k+1,n)
          enddo
        enddo
        do k=1,levs
          do i=1,lons_lat
            zadv(i,k,1+n)=-rdelp2(i,k)*
     &               (wflx(i,k)*dum(i,k)+wflx(i,k+1)*dup(i,k))
          enddo
        enddo
!       print *,' vadv zadv ',n,(zadv(i,1,1+n),i=1,4)
      enddo

! add vertical filterd advection
      do k=1,levs
        do i=1,lons_lat
          dhdt(i,k)=dhdt(i,k)+zadv(i,k,1)
        enddo
      enddo
      do n=1,ntrac
        do k=1,levs
          do i=1,lons_lat
            drqdt(i,k,n)=drqdt(i,k,n)+zadv(i,k,1+n)
          enddo
        enddo
      enddo
      do k=1,levs+1
        do i=1,lons_lat
          dpti(i,k) = -fb(i,k) - wflx(i,k)
        enddo
      enddo

! compute dcpdt and drdt for dkappadt
      drdt  = cons0
      dcpdt = cons0
      sumdt = cons0
!
      do n=1,ntrac
        if( ri(n) .ne. cons0 .and. cpi(n) .ne. cons0 ) then
          do k=1,levs
            do i=1,lons_lat
              drdt (i,k) = drdt (i,k) + drqdt(i,k,n)*ri(n)
              dcpdt(i,k) = dcpdt(i,k) + drqdt(i,k,n)*cpi(n)
              sumdt(i,k) = sumdt(i,k) + drqdt(i,k,n)
            enddo
          enddo
        endif
      enddo
      do k=1,levs
        do i=1,lons_lat
          drdt (i,k) = drdt (i,k) - ri(0) *sumdt(i,k)
          dcpdt(i,k) = dcpdt(i,k) - cpi(0)*sumdt(i,k)
        enddo
      enddo
!
! hydrostatic to get geopotential height local derivatives wrt x, y, t
!
      do k=1,levs
        do i=1,lons_lat
          dkhdl(i,k)=( xr(i,k)*dhdl(i,k)+hg(i,k)*drdl(i,k)
     &                -xkappa(i,k)*hg(i,k)*dcpdl(i,k) )/xcp(i,k)
          dkhdf(i,k)=( xr(i,k)*dhdf(i,k)+hg(i,k)*drdf(i,k)
     &                -xkappa(i,k)*hg(i,k)*dcpdf(i,k) )/xcp(i,k)
          dkhdt(i,k)=( xr(i,k)*dhdt(i,k)+hg(i,k)*drdt(i,k)
     &                -xkappa(i,k)*hg(i,k)*dcpdt(i,k) )/xcp(i,k)
        enddo
      enddo
      do k=1,levs
        do i=1,lons_lat
          dphixk(i,k)= rpp(i,k)*( dpp(i,k)*dkhdl(i,k)
     &                  +xkappa(i,k)*hg(i,k)*( (dpxi(i,k)-dpxi(i,k+1))
     &                  -rpp(i,k)*dpp(i,k)*(dpxi(i,k)+dpxi(i,k+1)) ) )
          dphiyk(i,k)= rpp(i,k)*( dpp(i,k)*dkhdf(i,k)
     &                  +xkappa(i,k)*hg(i,k)*( (dpyi(i,k)-dpyi(i,k+1))
     &                  -rpp(i,k)*dpp(i,k)*(dpyi(i,k)+dpyi(i,k+1)) ) )
          dphitk(i,k)= rpp(i,k)*( dpp(i,k)*dkhdt(i,k)
     &                  +xkappa(i,k)*hg(i,k)*( (dpti(i,k)-dpti(i,k+1))
     &                  -rpp(i,k)*dpp(i,k)*(dpti(i,k)+dpti(i,k+1)) ) )
           phiik(i,k)= rpp(i,k)*xkappa(i,k)*hg(i,k)*dpp(i,k)
        enddo
      enddo
      do i=1,lons_lat
        dphix(i,1)= cons0
        dphiy(i,1)= cons0
        dphit(i,1)= cons0
         phii(i,1)= cons0
      enddo
      do k=1,levs-1
        do i=1,lons_lat
          dphix(i,k  )= dphix(i,k)+dphixk(i,k)
          dphix(i,k+1)= dphix(i,k)+dphixk(i,k)
          dphiy(i,k  )= dphiy(i,k)+dphiyk(i,k)
          dphiy(i,k+1)= dphiy(i,k)+dphiyk(i,k)
          dphit(i,k  )= dphit(i,k)+dphitk(i,k)
          dphit(i,k+1)= dphit(i,k)+dphitk(i,k)
           phii(i,k  )=  phii(i,k)+ phiik(i,k)
           phii(i,k+1)=  phii(i,k)+ phiik(i,k)
        enddo
      enddo
      do i=1,lons_lat
        dphix(i,levs)= dphix(i,levs)+dphixk(i,levs)
        dphiy(i,levs)= dphiy(i,levs)+dphiyk(i,levs)
        dphit(i,levs)= dphit(i,levs)+dphitk(i,levs)
         phii(i,levs)=  phii(i,levs)+ phiik(i,levs)
      enddo
      do k=1,levs
        do i=1,lons_lat
          dphix(i,k)= dphix(i,k) / rcl + grav*zslam(i)
          dphiy(i,k)= dphiy(i,k) / rcl + grav*zsphi(i)
           phii(i,k)=  phii(i,k)       + grav*zs   (i)
        enddo
      enddo
!
! compute local derivative of phi
      do k=1,levs
        do i=1,lons_lat
          dum(i,k) = rcl*(ug(i,k)*dphix(i,k)+vg(i,k)*dphiy(i,k))
          dup(i,k) = xkappa(i,k)*hg(i,k)/ppl(i,k)
        enddo
      enddo

      do k=1,levs
        do i=1,lons_lat
          wg(i,k)=-0.5*(wflx(i,k)+wflx(i,k+1))*dup(i,k)
     &            +dum(i,k)
     &            +dphit(i,k)
          wg(i,k) = wg(i,k)/grav
          zg(i,k) = phii(i,k)/grav
        enddo
      enddo
       
!      if( me.eq.0 ) then
!        print *,' gfidi_hyb_gc_h_w_z w: ',( wg(1,k),k=1,levs)
!        print *,' gfidi_hyb_gc_h_w_z z: ',( zg(1,k),k=1,levs)
!      endif
        
!     print *,' end of gfidi_hyb_gc_h_w_z. '
!!

      return
      end

