      subroutine omega_hyb(njeff,nsize_ar,rcl,
     &                     expq,dphi,dlam,dg,ug,vg,vvel)
 
      use gfs_dyn_machine , only : kind_grid
 
      use gfs_dyn_resol_def
      use gfs_dyn_coordinate_def
!     use gfs_dyn_layout1 , only : me
      implicit none
 
      integer njeff,nsize_ar
 
      real(kind=kind_grid) rcl
      real(kind=kind_grid), dimension(nsize_ar,levs) :: dg, ug, vg, vvel
      real(kind=kind_grid), dimension(njeff)         :: dphi, dlam, expq

!     real(kind=kind_grid) dg(nsize_ar,levs), ug(nsize_ar,levs),
!    &                     vg(nsize_ar,levs),
!    &                     dphi(nsize_ar), dlam(nsize_ar),
!    &                     vvel(nsize_ar,levs),
!    &                     expq(nsize_ar)
 
 
       real(kind=kind_grid), dimension(njeff,levs) :: dpk, cg, cb, db
     &,                            workb, workc, prs, alfa, rlnp, rdel
       real(kind=kind_grid), dimension(levp1)      :: dot, dotinv
       real(kind=kind_grid)                        :: pk5(njeff,levp1)

!      real(kind=kind_grid)
!    &      pk5(njeff,levp1), dot(levp1), dotinv(levp1),
!    &      dpk(njeff,levs),     cg(njeff,levs),
!    &       cb(njeff,levs),     db(njeff,levs),
!    &    workb(njeff,levs),  workc(njeff,levs),  prs(njeff,levs),
!    &     alfa(njeff,levs),   rlnp(njeff,levs), rdel(njeff,levs)
 
      real(kind=kind_grid), parameter ::  cons0=0.d0, cons0p5=0.5d0
     &,                                   cons1=1.d0, cons2=2.d0
     &,                                   clog2=log(cons2)
      real(kind=kind_grid) rmin,rmax
      integer i,k,n,ifirst,il,ilat
 
!     print *,' enter omegtes_fd ' 		! hmhj

      do i=1,njeff
        expq(i) = exp(expq(i))
      enddo
      do k=1,levp1
        do i=1,njeff
          pk5(i,k) = ak5(k) + bk5(k)*expq(i)
        enddo
      enddo
 
      do k=1,levs
        do i=1,njeff
          prs(i,k)  = (pk5(i,k+1) + pk5(i,k) )*cons0p5
          dpk(i,k)  =  pk5(i,k+1) - pk5(i,k)
          rdel(i,k) =  cons1/dpk(i,k)
        enddo
      enddo
 
      k=1
      do i=1,njeff
        alfa(i,1) = clog2
        rlnp(i,1) = 99999.99
      enddo
 
      do k=2,levs
        do i=1,njeff
          rlnp(i,k) = log( pk5(i,k+1)/pk5(i,k) )
          alfa(i,k) = cons1 - ( pk5(i,k)/dpk(i,k) )*rlnp(i,k)
        enddo
      enddo
 
      do k=1,levs
        do i=1,njeff
          cg(i,k) = (ug(i,levs+1-k)*dlam(i)+vg(i,levs+1-k)*dphi(i))*rcl
        enddo
      enddo
 
      k=1
      do i=1,njeff
       db(i,1) = dg(i,levs)*dpk(i,1)
       cb(i,1) = cg(i,1)*dbk(1)
      enddo
 
      do k=1,levm1
        do i=1,njeff
          db(i,k+1) = db(i,k) + dg(i,levs-k)*dpk(i,k+1)
          cb(i,k+1) = cb(i,k) + cg(i,k+1)*dbk(k+1)
        enddo
      enddo
 
 
998   format('ilat=',i3)
999   format('k vv(k)=',i3,e13.3,' il=',i3,' slat=',f5.2,'  p=',f8.3)
 
 
      k=1
      do i=1,njeff
        workb(i,1) = alfa(i,1)*
     &            ( dg(i,levs)*dpk(i,1)+expq(i)*cb(i,1)*dbk(1) )
      enddo
 
      do k=2,levs
        do i=1,njeff
          workb(i,k) = rlnp(i,k)*( db(i,k-1)+expq(i)*cb(i,k-1) )
     &  +alfa(i,k)*( dg(i,levs+1-k)*dpk(i,k)+expq(i)*cg(i,k)*dbk(k) )
        enddo
      enddo
 
      k=1
      do i=1,njeff
        workc(i,1) = expq(i)*cg(i,1)*dbk(1)
      enddo
 
      do k=2,levs
        do i=1,njeff
          workc(i,k) = expq(i)*cg(i,k)
     &               * ( dbk(k)+ck(k)*rlnp(i,k)*rdel(i,k) )
        enddo
      enddo
 
      do k=1,levs
        do i=1,njeff
          vvel(i,levs+1-k) = rdel(i,k) * ( -workb(i,k) + workc(i,k))
     &                                 * prs(i,k)
        enddo
      enddo
 
!     print *,' leave omegtes_fd ' 		! hmhj

      return
      end
