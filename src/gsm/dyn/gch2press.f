      subroutine gch2press(njeff,nsize_ar,pgr,thgr,prsl,dprs)
!
! hmhj : this is modified hybrid by finite difference from henry juang
!        thgr can be t or h
! 2011 02 20 : henry jaung, add options for NDSL
!
 
      use gfs_dyn_machine , only : kind_grid
 
      use gfs_dyn_resol_def
      use namelist_dynamics_def
      use gfs_dyn_vert_def
      use gfs_dyn_coordinate_def
      use gfs_dyn_physcons, cp => con_cp , rd => con_rd, rk => con_rocp
      implicit none

      real(kind=kind_grid), parameter :: PT01=0.01, rkappa = cp / rd
      real(kind=kind_grid) prsl(nsize_ar,levs)
      real(kind=kind_grid) dprs(nsize_ar,levs)
      real(kind=kind_grid) pgr(nsize_ar)
      real(kind=kind_grid) thgr(nsize_ar,levs)

      real(kind=kind_grid) ppi(njeff,levs+1)
      real(kind=kind_grid) tki(njeff,levs+1)
      real(kind=kind_grid) tkrt0
 
      logical adjusted
      integer njeff,nsize_ar
      integer i,k
      real dpp,dpfact,pl,ph,pp
      integer kl,kh,kk,kl1,kh1
 
!     print *,' enter  gch2press '

      tki = 0.0
      do k=2,levs
        do i=1,njeff
          tkrt0 = (thgr(i,k-1)+thgr(i,k))/(thref(k-1)+thref(k))
          tki (i,k)=ck5(k)*tkrt0**rkappa
        enddo
      enddo
      do k=1,levp1
        do i=1,njeff
          ppi(i,k)  = ak5(k)+bk5(k)*pgr(i)+tki(i,k)
        enddo
      enddo
!
      do k=1,levs
        do i=1,njeff
          prsl(i,k) = (ppi(i,k)+ppi(i,k+1))*0.5
          dprs(i,k) = (ppi(i,k)-ppi(i,k+1))
        enddo
      enddo

!     print *,' leave gch2press. '

      return
      end
