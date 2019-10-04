      subroutine hyb2press(njeff,nsize_ar, pgr,prsl,dprs)
 
      use gfs_dyn_machine , only : kind_grid
 
 
      use gfs_dyn_resol_def
      use gfs_dyn_coordinate_def
      use gfs_dyn_physcons, rk => con_rocp
      implicit none

 
      real (kind=kind_grid), parameter :: rk1 = rk + 1.0, rkr = 1.0/rk
     &,                                   R100=100.0, PT01=0.01
      integer njeff,nsize_ar
      real(kind=kind_grid) pgr(nsize_ar)
      real(kind=kind_grid) prsl(nsize_ar,levs), dprs(nsize_ar,levs)
      real(kind=kind_grid) tem,pressfc(nsize_ar)
 
      real(kind=kind_grid) ppi(njeff,levs+1),ppik(njeff,levs+1)
      integer iq,ilat,me
      integer i,k
 
      do i=1,njeff
        pressfc(i) = exp(pgr(i))
      enddo
      do k=1,levp1
        do i=1,njeff
          ppi(i,levs+2-k) = ak5(k) + bk5(k)*pressfc(i) ! prsi are now pressures
        enddo
      enddo
 
      do i=1,njeff
        ppik(i,1) = (ppi(i,1)*PT01) ** rk
      enddo
      do k=1,levs
        do i=1,njeff
          ppik(i,k+1) = (ppi(i,k+1)*PT01) ** rk
          tem         = rk1 * (ppi(i,k) - ppi(i,k+1))
          ppik(i,k)   = (ppik(i,k)*ppi(i,k)-ppik(i,k+1)*ppi(i,k+1))
     &                 / tem
          prsl(i,k)    = R100 * ppik(i,k) ** rkr
          dprs(i,k)    = ppi(i,k) - ppi(i,k+1)
        enddo
      enddo
 
      return
      end
