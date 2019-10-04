      subroutine omega_gcdp(njeff,nsize_ar,rcl,
     &                      dpg,dpphi,dplam,dg,ug,vg,vvel) 
 
!
! hmhj : this is modified hybrid by finite difference from henry juang
! 20110220    henry jaung modified code to fit mass_dp and ndslfv
!
      use gfs_dyn_machine , only : kind_grid
 
      use gfs_dyn_resol_def
      use gfs_dyn_coordinate_def
      use gfs_dyn_physcons, rd => con_rd, cp => con_cp
      implicit none
 
      integer njeff,nsize_ar
      integer i,k
 
      real(kind=kind_grid) rcl
      real(kind=kind_grid) dg   (nsize_ar,levs),    ug(nsize_ar,levs),
     &                     vg   (nsize_ar,levs),  vvel(nsize_ar,levs),
     &                     dpphi(nsize_ar,levs), dplam(nsize_ar,levs),
     &                     dpg  (nsize_ar,levs)
 
       real(kind=kind_grid)
     &      pp (njeff,levp1),  dpx(njeff,levp1), dpy(njeff,levp1),
     &      dpp(njeff,levs),   dppx(njeff,levs),dppy(njeff,levs),
     &       db(njeff,levp1),  appx(njeff,levs),appy(njeff,levs) 
 
!     print *,' enter  omegtes_gc_fd '				! hmhj
!
      pp (:,levs+1) = 0.0
      dpx(:,levs+1) = 0.0
      dpy(:,levs+1) = 0.0
      do k=levs,1,-1
        do i=1,njeff
          pp (i,k)=pp (i,k+1)+dpg  (i,k)
          dpx(i,k)=dpx(i,k+1)+dplam(i,k)
          dpy(i,k)=dpy(i,k+1)+dpphi(i,k)
        enddo
      enddo
      do k=1,levs
        do i=1,njeff
          dpx(i,k)=dpx(i,k)*rcl
          dpy(i,k)=dpy(i,k)*rcl
        enddo
      enddo
      do k=1,levs
        do i=1,njeff
          dpp (i,k) = pp (i,k) - pp (i,k+1)
          dppx(i,k) = dpx(i,k) - dpx(i,k+1)
          dppy(i,k) = dpy(i,k) - dpy(i,k+1)
          appx(i,k) = dpx(i,k) + dpx(i,k+1)
          appy(i,k) = dpy(i,k) + dpy(i,k+1)
        enddo
      enddo
      do i=1,njeff
        db(i,levs+1) = 0.e0
      enddo
      do k=levs,1,-1
        do i=1,njeff
          db(i,k)=db(i,k+1)+dpp(i,k)*dg(i,k)
     &                     +ug(i,k)*dppx(i,k)+vg(i,k)*dppy(i,k)
        enddo
      enddo
      do k=1,levs
        do i=1,njeff
          vvel(i,k)= ug(i,k)*appx(i,k)+vg(i,k)*appy(i,k)
     &              -db(i,k)-db(i,k+1)
          vvel(i,k)= 0.5 * vvel(i,k)
        enddo
      enddo

!     print *,' leave omegtes_gc_h. '				! hmhj
 
      return
      end
