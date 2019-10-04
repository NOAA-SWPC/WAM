module module_dissip
use findmaxmin2
use stencilprint
contains

!*********************************************************************
!     dissip
!	Lateral momentum dissipation
!       R.Bleck                 September 2011
!*********************************************************************

  subroutine dissip (u_vel, v_vel, delp, dfflen, biharm)
  use module_control  ,only: npp,nip,nvl,nvlp1,PrintIpnDiag
  use module_constants,only: nprox,prox,rarea,sideln,rprox_ln,		&
                             deg_lat,cs,sn,actual,nedge,permedge
  implicit none
!  External type and dimension:

! --- 'diffusion length' dfflen = (diffusivity) * (time step) / (mesh size)
! ---                           = (diffusion velocity) x (time step)

  real   ,intent (IN)    :: dfflen(nvl)		! diffusion length scale (m)
  logical,intent (IN), optional :: biharm	! if yes, do biharmonic smoothg
!SMS$DISTRIBUTE (dh,nip) BEGIN
  real   ,intent (INOUT) :: u_vel (nvl,nip)	! field to be diffused
  real   ,intent (INOUT) :: v_vel (nvl,nip)	! field to be diffused
  real   ,intent (IN)    :: delp  (nvl,nip)	! lyr thknss, Pa

! Local variables:
  real    :: ufxrot(nvl,npp,nip)	! u momentum flux across edges
  real    :: vfxrot(nvl,npp,nip)	! v momentum flux across edges
!SMS$DISTRIBUTE END

! Local variables
  integer   :: k		! layer index
  integer   :: ipn		! icos point index
  integer   :: ipx		! neighbor across joint edge
  integer   :: edg,edgcount	! icos edge index
  real      :: factor
  real      :: uxy1,uxy2,vxy1,vxy2
  real      :: uflux,vflux,uold,vold
  character :: string*5
  logical   :: vrbos
  real,parameter   :: thshld = 1.e-11
  real      :: hfharm,a,b	! 0.5 * harmonic average
  hfharm(a,b) = a*b/(a+b)	! (see Appx.D, 1992 MICOM paper)

! print '(a/(10f8.1))','entering subr.dissip with dfflen =',dfflen

  if (maxval(dfflen).eq.0.) return

  call stencl(u_vel,nvl,1.,'(dissip) -u- input')
  call stencl(v_vel,nvl,1.,'(dissip) -v- input')

! do k=1,nvl,7
!  write (string,'(a,i3)') 'k=',k
!  call findmxmn2(u_vel,nvl,nip,k,'(dissip) u-in '//string)
!  call findmxmn2(v_vel,nvl,nip,k,'(dissip) v-in '//string)
! end do

!SMS$PARALLEL (dh,ipn) BEGIN
!SMS$EXCHANGE(delp,u_vel,v_vel)

!SMS$HALO_COMP(<1,1>) BEGIN
  do ipn = 1,nip			! horizontal loop
   vrbos=ipn.eq.PrintIpnDiag
   do edgcount = 1,nedge(ipn)		! loop through edges
    edg = permedge(edgcount,ipn)
    ipx = prox(edg,ipn)			! neighbor across shared edge
    do k = 1,nvl			! loop through layers
     if (dfflen(k).gt.0.) then

! --- Transform u,v at neighboring icos pt to local coord.system.
! --- cs and sn are coordinate transformation constants.
! --- uxy,vxy are values of u and v rotated into local system.

      uxy1= cs(1,edg,ipn)*u_vel(k,ipn)+sn(1,edg,ipn)*v_vel(k,ipn)
      vxy1=-sn(1,edg,ipn)*u_vel(k,ipn)+cs(1,edg,ipn)*v_vel(k,ipn)
      uxy2= cs(2,edg,ipn)*u_vel(k,ipx)+sn(2,edg,ipn)*v_vel(k,ipx)
      vxy2=-sn(2,edg,ipn)*u_vel(k,ipx)+cs(2,edg,ipn)*v_vel(k,ipx)

! --- momentum fluxes (pos.inward) in local (rotated) coord.system:
      ufxrot(k,edg,ipn) = (uxy2-uxy1)*dfflen(k	)			&
        *sideln(edg,ipn)*2.*hfharm(max(delp(k,ipn),thshld)		&
                                  ,max(delp(k,ipx),thshld))
      vfxrot(k,edg,ipn) = (vxy2-vxy1)*dfflen(k)				&
        *sideln(edg,ipn)*2.*hfharm(max(delp(k,ipn),thshld)		&
                                  ,max(delp(k,ipx),thshld))

      if (vrbos .and. mod(k,7).eq.1) then
!SMS$IGNORE BEGIN
       print 101,'orig u,v at',ipn,actual(ipx),k,edg,   &
         u_vel(k,ipn),v_vel(k,ipn),u_vel(k,ipx),v_vel(k,ipx)
       print 101,' rot u,v at',ipn,actual(ipx),k,edg,   &
         uxy1,vxy1,uxy2,vxy2
 101   format ('(dissip) ',a,2i7,2i3,3(f10.2,f8.2))
       factor=rarea(ipn)/max(delp(k,ipn),thshld)
       print 102,' u/vflx in rotated system',ipn,			&
         actual(ipx),k,edg,ufxrot(k,edg,ipn)*factor,			&
         vfxrot(k,edg,ipn)*factor
 102   format ('(dissip) ',a,2i7,2i3,2f7.2)
!SMS$IGNORE END
      end if

     end if				! dfflen > 0
    end do				! loop through layers
   end do				! loop through edges
  end do				! horizontal loop
!SMS$HALO_COMP END

  do ipn = 1,nip			! horizontal loop
   vrbos=ipn.eq.PrintIpnDiag
   do edg = 1,nprox(ipn)
    ipx = prox(edg,ipn)			! neighbor across shared edge
    do k = 1,nvl			! loop through layers
     if (dfflen(k).gt.0.) then
      factor=rarea(ipn)/max(delp(k,ipn),thshld)
! --- rotate momentum fluxes back to lat/lon coord.system
      uflux= cs(1,edg,ipn)*ufxrot(k,edg,ipn)-sn(1,edg,ipn)*vfxrot(k,edg,ipn)
      vflux= sn(1,edg,ipn)*ufxrot(k,edg,ipn)+cs(1,edg,ipn)*vfxrot(k,edg,ipn)
      uold=u_vel(k,ipn)
      vold=v_vel(k,ipn)
      u_vel(k,ipn)=u_vel(k,ipn)+uflux*factor
      v_vel(k,ipn)=v_vel(k,ipn)+vflux*factor

      if (vrbos .and. mod(k,7).eq.1) then
!SMS$IGNORE BEGIN
       print 102,'u/vflx in lat/lon system',ipn,			&
         actual(ipx),k,edg,uflux*factor,vflux*factor
       print 101,'old/new u,v',ipn,actual(ipx),k,edg,			&
        uold,vold,u_vel(k,ipn),v_vel(k,ipn)
!SMS$IGNORE END
      end if

     end if				! dfflen > 0
    end do				! loop through layers
   end do				! loop through edges
  end do				! horizontal loop

!SMS$PARALLEL END

  call stencl(u_vel,nvl,1.,'(dissip) -u- output')
  call stencl(v_vel,nvl,1.,'(dissip) -v- output')

! do k=1,nvl,7
!  write (string,'(a,i3)') 'k=',k
!  call findmxmn2(u_vel,nvl,nip,k,'(dissip) u-out '//string)
!  call findmxmn2(v_vel,nvl,nip,k,'(dissip) v-out '//string)
! end do

  return
  end subroutine dissip
end module module_dissip
