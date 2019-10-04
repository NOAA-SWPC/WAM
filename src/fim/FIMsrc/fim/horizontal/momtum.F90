module module_momtum
use stencilprint
use findmaxmin2
use findmaxmin3
contains
!*********************************************************************
!     momtum
!	Solves momentum equations
!	Alexander E. MacDonald			12/22/2004
!	J. Lee					September  2005
!	R. Bleck   major rewrite		April      2008
!	R. Bleck   removed omega diagnostics	August     2009
!*********************************************************************

  subroutine momtum (its,nf,of,vof,adbash1,adbash2,adbash3,	&
    u_velo,v_velo,exner,relvort,				&
    u_edg,v_edg,trcr_edg,bnll_edg,u_tndcy,v_tndcy,dp3d)

  use module_control  ,only: nvl,nvlp1,npp,nip,nabl,dt,ntra,	&
                             rleigh_light,rleigh_heavy,		&
                             PrintIpnDiag,veldff_bkgnd,veldff_boost
  use module_constants,only: nprox,rarea,sidevec_c,sidevec_e,corio
  use module_dissip

  implicit none

!  External type and dimension:
  integer,intent (IN)  :: its			! model time step
  integer,intent (IN)  :: nf,of,vof		! time slots: new,old,very old
  real   ,intent (IN)  :: adbash1,adbash2,adbash3
!SMS$DISTRIBUTE (dh,nip) BEGIN
  real   ,intent (IN)    :: u_edg    (nvl,npp,nip)
  real   ,intent (IN)    :: v_edg    (nvl,npp,nip)
  real   ,intent (INOUT) :: u_velo   (nvl,nip)
  real   ,intent (INOUT) :: v_velo   (nvl,nip)
  real   ,intent (IN)    :: exner    (nvlp1,nip)
  real   ,intent (IN)    :: trcr_edg (nvl,npp,nip,ntra)
  real   ,intent (IN)    :: bnll_edg (nvl,npp,nip)
  real   ,intent (INOUT) :: u_tndcy  (nvl,nip,nabl)
  real   ,intent (INOUT) :: v_tndcy  (nvl,nip,nabl)
  real   ,intent (OUT)   :: relvort  (nvl,nip)
  real   ,intent (IN)    :: dp3d     (nvl,nip)     ! layer thickness

! Local variables
  real pgfx(nvl,nip),pgfy(nvl,nip)
!SMS$DISTRIBUTE END

  integer   :: k		! layer index
  integer   :: ipn		! icos point index
  integer   :: edg		! icos edge index
  integer   :: ndamp		! number of layer subjected to dissipation
  character :: string*24
  logical   :: vrbos
  real      :: factor,wgt
  real      :: uzeta,vzeta,uold,vold
  real      :: dfflen(nvl)

!!  do k=1,nvl,7
!!   write (string,'(a,i3,a)') 'k=',k,' u_velo'
!!   call findmxmn2(u_velo,nvl,nip,k,string)
!!
!!   write (string,'(a,i3,a)') 'k=',k,' v_velo'
!!   call findmxmn2(v_velo,nvl,nip,k,string)
!!  end do
!!  print *

! --- dampen gravity waves near model top by mixing momentum laterally

! --- 'diffusion length' dfflen = (diffusivity) * (time step) / (mesh size)
! ---                           = (diffusion velocity) x (time step)

   ndamp=0.1*nvl		! number of layers subjected to dissipation
   do k=1,nvl
    wgt=max(0.,float(k-nvl+ndamp)/float(ndamp))
    dfflen(k)=dt*(veldff_bkgnd*(1.-wgt)+veldff_boost*wgt)
   end do

   call dissip(u_velo,v_velo,dp3d,dfflen,.false.)

!SMS$PARALLEL (dh,ipn) BEGIN
!  Initialize line integrals:
  pgfx   (:,:)=0.
  pgfy   (:,:)=0.
  relvort(:,:)=0.

!sms$compare_var(u_tndcy , "momtum.F90 - u_tndcy5 ")
!sms$compare_var(v_tndcy , "momtum.F90 - v_tndcy5 ")
!sms$compare_var(exner   , "momtum.F90 - exner5 ")
!sms$compare_var(bnll_edg, "momtum.F90 - bnll_edg5 ")

  do ipn=1,nip	                ! horizontal loop
   vrbos=ipn.eq.PrintIpnDiag

   ! loop through edges and compute line integrals of bernoulli function,
   ! potential temperature, and pressure gradient

   do edg=1,nprox(ipn)		! loop through edges
    do k=1,nvl			! loop through layers

     pgfx(k,ipn) = pgfx(k,ipn)					&
              -bnll_edg(k,edg,ipn)*sidevec_c(2,edg,ipn)		&
              +.5*(exner(k,ipn)+exner(k+1,ipn))			&
              *trcr_edg(k,edg,ipn,1)*sidevec_c(2,edg,ipn)
     pgfy(k,ipn) = pgfy(k,ipn)					&
              +bnll_edg(k,edg,ipn)*sidevec_c(1,edg,ipn)		&
              -.5*(exner(k,ipn)+exner(k+1,ipn))			&
              *trcr_edg(k,edg,ipn,1)*sidevec_c(1,edg,ipn)

   !  Vorticity is calculated as a line integral of the tangential wind
   !  component given by the dot product of the wind vector with sidevec.
   !  Sidevec is a vectorial representation of the edge.
 
     relvort(k,ipn) = relvort(k,ipn)				&
            + ((sidevec_e(1,edg,ipn)*u_edg(k,edg,ipn)		&
            +   sidevec_e(2,edg,ipn)*v_edg(k,edg,ipn))) 

    enddo	                ! loop through layers
   enddo                        ! loop through edges

!SMS$ignore begin
  if (vrbos) write (*,100) its,ipn
 100 format ('m o m t u m   time step',i6,'  ipn =',i8/		&
      4x,'   uold   unew  utdcy  gradp  corio',			&
      4x,'   vold   vnew  vtdcy  gradp  corio')
!SMS$ignore end

   do k=1,nvl			! loop through layers
    pgfx   (k,ipn) = pgfx   (k,ipn)*rarea(ipn)
    pgfy   (k,ipn) = pgfy   (k,ipn)*rarea(ipn)
    relvort(k,ipn) = relvort(k,ipn)*rarea(ipn)

    uzeta=(corio(ipn)+relvort(k,ipn))*u_velo(k,ipn)
    vzeta=(corio(ipn)+relvort(k,ipn))*v_velo(k,ipn)

    !  u/v tendcy is the sum of bernoulli fct. gradient and coriolis term
    u_tndcy (k,ipn,nf) = pgfx(k,ipn) + vzeta
    v_tndcy (k,ipn,nf) = pgfy(k,ipn) - uzeta

    ! advance velocity field in time
    uold=u_velo(k,ipn)
    vold=v_velo(k,ipn)

    u_velo(k,ipn) = u_velo(k,ipn)			&
            +adbash1*u_tndcy(k,ipn, nf)			&
            +adbash2*u_tndcy(k,ipn, of)			&
            +adbash3*u_tndcy(k,ipn,vof)
    v_velo(k,ipn) = v_velo(k,ipn)			&
            +adbash1*v_tndcy(k,ipn, nf)			&
            +adbash2*v_tndcy(k,ipn, of)			&
            +adbash3*v_tndcy(k,ipn,vof)

!SMS$ignore begin
   if (vrbos .and. mod(k,7).eq.1) write (*,101) k,			&
    uold,u_velo(k,ipn),u_tndcy(k,ipn,nf)*dt,pgfx(k,ipn)*dt, vzeta*dt,	&
    vold,v_velo(k,ipn),v_tndcy(k,ipn,nf)*dt,pgfy(k,ipn)*dt,-uzeta*dt
 101 format (i4,2f7.1,3f7.2,4x,2f7.1,3f7.2)
!SMS$ignore end
   enddo	                ! loop through layers

   ! Rayleigh damping of u,v near model top

   ndamp=0.25*nvl

!SMS$ignore begin
  if (vrbos) write (*,107) ipn,'u,v bfore Rayleigh damping:',	&
     (k,u_velo(k,ipn),v_velo(k,ipn),k=nvl-ndamp,nvl)
 107 format ('ipn=',i8,3x,a/(i15,2f9.2))
!SMS$ignore end

    if (u_velo(nvl,ipn)**2+v_velo(nvl,ipn)**2 .gt. 1.e4) then 
      do k=nvl-ndamp,nvl
        factor=1.-(dt/86400.)*rleigh_heavy*2.**(30.*(k-nvl)/nvl)
        u_velo(k,ipn)=u_velo(k,ipn)*factor
        v_velo(k,ipn)=v_velo(k,ipn)*factor
      end do
    else 
      do k=nvl-ndamp,nvl
        factor=1.-(dt/86400.)*rleigh_light*2.**(30.*(k-nvl)/nvl)
        u_velo(k,ipn)=u_velo(k,ipn)*factor
        v_velo(k,ipn)=v_velo(k,ipn)*factor
      end do
    end if

!SMS$ignore begin
  if (vrbos) write (*,107) ipn,'u,v after Rayleigh damping:',	&
     (k,u_velo(k,ipn),v_velo(k,ipn),k=nvl-ndamp,nvl)
!SMS$ignore end

  end do		                ! horizontal loop
!SMS$PARALLEL END

  write (string,'(a,i6,2x)') 'step',its
  call stencl(u_velo,nvl,1.,string(1:12)//'(atm momtum) new u')
  call stencl(v_velo,nvl,1.,string(1:12)//'(atm momtum) new v')
! call stencl(pgfx,nvl,1.e3,string(1:12)//'(atm momtum) pgfx x 1000')
! call stencl(pgfy,nvl,1.e3,string(1:12)//'(atm momtum) pgfy x 1000')

!sms$compare_var(u_tndcy, "momtum.F90 - u_tndcy6 ")
!sms$compare_var(v_tndcy, "momtum.F90 - v_tndcy6 ")
!sms$compare_var(u_velo , "momtum.F90 - u_velo6 ")
!sms$compare_var(v_velo , "momtum.F90 - v_velo6 ")

!!  do k=1,nvl,7
!!   write (string,'(a,i3,a)') 'k=',k,' u_tndcy'
!!   call findmxmn3(u_tndcy,nvl,nip,3,k,nf,string)
!!
!!   write (string,'(a,i3,a)') 'k=',k,' v_tndcy'
!!   call findmxmn3(v_tndcy,nvl,nip,3,k,nf,string)
!!  end do
!!  print *

  return
  end subroutine momtum
  end module module_momtum
