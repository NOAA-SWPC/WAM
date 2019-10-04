module  module_transp3d
!*********************************************************************
!	3-dimensional tracer transport module designed for intermittent
!	(i.e., long time step) execution. there are 3 entries:
!
!	transp0 - initializes mass flux arrays and saves initial -dp-
!	transp1 - builds up time integral of horizontal mass fluxes
!	transp2 - performs the actual transport operation

!	R. Bleck	     March 2009
!*********************************************************************
  use findmaxmin1
  use findmaxmin2
  use stencilprint

contains
  subroutine transp0(its,cumufx,dp3d,dpinit)
  use module_control ,only: nvl,npp,nip

  implicit none
  integer,intent(IN)  :: its		! model time step
!SMS$DISTRIBUTE (dh,nip) BEGIN
  real,intent(IN)     :: dp3d  (nvl,    nip)	! layer thickness
  real,intent(OUT)    :: cumufx(nvl,npp,nip)	! time-integrated mass flux
  real,intent(OUT)    :: dpinit(nvl,    nip)	! dp at start of time integr.
!SMS$DISTRIBUTE END


!SMS$PARALLEL (dh,ico) BEGIN
  cumufx(:,:,:)=0.
  dpinit(:,:)=dp3d(:,:)
!SMS$PARALLEL END

  print *,'tracer transport arrays initialized, time step',its
  return
  end subroutine transp0


  subroutine transp1(its,nf,of,vof,adbash1,adbash2,adbash3,	&
                     cumufx,massfx)
  use module_control  ,only: nvl,npp,nip
  use module_constants,only: nprox,rarea

  implicit none
  integer,intent(IN)  :: its		! model time step
  integer,intent(IN)  :: nf,of,vof	! time slots: new,old,very old
  real   ,intent(IN)  :: adbash1,adbash2,adbash3

!SMS$DISTRIBUTE (dh,nip) BEGIN
  real,intent(IN)     :: massfx(nvl,npp,nip,3)
  real,intent(INOUT)  :: cumufx(nvl,npp,nip)
!SMS$DISTRIBUTE END
  integer    :: ico		! Index for icos point number
  integer    :: edg,k
  logical    :: vrbos

!SMS$PARALLEL (dh,ico) BEGIN
  do ico=1,nip			! horizontal loop
   do edg=1,nprox(ico)
    do k=1,nvl
     cumufx(k,edg,ico)=cumufx(k,edg,ico)			&
       +adbash1*massfx(k,edg,ico, nf)				&
       +adbash2*massfx(k,edg,ico, of)				&
       +adbash3*massfx(k,edg,ico,vof)
    end do
   end do
  end do
!SMS$PARALLEL END

! print *,'mass fluxes added to time integral,  time step',its
  return
  end subroutine transp1


  subroutine transp2 (its,						&
      tracr,cumufx,dpinit,dpfinl,					&
      ttransp3d,ttransp3dEx,ttransp3dBa,TimingBarriers )

  use module_control  ,only: nvl,nvlp1,npp,nip,dt,ntra,ntrb,PrintIpnDiag
  use module_constants,only: nprox,prox,area,rarea
  use module_fct3d

  implicit none
! External variables:
  integer,intent (IN)    :: its			! model time step
! integer,intent (IN)    :: ntr			! number of tracer fields
  real*8 ,intent (INOUT) :: ttransp3d		! computation timer
  real*8 ,intent (INOUT) :: ttransp3dEx		! halo communication timer
  real*8 ,intent (INOUT) :: ttransp3dBa		! barrier timer for task skew
  logical,intent (IN)    :: TimingBarriers	! measure task skew when .true.
!SMS$DISTRIBUTE (dh,nip) BEGIN
  real   ,intent (INOUT) :: tracr (nvl,nip,ntra+ntrb) ! tracer
  real   ,intent (IN)    :: cumufx(nvl,npp,nip)       ! time-integr. mass flx
  real   ,intent (IN)    :: dpinit(nvl,nip)           ! init'l lyr thknss
  real   ,intent (IN)    :: dpfinl(nvl,nip)           ! final lyr thknss
! Local variables:
  real :: vertfx(nvl,nip),col_xpand(nip)
  real :: field(nvl,nip)
!SMS$DISTRIBUTE END
  integer    :: k		! layer index
  integer    :: ico		! Index for icos point number
  integer    :: edg		! icos edge number index
  integer    :: type		! tracer index
  logical    :: vrbos
  character  :: string*32
  real       :: hordiv(nvl),vertdv(nvl)
  real*8  :: t1,tstart,tstop,valmin,valmax

!sms$compare_var(tracr  , "transp3d.F90 - tracr1   ")
!sms$compare_var(cumufx , "transp3d.F90 - cumufx1  ")

  call StartTimer(t1)

! --- compute the various terms in the continuity equation integrated
! --- over time interval since last call to -transp0-
! --- the continuity eqn is split into horiz. and vert. terms as follows:
! ---        (dpfinl-dpinit) + hordiv + verdiv = 0

!SMS$EXCHANGE(cumufx,dpinit,dpfinl)

!SMS$PARALLEL (dh,ico) BEGIN
!SMS$HALO_COMP(<1,1>) BEGIN
  do ico=1,nip				! horizontal loop
   vrbos=ico.eq.PrintIpnDiag
   hordiv   (:    )=0.
   col_xpand(  ico)=0.
   vertfx   (:,ico)=0.

   do edg=1,nprox(ico)			! loop through edges
    do k=1,nvl				! loop through layers
     hordiv(k)=hordiv(k)+cumufx(k,edg,ico)
    end do
   end do

   do k=1,nvl				! loop through layers
    col_xpand(ico)=col_xpand(ico)+hordiv(k)
    vertdv(k)=(dpinit(k,ico)-dpfinl(k,ico))-hordiv(k)*rarea(ico)
    if (k.eq.1) then
     vertfx(k,ico)=                vertdv(k)
    else
     vertfx(k,ico)=vertfx(k-1,ico)+vertdv(k)
    end if

    if (vrbos) then
!SMS$ignore begin
    write (*,'(i7,i3,a,3es12.4)') ico,k,			&
    ' transp2 hordiv,col_xpand,vertfx:',hordiv(k)*rarea(ico),	&
      col_xpand(ico)*rarea(ico),vertfx(k,ico)
    call flush(6)
!SMS$ignore end
    end if

   end do
  end do				! horizontal loop
!SMS$HALO_COMP END
!SMS$PARALLEL END

! call findmxmn1(col_xpand,nip,'col_xpand')
! do k=1,nvl,7
!  write (string,'(a,i3,a)') 'lyr',k,' vertfx'
!  call findmxmn2(vertfx,nvl,nip,k,string)
! end do
! print *
! call flush(6)

! --- having determined the vertical flux term in the time-integrated
! --- continuity eqn, we can now perform the actual tracer transport

  if (ntrb.gt.0) then

   do type=1,ntrb	! loop over class B tracers
!SMS$PARALLEL (dh,ico) BEGIN
    do ico=1,nip
     do k=1,nvl
      field(k,ico)=tracr(k,ico,ntra+type)
     end do
    end do
!SMS$PARALLEL END

    write (string,'(a,i6,a,i2,a)') '(transp2) step',its,'  trcr',type,' in'
    call stencl(field,nvl,1.,trim(string))
   end do		! loop over tracers

   call StartTimer(tstart)

   call fct3d(tracr,ntra+ntrb,ntra+1,ntra+ntrb,cumufx,vertfx,		&
     area,rarea,dpinit,dpfinl,.false.)

   tstop=0
   call IncrementTimer(tstart,tstop)
   valmin=tstop*1.e3
   valmax=valmin
!SMS$REDUCE(valmin,MIN)
!SMS$REDUCE(valmax,MAX)
   print '(a,i2,a,2(i7,a))','time spent in subr fct3d (',ntrb,		&
    ' tracers)',nint(valmin),' -',nint(valmax),' msec'

   do type=1,ntrb	! loop over class B tracers
!SMS$PARALLEL (dh,ico) BEGIN
    do ico=1,nip
     do k=1,nvl
      field(k,ico)=tracr(k,ico,ntra+type)
     end do
    end do
!SMS$PARALLEL END

    write (string,'(a,i6,a,i2,a)') '(transp2) step',its,'  trcr',type,' out'
    call stencl(field,nvl,1.,trim(string))
   end do		! loop over tracers

  end if		!  ntrb > 0

!sms$compare_var(tracr  , "transp3d.F90 - tracr2   ")
!sms$compare_var(cumufx , "transp3d.F90 - cumufx2  ")

  call IncrementTimer(t1,ttransp3d)

  return
  end subroutine transp2
end module module_transp3d
