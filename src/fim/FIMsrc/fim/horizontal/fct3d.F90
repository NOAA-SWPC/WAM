module module_fct3d
contains

  subroutine fct3d(fields,numfld,frst,last,u,w,area,rarea,thko,thkn,diagno)

! --- 3-d transport routine adapted from HYCOM

! --- fld   - field to be transported
! --- u,w   - mass fluxes (x time step) satisfying continuity equation
! ---         (w(k) > 0 means mass flows from layer k to layer k+1)
! --- area  - grid cell size
! --- rarea - inverse of area
! --- thko,thkn - layer thickness at previous and new time step
  
  use module_control  ,only: npp,nip,nvl,PrintIpnDiag,ntra,ntrb
  use module_constants,only: nprox,prox,proxs
  use findmaxmin2
  use findmaxmin3

  implicit none
  integer,intent(IN) :: numfld		! total number of fields in 'fields'
  integer,intent(IN) :: frst,last	! range of fields to be advected
  logical,intent(IN) :: diagno		! activate diagnostic calculations
!SMS$DISTRIBUTE (dh,nip) BEGIN
  real,intent(INOUT) :: fields(nvl,nip,numfld)
  real   ,intent(IN) :: u(nvl,npp,nip),w(nvl,nip),		&
                        thko(nvl,nip),thkn(nvl,nip),		&
                        area(nip),rarea(nip)
! Local variables:
  real fld(nvl,nip),vertfx(nvl,nip),vertdv(nvl,nip),		&
       fmx(nvl,nip),fmn(nvl,nip),flp(nvl,nip),fln(nvl,nip),	&
       flx(nvl,npp,nip),uan(nvl,npp,nip),hordiv(nvl,nip)
!SMS$DISTRIBUTE END
  real a(nvl),b(nvl),c(nvl),dx,fcdx,yl,yr,totin,totou
  real q,clip,vlume(nvl),drift(nvl),bfore,after,slab,dslab,	&
       thkchg,dxlft,dxmid,dxrgt,bforek(nvl),afterk(nvl),	&
       amount,var1,var2
  integer i,edg,ix,edx,k,kp,n
  character string*24
  logical vrbos
  logical,parameter :: recovr=.false.
  real   ,parameter :: athird=1./3.,epsil=1.e-11,onemu=1.e-6

! print *,'entering fct3d...'

!SMS$PARALLEL (dh,i) BEGIN
  do i=1,nip
   vrbos=i.eq.PrintIpnDiag

   if (vrbos) then
! --- check mass conservation in test column
!SMS$IGNORE BEGIN
    write (*,'(i8,a/a)') i,					&
     '  fct3d -- time-integrated continuity eqn diagnostics:',	&
     '     thknss_tndcy  horiz.flxdiv   vert.flxdiv      residuum'
    do k=1,nvl			! vertical loop
     thkchg=thkn(k,i)-thko(k,i)
     hordiv(k,i)=0.
     do edg=1,nprox(i)
      hordiv(k,i)=hordiv(k,i)+u(k,edg,i)
     end do
     hordiv(k,i)=hordiv(k,i)*rarea(i)
     if (k.eq.1) then
      write (*,103) k,thkchg,hordiv(k,i),w(k,i),		&
       thkchg+hordiv(k,i)+w(k,i)
     else
      write (*,103) k,thkchg,hordiv(k,i),w(k,i)-w(k-1,i),	&
       thkchg+hordiv(k,i)+w(k,i)-w(k-1,i)
     end if
    end do			! vertical loop
103 format (i3,4es14.4)
!SMS$IGNORE END
   end if			! vrbos
  end do			! horiz. loop
!SMS$PARALLEL END

! --- optional: check mass conservation globally in select layers
  if (diagno) then

!SMS$PARALLEL (dh,i) BEGIN
   do i=1,nip			! horiz. loop
    do k=1,nvl,7		! vertical loop
     hordiv(k,i)=0.
     do edg=1,nprox(i)
      hordiv(k,i)=hordiv(k,i)+u(k,edg,i)
     end do
     hordiv(k,i)=hordiv(k,i)*rarea(i)
     if (k.eq.1) then
       hordiv(k,i)=hordiv(k,i)					&
         +w(k,i)         +thkn(k,i)-thko(k,i)
     else
       hordiv(k,i)=hordiv(k,i)					&
         +w(k,i)-w(k-1,i)+thkn(k,i)-thko(k,i)
     end if
    end do			! vertical loop
   end do			! horiz. loop
!SMS$PARALLEL END

   do k=1,nvl,7
    write (string,'(a,i3)') 'fct3d hordiv lyr',k
    call findmxmn2(hordiv,nvl,nip,k,string)
   end do			! vertical loop
  end if			! diagno

!SMS$EXCHANGE(fields,thko,thkn)

  do 1 n=frst,last		! loop over fields to be transported

  if (diagno) then
   do k=1,nvl,7
    write (string,'(a,i2.2,a,i3)') 'fct3d  in: fld',n,' lyr',k
    call findmxmn3(fields,nvl,nip,numfld,k,n,trim(string))
   end do
  end if

! --- get vertical flux by summing -fld- over upstream slab of thickness -w-

!SMS$PARALLEL (dh,i) BEGIN
!SMS$HALO_COMP(<1,1>) BEGIN
  do i=1,nip			! horiz. loop
   fld(:,i)=fields(:,i,n)

! --- fill massless cells with data from layer above or below
   do k=nvl-1,1,-1
    fld(k,i)=(fld(k,i)*thko(k,i)+fld(k+1,i)*onemu)	&
            /(         thko(k,i)+           onemu)
   end do
   do k=2,nvl
    fld(k,i)=(fld(k,i)*thko(k,i)+fld(k-1,i)*onemu)	&
            /(         thko(k,i)+           onemu)
   end do
  end do			! horiz. loop
!SMS$HALO_COMP END
!SMS$PARALLEL END

!SMS$PARALLEL (dh,i) BEGIN
  do i=1,nip			! horiz. loop
! --- fit 0th, 1st, or 2nd deg. polynomial to tracer in each cell
   a(1  )=fld(1  ,i)
   b(1  )=0.
   c(1  )=0.
   a(nvl)=fld(nvl,i)
   b(nvl)=0.
   c(nvl)=0.
   do k=2,nvl-1			! vertical loop
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! --- piecewise constant method:
!   a(k)=fld(k,i)
!   b(k)=0.
!   c(k)=0.
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! --- piecewise linear method:
! --- fit linear function a+bx to tracer in each cell (-.5 < x < +.5)
!   a(k)=fld(k,i)
!   b(k)=0.
!   if (fld(k,i).le.min(fld(k-1,i),fld(k+1,i)) .or.		&
!       fld(k,i).ge.max(fld(k-1,i),fld(k+1,i))) then
!     b(k)=0.
!   else if ((fld(k+1,i)-fld(k-1,i))*(fld(k-1,i)+fld(k+1,i)	&
!     -2.*fld(k,i)).gt.0.) then
!     b(k)=fld(k,i)-fld(k-1,i)
!   else
!     b(k)=fld(k+1,i)-fld(k,i)
!   end if
!   c(k)=0.
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! --- piecewise parabolic method:
! --- construct parabola  a+bx+cx^2  whose integral over [-.5,+.5] equals
! --- fld(k) and which passes though points yl,yr at [-.5,+.5] resp.
    dxlft=max(epsil,thko(k-1,i))
    dxmid=max(epsil,thko(k  ,i))
    dxrgt=max(epsil,thko(k+1,i))
    yl=(dxlft*fld(k,i)+dxmid*fld(k-1,i))/(dxlft+dxmid)
    yr=(dxrgt*fld(k,i)+dxmid*fld(k+1,i))/(dxrgt+dxmid)

    a(k)=1.5*fld(k,i)-.25*(yl+yr)
    b(k)=yr-yl
    c(k)=6.*(.5*(yl+yr)-fld(k,i))
    if (abs(yr-yl) .lt. 6.*abs(.5*(yl+yr)-fld(k,i))) then
! --- apex of parabola o !urs inside interval [-.5,+.5], implying an over-
! --- or undershoot situation. change curve to prevent over/undershoots.
     if (abs(yr-yl) .gt. 2.*abs(.5*(yl+yr)-fld(k,i))) then
! --- put apex of parabola on edge of interval [-.5,+.5]
      if ((yr-yl)*(.5*(yl+yr)-fld(k,i)) .gt. 0.) then
! --- apex at x=-.5
       a(k)=.25*(3.*fld(k,i)+yl)
       c(k)=3.*(fld(k,i)-yl)
       b(k)=c(k)
      else
! --- apex at x=+.5
       a(k)=.25*(3.*fld(k,i)+yr)
       c(k)=3.*(fld(k,i)-yr)
       b(k)=-c(k)
      end if
     else			!  -1/6 < x < +1/6
! --- moving apex won't help. replace parabola by constant.
      a(k)=fld(k,i)
      b(k)=0.
      c(k)=0.
     end if
    end if
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   end do			! vertical loop

   do k=1,nvl-1			! vertical loop
    slab=onemu
    if (w(k,i).lt.0.) then			! interface moves up in atm.
     amount=slab*fld(k+1,i)
     kp=k
 24  kp=kp+1
     if (slab.ge.-w(k,i)) goto 23
     if (thko(kp,i).gt.0.) then
      dslab=min(slab+thko(kp,i),-w(k,i))	&
           -min(slab           ,-w(k,i))
      dx=dslab/thko(kp,i)
      fcdx=a(kp)				&
          +b(kp)*.5*(dx-1.)			& !  not needed in pcm
          +c(kp)*(.25-dx*(.5-dx*athird))	  !  not needed in pcm,plm
      amount=amount+fcdx*dslab
      slab=slab+dslab
     end if
     if (kp.lt.nvl) go to 24
    else if (w(k,i).gt.0.) then			! interface moves down in atm.
     amount=slab*fld(k,i)
     kp=k+1
 25  kp=kp-1
     if (slab.ge.w(k,i)) goto 23
     if (thko(kp,i).gt.0.) then
      dslab=min(slab+thko(kp,i), w(k,i))	&
           -min(slab           , w(k,i))
      dx=dslab/thko(kp,i)
      fcdx=a(kp)				&
          +b(kp)*.5*(1.-dx)			& !  not needed in pcm
          +c(kp)*(.25-dx*(.5-dx*athird))	  !  not needed in pcm,plm
      amount=amount+fcdx*dslab
      slab=slab+dslab
     end if
     if (kp.gt.2) go to 25
    end if
 23  vertfx(k,i)=w(k,i)*amount/slab
   end do			! vertical loop

   vertfx(nvl,i)=0.		!  don't allow flux thru top
   vertdv(1,i)=vertfx(1,i)
   do k=2,nvl
    vertdv(k,i)=vertfx(k,i)-vertfx(k-1,i)
   end do
  end do			! horiz. loop
!SMS$PARALLEL END

  bfore=0.
  after=0.
  bforek(:)=0.
  afterk(:)=0.

  if (diagno) then
!SMS$PARALLEL (dh,i) BEGIN
   do i=1,nip
    do k=1,nvl
     bforek(k)=bforek(k)+fld(k,i)*thko(k,i)*area(i)
    end do
   end do
!SMS$PARALLEL END
  end if			! diagno

! --- compute low-order & antidiffusive (high- minus low-order) fluxes

  vlume(:)=0.
  drift(:)=0.

!SMS$PARALLEL (dh,i) BEGIN
  do i=1,nip			! horiz. loop
   do k=1,nvl			! vertical loop
    fmx(k,i)=fld(k,i)
    fmn(k,i)=fld(k,i)

    do edg=1,nprox(i)
     ix=prox(edg,i)
! +-+ +-+ +-+ +-+ +-+ +-+ +-+ +-+ +-+ +-+ +-+ +-+ +-+ +-+ +-+ +-+ +-+ +-+
!!   if (u(k,edg,i).ge.0.) then		! out-of cell > 0
!!    q=fld(k,i )
!!   else
!!    q=fld(k,ix)
!!   end if
!!   flx(k,edg,i)=u(k,edg,i)*q		! low order
! +-+ +-+ +-+ +-+ +-+ +-+ +-+ +-+ +-+ +-+ +-+ +-+ +-+ +-+ +-+ +-+ +-+ +-+
     edx=proxs(edg,i)		! index of joint edge as seen by neighbr
     flx(k,edg,i)=0.5*(		&	! low-order (out-of cell > 0)
         (u(k,edg,i )+abs(u(k,edg,i )))*fld(k,i )		&
       - (u(k,edx,ix)+abs(u(k,edx,ix)))*fld(k,ix))
! +-+ +-+ +-+ +-+ +-+ +-+ +-+ +-+ +-+ +-+ +-+ +-+ +-+ +-+ +-+ +-+ +-+ +-+
     q=.5*(fld(k,i)+fld(k,ix))		! high (2nd) order
     uan(k,edg,i)=q*u(k,edg,i)-flx(k,edg,i)	! antidiffusive

     fmx(k,i)=max(fmx(k,i),fld(k,ix))
     fmn(k,i)=min(fmn(k,i),fld(k,ix))
    end do			! loop over edges
    if (k.lt.nvl) then
     if (w(k  ,i).lt.0.) then
      fmx(k,i)=max(fmx(k,i),vertfx(k,i)/w(k,i))
      fmn(k,i)=min(fmn(k,i),vertfx(k,i)/w(k,i))
     end if
    end if
    if (k.gt.1) then
     if (w(k-1,i).gt.0.) then
      fmx(k,i)=max(fmx(k,i),vertfx(k-1,i)/w(k-1,i))
      fmn(k,i)=min(fmn(k,i),vertfx(k-1,i)/w(k-1,i))
     end if
    end if
   end do			! vertical loop
  end do			! horiz. loop
!SMS$PARALLEL END

!SMS$PARALLEL (dh,i) BEGIN
  do i=1,nip			! horiz. loop
   do k=1,nvl			! vertical loop
    hordiv(k,i)=0.
    do edg=1,nprox(i)
     hordiv(k,i)=hordiv(k,i)+flx(k,edg,i)
    end do
    hordiv(k,i)=hordiv(k,i)*rarea(i)

    q=fld(k,i)*thko(k,i)-hordiv(k,i)-vertdv(k,i)
    amount=max(0.,fmn(k,i)*thkn(k,i),min(q,fmx(k,i)*thkn(k,i)))
    if (recovr) then
     vlume(k)=vlume(k)+area(i)*thkn(k,i)
     drift(k)=drift(k)+(q-amount)*area(i)
    end if
    fld(k,i)=(fld(k,i)*onemu+amount)/(onemu+thkn(k,i))
   end do			! vertical loop
  end do			! horiz. loop

! --- at each grid point, determine the ratio of the largest permissible
! --- pos. (neg.) change in -fld- to the sum of all incoming (outgoing) fluxes

  do i=1,nip			! horiz. loop
   do k=1,nvl			! vertical loop
    totin=0.
    totou=0.
    do edg=1,nprox(i)
     totin=totin-min(0.,uan(k,edg,i))
     totou=totou+max(0.,uan(k,edg,i))
    end do
    flp(k,i)=(fmx(k,i)-fld(k,i))*thkn(k,i)/(totin+epsil)*rarea(i)
    fln(k,i)=(fld(k,i)-fmn(k,i))*thkn(k,i)/(totou+epsil)*rarea(i)
   end do			! vertical loop
  end do			! horiz. loop
!SMS$PARALLEL END

! --- limit antidiffusive fluxes

!SMS$EXCHANGE(flp,fln)

!SMS$PARALLEL (dh,i) BEGIN
  do i=1,nip			! horiz. loop
   do k=1,nvl			! vertical loop
    do edg=1,nprox(i)
     ix=prox(edg,i)
     clip=1.
     if (uan(k,edg,i).ge.0.) then		! out-of cell > 0
      clip=min(1.,fln(k,i),flp(k,ix))
     else
      clip=min(1.,flp(k,i),fln(k,ix))
     end if
     flx(k,edg,i)=uan(k,edg,i)*clip
    end do

    hordiv(k,i)=0.
    do edg=1,nprox(i)
     hordiv(k,i)=hordiv(k,i)+flx(k,edg,i)
    end do
    hordiv(k,i)=hordiv(k,i)*rarea(i)
    q=fld(k,i)*thkn(k,i)-hordiv(k,i)
    amount=max(0.,fmn(k,i)*thkn(k,i),min(q,fmx(k,i)*thkn(k,i)))
    if (recovr) drift(k)=drift(k)+(q-amount)*area(i)
    fld(k,i)=(fld(k,i)*onemu+amount)/(onemu+thkn(k,i))
   end do			! vertical loop
  end do			! horiz. loop
!SMS$PARALLEL END

  if (recovr) then

! --- recover 'clipped' amount and return to field layer by layer
   do k=1,nvl			! vertical loop
    var1=vlume(k)
    var2=drift(k)
!SMS$REDUCE(var1,var2,SUM)
    vlume(k)=var1
    drift(k)=var2

    if (vlume(k).ne.0.) then
     drift(k)=drift(k)/vlume(k)
     write (*,'(i2,a,es11.3)') k,'  tracer drift in fct3d',-drift(k)
    end if
   end do			! vertical loop

!SMS$PARALLEL (dh,i) BEGIN
   do i=1,nip
    do k=1,nvl
     fld(k,i)=fld(k,i)+drift(k)
    end do
   end do
!SMS$PARALLEL END
  end if			! recovr

  if (diagno) then
!SMS$PARALLEL (dh,i) BEGIN
   do i=1,nip
    do k=1,nvl
     afterk(k)=afterk(k)+fld(k,i)*thkn(k,i)*area(i)
    end do
   end do
!SMS$PARALLEL END

   var1=bforek(k)
   var2=afterk(k)
!SMS$REDUCE(var1,var2,SUM)
   bforek(k)=var1
   afterk(k)=var2

   do k=1,nvl
    bfore=bfore+bforek(k)
    after=after+afterk(k)
    print *,'bforek(k),afterk(k) in lyr',k,bforek(k),afterk(k)
   end do			! vertical loop
  end if			! diagno

  if (diagno) then
   if (bfore.ne.0.)						&
    write (*,'(a,1p,3e14.6,e11.1)') 'fct3d conservation:',	&
     bfore,after,after-bfore,(after-bfore)/bfore
  end if			! diagno

  if (recovr) then
   q=1.
   if (after.ne.0.) q=bfore/after
   write (*,'(a,f11.6)') 'fct3d: multiply tracer field by',q
   if (q.gt.1.1 .or. q.lt..9) stop '(excessive nonconservation)'

!SMS$PARALLEL (dh,i) BEGIN
   do i=1,nip
    do k=1,nvl
     fld(k,i)=fld(k,i)*q
    end do
   end do
!SMS$PARALLEL END
  end if			! recovr

  fields(:,:,n)=fld(:,:)

  if (diagno) then
   do k=1,nvl,7
    write (string,'(a,i2.2,a,i3)') 'fct3d out: fld',n,' lyr',k
    call findmxmn3(fields,nvl,nip,numfld,k,n,trim(string))
   end do
  end if

1 continue			! loop over fields
! print *,'exiting fct3d...'

  return
  end subroutine fct3d
end module module_fct3d
