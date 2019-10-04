module module_diagnoise
use findmaxmin1
contains
!*********************************************************************
!       diagnoise
!	A routine to diagnose external gravity wave noise in FIM
!       (expressed as rms of 2nd time derivative of surface pressure)
!	R. Bleck                July 2006
!*********************************************************************

subroutine diagnoise    &
(its,                   &	! model time step
ptdcy                   )       ! sfc.pres. tendency at 2 consec. time levels

use module_control,only: nvl,nip
implicit none

! Dimension and type external variables:
integer,intent (IN)  :: its		! model time step
!SMS$DISTRIBUTE(dh,NIP) BEGIN
real,   intent (IN)    :: ptdcy (nip,2) ! d(surf.prs.)/dt

! Local variables:
real    :: work(nip)
!SMS$DISTRIBUTE END
real    :: deriv2		! rms of 2nd derivative of surface pressure
real    :: derivma		! mean abs of 2nd derivative of surface pressure
real    :: scale = 5.		! arbitrary scale for plotting noise parameter
integer :: ivl			! layer index
integer :: ipn			! index for icosahedral grid
real    :: valmin,valmax

if (its.lt.2) return

deriv2=0.
derivma=0.

!SMS$PARALLEL (dh,ipn) BEGIN

do ipn=1,nip
  deriv2=deriv2+(ptdcy(ipn,1)-ptdcy(ipn,2))**2
  derivma=derivma+abs(ptdcy(ipn,1)-ptdcy(ipn,2))
end do
!sms$reduce(deriv2,derivma,SUM)

!SMS$PARALLEL END

! --- plot time series of noise index in stdout
! --- (type 'grep =+= stdout' to display the time series)
deriv2=sqrt(deriv2/nip)
derivma=(derivma/nip)
call linout(scale * deriv2,'x',its)
write (6,*)'rms-d(psfc)**2/d**2t) =',its,deriv2
write (6,*)'deriv-mean abs        =',its,derivma

!! valmin=minval(ptdcy(:,1))
!! valmax=maxval(ptdcy(:,1))
!! !SMS$REDUCE(valmin,min)
!! !SMS$REDUCE(valmax,max)
!! print *,'min/max of ptdcy(:,1):',valmin,valmax
!! 
!! valmin=minval(ptdcy(:,2))
!! valmax=maxval(ptdcy(:,2))
!! !SMS$REDUCE(valmin,min)
!! !SMS$REDUCE(valmax,max)
!! print *,'min/max of ptdcy(:,2):',valmin,valmax

work(:)=ptdcy(:,1)
call findmxmn1(work,nip,'ptdcy(:,1)')
work(:)=ptdcy(:,2)
call findmxmn1(work,nip,'ptdcy(:,2)')
 
return
end subroutine diagnoise


subroutine linout(value,char,labl)
!
! --- print single characters in a manner mimicking a curve plot in x,y space
! --- abscissa: down the page; ordinate: across the page
!
implicit none
real,intent(IN)        :: value
character*1,intent(IN) :: char
integer,intent(IN)     :: labl

integer,parameter      :: length=72
character*1 line(length)
integer l,n
 
! --- replace n-th element of array 'line' by character 'char', where
! --- n = 'value' modulo 'length'
! --- initialize 'line' by blanks before adding 'char'
! --- output 'line' after adding 'char'
! --- labl       -- abscissa value (integer), added to output line
!
do l=1,length
  line(l)=' '
end do
if (value.gt.0.) then
  n=int(mod(value+float(length-1),float(length)))+1
  line(n)=char
end if
write (*,'(''=+='',i6,80a1)') labl,(line(l),l=1,length)
return
end subroutine linout

end module module_diagnoise
