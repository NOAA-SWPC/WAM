  module module_gfs_tropp
!
! Abstract: This module contains routines needed for compute 
!           tropopause level fields
!
  use MODULE_gfs_machine, only:kint_mpi
  use MODULE_gfs_physcons, only: con_rog

  implicit none

  private
  public :: tpause

  contains

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
subroutine rsearch1(km1,z1,km2,z2,l2)
!$$$  subprogram documentation block
!
! subprogram:    rsearch1    search for a surrounding real interval
!   prgmmr: iredell    org: w/nmc23     date: 98-05-01
!
! abstract: this subprogram searches a monotonic sequences of real numbers
!   for intervals that surround a given search set of real numbers.
!   the sequences may be monotonic in either direction; the real numbers
!   may be single or double precision.
!
! program history log:
! 1999-01-05  mark iredell
!
! usage:    call rsearch1(km1,z1,km2,z2,l2)
!   input argument list:
!     km1    integer number of points in the sequence
!     z1     real (km1) sequence values to search
!            (z1 must be monotonic in either direction)
!     km2    integer number of points to search for
!     z2     real (km2) set of values to search for
!            (z2 need not be monotonic)
!
!   output argument list:
!     l2     integer (km2) interval locations from 0 to km1
!            (z2 will be between z1(l2) and z1(l2+1))
!
! subprograms called:
!   sbsrch essl binary search
!   dbsrch essl binary search
!
! remarks:
!   returned values of 0 or km1 indicate that the given search value
!   is outside the range of the sequence.
!
!   if a search value is identical to one of the sequence values
!   then the location returned points to the identical value.
!   if the sequence is not strictly monotonic and a search value is
!   identical to more than one of the sequence values, then the
!   location returned may point to any of the identical values.
!
!   if l2(k)=0, then z2(k) is less than the start point z1(1)
!   for ascending sequences (or greater than for descending sequences).
!   if l2(k)=km1, then z2(k) is greater than or equal to the end point
!   z1(km1) for ascending sequences (or less than or equal to for
!   descending sequences).  otherwise z2(k) is between the values
!   z1(l2(k)) and z1(l2(k+1)) and may equal the former.
!
! attributes:
!   language: fortran
!
!$$$
  implicit none
  integer,intent(in):: km1,km2
  real,intent(in):: z1(km1),z2(km2)
  integer,intent(out):: l2(km2)
  integer k1,k2
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  find the surrounding input interval for each output point.
  if(z1(1).le.z1(km1)) then
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  input coordinate is monotonically ascending.
    do k2=1,km2
      l2(k2)=km1
      do k1=1,km1-1
        if(z1(k1)<=z2(k2).and.z1(k1+1)>z2(k2)) then
          l2(k2)=k1
          exit
        endif
      enddo
    enddo
  else
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  input coordinate is monotonically descending.
    do k2=1,km2
      l2(k2)=km1
!      do k1=1,km1-1
!        if(z1(k1)<z2(k2)) then
      do k1=km1,2,-1
        if(z2(k2)>=z1(k1).and.z2(k2)<z1(k1-1))then
          l2(k2)=k1-1
          exit
        endif
      enddo
    enddo
  endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
end subroutine
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  subroutine tpause(km,p,u,v,t,h,ptp,utp,vtp,ttp,htp,shrtp)
!$$$  Subprogram documentation block
!
! Subprogram: tpause     Compute tropopause level fields
!   Prgmmr: Iredell      Org: np23        Date: 1999-10-18
!
! Abstract: This subprogram finds the tropopause level and computes fields 
!   at the tropopause level.  The tropopause is defined as the lowest level
!   above 500 mb which has a temperature lapse rate of less than 2 K/km.
!   The lapse rate must average less than 2 K/km over a 2 km depth.
!   If no such level is found below 50 mb, the tropopause is set to 50 mb.
!   The tropopause fields are interpolated linearly in lapse rate.
!   The tropopause pressure is found hydrostatically.
!   The tropopause wind shear is computed as the partial derivative
!   of wind speed with respect to height at the tropopause level.
!
! Program history log:
!   1999-10-18  Mark Iredell
!
! Usage:  call tpause(km,p,u,v,t,h,ptp,utp,vtp,ttp,htp,shrtp)
!   Input argument list:
!     km       integer number of levels
!     p        real (km) pressure (Pa)
!     u        real (km) x-component wind (m/s)
!     v        real (km) y-component wind (m/s)
!     t        real (km) temperature (K)
!     h        real (km) height (m)
!   Output argument list:
!     ptp      real tropopause pressure (Pa)
!     utp      real tropopause x-component wind (m/s)
!     vtp      real tropopause y-component wind (m/s)
!     ttp      real tropopause temperature (K)
!     htp      real tropopause height (m)
!     shrtp    real tropopause wind shear (1/s)
!
! Files included:
!   physcons.h     Physical constants
!
! Subprograms called:
!   rsearch1       search for a surrounding real interval
!
! Attributes:
!   Language: Fortran 90
!
!$$$
    implicit none
    integer,intent(in):: km
    real,intent(in),dimension(km):: p,u,v,t,h
    real,intent(out):: ptp,utp,vtp,ttp,htp,shrtp
    real,parameter:: ptplim(2)=(/500.e+2,50.e+2/),gamtp=2.e-3,hd=2.e+3
    real gamu,gamd,td,gami,wtp,spdu,spdd
    integer klim(2),k,kd,ktp, kd_array(1)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  find tropopause level
    call rsearch1(km-2,p(2),2,ptplim(1),klim(1))
!    klim(1)=klim(1)+2
    klim(1)=klim(1)+1
!    klim(2)=klim(2)+1
    klim(2)=klim(2)+2
    gamd=1.e+9
    ktp=klim(2)
    wtp=0
!    do k=klim(1),klim(2)
    do k=klim(1),klim(2),-1
!      gamu=(t(k-1)-t(k+1))/(h(k+1)-h(k-1))
      gamu=(t(k+1)-t(k-1))/(h(k-1)-h(k+1))
      if(gamu.le.gamtp) then
!        call rsearch1(km-k-1,h(k+1),1,h(k)+hd,kd)
	call rsearch1(k-2,h(2),1,(/h(k)+hd/),kd_array)
        kd = kd_array(1)
!        td=t(k+kd)+(h(k)+hd-h(k+kd))/(h(k+kd+1)-h(k+kd))*(t(k+kd+1)-t(k+kd))
	td=t(kd+2)+(h(k)+hd-h(2+kd))/(h(kd+1)-h(2+kd))*(t(kd+1)-t(2+kd))
        gami=(t(k)-td)/hd
        if(gami.le.gamtp) then
          ktp=k
          wtp=(gamtp-gamu)/(max(gamd,gamtp+0.1e-3)-gamu)
          exit
        endif
      endif
      gamd=gamu
    enddo
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  compute tropopause level fields
    utp=u(ktp)-wtp*(u(ktp)-u(ktp-1))
    vtp=v(ktp)-wtp*(v(ktp)-v(ktp-1))
    ttp=t(ktp)-wtp*(t(ktp)-t(ktp-1))
    htp=h(ktp)-wtp*(h(ktp)-h(ktp-1))
    ptp=p(ktp)*exp((h(ktp)-htp)*(1-0.5*(ttp/t(ktp)-1))/(con_rog*t(ktp)))
    spdu=sqrt(u(ktp)**2+v(ktp)**2)
    spdd=sqrt(u(ktp-1)**2+v(ktp-1)**2)
    shrtp=(spdu-spdd)/(h(ktp)-h(ktp-1))
    
    utp=u(ktp)-wtp*(u(ktp)-u(ktp+1))
    vtp=v(ktp)-wtp*(v(ktp)-v(ktp+1))
    ttp=t(ktp)-wtp*(t(ktp)-t(ktp+1))
    htp=h(ktp)-wtp*(h(ktp)-h(ktp+1))
    ptp=p(ktp)*exp((h(ktp)-htp)*(1-0.5*(ttp/t(ktp)-1))/(con_rog*t(ktp)))
    spdu=sqrt(u(ktp)**2+v(ktp)**2)
    spdd=sqrt(u(ktp+1)**2+v(ktp+1)**2)
    shrtp=(spdu-spdd)/(h(ktp)-h(ktp+1))
  end subroutine

  end module module_gfs_tropp
