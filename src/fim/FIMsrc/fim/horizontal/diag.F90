module module_diag
contains
!*********************************************************************
!     diag
!	Diagnostic program after prognostic variables are calculated
!	Alexander E. MacDonald  11/14/2005
!	J.Lee                   01/04/2006
!*********************************************************************

subroutine diag(its,  &
  ph3d,us3d,vs3d,     &     ! phi (=g*z), west wind and south wind on s
  ex3d,mp3d,dp3d,     &     ! pres,exner,mont pot,kin energy
  tr,trdp        )     ! specific humidity
use module_control  ,only: nvl,nvlp1,nip,ptop,ntra,ntrb
use module_constants,only: p1000,cp,rd,qvmin,qwmin
implicit none

! Dimension and type external varialbles:
integer,intent (IN) :: its
!SMS$DISTRIBUTE (dh,nip) BEGIN
real,intent (IN)    :: us3d(nvl,nip)	! west wind
real,intent (IN)    :: vs3d(nvl,nip)	! south wind
real,intent (INOUT) :: ph3d(nvlp1,nip)  ! phi (=gz)
real,intent (IN)    :: ex3d(nvlp1,nip)	! exner
real,intent (OUT)   :: mp3d(nvl,nip)	! montgomery potential
real,intent (IN)    :: dp3d(nvl,nip)     ! specific humidity
real,intent (INOUT) :: tr  (nvl,nip,ntra+ntrb)     ! specific humidity
real,intent (INOUT) :: trdp(nvl,nip,ntra+ntrb)     ! specific humidity
!SMS$DISTRIBUTE END

! Declare local variables:
integer    :: ipn  ! Index for icos point number
integer    :: ivl  ! Index vertical level
real       :: totp ! Total pressure summation variable
real       :: temp1(nvl),temp2(nvl)

!  Note that tr (tracers), velocity, us3d,vs3d and montgomery potential, mp3d
!  are all constant through the layers.  Phi (ph3d) and pressure (dp3d)
!  vary through the layer.

!  Layer variables:  us3d,vs3d,tr (tracers),mp3d
!  Level variables:  ex3d,ph3d

!SMS$PARALLEL (dh,ipn) BEGIN
!sms$compare_var(ex3d, "diag.F90 - ex3d5 ")
!sms$compare_var(ph3d, "diag.F90 - ph3d5 ")
do ipn=1,nip	!  global icos loop

!.........................................................
!   Determine bottom layer values
!.........................................................

  mp3d(1,ipn)=ex3d(1,ipn)*tr(1,ipn,1) + ph3d(1,ipn) ! mp at surface

  do ivl=2,nvl	!  vertical loop

    ! Hydrostatic eqn:  d mp/d theta = exner
    temp1(ivl)=ex3d(ivl,ipn)*(tr(ivl,ipn,1)-tr(ivl-1,ipn,1))

    ! Hydrostatic eqn:  d phi/d exner = - theta
    temp2(ivl)=tr(ivl,ipn,1)*(ex3d(ivl,ipn)-ex3d(ivl-1,ipn))

  enddo		! vertical loop
  do ivl=2,nvl	!  vertical loop

    ! Hydrostatic eqn:  d mp/d theta = exner
    mp3d(ivl,ipn)=mp3d(ivl-1,ipn)+temp1(ivl)

    ! Hydrostatic eqn:  d phi/d exner = - theta
    ph3d(ivl,ipn)=ph3d(ivl-1,ipn)-temp2(ivl)

  enddo		! vertical loop

  do ivl=1,nvl
    tr(ivl,ipn,2)=max(qvmin,tr(ivl,ipn,2))
    trdp(ivl,ipn,2)=tr(ivl,ipn,2)*dp3d(ivl,ipn)
    tr(ivl,ipn,3)=max(qwmin,tr(ivl,ipn,3))
    trdp(ivl,ipn,3)=tr(ivl,ipn,3)*dp3d(ivl,ipn)
  end do

enddo		! horizontal loop 

!SMS$PARALLEL END

return
end subroutine diag
end module module_diag
