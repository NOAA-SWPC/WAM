module module_globsum
  use module_control,   only: nvl, nvlp1, nip, dt, ntra, ntrb
  use module_constants, only: grvity, area, cp, rd

  implicit none

  real    :: qmstrold   ! save previous value
  real    :: qmstrcold  ! save previous value
  real    :: qmstrnold  ! save previous value
  logical :: qdtr_set = .false.

!JR These things were moved from dyn_run to here for restart capability.

  real :: qmass
  real :: qmsqv
  real :: qmsqw
  real :: qmste
  real :: qmstr  = 0.
  real :: qmstrn = 0.
  real :: qmstrc = 0.
  real :: qdtr
  real :: qdtrn
  real :: qdtrc

contains

!*********************************************************************
!     globsum
!       Calculate global sums of various useful quantities
!       Alexander E. MacDonald  11/14/2005
!       J. Lee                  01/04/2006
!       J. Rosinski             10/03/2011
!         Modified for restart capability
!*********************************************************************

  subroutine globsum (its, dp3d, tr, rn2d, rc2d, pr3d, ex3d, qf2d, qtrcr)
    integer, intent(in) :: its
!SMS$DISTRIBUTE (dh,nip) BEGIN
    real, intent(in) :: dp3d(nvl,nip)
    real, intent(in) :: tr(nvl,nip,ntra+ntrb)
    real, intent(in) :: rn2d(nip)         ! from PHY via CPL
    real, intent(in) :: rc2d(nip)         ! from PHY via CPL
    real, intent(in) :: pr3d(nvlp1,nip)
    real, intent(in) :: ex3d(nvlp1,nip)
    real, intent(in) :: qf2d(nip)         ! from PHY via CPL
!SMS$DISTRIBUTE END
    real, intent(out) :: qtrcr(ntra+ntrb)

! Local variables
    integer :: ipn  ! Index for icos point number
    integer :: k    ! Index vertical level
    integer :: t    ! Index for tracer type
    real    :: den
    real*8  :: summ,smqv,smqw,sumr,sume,sumrc,sumrn,sumtr(ntra+ntrb)

    summ  = 0.
    smqv  = 0.
    smqw  = 0.
    sumr  = 0.
    sume  = 0.
    sumrc = 0.
    sumrn = 0.
    sumtr(:) = 0.
!SMS$PARALLEL (dh,ipn) BEGIN
    do ipn=1,nip
      do k=1,nvl
        summ = summ + area(ipn)*dp3d(k,ipn)               ! integrated mass
        smqv = smqv + area(ipn)*dp3d(k,ipn)*tr(k,ipn,2)   ! integrated water vapor
        smqw = smqw + area(ipn)*dp3d(k,ipn)*tr(k,ipn,3)   ! integrated condensate
      end do

      sumr  = sumr + area(ipn)*rn2d(ipn)                  ! integrated precipitation
      sumrc = sumrc + area(ipn)*rc2d(ipn)                 ! integrated sub-gridscale (conv) precip
      sumrn = sumrn + area(ipn)*(rn2d(ipn)-rc2d(ipn))     ! integrated resolved precip
      den   = pr3d(1,ipn)/((rd/cp)*ex3d(1,ipn)*tr(1,ipn,1))
      sume  = sume+area(ipn)*den*qf2d(ipn)*dt             ! integrated evaporation

      do t=1,ntra+ntrb
        do k=1,nvl
          sumtr(t) = sumtr(t) + area(ipn)*dp3d(k,ipn)*tr(k,ipn,t) ! tracer
        end do
      end do
    end do           ! horizontal loop 

    summ = summ/grvity
    smqv = smqv/grvity
    smqw = smqw/grvity
!sms$reduce(summ,smqv,smqw,sumr,sume,sumrc,sumrn,SUM)
!SMS$PARALLEL END

    qmass = summ 
    qmsqv = smqv
    qmsqw = smqw
    qmste = sume
    qtrcr(:) = sumtr(:)
! save previous values
    qmstrold  = qmstr
    qmstrcold = qmstrc
    qmstrnold = qmstrn
! store new values
    qmstr  = sumr
    qmstrc = sumrc
    qmstrn = sumrn

    if (qdtr_set) then
      qdtr  = qmstr  - qmstrold
      qdtrn = qmstrn - qmstrnold
      qdtrc = qmstrc - qmstrcold
    else
      qdtr  = 0.
      qdtrn = 0.
      qdtrc = 0.
      qdtr_set = .true.
    end if

    return
  end subroutine globsum
end module module_globsum
