! Routines shared between WRF physics and chemistry.  

module module_wrf_share

contains


!*********************************************************************
subroutine wrf_set_array_bounds(nvl,nip,                 &
                                ids,ide,jds,jde,kds,kde, &
                                ims,ime,jms,jme,kms,kme, & 
                                its,ite,jts,jte,kts,kte)
!       Sets up WRF domain, memory, and tile bounds.  
!       Tom Henderson           July, 2009
!*********************************************************************

  implicit none

  ! Dummy arguments
  integer, intent(in   ) :: nvl,nip
  integer, intent(  out) :: ids,ide,jds,jde,kds,kde
  integer, intent(  out) :: ims,ime,jms,jme,kms,kme
  integer, intent(  out) :: its,ite,jts,jte,kts,kte

  ! locals
!SMS$DISTRIBUTE (dh,nip) BEGIN
  integer :: itmp(nip)
!SMS$DISTRIBUTE END

  ! set up trivial indices
  ims = 1
  ime = 1
  ids = 1
  ide = 1
  its = 1
  ite = 1
  jds = 1
  kms = 1
  kds = 1
  kts = 1

  ! set up non-parameter indices
  jde = nip    ! domain end in decomposed dimension
!SMS$PARALLEL (dh,j) BEGIN
!SMS$TO_LOCAL (<1,jds:lbound>, <1,jde:ubound>) BEGIN
  jts = jds
  jte = jde
!SMS$TO_LOCAL END
!SMS$PARALLEL END
  jms = LBOUND(itmp,1)
  jme = UBOUND(itmp,1)
  kme = nvl+1
  kde = nvl+1
  kte = nvl

  return

end subroutine wrf_set_array_bounds


end module module_wrf_share

