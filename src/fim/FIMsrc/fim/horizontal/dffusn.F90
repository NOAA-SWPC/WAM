module module_dffusn_lev
contains
!*********************************************************************
!     dffusn_lev
!	Diffuse level variable (no thickness weighting)
!	S. Sun                           September 2009
!*********************************************************************

  subroutine dffusn_lev (fld, dfflen, kdim, k1, k2)
  use module_control  ,only: npp,nip,PrintIpnDiag
  use module_constants,only: nprox,prox,rarea,sideln

  implicit none
!  External type and dimension:

  integer,intent (IN)    :: kdim		! vert.dim. of fld and dfflen
  integer,intent (IN)    :: k1,k2		! operate on levels k1 ... k2
  real   ,intent (IN)    :: dfflen (kdim)	! diffusion length scale (m)
!SMS$DISTRIBUTE (dh,nip) BEGIN
  real   ,intent (INOUT) :: fld (kdim,nip)	! field(s) to be diffused

! Local variables:
  real    :: flxdv(kdim,nip)	! line integral of dffus.flux across the edge
!SMS$DISTRIBUTE END

! Local variables
  integer   :: k		! layer index
  integer   :: ipn		! icos point index
  integer   :: ipx		! neighbor across joint edge
  integer   :: edg		! icos edge index
  real      :: factor

!SMS$PARALLEL (dh,ipn) BEGIN
!SMS$EXCHANGE(fld)

  do ipn = 1,nip			! horizontal loop
   flxdv(:,ipn) = 0.
   do edg = 1,nprox(ipn)
    ipx = prox(edg,ipn)			! neighbor across shared edge
    do k = k1,k2			! loop through levels
     flxdv(k,ipn) = flxdv(k,ipn)+(fld(k,ipn)-fld(k,ipx))		&
        *sideln(edg,ipn)

     if (ipn.eq.PrintIpnDiag .and. mod(k,7).eq.k1) then
!SMS$IGNORE BEGIN
      write (*,'(a,i8,i4,a,2f12.1,a,es11.2)') 'ipn,k=',ipn,k,		&
      '  (dffusn_lev) fld=',fld(k,ipn),fld(k,prox(edg,ipn)),'  flx=',	&
       fld(k,ipn)-fld(k,prox(edg,ipn))
!SMS$IGNORE END
     end if

    end do				! loop through levels
   end do				! loop through edges
  end do				! horizontal loop

  do ipn = 1,nip			! horizontal loop
   do k = k1,k2				! loop through levels
   factor = -dfflen(k)*rarea(ipn)
    fld(k,ipn) = fld(k,ipn) + flxdv(k,ipn)*factor

    if (ipn.eq.PrintIpnDiag .and. mod(k,7).eq.k1) then
!SMS$IGNORE BEGIN
     write (*,'(i8,i4,a,2es11.2,f12.1)') ipn,k,				&
      '  (dffusn_lev) flxdv,fac,fld=',flxdv(k,ipn),factor,fld(k,ipn)
!SMS$IGNORE END
    end if

   end do				! loop through levels
  end do				! horizontal loop

!!!SMS$EXCHANGE(fld)
!SMS$PARALLEL END

  return
  end subroutine dffusn_lev
end module module_dffusn_lev	! SMS doesn't like multiple routines in module


module module_dffusn_lyr
contains
!*********************************************************************
!     dffusn_lyr
!	Diffuse layer variable (thickness-weighted for conservation)
!	S. Sun                           September 2009
!*********************************************************************

  subroutine dffusn_lyr (fld, delp, dfflen)
  use module_control  ,only: npp,nip,nvl,PrintIpnDiag
  use module_constants,only: nprox,prox,rarea,sideln

  implicit none
!  External type and dimension:
  real   ,intent (IN)    :: dfflen              ! diffusion length scale (m)
!SMS$DISTRIBUTE (dh,nip) BEGIN
  real   ,intent (INOUT) :: fld   (nvl,nip)     ! field to be diffused
  real   ,intent (IN)    :: delp  (nvl,nip)     ! lyr thknss, Pa

! Local variables:
  real    :: flxdv(nvl,nip)	! line integral of dffus.flux across the edge
!SMS$DISTRIBUTE END

! Local variables
  integer   :: k		! layer index
  integer   :: ipn		! icos point index
  integer   :: ipx		! neighbor across joint edge
  integer   :: edg		! icos edge index
  real      :: factor
  real,parameter   :: thshld = 1.e-11
  real      :: hfharm,a,b
  hfharm(a,b) = a*b/(a+b)	! harmonic average x 0.5

!SMS$PARALLEL (dh,ipn) BEGIN
!SMS$EXCHANGE(delp,fld)

  do ipn = 1,nip			! horizontal loop
   flxdv(:,ipn) = 0.
   do edg = 1,nprox(ipn)
    ipx = prox(edg,ipn)			! neighbor across shared edge
    do k = 1,nvl			! loop through layers
     flxdv(k,ipn) = flxdv(k,ipn)+(fld(k,ipn)-fld(k,ipx))		&
        *sideln(edg,ipn)*2.*hfharm(max(delp(k,ipn),thshld)		&
                                  ,max(delp(k,ipx),thshld))

     if (ipn.eq.PrintIpnDiag .and. mod(k,7).eq.1) then
!SMS$IGNORE BEGIN
      write (*,'(a,i8,i4,a,2f12.1,a,es11.2)') 'ipn,k=',ipn,k,		&
      '  (dffusn_lyr) fld=',fld(k,ipn),fld(k,prox(edg,ipn)),'  flx=',	&
       fld(k,ipn)-fld(k,prox(edg,ipn))
!SMS$IGNORE END
     end if

    end do				! loop through layers
   end do				! loop through edges
  end do				! horizontal loop

  do ipn = 1,nip			! horizontal loop
   do k = 1,nvl				! loop through layers
    factor = -dfflen*rarea(ipn)/max(delp(k,ipn),thshld)
    fld(k,ipn) = fld(k,ipn) + flxdv(k,ipn)*factor

    if (ipn.eq.PrintIpnDiag .and. mod(k,7).eq.1) then
!SMS$IGNORE BEGIN
     write (*,'(i8,i4,a,2es11.2,f12.1)')ipn,k,				&
      '  (dffusn_lyr) flxdv,fac,fld=',flxdv(k,ipn),factor,fld(k,ipn)
!SMS$IGNORE END
    end if

   end do				! loop through layers
  end do				! horizontal loop

!!!SMS$EXCHANGE(fld)
!SMS$PARALLEL END

  return
  end subroutine dffusn_lyr
end module module_dffusn_lyr	! SMS doesn't like multiple routines in module
