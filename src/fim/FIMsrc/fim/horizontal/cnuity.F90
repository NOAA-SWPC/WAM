module module_cnuity
use stencilprint
use stenedgprint
use findmaxmin2
contains
!*********************************************************************
!     cnuity
! 	based on fct = flux corrected transport
!	J. Lee		Author				September, 2005
!	A. E. MacDonald Documentor			November,  2005
!	R. Bleck	major rewrite			April      2008
!	R. Bleck	revised omega diagnostics	August     2009
!	R. Bleck	discarded-flux recovery		November   2009
!
! 	This routine is based on Zalesak, JOURNAL OF COMPUTATIONAL
!       PHYSICS, 31, 335-362, 1979.  Dale Durran provides an
!       excellent discussion of flux corrected transport in his book
!     NUMERICAL METHODS FOR WAVE EQUATIONS IN GEOPHYSICAL FLUID DYNAMICS.
!*********************************************************************

  subroutine cnuity(its,		&
    nf,of,vof,				&
    adbash1,adbash2,adbash3,		&
    u_velo,v_velo,			&
    u_edg,v_edg,			&
    dp_edg,lp_edg,			&
    delp,pres,exner,			&
    dp_tndcy,dplo_tndcy,		&
    massfx,omega,			&
    tcnuity,tcnuityEx,tcnuityBa,	&
    TimingBarriers )

use module_control  ,only: nvl,nvlp1,npp,nip,nabl,dt,nd,PrintIpnDiag
use module_constants,only: nprox,prox,proxs,sidevec_c,		&
                           sidevec_e,rarea,area,p1000,rd,cp,    &
                           nedge,permedge
implicit none

!..............................................................
!	Sec. 0  Dimension and Type
!..............................................................

! External variables:
  integer,intent (IN)    :: its              ! model time step
  integer,intent (IN)    :: nf,of,vof	     ! time slots: new,old,very old
  real   ,intent (IN)    :: adbash1,adbash2,adbash3
  real*8 ,intent (INOUT) :: tcnuity          ! computation timer
  real*8 ,intent (INOUT) :: tcnuityEx        ! halo update communication timer
  real*8 ,intent (INOUT) :: tcnuityBa        ! barrier timer for task skew
  logical,intent (IN)    :: TimingBarriers   ! measure task skew when .true.
!SMS$DISTRIBUTE (dh,nip) BEGIN
  real   ,intent (IN)    :: u_velo    (nvl,nip)		! west wind, m/sec
  real   ,intent (IN)    :: v_velo    (nvl,nip)		! south wind,m/sec
  real   ,intent (IN)    :: u_edg     (nvl,npp,nip)	! u on edges, m/sec
  real   ,intent (IN)    :: v_edg     (nvl,npp,nip)	! v on edges, m/sec
  real   ,intent (IN)    :: dp_edg    (nvl,npp,nip)	! dp on edges, Pa
  real   ,intent (IN)    :: lp_edg    (nvl,npp,nip)	! midlyr p on edges, Pa
  real   ,intent (INOUT) :: delp      (nvl,nip)		! lyr thknss, Pa
  real   ,intent (INOUT) :: pres      (nvlp1,nip)	! prs on intfc, Pa
  real   ,intent (INOUT) :: exner     (nvlp1,nip)	! Exner fct
  real   ,intent (INOUT) :: dp_tndcy  (nvl,nip,nabl)	! Pa/sec
  real   ,intent (INOUT) :: dplo_tndcy(nvl,nip,nabl)	! Pa/sec
  real   ,intent (INOUT) :: massfx    (nvl,npp,nip,3)	! N/sec
  real   ,intent (OUT)   :: omega     (nvl,nip)		! N/sec

! Local variables:
  real    :: p_plus(nvl,nip)	! Zalesak's p_plus, N/sec
  real    :: p_mnus(nvl,nip)	! Zalesak's p_minus, N/sec
  real    :: r_plus(nvl,nip)	! Zalesak's r_plus, dimensionless
  real    :: r_mnus(nvl,nip)	! Zalesak's r_minus, dimensionless
  real    :: vnorm(nvl,npp,nip)	! outward-directed velocity on edge x edge lngth
  real    :: antifx(nvl,npp,nip)	! N/sec
  real    :: delp_lo(nvl,nip)		! Pa
  real    :: anti_tndcy(nvl,nip)	! Pa/sec
  real    :: psurf(nip)			! surface pressure
  real    :: recovr(npp,nip)		! fluxes discarded in clipping process
!SMS$DISTRIBUTE END

  integer :: k	 		! layer index
  integer :: ipn		! Index for icos cell number
  integer :: edg		! Index for icos edge number
  real    :: flxhi              ! thickness flux, high-order
  real    :: dpmax,dpmin	! Zalesak's phi_max,phi_min
  real    :: q_plus,q_mnus      ! Zalesak's q_plus,q_mnus (N/sec)
  real    :: dpdx,dpdy		! pressure gradient in x,y direction
  real    :: coltend,lyrtend	! column & layer pressure tendency
  real    :: old,clip
  integer :: ipx		! neighbor across joint edge
  integer :: edx		! joint edge index as seen by neighbor
  character :: string*32
  real,   parameter :: thshld = 1.e-11
  logical,parameter :: low_ord = .false. ! if true, skip antidiffusive part
  integer :: edgcount           ! count of icos edges
  real*8  :: t1

!.............................................................
!  Sec. 1 Calculate Anti Diffusive Flux, Low order forcing
!.............................................................

! Calculates the FCT low and high order fluxes, and defines
! the antidiffusive flux as the difference between the high and
! low-order fluxes.  The low order flux is computed based on 
! the assumption of piecewise continuity, with a constant value 
! in each cell.  The higher order uses the "gazebo" with a
! sloped line used for the flux integral.

!! do k=1,nvl,7
!!  write (string,'(a,i2)') 'cnuity: old dp k=',k
!!  call findmxmn2(delp,nvl,nip,k,string)
!! end do

!SMS$PARALLEL (dh,ipn) BEGIN

!sms$compare_var(sidevec_e, "cnuity.F90 - sidevec_e6 ")
!sms$compare_var(u_edg    , "cnuity.F90 - u_edg6 ")
!sms$compare_var(v_edg    , "cnuity.F90 - v_edg6 ")

  call StartTimer(t1)

  ! Initialize these local and INTENT(OUT) arrays so COMPARE_VAR does not get 
  ! confused by uninitialized edg = 6 edges of pentagonal grid cells.  This 
  ! code could be safely omitted if COMPARE_VAR were not used.  

!SMS$HALO_COMP(<1,1>) BEGIN
  do ipn = 1,nip		! horizontal loop
   vnorm  (:,npp,ipn)    = 0.
   antifx (:,npp,ipn)    = 0.
   massfx (:,npp,ipn,nf) = 0.

   do edgcount = 1,nedge(ipn)	! loop through edges
    edg = permedge(edgcount,ipn)
    do k = 1,nvl		! loop through layers
     vnorm(k,edg,ipn) = sidevec_e(2,edg,ipn)*u_edg(k,edg,ipn)	&
                      - sidevec_e(1,edg,ipn)*v_edg(k,edg,ipn)  
    end do			! loop through layer
   end do			! loop through edges
  end do			! horizontal loop
!SMS$HALO_COMP END

! Avoid exchange via HALO_COMP in previous loop and in edgvar.F90
!!!SMS$EXCHANGE(vnorm,dp_edg)

!sms$compare_var(vnorm   , "cnuity.F90 - vnorm7 ")
!sms$compare_var(dp_edg  , "cnuity.F90 - dp_edg7 ")
!sms$compare_var(delp    , "cnuity.F90 - delp7 ")

!SMS$HALO_COMP(<1,1>) BEGIN
  do ipn = 1,nip		! horizontal loop
   do edgcount = 1,nedge(ipn)	! loop through edges
    edg = permedge(edgcount,ipn)
    ipx = prox(edg,ipn)		! neighbor across shared edge
    edx = proxs(edg,ipn)	! neighbor's index of shared edge

    do k = 1,nvl		! loop through layers
  
     ! high-order mass flux (2nd order centered, out-of cell > 0)
     flxhi = vnorm(k,edg,ipn)*dp_edg(k,edg,ipn)

     ! low-order mass flux (donor-cell, out-of cell > 0)
     massfx(k,edg,ipn,nf) = 0.5*(					&
            (vnorm(k,edg,ipn)+abs(vnorm(k,edg,ipn)))*delp(k,ipn)	&
          - (vnorm(k,edx,ipx)+abs(vnorm(k,edx,ipx)))*delp(k,ipx)  )

     ! anti-diffusive flux = difference between high order flux (flxhi)
     ! and low-order flux (massfx), from Zalesek, p336, Eqn (3):

     antifx(k,edg,ipn) = flxhi-massfx(k,edg,ipn,nf)   	! N/sec
  
    end do			! loop through layers
   end do			! loop through edges
  end do			! horizontal loop
!SMS$HALO_COMP END

  write (string,'(a,i6,2x)') 'step',its
  call stencl(delp,nvl,.01,string(1:12)//'(cnuity) old dp')
! call stenedg(delp,massfx(1,1,1,nf),nvl,				&
!              string(1:12)//'(cnuity) old dp, low-ord flx')

  do ipn = 1,nip		! horizontal loop
   dplo_tndcy(:,ipn,nf) = 0.
   do edg = 1,nprox(ipn)	! loop through edges
    do k = 1,nvl		! loop through layers

     ! sum up edge fluxes to get low-order tendency term 
     dplo_tndcy(k,ipn,nf) = dplo_tndcy(k,ipn,nf)+massfx(k,edg,ipn,nf)

    end do			! loop through layers
   end do			! loop through edges

!.............................................................
!  Sec. 2.  Calculate new low-order dp using full Adams Bashforth
!.............................................................

   do k = 1,nvl			! loop through layers

    ! divide tendency by cell area to convert to  Pa/sec
    dplo_tndcy(k,ipn,nf) = -dplo_tndcy(k,ipn,nf)*rarea(ipn)

    ! get new value for the low-order delp field using Adams Bashforth
    ! 3 time levels, the one just calculated (nf), and the two prev ones,
    ! marked of (old field) and vof (very old field):

    delp_lo(k,ipn) = delp(k,ipn)				&
         +adbash1*dplo_tndcy(k,ipn, nf)				&
         +adbash2*dplo_tndcy(k,ipn, of)				&
         +adbash3*dplo_tndcy(k,ipn,vof)

   end do			! loop through layers
  end do			! horizontal loop

  if (low_ord) then		! use low-order mass fluxes only

   do ipn = 1,nip		! horizontal loop
    do k = 1,nvl		! loop through layers
     delp(k,ipn) = delp_lo(k,ipn)
     dp_tndcy(k,ipn,nf) = dplo_tndcy(k,ipn,nf)
    end do
   end do

  else				! evaluate and apply antidiffusive fluxes

! Dale Durrans book indicates that condition Zal (14) should be
! satisfied, although Zalesak says its cosmetic.  We believe Dale:
! Also calculate p_plus Zalesak (7) and p_minus Zalesak (10)
! (aggregate of incoming and outgoing fluxes, N/sec) 

   call IncrementTimer(t1,tcnuity)

   if (TimingBarriers) then
    call StartTimer(t1)
!SMS$BARRIER
    call IncrementTimer(t1,tcnuityBa)
   endif

   call StartTimer(t1)
!SMS$EXCHANGE(delp_lo)
!!!SMS$EXCHANGE(delp) exchanged in edgvar
   call IncrementTimer(t1,tcnuityEx)

!sms$compare_var(antifx, "cnuity.F90 - antifx8 ")
!sms$compare_var(delp_lo,"cnuity.F90 - delp_lo8 ")

   call StartTimer(t1)

   do ipn = 1,nip		! horizontal loop
    p_plus(:,ipn) = 0. 
    p_mnus(:,ipn) = 0.
    do edg = 1,nprox(ipn)	! loop through edges
     do k = 1,nvl		! loop through layers
      if(antifx(k,edg,ipn).le.0.) then			! flux into cell
       p_plus(k,ipn) = p_plus(k,ipn)-antifx(k,edg,ipn)
      else						! flux out-off cell 
       p_mnus(k,ipn) = p_mnus(k,ipn)+antifx(k,edg,ipn)    
      endif
     end do			! loop through layers
    end do			! loop through edges

!   At this stage we have calculated the low-order delp_lo, and the 
!   unclipped antidiffusive flux for the entire icos global grid.

!............................................................
!   Sec 3.  Monotonicity Limit on Fluxes
!............................................................

!   Determine the amount of antidiffusive fluxes that can be
!   added to the low-order solution without violating monotonicity.

    do k = 1,nvl		! loop through layers

     ! According to Zal (17), we must limit according to max of
     ! dp from any gazebo direction, in either current or prev time
     ! step.  Likewise for the minimum.
     ! For the pentagons prox(6,ipn) is set to prox(5,ipn) in init.F90.
     
     dpmax = max(delp_lo(k,ipn),delp(k,ipn),			&
           delp_lo(k,prox(1,ipn)),delp_lo(k,prox(2,ipn)),	&
           delp_lo(k,prox(3,ipn)),delp_lo(k,prox(4,ipn)),	&
           delp_lo(k,prox(5,ipn)),delp_lo(k,prox(6,ipn)),	&
           delp(k,prox(1,ipn)),delp(k,prox(2,ipn)),		&
           delp(k,prox(3,ipn)),delp(k,prox(4,ipn)),		&
           delp(k,prox(5,ipn)),delp(k,prox(6,ipn))  )

     dpmin = min(delp_lo(k,ipn),delp(k,ipn),			&
           delp_lo(k,prox(1,ipn)),delp_lo(k,prox(2,ipn)),	&
           delp_lo(k,prox(3,ipn)),delp_lo(k,prox(4,ipn)),	&
           delp_lo(k,prox(5,ipn)),delp_lo(k,prox(6,ipn)),  	&
           delp(k,prox(1,ipn)),delp(k,prox(2,ipn)),		&
           delp(k,prox(3,ipn)),delp(k,prox(4,ipn)),		&
           delp(k,prox(5,ipn)),delp(k,prox(6,ipn))  )

     dpmax = max(0.,dpmax)	! cannot allow negative dpmax
     dpmin = max(0.,dpmin)	! cannot allow negative dpmin	

     ! q_plus/q_mnus are the upper/lower limits on antidiffusive dp tendencies
     ! q_plus is Zalesak (8):
     q_plus = (dpmax-delp_lo(k,ipn)				& ! N/sec
        -(adbash2*min(0.,dp_tndcy(k,ipn, of))			&
        + adbash3*max(0.,dp_tndcy(k,ipn,vof))))			&
         /adbash1*area(ipn)

     ! q_mnus is Zalesak (11):
     q_mnus = (delp_lo(k,ipn)-dpmin				& ! N/sec
        +(adbash2*max(0.,dp_tndcy(k,ipn, of))			&
        + adbash3*min(0.,dp_tndcy(k,ipn,vof))))			&
         /adbash1*area(ipn)

     !  Having p_plus(k,ipn) and q_plus, we can calc r_plus, Zal (9):

     !  reduce fluxes to stay within limits posed by q_plus,q_mnus.
     !  r_plus,r_mnus are dimensionless

     r_plus(k,ipn) = max(0.,min(1.,q_plus/max(p_plus(k,ipn),thshld)))
     r_mnus(k,ipn) = max(0.,min(1.,q_mnus/max(p_mnus(k,ipn),thshld)))

    end do			! loop through layers
   end do			! horizontal loop
!
!   Now we have chosen flux limiters that will assure monotonicity.
!   Next, we do the clipping.

!.......................................................
!   Sec. 4. Flux Clipping
!.......................................................

!   As explained by Durran, once you have the r_plus and
!   r_mnus over the whole grid, you can assure that the clipping
!   can be done so that it does not cause a problem in the center
!   cell, NOR IN ANY OF THE NEIGHBORING CELLS.  The clipping
!   is from Zalesek (13):

   call IncrementTimer(t1,tcnuity)

if (TimingBarriers) then
   call StartTimer(t1)
!SMS$BARRIER
   call IncrementTimer(t1,tcnuityBa)
endif

   call StartTimer(t1)
!SMS$EXCHANGE(r_plus,r_mnus)
!sms$compare_var(antifx    , "cnuity.F90 - antifx9 ")
   call IncrementTimer(t1,tcnuityEx)

   call StartTimer(t1)

!SMS$HALO_COMP(<1,1>) BEGIN
!DIR$ vector always
   do ipn = 1,nip		! horizontal loop

    psurf(ipn) = 0.
    do k = nvl,1,-1		! loop through layers (top-down for psurf)
     psurf(ipn) = psurf(ipn)+delp   (k,ipn)
!    psurf(ipn) = psurf(ipn)+delp_lo(k,ipn)
    end do

    do edgcount = 1,nedge(ipn)	! loop through edges
     edg = permedge(edgcount,ipn)
     recovr(edg,ipn) = 0.

     do k = 1,nvl		! loop through layers
      if (antifx(k,edg,ipn).ge.0.) then			! outgoing
       clip = min(r_mnus(k,ipn),r_plus(k,prox(edg,ipn)))
      else						! incoming
       clip = min(r_plus(k,ipn),r_mnus(k,prox(edg,ipn)))
      end if

      ! limit antidiffusive fluxes
      ! set aside vertical integral of discarded fluxes for later recovery
      recovr(edg,ipn) = recovr(edg,ipn)+antifx(k,edg,ipn)*(1.-clip)
      antifx(k,edg,ipn) = antifx(k,edg,ipn)*clip

     ! add clipped antidiff to low-order flux to obtain total mass flux

      massfx(k,edg,ipn,nf) = massfx(k,edg,ipn,nf)+antifx(k,edg,ipn)
     end do			! loop through layers
    end do			! loop through edges
   end do			! horizontal loop

!ss   ! include fluxes discarded during clipping process as barotropic
!ss   !  corrections to antidiffusive fluxes
!ss
!ss   do ipn = 1,nip			! horizontal loop
!ss    do edgcount = 1,nedge(ipn)	! loop through edges
!ss     edg = permedge(edgcount,ipn)
!ss     if (recovr(edg,ipn).ge.0.) then		! outgoing
!ss      do k = 1,nvl
!ss       antifx(k,edg,ipn) = antifx(k,edg,ipn)			&
!ss         +recovr(edg,ipn)*delp   (k,ipn )/psurf(ipn )
!ss!        +recovr(edg,ipn)*delp_lo(k,ipn )/psurf(ipn )
!ss      end do
!ss     else					! incoming
!ss      ipx = prox(edg,ipn)
!ss      do k = 1,nvl
!ss       antifx(k,edg,ipn) = antifx(k,edg,ipn)			&
!ss         +recovr(edg,ipn)*delp   (k,ipx)/psurf(ipx)
!ss!        +recovr(edg,ipn)*delp_lo(k,ipx)/psurf(ipx)
!ss      end do
!ss     end if			! recovr > or < 0
!ss    end do			! loop through edges
!ss   end do			! horizontal loop

!SMS$HALO_COMP END

!DIR$ vector always
   do ipn = 1,nip		! horizontal loop
    do k = 1,nvl		! loop through layers
     anti_tndcy(k,ipn) = 0.

     ! sum up edge fluxes to get antidiffusive tendency term
     do edg = 1,nprox(ipn)	! loop through edges
      anti_tndcy(k,ipn) = anti_tndcy(k,ipn)+antifx(k,edg,ipn)
     end do			! loop through edges

     ! divide antidiff tendency by cell area to convert to  Pa/sec
     anti_tndcy(k,ipn) = -anti_tndcy(k,ipn)*rarea(ipn)

     ! combine low-order with clipped antidiffusive tendency
     dp_tndcy(k,ipn,nf) = dplo_tndcy(k,ipn,nf)+anti_tndcy(k,ipn)

     ! advance delp to new time step

     delp(k,ipn) = delp(k,ipn)				&
          +adbash1*dp_tndcy(k,ipn, nf)			&
          +adbash2*dp_tndcy(k,ipn, of)			&
          +adbash3*dp_tndcy(k,ipn,vof)
    end do			! loop through layers
   end do			! horizontal loop

  end if			! low_ord = true or false

!! anti_tndcy(:,:)=dp_tndcy(:,:,nf)
!! do k=1,nvl,7
!!  write (string,'(a,i2)') 'dp_tndcy (nf), k=',k
!!  call findmxmn2(anti_tndcy,nvl,nip,k,string)
!!  write (string,'(a,i2)') 'cnuity: new dp k=',k
!!  call findmxmn2(delp,nvl,nip,k,string)
!! end do

  write (string,'(a,i6,2x)') 'step',its                        
  call stenedg(delp_lo,antifx,nvl,				&
               string(1:12)//'(cnuity) low-ord dp, antidiff flx')
  call stencl(delp,nvl,.01,string(1:12)//'(cnuity) new dp')

  do ipn = 1,nip		! horizontal loop
   coltend = 0.			! integrated mass flux convergence
   do k = nvl,1,-1		! loop through layers (top-down for p,omega)

    ! update pressure and Exner fct by vertically summing up thickness values
    pres(k,ipn) = pres(k+1,ipn)+delp(k,ipn)
    exner(k,ipn) = cp*(pres(k,ipn)/p1000)**(rd/cp)

    ! evaluate omega = dp/dt as
    ! (partial_p/partial_t) + (v_dot_grad_p) + (s_dot partial_p/ partial_s)

    lyrtend = (adbash1*dp_tndcy(k,ipn, nf)			&
              +adbash2*dp_tndcy(k,ipn, of)			&
              +adbash3*dp_tndcy(k,ipn,vof))/dt
    omega(k,ipn) = coltend+.5*lyrtend		! evaluate at mid depth
    coltend = coltend+lyrtend			! flux conv. intgral

    ! pressure gradient
    dpdx = 0.
    dpdy = 0.
    do edgcount = 1,nedge(ipn)	! loop through edges
     edg = permedge(edgcount,ipn)
     dpdx = dpdx+lp_edg(k,edg,ipn)*sidevec_c(2,edg,ipn)
     dpdy = dpdy-lp_edg(k,edg,ipn)*sidevec_c(1,edg,ipn)
    end do

    old = omega(k,ipn)
    omega(k,ipn) = omega(k,ipn)					&
                  +(u_velo(k,ipn)*dpdx				&
                   +v_velo(k,ipn)*dpdy)*rarea(ipn)

    if (ipn.eq.PrintIpnDiag .and. mod(k,7).eq.1) then
!SMS$IGNORE BEGIN
     print '(a,i8,i4,a,3f9.3)','ipn,k  = ',ipn,k,			&
      '  omega terms 1+3,term 2,total:',old,(u_velo(k,ipn)*dpdx		&
         + v_velo(k,ipn)*dpdy)*rarea(ipn),omega(k,ipn)
!SMS$IGNORE END
    end if

   end do			! loop through layers
  end do			! horizontal loop

! do k = 1,nvl,7
!  write (string,'(a,i3,a)') 'k',k,' cnuity:omega'
!  call findmxmn2(omega,nvl,nip,k,string)
! end do
! print *

  call IncrementTimer(t1,tcnuity)

if (.not.low_ord) then
 !sms$compare_var(r_plus   , "cnuity.F90 - r_plus10 ")
 !sms$compare_var(r_mnus   , "cnuity.F90 - r_mnus10 ")
 !sms$compare_var(antifx   , "cnuity.F90 - antifx10 ")
end if
!sms$compare_var(massfx    , "cnuity.F90 - massfx10")
!sms$compare_var(dplo_tndcy, "cnuity.F90 - dplo_tndcy10")
!sms$compare_var(delp_lo   , "cnuity.F90 -  delp_lo10")

!SMS$PARALLEL END

  return
  end subroutine cnuity
  end module module_cnuity
