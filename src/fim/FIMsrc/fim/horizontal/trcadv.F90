module  module_trcadv
use stencilprint
use findmaxmin2
use findmaxmin1
contains
!*********************************************************************
!     trcadv
! 	trcadv = flux corrected transport for mass field tracers
!	J. Lee		Author				September 2005
!	A. E. MacDonald Documentor			November  2005
!	R. Bleck	major rewrite			April     2008
!       R. Bleck        fixed bug in high-order flux	April     2011
!
! 	This routine is based on Zalesak, JOURNAL OF COMPUTATIONAL
!       PHYSICS, 31, 335-362, 1979.  Dale Durran provides an
!       excellent discussion of flux corrected transport in his book
!     NUMERICAL METHODS FOR WAVE EQUATIONS IN GEOPHYSICAL FLUID DYNAMICS.
!*********************************************************************

  subroutine trcadv (its,nf,of,vof,adbash1,adbash2,adbash3,	&
    trcr_edg,tracr,trcdp,trc_tdcy,trclo_tdcy,massfx,delp,       &
    ttrcadv,ttrcadvEx,ttrcadvBa,TimingBarriers )

use module_control  ,only: nvl,npp,nip,nabl,dt,ntra,nd,PrintIpnDiag
use module_constants,only: nprox,prox,proxs,rarea,area,		&
                           nedge,permedge
implicit none

!..............................................................
!	Sec. 0  Dimension and Type
!..............................................................

! External variables:
  integer,intent (IN)    :: its		! model time step
  integer,intent (IN)    :: nf,of,vof	! time slots: new,old,very old
  real   ,intent (IN)    :: adbash1,adbash2,adbash3
  real*8 ,intent (INOUT) :: ttrcadv          ! computation timer
  real*8 ,intent (INOUT) :: ttrcadvEx        ! halo update communication timer
  real*8 ,intent (INOUT) :: ttrcadvBa        ! barrier timer for task skew
  logical,intent (IN)    :: TimingBarriers   ! measure task skew when .true.
!SMS$DISTRIBUTE (dh,nip) BEGIN
  real   ,intent (IN)    :: trcr_edg   (nvl,npp,nip,ntra)
  real   ,intent (INOUT) :: tracr      (nvl,nip,ntra)
  real   ,intent (INOUT) :: trcdp      (nvl,nip,ntra)
  real   ,intent (INOUT) :: trc_tdcy   (nvl,nip,nabl,ntra)
  real   ,intent (INOUT) :: trclo_tdcy (nvl,nip,nabl,ntra)
  real   ,intent (IN)    :: massfx     (nvl,npp,nip,3)
  real   ,intent (IN)    :: delp       (nvl,nip)
! Local variables:
  real  :: s_plus(nvl,nip)	! Zalesak's p_plus
  real  :: s_mnus(nvl,nip)	! Zalesak's p_minus
  real  :: r_plus(nvl,nip)	! Zalesak's r_plus
  real  :: r_mnus(nvl,nip)	! Zalesak's r_minus 
  real  :: trcr_lo(nvl,nip)	! tracer after low-order transport
  real  :: antiflx(nvl,npp,nip) ! antidiffusive tracer flux
  real  :: trcdp_lo(nvl,nip)	! tracer*dp after low-order transport
  real  :: trmax(nvl,nip)	! regional tracer max for flux clipping
  real  :: trmin(nvl,nip)	! regional tracer min for flux clipping
  real  :: anti_tdcy(nvl,nip,ntra)! antidiff trcr tendency
  real  :: flxlo(nvl,npp,nip)   ! tracer flux, low order
  real  :: q_plus(nip)		! Zalesak's q_plus
  real  :: q_mnus(nip)		! Zalesak's q_minus
!SMS$DISTRIBUTE END

  integer           :: k	! layer index
  integer           :: ipn	! icos point number index
  integer           :: edg	! icos edge number index
  integer           :: type	! tracer index (1=theta; 2=specif.hum., ...)
  real              :: del_dp	! delta_p used for upstream integral of flux
  real              :: flxhi	! tracer flux, high order
  integer           :: ipx	! neighbor across joint edge
  integer           :: edx	! neighbor's index of joint edge
  integer           :: edgcount ! count of icos edges
  character         :: string*32
  real,parameter    :: thshld = 1.e-11
  logical,parameter :: low_ord = .false. ! if true, skip antidiffusive part
  real*8  :: t1

!.............................................................
!  Sec. 1 Calculate Low and Antidiffusive Flux
!.............................................................

! Calculates the FCT low and high order fluxes, and defines
! the antidiffusive flux as the difference between the high and
! low order fluxes.  The low order flux is computed based on 
! the assumption of piecewise continuity, with a constant value 
! in each cell.  The higher order uses the "gazebo" with a
! sloped line used for the flux integral.

! Avoid exchange via HALO_COMP in cnuity.F90 and edgvar.F90
!!!SMS$EXCHANGE(massfx,trcr_edg)
!!SMS$EXCHANGE(tracr) exchanged in edgvar

!sms$compare_var(tracr  , "trcadv.F90 - tracr1   ")
!sms$compare_var(massfx , "trcadv.F90 - massfx1   ")
!sms$compare_var(trc_tdcy, "trcadv.F90 - trc_tdcy1  ")
!sms$compare_var(trcr_edg, "trcadv.F90 - trcr_edg1  ")

  call StartTimer(t1)

!!  do k=1,nvl,7
!!   write (string,'(a,i3,a)') 'k',k,' trcadv:theta'
!!   call findmxmn3(tracr,nvl,nip,ntra,k,1,string)
!!  end do
!!  print *

  write (string,'(a,i7)') '(atm trcadv) step',its
  call stencl(tracr,nvl,1.,trim(string)//', old theta')

  trclo_tdcy(:,:,nf,:)=0.
  trc_tdcy  (:,:,nf,:)=0.
  anti_tdcy (:,:,   :)=0.

  ! Initialize these local arrays so COMPARE_VAR does not get 
  ! confused by uninitialized edg=6 edges of pentagonal grid cells.  This 
  ! code could be safely omitted if COMPARE_VAR were not used.  
  antiflx(:,npp,:) = 0.

  do type=1,ntra			! loop through tracers

!SMS$PARALLEL (dh,ipn) BEGIN
!SMS$HALO_COMP(<1,1>) BEGIN
   do ipn=1,nip				! horizontal loop
    do edgcount=1,nedge(ipn)		! loop through edges
     edg = permedge(edgcount,ipn)
     ipx=prox(edg,ipn)			! neighbor across edge
     edx=proxs(edg,ipn)		! index of joint edge as seen by neighbor

     do k=1,nvl				! loop through layers

! --- high-order tracer flux (2nd order centered, out of cell > 0):
      flxhi=massfx(k,edg,ipn,nf)*trcr_edg(k,edg,ipn,type)	! trcdim x N/s

! --- low-order tracer flux (donor-cell, out of cell > 0):
      flxlo(k,edg,ipn)=0.5*(						&
           (massfx(k,edg,ipn,nf)+abs(massfx(k,edg,ipn,nf)))		&
           *tracr (k,ipn,type)						&
         - (massfx(k,edx,ipx,nf)+abs(massfx(k,edx,ipx,nf)))		&
           *tracr (k,ipx,type) )
    
      ! get anti-diffusive flux as the difference between high order
      ! (flxhi) and low order flux (flxlo), from Zalesak, p336, Eqn (3):

      antiflx(k,edg,ipn)=flxhi-flxlo(k,edg,ipn)   ! trcdim x N/sec

     end do			! loop through layers
    end do			! loop through edges
   end do			! horizontal loop
!SMS$HALO_COMP END

   do ipn=1,nip			! horizontal loop
    do edg=1,nprox(ipn)		! loop through edges
     do k=1,nvl			! loop through layers

      ! sum up edge fluxes to get low-order tracer tendency
      trclo_tdcy(k,ipn,nf,type)=trclo_tdcy(k,ipn,nf,type)+flxlo(k,edg,ipn)

     end do			! loop through layers
    end do			! loop through edges

  !.....................................................................
  !  Sec. 2.  Calculate new low order trcr*dp using full Adams Bashforth
  !.....................................................................

    do k=1,nvl			! loop through layers

  !  divide tendency by cell area to convert to (trcdim x Pa/sec)
     trclo_tdcy(k,ipn,nf,type)=-trclo_tdcy(k,ipn,nf,type)*rarea(ipn)

     ! get new value for the low order tracer*dp field using Adams Bashforth
     ! 3 time levels, the one just calculated (nf), and the two prev ones,
     ! marked of (old field) and vof (very old field):

     trcdp_lo(k,ipn) = trcdp(k,ipn,type)				&
             +adbash1*trclo_tdcy(k,ipn, nf,type)			&
             +adbash2*trclo_tdcy(k,ipn, of,type)			&
             +adbash3*trclo_tdcy(k,ipn,vof,type)

     ! set upper/lower bounds for ratio of trcdp_lo and delp
     trmax(k,ipn)=max(tracr(k,ipn,type),				&
       tracr(k,prox(1,ipn),type),tracr(k,prox(2,ipn),type),		&
       tracr(k,prox(3,ipn),type),tracr(k,prox(4,ipn),type),		&
       tracr(k,prox(5,ipn),type),tracr(k,prox(6,ipn),type)  )

     trmin(k,ipn)=min(tracr(k,ipn,type),				&
       tracr(k,prox(1,ipn),type),tracr(k,prox(2,ipn),type),		&
       tracr(k,prox(3,ipn),type),tracr(k,prox(4,ipn),type),		&
       tracr(k,prox(5,ipn),type),tracr(k,prox(6,ipn),type)  )

     trmax(k,ipn)=max(0.,trmax(k,ipn))	! cannot allow negative trmax
     trmin(k,ipn)=max(0.,trmin(k,ipn))	! cannot allow negative trmin	

     ! now divide trcdp_lo by delp to get new low order tracer field 
     trcr_lo(k,ipn)=max(trmin(k,ipn),min(trmax(k,ipn),			&
           trcdp_lo(k,ipn)/max(thshld,delp(k,ipn)) ))

    end do			! loop through layers
   end do			! horizontal loop

   if (low_ord) then

    tracr(:,:,type) = trcr_lo(:,:)

   else				! evaluate and apply antidiffusive fluxes

  ! Dale Durran's book indicates that condition Zal (14) should be
  ! satisfied, although Zalesak says it's cosmetic.  We believe Dale:
  ! Also Calculate s_plus Zalesak (7) and s_minus Zalesak (10).
  ! (aggregate of incoming and outgoing fluxes, trcdim x N/sec)

   s_plus=0.
   s_mnus=0.

   call IncrementTimer(t1,ttrcadv)

   if (TimingBarriers) then
    call StartTimer(t1)
!SMS$BARRIER
    call IncrementTimer(t1,ttrcadvBa)
   endif

   call StartTimer(t1)
!SMS$EXCHANGE(trcr_lo)
   call IncrementTimer(t1,ttrcadvEx)

   call StartTimer(t1)

   do ipn=1,nip			! horizontal loop
    do edg=1,nprox(ipn)		! loop through edges
     do k=1,nvl			! loop through layers
      if(antiflx(k,edg,ipn).le.0.) then			! flux into cell
       s_plus(k,ipn)=s_plus(k,ipn)-antiflx(k,edg,ipn)
      else						! flux out-of cell
       s_mnus(k,ipn)=s_mnus(k,ipn)+antiflx(k,edg,ipn)
      endif
     end do			! loop through layers
    end do			! loop through edges

  !............................................................
  !  Sec 3.  Monotonicity Limit on Fluxes
  !............................................................

  !  Determine the amount of antidiffusive fluxes that can be
  !  added to the low order solution without violating monotonicity.

    do k=1,nvl			! loop through layers

     ! According to Zal (17), we must limit according to max of
     ! tracer from any gazebo direction, in either current or prev time
     ! step.  Likewise for the minimum.
     ! For the pentagons prox(6,ipn) is set to prox(5,ipn) in init.F90.

     ! relax upper/lower bounds by incorporating new low-order field
     trmax(k,ipn)=max(trmax(k,ipn),trcr_lo(k,ipn),		&
       trcr_lo(k,prox(1,ipn)),trcr_lo(k,prox(2,ipn)),		&
       trcr_lo(k,prox(3,ipn)),trcr_lo(k,prox(4,ipn)),		&
       trcr_lo(k,prox(5,ipn)),trcr_lo(k,prox(6,ipn))  )

     trmin(k,ipn)=min(trmin(k,ipn),trcr_lo(k,ipn),		&
       trcr_lo(k,prox(1,ipn)),trcr_lo(k,prox(2,ipn)),		&
       trcr_lo(k,prox(3,ipn)),trcr_lo(k,prox(4,ipn)),		&
       trcr_lo(k,prox(5,ipn)),trcr_lo(k,prox(6,ipn))  )

     trmax(k,ipn)=max(0.,trmax(k,ipn))	! cannot allow negative trmax
     trmin(k,ipn)=max(0.,trmin(k,ipn))	! cannot allow negative trmin	

     ! q_plus/q_mnus are the upper/lower limits on antidiff trcr tendencies,
     ! dimensioned trdim x N/sec

     ! q_plus is Zalesak (8):
     q_plus(ipn) = ((trmax(k,ipn)-trcr_lo(k,ipn))*delp(k,ipn)	&
            -(adbash2*min(0.,trc_tdcy(k,ipn, of,type))		&
            + adbash3*max(0.,trc_tdcy(k,ipn,vof,type))))	&
             /adbash1*(area(ipn))

     ! q_mnus is Zalesak (11):
     q_mnus(ipn) = ((trcr_lo(k,ipn)-trmin(k,ipn))*delp(k,ipn)	&
            +(adbash2*max(0.,trc_tdcy(k,ipn, of,type))		&
            + adbash3*min(0.,trc_tdcy(k,ipn,vof,type))))	&
             /adbash1*(area(ipn))

     ! Having s_plus(k,ipn) and q_plus, we can calc r_plus, Zal (9):

     !  reduce fluxes to stay within limits posed by q_plus/q_mnus
     !  r_plus/r_mnus are dimensionless

     r_plus(k,ipn)=max(0.,min(1.,q_plus(ipn)/max(s_plus(k,ipn),thshld)))
     r_mnus(k,ipn)=max(0.,min(1.,q_mnus(ipn)/max(s_mnus(k,ipn),thshld)))

!!     if (ipn.eq.PrintIpnDiag) then
!!!SMS$IGNORE BEGIN
!!      print '(2i7,2i3,a/7es11.3)',its,ipn,k,type,' q_plus etc:',	&
!!            q_plus(ipn),trmax(k,ipn),trcr_lo(k,ipn),delp(k,ipn),	&
!!            - adbash2*min(0.,trc_tdcy(k,ipn, of,type)),		&
!!            - adbash3*max(0.,trc_tdcy(k,ipn,vof,type)),		&
!!             adbash1
!!      print '(2i7,2i3,a/7es11.3)',its,ipn,k,type,' q_mnus etc:',	&
!!            q_mnus(ipn),trcr_lo(k,ipn),trmin(k,ipn),delp(k,ipn),	&
!!            + adbash2*max(0.,trc_tdcy(k,ipn, of,type)),		&
!!            + adbash3*min(0.,trc_tdcy(k,ipn,vof,type)),		&
!!             adbash1
!!!SMS$IGNORE END
!!     end if

    end do			! loop through layers
   end do			! horizontal loop

!  do k=1,nvl,7
!   write (string,'(a3,i3)') ' k=',k
!   call findmxmn2(r_plus,nvl,nip,k,'(trcadv) r_plus'//string(1:6))
!   call findmxmn2(r_mnus,nvl,nip,k,'(trcadv) r_mnus'//string(1:6))
!  end do
!  call stencl(r_plus,nvl,1000.,'(trcadv) r_plus x 1000')
!  call stencl(r_mnus,nvl,1000.,'(trcadv) r_mnus x 1000')

  ! Now we have a flux limiter that will assure monotonicity.
  ! Next, we do the clipping.

  !.......................................................
  !  Sec. 4. Flux Clipping
  !.......................................................

  ! As explained by Durran, once you have the r_plus and
  ! r_mnus over the whole grid, you can assure that the clipping
  ! can be done so that it does not cause a problem in the center
  ! cell, NOR IN ANY OF THE NEIGHBORING CELLS.  The clipping
  ! is from Zalesak (13):

   call IncrementTimer(t1,ttrcadv)

   if (TimingBarriers) then
    call StartTimer(t1)
!SMS$BARRIER
    call IncrementTimer(t1,ttrcadvBa)
   endif

   call StartTimer(t1)
!SMS$EXCHANGE(r_plus,r_mnus)
   call IncrementTimer(t1,ttrcadvEx)

!sms$compare_var(s_plus   , "trcadv.F90 - s_plus4 ")
!sms$compare_var(s_mnus   , "trcadv.F90 - s_mnus4 ")
!sms$compare_var(r_plus   , "trcadv.F90 - r_plus4 ")
!sms$compare_var(r_mnus   , "trcadv.F90 - r_mnus4 ")
!sms$compare_var(antiflx  , "trcadv.F90 - antiflx4 ")
!sms$compare_var(trcr_lo  , "trcadv.F90 - trcr_lo3 ")

   call StartTimer(t1)

   do ipn=1,nip			! horizontal loop
    do edg=1,nprox(ipn)		! loop through edges
     do k=1,nvl			! loop through layers

      if (antiflx(k,edg,ipn).ge.0.) then		! outgoing
       antiflx(k,edg,ipn) = antiflx(k,edg,ipn)			&
         * min(r_mnus(k,ipn),r_plus(k,prox(edg,ipn)))
      else						! incoming
       antiflx(k,edg,ipn) = antiflx(k,edg,ipn)			&
         * min(r_plus(k,ipn),r_mnus(k,prox(edg,ipn)))
      end if

      ! sum up edge fluxes to get antidiff tracer tendency
      anti_tdcy(k,ipn,type)=anti_tdcy(k,ipn,type)+antiflx(k,edg,ipn)

     end do			! loop through layers
    end do			! loop through edges

    do k=1,nvl			! loop through layers

     ! divide tendency by cell area to convert to  trcdim x Pa/sec
     anti_tdcy(k,ipn,  type)=-anti_tdcy(k,ipn,  type)*rarea(ipn)

     ! combine low order with clipped antidiff tendency
     trc_tdcy(k,ipn,nf,type)=trclo_tdcy(k,ipn,nf,type)		&
                            + anti_tdcy(k,ipn,   type)

     ! advance (tracer x thickness) field in time
     trcdp(k,ipn,type) = trcdp(k,ipn,type)			&
           +adbash1*trc_tdcy(k,ipn, nf,type)			&
           +adbash2*trc_tdcy(k,ipn, of,type)			&
           +adbash3*trc_tdcy(k,ipn,vof,type)

     ! finally divide new trcr*dp by new delp to get new tracer field 
     tracr(k,ipn,type)=max(trmin(k,ipn),min(trmax(k,ipn),	&
           trcdp(k,ipn,type)/max(thshld,delp(k,ipn)) ))

    end do			! loop through layers
   end do			! horizontal loop

   end if			! low_ord = true or false
!SMS$PARALLEL END
  end do                	! loop through tracers

  write (string,'(a,i7)') '(atm trcadv) step',its
  call stencl(tracr,nvl,1.,trim(string)//', new theta')

  call IncrementTimer(t1,ttrcadv)

!sms$compare_var(antiflx   , "trcadv.F90 - antiflx5 ")
!sms$compare_var(trcdp_lo  , "trcadv.F90 - trcdp_lo5")
!sms$compare_var(trcdp     , "trcadv.F90 - trcdp5 ")
!sms$compare_var(trclo_tdcy, "trcadv.F90 - trclo_tdcy5  ")

  return
end subroutine trcadv
end module module_trcadv
