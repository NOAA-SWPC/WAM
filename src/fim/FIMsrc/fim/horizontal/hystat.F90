module module_hystat
use findmaxmin1
contains
!*********************************************************************
!     hystat
!	Hydrostatic equation
!	Alexander E. MacDonald  11/14/2005
!	J.Lee                   01/04/2006
!	R.Bleck                 12/06/2007
!*********************************************************************

   subroutine hystat (its,   &
     ph3d,                   &	! geopotential (=g*z)
     ex3d,mp3d,dp3d,         &	! ,exner fct, mont pot, layer thickness
     tr3d,trdp,psrf,ptdcy     )	! tracer, tracer x thickness, srf.prs tndcy

   use module_control  ,only: nvl,nvlp1,nip,ptop,ntra,PrintIpnDiag
   use module_constants,only: p1000,cp,rd,qvmin,qwmin
   implicit none

! Dimension and type external variables:
   integer,intent (IN) :: its			! model time step
!SMS$DISTRIBUTE (dh,nip) BEGIN
   real,intent (INOUT) :: ph3d(nvlp1,nip)	! geopotential
   real,intent (IN)    :: ex3d(nvlp1,nip)	! exner fct
   real,intent (OUT)   :: mp3d(nvl,nip)		! montgomery potential
   real,intent (IN)    :: dp3d(nvl,nip)		! layer thickness
   real,intent (INOUT) :: tr3d(nvl,nip,ntra)	! mass field tracers
   real,intent (OUT)   :: trdp(nvl,nip,ntra)	! tracer x thickness
   real,intent (INOUT) :: psrf(nip)		! surface pressure
!JR Cannot have intent(out) for ptdcy because this routine only sets half of the array, and Lahey
!JR assigns a "bad sequence of bits" to all intent(out) variables.
   real,intent (INOUT) :: ptdcy(nip,2)		! srf.pres.tdcy, 2 time lvls

! Declare local variables:
   real       :: work(nip)
!SMS$DISTRIBUTE END
   integer    :: ipn	! icos index
   integer    :: k	! layer index
   integer    :: ns	! tracer index
   logical    :: vrbos	! switch for 'verbose' mode

!  Note that tracers, velocity, and montgomery potential, mp3d
!  are all constant through the layers.  Phi (ph3d) and pressure (dp3d)
!  vary through the layer.

!  Layer variables:  tr3d,mp3d
!  Level variables:  ex3d,ph3d

!SMS$PARALLEL (dh,ipn) BEGIN
!sms$compare_var(ex3d, "diag.F90 - ex3d5 ")
!sms$compare_var(ph3d, "diag.F90 - ph3d5 ")

   do ipn=1,nip		!  global icos loop
    vrbos=ipn.eq.PrintIpnDiag

! --- srf.prs tendency is needed to evaluate model noise (diagnoise.F90)
     work(ipn)=psrf(ipn)
     psrf(ipn)=dp3d(nvl,ipn)
     do k=nvl-1,1,-1
      psrf(ipn)=psrf(ipn)+dp3d(k,ipn)
     end do
     ns=mod(its,2)+1
     if (its.gt.0) ptdcy(ipn,ns)=psrf(ipn)-work(ipn)

!.........................................................
!   Determine bottom layer values
!.........................................................

     mp3d(1,ipn)=ex3d(1,ipn)*tr3d(1,ipn,1)+ph3d(1,ipn)	! mont pot, layer 1
     ph3d(2,ipn)=mp3d(1,ipn)-tr3d(1,ipn,1)*ex3d(2,ipn)	! geopot, level 2

     do k=2,nvl		!  vertical loop

! Hydrostatic eqn:  d mp/d theta = exner fct, tr3d(.,.,1) = theta
       mp3d(k,ipn)=mp3d(k-1,ipn)+ex3d(k,ipn)			&
                   *(tr3d(k,ipn,1)-tr3d(k-1,ipn,1))

! get geopotential from identity  montg pot = geopot + theta*exner
       ph3d(k+1,ipn)=mp3d(k,ipn)-tr3d(k,ipn,1)*ex3d(k+1,ipn)

     enddo		! vertical loop

! keep vapor mixing ratio and water content above prescribed lower limits
     do k=1,nvl
       tr3d(k,ipn,2)=max(qvmin,tr3d(k,ipn,2))
       tr3d(k,ipn,3)=max(qwmin,tr3d(k,ipn,3))
     end do

! compute tracer amount per unit area (tracer concentration x layer thickness)
     do ns=1,ntra
       do k=1,nvl
         trdp(k,ipn,ns)=tr3d(k,ipn,ns)*dp3d(k,ipn)
       end do
     end do

     if (vrbos) then
!SMS$IGNORE BEGIN
      write (6,100) its,ipn,(k,1000.*(ex3d(k,ipn)/cp)**(cp/rd),	&
       ex3d(k,ipn),ph3d(k,ipn),mp3d(k,ipn),k=1,nvl),			&
        nvlp1,1000.*(ex3d(nvlp1,ipn)/cp)**(cp/rd),ex3d(nvlp1,ipn),	&
         ph3d(nvlp1,ipn)
 100 format ('its,ipn=',i6,i8,                                          &
      '  HYSTAT    pres   exn.fct    geopot    montg'/(i28,4f10.1))
!SMS$IGNORE END
     end if

   enddo		! horizontal loop 
!SMS$PARALLEL END

!! call findmxmn1(work,nip,'old srf.pres.')
!! call findmxmn1(psrf,nip,'new srf.pres.')
!! ns=mod(its+1,2)+1
!! work(:)=ptdcy(:,ns)
!! call findmxmn1(work,nip,'old srf.pres.tdcy')
!! ns=mod(its  ,2)+1
!! work(:)=ptdcy(:,ns)
!! call findmxmn1(work,nip,'new srf.pres.tdcy')

   return
   end subroutine hystat
end module module_hystat
