module module_edgvar
use stenedgprint
contains
!*********************************************************************
!	edgvar
!	Interpolates data to cell edges in the local stereographic grid
!	J. Lee                    Sep 2005
!	A. E. MacDonald		  Nov 2005  fim conversion
!	R. Bleck                  Apr 2008  cosmetic changes
!*********************************************************************
  
  subroutine edgvar (its,		&
    u_vel,v_vel,			& ! west, south wind on s level
    delp,pres,				& ! thickness, intfc pressure
    montg,tracr,			& ! montg.pot., mass field tracers
    u_edg,v_edg,dp_edg,lp_edg,		& ! u v dp and layer-prs edge variables
    bnll_edg,trc_edg,			& ! bernoulli fct, trcr edge variables
    tedgvar,tedgvarEx,tedgvarBa,	& ! timers
    TimingBarriers )
  
  use module_control  ,only: nvl,nvlp1,npp,nip,ntra,PrintIpnDiag
  use module_constants,only: nprox,proxs,prox,cs,sn,nedge,permedge,actual
  
  implicit none
  
!   Type and dimension external variables:
  
  integer,intent(IN)    :: its				! model time step
!SMS$DISTRIBUTE (dh,nip) BEGIN
!TBH:  Note that u_vel, v_vel, delp, pres, montg, and tracr are INTENT(IN) 
!TBH:  *except* for the SMS EXCHANGE directive which modifies their halos.  
!TBH:  Lahey is smart enough to catch this so we must specify INTENT(INOUT) 
!TBH:  for these dummy arguments.  
  real   ,intent(INOUT) :: u_vel (nvl,nip)		! west wind
  real   ,intent(INOUT) :: v_vel (nvl,nip)		! south wind
  real   ,intent(INOUT) :: delp  (nvl,nip)		! layer thickness
  real   ,intent(INOUT) :: pres  (nvlp1,nip)		! pres on interfaces
  real   ,intent(INOUT) :: montg (nvl,nip)		! montgomery potential
  real   ,intent(INOUT) :: tracr (nvl,nip,ntra)		! mass field tracers
  real   ,intent(OUT)   :: u_edg   (nvl,npp,nip)	! west wind on edges
  real   ,intent(OUT)   :: v_edg   (nvl,npp,nip)	! south wind on edges
  real   ,intent(OUT)   :: dp_edg  (nvl,npp,nip)	! layer thknss on edges
  real   ,intent(OUT)   :: lp_edg  (nvl,npp,nip)	! midlyr prs on edges
  real   ,intent(OUT)   :: trc_edg (nvl,npp,nip,ntra)	! tracers on edges
  real   ,intent(OUT)   :: bnll_edg(nvl,npp,nip)	! bernoulli fct on edges

  real work1d(nip),work2d(npp,nip)
!SMS$DISTRIBUTE END
  real*8 ,intent(INOUT) :: tedgvar             ! computation timer
  real*8 ,intent(INOUT) :: tedgvarEx           ! halo update communication timer
  real*8 ,intent(INOUT) :: tedgvarBa           ! barrier timer for task skew
  logical,intent(IN)    :: TimingBarriers      ! measure task skew when .true.
  
!   Type and dimension local variables:
  
  integer :: edg                   !   icos edge index
  integer :: ipx                   !   neighbor across edge with index 'edg'
  integer :: im1,ip1               !   neighbors across edg-1 and edg+1
  integer :: ipn                   !   icos index
  integer :: k 		           !   layer Index
  integer :: type                  !   tracer index: 1=theta, 2=qv, ....
  integer :: edgcount              !   count of icos edges
!    These are u and v at neighboring icos points, NOT on edge
  real    :: u_xy1,u_xy2,u_xy3,u_xy4 ! u on the xy local grid (m/s), at prox pt
  real    :: v_xy1,v_xy2,v_xy3,v_xy4 ! v on the xy local grid (m/s), at prox pt
  real    :: mpe  	           !   montg.pot. on edges
  real    :: kee                   !   kinetic energy on edges
  real*8  :: t1
  logical :: vrbos
  character      :: string*32
! real,parameter :: divby6=1./6.   !   use in trapezoidal rule
  real,parameter :: divby18=1./18. !   use in simpson's rule
  
  if (TimingBarriers) then
   call StartTimer(t1)
!SMS$BARRIER
   call IncrementTimer(t1,tedgvarBa)
  endif
  
!SMS$PARALLEL (dh,ipn) BEGIN
  call StartTimer(t1)
!SMS$EXCHANGE(u_vel,v_vel,delp,tracr,montg,pres)
  call IncrementTimer(t1,tedgvarEx)
! sms$compare_var(u_vel, "load_ls.F90 - u_vel5 ")
! sms$compare_var(v_vel, "load_ls.F90 - v_vel5 ")
! sms$compare_var(montg, "load_ls.F90 - montg5 ")
! sms$compare_var(cs  , "load_ls.F90 - cs5 ")
! sms$compare_var(sn  , "load_ls.F90 - sn5 ")
  
  call StartTimer(t1)
!SMS$HALO_COMP(<1,1>) BEGIN
  do ipn=1,nip				! horizontal loop
  vrbos=ipn.eq.PrintIpnDiag
   do edgcount=1,nedge(ipn)		! loop through edges
  
! --- edge quantities are interpolated from 4 icos cells -- the cells on
! --- either side of the edge plus the 2 immediate neighbors of this pair.
  
    edg = permedge(edgcount,ipn)
    ipx = prox(edg,ipn)
    im1=mod(edg-2+nprox(ipn),nprox(ipn))+1
    ip1=mod(edg             ,nprox(ipn))+1
    im1 = prox(im1,ipn)
    ip1 = prox(ip1,ipn)
  
    if (vrbos) then
!SMS$IGNORE BEGIN
     print 100,ipn,edg,ipn,actual(ipx),actual(im1),actual(ip1)
 100 format (i8,' (edgvar) edge',i2,' intpol based on icos cells'	&
      2i8,' (lrg wgt),'/51x,2i8,' (sml wgt)')
!SMS$IGNORE END
    end if
  
!DIR$ vector always
    do k=1,nvl				! Loop over layers
  
  !    Transform u,v at neighboring icos pt to local coord.system.
  !    cs and sn are coordinate transformation constants.
  !    u_xy,v_xy are values of u and v rotated into local system.
  !    (Unrolled to allow vectorization of the k loop)
  
     u_xy1 =  cs(1,edg,ipn)*u_vel(k,ipn)+sn(1,edg,ipn)*v_vel(k,ipn)
     v_xy1 = -sn(1,edg,ipn)*u_vel(k,ipn)+cs(1,edg,ipn)*v_vel(k,ipn)
     u_xy2 =  cs(2,edg,ipn)*u_vel(k,ipx)+sn(2,edg,ipn)*v_vel(k,ipx)
     v_xy2 = -sn(2,edg,ipn)*u_vel(k,ipx)+cs(2,edg,ipn)*v_vel(k,ipx)
     u_xy3 =  cs(3,edg,ipn)*u_vel(k,im1)+sn(3,edg,ipn)*v_vel(k,im1)
     v_xy3 = -sn(3,edg,ipn)*u_vel(k,im1)+cs(3,edg,ipn)*v_vel(k,im1)
     u_xy4 =  cs(4,edg,ipn)*u_vel(k,ip1)+sn(4,edg,ipn)*v_vel(k,ip1)
     v_xy4 = -sn(4,edg,ipn)*u_vel(k,ip1)+cs(4,edg,ipn)*v_vel(k,ip1)
  
  !   interpolate rotated wind components to edges
  
!    u_edg(k,edg,ipn) = (2.*(u_xy1+u_xy2)+u_xy3+u_xy4)*divby6
!    v_edg(k,edg,ipn) = (2.*(v_xy1+v_xy2)+v_xy3+v_xy4)*divby6
     u_edg(k,edg,ipn) = (8.*(u_xy1+u_xy2)+u_xy3+u_xy4)*divby18
     v_edg(k,edg,ipn) = (8.*(v_xy1+v_xy2)+v_xy3+v_xy4)*divby18
  
    end do			!  loop over layers

  !   interpolate layer thickness to edges
  
    do k=1,nvl	! Loop over layers
!    dp_edg (k,edg,ipn) = (2.*(delp(k,ipn)+delp(k,ipx))		&
!                       +      delp(k,im1)+delp(k,ip1))*divby6
     dp_edg (k,edg,ipn) = (8.*(delp(k,ipn)+delp(k,ipx))		&
                        +      delp(k,im1)+delp(k,ip1))*divby18
    end do			!  loop over layers
  
  !   interpolate mid-layer pressure to edges
  
    do k=1,nvl			! Loop over layers
     lp_edg (k,edg,ipn) =					&
!       (   (pres(k  ,ipn)+pres(k  ,ipx)			&
!           +pres(k+1,ipn)+pres(k+1,ipx))			&
!     + .5* (pres(k  ,im1)+pres(k  ,ip1)			&
!           +pres(k+1,im1)+pres(k+1,ip1)) )*divby6
        (4.*(pres(k  ,ipn)+pres(k  ,ipx)			&
           + pres(k+1,ipn)+pres(k+1,ipx))			&
      + .5* (pres(k  ,im1)+pres(k  ,ip1) 			&
            +pres(k+1,im1)+pres(k+1,ip1)) )*divby18
    end do			!  loop over layers
  
    do k=1,nvl			! Loop over layers
  
  !   interpolate montg.pot. to edges
  
!    mpe = (2.*(montg(k,ipn)+montg(k,ipx))			&
!             + montg(k,im1)+montg(k,ip1))*divby6
     mpe = (8.*(montg(k,ipn)+montg(k,ipx))			&
              + montg(k,im1)+montg(k,ip1))*divby18
  
  !   interpolate kinetic energy to edges
  
     kee = .5*(u_edg(k,edg,ipn)**2+v_edg(k,edg,ipn)**2)
  
  !   bernoulli function = mont.pot. (mp) + kinetic energy (ke):
     bnll_edg(k,edg,ipn) = mpe + kee
    end do			! loop over layers
   end do			! loop over edges
  end do			! horizontal loop
!SMS$HALO_COMP END
  
  !   interpolate tracer to edges
  
  do type=1,ntra			! loop through tracers
!SMS$HALO_COMP(<1,1>) BEGIN
   do ipn=1,nip			! horizontal loop
    trc_edg(:,npp,ipn,type)=0.
    do edgcount=1,nedge(ipn)	! loop through edges
  
     edg = permedge(edgcount,ipn)
     ipx = prox(edg,ipn)
     im1=mod(edg-2+nprox(ipn),nprox(ipn))+1
     ip1=mod(edg             ,nprox(ipn))+1
     im1 = prox(im1,ipn)
     ip1 = prox(ip1,ipn)
  
     do k=1,nvl    		! Loop over layers
      trc_edg (k,edg,ipn,type) =				&
!           (2.*(tracr(k,ipn,type)+tracr(k,ipx,type))		&
!              + tracr(k,im1,type)+tracr(k,ip1,type))*divby6
            (8.*(tracr(k,ipn,type)+tracr(k,ipx,type))		&
               + tracr(k,im1,type)+tracr(k,ip1,type))*divby18
     end do			!  loop over layers
    end do			!  loop over edges
   end do			!  horizontal loop
!SMS$HALO_COMP END
  
!  if (type.eq.1) then
!   write (string,'(a,i2)') '(atm edgvar) tracer',type
!   do k=1,nvl,7
!    work1d(:)=tracr(k,:,type)
!    work2d(:,:)=trc_edg(k,:,:,type)
!    call stenedg(work1d,work2d,1,trim(string)//', cell & edge')
!   end do
!  end if
  
  end do			!  loop through tracers
  
!  sms$compare_var(u_edg, "load_ls.F90 - u_edg ")
!SMS$PARALLEL END
  call IncrementTimer(t1,tedgvar)
  
  return
  end subroutine edgvar
end module module_edgvar
