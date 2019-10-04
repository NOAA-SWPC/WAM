module module_fim_dyn_run

contains

subroutine dyn_run(its)
!*********************************************************************
!	"Run" method for fim global model dynamics
!	Alexander E. MacDonald  12/24/2005
!	J. LEE                  12/28/2005
!*********************************************************************

use module_constants
use module_control  ,only: nts,nvl,nvlp1,nip,ntra,ntrb,itsStart,nabl,    &
                           PrintIpnDiag,TestDiagProgVars,TestDiagNoise,  &
                           PrintDiagProgVars,TimingBarriers,dt,dtratio,  &
                           ArchvTimeUnit,ArchvStep
                           
use module_variables,only: u_tdcy,v_tdcy,dp_tdcy,dpl_tdcy,               &
                           massfx,cumufx,dpinit,nf,of,vof,               &
                           us3d,vs3d,dp3d,pr3d,ph3d,ex3d,mp3d,           &
                           tk3d,u_edg,v_edg,dp_edg,diaga,diagb,          &
                           lp_edg,bnll_edg,adbash1,adbash2,adbash3,      &
                           tr3d,trc_edg,trdp,trl_tdcy,trc_tdcy,          &
                           ws3d,sdot,rh3d,pw2d,vor,psrf,ptdcy,worka,     &
                           curr_write_time
                                       
!SMS$IGNORE BEGIN
use module_initial_chem_namelists, only: chem_opt,cu_physics,mp_physics
use postdata                      ,only: gribout
USE gfs_physics_internal_state_mod, only:gis_phy


!SMS$IGNORE END
use module_wrf_variables,only: phys3dwrf,phys2dwrf,exch
use module_sfc_variables   ,only: rn2d,rc2d,ts2d,us2d,hf2d,qf2d,   &
                                  sw2d,lw2d,flxlwtoa2d,st3d,sm3d,t2m2d,q2m2d, &
                                  canopy2d, fice2d, hice2d, sheleg2d, slmsk2d, &
                                  rn2d0, rc2d0, rg2d0, flxswavg2d, flxlwavg2d
use module_abstart         ,only: abstart
use module_chem_output     ,only: chem_output
use module_wrf_output      ,only: wrf_output
use module_output          ,only: output
use module_edgvar          ,only: edgvar
use module_cnuity          ,only: cnuity
use module_trcadv          ,only: trcadv
use module_transp3d        ,only: transp0,transp1,transp2
use module_momtum          ,only: momtum
use module_hystat          ,only: hystat
use module_hybgen          ,only: hybgen 
use module_diagnoise       ,only: diagnoise
use module_globsum         ,only: globsum, qmass, qmsqv, qmsqw, qmste, qmstr, qmstrn, qmstrc, &
                                  qdtr, qdtrn, qdtrc
use module_profout         ,only: profout
use module_outDiags        ,only: outDiags
use module_outtime_dyn     ,only: tdyn,tout,thystat,tedgvar,tcnuity,    &
                                  ttrcadv,thybgen,tprofout,tabstart,    &
                                  tmomtum,ttransp,tedgvarEx,tcnuityEx,  &
                                  ttrcadvEx,ttranspEx,tedgvarBa,        &
                                  tcnuityBa,ttrcadvBa,toutputBa,        &
                                  ttranspBa
use module_chem_variables  ,only: aod2d,p10,pm25,trfall,tr1_tavg
use findmaxmin3
use module_core_setup      ,only: iam_compute_root,use_write_tasks

implicit none

!  Declare dummy arguments
integer, intent(in   ) :: its

!  Declare local variables:
integer :: itsm1        ! its - 1
integer :: ivl,ipn,type
real*8  :: t0,t1
! diagnostic sums
real    :: qtrcr(ntra+ntrb)
logical, save :: qdtr_print=.false.

real :: time
integer::ret
integer, external :: its2time
logical::post_init_file_called=.false.

!SMS$DISTRIBUTE(dh,1) BEGIN
real  :: gu10m(nip),gv10m(nip)
!SMS$DISTRIBUTE END

call StartTimer(t0)
itsm1 = its - 1

  !...........................................................
  ! Finish up previous dynamics time step unless this is the 
  ! first time step.  
  ! This complexity is required for the NCEP ESMF approach 
  ! in which single-phase DYN and PHY components alternate 
  ! execution during each time step.  
  !
!TBH:  Restore if statement label once Mark fixes PPP
!TBH skip_first_time1: if (its > itsStart ) then
if (its > 1 ) then

!sms$compare_var(st3d   , "dyn_run begin second half of iteration ")
!sms$compare_var(sm3d   , "dyn_run begin second half of iteration ")
!sms$compare_var(rn2d   , "dyn_run begin second half of iteration ")
!sms$compare_var(rc2d   , "dyn_run begin second half of iteration ")
!sms$compare_var(ts2d   , "dyn_run begin second half of iteration ")
!sms$compare_var(us2d   , "dyn_run begin second half of iteration ")
!sms$compare_var(hf2d   , "dyn_run begin second half of iteration ")
!sms$compare_var(sw2d   , "dyn_run begin second half of iteration ")

  !...........................................................
  ! Hybrid sigma-theta grid maintenance
  !
  call StartTimer(t1)
  call hybgen (itsm1,		&
   thetac,			& ! target pot.temp.
   us3d,vs3d,tr3d,		& ! zonal & merid wind, mass field tracers
   sdot,ex3d,dp3d,pr3d  )	  ! interface displ, exner, layer thknss, prs

  call IncrementTimer(t1,thybgen)

  !...........................................................
  ! Hydrostatic calculations
  call StartTimer(t1)
  call hystat (itsm1,		&
   ph3d,			& ! geopotential
   ex3d,mp3d,dp3d,		& ! exner fct, montg.pot., layer thickness
   tr3d,trdp,			& ! mass field tracers, tracer x thickness
   psrf,ptdcy    )		  ! surface pressure, srf.pres.tndcy
  call IncrementTimer(t1,thystat)

  !...........................................................
  ! Evaluate noise parameter
  !
  if( mod(itsm1,TestDiagNoise)==0.or.itsm1==itsStart ) then
    call diagnoise (itsm1,	&
     ptdcy                      ) ! sfc.pres.tendency at 2 consec. time levels
  endif
  if (mod (itsm1,TestDiagProgVars) == 0 .or. itsm1 == itsStart .or. its == itsStart+nts) then
    call globsum (itsm1, dp3d, tr3d, rn2d, rc2d, pr3d, ex3d, qf2d, qtrcr)
  endif 

  if (mod(itsm1,TestDiagProgVars) == 0 .or. (itsm1 == itsStart .and. PrintDiagProgVars /= -2)) then
    call outDiags (itsm1, nvl, nvlp1, nip, ntra+ntrb, pr3d, ph3d, tr3d, rn2d, &
                   rc2d, sdot, dp3d, us3d, vs3d, rh3d, lw2d, sw2d, hf2d, qf2d)
  endif
!sms$compare_var(mp3d, "dyn_run.F90 - mp3d10 ")

!sms$compare_var(st3d   , "dyn_run end iteration ")
!sms$compare_var(sm3d   , "dyn_run end iteration ")
!sms$compare_var(rn2d   , "dyn_run end iteration ")
!sms$compare_var(rc2d   , "dyn_run end iteration ")
!sms$compare_var(ts2d   , "dyn_run end iteration ")
!sms$compare_var(us2d   , "dyn_run end iteration ")
!sms$compare_var(hf2d   , "dyn_run end iteration ")
!sms$compare_var(sw2d   , "dyn_run end iteration ")

endif		! its > 1

!TBH:  Restore if statement label once Mark fixes PPP
!TBH endif skip_first_time1

call IncrementTimer(t0,tdyn)
call StartTimer(t0)
! Always call output routine, even for its-1 == 0
!............................................................

! If gribout is enabled and I am the compute root, prepare the grib output file
! for writing. The unsupported case where gribout is enabled and more than one
! write task is specified has already been handled in dyn_init. Note that it is
! assumed here that all routines calling icosio_out are doing so on the same
! time step. If any routine needs to write history to disk on a different
! schedule, the first "if" conditional below will need to change.

if (mod(itsm1,ArchvStep) == 0) then
  if (gribout.and..not.use_write_tasks.and.iam_compute_root()) then
    call post_init_file(its2time(itsm1),ret)
    if (ret /= 0) then
      write(6,*) 'output: bad return from post_init_file: stopping'
      call flush(6)
      stop
    endif
    post_init_file_called = .true.
  endif
endif

!SMS$PARALLEL(dh, ipn) BEGIN
if (itsm1==0) then
 ! physics has not yet been called
 do ipn=1,nip
   gu10m(ipn) = us3d(1,ipn)
   gv10m(ipn) = vs3d(1,ipn)
 end do
else
 do ipn=1,nip
   gu10m(ipn) = gis_phy%flx_fld%u10m(ipn,1)
   gv10m(ipn) = gis_phy%flx_fld%v10m(ipn,1)
 end do
endif
!SMS$PARALLEL END

! Generate the output field of main variables:
call output (itsm1, nts,                     &
             us3d, vs3d, dp3d,               & ! west & south wind, layer thickness
             pr3d, ex3d, mp3d,               & ! pressure, exner fct, montg.pot.,
             tr3d, rh3d, vor, ws3d,          & ! mass field tracers, relative humidity
             chem_opt,diaga, diagb,          & ! diagnostic arrays
             ph3d, tk3d, rn2d, rc2d, pw2d,   &
             ts2d, us2d, hf2d, qf2d, sw2d,   &
             lw2d, st3d, sm3d, t2m2d, q2m2d, &
             canopy2d, fice2d, hice2d,       &
             sheleg2d, slmsk2d,              &
             gu10m, gv10m, flxlwtoa2d,       &
             rn2d0, rc2d0, rg2d0, flxswavg2d, flxlwavg2d, & ! accumulated things
             toutputBa, TimingBarriers, curr_write_time)

! if done like this it should be separated into wrfphys and chem
if (chem_opt.gt.0) then ! output chem variables
  call chem_output(itsm1,nts,aod2d,exch,p10,pm25,pr3d,tk3d,tr3d,trfall,phys2dwrf,tr1_tavg)
endif
! if done like this it should be separated into wrfphys and chem
if ( mp_physics .gt. 0 .or. cu_physics .gt. 0) then ! output wrfphys variables
  call wrf_output(itsm1,nts,pr3d,tk3d,tr3d,phys2dwrf)
endif
  
! See the comment accompanying the post_init_file() call.
        
if (post_init_file_called) then
  call post_finalize_file(ret)
  if (ret.ne.0) then
    write(6,*) 'output: bad return from post_finalize_file: stopping'
    stop
  endif
  post_init_file_called = .false.
endif
        
call IncrementTimer(t0,tout)
call StartTimer(t0)

!TBH:  Restore if statement label once Mark fixes PPP
!TBH skip_first_time2: if (its > itsStart ) then
if (its > 1) then
  call StartTimer(t1)
  call profout(                 &
   itsm1,PrintIpnDiag,          & ! index time step
   us3d,vs3d,dp3d,              & ! west & south wind, layer thickness
   pr3d,tr3d(:,:,1),mp3d,tk3d,  & ! pressure, pot.temp., montg.pot, temp (k)
   tr3d(:,:,2),rh3d,            & ! specific and relative humidity
   ph3d,                        & ! geopot.
   ts2d,us2d,hf2d,qf2d  )         ! skin temp., ustar, snsbl heatflx, vapor flx
  call IncrementTimer(t1,tprofout)

  if (mod (itsm1,TestDiagProgVars) == 0 .or. itsm1 == 1 .or. its == itsStart+nts) then
    time=its2time(itsm1)
    write (6,*)'precip - total/nonconv/conv=',qdtr,qdtrn,qdtrc
    write (6,80)'Global 3D mass          =',qmass,' at ',time,ArchvTimeUnit,', time step=',itsm1
    write (6,80)'Global 3D water vapor   =',qmsqv,' at ',time,ArchvTimeUnit,', time step=',itsm1
    write (6,80)'Global 3D cloud water   =',qmsqw,' at ',time,ArchvTimeUnit,', time step=',itsm1
    write (6,80)'Global integ acc precip =',qmstr,' at ',time,ArchvTimeUnit,', time step=',itsm1

    if(qdtr_print) then
      write (6,100) PrintDiagProgVars,qdtr,time,ArchvTimeUnit,itsm1
100   format ('Global precip, last',i3,'h =',es14.7,' at ',f5.1,1x,a,', time step=',i8)
    else
      qdtr_print = .true.
    endif
    write (6,*)'Global integ evaporation=',qmste,' at time step=',itsm1
    do type=1,ntrb
      write (6,'(a,i2.2,a,es15.7,a,i12)') 'Tracer B',type,        &
       ' global amount =',qtrcr(ntra+type),'  at time step=',itsm1
    end do
  endif 
80  format (a,es14.7,a,f5.1,1x,2a,i8)

!TODO:  move this down into a subroutine...  
  if (mp_physics /= 0 .or. cu_physics /= 0) then
!
! store this time level t,qv
!
!SMS$PARALLEL(dh, ipn) BEGIN
    do ipn=1,nip
      do ivl=1,nvl
        phys3dwrf(ivl,ipn,3) = tr3d(ivl,ipn,2)
        phys3dwrf(ivl,ipn,7) = tr3d(ivl,ipn,1)*ex3d(ivl,ipn)/(1.+.6078*tr3d(ivl,ipn,2))/cp
!     if(ipn.eq.5000)then
!        write(6,*)'dyn0',ivl,phys3dwrf(ivl,ipn,7),ex3d(ivl,ipn)
!     endif
      enddo
    enddo
!SMS$PARALLEL END
  endif
endif
!TBH:  Restore if statement label once Mark fixes PPP
!TBH endif skip_first_time2

  !...........................................................
  ! Begin current dynamics time step unless this is the last 
  ! (itsStart+nts+1) iteration.  
  ! This complexity is required for the NCEP ESMF approach 
  ! in which single-phase DYN and PHY components alternate 
  ! execution during each time step.  
  !
!TBH:  Restore if statement label once Mark fixes PPP
!TBH skip_last_iteration: if (its <= nts ) then
if (its < itsStart+nts) then
  print "('its=',i7)",its

!sms$compare_var(st3d   , "dyn_run begin iteration ")
!sms$compare_var(sm3d   , "dyn_run begin iteration ")
!sms$compare_var(rn2d   , "dyn_run begin iteration ")
!sms$compare_var(rc2d   , "dyn_run begin iteration ")
!sms$compare_var(ts2d   , "dyn_run begin iteration ")
!sms$compare_var(us2d   , "dyn_run begin iteration ")
!sms$compare_var(hf2d   , "dyn_run begin iteration ")
!sms$compare_var(sw2d   , "dyn_run begin iteration ")

  !.............................................................
  ! Start the Adams Bashforth 3rd order time diff:
  call StartTimer(t1)
  call abstart (its,itsStart, &
    nf,of,vof,			& ! time slots for Adams Bashforth
    adbash1,adbash2,adbash3	) ! Adams Bashforth time dif. weights
  call IncrementTimer(t1,tabstart)

  !.........................................................
  ! interpolate variables to cell edges
  call edgvar (its,		&
    us3d,vs3d,			& ! west & south wind
    dp3d,pr3d,			& ! layer thknss, pressure
    mp3d,tr3d,			& ! montg.pot, tracer
    u_edg,v_edg,dp_edg,lp_edg,	& ! edge variables: u,v,thknss,midlyr-prs
    bnll_edg,trc_edg,		& ! edge variables: bernoulli-fct, tracer
    tedgvar,tedgvarEx,tedgvarBa,& ! timers
    TimingBarriers )              ! turn on timing barriers when .true.

  !........................................................
  ! Solve continuity equation
  call cnuity (its,		&
    nf,of,vof,			& ! time slots for Adams Bashforth
    adbash1,adbash2,adbash3,	& ! Adams Bashforth time dif. weights
    us3d,vs3d,			& ! west & south wind
    u_edg,v_edg,		& ! west & south wind on edges
    dp_edg,lp_edg,		& ! lyr thknss & mid-lyr pressure on edges
    dp3d,pr3d,ex3d,		& ! layer thickness, pressure, Exner fct.
    dp_tdcy,dpl_tdcy,		& ! forcing for dp and for low ord dp
    massfx,ws3d,		& ! mass flux across edges, omega=dp/dt
    tcnuity,tcnuityEx,tcnuityBa,& ! timers
    TimingBarriers )              ! turn on timing barriers when .true.

  !........................................................
  ! Build up time integral of mass fluxes
  if (ntrb > 0) then
    call transp1(its,             &
         nf,of,vof,                  & ! time slots for Adams Bashforth
         adbash1,adbash2,adbash3,    & ! Adams Bashforth time dif. weights
         cumufx,massfx)              ! time-integrated and instant. massfx

  !........................................................
  ! Solve transport equation for class B tracers
    if (mod(its,dtratio) == 0) then
      call transp2(its,            &
           tr3d,cumufx,               & ! tracer, time-integrated mass flux
           dpinit,dp3d,               & ! initial & final thickness
           ttransp,ttranspEx,ttranspBa,TimingBarriers )

  ! re-initialize time integrals
      call transp0(its,cumufx,dp3d,dpinit)
    end if
  end if

  !........................................................
  ! Solve tracer transport equation for class A tracers
  call trcadv (its,		&
    nf,of,vof,			& ! time slots for Adams Bashforth
    adbash1,adbash2,adbash3,	& ! Adams Bashforth time dif. weights
    trc_edg,			& ! tracer on edges
    tr3d,trdp,			& ! trcr, trcr*dp, trcr*dp low ord
    trc_tdcy,trl_tdcy,		& ! forcing, low ord forcing for trcr
    massfx,dp3d,                & ! mass flux across edges, layer thickness
    ttrcadv,ttrcadvEx,ttrcadvBa,& ! timers
    TimingBarriers )              ! turn on timing barriers when .true.

  !.........................................................   
  ! Solve momentum equation
  call StartTimer(t1)
  call momtum (its,		&
    nf,of,vof,			& ! time slots for Adams Bashforth
    adbash1,adbash2,adbash3,	& ! Adams Bashforth time dif. weights
    us3d,vs3d,ex3d,vor,		& ! west & south wind, exner fct, vorticity
    u_edg,v_edg,trc_edg,	& ! u,v,trcr on edges
    bnll_edg,			& ! bernoulli fct on edges
    u_tdcy,v_tdcy,dp3d	)	  ! forcing for u,v; layer thknss
  call IncrementTimer(t1,tmomtum)

!!!sms$compare_var(u_tdcy, "dyn_run.F90 - u_tdcy5 ")
!sms$compare_var(us3d, "dyn_run.F90 - us3d5 ")
!sms$compare_var(trdp, "dyn_run.F90 - trdp5 ")
!sms$compare_var(trc_tdcy, "dyn_run.F90 - trc_tdcy5 ")

  !.........................................................
  ! Save theta at end of dyn_run but before physics_run
  worka(:,:)=tr3d(:,:,1)
  
!sms$compare_var(st3d   , "dyn_run end first half of iteration ")
!sms$compare_var(sm3d   , "dyn_run end first half of iteration ")
!sms$compare_var(rn2d   , "dyn_run end first half of iteration ")
!sms$compare_var(rc2d   , "dyn_run end first half of iteration ")
!sms$compare_var(ts2d   , "dyn_run end first half of iteration ")
!sms$compare_var(us2d   , "dyn_run end first half of iteration ")
!sms$compare_var(hf2d   , "dyn_run end first half of iteration ")
!sms$compare_var(sw2d   , "dyn_run end first half of iteration ")

endif
!TBH:  Restore if statement label once Mark fixes PPP
!TBH endif skip_last_iteration

!TODO:  move this down into a subroutine...  
if (mp_physics /=0 .or. cu_physics /= 0) then
!
! get dynamic tendencies for th and qv
!SMS$PARALLEL(dh, ipn) BEGIN
  do ipn=1,nip
    do ivl=1,nvl
!     if(ipn.eq.5000)then
!        write(6,*)'dyn',dt,tr3d(ivl,ipn,1)* (ex3d(ivl,ipn)/(1.+.6078*tr3d(ivl,ipn,2))/cp),phys3dwrf(ivl,ipn,7)
!        write(6,*)'dyn2',ivl,tr3d(ivl,ipn,2),phys3dwrf(ivl,ipn,3),ex3d(ivl,ipn)
!     endif
      phys3dwrf(ivl,ipn,3) = (tr3d(ivl,ipn,2) - phys3dwrf(ivl,ipn,3))/dt
      phys3dwrf(ivl,ipn,7) = (tr3d(ivl,ipn,1)*(ex3d(ivl,ipn)/(1. + .6078*tr3d(ivl,ipn,2)))/cp - &
                              phys3dwrf(ivl,ipn,7))/dt
    enddo
  enddo
!SMS$PARALLEL END
endif

call IncrementTimer(t0,tdyn)

return
end subroutine dyn_run
end module module_fim_dyn_run
