module module_output
  implicit none
!*********************************************************************
!	output
!	Output program for fim global model
!	Alexander E. MacDonald  12/27/2004
!	J. Lee                  September, 2005
!*********************************************************************

contains

subroutine output (its, nts, &                       ! index time step, final timestep
                   us3d, vs3d, dp3d,               & ! west wind, south wind, delta pres 
                   pr3d, ex3d, mp3d,               & ! pressure, Exner, mont pot, 
                   tr, rh3d, vor, ws3d,            & ! tracers, specific and relative humidity, etc.
                   chem_opt,diaga, diagb,          & ! diagnostic arrays
                   ph3d, tk3d, rn2d, rc2d, pw2d,   &
                   ts2d, us2d, hf2d, qf2d, sw2d,   &
                   lw2d, st3d, sm3d, t2m2d, q2m2d, & ! geopotential, accumulated precip/rainfall
                   canopy2d, fice2d, hice2d,       &
                   sheleg2d, slmsk2d,              &
                   u10m, v10m, flxlwtoa2d,         &
                   rn2d0, rc2d0, rg2d0, flxswavg2d, flxlwavg2d, & ! accumulated things
                   toutputBa, TimingBarriers, curr_write_time)
  use module_constants
  use module_control,           only: dt,filename_len,FixedGridOrder,nip,ntra,ntrb,&
    nvar2d,nvarp,nvl,nvlp,nvlp1,ArchvStep,ArchvTimeUnit,PrintDiags, restart_freq,&
    yyyymmddhhmm, hrs_in_month,EnKFIO, readrestart, itsstart
  use module_core_setup,        only: use_write_tasks
  use module_op_diag,           only: op_diag
  use module_outFMTed,          only: outFMTed
  use module_outqv,             only: outqv
  use module_outqv_mn_lat,      only: outqv_mn_lat
  use module_outqv_mn_lat_land, only: outqv_mn_lat_land
  use module_outvar_enkf,       only: outvar_enkf
  use module_printMAXMIN,       only: printMAXMIN
  use restart,                  only: write_restart
!SMS$IGNORE BEGIN
  use icosio,                   only: icosio_out
!SMS$IGNORE END
  use findmaxmin2
  use module_header,            only: header

  implicit none

  ! External variable declarations:
  integer, intent(in) :: its      ! current time step
  integer, intent(in) :: nts      ! final time step
  integer, intent(in) :: chem_opt ! for chem pressure level files we need to know chem_opt
!SMS$DISTRIBUTE (dh,nip) BEGIN
  real   ,intent(IN)    :: us3d(nvl  ,nip),vs3d(nvl,nip),dp3d(nvl,nip)
  real   ,intent(IN)    :: pr3d(nvlp1,nip)
  real   ,intent(IN)    :: mp3d(nvl  ,nip)
  real   ,intent(IN)    :: diaga(nvl  ,nip),diagb(nvl,nip)
  real   ,intent(IN)    :: vor (nvl  ,nip)
  real   ,intent(IN)    :: ws3d(nvl  ,nip)
  real   ,intent(IN)    :: ph3d(nvlp1,nip), ex3d(nvlp1,nip)
  real   ,intent(IN)    :: tr(nvl    ,nip,ntra+ntrb)
  real   ,intent(INOUT) :: rh3d(nvl  ,nip)
  real   ,intent(INOUT) :: tk3d(nvl  ,nip)
  real   ,intent(IN)    :: rn2d(nip),rc2d(nip)
  real   ,intent(IN)    :: u10m(nip),v10m(nip)
!JR Moved these 5 things to arg list so they can be written to the restart file.
  real   ,intent(inout) :: rn2d0(nip),rc2d0(nip),rg2d0(nip),flxswavg2d(nip),flxlwavg2d(nip)
  real   ,intent(INOUT) :: pw2d(nip)
  real   ,intent(IN)    :: ts2d(nip),us2d(nip),hf2d(nip),qf2d(nip),sw2d(nip),&
    lw2d(nip),st3d(4,nip),sm3d(4,nip)
  real   ,intent(IN)    :: flxlwtoa2d(nip)
  real   ,intent(IN)    :: t2m2d(nip),q2m2d(nip)
!JR Added items from module_sfc_variables so everything comes from input arg list
!JR rather than some used from module
  real, intent(in) :: canopy2d(nip)
  real, intent(in) :: fice2d(nip)
  real, intent(in) :: hice2d(nip)
  real, intent(in) :: sheleg2d(nip)
  real, intent(in) :: slmsk2d(nip)

  real*8 ,intent(INOUT) :: toutputBa
  logical,intent(IN)    :: TimingBarriers
  integer, intent(inout) :: curr_write_time ! most recent time vars. were written

  real                  :: th3d(nvl,nip), qv3d(nvl,nip), qw3d(nvl,nip), hfop(nip), qfop(nip)
  real                  :: oz3d(nvl,nip)
  real                  :: td3d(nvl,nip),mslp(nip),rg2d(nip)
  real                  :: rn_xh(nip),rc_xh(nip),rg_xh(nip)
  real                  :: g3p (nvlp,nip,nvarp)
  !                        (nvarp=5)
  !                        1=height,2=temp,3=RH (w.r.t. water),4=u wind,5=v wind
  real                  :: g2d (nip,nvar2d)
  !                     ! additional diagnostic 2d variables from op_diag.F90
  real                  :: spd10m_dif(nip)
!SMS$DISTRIBUTE END

  integer :: LB, time, ipn
  integer :: accum_start     ! value of "time" from previous output call
  real*8  :: t0, t1=0.0d0, t2=0.0d0, t3=0.0d0
  character(len=filename_len), external :: filename

  integer, external :: its2time

  time = its2time(its)

  if (its == 0) then
    flxswavg2d(:) = 0.
    flxlwavg2d(:) = 0.
    rn2d0(:) = 0.0
    rc2d0(:) = 0.0
    rg2d0(:) = 0.0
  end if

! --- accumulate fields averaged over "ArchvIntvl"
!SMS$PARALLEL (dh,ipn) BEGIN
  do ipn=1,nip
    flxswavg2d(ipn) = flxswavg2d(ipn) + sw2d(ipn)
    flxlwavg2d(ipn) = flxlwavg2d(ipn) + lw2d(ipn)
  end do
!SMS$PARALLEL END

!JR Write history info every so often (mod(its,archvstep) == 0)
  if (mod(its,ArchvStep) == 0) then
    accum_start = curr_write_time
    curr_write_time = time
    if (TimingBarriers) then
      call StartTimer(t0)
!SMS$BARRIER
      call IncrementTimer(t0,toutputBa)
    endif

    call StartTimer(t0)
    !there is an SMS problem writing tr(:,:,1)

    th3d(:,:) = tr(:,:,1)
    qv3d(:,:) = tr(:,:,2)
    qw3d(:,:) = tr(:,:,3)
    oz3d(:,:) = tr(:,:,4)
    rg2d(:)   = rn2d(:) - rc2d(:)
    hfop(:)   = hf2d(:) * 1.25 * 1004.0
    qfop(:)   = qf2d(:) * 1.25 * 2.5e6
    ! Calculate x-hour interval precip (difference from total precip at the last
    ! output time)
    rn_xh(:)  = rn2d(:) - rn2d0(:)
    rc_xh(:)  = rc2d(:) - rc2d0(:)
    rg_xh(:)  = rg2d(:) - rg2d0(:)
    rn2d0(:)  = rn2d(:)
    rc2d0(:)  = rc2d(:)
    rg2d0(:)  = rg2d(:)
    spd10m_dif(:) = sqrt(u10m(:)**2+v10m(:)**2) - sqrt(us3d(1,:)**2 + vs3d(1,:)**2)
    flxswavg2d(:) = flxswavg2d(:) / ArchvStep
    flxlwavg2d(:) = flxlwavg2d(:) / ArchvStep

    ! Calculate various 3-d and 2-d diagnostic variables for outputting below.
    ! These are generally multivariate diagnostics, thus not do-able by the
    ! scalar FIMpost, which can only do horizontal interpolation (icos to
    ! lat/lon) one variable at a time.
    call op_diag(              &
      its,nts,                 & ! index time step, final timestep
      us3d,vs3d,dp3d,          & ! west wind, south wind, delta pres 
      pr3d,ex3d,mp3d,          & ! pressure, Exner, mont pot, 
      tr,vor,ws3d,             & ! tracers, etc.
      ph3d,rn2d,rc2d,          & ! geopotential, accumulated precip/rainfall
      ts2d,us2d,hf2d,qf2d,sw2d,lw2d,st3d,sm3d,&
      ! Below are output variables from op_diag
      tk3d,rh3d,td3d,pw2d,mslp,g3p,g2d,time&
      )   

    write (6,*) ' U wind comp at k=1 - max/min'
    call outqv(us3d ,deg_lat,deg_lon,nip,1,1.)
    call outqv_mn_lat(us3d,deg_lat,deg_lon,50.,1,nip,its,1.)
    write (6,*) ' U10M - max/min'
    call outqv(u10m ,deg_lat,deg_lon,nip,1,1.)
    call outqv_mn_lat(u10m,deg_lat,deg_lon,50.,1,nip,its,1.)
    write (6,*) ' V wind comp at k=1 - max/min'
    call outqv(vs3d ,deg_lat,deg_lon,nip,1,1.)
    call outqv_mn_lat(vs3d,deg_lat,deg_lon,50.,1,nip,its,1.)
    write (6,*) ' V10M - max/min'
    call outqv(v10m ,deg_lat,deg_lon,nip,1,1.)
    call outqv_mn_lat(v10m,deg_lat,deg_lon,50.,1,nip,its,1.)
    write (6,*) ' spd_dif(wind-10M- wind(k=1)) - max/min'
    call outqv(spd10m_dif ,deg_lat,deg_lon,nip,1,1.)
    call outqv_mn_lat(spd10m_dif,deg_lat,deg_lon,50.,1,nip,its,1.)

    write (6,*) ' Snow water equivalent - max/min'
    call outqv(sheleg2d ,deg_lat,deg_lon,nip,1,1.)
    call outqv_mn_lat(sheleg2d,deg_lat,deg_lon,50.,1,nip,its,1.)
    write (6,*) ' Canopy water - max/min'
    call outqv(canopy2d ,deg_lat,deg_lon,nip,1,1.)
    call outqv_mn_lat(canopy2d,deg_lat,deg_lon,30.,1,nip,its,1.)
    write (6,*) ' Canopy water - mean - land only'
    call outqv_mn_lat_land(canopy2d,deg_lat,deg_lon,30.,slmsk2d,1,1,nip,its,1.)
    write (6,*) ' Rain - since last output - max/min'
    call outqv(rn_xh,deg_lat,deg_lon,nip,1,1.)
    write (6,*) ' Rain-conv - since last output - max/min'
    call outqv(rc_xh,deg_lat,deg_lon,nip,1,1.)
    write (6,*) ' Rain-grid-scale - since last output - max/min'
    call outqv(rg_xh,deg_lat,deg_lon,nip,1,1.)
    ! write (6,*)' Rain0-  - max/min'
    ! call outqv       (  rn2d0 ,deg_lat,deg_lon,nip,1,1.)
    ! write (6,*)' Rain0-conv  - max/min'
    ! call outqv       (  rc2d0 ,deg_lat,deg_lon,nip,1,1.)
    if (its == itsStart-1) then
      write (6,*) ' Sea ice-hice - max/min'
      call outqv(hice2d,deg_lat,deg_lon,nip,1,1.)
      call outqv_mn_lat(hice2d,deg_lat,deg_lon,50.,1,nip,its,1.)
      write (6,*) ' Sea ice-fice - max/min'
      call outqv(fice2d,deg_lat,deg_lon,nip,1,1.)
      call outqv_mn_lat(fice2d,deg_lat,deg_lon,50.,1,nip,its,1.)
    end if
    write (6,*) ' Outgoing longwave radiation - max/min'
    call outqv(flxlwtoa2d,deg_lat,deg_lon,nip,1,1.)

    if (PrintDiags) then
      call outFMTed (its,pr3d,ph3d,us3d,vs3d,dp3d,mp3d,vor,tr,rh3d,tk3d,ws3d,&
                     st3d,sm3d,rn2d,pw2d,ts2d,us2d,hf2d,qf2d,sw2d,lw2d,time)
    else
      if (EnKFIO .and. .not. use_write_tasks .and. .not. FixedGridOrder) then
       ! write only fields needed for EnKF DA cycle to two files.
       ! 3d dynamical and surface fields. Only supported if no write task, and FixedGridOrder
       ! = .false.
        call outvar_enkf(time,pr3d,ex3d,us3d,vs3d,tr,ph3d)

      else

!JR Redo order for testing convenience to match specification in FIMnamelist ("wgrib" will give 
!JR same order). Optional argument "scalefactor" is scaling factor to be applied when GRIB file 
!JR is written. Optional argument "accum_start" specifies an accumulation start time different 
!JR from default when GRIB file is written.

        LB = LBOUND(g3p,2)
        call icosio_out (its, time, 'hgtP', g3p(1,LB,1), nvlp, filename('hgtP',its), header('hgtP',nvlp,its))
        call icosio_out (its, time, 'tmpP', g3p(1,LB,2), nvlp, filename('tmpP',its), header('tmpP',nvlp,its))
        call icosio_out (its, time, 'rp3P', g3p(1,LB,3), nvlp, filename('rp3P',its), header('rp3P',nvlp,its))
        call icosio_out (its, time, 'up3P', g3p(1,LB,4), nvlp, filename('up3P',its), header('up3P',nvlp,its))
        call icosio_out (its, time, 'vp3P', g3p(1,LB,5), nvlp, filename('vp3P',its), header('vp3P',nvlp,its))

        LB=LBOUND(g2d,1)
        call icosio_out (its, time, 'pr3D', pr3d, nvlp1,  filename('pr3D',its), header('pr3D',nvlp1,its))
        call icosio_out (its, time, 'ph3D', ph3d, nvlp1,  filename('ph3D',its), header('ph3D',nvlp1,its), scalefactor=1./9.8)
        call icosio_out (its, time, 'tk3D', tk3d, nvl,    filename('tk3D',its), header('tk3D',nvl,its))
        call icosio_out (its, time, 'td3D', td3d, nvl,    filename('td3D',its), header('td3D',nvl,its))
        call icosio_out (its, time, 'ws3D', ws3d, nvl,    filename('ws3D',its), header('ws3D',nvl,its))
        call icosio_out (its, time, 'rh3D', rh3d, nvl,    filename('rh3D',its), header('rh3D',nvl,its))
        call icosio_out (its, time, 'us3D', us3d, nvl,    filename('us3D',its), header('us3D',nvl,its))
        call icosio_out (its, time, 'vs3D', vs3d, nvl,    filename('vs3D',its), header('vs3D',nvl,its))
        call icosio_out (its, time, 'rn2D', rn2d, 1,      filename('2D__',its), header('rn2D',1,its), accum_start=accum_start)
        call icosio_out (its, time, 'rc2D', rc2d, 1,      filename('2D__',its), header('rc2D',1,its), accum_start=accum_start)
        call icosio_out (its, time, 'r12D', rn_xh,1,      filename('2D__',its), header('r12D',1,its))
        call icosio_out (its, time, 'r22D', rc_xh,1,      filename('2D__',its), header('r22D',1,its))
        call icosio_out (its, time, 'rg2D', rg2d, 1,      filename('2D__',its), header('rg2D',1,its), accum_start=accum_start)
        call icosio_out (its, time, 'pw2D', pw2d, 1,      filename('2D__',its), header('pw2D',1,its))
        call icosio_out (its, time, 'ts2D', ts2d, 1,      filename('2D__',its), header('ts2D',1,its))
        call icosio_out (its, time, 'us2D', us2d, 1,      filename('2D__',its), header('us2D',1,its))
        call icosio_out (its, time, 'hf2D', hfop, 1,      filename('2D__',its), header('hf2D',1,its))
        call icosio_out (its, time, 'qf2D', qfop, 1,      filename('2D__',its), header('qf2D',1,its))
        call icosio_out (its, time, 'sw2D', sw2d, 1,      filename('2D__',its), header('sw2D',1,its))
        call icosio_out (its, time, 'lw2D', lw2d, 1,      filename('2D__',its), header('lw2D',1,its))
        call icosio_out (its, time, 'ms2D', mslp, 1,      filename('2D__',its), header('ms2D',1,its))
        call icosio_out (its, time, 'sn2D', sheleg2d,  1, filename('2D__',its), header('sn2D',1,its))
        call icosio_out (its, time, 'cb2D', g2d(LB,1), 1, filename('2D__',its), header('cb2D',1,its))
        call icosio_out (its, time, 'ct2D', g2d(LB,3), 1, filename('2D__',its), header('ct2D',1,its))
        call icosio_out (its, time, 'u12D', u10m, 1,      filename('2D__',its), header('u12D',1,its))
        call icosio_out (its, time, 'v12D', v10m, 1,      filename('2D__',its), header('v12D',1,its))
        call icosio_out (its, time, 'rp2D', g2d(LB,5), 1, filename('2D__',its), header('rp2D',1,its))


!JR These aren't specified in default FIMnamelist

        call icosio_out (its, time, 'dp3D', dp3d,    nvl, filename('dp3D',its), header('dp3D',nvl,its))
        call icosio_out (its, time, 'mp3D', mp3d,    nvl, filename('mp3D',its), header('mp3D',nvl,its))
        call icosio_out (its, time, 'th3D', th3d,    nvl, filename('th3D',its), header('th3D',nvl,its))
        call icosio_out (its, time, 'qv3D', qv3d,    nvl, filename('qv3D',its), header('qv3D',nvl,its), scalefactor=1000.)
        call icosio_out (its, time, 'qw3D', qw3d,    nvl, filename('qw3D',its), header('qw3D',nvl,its), scalefactor=1000.)
        call icosio_out (its, time, 'oz3D', oz3d,    nvl, filename('oz3D',its), header('oz3D',nvl,its), scalefactor=1000.)
        call icosio_out (its, time, 'vo3D', vor,     nvl, filename('vo3D',its), header('vo3D',nvl,its))

        call icosio_out (its, time, 'da3D', diaga,   nvl, filename('da3D',its), header('da3D',nvl,its))
        call icosio_out (its, time, 'db3D', diagb,   nvl, filename('db3D',its), header('db3D',nvl,its))
        call icosio_out (its, time, 'cn2D', canopy2d,1,   filename('2D__',its), header('cn2D',1,its))
        call icosio_out (its, time, 'st3D', st3d,    4,   filename('2D__',its), header('st3D',4,its))
        call icosio_out (its, time, 'sm3D', sm3d,    4,   filename('2D__',its), header('sm3D',4,its))
        call icosio_out (its, time, 't22D', t2m2d,   1,   filename('2D__',its), header('t22D',1,its))
        call icosio_out (its, time, 'q22D', q2m2d,   1,   filename('2D__',its), header('q22D',1,its))
        call icosio_out (its, time, 'r32D', rg_xh,   1,   filename('2D__',its), header('r32D',1,its))
        call icosio_out (its, time, 'sa2d', flxswavg2d, 1,filename('2D__',its), header('sa2d',1,its))
        call icosio_out (its, time, 'la2d', flxlwavg2d, 1,filename('2D__',its), header('la2d',1,its))
        call icosio_out (its, time, 'ol2d', flxlwtoa2d, 1,filename('2D__',its), header('ol2d',1,its))
        flxswavg2d(:) = 0.
        flxlwavg2d(:) = 0.

      if(ntrb.gt.0 .and. chem_opt == 300)then
      LB=LBOUND(g3p,2)
      call icosio_out (its, time, 'so2P', g3p(1,LB,6), nvlp, filename('so2P',its), header('so2P',nvlp,its))
      call icosio_out (its, time, 'slfP', g3p(1,LB,7), nvlp, filename('slfP',its), header('slfP',nvlp,its))
      call icosio_out (its, time, 'dmsP', g3p(1,LB,8), nvlp, filename('dmsP',its), header('dmsP',nvlp,its))
      call icosio_out (its, time, 'msaP', g3p(1,LB,9), nvlp, filename('msaP',its), header('msaP',nvlp,its))
      call icosio_out (its, time, 'p25P', g3p(1,LB,10), nvlp, filename('p25P',its), header('p25P',nvlp,its))
      call icosio_out (its, time, 'bc1P', g3p(1,LB,11), nvlp, filename('bc1P',its), header('bc1P',nvlp,its))
      call icosio_out (its, time, 'bc2P', g3p(1,LB,12), nvlp, filename('bc2P',its), header('bc2P',nvlp,its))
      call icosio_out (its, time, 'oc1P', g3p(1,LB,13), nvlp, filename('oc1P',its), header('oc1P',nvlp,its))
      call icosio_out (its, time, 'oc2P', g3p(1,LB,14), nvlp, filename('oc2P',its), header('oc2P',nvlp,its))
      call icosio_out (its, time, 'd1sP', g3p(1,LB,15), nvlp, filename('d1sP',its), header('d1sP',nvlp,its))
      call icosio_out (its, time, 'd2sP', g3p(1,LB,16), nvlp, filename('d2sP',its), header('d2sP',nvlp,its))
      call icosio_out (its, time, 'd3sP', g3p(1,LB,17), nvlp, filename('d3sP',its), header('d3sP',nvlp,its))
      call icosio_out (its, time, 'd4sP', g3p(1,LB,18), nvlp, filename('d4sP',its), header('d4sP',nvlp,its))
      call icosio_out (its, time, 'd5sP', g3p(1,LB,19), nvlp, filename('d5sP',its), header('d5sP',nvlp,its))
      call icosio_out (its, time, 's1sP', g3p(1,LB,20), nvlp, filename('s1sP',its), header('s1sP',nvlp,its))
      call icosio_out (its, time, 's2sP', g3p(1,LB,21), nvlp, filename('s2sP',its), header('s2sP',nvlp,its))
      call icosio_out (its, time, 's3sP', g3p(1,LB,22), nvlp, filename('s3sP',its), header('s3sP',nvlp,its))
      call icosio_out (its, time, 's4sP', g3p(1,LB,23), nvlp, filename('s4sP',its), header('s4sP',nvlp,its))
      call icosio_out (its, time, 'p10P', g3p(1,LB,24), nvlp, filename('p10P',its), header('p10P',nvlp,its))
      endif
     endif
    endif

    call IncrementTimer(t0,t1)
    call StartTimer(t0)
    call printMAXMIN(its,nvl,nip,ntra+ntrb,tr,dp3d)
    call IncrementTimer(t0,t2)

    if (TimingBarriers) then
      call StartTimer(t0)
!SMS$BARRIER
      call IncrementTimer(t0,toutputBa)
    endif
  endif

  if (its == itsStart+nts-1) then
    print "(' OUTPUT time, maxmin time, restart time:',2F10.1)", t1, t2
  endif

  return
end subroutine output
end module module_output
