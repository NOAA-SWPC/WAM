module module_fim_dyn_init
  use findmaxmin2
  use stencilprint
  use stenedgprint
  use module_control  ,only: nvl,nvlp1,nip,nts,dt,numphr,TotalTime,ArchvIntvl, &
                             ArchvTimeUnit,ArchvStep,itsStart,                 &
                             readrestart,ntra,    &
                             npp,nd,glvl,PrintDiags,curve,ntrb,                &
                             PrintIpnDiag,PrintMAXMINtimes,FixedGridOrder,     &
                             NumCacheBlocksPerPE,ipnDiagLocal,ipnDiagPE,       &
                             yyyymmddhhmm,PrintDiagProgVars,TestDiagProgVars,  &
                             PrintDiagNoise,TestDiagNoise,control,             &
                             TimingBarriers,gfsltln_file,dt_reducer_numerator, &
                             dt_reducer_denominator,EnKFAnl,EnKFFileName,      &
                             FGFileName,FGFileNameSig,filename_len,varnamelen
  use module_constants,only: lat,lon,nprox,proxs,prox,area,cs,sn,sidevec_c,    &
                             sidevec_e,sideln,rprox_ln,deg_lat,deg_lon,rarea,  &
                             omegx2,raddeg,rsideln,corio,thetac,   &
                             dpsig,p1000,cp,rd,qvmin,qwmin,inv_perm,nedge,     &
                             permedge, actual
  use module_variables,only: u_tdcy,v_tdcy,dp_tdcy,trc_tdcy,trl_tdcy,          &
                             u_tdcy_phy,v_tdcy_phy,trc_tdcy_phy,               &
                             massfx,nf,of,vof,us3d,vs3d,ws3d,dp3d,tk3d,        &
                             pr3d,ex3d,tr3d,mp3d,ph3d,rh3d,trdp,vor,           &
                             sdot,dpl_tdcy,work2d,iwork2d,pw2d,psrf,ptdcy,     &
                             dpinit,cumufx
  use module_sfc_variables,only: rn2d,rc2d,qf2d,flxlwtoa2d
  use module_core_setup   ,only: iam_compute_root,iam_fim_task,iam_write_task, &
                                 my_comm,nct,nwt,output_intercomm,use_write_tasks
  use module_dyn_alloc    ,only: dyn_alloc
  use module_hybgen       ,only: hybgen
  use module_hystat       ,only: hystat
  use module_transp3d     ,only: transp0
  use module_header       ,only: header_size
  use module_wtinfo       ,only: wtinfo
!SMS$IGNORE BEGIN
  use icosio,   only: icosio_prep,icosio_set_inv_perm,icosio_setup
  use post,     only: post_init_readnl, post_init_slint, post_init_readnl_called, post_init_slint_called
  use postdata, only: fimout,gribout
  use sigio_module
  use units,    only: getunit, returnunit
!SMS$IGNORE END

  implicit none

contains

  subroutine dyn_init (client_server_io_in)
!*********************************************************************
!       Loads the initial variables and constants for the dynamics 
!       component.  
!       Alexander E. MacDonald  11/27/04
!       J. Lee                  September, 2005
!*********************************************************************

    logical, intent(in), optional :: client_server_io_in

! Local variables

!SMS$DISTRIBUTE(dh,1) BEGIN
    logical,allocatable :: in_halo(:)
    integer,allocatable :: itmp(:) ! For bounds-finding: See TODO below.
!SMS$DISTRIBUTE END

!SMS$DISTRIBUTE(dh,3) BEGIN
    real,allocatable :: work_edg(:,:,:)  ! work array for printing sidevec stencils
!SMS$DISTRIBUTE END

    character*16  :: string
    integer       :: ipn               ! Index for icos point number
    integer       :: isn               ! Index for icos edge number
    integer       :: ivl               ! Index for vertical level
    integer       :: i,its,idx
    integer       :: mype=-1           ! Set invalid default value
    integer       :: ipnGlobal
    logical       :: DiagPrint

    integer :: lusig
    integer(sigio_intkind) :: iret
    integer                :: nvp
    type(sigio_head)       :: head
    type(sigio_data)       :: data
    CHARACTER(len=9 )      :: jdate 
    CHARACTER(len=2 )      :: hh
    CHARACTER(len=80)      :: sanlFile

! IBM follows the Fortran standard precisely here requiring all 
! literal constants to have the same number of characters.  So 
! whitespace is significant.  
    character(11),dimension(0:3) :: CurveType = (/'IJ         ', &
                                                  'Hilbert    ', &
                                                  'IJ block   ', &
                                                  'SquareBlock'/)
    integer :: ipnx,isnx,ip1,im1,np1,nm1
    integer :: isncount
    integer :: nprocs = 1
    integer :: ios
!SMS$insert integer              :: HaloSize
!SMS$insert integer,allocatable  :: RegionSize(:)
!SMS$insert CHARACTER(len=80)    :: DecompInfoFile

    real*8 :: t0,t1=0.0D0

    integer :: ims,ips,ipe,ime    ! bounds: lower outer & inner, upper outer & inner
    integer :: icosio_comm        ! an MPI communicator to pass to icosio
    integer :: idum1, idum2       ! dummies
    logical :: client_server_io   ! run icosio in client-server mode?
    logical :: debugmsg_on        ! write-task debug message control
    logical :: ldum1, ldum2       ! dummies
    integer :: unitno             ! Unit number for I/O

    print *,'entering dyn_init ...'

!SMS$insert call nnt_me(mype)
!SMS$insert call NNT_NPROCS(nprocs)
    if (iam_write_task) nprocs = nct

    call control ()
    call wtinfo(idum1,nwt,idum2,ldum1,ldum2,debugmsg_on,my_comm)

! Initialize post. If gribout is enabled, ensure there is a max of 1 write task.

    call post_init_readnl(iret)
    if (iret /= 0) then
      write(6,*)'dyn_init: bad return from post_init: stopping'
      stop
    end if

    if (.not. post_init_readnl_called) then
      write(6,*)'dyn_init: logic error: cannot test gribout until post_init has been called'
      stop
    end if

    if (gribout .and. nwt > 1) then
      write(6,*)'dyn_init: Cannot have more than 1 write task when gribout enabled'
      stop
    end if

if (iam_write_task) then
  call post_init_slint (iret)
  if (iret /= 0) then
    write(6,*)'dyn_init: bad return from post_init_slint: stopping'
    stop
  endif
  if (.not. post_init_slint_called) then
    write(6,*)'dyn_init: logic error: cannot test gribout until post_init_slint has been called'
    stop
  endif
endif

!TODO:  remove dependence of serial prep code on ComputeTasks!  
!SMS$insert ios = 0 ! avoid Lahey complaint on non-root tasks
!SMS$insert allocate(RegionSize(nprocs))
!SMS$insert write(DecompInfoFile,"('DecompInfo_',i0,'.dat')") nprocs
!SMS$insert unitno = getunit ()
!SMS$insert if (unitno < 0) then
!SMS$insert   print*,'dyn_init: getunit failed for file=', trim(DecompInfoFile),'. Stopping'
!SMS$insert   stop
!SMS$insert end if
!SMS$insert open (unitno,file=TRIM(DecompInfoFile),status='old',iostat=ios)
!SMS$insert if (ios /= 0) then
!SMS$insert   print*,'ERROR:  Cannot find ',TRIM(DecompInfoFile),', prep and fim must be run with the same setting for ComputeTasks in FIMnamelist.'
!SMS$insert   stop
!SMS$insert endif
!SMS$insert read (unitno,*) HaloSize
!SMS$insert read (unitno,*) RegionSize
!SMS$insert close(unitno)
!SMS$insert call returnunit (unitno)
!SMS$CREATE_DECOMP(dh,<nip>,<HaloSize>:regionsize=RegionSize)

! Set icosio values.

    client_server_io = .true.
    if (present(client_server_io_in)) client_server_io = client_server_io_in

    icosio_comm = my_comm
    if (use_write_tasks) icosio_comm = output_intercomm

! TODO Similar bounds-finding logic involving a distributed itmp array already
! TODO exists in wrf_set_array_bounds(). Generalize that routine to allow finding
! TODO arbitrary -- memory, domain, patch or tile -- sets of bounds, and call it
! TODO here to find im[s|e] and ip[s|e].

    allocate(itmp(nip))
    ims = lbound(itmp,1)
    ime = ubound(itmp,1)
!sms$to_local(dh:<1,ips:lbound>,<1,nip:ubound>) begin
    ips = 1
    ipe = nip
!sms$to_local end

!JR Tell icosio to use unit 50. Since getunit() with no arguments may return 
!JR different numbers for different MPI tasks, want to be safe.
    unitno = getunit (50)
    if (unitno < 0) then
      print*,'dyn_init: getunit failed for icosio_setup. Stopping'
      stop
    end if

! Send icosio required values. If client/server io is enabled, write tasks will
! not return from the icosio_setup() call.

    call icosio_setup (binout_in=fimout,                    &
                       client_server_io_in=client_server_io,&
                       comm_in=icosio_comm,                 &
                       debugmsg_on_in=debugmsg_on,          &
                       dt_in=dt,                            &
                       filename_len_in=filename_len,        &
                       glvl_in=glvl,                        &
                       gribout_in=gribout,                  &
                       header_size_in=header_size,          &
                       i_am_write_task_in=iam_write_task,   &
                       ims_in=ims,                          &
                       ime_in=ime,                          &
                       ips_in=ips,                          &
                       ipe_in=ipe,                          &
                       lunout_in=unitno,                    &
                       nip_in=nip,                          &
                       nts_in=itsStart+nts-1,               &  ! terminal time step
                       nvl_in=nvl,                          &
                       print_diags_in=printdiags,           &
                       using_write_tasks_in=use_write_tasks,&
                       varname_len_in=varnamelen,           &
                       yyyymmddhhmm_in=yyyymmddhhmm)

! If client/server IO is disabled, write tasks return from the icosio_setup()
! call when they receive a shutdown command from the compute root after the last
! output frame has been processed, and return from the present subroutine at
! the following statement. FIM uses client/server IO by default, and disables it
! only for NEMS-enabled runs.

    if (iam_write_task) return

    call dyn_alloc ()
    allocate(in_halo(nip))

! Changes in this block may require changes in output.F90 as well
    ArchvIntvl = min (ArchvIntvl, TotalTime)

    print"('FIM Global Model')"
    print *, " "
    call datetime ()
    print *, " "
    print"(' Curve:    '                           ,A44              )",CurveType(curve)
    print"(' Number of cache blocks per processor:',I8,' blocks'     )",NumCacheBlocksPerPE
    print"(' Grid Level'                           ,I35              )",glvl
    print"(' Number of Processors:'                ,I24,' processors')",nprocs
    print"(' Global size:     '                    ,I28,' points'    )",nip
!SMS$insert print"(' Halo size:       '                    ,I28,' points'    )",HaloSize
    print"(' Forecast duration (',a2,'):'          ,I23)",ArchvTimeUnit,TotalTime
    print"(' Vertical resolution:'                 ,I25,' levels'    )",nvl

    if (glvl >=7 ) then
      print"(' Default time step:'                   ,I27,' seconds'   )",nint(dt*dt_reducer_denominator/dt_reducer_numerator)
      print"(' Time step reduced to '                ,I24,' seconds'   )",nint(dt)
    else
      print"(' Length of Time step: '                ,I24,' seconds'   )",nint(dt)
    end if

    print"(' Number of time steps:'                ,I24,' timesteps' )",nts
    print "(' Output every',I33,' timesteps')",ArchvStep
    print  "(' Print MAX,MIN routine times',L18)",PrintMAXMINtimes
    print  "(' Add timed barriers to measure task skew',L6)",TimingBarriers
    print  "(' Output in fixed (IJ) grid ordering',L11)",FixedGridOrder

    if (curve == 0) then !The grid order is already IJ for curve=0.
      FixedGridOrder = .false.
    end if

    if      (PrintDiagProgVars == 0) then
      TestDiagProgVars = 1
    else if (PrintDiagProgVars <  0) then
      TestDiagProgVars = 2*nts
    else
      TestDiagProgVars = PrintDiagProgVars*numphr
    end if

    if    (PrintDiagNoise == 0) then
      TestDiagNoise = 1
    elseif(PrintDiagNoise <  0) then
      TestDiagNoise = 2*nts
    else
      TestDiagNoise = PrintDiagNoise*numphr
    end if

    print "(' Forecast initial time        ',A16,' YYYYMMDDHHMM')",yyyymmddhhmm
    print "(' Diagnostic prints at ipn           ',I10)",PrintIpnDiag
    print "(' Print diagnostic messages          ',L10)",PrintDiags
    print "(' Print diagnostic prognosis vars    ',I10)",PrintDiagProgVars
    print "(' Print diagnostic gravity wave noise',I10)",PrintDiagNoise
    print *,' '
! TODO:  Remove the print of CompilerMPIstring ASAP.  It is no longer used by 
! TODO:  FIMrun or FIMwfm scripts.  (Should be OK to remove if no one has 
! TODO:  complained within a few weeks after this commit!)  
    print*,' CompilerMPIstring=UPDATE_SCRIPTS_TO_READ_BUILD_CONFIG_FROM_SRCDIR='
    print *,' '

    if (PrintIpnDiag > nip) then
      print*,'Fatal error: PrintIpnDiag must be <= nip',PrintIpnDiag,nip
      stop
    end if

! ..................................................................
! Sec. 1. Calculate icosahedral grid description data
! ..................................................................
! The variables to be read in are:

! lat(nip),lon(nip) - lat and lon of icos pts
! nprox(nip) = number of proximity points (6 or 5)
! proxs(npp,nip) = array holding icos side number (1 to 6)
! prox(npp,nip) = array holding index of icos cell across side
! area(nip) = area of cell
! sidevec_c(nd,npp,nip)= side vector projected from cell center
! sidevec_e(nd,npp,nip)= side vector projected from cell edge
! sideln(npp,nip) = length of side (edge)- local projection invarient
! rprox_ln(npp,nip) = reciprical of length between icos pts
! inv_perm = inverse permutation of the grid
!SMS$PARALLEL(dh, ipn) BEGIN
!SMS$SERIAL BEGIN
    unitno = getunit ()
    if (unitno < 0) then
      print*,'dyn_init: getunit failed for file=glvl.dat. Stopping'
      stop
    end if

    open (unitno, file='glvl.dat', form='unformatted', action='read', iostat=ios)
    if (ios /= 0) then
      print*,'dyn_init: cannot open glvl.dat for reading. Stopping'
      stop
    end if

    call TestGlvlHeader (unitno, 'glvl.dat', 'dyn_init', glvl)
    call TestCurveHeader(unitno, 'glvl.dat', 'dyn_init', curve)

!SMS$SERIAL END
!SMS$SERIAL BEGIN
    read (unitno,err=90) lat
!SMS$SERIAL END
!SMS$SERIAL BEGIN
    read (unitno, err=90) lon
!SMS$SERIAL END
!SMS$SERIAL BEGIN
    read (unitno, err=90) nprox
!SMS$SERIAL END
    do isn=1,size(proxs,1)
!SMS$SERIAL BEGIN
      read (unitno, err=90) iwork2d
!SMS$SERIAL END
      do ipn=1,nip
        proxs(isn,ipn) = iwork2d(ipn)
      end do
    end do

    do isn=1,size(prox,1)
!SMS$SERIAL BEGIN
      read (unitno, err=90) iwork2d
!SMS$SERIAL END
      do ipn=1,nip
        prox(isn,ipn) = iwork2d(ipn)
      end do
    end do
!SMS$SERIAL BEGIN
    read (unitno, err=90) area
!SMS$SERIAL END
    do isn=1,size(cs,2)
      do idx=1,size(cs,1)
!SMS$SERIAL BEGIN
        read (unitno, err=90) work2d
!SMS$SERIAL END
        do ipn=1,nip
          cs(idx,isn,ipn) = work2d(ipn)
        end do
      end do
    end do

    do isn=1,size(sn,2)
      do idx=1,size(sn,1)
!SMS$SERIAL BEGIN
        read (unitno, err=90) work2d
!SMS$SERIAL END
        do ipn=1,nip
          sn(idx,isn,ipn) = work2d(ipn)
        end do
      end do
    end do

    do isn=1,size(sidevec_c,2)
      do idx=1,size(sidevec_c,1)
!SMS$SERIAL BEGIN
        read (unitno, err=90) work2d
!SMS$SERIAL END
        do ipn=1,nip
          sidevec_c(idx,isn,ipn) = work2d(ipn)
        end do
      end do
    end do

    do isn=1,size(sidevec_e,2)
      do idx=1,size(sidevec_e,1)
!SMS$SERIAL BEGIN
        read (unitno, err=90) work2d
!SMS$SERIAL END
        do ipn=1,nip
          sidevec_e(idx,isn,ipn) = work2d(ipn)
        end do
      end do
    end do

    do isn=1,size(sideln,1)
!SMS$SERIAL BEGIN
      read (unitno, err=90) work2d
!SMS$SERIAL END
      do ipn=1,nip
        sideln(isn,ipn) = work2d(ipn)
      end do
    end do

    do isn=1,size(rprox_ln,1)
!SMS$SERIAL BEGIN
      read (unitno, err=90) work2d
!SMS$SERIAL END
      do ipn=1,nip
        rprox_ln(isn,ipn) = work2d(ipn)
      end do
    end do
!SMS$SERIAL BEGIN
    read (unitno, err=90) inv_perm
    close(unitno)
    call returnunit (unitno)
!SMS$SERIAL END
    do ipn=1,nip
      if(nprox(ipn)==5)then
        prox(6,ipn) = prox(5,ipn)
      end if
    end do

!SMS$PARALLEL END
!SMS$UNSTRUCTURED_GRID(PROX)

! Update halos of constant arrays.  Cannot rely on automatic halo update from 
! file read due to prox not yet being set up prior to !SMS$UNSTRUCTURED_GRID.  
!SMS$EXCHANGE(lat,lon,nprox,proxs,area,cs,sn,sidevec_c,sidevec_e,sideln,rprox_ln)

  ! set up nedge and permedge for computations in the halo (HALO_COMP)
!SMS$IGNORE BEGIN
    in_halo(:) = .true.
    nedge(:) = 0
!SMS$IGNORE END
!SMS$PARALLEL(dh, ipn) BEGIN
! set up nedge and permedge for interior cells
    do ipn=1,nip
      in_halo(ipn) = .false.
      nedge(ipn) = nprox(ipn)     ! NOOP for interior cells
      do isn=1,nprox(ipn)
        permedge(isn,ipn) = isn   ! NOOP for interior cells
      end do
    end do
! now set up nedge and permedge for halo cells
! NOTE:  *not* owner-computes!
    do ipn=1,nip
      do isn=1,nprox(ipn)
        ipnx = prox(isn,ipn)
        if (in_halo(ipnx)) then
          nedge(ipnx) = nedge(ipnx) + 1
          permedge(nedge(ipnx),ipnx) = proxs(isn,ipn)
        end if
      end do
    end do
!SMS$PARALLEL END

#define DEBUG_HALO_COMP
!TBH:  lots of error checks for debugging
#ifdef DEBUG_HALO_COMP
!SMS$PARALLEL(dh, ipn) BEGIN
! verify that in_halo is not screwed up in the interior
    do ipn=1,nip
      if (in_halo(ipn)) then
!SMS$IGNORE BEGIN
        print *,'ERROR C:  in_halo(',ipn,') = ',in_halo(ipn),'  but ipn is an interior cell!'
!SMS$IGNORE END
        stop
      end if
    end do
! verify that prox(isnx,ipnx) == ipn in the interior
    do ipn=1,nip
      do isn=1,nprox(ipn)
        ipnx = prox(isn,ipn)
        isnx = proxs(isn,ipn)
        if (.not.in_halo(ipnx)) then  ! avoid halo points which have not yet been set up
          if (prox(isnx,ipnx) /= ipn) then
!SMS$IGNORE BEGIN
            print *,'ERROR C:  prox(',isnx,',',ipnx,') /= ',ipn,'  me = ',mype,prox(isnx,ipnx)
!SMS$IGNORE END
            stop
          end if
        end if
      end do
    end do
! verify that nedge is OK in the interior
    do ipn=1,nip
      if (nedge(ipn) /= nprox(ipn)) then
!SMS$IGNORE BEGIN
        print *,'ERROR C:  nedge(',ipn,') /= nprox(',ipn,')  [',nedge(ipn),' /= ',nprox(ipn),']  me = ',mype
!SMS$IGNORE END
        stop
      end if
    end do
! verify that permedge is OK in the interior
    do ipn=1,nip
      do isn=1,nedge(ipn)
        if (permedge(isn,ipn) /= isn) then
!SMS$IGNORE BEGIN
          print *,'ERROR C:  permedge(',isn,',',ipn,') /= ',isn,')  [',permedge(isn,ipn),' /= ',isn,']  me = ',mype
!SMS$IGNORE END
          stop
        end if
      end do
    end do
!SMS$PARALLEL END
#endif
#undef DEBUG_HALO_COMP

!SMS$PARALLEL(dh, ipn) BEGIN
! fix prox for halo cells
! NOTE:  *not* owner-computes!
! TODO:  Move this loop into SMS_UnstructuredGrid along with the exchange of proxs above ??  
! TODO:  Would need to pass proxs into !SMS$UNSTRUCTURED_GRID too of course... 
    do ipn=1,nip
      do isn=1,nprox(ipn)
        ipnx = prox(isn,ipn)
        isnx = proxs(isn,ipn)
        if (in_halo(ipnx)) then
          ! ipnx is a halo cell
          ! point prox for halo cell at ipnx back to me (ipn)
          prox(isnx,ipnx) = ipn
          ! handle pointing prox for halo cell at ipnx to any halo cells that 
          ! are adjacent to itself (ipnx) and I (ipn)
          ! im1 and ip1 are my edges adjacent to cells that are also adjacent 
          ! to the halo cell at ipnx
          im1 = mod(isn-2+nprox(ipn),nprox(ipn)) + 1
          ip1 = mod(isn             ,nprox(ipn)) + 1
          ! nm1 and np1 are edges of halo cell at ipnx that are also adjacent 
          ! to me (ipn)
          nm1 = mod(isnx-2+nprox(ipnx),nprox(ipnx)) + 1
          np1 = mod(isnx              ,nprox(ipnx)) + 1
          ! note that nm1 and ip1 adjoin the same cell
          prox(nm1,ipnx) = prox(ip1,ipn)
          ! note that np1 and im1 adjoin the same cell
          prox(np1,ipnx) = prox(im1,ipn)
          ! note that prox(nm1,ipnx) and prox(np1,ipnx) are computed redundantly
        end if
      end do
    end do
!SMS$PARALLEL END

!SMS$IGNORE BEGIN
!print *,'DEBUG:  setting halo_comp values in first layer of dh__S1 and dh__E1 by brute force'
!print *,'DEBUG:  collapsed_halo_size = ',collapsed_halo_size
!print *,'DEBUG:  dh__NestLevel = ',dh__NestLevel
!print *,'DEBUG:  dh__S1(  1,0,1) = ',dh__S1(  1,0,1)
!print *,'DEBUG:  dh__S1(  1,1,1) = ',dh__S1(  1,1,1)
!print *,'DEBUG:  dh__E1(nip,0,1) = ',dh__E1(nip,0,1)
!print *,'DEBUG:  dh__E1(nip,1,1) = ',dh__E1(nip,1,1)
!SMS$IGNORE END

#define DEBUG_HALO_COMP
!TBH:  lots of error checks for debugging
#ifdef DEBUG_HALO_COMP
!SMS$PARALLEL(dh, ipn) BEGIN
! verify that in_halo is not screwed up in the interior
    do ipn=1,nip
      if (in_halo(ipn)) then
!SMS$IGNORE BEGIN
        print *,'ERROR X:  in_halo(',ipn,') = ',in_halo(ipn),'  but ipn is an interior cell!'
!SMS$IGNORE END
        stop
      end if
    end do
! verify that prox(isnx,ipnx) == ipn is still true in the interior
    do ipn=1,nip
      do isn=1,nprox(ipn)
        ipnx = prox(isn,ipn)
        isnx = proxs(isn,ipn)
        if (.not.in_halo(ipnx)) then  ! avoid halo points
          if (prox(isnx,ipnx) /= ipn) then
!SMS$IGNORE BEGIN
            print *,'ERROR X:  prox(',isnx,',',ipnx,') /= ',ipn,'  me = ',mype
!SMS$IGNORE END
            stop
          end if
        end if
      end do
    end do
!SMS$HALO_COMP(<1,1>) BEGIN
! verify that prox(isnx,ipnx) == ipn everywhere that we care
    do ipn=1,nip
      do isncount=1,nedge(ipn)
        isn = permedge(isncount,ipn)
        ipnx = prox(isn,ipn)
        isnx = proxs(isn,ipn)
        if (prox(isnx,ipnx) /= ipn) then
!SMS$IGNORE BEGIN
          print *,'ERROR X:  prox(',isnx,',',ipnx,') /= ',ipn,'  me = ',mype,'  in_halo(',ipn,') = ', &
                  in_halo(ipn),'  in_halo(',ipnx,') = ',in_halo(ipnx)
!SMS$IGNORE END
          stop
        end if
      end do
    end do
!SMS$HALO_COMP END
!SMS$PARALLEL END
#endif
#undef DEBUG_HALO_COMP

!SMS$PARALLEL(dh, ipn) BEGIN
    do ipn=1,nip
      corio(ipn) = omegx2*sin(lat(ipn))	! coriolis acceleration
      deg_lat(ipn) = raddeg*lat(ipn)	! latitude in degrees
      deg_lon(ipn) = raddeg*lon(ipn)	! longitude in degrees
      rarea(ipn) = 1./area(ipn)     	! reciprical of area
      call GetIpnGlobalMype (ipn, ipnGlobal, mype, DiagPrint)
      if(DiagPrint) then
        ipnDiagLocal = ipn
        ipnDiagPE    = mype
      end if
      actual(ipn) = ipn
    end do

    do ipn=1,nip
      do isn=1,nprox(ipn)
        rsideln(isn,ipn) = 1./sideln(isn,ipn)
      end do
    end do

!SMS$EXCHANGE(corio,deg_lat,deg_lon,rarea,rsideln,actual)

! ...............................................................
! Sec. 2. Load initial state
! ...............................................................

    call StartTimer(t0)
!SMS$SERIAL (<thetac,OUT> : default=ignore)  BEGIN
    unitno = getunit ()
    if (unitno < 0) then
      print*,'dyn_init: getunit failed for file=theta_coor.txt. Stopping'
      stop
    end if

    open (unitno, file='theta_coor.txt', form='formatted', action='read', iostat=ios)
    if (ios /= 0) then
      print*,'dyn_init: cannot open theta_coor.txt for reading. Stopping'
      stop
    end if
  
    read (unitno, *, iostat=ios) thetac
    if (ios /= 0) then
      print*,'dyn_init: error reading theta_coor.txt. Stopping'
      stop
    end if
  
!do ivl = 1, nvl
!  read (26,*) thetac(ivl)
!end do
    close (unitno)
!SMS$SERIAL END
    print '(a/(10f8.1))','thetac (deg K):', thetac
!SMS$SERIAL (<dpsig,OUT> : default=ignore)  BEGIN
! Re-use the same unit number that was just closed
    open (unitno, file='dpsig.txt', form='formatted', action='read', iostat=ios)
    if (ios /= 0) then
      print*,'dyn_init: cannot open dpsig.txt for reading. Stopping'
      stop
    end if

    read (unitno, *, iostat=ios) dpsig
    if (ios /= 0) then
      print*,'dyn_init: error reading dpsig.txt. Stopping'
      stop
    end if

!do ivl = 1, nvl
!  read(unitno,*) dpsig(ivl)
!end do
    close (unitno)
    call returnunit (unitno)
!SMS$SERIAL END
    print '(a/(10f8.1))','dpsig (Pa):',dpsig

    call GetJdate(yyyymmddhhmm,jdate)		! Julian date conversion
    hh = yyyymmddhhmm(9:10)

    if (EnKFAnl) then
! read from EnKF analysis
      call readenkfanal(nvl,FGFileName,EnKFFileName,FGFileNAmeSig,us3d,vs3d,dp3d,mp3d,pr3d,ex3d,ph3d,tr3d)
!    print*,'retured from readenkf'
!    DO i=1,nvl,6
!       print*,i,pr3d(i,400),ph3d(i,400),ex3d(i,400),tr3d(i,400,1)
!    END DO
    else if (.not. readrestart) then
      sanlFile = jdate // ".gfs.t" // hh // "z.sanl"
      print *,' get initial state from ',sanlFile

! get control info on input data
! "82" is a magic unit that might have to be byte-swapped

!SMS$SERIAL (<nvp,OUT> : default=ignore)  BEGIN
      lusig = getunit(82)
      if (lusig < 0) then
        print*,'dyn_init: getunit failed for unit=82. Stopping'
        stop
      end if

      call sigio_srohdc(lusig,sanlFile,head,data,iret)
      if (iret .ne. 0) then
        print '(a)','dyn_init: error reading '//sanlFile
        STOP
!TBH:  this code hangs because errexit does not call MPI_ABORT
!    call errmsg('dyn_init: error reading '//sanlFile)
!    call errexit(2)
      end if
      call returnunit (lusig)
      nvp = head%levs			! # of layers in input grid
!SMS$SERIAL END

      call ss2icos (nvp, sanlFile, us3d, vs3d, dp3d, mp3d, pr3d, ex3d, ph3d, tr3d, gfsltln_file)
    end if         ! enkfanl
!TBH:  ss2icos only sets tr3d(:,:,1:4)
!TBH:  so set remaining values to zero here
!TODO:  This may not be the correct approach.  Need to get 
!TODO:  Georg and Rainer together to discuss.  At the moment Georg 
!TODO:  believes that tr3d(:,:,5:ntra) are all zero anyway...  

    if (.not. readrestart) then
!SMS$IGNORE BEGIN
! put initialization inside SMS "ignore" so all values are set to 
! zero (including halo)
      if (ntra+ntrb > 4) then
        tr3d(:,:,5:ntra+ntrb) = 0.0
      end if
!SMS$IGNORE END

! Initialize post. If gribout is enabled, ensure there is a max of 1 write task.

!SMS$SERIAL BEGIN
  call post_init_slint (iret)
  if (iret /= 0) then
    write(6,*)'dyn_init: bad return from post_init_slint: stopping'
    stop
  endif
  if (.not. post_init_slint_called) then
    write(6,*)'dyn_init: logic error: cannot test gribout until post_init_slint has been called'
    stop
  endif
!SMS$SERIAL END

      call findmxmn2(ph3d,nvlp1,nip,1,'surf.geopot.')
      call IncrementTimer(t0,t1)
      print"(' DYNAMICS INPUT time:',F10.0)",t1

      do ipn=1,nip
        do ivl=1,nvl  					! vertical loop
          do i=1,ntra+ntrb
            trdp(ivl,ipn,i) = tr3d(ivl,ipn,i)*dp3d(ivl,ipn)
          end do
        end do
      end do
 
      its = itsStart - 1

      print *,'dyn_init calling hybgen ...'

      call hybgen(its,         &
           thetac,                 & ! target pot.temperature
           us3d,vs3d,tr3d,         & ! zonal, meridional wind, mass field tracers
           sdot,ex3d,dp3d,pr3d     ) ! intfc displ., exner, lyr thknss, pressure
!    print*,'retured from hybgen'
!    DO i=1,nvl,6
!       print*,i,pr3d(i,400),ph3d(i,400),ex3d(i,400),tr3d(i,400,1)
!    END DO

      do ipn=1,nip
        pr3d(1,ipn) = p1000*(ex3d(1,ipn)/cp)**(cp/rd)
        do ivl=1,nvl  !  vertical loop
          pr3d(ivl+1,ipn) = p1000*(ex3d(ivl+1,ipn)/cp)**(cp/rd)
          dp3d(ivl  ,ipn) = pr3d(ivl,ipn) - pr3d(ivl+1,ipn)
        end do
      end do           ! horizontal loop

      psrf  = 0.
      ptdcy = 0.

      print *,'dyn_init calling hystat ...'
      call hystat(its,              &
           ph3d,ex3d,mp3d,            &  ! geopotential, exner fct., montg.pot.
           dp3d,tr3d,trdp,            &  ! layer thknss, tracer, tracer x thknss
           psrf,ptdcy     )              ! srf.pressure, srf.prs. tendency
      
! ...............................................................
! Sec. 3. Initialize misc. variables
! ...............................................................

! Initialize tendency arrays
      u_tdcy   = 0.		! u tendency
      v_tdcy   = 0.		! v tendency
      dp_tdcy  = 0.		! dp tendency
      dpl_tdcy = 0.		! dp tendency, low order
      trc_tdcy = 0.		! tracer tendency
      trl_tdcy = 0.		! tracer tendency, low order
      massfx   = 0.		! mass flux (3 time levels)
      u_tdcy_phy   = 0.         ! physics u tendency
      v_tdcy_phy   = 0.         ! physics v tendency
      trc_tdcy_phy = 0.         ! physics tracer tendency
!SMS$PARALLEL END

! initial Adams-Bashforth indices
      nf  = 0	! "new field" index 
      of  = 2	! "old field" index
      vof = 1	! "very old field" index

      ws3d = 0.
      tk3d = 0.0
      pw2d = 0.     ! precipitable water
      rn2d = 0.     ! accumulated precipitation/rainfall
      rc2d = 0.0
      flxlwtoa2d = 0.0

      if (ntrb > 0) then
        call transp0(its,cumufx,dp3d,dpinit)	! initialize class B tracer transport
      end if
    end if    ! .not. readrestart

    if (readrestart) then
!SMS$SERIAL BEGIN
      call post_init_slint (iret)
      if (iret /= 0) then
        write(6,*)'dyn_init: bad return from post_init_slint: stopping'
        stop
      endif
      if (.not. post_init_slint_called) then
        write(6,*)'dyn_init: logic error: cannot test gribout until post_init_slint has been called'
       stop
      endif
!SMS$SERIAL END
    end if ! readrestart: initialize slint

    vor  = 0.
    rh3d = 0.

    deallocate(in_halo)

    if (fixedgridorder) call icosio_set_inv_perm(inv_perm)
    if (use_write_tasks) call icosio_prep

! --- exercising stencl routine
    call stencl(deg_lat,1,1.,'latitude')
    call stencl(deg_lon,1,1.,'longitude')

! All the work_edg stuff need only be done in the initial run
    if (readrestart) then
      return
    end if
!-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
! --- special diagnostics to figure out what goes on at the poles

!    if (PrintIpnDiag == 1 .or. PrintIpnDiag == nip) then  
    allocate (work_edg(1,npp,nip))
    do idx=1,2
!SMS$PARALLEL(dh, ipn) BEGIN
      do ipn=1,nip
        do isn=1,nprox(ipn)
          work_edg(1,isn,ipn) = sidevec_e(idx,isn,ipn)
!!SMS$ignore begin
!      if (ipn==PrintIpnDiag) write (*,'(2i7,a,f8.1))')		&
!       ipn,prox(isn,ipn),string,work_edg(1,isn,ipn)
!!SMS$ignore end
        end do
      end do
!SMS$EXCHANGE(work_edg)
!SMS$PARALLEL END
      write (string,'(a,i1,a)') ' sidevec_e(',idx,')'
      call stenedg(deg_lat,work_edg,1,string)
!SMS$PARALLEL(dh, ipn) BEGIN
      do ipn=1,nip
        do isn=1,nprox(ipn)
          work_edg(1,isn,ipn) = sidevec_c(idx,isn,ipn)
!!SMS$ignore begin
!      if (ipn==PrintIpnDiag) write (*,'(2i7,a,f8.1))')		&
!       ipn,prox(isn,ipn),string,work_edg(1,isn,ipn)
!!SMS$ignore end
        end do
      end do
!SMS$EXCHANGE(work_edg)
!SMS$PARALLEL END
      write (string,'(a,i1,a)') ' sidevec_c(',idx,')'
      call stenedg(deg_lat,work_edg,1,string)
    end do

    do idx=1,4
!SMS$PARALLEL(dh, ipn) BEGIN
      do ipn=1,nip
        do isn=1,nprox(ipn)
          work_edg(1,isn,ipn) = sn(idx,isn,ipn)
        end do
!!SMS$ignore begin
!      if (ipn==PrintIpnDiag) write (*,'(2i7,a,f8.1))')		&
!       ipn,prox(isn,ipn),string,work_edg(1,isn,ipn)
!!SMS$ignore end
      end do
!SMS$EXCHANGE(work_edg)
!SMS$PARALLEL END
      write (string,'(a,i1,a)') ' sn(',idx,')'
      call stenedg(deg_lat,work_edg,1,string)
!SMS$PARALLEL(dh, ipn) BEGIN
      do ipn=1,nip
        do isn=1,nprox(ipn)
          work_edg(1,isn,ipn) = cs(idx,isn,ipn)
        end do
!!SMS$ignore begin
!      if (ipn==PrintIpnDiag) write (*,'(2i7,a,f8.1))')		&
!       ipn,prox(isn,ipn),string,work_edg(1,isn,ipn)
!!SMS$ignore end
      end do
!SMS$EXCHANGE(work_edg)
!SMS$PARALLEL END
      write (string,'(a,i1,a)') ' cs(',idx,')'
      call stenedg(deg_lat,work_edg,1,string)
    end do
    deallocate (work_edg)
!-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

    print *,'... exiting dyn_init'
    return

90  write(6,*)'dyn_init: error reading a file'
    call flush(6)
    stop
  end subroutine dyn_init
end module module_fim_dyn_init
