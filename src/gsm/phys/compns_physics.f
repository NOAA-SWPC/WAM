!-----------------------------------------------------------------------
      subroutine compns_physics(deltim, iret,     ntrac,   nxpt,   nypt,
     &                          jintmx, jcap,     levs,    levr,   lonr,
     &                          latr,   ntoz,     ntcw,    ncld,   ntke,
     &                          lsoil,  nmtvr,    num_p3d, num_p2d,
     &                          npdf3d, ncnvcld3d,nshoc_3d, nshoc_2d,
     &                          ntot3d, ntot2d,
     &                          thermodyn_id, sfcpress_id,
     &                          nlunit,   me, gfs_phy_namelist)
!
!$$$  Subprogram Documentation Block
!
! Subprogram:  compns     Check and compute namelist frequencies
!   Prgmmr: Iredell       Org: NP23          Date: 1999-01-26
!
! Abstract: This subprogram checks global spectral model namelist
!           frequencies in hour units for validity.  If they are valid,
!           then the frequencies are computed in timestep units.
!           The following rules are applied:
!             1. the timestep must be positive;
!             2. the output frequency must be positive and
!                a multiple of the timestep to within tolerance;
!             3. the shortwave frequency must be positive and
!                a multiple of the timestep to within tolerance;
!             4. the longwave frequency must be positive and
!                a multiple of the timestep to within tolerance and
!                a multiple of the shortwave frequency;
!             5. the zeroing frequency must be positive and
!                a multiple of the timestep to within tolerance and
!                a multiple of the output frequency;
!             6. the restart frequency must be positive and
!                a multiple of the timestep to within tolerance and
!                a multiple of the longwave frequency and
!                a multiple of the zeroing frequency;
!             7. the initialization window must be non-negative and
!                a multiple of the timestep to within tolerance and
!                a multiple of the longwave frequency and
!                no longer than the restart frequency;
!             8. the cycling frequency must be non-negative and
!                a multiple of the timestep to within tolerance and
!                a multiple of the longwave frequency.
!
! Program History Log:
!   1999-01-26  Iredell
!   2009-05-04  Moorthi
!   2009-10-12  Sarah Lu, add grid_aldata (default to F)
!   2010-01-12  Sarah Lu, add fdaer (default to 0)
!   2010-08-03  Jun Wang, add fhdfi (default to 0)
!   2010-jul/augMoorthi  added many physics options including nst
!   2012-04-06  Henry Juang, add idea 
!   2012-10-18  S. Moorthi add use_ufo, dtphys,  and remove hdif_fac 
!   2010-08-03  Jun Wang, add sfcio option
!   2013-04-09  Jun Wang, set ivegsrc and cu_physics
!   2013-10-30  Xingren Wu, add a2oi_out and ngrid_a2oi for A/O/I coupling
!   2013-01     Y-T Hou    - 1. changed radiation calling time control
!                  variables, fhswr/fhlwr, from previous in unit of hours
!                  to in unit of seconds. thus radiation callinginterval
!                  can be made less than the old one-hour limit.
!                            2. added aerosol model selection flag, iaer_mdl
!                  for choices of opac-clim, gocart-clim, and gocart-prog...
!  2013-08-14 S. Moorthi - Merging new radiation and stochastic physics to NEMS
!  2014-03-18 Sarah Lu - Remove iaer_mdl
!  2014-04-19 Xingren Wu, add cplflx for A/O/I coupling
!  2014-04-28 Jun Wang - add cgwf,prslrd0
!  2014-05-02 Philip Pegion - added stocast physics (shum,skeb,sppt,vc)
!  2014-06-27 S. Moorthi - add semilag and turn off shuff_lats_r
!  2014-08-01 s. moorthi - add fixtrc - tracer fixer option
!  2015-05-28 Xingren Wu,  change ngrid_a2oi=50
!  2015-08-25 Xu Li      - change nst_fcst to be nstf_name
!  2015-10-02 s moorthi  - add grboro logical - input orog is binary, if false
!  2015-10-19 Weiyu Yang - add the f107_kp_size, f107_kp_interval.
!                          number of inputted f10.7 and kp data and the time
!                          interval of them.
!  2016-03-03 J. Han  - add cnvcld logical ncnvcld3d integer
!                       for convective cloudiness enhancement
!  2016-03-04 J. Han  - change newsas & sashal to imfdeepcnv 
!                        & imfshalcnv, respectively
!  2016-03-07 Weiyu Yang - add the f107_kp_skip_size, the skip size
!                          before the forecast starting.
!
!  2017-02-21  VAY       call wam_control_default & call idea_wamcontrol_init(me)
!  2017-03-12  VAY       commentary changes
!  2017-09-10  Weiyu Yang - Add IPE back coupling to WAM code.
! Usage:    call compns(deltim,
!    &                  fhout,fhswr,fhlwr,fhzer,fhres,fhcyc,
!    &                  nsout,nsswr,nslwr,nszer,nsres,nscyc,
!    &                  iret)
!   Input Arguments:
!     tol      - real error tolerance allowed for input frequencies
!                (e.g. 0.01 for 1% of timestep maximum error allowed)
!     deltim   - real timestep in seconds
!     fhout    - real output frequency in hours
!     fhswr    - real shortwave frequency in hours
!     fhlwr    - real longwave frequency in hours
!     fhzer    - real zeroing frequency in hours
!     fhres    - real restart frequency in hours
!     fhcyc    - real cycling frequency in hours
!   Output Arguments:
!     nsout    - integer output frequency in timesteps
!     nsswr    - integer shortwave frequency in timesteps
!     nslwr    - integer longwave frequency in timesteps
!     nszer    - integer zeroing frequency in timesteps
!     nsres    - integer restart frequency in timesteps
!     nscyc    - integer cycling frequency in timesteps
!     iret     - integer return code (0 if successful or
!                between 1 and 8 for which rule above was broken)
!     LDIAG3D  - switch for 3D diagnostic- (default = false)
!
! Attributes:
!   Language: Fortran 90
!
!$$$

      
      use namelist_physics_def
      use module_IPE_to_WAM, ONLY: ipe_to_wam_coupling
!
!vay-2017
!
      use namelist_wamphysics_def
!cmy mpi_def holds liope
      use mpi_def, only : liope
      use wam_ifp_mod, only: params

      implicit none

      real tol
 
      character (len=*), intent(in) :: gfs_phy_namelist
      integer, intent(in)           :: me, nlunit
      real,intent(inout)            :: deltim
      integer,intent(out)           :: iret
      integer ntrac,nxpt,nypt,jintmx,jcap,levs,lonr,latr
!     integer ntrac,jcap,levs,lonr,latr
      integer levr
      integer ntoz,ntcw,ncld,lsoil,nmtvr,num_p3d,num_p2d,member_num
     &,       npdf3d, ncnvcld3d, nshoc_3d, nshoc_2d, ntke 
     &,       ntot3d, ntot2d
      integer thermodyn_id, sfcpress_id, ifp_realtime_interval
      real    tfiltc
      logical lgoc3d

!sela - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     if output (fhout) more frequently than zeroing ,get partial rains
 
      namelist /nam_phy/FHMAX,FHOUT,FHRES,FHZER,FHSEG,FHROT,IGEN,
     & NGPTC,fhswr,fhlwr,fhcyc,ras,LGOC3D,FHGOC3D,LDIAG3D,reduced_grid,
     & shuff_lats_r,thermodyn_id,sfcpress_id,fhdfi,semilag,gg_tracers,
     & pre_rad,hybrid,gen_coord_hybrid,random_clds,liope,
     & ntrac,nxpt,nypt,jintmx,jcap,levs,lonr,latr,levr,
     & ntoz,ntcw,ncld,ntke,lsoil,nmtvr,zhao_mic,nsout,lsm,tfiltc,
     & isol, ico2, ialb, iems, iaer, iovr_sw, iovr_lw,ictm,
     & isubc_lw, isubc_sw, fdaer, lsidea, weimer_model,
     & ipe_to_wam_coupling, ifp_realtime_interval,
     & ncw, crtrh,old_monin,flgmin,cnvgwd,cgwf,prslrd0,ral_ts,fixtrc,
!    & ncw, crtrh,old_monin,flgmin,gfsio_in,gfsio_out,cnvgwd,
     & ccwf,shal_cnv,imfshalcnv,imfdeepcnv,
     & crick_proof,ccnorm,ctei_rm,mom4ice,
     & norad_precip,mstrat,trans_trac,dlqf,moist_adj,cscnv,
     & do_shoc,shocaftcnv,shoc_cld,
     & nstf_name,nst_anl,lsea,cal_pre,psautco,prautco,evpco,wminco,
     & pdfcld,shcnvcw,sup,cnvcld,redrag,hybedmf,dspheat,
     & fhout_hf,fhmax_hf,cdmbgwd,bkgd_vdif_m,bkgd_vdif_h,
     & grid_aldata,bkgd_vdif_s,use_ufo,dtphys,nemsio_in,nemsio_out,
     & sfcio_out,climate,a2oi_out,cplflx,ngrid_a2oi,griboro,
     & sppt,skeb,shum,vc,vcamp,
! iau parameters
     & iau,iaufiles_fg,iaufiles_anl,iaufhrs,iau_delthrs,iaulnp
!
      shuff_lats_r   = .false.
!     reshuff_lats_r = .false.
!
      fhmax    = 0
      fhout    = 0
      fhzer    = 0
      fhseg    = 0
      fhrot    = 0
      fhout_hf = 1
      fhmax_hf = 0
      dtphys   = 3600
      igen     = 0
      fhswr    = 0
      fhlwr    = 0
      fhcyc    = 0
      fhdfi    = 0
      tfiltc   = 0.85
      ccwf     = 1.0
      dlqf     = 0.0
      cgwf     = 0.1      ! cloud top fraction for convective gwd scheme
      ctei_rm  = 10.0
      NGPTC    = lonr
      ntrac    = 3
!
      bkgd_vdif_m = 1.0
      bkgd_vdif_h = 1.0
      bkgd_vdif_s = 1.0
      prslrd0     = 200.  ! pressure(pa) above which Raleigh damping applied
      ral_ts      = 10.0  ! time scale in days for Raleigh damping
!
      ras              = .false.
      cscnv            = .false.
      do_shoc          = .false.
      shocaftcnv       = .false.
      shoc_cld         = .false.
      zhao_mic         = .true.
      LDIAG3D          = .false.
      LGOC3D           = .false.
      fhgoc3d          =  0.0
      shal_cnv         = .true.
      imfshalcnv       = 1     
      crick_proof      = .false.
      ccnorm           = .false.
      imfdeepcnv       = 1     
      redrag           = .false.
      hybedmf          = .false.
      dspheat          = .false.
      norad_precip     = .false.   ! This is effective only for Ferrier/Moorthi
      mom4ice          = .false.   ! True when coupled to MOM4 OM
      mstrat           = .false.
      trans_trac       = .true.    ! This is effective only for RAS
      moist_adj        = .false.   ! Must be true to turn on moist convective
      cal_pre          = .false.   ! true for huiya's precip type algorithm
!
      reduced_grid     = .true.
!
      pre_rad          = .false.
      hybrid           = .true.
      gen_coord_hybrid = .false.                                     !hmhj
      random_clds      = .false.
      liope            = .false.
!
      old_monin        = .false.
      cnvgwd           = .false.
      use_ufo          = .false.

      fixtrc(:)        = .false.
      semilag          = .false.
      gg_tracers       = .false.

! idea add
      lsidea           = .false.
      weimer_model     = 'epot'
      ipe_to_wam_coupling = .false.
      ifp_realtime_interval = -1

! Add: a2oi_out, cplflx & ngrid_a2oi
      a2oi_out         = .false.     ! default not writing out Atm fields for Ocn/Ice
      cplflx           = .false.     ! default no cplflx collection
      ngrid_a2oi       = 50          ! # of Atm fields for the coupling use
!
      thermodyn_id     = 1
      sfcpress_id      = 1
!
!     ncw(1)           = 75
      ncw(1)           = 50
      ncw(2)           = 150
      crtrh(:)         = 0.85
      flgmin(:)        = 0.20
!
      psautco(:)       = 4.0E-4    ! Zhao scheme default opr value
      prautco(:)       = 1.0E-4    ! Zhao scheme default opr value
      evpco            = 2.0E-5    ! Zaho scheme evaporation coefficient
      wminco(:)        = 1.0E-5    ! Zhao scheme default water and ice floor value
      cdmbgwd(:)       = 1.0       ! Mtn Blking and GWD tuning factors

      sup              = 1.1
      pdfcld           = .false.
      shcnvcw          = .false.
      cnvcld           = .false.
      npdf3d           = 0
      ncnvcld3d        = 0

                                   ! For NST model
      nstf_name(:)     = 0         ! initialized as zeros
      nst_anl          = .false.   ! if true, using nst analysis over lakes
      lsea             = 0
!
      nemsio_in        = .true.
      nemsio_out       = .true.
      sfcio_out        = .false.
      climate          = .false.
      griboro          = .true.

!
      ! can specify up to 5 values for the stochastic physics parameters 
      ! (each is an array of length 5)
      sppt             = -999.  ! stochastic physics tendency amplitude
      shum             = -999.  ! stochastic boundary layer spf hum amp   
      skeb             = -999.  ! stochastic KE backscatter amplitude
      vcamp            = -999.  ! stochastic vorticity confinment amp
      vc              = 0.      ! deterministic vorticity confinement parameter.
! !  additions for stochastic perturbations... (details are in dynamics)
      do_sppt = .false.
      do_shum = .false.
      do_skeb = .false.
      do_vc   = .false.
! gaussian or power law variance spectrum for skeb (0: gaussian, 1:
! power law). If power law, skeb_lscale interpreted as a power not a
! length scale.
! additions for iau
      iau              = .false.
      iaufhrs          = -1
      iaufiles_fg      = ''
      iaufiles_anl     = ''
!
      nsout   = 0
      nsout_hf = 0
      lsm     = 1         ! NOAH LSM is the default when lsm=1
      levr    = 0
!           Default values for some radiation controls
!           ------------------------------------------
      isol     = 0        ! use prescribed solar constant
      ico2     = 0        ! prescribed global mean value (old opernl)
      ialb     = 0        ! use climatology alb, based on sfc type
!     ialb     = 1        ! use modis based alb
      iems     = 0        ! use fixed value of 1.0
      iaer     = 1        ! default aerosol effect in sw only
      iovr_sw  = 1        ! sw: max-random overlap clouds
      iovr_lw  = 1        ! lw: max-random overlap clouds
      ictm     = 1        ! ictm=0 => use data at initial cond time, if not
                          ! available, use latest, no extrapolation.
                          ! ictm=1 => use data at the forecast time, if
                          ! not
                          ! available, use latest and extrapolation.
                          ! ictm=yyyy0 => use yyyy data for the forecast
                          ! time,
                          ! no further data extrapolation.
                          ! ictm=yyyy1 = > use yyyy data for the fcst.
                          ! if needed, do extrapolation to match the
                          ! fcst time.
                          ! ictm=-1 => use user provided external data
                          ! for
                          ! the fcst time, no extrapolation.
                          ! ictm=-2 => same as ictm=0, but add seasonal
                          ! cycle
                          ! from climatology. no extrapolation.

      isubc_sw = 0        ! sw clouds without sub-grid approximation
      isubc_lw = 0        ! lw clouds without sub-grid approximation
                          ! =1 => sub-grid cloud with prescribed seeds
                          ! =2 => sub-grid cloud with randomly generated
                          ! seeds
!
      ivegsrc = 2         ! ivegsrc=0:USGS, 1:IGBP, 2:WMD
!
! The copy/pointer option (Sarah Lu) - 3D fields are allocated only in DYN;
! pointer is used to associate 3D fields in PHY back to DYN 

      grid_aldata = .false.
!
! The relaxation time in days to gocart forecast (Sarah Lu)
! default is 0 (use clim/anal fields)

      fdaer = 0.0
!
      if (me == 0) print *,' nlunit=',nlunit,' gfs_phy_namelist=',
     &                      gfs_phy_namelist

!$$$      read(5,nam_phy)
      open(unit=nlunit,file=gfs_phy_namelist)
      rewind (nlunit)
!     read(nlunit,nam_phy,err=999)
      read(nlunit,nam_phy)
!
!VAY-2017 ============================= lsidea
!   add copies of wam_control_in to RT-scripts
!
      if (lsidea) then
!
      call wam_control_default
      params % ifp_realtime_interval = ifp_realtime_interval
!
!       read_in SW-drivers WAM switches from "nam_wam_control"
!
      call idea_wamcontrol_init(me)

      ENDIF   !lsidea
!============================= lsidea

!
      if (cal_pre) then
        random_clds = .true.
        if (me == 0) print *,' cal_pre=',cal_pre,' random_clds=',
     &                       random_clds
      endif
!
      if (me == 0) print *,' fhmax=',fhmax,' nstf_name =',nstf_name
     &,'lsea =',lsea,' fhout_hf=',fhout_hf,' fhmax_hf=',fhmax_hf

      if (do_shoc) then
        nshoc_3d = 3
        nshoc_2d = 0
        shal_cnv = .false.
        imfshalcnv  = -1  
        hybedmf  = .false.
        if (me == 0)
     &  print *,' Simplified Higher Order Closure Model used for'
     &,         ' Boundary layer and Shallow Convection'
     &,         '  nshoc_3d=',nshoc_3d,'nshoc_2d=',nshoc_2d
     &,         ' ntke=',ntke
      else
        nshoc_3d = 0
        nshoc_2d = 0
      endif
!
!cup_physics category: 1:bmj, 2:kf, 3:gd, 4:sas, 5:ras
      if (ras) then
        cu_physics=5
      else
!sas
        cu_physics=4
      endif
!
      if (me == 0) then
        write(6,nam_phy)
        if (lsm == 1) then
          print *,' NOAH Land Surface Model used'
        elseif (lsm == 0) then
          print *,' OSU Land Surface Model used'
        else
          print *,' Unsupported LSM type - job aborted'
     &,                        ' - lsm=',lsm
          call mpi_quit(2222)
        endif
!
        print *,'in phy init, nemsio_in=',nemsio_in,'nemsio_out=',
     &           nemsio_out,'sfcio_out=',sfcio_out
!
        print *,' in compns_physics use_ufo=',use_ufo,
     &                            ' nst_anl=',nst_anl
        if (.not. cscnv) then 
          if (ras) then
            print *,' RAS Convection scheme used with ccwf=',ccwf
          else
            if (imfdeepcnv==0) then
              print *,' old SAS Convection scheme before July 2010 used'
            elseif(imfdeepcnv==1) then
              print *,' July 2010 version of SAS conv scheme used'
            elseif(imfdeepcnv==2) then
            print *,' scale- & aerosol-aware mass-flux deep conv scheme'
            endif
          endif
        else
          print *,'Chikira-Sugiyama Convection scheme used'
        endif
        print *,' The time filter coefficient tfiltc=',tfiltc,
     &   'cu_physics=',cu_physics
        if (.not. old_monin .and. .not. do_shoc) 
     &    print *,' New PBL scheme used'
        if (.not. shal_cnv) then
          imfshalcnv = -1      
          print *,' No shallow convection used'
        else
          if (imfshalcnv>0) print *,
     &     ' New Massflux based shallow convection used'
        endif
        if (cnvgwd) print *,' Convective GWD parameterization used'
          if (crick_proof) print *,' CRICK-Proof cloud water used in'
     &,                            ' radiation '
          if (ccnorm) print *,' Cloud condensate normalized by cloud'
     &,                       ' cover for radiation'
      endif
      if (.not. semilag) then
        fixtrc(:)  = .false.
        gg_tracers = .false.
      endif
!
      if (levr == 0) then
        levr = levs
      endif
      if (me == 0) then
        print *,' Radiative heating calculated at',levr, ' layers'
        if (iovr_sw == 0) then
          print *,' random cloud overlap for Shortwave IOVR_SW='
     &,           iovr_sw
        else
          print *,' max-random cloud overlap for Shortwave IOVR_SW='
     &,             iovr_sw
        endif
        if (iovr_lw == 0) then
          print *,' random cloud overlap for Longwave IOVR_LW='
     &,           iovr_lw
        else
          print *,' max-random cloud overlap for Longwave IOVR_LW='
     &,             iovr_lw
        endif
        if (isubc_sw == 0) then
          print *,' no sub-grid cloud for Shortwave ISUBC_SW='
     &,           isubc_sw
        else
          print *,' sub-grid cloud for Shortwave ISUBC_SW='
     &,           isubc_sw
        endif
        if (isubc_lw == 0) then
          print *,' no sub-grid cloud for Longwave ISUBC_LW='
     &,           isubc_lw
        else
          print *,' sub-grid cloud for Longwave ISUBC_LW='
     &,           isubc_lw
        endif
      endif
!
      if (zhao_mic) then        ! default setup for Zhao Microphysics
        num_p3d = 4
        num_p2d = 3
        if (pdfcld) then
          npdf3d = 3
        else
          shcnvcw = .false.
        endif
        if (me == 0) print *,' Using Zhao/Carr/Sundqvist  Microphysics'
      else                      ! Brad Ferrier's Microphysics
        num_p3d = 3
        num_p2d = 1
        pdfcld  = .false.
        shcnvcw = .false.

        if (me .eq. 0) print *,' Using Ferrier Microphysics'
      endif
!
      if(ras .or. cscnv) cnvcld = .false.
      if(do_shoc .or. pdfcld) cnvcld = .false.
      if(cnvcld) ncnvcld3d = 1
!
      ntot3d = num_p3d + nshoc_3d + npdf3d + ncnvcld3d
      ntot2d = num_p2d + nshoc_2d
      if (me == 0) print *,' num_p3d=',num_p3d,' num_p2d=',num_p2d,
     &                     ' crtrh=',crtrh,' npdf3d=',npdf3d,
     &                     ' pdfcld=',pdfcld,' shcnvcw=',shcnvcw,
     &                     ' cnvcld=',cnvcld,' ncnvcld3d=',ncnvcld3d,
     &                     ' do_shoc=',do_shoc,' nshoc3d=',nshoc_3d,
     &                     ' nshoc_2d=',nshoc_2d,' shoc_cld=',shoc_cld
     &,                    ' ntot3d=',ntot3d,' ntot2d=',ntot2d
     &,                    'shocaftcnv=',shocaftcnv
!
!sela - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      tol = 0.01
!  Compute nsout and check rule 2.
      if(nsout > 0) fhout = nsout*deltim/3600.
      nsout=nint(fhout*3600./deltim)
      if(nsout <= 0.or.abs(nsout-fhout*3600./deltim) > tol) then
        iret = 2
        return
      endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Compute nsout_hf and check rule 21.
!     if(nsout_hf.gt.0) fhout=nsout_hf*deltim/3600.
      nsout_hf = nint(fhout_hf*3600./deltim)
      if(nsout_hf <= 0.or.abs(nsout_hf-fhout_hf*3600./deltim)>tol) then
        iret = 9
        return
      endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Compute nsswr and check rule 3.
      nsswr=nint(fhswr/deltim)
      if(nsswr <= 0.or.abs(nsswr-fhswr/deltim).gt.tol) then
        iret=3
        return
      endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Compute nslwr and check rule 4.
      nslwr = nint(fhlwr/deltim)
      if(nslwr <= 0.or.abs(nslwr-fhlwr/deltim) > tol .or.
     &   mod(nslwr,nsswr).ne.0) then
        iret = 4
        return
      endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Compute nszer and check rule 5.
      nszer = nint(fhzer*3600./deltim)
      if(nszer <= 0.or.abs(nszer-fhzer*3600./deltim).gt.tol.or.
     &   mod(nszer,nsout).ne.0) then
        iret=5
        return
      endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Compute nsres and check rule 6.
      nsres=nint(fhres*3600./deltim)
      if(nsres <= 0.or.abs(nsres-fhres*3600./deltim).gt.tol.or.
     &   mod(nsres,nslwr).ne.0.or.mod(nsres,nszer).ne.0) then
        iret = 6
        return
      endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Compute ndfi and check rule 7.
      if(fhdfi == 0.) then
        ndfi = 0
        ldfi = .false.
      else
        ndfi = nint(2*fhdfi*3600./deltim)
        ldfi = .true.
        if(ndfi <= 0.or.abs(ndfi-2*fhdfi*3600./deltim) > tol .or.
     &     ndfi > nsres) then
           print *,'ndfi=',ndfi,'is not equal to',2*fhdfi*3600./deltim
          iret = 7
          return
        endif
      endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Compute nscyc and check rule 8.
      if(fhcyc == 0.) then
        nscyc = 0
      else
        nscyc = nint(fhcyc*3600./deltim)
        if(nscyc <= 0.or.abs(nscyc-fhcyc*3600./deltim) > tol .or.
     &     mod(nscyc,nslwr) /= 0) then
          iret=8
          return
        endif
      endif
!!
      IF (NGPTC > LONR) THEN
         NGPTC = LONR
         WRITE(0,*) "NGPTC IS TOO BIG, RESET NGPTC TO LONR",NGPTC
      ENDIF
      IF (ME == 0)   WRITE(0,*) "NGPTC IS SET TO NGPTC :",NGPTC
      IF (sppt(1) > 0 ) THEN
        do_sppt=.true.
      ENDIF
      IF (shum(1) > 0 ) THEN
        do_shum=.true.
      ENDIF
      IF (skeb(1) > 0 ) THEN
        do_skeb=.true.
      ENDIF
      IF (vc > tiny(vc) .or. vcamp(1) > 0 ) THEN
        do_vc=.true.
      ENDIF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  All checks are successful.
      iret = 0
      return
  999 print *,' error reading  namelist - execution terminated by user'
      call mpi_quit(999)
!
      end
