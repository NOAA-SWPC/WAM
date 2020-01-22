!-----------------------------------------------------------------------
      subroutine compns_dynamics (deltim,iret,
     &               ntrac,nxpt,nypt,jintmx,jcap,
     &               levs,levr,lonf,latg, ntoz,
     &               ntcw,ncld, ntke, spectral_loop, me,
     &               thermodyn_id,sfcpress_id,
     &               nlunit, gfs_dyn_namelist,ndfi)
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
!             9. the difgital filter  must be non-negative and
!
! Program History Log:
!   1999-01-26  Iredell
!   2007-02-01  H.-M. H. Juang modify to be used for dynamics only
!   2009-11-09  Jun Wang       added ndfi 
!   2010-09-08  Jun Wang       change gfsio to nemsio
!   2011-02-11  Henry Juang    add codes to fit mass_dp and ndslfv
!   2011-02-28  Sarah Lu       add thermodyn_id and sfcpress_id
!   2012-04-06  Henry Juang    add idea for lsidea
!   2012-10-05  Jun Wang       add sigio_out
!   2012-10-18  S. Morthi      add hdif_fac, hdif_fac2, slrd0
!   2013-02-02  Henry Juang    revised reduced grid and add x number
!   2013-04-02  Jun Wang       add dfilevs for digital filer
!   2014-05-02  Phil Pegion    add stochastic physics parameterization
!   2014-07-17  S  Moorthi     updating to latest semi-lagrangian
!   2014-07-21  S  Moorthi     removed num_reduce (lonsperlat read from a file now)
!                              added cdamp and k2o
!  2016-03-07   Weiyu Yang - add the wam_ipe_coupling, and height_dependent_g:
!                            the WAM IPE model coupling flag and the flag that
!                            flag for using the height dependent g in that coupling.
!  2017-11-06   Weiyu Yang -
!               i)   add wam_ipe_cpl_rst_input for WAM-IPE coupling restart run,
!               ii)  add wam_ipe_cpl_rst_output for WAM-IPE coupling restart run,
!               iii) add NC_output and DELOUT_NC for outputing the
!                    NetCDF diagnostic files.
!
! Usage:    call compns(deltim,
!    &                  fhout,fhres,
!    &                  nsout,nsres,
!    &                  iret)
!   Input Arguments:
!     tol      - real error tolerance allowed for input frequencies
!                (e.g. 0.01 for 1% of timestep maximum error allowed)
!     deltim   - real timestep in seconds
!     fhout    - real output frequency in hours
!     fhres    - real restart frequency in hours
!   Output Arguments:
!     nsout    - integer output frequency in timesteps
!     nsskeb   - integer skeb update frequency in timesteps
!     nsres    - integer restart frequency in timesteps
!     iret     - integer return code (0 if successful or
!                between 1 and 8 for which rule above was broken)
!
! Attributes:
!   Language: Fortran 90
!
!$$$

      
      use namelist_dynamics_def
      use pmgrid         , only : wgt_slg, quamon
      use gfs_dyn_mpi_def, only : liope
      use layout_grid_tracers , only : xhalo, yhalo
      implicit none

      real tol
 
      character (len=*), intent(in) :: gfs_dyn_namelist
      integer, intent(in)           :: me, nlunit
      real,intent(inout)            :: deltim
      integer,intent(out)           :: iret
      integer ntrac,nxpt,nypt,jintmx,jcap,levs,lonf,latg
      integer levr,lsoil,nmtvr,lonr,latr,ndfi,k
      integer ntoz,ntcw,ncld,ntke,spectral_loop,member_num
      integer thermodyn_id, sfcpress_id
      real    tfiltc,tem

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
      namelist /nam_dyn/FHMAX,FHOUT,FHRES,FHROT,FHDFI,IGEN,
     & NGPTC,shuff_lats_a,reshuff_lats_a,
     & nxpt,nypt,jintmx,jcap,levs,lonf,latg,levr,
     & ntrac,ntoz,ntcw,ncld,ntke,nsout,tfiltc,
     & nemsio_in,nemsio_out,liope,ref_temp,ref_pres,lsidea,
     & explicit,hybrid,gen_coord_hybrid,process_split, 
     & spectral_loop,ndslfv,mass_dp,semi_implicit_temp_profile,
     & reduced_grid,hdif_fac,hdif_fac2,slrd0,cdamp,k2o,
     & thermodyn_id, sfcpress_id,gg_tracers,phigs_d,
     & zflxtvd,sigio_out,ldfi_spect,redgg_a,
     & sl_epsln,settls_dep3ds,settls_dep3dg,
     & herm_x,herm_y,herm_z,lin_xyz,wgt_cub_lin_xyz,lin_xy,
     & wgt_cub_lin_xyz_trc,
     & fhout_hf,fhmax_hf,
     & cont_eq_opt1,quamon,wgtm,time_extrap_etadot,iter_one_no_interp,
     & opt1_3d_qcubic,dfilevs,yhalo,
     & iau,iaufiles_fg,iaufiles_anl,iaufhrs,iau_delthrs,yhalo,iaulnp,
     & sppt,sppt_tau,sppt_lscale,sppt_logit,
     & vcamp,vc_lscale,vc_tau,vc_logit,
     & iseed_shum,iseed_sppt,shum,shum_tau,shum_lscale,iseed_vc,
     & fhstoch,stochini,vc,skeb_varspect_opt,sppt_sfclimit,
     & skeb,skeb_tau,skeb_lscale,iseed_skeb,skeb_vfilt,skeb_diss_smooth,
     & skeb_sigtop1,skeb_sigtop2,sppt_sigtop1,sppt_sigtop2,
     & shum_sigefold,vc_sigtop1,vc_sigtop2,skebint,
! WAM IPE coupling flags
     & wam_ipe_coupling, height_dependent_g,
! For WAM IPE coupling restart run.
!----------------------------------
     & wam_ipe_cpl_rst_output, wam_ipe_cpl_rst_input,
! For outputing the NetCDF diagnostic files.
!-------------------------------------------
     & NC_output, DELOUT_NC, nc_fields

!
      fhmax      = 0
      fhout      = 0
      fhres      = 0
      fhrot      = 0
      fhout_hf   = 1
      fhmax_hf   = 0
      fhdfi      = 0
      dfilevs    = levs
      igen       = 0
      tfiltc     = 0.85
      ngptc      = lonf

      cdamp(1)   = 0.0
      cdamp(2)   = 1.0
      k2o        = -1
      gg_tracers = .false.

      if (semilag) then
        xhalo      = 1
        yhalo      = 10
        redgg_a    = .true.
        sl_epsln   = 0.05
        phigs_d    = 60.0
        cdamp(1)   = 5.0e4
        cdamp(2)   = 1.5

        herm_x     = .true.
        herm_y     = .true.
        herm_z     = .true.
        lin_xyz    = .false.
        lin_xy     = .false.
        levwgt(1)  = 15
        levwgt(2)  = 35
        wgtm(1)    = 0.0
        wgtm(2)    = 1.0
        quamon     = .false.

        wgt_cub_lin_xyz     = .false.
        wgt_cub_lin_xyz_trc = .false.
        cont_eq_opt1        = .false.
        time_extrap_etadot  = .false.
        iter_one_no_interp  = .false.
        opt1_3d_qcubic      = .false.

        settls_dep3ds       = .false.
        settls_dep3dg       = .false.
      endif
!
      hdif_fac    = 1.0
      hdif_fac2   = 1.0
      slrd0       = 0.002 ! sigma level above which Raleigh damping applied
!
      shuff_lats_a     = .true.
      reshuff_lats_a   = .false.
!
      reduced_grid     = .true.
!
      explicit         = .false.
      hybrid           = .false.
      gen_coord_hybrid = .false.                                     !hmhj
      mass_dp          = .false.                                     !hmhj
      semi_implicit_temp_profile    = .false.                        !hmhj
      process_split    = .false.                                     !hmhj
      liope            = .true.
!
      thermodyn_id     = 1
      sfcpress_id      = 1
!
      zflxtvd          = .false.
      ndslfv           = .false. ! non_iteration semi_Lagrangian finite volume
      ldfi_spect       = .false. ! digital filter on spectral coefficients
!
      nemsio_in         = .true.
      nemsio_out        = .true.
      sigio_out         = .false.
!
!
      ref_temp          = 300.0
      ref_pres          = 101.325 ! in centibars - good value for 2 time level scheme
!
      nsout             = 0
      nsout_hf          = 0
      levr              = 0
      spectral_loop     = 2       ! 1 for one-loop or 2 for two-loop
! !  additions for stochastic perturbations...
!------------------------------------------
      ! can specify up to 5 values for the stochastic physics parameters 
      ! (each is an array of length 5)
      sppt             = -999.  ! stochastic physics tendency amplitude
      shum             = -999.  ! stochastic boundary layer spf hum amp   
      skeb             = -999.  ! stochastic KE backscatter amplitude
      vcamp            = -999.  ! stochastic vorticity confinment amp
! logicals
      do_sppt = .false.
      do_shum = .false.
      do_skeb = .false.
      do_vc   = .false.
! number of passes of 1-2-1 vertical filter
! for SKEB random patterns.
      skeb_vfilt       = 0
      skebint          = 0
      sppt_tau         = -999.  ! time scales
      shum_tau         = -999.
      skeb_tau         = -999.
      vc_tau           = -999.
      sppt_lscale      = -999.  ! length scales
      shum_lscale      = -999.
      skeb_lscale      = -999.
      vc_lscale        = -999.
      iseed_sppt       = 0      ! random seeds (if 0 use system clock)
      iseed_shum       = 0
      iseed_skeb       = 0
      iseed_vc         = 0
! parameters to control vertical tapering of stochastic physics with
! height
      sppt_sigtop1 = 0.1
      sppt_sigtop2 = 0.025
      skeb_sigtop1 = 0.1
      skeb_sigtop2 = 0.025
      vc_sigtop1 = 0.1
      vc_sigtop2 = 0.025
      shum_sigefold = 0.2
! reduce amplitude of sppt near surface (lowest 2 levels)
      sppt_sfclimit = .false.
! wavenumber defining smoothing of skeb dissipation estimate.
      skeb_diss_smooth = 12.
! gaussian or power law variance spectrum for skeb (0: gaussian, 1:
! power law). If power law, skeb_lscale interpreted as a power not a
! length scale.
      skeb_varspect_opt = 0
      vc                = 0.      ! deterministic vorticity confinement parameter.
      sppt_logit        = .false. ! logit transform for sppt to bounded interval [-1,+1]
      vc_logit          = .false. ! logit transform for vc to bounded interval [-1,+1]
      fhstoch           = -999.0  ! forecast hour to dump random patterns
      stochini          = .false. ! true= read in pattern, false=initialize from seed

! idea add
      lsidea              = .false. ! idea add (WAM) logical variablei 
      wam_ipe_coupling    = .false.
      height_dependent_g  = .false.

      wam_ipe_cpl_rst_input  = .false.
      wam_ipe_cpl_rst_output = .false.

      NC_output           = .false.
      nc_fields           = .true.
      DELOUT_NC           = 3600
!
!  iau parameters
       iau              = .false.
       iaufhrs          = -1
       iaufiles_fg      = ''
       iaufiles_anl     = ''
      if (me == 0) print *,' nlunit=',nlunit,' gfs_dyn_namelist=',
     &                      gfs_dyn_namelist
!$$$      read(5,nam_dyn)

      open(unit=nlunit,file=gfs_dyn_namelist)
      rewind (nlunit)
      read(nlunit,nam_dyn)

      if (me == 0) print *,' in compns_dynamics semilag=',semilag

      if (k2o < 0) k2o = levs/2

      if (semilag) then
        if (herm_x .or. herm_y .or. herm_z) then
          if (me == 0)
     &      print*,'hermite interpolation used - turning off lin_xyz,'
     &,            'wgt_cub_lin_xyz,quamon,lin_xy'
            lin_xyz             = .false.
            wgt_cub_lin_xyz     = .false.
            wgt_cub_lin_xyz_trc = .false.
            lin_xy              = .false.
            quamon              = .false.
            cont_eq_opt1        = .false.
        endif

        if (settls_dep3dg .neqv. settls_dep3ds) then
          if (me == 0)
     &      print*,'**error: settls_dep3dg should be same as'
     &,            ' settls_dep3ds  - setting both to true '
           settls_dep3dg  = .true.
           settls_dep3ds  = .true.
        endif
!
!  Allocate and set up weighting function for two time level semi-Lagangian
!  Here levels ae from top to bottom
!
        if (wgt_cub_lin_xyz .or. wgt_cub_lin_xyz_trc) then
          allocate (wgt_slg(levs))
          do k = 1,levwgt(1)
             wgt_slg(k) = wgtm(1)
          enddo
          tem = (wgtm(2) - wgtm(1)) / (levwgt(2) - levwgt(1) + 1)
          do k = levwgt(1)+1,levwgt(2)
             wgt_slg(k) = wgtm(1) + tem * (k-levwgt(1))
          enddo
          do k = levwgt(2)+1,levs
             wgt_slg(k) = wgtm(2)
          enddo
!         wgt_slg(:) = 1.0
          if (me == 0) print *,' wgt_slg=',wgt_slg
        endif
      endif               ! end of if(semilag)
!
      if( lsidea ) then       ! idea add (WAM) case needs larger ref_temp
        if (levs > 100 .and. ref_temp < 400.0) ref_temp = 2500.0
      endif

      if (me == 0) write(6,nam_dyn)
      filta = tfiltc
      if(semilag) then
        phigs = phigs_d * (atan(1.0)/45.0)
        if (me == 0) then
          write(0,*)' phigs_d=',phigs_d,' phigs=',phigs
        endif
      endif
!
      if (levr == 0) then
        levr = levs
      endif

!
!sela - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      tol = 0.01
!  Check rule 1.
      if (me == 0)
     & print *,'lver=',levr,'deltim=',deltim,'nsout=',nsout,'fhout=',
     & fhout,'fhres=',fhres,'gen_coord_hybrid=',gen_coord_hybrid, 
     & 'ntoz=',ntoz
     &,' fhout_hf=',fhout_hf,' fhmax_hf=',fhmax_hf,' in dynamics'

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Compute nsout and check rule 2.
      if(nsout > 0) fhout = nsout*deltim/3600.
      nsout = nint(fhout*3600./deltim)
      if(nsout <= 0 .or. abs(nsout-fhout*3600./deltim) > tol) then
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
!  Compute nsres and check rule 6.
      nsres=nint(fhres*3600./deltim)
      if(nsres <= 0 .or. abs(nsres-fhres*3600./deltim) > tol) then
        iret = 6
        return
      endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Compute ndfi and check rule 7.
      if(fhdfi == 0.) then
        ndfi = 0
        ldfi_spect = .false.
      else
        ndfi = nint(2*fhdfi*3600./deltim)
        if(ndfi <= 0 .or. abs(ndfi-2*fhdfi*3600./deltim) > tol .or.
     &     ndfi > nsres) then
           print *,'ndfi=',ndfi,'is not equal to',2*fhdfi*3600./deltim
          iret = 7
          return
        endif
      endif
! PJP stochastic physics additions
      IF (sppt(1) > 0 ) THEN
        do_sppt=.true.
      ENDIF
      IF (shum(1) > 0 ) THEN
        do_shum=.true.
!     shum parameter has units of 1/hour, to remove time step
!     dependence.
!     change shum parameter units from per hour to per timestep
         DO k=1,5
            IF (shum(k) .gt. 0.0) shum(k)=shum(k)*deltim/3600.0
         ENDDO
      ENDIF
      IF (skeb(1) > 0 ) THEN
         do_skeb=.true.
!      adjust skeb values for resolution.
!      scaling is such that a value of 1.0 at T574 with a 900 second
!      timestep produces well-calibrated values of forecast spread.
         DO k=1,5
            IF  (skeb(k) .gt. 0.0) THEN
               skeb(k)=skeb(k)*deltim/(jcap*(jcap+1))*365765.0  ! 365765 is a scale factor so the base SKEB value in the namelist is 1.0
            ENDIF
         ENDDO
      ENDIF
!    compute frequencty to estimate dissipation timescale
      IF (skebint == 0.) skebint=deltim
      nsskeb=nint(skebint/deltim)                              ! skebint in seconds
      IF(nsskeb<=0 .or. abs(nsskeb-skebint/deltim)>tol) THEN
         WRITE(0,*) "SKEB interval is invalid",skebint
        iret=9
        return
      ENDIF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF (vc > tiny(vc) .or. vcamp(1) > 0 ) THEN
        do_vc=.true.
      ENDIF
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
      if (ngptc > lonf) then
         ngptc = lonf
         WRITE(0,*) "NGPTC IS TOO BIG, RESET NGPTC TO lonf",NGPTC
      endif
      IF (ME == 0)   WRITE(0,*) "NGPTC IS SET TO NGPTC :",NGPTC

!  all checks are successful.
      iret=0
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  All checks are successful.
!
      iret = 0
!
      return
      end subroutine compns_dynamics
