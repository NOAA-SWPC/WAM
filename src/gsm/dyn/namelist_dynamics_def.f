      module namelist_dynamics_def

!
! program lot
! 06 Apr 2012:    Henry Juang add some options for NDSL
! 05 Oct 2012:    Jun Wang    add sigio_out
! 02 Apr 2014:    Jun Wang    add dfilevs
! Mar 07 2014     Weiyu Yang - add wam_ipe_coupling, height_dependent_g
! 02 May 2014:    Philip Pegion add stochastic variables
! 06 NEV 2017:    Weiyu Yang -
!                 i)   add wam_ipe_cpl_rst_input for WAM-IPE coupling
!                 restart run,
!                 ii)  add wam_ipe_cpl_rst_output for WAM-IPE coupling
!                 restart run,
!                 iv)  add NC_output and FHOUT_NC for outputing the
!                      NetCDF diagnostic files.
!
      use gfs_dyn_machine
      implicit none
      
      integer nsres,nsout,igen,ngptc,levwgt(2),k2o,nsout_hf
      integer dfilevs,nsskeb
      integer FHOUT_NC
      real(kind=kind_evod) fhrot,fhmax,fhout,fhres,fhini,fhdfi
      real(kind=kind_evod) filta,ref_temp,sl_epsln,cdamp(2)
     &,                    ref_pres,fhout_hf,fhmax_hf
      real(kind=kind_evod) hdif_fac,hdif_fac2,slrd0,wgtm(2)
      REAL(KIND = kind_evod) :: phigs,phigs_d
      logical lsfwd,ldfi_spect, shuff_lats_a,reshuff_lats_a
      logical,target :: hybrid,gen_coord_hybrid
      logical zflxtvd,explicit

      logical nemsio_in, nemsio_out, sigio_out
      logical reduced_grid, semi_implicit_temp_profile
      logical mass_dp, process_split
!
      logical herm_x,  herm_y,  herm_z,  lin_xyz, wgt_cub_lin_xyz,lin_xy
     &,       semilag, redgg_a, gg_tracers
     &,       wgt_cub_lin_xyz_trc
     &,       time_extrap_etadot,settls_dep3ds,settls_dep3dg
     &,       iter_one_no_interp,cont_eq_opt1,opt1_3d_qcubic

      logical ndslfv
! hmhj idea add
      logical lsidea,iau
      integer    ::iaulnp=2 ! 0: no pressure update, 1: log pressure, 2: pressure
      integer    :: iau_delthrs = 6 ! iau time interval (to scale increments)
      character(len=120), dimension(7) ::  iaufiles_fg,iaufiles_anl
      real(kind=kind_evod), dimension(7) :: iaufhrs
! WAM IPE coupling flags.
!------------------------
      logical :: wam_ipe_coupling, height_dependent_g, NC_output
     &           wam_ipe_cpl_rst_input, wam_ipe_cpl_rst_output
! pjp stochastic phyics
      integer skeb_varspect_opt
      logical sppt_sfclimit

      real(kind=kind_evod) :: skeb_sigtop1,skeb_sigtop2,
     &                   sppt_sigtop1,sppt_sigtop2,shum_sigefold,
     &                   vc_sigtop1,vc_sigtop2
      real(kind=kind_evod) fhstoch,vc,skeb_diss_smooth,skebint
      real(kind=kind_evod), dimension(5) ::
     &   skeb,skeb_lscale,skeb_tau
      real(kind=kind_evod), dimension(5) ::
     &   sppt,sppt_lscale,sppt_tau
      real(kind=kind_evod), dimension(5) ::
     &   vcamp,vc_lscale,vc_tau
      real(kind=kind_evod), dimension(5) ::
     &   shum,shum_lscale,shum_tau
      integer(8),dimension(5) ::
     &   skeb_vfilt,iseed_sppt,iseed_vc,iseed_shum,iseed_skeb
      logical stochini,vc_logit,sppt_logit
      logical do_shum,do_sppt,do_skeb,do_vc

      character*20 ens_nam
!
      end module namelist_dynamics_def
