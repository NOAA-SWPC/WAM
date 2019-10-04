!
! !module: gfs_physics_run_mod --- run module of the grided
!                              component of the gfs physics.
!
! !description: gfs run module.
!
! !revision history:
!
!  november 2004      weiyu yang initial code.
!  may      2005      weiyu yang, updated to the new version of gfs.
!  janusry  2007      hann-ming henry juang for gfs dynamics only
!  july     2007      shrinivas moorthi for gfs physics only
!  november 2007      hann-ming henry juang continue for gfs physics
!  october  2009      jun wang add nsout option
!  oct 11   2009      sarah lu, grid_gr replaced by grid_fld
!  oct 17   2009      sarah lu, q is replaced by tracers(1)
!  dec 08   2009      sarah lu, add g3d_fld to do_physics_one_step
!                     calling argument
!  July     2010      Shrinivas Moorthi - Updated for new physics and added nst
!                     eliminated calls to common_vars
!  jul 21  2010       sarah lu, add g2d_fld to do_physics_one_step
!                     calling argument
!  Aug 03  2010       jun wang, fix lsout for dfi
!  Aug 25  2010       Jun Wang, add zhour_dfi for filtered dfi fields output
!  Oct 18  2010       s. moorthi added fscav to do tstep
!  Dec 23  2010       Sarah Lu, setup fscav from gfs_phy_tracer 
!  Nov 27  2011       Sarah Lu, zerout fcld, dqdt, and wet1
!  Apr 06  2012       Henry Juang, add idea
!  Apr 09  2012       Jun Wang save phys state at 3hr and set back to 
!                     3hr phys state aft dfi
!  Mar 09  2013       Jun Wang add restart step for idea
!  Mar 25  2014       Xingren Wu add aoi_fld for A/O/I coupling
!  Mar 31  2014       S Moorthi  Add sstForGSM to do_physics_onestep argument
!  Jul --  2014       S Moorthi  update for new physics and semilag
!  Sep 18  2014       S Moorthi  simplified argments for dfi_fixwr
!  Sep 30  2014       Sarah Lu,  Remove fscav array
!  May --  2015       S Moorthi - added SHOC option
!  Jun 09  2015       G Theurich, Generalize importData handling
!  Jul 29  2015       S Moorthi - added high frequency output capability
!  Jan     2016       P Tripp - NUOPC/GSM merge - importData
!  feb     2016       S Moorthi - grid-point digital filter fix to filter at initial
!                                 time and during restart
!  mar 05  2016       S Moorthi - add logic to change phour back to reset time after
!                                 clock was reset after digital filter initialization
!
! !interface:
!
      module gfs_physics_run_mod
!
!!uses:
!
      use gfs_physics_internal_state_mod, ONLY: gfs_physics_internal_state
      USE date_def,                       ONLY: fhour, zhour
      USE namelist_physics_def,           ONLY: nsout,ldfi,ndfi,nsout_hf, fhmax_hf
      use gfs_phy_tracer_config,          ONLY: gfs_phy_tracer
      use layout1,                        only: me
      use resol_def,                      only: kdt_start

      implicit none

      contains

      subroutine gfs_physics_run(gis_phy, rc)

      type(gfs_physics_internal_state), pointer, intent(inout) :: gis_phy
      integer, optional,                         intent(out)   :: rc

      real , save      :: timestep=0.0
      integer             rc1, k , i1, i2, i3, kdt_dif
      logical             lsout1
      real                tem

!***********************************************************************
!
!     lsfwd      logical true during a first forward step
!     lssav      logical true during a step for which
!                diagnostics are accumulated
!     lscca      logical true during a step for which convective clouds
!                are calculated from convective precipitation rates
!     phour      real forecast hour at the end of the previous timestep
!
!     lsout controls freq. of output
!     nsout controls freq. of output in time steps
!     fhout controls freq. of output in hours
!     nszer time steps between zeroing fluxes
!
!***********************************************************************

!       print *,' enter gfs_physics_run '
!
       rc1 = 0
!
!      print *,' uug=',gis_phy%grid_gr(1,gis_phy%g_u:gis_phy%g_u+gis_phy%levs-1)
!      print *,' pg=',gis_phy%grid_gr(1,gis_phy%g_p:gis_phy%g_p+gis_phy%levs-1)
!      print *,' dpg=',gis_phy%grid_gr(1,gis_phy%g_dp:gis_phy%g_dp+gis_phy%levs-1)
!
! ---------------------------------------------------------------------
! ======================================================================
!                     do one physics time step
! ---------------------------------------------------------------------
!     write(0,*)' in gfs_physics_run kdt=',gis_phy%kdt,' nsout=',nsout  &
!              ,' kdt_start=',kdt_start

       if (fhmax_hf > 0 .and. nsout_hf > 0 .and. gis_phy%phour <= fhmax_hf    &
           .and. gis_phy%kdt*gis_phy%deltim <= fhmax_hf*3600.) then
         lsout1 = mod(gis_phy%kdt ,nsout_hf) == 0
       else
         lsout1 = mod(gis_phy%kdt ,nsout) == 0
       endif

       kdt_dif = gis_phy%kdt - kdt_start
       if (.not. ldfi) then
         gis_phy%lsout = lsout1 .or. gis_phy%kdt == 1
       else
         gis_phy%lsout = lsout1 .and.                                   &
           (kdt_dif <= ndfi/2 .or. kdt_dif > ndfi) .or. gis_phy%kdt == 1
       endif

       if( ndfi > 0 .and. kdt_dif == ndfi/2+1 .and. .not. ldfi )  then
         call dfi_fixwr(2, gis_phy%sfc_fld,  gis_phy%nst_fld)
       endif
        

        if ( gis_phy%lgocart ) then
!          i1 = gis_phy%gfs_phy_tracer%ntrac_met+1 ! 1st chemical tracer (excluding o3)
!          i2 = gis_phy%gfs_phy_tracer%ntrac       ! last chemical tracer (excluding o3)
           i1 = gfs_phy_tracer%ntrac_met+1         ! 1st chemical tracer (excluding o3)
           i2 = gfs_phy_tracer%ntrac               ! last chemical tracer (excluding o3)
!
           i1 = gis_phy%lonr
           i2 = gis_phy%lats_node_r_max
           i3 = gis_phy%levs
           gis_phy%flx_fld%wet1(1:i1,1:i2)         = 0.
           gis_phy%g3d_fld%fcld(1:i1,1:i2,1:i3)    = 0.
           gis_phy%g3d_fld%dqdt(1:i1,1:i2,1:i3)    = 0.
           gis_phy%g3d_fld%cnv_mfc(1:i1,1:i2,1:i3) = 0.
           gis_phy%g3d_fld%cnv_mfd(1:i1,1:i2,1:i3) = 0.
           gis_phy%g3d_fld%cnv_qc(1:i1,1:i2,1:i3)  = 0.
        endif

        tem = (gis_phy%kdt-1) * gis_phy%deltim / 3600.0
        if (gis_phy%phour > tem) gis_phy%phour = tem 

!       gis_phy%phour = gis_phy%kdt * gis_phy%deltim / 3600.
!     if (me ==0) write(0,*)'gis_phy%phour=',gis_phy%phour,'gis_phy%kdt=',&
!         gis_phy%kdt,' gis_phy%deltim',gis_phy%deltim
!
! ======================================================================
        call do_physics_one_step(                                         &
                 gis_phy%deltim,   gis_phy%kdt,     gis_phy%phour,        &
                 gis_phy%grid_fld, gis_phy%sfc_fld, gis_phy%flx_fld,      &
                 gis_phy%nst_fld,  gis_phy%g3d_fld, gis_phy%g2d_fld,      &
                 gis_phy%aoi_fld,  gis_phy%importData,                    &
                 gis_phy%lats_nodes_r,   gis_phy%global_lats_r,           &
                 gis_phy%lonsperlar,                                      &
                 gis_phy%XLON,    gis_phy%XLAT,    gis_phy%COSZDG,        &
                 gis_phy%HPRIME,  gis_phy%SWH,     gis_phy%swhc,          &
                 gis_phy%HLW,     gis_phy%hlwc,                           &
! idea add by hmhj - commented by moorthi since unused
!                gis_phy%HTRSWB,  gis_phy%HTRLWB,                         &
                 gis_phy%FLUXR,   gis_phy%SFALB,                          &
                 gis_phy%SLAG,    gis_phy%SDEC,    gis_phy%CDEC,          &
                 gis_phy%OZPLIN,  gis_phy%JINDX1,  gis_phy%JINDX2,        &
                 gis_phy%DDY,                                             &
                 gis_phy%phy_f3d, gis_phy%phy_f2d, gis_phy%phy_fctd,      &
                 gis_phy%num_ctp,                                         &
                 gis_phy%NBLCK,                    gis_phy%ZHOUR_DFI,     &
!                gis_phy%NBLCK,   gis_phy%ZHOUR,   gis_phy%ZHOUR_DFI,     &
                 gis_phy%N3,      gis_phy%N4,                             &
                 gis_phy%LSOUT,   gis_phy%COLAT1,  gis_phy%CFHOUR1,       &
                 gis_phy%restart_step, gis_phy%mdl_parm)

!                        
! =======================================================================
!
! save phys fields for digital filter
!
       if( ldfi .and. kdt_dif == ndfi/2 )  then
!         write(0,*)'save phys state, at gis_phy%kdt=',gis_phy%kdt,'ldfi=',ldfi
         call dfi_fixwr(1, gis_phy%sfc_fld,  gis_phy%nst_fld)
       endif
!
      gis_phy%phour = fhour                        ! update hour
      gis_phy%zhour = zhour                        ! update hour
!
      if(present(rc)) then
          rc = rc1
      end if

      end subroutine gfs_physics_run

      end module gfs_physics_run_mod
