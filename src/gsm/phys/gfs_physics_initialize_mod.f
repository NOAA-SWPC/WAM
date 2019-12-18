!#include "../../../ESMFVersionDefine.h"

! !module: gfs_physics_initialize_mod 
!          --- initialization module of the gridded component of gfs physics.
!
! !description: gfs physics gridded component initialize module.
!
! !revision history:
!
!  november 2004  weiyu yang    initial code.
!  january  2006  s. moorthi    update to the new gfs version
!  august   2006  h. juang      add option to run generalized coordinates
!  january  2007  h. juang      change for dynamics only
!  july     2007  s. moorthi    change for physics only
!  november 2007  h. juang      continue for physics
!  may      2009  j. wang       change for quilt
!  oct 09   2009  Sarah Lu      coord def initialized (lats_nodes_r_fix,
!                               lats_node_r, ipt_lats_node_r)
!  oct 11   2009  Sarah Lu      grid_gr is replaced by grid_fld
!  oct 12   2009  Sarah Lu      initialize start_step
!  oct 16   2009  Sarah Lu      initialize gfs_phy_tracer
!  nov 14   2009  Sarah Lu      flx_fld and sfc_fld allocation modified
!  dec 10   2009  Sarah Lu      initialize lgocart and g3d_fld
!  jan 22   2010  Sarah Lu      increase ngrids_flx and nfxr to include aod 
!  feb 09   2010  Sarah Lu      set tracer_const (ri,cpi) from import state
!  feb 05   2010  J. Wang       put phy_f3d and phy_f2d into restart file
!  apr 09   2010  Sarah Lu      initialize global_lats_r, lonsperlar
!  July     2010  S. Moorthi    Upgrade to new physics + nst model
!  jul 14   2010  Sarah Lu      initialize g2d_fld
!  jul 23   2010  Sarah Lu      initialize ngrids_aer and buff_mult_pieceg
!  July 30  2010  S. Moorthi    Removed allocation of cldcov
!  Aug 19   2010  S. Moorthi    Updated for T574 + added num_reduce to namelist
!  Oct 18   2010  S. Moorthi    Added fscav initialization
!  Dec 23   2010  Sarah Lu      initialize scatter_lats, scatter_lons, g2d_fld%met
!  Mar 27   2010  J. Wang       add zsoil to sfc file
!  Oct 03   2011  W. Yang       Modified for using the ESMF 5.2.0r library.
!  Apr 06   2012  H. Juang      add idea
!  Aug 20   2013  S. Moorthi    Updating for YuTai's new radiation package in gfs
!  Nov 23   2013  Sarah Lu      specify climate in internal state from namelist
!  Mar 18   2014  Sarah Lu      Remove iaer_mdl in rad_init
!  Mar 25   2014  Xingren Wu    initialize aoi_fld
!  Mar 31   2014  S. Moorthi    Allocate and Initialize the array sstForGSM
!  Jun 26   2014  S. Moorthi    Modified to read lonsperlar from a file
!  Jul 11   2014  S. Moorthi    Add npdf3d for pdf clouds
!  Jul 18   2014  S. Moorthi    removed num_reduce
!  Sep 30   2014  Sarah Lu      Remove fscav initialization
!  Nov 17   2014  S. Moorthi    Add missing dxmin, dxmax and dxinv definitions (my error)
!  Feb 05   2015  s. moorthi    several changes including stochphys switch to reduce memory
!  Apr 09   2015  s. moorthi    cs convection related change
!  Jun 09   2015  G Theurich    Generalize importData handling
!  aug      2015  s. moorthi    implement SHOC related additions
!  Aug 10   2015  s. moorthi    add call to sfc_init
!  Aug 21   2015  Xu Li         change nst_fcst to be nstf_name
!  Jan      2016  P. Tripp      NUOPC/GSM merge - importData
!  March    2016  J. Han        Add ncnvcld3d for enhancing conv clouds
!  March    2016  Hang Lei      Initialize the physics variable for mdl_parm
!  Sept     2017  W. Yang       Add the IPE back coupling to WAM code,
! !interface:
!
      module gfs_physics_initialize_mod

!
!!uses:
!
      USE esmf
      USE gfs_physics_internal_state_mod, ONLY: gfs_physics_internal_state
      use nuopc_physics,  &
        only: model_parameters, use_nuopc, nuopc_phys_init,  &
              phys_init_savein, phys_init_readin
!     USE mpi_def,        ONLY: liope
      USE mpi_def,        ONLY: mc_comp, mc_io, mpi_comm_all_dup, mpi_comm_all
      USE resol_def,      ONLY: g_dpdt, lotgr, g_p, g_dp, nrcm, g_q,            &
                                g_gz, g_u, g_v, g_ps, g_t, ntoz,                &
                                ntcw, lonr, latr, ncld, num_p3d, num_p2d,npdf3d,&
                                ncnvcld3d,lsoil,nmtvr,levr,nlunit, ntrac, nxpt, &
                                jcap, levs, nypt, jintmx, ngrids_sfcc,          &
!jw
                                ngrids_sfcc2d,ngrids_sfcc3d,                    & !jwang
                                ngrids_aer,                                     & !sarah lu
                                ngrids_aoi,                                     & !xwu
                                ngrids_flx, levp1, lonrx, nfxr, ngrids_gg,      &
                                levm1, ivssfc, thermodyn_id, sfcpress_id,       &
                                ivssfc_restart, latrd, latr2, ivsupa, levh,     &
                                lgocart, scatter_lats, scatter_lons,            &
                                nr_nst, nf_nst, ngrids_nst, ivsnst,             &
                                nshoc_3d, nshoc_2d, ntke, ntot3d, ntot2d
!jw
      use mod_state,      ONLY: buff_mult_piecea2d,buff_mult_piecea3d,          &  !jwang
                                buff_mult_piecef,buff_mult_pieceg,              &
                                buff_mult_piecenst
      use coordinate_def, ONLY: ak5,bk5,ck5                                        !jwang
      use vert_def,       only: sl,si                                              !  Moorthi

      USE ozne_def,       ONLY: levozc, latsozp, blatc, timeozc, timeoz,        &
                                kozpl, levozp, pl_time, pl_lat, pl_pres,        &
                                kozc, dphiozc, latsozc, pl_coeff
      USE namelist_physics_def,                                                 &
                          ONLY: ras, cscnv, jo3, ldiag3d, ngptc, ens_nam,       &
                                reduced_grid, grid_aldata, nstf_name,           &
                                isol, ico2, ialb, iems, iaer, iovr_sw,          &
                                iovr_lw,ictm, isubc_sw,isubc_lw,                &
                                crick_proof, ccnorm, norad_precip,              &
                                hybrid, gen_coord_hybrid, climate,              &
                                a2oi_out, cplflx, stochphys, do_shoc,           &
! Add some things needed by physics wrapper
                                flipv, pre_rad, lsm, imfshalcnv,imfdeepcnv,     &
                                ncw, crtrh, cdmbgwd, ccwf, dlqf, ctei_rm, cgwf, &
                                prslrd0, ral_ts,old_monin, cnvgwd, shal_cnv,    &
                                cal_pre, mom4ice, mstrat, trans_trac, moist_adj,&
                                lsidea, pdfcld, shcnvcw, redrag, hybedmf,       &
                                dspheat, shoc_cld, shocaftcnv,nemsio_in

      USE module_ras,     ONLY: nrcmax, fix_ncld_hr
      use cs_conv       , only: nctp
      use physcons,       only: max_lon, max_lat, min_lon, min_lat, dxmin, dxmax, dxinv

      USE gg_def,         ONLY: sinlat_r, coslat_r, wgtcs_r, rcs2_r, wgt_r,   &
                                colrad_r
      USE layout1,        ONLY: nodes_comp, lats_node_r_max, lats_node_r,     &
                                nodes, lon_dims_ext, lon_dims_r,idrt,         &
                                ipt_lats_node_r, me, comp_task
      USE date_def,       ONLY: idate, fhour
!*    USE tracer_const,   ONLY: set_tracer_const

      USE gfs_physics_sfc_flx_set_mod,                                  &
                          ONLY: sfcvar_aldata, flxvar_aldata, flx_init, sfc_init
      USE d3d_def,        ONLY: d3d_init, d3d_zero
      use wam_jh_integral,only: jh_integral_init, jh_integral_zero
      use machine,        ONLY : kind_io4
      USE sfcio_module,   ONLY: sfcio_axdbta

      USE module_IPE_to_WAM, ONLY:  lowst_ipe_level,                   &
                                    ZMT, MMT, JHR, SHR, O2DR,          &
                                    ipe_to_wam_coupling

      USE gfs_physics_gridgr_mod,  ONLY: gridvar_aldata  
      USE gfs_physics_g3d_mod,     ONLY: g3d_aldata
      USE gfs_physics_g2d_mod,     ONLY: g2d_aldata
      use gfs_phy_tracer_config,   ONLY: gfs_phy_tracer, tracer_config_init
      use gfs_physics_nst_var_mod
      use module_radsw_parameters, only: nbdsw
      use module_radlw_parameters, only: nbdlw
      use gfs_physics_aoi_var_mod
      use module_CPLFIELDS, only: NImportFields
!#ifndef IBM
!     USE omp_lib
!#endif
      implicit none

      include 'mpif.h'

      contains

      subroutine gfs_physics_initialize(gis_phy, rc)

! this subroutine sets up the internal state variables,
! allocate internal state arrays for initializing the gfs system.
!----------------------------------------------------------------

      integer, parameter :: iunit=101
      type(gfs_physics_internal_state), pointer, intent(inout) :: gis_phy
      integer,                                   intent(out)   :: rc

      integer                            :: ierr

      integer                            :: i, j, k, l, n, ilat, locl, ikey,&
                                            nrank_all,nblck,latrhf,iret
      integer                            :: num_parthds
      logical                            :: file_exists=.false.
!!
      real (kind=kind_io4) blatc4
      real (kind=kind_io4), allocatable  :: pl_lat4(:)
      real (kind=kind_io4), allocatable  :: pl_pres4(:)
      real (kind=kind_io4), allocatable  :: pl_time4(:)
!
      real (kind=kind_phys), allocatable :: si_loc(:)
      real (kind=kind_phys), parameter   :: typical_pgr=95.0
!     integer, parameter :: iflip = 0      ! virtical profile index from toa
                                           ! to sfc
      integer, parameter :: iflip = 1      ! virtical profile index from sfc
                                           ! to toa
      real (kind=8) :: phys_ini_time=0, btime, timef
!!
!     include 'function2'
!!
      btime = timef()

! set up gfs internal state (gis_phy) dimension and values for physics etc
!-------------------------------------------------------------------
      me     = gis_phy%me
      nodes  = gis_phy%nodes
!      CALL ESMF_VMGetCurrent(vm, rc = ierr)

      call compns_physics(gis_phy%nam_gfs_phy%deltim, gis_phy%iret,     &
                          gis_phy%ntrac,                                &
                          gis_phy%nxpt,   gis_phy%nypt,  gis_phy%jintmx,&
                          gis_phy%jcap,   gis_phy%levs,  gis_phy%levr,  &
                          gis_phy%lonr,   gis_phy%latr,                 &
                          gis_phy%ntoz,   gis_phy%ntcw,  gis_phy%ncld,  &
                          gis_phy%ntke,   gis_phy%lsoil,  gis_phy%nmtvr,&
                          gis_phy%num_p3d, gis_phy%num_p2d,             &
                          gis_phy%npdf3d, gis_phy%ncnvcld3d,            &
                          gis_phy%nshoc_3d, gis_phy%nshoc_2d,           &
                          gis_phy%ntot3d, gis_phy%ntot2d,               &
                          gis_phy%thermodyn_id, gis_phy%sfcpress_id,    &
                          gis_phy%nam_gfs_phy%nlunit, me,               &
                          gis_phy%nam_gfs_phy%gfs_phy_namelist)

      if (me == 0) write(0,*)' after compns_physics ntke=', gis_phy%ntke
!
!     This is also called in nuopc_phys_init, don't run twice if using the wrapper
!       if ( .not. use_nuopc ) then
      CALL set_soilveg(me,gis_phy%nam_gfs_phy%nlunit)
!       end if
!* ri/cpi is filled from dyn export state attributes (Sarah Lu)
!*      call set_tracer_const(gis_phy%ntrac,me,gis_phy%nam_gfs_phy%nlunit)
!
! met+chem tracer specification (Sarah Lu)
! NOTE: This config_init call repeats the init routine in dyc gc.  
!       The redundant calls will be removed in a later revision.  
!       The tracer specification will then be passed in from dyn dc
!
!     call tracer_config_init( gis_phy%gfs_phy_tracer, gis_phy%ntrac,   &
      call tracer_config_init( gis_phy%ntrac, gis_phy%ntoz,             &
                               gis_phy%ntcw, gis_phy%ncld,              &
                               gis_phy%ntke,  me )
!      gfs_phy_tracer = gis_phy%gfs_phy_tracer
      gis_phy%lgocart = gfs_phy_tracer%doing_GOCART     ! for internal state
      lgocart = gis_phy%lgocart                         ! for resol_def module
      if( me == 0) then
       write(0,*)'LU_TRC, ntrac     =',gfs_phy_tracer%ntrac, gis_phy%ntrac
       write(0,*)'LU_TRC, ntrac_met =',gfs_phy_tracer%ntrac_met
       write(0,*)'LU_TRC, ntrac_chem=',gfs_phy_tracer%ntrac_chem
       write(0,*)'LU_TRC, lgocart   =',gis_phy%lgocart,lgocart
       do n = 1, gfs_phy_tracer%ntrac
         write(0,*)'LU_TRC, tracer_vname=',gfs_phy_tracer%vname(n)
       enddo
      endif

!
      nlunit   = gis_phy%nam_gfs_phy%nlunit
      ntrac    = gis_phy%ntrac
      nxpt     = gis_phy%nxpt
      nypt     = gis_phy%nypt
      jintmx   = gis_phy%jintmx
      jcap     = gis_phy%jcap
      levs     = gis_phy%levs
      levr     = gis_phy%levr
      lonr     = gis_phy%lonr
      latr     = gis_phy%latr
      ntoz     = gis_phy%ntoz
      ntcw     = gis_phy%ntcw
      ncld     = gis_phy%ncld
      ntke     = gis_phy%ntke
      lsoil    = gis_phy%lsoil
      nmtvr    = gis_phy%nmtvr
      num_p3d  = gis_phy%num_p3d
      num_p2d  = gis_phy%num_p2d
      npdf3d   = gis_phy%npdf3d
      ncnvcld3d= gis_phy%ncnvcld3d
      nshoc_3d = gis_phy%nshoc_3d
      nshoc_2d = gis_phy%nshoc_2d
      ntot3d   = gis_phy%ntot3d
      ntot2d   = gis_phy%ntot2d
      thermodyn_id = gis_phy%thermodyn_id
      sfcpress_id  = gis_phy%sfcpress_id
      if (gis_phy%nam_gfs_phy%Total_Member <= 1) then
        ens_nam=' '
      else
        write(ens_nam,'("_",I2.2)') gis_phy%nam_gfs_phy%Member_Id
      endif
!
!     ivssfc  = 200501
      ivssfc  = 200509
      ivssfc_restart  = 200509
      if (ivssfc > ivssfc_restart) ivssfc_restart = ivssfc
      ivsnst  = 200907
!
      allocate(gis_phy%zsoil(lsoil))
      if(lsoil == 2) then
        gis_phy%zsoil = (/-0.1,-2.0/)   
      elseif(lsoil == 4) then
        gis_phy%zsoil = (/-0.1,-0.4,-1.0,-2.0/)
      endif

      ivsupa  = 0
      if (levs > 99) ivsupa  = 200509
!
      levh   = ntrac*levs
      latrd  = latr + 2*jintmx
      latr2  = latr/2
      levm1  = levs-1 
      levp1  = levs+1 
      lonrx  = lonr + 1 + 2*nxpt + 1
!
      ngrids_sfcc = 32+LSOIL*3   ! No CV, CVB, CVT! includes T2M, Q2M, TISFC
!jw
      ngrids_sfcc2d = 32        ! No CV, CVB, CVT! includes T2M, Q2M, TISFC
      ngrids_sfcc3d = LSOIL*3   ! for smc,stc,slc

     gis_phy%climate = climate
     if (climate) then
!      ngrids_flx  = 66+36+8  ! additional 8 gocart avg output fields
!      ngrids_flx  = 66+36+8+4! additional 8 gocart, 4 sw fluxes
       ngrids_flx  = 66+36+5  ! additional 4 sw fluxes + frozen precip fraction
     else
!      ngrids_flx  = 66+43+8  ! additional 8 gocart avg output fields
!      ngrids_flx  = 66+43+8+4! additional 8 gocart, 4 sw fluxes
       ngrids_flx  = 66+43+5  ! additional 4 sw fluxes + frozen precip fraction
     endif
 
      nfxr        = 39               ! Add AOD
      ngrids_gg   = 2+LEVS*(4+ntrac)
!
      ngrids_aoi  = 40
!
      if (nstf_name(1) > 0) then         ! For NST model
!       ngrids_nst = 19              ! oceanic fields (for diurnal warming and sub-layer)
        nr_nst = 10                  ! oceanic fields: for diurnal warming modelrun
        nf_nst = 9                   ! oceanic fields: for GSI analysis
        ngrids_nst = nr_nst + nf_nst ! oceanic fields (for diurnal warming and sub-layer)
      else
        ngrids_nst = 0
      endif
!
      allocate(ak5(levp1))
      allocate(bk5(levp1))
      allocate(ck5(levp1))
      allocate(si(levp1))
      allocate(sl(levs))

      allocate ( lon_dims_r(latr),    stat = ierr )
      allocate ( lon_dims_ext(latrd), stat = ierr )
!
      allocate(colrad_r(latr), stat = ierr)
      allocate(wgt_r(latr2),   stat = ierr)
      allocate(wgtcs_r(latr2), stat = ierr)
      allocate(rcs2_r(latr2),  stat = ierr)
      allocate(sinlat_r(latr), stat = ierr)
      allocate(coslat_r(latr), stat = ierr)
!
      allocate ( gis_phy%lats_nodes_r(nodes),     stat = ierr )
      allocate ( gis_phy%lats_nodes_ext(nodes),   stat = ierr )
      allocate ( gis_phy%global_lats_r(latr),     stat = ierr )
      allocate ( gis_phy%global_lats_ext(latr),   stat = ierr )
      allocate ( gis_phy%lonsperlar(latr),        stat = ierr)
      allocate ( gis_phy%lats_nodes_r_fix(nodes), stat = ierr )    !added for mGrid
!     allocate ( gis_phy%global_lats_ext(latr+2*jintmx+2*nypt*(nodes-1)), stat = ierr )

      if (me == 0) write(0,*)' in gfs_physics_initialize reduced_grid=',reduced_grid

      if( reduced_grid ) then
        if (me == 0) print *,' run with reduced gaussian grid '
        inquire (file="lonsperlar.dat", exist=file_exists)
        if ( .not. file_exists ) then
          if ( me == 0 ) then
            print *,'   Requested lonsperlar.dat  data file does not exist'
            print *,'   *** Stopped in subroutine GFS_Init !!'
          endif
          call mpi_quit(1111)
        else
          open (iunit,file='lonsperlar.dat',status='old',form='formatted',     &
                                            action='read',iostat=iret)
          if (iret /= 0) then
            write(0,*)' iret while reading lonsperlar.dat ',iret
            call mpi_quit(1112)
          endif
          rewind iunit
          read (iunit,*,iostat=iret) latrhf,(gis_phy%lonsperlar(i),i=1,latrhf)
          if (latrhf+latrhf /= latr) then
             write(0,*)' latrhf=',latrhf,' not equal to latr/2=',latr/2
             call mpi_quit(1113)
          endif
          do i=1,latrhf
            gis_phy%lonsperlar(latr-i+1) = gis_phy%lonsperlar(i)
          enddo
          close(iunit)
        endif
!         write(0,*)' gis_phy%lonsperlar=',gis_phy%lonsperlar
      else
        if (me == 0) print *,' run with full gaussian grid '
        do j=1,latr
          gis_phy%lonsperlar(j) = lonr
        enddo
      endif
!
      g_gz   = 1
      g_ps   = g_gz  + 1
      g_t    = g_ps  + 1     
      g_u    = g_t   + levs
      g_v    = g_u   + levs
      g_q    = g_v   + levs
      g_p    = g_q   + levh
      g_dp   = g_p   + levs
      g_dpdt = g_dp  + levs
       
      lotgr  = g_dpdt+ levs - 1

      gis_phy%g_gz    = g_gz  
      gis_phy%g_ps    = g_ps  
      gis_phy%g_t     = g_t  
      gis_phy%g_u     = g_u  
      gis_phy%g_v     = g_v  
      gis_phy%g_q     = g_q  
      gis_phy%g_p     = g_p  
      gis_phy%g_dp    = g_dp 
      gis_phy%g_dpdt  = g_dpdt 
!
      gis_phy%lotgr = lotgr


      if (ras) then
        if (fix_ncld_hr) then
          nrcm = min(nrcmax, levs-1) * (gis_phy%deltim/1200) + 0.10001
        else
          nrcm = min(nrcmax, levs-1)
        endif
      else
        nrcm = 2
      endif
!
!     write(0,*)' before ozone read'
      if (ntoz < 0) then      ! Diagnostic ozone
        rewind (kozc)
        read (kozc,end=101) latsozc, levozc, timeozc, blatc4
  101   if (levozc  < 10 .or. levozc > 100) then
          rewind (kozc)
          levozc  = 17
          latsozc = 18
          blatc   = -85.0
        else
          blatc   = blatc4
        endif
        latsozp   = 2
        levozp    = 1
        timeoz    = 1
        pl_coeff  = 0
      else                       ! Prognostic Ozone
        rewind (kozpl)
        read (kozpl) pl_coeff, latsozp, levozp, timeoz
        allocate (pl_lat(latsozp), pl_pres(levozp),pl_time(timeoz+1), stat = ierr)
        allocate (pl_lat4(latsozp), pl_pres4(levozp),pl_time4(timeoz+1), stat = ierr)
        rewind (kozpl)
        read (kozpl) pl_coeff, latsozp, levozp, timeoz, pl_lat4, pl_pres4,  &
                     pl_time4
        pl_pres(:) = pl_pres4(:)
        pl_lat(:)  = pl_lat4(:)
        pl_time(:) = pl_time4(:)
        latsozc = 2
        blatc   = 0.0
        if (allocated(pl_lat4))  deallocate(pl_lat4)
        if (allocated(pl_pres4)) deallocate(pl_pres4)
        if (allocated(pl_time4)) deallocate(pl_time4)
      endif
      dphiozc = -(blatc+blatc)/(latsozc-1)
!
      if (me  == 0) then

!       print *,' g_gz ',g_gz
!       print *,' g_ps ',g_ps
!       print *,' g_t  ',g_t 
!       print *,' g_u  ',g_u 
!       print *,' g_v  ',g_v 
!       print *,' g_q  ',g_q 
!       print *,' g_p  ',g_p 
!       print *,' g_dp ',g_dp
!       print *,' g_dpdt ',g_dpdt
!       print *,' lotgr ',lotgr
        print *,' latsozp=',latsozp,' levozp=',levozp,' timeoz=',timeoz
        print *,' latsozc=',latsozc,' levozc=',levozc,' timeozc=',        &
                  timeozc, 'dphiozc=',dphiozc
!       print *,' pl_lat=',pl_lat
!       print *,' pl_pres=',pl_pres
!       print *,' pl_time=',pl_time

      endif

!     pl_pres(:) = log(0.1*pl_pres(:))       ! Natural log of pres in cbars
      pl_pres(:) = log(100.0*pl_pres(:))     ! Natural log of pres in Pa
!
      allocate(gis_phy%OZPLIN(LATSOZP,LEVOZP,pl_coeff,timeoz), stat = ierr) !OZONE P-L coeffcients

!
!   Compute dxmax, dxmin, and dxinv used in gbphys and grrad
!
      dxmax = log(1.0/(max_lon*max_lat))
      dxmin = log(1.0/(min_lon*min_lat))
      dxinv = 1.0 / (dxmax-dxmin)
      if (me  == 0) write(0,*)' dxmax=',dxmax,' dxmin=',dxmin,' dxinv=',dxinv

!      write(0,*)' finished array allocation in gfs_physics_initialize '

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!!
!!      create io communicator and comp communicator
!!
!      if (me == 0) write(*,*) 'io option ,liope :',liope
!
!jw      call mpi_comm_dup(mpi_comm_all, mpi_comm_all_dup, ierr)
!     call mpi_barrier (mpi_comm_all_dup,               ierr)
!
!      if (nodes == 1) liope=.false.
!jw      if (liope) then
!jw        call mpi_comm_rank(mpi_comm_all_dup,nrank_all,ierr)
!jw        icolor=1
!jw        ikey=1
!jw        nodes_comp=nodes-1
!jw        if (nrank_all.eq.nodes-1) then
!!  io server
!jw          write(*,*) 'io server task'
!jw          icolor=2
!jw          gis_phy%kcolor=mpi_undefined
!jw          call mpi_comm_split(mpi_comm_all_dup,icolor,ikey,mc_io,ierr)
!jw          call mpi_comm_split(mpi_comm_all_dup,gis_phy%kcolor,ikey,mc_comp,ierr)
!jw        else
!sela     write(*,*) 'compute server task '
!jw          icolor=mpi_undefined
!jw          gis_phy%kcolor=1
!jw          call mpi_comm_split(mpi_comm_all_dup,gis_phy%kcolor,ikey,mc_comp,ierr)
!jw          call mpi_comm_split(mpi_comm_all_dup,icolor,ikey,mc_io,ierr)
!jw          call mpi_comm_size(mc_comp,nodes,ierr)
!jw        endif
!jw      else
!jw        icolor=2
!jw        mc_comp=mpi_comm_all_dup
!jw      endif
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!!
        nodes_comp = nodes
        comp_task  = me < nodes_comp
!c
!$$$      time0=timer()
!
      if (me == 0) then
!       print 100, jcap,levs
!100     format (' smf ',i3,i3,' in gfs physics initialize ')
!#ifdef IBM
        print*,'number of threads is ',num_parthds()
!#else
!       print*,'number of threads is ',omp_get_num_threads()
!#endif
!jw        if (liope) then
!jw          print*,'number of mpi procs is ',nodes
!jw          print*,'number of mpi io procs is 1 (nodes)'
!jw        else
          print*,'number of mpi procs is ',nodes
!jw        endif
      endif
!c
      gis_phy%cons0    =    0.0d0
      gis_phy%cons0p5  =    0.5d0
      gis_phy%cons1200 = 1200.d0
      gis_phy%cons3600 = 3600.d0
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
      if(gis_phy%iret.ne.0) then
        if(me == 0) print *,' incompatible physics namelist -',         &
                            ' aborted in gfs_phy_initilize ',gis_phy%iret 
!*                          ' aborted in gfs_phy_initilize'
        call mpi_quit(13)
      endif
!!
!     if predicted ozon is desired set jo3=2
      jo3 = 2          !using predicted ozone in radiation.
! 
      gis_phy%lats_nodes_ext = 0

!       write(0,*)' gis_phy%lonsperlar2d=',gis_phy%lonsperlar
      call getcon_physics(gis_phy%n3,gis_phy%n4,                         &
                          gis_phy%lats_nodes_r,gis_phy%global_lats_r,    &
                          gis_phy%lonsperlar,                            &
                          gis_phy%lats_nodes_ext,gis_phy%global_lats_ext,&
                          gis_phy%colat1,gis_phy%idrt)

      idrt = gis_phy%idrt

!     write(0,*)' gis_phy%lonsperlar2c=',gis_phy%lonsperlar

      gis_phy%lats_node_r_max     = lats_node_r_max
      gis_phy%lats_nodes_r_fix(:) = gis_phy%lats_node_r_max 

!* set up scatter_lats and scatter_lons for simple scatter (Sarah Lu)
      allocate(scatter_lats (latr))
      allocate(scatter_lons (latr))
      scatter_lats(1:latr) = gis_phy%global_lats_r(1:latr)
      scatter_lons(1:latr) = gis_phy%lonsperlar(1:latr)

!* change lats_node_r to lats_node_r_max to allow the pointer option
!*    call sfcvar_aldata(lonr, lats_node_r, lsoil, gis_phy%sfc_fld, ierr)
!*    call flxvar_aldata(lonr, lats_node_r, gis_phy%flx_fld, ierr)

      call sfcvar_aldata(lonr, lats_node_r_max, lsoil, gis_phy%sfc_fld, ierr)
      call sfc_init(gis_phy%sfc_fld, ierr)
      call flxvar_aldata(lonr, lats_node_r_max, gis_phy%flx_fld, ierr)

!      print *,' check after sfc flx var_aldata ' 
      IF (me == 0) write(*,*) ' in "GFS_Initialize_ESMFMod,lonr,lats_node_r,&
                   nr_nst,nf_nst : ',lonr,lats_node_r,nr_nst,nf_nst,        &
                  'lats_node_r_max=',lats_node_r_max

!      write(0,*)'in GFS_Initialize_ESMFMod,lonr=',lonr,'lats_node_r=',    &
!                lats_node_r,'nr_nst=',nr_nst,' nf_nst=',nf_nst,           &
!                'lats_node_r_max=',lats_node_r_max

!    Modified by Moorthi

      if (nstf_name(1) > 0) then
        call nstvar_aldata(lonr,lats_node_r_max,gis_phy%nst_fld,ierr)
      endif

!    Add (Xingren Wu)
      if (a2oi_out .or. cplflx) then
        call aoivar_aldata(lonr,lats_node_r_max,gis_phy%aoi_fld,ierr)
      endif

!! allocate grid_fld                      --- Sarah Lu
      gis_phy%grid_aldata = grid_aldata
      if ( gis_phy%grid_aldata ) then
        if( me == 0) print *,'LU_PHY: grid_fld allocated ; copy is used' &
                            ,' gis_phy%ntrac=',gis_phy%ntrac
        call gridvar_aldata (lonr, lats_node_r_max, levs,               &
                             gis_phy%ntrac, gis_phy%grid_fld, ierr)  
      else                                                            
        allocate (gis_phy%grid_fld%rqtk(lonr,lats_node_r_max))
        if( me == 0) print *,'LU_PHY: grid_fld not allocated ; pointer is used'
      endif                                                     

!*    allocate (gis_phy%grid_gr(lonr*lats_node_r_max,lotgr), stat = ierr )
!*    ALLOCATE (gis_phy%CLDCOV(LEVS,LONR,LATS_NODE_R),  stat = ierr)
!     ALLOCATE (gis_phy%CLDCOV(LEVS,LONR,LATS_NODE_R_MAX), stat = ierr)

      ALLOCATE (gis_phy%XLON(LONR,LATS_NODE_R),         stat = ierr)
      ALLOCATE (gis_phy%XLAT(LONR,LATS_NODE_R),         stat = ierr)
      ALLOCATE (gis_phy%COSZDG(LONR,LATS_NODE_R),       stat = ierr)
      ALLOCATE (gis_phy%SFALB(LONR,LATS_NODE_R),        stat = ierr)
      ALLOCATE (gis_phy%HPRIME(NMTVR,LONR,LATS_NODE_R), stat = ierr)
      ALLOCATE (gis_phy%FLUXR(nfxr,LONR,LATS_NODE_R),   stat = ierr)

      ALLOCATE (gis_phy%importData(LONR,LATS_NODE_R,NImportFields),stat = ierr)
      gis_phy%importData = -99999.0

      nblck         = LONR/NGPTC + 1
      gis_phy%NBLCK = nblck

      ALLOCATE (gis_phy%SWH(NGPTC,LEVS,NBLCK,LATS_NODE_R), stat = ierr)
      ALLOCATE (gis_phy%HLW(NGPTC,LEVS,NBLCK,LATS_NODE_R), stat = ierr)
!     if (stochphys) then
      ALLOCATE (gis_phy%SWHC(NGPTC,LEVS,NBLCK,LATS_NODE_R), stat = ierr)
      ALLOCATE (gis_phy%HLWC(NGPTC,LEVS,NBLCK,LATS_NODE_R), stat = ierr)
!     endif

! idea add by hmhj - commented by moorthi since this is unused
!     ALLOCATE (gis_phy%HTRSWB(NGPTC,LEVS,NBDSW,NBLCK,LATS_NODE_R),stat=ierr)
!     ALLOCATE (gis_phy%HTRLWB(NGPTC,LEVS,NBDLW,NBLCK,LATS_NODE_R),stat=ierr)

!
      ALLOCATE (gis_phy%JINDX1(LATS_NODE_R), stat = ierr)
      ALLOCATE (gis_phy%JINDX2(LATS_NODE_R), stat = ierr)
      ALLOCATE (gis_phy%DDY(LATS_NODE_R),    stat = ierr)
!
      allocate (gis_phy%phy_f3d(NGPTC,LEVS,ntot3d,NBLCK,lats_node_r), stat = ierr)
      allocate (gis_phy%phy_f2d(lonr,lats_node_r,ntot2d),             stat = ierr)
      if (cscnv) then
        gis_phy%num_ctp = nctp
      else
        gis_phy%num_ctp = 1
      endif
      allocate (gis_phy%phy_fctd(lonr,lats_node_r,gis_phy%num_ctp),  stat = ierr)
!
      allocate (gis_phy%fhour_idate(1,5), stat = ierr )
!
!      write(0,*) ' check after lots allocates,size(fhour_idate)= ' ,   &
!        size(gis_phy%fhour_idate,1),size(gis_phy%fhour_idate,2),ierr

      if (ldiag3d) then
!* change lats_node_r to lats_node_r_max for consistency
!*      call d3d_init(ngptc,nblck,lonr,lats_node_r,levs,pl_coeff)
        call d3d_init(ngptc,nblck,lonr,lats_node_r_max,levs,pl_coeff)
      endif

     if (lsidea) then
        call jh_integral_init(ngptc,nblck)
      endif

!* allocate g3d_fld and g2d_fld
      if (lgocart) then
        call g3d_aldata (lonr, lats_node_r_max, levs,                   &
                         gis_phy%g3d_fld, ierr)  
!       call g2d_aldata (lonr, lats_node_r_max, gis_phy%gfs_phy_tracer, &
        call g2d_aldata (lonr, lats_node_r_max, gfs_phy_tracer,         &
                         gis_phy%g2d_fld, ierr)

        ngrids_aer = 0
        if ( gis_phy%g2d_fld%du%nfld > 0 )                              &
          ngrids_aer = ngrids_aer+gis_phy%g2d_fld%du%nfld

        if ( gis_phy%g2d_fld%su%nfld > 0 )                              &
          ngrids_aer = ngrids_aer+gis_phy%g2d_fld%su%nfld

        if ( gis_phy%g2d_fld%oc%nfld > 0 )                              &
          ngrids_aer = ngrids_aer+gis_phy%g2d_fld%oc%nfld

        if ( gis_phy%g2d_fld%bc%nfld > 0 )                              &
          ngrids_aer = ngrids_aer+gis_phy%g2d_fld%bc%nfld

        if ( gis_phy%g2d_fld%ss%nfld > 0 )                              &
          ngrids_aer = ngrids_aer+gis_phy%g2d_fld%ss%nfld

        if ( gis_phy%g2d_fld%met%nfld > 0 )                             &
          ngrids_aer = ngrids_aer+gis_phy%g2d_fld%met%nfld

!       print *, 'INIT_g2d_fld ngrids_aer = ',ngrids_aer
      endif

      if (ntot3d > 0) gis_phy%phy_f3d  = 0.0
      if (ntot2d > 0) gis_phy%phy_f2d  = 0.0
      if (cscnv)      gis_phy%phy_fctd = 0.0
!jw
!jw    allocate(buff_mult_piecea2d(lonr,lats_node_r_max,1:ngrids_sfcc2d+ngrids_nst+1))
       allocate(buff_mult_piecea2d(lonr,lats_node_r_max,ngrids_sfcc2d+1))
       allocate(buff_mult_piecea3d(lonr,lats_node_r_max,ngrids_sfcc3d+1))
       allocate(buff_mult_piecef  (lonr,lats_node_r_max,ngrids_flx+1))
       allocate(buff_mult_pieceg  (lonr,lats_node_r_max,ngrids_aer+1))
       if (nstf_name(1) > 0) then
         allocate(buff_mult_piecenst(lonr,lats_node_r_max,ngrids_nst+1))
       endif

!      write(0,*)' lonr=',lonr,' lats_node_r_max=',lats_node_r_max,       &
!      'ngrids_sfcc2d=',ngrids_sfcc2d,' ngrids_nst=',ngrids_nst,          &
!      'ngrids_flx=',ngrids_flx

       buff_mult_piecea2d(1:lonr,1:lats_node_r_max,1:ngrids_sfcc2d+1) = 0.
       buff_mult_piecea3d(1:lonr,1:lats_node_r_max,1:ngrids_sfcc3d+1) = 0.
       buff_mult_piecef(1:lonr,1:lats_node_r_max,1:ngrids_flx+1)      = 0.
       buff_mult_pieceg(1:lonr,1:lats_node_r_max,1:ngrids_aer+1)      = 0.
       if (nstf_name(1) > 0) then
         buff_mult_piecenst(1:lonr,1:lats_node_r_max,1:ngrids_nst+1)    = 0.
       endif
!!
!!
!       write(0,*)' gis_phy%lonsperlar2b=',gis_phy%lonsperlar
!       write(0,*)' before fix_fields'

      PRINT*, 'in phys initialize, lsidea, ipe_to_wam_coupling=', &
            lsidea, ipe_to_wam_coupling
      IF(lsidea .AND. ipe_to_wam_coupling) THEN
        lowst_ipe_level = 80

        IF(.NOT. ASSOCIATED(ZMT))  &
          ALLOCATE(ZMT(lonr, lats_node_r_max, lowst_ipe_level:levs)) 
        IF(.NOT. ASSOCIATED(MMT))  &
          ALLOCATE(MMT(lonr, lats_node_r_max, lowst_ipe_level:levs)) 
        IF(.NOT. ASSOCIATED(JHR))  &
          ALLOCATE(JHR(lonr, lats_node_r_max, lowst_ipe_level:levs)) 
        IF(.NOT. ASSOCIATED(SHR))  &
          ALLOCATE(SHR(lonr, lats_node_r_max, lowst_ipe_level:levs)) 
        IF(.NOT. ASSOCIATED(O2DR)) &
          ALLOCATE(O2DR(lonr, lats_node_r_max, lowst_ipe_level:levs)) 
! Just for simple testing, when run the model, should commenting the lines.
!--------------------------------------------------------------------------
!      ZMT(:,:,81:150)=1.0E-7
!      MMT(:,:,81:150)=1.0E-7
!      JHR(:,:,81:150)=1.0E-7
!      SHR(:,:,81:150)=1.0E-7
!      O2DR(:,:,81:150)=1.0E-7

!      ZMT(:,:,80)=-1.0E30
!      MMT(:,:,80)=-1.0E30
!      JHR(:,:,80)=-1.0E30
!      SHR(:,:,80)=-1.0E30
!      O2DR(:,:,80)=-1.0E30
! End of the simple testing.
!---------------------------

      END IF

      IF(lsidea .AND. gis_phy%restart_run) THEN
        call fix_fields_idea_rst(gis_phy%LONSPERLAR, gis_phy%GLOBAL_LATS_R,    &
                        gis_phy%XLON,       gis_phy%XLAT,   gis_phy%sfc_fld,   &
                        gis_phy%nst_fld,    gis_phy%HPRIME, gis_phy%JINDX1,    &
                        gis_phy%JINDX2,     gis_phy%DDY,    gis_phy%OZPLIN,    &
                        gis_phy%nam_gfs_phy%sfc_ini,                           &
                        gis_phy%nam_gfs_phy%grd_ini,                           &
                        gis_phy%nam_gfs_phy%nst_ini,                           &
                        nblck, gis_phy%phy_f3d, gis_phy%phy_f2d )
      ELSE
        call fix_fields(gis_phy%LONSPERLAR, gis_phy%GLOBAL_LATS_R,             &
                        gis_phy%XLON,       gis_phy%XLAT,   gis_phy%sfc_fld,   &
                        gis_phy%nst_fld,    gis_phy%HPRIME, gis_phy%JINDX1,    &
                        gis_phy%JINDX2,     gis_phy%DDY,    gis_phy%OZPLIN,    &
                        gis_phy%nam_gfs_phy%sfc_ini,                           &
                        gis_phy%nam_gfs_phy%grd_ini,                           &
                        gis_phy%nam_gfs_phy%nst_ini,                           &
                        nblck, gis_phy%phy_f3d, gis_phy%phy_f2d )
      END IF
!      print *,' GISXLAT=',gis_phy%XLAT(1,:)
!       write(0,*)' after fix_fields'
!!
! coord def (lats_node_r, ipt_lats_node_r, and lats_nodes_a_fix)     
      gis_phy%lats_node_r     = lats_node_r                      
      gis_phy%ipt_lats_node_r = ipt_lats_node_r              

!     call mpi_quit(3333)
!!
!! debug print (Sarah Lu)
!      if(me==0) then                                                  
!         do j=1,latr                                                
!         print *, 'PHY: lonsperlar=',j,gis_phy%lonsperlar(j)     
!         enddo                                                        
!         print *, 'PHY: lats_node_r_max=',gis_phy%lats_node_r_max 
!         print *, 'PHY: lats_nodes_r=',gis_phy%lats_nodes_r(:)   
!         print *, 'PHY: global_lats_r=',gis_phy%global_lats_r(:)
!      endif                                                          
!      print *, 'PHY:  lats_node_r=',me,gis_phy%lats_node_r      
!      n = 0                                                     
!      do j = 1, gis_phy%lats_node_r                              
!         ilat = gis_phy%global_lats_r(gis_phy%ipt_lats_node_r-1+j)  
!         l =  gis_phy%lonsperlar(ilat)                            
!         n = n + l                                                 
!         print *, 'PHY:, xlat=',me,n,gis_phy%ipt_lats_node_r-1+j,& 
!                   j, ilat, l, 57.29578*gis_phy%xlat(1,j)       
!      enddo                                             
!!
!!

!  ---  set up sigma levels before radiation initialization

        allocate (si_loc(levr+1))
        if ( hybrid .and. (.not.gen_coord_hybrid) ) then

!  ---  following sela: si(k)=(ak5(k)+bk5(k)*Typical_pgr)/Typical_pgr
!       ak(k) and bk(k) go from top to bottom, adjust si_loc direction
!       according model direction flag iflip.

!         if (iflip == 1) then        ! vertical from sfc to top
!*********************************************************************
!!!! COMMENTED WHILE TESTING OTHER PARTS OF CODE !!!!!!!!!!!!!!!!!!!!!!
!           si_loc(levr+1) = ak5(1)/typical_pgr+bk5(1)
            si_loc(levr+1) = 0.0
            if(nemsio_in) then
              do k = 1, levr
                si_loc(levr+1-k) = ak5(levp1-levr+k)/typical_pgr          &
     &                           + bk5(levp1-levr+k)
              enddo
            else
              do k = 1, levr
                si_loc(k) = float(levr+1-k) / float(levr)
              enddo
            endif
            if (me == 0) write(0,*)' si_loc=',si_loc(1:levr+1)
!*********************************************************************
!         else                        ! vertical from top to sfc
!           si_loc(1) = ak5(1)/typical_pgr+bk5(1)
!           do k = 1, levr
!             si_loc(k+1)      = ak5(levp1-levr+k)/typical_pgr
!             &
!    &                         + bk5(levp1-levr+k)
!           enddo
!         endif
        else
!  ---  using the model's native sigama coordinate, si
!       si(k) goes from bottom to top, adjust si_loc direction
!       according model direction flag iflip.

!         if (iflip == 1) then        ! vertical from sfc to top
            si = 0.
            do k = 1, levr
              si_loc(k) = si(k)
            enddo
            si_loc(levr+1) = si(levp1)
!         else                        ! vertical from top to sfc
!           do k = 1, levr
!             si_loc(k+1) = si(levr-k+1)
!           enddo
!           si_loc(1) = si(levp1)
!         endif
        endif    ! end_if_hybrid
!     This is called in nuopc_phys_init, don't run twice if using the wrapper
      if ( .not. use_nuopc ) then
        call rad_initialize                                             &
!  ---  inputs:
     &     ( si_loc,levr,ictm,isol,ico2,iaer,ialb,iems,ntcw,            &
     &       num_p3d,npdf3d,ntoz,iovr_sw,iovr_lw,isubc_sw,isubc_lw,     &
     &       crick_proof,ccnorm,norad_precip,idate,iflip,me )
      endif
!  ---  outputs: ( none )
!
!     zero fluxes and diagnostics
!
      gis_phy%zhour = gis_phy%phour

      if (me == 0)                                                      &
      write(0,*)' in physics initialize phour=',gis_phy%phour,          &
                ' fhour=',fhour,' zhour=',gis_phy%zhour

      gis_phy%FLUXR = 0.
!
      call flx_init(gis_phy%flx_fld, ierr)

!    Add (Xingren Wu)
      if (a2oi_out .or. cplflx) then
        call aoivar_init(gis_phy%aoi_fld, ierr)
      endif

      if (lsidea) then
        call jh_integral_zero
      endif

      if (ldiag3d) then
        call d3d_zero
      endif

      if (use_nuopc) then
        associate(mdl_parm => gis_phy%mdl_parm)

        call nuopc_phys_init ( &
       mdl_parm, ntcw, ncld, ntoz, NTRAC, levs, me, lsoil, &
       lsm, nmtvr, nrcm, levozp, lonr, latr, jcap, num_p3d, num_p2d, &
       npdf3d, ncnvcld3d, pl_coeff, ncw, crtrh, cdmbgwd, ccwf, dlqf,ctei_rm, &
       cgwf, prslrd0, ral_ts, ras, pre_rad, ldiag3d, lgocart,cplflx,flipv, &
       old_monin,cnvgwd, shal_cnv, imfshalcnv,imfdeepcnv, cal_pre, mom4ice, &
       mstrat,trans_trac,nstf_name,moist_adj,thermodyn_id,sfcpress_id, &
       gen_coord_hybrid, levr, lsidea, pdfcld, shcnvcw, redrag, &
       hybedmf, dspheat, dxmax, dxmin, dxinv, &
       cscnv, gis_phy%num_ctp, ntke, do_shoc, shocaftcnv, ntot3d, ntot2d, &
       si_loc, ictm, isol, ico2, iaer, ialb, iems, &
       iovr_sw,iovr_lw,isubc_sw,isubc_lw, shoc_cld, &
       crick_proof,ccnorm,norad_precip,idate,iflip, nlunit)

       call mdl_parm%print("In gfs_physics_initialize")
        end associate
      end if

!
! initialize start_step (Sarah Lu)

      gis_phy% start_step  = .true.

      rc = 0
!
      if (allocated(si_loc)) deallocate (si_loc)
      phys_ini_time = phys_ini_time + (timef() - btime)
      write(0,*)' phys_ini_time=',phys_ini_time*1.0e-3,' me=',me
!
!       write(0,*)' gis_phy%lonsperlar2=',gis_phy%lonsperlar
!       write(0,*)' returning from gfs_physics_initialize'
!
      end subroutine gfs_physics_initialize
!
! =========================================================================
!
      end module gfs_physics_initialize_mod
