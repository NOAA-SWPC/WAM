
! !module: gfs_dynamics_initialize_slg_mod 
!          --- initialize module of the
!              gridded component of the gfs dynamics system.
!              gfs dynamics related initilization 
!
! !description: gfs dynamics gridded component initialize module.
!
! !revision history:
!
!  november 2004  weiyu yang     initial code.
!  january 2006  s. moorthi      update to the new gfs version
!  august 2006   h. juang        add option to run generalized coordinates
!  december 2006 s. moorthi      gfsio included
!  january 2007 h. juang         change for dynamics only
!  May     2008 j. wang          change for gfs wrt grid component
!  Oct 04  2009 sarah lu         init xlon, xlat, lats_nodes_a_fix
!  Oct 05  2009 sarah lu         grid_gr unfolded from 2D to 3D
!  Oct 16  2009 sarah lu         initialize gfs_dyn_tracer
!  november 2009 j. wang         grid_gr_dfi for digital filter
!  Feb 05  2010 j. wang          add option to read in  restart file
!  Aug 19  2010 S. Moorthi       Updated for T574 + added num_reduce to namelist
!  Aug 25  2010 sarah lu         add option to compute tracer global sum
!  Sep 08  2010 J. Wang          changed gfsio file to nemsio file
!  Nov 01  2010 H. Juang         add non-iteration dimensional-split
!                                semi-Lagrangian dynamics (ndslfv)
!                                with mass_dp, process_split options.
!  Dec 16  2010 J. Wang          changed to nemsio library
!  Feb 20  2011 H. Juang         implement into nems for mass_dp and ndsl
!  Feb 28  2011 Sarah Lu         add thermodyn_id and sfcpress_id
!  Apr 06  2012 H. Juang         add idea
!  Sep 20  2012 J. Wang          add sigio option
!  Feb 04  2013 W. Yang          modified for the slg version.
!
!  Jun 26  2014 S. Moorthi       Modified to read lonsperlat from a file
!  Jul 21  2014 S. Moorthi       removed num_reduce
!  feb 14  2015 J. Wang          add option for read in nemsio
!
! !interface:
!
      module gfs_dynamics_initialize_slg_mod
!
!!uses:
!
      use gfs_dynamics_getcf_mod
      use gfs_dyn_machine, only : kind_io4
      use nemsio_module ,  only : nemsio_init
      use gfs_dyn_iau_module, only: init_iau
!
      use gfs_dyn_write_state,   only : buff_mult_pieceg
      use gfs_dyn_layout1,       only : ipt_lats_node_a, lats_node_a_max,lon_dim_a
      use gfs_dyn_resol_def,     only : adiabatic, thermodyn_id, sfcpress_id
      use namelist_dynamics_def, only : fhrot,fhini,nemsio_in, redgg_a, semilag, &
                                        do_sppt,do_skeb,do_shum,do_vc,gg_tracers
      use gfs_dyn_tracer_config, only: gfs_dyn_tracer, tracer_config_init,gfs_dyn_tracer
      use gfs_dyn_io_header,     only: z_r,z
!     use gfs_dyn_io_header,     only: z_r,z,gz_grid
      use gfs_dyn_coordinate_def
      use layout_lag          ,  only : lat1s_h,lats_dim_h,lon_dim_h
      use gfs_dyn_gg_def      ,  only : lats_nodes_h,global_lats_h
      use layout_grid_tracers ,  only : rgt_a,rgt_h,xhalo,yhalo
#ifndef IBM
      USE omp_lib
#endif

      implicit none

      contains

      subroutine gfs_dynamics_initialize_slg(gis_dyn, rc)

! this subroutine set up the internal state variables,
! allocate internal state arrays for initializing the gfs system.
!----------------------------------------------------------------
!
      implicit none
!
      integer, parameter :: iunit=101
      type(gfs_dynamics_internal_state), pointer, intent(inout) :: gis_dyn
      integer,                                    intent(out)   :: rc

      character*20         cfile, cfile2
      logical           :: file_exists=.false.
      integer           :: ierr, npe_single_member, n1, n2, latghf, iret

      integer           :: num_parthds
      integer           :: i, j, k, l, n, locl

      real(kind=8)      :: dyn_ini_time=0.0, btime, timef
!-------------------------------------------------------------------

! set up gfs internal state dimension and values for dynamics etc
!-------------------------------------------------------------------
      btime  = timef()
      me     = gis_dyn%me
      if (me == 0) write(0,*)'in initial,nbefore allocate lonsperlat,',&
                   allocated(gis_dyn%lonsperlat),'latg=',latg
!
      nodes  = gis_dyn%nodes
      nlunit = gis_dyn%nam_gfs_dyn%nlunit
      npe_single_member = gis_dyn%npe_single_member

      semilag = gis_dyn%SLG_FLAG
      call compns_dynamics(gis_dyn%deltim, gis_dyn%iret, gis_dyn%ntrac,	&
                           gis_dyn%nxpt,   gis_dyn%nypt, gis_dyn%jintmx,&
                           gis_dyn%jcap,                              	&
                           gis_dyn%levs,   gis_dyn%levr, 		&
                           gis_dyn%lonf,   gis_dyn%latg,          	&
                           gis_dyn%ntoz,   gis_dyn%ntcw, gis_dyn%ncld, 	&
                           gis_dyn%ntke,   gis_dyn%spectral_loop,      	&
                           me, gis_dyn%thermodyn_id,gis_dyn%sfcpress_id,&
                           gis_dyn%nam_gfs_dyn%nlunit, 			&
                           gis_dyn%nam_gfs_dyn%gfs_dyn_namelist,        &
                           gis_dyn%ndfi)                           ! jw
!
      call get_tracer_const(gis_dyn%ntrac,me,gis_dyn%nam_gfs_dyn%nlunit)
!
! met+chem tracer specification (Sarah Lu)
!
!     call tracer_config_init( gis_dyn%gfs_dyn_tracer, gis_dyn%ntrac,     &
      call tracer_config_init( gis_dyn%ntrac, gis_dyn%ntoz, gis_dyn%ntcw, &
                               gis_dyn%ncld,  gis_dyn%ntke, me )
!      gfs_dyn_tracer = gis_dyn%gfs_dyn_tracer
      if( me == 0) then
        write(0,*)'LU_TRC, exit tracer_config_init in dyn'
        write(0,*)'LU_TRC, ntrac=     ',gfs_dyn_tracer%ntrac,ntrac
        write(0,*)'LU_TRC, ntrac_met =',gfs_dyn_tracer%ntrac_met
        write(0,*)'LU_TRC, ntrac_chem=',gfs_dyn_tracer%ntrac_chem
        do n = 1, gfs_dyn_tracer%ntrac
          write(0,*)'LU_TRC, tracer_vname=',gfs_dyn_tracer%vname(n, :)
        enddo
      endif
!
      ntrac   = gis_dyn%ntrac
      nxpt    = gis_dyn%nxpt
      nypt    = gis_dyn%nypt
      jintmx  = gis_dyn%jintmx
      jcap    = gis_dyn%jcap
      levs    = gis_dyn%levs
      levr    = gis_dyn%levr
      lonf    = gis_dyn%lonf
      latg    = gis_dyn%latg
      ntoz    = gis_dyn%ntoz
      ntcw    = gis_dyn%ntcw
      ncld    = gis_dyn%ncld
      thermodyn_id = gis_dyn%thermodyn_id
      sfcpress_id  = gis_dyn%sfcpress_id
      nemsio_in    = gis_dyn%nemsio_in
      lon_dim_a = lonf + 2

      if (gis_dyn%nam_gfs_dyn%total_member <= 1) then
        ens_nam=' '
      else
        write(ens_nam,'("_",i2.2)') gis_dyn%nam_gfs_dyn%member_id
      endif
!
      levh   = ntrac*levs
!M    latgd  = latg+ 2*jintmx 
      jcap1  = jcap+1 
      jcap2  = jcap+2 
      latg2  = latg/2 
      levm1  = levs-1 
      levp1  = levs+1 
!M    lonfx  = lonf + 1 + 2*nxpt+1 
      lnt    = jcap2*jcap1/2 
      lnuv   = jcap2*jcap1 
      lnt2   = 2*lnt 
      lnt22  = 2*lnt+1 
      lnte   = (jcap2/2)*((jcap2/2)+1)-1 
      lnto   = (jcap2/2)*((jcap2/2)+1)-(jcap2/2) 
      lnted  = lnte 
      lntod  = lnto 

!jw   ngrids_gg       = 2+levs*(4+ntrac)
      ngrids_gg       = 2+levs*(5+ntrac)
      gis_dyn%lnt2    = lnt2

      allocate(lat1s_a(0:jcap))
      allocate(lon_dims_a(latg))
!M    allocate(lon_dims_a(latgd))
!M    allocate(lon_dims_ext(latgd))

      allocate(colrad_a(latg2))
      allocate(wgt_a(latg2))
      allocate(wgtcs_a(latg2))
      allocate(rcs2_a(latg2))
      allocate(sinlat_a(latg))
      allocate(coslat_a(latg))
      coslat_a = 0

      allocate(am(levs,levs))
      allocate(bm(levs,levs))
      allocate(cm(levs,levs))
      allocate(dm(levs,levs,jcap1))
      allocate(tor(levs))
      allocate(si(levp1))
      allocate(sik(levp1))
      allocate(sl(levs))
      allocate(slk(levs))
      allocate(del(levs))
      allocate(rdel2(levs))
      allocate(ci(levp1))
      allocate(cl(levs))
      allocate(tov(levs))
      allocate(sv(levs))
      allocate(tor_slg(LEVS))
      allocate(y_ecm(LEVS,LEVS))
      allocate(t_ecm(LEVS,LEVS))
      allocate(am_slg(LEVS,LEVS))
      allocate(bm_slg(LEVS,LEVS))
      allocate(sv_ecm(LEVS))
      allocate(sv_slg(LEVS))
      allocate(D_slg_m(levs,levs,jcap1))

      allocate(ak5(levp1))
      allocate(bk5(levp1))
      allocate(ck5(levp1)) 
      allocate(thref(levp1))
      allocate(ck(levs))
      allocate(dbk(levs))
      allocate(bkl(levs))
      allocate(amhyb(levs,levs))
      allocate(bmhyb(levs,levs))
      allocate(smhyb(levs,levs))
      allocate(hmhyb(levs,levs))
      allocate(svhyb(levs))
      allocate(tor_hyb(levs))
      allocate(d_hyb_m(levs,levs,jcap1))
      allocate(dm205_hyb(jcap1,levs,levs))

      allocate(spdmax(levs))

      if (me == 0)                                                       &
      write(0,*)'before allocate lonsperlat,',allocated(gis_dyn%lonsperlat),'latg=',latg

      allocate(gis_dyn%lonsperlat(latg))

      inquire (file="lonsperlat.dat", exist=file_exists)
      if ( .not. file_exists ) then
        if ( me == 0 ) then
          write(0,*)' Requested lonsperlat.dat  data file does not exist'
          write(0,*)'   *** Stopped in subroutine GFS_Init !!'
        endif
        call mpi_quit(1111)
      else
        open (iunit,file='lonsperlat.dat',status='old',form='formatted',      &
                                          action='read',iostat=iret)
        if (iret /= 0) then
          write(0,*)' iret while reading lonsperlat.dat ',iret
          call mpi_quit(1112)
        endif
        rewind iunit
        read (iunit,*,iostat=iret) latghf,(gis_dyn%lonsperlat(i),i=1,latghf)
        if (latghf+latghf /= latg) then
           write(0,*)' ltghf=',latghf,' not equal to latg/2=',latg/2
           call mpi_quit(1113)
        endif
        do i=1,latghf
          gis_dyn%lonsperlat(latg-i+1) = gis_dyn%lonsperlat(i)
        enddo
        close(iunit)
      endif
!
! spectral location
      P_GZ  = 0*LEVS+0*LEVH+1  !      GZE/O(LNTE/OD,2),
      P_ZEM = 0*LEVS+0*LEVH+2  !     ZEME/O(LNTE/OD,2,LEVS),
      P_DIM = 1*LEVS+0*LEVH+2  !     DIME/O(LNTE/OD,2,LEVS),
      P_TEM = 2*LEVS+0*LEVH+2  !     TEME/O(LNTE/OD,2,LEVS),
      P_QM  = 3*LEVS+0*LEVH+2  !      QME/O(LNTE/OD,2),
      P_ZE  = 3*LEVS+0*LEVH+3  !      ZEE/O(LNTE/OD,2,LEVS),
      P_DI  = 4*LEVS+0*LEVH+3  !      DIE/O(LNTE/OD,2,LEVS),
      P_TE  = 5*LEVS+0*LEVH+3  !      TEE/O(LNTE/OD,2,LEVS),
      P_Q   = 6*LEVS+0*LEVH+3  !       QE/O(LNTE/OD,2),
      P_DLAM= 6*LEVS+0*LEVH+4  !  DPDLAME/O(LNTE/OD,2),
      P_DPHI= 6*LEVS+0*LEVH+5  !  DPDPHIE/O(LNTE/OD,2),
      P_ULN = 6*LEVS+0*LEVH+6  !     ULNE/O(LNTE/OD,2,LEVS),
      P_VLN = 7*LEVS+0*LEVH+6  !     VLNE/O(LNTE/OD,2,LEVS),
      P_W   = 8*LEVS+0*LEVH+6  !       WE/O(LNTE/OD,2,LEVS),
      P_X   = 9*LEVS+0*LEVH+6  !       XE/O(LNTE/OD,2,LEVS),
      P_Y   =10*LEVS+0*LEVH+6  !       YE/O(LNTE/OD,2,LEVS),
      P_ZQ  =11*LEVS+0*LEVH+6  !      ZQE/O(LNTE/OD,2)
      P_RT  =11*LEVS+0*LEVH+7  !      RTE/O(LNTE/OD,2,LEVH),
      P_RM  =11*LEVS+1*LEVH+7  !      RME/O(LNTE/OD,2,LEVH),
      P_RQ  =11*LEVS+2*LEVH+7  !      RQE/O(LNTE/OD,2,LEVH),

      lotls  = 11*LEVS+3*LEVH+6

      g_gz   = 1
      g_uum  = g_gz  + 1        !  for grid point
      g_vvm  = g_uum + levs     !  for grid point
      g_ttm  = g_vvm + levs     !  for grid point
!     g_qm   = g_ttm + levs     !  for grid point
      g_rm   = g_ttm + levs     !  for grid point
      g_dpm  = g_rm  + levh     !  for grid point
      g_qm   = g_dpm + levs     !  for grid point

      g_uu   = g_qm  + 1        !  for grid point
      g_vv   = g_uu  + levs     !  for grid point
      g_tt   = g_vv  + levs     !  for grid point
!     g_q    = g_tt  + levs     !  for grid point

!     g_rm   = g_ttm + levs     !  for grid point
!     g_dpm  = g_rm  + levh     !  for grid point

      g_rq   = g_tt  + levs     !  for grid point
      g_dp   = g_rq  + levh     !  for grid point
      g_q    = g_dp  + levs     !  for grid point

      g_u    = g_q   + 1        !  for grid point
      g_v    = g_u   + levs     !  for grid point
      g_t    = g_v   + levs     !  for grid point
      g_rt   = g_t   + levs     !  for grid point
      g_dpn  = g_rt  + levh     !  for grid point
      g_zq   = g_dpn + levs     !  for grid point

      g_p    = g_zq  + 1        !  for grid point 
      g_dpdt = g_p   + levs     !  for grid point 
      g_zz   = g_dpdt+ levs     !  for grid point

      g_uup   = g_zz  + levs    !  for grid point
      g_vvp   = g_uup + levs    !  for grid point
      g_ttp   = g_vvp + levs    !  for grid point
      g_rqp   = g_ttp + levs    !  for grid point
      g_dpp   = g_rqp + levh    !  for grid point
      g_zqp   = g_dpp + levs    !  for grid point

      g_rqtk  = g_zqp + 1

!     lotgr  = g_zqp
      lotgr  = g_rqtk                ! Added by Moorthi
      lotgr6 = 4*levs + 1*levh + 1

!     write(0,*)' g_rqtk=',g_rqtk
!
        lots     = 6*levs + 1*levh + 5
        lots_slg = 8*levs + 1*levh + 4
        lotd     = 6*levs + 2*levh+0
!       lota     = 5*levs + 1*levh+2
        lota     = 5*levs + 1*levh + 3
        lotp     = 4*levs
!
        ksz      = 1
        ksd      = ksz  + levs
        kst      = ksd  + levs
        ksr      = kst  + levs
        ksdp     = ksr  + levh
        ksq      = ksdp + levs
        ksplam   = ksq+1
        kspphi   = ksplam + 1
        ksu      = kspphi + 1
        ksv      = ksu    + levs
        kzslam   = ksv    + levs
        kzsphi   = kzslam + 1
!
        kau      = 1
        kav      = kau  + levs
        kat      = kav  + levs
        kar      = kat  + levs
        kadp     = kar  + levh
        kaps     = kadp + levs
        kazs     = kaps + 1
        kap2     = kazs + 1
!
        kdpphi   = 1
        kzzphi   = kdpphi + levs
        kdplam   = kzzphi + levs
        kzzlam   = kdplam + levs
!
        kdtphi   = 1
        kdrphi   = kdtphi + levs
        kdtlam   = kdrphi + levh
        kdrlam   = kdtlam + levs
        kdulam   = kdrlam + levh
        kdvlam   = kdulam + levs
        kduphi   = kdvlam + levs
        kdvphi   = kduphi + levs
!
! point to internal state
        gis_dyn%p_zem   = p_zem            !     zeme/o(lnte/od,2,levs),
        gis_dyn%p_dim   = p_dim            !     dime/o(lnte/od,2,levs),
        gis_dyn%p_tem   = p_tem            !     teme/o(lnte/od,2,levs),
        gis_dyn%p_rm    = p_rm             !      rme/o(lnte/od,2,levh),
        gis_dyn%p_qm    = p_qm             !      qme/o(lnte/od,2),
        gis_dyn%p_ze    = p_ze             !      zee/o(lnte/od,2,levs),
        gis_dyn%p_di    = p_di             !      die/o(lnte/od,2,levs),
        gis_dyn%p_te    = p_te             !      tee/o(lnte/od,2,levs),
        gis_dyn%p_rq    = p_rq             !      rqe/o(lnte/od,2,levh),
        gis_dyn%p_q     = p_q              !       qe/o(lnte/od,2),
        gis_dyn%p_dlam  = p_dlam           !  dpdlame/o(lnte/od,2),
        gis_dyn%p_dphi  = p_dphi           !  dpdphie/o(lnte/od,2),
        gis_dyn%p_uln   = p_uln            !     ulne/o(lnte/od,2,levs),
        gis_dyn%p_vln   = p_vln            !     vlne/o(lnte/od,2,levs),
        gis_dyn%p_w     = p_w              !       we/o(lnte/od,2,levs),
        gis_dyn%p_x     = p_x              !       xe/o(lnte/od,2,levs),
        gis_dyn%p_y     = p_y              !       ye/o(lnte/od,2,levs),
        gis_dyn%p_rt    = p_rt             !      rte/o(lnte/od,2,levh),
        gis_dyn%p_zq    = p_zq             !      zqe/o(lnte/od,2)
        gis_dyn%p_gz    = p_gz             !      gze/o(lnte/od,2),
!
        gis_dyn%g_gz    = g_gz
        gis_dyn%g_uum   = g_uum
        gis_dyn%g_vvm   = g_vvm
        gis_dyn%g_ttm   = g_ttm
        gis_dyn%g_rm    = g_rm
        gis_dyn%g_dpm   = g_dpm
        gis_dyn%g_qm    = g_qm
        gis_dyn%g_uu    = g_uu
        gis_dyn%g_vv    = g_vv
        gis_dyn%g_tt    = g_tt
        gis_dyn%g_rq    = g_rq
        gis_dyn%g_dp    = g_dp
        gis_dyn%g_q     = g_q
        gis_dyn%g_u     = g_u
        gis_dyn%g_v     = g_v
        gis_dyn%g_t     = g_t
        gis_dyn%g_rt    = g_rt
        gis_dyn%g_dpn   = g_dpn
        gis_dyn%g_zq    = g_zq
        gis_dyn%g_p     = g_p
        gis_dyn%g_dpdt  = g_dpdt
        gis_dyn%g_zz    = g_zz
        gis_dyn%g_uup   = g_uup
        gis_dyn%g_vvp   = g_vvp
        gis_dyn%g_ttp   = g_ttp
        gis_dyn%g_rqp    = g_rqp
        gis_dyn%g_dpp    = g_dpp
        gis_dyn%g_zqp    = g_zqp
        gis_dyn%g_rqtk   = g_rqtk
!
        gis_dyn%lotls    = lotls
        gis_dyn%lotgr    = lotgr
        gis_dyn%lots     = lots 
        gis_dyn%lots_slg = lots_slg
        gis_dyn%lotd     = lotd
        gis_dyn%lota     = lota
        gis_dyn%lotp     = lotp
!
        allocate(gis_dyn%tee1(levs))

!       print *,' finish dimension in gfs_dynamics_initialize '

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!!
!!      create io communicator and comp communicator
!!
      nodes_comp=nodes
!
!$$$      time0=timer()
!
      if (me == 0) then
!       print 100, jcap,levs
!100   format (' smf ',i3,i3,' created august 2000 ev od ri ')
        print*,'number of threads is',num_parthds()
        print*,'number of mpi procs is',nodes
      endif
!
      gis_dyn%cons0    =    0.0d0
      gis_dyn%cons0p5  =    0.5d0
      gis_dyn%cons1200 = 1200.d0
      gis_dyn%cons3600 = 3600.d0
!
      ls_dim = (jcap1-1)/nodes+1
!!
!cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
      if (me == 0)                                                       &
       write(0,*)'before allocate ls_nodes,',allocated(gis_dyn%ls_nodes),&
       'ls_dim=', ls_dim,'nodes=',nodes
!
      allocate (      gis_dyn%ls_node (ls_dim*3) )
      allocate (      gis_dyn%ls_nodes(ls_dim,nodes) )
      allocate (  gis_dyn%max_ls_nodes(nodes) )
!
      allocate (  gis_dyn%lats_nodes_a_fix(nodes))     ! added for mGrid
!
      allocate (  gis_dyn%lats_nodes_a(nodes) )
      allocate ( gis_dyn%global_lats_a(latg) )
!
      allocate (   gis_dyn%lats_nodes_ext(nodes) )
      allocate ( gis_dyn%global_lats_ext(latg+2*jintmx+2*nypt*(nodes-1)) )

! For creating the ESMF interface state with the GFS
! internal parallel structure.   Weiyu.
!---------------------------------------------------
      ALLOCATE(gis_dyn%TRIE_LS_SIZE      (npe_single_member))
      ALLOCATE(gis_dyn%TRIO_LS_SIZE      (npe_single_member))
      ALLOCATE(gis_dyn%TRIEO_LS_SIZE     (npe_single_member))
      ALLOCATE(gis_dyn%LS_MAX_NODE_GLOBAL(npe_single_member))
      ALLOCATE(gis_dyn%LS_NODE_GLOBAL    (LS_DIM*3, npe_single_member))

      gis_dyn%LS_NODE_GLOBAL     = 0
      gis_dyn%LS_MAX_NODE_GLOBAL = 0
      gis_dyn%TRIEO_TOTAL_SIZE   = 0

      DO i = 1, npe_single_member
          CALL GET_LS_NODE(i-1, gis_dyn%LS_NODE_GLOBAL(1, i),               &
                            gis_dyn%LS_MAX_NODE_GLOBAL(i), gis_dyn%IPRINT)
          gis_dyn%TRIE_LS_SIZE(i) = 0
          gis_dyn%TRIO_LS_SIZE(i) = 0
          DO LOCL = 1, gis_dyn%LS_MAX_NODE_GLOBAL(i)
              gis_dyn%LS_NODE_GLOBAL(LOCL+  LS_DIM, i)   = gis_dyn%TRIE_LS_SIZE(i)
              gis_dyn%LS_NODE_GLOBAL(LOCL+  2*LS_DIM, i) = gis_dyn%TRIO_LS_SIZE(i)

              L = gis_dyn%LS_NODE_GLOBAL(LOCL, i)

              gis_dyn%TRIE_LS_SIZE(i) = gis_dyn%TRIE_LS_SIZE(i) + (JCAP+3-L)/2
              gis_dyn%TRIO_LS_SIZE(i) = gis_dyn%TRIO_LS_SIZE(i) + (JCAP+2-L)/2
          END DO
          gis_dyn%TRIEO_LS_SIZE(i) = gis_dyn%TRIE_LS_SIZE(i)  + gis_dyn%TRIO_LS_SIZE(i) + 3
          gis_dyn%TRIEO_TOTAL_SIZE = gis_dyn%TRIEO_TOTAL_SIZE + gis_dyn%TRIEO_LS_SIZE(i)
      END DO


!---------------------------------------------------
!
      gis_dyn%iprint = 0
      call get_ls_node( me, gis_dyn%ls_node, ls_max_node, gis_dyn%iprint )
!
!
      len_trie_ls = 0
      len_trio_ls = 0
      do locl=1,ls_max_node
         gis_dyn%ls_node(locl+  ls_dim) = len_trie_ls
         gis_dyn%ls_node(locl+2*ls_dim) = len_trio_ls
         l = gis_dyn%ls_node(locl)
         len_trie_ls = len_trie_ls+(jcap+3-l)/2
         len_trio_ls = len_trio_ls+(jcap+2-l)/2
      enddo
      if (me == 0) print *,'ls_node=',gis_dyn%ls_node(1:ls_dim),'2dim=',  &
         gis_dyn%ls_node(ls_dim+1:2*ls_dim),'3dim=',  &
         gis_dyn%ls_node(2*ls_dim+1:3*ls_dim)
!
!
      allocate ( gis_dyn%epse  (len_trie_ls) )
      allocate ( gis_dyn%epso  (len_trio_ls) )
      allocate ( gis_dyn%epsedn(len_trie_ls) )
      allocate ( gis_dyn%epsodn(len_trio_ls) )
!
      allocate ( gis_dyn%snnp1ev(len_trie_ls) )
      allocate ( gis_dyn%snnp1od(len_trio_ls) )
!
      allocate ( gis_dyn%ndexev(len_trie_ls) )
      allocate ( gis_dyn%ndexod(len_trio_ls) )
!
      allocate ( gis_dyn%plnev_a(len_trie_ls,latg2) )
      allocate ( gis_dyn%plnod_a(len_trio_ls,latg2) )
      allocate ( gis_dyn%pddev_a(len_trie_ls,latg2) )
      allocate ( gis_dyn%pddod_a(len_trio_ls,latg2) )
      allocate ( gis_dyn%plnew_a(len_trie_ls,latg2) )
      allocate ( gis_dyn%plnow_a(len_trio_ls,latg2) )
!
      gis_dyn%maxstp = 36

 
      if(me == 0) print*,'from compns_dynamics : iret=',gis_dyn%iret		&
                        ,' nsout=',nsout,' nsres=',nsres
      if(gis_dyn%iret/=0) then
        if(me == 0) print *,' incompatible namelist - aborted in main'
        call mpi_quit(13)
      endif
!!
      gis_dyn%lats_nodes_ext = 0
      call getcon_dynamics(gis_dyn%n3,              gis_dyn%n4,			&
                           gis_dyn%ls_node,         gis_dyn%ls_nodes,           &
                           gis_dyn%max_ls_nodes,    gis_dyn%lats_nodes_a,       &
                           gis_dyn%global_lats_a,   gis_dyn%lonsperlat,         &
                           gis_dyn%lats_node_a_max, gis_dyn%lats_nodes_ext,     &
                           gis_dyn%global_lats_ext, gis_dyn%epse,               &
                           gis_dyn%epso,            gis_dyn%epsedn,             &
                           gis_dyn%epsodn,          gis_dyn%snnp1ev,            &
                           gis_dyn%snnp1od,         gis_dyn%ndexev,             &
                           gis_dyn%ndexod,          gis_dyn%plnev_a,            &
                           gis_dyn%plnod_a,         gis_dyn%pddev_a,            &
                           gis_dyn%pddod_a,         gis_dyn%plnew_a,            &
                           gis_dyn%plnow_a,         gis_dyn%colat1)
!
      gis_dyn%lats_node_a     = gis_dyn%lats_nodes_a(me+1)
      gis_dyn%ipt_lats_node_a = ipt_lats_node_a

      if (me == 0)                                                        &
       write(0,*)'after getcon_dynamics,lats_node_a=',gis_dyn%lats_node_a &
         ,'ipt_lats_node_a=',gis_dyn%ipt_lats_node_a,'ngptc=',ngptc
!
!     if (gg_tracers) then
        if (.not. allocated(lats_nodes_h))  allocate (lats_nodes_h(nodes))
        if (.not. allocated(lat1s_h))       allocate (lat1s_h(0:jcap))
        if (.not. allocated(global_lats_h)) allocate (global_lats_h(latg+2*yhalo*nodes))

        call getcon_lag(gis_dyn%lats_nodes_a,gis_dyn%global_lats_a,        &
                        lats_nodes_h, global_lats_h,                       &
                        gis_dyn%lonsperlat,xhalo,yhalo)

        if (.not.allocated(rgt_a))          allocate (rgt_a(lonf,levs,lats_dim_a,ntrac))
        if (.not.allocated(rgt_h))          allocate (rgt_h(lon_dim_h,levs,lats_dim_h,ntrac))
        rgt_h = 0.0
        rgt_a = 0.0
        if (me == 0) write(0,*)' in initialize lon_dim_h-',lon_dim_h,' levs=',levs  &
                  ,' lats_dim_h=',lats_dim_h,' ntrac=',ntrac, ' size_rgt=',size(rgt_h,dim=3)
!     endif
!
!     print *,'ls_nodes=',gis_dyn%ls_nodes(1:ls_dim,me+1)
!!
      gis_dyn%nblck = lonf/ngptc + 1

!
! initialize coord def (xlon,xlat) and lats_nodes_a_fix
!
      gis_dyn%lats_nodes_a_fix(:) = gis_dyn%lats_node_a_max
      allocate ( gis_dyn%XLON(lonf,gis_dyn%lats_node_a) )
      allocate ( gis_dyn%XLAT(lonf,gis_dyn%lats_node_a) )

      call gfs_dyn_lonlat_para(gis_dyn%global_lats_a,        &
              gis_dyn%xlon, gis_dyn%xlat, gis_dyn%lonsperlat)

      allocate ( gis_dyn%trie_ls (len_trie_ls,2,lotls) )
      allocate ( gis_dyn%trio_ls (len_trio_ls,2,lotls) )
      allocate ( gis_dyn%grid_gr (lonf,lats_node_a_max,lotgr    ) )
      allocate ( gis_dyn%grid_gr6(lonf,lats_node_a_max,lotgr6*2 ) )
      allocate ( gis_dyn%pwat    (lonf,lats_node_a) )
      allocate ( gis_dyn%ptot    (lonf,lats_node_a) )
      allocate ( gis_dyn%ptrc    (lonf,lats_node_a,ntrac) )         !glbsum
      allocate ( z(lnt2) )
      allocate ( z_r(lnt2) )
!     allocate ( gz_grid(lon_dim_a,lats_dim_a) )
!
!     allocate ( gis_dyn%syn_ls_a(4*ls_dim,gis_dyn%lots,latg2) )
!     allocate ( gis_dyn%dyn_ls_a(4*ls_dim,gis_dyn%lotd,latg2) )
!
!     allocate ( gis_dyn%syn_gr_a_1(lonfx*gis_dyn%lots,lats_dim_ext) )
!     allocate ( gis_dyn%syn_gr_a_2(lonfx*gis_dyn%lots,lats_dim_ext) )
!     allocate ( gis_dyn%pyn_gr_a_1(lonfx*gis_dyn%lotp,lats_dim_ext) )
!     allocate ( gis_dyn%pyn_gr_a_2(lonfx*gis_dyn%lotp,lats_dim_ext) )
!     allocate ( gis_dyn%dyn_gr_a_1(lonfx*gis_dyn%lotd,lats_dim_ext) )
!     allocate ( gis_dyn%dyn_gr_a_2(lonfx*gis_dyn%lotd,lats_dim_ext) )
!     allocate ( gis_dyn%anl_gr_a_1(lonfx*gis_dyn%lota,lats_dim_ext) )
!     allocate ( gis_dyn%anl_gr_a_2(lonfx*gis_dyn%lota,lats_dim_ext) )
! stochastic perturbations
      if (do_shum) allocate (gis_dyn%shum_wts(lonf,lats_node_a_max,levs) )
      if (do_sppt) allocate (gis_dyn%sppt_wts(lonf,lats_node_a_max,levs) )
      if (do_skeb) allocate (gis_dyn%skebu_wts(lonf,lats_node_a_max,levs) )
      if (do_skeb) allocate (gis_dyn%skebv_wts(lonf,lats_node_a_max,levs) )
      if (do_vc) allocate (gis_dyn%vcu_wts(lonf,lats_node_a_max,levs) )
      if (do_vc) allocate (gis_dyn%vcv_wts(lonf,lats_node_a_max,levs) )
!!
!** allocate digital filter vars

      gis_dyn%grid_gr_dfi%z_imp    = 0 ; gis_dyn%grid_gr_dfi%ps_imp     = 0
      gis_dyn%grid_gr_dfi%z_imp    = 0 ; gis_dyn%grid_gr_dfi%ps_imp     = 0
      gis_dyn%grid_gr_dfi%temp_imp = 0 ; gis_dyn%grid_gr_dfi%u_imp      = 0
      gis_dyn%grid_gr_dfi%v_imp    = 0 ; gis_dyn%grid_gr_dfi%tracer_imp = 0
      gis_dyn%grid_gr_dfi%p_imp    = 0 ; gis_dyn%grid_gr_dfi%dp_imp     = 0
      gis_dyn%grid_gr_dfi%dpdt_imp = 0


      if(gis_dyn%ndfi > 0 ) then
        if(gis_dyn%esmf_sta_list%z_import == 1) then
          gis_dyn%grid_gr_dfi%z_imp = 1
          allocate( gis_dyn%grid_gr_dfi%hs(lonf,lats_node_a_max,1))
!$omp parallel do private(i,j)
          do j=1,lats_node_a_max
            do i=1,lonf
              gis_dyn%grid_gr_dfi%hs(i,j,1) = 0.
            enddo
          enddo
        endif
        if(gis_dyn%esmf_sta_list%ps_import == 1) then
          gis_dyn%grid_gr_dfi%ps_imp = 1
          allocate( gis_dyn%grid_gr_dfi%ps(lonf,lats_node_a_max,1))
!$omp parallel do private(i,j)
          do j=1,lats_node_a_max
            do i=1,lonf
              gis_dyn%grid_gr_dfi%ps(i,j,1) = 0.
            enddo
          enddo
        endif
        if(gis_dyn%esmf_sta_list%temp_import == 1) then
          gis_dyn%grid_gr_dfi%temp_imp = 1
          allocate( gis_dyn%grid_gr_dfi%t(lonf,lats_node_a_max,levs))
!$omp parallel do private(i,j,k)
          do k=1,levs
            do j=1,lats_node_a_max
              do i=1,lonf
                gis_dyn%grid_gr_dfi%t(i,j,k) = 0.
              enddo
            enddo
          enddo
        endif
        if(gis_dyn%esmf_sta_list%u_import == 1) then
          gis_dyn%grid_gr_dfi%u_imp = 1
          allocate( gis_dyn%grid_gr_dfi%u(lonf,lats_node_a_max,levs))
!$omp parallel do private(i,j,k)
          do k=1,levs
            do j=1,lats_node_a_max
              do i=1,lonf
                gis_dyn%grid_gr_dfi%u(i,j,k) = 0.
              enddo
            enddo
          enddo
        endif
        if(gis_dyn%esmf_sta_list%v_import == 1) then
          gis_dyn%grid_gr_dfi%v_imp = 1
          allocate( gis_dyn%grid_gr_dfi%v(lonf,lats_node_a_max,levs))
!$omp parallel do private(i,j,k)
          do k=1,levs
            do j=1,lats_node_a_max
              do i=1,lonf
                gis_dyn%grid_gr_dfi%v(i,j,k) = 0.
              enddo
            enddo
          enddo
        endif
        if(gis_dyn%esmf_sta_list%tracer_import == 1) then
          gis_dyn%grid_gr_dfi%tracer_imp = 1
          allocate( gis_dyn%grid_gr_dfi%tracer(lonf,lats_node_a_max,levh))
!$omp parallel do private(i,j,k)
          do k=1,levh
            do j=1,lats_node_a_max
              do i=1,lonf
                gis_dyn%grid_gr_dfi%tracer(i,j,k) = 0.
              enddo
            enddo
          enddo
        endif
        if(gis_dyn%esmf_sta_list%p_import == 1) then
          gis_dyn%grid_gr_dfi%p_imp = 1
          allocate( gis_dyn%grid_gr_dfi%p(lonf,lats_node_a_max,levs))
!$omp parallel do private(i,j,k)
          do k=1,levs
            do j=1,lats_node_a_max
              do i=1,lonf
                gis_dyn%grid_gr_dfi%p(i,j,k) = 0.
              enddo
            enddo
          enddo
        endif
        if(gis_dyn%esmf_sta_list%dp_import == 1) then
          gis_dyn%grid_gr_dfi%dp_imp = 1
          allocate( gis_dyn%grid_gr_dfi%dp(lonf,lats_node_a_max,levs))
!!$omp parallel do private(i,j,k)
!         do k=1,levs
!           do j=1,lats_node_a_max
!             do i=1,lonf
!               gis_dyn%grid_gr_dfi%dp(i,j,k) = 0.
!             enddo
!           enddo
!         enddo
        endif
        if(gis_dyn%esmf_sta_list%dpdt_import == 1) then
          gis_dyn%grid_gr_dfi%dpdt_imp = 1
          allocate( gis_dyn%grid_gr_dfi%dpdt(lonf,lats_node_a_max,levs))
!!$omp parallel do private(i,j,k)
!         do k=1,levs
!           do j=1,lats_node_a_max
!             do i=1,lonf
!               gis_dyn%grid_gr_dfi%dpdt(i,j,k) = 0.
!             enddo
!           enddo
!         enddo
        endif
      endif
!
!## allocate output vars
      allocate(buff_mult_pieceg(lonf,lats_node_a_max,ngrids_gg))
!$omp parallel do private(i,j,k)
      do k=1,ngrids_gg
        do j=1,lats_node_a_max
          do i=1,lonf
            buff_mult_pieceg(i,j,k) = 0.
          enddo
        enddo
      enddo
!
      adiabatic = gis_dyn%adiabatic
!##

!!
      allocate ( gis_dyn%fhour_idate(1,5) )
!      write(0,*)'after allocate fhour_idate'
!
      if (me == 0) then
        print*, ' lats_dim_a=', lats_dim_a, ' lats_node_a=', lats_node_a
!       print*, ' lats_dim_ext=', lats_dim_ext,                           &
!               ' lats_node_ext=', lats_node_ext
      endif
!c
!$omp parallel do private(i,j,k)
      do k=1,lotgr
        do j=1,lats_node_a_max
          do i=1,lonf
            gis_dyn%grid_gr(i,j,k)  = 0.0
          enddo
        enddo
      enddo
!$omp parallel do private(i,j,k)
      do k=1,lotgr6*2
        do j=1,lats_node_a_max
          do i=1,lonf
            gis_dyn%grid_gr6(i,j,k) = 0.0
          enddo
        enddo
      enddo
!$omp parallel do private(i,j)
      do j=1,lats_node_a
        do i=1,lonf
          gis_dyn%ptot(i,j) = 0.0
          gis_dyn%pwat(i,j) = 0.0
        enddo
      enddo
      do k=1,ntrac
!$omp parallel do private(i,j)
        do j=1,lats_node_a
          do i=1,lonf
            gis_dyn%ptrc(i,j,k) = 0.0                            !glbsum
          enddo
        enddo
      enddo

!......................................................................
!
!      write(0,*) 'number of latitudes ext. :',lats_node_ext,              &
!                  lats_dim_ext,lats_node_a
!!
!!
!      print *,' sig_ini=',gis_dyn%nam_gfs_dyn%sig_ini,			  &
!              ' sig_ini2=',gis_dyn%nam_gfs_dyn%sig_ini2 

      gis_dyn%pdryini = 0.0

      if (me == 0) then
        print *,' grid_ini=',trim(gis_dyn%nam_gfs_dyn%grid_ini),'fhrot=',fhrot,    &
        'fhini=',fhini,'restart_run=',gis_dyn%restart_run
      endif

      cfile  = gis_dyn%nam_gfs_dyn%sig_ini
      cfile2 = gis_dyn%nam_gfs_dyn%sig_ini2
      n1     = 11
      n2     = 12
!
      if (me == 0) write(0,*)' before input_fields_slg lotls=',lotls
      CALL input_fields_slg(n1, n2, gis_dyn%pdryini,                                 &
                            gis_dyn%trie_ls, gis_dyn%trio_ls, gis_dyn%grid_gr ,      &
                            gis_dyn%ls_node, gis_dyn%ls_nodes, gis_dyn%max_ls_nodes, &
                            gis_dyn%snnp1ev, gis_dyn%snnp1od,                        &
                            gis_dyn%epse,    gis_dyn%epso,                           &
                            cfile, cfile2,   gis_dyn%restart_run,                    &
                            gis_dyn%global_lats_a, gis_dyn%lonsperlat,               &

!                           gis_dyn%global_lats_a, gis_dyn%lonsperlat,gis_dyn%epsedn,&
!                           gis_dyn%epsodn,gis_dyn%plnev_a, gis_dyn%plnod_a,         &

                            gis_dyn%plnew_a, gis_dyn%plnow_a, gis_dyn%lats_nodes_a)

!                           gis_dyn%pwat, gis_dyn%ptot, gis_dyn%ptrc)
!
      if(.not. gis_dyn%restart_run .or. fhrot == 0.0 ) then
          gis_dyn%start_step    = .true.
          gis_dyn%reset_step    = .false.
          gis_dyn%restart_step  = .false.
      ELSE
        if (me == 0) print *,'restart,filenames=', TRIM(cfile),', ', TRIM(cfile2)
        gis_dyn% start_step    = .false.
        gis_dyn% reset_step    = .false.
        gis_dyn% restart_step  = .true.
        gis_dyn% ndfi = gis_dyn% ndfi + gis_dyn% kdt
        print *,'in gfs dyn init, ndfi=',gis_dyn% ndfi,'kdt=',gis_dyn%kdt
      END IF
!!
!!
      gis_dyn%zhour = fhour
      if (iau) then
         call init_iau(gis_dyn)
      endif
      rc            = 0
      dyn_ini_time  = dyn_ini_time + (timef() - btime)

      write(0,*)' dyn_ini_time=',dyn_ini_time*1.0e-3,' me=',me
!
!
      end subroutine gfs_dynamics_initialize_slg
!
      end module gfs_dynamics_initialize_slg_mod
