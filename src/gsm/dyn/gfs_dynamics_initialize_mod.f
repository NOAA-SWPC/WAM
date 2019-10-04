
! !module: gfs_dynamics_initialize_mod 
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
!  Jan 20  2013 J. Wang          idea diffusion init interface change
!  Jun 26  2014 S. Moorthi       modified to read lonsperlat from a file
!  Jul 21  2014 S. Moorthi       removed num_reduce ; some cleanup
!  Aug 17  2016 P. Pegion        add call for iau update init
!  Nov 17  2017 W. Yang          modify the code for WAM-IPE coupling restart run.
!
! !interface:
!
      module gfs_dynamics_initialize_mod
!
!!uses:
!
      use gfs_dynamics_getcf_mod
      use gfs_dyn_machine, only : kind_io4
      use nemsio_module , only : nemsio_init
!
      use gfs_dyn_write_state, only : buff_mult_pieceg
      use gfs_dyn_layout1, only : ipt_lats_node_a, lats_node_a_max
      use gfs_dyn_resol_def, only : adiabatic, thermodyn_id, sfcpress_id
      use namelist_dynamics_def, only : fhrot,fhini,nemsio_in,do_shum,&
                                        do_sppt,do_skeb,do_vc,iau,    &
                                        wam_ipe_cpl_rst_input,        &
                                        wam_ipe_cpl_rst_output
      use gfs_dyn_tracer_config, only: gfs_dyn_tracer, tracer_config_init,gfs_dyn_tracer
      use gfs_dyn_io_header, only: z_r,z
!  stochastic perturbations
          use gfs_dyn_stoch_data, only: init_stochdata
!   iau forcing
       use gfs_dyn_iau_module, only: init_iau
!#ifndef IBM
!     USE omp_lib
!#endif

      implicit none

      contains

      subroutine gfs_dynamics_initialize(gis_dyn, rc)

! this subroutine set up the internal state variables,
! allocate internal state arrays for initializing the gfs system.
!----------------------------------------------------------------
!
      implicit none
!
      integer, parameter :: iunit=101
      type(gfs_dynamics_internal_state), pointer, intent(inout) :: gis_dyn
      integer,                                    intent(out)   :: rc

      logical           :: file_exists=.false.

      integer 		:: i, j, l, n, ilat, locl, ikey, nrank_all,     &
                           num_parthds, ierr, latghf, iret
!!
      character(20) cfile


! set up gfs internal state dimension and values for dynamics etc
!-------------------------------------------------------------------
      me     = gis_dyn%me
      if (me == 0)                                                      &
      write(0,*)'in initial,nbefore allocate lonsperlat,',allocated(gis_dyn%lonsperlat),'latg=',latg
!
      nodes  = gis_dyn%nodes
      nlunit = gis_dyn%nam_gfs_dyn%nlunit

      call compns_dynamics(gis_dyn%nam_gfs_dyn%deltim, gis_dyn%iret,    &
                           gis_dyn%ntrac,                               &
                           gis_dyn%nxpt,   gis_dyn%nypt, gis_dyn%jintmx,&
                           gis_dyn%jcap,                              	&
                           gis_dyn%levs,   gis_dyn%levr, 		&
                           gis_dyn%lonf,   gis_dyn%latg,          	&
                           gis_dyn%ntoz,   gis_dyn%ntcw, gis_dyn%ncld, 	&
                           gis_dyn%ntke,   gis_dyn%spectral_loop, me,  	&
                           gis_dyn%thermodyn_id,gis_dyn%sfcpress_id,    &
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

      if (gis_dyn%nam_gfs_dyn%total_member <= 1) then
        ens_nam=' '
      else
        write(ens_nam,'("_",i2.2)') gis_dyn%nam_gfs_dyn%member_id
      endif
!
      levh   = ntrac*levs
      latgd  = latg+ 2*jintmx 
      jcap1  = jcap+1 
      jcap2  = jcap+2 
      latg2  = latg/2 
      levm1  = levs-1 
      levp1  = levs+1 
      lonfx  = lonf + 1 + 2*nxpt+1 
      lnt    = jcap2*jcap1/2 
      lnuv   = jcap2*jcap1 
      lnt2   = 2*lnt 
      lnt22  = 2*lnt+1 
      lnte   = (jcap2/2)*((jcap2/2)+1)-1 
      lnto   = (jcap2/2)*((jcap2/2)+1)-(jcap2/2) 
      lnted  = lnte 
      lntod  = lnto 

!jw      ngrids_gg       = 2+levs*(4+ntrac)
      ngrids_gg       = 2+levs*(5+ntrac)
      gis_dyn%lnt2    = lnt2

      allocate(lat1s_a(0:jcap))
      allocate(lon_dims_a(latgd))
      allocate(lon_dims_ext(latgd))

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

      if (nemsio_in) then
        call nemsio_init(ierr)
      endif
!
      if (me == 0)                                                       &
      write(0,*)'before allocate lonsperlat,',allocated(gis_dyn%lonsperlat),'latg=',latg
      allocate(gis_dyn%lonsperlat(latg))

      inquire (file="lonsperlat.dat", exist=file_exists)
      if ( .not. file_exists ) then
        if ( me == 0 ) then
          print *,'   Requested lonsperlat.dat  data file does not exist'
          print *,'   *** Stopped in subroutine GFS_Init !!'
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
        p_zem  =         1     !     zeme/o(lnte/od,2,levs),
        p_dim  = p_zem  +levs  !     dime/o(lnte/od,2,levs),
        p_tem  = p_dim  +levs  !     teme/o(lnte/od,2,levs),
      if( .not. ndslfv ) then
        p_rm   = p_tem  +levs  !      rme/o(lnte/od,2,levh),
        p_dpm  = p_rm   +levh  !      qme/o(lnte/od,2),
      else
        p_dpm  = p_tem  +levs  !      rme/o(lnte/od,2,levh),
      endif
        p_qm   = p_dpm  +levs  !      qme/o(lnte/od,2),

        p_ze   = p_qm   +1     !      zee/o(lnte/od,2,levs),
        p_di   = p_ze   +levs  !      die/o(lnte/od,2,levs),
        p_te   = p_di   +levs  !      tee/o(lnte/od,2,levs),
      if( .not. ndslfv ) then
        p_rq   = p_te   +levs  !      rqe/o(lnte/od,2,levh),
        p_dp   = p_rq   +levh  !       qe/o(lnte/od,2),
      else
        p_dp   = p_te   +levs  !      rqe/o(lnte/od,2,levh),
      endif
        p_q    = p_dp   +levs  !       qe/o(lnte/od,2),
        p_dlam = p_q    +1     !  dpdlame/o(lnte/od,2),
        p_dphi = p_dlam +1     !  dpdphie/o(lnte/od,2),
        p_uln  = p_dphi +1     !     ulne/o(lnte/od,2,levs),
        p_vln  = p_uln  +levs  !     vlne/o(lnte/od,2,levs),
        p_zslam= p_vln  +levs  !    zslam/o(lnte/od,2),
        p_zsphi= p_zslam+1     !    zsphi/o(lnte/od,2),
                                                                                
        p_zz   = p_zsphi+1
        p_dpphi= p_zz   +levs
        p_zzphi= p_dpphi+levs
        p_dplam= p_zzphi+levs
        p_zzlam= p_dplam+levs

        p_w    = p_zzlam+levs  !       we/o(lnte/od,2,levs),
        p_x    = p_w    +levs  !       xe/o(lnte/od,2,levs),
        p_y    = p_x    +levs  !       ye/o(lnte/od,2,levs),
      if( .not. ndslfv ) then
        p_rt   = p_y    +levs  !      rte/o(lnte/od,2,levh),
        p_dpn  = p_rt   +levh  !      dpe/o(lnte/od,2,levs),
      else
        p_dpn  = p_y    +levs  !      rte/o(lnte/od,2,levh),
      endif
        p_zq   = p_dpn  +levs  !      zqe/o(lnte/od,2),
        p_gz   = p_zq   +1     !      gze/o(lnte/od,2),
        p_lapgz= p_gz   +1     !      gze/o(lnte/od,2),

      lotls  = p_lapgz

      g_gz   = 1
      g_uum  = g_gz  + 1        !  for grid point
      g_vvm  = g_uum + levs     !  for grid point
      g_ttm  = g_vvm + levs     !  for grid point
      g_rm   = g_ttm + levs     !  for grid point
      g_dpm  = g_rm  + levh     !  for grid point
      g_qm   = g_dpm + levs     !  for grid point

      g_uu   = g_qm  + 1        !  for grid point
      g_vv   = g_uu  + levs     !  for grid point
      g_tt   = g_vv  + levs     !  for grid point
      g_rq   = g_tt  + levs     !  for grid point
      g_dp   = g_rq  + levh     !  for grid point
      g_q    = g_dp  + levs     !  for grid point

      g_u    = g_q   + 1        !  for grid point
      g_v    = g_u   + levs     !  for grid point
      g_t    = g_v   + levs     !  for grid point
      g_rt   = g_t   + levs     !  for grid point
      g_dpn  = g_rt  + levh     !  for grid point
      g_zq   = g_dpn + levs     !  for grid point

      g_p    = g_zq  + 1   	!  for grid point 
      g_dpdt = g_p   + levs   	!  for grid point 
      g_zz   = g_dpdt+ levs     !  for grid point

      g_uup   = g_zz  + levs    !  for grid point
      g_vvp   = g_uup + levs    !  for grid point
      g_ttp   = g_vvp + levs    !  for grid point
      g_rqp   = g_ttp + levs    !  for grid point
      g_dpp   = g_rqp + levh    !  for grid point
      g_zqp   = g_dpp + levs    !  for grid point

      lotgr  = g_zqp
      lotgr6 = 4*levs+1*levh+1
!c
      if( .not. ndslfv ) then
        lots = 6*levs+1*levh+5
        lotd = 6*levs+2*levh+0
        lota = 5*levs+1*levh+2
      else
        lots = 6*levs+5
        lotd = 6*levs
        lota = 5*levs+2
      endif
      lotp = 4*levs
!
        ksz     = 1
        ksd     = ksz    + levs
        kst     = ksd    + levs
      if( .not. ndslfv ) then
        ksr     = kst    + levs
        ksdp    = ksr    + levh
      else
        ksdp    = kst    + levs
      endif
        ksq     = ksdp   + levs
        ksplam  = ksq    + 1
        kspphi  = ksplam + 1
        ksu     = kspphi + 1
        ksv     = ksu    + levs
        kzslam  = ksv    + levs
        kzsphi  = kzslam + 1
!
        kau   = 1
        kav   = kau + levs
        kat   = kav + levs
      if( .not. ndslfv ) then
        kar   = kat + levs
        kadp  = kar + levh
      else
        kadp  = kat + levs
      endif
        kaps  = kadp + levs
        kazs  = kaps + 1
        kap2  = kazs + 1
!
      kdpphi  = 1
      kzzphi  = kdpphi + levs
      kdplam  = kzzphi + levs
      kzzlam  = kdplam + levs
!
        kdtphi  = 1
      if( .not. ndslfv ) then
        kdrphi  = kdtphi + levs
        kdtlam  = kdrphi + levh
        kdrlam  = kdtlam + levs
        kdulam  = kdrlam + levh
      else
        kdtlam  = kdtphi + levs
        kdulam  = kdtlam + levs
      endif
        kdvlam  = kdulam + levs
        kduphi  = kdvlam + levs
        kdvphi  = kduphi + levs
!
! point to internal state
        gis_dyn%p_zem   = p_zem            !     zeme/o(lnte/od,2,levs),
        gis_dyn%p_dim   = p_dim            !     dime/o(lnte/od,2,levs),
        gis_dyn%p_tem   = p_tem            !     teme/o(lnte/od,2,levs),
      if( .not. ndslfv ) then
        gis_dyn%p_rm    = p_rm             !      rme/o(lnte/od,2,levh),
      endif
        gis_dyn%p_dpm   = p_dpm            !     dpme/o(lnte/od,2),
        gis_dyn%p_qm    = p_qm             !      qme/o(lnte/od,2),
        gis_dyn%p_zslam = p_zslam          ! hmhj
        gis_dyn%p_zsphi = p_zsphi          ! hmhj
        gis_dyn%p_ze    = p_ze             !      zee/o(lnte/od,2,levs),
        gis_dyn%p_di    = p_di             !      die/o(lnte/od,2,levs),
        gis_dyn%p_te    = p_te             !      tee/o(lnte/od,2,levs),
      if( .not. ndslfv ) then
        gis_dyn%p_rq    = p_rq             !      rqe/o(lnte/od,2,levh),
      endif
        gis_dyn%p_dp    = p_dp             !      dpe/o(lnte/od,2),
        gis_dyn%p_q     = p_q              !       qe/o(lnte/od,2),
        gis_dyn%p_dlam  = p_dlam           !  dpdlame/o(lnte/od,2),
        gis_dyn%p_dphi  = p_dphi           !  dpdphie/o(lnte/od,2),
        gis_dyn%p_uln   = p_uln            !     ulne/o(lnte/od,2,levs),
        gis_dyn%p_vln   = p_vln            !     vlne/o(lnte/od,2,levs),
        gis_dyn%p_w     = p_w              !       we/o(lnte/od,2,levs),
        gis_dyn%p_x     = p_x              !       xe/o(lnte/od,2,levs),
        gis_dyn%p_y     = p_y              !       ye/o(lnte/od,2,levs),
      if( .not. ndslfv ) then
        gis_dyn%p_rt    = p_rt             !      rte/o(lnte/od,2,levh),
      endif
        gis_dyn%p_dpn   = p_dpn            !     dpne/o(lnte/od,2)
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
!
      gis_dyn%lotls = lotls
      gis_dyn%lotgr = lotgr
      gis_dyn%lots = lots 
      gis_dyn%lotd = lotd
      gis_dyn%lota = lota
      gis_dyn%lotp = lotp
!
      allocate(gis_dyn%tee1(levs))

!     print *,' finish dimension in gfs_dynamics_initialize '

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!!
!!      create io communicator and comp communicator
!!
      if (me == 0) write(0,*) 'io option ,liope :',liope
!!
      nodes_comp = nodes
!
      if (me == 0) then
        print 100, jcap,levs
100   format (' smf ',i3,i3,' created august 2000 ev od ri ')
!#ifdef IBM
        print*,'number of threads is',num_parthds()
!#else
!       print*,'number of threads is',omp_get_num_threads()
!#endif
        print*,'number of mpi procs is',nodes
      endif
!
      gis_dyn%cons0    =    0.0d0     !constant
      gis_dyn%cons0p5  =    0.5d0     !constant
      gis_dyn%cons1200 = 1200.d0      !constant
      gis_dyn%cons3600 = 3600.d0      !constant
!
      ls_dim = (jcap1-1)/nodes+1
!!
!cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
      if (me == 0)                                                       &
       write(0,*)'before allocate ls_nodes,',allocated(gis_dyn%ls_nodes),&
       'ls_dim=', ls_dim,'nodes=',nodes
!
      allocate ( gis_dyn%ls_node (ls_dim*3) )
      allocate ( gis_dyn%ls_nodes(ls_dim,nodes) )
      allocate ( gis_dyn%max_ls_nodes(nodes) )
      gis_dyn%ls_node      = 0
      gis_dyn%ls_nodes     = 0
      gis_dyn%max_ls_nodes = 0

!
      allocate ( gis_dyn%lats_nodes_a_fix(nodes))     ! added for mGrid
!
      allocate ( gis_dyn%lats_nodes_a(nodes) )
      allocate ( gis_dyn%global_lats_a(latg) )
!
      allocate ( gis_dyn%lats_nodes_ext(nodes) )
      allocate ( gis_dyn%global_lats_ext(latg+2*jintmx+2*nypt*(nodes-1)) )
!
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
                            gis_dyn%ls_node(ls_dim+1:2*ls_dim),'3dim=',   &
                            gis_dyn%ls_node(2*ls_dim+1:3*ls_dim)
!c
!c
      allocate ( gis_dyn%epse  (len_trie_ls) )
      allocate ( gis_dyn%epso  (len_trio_ls) )
      allocate ( gis_dyn%epsedn(len_trie_ls) )
      allocate ( gis_dyn%epsodn(len_trio_ls) )
!c
      allocate ( gis_dyn%snnp1ev(len_trie_ls) )
      allocate ( gis_dyn%snnp1od(len_trio_ls) )
!c
      allocate ( gis_dyn%ndexev(len_trie_ls) )
      allocate ( gis_dyn%ndexod(len_trio_ls) )
!c
      allocate ( gis_dyn%plnev_a(len_trie_ls,latg2) )
      allocate ( gis_dyn%plnod_a(len_trio_ls,latg2) )
      allocate ( gis_dyn%pddev_a(len_trie_ls,latg2) )
      allocate ( gis_dyn%pddod_a(len_trio_ls,latg2) )
      allocate ( gis_dyn%plnew_a(len_trie_ls,latg2) )
      allocate ( gis_dyn%plnow_a(len_trio_ls,latg2) )
!c
      gis_dyn%maxstp = 36

 
      if(me == 0) 							&
        print*,'from compns_dynamics : iret=',gis_dyn%iret		&
       ,' nsout=',nsout,' nsres=',nsres

      if(gis_dyn%iret /= 0) then
        if(me == 0) print *,' incompatible namelist - aborted in main'
        call mpi_quit(13)
      endif
!!
      gis_dyn%lats_nodes_ext = 0
      call getcon_dynamics(gis_dyn%n3,gis_dyn%n4,			&
           gis_dyn%ls_node,gis_dyn%ls_nodes,gis_dyn%max_ls_nodes,       &
           gis_dyn%lats_nodes_a,gis_dyn%global_lats_a,                  &
           gis_dyn%lonsperlat,gis_dyn%lats_node_a_max,                  &
           gis_dyn%lats_nodes_ext,gis_dyn%global_lats_ext,              &
           gis_dyn%epse,gis_dyn%epso,gis_dyn%epsedn,gis_dyn%epsodn,     &
           gis_dyn%snnp1ev,gis_dyn%snnp1od,				&
           gis_dyn%ndexev,gis_dyn%ndexod,  				&
           gis_dyn%plnev_a,gis_dyn%plnod_a,				&
           gis_dyn%pddev_a,gis_dyn%pddod_a, 				&
           gis_dyn%plnew_a,gis_dyn%plnow_a,gis_dyn%colat1)
!
      gis_dyn%lats_node_a     = gis_dyn%lats_nodes_a(me+1)
      gis_dyn%ipt_lats_node_a = ipt_lats_node_a

      if (me == 0)                                                       &
      write(0,*)'after getcon_dynamics,lats_node_a=',gis_dyn%lats_node_a &
       ,'ipt_lats_node_a=',gis_dyn%ipt_lats_node_a,'ngptc=',ngptc
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

      call gfs_dyn_lonlat_para(gis_dyn%global_lats_a, gis_dyn%xlon,     &
                               gis_dyn%xlat, gis_dyn%lonsperlat)
!!
      allocate ( gis_dyn%trie_ls (len_trie_ls,2,lotls) )
      allocate ( gis_dyn%trio_ls (len_trio_ls,2,lotls) )
      allocate ( gis_dyn%grid_gr (lonf,lats_node_a_max,lotgr    ) )
      allocate ( gis_dyn%grid_gr6(lonf,lats_node_a_max,lotgr6*2 ) )
      allocate ( gis_dyn%pwat    (lonf,lats_node_a) )
      allocate ( gis_dyn%ptot    (lonf,lats_node_a) )
      allocate ( gis_dyn%ptrc    (lonf,lats_node_a,ntrac) )         !glbsum
      allocate ( z(lnt2) )
      allocate ( z_r(lnt2) )
!c
      allocate (   gis_dyn%syn_ls_a(4*ls_dim,gis_dyn%lots,latg2) )
      allocate (   gis_dyn%dyn_ls_a(4*ls_dim,gis_dyn%lotd,latg2) )
!c
      allocate (   gis_dyn%syn_gr_a_1(lonfx*gis_dyn%lots,lats_dim_ext) )
      allocate (   gis_dyn%syn_gr_a_2(lonfx*gis_dyn%lots,lats_dim_ext) )
      allocate (   gis_dyn%pyn_gr_a_1(lonfx*gis_dyn%lotp,lats_dim_ext) )
      allocate (   gis_dyn%pyn_gr_a_2(lonfx*gis_dyn%lotp,lats_dim_ext) )
      allocate (   gis_dyn%dyn_gr_a_1(lonfx*gis_dyn%lotd,lats_dim_ext) )
      allocate (   gis_dyn%dyn_gr_a_2(lonfx*gis_dyn%lotd,lats_dim_ext) )
      allocate (   gis_dyn%anl_gr_a_1(lonfx*gis_dyn%lota,lats_dim_ext) )
      allocate (   gis_dyn%anl_gr_a_2(lonfx*gis_dyn%lota,lats_dim_ext) )
! stochastic perturbations
      if (do_shum) allocate (gis_dyn%shum_wts(lonf,lats_node_a_max,levs) )
      if (do_sppt) allocate (gis_dyn%sppt_wts(lonf,lats_node_a_max,levs) )
      if (do_skeb) allocate (gis_dyn%skebu_wts(lonf,lats_node_a_max,levs) )
      if (do_skeb) allocate (gis_dyn%skebv_wts(lonf,lats_node_a_max,levs) )
      if (do_vc) allocate (gis_dyn%vcu_wts(lonf,lats_node_a_max,levs) )
      if (do_vc) allocate (gis_dyn%vcv_wts(lonf,lats_node_a_max,levs) )
!
!** allocate digital filter vars
      gis_dyn%grid_gr_dfi%z_imp    = 0 ; gis_dyn%grid_gr_dfi%ps_imp     = 0
      gis_dyn%grid_gr_dfi%z_imp    = 0 ; gis_dyn%grid_gr_dfi%ps_imp     = 0
      gis_dyn%grid_gr_dfi%temp_imp = 0 ; gis_dyn%grid_gr_dfi%u_imp      = 0
      gis_dyn%grid_gr_dfi%v_imp    = 0 ; gis_dyn%grid_gr_dfi%tracer_imp = 0
      gis_dyn%grid_gr_dfi%p_imp    = 0 ; gis_dyn%grid_gr_dfi%dp_imp     = 0
      gis_dyn%grid_gr_dfi%dpdt_imp = 0
      if(gis_dyn%ndfi>0.and. gis_dyn%esmf_sta_list%z_import==1) then
        gis_dyn%grid_gr_dfi%z_imp=1
        allocate( gis_dyn%grid_gr_dfi%hs(lonf,lats_node_a_max,1))
      endif
      if(gis_dyn%ndfi>0.and. gis_dyn%esmf_sta_list%ps_import==1) then
        gis_dyn%grid_gr_dfi%ps_imp=1
        allocate( gis_dyn%grid_gr_dfi%ps(lonf,lats_node_a_max,1))
      endif
      if(gis_dyn%ndfi>0.and. gis_dyn%esmf_sta_list%temp_import==1) then
        gis_dyn%grid_gr_dfi%temp_imp=1
        allocate( gis_dyn%grid_gr_dfi%t(lonf,lats_node_a_max,levs))
      endif
      if(gis_dyn%ndfi>0.and. gis_dyn%esmf_sta_list%u_import==1) then
        gis_dyn%grid_gr_dfi%u_imp=1
        allocate( gis_dyn%grid_gr_dfi%u(lonf,lats_node_a_max,levs))
      endif
      if(gis_dyn%ndfi>0.and. gis_dyn%esmf_sta_list%v_import==1) then
        gis_dyn%grid_gr_dfi%v_imp=1
        allocate( gis_dyn%grid_gr_dfi%v(lonf,lats_node_a_max,levs))
      endif
      if(gis_dyn%ndfi>0.and. gis_dyn%esmf_sta_list%tracer_import==1) then
        gis_dyn%grid_gr_dfi%tracer_imp=1
        allocate( gis_dyn%grid_gr_dfi%tracer(lonf,lats_node_a_max,ntrac*levs))
      endif
      if(gis_dyn%ndfi>0.and. gis_dyn%esmf_sta_list%p_import==1) then
        gis_dyn%grid_gr_dfi%p_imp=1
        allocate( gis_dyn%grid_gr_dfi%p(lonf,lats_node_a_max,levs))
      endif
      if(gis_dyn%ndfi>0.and. gis_dyn%esmf_sta_list%dp_import==1) then
        gis_dyn%grid_gr_dfi%dp_imp=1
        allocate( gis_dyn%grid_gr_dfi%dp(lonf,lats_node_a_max,levs))
      endif
      if(gis_dyn%ndfi>0.and. gis_dyn%esmf_sta_list%dpdt_import==1) then
        gis_dyn%grid_gr_dfi%dpdt_imp=1
        allocate( gis_dyn%grid_gr_dfi%dpdt(lonf,lats_node_a_max,levs))
      endif
!
!## allocate output vars
      allocate(buff_mult_pieceg(lonf,lats_node_a_max,ngrids_gg))
      buff_mult_pieceg = 0.
      adiabatic = gis_dyn%adiabatic
!##

!!
      allocate (      gis_dyn%fhour_idate(1,5) )
!      write(0,*)'after allocate fhour_idate'
!
      if (me == 0) then
        print*, ' lats_dim_a=', lats_dim_a, ' lats_node_a=', lats_node_a
        print*, ' lats_dim_ext=', lats_dim_ext,                           &
                ' lats_node_ext=', lats_node_ext
      endif
!     write(150+me,*)' lats_dim_a=', lats_dim_a, ' lats_node_a=',lats_node_a
!c
      gis_dyn%grid_gr  = 0.0
      gis_dyn%grid_gr6 = 0.0
      gis_dyn%ptot     = 0.0
      gis_dyn%pwat     = 0.0
      gis_dyn%ptrc     = 0.0                            !glbsum

        ilat=lats_node_a
!c......................................................................
!c
!      write(0,*) 'number of latitudes ext. :',lats_node_ext,              &
!                  lats_dim_ext,lats_node_a
!!
!      print *,' sig_ini=',gis_dyn%nam_gfs_dyn%sig_ini,			  &
!              ' sig_ini2=',gis_dyn%nam_gfs_dyn%sig_ini2 
      gis_dyn%pdryini = 0.0

      if(lsidea) then
        if(wam_ipe_cpl_rst_input ) then
! Add the open lines for the WAM-IPE coupling restart file. WY.
!--------------------------------------------------------------
! For reading in restart file at the begining of the restart run.
!----------------------------------------------------------------
          open(unit=180, file='WAM_IPE_RST_rd', form='unformatted')
          rewind 180
        end if

        if(wam_ipe_cpl_rst_output ) then
! For writing out the restart file for the next restart WAM-IPE run.
!-------------------------------------------------------------------
          open(unit=181, file='WAM_IPE_RST_wrt', form='unformatted')
          rewind 181
        end if
      end if

      if (me == 0) then
        print *,' grid_ini=',trim(gis_dyn%nam_gfs_dyn%grid_ini),'fhrot=',fhrot, &
        'fhini=',fhini,'restart_run=',gis_dyn%restart_run
      endif

      if( .not. gis_dyn%restart_run) then
        if(nemsio_in) then
          cfile = gis_dyn%nam_gfs_dyn%grid_ini
        else
          cfile = gis_dyn%nam_gfs_dyn%sig_ini
        endif
!
        call input_fields(cfile, gis_dyn%pdryini,                            &
          gis_dyn%trie_ls, gis_dyn%trio_ls,  gis_dyn%grid_gr ,               &
          gis_dyn%ls_node, gis_dyn%ls_nodes, gis_dyn%max_ls_nodes,           &
          gis_dyn%global_lats_a, gis_dyn%lonsperlat,                         &
          gis_dyn%epse, gis_dyn%epso, gis_dyn%epsedn, gis_dyn%epsodn,        &
          gis_dyn%plnev_a, gis_dyn%plnod_a,                                  &
          gis_dyn%plnew_a, gis_dyn%plnow_a, gis_dyn%lats_nodes_a,            &
          gis_dyn%pwat, gis_dyn%ptot, gis_dyn%ptrc, gis_dyn%snnp1ev,         &
          gis_dyn%snnp1od)
!
          gis_dyn% start_step   = .true.
          gis_dyn% reset_step   = .false.
          gis_dyn% restart_step = .false.
      else
        if (me == 0) then
          print *,'restart,filename=',trim(gis_dyn%nam_gfs_dyn%grid_ini),       &
          trim(gis_dyn%nam_gfs_dyn%grid_ini2),trim(gis_dyn%nam_gfs_dyn%sig_ini),&
          trim(gis_dyn%nam_gfs_dyn%sig_ini2)
        endif
        call input_fields_rst(gis_dyn%nam_gfs_dyn%grid_ini,                  &
          gis_dyn%nam_gfs_dyn%grid_ini2,gis_dyn%nam_gfs_dyn%sig_ini,         &
          gis_dyn%nam_gfs_dyn%sig_ini2, gis_dyn%pdryini,                     &
          gis_dyn%trie_ls, gis_dyn%trio_ls,  gis_dyn%grid_gr ,               &
          gis_dyn%ls_node, gis_dyn%ls_nodes, gis_dyn%max_ls_nodes,           &
          gis_dyn%global_lats_a,  gis_dyn%lonsperlat,       &
          gis_dyn%epse, gis_dyn%epso, gis_dyn%plnev_a, gis_dyn%plnod_a,      &
          gis_dyn%plnew_a, gis_dyn%plnow_a, gis_dyn%snnp1ev,gis_dyn%snnp1od, &
          gis_dyn%lats_nodes_a)
!
          gis_dyn% start_step   = .false.
          gis_dyn% reset_step   = .false.
          gis_dyn% restart_step = .true.

! read in the restart file for WAM-IPE coupling restart run.
!-----------------------------------------------------------
          if(lsidea .and. wam_ipe_cpl_rst_input) then
            call input_for_wam_ipe_rst(gis_dyn%global_lats_a,   &
                                       gis_dyn%lonsperlat,      &
                                       gis_dyn%lats_nodes_a)
          end if
      endif
!!
      tov = 0.0
      if (.not. (hybrid.or.gen_coord_hybrid) ) then                   ! hmhj
       call setsig(si,ci,del,sl,cl,rdel2,tov,me)
       am=-8888888.
       bm=-7777777.
       call amhmtm(del,sv,am)
       call bmdi_sig(ci,bm)
      endif
      if( ndslfv ) then
        if(lsidea ) then
          call idea_deldifs_init                                        &
                (SL,gis_dyn%LS_NODE,hybrid,gen_coord_hybrid)
        else
          call deldifs_noq                                              &
                  (gis_dyn%epse,gis_dyn%epse,gis_dyn%epse,              &
                   gis_dyn%epse,gis_dyn%epse,gis_dyn%epse,              &
                   gis_dyn%epso,gis_dyn%epso,                           &
                   gis_dyn%epso,gis_dyn%epso,                           &
                   gis_dyn%cons0,sl,gis_dyn%ls_node,gis_dyn%epse,       &
                   0,hybrid,gen_coord_hybrid)
        endif
      else
        if(lsidea ) then
          CALL idea_deldifs_init                                        &
                (SL,gis_dyn%LS_NODE,hybrid,gen_coord_hybrid)
        else
          call deldifs(gis_dyn%epse,gis_dyn%epse,gis_dyn%epse,		&
                   gis_dyn%epse,gis_dyn%epse,gis_dyn%epse,      	&
                   gis_dyn%epso,gis_dyn%epso,gis_dyn%epso,		&
                   gis_dyn%epso,gis_dyn%epso,gis_dyn%epso,      	& 
                   gis_dyn%cons0,sl,gis_dyn%ls_node,gis_dyn%epse,	&
                   0,hybrid,gen_coord_hybrid)  
        endif
      endif
!!
      if (iau) then
          call init_iau(gis_dyn)
      endif
      gis_dyn%zhour = fhour
      rc = 0
!
! initialize stochastic module
           IF (do_shum .OR. do_sppt .OR. do_skeb .OR. do_vc) THEN
              call init_stochdata(gis_dyn%deltim,gis_dyn%LS_NODE)
           ENDIF
      end subroutine gfs_dynamics_initialize
!
      end module gfs_dynamics_initialize_mod
