!#include "../../../ESMFVersionDefine.h"

!
! !module: gfs_dynamics_run_mod --- run module of the grided
!                              component of the gfs dynamics system.
!
! !description: gfs run module.
!
! !revision history:
!
!  november 2004      weiyu yang initial code.
!  may      2005      weiyu yang, updated to the new version of gfs.
!  january  2007      hann-ming henry juang for gfs dynamics only
!  March    2009      Weiyu Yang, modified for the ensmeble NEMS run.
!  oct 05   2009      sarah lu, grid_gr unfolded from 2D to 3D
!  october  2009      jun wang, output data every time step
!  november 2009      jun wang, set digital filter variables
!  febrary  2010      jun wang, add restart step
!  aug  2010          sarah lu, compute tracer global sum
!  febrary  2011      henry juang add non-ierating dimensional-splitting semi_lagrangian
!                     (NDSL) advection with options of MASS_DP and NDSLFV
!  Feb      2011      Weiyu Yang, Updated to use both the ESMF 4.0.0rp2 library,
!                     ESMF 5 library and the the ESMF 3.1.0rp2 library.
!  Sep      2012      Jun Wang, add sigio option
!  dec      2014      Weiyu Yang, Add NEMSIO output.
!  feb     2016       S Moorthi - grid-point digital filter fix to filter at initial
!                                 time and during restart
!  aug      2016      Phil Pegion - add logic for iau

!
!
! !interface:
!
      module gfs_dynamics_run_mod
!
!!uses:
!
      USE ESMF
      use gfs_dynamics_internal_state_mod
      use namelist_dynamics_def, only : nemsio_in, semilag,do_shum,do_sppt, &
                                        do_skeb,do_vc,fhmax_hf,nsout_hf,nsskeb,skebint
      use gfs_get_pattern_mod
      use gfs_dyn_resol_def,     only : kdt_start
      use do_dynamics_mod,       only : do_dynamics_gridn2p
      use gfs_dyn_layout1,       only : ipt_lats_node_a,me,lats_node_a

      implicit none

      contains

      subroutine gfs_dynamics_run(gis_dyn, imp_gfs_dyn, rc)

      type(gfs_dynamics_internal_state), pointer, intent(inout) :: gis_dyn
      TYPE(ESMF_State),                           INTENT(inout) :: imp_gfs_dyn
      integer, optional,                          intent(out)   :: rc

!     real , save	:: timestep=0.0
      integer, save     :: kdt_save=0
      integer	           rc1, k, i, kdt_loc
      logical              lcom2mdl

      LOGICAL            :: Cpl_flag1
!***********************************************************************
!
!     lsfwd      logical true during a first forward step
!     phour      real forecast hour at the end of the previous timestep
!     fhrot      if =0 read one time level sigmas and te=tem zem=ze ...
!
!     lsout controls freq. of output
!     nsout controls freq. of output in time steps
!     fhout controls freq. of output in hours
!
!     nsres time steps between writing restart files
!
!***********************************************************************
       rc1 = 0
! ---------------------------------------------------------------------
! check if need do digital filter data collect
!
      kdt_loc = gis_dyn%kdt - kdt_start + 1
      gis_dyn%ldfi = gis_dyn%ndfi > 0 .and. kdt_loc  <= gis_dyn%ndfi   &
                                      .and. kdt_save <= gis_dyn%ndfi
! ---------------------------------------------------------------------
! change temperature and pressure back to model grid value at n+1  slot.
!
      if (.not. gis_dyn%start_step .and. .not. gis_dyn%reset_step       &
          .and. .not. gis_dyn%restart_step ) then
!
        if (fhmax_hf > 0 .and. nsout_hf > 0 .and. gis_dyn%phour <= fhmax_hf &
            .and. gis_dyn%kdt*gis_dyn%deltim <= fhmax_hf*3600.) then
          gis_dyn%lsout = mod(gis_dyn%kdt,nsout_hf) == 0 .or. gis_dyn%kdt == 1
        else
          gis_dyn%lsout = mod(gis_dyn%kdt,nsout) == 0 .or. gis_dyn%kdt == 1
        endif


        lcom2mdl = .not. (gis_dyn%ens .and. gis_dyn%cpl_flag) .and.         &
                   .not. (gis_dyn%kdt==0.and.(.not.nemsio_in) )


        IF(lcom2mdl) then
            call common_to_model_vars (gis_dyn%grid_gr(1,1,gis_dyn%g_zqp),  &
                                       gis_dyn%grid_gr(1,1,gis_dyn%g_ttp),  &
                                       gis_dyn%grid_gr(1,1,gis_dyn%g_rqp),  &
                                       gis_dyn%grid_gr(1,1,gis_dyn%g_uup),  &
                                       gis_dyn%grid_gr(1,1,gis_dyn%g_vvp),  &
                                       gis_dyn%grid_gr(1,1,gis_dyn%g_dpp),  &
                                       gis_dyn%global_lats_a,               &
                                       gis_dyn%lonsperlat)
        endif

!       lsfwd=.false.
!     else
!       print *,' the first time step, model running not from internal state.'
!       lsfwd=.true.
      endif
!
! ---------------------------------------------------------------------
!! get pwat and ptot
      if (.not. gis_dyn%reset_step ) then
        call getpwatptot(gis_dyn%grid_gr(1,1,gis_dyn%g_zqp),      &
                         gis_dyn%grid_gr(1,1,gis_dyn%g_ttp),      &
                         gis_dyn%grid_gr(1,1,gis_dyn%g_rqp),      &
                         gis_dyn%global_lats_a,                   &
                         gis_dyn%lonsperlat,                      &
                         gis_dyn%pwat, gis_dyn%ptot, gis_dyn%ptrc )  !glbsum
      endif

!
! ---------------------------------------------------------------------
! check whether reset step due to dfi
!
      kdt_save = kdt_save + 1
! get stochastic pertrubation patterns 
       IF ( do_shum) THEN
          call get_pattern_shum(gis_dyn%LS_NODE,gis_dyn%LS_NODES,gis_dyn%MAX_LS_NODES,  &
           gis_dyn%LATS_NODES_A,gis_dyn%GLOBAL_LATS_A,gis_dyn%LONSPERLAT,               &
           gis_dyn%PLNEV_A,gis_dyn%PLNOD_A,gis_dyn%shum_wts,gis_dyn%deltim)
       ENDIF
       IF ( do_sppt) THEN
          call get_pattern_sppt(gis_dyn%LS_NODE,gis_dyn%LS_NODES,gis_dyn%MAX_LS_NODES,  &
           gis_dyn%LATS_NODES_A,gis_dyn%GLOBAL_LATS_A,gis_dyn%LONSPERLAT,               &
           gis_dyn%PLNEV_A,gis_dyn%PLNOD_A,gis_dyn%sppt_wts,gis_dyn%deltim)
       ENDIF
       IF ( do_skeb) THEN
         if (me .EQ. 0 .OR. me .EQ. 8) then
             print*,'calling skeb',gis_dyn%trie_ls(1,1,p_ze)
          endif
         IF (mod(gis_dyn%kdt,nsskeb) == 1 .OR. nsskeb == 1) THEN
         call get_pattern_skeb(gis_dyn%trie_ls(1,1,gis_dyn%p_ze), gis_dyn%trio_ls(1,1,gis_dyn%p_ze),&
               gis_dyn%trie_ls(1,1,gis_dyn%p_di),  gis_dyn%trio_ls(1,1,gis_dyn%p_di),               &
               gis_dyn%LS_NODE,gis_dyn%LS_NODES,gis_dyn%MAX_LS_NODES,                   &
               gis_dyn%LATS_NODES_A,gis_dyn%GLOBAL_LATS_A,gis_dyn%LONSPERLAT,           &
               gis_dyn%epsedn,gis_dyn%epsodn,gis_dyn%snnp1ev,gis_dyn%snnp1od,           &
               gis_dyn%PLNEV_A,gis_dyn%PLNOD_A,gis_dyn%PLNEW_A,gis_dyn%PLNOW_A,         &
               gis_dyn%skebu_wts,gis_dyn%skebv_wts,skebint)
         ENDIF
       ENDIF
       IF ( do_vc) THEN
         call get_pattern_vc(gis_dyn%trie_ls(1,1,gis_dyn%p_ze), gis_dyn%trio_ls(1,1,gis_dyn%p_ze),&
               gis_dyn%trie_ls(1,1,gis_dyn%p_di), gis_dyn%trio_ls(1,1,gis_dyn%p_di),              &
               gis_dyn%LS_NODE,gis_dyn%LS_NODES,gis_dyn%MAX_LS_NODES,                   &
               gis_dyn%LATS_NODES_A,gis_dyn%GLOBAL_LATS_A,gis_dyn%LONSPERLAT,           &
               gis_dyn%epsedn,gis_dyn%epsodn,gis_dyn%snnp1ev,gis_dyn%snnp1od,           &
               gis_dyn%PLNEV_A,gis_dyn%PLNOD_A,gis_dyn%PLNEW_A,gis_dyn%PLNOW_A,         &
               gis_dyn%vcu_wts, gis_dyn%vcv_wts,gis_dyn%deltim)
       ENDIF


! ======================================================================
!                     do one time step with one-loop
! ---------------------------------------------------------------------
      IF(semilag) THEN
        CALL do_dynamics_slg_loop( gis_dyn%deltim,        gis_dyn%kdt,          &
                                   gis_dyn%PHOUR,         gis_dyn%TRIE_LS,      &
                                   gis_dyn%TRIO_LS,       gis_dyn%GRID_GR,      &
!                                  gis_dyn%grid_gr_dfi,   gis_dyn%GRID_GR6,     &
                                   gis_dyn%LS_NODE,       gis_dyn%LS_NODES,     &
                                   gis_dyn%MAX_LS_NODES,  gis_dyn%LATS_NODES_A, &
                                   gis_dyn%GLOBAL_LATS_A, gis_dyn%LONSPERLAT,   &
                                   gis_dyn%ldfi,          gis_dyn%grid_gr_dfi,  &
                                   gis_dyn%EPSE,          gis_dyn%EPSO,         &
                                   gis_dyn%EPSEDN,        gis_dyn%EPSODN,       &
                                   gis_dyn%SNNP1EV,       gis_dyn%SNNP1OD,      &
                                   gis_dyn%NDEXEV,        gis_dyn%NDEXOD,       &
                                   gis_dyn%PLNEV_A,       gis_dyn%PLNOD_A,      &
                                   gis_dyn%PDDEV_A,       gis_dyn%PDDOD_A,      &
                                   gis_dyn%PLNEW_A,       gis_dyn%PLNOW_A,      &
                                   gis_dyn%LSOUT,                               &
                                   gis_dyn%CFHOUR1,       gis_dyn%ENS,          &
                                   gis_dyn%start_step,    gis_dyn%restart_step, &
                                   gis_dyn%reset_step,    gis_dyn%end_step,     &
                                   gis_dyn%dfiend_step,   gis_dyn%pdryini,      &
                                   gis_dyn%nblck,         gis_dyn%zhour,        &
                                   gis_dyn%pwat,          gis_dyn%ptot,         &
                                   gis_dyn%ptrc,          gis_dyn%nfcstdate7,   &
                                   gis_dyn%iniauinterval)

        gis_dyn%phour = fhour

      ELSE ! for non-Semi-Lagrangian run.

        if( gis_dyn% spectral_loop == 1 ) then
!
!! grid_gr unfolded from 2D to 3D (Sarah Lu)
!      if(gis_dyn% kdt<=1) then
!          print *,'in gfs dyn run,bf one kdt=',gis_dyn%kdt,'ps=',maxval(gis_dyn%grid_gr(:,:,g_q)),   &
!        minval(gis_dyn%grid_gr(:,:,g_q)),'u=',maxval(gis_dyn%grid_gr(:,:,g_uu)),   &
!        minval(gis_dyn%grid_gr(:,:,g_uu)),'v=',maxval(gis_dyn%grid_gr(:,:,g_vv)),   &
!        minval(gis_dyn%grid_gr(:,:,g_vv)),'t=',maxval(gis_dyn%grid_gr(:,:,g_tt)),   &
!        minval(gis_dyn%grid_gr(:,:,g_tt)),'rq=',maxval(gis_dyn%grid_gr(:,:,g_rq)),   &
!        minval(gis_dyn%grid_gr(:,:,g_rq)),'psg=',maxval(gis_dyn%grid_gr(:,:,g_zq)),   &
!        minval(gis_dyn%grid_gr(:,:,g_zq)),'psm=',minval(gis_dyn%grid_gr(:,:,g_qm)),   &
!        minval(gis_dyn%grid_gr(:,:,g_qm))
!      endif

          call  do_dynamics_one_loop(    gis_dyn% deltim         ,	&
                   gis_dyn% kdt           ,gis_dyn% phour          ,	&
                   gis_dyn% trie_ls       ,gis_dyn% trio_ls        ,	&
                   gis_dyn% grid_gr       ,gis_dyn% grid_gr6       ,    &
                   gis_dyn% grid_gr_dfi   ,                             &
                   gis_dyn% ls_node       ,gis_dyn% ls_nodes       ,	&
                   gis_dyn% max_ls_nodes  ,				&
                   gis_dyn% lats_nodes_a  ,gis_dyn% global_lats_a  ,	&
                   gis_dyn% lonsperlat    ,				&
                   gis_dyn% lats_nodes_ext,gis_dyn% global_lats_ext,	&
                   gis_dyn% epse          ,gis_dyn% epso           ,	&
                   gis_dyn% epsedn        ,gis_dyn% epsodn         ,	&
                   gis_dyn% snnp1ev       ,gis_dyn% snnp1od        ,	&
                   gis_dyn% ndexev        ,gis_dyn% ndexod         ,	&
                   gis_dyn% plnev_a       ,gis_dyn% plnod_a        ,	&
                   gis_dyn% pddev_a       ,gis_dyn% pddod_a        ,	&
                   gis_dyn% plnew_a       ,gis_dyn% plnow_a        ,	&
                   gis_dyn% syn_ls_a      ,gis_dyn% dyn_ls_a       ,	&
                   gis_dyn% syn_gr_a_1    ,gis_dyn% pyn_gr_a_1     ,    &
                   gis_dyn% dyn_gr_a_1    ,gis_dyn% anl_gr_a_1     ,    &
                   gis_dyn% syn_gr_a_2    ,gis_dyn% pyn_gr_a_2     ,    &
                   gis_dyn% dyn_gr_a_2    ,gis_dyn% anl_gr_a_2     ,    &
                   gis_dyn% pwat          ,gis_dyn% ptot           ,	&
                   gis_dyn% ptrc          ,                             &
                   gis_dyn% pdryini       ,gis_dyn% nblck          ,	&
                   gis_dyn% zhour         ,				&
                   gis_dyn% n1            ,gis_dyn% n4             ,	&
                   gis_dyn% lsout	      ,gis_dyn% ldfi           ,&
                   gis_dyn% colat1        ,gis_dyn% cfhour1	       ,&
                   gis_dyn% start_step    ,gis_dyn% restart_step   ,    &
                   gis_dyn% reset_step    ,gis_dyn% end_step       ,    &
                   gis_dyn% dfiend_step   ,                             &
                   gis_dyn% nfcstdate7,                                 &
                   gis_dyn%Cpl_flag)

          IF(gis_dyn%end_step) THEN
             Cpl_flag1 = .true.
             CALL ESMF_AttributeSet(imp_gfs_dyn, 'Cpl_flag',            &
                                    Cpl_flag1, rc = rc1)
          ELSE
             Cpl_flag1 = .false.
             CALL ESMF_AttributeSet(imp_gfs_dyn, 'Cpl_flag',            &
                                    Cpl_flag1, rc = rc1)
          END IF
!
! ======================================================================
!                     do one time step with two-loop
! ---------------------------------------------------------------------
        else if( gis_dyn% spectral_loop == 2 ) then
!
          call  do_dynamics_two_loop(    gis_dyn% deltim         ,	&
                   gis_dyn% kdt           ,gis_dyn% phour          ,	&
                   gis_dyn% trie_ls       ,gis_dyn% trio_ls        ,	&
                   gis_dyn% grid_gr(1,1,1),gis_dyn% grid_gr_dfi    ,    &
                   gis_dyn% ls_node       ,gis_dyn% ls_nodes       ,	&
                   gis_dyn% max_ls_nodes  ,				&
                   gis_dyn% lats_nodes_a  ,gis_dyn% global_lats_a  ,	&
                   gis_dyn% lonsperlat    ,				&
                   gis_dyn% lats_nodes_ext,gis_dyn% global_lats_ext,	&
                   gis_dyn% epse          ,gis_dyn% epso           ,	&
                   gis_dyn% epsedn        ,gis_dyn% epsodn     ,	&
                   gis_dyn% snnp1ev       ,gis_dyn% snnp1od        ,	&
                   gis_dyn% ndexev        ,gis_dyn% ndexod         ,	&
                   gis_dyn% plnev_a       ,gis_dyn% plnod_a        ,	&
                   gis_dyn% pddev_a       ,gis_dyn% pddod_a        ,	&
                   gis_dyn% plnew_a       ,gis_dyn% plnow_a        ,	&
                   gis_dyn% syn_ls_a      ,gis_dyn% dyn_ls_a       ,	&
                   gis_dyn% syn_gr_a_1    ,gis_dyn% pyn_gr_a_1     ,    &
                   gis_dyn% dyn_gr_a_1    ,gis_dyn% anl_gr_a_1     ,    &
                   gis_dyn% syn_gr_a_2    ,gis_dyn% pyn_gr_a_2     ,    &
                   gis_dyn% dyn_gr_a_2    ,gis_dyn% anl_gr_a_2     ,    &
                   gis_dyn% pwat          ,gis_dyn% ptot           ,	&
                   gis_dyn% ptrc          ,                             & !glbsum
                   gis_dyn% pdryini       ,gis_dyn% nblck          ,	&
                   gis_dyn% zhour         ,				&
                   gis_dyn% n1            ,gis_dyn% n4             ,	&
                   gis_dyn% lsout         ,gis_dyn% ldfi           ,    & ! jw
                   gis_dyn% colat1        ,gis_dyn% cfhour1	,	&
                   gis_dyn%start_step     ,gis_dyn% restart_step   ,    &
                   gis_dyn%reset_step     ,gis_dyn% end_step       ,    &
                   gis_dyn% dfiend_step   ,                             & ! jw
                   gis_dyn% nfcstdate7    ,gis_dyn%iniauinterval)

        else
          print *,' number of spectral loop is wrong. it is ',	&
                    gis_dyn%spectral_loop
          call mpi_quit(99)
          stop

        endif

      END IF ! for the semilag if.
!
! =======================================================================
!                   end of one- or two-loop computation
! =======================================================================
!
! check whether wind speed is out of bound.
!
        do k=1,levs
          if(spdmax(k) < 0. .or. spdmax(k) > 2000.) then
            print *,'unphysical maximum speed',spdmax(k),' me=',me,' k=',k
            call mpi_quit(7)
            stop
          endif
        enddo

! =======================================================================
!   write out random patterns at middle of DFI window
       if( gis_dyn%ldfi .and. gis_dyn%kdt == gis_dyn%ndfi/2 )  then
          if (do_sppt .OR. do_shum .OR. do_skeb .OR. do_vc ) then
             call dump_patterns('stoch_dfi')
          endif
       endif
! update hour
       if( gis_dyn%ndfi>0 .and. gis_dyn%kdt==gis_dyn%ndfi/2+1 .and. .not.gis_dyn%ldfi )  then
          if (do_sppt .OR. do_shum .OR. do_skeb .OR. do_vc ) then
             call restore_patterns('stoch_dfi')
          endif
       endif
!
      gis_dyn%phour = fhour

! ---------------------------------------------------------------------
! prepare n+1 grid value back to common temperature and pressure
!
     IF(gis_dyn%ENS .AND. gis_dyn%end_step) THEN
          call model_to_common_vars (gis_dyn% grid_gr(1,1,gis_dyn%g_qm ),   &
                                     gis_dyn% grid_gr(1,1,gis_dyn%g_ttm),   &
                                     gis_dyn% grid_gr(1,1,gis_dyn%g_rm ),   &
                                     gis_dyn% grid_gr(1,1,gis_dyn%g_uum),   &
                                     gis_dyn% grid_gr(1,1,gis_dyn%g_vvm),   &
                                     gis_dyn% grid_gr(1,1,gis_dyn%g_dpm),   &
                                     gis_dyn% grid_gr(1,1,gis_dyn%g_p  ),   &
                                     gis_dyn% grid_gr(1,1,gis_dyn%g_dpdt ), &
                                     gis_dyn%global_lats_a,                 &
                                     gis_dyn%lonsperlat,0)
          call model_to_common_vars (gis_dyn% grid_gr(1,1,gis_dyn%g_q  ),   &
                                     gis_dyn% grid_gr(1,1,gis_dyn%g_tt ),   &
                                     gis_dyn% grid_gr(1,1,gis_dyn%g_rq ),   &
                                     gis_dyn% grid_gr(1,1,gis_dyn%g_uu ),   &
                                     gis_dyn% grid_gr(1,1,gis_dyn%g_vv ),   &
                                     gis_dyn% grid_gr(1,1,gis_dyn%g_dp ),   &
                                     gis_dyn% grid_gr(1,1,gis_dyn%g_p  ),   &
                                     gis_dyn% grid_gr(1,1,gis_dyn%g_dpdt ), &
                                     gis_dyn%global_lats_a,                 &
                                     gis_dyn%lonsperlat,1)
     ELSE

           call model_to_common_vars (gis_dyn% grid_gr(1,1,gis_dyn%g_zqp),   &
                                      gis_dyn% grid_gr(1,1,gis_dyn%g_ttp),   &
                                      gis_dyn% grid_gr(1,1,gis_dyn%g_rqp),   &
                                      gis_dyn% grid_gr(1,1,gis_dyn%g_uup),   &
                                      gis_dyn% grid_gr(1,1,gis_dyn%g_vvp),   &
                                      gis_dyn% grid_gr(1,1,gis_dyn%g_dpp),   &
                                      gis_dyn% grid_gr(1,1,gis_dyn%g_p  ),   &
                                      gis_dyn% grid_gr(1,1,gis_dyn%g_dpdt ), &
                                      gis_dyn%global_lats_a,                 &
                                      gis_dyn%lonsperlat,1)

      END IF
!cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

      if(present(rc)) then
          rc = rc1
      end if

      end subroutine gfs_dynamics_run

      end module gfs_dynamics_run_mod
