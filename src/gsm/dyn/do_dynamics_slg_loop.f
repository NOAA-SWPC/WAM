      SUBROUTINE do_dynamics_slg_loop(deltim,kdt,PHOUR,
     &                 TRIE_LS,TRIO_LS,GRID_GR,
     &                 LS_NODE,LS_NODES,MAX_LS_NODES,
     &                 LATS_NODES_A,GLOBAL_LATS_A,
     &                 LONSPERLAT,ldfi,grid_gr_dfi,
     &                 EPSE,EPSO,EPSEDN,EPSODN,
     &                 SNNP1EV,SNNP1OD,NDEXEV,NDEXOD,
     &                 PLNEV_A,PLNOD_A,PDDEV_A,PDDOD_A,
     &                 PLNEW_A,PLNOW_A,
     &                 LSOUT,CFHOUR1,SPS,
!    &                 SYN_GR_A_1, SYN_GR_A_2,
!    &                 ANL_GR_A_1, ANL_GR_A_2,
     &                 start_step,restart_step,reset_step,end_step,
     &                 dfiend_step, pdryini,nblck,ZHOUR,
     &                 pwat,ptot,ptrc,nfcstdate7,iniauinterval)
!!
      use gfs_dyn_machine     , only : kind_evod, kind_grid
      use gfs_dyn_resol_def   , only : latg,latg2,levh,levs,levp1,
     &                                 ntrac,ncld,lonf,
     &                                 lotgr,lotgr6,lotls,
     &                                 p_di,p_dim,p_q,p_qm,p_rm,p_rq,
     &                                 p_rt,p_te,p_tem,p_uln,p_vln,
     &                                 p_w,p_x,p_y,p_ze,p_zem,p_zq,
     &                                 adiabatic,        p_dpn, p_dp, 
!    &                                 adiabatic, lonfx, p_dpn, p_dp, 
     &                                 p_dpm, lots, lots_slg, lota,
     &                                 g_q,g_u,g_v,g_t,g_rt

      use gfs_dyn_layout1     , only : len_trie_ls,len_trio_ls,
     &                                 ls_dim,ls_max_node,
     &                                 me,me_l_0,nodes,lats_dim_a,
     &                                 lon_dim_a,
     &                                 ipt_lats_node_a,lats_node_a,
     &                                 lats_node_a_max, lats_dim_ext
      use gfs_dyn_vert_def,     only : am,bm,si,sl,sv,tov,del
      use gfs_dyn_date_def    , only : fhour,idate,shour,spdmax
      use gfs_dyn_gg_def      , only : rcs2_a
      use namelist_dynamics_def,only : ens_nam, mass_dp,
     &                                 gen_coord_hybrid, gg_tracers,
     &                                 hybrid, igen, explicit,
     &                                 nsres, fhdfi, ldfi_spect,
     &                                 sl_epsln, nsout, filta
     &,                                process_split, nemsio_out, fhrot

      use gfs_dyn_mpi_def     , only : kind_mpi,
     &                                 mc_comp,mpi_r_mpi,num_pes_fcst
!     USE module_gfs_mpi_def,     ONLY: num_pes_fcst
      USE gfs_dyn_dfi_mod,        ONLY: gfs_dfi_grid_gr
      USE gfs_dyn_coordinate_def, ONLY: ak5, bk5
      USE gfs_dyn_gg_def,         ONLY: wgt_a
      USE gfs_dyn_bfilt_def,      ONLY: bfilte, bfilto
      USE gfs_dyn_tracer_config,  ONLY: glbsum
      USE do_dynamics_mod,        ONLY: do_dynamics_spectdfini_slg,
     &                                  do_dynamics_syn2gridn_slg,
     &                                  do_dynamics_gridomega_slg,
     &                                  do_dynamics_gridp2n,
     &                                  do_dynamics_spectc2n,
     &                                  do_dynamics_spectn2m,
     &                                  do_dynamics_gridpdp,
     &                                  do_dynamics_gridn2p,
     &                                  do_dynamics_gridn2anl_slg,
     &                                  do_dynamics_gridn2anl_slg_dfi
      use layout_grid_tracers ,  only : rgt_a,rgt_h,xhalo,yhalo

      IMPLICIT NONE
!!     
      integer lat
      INTEGER, SAVE                        :: ndfih, kdtdfi
      CHARACTER(16)                        :: CFHOUR1
      INTEGER,INTENT(IN), dimension(latg)  :: LONSPERLAT, GLOBAL_LATS_A
      LOGICAL,INTENT(IN)                   :: iniauinterval
!!     
      INTEGER,INTENT(IN)                   :: nfcstdate7(7)
      REAL(KIND=KIND_EVOD),  INTENT(in)    :: deltim,PHOUR
      REAL(KIND=KIND_EVOD),INTENT(INOUT)   :: ZHOUR
!
      type(gfs_dfi_grid_gr),intent(inout) :: grid_gr_dfi
      logical,intent(in) :: ldfi

      REAL(KIND=KIND_GRID) GRID_GR(lonf*lats_node_a_max,lotgr)

      REAL(KIND = kind_evod) :: trie_ls_rqt(len_trie_ls,2,levs)
      REAL(KIND = kind_evod) :: trio_ls_rqt(len_trio_ls,2,levs)
      REAL(KIND = kind_evod) :: trie_ls_sfc(len_trie_ls,2)
      REAL(KIND = kind_evod) :: trio_ls_sfc(len_trio_ls,2)
      REAL(KIND = kind_grid) :: typdel(levs)
!     REAL(KIND = kind_grid) :: SYN_GR_A_1(LONFX*LOTS, LATS_DIM_EXT)
!     REAL(KIND = kind_grid) :: SYN_GR_A_2(LONFX*LOTS, LATS_DIM_EXT)
      real(kind=kind_grid), dimension(lon_dim_a,lota,lats_dim_a) ::
     &                          syn_gr_a_1,anl_gr_a_1
      real(kind=kind_grid), dimension(lonf,lota,lats_dim_a) ::
     &                          syn_gr_a_2,anl_gr_a_2
!     REAL(KIND = kind_grid) :: ANL_GR_A_1(LONFX*LOTA, LATS_DIM_EXT)
!     REAL(KIND = kind_grid) :: ANL_GR_A_2(LONFX*LOTA, LATS_DIM_EXT)

      REAL(KIND = kind_grid) :: ptotj(lats_node_a), pwatj(lats_node_a)
      REAL(KIND = kind_grid) :: ptotg(latg), pwatg(latg), typical_pgr
      REAL(KIND = kind_grid) :: sumwa, sumto, ptotp, pwatp, pdryg, 
     &                          pdryini, pcorr
!    &                          pdryini, pcorr, filtb

      REAL(KIND = kind_grid) :: ptrc(lonf,lats_node_a,ntrac)                !glbsum
      REAL(KIND = kind_grid) :: ptrcj(lats_node_a,ntrac)                    !glbsum
      REAL(KIND = kind_grid) :: tmpj(lats_node_a)                           !glbsum
      REAL(KIND = kind_grid) :: ptrcg(latg),sumtrc(ntrac),ptrcp(ntrac)      !glbsum

      REAL(KIND = kind_grid), dimension(lonf,lats_node_a) :: ptot, pwat

      integer ifirst
      data ifirst /1/
      save ifirst
!
      real, allocatable, save   :: gzie_ln(:,:),gzio_ln(:,:)

!     real, allocatable   :: gzie_ln(:,:),gzio_ln(:,:),factor_b2t_ref(:)
!     save gzie_ln,gzio_ln,factor_b2t_ref

      INTEGER NBLCK
      REAL(KIND=KIND_EVOD) TRIE_LS(LEN_TRIE_LS,2,lotls)
     &,                    TRIO_LS(LEN_TRIO_LS,2,lotls)
     &,                    sum_k_rqchange_ls(len_trie_ls,2)
     &,                    sum_k_rqchango_ls(len_trio_ls,2)

!!
      real(kind=kind_evod), allocatable, save :: save_qe_ls(:,:)
     &,                                          save_qo_ls(:,:)
!!
      integer              ls_node(ls_dim,3)
!!
!     ls_node(1,1) ... ls_node(ls_max_node,1) : values of L
!     ls_node(1,2) ... ls_node(ls_max_node,2) : values of jbasev
!     ls_node(1,3) ... ls_node(ls_max_node,3) : values of jbasod
!!
      INTEGER              LS_NODES(LS_DIM,NODES)
      INTEGER, dimension(nodes)  :: MAX_LS_NODES, LATS_NODES_A
!
      real(kind=kind_evod), dimension(len_trie_ls) :: epse, epsedn,
     &                                                snnp1ev
      real(kind=kind_evod), dimension(len_trio_ls) :: epso, epsodn,
     &                                                snnp1od
      real(kind=kind_evod), dimension(len_trie_ls,latg2) :: plnev_a,
     &                                                pddev_a, plnew_a
      real(kind=kind_evod), dimension(len_trie_ls,latg2) :: plnod_a,
     &                                                pddod_a, plnow_a
      INTEGER NDEXEV(LEN_TRIE_LS), NDEXOD(LEN_TRIO_LS)
!
      INTEGER kdt, lan, IERR, I, J, K, L, LOCL, N, item, jtem, 
     &        ltem, ktem, lons_lat, mtem, jt, kt,kk
     &,       kp_t,kp_y,kp_d,kp_z,kp_x,kp_w, stp, nt, iprint
     &,       INDLSEV,JBASEV,INDLSOD,JBASOD
      real(kind=kind_evod)  batah, tem
      REAL(KIND=kind_mpi)  coef00m(LEVS,ntrac)! temp. ozone clwater  
      REAL(KIND=kind_evod) coef00(LEVS,ntrac) ! temp. ozone clwater  
      REAL(KIND=kind_evod), parameter :: CONS0P5 = 0.5, CONS1=1.0, 
     &                                   CONS2=2.0
      include 'function_indlsev'

      LOGICAL LSOUT, SPS
!     LOGICAL, SAVE      :: fwd_step
      LOGICAL            :: start_step,   reset_step, end_step,
     &                      restart_step, dfiend_step
!     LOGICAL, PARAMETER :: ladj = .false.
      LOGICAL, PARAMETER :: ladj = .true.
!
!
! timings
!     real(kind=kind_evod) global_times_a(latg,nodes)

!     real*8               rtc ,timer1,timer2
!
      shour = kdt * deltim
      fhour = shour/3600.
      batah = 1.0 + sl_epsln       ! Moorthi

!     if (me == 0) write(0,*)' in do_dynamics_slg_loop shour=',shour

!     filtb = (cons1-filta)*cons0p5
!
!----------------------------------------------------------
!**********************************************************ME<NUM_PES_FCST*
!      if (me < num_pes_fcst) then
!----------------------------------------------------------
!
        if (.not. allocated(save_qe_ls))
     &        allocate (save_qe_ls(len_trie_ls,2))
        if (.not. allocated(save_qo_ls))
     &        allocate (save_qo_ls(len_trio_ls,2))
!        if(me==0) print *,'kdt=',kdt,'kdtdfi=',kdtdfi,
!     &  'reset_step=',reset_step

        if( start_step ) then

          call get_cd_hyb_slg(deltim,batah)

! --------------------------------------------------------------
! if the first step, from internal state, prepare syn for nonlinear
! -------- this section is called once only ---------
! --------------------------------------------------------------

!         ndfih  = nint(fhdfi*3600./deltim)
!         kdtdfi = kdt + ndfih

          if (me == 0) 
     &     print *,' start step from internal state (grid and spectral)'
     &,            ' ldfi_spect=',ldfi_spect

! dfi.
          if (ldfi_spect) then
            ndfih  = nint(fhdfi*3600./deltim)
            kdtdfi = kdt + ndfih
 
            if ( me == 0 ) print *,' calling spectdfini with ndfih=',
     &                             ndfih
            if (ndfih /= 0 ) then
              call do_dynamics_spectdfini_slg(-ndfih-1,ndfih,trie_ls,
     &                                     trio_ls)
            endif
          endif 

          start_step = .false.
! -------------------------------------------------------
        else if( reset_step ) then
! --------------------------------------------------------------
! if it is reset step to reset internal state to be import state
! -------- this section is called once only ---------
! --------------------------------------------------------------
!         print *,' reset internal values by import for all '
!         if (me == 0) write(0,*)' in reset_step option'

          if (ldfi_spect) then
            if( me == 0 ) print *,' finialize spectdfini '
            call do_dynamics_spectdfini_slg(ndfih+1,ndfih,
     &                                      trie_ls,trio_ls)
            call do_dynamics_spectc2n(trie_ls,trio_ls)
            call do_dynamics_spectn2m(trie_ls,trio_ls)
            ldfi_spect = .false.
          else
! move data from physics to n+1
            lsout      = .true.
            call do_dynamics_gridp2n(grid_gr,global_lats_a,lonsperlat)

            call do_dynamics_gridn2anl_slg_dfi(grid_gr,anl_gr_a_2,rcs2_a
     &,                                        global_lats_a,lonsperlat)
            call grid_to_spect_slg(anl_gr_a_1,anl_gr_a_2
     &,                            trie_ls,trio_ls,lsout
     &,                            ls_node,ls_nodes,max_ls_nodes
     &,                            lats_nodes_a,global_lats_a,lonsperlat
     &,                            epse,epso,plnew_a,plnow_a)
!
            do i=1,len_trie_ls
!             trie_ls(i,1,p_q)  = trie_ls(i,1,p_zq)
!             trie_ls(i,2,p_q)  = trie_ls(i,2,p_zq)
              trie_ls(i,1,p_qm) = trie_ls(i,1,p_q)
              trie_ls(i,2,p_qm) = trie_ls(i,2,p_q)
            enddo
            do i=1,len_trio_ls
!             trio_ls(i,1,p_q)  = trio_ls(i,1,p_zq)
!             trio_ls(i,2,p_q)  = trio_ls(i,2,p_zq)
              trio_ls(i,1,p_qm) = trio_ls(i,1,p_q)
              trio_ls(i,2,p_qm) = trio_ls(i,2,p_q)
            enddo
!
!$omp parallel do private(k,i)
            do k=1,levs
              do i=1,len_trie_ls
!               trie_ls(i,1,p_di+k-1) = trie_ls(i,1,p_x+k-1)
!               trie_ls(i,1,p_ze+k-1) = trie_ls(i,1,p_w+k-1)
!               trie_ls(i,1,p_te+k-1) = trie_ls(i,1,p_y+k-1)

!               trie_ls(i,2,p_di+k-1) = trie_ls(i,2,p_x+k-1)
!               trie_ls(i,2,p_ze+k-1) = trie_ls(i,2,p_w+k-1)
!               trie_ls(i,2,p_te+k-1) = trie_ls(i,2,p_y+k-1)

                trie_ls(i,1,p_dim+k-1) = trie_ls(i,1,p_di+k-1)
                trie_ls(i,1,p_zem+k-1) = trie_ls(i,1,p_ze+k-1)
                trie_ls(i,1,p_tem+k-1) = trie_ls(i,1,p_te+k-1)

                trie_ls(i,2,p_dim+k-1) = trie_ls(i,2,p_di+k-1)
                trie_ls(i,2,p_zem+k-1) = trie_ls(i,2,p_ze+k-1)
                trie_ls(i,2,p_tem+k-1) = trie_ls(i,2,p_te+k-1)
              enddo
              do i=1,len_trio_ls
!               trio_ls(i,1,p_di+k-1) = trio_ls(i,1,p_x+k-1)
!               trio_ls(i,1,p_ze+k-1) = trio_ls(i,1,p_w+k-1)
!               trio_ls(i,1,p_te+k-1) = trio_ls(i,1,p_y+k-1)

!               trio_ls(i,2,p_di+k-1) = trio_ls(i,2,p_x+k-1)
!               trio_ls(i,2,p_ze+k-1) = trio_ls(i,2,p_w+k-1)
!               trio_ls(i,2,p_te+k-1) = trio_ls(i,2,p_y+k-1)

                trio_ls(i,1,p_dim+k-1) = trio_ls(i,1,p_di+k-1)
                trio_ls(i,1,p_zem+k-1) = trio_ls(i,1,p_ze+k-1)
                trio_ls(i,1,p_tem+k-1) = trio_ls(i,1,p_te+k-1)

                trio_ls(i,2,p_dim+k-1) = trio_ls(i,2,p_di+k-1)
                trio_ls(i,2,p_zem+k-1) = trio_ls(i,2,p_ze+k-1)
                trio_ls(i,2,p_tem+k-1) = trio_ls(i,2,p_te+k-1)
              enddo
            enddo
!!!!!!!!!!!!!!!!!!!!!!!
!           do k=1,levs
!             trie_ls(:,:,p_w+k-1)   = trie_ls(:,:,p_ze+k-1)
!             trie_ls(:,:,p_x+k-1)   = trie_ls(:,:,p_di+k-1)
!             trie_ls(:,:,p_y+k-1)   = trie_ls(:,:,p_te+k-1)
!             trie_ls(:,:,p_dpn+k-1) = trie_ls(:,:,p_dp+k-1)
!             trio_ls(:,:,p_w+k-1)   = trio_ls(:,:,p_ze+k-1)
!             trio_ls(:,:,p_x+k-1)   = trio_ls(:,:,p_di+k-1)
!             trio_ls(:,:,p_y+k-1)   = trio_ls(:,:,p_te+k-1)
!             trio_ls(:,:,p_dpn+k-1) = trio_ls(:,:,p_dp+k-1)
!           enddo
!!!!!!!!!!!!!!!!!!!!!!!
            if (gg_tracers) then
              do lan=1,lats_node_a
                lat      = global_lats_a(ipt_lats_node_a-1+lan)
                lons_lat = lonsperlat(lat)
                do n=1,ntrac
!$omp parallel do private(i,k,kk)
                  do k=1,levs
                    kk = levs + 1 - k
                    do i=1,min(lonf,lons_lat)
                      rgt_a(i,kk,lan,n) = rgt_h(i+xhalo,k,lan+yhalo,n)
                    enddo
                  enddo
                enddo
              enddo
            else
!$omp parallel do private(k,i)
              do k=1,levh
                do i=1,len_trie_ls
                  trie_ls(i,1,p_rm+k-1) = trie_ls(i,1,p_rq+k-1)
                  trie_ls(i,2,p_rm+k-1) = trie_ls(i,2,p_rq+k-1)
                enddo
                do i=1,len_trio_ls
                  trio_ls(i,1,p_rm+k-1) = trio_ls(i,1,p_rq+k-1)
                  trio_ls(i,2,p_rm+k-1) = trio_ls(i,2,p_rq+k-1)
                enddo
              enddo
            endif
!----------------------------------------------------------
!
            if (.not. iniauinterval) then ! skip if in IAU forcing interval
              do i=1,len_trie_ls
                sum_k_rqchange_ls(i,1) = trie_ls(i,1,p_q)
                sum_k_rqchange_ls(i,2) = trie_ls(i,2,p_q)
              enddo
              do i=1,len_trio_ls
                sum_k_rqchango_ls(i,1) = trio_ls(i,1,p_q)
                sum_k_rqchango_ls(i,2) = trio_ls(i,2,p_q)
              enddo
  
              do i=1,len_trie_ls
                trie_ls(i,1,p_q) = save_qe_ls(i,1)
                trie_ls(i,2,p_q) = save_qe_ls(i,2)
              enddo
              do i=1,len_trio_ls
                trio_ls(i,1,p_q) = save_qo_ls(i,1)
                trio_ls(i,2,p_q) = save_qo_ls(i,2)
              enddo
            endif
!

!$omp parallel do private(k,kp_d,kp_z,kp_t,kp_x,kp_w,kp_y,j)
            do k=1,levs
              kp_d = P_DI + K - 1
              kp_z = P_ZE + K - 1
              kp_t = P_TE + K - 1
              kp_x = P_X  + K - 1
              kp_w = P_W  + K - 1
              kp_y = P_Y  + K - 1
              do j=1,len_trie_ls
                trie_ls(j,1,kp_d) = trie_ls(j,1,kp_x)
                trie_ls(j,2,kp_d) = trie_ls(j,2,kp_x)
                trie_ls(j,1,kp_z) = trie_ls(j,1,kp_w)
                trie_ls(j,2,kp_z) = trie_ls(j,2,kp_w)
                trie_ls(j,1,kp_t) = trie_ls(j,1,kp_y)
                trie_ls(j,2,kp_t) = trie_ls(j,2,kp_y)
              enddo
              do j=1,len_trio_ls
                trio_ls(j,1,kp_d) = trio_ls(j,1,kp_x)
                trio_ls(j,2,kp_d) = trio_ls(j,2,kp_x)
                trio_ls(j,1,kp_z) = trio_ls(j,1,kp_w)
                trio_ls(j,2,kp_z) = trio_ls(j,2,kp_w)
                trio_ls(j,1,kp_t) = trio_ls(j,1,kp_y)
                trio_ls(j,2,kp_t) = trio_ls(j,2,kp_y)
              enddo
            enddo

            if(.not. gg_tracers)then
!$omp parallel do private(k,kp_x,kp_t,j)
              do k=1,levh
                kp_x = P_RQ + K - 1
                kp_t = P_RT + K - 1
                do j=1,len_trie_ls
                  trie_ls(j,1,kp_x) = trie_ls(j,1,kp_t)
                  trie_ls(j,2,kp_x) = trie_ls(j,2,kp_t)
                enddo
                do j=1,len_trio_ls
                  trio_ls(j,1,kp_x) = trio_ls(j,1,kp_t)
                  trio_ls(j,2,kp_x) = trio_ls(j,2,kp_t)
                enddo
              enddo
            endif !  if(.not.gg_tracers

          endif

          reset_step = .false.
! --------------------------------------------------------------
        else if( restart_step ) then   ! restart from history file
! --------------------------------------------------------------
          if(me == 0) print *,'in restart step'

!         call get_cd_hyb(deltim)
          call get_cd_hyb_slg(deltim,batah)

          if (ldfi_spect) then
            ndfih  = nint(fhdfi*3600./deltim)
            kdtdfi = kdt + ndfih

            if ( me == 0 ) print *,' calling spectdfini with ndfih=',
     &        ndfih
            if (ndfih /= 0 ) then
              call do_dynamics_spectdfini_slg(-ndfih-1,ndfih,trie_ls,
     &                                     trio_ls)
            endif
          endif


!    THIS PART OF THE CODE IS IMCOMPLETE - RESTART needs UPDATING
!
!?         call spect_to_grid(trie_ls,trio_ls,
!?    &                       syn_gr_a_1,syn_gr_a_2,
!?    &                       ls_node,ls_nodes,max_ls_nodes,
!?    &                       lats_nodes_a,global_lats_a,lonsperlat,
!?    &                       epse,epso,epsedn,epsodn,
!?    &                       snnp1ev,snnp1od,plnev_a,plnod_a)

! ------------------------------------------------------
        else    ! end start_step, begin not start_step
! ------------------------------------------------------
! start linear computation
! -----------------------------------------------------------
!
!******************************************************************ADIABATIC
          if( .not. adiabatic ) then    ! logical variable in gfs_dyn_resol_def
! --------------------------------
! do after physics, not from input
! -----------------------------------------------------------
! update after physics

!         if(me==0) print *,'normal dyn step, kdt=',kdt
! move data from physics to n+1
            call do_dynamics_gridp2n(grid_gr,global_lats_a,lonsperlat)
!
            call do_dynamics_gridn2anl_slg(grid_gr,anl_gr_a_2
     &,                                    rcs2_a
     &,                                    global_lats_a,lonsperlat
     &,                                    iniauinterval)
!
! transform values in grid to spectral
!
            call grid_to_spect_slg(anl_gr_a_1,anl_gr_a_2
     &,                            trie_ls,trio_ls,lsout
     &,                            ls_node,ls_nodes,max_ls_nodes
     &,                            lats_nodes_a,global_lats_a,lonsperlat
     &,                            epse,epso,plnew_a,plnow_a)
!
!----------------------------------------------------------
!
          if (.not. iniauinterval) then ! skip if in IAU forcing interval
            do i=1,len_trie_ls
              sum_k_rqchange_ls(i,1) = trie_ls(i,1,p_q)
              sum_k_rqchange_ls(i,2) = trie_ls(i,2,p_q)
            enddo
            do i=1,len_trio_ls
              sum_k_rqchango_ls(i,1) = trio_ls(i,1,p_q)
              sum_k_rqchango_ls(i,2) = trio_ls(i,2,p_q)
            enddo

            do i=1,len_trie_ls
              trie_ls(i,1,p_q)       = save_qe_ls(i,1)
              trie_ls(i,2,p_q)       = save_qe_ls(i,2)
            enddo
            do i=1,len_trio_ls
              trio_ls(i,1,p_q)       = save_qo_ls(i,1)
              trio_ls(i,2,p_q)       = save_qo_ls(i,2)
            enddo
         endif


!     print *,' ----- do bfilter ------ '

!$omp parallel do private(k,i,jtem,ktem)
            do k=1,levs
              ktem = p_w   + k - 1
              jtem = p_vln + k - 1
              do i=1,len_trie_ls
                trie_ls(i,1,ktem) = trie_ls(i,1,ktem) + bfilte(i)*
     &                             (trie_ls(i,1,jtem)-trie_ls(i,1,ktem))

                trie_ls(i,2,ktem) = trie_ls(i,2,ktem) + bfilte(i)*
     &                             (trie_ls(i,2,jtem)-trie_ls(i,2,ktem))
              enddo
              do i=1,len_trio_ls
                trio_ls(i,1,ktem) = trio_ls(i,1,ktem) + bfilto(i)*
     &                             (trio_ls(i,1,jtem)-trio_ls(i,1,ktem))

                trio_ls(i,2,ktem) = trio_ls(i,2,ktem) + bfilto(i)*
     &                             (trio_ls(i,2,jtem)-trio_ls(i,2,ktem))
              enddo
            enddo

            if(.not. gg_tracers)then
!$omp parallel do private(k,i,jtem,ktem,tem)
              do k=1,levs
                ktem = p_rt + k - 1
                jtem = p_rq + k - 1

                do i=1,len_trie_ls
                  tem = bfilte(i)*(trie_ls(i,1,jtem)-trie_ls(i,1,ktem))
                  trie_ls_rqt(i,1,k) = tem
                  trie_ls(i,1,ktem)  = trie_ls(i,1,ktem) + tem
!
                  tem = bfilte(i)*(trie_ls(i,2,jtem)-trie_ls(i,2,ktem))
                  trie_ls_rqt(i,2,k) = tem
                  trie_ls(i,2,ktem)  = trie_ls(i,2,ktem) + tem

!                 trie_ls_rqt(i,1,k) = bfilte(i)*
!    &                               (trie_ls(i,1,item)-trie_ls(i,1,jtem))
!                 trie_ls_rqt(i,2,k) = bfilte(i)*
!    &                            (trie_ls(i,2,item)-trie_ls(i,2,jtem))
!                 trie_ls(i,1,ktem)  = trie_ls(i,1,ltem) +
!    &                     bfilte(i) *(trie_ls(i,1,ktem)-trie_ls(i,1,ltem))
!                 trie_ls(i,2,ktem)  = trie_ls(i,2,ltem) +
!    &                     bfilte(i) *(trie_ls(i,2,ktem)-trie_ls(i,2,ltem))
                enddo
!!
                do i=1,len_trio_ls
                  tem = bfilto(i)*(trio_ls(i,1,jtem)-trio_ls(i,1,ktem))
                  trio_ls_rqt(i,1,k)   = tem
                  trio_ls(i,1,ktem)     = trio_ls(i,1,ktem) + tem
! 
                  tem = bfilto(i)*(trio_ls(i,2,jtem)-trio_ls(i,2,ktem))
                  trio_ls_rqt(i,2,k)    = tem
                  trio_ls(i,2,p_rt+k-1) = trio_ls(i,2,ktem) + tem

!                 trio_ls_rqt(i,1,k) = bfilto(i)*
!    &                                (trio_ls(i,1,item)-trio_ls(i,1,jtem))
!                 trio_ls_rqt(i,2,k) = bfilto(i)*
!    &                                (trio_ls(i,2,item)-trio_ls(i,2,jtem))
!                 trio_ls(i,1,ktem)  = trio_ls(i,1,ltem) +
!    &                     bfilto(i) *(trio_ls(i,1,ktem)-trio_ls(i,1,ltem))
!                 trio_ls(i,2,ktem)  = trio_ls(i,2,ltem) +
!    &                     bfilto(i) *(trio_ls(i,2,ktem)-trio_ls(i,2,ltem))
                enddo
              enddo
!!
              do nt=2,ntrac
                ktem = p_rt + (nt-1)*levs - 1
                jtem = p_rq + (nt-1)*levs - 1
!             do k=levs*(nt-2)+1,levs*(nt-1)
!$omp parallel do private(k,i,jt,kt)
                do k=1,levs
                  kt = ktem + k
                  jt = jtem + k
                  do i=1,len_trie_ls
                    trie_ls(i,1,kt) = trie_ls(i,1,kt) + bfilte(i)*
     &                               (trie_ls(i,1,jt)-trie_ls(i,1,kt))
                    trie_ls(i,2,kt) = trie_ls(i,2,kt) + bfilte(i)*
     &                               (trie_ls(i,2,jt)-trie_ls(i,2,kt))
                 enddo
                 do i=1,len_trio_ls
                   trio_ls(i,1,kt) = trio_ls(i,1,kt) + bfilto(i)*
     &                              (trio_ls(i,1,jt)-trio_ls(i,1,kt))
                   trio_ls(i,2,kt) = trio_ls(i,2,kt) + bfilto(i)*
     &                              (trio_ls(i,2,jt)-trio_ls(i,2,kt))
                 enddo
                enddo
              enddo

!             do k=1,levh
!               item = p_rt+k-1
!               jtem = p_rq+k-1
!               do i=1,len_trie_ls
!                trie_ls(i,1,item) =  trie_ls(i,1,jtem) +
!    &                   bfilte(i) * (trie_ls(i,1,item)-trie_ls(i,1,jtem))
!                trie_ls(i,2,item) =  trie_ls(i,2,jtem) +
!    &                   bfilte(i) * (trie_ls(i,2,item)-trie_ls(i,2,jtem))
!               enddo
!               do i=1,len_trio_ls
!                trio_ls(i,1,item) =  trio_ls(i,1,jtem) +
!    &                   bfilto(i) * (trio_ls(i,1,item)-trio_ls(i,1,jtem))
!                trio_ls(i,2,item) =  trio_ls(i,2,jtem) +
!    &                   bfilto(i) * (trio_ls(i,2,item)-trio_ls(i,2,jtem))
!               enddo
!             enddo

            endif  ! if(.not.gg_tracers)
!!
!----------------------------------------------------------------------
!         print *,' ----- do pdry correction ------ '

             if (hybrid) then
               typical_pgr = 85.
               do k=1,levp1
                 si(levs+2-k) = ak5(k)/typical_pgr + bk5(k)
               enddo
             endif

             do k=1,levs
               typdel(k) = si(k)-si(k+1)
             enddo

!----------------------------------------------------------------------
! adjust moisture changes to the total mass to conserve dry mass
!
             do lan=1,lats_node_a
               lat      = global_lats_a(ipt_lats_node_a-1+lan)
               lons_lat = lonsperlat(lat)
               ptotp    = 0.
               pwatp    = 0.
               tem      = 0.5 / lons_lat
               do i=1,lons_lat
                  ptotp = ptotp + ptot(i,lan)
                  pwatp = pwatp + pwat(i,lan)
               enddo
               pwatj(lan) = pwatp * tem
               ptotj(lan) = ptotp * tem
             enddo
             call excha(lats_nodes_a,global_lats_a,ptotj,pwatj,ptotg,
     &                                                         pwatg)
             sumwa = 0.
             sumto = 0.
             do lat=1,latg
               sumto = sumto + wgt_a(min(lat,latg-lat+1))*ptotg(lat)
               sumwa = sumwa + wgt_a(min(lat,latg-lat+1))*pwatg(lat)
             enddo

             pdryg = sumto - sumwa
             if(pdryini <= 0.) pdryini = pdryg

             if ( glbsum ) then                                             !glbsum
               do lan=1,lats_node_a
                 lat      = global_lats_a(ipt_lats_node_a-1+lan)
                 lons_lat = lonsperlat(lat)
                 ptrcp(:) = 0.
                 tem      = 0.5 / lons_lat
                 do n = 1, ntrac
                   do i=1,lons_lat
                     ptrcp(n)   = ptrcp(n) + ptrc(i,lan,n)
                   enddo
                   ptrcj(lan,n) = ptrcp(n) * tem
                 enddo
               enddo
               do n = 1, ntrac
                 sumtrc(n) = 0.
                 tmpj(:)   = ptrcj(:,n)
                 call excha(lats_nodes_a,global_lats_a,ptotj,tmpj,
     &                                                ptotg,ptrcg)
                 do lat=1,latg
                  sumtrc(n) = sumtrc(n)
     &                      + wgt_a(min(lat,latg-lat+1))*ptrcg(lat)
                 enddo
               enddo
             endif                                                          !glbsum

             pcorr = (pdryini-pdryg)/sumto*sqrt(2.)
!
             if (me == 0) write(0,*)'pcorr pdryini pdryg ',pcorr,
     &                    pdryini,pdryg,' kdt=',kdt,' fhour=',fhour

             if (glbsum .and. me == me_l_0) then                            !glbsum
               write(70,111) kdt,fhour,idate                                !glbsum
               write(71,*)   kdt,(sumtrc(n),n=4,ntrac)                      !glbsum
             endif                                                          !glbsum
111          format ('kdt, fhour, idate=',i6,1x,f10.3,2x,4(i4,2x))          !glbsum


!*********************************************************************LADJ

             if (ladj) then
!
               do i=1,len_trie_ls
                 trie_ls(i,1,p_zq) = 0.
                 trie_ls(i,2,p_zq) = 0.
               enddo
               do i=1,len_trio_ls
                 trio_ls(i,1,p_zq) = 0.
                 trio_ls(i,2,p_zq) = 0.
               enddo
               if (me == me_l_0) then
!                trie_ls(1,1,p_zq) = pcorr
                 trie_ls(1,1,p_zq) = trie_ls(1,1,p_zq) + pcorr
               endif

               if(gg_tracers)then
                 do i=1,len_trie_ls
                   trie_ls(i,1,p_zq) = trie_ls(i,1,p_zq)
     &                               + sum_k_rqchange_ls(i,1)
                   trie_ls(i,2,p_zq) = trie_ls(i,2,p_zq)
     &                               + sum_k_rqchange_ls(i,2)
                 enddo
                 do i=1,len_trio_ls
                   trio_ls(i,1,p_zq) = trio_ls(i,1,p_zq)
     &                               + sum_k_rqchango_ls(i,1)
                   trio_ls(i,2,p_zq) = trio_ls(i,2,p_zq)
     &                               + sum_k_rqchango_ls(i,2)
                 enddo
               else
                 do k=1,levs
                   do i=1,len_trie_ls
                     trie_ls(i,1,p_zq) = trie_ls(i,1,p_zq)
     &                                 + typdel(k)*trie_ls_rqt(i,1,k)
                     trie_ls(i,2,p_zq) = trie_ls(i,2,p_zq)
     &                                 + typdel(k)*trie_ls_rqt(i,2,k)
                   enddo
                   do i=1,len_trio_ls
                     trio_ls(i,1,p_zq) = trio_ls(i,1,p_zq)
     &                                 + typdel(k)*trio_ls_rqt(i,1,k)
                     trio_ls(i,2,p_zq) = trio_ls(i,2,p_zq)
     &                                 + typdel(k)*trio_ls_rqt(i,2,k)
                   enddo
                 enddo
               endif               !fin if(gg_tracers)

!!
!$omp parallel do private(k,i,item,jtem,ktem,ltem,mtem)
               do k=1,levs
                 item = p_di+k-1
                 jtem = p_uln+k-1
                 ktem = p_x+k-1
                 ltem = p_te+k-1
                 mtem = p_y+k-1
 
                 do i=1,len_trie_ls
                  trie_ls(i,1,item) = bfilte(i)
     &                          * (trie_ls(i,1,jtem)-trie_ls(i,1,ktem))
                  trie_ls(i,1,ltem) = bfilte(i)
     &                          * (trie_ls(i,1,ltem)-trie_ls(i,1,mtem))
                  trie_ls(i,2,item) = bfilte(i)
     &                          * (trie_ls(i,2,jtem)-trie_ls(i,2,ktem))
                  trie_ls(i,2,ltem) = bfilte(i)
     &                          * (trie_ls(i,2,ltem)-trie_ls(i,2,mtem))
                 enddo
                 do i=1,len_trio_ls
                   trio_ls(i,1,item) = bfilto(i)
     &                           * (trio_ls(i,1,jtem)-trio_ls(i,1,ktem))
                   trio_ls(i,1,ltem) = bfilto(i)
     &                           * (trio_ls(i,1,ltem)-trio_ls(i,1,mtem))
                   trio_ls(i,2,item) = bfilto(i)
     &                           * (trio_ls(i,2,jtem)-trio_ls(i,2,ktem))
                   trio_ls(i,2,ltem) = bfilto(i)
     &                           * (trio_ls(i,2,ltem)-trio_ls(i,2,mtem))
                 enddo
               enddo
!!
!           do k=1,levs
!             item = p_x +k-1
!             jtem = p_di+k-1
!             ktem = p_y +k-1
!             ltem = p_te+k-1
!             do i=1,len_trie_ls
!              trie_ls(i,1,item) =  trie_ls(i,1,jtem) +
!    &                 bfilte(i) * (trie_ls(i,1,item)-trie_ls(i,1,jtem))
!              trie_ls(i,2,item) =  trie_ls(i,2,jtem) +
!    &                 bfilte(i) * (trie_ls(i,2,item)-trie_ls(i,2,jtem))
!              trie_ls(i,1,ktem) =  trie_ls(i,1,ltem) +
!    &                 bfilte(i) * (trie_ls(i,1,ktem)-trie_ls(i,1,ltem))
!              trie_ls(i,2,ktem) =  trie_ls(i,2,ltem) +
!    &                 bfilte(i) * (trie_ls(i,2,ktem)-trie_ls(i,2,ltem))
!             enddo
!             do i=1,len_trio_ls
!              trio_ls(i,1,item) =  trio_ls(i,1,jtem) +
!    &                 bfilto(i) * (trio_ls(i,1,item)-trio_ls(i,1,jtem))
!              trio_ls(i,2,item) =  trio_ls(i,2,jtem) +
!    &                 bfilto(i) * (trio_ls(i,2,item)-trio_ls(i,2,jtem))
!              trio_ls(i,1,ktem) =  trio_ls(i,1,ltem) +
!    &                 bfilto(i) * (trio_ls(i,1,ktem)-trio_ls(i,1,ltem))
!              trio_ls(i,2,ktem) =  trio_ls(i,2,ltem) +
!    &                 bfilto(i) * (trio_ls(i,2,ktem)-trio_ls(i,2,ltem))
!             enddo
!           enddo

!---------------------------------------------------------
!$OMP parallel do private(locl)
               do locl=1,ls_max_node


                 call impadje_slg(trie_ls(1,1,p_x),trie_ls(1,1,p_y),
     &                            trie_ls(1,1,p_q),trie_ls(1,1,p_di),
     &                            trie_ls(1,1,p_te),trie_ls(1,1,p_zq),
     &                            deltim,
     &                            trie_ls(1,1,p_uln),trie_ls(1,1,p_vln),
     &                            snnp1ev,ndexev,ls_node,locl,batah)
!!
                 call impadjo_slg(trio_ls(1,1,p_x),trio_ls(1,1,p_y),
     &                            trio_ls(1,1,p_q),trio_ls(1,1,p_di),
     &                            trio_ls(1,1,p_te),trio_ls(1,1,p_zq),
     &                            deltim,
     &                            trio_ls(1,1,p_uln),trio_ls(1,1,p_vln),
     &                            snnp1od,ndexod,ls_node,locl,batah)
               enddo
!
!!$OMP parallel do private(locl)
!              do locl=1,ls_max_node
!                call impadje_hyb(trie_ls(1,1,p_x),trie_ls(1,1,p_y),
!    &                            trie_ls(1,1,p_zq),trie_ls(1,1,p_di),
!    &                            trie_ls(1,1,p_te),trie_ls(1,1,p_q),
!    &                            deltim ,
!    &                            trie_ls(1,1,p_uln),
!    &                            trie_ls(1,1,p_vln),
!    &                            snnp1ev,ndexev,ls_node,locl)
!!
!                call impadjo_hyb(trio_ls(1,1,p_x),trio_ls(1,1,p_y),
!    &                            trio_ls(1,1,p_zq),trio_ls(1,1,p_di),
!    &                            trio_ls(1,1,p_te),trio_ls(1,1,p_q),
!    &                            deltim ,
!    &                            trio_ls(1,1,p_uln),
!    &                            trio_ls(1,1,p_vln),
!    &                            snnp1od,ndexod,ls_node,locl)
!              enddo

!---------------------------------------------------------

! -----------------------------------------------------------
             else  ! fin impadj,    following is with no impadj
! -----------------------------------------------------------
               DO k=1,LEVS
                del(k) = typDEL(k)                 ! sela 4.5.07
               ENDDO
               if (me == me_l_0) then
                 trie_ls(1,1,p_q) = trie_ls(1,1,p_q) + pcorr
               endif
!
! testing mass correction on sep 25
!!
              if (.not. iniauinterval) then
               if(gg_tracers)then
                 do i=1,len_trie_ls
                   trie_ls(i,1,p_q) = trie_ls(i,1,p_q)
     &                              + sum_k_rqchange_ls(i,1)
                   trie_ls(i,2,p_q) = trie_ls(i,2,p_q)
     &                              + sum_k_rqchange_ls(i,2)
                 enddo
                 do i=1,len_trio_ls
                   trio_ls(i,1,p_q) = trio_ls(i,1,p_q)
     &                              + sum_k_rqchango_ls(i,1)
                   trio_ls(i,2,p_q) = trio_ls(i,2,p_q)
     &                              + sum_k_rqchango_ls(i,2)
                 enddo
               else
                 do k=1,levs
                   do i=1,len_trie_ls
                    trie_ls(i,1,p_q) = trie_ls(i,1,p_q)
     &                               + del(k)*trie_ls_rqt(i,1,k)
                    trie_ls(i,2,p_q) = trie_ls(i,2,p_q)
     &                               + del(k)*trie_ls_rqt(i,2,k)
                   enddo
                   do i=1,len_trio_ls
                      trio_ls(i,1,p_q) = trio_ls(i,1,p_q)
     &                                 + del(k)*trio_ls_rqt(i,1,k)
                      trio_ls(i,2,p_q) = trio_ls(i,2,p_q)
     &                                 + del(k)*trio_ls_rqt(i,2,k)
                   enddo
                 enddo
               endif
              endif
! testing mass correction on sep 25
!
!$omp parallel do private(k,i,item,jtem,ktem,ltem,mtem)
               do k=1,levs
                 item = p_di+k-1
                 jtem = p_uln+k-1
                 ktem = p_x+k-1
                 ltem = p_te+k-1
                 mtem = p_y+k-1
                 do i=1,len_trie_ls
                   trie_ls(i,1,ktem) = trie_ls(i,1,ktem) + bfilte(i)
     &                            *(trie_ls(i,1,jtem)-trie_ls(i,1,ktem))
                   trie_ls(i,2,ktem) = trie_ls(i,2,ktem) + bfilte(i)
     &                            *(trie_ls(i,2,jtem)-trie_ls(i,2,ktem))
                   trie_ls(i,1,mtem) = trie_ls(i,1,mtem) + bfilte(i)
     &                            *(trie_ls(i,1,ltem)-trie_ls(i,1,mtem))
                   trie_ls(i,2,mtem) = trie_ls(i,2,mtem) + bfilte(i)
     &                            *(trie_ls(i,2,ltem)-trie_ls(i,2,mtem))
                 enddo

                 do i=1,len_trio_ls
                   trio_ls(i,1,ktem) = trio_ls(i,1,ktem)+bfilto(i)
     &                            *(trio_ls(i,1,jtem)-trio_ls(i,1,ktem))
                   trio_ls(i,2,ktem) = trio_ls(i,2,ktem) + bfilto(i)
     &                            *(trio_ls(i,2,jtem)-trio_ls(i,2,ktem))
                   trio_ls(i,1,mtem) = trio_ls(i,1,mtem) + bfilto(i)
     &                            *(trio_ls(i,1,ltem)-trio_ls(i,1,mtem))
                   trio_ls(i,2,mtem) = trio_ls(i,2,mtem) + bfilto(i)
     &                            *(trio_ls(i,2,ltem)-trio_ls(i,2,mtem))
                 enddo
               enddo
             endif   ! fin no ladj (i.e. no massadj)


!           if (me.eq.me_l_0) then
!             trie_ls(1,1,p_zq)=trie_ls(1,1,p_zq)+pcorr
!           endif
!!
!           do k=1,levs
!             do i=1,len_trie_ls
!               trie_ls(i,1,p_zq) = trie_ls(i,1,p_zq)
!    &                            + typdel(k)*trie_ls_rqt(i,1,k)
!               trie_ls(i,2,p_zq) = trie_ls(i,2,p_zq)
!    &                            + typdel(k)*trie_ls_rqt(i,2,k)
!             enddo
!             do i=1,len_trio_ls
!               trio_ls(i,1,p_zq) = trio_ls(i,1,p_zq)
!    &                            + typdel(k)*trio_ls_rqt(i,1,k)
!               trio_ls(i,2,p_zq) = trio_ls(i,2,p_zq)
!    &                            + typdel(k)*trio_ls_rqt(i,2,k)
!             enddo
!           enddo
!!

!           do k=1,levs
!             item = p_x +k-1
!             jtem = p_di+k-1
!             ktem = p_y +k-1
!             ltem = p_te+k-1
!             do i=1,len_trie_ls
!              trie_ls(i,1,item) =  trie_ls(i,1,jtem) +
!    &                 bfilte(i) * (trie_ls(i,1,item)-trie_ls(i,1,jtem))
!              trie_ls(i,2,item) =  trie_ls(i,2,jtem) +
!    &                 bfilte(i) * (trie_ls(i,2,item)-trie_ls(i,2,jtem))
!              trie_ls(i,1,ktem) =  trie_ls(i,1,ltem) +
!    &                 bfilte(i) * (trie_ls(i,1,ktem)-trie_ls(i,1,ltem))
!              trie_ls(i,2,ktem) =  trie_ls(i,2,ltem) +
!    &                 bfilte(i) * (trie_ls(i,2,ktem)-trie_ls(i,2,ltem))
!        trie_ls(i,1,p_dpn+k-1)=trie_ls(i,1,p_dp+k-1)+bfilte(i)*
!    &  (trie_ls(i,1,p_dpn+k-1)-trie_ls(i,1,p_dp+k-1))
!        trie_ls(i,2,p_dpn+k-1)=trie_ls(i,2,p_dp+k-1)+bfilte(i)*
!    &  (trie_ls(i,2,p_dpn+k-1)-trie_ls(i,2,p_dp+k-1))
!             enddo
!             do i=1,len_trio_ls
!              trio_ls(i,1,item) =  trio_ls(i,1,jtem) +
!    &                 bfilto(i) * (trio_ls(i,1,item)-trio_ls(i,1,jtem))
!              trio_ls(i,2,item) =  trio_ls(i,2,jtem) +
!    &                 bfilto(i) * (trio_ls(i,2,item)-trio_ls(i,2,jtem))
!              trio_ls(i,1,ktem) =  trio_ls(i,1,ltem) +
!    &                 bfilto(i) * (trio_ls(i,1,ktem)-trio_ls(i,1,ltem))
!              trio_ls(i,2,ktem) =  trio_ls(i,2,ltem) +
!    &                 bfilto(i) * (trio_ls(i,2,ktem)-trio_ls(i,2,ltem))
!        trio_ls(i,1,p_dpn+k-1)=trio_ls(i,1,p_dp+k-1)+bfilto(i)*
!    &  (trio_ls(i,1,p_dpn+k-1)-trio_ls(i,1,p_dp+k-1))
!        trio_ls(i,2,p_dpn+k-1)=trio_ls(i,2,p_dp+k-1)+bfilto(i)*
!    &  (trio_ls(i,2,p_dpn+k-1)-trio_ls(i,2,p_dp+k-1))
!             enddo
!           enddo

          endif    ! fin no adiabatic
!******************************************************************ADIABATIC


!$omp parallel do private(k,kp_d,kp_z,kp_t,kp_x,kp_w,kp_y,j)
           do k=1,levs
             kp_d = P_DI + K - 1
             kp_z = P_ZE + K - 1
             kp_t = P_TE + K - 1
             kp_x = P_X  + K - 1
             kp_w = P_W  + K - 1
             kp_y = P_Y  + K - 1
             do j=1,len_trie_ls
               trie_ls(j,1,kp_d) = trie_ls(j,1,kp_x)
               trie_ls(j,2,kp_d) = trie_ls(j,2,kp_x)
               trie_ls(j,1,kp_z) = trie_ls(j,1,kp_w)
               trie_ls(j,2,kp_z) = trie_ls(j,2,kp_w)
               trie_ls(j,1,kp_t) = trie_ls(j,1,kp_y)
               trie_ls(j,2,kp_t) = trie_ls(j,2,kp_y)
             enddo
             do j=1,len_trio_ls
               trio_ls(j,1,kp_d) = trio_ls(j,1,kp_x)
               trio_ls(j,2,kp_d) = trio_ls(j,2,kp_x)
               trio_ls(j,1,kp_z) = trio_ls(j,1,kp_w)
               trio_ls(j,2,kp_z) = trio_ls(j,2,kp_w)
               trio_ls(j,1,kp_t) = trio_ls(j,1,kp_y)
               trio_ls(j,2,kp_t) = trio_ls(j,2,kp_y)
             enddo
           enddo
           if(.not. gg_tracers)then
!$omp parallel do private(k,kp_x,kp_t,j)
             do k=1,levh
               kp_x = P_RQ + K - 1
               kp_t = P_RT + K - 1
               do j=1,len_trie_ls
                 trie_ls(j,1,kp_x) = trie_ls(j,1,kp_t)
                 trie_ls(j,2,kp_x) = trie_ls(j,2,kp_t)
               enddo
               do j=1,len_trio_ls
                 trio_ls(j,1,kp_x) = trio_ls(j,1,kp_t)
                 trio_ls(j,2,kp_x) = trio_ls(j,2,kp_t)
               enddo
             enddo
           endif !  if(.not.gg_tracers)
! ----------------------------------------
!         endif   ! fin no impadj
!*********************************************************************LADJ
! ----------------------------------------

        endif   ! end not start_step
! ------------------------------------------------
!      endif     ! only for fcst node

!      if (me == 0) write(0,*)' after start step if'
!**********************************************************ME<NUM_PES_FCST*
!--------------------------------------------
! =====================================================================
!--------------------------------------------
       IF(ldfi_spect) THEN
         call do_dynamics_spectdfini_slg
     &                      (kdt-kdtdfi,ndfih,trie_ls,trio_ls)
         if( kdt-kdtdfi == ndfih ) reset_step=.true.
         if( me == 0 ) print *,' do spectdfini at kdt=',kdt,
     &    'reset_step=',reset_step,kdt,kdtdfi,ndfih
       END IF
!
! =====================================================================
       IF(.not. restart_step) THEN
!
!--------------------------------------------
         IF (lsout .and. kdt /= 0 .and. .not. nemsio_out) THEN
!--------------------------------------------
!
           CALL WRTOUT_dynamics(PHOUR,FHOUR,ZHOUR,IDATE,
     &                 TRIE_LS,TRIO_LS,grid_gr,
     &                 SL,SI,
     &                 ls_node,LS_NODES,MAX_LS_NODES,
     &                 lats_nodes_a,global_lats_a,lonsperlat,
     &                 snnp1ev,snnp1od,
     &                 epsedn,epsodn,plnev_a,plnod_a,
     &                 epse  ,epso  ,plnew_a,plnow_a,
     &                 pdryini,'SIG.F')
!
! ----------------------------------
         ENDIF ! if ls_out
! ----------------------------------
!
         IF (mod(kdt,nsres) == 0 .and. kdt /= 0) THEN
!!
           CALL wrt_restart_dynamics(TRIE_LS,TRIO_LS,grid_gr,
     &          SI,fhour,idate,
     &          igen,pdryini,
     &          ls_node,ls_nodes,max_ls_nodes,
     &          global_lats_a,lonsperlat,lats_nodes_a,
     &          epse,epso,plnew_a,plnow_a,
     &          ens_nam,kdt,nfcstdate7)

         ENDIF
!
!-- end of restart step
       ELSE
          restart_step = .false.
       ENDIF

! =====================================================================
!
!       if (me == 0) write(0,*)' in dyn end_step=',end_step,
!    &' dfiend_step=',dfiend_step,' reset_step=',reset_step
!    &,' ldfi_spect=',ldfi_spect
          if(end_step .and. .not. nemsio_out) then
        RETURN

      elseif(dfiend_step) then !  dfi end step, return to dfi routine
!       if (me == 0) write(0,*)' returning from do tstep at kdt=',kdt,
!    &' reset_step=',reset_step,' dfiend_step=',dfiend_step
        return
!
      end if
!
      kdt = kdt + 1

!*********************************************** following is the tdostep.
      SHOUR = SHOUR + deltim

!     if (me == 0) write(0,*)' in slg dynamics shour=',shour,' kdt=',kdt
!      if (me < num_pes_fcst) then
!
!Start do_tstep slg part.
!------------------------

      if(ifirst == 1) then
!       allocate ( factor_b2t_ref(levs), gzie_ln(len_trie_ls,2),
!    &             gzio_ln(len_trio_ls,2) )

        allocate ( gzie_ln(len_trie_ls,2), gzio_ln(len_trio_ls,2) )

        call get_cd_hyb_slg(deltim,batah)

        k = 0
        CALL deldifs(
     &            TRIE_LS(1,1,P_RT+k-1), TRIE_LS(1,1,P_W+k-1),
     &            TRIE_LS(1,1,P_QM    ), TRIE_LS(1,1,P_X+k-1),
     &            TRIE_LS(1,1,P_Y +k-1), TRIE_LS(1,1,P_TEM+k-1),
     &            TRIO_LS(1,1,P_RT+k-1), TRIO_LS(1,1,P_W+k-1),
     &            TRIO_LS(1,1,P_QM    ), TRIO_LS(1,1,P_X+k-1),
     &            TRIO_LS(1,1,P_Y +k-1), TRIO_LS(1,1,P_TEM+k-1),
     &            deltim,SL,LS_NODE,coef00,k,hybrid,
     &            gen_coord_hybrid)

        ifirst=0
      endif

!     global_times_a = 0.
!     timer1         = rtc()

      call gloopa_hyb_slg (deltim,trie_ls,trio_ls,gzie_ln,gzio_ln,
     &                     ls_node,ls_nodes,max_ls_nodes,
     &                     lats_nodes_a,global_lats_a,
     &                     lonsperlat,ldfi,grid_gr_dfi,
     &                     epse,epso,epsedn,epsodn,
     &                     snnp1ev,snnp1od,
!    &                     snnp1ev,snnp1od,ndexev,ndexod,
     &                     plnev_a,plnod_a,pddev_a,pddod_a,
     &                     plnew_a,plnow_a,
     &                     kdt,batah,lsout, end_step)
!                  &       global_times_a,kdt,batah,lsout, end_step)

!     timer2 = rtc()

!$omp parallel do private(locl)
      do locl=1,ls_max_node
        call sicdife_hyb_slg(trie_ls(1,1,p_x  ), trie_ls(1,1,p_y ),
     &                       trie_ls(1,1,p_zq ), deltim/2.,
     &                       trie_ls(1,1,p_uln), trie_ls(1,1,p_vln),
     &                       ls_node,snnp1ev,ndexev,locl,batah)
        call sicdifo_hyb_slg(trio_ls(1,1,p_x  ), trio_ls(1,1,p_y ),
     &                       trio_ls(1,1,p_zq ), deltim/2.,
     &                       trio_ls(1,1,p_uln), trio_ls(1,1,p_vln),
     &                       ls_node,snnp1od,ndexod,locl,batah)
      enddo
      do j=1,len_trie_ls
        trie_ls(j,1,p_zq ) = trie_ls(j,1,p_zq) - gzie_ln(j,1)
        trie_ls(j,2,p_zq ) = trie_ls(j,2,p_zq) - gzie_ln(j,2)
      enddo
      do j=1,len_trio_ls
        trio_ls(j,1,p_zq ) = trio_ls(j,1,p_zq) - gzio_ln(j,1)
        trio_ls(j,2,p_zq ) = trio_ls(j,2,p_zq) - gzio_ln(j,2)
      enddo
!save n-1 values for diffusion, not really part of samilag scheme
      do j=1,len_trie_ls
        trie_ls(j,1,p_qm ) = trie_ls(j,1,p_zq)
        trie_ls(j,2,p_qm ) = trie_ls(j,2,p_zq)
      enddo
      do j=1,len_trio_ls
        trio_ls(j,1,p_qm ) = trio_ls(j,1,p_zq)
        trio_ls(j,2,p_qm ) = trio_ls(j,2,p_zq)
      enddo

!$omp parallel do private(k,kp_t,kp_y,j)
      do k=1,levs
        kp_t = p_tem + k - 1
        kp_y = p_y   + k - 1
        do j=1,len_trie_ls
          trie_ls(j,1,kp_t) = trie_ls(j,1,kp_y)
          trie_ls(j,2,kp_t) = trie_ls(j,2,kp_y)
        enddo
        do j=1,len_trio_ls
          trio_ls(j,1,kp_t) = trio_ls(j,1,kp_y)
          trio_ls(j,2,kp_t) = trio_ls(j,2,kp_y)
        enddo
      enddo

!--------------------------------------------------------
      if ( me  == me_l_0 ) then
        coef00(:,:) = 0.0
        do locl=1,ls_max_node
          l      = ls_node(locl,1)
          jbasev = ls_node(locl,2)
          if (l == 0) then
            n = 0
! 1 corresponds to temperature,  2 corresponds to ozon, 3 to clwater
            do k=1,levs
              coef00(k,1) = trie_ls(indlsev(n,l),1,p_y +k-1)
            enddo
          endif
        enddo
        coef00m = coef00
      endif

      call mpi_bcast(coef00m,levs*ntrac,mpi_r_mpi,me_l_0,mc_comp,
     &                                                      ierr)
      coef00 = coef00m
      call updown(sl,coef00(1,1))
!
      if (gg_tracers) then
!$omp parallel do shared(trie_ls,trio_ls)
!$omp+shared(deltim,sl,ls_node,coef00,hybrid,gen_coord_hybrid)
!$omp+private(k)
        do k=1,levs
          call deldifs_tracers(
     &               trie_ls(1,1,p_rt+k-1), trie_ls(1,1,p_w+k-1),
     &               trie_ls(1,1,p_qm    ), trie_ls(1,1,p_x+k-1),
     &               trie_ls(1,1,p_y +k-1), trie_ls(1,1,p_tem+k-1),
     &               trio_ls(1,1,p_rt+k-1), trio_ls(1,1,p_w+k-1),
     &               trio_ls(1,1,p_qm    ), trio_ls(1,1,p_x+k-1),
     &               trio_ls(1,1,p_y +k-1), trio_ls(1,1,p_tem+k-1),
     &               deltim,sl,ls_node,coef00,k,hybrid,
     &               gen_coord_hybrid)
        enddo
      else
!
!$omp parallel do shared(trie_ls,trio_ls)
!$omp+shared(deltim,sl,ls_node,coef00,hybrid,gen_coord_hybrid)
!$omp+private(k)
        do k=1,levs
          call deldifs(
     &               trie_ls(1,1,p_rt+k-1), trie_ls(1,1,p_w+k-1),
     &               trie_ls(1,1,p_qm    ), trie_ls(1,1,p_x+k-1),
     &               trie_ls(1,1,p_y +k-1), trie_ls(1,1,p_tem+k-1),
     &               trio_ls(1,1,p_rt+k-1), trio_ls(1,1,p_w+k-1),
     &               trio_ls(1,1,p_qm    ), trio_ls(1,1,p_x+k-1),
     &               trio_ls(1,1,p_y +k-1), trio_ls(1,1,p_tem+k-1),
     &               deltim,sl,ls_node,coef00,k,hybrid,
     &               gen_coord_hybrid)
        enddo
      endif
!--------------------------------------------------------
      do j=1,len_trie_ls
        trie_ls(j,1,p_q ) = trie_ls(j,1,p_zq)
        trie_ls(j,2,p_q ) = trie_ls(j,2,p_zq)
      enddo
      do j=1,len_trio_ls
        trio_ls(j,1,p_q ) = trio_ls(j,1,p_zq)
        trio_ls(j,2,p_q ) = trio_ls(j,2,p_zq)
      enddo
!         if (iprint .eq. 1) print*,' me = beg gloopb ',me
!         timer1 = rtc()

!!
      do i=1,len_trie_ls
        save_qe_ls(i,1) = trie_ls(i,1,p_q)
        save_qe_ls(i,2) = trie_ls(i,2,p_q)
      enddo
      do i=1,len_trio_ls
        save_qo_ls(i,1) = trio_ls(i,1,p_q)
        save_qo_ls(i,2) = trio_ls(i,2,p_q)
      enddo
!
!--------------------------------------------
! do transform from new spectral to grid  before for exit

          call spect_to_grid_slg(trie_ls,trio_ls,
     &                           syn_gr_a_1,syn_gr_a_2,
     &                           ls_node,ls_nodes,max_ls_nodes,
     &                           lats_nodes_a,global_lats_a,lonsperlat,
     &                           epsedn,epsodn,
     &                           snnp1ev,snnp1od,plnev_a,plnod_a)

! -------------------------------------------------------------------
!  get dpdt in grid point values for export state
!
        call do_dynamics_gridomega_slg(syn_gr_a_2,grid_gr,rcs2_a,
     &                                 global_lats_a,lonsperlat)

        call do_dynamics_syn2gridn_slg(syn_gr_a_2,grid_gr,
     &                                 ak5,bk5,kdt,
!    &                                 lon_dim_a,lots_sl,lats_dim_a,
     &                                 global_lats_a,lonsperlat)

!
        stp = 1
        call do_dynamics_gridpdp(grid_gr,global_lats_a,lonsperlat,stp)

!---------------------------------------------------------------------
! move data from n+1 to physics
        call do_dynamics_gridn2p(grid_gr,global_lats_a,lonsperlat)

!
! ---------------------------------------------------------------------
!
!################################################################################
!       spdmax_nodem = spdmax_node
!       call mpi_gather(spdmax_nodem,levs,MPI_R_MPI,
!    &                  spdmax_nodesm,levs,MPI_R_MPI,
!    &                  0,MC_COMP,ierr)
!!      spdmax_nodes = spdmax_nodesm
!
!sela call mpi_barrier (mpi_comm_world,ierr)
!
!--------------------------------------------
!       if ( me == 0 ) then
!--------------------------------------------
!
!        spdmax_nodes = spdmax_nodesm
!        do k=1,levs
!          spdmax(k) = cons0     !constant
!          do node=1,nodes
!            spdmax(k) = max(spdmax(k),spdmax_nodes(k,node))
!          enddo
!          spdmax(k) = sqrt(spdmax(k))
!        enddo
!
!        print*,'in do_dynamics_two_loop for spdmx at kdt=',kdt
!        print 100,(spdmax(k),k=1,levs)
!100      format(' spdmx(001:010)=',10f5.0,:/' spdmx(011:020)=',10f5.0,
!    &        :/' spdmx(021:030)=',10f5.0,:/' spdmx(031:040)=',10f5.0,
!    &        :/' spdmx(041:050)=',10f5.0,:/' spdmx(051:060)=',10f5.0,
!    &        :/' spdmx(061:070)=',10f5.0,:/' spdmx(071:080)=',10f5.0,
!    &        :/' spdmx(081:090)=',10f5.0,:/' spdmx(091:100)=',10f5.0,
!    &        :/' spdmx(101:110)=',10f5.0,:/' spdmx(111:120)=',10f5.0,
!    &        :/' spdmx(121:130)=',10f5.0,:/' spdmx(131:140)=',10f5.0,
!    &        :/' spdmx(141:150)=',10f5.0,:/' spdmx(151:160)=',10f5.0)
!
!--------------------------------------------
!       endif

!--------------------------------------------
!      call mpi_bcast(spdmax,levs,mpi_real8,0,MC_COMP,ierr)
!
!       deallocate ( spdlat )
!
!--------------------------------------------
!     endif ! only for fcst nodes
!--------------------------------------------
!################################################################################
      END SUBROUTINE do_dynamics_slg_loop
