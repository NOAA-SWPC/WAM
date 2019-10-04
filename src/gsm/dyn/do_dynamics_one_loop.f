      SUBROUTINE do_dynamics_one_loop(deltim,kdt,PHOUR,
     &                 TRIE_LS,TRIO_LS,GRID_GR, GRID_GR6,
     &                 grid_gr_dfi,   
     &                 LS_NODE,LS_NODES,MAX_LS_NODES,
     &                 LATS_NODES_A,GLOBAL_LATS_A,
     &                 LONSPERLAT,
     &                 LATS_NODES_EXT,GLOBAL_LATS_EXT,
     &                 EPSE,EPSO,EPSEDN,EPSODN,
     &                 SNNP1EV,SNNP1OD,NDEXEV,NDEXOD,
     &                 PLNEV_A,PLNOD_A,PDDEV_A,PDDOD_A,
     &                 PLNEW_A,PLNOW_A,
     &                 SYN_LS_A,DYN_LS_A,
     &                 SYN_GR_A_1,PYN_GR_A_1,DYN_GR_A_1,ANL_GR_A_1,
     &                 SYN_GR_A_2,PYN_GR_A_2,DYN_GR_A_2,ANL_GR_A_2,
     &                 pwat,ptot,ptrc,
     &                 pdryini,nblck,ZHOUR,N1,N4,
     &                 LSOUT,ldfi,COLAT1,CFHOUR1,
     &                 start_step,restart_step,reset_step,end_step,
     &                 dfiend_step,nfcstdate7,
     &                 Cpl_flag)
cc

! March 2009, Weiyu Yang modified for GEFS run.
! Aug 2010    Sarah Lu modified to compute tracer global sum
! Feb 2011    S Moorthi rearranged glbsu
! February 2010: Hann-Ming Henry Juang add non-iteration dimensional-split
!                Semi-Lagrangain (NDSL).
! Apr 2012    Henry Juang, add idea
! Apr 2012    Jun Wang, add gridp2n at reset step
! Jun 2013    Henry Juang, add option with grid moisture where requires spectral
!                          for NDSL
! Oct 2013    Henry Juang, modify wrtout passing arguments and mass_dp
!----------------------------------------------

      use gfs_dyn_resol_def
      use gfs_dyn_layout1
      use gfs_dyn_gg_def
      use gfs_dyn_vert_def
      use gfs_dyn_coordinate_def
      use gfs_dyn_date_def
      use namelist_dynamics_def
      use gfs_dyn_mpi_def
      use gfs_dyn_dfi_mod, only : gfs_dfi_grid_gr
      use gfs_dyn_tracer_config, only: glbsum                   !glbsum

      use gfs_dyn_physcons, only: p0 => con_p0
      use do_dynamics_mod

      IMPLICIT NONE
!!     
      CHARACTER(16)                     :: CFHOUR1
      INTEGER,INTENT(IN):: LONSPERLAT(LATG),N1,N4,nfcstdate7(7)
      REAL(KIND=KIND_EVOD),INTENT(IN):: deltim,PHOUR
      REAL(KIND=KIND_EVOD),INTENT(INOUT):: ZHOUR
!
      type(gfs_dfi_grid_gr),intent(inout) :: grid_gr_dfi
      logical,intent(in)  :: ldfi
!!     
      INTEGER NBLCK
!!!   LOTALL=13*LEVS+3*LEVH+8
      REAL(KIND=KIND_EVOD) TRIE_LS(LEN_TRIE_LS,2,LOTls)
      REAL(KIND=KIND_EVOD) TRIO_LS(LEN_TRIO_LS,2,LOTls)
      REAL(KIND=KIND_EVOD) GRID_GR (lonf*lats_node_a_max,lotgr    )
      REAL(KIND=KIND_EVOD) GRID_GR6(lonf*lats_node_a_max,lotgr6 * 2 )
!
      integer          ls_node(ls_dim,3)
!
      INTEGER          LS_NODES(LS_DIM,NODES)
      INTEGER          MAX_LS_NODES   (NODES)
      INTEGER          LATS_NODES_A   (NODES)
      INTEGER          LATS_NODES_EXT (NODES)
      INTEGER          GLOBAL_LATS_A(LATG)
      INTEGER          GLOBAL_LATS_EXT(LATG+2*JINTMX+2*NYPT*(NODES-1))
!
      real(kind=kind_evod) colat1

      REAL(KIND=KIND_EVOD)      EPSE(LEN_TRIE_LS)
      REAL(KIND=KIND_EVOD)      EPSO(LEN_TRIO_LS)
      REAL(KIND=KIND_EVOD)    EPSEDN(LEN_TRIE_LS)
      REAL(KIND=KIND_EVOD)    EPSODN(LEN_TRIO_LS)
      REAL(KIND=KIND_EVOD)   SNNP1EV(LEN_TRIE_LS)
      REAL(KIND=KIND_EVOD)   SNNP1OD(LEN_TRIO_LS)
      INTEGER                 NDEXEV(LEN_TRIE_LS)
      INTEGER                 NDEXOD(LEN_TRIO_LS)

      REAL(KIND=KIND_EVOD)   PLNEV_A(LEN_TRIE_LS,LATG2)
      REAL(KIND=KIND_EVOD)   PLNOD_A(LEN_TRIO_LS,LATG2)
      REAL(KIND=KIND_EVOD)   PDDEV_A(LEN_TRIE_LS,LATG2)
      REAL(KIND=KIND_EVOD)   PDDOD_A(LEN_TRIO_LS,LATG2)
      REAL(KIND=KIND_EVOD)   PLNEW_A(LEN_TRIE_LS,LATG2)
      REAL(KIND=KIND_EVOD)   PLNOW_A(LEN_TRIO_LS,LATG2)

      REAL(KIND=KIND_EVOD) SYN_LS_A(4*LS_DIM,LOTS,LATG2)
      REAL(KIND=KIND_EVOD) DYN_LS_A(4*LS_DIM,LOTD,LATG2)

      REAL(KIND=KIND_EVOD) SYN_GR_A_1(LONFX*LOTS,LATS_DIM_EXT)
      REAL(KIND=KIND_EVOD) PYN_GR_A_1(LONFX*LOTP,LATS_DIM_EXT)
      REAL(KIND=KIND_EVOD) DYN_GR_A_1(LONFX*LOTD,LATS_DIM_EXT)
      REAL(KIND=KIND_EVOD) ANL_GR_A_1(LONFX*LOTA,LATS_DIM_EXT)

      REAL(KIND=KIND_EVOD) SYN_GR_syq(LONFX*levh,LATS_DIM_EXT)

      REAL(KIND=KIND_EVOD) SYN_GR_A_2(LONFX*LOTS,LATS_DIM_EXT)
      REAL(KIND=KIND_EVOD) PYN_GR_A_2(LONFX*LOTP,LATS_DIM_EXT)
      REAL(KIND=KIND_EVOD) DYN_GR_A_2(LONFX*LOTD,LATS_DIM_EXT)
      REAL(KIND=KIND_EVOD) ANL_GR_A_2(LONFX*LOTA,LATS_DIM_EXT)

      real(kind=kind_evod) pdot(lonf,levs+1,lats_dim_ext)
!
      INTEGER LEV,LEVMAX
      real (kind=kind_grid) pdryini,pcorr
      real (kind=kind_grid) ptot(lonf,lats_node_a)
      real (kind=kind_grid) pwat(lonf,lats_node_a)
      real (kind=kind_grid) ptotj(lats_node_a),pwatj(lats_node_a)
      real (kind=kind_grid) ptotg(latg),pwatg(latg)
      real (kind=kind_grid) sumwa,sumto,ptotp,pwatp,pdryg

! For tracer gloabl sum (Sarah Lu)
      real (kind=kind_grid) ptrc(lonf,lats_node_a,ntrac)               !glbsum
      real (kind=kind_grid) ptrcj(lats_node_a,ntrac),tmpj(lats_node_a) !glbsum
      real (kind=kind_grid) ptrcg(latg),sumtrc(ntrac),ptrcp(ntrac)     !glbsum
!
      REAL (KIND=KIND_grid) filtb
      REAL (KIND=KIND_grid) cons0,cons1,cons2,cons0p5
      INTEGER               kdt,IERR,I,J,K,L,LOCL,N
      real(kind=kind_evod)  ye1(levs)
      REAL(KIND=kind_mpi)   coef00m(LEVS,ntrac) ! temp. ozone clwater  
      REAL(KIND=kind_evod)  coef00(LEVS,ntrac) ! temp. ozone clwater  
      INTEGER               INDLSEV,JBASEV
      INTEGER               INDLSOD,JBASOD

c idea-related changes
c Introduced arrays cvd00 and cvd00m (global mean temperature,h2o,o3,cld,
c  O and O2) similar to coef00 and coef00m and added diffusion coeffs.
c idea add 1
      REAL(KIND=kind_mpi)  cvd00m(LEVS,0:ntrac) ! temp, h2o,o3,water,o,o2
      REAL(KIND=kind_evod) cvd00(LEVS,0:ntrac) ! temp, h2o,o3,water,o,o2
      REAL(KIND=kind_evod),dimension(levs):: visc,cond,diff,tem,plyr

      integer               lan,lat
      integer               lon_dim,lons_lat,node,nvcn
      integer               iblk,njeff,lon,stp
      integer , parameter :: ngptcd = 12
      include 'function2'
      LOGICAL               LSOUT,ex_out
      LOGICAL               start_step,reset_step,end_step,dfiend_step
      LOGICAL               restart_step

      LOGICAL,      INTENT(inout) :: Cpl_flag
!!     
      LOGICAL, save               :: fwd_step = .true.
      LOGICAL, save               :: firstdfi = .true.
      REAL (KIND=KIND_grid), save :: dt,dt2,rdt2
!
      real(kind=kind_evod) xvcn
!
!     saved vertical advection of tracers from time step n-1
      real(kind=kind_evod),allocatable, save:: szdrdt(:,:)
      logical,save:: zfirst=.true.
!
      integer              iter_max
      integer              rc1
!
      real(kind=kind_evod) , allocatable :: spdlat(:,:)
!
      real(kind=kind_evod) spdmax_node (levs)
      real(kind=kind_mpi) spdmax_nodem (levs)
      real(kind=kind_evod) spdmax_nodes(levs,nodes)
      real(kind=kind_mpi) spdmax_nodesm(levs,nodes)
      real(kind=kind_mpi) thref_mpi (levs)
!
!
! timings
      real(kind=kind_evod) global_times_a(latg,nodes)
      integer tag,ireq1,ireq2
      real*8 rtc ,timer1,timer2
!
!
!     print *,' ----------- run in one loop ------------ '
!
!     shour=shour+deltim
      shour= kdt * deltim
      fhour = shour/3600.
      cons0=0.0d0
      cons1=1.0d0
      cons2=2.0d0
      cons0p5 = 0.5d0                        !constant
      filtb = (cons1-filta)*cons0p5          !constant

      if (me<num_pes_fcst) then                                         
      if( zfirst .and. .not.ndslfv ) then
        allocate (szdrdt(lonfx*levh,lats_dim_a))
        szdrdt=0.
      endif
      end if

! If input from the ensemble coupler ESMF state, skip the first linear
! computation part.
!---------------------------------------------------------------------
      IF(.NOT. CPl_flag) THEN
!
!----------------------------------------------------------
      if (me<num_pes_fcst) then                                          
!----------------------------------------------------------
!
!!
! -----------------------------------
      if( start_step ) then
! --------------------------------------------------------------
! if the first step from internal state, prepare syn for nonlinear 
! ----- this section is called once only -------
! --------------------------------------------------------------
!      print *,' start step from internal state '

        fwd_step = .true.
        dt   = deltim*0.5
        dt2  = cons2*dt
        rdt2 = 1./dt2

! prepare temp profile for semi_implicit matrix if need
        if ( semi_implicit_temp_profile ) then
          if (me.eq.me_l_0) then
            do k=1,levs
              n=p_te+k-1
              thref_mpi(k)=trie_ls(1,1,n)/sqrt(2.)
            enddo
!           print *,' thref from te is ',(thref_mpi(k),k=1,levs)
          endif
          call mpi_bcast(thref_mpi,levs,MPI_R_MPI,me_l_0,MC_COMP,ierr)
        endif

        if( gen_coord_hybrid ) then
          if( mass_dp ) then
            call get_cd_hyb_gcdp(dt)
          else
            call get_cd_hyb_gc(dt)
          endif
        else if( hybrid ) then
          call get_cd_hyb(dt)
        else
          call get_cd_sig(am,bm,dt,tov,sv)
        endif

        call spect_to_grid
     &      (trie_ls,trio_ls, 
     &       syn_gr_a_1,syn_gr_a_2,
     &       ls_node,ls_nodes,max_ls_nodes,
     &       lats_nodes_a,global_lats_a,lonsperlat,
     &       epse,epso,epsedn,epsodn,
     &       snnp1ev,snnp1od,plnev_a,plnod_a)
!      print *,' finish spect_to_grid and reset start_step to be false'
        start_step = .false.
! -------------------------------------------------------
      else if( reset_step ) then
! --------------------------------------------------------------
! if it is reset step to reset internal state to be import state
! ----- this section is called once only -------
! --------------------------------------------------------------

        fwd_step = .true.
        dt   = deltim*0.5
        dt2  = cons2*dt
        rdt2 = 1./dt2

        if( gen_coord_hybrid ) then
          if( mass_dp ) then
            call get_cd_hyb_gcdp(dt)
          else
            call get_cd_hyb_gc(dt)
          endif
        else if( hybrid ) then
          call get_cd_hyb(dt)
        else
          call get_cd_sig(am,bm,dt,tov,sv)
        endif
!
! move data from physics to n+1
        call do_dynamics_gridp2n(grid_gr,global_lats_a,lonsperlat)

! after digital filter, all updated in n

        call do_dynamics_gridn2anl(grid_gr,anl_gr_a_2,
     &                             global_lats_a,lonsperlat)

        call grid_to_spect(anl_gr_a_1,anl_gr_a_2,
     &       trie_ls,trio_ls,
     &       ls_node,ls_nodes,max_ls_nodes,
     &       lats_nodes_a,global_lats_a,lonsperlat,
     &       epse,epso,plnew_a,plnow_a)
!
!set n and n-1 time level values as import
! spectral
!       print *,' set time level n to time level n-1 '
        call do_dynamics_spectn2c(trie_ls,trio_ls)
        call do_dynamics_spectn2m(trie_ls,trio_ls)
! grid
        call do_dynamics_gridn2c(grid_gr,
     &                           global_lats_a,lonsperlat)
        call do_dynamics_gridn2m(grid_gr,
     &                           global_lats_a,lonsperlat)
!
! after digital filter, prepare syn for nonlinear forcing
!
        call spect_to_grid
     &      (trie_ls,trio_ls, 
     &       syn_gr_a_1,syn_gr_a_2,
     &       ls_node,ls_nodes,max_ls_nodes,
     &       lats_nodes_a,global_lats_a,lonsperlat,
     &       epse,epso,epsedn,epsodn,
     &       snnp1ev,snnp1od,plnev_a,plnod_a)
!
        reset_step = .false.
!
! -------------------------------------------------------
!### for restart step
      elseif(restart_step) then
!
        fwd_step = .false.
        dt   = deltim
        dt2  = cons2*dt
        rdt2 = 1./dt2
!
        if( gen_coord_hybrid ) then
          if( mass_dp ) then
            call get_cd_hyb_gcdp(dt)
          else
            call get_cd_hyb_gc(dt)
          endif
        else if( hybrid ) then
          call get_cd_hyb(dt)
        else
          call get_cd_sig(am,bm,dt,tov,sv)
        endif
!
        zfirst=.false.

        call spect_to_grid
     &      (trie_ls,trio_ls,
     &       syn_gr_a_1,syn_gr_a_2,
     &       ls_node,ls_nodes,max_ls_nodes,
     &       lats_nodes_a,global_lats_a,lonsperlat,
     &       epse,epso,epsedn,epsodn,
     &       snnp1ev,snnp1od,plnev_a,plnod_a)

! -------------------------------------------------------
      else	! end start_step, begin not start_step 
! ------------------------------------------------------
! start linear computation for model dynamics
! do after physics, not from input 
! -----------------------------------------------------------
! compute total tendency (linear and nonlinea, after physics if any)

! move data from physics to n+1
        if( process_split ) then
          call do_dynamics_gridap2n(grid_gr,global_lats_a,lonsperlat)
        else
          call do_dynamics_gridp2n(grid_gr,global_lats_a,lonsperlat)
        endif

        IF(kdt == 1) THEN
            grid_gr6 = grid_gr(:, 2 : lotgr6 * 2 + 1)
            call model_to_common_vars
     &          (grid_gr6(1,g_qm  - 1),
     &           grid_gr6(1,g_ttm - 1),
     &           grid_gr6(1,g_rm  - 1),
     &           grid_gr6(1,g_uum - 1),
     &           grid_gr6(1,g_vvm - 1),
     &           grid_gr6(1,g_dpm - 1),
     &           grid_gr (1,g_p      ),
     &           grid_gr (1,g_dpdt   ),
     &           global_lats_a, lonsperlat, 0)

            call model_to_common_vars
     &          (grid_gr6(1,g_q  - 1),
     &           grid_gr6(1,g_tt - 1),
     &           grid_gr6(1,g_rq - 1),
     &           grid_gr6(1,g_uu - 1),
     &           grid_gr6(1,g_vv - 1),
     &           grid_gr6(1,g_dp - 1),
     &           grid_gr (1,g_p      ),
     &           grid_gr (1,g_dpdt   ),
     &           global_lats_a, lonsperlat, 0)
        END IF

        call do_dynamics_gridt2anl(grid_gr,anl_gr_a_2,rdt2,
     &                             global_lats_a,lonsperlat)
!
!----------------------------------------------------------
! transform total tendency in grid to spectral
!
        call grid_to_spect(anl_gr_a_1,anl_gr_a_2,
     &     trie_ls,trio_ls,
     &     ls_node,ls_nodes,max_ls_nodes,
     &     lats_nodes_a,global_lats_a,lonsperlat,
     &     epse,epso,plnew_a,plnow_a)
!
!-----------------------------------------------------------------
! adjust moisture changes to the total mass to conserve dry mass
!
!     write(0,*)' pdryini=',pdryini

        do lan=1,lats_node_a
          lat = global_lats_a(ipt_lats_node_a-1+lan)
          lons_lat = lonsperlat(lat)
          ptotp=0.
          pwatp=0.
          do i=1,lons_lat
             ptotp     = ptotp + ptot(i,lan)
             pwatp     = pwatp + pwat(i,lan)
          enddo
          pwatj(lan)=pwatp/(2.*lonsperlat(lat))
          ptotj(lan)=ptotp/(2.*lonsperlat(lat))
        enddo
        call excha(lats_nodes_a,global_lats_a,ptotj,pwatj,ptotg,pwatg)
        sumwa=0.
        sumto=0.
        do lat=1,latg
           sumto=sumto+wgt_a(min(lat,latg-lat+1))*ptotg(lat)
           sumwa=sumwa+wgt_a(min(lat,latg-lat+1))*pwatg(lat)
        enddo
!        if(kdt==5) then
!      print *,' sumto=',sumto,' sumwa=',sumwa,' ptotg=',ptotg(1:latg)
!      print *,' sumwa=',sumwa,'pwatg=',pwatg(1:latg)
!        endif

        if ( glbsum ) then                                              !glbsum
          do lan=1,lats_node_a
            lat      = global_lats_a(ipt_lats_node_a-1+lan)
            lons_lat = lonsperlat(lat)
            ptrcp(:) = 0.
            do n = 1, ntrac
              do i=1,lons_lat
                ptrcp(n)   = ptrcp(n) + ptrc(i,lan,n)
              enddo
              ptrcj(lan,n) = ptrcp(n) / (2.*lons_lat)
            enddo
          enddo
!
          do n = 1, ntrac                                               !glbsum
           sumtrc(n)=0.                                                 !glbsum
           tmpj(:) = ptrcj(:,n)                                         !glbsum
           call excha(lats_nodes_a,global_lats_a,ptotj,tmpj,ptotg,ptrcg)!glbsum
           do lat=1,latg                                                !glbsum
             sumtrc(n)=sumtrc(n)+wgt_a(min(lat,latg-lat+1))*ptrcg(lat)  !glbsum
           enddo                                                        !glbsum
          enddo                                                         !glbsum
        endif                                                           !glbsum

        pdryg=sumto-sumwa
        if(pdryini.le.0.) pdryini=pdryg
            
        if( gen_coord_hybrid ) then                               ! hmhj
          pcorr=(pdryini-pdryg)      *sqrt(2.)                    ! hmhj
        else                                                      ! hmhj
          pcorr=(pdryini-pdryg)/sumto*sqrt(2.)
        endif                                                     ! hmhj
!
        if (me.eq.me_l_0) then
          print *,' pdryini pdryg pcorr = ',pdryini,pdryg,pcorr
          trie_ls(1,1,p_zq)=trie_ls(1,1,p_zq)+pcorr/dt2
          do k=p_dpn,p_dpn+levs-1
            n=k-p_dpn+1
            trie_ls(1,1,k)=trie_ls(1,1,k)+pcorr/dt2*dbk(n)
          enddo
        endif

        if (glbsum .and. me.eq.me_l_0) then                             !glbsum
          write(70,111) kdt,fhour,idate                                 !glbsum
          write(71,*)   kdt,(sumtrc(n),n=4,ntrac)                       !glbsum
        endif                                                           !glbsum
111     format ('kdt, fhour, idate=',i6,1x,f10.3,2x,4(i4,2x))           !glbsum

!----------------------------------------------------------
! update in spectral for vorticity and tracers in explicit way       
!
        call do_dynamics_spectupdatewrt(trie_ls,trio_ls,dt2)
        do locl=1,ls_max_node
                l=ls_node(locl,1)
           jbasev=ls_node(locl,2)
           if ( l .eq. 0 ) then
              n = 0
              do k=1,levs
                 trie_ls(indlsev(n,l),1,P_w+k-1)=cons0     !constant
                 trie_ls(indlsev(n,l),2,P_w+k-1)=cons0     !constant
              enddo
           endif
        end do
!
!----------------------------------------------------------
! update in spectral for div, temp and ps in semi-implicit or explicit
!
!----------------------------------------------------------
        if( explicit ) then						
!----------------------------------------------------------
!
          call do_dynamics_spectupdatexydpnzq(trie_ls,trio_ls,dt2)
!
!----------------------------------------------------------
        else	! do semi-implicit
!----------------------------------------------------------
!
!
          if( gen_coord_hybrid ) then                                       

            if( mass_dp ) then

!$omp parallel do private(locl)
              do locl=1,ls_max_node

       call sicdife_hyb_gcdp(trie_ls(1,1,P_dim), trie_ls(1,1,P_tem),
     &                       trie_ls(1,1,P_dpm), trie_ls(1,1,P_x  ),
     &                       trie_ls(1,1,P_y  ), trie_ls(1,1,P_dpn),
     &                       trie_ls(1,1,P_di ), trie_ls(1,1,P_te ),
     &                       trie_ls(1,1,P_dp ),dt,
     &                       trie_ls(1,1,P_uln), trie_ls(1,1,P_vln),
     &                       ls_node,snnp1ev,ndexev,locl)

       call sicdifo_hyb_gcdp(trio_ls(1,1,P_dim), trio_ls(1,1,P_tem),
     &                       trio_ls(1,1,P_dpm), trio_ls(1,1,P_x  ),
     &                       trio_ls(1,1,P_y  ), trio_ls(1,1,P_dpn),
     &                       trio_ls(1,1,P_di ), trio_ls(1,1,P_te ),
     &                       trio_ls(1,1,P_dp ),dt,
     &                       trio_ls(1,1,P_uln), trio_ls(1,1,P_vln),
     &                       ls_node,snnp1od,ndexod,locl)
              enddo
! update zq
              trie_ls(:,:,p_zq) = 0.0
              trio_ls(:,:,p_zq) = 0.0
              do k=levs,1,-1
                do i=1,len_trie_ls
                   trie_ls(i,1,P_zq)=
     x             trie_ls(i,1,P_zq)+ trie_ls(i,1,P_dpn+k-1)
                   trie_ls(i,2,P_zq)=
     x             trie_ls(i,2,P_zq)+ trie_ls(i,2,P_dpn+k-1)
                enddo
                do i=1,len_trio_ls
                   trio_ls(i,1,P_zq)=
     x             trio_ls(i,1,P_zq)+ trio_ls(i,1,P_dpn+k-1)
                   trio_ls(i,2,P_zq)=
     x             trio_ls(i,2,P_zq)+ trio_ls(i,2,P_dpn+k-1)
                enddo
              enddo

            else        ! end of mass_dp

!$omp parallel do private(locl)
              do locl=1,ls_max_node                                             

         call sicdife_hyb_gc(trie_ls(1,1,P_dim), trie_ls(1,1,P_tem),    
     &                       trie_ls(1,1,P_qm ), trie_ls(1,1,P_x  ),    
     &                       trie_ls(1,1,P_y  ), trie_ls(1,1,P_zq ),    
     &                       trie_ls(1,1,P_di ), trie_ls(1,1,P_te ),    
     &                       trie_ls(1,1,P_q  ),dt,                 
     &                       trie_ls(1,1,P_uln), trie_ls(1,1,P_vln),    
     &                       ls_node,snnp1ev,ndexev,locl)               

         call sicdifo_hyb_gc(trio_ls(1,1,P_dim), trio_ls(1,1,P_tem),    
     &                       trio_ls(1,1,P_qm ), trio_ls(1,1,P_x  ),    
     &                       trio_ls(1,1,P_y  ), trio_ls(1,1,P_zq ),    
     &                       trio_ls(1,1,P_di ), trio_ls(1,1,P_te ),    
     &                       trio_ls(1,1,P_q  ),dt,                 
     &                       trio_ls(1,1,P_uln), trio_ls(1,1,P_vln),    
     &                       ls_node,snnp1od,ndexod,locl)               
              enddo                                                             

            endif       ! end of gc

          else if(hybrid)then

!$omp parallel do private(locl)
            do locl=1,ls_max_node

         call sicdife_hyb(trie_ls(1,1,P_dim), trie_ls(1,1,P_tem),
     &                    trie_ls(1,1,P_qm ), trie_ls(1,1,P_x  ),
     &                    trie_ls(1,1,P_y  ), trie_ls(1,1,P_zq ),
     &                    trie_ls(1,1,P_di ), trie_ls(1,1,P_te ),
     &                    trie_ls(1,1,P_q  ),dt,
     &                    trie_ls(1,1,P_uln), trie_ls(1,1,P_vln),
     &                    ls_node,snnp1ev,ndexev,locl)

         call sicdifo_hyb(trio_ls(1,1,P_dim), trio_ls(1,1,P_tem),
     &                    trio_ls(1,1,P_qm ), trio_ls(1,1,P_x  ),
     &                    trio_ls(1,1,P_y  ), trio_ls(1,1,P_zq ),
     &                    trio_ls(1,1,P_di ), trio_ls(1,1,P_te ),
     &                    trio_ls(1,1,P_q  ),dt,
     &                    trio_ls(1,1,P_uln), trio_ls(1,1,P_vln),
     &                    ls_node,snnp1od,ndexod,locl)
            enddo

          else

!$omp parallel do private(locl)
            do locl=1,ls_max_node

         CALL SICDIFE_sig(TRIE_LS(1,1,P_DIM), TRIE_LS(1,1,P_TEM),
     &                    TRIE_LS(1,1,P_QM ), TRIE_LS(1,1,P_X  ),
     &                    TRIE_LS(1,1,P_Y  ), TRIE_LS(1,1,P_ZQ ),
     &                    AM,BM,TOV,SV,dt,
     &                    TRIE_LS(1,1,P_ULN), TRIE_LS(1,1,P_VLN),
     &                    LS_NODE,SNNP1EV,NDEXEV,locl,TRIE_LS(1,1,P_DI))

         CALL SICDIFO_sig(TRIO_LS(1,1,P_DIM), TRIO_LS(1,1,P_TEM),
     &                    TRIO_LS(1,1,P_QM ), TRIO_LS(1,1,P_X  ),  
     &                    TRIO_LS(1,1,P_Y  ), TRIO_LS(1,1,P_ZQ ),
     &                    AM,BM,TOV,SV,dt,
     &                    TRIO_LS(1,1,P_ULN), TRIO_LS(1,1,P_VLN),
     &                    LS_NODE,SNNP1OD,NDEXOD,locl,TRIO_LS(1,1,P_DI))
            enddo

          endif ! end of all sicdif

!----------------------------------------------------------
        endif 		! not explicit					
!----------------------------------------------------------
!----------------------------------------------------------
! compute coef00 for all, even for hybrid mode
      coef00(:,:) = 0.0
! idea add one line
      if( lsidea ) cvd00(:,:) = 0.0
      IF ( ME .EQ. ME_L_0 ) THEN
        DO LOCL=1,LS_MAX_NODE
          l=ls_node(locl,1)
          jbasev=ls_node(locl,2)
          IF ( L .EQ. 0 ) THEN
            N=0
! 1 Corresponds to temperature,  2 corresponds to ozon, 3 to clwater
            DO K=1,LEVS
              coef00(K,1) = TRIE_LS(INDLSEV(N,L),1,P_Y +K-1)
              if (ntoz .gt. 1 .and.                                     
     &            .not. (hybrid.or.gen_coord_hybrid)) then              
                if( .not.ndslfv ) then
                  coef00(K,ntoz) = TRIE_LS(INDLSEV(N,L),1,
     &                                   (ntoz-1)*levs+P_rt+K-1)
                else
                  coef00(k,ntoz) = 0.0
                endif
              endif
c idea add 2
              if( lsidea ) then
              cvd00(K,0) = TRIE_LS(INDLSEV(N,L),1,P_te +K-1)
              if( .not.ndslfv ) then
               do i=1,ntrac
                cvd00(K,i) = TRIE_LS(INDLSEV(N,L),1,(i-1)*levs+P_rq+K-1)
               enddo
              else
               do i=1,ntrac
                cvd00(K,i) = grid_gr(1,(i-1)*levs+g_rq+k-1)
                cvd00(K,i) = cvd00(k,i)/(.5*sqrt(2.))
               enddo
              endif
              endif
c idea add end
            ENDDO
          ENDIF
        END DO
      END IF
      coef00m = coef00
      CALL MPI_BCAST(coef00m,levs*ntrac,MPI_R_MPI,ME_L_0,MC_COMP,IERR)
      coef00=coef00m
c idea add 3
      if( lsidea ) then
      cvd00m(:,:)=cvd00(:,:)
      CALL MPI_BCAST(cvd00m,levs*(ntrac+1),MPI_R_MPI,ME_L_0,
     &MC_COMP,IERR)
      endif
c idea add end
      if( gen_coord_hybrid ) then                                       
        call updown_gc(sl,coef00(1,1))                                  
      else                                                              
        call updown(sl,coef00(1,1))
      endif                                                             
      if (ntoz .gt. 1 .and. .not. (hybrid.or.gen_coord_hybrid)) then    
               call updown(sl,coef00(1,ntoz))
      endif
!
      if( lsidea ) then

c idea add 4
c Get cvd00 back and normalize
        cvd00(:,:)=(.5*sqrt(2.))*cvd00m(:,:)
c Calculate global mean viscosity, conductivity, and diffusion
c coefficients
        visc=0.
        cond=0.
        diff=0.
        do k=1,levs
          plyr(k)=(ak5(k)+ak5(k+1)+p0*1.e-3*(bk5(k)+bk5(k+1)))*500.
        enddo
!hmhj debug
!       do k=0,ntrac
!       call mymaxmin(cvd00(1,k),levs,levs,1,' cvd00 ')
!       enddo
!       call mymaxmin(plyr,levs,levs,1,' plyr ')

        call idea_getcoef(levs,ntrac,cvd00,plyr,visc,cond,diff)
!hmhj debug
!       call mymaxmin(visc,levs,levs,1,' visc ')
!       call mymaxmin(cond,levs,levs,1,' cond ')
!       call mymaxmin(diff,levs,levs,1,' diff ')

      endif
!
! -----------------------------------------------------------------
      if( ndslfv ) then

        if( lsidea ) then

        print *,' ndslfv and lsidea are on '
!
!$omp parallel do shared(TRIE_LS,TRIO_LS)
!$omp+shared(dt,SL,LS_NODE,coef00,hybrid,gen_coord_hybrid)
!$omp+shared(visc,cond,diff)
        do k=1,levs
         CALL idea_deldifs_noq
     &               (TRIE_LS(1,1,P_W+k-1),
     &                TRIE_LS(1,1,P_QM    ), TRIE_LS(1,1,P_X+k-1),
     &                TRIE_LS(1,1,P_Y +k-1), TRIE_LS(1,1,P_TEM+k-1),
     &                TRIO_LS(1,1,P_W+k-1),
     &                TRIO_LS(1,1,P_QM    ), TRIO_LS(1,1,P_X+k-1),
     &                TRIO_LS(1,1,P_Y +k-1), TRIO_LS(1,1,P_TEM+k-1),
     &                dt,SL,LS_NODE,coef00,k,hybrid,
     &                gen_coord_hybrid,visc,cond,diff)
        enddo
        else
        print *,' ndslfv is on and lsidea is off '
!
!$omp parallel do shared(TRIE_LS,TRIO_LS)
!$omp+shared(dt,SL,LS_NODE,coef00,hybrid,gen_coord_hybrid)
        do k=1,levs
          CALL deldifs_noq
     &               (TRIE_LS(1,1,P_W+k-1),
     &                TRIE_LS(1,1,P_QM    ), TRIE_LS(1,1,P_X+k-1),
     &                TRIE_LS(1,1,P_Y +k-1), TRIE_LS(1,1,P_TEM+k-1),
     &                TRIO_LS(1,1,P_W+k-1),
     &                TRIO_LS(1,1,P_QM    ), TRIO_LS(1,1,P_X+k-1),
     &                TRIO_LS(1,1,P_Y +k-1), TRIO_LS(1,1,P_TEM+k-1),
     &                dt,SL,LS_NODE,coef00,k,hybrid,
     &                gen_coord_hybrid)
        enddo

        endif

      else


        if( lsidea ) then
        print *,' ndslfv is off and lsidea is on '
!
!$omp parallel do shared(TRIE_LS,TRIO_LS)
!$omp+shared(dt,SL,LS_NODE,coef00,hybrid,gen_coord_hybrid)
!$omp+shared(visc,cond,diff)
        do k=1,levs
         CALL idea_deldifs(TRIE_LS(1,1,P_RT+k-1), TRIE_LS(1,1,P_W+k-1),
     &                TRIE_LS(1,1,P_QM    ), TRIE_LS(1,1,P_X+k-1),
     &                TRIE_LS(1,1,P_Y +k-1), TRIE_LS(1,1,P_TEM+k-1),
     &                TRIO_LS(1,1,P_RT+k-1), TRIO_LS(1,1,P_W+k-1),
     &                TRIO_LS(1,1,P_QM    ), TRIO_LS(1,1,P_X+k-1),
     &                TRIO_LS(1,1,P_Y +k-1), TRIO_LS(1,1,P_TEM+k-1),
     &                dt,SL,LS_NODE,coef00,k,hybrid,
     &                gen_coord_hybrid,visc,cond,diff)
        enddo

        else
        print *,' ndslfv is off and lsidea is off '


!$omp parallel do shared(TRIE_LS,TRIO_LS)
!$omp+shared(dt,SL,LS_NODE,coef00,hybrid)
        do k=1,levs
         CALL deldifs(TRIE_LS(1,1,P_RT+k-1), TRIE_LS(1,1,P_W+k-1),
     &                TRIE_LS(1,1,P_QM    ), TRIE_LS(1,1,P_X+k-1),
     &                TRIE_LS(1,1,P_Y +k-1), TRIE_LS(1,1,P_TEM+k-1),    
     &                TRIO_LS(1,1,P_RT+k-1), TRIO_LS(1,1,P_W+k-1),
     &                TRIO_LS(1,1,P_QM    ), TRIO_LS(1,1,P_X+k-1),
     &                TRIO_LS(1,1,P_Y +k-1), TRIO_LS(1,1,P_TEM+k-1),    
     &                dt,SL,LS_NODE,coef00,k,hybrid,                
     &                gen_coord_hybrid)                                 
        enddo

        endif

      endif

! ----------------------------------------------------------------------
      if( ndslfv ) then

!$omp parallel do shared(TRIE_LS,NDEXEV,TRIO_LS,NDEXOD)
!$omp+shared(SL,SPDMAX,dt,LS_NODE)
        do k=1,levs
         CALL damp_speed_noq
     &                  (TRIE_LS(1,1,P_X+k-1), TRIE_LS(1,1,P_W +k-1),
     &                   TRIE_LS(1,1,P_Y+k-1),
     &                   NDEXEV,
     &                   TRIO_LS(1,1,P_X+k-1), TRIO_LS(1,1,P_W +k-1),
     &                   TRIO_LS(1,1,P_Y+k-1),
     &                   NDEXOD,
     &                   SL,SPDMAX(k),dt,LS_NODE)
        enddo

      else

!$omp parallel do shared(TRIE_LS,NDEXEV,TRIO_LS,NDEXOD)
!$omp+shared(SL,SPDMAX,dt,LS_NODE)
        do k=1,levs
         CALL damp_speed(TRIE_LS(1,1,P_X+k-1), TRIE_LS(1,1,P_W +k-1),
     &                   TRIE_LS(1,1,P_Y+k-1), TRIE_LS(1,1,P_RT+k-1),
     &                   NDEXEV,
     &                   TRIO_LS(1,1,P_X+k-1), TRIO_LS(1,1,P_W +k-1),
     &                   TRIO_LS(1,1,P_Y+k-1), TRIO_LS(1,1,P_RT+k-1),
     &                   NDEXOD,
     &                   SL,SPDMAX(k),dt,LS_NODE)
        enddo

      endif
!
!-------------------------------------------
      if(.not.fwd_step)then
!-------------------------------------------
      if( ndslfv ) then
      CALL filtereo_noq
     &             (TRIE_LS(1,1,P_TEM), TRIE_LS(1,1,P_TE ),
     &              TRIE_LS(1,1,P_Y  ), TRIE_LS(1,1,P_DIM),
     &              TRIE_LS(1,1,P_DI ), TRIE_LS(1,1,P_X  ),
     &              TRIE_LS(1,1,P_ZEM), TRIE_LS(1,1,P_ZE ),
     &              TRIE_LS(1,1,P_W  ),
     &              TRIE_LS(1,1,P_dpm), TRIE_LS(1,1,P_dp ),
     &              TRIE_LS(1,1,P_dpn),
     &              TRIO_LS(1,1,P_TEM), TRIO_LS(1,1,P_TE ),
     &              TRIO_LS(1,1,P_Y  ), TRIO_LS(1,1,P_DIM),
     &              TRIO_LS(1,1,P_DI ), TRIO_LS(1,1,P_X  ),
     &              TRIO_LS(1,1,P_ZEM), TRIO_LS(1,1,P_ZE ),
     &              TRIO_LS(1,1,P_W  ),
     &              TRIO_LS(1,1,P_dpm), TRIO_LS(1,1,P_dp ),
     &              TRIO_LS(1,1,P_dpn),
     &              FILTA,LS_NODE)
      else
      CALL filtereo(TRIE_LS(1,1,P_TEM), TRIE_LS(1,1,P_TE ),
     &              TRIE_LS(1,1,P_Y  ), TRIE_LS(1,1,P_DIM),
     &              TRIE_LS(1,1,P_DI ), TRIE_LS(1,1,P_X  ),
     &              TRIE_LS(1,1,P_ZEM), TRIE_LS(1,1,P_ZE ),
     &              TRIE_LS(1,1,P_W  ), TRIE_LS(1,1,P_RM ),
     &              TRIE_LS(1,1,P_RQ ), TRIE_LS(1,1,P_RT ),
     &              TRIE_LS(1,1,P_dpm), TRIE_LS(1,1,P_dp ),
     &              TRIE_LS(1,1,P_dpn),
     &              TRIO_LS(1,1,P_TEM), TRIO_LS(1,1,P_TE ),
     &              TRIO_LS(1,1,P_Y  ), TRIO_LS(1,1,P_DIM),
     &              TRIO_LS(1,1,P_DI ), TRIO_LS(1,1,P_X  ),
     &              TRIO_LS(1,1,P_ZEM), TRIO_LS(1,1,P_ZE ),
     &              TRIO_LS(1,1,P_W  ), TRIO_LS(1,1,P_RM ),
     &              TRIO_LS(1,1,P_RQ ), TRIO_LS(1,1,P_RT ),
     &              TRIO_LS(1,1,P_dpm), TRIO_LS(1,1,P_dp ),
     &              TRIO_LS(1,1,P_dpn),
     &              FILTA,LS_NODE)
      endif
!
      DO J=1,LEN_TRIE_LS
         TRIE_LS(J,1,P_QM)=TRIE_LS(J,1,P_Q )
         TRIE_LS(J,2,P_QM)=TRIE_LS(J,2,P_Q )
         TRIE_LS(J,1,P_Q )=TRIE_LS(J,1,P_ZQ)
         TRIE_LS(J,2,P_Q )=TRIE_LS(J,2,P_ZQ)
      ENDDO
      DO J=1,LEN_TRIO_LS
         TRIO_LS(J,1,P_QM)=TRIO_LS(J,1,P_Q )
         TRIO_LS(J,2,P_QM)=TRIO_LS(J,2,P_Q )
         TRIO_LS(J,1,P_Q )=TRIO_LS(J,1,P_ZQ)
         TRIO_LS(J,2,P_Q )=TRIO_LS(J,2,P_ZQ)
      ENDDO
!--------------------------------------------
      else	! fwd_step next
!--------------------------------------------
        call do_dynamics_spectn2c(trie_ls,trio_ls)
!--------------------------------------------
      endif	! end fwd_step

!--------------------------------------------
! ---------------------------------------------------------------------
! transform new spectral into grid, then do filter in grid-point values
! ---------------------------------------------------------------------
!
      call spect_to_grid
     &      (trie_ls,trio_ls, 
     &       syn_gr_a_1,syn_gr_a_2,
     &       ls_node,ls_nodes,max_ls_nodes,
     &       lats_nodes_a,global_lats_a,lonsperlat,
     &       epse,epso,epsedn,epsodn,
     &       snnp1ev,snnp1od,plnev_a,plnod_a)
      call do_dynamics_syn2gridn(syn_gr_a_2,grid_gr,
     &                           global_lats_a,lonsperlat)
!
! do filter in the grid point values and advance time with update
! ---------------------------
      if(.not.fwd_step)then
        call do_dynamics_gridfilter(grid_gr,filta,filtb,
     &                              global_lats_a,lonsperlat)
      else	! fwd_step next
        call do_dynamics_gridn2c(grid_gr,
     &                           global_lats_a,lonsperlat)
      endif	! end fwd_step

      if( fwd_step ) then
        fwd_step = .false.
        dt   = deltim
        dt2  = cons2*dt
        rdt2 = 1./dt2

        if( gen_coord_hybrid ) then
          if( mass_dp ) then
            call get_cd_hyb_gcdp(dt)
          else
            call get_cd_hyb_gc(dt)
          endif
        else if( hybrid ) then
          call get_cd_hyb(dt)
        else
          call get_cd_sig(am,bm,dt,tov,sv)
        endif
      endif

! ------------------------------------------------
      endif 	! end not start_step
! ------------------------------------------------
      endif 	! only for fcst pes
!
!--------------------------------------------
!-- digital filter state collect
!--------------------------------------------
!      if (me == 0)
!     &print *,'in one loop,call gfs_dfi_coll,ldfi=',ldfi,'kdt=',kdt
      IF (ldfi) THEN
        call gfs_dficoll_dynamics(grid_gr,grid_gr_dfi)
      ENDIF
!
!
! =====================================================================
      if(.not.restart_step) then
!
!--------------------------------------------
      IF (lsout) THEN
!--------------------------------------------
CC
        CALL f_hpmstart(32,"wrtout_dynamics")
CC
        CALL countperf(0,18,0.)
c
        CALL WRTOUT_dynamics(PHOUR,FHOUR,ZHOUR,IDATE,
     &              TRIE_LS,TRIO_LS,grid_gr,
     &              SL,SI,
     &              ls_node,LS_NODES,MAX_LS_NODES,
     &              lats_nodes_a,global_lats_a,lonsperlat,
     &              snnp1ev,snnp1od,
!    &              CFHOUR1,snnp1ev,snnp1od,
     &              epsedn,epsodn,plnev_a,plnod_a,
     &              epse  ,epso  ,plnew_a,plnow_a,
     &              pdryini,'SIG.F')
!
        CALL f_hpmstop(32)
!!
        CALL countperf(1,18,0.)
!!
! ----------------------------------
      ENDIF ! if ls_out
! ----------------------------------
!!
!       print *,'kdt=',kdt,'nsres=',nsres,'write=',mod(kdt,nsres),
!     &   'lonsperlat=',lonsperlat
!
       IF (mod(kdt,nsres).eq.0.and.kdt.ne.0) THEN
!!
         CALL wrt_restart_dynamics(TRIE_LS,TRIO_LS,grid_gr,
     &        SI,fhour,idate,igen,pdryini,
     &        ls_node,ls_nodes,max_ls_nodes,
     &        global_lats_a,lonsperlat,lats_nodes_a,
     &        epse,epso,plnew_a,plnow_a,
     &        ens_nam,kdt,nfcstdate7)
!
       ENDIF
!
!-- end of restart step
       else
          restart_step=.false.
       endif
! =====================================================================
!
!  finish integration, then return
! -----------------------------------
!      print *,'in dyn one loop,bf endstep,end_step=',end_step,
!     &   'dfiend_step=',dfiend_step,'kdt=',kdt
      IF(end_step) THEN
          Cpl_flag = .true.
          RETURN
!
      else if(dfiend_step) then
!
!  dfi end step, return to dfi routine
!------------------------------------
       if(me==0) 
     &     print *,'in dyn one step, return aft dfi,kdt=',kdt
          RETURN
      END IF
! run for Cpl_flag == .true.
!---------------------------
      ELSE 
          grid_gr6 = grid_gr(:, 2 : lotgr6*2 + 1)

          CALL grid_to_spect_inp
     &        (grid_gr(1, g_gz ),    grid_gr(1, g_qm),
     &         grid_gr(1, g_uum),    grid_gr(1, g_vvm),
     &         grid_gr(1, g_ttm),    grid_gr(1, g_rm),
     &         grid_gr(1, g_dpm),
     &         trie_ls,  trio_ls,
     &         ls_node,      ls_nodes,      max_ls_nodes,
     &         lats_nodes_a, global_lats_a, lonsperlat,
     &         epse, epso, plnew_a, plnow_a, plnev_a, plnod_a,
     &         pwat, ptot, ptrc)

          call do_dynamics_spectn2m(trie_ls,trio_ls)

          CALL grid_to_spect_inp
     &        (grid_gr(1, g_gz ),    grid_gr(1, g_q),
     &         grid_gr(1, g_uu),    grid_gr(1, g_vv),
     &         grid_gr(1, g_tt),    grid_gr(1, g_rq),
     &         grid_gr(1, g_dp),
     &         trie_ls,  trio_ls,
     &         ls_node,      ls_nodes,      max_ls_nodes,
     &         lats_nodes_a, global_lats_a, lonsperlat,
     &         epse, epso, plnew_a, plnow_a, plnev_a, plnod_a,
     &         pwat, ptot, ptrc)

          call do_dynamics_spectn2c(trie_ls,trio_ls)

          call do_dynamics_gridc2n(grid_gr, global_lats_a, lonsperlat)

          Cpl_flag = .false.
! ---------------------
      END IF   ! cpl_flag
! ---------------------
      kdt = kdt + 1
! ===================================================================
! ===================================================================
! start nonlinear computation for dynamics 
! ===================================================================
! ===================================================================
      iter_max=0
      do lan=1,lats_node_a
         lat = global_lats_a(ipt_lats_node_a-1+lan)
         lons_lat = lonsperlat(lat)
         iter_max = max ( iter_max , (lons_lat+ngptcd-1)/ngptcd )
      enddo
!
      allocate ( spdlat(levs,iter_max ) )
      do k=1,levs
         spdmax_node(k) = cons0     !constant
      enddo
!
! =====================================================================
!----------------------------------------------------------
      if (me<num_pes_fcst) then
!----------------------------------------------------------
!
! transform total tendency in grid to spectral
!
      call spect_to_gridxy
     &    (trie_ls,trio_ls,
     &     syn_gr_a_1,syn_gr_a_2,
     &     dyn_gr_a_1,dyn_gr_a_2,
     &     ls_node,ls_nodes,max_ls_nodes,
     &     lats_nodes_a,global_lats_a,lonsperlat,
     &     pddev_a,pddod_a)
!
! -------------------------------------------------------------------
! do pressure gradient related derivatives.

      if( mass_dp ) then

        call do_dynamics_gridzz(grid_gr,global_lats_a,lonsperlat)

        call gridzz_to_spect
     &      (grid_gr,trie_ls,trio_ls,
     &       ls_node,ls_nodes,max_ls_nodes,
     &       lats_nodes_a,global_lats_a,lonsperlat,
     &       epse,epso,plnew_a,plnow_a)

        call spectpz_to_gridxy
     &      (trie_ls,trio_ls,
     &       pyn_gr_a_1,pyn_gr_a_2,
     &       ls_node,ls_nodes,max_ls_nodes,
     &       lats_nodes_a,global_lats_a,lonsperlat,
     &       epse,epso,epsedn,epsodn,
     &       snnp1ev,snnp1od,plnev_a,plnod_a)

      endif
! -------------------------------------------------------------------
!  get dpdt in grid point values for export state
!
      call do_dynamics_gridomega(syn_gr_a_2,dyn_gr_a_2,pyn_gr_a_2,
     &                     grid_gr,rcs2_a,global_lats_a,lonsperlat)
!
! =============================
      if( .not. ndslfv ) then
! =============================

!----------------------------------------------------------
      do lan=1,lats_node_a  
!
        lat = global_lats_a(ipt_lats_node_a-1+lan)
!
        lon_dim = lon_dims_a(lan)
        lons_lat = lonsperlat(lat)

        if( .not. gen_coord_hybrid ) then        
!$omp parallel do private(k,j)
          do k=1,levs
            do j=1,lons_lat
              syn_gr_a_2(j+(kst-2+k)*lon_dim,lan)=
     &        syn_gr_a_2(j+(kst-2+k)*lon_dim,lan)-tov(k)
            enddo
          enddo
        endif                                   
! --------------------------------------------------------------------

!$omp parallel do schedule(dynamic,1) private(lon)
!$omp+private(njeff,iblk)
!$omp+private(nvcn,xvcn)
        do lon=1,lons_lat,ngptcd
!!
          njeff = min(ngptcd,lons_lat-lon+1)
          iblk  = (lon-1)/ngptcd + 1
!
          CALL countperf(0,10,0.)

          if( gen_coord_hybrid ) then                                    ! hmhj
            if( thermodyn_id.eq.3 ) then                                 ! hmhj

              if( mass_dp ) then

                call gfidi_hyb_gchdp(lon_dim, njeff, lat,             ! hmhj 
     &               syn_gr_a_2(lon+(ksd   -1)*lon_dim,lan),          ! hmhj
     &               syn_gr_a_2(lon+(ksz   -1)*lon_dim,lan),          ! hmhj
     &               syn_gr_a_2(lon+(kst   -1)*lon_dim,lan),          ! hmhj
     &               syn_gr_a_2(lon+(ksu   -1)*lon_dim,lan),          ! hmhj
     &               syn_gr_a_2(lon+(ksv   -1)*lon_dim,lan),          ! hmhj
     &               syn_gr_a_2(lon+(ksr   -1)*lon_dim,lan),          ! hmhj
     &               syn_gr_a_2(lon+(ksdp  -1)*lon_dim,lan),          ! hmhj
     &               pyn_gr_a_2(lon+(kdpphi-1)*lon_dim,lan),          ! hmhj
     &               pyn_gr_a_2(lon+(kdplam-1)*lon_dim,lan),          ! hmhj
     &               pyn_gr_a_2(lon+(kzzphi-1)*lon_dim,lan),          ! hmhj
     &               pyn_gr_a_2(lon+(kzzlam-1)*lon_dim,lan),          ! hmhj
     &               rcs2_a(min(lat,latg-lat+1)),                       ! hmhj
     &               spdlat(1,iblk),                                    ! hmhj
     &               dt,nvcn,xvcn,                                      ! hmhj
     &               dyn_gr_a_2(lon+(kdtphi-1)*lon_dim,lan),          ! hmhj
     &               dyn_gr_a_2(lon+(kdtlam-1)*lon_dim,lan),          ! hmhj
     &               dyn_gr_a_2(lon+(kdrphi-1)*lon_dim,lan),          ! hmhj
     &               dyn_gr_a_2(lon+(kdrlam-1)*lon_dim,lan),          ! hmhj
     &               dyn_gr_a_2(lon+(kdulam-1)*lon_dim,lan),          ! hmhj
     &               dyn_gr_a_2(lon+(kdvlam-1)*lon_dim,lan),          ! hmhj
     &               dyn_gr_a_2(lon+(kduphi-1)*lon_dim,lan),          ! hmhj
     &               dyn_gr_a_2(lon+(kdvphi-1)*lon_dim,lan),          ! hmhj
     &               anl_gr_a_2(lon+(kaps  -1)*lon_dim,lan),          ! hmhj
     &               anl_gr_a_2(lon+(kadp  -1)*lon_dim,lan),          ! hmhj
     &               anl_gr_a_2(lon+(kat   -1)*lon_dim,lan),          ! hmhj
     &               anl_gr_a_2(lon+(kar   -1)*lon_dim,lan),          ! hmhj
     &               anl_gr_a_2(lon+(kau   -1)*lon_dim,lan),          ! hmhj
     &               anl_gr_a_2(lon+(kav   -1)*lon_dim,lan),          ! hmhj
     &               szdrdt(lon,lan),zfirst)                          ! fyang

              else      ! not-mass_dp

              call gfidi_hyb_gc_h(lon_dim, njeff, lat,                   ! hmhj
     &               syn_gr_a_2(lon+(ksd   -1)*lon_dim,lan),          ! hmhj
     &               syn_gr_a_2(lon+(kst   -1)*lon_dim,lan),          ! hmhj
     &               syn_gr_a_2(lon+(ksz   -1)*lon_dim,lan),          ! hmhj
     &               syn_gr_a_2(lon+(ksu   -1)*lon_dim,lan),          ! hmhj
     &               syn_gr_a_2(lon+(ksv   -1)*lon_dim,lan),          ! hmhj
     &               syn_gr_a_2(lon+(ksr   -1)*lon_dim,lan),          ! hmhj
     &               syn_gr_a_2(lon+(kspphi-1)*lon_dim,lan),          ! hmhj
     &               syn_gr_a_2(lon+(ksplam-1)*lon_dim,lan),          ! hmhj
     &               syn_gr_a_2(lon+(ksq   -1)*lon_dim,lan),          ! hmhj
     &               syn_gr_a_2(lon+(kzsphi-1)*lon_dim,lan),          ! hmhj
     &               syn_gr_a_2(lon+(kzslam-1)*lon_dim,lan),          ! hmhj
     &               rcs2_a(min(lat,latg-lat+1)),                       ! hmhj
     &               spdlat(1,iblk),                                    ! hmhj
     &               dt,nvcn,xvcn,                                  	! hmhj
     &               dyn_gr_a_2(lon+(kdtphi-1)*lon_dim,lan),          ! hmhj
     &               dyn_gr_a_2(lon+(kdtlam-1)*lon_dim,lan),          ! hmhj
     &               dyn_gr_a_2(lon+(kdrphi-1)*lon_dim,lan),          ! hmhj
     &               dyn_gr_a_2(lon+(kdrlam-1)*lon_dim,lan),          ! hmhj
     &               dyn_gr_a_2(lon+(kdulam-1)*lon_dim,lan),          ! hmhj
     &               dyn_gr_a_2(lon+(kdvlam-1)*lon_dim,lan),          ! hmhj
     &               dyn_gr_a_2(lon+(kduphi-1)*lon_dim,lan),          ! hmhj
     &               dyn_gr_a_2(lon+(kdvphi-1)*lon_dim,lan),          ! hmhj
     &               anl_gr_a_2(lon+(kaps  -1)*lon_dim,lan),          ! hmhj
     &               anl_gr_a_2(lon+(kat   -1)*lon_dim,lan),          ! hmhj
     &               anl_gr_a_2(lon+(kar   -1)*lon_dim,lan),          ! hmhj
     &               anl_gr_a_2(lon+(kau   -1)*lon_dim,lan),          ! hmhj
     &               anl_gr_a_2(lon+(kav   -1)*lon_dim,lan),          ! hmhj
     &               szdrdt(lon,lan),zfirst)                          ! fyang

              endif 	! end of mass_dp 

            else  ! not enthalpy next

              call gfidi_hyb_gc(lon_dim, njeff, lat,
     &               syn_gr_a_2(lon+(ksd   -1)*lon_dim,lan),
     &               syn_gr_a_2(lon+(kst   -1)*lon_dim,lan),
     &               syn_gr_a_2(lon+(ksz   -1)*lon_dim,lan),
     &               syn_gr_a_2(lon+(ksu   -1)*lon_dim,lan),
     &               syn_gr_a_2(lon+(ksv   -1)*lon_dim,lan),
     &               syn_gr_a_2(lon+(ksr   -1)*lon_dim,lan),
     &               syn_gr_a_2(lon+(kspphi-1)*lon_dim,lan),
     &               syn_gr_a_2(lon+(ksplam-1)*lon_dim,lan),
     &               syn_gr_a_2(lon+(ksq   -1)*lon_dim,lan),
     &               syn_gr_a_2(lon+(kzsphi-1)*lon_dim,lan),
     &               syn_gr_a_2(lon+(kzslam-1)*lon_dim,lan),
     x               rcs2_a(min(lat,latg-lat+1)),
     &               spdlat(1,iblk),
     &               dt,nvcn,xvcn,
     &               dyn_gr_a_2(lon+(kdtphi-1)*lon_dim,lan),
     &               dyn_gr_a_2(lon+(kdtlam-1)*lon_dim,lan),
     &               dyn_gr_a_2(lon+(kdrphi-1)*lon_dim,lan),
     &               dyn_gr_a_2(lon+(kdrlam-1)*lon_dim,lan),
     &               dyn_gr_a_2(lon+(kdulam-1)*lon_dim,lan),
     &               dyn_gr_a_2(lon+(kdvlam-1)*lon_dim,lan),
     &               dyn_gr_a_2(lon+(kduphi-1)*lon_dim,lan),
     &               dyn_gr_a_2(lon+(kdvphi-1)*lon_dim,lan),
     &               anl_gr_a_2(lon+(kaps  -1)*lon_dim,lan),
     &               anl_gr_a_2(lon+(kat   -1)*lon_dim,lan),
     &               anl_gr_a_2(lon+(kar   -1)*lon_dim,lan),
     &               anl_gr_a_2(lon+(kau   -1)*lon_dim,lan),
     &               anl_gr_a_2(lon+(kav   -1)*lon_dim,lan),
     &               szdrdt(lon,lan),zfirst)

            endif                                                        ! hmhj

          else if( hybrid ) then                                         ! hmhj

            call gfidi_hyb(lon_dim, njeff, lat,
     &               syn_gr_a_2(lon+(ksd   -1)*lon_dim,lan),
     &               syn_gr_a_2(lon+(kst   -1)*lon_dim,lan),
     &               syn_gr_a_2(lon+(ksz   -1)*lon_dim,lan),
     &               syn_gr_a_2(lon+(ksu   -1)*lon_dim,lan),
     &               syn_gr_a_2(lon+(ksv   -1)*lon_dim,lan),
     &               syn_gr_a_2(lon+(ksr   -1)*lon_dim,lan),
     &               syn_gr_a_2(lon+(kspphi-1)*lon_dim,lan),
     &               syn_gr_a_2(lon+(ksplam-1)*lon_dim,lan),
     &               syn_gr_a_2(lon+(ksq   -1)*lon_dim,lan),
     &               syn_gr_a_2(lon+(kzsphi-1)*lon_dim,lan),
     &               syn_gr_a_2(lon+(kzslam-1)*lon_dim,lan),
     &               rcs2_a(min(lat,latg-lat+1)),
     &               spdlat(1,iblk),
     &               dt,nvcn,xvcn,
     &               dyn_gr_a_2(lon+(kdtphi-1)*lon_dim,lan),
     &               dyn_gr_a_2(lon+(kdtlam-1)*lon_dim,lan),
     &               dyn_gr_a_2(lon+(kdrphi-1)*lon_dim,lan),
     &               dyn_gr_a_2(lon+(kdrlam-1)*lon_dim,lan),
     &               dyn_gr_a_2(lon+(kdulam-1)*lon_dim,lan),
     &               dyn_gr_a_2(lon+(kdvlam-1)*lon_dim,lan),
     &               dyn_gr_a_2(lon+(kduphi-1)*lon_dim,lan),
     &               dyn_gr_a_2(lon+(kdvphi-1)*lon_dim,lan),
     &               anl_gr_a_2(lon+(kaps  -1)*lon_dim,lan),
     &               anl_gr_a_2(lon+(kat   -1)*lon_dim,lan),
     &               anl_gr_a_2(lon+(kar   -1)*lon_dim,lan),
     &               anl_gr_a_2(lon+(kau   -1)*lon_dim,lan),
     &               anl_gr_a_2(lon+(kav   -1)*lon_dim,lan),
     &               szdrdt(lon,lan),zfirst)
        
          else

          call gfidi_sig(lon_dim, njeff, lat,
     &               syn_gr_a_2(lon+(ksd   -1)*lon_dim,lan),
     &               syn_gr_a_2(lon+(kst   -1)*lon_dim,lan),
     &               syn_gr_a_2(lon+(ksz   -1)*lon_dim,lan),
     &               syn_gr_a_2(lon+(ksu   -1)*lon_dim,lan),
     &               syn_gr_a_2(lon+(ksv   -1)*lon_dim,lan),
     &               syn_gr_a_2(lon+(ksr   -1)*lon_dim,lan),
     &               syn_gr_a_2(lon+(kspphi-1)*lon_dim,lan),
     &               syn_gr_a_2(lon+(ksplam-1)*lon_dim,lan),
     &               syn_gr_a_2(lon+(ksq   -1)*lon_dim,lan),
     &               syn_gr_a_2(lon+(kzsphi-1)*lon_dim,lan),
     &               syn_gr_a_2(lon+(kzslam-1)*lon_dim,lan),
     &               rcs2_a(min(lat,latg-lat+1)),
     &               del,rdel2,ci,tov,spdlat(1,iblk),
     &               dt,sl,nvcn,xvcn,
     &               dyn_gr_a_2(lon+(kdtphi-1)*lon_dim,lan),
     &               dyn_gr_a_2(lon+(kdtlam-1)*lon_dim,lan),
     &               dyn_gr_a_2(lon+(kdrphi-1)*lon_dim,lan),
     &               dyn_gr_a_2(lon+(kdrlam-1)*lon_dim,lan),
     &               dyn_gr_a_2(lon+(kdulam-1)*lon_dim,lan),
     &               dyn_gr_a_2(lon+(kdvlam-1)*lon_dim,lan),
     &               dyn_gr_a_2(lon+(kduphi-1)*lon_dim,lan),
     &               dyn_gr_a_2(lon+(kdvphi-1)*lon_dim,lan),
     &               anl_gr_a_2(lon+(kaps  -1)*lon_dim,lan),
     &               anl_gr_a_2(lon+(kat   -1)*lon_dim,lan),
     &               anl_gr_a_2(lon+(kar   -1)*lon_dim,lan),
     &               anl_gr_a_2(lon+(kau   -1)*lon_dim,lan),
     &               anl_gr_a_2(lon+(kav   -1)*lon_dim,lan))

          endif                                                          ! hmhj

          CALL countperf(1,10,0.)
!
        enddo   !lon

! ---------------------------------------------------------------
!
        iblk=1
        do lon=1,lons_lat,ngptcd
          do k=1,levs
            spdmax_node(k)=max(spdmax_node(k),spdlat(k,iblk))
          enddo
          iblk=iblk+1
        enddo
!
! --------------------------------------------------------
      enddo   ! end of lan
!
! -------------------------------------------------------------------
!  update the grid point values
!
      if(.not.mass_dp) then
        call do_dynamics_gridaddlapgz(anl_gr_a_2,syn_gr_a_2,
     &                                  global_lats_a,lonsperlat)
      endif
      call do_dynamics_gridupdate(grid_gr,anl_gr_a_2,dt2,
     &                                  global_lats_a,lonsperlat)
!
! ===================================================
      else     ! ndslfv next
! ===================================================

!     call ndslfv_monoadvh (grid_gr,global_lats_a,lonsperlat,dt,kdt)
      call ndslfv_monoadvh2(grid_gr,global_lats_a,lonsperlat,dt)

      call do_dynamics_advhn2anl(grid_gr,anl_gr_a_2,
     &                                  global_lats_a,lonsperlat)

      call do_dynamics_gridc2syq(grid_gr,syn_gr_syq,
     &                             global_lats_a,lonsperlat)

      do lan=1,lats_node_a
!
        lat = global_lats_a(ipt_lats_node_a-1+lan)
!
        lon_dim = lon_dims_a(lan)
        lons_lat = lonsperlat(lat)

!$omp parallel do schedule(dynamic,1) private(lon)
!$omp+private(njeff,iblk)
!$omp+private(nvcn,xvcn)
        do lon=1,lons_lat,ngptcd
!!
          njeff=min(ngptcd,lons_lat-lon+1)
          iblk  = (lon-1)/ngptcd + 1
!
! mass_dp and thermodyn_id=3
!
                call gfidi_gchdp_noadv_noq (lon_dim, njeff, lat,        ! hmhj
     &               syn_gr_a_2(lon+(ksd   -1)*lon_dim,lan),          ! hmhj
     &               syn_gr_a_2(lon+(ksz   -1)*lon_dim,lan),          ! hmhj
     &               syn_gr_a_2(lon+(kst   -1)*lon_dim,lan),          ! hmhj
     &               syn_gr_a_2(lon+(ksu   -1)*lon_dim,lan),          ! hmhj
     &               syn_gr_a_2(lon+(ksv   -1)*lon_dim,lan),          ! hmhj
     &               syn_gr_syq(lon                   ,lan),          ! hmhj
     &               syn_gr_a_2(lon+(ksdp  -1)*lon_dim,lan),          ! hmhj
     &               syn_gr_a_2(lon+(ksq   -1)*lon_dim,lan),
     &               syn_gr_a_2(lon+(kzsphi-1)*lon_dim,lan),          ! hmhj
     &               syn_gr_a_2(lon+(kzslam-1)*lon_dim,lan),          ! hmhj
     &               pyn_gr_a_2(lon+(kdpphi-1)*lon_dim,lan),          ! hmhj
     &               pyn_gr_a_2(lon+(kdplam-1)*lon_dim,lan),          ! hmhj
     &               pyn_gr_a_2(lon+(kzzphi-1)*lon_dim,lan),          ! hmhj
     &               pyn_gr_a_2(lon+(kzzlam-1)*lon_dim,lan),          ! hmhj
     &               rcs2_a(min(lat,latg-lat+1)),                       ! hmhj
     &               spdlat(1,iblk),dt,                                 ! hmhj
     &               dyn_gr_a_2(lon+(kdtphi-1)*lon_dim,lan),          ! hmhj
     &               dyn_gr_a_2(lon+(kdtlam-1)*lon_dim,lan),          ! hmhj
     &               dyn_gr_a_2(lon+(kdulam-1)*lon_dim,lan),          ! hmhj
     &               dyn_gr_a_2(lon+(kdvlam-1)*lon_dim,lan),          ! hmhj
     &               dyn_gr_a_2(lon+(kduphi-1)*lon_dim,lan),          ! hmhj
     &               dyn_gr_a_2(lon+(kdvphi-1)*lon_dim,lan),          ! hmhj
     &               anl_gr_a_2(lon+(kaps  -1)*lon_dim,lan),          ! hmhj
     &               anl_gr_a_2(lon+(kadp  -1)*lon_dim,lan),          ! hmhj
     &               anl_gr_a_2(lon+(kat   -1)*lon_dim,lan),          ! hmhj
     &               anl_gr_a_2(lon+(kau   -1)*lon_dim,lan),          ! hmhj
     &               anl_gr_a_2(lon+(kav   -1)*lon_dim,lan),          ! hmhj
     &               anl_gr_a_2(lon+(kap2  -1)*lon_dim,lan),          ! hmhj
     &               pdot(lon,1,lan))                                 ! hmhj

        enddo   ! end of lon
! --------------------------------------------------------
        iblk=1
        do lon=1,lons_lat,ngptcd
          do k=1,levs
            spdmax_node(k)=max(spdmax_node(k),spdlat(k,iblk))
          enddo
          iblk=iblk+1
        enddo
! --------------------------------------------------------
      enddo   ! end of lan


!  update the grid point values without vertical advection
!
      call do_dynamics_gridupdate(grid_gr,anl_gr_a_2,dt2,
     &                                  global_lats_a,lonsperlat)

!  update the grid point with vertical advection
!
      call ndslfv_monoadvv (grid_gr,pdot,
     &                        global_lats_a,lonsperlat,dt)

! ================
      endif     ! end of ndslfv
! ================
! -------------------------------------------------------------------
!  update the grid point values of p and dp
!
!hmhj redefine dp due to physics ....
!hmhj if( mass_dp ) then
!hmhj   if( process_split ) then
!hmhj     stp = 0
!hmhj   else
!hmhj     stp = 1
!hmhj   endif
!hmhj   call do_dynamics_gridp   (grid_gr,global_lats_a,lonsperlat,stp)
!hmhj else
!hmhj   call do_dynamics_gridpdpn(grid_gr,global_lats_a,lonsperlat)
!hmhj endif
      stp=1
      if( process_split ) stp=0
      call do_dynamics_gridpdp(grid_gr,global_lats_a,lonsperlat,stp)
!
! ---------------------------------------------------------------------
! gather all speed max
!
      spdmax_nodem=spdmax_node
      call mpi_gather(spdmax_nodem,levs,MPI_R_MPI,
     &                spdmax_nodesm,levs,MPI_R_MPI,
     &                0,MC_COMP,ierr)
!WY bug fix.
!      spdmax_nodes=spdmax_nodesm
      IF(me == 0) THEN
          spdmax_nodes(:,:)=spdmax_nodesm(:, :)
      END IF
!
!sela call mpi_barrier (mpi_comm_world,ierr)
!
!--------------------------------------------
      if ( me .eq. 0 ) then
!--------------------------------------------
!
         do k=1,levs
            spdmax(k) = cons0     !constant
            do node=1,nodes
               spdmax(k)=max(spdmax(k),spdmax_nodes(k,node))
            enddo
            spdmax(k)=sqrt(spdmax(k))
         enddo
!
         print*,'in do_dynamics_one_loop for spdmx at kdt=',kdt
         print 100,(spdmax(k),k=1,levs)
100      format(' spdmx(001:010)=',10f5.0,:/' spdmx(011:020)=',10f5.0,
     &        :/' spdmx(021:030)=',10f5.0,:/' spdmx(031:040)=',10f5.0,
     &        :/' spdmx(041:050)=',10f5.0,:/' spdmx(051:060)=',10f5.0,
     &        :/' spdmx(061:070)=',10f5.0,:/' spdmx(071:080)=',10f5.0,
     &        :/' spdmx(081:090)=',10f5.0,:/' spdmx(091:100)=',10f5.0,
     &        :/' spdmx(101:110)=',10f5.0,:/' spdmx(111:120)=',10f5.0,
     &        :/' spdmx(121:130)=',10f5.0,:/' spdmx(131:140)=',10f5.0,
     &        :/' spdmx(141:150)=',10f5.0,:/' spdmx(151:160)=',10f5.0)
!
!--------------------------------------------
      endif
!--------------------------------------------
!
      call mpi_bcast(spdmax,levs,mpi_real8,
     &               0,MC_COMP,ierr)
!
      deallocate ( spdlat )
!
      if(zfirst) zfirst=.false.
!!
      if( process_split ) then
        call do_dynamics_gridm2p(grid_gr,global_lats_a,lonsperlat)
      else
        call do_dynamics_gridn2p(grid_gr,global_lats_a,lonsperlat)
      endif

!!
!--------------------------------------------
      endif ! only for forecast pes
!--------------------------------------------
      RETURN
      END
