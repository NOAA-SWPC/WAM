!-----------------------------------------------------------------------
!
      MODULE MODULE_CLOCKTIMES
!
!-----------------------------------------------------------------------
!
!***  List the clocktime counters for the various parts of
!***  the integration and print them as desired.
!
!-----------------------------------------------------------------------
!
      USE MODULE_KINDS
!
!-----------------------------------------------------------------------
!
      PRIVATE
!
      PUBLIC :: INTEGRATION_TIMERS                                      &
               ,PRINT_CLOCKTIMES                                        &
               ,TIMERS                                                  &
               ,cbcst_tim,pbcst_tim
!
!-----------------------------------------------------------------------
!
      TYPE INTEGRATION_TIMERS
!
        REAL(kind=KDBL) :: total_integ_tim,totalsum_tim
!
        REAL(kind=KDBL) :: solver_dyn_tim,solver_phy_tim
!
        REAL(kind=KDBL) :: adv1_tim,adv2_tim,bocoh_tim,bocov_tim        &
                          ,cdwdt_tim,cdzdt_tim,consts_tim               &
                          ,ddamp_tim,dht_tim                            &
                          ,exch_dyn,exch_phy                            &
                          ,exch_tim                                     &
                          ,fftfhn_tim,fftfwn_tim                        &
                          ,hdiff_tim,mono_tim                           &
                          ,pdtsdt_tim,pgforce_tim,poavhn_tim            &
                          ,polehn_tim,polewn_tim                        &
                          ,prefft_tim,presmud_tim                       &
                          ,solver_init_tim                              &
                          ,swaphn_tim,swapwn_tim                        &
                          ,updatet_tim                                  &
                          ,updateuv_tim                                 &
                          ,updates_tim                                  &
                          ,vsound_tim,vtoa_tim
!
        REAL(kind=KDBL) :: adjppt_tim,cucnvc_tim                        &
                          ,gsmdrive_tim,h_to_v_tim,gfs_phy_tim          &
                          ,phy_sum_tim                                  &
                          ,pole_swap_tim,radiation_tim,rdtemp_tim       &
                          ,turbl_tim                                    &
                          ,cltend_tim,rfupdate_tim,tqadjust_tim
!
        REAL(kind=KDBL) :: domain_run_1                                 &
                          ,domain_run_2                                 &
                          ,domain_run_3                                 &
                          ,pc_cpl_run_cpl1                              &
                          ,pc_cpl_run_cpl2                              &
                          ,pc_cpl_run_cpl3                              &
                          ,pc_cpl_run_cpl4                              &
                          ,pc_cpl_run_cpl5                              &
                          ,cpl1_recv_tim                                &
                          ,cpl2_send_tim                                &
                          ,cpl2_comp_tim                                &
                          ,cpl2_wait_tim                                &
                          ,parent_bookkeep_moving_tim                   &
                          ,parent_update_moving_tim                     &
                          ,t0_recv_move_tim                             &
                          ,read_moving_child_topo_tim                   &
                          ,barrier_move_tim,pscd_tim,pscd1_tim          &
                          ,pscd2_tim,pscd3_tim,pscd4_tim

!
!-----------------------------------------------------------------------
!***  Associated with moving nests
!-----------------------------------------------------------------------
!
        REAL(kind=KDBL) :: update_interior_from_nest_tim                &
                          ,update_interior_from_parent_tim
!
      END TYPE INTEGRATION_TIMERS
!
!-----------------------------------------------------------------------
!
      TYPE(INTEGRATION_TIMERS),DIMENSION(:),ALLOCATABLE,TARGET :: TIMERS   !<-- Timers for each domain
!
      REAL(kind=KDBL),DIMENSION(99) :: cbcst_tim,pbcst_tim
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!
      SUBROUTINE PRINT_CLOCKTIMES(NTIMESTEP                             &
                                 ,MY_DOMAIN_ID                          &
                                 ,MYPE                                  &
                                 ,NPE_PRINT                             &
                                 ,TIMERS_DOMAIN)
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument variables
!------------------------
!
      INTEGER(kind=KINT),INTENT(IN) :: NTIMESTEP                        &  !<-- Forecast timestep
                                      ,MY_DOMAIN_ID                     &  !<-- The domain's ID
                                      ,MYPE                             &  !<-- The task ID
                                      ,NPE_PRINT                           !<-- ID of task providing clocktime diagnostics
!
      TYPE(INTEGRATION_TIMERS),TARGET,INTENT(INOUT) :: TIMERS_DOMAIN       !<-- Assorted clocktime timers for current domain
!
!---------------------
!***  Local variables
!---------------------
!
      TYPE(INTEGRATION_TIMERS),POINTER :: TD
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      TD=>TIMERS_DOMAIN                                                    !<-- Abbreviate name of this domain's timers
!
!-----------------------------------------------------------------------
!
      td%totalsum_tim=td%adv1_tim                                       &
                  +td%bocoh_tim                                         &
                  +td%bocov_tim                                         &
                  +td%cdwdt_tim                                         &
                  +td%cdzdt_tim                                         &
                  +td%dht_tim                                           &
                  +td%ddamp_tim                                         &
                  +td%exch_dyn                                          &
                  +td%fftfhn_tim                                        &
                  +td%fftfwn_tim                                        &
                  +td%hdiff_tim                                         &
                  +td%pdtsdt_tim                                        &
                  +td%pgforce_tim                                       &
                  +td%poavhn_tim                                        &
                  +td%polehn_tim                                        &
                  +td%polewn_tim                                        &
                  +td%swaphn_tim                                        &
                  +td%swapwn_tim                                        &
                  +td%updatet_tim                                       &
                  +td%updateuv_tim                                      &
                  +td%updates_tim                                       &
                  +td%vsound_tim                                        &
                  +td%vtoa_tim
!
      td%totalsum_tim=td%totalsum_tim                                   &
                  +td%cucnvc_tim                                        &
                  +td%exch_phy                                          &
                  +td%gsmdrive_tim                                      &
                  +td%h_to_v_tim                                        &
                  +td%pole_swap_tim                                     &
                  +td%radiation_tim                                     &
                  +td%rdtemp_tim                                        &
                  +td%turbl_tim                                         &
                  +td%cltend_tim                                        &
                  +td%rfupdate_tim                                      &
                  +td%tqadjust_tim
!
      td%totalsum_tim=td%totalsum_tim                                   &
                  +td%solver_init_tim 
!
!-----------------------------------------------------------------------
!***  The designated MPI task writes clocktimes for its work.
!-----------------------------------------------------------------------
!
      IF(MYPE==NPE_PRINT)THEN
!
        write(0,*)' '
        write(0,FMT='(" Clocktimes for domain #",I2.2)') my_domain_id
!
        write(0,FMT='(" ntsd= ",I6," total_integration_tim=  ",g12.5)') ntimestep,td%total_integ_tim
!
        write(0,FMT='("   solver_init=          ",g12.5," pct= ",f7.2)') td%solver_init_tim &
                 ,td%solver_init_tim/td%total_integ_tim*100.
        write(0,FMT='("   consts=               ",g12.5," pct= ",f7.2)') td%consts_tim &
                 ,td%consts_tim/td%total_integ_tim*100.
!
        write(0,*)' DYNAMICS'
!
        write(0,FMT='("   solver_dyn_total=     ",g12.5," pct= ",f7.2)') td%solver_dyn_tim &
                 ,td%solver_dyn_tim/td%total_integ_tim*100.
        write(0,FMT='("   solver_dyn_w/o_exch=  ",g12.5," pct= ",f7.2)') td%solver_dyn_tim-td%exch_dyn &
                 ,(td%solver_dyn_tim-td%exch_dyn)/td%total_integ_tim*100.
        write(0,FMT='("   pgforce=              ",g12.5," pct= ",f7.2)') td%pgforce_tim &
                 ,td%pgforce_tim/td%total_integ_tim*100.
        write(0,FMT='("   dht=                  ",g12.5," pct= ",f7.2)') td%dht_tim &
                 ,td%dht_tim/td%total_integ_tim*100.
        write(0,FMT='("   ddamp=                ",g12.5," pct= ",f7.2)') td%ddamp_tim &
                ,td%ddamp_tim/td%total_integ_tim*100.
        write(0,FMT='("   pdtsdt=               ",g12.5," pct= ",f7.2)') td%pdtsdt_tim &
                ,td%pdtsdt_tim/td%total_integ_tim*100.
        write(0,FMT='("   vtoa=                 ",g12.5," pct= ",f7.2)') td%vtoa_tim &
                ,td%vtoa_tim/td%total_integ_tim*100.
        write(0,FMT='("   adv1=                 ",g12.5," pct= ",f7.2)') td%adv1_tim &
                ,td%adv1_tim/td%total_integ_tim*100.
!
        if(td%adv2_tim/=0.)then
          write(0,FMT='("   adv2=                 ",g12.5," pct= ",f7.2)') td%adv2_tim &
               ,td%adv2_tim/td%total_integ_tim*100.
        endif
!
        if(td%mono_tim/=0.)then
          write(0,FMT='("   mono=                 ",g12.5," pct= ",f7.2)') td%mono_tim &
               ,td%mono_tim/td%total_integ_tim*100.
        endif
!
        write(0,FMT='("   cdzdt=                ",g12.5," pct= ",f7.2)') td%cdzdt_tim &
                ,td%cdzdt_tim/td%total_integ_tim*100.
        write(0,FMT='("   cdwdt=                ",g12.5," pct= ",f7.2)') td%cdwdt_tim &
                ,td%cdwdt_tim/td%total_integ_tim*100.
        write(0,FMT='("   vsound=               ",g12.5," pct= ",f7.2)') td%vsound_tim &
                ,td%vsound_tim/td%total_integ_tim*100. 
        write(0,FMT='("   hdiff=                ",g12.5," pct= ",f7.2)') td%hdiff_tim &
                ,td%hdiff_tim/td%total_integ_tim*100.
        write(0,FMT='("   bocoh=                ",g12.5," pct= ",f7.2)') td%bocoh_tim &
                ,td%bocoh_tim/td%total_integ_tim*100.
        write(0,FMT='("   bocov=                ",g12.5," pct= ",f7.2)') td%bocov_tim &
                ,td%bocov_tim/td%total_integ_tim*100.
        write(0,FMT='("   updatet=              ",g12.5," pct= ",f7.2)') td%updatet_tim &
                ,td%updatet_tim/td%total_integ_tim*100.
        write(0,FMT='("   updateuv=             ",g12.5," pct= ",f7.2)') td%updateuv_tim &
                ,td%updateuv_tim/td%total_integ_tim*100.
        write(0,FMT='("   updates=              ",g12.5," pct= ",f7.2)') td%updates_tim &
                ,td%updates_tim/td%total_integ_tim*100.

        if(td%prefft_tim/=0.)then
          write(0,FMT='("   prefft=               ",g12.5," pct= ",f7.2)') td%prefft_tim &
                ,td%prefft_tim/td%total_integ_tim*100.
          write(0,FMT='("   fftfhn=               ",g12.5," pct= ",f7.2)') td%fftfhn_tim &
                ,td%fftfhn_tim/td%total_integ_tim*100.
          write(0,FMT='("   fftfwn=               ",g12.5," pct= ",f7.2)') td%fftfwn_tim &
                ,td%fftfwn_tim/td%total_integ_tim*100.
          write(0,FMT='("   polewn=               ",g12.5," pct= ",f7.2)') td%polewn_tim &
                ,td%polewn_tim/td%total_integ_tim*100.
          write(0,FMT='("   poavhn=               ",g12.5," pct= ",f7.2)') td%poavhn_tim &
                ,td%poavhn_tim/td%total_integ_tim*100.
        endif
!
        if(td%presmud_tim/=0.)then
          write(0,FMT='("   presmud=              ",g12.5," pct= ",f7.2)') td%presmud_tim &
                ,td%presmud_tim/td%total_integ_tim*100.
        endif
        write(0,*)' PHYSICS '
!
        write(0,FMT='("   solver_phy_total=     ",g12.5," pct= ",f7.2)') td%solver_phy_tim &
                 ,td%solver_phy_tim/td%total_integ_tim*100.
        write(0,FMT='("   solver_phy_w/o_exch=  ",g12.5," pct= ",f7.2)') td%solver_phy_tim-td%exch_phy &
                 ,(td%solver_phy_tim-td%exch_phy)/td%total_integ_tim*100.
        write(0,FMT='("   cucnvc=               ",g12.5," pct= ",f7.2)') td%cucnvc_tim &
                ,td%cucnvc_tim/td%total_integ_tim*100.
        write(0,FMT='("   gsmdrive=             ",g12.5," pct= ",f7.2)') td%gsmdrive_tim &
                ,td%gsmdrive_tim/td%total_integ_tim*100.
        write(0,FMT='("   cltend=               ",g12.5," pct= ",f7.2)') td%cltend_tim &
                ,td%cltend_tim/td%total_integ_tim*100.
        write(0,FMT='("   rime_factor_update=   ",g12.5," pct= ",f7.2)') td%rfupdate_tim &
                ,td%rfupdate_tim/td%total_integ_tim*100.
        write(0,FMT='("   tqadjust=             ",g12.5," pct= ",f7.2)') td%tqadjust_tim &
                ,td%tqadjust_tim/td%total_integ_tim*100.
        write(0,FMT='("   radiation=            ",g12.5," pct= ",f7.2)') td%radiation_tim &
                ,td%radiation_tim/td%total_integ_tim*100.
        write(0,FMT='("   rdtemp=               ",g12.5," pct= ",f7.2)') td%rdtemp_tim &
                ,td%rdtemp_tim/td%total_integ_tim*100.
        write(0,FMT='("   turbl=                ",g12.5," pct= ",f7.2)') td%turbl_tim &
                ,td%turbl_tim/td%total_integ_tim*100.
        write(0,FMT='("   h_to_v=               ",g12.5," pct= ",f7.2)') td%h_to_v_tim &
                ,td%h_to_v_tim/td%total_integ_tim*100.
!
        if(td%pole_swap_tim/=0.)then
          write(0,FMT='("   pole_swap=            ",g12.5," pct= ",f7.2)') td%pole_swap_tim &
                ,td%pole_swap_tim/td%total_integ_tim*100.
        endif
!
        write(0,*)' EXCHANGE TIMES '
!
        write(0,FMT='("   exch_dyn=             ",g12.5," pct= ",f7.2)') td%exch_dyn &
                ,td%exch_dyn/td%total_integ_tim*100.
!
        write(0,FMT='("   exch_phy=             ",g12.5," pct= ",f7.2)') td%exch_phy &
                ,td%exch_phy/td%total_integ_tim*100.
!
        td%exch_tim=td%exch_dyn+td%exch_phy
        write(0,FMT='("   exch_tim=             ",g12.5," pct= ",f7.2)') td%exch_tim &
                ,td%exch_tim/td%total_integ_tim*100.
!
        if(td%swaphn_tim/=0.)then
          write(0,FMT='("   swaphn=               ",g12.5," pct= ",f7.2)') td%swaphn_tim &
                ,td%swaphn_tim/td%total_integ_tim*100.
        endif
!
        if(td%swapwn_tim/=0.)then
          write(0,FMT='("   swapwn=               ",g12.5," pct= ",f7.2)') td%swapwn_tim &
                ,td%swapwn_tim/td%total_integ_tim*100.
        endif
!
        if(td%polehn_tim/=0.)then
          write(0,FMT='("   polehn=               ",g12.5," pct= ",f7.2)') td%polehn_tim &
                ,td%polehn_tim/td%total_integ_tim*100.
        endif
!
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE PRINT_CLOCKTIMES
!
!-----------------------------------------------------------------------
!
      END MODULE MODULE_CLOCKTIMES
!
!-----------------------------------------------------------------------
