      SUBROUTINE input_fields_slg(n1,n2, PDRYINI,TRIE_LS,TRIO_LS,
     &                 grid_gr,
     &                 LS_NODE,LS_NODES,MAX_LS_NODES,SNNP1EV,SNNP1OD,
     &                 epse,epso, cread, cread2, restart_run,
     &                 global_lats_a,lonsperlat,
!    &                 global_lats_a,lonsperlat,epsedn,epsodn,
     &                 plnew_a,plnow_a,lats_nodes_a)
!    &                 plnev_a,plnod_a,plnew_a,plnow_a,lats_nodes_a)
!    &                 pwat,ptot,ptrc)
!
!-----
!  revision history
!
! Feb, 2015 Jun Wang    Add option for nemsio
! Feb  2015 S Moorthi - changed jun's modifications and fixed for reading in nemsio
!                       to include corect stack etc for semi-Lagrangian model
! Nov  2015 S Moorthi-  logic to initialize from nemsio history file
!
      USE gfs_dyn_machine,       ONLY : KIND_EVOD, kind_grid
      USE gfs_dyn_resol_def,     ONLY : latg, latg2, levh, levs, 
     &                                  p_rm, p_q, p_ze, p_rq, p_te, 
     &                                  p_di, p_zem, p_tem, p_dim, 
     &                                  p_gz, p_qm, lonf, ntrac, lotls
     &,                                 lotgr
      USE gfs_dyn_layout1,       ONLY : me, nodes, ls_dim, len_trie_ls, 
     &                                  len_trio_ls, ipt_lats_node_a
     &,                                 lats_node_a_max

      USE gfs_dyn_date_def,      ONLY : fhour, idate
      USE namelist_dynamics_def, ONLY : fhini, nemsio_in, gg_tracers
     &,                                 fhrot
      USE gfs_dyn_layout1      , ONLY : lats_node_a
      use layout_grid_tracers ,  only : rgt_a,rgt_h,xhalo,yhalo
!
      IMPLICIT NONE
!!
 
      REAL(KIND=KIND_EVOD) PDRYINI
      INTEGER              N1,N2
      CHARACTER (len=*)   :: CREAD, CREAD2
!     REAL(KIND=KIND_EVOD) TRIE_LS(LEN_TRIE_LS,2,11*LEVS+3*LEVH+6)
!     REAL(KIND=KIND_EVOD) TRIO_LS(LEN_TRIO_LS,2,11*LEVS+3*LEVH+6)

      REAL(KIND=KIND_EVOD) TRIE_LS(LEN_TRIE_LS,2,lotls)
      REAL(KIND=KIND_EVOD) TRIO_LS(LEN_TRIO_LS,2,lotls)
      real(kind=kind_evod), dimension(len_trie_ls) :: epse, snnp1ev
      real(kind=kind_evod), dimension(len_trio_ls) :: epso, snnp1od
      REAL(KIND=KIND_GRID) GRID_GR(lonf,lats_node_a_max,lotgr)

!     REAL(KIND=KIND_EVOD) EPSEDN (LEN_TRIE_LS)
!     REAL(KIND=KIND_EVOD) EPSODN (LEN_TRIO_LS)
!     REAL(KIND=KIND_EVOD) PLNEV_a(LEN_TRIE_LS,latg2)
!     REAL(KIND=KIND_EVOD) PLNOD_a(LEN_TRIO_LS,latg2)

      REAL(KIND=KIND_EVOD) PLNEW_a(LEN_TRIE_LS,latg2)
      REAL(KIND=KIND_EVOD) PLNOW_a(LEN_TRIO_LS,latg2)
!
!     real(kind=kind_evod) zsg(lonf,lats_node_a)
!     real(kind=kind_evod) psg(lonf,lats_node_a)
!     real(kind=kind_evod) dpg(lonf,lats_node_a,levs)
!     real(kind=kind_evod) ttg(lonf,lats_node_a,levs)
!     real(kind=kind_evod) uug(lonf,lats_node_a,levs)
!     real(kind=kind_evod) vvg(lonf,lats_node_a,levs)
!     real(kind=kind_evod) rqg(lonf,lats_node_a,levh)

!     REAL(KIND=KIND_evod) pwat   (lonf,lats_node_a)
!     REAL(KIND=KIND_evod) ptot   (lonf,lats_node_a)
!     REAL(KIND=KIND_evod) ptrc   (lonf,lats_node_a,ntrac)          !glbsum
!
      integer, dimension(latg)  :: global_lats_a, lonsperlat
      integer, dimension(nodes) :: max_ls_nodes, lats_nodes_a
!
      INTEGER              LS_NODE (LS_DIM*3), LS_NODES(LS_DIM,NODES)
!
c$$$      INTEGER                LOTS,LOTD,LOTA
c$$$      PARAMETER            ( LOTS = 5*LEVS+1*LEVH+3 )
c$$$      PARAMETER            ( LOTD = 6*LEVS+2*LEVH+0 )
c$$$      PARAMETER            ( LOTA = 3*LEVS+1*LEVH+1 )
c$$$      INTEGER   P_GZ,P_ZEM,P_DIM,P_TEM,P_RM,P_QM
c$$$      INTEGER   P_ZE,P_DI,P_TE,P_RQ,P_Q,P_DLAM,P_DPHI,P_ULN,P_VLN
c$$$      INTEGER   P_W,P_X,P_Y,P_RT,P_ZQ
c$$$      PARAMETER(P_GZ  = 0*LEVS+0*LEVH+1,  !      GZE/O(LNTE/OD,2),
c$$$     X          P_ZEM = 0*LEVS+0*LEVH+2,  !     ZEME/O(LNTE/OD,2,LEVS),
c$$$     X          P_DIM = 1*LEVS+0*LEVH+2,  !     DIME/O(LNTE/OD,2,LEVS),
c$$$     X          P_TEM = 2*LEVS+0*LEVH+2,  !     TEME/O(LNTE/OD,2,LEVS),
c$$$     X          P_RM  = 3*LEVS+0*LEVH+2,  !      RME/O(LNTE/OD,2,LEVH),
c$$$     X          P_QM  = 3*LEVS+1*LEVH+2,  !      QME/O(LNTE/OD,2),
c$$$     X          P_ZE  = 3*LEVS+1*LEVH+3,  !      ZEE/O(LNTE/OD,2,LEVS),
c$$$     X          P_DI  = 4*LEVS+1*LEVH+3,  !      DIE/O(LNTE/OD,2,LEVS),
c$$$     X          P_TE  = 5*LEVS+1*LEVH+3,  !      TEE/O(LNTE/OD,2,LEVS),
c$$$     X          P_RQ  = 6*LEVS+1*LEVH+3,  !      RQE/O(LNTE/OD,2,LEVH),
c$$$     X          P_Q   = 6*LEVS+2*LEVH+3,  !       QE/O(LNTE/OD,2),
c$$$     X          P_DLAM= 6*LEVS+2*LEVH+4,  !  DPDLAME/O(LNTE/OD,2),
c$$$     X          P_DPHI= 6*LEVS+2*LEVH+5,  !  DPDPHIE/O(LNTE/OD,2),
c$$$     X          P_ULN = 6*LEVS+2*LEVH+6,  !     ULNE/O(LNTE/OD,2,LEVS),
c$$$     X          P_VLN = 7*LEVS+2*LEVH+6,  !     VLNE/O(LNTE/OD,2,LEVS),
c$$$     X          P_W   = 8*LEVS+2*LEVH+6,  !       WE/O(LNTE/OD,2,LEVS),
c$$$     X          P_X   = 9*LEVS+2*LEVH+6,  !       XE/O(LNTE/OD,2,LEVS),
c$$$     X          P_Y   =10*LEVS+2*LEVH+6,  !       YE/O(LNTE/OD,2,LEVS),
c$$$     X          P_RT  =11*LEVS+2*LEVH+6,  !      RTE/O(LNTE/OD,2,LEVH),
c$$$     X          P_ZQ  =11*LEVS+3*LEVH+6)  !      ZQE/O(LNTE/OD,2)
!     REAL(KIND=KIND_EVOD) TEE1(LEVS)
!     REAL(KIND=KIND_EVOD)  YE1(LEVS)

      INTEGER              IERR,IPRINT,J,JDT,K,L,LOCL,N,i,kk
     &,                    INDLSEV,JBASEV, INDLSOD,JBASOD, lan, lat
     &,                    lons_lat

      REAL(KIND=KIND_EVOD), parameter :: CONS0=0.0, CONS2=2.0,
     &                                   CONS600=600.0
      LOGICAL LSLAG, restart_run
!
      if(me.eq.0) PRINT  9876,N1,N2,FHOUR,idate
 9876 FORMAT(1H ,'N1,N2,FHOUR IN spect_fields ',2(I4,1X),F6.2,
     & ' idate no yet read in',4(1x,i4))
      IPRINT = 0
!$$$  IF ( ME .EQ. 0 ) IPRINT = 1
!
      if (me == 0) write(0,*)' cread=',cread

      if(nemsio_in) then
        if (me == 0) write(0,*)'read in nemsio file,cread=',trim(cread)

        CALL TREADEO_nemsio_slg(cread,IDATE,1,restart_run,
     &               TRIE_LS(1,1,P_GZ ), TRIE_LS(1,1,P_QM ),
     &               TRIE_LS(1,1,P_TEM), TRIE_LS(1,1,P_DIM),
     &               TRIE_LS(1,1,P_ZEM), TRIE_LS(1,1,P_RM ),
     &               TRIO_LS(1,1,P_GZ ), TRIO_LS(1,1,P_QM ),
     &               TRIO_LS(1,1,P_TEM), TRIO_LS(1,1,P_DIM),
     &               TRIO_LS(1,1,P_ZEM), TRIO_LS(1,1,P_RM ),
     &               grid_gr,
     &               LS_NODE,LS_NODES,MAX_LS_NODES,
     &               PDRYINI,IPRINT,
     &               global_lats_a,lats_nodes_a,lonsperlat,
     &               epse, epso, plnew_a, plnow_a,
     &               snnp1ev, snnp1od)

!    &               plnev_a, plnod_a, snnp1ev, snnp1od)
!    &               pwat, ptot, ptrc, slg_flag)            !glbsum

!       write(0,*)'finished reading in nemsio file,cread=',
!    &  trim(cread)
!    &,' p_gz=',TRIE_LS(1,1,P_GZ ),' p_qm=',TRIE_LS(1,1,P_QM ),
!    &' p_tem=',TRIE_LS(1,1,P_TEM)
      else
        CALL TREADEO_slg(N1,FHOUR,IDATE,
     &               TRIE_LS(1,1,P_GZ ), TRIE_LS(1,1,P_QM ),
     &               TRIE_LS(1,1,P_TEM), TRIE_LS(1,1,P_DIM),
     &               TRIE_LS(1,1,P_ZEM), TRIE_LS(1,1,P_RM ),
     &               TRIO_LS(1,1,P_GZ ), TRIO_LS(1,1,P_QM ),
     &               TRIO_LS(1,1,P_TEM), TRIO_LS(1,1,P_DIM),
     &               TRIO_LS(1,1,P_ZEM), TRIO_LS(1,1,P_RM ),
     &               LS_NODE,
     &               SNNP1EV,SNNP1OD,PDRYINI,IPRINT,
     &               cread)
      endif
 
      if(me == 0) PRINT 9877, N1,FHOUR
 9877 FORMAT(1H ,'N1,FHOUR AFTER TREAD',1(I4,1X),F6.2)
 
!     fhini = fhour
      if (me == 0) write(0,*)' fhini=',fhini,' restart_run=',
     &                       restart_run
      CALL RMS_spect(TRIE_LS(1,1,P_QM ), TRIE_LS(1,1,P_DIM),
     &               TRIE_LS(1,1,P_TEM), TRIE_LS(1,1,P_ZEM),
     &               TRIE_LS(1,1,P_RM ),
     &               TRIO_LS(1,1,P_QM ), TRIO_LS(1,1,P_DIM),
     &               TRIO_LS(1,1,P_TEM), TRIO_LS(1,1,P_ZEM),
     &               TRIO_LS(1,1,P_RM ),
     &               LS_NODES,MAX_LS_NODES)
!---------------------------------------------------------------
      if(.NOT. restart_run .and. fhini == fhrot) THEN  !set n time level values to n-1 time
        do i=1,len_trie_ls
           trie_ls(i,1,P_q )=trie_ls(i,1,P_qm )
           trie_ls(i,2,P_q )=trie_ls(i,2,P_qm )
        enddo
        do i=1,len_trio_ls
           trio_ls(i,1,P_q )=trio_ls(i,1,P_qm )
           trio_ls(i,2,P_q )=trio_ls(i,2,P_qm )
        enddo
 
!$omp parallel do private(i,k)
        do k=1,levs
          do i=1,len_trie_ls
            trie_ls(i,1,P_te +k-1)=trie_ls(i,1,P_tem +k-1)
            trie_ls(i,2,P_te +k-1)=trie_ls(i,2,P_tem +k-1)
 
            trie_ls(i,1,P_di +k-1)=trie_ls(i,1,P_dim +k-1)
            trie_ls(i,2,P_di +k-1)=trie_ls(i,2,P_dim +k-1)
 
            trie_ls(i,1,P_ze +k-1)=trie_ls(i,1,P_zem +k-1)
            trie_ls(i,2,P_ze +k-1)=trie_ls(i,2,P_zem +k-1)
          enddo
          do i=1,len_trio_ls
            trio_ls(i,1,P_te +k-1) = trio_ls(i,1,P_tem+k-1)
            trio_ls(i,2,P_te +k-1) = trio_ls(i,2,P_tem+k-1)
 
            trio_ls(i,1,P_di +k-1) = trio_ls(i,1,P_dim+k-1)
            trio_ls(i,2,P_di +k-1) = trio_ls(i,2,P_dim+k-1)
 
            trio_ls(i,1,P_ze +k-1) = trio_ls(i,1,P_zem+k-1)
            trio_ls(i,2,P_ze +k-1) = trio_ls(i,2,P_zem+k-1)
          enddo
        enddo
        if (gg_tracers .and. nemsio_in) then
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
!$omp parallel do private(i,k)
          do k=1,levh
            do i=1,len_trie_ls
              trie_ls(i,1,P_rq +k-1) = trie_ls(i,1,P_rm +k-1)
              trie_ls(i,2,P_rq +k-1) = trie_ls(i,2,P_rm +k-1)
            enddo
            do i=1,len_trio_ls
              trio_ls(i,1,P_rq +k-1) = trio_ls(i,1,P_rm+k-1)
              trio_ls(i,2,P_rq +k-1) = trio_ls(i,2,P_rm+k-1)
            enddo
          enddo
        endif
!--------------------------------------------------------
      else
!--------------------------------------------------------
        IPRINT = 0
!$$$      IF ( ME == 0 ) IPRINT = 1
      if (me == 0) write(0,*)' cread2=',cread2
        if (nemsio_in) then
          CALL TREADEO_nemsio_slg(cread2,IDATE,2,restart_run,
     &                 TRIE_LS(1,1,P_GZ ), TRIE_LS(1,1,P_Q ),
     &                 TRIE_LS(1,1,P_TE), TRIE_LS(1,1,P_DI),
     &                 TRIE_LS(1,1,P_ZE), TRIE_LS(1,1,P_RQ ),
     &                 TRIO_LS(1,1,P_GZ), TRIO_LS(1,1,P_Q ),
     &                 TRIO_LS(1,1,P_TE), TRIO_LS(1,1,P_DI),
     &                 TRIO_LS(1,1,P_ZE), TRIO_LS(1,1,P_RQ ),
     &                 grid_gr,
     &                 LS_NODE,LS_NODES,MAX_LS_NODES,
     &                 PDRYINI,IPRINT,
     &                 global_lats_a,lats_nodes_a,lonsperlat,
     &                 epse, epso, plnew_a, plnow_a,
     &                 snnp1ev, snnp1od)
        else
          CALL TREADEO_slg(N2,FHOUR,IDATE,
     &                 TRIE_LS(1,1,P_GZ), TRIE_LS(1,1,P_Q ),
     &                 TRIE_LS(1,1,P_TE), TRIE_LS(1,1,P_DI),
     &                 TRIE_LS(1,1,P_ZE), TRIE_LS(1,1,P_RQ),
     &                 TRIO_LS(1,1,P_GZ), TRIO_LS(1,1,P_Q ),
     &                 TRIO_LS(1,1,P_TE), TRIO_LS(1,1,P_DI),
     &                 TRIO_LS(1,1,P_ZE), TRIO_LS(1,1,P_RQ),
     &                 LS_NODE,
     &                 SNNP1EV,SNNP1OD,PDRYINI,IPRINT,
     &                 cread2)
          if(me == 0) PRINT 9878, N2,FHOUR
 9878     FORMAT(1H ,'N2,FHOUR AFTER TREAD',1(I4,1X),F6.2)
        endif
      endif

      END SUBROUTINE input_fields_slg
