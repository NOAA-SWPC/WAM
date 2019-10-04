      SUBROUTINE input_fields(cread, PDRYINI,TRIE_LS,TRIO_LS,grid_gr,
     &                        LS_NODE,LS_NODES,MAX_LS_NODES,
     &                        global_lats_a,lonsperlat,
     &                        epse,epso,epsedn,epsodn,plnev_a,plnod_a,
     &                        plnew_a,plnow_a,lats_nodes_a,
     &                        pwat,ptot,ptrc,SNNP1EV,SNNP1OD)
!!
!! Aug 2010    Sarah Lu, modified to compute tracer global sum
!! 20100908    J. WANG   remove gfsio module
!! 20110220    H. Juang  change to have mass dp and NDSL
!! 20120913    J. Wang   add sigio option
!! 20130930    H. Juang  change the assignment of grid to n then n-1

      use gfs_dyn_resol_def
      use gfs_dyn_layout1
      use gfs_dyn_gg_def
      use gfs_dyn_vert_def
      use gfs_dyn_date_def
      use namelist_dynamics_def
      use gfs_dyn_mpi_def
      IMPLICIT NONE
!!
 
cmy fix pdryini type
cmy      REAL(KIND=KIND_EVOD) PDRYINI
      REAL(KIND=kind_grid) PDRYINI
      CHARACTER (len=*)   :: CREAD
      REAL(KIND=KIND_EVOD) TRIE_LS(LEN_TRIE_LS,2,LOTLS)
      REAL(KIND=KIND_EVOD) TRIO_LS(LEN_TRIO_LS,2,LOTLS)
      REAL(KIND=KIND_GRID) GRID_GR(lonf,lats_node_a_max,lotgr)
      REAL(KIND=KIND_EVOD) EPSE   (LEN_TRIE_LS)
      REAL(KIND=KIND_EVOD) EPSO   (LEN_TRIO_LS)
      REAL(KIND=KIND_EVOD) EPSEDN (LEN_TRIE_LS)
      REAL(KIND=KIND_EVOD) EPSODN (LEN_TRIO_LS)
      REAL(KIND=KIND_EVOD) PLNEV_a(LEN_TRIE_LS,latg2)
      REAL(KIND=KIND_EVOD) PLNOD_a(LEN_TRIO_LS,latg2)
      REAL(KIND=KIND_EVOD) PLNEW_a(LEN_TRIE_LS,latg2)
      REAL(KIND=KIND_EVOD) PLNOW_a(LEN_TRIO_LS,latg2)
      REAL(KIND=KIND_EVOD) snnp1ev(LEN_TRIE_LS)
      REAL(KIND=KIND_EVOD) snnp1od(LEN_TRIO_LS)
!
      real(kind=kind_grid) zsg(lonf,lats_node_a)
      real(kind=kind_grid) psg(lonf,lats_node_a)
      real(kind=kind_grid) dpg(lonf,lats_node_a,levs)
      real(kind=kind_grid) ttg(lonf,lats_node_a,levs)
      real(kind=kind_grid) uug(lonf,lats_node_a,levs)
      real(kind=kind_grid) vvg(lonf,lats_node_a,levs)
      real(kind=kind_grid) rqg(lonf,lats_node_a,levh)

      REAL(KIND=KIND_GRID) pwat   (lonf,lats_node_a)
      REAL(KIND=KIND_GRID) ptot   (lonf,lats_node_a)
      REAL(KIND=KIND_GRID) ptrc   (lonf,lats_node_a,ntrac)          !glbsum
!
      integer global_lats_a(latg), lonsperlat(latg)
 
cmy bug fix on dimension of ls_node
      INTEGER              LS_NODE (LS_DIM*3)
      INTEGER              LS_NODES(LS_DIM,NODES)
      INTEGER          MAX_LS_NODES(NODES)
      integer            lats_nodes_a(nodes)
!
      INTEGER              IERR,IPRINT,J,JDT,K,L,LOCL,N,i
      REAL(KIND=KIND_EVOD) TEE1(LEVS)
      REAL(KIND=KIND_EVOD)  YE1(LEVS)
      INTEGER              INDLSEV,JBASEV
      INTEGER              INDLSOD,JBASOD
      REAL(KIND=KIND_EVOD), parameter :: CONS0=0.0, CONS2=2.0,
     &                                   CONS600=600.0
      integer		lan, lat, lons_lat, jlonf, step
!
      if(me.eq.0) PRINT  9876,FHOUR,idate
 9876 FORMAT(1H ,'FHOUR IN input_fields ',F6.2,
     & ' idate no yet read in',4(1x,i4))
      IPRINT = 0
c$$$  IF ( ME .EQ. 0 ) IPRINT = 1
!
      if (me == 0) write(0,*)'input field, cread=',cread,'ntoz=',ntoz,
     &  'nemsio_in=',nemsio_in

      if(nemsio_in) then
        if (me == 0) write(0,*)'read in nemsio file,cread=',trim(cread)

        CALL TREADEO_nemsio(cread,IDATE,
     &                      trie_ls,trio_ls,
     &                      zsg, psg, ttg, uug, vvg, rqg, dpg,
     X                      LS_NODE,LS_NODES,MAX_LS_NODES,
     X                      PDRYINI,IPRINT,
     &                      global_lats_a,lats_nodes_a,lonsperlat,
     &                      epse, epso, plnew_a, plnow_a, 
     &                      plnev_a, plnod_a, snnp1ev, snnp1od,
     &                      pwat, ptot, ptrc)            !glbsum
      
        do j = 1, lats_node_a
          grid_gr(:,j,g_gz) = zsg(:,j)
          grid_gr(:,j,g_q ) = psg(:,j)
          grid_gr(:,j,g_dp:g_dp+levs-1) = dpg(:,j,1:levs)
          grid_gr(:,j,g_tt:g_tt+levs-1) = ttg(:,j,1:levs)
          grid_gr(:,j,g_uu:g_uu+levs-1) = uug(:,j,1:levs)
          grid_gr(:,j,g_vv:g_vv+levs-1) = vvg(:,j,1:levs)
          grid_gr(:,j,g_rq:g_rq+levh-1) = rqg(:,j,1:levh)
        enddo

      else
        print *,'read in sigio file,cread=',trim(cread)
        CALL TREADEO(IDATE,trie_ls,trio_ls,
     &               grid_gr,
     X               LS_NODE,LS_NODES,MAX_LS_NODES,
     X               PDRYINI,IPRINT,
     &               global_lats_a,lats_nodes_a,lonsperlat, cread,
     &               epse, epso, epsedn,epsodn,plnew_a, plnow_a,
     &               plnev_a, plnod_a,snnp1ev,snnp1od)
      endif

      fhini=fhour
      if(me.eq.0) PRINT 9877, FHOUR
 9877 FORMAT(1H ,'FHOUR AFTER TREAD',F6.2)
 
      if (me .eq. 0) write(0,*)' fhini=',fhini,'last_fcst_pe=',
     &     last_fcst_pe,'fhrot=',fhrot

!---------------------------------------------------------------
      if(fhini.eq.fhrot) THEN
!set n time level values to n-1 time
! spectral
!       print *,' set time level n to time level n-1 '

       trie_ls(:,:,P_qm )=trie_ls(:,:,P_q )
       trie_ls(:,:,p_dpm:p_dpm+levs-1)=trie_ls(:,:,p_dp:p_dp+levs-1)
       trie_ls(:,:,p_tem:p_tem+levs-1)=trie_ls(:,:,p_te:p_te+levs-1)
       trie_ls(:,:,p_dim:p_dim+levs-1)=trie_ls(:,:,p_di:p_di+levs-1)
       trie_ls(:,:,p_zem:p_zem+levs-1)=trie_ls(:,:,p_ze:p_ze+levs-1)
       trio_ls(:,:,P_qm )=trio_ls(:,:,P_q )
       trio_ls(:,:,p_dpm:p_dpm+levs-1)=trio_ls(:,:,p_dp:p_dp+levs-1)
       trio_ls(:,:,p_tem:p_tem+levs-1)=trio_ls(:,:,p_te:p_te+levs-1)
       trio_ls(:,:,p_dim:p_dim+levs-1)=trio_ls(:,:,p_di:p_di+levs-1)
       trio_ls(:,:,p_zem:p_zem+levs-1)=trio_ls(:,:,p_ze:p_ze+levs-1)

       if( .not. ndslfv ) then
         trie_ls(:,:,p_rm:p_rm+levh-1)=trie_ls(:,:,p_rq:p_rq+levh-1)
         trio_ls(:,:,p_rm:p_rm+levh-1)=trio_ls(:,:,p_rq:p_rq+levh-1)
       endif

! grid
        grid_gr(:,:,g_qm )=grid_gr(:,:,g_q )
        grid_gr(:,:,g_dpm:g_dpm+levs-1)=grid_gr(:,:,g_dp:g_dp+levs-1)
        grid_gr(:,:,g_ttm:g_ttm+levs-1)=grid_gr(:,:,g_tt:g_tt+levs-1)
        grid_gr(:,:,g_uum:g_uum+levs-1)=grid_gr(:,:,g_uu:g_uu+levs-1)
        grid_gr(:,:,g_vvm:g_vvm+levs-1)=grid_gr(:,:,g_vv:g_vv+levs-1)
        grid_gr(:,:,g_rm :g_rm +levh-1)=grid_gr(:,:,g_rq:g_rq+levh-1)

      endif
!
! fill up n+1 grid_gr in case of internal2export used.
!
        grid_gr(:,:,g_zq)=grid_gr(:,:,g_q )
        grid_gr(:,:,g_dpn:g_dpn+levs-1)=grid_gr(:,:,g_dp:g_dp+levs-1)
        grid_gr(:,:,g_t  :g_t  +levs-1)=grid_gr(:,:,g_tt:g_tt+levs-1)
        grid_gr(:,:,g_u  :g_u  +levs-1)=grid_gr(:,:,g_uu:g_uu+levs-1)
        grid_gr(:,:,g_v  :g_v  +levs-1)=grid_gr(:,:,g_vv:g_vv+levs-1)
        grid_gr(:,:,g_rt :g_rt +levh-1)=grid_gr(:,:,g_rq:g_rq+levh-1)
!
! fill up p grid_gr
!
        grid_gr(:,:,g_zqp)=grid_gr(:,:,g_q )
        grid_gr(:,:,g_dpp:g_dpp+levs-1)=grid_gr(:,:,g_dp:g_dp+levs-1)
        grid_gr(:,:,g_ttp:g_ttp+levs-1)=grid_gr(:,:,g_tt:g_tt+levs-1)
        grid_gr(:,:,g_uup:g_uup+levs-1)=grid_gr(:,:,g_uu:g_uu+levs-1)
        grid_gr(:,:,g_vvp:g_vvp+levs-1)=grid_gr(:,:,g_vv:g_vv+levs-1)
        grid_gr(:,:,g_rqp:g_rqp+levh-1)=grid_gr(:,:,g_rq:g_rq+levh-1)
!     
!--------------------------------------------------------
!!
      RETURN
      END
