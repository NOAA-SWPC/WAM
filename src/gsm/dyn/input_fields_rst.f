      SUBROUTINE input_fields_rst(gread,gread2,cread, cread2, 
     &                 PDRYINI,TRIE_LS,TRIO_LS,
     &                 grid_gr,LS_NODE,LS_NODES,MAX_LS_NODES,
     &                 global_lats_a,lonsperlat,
     &                 epse,epso,plnev_a,plnod_a,plnew_a,plnow_a,
     &                 snnp1ev,snnp1od,lats_nodes_a)
!!
!program log
!  20100205  J. WANG     Read in input restart files without computing 
!                        pwat nad ptot
!  20100908  J. WANG     remove gfsio module
!  20110220  H. Juang    remove some un-necessary name in treads_nemsio
!  20170510  W. YANG     modified for the WAM restart.
!
      use gfs_dyn_resol_def
      use gfs_dyn_layout1
      use gfs_dyn_gg_def
      use gfs_dyn_vert_def
      use gfs_dyn_date_def
      use namelist_dynamics_def
      use gfs_dyn_mpi_def
      use gfs_dyn_physcons, cp => con_cp , rd => con_rd,
     &    rerth => con_rerth, grav => con_g
      use gfs_dyn_coordinate_def
      IMPLICIT NONE
!!
      REAL(KIND=kind_grid) PDRYINI
      CHARACTER (len=*)   :: CREAD, CREAD2,gread,gread2
      REAL(KIND=KIND_EVOD) TRIE_LS(LEN_TRIE_LS,2,LOTLS)
      REAL(KIND=KIND_EVOD) TRIO_LS(LEN_TRIO_LS,2,LOTLS)
      REAL(KIND=KIND_GRID) GRID_GR(lonf,lats_node_a_max,lotgr)
      REAL(KIND=KIND_EVOD) EPSE   (LEN_TRIE_LS)
      REAL(KIND=KIND_EVOD) EPSO   (LEN_TRIO_LS)
      REAL(KIND=KIND_EVOD) PLNEV_a(LEN_TRIE_LS,latg2)
      REAL(KIND=KIND_EVOD) PLNOD_a(LEN_TRIO_LS,latg2)
      REAL(KIND=KIND_EVOD) PLNEW_a(LEN_TRIE_LS,latg2)
      REAL(KIND=KIND_EVOD) PLNOW_a(LEN_TRIO_LS,latg2)
      real(kind=kind_evod),intent(in) :: snnp1ev(len_trie_ls)
      real(kind=kind_evod),intent(in) :: snnp1od(len_trio_ls)
!
      real(kind=kind_grid) zsg(lonf,lats_node_a)
      real(kind=kind_grid) psg(lonf,lats_node_a)
      real(kind=kind_grid) dpg(lonf,lats_node_a,levs)
      real(kind=kind_grid) ttg(lonf,lats_node_a,levs)
      real(kind=kind_grid) uug(lonf,lats_node_a,levs)
      real(kind=kind_grid) vvg(lonf,lats_node_a,levs)
      real(kind=kind_grid) rqg(lonf,lats_node_a,levh)
!
      integer nblck
      integer global_lats_a(latg), lonsperlat(latg)
!
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
      integer		lan, lat, lons_lat, jlonf,nnl,nn,kk,lon
      integer          indev1,indev2,indev,indod1,indod2,indod
      REAL(KIND=KIND_EVOD)  ga2
      real(kind=kind_mpi_r),allocatable :: trieo_ls_nodes(:,:,:,:)
      real(kind=kind_mpi_r),allocatable :: trieo_ls_node(:,:,:)
      integer ioproc,lenrec,ii
!
      include 'function2'
!
!------------------------------------------------------------------
!
      if(me.eq.0) PRINT  9876,FHOUR,idate
 9876 FORMAT(1H ,'FHOUR IN input_fields ',F6.2,
     & ' idate no yet read in',4(1x,i4))
      IPRINT = 0
!$$$  IF ( ME .EQ. 0 ) IPRINT = 1
!
!--- n-1 time step grid 
! 
!     print *,'in inputfile_rst,gread=',gread,'gread2=',gread2,
!    & 'cread=',cread,'cread2=',cread2
!     if (me .eq. 0) write(0,*)' gread=',gread,'ntoz=',ntoz
        CALL TREADG_nemsio(gread,FHOUR,IDATE,
     &               zsg, psg, ttg, uug, vvg, rqg, dpg,
     X               PDRYINI,IPRINT,
     &               global_lats_a,lats_nodes_a,lonsperlat) 

       do j=1,lats_node_a
       grid_gr(:,j,g_gz) = zsg(:,j)
       grid_gr(:,j,g_qm) = psg(:,j)
       grid_gr(:,j,g_dpm:g_dpm+levs-1) = dpg(:,j,1:levs)
       grid_gr(:,j,g_ttm:g_ttm+levs-1) = ttg(:,j,1:levs)
       grid_gr(:,j,g_uum:g_uum+levs-1) = uug(:,j,1:levs)
       grid_gr(:,j,g_vvm:g_vvm+levs-1) = vvg(:,j,1:levs)
       grid_gr(:,j,g_rm :g_rm +levh-1) = rqg(:,j,1:levh)
       enddo
!
!---------------------------------------------------------------
!--------------------------------------------------------
!
!--- n time step grid 
! 
      if (me .eq. 0) write(0,*)' gread2=',gread2
          CALL TREADG_nemsio(gread2,fhour,idate,
     &                 zsg, psg, ttg, uug, vvg, rqg, dpg,
     X                 PDRYINI,IPRINT,
     &                 global_lats_a,lats_nodes_a,lonsperlat) 
!
       do j=1,lats_node_a
!      grid_gr(:,j,g_gz) = zsg(:,j)
       grid_gr(:,j,g_q ) = psg(:,j)
       grid_gr(:,j,g_dp :g_dp +levs-1) = dpg(:,j,1:levs)
       grid_gr(:,j,g_tt :g_tt +levs-1) = ttg(:,j,1:levs)
       grid_gr(:,j,g_uu :g_uu +levs-1) = uug(:,j,1:levs)
       grid_gr(:,j,g_vv :g_vv +levs-1) = vvg(:,j,1:levs)
       grid_gr(:,j,g_rq :g_rq +levh-1) = rqg(:,j,1:levh)
       enddo
!
! =======================================================================
!
!-------------------------------------------------------------------------
!---read in n-1 time step spectral file
      if (me .eq. 0) write(0,*)'in input, sread1, cread=',trim(cread)
          CALL TREADS_nemsio(cread,FHOUR,IDATE,
     &                       trie_ls, trio_ls,
     &                       LS_NODE)

        trie_ls(:,:,p_qm)=trie_ls(:,:,p_q )
        trie_ls(:,:,p_dpm:p_dpm+levs-1)=trie_ls(:,:,p_dp:p_dp+levs-1)
        trie_ls(:,:,p_tem:p_tem+levs-1)=trie_ls(:,:,p_te:p_te+levs-1)
        trie_ls(:,:,p_dim:p_dim+levs-1)=trie_ls(:,:,p_di:p_di+levs-1)
        trie_ls(:,:,p_zem:p_zem+levs-1)=trie_ls(:,:,p_ze:p_ze+levs-1)
        trio_ls(:,:,p_qm)=trio_ls(:,:,p_q )
        trio_ls(:,:,p_dpm:p_dpm+levs-1)=trio_ls(:,:,p_dp:p_dp+levs-1)
        trio_ls(:,:,p_tem:p_tem+levs-1)=trio_ls(:,:,p_te:p_te+levs-1)
        trio_ls(:,:,p_dim:p_dim+levs-1)=trio_ls(:,:,p_di:p_di+levs-1)
        trio_ls(:,:,p_zem:p_zem+levs-1)=trio_ls(:,:,p_ze:p_ze+levs-1)
        if ( .not. ndslfv ) then
          trie_ls(:,:,p_rm :p_rm +levh-1)=trie_ls(:,:,p_rq:p_rq+levh-1)
          trio_ls(:,:,p_rm :p_rm +levh-1)=trio_ls(:,:,p_rq:p_rq+levh-1)
        endif
!
!---
      fhini=fhour
      if(me.eq.0) PRINT 9877, FHOUR
 9877 FORMAT(1H ,'FHOUR AFTER TREAD',F6.2)

!     if (me .eq. 0) write(0,*)' fhini=',fhini,'last_fcst_pe=',
!    &     last_fcst_pe,'fhrot=',fhrot
!     if (me<=last_fcst_pe) then
!       CALL RMS_spect(TRIE_LS(1,1,P_QM ), TRIE_LS(1,1,P_DIM),
!    X             TRIE_LS(1,1,P_TEM), TRIE_LS(1,1,P_ZEM),
!    X             TRIE_LS(1,1,P_RM ),
!    X             TRIO_LS(1,1,P_QM ), TRIO_LS(1,1,P_DIM),
!    X             TRIO_LS(1,1,P_TEM), TRIO_LS(1,1,P_ZEM),
!    X             TRIO_LS(1,1,P_RM ),
!    X             LS_NODES,MAX_LS_NODES)
!     endif

!
!--- N time step spectral 
      if (me .eq. 0) write(0,*)'in input, sread2, cread=',trim(cread2)
          CALL TREADS_nemsio(cread2,FHOUR,IDATE,
     &                       trie_ls, trio_ls,
     &                       LS_NODE)
!
        if(me.eq.0) PRINT 9878, FHOUR
 9878   FORMAT(1H ,'FHOUR AFTER TREAD',F6.2)
!
!--------------------------------------------------------------
! laplacian terrain for divergence
!
        ga2=grav/(rerth*rerth)
        do locl=1,ls_max_node
                L=ls_node(locl)
           jbasev=ls_node(locl+ls_dim)
           indev1 = indlsev(L,L)
           if (mod(L,2).eq.mod(jcap+1,2)) then
              indev2 = indlsev(jcap+1,L)
           else
              indev2 = indlsev(jcap  ,L)
           endif
           do indev = indev1 , indev2
              trie_ls(indev,1,p_lapgz)=
     X        trie_ls(indev,1,p_gz)*snnp1ev(indev)*ga2
              trie_ls(indev,2,p_lapgz)=
     X        trie_ls(indev,2,p_gz)*snnp1ev(indev)*ga2
           end do
        end do
        do locl=1,ls_max_node
                L=ls_node(locl)
           jbasod=ls_node(locl+2*ls_dim)
           indod1 = indlsod(L+1,L)
           if (mod(L,2).eq.mod(jcap+1,2)) then
              indod2 = indlsod(jcap  ,L)
           else
              indod2 = indlsod(jcap+1,L)
           endif
           do indod = indod1 , indod2
              trio_ls(indod,1,p_lapgz)=
     X        trio_ls(indod,1,p_gz)*snnp1od(indod)*ga2
              trio_ls(indod,2,p_lapgz)=
     X        trio_ls(indod,2,p_gz)*snnp1od(indod)*ga2
           end do
        end do
!
!
!------------------------------------------------------------------
!
      ioproc=0
      trie_ls=0.0
      trio_ls=0.0
         allocate ( trieo_ls_node  ( len_trie_ls_max+len_trio_ls_max,
     &                               2, lotls ),
     &              stat=ierr )
         trieo_ls_node = 0.0
      if ( me .eq. ioproc ) then
         allocate ( trieo_ls_nodes ( len_trie_ls_max+len_trio_ls_max,
     &                               2, lotls, nodes ),
     &              stat=ierr )
         trieo_ls_nodes = 0.0
         READ(1051) trieo_ls_nodes
       else
         allocate (trieo_ls_nodes(1, 1, 1, 1), stat = ierr)
      endif
      lenrec = (len_trie_ls_max+len_trio_ls_max) * 2 * lotls

      call mpi_scatter(trieo_ls_nodes, lenrec, MPI_R_MPI_R,
     &          trieo_ls_node, lenrec, MPI_R_MPI_R,
     &          ioproc, MPI_COMM_ALL, ierr)
      DO k = 1, lotls
        DO j = 1, 2
          DO i = 1, len_trie_ls
            trie_ls(i, j, k) = trieo_ls_node(i, j, k)
          END DO
          DO i = 1, len_trio_ls
            ii = i + len_trie_ls_max
            trio_ls(i, j, k) = trieo_ls_node(ii, j, k)
          END DO
        END DO
      END DO

      call mpi_barrier(MPI_COMM_ALL,ierr)
      DEALLOCATE(trieo_ls_nodes)
      DEALLOCATE(trieo_ls_node)

!
!--------------------------------------------------------------
! fill up n+1 grid_gr in case of internal2export used.
!
        trie_ls(:,:,p_zq)=trie_ls(:,:,p_q )
        trie_ls(:,:,p_y :p_y +levs-1)=trie_ls(:,:,p_te:p_te+levs-1)
        trie_ls(:,:,p_x :p_x +levs-1)=trie_ls(:,:,p_di:p_di+levs-1)
        trie_ls(:,:,p_w :p_w +levs-1)=trie_ls(:,:,p_ze:p_ze+levs-1)
        trio_ls(:,:,p_zq)=trio_ls(:,:,p_q )
        trio_ls(:,:,p_y :p_y +levs-1)=trio_ls(:,:,p_te:p_te+levs-1)
        trio_ls(:,:,p_x :p_x +levs-1)=trio_ls(:,:,p_di:p_di+levs-1)
        trio_ls(:,:,p_w :p_w +levs-1)=trio_ls(:,:,p_ze:p_ze+levs-1)
        if ( .not. ndslfv ) then
          trie_ls(:,:,p_rt :p_rt +levh-1)=trie_ls(:,:,p_rq:p_rq+levh-1)
          trio_ls(:,:,p_rt :p_rt +levh-1)=trio_ls(:,:,p_rq:p_rq+levh-1)
        endif
!
!--------------------------------------------------------------
! fill up n+1 grid_gr in case of internal2export used.
!
        grid_gr(:,:,g_zq)=grid_gr(:,:,g_q )
        grid_gr(:,:,g_t :g_t +levs-1)=grid_gr(:,:,g_tt:g_tt+levs-1)
        grid_gr(:,:,g_u :g_u +levs-1)=grid_gr(:,:,g_uu:g_uu+levs-1)
        grid_gr(:,:,g_v :g_v +levs-1)=grid_gr(:,:,g_vv:g_vv+levs-1)
        grid_gr(:,:,g_rt:g_rt+levh-1)=grid_gr(:,:,g_rq:g_rq+levh-1)

!     
!--------------------------------------------------------
!!
      RETURN
      END
