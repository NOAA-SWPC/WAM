       SUBROUTINE TREADEO_IAU(IDATE,grid_gr,
     X               LS_NODE,LS_NODES,MAX_LS_NODES,
     X               IPRINT,
     &               global_lats_a,lats_nodes_a,lonsperlat, cfile,
     &               epse, epso, epsedn,epsodn, plnew_a, plnow_a,
     &               plnev_a, plnod_a,snnp1ev, snnp1od)
!
!-- revision history
! 
!   Sep 2012, J. Wang: change interface for nems gfs
!   Oct 2012, H. Juang: recover terrain laplacian preparation
!   May 2013  S. Moorthi - correct definition of idvt for WAM
!   Jun 2013, H. Juang: add dp and spectral moisture in case of ndsl 
!   Sep 2013, H. Juang: add computation of dp, not only give dp space
!
      use gfs_dyn_resol_def
      use gfs_dyn_layout1
      use gfs_dyn_coordinate_def					! hmhj
!      use sig_io
      use gfs_dyn_io_header
      use gfs_dyn_iau_module, only: gq_iau,gu_iau,gv_iau,gt_iau,grq_iau
      use namelist_dynamics_def
!      use gfs_dyn_namelist_def
      use gfs_dyn_vert_def
      use gfs_dyn_mpi_def
      use gfs_dyn_physcons, rerth => con_rerth, grav => con_g, 
     &    rkap => con_rocp, cpd => con_cp
      USE gfs_dyn_date_def, ONLY: FHOUR
      use gfs_dyn_tracer_config, ONLY : gfs_dyn_tracer,glbsum   ! generalized tracer
      use sigio_module
      use sigio_r_module

!
      IMPLICIT NONE
!
      INTEGER,intent(inout) :: IDATE(4)
      REAL(KIND=KIND_GRID),dimension(lonf,lats_node_a_max,gq_iau),
     &   intent(out) :: grid_gr
      INTEGER,intent(in) :: LS_NODE(LS_DIM,3)
      INTEGER,intent(in) :: LS_NODES(LS_DIM,NODES)
      INTEGER,intent(in) :: MAX_LS_NODES(NODES)
      integer,intent(in) :: lats_nodes_a(nodes)
      INTEGER,intent(inout) :: IPRINT
      integer,intent(in) :: global_lats_a(latg),lonsperlat(latg)
      character*(*),intent(in) :: cfile
      real(kind=kind_evod)   :: trie_ls(len_trie_ls,2,gq_iau)
      real(kind=kind_evod)   :: trio_ls(len_trio_ls,2,gq_iau)
      real(kind=kind_evod),dimension(len_trie_ls),intent(in) :: epse,
     &  epsedn
      real(kind=kind_evod),dimension(len_trio_ls),intent(in) :: epso,
     &  epsodn
!
      real(kind=kind_evod),dimension(len_trie_ls,latg2),intent(in) 
     &   :: plnew_a,plnev_a
      real(kind=kind_evod),dimension(len_trio_ls,latg2),intent(in) 
     &   :: plnow_a,plnod_a
!
      REAL(KIND=KIND_EVOD),intent(in) :: snnp1ev(LEN_TRIE_LS)
      REAL(KIND=KIND_EVOD),intent(in) :: snnp1od(LEN_TRIO_LS)
      integer            lats_nodes_r(nodes)
!
!------local vars
!
      REAL(KIND=KIND_EVOD)  QE(LEN_TRIE_LS,2,1)
     &,                    TEE(LEN_TRIE_LS,2,LEVS)
     &,                    DIE(LEN_TRIE_LS,2,LEVS)
     &,                    ZEE(LEN_TRIE_LS,2,LEVS)
     &,                    UUE(LEN_TRIE_LS,2,LEVS)
     &,                    VVE(LEN_TRIE_LS,2,LEVS)
     &,                    RQE(LEN_TRIE_LS,2,LEVS)
     &,                     QO(LEN_TRIO_LS,2,1)
     &,                    TEO(LEN_TRIO_LS,2,LEVS)
     &,                     DIO(LEN_TRIO_LS,2,LEVS)
     &,                    ZEO(LEN_TRIO_LS,2,LEVS)
     &,                    UUO(LEN_TRIO_LS,2,LEVS)
     &,                    VVO(LEN_TRIO_LS,2,LEVS)
     &,                    RQO(LEN_TRIO_LS,2,LEVS)
     &,                    Z00
!
      real(kind=kind_evod), allocatable  :: syn_gr_a_2(:,:,:)
!
!
!cmr  ls_node(1,1) ... ls_node(ls_max_node,1) : values of L
!cmr  ls_node(1,2) ... ls_node(ls_max_node,2) : values of jbasev
!cmr  ls_node(1,3) ... ls_node(ls_max_node,3) : values of jbasod
!
      INTEGER              J,K,L,LOCL,N,kk,ilan,itrac
      integer              i,lan,lat,lons_lat,lon,njeff,nn,lon_dim
      integer              indev
      integer              indod
      integer              indev1,indev2
      integer              indod1,indod2
      REAL(KIND=KIND_EVOD) WAVES
      REAL(KIND=KIND_EVOD), target ::  TRISCA(LNT2)
      REAL(KIND=KIND_EVOD) sikp1(levp1)
      INTEGER              INDLSEV,JBASEV
      INTEGER              INDLSOD,JBASOD
c$$$      REAL(KIND=KIND_IO4) Z(lnt2)
!
      type(sigio_head) head
      type(sigio_dbti) dati
!
      integer              iret, num_dta,nft,lon1,lon2
      real(kind=kind_evod) psurfff
      real(kind=kind_evod) pressk, tem, rkapi, rkapp1
!
      real(kind=kind_evod) teref(levp1),ck5p(levp1)			! hmhj
 
 
      INCLUDE 'function2'
!!
!!
      print *,' enter treadeo.io_fd '					! hmhj
      if (semilag) then
         lon1=lon_dim_a
         lon2=lonf
      else
         lon1=lonfx
         lon2=lonfx
      endif
      allocate(syn_gr_a_2(lon2,lats_node_a,levs))

      nft=152
      call sigio_rropen(nft,cfile,iret)
      call sigio_alhead(head,iret)
      call sigio_rrhead(nft,head,iret)
!
      ivsinp = head%ivs
      IF (me .eq. 0) THEN
        print *,' In treadeo iret=',iret,' cfile=',cfile
     &,' ivs=',head%ivs,' levs=',head%levs
      ENDIF
      IF (iret .ne. 0) THEN
        print *,' unable to read from unit ',nft,' Job Aborted'
     &,' iret=',iret,' me=',me
        call mpi_quit(7777)
      ENDIF
!
      idvc = head%idvc  !idvc=3:sigma-theta and/or p, idvc=2:sigma-p, idvc=1:sigma files
      idsl = head%idsl
!
      dati%i = 2                               ! Surface pressure
      dati%f => TRISCA                       
      call sigio_rrdbti(nft,head,dati,iret)
!     IF (me == 0) print *,' SFCPRES=',trisca(1:10)
      CALL TRISEORI(TRISCA,QE(1,1,1),QO(1,1,1),1,LS_NODE)
      DO K=1,LEVS
        dati%i = k + 2                        ! Virtual Temperature or CpT
        dati%f => TRISCA                       

        call sigio_rrdbti(nft,head,dati,iret)

        CALL TRISEORI(TRISCA,TEE(1,1,K),TEO(1,1,K),1,LS_NODE)
      ENDDO
!
      DO K=1,LEVS
         dati%i = levs + 2 + (k-1) * 2 + 1     ! Divergence
         dati%f => TRISCA
         call sigio_rrdbti(nft,head,dati,iret)
         CALL TRISEORI(TRISCA,DIE(1,1,K),DIO(1,1,K),1,LS_NODE)
!
         dati%i = levs + 2 + (k-1) * 2 + 2     ! Vorticity
         dati%f => TRISCA
         call sigio_rrdbti(nft,head,dati,iret)
         CALL TRISEORI(TRISCA,ZEE(1,1,K),ZEO(1,1,K),1,LS_NODE)
      END DO
csela print*,' levh=',levh
!
!
!
      DO itrac=1,ntrac
        RQE=0.
        RQO=0.
        kk = 0
        IF (itrac == 1) THEN
          kk = 1
        ELSEIF (itrac == ntoz) THEN
          kk = ntoz
        ELSEIF (itrac >= ntcw .and. itrac < ntcw+ncld-1) THEN
          DO n=1,ncld
            IF (itrac == ntcw+n-1) kk = ntcw+n-1
          ENDDO
        ELSE
          kk = itrac
        ENDIF
!
         DO k=1,levs
           dati%i = levs * (2+itrac) + 2 + k             ! Tracers starting with q
           dati%f => TRISCA
           call sigio_rrdbti(nft,head,dati,iret)
           CALL TRISEORI(TRISCA,RQE(1,1,k),RQO(1,1,k),1,LS_NODE)
         END DO
         call spect_to_grid_iau(RQE,RQO,
     &         syn_gr_a_2,
     &         ls_node,ls_nodes,max_ls_nodes,
     &         lats_nodes_a,global_lats_a,lonsperlat,
     &         plnev_a,plnod_a,levs,lon1,lon2)
!jw!$omp parallel do private(lan)
!jw!$omp+private(lat,lon_dim,lons_lat,i,itrac)
      DO lan=1,lats_node_a
        lat = global_lats_a(ipt_lats_node_a-1+lan)
        lon_dim = lon_dims_a(lan)
        lons_lat = lonsperlat(lat)
        DO k=1,levs
          DO i=1,lons_lat
            grid_gr(i,lan,grq_iau+(kk-1)*levs+k-1)=syn_gr_a_2(i,lan,k)  ! make sure I have them in the right slots
          ENDDO
        ENDDO
      ENDDO

      END DO
 
!
      DO k=1,levs
         call dezouv(DIE(1,1,k), ZEO(1,1,k),
     x               UUE(1,1,k), VVO(1,1,k),
     x               epsedn,epsodn,snnp1ev,snnp1od,ls_node)
!
         call dozeuv(DIO(1,1,k), ZEE(1,1,k),
     x               UUO(1,1,k), VVE(1,1,k),
     x               epsedn,epsodn,snnp1ev,snnp1od,ls_node)
      ENDDO
      call spect_to_grid_iau(QE,QO,
     &         syn_gr_a_2(:,:,1),
     &         ls_node,ls_nodes,max_ls_nodes,
     &         lats_nodes_a,global_lats_a,lonsperlat,
     &         plnev_a,plnod_a,1,lon1,lon2)
!jw!$omp parallel do private(lan)
!jw!$omp+private(lat,lon_dim,lons_lat,i)
      DO lan=1,lats_node_a
        lat = global_lats_a(ipt_lats_node_a-1+lan)
        lon_dim = lon_dims_a(lan)
        lons_lat = lonsperlat(lat)
        DO i=1,lons_lat
          grid_gr(i,lan,gq_iau)=syn_gr_a_2(i,lan,1)
        ENDDO
      ENDDO

      call spect_to_grid_iau(UUE,UUO,
     &         syn_gr_a_2,
     &         ls_node,ls_nodes,max_ls_nodes,
     &         lats_nodes_a,global_lats_a,lonsperlat,
     &         plnev_a,plnod_a,levs,lon1,lon2)
!jw!$omp parallel do private(lan)
!jw!$omp+private(lat,lon_dim,lons_lat,i,k)
      DO lan=1,lats_node_a
        lat = global_lats_a(ipt_lats_node_a-1+lan)
        lon_dim = lon_dims_a(lan)
        lons_lat = lonsperlat(lat)
        DO k=1,levs
          DO i=1,lons_lat
            grid_gr(i,lan,gu_iau+k-1)=syn_gr_a_2(i,lan,k)
          ENDDO
        ENDDO
      ENDDO

      call spect_to_grid_iau(VVE,VVO,
     &         syn_gr_a_2,
     &         ls_node,ls_nodes,max_ls_nodes,
     &         lats_nodes_a,global_lats_a,lonsperlat,
     &         plnev_a,plnod_a,levs,lon1,lon2)
!jw!$omp parallel do private(lan)
!jw!$omp+private(lat,lon_dim,lons_lat,i,k)
      DO lan=1,lats_node_a
        lat = global_lats_a(ipt_lats_node_a-1+lan)
        lon_dim = lon_dims_a(lan)
        lons_lat = lonsperlat(lat)
        DO k=1,levs
          DO i=1,lons_lat
            grid_gr(i,lan,gv_iau+k-1)=syn_gr_a_2(i,lan,k)
          ENDDO
        ENDDO
      ENDDO

      call spect_to_grid_iau(TEE,TEO,
     &         syn_gr_a_2,
     &         ls_node,ls_nodes,max_ls_nodes,
     &         lats_nodes_a,global_lats_a,lonsperlat,
     &         plnev_a,plnod_a,levs,lon1,lon2)
!jw!$omp parallel do private(lan)
!jw!$omp+private(lat,lon_dim,lons_lat,i,k)
      DO lan=1,lats_node_a
        lat = global_lats_a(ipt_lats_node_a-1+lan)
        lon_dim = lon_dims_a(lan)
        lons_lat = lonsperlat(lat)
        DO k=1,levs
          DO i=1,lons_lat
            grid_gr(i,lan,gt_iau+k-1)=syn_gr_a_2(i,lan,k)
          ENDDO
        ENDDO
      ENDDO

!      IF(me==0) THEN
!        open (991,file='ps_grid',form='unformatted') 
!        write(991)real(grid_gr(1:lonf,1:lats_node_a,gq_iau),kind=4)
!        close(991)
!      ENDIF
!      call mpi_barrier(MPI_COMM_ALL,i)
!      call mpi_abort(mpi_comm_all,j,i)
!      print *,'in treadeo.io, ps=',
!     &   maxval(grid_gr(:,1:lats_node_a,gq_iau)),
!     &   minval(grid_gr(:,1:lats_node_a,gq_iau))
!      print*,'in treadeo.io, t=',
!     &   maxval(grid_gr(:,1:lats_node_a,gt_iau:gt_iau+levs-1)),
!     &   minval(grid_gr(:,1:lats_node_a,gt_iau:gt_iau+levs-1))
!      print*,'in treadeo.io, tracer=',
!     &   maxval(grid_gr(:,1:lats_node_a,grq_iau:grq_iau+levh-1)),
!     &   minval(grid_gr(:,1:lats_node_a,grq_iau:grq_iau+levh-1))
!      print*,'in treadeo.io, u=',
!     &   maxval(grid_gr(:,1:lats_node_a,gu_iau:gu_iau+levs-1)),
!     &   minval(grid_gr(:,1:lats_node_a,gu_iau:gu_iau+levs-1))
!      print*,'in treadeo.io, v=',
!     &   maxval(grid_gr(:,1:lats_node_a,gv_iau:gv_iau+levs-1)),
!     &   minval(grid_gr(:,1:lats_node_a,gv_iau:gv_iau+levs-1))
!  
!      print *,' leave treadeo.io_fd '

      RETURN
      END
