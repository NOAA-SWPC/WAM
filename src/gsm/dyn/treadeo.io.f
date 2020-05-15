       SUBROUTINE TREADEO(IDATE,trie_ls,trio_ls,
     &               grid_gr,
     X               LS_NODE,LS_NODES,MAX_LS_NODES,
     X               PDRYINI,IPRINT,
     &               global_lats_a,lats_nodes_a,lonsperlat, cfile,
     &               epse, epso, epsedn,epsodn, plnew_a, plnow_a,
     &               plnev_a, plnod_a,snnp1ev, snnp1od)
!
!      SUBROUTINE TREADEO(NFT,FHOUR,IDATE,
!     X                   GZE,QE,TEE,DIE,ZEE,RQE,
!     X                   GZO,QO,TEO,DIO,ZEO,RQO,
!     X                   LS_NODE,LS_NODES,MAX_LS_NODES,
!     &                   plnev_r,plnod_r,plnew_r,plnow_r,lats_nodes_r,
!     X                   SNNP1EV,SNNP1OD,pdryini,IPRINT,
!     &                   phy_f3d, phy_f2d, global_lats_a,
!     &                   lonsperlar,cfile)
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
      real(kind=kind_evod),intent(out) :: trie_ls(len_trie_ls,2,lotls)
      real(kind=kind_evod),intent(out) :: trio_ls(len_trio_ls,2,lotls)
      REAL(KIND=KIND_GRID),dimension(lonf,lats_node_a_max,lotgr),
     &   intent(out) :: grid_gr
      INTEGER,intent(in) :: LS_NODE(LS_DIM,3)
      INTEGER,intent(in) :: LS_NODES(LS_DIM,NODES)
      INTEGER,intent(in) :: MAX_LS_NODES(NODES)
      integer,intent(in) :: lats_nodes_a(nodes)
      INTEGER,intent(inout) :: IPRINT
      integer,intent(in) :: global_lats_a(latg),lonsperlat(latg)
      character*(*),intent(in) :: cfile
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
      INTEGER              NTRACI, ntozi, ntcwi, ncldi
     &,                    nt0, direction,levsi
!
      REAL(KIND=KIND_EVOD) GZE(LEN_TRIE_LS,2)
     &,                 LAPGZE(LEN_TRIE_LS,2)
     &,                     QE(LEN_TRIE_LS,2)
     &,                    DPE(LEN_TRIE_LS,2,LEVS)
     &,                    TEE(LEN_TRIE_LS,2,LEVS)
     &,                    DIE(LEN_TRIE_LS,2,LEVS)
     &,                    ZEE(LEN_TRIE_LS,2,LEVS)
     &,                    RQE(LEN_TRIE_LS,2,LEVS,ntrac)
     &,                    GZO(LEN_TRIO_LS,2)
     &,                 LAPGZO(LEN_TRIO_LS,2)
     &,                     QO(LEN_TRIO_LS,2)
     &,                    DPO(LEN_TRIO_LS,2,LEVS)
     &,                    TEO(LEN_TRIO_LS,2,LEVS)
     &,                     DIO(LEN_TRIO_LS,2,LEVS)
     &,                    ZEO(LEN_TRIO_LS,2,LEVS)
     &,                    RQO(LEN_TRIO_LS,2,LEVS,ntrac)
!    &,                    SL(LEVS)
!    &,                    SI(LEVP1)
     &,                    Z00
!
      real(kind=kind_evod) syn_gr_a_1(lonfx*(lots+1),lats_dim_a)
      real(kind=kind_evod) syn_gr_a_2(lonfx*(lots+1),lats_dim_a)
      real(kind=kind_evod)  rq_gr_a_1(lonfx*levh,lats_dim_a)
      real(kind=kind_evod)  rq_gr_a_2(lonfx*levh,lats_dim_a)
      real(kind=kind_evod)  dp_gr_a_1(lonfx*levs,lats_dim_a)
      real(kind=kind_evod)  dp_gr_a_2(lonfx*levs,lats_dim_a)
      real(kind=kind_evod)   gtv (ngptc,levs), gq  (ngptc)
      real(kind=kind_evod)   prsl(ngptc,levs), dprs(ngptc,levs)
!
!
!cmr  ls_node(1,1) ... ls_node(ls_max_node,1) : values of L
!cmr  ls_node(1,2) ... ls_node(ls_max_node,2) : values of jbasev
!cmr  ls_node(1,3) ... ls_node(ls_max_node,3) : values of jbasod
!
      INTEGER              J,K,L,LOCL,N,lv,kk,ilan
      integer              i,lan,lat,lons_lat,lon,njeff,nn,lon_dim
      integer              indev
      integer              indod
      integer              indev1,indev2
      integer              indod1,indod2
      REAL(KIND=KIND_EVOD) GA2,GENCODE,GZBAR,ORDER,REALFORM
      REAL(KIND=KIND_EVOD) TRUN,WAVES,XLAYERS
      REAL(KIND=KIND_EVOD) XI(LEVP1),XL(LEVS)
      REAL(KIND=KIND_EVOD), target ::  TRISCA(LNT2)
      REAL(KIND=KIND_EVOD) sikp1(levp1)
      REAL(KIND=KIND_IO4)   VTID,RUNID4,USRID4,pdryini4,XNCLD,xgf
      REAL(KIND=KIND_IO8)   PDRYINI
      INTEGER              INDLSEV,JBASEV
      INTEGER              INDLSOD,JBASOD
c$$$      REAL(KIND=KIND_IO4) Z(lnt2)
!
      type(sigio_head) head
      type(sigio_dbti) dati
!
      integer              iret, num_dta,nft
      real(kind=kind_evod) psurfff
      real(kind=kind_evod) pressk, tem, rkapi, rkapp1
!
      real(kind=kind_evod) teref(levp1),ck5p(levp1)			! hmhj
 
 
      INCLUDE 'function2'
!!
!!
      nft=152
      call sigio_rropen(nft,cfile,iret)
      call sigio_alhead(head,iret)
      call sigio_rrhead(nft,head,iret)
!
      ivsinp = head%ivs
      if (me .eq. 0) then
        print *,' In treadeo iret=',iret,' cfile=',cfile
     &,' ivs=',head%ivs,' levs=',head%levs
      endif
      if (iret .ne. 0) then
        print *,' unable to read from unit ',nft,' Job Aborted'
     &,' iret=',iret,' me=',me
        call mpi_quit(7777)
      endif
!
      idvc = head%idvc  !idvc=3:sigma-theta and/or p, idvc=2:sigma-p, idvc=1:sigma files
      idsl = head%idsl
!     rewind(nft)
!     READ(NFT)

      RKAPI  = 1.0 / RKAP
      RKAPP1 = 1.0 + RKAP

      if (me .eq. 0) then
        print *,' gen_coord_hybrid=',gen_coord_hybrid,
     &' idvc=',head%idvc,' idvm=',head%idvm,' hybrid=',hybrid
      endif
      if (gen_coord_hybrid) then					! hmhj

        sfcpress_id = mod ( head%idvm , 10 )				! hmhj
        thermodyn_id = mod ( head%idvm / 10 , 10 )			! hmhj
        if (me .eq. 0) then
          print *,' sfcpress_id thermodyn_id ',sfcpress_id,thermodyn_id	! hmhj
        endif
!   ak bk ck in file have the same order as model			! hmhj
        do k=1,levp1							! hmhj
          ak5(k) = head%vcoord(k,1)/1000.				! hmhj
          bk5(k) = head%vcoord(k,2)					! hmhj
          ck5(k) = head%vcoord(k,3)/1000.				! hmhj
        enddo								! hmhj
        vertcoord_id=0							! hmhj
        do k=1,levp1							! hmhj
          if( ck5(k).ne.0.0 ) vertcoord_id=3				! hmhj
        enddo								! hmhj
! provide better estimated press 					! hmhj
        psurfff = 101.3							! hmhj
        if( thermodyn_id.eq.3 ) then					! hmhj
         do k=1,levp1							! hmhj
          thref(k)=300.*cpd						! hmhj
          teref(k)=255.*cpd						! hmhj
         enddo								! hmhj
        else								! hmhj
         do k=1,levp1							! hmhj
          thref(k)=300.							! hmhj
          teref(k)=255.							! hmhj
         enddo								! hmhj
        endif								! hmhj
        ck5p(levp1)=ck5(levp1)						! hmhj
! modify.Weiyu.
!--------------
!        do k=1,levp1							! hmhj
        do k=1,levp1 - 1							! hmhj
          ck5p(k)=ck5(k)*(teref(k)/thref(k))**rkapi			! hmhj
        enddo								! hmhj
        if( me.eq.0 ) then						! hmhj
          do k=1,levp1							! hmhj
            pressk=ak5(k)+bk5(k)*psurfff+ck5p(k)			! hmhj
            print 180,k,ak5(k),bk5(k),ck5(k),pressk			! hmhj
180         format('k=',i3,'  ak5=',f13.9,'  bk5=',e15.8,		! hmhj
     &            '   ck5=',f13.9,'  closed pressk=',f10.6)		! hmhj
          enddo								! hmhj
        endif								! hmhj
        do k=1,levs
          dbk(k) = bk5(k)-bk5(k+1)
          bkl(k) = (bk5(k)+bk5(k+1))*0.5
          ck(k)  = ak5(k)*bk5(k+1)-ak5(k+1)*bk5(k)
        enddo
        do k=1,levp1							! hmhj
          si(k)=ak5(k)/psurfff+bk5(k)+ck5p(k)/psurfff			! hmhj
        enddo								! hmhj
        do k=1,levs							! hmhj
          sl(k)=0.5*(si(k)+si(k+1))					! hmhj
        enddo								! hmhj

      else if (hybrid .and. idvc .eq. 2) then				! hmhj
!       idsl=slid  !=2,pk=0.5*(p(k+1/2)+p(k-1/2)) check alfa(1)  am_bm
!   ak bk order in "sigma" file is bottom to top !!!!!!!!!!!!!!!!!!
        psurfff = 101.3
        do k=1,levp1
          ak5(k) = head%vcoord(levp1+1-k,1)/1000.
          bk5(k) = head%vcoord(levp1+1-k,2)
          pressk = ak5(k)+bk5(k)*psurfff
          
          if(me.eq.0)print 190,k,ak5(k),bk5(k),pressk
190       format('k=',i3,'  ak5=',E15.8,'  bk5=',e15.8,
     &           '  pressk=',E14.6)
           
        enddo
        do k=1,levs
          dbk(k) = bk5(k+1)-bk5(k)
          bkl(k) = (bk5(k+1)+bk5(k))*0.5
          ck(k)  = ak5(k+1)*bk5(k)-ak5(k)*bk5(k+1)
          if(me.eq.0)print 200,k,dbk(k),ck(k)
200       format('k=',i3,'  dbk=',f9.6,'  ck=',e13.5)
        enddo
!
! hmhj give an estimated si and sl for dynamics
        do k=1,levs+1
          si(levs+2-k)=ak5(k)/psurfff+bk5(k) !ak(k) bk(k) go top to bottom
        enddo
        do k=1,levs
          sl(k)=0.5*(si(k)+si(k+1))
        enddo
!
!     elseif (head%idvc .eq. 1) then
      elseif (head%idvc .le. 1) then
        si(:)    = head%vcoord(:,1)
        sik(:)   = si(:) ** rkap
        sikp1(:) = si(:) ** rkapp1
        do k=1,levs
          tem      = rkapp1 * (si(k) - si(k+1))
          slk(k)   = (sikp1(k)-sikp1(k+1))/tem
          sl(k)    = slk(k) ** rkapi
!         sl(k)    = ((sikp1(k)-sikp1(k+1))/tem)**rkapi
          if (me .eq. 0) print 250, k, si(k), sl(k)
250       format('k=',i2,'  si=',f9.6,'  sl=',e13.5)
        enddo
      else
        print *,' Non compatible Initial state IDVC=',head%idvc
     &,' iret=',iret
        call MPI_QUIT(333)
      endif
!
csela print*,' read second record successfully '
      FHOUR       = head%fhour
      idate       = head%idate
      WAVES       = head%jcap
      lonr        = head%lonr
      latr        = head%latr
!      XLAYERS     = head%levs
      levsi       = head%levs
      itrun       = head%itrun
      iorder       = head%iorder
      irealf       = head%irealf
!      REALFORM    = head%irealf
      icen        = 7
      icen2       = head%icen2
      igen        = head%igen
      idpp        = head%idpp
      idrun       = head%idrun
      itrun       = head%itrun
      idusr       = head%idusr
      ienst       = head%iens(1)
      iensi       = head%iens(2)
      if (fhour .gt. 0.0 .and. head%nxss .eq. 0 .and.
     &    head%pdryini > 0.0 ) then
        if (pdryini .eq. 0.0) pdryini = head%pdryini
      endif
      if (me == 0) print *,' IN TREAD PDRYINI=',pdryini,
     &                     ' head=',head%pdryini
!sela ntraci = nint(tracers-1)
      ntraci = head%ntrac
      if (head%idvt .gt. 0.0) then
        if (head%idvt == 200) then
          ntozi = 2
          ntcwi = 3
          ncldi = 1
        else
          nt0   = mod(head%idvt,10)
          if (nt0 > 0) then
            ntcwi = head%idvt / 10
            ntozi = head%idvt - ntcwi * 10 + 1
            ntcwi = ntcwi + 1
          else
            ntcwi = ntcw
            ntozi = ntoz
          endif
          ncldi = head%ncldt
        endif
      elseif(ntraci .eq. 2) then
        ntozi = 2
        ntcwi = 0
        ncldi = 0
      elseif(ntraci .eq. 3) then
        ntozi = 2
        ntcwi = 3
        ncldi = 1
      else
        ntozi = 0
        ntcwi = 0
        ncldi = 0
      endif
      ixgr = head%ixgr
!
      if (ntrac <= 3) then
        idvt = (ntcw-1)*10 + ntoz - 1
      else
        idvt = head%idvt
      endif
!
      IF (me.eq.0) THEN
        write(*,*)'nfile,in treadeo fhour,idate=',nft,fhour,idate
     &, ' ntozi=',ntozi,' ntcwi=',ntcwi,' ncldi=',ncldi
     &, ' ntraci=',ntraci,' tracers=',head%ntrac,' vtid=',head%idvt
     &,  head%ncldt,' idvc=',head%idvc,' jcap=',head%jcap
     &, ' ixgr=',ixgr,' pdryini=',pdryini
      ENDIF
!cjfe
      IF(IPRINT.EQ.1)
     X PRINT *,'TREAD UNIT,FHOUR,IDATE=',NFT,FHOUR,IDATE
 
!
      dati%i = 1                                           ! hs
      dati%f => TRISCA
      call sigio_rrdbti(nft,head,dati,iret)
!     if (me == 0) print *,' Z_R=',trisca(1:10),' iret=',iret
      Z   = TRISCA
      Z_R = TRISCA
      CALL TRISEORI(TRISCA,GZE,GZO,1,LS_NODE)
      Z00=TRISCA(1)
      GA2=GRAV/(RERTH*RERTH)
! hmhj laplacian terrain for divergence
      DO LOCL=1,LS_MAX_NODE
              L=ls_node(locl,1)
         jbasev=ls_node(locl,2)
         indev1 = indlsev(L,L)
         if (mod(L,2).eq.mod(jcap+1,2)) then
            indev2 = indlsev(jcap+1,L)
         else
            indev2 = indlsev(jcap  ,L)
         endif
         do indev = indev1 , indev2
            LAPGZE(INDEV,1)=
     X         GZE(INDEV,1)*SNNP1EV(INDEV)*GA2
            LAPGZE(INDEV,2)=
     X         GZE(INDEV,2)*SNNP1EV(INDEV)*GA2
         END DO
      END DO
!     if (me == 0) print *,' GZE=',maxval(GZe(:,1:2)),
!    &  minval(GZe(:,1:2))
      DO LOCL=1,LS_MAX_NODE
              L=ls_node(locl,1)
         jbasod=ls_node(locl,3)
         indod1 = indlsod(L+1,L)
         if (mod(L,2).eq.mod(jcap+1,2)) then
            indod2 = indlsod(jcap  ,L)
         else
            indod2 = indlsod(jcap+1,L)
         endif
         do indod = indod1 , indod2
            LAPGZO(INDOD,1)=
     X         GZO(INDOD,1)*SNNP1OD(INDOD)*GA2
            LAPGZO(INDOD,2)=
     X         GZO(INDOD,2)*SNNP1OD(INDOD)*GA2
         END DO
      END DO
!     if (me == 0) print *,' GZE=',maxval(GZo(:,1:2)),
!    &  minval(GZo(:,1:2))
 
      if (mod(head%idvm/10,10) == 3 .and. me == 0)then
        print *,' CPI=',head%cpi(1:ntraci+1)
        print *,' RI=',head%ri(1:ntraci+1)
      endif
      dati%i = 2                               ! Surface pressure
      dati%f => TRISCA                       
      call sigio_rrdbti(nft,head,dati,iret)
!     if (me == 0) print *,' SFCPRES=',trisca(1:10)
      CALL TRISEORI(TRISCA,QE,QO,1,LS_NODE)
      DO K=1,LEVS
        dati%i = k + 2                        ! Virtual Temperature or CpT
        dati%f => TRISCA                       

        call sigio_rrdbti(nft,head,dati,iret)

        CALL TRISEORI(TRISCA,TEE(1,1,K),TEO(1,1,K),1,LS_NODE)
      enddo
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
      RQE=0.
      RQO=0.
      DO K=1,ntraci
        kk = 0
        if (k .eq. 1) then
          kk = 1
        elseif (k .eq. ntozi) then
          kk = ntoz
        elseif (k .ge. ntcwi .and. k .lt. ntcwi+ncldi-1) then
          do n=1,ncldi
            if (k .eq. ntcwi+n-1) kk = ntcw+n-1
          enddo
        else
          kk = k
        endif
!
        DO lv=1,levs
          dati%i = levs * (2+k) + 2 + lv             ! Tracers starting with q
          dati%f => TRISCA
          call sigio_rrdbti(nft,head,dati,iret)
          CALL TRISEORI(TRISCA,RQE(1,1,lv,KK),RQO(1,1,lv,KK),1,LS_NODE)
        END DO
      END DO
!
!jw      if (((ixgr .eq. 4 .and. num_p3d .eq. 4) .or.   ! Zhao Scheme!
!jw     &     (ixgr .eq. 5 .and. num_p3d .eq. 3))       ! Ferrier Scheme!
!jw     &     .and. fhour .gt. 0.1) then
!jw          kmsk(:,:)  = 0
!jw          num_dta = (ntraci+3)*levs + 2
!jw          do nn=1,num_p3d
!jw            do k=1,levs
!jw              dati%i = num_dta + (nn-1)*levs + k      ! physics 3D grid fields
!jw              dati%f => buff1
!jw              call sigio_rrdbti(nft,head,dati,iret)
!jw              call split2d_r(buff1(1),buffo,global_lats_a)
!jw              CALL interpred(1,kmsk,buffo,buff2,global_lats_a,
!jw     &                lonsperlar)
!
!jw              do lan=1,LATS_NODE_R
!jw                lat = global_lats_a(ipt_lats_node_a-1+lan)
!jw                lons_lat = lonsperlar(lat)
!jw                DO i=1,lons_lat
!jw                  phy_f3d(i,k,nn,lan) = buff2(i,lan)
!jw                enddo
!jw              enddo
!jw            enddo
!jw          enddo
!jw          do nn=1,num_p2d
!jw            dati%i = num_dta + num_p3d*levs + nn      ! physics 2D grid fields
!jw            dati%f => buff1
!jw            call sigio_rrdbti(nft,head,dati,iret)
!jw            call split2d_r(buff1,buffo,global_lats_a)
!jw            CALL interpred(1,kmsk,buffo,buff2,global_lats_a,
!jw     &              lonsperlar)
!jw            do j=1,lats_node_r
!jw              do i=1,lonr
!jw                phy_f2d(i,nn,j) = buff2(i,j)
!jw              enddo
!jw            enddo
!jw          enddo
!       endif
!jw      endif
!jw      if (head%nxss .gt. 0) then
!jw        dati%i = num_dta + num_p3d*levs + num_p2d + 1    ! pdryini
!jw        dati%f => buff1
!jw        call sigio_rrdbti(nft,head,dati,iret)
!jw        pdryini = buff1(1)
!jw      endif
!
      iprint=0
 
!sela IF(IPRINT.NE.1) RETURN
!     DO K=1,LEVS
!sela    XL(K)=XL(K)-SL(K)
!        sL(K)=XL(K)
!     END DO
!Moor PRINT 100,(XL(K),K=1,LEVS)
!     DO K=1,LEVP1
!sela    XI(K)=XI(K)-SI(K)
!        sI(K)=XI(K)
!     END DO
!Moor PRINT 100,(XI(K),K=1,LEVP1)
100   FORMAT(1H0, 12   (E10.3))
!Moor PRINT 101,NFT,FHOUR,IDATE,Z00
101   FORMAT (1H0, 'IF ABOVE TWO ROWS NOT ZERO,INCONSISTENCY IN SIG.DEF'
     1,'ON NFT=',I2,2X,F6.1,2X,4(I4),'Z00=',E12.4)
!!!!
!
!   Convert from virtual temperature to enthalpy if need
!
!jw      if( thermodyn_id.le.1 .and. sfcpress_id.le.1
!jw     &    .and. gen_coord_hybrid ) then

!
!jw        if (.NOT.LIOPE.or.icolor.ne.2) then
!
!jw          direction=1	! from (tv,lnps) to (enthalpy,ps)
!jw          call spect_tv_enthalpy_ps
!jw     &       (direction,
!jw     &        QE,QO,TEE,TEO,RQE,RQO,
!jw     &        ls_node,ls_nodes,max_ls_nodes,
!jw     &        lats_nodes_a,global_lats_a,lonsperlat,
!jw     &        plnev_a,plnod_a,plnew_a,plnow_a)

!jw        endif	! .NOT.LIOPE.or.icolor.ne.2
!
!jw        if( thermodyn_id==3 ) then
!jw          do k=1,levp1
!jw            thref(k)=300.*cpd	
!jw            teref(k)=255.*cpd
!jw          enddo		
!jw          thermodyn_id = 3
!jw          sfcpress_id  = 2
!jw        else
!jw          thermodyn_id = 1
!jw          sfcpress_id  = 2
!jw        endif

!jw      endif
!
!set to spect coef

      DO i=1,len_trie_ls
        trie_ls(i,1,P_LAPGZ)=LAPGZE(i,1)
        trie_ls(i,2,P_LAPGZ)=LAPGZE(i,2)
        trie_ls(i,1,P_GZ)=GZE(i,1)
        trie_ls(i,2,P_GZ)=GZE(i,2)
        trie_ls(i,1,P_Q)=QE(i,1)
        trie_ls(i,2,P_Q)=QE(i,2)
      ENDDO
      DO i=1,len_trio_ls
        trio_ls(i,1,P_LAPGZ)=LAPGZO(i,1)
        trio_ls(i,2,P_LAPGZ)=LAPGZO(i,2)
        trio_ls(i,1,P_GZ)=GZO(i,1)
        trio_ls(i,2,P_GZ)=GZO(i,2)
        trio_ls(i,1,P_Q)=QO(i,1)
        trio_ls(i,2,P_Q)=QO(i,2)
      ENDDO

! 
      DO K=1,LEVS
!
        DO i=1,len_trie_ls
          trie_ls(i,1,P_TE+k-1)=TEE(i,1,k)
          trie_ls(i,2,P_TE+k-1)=TEE(i,2,k)
          trie_ls(i,1,P_DI+k-1)=DIE(i,1,k)
          trie_ls(i,2,P_DI+k-1)=DIE(i,2,k)
          trie_ls(i,1,P_ZE+k-1)=ZEE(i,1,k)
          trie_ls(i,2,P_ZE+k-1)=ZEE(i,2,k)
        ENDDO
        DO i=1,len_trio_ls
          trio_ls(i,1,P_TE+k-1)=TEO(i,1,k)
          trio_ls(i,2,P_TE+k-1)=TEO(i,2,k)
          trio_ls(i,1,P_DI+k-1)=DIO(i,1,k)
          trio_ls(i,2,P_DI+k-1)=DIO(i,2,k)
          trio_ls(i,1,P_ZE+k-1)=ZEO(i,1,k)
          trio_ls(i,2,P_ZE+k-1)=ZEO(i,2,k)
        ENDDO

      ENDDO
!
      if( .not. ndslfv ) then
!     print *,'in treadeo, readin rh,levh=',levh,'ntraci=',
!    &  ntraci,'levs=',levs,'len_trio_ls=',len_trio_ls
      DO K=1,ntraci
        DO L=1,levs

          kk=(k-1)*levs+L-1
          DO i=1,len_trie_ls
            trie_ls(i,1,P_RQ+kk)=RQE(i,1,l,k)
            trie_ls(i,2,P_RQ+kk)=RQE(i,2,l,k)
          ENDDO
          DO i=1,len_trio_ls
            trio_ls(i,1,P_RQ+kk)=RQO(i,1,l,k)
            trio_ls(i,2,P_RQ+kk)=RQO(i,2,l,k)
          ENDDO
!
        ENDDO
      ENDDO
      endif

!set trie/o at n+1
      do k=1,levs
        trie_ls(:,:,p_w+k-1)=trie_ls(:,:,p_ze+k-1)
        trie_ls(:,:,p_x+k-1)=trie_ls(:,:,p_di+k-1)
        trie_ls(:,:,p_y+k-1)=trie_ls(:,:,p_te+k-1)
        trio_ls(:,:,p_w+k-1)=trio_ls(:,:,p_ze+k-1)
        trio_ls(:,:,p_x+k-1)=trio_ls(:,:,p_di+k-1)
        trio_ls(:,:,p_y+k-1)=trio_ls(:,:,p_te+k-1)
      enddo
      trie_ls(:,:,p_zq) = trie_ls(:,:,p_q)
      trio_ls(:,:,p_zq) = trio_ls(:,:,p_q)
      if( .not. ndslfv ) then
        do k=1,levh
          trie_ls(:,:,p_rt+k-1)=trie_ls(:,:,p_rq+k-1)
          trio_ls(:,:,p_rt+k-1)=trio_ls(:,:,p_rq+k-1)
        enddo
      endif
!
!     if(me==0) then
!       open (991,file='di_in',form='unformatted') 
!       write(991)trie_ls(1:len_trie_ls,1,P_di)
!       write(991)trie_ls(1:len_trie_ls,2,P_di)
!       write(991)trio_ls(1:len_trio_ls,1,P_di)
!       write(991)trio_ls(1:len_trio_ls,2,P_di)
!       close(991)
!       print *,'in treadeo,len_trie_ls=',len_trie_ls,
!    &   'len_trio_ls=',len_trio_ls
!     endif
!
! transfer to grid space
!
      call spect_to_grid_gz
     &        (trie_ls,trio_ls,
     &         syn_gr_a_1,syn_gr_a_2,
     &         ls_node,ls_nodes,max_ls_nodes,
     &         lats_nodes_a,global_lats_a,lonsperlat,
     &         epse,epso,epsedn,epsodn,
     &         snnp1ev,snnp1od,plnev_a,plnod_a)
!
!set data from syn back to grid_gr
!
!jw!$omp parallel do private(lan)
!jw!$omp+private(lat,lon_dim,lons_lat,i,k)
!
      do lan=1,lats_node_a
        lat = global_lats_a(ipt_lats_node_a-1+lan)
        lon_dim = lon_dims_a(lan)
        lons_lat = lonsperlat(lat)
        do k=1,levs
          do i=1,lons_lat
            grid_gr(i,lan,G_uu+k-1)=
     &         syn_gr_a_2(i+(ksu -2+k)*lon_dim,lan)
            grid_gr(i,lan,G_vv+k-1)=
     &         syn_gr_a_2(i+(ksv -2+k)*lon_dim,lan)
            grid_gr(i,lan,G_tt+k-1)=
     &         syn_gr_a_2(i+(kst -2+k)*lon_dim,lan)
          enddo
        enddo
        do i=1,lons_lat
            grid_gr(i,lan,G_q)= 
     &         syn_gr_a_2(i+(ksq-1)*lon_dim,lan)
            grid_gr(i,lan,G_gz)= 
     &         syn_gr_a_2(i+kzsphi*lon_dim,lan)
!      print *,'in treadeo.io, lan',lan,'ps=',
!     &  maxval(grid_gr(1:lons_lat,lan,g_qm)),
!     &  minval(grid_gr(1:lons_lat,lan,g_qm))
        enddo
!
        if( .not. ndslfv ) then
        do k=1,levh
          do i=1,lons_lat
            grid_gr(i,lan,G_rq+k-1)= 
     &         syn_gr_a_2(i+(ksr-2+k)*lon_dim,lan)
          enddo
        enddo
        endif
!
      enddo

      if( ndslfv ) then
        call spect_to_grid_rqt
     &    (rqe,rqo, rq_gr_a_1, rq_gr_a_2,levh,
     &     ls_node,ls_nodes,max_ls_nodes,
     &     lats_nodes_a,global_lats_a,lonsperlat,
     &     plnev_a,plnod_a)
        do lan=1,lats_node_a
          lat = global_lats_a(ipt_lats_node_a-1+lan)
          lon_dim = lon_dims_a(lan)
          lons_lat = lonsperlat(lat)
          do k=1,levh
            do i=1,lons_lat
              grid_gr(i,lan,G_rq+k-1)=
     &            rq_gr_a_2(i+(k-1)*lon_dim,lan)
            enddo
          enddo
        enddo
      endif

! ============ add dp =============
      if( ndslfv ) then
! prepare dp from sigio
        do lan=1,lats_node_a
          lat      = global_lats_a(ipt_lats_node_a-1+lan)
          lon_dim  = lon_dims_a(lan)
          lons_lat = lonsperlat(lat)
! !$omp parallel do schedule(dynamic,1) private(lon)
! !$omp+private(i,k,ilan,njeff)
! !$omp+private(gq,gtv,prsl,dprs)
          do lon=1,lons_lat,ngptc
            njeff = min(ngptc,lons_lat-lon+1)
            ilan  = lon - 1 
            do k=1,levs
              do i=1,njeff
                gtv(i,k) = grid_gr(ilan+i,lan,G_tt+k-1)
              enddo
            enddo
            do i=1,njeff
              gq(i) = grid_gr(ilan+i,lan,G_q)
            enddo
            if( gen_coord_hybrid ) then
              call gch2press(njeff,ngptc,gq, gtv, prsl, dprs)
            else if( hybrid ) then
              call hyb2press(njeff,ngptc,gq, prsl, dprs)
            else
              call sig2press(njeff,ngptc,gq, prsl, dprs)
            endif
            do k=1,levs
              do i=1,njeff
                dp_gr_a_2(ilan+i+(k-1)*lon_dim,lan) = dprs(i,k)
                grid_gr(ilan+i,lan,G_dp+k-1) = dprs(i,k)
              enddo
            enddo
!           do k=1,levs
!             print *,' k=',k
!             call mymaxmin(dprs(1,k),njeff,1,1,' dprs ' )
!           enddo
          enddo
        enddo
! obtain dp in spectral space
        call grid_to_spect_rqt(dp_gr_a_1,dp_gr_a_2,
     &                         dpe,dpo,levs,
     &                         ls_node,ls_nodes,max_ls_nodes,
     &                 lats_nodes_a,global_lats_a,lonsperlat,
     &                   epse,epso,plnew_a,plnow_a) 
        DO K=1,LEVS
          DO i=1,len_trie_ls
            trie_ls(i,1,P_DP+k-1)=DPE(i,1,k)
            trie_ls(i,2,P_DP+k-1)=DPE(i,2,k)
          ENDDO
          DO i=1,len_trio_ls
            trio_ls(i,1,P_DP+k-1)=DPO(i,1,k)
            trio_ls(i,2,P_DP+k-1)=DPO(i,2,k)
          ENDDO
        ENDDO

        do k=1,levs
          trie_ls(:,:,p_dpn+k-1)=trie_ls(:,:,p_dp+k-1)
          trio_ls(:,:,p_dpn+k-1)=trio_ls(:,:,p_dp+k-1)
        enddo

      endif
!
!      if(me==0) then
!        open (991,file='ps_grid',form='unformatted') 
!        write(991)grid_gr(1:lonf,1:lats_node_a,g_qm)
!        close(991)
!      endif
!  

      RETURN
      END
