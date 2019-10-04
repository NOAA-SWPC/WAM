      SUBROUTINE TREADEO_slg(NFT,FHOUR,IDATE,
     &                       GZE,QE,TEE,DIE,ZEE,RQE,
     &                       GZO,QO,TEO,DIO,ZEO,RQO,
     &                       LS_NODE,
     &                       SNNP1EV,SNNP1OD,pdryini,IPRINT,
     &                       cfile)
 
      use gfs_dyn_machine,   only: kind_evod, kind_io8

      use gfs_dyn_resol_def, only: latg, lonf, levs, levp1, latg2, 
     &                             ntcw, jcap, num_p2d, ntke,
     &                             num_p3d, ntoz, ivsinp, lnt2, ntrac
      use gfs_dyn_layout1,   only: me, len_trie_ls, len_trio_ls, 
     &                             lats_node_a, ls_dim, ls_max_node
      use namelist_dynamics_def, ONLY: igen, hybrid
      USE gfs_dyn_io_header,     ONLY: ienst, iensi, idvt, z, z_r, 
     &                                 itrun, icen2
      USE gfs_dyn_vert_def,      ONLY: si, sl
      USE sigio_module,          ONLY: sigio_head, sigio_alhead
      USE sigio_r_module,        ONLY: sigio_dbti, sigio_rropen,
     &                                 sigio_rrhead, sigio_rrdbti
      USE gfs_dyn_coordinate_def, ONLY: ck, bkl, dbk, ak5, bk5,
     &                                  idsl, idvc
      use gfs_dyn_physcons, ONLY: rerth => con_rerth, grav => con_g
!     use gfs_dyn_physcons, ONLY: rerth => con_rerth, grav => con_g, 
!    &                            rkap => con_rocp
!
      implicit none
      character*(*) cfile
      INTEGER              NFT
      REAL(KIND=KIND_EVOD) FHOUR
      INTEGER              IDATE(4),NTRACI, ntozi, ntcwi, ncldi, ixgr
     &,                    nt0
!
      REAL(KIND=KIND_EVOD) GZE(LEN_TRIE_LS,2)
     &,                     QE(LEN_TRIE_LS,2)
     &,                    TEE(LEN_TRIE_LS,2,LEVS)
     &,                    DIE(LEN_TRIE_LS,2,LEVS)
     &,                    ZEE(LEN_TRIE_LS,2,LEVS)
     &,                    RQE(LEN_TRIE_LS,2,LEVS,ntrac)
     &,                    GZO(LEN_TRIO_LS,2)
     &,                     QO(LEN_TRIO_LS,2)
     &,                    TEO(LEN_TRIO_LS,2,LEVS)
     &,                     DIO(LEN_TRIO_LS,2,LEVS)
     &,                    ZEO(LEN_TRIO_LS,2,LEVS)
     &,                    RQO(LEN_TRIO_LS,2,LEVS,ntrac)
 
      integer              ls_node(ls_dim,3)
!
!cmr  ls_node(1,1) ... ls_node(ls_max_node,1) : values of L
!cmr  ls_node(1,2) ... ls_node(ls_max_node,2) : values of jbasev
!cmr  ls_node(1,3) ... ls_node(ls_max_node,3) : values of jbasod
!
      REAL(KIND=KIND_EVOD) SNNP1EV(LEN_TRIE_LS), SNNP1OD(LEN_TRIO_LS)
      integer              i,j,k,l,n,lv,kk,lan,lat,lons_lat,lon,nn,locl
     &,                    iprint
      integer              indev, indod
      integer              indev1,indev2,indod1,indod2
     &,                    indlsev,indlsod,jbasev,jbasod
!
!     REAL(KIND=KIND_EVOD), target ::  TRISCA(LNT2)
      real(kind=kind_io8), allocatable, target :: trisca(:)
      REAL(KIND=KIND_IO8)   PDRYINI
!
      type(sigio_head) head
      type(sigio_dbti) dati
!
      integer              iret, num_dta
      real(kind=kind_evod) ga2, psurfff, pressk
!
!     integer kmsk(lonf,latg)
!     real(kind=kind_io8), target ::  buff1(lonf*latg)
!     real(kind=kind_io8), dimension(lonf,lats_node_a) ::  buffo, buff2

      integer, allocatable                     :: kmsk(:,:)
      real(kind=kind_io8), allocatable, target :: buff1(:)
      real(kind=kind_io8), allocatable         :: buffo(:,:), buff2(:,:)
 
      INCLUDE 'function2'

      call sigio_rropen(nft,cfile,iret)
      call sigio_alhead(head,iret)
      call sigio_rrhead(nft,head,iret)
!
      ivsinp = head%ivs
      if (me == 0) then
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

      if (hybrid .and. idvc .eq. 2) then				! hmhj
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
          if(me == 0)print 200,k,dbk(k),ck(k)
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
      else
        print *,' Non compatible Initial state IDVC=',head%idvc
     &,' iret=',iret
        call MPI_QUIT(333)
      endif
!
      FHOUR       = head%fhour
      idate       = head%idate
      itrun       = head%itrun
      icen2       = head%icen2
      igen        = head%igen
      ienst       = head%iens(1)
      iensi       = head%iens(2)
!     runid       = head%idrun
!     usrid       = head%idusr
      if (fhour > 0.0 .and. head%nxss == 0 .and.
     &    head%pdryini > 0.0 ) then
        if (pdryini == 0.0) pdryini = head%pdryini
      endif
      if (me == 0) print *,' IN TREAD PDRYINI=',pdryini,
     &                     ' head=',head%pdryini
      ntraci = head%ntrac
      if (head%idvt > 0.0) then
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
      elseif(ntraci == 2) then
        ntozi = 2
        ntcwi = 0
        ncldi = 0
      elseif(ntraci == 3) then
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
      IF (me == 0) THEN
        write(*,*)'nfile,in treadeo fhour,idate=',nft,fhour,idate
     &, ' ntozi=',ntozi,' ntcwi=',ntcwi,' ncldi=',ncldi
     &, ' ntraci=',ntraci,' tracers=',head%ntrac,' vtid=',head%idvt
     &,  head%ncldt,' idvc=',head%idvc,' jcap=',head%jcap
     &, ' ixgr=',ixgr,' pdryini=',pdryini
      ENDIF
!jfe
      IF(IPRINT == 1) PRINT *,'TREAD UNIT,FHOUR,IDATE=',NFT,FHOUR,IDATE
 
!
      allocate(trisca(lnt2))

      dati%i = 1                                           ! hs
      dati%f => TRISCA
      call sigio_rrdbti(nft,head,dati,iret)
      if (me == 0) print *,' Z_R=',trisca(1:10),' iret=',iret
      Z   = TRISCA
      Z_R = TRISCA
      CALL TRISEORI(TRISCA,GZE,GZO,1,LS_NODE)
      GA2 = GRAV/(RERTH*RERTH)
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
            GZE(INDEV,1) = GZE(INDEV,1)*SNNP1EV(INDEV)*GA2
            GZE(INDEV,2) = GZE(INDEV,2)*SNNP1EV(INDEV)*GA2
         END DO
      END DO
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
            GZO(INDOD,1) = GZO(INDOD,1)*SNNP1OD(INDOD)*GA2
            GZO(INDOD,2) = GZO(INDOD,2)*SNNP1OD(INDOD)*GA2
         END DO
      END DO
 
      if (mod(head%idvm/10,10) == 3 .and. me == 0)then
        print *,' CPI=',head%cpi(1:ntraci+1)
        print *,' RI=',head%ri(1:ntraci+1)
      endif
      dati%i = 2                               ! Surface pressure
      dati%f => TRISCA                       
      call sigio_rrdbti(nft,head,dati,iret)
      if (me == 0) print *,' SFCPRES=',trisca(1:10)
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
!
!
      do k=1,ntrac
!$omp parallel do private(i,lv)
        do lv=1,levs
          do i=1,len_trie_ls
            rqe(i,1,lv,k) = 0.0
            rqe(i,2,lv,k) = 0.0
          enddo
          do i=1,len_trio_ls
            rqo(i,1,lv,k) = 0.0
            rqo(i,2,lv,k) = 0.0
          enddo
        enddo
      enddo
      if (ntke > 0) then
        rqe(1,1,:,ntke) = 1.0e-7
      endif
!
      DO K=1,ntraci
        kk = 0
        if (k == 1) then
          kk = 1
        elseif (k == ntozi) then
          kk = ntoz
        elseif (k >= ntcwi .and. k < ntcwi+ncldi-1) then
          do n=1,ncldi
            if (k == ntcwi+n-1) kk = ntcw+n-1
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
!     if (((ixgr .eq. 4 .and. num_p3d .eq. 4) .or.   ! Zhao Scheme!
!    &     (ixgr .eq. 5 .and. num_p3d .eq. 3))       ! Ferrier Scheme!
!    &     .and. fhour .gt. 0.1) then
!       allocate(kmsk(lonf,latg), buff1(lonf,latg))
!       allocate(buff2(lonf,lats_node_a),buffo(lonf,lats_node_a))
!       kmsk(:,:)  = 0
!       num_dta = (ntraci+3)*levs + 2
!       do nn=1,num_p3d
!         do k=1,levs
!           dati%i = num_dta + (nn-1)*levs + k      ! physics 3D grid fields
!           dati%f => buff1
!           call sigio_rrdbti(nft,head,dati,iret)
!            call split2d_r(buff1(1),buffo,global_lats_a)
!            CALL interpred(1,kmsk,buffo,buff2,global_lats_a,
!     &                     lonsperlat)
!
!           do lan=1,LATS_NODE_A
!             lat = global_lats_a(ipt_lats_node_a-1+lan)
!             lons_lat = lonsperlat(lat)
!             DO i=1,lons_lat
!               phy_f3d(i,k,nn,lan) = buff2(i,lan)
!             enddo
!           enddo
!         enddo
!       enddo
!       do nn=1,num_p2d
!         dati%i = num_dta + num_p3d*levs + nn      ! physics 2D grid fields
!         dati%f => buff1
!         call sigio_rrdbti(nft,head,dati,iret)
!          call split2d_r(buff1,buffo,global_lats_a)
!          CALL interpred(1,kmsk,buffo,buff2,global_lats_a,
!     &                   lonsperlat)
!         do j=1,lats_node_a
!           do i=1,lonf
!             phy_f2d(i,nn,j) = buff2(i,j)
!           enddo
!         enddo
!       enddo
!       deallocate (kmsk, buffo, buff2)
!     endif
      if (head%nxss > 0) then
        if (.not. allocated(buff1)) allocate(buff1(1))
        dati%i = num_dta + num_p3d*levs + num_p2d + 1    ! pdryini
        dati%f => buff1
        call sigio_rrdbti(nft,head,dati,iret)
        pdryini = buff1(1)
        deallocate(buff1)
      endif
!
      if (allocated(trisca)) deallocate (trisca)
      iprint = 0
 
      END SUBROUTINE TREADEO_slg
