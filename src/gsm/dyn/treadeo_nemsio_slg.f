      subroutine treadeo_nemsio_slg(cfile,IDATE,tlev,restart_run
     &,                             GZE,QE,TEE,DIE,ZEE,RQE
     &,                             GZO,QO,TEO,DIO,ZEO,RQO
     &,                             grid_gr
     &,                             LS_NODE,LS_NODES,MAX_LS_NODES
     &,                             pdryini,IPRINT
     &,                             global_lats_a,lats_nodes_a
     &,                             lonsperlat
     &,                             epse,epso,plnew_a,plnow_a
     &,                             snnp1ev,snnp1od)
!    &,                             plnev_a,plnod_a,snnp1ev,snnp1od)
!    &                              pwat,ptot,ptrc,slg_flag)

!!
!! Revision history:
!         2008    Henry Juang, original code
!  Nov 23 2009    Sarah Lu, tracer read-in is generalized (loop through ntrac, with
!                 tracer name specified in gfs_dyn_tracer_config)
!  Apr 09 2010    Sarah Lu, set rqg initial value to 1.e-15
!  Aug 17 2010    Sarah Lu, clean debug print
!  Aug 25 2010    Sarah Lu, modified to compute tracer global sum
!  Sep 08 2010    Jun Wang, change to nemsio format file
!  Dec 16 2010    Jun Wang, change to nemsio library
!  Feb 20 2011    Henry Juang, change to have mass dp and ndslfv
!  Feb 26 2011    Sarah Lu, modify to read both cold-start (from chgres) and 
!                 warm-start (from replay) ICs
!  Jun    2011    Jun Wang   reading fields from grib file with either w3_d 
!                 or w3_4 lib
!  Nov 11 2011    Sarah Lu, change floor value for tracer initial values;
!                 remove nvcoord read-in and check              
!  Sep 20 2012    Jun Wang, set n time step trie/o to be consistnet with sigio input
!  Feb 18 2015    S Moorthi - adapted to two time level semi_Lagrangian model
!  Nov ?? 2015    S Moorthi - fix for restarting from nemsio history file
!  
 
      use gfs_dyn_resol_def
      use gfs_dyn_layout1
      use gfs_dyn_coordinate_def
      use gfs_dyn_io_header
      use namelist_dynamics_def
      use gfs_dyn_vert_def
      use gfs_dyn_mpi_def
      use gfs_dyn_physcons, rerth => con_rerth, grav => con_g,
     &                      rkap  => con_rocp,  cpd  => con_cp
      use nemsio_module
      use nemsio_def
      USE gfs_dyn_date_def,      ONLY : FHOUR
      use gfs_dyn_tracer_config, ONLY : gfs_dyn_tracer   ! generalized tracer
!     use gfs_dyn_io_header, only : z, z_r

!
      IMPLICIT NONE
!
      real(kind=kind_evod), parameter :: pa2cb=0.001
      character*(*) cfile
      INTEGER       IDATE(4), idate7(7)
     &,             levsi, jcapi, latgi, lonfi, tlev
      logical       restart_run
!     logical       slg_flag
!
      REAL(KIND=KIND_GRID) GRID_GR(lonf,lats_node_a_max,lotgr)
      real(kind=kind_evod), dimension(len_trie_ls,2) :: gze, qe
      real(kind=kind_evod), dimension(len_trio_ls,2) :: gzo, qo
      real(kind=kind_evod), dimension(len_trie_ls,2,levs) :: tee,die,zee
      real(kind=kind_evod), dimension(len_trio_ls,2,levs) :: teo,dio,zeo
      real(kind=kind_evod), dimension(len_trie_ls,2,levs,ntrac) :: rqe
      real(kind=kind_evod), dimension(len_trio_ls,2,levs,ntrac) :: rqo

!!
      real(kind=kind_evod), dimension(len_trie_ls) :: epse, snnp1ev
      real(kind=kind_evod), dimension(len_trio_ls) :: epso, snnp1od
!!
      real(kind=kind_evod), dimension(len_trie_ls,latg2) :: plnew_a
!    &,                     plnev_a
      real(kind=kind_evod), dimension(len_trio_ls,latg2) :: plnow_a
!    &,                     plnod_a
!!
!     real(kind=kind_evod) trie_ls(len_trie_ls,2,lotls)
!     real(kind=kind_evod) trio_ls(len_trio_ls,2,lotls)
!
      integer              ls_node(ls_dim,3)
!
!cmr  ls_node(1,1) ... ls_node(ls_max_node,1) : values of L
!cmr  ls_node(1,2) ... ls_node(ls_max_node,2) : values of jbasev
!cmr  ls_node(1,3) ... ls_node(ls_max_node,3) : values of jbasod
!
      INTEGER, dimension(nodes) :: MAX_LS_NODES, lats_nodes_a
      INTEGER              LS_NODES(LS_DIM,NODES)
      INTEGER              J,K,L,LOCL,N,lv,kk,w3rlkind,w3ikind,iprint
     &,                    i,lan,lat,iblk,lons_lat,il,lon,njeff,nn
     &,                    indev,indod,indev1,indev2,indod1,indod2
     &,                    nfhour,nfminute,nfsecondn,nfsecondd 
!
      REAL(KIND=KIND_EVOD) XI(LEVP1),sikp1(levp1),XL(LEVS)
      REAL(KIND=KIND_IO4)  VTID,RUNID4,pdryini4,XNCLD,xgf
     &,                    TRUN,WAVES,XLAYERS

      REAL(KIND=KIND_grid)  PDRYINI
      real(kind=kind_io4), allocatable :: vcoord4(:,:,:)
! for generalized tracers
      integer                      nreci
      character*8               :: vname
      character*8, allocatable  :: recnamei(:), reclevtypi(:)
      integer,     allocatable  :: reclevi(:)
!
      integer              iret, num_dta, ijm, tlmeta
      real(kind=kind_evod) ga2, psurfff, pressk, tem
      real(kind=kind_evod), parameter :: rkapi=1.0/rkap,
     &                                   rkapp1=1.0+rkap
!
      integer kmsk(lonf,latg), global_lats_a(latg), lonsperlat(latg)
      real(kind=kind_io8), dimension(lonf,lats_node_a) ::  buffo, buff2
      real(kind=kind_evod) teref(levp1),ck5p(levp1)
!    &,                    ttref(levp1)
!
      real (kind=kind_io4), allocatable ::  nemsio_data(:)
!!
      real(kind=kind_grid), dimension(lonf,lats_node_a) :: zsg, psg
     &,                       pwat, ptot
      real(kind=kind_grid), dimension(lonf,lats_node_a,levs) :: uug
     &,                       vvg, ttg
!    &,                       vvg, ttg, dpg
      real(kind=kind_grid) rqg(lonf,lats_node_a,levh)
!
!     REAL(KIND=KIND_GRID) ptrc   (lonf,lats_node_a,ntrac)                !glbsum
!
      real(8) timef,stime,etime
!
      INTEGER              INDLSEV,JBASEV,INDLSOD,JBASOD
      INCLUDE 'function2'
!------------------------------------------------------------------------
!       Input file is in grid-point space - use gfs_io package
!
      if (me == 0) write(0,*)' before nemsio_open cfile=',cfile

      stime = timef()
      call nemsio_open(gfile_in,trim(cfile),'read',iret)
      etime = timef()

      if (me == 0) write(0,*)'in read nemsio file, open time='
     &,                       timef()-stime
     &,     ' lonf=',lonf,' lats_node_a=',lats_node_a
!
      call nemsio_getfilehead(gfile_in,iret=iret,
     &                        version=ivsupa,idate=idate7,
     &                        nfhour=nfhour,nfminute=nfminute,
     &                        nfsecondn=nfsecondn,nfsecondd=nfsecondd,
     &                        dimy=latgi,dimx=lonfi,dimz=levsi,
     &                        jcap=jcapi,idvc=idvc,
     &                        ncldt=ncldt,tlmeta=tlmeta)

!     write(0,*)' after nemsio_getfilehead levsi=',levsi,latgi,lonfi
!    &,' me=',me
       idate(1)   = idate7(4)
       idate(2:3) = idate7(2:3)
       idate(4)   = idate7(1)
!
      call nemsio_getheadvar(gfile_in,'iorder',iorder,iret=iret)
      call nemsio_getheadvar(gfile_in,'irealf',irealf,iret=iret)
      call nemsio_getheadvar(gfile_in,'igen',igen,iret=iret)
      call nemsio_getheadvar(gfile_in,'dimx',lonb,iret=iret)
      call nemsio_getheadvar(gfile_in,'dimy',latb,iret=iret)
      call nemsio_getheadvar(gfile_in,'icen2',icen2,iret=iret)
      call nemsio_getheadvar(gfile_in,'iens',iens,iret=iret)
      call nemsio_getheadvar(gfile_in,'idpp',idpp,iret=iret)
      call nemsio_getheadvar(gfile_in,'idrun',idrun,iret=iret)
      call nemsio_getheadvar(gfile_in,'itrun',itrun,iret=iret)
      call nemsio_getheadvar(gfile_in,'idusr',idusr,iret=iret)
      call nemsio_getheadvar(gfile_in,'pdryini',pdryini4,iret=iret)

!      call gfsio_getfilehead(gfile_in,iret=iret,
!     &  version=ivsupa,fhour=fhour4,idate=idate,
!     &  latb=latb,lonb=lonb,levs=levsi,jcap=jcapi,itrun=itrun,
!     &  iorder=iorder,irealf=irealf,igen=igen,latf=latgi,lonf=lonfi,
!     &  latr=latri,lonr=lonfi,ntrac=ntraci,icen2=icen2,iens=iens,
!     &  idpp=idpp,idsl=idsl,idvc=idvc,idvm=idvm,idvt=idvt,idrun=idrun,
!     &  idusr=idusr,pdryini=pdryini4,ncldt=ncldt,nvcoord=nvcoord)
!
      if (me == 0) then
        write(0,*)'iret=',iret,
     &     ' levsi=',levsi,' idate=',idate,
     &   'lonf=',lonf,'lonfi=',lonfi,'latg=',latg,'latgi=',latgi,
     &   'jcap=',jcap,'jcapi=',jcapi,'levs=',levs,'levsi=',levsi,
     &   'idvc=',idvc,'tlmeta=',tlmeta,
     &   'gen_coord_hybrid=',gen_coord_hybrid,'pdryini=',pdryini
        if(lonf .ne. lonfi .or. latg .ne. latgi .or.
     &     jcap .ne. jcapi .or. levs .ne. levsi) then
          print *,' Input resolution and the model resolutions are'
     &,  ' different- run aborted'
          call mpi_quit(555)
        endif
        if ( gen_coord_hybrid ) then
          if(me==0) print *, ' Use sigma-theta-p hybrid coordinate'
          if (idvc == 3 ) then
           if(me==0)   
     &       print *, ' Cold_start input is consistent, run continues'
          else 
           if(me==0)
     &       print *, ' Cold_start input is different, run aborted'
           call mpi_quit(556)
          endif
        endif   
        if ( hybrid ) then
          if(me==0)print *, ' Use sigma-p hybrid coordinate'
          if (idvc == 2 ) then
           if(me==0)
     &      print *, ' Cold_start input is consistent, run continues'
          else 
           if(me==0)
     &       print *, ' Cold_start input is different, run aborted'
           call mpi_quit(557)
          endif
        endif   
      endif
!
      allocate (vcoord4(levsi+1,3,2))
      allocate (vcoord(levsi+1,3))
      call nemsio_getfilehead(gfile_in,iret=iret,vcoord=vcoord4)
!
!     if (me == 0) then
!     print *,' nvcoord=',nvcoord,' vcoord4=',vcoord4(:,1:nvcoord)
!    &,' iret=',iret
!     endif
!
      vcoord(:,1:3) = vcoord4(:,1:3,1)
!     if (me .eq. 0) print *,' vcoord=',vcoord(:,1:nvcoord)
      deallocate (vcoord4)
!
! for generalized tracers
! retrieve nreci, recnamei, reclevtypi, and reclevi
      call nemsio_getfilehead(gfile_in,iret=iret,nrec=nreci)
      if (me == 0) then
        print *, 'LU_TRC: nreci =', nreci, iret
      endif

      allocate (recnamei(nreci))
      allocate (reclevtypi(nreci))
      allocate (reclevi(nreci))
      call nemsio_getfilehead(gfile_in,iret=iret,recname=recnamei,
     &                       reclevtyp=reclevtypi,reclev=reclevi)
       stime=timef()
!        print *,'after nemsioheader,time=',stime-etime

!
      if (gen_coord_hybrid) then

!   ak bk ck in file have the same order as model
        do k=1,levp1
          ak5(k) = vcoord(k,1)/1000.
          bk5(k) = vcoord(k,2)
          ck5(k) = vcoord(k,3)/1000.
        enddo
        vertcoord_id=0
        do k=1,levp1
          if( ck5(k).ne.0.0 ) vertcoord_id=3
        enddo
! provide better estimated press
        psurfff = 101.3
        if( thermodyn_id.eq.3 ) then
          do k=1,levs
            thref(k) = 300.0*cpd
            teref(k) = 255.0*cpd
          enddo
        else
         do k=1,levp1
          thref(k) = 300.0
          teref(k) = 255.0
         enddo
        endif
        ck5p(levp1) = ck5(levp1)
        do k=1,levs
          ck5p(k) = ck5(k)*(teref(k)/thref(k))**rkapi
        enddo
        if( me.eq.0 ) then
          do k=1,levp1
            pressk=ak5(k)+bk5(k)*psurfff+ck5p(k)
            print 180,k,ak5(k),bk5(k),ck5(k),pressk
180         format('k=',i3,'  ak5=',f13.6,'  bk5=',e13.5,
     &            '   ck5=',f13.6,'  closed pressk=',f10.6)
          enddo
        endif
        do k=1,levs
          dbk(k) = bk5(k)-bk5(k+1)
          bkl(k) = (bk5(k)+bk5(k+1))*0.5
          ck(k)  = ak5(k)*bk5(k+1)-ak5(k+1)*bk5(k)
        enddo
        do k=1,levp1
          si(k) = ak5(k)/psurfff + bk5(k) + ck5p(k)/psurfff
        enddo
        do k=1,levs
          sl(k) = 0.5*(si(k)+si(k+1))
        enddo

      else if (hybrid .and. idvc .eq. 2) then
!       idsl=slid  !=2,pk=0.5*(p(k+1/2)+p(k-1/2)) check alfa(1)  am_bm
!   ak bk order in "sigma" file is bottom to top !!!!!!!!!!!!!!!!!!
        psurfff = 101.3
        ck5=0.
        do k=1,levp1
          ak5(k) = vcoord(levp1+1-k,1)/1000.
          bk5(k) = vcoord(levp1+1-k,2)
          pressk = ak5(k) + bk5(k)*psurfff

          if(me.eq.0)print 190,k,ak5(k),bk5(k),pressk
190       format('k=',i3,'  ak5=',E14.6,'  bk5=',e14.6,
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
          si(levs+2-k) = ak5(k)/psurfff + bk5(k) !ak(k) bk(k) go top to bottom
        enddo
        do k=1,levs
          sl(k) = 0.5*(si(k)+si(k+1))
        enddo
!
      elseif (idvc .le. 1) then
        si(:)    = vcoord(:,1)
        sik(:)   = si(:) ** rkap
        sikp1(:) = si(:) ** rkapp1
        do k=1,levs
          tem      = rkapp1 * (si(k) - si(k+1))
          slk(k)   = (sikp1(k)-sikp1(k+1))/tem
          sl(k)    = slk(k) ** rkapi
!         sl(k)    = ((sikp1(k)-sikp1(k+1))/tem)**rkapi
          if (me .eq. 0) print 250, k, si(k), sl(k)
250       format('k=',i3,'  si=',f9.6,'  sl=',e13.5)
        enddo
      else
        print *,' Non compatible Initial state IDVC=',idvc
     &,' iret=',iret
        call MPI_QUIT(560)
      endif
!
      FHOUR       = real(nfhour,8)+real(nfminute,8)/60.+          
     &              real(nfsecondn,8)/(real(nfsecondd,8)*3600.)
      if (tlev == 1) fhini = fhour
!     WAVES       = jcap
!     XLAYERS     = levs
      itrun       = itrun
      icen        = 7
      icen2       = icen2
      igen        = igen
      ienst       = iens(1)
      iensi       = iens(2)
!Weiyu, pdryini4 did not initiaqlized, temp give value zero.
!-----------------------------------------------------------
      pdryini4 = 0.0
      if (pdryini == 0.0) pdryini = pdryini4
!
!
      IF (me == 0) THEN
        write(0,*)'cfile,in treadeo fhour,idate=',cfile,fhour,idate
     &, ' idvc=',idvc,' jcap=',jcap, ' pdryini=',pdryini
      ENDIF
!
      allocate (nemsio_data(lonf*latg))
!  Read orog
       stime=timef()

      call w3kind(w3rlkind,w3ikind)
!
      if(w3rlkind==8) then
        call nemsio_readrecv(gfile_in,'hgt','sfc',1,nemsio_data,
     &                       iret=iret)
      elseif(w3rlkind==4) then
        call nemsio_readrecvw34(gfile_in,'hgt','sfc',1,nemsio_data, 
     &                          iret=iret)
      endif

!      print *,'in treadeo,time=',timef()-stime,'hgt=',
!     &  maxval(nemsio_data),minval(nemsio_data), 'iret=',iret
      call split2d(nemsio_data,buffo,global_lats_a)
!      print *,'in treadeo,buffo=',maxval(buffo),minval(buffo)
      CALL interpred(1,kmsk,buffo,zsg,global_lats_a,lonsperlat)
!      print *,'in treadeo,zgs=',maxval(zsg),minval(zsg)
      ijm = lonf*lats_node_a
!
!  Read ps
      if(w3rlkind==8) then
        call nemsio_readrecv(gfile_in,'pres','sfc',1,nemsio_data,
     &                       iret=iret)
      elseif(w3rlkind==4) then
        call nemsio_readrecvw34(gfile_in,'pres','sfc',1,nemsio_data,
     &                          iret=iret)
      endif
      call split2d(nemsio_data,buffo,global_lats_a)
      CALL interpred(1,kmsk,buffo,psg,global_lats_a,lonsperlat)
!
!  Read u
      do k=1,levs
       if(w3rlkind==8) then
         call nemsio_readrecv(gfile_in,'ugrd','mid layer',k,nemsio_data,
     &                        iret=iret)
       elseif(w3rlkind==4) then
         call nemsio_readrecvw34(gfile_in,'ugrd','mid layer',k,
     &                           nemsio_data,iret=iret)
       endif
!      print *,'in treadeo,ugrd=',maxval(nemsio_data),
!    & minval(nemsio_data),'iret=',iret,'k=',k,' me=',me
!    &,' xdim=',size(uug,dim=1),' ydim=',size(uug,dim=2)
!    &,' zdim=',size(uug,dim=3)
       call split2d(nemsio_data,buffo,global_lats_a)
       CALL interpred(1,kmsk,buffo,uug(1,1,k),global_lats_a,lonsperlat)
      enddo
!  Read v
      do k=1,levs
       if(w3rlkind==8) then
        call nemsio_readrecv(gfile_in,'vgrd','mid layer',k,nemsio_data,
     &                       iret=iret)
       elseif(w3rlkind==4) then
        call nemsio_readrecvw34(gfile_in,'vgrd','mid layer',k,
     &                          nemsio_data,iret=iret)
       endif
       call split2d(nemsio_data,buffo,global_lats_a)
       CALL interpred(1,kmsk,buffo,vvg(1,1,k),global_lats_a,lonsperlat)
      enddo
!  Read T   -- this is real temperature
      do k=1,levs
       if(w3rlkind==8) then
        call nemsio_readrecv(gfile_in,'tmp','mid layer',k,nemsio_data,
     &                       iret=iret)
       elseif(w3rlkind==4) then
        call nemsio_readrecvw34(gfile_in,'tmp','mid layer',k,
     &                          nemsio_data,iret=iret)
       endif
        call split2d(nemsio_data,buffo,global_lats_a)
        CALL interpred(1,kmsk,buffo,ttg(1,1,k),global_lats_a,lonsperlat)
      enddo

!  Read dp 
!     do k=1,levs
!      if(w3rlkind==8) then
!       call nemsio_readrecv(gfile_in,'dpres','mid layer',k,nemsio_data,
!    &     iret=iret)
!      elseif(w3rlkind==4) then
!       call nemsio_readrecvw34(gfile_in,'dpres','mid layer',k,
!    &     nemsio_data,iret=iret)
!      endif
!       call split2d(nemsio_data,buffo,global_lats_a)
!       CALL interpred(1,kmsk,buffo,dpg(1,1,k),global_lats_a,lonsperlat)
!     enddo
!     if( me.eq.0 ) print *,' read dpg ',(dpg(1,1,k),k=1,levs)
!
!  Initial Tracers with zero
!
!*    rqg(:,:,:) = 0.0
      rqg(:,:,:) = 0.0

!! Generalized tracers: 
!! Loop through ntrac to read in met + chem tracers
!*
      do n = 1, ntrac
        vname = trim(gfs_dyn_tracer%vname(n, 1))
        if(me==0) print *,'LU_TRC: initialize ',n,vname
     &,' w3rlkind=',w3rlkind
        do k=1,levs
         if(w3rlkind==8) then
          call nemsio_readrecv(gfile_in,trim(vname),
     &                        'mid layer',k,nemsio_data,iret=iret)
         elseif(w3rlkind==4) then
          call nemsio_readrecvw34(gfile_in,trim(vname),
     &                           'mid layer',k,nemsio_data,iret=iret)
         endif
          if(iret == 0) then
!*          if(me==0) print *,'LU_TRC: tracer read in ok -',
!*   &                gfs_dyn_tracer%vname(n, 1),k
            if(me==0 .and. k==1) 
     &                print *,'LU_TRC: tracer read in ok '
!           if (me == 0 .and. n == 1) then
!             write(0,*)' k=',k,' rqg=',minval(nemsio_data)
!    &,maxval(nemsio_data)
!           endif

            call split2d(nemsio_data,buffo,global_lats_a)
            CALL interpred(1,kmsk,buffo,rqg(1,1,k+(n-1)*levs),
     &                     global_lats_a,lonsperlat)
          else
!*          if(me==0) print *,'LU_TRC: tracer not found in input; ',
!*   &         'set chem tracer to default values',me,k
            if(me==0 .and. k==1)  print *,
     &         'LU_TRC: tracer not found in input; set to default'
          endif
        enddo
      enddo       
      etime = timef()
!        print *,'after nemsioheader,time=',etime-stime
!
!   Convert from Gaussian grid to spectral space
!   including converting to model_uvtp if necessary
!
!       write(0,*)' num_pes_fcst=',num_pes_fcst,' me=',me
!    &,' tlev=',tlev
      if(me < num_pes_fcst) then

!       write(0,*)' before  grid_to_spect_inp_slg'
!    &,' zsg=',maxval(zsg),minval(zsg),' psg=',maxval(psg)
!    &,minval(psg),' tlev=',tlev

!       do j=1,lats_node_a          ! save height for output - Moorthi
!         do i=1,lonf
!           gz_grid(i,j) = zsg(i,j)
!         enddo
!       enddo
!
       if (tlev == 1) then
!$omp parallel do private(i,j,lat,lons_lat)
          do j = 1, lats_node_a
            lat      = global_lats_a(ipt_lats_node_a-1+j)
            lons_lat = lonsperlat(lat)
            do i=1,lons_lat
              grid_gr(i,j,g_gz)  = zsg(i,j)
              grid_gr(i,j,g_qm ) = log(psg(i,j)*pa2cb)
            enddo
          enddo
!$omp parallel do private(i,j,k,lat,lons_lat)
          do k=1,levs
            do j = 1, lats_node_a
              lat      = global_lats_a(ipt_lats_node_a-1+j)
              lons_lat = lonsperlat(lat)
              do i=1,lons_lat
                grid_gr(i,j,g_ttm+k-1) = ttg(i,j,k)
                grid_gr(i,j,g_uum+k-1) = uug(i,j,k)
                grid_gr(i,j,g_vvm+k-1) = vvg(i,j,k)
!               grid_gr(i,j,g_dpm+k-1) = dpg(i,j,k)
              enddo
            enddo
!          write(1000+me,*)'tmin=',minval(ttg(1:lons_lat,:,k)),' tmax=',
!    &maxval(ttg(1:lons_lat,:,k)),' umin=',minval(uug(1:lons_lat,:,k)),
!    &' umax=', maxval(uug(1:lons_lat,:,k)),' vmin=',
!    & minval(vvg(1:lons_lat,:,k)),' vmax=',
!    & maxval(vvg(1:lons_lat,:,k)),' k=',k
          enddo
!$omp parallel do private(lan,lat,lons_lat)
          do lan=1,lats_node_a
            lat      = global_lats_a(ipt_lats_node_a-1+lan)
            lons_lat = lonsperlat(lat)
            call hyb2press(lons_lat, lonf
     &,                    grid_gr(1:lonf,lan,g_qm)
     &,                    grid_gr(1:lonf,lan,g_p:g_p+levs-1)
     &,                    grid_gr(1:lonf,lan,g_dpm:g_dpm+levs-1))
          enddo
!$omp parallel do private(i,j,k)
          do k=1,levh
            do j = 1, lats_node_a
              do i=1,lonf
                grid_gr(i,j,g_rm+k-1) = rqg(i,j,k)
              enddo
            enddo
          enddo
        endif
!
        if (.not. restart_run .and. fhini == fhrot) then
!$omp parallel do private(i,j)
          do j = 1, lats_node_a
            do i=1,lonf
!             grid_gr(i,j,g_q ) = log(psg(i,j)*pa2cb)
              grid_gr(i,j,g_q ) = grid_gr(i,j,g_qm)
            enddo
          enddo
!$omp parallel do private(i,j,k)
          do k=1,levs
            do j = 1, lats_node_a
              do i=1,lonf
                grid_gr(i,j,g_tt+k-1) = ttg(i,j,k)
                grid_gr(i,j,g_uu+k-1) = uug(i,j,k)
                grid_gr(i,j,g_vv+k-1) = vvg(i,j,k)
!               grid_gr(i,j,g_dp+k-1) = grid_gr(i,j,g_dpm+k-1)
              enddo
            enddo
          enddo
!$omp parallel do private(i,j,k)
          do k=1,levh
            do j = 1, lats_node_a
              do i=1,lonf
                grid_gr(i,j,g_rq+k-1) = rqg(i,j,k)
              enddo
            enddo
          enddo
        elseif (tlev == 2) then
!$omp parallel do private(i,j)
          do j = 1, lats_node_a
            do i=1,lonf
              grid_gr(i,j,g_gz) = zsg(i,j)
              grid_gr(i,j,g_q ) = log(psg(i,j)*pa2cb)
            enddo
          enddo
!$omp parallel do private(i,j,k)
          do k=1,levs
            do j = 1, lats_node_a
              do i=1,lonf
                grid_gr(i,j,g_tt+k-1) = ttg(i,j,k)
                grid_gr(i,j,g_uu+k-1) = uug(i,j,k)
                grid_gr(i,j,g_vv+k-1) = vvg(i,j,k)
!               grid_gr(i,j,g_dp+k-1) = dpg(i,j,k)
              enddo
            enddo
          enddo
!$omp parallel do private(i,j,k)
          do k=1,levh
            do j = 1, lats_node_a
              do i=1,lonf
                grid_gr(i,j,g_rq+k-1) = rqg(i,j,k)
              enddo
            enddo
          enddo
        endif

!     write(0,*)' me=',me,' calling grid_to_spect_inp_slg'
        call grid_to_spect_inp_slg
     &                    (zsg,psg,uug,vvg,ttg,rqg,
     &                     GZE,QE,TEE,DIE,ZEE,RQE,
     &                     GZO,QO,TEO,DIO,ZEO,RQO,
     &                     ls_node,ls_nodes,max_ls_nodes,
     &                     lats_nodes_a,global_lats_a,lonsperlat,
     &                     epse,epso,plnew_a,plnow_a)

        call triseof(gze,gzo,z_r,1,ls_node)
        z = z_r
!       write(0,*)' after  grid_to_spect_inp_slg '
!       write(0,*)' tee=',tee(1,1,:)
!       write(0,*)' qe=',qe(1,1),' gze=',gze(1,1)

!  convert to laplacian of terrain for divergence
        ga2 = grav/(rerth*rerth)
        do locl=1,ls_max_node
                L = ls_node(locl,1)
           jbasev = ls_node(locl,2)
           indev1 = indlsev(L,L)
           if (mod(L,2) == mod(jcap+1,2)) then
              indev2 = indlsev(jcap+1,L)
           else
              indev2 = indlsev(jcap  ,L)
           endif
           do indev = indev1 , indev2
             GZE(INDEV,1) = GZE(INDEV,1)*SNNP1EV(INDEV)*GA2
             GZE(INDEV,2) = GZE(INDEV,2)*SNNP1EV(INDEV)*GA2
           end do
        end do
        do locl=1,ls_max_node
                L = ls_node(locl,1)
           jbasod = ls_node(locl,3)
           indod1 = indlsod(L+1,L)
           if (mod(L,2) == mod(jcap+1,2)) then
              indod2 = indlsod(jcap  ,L)
           else
              indod2 = indlsod(jcap+1,L)
           endif
           do indod = indod1 , indod2
             GZO(INDOD,1) = GZO(INDOD,1)*SNNP1OD(INDOD)*GA2
             GZO(INDOD,2) = GZO(INDOD,2)*SNNP1OD(INDOD)*GA2
           enddo
        enddo
      endif
!
      call nemsio_close(gfile_in,iret)
!
      iprint = 0
 
!     write(0,*)' returning from treadei_nemsio_slgr tlev', tlev
!!!!
      RETURN
      END
