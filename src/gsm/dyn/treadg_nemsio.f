      SUBROUTINE TREADG_nemsio(gfilename,FHOUR,IDATE,
     &                   zsg,psg,ttg,uug,vvg,rqg,dpg,
     &                   pdryini,IPRINT,
     &                   global_lats_a,lats_nodes_a,lonsperlat)
!
!-------------------------------------------------------------------
!*** program log
!*** Dec, 2009 Jun Wang:  read grid point variables for restart
!*** Dec, 2010 Jun Wang:  change to nemsio library 
!*** Feb, 2011 Henry Juang: add dpg for mass_dp and ndsl
!*** Feb, 2011 Sarah Lu: change for new nemsio file header
!*** Nov, 2011 Jun Wang: remove nvccord from restart file
!*** Nov, 2011 Sarah Lu: change floor value for tracer initial values
!*** Mar, 2012 Jun Wang: restart for spectral his output
!-------------------------------------------------------------------
!
!!
      use gfs_dyn_resol_def
      use gfs_dyn_layout1
      use gfs_dyn_coordinate_def					
      use gfs_dyn_io_header, only : idvt,lonr,latr,icen2,ienst,iensi,
     &                              itrun,idpp,idrun,idusr,ncldt,
     &                              irealf,iorder,icen,iens
      use namelist_dynamics_def
      use gfs_dyn_vert_def
      use gfs_dyn_mpi_def
      use gfs_dyn_physcons, rerth => con_rerth
     &,             grav => con_g, rkap => con_rocp
     &,             cpd => con_cp
      use nemsio_module
      use gfs_dyn_tracer_config, ONLY : gfs_dyn_tracer   ! generalized tracer
!
      IMPLICIT NONE
      character(*) gfilename
      REAL(KIND=KIND_EVOD) FHOUR
      INTEGER              IDATE(4)
     &,                    latbi, lonbi, levsi, jcapi,
     &                     latgi, lonfi, latri, lonri,idate7(7)
!
      integer              lats_nodes_a(nodes)
      INTEGER              IPRINT
      INTEGER              J,K,L,LOCL,N,lv,kk
      integer              i,lan,lat,iblk,lons_lat,il,lon,njeff,nn
      REAL(KIND=KIND_EVOD) TRUN,WAVES,XLAYERS
      REAL(KIND=KIND_EVOD) XI(LEVP1),XL(LEVS)
      REAL(KIND=KIND_EVOD) sikp1(levp1)
      REAL(KIND=KIND_IO4)   VTID,RUNID4,fhour4,XNCLD,xgf
      REAL(KIND=KIND_grid)  PDRYINI
      real(kind=kind_io4), allocatable ::  vcoord4(:,:,:)
! for generalized tracers
      integer                      nreci
      character*16               :: vname
!
      type (nemsio_gfile) gfile_in
!
      integer              iret, num_dta
      integer              im,jm,fieldsize
      real(kind=kind_evod) psurfff
      real(kind=kind_evod) pressk, tem
      real(kind=kind_evod), parameter :: rkapi=1.0/rkap,
     &                                   rkapp1=1.0+rkap
!
      integer kmsk(lonf,latg), global_lats_a(latg), lonsperlat(latg)
      real(kind=kind_io8) buffo(lonf,lats_node_a)
     &,                   buff2(lonf,lats_node_a)
      real(kind=kind_evod) teref(levp1),ck5p(levp1)			! hmhj
!
      real (kind=kind_io8), allocatable ::  nemsio_data(:)
!!
      real(kind=kind_grid) zsg(lonf,lats_node_a)
      real(kind=kind_grid) psg(lonf,lats_node_a)
      real(kind=kind_grid) uug(lonf,lats_node_a,levs)
      real(kind=kind_grid) vvg(lonf,lats_node_a,levs)
      real(kind=kind_grid) ttg(lonf,lats_node_a,levs)
      real(kind=kind_grid) dpg(lonf,lats_node_a,levs)
      real(kind=kind_grid) rqg(lonf,lats_node_a,levh)
!
!------------------------------------------------------------------
      print *,'in treadg_nemsio,gfilename=',trim(gfilename)
!
!--- Input file is in grid-point space - use nems_io package
!
      call nemsio_init()
!
      call nemsio_open(gfile_in,trim(gfilename),'read',iret)
      print *,'after nemsio_open,iret=',iret
!
      call nemsio_getfilehead(gfile_in,iret=iret,
     &  idate=idate7,
     &  dimy=im,dimx=jm,dimz=levsi,jcap=jcapi,
     &  idvc=idvc)
!      print *,'after nemsio_getfilehead,iret=',iret
!!
      call nemsio_getheadvar(gfile_in,'fhour',fhour,iret=iret)
      call nemsio_getheadvar(gfile_in,'iorder',iorder,iret=iret)
      call nemsio_getheadvar(gfile_in,'irealf',irealf,iret=iret)
      call nemsio_getheadvar(gfile_in,'igen',igen,iret=iret)
      call nemsio_getheadvar(gfile_in,'latf',latgi,iret=iret)
      call nemsio_getheadvar(gfile_in,'lonf',lonfi,iret=iret)
      call nemsio_getheadvar(gfile_in,'icen2',icen2,iret=iret)
      call nemsio_getheadvar(gfile_in,'iens',iens,iret=iret)
      call nemsio_getheadvar(gfile_in,'idpp',idpp,iret=iret)
      call nemsio_getheadvar(gfile_in,'idrun',idrun,iret=iret)
      call nemsio_getheadvar(gfile_in,'itrun',itrun,iret=iret)
      call nemsio_getheadvar(gfile_in,'idvt',idvt,iret=iret)
      call nemsio_getheadvar(gfile_in,'pdryini',pdryini,iret=iret)
      lonr=lonfi
      latr=latgi

!
      if (me == 0) then
        print *,'iret=',iret,
     &     ' levsi=',levsi,'idate7=',idate7,
     &   'lonf=',lonf,'lonfi=',lonfi,'latg=',latg,'latgi=',latgi,
     &   'jcap=',jcap,'jcapi=',jcapi,'levs=',levs,'levsi=',levsi,
     &   'idvc=',idvc,'itrun=',itrun,
     &   'gen_coord_hybrid=',gen_coord_hybrid,'pdryini=',pdryini
        if(lonf .ne. lonfi .or. latg .ne. latgi .or.
     &     jcap .ne. jcapi .or. levs .ne. levsi) then
          print *,' Input resolution and the model resolutions are'
     &,  ' different- run aborted'
          call mpi_quit(660)
        endif
        if ( gen_coord_hybrid ) then
          print *, ' Use generalized hybrid coordinate'
          if (idvc == 3 ) then
           print *, ' Restart input is consistent, run continues'
          else
           print *, ' Restart input is different, run aborted'
           call mpi_quit(661)
          endif
        endif
        if ( hybrid ) then
          print *, ' Use sigma-p hybrid coordinate'
          if (idvc == 2 ) then
           print *, ' Restart input is consistent, run continues'
          else
           print *, ' Restart input is different, run aborted'
           call mpi_quit(662)
          endif
        endif


      endif
!
      allocate (vcoord4(levsi+1,3,2))
      if(.not.allocated(vcoord)) allocate (vcoord(levsi+1,3))
      call nemsio_getfilehead(gfile_in,iret=iret,vcoord=vcoord4)
!
      vcoord(:,1:3) = vcoord4(:,1:3,1)
!      if (me .eq. 0) print *,' vcoord1=',vcoord(1:,1:3)
      deallocate (vcoord4)
!
!---  for generalized tracers
! retrieve nreci, recnamei, reclevtypi, and reclevi
      call nemsio_getfilehead(gfile_in,iret=iret,nrec=nreci)
      if (me == 0) then
        print *, 'LU_TRC: nreci =', nreci, iret
      endif
!
!---
      if (gen_coord_hybrid) then                                        ! hmhj

!        sfcpress_id  = mod(idvm , 10)
!        thermodyn_id = mod(idvm/10 , 10)
!   ak bk ck in file have the same order as model                       ! hmhj
        do k=1,levp1                                                    ! hmhj
          ak5(k) = vcoord(k,1)/1000.                                    ! hmhj
          bk5(k) = vcoord(k,2)                                          ! hmhj
          ck5(k) = vcoord(k,3)/1000.                                    ! hmhj
        enddo                                                           ! hmhj
        vertcoord_id=0                                                  ! hmhj
        do k=1,levp1                                                    ! hmhj
          if( ck5(k).ne.0.0 ) vertcoord_id=3                            ! hmhj
        enddo
! provide better estimated press                                        ! hmhj
        psurfff = 101.3d0
        if( thermodyn_id.eq.3 ) then                                    ! hmhj
          do k=1,levs                                                   ! hmhj
            thref(k) = 300.0*cpd                                        ! hmhj
            teref(k) = 255.0*cpd                                        ! hmhj
          enddo                                                         ! hmhj
        else                                                            ! hmhj
         do k=1,levp1                                                   ! hmhj
          thref(k) = 300.0                                              ! hmhj
          teref(k) = 255.0                                              ! hmhj
         enddo                                                          ! hmhj
        endif
        ck5p(levp1) = ck5(levp1)                                        ! hmhj
        do k=1,levs                                                     ! hmhj
          ck5p(k) = ck5(k)*(teref(k)/thref(k))**rkapi                   ! hmhj
        enddo
        if( me.eq.0 ) then                                              ! hmhj
          do k=1,levp1                                                  ! hmhj
            pressk=ak5(k)+bk5(k)*psurfff+ck5p(k)                        ! hmhj
            print 180,k,ak5(k),bk5(k),ck5(k),pressk                     ! hmhj
 180     format('k=',i2,'  ak5=',f13.6,'  bk5=',e13.5,                  ! hmhj
     &            '   ck5=',f13.6,'  closed pressk=',f10.6)             ! hmhj
          enddo                                                         ! hmhj
        endif                                                           ! hmhj
        do k=1,levs
          dbk(k) = bk5(k)-bk5(k+1)
          bkl(k) = (bk5(k)+bk5(k+1))*0.5
          ck(k)  = ak5(k)*bk5(k+1)-ak5(k+1)*bk5(k)
        enddo
        do k=1,levp1                                                    ! hmhj
          si(k) = ak5(k)/psurfff + bk5(k) + ck5p(k)/psurfff             ! hmhj
        enddo                                                           ! hmhj
        do k=1,levs                                                     ! hmhj
          sl(k) = 0.5*(si(k)+si(k+1))                                   ! hmhj
        enddo                                                           ! hmhj

      else if (hybrid .and. idvc .eq. 2) then
!       idsl=slid  !=2,pk=0.5*(p(k+1/2)+p(k-1/2)) check alfa(1)  am_bm
!   ak bk order in "sigma" file is bottom to top !!!!!!!!!!!!!!!!!!
        psurfff = 101.3
        do k=1,levp1
          ak5(k) = vcoord(levp1+1-k,1)/1000.
          bk5(k) = vcoord(levp1+1-k,2)
          pressk = ak5(k) + bk5(k)*psurfff

          if(me.eq.0)print 190,k,ak5(k),bk5(k),pressk
190       format('k=',i2,'  ak5=',E14.6,'  bk5=',e14.6,
     &           '  pressk=',E14.6)

        enddo
        do k=1,levs
          dbk(k) = bk5(k+1)-bk5(k)
          bkl(k) = (bk5(k+1)+bk5(k))*0.5
          ck(k)  = ak5(k+1)*bk5(k)-ak5(k)*bk5(k+1)
          if(me.eq.0)print 200,k,dbk(k),ck(k)
200       format('k=',i2,'  dbk=',f8.6,'  ck=',e13.5)
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
250       format('k=',i2,'  si=',f8.6,'  sl=',e13.5)
        enddo

      else
        print *,' Non compatible Initial state IDVC=',idvc
     &,' iret=',iret
        call MPI_QUIT(670)
      endif
!
      idate(1)    = idate7(4)
      idate(2:3)  = idate7(2:3)
      idate(4)    = idate7(1)
      WAVES       = jcap
      XLAYERS     = levs
      icen        = 7
      ienst       = iens(1)
      iensi       = iens(2)
!
!
      IF (me.eq.0) THEN
        write(0,*)'gfile,in treadeo fhour,idate=',gfilename,fhour,idate
     &, ' idvc=',idvc,' jcap=',jcap 
     &, ' pdryini=',pdryini,'idate=',idate
      ENDIF
!
      fieldsize=im*jm
      allocate (nemsio_data(fieldsize))
!  Read orog
      call nemsio_readrecv(gfile_in,'hgt','sfc',1,nemsio_data,
     &     iret=iret)
      call split2d_rdGRD(nemsio_data,zsg,fieldsize,global_lats_a,
     &     lonsperlat)
!
!  Read ps
      call nemsio_readrecv(gfile_in,'pres','sfc',1,nemsio_data,
     &     iret=iret)
!      print *,'aftr read in pressfc=',maxval(nemsio_data),
!     &    minval(nemsio_data)
      call split2d_rdGRD(nemsio_data,psg,fieldsize,global_lats_a,
     &  lonsperlat)
!      print *,'in treadg_nemsio,psg=',maxval(psg(1:lonf,1:lats_node_a)),
!     &    minval(psg(1:lonf,1:lats_node_a)),'psg(1:5,1)=',psg(1:5,1),
!     & 'psg(1:5,2)=',psg(1:5,2),'psg(1:5,3)=',psg(1:5,3)
!
!  Read u
      do k=1,levs
        call nemsio_readrecv(gfile_in,'ugrd','mid layer',k,nemsio_data,
     &       iret=iret)
        call split2d_rdGRD(nemsio_data,uug(:,:,k),fieldsize,
     &    global_lats_a,lonsperlat)
      enddo
!  Read v
      do k=1,levs
        call nemsio_readrecv(gfile_in,'vgrd','mid layer',k,nemsio_data,
     &    iret=iret)
        call split2d_rdGRD(nemsio_data,vvg(:,:,k),fieldsize,
     &    global_lats_a,lonsperlat)
      enddo
!  Read T   -- this is real temperature
      do k=1,levs
        call nemsio_readrecv(gfile_in,'tmp','mid layer',k,nemsio_data,
     &    iret=iret)
        call split2d_rdGRD(nemsio_data,ttg(:,:,k),fieldsize,
     &    global_lats_a,lonsperlat)
      enddo
!  Read dp  
      do k=1,levs
        call nemsio_readrecv(gfile_in,'dpres','mid layer',k,nemsio_data,
     &    iret=iret)
        call split2d_rdGRD(nemsio_data,dpg(:,:,k),fieldsize,
     &    global_lats_a,lonsperlat)
      enddo
!
!  Initial Tracers with zero
!
!     rqg(:,:,:) = 0.0
      rqg(:,:,:) = 1.0E-20

!! Generalized tracers: 
!! Loop through ntrac to read in met + chem tracers
!*
      do k=1,levh
          call nemsio_readrecv(gfile_in,'tracer',
     &                   'tracer layer',k,nemsio_data,iret=iret)
          if(iret == 0) then
!            if(me==0) print *,'TRACER read in ok -',
!     &                gfs_dyn_tracer%vname(n),k
!      print *,'aftr read in tmp=',maxval(nemsio_data),
!     &    minval(nemsio_data),'k=',k
            call split2d_rdGRD(nemsio_data,rqg(:,:,k),fieldsize,
     &      global_lats_a,lonsperlat)
!      print *,'in treadg_nemsio,k=',k,'rqg=',
!     &    maxval(rqg(1:lonf,1:lats_node_a,k)),
!     &    minval(rqg(1:lonf,1:lats_node_a,k))
          else
            if(me==0) print *,'TRACER: tracer not found in input; ',
     &         'set chem tracer to default values',me,k
          endif
       enddo
!
      call nemsio_close(gfile_in,iret)
      call nemsio_finalize()
!
      deallocate(nemsio_data)
      iprint=0
 
!!!!
      RETURN
      END
