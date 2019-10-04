      subroutine treadeo_nemsio_iau(cfile,IDATE
     &,                             psg,uug,vvg,ttg,rqg 
     &,                             LS_NODE,LS_NODES,MAX_LS_NODES
     &,                             IPRINT
     &,                             global_lats_a,lats_nodes_a
     &,                             lonsperlat)

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
!  
 
      use gfs_dyn_resol_def
      use gfs_dyn_layout1
      use gfs_dyn_coordinate_def
      use namelist_dynamics_def
      use gfs_dyn_vert_def
      use gfs_dyn_mpi_def
      use gfs_dyn_physcons, rkap => con_rocp
     &,             cpd => con_cp
      use nemsio_module
      use nemsio_def
      use gfs_dyn_tracer_config, ONLY : gfs_dyn_tracer   ! generalized tracer
      use gfs_dyn_io_header

!
      IMPLICIT NONE
      character*(*) cfile
      INTEGER             IDATE(4), idate7(7)
     &,                   levsi, jcapi, latgi, lonfi
      logical       fromiau 
      integer              ls_node(ls_dim,3)
!
      INTEGER, dimension(nodes) :: MAX_LS_NODES, lats_nodes_a
      INTEGER              LS_NODES(LS_DIM,NODES)
      INTEGER              J,K,L,LOCL,N,lv,kk,w3rlkind,w3ikind,iprint
     &,                    i,lan,lat,iblk,lons_lat,il,lon,njeff,nn
     &,                    indev,indod,indev1,indev2,indod1,indod2
     &,                    nfhour,nfminute,nfsecondn,nfsecondd 
!

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
!
      real (kind=kind_io4), allocatable ::  nemsio_data(:)
!!
      real(kind=kind_grid), dimension(lonf,lats_node_a) :: psg
      real(kind=kind_grid), dimension(lonf,lats_node_a,levs) :: uug
     &,                       vvg, ttg
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
      call nemsio_getheadvar(gfile_in,'iens',iens,iret=iret)
      call nemsio_getheadvar(gfile_in,'idpp',idpp,iret=iret)
      call nemsio_getheadvar(gfile_in,'idrun',idrun,iret=iret)
      call nemsio_getheadvar(gfile_in,'idusr',idusr,iret=iret)

      if (me == 0) then
        write(0,*)'iret=',iret,
     &     ' levsi=',levsi,' idate=',idate,
     &   'lonf=',lonf,'lonfi=',lonfi,'latg=',latg,'latgi=',latgi,
     &   'jcap=',jcap,'jcapi=',jcapi,'levs=',levs,'levsi=',levsi,
     &   'idvc=',idvc,'tlmeta=',tlmeta,
     &   'gen_coord_hybrid=',gen_coord_hybrid 
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
      igen        = igen
      ienst       = iens(1)
      iensi       = iens(2)
!Weiyu, pdryini4 did not initiaqlized, temp give value zero.
!-----------------------------------------------------------
!
!
      IF (me.eq.0) THEN
        write(0,*)'cfile,in treadeo idate=',cfile,idate
     &, ' idvc=',idvc,' jcap=',jcap
      ENDIF
!
      allocate (nemsio_data(lonf*latg))
       stime=timef()

      call w3kind(w3rlkind,w3ikind)
!
!  Read ps
      if(w3rlkind==8) then
        call nemsio_readrecv(gfile_in,'pres','sfc',1,nemsio_data,
     &    iret=iret)
      elseif(w3rlkind==4) then
        call nemsio_readrecvw34(gfile_in,'pres','sfc',1,nemsio_data,
     &    iret=iret)
      endif
      call split2d(nemsio_data,buffo,global_lats_a)
      CALL interpred(1,kmsk,buffo,psg,global_lats_a,lonsperlat)
!
!  Read u
      do k=1,levs
       if(w3rlkind==8) then
         call nemsio_readrecv(gfile_in,'ugrd','mid layer',k,nemsio_data,
     &     iret=iret)
       elseif(w3rlkind==4) then
         call nemsio_readrecvw34(gfile_in,'ugrd','mid layer',k,
     &     nemsio_data,iret=iret)
       endif
!      print *,'in treadeo,ugrd=',maxval(nemsio_data),
!     & minval(nemsio_data),'iret=',iret,'k=',k
       call split2d(nemsio_data,buffo,global_lats_a)
       CALL interpred(1,kmsk,buffo,uug(1,1,k),global_lats_a,lonsperlat)
      enddo
!  Read v
      do k=1,levs
       if(w3rlkind==8) then
        call nemsio_readrecv(gfile_in,'vgrd','mid layer',k,nemsio_data,
     &   iret=iret)
       elseif(w3rlkind==4) then
        call nemsio_readrecvw34(gfile_in,'vgrd','mid layer',k,
     &   nemsio_data,iret=iret)
       endif
       call split2d(nemsio_data,buffo,global_lats_a)
       CALL interpred(1,kmsk,buffo,vvg(1,1,k),global_lats_a,lonsperlat)
      enddo
!  Read T   -- this is real temperature
      do k=1,levs
       if(w3rlkind==8) then
        call nemsio_readrecv(gfile_in,'tmp','mid layer',k,nemsio_data,
     &     iret=iret)
       elseif(w3rlkind==4) then
        call nemsio_readrecvw34(gfile_in,'tmp','mid layer',k,
     &     nemsio_data,iret=iret)
       endif
        call split2d(nemsio_data,buffo,global_lats_a)
        CALL interpred(1,kmsk,buffo,ttg(1,1,k),global_lats_a,lonsperlat)
      enddo

!  Initial Tracers with zero
!
!*    rqg(:,:,:) = 0.0
      rqg(:,:,:) = 1.e-20

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
     &                       'mid layer',k,nemsio_data,iret=iret)
         elseif(w3rlkind==4) then
          call nemsio_readrecvw34(gfile_in,trim(vname),
     &                       'mid layer',k,nemsio_data,iret=iret)
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
      call nemsio_close(gfile_in,iret)
      iprint = 0
!!!!
!!!! Added by PJP
      deallocate (recnamei)
      deallocate (reclevtypi)
      deallocate (reclevi)
      deallocate (nemsio_data)
      RETURN
      END
