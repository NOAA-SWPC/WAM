      subroutine twriteg_rst(fname,IOPROC,fhour,idate,
     x           si,pdryini,global_lats_a,lonsperlat,lats_nodes_a,
     x           psg,dpg,ttg,uug,vvg,rqg,zsg,
     &           kdt,nfcstdate7)
!
!-------------------------------------------------------------------
!*** program log
!*** Dec, 2009 Jun Wang:  write spectral variables for restart
!*** Dec, 2010 Jun Wang:  change to nemsio library
!*** Feb, 2011 Henry Juang: add option for mass_dp and  NDSL
!*** Nov, 2011 Jun Wang:  remove nvcoord in restart file (not used)
!*** May, 2013 Shrinivas Moorthi: correct definition of idvt for WAM
!-------------------------------------------------------------------
!
      use gfs_dyn_resol_def
      use gfs_dyn_layout1
      use gfs_dyn_io_header
      use gfs_dyn_coordinate_def	
      use namelist_dynamics_def
      use gfs_dyn_mpi_def
      use nemsio_module
!
      implicit none
!
      character(*),intent(in) :: fname
      integer,intent(in) :: ioproc
      real(kind=kind_evod),intent(in) :: fhour
      integer,intent(in) :: idate(4)
      integer,intent(in) :: kdt,nfcstdate7(7)
!
      real(kind=kind_evod) si(levp1)
      REAL(KIND=KIND_IO8),intent(in) :: pdryini
      integer,intent(in) ::          global_lats_a(latg)
      integer,intent(in) ::          lonsperlat(latg)
      integer,intent(in) ::          lats_nodes_a(nodes)
!
      REAL(KIND=KIND_GRID),intent(in) :: zsg(lonf,lats_node_a_max)
      REAL(KIND=KIND_GRID),intent(in) :: psg(lonf,lats_node_a_max)
      REAL(KIND=KIND_GRID),intent(in) :: dpg(lonf,lats_node_a_max,levs)
      REAL(KIND=KIND_GRID),intent(in) :: ttg(lonf,lats_node_a_max,levs)
      REAL(KIND=KIND_GRID),intent(in) :: uug(lonf,lats_node_a_max,levs)
      REAL(KIND=KIND_GRID),intent(in) :: vvg(lonf,lats_node_a_max,levs)
      REAL(KIND=KIND_GRID),intent(in) :: rqg(lonf,lats_node_a_max,levh)
!
!local variables:
      REAL(kind=8) t1,t2,t3,t4,t5,t6,ta,tb,rtc
!
      integer              ierr,j,k,l,lenrec,locl,n,node,jrec,idate7(7)
!
      real(kind=kind_ior),allocatable ::  tmp8(:,:),buf(:)
      real(kind=kind_ior),allocatable ::  GZ(:)
!
      type(nemsio_gfile) gfile
      integer kps,ktt,kuu,kvv,krq,kdp
      integer nfhour,nfminute,nfsecondn,nfsecondd,nrec,nmeta
      integer nmetavari,nmetavarr8,nmetaaryi,fieldsize
      character(16),allocatable :: recname(:),reclevtyp(:)
      character(16),allocatable :: variname(:),varr8name(:),aryiname(:)
      integer,allocatable :: reclev(:)
      integer,allocatable :: varival(:),aryilen(:),aryival(:,:)
      real(kind=kind_io4),allocatable :: vcoord4(:,:,:)
      real(8),allocatable :: varr8val(:)
      integer iret, ipt_lats
      integer  il,ilen,i,msgtag,ls_diml,nodesl,ij
      logical first
      save first,nmetavari,nmetavarr8,nmetaaryi,recname,reclevtyp, 
     &     reclev,variname,varr8name,aryiname,varival,varr8val,
     &     aryilen,aryival,vcoord4,nmeta,nrec,GZ
      data first /.true./
!
!
      real(kind=kind_mpi_r),allocatable :: grid_node (:,:,:)
      real(kind=kind_mpi_r),allocatable :: grid_nodes(:,:,:,:)
!
      real(kind=kind_mpi_r),allocatable :: grid_gz_nodes(:,:,:)
!
      integer      lan,lat,iblk,lons_lat,lon,NJEFF,nn,lv,lotg_rst
!
!---------------------------------------------------------------------
!
!       print *,' enter twriteg_rst ' 
!       print *,'lonf=',lonf,'lats_node_a_max=',lats_node_a_max,
!     &  'total_levels=',3*levs+1*levh+1,'lonsperlat=',lonsperlat

      lotg_rst = 4*levs+1*levh+1	! ps, dp, tt, uu, vv, rq

      allocate ( grid_node ( lonf,lats_node_a_max,lotg_rst ) )
!
      fieldsize=sum(lonsperlat)
!      print *,'fieldsize=',fieldsize
!
!collect data 
      kps=1
      kdp=kps+1
      ktt=kdp+levs
      kuu=ktt+levs
      kvv=kuu+levs
      krq=kvv+levs
!
      do j=1,lats_node_a_max
        do i=1,lonf
          grid_node(i,j,kps) = psg(i,j)
        enddo
      enddo
!
      do k=1,levs
        do j=1,lats_node_a_max
          do i=1,lonf
            grid_node(i,j,kdp+k-1) = dpg(i,j,k)
            grid_node(i,j,ktt+k-1) = ttg(i,j,k)
            grid_node(i,j,kuu+k-1) = uug(i,j,k)
            grid_node(i,j,kvv+k-1) = vvg(i,j,k)
          enddo
        enddo
      enddo
!
      do k=1,levh
        do j=1,lats_node_a_max
          do i=1,lonf
            grid_node(i,j,krq+k-1) = rqg(i,j,k)
          enddo
        enddo
      enddo
!
!      print *,'pdryini=',pdryini,'lats_nodes_a=',lats_nodes_a, 
!     &  'nodes=',nodes
!
!WY bug fix.
!-----------
      if ( me .eq. ioproc ) then
         allocate ( grid_nodes ( lonf,lats_node_a_max,
     &                           lotg_rst, nodes ),stat=ierr )
         if(first) then
           allocate ( grid_gz_nodes (lonf,lats_node_a_max,nodes),
     &                stat=ierr )
         endif
      else
         allocate (grid_nodes(1, 1, 1, 1), stat = ierr)
         if(first) then
           allocate (grid_gz_nodes(1, 1, 1), stat = ierr )
         endif
      endif
      if (ierr .ne. 0) then
        call mpi_abort(mpi_comm_all,ierr,i)
      endif
!
      if(nodes>1) then
        lenrec = lonf*lats_node_a_max * lotg_rst
!
        t1=rtc()
!        print *,'after allocate grid_nodes,lenrec=',lenrec
        call mpi_gather( grid_node , lenrec, MPI_R_MPI_R,
     x                 grid_nodes, lenrec, MPI_R_MPI_R,
     x                 ioproc, MPI_COMM_ALL, ierr)
      else
        grid_nodes(:,:,:,1)=grid_node(:,:,:)
      endif
      deallocate(grid_node)
      if(first) then
        if(nodes>1) then
          lenrec=lonf*lats_node_a_max 
          call mpi_gather( zsg , lenrec, MPI_R_MPI_R,
     x                 grid_gz_nodes, lenrec, MPI_R_MPI_R,
     x                 ioproc, MPI_COMM_ALL, ierr)
        else
          grid_gz_nodes(:,:,1)=zsg(:,:)
        endif
        if ( me .eq. ioproc ) then
!
        allocate ( tmp8 ( lonf,latg ) )
        allocate ( gz(fieldsize) )
        ipt_lats=1
        do node=1,nodes
           do j=1,lats_nodes_a(node)
            lat=global_lats_a(ipt_lats-1+j)
            do i=1,lonf
              tmp8(i,lat)=grid_gz_nodes(i,j,node)
            enddo
          enddo
          ipt_lats=ipt_lats+lats_nodes_a(node)
        enddo
        ij=0
        do j=1,latg
            do i=1,lonsperlat(j)
              ij=ij+1
              GZ(ij) = tmp8(i,j)
            enddo
        enddo
        deallocate(tmp8)
!WY bug fix.
!-----------
!        deallocate ( grid_gz_nodes )
!         print *,'gz=',maxval(gz),minval(gz)
       endif
       deallocate ( grid_gz_nodes )
      endif
!
      t2=rtc()
!
      IF (me.eq.ioproc) THEN
 
!        print *,' in TWRITEG fhour=',fhour,'ij=',ij,'fieldsize=',
!     &    fieldsize
!
        if (first) then
!
          nmeta=6
          nrec=lotg_rst+1
          allocate(vcoord4(levp1,3,2))
          vcoord4 = 0.0
          if (lsidea) then
            idvt  = 200
          else
            idvt  = (ntoz-1) + 10 * (ntcw-1)
          endif

          nmetavari=15
          allocate(variname(nmetavari),varival(nmetavari))
          variname(1:nmetavari)=(/"latb     ","lonb     ","itrun    ", 
     &   "iorder   ","irealf   ","igen     ","latf     ","lonf     ",
     &   "icen2    ","idpp     ","idvt     ","idrun    ","idusr    ",
     &   "ivs      ","NTIMESTEP"/)
          varival(1:nmetavari-1)=(/latb,lonb,itrun,2,2,igen,
     &    latg,lonf,icen2,0,idvt,0,0,200509 /)

          nmetaaryi=2
          allocate(aryiname(nmetaaryi),aryilen(nmetaaryi))
          aryiname(1)='iens'
          aryilen(1)=2
          aryiname(2)='FCSTDATE'
          aryilen(2)=7
          allocate(aryival(maxval(aryilen),nmetaaryi))
          aryival(1,1)=ienst
          aryival(2,1)=iensi
!!
          if (gen_coord_hybrid) then					! hmhj

!           idvc    = vctype						! hmhj
            idvc    = 3      					! hmhj
            idvm    = 32    ! 1: ln(ps) 2:ps				! hmhj
            idsl    = 2    ! idsl=2 for middle of layer		! hmhj
            do k=1,levp1							! hmhj
              vcoord4(k,1,1)=ak5(k)*1000.				! hmhj
              vcoord4(k,2,1)=bk5(k)					! hmhj
              vcoord4(k,3,1)=ck5(k)*1000.				! hmhj
            enddo								! hmhj

          else if ( hybrid ) then
            idvc    = 2    ! for hybrid vertical coord.
            do k=1,levp1
              vcoord4(k,1,1)=ak5(levp1+1-k)*1000.
              vcoord4(k,2,1)=bk5(levp1+1-k)
!           if(me.eq.0)print 190,k,head%vcoord(k,1),head%vcoord(k,2)
190         format('in twrite k=',i2,'  ak5r4=',f13.6,'  bk5r4=',e13.5)
            enddo
          else
            idvc    = 1    ! for sigma vertical coord. (default)
            vcoord4(:,1,1) = si (:)
          endif
!
!-- field infomation
          allocate(recname(nrec),reclevtyp(nrec),reclev(nrec))
          recname(1)='hgt'
          recname(2)='pres'
          recname(       3:  levs+2)='dpres'
          recname(  levs+3:2*levs+2)='tmp'
          recname(2*levs+3:3*levs+2)='ugrd'
          recname(3*levs+3:4*levs+2)='vgrd'
          recname(4*levs+3:4*levs+levh+2)='tracer'
          reclevtyp(1)='sfc'
          reclevtyp(2)='sfc'
          reclevtyp(3:4*levs+2)='mid layer'
          reclevtyp(4*levs+3:4*levs+levh+2)='tracer layer'
          reclev(1)=1
          reclev(2)=1
          do i=1,levs
            reclev(i+2)=i
            reclev(i+2+  levs)=i
            reclev(i+2+2*levs)=i
            reclev(i+2+3*levs)=i
          enddo
          do i=1,levh
            reclev(i+2+4*levs)=i
          enddo
!
          nmetavarr8=2
          allocate(varr8name(nmetavarr8),varr8val(nmetavarr8))
          varr8name(1:nmetavarr8)=(/'pdryini','fhour  '/)
!
!endof first
        endif
!
        idate7=0;idate7(7)=1
        idate7(1)=idate(4)
        idate7(2:3)=idate(2:3)
        idate7(4)=idate(1)
!
        nfhour=int(fhour)
        nfminute=int((fhour-nfhour)*60)
        nfsecondn=int(((fhour-nfhour)*60-nfminute)*60)
        nfsecondd=1
!
        varr8val(1)=pdryini
        varr8val(2)=fhour
        varival(nmetavari)=kdt
        aryival(1:7,2)=nfcstdate7(1:7)
!
        call nemsio_open(gfile,fname,'write',iret,modelname='GFS',   
     &    gdatatype='bin8',idate=idate7,nfhour=nfhour,nfminute=nfminute,
     &    nfsecondn=nfsecondn,nfsecondd=nfsecondd,dimx=fieldsize,dimy=1,
     &    dimz=levs,nmeta=nmeta,jcap=jcap,idsl=idsl,
     &    idvm=idvm,idvc=idvc,ntrac=ntrac,nrec=nrec,ncldt=ncld,
     &    recname=recname,reclevtyp=reclevtyp,reclev=reclev,
     &    vcoord=vcoord4,extrameta=.true.,nmetavari=nmetavari,
     &    nmetavarr8=nmetavarr8,nmetaaryi=nmetaaryi,
     &    variname=variname,varival=varival,varr8name=varr8name,
     &    varr8val=varr8val,aryiname=aryiname,aryilen=aryilen,
     &    aryival=aryival)
!
!--- write out data
!
       jrec=1
       call nemsio_writerec(gfile,1,GZ,iret=iret)
!
       allocate ( tmp8 ( lonf,latg ) )
       allocate ( buf(fieldsize) )
       do k=1,lotg_rst
         jrec=k+1
         ipt_lats=1
         do node=1,nodes
           do j=1,lats_nodes_a(node)
            lat=global_lats_a(ipt_lats-1+j)
            do i=1,lonf
              tmp8(i,lat)=grid_nodes(i,j,k,node)
            enddo
           enddo
           ipt_lats=ipt_lats+lats_nodes_a(node)
         enddo
         ij=0
         do j=1,latg
            do i=1,lonsperlat(j)
              ij=ij+1
              buf(ij) = tmp8(i,j)
            enddo
         enddo
         call nemsio_writerec(gfile,jrec,buf,iret=iret)
       end do
       deallocate(buf,tmp8)
!WY bug fix.
!-----------
!       deallocate(grid_nodes)
       call nemsio_close(gfile,iret=iret)
!
        t4=rtc  ()
!sela print *, ' DISK TIME FOR SIG TWRITEO WRT ',t4-t3
!
      endif   !me.eq.ioproc
      deallocate(grid_nodes)
!!
      if(first) then
          first = .false.
      endif
!
      call mpi_barrier(MPI_COMM_ALL,ierr)
!      print *,' leave twrites_nemsio ' 

      return
      end
