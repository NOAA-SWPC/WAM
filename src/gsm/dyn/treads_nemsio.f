      SUBROUTINE TREADS_nemsio(sfile,FHOUR,IDATE,
     &                         trie_ls,trio_ls,
     &                         LS_NODE)
!!
!-------------------------------------------------------------------
!*** program log
!*** Dec, 2009 Jun Wang:  read spectral variables for restart
!*** Dec, 2010 Jun Wang:  change to nemsio library
!*** Feb, 2011 Henry Juang: change argument to be generic one to fit NDSL
!*** Mar, 2013 Jun Wang:  save orog from restart
!-------------------------------------------------------------------
!
!  
      use gfs_dyn_resol_def
      use gfs_dyn_layout1
      use gfs_dyn_coordinate_def					! hmhj
      use gfs_dyn_io_header
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
      character*(*) sfile
      REAL(KIND=KIND_EVOD) FHOUR
      INTEGER              IDATE(4),NTRACI, nreci,jcapi,idate7(7),levsi
!!
      REAL(KIND=KIND_EVOD) TRIE_LS(LEN_TRIE_LS,2,LOTLS)
      REAL(KIND=KIND_EVOD) TRIO_LS(LEN_TRIO_LS,2,LOTLS)
!
      integer              ls_node(ls_dim,3)
!
!cmr  ls_node(1,1) ... ls_node(ls_max_node,1) : values of L
!cmr  ls_node(1,2) ... ls_node(ls_max_node,2) : values of jbasev
!cmr  ls_node(1,3) ... ls_node(ls_max_node,3) : values of jbasod
!
      INTEGER              J,K,L,LOCL,N,lv,kk
      integer              i,lan,lat,lons_lat,il,lon,njeff,nn
      integer              indev
      integer              indod
      integer              indev1,indev2
      integer              indod1,indod2
! for generalized tracers
!      character*16, allocatable  :: recnamei(:)
!      character*16, allocatable  :: reclevtypi(:)
!      integer,      allocatable  :: reclevi(:)
!
      REAL(KIND=KIND_IO4)  fhour4
      type (nemsio_gfile)  gfile_in
      integer              iret,num_dta,im,jm,fieldsize
!
      real (kind=kind_io8), allocatable ::  nemsio_sdata(:)
!!
!------------------------------------------------------------------
      print *,'in treads_nemsio,sfile=',trim(sfile)
!
!--- Input file is in grid-point space - use nems_io package
!
      call nemsio_init()
!
      call nemsio_open(gfile_in,trim(sfile),'read',iret)
!
      call nemsio_getfilehead(gfile_in,iret=iret,
     &  idate=idate7,dimx=im,dimy=jm,dimz=levsi,jcap=jcapi,
     &  nrec=nreci)
!
      call nemsio_getheadvar(gfile_in,'fhour',fhour,iret=iret)
!
      if(me==0) print *,'idate=',idate,'im=',im,'jm=',jm,
     &  'levsi=',levsi,'jcap=',jcap,'nrec=',nreci,'fhour=',fhour
!
!---  for generalized tracers
! retrieve nreci, recnamei, reclevtypi, and reclevi
!
!      allocate (recnamei(nreci))
!      allocate (reclevtypi(nreci))
!      allocate (reclevi(nreci))
!      call nemsio_getfilehead(gfile_in,iret=iret,recname=recnamei,
!     &                       reclevtyp=reclevtypi,reclev=reclevi)
!
!---
      fieldsize=im*jm
      allocate (nemsio_sdata(fieldsize))
!  Read orog
      call nemsio_readrecv(gfile_in,'gz','sfc',1,nemsio_sdata,iret=iret)
      z=nemsio_sdata
!      print *,'in treads,read gze=',maxval(nemsio_sdata),
!     &   minval(nemsio_sdata)
      call triseori(nemsio_sdata,                                      
     &              trie_ls(1,1,p_gz),trio_ls(1,1,p_gz),1,ls_node)
!      print *,'in treads,gze=',maxval(gze(:,1)),minval(gze(:,1)),
!     & maxval(gze(:,2)),minval(gze(:,2)),'gzo=',maxval(gzo(:,1)),
!     & minval(gzo(:,1)),maxval(gzo(:,2)),minval(gzo(:,2)),'gze(1:5,2)='
!     &,gze(1:5,1),gze(1:5,2),'gzo(1:5,2=',gzo(1:5,1),gzo(1:5,2)
!
!  Read ps
      call nemsio_readrecv(gfile_in,'pres','sfc',1,nemsio_sdata,
     &   iret=iret)
      call triseori(nemsio_sdata,
     &              trie_ls(1,1,p_q),trio_ls(1,1,p_q),1,ls_node)
!      print *,'in treads,qe=',maxval(qe),minval(qe),'qo=',
!     &  maxval(qo),minval(qo),'qe(1:5,1:2)=',qe(1:5,1),qe(1:5,2),
!     &  'qo(1:5,1:2)=',qo(1:5,1),qo(1:5,2)
!
!  Read dp
      do k=1,levs
        kk = p_dp + k - 1
        call nemsio_readrecv(gfile_in,'dpres','mid layer',k,
     &       nemsio_sdata,iret=iret)
        call triseori(nemsio_sdata,
     &                trie_ls(1,1,kk),trio_ls(1,1,kk),1,ls_node)
      enddo
!
!  Read t
      do k=1,levs
        kk = p_te + k - 1
        call nemsio_readrecv(gfile_in,'tmp','mid layer',k,nemsio_sdata,
     &    iret=iret)
!      print *,'in treads,read tmp=',maxval(nemsio_sdata),
!     &   minval(nemsio_sdata),'k=',k
        call triseori(nemsio_sdata,   
     &                trie_ls(1,1,kk),trio_ls(1,1,kk),1,ls_node)
!      print *,'in treads,',k,'tee=',maxval(tee(:,:,k)),
!     &  minval(tee(:,:,k)),'teo=',maxval(teo(:,:,k)),minval(teo(:,:,k))
      enddo
!
!  Read di
      do k=1,levs
        kk = p_di + k - 1
        call nemsio_readrecv(gfile_in,'di','mid layer',k,nemsio_sdata,
     &    iret=iret)
        call triseori(nemsio_sdata,  
     &                trie_ls(1,1,kk),trio_ls(1,1,kk),1,ls_node)
      enddo
!  Read ze 
      do k=1,levs
        kk = p_ze + k - 1
        call nemsio_readrecv(gfile_in,'ze','mid layer',k,nemsio_sdata,
     &    iret=iret)
        call triseori(nemsio_sdata, 
     &                trie_ls(1,1,kk),trio_ls(1,1,kk),1,ls_node)
      enddo
!

      if( .not. ndslfv ) then

      trie_ls(:,:,p_rq:p_rq+levh-1) = 0.0
      trio_ls(:,:,p_rq:p_rq+levh-1) = 0.0

!! Generalized tracers: 
!! Loop through ntrac to read in met + chem tracers
!*
      do k=1,levh
        kk = p_rq + k - 1
        call nemsio_readrecv(gfile_in,'rq',
     &     'tracer layer',k,nemsio_sdata,iret=iret)
        if(iret == 0) then
          call triseori(nemsio_sdata,   
     &                  trie_ls(1,1,kk),trio_ls(1,1,kk),1,ls_node)
        else
          if(me==0) print *,'TRACER: tracer not found in input; ',
     &         'set chem tracer to default values',me,k
        endif
      enddo       
!
      endif
!
      call nemsio_close(gfile_in,iret)
      call nemsio_finalize()
!
      deallocate(nemsio_sdata)
!      print *,'endof treads_nemsio'
!!!!
      RETURN
      END
