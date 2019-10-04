       SUBROUTINE PARA_NST_W(IOPROC,nst_fld, cfile,xhour,idate,
     &         lats_nodes_r,global_lats_r,lonsperlar,
     &         ens_nam)
!!
!! revision history
!  Dec      2010   Jun Wang  change to nemsio library
!
      use resol_def,     ONLY: latr, lonr, levs
      use layout1,       ONLY: lats_node_r, me,ipt_lats_node_r,nodes
      use nemsio_module
      use gfs_physics_nst_var_mod, ONLY: NST_Var_Data
      USE machine,   ONLY: kind_io4, kind_ior, kind_io8,kind_phys
      implicit none
!!
      TYPE(NST_Var_Data),intent(in)    :: nst_fld
      integer,intent(in)               :: idate(4),ioproc
      real(kind=kind_io8),intent(in)   :: xhour
      character*(*),intent(in)         :: cfile,ens_nam
      INTEGER,intent(in)               :: lats_nodes_r(nodes)
      INTEGER,intent(in)               :: GLOBAL_LATS_R(latr)
      INTEGER,intent(in)               :: lonsperlar(latr)
!
!!
      real(kind=kind_io8),allocatable:: bfo(:,:)
      integer k,lan,i,nphyfld,fieldsize
!!
      integer,save:: version
!
      type(nemsio_gfile)  :: gfile
!
      integer iret,nrec,nmetavari,nmetavarr8,nmetaaryi,nmetaaryr8,
     &    idate7(7),nmeta,jrec,l,nsfcrec,nnstrec
      integer nfhour,nfminute,nfsecondn,nfsecondd
      character(16),allocatable :: variname(:),varr8name(:)
!      character(16),allocatable :: aryiname(:),aryr8name(:)
!      integer,allocatable :: varival(:),aryilen(:),aryr8len(:),
!     &     aryival(:,:)
      real(kind=kind_io8),allocatable :: varr8val(:),aryr8val(:,:)
      character(16),allocatable :: recname(:),reclevtyp(:)
      integer,allocatable :: reclev(:)
      logical first
      sAve first,nmetavarr8,varr8name,varr8val,
     &     nmeta,nrec,nnstrec,recname,reclevtyp,reclev
!     &     nmetaaryi,nmetavari,variname,varival,aryiname,aryilen,
!     &     aryival,nmetaaryr8,aryr8name,aryr8len,aryr8val,
      data first /.true./
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!      print *,'in para_nst_w'
      nnstrec=19
      fieldsize=sum(lonsperlar)
      allocate(bfo(fieldsize,nnstrec))
!
      call fld_collect_nst(nst_fld,
     &  fieldsize,nnstrec,bfo,lats_nodes_r,global_lats_r,lonsperlar,
     &  ioproc,nnstrec)
!      write(0,*)'after fld_collect'
!
      if (me.eq.ioproc) then
        if (first) then
          nmeta=5
!
          nmetavarr8=1
          allocate(varr8name(nmetavarr8),varr8val(nmetavarr8))
          varr8name=(/'fhour'/)
!
          nrec=nnstrec
          allocate(recname(nrec),reclevtyp(nrec),reclev(nrec))
!record name
          RECNAME(1)='slmsk'
          RECNAME(2)='xt'
          RECNAME(3)='xs'
          RECNAME(4)='xu'
          RECNAME(5)='xv'
          RECNAME(6)='xz'
          RECNAME(7)='zm'
          RECNAME(8)='xtts'
          RECNAME(9)='xzts'
          RECNAME(10)='dtcool'
          RECNAME(11)='zc'
          RECNAME(12)='c0'
          RECNAME(13)='cd'
          RECNAME(14)='w0'
          RECNAME(15)='wd'
          RECNAME(16)='dconv'
          RECNAME(17)='ifd'
          RECNAME(18)='tref'
          RECNAME(19)='qrain'
!
          RECLEVTYP(1:nrec)='sfc'
          RECLEV(1:nrec)=1
!
!endif first 
        endif
!
        varr8val(1)=xhour
!
        idate7(1:6)=0;idate7(7)=1
        idate7(1)=idate(4)
        idate7(2:3)=idate(2:3)
        idate7(4)=idate(1)
        nfhour=int(xhour)
        nfminute=int((xhour-nfhour)*60.)
        nfsecondn=int((xhour-nfhour)*3600.-nfminute*60.)
        nfsecondd=1
!
        PRINT 99,xhour,IDATE
99      FORMAT(1H ,'in fixio HOUR=',f8.2,3x,'IDATE=',
     &  4(1X,I4))
!!
! open nemsio sfc restart file
        call nemsio_init()
!
        write(0,*)'before nemsio_open for restart nst file'
        call nemsio_open(gfile,trim(cfile),'write',iret,   
     & modelname='GFS',gdatatype='bin8',idate=idate7,nfhour=nfhour, 
     & nfminute=nfminute,nfsecondn=nfsecondn,nfsecondd=nfsecondd,   
     & dimx=fieldsize,dimy=1,dimz=levs,nrec=nrec,         
     & nmeta=5,recname=recname,reclevtyp=reclevtyp,reclev=reclev,   
     & extrameta=.true., nmetavarr8=nmetavarr8,
     & varr8name=varr8name,varr8val=varr8val)         

!       print *,'after restart nst nemsio_open, iret=',iret

        do jrec=1,nrec
          call nemsio_writerec(gfile,jrec,bfo(:,jrec),iret=iret)
        enddo
        deallocate(bfo)
!
        call nemsio_close(gfile)
        call nemsio_finalize()
!
! endof ioproc
      ENDIF
!
      if(first) first=.false.
         
      return
      end
