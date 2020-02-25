      SUBROUTINE read_vcoord(nread,cfile)
!
      USE machine,        only : kind_io4
      use coordinate_def, only : ak5,bk5,ck5
      use layout1,        ONLY: me
      use nemsio_module
      use sigio_module
      use sigio_r_module
!
      implicit none
      integer       nread
      character*(*) cfile
!      nft=152
      integer :: idvc, iret, levp1, k
      call sigio_rropen(nread,cfile,iret)
      call sigio_alhead(nread,iret)
      call sigio_rrhead(nread,head,iret)
      idvc = head%idvc  !idvc=3:sigma-theta and/or p, idvc=2:sigma-p, idvc=1:sigma file
      levp1 =head%levs+1    
        do k=1,levp1					
          ak5(k) = head%vcoord(k,1)/1000.
          bk5(k) = head%vcoord(k,2)				
          ck5(k) = head%vcoord(k,3)/1000.			
        enddo
      call sigio_rclose(nread,iret)
!
      end subroutine read_vcoord_sigvay

      SUBROUTINE read_vcoord(nread,cfile)
!
!***********************************************************************
!
! !REVISION HISTORY:
!
! 2016/03/01 Jun Wang read in vertical coordinate info from grid point 
!                     model level file
!
      USE machine,        only : kind_io4
      use coordinate_def, only : ak5,bk5,ck5
      use layout1,        ONLY: me
      use nemsio_module
!
      integer       nread
      character*(*) cfile
!
      integer idvc,k,dimz,iret
      logical gen_coord_hybrid, hybrid
      type(nemsio_gfile) gfile_in
      real(kind=kind_io4),allocatable :: vcoord(:,:,:)
!
!!read grid point model level file
      if (me == 0) print *,' nread=',nread,' cfile=',cfile

      call nemsio_init()
!
      call nemsio_open(gfile_in,trim(cfile),'read',iret=iret)
!
      call nemsio_getfilehead(gfile_in,iret=iret,
     &     dimz=dimz,idvc=idvc)
!
      gen_coord_hybrid=.false.
      call nemsio_getheadvar(gfile_in,'GEN_COOR', 
     &     gen_coord_hybrid, iret=iret)
      hybrid=.false.
      call nemsio_getheadvar(gfile_in,'HYBRID',hybrid,iret=iret)  
!      print *,' in phys, read vcoord,hybrid=',hybrid,'dimz=',
!     &  dimz,'idvc=',idvc
      if (iret /= 0) then
        if (idvc==2) hybrid=.true.
      endif
!
      allocate(vcoord(dimz+1,3,2))
      call nemsio_getfilehead(gfile_in,iret=iret,vcoord=vcoord)
!
      call nemsio_close(gfile_in,iret)
! 
      if (gen_coord_hybrid) then
!   ak bk ck in file have the same order as model
        do k=1,dimz+1
          ak5(k) = vcoord(k,1,1)/1000.
          bk5(k) = vcoord(k,2,1)
          ck5(k) = vcoord(k,3,1)/1000.
        enddo
      else if (hybrid .and. idvc .eq. 2) then
        ck5=0.
        do k=1,dimz+1
          ak5(k) = vcoord(dimz+2-k,1,1)/1000.
          bk5(k) = vcoord(dimz+2-k,2,1)
        enddo
      endif
!
      deallocate(vcoord)

      end subroutine read_vcoord
!
