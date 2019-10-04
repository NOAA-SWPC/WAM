      subroutine fld_collect_nst(nst_fld,
     &    fldsz,nfld,bfo,lats_nodes_r,global_lats_r,lonsperlar,ioproc,
     &    nrecs)
!
!-------------------------------------------------------------------
!*** program log
!*** Dec, 2009 Jun Wang: ioproc collects partial restart file fields 
!*** Aug, 2010 Jun Wang: changed for nst restart file
!-------------------------------------------------------------------
!
      use resol_def,     ONLY: latr, lonr, levs
      use layout1,       ONLY: lats_node_r, me,ipt_lats_node_r,
     &                         lats_node_r_max,nodes
      use gfs_physics_nst_var_mod, ONLY: NST_Var_Data
      USE machine,   ONLY: kind_io4, kind_ior, kind_io8,kind_phys
      use mpi_def,   only: mpi_comm_all,MPI_R_MPI_R
      implicit none

      TYPE(Nst_Var_Data),intent(in)    :: Nst_fld
      INTEGER,intent(in)               :: lats_nodes_r(nodes)
      INTEGER,intent(in)               :: GLOBAL_LATS_R(latr)
      INTEGER,intent(in)               :: lonsperlar(latr)
!
      integer,intent(in) :: fldsz,nfld,ioproc,nrecs
      real(kind=kind_io8),intent(inout) :: bfo(fldsz,nfld)
!
!local variables:
      REAL(kind=8) t1,t2,t3,t4,t5,t6,ta,tb,rtc
!
      integer      ierr,j,k,l,lenrec,locl,n,node,nf,ipt_lats,nrecs1
!
      real(kind=kind_ior),allocatable ::  tmp8(:,:)
!
      real(kind=kind_ior),allocatable :: grid_node(:,:,:)
      real(kind=kind_ior),allocatable :: grid_nodes(:,:,:,:)
cc
      integer      lan,lat,lons_lat,lon,i,il
cc
!---------------------------------------------------------------------
!
!      Print *,'in fld_collect,nst,lonr=',lonr,'lats_node_r_max=',
!     &   lats_node_r_max,
!     &  'total_levels=',nfld,'lonsperlar=',lonsperlar

      allocate ( grid_node ( lonr,lats_node_r_max,nfld ) )
      grid_node=0.
!
!collect data 
      do j=1,lats_node_r_max
        do i=1,lonr
          grid_node(i,j,1) = nst_fld%slmsk(i,j)
          grid_node(i,j,2) = nst_fld%xt(i,j)
          grid_node(i,j,3) = nst_fld%xs(i,j)
          grid_node(i,j,4) = nst_fld%xu(i,j)
          grid_node(i,j,5) = nst_fld%xv(i,j)
          grid_node(i,j,6) = nst_fld%xz(i,j)
          grid_node(i,j,7) = nst_fld%zm(i,j)
          grid_node(i,j,8) = nst_fld%xtts(i,j)
          grid_node(i,j,9) = nst_fld%xzts(i,j)
          grid_node(i,j,10) = nst_fld%dt_cool(i,j)
          grid_node(i,j,11) = nst_fld%z_c(i,j)
          grid_node(i,j,12) = nst_fld%c_0(i,j)
          grid_node(i,j,13) = nst_fld%c_d(i,j)
          grid_node(i,j,14) = nst_fld%w_0(i,j)
          grid_node(i,j,15) = nst_fld%w_d(i,j)
          grid_node(i,j,16) = nst_fld%d_conv(i,j)
          grid_node(i,j,17) = nst_fld%ifd(i,j)
          grid_node(i,j,18) = nst_fld%tref(i,j)
          grid_node(i,j,19) = nst_fld%qrain(i,j)
        enddo
      enddo
!
!      print *,'in fld_collect,nfld=',nfld
!
!WY bug fix.
!-----------
      if ( me .eq. ioproc ) then
         allocate ( grid_nodes ( lonr,lats_node_r_max,
     &                           nfld, nodes ),stat=ierr )
      else
         allocate(grid_nodes(1, 1, 1, 1), stat = ierr)
      endif
      if (ierr .ne. 0) then
        call mpi_abort(mpi_comm_all,ierr,i)
      endif
!
      lenrec = lonr*lats_node_r_max * nfld
!
      t1=rtc()
      if(nodes>1) then
        call mpi_gather( grid_node , lenrec, MPI_R_MPI_R,
     x                 grid_nodes, lenrec, MPI_R_MPI_R,
     x                 ioproc, MPI_COMM_ALL, ierr)
      elseif(nodes==1) then
        grid_nodes(:,:,:,1)=grid_node(:,:,:)
      else
        print *,'ERROR: the totals nodes=',nodes
        call mpi_abort(mpi_comm_all,ierr,i)
      endif
!
      deallocate(grid_node)
!      print *,'after gather grid_nodes, ierr=',ierr
!
      IF (me.eq.ioproc) THEN
 
!        print *,' in fld_collect,nst, compact data'
!
!compact the data
!
       allocate ( tmp8 ( lonr,latr ) )
       do k=1,nfld
         ipt_lats=1
         do node=1,nodes
           do j=1,lats_nodes_r(node)
            lat=global_lats_r(ipt_lats-1+j)
            do i=1,lonr
              tmp8(i,lat)=grid_nodes(i,j,k,node)
            enddo
           enddo
           ipt_lats=ipt_lats+lats_nodes_r(node)
         enddo
         il=0
         do j=1,latr
            do i=1,lonsperlar(j)
              il=il+1
              bfo(il,k) = tmp8(i,j)
            enddo
         enddo
       end do
       deallocate(tmp8)
!WY bug fix.
!-----------
!       deallocate(grid_nodes)
!
        t4=rtc  ()
      endif   !me.eq.ioproc
      deallocate(grid_nodes)
!!
!
      call mpi_barrier(MPI_COMM_ALL,ierr)
!      print *,' leave fld_collect_nst ' 

      return
      end
