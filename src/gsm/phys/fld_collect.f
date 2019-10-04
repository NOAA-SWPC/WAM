      subroutine fld_collect(sfc_fld,phy_f2d,phy_f3d,ngptc,nblck,
     &                       fldsz,nfld,bfo,lats_nodes_r,global_lats_r,
     &                       lonsperlar,ioproc, nrecs)
!
!-------------------------------------------------------------------
!*** program log
!*** Dec, 2009 Jun Wang: ioproc collects partial restart file fields 
!*** Sep, 2011 Jun Wang: add cv/cvt/cvb to restart file fields 
!-------------------------------------------------------------------
!
      use resol_def,     ONLY: latr, lonr, levs, lsoil, ivssfc_restart,
     &                         ntot3d, ntot2d
      use layout1,       ONLY: lats_node_r, me,ipt_lats_node_r,
     &                         lats_node_r_max,nodes
      USE machine,       ONLY: kind_io4, kind_ior, kind_io8,kind_phys
      use mpi_def,       only: mpi_comm_all,MPI_R_MPI_R
      use gfs_physics_sfc_flx_mod, ONLY: Sfc_Var_Data
      implicit none

      TYPE(Sfc_Var_Data),intent(in)    :: sfc_fld
      INTEGER,intent(in)               :: lats_nodes_r(nodes)
      INTEGER,intent(in)               :: GLOBAL_LATS_R(latr)
      INTEGER,intent(in)               :: lonsperlar(latr)
!
      integer,intent(in) :: ngptc, nblck,fldsz,nfld,ioproc,nrecs
      REAL (KIND=KIND_phys),intent(in) ::
     &     phy_f3d(ngptc,levs,ntot3d,nblck,LATS_NODE_R)
     &,    phy_f2d(LONR,LATS_NODE_R,ntot2d)
      real(kind=kind_io8) :: bfo(fldsz,nfld)
!
!local variables:
      REAL(kind=8) t1,t2,t3,t4,t5,t6,ta,tb,rtc
!
      integer      ierr,j,k,l,lenrec,n,node,nf,ipt_lats,nrecs1
!
      real(kind=kind_ior),allocatable ::  tmp8(:,:)
!
      real(kind=kind_ior),allocatable :: grid_node(:,:,:)
      real(kind=kind_ior),allocatable :: grid_nodes(:,:,:,:)
!
      integer      lan,lat,iblk,lons_lat,lon,NJEFF,i,il
!
!---------------------------------------------------------------------
!
      if ( me  ==  ioproc ) then
!       write(0,*)' enter para_fix_w '            		! hmhj
        write(0,*)'lonr=',lonr,'lats_node_r_max=',lats_node_r_max,
     &    'total_levels=',nfld,'lonsperlar=',lonsperlar
      endif

      allocate ( grid_node ( lonr,lats_node_r_max,nfld ) )

      grid_node = 0.
      nrecs1 = nrecs+1
!
!collect data 
      do j=1,lats_node_r_max
        do i=1,lonr
          grid_node(i,j,1) = sfc_fld%tsea(i,j)
          grid_node(i,j,2) = sfc_fld%weasd(i,j)
          grid_node(i,j,3) = sfc_fld%tg3(i,j)
          grid_node(i,j,4) = sfc_fld%zorl(i,j)
          grid_node(i,j,5) = sfc_fld%cv(i,j)
          grid_node(i,j,6) = sfc_fld%cvt(i,j)
          grid_node(i,j,7) = sfc_fld%cvb(i,j)
          grid_node(i,j,8) = sfc_fld%alvsf(i,j)
          grid_node(i,j,9) = sfc_fld%alvwf(i,j)
          grid_node(i,j,10) = sfc_fld%alnsf(i,j)
          grid_node(i,j,11) = sfc_fld%alnwf(i,j)
          grid_node(i,j,12) = sfc_fld%slmsk(i,j)
          grid_node(i,j,13) = sfc_fld%vfrac(i,j)
          grid_node(i,j,14) = sfc_fld%canopy(i,j)
          grid_node(i,j,15) = sfc_fld%f10m(i,j)
!         grid_node(i,j,16) = sfc_fld%t2m(i,j)
!         grid_node(i,j,17) = sfc_fld%q2m(i,j)
          grid_node(i,j,16) = sfc_fld%vtype(i,j)
          grid_node(i,j,17) = sfc_fld%stype(i,j)
          grid_node(i,j,18) = sfc_fld%facsf(i,j)
          grid_node(i,j,19) = sfc_fld%facwf(i,j)
          grid_node(i,j,20) = sfc_fld%uustar(i,j)
          grid_node(i,j,21) = sfc_fld%ffmm(i,j)
          grid_node(i,j,22) = sfc_fld%ffhh(i,j)
          grid_node(i,j,23) = sfc_fld%hice(i,j)
          grid_node(i,j,24) = sfc_fld%fice(i,j)
          grid_node(i,j,25) = sfc_fld%tisfc(i,j)
          grid_node(i,j,26) = sfc_fld%tprcp(i,j)
          grid_node(i,j,27) = sfc_fld%srflag(i,j)
          grid_node(i,j,28) = sfc_fld%snwdph(i,j)
          grid_node(i,j,29) = sfc_fld%shdmin(i,j)
          grid_node(i,j,30) = sfc_fld%shdmax(i,j)
          grid_node(i,j,31) = sfc_fld%slope(i,j)
          grid_node(i,j,32) = sfc_fld%snoalb(i,j)
          grid_node(i,j,33) = sfc_fld%oro(i,j)
          grid_node(i,j,34) = sfc_fld%sncovr(i,j)
        enddo
      enddo
      do l=1,lsoil
        do j=1,lats_node_r_max
          do i=1,lonr
            grid_node(i,j,nrecs+l)         = sfc_fld%smc(l,i,j)
            grid_node(i,j,nrecs+l+lsoil)   = sfc_fld%stc(l,i,j)
            grid_node(i,j,nrecs+l+2*lsoil) = sfc_fld%slc(l,i,j)
          enddo
        enddo
      enddo
!
!phy_f2d
      do k=1,ntot2d
        do j=1,lats_node_r
          do i=1,lonr
            grid_node(i,j,nrecs+3*lsoil+k) = phy_f2d(i,j,k)
          enddo
        enddo
      enddo
!
!phy_f3d
      do l=1,ntot3d
        do k=1,levs
          do j=1,lats_node_r
            lat = global_lats_r(ipt_lats_node_r-1+j)
            lons_lat = lonsperlar(lat)
            iblk = 0
            il   = 1
            do lon=1,lons_lat,NGPTC
              NJEFF = MIN(NGPTC,lons_lat-lon+1)
              iblk  = iblk+1
              do i=1,NJEFF
                grid_node(il,j,nrecs+3*lsoil+ntot2d+(l-1)*levs+k) = 
     &            phy_f3d(i,k,l,iblk,j)
                il = il+1
              enddo
            enddo
          enddo
        enddo
      enddo
!     write(0,*) 'in fld_collect,nfld=',nfld,'total field=',
!    &   nrecs+3*lsoil+ntot2d+ntot3d*levs,'nodes=',nodes
!
!WY bug fix.
!-----------
      if ( me  == ioproc ) then
         allocate ( grid_nodes ( lonr,lats_node_r_max,
     &                           nfld, nodes ),stat=ierr )
      else
         allocate(grid_nodes(1, 1, 1, 1), stat = ierr)
      endif
      if (ierr /= 0) then
        call mpi_abort(mpi_comm_all,ierr,i)
      endif
!
      lenrec = lonr*lats_node_r_max * nfld
!
      t1 = rtc()
!      print *,'after allocate grid_nodes,lenrec=',lenrec
      if (nodes > 1) then
        call mpi_gather( grid_node , lenrec, MPI_R_MPI_R,
     x                   grid_nodes, lenrec, MPI_R_MPI_R,
     x                   ioproc, MPI_COMM_ALL, ierr)
      elseif(nodes == 1) then
        grid_nodes(:,:,:,1) = grid_node(:,:,:)
      else
        print *,'ERROR: the totals nodes=',nodes
        call mpi_abort(mpi_comm_all,ierr,i)
      endif
!
      deallocate(grid_node)
!
      IF (me == ioproc) THEN
 
!
!compact the data
!
        allocate ( tmp8 (lonr,latr) )
        do k=1,nfld
          ipt_lats=1
          do node=1,nodes
            do j=1,lats_nodes_r(node)
              lat = global_lats_r(ipt_lats-1+j)
              do i=1,lonr
                tmp8(i,lat) = grid_nodes(i,j,k,node)
              enddo
            enddo
            ipt_lats = ipt_lats+lats_nodes_r(node)
          enddo
          il = 0
          do j=1,latr
             do i=1,lonsperlar(j)
               il = il+1
               bfo(il,k) = tmp8(i,j)
             enddo
          enddo
        end do
        deallocate(tmp8)
!WY bug fix.
!-----------
!       deallocate(grid_nodes)
!
        t4 = rtc  ()
      endif                                    ! me == ioproc

      deallocate(grid_nodes)
!!
!
      call mpi_barrier(MPI_COMM_ALL,ierr)
!      print *,' leave fld_collect ' 

      return
      end
