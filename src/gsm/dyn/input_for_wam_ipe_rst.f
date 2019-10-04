      SUBROUTINE input_for_wam_ipe_rst(global_lats_a,lonsperlat, 
     &                                 lats_nodes_a)
!!
!  Nov 14 2017    Weiyu Yang - create this subroutine to read in the
!                              WAM-IPE coupling interface fields, which will 
!                              send to the IPE model for the first step 
!                              WAM-IPE restart run.
!! Revision history:
!

      use get_variables_for_WAM_IPE_coupling, ONLY: wwg, zzg, uug,
     &                                         vvg, ttg, rqg, n2g
      use gfs_dyn_resol_def
      use gfs_dyn_layout1
      use gfs_dyn_mpi_def

      implicit none

      integer, dimension(latg) :: global_lats_a, lonsperlat
      integer, dimension(nodes_comp) :: lats_nodes_a
!
      real(kind=kind_grid), dimension(:, :, :), allocatable :: 
     &                                           buff_mult_pieceg_ipe
      real(kind=kind_grid),dimension(:,:,:),  allocatable :: buff_final
      real(kind=kind_grid),dimension(:,:,:),  allocatable :: grid_node
      real(kind=kind_grid),dimension(:,:,:,:),allocatable :: grid_nodes
!
      integer i, j, k, ngrids_gg_ipe, ioproc, lenrec
      integer lat, ipt_lats, ierr, node
!
      ngrids_gg_ipe = 6*levs+levh

      ioproc = nodes_comp - 1

      if(.not. allocated(grid_node)) then
         allocate(grid_node(lonf,lats_node_a_max,ngrids_gg_ipe))
      end if

      if(me == ioproc) then
        if(.not. allocated(buff_final)) then
           allocate(buff_final(lonf, latg, ngrids_gg_ipe))
        endif
        if(.not. allocated(grid_nodes)) then
           allocate(grid_nodes(lonf, lats_node_a_max, ngrids_gg_ipe, 
     &                         nodes_comp))
           grid_nodes = 0.0
        endif

! The following is to read in the interface restart file for WAM-IPE
! coupling restart run.
!-------------------------------------------------------------------
        rewind 180
        READ(180) buff_final
        
        DO k = 1, ngrids_gg_ipe
          ipt_lats = 1
          DO node = 1, nodes_comp
            DO j = 1, lats_nodes_a(node)
              lat = global_lats_a(ipt_lats-1+j)
              DO i = 1, lonf
                grid_nodes(i, j, k, node) = buff_final(i, lat, k)
              END DO
            END DO
            ipt_lats = ipt_lats+lats_nodes_a(node)
          END DO
        END DO
      else
        if(.not. allocated(grid_nodes)) then
           allocate(grid_nodes(1, 1, 1, 1))
        endif
      endif

      IF(.NOT. ALLOCATED(wwg)) THEN
         ALLOCATE(wwg(lonf,lats_node_a,levs))
         wwg = 0.0
      endif
      IF(.NOT. ALLOCATED(zzg)) then
         ALLOCATE(zzg(lonf,lats_node_a,levs))
         zzg = 0.0
      endif
      IF(.NOT. ALLOCATED(uug)) then
         ALLOCATE(uug(lonf,lats_node_a,levs))
         uug = 0.0
      endif
      IF(.NOT. ALLOCATED(vvg)) then
         ALLOCATE(vvg(lonf,lats_node_a,levs))
         vvg=0.0
      endif
      IF(.NOT. ALLOCATED(ttg)) then
         ALLOCATE(ttg(lonf,lats_node_a,levs))
         ttg = 0.0
      endif
      IF(.NOT. ALLOCATED(rqg)) then
         ALLOCATE(rqg(lonf,lats_node_a,levh))
         rqg = 0.0
      endif
      IF(.NOT. ALLOCATED(n2g)) then
         ALLOCATE(n2g(lonf,lats_node_a,levs))
         n2g = 0.0
      endif

      if(nodes_comp>1) then
        lenrec = lonf * lats_node_a_max * ngrids_gg_ipe
        call mpi_scatter(grid_nodes, lenrec, mpi_real8,
     x                 grid_node, lenrec, mpi_real8,
     x                 ioproc, MPI_COMM_ALL, ierr)
      else
        grid_node(:,:,:)=grid_nodes(:,:,:,1)
      end if

      if(.not. allocated(buff_mult_pieceg_ipe)) then
         allocate(buff_mult_pieceg_ipe(lonf,lats_node_a,ngrids_gg_ipe))
      endif
!
      CALL mpi_barrier(mpi_comm_all, ierr)

      DO k = 1, ngrids_gg_ipe
        DO j = 1, lats_node_a
          DO i = 1, lonf
            buff_mult_pieceg_ipe(i, j, k) = grid_node(i, j, k)
          END DO
        END DO
      END DO

      do k=1,levs
        do j=1,lats_node_a
          do i=1,lonf
            wwg(i,j,k) = buff_mult_pieceg_ipe(i, j, k)
          enddo
        enddo
      enddo
!!
      do k=1,levs
        do j=1,lats_node_a
          do i=1,lonf
            zzg(i,j,k) = buff_mult_pieceg_ipe(i, j, levs + k)
          enddo
        enddo
      enddo
!
      do k=1,levs
        do j=1,lats_node_a
          do i=1,lonf
            uug(i,j,k) = buff_mult_pieceg_ipe(i, j, 2 * levs + k)
          enddo
        enddo
      enddo
!
      do k=1,levs
        do j=1,lats_node_a
          do i=1,lonf
            vvg(i,j,k) = buff_mult_pieceg_ipe(i, j, 3 * levs + k)
          enddo
        enddo
      enddo
!
      do k=1,levs
        do j=1,lats_node_a
          do i=1,lonf
            ttg(i,j,k) = buff_mult_pieceg_ipe(i, j, 4 * levs + k)
          enddo
        enddo
      enddo
!
      if (levh > 0) then
        do k=1,levh
          do j=1,lats_node_a
            do i=1,lonf
              rqg(i,j,k) = buff_mult_pieceg_ipe(i, j, 5 * levs + k)
            enddo
          enddo
        enddo
      endif
!
      do k=1,levs
        do j=1,lats_node_a
          do i=1,lonf
            n2g(i,j,k) = buff_mult_pieceg_ipe(i, j, 5 * levs + levh + k)
          enddo
        enddo
      enddo

! For test.
!----------
!      print*, 'read wam-ipe cpl rst file, me, uug=', me, uug(3,1,100)
!      print*, 'read wam-ipe cpl rst file, me, vvg=', me, vvg(3,1,100)
!      print*, 'read wam-ipe cpl rst file, me, ttg=', me, ttg(3,1,100)
!      print*, 'read wam-ipe cpl rst file, me, wwg=', me, wwg(3,1,100)
!      print*, 'read wam-ipe cpl rst file, me, zzg=', me, zzg(3,1,100)
!      print*, 'read wam-ipe cpl rst file, me, n2g=', me, n2g(3,1,100)
!      print*, 'read wam-ipe cpl rst file, me, rqg=', me, rqg(3,1,100),
!     &          rqg(3,1,250),rqg(3,1,400)

      END SUBROUTINE input_for_wam_ipe_rst
