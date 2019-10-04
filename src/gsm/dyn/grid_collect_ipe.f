      SUBROUTINE grid_collect_ipe(wwg,zzg,uug,vvg,ttg,rqg,n2g,
     &              global_lats_a,lonsperlat, lats_nodes_a, kdt, 
     &              deltim)
!!
!! Revision history:
!  2007           Henry Juang, original code
!  2008           Jun Wang  modified buff for write grid component
!  Nov 23 2009    Sarah Lu, comment out 4D tracer
!  Sep 08 2010    Jun Wang  change gfsio to nemsio
!  Dec 16 2010    Jun Wang  change to nemsio library
!  Feb 20 2011    Hann-Ming Henry Juang add code for NDSL
!  Sep 24 2014    S Moorthi - some cleanup and optimization
!  Feb 04 2015    S. Moorthi - threading and optimization
!  May 17 2016    Weiyu Yang - modified from grid_collect.f for WAM-IPE
!                              coupling outputs.
!  Nov 08 2017    Weiyu Yang - modified for adding the NetCDF diagnostic 
!                              outputs and the WAM-IPE coupling restart
!                              functions.
!

      use gfs_dyn_resol_def
      use gfs_dyn_layout1
      use gfs_dyn_mpi_def
      use namelist_dynamics_def, ONLY: wam_ipe_cpl_rst_output,
     &                                 grads_output, NC_output

      implicit none

      integer, dimension(latg) :: global_lats_a, lonsperlat
      integer, dimension(nodes_comp) :: lats_nodes_a
!
      real(kind=kind_grid), dimension(lonf,lats_node_a,levs) ::  uug,vvg
     &,                                                          ttg,wwg
     &,                                                          zzg,n2g
      real(kind=kind_grid), dimension(lonf,lats_node_a,levh) ::  rqg
      real(kind=kind_grid), dimension(:, :, :), allocatable :: 
     &                                              buff_mult_pieceg_ipe
!
      real(kind=kind_grid), dimension(lonf,lats_node_a) :: buffi
      integer, dimension(lonf,lats_node_a)             :: kmsk
      integer i, j, k, ngrids_gg_ipe, kdt
      real   deltim
!
! For test.
!----------
!      print*, 'me, uug=', me, uug(3,1,100)
!      print*, 'me, vvg=', me, vvg(3,1,100)
!      print*, 'me, ttg=', me, ttg(3,1,100)
!      print*, 'me, wwg=', me, wwg(3,1,100)
!      print*, 'me, zzg=', me, zzg(3,1,100)
!      print*, 'me, n2g=', me, n2g(3,1,100)
!      print*, 'me, rqg=', me, rqg(3,1,100),rqg(3,1,250),rqg(3,1,400)

      ngrids_gg_ipe = 6*levs+levh
      kmsk = 0

      if(.not. allocated(buff_mult_pieceg_ipe)) then
         allocate(buff_mult_pieceg_ipe(lonf,lats_node_a,ngrids_gg_ipe))
         buff_mult_pieceg_ipe = 0.0
      endif
!
      do k=1,levs
!$omp parallel do private(i,j)
        do j=1,lats_node_a
          do i=1,lonf
            buffi(i,j) = wwg(i,j,k)
          enddo
        enddo
        IF(grads_output .OR. NC_output) THEN
          CALL uninterpred_dyn(1,kmsk,buff_mult_pieceg_ipe(1,1,k),buffi,
     &                     global_lats_a,lonsperlat)
        ELSE IF(wam_ipe_cpl_rst_output) THEN
          do j=1,lats_node_a
            do i=1,lonf
              buff_mult_pieceg_ipe(i, j, k) = buffi(i,j)
            enddo
          enddo
        END IF

!      write(0,*)'in grid collect, buff_wwg=',' me=',me,
!    & maxval(buff_mult_pieceg_ipe(1:lonf,1:lats_node_a,k)),
!    & minval(buff_mult_pieceg_ipe(1:lonf,1:lats_node_a,k))
      enddo
!!
      do k=1,levs
!$omp parallel do private(i,j)
        do j=1,lats_node_a
          do i=1,lonf
            buffi(i,j) = zzg(i,j,k)
          enddo
        enddo
        IF(grads_output .OR. NC_output) THEN
          CALL uninterpred_dyn(1,kmsk,
     &       buff_mult_pieceg_ipe(1,1,levs+k),buffi,
     &       global_lats_a,lonsperlat)
        ELSE IF(wam_ipe_cpl_rst_output) THEN
          do j=1,lats_node_a
            do i=1,lonf
              buff_mult_pieceg_ipe(i, j, levs+k) = buffi(i,j) 
            enddo
          enddo
        END IF

!      write(0,*)'in grid collect, buff_zzg=',' me=',me,
!    & maxval(buff_mult_pieceg_ipe(1:lonf,1:lats_node_a,levs+k)),
!    & minval(buff_mult_pieceg_ipe(1:lonf,1:lats_node_a,levs+k))
      enddo
!
      do k=1,levs
!$omp parallel do private(i,j)
        do j=1,lats_node_a
          do i=1,lonf
            buffi(i,j) = uug(i,j,k)
          enddo
        enddo
        IF(grads_output .OR. NC_output) THEN
          CALL uninterpred_dyn(1,kmsk,
     &       buff_mult_pieceg_ipe(1,1,2*levs+k),buffi,
     &       global_lats_a,lonsperlat)
        ELSE IF(wam_ipe_cpl_rst_output) THEN
          do j=1,lats_node_a
            do i=1,lonf
              buff_mult_pieceg_ipe(i, j, 2*levs+k) = buffi(i,j)
            enddo
          enddo
        END IF

!      write(0,*)'in grid collect, buff_uug=',' me=',me,
!    & maxval(buff_mult_pieceg_ipe(1:lonf,1:lats_node_a,2*levs+k)),
!    & minval(buff_mult_pieceg_ipe(1:lonf,1:lats_node_a,2*levs+k))
      enddo
!
      do k=1,levs
!$omp parallel do private(i,j)
        do j=1,lats_node_a
          do i=1,lonf
            buffi(i,j) = vvg(i,j,k)
          enddo
        enddo
        IF(grads_output .OR. NC_output) THEN
          CALL uninterpred_dyn(1,kmsk,
     &       buff_mult_pieceg_ipe(1,1,3*levs+k),buffi,
     &       global_lats_a,lonsperlat)
        ELSE IF(wam_ipe_cpl_rst_output) THEN
          do j=1,lats_node_a
            do i=1,lonf
              buff_mult_pieceg_ipe(i, j, 3*levs+k) = buffi(i,j)
            enddo
          enddo
        END IF

!      write(0,*)'in grid collect, buff_vvg=',' me=',me,
!    & maxval(buff_mult_pieceg_ipe(1:lonf,1:lats_node_a,3*levs+k)),
!    & minval(buff_mult_pieceg_ipe(1:lonf,1:lats_node_a,3*levs+k))
      enddo
!
      do k=1,levs
!$omp parallel do private(i,j)
        do j=1,lats_node_a
          do i=1,lonf
            buffi(i,j) = ttg(i,j,k)
          enddo
        enddo
        IF(grads_output .OR. NC_output) THEN
          CALL uninterpred_dyn(1,kmsk,
     &       buff_mult_pieceg_ipe(1,1,4*levs+k),buffi,
     &       global_lats_a,lonsperlat)
        ELSE IF(wam_ipe_cpl_rst_output) THEN
          do j=1,lats_node_a
            do i=1,lonf
              buff_mult_pieceg_ipe(i, j, 4*levs+k) = buffi(i,j)
            enddo
          enddo
        END IF

!      write(0,*)'in grid collect, buff_ttg=',' me=',me,
!    & maxval(buff_mult_pieceg_ipe(1:lonf,1:lats_node_a,4*levs+k)),
!    & minval(buff_mult_pieceg_ipe(1:lonf,1:lats_node_a,4*levs+k))
      enddo
!
      if (levh > 0) then
        do k=1,levh
!$omp parallel do private(i,j)
          do j=1,lats_node_a
            do i=1,lonf
              buffi(i,j) = rqg(i,j,k)
            enddo
          enddo
          IF(grads_output .OR. NC_output) THEN
          CALL uninterpred_dyn(1,kmsk,
     &       buff_mult_pieceg_ipe(1,1,5*levs+k),buffi,
     &       global_lats_a,lonsperlat)
          ELSE IF(wam_ipe_cpl_rst_output) THEN
            do j=1,lats_node_a
              do i=1,lonf
                buff_mult_pieceg_ipe(i, j, 5*levs+k) = buffi(i,j)
              enddo
            enddo
          END IF

!      write(0,*)'in grid collect, buff_rqg=',' me=',me,
!    & maxval(buff_mult_pieceg_ipe(1:lonf,1:lats_node_a,5*levs+k)),
!    & minval(buff_mult_pieceg_ipe(1:lonf,1:lats_node_a,5*levs+k))
        enddo
      endif
!
      do k=1,levs
!$omp parallel do private(i,j)
        do j=1,lats_node_a
          do i=1,lonf
            buffi(i,j) = n2g(i,j,k)
          enddo
        enddo
        IF(grads_output .OR. NC_output) THEN
          CALL uninterpred_dyn(1,kmsk,
     &       buff_mult_pieceg_ipe(1,1,5*levs+levh+k),buffi,
     &       global_lats_a,lonsperlat)
        ELSE IF(wam_ipe_cpl_rst_output) THEN
          do j=1,lats_node_a
            do i=1,lonf
              buff_mult_pieceg_ipe(i, j, 5*levs+levh+k) = buffi(i,j)
            enddo
          enddo
        END IF

!      write(0,*)'in grid collect, buff_n2g=',' me=',me,
!    & maxval(buff_mult_pieceg_ipe(1:lonf,1:lats_node_a,5*levs+levh+k)),
!    & minval(buff_mult_pieceg_ipe(1:lonf,1:lats_node_a,5*levs+levh+k))
      enddo

      CALL atmgg_move_ipe(buff_mult_pieceg_ipe,ngrids_gg_ipe, kdt,
     &                    global_lats_a, lats_nodes_a, deltim)

      return
      end

      subroutine atmgg_move_ipe(buff_mult_pieceg_ipe,ngrids_gg_ipe, kdt,
     &                          global_lats_a, lats_nodes_a, deltim)
c
c***********************************************************************
c
      use gfs_dyn_resol_def
      use gfs_dyn_write_state
      use gfs_dyn_layout1
      use gfs_dyn_mpi_def
      use namelist_dynamics_def, ONLY: wam_ipe_cpl_rst_output,
     &                                 grads_output,
     &                                 FHOUT_grads, NC_output,
     &                                 FHOUT_NC, FHRES
      implicit none
!
      integer ngrids_gg_ipe
!      real(kind=kind_io4), dimension(lonf,lats_node_a,ngrids_gg_ipe)
!     &                                           :: buff_mult_pieceg_ipe
!      real(kind=kind_io4), dimension(lonf,lats_node_a_max,ngrids_gg_ipe)
!     &                                           :: grid_node
!      real(kind=kind_io4),dimension(:,:,:),  allocatable :: buff_final
!      real(kind=kind_io4),dimension(:,:,:,:),allocatable :: grid_nodes
      real(kind=kind_grid), dimension(lonf,lats_node_a,ngrids_gg_ipe)
     &                                           :: buff_mult_pieceg_ipe
      real(kind=kind_grid),dimension(lonf,lats_node_a_max,ngrids_gg_ipe)
     &                                           :: grid_node
      real(kind=kind_grid),dimension(:,:,:),  allocatable :: buff_final
      real(kind=kind_grid),dimension(:,:,:,:),allocatable :: grid_nodes
      integer, dimension(latg)       :: global_lats_a
      integer, dimension(nodes_comp) :: lats_nodes_a
      integer ioproc, kdt, lat, ipt_lats
      integer j,k,i,ierr, node
      integer lenrec
      real    deltim
!
      ioproc = nodes_comp - 1

      grid_node = 0.0
      DO k = 1, ngrids_gg_ipe
        DO j = 1, lats_node_a
          DO i = 1, lonf
            grid_node(i, j, k) = buff_mult_pieceg_ipe(i, j, k)
          END DO
        END DO
      END DO
!!
      if(me == ioproc) then
        if(.not. allocated(buff_final)) then
           allocate(buff_final(lonf, latg, ngrids_gg_ipe))
        endif
        if(.not. allocated(grid_nodes)) then
           allocate(grid_nodes(lonf, lats_node_a_max, ngrids_gg_ipe, 
     &                         nodes_comp))
        endif
      else
        if(.not. allocated(grid_nodes)) then
           allocate(grid_nodes(1, 1, 1, 1))
        endif
      endif
!
      if(nodes_comp>1) then
        lenrec = lonf * lats_node_a_max * ngrids_gg_ipe
!
        call mpi_gather( grid_node , lenrec, mpi_real8,
     x                 grid_nodes, lenrec, mpi_real8,
     x                 ioproc, MPI_COMM_ALL, ierr)
      else
        grid_nodes(:,:,:,1)=grid_node(:,:,:)
      endif

      IF(me == ioproc) THEN
        DO k = 1, ngrids_gg_ipe
          ipt_lats = 1
          DO node = 1, nodes_comp
            DO j = 1, lats_nodes_a(node)
              lat = global_lats_a(ipt_lats-1+j)
              DO i = 1, lonf
                buff_final(i, lat, k) = grid_nodes(i, j, k, node)
              END DO
            END DO
            ipt_lats = ipt_lats+lats_nodes_a(node)
          END DO
        END DO
      END IF

      call mpi_barrier(mpi_comm_all,ierr)
      deallocate(grid_nodes)

! Write out the wwg, zzg, uug, vvg, ttg, rqg, n2g full grid fields to
! disk.
!--------------------------------------------------------------------
! buff_final contains wwg, zzg, uug, vvg, ttg, rqg, n2g.
!-------------------------------------------------------
      if(me == ioproc) then

! The following is only to output for making figures and comparisons.
! WY.
!------------------------------------------------------------------------
        IF(grads_output .AND.
     &    MOD(NINT(deltim) * kdt, FHOUT_grads * 3600) == 0) THEN
          PRINT*,'Output the WAM-IPE coupling fields at kdt=', kdt
          print*, 'kdt, lonf, latg, ngrids_gg_ipe=',kdt, lonf, latg,
     &             ngrids_gg_ipe
          write(178) kdt, lonf, latg, ngrids_gg_ipe
          write(178) buff_final
        END IF

! The following is to the NetCDF diagnostic files.
!-------------------------------------------------
!        IF(NC_output .AND.
!     &    MOD(NINT(deltim) * kdt, FHOUT_NC * 3600) == 0) THEN
!          CALL write_NC(kdt, lonf, latg, ngrids_gg_ipe,
!     &                  buff_final)
!        END IF

! The following is to output the interface restart file for WAM-IPE
! coupling restart run.
!------------------------------------------------------------------
! Joe Schoonover ( Jan. 31, 2018 ) : Needed to change kdt to
! kdt+1 to write the neutrals at a multiple of FHRES, instead
! of a FHRES+dt
        IF(wam_ipe_cpl_rst_output .AND.
     &    MOD(NINT(deltim) * (kdt+1), NINT(FHRES) * 3600) == 0) THEN
! restart file will keep the last output file.
!---------------------------------------------
          rewind 181
          WRITE(181) buff_final
        END IF

        deallocate(buff_final)
      end if   ! if(me == ioproc). 

      return
      end

      subroutine uninterpred_dyn(iord,kmsk,f,fi,global_lats_a,
     &     lonsperlat)
!!
      use gfs_dyn_resol_def,   ONLY: latg, lonf
      use gfs_dyn_layout1,     ONLY: lats_node_a, ipt_lats_node_a
      USE machine,     ONLY: kind_io8
      implicit none
!!
      integer              global_lats_a(latg)
      integer,intent(in):: iord
      integer,intent(in):: kmsk(lonf,lats_node_a)
      integer,intent(in):: lonsperlat(latg)
      real(kind=kind_io8),intent(out):: f(lonf,lats_node_a)
      real(kind=kind_io8),intent(in) :: fi(lonf,lats_node_a)
      integer j,lons,lat
!!
      do j=1,lats_node_a
        lat  = global_lats_a(ipt_lats_node_a-1+j)
        lons = lonsperlat(lat)
        if(lons .ne. lonf) then
          call intlon_dyn(iord,1,1,lons,lonf,
     &                     kmsk(1,j),fi(1,j),f(1,j))
        else
          f(:,j) = fi(:,j)
        endif
      enddo
      end subroutine

      subroutine intlon_dyn(iord,imon,imsk,m1,m2,k1,f1,f2)
      use machine, ONLY: kind_io8
      implicit none
      integer,intent(in):: iord,imon,imsk,m1,m2
      integer,intent(in):: k1(m1)
      real (kind=kind_io8),intent(in):: f1(m1)
      real (kind=kind_io8),intent(out):: f2(m2)
      integer i2,in,il,ir
      real (kind=kind_io8) r,x1
      r=real(m1)/real(m2)
      do i2=1,m2
         x1=(i2-1)*r
         il=int(x1)+1
         ir=mod(il,m1)+1
          if(iord.eq.2.and.(imsk.eq.0.or.k1(il).eq.k1(ir))) then
            f2(i2)=f1(il)*(il-x1)+f1(ir)*(x1-il+1)
          else
            in=mod(nint(x1),m1)+1
            f2(i2)=f1(in)
          endif
      enddo
      end subroutine
