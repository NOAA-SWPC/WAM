      module gci
      implicit none

      contains
      SUBROUTINE grid_collect_ipe(wwg,zzg,uug,vvg,ttg,rqg,n2g,
     &             global_lats_a,lonsperlat, lats_nodes_a, kdt, deltim,
     &             restart_step, den, gmol)
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
!  Dec 09 2019    A. Kubaryk - near-total rewrite

      use gfs_dyn_resol_def
      use gfs_dyn_layout1
      use gfs_dyn_mpi_def

      implicit none

      integer, dimension(latg) :: global_lats_a, lonsperlat
      integer, dimension(nodes_comp) :: lats_nodes_a
!
      real(kind=kind_grid), dimension(lonf,lats_node_a,levs) ::  uug,vvg
     &,                                                          ttg,wwg
     &,                                                          zzg,n2g
      real(kind=kind_grid), dimension(lonf,lats_node_a,levh) ::  rqg
      logical, intent(in) :: restart_step
      real(kind=kind_grid), dimension(lonf,lats_node_a,levs),
     &                      optional                         :: den,gmol
!
      real(kind=kind_grid), dimension(:, :, :), allocatable ::
     &                                              buff_mult_pieceg_ipe
!
      integer i, j, k, ngrids_gg_ipe, kdt
      integer, dimension(lonf,lats_node_a)             :: kmsk
      logical nc_out
      real deltim
!
      if (PRESENT(den) .and. PRESENT(gmol)) then
        nc_out = .true.
      else
        nc_out = .false.
      end if

      if (nc_out) then
        ngrids_gg_ipe = 8*levs+levh
      else
        ngrids_gg_ipe = 6*levs+levh
      end if

      if(.not. allocated(buff_mult_pieceg_ipe)) then
         allocate(buff_mult_pieceg_ipe(lonf,lats_node_a,ngrids_gg_ipe))
      endif

      if (nc_out) then
        do k=1,levs
          call uninterpred_dyn(1,kmsk,
     &                         buff_mult_pieceg_ipe(1,1,k),
     &                         wwg(:,:,k),global_lats_a,lonsperlat)
          call uninterpred_dyn(1,kmsk,
     &                         buff_mult_pieceg_ipe(1,1,k+levs),
     &                         zzg(:,:,k),global_lats_a,lonsperlat)
          call uninterpred_dyn(1,kmsk,
     &                         buff_mult_pieceg_ipe(1,1,k+levs*2),
     &                         uug(:,:,k),global_lats_a,lonsperlat)
          call uninterpred_dyn(1,kmsk,
     &                         buff_mult_pieceg_ipe(1,1,k+levs*3),
     &                         vvg(:,:,k),global_lats_a,lonsperlat)
          call uninterpred_dyn(1,kmsk,
     &                         buff_mult_pieceg_ipe(1,1,k+levs*4),
     &                         ttg(:,:,k),global_lats_a,lonsperlat)
          call uninterpred_dyn(1,kmsk,
     &                         buff_mult_pieceg_ipe(1,1,k+levs*5+levh),
     &                         n2g(:,:,k),global_lats_a,lonsperlat)
          call uninterpred_dyn(1,kmsk,
     &                         buff_mult_pieceg_ipe(1,1,k+levs*6+levh),
     &                         den(:,:,k),global_lats_a,lonsperlat)
          call uninterpred_dyn(1,kmsk,
     &                         buff_mult_pieceg_ipe(1,1,k+levs*7+levh),
     &                         gmol(:,:,k),global_lats_a,lonsperlat)
        end do
        do k=1,levh
          call uninterpred_dyn(1,kmsk,
     &                         buff_mult_pieceg_ipe(1,1,k+levs*5),
     &                         rqg(:,:,k),global_lats_a,lonsperlat)
        end do
      else
            buff_mult_pieceg_ipe(:,:,       1:  levs)      = wwg
            buff_mult_pieceg_ipe(:,:,  levs+1:2*levs)      = zzg
            buff_mult_pieceg_ipe(:,:,2*levs+1:3*levs)      = uug
            buff_mult_pieceg_ipe(:,:,3*levs+1:4*levs)      = vvg
            buff_mult_pieceg_ipe(:,:,4*levs+1:5*levs)      = ttg
            buff_mult_pieceg_ipe(:,:,5*levs+1:5*levs+levh) = rqg
            buff_mult_pieceg_ipe(:,:,5*levs+levh+1:6*levs+levh) = n2g
      end if

      CALL atmgg_move_ipe(buff_mult_pieceg_ipe,ngrids_gg_ipe, kdt,
     &                    global_lats_a, lats_nodes_a, deltim,
     &                    nc_out, restart_step)

      END SUBROUTINE grid_collect_ipe

      subroutine atmgg_move_ipe(buff_mult_pieceg_ipe,ngrids_gg_ipe, kdt,
     &                       global_lats_a, lats_nodes_a, deltim,
     &                       nc_out, restart_step)
!!!
      use gfs_dyn_resol_def
      use gfs_dyn_write_state
      use gfs_dyn_layout1
      use gfs_dyn_mpi_def
      use namelist_dynamics_def, ONLY: wam_ipe_cpl_rst_output,
     &                                 NC_output, FHRES, ens_nam
      implicit none
!
      integer ngrids_gg_ipe
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
      integer lenrec, ndig, nfill
      real    deltim
      logical nc_out, restart_step
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

! Write out the wwg, zzg, uug, vvg, ttg, rqg, n2g,(den, gmol), full 
! grid fields to disk.
!------------------------------------------------------------------
! buff_final contains wwg, zzg, uug, vvg, ttg, rqg, n2g (den, gmol).
!-------------------------------------------------------------------
      if(me == ioproc) then
        if (nc_out) then
          call write_nc(kdt.eq.0 .or. restart_step,
     &                   lonf, latg, ngrids_gg_ipe,buff_final)
        else
          PRINT*, 'write out WAM IPE CPL RST file, kdt=', kdt
          rewind 181
          WRITE(181) buff_final
        END IF

        deallocate(buff_final)
      end if   ! if(me == ioproc). 

      end subroutine atmgg_move_ipe


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
      end subroutine uninterpred_dyn

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

      subroutine write_nc(first_step,lonf,latg,ngrids_gg_ipe,buff_final)
      use netcdf
      use gfs_dyn_gg_def, only:    sinlat_a
      use gfs_dyn_resol_def, only: levs
      use namelist_dynamics_def, ONLY: nc_fields
      implicit none

      logical, intent(in) :: first_step
      integer, intent(in) :: lonf
      integer, intent(in) :: latg
      integer, intent(in) :: ngrids_gg_ipe
      real, dimension(lonf,latg,ngrids_gg_ipe), intent(in) :: buff_final

      ! Local Variables
      integer, parameter :: fields = 13
      character(2), dimension(fields),parameter::var=(/'w ', 'z ', 'u ',&
     &                                           'v ', 't ', 'qr', 'o3',&
     &                                           'cw', 'o ', 'o2', 'n2',&
     &                                           'dn', 'gm'/)
      character(7),dimension(fields),parameter::units=(/'m/s','m','m/s',&
     &                                        'm/s','K','kg/kg','kg/kg',&
     &                                     'kg/kg','m^-3','m^-3','m^-3',&
     &                                         'kg*m^-3','kg/mol'/)
      real, dimension(latg) :: lats
      real, dimension(lonf) :: lons
      integer, dimension(levs) :: levels
      integer :: ncstatus, ncid, x_dimid, y_dimid, z_dimid, time_dimid
      integer :: xt_dimid, yt_dimid, zt_dimid
      integer :: varid, i, k, time
      integer :: start(4), count(4)
      character(13), parameter :: filename='wam_fields.nc'
      real :: pi
      logical :: file_exists

      pi = atan(1.0)*4.0
      count = (/lonf,latg,levs,1/)

      lats = asin(sinlat_a)*180. / pi
      do i = 1,lonf
        lons(i) = (i-1) * 360. / lonf
      enddo
      do i = 1,levs
        levels(i) = i
      enddo

      INQUIRE(FILE=filename, EXIST=file_exists)

      if (.not. file_exists .or. first_step) then
!        create
         ncstatus=nf90_create(filename, NF90_NETCDF4, ncid )

!        dimensions
         ncstatus=nf90_def_dim(ncid, "lon",  lonf, x_dimid)
         ncstatus=nf90_def_dim(ncid, "lat",  latg, y_dimid)
         ncstatus=nf90_def_dim(ncid, "levs", levs, z_dimid)
         ncstatus=nf90_def_dim(ncid, "time", nf90_unlimited,time_dimid)
         ncstatus=nf90_put_att(ncid, time_dimid, "axis", "T")

         ncstatus=nf90_def_var(ncid, "lon", NF90_FLOAT, (/x_dimid/),    &
     &                         xt_dimid)
         ncstatus=nf90_put_att(ncid, xt_dimid, "axis", "X")
         ncstatus=nf90_put_att(ncid, xt_dimid, "long_name", "longitude")
         ncstatus=nf90_put_att(ncid, xt_dimid, "units", "degrees_east")
         ncstatus=nf90_put_var(ncid, xt_dimid, lons)

         ncstatus=nf90_def_var(ncid, "lat", NF90_FLOAT, (/y_dimid/),    &
     &                         yt_dimid)
         ncstatus=nf90_put_att(ncid, yt_dimid, "axis", "Y")
         ncstatus=nf90_put_att(ncid, yt_dimid, "long_name", "latitude")
         ncstatus=nf90_put_att(ncid, yt_dimid, "units", "degrees_north")
         ncstatus=nf90_put_var(ncid, yt_dimid, lats)

         ncstatus=nf90_def_var(ncid, "levs", NF90_FLOAT, (/z_dimid/),   &
     &                         zt_dimid)
         ncstatus=nf90_put_att(ncid, zt_dimid, "axis", "Z")
         ncstatus=nf90_put_att(ncid, zt_dimid, "long_name", "level")
         ncstatus=nf90_put_var(ncid, zt_dimid, levels)

         do i=1,fields
            if (nc_fields(i)) then
            ncstatus=nf90_def_var(ncid,trim(var(i)), NF90_FLOAT, (/     &
     &                 x_dimid, y_dimid, z_dimid, time_dimid /), varid)
            ncstatus=nf90_put_att(ncid, varid, "units", trim(units(i)))
            ncstatus=nf90_def_var_deflate(ncid, varid, shuffle = 1,     &
     &                                   deflate = 1, deflate_level = 5)
            endif
         enddo
         ncstatus=nf90_enddef(ncid)
         ncstatus=nf90_close(ncid)
      else
         ncstatus=nf90_open(filename,NF90_WRITE, ncid)
         ncstatus=nf90_inq_dimid(ncid, "time", time_dimid)
         ncstatus=nf90_inquire_dimension(ncid,time_dimid,len=time)
         start = (/1,1,1,time+1/)
         do i=1,fields
            if (nc_fields(i)) then
            ncstatus=nf90_inq_varid(ncid, trim(var(i)), varid)
               ncstatus=nf90_put_var(ncid, varid,                       &
     &                  buff_final(:,:,(i-1)*levs+1:i*levs),            &
     &                  start=start,count=count)
            endif
        enddo !i
      endif ! kdt
      ncstatus=nf90_close(ncid)
      end subroutine
      end module
