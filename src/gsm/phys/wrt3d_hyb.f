      subroutine wrt3d_hyb(dt6dt,kdt,global_lats_r,lonsperlar,nblck,
     &                     idate)
      use resol_def
      use layout1
      use namelist_physics_def
      use machine
      implicit none
!!
      integer, parameter :: numv = 2
      integer, parameter :: ioproc=0
      integer, parameter :: nwave = 37
!!
      integer, intent(in) :: kdt,nblck
      real (kind=kind_rad),dimension(nwave,2),
     &                     intent(in) :: dt6dt
      integer, intent(in), dimension(latr) :: global_lats_r, lonsperlar
      integer, dimension(4), intent(in) :: idate
! Local variables
      integer :: nv,k,lan,lon,i,iblk,lons_lat,njeff,il,lat
      real (kind=kind_io4) :: wrkga(lonr,latr)
      real (kind=kind_io8) :: rtime
      real (kind=kind_io8) :: glolal(lonr,lats_node_r)
      real (kind=kind_io8) :: buffo(lonr,lats_node_r)
      integer:: kmsk(lonr*latr)
      integer, dimension(nodes) :: lats_nodes_r
!..........................................................
!     temperature tendencies
!
      kmsk = 0
      do nv=1,numv
        do k=1,nwave
          if(me.eq.ioproc) call ncwrite(nv, k, kdt, dt6dt(k,nv), idate)
        enddo
      enddo

      return
      end subroutine wrt3d_hyb

!!

      subroutine ncwrite(nv, k, kdt, workga, idate)
      use netcdf
      use layout1
      use gg_def
      use machine
      use resol_def
      implicit none
!
      integer, intent(in) :: nv, k, kdt
      integer, intent(in), dimension(4) :: idate
      real(kind=kind_rad), intent(in) :: workga
! Local variables
      real :: pi

      integer :: id, x_dimid, x_varid, y_dimid, y_varid, t_dimid,
     &           t_varid, varid, t_len, i, adj, z_dimid, z_varid
      integer, parameter :: numv = 2
      integer, parameter :: num_variables = 2
      integer, dimension(2) :: dimids
      integer, dimension(levs) :: nlevs

      character(len=*), dimension(num_variables), parameter :: vars  =
     &           (/"flux","rlmeuv"/)
      character(len=*), dimension(num_variables), parameter :: lvars =
     &           (/"Parameterized EUV Flux", "Wavelength"/)
      character(len=*), dimension(num_variables), parameter :: units =
     &           (/"unknown","cm"/)
      character(len=8), parameter :: fname = "dt6dt.nc"

      character(len=33) :: time_units

      logical :: exists
!
      inquire( file = fname, exist=exists )

      if (.not. exists) then
        write(time_units, "(A14,I0.4,A1,I0.2,A1,I0.2,A1,I0.2,A6)"),
     &   "minutes since ", idate(4), "-", idate(2), "-", idate(3),
     &   " ", idate(1),":00:00"
        call check(nf90_create(fname, nf90_clobber, id))

        call check(nf90_def_dim(id,"nwave",37, z_dimid))
        call check(nf90_def_dim(id,"time",NF90_UNLIMITED,t_dimid))

        call check(nf90_def_var(id,"time",NF90_INT,t_dimid,t_varid))
        call check(nf90_put_att(id,t_varid,"long_name","time"))
        call check(nf90_put_att(id,t_varid,"units",time_units))
        call check(nf90_put_att(id,x_varid,"axis","T"))

        dimids = (/ z_dimid, t_dimid /)

        do i=1,numv
          call check(nf90_def_var(id,vars(i),NF90_REAL,dimids,varid))
          call check(nf90_put_att(id,varid,"units",units(i)))
          call check(nf90_put_att(id,varid,"long_name", lvars(i)))
        end do

        call check(nf90_enddef(id))

        call check(nf90_close(id))
      end if

      call check(nf90_open(fname, NF90_WRITE, id))

      call check(nf90_inq_varid(id, vars(nv), varid))
      call check(nf90_put_var(id,varid,(/workga/),start=(/k,kdt+1/),
     &                              count=(/1,1/)))

      if (nv .eq. 1) then
        call check(nf90_inq_varid(id, "time", varid))
        call check(nf90_put_var(id,varid,kdt,start=(/kdt+1/)))
      end if

      call check(nf90_close(id))

      end subroutine ncwrite

!!

      subroutine check(istatus)

      use netcdf

      implicit none

      integer, intent(in) :: istatus

      if (istatus /= nf90_noerr) then
        write(*,*) TRIM(ADJUSTL(nf90_strerror(istatus)))
      end if

      end subroutine check
