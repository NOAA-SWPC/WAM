! input forcing parameter module
      module wam_ifp_mod

      use comio
      use wam_ifp_class
      use mpi_def, only: MPI_COMM_ALL, MPI_INFO_NULL
      use layout1, only: me

      implicit none

      type(param_t) :: params
      type(farr_t)  :: farr
! Legacy variables
      real    :: f107_fix, f107d_fix, kp_fix
      real    :: swpcf107_fix, swpcf107d_fix, swpckp_fix
!

      contains

      subroutine read_ifp

      character(len=19), parameter :: filename = "input_parameters.nc"
      class(COMIO_T), pointer :: io  => null()
      integer, parameter      :: fmt =  COMIO_FMT_PNETCDF
      integer, pointer :: dims(:)

      io => COMIO_T(fmt=fmt, comm=MPI_COMM_ALL, info=MPI_INFO_NULL)

      call dealloc()
      call check_write_lock()
      call manage_read_lock(.true.)

      call io % open(filename, "r")
      call io % description("skip", params % skip)
      call io % description("ifp_interval", params % ifp_interval)

      call io % domain("f107", dims)
      call alloc(dims(1))
      call io % read("f107",  farr % f107)
      call io % read("f107d", farr % f107d)
      call io % read("kp",    farr % kp)
      call io % read("kpa",   farr % kpa)
      call io % read("nhp",   farr % nhp)
      call io % read("nhpi",  farr % nhpi)
      call io % read("shp",   farr % shp)
      call io % read("shpi",  farr % shpi)
      call io % read("swden", farr % swden)
      call io % read("swang", farr % swang)
      call io % read("swvel", farr % swvel)
      call io % read("swbz",  farr % swbz)
      call io % read("swbt",  farr % swbt)
      if ( me .eq. 0 ) then
          write(6,*) farr % f107(1)
          write(6,*) farr % f107d(1)
          write(6,*) farr % kp(1)
          write(6,*) farr % kpa(1)
          write(6,*) farr % nhp(1)
          write(6,*) farr % nhpi(1)
          write(6,*) farr % shp(1)
          write(6,*) farr % shpi(1)
          write(6,*) farr % swden(1)
          write(6,*) farr % swang(1)
          write(6,*) farr % swvel(1)
          write(6,*) farr % swbz(1)
          write(6,*) farr % swbt(1)
      end if
      call io % close()

      call manage_read_lock(.false.)

      end subroutine read_ifp

      subroutine check_write_lock
      character(len=22), parameter :: filename = "input_parameters.wlock"
      logical :: not_ready
      integer :: iostat

      not_ready = .true.
      do while (not_ready)
        inquire(file=filename, exist=not_ready, iostat=iostat)
        if (not_ready) call sleep(1)
      end do

      end subroutine check_write_lock

      subroutine manage_read_lock(create)
      logical, intent(in) :: create

      character(len=26), parameter :: filename = "input_parameters.rlock.wam"
      character(len=29) :: lockfile
      integer, parameter :: unit = 79

      write (lockfile, "(A22,I0.3)") filename,me
      open(unit, file=lockfile, status="replace", action="write")
      if (create) then
        close(unit)
      else
        close(unit, status="delete")
      end if

      end subroutine manage_read_lock

      subroutine alloc(dim)
        integer, intent(in) :: dim

        if (.not.allocated(farr%f107))  allocate(farr%f107 (dim))
        if (.not.allocated(farr%f107d)) allocate(farr%f107d(dim))
        if (.not.allocated(farr%kp))    allocate(farr%kp   (dim))
        if (.not.allocated(farr%kpa))   allocate(farr%kpa  (dim))
        if (.not.allocated(farr%nhp))   allocate(farr%nhp  (dim))
        if (.not.allocated(farr%nhpi))  allocate(farr%nhpi (dim))
        if (.not.allocated(farr%shp))   allocate(farr%shp  (dim))
        if (.not.allocated(farr%shpi))  allocate(farr%shpi (dim))
        if (.not.allocated(farr%swden)) allocate(farr%swden(dim))
        if (.not.allocated(farr%swvel)) allocate(farr%swvel(dim))
        if (.not.allocated(farr%swang)) allocate(farr%swang(dim))
        if (.not.allocated(farr%swbz))  allocate(farr%swbz (dim))
        if (.not.allocated(farr%swbt))  allocate(farr%swbt (dim))

      end subroutine alloc

      subroutine dealloc()

        if (allocated(farr%f107))  deallocate(farr%f107)
        if (allocated(farr%f107d)) deallocate(farr%f107d)
        if (allocated(farr%kp))    deallocate(farr%kp)
        if (allocated(farr%kpa))   deallocate(farr%kpa)
        if (allocated(farr%nhp))   deallocate(farr%nhp)
        if (allocated(farr%nhpi))  deallocate(farr%nhpi)
        if (allocated(farr%shp))   deallocate(farr%shp)
        if (allocated(farr%shpi))  deallocate(farr%shpi)
        if (allocated(farr%swden)) deallocate(farr%swden)
        if (allocated(farr%swvel)) deallocate(farr%swvel)
        if (allocated(farr%swang)) deallocate(farr%swang)
        if (allocated(farr%swbz))  deallocate(farr%swbz)
        if (allocated(farr%swbt))  deallocate(farr%swbt)

      end subroutine dealloc
! legacy code below, not sure this is still needed
!==========================================================
! Below two service subs to disable "read_wam_f107_kp_txt"
! during model tune-ups
!==========================================================
      subroutine fix_spweather_data
!=======================================================================
!VAY 2016: This is temporal "substitue" for "read_wam_f107_kp_txt"
!    with fixed Kp and F107 data to work with long-term WAM run
! TO DO "advance_solar" KP-F107 drivers with WAM calendar
!=======================================================================
      swpcf107_fix = 100.
      swpckp_fix   = 1.
      swpcf107d_fix = swpcf107_fix

      f107_fix = 100.
      kp_fix   = 1.
      f107d_fix = f107_fix
      end subroutine fix_spweather_data
!
      subroutine read_spweather_real_data
!=======================================================================
!VAY 2016: This is temporal "substitue" for "read_wam_f107_kp_txt"
!    with fixed Kp and F107 data to work with long-term WAM run
! TO do "advance_solar" KP-F107 drivers with WAM calendar
!=======================================================================
      f107_fix = 100.
      kp_fix   = 1.
      f107d_fix = f107_fix

      end subroutine read_spweather_real_data
!

      end module wam_ifp_mod
