      module wam_jh_integral
!     Contains code for storing JH values from IDEA physics
!     and summing the array for a global volumetric integral.
      use layout1,   only: lats_node_r, ipt_lats_node_r
      use gg_def,    only: coslat_r, dlat_r, sinlat_r
      use resol_def, only: latr
      use namelist_physics_def, only: output_jh_integral
      use netcdf
      use mpi_def
      implicit none

      real                :: jh_nh_integral, jh_sh_integral
      real, allocatable   :: jh_global(:,:,:)
      character(20), parameter :: filename="jh.nc4"

      public  :: jh_nh_integral, jh_sh_integral, jh_global
      private :: filename

      contains
!
!
!
      subroutine jh_integral_init(ngptc,nblck)
      implicit none

      integer, intent(in) :: ngptc, nblck

      if (.not.allocated(jh_global)) then
        allocate (jh_global(ngptc,nblck,lats_node_r))
      endif

      if (output_jh_integral) call start_jh_output

      end subroutine jh_integral_init
!
!
!
      subroutine jh_integral_zero()
      implicit none

! Local variables
      real, parameter     :: zero=0.0

      jh_global = zero

      end subroutine jh_integral_zero
!
!
!
      subroutine do_jh_integral(global_lats_r,lonsperlar,nblck,ngptc,me)
      implicit none

      integer, intent(in) :: global_lats_r(latr)
      integer, intent(in) :: lonsperlar(latr)
      integer, intent(in) :: nblck, ngptc, me
! Local variables
      integer             :: lan, lat, lons_lat
      integer             :: i, n, ierr
      real                :: sh_integral, nh_integral
      real                :: sh_integral_sum, nh_integral_sum
      real                :: pi, volfac, dlat

      pi = atan(1.0)*4.0
      sh_integral = 0.0
      nh_integral = 0.0

      do lan=1,lats_node_r
         lat = global_lats_r(ipt_lats_node_r-1+lan)
         lons_lat = lonsperlar(lat)
         if (lat .lt. latr) then
            dlat = abs(asin(sinlat_r(lat))-asin(sinlat_r(lat+1)))
         else
            dlat = abs(asin(sinlat_r(lat))+pi*0.5)
         endif
         volfac = (2*pi/lons_lat) * dlat * coslat_r(lat) * 1.0e-09

         if (sinlat_r(lat).lt.0) then
            do n=1,nblck
               do i=1,ngptc
                  sh_integral = sh_integral + jh_global(i,n,lan) *      &
     &                           volfac
               end do !i
            end do !n
         else
            do n=1,nblck
               do i=1,ngptc
                  nh_integral = nh_integral + jh_global(i,n,lan) *      &
     &                           volfac
               end do !i
            end do !n
         end if !hemispheric
      end do !lan

      call mpi_allreduce(sh_integral, jh_sh_integral, 1, mpi_real8,     &
     &                    mpi_sum, mpi_comm_all, ierr)

      call mpi_allreduce(nh_integral, jh_nh_integral, 1, mpi_real8,     &
     &                    mpi_sum, mpi_comm_all, ierr)
      call mpi_barrier(mpi_comm_all, ierr)

      end subroutine do_jh_integral
!
      subroutine start_jh_output()
      implicit none

! Local variables
      integer           :: ncstatus, ncid, time_dimid
      integer           :: jh_sh_varid, jh_nh_varid

! create
      ncstatus=nf90_create(filename, NF90_NETCDF4, ncid )

! dimensions
      ncstatus=nf90_def_dim( ncid, "time", nf90_unlimited, time_dimid )
      ncstatus=nf90_put_att( ncid, time_dimid, "axis", "T")

! variables
      ncstatus=nf90_def_var( ncid, "jh_sh", NF90_FLOAT, (/time_dimid/), &
     &                        jh_sh_varid)
      ncstatus=nf90_put_att( ncid, jh_sh_varid,"long_name",             &
     &                        "SH Joule Heating Integral" )
      ncstatus=nf90_put_att( ncid, jh_sh_varid, "units", "J/s" )

      ncstatus=nf90_def_var( ncid, "jh_nh", NF90_FLOAT, (/time_dimid/), &
     &                        jh_nh_varid)
      ncstatus=nf90_put_att( ncid, jh_nh_varid,"long_name",             &
     &                        "NH Joule Heating Integral" )
      ncstatus=nf90_put_att( ncid, jh_nh_varid, "units", "J/s" )

! close
      ncstatus=nf90_enddef(ncid)
      ncstatus=nf90_close(ncid)

      end subroutine start_jh_output
!
!
!
      subroutine write_jh_output(kdt)
      implicit none

      integer, intent(in) :: kdt
! Local variables
      integer             :: ncid, ncstatus, varid
      integer             :: start(1)

      start = (/kdt/)

      ncstatus=nf90_open(filename, NF90_WRITE, ncid)

      ncstatus=nf90_inq_varid(ncid, "jh_sh", varid)
      ncstatus=nf90_put_var(ncid, varid, jh_sh_integral, start=start)

      ncstatus=nf90_inq_varid(ncid, "jh_nh", varid)
      ncstatus=nf90_put_var(ncid, varid, jh_nh_integral, start=start)

      ncstatus=nf90_close(ncid)

      end subroutine write_jh_output
!
      end module wam_jh_integral
