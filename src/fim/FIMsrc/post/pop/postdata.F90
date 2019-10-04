! The purpose of module postdata is to hold all the namelist values for post and pop. An advantage to keeping
! this information in a separate module is that executables such as get_gribout also need to "use" this
! module. If it were contained in post.F90, tons of extra baggage would need to be included at link time 
! (e.g. netcdf, module_control.
!
! Unfortunately the default data spec needs to be "public", because post.F90 and pop.F90 both need to "use" 
! all of it.

module postdata
  implicit none

  public
  save

!JR postnamelist
  integer, parameter :: max_vars=200        ! max number of output vars
  integer, parameter :: max_pathlen=512     ! max length of pathname
  integer, parameter :: max_filelen=32      ! max length of filename
!JR For some reason, the FIM fields created all have a length of 4
!JR The choice of 16 for max length is arbitrary
  integer, parameter :: max_varnamelen=16   ! max length of variable name

  character(len=max_pathlen) :: datadir = ''! path to input FIM data
  character(len=max_pathlen) :: outputdir = '' ! path to output GRIB or netCDF files
  character(len=1) :: input = ''            !JR unknown functionality not ported from pop
  character(len=1) :: output = ''           !JR unknown functionality not ported from pop
  character(len=max_varnamelen) :: mode     ! Unused namelist placebo (needed by pop)
  character(len=16) :: output_fmt           ! 'grib' or 'nc'
  logical :: multiple_output_files          ! placebo: In pop a flag indicating multiple output files
  character(len=max_filelen) :: gribtable = '' ! name of grib table
  integer :: grid_id                        ! GRIB-specific var. related to horizontal resolution
  integer :: mx = -999                      ! number of points in x on output lat-lon grid
  integer :: my = -999                      ! number of points in y on output lat-lon grid
  logical :: latlonfld = .true.             ! flag indicates output grib file is lat-lon (default true)
  integer :: is = -999                      ! Interpolation scheme. Only valid value is 1. Init to invalid
  integer :: vres = -999                    ! number of vertical levels (init to bad value)
  character(len=max_varnamelen) :: var_list(max_vars) = ' ' ! list of GRIB output variables
  integer :: nsmooth_var(max_vars) = 0      ! number of smoothing invocations
!JR t1, t2, and delta_t are unused placebos by FIM when gribout=.true.
  integer :: t1                         ! used only in pop.F90
  integer :: t2, delta_t                ! used only in pop.F90
!JR fimout and gribout are unused placebos by pop when run as a separate executable.
  logical :: fimout = .true.            ! whether to write FIM-style output files directly (default true)
  logical :: gribout = .false.          ! whether to write GRIB files directly (default false)

  namelist /postnamelist/ datadir, outputdir, input, output, output_fmt, multiple_output_files, &
                          gribtable, grid_id, mx, my, latlonfld, is, vres, mode, var_list, &
                          nsmooth_var, t1, t2, delta_t, gribout, fimout

contains

  subroutine post_read_namelist (retcode)
    integer, intent(out) :: retcode

    integer :: j
    integer :: ioerr=0
    integer :: me=0
    logical, save :: post_namelist_read = .false.
    
    retcode = 0

    if (.not. post_namelist_read) then
!sms$insert call nnt_me(me)
      open (11, file="./FIMnamelist", status='old', action='read', iostat=ioerr)
      if (ioerr /= 0) then
        write(*,'(a,i0)') 'post_read_namelist: error opening FIMnamelist by task ',me
        call flush(6)
        retcode = -1
        return
      end if

      read (11, nml=postnamelist, iostat=ioerr)
      if (ioerr /= 0) then
        write(*,'(a,i0)') 'post_init: error reading POSTnamelist from FIMnamelist by task ',me
        call flush(6)
        retcode = -1
        return
      end if
      close (11)

      ! Ensure against old-style input fields by checking for embedded spaces in the field name

      do j=1,max_vars
        if (index (trim (var_list(j)), ' ') /= 0) then
          write(*,'(a,i0,a,a)') 'post_read_namelist: var_list element ',j, ' contains an embedded blank character ', &
                     'which is not allowed.'
          write(*,'(a,a)') 'The first bad input field is:', trim (var_list(j))
          write(*,'(a,a)') 'Perhaps you are using the old-style syntax of "var1 var2 var3 ..." instead of ', &
                     'the new syntax of "var1", "var2", "var3", ...?'
          retcode = -1
          return
        end if
      end do
            
      post_namelist_read = .true.
    end if
    return
  end subroutine post_read_namelist
end module postdata
