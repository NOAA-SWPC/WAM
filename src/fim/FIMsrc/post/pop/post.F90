! The contents of this file were taken from pop.F90, and broken into 4 parts due to needing
! to be called at different times from FIM: 
! 1) things that are done once during the run (post_init)
! 2) things that are done every output time slice (post_init_file)
! 3) things that are done for each output field (post_write_field)
! 4) closing the output GRIB file (post_finalize_file)
!
! Jim Rosinski, Feb, 2011
! Split the subroutine post_init() into two subroutines: post_init_readnl() and post_init_slint(),
! in order to accommodate the initialization logic in dyn_init() subroutine. 
! N. Wang, Sep, 2011

module post
  use module_control, only: curve, fixedgridorder, nvlp, glvl, nip, nvl, nvlp1, &
                            yyyymmddhhmm, ArchvIntvl, totaltime
  use fimnc, only: var_info, set_model_nlevels
  use slint, only: bilinear_init_i2r, bilinear_interp_i2r, tgt_grid
  use postdata, only: post_read_namelist, datadir, outputdir, input, output, output_fmt, max_vars, var_list, &
                      multiple_output_files, gribtable, grid_id, mx, my, latlonfld, is, vres, &
                      mode, nsmooth_var, &
                      max_pathlen, max_varnamelen, gribout, fimout
  use units, only: getunit, returnunit

  implicit none

  private
  save

!JR For some reason, the FIM fields created all have a length of 4
!JR The choice of 16 for max length is arbitrary
  integer :: file_handle = -1               ! File handle for netcdf output: init to bad value
  integer :: nvars                          ! Number of variables to output
  integer :: time = -999                    ! needed by grib routines
  integer :: timebeg = -999                 ! needed by grib routines
  integer :: maxlevs                        ! max number of levels to allocate for data_xyz

  character(len=10) :: jdate                ! julian date (?) used to construct name of grib file

  logical :: post_init_readnl_called = .false.  ! flag indicates part of the post init has been called
  logical :: post_init_slint_called = .false.   ! flag indicates part of the post init has been called

  real, allocatable :: data_xyz(:,:,:)      ! data array sent to grib file (allocated in post_init)
  real :: r2d                               ! convert radians to degrees
  real :: pi                                ! 3.14159...

  logical :: gribfile_is_open = .false. ! indicates a GRIB file is currently open (thus writable)

  ! Public variables

  public :: data_xyz
  public :: fimout
  public :: gribfile_is_open
  public :: gribout
  public :: jdate
  public :: post_init_readnl_called
  public :: post_init_slint_called
  public :: time
  public :: timebeg
  public :: maxlevs

  ! Public routines

  public :: post_init_readnl       ! Called once per run: reads namelist
  public :: post_init_slint       ! Called once per run: initialize the interpolation package
  public :: varindex

contains

!
! post_init must be called once during the model run. It reads the namelist, determines
! some resolution-dependent variables and allocates requisite one-time memory
!
  subroutine post_init_readnl (retcode)
    integer, intent(out) :: retcode         ! return code to caller

    integer :: ret                          ! return code from called routines
    retcode = 0

! Read the postnamelist, then check validity of input
    call post_read_namelist (ret)

    if (ret < 0) then
      write(6,*) 'post_init_readnl: bad return from post_read_namelist'
      call flush(6)
      retcode = -1
      return
    end if

!JR Only is=1 works for now
    if (is /= 1) then
      write(6,*) 'post_init_readnl: only is=1 is currently supported'
      call flush(6)
      retcode = -1
      return
    end if

!JR Ensure "input" and "output" are empty--handle later if needed
    if (input /= '') then
      write(6,*) 'post_init_readnl: input must be empty'
      call flush(6)
      retcode = -1
      return
    end if

    if (output /= '') then
      write(6,*) 'post_init_readnl: output must be empty'
      call flush(6)
      retcode = -1
      return
    end if

!JR Determine number of variables input
    nvars = 0
    do while (var_list(nvars+1) /= ' ' .and. nvars < max_vars)
      nvars = nvars + 1
    end do
    
!JR FIM default is "grib", grid_id=228 => mx=144, my=73
    if (output_fmt == 'grib') then
      call gridid2mxmy (grid_id, mx, my)
    else
      write(0,*) 'post_init_readnl: Only output_fmt="grib" is currently supported'
      call flush(6)
      retcode = -1
      return
    end if
    post_init_readnl_called = .true.

  end subroutine post_init_readnl


  subroutine post_init_slint (retcode)
    integer, intent(out) :: retcode         ! return code to caller

    integer :: ioerr                        ! return value from IO routines
    integer :: ret                          ! return code from called routines
    integer :: unitno                       ! unit number
    character(len=max_pathlen) :: init_file ! full path to file containing icos grid info
    real, allocatable :: llpoints(:,:)

    retcode = 0
    if (post_init_slint_called) then
      write(6,*) 'post_init_slint: must only be invoked once during FIM execution'
      call flush(6)
      retcode = -1
      return
    endif

    maxlevs = max (nvlp1,nvlp)
    allocate (data_xyz(mx, my, maxlevs))
    datadir = datadir(1:len_trim(datadir)) // '/'

    init_file = './icos_grid_info_level.dat'
    allocate (llpoints(nip, 2))

    unitno = getunit ()
    if (unitno < 0) then
      print*,'post_init: getunit failed for ', trim(init_file), '  Stopping'
      stop
    end if

    open (unitno, file=init_file, status='old', action='read', form='unformatted', iostat=ioerr)
    if (ioerr /= 0) then
      write(6,*)'post_init_slint: bad attempt to open ', trim (init_file)
      call flush(6)
      retcode = -1
      return
    end if

    call testglvlheader (unitno, init_file, 'post_init_slint', glvl)
!JR Again, skip fixedgridorder test because we're always reading icos_grid_info_level.dat
!    if (.not. fixedgridorder) then
      call testcurveheader (unitno, init_file, 'post_init_slint', curve)
!    end if
    read (unitno, iostat=ioerr) llpoints(:, 1), llpoints(:, 2)
    if (ioerr == 0) then
!      write(6,*)'post_init_slint: successfully read llpoints from ', trim (init_file)
    else
      write(6,*)'post_init_slint: bad attempt to read ', trim (init_file), ' nelem=', &
                ubound(llpoints,1), ' iostat=', ioerr
      call flush(6)
      retcode = -1
      return
    end if
    close (unitno)
    call returnunit (unitno)
    call bilinear_init_i2r (mx, my, llpoints, nip)

    deallocate (llpoints)

!JR fix later for netcdf
    if (output_fmt == 'nc') then
!JR      call init_cdf_vars (output, file_handle, mx, my, nvlp1, &
!JR                          nfct + 1, 1, nvars, var_list, date_str)
    else if (output_fmt == 'grib') then
!JR The call to set_model_nlevels sets nvl in fimnc.F90. It is needed by the call to var_info
      call set_model_nlevels (nvl)
! Read the gribtable
      call initgrib (gribtable)
    end if

    post_init_slint_called = .true.
    return
  end subroutine post_init_slint

  subroutine gridid2mxmy (gridid, mx, my)
    integer, intent(in) :: gridid
    integer, intent(out) :: mx, my

    if (gridid == 228) then
      mx = 144
      my = 73
    else if (gridid == 45) then
      mx = 288
      my = 145
    else if (gridid == 3) then
      mx = 360
      my = 181
    else if (gridid == 4) then
      mx = 720
      my = 361
    end if
    return
  end subroutine gridid2mxmy

! varindex returns the index of the field in var_list, or -1 if not found
! TODO: devise a more elegant approach
! TODO: ensure this function doesn't cost a lot

  integer function varindex (varname)
    character(len=*), intent(in) :: varname
    integer :: i

    varindex = -1  ! return value if not found
    do i=1,max_vars
      if (varname == var_list(i)) then
        varindex = i
      end if
    end do
    return
  end function varindex
end module post

subroutine post_finalize_file (retcode)

! post_finalize_file closes the open GRIB file

  use post,     only: gribfile_is_open
  use postdata, only: output_fmt

  implicit none
  
  integer, intent(out) :: retcode

  retcode = 0

  if (.not. gribfile_is_open) then
    write(6,*)'post_finalize_file: gribfile already closed'
    call flush(6)
    retcode = -1
    return
  end if
!JR Looks like closegrib knows some "currently open unit number" under the covers
  if (output_fmt == 'grib') then
    call closegrib ()
    gribfile_is_open = .false.
  end if

end subroutine post_finalize_file

subroutine post_init_file (newtime, retcode)

! post_init_file opens the output GRIB or (for future) netCDF file

  use module_control, only: yyyymmddhhmm
  use post,           only: gribfile_is_open, jdate, post_init_slint_called, time, timebeg
  use postdata,       only: max_pathlen, outputdir, output_fmt

  implicit none

  integer, intent(in)        :: newtime
  integer, intent(out)       :: retcode ! return code to caller

  character(len=19)          :: date_str ! unknown: used only for netcdf files
  character(len=6)           :: ahr 
  character(len=max_pathlen) :: gribfile ! full path of grib file to open
  integer                    :: nfct = -1 ! number of forecast times (netcdf only)
  integer                    :: year, month, day, hour, minute, jday
  integer, external          :: iw3jdn

  retcode = 0
  write(6,*) 'post_init_file: newtime=',newtime
! commented on 10/03/2011
!  if (.not. post_init_slint_called) then
!    write(6,*) 'post_init_file: post_init_slint has not yet been called'
!    call flush(6)
!    retcode = -1
!    return
!  end if

  if (gribfile_is_open) then
    write(6,*) 'post_init_file: gribfile is already open--cannot open a new one'
    call flush(6)
    retcode = -1
    return
  end if

  time = newtime
  timebeg = max (time-1, 0)
  write(6,*)'post_init_file: initializing GRIB file for time=', time

! open and init the netCDF or GRIB file

  if (output_fmt == 'grib') then
! get date info from the date string
    read (unit=yyyymmddhhmm(1:4), fmt='(i4)') year
    read (unit=yyyymmddhhmm(5:6), fmt='(i2)') month
    read (unit=yyyymmddhhmm(7:8), fmt='(i2)') day
    read (unit=yyyymmddhhmm(9:10), fmt='(i2)') hour
    read (unit=yyyymmddhhmm(11:12), fmt='(i2)') minute

! create a year 'month-date-hour-minute' date string
    date_str = yyyymmddhhmm(1:4) // '-' // yyyymmddhhmm(5:6) // '-' // yyyymmddhhmm(7:8) // &
      '-' // yyyymmddhhmm(9:10) // ':' // yyyymmddhhmm(11:12) // ':00'

! create the jdate string
    jday = iw3jdn (year, month, day) - iw3jdn (year, 1, 1) + 1
    write (unit=jdate(1:2), fmt='(i2.2)') mod (year, 100) 
    write (unit=jdate(3:5), fmt='(i3.3)') jday 
    write (unit=jdate(6:7), fmt='(i2.2)') hour 

! changed to write i6.6 instead of i3.3 so we can output time resolutions > 3 digits (ex., minutes) as part of the filename
    jdate = jdate(1:7) // '000'
    WRITE(gribfile,'(a,"/",a7,i6.6)')  outputdir(1:len_trim(outputdir)), jdate, time
    write(6,*) 'post_init_file: gribfile: ', trim(gribfile)
    call opengrib (gribfile)
    gribfile_is_open = .true.
  end if

end subroutine post_init_file


subroutine post_write_field (vardata, varname, scalefactor, accum_start, retcode)

! post_write_field writes a single field to the already-open GRIB file

  use fimnc,          only: var_info
  use module_control, only: nip, nvlp1, nvlp, pres_hpa
  use post,           only: data_xyz, gribfile_is_open, jdate, time, timebeg, varindex, maxlevs
  use postdata,       only: max_varnamelen, mx, my, nsmooth_var, output_fmt, var_list, outputdir
  use slint,          only: bilinear_interp_i2r
  
  implicit none
  
  real, intent(in) ::         vardata(*)  ! Data on FIM grid to be written.
! Can be as big as nip*nvlp1. Note collapsing of dims to 1
  character(len=*), intent(in) :: varname ! FIM variable name
  real, intent(in) :: scalefactor         ! scaling factor for GRIB data
  integer, intent(in) :: accum_start      ! beginning time for GRIB data (-1 means use default)
  integer, intent(out) :: retcode         ! return code to send caller

  integer, parameter :: unknown = -999 ! placeholder in netcdf output
  integer :: nlevels   ! number of levels in field
  integer :: v         ! index of varname in var_list
  integer :: i, k      ! loop indices
  character(len=max_varnamelen+2) :: var_name  ! varname maybe with '_B' appended
  character(len=max_varnamelen)   :: units     ! needed only for netcdf output
  character(len=max_varnamelen)   :: var_desc  ! variable description
  real, allocatable :: vardata_scaled(:)

  character(len=128) :: intbinfile
  character(len=6)        :: ahr 
  character(len=8)        :: FMT='(I3.3)'

  retcode = 0
  if (.not. gribfile_is_open) then
    write(6,*)'post_write_field: trying to write field ', trim(varname), ' but gribfile not open'
    call flush(6)
    retcode = -1
    return
  end if

  v = varindex (varname)
  if (v < 0) then
    write(6,*)'post_write_field: no grib equivalent found for ', trim(varname), ' so skipping'
    call flush(6)
    return
  end if

  call var_info (varname, var_desc, units, nlevels, nvlp)

! Scale data for GRIB output if required
  if (scalefactor == 1.) then
    do k=1,nlevels 
      call bilinear_interp_i2r (k, nlevels, vardata, data_xyz) 
    end do
  else
!TODO modify bilinear_interp_i2r_post to do 1 level at a time, avoiding the necessity of
!     allocating a full 3-d temporary array here
    allocate (vardata_scaled(1:nip*nlevels))
    vardata_scaled(1:nip*nlevels) = vardata(1:nip*nlevels) * scalefactor
    do k=1,nlevels 
      call bilinear_interp_i2r (k, nlevels, vardata_scaled, data_xyz) 
    end do
    deallocate (vardata_scaled)
  end if

  do i=1,nsmooth_var(v)
    call smooth (data_xyz, mx, my, nlevels, 0.2)
  end do
!JR Zero out output array above top level needed
  if (nlevels < maxlevs) then
    do k=nlevels+1,maxlevs
      data_xyz(:,:,k) = 0.0
    end do
  end if

  if (output_fmt == 'nc') then
!JR      call write_data (file_handle, date_str, var_list(v), & 
!JR                       var_desc, data_xyz, units, unknown) 
  else if (output_fmt == 'grib') then
!JR Why the test on nlevels > 2???
   if (nlevels > 2) then
      var_name = var_list(v)(1:len_trim(var_list(v)))
      IF(var_name /= "hgtP" .AND. var_name /= "tmpP" .AND. &
         var_name /= "up3P" .AND. var_name /= "vp3P" .AND. &
         var_name /= "oc1P" .AND. var_name /= "oc2P" .AND. &
         var_name /= "bc1P" .AND. var_name /= "bc2P" .AND. &
         var_name /= "so2P" .AND. var_name /= "slfP" .AND. &
         var_name /= "d1sP" .AND. var_name /= "d2sP" .AND. &
         var_name /= "d3sP" .AND. var_name /= "d4sP" .AND. &
         var_name /= "d5sP" .AND. var_name /= "s1sP" .AND. &
         var_name /= "s2sP" .AND. var_name /= "s3sP" .AND. &
         var_name /= "s4sP" .AND. var_name /= "dmsP" .AND. &
         var_name /= "msaP" .AND. var_name /= "p25P" .AND. &
         var_name /= "rh3P" .AND. var_name /= "p10P") THEN
         var_name = var_list(v)(1:len_trim(var_list(v))) // '_B'
      end if
    else
      var_name = var_list(v)(1:len_trim(var_list(v)))
    end if
    write(6,*)'post_write_field: writing field ', trim(var_name)

!JR Some variables are accumulated from a time other than "timebeg"
!JR Flag value -1 says use default

    if (accum_start == -1) then  ! flag value meaning use default
      call writegrib (var_name, mx, my, nlevels, 0, &
        0, data_xyz, jdate, time, timebeg, nvlp, pres_hpa)
    else
      call writegrib (var_name, mx, my, nlevels, 0, &
        0, data_xyz, jdate, time, accum_start, nvlp, pres_hpa)
    end if
    write(6,*)'post_write_field: done writing GRIB field ', trim(varname), ' index=', v
  end if
  return
end subroutine post_write_field
