!=============================================================
! Post processor utility
!
!
! Main packages: slint (spherical linear interpolation),
!                vlint (vertical linear interpolation),
!                fimnc (FIM netCDF utility routines),
!                gribio (FIM grib utility routines) 
!
! Ning Wang, March 2007 
!
!
!=============================================================

      program pop
      USE fimnc
      use module_control,only: TotalTime,ArchvIntvl,curve,NumCacheBlocksPerPE,&
                               PrintMAXMINtimes,FixedGridOrder,nvlp,pres_hpa,&
                               TimingBarriers,glvl,nip,nvl,control,&
                               yyyymmddhhmm,ArchvTimeUnit,numphr
      use postdata, only: post_read_namelist, datadir, outputdir, input, output, output_fmt, max_vars, var_list, &
                          multiple_output_files, gribtable, grid_id, mx, my, latlonfld, is, vres, &
                          mode, nsmooth_var, &
                          t1, t2, delta_t
      USE slint, ONLY: bilinear_init_i2r, bilinear_interp_i2r, tgt_grid

      IMPLICIT NONE

! define a derived type for variable data
      TYPE array_pointer
         REAL, DIMENSION(:,:), POINTER :: p
      END TYPE array_pointer

      INTEGER                 :: nverlvs  !The number of vertical levels defined in the Makefile
      INTEGER                 :: iargc    ! command line argument count
      INTEGER                 :: nvars, nfct
      CHARACTER(len=char_len) :: init_file
      CHARACTER(len=6)        :: ahr 
      CHARACTER(len=8)        :: FMT='(I3.3)'
      CHARACTER(len=80)       :: FMT2
      INTEGER                 :: file_handle 
      INTEGER                 :: var_num, nlevels
      INTEGER                 :: year, month, day, hour, minute, jday, IW3JDN, is2Dvar
      INTEGER                 :: nt, i, j, k , idx, ierr 
      REAL, ALLOCATABLE       :: data_xyz(:,:,:), vardata(:), vardata_n(:)
      REAL, ALLOCATABLE       :: data_xyz_var(:,:,:)
      REAL, ALLOCATABLE       :: data_xyz_pr(:,:,:)
      REAL, ALLOCATABLE       :: llpoints(:,:)
      CHARACTER(len=char_len) :: var_name, var_description, units
      CHARACTER(len=19      ) :: date_str, date_str2
      CHARACTER(len=10      ) :: jdate 
      CHARACTER(len=char_len) :: gribfile
      REAL                    :: missing_value, r2d, pi
      integer :: ioerr
      integer :: ret      ! returned from subroutine call
      integer :: maxlevs  ! max number of levels for array allocation

! check and get the commandline arguments
      IF (iargc() .NE. 0) THEN
        WRITE(0,*) 'Usage: pop'
        WRITE(0,*) 'Note: Make sure that namelist file pop.nl is in the current directory'
        STOP
      END IF

! set up control variables
!JR Subroutine "control" reads in postnamelist, but it is opened and read in again here
!TODO: Only read it once

      call control ()
      nverlvs = nvl
      delta_t=ArchvIntvl
      t2=TotalTime
      nsmooth_var=0
      t1=0
      call post_read_namelist (ret)
      if (ret < 0) then
        write(6,*) 'pop: bad return from post_read_namelist'
        stop
      end if

!JR Determine number of variables input
      nvars = 0
      do while (var_list(nvars+1) /= ' ' .and. nvars < max_vars)
        nvars = nvars + 1
      end do

!JR is=1 in default FIM => horiz. interp and no vert. interp
      IF (is == 3 .AND. nvars > 3) THEN
        WRITE(0,*) 'Only allow maximum 3 variables for vertical cross sections'
        STOP
      ENDIF

! get date info from the date string
      READ(UNIT=yyyymmddhhmm(1:4), FMT='(I4)') year
      READ(UNIT=yyyymmddhhmm(5:6), FMT='(I2)') month
      READ(UNIT=yyyymmddhhmm(7:8), FMT='(I2)') day
      READ(UNIT=yyyymmddhhmm(9:10), FMT='(I2)') hour
      READ(UNIT=yyyymmddhhmm(11:12), FMT='(I2)') minute

! create a year 'month-date-hour-minute' date string
      date_str = yyyymmddhhmm(1:4) // "-" // yyyymmddhhmm(5:6) // "-" // yyyymmddhhmm(7:8) // "-" // &
                 yyyymmddhhmm(9:10) // ":" // yyyymmddhhmm(11:12) // ":00"
      date_str2 = date_str
     
! create the jdate string
      jday = IW3JDN(year,month,day) - IW3JDN(year,1, 1) + 1
      WRITE(UNIT=jdate(1:2), FMT='(I2.2)') MOD (year, 100) 
      WRITE(UNIT=jdate(3:5), FMT='(I3.3)') jday 
      WRITE(UNIT=jdate(6:7), FMT='(I2.2)') hour 
      jdate = jdate(1:7) // '000'

! compute the number of forecast time and number of icosahedral grid
      nfct = (t2 - t1) / delta_t 

!JR FIM default is "grib"
!JR grid_id=228 => mx=144, my=73
      IF (output_fmt == "grib") THEN
        CALL gridid2mxmy(grid_id, mx, my)
      ENDIF

      IF (is == 0) THEN
        mx = 1024
        my = nip / 1024 + 1
      ENDIF

      IF (is == 0) THEN
        ALLOCATE(vardata(nip * (nverlvs+1)))
        ALLOCATE(vardata_n(mx * my * (nverlvs+1)))
      ELSE IF (is == 1) THEN
!JR The following is a HACK to get pop to work when number of pressure levels exceeds number of model levels
        maxlevs = max (nverlvs+1,nvlp)
        ALLOCATE(vardata(nip * maxlevs))
        ALLOCATE(data_xyz(mx, my, maxlevs))
      ELSE IF (is == 2) THEN
        ALLOCATE(vardata(nip * (nverlvs+1)))
        ALLOCATE(data_xyz_var(mx,my,nverlvs+1))
        ALLOCATE(data_xyz_pr(mx,my,nverlvs+1))
        ALLOCATE(data_xyz(mx, my, nvlp))
      ELSE IF (is == 3) THEN
        ALLOCATE(vardata(nip * (nverlvs+1)))
        ALLOCATE(data_xyz_var(mx,my,nverlvs+1))
        ALLOCATE(data_xyz_pr(mx,my,nverlvs+1))
        ALLOCATE(data_xyz(mx, my, vres))
      ENDIF

      missing_value = -99.0

! open and init the netCDF or GRIB file
      IF(is == 0) THEN
        CALL set_model_nlevels(nverlvs)
        CALL initgrib(gribtable)
      ELSE IF(is == 1) THEN
        IF (output_fmt == "nc") THEN
          CALL init_cdf_vars(output, file_handle, mx, my, nverlvs+1, nfct + 1, 1, nvars, var_list, date_str, nvlp)
        ELSEIF (output_fmt == "grib") THEN
          CALL set_model_nlevels(nverlvs)
          CALL initgrib(gribtable)
          IF (.not. multiple_output_files) CALL opengrib(output)
        END IF
      ELSE IF (is == 2) THEN
        IF (output_fmt == "nc") THEN
          CALL init_cdf_vars_v(output, file_handle, mx, my, nvlp, nfct + 1, 1, nverlvs, nvars, var_list, date_str, nvlp)
        ELSEIF (output_fmt == "grib") THEN
          CALL set_model_nlevels(nverlvs)
          CALL initgrib(gribtable)
          IF (.not. multiple_output_files) CALL opengrib(output)
        END IF
      ELSE IF (is == 3) THEN
        IF (output_fmt == "nc") THEN
          CALL init_cdf_vars_v(output, file_handle, mx, my, vres, nfct + 1, 1, nverlvs, nvars, var_list, date_str, nvlp)
        ELSEIF (output_fmt == "grib") THEN
          CALL set_model_nlevels(nverlvs)
          CALL initgrib(gribtable)
          IF (.not. multiple_output_files) CALL opengrib(output)
        END IF
      ENDIF
    
      datadir = datadir(1:LEN_TRIM(datadir)) // '/'

! if interpolation scheme is not 0, init the horizontal interpolation.
      IF (is /= 0) THEN
        IF (FixedGridOrder) THEN
          FMT2 = '(a,"latlonIJ.dat")'
          write(init_file,FMT2) datadir(1:LEN_TRIM(datadir))
          ALLOCATE(llpoints(nip, 2))
          OPEN (10, file=init_file, action='read', form='unformatted', iostat=ioerr)
          if (ioerr == 0) then
            write(6,*)'pop: successfully opened init_file=', trim(init_file)
          else
            write(6,*)'pop: failed to open init_file=', trim(init_file)
          end if
          call TestGlvlHeader(10,init_file,'pop',glvl)
          READ (10, iostat=ioerr) llpoints(:, 1), llpoints(:, 2)
          if (ioerr /= 0) then
            write(6,*)'pop: bad attempt to read ', trim (init_file), 'nelem=', &
                      ubound(llpoints,1), ' iostat=', ioerr
          end if
          CLOSE(10)
        ELSE
          init_file = './icos_grid_info_level.dat'
          ALLOCATE(llpoints(nip, 2))
          OPEN (10, file=init_file, action='read', form='unformatted', iostat=ioerr)
          if (ioerr == 0) then
            write(6,*)'pop: successfully opened init_file=', trim(init_file)
          else
            write(6,*)'pop: failed to open init_file=', trim(init_file)
          end if
          call TestGlvlHeader (10,init_file,'pop',glvl )
          call TestCurveHeader(10,init_file,'pop',curve)
          READ (10, iostat=ioerr) llpoints(:, 1), llpoints(:, 2)
          if (ioerr == 0) then
!            write(6,*)'pop: successfully read llpoints nelem=',ubound(llpoints,1)
          else
            write(6,*)'pop: bad attempt to read ', trim (init_file), 'nelem=', &
                      ubound(llpoints,1), ' iostat=', ioerr
          end if
          CLOSE(10)
        ENDIF
!        write(6,*)'JR pop: calling bilinear_init_i2r llpoints=', llpoints(:,1)
        CALL bilinear_init_i2r(mx, my, llpoints, nip)
        
        DEALLOCATE(llpoints)
      ENDIF
    
      IF (latlonfld) THEN
        pi = 4.0*ATAN(1.0)
        r2d = 180.0 / pi 
        ALLOCATE(llpoints(nip, 2))
        FMT2 = '(a,"latlonIJ.dat")'
        write(init_file,FMT2) datadir(1:LEN_TRIM(datadir))
        OPEN(10,file=init_file,status='old',form='unformatted')
        call TestGlvlHeader(10,'latlonIJ.dat','pop',glvl)
        READ(10) llpoints(:, 1), llpoints(:, 2)
        CLOSE(10)
        llpoints(1:nip, 1) = llpoints(1:nip, 1) * r2d
        llpoints(1:nip, 2) = (llpoints(1:nip, 2) - pi) * r2d
      END IF

      IF (is == 0) THEN  ! no interpolation, native grid can only be saved in GRIB file format at this point.
        DO nt = t1, t2, delta_t
          IF (output_fmt.eq."grib" .and. multiple_output_files) THEN
            WRITE(ahr, FMT) nt
            gribfile = outputdir(1:LEN_TRIM(outputdir)) // '/' // jdate  // ahr
            CALL opengrib(gribfile)
          ELSE
             PRINT*, 'Native grid can only be saved in GRIB file format.'
             STOP
          ENDIF
          DO var_num = 1, nvars ! for each variable
            var_name = var_list(var_num)(1:LEN_TRIM(var_list(var_num)))
            CALL read_1var(nt, delta_t, datadir, vardata, var_name, nip, nverlvs + 1,ArchvTimeUnit, nvlp)
            CALL var_info(var_list(var_num), var_description, units, nlevels, nvlp)
            DO i = 1, nlevels  
              DO j = 1, nip
                vardata_n(mx*my*(i-1)+j) = vardata(i + (j - 1) * (nverlvs + 1))
              END DO
              Do j = nip + 1, mx * my
                vardata_n(mx*my*(i-1)+j) = vardata(i + (nip - 1) * (nverlvs + 1))
              END DO
            ENDDO
            IF (nlevels > 2) THEN
              var_name = var_list(var_num)(1:LEN_TRIM(var_list(var_num))) // '_B'
            ELSE
              var_name = var_list(var_num)(1:LEN_TRIM(var_list(var_num)))
            ENDIF
            CALL writegrib(var_name, nip, 1, nlevels, glvl, curve, vardata_n, jdate, nt, &
                           tba(nt,t1,delta_t,var_name), nvlp, pres_hpa)
          ENDDO
          IF (latlonfld) THEN
            vardata_n(1:nip) = llpoints(1:nip,1) 
            vardata_n(nip+1:mx * my) = 0.0 
            var_name = "LAT"
            CALL writegrib(var_name, nip, 1, 1, glvl, curve, vardata_n, jdate, t1, 0, nvlp, pres_hpa)
            vardata_n(1:nip) = llpoints(1:nip,2) 
            vardata_n(nip+1:mx * my) = 0.0 
            var_name = "LON"
            CALL writegrib(var_name, nip, 1, 1, glvl, curve, vardata_n, jdate, t1, 0, nvlp, pres_hpa)
          ENDIF
          IF (output_fmt.eq."grib" .and. multiple_output_files) CALL closegrib()
          PRINT "('*'$)"  ! progress '*'
        ENDDO
      ELSE IF (is == 1) THEN ! Horizontal interpolation
        DO nt = t1, t2, delta_t
          IF (output_fmt.eq."grib" .and. multiple_output_files) THEN
            WRITE(ahr, FMT) nt
            gribfile = outputdir(1:LEN_TRIM(outputdir)) // '/' // jdate // ahr
            open (unit=1, file=gribfile, status='old', action='read', iostat=ioerr)
            close (1)
            if (ioerr == 0) then
              write(6,*)'pop: GRIB file ', trim(gribfile), ' already exists: stopping...'
              stop 999
            end if
            CALL opengrib(gribfile)
          ENDIF
          DO var_num = 1, nvars ! for each variable
            IF (input /= "") THEN
              CALL read_1var_direct(input, nip, nverlvs, vardata,ArchvTimeUnit, nvlp)
              nlevels = nverlvs
            ELSE
              var_name = var_list(var_num)(1:LEN_TRIM(var_list(var_num)))
              CALL read_1var(nt, delta_t, datadir, vardata, var_name, nip, nverlvs + 1,ArchvTimeUnit, nvlp)
              CALL var_info(var_list(var_num), var_description, units, nlevels, nvlp)
            ENDIF
!            write(6,*)'JR pop: calling bilinear_interp_i2r for field ', trim(var_name)
            DO k=1,nlevels 
              CALL bilinear_interp_i2r(k, nlevels, vardata, data_xyz) 
            END DO !level loop
!            write(6,*)'JR pop after bilinear_interp_i2r: nlevels,vardata,data_xyz=', &
!                      nlevels, vardata(1), data_xyz(1,1,1)
            DO i = 1, nsmooth_var(var_num)
              CALL smooth(data_xyz,mx,my,nlevels,0.2)
            END DO

              
            IF (nlevels < nverlvs + 1) THEN
              data_xyz(:,:,nlevels + 1:) = 0.0
            END IF
            IF (output_fmt == "nc") THEN
              CALL write_data (file_handle, date_str, var_list(var_num), & 
              var_description, data_xyz, units, nt / delta_t) 
            ELSEIF (output_fmt == "grib") THEN
              IF (nlevels > 2) THEN
                var_name = var_list(var_num)(1:LEN_TRIM(var_list(var_num)))
!JR The "_B" means the data are on native levels rather than pressure levels
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
                  var_name = var_list(var_num)(1:LEN_TRIM(var_list(var_num))) // '_B'
                ENDIF
              ELSE
                var_name = var_list(var_num)(1:LEN_TRIM(var_list(var_num)))
              ENDIF
              CALL writegrib(var_name, mx, my, nlevels, 0, 0, data_xyz, jdate, nt, &
                             tba(nt,t1,delta_t,var_name),nvlp, pres_hpa)
            ENDIF
          END DO 
          CALL mmdddateadj(date_str, delta_t)  
          IF (output_fmt.eq."grib" .and. multiple_output_files) CALL closegrib()
          PRINT "('*'$)"  ! progress '*'
        END DO 
      ELSE ! vertical interpolation
        var_list(nvars+1) = "pr3D"
        DO nt = t1, t2, delta_t
          IF (output_fmt.eq."grib" .and. multiple_output_files) THEN
            WRITE(ahr, FMT) nt
            gribfile = outputdir(1:LEN_TRIM(outputdir)) // '/' // jdate // ahr
            CALL opengrib(gribfile)
          ENDIF

          ! interpolate pressure to latlon grid 
          var_name = var_list(var_num)(1:LEN_TRIM(var_list(var_num)))
          CALL read_1var(nt, delta_t, datadir, vardata, var_name, nip, nverlvs + 1,ArchvTimeUnit,nvlp)
          CALL var_info(var_list(nvars+1), var_description, units, nlevels,nvlp)
          DO k = 1, nlevels  
            CALL bilinear_interp_i2r (k, nlevels, vardata, data_xyz_pr) 
          END DO

          ! interpolate the specified variables to latlon grid 
          DO var_num = 1, nvars ! for each variable
            var_name = var_list(var_num)(1:LEN_TRIM(var_list(var_num)))
            CALL read_1var(nt, delta_t, datadir, vardata, var_name, nip, nverlvs + 1,ArchvTimeUnit, nvlp)
            CALL var_info(var_list(var_num), var_description, units, nlevels, nvlp)
            DO k = 1, nlevels 
              CALL bilinear_interp_i2r (k, nlevels, vardata, data_xyz_var) 
            END DO 
            ! interpolate to vertical plane
            data_xyz = 0.0
            IF (is == 2) THEN
              IF (is2Dvar(var_list(var_num)) == 1) THEN
                var_name = var_list(var_num)(1:LEN_TRIM(var_list(var_num)))
                DO i = 1, nsmooth_var(var_num)
                  CALL smooth(data_xyz_var,mx,my,1,0.2)
                END DO
                CALL writegrib(var_name, mx, my, 1, 0, 0, data_xyz_var, jdate, nt, &
                               tba(nt,t1,delta_t,var_name),nvlp,pres_hpa) 
                CYCLE
              ENDIF
              CALL vlint2coor(mx, my, nlevels, nverlvs + 1, data_xyz_pr, data_xyz_var, data_xyz, pres_hpa, nvlp)  
            ELSE IF (nlevels == nverlvs) THEN
              IF (mode == "step") THEN
                CALL v_interp(mx, my, nlevels, data_xyz_pr, data_xyz_var, vres, data_xyz, missing_value, 1)  
              ELSE IF (mode == "linear") THEN
                CALL v_interp(mx, my, nlevels, data_xyz_pr, data_xyz_var, vres, data_xyz, missing_value, 0)  
              ENDIF
            ELSEIF (nlevels == nverlvs + 1) THEN ! variables defined on interface
              CALL v_interp_lvlvar(mx, my, nlevels, data_xyz_pr, data_xyz_var, vres, data_xyz, missing_value, 0)
            ENDIF
            DO i = 1, nsmooth_var(var_num)
              IF (is == 2) THEN
                CALL smooth(data_xyz,mx,my,nvlp,0.2)
              ELSE if (is == 3) THEN
                CALL smooth(data_xyz,mx,my,vres,0.2)
              END IF
            END DO

! write data to the netCDF file, the variable interpolated
            IF (output_fmt == "nc") THEN
              CALL write_data (file_handle, date_str, var_list(var_num), & 
              var_description, data_xyz, units, nt / delta_t) 
            ELSEIF (output_fmt == "grib") THEN
              var_name = var_list(var_num)(1:LEN_TRIM(var_list(var_num))) // '_P'
              IF (is == 2) THEN
                CALL writegrib(var_name, mx, my, nvlp, 0, 0, data_xyz, jdate, nt, &
                               tba(nt,t1,delta_t,var_name),nvlp,pres_hpa)
              ELSE
                CALL writegrib(var_name, mx, my, vres, 0, 0, data_xyz, jdate, nt, &
                               tba(nt,t1,delta_t,var_name),nvlp,pres_hpa)
              ENDIF
            ENDIF
          END DO
          ! write pressure field to the nc file
          IF (is == 3) THEN
            DO k=1,vres 
              IF (k <= nverlvs + 1) THEN
                data_xyz(1:mx, 1:my, k) = data_xyz_pr(1:mx, 1:my, k)
              ELSE
                data_xyz(1:mx, 1:my, k) = 0
              ENDIF
            END DO !level loop
            CALL var_info(var_list(nvars+1), var_description, units, nlevels, nvlp)
            IF (output_fmt == "nc") THEN
              CALL write_data (file_handle, date_str, var_list(nvars+1), & 
              var_description, data_xyz, units, nt / delta_t) 
            ELSEIF (output_fmt == "grib") THEN
              var_name = var_list(nvars+1)(1:LEN_TRIM(var_list(nvars+1))) // '_B'
              CALL writegrib(var_name, mx, my, nlevels, 0, 0, data_xyz, jdate, nt, &
                             tba(nt,t1,delta_t,var_name), nvlp, pres_hpa)
            ENDIF
          ENDIF
          CALL mmdddateadj(date_str, delta_t)  
          IF (output_fmt.eq."grib" .and. multiple_output_files) CALL closegrib()
          PRINT "('*'$)"  ! progress '*'
! time iteration
        END DO 
      ENDIF

      IF (output_fmt == "nc") THEN
        CALL close_cdf(file_handle)
        IF (is == 1) THEN
          CALL write_cdf_cdfapi(output, mx, my, nverlvs+1, nfct + 1, date_str2, delta_t, 0, nvlp, pres_hpa)
        ELSEIF (is == 2) THEN
          CALL write_cdf_cdfapi(output, mx, my, nvlp, nfct + 1, date_str2, delta_t, 1, nvlp, pres_hpa)
        ELSEIF (is == 3) THEN
          CALL write_cdf_cdfapi(output, mx, my, vres, nfct + 1, date_str2, delta_t, 2, nvlp, pres_hpa)
        ENDIF
      ELSEIF (output_fmt == "grib") THEN
        CALL endgrib()
      ENDIF

      PRINT*, '  '

      IF (is == 0) THEN
        DEALLOCATE(vardata, vardata_n)
      ELSE IF (is == 1) THEN
        DEALLOCATE(vardata, data_xyz)
      ELSE
        DEALLOCATE(data_xyz_var, data_xyz_pr, vardata, data_xyz)
      ENDIF

      STOP

    contains

      integer function tba (nt, t1, delta_t, varname)
        implicit none

        integer, intent(in) :: nt, t1, delta_t
        character(len=*), intent(in) :: varname
!JR Changed logic to give correct accumulation behavior for precip. variables.
!JR The equivalent behavior for gribout=.true. is handled in horizontal/output.F90 
!JR        tba=nt-delta_t
!JR        if (tba.lt.0.or.varname.eq.'rn2D'.or.varname.eq.'rc2D'.or.varname.eq.'rg2D') tba=0
        if (nt > delta_t .and. (varname == 'rn2D' .or. varname == 'rc2D' .or. varname == 'rg2D')) then
          tba = nt - delta_t
        else
          tba = 0
        end if
        return
      end function tba
    end program pop

SUBROUTINE read_1var(nt, delta_t, datadir, vardata, var_name, nip, nverlvs,ArchvTimeUnit, nvlp)
      USE fimnc, ONLY: var_info,char_len
      IMPLICIT NONE
      
      INTEGER nt, delta_t, nip, nverlvs, nvlp
      CHARACTER *(*) :: datadir
      CHARACTER *(*) :: var_name
      CHARACTER(80) :: header(10)
!JR Don't know the actual dimension of vardata
      REAL vardata(*)
      character(2),intent(in)::ArchvTimeUnit

      CHARACTER (len=char_len)  :: var_description, units
      INTEGER :: nlevels

      INTEGER time, lunout, is2Dvar, ioerr
      CHARACTER (len=char_len) :: filename

! when all 2d vars, sm3d, and st3d are in one file, comment the following statements
! if these variables are in seperate files.
      IF (is2Dvar(var_name) == 1) THEN
        CALL read_2Dvar(nt, datadir, var_name, nip, nverlvs, vardata,ArchvTimeUnit, nvlp)
!print*, 'var =', var_name, minval(vardata(1:nip)), maxval(vardata(1:nip))
        RETURN
      ENDIF

      lunout = 29
      time = nt
      WRITE(filename,"('fim_out_',a4,i6.6,a2)") var_name,time,ArchvTimeUnit
      filename = datadir(1:LEN_TRIM(datadir)) // filename
      CALL var_info(var_name, var_description, units, nlevels, nvlp)
      OPEN (lunout, file=filename, form="unformatted", action='read', iostat=ioerr)
      if (ioerr /= 0) then
        write(6,*)'read_1var: Failed to open file ', trim(filename), '. Stopping'
        stop
      end if

      READ (lunout, iostat=ioerr) header
      if (ioerr /= 0) then
        write(6,*)'read_1var: Failed to read header from file ',trim(filename), '. Stopping'
        stop
      end if

      READ (lunout, iostat=ioerr) vardata(1:nip*nlevels)
      if (ioerr /= 0) then
        write(6,*)'read_1var: Failed to read vardata from file ',trim(filename), '. ioerr=', ioerr
        write(6,*)'var_name=', var_name, 'nlevels=', nlevels
        write(6,*)'Stopping.'
        stop
      end if

      CLOSE(lunout)
      ! special treament for qv3d and qw3d
      IF (var_name == "qv3D" .OR. var_name == "qw3D") THEN
        vardata(1:nip*nlevels) = vardata(1:nip*nlevels) * 1000.0
      ENDIF
      ! special treament for ph3d 
      IF (var_name == "ph3D") THEN
!JR Changed to multiply by reciprocal to match scalefactor passed in when gribout=true
        vardata(1:nip*nlevels)  = vardata(1:nip*nlevels) * (1./9.8)
      ENDIF
      ! special treament for oz3d
      IF (var_name == "oz3D") THEN
          vardata(1:nip*nlevels)  = vardata(1:nip*nlevels) * 1000.0
      ENDIF
!print*, 'var ', var_name, minval(vardata(1:nip*nlevels)), maxval(vardata(1:nip*nlevels))
END SUBROUTINE read_1var

SUBROUTINE read_2Dvar(nt, datadir, var_name, nip, nverlvs, vardata,ArchvTimeUnit, nvlp)
      USE fimnc, ONLY: var_info, char_len
      IMPLICIT NONE
      
      CHARACTER *(*) :: datadir, var_name
      INTEGER  nt, nip, nverlvs

!JR Don't know the actual dimension of vardata
      REAL vardata(*)

      INTEGER i, j, lunout, nlevels, nvlp
      CHARACTER(len=char_len) :: filename

      CHARACTER(len=80) :: header(10)
      CHARACTER(len=char_len) :: varname, varname_uc
      CHARACTER(len=char_len)  :: var_description, units
      character(2),intent(in)::ArchvTimeUnit
      
      lunout = 29
      varname_uc = var_name(1:LEN_TRIM(var_name))
      CALL tolowercase(var_name)
      WRITE(filename,"('fim_out_2D__',i6.6,a2)") nt,ArchvTimeUnit
      filename = datadir(1:LEN_TRIM(datadir)) // filename
      OPEN(lunout,file=filename,form="unformatted")

      READ(lunout) header
      READ(header,FMT="(4X A4)") varname  ! for now 10:00 am
      CALL tolowercase(varname)

      DO WHILE (var_name /= varname) 
          IF (varname == "sm3d" .OR. varname == "st3d") THEN
            READ(lunout) vardata(1:4*nip)
          ELSE
            READ(lunout) vardata(1:nip)
          ENDIF
          READ(lunout) header
          READ(header,FMT="(4X A4)") varname 
          CALL tolowercase(varname)
      END DO

      CALL var_info(varname_uc, var_description, units, nlevels, nvlp)

      ! special treatment for rainfall amount
      IF (var_name == "rn2d" .OR. var_name == "rc2d" .OR. &
          var_name == "rg2d" ) THEN
        nlevels = 1
      END IF

      READ (lunout) vardata(1:nip*nlevels)
      CLOSE(lunout)

END SUBROUTINE read_2Dvar

SUBROUTINE read_1var_direct(filename, nip, nverlvs, vardata)
      USE fimnc
      IMPLICIT NONE
      
      CHARACTER *(*) :: filename
      INTEGER  nip, nverlvs
      REAL vardata(nip * nverlvs)

      INTEGER :: lunout=30

      OPEN (lunout,file=filename,form="unformatted")
      READ(lunout) vardata(1:nip*nverlvs)
      CLOSE(lunout)
END SUBROUTINE read_1var_direct

INTEGER FUNCTION is2Dvar(var_name)
      IMPLICIT NONE

      CHARACTER *(*) :: var_name

      IF(var_name(3:4) == "2D" .OR. var_name(3:4) == "2d" .or. var_name(1:4) == "iash") THEN
        is2Dvar = 1
      ELSE
        is2Dvar = 0
      ENDIF

      RETURN
END FUNCTION is2Dvar

SUBROUTINE gridid2mxmy(gridid, mx, my)
      IMPLICIT NONE
      INTEGER gridid, mx, my
      IF (gridid == 228) THEN
        mx = 144
        my = 73
      ELSEIF (gridid == 45) THEN
        mx = 288
        my = 145
      ELSEIF (gridid == 3) THEN
        mx = 360
        my = 181
      ELSEIF (gridid == 4) THEN
        mx = 720
        my = 361
      ENDIF
END SUBROUTINE gridid2mxmy

SUBROUTINE mmdddateadj(mmdddate, delta_t)
      IMPLICIT NONE
      CHARACTER*19 mmdddate
      INTEGER delta_t

      INTEGER year, month, day, hour, minute,second,ly 
      INTEGER month_ny(12), month_ly(12)
      CHARACTER dash, col
      DATA month_ny/31,28,31,30,31,30,31,31,30,31,30,31/
      DATA month_ly/31,29,31,30,31,30,31,31,30,31,30,31/
      dash = '-'
      col = ':'

      READ(mmdddate,'(i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,i2)') year,month,day,hour,minute,second
!print*, year,month,day,hour,minute,second
      IF (mod(year, 4) == 0 .AND. mod(year, 100) /= 0) THEN
        ly = 1
      ELSE
        ly = 0 
      ENDIF
      IF (hour + delta_t .GE. 24) THEN
        day = day + 1
        IF ((ly == 0 .AND. day > month_ny(month)) .OR. (ly == 1 .AND. day > month_ly(month))) THEN
          day = 1
          month = month + 1
        ENDIF 
        IF (month > 12) THEN
          month = 1
          year = year + 1
        ENDIF
      ENDIF
      hour = MOD(hour + delta_t, 24)
!print*, year,month,day,hour,minute,second
      WRITE(mmdddate,'(i4.4,A1,i2.2,A1,i2.2,A1,i2.2,A1,i2.2,A1,i2.2)') year,dash,month,dash,day,dash,hour,col,minute,col,second 
     
END SUBROUTINE mmdddateadj

SUBROUTINE jdateadj(jdate, delta_t)
      IMPLICIT NONE
      CHARACTER*7 jdate
      INTEGER delta_t

      INTEGER yr, jday, hr, ly 

      READ(jdate,'(i2,i3,i2)') yr,jday,hr 
      IF (hr + delta_t .GE. 24) THEN
        jday = jday + 1
        IF (mod(yr, 4) == 0 .AND. mod(yr, 100) /= 0) THEN
          ly = 1
        ELSE
          ly = 0 
        ENDIF
        IF (jday > 365 + ly) THEN
          jday = 1
          yr = yr + 1
        ENDIF
      ENDIF
      hr = MOD(hr + delta_t, 24)
      WRITE(jdate,'(i2.2,i3.3,i2.2)') yr,jday,hr
     
END SUBROUTINE jdateadj

