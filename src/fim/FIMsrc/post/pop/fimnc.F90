!-------------------------------------------------------------------------
! This file contsins the subroutines to create FIM output
! CF compatible netCDF file.
!
! These subroutines use WRF I/O API to create and write
! the netCDF file. At this point we have to use direct
! netCDF calls in addition to the WRF I/O API defined
! subroutines to make the netCDF file  CF compatible.
!
! Ning Wang   Feb. 2006   Apapted partially from Dan Shaffer's
!                         post.F90;
!
! Ning Wang   May  2006   Added vertically interpolated datasets.
!
!-------------------------------------------------------------------------
!#define VARIABLE_LEVELS
#define ERROR_CHECK(a) call IO_DEBUG(a, __LINE__, __FILE__)
      MODULE fimnc
        IMPLICIT NONE 

        INTEGER, PARAMETER :: dummy_com = 999 ! Dummy MPI communicator
        INTEGER, PARAMETER :: WRF_REAL = 104
        INTEGER, PARAMETER :: char_len = 512

        INTEGER :: nvls, STATUS
        CHARACTER(len=char_len), dimension(3) :: dim_names
        INTEGER, DIMENSION(4) :: domain_start = (/1, 1, 1, 1/)
        INTEGER, DIMENSION(4) :: domain_end
        CHARACTER(len=char_len)  :: sysdepinfo = " "

        INTEGER :: dim_allo_set = 0
        REAL, ALLOCATABLE :: lons(:)
        REAL, ALLOCATABLE :: lats(:)
        REAL, ALLOCATABLE :: levels(:)
        REAL, ALLOCATABLE :: times(:)

      CONTAINS

! subroutine to init and create the netcdf file with specified variables
        SUBROUTINE init_cdf_vars(output_file, file_handle, lon_dim, lat_dim, &
          level_dim, num_fcst_times, training_commit, nvars, var_names, date, nvlp) 

          CHARACTER(len=*) :: output_file
          INTEGER :: file_handle
          INTEGER :: nvlp
          INTEGER :: lon_dim, lat_dim, level_dim, num_fcst_times
          INTEGER :: training_commit, nvars
          CHARACTER(len=*) var_names(nvars)
          CHARACTER(len=19)  :: date

          CHARACTER(len=char_len)  :: var_name, var_description, units

          REAL :: dummy_array(1)
          INTEGER :: i, nlevels
          
          nvls = level_dim - 1

          domain_end(1) = lon_dim
          domain_end(2) = lat_dim
          domain_end(3) = level_dim
          domain_end(4) = num_fcst_times
          
          CALL ext_ncd_ioinit(sysdepinfo, STATUS)

          CALL ext_ncd_open_for_write_begin (output_file , &
            dummy_com, dummy_com, sysdepinfo, file_handle, STATUS )

          CALL ext_ncd_put_dom_ti_char (file_handle, "Title", &
            "FIM_OUTPUT", STATUS)

          ERROR_CHECK(STATUS)

          IF(training_commit == 1) THEN
            DO i = 1, nvars
              CALL var_info(var_names(i), var_description, units, nlevels, nvlp) 
              CALL write_meta(file_handle, date, var_names(i), & 
                var_description, units, nlevels, dummy_array)
            END DO
            CALL ext_ncd_open_for_write_commit (file_handle, STATUS)
            ERROR_CHECK(STATUS)
          ENDIF

        END SUBROUTINE init_cdf_vars

! subroutine to init and create the netcdf file for vertical cross sections 
        SUBROUTINE init_cdf_vars_v(output_file, file_handle, lon_dim, lat_dim, &
          level_dim, num_fcst_times, training_commit, nverlvs, nvars, var_names, date, nvlp) 

          CHARACTER *(*) :: output_file
          INTEGER :: file_handle
          INTEGER :: nvlp
          INTEGER :: lon_dim, lat_dim, level_dim, num_fcst_times
          INTEGER :: training_commit, nverlvs, nvars
          CHARACTER(len=19)  :: date

          CHARACTER(len=char_len)  :: var_name, var_description, units
          CHARACTER(len=*)  :: var_names(nvars)
          CHARACTER(len=char_len)  :: pr3d_name

          INTEGER :: nlevels, i
          REAL :: dummy_array(1)
          
          pr3d_name = "pr3d"
          nvls = nverlvs 
          domain_end(1) = lon_dim
          domain_end(2) = lat_dim
          domain_end(3) = level_dim
          domain_end(4) = num_fcst_times
          
          CALL ext_ncd_ioinit(sysdepinfo, STATUS)

          CALL ext_ncd_open_for_write_begin (output_file , &
            dummy_com, dummy_com, sysdepinfo, file_handle, STATUS )

          CALL ext_ncd_put_dom_ti_char (file_handle, "Title", &
            "FIM_OUTPUT", STATUS)

          ERROR_CHECK(STATUS)

          IF(training_commit == 1) THEN
            DO i = 1, nvars
              CALL var_info(var_names(i), var_description, units, nlevels, nvlp) 
              CALL write_meta(file_handle, date, var_names(i), & 
              var_description, units, level_dim, dummy_array)
            END DO
            CALL var_info(pr3d_name, var_description, units, nlevels, nvlp) 
            CALL write_meta(file_handle, date, pr3d_name, & 
              var_description, units, nlevels, dummy_array)
            CALL ext_ncd_open_for_write_commit (file_handle, STATUS)
            ERROR_CHECK(STATUS)
          ENDIF

        END SUBROUTINE init_cdf_vars_v

! subroutine to init and create the netcdf file with specified variables of
! the original icosahedral grid.
        SUBROUTINE init_cdf_icos(output_file, file_handle, nip, &
          level_dim, num_fcst_times, training_commit, nvars, var_names, date, nvlp) 

          CHARACTER *(*) :: output_file
          INTEGER :: file_handle
          INTEGER :: nvlp
          INTEGER :: nip, level_dim, num_fcst_times
          INTEGER :: training_commit, nvars
          CHARACTER(len=char_len) var_names(nvars)
          CHARACTER(len=19)  :: date

          CHARACTER(len=char_len)  :: var_name, var_description, units

          REAL :: dummy_array(1)
          INTEGER :: i, nlevels
          
          nvls = level_dim - 1

          domain_end(1) = nip
          domain_end(2) = 1
          domain_end(3) = level_dim
          domain_end(4) = num_fcst_times
          
          CALL ext_ncd_ioinit(sysdepinfo, STATUS)

          CALL ext_ncd_open_for_write_begin (output_file , &
            dummy_com, dummy_com, sysdepinfo, file_handle, STATUS )

          CALL ext_ncd_put_dom_ti_char (file_handle, "Title", &
            "FIM_OUTPUT", STATUS)

          ERROR_CHECK(STATUS)

          IF(training_commit == 1) THEN
            DO i = 1, nvars
              CALL var_info(var_names(i), var_description, units, nlevels, nvlp) 
              CALL write_meta(file_handle, date, var_names(i), & 
                var_description, units, nlevels, dummy_array)
            END DO
            CALL ext_ncd_open_for_write_commit (file_handle, STATUS)
            ERROR_CHECK(STATUS)
          ENDIF

        END SUBROUTINE init_cdf_icos

! subroutine to close the netcdf file
        SUBROUTINE close_cdf (file_handle) 
          INTEGER :: file_handle
          call ext_ncd_ioclose(file_handle, STATUS )
          ERROR_CHECK(STATUS)

          call ext_ncd_ioexit(   STATUS )
          ERROR_CHECK(STATUS)

        END SUBROUTINE close_cdf

! subroutine to set dim_names once before first use
        SUBROUTINE set_dim_names
          LOGICAL,SAVE :: first_time = .true.
          IF (first_time) THEN
            dim_names(1) = "lon"
            dim_names(2) = "lat"
            dim_names(3) = "levels"
            first_time = .false.
          ENDIF
        END SUBROUTINE set_dim_names

! subroutines to write a variable to the netcdf file
        SUBROUTINE write_data (file_handle, date, var_name, & 
          var_description, var_data, units, time_step) 
          INTEGER :: file_handle, time_step
          CHARACTER(*)  :: date 
          REAL :: var_data(*)
          CHARACTER(len=*)  :: var_name, var_description, units

          CALL set_dim_names
          domain_start(4) = time_step
          domain_end(4) = time_step 
          CALL ext_ncd_write_field (file_handle, date, var_name, &
            var_data, WRF_REAL, dummy_com, dummy_com,  &
            0, &            ! Domain Descriptor
            "XYZ", &        ! MemoryOrder
            "", &           ! Stagger
            dim_names, &
            domain_start, &        ! Domain Start
            domain_end, &          ! Domain End
            domain_start, &        ! Memory Start
            domain_end, &          ! Memory End
            domain_start, &        ! Patch Start
            domain_end, &          ! Patch End
            STATUS )
          ERROR_CHECK(STATUS)

        END SUBROUTINE write_data


! subroutine to write meta info
        SUBROUTINE write_meta (file_handle, date, var_name, & 
          var_description, units, nlevels, var_data) 
          INTEGER :: file_handle, nlevels
          CHARACTER(*)  :: date 
          REAL :: var_data(*)
          CHARACTER(len=char_len) :: var_name, var_description, units

          CALL set_dim_names
          domain_start(4) = 1
          domain_end(4) = 1 
          CALL ext_ncd_write_field (file_handle, date, var_name, &
            var_data, WRF_REAL, dummy_com, dummy_com,  &
            0, &            ! Domain Descriptor
            "XYZ", &        ! MemoryOrder
            "", &           ! Stagger
            dim_names, &
            domain_start, &        ! Domain Start
            domain_end, &          ! Domain End
            domain_start, &        ! Memory Start
            domain_end, &          ! Memory End
            domain_start, &        ! Patch Start
            domain_end, &          ! Patch End
            STATUS )
          ERROR_CHECK(STATUS)

        ! Add metadata for the field 
          CALL ext_ncd_put_var_ti_char (file_handle ,"description", & 
            var_name, var_description, STATUS)
          ERROR_CHECK(STATUS)

          CALL ext_ncd_put_var_ti_char (file_handle,"units", var_name, &
            units, STATUS )
          ERROR_CHECK(STATUS)

        END SUBROUTINE write_meta

! subroutine needed by wrfio API
        SUBROUTINE io_debug (message_level, line_number, file_name) 
          INTEGER, intent(in) :: message_level, line_number 
          CHARACTER(len=*), intent(in) :: file_name

          IF (message_level > 0) THEN
            PRINT *, "Error at ", file_name, line_number
          END IF
          END SUBROUTINE io_debug

! subroutine to obtain the variable info through variable number
        SUBROUTINE var_info(var_name, var_description, units, nlevels, nvlp)
          CHARACTER(len=*)  :: var_name, var_description, units
          INTEGER :: nlevels, nvlp

          SELECT CASE (var_name)
      
! 3D variables for dynamics 
          CASE ("us3D")
          var_description = "U wind"
          units = "meter.second"
          nlevels = nvls 

          CASE ("vs3D")
          var_description = "V wind"
          units = "meter.second"
          nlevels = nvls

          CASE ("dp3D")
          var_description = "3-D Delta Pressure"
          units = "pascals"
          nlevels = nvls

          CASE ("pr3D")
          var_description = "3-D Pressure"
          units = "pascals"
          nlevels = nvls + 1

          CASE ("mp3D")
          var_description = "Montgomery Potential"
          units = "meter2.second2"
          nlevels = nvls

          CASE ("th3D")
          var_description = "Potential Temperature"
          units = "kelvin"
          nlevels = nvls

          CASE ("ph3D")
          var_description = "Geo Potential"
          units = "meter2/second2"
          nlevels = nvls + 1

          CASE ("qv3D")
          var_description = "Specific humidity"
          units = "non-dimensional"
          nlevels = nvls

          CASE ("rh3D")
          var_description = "Relative humidity"
          units = "percentage"
          nlevels = nvls

          CASE ("vr3D")
          var_description = "Vorticity"
          units = "1/second"
          nlevels = nvls

          CASE ("ws3D")
          var_description = "Omega"
          units = "Micro-bar/second"
          nlevels = nvls

          CASE ("tk3D")
          var_description = "Temperature"
          units = "Kelvin"
          nlevels = nvls

          CASE ("td3D")
          var_description = "Dew-point temperature"
          units = "Kelvin"
          nlevels = nvls

! 3D variables for chemistry
          CASE ("d2st")
          var_description = "coarse dust particles"
          units = "ug/kg"
          nlevels = nvls

          CASE ("d1st")
          var_description = "fine dust particles"
          units = "ug/kg"
          nlevels = nvls

          CASE ("dms1")
          var_description = "dms"
          units = "ppm"
          nlevels = nvls

          CASE ("pso2")
          var_description = "so2"
          units = "ppm"
          nlevels = nvls

          CASE ("sulf")
          var_description = "sulfate"
          units = "ppm"
          nlevels = nvls

          CASE ("pp25")
          var_description = "other primary pm25 +volcanic ash "
          units = "ug/Kg"
          nlevels = nvls

          CASE ("pp10")
          var_description = "other primary pm10 +volcanic ash"
          units = "ug/Kg"
          nlevels = nvls

          CASE ("obc1")
          var_description = "hydrophobic organic carbon"
          units = "ug/Kg"
          nlevels = nvls

          CASE ("obc2")
          var_description = "hydrophillic organic carbon"
          units = "ug/Kg"
          nlevels = nvls

          CASE ("pbc1")
          var_description = "hydrophobic black carbon"
          units = "ug/Kg"
          nlevels = nvls

          CASE ("pbc2")
          var_description = "hydrophillic black carbon"
          units = "ug/Kg"
          nlevels = nvls

          CASE ("ash1")
          var_description = "volcanic ash size bin 1"
          units = "ug/Kg"
          nlevels = nvls

          CASE ("ash2")
          var_description = "volcanic ash size bin 2"
          units = "ug/Kg"
          nlevels = nvls

          CASE ("ash3")
          var_description = "volcanic ash size bin 3"
          units = "ug/Kg"
          nlevels = nvls

          CASE ("ash4")
          var_description = "volcanic ash size bin 4"
          units = "ug/Kg"
          nlevels = nvls

          CASE ("c13D")
          var_description = "radioactive tracer 1, explosive emissions with height"
          units = "?"
          nlevels = nvls

          CASE ("c23D")
          var_description = "radioactive tracer 2, linear emissions with height"
          units = "?"
          nlevels = nvls

          CASE ("oc1P")
          var_description = "hydrophobic organic carbon"
          units = "ug/Kg"
          nlevels = 40

          CASE ("oc2P")
          var_description = "hydrophobic organic carbon"
          units = "ug/Kg"
          nlevels = 40

          CASE ("bc1P")
          var_description = "hydrophobic black carbon"
          units = "ug/Kg"
          nlevels = 40

          CASE ("bc2P")
          var_description = "hydrophobic black carbon"
          units = "ug/Kg"
          nlevels = 40

          CASE ("so2P")
          var_description = "so2"
          units = "ppm"
          nlevels = 40

          CASE ("slfP")
          var_description = "sulfate"
          units = "ppm"
          nlevels = 40

          CASE ("d1sP")
          var_description = "dust particles"
          units = "ppm"
          nlevels = 40

          CASE ("d2sP")
          var_description = "dust particles"
          units = "ppm"
          nlevels = 40

          CASE ("d3sP")
          var_description = "dust particles"
          units = "ppm"
          nlevels = 40

          CASE ("d4sP")
          var_description = "dust particles"
          units = "ppm"
          nlevels = 40

          CASE ("d5sP")
          var_description = "dust particles"
          units = "ppm"
          nlevels = 40

          CASE ("s1sP")
          var_description = "seasalt"
          units = "ppm"
          nlevels = 40

          CASE ("s2sP")
          var_description = "seasalt"
          units = "ppm"
          nlevels = 40

          CASE ("s3sP")
          var_description = "seasalt"
          units = "ppm"
          nlevels = 40

          CASE ("s4sP")
          var_description = "seasalt"
          units = "ppm"
          nlevels = 40

          CASE ("dmsP")
          var_description = "dms"
          units = "ppm"
          nlevels = 40

          CASE ("msaP")
          var_description = "msa"
          units = "ppm"
          nlevels = 40

          CASE ("p25P")
          var_description = "other primary pm25"
          units = "ug/Kg"
          nlevels = 40

          CASE ("p10P")
          var_description = "other primary pm25"
          units = "ug/Kg"
          nlevels = 40

! 3D variables for physics 
          CASE ("qw3D")
          var_description = "liquid cloud mixing ratio"
          units = "Kg/Kg"
          nlevels = nvls

          CASE ("hl3D")
          var_description = "long-wave heating rate"
          units = "Kelvin/second"
          nlevels = nvls

          CASE ("hs3D")
          var_description = "short-wave heating rate"
          units = "Kelvin/second"
          nlevels = nvls

          CASE ("oz3D")
          var_description = "ozone"
          units = "Kg/Kg"
          nlevels = nvls

          CASE ("ar3D")
          var_description = "aerosol"
          units = "Kg/Kg"
          nlevels = nvls

          CASE ("cf3D")
          var_description = "cloud fraction"
          units = "percent(%)"
          nlevels = nvls

          CASE ("st3D")
          var_description = "st"
          units = "kg/meter^2"
          nlevels = 4

          CASE ("sm3D")
          var_description = "sm"
          units = "kg/meter^2"
          nlevels = 4

!2D variables for physics
          CASE ("sn2D")
          var_description = "snow water equivalent"
          units = "meter"
          nlevels = 1

          CASE ("rn2D")
          var_description = "rainfall(accumulated total)"
          units = "millimeter"
          nlevels = 1

          CASE ("rc2D")
          var_description = "rainfall(accumulated conv.)"
          units = "millimeter"
          nlevels = 1

          CASE ("rg2D")
          var_description = "rainfall(accumulated large-scale)"
          units = "millimeter"
          nlevels = 1

          CASE ("r12D")
          var_description = "precipitation (total, since last output)"
          units = "millimeter"
          nlevels = 1

          CASE ("r22D")
          var_description = "precipitation (conv, since last output)"
          units = "millimeter"
          nlevels = 1

          CASE ("r32D")
          var_description = "precipitation (large-scale, since last output)"
          units = "millimeter"
          nlevels = 1

          CASE ("pw2D")
          var_description = "precipitable water"
          units = "millimeter"
          nlevels = 1

          CASE ("ts2D")
          var_description = "skin temperature"
          units = "deg. Kelvin"
          nlevels = 1

          CASE ("us2D")
          var_description = "friction velocity"
          units = "meter/sec."
          nlevels = 1

          CASE ("u12D")
          var_description = "u-component of wind"
          units = "meter/sec."
          nlevels = 1

          CASE ("v12D")
          var_description = "v-component of wind"
          units = "meter/sec."
          nlevels = 1

          CASE ("hf2D")
          var_description = "sensible heat flux"
          units = "watt/meter^2"
          nlevels = 1

          CASE ("qf2D")
          var_description = "water vapor flux"
          units = "kg/meter^2"
          nlevels = 1

          CASE ("sw2D")
          var_description = "sw"
          units = "kg/meter^2"
          nlevels = 1

          CASE ("lw2D")
          var_description = "lw"
          units = "kg/meter^2"
          nlevels = 1

          CASE ("ms2D")
          var_description = "mean sea level pressure"
          units = "Pa"
          nlevels = 1

          CASE ("ct2D")
          var_description = "cloud top height"
          units = "m"
          nlevels = 1

          CASE ("cb2D")
          var_description = "cloud base height"
          units = "m"
          nlevels = 1

          CASE ("io2D")
          var_description = "integrated organic carbon"
          units = "ug/kg"
          nlevels = 1

          CASE ("ib2D")
          var_description = "integrated black carbon"
          units = "ug/kg"
          nlevels = 1

          CASE ("id2D")
          var_description = "integrated fine dust"
          units = "ug/kg"
          nlevels = 1

          CASE ("is2D")
          var_description = "integrated sulfate"
          units = "ppm"
          nlevels = 1

          CASE ("ia2D")
          var_description = "integrated PM25"
          units = "ug/m3"
          nlevels = 1

          CASE ("ao2D")
          var_description = "Aerosol Optical Depth"
          units = "unitless"
          nlevels = 1

          CASE ("iash")
          var_description = "Vertically integrated volcanic ash"
          units = "ug/kg"
          nlevels = 1

          CASE ("fl2D")
          var_description = "fallout"
          units = "?"
          nlevels = 1

          CASE ("rp2D")
          var_description = "relative humidity with respect to precipitable water"
          units = "%"
          nlevels = 1

          CASE ("ol2D")
          var_description = "outgoing LW radiation at top of atmosphere"
          units = "W/m**2"
          nlevels = 1


! FIM output variables at standard pressure levels
          CASE ("hgtP")
          var_description = "height at pressure levels"
          units = "meter2/second2"
          nlevels = nvlp

          CASE ("tmpP")
          var_description = "temperature at pressure levels"
          units = "deg"
          nlevels = nvlp

          CASE ("rp3P")
          var_description = "Relative humidity"
          units = "percentage"
          nlevels = nvlp

          CASE ("up3P")
          var_description = "U wind"
          units = "meter.second"
          nlevels = nvlp 

          CASE ("vp3P")
          var_description = "V wind"
          units = "meter.second"
          nlevels = nvlp

!temp variables 
          CASE ("t1xx")
          var_description = "temporary variable 1"
          units = "unit"
          nlevels = nvls

          CASE ("t2xx")
          var_description = "temporary variable 2"
          units = "unit"
          nlevels = nvls

          CASE ("t3xx")
          var_description = "temporary variable 3"
          units = "unit"
          nlevels = nvls + 1

          CASE ("t4xx")
          var_description = "temporary variable 4"
          units = "unit"
          nlevels = 4 

          CASE ("t5xx")
          var_description = "temporary variable 5"
          units = "unit"
          nlevels = 1 

          CASE default
          write(6,*)'var_info: unknown input variable:', var_name(1:len_trim(var_name))
          STOP

          END SELECT  

        END SUBROUTINE var_info

! subroutine to set number of levels for the model
        SUBROUTINE set_model_nlevels(nlevels)
          INTEGER nlevels
          nvls = nlevels
        END SUBROUTINE set_model_nlevels

! subroutine to specificaly train and commit to write
        SUBROUTINE train_commite (file_handle, var_name, &
          var_description, units, nlevels, date)
          INTEGER :: file_handle, nlevels
          CHARACTER(len=char_len)  :: var_name, var_description, units
          CHARACTER(len=19) date 
          REAL :: dummy_array(1)

            CALL write_meta(file_handle, date, var_name, & 
                var_description, units, nlevels, dummy_array)

            CALL ext_ncd_open_for_write_commit (file_handle, STATUS)
            ERROR_CHECK(STATUS)

        END SUBROUTINE train_commite

! subroutine to write all dimension variables using direct netCDF API
        SUBROUTINE write_cdf_cdfapi(nc_file, lon_dim, lat_dim, & 
          level_dim, num_fcst_times, date, delta_t, vip_flg, nvlp, pres_hpa)
          include 'netcdf.inc'
          CHARACTER *(*) :: nc_file
          CHARACTER(len = 19) date
          INTEGER :: lon_dim, lat_dim, level_dim, num_fcst_times, delta_t, vip_flg, nvlp
          INTEGER :: pres_hpa(*)

          INTEGER :: STATUS, NCID, its, i
          INTEGER :: LAT_VAR_ID, LAT_DIM_ID, LON_VAR_ID, LON_DIM_ID
          INTEGER :: LEVELS_VAR_ID, LEVELS_DIM_ID
          INTEGER :: TIME_VAR_ID, TIME_DIM_ID
          REAL :: D2R 
          CHARACTER(len=50) time_att_str

          D2R = atan(1.0) * 4.0 / 180.0

          IF (dim_allo_set == 0) THEN
            ALLOCATE(lons(lon_dim))
            DO its = 1, lon_dim
              lons(its) = 2.0 * 180.0 * REAL(its - 1.0)  &
                / REAL(lon_dim) 
            END DO

            ALLOCATE(lats(lat_dim))
            DO its = 1, lat_dim
              lats(its) = 180.0 * (REAL(its - 1) - REAL(lat_dim - 1) / 2.0) &
                / REAL(lat_dim - 1.0)  
              lats(its) = -lats(its)
            END DO

            ALLOCATE(levels(level_dim))
            IF (vip_flg == 0) THEN  
              DO its = 1, level_dim
                levels(its) = REAL(its)
              END DO
            ELSE IF (vip_flg == 1) THEN 
              DO its = 1, level_dim
                levels(its) = REAL(pres_hpa(its))
              END DO
            ELSE IF (vip_flg == 2) THEN 
              DO its = 1, level_dim
                levels(its) = REAL(level_dim - its) * 10  ! 10 mb per vertical grid
              END DO
            END IF

            ALLOCATE(times(num_fcst_times))
            DO its = 1, num_fcst_times
              times(its) = REAL(its - 1) * REAL(delta_t)
            END DO
          END IF

          STATUS = nf_open(nc_file, NF_WRITE, NCID)
          ERROR_CHECK(STATUS)

          STATUS = nf_redef(NCID)
          ERROR_CHECK(STATUS)

          STATUS = nf_inq_dimid(NCID, 'lat', LAT_DIM_ID)
          ERROR_CHECK(STATUS)

          STATUS = nf_def_var(NCID, 'lat', nf_float, 1, LAT_DIM_ID, & 
            LAT_VAR_ID)
          ERROR_CHECK(STATUS)

          STATUS = nf_put_att_text(NCID, LAT_VAR_ID, "units", 13, &
            "degrees_north")
          ERROR_CHECK(STATUS)

          STATUS = nf_inq_dimid(NCID, 'lon', LON_DIM_ID)
          ERROR_CHECK(STATUS)

          STATUS = nf_def_var(NCID, 'lon', nf_float, 1, LON_DIM_ID,&
            LON_VAR_ID)
          ERROR_CHECK(STATUS)

          STATUS = nf_put_att_text(NCID, LON_VAR_ID, "units", 12, &
            "degrees_east")
          ERROR_CHECK(STATUS)

          STATUS = nf_inq_dimid(NCID, 'levels', LEVELS_DIM_ID)
          ERROR_CHECK(STATUS)

          STATUS = nf_def_var(NCID, 'levels', nf_float, 1, &
            LEVELS_DIM_ID, LEVELS_VAR_ID)
          ERROR_CHECK(STATUS)

          STATUS = nf_put_att_text(NCID, LEVELS_VAR_ID, "units", 5, &
             "meter")
          ERROR_CHECK(STATUS)

          STATUS = nf_inq_dimid(NCID, 'Time', TIME_DIM_ID)
          ERROR_CHECK(STATUS)

          STATUS = nf_def_var(NCID, 'Time', nf_float, 1, TIME_DIM_ID, &
            TIME_VAR_ID)
          ERROR_CHECK(STATUS)

!          time_att_str = "Fcst intv from " //  date
          time_att_str = "hours since " //  date
          STATUS = nf_put_att_text(NCID, TIME_VAR_ID, "units", 31, &
            time_att_str)
!            "hours since 2006-01-01 00:00:00")
          ERROR_CHECK(STATUS)

          STATUS = nf_enddef(NCID)
          ERROR_CHECK(STATUS)

          STATUS = nf_put_var_real(NCID, LAT_VAR_ID, lats)
          ERROR_CHECK(STATUS)

          STATUS = nf_put_var_real(NCID, LON_VAR_ID, lons)
          ERROR_CHECK(STATUS)

          STATUS = nf_put_var_real(NCID, LEVELS_VAR_ID, LEVELS)
          ERROR_CHECK(STATUS)

          STATUS = nf_put_var_real(NCID, TIME_VAR_ID, TIMES)
          ERROR_CHECK(STATUS)

          STATUS = nf_close(NCID)
          ERROR_CHECK(STATUS)

        END SUBROUTINE write_cdf_cdfapi

      END MODULE fimnc

        SUBROUTINE wrf_debug (message_level, message) 
          INTEGER, INTENT(in) :: message_level 
          CHARACTER(len=*), INTENT(in) :: message

          IF (message_level > 0) THEN
            PRINT *, "Error ", message
          END IF
        END SUBROUTINE wrf_debug
