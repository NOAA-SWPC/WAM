    PROGRAM slintest
    USE slint, only:bilinear_init_fn

    IMPLICIT NONE

    INTEGER mx, my, i, nip
    CHARACTER (len=128) grid_file1, grid_file2
    REAL llpoints(10242, 2)
    REAL time_beg, time_end

    nip = 163842
    grid_file1 = "gfsltln_t382.dat"
    grid_file2 = "icos_grid_info_level.dat"
    mx = 1152
    my = 576
    CALL cpu_time(time_beg)
    CALL bilinear_init_fn(grid_file1, mx*my, grid_file2, nip)
    CALL cpu_time(time_end)
    PRINT*, 'It took ', time_end - time_beg, 'seconds to init slint'

!    nip = 10242
!    grid_file = "grid_info.dat"
!    mx = 128
!    my = 64
!    OPEN (10,file=grid_file,status='old',form='unformatted')
!    READ(10)
!    READ(10)
!    READ (10) llpoints(:, 1), llpoints(:, 2)
!    CLOSE(10)

!    print*, llpoints(1:32, 1)
!    print*, llpoints(1:32, 2)
!    DO i = 1, 100
!    print*, i
!      CALL bilinear_interp_i2r() 
!    END DO

    END PROGRAM slintest
