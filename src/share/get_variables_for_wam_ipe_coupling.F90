module get_variables_for_WAM_IPE_coupling

  IMPLICIT NONE

  INTEGER, PARAMETER, PRIVATE :: kind_grid = 8

  REAL(kind=kind_grid), DIMENSION(:, :, :), ALLOCATABLE, TARGET :: wwg, zzg, &
                                                                   uug, vvg, &
                                                                   ttg, rqg, &
                                                                   n2g, den, &
                                                                   gmol
  REAL(kind=kind_grid), DIMENSION(:, :, :), ALLOCATABLE :: ppg
  REAL(kind=kind_grid), DIMENSION(:, :),    ALLOCATABLE :: ps

END module get_variables_for_WAM_IPE_coupling

