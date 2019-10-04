      module get_variables_for_WAM_IPE_coupling

      USE ESMF
      IMPLICIT NONE

      REAL(ESMF_KIND_R8), DIMENSION(:, :, :), ALLOCATABLE :: wwg, zzg, uug,vvg, ttg, rqg, n2g, ppg
      REAL(ESMF_KIND_R8), DIMENSION(:, :),    ALLOCATABLE :: ps

      END module get_variables_for_WAM_IPE_coupling

