      module get_variables_for_WAM_IPE_coupling

      use gfs_dyn_resol_def, ONLY: lonf, levs, levh
      use gfs_dyn_layout1,   ONLY: lats_node_a
      use gfs_dyn_machine

      IMPLICIT NONE

      REAL, DIMENSION(:, :, :), ALLOCATABLE :: wwg, zzg, uug, 
     &                                         vvg, ttg, rqg,
     &                                         n2g, ppg
      REAL, DIMENSION(:, :),    ALLOCATABLE :: ps

      END module get_variables_for_WAM_IPE_coupling

