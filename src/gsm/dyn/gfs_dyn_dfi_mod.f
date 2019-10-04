       module gfs_dyn_dfi_mod
!
!*** jw: set up pointers for digital filter variables
!
        use gfs_dyn_machine, only:kind_io8
        implicit none
!
        type gfs_dfi_grid_gr

          integer z_imp, ps_imp, temp_imp, u_imp, v_imp, tracer_imp,    &
     &            p_imp, dp_imp, dpdt_imp
          real(kind=kind_io8),pointer:: hs(:,:,:)     => null()
          real(kind=kind_io8),pointer:: ps(:,:,:)     => null()
          real(kind=kind_io8),pointer:: t(:,:,:)      => null()
          real(kind=kind_io8),pointer:: u(:,:,:)      => null()
          real(kind=kind_io8),pointer:: v(:,:,:)      => null()
          real(kind=kind_io8),pointer:: tracer(:,:,:) => null()
          real(kind=kind_io8),pointer:: p(:,:,:)      => null()
          real(kind=kind_io8),pointer:: dp(:,:,:)     => null()
          real(kind=kind_io8),pointer:: dpdt(:,:,:)   => null()

        end type gfs_dfi_grid_gr
!
       end module gfs_dyn_dfi_mod
