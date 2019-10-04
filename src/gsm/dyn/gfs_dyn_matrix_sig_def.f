      module gfs_dyn_matrix_sig_def
      use gfs_dyn_machine
      implicit none
      
       real(kind=kind_evod) , allocatable ::
     .  tor_sig(:), d_m(:,:,:),dm205(:,:,:)
      end module gfs_dyn_matrix_sig_def
