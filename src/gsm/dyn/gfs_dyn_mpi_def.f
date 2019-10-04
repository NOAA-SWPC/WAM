      module gfs_dyn_mpi_def
      use gfs_dyn_machine
      use module_gfs_mpi_def
!    intro of implicit none > >
      implicit none
!	integer :: mpi_real4
!	integer :: mpi_real8
!	^^^^^^^^^^^^^^^^^^^^^

      integer MPI_R_IO, MPI_R_MPI, MPI_R_DEF, MPI_A_DEF
     &,       MPI_R_IO_R,MPI_R_MPI_R
      PARAMETER (MPI_R_IO =MPI_REAL4)
      PARAMETER (MPI_R_IO_R=MPI_REAL8)

ccmr  PARAMETER (MPI_R_MPI=MPI_REAL8)
      PARAMETER (MPI_R_MPI=MPI_REAL4)
      PARAMETER (MPI_R_MPI_R=MPI_REAL8)

      PARAMETER (MPI_R_DEF=MPI_REAL8)
      PARAMETER (MPI_A_DEF=MPI_REAL8)

      integer kind_mpi,kind_sum,kind_mpi_r
ccmr  PARAMETER (kind_mpi=8,kind_sum=8)
      PARAMETER (kind_mpi=4,kind_sum=4,kind_mpi_r=8)
      integer ngrids_sfc,ngrid_global
      parameter(ngrids_sfc=100)

!     REAL(KIND=KIND_io4) ,ALLOCATABLE, target :: buf_sig(:,:)
!    &,                                           buff_grid(:,:)
!    &,                                           buf_sig_n(:,:,:)
!     REAL(KIND=KIND_ior) ,ALLOCATABLE, target :: buf_sig_r(:,:)
!     REAL(KIND=KIND_ior) ,ALLOCATABLE, target :: buf_sig_r(:,:)
!    &,                                           buf_sig_rn(:,:,:)
!    &,                                           buf_grd_r(:,:,:)
      REAL(KIND=KIND_io4) ,POINTER ::  buff_mult(:,:,:)
!     REAL(KIND=KIND_io4) ,POINTER ::  buff_multg(:,:)
      real tmm(10,10)
      end module gfs_dyn_mpi_def
