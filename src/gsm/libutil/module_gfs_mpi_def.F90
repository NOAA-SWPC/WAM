      module module_gfs_mpi_def
!
	implicit none
!
      include 'mpif.h'
!
      integer,parameter :: max_inter_groups=100     !max number of quilt server groups
!
      integer stat(MPI_STATUS_SIZE),info
      INTEGER :: MC_COMP, MC_IO, MPI_COMM_ALL, MPI_COMM_ALL_DUP
      logical LIOPE
!
      INTEGER :: num_pes,num_pes_fcst,first_fcst_pe,last_fcst_pe 
      INTEGER :: mpi_comm_inter,mpi_comm_comp,mpi_inter_b
      INTEGER,allocatable :: petlist_fcst(:)
      logical QUILTING

      end module module_gfs_mpi_def
