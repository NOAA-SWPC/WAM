      module n_module_gfs_mpi_def
!
      include 'mpif.h'
!
      integer,parameter :: max_inter_groups=100     !max number of quilt server groups
!
      integer stat(MPI_STATUS_SIZE),info
      INTEGER :: MC_COMP, MC_IO, MPI_COMM_ALL, MPI_COMM_ALL_DUP
      INTEGER :: N_GROUP
      logical LIOPE
!
      INTEGER :: num_pes,num_pes_fcst,first_fcst_pe,last_fcst_pe 
      INTEGER :: write_tasks_per_group, write_groups 
      INTEGER :: mpi_comm_inter,mpi_comm_comp,mpi_inter_b
      INTEGER,allocatable :: petlist_fcst(:),petlist_write(:,:)
      INTEGER,dimension(max_inter_groups) :: mpi_comm_inter_array
      logical QUILTING
      character*20 ensmem_name

      end module n_module_gfs_mpi_def
