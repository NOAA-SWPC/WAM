      module module_io_mpi_def
!
      include 'mpif.h'
!
      integer,parameter :: max_inter_groups=100     !max number of quilt server groups
!
      INTEGER :: N_GROUP
      INTEGER :: num_pes_fcst,first_fcst_pe,last_fcst_pe 
      INTEGER :: num_pes_wrt,write_tasks_per_group, write_groups 
      INTEGER :: mpi_comm_comp
      INTEGER,allocatable :: petlist_write(:,:)
      INTEGER,dimension(max_inter_groups) :: mpi_comm_inter_array
      logical QUILTING
      character*20 ensmem_name

      end module module_io_mpi_def
