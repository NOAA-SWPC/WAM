program get_num_cores

  use module_wtinfo,only:wtinfo
  use read_queue_namelist,only:getnprocs

  implicit none

  integer :: cpn                       ! number of cores per node
  integer :: leftover                  ! left over tasks for a mod() calculation
  integer :: mwtpn                     ! max write tasks per node
  integer :: nct                       ! number of compute tasks
  integer :: num_nodes_wt              ! number of nodes filled w/ write tasks
  integer :: numcores_batch            ! numcores_mpirun modified to fill all nodes
  integer :: numcores_donothing        ! number of cores which will be MPI do-nothing tasks
  integer :: numcores_mpirun           ! cores needed by mpirun cmd
  integer :: nwt                       ! number of write tasks
  integer :: tot_nodes                 ! number of nodes needed by batch environment
  logical :: debugmsg_on               ! write-task debug message control

  logical :: root_own_node             ! whether root has a node to himself
  logical :: abort_on_bad_task_distrib ! ignored in this program
  logical :: compute_tasks_after_write_tasks ! whether there are compute tasks after write tasks

  call wtinfo(cpn,nwt,mwtpn,root_own_node,abort_on_bad_task_distrib,debugmsg_on)
  call getnprocs (nct)

  if (mwtpn > cpn) then
    write(6,*) 'get_num_cores: Max write tasks per node=', mwtpn, &
               ' exceeds cores_per_node=', cpn
    stop 999
  end if

!JR TODO: FIX THIS!

  if (nct < cpn .and. .not. root_own_node) then
    write(6,*) 'get_num_cores: root_own_node false and nct < cpn is not allowed.'
    write(6,*) 'This is because core_setup_fim doesnt know the number of compute tasks a-priori.'
    write(6,*) 'It could easily call GetNprocs to get the value, but that routine is INSANELY expensive because'
    write(6,*) 'all MPI tasks open 2 files and read the namelists contained therein in order to get the value.'
    stop 999
  end if

  numcores_donothing = 0

! Initialize numcores_mpirun to number of compute tasks, then add for root being
! on his own node, and write tasks

  compute_tasks_after_write_tasks = .false.
  numcores_mpirun = nct
  if (root_own_node) then
    if (nwt > 0 .or. nct > 1) then ! fill rest of 1st node with idle tasks
      numcores_mpirun = numcores_mpirun + cpn - 1
      numcores_donothing = numcores_donothing + cpn - 1
    end if
    if (nct > 1) then
      compute_tasks_after_write_tasks = .true.
    end if
  else
    if (nct > cpn) then    ! more than root node needed for compute tasks
      compute_tasks_after_write_tasks = .true.
    else if (nwt > 0) then ! fill rest of 1st node with idle tasks
      numcores_mpirun = numcores_mpirun + (cpn - nct)
      numcores_donothing = numcores_donothing + (cpn - nct)
    end if
  end if

! Write tasks: first handle nodes filled with write tasks

  num_nodes_wt = nwt / mwtpn
  numcores_mpirun = numcores_mpirun + (num_nodes_wt * cpn)
  numcores_donothing = numcores_donothing + num_nodes_wt * (cpn - mwtpn)

! Last node with write tasks may have less than mwtpn

  leftover = mod (nwt, mwtpn)
  if (leftover > 0) then
    if (compute_tasks_after_write_tasks) then   
      ! Allocate the full node since compute tasks start at next empty node
      numcores_mpirun = numcores_mpirun + cpn
      numcores_donothing = numcores_donothing + (cpn - leftover)
    else
      ! No compute tasks to follow: Just account for the remaining write tasks
      numcores_mpirun = numcores_mpirun + leftover
      num_nodes_wt = num_nodes_wt + 1
    end if
  end if

! Check to ensure numbers were calculated correctly

  if (nct + nwt + numcores_donothing /= numcores_mpirun) then
    write(*,*) 'get_num_cores failure:', nct, nwt, numcores_donothing, numcores_mpirun
    stop 999
  end if

  write (*,'(a,i0)') 'num_cores_mpirun:', numcores_mpirun

! numcores_batch considers all nodes to be full.
! This is critical for jaguar* because the core count
! requested on the #PBS line must be a multiple of cpn

  numcores_batch = numcores_mpirun
  leftover = mod (numcores_batch, cpn)
  if (leftover /= 0) then
    numcores_batch = numcores_batch + (cpn - leftover)
  end if
  write (*,'(a,i0)') 'num_cores_batch:', numcores_batch

  tot_nodes = numcores_batch / cpn
  write (*,'(a,i0)') 'tot_nodes:', tot_nodes
  write (*,'(a,i0)') 'num_nodes_wt:', num_nodes_wt
  write (*,'(a,i0)') 'num_cores_donothing:', numcores_donothing

! Print number of cores left over due to final allocated node being less than full

  write (*,'(a,i0)') 'num_cores_notattached:', numcores_batch - numcores_mpirun

! root_own_node is needed in some of the scripts: print in a way that grep 
! will be sure to find it

  if (root_own_node) then
    write(*,'(a)')'root_own_node:TRUE'
  else
    write(*,'(a)')'root_own_node:FALSE'
  end if
end program get_num_cores
