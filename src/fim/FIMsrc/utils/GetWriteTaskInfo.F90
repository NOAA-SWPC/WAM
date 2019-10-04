program GetWriteTaskInfo

  use module_wtinfo,only:wtinfo

  implicit none

  integer :: cpn                       ! cores per node (ignored in this pgm)
  integer :: nwt                       ! number of write tasks
  integer :: mwtpn                     ! max write tasks per node
  logical :: root_own_node             ! ignored in this pgm
  logical :: abort_on_bad_task_distrib ! ignored in this pgm
  logical :: debugmsg_on               ! ignored in this pgm

  call wtinfo(cpn,nwt,mwtpn,root_own_node,abort_on_bad_task_distrib,debugmsg_on)

  write (*,'(a,i0)') 'num_write_tasks:',nwt
  write (*,'(a,i0)') 'max_write_tasks_per_node:',mwtpn

end program GetWriteTaskInfo
