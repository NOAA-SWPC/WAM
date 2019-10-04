program GetComputeTasks
use read_queue_namelist,only: GetNprocs
implicit none
integer ::      nprocs
call GetNprocs (nprocs)
write(6,"(i0)") nprocs
end program GetComputeTasks
