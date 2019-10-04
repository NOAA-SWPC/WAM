program GetQueueTime
use read_queue_namelist,only: GetMaxQueueTime
implicit none
character(8) :: QueueTime
call GetMaxQueueTime (QueueTime)
write(*,"(A)") QueueTime
end program GetQueueTime
