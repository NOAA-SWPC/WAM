program GetParallelism
use read_queue_namelist,only: ReadQUEUEnamelist,ComputeTasks
call ReadQUEUEnamelist
if(ComputeTasks=='S'.or.ComputeTasks=='s') then
  write(6,"(a6)") 'serial'
else
  write(6,"(a8)") 'parallel'
endif
end program GetParallelism
