program GetFIMDIR
use read_queue_namelist,only: ReadQUEUEnamelist,FIMDIR
implicit none
call ReadQUEUEnamelist
write(*,"(A)") FIMDIR
end program GetFIMDIR
