program GetPREPDIR
use read_queue_namelist,only: ReadQUEUEnamelist,PREPDIR
implicit none
call ReadQUEUEnamelist
write(*,"(A)") PREPDIR
end program GetPREPDIR
