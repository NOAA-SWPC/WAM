program GetDATADIR
use read_queue_namelist,only: ReadQUEUEnamelist,DATADIR
implicit none
call ReadQUEUEnamelist
write(*,"(A)") DATADIR
end program GetDATADIR
