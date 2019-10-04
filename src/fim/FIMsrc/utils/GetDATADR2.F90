program GetDATADR2
use read_queue_namelist,only: ReadQUEUEnamelist,DATADR2
implicit none
call ReadQUEUEnamelist
write(*,"(A)") DATADR2
end program GetDATADR2
