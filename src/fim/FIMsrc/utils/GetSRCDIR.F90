program GetSRCDIR
use read_queue_namelist,only: ReadQUEUEnamelist,SRCDIR
implicit none
call ReadQUEUEnamelist
write(*,"(A)") SRCDIR
end program GetSRCDIR
