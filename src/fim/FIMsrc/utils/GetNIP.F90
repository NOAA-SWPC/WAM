program GetNIP
use read_queue_namelist,only: ReturnNIP
implicit none
integer      :: nip
call ReturnNIP(nip)
write(6,"(i0)") nip
end program GetNIP

