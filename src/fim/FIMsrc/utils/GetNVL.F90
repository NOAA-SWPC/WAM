program GetNVL
use read_queue_namelist,only: ReturnNVL
implicit none
integer      :: glvl
call ReturnNVL (glvl)
write(6,"(i0)") glvl
end program GetNVL

