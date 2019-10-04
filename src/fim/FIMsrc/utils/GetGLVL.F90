program GetGLVL
use read_queue_namelist,only: ReturnGLVL
implicit none
integer      :: glvl
call ReturnGLVL(glvl)
write(6,"(i0)") glvl
end program GetGLVL

