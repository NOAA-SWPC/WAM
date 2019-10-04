program GetWRFmpEnabled
use read_queue_namelist,only: GetWRFmpOn
implicit none
logical :: WRFmp_enabled
call GetWRFmpOn(WRFmp_enabled)
if (WRFmp_enabled) then
  write(6,*) 'true'
else
  write(6,*) 'false'
endif
end program GetWRFmpEnabled
