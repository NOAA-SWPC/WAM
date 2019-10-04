program GetWRFcuEnabled
use read_queue_namelist,only: GetWRFcuOn
implicit none
logical :: WRFcu_enabled
call GetWRFcuOn(WRFcu_enabled)
if (WRFcu_enabled) then
  write(6,*) 'true'
else
  write(6,*) 'false'
endif
end program GetWRFcuEnabled
