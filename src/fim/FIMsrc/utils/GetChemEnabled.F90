program GetChemEnabled
use read_queue_namelist,only: GetChemOn
implicit none
logical :: chem_enabled
call GetChemOn(chem_enabled)
if (chem_enabled) then
  write(6,*) 'true'
else
  write(6,*) 'false'
endif
end program GetChemEnabled
