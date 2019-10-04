module module_out2D
contains
subroutine out2D(its,VarName,var,lun2d,time)
use module_control,only: nvl,nvlp1,nip,dt,glvl,curve,yyyymmddhhmm

integer      ,intent(IN) :: its
character*(*),intent(IN) :: VarName
real         ,intent(IN) :: var(size(var))
integer      ,intent(IN) :: lun2d
integer      ,intent(IN) :: time

write(header,100) VarName,yyyymmddhhmm,nvl,glvl,curve,its,time
write(lun2d) header
write(lun2d) var

100 format('FIM ',A,' Forecast initial time YYYYMMDDHHMM: ',A12,/,&
           'Level ',I0,', GLVL= ',I0,', Memory Layout ',I0,', Step ',I0,', ',I0,' hours',/,&
           '3',/,&
           '4',/,&
           '5',/,&
           '6',/,&
           '7',/,&
           '8',/,&
           '9',/,&
           '10')
return
end subroutine out2D
end module module_out2D
