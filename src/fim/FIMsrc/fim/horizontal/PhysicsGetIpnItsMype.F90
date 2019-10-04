subroutine PhysicsGetIpnItsMype(ipnGlobal,itsOut,mype,DiagPrint)
!This routine returns the global ipn (ipnGLobal), its (itsOut),the processor number (mype) and DiagPrint.
!DiagPrint=T means that the current ipn in physics.F90 matches PrintIpnDiag from the namelist.
!This routine only works for routines that are in the calling tree whose base is in a 
!parallel DO IPN=1,NIP loop in physics.F90.

use module_control,only: PrintIpnDiag
use module_physics,only: ipn,itsP
!SMS$insert use module_physics,only: dh
implicit none
integer,intent(out):: ipnGlobal,itsOut,mype
logical,intent(out):: DiagPrint

ipnGlobal = ipn
!Not needed for dynamic memory!SMS$insert call nnt_UnsToGlobal(dh,ipn,ipnGlobal)
DiagPrint = ipnGlobal==PrintIpnDiag
itsOut = itsP
mype=0
!SMS$insert call nnt_me(mype)
return
end subroutine PhysicsGetIpnItsMype
