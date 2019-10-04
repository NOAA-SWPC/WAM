subroutine GetIpnGlobalMype(ipn,ipnGlobal,mype,DiagPrint)
!For an input ipn this routine returns the global ipn (ipnGLobal), the processor number (mype) and DiagPrint.
!DiagPrint=T means that the input ipn matches PrintIpnDiag from the namelist.

use module_control,only: PrintIpnDiag
!SMS$insert use module_decomp
implicit none
integer,intent(IN) :: ipn
integer,intent(out):: ipnGlobal,mype
logical,intent(out):: DiagPrint

ipnGlobal = ipn
!Not needed for dynamic memory!SMS$insert call nnt_UnsToGlobal(dh,ipn,ipnGlobal)
DiagPrint = ipnGlobal==PrintIpnDiag
mype=0
!SMS$insert call nnt_me(mype)
return
end subroutine GetIpnGlobalMype
