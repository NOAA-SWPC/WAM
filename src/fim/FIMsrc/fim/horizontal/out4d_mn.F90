module module_out4d_mn
contains
!
!$$$   SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .                                       .
! SUBPROGRAM:    OUTQV_MN    PRINT MEAN VALUE FOR EACH LAYER FROM ICOS-ARRAY
!   PRGMMR:  JIN, ADAPTED FROM OUTQV ORIGINALLY BY S.BENJAMIN  DATE: 07-06-20
!
! ABSTRACT:  PRINT MEAN LAYER VALUE OF 2-D ARRAY
!
! PROGRAM HISTORY LOG:
!    2007/06        J. Lee          adapted codes from outqv.F90
!
! USAGE:   CALL OUTQV_MN(QVA,NVL,NIP,ITS)
!
!   INPUT ARGUMENT LIST:
!     QVA      - REAL     2-D ICOSAHEDRAL TRACER ARRAY
!     NVL      - INTEGER  NO. OF POINTS IN VERTICAL DIRECTION
!     NIP      - INTEGER  NO. OF ICOS POINTS
!     ITS      - INTEGER  NO. OF TIME STEP
!
! OUTPUT ARGUMENT LIST: none
!
! REMARKS: NONE
!

SUBROUTINE OUT4D_MN(QVA,NVL,NIP,NTR,ITS,factor,n)
implicit none
INTEGER,intent(IN) ::  NIP,NVL,NTR,n,ITS
!SMS$DISTRIBUTE (dh,nip) BEGIN
REAL   ,intent(IN) ::  QVA(Nvl,Nip,Ntr)
!SMS$DISTRIBUTE END
REAL*8 sum
REAL qvamn, factor
INTEGER IPN,IVL

!SMS$PARALLEL (dh,ipn) BEGIN
DO IVL=1,NVL
  sum=0.d0
  DO IPN=1,NIP
    sum = sum + QVA(ivl,ipn,n)
  ENDDO
  qvamn=sum/float(nip)
!SMS$reduce(qvamn,SUM)
!
write (6,120)ITS,IVL,qvamn*factor
120   format ('ITS=',i6,'  K=',i3,'   MeanVal=',f12.4 )
ENDDO

!SMS$PARALLEL END
RETURN
end subroutine out4d_mn
end module module_out4d_mn
