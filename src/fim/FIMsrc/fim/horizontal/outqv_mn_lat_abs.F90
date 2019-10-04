module module_outqv_mn_lat_abs
contains
!
!$$$   SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .                                       .
! SUBPROGRAM:    OUTQV_MN_lat_abs    PRINT MEAN abs VALUE FOR EACH LAYER FROM ICOS-ARRAY
!                                      for two latitude belts 
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

SUBROUTINE OUTQV_MN_lat_abs (QVA,lat,lon,xlim_lat,NVL,NIP,ITS,factor)
implicit none
INTEGER,intent(IN) ::  NIP,NVL,ITS
!SMS$DISTRIBUTE (dh,nip) BEGIN
REAL   ,intent(IN) ::  QVA(Nvl,Nip)
real   ,intent(IN) :: lat(nip),lon(nip)

!SMS$DISTRIBUTE END
REAL*8 sum, sum1
REAL qvamn, qvamn1, factor, xlim_lat
INTEGER IPN,IVL, isum, isum1

write (6,118)ITS,xlim_lat,factor
118   format ('ITS=',i6,'  Latitude-limit=',f6.1,  &
     '  Scaling-factor=',G10.2 )

!SMS$PARALLEL (dh,ipn) BEGIN
DO IVL=1,NVL
  sum=0.d0
  sum1=0.d0
  isum = 0
  isum1 = 0
  DO IPN=1,NIP
   if (abs(lat(ipn)).gt.xlim_lat) then
    sum = sum + abs(QVA(ivl,ipn) )
    isum = isum + 1 
   else 
    sum1 = sum1 + abs( QVA(ivl,ipn) )
    isum1= isum1+ 1
   end if
  ENDDO
!SMS$reduce(isum,isum1,sum,sum1,SUM)
  qvamn=sum  /float(isum)
  qvamn1=sum1/float(isum1)
!
write (6,120)ITS,IVL,xlim_lat,isum        &
       ,qvamn*factor,isum1,qvamn1*factor
120   format ('ITS=',i6,'  K=',i3,        &
    '  lat GT/LE ',f6.2, 2(' npts=',i8,   &
    '   MeanVal=',f12.4 ) )
ENDDO

!SMS$PARALLEL END
RETURN
end subroutine outqv_mn_lat_abs
end module module_outqv_mn_lat_abs
