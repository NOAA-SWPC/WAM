module module_outqv
contains
!
!$$$   SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .                                       .
! SUBPROGRAM:    OUTQV       PRINT MAX VALUE OF ARRAY
!   PRGMMR:  BENJAMIN, STAN ORG: ERL/PROFS      DATE: 93-01-18
!
! ABSTRACT:  PRINT MAX VALUE AND IPN,IVL OF 2-D ARRAY
!
! PROGRAM HISTORY LOG:
!    88/05/31       S. BENJAMIN     ORIGINAL VERSION
!    2006/02        J. Lee          convert from F77 to F90
!    2006/02        J. Middlecoff   converted to icos notation and parallelized
!
! USAGE:   CALL OUTQV(QVA,lat,lon,NIP,NVL)
!
!   INPUT ARGUMENT LIST:
!     QVA      - REAL     2-D ARRAY
!     lat      - REAL     1-D array
!     lon      - REAL     1-D array
!     NIP      - INTEGER  NO. OF ICOS POINTS
!     NVL      - INTEGER  NO. OF POINTS IN VERTICAL DIRECTION
!
! OUTPUT ARGUMENT LIST: none
!
! REMARKS: NONE
!

SUBROUTINE OUTQV(QVA,lat,lon,NIP,NVL,factor)
implicit none
INTEGER,intent(IN) ::  NIP,NVL
!SMS$DISTRIBUTE (dh,nip) BEGIN
REAL   ,intent(IN) ::  QVA(Nvl,Nip)
real   ,intent(IN) :: lat(nip),lon(nip)
!SMS$DISTRIBUTE END
REAL QVAMAX,qvamin,latIMAX,latIMIN,lonIMAX,lonIMIN,MyQvamax,MyQvamin
REAL factor
INTEGER imax,imin,IPN,IVL,mype,im,MyStartG,MyEndG
logical DiagPrint

  imin = 0
  imax = 0
!SMS$PARALLEL (dh,ipn) BEGIN
print"('    IVL    IPN    LAT     LON    MAX VALUE       IPN    LAT     LON    MIN VALUE')"
DO IVL=1,NVL
  QVAMAX=-1.E30
  QVAMin= 1.E30
  DO IPN=1,NIP
    IF(QVA(ivl,IPN).GE.QVAMAX)THEN
      QVAMAX=QVA(ivl,IPN)
      IMAX=IPN
    ENDIF
    IF(QVA(IVL,ipn).LE.QVAMin)THEN
      QVAMIN=QVA(ivl,IPN)
      IMIN=IPN
    ENDIF
  ENDDO
  MyQvamax = qvamax 
  MyQvamin = qvamin 
!SMS$reduce(qvamax,max)
!SMS$reduce(qvamin,min)
  if(MyQvamax == qvamax) then
    im = imax
    call GetIpnGlobalMype(im,imax,mype,DiagPrint)
  else
    imax    = 0
  endif
  if(MyQvamin == qvamin) then
    im = imin
    call GetIpnGlobalMype(im,imin,mype,DiagPrint)
  else
    imin    = 0
  endif
!SMS$reduce(imax,imin,max)
  ! exit if data is totally hosed
  IF (imin == 0) THEN
    PRINT *,'ERROR OUTQV:  imin==0 at ivl = ',ivl
    STOP
  ENDIF
  IF (imax == 0) THEN
    PRINT *,'ERROR OUTQV:  imax==0 at ivl = ',ivl
    STOP
  ENDIF
  latIMAX = -9999.
  lonIMAX = -9999.
  latIMIN = -9999.
  lonIMIN = -9999.
  DO IPN=1,nip
    MyStartG = ipn
    exit
  enddo
  DO IPN=1,nip
    MyEndG = ipn
  enddo
  if(imin >= MyStartG .and. imin <= MyEndG) then
    latIMIN = lat(imin)
    lonIMIN = lon(imin)
  endif
  if(imax >= MyStartG .and. imax <= MyEndG) then
    latIMAX = lat(imax)
    lonIMAX = lon(imax)
  endif
  if (lonIMAX .gt. 180.) lonimax = lonimax-360.
  if (lonIMin .gt. 180.) lonimin = lonimin-360.
!SMS$reduce(latIMAX,lonIMAX,latIMIN,lonIMIN,max)
  Qvamax = qvamax *factor
  Qvamin = qvamin *factor

  print"(i6,i8,2f8.1,1pe12.3,i10,0p2f8.1,1pe12.3)",IVL,IMAX,latIMAX,lonIMAX,QVAMAX, &
                                                       IMIN,latIMIN,lonIMIN,QVAMIN
ENDDO
!SMS$PARALLEL END
RETURN
end subroutine outqv
end module module_outqv
