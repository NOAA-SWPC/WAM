      SUBROUTINE GETOZN(LENT, LM, O3B, K1, K2, FAC, PRSLK, iflip
     &,                 XLAT,ko3)
!
!     This code is written By Shrinivas Moorthi
!
      USE MACHINE , ONLY : kind_phys,kind_rad,kind_io4
      USE FUNCPHYS , ONLY : fpkap
      USE PHYSCONS, ROCP => con_ROCP, PI => con_PI
      implicit none
!     include 'constant.h'
!
      integer jmr, blte, dlte, loz
!  4X5 ozone data
!     parameter (jmr=45, BLTE=-86.0, DLTE=4.0)
! GEOS ozone data
      parameter (jmr=18, BLTE=-85.0, DLTE=10.0, LOZ=17)
!
      real (kind=kind_rad) p00, daylen
      PARAMETER (p00=1000.0)
      logical geosoz
      parameter (geosoz=.true., daylen=86400.0)
!
      integer lent, lm, k1, k2, iflip, ko3
!
!     LOCALS
!
      real (kind=kind_rad) O3r(JMr,LOZ,12), O3B(lent,LM)
     &,                    PRSLK(LENT,LM),  O3I(lent,LOZ), PKSTR(LOZ)
     &,                    PSTR(LOZ),       xlat(lent)
     &,                    wk1(lent)
!
      integer imond(12), ilat(jmr,12)
      real (kind=kind_io4) pstr4(loz), o3clim4(jmr,loz,12)
!
      LOGICAL     FIRST
      DATA  FIRST/.TRUE./
      DATA PKSTR/LOZ*0.0/, PSTR/LOZ*0.0/
!     DATA O3Z/JMLZ13*0.0/
!
      real (kind=kind_rad) tem, tem1, tem2, tem3, tem4
     &,                    temp, rdg, fac, deglat, elte
!
      SAVE first, pkstr, pstr, o3r, elte
!
      integer i, j, k, l, nm, j1, j2, ll
!
!     if(me .eq. 0) WRITE (6,989) jdat(1),jdat(2),jdat(3),jdat(5)
! 989 FORMAT(' UPDATING OZONE FOR ', I4,I3,I3,I3)
!
      IF (FIRST) THEN
         REWIND KO3
         elte = blte + (jmr-1) * dlte
         do l=1,loz
           READ (KO3,15) pstr4(l)
         enddo
         pstr = pstr4
         DO nm=1,12
           do j=1,jmr
             READ (KO3,19) imond(nm),ilat(j,nm),
     &                     (o3clim4(j,l,nm), l=1,10)
             READ (KO3,20) (o3clim4(j,l,nm), l=11,loz)
           ENDDO
         ENDDO
         O3R = o3clim4
         do  nm=1,12
           do l=1,loz
             do j=1,jmr
               o3r(j,l,nm) = o3r(j,l,nm) * 1.655e-6
             enddo
           enddo
         enddo
!
   15    format(f10.3)
   19    format(i2,i4,10f6.2)
   20    format(6x,10f6.2)
!
         PRINT *,' FOUND OZONE DATA FOR LEVELS PSTR=',(PSTR(L),L=1,LOZ)
!        print *,' O3=',(o3r(15,l,1),l=1,loz)
!
         DO L=1,LOZ
           PKSTR(L) = fpkap(PSTR(L)*100.0)
!          PKSTR(L) = (PSTR(L)/P00) ** ROCP
         ENDDO
!
         FIRST  = .FALSE.
      ENDIF
!
      RDG = 180.0 / PI
      DO I=1,LENT
        deglat = xlat(i)*rdg
        if (deglat .gt. blte .and. deglat .lt. elte) then
          tem1 = (deglat - BLTE)/DLTE + 1
          J1   = tem1
          J2   = J1 + 1
          tem1 = tem1 - J1
        elseif (deglat .le. blte) then
          j1 = 1
          j2 = 1
          tem1 = 1.0
        elseif (deglat .ge. elte) then
          j1 = jmr
          j2 = jmr
          tem1 = 1.0
        endif
        tem2 = 1.0 - tem1
        DO J=1,LOZ
          tem3     = tem2*o3r(j1,J,k1) + tem1*o3r(j2,J,k1)
          tem4     = tem2*o3r(j1,J,k2) + tem1*o3r(j2,J,k2)
          O3I(I,J) = tem4*fac          + tem3*(1.0-fac)
        ENDDO
      ENDDO
!     DO I=1,LENT
!       PIK(I) = (PS(I)*10.0/P00) ** ROCP
!     ENDDO
!
      DO L=1,LM
        LL = L
        if (iflip .eq. 1) LL = LM + 1 -L
        DO I=1,LENT
          WK1(I) = PRSLK(I,LL)
        ENDDO
        DO K=1,LOZ-1
          temp = 1.0 / (PKSTR(K+1) - PKSTR(K))
          DO I=1,LENT
            IF (WK1(I).GT.PKSTR(K) .AND. WK1(I).LE.PKSTR(K+1)) THEN
              TEM      = (PKSTR(K+1) - WK1(I)) * TEMP
              O3B(I,L) = TEM * O3I(I,K) + (1.0 - TEM) * O3I(I,K+1)
            ENDIF
          ENDDO
        ENDDO
        DO I=1,LENT
          IF (WK1(I) .GT. PKSTR(LOZ)) O3B(I,L) = O3I(I,LOZ)
          IF (WK1(I) .LT. PKSTR(1))   O3B(I,L) = O3I(I,1)
        ENDDO
      ENDDO
!
      RETURN
      END
