      SUBROUTINE CLO89(CLDFAC,CAMT,NCLDS,KBTM,KTOP
     &,                L, LP1, IMAX, NVECT)
!
!
!     SUBROUTINE CLO88 COMPUTES CLOUD TRANSMISSION FUNCTIONS FOR THE
!  LONGWAVE CODE,USING CODE WRITTEN BY BERT KATZ (301-763-8161).
!  AND MODIFIED BY DAN SCHWARZKOPF IN DECEMBER,1988.
!                INPUTS:          (MODULE BLOCK)
!      CAMT,KTOP,KBTM,NCLDS         RADISW
!                OUTPUT:
!      CLDFAC                       CLDCOM
!
!          CALLED BY:      RADMN OR MODEL ROUTINE
!          CALLS    :
!
!
      USE MACHINE , ONLY : kind_rad
      implicit none
!
      integer IMAX,KP,K,LP1,L,NVECT
      integer NCLDS(IMAX),KTOP(IMAX,LP1),KBTM(IMAX,LP1)
      real (kind=kind_rad) CAMT(IMAX,LP1),CLDFAC(IMAX,LP1,LP1)
     &,                     CLDROW(LP1)
!
       real (kind=kind_rad) CLDIPT(LP1,LP1,NVECT)
!
      real (kind=kind_rad) xcld
      integer iq, itop, jtop, ip, ir, i, j, k1, k2, kt, nc, kb
!
      DO IQ=1,IMAX,NVECT
!
         ITOP = IQ + (NVECT-1)
         IF (ITOP .GT. IMAX) ITOP = IMAX
         JTOP = ITOP - IQ + 1
!
         DO IP=1,JTOP
           IR = IQ + IP - 1
           IF (NCLDS(IR).EQ.0) THEN
             DO J=1,LP1
               DO I=1,LP1
                 CLDIPT(I,J,IP) = 1.
               ENDDO
             ENDDO
           ENDIF
           IF (NCLDS(IR).GE.1) THEN
             XCLD = 1.-CAMT(IR,2)
             K1   = KTOP(IR,2) + 1
             K2   = KBTM(IR,2)
             DO J=1,LP1
               CLDROW(J) = 1.
             ENDDO
             DO J=1,K2
               CLDROW(J) = XCLD
             ENDDO
             KB = MAX(K1,K2+1)
             DO K=KB,LP1
               DO KP=1,LP1
                 CLDIPT(KP,K,IP) = CLDROW(KP)
               ENDDO
             ENDDO
             DO J=1,LP1
               CLDROW(J) = 1.
             ENDDO
             DO J=K1,LP1
               CLDROW(J) = XCLD
             ENDDO
             KT = MIN(K1-1,K2)
             DO K=1,KT
               DO KP=1,LP1
                 CLDIPT(KP,K,IP) = CLDROW(KP)
               ENDDO
             ENDDO
             IF (K2+1 .LE. K1-1) THEN
               DO J=K2+1,K1-1
                 DO I=1,LP1
                   CLDIPT(I,J,IP) = 1.
                 ENDDO
               ENDDO
             ELSE IF(K1.LE.K2) THEN
               DO J=K1,K2
                 DO I=1,LP1
                   CLDIPT(I,J,IP) = XCLD
                 ENDDO
               ENDDO
             ENDIF
           ENDIF
           IF (NCLDS(IR).GE.2) THEN
             DO NC=2,NCLDS(IR)
               XCLD = 1. - CAMT(IR,NC+1)
               K1   = KTOP(IR,NC+1)+1
               K2   = KBTM(IR,NC+1)
               DO J=1,LP1
                 CLDROW(J) = 1.
               ENDDO
               DO J=1,K2
                 CLDROW(J) = XCLD
               ENDDO
               KB = MAX(K1,K2+1)
               DO K=KB,LP1
                 DO KP=1,LP1
                    CLDIPT(KP,K,IP) = CLDIPT(KP,K,IP)*CLDROW(KP)
!                   CLDFIP(KP,K)    = CLDROW(KP)
                 ENDDO
               ENDDO
               DO J=1,LP1
                 CLDROW(J) = 1.
               ENDDO
               DO J=K1,LP1
                 CLDROW(J) = XCLD
               ENDDO
               KT = MIN(K1-1,K2)
               DO K=1,KT
                 DO KP=1,LP1
                   CLDIPT(KP,K,IP) = CLDIPT(KP,K,IP)*CLDROW(KP)
!                  CLDFIP(KP,K)    = CLDROW(KP)
                 ENDDO
               ENDDO
!              IF(K2+1.LE.K1-1) THEN
!                DO J=K2+1,K1-1
!                  DO I=1,LP1
!                    CLDIPT(I,J,IP) = 1.
!                  ENDDO
!                ENDDO
               IF (K1 .LE. K2) THEN
                 DO J=K1,K2
                   DO I=1,LP1
                     CLDIPT(I,J,IP) = CLDIPT(I,J,IP)*XCLD
                   ENDDO
                 ENDDO
               ENDIF
!              DO J=1,LP1
!                DO I=1,LP1
!                  CLDIPT(I,J,IP) = CLDIPT(I,J,IP)*CLDFIP(I,J)
!                ENDDO
!              ENDDO
             ENDDO
           ENDIF
!
         ENDDO
!
         DO J=1,LP1
           DO I=1,LP1
             DO IP=1,JTOP
               IR             = IQ + IP - 1
               CLDFAC(IR,I,J) = CLDIPT(I,J,IP)
             ENDDO
           ENDDO
         ENDDO
!
      ENDDO
!
      RETURN
      END
!
      SUBROUTINE E1E290(G1,G2,G3,G4,G5,EMISS,FXOE1,DTE1,FXOE2,DTE2,
     &                  AVEPHI, L, LP1, IMAX, EM1V, EM1VW, T1, T2, T4)
CFPP$ NOCONCUR R
!
!     SUBROUTINE E1E290 COMPUTES THE EXCHANGE TERMS IN THE FLUX EQUATION
!  FOR LONGWAVE RADIATION FOR ALL TERMS EXCEPT THE EXCHANGE WITH THE
!  TOP OF THE ATMOSPHERE. THE METHOD IS A TABLE LOOKUP ON A PRE-
!  COMPUTED E2 FUNCTION (DEFINED IN REF. (4)).
!      THE E1 FUNCTION  CALCULATIONS (FORMERLY DONE IN SUBROUTINE
!  E1V88 COMPUTE THE FLUX RESULTING FROM THE EXCHANGE OF PHOTONS
!  BETWEEN A LAYER AND THE TOP OF THE ATMOSPHERE.  THE METHOD IS A
!  TABLE LOOKUP ON A PRE-COMPUTED E1 FUNCTION.
!     CALCULATIONS ARE DONE IN TWO FREQUENCY RANGES:
!       1) 0-560,1200-2200 CM-1   FOR Q(APPROX)
!       2) 160-560 CM-1           FOR Q(APPROX,CTS).
!  MOTIVATION FOR THESE CALCULATIONS IS IN REFERENCES (1) AND (4).
!       INPUTS:                    (MODULE BLOCKS)
!     EM1V,EM1VW,T1,T2,T4              TABCOM
!     AVEPHI                           TFCOM
!     TEMP                             RADISW
!     T                                KDACOM
!     FXOE1,DTE1                ARGUMENT LIST
!     FXOE2,DTE2                ARGUMENT LIST
!       OUTPUTS:
!     EMISS                            TFCOM
!     G1,G2,G3                  ARGUMENT LIST,FOR 1ST FREQ. RANGE
!     G4,G5                     ARGUMENT LIST,FOR 2ND FREQ. RANGE
!
!        CALLED BY :     FST88
!        CALLS     :
!
!
c$$$      USE MACHINE , ONLY : kind_rad
      USE HCON
      implicit none
!
      integer  L, LP1, IMAX
      real (kind=kind_rad) AVEPHI(IMAX,LP1), EMISS(IMAX,LP1)
!
      integer   IT1(IMAX,3*L+2)
      real (kind=kind_rad) FYO(IMAX,LP1), DU(IMAX,LP1),
     &                      WW1(IMAX,LP1), WW2(IMAX,LP1)
!---VARIABLES IN THE ARGUMENT LIST
      real (kind=kind_rad) T1(5040),        T2(5040), T4(5040)
     &,                     EM1V(5040),      EM1VW(5040)
     &,                     FXOE1(IMAX,LP1), DTE1(IMAX,LP1)
     &,                     FXOE2(IMAX,LP1), DTE2(IMAX,LP1),
     &                      G1(IMAX,LP1),    G2(IMAX,L)
     &,                     G3(IMAX,LP1),    G4(IMAX,LP1), G5(IMAX,L)
!
      real (kind=kind_rad) TMP3, tem1, tem2, tem3, tem4
      integer LL, LLP1, KP, I, IVAL, K1, K2, item
!
!---FIRST WE OBTAIN THE EMISSIVITIES AS A FUNCTION OF TEMPERATURE
!   (INDEX FXO) AND WATER AMOUNT (INDEX FYO). THIS PART OF THE CODE
!   THUS GENERATES THE E2 FUNCTION. THE FXO INDICES HAVE BEEN
!   OBTAINED IN FST88, FOR CONVENIENCE.
!
!---THIS SUBROUTINE EVALUATES THE K=1 CASE ONLY--
!
!---THIS LOOP REPLACES LOOPS GOING FROMI=1,IMAX AND KP=2,LP1 PLUS
!   THE SPECIAL CASE FOR THE LP1TH LAYER.
!     LP2  = L + 2
      LL   = L + L
      LLP1 = LL + 1
!
      DO KP=1,LP1
        DO I=1,IMAX
          TMP3        = LOG10(AVEPHI(I,KP)) + H16E1
          FYO(I,KP)   = AINT(TMP3*TEN)
          DU(I,KP)    = TMP3 - HP1*FYO(I,KP)
          FYO(I,KP)   = H28E1 * FYO(I,KP)
          IVAL        = FYO(I,KP) + FXOE2(I,KP)
          EMISS(I,KP) = T1(IVAL ) + DU(I,KP)   * T2(IVAL)
     &                            + DTE2(I,KP) * T4(IVAL)
        ENDDO
      ENDDO
!
!---THE SPECIAL CASE EMISS(I,L) (LAYER KP) IS OBTAINED NOW
!   BY AVERAGING THE VALUES FOR L AND LP1:
      DO I=1,IMAX
        EMISS(I,L) = HAF*(EMISS(I,L) + EMISS(I,LP1))
      ENDDO
!
!   CALCULATIONS FOR THE KP=1 LAYER ARE NOT PERFORMED, AS
!   THE RADIATION CODE ASSUMES THAT THE TOP FLUX LAYER (ABOVE THE
!   TOP DATA LEVEL) IS ISOTHERMAL, AND HENCE CONTRIBUTES NOTHING
!   TO THE FLUXES AT OTHER LEVELS.
!
!***THE FOLLOWING IS THE CALCULATION FOR THE E1 FUNCTION, FORMERLY
!    DONE IN SUBROUTINE E1V88. THE MOVE TO E1E288 IS DUE TO THE
!    SAVINGS IN OBTAINING INDEX VALUES (THE TEMP. INDICES HAVE
!    BEEN OBTAINED IN FST88, WHILE THE U-INDICES ARE OBTAINED
!    IN THE E2 CALCS.,WITH K=1).
!
!   FOR TERMS INVOLVING TOP LAYER, DU IS NOT KNOWN; IN FACT, WE
!   USE INDEX 2 TO REPERSENT INDEX 1 IN PREV. CODE. THIS MEANS THAT
!    THE IT1 INDEX 1 AND LLP1 HAS TO BE CALCULATED SEPARATELY. THE
!   INDEX LLP2 GIVES THE SAME VALUE AS 1; IT CAN BE OMITTED.
!
      DO I=1,IMAX
        IT1(I,1) = FXOE1(I,1)
        WW1(I,1) = TEN - DTE1(I,1)
        WW2(I,1) = HP1
      ENDDO
CDIR$ IVDEP
      DO KP=1,L
        K1 = KP + 1
        K2 = KP + LP1
        DO I=1,IMAX
          IT1(I,K1) = FYO(I,KP) + FXOE1(I,K1)
          IT1(I,K2) = FYO(I,KP) + FXOE1(I,KP)
          WW1(I,K1) = TEN       - DTE1(I,K1)
          WW2(I,K1) = HP1       - DU(I,KP)
        ENDDO
      ENDDO
      DO KP=1,L
        DO I=1,IMAX
          IT1(I,KP+LLP1) = FYO(I,KP) + FXOE1(I,1)
        ENDDO
      ENDDO
!
!  G3(I,1) HAS THE SAME VALUES AS G1 (AND DID ALL ALONG)
      DO I=1,IMAX
        TEM1 = WW1(I,1) * WW2(I,1)
        TEM2 = WW2(I,1) * DTE1(I,1)
        ITEM = IT1(I,1)
        G1(I,1) = TEM1 * EM1V(ITEM)  + TEM2 * EM1V(ITEM+1)
        G4(I,1) = TEM1 * EM1VW(ITEM) + TEM2 * EM1VW(ITEM+1)
        G3(I,1) = G1(I,1)
      ENDDO
      DO KP=1,L
        K1 = KP + 1
        DO I=1,IMAX
          TEM1 = WW1(I,K1)  * WW2(I,K1)
          TEM2 = WW2(I,K1)  * DTE1(I,K1)
          TEM3 = WW1(I,K1)  * DU(I,KP)
          TEM4 = DTE1(I,K1) * DU(I,KP)
          ITEM = IT1(I,K1)
!
          G1(I,K1) = TEM1 * EM1V(ITEM)     + TEM2 * EM1V(ITEM+1)+
     &               TEM3 * EM1V(ITEM+28)  + TEM4 * EM1V(ITEM+29)
          G4(I,K1) = TEM1 * EM1VW(ITEM)    + TEM2 * EM1VW(ITEM+1)+
     &               TEM3 * EM1VW(ITEM+28) + TEM4 * EM1VW(ITEM+29)
!
          TEM1 = WW1(I,KP)  * WW2(I,K1)
          TEM2 = WW2(I,K1)  * DTE1(I,KP)
          TEM3 = WW1(I,KP)  * DU(I,KP)
          TEM4 = DTE1(I,KP) * DU(I,KP)
          ITEM = IT1(I,LP1+KP)
!
          G2(I,KP) = TEM1 * EM1V(ITEM)     + TEM2 * EM1V(ITEM+1)+
     &               TEM3 * EM1V(ITEM+28)  + TEM4 * EM1V(ITEM+29)
          G5(I,KP) = TEM1 * EM1VW(ITEM)    + TEM2 * EM1VW(ITEM+1)+
     &               TEM3 * EM1VW(ITEM+28) + TEM4 * EM1VW(ITEM+29)
        ENDDO
      ENDDO
      DO KP=2,LP1
        DO I=1,IMAX
          ITEM = IT1(I,LL+KP)
          G3(I,KP) = WW1(I,1)  * WW2(I,KP)  * EM1V(ITEM)+
     &               WW2(I,KP) * DTE1(I,1)  * EM1V(ITEM+1)+
     &               WW1(I,1)  * DU(I,KP-1) * EM1V(ITEM+28)+
     &               DTE1(I,1) * DU(I,KP-1) * EM1V(ITEM+29)
        ENDDO
      ENDDO
!
      RETURN
      END
      SUBROUTINE E290(EMISSB,EMISS,AVEPHI,KLEN,FXOE2,DTE2
     &,               L, LP1, IMAX, T1, T2, T4)
CFPP$ NOCONCUR R
!
!     SUBROUTINE E290 COMPUTES THE EXCHANGE TERMS IN THE FLUX EQUATION
!  FOR LONGWAVE RADIATION FOR ALL TERMS EXCEPT THE EXCHANGE WITH THE
!  TOP OF THE ATMOSPHERE. THE METHOD IS A TABLE LOOKUP ON A PRE-
!  COMPUTED E2 FUNCTION (DEFINED IN REF. (4)).
!     CALCULATIONS ARE DONE IN THE FREQUENCY RANGE:
!       1) 0-560,1200-2200 CM-1   FOR Q(APPROX)
!  MOTIVATION FOR THESE CALCULATIONS IS IN REFERENCES (1) AND (4).
!       INPUTS:                    (MODULE BLOCKS)
!     T1,T2,T4,                  TABCOM
!     AVEPHI                           TFCOM
!     FXOE2,DTE2,KLEN           ARGUMENT LIST
!       OUTPUTS:
!     EMISS,EMISSB                     TFCOM
!
!        CALLED BY :     FST88
!        CALLS     :
!
!
c$$$      USE MACHINE , ONLY : kind_rad
      USE HCON
      implicit none
!
      integer L, LP1, IMAX, KLEN
      real (kind=kind_rad) EMISSB(IMAX,LP1), EMISS(IMAX,LP1)
     &,                     AVEPHI(IMAX,LP1)
     &,                     DT(IMAX,LP1),     FYO(IMAX,LP1)
     &,                     DU(IMAX,LP1)
      integer               IVAL(IMAX,LP1)
!---VARIABLES IN THE ARGUMENT LIST
      real (kind=kind_rad) T1(5040),        T2(5040),   T4(5040)
     &,                     FXOE2(IMAX,LP1), DTE2(IMAX,LP1)
!
      real (kind=kind_rad) tmp3
      integer lp2,i, k, k1, item
!
!---FIRST WE OBTAIN THE EMISSIVITIES AS A FUNCTION OF TEMPERATURE
!   (INDEX FXO) AND WATER AMOUNT (INDEX FYO). THIS PART OF THE CODE
!   THUS GENERATES THE E2 FUNCTION.
!
!---CALCULATIONS FOR VARYING KP (FROM KP=K+1 TO LP1, INCLUDING SPECIAL
!   CASE: RESULTS ARE IN EMISS
 
      LP2 = L + 2
      DO K=1,LP2-KLEN
        K1 = K + KLEN - 1
        DO I=1,IMAX
          TMP3        = LOG10(AVEPHI(I,K1)) + H16E1
          FYO(I,K)    = AINT(TMP3*TEN)
          DU(I,K)     = TMP3  - HP1*FYO(I,K)
          FYO(I,K)    = H28E1 * FYO(I,K)
          ITEM        = FYO(I,K) + FXOE2(I,K1)
          EMISS(I,K1) = T1(ITEM) + DU(I,K)    * T2(ITEM)
     &                           + DTE2(I,K1) * T4(ITEM)
        ENDDO
      ENDDO
 
!---THE SPECIAL CASE EMISS(I,L) (LAYER KP) IS OBTAINED NOW
!   BY AVERAGING THE VALUES FOR L AND LP1:
 
      DO I=1,IMAX
        EMISS(I,L) = HAF*(EMISS(I,L) + EMISS(I,LP1))
      ENDDO
 
!---NOTE THAT EMISS(I,LP1) IS NOT USEFUL AFTER THIS POINT.
!
!---CALCULATIONS FOR KP=KLEN AND VARYING K; RESULTS ARE IN EMISSB.
!  IN THIS CASE, THE TEMPERATURE INDEX IS UNCHANGED, ALWAYS BEING
!  FXO(I,KLEN-1); THE WATER INDEX CHANGES, BUT IS SYMMETRICAL WITH
!  THAT FOR THE VARYING KP CASE.NOTE THAT THE SPECIAL CASE IS NOT
!  INVOLVED HERE.
!     (FIXED LEVEL) K VARIES FROM (KLEN+1) TO LP1; RESULTS ARE IN
!   EMISSB(I,(KLEN) TO L)
!
      DO K=1,LP1-KLEN
        DO I=1,IMAX
          DT(I,K)   = DTE2(I,KLEN-1)
          IVAL(I,K) = FYO(I,K) + FXOE2(I,KLEN-1)
        ENDDO
      ENDDO
!
      DO K=1,LP1-KLEN
        K1 = K + KLEN - 1
        DO I=1,IMAX
          ITEM = IVAL(I,K)
          EMISSB(I,K1) = T1(ITEM) + DU(I,K) * T2(ITEM)
     &                            + DT(I,K) * T4(ITEM)
        ENDDO
      ENDDO
!
      RETURN
      END
      SUBROUTINE E2SPEC(EMISS,AVEPHI,FXOSP,DTSP
     &,                 LP1, IMAX, T1, T2, T4)
CFPP$ NOCONCUR R
!
!     SUBROUTINE E2SPEC COMPUTES THE EXCHANGE TERMS IN THE FLUX EQUATION
!  FOR LONGWAVE RADIATION FOR 2 TERMS USED FOR NEARBY LAYER COMPU-
!  TATIONS. THE METHOD IS A TABLE LOOKUP ON A PRE-
!  COMPUTED E2 FUNCTION (DEFINED IN REF. (4)).
!     CALCULATIONS ARE DONE IN THE FREQUENCY RANGE:
!        0-560,1200-2200 CM-1
!  MOTIVATION FOR THESE CALCULATIONS IS IN REFERENCES (1) AND (4).
!       INPUTS:                    (MODULE BLOCKS)
!     T1,T2,T4,                  TABCOM
!     AVEPHI                           TFCOM
!     FXOSP,DTSP                ARGUMENT LIST
!       OUTPUTS:
!     EMISS                            TFCOM
!
!        CALLED BY :     FST88
!        CALLS     :
!
!
c$$$      USE MACHINE , ONLY : kind_rad
      USE HCON
      implicit none
!
      integer LP1, IMAX
      real (kind=kind_rad) AVEPHI(IMAX,LP1), EMISS(IMAX,LP1)
!---VARIABLES IN THE ARGUMENT LIST
      real (kind=kind_rad) T1(5040),      T2(5040), T4(5040)
     &,                     FXOSP(IMAX,2), DTSP(IMAX,2)
      real (kind=kind_rad) TMP3, FYO, DU
      integer k, i, ival
!
!---FIRST WE OBTAIN THE EMISSIVITIES AS A FUNCTION OF TEMPERATURE
!   (INDEX FXO) AND WATER AMOUNT (INDEX FYO). THIS PART OF THE CODE
!   THUS GENERATES THE E2 FUNCTION.
!
      DO  K=1,2
        DO  I=1,IMAX
          TMP3       = LOG10(AVEPHI(I,K)) + H16E1
          FYO        = AINT(TMP3*TEN)
          DU         = TMP3 - HP1*FYO
          IVAL       = H28E1*FYO + FXOSP(I,K)
          EMISS(I,K) = T1(IVAL) + DU*T2(IVAL) + DTSP(I,K)*T4(IVAL)
      ENDDO
      ENDDO
!
      RETURN
      END
!     SUBROUTINE E3V88 COMPUTES NEARBY LAYER TRANSMISSIVITIES FOR
!  H2O USING A TABLE LOOKUP OF THE PRE-COMPUTED E3 FUNCTION
! ( DESCRIBED IN REF. (4)).
!         INPUTS:                 (MODULE BLOCKS,ARGS.)
!       TV,AV                      ARGUMENT LIST
!       EM3                        TABCOM
!          OUTPUTS:
!       EMV                        ARGUMENT LIST
!
!       CALLED BY  :    FST88
!       CALLS      :    ALOG10H,ALOG10V
!
      SUBROUTINE E3V88(EMV,TV,AV, LLP1, IMAX, EM3V)
CFPP$ NOCONCUR R
!
c$$$      USE MACHINE , ONLY : kind_rad
      USE HCON
      implicit none
!
      integer  LLP1,IMAX
!
!    DIMENSIONS OF ARRAYS IN ARGUMENT LIST
      real (kind=kind_rad) EMV(IMAX,LLP1), TV(IMAX,LLP1)
     &,                     AV(IMAX,LLP1),EM3V(5040)
!
      real (kind=kind_rad) FXO, TMP3, DT, FYO, DU, WW1, WW2
      integer k, i, it
!
!---THE FOLLOWING LOOP REPLACES A DOUBLE LOOP OVER I (1-IMAX) AND
!   K (1-LLP1)
!
      DO K=1,LLP1
        DO I=1,IMAX
          FXO  = AINT(TV(I,K)*HP1)
          TMP3 = LOG10(AV(I,K)) + H16E1
          DT   = TV(I,K) - TEN*FXO
          FYO  = AINT(TMP3*TEN)
          DU   = TMP3 - HP1*FYO
!---OBTAIN INDEX FOR TABLE LOOKUP; THIS VALUE WILL HAVE TO BE
!   DECREMENTED BY 9 TO ACCOUNT FOR TABLE TEMPS STARTING AT 100K.
          IT   = FXO + FYO*H28E1
          WW1  = TEN - DT
          WW2  = HP1 - DU
          EMV(I,K) = WW2 * (WW1*EM3V(IT-9)  + DT*EM3V(IT-8))
     &             + DU  * (WW1*EM3V(IT+19) + DT*EM3V(IT+20))
        ENDDO
      ENDDO
!
      RETURN
      END
!     *****************************************************************
!          SUBROUTINE FST88 IS THE MAIN COMPUTATION MODULE OF THE
!     LONG-WAVE RADIATION CODE. IN IT ALL "EMISSIVITY" CALCULATIONS,
!     INCLUDING CALLS TO TABLE LOOKUP SUBROUTINES. ALSO,AFTER CALLING
!     SUBROUTINE "SPA88", FINAL COMBINED HEATING RATES AND GROUND
!     FLUX ARE OBTAINED.
!     *****************************************************************
!              INPUTS:
!        BETINW,BETAWD,AB15WD              BDWIDE
!        BETAD,BO3RND,AO3RND               BANDTA
!        CLDFAC                            CLDCOM
!        QH2O,P,DELP2,DELP,T,VAR1,VAR2,    KDACOM
!        VAR3,VAR4,CNTVAL                  KDACOM
!        TOTVO2,TOTO3,TOTPHI,EMPL,EMX1     KDACOM
!        TPHIO3,EMX2                       KDACOM
!        TEMP,PRESS                        RADISW
!        NCLDS,KTOP,KBTM,CAMT              RADISW
!        IND,INDX2,KMAXV,SOURCE,DSRCE      TABCOM
!        SKC1R,SKC3R,KMAXVM,NREP1,NREP2    TABCOM
!        NST1,NST2,NRP1,NRP2               TABCOM
!        CO2NBL,CO21                       TFCOM
!        CO2SP1,CO2SP2                     TFCOM
!              OUTPUTS:
!        HEATRA,GRNFLX,TOPFLX              LWOUT
!
!          CALLED BY  :    RADMN OR MAIN PGM
!          CALLS      :    CLO88,E1E288,E3V88,SPA88,NLTE
!
!        PASSED VARIABLES:
!              IN E3V88:
!        EMD     =  E3 FUNCTION FOR H2O LINES (0-560,1200-2200 CM-1)
!                     COMPUTED IN E3V88
!        TPL     =  TEMPERATURE INPUT FOR E3 CALCULATION IN E3V88
!        EMPL    =  H2O AMOUNT,INPUT FOR E3 CALCULATION IN E3V88
!                   (COMPUTED IN LWR88; STORED IN KDACOM.H)
!              IN E1E288:
!        E1CTS1  =  E1 FUNCTION FOR THE (I+1)TH LEVEL USING THE
!                   TEMPERATURE OF THE ITH DATA LEVEL,COMPUTED OVER
!                   THE FREQUENCY RANGE 0-560,1200-2200 CM-1. (E1CTS1-
!                   E1CTW1) IS USED IN OBTAINING THE FLUX AT THE TOP
!                   IN THE 0-160,1200-2200 CM-1 RANGE (FLX1E1).
!        E1CTS2  =  E1 FUNCTION FOR THE ITH LEVEL, USING THE TEMP. OF
!                   THE ITH DATA LEVEL,COMPUTED OVER THE FREQUENCY RANGE
!                   0-560,1200-2200 CM-1. (E1CTS2-E1CTW2) IS ALSO USED
!                   IN OBTAINING THE FLUX AT THE TOP IN THE 0-160,.
!                   1200-2200 CM-1 RANGE.
!        E1FLX   =  E1 FCTN. FOR THE ITH LEVEL,USING THE TEMPERATURE AT
!                   THE TOP OF THE ATMOSPHERE. COMPUTED OVER THE FREQ.
!                   RANGE 0-560,1200-2200 CM-1. USED FOR Q(APPROX) TERM.
!                   (IN MODULE BLOCK TFCOM)
!        E1CTW1  =  LIKE E1CTS1,BUT COMPUTED OVER THE 160-560 CM-1 RANGE
!                   AND USED FOR Q(APPROX,CTS) CALCULATION
!        E1CTW2  =  LIKE E1CTS2,BUT COMPUTED OVER THE 160-560 CM-1 RANGE
!                   AND USED FOR Q(APPROX,CTS) CALCULATION
!        FXO     =  TEMPERATURE INDEX USED FOR E1 FUNCTION AND ALSO
!                   USED FOR SOURCE FUNCTION CALC. IN FST88.
!        DT      =  TEMP. DIFF.BETWEEN MODEL TEMPS. AND TEMPS. AT
!                   TABULAR VALUES OF E1 AND SOURCE FCTNS. USED IN
!                   FST88 AND IN E1 FUNCTION CALC.
!        FXOE2   =  TEMPERATURE INDEX USED FOR E2 FUNCTION
!        DTE2    =  TEMP. DIFF. BETWEEN MODEL TEMP. AND TEMPS. AT
!                   TABULAR VALUES OF E2 FUNCTION.
      SUBROUTINE FST88(HEATRA,GRNFLX,TOPFLX,
     &                 QH2O,PRESS,P,DELP,DELP2,TEMP,T,
     &                 CLDFAC,
!    &                 CLDFAC,NCLDS,KTOP,KBTM,CAMT,
     &                 CO21,CO2NBL,CO2SP1,CO2SP2,
     &                 VAR1,VAR2,VAR3,VAR4,CNTVAL,
     &                 TOTO3,TPHIO3,TOTPHI,TOTVO2,
     &                 EMX1,EMX2,EMPL
     &,                L, LP1, LP1V, LLP1, IMAX
     &,                SOURCE,DSRCE)
!
CFPP$ NOCONCUR R
!
c$$$      USE MACHINE , ONLY : kind_rad
      USE HCON
      USE RNDDTA
      implicit none
!
      integer L, LP1, LP1V, LLP1, IMAX
!
      real (kind=kind_rad) SOURCE(28,NBLY), DSRCE(28,NBLY)
!
      real (kind=kind_rad) QH2O(IMAX,LP1), PRESS(IMAX,LP1)
     &,                     P(IMAX,LP1),    DELP(IMAX,L),  DELP2(IMAX,L)
     &,                     TEMP(IMAX,LP1), T(IMAX,LP1)
     &,                     CLDFAC(IMAX,LP1,LP1)
     &,                     CO21(IMAX,LP1,LP1), CO2NBL(IMAX,L)
     &,                     CO2SP1(IMAX,LP1),   CO2SP2(IMAX,LP1)
     &,                     VAR1(IMAX,L),       VAR2(IMAX,L)
     &,                     VAR3(IMAX,L),       VAR4(IMAX,L)
     &,                     CNTVAL(IMAX,LP1),   TOPFLX(IMAX)
     &,                     HEATRA(IMAX,L),     GRNFLX(IMAX)
!
      real (kind=kind_rad) GXCTS(IMAX),      FLX1E1(IMAX)
     &,                     AVEPHI(IMAX,LP1), EMISS(IMAX,LP1)
     &,                     EMISSB(IMAX,LP1)
     &,                    TOTO3(IMAX,LP1),   TPHIO3(IMAX,LP1)
     &,                    TOTPHI(IMAX,LP1),  TOTVO2(IMAX,LP1)
     &,                    EMX1(IMAX),        EMX2(IMAX)
     &,                    EMPL(IMAX,LLP1)
!
      real (kind=kind_rad) EXCTS(IMAX,L),   CTSO3(IMAX,L), CTS(IMAX,L)
     &,                     E1FLX(IMAX,LP1), CO2SP(IMAX,LP1)
     &,                     TO3SPC(IMAX,L),  TO3SP(IMAX,LP1)
      real (kind=kind_rad) OSS(IMAX,LP1),   CSS(IMAX,LP1)
     &,                     SS1(IMAX,LP1),   SS2(IMAX,LP1)
     &,                     TC(IMAX,LP1),    DTC(IMAX,LP1)
     &,                     SORC(IMAX,LP1,NBLY), CSOUR(IMAX,LP1)
!CC
      real (kind=kind_rad) AVVO2(IMAX,LP1), OVER1D(IMAX,LP1)
     &,                     TO31D(IMAX,LP1), CONT1D(IMAX,LP1)
     &,                     AVMO3(IMAX,LP1), AVPHO3(IMAX,LP1)
     &,                     C(IMAX,LLP1),    C2(IMAX,LLP1)
     &,                     EMSPEC(IMAX,2)
!
!---DIMENSION OF VARIABLES EQUIVALENCED TO THOSE IN VTEMP---
      real (kind=kind_rad) VTMP3(IMAX,LP1)
     &,                     ALP(IMAX,LLP1)     ! TPL used in place of CSUB
     &,                     DELPR1(IMAX,LP1), DELPR2(IMAX,LP1)
     &,                     EMISDG(IMAX,LP1), CONTDG(IMAX,LP1)
     &,                     TO3DG(IMAX,LP1),  FLXNET(IMAX,LP1)
     &,                     VSUM1(IMAX,LP1)
!
!---DIMENSION OF VARIABLES PASSED TO OTHER SUBROUTINES---
!   (AND NOT FOUND IN MODULE BLOCKS)
      real (kind=kind_rad) E1CTS1(IMAX,LP1), E1CTS2(IMAX,L)
     &,                     E1CTW1(IMAX,LP1), E1CTW2(IMAX,L)
     &,                     TPL(IMAX,LLP1)        ! TPL used as EMD as well
!   IT IS POSSIBLE TO EQUIVALENCE EMD,TPL TO THE ABOVE VARIABLES,
!   AS THEY GET CALLED AT DIFFERENT TIMES
      real (kind=kind_rad) FXO(IMAX,LP1),   DT(IMAX,LP1)
     &,                     FXOE2(IMAX,LP1), DTE2(IMAX,LP1)
     &,                     FXOSP(IMAX,2),   DTSP(IMAX,2)
!
!     DIMENSION OF LOCAL VARIABLES
      real (kind=kind_rad) RLOG(IMAX,L),     FLX(IMAX,LP1)
     &,                     TOTEVV(IMAX,LP1), CNTTAU(IMAX,LP1)
!
      real (kind=kind_rad) vtmp, fac1, tem, tmp3, du, fyo, dt3
     &,                     ww1, ww2, fxo3, csub2
      integer lm1, ll, llm1, llm2, i, k, item, k1, kk, klen
     &,       kk1, kkk, kp, ival, it
cc
cc--------------------------------------------------------------------
cc
      real(kind=kind_rad) cons_1pdm25          !constant
cc
      cons_1pdm25      =       1.d-25          !constant
cc
cc--------------------------------------------------------------------
cc
!
!          FIRST SECTION IS TABLE LOOKUP FOR SOURCE FUNCTION AND
!     DERIVATIVE (B AND DB/DT).ALSO,THE NLTE CO2 SOURCE FUNCTION
!     IS OBTAINED
!
!---IN CALCS. BELOW, DECREMENTING THE INDEX BY 9
!   ACCOUNTS FOR THE TABLES BEGINNING AT T=100K.
!   AT T=100K.
!
      LM1  = L - 1
      LL   = LLP1 - 1
      LLM1 = LL - 1
      LLM2 = LL - 2
!
!                                 ******* E1 SOURCE *******
      DO K=1,LP1
        DO I=1,IMAX
          VTMP         = AINT(TEMP(I,K)*HP1)
          FXO(I,K)     = VTMP      - 9.0
          DT(I,K)      = TEMP(I,K) - TEN*VTMP
!
          ITEM         = FXO(I,K)
!
!       SOURCE FUNCTION FOR 14 COMBINED BANDS
!         BAND 9  - (560-670 CM-1)  BAND 10 - (670-800 CM-1)
!         BAND 11 - (800-900 CM-1)  BAND 12 - (900-990 CM-1)
!         BAND 13 - (990-1070 CM-1) BAND 14 - (1070-1200 CM-1)
!
          SORC(I,K,1)  = SOURCE(ITEM,1)  + DT(I,K)*DSRCE(ITEM,1)  ! Band 1
          SORC(I,K,2)  = SOURCE(ITEM,2)  + DT(I,K)*DSRCE(ITEM,2)  ! Band 2
          SORC(I,K,3)  = SOURCE(ITEM,3)  + DT(I,K)*DSRCE(ITEM,3)  ! Band 3
          SORC(I,K,4)  = SOURCE(ITEM,4)  + DT(I,K)*DSRCE(ITEM,4)  ! Band 4
          SORC(I,K,5)  = SOURCE(ITEM,5)  + DT(I,K)*DSRCE(ITEM,5)  ! Band 5
          SORC(I,K,6)  = SOURCE(ITEM,6)  + DT(I,K)*DSRCE(ITEM,6)  ! Band 6
          SORC(I,K,7)  = SOURCE(ITEM,7)  + DT(I,K)*DSRCE(ITEM,7)  ! Band 7
          SORC(I,K,8)  = SOURCE(ITEM,8)  + DT(I,K)*DSRCE(ITEM,8)  ! Band 8
          SORC(I,K,9)  = SOURCE(ITEM,9)  + DT(I,K)*DSRCE(ITEM,9)  ! Band 9
          SORC(I,K,10) = SOURCE(ITEM,10) + DT(I,K)*DSRCE(ITEM,10) ! Band 10
          SORC(I,K,11) = SOURCE(ITEM,11) + DT(I,K)*DSRCE(ITEM,11) ! Band 11
          SORC(I,K,12) = SOURCE(ITEM,12) + DT(I,K)*DSRCE(ITEM,12) ! Band 12
          SORC(I,K,13) = SOURCE(ITEM,13) + DT(I,K)*DSRCE(ITEM,13) ! Band 13
          SORC(I,K,14) = SOURCE(ITEM,14) + DT(I,K)*DSRCE(ITEM,14) ! Band 14
        ENDDO
      ENDDO
!
!---TEMP. INDICES FOR E2 (KP=1 LAYER NOT USED IN FLUX CALCULATIONS)
      DO K=1,L
        DO I=1,IMAX
          VTMP       = AINT(T(I,K+1)*HP1)
          FXOE2(I,K) = VTMP   - 9.0
          DTE2(I,K)  = T(I,K+1) - TEN*VTMP
        ENDDO
      ENDDO
!---SPECIAL CASE TO HANDLE KP=LP1 LAYER AND SPECIAL E2 CALCS.
      DO I=1,IMAX
        FXOE2(I,LP1) = FXO(I,L)
        DTE2(I,LP1)  = DT(I,L)
        FXOSP(I,1)   = FXOE2(I,LM1)
        FXOSP(I,2)   = FXO(I,LM1)
        DTSP(I,1)    = DTE2(I,LM1)
        DTSP(I,2)    = DT(I,LM1)
      ENDDO
!
!        THE FOLLOWING SUBROUTINE OBTAINS NLTE SOURCE FUNCTION FOR CO2
!
!     CALL NLTE
!
!---OBTAIN SPECIAL SOURCE FUNCTIONS FOR THE 15 UM BAND (CSOUR)
!   AND THE WINDOW REGION (SS1).  ALSO
!---COMPUTE TEMP**4 (TC) AND VERTICAL TEMPERATURE DIFFERENCES
!   (OSS,CSS,SS2,DTC). ALL THESE WILL BE USED LATER IN FLUX COMPUTATIONS.
!
      DO K=1,LP1
        DO I=1,IMAX
          SS1(I,K)   = SORC(I,K,11) + SORC(I,K,12) + SORC(I,K,14)
          CSOUR(I,K) = SORC(I,K,9) + SORC(I,K,10)
          VTMP       = TEMP(I,K) * TEMP(I,K)
          TC(I,K)    = VTMP * VTMP
        ENDDO
      ENDDO
      DO K=1,L
        K1 = K + 1
        DO I=1,IMAX
          OSS(I,K1) = SORC(I,K1,13) - SORC(I,K,13)
          CSS(I,K1) = CSOUR(I,K1)   - CSOUR(I,K)
          DTC(I,K1) = TC(I,K1)      - TC(I,K)
          SS2(I,K1) = SS1(I,K1)     - SS1(I,K)
        ENDDO
      ENDDO
!
!
!---THE FOLLOWIMG IS A DRASTIC REWRITE OF THE RADIATION CODE TO
!    (LARGELY) ELIMINATE THREE-DIMENSIONAL ARRAYS. THE CODE WORKS
!    ON THE FOLLOWING PRINCIPLES:
!
!          LET K = FIXED FLUX LEVEL, KP = VARYING FLUX LEVEL
!          THEN FLUX(K)=SUM OVER KP : (DELTAB(KP)*TAU(KP,K))
!               OVER ALL KP'S, FROM 1 TO LP1.
!
!          WE CAN BREAK DOWN THE CALCULATIONS FOR ALL K'S AS FOLLOWS:
!
!          FOR ALL K'S K=1 TO LP1:
!              FLUX(K)=SUM OVER KP : (DELTAB(KP)*TAU(KP,K))  (1)
!                      OVER ALL KP'S, FROM K+1 TO LP1
!          AND
!              FOR KP FROM K+1 TO LP1:
!                 FLUX(KP) = DELTAB(K)*TAU(K,KP)              (2)
!
!          NOW IF TAU(K,KP)=TAU(KP,K) (SYMMETRICAL ARRAYS)
!          WE CAN COMPUTE A 1-DIMENSIONAL ARRAY TAU1D(KP) FROM
!          K+1 TO LP1, EACH TIME K IS INCREMENTED.
!          EQUATIONS (1) AND (2) THEN BECOME:
!
!             TAU1D(KP) = (VALUES FOR TAU(KP,K) AT THE PARTICULAR K)
!             FLUX(K) = SUM OVER KP : (DELTAB(KP)*TAU1D(KP))   (3)
!             FLUX(KP) = DELTAB(K)*TAU1D(KP)                   (4)
!
!         THE TERMS FOR TAU (K,K) AND OTHER SPECIAL TERMS (FOR
!         NEARBY LAYERS) MUST, OF COURSE, BE HANDLED SEPARATELY, AND
!         WITH CARE.
!
!      COMPUTE "UPPER TRIANGLE" TRANSMISSION FUNCTIONS FOR
!      THE 9.6 UM BAND (TO3SP) AND THE 15 UM BAND (OVER1D). ALSO,
!      THE
!      STAGE 1...COMPUTE O3 ,OVER TRANSMISSION FCTNS AND AVEPHI
!---DO K=1 CALCULATION (FROM FLUX LAYER KK TO THE TOP) SEPARATELY
!   AS VECTORIZATION IS IMPROVED,AND OZONE CTS TRANSMISSIVITY
!   MAY BE EXTRACTED HERE.
!
      DO K=1,L
        DO I=1,IMAX
          AVEPHI(I,K) = TOTPHI(I,K+1)
        ENDDO
      ENDDO
!
!---IN ORDER TO PROPERLY EVALUATE EMISS INTEGRATED OVER THE (LP1)
!   LAYER, A SPECIAL EVALUATION OF EMISS IS DONE. THIS REQUIRES
!   A SPECIAL COMPUTATION OF AVEPHI, AND IT IS STORED IN THE
!   (OTHERWISE VACANT) LP1'TH POSITION
!
      DO I=1,IMAX
        AVEPHI(I,LP1) = AVEPHI(I,LM1) + EMX1(I)
      ENDDO
!
!   COMPUTE FLUXES FOR K=1
!
      CALL E1E290(E1CTS1,E1CTS2,E1FLX,E1CTW1,E1CTW2,EMISS,
     &            FXO,DT,FXOE2,DTE2,AVEPHI
     &,           L, LP1, IMAX, EM1V, EM1VW, T1, T2, T4)
!
      DO K=1,L
        K1 = K + 1
        DO I=1,IMAX
          FAC1        = BO3RND(2)*TPHIO3(I,K1)/TOTO3(I,K1)
          TO3SPC(I,K) = HAF*(FAC1*
     &        (SQRT(ONE+(FOUR*AO3RND(2)*TOTO3(I,K1))/FAC1)-ONE))
!
!   FOR K=1, TO3SP IS USED INSTEAD OF TO31D (THEY ARE EQUAL IN THIS
!   CASE); TO3SP IS PASSED TO SPA90, WHILE TO31D IS A WORK-ARRAY.
!
          TO3SP(I,K)  = EXP(HM1EZ*(TO3SPC(I,K)+SKO3R*TOTVO2(I,K1)))
          OVER1D(I,K) = EXP(HM1EZ*(SQRT(AB15WD*TOTPHI(I,K1))+
     &                                         SKC1R*TOTVO2(I,K1)))
!
!---BECAUSE ALL CONTINUUM TRANSMISSIVITIES ARE OBTAINED FROM THE
!  2-D QUANTITY CNTTAU (AND ITS RECIPROCAL TOTEVV) WE STORE BOTH
!  OF THESE HERE. FOR K=1, CONT1D EQUALS CNTTAU
!
          CNTTAU(I,K) = EXP(HM1EZ*TOTVO2(I,K1))
          TOTEVV(I,K) = 1. / max(CNTTAU(I,K),cons_1pdm25)     !constant
!         TOTEVV(I,K) = 1. / CNTTAU(I,K) ! commenteed by Moorthi 03/08/2000
        ENDDO
      ENDDO
      DO K=1,L
        DO I=1,IMAX
          CO2SP(I,K+1) = OVER1D(I,K)*CO21(I,1,K+1)
        ENDDO
      ENDDO
      DO K=1,L
        K1 = K + 1
        DO I=1,IMAX
          CO21(I,K1,1) = CO21(I,K1,1)*OVER1D(I,K)
        ENDDO
      ENDDO
!---RLOG IS THE NBL AMOUNT FOR THE 15 UM BAND CALCULATION
      DO I=1,IMAX
        RLOG(I,1) = OVER1D(I,1)*CO2NBL(I,1)
      ENDDO
!---THE TERMS WHEN KP=1 FOR ALL K ARE THE PHOTON EXCHANGE WITH
!   THE TOP OF THE ATMOSPHERE, AND ARE OBTAINED DIFFERENTLY THAN
!   THE OTHER CALCULATIONS
      DO K=2,LP1
        DO I=1,IMAX
          TEM  = TC(I,1)*E1FLX(I,K)        + SS1(I,1)*CNTTAU(I,K-1)
     &         + SORC(I,1,13)*TO3SP(I,K-1) + CSOUR(I,1)*CO2SP(I,K)
          FLX(I,K)  = TEM * CLDFAC(I,1,K)
        ENDDO
      ENDDO
 
      DO I=1,IMAX
        FLX(I,1)  = TC(I,1)*E1FLX(I,1) + SS1(I,1) + SORC(I,1,13)
     &            + CSOUR(I,1)
      ENDDO
!---THE KP TERMS FOR K=1...
      DO KP=2,LP1
        DO I=1,IMAX
          TEM      = OSS(I,KP)*TO3SP(I,KP-1) + SS2(I,KP)*CNTTAU(I,KP-1)
     &             + CSS(I,KP)*CO21(I,KP,1)  + DTC(I,KP)*EMISS(I,KP-1)
          FLX(I,1)  = FLX(I,1)  + TEM * CLDFAC(I,KP,1)
        ENDDO
      ENDDO
!
!          SUBROUTINE SPA88 IS CALLED TO OBTAIN EXACT CTS FOR WATER
!     CO2 AND O3, AND APPROXIMATE CTS CO2 AND O3 CALCULATIONS.
!
      CALL SPA88(EXCTS,CTSO3,GXCTS,SORC,CSOUR,
     &           CLDFAC,TEMP,PRESS,VAR1,VAR2,
     &           P,DELP,DELP2,TOTVO2,TO3SP,TO3SPC,
     &           CO2SP1,CO2SP2,CO2SP
     &,          L, LP1, IMAX)
!
!    THIS SECTION COMPUTES THE EMISSIVITY CTS HEATING RATES FOR 2
!    EMISSIVITY BANDS: THE 0-160,1200-2200 CM-1 BAND AND THE 800-
!    990,1070-1200 CM-1 BAND. THE REMAINING CTS COMTRIBUTIONS ARE
!    CONTAINED IN CTSO3, COMPUTED IN SPA88.
!
      DO I=1,IMAX
        VTMP3(I,1) = 1.
      ENDDO
      DO K=1,L
        DO I=1,IMAX
          VTMP3(I,K+1) = CNTTAU(I,K)*CLDFAC(I,K+1,1)
        ENDDO
      ENDDO
      DO K=1,L
        DO I=1,IMAX
          CTS(I,K) = TC(I,K) * (E1CTW2(I,K)*CLDFAC(I,K+1,1)
     &                        - E1CTW1(I,K)*CLDFAC(I,K,1)) +
     &                          SS1(I,K)*(VTMP3(I,K+1)-VTMP3(I,K))
!
        ENDDO
      ENDDO
!
      DO K=1,L
        DO I=1,IMAX
          VTMP3(I,K)=TC(I,K)*(CLDFAC(I,K,1)*(E1CTS1(I,K)-E1CTW1(I,K)) -
     &                        CLDFAC(I,K+1,1)*(E1CTS2(I,K)-E1CTW2(I,K)))
!
        ENDDO
      ENDDO
      DO I=1,IMAX
        TEM = TC(I,LP1) * (E1CTS1(I,LP1)-E1CTW1(I,LP1))
        FLX1E1(I) = TEM * CLDFAC(I,LP1,1)
      ENDDO
      DO K=1,L
        DO I=1,IMAX
          FLX1E1(I) = FLX1E1(I) + VTMP3(I,K)
        ENDDO
      ENDDO
!
!---NOW REPEAT FLUX CALCULATIONS FOR THE K=2..LM1  CASES.
!   CALCULATIONS FOR FLUX LEVEL L AND LP1 ARE DONE SEPARATELY, AS ALL
!   EMISSIVITY AND CO2 CALCULATIONS ARE SPECIAL CASES OR NEARBY LAYERS.
!
      DO K=2,LM1
        KLEN = K
        DO KK=1,LP1-K
          DO I=1,IMAX
            AVEPHI(I,KK+K-1) = TOTPHI(I,KK+K) - TOTPHI(I,K)
          ENDDO
        ENDDO
        DO I=1,IMAX
          AVEPHI(I,LP1) = AVEPHI(I,LM1) + EMX1(I)
        ENDDO
!
!---COMPUTE EMISSIVITY FLUXES (E2) FOR THIS CASE. NOTE THAT
!   WE HAVE OMITTED THE NEARBY LATER CASE (EMISS(I,K,K)) AS WELL
!   AS ALL CASES WITH K=L OR LP1. BUT THESE CASES HAVE ALWAYS
!   BEEN HANDLED AS SPECIAL CASES, SO WE MAY AS WELL COMPUTE
!    THEIR FLUXES SEPARASTELY.
!
        CALL E290(EMISSB,EMISS,AVEPHI,KLEN,FXOE2,DTE2
     &,                 L, LP1, IMAX, T1, T2, T4)
!
        DO KK=1,LP1-K
          KKK = KK + K
          KK1 = KKK - 1
          DO I=1,IMAX
            AVMO3(I,KK1)  = TOTO3(I,KKK)  - TOTO3(I,K)
            AVPHO3(I,KK1) = TPHIO3(I,KKK) - TPHIO3(I,K)
            AVVO2(I,KK1)  = TOTVO2(I,KKK) - TOTVO2(I,K)
            CONT1D(I,KK1) = CNTTAU(I,KK1) * TOTEVV(I,K-1)
          ENDDO
        ENDDO
!
CDIR$ IVDEP
        DO KK=1,LP1-K
          KKK = KK  + K
          KK1 = KKK - 1
          DO I=1,IMAX
            FAC1 = BO3RND(2)*AVPHO3(I,KK1)/AVMO3(I,KK1)
            VTMP = HAF*(FAC1*(SQRT(ONE+(FOUR*AO3RND(2)*AVMO3(I,KK1))
     &                                  /FAC1)-ONE))
            TO31D(I,KK1)  = EXP(HM1EZ*(VTMP+SKO3R*AVVO2(I,KK1)))
            OVER1D(I,KK1) = EXP(HM1EZ*(SQRT(AB15WD*AVEPHI(I,KK1))+
     &                                 SKC1R*AVVO2(I,KK1)))
            CO21(I,KKK,K) = OVER1D(I,KK1)*CO21(I,KKK,K)
          ENDDO
        ENDDO
        DO KP=K+1,LP1
          DO I=1,IMAX
            CO21(I,K,KP) = OVER1D(I,KP-1)*CO21(I,K,KP)
          ENDDO
        ENDDO
!---RLOG IS THE NBL AMOUNT FOR THE 15 UM BAND CALCULATION
        DO I=1,IMAX
          RLOG(I,K) = OVER1D(I,K)*CO2NBL(I,K)
        ENDDO
!---THE KP TERMS FOR ARBIRRARY K..
        DO KP=K+1,LP1
          DO I=1,IMAX
            TEM       = OSS(I,KP)*TO31D(I,KP-1)+SS2(I,KP)*CONT1D(I,KP-1)
     &                + CSS(I,KP)*CO21(I,KP,K) +DTC(I,KP)*EMISS(I,KP-1)
            FLX(I,K)  = FLX(I,K)  + TEM * CLDFAC(I,KP,K)
          ENDDO
        ENDDO
        DO KP=K+1,LP1
          DO I=1,IMAX
            TEM        = OSS(I,K)*TO31D(I,KP-1)+SS2(I,K)*CONT1D(I,KP-1)
     &                 + CSS(I,K)*CO21(I,K,KP) +DTC(I,K)*EMISSB(I,KP-1)
            FLX(I,KP)  = FLX(I,KP)  + TEM * CLDFAC(I,K,KP)
          ENDDO
        ENDDO
      ENDDO
!
!   NOW DO K=L CASE. SINCE THE KP LOOP IS LENGTH 1, MANY SIMPLIFI-
!   CATIONS OCCUR. ALSO, THE CO2 QUANTITIES (AS WELL AS THE EMISS
!  QUANTITIES) ARE COMPUTED IN THE NBL SEDCTION; THEREFORE, WE WANT
!  ONLY OVER,TO3 AND CONT1D (OVER(I,L),TO31D(I,L) AND CONT1D(I,L)
!  ACCORDING TO THE NOTATION. THUS NO CALL IS MADE TO THE E290
!  SUBROUTINE.
!         THE THIRD SECTION CALCULATES BOUNDARY LAYER AND NEARBY LAYER
!     CORRECTIONS TO THE TRANSMISSION FUNCTIONS OBTAINED ABOVE. METHODS
!     ARE GIVEN IN REF. (4).
!          THE FOLLOWING RATIOS ARE USED IN VARIOUS NBL CALCULATIONS:
!
!   THE REMAINING CALCULATIONS ARE FOR :
!                        1) THE (K,K) TERMS, K=2,LM1;
!                        2) THE (L,L) TERM
!                        3) THE (L,LP1) TERM
!                        4) THE (LP1,L) TERM
!                        5) THE (LP1,LP1) TERM.
!     EACH IS UNIQUELY HANDLED; DIFFERENT FLUX TERMS ARE COMPUTED
!     DIFFERENTLY
!
!
!          FOURTH SECTION OBTAINS WATER TRANSMISSION FUNCTIONS
!     USED IN Q(APPROX) CALCULATIONS AND ALSO MAKES NBL CORRECTIONS:
!     1) EMISS (I,J) IS THE TRANSMISSION FUNCTION MATRIX OBTAINED
!     BY CALLING SUBROUTINE E1E288;
!     2) "NEARBY LAYER" CORRECTIONS (EMISS(I,I)) ARE OBTAINED
!     USING SUBROUTINE E3V88;
!     3) SPECIAL VALUES AT THE SURFACE (EMISS(L,LP1),EMISS(LP1,L),
!     EMISS(LP1,LP1)) ARE CALCULATED.
!
!
!      OBTAIN ARGUMENTS FOR E1E288 AND E3V88:
!
      DO I=1,IMAX
        TPL(I,1)    = TEMP(I,L)
        TPL(I,LP1)  = HAF*(T(I,LP1) + TEMP(I,L))
        TPL(I,LLP1) = HAF*(T(I,L)   + TEMP(I,L))
!
!---E2 FUNCTIONS ARE REQUIRED IN THE NBL CALCULATIONS FOR 2 CASES,
!   DENOTED (IN OLD CODE) AS (L,LP1) AND (LP1,LP1)
        AVEPHI(I,1) = VAR2(I,L)
        AVEPHI(I,2) = VAR2(I,L) + EMPL(I,L)
      ENDDO
      DO K=2,L
        DO I=1,IMAX
          TPL(I,K)   = T(I,K)
          TPL(I,K+L) = T(I,K)
        ENDDO
      ENDDO
!
!     Inlining of E2SPEC
!
!     SUBROUTINE E2SPEC COMPUTES THE EXCHANGE TERMS IN THE FLUX EQUATION
!  FOR LONGWAVE RADIATION FOR 2 TERMS USED FOR NEARBY LAYER COMPU-
!  TATIONS. THE METHOD IS A TABLE LOOKUP ON A PRE-
!  COMPUTED E2 FUNCTION (DEFINED IN REF. (4)).
!
      DO  K=1,2
        DO  I=1,IMAX
          TMP3       = LOG10(AVEPHI(I,K)) + H16E1
          FYO        = AINT(TMP3*TEN)
          DU         = TMP3 - HP1*FYO
          IVAL       = H28E1*FYO + FXOSP(I,K)
          EMISS(I,K) = T1(IVAL)  + DU*T2(IVAL) + DTSP(I,K)*T4(IVAL)
        ENDDO
      ENDDO
!
!     CALL E3V88 FOR NBL H2O TRANSMISSIVITIES
!     SUBROUTINE E3V88 COMPUTES NEARBY LAYER TRANSMISSIVITIES FOR
!  H2O USING A TABLE LOOKUP OF THE PRE-COMPUTED E3 FUNCTION
! ( DESCRIBED IN REF. (4)).
!
!     Inlining of E3V88
!
      DO K=1,LLP1
        DO I=1,IMAX
          FXO3 = AINT(TPL(I,K)*HP1)
          TMP3 = LOG10(EMPL(I,K)) + H16E1
          DT3  = TPL(I,K) - TEN*FXO3
          FYO  = AINT(TMP3*TEN)
          DU   = TMP3 - HP1*FYO
!---OBTAIN INDEX FOR TABLE LOOKUP; THIS VALUE WILL HAVE TO BE
!   DECREMENTED BY 9 TO ACCOUNT FOR TABLE TEMPS STARTING AT 100K.
          IT   = FXO3 + FYO*H28E1
          WW1  = TEN - DT3
          WW2  = HP1 - DU
          TPL(I,K) = WW2 * (WW1*EM3V(IT-9)  + DT3*EM3V(IT-8))+
     &               DU  * (WW1*EM3V(IT+19) + DT3*EM3V(IT+20))
        ENDDO
      ENDDO
!
!   COMPUTE NEARBY LAYER AND SPECIAL-CASE TRANSMISSIVITIES FOR EMISS
!    USING METHODS FOR H2O GIVEN IN REF. (4)
CDIR$ IVDEP
      DO K=2,L
        DO I=1,IMAX
          EMISDG(I,K) = TPL(I,K+L) + TPL(I,K)
        ENDDO
      ENDDO
!
!   NOTE THAT EMX1/2 (PRESSURE SCALED PATHS) ARE NOW COMPUTED IN
!   LWR88
      DO I=1,IMAX
        EMSPEC(I,1) = (TPL(I,1)*EMPL(I,1)-TPL(I,LP1)*EMPL(I,LP1))/
     &                 EMX1(I) + QUARTR*(EMISS(I,1)+EMISS(I,2))
        EMISDG(I,LP1)=TWO*TPL(I,LP1)
        EMSPEC(I,2) = TWO*(TPL(I,1)*EMPL(I,1)-TPL(I,LLP1)*EMPL(I,LLP1))/
     &                EMX2(I)
      ENDDO
      DO I=1,IMAX
        FAC1 = BO3RND(2)*VAR4(I,L)/VAR3(I,L)
        VTMP = HAF*(FAC1*(SQRT(ONE+(FOUR*AO3RND(2)*VAR3(I,L))/FAC1)
     &                                                        -ONE))
        TO31D(I,L)  = EXP(HM1EZ*(VTMP+SKO3R*CNTVAL(I,L)))
        OVER1D(I,L) = EXP(HM1EZ*(SQRT(AB15WD*VAR2(I,L))+
     &                           SKC1R*CNTVAL(I,L)))
        CONT1D(I,L) = CNTTAU(I,L)*TOTEVV(I,LM1)
        RLOG(I,L)   = OVER1D(I,L)*CO2NBL(I,L)
      ENDDO
      DO K=1,L
        K1 = K + 1
        DO I=1,IMAX
          RLOG(I,K)    = LOG(RLOG(I,K))
          DELPR2(I,K1) = DELP(I,K)*(P(I,K1)-PRESS(I,K))
          TPL(I,K)     = -SQRT(DELPR2(I,K1))*RLOG(I,K)
        ENDDO
      ENDDO
      DO K=1,LM1
        K1 = K + 1
        DO I=1,IMAX
          DELPR1(I,K1) = DELP(I,K1)*(PRESS(I,K1)-P(I,K1))
          TPL(I,K+L)   = -SQRT(DELPR1(I,K1))*RLOG(I,K1)
        ENDDO
      ENDDO
      DO I=1,IMAX
        TPL(I,LL)   = -RLOG(I,L)
        TPL(I,LLP1) = -RLOG(I,L)*SQRT(DELP(I,L)*(P(I,LP1)-PRESS(I,LM1)))
      ENDDO
!        THE FIRST COMPUTATION IS FOR THE 15 UM BAND,WITH THE
!     FOR THE COMBINED H2O AND CO2 TRANSMISSION FUNCTION.
!
!       PERFORM NBL COMPUTATIONS FOR THE 15 UM BAND
!***THE STATEMENT FUNCTION SF IN PREV. VERSIONS IS NOW EXPLICITLY
!   EVALUATED.
      DO K=1,LLP1
        DO I=1,IMAX
          C(I,K)=TPL(I,K)*(HMP66667+TPL(I,K)*(QUARTR+TPL(I,K)*HM6666M2))
        ENDDO
      ENDDO
      DO I=1,IMAX
        CO21(I,LP1,LP1) = ONE+C(I,L)
        CO21(I,LP1,L)   = ONE+(DELP2(I,L)*C(I,LL)-(PRESS(I,L)-P(I,L))*
     &                    C(I,LLM1))/(P(I,LP1)-PRESS(I,L))
        CO21(I,L,LP1)   = ONE+((P(I,LP1)-PRESS(I,LM1))*C(I,LLP1)-
     &          (P(I,LP1)-PRESS(I,L))*C(I,L))/(PRESS(I,L)-PRESS(I,LM1))
      ENDDO
      DO K=2,L
        DO I=1,IMAX
          CO21(I,K,K) = ONE + HAF*(C(I,LM1+K)+C(I,K-1))
        ENDDO
      ENDDO
!
!    COMPUTE NEARBY-LAYER TRANSMISSIVITIES FOR THE O3 BAND AND FOR THE
!    ONE-BAND CONTINUUM BAND (TO3 AND EMISS2). THE SF2 FUNCTION IS
!    USED. THE METHOD IS THE SAME AS DESCRIBED FOR CO2 IN REF (4).
CDIR$ IVDEP
      DO K=1,LM1
        K1 = K + 1
        DO I=1,IMAX
          TPL(I,K1)  = CNTVAL(I,K1) * DELPR1(I,K1)
          TPL(I,K+L) = CNTVAL(I,K)  * DELPR2(I,K1)
        ENDDO
      ENDDO
!---THE SF2 FUNCTION IN PREV. VERSIONS IS NOW EXPLICITLY EVALUATED
      DO K=1,LLM2
        DO I=1,IMAX
          TEM       = TPL(I,K+1)
          CSUB2     = SKO3R*TEM
          C(I,K+1)  = TEM  *(HMP5+TEM  *(HP166666-TEM*H41666M2))
          C2(I,K+1) = CSUB2*(HMP5+CSUB2*(HP166666-CSUB2*H41666M2))
        ENDDO
      ENDDO
      DO I=1,IMAX
        CONTDG(I,LP1) = 1. + C(I,LLM1)
        TO3DG(I,LP1)  = 1. + C2(I,LLM1)
      ENDDO
      DO K=2,L
        DO I=1,IMAX
          CONTDG(I,K) = ONE + HAF * (C(I,K)  + C(I,LM1+K))
          TO3DG(I,K)  = ONE + HAF * (C2(I,K) + C2(I,LM1+K))
        ENDDO
      ENDDO
!
!---NOW OBTAIN FLUXES
!
!    FOR THE DIAGONAL TERMS...
      DO K=2,LP1
        DO I=1,IMAX
          TEM       = DTC(I,K)*EMISDG(I,K) + SS2(I,K)*CONTDG(I,K)
     &              + OSS(I,K)*TO3DG(I,K)  + CSS(I,K)*CO21(I,K,K)
          FLX(I,K)  = FLX(I,K)  + TEM * CLDFAC(I,K,K)
        ENDDO
      ENDDO
!     FOR THE TWO OFF-DIAGONAL TERMS...
      DO I=1,IMAX
        TEM         = CSS(I,LP1)*CO21(I,LP1,L) + DTC(I,LP1)*EMSPEC(I,2)
     &              + OSS(I,LP1)*TO31D(I,L)    + SS2(I,LP1)*CONT1D(I,L)
        FLX(I,L)    = FLX(I,L)  + TEM  * CLDFAC(I,LP1,L)
!
        TEM         = CSS(I,L)*CO21(I,L,LP1) + OSS(I,L)*TO31D(I,L)
     &              + SS2(I,L)*CONT1D(I,L)   + DTC(I,L)*EMSPEC(I,1)
        FLX(I,LP1)  = FLX(I,LP1)  + TEM * CLDFAC(I,L,LP1)
      ENDDO
!
!     FINAL SECTION OBTAINS EMISSIVITY HEATING RATES,
!     TOTAL HEATING RATES AND THE FLUX AT THE GROUND
!
      DO K=1,L
        DO I=1,IMAX
!     .....CALCULATE THE TOTAL HEATING RATES
          TEM         = RADCON * DELP(I,K)
          VSUM1(I,K)  = FLX(I,K+1) - FLX(I,K) - CTS(I,K)
     &                - CTSO3(I,K) + EXCTS(I,K)
          HEATRA(I,K) = TEM * VSUM1(I,K)
!
!        print *,' heatra=',heatra(i,k),' flx=', flx(i,K+1),flx(i,k)
!    &,' cts=',cts(i,K),' ctso3=',ctso3(i,k)
!    &,' excts=',excts(i,k),' k=',k,' tem=',radcon*delp(i,k)
!
!
        ENDDO
      ENDDO
!
!     .....CALCULATE THE FLUX AT EACH FLUX LEVEL USING THE FLUX AT THE
!    TOP (FLX1E1+GXCTS) AND THE INTEGRAL OF THE HEATING RATES (VSUM1)
      DO I=1,IMAX
        TOPFLX(I)   = FLX1E1(I) + GXCTS(I)
        FLXNET(I,1) = TOPFLX(I)
      ENDDO
!---ONLY THE SURFACE VALUE OF FLUX (GRNFLX) IS NEEDED UNLESS
!    THE THICK CLOUD SECTION IS INVOKED.
!
      DO K=1,L
        DO I=1,IMAX
          FLXNET(I,K+1) = FLXNET(I,K) + VSUM1(I,K)
        ENDDO
      ENDDO
      DO I=1,IMAX
        GRNFLX(I) = FLXNET(I,LP1)
      ENDDO
!
!  ***************************************************************
!  *   THICK CLOUD SECTION NO LONGER USED ....K.A.C. SEP96       *
!  ***************************************************************
!     THIS IS THE THICK CLOUD SECTION.OPTIONALLY,IF THICK CLOUD
!     FLUXES ARE TO BE "CONVECTIVELY ADJUSTED",IE,DF/DP IS CONSTANT,
!     FOR CLOUDY PART OF GRID POINT, THE FOLLOWING CODE IS EXECUTED.
!***FIRST,COUNT THE NUMBER OF CLOUDS ALONG THE LAT. ROW. SKIP THE
!   ENTIRE THICK CLOUD COMPUTATION IF THERE ARE NO CLOUDS.
!KC      ICNT=0
!KC      DO 1301 I=1,IMAX
!KC      ICNT=ICNT+NCLDS(I)
!KC1301  CONTINUE
!KC      IF (ICNT.EQ.0) GO TO 6999
!---FIND THE MAXIMUM NUMBER OF CLOUDS IN THE LATITUDE ROW
!KC      KCLDS=NCLDS(1)
!KC      DO 2106 I=2,IMAX
!KC      KCLDS=MAX(NCLDS(I),KCLDS)
!KC2106  CONTINUE
!
!***OBTAIN THE PRESSURES AND FLUXES OF THE TOP AND BOTTOM OF
!   THE NC'TH CLOUD (IT IS ASSUMED THAT ALL KTOP AND KBTM'S HAVE
!   BEEN DEFINED!).
!KC      DO 1361 KK=1,KCLDS
!KC      KMIN=LP1
!KC      KMAX=0
!KC      DO 1362 I=1,IMAX
!KC        J1=KTOP(I,KK+1)
!       IF (J1.EQ.1) GO TO 1362
!KC        J3=KBTM(I,KK+1)
!KC        IF (J3.GT.J1) THEN
!KC          PTOP(I)=P(I,J1)
!KC          PBOT(I)=P(I,J3+1)
!KC          FTOP(I)=FLXNET(I,J1)
!KC          FBOT(I)=FLXNET(I,J3+1)
!***OBTAIN THE "FLUX DERIVATIVE" DF/DP (DELPTC)
!KC          DELPTC(I)=(FTOP(I)-FBOT(I))/(PTOP(I)-PBOT(I))
!KC          KMIN=MIN(KMIN,J1)
!KC          KMAX=MAX(KMAX,J3)
!KC        ENDIF
!KC1362  CONTINUE
!KC      KMIN=KMIN+1
!***CALCULATE THE TOT. FLUX CHG. FROM THE TOP OF THE CLOUD, FOR
!   ALL LEVELS.
!KC      DO 1365 K=KMIN,KMAX
!KC      DO 1363 I=1,IMAX
!       IF (KTOP(I,KK+1).EQ.1) GO TO 1363
!KC        IF(KTOP(I,KK+1).LT.K .AND. K.LE.KBTM(I,KK+1)) THEN
!KC          Z1(I,K)=(P(I,K)-PTOP(I))*DELPTC(I)+FTOP(I)
!ORIGINAL FLXNET(I,K)=FLXNET(I,K)*(ONE-CAMT(I,KK+1)) +
!ORIGINAL1            Z1(I,K)*CAMT(I,KK+1)
!KC          FLXNET(I,K)=Z1(I,K)
!KC        ENDIF
!KC1363  CONTINUE
!KC1365  CONTINUE
!KC1361  CONTINUE
!***USING THIS FLUX CHG. IN THE CLOUDY PART OF THE GRID BOX, OBTAIN
!   THE NEW FLUXES, WEIGHTING THE CLEAR AND CLOUDY FLUXES:AGAIN, ONLY
!    THE FLUXES IN THICK-CLOUD LEVELS WILL EVENTUALLY BE USED.
!     DO 6051 K=1,LP1
!     DO 6051 I=1,IMAX
!     FLXNET(I,K)=FLXNET(I,K)*(ONE-CAMT(I,NC)) +
!    1            Z1(I,K)*CAMT(I,NC)
!051  CONTINUE
!***MERGE FLXTHK INTO FLXNET FOR APPROPRIATE LEVELS.
!     DO 1401 K=1,LP1
!     DO 1401 I=1,IMAX
!     IF (K.GT.ITOP(I) .AND. K.LE.IBOT(I)
!    1  .AND.  (NC-1).LE.NCLDS(I))  THEN
!          FLXNET(I,K)=FLXTHK(I,K)
!     ENDIF
!401  CONTINUE
!
!******END OF CLOUD LOOP*****
!KC6001  CONTINUE
!KC6999  CONTINUE
!***THE FINAL STEP IS TO RECOMPUTE THE HEATING RATES BASED ON THE
!   REVISED FLUXES:
!KC      DO 6101 K=1,L
!KC      DO 6101 I=1,IMAX
!KC      HEATRA(I,K)=RADCON*(FLXNET(I,K+1)-FLXNET(I,K))*DELP(I,K)
!KC6101  CONTINUE
!     THE THICK CLOUD SECTION ENDS HERE.
!  ***************************************************************
!  *   THICK CLOUD SECTION NO LONGER USED ....K.A.C. SEP96       *
!  ***************************************************************
      RETURN
      END
      SUBROUTINE HCONST
CFPP$ NOCONCUR R
!
      USE HCON
      implicit none
!     SUBROUTINE HCONST DEFINES VARIABLES TO REPRESENT FLOATING-
!       POINT CONSTANTS.
!
!     COMDECK HCON CONTAINS THE MODULE BLOCK FOR THESE FLOATING-
!       POINT CONSTANTS.
!
!     THE NAMING CONVENTIONS FOR THE FLOATING-POINT VARIABLES ARE
!       AS FOLLOWS:
!
!   1) PHYSICAL AND MATHEMATICAL CONSTANTS WILL BE GIVEN NAMES
!     RELEVANT TO THEIR MEANING
!   2) OTHER CONSTANTS WILL BE GIVEN NAMES RELEVANT TO THEIR VALUE
!      AND ADHERING TO THE FOLLOWING CONVENTIONS:
!       A) THE FIRST LETTER WILL BE REPRESENTED WITH AN 'H' EXCEPT
!          FOR I) AND J) BELOW
!       B) A DECIMAL POINT WILL BE REPRESENTED WITH A 'P'
!       C) THERE WILL BE NO EMBEDDED '0'(ZERO); ALL 0'S WILL
!          BE REPRESENTED WITH A 'Z'
!       D) A MINUS SIGN WILL BE REPRESENTED WITH AN 'M'
!       E) THE DECIMAL POINT IS ASSUMED AFTER THE FIRST DIGIT FOR
!          NUMBERS WITH EXPONENTS
!       F) POSITIVE EXPONENTS ARE INDICATED WITH 'E';NEGATIVE
!          EXPONENTS WITH 'M'
!       G) DIGITS ARE TRUNCATED IN ORDER TO HAVE NO MORE THAN 8
!          CHARACTERS PER NAME
!       H) NUMBERS LESS THAN 0.1 AND GREATER THAN 10. WILL BE
!          REPRESENTED IN EXPONENT FORMAT (EXCEPT A FEW SPECIAL CASES)
!       I) THE WHOLE NUMBERS FROM 0.0 THROUGH 10.,AND 20.,30.,40.,50.,
!          60.,70.,80.,90.,100.,WILL BE SPELLED OUT
!       J) GOOD JUDGMENT WILL PREVAIL OVER ALL CONVENTIONS
!
!       EXAMPLES
!     CONSTANT           VARIABLE NAME             CONVENTION
!      600.                 LHEATC                  1)
!      680.                 LHEATS                  1)
!     1.4142               SQROOT2                  1)
!     2.0                    TWO                    2)-(I)
!    -3.0                  HM3PZ                    2)-(A,B,D)
!    310.                  C31E2                    2)-(A,E,F,H)
!   -0.7239E-9             HM723M1Z                 2)-(A,C,D,E,F,G,H)
!     0.0                   ZERO                    2)-(I)
!     0.1                   HP1                     2)-(A,B,H)
!     0.01                 H1M2                     2)-(A,E,F,H)
!     30.                  THIRTY                   2)-(H,I)
!     0.5                  HAF                      2)-(J)
!     9.0                  HNINE                    2)-(J)
!
!
!******THE FOLLOWING ARE PHYSICAL CONSTANTS*****
!        ARRANGED IN ALPHABETICAL ORDER
!
      AMOLWT   = 28.9644
      CSUBP    = 1.00484E7
      DIFFCTR  = 1.66
      G        = 980.665
      GINV     = 1./G
      GRAVDR   = 980.0
      O3DIFCTR = 1.90
      P0       = 1013250.
      P0INV    = 1./P0
      GP0INV   = GINV*P0INV
      P0XZP2   = 202649.902
      P0XZP8   = 810600.098
      P0X2     = 2.*1013250.
      RADCON   = 8.427
      RADCON1  = 1./8.427
      RATCO2MW = 1.519449738
      RATH2OMW = .622
      RGASK    = 8.3142E7
      RGASSP   = 8.31432E7
      SECPDA   = 8.64E4
!
!******THE FOLLOWING ARE MATHEMATICAL CONSTANTS*******
!        ARRANGED IN DECREASING ORDER
      HUNDRED = 100.
      HNINETY = 90.
      SIXTY   = 60.
      FIFTY   = 50.
      TEN     = 10.
      EIGHT   = 8.
      FIVE    = 5.
      FOUR    = 4.
      THREE   = 3.
      TWO     = 2.
      ONE     = 1.
      HAF     = 0.5
      QUARTR  = 0.25
      ZERO    = 0.
!
!******FOLLOWING ARE POSITIVE FLOATING POINT CONSTANTS(H'S)
!       ARRANGED IN DECREASING ORDER
      H83E26   = 8.3E26
      H71E26   = 7.1E26
      H1E15    = 1.E15
      H1E13    = 1.E13
      H1E11    = 1.E11
      H1E8     = 1.E8
      H2E6     = 2.0E6
      H1E6     = 1.0E6
      H69766E5 = 6.97667E5
      H4E5     = 4.E5
      H165E5   = 1.65E5
      H5725E4  = 57250.
      H488E4   = 48800.
      H1E4     = 1.E4
      H24E3    = 2400.
      H20788E3 = 2078.8
      H2075E3  = 2075.
      H18E3    = 1800.
      H1224E3  = 1224.
      H67390E2 = 673.9057
      H5E2     = 500.
      H3082E2  = 308.2
      H3E2     = 300.
      H2945E2  = 294.5
      H29316E2 = 293.16
      H26E2    = 260.0
      H25E2    = 250.
      H23E2    = 230.
      H2E2     = 200.0
      H15E2    = 150.
      H1386E2  = 138.6
      H1036E2  = 103.6
      H8121E1  = 81.21
      H35E1    = 35.
      H3116E1  = 31.16
      H28E1    = 28.
      H181E1   = 18.1
      H18E1    = 18.
      H161E1   = 16.1
      H16E1    = 16.
      H1226E1  = 12.26
      H9P94    = 9.94
      H6P08108 = 6.081081081
      H3P6     = 3.6
      H3P5     = 3.5
      H2P9     = 2.9
      H2P8     = 2.8
      H2P5     = 2.5
      H1P8     = 1.8
      H1P4387  = 1.4387
      H1P41819 = 1.418191
      H1P4     = 1.4
      H1P25892 = 1.258925411
      H1P082   = 1.082
      HP816    = 0.816
      HP805    = 0.805
      HP8      = 0.8
      HP60241  = 0.60241
      HP602409 = 0.60240964
      HP6      = 0.6
      HP526315 = 0.52631579
      HP518    = 0.518
      HP5048   = 0.5048
      HP3795   = 0.3795
      HP369    = 0.369
      HP26     = 0.26
      HP228    = 0.228
      HP219    = 0.219
      HP166666 = .166666
      HP144    = 0.144
      HP118666 = 0.118666192
      HP1=0.1
!        (NEGATIVE EXPONENTIALS BEGIN HERE)
      H658M2   = 0.0658
      H625M2   = 0.0625
      H44871M2 = 4.4871E-2
      H44194M2 = .044194
      H42M2    = 0.042
      H41666M2 = 0.0416666
      H28571M2 = .02857142857
      H2118M2  = 0.02118
      H129M2   = 0.0129
      H1M2     = .01
      H559M3   = 5.59E-3
      H3M3     = 0.003
      H235M3   = 2.35E-3
      H1M3     = 1.0E-3
      H987M4   = 9.87E-4
      H323M4   = 0.000323
      H3M4     = 0.0003
      H285M4   = 2.85E-4
      H1M4     = 0.0001
      H75826M4 = 7.58265E-4
      H6938M5  = 6.938E-5
      H394M5   = 3.94E-5
      H37412M5 = 3.7412E-5
      H15M5    = 1.5E-5
      H1439M5  = 1.439E-5
      H128M5   = 1.28E-5
      H102M5   = 1.02E-5
      H1M5     = 1.0E-5
      H7M6     = 7.E-6
      H4999M6  = 4.999E-6
      H451M6   = 4.51E-6
      H25452M6 = 2.5452E-6
      H1M6     = 1.E-6
      H391M7   = 3.91E-7
      H1174M7  = 1.174E-7
      H8725M8  = 8.725E-8
      H327M8   = 3.27E-8
      H257M8   = 2.57E-8
      H1M8     = 1.0E-8
      H23M10   = 2.3E-10
      H14M10   = 1.4E-10
      H11M10   = 1.1E-10
      H1M10    = 1.E-10
      H83M11   = 8.3E-11
      H82M11   = 8.2E-11
      H8M11    = 8.E-11
      H77M11   = 7.7E-11
      H72M11   = 7.2E-11
      H53M11   = 5.3E-11
      H48M11   = 4.8E-11
      H44M11   = 4.4E-11
      H42M11   = 4.2E-11
      H37M11   = 3.7E-11
      H35M11   = 3.5E-11
      H32M11   = 3.2E-11
      H3M11    = 3.0E-11
      H28M11   = 2.8E-11
      H24M11   = 2.4E-11
      H23M11   = 2.3E-11
      H2M11    = 2.E-11
      H18M11   = 1.8E-11
      H15M11   = 1.5E-11
      H14M11   = 1.4E-11
      H114M11  = 1.14E-11
      H11M11   = 1.1E-11
      H1M11    = 1.E-11
      H96M12   = 9.6E-12
      H93M12   = 9.3E-12
      H77M12   = 7.7E-12
      H74M12   = 7.4E-12
      H65M12   = 6.5E-12
      H62M12   = 6.2E-12
      H6M12    = 6.E-12
      H45M12   = 4.5E-12
      H44M12   = 4.4E-12
      H4M12    = 4.E-12
      H38M12   = 3.8E-12
      H37M12   = 3.7E-12
      H3M12    = 3.E-12
      H29M12   = 2.9E-12
      H28M12   = 2.8E-12
      H24M12   = 2.4E-12
      H21M12   = 2.1E-12
      H16M12   = 1.6E-12
      H14M12   = 1.4E-12
      H12M12   = 1.2E-12
      H8M13    = 8.E-13
      H46M13   = 4.6E-13
      H36M13   = 3.6E-13
      H135M13  = 1.35E-13
      H12M13   = 1.2E-13
      H1M13    = 1.E-13
      H3M14    = 3.E-14
      H15M14   = 1.5E-14
      H14M14   = 1.4E-14
      H101M16  = 1.01E-16
      H1M16    = 1.0E-16
      H1M17    = 1.E-17
      H1M18    = 1.E-18
      H1M19    = 1.E-19
      H1M20    = 1.E-20
      H1M21    = 1.E-21
      H1M22    = 1.E-22
      H1M23    = 1.E-23
      H1M24    = 1.E-24
      H26M30   = 2.6E-30
      H14M30   = 1.4E-30
      H25M31   = 2.5E-31
      H21M31   = 2.1E-31
      H12M31   = 1.2E-31
      H9M32    = 9.E-32
      H55M32   = 5.5E-32
      H45M32   = 4.5E-32
      H4M33    = 4.E-33
      H62M34   = 6.2E-34
!
!******FOLLOWING ARE NEGATIVE FLOATING POINT CONSTANTS (HM'S)
!          ARRANGED IN DESCENDING ORDER
      HM2M2    = -.02
      HM6666M2 = -.066667
      HMP5     = -0.5
      HMP575   = -0.575
      HMP66667 = -.66667
      HMP805   = -0.805
      HM1EZ    = -1.
      HM13EZ   = -1.3
      HM19EZ   = -1.9
      HM1E1    = -10.
      HM1597E1 = -15.97469413
      HM161E1  = -16.1
      HM1797E1 = -17.97469413
      HM181E1  = -18.1
      HM8E1    = -80.
      HM1E2    = -100.
!
      RETURN
      END
      SUBROUTINE LWR88(HEATRA,GRNFLX,TOPFLX,
     &                 PRESS,TEMP,RH2O,QO3,CLDFAC,
     &                 CAMT,NCLDS,KTOP,KBTM
!
     &,                L, LP1, LP1V, LLP1, IMAX
!
     &,                STEMP, GTEMP
     &,                CDTM51, CO2M51, C2DM51, CDTM58, CO2M58, C2DM58
     &,                CDT51,  CO251,  C2D51,  CDT58,  CO258,  C2D58
     &,                CDT31,  CO231,  C2D31,  CDT38,  CO238,  C2D38
     &,                CDT71,  CO271,  C2D71,  CDT78,  CO278,  C2D78
!
!    &,                CO211,CO218   ! Not used!!!
!
     &,                SOURCE,DSRCE)
!
CFPP$ NOCONCUR R
!     SUBROUTINE LWR88 COMPUTES TEMPERATURE-CORRECTED CO2 TRANSMISSION
!   FUNCTIONS AND ALSO COMPUTES THE PRESSURE GRID AND LAYER OPTICAL
!   PATHS.
!          INPUTS:                (MODULE BLOCKS)
!      CLDFAC                          CLDCOM
!      PRESS,TEMP,RH2O,QO3             RADISW
!      CAMT,NCLDS,KTOP,KBTM            RADISW
!      CO251,CO258,CDT51,CDT58         CO2BD3
!      C2D51,C2D58,CO2M51,CO2M58       CO2BD3
!      CDTM51,CDTM58,C2DM51,C2DM58     CO2BD3
!      STEMP,GTEMP                     CO2BD3
!      CO231,CO238,CDT31,CDT38         CO2BD2
!      C2D31,C2D38                     CO2BD2
!      CO271,CO278,CDT71,CDT78         CO2BD4
!      C2D71,C2D78                     CO2BD4
!      BETINW                          BDWIDE
!          OUTPUTS:
!      HEATRA,GRNFLX,TOPFLX            LWOUT
!          CALLED BY:
!      RADMN OR INPUT ROUTINE OF MODEL
!          CALLS:
!      FST88
!
c$$$      USE MACHINE , ONLY : kind_rad
      USE HCON
      USE RNDDTA
      implicit none
!
!   B0,B1,B2,B3 ARE COEFFICIENTS USED TO CORRECT FOR THE USE OF 250K IN
!   THE PLANCK FUNCTION USED IN EVALUATING PLANCK-WEIGHTED CO2
!   TRANSMISSION FUNCTIONS. (SEE REF. 4)
!
      real (kind=kind_rad) B0, B1, B2, B3
      PARAMETER (B0=-.51926410E-4, B1=-.18113332E-3,
     &           B2=-.10680132E-5, B3=-.67303519E-7)
!
      integer L, LP1, LP1V, LLP1, IMAX
      integer NCLDS(IMAX), KTOP(IMAX,LP1), KBTM(IMAX,LP1)
!
      real (kind=kind_rad) CO251(LP1,LP1), CO258(LP1,LP1)
     &,                     CDT51(LP1,LP1), CDT58(LP1,LP1)
     &,                     C2D51(LP1,LP1), C2D58(LP1,LP1)
     &,                     CO2M51(L),      CO2M58(L)
     &,                     CDTM51(L),      CDTM58(L)
     &,                     C2DM51(L),      C2DM58(L)
     &,                     STEMP(LP1),     GTEMP(LP1)
!
      real (kind=kind_rad) CO231(LP1), CO238(LP1), CDT31(LP1)
     &,                     CDT38(LP1), C2D31(LP1), C2D38(LP1)
!
      real (kind=kind_rad) CO271(LP1), CO278(LP1), CDT71(LP1)
     &,                     CDT78(LP1), C2D71(LP1), C2D78(LP1)
!
      real (kind=kind_rad) SOURCE(28,NBLY), DSRCE(28,NBLY)
!
      real (kind=kind_rad) PRESS(IMAX,LP1),      TEMP(IMAX,LP1)
     &,                     RH2O(IMAX,L),         QO3(IMAX,L)
     &,                     CLDFAC(IMAX,LP1,LP1), CAMT(IMAX,LP1)
     &,                     HEATRA(IMAX,L),       GRNFLX(IMAX)
     &,                     TOPFLX(IMAX),         DELP2(IMAX,L)
!
      real (kind=kind_rad) QH2O(IMAX,L),       T(IMAX,LP1)
     &,                     P(IMAX,LP1),        DELP(IMAX,L)
     &,                     CO21(IMAX,LP1,LP1), CO2NBL(IMAX,L)
     &,                     CO2SP1(IMAX,LP1),   CO2SP2(IMAX,LP1)
     &,                     VAR1(IMAX,L),       VAR2(IMAX,L)
     &,                     VAR3(IMAX,L),       VAR4(IMAX,L)
     &,                     TOTO3(IMAX,LP1),    TPHIO3(IMAX,LP1)
     &,                     TOTPHI(IMAX,LP1),   TOTVO2(IMAX,LP1)
     &,                     EMX1(IMAX),         EMX2(IMAX)
     &,                     EMPL(IMAX,LLP1),    CNTVAL(IMAX,LP1)
!
      real (kind=kind_rad) CO2R1(IMAX,LP1),    DCO2D1(IMAX,LP1)
     &,                     D2CD21(IMAX,LP1),   D2CD22(IMAX,LP1)
     &,                     CO2R2(IMAX,LP1),    DCO2D2(IMAX,LP1)
     &,                     CO2MR(IMAX,L),      CO2MD(IMAX,L)
     &,                     CO2M2D(IMAX,L),     TDAV(IMAX,LP1)
     &,                     TSTDAV(IMAX,LP1)  , VV(IMAX,L)
     &,                     VSUM3(IMAX,LP1),    DIFT(IMAX,LP1)
     &,                     A1(IMAX),           A2(IMAX)
     &,                     TLSQU(IMAX,LP1)
!
!
      real (kind=kind_rad) texpsl, tem, vsum2, CO2R, DCO2DT, D2CDT2
      integer LL, LM1, LP2, K, I, KP, K1, KK
!
      LL  = LLP1 - 1
      LM1 = L - 1
      LP2 = L + 2
!
!****COMPUTE FLUX PRESSURES (P) AND DIFFERENCES (DELP2,DELP)
!****COMPUTE FLUX LEVEL TEMPERATURES (T) AND CONTINUUM TEMPERATURE
!    CORRECTIONS (TEXPSL)
!
      DO K=2,L
        DO I=1,IMAX
          P(I,K) = HAF*(PRESS(I,K-1)+PRESS(I,K))
          T(I,K) = HAF*(TEMP(I,K-1)+TEMP(I,K))
        ENDDO
      ENDDO
      DO I=1,IMAX
        P(I,1)   = ZERO
        P(I,LP1) = PRESS(I,LP1)
        T(I,1)   = TEMP(I,1)
        T(I,LP1) = TEMP(I,LP1)
      ENDDO
      DO K=1,L
        DO I=1,IMAX
          DELP2(I,K) = P(I,K+1) - P(I,K)
          DELP(I,K)  = ONE / DELP2(I,K)
        ENDDO
      ENDDO
!****COMPUTE ARGUMENT FOR CONT.TEMP.COEFF.
!    (THIS IS 1800.(1./TEMP-1./296.))..THEN TAKE EXPONENTIAL
!     DO  I=1,IMAX*LP1
!       TEXPSL(I,1) = EXP(H18E3/TEMP(I,1)-H6P08108)
!     ENDDO
!***COMPUTE OPTICAL PATHS FOR H2O AND O3, USING THE DIFFUSIVITY
!   APPROXIMATION FOR THE ANGULAR INTEGRATION (1.66). OBTAIN THE
!   UNWEIGHTED VALUES(VAR1,VAR3) AND THE WEIGHTED VALUES(VAR2,VAR4).
!   THE QUANTITIES H3M4(.0003) AND H3M3(.003) APPEARING IN THE VAR2 AND
!   VAR4 EXPRESSIONS ARE THE APPROXIMATE VOIGT CORRECTIONS FOR H2O AND
!   O3,RESPECTIVELY.
!
      DO K=1,L
        DO I=1,IMAX
          QH2O(I,K) = RH2O(I,K)*DIFFCTR
!
!---VV IS THE LAYER-MEAN PRESSURE (IN ATM),WHICH IS NOT THE SAME AS
!   THE LEVEL PRESSURE (PRESS)
!
          VV(I,K)   = HAF*(P(I,K+1)+P(I,K))*P0INV
          VAR1(I,K) = DELP2(I,K) * QH2O(I,K)*GINV
          VAR3(I,K) = DELP2(I,K) * QO3(I,K)*DIFFCTR*GINV
          VAR2(I,K) = VAR1(I,K)  * (VV(I,K)+H3M4)
          VAR4(I,K) = VAR3(I,K)  * (VV(I,K)+H3M3)
!
!  COMPUTE OPTICAL PATH FOR THE H2O CONTINUUM, USING ROBERTS COEFFS.
!  (BETINW),AND TEMP. CORRECTION (TEXPSL). THE DIFFUSIVITY FACTOR
!  (WHICH CANCELS OUT IN THIS EXPRESSION) IS ASSUMED TO BE 1.66. THE
!  USE OF THE DIFFUSIVITY FACTOR HAS BEEN SHOWN TO BE A SIGNIFICANT
!  SOURCE OF ERROR IN THE CONTINUUM CALCS.,BUT THE TIME PENALTY OF
!  AN ANGULAR INTEGRATION IS SEVERE.
!
          TEXPSL      = EXP(H18E3/TEMP(I,K)-H6P08108)
          CNTVAL(I,K) = TEXPSL*RH2O(I,K)*VAR2(I,K)*BETINW/
     &                 (RH2O(I,K)+RATH2OMW)
        ENDDO
      ENDDO
!   COMPUTE SUMMED OPTICAL PATHS FOR H2O,O3 AND CONTINUUM
      DO I=1,IMAX
        TOTPHI(I,1) = ZERO
        TOTO3(I,1)  = ZERO
        TPHIO3(I,1) = ZERO
        TOTVO2(I,1) = ZERO
      ENDDO
      DO K=2,LP1
        DO I=1,IMAX
          TOTPHI(I,K) = TOTPHI(I,K-1) + VAR2(I,K-1)
          TOTO3(I,K)  = TOTO3(I,K-1)  + VAR3(I,K-1)
          TPHIO3(I,K) = TPHIO3(I,K-1) + VAR4(I,K-1)
          TOTVO2(I,K) = TOTVO2(I,K-1) + CNTVAL(I,K-1)
        ENDDO
      ENDDO
!
!---EMX1 IS THE ADDITIONAL PRESSURE-SCALED MASS FROM PRESS(L) TO
!   P(L). IT IS USED IN NEARBY LAYER AND EMISS CALCULATIONS.
!---EMX2 IS THE ADDITIONAL PRESSURE-SCALED MASS FROM PRESS(L) TO
!   P(LP1). IT IS USED IN CALCULATIONS BETWEEN FLUX LEVELS L AND LP1.
!
      DO I=1,IMAX
        TEM     = QH2O(I,L)*PRESS(I,L)*GP0INV
        EMX1(I) = TEM * (PRESS(I,L)-P(I,L))
        EMX2(I) = TEM * (P(I,LP1)-PRESS(I,L))
      ENDDO
!---EMPL IS THE PRESSURE SCALED MASS FROM P(K) TO PRESS(K) (INDEX 2-LP1)
!   OR TO PRESS(K+1) (INDEX LP2-LL)
      DO K=1,L
        DO I=1,IMAX
          EMPL(I,K+1)=QH2O(I,K)*P(I,K+1)*(P(I,K+1)-PRESS(I,K))*GP0INV
        ENDDO
      ENDDO
      DO K=1,LM1
        KK = K + LP1
        K1 = K + 1
        DO I=1,IMAX
          EMPL(I,KK)=QH2O(I,K1)*P(I,K1)*(PRESS(I,K1)-P(I,K1))*GP0INV
        ENDDO
      ENDDO
      DO I=1,IMAX
        EMPL(I,1)    = VAR2(I,L)
        EMPL(I,LLP1) = EMPL(I,LL)
      ENDDO
!***COMPUTE WEIGHTED TEMPERATURE (TDAV) AND PRESSURE (TSTDAV) INTEGRALS
!   FOR USE IN OBTAINING TEMP. DIFFERENCE BET. SOUNDING AND STD.
!   TEMP. SOUNDING (DIFT)
      DO I=1,IMAX
        TSTDAV(I,1) = ZERO
        TDAV(I,1)   = ZERO
      ENDDO
      DO K=1,LP1
        DO I=1,IMAX
          VSUM3(I,K) = TEMP(I,K) - STEMP(K)
        ENDDO
      ENDDO
      DO K=1,L
        DO I=1,IMAX
          VSUM2         = GTEMP(K)    * DELP2(I,K)
          TSTDAV(I,K+1) = TSTDAV(I,K) + VSUM2
          TDAV(I,K+1)   = TDAV(I,K)   + VSUM2 * VSUM3(I,K)
        ENDDO
      ENDDO
!
!****EVALUATE COEFFICIENTS FOR CO2 PRESSURE INTERPOLATION (A1,A2)
      TEM = 1.0 / P0XZP2
      DO I=1,IMAX
        A1(I) = (PRESS(I,LP1)-P0XZP8)*TEM
        A2(I) = (P0-PRESS(I,LP1))*TEM
      ENDDO
!
!***PERFORM CO2 PRESSURE INTERPOLATION ON ALL INPUTTED TRANSMISSION
!   FUNCTIONS AND TEMP. DERIVATIVES
!---SUCCESSIVELY COMPUTING CO2R,DCO2DT AND D2CDT2 IS DONE TO SAVE
!   STORAGE (AT A SLIGHT LOSS IN COMPUTATION TIME)
      DO K=1,LP1
        DO I=1,IMAX
          CO2R1(I,K)  =       A1(I)*CO231(K)+A2(I)*CO238(K)
          D2CD21(I,K) = H1M3*(A1(I)*C2D31(K)+A2(I)*C2D38(K))
          DCO2D1(I,K) = H1M2*(A1(I)*CDT31(K)+A2(I)*CDT38(K))
          CO2R2(I,K)  =       A1(I)*CO271(K)+A2(I)*CO278(K)
          D2CD22(I,K) = H1M3*(A1(I)*C2D71(K)+A2(I)*C2D78(K))
          DCO2D2(I,K) = H1M2*(A1(I)*CDT71(K)+A2(I)*CDT78(K))
        ENDDO
      ENDDO
      DO K=1,L
        DO I=1,IMAX
          CO2MR(I,K)  =       A1(I)*CO2M51(K)+A2(I)*CO2M58(K)
          CO2MD(I,K)  = H1M2*(A1(I)*CDTM51(K)+A2(I)*CDTM58(K))
          CO2M2D(I,K) = H1M3*(A1(I)*C2DM51(K)+A2(I)*C2DM58(K))
        ENDDO
      ENDDO
!
!***COMPUTE CO2 TEMPERATURE INTERPOLATIONS FOR ALL BANDS,USING DIFT
!
!   THE CASE WHERE K=1 IS HANDLED FIRST. WE ARE NOW REPLACING
!   3-DIMENSIONAL ARRAYS BY 2-D ARRAYS, TO SAVE SPACE. THUS THIS
!   CALCULATION IS FOR (I,KP,1)
!
      DO KP=2,LP1
        DO I=1,IMAX
          DIFT(I,KP) = TDAV(I,KP) / TSTDAV(I,KP)
        ENDDO
      ENDDO
      DO I=1,IMAX
        CO21(I,1,1) = 1.0
        CO2SP1(I,1) = 1.0
        CO2SP2(I,1) = 1.0
      ENDDO
      DO KP=2,LP1
        DO I=1,IMAX
!---CALCULATIONS FOR KP>1 FOR K=1
          CO2R         =       A1(I)*CO251(KP,1)+A2(I)*CO258(KP,1)
          DCO2DT       = H1M2*(A1(I)*CDT51(KP,1)+A2(I)*CDT58(KP,1))
          D2CDT2       = H1M3*(A1(I)*C2D51(KP,1)+A2(I)*C2D58(KP,1))
          CO21(I,KP,1) = CO2R+DIFT(I,KP)*(DCO2DT+HAF*DIFT(I,KP)*D2CDT2)
!---CALCULATIONS FOR (EFFECTIVELY) KP=1,K>KP. THESE USE THE
!   SAME VALUE OF DIFT DUE TO SYMMETRY
          CO2R         =       A1(I)*CO251(1,KP)+A2(I)*CO258(1,KP)
          DCO2DT       = H1M2*(A1(I)*CDT51(1,KP)+A2(I)*CDT58(1,KP))
          D2CDT2       = H1M3*(A1(I)*C2D51(1,KP)+A2(I)*C2D58(1,KP))
          CO21(I,1,KP) = CO2R+DIFT(I,KP)*(DCO2DT+HAF*DIFT(I,KP)*D2CDT2)
        ENDDO
      ENDDO
!
!   THE TRANSMISSION FUNCTIONS USED IN SPA88 MAY BE COMPUTED NOW.
!---(IN THE 250 LOOP,DIFT REALLY SHOULD BE (I,1,K), BUT DIFT IS
!    INVARIANT WITH RESPECT TO K,KP,AND SO (I,1,K)=(I,K,1))
      DO K=2,LP1
        DO I=1,IMAX
          CO2SP1(I,K) = CO2R1(I,K)+DIFT(I,K)*(DCO2D1(I,K)+HAF*DIFT(I,K)*
     &                                                     D2CD21(I,K))
          CO2SP2(I,K) = CO2R2(I,K)+DIFT(I,K)*(DCO2D2(I,K)+HAF*DIFT(I,K)*
     &                                                     D2CD22(I,K))
        ENDDO
      ENDDO
!
!   NEXT THE CASE WHEN K=2...L
      DO K=2,L
        DO KP=K+1,LP1
          DO I=1,IMAX
           DIFT(I,KP)   = (TDAV(I,KP)-TDAV(I,K))/
     &                    (TSTDAV(I,KP)-TSTDAV(I,K))
!
           CO2R         =       A1(I)*CO251(KP,K)+A2(I)*CO258(KP,K)
           DCO2DT       = H1M2*(A1(I)*CDT51(KP,K)+A2(I)*CDT58(KP,K))
           D2CDT2       = H1M3*(A1(I)*C2D51(KP,K)+A2(I)*C2D58(KP,K))
           CO21(I,KP,K) = CO2R+DIFT(I,KP)*(DCO2DT+HAF*DIFT(I,KP)*D2CDT2)
!
           CO2R         =       A1(I)*CO251(K,KP)+A2(I)*CO258(K,KP)
           DCO2DT       = H1M2*(A1(I)*CDT51(K,KP)+A2(I)*CDT58(K,KP))
           D2CDT2       = H1M3*(A1(I)*C2D51(K,KP)+A2(I)*C2D58(K,KP))
           CO21(I,K,KP) = CO2R+DIFT(I,KP)*(DCO2DT+HAF*DIFT(I,KP)*D2CDT2)
          ENDDO
        ENDDO
      ENDDO
!   FINALLY THE CASE WHEN K=KP,K=2..LP1
      DO K=2,LP1
        DO I=1,IMAX
          DIFT(I,K)   = HAF*(VSUM3(I,K)+VSUM3(I,K-1))
          CO2R        =       A1(I)*CO251(K,K)+A2(I)*CO258(K,K)
          DCO2DT      = H1M2*(A1(I)*CDT51(K,K)+A2(I)*CDT58(K,K))
          D2CDT2      = H1M3*(A1(I)*C2D51(K,K)+A2(I)*C2D58(K,K))
          CO21(I,K,K) = CO2R+DIFT(I,K)*(DCO2DT+HAF*DIFT(I,K)*D2CDT2)
        ENDDO
      ENDDO
!--- WE AREN'T DOING NBL TFS ON THE 100 CM-1 BANDS .
      DO K=1,L
        DO I=1,IMAX
          CO2NBL(I,K) = CO2MR(I,K) + VSUM3(I,K) *
     &                  (CO2MD(I,K)+HAF*VSUM3(I,K)*CO2M2D(I,K))
        ENDDO
      ENDDO
!***COMPUTE TEMP. COEFFICIENT BASED ON T(K) (SEE REF.2)
      DO K=1,LP1
        DO I=1,IMAX
          IF (T(I,K).LE.H25E2) THEN
            TEM        = T(I,K) - H25E2
            TLSQU(I,K) = B0 + TEM * (B1 + TEM * (B2 + B3*TEM))
          ELSE
            TLSQU(I,K) = B0
          ENDIF
        ENDDO
      ENDDO
!***APPLY TO ALL CO2 TFS
      DO K=1,LP1
        DO KP=1,LP1
          DO I=1,IMAX
            CO21(I,KP,K) = CO21(I,KP,K)*(ONE-TLSQU(I,KP)) + TLSQU(I,KP)
          ENDDO
        ENDDO
        DO I=1,IMAX
          CO2SP1(I,K) = CO2SP1(I,K)*(ONE-TLSQU(I,1)) + TLSQU(I,1)
          CO2SP2(I,K) = CO2SP2(I,K)*(ONE-TLSQU(I,1)) + TLSQU(I,1)
        ENDDO
      ENDDO
      DO K=1,L
        DO I=1,IMAX
          CO2NBL(I,K) = CO2NBL(I,K)*(ONE-TLSQU(I,K))+TLSQU(I,K)
        ENDDO
      ENDDO
      CALL FST88(HEATRA,GRNFLX,TOPFLX,
     &           QH2O,PRESS,P,DELP,DELP2,TEMP,T,
     &           CLDFAC,
!    &           CLDFAC,NCLDS,KTOP,KBTM,CAMT,
     &           CO21,CO2NBL,CO2SP1,CO2SP2,
     &           VAR1,VAR2,VAR3,VAR4,CNTVAL,
     &           TOTO3,TPHIO3,TOTPHI,TOTVO2,
     &           EMX1,EMX2,EMPL
     &,          L, LP1, LP1V, LLP1, IMAX
     &,          SOURCE,DSRCE)
      RETURN
      END
!     SUBROUTINE SPA88 COMPUTES EXACT CTS HEATING RATES AND FLUXES AND
!  CORRESPONDING CTS EMISSIVITY QUANTITIES FOR H2O,CO2 AND O3.
!          INPUTS:                (MODULE BLOCKS)
!       ACOMB,BCOMB,APCM,BPCM                  BDCOMB
!       ATPCM,BTPCM,BETACM                     BDCOMB
!       BETINW                                 BDWIDE
!       TEMP,PRESS                             RADISW
!       VAR1,VAR2,P,DELP,DELP2                 KDACOM
!       TOTVO2,TO3SP,TO3SPC                    TFCOM
!       CO2SP1,CO2SP2,CO2SP                    TFCOM
!       CLDFAC                                 CLDCOM
!       SKO2D                                  TABCOM
!       SORC,CSOUR                             SRCCOM
!           OUTPUTS:
!       EXCTS,CTSO3                            TFCOM
!       GXCTS                                  RDFLUX
!           CALLED BY:
!       FST88
!            CALLS:
!
      SUBROUTINE SPA88(EXCTS,CTSO3,GXCTS,SORC,CSOUR,
     &                 CLDFAC,TEMP,PRESS,VAR1,VAR2,
     &                 P,DELP,DELP2,TOTVO2,TO3SP,TO3SPC,
     &                 CO2SP1,CO2SP2,CO2SP
     &,                L, LP1, IMAX)
CFPP$ NOCONCUR R
!
c$$$      USE MACHINE , ONLY : kind_rad
      USE HCON
      USE RNDDTA
      implicit none
!
      integer L, LP1, IMAX
!
      real (kind=kind_rad) SORC(IMAX,LP1,NBLY), CSOUR(IMAX,LP1)
     &,                     CLDFAC(IMAX,LP1,LP1)
     &,                     TEMP(IMAX,LP1),      PRESS(IMAX,LP1)
     &,                     VAR1(IMAX,L),        VAR2(IMAX,L)
     &,                     P(IMAX,LP1),         DELP(IMAX,L)
     &,                     DELP2(IMAX,L),       TOTVO2(IMAX,LP1)
     &,                     TO3SPC(IMAX,L),      TO3SP(IMAX,LP1)
     &,                     CO2SP1(IMAX,LP1),    CO2SP2(IMAX,LP1)
     &,                     CO2SP(IMAX,LP1),     EXCTS(IMAX,L)
     &,                     CTSO3(IMAX,L),       GXCTS(IMAX)
!
      real (kind=kind_rad) PHITMP(IMAX,L),      PSITMP(IMAX,L)
     &,                     TT(IMAX,L),          CTMP(IMAX,LP1)
     &,                     X(IMAX,L),           Y(IMAX,L)
     &,                     TOPM(IMAX,L),        TOPPHI(IMAX,L)
     &,                     CTMP3(IMAX,LP1),     CTMP2(IMAX,LP1)
!
      real (kind=kind_rad) F, FF, AG, AGG, FAC1, FAC2, TEM
      integer lm1, k, i, ib
!
      LM1 = L - 1
!
!---COMPUTE TEMPERATURE QUANTITIES FOR USE IN PROGRAM
      DO K=1,L
        DO I=1,IMAX
          X(I,K) = TEMP(I,K) - H25E2
          Y(I,K) = X(I,K) * X(I,K)
!
!     Initialize some arrays
!
          EXCTS(I,K) =  0.0
        ENDDO
      ENDDO
!---INITIALIZE CTMP(I,1),CTMP2(I,1),CTMP3(I,1) TO UNITY; THESE ARE
!   TRANSMISSION FCTNS AT THE TOP.
      DO I=1,IMAX
        CTMP(I,1)  = ONE
        CTMP2(I,1) = 1.
        CTMP3(I,1) = 1.
        GXCTS(I)   = 0.0
 
!       For Clear Sky
      ENDDO
!
!***BEGIN LOOP ON FREQUENCY BANDS ***
!
!-----CALCULATION FOR BAND 1 (COMBINED BAND 1)
!-----CALCULATION FOR BAND 2 (COMBINED BAND 2)
!-----CALCULATION FOR BAND 3 (COMBINED BAND 3)
!-----CALCULATION FOR BAND 4 (COMBINED BAND 4)
!-----CALCULATION FOR BAND 5 (COMBINED BAND 5)
!-----CALCULATION FOR BAND 6 (COMBINED BAND 6)
!-----CALCULATION FOR BAND 7 (COMBINED BAND 7)
!-----CALCULATION FOR BAND 8 (COMBINED BAND 8)
!-----CALCULATION FOR BAND 9 ( 560-670 CM-1; INCLUDES CO2)
!-----CALCULATION FOR BAND 10 (670-800 CM-1; INCLUDES CO2)
!-----CALCULATION FOR BAND 11 (800-900 CM-1)
!-----CALCULATION FOR BAND 12 (900-990 CM-1)
!-----CALCULATION FOR BAND 13 (990-1070 CM-1; INCLUDES O3))
!-----CALCULATION FOR BAND 14 (1070-1200 CM-1)
!
      DO IB=1,14
!
!
!---    CALCULATION FOR SINGLE BAND (COMBINED BAND)
!
!---    OBTAIN TEMPERATURE CORRECTION (CAPPHI,CAPPSI),THEN MULTIPLY
!       BY OPTICAL PATH (VAR1,VAR2) TO COMPUTE TEMPERATURE-CORRECTED
!       OPTICAL PATH AND MEAN PRESSURE FOR A LAYER (PHITMP,PSITMP)
!
        DO K=1,L
          DO I=1,IMAX
            F           = H44194M2*(APCM(IB)*X(I,K)+BPCM(IB)*Y(I,K))
            FF          = H44194M2*(ATPCM(IB)*X(I,K)+BTPCM(IB)*Y(I,K))
            AG          = (H1P41819+F)*F   + ONE
            AGG         = (H1P41819+FF)*FF + ONE
!
            AG          = AG * AG      !  AG ** 2
            AG          = AG * AG      !  AG ** 4
            AG          = AG * AG      !  AG ** 8
            AGG         = AGG * AGG
            AGG         = AGG * AGG
            AGG         = AGG * AGG
!
            PHITMP(I,K) = VAR1(I,K) * (AG*AG)  ! AG ** 16
            PSITMP(I,K) = VAR2(I,K) * (AGG*AGG)
          ENDDO
        ENDDO
!---    OBTAIN OPTICAL PATH,MEAN PRESSURE FROM THE TOP TO THE PRESSURE
!       P(K) (TOPM,TOPPHI)
        DO I=1,IMAX
          TOPM(I,1)   = PHITMP(I,1)
          TOPPHI(I,1) = PSITMP(I,1)
        ENDDO
        DO K=2,L
          DO I=1,IMAX
            TOPM(I,K)   = TOPM(I,K-1)   + PHITMP(I,K)
            TOPPHI(I,K) = TOPPHI(I,K-1) + PSITMP(I,K)
          ENDDO
        ENDDO
!---    TT IS THE CLOUD-FREE CTS TRANSMISSION FUNCTION
        IF (IB .LT. 5) THEN
          DO K=1,L
            DO I=1,IMAX
              FAC1      = ACOMB(IB)*TOPM(I,K)
              FAC2      = FAC1*TOPM(I,K)/(BCOMB(IB)*TOPPHI(I,K))
              TT(I,K)   = EXP(HM1EZ*FAC1/SQRT(1.+FAC2))
            ENDDO
          ENDDO
        ELSEIF (IB .LT. 9 .OR. (IB .GT. 10 .AND. IB .NE. 13)) THEN
          DO K=1,L
            DO I=1,IMAX
              FAC1      = ACOMB(IB)*TOPM(I,K)
              FAC2      = FAC1*TOPM(I,K)/(BCOMB(IB)*TOPPHI(I,K))
              TT(I,K)   = EXP(HM1EZ*(FAC1/SQRT(1.+FAC2)+
     &                    BETACM(IB)*TOTVO2(I,K+1)*SKO2D))
          ENDDO
          ENDDO
        ELSEIF (IB .EQ. 9) THEN
          DO K=1,L
            DO I=1,IMAX
              FAC1      = ACOMB(IB)*TOPM(I,K)
              FAC2      = FAC1*TOPM(I,K)/(BCOMB(IB)*TOPPHI(I,K))
              TT(I,K)   = EXP(HM1EZ*(FAC1/SQRT(1.+FAC2)+
     &                    BETACM(IB)*TOTVO2(I,K+1)*SKO2D))*CO2SP1(I,K+1)
            ENDDO
          ENDDO
        ELSEIF (IB .EQ. 10) THEN
          DO K=1,L
            DO I=1,IMAX
              FAC1      = ACOMB(IB)*TOPM(I,K)
              FAC2      = FAC1*TOPM(I,K)/(BCOMB(IB)*TOPPHI(I,K))
              TT(I,K)   = EXP(HM1EZ*(FAC1/SQRT(1.+FAC2)+
     &                    BETACM(IB)*TOTVO2(I,K+1)*SKO2D))*CO2SP2(I,K+1)
            ENDDO
          ENDDO
        ELSEIF (IB .EQ. 13) THEN
          DO K=1,L
            DO I=1,IMAX
              FAC1      = ACOMB(IB)*TOPM(I,K)
              FAC2      = FAC1*TOPM(I,K)/(BCOMB(IB)*TOPPHI(I,K))
              TT(I,K)   = EXP(HM1EZ*(FAC1/SQRT(1.+FAC2)+
     &                    BETACM(IB)*TOTVO2(I,K+1)*SKO2D +TO3SPC(I,K)))
            ENDDO
          ENDDO
        ENDIF
        DO K=1,L
          DO I=1,IMAX
            CTMP(I,K+1)  = TT(I,K)*CLDFAC(I,K+1,1)
          ENDDO
        ENDDO
!
!---    EXCTS IS THE CTS COOLING RATE ACCUMULATED OVER FREQUENCY BANDS
        DO K=1,L
          DO I=1,IMAX
            EXCTS(I,K)  = EXCTS(I,K)
     &                  + SORC(I,K,IB) * (CTMP(I,K+1)-CTMP(I,K))
          ENDDO
        ENDDO
!---    GXCTS IS THE EXACT CTS TOP FLUX ACCUMULATED OVER FREQUENCY BANDS
        DO I=1,IMAX
          TEM = TT(I,L)*SORC(I,L,IB)+
     &       (HAF*DELP(I,L)*(TT(I,LM1)*(P(I,LP1)-PRESS(I,L)) +
     &       TT(I,L)*(P(I,LP1)+PRESS(I,L)-TWO*P(I,L)))) *
     &       (SORC(I,LP1,IB)-SORC(I,L,IB))
          GXCTS(I)  = GXCTS(I)  + TEM * CLDFAC(I,LP1,1)
        ENDDO
!
      ENDDO                         ! Band Loop Ends here!
!
!
!   OBTAIN CTS FLUX AT THE TOP BY INTEGRATION OF HEATING RATES AND
!   USING CTS FLUX AT THE BOTTOM (CURRENT VALUE OF GXCTS). NOTE
!   THAT THE PRESSURE QUANTITIES AND CONVERSION FACTORS HAVE NOT
!   BEEN INCLUDED EITHER IN EXCTS OR IN GXCTS. THESE CANCEL OUT, THUS
!   REDUCING COMPUTATIONS!
!
      DO K=1,L
        DO I=1,IMAX
          GXCTS(I)  = GXCTS(I)  - EXCTS(I,K)
        ENDDO
      ENDDO
!
!   NOW SCALE THE COOLING RATE (EXCTS) BY INCLUDING THE PRESSURE
!   FACTOR (DELP) AND THE CONVERSION FACTOR (RADCON)
!
!     DO I=1,IMAX*L
!       EXCTS(I,1)  = EXCTS(I,1) *RADCON*DELP(I,1)
!     ENDDO
!---THIS IS THE END OF THE EXACT CTS COMPUTATIONS; AT THIS POINT
!   EXCTS HAS ITS APPROPRIATE VALUE.
!
!*** COMPUTE APPROXIMATE CTS HEATING RATES FOR 15UM AND 9.6 UM BANDS
!     (CTSO3)
      DO K=1,L
        DO I=1,IMAX
          CTMP2(I,K+1)  = CO2SP(I,K+1) * CLDFAC(I,K+1,1)
          CTMP3(I,K+1)  = TO3SP(I,K) * CLDFAC(I,K+1,1)
        ENDDO
      ENDDO
      DO K=1,L
        DO I=1,IMAX
!         CTSO3(I,K) = RADCON*DELP(I,K)*
          CTSO3(I,K) = CSOUR(I,K)*(CTMP2(I,K+1)-CTMP2(I,K)) +
     &                 SORC(I,K,13)*(CTMP3(I,K+1)-CTMP3(I,K))
!
        ENDDO
      ENDDO
!
      RETURN
      END
      SUBROUTINE LWTABLE(LP1,LP1V, SOURCE,DSRCE)
CFPP$ NOCONCUR R
!     SUBROUTINE TABLE COMPUTES TABLE ENTRIES USED IN THE LONGWAVE RADIA
!     PROGRAM. ALSO CALCULATED ARE INDICES USED IN STRIP-MINING AND FOR
!     SOME PRE-COMPUTABLE FUNCTIONS.
!         INPUTS:
!         OUTPUTS:
!       EM1V,EM1VW,T1,T2,T4                     TABCOM
!       EM3,SOURCE,DSRCE,IND,INDX2,KMAXV        TABCOM
!       KMAXVM,                                 TABCOM
!       AO3RND,BO3RND,AB15                      BANDTA
!       AB15WD,SKC1R,SKO3R,SKO2D                BDWIDE
!
c$$$      USE MACHINE , ONLY : kind_rad
      USE HCON
      USE RNDDTA
      implicit none
!
      integer lp1, lp1v,i1
!
      real (kind=kind_rad) SOURCE(28,NBLY), DSRCE(28,NBLY)
!
      real (kind=kind_rad) SUM(28,180),     PERTSM(28,180)
     &,                     SUM3(28,180),    SUMWDE(28,180)
     &,                     SRCWD(28,NBLX),  SRC1NB(28,NBLW)
     &,                     DBDTNB(28,NBLW), ZMASS(181)
     &,                     ZROOT(181),      SC(28),     DSC(28)
     &,                     XTEMV(28),       TFOUR(28),  FORTCU(28)
     &,                     X(28),           X1(28),     X2(180)
     &,                     SRCS(28),        SUM4(28),   SUM6(28)
     &,                     SUM7(28),        SUM8(28),   SUM4WD(28)
     &,                     R1(28),          R2(28),     S2(28)
     &,                     T3(28),          R1WD(28)
     &,                     EXPO(180),       FAC(180)
     &,                     CNUSB(30),       DNUSB(30)
     &,                     ALFANB(NBLW),    AROTNB(NBLW)
     &,                     ANB(NBLW),       BNB(NBLW),  CENTNB(NBLW)
     &,                     DELNB(NBLW),     BETANB(NBLW)
!
      real (kind=kind_rad) CENT, DEL, BDHI, BDLO, ANU, C1
      integer L, LP2, N, J, JP, I, IA, NSUBDS, NSB
      real (kind=kind_rad) ARNDM1(64), ARNDM2(64), ARNDM3(35)
     &,                     BRNDM1(64), BRNDM2(64), BRNDM3(35)
     &,                     AP1(64),    AP2(64),    AP3(35)
     &,                     BP1(64),    BP2(64),    BP3(35)
     &,                     ATP1(64),   ATP2(64),   ATP3(35)
     &,                     BTP1(64),   BTP2(64),   BTP3(35)
     &,                     BETAD1(64), BETAD2(64), BETAD3(35)
     &,                     BANDL1(64), BANDL2(64), BANDL3(35)
     &,                     BANDH1(64), BANDH2(64), BANDH3(35)
!
!***THE FOLLOWING DATA STATEMENTS ARE BAND PARAMETERS OBTAINED USING
!   THE 1982 AFGL CATALOG ON THE SPECIFIED BANDS
!
      DATA ARNDM1  /
     &   0.354693E+00,  0.269857E+03,  0.167062E+03,  0.201314E+04,
     &   0.964533E+03,  0.547971E+04,  0.152933E+04,  0.599429E+04,
     &   0.699329E+04,  0.856721E+04,  0.962489E+04,  0.233348E+04,
     &   0.127091E+05,  0.104383E+05,  0.504249E+04,  0.181227E+05,
     &   0.856480E+03,  0.136354E+05,  0.288635E+04,  0.170200E+04,
     &   0.209761E+05,  0.126797E+04,  0.110096E+05,  0.336436E+03,
     &   0.491663E+04,  0.863701E+04,  0.540389E+03,  0.439786E+04,
     &   0.347836E+04,  0.130557E+03,  0.465332E+04,  0.253086E+03,
     &   0.257387E+04,  0.488041E+03,  0.892991E+03,  0.117148E+04,
     &   0.125880E+03,  0.458852E+03,  0.142975E+03,  0.446355E+03,
     &   0.302887E+02,  0.394451E+03,  0.438112E+02,  0.348811E+02,
     &   0.615503E+02,  0.143165E+03,  0.103958E+02,  0.725108E+02,
     &   0.316628E+02,  0.946456E+01,  0.542675E+02,  0.351557E+02,
     &   0.301797E+02,  0.381010E+01,  0.126319E+02,  0.548010E+01,
     &   0.600199E+01,  0.640803E+00,  0.501549E-01,  0.167961E-01,
     &   0.178110E-01,  0.170166E+00,  0.273514E-01,  0.983767E+00/
      DATA ARNDM2  /
     &   0.753946E+00,  0.941763E-01,  0.970547E+00,  0.268862E+00,
     &   0.564373E+01,  0.389794E+01,  0.310955E+01,  0.128235E+01,
     &   0.196414E+01,  0.247113E+02,  0.593435E+01,  0.377552E+02,
     &   0.305173E+02,  0.852479E+01,  0.116780E+03,  0.101490E+03,
     &   0.138939E+03,  0.324228E+03,  0.683729E+02,  0.471304E+03,
     &   0.159684E+03,  0.427101E+03,  0.114716E+03,  0.106190E+04,
     &   0.294607E+03,  0.762948E+03,  0.333199E+03,  0.830645E+03,
     &   0.162512E+04,  0.525676E+03,  0.137739E+04,  0.136252E+04,
     &   0.147164E+04,  0.187196E+04,  0.131118E+04,  0.103975E+04,
     &   0.621637E+01,  0.399459E+02,  0.950648E+02,  0.943161E+03,
     &   0.526821E+03,  0.104150E+04,  0.905610E+03,  0.228142E+04,
     &   0.806270E+03,  0.691845E+03,  0.155237E+04,  0.192241E+04,
     &   0.991871E+03,  0.123907E+04,  0.457289E+02,  0.146146E+04,
     &   0.319382E+03,  0.436074E+03,  0.374214E+03,  0.778217E+03,
     &   0.140227E+03,  0.562540E+03,  0.682685E+02,  0.820292E+02,
     &   0.178779E+03,  0.186150E+03,  0.383864E+03,  0.567416E+01/
      DATA ARNDM3  /
     &   0.225129E+03,  0.473099E+01,  0.753149E+02,  0.233689E+02,
     &   0.339802E+02,  0.108855E+03,  0.380016E+02,  0.151039E+01,
     &   0.660346E+02,  0.370165E+01,  0.234169E+02,  0.440206E+00,
     &   0.615283E+01,  0.304077E+02,  0.117769E+01,  0.125248E+02,
     &   0.142652E+01,  0.241831E+00,  0.483721E+01,  0.226357E-01,
     &   0.549835E+01,  0.597067E+00,  0.404553E+00,  0.143584E+01,
     &   0.294291E+00,  0.466273E+00,  0.156048E+00,  0.656185E+00,
     &   0.172727E+00,  0.118349E+00,  0.141598E+00,  0.588581E-01,
     &   0.919409E-01,  0.155521E-01,  0.537083E-02/
      DATA BRNDM1  /
     &   0.789571E-01,  0.920256E-01,  0.696960E-01,  0.245544E+00,
     &   0.188503E+00,  0.266127E+00,  0.271371E+00,  0.330917E+00,
     &   0.190424E+00,  0.224498E+00,  0.282517E+00,  0.130675E+00,
     &   0.212579E+00,  0.227298E+00,  0.138585E+00,  0.187106E+00,
     &   0.194527E+00,  0.177034E+00,  0.115902E+00,  0.118499E+00,
     &   0.142848E+00,  0.216869E+00,  0.149848E+00,  0.971585E-01,
     &   0.151532E+00,  0.865628E-01,  0.764246E-01,  0.100035E+00,
     &   0.171133E+00,  0.134737E+00,  0.105173E+00,  0.860832E-01,
     &   0.148921E+00,  0.869234E-01,  0.106018E+00,  0.184865E+00,
     &   0.767454E-01,  0.108981E+00,  0.123094E+00,  0.177287E+00,
     &   0.848146E-01,  0.119356E+00,  0.133829E+00,  0.954505E-01,
     &   0.155405E+00,  0.164167E+00,  0.161390E+00,  0.113287E+00,
     &   0.714720E-01,  0.741598E-01,  0.719590E-01,  0.140616E+00,
     &   0.355356E-01,  0.832779E-01,  0.128680E+00,  0.983013E-01,
     &   0.629660E-01,  0.643346E-01,  0.717082E-01,  0.629730E-01,
     &   0.875182E-01,  0.857907E-01,  0.358808E+00,  0.178840E+00/
      DATA BRNDM2  /
     &   0.254265E+00,  0.297901E+00,  0.153916E+00,  0.537774E+00,
     &   0.267906E+00,  0.104254E+00,  0.400723E+00,  0.389670E+00,
     &   0.263701E+00,  0.338116E+00,  0.351528E+00,  0.267764E+00,
     &   0.186419E+00,  0.238237E+00,  0.210408E+00,  0.176869E+00,
     &   0.114715E+00,  0.173299E+00,  0.967770E-01,  0.172565E+00,
     &   0.162085E+00,  0.157782E+00,  0.886832E-01,  0.242999E+00,
     &   0.760298E-01,  0.164248E+00,  0.221428E+00,  0.166799E+00,
     &   0.312514E+00,  0.380600E+00,  0.353828E+00,  0.269500E+00,
     &   0.254759E+00,  0.285408E+00,  0.159764E+00,  0.721058E-01,
     &   0.170528E+00,  0.231595E+00,  0.307184E+00,  0.564136E-01,
     &   0.159884E+00,  0.147907E+00,  0.185666E+00,  0.183567E+00,
     &   0.182482E+00,  0.230650E+00,  0.175348E+00,  0.195978E+00,
     &   0.255323E+00,  0.198517E+00,  0.195500E+00,  0.208356E+00,
     &   0.309603E+00,  0.112011E+00,  0.102570E+00,  0.128276E+00,
     &   0.168100E+00,  0.177836E+00,  0.105533E+00,  0.903330E-01,
     &   0.126036E+00,  0.101430E+00,  0.124546E+00,  0.221406E+00/
      DATA BRNDM3  /
     &   0.137509E+00,  0.911365E-01,  0.724508E-01,  0.795788E-01,
     &   0.137411E+00,  0.549175E-01,  0.787714E-01,  0.165544E+00,
     &   0.136484E+00,  0.146729E+00,  0.820496E-01,  0.846211E-01,
     &   0.785821E-01,  0.122527E+00,  0.125359E+00,  0.101589E+00,
     &   0.155756E+00,  0.189239E+00,  0.999086E-01,  0.480993E+00,
     &   0.100233E+00,  0.153754E+00,  0.130780E+00,  0.136136E+00,
     &   0.159353E+00,  0.156634E+00,  0.272265E+00,  0.186874E+00,
     &   0.192090E+00,  0.135397E+00,  0.131497E+00,  0.127463E+00,
     &   0.227233E+00,  0.190562E+00,  0.214005E+00/
      DATA AP1     /
     &  -0.675950E-02, -0.909459E-02, -0.800214E-02, -0.658673E-02,
     &  -0.245580E-02, -0.710464E-02, -0.205565E-02, -0.446529E-02,
     &  -0.440265E-02, -0.593625E-02, -0.201913E-02, -0.349169E-02,
     &  -0.209324E-02, -0.127980E-02, -0.388007E-02, -0.140542E-02,
     &   0.518346E-02, -0.159375E-02,  0.250508E-02,  0.132182E-01,
     &  -0.903779E-03,  0.110959E-01,  0.924528E-03,  0.207428E-01,
     &   0.364166E-02,  0.365229E-02,  0.884367E-02,  0.617260E-02,
     &   0.701340E-02,  0.184265E-01,  0.992822E-02,  0.908582E-02,
     &   0.106581E-01,  0.276268E-02,  0.158414E-01,  0.145747E-01,
     &   0.453080E-02,  0.214767E-01,  0.553895E-02,  0.195031E-01,
     &   0.237016E-01,  0.112371E-01,  0.275977E-01,  0.188833E-01,
     &   0.131079E-01,  0.130019E-01,  0.385122E-01,  0.111768E-01,
     &   0.622620E-02,  0.194397E-01,  0.134360E-01,  0.207829E-01,
     &   0.147960E-01,  0.744479E-02,  0.107564E-01,  0.181562E-01,
     &   0.170062E-01,  0.233303E-01,  0.256735E-01,  0.274745E-01,
     &   0.279259E-01,  0.197002E-01,  0.140268E-01,  0.185933E-01/
      DATA AP2     /
     &   0.169525E-01,  0.214410E-01,  0.136577E-01,  0.169510E-01,
     &   0.173025E-01,  0.958346E-02,  0.255024E-01,  0.308943E-01,
     &   0.196031E-01,  0.183608E-01,  0.149419E-01,  0.206358E-01,
     &   0.140654E-01,  0.172797E-01,  0.145470E-01,  0.982987E-02,
     &   0.116695E-01,  0.811333E-02,  0.965823E-02,  0.649977E-02,
     &   0.462192E-02,  0.545929E-02,  0.680407E-02,  0.291235E-02,
     &  -0.974773E-03,  0.341591E-02,  0.376198E-02,  0.770610E-03,
     &  -0.940864E-04,  0.514532E-02,  0.232371E-02, -0.177741E-02,
     &  -0.374892E-03, -0.370485E-03, -0.221435E-02, -0.490000E-02,
     &   0.588664E-02,  0.931411E-03, -0.456043E-03, -0.545576E-02,
     &  -0.421136E-02, -0.353742E-02, -0.174276E-02, -0.361246E-02,
     &  -0.337822E-02, -0.867030E-03, -0.118001E-02, -0.222405E-02,
     &  -0.725144E-03,  0.118483E-02,  0.995087E-02,  0.273812E-03,
     &   0.417298E-02,  0.764294E-02,  0.631568E-02, -0.213528E-02,
     &   0.746130E-02,  0.110337E-02,  0.153157E-01,  0.504532E-02,
     &   0.406047E-02,  0.192895E-02,  0.202058E-02,  0.126420E-01/
      DATA AP3     /
     &   0.310028E-02,  0.214779E-01,  0.560165E-02,  0.661070E-02,
     &   0.694966E-02,  0.539194E-02,  0.103745E-01,  0.180150E-01,
     &   0.747133E-02,  0.114927E-01,  0.115213E-01,  0.160709E-02,
     &   0.154278E-01,  0.112067E-01,  0.148690E-01,  0.154442E-01,
     &   0.123977E-01,  0.237539E-01,  0.162820E-01,  0.269484E-01,
     &   0.178081E-01,  0.143221E-01,  0.262468E-01,  0.217065E-01,
     &   0.107083E-01,  0.281220E-01,  0.115565E-01,  0.231244E-01,
     &   0.225197E-01,  0.178624E-01,  0.327708E-01,  0.116657E-01,
     &   0.277452E-01,  0.301647E-01,  0.349782E-01/
      DATA BP1     /
     &   0.717848E-05,  0.169280E-04,  0.126710E-04,  0.758397E-05,
     &  -0.533900E-05,  0.143490E-04, -0.595854E-05,  0.296465E-05,
     &   0.323446E-05,  0.115359E-04, -0.692861E-05,  0.131477E-04,
     &  -0.624945E-05, -0.756955E-06,  0.107458E-05, -0.159796E-05,
     &  -0.290529E-04, -0.170918E-05, -0.193934E-04, -0.707209E-04,
     &  -0.148154E-04, -0.383162E-04, -0.186050E-04, -0.951796E-04,
     &  -0.210944E-04, -0.330590E-04, -0.373087E-04, -0.408972E-04,
     &  -0.396759E-04, -0.827756E-04, -0.573773E-04, -0.325384E-04,
     &  -0.449411E-04, -0.271450E-04, -0.752791E-04, -0.549699E-04,
     &  -0.225655E-04, -0.102034E-03, -0.740322E-05, -0.668846E-04,
     &  -0.106063E-03, -0.304840E-04, -0.796023E-04,  0.504880E-04,
     &   0.486384E-04, -0.531946E-04, -0.147771E-03, -0.406785E-04,
     &   0.615750E-05, -0.486264E-04, -0.419335E-04, -0.819467E-04,
     &  -0.709498E-04,  0.326984E-05, -0.369743E-04, -0.526848E-04,
     &  -0.550050E-04, -0.684057E-04, -0.447093E-04, -0.778390E-04,
     &  -0.982953E-04, -0.772497E-04, -0.119430E-05, -0.655187E-04/
      DATA BP2     /
     &  -0.339078E-04,  0.716657E-04, -0.335893E-04,  0.220239E-04,
     &  -0.491012E-04, -0.393325E-04, -0.626461E-04, -0.795479E-04,
     &  -0.599181E-04, -0.578153E-04, -0.597559E-05, -0.866750E-04,
     &  -0.486783E-04, -0.580912E-04, -0.647368E-04, -0.350643E-04,
     &  -0.566635E-04, -0.385738E-04, -0.463782E-04, -0.321485E-04,
     &  -0.177300E-04, -0.250201E-04, -0.365492E-04, -0.165218E-04,
     &  -0.649177E-05, -0.218458E-04, -0.984604E-05, -0.120034E-04,
     &  -0.110119E-06, -0.164405E-04, -0.141396E-04,  0.315347E-05,
     &  -0.141544E-05, -0.297320E-05, -0.216248E-05,  0.839264E-05,
     &  -0.178197E-04, -0.106225E-04, -0.468195E-05,  0.997043E-05,
     &   0.679709E-05,  0.324610E-05, -0.367325E-05,  0.671058E-05,
     &   0.509293E-05, -0.437392E-05, -0.787922E-06, -0.271503E-06,
     &  -0.437940E-05, -0.128205E-04, -0.417830E-04, -0.561134E-05,
     &  -0.209940E-04, -0.414366E-04, -0.289765E-04,  0.680406E-06,
     &  -0.558644E-05, -0.530395E-05, -0.622242E-04, -0.159979E-05,
     &  -0.140286E-04, -0.128463E-04, -0.929499E-05, -0.327886E-04/
      DATA BP3     /
     &  -0.189353E-04, -0.737589E-04, -0.323471E-04, -0.272502E-04,
     &  -0.321731E-04, -0.326958E-04, -0.509157E-04, -0.681890E-04,
     &  -0.362182E-04, -0.354405E-04, -0.578392E-04,  0.238627E-05,
     &  -0.709028E-04, -0.518717E-04, -0.491859E-04, -0.718017E-04,
     &  -0.418978E-05, -0.940819E-04, -0.630375E-04, -0.478469E-04,
     &  -0.751896E-04, -0.267113E-04, -0.109019E-03, -0.890983E-04,
     &  -0.177301E-04, -0.120216E-03,  0.220464E-04, -0.734277E-04,
     &  -0.868068E-04, -0.652319E-04, -0.136982E-03, -0.279933E-06,
     &  -0.791824E-04, -0.111781E-03, -0.748263E-04/
      DATA ATP1    /
     &  -0.722782E-02, -0.901531E-02, -0.821263E-02, -0.808024E-02,
     &  -0.320169E-02, -0.661305E-02, -0.287272E-02, -0.486143E-02,
     &  -0.242857E-02, -0.530288E-02, -0.146813E-02, -0.566474E-03,
     &  -0.102192E-02,  0.300643E-03, -0.331655E-02,  0.648220E-03,
     &   0.552446E-02, -0.933046E-03,  0.205703E-02,  0.130638E-01,
     &  -0.229828E-02,  0.715648E-02,  0.444446E-03,  0.193500E-01,
     &   0.364119E-02,  0.252713E-02,  0.102420E-01,  0.494224E-02,
     &   0.584934E-02,  0.146255E-01,  0.921986E-02,  0.768012E-02,
     &   0.916105E-02,  0.276223E-02,  0.125245E-01,  0.131146E-01,
     &   0.793016E-02,  0.201536E-01,  0.658631E-02,  0.171711E-01,
     &   0.228470E-01,  0.131306E-01,  0.226658E-01,  0.176086E-01,
     &   0.149987E-01,  0.143060E-01,  0.313189E-01,  0.117070E-01,
     &   0.133522E-01,  0.244259E-01,  0.148393E-01,  0.223982E-01,
     &   0.151792E-01,  0.180474E-01,  0.106299E-01,  0.191016E-01,
     &   0.171776E-01,  0.229724E-01,  0.275530E-01,  0.302731E-01,
     &   0.281662E-01,  0.199525E-01,  0.192588E-01,  0.173220E-01/
      DATA ATP2    /
     &   0.195220E-01,  0.169371E-01,  0.193212E-01,  0.145558E-01,
     &   0.189654E-01,  0.122030E-01,  0.186206E-01,  0.228842E-01,
     &   0.139343E-01,  0.164006E-01,  0.137276E-01,  0.154005E-01,
     &   0.114575E-01,  0.129956E-01,  0.115305E-01,  0.929260E-02,
     &   0.106359E-01,  0.771623E-02,  0.106075E-01,  0.597630E-02,
     &   0.493960E-02,  0.532554E-02,  0.646175E-02,  0.302693E-02,
     &   0.150899E-02,  0.310333E-02,  0.533734E-02,  0.239094E-03,
     &   0.356782E-02,  0.707574E-02,  0.215758E-02, -0.527589E-03,
     &   0.643893E-03, -0.101916E-02, -0.383336E-02, -0.445966E-02,
     &   0.880190E-02,  0.245662E-02, -0.560923E-03, -0.582201E-02,
     &  -0.323233E-02, -0.454197E-02, -0.240905E-02, -0.343160E-02,
     &  -0.335156E-02, -0.623846E-03,  0.393633E-03, -0.271593E-02,
     &  -0.675874E-03,  0.920642E-03,  0.102168E-01, -0.250663E-03,
     &   0.437126E-02,  0.767434E-02,  0.569931E-02, -0.929326E-03,
     &   0.659414E-02,  0.280687E-02,  0.127614E-01,  0.780789E-02,
     &   0.374807E-02,  0.274288E-02,  0.534940E-02,  0.104349E-01/
       DATA ATP3   /
     &   0.294379E-02,  0.177846E-01,  0.523249E-02,  0.125339E-01,
     &   0.548538E-02,  0.577403E-02,  0.101532E-01,  0.170375E-01,
     &   0.758396E-02,  0.113402E-01,  0.106960E-01,  0.107782E-01,
     &   0.136148E-01,  0.992064E-02,  0.167276E-01,  0.149603E-01,
     &   0.136259E-01,  0.234521E-01,  0.166806E-01,  0.298505E-01,
     &   0.167592E-01,  0.186679E-01,  0.233062E-01,  0.228467E-01,
     &   0.128947E-01,  0.293979E-01,  0.219815E-01,  0.220663E-01,
     &   0.272710E-01,  0.237139E-01,  0.331743E-01,  0.208799E-01,
     &   0.281472E-01,  0.318440E-01,  0.370962E-01/
      DATA BTP1    /
     &   0.149748E-04,  0.188007E-04,  0.196530E-04,  0.124747E-04,
     &  -0.215751E-07,  0.128357E-04, -0.265798E-05,  0.606262E-05,
     &   0.287668E-05,  0.974612E-05, -0.833451E-05,  0.584410E-05,
     &  -0.452879E-05, -0.782537E-05,  0.786165E-05, -0.768351E-05,
     &  -0.196168E-04,  0.177297E-06, -0.129258E-04, -0.642798E-04,
     &  -0.986297E-05, -0.257145E-04, -0.141996E-04, -0.865089E-04,
     &  -0.141691E-04, -0.272578E-04, -0.295198E-04, -0.308878E-04,
     &  -0.313193E-04, -0.669272E-04, -0.475777E-04, -0.221332E-04,
     &  -0.419930E-04, -0.102519E-04, -0.590184E-04, -0.574771E-04,
     &  -0.240809E-04, -0.913994E-04, -0.908886E-05, -0.721074E-04,
     &  -0.902837E-04, -0.447582E-04, -0.664544E-04, -0.143150E-04,
     &  -0.511866E-05, -0.559352E-04, -0.104734E-03, -0.305206E-04,
     &   0.103303E-04, -0.613019E-04, -0.320040E-04, -0.738909E-04,
     &  -0.388263E-04,  0.306515E-04, -0.352214E-04, -0.253940E-04,
     &  -0.521369E-04, -0.746260E-04, -0.744124E-04, -0.881905E-04,
     &  -0.933645E-04, -0.664045E-04, -0.570712E-05, -0.566312E-04/
      DATA BTP2    /
     &  -0.364967E-04,  0.393501E-06, -0.234050E-04, -0.141317E-04,
     &  -0.525480E-04, -0.172241E-04, -0.410843E-04, -0.358348E-04,
     &  -0.256168E-04, -0.509482E-04, -0.180570E-04, -0.555356E-04,
     &  -0.271464E-04, -0.274040E-04, -0.480889E-04, -0.275751E-04,
     &  -0.415681E-04, -0.383770E-04, -0.280139E-04, -0.287919E-04,
     &  -0.125865E-04, -0.265467E-04, -0.172765E-04, -0.164611E-04,
     &   0.189183E-04, -0.171219E-04, -0.132766E-04, -0.344611E-05,
     &  -0.442832E-05, -0.185779E-04, -0.139755E-04,  0.168083E-05,
     &  -0.395287E-05, -0.297871E-05,  0.434383E-05,  0.131741E-04,
     &  -0.192637E-04, -0.549551E-05,  0.122553E-05,  0.204627E-04,
     &   0.154027E-04,  0.953462E-05,  0.131125E-05,  0.732839E-05,
     &   0.755405E-05, -0.305552E-05, -0.434858E-05,  0.308409E-05,
     &  -0.164787E-05, -0.818533E-05, -0.355041E-04, -0.504696E-05,
     &  -0.229022E-04, -0.356891E-04, -0.230346E-04,  0.518835E-05,
     &  -0.160187E-04, -0.104617E-04, -0.464754E-04, -0.115807E-04,
     &  -0.130230E-04, -0.603491E-05, -0.125324E-04, -0.165516E-04/
      DATA BTP3    /
     &  -0.991679E-05, -0.529432E-04, -0.200199E-04, -0.181977E-04,
     &  -0.220940E-04, -0.204483E-04, -0.432584E-04, -0.449109E-04,
     &  -0.247305E-04, -0.174253E-04, -0.484446E-04,  0.354150E-04,
     &  -0.425581E-04, -0.406562E-04, -0.505495E-04, -0.651856E-04,
     &  -0.153953E-04, -0.894294E-04, -0.616551E-04, -0.846504E-04,
     &  -0.699414E-04, -0.376203E-04, -0.940985E-04, -0.753050E-04,
     &  -0.183710E-04, -0.123907E-03, -0.279347E-04, -0.736381E-04,
     &  -0.103588E-03, -0.754117E-04, -0.140991E-03, -0.366687E-04,
     &  -0.927785E-04, -0.125321E-03, -0.115290E-03/
      DATA BETAD1  /
     &   0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &   0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &   0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &   0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &   0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &   0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &   0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &   0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &   0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &   0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &   0.234879E+03,  0.217419E+03,  0.201281E+03,  0.186364E+03,
     &   0.172576E+03,  0.159831E+03,  0.148051E+03,  0.137163E+03,
     &   0.127099E+03,  0.117796E+03,  0.109197E+03,  0.101249E+03,
     &   0.939031E+02,  0.871127E+02,  0.808363E+02,  0.750349E+02,
     &   0.497489E+02,  0.221212E+02,  0.113124E+02,  0.754174E+01,
     &   0.589554E+01,  0.495227E+01,  0.000000E+00,  0.000000E+00/
      DATA BETAD2  /
     &   0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &   0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &   0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &   0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &   0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &   0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &   0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &   0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &   0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &   0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &   0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &   0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &   0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &   0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &   0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &   0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00/
      DATA BETAD3  /
     &   0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &   0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &   0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &   0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &   0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &   0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &   0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &   0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,
     &   0.000000E+00,  0.000000E+00,  0.000000E+00/
      DATA BANDL1 /
     &   0.000000E+00,  0.100000E+02,  0.200000E+02,  0.300000E+02,
     &   0.400000E+02,  0.500000E+02,  0.600000E+02,  0.700000E+02,
     &   0.800000E+02,  0.900000E+02,  0.100000E+03,  0.110000E+03,
     &   0.120000E+03,  0.130000E+03,  0.140000E+03,  0.150000E+03,
     &   0.160000E+03,  0.170000E+03,  0.180000E+03,  0.190000E+03,
     &   0.200000E+03,  0.210000E+03,  0.220000E+03,  0.230000E+03,
     &   0.240000E+03,  0.250000E+03,  0.260000E+03,  0.270000E+03,
     &   0.280000E+03,  0.290000E+03,  0.300000E+03,  0.310000E+03,
     &   0.320000E+03,  0.330000E+03,  0.340000E+03,  0.350000E+03,
     &   0.360000E+03,  0.370000E+03,  0.380000E+03,  0.390000E+03,
     &   0.400000E+03,  0.410000E+03,  0.420000E+03,  0.430000E+03,
     &   0.440000E+03,  0.450000E+03,  0.460000E+03,  0.470000E+03,
     &   0.480000E+03,  0.490000E+03,  0.500000E+03,  0.510000E+03,
     &   0.520000E+03,  0.530000E+03,  0.540000E+03,  0.550000E+03,
     &   0.560000E+03,  0.670000E+03,  0.800000E+03,  0.900000E+03,
     &   0.990000E+03,  0.107000E+04,  0.120000E+04,  0.121000E+04/
      DATA BANDL2 /
     &   0.122000E+04,  0.123000E+04,  0.124000E+04,  0.125000E+04,
     &   0.126000E+04,  0.127000E+04,  0.128000E+04,  0.129000E+04,
     &   0.130000E+04,  0.131000E+04,  0.132000E+04,  0.133000E+04,
     &   0.134000E+04,  0.135000E+04,  0.136000E+04,  0.137000E+04,
     &   0.138000E+04,  0.139000E+04,  0.140000E+04,  0.141000E+04,
     &   0.142000E+04,  0.143000E+04,  0.144000E+04,  0.145000E+04,
     &   0.146000E+04,  0.147000E+04,  0.148000E+04,  0.149000E+04,
     &   0.150000E+04,  0.151000E+04,  0.152000E+04,  0.153000E+04,
     &   0.154000E+04,  0.155000E+04,  0.156000E+04,  0.157000E+04,
     &   0.158000E+04,  0.159000E+04,  0.160000E+04,  0.161000E+04,
     &   0.162000E+04,  0.163000E+04,  0.164000E+04,  0.165000E+04,
     &   0.166000E+04,  0.167000E+04,  0.168000E+04,  0.169000E+04,
     &   0.170000E+04,  0.171000E+04,  0.172000E+04,  0.173000E+04,
     &   0.174000E+04,  0.175000E+04,  0.176000E+04,  0.177000E+04,
     &   0.178000E+04,  0.179000E+04,  0.180000E+04,  0.181000E+04,
     &   0.182000E+04,  0.183000E+04,  0.184000E+04,  0.185000E+04/
      DATA BANDL3 /
     &   0.186000E+04,  0.187000E+04,  0.188000E+04,  0.189000E+04,
     &   0.190000E+04,  0.191000E+04,  0.192000E+04,  0.193000E+04,
     &   0.194000E+04,  0.195000E+04,  0.196000E+04,  0.197000E+04,
     &   0.198000E+04,  0.199000E+04,  0.200000E+04,  0.201000E+04,
     &   0.202000E+04,  0.203000E+04,  0.204000E+04,  0.205000E+04,
     &   0.206000E+04,  0.207000E+04,  0.208000E+04,  0.209000E+04,
     &   0.210000E+04,  0.211000E+04,  0.212000E+04,  0.213000E+04,
     &   0.214000E+04,  0.215000E+04,  0.216000E+04,  0.217000E+04,
     &   0.218000E+04,  0.219000E+04,  0.227000E+04/
      DATA BANDH1 /
     &   0.100000E+02,  0.200000E+02,  0.300000E+02,  0.400000E+02,
     &   0.500000E+02,  0.600000E+02,  0.700000E+02,  0.800000E+02,
     &   0.900000E+02,  0.100000E+03,  0.110000E+03,  0.120000E+03,
     &   0.130000E+03,  0.140000E+03,  0.150000E+03,  0.160000E+03,
     &   0.170000E+03,  0.180000E+03,  0.190000E+03,  0.200000E+03,
     &   0.210000E+03,  0.220000E+03,  0.230000E+03,  0.240000E+03,
     &   0.250000E+03,  0.260000E+03,  0.270000E+03,  0.280000E+03,
     &   0.290000E+03,  0.300000E+03,  0.310000E+03,  0.320000E+03,
     &   0.330000E+03,  0.340000E+03,  0.350000E+03,  0.360000E+03,
     &   0.370000E+03,  0.380000E+03,  0.390000E+03,  0.400000E+03,
     &   0.410000E+03,  0.420000E+03,  0.430000E+03,  0.440000E+03,
     &   0.450000E+03,  0.460000E+03,  0.470000E+03,  0.480000E+03,
     &   0.490000E+03,  0.500000E+03,  0.510000E+03,  0.520000E+03,
     &   0.530000E+03,  0.540000E+03,  0.550000E+03,  0.560000E+03,
     &   0.670000E+03,  0.800000E+03,  0.900000E+03,  0.990000E+03,
     &   0.107000E+04,  0.120000E+04,  0.121000E+04,  0.122000E+04/
      DATA BANDH2 /
     &   0.123000E+04,  0.124000E+04,  0.125000E+04,  0.126000E+04,
     &   0.127000E+04,  0.128000E+04,  0.129000E+04,  0.130000E+04,
     &   0.131000E+04,  0.132000E+04,  0.133000E+04,  0.134000E+04,
     &   0.135000E+04,  0.136000E+04,  0.137000E+04,  0.138000E+04,
     &   0.139000E+04,  0.140000E+04,  0.141000E+04,  0.142000E+04,
     &   0.143000E+04,  0.144000E+04,  0.145000E+04,  0.146000E+04,
     &   0.147000E+04,  0.148000E+04,  0.149000E+04,  0.150000E+04,
     &   0.151000E+04,  0.152000E+04,  0.153000E+04,  0.154000E+04,
     &   0.155000E+04,  0.156000E+04,  0.157000E+04,  0.158000E+04,
     &   0.159000E+04,  0.160000E+04,  0.161000E+04,  0.162000E+04,
     &   0.163000E+04,  0.164000E+04,  0.165000E+04,  0.166000E+04,
     &   0.167000E+04,  0.168000E+04,  0.169000E+04,  0.170000E+04,
     &   0.171000E+04,  0.172000E+04,  0.173000E+04,  0.174000E+04,
     &   0.175000E+04,  0.176000E+04,  0.177000E+04,  0.178000E+04,
     &   0.179000E+04,  0.180000E+04,  0.181000E+04,  0.182000E+04,
     &   0.183000E+04,  0.184000E+04,  0.185000E+04,  0.186000E+04/
      DATA BANDH3 /
     &   0.187000E+04,  0.188000E+04,  0.189000E+04,  0.190000E+04,
     &   0.191000E+04,  0.192000E+04,  0.193000E+04,  0.194000E+04,
     &   0.195000E+04,  0.196000E+04,  0.197000E+04,  0.198000E+04,
     &   0.199000E+04,  0.200000E+04,  0.201000E+04,  0.202000E+04,
     &   0.203000E+04,  0.204000E+04,  0.205000E+04,  0.206000E+04,
     &   0.207000E+04,  0.208000E+04,  0.209000E+04,  0.210000E+04,
     &   0.211000E+04,  0.212000E+04,  0.213000E+04,  0.214000E+04,
     &   0.215000E+04,  0.216000E+04,  0.217000E+04,  0.218000E+04,
     &   0.219000E+04,  0.220000E+04,  0.238000E+04/
      real (kind=kind_rad) ALB1(21,7),ALB2(21,7),ALB3(21,6)
!
      DATA ALB1/ .061,.062,.072,.087,.115,.163,.235,.318,.395,.472,.542,
     & .604,.655,.693,.719,.732,.730,.681,.581,.453,.425,.061,.062,.070,
     & .083,.108,.145,.198,.263,.336,.415,.487,.547,.595,.631,.656,.670,
     & .652,.602,.494,.398,.370,.061,.061,.068,.079,.098,.130,.174,.228,
     & .290,.357,.424,.498,.556,.588,.603,.592,.556,.488,.393,.342,.325,
     & .061,.061,.065,.073,.086,.110,.150,.192,.248,.306,.360,.407,.444,
     & .469,.480,.474,.444,.386,.333,.301,.290,.061,.061,.065,.070,.082,
     & .101,.131,.168,.208,.252,.295,.331,.358,.375,.385,.377,.356,.320,
     & .288,.266,.255,.061,.061,.063,.068,.077,.092,.114,.143,.176,.210,
     & .242,.272,.288,.296,.300,.291,.273,.252,.237,.266,.220,.061,.061,
     & .062,.066,.072,.084,.103,.127,.151,.176,.198,.219,.236,.245,.250,
     & .246,.235,.222,.211,.205,.200/
      DATA ALB2/ .061,.061,.061,.065,.071,.079,.094,.113,.134,.154,.173,
     & .185,.190,.193,.193,.190,.188,.185,.182,.180,.178,.061,.061,.061,
     & .064,.067,.072,.083,.099,.117,.135,.150,.160,.164,.165,.164,.162,
     & .160,.159,.158,.157,.157,.061,.061,.061,.062,.065,.068,.074,.084,
     & .097,.111,.121,.127,.130,.131,.131,.130,.129,.127,.126,.125,.122,
     & .061,.061,.061,.061,.062,.064,.070,.076,.085,.094,.101,.105,.107,
     & .106,.103,.100,.097,.096,.095,.095,.095,.061,.061,.061,.060,.061,
     & .062,.065,.070,.075,.081,.086,.089,.090,.088,.084,.080,.077,.075,
     & .074,.074,.074,.061,.061,.060,.060,.060,.061,.063,.065,.068,.072,
     & .076,.077,.076,.074,.071,.067,.064,.062,.061,.061,.061,.061,.061,
     & .060,.060,.060,.060,.061,.062,.065,.068,.069,.069,.068,.065,.061,
     & .058,.055,.054,.053,.052,.052/
      DATA ALB3/ .061,.061,.060,.060,.060,.060,.060,.060,.062,.065,.065,
     & .063,.060,.057,.054,.050,.047,.046,.045,.044,.044,.061,.061,.060,
     & .060,.060,.059,.059,.059,.059,.059,.058,.055,.051,.047,.043,.039,
     & .035,.033,.032,.031,.031,.061,.061,.060,.060,.060,.059,.059,.058,
     & .057,.056,.054,.051,.047,.043,.039,.036,.033,.030,.028,.027,.026,
     & .061,.061,.060,.060,.060,.059,.059,.058,.057,.055,.052,.049,.045,
     & .040,.036,.032,.029,.027,.026,.025,.025,.061,.061,.060,.060,.060,
     & .059,.059,.058,.056,.053,.050,.046,.042,.038,.034,.031,.028,.026,
     & .025,.025,.025,.061,.061,.060,.060,.059,.058,.058,.057,.055,.053,
     & .050,.046,.042,.038,.034,.030,.028,.029,.025,.025,.025/
!
!   INITIALISATION
      ARNDM(1:64)=ARNDM1(1:64)
      ARNDM(65:128)=ARNDM2(1:64)
      ARNDM(129:163)=ARNDM3(1:35)
      BRNDM(1:64)=BRNDM1(1:64)
      BRNDM(65:128)=BRNDM2(1:64)
      BRNDM(129:163)=BRNDM3(1:35)
      AP(1:64)=AP1(1:64)
      AP(65:128)=AP2(1:64)
      AP(129:163)=AP3(1:35)
      BP(1:64)=BP1(1:64)
      BP(65:128)=BP2(1:64)
      BP(129:163)=BP3(1:35)
      ATP(1:64)=ATP1(1:64)
      ATP(65:128)=ATP2(1:64)
      ATP(129:163)=ATP3(1:35)
      BTP(1:64)=BTP1(1:64)
      BTP(65:128)=BTP2(1:64)
      BTP(129:163)=BTP3(1:35)
      BETAD(1:64)=BETAD1(1:64)
      BETAD(65:128)=BETAD2(1:64)
      BETAD(129:163)=BETAD3(1:35)
      BANDLO(1:64)=BANDL1(1:64)
      BANDLO(65:128)=BANDL2(1:64)
      BANDLO(129:163)=BANDL3(1:35)
      BANDHI(1:64)=BANDH1(1:64)
      BANDHI(65:128)=BANDH2(1:64)
      BANDHI(129:163)=BANDH3(1:35)
      ALBD(1:21,1:7)=ALB1(1:21,1:7)
      ALBD(1:21,8:14)=ALB2(1:21,1:7)
      ALBD(1:21,15:20)=ALB3(1:21,1:6)
!****************************************
!***COMPUTE LOCAL QUANTITIES AND AO3,BO3,AB15
!....FOR NARROW-BANDS...
!
      L   = LP1 - 1
      LP2 = LP1 + 1
!
      DO N=1,NBLW
        ANB(N)    = ARNDM(N)
        BNB(N)    = BRNDM(N)
        CENTNB(N) = HAF * (BANDLO(N) + BANDHI(N))
        DELNB(N)  = BANDHI(N) - BANDLO(N)
        BETANB(N) = BETAD(N)
      ENDDO
      AB15(1) = ANB(57)*BNB(57)
      AB15(2) = ANB(58)*BNB(58)
!....FOR WIDE BANDS...
      AB15WD = AWIDE*BWIDE
!
!***COMPUTE RATIOS OF CONT. COEFFS
      SKC1R = BETAWD    / BETINW
      SKO3R = BETAD(61) / BETINW
      SKO2D = ONE       / BETINW
!
!****BEGIN TABLE COMPUTATIONS HERE***
!***COMPUTE TEMPS, MASSES FOR TABLE ENTRIES
!---NOTE: THE DIMENSIONING AND INITIALIZATION OF XTEMV AND OTHER ARRAYS
!   WITH DIMENSION OF 28 IMPLY A RESTRICTION OF MODEL TEMPERATURES FROM
!   100K TO 370K.
!---THE DIMENSIONING OF ZMASS,ZROOT AND OTHER ARRAYS WITH DIMENSION OF
!   180 IMPLY A RESTRICTION OF MODEL H2O AMOUNTS SUCH THAT OPTICAL PATHS
!   ARE BETWEEN 10**-16 AND 10**2, IN CGS UNITS.
      ZMASS(1) = H1M16
      DO J=1,180
        JP        = J + 1
        ZROOT(J)  = SQRT(ZMASS(J))
        ZMASS(JP) = ZMASS(J)*H1P25892
      ENDDO
      DO I=1,28
        XTEMV(I)  = HNINETY + TEN*I
        TFOUR(I)  = 1.0 / (XTEMV(I)*XTEMV(I)*XTEMV(I)*XTEMV(I))
        FORTCU(I) = 1.0 / (FOUR*XTEMV(I)*XTEMV(I)*XTEMV(I))
      ENDDO
!******THE COMPUTATION OF SOURCE,DSRCE IS  NEEDED ONLY
!   FOR THE COMBINED WIDE-BAND CASE.TO OBTAIN THEM,THE SOURCE
!   MUST BE COMPUTED FOR EACH OF THE (NBLX) WIDE BANDS(=SRCWD)
!   THEN COMBINED (USING IBAND) INTO SOURCE.
      DO N=1,NBLY
        DO I=1,28
          SOURCE(I,N) = ZERO
        ENDDO
      ENDDO
      DO N=1,NBLX
        DO I=1,28
          SRCWD(I,N)=ZERO
        ENDDO
      ENDDO
!---BEGIN FREQ. LOOP (ON N)
      DO N=1,NBLX
        IF (N.LE.46) THEN
!***THE 160-1200 BAND CASES
          CENT = CENTNB(N+16)
          DEL  = DELNB(N+16)
          BDLO = BANDLO(N+16)
          BDHI = BANDHI(N+16)
        ENDIF
        IF (N.EQ.NBLX) THEN
!***THE 2270-2380 BAND CASE
          CENT = CENTNB(NBLW)
          DEL  = DELNB(NBLW)
          BDLO = BANDLO(NBLW)
          BDHI = BANDHI(NBLW)
        ENDIF
!***FOR PURPOSES OF ACCURACY, ALL EVALUATIONS OF PLANCK FCTNS ARE MADE
!  ON 10 CM-1 INTERVALS, THEN SUMMED INTO THE (NBLX) WIDE BANDS.
        NSUBDS = (DEL-H1M3)/10+1
        DO NSB=1,NSUBDS
          IF (NSB.NE.NSUBDS) THEN
            CNUSB(NSB) = TEN*(NSB-1) + BDLO + FIVE
            DNUSB(NSB) = TEN
          ELSE
            CNUSB(NSB) = HAF * (TEN*(NSB-1)+BDLO+BDHI)
            DNUSB(NSB) = BDHI - (TEN*(NSB-1)+BDLO)
          ENDIF
          C1 = (H37412M5)*CNUSB(NSB)**3
!---BEGIN TEMP. LOOP (ON I)
          DO I=1,28
            X(I)       = H1P4387*CNUSB(NSB) / XTEMV(I)
            X1(I)      = EXP(X(I))
            SRCS(I)    = C1 / (X1(I)-ONE)
            SRCWD(I,N) = SRCWD(I,N)+SRCS(I)*DNUSB(NSB)
          ENDDO
        ENDDO
      ENDDO                          ! End of N Loop!
!***THE FOLLOWING LOOPS CREATE THE COMBINED WIDE BAND QUANTITIES SOURCE
!   AND DSRCE
      DO N=1,40
        DO I=1,28
          SOURCE(I,IBAND(N)) = SOURCE(I,IBAND(N)) + SRCWD(I,N)
        ENDDO
      ENDDO
      DO N=9,NBLY
        DO I=1,28
          SOURCE(I,N) = SRCWD(I,N+32)
        ENDDO
      ENDDO
      DO N=1,NBLY
        DO I=1,27
          DSRCE(I,N) = (SOURCE(I+1,N)-SOURCE(I,N))*HP1
        ENDDO
      ENDDO
      DO N=1,NBLW
        ALFANB(N) = BNB(N)*ANB(N)
        AROTNB(N) = SQRT(ALFANB(N))
      ENDDO
!***FIRST COMPUTE PLANCK FCTNS (SRC1NB) AND DERIVATIVES (DBDTNB) FOR
!   USE IN TABLE EVALUATIONS. THESE ARE DIFFERENT FROM SOURCE,DSRCE
!   BECAUSE DIFFERENT FREQUENCY PTS ARE USED IN EVALUATION, THE FREQ.
!   RANGES ARE DIFFERENT, AND THE DERIVATIVE ALGORITHM IS DIFFERENT.
!
      DO N=1,NBLW
        CENT = CENTNB(N)
        DEL  = DELNB(N)
!---NOTE: AT PRESENT, THE IA LOOP IS ONLY USED FOR IA=2. THE LOOP STRUCT
!   IS KEPT SO THAT IN THE FUTURE, WE MAY USE A QUADRATURE SCHEME FOR
!   THE PLANCK FCTN EVALUATION, RATHER THAN USE THE MID-BAND FREQUENCY.
        DO IA=1,3
          ANU = CENT + HAF*(IA-2)*DEL
          C1  = (H37412M5)*ANU*ANU*ANU + H1M20
!---TEMPERATURE LOOP---
          DO I=1,28
            X(I)   = H1P4387 * ANU / XTEMV(I)
            X1(I)  = EXP(X(I))
            SC(I)  = C1 / ((X1(I)-ONE)+H1M20)
            DSC(I) = SC(I)*SC(I)*X(I)*X1(I) / (XTEMV(I)*C1)
          ENDDO
          IF (IA.EQ.2) THEN
            DO I=1,28
              SRC1NB(I,N ) =DEL * SC(I)
              DBDTNB(I,N) = DEL * DSC(I)
            ENDDO
          ENDIF
        ENDDO
      ENDDO
!***NEXT COMPUTE R1,R2,S2,AND T3- COEFFICIENTS USED FOR E3 FUNCTION
!   WHEN THE OPTICAL PATH IS LESS THAN 10-4. IN THIS CASE, WE ASSUME A
!   DIFFERENT DEPENDENCE ON (ZMASS).
!---ALSO OBTAIN R1WD, WHICH IS R1 SUMMED OVER THE 160-560 CM-1 RANGE
      DO I=1,28
        SUM4(I)   = ZERO
        SUM6(I)   = ZERO
        SUM7(I)   = ZERO
        SUM8(I)   = ZERO
        SUM4WD(I) = ZERO
      ENDDO
      DO N=1,NBLW
        CENT=CENTNB(N)
!***PERFORM SUMMATIONS FOR FREQ. RANGES OF 0-560,1200-2200 CM-1 FOR SUM4
!   SUM6,SUM7,SUM8
        IF (CENT.LT.560. .OR. CENT.GT.1200..AND.CENT.LE.2200.) THEN
          DO I=1,28
            SUM4(I) = SUM4(I) + SRC1NB(I,N)
            SUM6(I) = SUM6(I) + DBDTNB(I,N)
            SUM7(I) = SUM7(I) + DBDTNB(I,N) * AROTNB(N)
            SUM8(I) = SUM8(I) + DBDTNB(I,N) * ALFANB(N)
          ENDDO
        ENDIF
!***PERFORM SUMMATIONS OVER 160-560 CM-1 FREQ RANGE FOR E1 CALCS (SUM4WD
        IF (CENT.GT.160. .AND. CENT.LT.560.) THEN
          DO I=1,28
            SUM4WD(I) = SUM4WD(I) + SRC1NB(I,N)
          ENDDO
        ENDIF
      ENDDO
      DO I=1,28
        R1(I)   = SUM4(I)   * TFOUR(I)
        R2(I)   = SUM6(I)   * FORTCU(I)
        S2(I)   = SUM7(I)   * FORTCU(I)
        T3(I)   = SUM8(I)   * FORTCU(I)
        R1WD(I) = SUM4WD(I) * TFOUR(I)
      ENDDO
      DO J=1,180
        DO I=1,28
          SUM(I,J)    = ZERO
          PERTSM(I,J) = ZERO
          SUM3(I,J)   = ZERO
          SUMWDE(I,J) = ZERO
        ENDDO
      ENDDO
!---FREQUENCY LOOP BEGINS---
      DO N=1,NBLW
        CENT = CENTNB(N)
!***PERFORM CALCULATIONS FOR FREQ. RANGES OF 0-560,1200-2200 CM-1
        IF (CENT.LT.560. .OR. CENT.GT.1200..AND.CENT.LE.2200.) THEN
          DO J=1,180
            X2(J)   = AROTNB(N) * ZROOT(J)
            EXPO(J) = EXP(-X2(J))
          ENDDO
          DO J=1,180
            IF (X2(J).GE.HUNDRED) THEN
              EXPO(J) = ZERO
            ENDIF
          ENDDO
          DO J=121,180
            FAC(J) = ZMASS(J)*(ONE-(ONE+X2(J))*EXPO(J))/(X2(J)*X2(J))
          ENDDO
          DO J=1,180
            DO I=1,28
              SUM(I,J)    = SUM(I,J)    + SRC1NB(I,N)*EXPO(J)
              PERTSM(I,J) = PERTSM(I,J) + DBDTNB(I,N)*EXPO(J)
            ENDDO
          ENDDO
          DO J=121,180
            DO I=1,28
              SUM3(I,J) = SUM3(I,J) + DBDTNB(I,N)*FAC(J)
            ENDDO
          ENDDO
        ENDIF
!---COMPUTE SUM OVER 160-560 CM-1 RANGE FOR USE IN E1 CALCS (SUMWDE)
        IF (CENT.GT.160. .AND. CENT.LT.560.) THEN
          DO J=1,180
            DO I=1,28
              SUMWDE(I,J) = SUMWDE(I,J) + SRC1NB(I,N)*EXPO(J)
            ENDDO
          ENDDO
        ENDIF
      ENDDO
      I1=0
      DO J=1,180
        DO I=1,28
          I1=I1+1
          EM1V(I1)    = SUM(I,J)    * TFOUR(I)
          T1(I1) = PERTSM(I,J) * FORTCU(I)
        ENDDO
      ENDDO
      I1=120*28
      DO J=121,180
        DO I=1,28
          I1=I1+1
          EM3V(I1) = SUM3(I,J) * FORTCU(I)
        ENDDO
      ENDDO
      I1=0
      DO J=1,179
        DO I=1,28
          I1=I1+1
          T2(I1)=(T1(I1+28)-T1(I1))*TEN
        ENDDO
      ENDDO
      I1=0
      DO J=1,180
        DO I=1,27
          T4(I1+I)=(T1(I1+I+1)-T1(I1+I))*HP1
        ENDDO
        I1 = I1 + 28
        T4(I1)=ZERO
      ENDDO
      I1=179*28
      DO I=1,28
        I1=I1+1
        T2(I1)=ZERO
      ENDDO
      I1=0
      DO J=1,2
        DO I=1,28
          I1=I1+1
          EM1V(I1)=R1(I)
        ENDDO
      ENDDO
      I1=0
      DO J=1,120
        DO I=1,28
          I1=I1+1
          EM3V(I1)=R2(I)/TWO-S2(I)*SQRT(ZMASS(J))/THREE
     &             + T3(I)*ZMASS(J)/EIGHT
        ENDDO
      ENDDO
      I1=120*28
      DO J=121,180
        DO I=1,28
          I1=I1+1
          EM3V(I1)=EM3V(I1)/ZMASS(J)
        ENDDO
      ENDDO
!***NOW COMPUTE E1 TABLES FOR 160-560 CM-1 BANDS ONLY.
!   WE USE R1WD AND SUMWDE OBTAINED ABOVE.
      I1=0
      DO J=1,180
        DO I=1,28
          I1=I1+1
          EM1VW(I1)=SUMWDE(I,J)*TFOUR(I)
        ENDDO
      ENDDO
      I1=0
      DO J=1,2
        DO I=1,28
          I1=I1+1
          EM1VW(I1)=R1WD(I)
        ENDDO
      ENDDO
!
      RETURN
      END
