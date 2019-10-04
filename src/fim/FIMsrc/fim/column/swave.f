      SUBROUTINE SWR95(S0,ISRC,PL,TA,WA,OA,CO2,COSZ,TAUCL,
     &            CCLY,CFAC,ICFC,ICWP,CWP,CIP,REW,REI,FICE,
     &            ALBUVB,ALBUVD,ALBIRB,ALBIRD,KPRF,IDXC,CMIX,DENN,RH,
     &            HTRC,TUPFXC,TDNFLX,SUPFXC,SDNFXC,
     &            TUPFX0,SUPFX0,SDNFX0,
     &            SDNFVB,SDNFVD,SDNFNB,SDNFND, NDAY, IR,
! --- FOR UV-B BAND FLUXES
     &            SUVBFC,SUVBF0,
! --- END UV-B
     &            L, LP1, IMAX, NSRC, NBD, NVB, NAE, NDM, NXC, NDN,
     &            HAER, IDM, DZ, HZ, TAUR, me, ix2)
!    &,           lprnt)
CFPP$ NOCONCUR R
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:    SWR95      COMPUTES SHORT-WAVE RADIATIVE HEATING
!   PROGRAMMER: YU-TAI HOU  ORG: W/NMC20      DATE: 95-02-09
!
! ABSTRACT: THIS CODE IS A MODIFIED VERSION OF M.D. CHOU'S SW
!   RADIATION CODE TO FIT NMC MRF AND CLIMATE MODELS.  IT COMPUTES
!   SW ATMOSPHERIC ABSORPTION AND SCATTERING EFFECTS DUE TO O3,
!   H2O,CO2,O2,CLOUDS, AND AEROSOLS, ETC.
!   IT HAS 8 UV+VIS BANDS AND 3 NIR BANDS (10 K-VALUES EACH).
!
! REFERENCES: CHOU (1986, J. CLIM. APPL.METEOR.)
!   CHOU (1990, J. CLIM.), AND CHOU (1992, J. ATMS. SCI.)
!   CHOU AND SUAREZ (1999, NASA/TM-1999-104606,VOL.15)
!
! PROGRAM HISTORY LOG:
!   94-06-12   M.D. CHOU, GLA.
!   95-02-09   YU-TAI HOU      - RECODE FOR NMC MODELS
!   98-08-03   YU-TAI HOU      - UPDATED CLOUD RADIATIVE PROPERTIES
!            CALCULATION. USE SLINGO'S METHOD (JAS 1989) ON WATER
!            CLOUD, EBERT AND CURRY'S METHOD (JGR 1992) ON ICE CLOUD.
!   99-03-25   YU-TAI HOU      - UPDATED CLOUD RADIATIVE PROPERTY
!            CALCULATIONS USE CHOU ET AL. NEW METHOD (J. CLIM 1998)
!   99-04-27   YU-TAI HOU      - UPDATED CLOUD RADIATIVE PROPERTY
!            CALCULATIONS USE LINEAR T-ADJUSTED METHOD.
!   99-09-13   YU-TAI HOU      - UPDATED TO CHOU'S JUNE,99 VERSION
!
! USAGE:       CALL SWR95
!
! ATTRIBUTES:
!   LANGUAGE:  FORTRAN 77 & Fortran 90
!   MACHINE:   CRAY C-90, IBM SP, SGI
!
! INPUT PARAMETERS:
!   S0     : SOLAR CONSTANT
!   ISRC   : FLAGS FOR SELECTING ABSORBERS
!            1:AEROSOLS, 2:O2, 3:CO2, 4:H2O, 5:O3
!            =0:WITHOUT IT,  =1: WITH IT.
!   PL     : MODEL LEVEL PRESSURE IN MB
!   TA     : MODEL LAYER TEMPERATURE IN K
!   WA     : LAYER SPECIFIC HUMIDITY IN GM/GM
!   OA     : LAYER OZONE CONCENTRATION IN GM/GM
!   CO2    : CO2 MIXING RATION BY VOLUME
!   COSZ   : COSINE OF SOLAR ZENITH ANGLE
!   TAUCL  : OPTICAL DEPTH OF CLOUD LAYERS
!   CCLY   : LAYER CLOUD FRACTION
!   CFAC   : FRACTION OF CLEAR SKY VIEW AT THE LAYER INTERFACE
!   ICFC   : =0 NO CLOUD FACTOR TO WEIGH CLEAR AND CLOUDY FLUXES
!            =1 USE CLOUD FACTOR TO WEIGH CLEAR AND CLOUDY FLUXES
!   ICWP   : FLAG INDICATES THE METHOD USED FOR CLOUD PROPERTIES
!            CALCULATIONS, =0 USE T-P; =1 USE CWC/CIC.
!   CWP    : LAYER CLOUD WATER PATH (G/M**2)
!   CIP    : LAYER CLOUD ICE PATH (G/M**2)
!   REW    : LAYER WATER CLOUD DROP EFFECTIVE RADIUS (MICRON)
!   REI    : LAYER ICE CLOUD DROP EFFECTIVE RADIUS
!   FICE   : FRACTION OF CLOUD ICE CONTENT
!   ALBUVB : UV+VIS SURF DIRECT ALBEDO
!   ALBUVD : UV+VIS SURF DIFFUSED ALBEDO
!   ALBIRB : NIR SURF DIRECT ALBEDO
!   ALBIRD : NIR SURF DIFFUSED ALBEDO
!   PAER   : AEROSOL PROFILES (FRACTION)
!
! OUTPUT PARAMETER:
!   HTRC   : HEATING RATES FOR CLOUDY SKY IN  K/DAY
!   TUPFXC : UPWARD FLUX AT TOA FOR CLOUDY SKY  W/M**2
!   TDNFLX : DNWARD FLUX AT TOA FOR ALL SKY  W/M**2
!   SUPFXC : UPWARD FLUX AT SFC FOR CLOUDY SKY  W/M**2
!   SDNFXC : DNWARD FLUX AT SFC FOR CLOUDY SKY  W/M**2
!   TUPFX0 : UPWARD FLUX AT TOA FOR CLEAR SKY   W/M**2
!   SUPFX0 : UPWARD FLUX AT SFC FOR CLEAR SKY   W/M**2
!   SDNFX0 : DNWARD FLUX AT SFC FOR CLEAR SKY   W/M**2
!   SDNFVB : DOWNWARD SURFACE VIS BEAM FLUX     W/M**2
!   SDNFNB : DOWNWARD SURFACE NIR BEAM FLUX     W/M**2
!   SDNFVD : DOWNWARD SURFACE VIS DIFF FLUX     W/M**2
!   SDNFND : DOWNWARD SURFACE NIR DIFF FLUX     W/M**2
!
! NOTE:
!   FOR ALL QUANTITIES, K=1 IS THE TOP LEVEL/LAYER, EXCEPT
!   SI AND SL, FOR WHICH K=1 IS THE SURFACE LEVEL/LAYER.
!
!$$$
!
      USE MACHINE , ONLY : kind_rad
      implicit none
!
      integer L, LP1, IMAX, NSRC, NBD, NVB, NAE, NDM, NXC, NDN
     &,       ICFC, ICWP, NDAY, IR(NDAY), me, ix2
!
      integer     IDM (ix2,L,NAE), IDXC(NXC,IMAX), KPRF(IMAX)
      integer     ISRC(NSRC)
      real (kind=kind_rad)  HAER(NDM,NAE)
     &,                     DZ(IMAX,L), HZ(IMAX,L+1), TAUR(ix2,L,NBD)
     &,                     CMIX(NXC,IMAX), DENN(NDN,IMAX)
! ---  INPUT
      real (kind=kind_rad) S0, CO2
     &, PL (IMAX,LP1), TA(IMAX,L),   WA(IMAX,L),     OA(IMAX,L)
     &, TAUCL(IMAX,L), CCLY(IMAX,L), CFAC(IMAX,LP1), COSZ(IMAX)
     &, ALBUVB(IMAX),  ALBUVD(IMAX), ALBIRB(IMAX),   ALBIRD(IMAX)
     &, RH(IMAX,L),    FICE(IMAX,L)
     &, CWP(IMAX,L),   CIP(IMAX,L),  REW(IMAX,L),    REI(IMAX,L)
 
! ---  OUTPUT
      real (kind=kind_rad)
     &  TUPFXC(IMAX), SUPFXC(IMAX), SDNFXC(IMAX), TDNFLX(IMAX)
     &, TUPFX0(IMAX), SUPFX0(IMAX), SDNFX0(IMAX), HTRC(IMAX,L)
     &, SDNFVB(IMAX), SDNFVD(IMAX), SDNFNB(IMAX), SDNFND(IMAX)
     &, SDN0VB(IMAX), SDN0VD(IMAX), SDN0NB(IMAX), SDN0ND(IMAX)
 
! --- OUTPUT FOR UV-B BAND DOWNWARD SURFACE AND TOP DOWNWARD FLUXES
      real (kind=kind_rad)
     &  SUVBF0(IMAX),     SUVBFC(IMAX)
 
! ---  INTERNAL ARRAY
      real (kind=kind_rad)  HTR0 (IMAX,L)
!
!     Locals
!
      real (kind=kind_rad)
     &  PL1 (NDAY,LP1), TA1(NDAY,L),   WA1(NDAY,L),     OA1(NDAY,L)
     &, TAUCL1(NDAY,L), CCLY1(NDAY,L), CFAC1(NDAY,LP1), COSZ1(NDAY)
     &, AL1UVB(NDAY),   AL1UVD(NDAY),  AL1IRB(NDAY),    AL1IRD(NDAY)
     &, RH1(NDAY,L),    FICE1(NDAY,L)
     &, CWP1(NDAY,L),   CIP1(NDAY,L), REW1(NDAY,L),     REI1(NDAY,L)
 
! --- OUTPUT FOR UV-B BAND DOWNWARD SURFACE AND TOP DOWNWARD FLUXES
      real (kind=kind_rad)
     &  SUV1F0(NDAY),     SUV1FC(NDAY)
 
      real (kind=kind_rad)
     &  TU1FXC(NDAY), SU1FXC(NDAY), SD1FXC(NDAY), TD1FLX(NDAY)
     &, TU1FX0(NDAY), SU1FX0(NDAY), SD1FX0(NDAY), HTRC1(NDAY,L)
     &, SD1FVB(NDAY), SD1FVD(NDAY), SD1FNB(NDAY), SD1FND(NDAY)
     &, SD10VB(NDAY), SD10VD(NDAY), SD10NB(NDAY), SD10ND(NDAY)
     &, HTR01(NDAY,L)
      real (kind=kind_rad) CMIX1(NXC,NDAY), DENN1(NDN,NDAY)
     &,                    DZ1(NDAY,L),     HZ1(NDAY,L+1)
     &,                    TAUR1(nday,L,NBD)
      integer IDXC1(NXC,NDAY), KPRF1(NDAY), idm1(nday,L,NAE)
!
      integer i, ii, k, j
!     logical lprnt
!     integer ipnGlobal,its,mype
!     logical DiagPrint
!     call PhysicsGetIpnItsMype(ipnGlobal,its,mype,DiagPrint)
!
!
!===> ... BEGIN HERE
!
      DO I=1,NDAY
        II         = IR(I)
        TDNFLX(II) = S0 * COSZ(II)
        SDN0VB(I)  = 0.0
        SDN0VD(I)  = 0.0
        SDN0NB(I)  = 0.0
        SDN0ND(I)  = 0.0
      ENDDO
!
!     Reduce the vectors
!
      DO K=1,L
        DO I=1,NDAY
          II          = IR(I)
          PL1(I,K)    = PL(II,K)
          TA1(I,K)    = TA(II,K)
          WA1(I,K)    = WA(II,K)
          OA1(I,K)    = OA(II,K)
          TAUCL1(I,K) = TAUCL(II,K)
          CFAC1(I,K)  = CFAC(II,K)
          CCLY1(I,K)  = CCLY(II,K)
          FICE1(I,K)  = FICE(II,K)
          CWP1(I,K)   = CWP(II,K)
          CIP1(I,K)   = CIP(II,K)
          REW1(I,K)   = REW(II,K)
          REI1(I,K)   = REI(II,K)
          RH1(I,K)    = RH(II,K)
!
          HZ1(I,K)    = HZ(II,K)
          DZ1(I,K)    = DZ(II,K)
        ENDDO
      ENDDO
      DO I=1,NDAY
        II           = IR(I)
        PL1(I,LP1)   = PL(II,LP1)
        CFAC1(I,LP1) = CFAC(II,LP1)
        AL1UVB(I)    = ALBUVB(II)
        AL1UVD(I)    = ALBUVD(II)
        AL1IRB(I)    = ALBIRB(II)
        AL1IRD(I)    = ALBIRD(II)
        COSZ1(I)     = COSZ(II)
        TD1FLX(I)    = TDNFLX(II)
!
        KPRF1(I)     = KPRF(II)
        HZ1(I,LP1)   = HZ(II,LP1)
      ENDDO
      DO K=1,NXC
        DO I=1,NDAY
          II         = IR(I)
          IDXC1(K,I) = IDXC(K,II)
          CMIX1(K,I) = CMIX(K,II)
        ENDDO
      ENDDO
      DO K=1,NDN
        DO I=1,NDAY
          II         = IR(I)
          DENN1(K,I) = DENN(K,II)
        ENDDO
      ENDDO
      DO J=1,NBD
        DO K=1,L
          DO I=1,NDAY
            II           = IR(I)
            TAUR1(I,K,J) = TAUR(II,K,J)
          ENDDO
        ENDDO
      ENDDO
      DO J=1,NAE
        DO K=1,L
          DO I=1,NDAY
            II          = IR(I)
            IDM1(I,K,J) = IDM(II,K,J)
          ENDDO
        ENDDO
      ENDDO
!
      CALL SWR95A(S0,ISRC,PL1,TA1,WA1,OA1,CO2,COSZ1,TAUCL1,
     &            CCLY1,CFAC1,ICFC,ICWP,CWP1,CIP1,REW1,REI1,FICE1,
     &            AL1UVB,AL1UVD,AL1IRB,AL1IRD,KPRF1,IDXC1,CMIX1,DENN1,
     &            RH1,
     &            HTRC1,TU1FXC,TD1FLX,SU1FXC,SD1FXC,
     &            TU1FX0,SU1FX0,SD1FX0,
     &            SD1FVB,SD1FVD,SD1FNB,SD1FND,
! --- FOR UV-B BAND FLUXES
     &            SUV1FC,SUV1F0,
! --- END UV-B
     &            L, LP1, NDAY, NSRC, NBD, NVB, NAE, NDM, NXC, NDN,
     &            HAER, IDM1, DZ1, HZ1, TAUR1, me)
!    &,           lprnt)
!
      DO I=1,NDAY
        II         = IR(I)
        SDNFNB(II) = SD1FNB(I)
        SDNFND(II) = SD1FND(I)
        SDNFVB(II) = SD1FVB(I)
        SDNFVD(II) = SD1FVD(I)
        TUPFX0(II) = TU1FX0(I)
        TUPFXC(II) = TU1FXC(I)
        SUPFX0(II) = SU1FX0(I)
        SUPFXC(II) = SU1FXC(I)
        SDNFX0(II) = SD1FX0(I)
        SDNFXC(II) = SD1FXC(I)

! --- FOR UV-B BAND FLUXES
        SUVBF0(II) = SUV1F0(I)
        SUVBFC(II) = SUV1FC(I)
! --- END UV-B
      ENDDO
!
      DO K=1,L
        DO I=1,NDAY
          HTRC(IR(I),K) = HTRC1(I,K)
        ENDDO
      ENDDO
!
      RETURN
      END
      SUBROUTINE SWR95A(S0,ISRC,PL,TA,WA,OA,CO2,COSZ,TAUCL,
     &            CCLY,CFAC,ICFC,ICWP,CWP,CIP,REW,REI,FICE,
     &            ALBUVB,ALBUVD,ALBIRB,ALBIRD,KPRF,IDXC,CMIX,DENN,RH,
     &            HTRC,TUPFXC,TDNFLX,SUPFXC,SDNFXC,
     &            TUPFX0,SUPFX0,SDNFX0,
     &            SDNFVB,SDNFVD,SDNFNB,SDNFND,
! --- FOR UV-B BAND FLUXES
     &            SUVBFC,SUVBF0,
! --- END UV-B
     &            L, LP1, IMAX, NSRC, NBD, NVB, NAE, NDM, NXC, NDN,
     &            HAER, IDM, DZ, HZ, TAUR, me)
!    &,           lprnt)
!FPP$ NOCONCUR R
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM:    SWR95      COMPUTES SHORT-WAVE RADIATIVE HEATING
!   PROGRAMMER: YU-TAI HOU  ORG: W/NMC20      DATE: 95-02-09
!
! ABSTRACT: THIS CODE IS A MODIFIED VERSION OF M.D. CHOU'S SW
!   RADIATION CODE TO FIT NMC MRF AND CLIMATE MODELS.  IT COMPUTES
!   SW ATMOSPHERIC ABSORPTION AND SCATTERING EFFECTS DUE TO O3,
!   H2O,CO2,O2,CLOUDS, AND AEROSOLS, ETC.
!   IT HAS 8 UV+VIS BANDS AND 3 NIR BANDS (10 K-VALUES EACH).
!
! REFERENCES: CHOU (1986, J. CLIM. APPL.METEOR.)
!   CHOU (1990, J. CLIM.), AND CHOU (1992, J. ATMS. SCI.)
!   CHOU AND SUAREZ (1999, NASA/TM-1999-104606,VOL.15)
!
! PROGRAM HISTORY LOG:
!   94-06-12   M.D. CHOU, GLA.
!   95-02-09   YU-TAI HOU      - RECODE FOR NMC MODELS
!   98-08-03   YU-TAI HOU      - UPDATED CLOUD RADIATIVE PROPERTIES
!            CALCULATION. USE SLINGO'S METHOD (JAS 1989) ON WATER
!            CLOUD, EBERT AND CURRY'S METHOD (JGR 1992) ON ICE CLOUD.
!   99-03-25   YU-TAI HOU      - UPDATED CLOUD RADIATIVE PROPERTY
!            CALCULATIONS USE CHOU ET AL. NEW METHOD (J. CLIM 1998)
!   99-04-27   YU-TAI HOU      - UPDATED CLOUD RADIATIVE PROPERTY
!            CALCULATIONS USE LINEAR T-ADJUSTED METHOD.
!   99-09-13   YU-TAI HOU      - UPDATED TO CHOU'S JUNE,99 VERSION
!
! USAGE:       CALL SWR95
!
! ATTRIBUTES:
!   LANGUAGE:  FORTRAN 77 & Fortran 90
!   MACHINE:   CRAY C-90, IBM SP, SGI
!
! INPUT PARAMETERS:
!   S0     : SOLAR CONSTANT
!   ISRC   : FLAGS FOR SELECTING ABSORBERS
!            1:AEROSOLS, 2:O2, 3:CO2, 4:H2O, 5:O3
!            =0:WITHOUT IT,  =1: WITH IT.
!   PL     : MODEL LEVEL PRESSURE IN MB
!   TA     : MODEL LAYER TEMPERATURE IN K
!   WA     : LAYER SPECIFIC HUMIDITY IN GM/GM
!   OA     : LAYER OZONE CONCENTRATION IN GM/GM
!   CO2    : CO2 MIXING RATION BY VOLUME
!   COSZ   : COSINE OF SOLAR ZENITH ANGLE
!   TAUCL  : OPTICAL DEPTH OF CLOUD LAYERS
!   CCLY   : LAYER CLOUD FRACTION
!   CFAC   : FRACTION OF CLEAR SKY VIEW AT THE LAYER INTERFACE
!   ICFC   : =0 NO CLOUD FACTOR TO WEIGH CLEAR AND CLOUDY FLUXES
!            =1 USE CLOUD FACTOR TO WEIGH CLEAR AND CLOUDY FLUXES
!   ICWP   : FLAG INDICATES THE METHOD USED FOR CLOUD PROPERTIES
!            CALCULATIONS, =0 USE T-P; =1 USE CWC/CIC.
!   CWP    : LAYER CLOUD WATER PATH (G/M**2)
!   CIP    : LAYER CLOUD ICE PATH (G/M**2)
!   REW    : LAYER WATER CLOUD DROP EFFECTIVE RADIUS (MICRON)
!   REI    : LAYER ICE CLOUD DROP EFFECTIVE RADIUS
!   FICE   : FRACTION OF CLOUD ICE CONTENT
!   ALBUVB : UV+VIS SURF DIRECT ALBEDO
!   ALBUVD : UV+VIS SURF DIFFUSED ALBEDO
!   ALBIRB : NIR SURF DIRECT ALBEDO
!   ALBIRD : NIR SURF DIFFUSED ALBEDO
!   PAER   : AEROSOL PROFILES (FRACTION)
!
! OUTPUT PARAMETER:
!   HTRC   : HEATING RATES FOR CLOUDY SKY IN  K/DAY
!   TUPFXC : UPWARD FLUX AT TOA FOR CLOUDY SKY  W/M**2
!   TDNFLX : DNWARD FLUX AT TOA FOR ALL SKY  W/M**2
!   SUPFXC : UPWARD FLUX AT SFC FOR CLOUDY SKY  W/M**2
!   SDNFXC : DNWARD FLUX AT SFC FOR CLOUDY SKY  W/M**2
!   TUPFX0 : UPWARD FLUX AT TOA FOR CLEAR SKY   W/M**2
!   SUPFX0 : UPWARD FLUX AT SFC FOR CLEAR SKY   W/M**2
!   SDNFX0 : DNWARD FLUX AT SFC FOR CLEAR SKY   W/M**2
!   SDNFVB : DOWNWARD SURFACE VIS BEAM FLUX     W/M**2
!   SDNFNB : DOWNWARD SURFACE NIR BEAM FLUX     W/M**2
!   SDNFVD : DOWNWARD SURFACE VIS DIFF FLUX     W/M**2
!   SDNFND : DOWNWARD SURFACE NIR DIFF FLUX     W/M**2
!
! NOTE:
!   FOR ALL QUANTITIES, K=1 IS THE TOP LEVEL/LAYER, EXCEPT
!   SI AND SL, FOR WHICH K=1 IS THE SURFACE LEVEL/LAYER.
!
!$$$
!
!
      USE MACHINE , ONLY : kind_rad,kind_phys
      implicit none
!
      integer L, LP1, IMAX, NSRC, NBD, NVB, NAE, NDM, NXC, NDN
     &,       ICFC, ICWP, me
      integer IDM (IMAX,L,NAE), IDXC(NXC,IMAX), KPRF(IMAX), ISRC(NSRC)
      real (kind=kind_rad) HAER(NDM,NAE)
     &,                     DZ(IMAX,L), HZ(IMAX,L+1), TAUR(IMAX,L,NBD)
     &,                     CMIX(NXC,IMAX), DENN(NDN,IMAX)
!
! ---  INPUT
      real (kind=kind_rad) S0, CO2
     &, PL (IMAX,LP1), TA(IMAX,L),   WA(IMAX,L),    OA(IMAX,L)
     &, TAUCL(IMAX,L), CCLY(IMAX,L), CFAC(IMAX,LP1),COSZ(IMAX)
     &, ALBUVB(IMAX),  ALBUVD(IMAX), ALBIRB(IMAX),  ALBIRD(IMAX)
     &, RH(IMAX,L),    FICE(IMAX,L)
     &, CWP(IMAX,L),   CIP(IMAX,L),  REW(IMAX,L),   REI(IMAX,L)
!
 
! ---  OUTPUT
      real (kind=kind_rad)
     &  TUPFXC(IMAX), SUPFXC(IMAX), SDNFXC(IMAX), TDNFLX(IMAX)
     &, TUPFX0(IMAX), SUPFX0(IMAX), SDNFX0(IMAX), HTRC(IMAX,L)
     &, SDNFVB(IMAX), SDNFVD(IMAX), SDNFNB(IMAX), SDNFND(IMAX)
     &, SDN0VB(IMAX), SDN0VD(IMAX), SDN0NB(IMAX), SDN0ND(IMAX)
 
! --- OUTPUT FOR UV-B BAND DOWNWARD SURFACE AND TOP DOWNWARD FLUXES
      real (kind=kind_rad)
     &  SUVBFC(IMAX),     SUVBF0(IMAX)
 
! ---  INTERNAL ARRAY
      real (kind=kind_rad)
     &  FNET0(IMAX,LP1), FNETC(IMAX,LP1), HTR0 (IMAX,LP1)
     &, DFLX0(IMAX,LP1), DFLXC(IMAX,LP1), DP   (IMAX,L)
     &, SCAL (IMAX,L),   SWH  (IMAX,LP1), SO2  (IMAX,LP1)
     &, WH   (IMAX,L),   OH   (IMAX,L),   SWU  (IMAX,LP1)
     &, CF0  (IMAX),     CF1  (IMAX),     SNT  (IMAX)
     &, CNT  (IMAX)
      real (kind=kind_rad) rewi(imax,L), reii(imax,L)
!     logical lprnt
!
      real (kind=kind_rad) ZTHIK(IMAX,L), CSMIK(IMAX,L)
!
      real (kind=kind_rad) taucrt
      integer IFPR, IBND
      DATA TAUCRT / 0.05 /, IFPR / 0 /
      DATA IBND / 1 /        !===> ... IBND=1:USE ONE NIR BAND
!     DATA IBND / 2 /        !===> ... IBND=2:USE TWO NIR BANDS
c$$$      SAVE TAUCRT, IFPR, IBND
!
      real (kind=kind_rad) tfac, tem, rcf1, ccc, xa, to2
     &,                     u1, du, w1, dw, fac
      integer i, k, jtop
!
      include 'co2tab_sw.h'
!
cc
cc--------------------------------------------------------------------
cc
      real(kind=kind_rad) cons_1pdm11          !constant
      integer ipnGlobal,its,mype
      logical DiagPrint
cc
      cons_1pdm11     =        1.d-11          !constant
cc
cc--------------------------------------------------------------------
cc
!===> ... BEGIN HERE
!     call PhysicsGetIpnItsMype(ipnGlobal,its,mype,DiagPrint)

! This AEROSOL print shuold be uncommented and printed when PrintDiags is true - JFM
!     IF (IFPR .EQ. 0) THEN
!       if (me.eq.0) WRITE(6,12) (ISRC(I),I=1,NSRC)
! 12    FORMAT(3X,'AEROSOL, O2, CO2, H2O, O3 =',5I3)
!       IFPR = 1
!     END IF
!
      DFLXC(:,:) = 0.0
      DO I=1,IMAX
        SWH (I,1) = 0.0
        SO2 (I,1) = 0.0
        TUPFXC(I) = 0.0
        TUPFX0(I) = 0.0
        SUPFXC(I) = 0.0
        SUPFX0(I) = 0.0
        SDNFXC(I) = 0.0
        SDNFX0(I) = 0.0
        CF0(I)    = CFAC(I,LP1)
        CF1(I)    = 1.0 - CF0(I)
        SNT(I)    = 1.0 / COSZ(I) ! SNT = SECANT OF SOLAR ZENITH ANGLE
!
        SDNFVB(I) = 0.0
        SDNFVD(I) = 0.0
        SDNFNB(I) = 0.0
        SDNFND(I) = 0.0
        SDN0VB(I) = 0.0
        SDN0VD(I) = 0.0
        SDN0NB(I) = 0.0
        SDN0ND(I) = 0.0
 
! --- FOR UV-B FLUXES
        SUVBFC(I) = 0.0
        SUVBF0(I) = 0.0
! --- END UV-B
 
      ENDDO
!
      TFAC = 0.5 / 300.0
      DO K=1,L
        DO I=1,IMAX
!===> ... LAYER THICKNESS AND PRESSURE SCALING FUNCTION FOR
!         WATER VAPOR ABSORPTION
          DP  (I,K) = PL(I,K+1) - PL(I,K)
          SCAL(I,K) = DP(I,K) * (TFAC*(PL(I,K)+PL(I,K+1)))**0.8
!===> ... SCALED ABSORBER AMOUNTS FOR H2O(WH,SWH), UNIT : G/CM**2
          TEM     = 0.00135*(TA(I,K)-240.0)
          WH(I,K) = 1.02 * WA(I,K) * SCAL(I,K)
!    &            * EXP(0.00135*(TA(I,K)-240.0))
     &            * (1.0 + TEM + 0.5*TEM*TEM) + 1.0E-11
        ENDDO
      ENDDO
      DO K=1,L
        DO I=1,IMAX
          SWH(I,K+1) = SWH(I,K) + WH(I,K)
          ZTHIK(I,K) = COSZ(I)
          CSMIK(I,K) = SNT(I)
        ENDDO
      ENDDO
!
!===> ... INITIALIZE FLUXES
!
      DO K=1,LP1
        DO I=1,IMAX
          FNET0(I,K) = 0.0
          FNETC(I,K) = 0.0
          DFLX0(I,K) = 0.0
        ENDDO
      ENDDO
!
      IF (ICFC .EQ. 1) THEN
        DO I=1,IMAX
          CFAC(I,LP1) = 0.0
        END DO
        DO K=1,L
          DO I=1,IMAX
            IF (CF1(I) .GT. 0.0) THEN
              RCF1 = 1.0 / CF1(I)
              CFAC(I,K) = (CFAC(I,K) - CF0(I)) * RCF1
              CCLY(I,K) = CCLY(I,K) * RCF1
            END IF
          END DO
        END DO
      END IF
!
      IF (ICWP.NE. 1) THEN
        DO K=1,L
          DO I=1,IMAX
!0900       TAUCL(I,K) = TAUCL(I,K) * CCLY(I,K)
            TAUCL(I,K) = TAUCL(I,K) * CCLY(I,K) * SQRT(CCLY(I,K))
          END DO
        END DO
      ELSE
        DO K=1,L
          DO I=1,IMAX
!0799       CCC = CCLY(I,K) * SQRT(CCLY(I,K))
            CCC = CCLY(I,K)
            CWP(I,K) = CWP(I,K) * CCC
            CIP(I,K) = CIP(I,K) * CCC
            REWI(I,K) = 1.0 / REW(I,K)
            REII(I,K) = 1.0 / REI(I,K)
          END DO
        END DO
      END IF
!
      IF (ISRC(4) .EQ. 1) THEN     !===> ... COMPUTE NIR FLUXES
!                                            ------------------
        CALL SOLIR(WH,TA,TAUCL,CSMIK,ZTHIK,IBND,FICE,
     &             ISRC(1),KPRF,IDXC,CMIX,DENN,RH,ALBIRB,ALBIRD,
     &             ICWP,CWP,CIP,CCLY,REW,REI,REWI,REII,
     &             TUPFXC,SUPFXC,SDNFXC,TUPFX0,SUPFX0,SDNFX0,
     &             FNET0,FNETC,SDN0NB,SDN0ND,SDNFNB,SDNFND,
     &             L, LP1, IMAX, NBD, NVB, NAE, NDM, NXC, NDN,
     &             HAER, IDM, DZ, HZ, TAUR)
!    &,            lprnt)
      END IF
!
      IF (ISRC(5) .EQ. 1) THEN     !===> ... COMPUTE UV+VISIBLE FLUXES
!                                            -------------------------
!               SCALED AMOUNTS FOR O3(WH), UNIT : (CM-AMT)STP FOR O3.
        XA = 1.02 * 466.7
        DO K=1,L
          DO I=1,IMAX
            OH(I,K) = XA * OA(I,K) * DP(I,K) + 1.0E-11
          ENDDO
        ENDDO
!
        CALL SOLUV(WH,OH,TA,TAUCL,CSMIK,ZTHIK,FICE,
     &             ISRC(1),KPRF,IDXC,CMIX,DENN,RH,ALBUVB,ALBUVD,
     &             ICWP,CWP,CIP,CCLY,REW,REI,REWI,REII,
     &             TUPFXC,SUPFXC,SDNFXC,TUPFX0,SUPFX0,SDNFX0,
     &             FNET0,FNETC,SDN0VB,SDN0VD,SDNFVB,SDNFVD,
! --- FOR UV-B BAND FLUXES
     &             SUVBFC,SUVBF0,
! --- END UV-B
     &             L, LP1, IMAX, NBD, NVB, NAE, NDM, NXC, NDN,
     &             HAER, IDM, DZ, HZ, TAUR)
!    &,            lprnt)

      END IF
!
!===> ... COMPUTE THE ABSORPTION DUE TO OXYGEN,CHOU(1990,J.CLIMATE,209-217)
!         PRESSURE SCALED AMOUNTS FOR O2(O2,SO2), UNIT IS (CM-ATM)STP FOR O2.
!         THE CONSTANT 165.22=(1000/980)*23.14%*(22400/32)
!
      IF (ISRC(2) .EQ. 1) THEN
        DO I=1,IMAX
          CNT(I) = 165.22 * SNT(I)
        END DO
        DO K=1,L
          DO I=1,IMAX
            SO2(I,K+1) = SO2(I,K) + CNT(I) * SCAL(I,K)
          ENDDO
        ENDDO
!===> ... COMPUTE FLUX REDUCTION DUE TO OXYGEN, THE CONSTANT 0.0633 IS
!         THE FRACTION OF INSOLATION CONTAINED IN THE OXYGEN BANDS.
!         TO2 IS THE BROADBAND TRANSMISSION FUNCTION FOR OXYGEN
        DO K=2,LP1
          DO I=1,IMAX
            TO2        = EXP(-0.145E-3 * SQRT(SO2(I,K)) )
            DFLX0(I,K) = 0.0633 * (1.0 - TO2)
          ENDDO
        ENDDO
      END IF
!
!===> ... TABLE LOOK-UP FOR THE ABSORPTION DUE TO CO2
!         COMPUTE SCALED AMOUNTS FOR CO2(WC,SO2).
!         THE CONSTANT 789=(1000/980)*(44/28.97)*(22400/44)
!
      IF (ISRC(3) .EQ. 1) THEN
        DO I=1,IMAX
          CNT(I)   = CO2 * SNT(I)
          SO2(I,1) = MAX(SO2(I,1), cons_1pdm11)     !constant
        END DO
        DO K=1,L
          DO I=1,IMAX
            SO2(I,K+1) = SO2(I,K) + 789.0 * CNT(I)*SCAL(I,K)
          ENDDO
        ENDDO
!
!===> ... FOR CO2 ABSORPTION IN SPECTRUM 1.220-2.270 MICRON
!         BOTH WATER VAPOR AND CO2 ABSORPTIONS ARE MODERATE
!         SO2 AND SWH ARE THE CO2 AND WATER VAPOR AMOUNTS
!         INTEGRATED FROM THE TOP OF THE ATMOSPHERE
!
        U1 = -3.0
        DU = 0.15
        W1 = -4.0
        DW = 0.15
        DO K=2,LP1
          DO I=1,IMAX
            SWU(I,K) = LOG10(SO2(I,K))
            SWH(I,K) = LOG10(SWH(I,K)*SNT(I))
          END DO
        END DO
!
!===> ... DFLX0 IS THE UPDATED FLUX REDUCTION
!
        CALL FLXCO2(SWU,U1,DU,NU,SWH,W1,DW,NW,CAH,DFLX0
     &,                                           L, IMAX)
!
!===> ... FOR CO2 ABSORPTION IN SPECTRUM 2.270-10.00 MICRON
!         WHERE THE CO2 ABSORPTION HAS A LARGE IMPACT ON THE
!         HEATING OF MIDDLE ATMOSPHERE
!
        U1 = 0.250E-3
        DU = 0.050E-3
        W1 = -2.0
        DW = 0.05
!===> ... CO2 MIXING RATIO IS INDEPENDENT OF SPACE
!         SWH IS THE LOGARITHM OF PRESSURE
        DO K=2,LP1
          DO I=1,IMAX
            SWU(I,K) = CNT(I)
            SWH(I,K) = LOG10(PL(I,K))
          END DO
        END DO
!===> ... DFLX0 IS THE UPDATED FLUX REDUCTION
!
        CALL FLXCO2(SWU,U1,DU,NX,SWH,W1,DW,NY,COA,DFLX0
     &,                                           L, IMAX)
!
      ENDIF
!
!===> ... ADJUST FOR THE EFFECT OF O2 AND CO2 ON CLEAR SKY NET FLUX
!
      IF (ISRC(2).EQ.1 .OR. ISRC(3).EQ.1) THEN
!       DO K=1,LP1
!         DO I=1,IMAX
!           FNET0(I,K) = FNET0(I,K) - DFLX0(I,K)
!         ENDDO
!       ENDDO
!
!===> ... ADJUST FOR THE EFFECT OF O2 AND CO2 ON CLOUD SKY NET FLUX
!
        DO I=1,IMAX
          JTOP = LP1
!===> ... ABOVE CLOUDS
          DO K=1,LP1
            DFLXC(I,K) = DFLX0(I,K)
            IF (CFAC(I,K) .LT. 1.0) THEN
              JTOP = K
              EXIT
            END IF
          END DO
!===> ... BELOW CLOUD TOP
          IF (JTOP .LT. LP1) THEN
            DO K=JTOP+1,LP1
              DFLXC(I,K) = DFLX0(I,K) * (FNETC(I,K)/FNET0(I,K))
            END DO
          END IF
          DO K=1,LP1
            FNET0(I,K) = FNET0(I,K) - DFLX0(I,K)
            FNETC(I,K) = FNETC(I,K) - DFLXC(I,K)
          END DO
        ENDDO
!
!===> ... ADJUST FOR OTHER FLUXES
!
        DO I=1,IMAX
          SDNFX0(I) = SDNFX0(I) - DFLX0(I,LP1)
          SDNFXC(I) = SDNFXC(I) - DFLXC(I,LP1)
          SDN0NB(I) = SDN0NB(I) - DFLX0(I,LP1)
          SDNFNB(I) = SDNFNB(I) - DFLXC(I,LP1)
        ENDDO
      END IF
!
      IF (ICFC .EQ. 1) THEN
!===> ... COMPUTE FINAL FLUXES AT TOP AND SURFACE
        DO I=1,IMAX
          SDNFVB(I) = CF0(I)*SDN0VB(I) + CF1(I)*SDNFVB(I)
          SDNFVD(I) = CF0(I)*SDN0VD(I) + CF1(I)*SDNFVD(I)
          SDNFNB(I) = CF0(I)*SDN0NB(I) + CF1(I)*SDNFNB(I)
          SDNFND(I) = CF0(I)*SDN0ND(I) + CF1(I)*SDNFND(I)
          TUPFXC(I) = CF0(I)*TUPFX0(I) + CF1(I)*TUPFXC(I)
          SUPFXC(I) = CF0(I)*SUPFX0(I) + CF1(I)*SUPFXC(I)
          SDNFXC(I) = CF0(I)*SDNFX0(I) + CF1(I)*SDNFXC(I)
! --- FOR UV-B FLUX
          SUVBFC(I) = CF0(I)*SUVBF0(I) + CF1(I)*SUVBFC(I)
! --- END UV-B
        ENDDO
        DO K=1,LP1
          DO I=1,IMAX
            FNETC (I,K) = CF0(I)*FNET0(I,K) + CF1(I)*FNETC(I,K)
          ENDDO
        ENDDO
      END IF
!
!===> ... CONVERT FLUX UNIT TO W/M**2
!
      DO K=1,LP1
        DO  I=1,IMAX
!CLEAR    FNET0 (I,K) = FNET0(I,K) * TDNFLX(I)
          FNETC (I,K) = FNETC(I,K) * TDNFLX(I)
        ENDDO
      ENDDO
      DO I=1,IMAX
        SDNFNB(I) = SDNFNB(I) * TDNFLX(I)
        SDNFND(I) = SDNFND(I) * TDNFLX(I)
        SDNFVB(I) = SDNFVB(I) * TDNFLX(I)
        SDNFVD(I) = SDNFVD(I) * TDNFLX(I)
        TUPFX0(I) = TUPFX0(I) * TDNFLX(I)
        TUPFXC(I) = TUPFXC(I) * TDNFLX(I)
        SUPFX0(I) = SUPFX0(I) * TDNFLX(I)
        SUPFXC(I) = SUPFXC(I) * TDNFLX(I)
        SDNFX0(I) = SDNFX0(I) * TDNFLX(I)
        SDNFXC(I) = SDNFXC(I) * TDNFLX(I)
 
! --- FOR UV-B BAND FLUXES
        SUVBF0(I) = SUVBF0(I) * TDNFLX(I)
        SUVBFC(I) = SUVBFC(I) * TDNFLX(I)
! --- END UV-B
 
      ENDDO
!
!===> ... FAC IS THE FACTOR FOR HEATING RATES (IN K/DAY)
!         IF USE K/SEC, RESULT SHOULD BE DEVIDED BY 86400.
!
!     FAC = 3.6*24./10.031*.98
      FAC = 8.4410328
!
      DO K=1,L
        DO I=1,IMAX
!CLEAR    HTR0(I,K) = (FNET0(I,K)-FNET0(I,K+1)) * FAC / DP(I,K)
          HTRC(I,K) = (FNETC(I,K)-FNETC(I,K+1)) * FAC / DP(I,K)
        ENDDO
      ENDDO
!
      RETURN
      END
      SUBROUTINE SOLUV(WZ,OZ,TA,TAUCL,CSM,ZTH,FICE,
     &                 KAER,KPRF,IDXC,CMIX,DENN,RH,ALBB,ALBD,
     &                 ICWP,CWP,CIP,CCLY,REW,REI,REWI,REII,
     &                 TUPFXC,SUPFXC,SDNFXC,TUPFX0,SUPFX0,SDNFX0,
     &                 FNET0,FNETC,DWSFB0,DWSFD0,DWSFBC,DWSFDC,
! --- FOR UV-B BAND FLUXES
     &                 SUVBFC,SUVBF0,
! --- END UV-B
     &                 L,LP1,IMAX,NBD,NVB,NAE,NDM,NXC,NDN,
     &                 HAER,IDM,DZ,HZ,TAUR)
!    &,                lprnt)
!FPP$ NOCONCUR R
!*******************************************************************
!  COMPUTE SOLAR FLUX IN THE UV+VISIBLE REGION
!  THE UV+VISIBLE REGION IS GROUPED INTO 8 BANDS:
!    UV-C     (.175-.225);(.225-.245,.260-.280);(.245-.260);
!    UV-B     (.280-.295);(.295-.310);(.310-.320);
!    UV-A     (.320-.400);
!    PAR      (.400-.700)
!
!  INPUT PARAMETERS:                            UNITS
!    WZ,OZ,TA,TAUCL,CSM,FICE,KAER,PAER,ALBB,ALBD
!    ICWP,CWP,CIP,CCLV,REW,REI
!
!  OUTPUT PARAMETERS:
!    FNET0  : CLEAR SKY NET FLUX
!    FNETC  : CLOUDY SKY NET FLUX
!    TUPFXC : CLOUDY SKY UPWARD FLUX AT TOA
!    SUPFXC : CLOUDY SKY UPWARD FLUX AT SFC
!    SDNFXC : CLOUDY SKY DOWNWARD FLUX AT SFC
!    TUPFX0 : CLEAR SKY UPWARD FLUX AT TOA
!    SUPFX0 : CLEAR SKY UPWARD FLUX AT SFC
!    SDNFX0 : CLEAR SKY DOWNWARD FLUX AT SFC
!    DWSFB0 : CLEAR SKY SFC DOWN DIR. FLUX
!    DWSFD0 : CLEAR SKY SFC DOWN DIF. FLUX
!    DWSFBC : CLOUDY SKY SFC DOWN DIR. FLUX
!    DWSFDC : CLOUDY SKY SFC DOWN DIF. FLUX
!
!  FIXED INPUT DATA:
!    FRACTION OF SOLAR FLUX CONTAINED
!       IN THE 8 BANDS (SS)                     FRACTION
!    RAYLEIGH OPTICAL THICKNESS (TAURAY)        /MB
!    OZONE ABSORPTION COEFFICIENT (AK)          /(CM-ATM)STP
!
!  THE FOLLOWING PARAMETERS MUST BE SPECIFIED BY USERS:
!    CLOUD ASYMMETRY FACTOR (ASYCL)             N/D
!  AEROSOL PARAMETERS ARE FROM SUBPROGRAM AEROS:
!
!  PROGRAM HISTORY LOG:
!   94-06-12   M.D. CHOU, GLA.
!   95-02-09   YU-TAI HOU      - RECODE FOR NMC MODELS
!   98-08-03   YU-TAI HOU      - UPDATED CLOUD RADIATIVE PROPERTIES
!            CALCULATION. USE SLINGO'S METHOD (JAS 1989) ON WATER
!            CLOUD, EBERT AND CURRY'S METHOD (JGR 1992) ON ICE CLOUD.
!   99-03-25   YU-TAI HOU      - UPDATED CLOUD PROPERTIES USE THE
!            MOST RECENT CHOU ET AL. DATA (J. CLIM 1998)
!   99-04-27   YU-TAI HOU      - UPDATED CLOUD RADIATIVE PROPERTY
!            CALCULATIONS USE LINEAR T-ADJUSTED METHOD.
!   99-09-13   YU-TAI HOU      - UPDATED TO CHOU'S JUNE,1999 VERSION
!
!********************************************************************
!
!
      USE MACHINE , ONLY : kind_rad
      implicit none
!
      integer nvbb
!     PARAMETER (NVBB=4)
      PARAMETER (NVBB=8)
!
      integer L, LP1, IMAX, NBD, NVB, NAE, NDM, NXC, NDN
     &,       KAER, ICWP
!
      integer     IDM (IMAX,L,NAE), IDXC(NXC,IMAX), KPRF(IMAX)
      real (kind=kind_rad) HAER(NDM,NAE)
     &,                     DZ(IMAX,L), HZ(IMAX,L+1), TAUR(IMAX,L,NBD)
     &,                     CMIX(NXC,IMAX), DENN(NDN,IMAX)
!
! --- INPUT
      real (kind=kind_rad)
     &  OZ(IMAX,L),   TAUCL(IMAX,L), ALBB(IMAX), ALBD(IMAX)
     &, CSM(IMAX,L),  ZTH(IMAX,L),   RH(IMAX,L)
     &, TA(IMAX,L),   FICE(IMAX,L)
     &, CWP(IMAX,L),  CIP(IMAX,L),   REW(IMAX,L), REI(IMAX,L)
     &, CCLY(IMAX,L), WZ(IMAX,L)
     &, REWI(IMAX,L), REII(IMAX,L)
! --- OUTPUT
      real (kind=kind_rad)
     &  FNET0 (IMAX,LP1), DWSFB0(IMAX), DWSFD0(IMAX)
     &, FNETC (IMAX,LP1), DWSFBC(IMAX), DWSFDC(IMAX)
     &, TUPFXC(IMAX),     SUPFXC(IMAX), SDNFXC(IMAX)
     &, TUPFX0(IMAX),     SUPFX0(IMAX), SDNFX0(IMAX)
! --- OUTPUT FOR UV-B BAND DOWNWARD SURFACE FLUXES
      real (kind=kind_rad)
     &  SUVBFC(IMAX),     SUVBF0(IMAX)
! --- TEMPORARY ARRAY
      real (kind=kind_rad)
     &  UPFLUX(IMAX,LP1), DWFLUX(IMAX,LP1)
     &, DWSFXB(IMAX),     DWSFXD(IMAX)
     &, TAUTO (IMAX,L),   SSATO (IMAX,L),   ASYTO (IMAX,L)
     &, TAURS (IMAX,L),   SSAT1 (IMAX,L),   ASYT1 (IMAX,L)
     &, TAUAER(IMAX,L),   SSAAER(IMAX,L),   ASYAER(IMAX,L)
     &, FFFCW (IMAX,L),   FFFT1 (IMAX,L),   FFFTO (IMAX,L)
     &, ASYCW (IMAX,L),   SSACW (IMAX,L)
! --- SOLAR FLUX AND ABSORPTION COEFFICIENTS
     &, SS(NVBB),         AK(NVBB),         WK(NVBB)
!0499
! --- T ADJUSTED CLD PROPERTY METHOD
     &, A0W(2), A1W(2), B0W(2), B1W(2), B0I(2), B1I(2), B2I(2)
     &, A0I(2), A1I(2), C0W(2), C1W(2), C0I(2), C1I(2), C2I(2)
     &, SSAW0(2), SSAI0(2), ASYW0(2), ASYI0(2)
     &, FACW(IMAX,L), FACI(IMAX,L)
     &, FFFRS0,       FPMIN, FPMAX
!
      logical cloudy(imax)
      Integer ncloud
!     logical lprnt
!
!
      DATA SS / 0.00057, 0.00367, 0.00083, 0.00417,
     &          0.00600, 0.00556, 0.05913, 0.39081 /
      DATA AK / 30.47, 187.2, 301.9, 42.83,
     &          7.090, 1.250, .0345, .0572 /
      DATA WK / 7*0.0E0, 0.75E-3 /
      DATA SSAW0 /.999998,.999998/, SSAI0 /.999994,.999995/
     &     ASYW0 / 0.853,  0.853 /, ASYI0 / 0.7991, 0.7998/
     &,    FFFRS0 / 0.1 /
      DATA FPMIN, FPMAX / 1.0E-8, 0.999999 /
!0898 - COEFF FOR WATER CLOUD
                         D A T A
! --- T ADJUSTED WATER/ICE CLOUD COEFF.
     &   A0W / 0.2807E-1,0.2798E-1 /, A1W / 0.1307E+1,0.1309E+1 /
     &,  B0W / -.1176E-6,-.1810E-6 /, C0W / 0.8276E+0,0.8272E+0 /
     &,  B1W / 0.1770E-6,0.1778E-6 /, C1W / 0.2541E-2,0.2565E-2 /
     &,  A0I / -.3011E-4,-.5975E-5 /, A1I / 0.2519E+1,0.2517E+1 /
     &,  B0I / 0.1688E-6,0.1721E-6 /, C0I / 0.7473E+0,0.7480E+0 /
     &,  B1I / 0.9936E-7,0.9177E-7 /, C1I / 0.1015E-2,0.1015E-2 /
     &,  B2I /-.1114E-10,-.1125E-10/, C2I / -.2524E-5,-.2531E-5 /
!
c$$$      SAVE SS, AK, WK, ASYW0, ASYI0, SSAW0, SSAI0, FFFRS0, FPMIN, FPMAX
c$$$      SAVE A0W,A1W,B0W,B1W,C0W,C1W,
c$$$     &     A0I,A1I,B0I,B1I,C0I,C1I,B2I,C2I
!
      real (kind=kind_rad) tau1,  tau2,  ssa1,  ssa2,  asy1,  asy2
     &,                     ssaw1, ssaw2, asyw1, asyw2, tauoz, tauwv
     &,                     tem
      integer i, k, iv
cc
cc--------------------------------------------------------------------
cc
      real(kind=kind_rad) cons_0               !constant
      real(kind=kind_rad) cons_10              !constant
      real(kind=kind_rad) cons_30              !constant
!     integer ipnGlobal,its,mype
!     logical DiagPrint
!     call PhysicsGetIpnItsMype(ipnGlobal,its,mype,DiagPrint)

cc
      cons_0          =         0.d0           !constant
      cons_10         =        10.d0           !constant
      cons_30         =        30.d0           !constant
cc
cc--------------------------------------------------------------------
cc
!
      DO K=1,L
        DO I=1,IMAX
          FACW(I,K) = MAX(cons_0, MIN(cons_10,273.15-TA(I,K)))*0.1     !constant
          FACI(I,K) = MAX(cons_0, MIN(cons_30,263.15-TA(I,K)))/30.0     !constant
        ENDDO
      ENDDO
      cloudy(:) = .false.
!
      IF (NVB .NE. NVBB) THEN
         PRINT *,' NVB=',NVB,' NVBB=',NVBB,' RUN STOPPED'
         STOP
      ENDIF
!
      IF (ICWP .NE. 1) THEN
        DO K=1,L
        DO I=1,IMAX
          IF (TAUCL(I,K) .GT. 0.0) THEN
            TAU2 = FICE(I,K) * TAUCL(I,K)
            TAU1 = TAUCL(I,K) - TAU2
C0499 - T-ADJ PROP FROM SPECIFIED SSA AND ASY
            SSA1 = FACW(I,K)*SSAW0(1) + (1.0-FACW(I,K))*SSAW0(2)
            SSA2 = FACI(I,K)*SSAI0(1) + (1.0-FACI(I,K))*SSAI0(2)
            SSAW1 = SSA1 * TAU1
            SSAW2 = SSA2 * TAU2
!           SSA1 = (1.0-FICE(I,K))*(FACW(I,K) *SSAW0(1)
!    &                       + (1.0-FACW(I,K))*SSAW0(2) )
!           SSA2 =        FICE(I,K) *(FACI(I,K) *SSAI0(1)
!    &                       + (1.0-FACI(I,K))*SSAI0(2) )
!           SSAW1 = SSA1 * TAUCL(I,K)
!           SSAW2 = SSA2 * TAUCL(I,K)
            SSACW(I,K) = SSAW1 + SSAW2
!           ASY1 = (1.0-FICE(I,K))*(FACW(I,K) *ASYW0(1)
!    &                       + (1.0-FACW(I,K))*ASYW0(2) )
!           ASY2 =        FICE(I,K) *(FACI(I,K) *ASYI0(1)
!    &                       + (1.0-FACI(I,K))*ASYI0(2) )
            ASY1 = FACW(I,K)*ASYW0(1) + (1.0-FACW(I,K))*ASYW0(2)
            ASY2 = FACI(I,K)*ASYI0(1) + (1.0-FACI(I,K))*ASYI0(2)
            ASYW1 = ASY1 * SSAW1
            ASYW2 = ASY2 * SSAW2
            ASYCW(I,K) = ASYW1 + ASYW2
            FFFCW(I,K) = ASY1*ASYW1 + ASY2*ASYW2
            cloudy(i) = .true.
          ELSE
            SSACW(I,K) = 1.0
            ASYCW(I,K) = 0.0
            FFFCW(I,K) = 0.0
          END IF
        ENDDO
        ENDDO
      ELSE
        DO K=1,L
        DO I=1,IMAX
          IF (CCLY(I,K) .GT. 0.01) THEN
!0499 --- T-ADJ PROP FROM ICE/WATER PATHS
            TAU1 = CWP(I,K)*(   FACW(I,K) *(A0W(1)+A1W(1)*REWI(I,K))
     &                     +(1.-FACW(I,K))*(A0W(2)+A1W(2)*REWI(I,K)))
            TAU2 = CIP(I,K)*(   FACI(I,K) *(A0I(1)+A1I(1)*REII(I,K))
     &                     +(1.-FACI(I,K))*(A0I(2)+A1I(2)*REII(I,K)))
            TAUCL(I,K) = TAU1 + TAU2
            SSA1 = 1.0 - (   FACW(I,K) *(B0W(1)+B1W(1)*REW(I,K))
     &                 + (1.-FACW(I,K))*(B0W(2)+B1W(2)*REW(I,K)) )
            SSA2 = 1.0 - (  FACI(I,K) *(B0I(1)
     &                 + (B1I(1)+B2I(1)*REI(I,K))*REI(I,K))
     &                 + (1.-FACI(I,K))*(B0I(2)
     &                 + (B1I(2)+B2I(2)*REI(I,K))*REI(I,K)) )
            SSAW1 = SSA1 * TAU1
            SSAW2 = SSA2 * TAU2
            SSACW(I,K) = SSAW1 + SSAW2
            ASY1 =     FACW(I,K) *(C0W(1)+C1W(1)*REW(I,K))
     &           + (1.-FACW(I,K))*(C0W(2)+C1W(2)*REW(I,K))
            ASY2 =     FACI(I,K) *(C0I(1)
     &           + (C1I(1)+C2I(1)*REI(I,K))*REI(I,K) )
     &           + (1.-FACI(I,K))*(C0I(2)
     &           + (C1I(2)+C2I(2)*REI(I,K))*REI(I,K) )
            ASYW1 = ASY1 * SSAW1
            ASYW2 = ASY2 * SSAW2
            ASYCW(I,K) = ASYW1 + ASYW2
            FFFCW(I,K) = ASY1*ASYW1 + ASY2*ASYW2
            cloudy(i)  = .true.
          ELSE
            TAUCL(I,K) = 0.0
            SSACW(I,K) = 1.0
            ASYCW(I,K) = 0.0
            FFFCW(I,K) = 0.0
          END IF
        ENDDO
        ENDDO
      END IF
!
      NCLOUD = 0
      DO I=1,IMAX
        if (cloudy(i)) ncloud = ncloud + 1
      ENDDO
!
!===> ... INTEGRATION OVER SPECTRAL BANDS
!
      DO IV=1,NVB
!
!===> ... LAYER OPTICAL DEPTH DUE TO RAYLEIGH SCATTERING
!
        DO K=1,L
          DO I=1,IMAX
            TAURS(I,K)  = TAUR(I,K,IV)
            SSAAER(I,K) = 0.0
            ASYAER(I,K) = 0.0
            TAUAER(I,K) = 0.0
          ENDDO
        ENDDO
!JFM commented out the call to aeros
!        if (kaer .ge. 1) then !==> AEROSOL OPTICAL PROPERTIES
!          CALL AEROS(IV,KPRF,IDXC,CMIX,DENN,RH
!     &,              TAUAER,SSAAER,ASYAER
!     &,              L,IMAX,NAE,NBD,NDM,NXC, NDN
!     &,              HAER,IDM,DZ,HZ)
!        endif
!
!===> ... COMPUTE TOTAL OPTICAL THICKNESS, SINGLE SCATTERING ALBEDO,
!         AND ASYMMETRY FACTOR FOR CLEAR SKY
!
        DO K=1,L
        DO I=1,IMAX
          TAUOZ = AK(IV)*OZ(I,K)
          TAUWV = WK(IV)*WZ(I,K)
          TAUTO(I,K) = MAX(FPMIN, TAUOZ+TAUWV+TAUAER(I,K)+TAURS(I,K))
          SSAT1(I,K) = SSAAER(I,K)*TAUAER(I,K) + TAURS(I,K)
          ASYT1(I,K) = ASYAER(I,K)*SSAAER(I,K)*TAUAER(I,K)
          FFFT1(I,K) = ASYAER(I,K)*ASYT1(I,K) + FFFRS0*TAURS(I,K)
!
          SSATO(I,K) = MIN(FPMAX, SSAT1(I,K)/TAUTO(I,K))
          TEM        = 1.0 / MAX(FPMIN, SSAT1(I,K))
          ASYTO(I,K) = ASYT1(I,K) * TEM
          FFFTO(I,K) = FFFT1(I,K) * TEM
        ENDDO
        ENDDO
!
!===> ... CLEAR SKY FLUXES CALCULATIONS
!
        CALL SWFLUX(TAUTO,SSATO,ASYTO,FFFTO,CSM,ZTH,ALBB,ALBD,
     &              UPFLUX,DWFLUX,DWSFXB,DWSFXD, L, LP1, IMAX)
!    &,             lprnt)
!
        DO K=1,LP1
          DO I=1,IMAX
!JFM and BAO limited FNET0 by .01
            FNET0(I,K) = max(.01,FNET0(I,K) + 
     &                   (DWFLUX(I,K) - UPFLUX(I,K))*SS(IV))
          ENDDO
        ENDDO
        DO I=1,IMAX
          TUPFX0(I) = TUPFX0(I) + UPFLUX(I,1)   * SS(IV)
          SUPFX0(I) = SUPFX0(I) + UPFLUX(I,LP1) * SS(IV)
          SDNFX0(I) = SDNFX0(I) + DWFLUX(I,LP1) * SS(IV)
          DWSFB0(I) = DWSFB0(I) + DWSFXB(I)     * SS(IV)
          DWSFD0(I) = DWSFD0(I) + DWSFXD(I)     * SS(IV)
        ENDDO
 
! --- FOR UV-B FLUX OUTPUT
       IF (IV.GE.5 .AND. IV.LE.6) THEN
        DO I=1,IMAX
          SUVBF0(I) = SUVBF0(I) + DWFLUX(I,LP1) * SS(IV)
        ENDDO
       END IF
! --- END UV-B
 
        IF (NCLOUD .GT. 0) THEN
!
!===> ... COMPUTE TOTAL OPTICAL THICKNESS, SINGLE SCATTERING ALBEDO,
!         AND ASYMMETRY FACTOR FOR CLOUDY SKY
!
        DO K=1,L
          DO I=1,IMAX
            IF (TAUCL(I,K) .GT. 0.0) THEN
              TAUTO(I,K) = TAUCL(I,K) + TAUTO(I,K)
              SSAT1(I,K) = SSACW(I,K) + SSAT1(I,K)
              SSATO(I,K) = MIN(FPMAX, SSAT1(I,K)/TAUTO(I,K))
              TEM        = 1.0  / MAX(FPMIN, SSAT1(I,K))
              ASYTO(I,K) = (ASYCW(I,K) + ASYT1(I,K)) * TEM
              FFFTO(I,K) = (FFFCW(I,K) + FFFT1(I,K)) * TEM
            END IF
          ENDDO
        ENDDO
C
C===> ... CLOUDY SKY FLUXES CALCULATIONS
C
        CALL SWFLUX(TAUTO,SSATO,ASYTO,FFFTO,CSM,ZTH,ALBB,ALBD,
     &              UPFLUX,DWFLUX,DWSFXB,DWSFXD, L, LP1, IMAX)
!    &,             lprnt)
!
        DO K=1,LP1
          DO I=1,IMAX
            FNETC(I,K) = FNETC(I,K) + (DWFLUX(I,K) - UPFLUX(I,K))*SS(IV)
          ENDDO
        ENDDO
        DO I=1,IMAX
          TUPFXC(I) = TUPFXC(I) + UPFLUX(I,1)   * SS(IV)
          SUPFXC(I) = SUPFXC(I) + UPFLUX(I,LP1) * SS(IV)
          SDNFXC(I) = SDNFXC(I) + DWFLUX(I,LP1) * SS(IV)
          DWSFBC(I) = DWSFBC(I) + DWSFXB(I)     * SS(IV)
          DWSFDC(I) = DWSFDC(I) + DWSFXD(I)     * SS(IV)
        ENDDO
 
! --- FOR UV-B FLUX OUTPUT
        IF (IV.GE.5 .AND. IV.LE.6) THEN
         DO I=1,IMAX
          SUVBFC(I) = SUVBFC(I) + DWFLUX(I,LP1) * SS(IV)
         ENDDO
        END IF
! --- END UV-B
 
        ELSE
          DO K=1,LP1
            DO I=1,IMAX
              FNETC(I,K) = FNET0(I,K)
            ENDDO
          ENDDO
          DO I=1,IMAX
            TUPFXC(I) = TUPFX0(I)
            SUPFXC(I) = SUPFX0(I)
            SDNFXC(I) = SDNFX0(I)
            DWSFBC(I) = DWSFB0(I)
            DWSFDC(I) = DWSFD0(I)
          ENDDO

        ENDIF
!
      ENDDO           !    INTEGRATION OVER SPECTRAL BANDS LOOP END
!
      RETURN
      END
 
      SUBROUTINE SOLIR(WH,TA,TAUCL,CSM,ZTH,IBND,FICE,
     &                 KAER,KPRF,IDXC,CMIX,DENN,RH,ALBB,ALBD,
     &                 ICWP,CWP,CIP,CCLY,REW,REI,REWI,REII,
     &                 TUPFXC,SUPFXC,SDNFXC,TUPFX0,SUPFX0,SDNFX0,
     &                 FNET0,FNETC,DWSFB0,DWSFD0,DWSFBC,DWSFDC,
     &                 L,LP1,IMAX,NBD,NVB,NAE,NDM,NXC,NDN,
     &                 HAER,IDM,DZ,HZ,TAUR)
!    &,                lprnt)
!FPP$ NOCONCUR R
!********************************************************************
!  COMPUTE SOLAR FLUX IN THE NIR REGION (3 BANDS, 10-K PER BAND)
!  THE NIR REGION HAS THREE WATER VAPOR BANDS, TEN K's FOR EACH BAND.
!    1.   1000-4400 (/cm)         2.27-10.0 (micron)
!    2.   4400-8200               1.22-2.27
!    3.   8200-14300              0.70-1.22
!
!  INPUT PARAMETERS:                           UNITS
!    WH,TA,TAUCL,CSM,IBND,FICE,KAER,PAER,ALBB,ALBD
!    ICWP,CWP,CIP,CCLV,REW,REI
!  FIXED INPUT DATA:
!    H2O ABSORPTION COEFFICIENT (XK)           CM**2/GM
!    K-DISTRIBUTION FUNCTION    (HK)           FRACTION
!
!  THE FOLLOWING PARAMETERS MUST SPECIFIED BY USERS:
!    CLOUD SINGLE SCATTERING ALBEDO (SACL)     N/D
!    CLOUD ASYMMETRY FACTOR (ASYCL)            N/D
!  AEROSOLS OPTICAL PARAMETERS ARE OBTAINED FROM CALLING
!    SUBPROGRAM AEROS
!
!  OUTPUT PARAMETERS:
!    FNET0  : CLEAR SKY NET FLUX
!    FNETC  : CLOUDY SKY NET FLUX
!    TUPFXC : CLOUDY SKY UPWARD FLUX AT TOA
!    SUPFXC : CLOUDY SKY UPWARD FLUX AT SFC
!    SDNFXC : CLOUDY SKY DOWNWARD FLUX AT SFC
!    TUPFX0 : CLEAR SKY UPWARD FLUX AT TOA
!    SUPFX0 : CLEAR SKY UPWARD FLUX AT SFC
!    SDNFX0 : CLEAR SKY DOWNWARD FLUX AT SFC
!    DWSFB0 : CLEAR SKY SFC DOWN DIR. FLUX
!    DWSFD0 : CLEAR SKY SFC DOWN DIF. FLUX
!    DWSFBC : CLOUDY SKY SFC DOWN DIR. FLUX
!    DWSFDC : CLOUDY SKY SFC DOWN DIF. FLUX
!
!  PROGRAM HISTORY LOG:
!   94-06-12   M.D. CHOU, GLA.
!   95-02-09   YU-TAI HOU      - RECODE FOR NMC MODELS
!   98-08-03   YU-TAI HOU      - UPDATED CLOUD RADIATIVE PROPERTIES
!            CALCULATION. USE SLINGO'S METHOD (JAS 1989) ON WATER
!            CLOUD, EBERT AND CURRY'S METHOD (JGR 1992) ON ICE CLOUD.
!   99-03-25   YU-TAI HOU      - UPDATED CLOUD PROPERTIES USE THE
!            MOST RECENT CHOU ET AL. DATA (J. CLIM 1998)
!   99-04-27   YU-TAI HOU      - UPDATED CLOUD RADIATIVE PROPERTY
!            CALCULATIONS USE LINEAR T-ADJUSTED METHOD.
!C   99-04-27   YU-TAI HOU      - UPDATED CLOUD RADIATIVE PROPERTY
!            CALCULATIONS USE LINEAR T-ADJUSTED METHOD.
!   99-09-13   YU-TAI HOU      - UPDATED TO CHOU'S JUNE,1999 VERSION
!
!********************************************************************
!
!
      USE MACHINE , ONLY : kind_rad
      implicit none
!
      integer NRB, NK0
      PARAMETER (NRB=4,NK0=10)
!
      integer L, LP1, IMAX, NBD, NVB, NAE, NDM, NXC, NDN
     &,       KAER, ICWP, IBND
!
      integer     IDM (IMAX,L,NAE), IDXC(NXC,IMAX), KPRF(IMAX)
      real (kind=kind_rad) HAER(NDM,NAE)
     &,                     DZ(IMAX,L), HZ(IMAX,L+1), TAUR(IMAX,L,NBD)
     &,                     CMIX(NXC,IMAX), DENN(NDN,IMAX)
!
! --- INPUT
      real (kind=kind_rad)
     &  WH(IMAX,L),   TAUCL(IMAX,L), CSM(IMAX,L), RH(IMAX,L)
     &, ALBB(IMAX),   ALBD(IMAX),    ZTH(IMAX,L)
     &, CWP(IMAX,L),  CIP(IMAX,L),   REW(IMAX,L), REI(IMAX,L)
     &, CCLY(IMAX,L), TA(IMAX,L),    FICE(IMAX,L)
     &, REWI(IMAX,L), REII(IMAX,L)
C --- OUTPUT
      real (kind=kind_rad)
     &  FNET0 (IMAX,LP1), DWSFB0(IMAX), DWSFD0(IMAX)
     &, FNETC (IMAX,LP1), DWSFBC(IMAX), DWSFDC(IMAX)
     &, TUPFXC(IMAX),     SUPFXC(IMAX), SDNFXC(IMAX)
     &, TUPFX0(IMAX),     SUPFX0(IMAX), SDNFX0(IMAX)
!
      Integer ncloud
      logical cloudy(imax)
!     logical lprnt
!
! --- TEMPORARY ARRAY
      real (kind=kind_rad)
     &  UPFLUX(IMAX,LP1), DWFLUX(IMAX,LP1)
     &, DWSFXB(IMAX),     DWSFXD(IMAX)
     &, TAUTO (IMAX,L),   SSATO (IMAX,L),     ASYTO (IMAX,L)
     &, TAURS (IMAX,L),   SSAT1 (IMAX,L),     ASYT1 (IMAX,L)
     &, TAUAER(IMAX,L),   SSAAER(IMAX,L),     ASYAER(IMAX,L)
     &, XK  (NK0),        HK  (NK0,NRB)
!0499 --- T ADJUSTED CLD PROPERTY METHOD
     &, SSAW0(NRB,2),   SSAI0(NRB,2),   ASYW0(NRB,2), ASYI0(NRB,2)
     &, FFFCW (IMAX,L), FFFT1 (IMAX,L), FFFTO (IMAX,L)
     &, ASYCW(IMAX,L),  SSACW(IMAX,L)
!
      real (kind=kind_rad)
     &  A0W(NRB,2), A1W(NRB,2), B0W(NRB,2), B1W(NRB,2)
     &, A0I(NRB,2), A1I(NRB,2), C0W(NRB,2), C1W(NRB,2)
     &, B0I(NRB,2), B1I(NRB,2), B2I(NRB,2), FACW(IMAX,L)
     &, C0I(NRB,2), C1I(NRB,2), C2I(NRB,2), FACI(IMAX,L)
     &, FFFRS0,     FPMIN,      FPMAX
!
      DATA XK / 0.0010, 0.0133, 0.0422, 0.1334, 0.4217,
     &          1.3340, 5.6230, 31.620, 177.80, 1000.0 /
      DATA HK / .01074, .00360, .00411, .00421, .00389,
     &          .00326, .00499, .00465, .00245, .00145,
     &          .08236, .01157, .01133, .01143, .01240,
     &          .01258, .01381, .00650, .00244, .00094,
     &          .20673, .03497, .03011, .02260, .01336,
     &          .00696, .00441, .00115, .00026, .00000,
     &          .29983, .05014, .04555, .03824, .02965,
     &          .02280, .02321, .01230, .00515, .00239 /
!
      DATA SSAW0/.7578,.9869,.9997,.9869, .7570,.9868,.9998,.9916/
     &,    ASYW0/.8678,.8185,.8354,.8315, .8723,.8182,.8354,.8311/
     &,    SSAI0/.7283,.9442,.9994,.9620, .7368,.9485,.9995,.9750/
     &,    ASYI0/.9058,.8322,.8068,.8220, .9070,.8304,.8067,.8174/
      DATA FFFRS0 / 0.1 /
      DATA FPMIN,FPMAX /1.0E-8, 0.999999/
!
                         D A T A
!0499 - T-ADJUSTED CLD PROP COEFF, WATER CLOUD
     &    A0W / 1.466E-2, 2.276E-2, 2.654E-2, 2.494E-2
     &,         1.528E-2, 2.286E-2, 2.642E-2, 2.517E-2 /
     &,   A1W / 1.617E+0, 1.451E+0, 1.351E+0, 1.392E+0
     &,         1.611E+0, 1.449E+0, 1.353E+0, 1.386E+0 /
     &,   B0W / 1.708E-1, 5.314E-4,-4.594E-6, 6.473E-3
     &,         1.674E-1, 5.427E-4,-3.306E-6, 3.218E-3 /
     &,   B1W / 7.142E-3, 1.258E-3, 2.588E-5, 6.649E-4
     &,         7.561E-3, 1.263E-3, 2.287E-5, 5.217E-4 /
     &,   C0W / 8.266E-1, 7.507E-1, 7.925E-1, 7.811E-1
     &,         8.344E-1, 7.501E-1, 7.922E-1, 7.808E-1 /
     &,   C1W / 4.119E-3, 6.770E-3, 4.297E-3, 5.034E-3
     &,         3.797E-3, 6.812E-3, 4.323E-3, 5.031E-3 /
!
                         D A T A
!0499 - T-ADJUSTED CLD PROP COEFF, ICE CLOUD
     &    A0I / 2.822E-4,-3.248E-5,-3.758E-5,-1.214E-5
     &,         2.712E-4,-4.308E-5,-3.917E-5,-2.456E-5 /
     &,   A1I / 2.491E+0, 2.522E+0, 2.522E+0, 2.520E00
     &,         2.489E+0, 2.523E+0, 2.522E+0, 2.521E00 /
     &,   B0I / 1.853E-1, 2.544E-3,-7.701E-7, 1.461E-2
     &,         1.738E-1, 2.461E-3,-8.979E-7, 7.083E-3 /
     &,   B1I / 1.841E-3, 1.023E-3, 9.849E-6, 4.612E-4
     &,         1.887E-3, 9.436E-4, 8.102E-6, 3.495E-4 /
     &,   B2I /-6.671E-6,-2.266E-6,-.3988E-9,-1.202E-6
     &,        -6.615E-6,-2.107E-6,-.1862E-9,-8.500E-7 /
     &,   C0I / 8.388E-1, 7.572E-1, 7.519E-1, 7.600E-1
     &,         8.414E-1, 7.566E-1, 7.519E-1, 7.566E-1 /
     &,   C1I / 1.519E-3, 1.563E-3, 1.099E-3, 1.275E-3
     &,         1.477E-3, 1.537E-3, 1.097E-3, 1.241E-3 /
     &,   C2I /-6.702E-6,-5.232E-6,-3.081E-6,-4.020E-6
     &,        -6.403E-6,-5.130E-6,-3.070E-6,-3.804E-6 /
!
c$$$      SAVE XK, HK, SSAW0,SSAI0, ASYW0,ASYI0, FFFRS0, FPMIN, FPMAX
c$$$      SAVE A0W,A1W,B0W,B1W,C0W,C1W,A0I,A1I,B0I,B1I,B2I,C0I,C1I,C2I
!
      real (kind=kind_rad) tau1,  tau2,  ssa1,  ssa2,  asy1,  asy2
     &,                     ssaw1, ssaw2, asyw1, asyw2, tauwv
     &,                     tem
      integer i, k, ibb1, ibb2, ib, ib1, ik
cc
cc--------------------------------------------------------------------
cc
      real(kind=kind_rad) cons_0               !constant
      real(kind=kind_rad) cons_10              !constant
      real(kind=kind_rad) cons_30              !constant
!     integer ipnGlobal,its,mype
!     logical DiagPrint
!     call PhysicsGetIpnItsMype(ipnGlobal,its,mype,DiagPrint)

cc
      cons_0          =         0.d0           !constant
      cons_10         =        10.d0           !constant
      cons_30         =        30.d0           !constant
cc
cc--------------------------------------------------------------------
cc
!
      DO K=1,L
        DO I=1,IMAX
          FACW(I,K) = MAX(cons_0, MIN(cons_10,273.15-TA(I,K)))*0.1      !constant
          FACI(I,K) = MAX(cons_0, MIN(cons_30,263.15-TA(I,K)))/30.0     !constant
        ENDDO
      ENDDO
!
!===> ... LOOP OVER THREE NIR BANDS
!
      IF (IBND .EQ. 1) THEN
        IBB1 = NRB
        IBB2 = NRB
      ELSE
        IBB1 = 1
        IBB2 = NRB - 1
      END IF
      DO IB=IBB1,IBB2
        IB1 = NVB + IB
!
!===> ... LAYER OPTICAL DEPTH DUE TO RAYLEIGH SCATTERING
!
        DO K=1,L
          DO I=1,IMAX
            TAURS(I,K)  = TAUR(I,K,IB1)
            SSAAER(I,K) = 0.0
            ASYAER(I,K) = 0.0
            TAUAER(I,K) = 0.0
          ENDDO
        ENDDO
!JFM commented out the call to aeros
!        if (kaer .ge. 1) then !==> AEROSOL OPTICAL PROPERTIES
!          CALL AEROS(IB1,KPRF,IDXC,CMIX,DENN,RH
!     &,              TAUAER,SSAAER,ASYAER
!     &,              L,IMAX,NAE,NBD,NDM,NXC, NDN
!     &,              HAER,IDM,DZ,HZ)
!        endif
        cloudy(:) = .false.
!
!0898 ... GET CLOUD PROPERTIES FROM CWP AND CIP
!
      IF (ICWP .EQ. 1) THEN
        DO K=1,L
        DO I=1,IMAX
          IF (CCLY(I,K) .GT. 0.0) THEN
! --- T-ADJ METHOD
            TAU1=CWP(I,K)*(  FACW(I,K) *(A0W(IB,1)+A1W(IB,1)*REWI(I,K))
     &                 +(1.0-FACW(I,K))*(A0W(IB,2)+A1W(IB,2)*REWI(I,K)))
            TAU2=CIP(I,K)*(  FACI(I,K) *(A0I(IB,1)+A1I(IB,1)*REII(I,K))
     &                 +(1.0-FACI(I,K))*(A0I(IB,2)+A1I(IB,2)*REII(I,K)))
            TAUCL(I,K) = TAU1 + TAU2
            SSA1 = 1.0 - (  FACW(I,K) *(B0W(IB,1)+B1W(IB,1)*REW(I,K))
     &               + (1.0-FACW(I,K))*(B0W(IB,2)+B1W(IB,2)*REW(I,K)))
            SSA2 = 1.0 - (  FACI(I,K) *(B0I(IB,1)
     &                 + (B1I(IB,1)+B2I(IB,1)*REI(I,K))*REI(I,K))
     &               + (1.0-FACI(I,K))*(B0I(IB,2)
     &                 + (B1I(IB,2)+B2I(IB,2)*REI(I,K))*REI(I,K)) )
            SSAW1 = SSA1 * TAU1
            SSAW2 = SSA2 * TAU2
            SSACW(I,K) = SSAW1 + SSAW2
            ASY1 =    FACW(I,K) *(C0W(IB,1)+C1W(IB,1)*REW(I,K))
     &         + (1.0-FACW(I,K))*(C0W(IB,2)+C1W(IB,2)*REW(I,K))
            ASY2 =    FACI(I,K) *(C0I(IB,1)
     &              + (C1I(IB,1)+C2I(IB,1)*REI(I,K))*REI(I,K))
     &         + (1.0-FACI(I,K))*(C0I(IB,2)
     &              + (C1I(IB,2)+C2I(IB,2)*REI(I,K))*REI(I,K))
            ASYW1 = ASY1 * SSAW1
            ASYW2 = ASY2 * SSAW2
            ASYCW(I,K) = ASYW1 + ASYW2
            FFFCW(I,K) = ASY1*ASYW1 + ASY2*ASYW2
            cloudy(i)  = .true.
          ELSE
            TAUCL(I,K) = 0.0
            SSACW(I,K) = 1.0
            ASYCW(I,K) = 0.0
            FFFCW(I,K) = 0.0
          END IF
        ENDDO
        ENDDO
      ELSE
        DO K=1,L
        DO I=1,IMAX
          IF (TAUCL(I,K) .GT. 0.0) THEN
            TAU2 = FICE(I,K) * TAUCL(I,K)
            TAU1 = TAUCL(I,K) - TAU2
            SSA1 = FACW(I,K)*SSAW0(IB,1) + (1.0-FACW(I,K))*SSAW0(IB,2)
            SSA2 = FACI(I,K)*SSAI0(IB,1) + (1.0-FACI(I,K))*SSAI0(IB,2)
            SSAW1 = SSA1 * TAU1
            SSAW2 = SSA2 * TAU2
!           SSA1 = (1.0-FICE(I,K)) * (FACW(I,K) * SSAW0(IB,1)
!    &                         + (1.0-FACW(I,K))* SSAW0(IB,2))
!           SSA2 =      FICE(I,K)  * (FACI(I,K) * SSAI0(IB,1)
!    &                         + (1.0-FACI(I,K))* SSAI0(IB,2))
!           SSAW1 = SSA1 * TAUCL(I,K)
!           SSAW2 = SSA2 * TAUCL(I,K)
            SSACW(I,K) = SSAW1 + SSAW2
!           ASY1 = (1.0-FICE(I,K)) * (FACW(I,K) * ASYW0(IB,1)
!    &                         + (1.0-FACW(I,K))* ASYW0(IB,2))
!           ASY2 =      FICE(I,K)  * (FACI(I,K) * ASYI0(IB,1)
!    &                         + (1.0-FACI(I,K))* ASYI0(IB,2))
            ASY1 = FACW(I,K)*ASYW0(IB,1) + (1.0-FACW(I,K))*ASYW0(IB,2)
            ASY2 = FACI(I,K)*ASYI0(IB,1) + (1.0-FACI(I,K))*ASYI0(IB,2)
            ASYW1 = ASY1 * SSAW1
            ASYW2 = ASY2 * SSAW2
            ASYCW(I,K) = ASYW1 + ASYW2
            FFFCW(I,K) = ASY1*ASYW1 + ASY2*ASYW2
            cloudy(i) = .true.
          ELSE
            SSACW(I,K) = 1.0
            ASYCW(I,K) = 0.0
            FFFCW(I,K) = 0.0
          END IF
        ENDDO
        ENDDO
      END IF
!
      NCLOUD = 0
      DO I=1,IMAX
        if (cloudy(i)) ncloud = ncloud + 1
      ENDDO
!
!===> ... IK IS THE INDEX FOR THE K-DISTRIBUTION FUNCTION (OR THE
!     ABSORPTION COEFFICIENT)
!
        DO IK=1,NK0
!
         IF (HK(IK,IB) .GE. 0.00001) THEN
!
!===> ... COMPUTE TATAL OPTICAL THICKNESS, SINGLE SCATTERING ALBEDO,
!         AND ASYMMETRY FACTOR FOR CLEAR SKY
!
           DO K=1,L
           DO I=1,IMAX
             TAUWV      = XK(IK)*WH(I,K)
             TAUTO(I,K) = MAX(FPMIN, TAUWV+TAUAER(I,K)+TAURS(I,K))
             SSAT1(I,K) = SSAAER(I,K)*TAUAER(I,K)+TAURS(I,K)
             ASYT1(I,K) = ASYAER(I,K)*SSAAER(I,K)*TAUAER(I,K)
             FFFT1(I,K) = ASYAER(I,K)*ASYT1(I,K) + FFFRS0*TAURS(I,K)
             SSATO(I,K) = MIN(FPMAX, SSAT1(I,K)/TAUTO(I,K))
             TEM        = 1.0 / MAX(FPMIN, SSAT1(I,K))
             ASYTO(I,K) = ASYT1(I,K) * TEM
             FFFTO(I,K) = FFFT1(I,K) * TEM
           ENDDO
           ENDDO
!
!===> ... CLEAR SKY FLUXES CALCULATIONS
!
           CALL SWFLUX(TAUTO,SSATO,ASYTO,FFFTO,CSM,ZTH,ALBB,ALBD,
     &                 UPFLUX,DWFLUX,DWSFXB,DWSFXD, L, LP1, IMAX)
!    &,                lprnt)
!
           DO K=1,LP1
             DO I=1,IMAX
!JFM and BAO limited FNET0 to .01
               FNET0 (I,K) = max(.01,FNET0 (I,K)
     &                     + (DWFLUX(I,K) - UPFLUX(I,K))*HK(IK,IB))
             ENDDO
           ENDDO
           DO I=1,IMAX
             TUPFX0(I) = TUPFX0(I) + UPFLUX(I,1)   * HK(IK,IB)
             SUPFX0(I) = SUPFX0(I) + UPFLUX(I,LP1) * HK(IK,IB)
             SDNFX0(I) = SDNFX0(I) + DWFLUX(I,LP1) * HK(IK,IB)
             DWSFB0(I) = DWSFB0(I) + DWSFXB(I)     * HK(IK,IB)
             DWSFD0(I) = DWSFD0(I) + DWSFXD(I)     * HK(IK,IB)
           ENDDO
           IF (NCLOUD .GT. 0) THEN
!
!===> ... COMPUTE TATAL OPTICAL THICKNESS, SINGLE SCATTERING ALBEDO,
!         AND ASYMMETRY FACTOR FOR CLOUDY SKY
!
           DO K=1,L
             DO I=1,IMAX
               IF (TAUCL(I,K) .GE. 0.001) THEN
                 TAUTO(I,K) = TAUCL(I,K) + TAUTO(I,K)
                 SSAT1(I,K) = SSACW(I,K) + SSAT1(I,K)
                 SSATO(I,K) = MIN(FPMAX, SSAT1(I,K)/TAUTO(I,K))
                 TEM        = 1.0 / MAX(FPMIN, SSAT1(I,K))
                 ASYTO(I,K) = (ASYCW(I,K) + ASYT1(I,K)) * TEM
                 FFFTO(I,K) = (FFFCW(I,K) + FFFT1(I,K)) * TEM
               END IF
             ENDDO
           ENDDO
!
!===> ... CLOUDY SKY FLUXES CALCULATIONS
!
           CALL SWFLUX(TAUTO,SSATO,ASYTO,FFFTO,CSM,ZTH,ALBB,ALBD,
     &                 UPFLUX,DWFLUX,DWSFXB,DWSFXD, L, LP1, IMAX)
!    &,                lprnt)
!
           DO K=1,LP1
             DO I=1,IMAX
               FNETC(I,K) = FNETC(I,K)
     &                    + (DWFLUX(I,K) - UPFLUX(I,K))*HK(IK,IB)
             ENDDO
           ENDDO
           DO I=1,IMAX
             TUPFXC(I) = TUPFXC(I) + UPFLUX(I,1)   * HK(IK,IB)
             SUPFXC(I) = SUPFXC(I) + UPFLUX(I,LP1) * HK(IK,IB)
             SDNFXC(I) = SDNFXC(I) + DWFLUX(I,LP1) * HK(IK,IB)
             DWSFBC(I) = DWSFBC(I) + DWSFXB(I)     * HK(IK,IB)
             DWSFDC(I) = DWSFDC(I) + DWSFXD(I)     * HK(IK,IB)
           ENDDO
           ELSE
             DO K=1,LP1
               DO I=1,IMAX
                 FNETC(I,K) = FNET0(I,K)
               ENDDO
             ENDDO
             DO I=1,IMAX
               TUPFXC(I) = TUPFX0(I)
               SUPFXC(I) = SUPFX0(I)
               SDNFXC(I) = SDNFX0(I)
               DWSFBC(I) = DWSFB0(I)
               DWSFDC(I) = DWSFD0(I)
             ENDDO
           ENDIF
!
         ENDIF
        ENDDO           ! K-distribution loop ends here
      ENDDO           ! Loop over NIR bands ends here
!
      RETURN
      END
 
      SUBROUTINE SWFLUX(TAU,SSC,G0,FF,CSM,ZTH,ALB,ALD,
     &                  UPFLUX,DWFLUX,DWSFCB,DWSFCD, L, LP1, IMAX)
!    &,                lprnt)
!FPP$ NOCONCUR R
!********************************************************************
!  USES THE DELTA-EDDINGTON APPROXIMATION TO COMPUTE THE BULK
!  SCATTERING PROPERTIES OF A SINGLE LAYER CODED FOLLOWING
!  COAKLEY ET AL.  (JAS, 1982)
!
!  INPUTS:
!    TAU: THE EFFECTIVE OPTICAL THICKNESS
!    SSC: THE EFFECTIVE SINGLE SCATTERING ALBEDO
!    G0:  THE EFFECTIVE ASYMMETRY FACTOR
!    FF:  THE EFFECTIVE FORWARD SCATTERING FACTOR
!    CSM: THE EFFECTIVE SECANT OF THE ZENITH ANGLE
!    ALB: SURFACE ALBEDO FOR DIRECT RADIATION
!    ALD: SURFACE ALBEDO FOR DIFFUSED RADIATION
!
!  OUTPUTS:
!    UPFLUX: UPWARD FLUXES
!    DWFLUX: DOWNWARD FLUXES
!    DWSFCB: DOWNWARD SURFACE FLUX DIRECT COMPONENT
!    DWSFCD: DOWNWARD SURFACE FLUX DIFFUSED COMPONENT
!********************************************************************
!
!
      USE MACHINE , ONLY : kind_rad
      implicit none
!
      integer L, LP1, IMAX
!
! --- INPUT
      real (kind=kind_rad)
     &  TAU(IMAX,L), SSC(IMAX,L), G0(IMAX,L), FF(IMAX,L)
     &, CSM(IMAX,L), ZTH(IMAX,L), ALB(IMAX),  ALD(IMAX)
! --- OUTPUT
      real (kind=kind_rad)
     &  UPFLUX(IMAX,LP1),DWFLUX(IMAX,LP1),DWSFCB(IMAX),DWSFCD(IMAX)
! --- TEMPORARY
      real (kind=kind_rad)
     &  TTB(IMAX,LP1),TDN(IMAX,LP1),RUP(IMAX,LP1), TT (IMAX,LP1,2)
     &, RFU(IMAX,LP1),RFD(IMAX,LP1),TB (IMAX,LP1), RR (IMAX,LP1,2)
!     logical lprnt
!
      real (kind=kind_rad) ZTHD, CSMD,  EPSLN, AA,   TAUP, SSCP, GP
     &,                     OMS1, OGS1,  TLAM,  U1,   U1P1, U1M1, E1
     &,                     U1E,  U1EPE, U1EME, DEN,  RF1,  TF1,  ZZTH
     &,                     ZZ,   DEN1,  GAMA,  ALFA, AMG,  APG,  ZA
     &,                     SLAM, BB
      integer I, K, KM1
cc
cc--------------------------------------------------------------------
cc
      real(kind=kind_rad) cons_0              !constant
      real(kind=kind_rad) cons_m30            !constant
      real(kind=kind_rad) cons_30             !constant
!     integer ipnGlobal,its,mype
!     logical DiagPrint
!     call PhysicsGetIpnItsMype(ipnGlobal,its,mype,DiagPrint)

cc
      cons_0         =        0.d0            !constant
      cons_m30       =      -30.d0            !constant
      cons_30        =       30.d0            !constant
cc
cc--------------------------------------------------------------------
cc
!
!===> ... DIFFUSE INCIDENT RADIATION IS APPROXIMATED BY BEAM RADIATION
!         WITH AN INCIDENT ANGLE OF 53 DEGREES. COS(53) = 0.602
      ZTHD = 0.602
      CSMD = 1.0 / ZTHD
      epsln = 1.0e-30
!
!===> ... DELTA-EDDINGTON SCALING OF SINGLE SCATTERING ALBEDO,
!         OPTICAL THICKNESS, AND ASYMMETRY FACTOR, K & H EQS(27-29)
!

      DO K=1,L
        DO I=1,IMAX
          AA    = 1.0 - FF(I,K)*SSC(I,K)
          TAUP  = TAU(I,K) * AA
          SSCP  = SSC(I,K) * (1.0 - FF(I,K)) / AA
          GP    = (G0(I,K) - FF(I,K)) / (1.0 - FF(I,K))
!
          OMS1  = 1.0 - SSCP
          OGS1  = 1.0 - SSCP*GP
          TLAM  = 3.0 * OMS1*OGS1
          SLAM  = SQRT(TLAM)
!
          U1    = 1.5 * OGS1 / SLAM
          U1P1  = U1 + 1.0
          U1M1  = U1 - 1.0
          E1    = EXP(max(-TAUP*SLAM, cons_m30))     !constant
          U1E   = U1 * E1
          U1EPE = U1E + E1
          U1EME = U1E - E1
          DEN   = 1.0 / ((U1P1 + U1EME)*(U1P1 - U1EME))
          RF1   = (U1P1 + U1EPE) * (U1M1 - U1EME) * DEN
          TF1   = 4.0 * U1E * DEN
!
!===> ... COMPUTE LAYER TRANSMISSIONS AND REFLECTIONS
!         (I,K,J) J=1,2 FOR LAYER K ILLUMINATED BY DIFFUSE AND
!                       DIRECT INCOMING RADIATION
!         RR   :  LAYER REFLECTION
!         TT   :  LAYER TOTAL TRANSMISSION
!         TB   :  LAYER DIRECT TRANSMISSION
!
!      Diffuse Radiation
!      -----------------
          ZZTH = ZTHD
          ZZ   = ZZTH * ZZTH
          DEN1 = 1.0 - TLAM*ZZ
          IF (ABS(DEN1) .LT. 1.0E-8) THEN    !===> ... SAFETY CHECK
            ZZTH = ZZTH + 0.001
            ZZ   = ZZTH * ZZTH
            DEN1 = 1.0 - TLAM*ZZ
          END IF
          DEN1   = SSCP / DEN1
!
          GAMA   = 0.50 * (1.0 + 3.0*GP*OMS1*ZZ)*DEN1
          ALFA   = 0.75 * ZTHD * (GP + OGS1)*DEN1
          AMG    = ALFA - GAMA
          APG    = ALFA + GAMA
!
          TB(I,K)   = EXP( -MIN(cons_30, TAUP*CSMD) )     !constant
          ZA        = AMG * TB(I,K)
          RR(I,K,1) = ZA*TF1 + APG*RF1 - AMG
          TT(I,K,1) = ZA*RF1 + APG*TF1 + (1.0-APG)*TB(I,K)
!
!      Direct Radiation
!      ----------------
          ZZTH = ZTH(I,K)
          ZZ   = ZZTH * ZZTH
          DEN1 = 1.0 - TLAM*ZZ
          IF (ABS(DEN1) .LT. 1.0E-8) THEN    !===> ... SAFETY CHECK
            ZZTH = ZZTH + 0.001
            ZZ   = ZZTH * ZZTH
            DEN1 = 1.0 - TLAM*ZZ
          END IF
          DEN1   = SSCP / DEN1
!
          GAMA   = 0.50 * (1.0 + 3.0*GP*OMS1*ZZ)*DEN1
          ALFA   = 0.75 * ZTH(I,K) * (GP + OGS1)*DEN1
          AMG    = ALFA - GAMA
          APG    = ALFA + GAMA
!
          TB(I,K)   = EXP( -MIN(cons_30, TAUP*CSM(I,K)) )     !constant
          ZA        = AMG * TB(I,K)
          RR(I,K,2) = ZA*TF1 + APG*RF1 - AMG
          TT(I,K,2) = ZA*RF1 + APG*TF1 + (1.0-APG)*TB(I,K)
!
          TB(I,K)   = MAX(cons_0, TB(I,K))       !constant
          RR(I,K,2) = MAX(cons_0, RR(I,K,2))     !constant
          TT(I,K,2) = MAX(cons_0, TT(I,K,2))     !constant
          RR(I,K,1) = MAX(cons_0, RR(I,K,1))     !constant
          TT(I,K,1) = MAX(cons_0, TT(I,K,1))     !constant
        ENDDO
      ENDDO
!
! --- AT THE SURFACE
!
      DO I=1,IMAX
        TB(I,LP1)   = 0.0
        RR(I,LP1,2) = ALB(I)
        TT(I,LP1,2) = 0.0
        RR(I,LP1,1) = ALD(I)
        TT(I,LP1,1) = 0.0
      END DO
!
      DO I=1,IMAX
        TTB(I,1) = TB(I,1)
        TDN(I,1) = TT(I,1,2)
        RFD(I,1) = RR(I,1,1)
        TTB(I,L) = 0.0
      ENDDO
!
!===> ... LAYERS ADDED DOWNWARD STARTING FROM TOP
!
      DO K=2,LP1
        DO I=1,IMAX
           DEN = TT(I,K,1) / (1.0 - RFD(I,K-1) * RR(I,K,1))
           TTB(I,K) = TTB(I,K-1) * TB(I,K)
           if (ttb(i,k) .lt. epsln) ttb(i,k) = 0.0
           TDN(I,K) = TTB(I,K-1)*TT(I,K,2)+(TDN(I,K-1)-TTB(I,K-1)
     1              + TTB(I,K-1)*RR(I,K,2)*RFD(I,K-1)) * DEN
           RFD(I,K) = RR(I,K,1) + TT(I,K,1)*RFD(I,K-1) * DEN
        ENDDO
      ENDDO
!
!===> ... LAYERS ADDED UPWARD STARTING FROM SURFACE
!
      DO I=1,IMAX
        RFU(I,LP1) = RR(I,LP1,1)
        RUP(I,LP1) = RR(I,LP1,2)
      ENDDO
      DO K=L,1,-1
        DO I=1,IMAX
           DEN = TT(I,K,1) / (1.0 - RFU(I,K+1) * RR(I,K,1))
           RUP(I,K) = RR(I,K,2) + ((TT(I,K,2)-TB(I,K))*RFU(I,K+1)
     1              + TB(I,K)*RUP(I,K+1)) * DEN
           RFU(I,K) = RR(I,K,1) + TT(I,K,1)*RFU(I,K+1) * DEN
        ENDDO
      ENDDO
!
!===> ... FIND UPWARD AND DOWNWARD FLUXES
!
      DO I=1,IMAX
          UPFLUX(I,1) = RUP(I,1)
          DWFLUX(I,1) = 1.0
      ENDDO
      DO K=2,LP1
        KM1 = K - 1
        DO I=1,IMAX
            DEN         = 1.0 / (1.0 - RFD(I,KM1)*RFU(I,K))
            AA          = TTB(I,KM1) * RUP(I,K)
            BB          = TDN(I,KM1) - TTB(I,KM1)
            UPFLUX(I,K) = (AA + BB*RFU(I,K)) * DEN
            DWFLUX(I,K) = TTB(I,KM1) + (AA*RFD(I,KM1) + BB) * DEN
        ENDDO
      ENDDO

!
!===> ... SURFACE DOWNWARD FLUXES
!
      DO I=1,IMAX
        DWSFCB(I) = TTB(I,L)
        DWSFCD(I) = DWFLUX(I,LP1)-DWSFCB(I)
      ENDDO
!
      RETURN
      END
 
      SUBROUTINE FLXCO2(SWC,U1,DU,NU,SWH,W1,DW,NW,TBL,DFLX
     &,                                               L, IMAX)
!FPP$ NOCONCUR R
!********************************************************************
!  COMPUTE THE ABSORPTION DUE TO CO2. REF: CHOU (J. CLIMATE, 1990,
!     209-217)
!     UPDATED SEP. 1999 BASED ON NASA/TM-1999-104606, VOL 15.
!  THE EFFECT OF CO2 ABSORPTION BELOW THE CLOUD TOP IS NEGLECTED.
!  INPUT VARIABLES:
!     SWC         : COLUMN AMOUNT OF CO2
!     SWH         : COLUMN AMOUNT OF WATER VAPOR
!     U1,DU,W1,DW : COEFFICIENTS
!     TBL         : LOOK UP CO2 ABSORPTION TABLE
!     NU,NW       : TABLE DIMENSIONS
!  OUTPUT VARIABLES:
!     DFLX        : ADDITIONAL FLUX REDUCTION DUE TO CO2 FOR CLEAR SKY
!
!********************************************************************
!
      USE MACHINE , ONLY : kind_rad
      implicit none
!
      integer imax, L, NU, NW
      real (kind=kind_rad) SWC(IMAX,L+1),   SWH(IMAX,L+1)
     &,                     DFLX(IMAX,L+1),  TBL(NU,NW)
      real (kind=kind_rad) u1, du, w1, dw, x1, y1, clog, wlog, dc, dd
     &,                     x2, y2
      integer i, k, ic, iw, ic1, iw1
!
! ... TABLE LOOK-UP FOR THE REDUCTION OF CLEAR-SKY SOLAR
!
      X1 = U1 - 0.5*DU
      Y1 = W1 - 0.5*DW
      DO K=2,L+1
        DO I=1,IMAX
            CLOG = SWC(I,K)
            WLOG = SWH(I,K)
            IC = INT( (CLOG - X1)/DU + 1.0)
            IW = INT( (WLOG - Y1)/Dw + 1.0)
            IC = MAX(2, MIN(NU, IC))
            IW = MAX(2, MIN(NW, IW))
            IC1 = IC - 1
            IW1 = IW - 1
            DC = CLOG - FLOAT(IC-2)*DU - U1
            DD = WLOG - FLOAT(IW-2)*DW - W1
            X2 = TBL(IC1,IW1) + (TBL(IC1,IW)-TBL(IC1,IW1))/DW * DD
            Y2 = X2 + (TBL(IC,IW1) - TBL(IC1,IW1))/DU * DC
            DFLX(I,K) = DFLX(I,K) + Y2
        ENDDO
      ENDDO
!
      RETURN
      END
      SUBROUTINE AEROS(IB,KPRF,IDXC,CMIX,DENN,RH
     &,                TAU,SSA,ASY
     &,                L, IMAX, NAE, NBD, NDM, NXC, NDN
     &,                HAER, IDM, DZ, HZ)
!FPP$ NOCONCUR R
!********************************************************************
!  COMPUTE AEROSOLS OPTICAL PROPERTIES IN EIGHT UV+VIS BANDS AND
!  FOUR NIR BANDS. THERE ARE SEVEN DIFFERENT VERTICAL PROFILE
!  STRUCTURES. IN THE TROPOSPHERE, AEROSOL DISTRIBUTION AT EACH
!  GRID POINT IS COMPOSED FROM UP TO SIX COMPONENTS OUT OF A TOTAL
!  OF TEN DIFFERENT SUBSTANCES.
!  REF: WMO REPORT WCP-112 (1986)
!
!  1999-10-13 Y.H. UPDATED TO OPAC DATA (1998)
!
!  BAND: 1. 0.175-0.225 (UV-C)     2. 0.225-0.245;0.260-0.280 (UV-C)
!        3. 0.245-0.260 (UV-C)     4. 0.280-0.295 (UV-B)
!        5. 0.295-0.310 (UV-B)     6. 0.310-0.320 (UV-B)
!        7. 0.320-0.400 (UV-A)     8. 0.400-0.700 (PAR)
!        9. 2.27 - 4.0  (NIR)     10. 1.22 - 2.27 (NIR)
!       11. 0.70 - 1.22 (NIR)     12. 0.70 - 4.0  (NIR)
!
!  INPUT VARIABLES:
!     IB   - SPECTRAL BAND INDEX                   -     1
!     KAER - =0 DO NOT COMPUTE AEROSOLS            -     1
!            =1 COMPUTE AEROSOL PROFILES
!     KPRF - INDECIES OF AEROSOL PROF STRUCTURES   -     IMAX
!     IDXC - INDECIES OF AEROSOL COMPONENTS        -     NXC*IMAX
!     CMIX - MIXING RATIOES OF AEROSOL COMPONENTS  -     NXC*IMAX
!     DENN - AEROSOL NUMBER DENSITIES              -     NDN*IMAX
!     RH   - RELATIVE HUMIDITY                     -     IMAX*L
!
!  OUTPUT VARIABLES:
!     TAU  - OPTICAL DEPTH                         -     IMAX*L
!     SSA  - SINGLE SCATTERING ALBEDO              -     IMAX*L
!     ASY  - ASYMMETRY PARAMETER                   -     IMAX*L
!     TAURS- RAYLEIGH SCATTERING OPTICAL DEPTH     -     IMAX*L
!
!  VARIALBES IN COMMON BLOCK:
!     HAER - SCALE HEIGHT OF AEROSOLS              KM    NDM*NAE
!     IDM  - AEROSOL DOMAIN INDEX                  -     L*NAE
!     DZ   - LAYER THICKNESS                       KM    L
!     HZ   - LEVEL HEIGHT                          KM    L+1
!     TAUR - RAYLEIGH SCATTERING OPTICAL DEPTH     -     L*NBD
!********************************************************************
!
      USE MACHINE , ONLY : kind_rad
      implicit none
!
      integer ncm1, ncm2, ncm, ncf, nbdd
      PARAMETER (NCM1=6, NCM2=4, NCM=NCM1+NCM2, NCF=3, NBDD=12)
!
! --- INPUT
      integer L, IMAX, NAE, NBD, NDM, NXC, NDN
      integer  IDXC(NXC,IMAX), KPRF(IMAX), IDM(IMAX,L,NAE)
!
      real (kind=kind_rad) CMIX(NXC,IMAX), DENN(NDN,IMAX), RH(IMAX,L)
! --- OUTPUT
      real (kind=kind_rad)
     &  TAU(IMAX,L),  SSA(IMAX,L),  ASY(IMAX,L)
! --- AEROSOL DATA
      real (kind=kind_rad)
     &  EXT0(NCM1,NBDD),    SCA0(NCM1,NBDD),    SSA0(NCM1,NBDD)
     &, AEXT(NCF,NCM2,NBDD),BSCA(NCF,NCM2,NBDD),ASF0(NCM1,NBDD)
     &, CSSA(NCF,NCM2,NBDD),DASF(NCF,NCM2,NBDD),ABPW(NCM2,NBDD)
     &, ESTR(NBDD)
!
      real (kind=kind_rad) HAER(NDM,NAE),  DZ(IMAX,L),  HZ(IMAX,L+1)
!
      real (kind=kind_rad) crt1, crt2
      DATA  CRT1,CRT2 / 30.0, 0.03333 /
      SAVE CRT1, CRT2
     &,    EXT0,SCA0,SSA0,AEXT,BSCA,ASF0,CSSA,DASF,ABPW,ESTR
!
      integer i, j, k, idom, kpf, icmp, ib, ic, ic1
!
! --- EXTINCTION COEFFS OF 6 NONE RH DEP COMPNTS IN 12 BANDS
!         INSO      SOOT      MINM      MIAM      MICM      MITR
      DATA  EXT0
     & /8.052E-03,1.356E-06,1.032E-04,2.821E-03,7.597E-02,5.445E-03,
     &  8.076E-03,1.313E-06,1.023E-04,2.835E-03,7.608E-02,5.465E-03,
     &  8.060E-03,1.342E-06,1.029E-04,2.826E-03,7.601E-02,5.452E-03,
     &  8.112E-03,1.248E-06,1.010E-04,2.857E-03,7.625E-02,5.496E-03,
     &  8.148E-03,1.173E-06,9.940E-05,2.879E-03,7.641E-02,5.527E-03,
     &  8.156E-03,1.154E-06,9.895E-05,2.884E-03,7.645E-02,5.534E-03,
     &  8.282E-03,8.612E-07,9.008E-05,2.969E-03,7.703E-02,5.650E-03,
     &  8.524E-03,5.432E-07,6.925E-05,3.140E-03,7.811E-02,5.876E-03,
     &  6.435E-03,7.664E-08,3.429E-06,2.452E-03,9.026E-02,6.150E-03,
     &  9.062E-03,1.471E-07,1.413E-05,3.365E-03,8.368E-02,6.691E-03,
     &  9.021E-03,2.626E-07,3.506E-05,3.411E-03,8.040E-02,6.348E-03,
     &  8.823E-03,2.132E-07,2.628E-05,3.325E-03,8.214E-02,6.442E-03/
!
! --- SCATTERING COEFFS OF 6 NONE RH DEP COMPNTS IN 12 BANDS
!         INSO      SOOT      MINM      MIAM      MICM      MITR
      DATA  SCA0
     & /4.447E-03,4.177E-07,8.264E-05,1.625E-03,4.142E-02,3.034E-03,
     &  4.723E-03,4.061E-07,8.314E-05,1.662E-03,4.150E-02,3.080E-03,
     &  4.539E-03,4.138E-07,8.281E-05,1.637E-03,4.145E-02,3.049E-03,
     &  5.136E-03,3.887E-07,8.389E-05,1.718E-03,4.162E-02,3.149E-03,
     &  5.404E-03,3.623E-07,8.453E-05,1.795E-03,4.184E-02,3.252E-03,
     &  5.423E-03,3.540E-07,8.464E-05,1.817E-03,4.192E-02,3.285E-03,
     &  5.729E-03,2.301E-07,8.277E-05,2.178E-03,4.367E-02,3.840E-03,
     &  6.255E-03,1.129E-07,6.675E-05,2.750E-03,5.242E-02,4.909E-03,
     &  5.553E-03,7.201E-10,2.423E-06,2.140E-03,6.357E-02,5.264E-03,
     &  7.757E-03,6.464E-09,1.348E-05,3.178E-03,6.539E-02,6.172E-03,
     &  7.229E-03,2.859E-08,3.417E-05,3.205E-03,6.186E-02,5.805E-03,
     &  7.253E-03,1.972E-08,2.546E-05,3.119E-03,6.317E-02,5.886E-03/
!
! --- SING SCTR ALBEDOES OF 6 NONE RH DEP COMPNTS IN 12 BANDS
!         INSO      SOOT      MINM      MIAM      MICM      MITR
      DATA  SSA0
     & /5.524E-01,3.081E-01,8.004E-01,5.760E-01,5.452E-01,5.571E-01,
     &  5.846E-01,3.095E-01,8.125E-01,5.860E-01,5.455E-01,5.634E-01,
     &  5.631E-01,3.086E-01,8.044E-01,5.793E-01,5.453E-01,5.592E-01,
     &  6.329E-01,3.116E-01,8.307E-01,6.011E-01,5.459E-01,5.728E-01,
     &  6.631E-01,3.084E-01,8.509E-01,6.230E-01,5.476E-01,5.884E-01,
     &  6.647E-01,3.061E-01,8.560E-01,6.298E-01,5.483E-01,5.935E-01,
     &  6.918E-01,2.672E-01,9.189E-01,7.336E-01,5.670E-01,6.797E-01,
     &  7.336E-01,2.025E-01,9.651E-01,8.748E-01,6.708E-01,8.345E-01,
     &  8.612E-01,8.972E-03,7.039E-01,8.739E-01,7.035E-01,8.567E-01,
     &  8.565E-01,4.103E-02,9.460E-01,9.438E-01,7.816E-01,9.224E-01,
     &  8.009E-01,1.024E-01,9.733E-01,9.395E-01,7.693E-01,9.142E-01,
     &  8.228E-01,7.671E-02,9.463E-01,9.367E-01,7.693E-01,9.133E-01/
!
! --- ASYMMETRY FACTORS OF 6 NONE RH DEP COMPNTS IN 12 BANDS
!         INSO      SOOT      MINM      MIAM      MICM      MITR
      DATA  ASF0
     & /9.390E-01,5.020E-01,7.320E-01,9.030E-01,9.460E-01,9.310E-01,
     &  9.219E-01,4.873E-01,7.272E-01,8.946E-01,9.469E-01,9.256E-01,
     &  9.333E-01,4.971E-01,7.304E-01,9.002E-01,9.463E-01,9.292E-01,
     &  8.963E-01,4.653E-01,7.200E-01,8.820E-01,9.482E-01,9.175E-01,
     &  8.798E-01,4.468E-01,7.126E-01,8.668E-01,9.484E-01,9.064E-01,
     &  8.787E-01,4.437E-01,7.109E-01,8.627E-01,9.481E-01,9.031E-01,
     &  8.600E-01,3.960E-01,6.880E-01,8.030E-01,9.400E-01,8.500E-01,
     &  8.280E-01,3.306E-01,6.616E-01,7.342E-01,8.923E-01,7.742E-01,
     &  8.989E-01,7.289E-02,3.909E-01,6.921E-01,8.067E-01,7.122E-01,
     &  8.136E-01,1.494E-01,5.319E-01,6.878E-01,8.103E-01,6.977E-01,
     &  7.864E-01,2.302E-01,6.087E-01,6.942E-01,8.406E-01,7.178E-01,
     &  8.045E-01,1.937E-01,5.688E-01,6.920E-01,8.282E-01,7.110E-01/
!
! --- FITTING COEFFS OF EXT OF 4 RH DEP COMPNTS IN 12 BANDS
!               WASO/SSCM                      SSAM/SUSO
      DATA (((AEXT(I,J,K),I=1,NCF),J=1,NCM2),K=1,6)
     & /1.595E-05,1.330E-05,4.237E-10, 2.304E-03,2.784E-03,1.272E-09,
     &  1.448E-01,1.770E-01,9.004E-08, 2.309E-04,2.389E-04,7.845E-10,
     &  1.528E-05,1.268E-05,4.174E-10, 2.320E-03,2.798E-03,1.275E-09,
     &  1.451E-01,1.773E-01,9.000E-08, 2.302E-04,2.409E-04,7.955E-10,
     &  1.573E-05,1.310E-05,4.216E-10, 2.309E-03,2.789E-03,1.273E-09,
     &  1.449E-01,1.771E-01,9.003E-08, 2.307E-04,2.396E-04,7.882E-10,
     &  1.426E-05,1.176E-05,4.079E-10, 2.344E-03,2.820E-03,1.279E-09,
     &  1.455E-01,1.778E-01,8.994E-08, 2.292E-04,2.440E-04,8.120E-10,
     &  1.331E-05,1.090E-05,3.976E-10, 2.368E-03,2.842E-03,1.282E-09,
     &  1.458E-01,1.782E-01,8.997E-08, 2.272E-04,2.461E-04,8.286E-10,
     &  1.312E-05,1.073E-05,3.951E-10, 2.374E-03,2.848E-03,1.283E-09,
     &  1.459E-01,1.782E-01,9.000E-08, 2.265E-04,2.462E-04,8.322E-10/
      DATA (((AEXT(I,J,K),I=1,NCF),J=1,NCM2),K=7,NBDD)
     & /1.011E-05,8.090E-06,3.512E-10, 2.468E-03,2.964E-03,1.291E-09,
     &  1.470E-01,1.805E-01,8.984E-08, 2.115E-04,2.436E-04,8.894E-10,
     &  6.646E-06,5.898E-06,6.175E-11, 2.609E-03,3.191E-03,1.322E-09,
     &  1.485E-01,1.805E-01,9.088E-08, 1.718E-04,2.079E-04,9.422E-10,
     &  4.815E-07,8.023E-07,3.772E-12, 1.482E-03,2.153E-03,1.589E-09,
     &  1.673E-01,1.966E-01,9.462E-08, 2.483E-05,3.554E-05,1.950E-11,
     &  9.357E-07,7.900E-07,5.441E-12, 2.246E-03,3.217E-03,1.599E-09,
     &  1.580E-01,1.893E-01,9.297E-08, 4.494E-05,6.263E-05,3.156E-11,
     &  2.732E-06,2.527E-06,1.241E-11, 2.666E-03,3.547E-03,1.433E-09,
     &  1.524E-01,1.844E-01,9.179E-08, 9.999E-05,1.327E-04,1.971E-10,
     &  2.007E-06,1.856E-06,9.589E-12, 2.439E-03,3.323E-03,1.496E-09,
     &  1.553E-01,1.868E-01,9.237E-08, 7.643E-05,1.000E-04,1.695E-10/
!
! --- FITTING COEFFS OF SCA OF 4 RH DEP COMPNTS IN 12 BANDS
!               WASO/SSCM                      SSAM/SUSO
      DATA (((BSCA(I,J,K),I=1,NCF),J=1,NCM2),K=1,6)
     & /1.431E-05,1.314E-05,4.229E-10, 2.303E-03,2.784E-03,1.272E-09,
     &  1.448E-01,1.771E-01,8.997E-08, 2.309E-04,2.389E-04,7.845E-10,
     &  1.400E-05,1.257E-05,4.168E-10, 2.319E-03,2.798E-03,1.275E-09,
     &  1.451E-01,1.774E-01,8.996E-08, 2.302E-04,2.409E-04,7.955E-10,
     &  1.421E-05,1.295E-05,4.209E-10, 2.309E-03,2.789E-03,1.273E-09,
     &  1.449E-01,1.772E-01,8.997E-08, 2.307E-04,2.396E-04,7.882E-10,
     &  1.355E-05,1.172E-05,4.077E-10, 2.344E-03,2.820E-03,1.279E-09,
     &  1.455E-01,1.779E-01,8.993E-08, 2.292E-04,2.440E-04,8.120E-10,
     &  1.294E-05,1.091E-05,3.976E-10, 2.368E-03,2.843E-03,1.282E-09,
     &  1.458E-01,1.782E-01,8.997E-08, 2.272E-04,2.461E-04,8.286E-10,
     &  1.276E-05,1.074E-05,3.951E-10, 2.374E-03,2.848E-03,1.283E-09,
     &  1.459E-01,1.782E-01,9.000E-08, 2.265E-04,2.462E-04,8.322E-10/
      DATA (((BSCA(I,J,K),I=1,NCF),J=1,NCM2),K=7,NBDD)
     & /9.919E-06,8.093E-06,3.512E-10, 2.468E-03,2.964E-03,1.291E-09,
     &  1.470E-01,1.805E-01,8.984E-08, 2.115E-04,2.436E-04,8.894E-10,
     &  6.505E-06,5.900E-06,6.173E-11, 2.609E-03,3.191E-03,1.322E-09,
     &  1.485E-01,1.805E-01,9.088E-08, 1.718E-04,2.079E-04,9.422E-10,
     &  1.224E-07,1.088E-07,1.425E-12, 1.088E-03,1.420E-03,1.191E-09,
     &  1.223E-01,1.253E-01,6.140E-08, 9.093E-06,1.152E-05,1.221E-11,
     &  8.152E-07,7.708E-07,5.416E-12, 2.240E-03,3.214E-03,1.588E-09,
     &  1.542E-01,1.867E-01,8.824E-08, 4.481E-05,6.249E-05,3.148E-11,
     &  2.570E-06,2.513E-06,1.240E-11, 2.665E-03,3.547E-03,1.432E-09,
     &  1.518E-01,1.844E-01,9.171E-08, 9.999E-05,1.327E-04,1.971E-10,
     &  1.848E-06,1.798E-06,9.430E-12, 2.411E-03,3.277E-03,1.467E-09,
     &  1.507E-01,1.814E-01,8.835E-08, 7.545E-05,9.867E-05,1.674E-10/
!
! --- FITTING COEFFS OF SSA OF 4 RH DEP COMPNTS IN 12 BANDS
!               WASO/SSCM                      SSAM/SUSO
      DATA (((CSSA(I,J,K),I=1,NCF),J=1,NCM2),K=1,6)
     & /8.820E-01,1.329E-01,8.925E-02, 9.999E-01,2.130E-04,4.523E-05,
     &  9.994E-01,1.319E-03,-8.368E-4, 1.000E+00,0.000E+00,0.000E+00,
     &  9.071E-01,1.059E-01,6.894E-02, 9.999E-01,1.776E-04,-2.168E-5,
     &  9.995E-01,1.062E-03,-6.559E-4, 1.000E+00,0.000E+00,0.000E+00,
     &  8.904E-01,1.239E-01,8.248E-02, 9.999E-01,2.012E-04,2.292E-05,
     &  9.994E-01,1.233E-03,-7.760E-4, 1.000E+00,0.000E+00,0.000E+00,
     &  9.449E-01,6.534E-02,3.848E-02, 1.000E+00,1.244E-04,-1.202E-4,
     &  9.997E-01,6.763E-04,-3.884E-4, 1.000E+00,0.000E+00,0.000E+00,
     &  9.684E-01,3.971E-02,1.988E-02, 1.000E+00,7.586E-05,-1.400E-4,
     &  9.998E-01,3.888E-04,-2.249E-4, 1.000E+00,0.000E+00,0.000E+00,
     &  9.697E-01,3.814E-02,1.904E-02, 1.000E+00,6.654E-05,-1.230E-4,
     &  9.999E-01,3.521E-04,-2.188E-4, 1.000E+00,0.000E+00,0.000E+00/
      DATA (((CSSA(I,J,K),I=1,NCF),J=1,NCM2),K=7,NBDD)
     & /9.790E-01,2.712E-02,1.340E-02, 1.000E+00,0.000E+00,0.000E+00,
     &  1.000E+00,0.000E+00,0.000E+00, 1.000E+00,0.000E+00,0.000E+00,
     &  9.742E-01,3.376E-02,1.754E-02, 1.000E+00,0.000E+00,0.000E+00,
     &  1.000E+00,0.000E+00,0.000E+00, 1.000E+00,0.000E+00,0.000E+00,
     &  4.741E-01,-2.284E-2,4.737E-01, 7.713E-01,-2.185E-1,3.434E-01,
     &  7.285E-01,-2.238E-1,2.247E-01, 4.767E-01,2.866E-01,2.205E-01,
     &  8.438E-01,1.973E-01,1.269E-01, 9.962E-01,5.396E-04,-7.588E-3,
     &  9.782E-01,-7.483E-3,-5.631E-2, 9.936E-01,4.980E-03,6.300E-04,
     &  9.198E-01,1.049E-01,6.285E-02, 9.995E-01,1.023E-03,-7.783E-4,
     &  9.960E-01,7.559E-03,-4.647E-3, 1.000E+00,-1.103E-7,-1.787E-6,
     &  8.678E-01,1.263E-01,1.084E-01, 9.842E-01,-1.268E-2,1.762E-02,
     &  9.728E-01,-1.375E-2,-1.121E-2, 9.651E-01,1.992E-02,1.374E-02/
!
! --- FITTING COEFFS OF ASF OF 4 RH DEP COMPNTS IN 12 BANDS
!               WASO/SSCM                      SSAM/SUSO
      DATA (((DASF(I,J,K),I=1,NCF),J=1,NCM2),K=1,6)
     & /7.240E-01,8.256E-02,3.256E-02, 7.758E-01,1.140E-01,1.584E-02,
     &  8.495E-01,2.280E-02,-1.158E-1, 7.431E-01,7.601E-02,-1.690E-2,
     &  7.178E-01,9.182E-02,3.936E-02, 7.739E-01,1.174E-01,1.304E-02,
     &  8.500E-01,2.590E-02,-1.115E-1, 7.470E-01,7.510E-02,-2.346E-2,
     &  7.220E-01,8.564E-02,3.483E-02, 7.752E-01,1.151E-01,1.491E-02,
     &  8.497E-01,2.383E-02,-1.144E-1, 7.444E-01,7.571E-02,-1.909E-2,
     &  7.085E-01,1.057E-01,4.957E-02, 7.710E-01,1.225E-01,8.826E-03,
     &  8.507E-01,3.054E-02,-1.052E-1, 7.527E-01,7.373E-02,-3.331E-2,
     &  7.013E-01,1.156E-01,5.817E-02, 7.681E-01,1.260E-01,7.093E-03,
     &  8.505E-01,3.608E-02,-9.727E-2, 7.574E-01,7.191E-02,-3.971E-2,
     &  7.003E-01,1.167E-01,5.963E-02, 7.674E-01,1.263E-01,7.395E-03,
     &  8.502E-01,3.756E-02,-9.510E-2, 7.581E-01,7.138E-02,-4.017E-2/
      DATA (((DASF(I,J,K),I=1,NCF),J=1,NCM2),K=7,NBDD)
     & /6.851E-01,1.307E-01,8.454E-02, 7.614E-01,1.239E-01,4.612E-03,
     &  8.486E-01,5.475E-02,-6.403E-2, 7.681E-01,6.680E-02,-4.495E-2,
     &  6.554E-01,1.454E-01,1.182E-01, 7.606E-01,1.143E-01,-2.758E-2,
     &  8.424E-01,6.848E-02,-4.633E-2, 7.632E-01,8.161E-02,-2.990E-2,
     &  3.633E-01,1.483E-01,2.637E-01, 7.303E-01,1.911E-01,8.386E-03,
     &  8.526E-01,1.819E-01,-1.123E-1, 5.014E-01,2.269E-01,2.352E-01,
     &  5.068E-01,1.733E-01,2.345E-01, 7.797E-01,1.068E-01,-9.434E-2,
     &  8.157E-01,1.042E-01,4.413E-03, 6.629E-01,1.697E-01,1.018E-01,
     &  5.899E-01,1.605E-01,1.771E-01, 7.728E-01,9.568E-02,-8.258E-2,
     &  8.304E-01,8.425E-02,-3.121E-2, 7.264E-01,1.228E-01,2.347E-02,
     &  5.476E-01,1.640E-01,2.019E-01, 7.725E-01,1.059E-01,-7.896E-2,
     &  8.272E-01,9.706E-02,-2.339E-2, 6.906E-01,1.454E-01,6.351E-02/
!
! --- POWER FACTOR FOR EXT, SCA FITTING COEFFS OF 4 RH DEP COMPNTS
!        WASO  SSAM  SSCM  SUSO   WASO  SSAM  SSCM  SUSO
      DATA ABPW
     & / 24.0, 33.0, 33.0, 28.0,  24.0, 33.0, 33.0, 28.0,
     &   24.0, 33.0, 33.0, 28.0,  24.0, 33.0, 33.0, 28.0,
     &   24.0, 33.0, 33.0, 28.0,  24.0, 33.0, 33.0, 28.0,
     &   24.0, 33.0, 33.0, 28.0,  27.0, 33.0, 33.0, 28.0,
     &   29.0, 33.0, 33.0, 34.0,  29.0, 33.0, 33.0, 34.0,
     &   29.0, 33.0, 33.0, 31.0,  29.0, 33.0, 33.0, 31.0 /
! --- EXTINGCTION COEFFS IN STRATOSPHERE FOR 12 BANDS
      DATA  ESTR
     & / 3.39E-4, 3.34E-4, 3.38E-4, 3.28E-4, 3.22E-4, 3.18E-4,
     &   3.01E-4, 2.09E-4, 1.70E-5, 5.01E-5, 1.03E-4, 7.72E-5 /
!
      real (kind=kind_rad) drh,  ext1, sca1, ssa1, asf1, drh1, drh2
     &,                     ex00, sc00, ss00, as00, ssa2, asf2, ext2
     &,                     ex01, sc01, ss01, as01
     &,                     ex02, sc02, ss02, as02
     &,                     ex03, sc03, ss03, as03, hd,   sig0u
     &,                     sig0l, ratio, hdi, tt
      if (NBD .NE. NBDD) then
         print *,' IN AEROS NBD =', NBD,' NBDD=',NBDD
         CALL ABORT
      endif
      DO I=1,IMAX
!
        KPF = KPRF(I)
        DO K=1,L
          IDOM = IDM(I,K,KPF)
          DRH = RH(I,K) - 0.5
!
          IF (IDOM .EQ. 1) THEN
! --- 1ST DOMAN - MIXING LAYER
            EXT1 = 0.0
            SCA1 = 0.0
            SSA1 = 0.0
            ASF1 = 0.0
            DO ICMP=1,NXC
              IC = IDXC(ICMP,I)
              IF (IC .GT. NCM1) THEN
                IC1 = IC - NCM1
                DRH1 = EXP(ABPW(IC1,IB)*DRH)
                DRH2 = DRH * DRH
                EX00 = AEXT(1,IC1,IB) + AEXT(2,IC1,IB)*DRH
     &               + AEXT(3,IC1,IB)*DRH1
                SC00 = BSCA(1,IC1,IB) + BSCA(2,IC1,IB)*DRH
     &               + BSCA(3,IC1,IB)*DRH1
                SS00 = CSSA(1,IC1,IB) + CSSA(2,IC1,IB)*DRH
     &               + CSSA(3,IC1,IB)*DRH2
                AS00 = DASF(1,IC1,IB) + DASF(2,IC1,IB)*DRH
     &               + DASF(3,IC1,IB)*DRH2
              ELSE IF (IC .GT. 0) THEN
                EX00 = EXT0(IC,IB)
                SC00 = SCA0(IC,IB)
                SS00 = SSA0(IC,IB)
                AS00 = ASF0(IC,IB)
              ELSE
                EX00 = 0.0
                SC00 = 0.0
                SS00 = 0.0
                AS00 = 0.0
              END IF
              EXT1 = EXT1 + CMIX(ICMP,I) * EX00
              SCA1 = SCA1 + CMIX(ICMP,I) * SC00
              SSA1 = SSA1 + CMIX(ICMP,I) * SS00 * EX00
              ASF1 = ASF1 + CMIX(ICMP,I) * AS00 * SC00
            END DO
            EXT2 = EXT1 * DENN(1,I)
            SSA2 = SSA1 / EXT1
            ASF2 = ASF1 / SCA1
          ELSE IF (IDOM .EQ. 2) THEN
! --- 2ND DOMAIN - MINERAL TRANSPORT LAYERS
            EXT2 = EXT0(6,IB) * DENN(2,I)
            SSA2 = SSA0(6,IB)
            ASF2 = ASF0(6,IB)
          ELSE IF (IDOM .EQ. 3) THEN
! --- 3RD DOMAIN - FREE TROPOSPHERIC LAYERS
!   1:INSO 0.17E-3; 2:SOOT 0.4; 7:WASO 0.59983; N:730
            DRH1 = EXP(ABPW(1,IB)*DRH)
            DRH2 = DRH * DRH
            EX01 = EXT0(1,IB)
            SC01 = SCA0(1,IB)
            SS01 = SSA0(1,IB)
            AS01 = ASF0(1,IB)
            EX02 = EXT0(2,IB)
            SC02 = SCA0(2,IB)
            SS02 = SSA0(2,IB)
            AS02 = ASF0(2,IB)
            EX03 = AEXT(1,1,IB) + AEXT(2,1,IB)*DRH + AEXT(3,1,IB)*DRH1
            SC03 = BSCA(1,1,IB) + BSCA(2,1,IB)*DRH + BSCA(3,1,IB)*DRH1
            SS03 = CSSA(1,1,IB) + CSSA(2,1,IB)*DRH + CSSA(3,1,IB)*DRH2
            AS03 = DASF(1,1,IB) + DASF(2,1,IB)*DRH + DASF(3,1,IB)*DRH2
            EXT1 = 0.17E-3*EX01 + 0.4*EX02 + 0.59983*EX03
            SCA1 = 0.17E-3*SC01 + 0.4*SC02 + 0.59983*SC03
            SSA1 = 0.17E-3*SS01*EX01 + 0.4*SS02*EX02 + 0.59983*SS03*EX03
            ASF1 = 0.17E-3*AS01*SC01 + 0.4*AS02*SC02 + 0.59983*AS03*SC03
            EXT2 = EXT1 * 730.0
            SSA2 = SSA1 / EXT1
            ASF2 = ASF1 / SCA1
          ELSE IF (IDOM .EQ. 4) THEN
! --- 4TH DOMAIN - STRATOSPHERIC LAYERS
            EXT2 = ESTR(IB)
            SSA2 = 0.9
            ASF2 = 0.6
          ELSE
! --- UPPER STRATOSPHERE ASSUME NO AEROSOL
            EXT2 = 0.0
            SSA2 = 1.0
            ASF2 = 0.0
          END IF
!
          HD = HAER(IDOM,KPF)
          IF (HD .GT. 0.0E0) THEN
            HDI      = 1.0 / HD
            SIG0U    = EXP(-HZ(I,K)  *HDI)
            SIG0L    = EXP(-HZ(I,K+1)*HDI)
            TAU(I,K) = EXT2 * HD*(SIG0L - SIG0U)
          ELSE
            TAU(I,K) = EXT2 * DZ(I,K)
!           TAU(I,K) = (EXT2-HD*HH(K)*ALOG(0.5*(SIG0U+SIG0L)))*DZ(K)
          END IF
          SSA(I,K) = SSA2
          ASY(I,K) = ASF2
!         write(6,112) IB,K,I,IDOM,HD,HH(K),DRH,DZ(K),SIG0U,SIG0L,
!    &                 DENN(1,I),DENN(2,I),EXT2,TAU(I,K),SSA2,ASF2
!112      format(1x,'IB,K,I=',3i3,' IDOM,HD,HH,DRH=',i2,3f5.2,
!    &           ' DZ,SIG0U,SIG0L=',3f6.3,/' DENN=',2f8.2,
!    &           ' EXT2,TAU,SSA,ASF=',4f6.3)
        END DO
      END DO
!
!===> ... SMOOTH PROFILE AT DOMAIN BOUNDARIES
!
      DO K=2,L
        DO I=1,IMAX
          RATIO = 1.0
          IF (TAU(I,K) .GT. 0.0) RATIO = TAU(I,K-1) / TAU(I,K)
          TT = TAU(I,K) + TAU(I,K-1)
          IF (RATIO .GT. CRT1) THEN
            TAU(I,K) = 0.2 * TT
            TAU(I,K-1) = TT - TAU(I,K)
          ELSE IF (RATIO .LT. CRT2) THEN
            TAU(I,K) = 0.8 * TT
            TAU(I,K-1) = TT - TAU(I,K)
          END IF
        ENDDO
      ENDDO
!
      RETURN
      END
