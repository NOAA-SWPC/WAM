 
      SUBROUTINE ASTRONOMY(lonl2,latd,nlats,lons_lar,sinlat,coslat,
     &  xlon,fhswr,idate,phour,lsswr,lslwr,
     &  SOLC,RSIN1,RCOS1,RCOS2,slag,sdec,cdec,COSZEN,coszdg,
     &  global_lats_r)
 
      USE MACHINE , ONLY :kind_rad

cc
      use resol_def
      use layout1
      implicit none
 
      integer  global_lats_r(latr)
      integer latd,nlats,lons_lar(latr),idate(4),lonl2
      integer JDNMC,kyear,jd,IMON,IDAY,IZTIM,IHR
      integer IM,ID,IYEAR
      integer                 lonsperlar(latr)
      logical lsswr,lslwr
      character*4 munth
c
      integer loz,jmr,jmout
      parameter (jmr=18,loz=17,jmout=37)
c
      real(kind=kind_rad) O3CLIM(JMR,LOZ,12),
     1                  o3out(jmout,loz),pstr(loz)
c
      real (kind=kind_rad) sinlat(latr),coslat(latr),xlon(LONL2,latd)
      real (kind=kind_rad) fhswr,phour
      real (kind=kind_rad) slag,sdec,cdec,solhr,fjdnmc,sc
      real (kind=kind_rad) COSZDG(LONL2,latd),COSZEN(LONL2,latd)
      real (kind=kind_rad) FJD,DLT,R1,ALF,XMIN
      real (kind=kind_rad) SOLC,RSIN1,RCOS1,RCOS2
cc
cc--------------------------------------------------------------------
cc
      real(kind=kind_rad), parameter ::  cons_24=24.0d0
cc
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
      SOLHR=MOD(PHOUR+IDATE(1),cons_24)
 
!sela if(me.eq.0) PRINT 1001, JCAP, LEVS
 1001 FORMAT (1H0,'GFDL/HOU REDCRAD',I2,I2,'became oper. june 15 1998')
C
C    ****************************************************************
C... * ASTRONOMY CALCULATIONS-ONCE FOR EACH NEW RADIATION interval  *
C    ****************************************************************
C..      GET 4 DIGIT YEAR FOR JULIAN DAY COMPUTATION
      kyear = IDATE(4)
      IMON  = IDATE(2)
      IDAY  = IDATE(3)
      IZTIM = IDATE(1)
      CALL COMPJD(KYEAR,IMON,IDAY,IZTIM,0,JDNMC,FJDNMC)
      CALL FCSTIM(PHOUR,IMON,IDAY,IZTIM,JDNMC,FJDNMC,
     1            RSIN1,RCOS1,RCOS2,JD,FJD)
C..**************************
      IF(lsswr) THEN
        CALL SOLAR(JD,FJD,R1,DLT,ALF,SLAG,SDEC,CDEC)
c        if(me.eq.0)print*,'in astronomy completed sr solar'
        CALL COSZMN(fhswr,SOLHR,SINLAT,COSLAT,SDEC,CDEC,SLAG,
     &    XLON,LONL2,latd,COSZEN,.TRUE.,COSZDG,nlats,lons_lar,
     &    global_lats_r)
c        if(me.eq.0)print*,'in astronomy completed sr coszmn'
C
C       CALCULATE SOLAR INPUT APPROPRIATE FOR DATE
        sc=2.
        SOLC=SC/(R1*R1)
      ENDIF
 
      CALL CDATE(JD,FJD,MUNTH,IM,ID,IYEAR,IHR,XMIN)
c        if(me.eq.0)print*,'in astronomy completed sr cdate'
!JFM  IF (me.eq.0) CALL PRTIME(ID,MUNTH,IYEAR,IHR,XMIN,
!JFM &                     JD,FJD,DLT,ALF,R1,SLAG,SOLC)
c        if(me.eq.0)print*,'in astronomy completed sr prtime'
c
!JFM  call o3intpnasa(phour,idate,o3clim,pstr,o3out)
c        if(me.eq.0)print*,'completed sr o3intpnasa and astronomy!!!'
 
      RETURN
      END
c
c***********************************************************************
c
      SUBROUTINE CDATE(JD,FJD,MUNTH,IM,ID,IYEAR,IHR,XMIN)
CFPP$ NOCONCUR R
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .                                       .
C SUBPROGRAM:    CDATE       COMPUTES DAY,MONTH,YR FROM JULIAN DAY
C   PRGMMR: KENNETH CAMPANA  ORG: W/NMC23    DATE: 89-07-07
C
C ABSTRACT: COMPUTES MONTH,DAY,YEAR FROM JULIAN DAY.
C
C PROGRAM HISTORY LOG:
C   77-06-07  ROBERT WHITE,GFDL
C   98-05-15  IREDELL   Y2K COMPLIANCE
C
C USAGE:    CALL CDATE(JD,FJD,MUNTH,IM,ID,IYEAR,IHR,XMIN)
C   INPUT ARGUMENT LIST:
C     JD       - JULIAN DAY FOR CURRENT FCST HOUR.
C     FJD      - FRACTION OF THE JULIAN DAY.
C   OUTPUT ARGUMENT LIST:
C     MUNTH    - MONTH (CHARACTER).
C     IM       - MONTH (INTEGER).
C     ID       - DAY OF THE MONTH.
C     IYEAR    - YEAR.
C     IHR      - HOUR OF THE DAY.
C     XMIN     - MINUTE OF THE HOUR.
C
C SUBPROGRAMS CALLED:
C   W3FS26     YEAR, MONTH, DAY FROM JULIAN DAY NUMBER
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN.
C
C$$$
      USE MACHINE , ONLY :kind_rad
      implicit none
 
      character*4 month(12),munth
      integer JD,IM,ID,IYEAR,IHR
      integer jda,mfjd,idaywk,idayyr
 
      real (kind=kind_rad) fjd,xmin
 
      DATA MONTH /'JAN.','FEB.','MAR.','APR.','MAY ','JUNE',
     &            'JULY','AUG.','SEP.','OCT.','NOV ','DEC.'/
      IF(FJD.GE.0.5) THEN
        JDA=JD+1
        MFJD=NINT(FJD*1440.)
        IHR=MFJD/60-12
        XMIN=MFJD-(IHR+12)*60
      ELSE
        JDA=JD
        MFJD=NINT(FJD*1440.)
        IHR=MFJD/60+12
        XMIN=MFJD-(IHR-12)*60
      ENDIF
      CALL W3FS26(JDA,IYEAR,IM,ID,IDAYWK,IDAYYR)
      MUNTH=MONTH(IM)
      END
c
c***********************************************************************
c
      SUBROUTINE COMPJD(JYR,JMNTH,JDAY,JHR,JMN,JD,FJD)
CFPP$ NOCONCUR R
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .                                       .
C SUBPROGRAM:    COMPJD      COMPUTES JULIAN DAY AND FRACTION
C   PRGMMR: KENNETH CAMPANA  ORG: W/NMC23    DATE: 89-07-07
C
C ABSTRACT: COMPUTES JULIAN DAY AND FRACTION
C   FROM YEAR, MONTH, DAY AND TIME UTC.
C
C PROGRAM HISTORY LOG:
C   77-05-06  RAY ORZOL,GFDL
C   98-05-15  IREDELL   Y2K COMPLIANCE
C
C USAGE:    CALL COMPJD(JYR,JMNTH,JDAY,JHR,JMN,JD,FJD)
C   INPUT ARGUMENT LIST:
C     JYR      - YEAR (4 DIGITS)-INTIAL FCST TIME.
C     JMNTH    - MONTH-INITIAL FCST TIME.
C     JDAY     - DAY-INITIAL FCST TIME.
C     JHR      - Z-TIME OF INITIAL FCST TIME.
C     JMN      - MINUTES (ZERO PASSED FROM CALLING PROGRAM).
C   OUTPUT ARGUMENT LIST:
C     JD       - JULIAN DAY.
C     FJD      - FRACTION OF THE JULIAN DAY.
C
C SUBPROGRAMS CALLED:
C   IW3JDN     COMPUTE JULIAN DAY NUMBER
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN.
C
C$$$
      USE MACHINE , ONLY :kind_rad
      implicit none
!
      integer JYR,JMNTH,JDAY,JHR,JMN,JD
      integer IW3JDN
 
      real (kind=kind_rad) FJD
 
      JD=IW3JDN(JYR,JMNTH,JDAY)
      IF(JHR.LT.12) THEN
        JD=JD-1
        FJD=0.5+JHR/24.+JMN/1440.
      ELSE
        FJD=(JHR-12)/24.+JMN/1440.
      ENDIF
      END
c
c***********************************************************************
c
      SUBROUTINE COSZMN(DTSWAV,SOLHR,SINLAT,COSLAT,SDEC,CDEC,SLAG,
     1             XLON,NLON2,latd,COSZEN,LDG,COSZDG,nlats,lons_lar,
     &             global_lats_r)
c
c***********************************************************************
c
C===>  COMPUTE MEAN COS SOLAR ZEN ANGL OVER DTSWAV HRS
C....   COSINE OF SOLAR ZEN ANGL FOR BOTH N. AND S. HEMISPHERES.
C        SOLHR=TIME(HRS) AFTER 00Z (GREENWICH TIME)..
C        XLON IS EAST LONG(RADIANS)..
C        SINLAT, COSLAT ARE SIN AND COS OF LATITUDE (N. HEMISPHERE)
C        SDEC, CDEC = THE SINE AND COSINE OF THE SOLAR DECLINATION.
C        SLAG = EQUATION OF TIME
C
      USE MACHINE , ONLY :kind_rad

      use resol_def
      use layout1
      implicit none 
      integer              global_lats_r(latr)
      integer nlats,lons_lar(latr)
      integer NLON2,latd
      integer ISTSUN(NLON2)
      integer nstp,istp
      integer i,it,j,nlon,lat
      LOGICAL LDG
      real (kind=kind_rad) DTSWAV,SOLHR,SDEC,CDEC,SLAG
      real (kind=kind_rad) XLON(NLON2,latd),COSZEN(NLON2,latd)
      real (kind=kind_rad) COSZDG(NLON2,latd)
      real (kind=kind_rad) SINLAT(latr),COSLAT(latr),COSZN(NLON2)
      real (kind=kind_rad) PID12,CNS,SS,CC
cc
cc--------------------------------------------------------------------
cc
      real(kind=kind_rad), parameter :: cons_0=0.0d0
cc
c @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
      NLON=NLON2/2
      NSTP = 6
      ISTP = 1 ! jbao NSTP*DTSWAV
      PID12 = (2.E0 * ASIN(1.E0)) / 12.E0
 
      DO j=1,nlats
        DO i=1,NLON2
          COSZEN(i,j) = 0.E0
          ISTSUN(i) = 0
        ENDDO
        DO IT=1,ISTP
          CNS = PID12 * (SOLHR-12.E0+(IT-1)*1.E0/NSTP) +SLAG
          lat = global_lats_r(ipt_lats_node_r-1+j)
          SS= SINLAT(lat)*SDEC
          CC= COSLAT(lat)*CDEC
cjfe      DO i=1,lonsinpe(0,j)
          DO i=1,lons_lar(lat)
            COSZN(i) = SS + CC * COS(CNS + XLON(i,j))
            COSZEN(i,j) = COSZEN(i,j) +  MAX (cons_0, COSZN(i))
            IF(COSZN(i).GT.0.E0) ISTSUN(i) = ISTSUN(i) + 1
          ENDDO
cjfe      SS=-SS
cjfe      DO I=NLON+1,NLON+lonsinpe(0,j)
cjfe        COSZN(i) = SS + CC * COS(CNS + XLON(i,j))
cjfe        COSZEN(i,j) = COSZEN(i,j) + AMAX1(0.E0, COSZN(i))
cjfe        IF(COSZN(i).GT.0.E0) ISTSUN(i) = ISTSUN(i) + 1
cjfe      ENDDO
        ENDDO
        DO i=1,NLON2
          IF(LDG) COSZDG(i,j) = COSZEN(i,j) / ISTP
          IF(ISTSUN(i).GT.0) COSZEN(i,j) = COSZEN(i,j) / ISTSUN(i)
        ENDDO
      ENDDO
 
      RETURN
      END
c
c***********************************************************************
c
      SUBROUTINE FCSTIM(FHOUR,IMON,IDAY,IZTIM,JDNMC,FJDNMC,
     1                  RSIN1,RCOS1,RCOS2,JD,FJD)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .                                       .
C SUBPROGRAM:    FCSTIM      SET FORECAST ORBIT PARMS AND JULIAN DAY.
C   PRGMMR: KENNETH CAMPANA  ORG: W/NMC23    DATE: 89-07-07
C
C ABSTRACT: FOR A GIVEN FORECAST HOUR AND INITIAL JULIAN DAY,
C   THREE ORBIT PARAMETERS AND THE FORECAST JULIAN DAY ARE COMPUTED.
C
C PROGRAM HISTORY LOG:
C   98-05-15  IREDELL   Y2K COMPLIANCE
C
C USAGE:    CALL FCSTIM(FHOUR,IMON,IDAY,IZTIM,JDNMC,FJDNMC,
C    1                  RSIN1,RCOS1,RCOS2,JD,FJD)
C   INPUT ARGUMENT LIST:
C     FHOUR    - FORECAST HOUR
C     IMON     - NOT USED
C     IDAY     - NOT USED
C     IZTIM    - NOT USED
C     JDNMC    - INITIAL JULIAN DAY.
C     FJDNMC   - INITIAL FRACTION OF THE JULIAN DAY.
C     RLAG     - DAY OF PERIHELION?
C     YEAR     - DAYS IN YEAR
C   OUTPUT ARGUMENT LIST:
C     RSIN1    - ORBIT PARAMETER
C     RCOS1    - ORBIT PARAMETER
C     RCOS2    - ORBIT PARAMETER
C     JD       - FORECAST JULIAN DAY.
C     FJD      - FORECAST FRACTION OF THE JULIAN DAY.
C
C SUBPROGRAMS CALLED:
C   W3FS26     YEAR, MONTH, DAY FROM JULIAN DAY NUMBER
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN.
C
C$$$
      USE MACHINE , ONLY :kind_rad
      implicit none
 
      integer IMON,IDAY,IZTIM,JDNMC,JD
      integer JDA,IYEAR,IM,ID,IDAYWK,IDAYYR
      real (kind=kind_rad) TPI
      PARAMETER (TPI=2.E0*3.141593E+0)
 
      real (kind=kind_rad) FHOUR,FJDNMC,RSIN1,RCOS1,RCOS2,FJD
      real (kind=kind_rad) XDA,DYINC,DYFCST,rlag,year,RANG
 
      IF(FJDNMC.GE.0.5) THEN
        JDA=JDNMC+1
        XDA=FJDNMC-0.5
      ELSE
        JDA=JDNMC
        XDA=FJDNMC+0.5
      ENDIF
      CALL W3FS26(JDA,IYEAR,IM,ID,IDAYWK,IDAYYR)
      DYINC=FHOUR/24
      DYFCST=IDAYYR+XDA+DYINC
      rlag=14.8125
      year=365.25
      RANG=TPI*(DYFCST-RLAG)/YEAR
      RSIN1=SIN(RANG)
      RCOS1=COS(RANG)
      RCOS2=COS(2*RANG)
      JD=JDNMC+FJDNMC+DYINC
      FJD=JDNMC+FJDNMC+DYINC-JD
      END
c
c***********************************************************************
c
      SUBROUTINE PRTIME(ID,MUNTH,IYEAR,IHR,XMIN,JD,FJD,
     1                  DLT,ALF,R1,SLAG,SOLC)
 
      use machine
      use physcons, pi => con_pi
      implicit none
 
      character*4 munth
      integer ID,IYEAR,IHR,JD
      integer LTD,LTM,IHALP,IYY
      integer SIGN,SIGB,DSIG
      real (kind=kind_rad) XMIN,FJD,DLT,ALF,R1,SLAG,SOLC
      real (kind=kind_rad) DEGRAD,HPI,ZERO,SIX,SIXTY,Q22855
      real (kind=kind_rad) DLTD,DLTM,DLTS,HALP,YMIN,ASEC,EQT,EQSEC
 
      PARAMETER (DEGRAD=180.E0/PI,HPI=0.5E0*PI)
      DATA       SIGN/1H-/,      SIGB/1H /
      DATA ZERO,SIX,SIXTY,Q22855/0.0,6.0,60.0,228.55735/
      SAVE SIGN,ZERO,SIX,SIXTY,Q22855
      DLTD=DEGRAD*DLT
      LTD=DLTD
      DLTM=SIXTY*(ABS(DLTD)-ABS(FLOAT(LTD)))
      LTM=DLTM
      DLTS=SIXTY*(DLTM-FLOAT(LTM))
      DSIG=SIGB
      IF((DLTD.LT.ZERO).AND.(LTD.EQ.0)) DSIG=SIGN
      HALP=SIX*ALF/HPI
      IHALP=HALP
      YMIN=ABS(HALP-FLOAT(IHALP))*SIXTY
      IYY=YMIN
      ASEC=(YMIN-FLOAT(IYY))*SIXTY
      EQT=Q22855*SLAG
      EQSEC=SIXTY*EQT
!jbao      PRINT 1004,    ID,MUNTH,IYEAR,IHR,XMIN,JD,FJD,R1,HALP,IHALP,
!jbao     1       IYY,ASEC,DLTD,DSIG,LTD,LTM,DLTS,EQT,EQSEC,SLAG,SOLC
 1004 FORMAT('0 FORECAST DATE',9X,I3,A5,I6,' AT',I3,' HRS',F6.2,' MINS'/
     1       '  JULIAN DAY',12X,I8,2X,'PLUS',F11.6/
     2       '  RADIUS VECTOR',9X,F10.7/
     3       '  RIGHT ASCENSION OF SUN',F12.7,' HRS, OR',I4,' HRS',I4,
     4                                 ' MINS',F6.1,' SECS'/
     5       '  DECLINATION OF THE SUN',F12.7,' DEGS, OR',A2,I3,
     6                                 ' DEGS',I4,' MINS',F6.1,' SECS'/
     7       '  EQUATION OF TIME',6X,F12.7,' MINS, OR',F10.2,' SECS, OR'
     8                           ,F9.6,' RADIANS'/
     9       '  SOLAR CONSTANT',8X,F12.7//)
      RETURN
      END
c
c***********************************************************************
c
      SUBROUTINE SOLAR(JD,FJD,R,DLT,ALP,SLAG,SDEC,CDEC)
C
C
C    *******************************************************************
C    *                            S O L A R                            *
C... *  PATTERNED AFTER ORIGINAL GFDL CODE---                          *
C... *     BUT NO CALCULATION OF LATITUDE MEAN COS SOLAR ZENITH ANGLE..*
C... *     ZENITH ANGLE CALCULATIONS DONE IN SR   ZENITH IN THIS CASE..*
C... *  HR ANGLE,MEAN COSZ,AND MEAN TAUDA CALC REMOVED--K.A.C. MAR 89  *
C    *  UPDATES BY HUALU PAN TO LIMIT ITERATIONS IN NEWTON METHOD AND  *
C    *  ALSO CCR REDUCED FROM(1.3E-7)--BOTH TO AVOID NONCONVERGENCE IN *
C    *  NMC S HALF PRECISION VERSION OF GFDL S CODE   ----  FALL 1988  *
C    *******************************************************************
C
C.....SOLAR COMPUTES RADIUS VECTOR, DECLINATION AND RIGHT ASCENSION OF
C.....SUN, EQUATION OF TIME
C
      USE MACHINE , ONLY :kind_rad,kind_phys
      use physcons, pi => con_pi
      implicit none
c 
      real(kind=kind_rad) CYEAR,SVT6,CCR,TPP
      integer JDOR
                                   D A T A
     1   CYEAR/365.25/,      CCR/1.3E-6/
C
C.....TPP = DAYS BETWEEN EPOCH AND PERIHELION PASSAGE OF 1900
C.....SVT6 = DAYS BETWEEN PERIHELION PASSAGE AND MARCH EQUINOX OF 1900
C.....JDOR = JD OF EPOCH WHICH IS JANUARY 0, 1900 AT 12 HOURS UT
C
                                   D A T A
     1   TPP/1.55/,          SVT6/78.035/,       JDOR/2415020/
C
C    *******************************************************************
 
 
      real(kind=kind_rad) TPI,HPI,RAD
      PARAMETER (TPI=2.0*PI,HPI=0.5*PI,RAD=180.0/PI)
      integer JD,JDOE,ITER
      real(kind=kind_rad) FJD,R,DLT,ALP,SLAG,SDEC,CDEC
      real(kind=kind_rad) DAT,T,YEAR,TYEAR,EC,ANGIN,ADOR,DELEQN
      real(kind=kind_rad) SNI,TINI,ER,QQ,E,EP,CD,HE,EQ,DATE
      real(kind=kind_rad) EM,CR,W,TST,SUN
C
      DAT=FLOAT(JD-JDOR)-TPP+FJD
C    COMPUTES TIME IN JULIAN CENTURIES AFTER EPOCH
      T=FLOAT(JD-JDOR)/36525.E0
C    COMPUTES LENGTH OF ANOMALISTIC AND TROPICAL YEARS (MINUS 365 DAYS)
      YEAR=.25964134E0+.304E-5*T
      TYEAR=.24219879E0-.614E-5*T
C    COMPUTES ORBIT ECCENTRICITY AND ANGLE OF EARTH'S INCLINATION FROM T
      EC=.01675104E0-(.418E-4+.126E-6*T)*T
      ANGIN=23.452294E0-(.0130125E0+.164E-5*T)*T
      ADOR=JDOR
      JDOE=ADOR+(SVT6*CYEAR)/(YEAR-TYEAR)
C    DELEQN=UPDATED SVT6 FOR CURRENT DATE
      DELEQN=FLOAT(JDOE-JD)*(YEAR-TYEAR)/CYEAR
      YEAR=YEAR+365.E0
      SNI=SIN(ANGIN/RAD)
      TINI=1.E0/TAN(ANGIN/RAD)
      ER=SQRT((1.E0+EC)/(1.E0-EC))
      QQ=DELEQN*TPI/YEAR
C    DETERMINE TRUE ANOMALY AT EQUINOX
      E=1.E0
      ITER = 0
 32   EP=E-(E-EC*SIN(E)-QQ)/(1.E0-EC*COS(E))
      CD=ABS(E-EP)
      E=EP
      ITER = ITER + 1
      IF(ITER.GT.10) THEN
        WRITE(6,*) ' ITERATION COUNT FOR LOOP 32 =', ITER
        WRITE(6,*) ' E, EP, CD =', E, EP, CD
      ENDIF
      IF(ITER.GT.10) GOTO 1032
      IF(CD.GT.CCR) GO TO 32
 1032 CONTINUE
      HE=.5E0*E
      EQ=2.E0*ATAN(ER*TAN(HE))
C    DATE=DAYS SINCE LAST PERIHELION PASSAGE
      DATE = MOD(DAT,YEAR)
C    SOLVE ORBIT EQUATIONS BY NEWTON'S METHOD
      EM=TPI*DATE/YEAR
      E=1.E0
      ITER = 0
 31   EP=E-(E-EC*SIN(E)-EM)/(1.E0-EC*COS(E))
      CR=ABS(E-EP)
      E=EP
      ITER = ITER + 1
      IF(ITER.GT.10) THEN
        WRITE(6,*) ' ITERATION COUNT FOR LOOP 31 =', ITER
      ENDIF
      IF(ITER.GT.10) GOTO 1031
      IF(CR.GT.CCR) GO TO 31
 1031 CONTINUE
      R=1.E0-EC*COS(E)
      HE=.5E0*E
      W=2.E0*ATAN(ER*TAN(HE))
C>YH  SIND=SNI*SIN(W-EQ)
C>YH  DLT=ASIN(SIND)
      SDEC=SNI*SIN(W-EQ)
      CDEC=SQRT(1.E0 - SDEC*SDEC)
      DLT=ASIN(SDEC)
      ALP=ASIN(TAN(DLT)*TINI)
      TST=COS(W-EQ)
      IF(TST.LT.0.E0) ALP=PI-ALP
      IF(ALP.LT.0.E0) ALP=ALP+TPI
      SUN=TPI*(DATE-DELEQN)/YEAR
      IF(SUN.LT.0.E0) SUN=SUN+TPI
      SLAG=SUN-ALP-.03255E0
      RETURN
      END
c
c***********************************************************************
c
      subroutine o3intpnasa(fhour,idate,o3clim,pstr,o3out)
c     ********************************************************
c     *  COMPUTES O3 CLIMO FROM 12 MONTH DATASET, LINEARLY   *
c     *   INTERPOLATED TO DAY,MON OF THE FCST.  THEN CREATE  *
c     *   A 5 DEG ARRAY FROM THE 10 DEG CLIMATOLOGY...FOR    *
c     *   EASE WHEN DOING A LATITUDINAL INTERPOLATION        *
c     *  THANKS TO S MOORTHI for NEW O3 CLIMO...KAC  DEC 1996*
c     * INPUT:                                               *
c     *   idate=NMC date-time                                *
c     *   fhour=forecast hour                                *
c     *   o3clim=10-deg O3 climo for each month(np->spole)   *
c     * OUTPUT :                                             *
c     *   o3out=5-deg O3 climo for forecast date(np->spole)  *
c     ********************************************************
c GEOS ozone data
      USE MACHINE , ONLY :kind_rad
      implicit none
      integer loz,jmr,jmout
      parameter (jmr=18,loz=17,jmout=37)
C
      INTEGER DAYS(12),idate(4)
      integer ida,imo,numdyz,imo1,jday,nmdtot,ndayr,mday
      integer monL,monC,monR,midL,midC,midR,jmr1,j1,j2
      integer l,j,ken
 
      real(kind=kind_rad) fhour
      real(kind=kind_rad) O3CLIM(JMR,LOZ,12),O3TMP(jmr,loz)
      real(kind=kind_rad) o3out(jmout,loz),pstr(loz)
      real(kind=kind_rad) difL,difR,delday
c
CCCC      common /o3nasaclim/o3clim,pstr,o3out
      DATA  DAYS/31,28,31,30,31,30,31,31,30,31,30,31/
C...     BEGIN HERE  ......
      ida=idate(3)
      imo=idate(2)
c...   FIND current day and month, initial values in ida,imo!
c       will not worry about leap year, since it will take a
c       120-year (WHAT?) forecast to be off by 1 month.  If this
c       is deemed a problem, need to redo this calculation.
c
      if (fhour.ge.24.) then
c... number of days into the forecast
       numdyz=fhour/24.0 + 0.01
c... get day-of-year, remember climate runs are for years
       imo1=imo-1
       jday = ida
       if (imo1.gt.0) then
        jday=0
        do 7 ken=1,imo1
         jday=jday+days(ken)
    7   continue
        jday=jday+ida
       end if
       nmdtot = jday+numdyz
       ndayr = mod(nmdtot,365)
       if (ndayr.eq.0) ndayr=365
c... now get month from day-of-year
       mday=0
       do 8 ken=1,11
        mday=mday+days(ken)
        imo=ken
        if (ndayr.le.mday) then
         ida=ndayr-(mday-days(imo))
         go to 9
        end if
    8  continue
       imo=12
    9  continue
ccc    print 66,fhour,numdyz,jday,nmdtot,ndayr
   66  format(' SBUVO3 climo hr=',f10.1,
     1        ' numdyz,jday,nmdtot,ndayr=',4i8)
      end if
C
c...   do a linear interpolation in time, where we assume that
c       the ozone data is valid for mid-month
c      monL is the preceeding month, monC for current mo, and
c      monR is the future month..
      monL=imo-1
      monC=imo
      monR=imo+1
      if (monL.lt.1) monL=12
      if (monR.gt.12) monR=1
c...    difL=number of days beteen mid-months of the current and
c            preceeding mo, difR=same for current and future mo..
c...    delL=number of days between current day and mon,
c       delR=same for current day and next month.
c       sign convention as if we were using day of year calculations.
      midL=days(monL)/2
      midC=days(monC)/2
      midR=days(monR)/2
      difL=-(days(monL)-midL+midC)
      difR= (days(monC)-midC+midR)
      delday=ida-midC
      if (ida.gt.midC) then
       do 60 j=1,jmr
        do 60 l=1,loz
          O3TMP(j,l)=o3clim(j,l,monC) +
     1        (o3clim(j,l,monR)-o3clim(j,l,monC)) * delday/difR
   60  continue
      else if (ida.lt.midC) then
       do 65 j=1,jmr
        do 65 l=1,loz
          O3TMP(j,l)=o3clim(j,l,monC) +
     1        (o3clim(j,l,monL)-o3clim(j,l,monC)) * delday/difL
   65  continue
      else if (ida.eq.midC) then
       do 70 j=1,jmr
        do 70 l=1,loz
          O3TMP(j,l)=o3clim(j,l,monC)
   70  continue
      end if
!cselaprint 200,imo,ida
c...   linearly interpolate to 5 deg zonal means
      jmr1=jmr-1
      do 80 j=1,jmr1
       j1=j*2
       j2=j1+1
       do 80 l=1,loz
        o3out(j1,l)=O3TMP(j,l)
        o3out(j2,l)=0.5*(O3TMP(j,l)+O3TMP(j+1,l))
   80 continue
      do 85 l=1,loz
       o3out(1,l)=O3TMP(1,l)
       o3out(jmout-1,l)=O3TMP(jmr,l)
       o3out(jmout,l)=O3TMP(jmr,l)
   85 continue
  200 format(1h1,'from o3intpnasa ozone climatology for mont,day=',2i4)
      return
      END
