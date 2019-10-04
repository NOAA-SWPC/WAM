!-----------------------------------------------------------------------
!
      MODULE MODULE_PRECIP_ADJUST
!
! This module contains 3 subroutines:
!     READPCP
!     CHKSNOW
!     ADJPPT
!-----------------------------------------------------------------------
!***
!***  Specify the diagnostic point here: (i,j) and the processor number.
!***  Remember that in WRF, local and global (i,j) are the same, so don't
!***  use the "local(i,j)" output from glb2loc.f; use the GLOBAL (I,J)
!***  and the PE_WRF.
!***
!
      USE MODULE_DM_PARALLEL,ONLY : DSTRB
      INTEGER,PRIVATE,SAVE :: ITS_B1,ITE_B1,JTS_B2,JTE_B2
      
      INTEGER :: ITEST=346,JTEST=256,TESTPE=53
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
      SUBROUTINE READPCP(MYPE,MPI_COMM_COMP                             &
                        ,PPTDAT,DDATA,LSPA,PCPHR                        &
                        ,MY_DOMAIN_ID                                   &
                        ,IDS,IDE,JDS,JDE,LM                             &
                        ,IMS,IME,JMS,JME                                &
                        ,ITS,ITE,JTS,JTE                                &
                        ,ITS_B1,ITE_B1,JTS_B2,JTE_B2)
!
!     ****************************************************************
!     *                                                              *
!     *   PRECIPITATION ASSIMILATION INITIALIZATION.                 *
!     *    READ IN PRECIP ANALYSIS AND DATA MASK AND SET UP ALL      *
!     *    APPROPRIATE VARIABLES.                                    *
!     *                   MIKE BALDWIN, MARCH 1994                   *
!     *                   Adapted to 2-D code, Ying Lin, Mar 1996    *
!     *                   For WRF/NMM: Y.Lin Mar 2005                *
!     *                                                              *
!     ****************************************************************
!-----------------------------------------------------------------------
!
! READ THE BINARY VERSION OF THE PRECIP ANALYSIS.
!

      IMPLICIT NONE
      INTEGER,INTENT(IN) :: MYPE,MPI_COMM_COMP                         &
                           ,MY_DOMAIN_ID                               &
                           ,IDS,IDE,JDS,JDE,LM                         &
                           ,IMS,IME,JMS,JME                            &
                           ,ITS,ITE,JTS,JTE                            &
                           ,ITS_B1,ITE_B1,JTS_B2,JTE_B2
      INTEGER,INTENT(IN) :: PCPHR
      REAL,DIMENSION(IDS:IDE,JDS:JDE) :: TEMPG
      REAL,DIMENSION(IMS:IME,JMS:JME) :: TEMPL
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(OUT) :: DDATA, LSPA
      REAL,DIMENSION(IMS:IME,JMS:JME,1:PCPHR),INTENT(OUT) :: PPTDAT
      INTEGER :: I, IER, IHR, J, N, NUNIT_PCP
      CHARACTER*256 :: MESSAGE
      CHARACTER(14) :: FILENAME
      CHARACTER(6),SAVE :: FMT_ID='(I2.2)'                              &
                          ,FMT_HR='(I1.1)'
      CHARACTER(2) :: CHAR_ID
      CHARACTER(1) :: CHAR_HR
      LOGICAL :: OPENED
!-----------------------------------------------------------------------
!
! Get the value of MYPE:
!
!
      TEMPG=999.
      IF(MYPE==0)THEN
      write(0,*)'PCPHR=',PCPHR
      write(0,*)'IDS,IDE,JDS,JDE in ADJPCP=',IDS,IDE,JDS,JDE
      ENDIF
!
      WRITE(CHAR_ID,FMT_ID)MY_DOMAIN_ID
!
      hours: DO IHR=1,PCPHR
!
        WRITE(CHAR_HR,FMT_HR)IHR
        FILENAME='pcp.hr'//CHAR_HR//'.'//CHAR_ID//'.bin'
!
        IF(MYPE==0)THEN
!
          DO N=51,99
            INQUIRE(N,opened=OPENED)
            IF(.NOT.OPENED)THEN
              NUNIT_PCP=N
              EXIT
            ENDIF
          ENDDO
!
          CLOSE(NUNIT_PCP)
!rv       OPEN(unit=NUNIT_PCP,file=FILENAME,form='UNFORMATTED'        &
!rv           ,STATUS='REPLACE',IOSTAT=IER)
          OPEN(unit=NUNIT_PCP,file=FILENAME,form='UNFORMATTED',IOSTAT=IER)
!rv
          IF(IER/=0)THEN
            WRITE(0,*)' Failed to open ',FILENAME,' in READPCP ier=',IER
          ENDIF
          READ(NUNIT_PCP) ((TEMPG(I,J),I=IDS,IDE),J=JDS,JDE)
!         WRITE(60+IHR,*)((TEMPG(I,J),I=IDS,IDE),J=JDS,JDE)
          WRITE(0,*) 'IHR=', IHR, ' FINISHED READING PCP TO TEMPG'
          CLOSE(NUNIT_PCP)
!
          DO J=JDS,JDE
            DO I=IDS,IDE
! In the binary version of the precip data, missing data are denoted as '999.'
! Convert the valid data from mm to m:

              IF (TEMPG(I,J).LT.900.) TEMPG(I,J)=TEMPG(I,J)*0.001
            ENDDO
          ENDDO
        ENDIF
!
! Distribute to local temp array:
        CALL DSTRB(TEMPG,TEMPL,1,1,1,1,1,MYPE,MPI_COMM_COMP)
!
! Place into correct hour slot in PPTDAT:
        
        DO J=JMS,JME
          DO I=IMS,IME
            PPTDAT(I,J,IHR)=TEMPL(I,J)
          ENDDO
        ENDDO
!
        IF(MYPE==0)THEN
          WRITE(0,*) 'ADJPPT-READPCP, IHR',IHR, 'PPTDAT=',        &
     &      PPTDAT(1,1,IHR)
        ENDIF

      ENDDO hours
!
! Give DDATA (hourly precipitation analysis partitioned into each physics
! timestep; partitioning done in ADJPPT) an initial value of 999, because
! TURBL/SURFCE is called before ADJPPT.  Also initialize LSPA to zero.
!
      DDATA=999.
      LSPA=0.
!
      END SUBROUTINE READPCP
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE CHKSNOW(MYPE,NTSD,DT,NPHS,SR,PPTDAT,PCPHR              &
                        ,IDS,IDE,JDS,JDE,LM                             &
                        ,IMS,IME,JMS,JME                                &
                        ,ITS,ITE,JTS,JTE                                &
                        ,ITS_B1,ITE_B1,JTS_B2,JTE_B2)
!
! AT THE FIRST PHYSICS TIME STEP AFTER THE TOP OF EACH HOUR, CHECK THE SNOW
! ARRAY AGAINST THE SR (SNOW/TOTAL PRECIP RATIO).  IF SR .GE. 0.9, SET THIS
! POINT TO MISSING (SO WE WON'T DO SNOW ADJUSTMENT HERE).
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      INTEGER,INTENT(IN) :: MYPE,NTSD,NPHS
      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,LM                          &
                           ,IMS,IME,JMS,JME                             &
                           ,ITS,ITE,JTS,JTE                             &
                           ,ITS_B1,ITE_B1,JTS_B2,JTE_B2
      INTEGER,INTENT(IN) :: PCPHR
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: SR
      REAL,DIMENSION(IMS:IME,JMS:JME,1:PCPHR),INTENT(INOUT) :: PPTDAT
      REAL,INTENT(IN) :: DT
      REAL :: TIMES
      INTEGER :: I, J, IHR
      CHARACTER*256 :: MESSAGE
!-----------------------------------------------------------------------
      TIMES=NTSD*DT
      IF (MOD(TIMES,3600.) < NPHS*DT) THEN
        IHR=INT(TIMES)/3600+1
        IF (IHR > PCPHR) GO TO 10
        DO J=JTS_B2,JTE_B2
        DO I=ITS_B1,ITE_B1
          IF (SR(I,J) >= 0.9) PPTDAT(I,J,IHR) = 999.
        ENDDO
        ENDDO
!
! Get the value of MYPE:
!
        IF (MYPE==0) THEN
          WRITE(0,1010) TIMES,SR(1,1)
 1010     FORMAT('ADJPPT-CHKSNOW: TIMES, SR=',F6.0,1X,F6.4)
        ENDIF
      ENDIF
 10   CONTINUE
      END SUBROUTINE CHKSNOW
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE ADJPPT(MYPE,NTSD,DT,NPHS,PREC,LSPA,PPTDAT,DDATA,PCPHR  &
                       ,IDS,IDE,JDS,JDE,LM                              &
                       ,IMS,IME,JMS,JME                                 &
                       ,ITS,ITE,JTS,JTE                                 &
                       ,ITS_B1,ITE_B1,JTS_B2,JTE_B2)

!***********************************************************************
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .     
! SUBPROGRAM:    ADJPPT    PRECIPITATION/CLOUD ADJUSTMENT
!    PRGRMMR:    Y. LIN    ORG: W/NP22     DATE: 2005/03/30
!     
! ABSTRACT:
!     ADJPPT  MAKES ADJUSTMENT TO MODEL'S TEMPERATURE, MOISTURE, HYDROMETEOR
!     FIELDS TO BE MORE CONSISTENT WITH THE OBSERVED PRECIPITATION AND CLOUD
!     TOP PRESSURE
!     
!     FOR NOW, AS A FIRST STEP, JUST PARTITION THE INPUT HOURLY PRECIPITATION
!     OBSERVATION INTO TIME STEPS, AND FEED IT INTO THE SOIL.
! PROGRAM HISTORY LOG:
!
!   2005/03/30  LIN      - BAREBONES PRECIPITATION PARTITION/FEEDING TO GROUND
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE : IBM 
!$$$  
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
      INTEGER,INTENT(IN) :: MYPE,NPHS, NTSD
      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,LM                          &
                           ,IMS,IME,JMS,JME                             &
                           ,ITS,ITE,JTS,JTE                             &
                           ,ITS_B1,ITE_B1,JTS_B2,JTE_B2
      REAL,INTENT(IN) :: DT
      INTEGER,INTENT(IN) :: PCPHR
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: PREC
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(OUT) :: DDATA, LSPA
      REAL,DIMENSION(IMS:IME,JMS:JME,1:PCPHR),INTENT(IN) :: PPTDAT
!-----------------------------------------------------------------------
!***
!***  LOCAL VARIABLES
!***
!-----------------------------------------------------------------------
      REAL :: DTPHS, FRACT, FRACT1, FRACT2, TIMES, TPHS1, TPHS2
      INTEGER :: I, J, IHR, IHR1, IHR2, NTSP
      CHARACTER*256 :: MESSAGE
!
! Get the value of MYPE:
!
!
      TIMES=NTSD*DT
      IHR=INT(TIMES)/3600+1
! Size of physics time step:
      DTPHS=NPHS*DT
!
! Compute the beginning and ending time of the current physics time step,
! TPHS1 and TPHS2:
!  
      NTSP=NTSD/NPHS+1
      TPHS1=(NTSP-1)*DTPHS
      TPHS2=NTSP*DTPHS
!
      IHR1=INT(TPHS1)/3600+1
      IHR2=INT(TPHS2)/3600+1
!
! Fraction of an hour that falls into IHR1 and IHR2.  Note that IHR1 and IHR2
! might be identical.
      IF (IHR1 > PCPHR) THEN 
        GO TO 200
      ELSEIF (IHR2 > PCPHR) THEN
        IHR2=3
        FRACT1=(3600.- MOD(INT(TPHS1),3600))/3600.
        FRACT2=0.
      ELSEIF (IHR1 .EQ. IHR2) THEN
         FRACT1=0.5*DTPHS/3600.
         FRACT2=FRACT1
      ELSE
         FRACT1=(3600.- MOD(INT(TPHS1),3600))/3600.
         FRACT2=FLOAT(MOD(INT(TPHS2),3600))/3600.
      ENDIF
!
      FRACT=FRACT1 + FRACT2
!
      IF (MYPE==0) THEN
         WRITE(0,1010) NTSD,NTSP,TIMES,IHR1,IHR2,TPHS1,TPHS2,      &
      &    FRACT1,FRACT2
 1010    FORMAT('ADJPPT: NTSD,NTSP,TIMES=',I4,1X,I4,1X,F6.0,' IHR1,IHR2=' &
      &   ,I1,1X,I1,' TPHS1,TPHS2=',F6.0,1X,F6.0,' FRACT1,FRACT2='        &
      &   ,2(1X,F6.4))
      ENDIF
!
!-----------------------------------------------------------------------
!   FRACT1/2 IS THE FRACTION OF IHR1/2'S PRECIP THAT WE WANT FOR
!   THIS ADJUSTMENT (assuming that the physics time step spans over
!   IHR1 and IHR2.  If not, then IHR1=IHR2).
!-----------------------------------------------------------------------
!   SET UP OBSERVED PRECIP FOR THIS TIMESTEP IN DDATA
!-----------------------------------------------------------------------
!
      DO J=JTS_B2,JTE_B2
      DO I=ITS_B1,ITE_B1
! Note sometimes IHR1=IHR2.  
        IF (PPTDAT(I,J,IHR1).GT.900..OR.PPTDAT(I,J,IHR2).GT.900.) THEN
          DDATA(I,J) = 999.
          LSPA(I,J) = LSPA(I,J) + PREC(I,J)
          GO TO 100
        ELSE
          IF (IHR2 .LE. PCPHR) then
            DDATA(I,J) = PPTDAT(I,J,IHR1)*FRACT1                        &
     &                 + PPTDAT(I,J,IHR2)*FRACT2
          ELSE
            DDATA(I,J) = PPTDAT(I,J,IHR1)*FRACT1 
          ENDIF
!
           LSPA(I,J) = LSPA(I,J) + DDATA(I,J)
        ENDIF
        IF (I.EQ.1 .AND. J.EQ.1 .AND. MYPE.EQ.0) THEN
          WRITE(0,1020) DDATA(I,J), PREC(I,J), LSPA(I,J)
 1020     FORMAT('ADJPPT: DDATA=',E12.6, ' PREC=',E12.6,' LSPA=',E12.6)
        ENDIF

 100    CONTINUE
      ENDDO
      ENDDO
!
 200  CONTINUE

      END SUBROUTINE ADJPPT
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      END MODULE MODULE_PRECIP_ADJUST
!
!-----------------------------------------------------------------------
