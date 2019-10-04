!-----------------------------------------------------------------------
!
      MODULE MODULE_MICROPHYSICS_NMM
!
!-----------------------------------------------------------------------
!
!***  THE MICROPHYSICS DRIVERS AND PACKAGES

!      11-06-2009 W. Wang put NAM micorphysics into a single module   
!      02-10-2010 W. Wang added wsm6 
!-----------------------------------------------------------------------
!
! HISTORY LOG:
!
!    11-06-2009 W. Wang - Put NAM/Ferrier microphysics into 
!                         a single module.  
!
!-----------------------------------------------------------------------
!
      USE MODULE_KINDS
!
      USE MODULE_CONSTANTS,ONLY : CICE,CLIQ,CPV,EP_1,EP_2,EPSILON,G     &
                                 ,P608,PSAT,R_D,R_V,RHOAIR0,RHOWATER    &
                                 ,SVPT0,XLF,XLV                         &
                                 ,CAPPA,CP,EPSQ
!
      USE MODULE_CONTROL,ONLY : NMMB_FINALIZE
! MP options
      USE MODULE_MP_ETANEW
      USE MODULE_MP_FER_HIRES
      USE MODULE_MP_WSM6
      USE MODULE_MP_THOMPSON
      USE MODULE_MP_GFS
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      PRIVATE
!
!      PUBLIC :: FERRIER_INIT,FPVS,GSMDRIVE,WSM3INIT
      PUBLIC :: GSMDRIVE
!
!-----------------------------------------------------------------------
!
      INTEGER :: MYPE
      REAL, PRIVATE,PARAMETER ::                                     &
!--- Physical constants follow:
           XLS=2.834E6,R_G=1./G
!
      INTEGER,PUBLIC,PARAMETER :: MICRO_RESTART=7501
!
      CONTAINS
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
      SUBROUTINE GSMDRIVE(ITIMESTEP,DT,NPHS                             &
                         ,DX,DY,SM,FIS                                  &
                         ,DSG2,SGML2,PDSG1,PSGML1,PT,PD                 &
                         ,T,Q,CWM,OMGALF                                &
                         ,TRAIN,SR                                      &
                         ,F_ICE,F_RAIN,F_RIMEF                          &
                         ,QC,QR,QI,QS,QG,NI,NR                          &
                         ,F_QC,F_QR,F_QI,F_QS,F_QG,F_NI,F_NR            &
                         ,PREC,ACPREC,AVRAIN,ACPREC_TOT                 &
                         ,acpcp_ra,acpcp_sn,acpcp_gr, refl_10cm         &
                         ,re_cloud,re_ice,re_snow                       &
                         ,has_reqc,has_reqi,has_reqs                    &
                         ,MP_RESTART_STATE                              &
                         ,TBPVS_STATE,TBPVS0_STATE                      &
                         ,SPECIFIED,NESTED                              &
                         ,MICROPHYSICS                                  &
                         ,RHGRD                                         &
                         ,TP1,QP1,PSP1                                  &
                         ,USE_RADAR                                     &
                         ,DFI_TTEN                                      &
                         ,IDS,IDE,JDS,JDE,LM                            &
                         ,IMS,IME,JMS,JME                               &
                         ,ITS,ITE,JTS,JTE                               &
                         ,ITS_B1,ITE_B1,JTS_B1,JTE_B1,MPRATES,D_SS)
!***********************************************************************
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .
! SUBPROGRAM:    GSMDRIVE    MICROPHYSICS OUTER DRIVER
!   PRGRMMR: BLACK           ORG: W/NP22     DATE: 02-03-26
!
! ABSTRACT:
!     RADIATION SERVES AS THE INTERFACE BETWEEN THE NMMB PHYSICS COMPONENT
!     AND THE WRF MICROPHYSICS DRIVER.
!
! PROGRAM HISTORY LOG:
!   02-03-26  BLACK      - ORIGINATOR
!   04-11-18  BLACK      - THREADED
!   06-07-31  BLACK      - BUILT INTO NMMB PHYSICS COMPONENT
!   08-08     JANJIC     - Synchronize WATER array and Q.
!
! USAGE: CALL GSMDRIVE FROM PHY_RUN
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE : IBM
!$$$
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      INTEGER,INTENT(IN) :: D_SS,ITIMESTEP,NPHS                         &
                           ,IDS,IDE,JDS,JDE,LM                          &
                           ,IMS,IME,JMS,JME                             &
                           ,ITS,ITE,JTS,JTE                             &
                           ,ITS_B1,ITE_B1,JTS_B1,JTE_B1
!
      LOGICAL,INTENT(IN) :: USE_RADAR
!
      REAL,INTENT(IN) :: DT,DX,DY,PT,RHGRD
!
      REAL,DIMENSION(IMS:IME,JMS:JME,1:LM,D_SS)  :: MPRATES
!
      REAL,DIMENSION(1:LM),INTENT(IN) :: DSG2,PDSG1,PSGML1,SGML2
!
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: FIS,PD,SM
!
      REAL,DIMENSION(IMS:IME,JMS:JME,1:LM),INTENT(IN) :: OMGALF         &
                                                        ,DFI_TTEN
!
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT) :: ACPREC,PREC      &
                                                      ,ACPREC_TOT       &
                                                      ,AVRAIN              !<-- Was a scalar
! G. Thompson added next 4 lines.
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT) :: acpcp_ra,acpcp_sn,acpcp_gr
      REAL,DIMENSION(IMS:IME,JMS:JME,1:LM),INTENT(INOUT) :: refl_10cm
      REAL,DIMENSION(IMS:IME,JMS:JME,1:LM),INTENT(INOUT) :: re_cloud, re_ice, re_snow
      INTEGER,INTENT(IN):: has_reqc, has_reqi, has_reqs
!
      REAL,DIMENSION(IMS:IME,JMS:JME,1:LM),INTENT(INOUT) :: CWM,Q,T     &
                                                           ,TRAIN
!
      REAL(kind=KFPT),DIMENSION(:,:,:),POINTER,INTENT(INOUT) ::         &
     &                                          QC,QI,QR,QS,QG,NI,NR
!
      REAL,DIMENSION(IMS:IME,JMS:JME,1:LM),INTENT(INOUT) :: F_ICE       &
                                                           ,F_RAIN      &
                                                           ,F_RIMEF
!
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(OUT) :: SR
!
      CHARACTER(99),INTENT(IN) :: MICROPHYSICS
!
      LOGICAL,INTENT(IN) :: NESTED,SPECIFIED
!
      LOGICAL,INTENT(IN) :: F_QC,F_QR,F_QI,F_QS,F_QG,F_NI,F_NR
!
!***  State Variables for ETAMPNEW Microphysics 
!
      REAL,DIMENSION(:),INTENT(INOUT) :: MP_RESTART_STATE               &
                                        ,TBPVS_STATE,TBPVS0_STATE
!*** GFS microphysics
      REAL, DIMENSION(IMS:IME,JMS:JME,1:LM), INTENT(INOUT) :: TP1,QP1
      REAL, DIMENSION(IMS:IME,JMS:JME), INTENT(INOUT)      :: PSP1
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER :: I,IJ,J,K,MP_PHYSICS,N,NTSD
      INTEGER :: ITSLOC,ITELOC,JTSLOC,JTELOC 
!
      INTEGER,DIMENSION(IMS:IME,JMS:JME) :: LOWLYR
!
      REAL :: DPL,DTPHS,PCPCOL,PDSL,PHMID,QW,RDTPHS,TNEW
      REAL :: MP_TTEN,mytten
!
      REAL,DIMENSION(1:LM) :: QL,TL
!
      REAL,DIMENSION(IMS:IME,JMS:JME) :: CUBOT,CUTOP,RAINNC,RAINNCV     &
                                        ,SNOWNC,SNOWNCV,XLAND           &
                                        ,graupelnc,graupelncv
!
      REAL,DIMENSION(IMS:IME,JMS:JME,1:LM) :: DZ                        &
                                             ,P_PHY,PI_PHY              &
                                             ,RR,TH_PHY,QV
!
      LOGICAL :: WARM_RAIN,F_QT,USE_QV
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      NTSD=ITIMESTEP
      DTPHS=NPHS*DT
      RDTPHS=1./DTPHS
!
!-- AVRAIN was a scalar but changed to a 2D array to allow for updates in ESMF
!
      DO J=JTS,JTE
      DO I=ITS,ITE
         AVRAIN(I,J)=AVRAIN(I,J)+1.
      ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
!***  NOTE:  THE NMMB HAS IJK STORAGE WITH LAYER 1 AT THE TOP.
!***         THE WRF PHYSICS DRIVERS HAVE IKJ STORAGE WITH LAYER 1
!***         AT THE BOTTOM.
!-----------------------------------------------------------------------
!
!!!   DO J=JTS_B1,JTE_B1
!!!   DO I=ITS_B1,ITE_B1
!.......................................................................
!$omp parallel do                                                       &
!$omp& private(j,i,k,pdsl,dpl,phmid,ql,tl)
!.......................................................................
      DO J=JTS,JTE
      DO I=ITS,ITE
!
        PDSL=PD(I,J)
        LOWLYR(I,J)=1
        XLAND(I,J)=SM(I,J)+1.
!
!-----------------------------------------------------------------------
!***   FILL RAINNC WITH ZERO (NORMALLY CONTAINS THE NONCONVECTIVE
!***                          ACCUMULATED RAIN BUT NOT YET USED BY NMM)
!***   COULD BE OBTAINED FROM ACPREC AND CUPREC (ACPREC-CUPREC)
!-----------------------------------------------------------------------
!..The NC variables were designed to hold simulation total accumulations
!.. whereas the NCV variables hold timestep only values, so change below
!.. to zero out only the timestep amount preparing to go into each
!.. micro routine while allowing NC vars to accumulate continually.
!.. But, the fact is, the total accum variables are local, never saved
!.. nor written so they go nowhere at the moment.
!
        RAINNC (I,J)=0. ! NOT YET USED BY NMM
        RAINNCv(I,J)=0.
        SNOWNCv(I,J)=0.
        graupelncv(i,j) = 0.0
!
!-----------------------------------------------------------------------
!***  FILL THE SINGLE-COLUMN INPUT
!-----------------------------------------------------------------------
!
        DO K=LM,1,-1   ! We are moving down from the top in the flipped arrays
!
          DPL=DSG2(K)*PDSL+PDSG1(K)
          PHMID=SGML2(K)*PDSL+PSGML1(K)
          TL(K)=T(I,J,K)
          QL(K)=AMAX1(Q(I,J,K),EPSQ)
!
          RR(I,J,K)=PHMID/(R_D*TL(K)*(P608*QL(K)+1.))
          PI_PHY(I,J,K)=(PHMID*1.E-5)**CAPPA
          TH_PHY(I,J,K)=TL(K)/PI_PHY(I,J,K)
          P_PHY(I,J,K)=PHMID
          DZ(I,J,K)=DPL*R_G/RR(I,J,K)
!
        ENDDO    !- DO K=LM,1,-1
!
      ENDDO    !- DO I=ITS,ITE
      ENDDO    !- DO J=JTS,JTE
!.......................................................................
!$omp end parallel do
!.......................................................................
!
!-----------------------------------------------------------------------
!***  IF NEEDED, UPDATE WATER VAPOR RATIO FROM SPECIFIC HUMIDITY.
!-----------------------------------------------------------------------
!
      IF(TRIM(MICROPHYSICS)=='wsm6' .OR. TRIM(MICROPHYSICS)=='thompson')THEN
        USE_QV=.TRUE.    !-- Initialize QV, update Q & CWM at the end
      ELSE    
        USE_QV=.FALSE.
      ENDIF   
!
      IF(USE_QV) THEN
!.......................................................................
!$omp parallel do                                                       &
!$omp& private(i,j,k)
!.......................................................................
        DO K=1,LM
          DO J=JTS_B1,JTE_B1
          DO I=ITS_B1,ITE_B1
            QV(I,J,K)=Q(I,J,K)/(1.-Q(I,J,K))
          ENDDO
          ENDDO
        ENDDO
!.......................................................................
!$omp end parallel do
!.......................................................................
      ENDIF
!
!-----------------------------------------------------------------------
!
!***  CALL MICROPHYSICS
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!***  TRANSLATE THE MICROPHYSICS OPTIONS IN THE CONFIG FILE TO THEIR
!***  ANALOGS IN THE WRF REGISTRY SO THAT THE WRF MICROPHYSICS DRIVER
!***  REMAINS UNTOUCHED.
!-----------------------------------------------------------------------
!
!     SELECT CASE (TRIM(MICROPHYSICS))
!       CASE ('fer')
!          MP_PHYSICS=5
!       CASE ('kes')
!         MP_PHYSICS=1
!       CASE ('lin')
!         MP_PHYSICS=2
!       CASE ('wsm3')
!         MP_PHYSICS=3
!       CASE ('tho')
!         MP_PHYSICS=8
!       CASE DEFAULT
!         WRITE(0,*)' User selected MICROPHYSICS=',MICROPHYSICS
!         WRITE(0,*)' Improper selection of Microphysics scheme in GSMDRIVE'
!!!       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
!         CALL NMMB_FINALIZE
!     END SELECT
!
!---------------------------------------------------------------------
!  Check for microphysics type.  We need a clean way to
!  specify these things!
!---------------------------------------------------------------------


        ITSLOC = MAX(ITS_B1,IDS)
        ITELOC = MIN(ITE_B1,IDE-1)
        JTSLOC = MAX(JTS_B1,JDS)
        JTELOC = MIN(JTE_B1,JDE-1)

        micro_select: SELECT CASE (TRIM(MICROPHYSICS))
!
          CASE ('fer')
            CALL ETAMP_NEW(                                                   &
                   ITIMESTEP=ntsd,DT=dtphs,DX=dx,DY=dy                        &
                  ,DZ8W=dz,RHO_PHY=rr,P_PHY=p_phy,PI_PHY=pi_phy,TH_PHY=th_phy &
                  ,Q=Q                                                        &
                  ,QC=QC                                                      &
                  ,QS=QS                                                      &
                  ,QR=QR                                                      &
                  ,QT=cwm                                                     &
                  ,LOWLYR=LOWLYR,SR=SR                                        &
                  ,F_ICE_PHY=F_ICE,F_RAIN_PHY=F_RAIN                          &
                  ,F_RIMEF_PHY=F_RIMEF                                        &
                  ,RAINNC=rainnc,RAINNCV=rainncv                              &
                  ,IDS=ids,IDE=ide, JDS=jds,JDE=jde, KDS=1,KDE=LM+1           &
                  ,IMS=ims,IME=ime, JMS=jms,JME=jme, KMS=1,KME=LM             &
                  ,ITS=itsloc,ITE=iteloc, JTS=jtsloc,JTE=jteloc, KTS=1,KTE=LM &
                  ,MP_RESTART_STATE=mp_restart_state                          &
                  ,TBPVS_STATE=tbpvs_state,TBPVS0_STATE=tbpvs0_state          &
                  ,D_SS=d_ss,MPRATES=mprates                                  &
                                                                            )
          CASE ('fer_hires')
            CALL FER_HIRES(                                                   &
                   ITIMESTEP=ntsd,DT=dtphs,DX=dx,DY=dy,RHgrd=RHGRD            &
                  ,DZ8W=dz,RHO_PHY=rr,P_PHY=p_phy,PI_PHY=pi_phy,TH_PHY=th_phy &
                  ,Q=Q                                                        &
                  ,QC=QC                                                      &
                  ,QS=QS                                                      &
                  ,QR=QR                                                      &
                  ,QT=cwm                                                     &
                  ,LOWLYR=LOWLYR,SR=SR                                        &
                  ,F_ICE_PHY=F_ICE,F_RAIN_PHY=F_RAIN                          &
                  ,F_RIMEF_PHY=F_RIMEF                                        &
                  ,RAINNC=rainnc,RAINNCV=rainncv                              &
                  ,IDS=ids,IDE=ide, JDS=jds,JDE=jde, KDS=1,KDE=LM+1           &
                  ,IMS=ims,IME=ime, JMS=jms,JME=jme, KMS=1,KME=LM             &
                  ,ITS=itsloc,ITE=iteloc, JTS=jtsloc,JTE=jteloc, KTS=1,KTE=LM &
                  ,MP_RESTART_STATE=mp_restart_state                          &
                  ,TBPVS_STATE=tbpvs_state,TBPVS0_STATE=tbpvs0_state          &
                  ,D_SS=d_ss,MPRATES=mprates                                  &
                  ,refl_10cm=refl_10cm                                        &
                                                                            )
          CASE ('gfs')
            CALL GFSMP(DT=dtphs,                                               &
                   dz8w=dz,rho_phy=rr,p_phy=p_phy,pi_phy=pi_phy,th_phy=th_phy, &
                   SR=SR,QT=CWM, F_ICE_PHY=F_ICE,                              &
                   RAINNC=RAINNC,RAINNCV=RAINNCV,                              &
                   Q=Q,QC=QC,QI=QI,                                            &
                   F_QC=F_QC,F_QI=F_QI,                                        &
                   TP1=TP1,QP1=QP1,PSP1=PSP1,                                  &
                   IDS=ids,IDE=ide, JDS=jds,JDE=jde, KDS=1,KDE=LM+1,           &
                   IMS=ims,IME=ime, JMS=jms,JME=jme, KMS=1,KME=LM ,            &
                   ITS=itsloc,ITE=iteloc, JTS=jtsloc,JTE=jteloc, KTS=1,KTE=LM )
          CASE ('wsm6')
             CALL wsm6(                                             &
                  TH=th_phy                                         &
                 ,Q=QV                                              &
                 ,QC=QC                                             &
                 ,QR=QR                                             &
                 ,QI=QI                                             &
                 ,QS=QS                                             &
                 ,QG=QG                                             &
                 ,DEN=rr,PII=pi_phy,P=p_phy,DELZ=dz                 &
                 ,DELT=dtphs,G=g,CPD=cp,CPV=cpv                     &
                 ,RD=r_d,RV=r_v,T0C=svpt0                           &
                 ,EP1=ep_1, EP2=ep_2, QMIN=epsilon                  &
                 ,XLS=xls, XLV0=xlv, XLF0=xlf                       &
                 ,DEN0=rhoair0, DENR=rhowater                       &
                 ,CLIQ=cliq,CICE=cice,PSAT=psat                     &
                 ,RAIN=rainnc ,RAINNCV=rainncv                      &
                 ,SNOW=snownc ,SNOWNCV=snowncv                      &
                 ,SR=sr                                             &
                 ,GRAUPEL=graupelnc ,GRAUPELNCV=graupelncv          &
                 ,IDS=ids,IDE=ide, JDS=jds,JDE=jde, KDS=1,KDE=LM+1  &
                 ,IMS=ims,IME=ime, JMS=jms,JME=jme, KMS=1,KME=LM    &
                 ,ITS=itsloc,ITE=iteloc, JTS=jtsloc,JTE=jteloc, KTS=1,KTE=LM &
                 ,D_SS=d_ss,MPRATES=mprates)
          CASE ('thompson')
!+---+-----------------------------------------------------------------+
!            write(6,*)'DEBUG-GT, calling mp_gt_driver'
             CALL mp_gt_driver(                                     &
                  qv=qv                                             &
                 ,qc=qc                                             &
                 ,qr=qr                                             &
                 ,qi=qi                                             &
                 ,qs=qs                                             &
                 ,qg=qg                                             &
                 ,ni=ni                                             &
                 ,nr=nr                                             &
                 ,TH=th_phy,PII=pi_phy,P=p_phy,dz=dz,dt_in=dtphs    &
                 ,itimestep=ntsd                                    &
                 ,RAINNC=rainnc ,RAINNCV=rainncv                    &
                 ,SNOWNC=snownc ,SNOWNCV=snowncv                    &
                 ,GRAUPELNC=graupelnc ,GRAUPELNCV=graupelncv        &
                 ,SR=sr                                             &
                 ,refl_10cm=refl_10cm(ims,jms,1)                    &
                 ,diagflag=.true.                                   &
                 ,do_radar_ref=1                                    &
                 ,re_cloud=re_cloud(ims,jms,1)                      &
                 ,re_ice=re_ice(ims,jms,1)                          &
                 ,re_snow=re_snow(ims,jms,1)                        &
                 ,has_reqc=has_reqc                                 &
                 ,has_reqi=has_reqi                                 &
                 ,has_reqs=has_reqs                                 &
                 ,IDS=ids,IDE=ide, JDS=jds,JDE=jde, KDS=1,KDE=LM+1  &
                 ,IMS=ims,IME=ime, JMS=jms,JME=jme, KMS=1,KME=LM    &
                 ,ITS=itsloc,ITE=iteloc, JTS=jtsloc,JTE=jteloc, KTS=1,KTE=LM &
                 ,D_SS=d_ss,MPRATES=mprates                         )
!
!..rainncv is actually all precip, so need to subtract snow/graupel to isolate rain only
!
                DO J=JMS,JME
                DO I=IMS,IME
                   acpcp_sn(I,J) = acpcp_sn(I,J) + snowncv(i,j)
                   acpcp_gr(I,J) = acpcp_gr(I,J) + graupelncv(i,j)
                   acpcp_ra(I,J) = acpcp_ra(I,J)                        &
                    + MAX(0., rainncv(i,j)-snowncv(i,j)-graupelncv(i,j))
                ENDDO
                ENDDO
!+---+-----------------------------------------------------------------+

          CASE DEFAULT
            WRITE(0,*)' The microphysics option does not exist: MICROPHYSICS = ',TRIM(MICROPHYSICS)
            CALL NMMB_FINALIZE

        END SELECT micro_select

!            
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  THE FOLLOWING MUST BE RECONCILED WHEN THREADING IS TURNED ON.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!jaa!$omp parallel do                                                       &
!jaa!$omp& private(ij)
!     DO IJ=1,NUM_TILES
!       CALL MICROPHYSICS_ZERO_OUT(                                     &
!                    WATER,N_MOIST,CONFIG_FLAGS                         &
!                   ,IDS,IDE,JDS,JDE,KDS,KDE                            &
!                   ,IMS,IME,JMS,JME,KMS,KME                            &
!                   ,GRID%I_START(IJ),GRID%I_END(IJ)                    &
!                   ,GRID%J_START(IJ),GRID%J_END(IJ)                    &
!                   ,KTS,KTE                                       )
!     ENDDO
!
!-----------------------------------------------------------------------
!
      IF(USE_QV) THEN    !-- Update Q & CWM for WSM6 & Thompson microphysics
!.......................................................................
!$omp parallel do                                                       &
!$omp& private(i,j,k)
!.......................................................................
        DO K=1,LM
          DO J=JTS_B1,JTE_B1
          DO I=ITS_B1,ITE_B1
            Q(I,J,K)=QV(I,J,K)/(1.+QV(I,J,K))
            CWM(I,J,K)=QC(i,j,k)+QR(i,j,k)+QI(i,j,k)+QS(i,j,k)+QG(i,j,k)
          ENDDO
          ENDDO
        ENDDO
!.......................................................................
!$omp end parallel do
!.......................................................................
      ENDIF
!
!.......................................................................
!$omp parallel do                                                       &
!$omp& private(i,j,k,TNEW,MP_TTEN)
!.......................................................................
      DO K=1,LM
        DO J=JTS_B1,JTE_B1
        DO I=ITS_B1,ITE_B1
!
!-----------------------------------------------------------------------
!***  UPDATE TEMPERATURE, SPECIFIC HUMIDITY, CLOUD WATER, AND HEATING.
!-----------------------------------------------------------------------
!
          TNEW=TH_PHY(I,J,K)*PI_PHY(I,J,K)
          TRAIN(I,J,K)=TRAIN(I,J,K)+(TNEW-T(I,J,K))*RDTPHS
          IF (USE_RADAR) THEN
            MP_TTEN=(TNEW-T(I,J,K))*RDTPHS              
            IF(DFI_TTEN(I,J,K)>MP_TTEN.AND.DFI_TTEN(I,J,K)<0.01        &
                                      .AND.MP_TTEN<0.0018)THEN 
              MP_TTEN=DFI_TTEN(I,J,K)                              
            END IF                                                
            T(I,J,K)=T(I,J,K)+MP_TTEN/RDTPHS
          ELSE
            T(I,J,K)=TNEW
          ENDIF
        ENDDO
        ENDDO
      ENDDO
!.......................................................................
!$omp end parallel do
!.......................................................................
!
!-----------------------------------------------------------------------
!***  UPDATE PRECIPITATION
!-----------------------------------------------------------------------
!
!jaa!$omp parallel do                                                       &
!jaa!$omp& private(i,j,pcpcol)
      DO J=JTS_B1,JTE_B1
      DO I=ITS_B1,ITE_B1
        PCPCOL=RAINNCV(I,J)*1.E-3
        PREC(I,J)=PREC(I,J)+PCPCOL
        ACPREC(I,J)=ACPREC(I,J)+PCPCOL
        ACPREC_TOT(I,J)=ACPREC_TOT(I,J)+PCPCOL
!
! NOTE: RAINNC IS ACCUMULATED INSIDE MICROPHYSICS BUT NMM ZEROES IT OUT ABOVE
!       SINCE IT IS ONLY A LOCAL ARRAY FOR NOW
!
      ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE GSMDRIVE
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------

      END MODULE MODULE_MICROPHYSICS_NMM

!-----------------------------------------------------------------------
