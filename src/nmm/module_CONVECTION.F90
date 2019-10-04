!-----------------------------------------------------------------------
!
      MODULE MODULE_CONVECTION
!
!-----------------------------------------------------------------------
!
!***  THE CONVECTION DRIVERS AND PACKAGES
!
!-----------------------------------------------------------------------
!
      USE MODULE_KINDS
!
      USE MODULE_CONTROL,ONLY : NMMB_FINALIZE

      USE MODULE_CU_BMJ
      USE MODULE_CU_SAS
      USE MODULE_CU_SASHUR
      USE MODULE_CU_SCALE
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      PRIVATE
!
      PUBLIC :: CUCNVC
      PUBLIC :: BMJSCHEME, SASSCHEME, SASHURSCHEME, SCALECUSCHEME
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  THE CONVECTION OPTIONS
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      INTEGER(KIND=KINT),PARAMETER :: KFETASCHEME=1                     &
                                     ,BMJSCHEME=2                       &
                                     ,GDSCHEME=3                        &
                                     ,SASSCHEME=4                       &
                                     ,SASHURSCHEME=84                   &
                                     ,SCALECUSCHEME=94
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
      SUBROUTINE CUCNVC(NTSD,DT,NCNVC,NRADS,NRADL,MINUTES_HISTORY       &
                       ,ENTRAIN,NEWALL,NEWSWAP,NEWUPUP,NODEEP           &
                       ,FRES,FR,FSL,FSS                                 &
                       ,DYH,RESTRT,HYDRO                                &
                       ,CLDEFI                                          &
                       ,F_ICE,F_RAIN                                    &
                       ,QC,QR,QI,QS,QG                                  &
                       ,F_QC,F_QR,F_QI,F_QS,F_QG                        &
                       ,DSG2,SGML2,SG2,PDSG1,PSGML1,PSG1                &
                       ,dxh                                             &
                       ,PT,PD,T,Q,CWM,TCUCN                             &
                       ,OMGALF,U,V                                      &
                       ,FIS,W0AVG                                       &
                       ,PREC,ACPREC,CUPREC,ACPREC_TOT,CUPPT,CPRATE      &
                       ,CNVBOT,CNVTOP,SM,LPBL                           &
                       ,HTOP,HTOPD,HTOPS                                &
                       ,HBOT,HBOTD,HBOTS                                &
                       ,AVCNVC,ACUTIM                                   &
                       ,RSWIN,RSWOUT                                    &
                       ,CONVECTION,CU_PHYSICS,MICROPHYSICS              &
!!!! added for SAS
                       ,SICE,QWBS,TWBS,PBLH,DUDT_PHY,DVDT_PHY           &
!!!
!!!  added for SAS-hurricane
                      ,MOMMIX,PGCON,SAS_MASS_FLUX   &   ! hwrf,namelist
                       ,SHALCONV,SHAL_PGCON          &   !hwrf,namelist
                       ,W_TOT,PSGDT                  &   ! test w from omgalf vs W_tot
!!
                       ,A2,A3,A4,CAPPA,CP,ELIV,ELWV,EPSQ,G              &
                       ,P608,PQ0,R_D,TIW                                &
                       ,IDS,IDE,JDS,JDE,LM                              &
                       ,IMS,IME,JMS,JME                                 &
                       ,ITS,ITE,JTS,JTE                                 &
                       ,ITS_B1,ITE_B1,JTS_B1,JTE_B1                     &
                                                    )
!***********************************************************************
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .
! SUBPROGRAM:    CUCNVC      CONVECTIVE PRECIPITATION OUTER DRIVER
!   PRGRMMR: BLACK           ORG: W/NP22     DATE: 02-03-21
!
! ABSTRACT:
!     CUCVNC DRIVES THE WRF CONVECTION SCHEMES
!
! PROGRAM HISTORY LOG:
!   02-03-21  BLACK      - ORIGINATOR
!   04-11-18  BLACK      - THREADED
!   06-10-11  BLACK      - BUILT INTO UMO PHYSICS COMPONENT
!   08-08     JANJIC     - Synchronize WATER array and Q.
!   10-10-26  WEIGUO WANG - add GFS SAS convection
!   14-06-19  WEIGUO WANG - add hurricane SAS (moved from hwrf)
!   16-08-29  WEIGUO WANG - add scale-aware convection schemes
! USAGE: CALL CUCNVC FROM PHY_RUN
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
      character(99),intent(in):: &
       convection,microphysics
!
      logical(kind=klog),intent(in):: &
       hydro,restrt &
      ,entrain,newall,newswap,newupup,nodeep &
      ,f_qc,f_qr,f_qi,f_qs,f_qg
!
      integer(kind=kint),intent(in):: &
       cu_physics &
      ,ids,ide,jds,jde,lm &
      ,ims,ime,jms,jme &
      ,its,ite,jts,jte &
      ,its_b1,ite_b1,jts_b1,jte_b1 &
      ,ncnvc,minutes_history &
      ,nrads,nradl,ntsd
!
      integer(kind=kint),dimension(ims:ime,jms:jme),intent(in):: &
       lpbl
!
      real(kind=kfpt),intent(in):: &
       a2,a3,a4,cappa,cp,dt,dyh,eliv,elwv,epsq &
      ,fres,fr,fsl,fss,g,p608,pq0,pt,r_d,tiw
!
      real(kind=kfpt),dimension(1:lm),intent(in):: &
       dsg2,pdsg1,psgml1,sgml2
!
      real(kind=kfpt),dimension(1:lm+1),intent(in):: &
       psg1,sg2
!
      real(kind=kfpt),dimension(jds:jde),intent(in):: &
       dxh
!
      real(kind=kfpt),dimension(ims:ime,jms:jme),intent(in):: &
       fis,pd &
      ,rswin,rswout,sm
!
      real(kind=kfpt),dimension(ims:ime,jms:jme),intent(inout):: &
       acprec,cldefi &
      ,acprec_tot &
      ,cnvbot,cnvtop &
      ,cuppt,cuprec &
      ,hbot,htop &
      ,hbotd,htopd &
      ,hbots,htops &
      ,prec,cprate &
      ,acutim,avcnvc  !<-- were scalars
!
      real(kind=kfpt),dimension(ims:ime,jms:jme),intent(in):: &
       sice,qwbs,twbs,pblh  !fOR SAS

      REAL(kind=kfpt), OPTIONAL, INTENT(IN) :: &
              PGCON,sas_mass_flux,shal_pgcon,mommix,shalconv         !sashur
!!      INTEGER(kind=kint), OPTIONAL, INTENT(IN) :: shalconv  !sashur
      real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(in)::W_TOT
      real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm-1),intent(in)::PSGDT !vertical mass flux
!
      real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(in):: &
       omgalf,u,v
!
      real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(out):: &
       dudt_phy,dvdt_phy
!
      real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(inout):: &
       q,t &
      ,f_ice &
      ,f_rain &
      ,cwm &
      ,tcucn
!
      real(kind=kfpt),dimension(ims:ime,1:lm+1,jms:jme),intent(inout):: &
       w0avg
!
      REAL(KIND=KFPT),DIMENSION(:,:,:),POINTER,INTENT(INOUT) ::         &
     &                                                QC,QI,QR,QS,QG
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      logical(kind=klog):: &
       restart,warm_rain,F_QGr
!
      logical(kind=klog),dimension(ims:ime,jms:jme):: &
       cu_act_flag
!
      integer(kind=kint):: &
       i,j &
      ,k &
      ,mnto &
      ,n,ncubot,ncutop,n_timstps_output
!
      integer(kind=kint),dimension(ims:ime,jms:jme):: &
       KPBL,LBOT,LTOP
!
      real(kind=kfpt):: &
       cf_hi,dtcnvc,dtdt,fice,frain,g_inv &
      ,pcpcol,pdsl,ql,ql_k,rdtcnvc &
      ,QCW,QCI,QRain,QSnow,QGraup &
      ,tl
!
      REAL(kind=kfpt),DIMENSION(IMS:IME,JMS:JME):: &
       CUBOT,CUTOP,NCA &
      ,RAINC,RAINCV,SFCZ,XLAND
!
      REAL(kind=kfpt),DIMENSION(IMS:IME,JMS:JME,1:LM):: &
       DZ,PHMID,exner &
      ,th,rr &
      ,RQCCUTEN,RQRCUTEN &
      ,RQICUTEN,RQSCUTEN &
      ,RQCUTEN,RTHCUTEN &
      ,RQGCUTEN &
      ,u_phy,v_phy

      REAL(kind=kfpt),DIMENSION(IMS:IME,JMS:JME,1:LM+1):: &
       PHINT
!-----------------------------------------------------------------------
!***  For temperature change check only.
!-----------------------------------------------------------------------
!zj      REAL(kind=kfpt) :: DTEMP_CHECK=1.0
      REAL(kind=kfpt) :: TCHANGE
!-----------------------------------------------------------------------
!***********************************************************************
!
!-----------------------------------------------------------------------
!***  RESET THE HBOT/HTOP CONVECTIVE CLOUD BOTTOM (BASE) AND TOP ARRAYS
!***  USED IN RADIATION.  THEY STORE THE MAXIMUM VERTICAL LIMITS OF
!***  CONVECTIVE CLOUD BETWEEN RADIATION CALLS.  THESE ARRAYS ARE OUT
!***  OF THE WRF PHYSICS AND THUS THEIR VALUES INCREASE UPWARD.
!***  CUPPT IS THE ACCUMULATED CONVECTIVE PRECIPITATION BETWEEN
!***  RADIATION CALLS.
!-----------------------------------------------------------------------
!
      IF(MOD(NTSD,NRADS)==0.OR.MOD(NTSD,NRADL)==0)THEN
         DO J=JMS,JME
         DO I=IMS,IME
           HTOP(I,J)=0.
           HBOT(I,J)=REAL(LM+1)
           CUPPT(I,J)=0.
         ENDDO
         ENDDO
      ENDIF
!
!-----------------------------------------------------------------------
      IF(MOD(NTSD,NCNVC)/=0.AND.CONVECTION=='bmj')RETURN
      IF(MOD(NTSD,NCNVC)/=0.AND.CONVECTION=='sas')RETURN
      IF(MOD(NTSD,NCNVC)/=0.AND.CONVECTION=='sashur')RETURN
      IF(MOD(NTSD,NCNVC)/=0.AND.CONVECTION=='scalecu')RETURN
!-----------------------------------------------------------------------
!
      RESTART=RESTRT
!
!-----------------------------------------------------------------------
!
      IF(MICROPHYSICS=='fer' .OR. MICROPHYSICS=='fer_hires') THEN
         F_QGr=.FALSE.
      ELSE
         F_QGr=F_QG
      ENDIF
      IF(CONVECTION=='kf')THEN
!
        IF(.NOT.RESTART.AND.NTSD==0)THEN
!jaa!zj$omp parallel do                                                       &
!jaa!zj$omp& private(i,j,k)
          DO J=JTS,JTE
          DO K=1,LM+1
          DO I=ITS,ITE
            W0AVG(I,K,J)=0.
          ENDDO
          ENDDO
          ENDDO
        ENDIF
!
      ENDIF
!
!-----------------------------------------------------------------------
!***  GENERAL PREPARATION
!-----------------------------------------------------------------------
!
!-- AVCNVC,ACUTIM were scalars but changed to 2D arrays to allow for updates in ESMF
!
      DO J=JTS,JTE
      DO I=ITS,ITE
         AVCNVC(I,J)=AVCNVC(I,J)+1.
         ACUTIM(I,J)=ACUTIM(I,J)+1.
      ENDDO
      ENDDO
!
      DTCNVC=NCNVC*DT
      RDTCNVC=1./DTCNVC
      G_INV=1./G
!
!.......................................................................
!zj$omp parallel do &
!zj$omp& private(j,i,k,pdsl,ql,tl)
!.......................................................................
      DO J=JTS,JTE
      DO I=ITS,ITE
!
        PDSL=PD(I,J)
        RAINCV(I,J)=0.
        RAINC(I,J)=0.
        PHINT(I,J,LM+1)=SG2(LM+1)*PDSL+PSG1(LM+1)
        XLAND(I,J)=SM(I,J)+1.
        NCA(I,J)=0.
        SFCZ(I,J)=FIS(I,J)*G_INV
!
        CUTOP(I,J)=999.
        CUBOT(I,J)=999.
!
!***  LPBL IS THE MODEL LAYER CONTAINING THE PBL TOP
!***  COUNTING DOWNWARD FROM THE TOP OF THE DOMAIN
!***  SO KPBL IS THE SAME LAYER COUNTING UPWARD FROM
!***  THE GROUND.
!
        KPBL(I,J)=LPBL(I,J)
!
!-----------------------------------------------------------------------
!***  FILL VERTICAL WORKING ARRAYS.
!-----------------------------------------------------------------------
!
        DO K=1,LM
!
          PHINT(I,J,K)=SG2(K)*PDSL+PSG1(K) !zj
          PHMID(I,J,K)=SGML2(K)*PDSL+PSGML1(K)

          QL=MAX(Q(I,J,K),EPSQ)
          TL=T(I,J,K)
          RR(I,J,K)=PHMID(I,J,K)/(R_D*TL*(.608*ql+1.))
          T(I,J,K)=TL
!
          EXNER(I,J,K)=(PHMID(I,J,K)*1.E-5)**CAPPA
          TH(I,J,K)=TL/EXNER(I,J,K)
!
        ENDDO
      ENDDO
      ENDDO
!.......................................................................
!zj$omp end parallel do
!.......................................................................
!
!-----------------------------------------------------------------------
!***  Compute velocity components at mass points.
!-----------------------------------------------------------------------
!
!.......................................................................
!zj$omp parallel do &
!zj$omp& private(j,i,k)
!.......................................................................
      do k=1,lm
        do j=jms,jme
          do i=ims,ime
            u_phy(i,j,k)=0.
            v_phy(i,j,k)=0.
!
            RTHCUTEN(I,J,K)=0.
            RQCUTEN(I,J,K)=0.
            RQCCUTEN(I,J,K)=0.
            RQRCUTEN(I,J,K)=0.
            RQICUTEN(I,J,K)=0.
            RQSCUTEN(I,J,K)=0.
            RQGCUTEN(I,J,K)=0.
            dudt_phy(i,j,k)=0.
            dvdt_phy(i,j,k)=0.
          enddo
        enddo
!
        do j=jts_b1,jte_b1
          do i=its_b1,ite_b1
            u_phy(i,j,k)=(u(i,j  ,k)+u(i-1,j  ,k) &
                         +u(i,j-1,k)+u(i-1,j-1,k))*0.25
            v_phy(i,j,k)=(v(i,j  ,k)+v(i-1,j  ,k) &
                         +v(i,j-1,k)+v(i-1,j-1,k))*0.25
          ENDDO
        ENDDO
      ENDDO
!.......................................................................
!zj$omp end parallel do
!.......................................................................
!-----------------------------------------------------------------------
!.......................................................................
!zj$omp parallel do                                                       &
!zj$omp private(i,j,k,ql_k)
!.......................................................................
      DO J=JTS,JTE
        DO I=ITS,ITE
          DZ(I,J,LM)=T(I,J,LM)*(.608*Q(I,J,LM)+1.)*R_D &
                    *(PHINT(I,J,LM+1)-PHINT(I,J,LM)) &
                    /(PHMID(I,J,LM)*G)
        ENDDO
!
        DO K=LM-1,1,-1
        DO I=ITS,ITE
          QL_K=MAX(Q(I,J,K),EPSQ)
          DZ(I,J,K)=T(I,J,K)*(.608*QL_K+1.)*R_D &
                    *(PHINT(I,J,K+1)-PHINT(I,J,K)) &
                    /(PHMID(I,J,K)*G)
        ENDDO
        ENDDO
!
      ENDDO
!.......................................................................
!zj$omp end parallel do
!.......................................................................
!
!write(0,*)'A2,A3,A4,cappa,CP,ELIV,ELWV,EPSQ,p608,PQ0,R_D,TIW' &
!,A2,A3,A4,cappa,CP,ELIV,ELWV,EPSQ,p608,PQ0,R_D,TIW

!
!-----------------------------------------------------------------------
!
!***  SINGLE-COLUMN CONVECTION
!
!-----------------------------------------------------------------------
      IF (CU_PHYSICS /= 0) THEN

          cps_select: SELECT CASE(cu_physics)

            CASE (BMJSCHEME)

            call  bmjdrv( &
                         ids,ide,jds,jde &
                        ,ims,ime,jms,jme &
                        ,its,ite,jts,jte,lm &
                        ,its_b1,ite_b1,jts_b1,jte_b1 &
                        ,entrain,newall,newswap,newupup,nodeep &
                        ,a2,a3,a4,cappa,cp,eliv,elwv,epsq,g &
                        ,p608,pq0,r_d,tiw &
                        ,fres,fr,fsl,fss &
                        ,dt,dyh,ntsd,ncnvc &
                        ,raincv,cutop,cubot,dxh,kpbl &
                        ,th,t,q,u_phy,v_phy,dudt_phy,dvdt_phy &
                        ,phint,phmid,exner &
                        ,cldefi,xland,cu_act_flag &
                      ! optional
                        ,rthcuten,rqcuten &
                        )
!-----------------------------------------------------------------------
           CASE (SASSCHEME)
           call sasdrv( &
                       ims,ime,jms,jme &
                      ,its,ite,jts,jte,lm &
                      ,dt,ntsd,ncnvc &
                      ,th,t,sice,omgalf,twbs,qwbs,pblh,u_phy,v_phy & !zj orig u&v
                      ,q,qc,qr,qi,qs,qg &
                      ,f_qc,f_qr,f_qi,f_qs,F_QGr &
                      ,phint,phmid,exner,rr,dz &
                      ,xland,cu_act_flag &
                      ,psgdt &
                      ,raincv,cutop,cubot &
                      ,dudt_phy,dvdt_phy &
                      ! optional
                      ,rthcuten, rqcuten &
                      ,rqccuten, rqrcuten &
                      ,rqicuten, rqscuten &
                      ,rqgcuten  &
                      )
!! 2014-06-19
!! Weiguo Wang added SAS version from HWRF
           CASE (SASHURSCHEME)
           call sasdrv_hur( &
                       ims,ime,jms,jme &
                      ,its,ite,jts,jte,lm &
                      ,dt,ntsd,ncnvc &
                      ,th,t,sice,omgalf,twbs,qwbs,pblh,u_phy,v_phy & !zj orig u&v
                      ,q,qc,qr,qi,qs,qg &
                      ,f_qc,f_qr,f_qi,f_qs,F_QGr &
                      ,phint,phmid,exner,rr,dz &
                      ,xland,cu_act_flag &
                      ,MOMMIX,PGCON,SAS_MASS_FLUX   &   ! hwrf,namelist
                      ,SHALCONV,SHAL_PGCON          &   ! hwrf,namelist
                      ,W_TOT,PSGDT                  &
!                     ,PRATEC                       &   ! hwrf, useful??
                      ,raincv,cutop,cubot &
                      ,dudt_phy,dvdt_phy &
                      ! optional
                      ,rthcuten, rqcuten &
                      ,rqccuten, rqrcuten &
                      ,rqicuten, rqscuten &
                      ,rqgcuten  &
                      )
!!2014-06-19
!!2016-08-29
!! Weiguo Wang added scale-aware SAS version, same as  HWRF
           CASE (SCALECUSCHEME)
           call scalecudrv( &
                       ids,jde         &
                      ,ims,ime,jms,jme &
                      ,its,ite,jts,jte,lm &
                      ,dt,ntsd,ncnvc &
                      ,th,t,sice,omgalf,twbs,qwbs,pblh,u_phy,v_phy & !zj orig u&v
                      ,q,qc,qr,qi,qs,qg &
                      ,f_qc,f_qr,f_qi,f_qs,F_QGr &
                      ,phint,phmid,exner,rr,dz &
                      ,xland,cu_act_flag &
                      ,dxh, dyh           &
                      ,MOMMIX,PGCON,SAS_MASS_FLUX   &   ! hwrf,namelist
                      ,SHALCONV,SHAL_PGCON          &   ! hwrf,namelist
                      ,W_TOT,PSGDT                  &
!                     ,PRATEC                       &   ! hwrf, useful??
                      ,raincv,cutop,cubot &
                      ,dudt_phy,dvdt_phy &
                      ! optional
                      ,rthcuten, rqcuten &
                      ,rqccuten, rqrcuten &
                      ,rqicuten, rqscuten &
                      ,rqgcuten  &
                      )
!!2016-08-29
            CASE DEFAULT

              WRITE( 0 , * ) 'The cumulus option does not exist: cu_physics = ', cu_physics

          END SELECT cps_select

      END IF
!
!-----------------------------------------------------------------------
!
!***  CNVTOP/CNVBOT HOLD THE MAXIMUM VERTICAL LIMITS OF CONVECTIVE CLOUD
!***  BETWEEN HISTORY OUTPUT TIMES.  HBOTS/HTOPS STORE SIMILIAR INFORMATION
!***  FOR SHALLOW (NONPRECIPITATING) CONVECTION, AND HBOTD/HTOPD ARE FOR
!***  DEEP (PRECIPITATING) CONVECTION.
!
      CF_HI=REAL(MINUTES_HISTORY)/60.
      N_TIMSTPS_OUTPUT=NINT(3600.*CF_HI/DT)
      MNTO=MOD(NTSD,N_TIMSTPS_OUTPUT)
!
      IF(MNTO>0.AND.MNTO<=NCNVC)THEN
        DO J=JTS,JTE
        DO I=ITS,ITE
          CNVBOT(I,J)=REAL(LM+1.)
          CNVTOP(I,J)=0.
          HBOTD(I,J)=REAL(LM+1.)
          HTOPD(I,J)=0.
          HBOTS(I,J)=REAL(LM+1.)
          HTOPS(I,J)=0.
        ENDDO
        ENDDO
      ENDIF
!
!-----------------------------------------------------------------------
!.......................................................................
!zj$omp parallel do                                                       &
!zj$omp& private(j,k,i,dtdt,tchange,pcpcol,ncubot,ncutop,QCW,QRain,QCI,QSnow,QGraup)
!.......................................................................
!-----------------------------------------------------------------------
      do j=jts_b1,jte_b1
      do i=its_b1,ite_b1
!-----------------------------------------------------------------------
!
!***  UPDATE TEMPERATURE, SPECIFIC HUMIDITY, AND HEATING.
!
        DO K=1,LM
!
!***  RTHCUTEN IN BMJDRV IS DTDT OVER exner.
!
          DTDT=RTHCUTEN(I,J,K)*exner(I,J,K)
          T(I,J,K)=T(I,J,K)+DTDT*DTCNVC
          Q(I,J,K)=Q(I,J,K)+RQCUTEN(I,J,K)*DTCNVC
          TCUCN(I,J,K)=TCUCN(I,J,K)+DTDT

!!! WANG, 11-2-2010 SAS convection; modified on 11-20-2014 by BSF
sas_test: IF(CONVECTION=='sas') THEN
            QC(I,J,K)=QC(I,J,K)+DTCNVC*RQCCUTEN(I,J,K)
            QCW=QC(I,J,K)
            QRain=0.
            QCI=0.
            QSnow=0.
            QGraup=0.
            IF(F_QR) THEN
              QR(I,J,K)=QR(I,J,K)+DTCNVC*RQRCUTEN(I,J,K)
              QRain=QR(I,J,K)
            ENDIF
            IF(F_QI) THEN
              QI(I,J,K)=QI(I,J,K)+DTCNVC*RQICUTEN(I,J,K)
              QCI=QI(I,J,K)
            ENDIF
            IF(F_QS) THEN
              QS(I,J,K)=QS(I,J,K)+DTCNVC*RQSCUTEN(I,J,K)
              QSnow=QS(I,J,K)
            ENDIF
            IF(F_QGr) THEN
              QG(I,J,K)=QG(I,J,K)+DTCNVC*RQGCUTEN(I,J,K)
              QGraup=QG(I,J,K)
            ENDIF
!-- Couple CWM, F_ice, & F_rain arrays
            CWM(I,J,K)=QCW+QRain+QCI+QSnow+QGraup
            F_ICE(I,J,K)=0.
            F_RAIN(I,J,K)=0.
            IF(CWM(I,J,K)>EPSQ) F_ICE(I,J,K)=(QCI+QSnow+QGraup)/CWM(I,J,K)
            IF(QRain>EPSQ) F_RAIN(I,J,K)=QRain/(QCW+QRain)
          ENDIF  sas_test
!!! wang, 11-2-2010; modified on 11-20-2014 by BSF
!
!zj          TCHANGE=DTDT*DTCNVC
!zj          IF(ABS(TCHANGE)>DTEMP_CHECK)THEN
!zj            WRITE(0,*)'BIG T CHANGE BY CONVECTION:',TCHANGE,' at (',I,',',J,',',K,')'
!zj	  ENDIF
!
        ENDDO

!write(0,*),'t',(rthcuten(i,j,k),k=1,lm)
!write(0,*),'q',(rqcuten(i,j,k),k=1,lm)
!write(0,*),'u',(dudt_phy(i,j,k),k=1,lm)
!write(0,*),'v',(dvdt_phy(i,j,k),k=1,lm)
!write(0,*),'exner',(exner(i,j,k),k=1,lm)


!
!***  UPDATE PRECIPITATION
!
        PCPCOL=RAINCV(I,J)*1.E-3*NCNVC
        PREC(I,J)=PREC(I,J)+PCPCOL
        ACPREC(I,J)=ACPREC(I,J)+PCPCOL
        ACPREC_TOT(I,J)=ACPREC_TOT(I,J)+PCPCOL
        CUPREC(I,J)=CUPREC(I,J)+PCPCOL
        CUPPT(I,J)=CUPPT(I,J)+PCPCOL
        CPRATE(I,J)=PCPCOL
!
!***  SAVE CLOUD TOP AND BOTTOM FOR RADIATION (HTOP/HBOT) AND
!***  FOR OUTPUT (CNVTOP/CNVBOT, HTOPS/HBOTS, HTOPD/HBOTD) ARRAYS.
!***  MUST BE TREATED SEPARATELY FROM EACH OTHER.
!
        NCUTOP=NINT(CUTOP(I,J))
        NCUBOT=NINT(CUBOT(I,J))
!
        IF(NCUTOP>1.AND.NCUTOP<LM+1)THEN
          HTOP(I,J)=MAX(CUTOP(I,J),HTOP(I,J))
          CNVTOP(I,J)=MAX(CUTOP(I,J),CNVTOP(I,J))
          IF(PCPCOL>0.)THEN
            HTOPD(I,J)=MAX(CUTOP(I,J),HTOPD(I,J))
          ELSE
            HTOPS(I,J)=MAX(CUTOP(I,J),HTOPS(I,J))
          ENDIF
        ENDIF
        IF(NCUBOT>0.AND.NCUBOT<LM+1)THEN
          HBOT(I,J)=MIN(CUBOT(I,J),HBOT(I,J))
          CNVBOT(I,J)=MIN(CUBOT(I,J),CNVBOT(I,J))
          IF(PCPCOL>0.)THEN
            HBOTD(I,J)=MIN(CUBOT(I,J),HBOTD(I,J))
          ELSE
            HBOTS(I,J)=MIN(CUBOT(I,J),HBOTS(I,J))
          ENDIF
        ENDIF
!
      ENDDO
      ENDDO
!.......................................................................
!zj$omp end parallel do
!.......................................................................
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE CUCNVC
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      END MODULE MODULE_CONVECTION
!
!-----------------------------------------------------------------------
