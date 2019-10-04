!
      MODULE MODULE_RA_RRTM
!
!-----------------------------------------------------------------------
!
!***  THE RADIATION DRIVERS AND PACKAGES
!
!-----------------------------------------------------------------------
!
      USE MODULE_INCLUDE
!
      use physparam,     only : icldflg, ioznflg, kind_phys, icmphys

      USE MODULE_CONSTANTS, ONLY : R,CP,PI,EPSQ,STBOLT,EP_2
      USE MODULE_MP_FER_HIRES, ONLY : FPVS

      use module_radiation_driver_nmmb,  only : grrad_nmmb,dayparts

      use module_radsw_parameters,  only : topfsw_type, sfcfsw_type
      use module_radlw_parameters,  only : topflw_type, sfcflw_type

      use module_radsw_main_nmmb, only : buggal, buggaloff, buggalon

!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      PRIVATE
!
      PUBLIC :: RRTM, RRTM_INIT
!
!-----------------------------------------------------------------------
!
!--- Used for Gaussian look up tables
!

!Moved to here from out of the RRTM routine to fix races found by Intel Inspector, jm 20131222,20131226
      LOGICAL       :: CNCLD=.TRUE.
      LOGICAL       :: OPER=.TRUE.
!$OMP THREADPRIVATE(CNCLD,OPER)

      REAL, PRIVATE,PARAMETER :: XSDmax=3.1, DXSD=.01
      INTEGER, PRIVATE,PARAMETER :: NXSD=XSDmax/DXSD
      REAL, DIMENSION(NXSD),PRIVATE,SAVE :: AXSD
      REAL, PRIVATE :: RSQR
      LOGICAL, PRIVATE, SAVE :: SDprint=.FALSE.

!-------------------------------
      INTEGER, SAVE, DIMENSION(3)     :: LTOP
      REAL,SAVE,DIMENSION(4) :: PTOPC
!--------------------------------
!
      REAL, PARAMETER ::         &
     &   RHgrd=1.00              & !--- RH (unitless) for onset of condensation
     &,  TRAD_ice=273.15-30.     & !--- Very tunable parameter
     &,  ABSCOEF_W=800.          & !--- Very tunable parameter
     &,  ABSCOEF_I=500.          & !--- Very tunable parameter
     &,  Qconv=0.1e-3            & !--- Very tunable parameter

     &,  CTauCW=ABSCOEF_W*Qconv  &
     &,  CTauCI=ABSCOEF_I*Qconv

!-- Set to TRUE to bogus in small amounts of convective clouds into the
!   input cloud calculations, but only if all of the following conditions 
!   are met:
!     (1) The maximum condensate mixing ratio is < QWmax.
!     (2) Only shallow convection is present, do not apply to deep convection.
!     (3) Only apply if the depth of shallow convection is between 
!         CU_DEEP_MIN (50 hPa) and CU_DEEP_MAX (200 hPa).
!     (4) Convective precipitation rate must be <0.01 mm/h.  
!
      LOGICAL, SAVE :: CUCLD=.FALSE.     &   ! was .TRUE.
     &                ,SUBGRID=.TRUE.
!
!-- After several tuning experiments, a value for QW_CU=0.003 g/kg should 
!   produce a cloud fraction of O(25%) and a SW reduction of O(100 W/m**2) 
!   for shallow convection with a maximum depth/thickness of O(200 hPa).
!-- QW_Cu=0.003 g/kg, which translates to a 3% cloud fraction in 
!   subroutine progcld2 at each model layer near line 960 in 
!   radiation_clouds.f, which translates to a O(25%) total cloud fraction
!   in the lower atmosphere (i.e., for "low-level" cloud fractions).
!
      REAL, PARAMETER :: QW_Cu=0.003E-3,QWmax=1.E-7,CUPPT_min=1.e-5   &
                        ,CU_DEEP_MIN=50.E2,CU_DEEP_MAX=200.E2

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  THE RADIATION PACKAGE OPTIONS
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
      SUBROUTINE RRTM (NTIMESTEP,DT,JDAT                                &
     &                    ,NPHS,GLAT,GLON                               &
     &                    ,NRADS,NRADL                                  &
     &                    ,DSG2,SGML2,PDSG1,PSGML1                      &
     &                    ,PT,PD                                        &
     &                    ,T,Q,CW,O3                                    &
     &                    ,ALBEDO                                       &
     &                    ,F_ICE,F_RAIN                                 &
     &                    ,QC,QI,QS,QR,QG,NI                            &
     &                    ,F_QC,F_QI,F_QS,F_QR,F_QG,F_NI                &
     &                    ,NUM_WATER                                    &
     &                    ,CLD_FRACTION                                 &
     &                    ,SM,CLDFRA                                    &
     &                    ,RLWTT,RSWTT                                  &
     &                    ,RLWIN,RSWIN                                  &
     &                    ,RSWINC,RSWOUT                                &
     &                    ,RLWTOA,RSWTOA                                &
     &                    ,CZMEAN,SIGT4                                 &
     &                    ,CFRACL,CFRACM,CFRACH                         &
     &                    ,ACFRST,NCFRST                                &
     &                    ,ACFRCV,NCFRCV                                &
     &                    ,CUPPT,SNOWC,SI                               & !was SNOW
     &                    ,HTOP,HBOT                                    &
     &                    ,TSKIN,Z0,SICE,F_RIMEF,MXSNAL,SGM,STDH,OMGALF &
     &                    ,IMS,IME,JMS,JME                              &
     &                    ,ITS,ITE,JTS,JTE                              &
     &                    ,LM                                           &
     &                    ,SOLCON                                       &
     &                    ,MYPE )
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      INTEGER,INTENT(IN) :: IME,IMS,ITE,ITS                             &
     &                     ,JME,JMS,JTE,JTS                             &
     &                     ,LM,MYPE                                     &
     &                     ,NTIMESTEP                                   &
     &                     ,NPHS,NRADL,NRADS                            &
     &                     ,CLD_FRACTION                                &
     &                     ,NUM_WATER           !-- not used any more
!
      INTEGER,INTENT(IN) :: JDAT(8)
!
      INTEGER,DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT) :: NCFRCV,NCFRST
!
      REAL,INTENT(IN) :: PT,DT

      real (kind=kind_phys), INTENT(IN) :: SOLCON
!
      REAL,DIMENSION(1:LM),INTENT(IN) :: DSG2,PDSG1,PSGML1,SGML2
!
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: CUPPT               &
                                                   ,GLAT,GLON           &
                                                   ,PD,SM,SNOWC,SI        !was SNOW
!
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: ALBEDO

      REAL,DIMENSION(IMS:IME,JMS:JME,LM),INTENT(IN) :: CW,O3,Q,T
!
      REAL,DIMENSION(IMS:IME,JMS:JME,1:LM),INTENT(IN) :: F_ICE,F_RAIN
!
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT) :: ACFRCV,ACFRST    &
                                                      ,RLWIN,RLWTOA     &
                                                      ,RSWIN,RSWOUT     &
                                                      ,HBOT,HTOP        &
                                                      ,RSWINC,RSWTOA
!
      REAL,DIMENSION(IMS:IME,JMS:JME,LM),INTENT(INOUT) :: RLWTT,RSWTT
!
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT) :: CFRACH,CFRACL    &
                                                      ,CFRACM,CZMEAN    &
                                                      ,SIGT4
!
      LOGICAL,INTENT(IN) :: F_QC,F_QS,F_QI,F_QR,F_QG,F_NI

      REAL,DIMENSION(IMS:IME,JMS:JME,1:LM),INTENT(INOUT) :: QC,QS       &
     &                     ,QI,QR,QG,NI 
!
      REAL,DIMENSION(IMS:IME,JMS:JME,1:LM),INTENT(INOUT) :: CLDFRA
!
       REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: TSKIN,Z0,SICE      &
                                                    ,MXSNAL,STDH        
!
       REAL,DIMENSION(1:LM+1),INTENT(IN) :: SGM
!
       REAL,DIMENSION(IMS:IME,JMS:JME,1:LM),INTENT(IN) :: F_RIMEF,OMGALF
!
       real(kind=kind_phys),DIMENSION(IMS:IME,JMS:JME) :: COSZDG    ! future output
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
#define LENIVEC (ITE-ITS+1)

!
      LOGICAL ::  LSSAV=.TRUE.
                ! logical flag for store 3-d cloud field
                ! ** need to be .TRUE. for non-zero FLUXR_V & CLDCOV_V off GRRAD

      LOGICAL :: LPRNT=.FALSE.
!
      LOGICAL :: LSLWR, LSSWR

      INTEGER,PARAMETER :: NFLUXR=39

!==========================================================================
!  Special for the lwrad to enhence the emissivity
!  it is similar to *CPATHFAC4LW to odcld in radlw  (Hsin-Mu Lin, 20140520)
!==========================================================================

      real(kind=kind_phys), PARAMETER :: CPATHFAC4LW=1.5

!
!-- WARNING: NTRAC must be large enough to account for 
!   different hydrometeor species +2, for ozone & aerosols
!
      INTEGER,PARAMETER :: NTRAC=9   ! GR1 dimension for ozone, aerosol, & clouds
      INTEGER,PARAMETER :: NTCW =3   ! ARRAY INDEX LOCATION FOR CLOUD CONDENSATE
      INTEGER,PARAMETER :: NCLDX=1   ! only used when ntcw .gt. 0

!
      INTEGER :: NUMX, NUMY, NFXR, KFLIP, I, L, J, K

      INTEGER :: ICWP, NTOZ
!
      real(kind=kind_phys),DIMENSION(LENIVEC) :: RTvR
      real(kind=kind_phys) ::  DTSW, DTLW, FHSWR, FHLWR, ARG_CW
!
      real(kind=kind_phys),DIMENSION(LENIVEC) ::                                    &
                             FLGMIN_L, CV, CVB, CVT, HPRIME_V, TSEA,      &
                             TISFC, FICE, ZORL, SLMSK, SNWDPH, SNCOVR,    &
                             SNOALB, ALVSF1, ALNSF1, ALVWF1, ALNWF1,      &
                             FACSF1, FACWF1, SFCNSW, SFCDSW, SFALB,       &
                             SFCDLW, TSFLW, TOAUSW, TOADSW, SFCCDSW,      &
                             TOAULW, SFCUSW, COSZEN_V, COSZDG_V,          &
                             SEMIS, XLAT, XLON, SINLAT, COSLAT

                           !===================================
                           ! SEMIS: surface lw emissivity
                           !        is intended output in GLOOPR
                           !        ** not NMMB in RRTM driver
                           !===================================

      INTEGER, DIMENSION(LENIVEC) :: ICSDSW, ICSDLW

!---  variables of instantaneous calculated toa/sfc radiation fluxes
!      ** IM=1 for the dimension
!
      type (topfsw_type), dimension(LENIVEC) :: TOPFSW
      type (sfcfsw_type), dimension(LENIVEC) :: SFCFSW

      type (topflw_type), dimension(LENIVEC) :: TOPFLW
      type (sfcflw_type), dimension(LENIVEC) :: SFCFLW

      real(kind=kind_phys),DIMENSION(LENIVEC,LM) ::                               &
     &                      CLDCOV_V,PRSL,PRSLK,GT,GQ, VVEL,F_ICEC,     &
                               F_RAINC,R_RIME,TAUCLOUDS,CLDF

      real(kind=kind_phys),DIMENSION(LENIVEC,LM+1) :: PRSI
!
      real(kind=kind_phys),DIMENSION(LENIVEC,5)  :: CLDSA_V
!
      real(kind=kind_phys),DIMENSION(LENIVEC,NFLUXR) :: FLUXR_V
!
      real(kind=kind_phys),DIMENSION(LENIVEC,LM,NTRAC) :: GR1   
!
      real(kind=kind_phys),DIMENSION(LENIVEC,LM) :: SWH, HLW
!
      REAL,DIMENSION(ITS:ITE,1:LM+1) :: P8W   !j dim removed, assume this is only called for 1 j at a time
!
      REAL,DIMENSION(ITS:ITE,1:LM)   :: P_PHY   !j dim removed, assume this is only called for 1 j at a time
!
      INTEGER :: JDOY, JDAY, JDOW, MMM, MMP, MM, IRET, MONEND, &
                 MON1, IS2, ISX, KPD9, IS1, NN, MON2, MON, IS, &  
                 LUGB, LEN, JMSK, IMSK       
!
      REAL :: WV,QICE,QCLD,CLFR,ESAT,QSAT,RHUM,RHtot,ARG,SDM,   &
               PMOD,CONVPRATE,CLSTP,P1,P2,CC1,CC2,CLDMAX,CL1,CL2, &
               CR1,DPCL,PRS1,PRS2,DELP,TCLD,CTau,CFSmax,CFCmax,  &
               CFRAVG,TDUM,CU_DEPTH
!
      INTEGER :: IXSD,NTSPH,NRADPP,NC,NMOD,LCNVT,LCNVB,NLVL,MALVL, &
                 LLTOP,LLBOT,KBT2,KTH1,KBT1,KTH2,KTOP1,LM1,LL
!
      REAL, PARAMETER :: EPSQ1=1.E-5,EPSQ2=1.E-8,EPSO3=1.E-10,H0=0., &
                         H1=1.,HALF=.5,CUPRATE=24.*1000., &
                         HPINC=HALF*1.E1, CLFRmin=0.01, TAUCmax=4.161, &
                         XSDmin=-XSDmax, DXSD1=-DXSD, STSDM=0.01, & 
                         CVSDM=.04,DXSD2=HALF*DXSD,DXSD2N=-DXSD2,PCLDY=0.25
!
      REAL,DIMENSION(10),SAVE :: CC,PPT

! moved out of routine into module and made threadprivate (jm 20131226, see above)
!jm      LOGICAL, SAVE :: CNCLD=.TRUE.
!jm      LOGICAL, SAVE :: OPER=.TRUE.

!
      DATA CC/0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0/
      DATA PPT/0.,.14,.31,.70,1.6,3.4,7.7,17.,38.,85./
!
      REAL,DIMENSION(0:LM)  :: CLDAMT
!
      LOGICAL :: BITX,BITY,BITZ,BITW,BIT1,BIT2,NEW_CLOUD,CU_cloud(LENIVEC)
!
      REAL :: CTHK(3)
      DATA CTHK/20000.0,20000.0,20000.0/
! 
      REAL,DIMENSION(ITS:ITE,JTS:JTE,3):: CLDCFR
      INTEGER,DIMENSION(ITS:ITE,JTS:JTE,3):: MBOT,MTOP

      REAL,DIMENSION(ITS:ITE,JTS:JTE):: CUTOP,CUBOT
      
      REAL,DIMENSION(ITS:ITE,JTS:JTE,LM) :: TauCI,CSMID,CCMID
!
      INTEGER,DIMENSION(ITS:ITE,JTS:JTE,LM+1) :: KTOP, KBTM

      REAL,DIMENSION(ITS:ITE,JTS:JTE,LM+1) :: CAMT
!
      INTEGER,DIMENSION(ITS:ITE,JTS:JTE) :: NCLDS, KCLD
!
      REAL,DIMENSION(ITS:ITE,JTS:JTE,LM) :: TAUTOTAL
!
      INTEGER :: NKTP,NBTM,NCLD,LML, ihr, imin, INDX_CW
!
      REAL :: CLFR1,TauC,DELPTOT

      
      real (kind=kind_phys), parameter :: f24 = 24.0     ! hours/day
      real (kind=kind_phys) :: fhr, solhr, SFCALBEDO(LENIVEC), SMX(LENIVEC)

      integer, dimension(CHNK_RRTM) :: dp_start,dp_len
      logical, dimension(CHNK_RRTM) :: dp_day
      integer                       :: ndayparts, idaypart,dps,dpl
      logical                       :: dpd
      integer im, ii,tid
      integer,external :: omp_get_thread_num
!--------------------------------------------------------------------------------------------------
!
!***THIS SUBROUTINE SELECTS AND PREPARES THE NECESSARY INPUTS FOR GRRAD (GFS RRTM DRIVER)
!
!   GRRAD IS CALLED COLUMN BY COLUMN
!
!INPUTS/OUPUTS OF GRRAD: 
!    INPUT VARIABLES:                                                   !
!      PRSI  (LM+1)    : MODEL LEVEL PRESSURE IN CB (KPA)               !
!      PRSL  (LM)      : MODEL LAYER MEAN PRESSURE IN CB (KPA)          !
!      PRSLK (LM)      : Exner function (dimensionless)                 !
!      GT    (LM)      : MODEL LAYER MEAN TEMPERATURE IN K              !
!      GQ    (LM)      : LAYER SPECIFIC HUMIDITY IN GM/GM               !
!      GR1   (LM,NTRAC): TRACER ARRAY (OZONE, AEROSOL, Various Hydrometeors) !
!      VVEL   (LM)      : LAYER MEAN VERTICAL VELOCITY IN CB/SEC         ! !not used
!      SLMSK (1)       : SEA/LAND MASK ARRAY (SEA:0,LAND:1,SEA-ICE:2)   !
!      XLON,XLAT       : GRID LONGITUDE/LATITUDE IN RADIANS             !
!      TSEA  (1)       : SURFACE TEMPERATURE IN K                       !
!      SNWDPH (1)       : SNOW DEPTH WATER EQUIVALENT IN MM              !
!      SNCOVR(1)       : SNOW COVER IN FRACTION                         !
!      SNOALB(1)       : MAXIMUM SNOW ALBEDO IN FRACTION                !
!      ZORL  (1)       : SURFACE ROUGHNESS IN CM                        !
!      HPRIM_V (1)       : TOPOGRAPHIC STANDARD DEVIATION IN M            !
!      ALVSF1 (1)       : MEAN VIS ALBEDO WITH STRONG COSZ DEPENDENCY    !
!      ALNSF1 (1)       : MEAN NIR ALBEDO WITH STRONG COSZ DEPENDENCY    !
!      ALVWF1 (1)       : MEAN VIS ALBEDO WITH WEAK COSZ DEPENDENCY      !
!      ALNWF1 (1)       : MEAN NIR ALBEDO WITH WEAK COSZ DEPENDENCY      !
!      FACSF1 (1)       : FRACTIONAL COVERAGE WITH STRONG COSZ DEPENDEN  !
!      FACWF1 (1)       : FRACTIONAL COVERAGE WITH WEAK COSZ DEPENDENCY  !
!      FICE  (1)       : ICE FRACTION OVER OPEN WATER GRID              !
!      TISFC (1)       : SURFACE TEMPERATURE OVER ICE FRACTION          !
!      SOLCON          : SOLAR CONSTANT (SUN-EARTH DISTANT ADJUSTED)    !
!
!-- Following 5 quantities are defined within a local 'tile':
!      SINLAT_t        : SINE OF LATITUDE                               !
!      COSLAT_t        : COSINE OF LATITUDE                             !
!      XLON_t          : LONGITUDE                                      !
!      COSZEN_t        : MEAN COS OF ZENITH ANGLE OVER RAD CALL PERIOD  !
!      COSZDG_t        : MEAN COS OF ZENITH ANGLE OVER RAD CALL PERIOD  !
!
!      CV    (1)       : FRACTION OF CONVECTIVE CLOUD                   ! !not used
!      CVT, CVB (1)    : CONVECTIVE CLOUD TOP/BOTTOM PRESSURE IN CB     ! !not used
!      IOVRSW/IOVRLW   : CONTROL FLAG FOR CLOUD OVERLAP (SW/LW RAD)     !
!                        =0 RANDOM OVERLAPPING CLOUDS                   !
!                        =1 MAX/RAN OVERLAPPING CLOUDS                  !
!      F_ICEC (LM)     : FRACTION OF CLOUD ICE  (IN FERRIER SCHEME)     !
!      F_RAINC(LM)     : FRACTION OF RAIN WATER (IN FERRIER SCHEME)     !
!      RRIME  (LM)     : MASS RATIO OF TOTAL TO UNRIMED ICE ( >= 1 )    !
!      FLGMIN_L(1)     : MINIMIM LARGE ICE FRACTION                     !
!                        =8 THOMPSON MICROPHYSICS SCHEME                ! G. Thompson 23Feb2013
!      NTCW            : =0 NO CLOUD CONDENSATE CALCULATED              !
!                        >0 ARRAY INDEX LOCATION FOR CLOUD CONDENSATE   !
!      NCLDX           : ONLY USED WHEN NTCW .GT. 0                     !
!      NTOZ            : =0 CLIMATOLOGICAL OZONE PROFILE                !
!                        >0 INTERACTIVE OZONE PROFILE                   ! !does not work currently
!      NTRAC           : DIMENSION VERIABLE FOR ARRAY GR1               !
!      NFXR            : SECOND DIMENSION OF INPUT/OUTPUT ARRAY FLUXR   !
!      DTLW, DTSW      : TIME DURATION FOR LW/SW RADIATION CALL IN SEC  !
!      LSSAV           : LOGICAL FLAG FOR STORE 3-D CLOUD FIELD         !
!      LM              : VERTICAL LAYER DIMENSION                       !
!      MYPE            : CONTROL FLAG FOR PARALLEL PROCESS              !
!      LPRNT           : CONTROL FLAG FOR DIAGNOSTIC PRINT OUT          !
!      TAUCLOUDS(LM)   : CLOUD OPTICAL DEPTH FROM NMMB (ferrier+bmj)    ! !new
!      CLDF(LM)        : CLOUD FRACTION FROM NMMB (ferrier+bmj)         ! !new
!                                                                       !
!    OUTPUT VARIABLES:                                                  !
!      SWH (LM)       : TOTAL SKY SW HEATING RATE IN K/SEC              !
!      SFCNSW(1)      : TOTAL SKY SURFACE NET SW FLUX IN W/M**2         !
!      SFCDSW(1)      : TOTAL SKY SURFACE DOWNWARD SW FLUX IN W/M**2    !
!      SFALB (1)      : MEAN SURFACE DIFFUSED ALBEDO                    !
!      HLW (LM)       : TOTAL SKY LW HEATING RATE IN K/SEC              !
!      SFCDLW(1)      : TOTAL SKY SURFACE DOWNWARD LW FLUX IN W/M**2    !
!      TSFLW (1)      : SURFACE AIR TEMP DURING LW CALCULATION IN K     !
!
!      TOAUSW (IM)    : TOTAL SKY TOA UPWARD SW FLUX IN W/M**2         ! !new
!      TOADSW (IM)    : TOTAL SKY TOA DOWNWARD SW FLUX IN W/M**2       ! !new
!      SFCCDSW(IM)    : CLEAR SKY SURFACE SW DOWNWARD FLUX IN W/M**2   ! !new
!      TOAULW (IM)    : TOTAL SKY TOA LW FLUX W/M**2                   ! !new
!      SFCUSW (IM)    : TOTAL SKY SURFACE SW UPWARD FLUX IN W/M**2     ! !new
!                                                                       !
!    INPUT AND OUTPUT VARIABLES:                                        !
!      FLUXR_V (IX,NFXR) : TO SAVE 2-D FIELDS                           !
!                          (bucket)                                     !
!      CLDSA_V(IX,5)     : TO SAVE 2-D CLOUD FRACTION. L/M/H/TOT/BL     !
!                          (instantaneous)                              !
!      CLDCOV_V(IX,LM)   : TO SAVE 3-D CLOUD FRACTION                   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!SELECT OPTIONS IN GRRAD
!
      NFXR=NFLUXR    ! second dimension of input/output array fluxr (FLUXR_V)

      ICSDSW(:)=0    ! auxiliary special cloud related array for SW
                     ! *** not used in this version of code ***
                     ! can be any value at this moment
      ICSDLW(:)=0    ! auxiliary special cloud related array for LW
                     ! *** not used in this version of code ***
                     ! can be any value at this moment

      ICWP = icldflg
      NTOZ = ioznflg

!------------------------------
! for np3d=5 (Lin, 20150601)
!------------------------------

      IF (ICMPHYS == 5 ) THEN
         ICWP = -1
      ENDIF
!
!=========================================================================
!
      IF (ICWP/=-1 .AND. CNCLD) THEN
         CNCLD=.FALSE.        !-- used when ICWP=1, 0
      ENDIF
!
!--- Cloud water index for the GR1 array
!
      IF (F_NI) THEN
        INDX_CW=4  !-- Thompson 
      ELSE
        INDX_CW=3  !-- All others (as of Oct 2014)
      ENDIF
!
!CLOUDS
!
!----------------------CONVECTION--------------------------------------
!  NRADPP IS THE NUMBER OF TIME STEPS TO ACCUMULATE CONVECTIVE PRECIP
!     FOR RADIATION
!   NOTE: THIS WILL NOT WORK IF NRADS AND NRADL ARE DIFFERENT UNLESS
!         THEY ARE INTEGER MULTIPLES OF EACH OTHER
!  CLSTP IS THE NUMBER OF HOURS OF THE ACCUMULATION PERIOD
!
      NTSPH=NINT(3600./DT)
      NRADPP=MIN(NRADS,NRADL)
      CLSTP=1.0*NRADPP/NTSPH
      CONVPRATE=CUPRATE/CLSTP

      IF (ICWP>0 .AND. CUCLD) CONVPRATE=1000./CLSTP    !-- convert to mm/h
!
      LM1=LM-1
!
      DO J=JTS,JTE
      DO I=ITS,ITE
!
        P8W(I,1)=PT
!
        DO K=1,LM
          P8W(I,K+1)=P8W(I,K)+PDSG1(K)+DSG2(K)*PD(I,J)
          P_PHY(I,K)=SGML2(K)*PD(I,J)+PSGML1(K)
          CCMID(I,J,K)=0.
          CSMID(I,J,K)=0.
        ENDDO
      ENDDO
      ENDDO
!
! --- initialize for non Thompson cloud fraction (used only in gfdl type)
!     for thompson cloud fraction, "CLDFRC" is direct INPUT
!     
      IF (CLD_FRACTION==0) THEN
         DO K=1,LM
         DO J=JTS,JTE
         DO I=ITS,ITE
            CLDFRA(I,J,K)=0.
         ENDDO
         ENDDO
         ENDDO
      ENDIF

! ---- 

      DO K=1,LM
      DO J=JTS,JTE
      DO I=ITS,ITE
         TAUTOTAL(I,J,K)=0.
      ENDDO
      ENDDO
      ENDDO
!
      DO J=JTS,JTE
      DO I=ITS,ITE
          CFRACH(I,J)=0.
          CFRACL(I,J)=0.
          CFRACM(I,J)=0.
          CZMEAN(I,J)=0.
          SIGT4(I,J)=0.
      ENDDO
      ENDDO
!
      DO K=1,3
      DO J=JTS,JTE
      DO I=ITS,ITE
        CLDCFR(I,J,K)=0.
        MTOP(I,J,K)=0
        MBOT(I,J,K)=0
      ENDDO
      ENDDO
      ENDDO
!
      DO J=JTS,JTE
      DO I=ITS,ITE
       CUTOP(I,J)=LM+1-HTOP(I,J)
       CUBOT(I,J)=LM+1-HBOT(I,J)
      ENDDO
      ENDDO
!      
!-----------------------------------------------------------------------
!---  COMPUTE GRID-SCALE CLOUD COVER FOR RADIATION  (Ferrier, Nov '04)
!
!--- Assumes Gaussian-distributed probability density functions (PDFs) for
!    total relative humidity (RHtot) within the grid for convective and
!    grid-scale cloud processes.  The standard deviation of RHtot is assumed
!    to be larger for convective clouds than grid-scale (stratiform) clouds.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
      ICWP_Test: IF (ICWP==-1) THEN   !-- *** Start of old NAM/GFDL cloud inputs ***
!-----------------------------------------------------------------------
       DO J=JTS,JTE
       DO I=ITS,ITE 
!
        DO 255 L=1,LM
!
          WV=MAX(EPSQ,Q(I,J,L))/(1.-MAX(EPSQ,Q(I,J,L)))   !-- Water vapor mixing ratio
          QICE=MAX(QS(I,J,L),0.)                          !-- Ice mixing ratio
          QCLD=QICE+MAX(QS(I,J,L),0.)                     !-- Total cloud water + ice mixing ratio
!rv------------------------------------
!rv   This should be temporary fix!!!!!
!rv   New (currently operational) calculation of cloud fraction is
!rv   causing different results with different decomposition
!rv   We should find cause of this!!!!!
!rv------------------------------------
          OPER_flag: IF (OPER) THEN
!rv------------------------------------
!-- From model tuning experiments vs CLAVR grid-to-grid verification:
!-- 100% cloud fractions at 0.01 g/kg (1.e-5 kg/kg) cloud mixing ratios
!-- 10% cloud fractions at 1.e-4 g/kg (1.e-7 kg/kg) cloud mixing ratios
!-- 1% cloud fractions at 1.e-6 g/kg (1.e-9 kg/kg) cloud mixing ratios
!
            CLFR=MIN(H1, MAX(H0,1.e5*QCLD))
            CLFR=SQRT(CLFR)
            IF (CLFR>=CLFRmin) CSMID(I,J,L)=CLFR
!rv------------------------------------
          else OPER_flag
!rv------------------------------------

!
            IF (QCLD .LE. EPSQ) GO TO 255                               !--- Skip if no condensate is present
            CLFR=H0
!
            WV=MAX(EPSQ,Q(I,J,L))/(1.-MAX(EPSQ,Q(I,J,L)))
!
!--- Saturation vapor pressure w/r/t water ( >=0C ) or ice ( <0C )
!
            ESAT=1000.*FPVS(T(I,J,L))                                   !--- Saturation vapor pressure (Pa)
            ESAT=MIN(ESAT, 0.99*P_PHY(I,L) )                            !--- Put limits on ESAT
            QSAT=EP_2*ESAT/(P_PHY(I,L)-ESAT)                            !--- Saturation mixing ratio
!
            RHUM=WV/QSAT                                                !--- Relative humidity
!
!--- Revised cloud cover parameterization (temporarily ignore rain)
!
            RHtot=(WV+QCLD)/QSAT                                        !--- Total relative humidity
!
            LCNVT=NINT(CUTOP(I,J))
            LCNVT=MIN(LM,LCNVT)
            LCNVB=NINT(CUBOT(I,J))
            LCNVB=MIN(LM,LCNVB)
            IF (L.GE.LCNVT .AND. L.LE.LCNVB) THEN
               SDM=CVSDM
            ELSE
               SDM=STSDM
            ENDIF
            ARG=(RHtot-RHgrd)/SDM
            IF (ARG.LE.DXSD2 .AND. ARG.GE.DXSD2N) THEN
               CLFR=HALF
            ELSE IF (ARG .GT. DXSD2) THEN
               IF (ARG .GE. XSDmax) THEN
                  CLFR=H1
               ELSE
                  IXSD=INT(ARG/DXSD+HALF)
                  IXSD=MIN(NXSD, MAX(IXSD,1))
                  CLFR=HALF+AXSD(IXSD)
               ENDIF              !--- End IF (ARG .GE. XSDmax)
            ELSE
               IF (ARG .LE. XSDmin) THEN
                  CLFR=H0
               ELSE
                  IXSD=INT(ARG/DXSD1+HALF)
                  IXSD=MIN(NXSD, MAX(IXSD,1))
                  CLFR=HALF-AXSD(IXSD)
                  IF (CLFR .LT. CLFRmin) CLFR=H0
               ENDIF        !--- End IF (ARG .LE. XSDmin)
            ENDIF           !--- IF (ARG.LE.DXSD2 .AND. ARG.GE.DXSD2N)
            CSMID(I,J,L)=CLFR
!rv------------------------------------
          endif  OPER_flag
!rv------------------------------------
!
255     CONTINUE         !--- End DO L=1,LM

       ENDDO ! End DO I=ITS,ITE
       ENDDO ! End DO J=JTS,JTE

!***********************************************************************
!******************  END OF GRID-SCALE CLOUD FRACTIONS  ****************


!***********************************************************************
!---  COMPUTE CONVECTIVE CLOUD COVER FOR RADIATION
!
!--- The parameterization of Slingo (1987, QJRMS, Table 1, p. 904) is
!    used for convective cloud fraction as a function of precipitation
!    rate.  Cloud fractions have been increased by 20% for each rainrate
!    interval so that shallow, nonprecipitating convection is ascribed a
!    constant cloud fraction of 0.1  (Ferrier, Feb '02).
!***********************************************************************
!
       GFDL_Conv: IF (CNCLD) THEN

        DO J=JTS,JTE
         DO I=ITS,ITE
!
!***  CLOUD TOPS AND BOTTOMS COME FROM CUCNVC
!     Convective clouds need to be at least 2 model layers thick
!
          IF (CUBOT(I,J)-CUTOP(I,J) .GT. 1.0) THEN
!--- Compute convective cloud fractions if appropriate  (Ferrier, Feb '02)
            CLFR=CC(1)
            PMOD=CUPPT(I,J)*CONVPRATE
            IF (PMOD .GT. PPT(1)) THEN
              DO NC=1,10
                IF(PMOD.GT.PPT(NC)) NMOD=NC
              ENDDO
              IF (NMOD .GE. 10) THEN
                CLFR=CC(10)
              ELSE
                CC1=CC(NMOD)
                CC2=CC(NMOD+1)
                P1=PPT(NMOD)
                P2=PPT(NMOD+1)
                CLFR=CC1+(CC2-CC1)*(PMOD-P1)/(P2-P1)
              ENDIF      !--- End IF (NMOD .GE. 10) ...
              CLFR=MIN(H1, CLFR)
            ENDIF        !--- End IF (PMOD .GT. PPT(1)) ...
!
!***  ADD LVL TO BE CONSISTENT WITH OTHER WORKING ARRAYS
!
            LCNVT=NINT(CUTOP(I,J))
            LCNVT=MIN(LM,LCNVT)
            LCNVB=NINT(CUBOT(I,J))
            LCNVB=MIN(LM,LCNVB)
!
!--- Build in small amounts of subgrid-scale convective condensate
!    (simple assumptions), but only if the convective cloud fraction
!    exceeds that of the grid-scale cloud fraction
!
            DO L=LCNVT,LCNVB
              ARG=MAX(H0, H1-CSMID(I,J,L))
              CCMID(I,J,L)=MIN(ARG,CLFR)
            ENDDO           !--- End DO LL=LCNVT,LCNVB
          ENDIF             !--- IF (CUBOT(I,J)-CUTOP(I,J) .GT. 1.0) ...
         ENDDO               ! End DO I=ITS,ITE
        ENDDO                ! End DO J=JTS,JTE
       ENDIF  GFDL_Conv      !--- End IF (CNCLD) ...
!
!*********************************************************************
!***************  END OF CONVECTIVE CLOUD FRACTIONS  *****************
!*********************************************************************
!***
!*** INITIALIZE ARRAYS FOR USES LATER
!***

       DO I=ITS,ITE
       DO J=JTS,JTE
!
         LML=LM
!***
!*** NOTE: LAYER=1 IS THE SURFACE, AND LAYER=2 IS THE FIRST CLOUD
!***       LAYER ABOVE THE SURFACE AND SO ON.
!***
         KTOP(I,J,1)=LM+1
         KBTM(I,J,1)=LM+1
         CAMT(I,J,1)=1.0
         KCLD(I,J)=2
!
         DO 510 L=2,LM+1
           CAMT(I,J,L)=0.0
           KTOP(I,J,L)=1
           KBTM(I,J,L)=1
  510    CONTINUE
!### End changes so far
!***
!*** NOW CALCULATE THE AMOUNT, TOP, BOTTOM AND TYPE OF EACH CLOUD LAYER
!*** CLOUD TYPE=1: STRATIFORM CLOUD
!***       TYPE=2: CONVECTIVE CLOUD
!*** WHEN BOTH CONVECTIVE AND STRATIFORM CLOUDS EXIST AT THE SAME POINT,
!*** SELECT CONVECTIVE CLOUD WITH THE HIGHER CLOUD FRACTION.
!*** CLOUD LAYERS ARE SEPARATED BY TOTAL ABSENCE OF CLOUDINESS.
!*** NOTE: THERE IS ONLY ONE CONVECTIVE CLOUD LAYER IN ONE COLUMN.
!*** KTOP AND KBTM ARE THE TOP AND BOTTOM OF EACH CLOUD LAYER IN TERMS
!*** OF MODEL LEVEL.
!***
         NEW_CLOUD=.TRUE.
!
      DO L=2,LML
        LL=LML-L+1                                  !-- Model layer
        CLFR=MAX(CCMID(I,J,LL),CSMID(I,J,LL))       !-- Cloud fraction in layer
        CLFR1=MAX(CCMID(I,J,LL+1),CSMID(I,J,LL+1))  !-- Cloud fraction in lower layer
!-------------------
        IF (CLFR .GE. CLFRMIN) THEN
!--- Cloud present at level
          IF (NEW_CLOUD) THEN
!--- New cloud layer
            IF(L==2.AND.CLFR1>=CLFRmin)THEN
              KBTM(I,J,KCLD(I,J))=LL+1
              CAMT(I,J,KCLD(I,J))=CLFR1
            ELSE
              KBTM(I,J,KCLD(I,J))=LL
              CAMT(I,J,KCLD(I,J))=CLFR
            ENDIF
            NEW_CLOUD=.FALSE.
          ELSE
!--- Existing cloud layer
            CAMT(I,J,KCLD(I,J))=AMAX1(CAMT(I,J,KCLD(I,J)), CLFR)
          ENDIF        ! End IF (NEW_CLOUD .EQ. 0) ...
        ELSE IF (CLFR1 .GE. CLFRMIN) THEN
!--- Cloud is not present at level but did exist at lower level, then ...
          IF (L .EQ. 2) THEN
!--- For the case of ground fog
           KBTM(I,J,KCLD(I,J))=LL+1
           CAMT(I,J,KCLD(I,J))=CLFR1
          ENDIF
          KTOP(I,J,KCLD(I,J))=LL+1
          NEW_CLOUD=.TRUE.
          KCLD(I,J)=KCLD(I,J)+1
          CAMT(I,J,KCLD(I,J))=0.0
        ENDIF
!-------------------
      ENDDO      !--- End DO L loop
!***
!*** THE REAL NUMBER OF CLOUD LAYERS IS (THE FIRST IS THE GROUND;
!*** THE LAST IS THE SKY):
!***
      NCLDS(I,J)=KCLD(I,J)-2
      NCLD=NCLDS(I,J)
!***
!***  NOW CALCULATE CLOUD RADIATIVE PROPERTIES
!***
      IF(NCLD.GE.1)THEN
!***
!*** NOTE: THE FOLLOWING CALCULATIONS, THE UNIT FOR PRESSURE IS MB!!!
!***
        DO NC=2,NCLD+1
!
        TauC=0.    !--- Total optical depth for each cloud layer (solar & longwave)
        NKTP=LM+1
        NBTM=0
        BITX=CAMT(I,J,NC).GE.CLFRMIN
        NKTP=MIN(NKTP,KTOP(I,J,NC))
        NBTM=MAX(NBTM,KBTM(I,J,NC))
!
        DO LL=NKTP,NBTM
          L=NBTM-LL+NKTP 
          IF(LL.GE.KTOP(I,J,NC).AND.LL.LE.KBTM(I,J,NC).AND.BITX)THEN
            PRS1=P8W(I,L)*0.01 
            PRS2=P8W(I,L+1)*0.01
            DELP=PRS2-PRS1
!
            CTau=0.
!-- For crude estimation of convective cloud optical depths
            IF (CCMID(I,J,L) .GE. CLFRmin) THEN
              IF (T(I,J,L) .GE. TRAD_ice) THEN
                CTau=CTauCW            !--- Convective cloud water
              ELSE
                CTau=CTauCI            !--- Convective ice
              ENDIF
            ENDIF
!
!-- For crude estimation of grid-scale cloud optical depths
!
!--   => The following 2 lines were intended to reduce cloud optical depths further
!        than what's parameterized in the NAM and what's theoretically justified
            CTau=CTau+ABSCOEF_W*QC(I,J,L)+ABSCOEF_I*QS(I,J,L)

            TAUTOTAL(I,J,L)=CTau*DELP                          !Total model level cloud optical depth
            CLDFRA(I,J,L)=MAX(CCMID(I,J,LL),CSMID(I,J,LL))     !Cloud fraction at model level           
            TauC=TauC+DELP*CTau                                !Total cloud optical depth as in GFDL
!
          ENDIF      !--- End IF(LL.GE.KTOP(I,NC) ....
        ENDDO        !--- End DO LL
!
      ENDDO
!
      ENDIF       ! NCLD.GE.1
!
      ENDDO  !  DO I=ITS,ITE
      ENDDO  !  DO J=JTS,JTE
!-----------------------------------------------------------------------
      ENDIF  ICWP_Test   !*** End of Old NAM/GFDL cloud inputs ***
!-----------------------------------------------------------------------

      FHSWR=(NRADS*DT)/3600.        ! [h]
      FHLWR=(NRADL*DT)/3600.        ! [h]
      DTLW =(NRADL*DT)              ! [s]
      DTSW =(NRADS*DT)              ! [s]
      LSSWR=MOD(NTIMESTEP,NRADS)==0
      LSLWR=MOD(NTIMESTEP,NRADL)==0

!==========================================================================
!  Similar to GFS "gloopr.f" line #370,  #413
!  The following block is from old "radiation_astronomy_nmmb.f"
!==========================================================================

      ihr   = JDAT(5)
      imin  = JDAT(6)

!  --- ...  hour of forecast time

      !  solhr = mod( float(ihr), f24 )    ! previous version

      !=== the new calculatuion will eliminate the time lag due to
      !    "jdate(5)" handled by ESMF  (201208)

      fhr = float(ihr)+float(imin)/60.
      solhr = mod( fhr, f24 )


!==========================================================================
! Main domain loop: calling grrad
!==========================================================================
!     

      DO J=JTS,JTE  !start grrad loop column by column

      XLON(1:LENIVEC) = GLON(its:ite,J)
      XLAT(1:LENIVEC)=GLAT(its:ite,J)

      SINLAT(1:LENIVEC) = SIN ( XLAT(1:LENIVEC) )
      COSLAT(1:LENIVEC) = COS ( XLAT(1:LENIVEC) )

      call dayparts( XLON,SINLAT,COSLAT,SOLHR,MYPE,FHSWR,NRADS,         &
                     LENIVEC,                                           &
                     dp_start,dp_len,dp_day,ndayparts )

      dayparts: DO idaypart=1,ndayparts
        dps = dp_start(idaypart)
        dpl = dp_len(idaypart)
        dpd = dp_day(idaypart)  !logical
        im = dpl

#define IBEG (its+dps-1)
#define IDEX (its+dps-1+I-1)
#define IEND (its+dps-1+im-1)
#define IRANGE IBEG:IEND 
#define IITER  IBEG,IEND

       ! if ( GLON(I,J) >= 0.0 ) then
           XLON(1:im) = GLON(IRANGE,J)
       ! else
       !    XLON(1) = GLON(I,J) + PI        ! if in -pi->+pi convert to 0->2pi
       ! endif

       XLAT(1:im)=GLAT(IRANGE,J)
       SINLAT(1:im) = SIN ( XLAT(1:im) )
       COSLAT(1:im) = COS ( XLAT(1:im) )

       TSEA(1:im)=TSKIN(IRANGE,J)
       TISFC(1:im)=TSKIN(IRANGE,J)                  ! change later if necessary
       ZORL(1:im)=Z0(IRANGE,J)*100.d0
       SNWDPH(1:im)=SI(IRANGE,J)                    ! snwdph[mm]
       SNCOVR(1:im)=SNOWC(IRANGE,J)
       SNOALB(1:im)=MXSNAL(IRANGE,J)
       HPRIME_V(1:im)=STDH(IRANGE,J)

       WHERE (SICE(IRANGE,J).GT.0.5)                 ! slmsk - ocean  - 0
         SLMSK(1:im)= 2.0d0                    !         land   - 1
         FICE(1:im)=SICE(IRANGE,J)                  ! change this later
       ELSEWHERE                            !         seaice - 2
         SLMSK(1:im)= 1.0d0-SM(IRANGE,J)            !
         FICE(1:im)= 0.0d0                     ! change this later
       ENDWHERE
!
!---
!
      FLGMIN_L(1:im)= 0.20d0 ! --- for ferrier

      CV (1:im)=0.d0         ! not in use
      CVB(1:im)=0.d0         ! not in use
      CVT(1:im)=0.d0         ! not in use

      PRSI(1:im,1)=P8W(IRANGE,1)/1000.                                ! [kPa]
!
      DO L=1,LM
        PRSI(1:im,L+1)=P8W(IRANGE,L+1)/1000.                          ! (pressure on interface) [kPa]
        PRSL(1:im,L)=P_PHY(IRANGE,L)/1000.                            ! (pressure on mid-layer) [kPa] 
        PRSLK(1:im,L)=(PRSL(1:im,L)*0.01d0)**(R/CP)
        RTvR(1:im)=1./(R*(Q(IRANGE,J,L)*0.608+1.- &
                   CW(IRANGE,J,L))*T(IRANGE,J,L))
        VVEL(1:im,L)=OMGALF(IRANGE,J,L)*1000.d0   &
                   *PRSL(1:im,L)*RTvR(1:im)            !not used
        GT(1:im,L)=T(IRANGE,J,L)
        GQ(1:im,L)=Q(IRANGE,J,L)
!
!--- GR1(:,:,1) - ozone
!    GR1(:,:,2) - reserved for prognostic aerosols in the future
!    GR1(:,:,3) - total condensate
!    GR1(:,:,4-9) - hydrometeor species from Thompson scheme
!
        DO ii=1,NTRAC
          GR1(1:im,L,ii)=0.d0
        ENDDO
!
        IF (NTOZ>0) GR1(1:im,l,1)=MAX(O3(IRANGE,J,L),EPSO3)
!
        CLDCOV_V(1:im,L)=0.d0                     ! used for prognostic cloud
        TAUCLOUDS(1:im,L)=TAUTOTAL(IRANGE,J,L)    ! CLOUD OPTICAL DEPTH (ICWP==-1)
        CLDF(1:im,L)=CLDFRA(IRANGE,J,L)           ! CLOUD FRACTION

        GR1(1:im,L,3)=CW(IRANGE,J,L)              ! total condensate
!
!----------
!
        thompson_test: IF (F_NI) THEN
!
!----------
!-- IF F_NI=true, then this means the microphysics is Thompson only.
!   Other progcld"X" drivers must be introduced into grrad_nmmb.f
!   in order to use GR1(:,:,N) where N>=4.  (BSF, Oct 2014)
!
!-- Warnings from Thompson:
!.. This is an awful way to deal with different physics having different
!.. number of species.  Something must eventually be done to resolve this
!.. section to be more flexible.  For now, we are directly passing each
!.. species in the water array into the GR1 array for use in the RRTM
!.. radiation scheme, but that requires a priori knowledge of which species
!.. is which index number at some later time.  This is far from optimal,
!.. but we proceed anyway.  Future developers be careful.
!.. If the WATER species include separate hydrometeor species, then
!.. fill in other elements even if unused.  Thompson microphysics
!.. will utilize certain elements when computing cloud optical depth.
!----------
!
           GR1(1:im,L,4)=QC(IRANGE,J,L)
           GR1(1:im,L,5)=QI(IRANGE,J,L)
           GR1(1:im,L,6)=QS(IRANGE,J,L)
           GR1(1:im,L,7)=QR(IRANGE,J,L)
           GR1(1:im,L,8)=QG(IRANGE,J,L)
           GR1(1:im,L,9)=NI(IRANGE,J,L)
!----------
        ELSE  thompson_test
!----------
           F_ICEC(1:im,L)=F_ICE(IRANGE,J,L)
           F_RAINC(1:im,L)=F_RAIN(IRANGE,J,L)
           R_RIME(1:im,L)=F_RIMEF(IRANGE,J,L)
!----------
        ENDIF  thompson_test
!----------

      ENDDO
!
!
      subgrid_cloud: IF (SUBGRID) THEN
!
!-- Build in tiny amounts of subgrid-scale cloud when no cloud is
!   present and RH > 95%.
!-- Note GR1(ii,L,3) is total condensate for all microphysics schemes
!
        thompson_testx: IF (F_NI) THEN
!
!-- Build in tiny amounts of subgrid-scale cloud for the Thompson scheme
!
          DO L=1,LM
          DO I=1,IM
            WV=GQ(I,L)/(1.-GQ(I,L))                !- Water vapor mixing ratio
            TCLD=REAL(GT(I,L))                     !- Temperature (deg K)
            ESAT=FPVS(TCLD)                        !- Saturation vapor pressure (kPa)
            P1=REAL(PRSL(I,L))                     !- Pressure (kPa)
            ESAT=MIN(ESAT, 0.99*P1)                !- Limit saturation vapor pressure

IF(P1<1.E-2) WRITE(6,"(a,3i4,2g11.4)") 'I,J,L,PRSL,E_sat=',I,J,L,P1,ESAT   !dbg

            QSAT=EP_2*ESAT/(P1-ESAT)               !- Saturation mixing ratio
            RHUM=WV/QSAT                           !- Relative humidity
            IF (GR1(I,L,3)<EPSQ .AND. RHUM>0.95) THEN
              ARG=MIN(0.01, RHUM-0.95)*QSAT
              GR1(I,L,3)=MIN(0.01E-3, ARG)
              IF (TCLD>TRAD_ICE) THEN
                GR1(I,L,4)=GR1(I,L,3)
              ELSE
                GR1(I,L,5)=GR1(I,L,3)
              ENDIF
            ENDIF      !- IF (GR1(ii,L,3)<EPSQ ...
          ENDDO        !- DO I
          ENDDO        !- DO L
!
        ELSE  thompson_testx
!
!-- Build in tiny amounts of subgrid-scale cloud for other microphysics schemes
!
          DO L=1,LM 
          DO I=1,IM
            WV=GQ(I,L)/(1.-GQ(I,L))                !- Water vapor mixing ratio
            TCLD=REAL(GT(I,L))                     !- Temperature (deg K)
            ESAT=FPVS(TCLD)                        !- Saturation vapor pressure (kPa)
            P1=REAL(PRSL(I,L))                     !- Pressure (kPa)
            ESAT=MIN(ESAT, 0.99*P1)                !- Limit saturation vapor pressure

IF(P1<1.E-2) WRITE(6,"(a,3i4,2g11.4)") 'I,J,L,PRSL,E_sat=',I,J,L,P1,ESAT   !dbg

            QSAT=EP_2*ESAT/(P1-ESAT)               !- Saturation mixing ratio
            RHUM=WV/QSAT                           !- Relative humidity
            IF (GR1(I,L,3)<EPSQ .AND. RHUM>0.95) THEN
              ARG=MIN(0.01, RHUM-0.95)*QSAT
              GR1(I,L,3)=MIN(0.01E-3, ARG)
              IF (TCLD>TRAD_ICE) THEN
                F_ICEC(I,L)=0.
              ELSE
                F_ICEC(I,L)=1.
              ENDIF
              F_RAINC(I,L)=0.
              R_RIME(I,L)=1.
            ENDIF
          ENDDO
          ENDDO

        ENDIF  thompson_testx

      ENDIF  subgrid_cloud
!
!-- Bogus in tiny amounts of shallow convection, but only if there are no
!   grid-scale clouds nor convective precipitation present.  Arrays CUTOP,
!   CUBOT are flipped with 1 at the top & LM at the surface (BSF, 7/18/2012)
!-- There are extra, nested IF statements to filter conditions as an extra
!   layer of caution.  
!
      CU_cloud=.FALSE.
      CU_Bogus1: IF (CUCLD) THEN
       DO I=IBEG,IEND     !-- to be consistent with IM form "daypart" (ver. 44371, 20140815)
         ii = i-IBEG+1
         LCNVT=MIN(LM, NINT(CUTOP(I,J)) )   !-- Convective cloud top
         LCNVB=MIN(LM, NINT(CUBOT(I,J)) )   !-- Convective cloud base
         CU_DEPTH=0.
         CU_Index: IF (LCNVB-LCNVT>1) THEN
            CU_DEPTH=1000.*(PRSL(ii,LCNVB)-PRSL(ii,LCNVT))   !- Pa
            CU_Deep: IF (CU_DEPTH>=CU_DEEP_MIN .AND. CU_DEPTH<=CU_DEEP_MAX) THEN
               QCLD=MAXVAL( GR1(ii,1:LM,3) )  !- Maximum *total condensate*
               PMOD=CUPPT(I,J)*CONVPRATE
               CU_Clds: IF (QCLD<QWmax .AND. PMOD<=CUPPT_min) THEN
                  CU_cloud(ii)=.TRUE.
                  DO L=LCNVT,LCNVB
                    GR1(ii,L,3)=GR1(ii,L,3)+QW_Cu  !- for *cloud water*
                    IF (F_NI) GR1(ii,L,INDX_CW)=GR1(ii,L,INDX_CW)+QW_Cu    !- *total condensate* in Thompson
                  ENDDO
               ENDIF CU_Clds
            ENDIF CU_Deep
         ENDIF CU_Index
       ENDDO
      ENDIF CU_Bogus1
!
      DO NC=1,5
        CLDSA_V(1:im,NC)=0.d0                 !used for prognostic cloud
      ENDDO
      DO NC=1,NFLUXR
        FLUXR_V(1:im,NC)=0.d0                 !used for prognostic cloud
      ENDDO
!
!---

      SFCALBEDO(1:im) = ALBEDO(IRANGE,J)
      SMX(1:im) = SM(IRANGE,J)

!  --- ...  calling radiation driver

!      write(0,'("ibeg,jts,LENIVEC,im,dpd ",5i7)')IBEG,jts,LENIVEC,im,dpd 
      call grrad_nmmb                                                   &
!!  ---  inputs:
           ( PRSI,PRSL,PRSLK,GT,GQ,GR1,VVEL,SLMSK,                      &
             XLON,XLAT,TSEA,SNWDPH,SNCOVR,SNOALB,ZORL,HPRIME_V,         &
           !  ALVSF1,ALNSF1,ALVWF1,ALNWF1,FACSF1,FACWF1,                &  ! processed inside grrad
             SFCALBEDO,SMX,                                             &  ! input for albedo cal
             FICE,TISFC,                                                &
             SINLAT,COSLAT,SOLHR, JDAT, SOLCON,                         &
             FHSWR ,NRADS,                                              &  ! extra input
             CV,CVT,CVB, F_ICEC, F_RAINC, R_RIME, FLGMIN_L,             &
             ICSDSW,ICSDLW,NTCW,NCLDX,NTOZ,NTRAC,NFXR,                  &  ! Use NTRAC instead of NUM_WATER
             CPATHFAC4LW,                                               &  ! enhance factor of cloud depth for LW
             DTLW,DTSW,LSSWR,LSLWR,LSSAV,                               &
             IBEG,jts, LENIVEC, im, LM, dpd, MYPE, LPRNT, 0, 0,         &  ! jm dpd is true for day, false for night
!  ---  additional inputs:                                                 ! GFDL type
             TAUCLOUDS,CLDF,                                            &  ! GFDL type
             CLD_FRACTION,                                              &  ! Thompson cloud fraction
!!  ---  outputs:
             SWH,TOPFSW,SFCFSW,SFALB,COSZEN_V,COSZDG_V,                 &
             HLW,TOPFLW,SFCFLW,TSFLW,SEMIS,CLDCOV_V,CLDSA_V,            &
!!  ---  input/output:
             FLUXR_V                                                    &
!! ---  optional outputs:
 !           ,HTRSWB,HTRLWB                                              &
           )
!      do i=1,im
!      write(0,'("zap,i,j,hlw ",2i7,e25.15)')IDEX,jts,HLW(i,1)
!      enddo

      COSZDG(IRANGE,J) = COSZDG_V (1:im)
      CZMEAN(IRANGE,J) = COSZEN_V (1:im)

      DO L=1,LM
        RLWTT(IRANGE,J,L)=HLW(1:im,L)
        RSWTT(IRANGE,J,L)=SWH(1:im,L)
      ENDDO

!=========================================================
! modify this section by using TOPFSW,SFCFSW,TOPFLW,SFCFLW
! instead of TOAUSW,TOADSW,SFCCDSW,TOAULW,SFCUSW
!=========================================================

      ! RLWIN(I,J)=SFCDLW(1)
      ! RSWIN(I,J)=SFCDSW(1)
      ! RSWINC(I,J)=SFCCDSW(1)
      ! RSWOUT(I,J)=RSWIN(I,J)*SFALB(1)
      ! RLWTOA(I,J)=TOAULW(1)
      ! RSWTOA(I,J)=TOAUSW(1)

      ! RSWOUT(I,J)=RSWIN(I,J)*SFALB(1)

      RLWIN(IRANGE,J) =SFCFLW(1:im)%dnfxc
      RSWIN(IRANGE,J) =SFCFSW(1:im)%dnfxc
      RSWOUT(IRANGE,J)=SFCFSW(1:im)%upfxc
      RSWINC(IRANGE,J)=SFCFSW(1:im)%dnfx0

      RLWTOA(IRANGE,J)=TOPFLW(1:im)%upfxc
      RSWTOA(IRANGE,J)=TOPFSW(1:im)%upfxc

!=================================================================
! For non GFDL type cloud (use cloud fields from outputs of GRRAD)
!=================================================================

      IF (ICWP /= -1) THEN
         IF ( LSSAV ) THEN
         !===========================================================
         ! Eliminate cloud fraction form GR1 & RH<95% (20140334, Lin)
         ! EPSQ2=1.e-8 (and not EPSQ=1.e-12) based on multiple tests
         !===========================================================
            DO I=1,IM
              ARG_CW = MAXVAL( CW(IDEX,J,1:LM) )  !- for *total condensate*

              IF (ARG_CW<EPSQ2) THEN  
                 DO L=1,LM
                    CLDCOV_V(I,L) = 0.d0
                 ENDDO
                 DO NC=1,5
                    CLDSA_V(I,NC) = 0.d0
                 ENDDO
              ENDIF
            ENDDO
         !===== end of eliminating extra cloud fraction =====
         !=========================================================

            DO L=1,LM
               CLDFRA(IRANGE,J,L)=CLDCOV_V(1:im,L)
               CSMID(IRANGE,J,L)=CLDCOV_V(1:im,L)
            ENDDO

            CFRACL(IRANGE,J)=CLDSA_V(1:im,1)
            CFRACM(IRANGE,J)=CLDSA_V(1:im,2)
            CFRACH(IRANGE,J)=CLDSA_V(1:im,3)
!
!@@@ To Do:  @@@
!@@@ Add CFRACT array to calculate the total cloud fraction, replace
!    the instantaneous cloud fraction in the post
!
            ACFRST(IRANGE,J)=ACFRST(IRANGE,J) + CLDSA_V(1:im,4)
            NCFRST(IRANGE,J)=NCFRST(IRANGE,J) + 1
!-- Added a time-averaged convective cloud fraction calculation
            WHERE (CU_cloud(1:im)) ACFRCV(IRANGE,J)=ACFRCV(IRANGE,J)+CLDSA_V(1:im,4)
            NCFRCV(IRANGE,J)=NCFRCV(IRANGE,J)+1
         ELSE
            PRINT *, '*** CLDFRA=0, need to set LSSAV=TRUE'
            STOP
         ENDIF
      ENDIF

!-----------------------------------------------------------------------
!
      ENDDO  dayparts
      ENDDO     ! --- END J LOOP for grrad


!
!-----------------------------------------------------------------------
!***  LONGWAVE
!-----------------------------------------------------------------------
!
      IF(MOD(NTIMESTEP,NRADL)==0)THEN
        DO J=JTS,JTE
          DO I=ITS,ITE
!
            TDUM=T(I,J,LM)
            SIGT4(I,J)=STBOLT*TDUM*TDUM*TDUM*TDUM
!
          ENDDO
        ENDDO
      ENDIF 
!-----------------------------------------------------------------------


!
!*** --------------------------------------------------------------------------
!***  DETERMINE THE FRACTIONAL CLOUD COVERAGE FOR HIGH, MID
!***  AND LOW OF CLOUDS FROM THE CLOUD COVERAGE AT EACH LEVEL
!***
!***  NOTE: THIS IS FOR DIAGNOSTICS ONLY!!!
!***
!***
!
!----------------------------------------------------------------------------
      ICWP_Test2: IF (ICWP==-1) THEN   !-- *** Start of old NAM/GFDL cloud ***
!----------------------------------------------------------------------------

       DO J=JTS,JTE
       DO I=ITS,ITE
!!
       DO L=0,LM
         CLDAMT(L)=0.
       ENDDO
!!
!!***  NOW GOES LOW, MIDDLE, HIGH
!!
       DO 480 NLVL=1,3
       CLDMAX=0.
       MALVL=LM
       LLTOP=LM+1-LTOP(NLVL)   !!!!COMES FROM GFDL INIT
!!***
!!***  GO TO THE NEXT CLOUD LAYER IF THE TOP OF THE CLOUD-TYPE IN
!!***  QUESTION IS BELOW GROUND OR IS IN THE LOWEST LAYER ABOVE GROUND.
!!***
       IF(LLTOP.GE.LM)GO TO 480
!!
       IF(NLVL.GT.1)THEN
         LLBOT=LM+1-LTOP(NLVL-1)-1
         LLBOT=MIN(LLBOT,LM1)
       ELSE
         LLBOT=LM1
       ENDIF
!!
       DO 435 L=LLTOP,LLBOT
       CLDAMT(L)=AMAX1(CSMID(I,J,L),CCMID(I,J,L))
       IF(CLDAMT(L).GT.CLDMAX)THEN
         MALVL=L
         CLDMAX=CLDAMT(L)
       ENDIF
   435 CONTINUE
!!*********************************************************************
!! NOW, CALCULATE THE TOTAL CLOUD FRACTION IN THIS PRESSURE DOMAIN
!! USING THE METHOD DEVELOPED BY Y.H., K.A.C. AND A.K. (NOV., 1992).
!! IN THIS METHOD, IT IS ASSUMED THAT SEPERATED CLOUD LAYERS ARE
!! RADOMLY OVERLAPPED AND ADJACENT CLOUD LAYERS ARE MAXIMUM OVERLAPPED.
!! VERTICAL LOCATION OF EACH TYPE OF CLOUD IS DETERMINED BY THE THICKEST
!! CONTINUING CLOUD LAYERS IN THE DOMAIN.
!!*********************************************************************
       CL1=0.0
       CL2=0.0
       KBT1=LLBOT
       KBT2=LLBOT
       KTH1=0
       KTH2=0
!!
       DO 450 LL=LLTOP,LLBOT
       L=LLBOT-LL+LLTOP
       BIT1=.FALSE.
       CR1=CLDAMT(L)
       BITX=(P8W(I,L).GE.PTOPC(NLVL+1)).AND.                           &
      &     (P8W(I,L).LT.PTOPC(NLVL)).AND.                             &
      &     (CLDAMT(L).GT.0.0)
       BIT1=BIT1.OR.BITX
       IF(.NOT.BIT1)GO TO 450
!!***
!!***  BITY=T: FIRST CLOUD LAYER; BITZ=T:CONSECUTIVE CLOUD LAYER
!!***  NOTE:  WE ASSUME THAT THE THICKNESS OF EACH CLOUD LAYER IN THE
!!***         DOMAIN IS LESS THAN 200 MB TO AVOID TOO MUCH COOLING OR
!!***         HEATING. SO WE SET CTHK(NLVL)=200*E2. BUT THIS LIMIT MAY
!!***         WORK WELL FOR CONVECTIVE CLOUDS. MODIFICATION MAY BE
!!***         NEEDED IN THE FUTURE.
!!***
       BITY=BITX.AND.(KTH2.LE.0)
       BITZ=BITX.AND.(KTH2.GT.0)
!!
       IF(BITY)THEN
         KBT2=L
         KTH2=1
       ENDIF
!!
       IF(BITZ)THEN
         KTOP1=KBT2-KTH2+1
         DPCL=P_PHY(I,KBT2)-P_PHY(I,KTOP1)
         IF(DPCL.LT.CTHK(NLVL))THEN
           KTH2=KTH2+1
         ELSE
           KBT2=KBT2-1
         ENDIF
       ENDIF
       IF(BITX)CL2=AMAX1(CL2,CR1)
!!***
!!*** AT THE DOMAIN BOUNDARY OR SEPARATED CLD LAYERS, RANDOM OVERLAP.
!!*** CHOOSE THE THICKEST OR THE LARGEST FRACTION AMT AS THE CLD
!!*** LAYER IN THAT DOMAIN.
!!***
       BIT2=.FALSE.
       BITY=BITX.AND.(CLDAMT(L-1).LE.0.0.OR. &
            P8W(I,L-1).LT.PTOPC(NLVL+1))
       BITZ=BITY.AND.CL1.GT.0.0
       BITW=BITY.AND.CL1.LE.0.0
       BIT2=BIT2.OR.BITY
       IF(.NOT.BIT2)GO TO 450
!!
!!
       IF(BITZ)THEN
         KBT1=INT((CL1*KBT1+CL2*KBT2)/(CL1+CL2))
         KTH1=INT((CL1*KTH1+CL2*KTH2)/(CL1+CL2))+1
         CL1=CL1+CL2-CL1*CL2
       ENDIF
!!
       IF(BITW)THEN
         KBT1=KBT2
         KTH1=KTH2
         CL1=CL2
       ENDIF
!!
       IF(BITY)THEN
         KBT2=LLBOT
         KTH2=0
         CL2=0.0
       ENDIF
  450 CONTINUE
!
        CLDCFR(I,J,NLVL)=AMIN1(1.0,CL1)
        MTOP(I,J,NLVL)=MIN(KBT1,KBT1-KTH1+1)
        MBOT(I,J,NLVL)=KBT1

  480 CONTINUE

      ENDDO ! End DO I=ITS,ITE
      ENDDO ! End DO J=ITS,JTE

!!
      DO J=JTS,JTE
      DO I=ITS,ITE

        CFRACL(I,J)=CLDCFR(I,J,1)
        CFRACM(I,J)=CLDCFR(I,J,2)
        CFRACH(I,J)=CLDCFR(I,J,3)

        IF(CNCLD)THEN
          CFSmax=0.   !-- Maximum cloud fraction (stratiform component)
          CFCmax=0.   !-- Maximum cloud fraction (convective component)
          DO L=1,LM
            CFSmax=MAX(CFSmax, CSMID(I,J,L) )
            CFCmax=MAX(CFCmax, CCMID(I,J,L) )
          ENDDO
          ACFRST(I,J)=ACFRST(I,J)+CFSmax
          NCFRST(I,J)=NCFRST(I,J)+1
          ACFRCV(I,J)=ACFRCV(I,J)+CFCmax
          NCFRCV(I,J)=NCFRCV(I,J)+1
        ELSE
  !--- Count only locations with grid-scale cloudiness, ignore convective clouds
  !    (option not used, but if so set to the total cloud fraction)
          CFRAVG=1.-(1.-CFRACL(I,J))*(1.-CFRACM(I,J))*(1.-CFRACH(I,J))
          ACFRST(I,J)=ACFRST(I,J)+CFRAVG
          NCFRST(I,J)=NCFRST(I,J)+1
        ENDIF

      ENDDO  !  DO I=ITS,ITE
      ENDDO  !  DO J=JTS,JTE

!-----------------------------------------------------------------------
      ENDIF  ICWP_Test2   !*** End of Old NAM/GFDL cloud ***
!-----------------------------------------------------------------------


      END SUBROUTINE RRTM


!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
      SUBROUTINE RRTM_INIT(EMISS,SFULL,SHALF,PPTOP,                     &
     &                     JULYR,MONTH,IDAY,GMT,                        &
     &                     CO2TF,                                       &
     &                     IDS, IDE, JDS, JDE, KDS, KDE,                &
     &                     IMS, IME, JMS, JME, KMS, KME,                &
     &                     ITS, ITE, JTS, JTE, KTS, KTE              )
!-----------------------------------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------------------------------
      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,KDS,KDE                     &
     &                     ,IMS,IME,JMS,JME,KMS,KME                     &
     &                     ,ITS,ITE,JTS,JTE,KTS,KTE
      INTEGER,INTENT(IN) :: JULYR,MONTH,IDAY,CO2TF
      REAL,INTENT(IN) :: GMT,PPTOP
      REAL,DIMENSION(KMS:KME),INTENT(IN) :: SFULL, SHALF
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT) :: EMISS
!
      INTEGER :: I,IHRST,J,N
      REAL :: PCLD,XSD,SQR2PI
      REAL :: SSLP=1013.25
      REAL, PARAMETER :: PTOP_HI=150.,PTOP_MID=350.,PTOP_LO=642.,       &
     &                   PLBTM=105000.
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!***  INITIALIZE DIAGNOSTIC LOW,MIDDLE,HIGH CLOUD LAYER PRESSURE LIMITS.
!
      LTOP(1)=0
      LTOP(2)=0
      LTOP(3)=0
!
      DO N=1,KTE
        PCLD=(SSLP-PPTOP*10.)*SHALF(N)+PPTOP*10.
        IF(PCLD>=PTOP_LO)LTOP(1)=N
        IF(PCLD>=PTOP_MID)LTOP(2)=N
        IF(PCLD>=PTOP_HI)LTOP(3)=N
!       PRINT *,N,PCLD,SHALF(N),PSTAR,PPTOP
      ENDDO
!***
!***  ASSIGN THE PRESSURES FOR CLOUD DOMAIN BOUNDARIES
!***
      PTOPC(1)=PLBTM
      PTOPC(2)=PTOP_LO*100.
      PTOPC(3)=PTOP_MID*100.
      PTOPC(4)=PTOP_HI*100.
!
!***  FOR NOW, GFDL RADIATION ASSUMES EMISSIVITY = 1.0
!
      DO J=JTS,JTE
      DO I=ITS,ITE
        EMISS(I,J) = 1.0
      ENDDO
      ENDDO
!
!---  Calculate the area under the Gaussian curve at the start of the
!---  model run and build the look up table AXSD
!
      SQR2PI=SQRT(2.*PI)
      RSQR=1./SQR2PI
      DO I=1,NXSD
        XSD=REAL(I)*DXSD
        AXSD(I)=GAUSIN(XSD)
      ENDDO
!
!-----------------------------------------------------------------------
      END SUBROUTINE RRTM_INIT
!-----------------------------------------------------------------------
!----------------------------------------------------------------------
!
      REAL FUNCTION GAUSIN(xsd)
      REAL, PARAMETER :: crit=1.e-3
      REAL A1,A2,RN,B1,B2,B3,SUM,xsd
!
!  This function calculate area under the Gaussian curve between mean
!  and xsd # of standard deviation (03/22/2004  Hsin-mu Lin)
!
      a1=xsd*RSQR
      a2=exp(-0.5*xsd**2)
      rn=1.
      b1=1.
      b2=1.
      b3=1.
      sum=1.
      do while (b2 .gt. crit)
         rn=rn+1.
         b2=xsd**2/(2.*rn-1.)
         b3=b1*b2
         sum=sum+b3
         b1=b3
      enddo
      GAUSIN=a1*a2*sum
      RETURN
      END FUNCTION GAUSIN
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      END MODULE MODULE_RA_RRTM
!
!-----------------------------------------------------------------------
