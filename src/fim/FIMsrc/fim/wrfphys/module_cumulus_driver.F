!WRF:MEDIATION_LAYER:PHYSICS
!

MODULE module_cumulus_driver
   USE module_constants,only: g=>grvity,cp,r_d=>rd,r_v=>rv
   REAL    , PARAMETER :: XLV0         = 3.15E6
   REAL    , PARAMETER :: XLV1         = 2370.
   REAL    , PARAMETER :: XLS0         = 2.905E6
   REAL    , PARAMETER :: XLS1         = 259.532

   REAL    , PARAMETER :: XLS          = 2.85E6
   REAL    , PARAMETER :: XLV          = 2.5E6
   REAL    , PARAMETER :: XLF          = 3.50E5

   REAL    , PARAMETER ::  SVP1=0.6112
   REAL    , PARAMETER ::  SVP2=17.67
   REAL    , PARAMETER ::  SVP3=29.65
   REAL    , PARAMETER ::  SVPT0=273.15
   REAL    , PARAMETER ::  EP_1=R_v/R_d-1.
   REAL    , PARAMETER ::  EP_2=R_d/R_v
CONTAINS
   SUBROUTINE cumulus_driver(                                     &
                 ! Order dependent args for domain, mem, and tile dims
                      ids,ide, jds,jde, kds,kde                       &
                     ,ims,ime, jms,jme, kms,kme                       &
                      ,its,ite,jts,jte,kts,kte                        &
                 ! Order independent args (use VAR= in call)
                 ! --Prognostic
                     ,u,v,th,t,w                                      &
                     ,p,pi,rho                                        &
                 ! --Other arguments
                     ,itimestep,dt,dx,cudt,curr_secs,adapt_step_flag  &
                     ,rainc,raincv,pratec                             &
                     ,dz8w,p8w,forcet,forceq                          &
                     ,stepcu,gsw                                &
                     ,xland,cu_act_flag       &
                     ,htop,hbot,kpbl,ht                               &  
                     ,ensdim,maxiens,maxens,maxens2,maxens3           &
                 ! Package selection variable
                     ,cu_physics                                      &
                 ! Optional moisture tracers 
                     ,qv_curr, qc_curr, qr_curr                       &
                     ,qi_curr, qs_curr, qg_curr                       & 
                     ,qv_prev, qc_prev, qr_prev                       & 
                     ,qi_prev, qs_prev, qg_prev                       &
                 ! Optional arguments for GD scheme
                     ,apr_gr,apr_w,apr_mc,apr_st,apr_as,apr_capma     &
                     ,apr_capme,apr_capmi,edt_out,clos_choice         &
                     ,mass_flux,xf_ens,pr_ens,cugd_avedx,imomentum    &
                     ,cugd_tten,cugd_qvten,cugd_qcten                 &
                     ,cugd_ttens,cugd_qvtens                          &
                     ,gd_cloud,gd_cloud2      &
                 ! Optional moisture and other tendencies
                     ,rqvcuten,rqccuten,rqrcuten                      &
                     ,rqicuten,rqscuten,rqgcuten                      &
                     ,rqvblten,rqvften                                &
                     ,rthcuten,rthraten,rthblten,rthften              &
                 ! Optional moisture tracer flags 
                     ,f_qv,f_qc,f_qr                                  &
                     ,f_qi,f_qs,f_qg                                  &
                                                                      )
!----------------------------------------------------------------------
   USE module_initial_chem_namelists , ONLY :                            &
                                           KFSCHEME,BMJSCHEME         &
                                          ,KFETASCHEME,GDSCHEME       &
                                          ,G3SCHEME                   &
                                          ,SASSCHEME

! *** add new modules of schemes here

!  USE module_cu_kf
!  USE module_cu_bmj
!  USE module_dm
!  USE module_domain, ONLY: domain
!  USE module_cu_kfeta
   USE module_cu_gd, ONLY : GRELLDRV
   USE module_cu_g3, ONLY : G3DRV,CONV_GRELL_SPREAD3D
!  USE module_cu_sas

   !  This driver calls subroutines for the cumulus parameterizations.
   !
   !  1. Kain & Fritsch (1993)
   !  2. Betts-Miller-Janjic (Janjic, 1994)
   !
!----------------------------------------------------------------------
   IMPLICIT NONE
!======================================================================
! Grid structure in physics part of WRF
!----------------------------------------------------------------------
! The horizontal velocities used in the physics are unstaggered
! relative to temperature/moisture variables. All predicted
! variables are carried at half levels except w, which is at full
! levels. Some arrays with names (*8w) are at w (full) levels.
!
!----------------------------------------------------------------------
! In WRF, kms (smallest number) is the bottom level and kme (largest
! number) is the top level.  In your scheme, if 1 is at the top level,
! then you have to reverse the order in the k direction.
!
!         kme      -   half level (no data at this level)
!         kme    ----- full level
!         kme-1    -   half level
!         kme-1  ----- full level
!         .
!         .
!         .
!         kms+2    -   half level
!         kms+2  ----- full level
!         kms+1    -   half level
!         kms+1  ----- full level
!         kms      -   half level
!         kms    ----- full level
!
!======================================================================
! Definitions
!-----------
! Rho_d      dry density (kg/m^3)
! Theta_m    moist potential temperature (K)
! Qv         water vapor mixing ratio (kg/kg)
! Qc         cloud water mixing ratio (kg/kg)
! Qr         rain water mixing ratio (kg/kg)
! Qi         cloud ice mixing ratio (kg/kg)
! Qs         snow mixing ratio (kg/kg)
!-----------------------------------------------------------------
!-- DT            time step (second)
!-- CUDT          cumulus time step (minute)
!-- curr_secs     current forecast time (seconds)
!-- itimestep     number of time step (integer)   
!-- DX            horizontal space interval (m)
!-- rr            dry air density (kg/m^3)
!
!-- RTHCUTEN      Theta tendency due to 
!                 cumulus scheme precipitation (K/s)
!-- RQVCUTEN      Qv tendency due to 
!                 cumulus scheme precipitation (kg/kg/s)
!-- RQRCUTEN      Qr tendency due to 
!                 cumulus scheme precipitation (kg/kg/s)
!-- RQCCUTEN      Qc tendency due to 
!                 cumulus scheme precipitation (kg/kg/s)
!-- RQSCUTEN      Qs tendency due to 
!                 cumulus scheme precipitation (kg/kg/s)
!-- RQICUTEN      Qi tendency due to 
!                 cumulus scheme precipitation (kg/kg/s)
!
!-- RAINC         accumulated total cumulus scheme precipitation (mm)
!-- RAINCV        cumulus scheme precipitation (mm)
!-- PRATEC        precipitiation rate from cumulus scheme (mm/s)
!                 time in KF cumulus scheme (integer)
!-- u_phy         u-velocity interpolated to theta points (m/s)
!-- v_phy         v-velocity interpolated to theta points (m/s)
!-- th_phy        potential temperature (K)
!-- t_phy         temperature (K)
!-- w             vertical velocity (m/s)
!-- moist         moisture array (4D - last index is species) (kg/kg)
!-- dz8w          dz between full levels (m)
!-- p8w           pressure at full levels (Pa)  
!-- p_phy         pressure (Pa)
!-- pi_phy        exner function (dimensionless)
!                 points (dimensionless)
!-- RTHRATEN      radiative temp forcing for Grell-Devenyi scheme
!-- RTHBLTEN      PBL temp forcing for Grell-Devenyi scheme
!-- RQVBLTEN      PBL moisture forcing for Grell-Devenyi scheme
!-- RTHFTEN
!-- RQVFTEN
!-- MASS_FLUX
!-- XF_ENS
!-- PR_ENS
!-- CU_ACT_FLAG
!-- rho           density (kg/m^3)
!-- STEPCU        # of fundamental timesteps between convection calls
!-- XLAND         land-sea mask (1.0 for land; 2.0 for water)
!-- XLV0          latent heat of vaporization constant 
!                 used in temperature dependent formula (J/kg)
!-- XLV1          latent heat of vaporization constant 
!                 used in temperature dependent formula (J/kg/K)
!-- XLS0          latent heat of sublimation constant 
!                 used in temperature dependent formula (J/kg)
!-- XLS1          latent heat of sublimation constant
!                 used in temperature dependent formula (J/kg/K)
!-- R_d           gas constant for dry air ( 287. J/kg/K)
!-- R_v           gas constant for water vapor (461 J/k/kg)
!-- Cp            specific heat at constant pressure (1004 J/k/kg)
!-- rvovrd        R_v divided by R_d (dimensionless)
!-- G             acceleration due to gravity (m/s^2)
!-- EP_1          constant for virtual temperature 
!                 (R_v/R_d - 1) (dimensionless)
!-- pi_phy        the exner function, (p/p0)**(R/Cp) (none unit)
!-- ids           start index for i in domain
!-- ide           end index for i in domain
!-- jds           start index for j in domain
!-- jde           end index for j in domain
!-- kds           start index for k in domain
!-- kde           end index for k in domain
!-- ims           start index for i in memory
!-- ime           end index for i in memory
!-- jms           start index for j in memory
!-- jme           end index for j in memory
!-- kms           start index for k in memory
!-- kme           end index for k in memory
!-- i_start       start indices for i in tile
!-- i_end         end indices for i in tile
!-- j_start       start indices for j in tile
!-- j_end         end indices for j in tile
!-- kts           start index for k in tile
!-- kte           end index for k in tile
!-- num_tiles     number of tiles
!-- HBOT          index of lowest model layer with convection
!-- HTOP          index of highest model layer with convection
!-- LBOT          index of lowest model layer with convection
!-- LTOP          index of highest model layer with convection
!-- KPBL          layer index of the PBL
!-- periodic_x    T/F this is using periodic lateral boundaries in the X direction
!-- periodic_y    T/F this is using periodic lateral boundaries in the Y-direction
!
!======================================================================

   INTEGER,      INTENT(IN   )    ::                             &
                                      ids,ide, jds,jde, kds,kde, &
                                      ims,ime, jms,jme, kms,kme, &
                                      its,ite, jts,jte, kts,kte, &        
                                      itimestep

   INTEGER,      INTENT(IN   )    ::                             &
                  ensdim,maxiens,maxens,maxens2,maxens3

   INTEGER, OPTIONAL,     INTENT(IN   )    ::                    &
                   cugd_avedx,clos_choice

   INTEGER,      INTENT(IN   )    ::   cu_physics
   INTEGER,      INTENT(IN   )    ::   STEPCU

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                 &
         INTENT(IN ) ::                                          &
                                                           dz8w  &
                                                      ,     p8w  &
                                                      ,       p  &
                                                      ,      pi  &
                                                      ,       u  &
                                                      ,       v  &
                                                      ,      th  &
                                                      ,       t  &
                                                      ,     rho  &
                                                      ,       w


   REAL, DIMENSION( ims:ime , jms:jme ), INTENT(IN) ::           &
                                                 GSW,HT,XLAND

   REAL, DIMENSION( ims:ime , jms:jme ),                         &
          INTENT(INOUT) ::                                RAINC  &
                                                    ,    RAINCV  &
                                                    ,      HTOP  & 
                                                    ,      HBOT 
 

   REAL, DIMENSION( ims:ime , jms:jme ),INTENT(INOUT),OPTIONAL :: &
        PRATEC
   REAL, DIMENSION( ims:ime , jms:jme ) :: tmppratec
                                                    
   INTEGER, DIMENSION( ims:ime , jms:jme ),                      &
                    INTENT(IN) ::                          KPBL


   LOGICAL, DIMENSION( ims:ime , jms:jme ),                      &
          INTENT(INOUT) :: CU_ACT_FLAG

   REAL,  INTENT(IN   ) :: DT, DX
   INTEGER,      INTENT(IN   ),OPTIONAL    ::                             &
                                      imomentum
   REAL,  INTENT(IN   ),OPTIONAL :: CUDT
   REAL,  INTENT(IN   ),OPTIONAL :: CURR_SECS
   LOGICAL,INTENT(IN   ),OPTIONAL    ::     adapt_step_flag
   REAL   :: cudt_pass, curr_secs_pass
   LOGICAL :: adapt_step_flag_pass,periodic_x,periodic_y

!
! optional arguments
!
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                 &
         OPTIONAL, INTENT(INOUT) ::                              &
                      ! optional moisture tracers
                      ! 2 time levels; if only one then use CURR
                      qv_curr, qc_curr, qr_curr                  &
                     ,qi_curr, qs_curr, qg_curr                  & 
                     ,qv_prev, qc_prev, qr_prev                  & 
                     ,qi_prev, qs_prev, qg_prev                  &
                      ! optional moisture and other tendencies
                     ,rqvcuten,rqccuten,rqrcuten                 &
                     ,rqicuten,rqscuten,rqgcuten                 &
                     ,rqvblten,rqvften                           &
                     ,rthraten,rthblten                          &
                     ,cugd_tten,cugd_qvten,cugd_qcten            &
                     ,cugd_ttens,cugd_qvtens                     &
                                                      ,   forcet &
                                                      ,   forceq &
                     ,rthften,rthcuten

   REAL, DIMENSION( ims:ime , jms:jme ),                         &
                    OPTIONAL,                                    &
                    INTENT(INOUT) ::                             &
                apr_gr,apr_w,apr_mc,apr_st,apr_as,apr_capma      &
               ,apr_capme,apr_capmi,edt_out                      &
                                                    , MASS_FLUX
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                 &
         OPTIONAL, INTENT(INOUT) ::                              &
                  GD_CLOUD,GD_CLOUD2
   REAL, DIMENSION( ims:ime , jms:jme , 1:ensdim ),              &
          OPTIONAL,                                              &
          INTENT(INOUT) ::                       XF_ENS, PR_ENS

!
! Flags relating to the optional tendency arrays declared above
! Models that carry the optional tendencies will provdide the
! optional arguments at compile time; these flags all the model
! to determine at run-time whether a particular tracer is in
! use or not.
!
   LOGICAL, INTENT(IN), OPTIONAL ::                             &
                                                      f_qv      &
                                                     ,f_qc      &
                                                     ,f_qr      &
                                                     ,f_qi      &
                                                     ,f_qs      &
                                                     ,f_qg


! LOCAL  VAR

   INTEGER :: i,j,k,ij

!-----------------------------------------------------------------

   if (.not. PRESENT(CURR_SECS)) then
      curr_secs_pass = -1
   else 
      curr_secs_pass = curr_secs
   endif

   if (.not. PRESENT(CUDT)) then
      cudt_pass = -1
   else
      cudt_pass = cudt
   endif

   if (.not. PRESENT(adapt_step_flag)) then
      adapt_step_flag_pass = .false.
   else
      adapt_step_flag_pass = adapt_step_flag
   endif
   periodic_x=.true.
   periodic_y=.true.

   ! Initialize tmppratec to pratec

   if ( PRESENT ( pratec ) ) then
      tmppratec(:,:) = pratec(:,:)
   else
      tmppratec(:,:) = 0.
   end if


   IF (cu_physics .eq. 0) return
      if(cu_physics .eq. 5 ) then
        do j=jts,min(jte,jde-1)
        do k=kts,kte
        do i=its,min(ite,ide-1)
           RTHFTEN(i,k,j)=(RTHFTEN(i,k,j)+RTHRATEN(i,k,j) &
                               +RTHBLTEN(i,k,j))*pi(i,k,j)
           RQVFTEN(i,k,j)=RQVFTEN(i,k,j)+RQVBLTEN(i,k,j)
       enddo
       enddo
       enddo
!     else if(cu_physics .eq. 3 ) then
!       do j=jts,min(jte,jde-1)
!       do k=kts,kte
!       do i=its,min(ite,ide-1)
!          RTHFTEN(i,k,j)=(RTHFTEN(i,k,j)+RTHRATEN(i,k,j) &
!                              +RTHBLTEN(i,k,j))
!          RQVFTEN(i,k,j)=RQVFTEN(i,k,j)+RQVBLTEN(i,k,j)
!      enddo
!      enddo
!      enddo
      endif

! DON'T JUDGE TIME STEP HERE, SINCE KF NEEDS ACCUMULATED W FIELD.
! DO IT INSIDE THE INDIVIDUAL CUMULUS SCHEME



   cps_select: SELECT CASE(cu_physics)

     CASE (GDSCHEME)
!TBH          print *,'call grelldrv'
!         CALL wrf_debug(100,'in grelldrv')
          CALL GRELLDRV(                                        &
                DT=dt, ITIMESTEP=itimestep, DX=dx               &
               ,U=u,V=v,T=t,W=w ,RHO=rho                        &
               ,P=p,PI=pi ,Q=qv_curr ,RAINCV=raincv             &
               ,DZ8W=dz8w,P8W=p8w,XLV=xlv,CP=cp,G=g,R_V=r_v     &
               ,PRATEC=tmppratec                                &
               ,APR_GR=apr_gr,APR_W=apr_w,APR_MC=apr_mc         &
               ,APR_ST=apr_st,APR_AS=apr_as                     &
               ,APR_CAPMA=apr_capma,APR_CAPME=apr_capme         &
               ,APR_CAPMI=apr_capmi,MASS_FLUX=mass_flux         &
               ,XF_ENS=xf_ens,PR_ENS=pr_ens,HT=ht               &
               ,xland=xland,gsw=gsw                             &
               ,GDC=gd_cloud,GDC2=gd_cloud2 &
               ,ENSDIM=ensdim,MAXIENS=maxiens,MAXENS=maxens     &
               ,MAXENS2=maxens2,MAXENS3=maxens3                 &
               ,STEPCU=STEPCU,htop=htop,hbot=hbot               &
               ,CU_ACT_FLAG=CU_ACT_FLAG     &
               ,IDS=ids,IDE=ide,JDS=jds,JDE=jde,KDS=kds,KDE=kde &
               ,IMS=ims,IME=ime,JMS=jms,JME=jme,KMS=kms,KME=kme &
               ,ITS=its,ITE=ite,JTS=jts,JTE=jte,KTS=kts,KTE=kte &
               ,PERIODIC_X=periodic_x,PERIODIC_Y=periodic_y     &
              ! optionals
               ,RTHCUTEN=RTHCUTEN ,RTHFTEN=RTHFTEN              &
               ,RQICUTEN=RQICUTEN ,RQVFTEN=RQVFTEN              &
               ,RTHRATEN=RTHRATEN,RTHBLTEN=RTHBLTEN             &
               ,RQVCUTEN=RQVCUTEN,RQCCUTEN=RQCCUTEN             &
               ,RQVBLTEN=RQVBLTEN                               &
               ,F_QV=f_qv,F_QC=f_qc,F_QR=f_qr                   &
               ,F_QI=f_qi,F_QS=f_qs                             &
                                                                )
!    CASE (SASSCHEME)
                                                                                                                                           
!         CALL wrf_debug(100,'in cu_sas')
!         CALL CU_SAS(                                          &
!               DT=dt,ITIMESTEP=itimestep,STEPCU=STEPCU         &
!              ,RAINCV=RAINCV,PRATEC=tmpPRATEC,HTOP=HTOP,HBOT=HBOT &
!              ,U3D=u,V3D=v,W=w,T3D=t,PI3D=pi,RHO3D=rho         &
!              ,QV3D=QV_CURR,QC3D=QC_CURR,QI3D=QI_CURR          &
!              ,DZ8W=dz8w,PCPS=p,P8W=p8w,XLAND=XLAND            &
!              ,CU_ACT_FLAG=CU_ACT_FLAG                         &
!              ,CUDT=cudt_pass                                  &
!              ,CURR_SECS=curr_secs_pass                        &
!              ,ADAPT_STEP_FLAG=adapt_step_flag_pass            &
!              ,IDS=ids,IDE=ide,JDS=jds,JDE=jde,KDS=kds,KDE=kde &
!              ,IMS=ims,IME=ime,JMS=jms,JME=jme,KMS=kms,KME=kme &
!              ,ITS=its,ITE=ite,JTS=jts,JTE=jte,KTS=kts,KTE=kte &
              ! optionals
!              ,RTHCUTEN=RTHCUTEN,RQVCUTEN=RQVCUTEN             &
!              ,RQCCUTEN=RQCCUTEN,RQICUTEN=RQICUTEN             &
!              ,F_QV=f_qv,F_QC=f_qc,F_QR=f_qr                   &
!              ,F_QI=f_qi,F_QS=f_qs                             &
!                                                               )
     CASE (G3SCHEME)
!         CALL wrf_debug(100,'in grelldrv')
          CALL G3DRV(                                           &
                DT=dt, ITIMESTEP=itimestep, DX=dx               &
               ,U=u,V=v,T=t,W=w ,RHO=rho                        &
               ,P=p,PI=pi,Q=qv_curr,RAINCV=raincv               &
               ,DZ8W=dz8w ,P8W=p8w,XLV=xlv,CP=cp,G=g,R_V=r_v    &
               ,APR_GR=apr_gr,APR_W=apr_w,APR_MC=apr_mc         &
               ,APR_ST=apr_st,APR_AS=apr_as,PRATEC=tmppratec    &
               ,APR_CAPMA=apr_capma,APR_CAPME=apr_capme         &
               ,APR_CAPMI=apr_capmi,MASS_FLUX=mass_flux         &
               ,XF_ENS=xf_ens,PR_ENS=pr_ens,HT=ht               &
               ,xland=xland,gsw=gsw,edt_out=edt_out             &
               ,GDC=gd_cloud,GDC2=gd_cloud2                     &
               ,cugd_tten=cugd_tten,cugd_qvten=cugd_qvten       &
               ,cugd_ttens=cugd_ttens,cugd_qvtens=cugd_qvtens   &
               ,cugd_qcten=cugd_qcten,cugd_avedx=cugd_avedx     &
               ,imomentum=imomentum                             &
               ,ENSDIM=ensdim,MAXIENS=maxiens,MAXENS=maxens     &
               ,MAXENS2=maxens2,MAXENS3=maxens3,ichoice=clos_choice &
               ,STEPCU=STEPCU,htop=htop,hbot=hbot               &
               ,CU_ACT_FLAG=CU_ACT_FLAG     &
               ,IDS=ids,IDE=ide,JDS=jds,JDE=jde,KDS=kds,KDE=kde &
               ,IMS=ims,IME=ime,JMS=jms,JME=jme,KMS=kms,KME=kme &
               ,ITS=its,ITE=ite,JTS=jts,JTE=jte,KTS=kts,KTE=kte &
               ,PERIODIC_X=periodic_x,PERIODIC_Y=periodic_y     &
              ! optionals
               ,RTHCUTEN=RTHCUTEN ,RTHFTEN=RTHFTEN              &
               ,RQICUTEN=RQICUTEN ,RQVFTEN=RQVFTEN              &
               ,RQVCUTEN=RQVCUTEN,RQCCUTEN=RQCCUTEN             &
               ,F_QV=f_qv,F_QC=f_qc,F_QR=f_qr                   &
               ,F_QI=f_qi,F_QS=f_qs                             &
                                                                )

     CASE DEFAULT 

!        WRITE( wrf_err_message , * ) 'The cumulus option does not exist: cu_physics = ', cu_physics
!        CALL wrf_error_fatal ( wrf_err_message )

   END SELECT cps_select

!   IF(cu_physics .eq. 5 )then
!!#ifdef DM_PARALLEL
!!#     include "HALO_CUP_G3_OUT.inc"
!!#endif
!        call conv_grell_spread3d(rthcuten=rthcuten,rqvcuten=rqvcuten                &
!     &            ,rqccuten=rqccuten,raincv=raincv,cugd_avedx=cugd_avedx            &
!     &            ,cugd_tten=cugd_tten,cugd_qvten=cugd_qvten,rqicuten=rqicuten      &
!     &            ,cugd_ttens=cugd_ttens,cugd_qvtens=cugd_qvtens                    &
!     &            ,cugd_qcten=cugd_qcten,pi_phy=pi,moist_qv=qv_curr                 &
!     &            ,PRATEC=tmppratec,dt=dt,num_tiles=num_tiles                       &
!     &            ,imomentum=imomentum                             &
!     &            ,F_QV=F_QV,F_QC=F_QC,F_QR=F_QR,F_QI=F_QI,F_QS=F_QS                &
!     &            ,ids=IDS,ide=IDE, jds=JDS,jde=JDE, kds=KDS,kde=KDE                &
!     &            ,ims=IMS,ime=IME, jms=JMS,jme=JME, kms=KMS,kme=KME                &
!     &            ,i_start=i_start,i_end=i_end                                      &
!     &            ,j_start=j_start,j_end=j_end                                      &
!     &            ,kts=kts, kte=kte)
!   endif

      !
      ! Copy pratec back to output array, if necessary.
      !
      if (PRESENT(PRATEC)) then
         pratec(:,:) = tmppratec(:,:)
      endif
   END SUBROUTINE cumulus_driver

END MODULE module_cumulus_driver
