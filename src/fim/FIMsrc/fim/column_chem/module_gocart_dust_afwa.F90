MODULE GOCART_DUST_AFWA
!
!  this module developed by Sandra Jones (AFWA and AER) and Glenn Creighton (AFWA)
!  for serious questions:q

!
!  this module developed by Sandra Jones (AFWA and AER) 
!  and Glenn Creighton (AFWA). For serious questions contact
!


  USE module_data_gocart_dust
  USE namelist_soilveg
  USE module_initial_chem_namelists

CONTAINS
  subroutine gocart_dust_afwa_driver(ktau,dt,alt,t_phy,moist,u_phy,        &
         v_phy,chem,rho_phy,dz8w,smois,u10,v10,p8w,erod,                   &
         ivgtyp,isltyp,vegfra,xland,xlat,xlong,gsw,area,g,emis_dust,       &
         dustin,ust,znt,clay,sand,alpha,gamma,                             &
         num_emis_dust,num_moist,num_chem,num_soil_layers,                 &
         ids,ide, jds,jde, kds,kde,                                        &
         ims,ime, jms,jme, kms,kme,                                        &
         its,ite, jts,jte, kts,kte                                         )
  IMPLICIT NONE

   INTEGER,      INTENT(IN   ) :: ktau,                            &
                                  ids,ide, jds,jde, kds,kde,               &
                                  ims,ime, jms,jme, kms,kme,               &
                                  its,ite, jts,jte, kts,kte,               &
                   num_emis_dust,num_moist,num_chem,num_soil_layers
   INTEGER,DIMENSION( ims:ime , jms:jme )                  ,               &
          INTENT(IN   ) ::                                                 &
                                                     ivgtyp,               &
                                                     isltyp
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_moist ),                &
         INTENT(IN ) ::                              moist
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),                 &
         INTENT(INOUT ) ::                           chem
   REAL, DIMENSION( ims:ime, 1, jms:jme,num_emis_dust),OPTIONAL,           &
         INTENT(INOUT ) ::                                                 &
         emis_dust
   REAL, DIMENSION( ims:ime, num_soil_layers, jms:jme ) ,     &
         INTENT(INOUT) ::                            smois
   REAL, DIMENSION( ims:ime , jms:jme, ndcls )             ,               &
         INTENT(IN   ) ::                            erod
!  REAL,  DIMENSION( ims:ime , jms:jme, 5 )                   ,               &
!         INTENT(INout   ) ::    dustin
   REAL, DIMENSION( ims:ime , jms:jme )                    ,               &
         INTENT(IN   ) ::                                                  &
                                                     u10,                  &
                                                     v10,                  &
                                                     gsw,                  &
                                                     vegfra,               &
                                                     xland,                &
                                                     xlat,                 &
                                                     xlong,area,           &
                                                     ust,                  &
                                                     znt,                  &
                                                     clay,                 &
                                                     sand,dustin
   REAL, DIMENSION( ims:ime , kms:kme , jms:jme ),                         &
         INTENT(IN   ) ::                                                  &
                                                     alt,                  &
                                                     t_phy,                &
                                                     dz8w,p8w,             &
                                                     u_phy,v_phy,rho_phy
  REAL, INTENT(IN   ) :: dt,g

! Local variables

  integer :: nmx,smx,i,j,k,imx,jmx,lmx
  integer,dimension (1,1) :: ilwi
  real*8, DIMENSION (1,1) :: erodtot
  REAL*8, DIMENSION (1,1) :: gravsm
  REAL*8, DIMENSION (1,1) :: drylimit
  real*8, DIMENSION (5)   :: tc,bems
!  real*8, dimension (1,1) :: w10m
  real*8, dimension (1,1) :: airden,airmas,ustar
  real*8, dimension (1) :: dxy
  real*8, dimension (3) :: massfrac
  real*8 :: conver,converi
  real, INTENT(IN   ) :: alpha, gamma

  conver=1.e-9
  converi=1.e9

! Number of dust bins

  imx=1
  jmx=1
  lmx=1
  nmx=ndust
  smx=nsalt

  k=kts
  do j=jts,jte
  do i=its,ite

! Don't do dust over water!!!

    ilwi(1,1)=0
    if(xland(i,j).lt.1.5)then
      ilwi(1,1)=1

! Total concentration at lowest model level. This is still hardcoded for 5 bins.

!    if(config_flags%chem_opt == 2 .or. config_flags%chem_opt == 11 ) then
!       tc(:)=1.e-16*conver
!    else
        tc(1)=chem(i,kts,j,p_dust_1)*conver
        tc(2)=chem(i,kts,j,p_dust_2)*conver
        tc(3)=chem(i,kts,j,p_dust_3)*conver
        tc(4)=chem(i,kts,j,p_dust_4)*conver
        tc(5)=chem(i,kts,j,p_dust_5)*conver
!    endif

!     tc(1)=chem(i,kts,j,p_dust_1)*conver
!     tc(2)=chem(i,kts,j,p_dust_2)*conver
!     tc(3)=chem(i,kts,j,p_dust_3)*conver
!     tc(4)=chem(i,kts,j,p_dust_4)*conver
!     tc(5)=chem(i,kts,j,p_dust_5)*conver

! Air mass and density at lowest model level.

      airmas(1,1)=-(p8w(i,kts+1,j)-p8w(i,kts,j))*area(i,j)/g
      airden(1,1)=rho_phy(i,kts,j)
      ustar(1,1)=ust(i,j)
      dxy(1)=area(i,j)
 
! Total erodibility.

      erodtot(1,1)=SUM(erod(i,j,:))

! Mass fractions of clay, silt, and sand.

      massfrac(1)=clay(i,j)
      massfrac(2)=1-(clay(i,j)+sand(i,j))
      massfrac(3)=sand(i,j)

! Don't allow roughness lengths greater than 20 cm to be lofted.
! This kludge accounts for land use types like urban areas and
! forests which would otherwise show up as high dust emitters.
! This is a placeholder for a more widely accepted kludge
! factor in the literature, which reduces lofting for rough areas.
! Forthcoming...

      IF (znt(i,j) .gt. 0.2) then
        ilwi(1,1)=0
      ENDIF

! Do not allow areas with bedrock, lava, or land-ice to loft

      IF (isltyp(i,j) .eq. 15 .or. isltyp(i,j) .eq. 16. .or. &
          isltyp(i,j) .eq. 18) then
        ilwi(1,1)=0
      ENDIF
      IF (isltyp(i,j) .eq. 0)then
            ilwi(1,1)=0
      endif
      if(ilwi(1,1) == 0 ) cycle

! Calculate gravimetric soil moisture and drylimit.

!      gravsm(1,1)=100*smois(i,1,j)/((1.-maxsmc(isltyp(i,j)))*(2.65*(1-clay(i,j))+2.50*clay(i,j)))
      gravsm(1,1)=100.*smois(i,1,j)/((1.-maxsmc(isltyp(i,j)))*(2.65*(1.-clay(i,j))+2.50*clay(i,j)))
      drylimit(1,1)=14.0*clay(i,j)*clay(i,j)+17.0*clay(i,j)
!     write(0,*) "gravsm(",i,",",j,")=",gravsm(1,1)," drylimit=",drylimit(1)
 
! Call dust emission routine.
! print *, "i,j=",i,j 
! print *, "ustar before call=",ustar(1,1)
      call source_dust(imx, jmx, lmx, nmx, smx, dt, tc, ustar, massfrac, &
                       erodtot, ilwi, dxy, gravsm, airden, airmas, &
                       bems, g, drylimit, alpha, gamma)

!     write(0,*)tc(1)
!     write(0,*)tc(2)
!     write(0,*)tc(3)
!     write(0,*)tc(4)
!     write(0,*)tc(5)
!    if(config_flags%chem_opt == 2 .or. config_flags%chem_opt == 11 ) then
!     dustin(i,j,1:5)=tc(1:5)*converi
!    else
     chem(i,kts,j,p_dust_1)=tc(1)*converi
     chem(i,kts,j,p_dust_2)=tc(2)*converi
     chem(i,kts,j,p_dust_3)=tc(3)*converi
     chem(i,kts,j,p_dust_4)=tc(4)*converi
     chem(i,kts,j,p_dust_5)=tc(5)*converi
!    endif

!     chem(i,kts,j,p_dust_1)=tc(1)*converi
!     chem(i,kts,j,p_dust_2)=tc(2)*converi
!     chem(i,kts,j,p_dust_3)=tc(3)*converi
!     chem(i,kts,j,p_dust_4)=tc(4)*converi
!     chem(i,kts,j,p_dust_5)=tc(5)*converi

! For output diagnostics

      emis_dust(i,1,j,p_edust1)=bems(1)
      emis_dust(i,1,j,p_edust2)=bems(2)
      emis_dust(i,1,j,p_edust3)=bems(3)
      emis_dust(i,1,j,p_edust4)=bems(4)
      emis_dust(i,1,j,p_edust5)=bems(5)
    endif
  enddo
  enddo
!

end subroutine gocart_dust_afwa_driver

  
  SUBROUTINE source_dust(imx, jmx, lmx, nmx, smx, dt1, tc, ustar, massfrac,&
                         erod, ilwi, dxy, gravsm, airden, airmas, &
                         bems, g0, drylimit, alpha, gamma)

! ****************************************************************************
! *  Evaluate the source of each dust particles size bin by soil emission  
! *
! *  Input:
! *         EROD      Fraction of erodible grid cell                (-)
! *         ILWI      Land/water flag                               (-)
! *         GRAVSM    Gravimetric soil moisture                     (g/g)
! *         DRYLIMIT  Upper GRAVSM limit for air-dry soil           (g/g)
! *         ALPHA     Constant to fudge the total emission of dust  (1/m)
! *         DXY       Surface of each grid cell                     (m2)
! *         AIRMAS    Mass of air for each grid box                 (kg)
! *         AIRDEN    Density of air for each grid box              (kg/m3)
! *         USTAR     Friction velocity                             (m/s)
! *         DT1       Time step                                     (s)
! *         NMX       Number of dust bins                           (-)
! *         SMX       Number of saltation bins                      (-)
! *         IMX       Number of I points                            (-)
! *         JMX       Number of J points                            (-)
! *         LMX       Number of L points                            (-)
! *
! *  Data:
! *         MASSFRAC  Fraction of mass in each of 3 soil classes    (-)
! *         SPOINT    Pointer to 3 soil classes                     (-)
! *         DEN_DUST  Dust density                                  (kg/m3)
! *         DEN_SALT  Saltation particle density                    (kg/m3)
! *         REFF_SALT Reference saltation particle diameter         (m)
! *         REFF_DUST Reference dust particle diameter              (m)
! *         LO_DUST   Lower diameter limits for dust bins           (m)
! *         UP_DUST   Upper diameter limits for dust bins           (m)
! *         FRAC_SALT Soil class mass fraction for saltation bins   (-)
! *
! *  Parameters:
! *         CMB       Constant of proportionality                   (-)
! *         MMD_DUST  Mass median diameter of dust                  (m)
! *         GSD_DUST  Geometric standard deviation of dust          (-)
! *         LAMBDA    Side crack propogation length                 (m)
! *         CV        Normalization constant                        (-)
! *         G0        Gravitational acceleration                    (m/s2)
! *         G         Gravitational acceleration in cgs             (cm/s2)
! *      
! *  Working:
! *         U_TS0     "Dry" threshold friction velocity             (m/s)
! *         U_TS      Moisture-adjusted threshold friction velocity (m/s)
! *         RHOA      Density of air in cgs                         (g/cm3)
! *         DEN       Dust density in cgs                           (g/cm3)
! *         DIAM      Dust diameter in cgs                          (cm)
! *         DMASS     Saltation mass distribution                   (-)
! *         DSURFACE  Saltation surface area per unit mass          (m2/kg)
! *         DS_REL    Saltation surface area distribution           (-)
! *         SALT      Saltation flux                                (kg/m/s)
! *         DLNDP     Dust bin width                                (-)
! *         EMIT      Total vertical mass flux                      (kg/m2/s)
! *         EMIT_VOL  Total vertical volume flux                    (m/s)
! *         DSRC      Mass of emitted dust               (kg/timestep/cell)
! *
! *  Output:
! *         TC        Total concentration of dust        (kg/kg/timestep/cell)
! *         BEMS      Source of each dust type           (kg/timestep/cell) 
! *
! ****************************************************************************

  INTEGER, INTENT(IN)   :: nmx,imx,jmx,lmx,smx
  INTEGER, INTENT(IN)   :: ilwi(imx,jmx)
  REAL*8, INTENT(IN)    :: erod(imx,jmx)
  REAL*8, INTENT(IN)    :: ustar(imx,jmx)
  REAL*8, INTENT(IN)    :: gravsm(imx,jmx)
  REAL*8, INTENT(IN)    :: drylimit(imx,jmx) 
  REAL*8, INTENT(IN)    :: dxy(jmx)
  REAL*8, INTENT(IN)    :: airden(imx,jmx,lmx), airmas(imx,jmx,lmx)
  REAL*8, INTENT(INOUT) :: tc(imx,jmx,lmx,nmx)
  REAL*8, INTENT(OUT)   :: bems(imx,jmx,nmx) 
  REAL, INTENT(IN)    :: g0,dt1

  REAL*8    :: den(smx), diam(smx)
  REAL*8    :: dvol(nmx), distr_dust(nmx), dlndp(nmx)
  REAL*8    :: dsurface(smx), ds_rel(smx)
  REAL*8    :: massfrac(3)
  REAL*8    :: u_ts0, u_ts, dsrc, srce, dmass, dvol_tot
  REAL*8    :: emit, emit_vol
  REAL      :: rhoa, g
  INTEGER   :: i, j, m, s

!! Sandblasting mass efficiency, aka "fudge factor" (based on Tegen et al, 
!! 2006 and Hemold et al, 2007)
!
!  REAL, PARAMETER :: alpha=1.8E-8  ! (m^-1)

! Global tuning constant, alpha.  Sandblasting mass efficiency, beta.
! Beta maxes out for clay fractions above 0.2 = betamax.

  REAL, INTENT(IN)  :: alpha
  REAL, PARAMETER :: betamax=5.25E-4
  REAL*8 :: beta

! Experimental optional exponential tuning constant for erodibility.
! 0 < gamma < 1 -> more relative impact by low erodibility regions.
  
  REAL, INTENT(IN) :: gamma

! Constant of proportionality from Marticorena et al, 1997 (unitless)
! Arguably more ~consistent~ fudge than alpha, which has many walnuts
! sprinkled throughout the literature. - GC

  REAL, PARAMETER :: cmb=1.0    
! REAL, PARAMETER :: cmb=2.61   ! from White,1979

! Parameters used in Kok distribution function. Advise not to play with 
! these without the expressed written consent of someone who knows what
! they're doing. - GC

  REAL, PARAMETER :: mmd_dust=3.4D-6  ! median mass diameter (m)
  REAL, PARAMETER :: gsd_dust=3.0     ! geom. std deviation
  REAL, PARAMETER :: lambda=12.0D-6   ! crack propogation length (m)
  REAL, PARAMETER :: cv=12.62D-6      ! normalization constant

! Calculate saltation surface area distribution from sand, silt, and clay
! mass fractions and saltation bin fraction. This will later become a 
! modifier to the total saltation flux.  The reasoning here is that the 
! size and availability of saltators affects saltation efficiency. Based
! on Eqn. (32) in Marticorena & Bergametti, 1995 (hereon, MB95).

  DO n=1,smx
    dmass=massfrac(spoint(n))*frac_salt(n)
    dsurface(n)=0.75*dmass/(den_salt(n)*reff_salt(n))  
  ENDDO
  
! The following equation yields relative surface area fraction.  It will only
! work if you are representing the "full range" of all three soil classes.
! For this reason alone, we have incorporated particle sizes that encompass
! the clay class, to account for the its relative area over the basal
! surface, even though these smaller bins would be unlikely to play any large
! role in the actual saltation process. - GC

  stotal=SUM(dsurface(:))
  DO n=1,smx
    ds_rel(n)=dsurface(n)/stotal
  ENDDO

! Calculate total dust emission due to saltation of sand sized particles.
! Begin by calculating DRY threshold friction velocity (u_ts0).  Next adjust
! u_ts0 for moisture to get threshold friction velocity (u_ts). Then
! calculate saltation flux (salt) where ustar has exceeded u_ts.  Finally, 
! calculate total dust emission (tot_emit), taking into account erodibility. 

 g = g0*1.0E2
 emit=0.0

 DO n = 1, smx
   den(n) = den_salt(n)*1.0D-3         ! (g cm^-3)
   diam(n) = 2.0*reff_salt(n)*1.0D2    ! (cm)
   DO i = 1,imx
     DO j = 1,jmx
       rhoa = airden(i,j,1)*1.0D-3       ! (g cm^-3)

     ! Threshold friction velocity as a function of the dust density and
     ! diameter from Bagnold (1941) (m s^-1).

       u_ts0 = 0.13*1.0D-2*SQRT(den(n)*g*diam(n)/rhoa)* &
               SQRT(1.0+0.006/den(n)/g/(diam(n))**2.5)/ &
               SQRT(1.928*(1331.0*(diam(n))**1.56+0.38)**0.092-1.0) 

     ! Friction velocity threshold correction function based on physical
     ! properties related to moisture tension. Soil moisture greater than
     ! dry limit serves to increase threshold friction velocity (making
     ! it more difficult to loft dust). When soil moisture has not reached
     ! dry limit, treat as dry (no correction to threshold friction
     ! velocity). GC

       IF (gravsm(i,j) > drylimit(i,j)) THEN
         u_ts = MAX(0.0D+0,u_ts0*(sqrt(1.0+1.21*(gravsm(i,j)-drylimit(i,j))**0.68)))
       ELSE
         u_ts = u_ts0
       END IF 

     ! Saltation flux from Marticorena & Bergametti 1995 (MB95). ds_rel is
     ! the relative surface area distribution

       IF (ustar(i,j) .gt. u_ts .and. erod(i,j) .gt. 0.0 .and. ilwi(i,j) == 1) THEN
         salt = cmb*ds_rel(n)*(airden(i,j,1)/g0)*(ustar(i,j)**3)* &
                (1. + u_ts/ustar(i,j))*(1. - (u_ts**2)/(ustar(i,j)**2))  ! (kg m^-1 s^-1)
       ELSE 
         salt = 0.0
       ENDIF

     ! Calculate total vertical mass flux (note beta has units of m^-1)
     ! Beta acts to tone down dust in areas with so few dust-sized particles that the
     ! lofting efficiency decreases.  Otherwise, super sandy zones would be huge dust
     ! producers, which is generally not the case.  Equation derived from wind-tunnel 
     ! experiments (see MB95).

       beta=10**(13.6*massfrac(1)-6.0)  ! (unitless)
       if (beta .gt. betamax) then
         beta=betamax
       endif
      ! emit=emit+salt*erod(i,j)*alpha*beta    ! (kg m^-2 s^-1)
       emit=emit+salt*(erod(i,j)**gamma)*alpha*beta    ! (kg m^-2 s^-1)
     END DO
   END DO
 END DO

! Now that we have the total dust emission, distribute into dust bins using 
! lognormal distribution (Dr. Jasper Kok, in press), and
! calculate total mass emitted over the grid box over the timestep. 
!
! In calculating the Kok distribution, we assume upper and lower limits to each bin.
! For reff_dust=(/0.73D-6,1.4D-6,2.4D-6,4.5D-6,8.0D-6/) (default),
! lower limits were ASSUMED at lo_dust=(/0.1D-6,1.0D-6,1.8D-6,3.0D-6,6.0D-6/)
! upper limits were ASSUMED at up_dust=(/1.0D-6,1.8D-6,3.0D-6,6.0D-6,10.0D-6/)
! These may be changed within module_data_gocart_dust.F, but make sure it is
! consistent with reff_dust values.  These values were taken from the original
! GOCART bin configuration. We use them here to calculate dust bin width, dlndp.
! dVol is the volume distribution. You know...if you were wondering. GC

 dvol_tot=0.
 DO n=1,nmx
   dlndp(n)=LOG(up_dust(n)/lo_dust(n))
   dvol(n)=(2.0*reff_dust(n)/cv)*(1.+ERF(LOG(2.0*reff_dust(n)/mmd_dust)/(SQRT(2.)*LOG(gsd_dust))))*&
         EXP(-(2.0*reff_dust(n)/lambda)**3.0)*dlndp(n)
   dvol_tot=dvol_tot+dvol(n)
  ! Convert mass flux to volume flux
   emit_vol=emit/den_dust(n) ! (m s^-1)
 END DO 
 DO n=1,nmx
   distr_dust(n)=dvol(n)/dvol_tot
   !print *,"distr_dust(",n,")=",distr_dust(n)
 END DO

! Now distribute total vertical emission into dust bins and update concentration.

 DO n=1,nmx 
   DO i=1,imx
     DO j=1,jmx
      ! Calculate total mass emitted
        dsrc = emit_vol*den_dust(n)*distr_dust(n)*dxy(j)*dt1  ! (kg)
        IF (dsrc < 0.0) dsrc = 0.0

      ! Update dust mixing ratio at first model level.
        tc(i,j,1,n) = tc(i,j,1,n) + dsrc / airmas(i,j,1) ! (kg/kg)
     !   bems(i,j,n) = dsrc  ! diagnostic
        bems(i,j,n) = 1000.*dsrc/(dxy(j)*dt1) ! diagnostic (g/m2/s)
     END DO
   END DO
 END DO

END SUBROUTINE source_dust


END MODULE GOCART_DUST_AFWA
