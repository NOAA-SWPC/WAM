MODULE GOCART_DUST
  

  USE module_data_gocart_dust
  USE namelist_soilveg
  USE module_initial_chem_namelists
! USE module_initial_chem_namelist_defaults !, only: set_species

CONTAINS
  subroutine gocart_dust_driver(ktau,dt,alt,t_phy,moist,u_phy,  &
         v_phy,chem,rho_phy,dz8w,smois,u10,v10,p8w,erod,                  &
         ivgtyp,isltyp,vegfra,xland,xlat,xlong,gsw,area,g,emis_dust,        &
         dusthelp,num_emis_dust,num_moist,num_chem,num_soil_layers,                 &
         start_month,                                               &
         ids,ide, jds,jde, kds,kde,                                        &
         ims,ime, jms,jme, kms,kme,                                        &
         its,ite, jts,jte, kts,kte                                         )
! USE module_configure
! USE module_state_description
  IMPLICIT NONE
!  TYPE(grid_config_rec_type),  INTENT(IN   )    :: config_flags

   INTEGER,      INTENT(IN   ) :: ktau,start_month,                  &
         num_emis_dust,num_moist,num_chem,num_soil_layers,                 &
                             ids,ide, jds,jde, kds,kde,               &
                                  ims,ime, jms,jme, kms,kme,               &
                                  its,ite, jts,jte, kts,kte
   INTEGER,DIMENSION( ims:ime , jms:jme )                  ,               &
          INTENT(IN   ) ::                                                 &
                                                     ivgtyp,               &
                                                     isltyp
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_moist ),                &
         INTENT(IN ) ::                                   moist
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),                 &
         INTENT(INOUT ) ::                                   chem
   REAL, DIMENSION( ims:ime, 1, jms:jme,num_emis_dust),OPTIONAL,&
         INTENT(INOUT ) ::                                                 &
         emis_dust
  REAL, DIMENSION( ims:ime, num_soil_layers, jms:jme ) ,      &
      INTENT(INOUT) ::                               smois
   REAL,  DIMENSION( ims:ime , jms:jme, 3 )                   ,               &
          INTENT(IN   ) ::    erod
   REAL,  DIMENSION( ims:ime , jms:jme )                   ,               &
          INTENT(IN   ) ::                                                 &
                                                     u10,                  &
                                                     v10,                  &
                                                     gsw,                  &
                                                  vegfra,                  &
                                                     xland,                &
                                                     xlat,                 &
                                                     xlong,area
   REAL,  DIMENSION( ims:ime , jms:jme ),                        &
          INTENT(OUT  ) :: dusthelp
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme ),                        &
          INTENT(IN   ) ::                                                 &
                                                        alt,               &
                                                      t_phy,               &
                                                     dz8w,p8w,             &
                                              u_phy,v_phy,rho_phy
 
  REAL, INTENT(IN   ) :: dt,g
!
! local variables
!
  integer :: nmx,i,j,k,imx,jmx,lmx,ipr
  integer,dimension (1,1) :: ilwi
  real*8, DIMENSION (1,1,3,1) :: erodin
  real*8, DIMENSION (5) :: tc,bems
  real*8, dimension (1,1) :: w10m,gwet,airden,airmas
  real*8, dimension (1) :: dxy
  real*8  tcs,conver,converi
  real dttt
  real*8,parameter::max_default=0.
! write(6,*)'in dust driver ',ktau,dt,start_month
! conver=1.e-9*mwdry
! converi=1.e9/mwdry
  conver=1.e-9
  converi=1.e9
!
! number of dust bins
!
  imx=1
  jmx=1
  lmx=1
  nmx=5 
  k=kts
  if(chem_opt == 304 .or. chem_opt == 316 .or. chem_opt == 317) Then
! print *,'chem_opt = ',chem_opt,'in gocart_dust',p_dust_1,p_dust_2
  dusthelp(:,:)=0.
  do j=jts,jte
  do i=its,ite
!
! 
!
     if(xland(i,j).lt.1.5 .and. xland(i,j).gt.0.5)then
     ilwi(1,1)=1
     tc(1)=chem(i,kts,j,p_dust_1)*conver
     tcs=tc(1)
     tc(2)=1.d-30
     tc(3)=chem(i,kts,j,p_dust_2)*conver
     tc(4)=1.d-30
     tc(5)=1.d-30
     w10m(1,1)=sqrt(u10(i,j)*u10(i,j)+v10(i,j)*v10(i,j))
     airmas(1,1)=-(p8w(i,kts+1,j)-p8w(i,kts,j))*area(i,j)/g
!    if(j.eq.681)then
!     write(6,*)chem(i,kts,j,p_dust_1),chem(i,kts,j,p_dust_2) 
!     write(6,*)tc(1),tc(3)
!     write(6,*)smois(i,1,j),maxsmc(isltyp(i,j))
!     write(6,*)p8w(i,kts+1,j),p8w(i,kts,j),area(i,j),u_phy(i,kts,j)
!     write(6,*)erod(i,j,1),u10(i,j),g,rho_phy(i,kts,j)
!    endif
!
! we donṫ trust the u10,v10 values, is model layers are very thin near surface
!
     if(dz8w(i,kts,j).lt.12.)w10m=sqrt(u_phy(i,kts,j)*u_phy(i,kts,j)+v_phy(i,kts,j)*v_phy(i,kts,j))
     erodin(1,1,1,1)=erod(i,j,1)!/area(i,j)
     erodin(1,1,2,1)=erod(i,j,2)!/area(i,j)
     erodin(1,1,3,1)=erod(i,j,3)!/area(i,j)
!
!  volumetric soil moisture over porosity
!
     if(isltyp(i,j).eq.0)then
      ilwi(1,1)=0
      gwet(1,1)=1.
     else
      gwet(1,1)=smois(i,1,j)/maxsmc(isltyp(i,j))

     endif
!    gwet(1,1)=.1
     airden(1,1)=rho_phy(i,kts,j)
     dxy(1)=area(i,j)
     ipr=0
!    if(j.eq.681)ipr=1
    call source_du( imx,jmx,lmx,nmx, dt, tc, &
                     erodin, ilwi, dxy, w10m, gwet, airden, airmas, &
                     bems,start_month,g,ipr)
!    if(tc(1).gt.tcs)then
!       print *,k,j,tc(1),tc(2)
!       print *,p_dust_1,p_dust_2,tc(3)
!    endif
     chem(i,kts,j,p_dust_1)=max(max_default,(tc(1)+.3125*tc(2))*converi)
     chem(i,kts,j,p_dust_2)=max(max_default,(.67*tc(2)+tc(3))*converi)
     dusthelp(i,j)=max(max_default,tc(2)*converi)
!    if(j.eq.681)then
!     write(6,*)chem(i,kts,j,p_dust_1),chem(i,kts,j,p_dust_2),dusthelp(i,j) 
!     write(6,*)dt,airmas(1,1),dusthelp(i,j)
!     write(6,*)tc(1),tc(2),tc(3)
!    endif
! for output diagnostics
     emis_dust(i,1,j,p_edust1)=bems(1)
     emis_dust(i,1,j,p_edust2)=bems(2)
     emis_dust(i,1,j,p_edust3)=bems(3)
     endif
  enddo
  enddo
  else
! print *,'chem_opt = ',chem_opt,'in gocart_dust2',p_dust_1,p_dust_2
  do j=jts,jte
  do i=its,ite
!
! 
!
!   write(6,*)i,j,xland(i,j)
     if(xland(i,j).lt.1.5 .and. xland(i,j).gt.0.5)then
     ilwi(1,1)=1
     tc(1)=chem(i,kts,j,p_dust_1)*conver
     tc(2)=chem(i,kts,j,p_dust_2)*conver
     tc(3)=chem(i,kts,j,p_dust_3)*conver
     tc(4)=chem(i,kts,j,p_dust_4)*conver
     tc(5)=chem(i,kts,j,p_dust_5)*conver
     w10m(1,1)=sqrt(u10(i,j)*u10(i,j)+v10(i,j)*v10(i,j))
     airmas(1,1)=-(p8w(i,kts+1,j)-p8w(i,kts,j))*area(i,j)/g
!
! we donṫ trust the u10,v10 values, is model layers are very thin near surface
!
     if(dz8w(i,kts,j).lt.12.)w10m=sqrt(u_phy(i,kts,j)*u_phy(i,kts,j)+v_phy(i,kts,j)*v_phy(i,kts,j))
     erodin(1,1,1,1)=erod(i,j,1)!/area(i,j)
     erodin(1,1,2,1)=erod(i,j,2)!/area(i,j)
     erodin(1,1,3,1)=erod(i,j,3)!/area(i,j)
!
!  volumetric soil moisture over porosity
!
     if(isltyp(i,j).eq.0)then
      ilwi(1,1)=0
      gwet(1,1)=1.
     else
      gwet(1,1)=smois(i,1,j)/maxsmc(isltyp(i,j))

     endif
!    gwet(1,1)=.1
     airden(1,1)=rho_phy(i,kts,j)
     dxy(1)=area(i,j)
     ipr=0
!    if(erod(i,j,1).gt.0. .and. gwet(1,1).lt.0.2 .and.j.lt.100)then
!    ipr=1
!    write(6,*)j,smois(i,1,j),maxsmc(isltyp(i,j)),erod(i,j,1),area(i,j)
!    write(6,*)w10m(1,1),airmas(1,1),airden(1,1),gwet(1,1)
!    write(6,*)g,dxy(1)
!    endif
!    dttt=3600.
!    if(erod(i,j,1).gt.0. .and. gwet(1,1).lt.0.2 .and.j.eq.9222 )then
!       ipr=1
!    endif
    call source_du( imx,jmx,lmx,nmx, dt, tc, &
                     erodin, ilwi, dxy, w10m, gwet, airden, airmas, &
                     bems,start_month,g,ipr)
!    write(0,*)tc(1)
!    write(0,*)tc(2)
!    write(0,*)tc(3)
!    write(0,*)tc(4)
!    write(0,*)tc(5)
!    if(erod(i,j,1).gt.0. .and. gwet(1,1).lt.0.2  .and.j.eq.9222)then
!    write(6,*)j,bems(1),bems(2),chem(i,kts,j,p_dust_1),tc(1)*converi
!    endif
     chem(i,kts,j,p_dust_1)=max(max_default,tc(1)*converi)
     chem(i,kts,j,p_dust_2)=max(max_default,tc(2)*converi)
     chem(i,kts,j,p_dust_3)=max(max_default,tc(3)*converi)
     chem(i,kts,j,p_dust_4)=max(max_default,tc(4)*converi)
     chem(i,kts,j,p_dust_5)=max(max_default,tc(5)*converi)
! for output diagnostics
     emis_dust(i,1,j,p_edust1)=bems(1)
     emis_dust(i,1,j,p_edust2)=bems(2)
     emis_dust(i,1,j,p_edust3)=bems(3)
     emis_dust(i,1,j,p_edust4)=bems(4)
     emis_dust(i,1,j,p_edust5)=bems(5)
     endif
  enddo
  enddo
  endif
!

end subroutine gocart_dust_driver

  
  SUBROUTINE source_du( imx,jmx,lmx,nmx, dt1, tc, &
                     erod, ilwi, dxy, w10m, gwet, airden, airmas, &
                     bems,month,g0,ipr)

! ****************************************************************************
! *  Evaluate the source of each dust particles size classes  (kg/m3)        
! *  by soil emission.
! *  Input:
! *         EROD      Fraction of erodible grid cell                (-)
! *                   for 1: Sand, 2: Silt, 3: Clay
! *         DUSTDEN   Dust density                                  (kg/m3)
! *         DXY       Surface of each grid cell                     (m2)
! *         AIRVOL    Volume occupy by each grid boxes              (m3)
! *         NDT1      Time step                                     (s)
! *         W10m      Velocity at the anemometer level (10meters)   (m/s)
! *         u_tresh   Threshold velocity for particule uplifting    (m/s)
! *         CH_dust   Constant to fudge the total emission of dust  (s2/m2)
! *      
! *  Output:
! *         DSRC      Source of each dust type           (kg/timestep/cell) 
! *
! *  Working:
! *         SRC       Potential source                   (kg/m/timestep/cell)
! *
! ****************************************************************************

! USE module_data_gocart
! USE module_data_gocart_dust

  

  INTEGER, INTENT(IN)    :: nmx,imx,jmx,lmx
  REAL*8,    INTENT(IN)    :: erod(imx,jmx,ndcls,ndsrc)
  INTEGER, INTENT(IN)    :: ilwi(imx,jmx),month

  REAL*8,    INTENT(IN)    :: w10m(imx,jmx), gwet(imx,jmx)
  REAL*8,    INTENT(IN)    :: dxy(jmx)
  REAL*8,    INTENT(IN)    :: airden(imx,jmx,lmx), airmas(imx,jmx,lmx)
  REAL*8,    INTENT(INOUT) :: tc(imx,jmx,lmx,nmx)
  REAL*8,    INTENT(OUT)   :: bems(imx,jmx,nmx) 

  REAL*8    :: den(nmx), diam(nmx)
  REAL*8    :: tsrc, u_ts0, cw, u_ts, dsrc, srce
  REAL, intent(in)    :: g0
  REAL    :: rhoa, g,dt1
  INTEGER :: i, j, n, m, k, ipr


  REAL*8                  :: tcmw(nmx), ar(nmx), tcvv(nmx)
  REAL*8                  :: ar_wetdep(nmx), kc(nmx)
  CHARACTER(LEN=20)     :: tcname(nmx), tcunits(nmx)
  LOGICAL               :: aerosol(nmx)


! REAL*8 :: tc1(imx,jmx,lmx,nmx)
! REAL*8, TARGET :: tcms(imx,jmx,lmx,nmx) ! tracer mass (kg; kgS for sulfur case)
! REAL*8, TARGET :: tcgm(imx,jmx,lmx,nmx) ! g/m3

  !-----------------------------------------------------------------------  
  ! sea salt specific
  !-----------------------------------------------------------------------  
! REAL*8, DIMENSION(nmx) :: ssaltden, ssaltreff, ra, rb
! REAL*8 :: ch_ss(nmx,12)

  !-----------------------------------------------------------------------  
  ! emissions (input)
  !-----------------------------------------------------------------------  
! REAL*8 :: e_an(imx,jmx,2,nmx), e_bb(imx,jmx,nmx), &
!         e_ac(imx,jmx,lmx,nmx)

  !-----------------------------------------------------------------------  
  ! diagnostics (budget)
  !-----------------------------------------------------------------------
!  ! tendencies per time step and process
!  REAL, TARGET :: bems(imx,jmx,nmx), bdry(imx,jmx,nmx), bstl(imx,jmx,nmx)
!  REAL, TARGET :: bwet(imx,jmx,nmx), bcnv(imx,jmx,nmx)
!
!  ! integrated tendencies per process
!  REAL, TARGET :: tems(imx,jmx,nmx), tstl(imx,jmx,nmx)
!  REAL, TARGET :: tdry(imx,jmx,nmx), twet(imx,jmx,nmx), tcnv(imx,jmx,nmx)

  ! global mass balance per time step 
  REAL*8 :: tmas0(nmx), tmas1(nmx)
  REAL*8 :: dtems(nmx), dttrp(nmx), dtdif(nmx), dtcnv(nmx)
  REAL*8 :: dtwet(nmx), dtdry(nmx), dtstl(nmx)
  REAL*8 :: dtems2(nmx), dttrp2(nmx), dtdif2(nmx), dtcnv2(nmx)
  REAL*8 :: dtwet2(nmx), dtdry2(nmx), dtstl2(nmx)
  real :: gthresh




  ! executable statemenst
  gthresh=.5

  DO n = 1, nmx
     ! Threshold velocity as a function of the dust density and the diameter
     ! from Bagnold (1941)
     den(n) = den_dust(n)*1.0D-3
     diam(n) = 2.0*reff_dust(n)*1.0D2
     g = g0*1.0E2
     ! Pointer to the 3 classes considered in the source data files
     m = ipoint(n)
     tsrc = 0.0
     DO k = 1, ndsrc
        ! No flux if wet soil 
        DO i = 1,imx
           DO j = 1,jmx
              rhoa = airden(i,j,1)*1.0D-3
              u_ts0 = 0.13*1.0D-2*SQRT(den(n)*g*diam(n)/rhoa)* &
                   SQRT(1.0+0.006/den(n)/g/(diam(n))**2.5)/ &
                   SQRT(1.928*(1331.0*(diam(n))**1.56+0.38)**0.092-1.0) 
!             write(0,*)u_ts0,den(n),diam(n),rhoa,g
              ! Fraction of emerged surfaces (subtract lakes, coastal ocean,..)
!              cw = 1.0 - water(i,j)
              
              ! Case of surface dry enough to erode
              IF (gwet(i,j) < gthresh) THEN
!              IF (gwet(i,j) < 0.5) THEN  !  Pete's modified value
                 u_ts = MAX(0.0D+0,u_ts0*(1.2D+0+2.0D-1*LOG10(MAX(1.0D-3, gwet(i,j)))))
              ELSE
                 ! Case of wet surface, no erosion
                 u_ts = 100.0
              END IF
              srce = frac_s(n)*erod(i,j,m,k)*dxy(j)  ! (m2)
              IF (ilwi(i,j) == 1 ) THEN
                 dsrc = ch_dust(n,month)*srce*w10m(i,j)**2 &
                      * (w10m(i,j) - u_ts)*dt1  ! (kg)
!                 IF (gwet(i,j) < 0.2 .and. ipr.eq.1)write(6,*)n,month,m,ch_dust(n,month),srce,w10m(i,j),u_ts,gwet(i,j)
!                 IF (gwet(i,j) < 0.2 .and. ipr.eq.1)write(6,*)ipoint(m),den_dust(n),erod(i,j,m,k),dxy(j)
!                 IF (gwet(i,j) < 0.2 .and. ipr.eq.1)write(6,*)srce,dsrc,frac_s(n)
!                 IF (gwet(i,j) < 0.2 .and. ipr.eq.1)write(6,*)airmas(i,j,1),dt1
              ELSE 
                 dsrc = 0.0
              END IF
!              dsrc = cw*ch_dust(k)*srce*w10m(i,j)**2 &
!                   * (w10m(i,j) - u_ts)*dt1  ! (kg)
!              dsrc = cw*ch_dust(n,dt(1)%mn)*srce*w10m(i,j)**2 &
!                   * (w10m(i,j) - u_ts)*dt1  ! (kg)
              IF (dsrc < 0.0) dsrc = 0.0
              
              ! Update dust mixing ratio at first model level.
! scale down dust by .6
              tc(i,j,1,n) = tc(i,j,1,n) + .7*dsrc / airmas(i,j,1)
              bems(i,j,n) = .7*dsrc
           END DO
        END DO
     END DO
  END DO
  
END SUBROUTINE source_du


END MODULE GOCART_DUST
