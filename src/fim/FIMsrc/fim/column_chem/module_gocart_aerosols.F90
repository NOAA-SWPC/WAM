MODULE module_gocart_aerosols

    USE module_initial_chem_namelists,only: p_bc1,p_bc2,p_oc1,p_oc2,      &
              p_dust_1,p_dust_2,p_dust_3,p_dust_4,p_dust_5,p_vash_1,     &
              p_seas_1,p_seas_2,p_seas_3,p_seas_4,p_sulf,p_p25,p_so2
    USE module_data_gocart_chem,only:airmw
!odule_species_decs
! USE module_model_constants, only: mwdry
    INTEGER, PARAMETER ::NBC1=1, NOC1=2, NBC2=3, NOC2=4

CONTAINS
  subroutine gocart_aerosols_driver(ktau,dt,t_phy,moist,                   &
         chem,rho_phy,dz8w,p8w,area,g,         &
         chem_opt,num_chem,num_moist,                                      &
         ids,ide, jds,jde, kds,kde,                                        &
         ims,ime, jms,jme, kms,kme,                                        &
         its,ite, jts,jte, kts,kte                                         )
! USE module_configure
! USE module_state_description
  IMPLICIT NONE
!  TYPE(grid_config_rec_type),  INTENT(IN   )    :: config_flags

   INTEGER,      INTENT(IN   ) :: ktau,                     &
                                  chem_opt,num_chem,num_moist,             &
                                  ids,ide, jds,jde, kds,kde,               &
                                  ims,ime, jms,jme, kms,kme,               &
                                  its,ite, jts,jte, kts,kte
    REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_moist ),                &
         INTENT(IN ) ::                                   moist
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),                 &
         INTENT(INOUT ) ::                                   chem
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme ),                        &
          INTENT(IN   ) ::                     t_phy,                      &
                                            dz8w,p8w,                      &
                                              rho_phy
   REAL,  DIMENSION( ims:ime ,  jms:jme ),                                 &
          INTENT(IN   ) ::     area

  REAL, INTENT(IN   ) :: dt,g
  integer :: ndt1,nmx,i,j,k,imx,jmx,lmx
  real*8, DIMENSION (1,1,1) :: tmp,airden,airmas
  REAL*8    :: chmlos(1,1,1,4)
  REAL*8    :: bchmlos(1,1,4)
  REAL*8 :: pc2(1,1,1,2)
  REAL*8 :: tc(4)
  real, parameter :: mw_c = 12.
  real mwdry,tt1,tt2
       mwdry=airmw
       imx=1
       jmx=1
       lmx=1
       nmx=4
       ndt1=ifix(dt)
!
!
       chmlos = 0.
       bchmlos = 0.
       do j=jts,jte
       do k=kts,kte-1
       do i=its,ite
          airmas(1,1,1)=-(p8w(i,k+1,j)-p8w(i,k,j))*area(i,j)/g
          pc2(1,1,1,1)=0.
          pc2(1,1,1,2)=0.
          tc(1)=chem(i,k,j,p_bc1)/mw_c*mwdry*1.d-9
          tc(2)=chem(i,k,j,p_oc1)/mw_c*mwdry*1.d-9
          tc(3)=chem(i,k,j,p_bc2)/mw_c*mwdry*1.d-9
          tc(4)=chem(i,k,j,p_oc2)/mw_c*mwdry*1.d-9
          tt1=tc(3)

         CALL chem_1(imx,jmx,lmx, nmx, ndt1, airmas, tc, &
              chmlos, bchmlos, pc2)
         CALL chem_2(imx,jmx,lmx, nmx, ndt1, airmas, tc, pc2)
         tt2 = tc(3) -tt1
          chem(i,k,j,p_bc1)=tc(1)/mwdry*mw_c*1.e9
          chem(i,k,j,p_oc1)=tc(2)/mwdry*mw_c*1.e9
          chem(i,k,j,p_bc2)=tc(3)/mwdry*mw_c*1.e9
          chem(i,k,j,p_oc2)=(tc(4)+8.*tt2)/mwdry*mw_c*1.e9

       enddo
       enddo
       enddo
end subroutine gocart_aerosols_driver
       subroutine sum_pm_gocart (                                      &
            alt, chem,pm2_5_dry, pm2_5_dry_ec, pm10,                   &
            num_chem,chem_opt,ids,ide, jds,jde, kds,kde,                        &
            ims,ime, jms,jme, kms,kme,                                 &
            its,ite, jts,jte, kts,kte                                  )
!  USE module_configure
!  USE module_state_description
   IMPLICIT NONE
   REAL, PARAMETER :: mwso4 = 96.066
   INTEGER,      INTENT(IN   ) :: chem_opt,ids,ide, jds,jde, kds,kde,               &
                                  ims,ime, jms,jme, kms,kme,               &
                                  its,ite, jts,jte, kts,kte,num_chem
    REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                          &
         INTENT(INOUT ) :: pm2_5_dry, pm2_5_dry_ec, pm10     
    REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                          &
         INTENT(IN ) :: alt
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),                 &
         INTENT(IN ) ::                                   chem
   real minv,maxv,d_2_5,s_2_5,d_10,sulfate,mwdry
   integer i,j,k,ii,jj,n,maxp,maxs,maxd
    mwdry=airmw
    d_2_5=0.3125
    s_2_5=0.75
    d_10=0.67

!
! sum up pm2_5 and pm10 output
!
!     maxv=maxval(chem(its:ite, kts:kte, jts:jte, :))
!     minv=minval(chem(its:ite, kts:kte, jts:jte, :))
!TBH commented out for now -- will not work in parallel without reduce
!TBH      write(6,*)'in sumpm ',minv,maxv
      pm2_5_dry(its:ite, kts:kte, jts:jte)    = 0.
      pm10(its:ite, kts:kte, jts:jte)    = 0.
      pm2_5_dry_ec(its:ite, kts:kte, jts:jte) = 0.
!TBH commented out for now while debugging
!TBH      write(6,*)'p_p25 = ',p_p25,p_dust_1
!TBH      write(6,*)kts,kte,kms,kme,jts,jte,jms,jme
!     do n=p_so2,p_dust_1
!     do j=jts,jte                    
!     do k=kts,kte
!     do i=its,ite
!       pm10(i,k,j)=chem(i,k,j,n)
!       if(n.eq.p_dust_1 .and. chem(i,k,j,n).lt.0.)write(6,*)'ERROR:',i,j,k,chem(i,k,j,n)
!     enddo
!     enddo
!     enddo
!     maxv=maxval(pm10(its:ite, kts:kte, jts:jte))
!     minv=minval(pm10(its:ite, kts:kte, jts:jte))
!     write(6,*)'in sumpm0 ',n,minv,maxv
!     enddo
!     minv=minval(alt(its:ite, kts:kte, jts:jte))
!     write(6,*)'alt',minv
      do j=jts,jte                    
         do k=kts,kte
         do i=its,ite
            sulfate=chem(i,k,j,p_sulf)*mwso4/mwdry*1.e3
            do n=p_p25,p_dust_1
               pm2_5_dry(i,k,j) = pm2_5_dry(i,k,j)+chem(i,k,j,n)        
            enddo
            if(chem_opt.eq.300)then
               pm2_5_dry(i,k,j) = pm2_5_dry(i,k,j)+chem(i,k,j,p_dust_2)*d_2_5     &
                                            +chem(i,k,j,p_seas_1)           &        
                                            +chem(i,k,j,p_seas_2)*s_2_5     &
                                            +sulfate
            else if(chem_opt .eq. 317) then
               pm2_5_dry(i,k,j) = pm2_5_dry(i,k,j)+chem(i,k,j,p_dust_2)*d_2_5     &
                                            +chem(i,k,j,p_seas_1)           &        
                                            +chem(i,k,j,p_seas_2)*s_2_5     &
                                            +chem(i,k,j,p_vash_1)           &
                                            +sulfate
            else
               pm2_5_dry(i,k,j) = pm2_5_dry(i,k,j)+chem(i,k,j,p_seas_1)           &        
                                            +sulfate
            endif
   
            !Convert the units from mixing ratio to concentration (ug m^-3)
            pm2_5_dry(i,k,j)    = pm2_5_dry(i,k,j) / alt(i,k,j)
      enddo
      enddo
      enddo
      maxd=max(p_dust_2,p_dust_3)
      maxp=max(p_dust_2,p_dust_4)
      maxs=max(p_seas_2,p_seas_3)
!     print *,'maxd,maxp,maxs = ',maxd,maxp,maxs
      do j=jts,jte
         do k=kts,kte
            do i=its,ite
               sulfate=chem(i,k,j,p_sulf)*mwso4/mwdry*1.e3
               do n=p_p25,maxd
                  pm10(i,k,j) = pm10(i,k,j)+chem(i,k,j,n)        
               enddo
               do n=p_seas_1,maxs
                  pm10(i,k,j) = pm10(i,k,j)+chem(i,k,j,n)        
               enddo
               pm10(i,k,j) = pm10(i,k,j) + sulfate                 &
                            +chem(i,k,j,maxp)*d_10
               pm10(i,k,j) = pm10(i,k,j)/ alt(i,k,j)
!              if(pm10(i,k,j).gt.10000)write(6,*)j,k,sulfate,pm10(i,k,j)
            enddo
         enddo
      enddo
!     maxv=maxval(pm10(its:ite, kts:kte, jts:jte))
!     minv=minval(pm10(its:ite, kts:kte, jts:jte))
!     write(6,*)'in sumpm2 ',maxv,minv

end subroutine sum_pm_gocart

SUBROUTINE chem_1(imx,jmx,lmx, nmx, &
     ndt1, airm, tc, chmlos, bchmlos, pc2)
! ****************************************************************************
! **                                                                        **
! **  For tracers with dry deposition, the loss rate of dry dep is combined **
! **  in chem loss term. Assuming a conversion time of 2.5 days (1-4 days,  **
! **  Lynn Russell) of conversion time from hydrophobic to hydrophilic.     **
! **     BC1 --> BC2         k = 4.63e-6                                    **
! **     OC1 --> OC2         k = 4.63e-6                                    **
! **     BC1 --> drydep      kd = DRYDf (sec-1)                             **
! **     OC1 --> drydep      kd = DRYDf (sec-1)                             **
! **                                                                        **
! ****************************************************************************


  IMPLICIT NONE

  INTEGER, INTENT(IN)    :: lmx, nmx,imx,jmx, ndt1
  REAL*8,    INTENT(IN)    :: airm(imx,jmx,lmx)
  REAL*8,    INTENT(INOUT) :: tc(imx,jmx,lmx,nmx)
  REAL*8,    INTENT(INOUT) :: chmlos(imx,jmx,lmx,nmx)
  REAL*8,    INTENT(INOUT) :: bchmlos(imx,jmx,nmx)
  REAL*8,    INTENT(OUT)   :: pc2(imx,jmx,lmx,2)

  REAL*8 :: r1, c0, r2, rkt, c1
  INTEGER :: np, n, i, j, l

  ! executable statements

  r1 = 4.63E-6
  
  DO n = 1,nmx
     IF (n == NBC1 .OR. n == NOC1) THEN
        IF (n == NBC1) np = 1
        IF (n == NOC1) np = 2
        DO l = 1,lmx
           DO j = 1,jmx
              DO i = 1,imx
                 
                 c0 = tc(i,j,l,n)
                 r2 = 0.0 ! used to be loss due to dry dep
                 rkt = (r1 + r2) * REAL(ndt1)
                 
                 c1 = c0 *  EXP(-rkt)
                 c1 = MAX(c1, 1.0D-32)
                 tc(i,j,l,n) = c1
                 
                 pc2(i,j,l,np) = (c0 - c1) * r1/(r1 + r2)
                 
                 !   Diagnostics:  
!                chmlos(i,j,l,np) = chmlos(i,j,l,np) + pc2(i,j,l,n)*airm(i,j,l)
                 chmlos(i,j,l,n) = chmlos(i,j,l,n) + pc2(i,j,l,np)*airm(i,j,l)
                 
              END DO
           END DO
        END DO

        DO j = 1,jmx
           DO i = 1,imx
              bchmlos(i,j,n) = bchmlos(i,j,n) + SUM(chmlos(i,j,:,n))
           END DO
        END DO

     END IF
  END DO
  
END SUBROUTINE chem_1

SUBROUTINE chem_2(imx,jmx,lmx, nmx, &
                  ndt1, airm, tc,  pc2)
! ****************************************************************************
! *                                                                          *
! *  C2 = C2_0 * exp(-kt) + PC2/kt * (1.-exp(-kt))                           *
! *    where k = dry deposition, C2 = BC2 or OC2.                            *
! *                                                                          *
! ****************************************************************************


  IMPLICIT NONE

  INTEGER, INTENT(IN)    :: lmx,imx,jmx, nmx, ndt1
  REAL*8,    INTENT(IN)    :: airm(imx,jmx,lmx)
  REAL*8,    INTENT(INOUT) :: tc(imx,jmx,lmx,nmx)
  REAL*8,    INTENT(IN)    :: pc2(imx,jmx,lmx,2)

  INTEGER :: np, n, i, j, l
  REAL*8  :: c0, pp, rkt, c1

  ! executable statements
  
  DO n = 1,nmx
     IF (n == NBC2 .OR. n == NOC2) THEN
        IF (n == NBC2) np = 1
        IF (n == NOC2) np = 2
!CMIC$ doall autoscope
        DO l = 1,lmx
           DO j = 1,jmx
              DO i = 1,imx
                 
                 c0 = tc(i,j,l,n)
!                 pp = pc2(i,j,l,n)
                 pp = pc2(i,j,l,np)
                 c1 = c0 + pp
                 
                 c1 = MAX(c1, 1.0D-32)
                 tc(i,j,l,n) = c1
                 
                 
              END DO
           END DO
        END DO
     END IF
  END DO
  
END SUBROUTINE chem_2

END MODULE module_gocart_aerosols
