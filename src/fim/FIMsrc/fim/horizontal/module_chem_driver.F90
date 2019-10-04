MODULE MODULE_CHEM_DRIVER
   IMPLICIT NONE
CONTAINS
     subroutine chem_driver(ktau)
!
! FIM version variables
!
  USE module_control, only: dt, yyyymmddhhmm, nvl, nvlp1,ntra,ntrb
  USE module_wrf_control, only: ims,ime,jms,jme,kms,kme,           &
                                ids,ide,jds,jde,kds,kde,           &
                                its,ite,jts,jte,kts,kte,           &
                                nvl_gocart,num_emis_ant,num_moist, &
                                num_ext_coef,num_bscat_coef,num_asym_par,&
                                num_chem,num_soil_layers,numgas,nbands,  &
                                num_emis_vol,CallChemistry,CallBiom
!TODO:  replace use-association of FIM dynamics variables with coupler
  USE module_variables
  USE module_chem_variables
  USE module_wrf_variables,only:exch,pb2d
  USE module_constants
  USE module_chem_constants ,only: p_gocart
  USE module_sfc_variables
! TBH:  Ignore these so PPP doesn't have to translate them
!SMS$IGNORE BEGIN
  USE module_plumerise1, only: plumerise_driver
  USE module_chem_prep_fim,only: chem_prep_fim
  USE module_gocart_seasalt,only: gocart_seasalt_driver
  USE gocart_dust,only: gocart_dust_driver
  USE gocart_dust_afwa,only: gocart_dust_afwa_driver
  USE module_gocart_settling,only: gocart_settling_driver
  USE module_gocart_opt,only: aero_opt
  USE module_vash_settling,only: vash_settling_driver,vashshort_settling_driver
  USE module_wetdep_ls,only: wetdep_ls
  USE module_dms_emis,only: gocart_dmsemis
  USE module_ctrans_grell,only:grelldrvct
  USE module_initial_chem_namelists,only:p_qc,p_dms,p_seas_1,chem_opt,biomass_burn_opt,p_oc1,   &
                                         p_bc1,p_p25,p_p10,p_sulf,p_so2,seas_opt,dust_opt,      &
                                         dmsemis_opt,chem_in_opt,aer_ra_feedback,kemit,p_tr1,p_tr2
  USE module_dry_dep_driver,only:dry_dep_driver
  USE module_gocart_aerosols,only:sum_pm_gocart,gocart_aerosols_driver
  USE module_gocart_chem,only:gocart_chem_driver 
  USE module_data_gocart_chem,only:airmw,mw_so4_aer
  USE module_optical_driver, only: optical_driver
  USE module_aer_opt_out, only: aer_opt_out
  USE module_aer_ra, only: aer_ra
!SMS$IGNORE END
  USE module_chemvars
  USE module_wrfphysvars
  USE module_outtime_chem,only: telapsed=>tchem

  IMPLICIT NONE

  INTEGER, INTENT (IN) :: KTAU

   real :: var_rmv(ims:ime, jms:jme ,num_chem)
   real :: swdown(ims:ime, jms:jme )
   real :: dep_vel_o3(ims:ime, jms:jme )
   real :: e_co(ims:ime, jms:jme )
   real :: raincv_b(ims:ime, jms:jme )
   real :: cu_co_ten(ims:ime, kms:kme,jms:jme )
   REAL :: dusthelp(ims:ime, jms:jme ),seashelp(ims:ime,jms:jme )
   real :: tr_fall(ims:ime, jms:jme,num_chem )

   integer,save :: current_month

   !shc end stuff for MEGAN v2.04

      REAL                     :: dtstep_plume,dtstep, dx, gmt
      REAL                     :: dust_alpha,dust_gamma
      REAL (kind=8) :: curr_secs
!
! Local variables...
!
      INTEGER :: ib,i, j, k, nv,nvv,ksub, dust_emiss_active, seasalt_emiss_active
      INTEGER :: call_plume,call_gocart,julday,current_gmt,curr_mins,current_day,current_year
      REAL :: factor,factor2,conv,mwdry,xlv,maxv,minv
      CHARACTER (LEN=80) :: message,filename 
      CHARACTER(len=9 )                :: jdate 
      real*8 :: t0
! ..
! ..
! .. Intrinsic Functions ..
      INTRINSIC max, min
     call StartTimer(t0)
    dust_alpha=1.0
    dust_gamma=1.6
     call_plume=0
    if(biomass_burn_opt > 0 ) then
       if ((mod(ktau,CallBiom)==0).or.(ktau==1)) then
         call_plume=1
         dtstep_plume=dt*CallBiom
       endif
    endif
    if(chem_opt == 0) then
     write(6,*)'Shouldnt be here with chem_opt = 0 '
    endif
    dusthelp(:,:)=0.
    seashelp(:,:)=0.


!sms$compare_var(st3d   , "begin chem_driver ")
!sms$compare_var(sm3d   , "begin chem_driver ")
!sms$compare_var(ts2d   , "begin chem_driver ")
!sms$compare_var(us2d   , "begin chem_driver ")
!sms$compare_var(sw2d   , "begin chem_driver ")

     call_gocart=0

     if(chem_opt >= 300 .and. chem_opt <  500 ) then
        if ((mod(ktau,CallChemistry)==0).or.(ktau==1)) then
          call_gocart=1
        endif
     endif
     xlv=2.5e6
     dtstep=dt*CallChemistry
     curr_secs=ktau*dt
     curr_mins=curr_secs/60.
     if(ktau.le.1)then
         dtstep=dt
         READ(UNIT=yyyymmddhhmm(5:6), FMT='(I2)') current_month
!sms$ignore begin
         rcav(jts:jte)=rc2d(jts:jte)
         rnav(jts:jte)=rn2d(jts:jte)-rc2d(jts:jte)
!sms$ignore end
     else
!sms$ignore begin
         rcav(jts:jte)=rc2d(jts:jte)-rcav(jts:jte)
         rcav(jts:jte)=max(0.,rcav(jts:jte))
         rnav(jts:jte)=rn2d(jts:jte)-rc2d(jts:jte)-rnav(jts:jte)
         rnav(jts:jte)=max(0.,rnav(jts:jte))
!sms$ignore end
     endif
! make the fim varialbles conform with the chem_driver
!    if(ktau.le.1)then
        READ(UNIT=yyyymmddhhmm(1:4), FMT='(I4)') current_year
        READ(UNIT=yyyymmddhhmm(5:6), FMT='(I2)') current_month
!    endif
     READ(UNIT=yyyymmddhhmm(7:8), FMT='(I2)') current_day
     READ(UNIT=yyyymmddhhmm(9:10), FMT='(I2)') current_gmt
     call GetJdate(yyyymmddhhmm,jdate)                ! Julian date conversion
     READ(UNIT=jdate(3:5), FMT='(I3)')julday
!
! we have to read co2 emissions every few hours
!
     if(chem_opt == 500 ) then
     if(ktau.le.1 .or. mod(curr_mins,180).eq.0)then
        write(filename,'("co2_flux.",I4,I2.2,I2.2,I2.2,".bin")')current_year,current_month,current_day,current_gmt
        write(6,*)filename
!SMS$SERIAL (<emiss_co2,OUT> : default=ignore)  BEGIN
        open(unit=28,file=filename, form="unformatted")
        read(28)emiss_co2
        close(28)
!       maxv=maxval(emiss_co2)
!SMS$SERIAL END
!       write(6,*)'maxv,p_tr1,p_tr2 on input for co2_ant = ',maxv,p_tr1,p_tr2
        call flush(6)
        do j=jts,jte
        emiss_ab(j,p_tr1)=emiss_co2(j)
        emiss_ab(j,p_tr2)=emiss_co2(j)
        enddo
     endif ! every 3 hours
     endif ! chem_opt = 500


     gmt=current_gmt
     call chem_prep_fim(ktau,dt,rh3d,tr3d,tk3d,st3d,sm3d,dp3d,mp3d,  &
             ts2d,us2d,sw2d,pr3d,emiss_ash_mass,emiss_ash_height,    &
             emiss_ash_dt,dm0,emiss_tr_mass,emiss_tr_height,         &
             emiss_tr_dt, VFRAC2d,VTYPE2d,STYPE2d,us3d,vs3d,ws3d,    &
             slmsk2d,zorl2d,exch,pb2d,hf2d,oh_backgd,h2o2_backgd,    &
             no3_backgd,backg_oh,backg_h2o2,backg_no3,p_gocart,      &
             nvl_gocart,ttday,tcosz,gmt,julday,dtstep,ph3d,area,ero1,&
             ero2,ero3,rcav,raincv_b,deg_lat,deg_lon,nvl,nvlp1,ntra, &
             relhum,rri,t_phy,moist,u_phy,v_phy,p_phy,chem,tsk,ntrb, &
             grvity,rd,p1000,cp,erod,emis_ant,emis_vol,e_co,dms_0,   &
             u10,v10,ivgtyp,isltyp,gsw,vegfra,rmol,ust,znt,xland,dxy,&
             t8w,p8w,exch_h,pbl,hfx,xlat,xlong,z_at_w,zmid,dz8w,vvel,&
             rho_phy,smois,num_soil_layers,num_chem,num_moist,       &
             emiss_abu,ebu_in_oc,ebu_in_bc,ebu_in_pm25,ebu_in_pm10,  &
             ebu_in_so2,ebu_in_sulf,emiss_ab,num_emis_ant,           &
             num_emis_vol,kemit,call_gocart,ids,ide, jds,jde, kds,   &
             kde,plumestuff,mean_fct_agtf,mean_fct_agef,             &
             mean_fct_agsv,mean_fct_aggr,firesize_agtf,firesize_agef,&
             firesize_agsv,firesize_aggr,chem_in_opt,                &
             ims,ime, jms,jme, kms,kme,                              &
             its,ite,jts,jte,kts,kte)
  mwdry=28.
!if(mod(ktau,CallChemistry)==0.or.ktau==1) then
! write(6,*)'in chem_driver, now do chemistry = ',ktau,dtstep,CallChemistry,kemit
  if(seas_opt == 1 )then
!    print *,'get seasalt emissions'
     call gocart_seasalt_driver(ktau,dt,rri,t_phy,moist,u_phy,  &
         v_phy,chem,rho_phy,dz8w,u10,v10,p8w,                  &
         xland,xlat,xlong,area,grvity,emis_seas, &
         seashelp,num_emis_seas,num_moist,num_chem,        &
         ids,ide, jds,jde, kds,kde,                                        &
         ims,ime, jms,jme, kms,kme,                                        &
         its,ite, jts,jte, kts,kte)
  endif
  if(dust_opt == 1 )then
     call gocart_dust_driver(ktau,dt,rri,t_phy,moist,u_phy,       &
         v_phy,chem,rho_phy,dz8w,smois,u10,v10,p8w,erod,ivgtyp,isltyp,&
         vegfra,xland,xlat,xlong,gsw,area,grvity,emis_dust,dusthelp,  &
         num_emis_dust,num_moist,num_chem,num_soil_layers,            &
         current_month,                                               &
         ids,ide, jds,jde, kds,kde,                                   &
         ims,ime, jms,jme, kms,kme,                                   &
         its,ite, jts,jte, kts,kte)
  endif
  if(dust_opt == 3 )then
     call gocart_dust_afwa_driver(ktau,dt,rri,t_phy,moist,u_phy,      &
         v_phy,chem,rho_phy,dz8w,smois,u10,v10,p8w,erod,ivgtyp,isltyp,&
         vegfra,xland,xlat,xlong,gsw,area,grvity,emis_dust,dusthelp,  &
         ust,znt,clayfrac,sandfrac,dust_alpha,dust_gamma,             &
         num_emis_dust,num_moist,num_chem,num_soil_layers,            &
         ids,ide, jds,jde, kds,kde,                                   &
         ims,ime, jms,jme, kms,kme,                                   &
         its,ite, jts,jte, kts,kte)
  endif
!endif
if(call_plume == 1 )  then
   call plumerise_driver (ktau,dtstep_plume,num_chem,num_moist,                &
           ebu_no,ebu_co,ebu_co2,ebu_eth,ebu_hc3,ebu_hc5,ebu_hc8,         &
           ebu_ete,ebu_olt,ebu_oli,ebu_pm25,ebu_pm10,ebu_dien,ebu_iso,    &
           ebu_api,ebu_lim,ebu_tol,ebu_xyl,ebu_csl,ebu_hcho,ebu_ald,      &
           ebu_ket,ebu_macr,ebu_ora1,ebu_ora2,ebu_bc,ebu_oc,ebu_so2,      &
           ebu_sulf,                                                      &
           ebu_in_no,ebu_in_co,ebu_in_co2,ebu_in_eth,ebu_in_hc3,ebu_in_hc5,      &
           ebu_in_hc8,ebu_in_ete,ebu_in_olt,ebu_in_oli,ebu_in_pm25,ebu_in_pm10,  &
           ebu_in_dien,ebu_in_iso,ebu_in_api,ebu_in_lim,ebu_in_tol,ebu_in_xyl,   &
           ebu_in_csl,ebu_in_hcho,ebu_in_ald,ebu_in_ket,ebu_in_macr,ebu_in_ora1, &
           ebu_in_ora2,ebu_in_bc,ebu_in_oc,ebu_in_so2,ebu_in_sulf,    &
           mean_fct_agtf,mean_fct_agef,mean_fct_agsv,mean_fct_aggr,              &
           firesize_agtf,firesize_agef,firesize_agsv,firesize_aggr,              &
           t_phy,moist,chem,rho_phy,vvel,u_phy,v_phy,p_phy,z_at_w,zmid,    &
           ids,ide, jds,jde, kds,kde,                                        &
           ims,ime, jms,jme, kms,kme,                                        &
           its,ite, jts,jte, kts,kte)
endif
!if(call_gocart == 1 )  then
if(dmsemis_opt == 1 ) then
     do j=jts,jte
     do i=its,ite
        diaga(5,j)=emis_dust(i,1,j,1)
        diaga(6,j)=chem(i,1,j,p_dms)
     enddo
     enddo
   call gocart_dmsemis(dt,rri,t_phy,u_phy,  &
         v_phy,chem,rho_phy,dz8w,u10,v10,p8w,dms_0,tsk,                  &
         ivgtyp,isltyp,xland,area,grvity,mwdry, &
         num_chem,p_dms,ids,ide, jds,jde, kds,kde,                                        &
         ims,ime, jms,jme, kms,kme,                                        &
         its,ite, jts,jte, kts,kte)
     do j=jts,jte
     do i=its,ite
        diaga(6,j)=-(diaga(6,j)-chem(i,1,j,p_dms))/dt
     enddo
     enddo
endif
if(chem_opt == 300 .and. (dust_opt == 1 .or. seas_opt == 1 .or. dust_opt == 3))then

   call gocart_settling_driver(dt,t_phy,moist,  &
         chem,rho_phy,dz8w,p8w,p_phy,         &
         dusthelp,seashelp,area,grvity,                         &
         num_moist,num_chem,                 &
         ids,ide, jds,jde, kds,kde,                                        &
         ims,ime, jms,jme, kms,kme,                                        &
         its,ite, jts,jte, kts,kte)
endif
!
! 10 volcanic size bins
!
if(CHEM_OPT == 316 )  then
    call vash_settling_driver(dt,t_phy,moist,                              &
         chem,rho_phy,dz8w,p8w,p_phy,area,                                 &
         ash_fall,grvity,num_moist,num_chem    ,                                &
         ids,ide, jds,jde, kds,kde,                                        &
         ims,ime, jms,jme, kms,kme,                                        &
         its,ite, jts,jte, kts,kte                                         )
     do j=jts,jte
     do i=its,ite
        ashfall(j)=ash_fall(i,j)
     enddo
     enddo
endif
!
! 4 volcanic size bins
!
if(CHEM_OPT == 317 .or. CHEM_OPT == 502 .or. CHEM_OPT == 300)  then
     do j=jts,jte
     do i=its,ite
        diaga(2,j)=chem(i,1,j,p_p25)
     enddo
     enddo
    call vashshort_settling_driver(dt,t_phy,moist,                              &
         chem,rho_phy,dz8w,p8w,p_phy,area,                                 &
         ash_fall,grvity,num_moist,num_chem    ,                                &
         ids,ide, jds,jde, kds,kde,                                        &
         ims,ime, jms,jme, kms,kme,                                        &
         its,ite, jts,jte, kts,kte                                         )
     do j=jts,jte
     do i=its,ite
        diaga(2,j)=(chem(i,1,j,p_p25)-diaga(2,j))/dt
     enddo
     enddo
     do j=jts,jte
     do i=its,ite
        ashfall(j)=ash_fall(i,j)
     enddo
     enddo
endif
!
! add biomass burning emissions at every timestep
!
if(BIOMASS_BURN_OPT == 1 )  then
  do i=its,ite
     do k=kts,kte-2
       do j=jts,jte
!factro for pm emissions, factor2 for burn emissions
         factor=dt*rri(i,k,j)/dz8w(i,k,j)
         factor2=4.828e-4*dt*rri(i,k,j)/(60.*dz8w(i,k,j))
         chem(i,k,j,p_oc1)=chem(i,k,j,p_oc1)+(ebu_oc(i,k,j))*factor
         chem(i,k,j,p_bc1)=chem(i,k,j,p_bc1)+(ebu_bc(i,k,j))*factor
         chem(i,k,j,p_p25)=chem(i,k,j,p_p25)+(ebu_pm25(i,k,j))*factor
!        chem(i,k,j,p_p10)=chem(i,k,j,p_p10)+(ebu_pm10(i,k,j))*factor

         chem(i,k,j,p_so2)=chem(i,k,j,p_so2)+ebu_so2(i,k,j)*factor2
       enddo
     enddo
  enddo
endif

!
! subgrid convective transport
!
   call grelldrvct(DT,ktau,                       &
              rho_phy,RAINCV_b,chem,tr_fall,              &
              U_phy,V_phy,t_phy,moist,dz8w,p_phy,                       &
              XLV,CP,grvity,rv,z_at_w,cu_co_ten,                         &
              numgas,chem_opt,                                   &
              num_chem,num_moist,                               &
              ids,ide, jds,jde, kds,kde,                        &
              ims,ime, jms,jme, kms,kme,                        &
              its,ite, jts,jte, kts,kte                         )
     do j=jts,jte
     do i=its,ite
        diaga(3,j)=chem(i,1,j,p_p25)
     enddo
     enddo
    call dry_dep_driver(ktau,dt,                &
               moist,p8w,rri,                                             &
               chem,rho_phy,dz8w,exch_h,hfx,                              &
               ivgtyp,tsk,pbl,ust,znt,zmid,z_at_w,                           &
               xland,dep_vel_o3,grvity,                                        &
               e_co,kemit,numgas,                                         &
               num_chem,num_moist,                                        &
               ids,ide, jds,jde, kds,kde,                                 &
               ims,ime, jms,jme, kms,kme,                                 &
               its,ite, jts,jte, kts,kte                                  )
     do j=jts,jte
     do i=its,ite
        diaga(3,j)=(chem(i,1,j,p_p25)-diaga(3,j))/dt
        diaga(4,j)=dep_vel_o3(i,j)
     enddo
     enddo
!
! ls wet deposition
!
     if(chem_opt .ne. 500) then
     call wetdep_ls(dt,chem,rnav,moist,rho_phy,var_rmv,num_moist, &
         num_chem,numgas,p_qc,dz8w,vvel,chem_opt,                 &
         ids,ide, jds,jde, kds,kde,                                        &
         ims,ime, jms,jme, kms,kme,                                        &
         its,ite, jts,jte, kts,kte                                         )
      endif

!if(mod(ktau,CallChemistry)==0.or.ktau==1) then
if(call_gocart == 1)then
! write(6,*)'in chem_driver, now do gocart chemistry and aod '
       call gocart_chem_driver(ktau,dtstep, gmt,julday,t_phy,moist,                                      &
         chem,rho_phy,dz8w,p8w,backg_oh,backg_h2o2,backg_no3,         &
         area,grvity,xlat,xlong,ttday,tcosz, &
         chem_opt,num_chem,num_moist,                                      &
         ids,ide, jds,jde, kds,kde,                                        &
         ims,ime, jms,jme, kms,kme,                                        &
         its,ite, jts,jte, kts,kte                                         )
      call gocart_aerosols_driver(ktau,dtstep,t_phy,moist,  &
         chem,rho_phy,dz8w,p8w,area,grvity,         &
         chem_opt,num_chem,num_moist,                                      &
         ids,ide, jds,jde, kds,kde,                                        &
         ims,ime, jms,jme, kms,kme,                                        &
         its,ite, jts,jte, kts,kte                                         )
      if(aer_ra_feedback == 2 )then
         call aero_opt('sw',dz8w,chem             &
                   ,rri,relhum,aod                       &
                   ,extt,ssca,asympar,num_chem                &
                   ,ids,ide, jds,jde, kds,kde                     &
                   ,ims,ime, jms,jme, kms,kme                     &
                   ,its,ite, jts,jte, kts,kte )

      endif
      if(aer_ra_feedback == 1 )then
              call optical_driver(curr_secs,dtstep,          &
               chem,dz8w,rri,relhum,                               &
!              h2oai,h2oaj,                                        &
               tauaersw,gaersw,waersw,bscoefsw,tauaerlw,           &
               l2aer,l3aer,l4aer,l5aer,l6aer,l7aer,                &
               num_chem,chem_opt,ids,ide, jds,jde, kds,kde,        &      
               ims,ime, jms,jme, kms,kme,                          &
               its,ite, jts,jte, kts,kte)
              call aer_opt_out(aod,dz8w,                           &
                    ext_coeff,bscat_coeff,asym_par,                &
                    tauaersw,gaersw,waersw,tauaerlw,               &
                    num_ext_coef,num_bscat_coef,num_asym_par,      &
                    ids,ide, jds,jde, kds,kde,                     &
                    ims,ime, jms,jme, kms,kme,                     &
                    its,ite, jts,jte, kts,kte )
              call aer_ra(dz8w                                     &
                   ,extt,ssca,asympar,nbands                       &
                   ,tauaersw,gaersw,waersw,tauaerlw                &
                   ,ids,ide, jds,jde, kds,kde                      &
                   ,ims,ime, jms,jme, kms,kme                      &
                   ,its,ite, jts,jte, kts,kte )


      endif
endif
!      vcsulf_old(its:ite,kts:kte,jts:jte) = &
!           max(chem(its:ite,kts:kte,jts:jte,p_sulf),epsilc)
!      vcso2_old(its:ite,kts:kte,jts:jte) = &
!           max(chem(its:ite,kts:kte,jts:jte,p_so2),epsilc)
!      vch2o2_old(its:ite,kts:kte,jts:jte) = &
!           max(chem(its:ite,kts:kte,jts:jte,p_h2o2),epsilc)

      if(chem_opt < 500 ) then
      call sum_pm_gocart (                                      &
            rri, chem,pm2_5_dry, pm2_5_dry_ec, pm10,                   &
            num_chem,chem_opt,ids,ide, jds,jde, kds,kde,                        &
            ims,ime, jms,jme, kms,kme,                                 &
            its,ite, jts,jte, kts,kte                                  )

       endif
!
!      store aerosol optical variables for feedback in radiation
!
       if(aer_ra_feedback >= 1 )then
          do ib=1,nbands
          do j=jts,jte
          do k=kts,kte
          do i=its,ite
           ext_cof(k,j,ib)=extt(i,k,j,ib)
           sscal(k,j,ib)=ssca(i,k,j,ib)
           asymp(k,j,ib)=asympar(i,k,j,ib)
          enddo
          enddo
          enddo
          enddo
!         print *,'in chem_driver ',nbands,maxval(extt),maxval(ext_cof)
! aod only output
          do j=jts,jte
          do i=its,ite
           aod2d(j)=aod(i,j)
!          aod2d(j)=dep_vel_o3(i,j)
          enddo
          enddo
!         print *,maxval(aod2d)
       endif   ! feedback to radiation
!
! pm25 and pm10 for output , not for tracer options
!
       if(chem_opt < 500) then
       do j=jts,jte
       do k=kts,kte
       do i=its,ite
        pm25(k,j)=pm2_5_dry(i,k,j)
        p10(k,j)=pm10(i,k,j)
       enddo
       enddo
       enddo
       endif
!
! put chem stuff back into tracer array
!
      do nv=1,num_chem
         nvv=ntra+nv
         do j=jts,jte
            do k=kts,kte
               do i=its,ite
                tr3d(k,j,nvv)=max(epsilc,chem(i,k,j,nv))
                trdp(k,j,nvv)=tr3d(k,j,nvv)*dp3d(k,j)
               enddo
            enddo
         enddo
         if(chem_opt == 501 ) then
            do j=jts,jte
            do i=its,ite
              trfall(j,nv)=trfall(j,nv)+tr_fall(i,j,nv)+var_rmv(i,j,nv)
            enddo
            enddo
         endif
      enddo
         if(chem_opt == 500 ) then
         do j=jts,jte
            do k=kts,kte
               do i=its,ite
                 tr1_tavg(k,j)=tr1_tavg(k,j)+chem(i,k,j,p_tr1)
               enddo
            enddo
         enddo
         endif

!sms$compare_var(st3d   , "end chem_driver ")
!sms$compare_var(sm3d   , "end chem_driver ")
!sms$compare_var(ts2d   , "end chem_driver ")
!sms$compare_var(us2d   , "end chem_driver ")
!sms$compare_var(sw2d   , "end chem_driver ")

      call IncrementTimer(t0,telapsed)

!endif  ! if(chem_opt >= 0) then

     END subroutine chem_driver

END MODULE MODULE_CHEM_DRIVER

