MODULE MODULE_CHEM_PREP_FIM
USE module_data_gocart_chem,only:airmw,mw_so4_aer
USE module_gocart_chem,only:szangle
!USE module_chem_namelist_defaults
  USE module_initial_chem_namelists ! ,only:p_bc1,p_oc1,p_sulf,p_e_bc,p_e_oc,p_e_sulf
CONTAINS
     subroutine chem_prep_fim(ktau,dtstep,rh3d,tr3d,tk3d,st3d,sm3d,dp3d,mp3d,ts2d,us2d,   &
             sw2d,pr3d,emiss_ash_mass,emiss_ash_height,emiss_ash_dt,dm0,  &
             emiss_tr_mass,emiss_tr_height,emiss_tr_dt, &
             VFRAC2d,VTYPE2d,STYPE2d,us3d,vs3d,ws3d,slmsk2d,zorl2d,exch,pb2d,hf2d,&
             oh_backgd,h2o2_backgd,no3_backgd,backg_oh,backg_h2o2,backg_no3,p_gocart,nvl_gocart, &
             ttday,tcosz,gmt,julday,dtstepc,                                                     &
             ph3d,area,ero1,ero2,ero3,rcav,raincv_b,deg_lat,deg_lon,nvl,nvlp1,ntra,      &
             relhum,rri,t_phy,moist,u_phy,v_phy,p_phy,chem,tsk,ntrb,g,rd,p1000,cp,           &
             erod,emis_ant,emis_vol,e_co,dms_0,&
             u10,v10,ivgtyp,isltyp,gsw,vegfra,rmol,ust,znt,xland,dxy,t8w,p8w,exch_h,pbl,hfx,   &
             xlat,xlong,z_at_w,zmid,dz8w,vvel,rho_phy,smois,num_soil_layers,num_chem,num_moist,&
             emiss_abu,ebu_in_oc,ebu_in_bc,ebu_in_pm25,ebu_in_pm10,ebu_in_so2,ebu_in_sulf,     &
             emiss_ab,num_emis_ant,num_emis_vol,kemit,call_gocart,ids,ide, jds,jde, kds,kde,                &
             plumestuff,mean_fct_agtf,mean_fct_agef,mean_fct_agsv,mean_fct_aggr,              &
             firesize_agtf,firesize_agef,firesize_agsv,firesize_aggr,              &
             chem_in_opt,ims,ime, jms,jme, kms,kme,                                      &
             its,ite,jts,jte,kts,kte)
!
! input fim variables
!
IMPLICIT NONE
INTEGER,      INTENT(IN   ) :: chem_in_opt,ktau,nvl,nvlp1,ntra,ntrb,nvl_gocart,call_gocart
REAL,         INTENT(IN   ) :: g,rd,p1000,cp,dtstep,dtstepc,gmt
real, intent(in) :: dp3d(nvl,jms:jme)      ! del p between coord levels (pascals)
real, intent(in) :: mp3d(nvl,jms:jme)      ! Montgomery Potential (m^2/s^2)
real, intent(inout) :: tk3d(nvl,jms:jme)      ! temperature, kelvin
real, intent(in) :: rh3d(nvl,jms:jme)      ! temperature, kelvin
real, intent(in) :: exch(nvl,jms:jme)      !
real, intent(in) :: oh_backgd(nvl_gocart,jms:jme)      !
real, intent(in) :: h2o2_backgd(nvl_gocart,jms:jme)      !
real, intent(in) :: no3_backgd(nvl_gocart,jms:jme)      !
real, intent(in) :: tr3d(nvl,jms:jme,ntra+ntrb)  ! 1=pot.temp, 2=water vapor, 3=cloud water, 4=ozone
real, intent(in) :: st3d(4,jms:jme)        ! soil temperature
real, intent(in) :: sm3d(4,jms:jme)        ! soil moisture
real, intent(in) :: ts2d(jms:jme)          ! skin temperature
real, intent(in) :: us2d(jms:jme)          ! friction velocity/equivalent momentum flux
real, intent(in) :: pb2d(jms:jme)          ! 
real, intent(in) :: rcav(jms:jme)          ! 
real, intent(in) :: hf2d(jms:jme)          ! 
real, intent(in) :: sw2d(jms:jme)          ! downward short-wave radiation flux
real, intent(in) :: pr3d(nvlp1,jms:jme)    ! pressure (pascal)
!real, intent(in) :: ex3d(nvlp1,jms:jme)    ! exner function
real, intent(in) :: ph3d(nvlp1,jms:jme)    ! geopotential (=gz), m^2/s^2
real, intent(in) :: emiss_ab(jms:jme,num_emis_ant)           ! 
real, intent(in) :: emiss_abu(jms:jme,num_emis_ant)           ! 
real, intent(in) :: plumestuff(jms:jme,8)           ! 
real, intent(in) :: ero1(jms:jme)           ! 
real, intent(in) :: ero2(jms:jme)           ! 
real, intent(in) :: ero3(jms:jme)           ! 
real, intent(inout) :: emiss_ash_mass(jms:jme)           ! 
real, intent(inout) :: emiss_ash_height(jms:jme)           ! 
real, intent(in) :: emiss_ash_dt(jms:jme)           ! 
real, intent(in) :: emiss_tr_mass(jms:jme)           ! 
real, intent(in) :: emiss_tr_height(jms:jme)           ! 
real, intent(in) :: emiss_tr_dt(jms:jme)           ! 
real, intent(in) :: dm0(jms:jme)           ! 
real, intent(in) :: p_gocart(56)           ! 
real, intent(in) :: area(jms:jme)             ! the area of cell polygon (m**2)
real, dimension (jms:jme), intent(in) :: vfrac2d,VTYPE2d,STYPE2d,zorl2d,slmsk2d
real, dimension (nvl,jms:jme), intent(in) :: us3d,vs3d,ws3d
real, intent(in) :: deg_lat(jms:jme),deg_lon(jms:jme)  ! lat and lon in degrees

   INTEGER,      INTENT(IN   ) :: num_soil_layers,num_chem,num_moist,julday,     &
                                  num_emis_vol,num_emis_ant,ids,ide, jds,jde, kds,kde,        &
                                  kemit,ims,ime, jms,jme, kms,kme,               &
                                  its,ite, jts,jte, kts,kte
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_moist ),                &
         INTENT(OUT ) ::                                   moist
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),                 &
         INTENT(OUT ) ::                                   chem
   REAL, DIMENSION( ims:ime, kms:kemit, jms:jme, num_emis_ant ),                 &
         INTENT(inout ) ::                                   emis_ant
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_emis_vol ),                 &
         INTENT(inout ) ::                                   emis_vol
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         ,               &
          INTENT(OUT   ) ::                                                 &
                                                        rri,               &
                                                      t_phy,               &
                                                      p_phy,               &
                                              relhum, dz8w,p8w,t8w,        &
                                              z_at_w , zmid ,exch_h,       &
                                              u_phy,v_phy,vvel,rho_phy
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         ,               &
          INTENT(OUT   ) ::                                                 &
                                   backg_oh,backg_h2o2,backg_no3
   REAL,  DIMENSION( ims:ime , jms:jme )         ,               &
          INTENT(OUT   ) ::                                                 &
                              ttday,tcosz
   REAL,DIMENSION( ims:ime , jms:jme )                  ,               &
          INTENT(OUT   ) ::                                                 &
         ebu_in_oc,ebu_in_bc,ebu_in_pm25,ebu_in_pm10,ebu_in_so2,ebu_in_sulf
   REAL,DIMENSION( ims:ime , jms:jme )                  ,               &
          INTENT(OUT   ) ::                                                 &
        mean_fct_agtf,mean_fct_agef,mean_fct_agsv,mean_fct_aggr,            &
        firesize_agtf,firesize_agef,firesize_agsv,firesize_aggr
   INTEGER,DIMENSION( ims:ime , jms:jme )                  ,               &
          INTENT(OUT   ) ::                                                 &
                                                     ivgtyp,               &
                                                     isltyp
   REAL, DIMENSION( ims:ime, jms:jme,3)::&
         erod

   REAL,  DIMENSION( ims:ime , jms:jme )                   ,               &
          INTENT(OUT   ) ::                                                 &
                                                     u10,                  &
                                                     v10,                  &
                                                     gsw,                  &
                                                  vegfra,                  &
                                                     rmol,                 &
                                                     ust,                  &
                                                     xland,                &
                                                     xlat,e_co,dms_0,      &
                                                     xlong,tsk,raincv_b,   &
                                                     dxy,znt,pbl,hfx
  REAL, DIMENSION( ims:ime, num_soil_layers, jms:jme ) ,      &
      INTENT(OUT) ::                               smois
   integer i,j,k,kk,nv,jmax,jmaxi,l,ll,n,ndystep,ixhour
   real maxv,factor,factor2,pu,pl,aln,pwant,rlat
   real thv,xhour,xmin,gmtp,xlonn,xtime,real_time
   real, DIMENSION (1,1) :: sza,cosszax
   real, DIMENSION (jms:jme) :: so2_mass
! volcanic stuff
   integer :: ko,k_final,k_initial,kl,kk4,curr_hours,curr_secs
   real :: percen_mass_umbrel,x1,base_umbrel,ashz_above_vent,base_umbrel2
   real, DIMENSION (kms:kme) :: vert_mass_dist
   real :: eh,h1,h2,h3,h4,h5,h6
! .. Intrinsic Functions ..
      INTRINSIC max, min, float
!
     so2_mass(:)=0.
!    h1=1-0
!    h2=7-0
!    h3=19-0
!    h4=69-0
!    h5=84-0
     h1=9-0
     h2=16-0
     h3=58-0
     h4=79-0
     h5=109-0
     h6=129-0
!
! use these values for real-time default (if volcano starts at h1 =0 )
! if h1 ne.0, then care has to be taken that ash_height is correct (below is hardwird for special case)
! ash_height can come in through FIMnamelist (read in in chem_init.F90)
!
     h1=240
     h2=240
     h3=240
     h4=240
     h5=240
     h6=240
       percen_mass_umbrel=.75 
       base_umbrel=.25    ! fraction
       base_umbrel2=1.    ! evenly distribution
       real_time=float(ktau)*dtstep/60.


       if(ktau.le.1)then
         emis_ant(:,:,:,:)=0.
         emis_vol(:,:,:,:)=0.
       endif
       e_co(:,:)=0.
       do i=its,ite
       do j=jts,jte
        z_at_w(i,kts,j)=max(0.,ph3d(kts,j)/g)
       enddo
       enddo
       do i=its,ite
       do k=kts,kte
       do j=jts,jte
        dz8w(i,k,j)=(ph3d(k+1,j)-ph3d(k,j))/g
        if(dz8w(i,k,j).lt.0.)dz8w(i,k,j)=-dz8w(i,k,j)
        z_at_w(i,k+1,j)=z_at_w(i,k,j)+dz8w(i,k,j)
       enddo
       enddo
       enddo
       do i=its,ite
       do k=kts,kte+1
       do j=jts,jte
!       z_at_w(i,k,j)=ph3d(k,j)/g
        p8w(i,k,j)=pr3d(k,j)
       enddo
       enddo
       enddo
       do i=its,ite
       do j=jts,jte
         raincv_b(i,j)=rcav(j)
         pbl(i,j)=pb2d(j)
         dms_0(i,j)=dm0(j)
         hfx(i,j)=hf2d(j)
         erod(i,j,1)=ero1(j)
         erod(i,j,2)=ero2(j)
         erod(i,j,3)=ero3(j)
         xlat(i,j)=deg_lat(j)
         xlong(i,j)=deg_lon(j)
         ust(i,j)=us2d(j)
         tsk(i,j)=ts2d(j)
         ivgtyp(i,j)=VTYPE2d(j)
         isltyp(i,j)=STYPE2d(j)
         gsw(i,j)=sw2d(j)
         vegfra(i,j)=VFRAC2d(j)
!         if(j.eq.681)write(6,*)ivgtyp(i,j),isltyp(i,j),vegfra(i,j),slmsk2d(j)
!        if(ivgtyp(i,j).ne.0)write(6,*)i,j,ivgtyp(i,j),isltyp(i,j),vegfra(i,j),pb2d(j),VTYPE2d(j)
         rmol(i,j)=0.
         znt(i,j)=zorl2d(j)*.01
!SLMSK   - SEA(0),LAND(1),ICE(2) MASK
         xland(i,j)=1.
         if(slmsk2d(j).eq.0)xland(i,j)=0.
         if(slmsk2d(j).eq.1)xland(i,j)=1.
         if(slmsk2d(j).eq.2)xland(i,j)=2.
!        if (slmsk2d(j).gt.0.)write(6,*)j,slmsk2d(j)
         dxy(i,j)=area(j)
         u10(i,j)=us3d(1,j)
         v10(i,j)=vs3d(1,j)
       enddo
       enddo
       factor=0.
       jmax=0
       jmaxi=0
       k=1
       if(p_bc2 .gt. 1)then  ! "regular" chem options
       do i=its,ite
       do j=jts,jte
         k=1
         emis_ant(i,k,j,p_e_bc)=emiss_ab(j,p_e_bc)
         emis_ant(i,k,j,p_e_oc)=emiss_ab(j,p_e_oc)
         emis_ant(i,k,j,p_e_sulf)=emiss_ab(j,p_e_sulf)
         emis_ant(i,k,j,p_e_so2)=emiss_ab(j,p_e_so2)
!
         ebu_in_oc(i,j)=emiss_abu(j,p_e_oc)
         ebu_in_bc(i,j)=emiss_abu(j,p_e_bc)
         ebu_in_pm25(i,j)=emiss_abu(j,p_e_pm_25)
         ebu_in_pm10(i,j)=emiss_abu(j,p_e_pm_10)
         ebu_in_so2(i,j)=emiss_abu(j,p_e_so2)
         ebu_in_sulf(i,j)=0. ! for now
         mean_fct_agtf(i,j)=plumestuff(j,1)
         mean_fct_agef(i,j)=plumestuff(j,2)
         mean_fct_agsv(i,j)=plumestuff(j,3)
         mean_fct_aggr(i,j)=plumestuff(j,4)
         firesize_agtf(i,j)=plumestuff(j,5)
         firesize_agef(i,j)=plumestuff(j,6)
         firesize_agsv(i,j)=plumestuff(j,7)
         firesize_aggr(i,j)=plumestuff(j,8)
       enddo
       enddo
       else if (p_tr2 .gt. 1) then  ! tracer options
! tracer run
          do i=its,ite
          do j=jts,jte
          k=kts
            emis_ant(i,k,j,p_e_tr1)=emiss_ab(j,p_e_tr1)
            emis_ant(i,k,j,p_e_tr2)=emiss_ab(j,p_e_tr2)
          enddo
          enddo
        else if ((p_tr2 .gt. 1) .and. (p_bc2 .gt. 1))then
          stop 'in chem_prep_fim, 111'
       endif
       do i=its,ite
       do k=kts,kte
       do j=jts,jte
         thv=tr3d(k,j,1)/(1.+0.6078*tr3d(k,j,2))
         tk3d(k,j)=thv*(.5*(p8w(i,k,j)+p8w(i,k+1,j))/p1000)**(rd/cp)
       enddo
       enddo
       enddo
       do i=its,ite
       do k=kts,kte+1
       kk=min(k,kte)
       do j=jts,jte
        zmid(i,k,j)=.5*(ph3d(kk+1,j)+ph3d(kk,j))/g
        dz8w(i,k,j)=z_at_w(i,kk+1,j)-z_at_w(i,kk,j)
        t_phy(i,k,j)=tk3d(kk,j)
!       relhum(i,k,j)=rh3d(kk,j)
        p_phy(i,k,j)=.5*(p8w(i,kk,j)+p8w(i,kk+1,j))
        u_phy(i,k,j)=us3d(kk,j)
        exch_h(i,k,j)=exch(kk,j)
        v_phy(i,k,j)=vs3d(kk,j)
        rho_phy(i,k,j)= p_phy(i,k,j)/(RD*T_phy(i,k,j)*(1.+.608*tr3d(kk,j,2)))
        rri(i,k,j)=1./rho_phy(i,k,j)
        vvel(i,k,j)=-ws3d(kk,j)*rri(i,k,j)/g
        moist(i,k,j,:)=0.
        moist(i,k,j,1)=tr3d(kk,j,2)
        if(t_phy(i,k,j).gt.265.)then
           moist(i,k,j,2)=tr3d(kk,j,3)
           moist(i,k,j,3)=0.
           if(moist(i,k,j,2).lt.1.e-8)moist(i,k,j,2)=0.
        else
           moist(i,k,j,2)=0.
           moist(i,k,j,3)=tr3d(kk,j,3)
           if(moist(i,k,j,3).lt.1.e-8)moist(i,k,j,3)=0.
        endif
        relhum(i,k,j) = .95
        relhum(i,k,j) = MIN( .95, moist(i,k,j,1) / &
             (3.80*exp(17.27*(t_phy(i,k,j)-273.)/ &
             (t_phy(i,k,j)-36.))/(.01*p_phy(i,k,j))))
        relhum(i,k,j)=max(0.1,relhum(i,k,j))

       enddo
       enddo
       enddo
       do nv=1,num_chem
       do i=its,ite
       do k=kts,kte+1
         kk=min(k,kte)
       do j=jts,jte
         chem(i,k,j,nv)=tr3d(kk,j,ntra+nv)
       enddo
       enddo
       enddo
       enddo
!
! gocart background fields only if gocart is called
!
       if(call_gocart.eq.1)then
       do i=its,ite
       do j=jts,jte
          do k=kts,kte
            do ll=2,nvl_gocart
              l=ll
              if(p_gocart(l).lt..01*p_phy(i,k,j))exit
            enddo
            pu=alog(p_gocart(l))
            pl=alog(p_gocart(l-1))
            pwant=alog(.01*p_phy(i,k,j))
            if(pwant.gt.pl)then
              backg_oh(i,k,j)=oh_backgd(l,j)
              backg_h2o2(i,k,j)=h2o2_backgd(l,j)
              backg_no3(i,k,j)=no3_backgd(l,j)
            else
              aln=(oh_backgd(l,j)*(pwant-pl)+            &
                  oh_backgd(l-1,j)*(pu-pwant))/(pu-pl)
               backg_oh(i,k,j)=aln
              aln=(h2o2_backgd(l,j)*(pwant-pl)+            &
                  h2o2_backgd(l-1,j)*(pu-pwant))/(pu-pl)
               backg_h2o2(i,k,j)=aln
              aln=(no3_backgd(l,j)*(pwant-pl)+            &
                  no3_backgd(l-1,j)*(pu-pwant))/(pu-pl)
               backg_no3(i,k,j)=aln
            endif
          enddo
       enddo
       enddo
       endif   ! end gocart stuff
       nv=1
       k=1
!      emis_ant=0.
!TBH       write(6,*)"airmw,mw_so4_aer,p_e_bc,p_e_oc=",airmw,mw_so4_aer,p_e_bc,p_e_oc
!TBH       write(6,*)"dtstep,rri(1,k,jmax),dz8w(1,k,jmax)=",dtstep,rri(1,k,jmax),dz8w(1,k,jmax)
       factor2=0.
       factor=0.
       if(p_bc2 .gt. 1)then
!!SMS$SERIAL BEGIN
!       write(0,*)'adding emissions for gocart',dtstep
!        write(0,*)'max ebc = ',maxval(emis_ant(:,:,:,p_e_bc))
!        write(0,*)'max eso2 = ',maxval(emis_ant(:,:,:,p_e_so2))
!        write(0,*)'max rri1 = ',maxval(rri(:,1,:))
!        write(0,*)'max dz1 = ',maxval(dz8w(:,1,:))
!!SMS$SERIAL END
       do i=its,ite
!TBH       write(6,*)"i = ",i,k
       do j=jts,jte
         factor=dtstep*rri(i,k,j)/dz8w(i,k,j)
         factor2=4.828e-4*dtstep*rri(i,k,j)/(60.*dz8w(i,k,j))
         chem(i,k,j,p_bc1)=chem(i,k,j,p_bc1)+emis_ant(i,k,j,p_e_bc)*factor
         chem(i,k,j,p_oc1)=chem(i,k,j,p_oc1)+emis_ant(i,k,j,p_e_oc)*factor
         chem(i,k,j,p_sulf)=chem(i,k,j,p_sulf)+emis_ant(i,k,j,p_e_sulf)*factor2
         chem(i,k,j,p_so2)=chem(i,k,j,p_so2)+emis_ant(i,k,j,p_e_so2)*factor2
       enddo
       enddo
       else if (p_tr2 .gt. 1)then    !co2 here
          do i=its,ite
          do j=jts,jte
!           factor2=dtstep*rri(i,k,j)/dz8w(i,k,j)
            factor2=4.828e-4*dtstep*rri(i,k,j)/(60.*dz8w(i,k,j))
            chem(i,k,j,p_tr1)=chem(i,k,j,p_tr1)+emis_ant(i,k,j,p_e_tr1)*factor2
            chem(i,k,j,p_tr2)=chem(i,k,j,p_tr2)+emis_ant(i,k,j,p_e_tr2)*factor2
          enddo
          enddo   
       else if ((p_tr2 .gt. 1) .and. (p_bc2 .gt. 1))then
          stop 'in chem_prep_fim, 112'
       endif
       do i=its,ite
       do j=jts,jte
       do nv=1,num_soil_layers
         smois(i,nv,j)=sm3d(nv,j)
       enddo
       enddo
       enddo
     curr_secs=ktau*ifix(dtstep)
     curr_hours=curr_secs/3600
!
!     do volcanoes if avaiable
!
!     if(chem_opt == 502 ) then
!       do j=jts,jte
!       if(emiss_ash_dt(j).le.0)CYCLE
!       emiss_ash_mass(j)=0.
!       emiss_ash_height(j)=0.
!       enddo
!
! default
!
        do j=jts,jte
        if(emiss_ash_dt(j).le.0)CYCLE
        so2_mass(j)=1.5e4*3600.*1.e9/64./area(j)
        eh=2600.*(emiss_ash_height(j)*.0005)**4.1494
        emiss_ash_mass(j)=eh*1.e9/area(j)
!       write(0,*)'h0 default ash mass = ',j,emiss_ash_dt(j),emiss_ash_height(j),emiss_ash_mass(j)
        enddo
! hard code for special retro case (set h1 - h6 properly
!
     if(curr_hours.ge.h1 .and. curr_hours.lt.h2)then
        do j=jts,jte
!       if(j.eq.jts)write(0,*)'curr_secs,curr_hours = ',curr_secs,curr_hours
        if(emiss_ash_dt(j).le.0)CYCLE
        emiss_ash_height(j)=5834.
!       eh=2600.*(emiss_ash_height(j)*.0005)**4.1494
        eh=3.11e5
        emiss_ash_mass(j)=eh*1.e9/area(j)
!       write(6,*)'h1 adjusted ash mass = ',emiss_ash_mass(j)
        enddo
     else if(curr_hours.ge.h2 .and. curr_hours.lt.h3)then
        do j=jts,jte
        if(emiss_ash_dt(j).le.0)CYCLE
        emiss_ash_height(j)=3834.
!       eh=2600.*(emiss_ash_height(j)*.0005)**4.1494
        eh=3.87e4
        emiss_ash_mass(j)=eh*1.e9/area(j)
!       write(6,*)'h2 adjusted ash mass = ',emiss_ash_mass(j)
        enddo
     else if(curr_hours.ge.h3 .and. curr_hours.lt.h4)then
        do j=jts,jte
        if(emiss_ash_dt(j).le.0)CYCLE
        emiss_ash_height(j)=5834.
!       eh=2600.*(emiss_ash_height(j)*.0005)**4.1494
        eh=3.11e5
        emiss_ash_mass(j)=eh*1.e9/area(j)
!       write(6,*)'h3 adjusted ash mass = ',emiss_ash_mass(j)
        enddo
     else if(curr_hours.ge.h4 .and. curr_hours.lt.h5)then
        do j=jts,jte
        if(emiss_ash_dt(j).le.0)CYCLE
        emiss_ash_height(j)=3334.
!       eh=2600.*(emiss_ash_height(j)*.0005)**4.1494
        eh=2.17e4
        emiss_ash_mass(j)=eh*1.e9/area(j)
!       write(6,*)'h4 adjusted ash mass = ',emiss_ash_mass(j)
        enddo
     else if(curr_hours.ge.h5 .and. curr_hours.lt.h6)then
        do j=jts,jte
        if(emiss_ash_dt(j).le.0)CYCLE
        emiss_ash_height(j)=3334.
!       eh=2600.*(emiss_ash_height(j)*.0005)**4.1494
        eh=2.17e4
        emiss_ash_mass(j)=eh*1.e9/area(j)
!       write(6,*)'h5 adjusted ash mass = ',emiss_ash_mass(j)
        enddo
     else if(curr_hours.ge.h6)then
        do j=jts,jte
        if(emiss_ash_dt(j).le.0)CYCLE
        emiss_ash_height(j)=2334.
        eh=4.93e3
        emiss_ash_mass(j)=eh*1.e9/area(j)
        enddo

     endif ! every o often
!     endif ! chem_opt = 502
! real-time application, keeping eruption constant
!
         if(ktau.le.2)then
            emis_vol(:,:,:,:)=0.
!      if(curr_hours.eq.h1 .or. curr_hours.eq.h2 .or. curr_hours.eq.h3 &
!         .or. curr_hours.eq.h4 .or. curr_hours.eq.h5 .or. curr_hours.eq.h6 .or. h1.gt.239)then
!         .or. curr_hours.eq.0)then
!         if(chem_opt == 316 .or. chem_opt == 317 .or. chem_opt == 502) then
! volcanic emissions
!
           do j=jts,jte
           do i=its,ite
             if(emiss_ash_dt(j).le.0)CYCLE
             if(emiss_ash_height(j).le.0.)CYCLE
             ashz_above_vent=emiss_ash_height(j) +z_at_w(i,kts,j)
!            write(0,*)'found and adjusted active volcano at j,kts,kpe = ',j,kts,kte
!            write(0,*)emiss_ash_height(j),emiss_ash_mass(j),emiss_ash_dt(j),ashz_above_vent
             do k=kte-1,kts,-1
                if(z_at_w(i,k,j) < ashz_above_vent)then
                  k_final=k+1
                  exit
                endif !inner
             enddo
             do k=kte-1,kts,-1
               if(z_at_w(i,k,j) < (1.-base_umbrel)*ashz_above_vent)then
                  k_initial=k
                  exit
                endif !inner
             enddo
             vert_mass_dist=0.
!              k_initial=int((k_final+k_initial)*0.5)
             
           !- parabolic vertical distribution between k_initial and k_final
               kk4 = k_final-k_initial+2
               do ko=1,kk4-1
                   kl=ko+k_initial-1
                   vert_mass_dist(kl) = 6.*percen_mass_umbrel* float(ko)/float(kk4)**2 * (1. - float(ko)/float(kk4))
               enddo
               if(sum(vert_mass_dist(kts:kte)) .ne. percen_mass_umbrel) then
                   x1= ( percen_mass_umbrel- sum(vert_mass_dist(kts:kte)) )/float(k_final-k_initial+1)
                   do ko=k_initial,k_final
                     vert_mass_dist(ko) = vert_mass_dist(ko)+ x1 !- values between 0 and 1.
                   enddo
                   !pause
               endif !inner
               !k_final > 0 .and. k_initial >
    
    !linear detrainment from vent to base of umbrella
               do ko=1,k_initial-1
                  vert_mass_dist(ko)=float(ko)/float(k_initial-1)
               enddo
               x1=sum(vert_mass_dist(1:k_initial-1))
    
               do ko=1,k_initial-1
                   vert_mass_dist(ko)=(1.-percen_mass_umbrel)*vert_mass_dist(ko)/x1
               enddo
               if(chem_opt == 316 ) then 
               do ko=1,k_final
                 emis_vol(i,ko,j,p_e_vash1)=.02*vert_mass_dist(ko)*emiss_ash_mass(j)
                 emis_vol(i,ko,j,p_e_vash2)=.04*vert_mass_dist(ko)*emiss_ash_mass(j)
                 emis_vol(i,ko,j,p_e_vash3)=.11*vert_mass_dist(ko)*emiss_ash_mass(j)
                 emis_vol(i,ko,j,p_e_vash4)=.09*vert_mass_dist(ko)*emiss_ash_mass(j)
                 emis_vol(i,ko,j,p_e_vash5)=.09*vert_mass_dist(ko)*emiss_ash_mass(j)
                 emis_vol(i,ko,j,p_e_vash6)=.13*vert_mass_dist(ko)*emiss_ash_mass(j)
                 emis_vol(i,ko,j,p_e_vash7)=.16*vert_mass_dist(ko)*emiss_ash_mass(j)
                 emis_vol(i,ko,j,p_e_vash8)=.16*vert_mass_dist(ko)*emiss_ash_mass(j)
                 emis_vol(i,ko,j,p_e_vash9)=.1*vert_mass_dist(ko)*emiss_ash_mass(j)
                 emis_vol(i,ko,j,p_e_vash10)=.1*vert_mass_dist(ko)*emiss_ash_mass(j)
               enddo
               do ko=k_final+1,kte
                 emis_vol(i,ko,j,p_e_vash1)=0.
                 emis_vol(i,ko,j,p_e_vash2)=0.
                 emis_vol(i,ko,j,p_e_vash3)=0.
                 emis_vol(i,ko,j,p_e_vash4)=0.
                 emis_vol(i,ko,j,p_e_vash5)=0.
                 emis_vol(i,ko,j,p_e_vash6)=0.
                 emis_vol(i,ko,j,p_e_vash7)=0.
                 emis_vol(i,ko,j,p_e_vash8)=0.
                 emis_vol(i,ko,j,p_e_vash9)=0.
                 emis_vol(i,ko,j,p_e_vash10)=0.
               enddo
               elseif (chem_opt == 317 .or. chem_opt == 502) then
!
! reduced vocanic ash transport
!
               do ko=1,k_final
                 emis_vol(i,ko,j,p_e_vash1)=.11*vert_mass_dist(ko)*emiss_ash_mass(j)
                 emis_vol(i,ko,j,p_e_vash2)=.08*vert_mass_dist(ko)*emiss_ash_mass(j)
                 emis_vol(i,ko,j,p_e_vash3)=.05*vert_mass_dist(ko)*emiss_ash_mass(j)
                 emis_vol(i,ko,j,p_e_vash4)=.035*vert_mass_dist(ko)*emiss_ash_mass(j)
               enddo
               elseif (chem_opt == 300) then
!
! if applied to gocart we only need finest ash bins, we use the coarse one for so2
!
               do ko=1,k_final
                 emis_vol(i,ko,j,p_e_vash1)=vert_mass_dist(ko)*so2_mass(j)
                 emis_vol(i,ko,j,p_e_vash2)=.08*vert_mass_dist(ko)*emiss_ash_mass(j)
                 emis_vol(i,ko,j,p_e_vash3)=.05*vert_mass_dist(ko)*emiss_ash_mass(j)
                 emis_vol(i,ko,j,p_e_vash4)=.035*vert_mass_dist(ko)*emiss_ash_mass(j)
               enddo
               endif !chem_opt==316 or 317,300,502
 
               do ko=k_final+1,kte
                 emis_vol(i,ko,j,p_e_vash1)=0.
                 emis_vol(i,ko,j,p_e_vash2)=0.
                 emis_vol(i,ko,j,p_e_vash3)=0.
                 emis_vol(i,ko,j,p_e_vash4)=0.
               enddo
       enddo
       enddo
!      endif ! chem_opt == 316 .or. chem_opt == 317 .or. chem_opt == 502
      endif ! curr_mins 
! initialy
!
       if(ktau.le.1)then

         if(chem_in_opt == 0 ) then
         if(chem_opt >= 300 .and. chem_opt <  500 ) then
           do j=jts,jte
              do k=kts,kte
                 do i=its,ite
                   do n=1,num_chem
                     chem(i,k,j,n)=1.e-12
                   enddo
                   chem(i,k,j,p_dms)=0.1e-6
                   chem(i,k,j,p_so2)=5.e-6
                   chem(i,k,j,p_sulf)=3.e-6
                   chem(i,k,j,p_msa)=0.1e-6
                   chem(i,k,j,p_bc1)=0.1e-3
                   chem(i,k,j,p_bc2)=0.1e-3
                   chem(i,k,j,p_oc1)=0.1e-3
                   chem(i,k,j,p_oc2)=0.1e-3
                   chem(i,k,j,p_p25)=1.
                   chem(i,k,j,p_p10)=1.
                 enddo
              enddo
           enddo
         endif ! chem_opt >= 300 .and. chem_opt <  500
         endif ! chem_in_opt == 0

!
! next is done to scale background oh and no3 in dependence on average zenith angle and day/night for no3
! this is done since background values are only available as average/month. It will not be necessary if other
! chemistry packages are used that provide oh,no3,h2o2
!
         if(chem_opt >= 300 .and. chem_opt <  500 ) then
           ndystep=86400/ifix(dtstepc)
           do i=its,ite
           do j=jts,jte
                   tcosz(i,j)=0.
                   ttday(i,j)=0.
                   rlat=xlat(i,j)*3.1415926535590/180.
                   xlonn=xlong(i,j)
!                 if(j.eq.681)then
!                   write(6,*)'szangle1',xlat(i,j),xlong(i,j)
!                   write(6,*)julday, sza, xlonn,rlat
!                   write(6,*)'ndystep,dtstepc,gmt = ',ndystep,dtstepc,gmt
!                 endif
                   do n=1,ndystep
                      xtime=n*dtstepc/60.
                      ixhour=ifix(gmt+.01)+ifix(xtime/60.)
                      xhour=float(ixhour)
                      xmin=60.*gmt+(xtime-xhour*60.)
                      gmtp=mod(xhour,24.)
                      gmtp=gmtp+xmin/60.
                      CALL szangle(1, 1, julday, gmtp, sza, cosszax,xlonn,rlat)
                      TCOSZ(i,j)=TCOSZ(I,J)+cosszax(1,1)
                      if(cosszax(1,1).gt.0.)ttday(i,j)=ttday(i,j)+dtstepc
                    enddo
!if (TCOSZ(i,j) == 0.0) then
!endif
!             if(j.eq.681)then
!                write(6,*)'szangle'
!                write(6,*)TCOSZ(i,j),ttday(i,j),julday, gmtp, sza, cosszax,xlonn,rlat
!             endif
           enddo
          enddo
         endif !chem_opt >= 300 .and. chem_opt <  500
         endif ! end ktau <=1
   if(CHEM_OPT == 316 )  then
       do i=its,ite
       do j=jts,jte
       if(emiss_ash_dt(j).le.0.)CYCLE
       do k=kts,kte-2
         factor2=dtstep*rri(i,k,j)/dz8w(i,k,j)
         chem(i,k,j,p_vash_1)=chem(i,k,j,p_vash_1)+emis_vol(i,k,j,p_e_vash1)*factor2
         chem(i,k,j,p_vash_2)=chem(i,k,j,p_vash_2)+emis_vol(i,k,j,p_e_vash2)*factor2
         chem(i,k,j,p_vash_3)=chem(i,k,j,p_vash_3)+emis_vol(i,k,j,p_e_vash3)*factor2
         chem(i,k,j,p_vash_4)=chem(i,k,j,p_vash_4)+emis_vol(i,k,j,p_e_vash4)*factor2
         chem(i,k,j,p_vash_5)=chem(i,k,j,p_vash_5)+emis_vol(i,k,j,p_e_vash5)*factor2
         chem(i,k,j,p_vash_6)=chem(i,k,j,p_vash_6)+emis_vol(i,k,j,p_e_vash6)*factor2
         chem(i,k,j,p_vash_7)=chem(i,k,j,p_vash_7)+emis_vol(i,k,j,p_e_vash7)*factor2
         chem(i,k,j,p_vash_8)=chem(i,k,j,p_vash_8)+emis_vol(i,k,j,p_e_vash8)*factor2
         chem(i,k,j,p_vash_9)=chem(i,k,j,p_vash_9)+emis_vol(i,k,j,p_e_vash9)*factor2
         chem(i,k,j,p_vash_10)=chem(i,k,j,p_vash_10)+emis_vol(i,k,j,p_e_vash10)*factor2
       enddo
       enddo
       enddo
     endif
   if(CHEM_OPT == 317 .or. CHEM_OPT == 502)  then
       do i=its,ite
       do j=jts,jte
       if(emiss_ash_dt(j).le.0.)CYCLE
!      print *,'adding volcanic ash emissions, dt = ',dtstep,rri(i,k_final,j),dz8w(i,k_final,j)
       do k=kts,kte-2
         factor2=dtstep*rri(i,k,j)/dz8w(i,k,j)
         chem(i,k,j,p_vash_1)=chem(i,k,j,p_vash_1)+emis_vol(i,k,j,p_e_vash1)*factor2
         chem(i,k,j,p_vash_2)=chem(i,k,j,p_vash_2)+emis_vol(i,k,j,p_e_vash2)*factor2
         chem(i,k,j,p_vash_3)=chem(i,k,j,p_vash_3)+emis_vol(i,k,j,p_e_vash3)*factor2
         chem(i,k,j,p_vash_4)=chem(i,k,j,p_vash_4)+emis_vol(i,k,j,p_e_vash4)*factor2
       enddo
       enddo
       enddo
     endif
   if(CHEM_OPT == 300 )  then
!
! for gocart only lump ash into p25 and p10
!
       do i=its,ite
       do j=jts,jte
       if(emiss_ash_dt(j).le.0.)CYCLE
       do k=kts,kte-2
         factor=4.828e-4*dtstep*rri(i,k,j)/(60.*dz8w(i,k,j))
         factor2=dtstep*rri(i,k,j)/dz8w(i,k,j)
         chem(i,k,j,p_p25)=chem(i,k,j,p_p25)                          &
                          +emis_vol(i,k,j,p_e_vash4)*factor2   
         chem(i,k,j,p_so2)=chem(i,k,j,p_so2)                          &
                          +emis_vol(i,k,j,p_e_vash1)*factor   
         chem(i,k,j,p_p10)=chem(i,k,j,p_p10)                          &
!                         +.5* emis_vol(i,k,j,p_e_vash4)*factor2      &  
                          +1.* emis_vol(i,k,j,p_e_vash3)*factor2      &  
                          +.5* emis_vol(i,k,j,p_e_vash2)*factor2
       enddo
       enddo
       enddo
     endif
!
! option 501 was only used for cesium ensemble - Japan 2010
!
         if(chem_opt == 501 ) then
! explosive tr emissions
!
           do j=jts,jte
           do i=its,ite
             if(emiss_tr_dt(j).le.0  .or. emiss_tr_height(j).le.0.)CYCLE
             ashz_above_vent=emiss_tr_height(j)+z_at_w(i,kts,j)
             do k=kte-1,kts,-1
                if(z_at_w(i,k,j) < ashz_above_vent)then
                  k_final=k+1
                  exit
                endif
             enddo
             do k=kte-1,kts,-1
               if(z_at_w(i,k,j) < (1.-base_umbrel)*ashz_above_vent)then
                  k_initial=k
                  exit
                endif
             enddo
             vert_mass_dist=0.
             k_initial=int((k_final+k_initial)*0.5)

           !- parabolic vertical distribution between k_initial and k_final
             kk4 = k_final-k_initial+2
             do ko=1,kk4-1
                   kl=ko+k_initial-1
                   vert_mass_dist(kl) = 6.*percen_mass_umbrel* float(ko)/float(kk4)**2 * (1. - float(ko)/float(kk4))
             enddo
             if(sum(vert_mass_dist(kts:kte)) .ne. percen_mass_umbrel) then
                   x1= ( percen_mass_umbrel- sum(vert_mass_dist(kts:kte)) )/float(k_final-k_initial+1)
                   do ko=k_initial,k_final
                     vert_mass_dist(ko) = vert_mass_dist(ko)+ x1 !- values between 0 and 1.
                   enddo
              endif

    !linear detrainment from vent to base of umbrella
               do ko=1,k_initial-1
                  vert_mass_dist(ko)=float(ko)/float(k_initial-1)
               enddo
               x1=sum(vert_mass_dist(1:k_initial-1))

               do ko=1,k_initial-1
                   vert_mass_dist(ko)=(1.-percen_mass_umbrel)*vert_mass_dist(ko)/x1
               enddo
! tr emissions for umbrella (explosive) type emissons
!
               do ko=1,k_final
                 emis_ant(i,ko,j,p_e_tr1)=vert_mass_dist(ko)*emiss_tr_mass(j)
                 emis_ant(i,ko,j,p_e_tr2)=1./float(k_final)*emiss_tr_mass(j)
               enddo
       if(emiss_tr_dt(j).le.0.)CYCLE
       do k=kts,kte-2
         factor2=dtstep*rri(i,k,j)/dz8w(i,k,j)
         chem(i,k,j,p_tr2)=chem(i,k,j,p_tr2)+emis_ant(i,k,j,p_e_tr2)*factor2
         if(real_time.gt.360.)chem(i,k,j,p_tr1)=chem(i,k,j,p_tr1)+emis_ant(i,k,j,p_e_tr2)*factor2
       enddo
       enddo 
       enddo       
       endif       



END  subroutine chem_prep_fim
END MODULE MODULE_CHEM_PREP_FIM
