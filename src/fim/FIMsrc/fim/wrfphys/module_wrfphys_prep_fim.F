MODULE MODULE_WRFPHYS_PREP_FIM
CONTAINS
     subroutine wrfphys_prep_fim(ktau,dtstep,tr3d,tk3d,st3d,sm3d,dp3d,mp3d,ts2d,us2d,sw2d,pr3d,  &
             VFRAC2d,VTYPE2d,STYPE2d,us3d,vs3d,ws3d,slmsk2d,zorl2d,exch,pb2d,hf2d,&
             ex3d,pi_phy,gmt,julday,ph3d,deg_lat,deg_lon,nvl,nvlp1,ntra,ntrb,      &
             th,rri,t_phy,moist,u_phy,v_phy,p_phy,tsk,g,rd,cp,&
             u10,v10,ivgtyp,isltyp,gsw,vegfra,rmol,ust,znt,xland,t8w,p8w,exch_h,pbl,hfx,ht,   &
             phys3dwrf,rqvblten,rqvften,rthraten,rthblten,rthften,                                      &
             xlat,xlong,z_at_w,zmid,dz8w,vvel,rho_phy,smois,num_soil_layers,num_moist,&
             ids,ide, jds,jde, kds,kde,                                      &
             ims,ime, jms,jme, kms,kme,                                      &
             its,ite,jts,jte,kts,kte)
!
! input fim variables
!
IMPLICIT NONE
INTEGER,      INTENT(IN   ) :: ktau,nvl,nvlp1,ntra,ntrb
REAL,         INTENT(IN   ) :: g,rd,dtstep,gmt,cp
real, intent(in) :: dp3d(nvl,jms:jme)      ! del p between coord levels (pascals)
real, intent(in) :: mp3d(nvl,jms:jme)      ! Montgomery Potential (m^2/s^2)
real, intent(in) :: tk3d(nvl,jms:jme)      ! temperature, kelvin
real, intent(in) :: exch(nvl,jms:jme)      !
real, intent(in) :: tr3d(nvl,jms:jme,ntra+ntrb)  ! 1=pot.temp, 2=water vapor, 3=cloud water, 4=ozone
real, intent(in) :: phys3dwrf(nvl,jms:jme,11)  ! 1=pot.temp, 2=water vapor, 3=cloud water, 4=ozone
real, intent(in) :: st3d(4,jms:jme)        ! soil temperature
real, intent(in) :: sm3d(4,jms:jme)        ! soil moisture
real, intent(in) :: ts2d(jms:jme)          ! skin temperature
real, intent(in) :: us2d(jms:jme)          ! friction velocity/equivalent momentum flux
real, intent(in) :: pb2d(jms:jme)          ! 
real, intent(in) :: hf2d(jms:jme)          ! 
real, intent(in) :: sw2d(jms:jme)          ! downward short-wave radiation flux
real, intent(in) :: pr3d(nvlp1,jms:jme)    ! pressure (pascal)
real, intent(in) :: ex3d(nvlp1,jms:jme)    ! exner function
real, intent(in) :: ph3d(nvlp1,jms:jme)    ! geopotential (=gz), m^2/s^2
real, dimension (jms:jme), intent(in) :: vfrac2d,VTYPE2d,STYPE2d,zorl2d,slmsk2d
real, dimension (nvl,jms:jme), intent(in) :: us3d,vs3d,ws3d
real, intent(in) :: deg_lat(jms:jme),deg_lon(jms:jme)  ! lat and lon in degrees

   INTEGER,      INTENT(IN   ) :: num_soil_layers,num_moist,julday,     &
                                  ids,ide, jds,jde, kds,kde,        &
                                  ims,ime, jms,jme, kms,kme,               &
                                  its,ite, jts,jte, kts,kte
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_moist ),                &
         INTENT(OUT ) ::                                   moist
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         ,               &
          INTENT(OUT   ) ::                                                 &
                                                        rri,               &
                                                      t_phy,th,            &
                                                      p_phy,pi_phy,        &
                                                      dz8w,p8w,t8w,        &
                        rqvblten,rqvften,rthraten,rthblten,rthften,        &
                                              z_at_w , zmid ,exch_h,       &
                                              u_phy,v_phy,vvel,rho_phy
   INTEGER,DIMENSION( ims:ime , jms:jme )                  ,               &
          INTENT(OUT   ) ::                                                 &
                                                     ivgtyp,               &
                                                     isltyp
   REAL,  DIMENSION( ims:ime , jms:jme )                   ,               &
          INTENT(OUT   ) ::                                                 &
                                                     u10,                  &
                                                     v10,                  &
                                                     gsw,                  &
                                                  vegfra,                  &
                                                     rmol,                 &
                                                     ust,                  &
                                                     xland,                &
                                                     xlat,      &
                                                     xlong,tsk,   &
                                                     znt,pbl,hfx,ht
  REAL, DIMENSION( ims:ime, num_soil_layers, jms:jme ) ,      &
      INTENT(OUT) ::                               smois
   integer i,j,k,kk,nv,jmax,jmaxi,l,ll,n,ndystep,ixhour
   real maxv,factor,factor2,pu,pl,aln,pwant,rlat
   real xhour,xmin,gmtp,xlonn,xtime
! .. Intrinsic Functions ..
      INTRINSIC max, min, float
!TBH       write(6,*)'in prep fim !!!!!!!!!!!!!!!!!!!'
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
        z_at_w(i,kts,j)=max(0.,ph3d(kts,j)/g)
        ht(i,j)=z_at_w(i,kts,j)
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
!      print *,maxval(vtype2d)
!      print *,maxval(stype2d)
!      print *,maxval(VFRAC2d)
       do i=its,ite
       do j=jts,jte
         pbl(i,j)=pb2d(j)
         hfx(i,j)=hf2d(j)
         xlat(i,j)=deg_lat(j)
         xlong(i,j)=deg_lon(j)
         ust(i,j)=us2d(j)
         tsk(i,j)=ts2d(j)
         ivgtyp(i,j)=VTYPE2d(j)
         isltyp(i,j)=STYPE2d(j)
         gsw(i,j)=sw2d(j)
         vegfra(i,j)=VFRAC2d(j)
!        if(j.eq.1000)write(6,*)'j1000 ',ivgtyp(i,j),isltyp(i,j),vegfra(i,j),slmsk2d(j)
!        if(ivgtyp(i,j).ne.0)write(6,*)i,j,ivgtyp(i,j),isltyp(i,j),vegfra(i,j),pb2d(j),VTYPE2d(j)
         rmol(i,j)=0.
         znt(i,j)=zorl2d(j)*.01
!SLMSK   - SEA(0),LAND(1),ICE(2) MASK
         xland(i,j)=1.
         if(slmsk2d(j).eq.0)xland(i,j)=0.
         if(slmsk2d(j).eq.1)xland(i,j)=1.
         if(slmsk2d(j).eq.2)xland(i,j)=2.
!        if (slmsk2d(j).gt.0.)write(6,*)j,slmsk2d(j)
         u10(i,j)=us3d(1,j)
         v10(i,j)=vs3d(1,j)
!        if(j.eq.1000)then
!          write(6,*)xlat(i,j),xlong(i,j),ph3d(1,j),xland(i,j)
!          write(6,*)vegfra(i,j),ust(i,j),tsk(i,j),ivgtyp(i,j)
!          write(6,*)isltyp(i,j),gsw(i,j),znt(i,j),z_at_w(i,kts,j)
!        endif
       enddo
       enddo
!       if(j.eq.3836)write(6,*)'in prep_fim !!!!!!!!!!!!!!!!!!!!!'
       do i=its,ite
       do k=kts,kte+1
       kk=min(k,kte)
       do j=jts,jte
        rqvblten(i,k,j)=phys3dwrf(kk,j,2)
        rthblten(i,k,j)=phys3dwrf(kk,j,5)
        rthraten(i,k,j)=phys3dwrf(kk,j,6)
        rthften(i,k,j)=phys3dwrf(kk,j,7)+rthblten(i,k,j)+rthraten(i,k,j)
        rqvften(i,k,j)=phys3dwrf(kk,j,3)+rqvblten(i,k,j)
!       t_phy(i,k,j)=tk3d(kk,j)
        th(i,k,j)=tr3d(kk,j,1)/(1.+.6078*tr3d(kk,j,2))
        zmid(i,k,j)=.5*(z_at_w(i,kk+1,j)+z_at_w(i,kk,j))
!       dz8w(i,k,j)=z_at_w(i,kk+1,j)-z_at_w(i,kk,j)
        p_phy(i,k,j)=.5*(p8w(i,kk,j)+p8w(i,kk+1,j))
        pi_phy(i,k,j)=.5*(ex3d(kk,j)+ex3d(kk+1,j))/cp
        t_phy(i,k,j)=th(i,k,j)*pi_phy(i,k,j)
        u_phy(i,k,j)=us3d(kk,j)
        exch_h(i,k,j)=exch(kk,j)
        v_phy(i,k,j)=vs3d(kk,j)
        rho_phy(i,k,j)= p_phy(i,k,j)/(RD*T_phy(i,k,j)*(1.+.608*tr3d(kk,j,2)))
        rri(i,k,j)=1./rho_phy(i,k,j)
        vvel(i,k,j)=-ws3d(kk,j)*rri(i,k,j)/g
!       if(j.eq.89)then
!       write(6,*)k,z_at_w(i,k,j),zmid(i,k,j),t_phy(i,k,j),p_phy(i,k,j)
!       endif
!       if(z_at_w(i,k,j).lt.0 .or. zmid(i,k,j).lt.0.) then
!        write(6,*)k,z_at_w(i,k,j),zmid(i,k,j),ph3d(kk,j),ph3d(kk+1,j)
!        write(6,*)i,j,t_phy(i,k,j),th(i,k,j),pi_phy(i,k,j),p_phy(i,k,j)
!        stop
!       endif
       enddo
       enddo
       enddo
       do i=its,ite
       do k=2,kte
       do j=jts,jte
        t8w(i,k,j)=.5*(t_phy(i,k,j)-t_phy(i,k-1,j)) ! .5*(tk3d(k-1,j)+tk3d(k,j))
       enddo
       enddo
       enddo
! do we know this here? do we need this?
       do i=its,ite
       do j=jts,jte
        t8w(i,1,j)=t_phy(i,1,j)
!       t8w(i,kte+1,j)=tk3d(kte,j)
       enddo
       enddo
! qv is in tr3d(2) !
!      write(6,*)'in wrfphysprep, kte,nvl,num_moist = ',kte,nvl,num_moist
       do nv=2,num_moist
       do i=its,ite
       do k=kts,kte+1
       kk=min(k,kte)
       do j=jts,jte
        moist(i,k,j,nv)=tr3d(kk,j,ntra+nv-1)
        if(moist(i,k,j,nv).lt.1.e-15)moist(i,k,j,nv)=0.
       enddo
       enddo
       enddo
       enddo
! p_qv is always 1!
       do i=its,ite
       do k=kts,kte+1
       kk=min(k,kte)
       do j=jts,jte
        moist(i,k,j,1)=tr3d(kk,j,2)
        if(j.eq.1000)then
!TBH        write(6,*)k,moist(i,k,j,p_qv),moist(i,k,j,p_qc),moist(i,k,j,p_qi)
        endif
       enddo
       enddo
       enddo
       do i=its,ite
       do j=jts,jte
       do nv=1,num_soil_layers
         smois(i,nv,j)=sm3d(nv,j)
       enddo
       enddo
       enddo
!TBH       write(6,*)'done prep fim !!!!!!!!!!!!!!!!!!!'
!      maxv=maxval(raincv_b)
!      write(6,*)maxv

END  subroutine wrfphys_prep_fim
END MODULE MODULE_WRFPHYS_PREP_FIM
