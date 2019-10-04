
!NOTE:  For the moment, WRF physics is tightly coupled to FIM dynamics 
!NOTE:  via use-association with modules that declare dynamics variables.  

module module_wrfphysics

implicit none

integer :: ipn          ! Index for icos point number
integer :: itsP         ! Public version of its

contains
!*********************************************************************
!     wrf_physics
!	Calculates column forcing for global fim using physics routines from WRF
!	04/30/2009 - Georg Grell                - original version
!*********************************************************************

subroutine wrf_physics (ktau)

use module_control ,only: nvl,nip,nvlp1,dt,ntra,ntrb,yyyymmddhhmm
use module_wrf_control, only: num_soil_layers, num_moist, & 
                              ids,ide, jds,jde, kds,kde,            &
                              ims,ime, jms,jme, kms,kme,            &
                              its,ite, jts,jte, kts,kte
use module_constants ,only: cp, rd, grvity, qvmin, deg_lat, deg_lon
! move to module_sfc_variables
!USE module_fim_phy_init,only : zorl2d,vfrac2d,vtype2d,stype2d
! BEGIN dynamics variables...  
USE module_sfc_variables ,only : rn2d,rc2d,ts2d,us2d,hf2d,sw2d,     &
                                 slmsk2d,st3d,sm3d,zorl2d,vfrac2d,vtype2d,stype2d
use module_variables,only: us3d,vs3d,dp3d,pr3d,ph3d,ex3d,mp3d,diaga,diagb,      &
                           tr=>tr3d,trdp,ws3d,tk3d
! END dynamics variables
use module_wrf_variables,only: phys3dwrf,phys2dwrf,exch,pb2d

! TBH:  Ignore these so PPP doesn't have to translate them
!SMS$IGNORE BEGIN
USE module_initial_chem_namelists
USE module_wrfphysvars
USE module_microphysics_driver, only:microphysics_driver
USE module_cumulus_driver, only: cumulus_driver
USE module_wrfphys_prep_fim, only:wrfphys_prep_fim
!USE module_chemvars, only: chem
!SMS$IGNORE END

implicit none

! Dimension and type external variables:
integer,intent (IN)   :: ktau           ! model time step

integer :: ivl,current_month,current_gmt,julday,i,j,k,nv,nvv,stepcu,kpbl(ims:ime,jms:jme)
real :: dx,dy,gmt
real :: mu(ims:ime,jms:jme),edt_out(ims:ime,jms:jme),cutop(ims:ime,jms:jme),cubot(ims:ime,jms:jme)
real :: pr_ens(ims:ime,jms:jme,ensdim),xf_ens(ims:ime,jms:jme,ensdim)
   LOGICAL, DIMENSION( ims:ime , jms:jme ) ::                       CU_ACT_FLAG
real :: maxtmp

CHARACTER(len=9 )                :: jdate

!TODO:  Need better initial values for these?  
  cutop = 0.0
  cubot = 1000.

!sms$compare_var(st3d   , "begin wrf_physics ")
!sms$compare_var(sm3d   , "begin wrf_physics ")
!sms$compare_var(rn2d   , "begin wrf_physics ")
!sms$compare_var(rc2d   , "begin wrf_physics ")
!sms$compare_var(ts2d   , "begin wrf_physics ")
!sms$compare_var(us2d   , "begin wrf_physics ")
!sms$compare_var(hf2d   , "begin wrf_physics ")
!sms$compare_var(sw2d   , "begin wrf_physics ")
!sms$compare_var(slmsk2d, "begin wrf_physics ")

! dx,dy not needed as of now (used to be needed for microphys/etampnew)
    dx=0.
    dy=0.
    stepcu=1
    mu(:,:)=1.
     if(ktau.le.1)then
        READ(UNIT=yyyymmddhhmm(5:6), FMT='(I2)') current_month
     endif
     READ(UNIT=yyyymmddhhmm(9:10), FMT='(I2)') current_gmt
     call GetJdate(yyyymmddhhmm,jdate)                ! Julian date conversion
     READ(UNIT=jdate(3:5), FMT='(I3)')julday

     gmt=current_gmt

!----------------------------------------------------------------------
!    Loop begins over all horizontal grid points
!----------------------------------------------------------------------

!SMS$PARALLEL (dh,ipn) BEGIN
  do ipn=1,nip
    do ivl=1,nvl
      tr(ivl,ipn,2) = max(qvmin, tr(ivl,ipn,2))
    enddo
  enddo
!SMS$PARALLEL END
!    print *,'call wrfphys_prep'
     call wrfphys_prep_fim(ktau,dt,tr,tk3d,st3d,sm3d,dp3d,mp3d,ts2d,us2d,sw2d,pr3d,  &
             VFRAC2d,VTYPE2d,STYPE2d,us3d,vs3d,ws3d,slmsk2d,zorl2d,exch,pb2d,hf2d,&
             ex3d,pi_phy,gmt,julday,ph3d,deg_lat,deg_lon,nvl,nvlp1,ntra,ntrb,      &
             th_phy,rri,t_phy,moist,u_phy,v_phy,p_phy,tsk,grvity,rd,cp,&
             u10,v10,ivgtyp,isltyp,gsw,vegfra,rmol,ust,znt,xland,t8w,p8w,exch_h,pbl,hfx,ht,   &
             phys3dwrf,rqvblten,rqvften,rthraten,rthblten,rthften,                                      &
             xlat,xlong,z_at_w,zmid,dz8w,vvel,rho_phy,smois,num_soil_layers,num_moist,&
             ids,ide, jds,jde, kds,kde,                                      &
             ims,ime, jms,jme, kms,kme,                                      &
             its,ite,jts,jte,kts,kte)
!
! some of the optional stuff related to chem have for now been commented out
!
!    print *,'call cumulus driver'
   if(cu_physics.gt.0)then
   raincv(:,:)=0.
   rainc(:,:)=0.
   CALL cumulus_driver(U=u_phy,V=v_phy,TH=th_phy,T=t_phy,W=vvel,P=p_phy,PI=pi_phy,RHO=rho_phy,ITIMESTEP=ktau,DT=dt,&
       DX=dx,RAINC=rainc,RAINCV=raincv,HTOP=cutop,HBOT=cubot,KPBL=kpbl,DZ8W=dz8w,P8W=p8w,STEPCU=stepcu,&
       XLAND=xland,APR_GR=apr_gr,APR_W=apr_w,APR_MC=apr_mc,APR_ST=apr_st,APR_AS=apr_as,APR_CAPMA=apr_capma, &
       APR_CAPME=apr_capme,APR_CAPMI=apr_capmi,MASS_FLUX=mass_flux,XF_ENS=xf_ens,PR_ENS=pr_ens,HT=ht,EDT_OUT=edt_out,&
       imomentum=imomentum,clos_choice=clos_choice,cugd_tten=cugd_tten,cugd_qvten=cugd_qvten,cugd_qcten=cugd_qcten,  &
       cugd_ttens=cugd_ttens,cugd_qvtens=cugd_qvtens,ENSDIM=ensdim,MAXIENS=maxiens,MAXENS=maxens,MAXENS2=maxens2,    &
       MAXENS3=maxens3,CU_ACT_FLAG=cu_act_flag,GSW=gsw,cugd_avedx=cugd_avedx,CU_PHYSICS=cu_physics,                  &
       IDS=ids,IDE=ide,JDS=jds,JDE=jde,KDS=kds,KDE=kde,&
       IMS=ims,IME=ime,JMS=jms,JME=jme,KMS=kms,KME=kme,             &
       ITS=its,ITE=ite,JTS=jts,JTE=jte,KTS=kts,KTE=kte,RQVCUTEN=rqvcuten,RQCCUTEN=rqccuten,RQICUTEN=rqicuten,RQVBLTEN=rqvblten,&
       RQVFTEN=rqvften,RTHRATEN=rthraten,RTHBLTEN=rthblten,RTHCUTEN=rthcuten,RTHFTEN=rthften,                        &
       QV_CURR=moist(ims,kms,jms,P_QV),F_QV=F_QV,QC_CURR=moist(ims,kms,jms,P_QC),F_QC=F_QC,                          &
       QR_CURR=moist(ims,kms,jms,P_QR),F_QR=F_QR,QI_CURR=moist(ims,kms,jms,P_QI),F_QI=F_QI,                          &
       QS_CURR=moist(ims,kms,jms,P_QS),F_QS=F_QS,QG_CURR=moist(ims,kms,jms,P_QG),F_QG=F_QG,GD_CLOUD=GD_CLOUD,GD_CLOUD2=GD_CLOUD2)

      do j=jts,jte
      do k=kts,kte
      do i=its,ite
       tr(k,j,2)=tr(k,j,2)+rqvcuten(i,k,j)*dt
       diagb(k,j)=rqvcuten(i,k,j)
       diaga(k,j)=rthcuten(i,k,j)*pi_phy(i,k,j)
       tr(k,j,1)=tr(k,j,1)+rthcuten(i,k,j)*dt*(1.+.6078*tr(k,j,2))
       if(mp_physics.gt.0)then
          tr(k,j,3)=tr(k,j,3)+rqccuten(i,k,j)*dt
       else
          tr(k,j,3)=tr(k,j,3)+(rqicuten(i,k,j)+rqccuten(i,k,j))*dt
       endif
      enddo
      enddo
      enddo
     endif
     if(mp_physics.gt.0)then
     print *,'call microphys driver'
      maxtmp = maxval(moist(its:ite,kts:kte,jts:jte,p_qv))
!SMS$REDUCE(maxtmp,MAX)
      write(6,*)'1max qv = ',maxtmp
      maxtmp = maxval(moist(its:ite,kts:kte,jts:jte,p_qc))
!SMS$REDUCE(maxtmp,MAX)
      write(6,*)'1max qc = ',maxtmp
      maxtmp = maxval(moist(its:ite,kts:kte,jts:jte,p_qr))
!SMS$REDUCE(maxtmp,MAX)
      write(6,*)'1max qr = ',maxtmp
      maxtmp = maxval(moist(its:ite,kts:kte,jts:jte,p_qi))
!SMS$REDUCE(maxtmp,MAX)
      write(6,*)'1max qi = ',maxtmp
      maxtmp = maxval(moist(its:ite,kts:kte,jts:jte,p_qs))
!SMS$REDUCE(maxtmp,MAX)
      write(6,*)'1max qs = ',maxtmp
      maxtmp = minval(moist(its:ite,kts:kte,jts:jte,p_qs))
!SMS$REDUCE(maxtmp,MIN)
      write(6,*)'1min qs = ',maxtmp
     rainnc(:,:)=0.
     rainncv(:,:)=0.
     snownc(:,:)=0.
     snowncv(:,:)=0.
     graupelnc(:,:)=0.
     graupelncv(:,:)=0.
     CALL microphysics_driver(                                            &
      &         DT=dt             ,DX=dx              ,DY=dy   &
      &        ,DZ8W=dz8w                      &
      &        ,ITIMESTEP=ktau             &
      &        ,P8W=p8w            ,P=p_phy            ,PI_PHY=pi_phy     &  
      &        ,RHO=rho_phy                   &
      &        ,SR=sr              ,TH=th_phy                        &
      &        ,T8W=t8w                                                   &
      &        ,CLDFRA=cldfra, EXCH_H=exch_h &
!#ifdef WRF_CHEM
!      &        ,QLSINK=qlsink,CLDFRA_OLD=cldfra_old             &
!      &        ,PRECR=precr, PRECI=preci, PRECS=precs, PRECG=grid
!%precg &
!      &        ,CHEM_OPT=config_flags%chem_opt, PROGN=config_flags%progn  &
!#endif
      &        ,XLAND=xland                                          &
      &        ,MP_PHYSICS=config_flags%mp_physics                        &
      &        ,ID=1                                                 &
      &        ,IDS=ids,IDE=ide, JDS=jds,JDE=jde, KDS=kds,KDE=kde         &
      &        ,ITS=its,ITE=ite, JTS=jts,JTE=jte, KTS=kts,KTE=kte         &
      &        ,IMS=ims,IME=ime, JMS=jms,JME=jme, KMS=kms,KME=kme         &
                 ! Optional
      &        , RAINNC=rainnc, RAINNCV=rainncv                 &
      &        , SNOWNC=snownc, SNOWNCV=snowncv                 &
      &        , GRAUPELNC=graupelnc, GRAUPELNCV=graupelncv     &
      &        , W=vvel, Z=zmid, HT=ht                         &
      &        , QV_CURR=moist(ims,kms,jms,P_QV), F_QV=F_QV               &
      &        , QC_CURR=moist(ims,kms,jms,P_QC), F_QC=F_QC               &
      &        , QR_CURR=moist(ims,kms,jms,P_QR), F_QR=F_QR               &
      &        , QI_CURR=moist(ims,kms,jms,P_QI), F_QI=F_QI               &
      &        , QS_CURR=moist(ims,kms,jms,P_QS), F_QS=F_QS               &
      &        , QG_CURR=moist(ims,kms,jms,P_QG), F_QG=F_QG               &
      &        , QNDROP_CURR=scalar(ims,kms,jms,P_QNDROP), F_QNDROP=F_QNDROP &
      &        , QNI_CURR=scalar(ims,kms,jms,P_QNI), F_QNI=F_QNI          &
      &        , QT_CURR=scalar(ims,kms,jms,P_QT), F_QT=F_QT              &
      &        , QNS_CURR=scalar(ims,kms,jms,P_QNS), F_QNS=F_QNS          &
      &        , QNR_CURR=scalar(ims,kms,jms,P_QNR), F_QNR=F_QNR          &
      &        , QNG_CURR=scalar(ims,kms,jms,P_QNG), F_QNG=F_QNG          &
      &        , QNN_CURR=scalar(ims,kms,jms,P_QNN), F_QNN=F_QNN          &
      &        , QNC_CURR=scalar(ims,kms,jms,P_QNC), F_QNC=F_QNC          &
      &        , qrcuten=rqrcuten, qscuten=rqscuten             &
      &        , qicuten=rqicuten,mu=mu                        &
      &        , HAIL=config_flags%gsfcgce_hail                           & ! for gsfcgce
      &        , ICE2=config_flags%gsfcgce_2ice                           & ! for gsfcgce
                                                                          )

      do nv=1,num_moist-1
      nvv=4+nv
      do j=jts,jte
      do k=kts,kte
      do i=its,ite
       if(moist(i,k,j,nv+1).lt.1.e-15)moist(i,k,j,nv+1)=0.
       tr(k,j,nvv)=moist(i,k,j,nv+1)
      enddo
      enddo
      enddo
      enddo

      do j=jts,jte
      do k=kts,kte
      do i=its,ite
       tr(k,j,2)=max(1.e-15,moist(i,k,j,1))
       tr(k,j,3)=moist(i,k,j,p_qc)+moist(i,k,j,p_qr)+moist(i,k,j,p_qs)+moist(i,k,j,p_qi)
       tr(k,j,1)=th_phy(i,k,j)*(1.+0.6078*tr(k,j,2))
      enddo
      enddo
      enddo
      do nv=1,ntra+ntrb
      do j=jts,jte
      do k=kts,kte
      do i=its,ite
       trdp(k,j,nv)=tr(k,j,nv)*dp3d(k,j)
      enddo
      enddo
      enddo
      enddo
      maxtmp = maxval(moist(its:ite,kts:kte,jts:jte,p_qv))
!SMS$REDUCE(maxtmp,MAX)
      write(6,*)'max qv = ',maxtmp
      maxtmp = maxval(moist(its:ite,kts:kte,jts:jte,p_qc))
!SMS$REDUCE(maxtmp,MAX)
      write(6,*)'max qc = ',maxtmp
      maxtmp = maxval(moist(its:ite,kts:kte,jts:jte,p_qr))
!SMS$REDUCE(maxtmp,MAX)
      write(6,*)'max qr = ',maxtmp
      maxtmp = maxval(moist(its:ite,kts:kte,jts:jte,p_qi))
!SMS$REDUCE(maxtmp,MAX)
      write(6,*)'max qi = ',maxtmp
      maxtmp = maxval(moist(its:ite,kts:kte,jts:jte,p_qs))
!SMS$REDUCE(maxtmp,MAX)
      write(6,*)'max qs = ',maxtmp
      maxtmp = minval(moist(its:ite,kts:kte,jts:jte,p_qs))
!SMS$REDUCE(maxtmp,MIN)
      write(6,*)'min qs = ',maxtmp
      maxtmp = maxval(rainncv(its:ite,jts:jte))
!SMS$REDUCE(maxtmp,MAX)
      write(6,*)'max rainnc = ',maxtmp
      endif ! mp_physics
!      maxtmp = maxval(raincv(its:ite,jts:jte))
!!SMS$REDUCE(maxtmp,MAX)
!      write(6,*)'max rainc = ',maxtmp
!
! phys2dwrf (6:8) for feedback to gbphys
!
      if(cu_physics .ne. 0) then
         do j=jts,jte
         do i=its,ite
              phys2dwrf(j,5)=phys2dwrf(j,5)+raincv(i,j)
              phys2dwrf(j,6)=raincv(i,j)*.001
              phys2dwrf(j,7)=cubot(i,j)
              phys2dwrf(j,8)=cutop(i,j)
!             rc2d(j)=rc2d(j)+raincv(i,j) ! phys2dwrf(j,5)
!             if(raincv(i,j).gt.0)write(6,*)'in wrfphys', phys2dwrf(j,6),cubot(i,j),cutop(i,j)
         enddo
         enddo
      endif
      if(mp_physics .ne. 0) then
      do j=jts,jte
      do i=its,ite
       phys2dwrf(j,1)=rainncv(i,j)+phys2dwrf(j,1)
       phys2dwrf(j,2)=rainncv(i,j)
!      if(rainncv(i,j).gt.0.)write(6,*)i,j,rainncv(i,j),moist(i,1,j,p_qr)
       phys2dwrf(j,3)=snowncv(i,j)+phys2dwrf(j,3)
       phys2dwrf(j,4)=snowncv(i,j)
!      rn2d(j)=phys2dwrf(j,5)+phys2dwrf(j,1)+phys2dwrf(j,3)
!      if(mp_physis .ne. 0) if(raincv(i,j).gt.0.)write(6,*)i,j,raincv(i,j)
      enddo
      enddo
      endif
!      maxtmp = maxval(phys2dwrf(jts:jte,5))
!!SMS$REDUCE(maxtmp,MAX)
!      write(6,*)'max rainct = ',maxtmp

!sms$compare_var(st3d   , "end wrf_physics ")
!sms$compare_var(sm3d   , "end wrf_physics ")
!sms$compare_var(rn2d   , "end wrf_physics ")
!sms$compare_var(rc2d   , "end wrf_physics ")
!sms$compare_var(ts2d   , "end wrf_physics ")
!sms$compare_var(us2d   , "end wrf_physics ")
!sms$compare_var(hf2d   , "end wrf_physics ")
!sms$compare_var(sw2d   , "end wrf_physics ")
!sms$compare_var(slmsk2d, "end wrf_physics ")

return
end subroutine wrf_physics
end module module_wrfphysics
