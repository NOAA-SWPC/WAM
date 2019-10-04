!WRF:MODEL_LAYER:PHYSICS
!

MODULE module_cu_gd

CONTAINS

!-------------------------------------------------------------
   SUBROUTINE GRELLDRV(                                         &
               DT,itimestep,DX                                  &
              ,rho,RAINCV,PRATEC                                &
              ,U,V,t,W,q,p,pi                                   &
              ,dz8w,p8w,XLV,CP,G,r_v                            &
              ,STEPCU,htop,hbot                                 &
              ,CU_ACT_FLAG                            &
              ,APR_GR,APR_W,APR_MC,APR_ST,APR_AS                &
              ,APR_CAPMA,APR_CAPME,APR_CAPMI                    &
              ,MASS_FLUX,XF_ENS,PR_ENS,HT,XLAND,gsw             &
              ,GDC,GDC2                                         &
              ,ensdim,maxiens,maxens,maxens2,maxens3            &
              ,ids,ide, jds,jde, kds,kde                        &
              ,ims,ime, jms,jme, kms,kme                        &
              ,its,ite, jts,jte, kts,kte                        &
              ,periodic_x,periodic_y                            &
              ,RQVCUTEN,RQCCUTEN,RQICUTEN                       &
              ,RQVFTEN,RQVBLTEN                                 &
              ,RTHFTEN,RTHCUTEN,RTHRATEN,RTHBLTEN               &
              ,F_QV    ,F_QC    ,F_QR    ,F_QI    ,F_QS         &
                                                                )
!-------------------------------------------------------------
   IMPLICIT NONE
!-------------------------------------------------------------
   INTEGER,      INTENT(IN   ) ::                               &
                                  ids,ide, jds,jde, kds,kde,    & 
                                  ims,ime, jms,jme, kms,kme,    & 
                                  its,ite, jts,jte, kts,kte
   LOGICAL periodic_x,periodic_y

   integer, intent (in   )              ::                      &
                       ensdim,maxiens,maxens,maxens2,maxens3
  
   INTEGER,      INTENT(IN   ) :: STEPCU, ITIMESTEP

   REAL,         INTENT(IN   ) :: XLV, R_v
   REAL,         INTENT(IN   ) :: CP,G

   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         ,    &
          INTENT(IN   ) ::                                      &
                                                          U,    &
                                                          V,    &
                                                          W,    &
                                                         pi,    &
                                                          t,    &
                                                          q,    &
                                                          p,    &
                                                       dz8w,    &
                                                       p8w,    &
                                                        rho
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         ,    &
          OPTIONAL                                         ,    &
          INTENT(INOUT   ) ::                                   &
               GDC,GDC2


   REAL, DIMENSION( ims:ime , jms:jme ),INTENT(IN) :: GSW,HT,XLAND
!
   REAL, INTENT(IN   ) :: DT, DX
!

   REAL, DIMENSION( ims:ime , jms:jme ),                        &
         INTENT(INOUT) ::         RAINCV, PRATEC, MASS_FLUX,    &
                          APR_GR,APR_W,APR_MC,APR_ST,APR_AS,    &
                              APR_CAPMA,APR_CAPME,APR_CAPMI,htop,hbot
!+lxz
!  REAL, DIMENSION( ims:ime , jms:jme ) :: & !, INTENT(INOUT) ::       &
!        HTOP,     &! highest model layer penetrated by cumulus since last reset in radiation_driver
!        HBOT       ! lowest  model layer penetrated by cumulus since last reset in radiation_driver
!                   ! HBOT>HTOP follow physics leveling convention

   LOGICAL, DIMENSION( ims:ime , jms:jme ),                     &
         INTENT(INOUT) ::                       CU_ACT_FLAG

!
! Optionals
!
   REAL, DIMENSION( ims:ime , kms:kme , jms:jme ),              &
         OPTIONAL,                                              &
         INTENT(INOUT) ::                           RTHFTEN,    &
                                                    RQVFTEN

   REAL, DIMENSION( ims:ime , kms:kme , jms:jme ),              &
         OPTIONAL,                                              &
         INTENT(IN   ) ::                                       &
                                                   RTHRATEN,    &
                                                   RTHBLTEN,    &
                                                   RQVBLTEN
   REAL, DIMENSION( ims:ime , kms:kme , jms:jme ),              &
         OPTIONAL,                                              &
         INTENT(INOUT) ::                                       &
                                                   RTHCUTEN,    &
                                                   RQVCUTEN,    &
                                                   RQCCUTEN,    &
                                                   RQICUTEN
!
! Flags relating to the optional tendency arrays declared above
! Models that carry the optional tendencies will provdide the
! optional arguments at compile time; these flags all the model
! to determine at run-time whether a particular tracer is in
! use or not.
!
   LOGICAL, OPTIONAL ::                                      &
                                                   F_QV      &
                                                  ,F_QC      &
                                                  ,F_QR      &
                                                  ,F_QI      &
                                                  ,F_QS



! LOCAL VARS
     real,    dimension ( ims:ime , jms:jme , 1:ensdim) ::      &
        massfln,xf_ens,pr_ens
     real,    dimension (its:ite,kts:kte+1) ::                    &
        OUTT,OUTQ,OUTQC,phh,cupclw
     real,    dimension (its:ite,kts:kte+1) ::  phf
     real,    dimension (its:ite)         ::                    &
        pret, ter11, aa0, fp
!+lxz
     integer, dimension (its:ite) ::                            &
        kbcon, ktop, ierr
!.lxz
     integer, dimension (its:ite,jts:jte) ::                    &
        iact_old_gr
     integer :: ichoice,iens,ibeg,iend,jbeg,jend,isum

!
! basic environmental input includes moisture convergence (mconv)
! omega (omeg), windspeed (us,vs), and a flag (aaeq) to turn off
! convection for this call only and at that particular gridpoint
!
     real,    dimension (its:ite,kts:kte+1) ::                    &
        T2d,TN,q2d,qo,PO,P2d,US,VS,omeg,pdot
     real, dimension (its:ite)            ::                    &
        Z1,PSUR,AAEQ,direction,mconv,cuten,umean,vmean,pmean

   INTEGER :: i,j,k,ICLDCK,ipr,jpr
   REAL    :: tcrit,dp,dq,dpm1,dpm0,dpp1,dpall
   INTEGER :: itf,jtf,ktf
   REAL    :: rkbcon,rktop        !-lxz

   ichoice=0
   iens=1
     ipr=1
     jpr=0
!    jpr=3077
!    ipr=0

      ibeg=its
      iend=ite
      jbeg=jts
      jend=jte

   tcrit=258.

   itf=ite
   ktf=kte
   jtf=jte
!  write(6,*)'ibeg,iend,jts,jtf,dt= ',ibeg,iend,jts,jtf,dt
!                                                                      
            isum=0
     DO 100 J = jts,jtf  
     DO I= its,itf
        ierr(i)=0.
        cuten(i)=0.
        iact_old_gr(i,j)=0
        mass_flux(i,j)=0.
        pratec(i,j) = 0.
        raincv(i,j)=0.
        CU_ACT_FLAG(i,j) = .true.
     ENDDO
     DO k=1,ensdim
     DO I= its,itf
        massfln(i,j,k)=0.
     ENDDO
     ENDDO
!j#if ( EM_CORE == 1 )
!    DO k= kts,ktf
!    DO I= its,itf
!       RTHFTEN(i,k,j)=(RTHFTEN(i,k,j)+RTHRATEN(i,k,j)+RTHBLTEN(i,k,j))*pi(i,k,j)
!       RQVFTEN(i,k,j)=RQVFTEN(i,k,j)+RQVBLTEN(i,k,j)
!    ENDDO
!    ENDDO
!    !  hydrostatic pressure, first on full levels
!    DO I=ITS,ITF
!        phf(i,1) = p8w(i,1,j)
!    ENDDO
!    !  integrate up, dp = -rho * g * dz
!    DO K=kts+1,ktf+1
!    DO I=ITS,ITF
!        phf(i,k) = phf(i,k-1) - rho(i,k-1,j) * g * dz8w(i,k-1,j)
!    ENDDO
!    ENDDO
!    !  scale factor so that pressure is not zero after integration
!    DO I=ITS,ITF
!        fp(i) = (p8w(i,kts,j)-p8w(i,kte,j))/(phf(i,kts)-phf(i,kte))
!    ENDDO
!    !  re-integrate up, dp = -rho * g * dz * scale_factor
!    DO K=kts+1,ktf+1
!    DO I=ITS,ITF
!        phf(i,k) = phf(i,k-1) - rho(i,k-1,j) * g * dz8w(i,k-1,j) * fp(i)
!    ENDDO
!    ENDDO
!    !  put hydrostatic pressure on half levels
!    DO K=kts,ktf
!    DO I=ITS,ITF
!        phh(i,k) = (phf(i,k) + phf(i,k+1))*0.5
!    ENDDO
!    ENDDO
!#endif
     DO I=ITS,ITF
!#if ( EM_CORE == 1 )
         PSUR(I)=p8w(I,1,J)*.01
!#endif
!#if ( NMM_CORE == 1 )
!         PSUR(I)=p(I,1,J)*.01
!#endif
         TER11(I)=HT(i,j)
         mconv(i)=0.
         aaeq(i)=0.
         direction(i)=0.
         pret(i)=0.
         umean(i)=0.
         vmean(i)=0.
         pmean(i)=0.
         ktop(i) = 0
         kbcon(i) = 1000
         htop(i,j)= 0.
         hbot(i,j)=1000.
     ENDDO
     DO K=kts,ktf
     DO I=ITS,ITF
         omeg(i,k)=0.
         cupclw(i,k)=0.
!#if ( EM_CORE == 1 )
!         po(i,k)=phh(i,k)*.01
!#endif

!#if ( NMM_CORE == 1 )
         po(i,k)=p(i,k,j)*.01
!#endif
         P2d(I,K)=PO(i,k)
         US(I,K) =u(i,k,j)
         VS(I,K) =v(i,k,j)
         T2d(I,K)=t(i,k,j)
         q2d(I,K)=q(i,k,j)
         omeg(I,K)= -g*rho(i,k,j)*w(i,k,j)
         pdot(i,k)=.01*omeg(i,k)
         TN(I,K)=t2d(i,k)+RTHFTEN(i,k,j)*dt
         IF(TN(I,K).LT.150.)TN(I,K)=T2d(I,K)
         QO(I,K)=q2d(i,k)+RQVFTEN(i,k,j)*dt
         IF(Q2d(I,K).LT.1.E-08)Q2d(I,K)=1.E-08
         IF(QO(I,K).LT.1.E-08)QO(I,K)=1.E-08
         OUTT(I,K)=0.
         OUTQ(I,K)=0.
         OUTQC(I,K)=0.
         if(j.eq.jpr)then
           write(6,*)k,p2d(i,k),t2d(i,k),tn(i,k),ht(i,j)
           write(6,*)'2',rqvften(i,k,j),q2d(i,k),QO(i,k)
           write(6,*)'3',rthblten(i,k,j),rthraten(i,k,j),RQVBLTEN(i,k,j)
         endif
!        RTHFTEN(i,k,j)=0.
!        RQVFTEN(i,k,j)=0.
     ENDDO
     ENDDO
      do k=  kts+1,ktf-1
      DO I = its,itf
         if((p2d(i,1)-p2d(i,k)).gt.150.and.p2d(i,k).gt.300)then
            dp=-.5*(p2d(i,k+1)-p2d(i,k-1))
            umean(i)=umean(i)+us(i,k)*dp
            vmean(i)=vmean(i)+vs(i,k)*dp
            pmean(i)=pmean(i)+dp
         endif
      enddo
      enddo
      DO I = its,itf
         umean(i)=umean(i)/pmean(i)
         vmean(i)=vmean(i)/pmean(i)
         direction(i)=(atan2(umean(i),vmean(i))+3.1415926)*57.29578
         if(direction(i).gt.360.)direction(i)=direction(i)-360.
      ENDDO
      DO K=kts,ktf-1
      DO I = its,itf
        dq=(q2d(i,k+1)-q2d(i,k))
        mconv(i)=mconv(i)+omeg(i,k)*dq/g
         if(j.eq.jpr)then
           write(6,*)'4',k,p2d(i,k),mconv(i),omeg(i,k)
         endif
      ENDDO
      ENDDO
      DO I = its,itf
        if(mconv(i).lt.0.)mconv(i)=0.
      ENDDO
!
!---- CALL CUMULUS PARAMETERIZATION
!
      CALL CUP_enss(pdot,outqc,j,AAEQ,T2d,Q2d,TER11,TN,QO,PO,PRET,     &
           P2d,OUTT,OUTQ,DT,PSUR,US,VS,tcrit,iens,                &
           mconv,massfln,iact_old_gr,omeg,direction,MASS_FLUX,    &
           maxiens,maxens,maxens2,maxens3,ensdim,                 &
           APR_GR,APR_W,APR_MC,APR_ST,APR_AS,                     &
           APR_CAPMA,APR_CAPME,APR_CAPMI,kbcon,ktop,              &
           ierr,xf_ens,pr_ens,XLAND,gsw,cupclw,                        &
           xlv,r_v,cp,g,ichoice,ipr,jpr,                          &
           ids,ide, jds,jde, kds,kde,                             &
           ims,ime, jms,jme, kms,kme,                             &
           its,ite, jts,jte, kts,kte                             )

      CALL neg_check(dt,q2d,outq,outt,outqc,pret,its,ite,kts,kte,itf,ktf)
      if(j.ge.jbeg.and.j.le.jend)then
            DO I=its,itf
!             cuten(i)=0.
              if(i.ge.ibeg.and.i.le.iend)then
!                write(6,*)'final pret ',ierr(i),j,pret(i),ht(i,j)
              if(pret(i).gt.0.)then
                 isum=isum+1
!                if(j.eq.jpr)write(6,*)'final pret ',i,j,pret(i)
                 pratec(i,j)=pret(i)
                 raincv(i,j)=pret(i)! *dt
                 cuten(i)=1.
                 rkbcon = kte+kts - kbcon(i)
                 rktop  = kte+kts -  ktop(i)
!                if (ktop(i)  > HTOP(i,j)) HTOP(i,j) = ktop(i)+.001
!                if (kbcon(i) < HBOT(i,j)) HBOT(i,j) = kbcon(i)+.001
                 HTOP(i,j) = ktop(i)+.001
                 HBOT(i,j) = kbcon(i)+.001
              endif
              else
               pret(i)=0.
              endif
            ENDDO
            K=kts
            DO K=kts,ktf-2
            DO I=its,itf
               RTHCUTEN(I,K,J)=outt(i,k)*cuten(i)/pi(i,k,j)
               RQVCUTEN(I,K,J)=outq(i,k)*cuten(i)
            ENDDO
            ENDDO
            DO K=kts+1,ktf-2
            DO I=its,itf
               dpm1=p8w(i,k-1,j)-p8w(i,k,j)
               dpm0=p8w(i,k,j)-p8w(i,k+1,j)
               dpp1=p8w(i,k+1,j)-p8w(i,k+2,j)
               dpall=.25*(p8w(i,k-1,j)-p8w(i,k+2,j))
!              RTHCUTEN(I,K,J)=(.25*dpm1*outt(i,k-1)/pi(i,k-1,j)   &
!                              +.5*dpm0*outt(i,k)/pi(i,k,j) &
!                              +.25*dpp1*outt(i,k+1)/pi(i,k+1,j))*cuten(i)/dpall
!              RQVCUTEN(I,K,J)=.25*(dpm1*outq(i,k-1)+2.*dpm0*outq(i,k)   &
!                              +dpp1*outq(i,k+1))*cuten(i)/dpall
               RTHCUTEN(I,K,J)=(.25*dpm1*outt(i,k-1)+.5*dpm0*outt(i,k) &
                               +.25*dpp1*outt(i,k+1))*cuten(i)/dpall/pi(i,k,j)
               RQVCUTEN(I,K,J)=.25*(dpm1*outq(i,k-1)+2.*dpm0*outq(i,k)   &
                               +dpp1*outq(i,k+1))*cuten(i)/dpall
!              RTHCUTEN(I,K,J)=(.25*outt(i,k-1)+.5*outt(i,k) &
!                              +.25*outt(i,k+1))*cuten(i)/pi(i,k,j)
!              RQVCUTEN(I,K,J)=.25*(outq(i,k-1)+2.*outq(i,k)   &
!                              +outq(i,k+1))*cuten(i)
!              if(j.eq.jpr)then
!               write(6,*)'dthdt,dqdt,dtdt = ',k,RTHCUTEN(I,K,J),RQVCUTEN(I,K,J),outt(i,k)
!              endif
            ENDDO
            ENDDO
!           if(j.eq.3077 .and. itimestep.eq.2)then
!           DO K=kts,ktf-1
!           DO I=its,itf
!               write(12,*)p2d(i,k),RTHCUTEN(I,K,J)*pi(i,k,j)
!               write(13,*)p2d(i,k),RQVCUTEN(I,K,J)
!           ENDDO
!           ENDDO
!           endif

            IF(PRESENT(RQCCUTEN)) THEN
              IF ( F_QC ) THEN
                DO K=kts,ktf
                DO I=its,itf
                   RQCCUTEN(I,K,J)=outqc(I,K)*cuten(i)
                   IF ( PRESENT( GDC ) ) GDC(I,K,J)=CUPCLW(I,K)*cuten(i)
                   IF ( PRESENT( GDC2 ) ) GDC2(I,K,J)=0.
                ENDDO
                ENDDO
              ENDIF
            ENDIF

!......     QSTEN STORES GRAUPEL TENDENCY IF IT EXISTS, OTHERISE SNOW (V2)     

            IF(PRESENT(RQICUTEN).AND.PRESENT(RQCCUTEN))THEN
              IF (F_QI) THEN
                DO K=kts,ktf
                DO I=its,itf
                   if(t2d(i,k).lt.258.)then
                      RQICUTEN(I,K,J)=outqc(I,K)*cuten(i)
                      RQCCUTEN(I,K,J)=0.
                      IF ( PRESENT( GDC2 ) ) GDC2(I,K,J)=CUPCLW(I,K)*cuten(i)
                   else
                      RQICUTEN(I,K,J)=0.
                      RQCCUTEN(I,K,J)=outqc(I,K)*cuten(i)
                      IF ( PRESENT( GDC ) ) GDC(I,K,J)=CUPCLW(I,K)*cuten(i)
                   endif
                ENDDO
                ENDDO
              ENDIF
            ENDIF
        endif !jbeg,jend

 100    continue
!                write(6,*)'!!!!!!!!!    isum = ',isum

   END SUBROUTINE GRELLDRV


   SUBROUTINE CUP_enss(pdot,OUTQC,J,AAEQ,T,Q,Z1,                            &
              TN,QO,PO,PRE,P,OUTT,OUTQ,DTIME,PSUR,US,VS,               &
              TCRIT,iens,mconv,massfln,iact,                           &
              omeg,direction,massflx,maxiens,                          &
              maxens,maxens2,maxens3,ensdim,                           &
              APR_GR,APR_W,APR_MC,APR_ST,APR_AS,                       &
              APR_CAPMA,APR_CAPME,APR_CAPMI,kbcon,ktop,                &   !-lxz
              ierr,xf_ens,pr_ens,xland,gsw,cupclw,                          &
              xl,rv,cp,g,ichoice,ipr,jpr,                              &
              ids,ide, jds,jde, kds,kde,                               &
              ims,ime, jms,jme, kms,kme,                               &
              its,ite, jts,jte, kts,kte                               )

   IMPLICIT NONE

      real :: cthk,cincrmax,cincrmin,dthk,val1,val2
      parameter(cthk=150.,cincrmax=180.,cincrmin=120.,dthk=25.)

     integer                                                           &
        ,intent (in   )                   ::                           &
        ids,ide, jds,jde, kds,kde,                                     &
        ims,ime, jms,jme, kms,kme,                                     &
        its,ite, jts,jte, kts,kte,ipr,jpr
     integer, intent (in   )              ::                           &
        j,ensdim,maxiens,maxens,maxens2,maxens3,ichoice,iens
  !
  ! 
  !
     real,    dimension (ims:ime,jms:jme,1:ensdim)                     &
        ,intent (inout)                   ::                           &
        massfln,xf_ens,pr_ens
     real,    dimension (ims:ime,jms:jme)                              &
        ,intent (inout )                  ::                           &
               APR_GR,APR_W,APR_MC,APR_ST,APR_AS,APR_CAPMA,     &
               APR_CAPME,APR_CAPMI,massflx
     real,    dimension (ims:ime,jms:jme)                              &
        ,intent (in   )                   ::                           &
               xland,gsw
     integer, dimension (its:ite,jts:jte)                              &
        ,intent (in   )                   ::                           &
        iact
  ! outtem = output temp tendency (per s)
  ! outq   = output q tendency (per s)
  ! outqc  = output qc tendency (per s)
  ! pre    = output precip
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
        OUTT,OUTQ,OUTQC,CUPCLW
     real,    dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        pre
!+lxz
     integer,    dimension (its:ite)                                   &
        ,intent (out  )                   ::                           &
        kbcon,ktop
!.lxz
  !
  ! basic environmental input includes moisture convergence (mconv)
  ! omega (omeg), windspeed (us,vs), and a flag (aaeq) to turn off
  ! convection for this call only and at that particular gridpoint
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        T,TN,PO,P,US,VS,omeg,pdot
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (inout)                   ::                           &
         Q,QO
     real, dimension (its:ite)                                         &
        ,intent (in   )                   ::                           &
        Z1,PSUR,AAEQ,direction,mconv

       
       real                                                            &
        ,intent (in   )                   ::                           &
        dtime,tcrit,xl,cp,rv,g


!
!  local ensemble dependent variables in this routine
!
     real,    dimension (its:ite,1:maxens)  ::                         &
        xaa0_ens
     real,    dimension (1:maxens)  ::                                 &
        mbdt_ens
     real,    dimension (1:maxens2) ::                                 &
        edt_ens
     real,    dimension (its:ite,1:maxens2) ::                         &
        edtc
     real,    dimension (its:ite,kts:kte,1:maxens2) ::                 &
        dellat_ens,dellaqc_ens,dellaq_ens,pwo_ens
!
!
!
!***************** the following are your basic environmental
!                  variables. They carry a "_cup" if they are
!                  on model cloud levels (staggered). They carry
!                  an "o"-ending (z becomes zo), if they are the forced
!                  variables. They are preceded by x (z becomes xz)
!                  to indicate modification by some typ of cloud
!
  ! z           = heights of model levels
  ! q           = environmental mixing ratio
  ! qes         = environmental saturation mixing ratio
  ! t           = environmental temp
  ! p           = environmental pressure
  ! he          = environmental moist static energy
  ! hes         = environmental saturation moist static energy
  ! z_cup       = heights of model cloud levels
  ! q_cup       = environmental q on model cloud levels
  ! qes_cup     = saturation q on model cloud levels
  ! t_cup       = temperature (Kelvin) on model cloud levels
  ! p_cup       = environmental pressure
  ! he_cup = moist static energy on model cloud levels
  ! hes_cup = saturation moist static energy on model cloud levels
  ! gamma_cup = gamma on model cloud levels
!
!
  ! hcd = moist static energy in downdraft
  ! zd normalized downdraft mass flux
  ! dby = buoancy term
  ! entr = entrainment rate
  ! zd   = downdraft normalized mass flux
  ! entr= entrainment rate
  ! hcd = h in model cloud
  ! bu = buoancy term
  ! zd = normalized downdraft mass flux
  ! gamma_cup = gamma on model cloud levels
  ! mentr_rate = entrainment rate
  ! qcd = cloud q (including liquid water) after entrainment
  ! qrch = saturation q in cloud
  ! pwd = evaporate at that level
  ! pwev = total normalized integrated evaoprate (I2)
  ! entr= entrainment rate
  ! z1 = terrain elevation
  ! entr = downdraft entrainment rate
  ! jmin = downdraft originating level
  ! kdet = level above ground where downdraft start detraining
  ! psur        = surface pressure
  ! z1          = terrain elevation
  ! pr_ens = precipitation ensemble
  ! xf_ens = mass flux ensembles
  ! massfln = downdraft mass flux ensembles used in next timestep
  ! omeg = omega from large scale model
  ! mconv = moisture convergence from large scale model
  ! zd      = downdraft normalized mass flux
  ! zu      = updraft normalized mass flux
  ! dir     = "storm motion"
  ! mbdt    = arbitrary numerical parameter
  ! dtime   = dt over which forcing is applied
  ! iact_gr_old = flag to tell where convection was active
  ! kbcon       = LFC of parcel from k22
  ! k22         = updraft originating level
  ! icoic       = flag if only want one closure (usually set to zero!)
  ! dby = buoancy term
  ! ktop = cloud top (output)
  ! xmb    = total base mass flux
  ! hc = cloud moist static energy
  ! hkb = moist static energy at originating level
  ! mentr_rate = entrainment rate

     real,    dimension (its:ite,kts:kte) ::                           &
        he,hes,qes,z,                                                  &
        heo,heso,qeso,zo,                                              &
        xhe,xhes,xqes,xz,xt,xq,                                        &

        qes_cup,q_cup,he_cup,hes_cup,z_cup,p_cup,gamma_cup,t_cup,      &
        qeso_cup,qo_cup,heo_cup,heso_cup,zo_cup,po_cup,gammao_cup,     &
        tn_cup,                                                        &
        xqes_cup,xq_cup,xhe_cup,xhes_cup,xz_cup,xp_cup,xgamma_cup,     &
        xt_cup,                                                        &

        dby,qc,qrcd,pwd,pw,hcd,qcd,dbyd,hc,qrc,zu,zd,clw_all,          &
        dbyo,qco,qrcdo,pwdo,pwo,hcdo,qcdo,dbydo,hco,qrco,zuo,zdo,      &
        xdby,xqc,xqrcd,xpwd,xpw,xhcd,xqcd,xhc,xqrc,xzu,xzd,            &

  ! cd  = detrainment function for updraft
  ! cdd = detrainment function for downdraft
  ! dellat = change of temperature per unit mass flux of cloud ensemble
  ! dellaq = change of q per unit mass flux of cloud ensemble
  ! dellaqc = change of qc per unit mass flux of cloud ensemble

        cd,cdd,scr1,DELLAH,DELLAQ,DELLAT,DELLAQC

  ! aa0 cloud work function for downdraft
  ! edt = epsilon
  ! aa0     = cloud work function without forcing effects
  ! aa1     = cloud work function with forcing effects
  ! xaa0    = cloud work function with cloud effects (ensemble dependent)
  ! edt     = epsilon
     real,    dimension (its:ite) ::                                   &
       edt,edto,edtx,AA1,AA0,XAA0,HKB,HKBO,aad,XHKB,QKB,QKBO,          &
       XMB,XPWAV,XPWEV,PWAV,PWEV,PWAVO,PWEVO,BU,BUO,cap_max,xland1,     &
       cap_max_increment,closure_n,xmbmax,pbcdif
     integer,    dimension (its:ite) ::                                &
       kzdown,KDET,K22,KB,JMIN,kstabi,kstabm,K22x,                     &   !-lxz
       KBCONx,KBx,KTOPx,ierr,ierr2,ierr3,KBMAX 

     integer                              ::                           &
       nall,iedt,nens,nens3,ki,I,K,KK,iresult
     real                                 ::                           &
      day,dz,mbdt,entr_rate,radius,entrd_rate,mentr_rate,mentrd_rate,  &
      zcutdown,edtmax,edtmin,depth_min,zkbmax,z_detr,zktop,            &
      massfld,dh,cap_maxs
! thorpex 20110506 - mods to implement SAS downdraft detrainment
     real                                 ::           &
       xlamdd,xlamde,beta,dp
! thorpex 20110511 - mods to implement SAS updraft (en/de)trainment
     real                                 ::           &
       clam,clam2,clam3,cxlamu,cxlamu2,cxlamu3,                        &
       fent1,fent2,pgcon,tem,tem1,frh,w1,w1l,w1s,w2,w2l,&
       w2s,w3,w3l,w3s,w4,w4l,w4s,dzmax
     real,    dimension (its:ite,kts:kte) ::           &
       xlamue,xlamue_sflx,xlamue3
     real,    dimension (its:ite) ::           &
       cincr,xlamud,xlamud_sflx,xlamud3


     integer :: itf,jtf,ktf,jprt
     integer :: jmini
     logical :: keep_going
     logical,    dimension (its:ite) :: flg
      w1l     = -8.e-3
      w2l     = -4.e-2
      w3l     = -5.e-3
      w4l     = -5.e-4
      w1s     = -2.e-4
      w2s     = -2.e-3
      w3s     = -1.e-3
      w4s     = -2.e-5
      w1=w1l
      w2=w2l
      w3=w3l
      w4=w4l
      pbcdif(:) = 0.
      cincr(:) = 0.


     itf=ite
     ktf=kte
     jtf=jte

!sms$distribute end
      day=86400.
!TBH      if(j.eq.jpr)write(6,*)'in cup_enss'
      do i=its,itf
        closure_n(i)=16.
        xland1(i)=1.
        if(xland(i,j).gt.1.5)xland1(i)=0.
        cap_max_increment(i)=25.
      enddo
!
!--- specify entrainmentrate and detrainmentrate
!
      if(iens.le.4)then
      radius=14000.-float(iens)*2000.
      else
      radius=12000.
      endif
!
!--- gross entrainment rate (these may be changed later on in the
!--- program, depending what your detrainment is!!)
!
      entr_rate=.2/radius

!
!--- entrainment of mass
!
! thorpex 20110506 - adopted from SAS 
!   entrainment coefficients
      cxlamu  = 1.0e-4
      cxlamu2  = 1.0e-4
      clam    = .1
      clam2    = .1
      clam3    = .1
      pgcon = 0.55       ! Zhang & Wu (2003,JAS)
!   detrainment coefficients
      xlamde  = 1.0e-4
      xlamdd  = 1.0e-4
      beta   = .05
      mentrd_rate=xlamde
      mentr_rate=entr_rate
      xlamue(:,:) = 0.
      xlamud(:) = 0.

!
!--- initial detrainmentrates
!
      do k=kts,ktf
      do i=its,itf
        cupclw(i,k)=0.
        cd(i,k)=0.1*entr_rate
        cdd(i,k)=0.
      enddo
      enddo
!
!--- max/min allowed value for epsilon (ratio downdraft base mass flux/updraft
!    base mass flux
!
      edtmax=.3
      edtmin=.0
!
!--- minimum depth (m), clouds must have
!
      depth_min=500.
!
!--- maximum depth (mb) of capping 
!--- inversion (larger cap = no convection)
!
      cap_maxs=100.
!sms$to_local(grid_dh: <1, mix :size>, <2, mjx :size>) begin
      DO 7 i=its,itf
        kbmax(i)=1
        aa0(i)=0.
        aa1(i)=0.
        aad(i)=0.
        edt(i)=0.
        kstabm(i)=ktf-1
        IERR(i)=0
        IERR2(i)=0
        IERR3(i)=0
        if(aaeq(i).lt.-1.)then
           ierr(i)=20
        endif
 7    CONTINUE
!
!--- first check for upstream convection
!
      do i=its,itf
          cap_max(i)=cap_maxs
!         if(tkmax(i,j).lt.2.)cap_max(i)=25.
!         if(gsw(i,j).lt.1.)cap_max(i)=25.

        iresult=0
!       massfld=0.
!     call cup_direction2(i,j,direction,iact, &
!          cu_mfx,iresult,0,massfld,  &
!          ids,ide, jds,jde, kds,kde, &
!          ims,ime, jms,jme, kms,kme, &
!          its,ite, jts,jte, kts,kte)
!         cap_max(i)=cap_maxs
       if(iresult.eq.1)then
          cap_max(i)=cap_maxs+20.
       endif
!      endif
      enddo
!
!--- max height(m) above ground where updraft air can originate
!
      zkbmax=4000.
!
!--- height(m) above which no downdrafts are allowed to originate
!
      zcutdown=3000.
!
!--- depth(m) over which downdraft detrains all its mass
!
      z_detr=1250.
!
      do nens=1,maxens
         mbdt_ens(nens)=(float(nens)-3.)*dtime*1.e-3+dtime*5.E-03
      enddo
      do nens=1,maxens2
         edt_ens(nens)=.95-float(nens)*.01
      enddo
!     if(j.eq.jpr)then
!       print *,'radius ensemble ',iens,radius
!       print *,mbdt_ens
!       print *,edt_ens
!     endif
!
!--- environmental conditions, FIRST HEIGHTS
!
      do i=its,itf
         if(ierr(i).ne.20)then
            do k=1,maxens*maxens2*maxens3
               xf_ens(i,j,(iens-1)*maxens*maxens2*maxens3+k)=0.
               pr_ens(i,j,(iens-1)*maxens*maxens2*maxens3+k)=0.
            enddo
         endif
      enddo
!
!--- calculate moist static energy, heights, qes
!
      call cup_env(z,qes,he,hes,t,q,p,z1, &
           psur,ierr,tcrit,0,xl,cp,   &
           ids,ide, jds,jde, kds,kde, &
           ims,ime, jms,jme, kms,kme, &
           its,ite, jts,jte, kts,kte)
      call cup_env(zo,qeso,heo,heso,tn,qo,po,z1, &
           psur,ierr,tcrit,0,xl,cp,   &
           ids,ide, jds,jde, kds,kde, &
           ims,ime, jms,jme, kms,kme, &
           its,ite, jts,jte, kts,kte)
!
!--- environmental values on cloud levels
!
      call cup_env_clev(t,qes,q,he,hes,z,p,qes_cup,q_cup,he_cup, &
           hes_cup,z_cup,p_cup,gamma_cup,t_cup,psur, &
           ierr,z1,xl,rv,cp,          &
           ids,ide, jds,jde, kds,kde, &
           ims,ime, jms,jme, kms,kme, &
           its,ite, jts,jte, kts,kte)
      call cup_env_clev(tn,qeso,qo,heo,heso,zo,po,qeso_cup,qo_cup, &
           heo_cup,heso_cup,zo_cup,po_cup,gammao_cup,tn_cup,psur,  &
           ierr,z1,xl,rv,cp,          &
           ids,ide, jds,jde, kds,kde, &
           ims,ime, jms,jme, kms,kme, &
           its,ite, jts,jte, kts,kte)
      do i=its,itf
      if(ierr(i).eq.0)then
!
      do k=kts,ktf-2
        if(zo_cup(i,k).gt.zkbmax+z1(i))then
          kbmax(i)=k
          go to 25
        endif
      enddo
 25   continue
!
!
!--- level where detrainment for downdraft starts
!
      do k=kts,ktf
        if(zo_cup(i,k).gt.z_detr+z1(i))then
          kdet(i)=k
          go to 26
        endif
      enddo
 26   continue
!
      endif
      enddo
!
!
!
!------- DETERMINE LEVEL WITH HIGHEST MOIST STATIC ENERGY CONTENT - K22
!-------- start with level2, first cup level is ground
!
      CALL cup_MAXIMI(HEO_CUP,2,KBMAX,K22,ierr, &
           ids,ide, jds,jde, kds,kde, &
           ims,ime, jms,jme, kms,kme, &
           its,ite, jts,jte, kts,kte)
       DO 36 i=its,itf
         IF(ierr(I).eq.0.)THEN
         IF(K22(I).GE.KBMAX(i))ierr(i)=2
         endif
 36   CONTINUE
!           
!  look for the level of free convection as cloud base -from SAS option
!           
      do i=its,itf
        kbcon(i) = ktf
        flg(i)=.true.
      enddo 
      do k = 2, ktf-1
        do i=its,itf
          if (flg(i) .and. ierr(i).eq.0 .and. k.le.kbmax(i)) then
            if(k.gt.k22(i).and.heo_cup(i,k22(i)).gt.heso_cup(i,k)) then
              kbcon(i) = k
              flg(i)   = .false.
            endif
          endif
        enddo
      enddo 
      do i=its,itf
        if(ierr(i).eq.0)then
           if(kbcon(i).eq.ktf) ierr(i) = 3
        endif
      enddo 

!
!--- DETERMINE THE LEVEL OF CONVECTIVE CLOUD BASE  - KBCON
!

! determine critical convective inhibition
!  as a function of vertical velocity at cloud base.
!
      do i=its,itf
        if(ierr(i).eq.0) then
          if(xland1(i).eq.0)then
             w1=w1s
             w2=w2s
             w3=w3s
             w4=w4s
          endif
          if(pdot(i,kbcon(i)).le.w4) then
            tem = (pdot(i,kbcon(i)) - w4) / (w3 - w4)
          elseif(pdot(i,kbcon(i)).ge.-w4) then
            tem = - (pdot(i,kbcon(i)) + w4) / (w4 - w3)
          else
            tem = 0.
          endif
          val1    =             -1.
          tem = max(tem,val1)
          val2    =             1.
          tem = min(tem,val2)
          tem = 1. - tem
          tem1= .5*(cincrmax-cincrmin)
          cincr(i) = cincrmax - tem * tem1
          pbcdif(i) = p_cup(i,k22(i)) - p_cup(i,kbcon(i))
          if(pbcdif(i).gt.cincr(i)) then
             ierr(i) = 7
          endif
!         if(pbcdif(i).gt.cincr(i)-25.) then
!            ierr2(i) = 7
!         endif
!         if(pbcdif(i).gt.cincr(i)-50.) then
!            ierr3(i) = 7
!         endif
        endif
      enddo

!     call cup_kbcon(cap_max_increment,1,k22,kbcon,heo_cup,heso_cup, &
!          ierr,kbmax,po_cup,cap_max, &
!          ids,ide, jds,jde, kds,kde, &
!          ims,ime, jms,jme, kms,kme, &
!          its,ite, jts,jte, kts,kte)
!     call cup_kbcon_cin(1,k22,kbcon,heo_cup,heso_cup,z,tn_cup, &
!          qeso_cup,ierr,kbmax,po_cup,cap_max,xl,cp,&
!          ids,ide, jds,jde, kds,kde, &
!          ims,ime, jms,jme, kms,kme, &
!          its,ite, jts,jte, kts,kte)
!
!--- increase detrainment in stable layers
!
!     CALL cup_minimi(HEso_cup,Kbcon,kstabm,kstabi,ierr,  &
!          ids,ide, jds,jde, kds,kde, &
!          ims,ime, jms,jme, kms,kme, &
!          its,ite, jts,jte, kts,kte)
! thorpex2 (SASRQ?) - Adaptation of SAS updraft entrainment
! Updraft entrainment rate is initially set to an inverse function of height 
      do i=its,itf
        do k=kts,ktf-1
          xlamue(i,k) = clam / z_cup(i,k)
        enddo
      enddo
! Assume updraft entrainment rate above cloud base is same as that
! at cloud base.
! Assume updraft detrainment rate equals updraft entrainment rate 
! at cloud base.
      do i=its,itf
        if (ierr(i).eq.0) then
          xlamud(i) = xlamue(i,kbcon(i))
          do k=kbcon(i)+1,ktf-1
            xlamue(i,k) = xlamue(i,kbcon(i))
          enddo
        endif
      enddo
!
!  final entrainment rate as the sum of turbulent part and organized entrainment
!  depending on the environmental relative humidity
!  functions rapidly decreasing with height, mimicking a cloud ensemble
!    (Bechtold et al., 2008)

      do i=its,itf
        if (ierr(i).eq.0) then
          dzmax=0.
          do k=kbcon(i),ktf-1
            dz=zo_cup(i,k)-zo_cup(i,k-1)
            dzmax=max(dzmax,dz)
            tem = qeso_cup(i,k)/qeso_cup(i,kbcon(i))
            frh = 1.-min(qo_cup(i,k)/qeso_cup(i,k),1.)
            fent1 = tem**2
            fent2 = tem**3
            xlamue(i,k) = xlamue(i,k)*fent1 + cxlamu*frh*fent2
            if(j.eq.jpr)print *,'k',k,fent1,fent2,xlamue(i,k),cxlamu*frh*fent2
          enddo
          tem = xlamud(i)-1./dzmax
          do k=kts,ktf-1
            xlamue(i,k)=max(xlamue(i,k),tem)
          enddo
        endif
      enddo

      do i=its,itf
        if (ierr(i).eq.0) then
          do k=kts+1,ktf-1
            cd(i,k)=xlamud(i)
          enddo
          cd(i,1)=xlamud(i)
        endif
      enddo

!     do i=its,itf
!     IF(ierr(I).eq.0.)THEN
!       if(kstabm(i)-1.gt.kstabi(i))then
!          do k=kstabi(i),kstabm(i)-1
!            cd(i,k)=cd(i,k-1)+1.5*entr_rate
!            if(cd(i,k).gt.10.0*entr_rate)cd(i,k)=10.0*entr_rate
!          enddo
!       ENDIF
!     ENDIF
!     ENDDO
!
!--- calculate incloud moist static energy
!
      call cup_up_he(k22,hkb,z_cup,cd,mentr_rate,he_cup,hc, &
           kbcon,ierr,dby,he,hes_cup, &
           xlamue,xlamud,             &
           ids,ide, jds,jde, kds,kde, &
           ims,ime, jms,jme, kms,kme, &
           its,ite, jts,jte, kts,kte)
      call cup_up_he(k22,hkbo,zo_cup,cd,mentr_rate,heo_cup,hco, &
           kbcon,ierr,dbyo,heo,heso_cup, &
           xlamue,xlamud,             &
           ids,ide, jds,jde, kds,kde, &
           ims,ime, jms,jme, kms,kme, &
           its,ite, jts,jte, kts,kte)

!--- DETERMINE CLOUD TOP - KTOP
!
      call cup_ktop(1,dbyo,kbcon,ktop,ierr, &
           ids,ide, jds,jde, kds,kde, &
           ims,ime, jms,jme, kms,kme, &
           its,ite, jts,jte, kts,kte)
      DO 37 i=its,itf
         kzdown(i)=0
         if(ierr(i).eq.0)then
            zktop=(zo_cup(i,ktop(i))-z1(i))*.6
            zktop=min(zktop+z1(i),zcutdown+z1(i))
            do k=kts,kte
              if(zo_cup(i,k).gt.zktop)then
                 kzdown(i)=k
                 go to 37
              endif
              enddo
         endif
 37   CONTINUE
!
!--- DOWNDRAFT ORIGINATING LEVEL - JMIN
!
      call cup_minimi(HEo_cup,Kbcon,kbmax,JMIN,ierr, &
           ids,ide, jds,jde, kds,kde, &
           ims,ime, jms,jme, kms,kme, &
           its,ite, jts,jte, kts,kte)
      DO 100 i=its,ite
        IF(ierr(I).eq.0.)THEN
!
!--- check whether it would have buoyancy, if there where
!--- no entrainment/detrainment
!
!jm begin 20061212:  the following code replaces code that
!   was too complex and causing problem for optimization.
!   Done in consultation with G. Grell.
        jmini = jmin(i)
        keep_going = .TRUE.
        DO WHILE ( keep_going )
          keep_going = .FALSE.
!         if ( jmini - 1 .lt. kdet(i)   ) kdet(i) = jmini-1
          if ( jmini     .ge. ktop(i)-1 ) jmini = ktop(i) - 2
          ki = jmini
          hcdo(i,ki)=heso_cup(i,ki)
          DZ=Zo_cup(i,Ki+1)-Zo_cup(i,Ki)
          dh=0.
          DO k=ki-1,1,-1
            hcdo(i,k)=heso_cup(i,jmini)
            DZ=Zo_cup(i,K+1)-Zo_cup(i,K)
            dh=dh+dz*(HCDo(i,K)-heso_cup(i,k))
            IF(dh.gt.0.)THEN
              jmini=jmini-1
              IF ( jmini .gt. kbcon(i) ) THEN
                keep_going = .TRUE.
              ELSE
                ierr(i) = 9
                EXIT
              ENDIF
            ENDIF
          ENDDO
        ENDDO
        jmin(i) = jmini
        IF ( jmini .le. kbcon(i) ) THEN
          ierr(i)=4
        ENDIF
!jm end 20061212
      ENDIF
100   CONTINUE
!
! - Must have at least depth_min m between cloud convective base
!     and cloud top.
!
      do i=its,itf
      IF(ierr(I).eq.0.)THEN
      IF(-zo_cup(I,KBCON(I))+zo_cup(I,KTOP(I)).LT.depth_min)then
            ierr(i)=6
      endif
      endif
      enddo

!
!c--- normalized updraft mass flux profile
!
      call cup_up_nms(zu,z_cup,mentr_rate,cd,kbcon,ktop,ierr,k22, &
           xlamue, xlamud,            &
           ids,ide, jds,jde, kds,kde, &
           ims,ime, jms,jme, kms,kme, &
           its,ite, jts,jte, kts,kte)
      call cup_up_nms(zuo,zo_cup,mentr_rate,cd,kbcon,ktop,ierr,k22, &
           xlamue, xlamud,            &
           ids,ide, jds,jde, kds,kde, &
           ims,ime, jms,jme, kms,kme, &
           its,ite, jts,jte, kts,kte)
!
!c--- normalized downdraft mass flux profile,also work on bottom detrainment
!--- in this routine
      jprt=0
!     if(j.eq.3077)jprt=1
!
      call cup_dd_nms(zd,z_cup,cdd,mentrd_rate,jmin,ierr, &
           0,kdet,z1,                 &
              kbcon,beta,xlamdd,xlamde,                      &
           ids,ide, jds,jde, kds,kde, &
           ims,ime, jms,jme, kms,kme, &
           its,ite, jts,jte, kts,kte,jprt)
      call cup_dd_nms(zdo,zo_cup,cdd,mentrd_rate,jmin,ierr, &
           1,kdet,z1,                 &
              kbcon,beta,xlamdd,xlamde,                      &
           ids,ide, jds,jde, kds,kde, &
           ims,ime, jms,jme, kms,kme, &
           its,ite, jts,jte, kts,kte ,0)
!
!--- downdraft moist static energy
!
      call cup_dd_he(hes_cup,zd,hcd,z_cup,cdd,mentrd_rate, &
           jmin,ierr,he,dbyd,he_cup,  &
           ids,ide, jds,jde, kds,kde, &
           ims,ime, jms,jme, kms,kme, &
           its,ite, jts,jte, kts,kte)
      call cup_dd_he(heso_cup,zdo,hcdo,zo_cup,cdd,mentrd_rate, &
           jmin,ierr,heo,dbydo,he_cup,&
           ids,ide, jds,jde, kds,kde, &
           ims,ime, jms,jme, kms,kme, &
           its,ite, jts,jte, kts,kte)
!
!--- calculate moisture properties of downdraft
!
      call cup_dd_moisture(zd,hcd,hes_cup,qcd,qes_cup, &
           pwd,q_cup,z_cup,cdd,mentrd_rate,jmin,ierr,gamma_cup, &
           pwev,bu,qrcd,q,he,t_cup,2,xl, &
           ids,ide, jds,jde, kds,kde, &
           ims,ime, jms,jme, kms,kme, &
           its,ite, jts,jte, kts,kte)
      call cup_dd_moisture(zdo,hcdo,heso_cup,qcdo,qeso_cup, &
           pwdo,qo_cup,zo_cup,cdd,mentrd_rate,jmin,ierr,gammao_cup, &
           pwevo,bu,qrcdo,qo,heo,tn_cup,1,xl, &
           ids,ide, jds,jde, kds,kde, &
           ims,ime, jms,jme, kms,kme, &
           its,ite, jts,jte, kts,kte)
!
!--- calculate moisture properties of updraft
!
      call cup_up_moisture(ierr,z_cup,qc,qrc,pw,pwav, &
           kbcon,ktop,cd,dby,mentr_rate,clw_all,      &
           q,GAMMA_cup,zu,qes_cup,k22,q_cup,xl, &
           xlamue, xlamud,            &
           ids,ide, jds,jde, kds,kde, &
           ims,ime, jms,jme, kms,kme, &
           its,ite, jts,jte, kts,kte)
      do k=kts,ktf
      do i=its,itf
         cupclw(i,k)=qrc(i,k)
      enddo
      enddo
      call cup_up_moisture(ierr,zo_cup,qco,qrco,pwo,pwavo, &
           kbcon,ktop,cd,dbyo,mentr_rate,clw_all, &
           qo,GAMMAo_cup,zuo,qeso_cup,k22,qo_cup,xl,&
           xlamue, xlamud,            &
           ids,ide, jds,jde, kds,kde, &
           ims,ime, jms,jme, kms,kme, &
           its,ite, jts,jte, kts,kte)
!
!--- calculate workfunctions for updrafts
!
      call cup_up_aa0(aa0,z,zu,dby,GAMMA_CUP,t_cup, &
           kbcon,ktop,ierr,           &
           ids,ide, jds,jde, kds,kde, &
           ims,ime, jms,jme, kms,kme, &
           its,ite, jts,jte, kts,kte)
      call cup_up_aa0(aa1,zo,zuo,dbyo,GAMMAo_CUP,tn_cup, &
           kbcon,ktop,ierr,           &
           ids,ide, jds,jde, kds,kde, &
           ims,ime, jms,jme, kms,kme, &
           its,ite, jts,jte, kts,kte)
      do i=its,itf
         if(ierr(i).eq.0)then
           if(aa1(i).eq.0.)then
               ierr(i)=17
           endif
         endif
      enddo
!
!--- DETERMINE DOWNDRAFT STRENGTH IN TERMS OF WINDSHEAR
!
      call cup_dd_edt(ierr,us,vs,zo,ktop,kbcon,edt,po,pwavo, &
           pwevo,edtmax,edtmin,maxens2,edtc, &
           ids,ide, jds,jde, kds,kde, &
           ims,ime, jms,jme, kms,kme, &
           its,ite, jts,jte, kts,kte)
      do 250 iedt=1,maxens2
        do i=its,itf
      if(j.eq.jpr.and.iedt.eq.1)then
!        i=ipr
         write(0,*)'250',k22(I),kbcon(i),ktop(i),jmin(i)
         write(0,*)'250',cincr(i),pbcdif(i),pdot(i,kbcon(i))
         write(0,*)edt(i),aa0(i),aa1(i)
         do k=kts,ktf
           write(0,*)k,z(i,k),he(i,k),hes(i,k)
         enddo
         write(0,*)'end 250 loop ',iedt,edt(ipr),ierr(ipr)
         do k=1,ktop(i)+1
           write(0,*)zu(i,k),zd(i,k),pw(i,k),pwd(i,k)
         enddo
      endif

         if(ierr(i).eq.0)then
         edt(i)=edtc(i,iedt)
         edto(i)=edtc(i,iedt)
         edtx(i)=edtc(i,iedt)
         endif
        enddo
        do k=kts,ktf
        do i=its,itf
           dellat_ens(i,k,iedt)=0.
           dellaq_ens(i,k,iedt)=0.
           dellaqc_ens(i,k,iedt)=0.
           pwo_ens(i,k,iedt)=0.
        enddo
        enddo
!
      do i=its,itf
        aad(i)=0.
      enddo
!     do i=its,itf
!       if(ierr(i).eq.0)then
!        eddt(i,j)=edt(i)
!        EDTX(I)=EDT(I)
!        BU(I)=0.
!        BUO(I)=0.
!       endif
!     enddo
!
!---downdraft workfunctions
!
!     call cup_dd_aa0(edt,ierr,aa0,jmin,gamma_cup,t_cup, &
!          hcd,hes_cup,z,zd,          &
!          ids,ide, jds,jde, kds,kde, &
!          ims,ime, jms,jme, kms,kme, &
!          its,ite, jts,jte, kts,kte)
!     call cup_dd_aa0(edto,ierr,aad,jmin,gammao_cup,tn_cup, &
!          hcdo,heso_cup,zo,zdo,      &
!          ids,ide, jds,jde, kds,kde, &
!          ims,ime, jms,jme, kms,kme, &
!          its,ite, jts,jte, kts,kte)
!
!--- change per unit mass that a model cloud would modify the environment
!
!--- 1. in bottom layer
!
      call cup_dellabot(ipr,jpr,heo_cup,ierr,zo_cup,po,hcdo,edto, &
           zdo,cdd,heo,dellah,j,mentrd_rate,zo,g, &
           ids,ide, jds,jde, kds,kde, &
           ims,ime, jms,jme, kms,kme, &
           its,ite, jts,jte, kts,kte)
      call cup_dellabot(ipr,jpr,qo_cup,ierr,zo_cup,po,qrcdo,edto, &
           zdo,cdd,qo,dellaq,j,mentrd_rate,zo,g,&
           ids,ide, jds,jde, kds,kde, &
           ims,ime, jms,jme, kms,kme, &
           its,ite, jts,jte, kts,kte)
!
!--- 2. everywhere else
!
      call cup_dellas(ierr,zo_cup,po_cup,hcdo,edto,zdo,cdd,    &
           heo,dellah,j,mentrd_rate,zuo,g,                     &
           cd,hco,ktop,k22,kbcon,mentr_rate,jmin,heo_cup,kdet, &
           k22,ipr,jpr,'deep',xlamue,xlamud,                                 &
           ids,ide, jds,jde, kds,kde, &
           ims,ime, jms,jme, kms,kme, &
           its,ite, jts,jte, kts,kte)
!
!-- take out cloud liquid water for detrainment
!
!??   do k=kts,ktf
      do k=kts,ktf-1
      do i=its,itf
       scr1(i,k)=0.
       dellaqc(i,k)=0.
       if(ierr(i).eq.0)then
!       print *,'in vupnewg, after della ',ierr(i),aa0(i),i,j
         scr1(i,k)=qco(i,k)-qrco(i,k)
         if(k.eq.ktop(i)-0)dellaqc(i,k)= &
                      .01*zuo(i,ktop(i))*qrco(i,ktop(i))* &
                      9.81/(po_cup(i,k)-po_cup(i,k+1))
         if(k.lt.ktop(i).and.k.gt.kbcon(i))then
           dz=zo_cup(i,k+1)-zo_cup(i,k)
           dellaqc(i,k)=.01*9.81*cd(i,k)*dz*zuo(i,k) &
                        *.5*(qrco(i,k)+qrco(i,k+1))/ &
                        (po_cup(i,k)-po_cup(i,k+1))
         endif
       endif
      enddo
      enddo
      do i=its,itf
      dellaqc(i,ktf)=0.
      enddo
      call cup_dellas(ierr,zo_cup,po_cup,qrcdo,edto,zdo,cdd, &
           qo,dellaq,j,mentrd_rate,zuo,g, &
           cd,scr1,ktop,k22,kbcon,mentr_rate,jmin,qo_cup,kdet, &
           k22,ipr,jpr,'deep', xlamue, xlamud,             &
              ids,ide, jds,jde, kds,kde,                     &
              ims,ime, jms,jme, kms,kme,                     &
              its,ite, jts,jte, kts,kte    )
!
!--- using dellas, calculate changed environmental profiles
!
!     do 200 nens=1,maxens
      mbdt=mbdt_ens(2)
      do i=its,itf
      xaa0_ens(i,1)=0.
      xaa0_ens(i,2)=0.
      xaa0_ens(i,3)=0.
      enddo

!     mbdt=mbdt_ens(nens)
!     do i=its,itf 
!     xaa0_ens(i,nens)=0.
!     enddo
      do k=kts,ktf
      do i=its,itf
         dellat(i,k)=0.
         if(ierr(i).eq.0)then
            XHE(I,K)=DELLAH(I,K)*MBDT+HEO(I,K)
            XQ(I,K)=DELLAQ(I,K)*MBDT+QO(I,K)
            DELLAT(I,K)=(1./cp)*(DELLAH(I,K)-xl*DELLAQ(I,K))
            XT(I,K)= DELLAT(I,K)*MBDT+TN(I,K)
            IF(XQ(I,K).LE.0.)XQ(I,K)=1.E-08
!            if(i.eq.ipr.and.j.eq.jpr)then
!              print *,'dellas',k,DELLAH(I,K),DELLAQ(I,K),DELLAT(I,K)
!            endif
         ENDIF
      enddo
      enddo
      do i=its,itf
      if(ierr(i).eq.0)then
      XHE(I,ktf)=HEO(I,ktf)
      XQ(I,ktf)=QO(I,ktf)
      XT(I,ktf)=TN(I,ktf)
      IF(XQ(I,ktf).LE.0.)XQ(I,ktf)=1.E-08
      endif
      enddo
!
!--- calculate moist static energy, heights, qes
!
      call cup_env(xz,xqes,xhe,xhes,xt,xq,po,z1, &
           psur,ierr,tcrit,2,xl,cp,   &
           ids,ide, jds,jde, kds,kde, &
           ims,ime, jms,jme, kms,kme, &
           its,ite, jts,jte, kts,kte)
!
!--- environmental values on cloud levels
!
      call cup_env_clev(xt,xqes,xq,xhe,xhes,xz,po,xqes_cup,xq_cup, &
           xhe_cup,xhes_cup,xz_cup,po_cup,gamma_cup,xt_cup,psur,   &
           ierr,z1,xl,rv,cp,          &
           ids,ide, jds,jde, kds,kde, &
           ims,ime, jms,jme, kms,kme, &
           its,ite, jts,jte, kts,kte)
!
!
!**************************** static control
!
!--- moist static energy inside cloud
!
      do i=its,itf
        if(ierr(i).eq.0)then
          xhkb(i)=xhe(i,k22(i))
        endif
      enddo
      call cup_up_he(k22,xhkb,xz_cup,cd,mentr_rate,xhe_cup,xhc, &
           kbcon,ierr,xdby,xhe,xhes_cup, &
           xlamue, xlamud,            &
           ids,ide, jds,jde, kds,kde, &
           ims,ime, jms,jme, kms,kme, &
           its,ite, jts,jte, kts,kte)
!
!c--- normalized mass flux profile
!
      call cup_up_nms(xzu,xz_cup,mentr_rate,cd,kbcon,ktop,ierr,k22, &
           xlamue, xlamud,            &
           ids,ide, jds,jde, kds,kde, &
           ims,ime, jms,jme, kms,kme, &
           its,ite, jts,jte, kts,kte)
!
!--- moisture downdraft
!
      call cup_dd_nms(xzd,xz_cup,cdd,mentrd_rate,jmin,ierr, &
           1,kdet,z1,                 &
              kbcon,beta,xlamdd,xlamde,                      &
           ids,ide, jds,jde, kds,kde, &
           ims,ime, jms,jme, kms,kme, &
           its,ite, jts,jte, kts,kte,0)
      call cup_dd_he(xhes_cup,xzd,xhcd,xz_cup,cdd,mentrd_rate, &
           jmin,ierr,xhe,dbyd,xhe_cup,&
           ids,ide, jds,jde, kds,kde, &
           ims,ime, jms,jme, kms,kme, &
           its,ite, jts,jte, kts,kte)
      call cup_dd_moisture(xzd,xhcd,xhes_cup,xqcd,xqes_cup, &
           xpwd,xq_cup,xz_cup,cdd,mentrd_rate,jmin,ierr,gamma_cup, &
           xpwev,bu,xqrcd,xq,xhe,xt_cup,3,xl, &
           ids,ide, jds,jde, kds,kde, &
           ims,ime, jms,jme, kms,kme, &
           its,ite, jts,jte, kts,kte)

!
!------- MOISTURE updraft
!
      call cup_up_moisture(ierr,xz_cup,xqc,xqrc,xpw,xpwav, &
           kbcon,ktop,cd,xdby,mentr_rate,clw_all, &
           xq,GAMMA_cup,xzu,xqes_cup,k22,xq_cup,xl, &
           xlamue, xlamud,            &
           ids,ide, jds,jde, kds,kde, &
           ims,ime, jms,jme, kms,kme, &
           its,ite, jts,jte, kts,kte)
!
!--- workfunctions for updraft
!
      call cup_up_aa0(xaa0,xz,xzu,xdby,GAMMA_CUP,xt_cup, &
           kbcon,ktop,ierr,           &
           ids,ide, jds,jde, kds,kde, &
           ims,ime, jms,jme, kms,kme, &
           its,ite, jts,jte, kts,kte)
!
!--- workfunctions for downdraft
!
!
!     call cup_dd_aa0(edtx,ierr,xaa0,jmin,gamma_cup,xt_cup, &
!          xhcd,xhes_cup,xz,xzd,      &
!          ids,ide, jds,jde, kds,kde, &
!          ims,ime, jms,jme, kms,kme, &
!          its,ite, jts,jte, kts,kte)
      do 200 nens=1,maxens
      do i=its,itf 
         if(ierr(i).eq.0)then
           xaa0_ens(i,nens)=xaa0(i)
           nall=(iens-1)*maxens3*maxens*maxens2 &
                +(iedt-1)*maxens*maxens3 &
                +(nens-1)*maxens3
           do k=kts,ktf
              if(k.le.ktop(i))then
                 do nens3=1,maxens3
                 if(nens3.eq.7)then
!--- b=0
                 pr_ens(i,j,nall+nens3)=pr_ens(i,j,nall+nens3)+ &
                                    pwo(i,k) 
!                                  +edto(i)*pwdo(i,k)
!--- b=beta
                 else if(nens3.eq.8)then
                 pr_ens(i,j,nall+nens3)=pr_ens(i,j,nall+nens3)+ &
                                    pwo(i,k)
!--- b=beta/2
                 else if(nens3.eq.9)then
                 pr_ens(i,j,nall+nens3)=pr_ens(i,j,nall+nens3)+ &
                                    pwo(i,k)
!                                  +.5*edto(i)*pwdo(i,k)
                 else
                 pr_ens(i,j,nall+nens3)=pr_ens(i,j,nall+nens3)+ &
                                    pwo(i,k)+edto(i)*pwdo(i,k)
                 endif
                 enddo
              endif
           enddo
         if(pr_ens(i,j,nall+7).lt.1.e-6)then
            ierr(i)=18
            do nens3=1,maxens3
               pr_ens(i,j,nall+nens3)=0.
            enddo
         endif
         do nens3=1,maxens3
           if(pr_ens(i,j,nall+nens3).lt.1.e-4)then
            pr_ens(i,j,nall+nens3)=0.
           endif
         enddo
         endif
      enddo
 200  continue
!
!--- LARGE SCALE FORCING
!
!
!------- CHECK wether aa0 should have been zero
!
!
!     CALL cup_MAXIMI(HEO_CUP,3,KBMAX,K22x,ierr, &
!          ids,ide, jds,jde, kds,kde, &
!          ims,ime, jms,jme, kms,kme, &
!          its,ite, jts,jte, kts,kte)
      do i=its,itf
         ierr2(i)=ierr(i)
         ierr3(i)=ierr(i)
      enddo
!     call cup_kbcon(cap_max_increment,2,k22x,kbconx,heo_cup, &
!          heso_cup,ierr2,kbmax,po_cup,cap_max, &
!          ids,ide, jds,jde, kds,kde, &
!          ims,ime, jms,jme, kms,kme, &
!          its,ite, jts,jte, kts,kte)
!     call cup_kbcon(cap_max_increment,3,k22x,kbconx,heo_cup, &
!          heso_cup,ierr3,kbmax,po_cup,cap_max, &
!          ids,ide, jds,jde, kds,kde, &
!          ims,ime, jms,jme, kms,kme, &
!          its,ite, jts,jte, kts,kte)
!
!--- DETERMINE THE LEVEL OF CONVECTIVE CLOUD BASE  - KBCON
!
      call cup_forcing_ens(closure_n,xland1,aa0,aa1,xaa0_ens,mbdt_ens,dtime,   &
           ierr,ierr2,ierr3,xf_ens,j,'deeps',                 &
           maxens,iens,iedt,maxens2,maxens3,mconv,            &
           po_cup,ktop,omeg,zdo,k22,zuo,pr_ens,edto,kbcon,    &
           massflx,iact,direction,ensdim,massfln,ichoice,     &
           ids,ide, jds,jde, kds,kde, &
           ims,ime, jms,jme, kms,kme, &
           its,ite, jts,jte, kts,kte)
!
      do k=kts,ktf
      do i=its,itf
        if(ierr(i).eq.0)then
           dellat_ens(i,k,iedt)=dellat(i,k)
           dellaq_ens(i,k,iedt)=dellaq(i,k)
           dellaqc_ens(i,k,iedt)=dellaqc(i,k)
           pwo_ens(i,k,iedt)=pwo(i,k)+edt(i)*pwdo(i,k)
        else 
           dellat_ens(i,k,iedt)=0.
           dellaq_ens(i,k,iedt)=0.
           dellaqc_ens(i,k,iedt)=0.
           pwo_ens(i,k,iedt)=0.
        endif
!      if(i.eq.ipr.and.j.eq.jpr)then
!        print *,iens,iedt,dellat(i,k),dellat_ens(i,k,iedt), &
!          dellaq(i,k), dellaqc(i,k)
!      endif
      enddo
      enddo
 250  continue
!
!--- FEEDBACK
!
      do i=its,itf
        if (ierr(i).eq.0)then
          k=kbcon(i)
          dp = 100.* (p_cup(i,k) - p_cup(i,k+1))
          xmbmax(i) = dp / (g * dtime)
!         if (j.eq.3077) print *,'xmbmax=',xmbmax(i),dp,g,dtime
        endif
      enddo

      call cup_output_ens(xf_ens,ierr,dellat_ens,dellaq_ens, &
           dellaqc_ens,outt,outq,outqc,pre,pwo_ens,xmb,ktop, &
           j,'deep',maxens2,maxens,iens,ierr2,ierr3,         &
           pr_ens,maxens3,ensdim,massfln,xmbmax,             &
           APR_GR,APR_W,APR_MC,APR_ST,APR_AS,                &
           APR_CAPMA,APR_CAPME,APR_CAPMI,closure_n,xland1,   &
           ids,ide, jds,jde, kds,kde, &
           ims,ime, jms,jme, kms,kme, &
           its,ite, jts,jte, kts,kte)
      do i=its,itf
           PRE(I)=MAX(PRE(I),0.)
           if(i.eq.ipr.and.j.eq.jpr)then
             write(0,*)'j,pre(i),aa0(i),aa1(i)'
             write(0,*)j,pre(i),aa0(i),aa1(i)
           endif

      enddo
!
!---------------------------done------------------------------
!

   END SUBROUTINE CUP_enss


   SUBROUTINE cup_dd_aa0(edt,ierr,aa0,jmin,gamma_cup,t_cup, &
              hcd,hes_cup,z,zd,                             &
              ids,ide, jds,jde, kds,kde,                    &
              ims,ime, jms,jme, kms,kme,                    &
              its,ite, jts,jte, kts,kte                    )

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
        ids,ide, jds,jde, kds,kde,                                     &
        ims,ime, jms,jme, kms,kme,                                     &
        its,ite, jts,jte, kts,kte
  ! aa0 cloud work function for downdraft
  ! gamma_cup = gamma on model cloud levels
  ! t_cup = temperature (Kelvin) on model cloud levels
  ! hes_cup = saturation moist static energy on model cloud levels
  ! hcd = moist static energy in downdraft
  ! edt = epsilon
  ! zd normalized downdraft mass flux
  ! z = heights of model levels 
  ! ierr error value, maybe modified in this routine
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        z,zd,gamma_cup,t_cup,hes_cup,hcd
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        edt
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        jmin
!
! input and output
!


     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
     real,    dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        aa0
!
!  local variables in this routine
!

     integer                              ::                           &
        i,k,kk
     real                                 ::                           &
        dz
!
     integer :: itf, ktf
!
     itf=ite
     ktf=kte
!
!??    DO k=kts,kte-1
       DO k=kts,ktf-1
       do i=its,itf
         IF(ierr(I).eq.0.and.k.lt.jmin(i))then
         KK=JMIN(I)-K
!
!--- ORIGINAL
!
         DZ=(Z(I,KK)-Z(I,KK+1))
         AA0(I)=AA0(I)+zd(i,kk)*EDT(I)*DZ*(9.81/(1004.*T_cup(I,KK))) &
            *((hcd(i,kk)-hes_cup(i,kk))/(1.+GAMMA_cup(i,kk)))
         endif
      enddo
      enddo

   END SUBROUTINE CUP_dd_aa0


   SUBROUTINE cup_dd_edt(ierr,us,vs,z,ktop,kbcon,edt,p,pwav, &
              pwev,edtmax,edtmin,maxens2,edtc,               &
              ids,ide, jds,jde, kds,kde,                     &
              ims,ime, jms,jme, kms,kme,                     &
              its,ite, jts,jte, kts,kte                     )

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        ids,ide, jds,jde, kds,kde,           &
        ims,ime, jms,jme, kms,kme,           &
        its,ite, jts,jte, kts,kte
     integer, intent (in   )              ::                           &
        maxens2
  !
  ! ierr error value, maybe modified in this routine
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        us,vs,z,p
     real,    dimension (its:ite,1:maxens2)                            &
        ,intent (out  )                   ::                           &
        edtc
     real,    dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        edt
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        pwav,pwev
     real                                                              &
        ,intent (in   )                   ::                           &
        edtmax,edtmin
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        ktop,kbcon
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
!
!  local variables in this routine
!

     integer i,k,kk
     real    einc,pef,pefb,prezk,zkbc
     real,    dimension (its:ite)         ::                           &
      vshear,sdp,vws

     integer :: itf, ktf

     itf=ite
     ktf=kte
!
!--- DETERMINE DOWNDRAFT STRENGTH IN TERMS OF WINDSHEAR
!
! */ calculate an average wind shear over the depth of the cloud
!
       do i=its,itf
        edt(i)=0.
        vws(i)=0.
        sdp(i)=0.
        vshear(i)=0.
       enddo
       do kk = kts,ktf-1
         do 62 i=its,itf
          IF(ierr(i).ne.0)GO TO 62
          if (kk .le. min0(ktop(i),ktf-1) .and. kk .ge. kbcon(i)) then
             vws(i) = vws(i)+ &
              (abs((us(i,kk+1)-us(i,kk))/(z(i,kk+1)-z(i,kk))) &
          +   abs((vs(i,kk+1)-vs(i,kk))/(z(i,kk+1)-z(i,kk)))) * &
              (p(i,kk) - p(i,kk+1))
            sdp(i) = sdp(i) + p(i,kk) - p(i,kk+1)
          endif
          if (kk .eq. ktf-1)vshear(i) = 1.e3 * vws(i) / sdp(i)
   62   continue
       end do
      do i=its,itf
         IF(ierr(i).eq.0)then
            pef=(1.591-.639*VSHEAR(I)+.0953*(VSHEAR(I)**2) &
               -.00496*(VSHEAR(I)**3))
            if(pef.gt.edtmax)pef=edtmax
            if(pef.lt.edtmin)pef=edtmin
!
!--- cloud base precip efficiency
!
            zkbc=z(i,kbcon(i))*3.281e-3
            prezk=.02
            if(zkbc.gt.3.)then
               prezk=.96729352+zkbc*(-.70034167+zkbc*(.162179896+zkbc &
               *(- 1.2569798E-2+zkbc*(4.2772E-4-zkbc*5.44E-6))))
            endif
            if(zkbc.gt.25)then
               prezk=2.4
            endif
            pefb=1./(1.+prezk)
            if(pefb.gt.edtmax)pefb=edtmax
            if(pefb.lt.edtmin)pefb=edtmin
            EDT(I)=1.-.5*(pefb+pef)
!--- edt here is 1-precipeff!
!           einc=(1.-edt(i))/float(maxens2)
!           einc=edt(i)/float(maxens2+1)
!--- 20 percent
            einc=.2*edt(i)
            do k=1,maxens2
                edtc(i,k)=edt(i)+float(k-2)*einc
            enddo
         endif
      enddo
      do i=its,itf
         IF(ierr(i).eq.0)then
            do k=1,maxens2
               EDTC(I,K)=-EDTC(I,K)*PWAV(I)/PWEV(I)
               IF(EDTC(I,K).GT.edtmax)EDTC(I,K)=edtmax
               IF(EDTC(I,K).LT.edtmin)EDTC(I,K)=edtmin
            enddo
         endif
      enddo

   END SUBROUTINE cup_dd_edt


   SUBROUTINE cup_dd_he(hes_cup,zd,hcd,z_cup,cdd,entr,       &
              jmin,ierr,he,dby,he_cup,                       &
              ids,ide, jds,jde, kds,kde,                     &
              ims,ime, jms,jme, kms,kme,                     &
              its,ite, jts,jte, kts,kte                     )

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
                                  ids,ide, jds,jde, kds,kde,           &
                                  ims,ime, jms,jme, kms,kme,           &
                                  its,ite, jts,jte, kts,kte
  ! hcd = downdraft moist static energy
  ! he = moist static energy on model levels
  ! he_cup = moist static energy on model cloud levels
  ! hes_cup = saturation moist static energy on model cloud levels
  ! dby = buoancy term
  ! cdd= detrainment function 
  ! z_cup = heights of model cloud levels 
  ! entr = entrainment rate
  ! zd   = downdraft normalized mass flux
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        he,he_cup,hes_cup,z_cup,cdd,zd
  ! entr= entrainment rate 
     real                                                              &
        ,intent (in   )                   ::                           &
        entr
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        jmin
!
! input and output
!

   ! ierr error value, maybe modified in this routine

     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr

     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
        hcd,dby
!
!  local variables in this routine
!

     integer                              ::                           &
        i,k,ki
     real                                 ::                           &
        dz

     integer :: itf, ktf

     itf=ite
     ktf=kte

      do k=kts+1,ktf
      do i=its,itf
      dby(i,k)=0.
      IF(ierr(I).eq.0)then
         hcd(i,k)=hes_cup(i,k)
      endif
      enddo
      enddo
!
      do 100 i=its,itf
      IF(ierr(I).eq.0)then
      k=jmin(i)
      hcd(i,k)=hes_cup(i,k)
      dby(i,k)=hcd(i,jmin(i))-hes_cup(i,k)
!
      do ki=jmin(i)-1,1,-1
         DZ=Z_cup(i,Ki+1)-Z_cup(i,Ki)
         HCD(i,Ki)=(HCD(i,Ki+1)*(1.-.5*CDD(i,Ki)*DZ) &
                  +entr*DZ*HE(i,Ki) &
                  )/(1.+entr*DZ-.5*CDD(i,Ki)*DZ)
         dby(i,ki)=HCD(i,Ki)-hes_cup(i,ki)
      enddo
!
      endif
!--- end loop over i
100    continue


   END SUBROUTINE cup_dd_he


   SUBROUTINE cup_dd_moisture(zd,hcd,hes_cup,qcd,qes_cup,    &
              pwd,q_cup,z_cup,cdd,entr,jmin,ierr,            &
              gamma_cup,pwev,bu,qrcd,                        &
              q,he,t_cup,iloop,xl,                           &
              ids,ide, jds,jde, kds,kde,                     &
              ims,ime, jms,jme, kms,kme,                     &
              its,ite, jts,jte, kts,kte                     )

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
                                  ids,ide, jds,jde, kds,kde,           &
                                  ims,ime, jms,jme, kms,kme,           &
                                  its,ite, jts,jte, kts,kte
  ! cdd= detrainment function 
  ! q = environmental q on model levels
  ! q_cup = environmental q on model cloud levels
  ! qes_cup = saturation q on model cloud levels
  ! hes_cup = saturation h on model cloud levels
  ! hcd = h in model cloud
  ! bu = buoancy term
  ! zd = normalized downdraft mass flux
  ! gamma_cup = gamma on model cloud levels
  ! mentr_rate = entrainment rate
  ! qcd = cloud q (including liquid water) after entrainment
  ! qrch = saturation q in cloud
  ! pwd = evaporate at that level
  ! pwev = total normalized integrated evaoprate (I2)
  ! entr= entrainment rate 
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        zd,t_cup,hes_cup,hcd,qes_cup,q_cup,z_cup,cdd,gamma_cup,q,he 
     real                                                              &
        ,intent (in   )                   ::                           &
        entr,xl
     integer                                                           &
        ,intent (in   )                   ::                           &
        iloop
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        jmin
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
        qcd,qrcd,pwd
     real,    dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        pwev,bu
!
!  local variables in this routine
!

     integer                              ::                           &
        i,k,ki
     real                                 ::                           &
        dh,dz,dqeva

     integer :: itf, ktf

     itf=ite
     ktf=kte

      do i=its,itf
         bu(i)=0.
         pwev(i)=0.
      enddo
      do k=kts,ktf
      do i=its,itf
         qcd(i,k)=0.
         qrcd(i,k)=0.
         pwd(i,k)=0.
      enddo
      enddo
!
!
!
      do 100 i=its,itf
      IF(ierr(I).eq.0)then
      k=jmin(i)
      DZ=Z_cup(i,K+1)-Z_cup(i,K)
      qcd(i,k)=q_cup(i,k)
!     qcd(i,k)=.5*(qes_cup(i,k)+q_cup(i,k))
      qrcd(i,k)=qes_cup(i,k)
      pwd(i,jmin(i))=min(0.,qcd(i,k)-qrcd(i,k))
      pwev(i)=pwev(i)+pwd(i,jmin(i))
      qcd(i,k)=qes_cup(i,k)
!
      DH=HCD(I,k)-HES_cup(I,K)
      bu(i)=dz*dh
      do ki=jmin(i)-1,1,-1
         DZ=Z_cup(i,Ki+1)-Z_cup(i,Ki)
         QCD(i,Ki)=(qCD(i,Ki+1)*(1.-.5*CDD(i,Ki)*DZ) &
                  +entr*DZ*q(i,Ki) &
                  )/(1.+entr*DZ-.5*CDD(i,Ki)*DZ)
!
!--- to be negatively buoyant, hcd should be smaller than hes!
!
         DH=HCD(I,ki)-HES_cup(I,Ki)
         bu(i)=bu(i)+dz*dh
         QRCD(I,Ki)=qes_cup(i,ki)+(1./XL)*(GAMMA_cup(i,ki) &
                  /(1.+GAMMA_cup(i,ki)))*DH
         dqeva=qcd(i,ki)-qrcd(i,ki)
         if(dqeva.gt.0.)dqeva=0.
         pwd(i,ki)=zd(i,ki)*dqeva
         qcd(i,ki)=qrcd(i,ki)
         pwev(i)=pwev(i)+pwd(i,ki)
!        if(iloop.eq.1.and.i.eq.102.and.j.eq.62)then
!         print *,'in cup_dd_moi ', hcd(i,ki),HES_cup(I,Ki),dh,dqeva
!        endif
      enddo
!
!--- end loop over i
       if(pwev(I).eq.0.and.iloop.eq.1)then
!        print *,'problem with buoy in cup_dd_moisture',i
         ierr(i)=7
       endif
       if(BU(I).GE.0.and.iloop.eq.1)then
!        print *,'problem with buoy in cup_dd_moisture',i
         ierr(i)=7
       endif
      endif
100    continue

   END SUBROUTINE cup_dd_moisture


   SUBROUTINE cup_dd_nms(zd,z_cup,cdd,entr,jmin,ierr,        &
              itest,kdet,z1,                                 &
              kbcon,beta,xlamdd,xlamde,                      &
              ids,ide, jds,jde, kds,kde,                     &
              ims,ime, jms,jme, kms,kme,                     &
              its,ite, jts,jte, kts,kte,jpr                  )

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
                                  ids,ide, jds,jde, kds,kde,           &
                                  ims,ime, jms,jme, kms,kme,           &
                                  its,ite, jts,jte, kts,kte,jpr
  ! z_cup = height of cloud model level
  ! z1 = terrain elevation
  ! entr = downdraft entrainment rate
  ! jmin = downdraft originating level
  ! kdet = level above ground where downdraft start detraining
  ! itest = flag to whether to calculate cdd
  
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        z_cup
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        z1 
     real                                                              &
        ,intent (in   )                   ::                           &
        entr,xlamdd,xlamde,beta 
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        jmin,kdet,kbcon
     integer                                                           &
        ,intent (in   )                   ::                           &
        itest
!
! input and output
!

   ! ierr error value, maybe modified in this routine

     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
                                                                 ierr
   ! zd is the normalized downdraft mass flux
   ! cdd is the downdraft detrainmen function

     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
                                                             zd
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (inout)                   ::                           &
                                                             cdd
!
!  local variables in this routine
!

     integer                              ::                           &
                                                  i,k,ki
     real                                 ::                           &
                                            a,perc,dz,ptem,tem
     real, dimension (its:ite)         ::   xlamd

     integer :: itf, ktf

     itf=ite
     ktf=kte
!
!--- perc is the percentage of mass left when hitting the ground
!
      perc=.03

      do k=kts,ktf
      do i=its,itf
         zd(i,k)=0.
      enddo
      enddo
      a=1.-perc
! only define this entrain/detrain nonsense once
      if (itest.eq.0) then

!       if(jpr.eq.1)write(6,*)'itest = 0 here '
        do k=kts,ktf
        do i=its,itf
          cdd(i,k)=0.
        enddo
        enddo

        do i=its,itf
          IF(ierr(I).eq.0)then
            dz = z_cup(i,kbcon(i))/float(kbcon(i))
            tem = 1./float(kbcon(i))
            xlamd(i) = (1.-beta**tem)/dz
            if(jpr.eq.1)write(6,*)'b',beta**tem,tem,beta,xlamd(i)

            do ki=jmin(i)-1,1,-1
              if (ki.ge.kbcon(i)) then
                cdd(i,ki) = xlamdd
              else
                cdd(i,ki) = xlamdd + xlamd(i)
              endif
!           if(jpr.eq.1)write(6,*)ki,cdd(i,ki),xlamdd,xlamd(i)
            end do

          ENDIF
        enddo
      endif ! itest=0
!
!
!
      do 100 i=its,itf
      IF(ierr(I).eq.0)then
      zd(i,jmin(i))=1.
!
!--- integrate downward, specify detrainment(cdd)!
!
      do ki=jmin(i)-1,1,-1
         DZ=Z_cup(i,Ki+1)-Z_cup(i,Ki)
! Thorpex 20110506  (SASRQ3)
            ptem = cdd(i,ki) - xlamde
            zd(i,ki) = zd(i,ki+1)*(1. - ptem * dz)
!           if(jpr.eq.1)write(6,*)ki,cdd(i,ki),zd(i,ki),xlamde

!        if(ki.le.kdet(i).and.itest.eq.0)then
!          cdd(i,ki)=entr+(1.- (a*(z_cup(i,ki)-z1(i)) &
!                    +perc*(z_cup(i,kdet(i))-z1(i)) ) &
!                        /(a*(z_cup(i,ki+1)-z1(i)) &
!                     +perc*(z_cup(i,kdet(i))-z1(i))))/dz
!        endif
!        zd(i,ki)=zd(i,ki+1)*(1.+(entr-cdd(i,ki))*dz)
      enddo
!
      endif
!--- end loop over i
100    continue

   END SUBROUTINE cup_dd_nms


   SUBROUTINE cup_dellabot(ipr,jpr,he_cup,ierr,z_cup,p_cup,  &
              hcd,edt,zd,cdd,he,della,j,mentrd_rate,z,g,     &
              ids,ide, jds,jde, kds,kde,                     &
              ims,ime, jms,jme, kms,kme,                     &
              its,ite, jts,jte, kts,kte                     )

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        ids,ide, jds,jde, kds,kde,           &
        ims,ime, jms,jme, kms,kme,           &
        its,ite, jts,jte, kts,kte
     integer, intent (in   )              ::                           &
        j,ipr,jpr
  !
  ! ierr error value, maybe modified in this routine
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
        della
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in  )                   ::                           &
        z_cup,p_cup,hcd,zd,cdd,he,z,he_cup
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        edt
     real                                                              &
        ,intent (in   )                   ::                           &
        g,mentrd_rate
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
!
!  local variables in this routine
!

      integer i
      real detdo,detdo1,detdo2,entdo,dp,dz,subin,                      &
      totmas
!
     integer :: itf, ktf

     itf=ite
     ktf=kte
!
!
!      if(j.eq.jpr)print *,'in cup dellabot '
      do 100 i=its,itf
      della(i,1)=0.
      if(ierr(i).ne.0)go to 100
      dz=z_cup(i,2)-z_cup(i,1)
      DP=100.*(p_cup(i,1)-P_cup(i,2))
      detdo1=edt(i)*zd(i,2)*CDD(i,1)*DZ
      detdo2=edt(i)*zd(i,1)
      entdo=edt(i)*zd(i,2)*mentrd_rate*dz
      subin=-EDT(I)*zd(i,2)
      detdo=detdo1+detdo2-entdo+subin
      if(detdo.gt.1.e-6)then
       write(6,*)'totmas = ',detdo
       write(6,*)detdo1,detdo2,entdo,subin
      endif
!     DELLA(I,1)=(detdo1*.5*(HCD(i,1)+HCD(i,2)) &
!                +detdo2*hcd(i,1) &
!                +subin*he_cup(i,2) &
!                -entdo*he(i,1))*g/dp
      DELLA(I,1)=detdo2*(hcd(i,2) -he_cup(i,2))*g/dp
 100  CONTINUE

   END SUBROUTINE cup_dellabot


   SUBROUTINE cup_dellas(ierr,z_cup,p_cup,hcd,edt,zd,cdd,              &
              he,della,j,mentrd_rate,zu,g,                             &
              cd,hc,ktop,k22,kbcon,mentr_rate,jmin,he_cup,kdet,kpbl,   &
              ipr,jpr,name, xlamue, xlamud,                                          &
              ids,ide, jds,jde, kds,kde,                               &
              ims,ime, jms,jme, kms,kme,                               &
              its,ite, jts,jte, kts,kte                               )

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        ids,ide, jds,jde, kds,kde,           &
        ims,ime, jms,jme, kms,kme,           &
        its,ite, jts,jte, kts,kte
     integer, intent (in   )              ::                           &
        j,ipr,jpr
  !
  ! ierr error value, maybe modified in this routine
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
        della
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in  )                   ::                           &
        z_cup,p_cup,hcd,zd,cdd,he,hc,cd,zu,he_cup,xlamue
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        edt,xlamud
     real                                                              &
        ,intent (in   )                   ::                           &
        g,mentrd_rate,mentr_rate
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        kbcon,ktop,k22,jmin,kdet,kpbl
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
      character *(*), intent (in)        ::                           &
       name
!
!  local variables in this routine
!

      integer i,k
      real detdo1,detdo2,entdo,dp,dz,subin,detdo,entup,                &
      detup,subdown,entdoj,entupk,detupk,totmas
!
     integer :: itf, ktf

     itf=ite
     ktf=kte
!
!
      i=ipr
!      if(j.eq.jpr)then
!        print *,'in dellas kpbl(i),k22(i),kbcon(i),ktop(i),jmin(i)'
!        print *,kpbl(i),k22(i),kbcon(i),ktop(i),jmin(i)
!      endif
       DO K=kts+1,ktf
       do i=its,itf
          della(i,k)=0.
       enddo
       enddo
!
! if k22=2, there is also a term due to updrafts at level1 !
! this will overwrite dellabot in this case
!
       DO i=its,ite
       IF(ierr(i).eq.0 .and. k22(i).eq.2)then
            detdo2=edt(i)*zd(i,1)
            entupk=zu(i,k22(i))
            dp=100.*(p_cup(i,1)-p_cup(i,2))
            DELLA(I,1)=(detdo2*(hcd(i,2) -he_cup(i,2))        &
                       -entupk*(he_cup(i,k22(i))-he_cup(i,2)) &
                       )*g/dp
       endif
       enddo
       DO 100 k=kts+1,ktf-1
       DO 100 i=its,ite
         IF(ierr(i).ne.0)GO TO 100
         IF(K.Gt.KTOP(I))GO TO 100
!
!--- SPECIFY DETRAINMENT OF DOWNDRAFT, HAS TO BE CONSISTENT
!--- WITH ZD CALCULATIONS IN SOUNDD.
!
! updraft starts at k22 (>1). thi will effect della(k22-1). At that level (k22-1) 
! resulting subin is active ... updraft entr/detr start going from kbcon  
!  to kbcon+1, which would effect kbcon for dellas
         DZ=Z_cup(I,K+1)-Z_cup(I,K)
         detdo=edt(i)*CDD(i,K)*DZ*ZD(i,k+1)
         entdo=edt(i)*mentrd_rate*dz*zd(i,k+1)
         subin=zu(i,k+1)-zd(i,k+1)*edt(i)
         subdown=(zu(i,k)-zd(i,k)*edt(i))
         entup=0.
         detup=0.
         entdoj=0.
         entupk=0.
         detupk=0.
         if(k.ge.kbcon(i).and.k.lt.ktop(i))then
            entup=xlamue(i,k+1)*dz*zu(i,k+0)
            detup=xlamud(i)*DZ*ZU(i,k+0)
         endif
!
         if(k.eq.jmin(i))then
         entdoj=edt(i)*zd(i,k)
         endif
! k is at least kts+1, if k22 happens to be 2, this is handled above loop 100
         if(k.eq.k22(i)-1)then
            entupk=zu(i,k22(i))
            dp=100.*(p_cup(i,k)-p_cup(i,k+1))
            della(i,k)=(subin*he_cup(i,k+1) &
                    -subdown*he_cup(i,k) &
                    +detdo*.5*(HCD(i,K+1)+HCD(i,K)) &
                    -entdo*he(i,k) &
                    -entupk*he_cup(i,k22(i)) &
                    -entdoj*he_cup(i,jmin(i)) &
                     )*g/dp
         endif

!        if(k.gt.kdet(i))then
!           detdo=0.
!        endif

         if(k.eq.ktop(i)-0)then
         detupk=zu(i,ktop(i))
         subin=0.
         endif
         if(k.le.k22(i))then
            detup=0.
         endif
!C
!C--- CHANGED DUE TO SUBSIDENCE AND ENTRAINMENT
!C
         totmas=subin-subdown+detup-entup-entdo+ &
                 detdo-entupk-entdoj+detupk
!         if(j.eq.jpr.and.i.eq.ipr)print *,'k,totmas,sui,sud = ',k,
!     1   totmas,subin,subdown
!         if(j.eq.jpr.and.i.eq.ipr)print *,'updr stuff = ',detup,
!     1      entup,entupk,detupk
!         if(j.eq.jpr.and.i.eq.ipr)print *,'dddr stuff = ',entdo,
!     1      detdo,entdoj
         if(abs(totmas).gt.1.e-6)then
!        if(j.eq.jpr)then
            print *,'*********************',i,j,k,totmas,name
            print *,'k22,kbcon,ktop,jmin= ',k22(i),kbcon(i),ktop(i),jmin(i)
          print *,'updr stuff = ',subin,subdown,detup,entup,entupk,detupk
          print *,'dddr stuff = ',entdo,detdo,entdoj
!        call wrf_error_fatal ( 'totmas .gt.1.e-6' )
         endif
!        dp=100.*(p_cup(i,k-1)-p_cup(i,k))
         dp=100.*(p_cup(i,k)-p_cup(i,k+1))
         della(i,k)=(subin*he_cup(i,k+1) &
                    -subdown*he_cup(i,k) &
                    +detup*.5*(HC(i,K+1)+HC(i,K)) &
                    +detdo*.5*(HCD(i,K+1)+HCD(i,K)) &
                    -entup*he(i,k) &
                    -entdo*he(i,k) &
                    -entupk*he_cup(i,k22(i)+1) &
                    -entdoj*he_cup(i,jmin(i)) &
                    +detupk*hc(i,ktop(i)) &
                     )*g/dp
!      if(i.eq.ipr.and.j.eq.jpr)then
!        print *,k,della(i,k),subin*he_cup(i,k+1),subdown*he_cup(i,k),
!     1            detdo*.5*(HCD(i,K+1)+HCD(i,K))
!        print *,k,detup*.5*(HC(i,K+1)+HC(i,K)),detupk*hc(i,ktop(i)),
!     1         entup*he(i,k),entdo*he(i,k)
!        print *,k,he_cup(i,k+1),he_cup(i,k),entupk*he_cup(i,k)
!      endif

 100  CONTINUE

   END SUBROUTINE cup_dellas


   SUBROUTINE cup_direction2(i,j,dir,id,massflx,             &
              iresult,imass,massfld,                         &
              ids,ide, jds,jde, kds,kde,                     &
              ims,ime, jms,jme, kms,kme,                     &
              its,ite, jts,jte, kts,kte                     )

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        ids,ide, jds,jde, kds,kde,           &
        ims,ime, jms,jme, kms,kme,           &
        its,ite, jts,jte, kts,kte
     integer, intent (in   )              ::                           &
        i,j,imass
     integer, intent (out  )              ::                           &
        iresult
  !
  ! ierr error value, maybe modified in this routine
  !
     integer,    dimension (ims:ime,jms:jme)                           &
        ,intent (in   )                   ::                           &
        id
     real,    dimension (ims:ime,jms:jme)                              &
        ,intent (in   )                   ::                           &
        massflx
     real,    dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        dir
     real                                                              &
        ,intent (out  )                   ::                           &
        massfld
!
!  local variables in this routine
!

       integer k,ia,ja,ib,jb
       real diff
!
!
!
       if(imass.eq.1)then
           massfld=massflx(i,j)
       endif
       iresult=0
!      return
       diff=22.5
       if(dir(i).lt.22.5)dir(i)=360.+dir(i)
       if(id(i,j).eq.1)iresult=1
!      ja=max(2,j-1)
!      ia=max(2,i-1)
!      jb=min(mjx-1,j+1)
!      ib=min(mix-1,i+1)
       ja=j-1
       ia=i-1
       jb=j+1
       ib=i+1
        if(dir(i).gt.90.-diff.and.dir(i).le.90.+diff)then
!--- steering flow from the east
          if(id(ib,j).eq.1)then
            iresult=1
            if(imass.eq.1)then
               massfld=max(massflx(ib,j),massflx(i,j))
            endif
            return
          endif
        else if(dir(i).gt.135.-diff.and.dir(i).le.135.+diff)then
!--- steering flow from the south-east
          if(id(ib,ja).eq.1)then
            iresult=1
            if(imass.eq.1)then
               massfld=max(massflx(ib,ja),massflx(i,j))
            endif
            return
          endif
!--- steering flow from the south
        else if(dir(i).gt.180.-diff.and.dir(i).le.180.+diff)then
          if(id(i,ja).eq.1)then
            iresult=1
            if(imass.eq.1)then
               massfld=max(massflx(i,ja),massflx(i,j))
            endif
            return
          endif
!--- steering flow from the south west
        else if(dir(i).gt.225.-diff.and.dir(i).le.225.+diff)then
          if(id(ia,ja).eq.1)then
            iresult=1
            if(imass.eq.1)then
               massfld=max(massflx(ia,ja),massflx(i,j))
            endif
            return
          endif
!--- steering flow from the west
        else if(dir(i).gt.270.-diff.and.dir(i).le.270.+diff)then
          if(id(ia,j).eq.1)then
            iresult=1
            if(imass.eq.1)then
               massfld=max(massflx(ia,j),massflx(i,j))
            endif
            return
          endif
!--- steering flow from the north-west
        else if(dir(i).gt.305.-diff.and.dir(i).le.305.+diff)then
          if(id(ia,jb).eq.1)then
            iresult=1
            if(imass.eq.1)then
               massfld=max(massflx(ia,jb),massflx(i,j))
            endif
            return
          endif
!--- steering flow from the north
        else if(dir(i).gt.360.-diff.and.dir(i).le.360.+diff)then
          if(id(i,jb).eq.1)then
            iresult=1
            if(imass.eq.1)then
               massfld=max(massflx(i,jb),massflx(i,j))
            endif
            return
          endif
!--- steering flow from the north-east
        else if(dir(i).gt.45.-diff.and.dir(i).le.45.+diff)then
          if(id(ib,jb).eq.1)then
            iresult=1
            if(imass.eq.1)then
               massfld=max(massflx(ib,jb),massflx(i,j))
            endif
            return
          endif
        endif

   END SUBROUTINE cup_direction2


   SUBROUTINE cup_env(z,qes,he,hes,t,q,p,z1,                 &
              psur,ierr,tcrit,itest,xl,cp,                   &
              ids,ide, jds,jde, kds,kde,                     &
              ims,ime, jms,jme, kms,kme,                     &
              its,ite, jts,jte, kts,kte                     )

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        ids,ide, jds,jde, kds,kde,           &
        ims,ime, jms,jme, kms,kme,           &
        its,ite, jts,jte, kts,kte
  !
  ! ierr error value, maybe modified in this routine
  ! q           = environmental mixing ratio
  ! qes         = environmental saturation mixing ratio
  ! t           = environmental temp
  ! tv          = environmental virtual temp
  ! p           = environmental pressure
  ! z           = environmental heights
  ! he          = environmental moist static energy
  ! hes         = environmental saturation moist static energy
  ! psur        = surface pressure
  ! z1          = terrain elevation
  ! 
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        p,t
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
        he,hes,qes
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (inout)                   ::                           &
        z,q
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        psur,z1
     real                                                              &
        ,intent (in   )                   ::                           &
        xl,cp
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
     integer                                                           &
        ,intent (in   )                   ::                           &
        itest
!
!  local variables in this routine
!

     integer                              ::                           &
       i,k,iph
      real, dimension (1:2) :: AE,BE,HT
      real, dimension (its:ite,kts:kte) :: tv
      real :: tcrit,e,tvbar

     integer :: itf, ktf

     itf=ite
     ktf=kte

      HT(1)=XL/CP
      HT(2)=2.834E6/CP
      BE(1)=.622*HT(1)/.286
      AE(1)=BE(1)/273.+ALOG(610.71)
      BE(2)=.622*HT(2)/.286
      AE(2)=BE(2)/273.+ALOG(610.71)
!      print *, 'TCRIT = ', tcrit,its,ite
      DO k=kts,ktf
      do i=its,itf
        if(ierr(i).eq.0)then
!Csgb - IPH is for phase, dependent on TCRIT (water or ice)
        IPH=1
        IF(T(I,K).LE.TCRIT)IPH=2
!       print *, 'AE(IPH),BE(IPH) = ',AE(IPH),BE(IPH),AE(IPH)-BE(IPH),T(i,k),i,k
        E=EXP(AE(IPH)-BE(IPH)/T(I,K))
!       print *, 'P, E = ', P(I,K), E
        QES(I,K)=.622*E/(100.*P(I,K)-E)
        IF(QES(I,K).LE.1.E-08)QES(I,K)=1.E-08
        IF(Q(I,K).GT.QES(I,K))Q(I,K)=QES(I,K)
        TV(I,K)=T(I,K)+.608*Q(I,K)*T(I,K)
        endif
      enddo
      enddo
!
!--- z's are calculated with changed h's and q's and t's
!--- if itest=2
!
      if(itest.ne.2)then
         do i=its,itf
           if(ierr(i).eq.0)then
             Z(I,1)=max(0.,Z1(I))-(ALOG(P(I,1))- &
                 ALOG(PSUR(I)))*287.*TV(I,1)/9.81
           endif
         enddo

! --- calculate heights
         DO K=kts+1,ktf
         do i=its,itf
           if(ierr(i).eq.0)then
              TVBAR=.5*TV(I,K)+.5*TV(I,K-1)
              Z(I,K)=Z(I,K-1)-(ALOG(P(I,K))- &
               ALOG(P(I,K-1)))*287.*TVBAR/9.81
           endif
         enddo
         enddo
      else
         do k=kts,ktf
         do i=its,itf
           if(ierr(i).eq.0)then
             z(i,k)=(he(i,k)-1004.*t(i,k)-2.5e6*q(i,k))/9.81
             z(i,k)=max(1.e-3,z(i,k))
           endif
         enddo
         enddo
      endif
!
!--- calculate moist static energy - HE
!    saturated moist static energy - HES
!
       DO k=kts,ktf
       do i=its,itf
         if(ierr(i).eq.0)then
         if(itest.eq.0)HE(I,K)=9.81*Z(I,K)+1004.*T(I,K)+2.5E06*Q(I,K)
         HES(I,K)=9.81*Z(I,K)+1004.*T(I,K)+2.5E06*QES(I,K)
         IF(HE(I,K).GE.HES(I,K))HE(I,K)=HES(I,K)
!         if(i.eq.2)then
!           print *,k,z(i,k),t(i,k),p(i,k),he(i,k),hes(i,k)
!         endif
         endif
      enddo
      enddo

   END SUBROUTINE cup_env


   SUBROUTINE cup_env_clev(t,qes,q,he,hes,z,p,qes_cup,q_cup,   &
              he_cup,hes_cup,z_cup,p_cup,gamma_cup,t_cup,psur, &
              ierr,z1,xl,rv,cp,                                &
              ids,ide, jds,jde, kds,kde,                       &
              ims,ime, jms,jme, kms,kme,                       &
              its,ite, jts,jte, kts,kte                       )

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        ids,ide, jds,jde, kds,kde,           &
        ims,ime, jms,jme, kms,kme,           &
        its,ite, jts,jte, kts,kte
  !
  ! ierr error value, maybe modified in this routine
  ! q           = environmental mixing ratio
  ! q_cup       = environmental mixing ratio on cloud levels
  ! qes         = environmental saturation mixing ratio
  ! qes_cup     = environmental saturation mixing ratio on cloud levels
  ! t           = environmental temp
  ! t_cup       = environmental temp on cloud levels
  ! p           = environmental pressure
  ! p_cup       = environmental pressure on cloud levels
  ! z           = environmental heights
  ! z_cup       = environmental heights on cloud levels
  ! he          = environmental moist static energy
  ! he_cup      = environmental moist static energy on cloud levels
  ! hes         = environmental saturation moist static energy
  ! hes_cup     = environmental saturation moist static energy on cloud levels
  ! gamma_cup   = gamma on cloud levels
  ! psur        = surface pressure
  ! z1          = terrain elevation
  ! 
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        qes,q,he,hes,z,p,t
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
        qes_cup,q_cup,he_cup,hes_cup,z_cup,p_cup,gamma_cup,t_cup
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        psur,z1
     real                                                              &
        ,intent (in   )                   ::                           &
        xl,rv,cp
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
!
!  local variables in this routine
!

     integer                              ::                           &
       i,k

     integer :: itf, ktf

     itf=ite
     ktf=kte

      do k=kts+1,ktf
      do i=its,itf
        if(ierr(i).eq.0)then
        qes_cup(i,k)=.5*(qes(i,k-1)+qes(i,k))
        q_cup(i,k)=.5*(q(i,k-1)+q(i,k))
        hes_cup(i,k)=.5*(hes(i,k-1)+hes(i,k))
        he_cup(i,k)=.5*(he(i,k-1)+he(i,k))
        if(he_cup(i,k).gt.hes_cup(i,k))he_cup(i,k)=hes_cup(i,k)
        z_cup(i,k)=.5*(z(i,k-1)+z(i,k))
        p_cup(i,k)=.5*(p(i,k-1)+p(i,k))
        t_cup(i,k)=.5*(t(i,k-1)+t(i,k))
        gamma_cup(i,k)=(xl/cp)*(xl/(rv*t_cup(i,k) &
                       *t_cup(i,k)))*qes_cup(i,k)
        endif
      enddo
      enddo
      do i=its,itf
        if(ierr(i).eq.0)then
        qes_cup(i,1)=qes(i,1)
        q_cup(i,1)=q(i,1)
        hes_cup(i,1)=hes(i,1)
        he_cup(i,1)=he(i,1)
        z_cup(i,1)=.5*(z(i,1)+z1(i))
        p_cup(i,1)=.5*(p(i,1)+psur(i))
        t_cup(i,1)=t(i,1)
        gamma_cup(i,1)=xl/cp*(xl/(rv*t_cup(i,1) &
                       *t_cup(i,1)))*qes_cup(i,1)
        endif
      enddo

   END SUBROUTINE cup_env_clev


   SUBROUTINE cup_forcing_ens(closure_n,xland,aa0,aa1,xaa0,mbdt,dtime,ierr,ierr2,ierr3,&
              xf_ens,j,name,maxens,iens,iedt,maxens2,maxens3,mconv,    &
              p_cup,ktop,omeg,zd,k22,zu,pr_ens,edt,kbcon,massflx,      &
              iact_old_gr,dir,ensdim,massfln,icoic,                    &
              ids,ide, jds,jde, kds,kde,                               &
              ims,ime, jms,jme, kms,kme,                               &
              its,ite, jts,jte, kts,kte                               )

   IMPLICIT NONE

     integer                                                           &
        ,intent (in   )                   ::                           &
        ids,ide, jds,jde, kds,kde,           &
        ims,ime, jms,jme, kms,kme,           &
        its,ite, jts,jte, kts,kte
     integer, intent (in   )              ::                           &
        j,ensdim,maxens,iens,iedt,maxens2,maxens3
  !
  ! ierr error value, maybe modified in this routine
  ! pr_ens = precipitation ensemble
  ! xf_ens = mass flux ensembles
  ! massfln = downdraft mass flux ensembles used in next timestep
  ! omeg = omega from large scale model
  ! mconv = moisture convergence from large scale model
  ! zd      = downdraft normalized mass flux
  ! zu      = updraft normalized mass flux
  ! aa0     = cloud work function without forcing effects
  ! aa1     = cloud work function with forcing effects
  ! xaa0    = cloud work function with cloud effects (ensemble dependent)
  ! edt     = epsilon
  ! dir     = "storm motion"
  ! mbdt    = arbitrary numerical parameter
  ! dtime   = dt over which forcing is applied
  ! iact_gr_old = flag to tell where convection was active
  ! kbcon       = LFC of parcel from k22
  ! k22         = updraft originating level
  ! icoic       = flag if only want one closure (usually set to zero!)
  ! name        = deep or shallow convection flag
  !
     real,    dimension (ims:ime,jms:jme,1:ensdim)                     &
        ,intent (inout)                   ::                           &
        pr_ens
     real,    dimension (ims:ime,jms:jme,1:ensdim)                     &
        ,intent (out  )                   ::                           &
        xf_ens,massfln
     real,    dimension (ims:ime,jms:jme)                              &
        ,intent (in   )                   ::                           &
        massflx
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        omeg,zd,zu,p_cup
     real,    dimension (its:ite,1:maxens)                             &
        ,intent (in   )                   ::                           &
        xaa0
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        aa1,edt,dir,mconv,xland
     real,    dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        aa0,closure_n
     real,    dimension (1:maxens)                                     &
        ,intent (in   )                   ::                           &
        mbdt
     real                                                              &
        ,intent (in   )                   ::                           &
        dtime
     integer, dimension (its:ite,jts:jte)                              &
        ,intent (in   )                   ::                           &
        iact_old_gr
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        k22,kbcon,ktop
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr,ierr2,ierr3
     integer                                                           &
        ,intent (in   )                   ::                           &
        icoic
      character *(*), intent (in)         ::                           &
       name
!
!  local variables in this routine
!

     real,    dimension (1:maxens3)       ::                           &
       xff_ens3
     real,    dimension (1:maxens)        ::                           &
       xk
     integer                              ::                           &
       i,k,nall,n,ne,nens,nens3,iresult,iresultd,iresulte,mkxcrt,kclim
     parameter (mkxcrt=15)
     real                                 ::                           &
       a1,massfld,xff0,xomg,aclim1,aclim2,aclim3,aclim4
     real,    dimension(1:mkxcrt)         ::                           &
       pcrit,acrit,acritt

     integer :: itf,nall2

     itf=ite

      DATA PCRIT/850.,800.,750.,700.,650.,600.,550.,500.,450.,400.,    &
                 350.,300.,250.,200.,150./
      DATA ACRIT/.0633,.0445,.0553,.0664,.075,.1082,.1521,.2216,       &
                 .3151,.3677,.41,.5255,.7663,1.1686,1.6851/
!  GDAS DERIVED ACRIT
      DATA ACRITT/.203,.515,.521,.566,.625,.665,.659,.688,             &
                  .743,.813,.886,.947,1.138,1.377,1.896/
!
       nens=0

!--- LARGE SCALE FORCING
!
       DO 100 i=its,itf
!       if(i.eq.ipr.and.j.eq.jpr)print *,'ierr = ',ierr(i)
          if(name.eq.'deeps'.and.ierr(i).gt.995)then
!          print *,i,j,ierr(i),aa0(i)
           aa0(i)=0.
           ierr(i)=0
          endif
          IF(ierr(i).eq.0)then
!           kclim=0
           do k=mkxcrt,1,-1
             if(p_cup(i,ktop(i)).lt.pcrit(k))then
               kclim=k
               go to 9
             endif
           enddo
           if(p_cup(i,ktop(i)).ge.pcrit(1))kclim=1
 9         continue
           kclim=max(kclim,1)
           k=max(kclim-1,1)
           aclim1=acrit(kclim)*1.e3
           aclim2=acrit(k)*1.e3
           aclim3=acritt(kclim)*1.e3
           aclim4=acritt(k)*1.e3
!           print *,'p_cup(ktop(i)),kclim,pcrit(kclim)'
!           print *,p_cup(i,ktop(i)),kclim,pcrit(kclim)
!           print *,'aclim1,aclim2,aclim3,aclim4'
!           print *,aclim1,aclim2,aclim3,aclim4
!           print *,dtime,name,ierr(i),aa1(i),aa0(i)
!          print *,dtime,name,ierr(i),aa1(i),aa0(i)
!
!--- treatment different for this closure
!
             if(name.eq.'deeps')then
!
                xff0= (AA1(I)-AA0(I))/DTIME
                xff_ens3(1)=(AA1(I)-AA0(I))/dtime
                xff_ens3(2)=.9*xff_ens3(1)
                xff_ens3(3)=1.1*xff_ens3(1)
!   
!--- more original Arakawa-Schubert (climatologic value of aa0)
!
!
!--- omeg is in bar/s, mconv done with omeg in Pa/s
!     more like Brown (1979), or Frank-Cohen (199?)
!
                xff_ens3(4)=-omeg(i,k22(i))/9.81
                xff_ens3(5)=-omeg(i,kbcon(i))/9.81
                xff_ens3(6)=-omeg(i,1)/9.81
                do k=2,kbcon(i)-1
                  xomg=-omeg(i,k)/9.81
                  if(xomg.gt.xff_ens3(6))xff_ens3(6)=xomg
                enddo
!
!--- more like Krishnamurti et al.
!
                xff_ens3(7)=mconv(i)
                xff_ens3(8)=mconv(i)
                xff_ens3(9)=mconv(i)
!
!--- more like Fritsch Chappel or Kain Fritsch (plus triggers)
!
                xff_ens3(10)=AA1(I)/(60.*20.)
                xff_ens3(11)=AA1(I)/(60.*30.)
                xff_ens3(12)=AA1(I)/(60.*40.)
!  
!--- more original Arakawa-Schubert (climatologic value of aa0)
!
                xff_ens3(13)=max(0.,(AA1(I)-aclim1)/dtime)
                xff_ens3(14)=max(0.,(AA1(I)-aclim2)/dtime)
                xff_ens3(15)=max(0.,(AA1(I)-aclim3)/dtime)
                xff_ens3(16)=max(0.,(AA1(I)-aclim4)/dtime)
!               if(ierr2(i).gt.0.or.ierr3(i).gt.0)then
!                 xff_ens3(10)=0.
!                 xff_ens3(11)=0.
!                 xff_ens3(12)=0.
!                 xff_ens3(13)=0.
!                 xff_ens3(14)=0.
!                 xff_ens3(15)=0.
!                 xff_ens3(16)=0.
!               endif

                do nens=1,maxens
                   XK(nens)=(XAA0(I,nens)-AA1(I))/MBDT(2)
                   if(xk(nens).le.0.and.xk(nens).gt.-1.e-6) &
                           xk(nens)=-1.e-6
                   if(xk(nens).gt.0.and.xk(nens).lt.1.e-6) &
                           xk(nens)=1.e-6
                enddo
!
!--- add up all ensembles
!
                do 350 ne=1,maxens
!
!--- for every xk, we have maxens3 xffs
!--- iens is from outermost ensemble (most expensive!
!
!--- iedt (maxens2 belongs to it)
!--- is from second, next outermost, not so expensive
!
!--- so, for every outermost loop, we have maxens*maxens2*3
!--- ensembles!!! nall would be 0, if everything is on first
!--- loop index, then ne would start counting, then iedt, then iens....
!
                   iresult=0
                   iresultd=0
                   iresulte=0
                   nall=(iens-1)*maxens3*maxens*maxens2 &
                        +(iedt-1)*maxens*maxens3 &
                        +(ne-1)*maxens3
!
! over water, enfor!e small cap for some of the closures
!
                if(xland(i).lt.-0.1)then
                 if(ierr2(i).gt.0.or.ierr3(i).gt.0)then
!       - ierr2 - 75 mb cap thickness, ierr3 - 125 cap thickness

! - for larger cap, set Grell closure to zero
                      xff_ens3(1) =0.
                      massfln(i,j,nall+1)=0.
                      xff_ens3(2) =0.
                      massfln(i,j,nall+2)=0.
                      xff_ens3(3) =0.
                      massfln(i,j,nall+3)=0.
                      closure_n(i)=closure_n(i)-1.

                      xff_ens3(7) =0.
                      massfln(i,j,nall+7)=0.
                      xff_ens3(8) =0.
                      massfln(i,j,nall+8)=0.
                      xff_ens3(9) =0.
!                     massfln(i,j,nall+9)=0.
                      closure_n(i)=closure_n(i)-1.
                endif
!
!   also take out some closures in general
!
                      xff_ens3(4) =0.
                      massfln(i,j,nall+4)=0.
                      xff_ens3(5) =0.
                      massfln(i,j,nall+5)=0.
                      xff_ens3(6) =0.
                      massfln(i,j,nall+6)=0.
                      closure_n(i)=closure_n(i)-3.

                      xff_ens3(10)=0.
                      massfln(i,j,nall+10)=0.
                      xff_ens3(11)=0.
                      massfln(i,j,nall+11)=0.
                      xff_ens3(12)=0.
                      massfln(i,j,nall+12)=0.
                      if(ne.eq.1)closure_n(i)=closure_n(i)-3
                      xff_ens3(13)=0.
                      massfln(i,j,nall+13)=0.
                      xff_ens3(14)=0.
                      massfln(i,j,nall+14)=0.
                      xff_ens3(15)=0.
                      massfln(i,j,nall+15)=0.
                      massfln(i,j,nall+16)=0.
                      if(ne.eq.1)closure_n(i)=closure_n(i)-4

                endif
!
! end water treatment
!
!--- check for upwind convection
!                  iresult=0
                   massfld=0.

!                  call cup_direction2(i,j,dir,iact_old_gr, &
!                       massflx,iresult,1,                  &
!                       massfld,                            &
!                       ids,ide, jds,jde, kds,kde,          &
!                       ims,ime, jms,jme, kms,kme,          &
!                       its,ite, jts,jte, kts,kte          )
!                  if(i.eq.ipr.and.j.eq.jpr.and.iedt.eq.1.and.ne.eq.1)then
!                  if(iedt.eq.1.and.ne.eq.1)then
!                   print *,massfld,ne,iedt,iens
!                   print *,xk(ne),xff_ens3(1),xff_ens3(2),xff_ens3(3)
!                  endif
!                  print *,i,j,massfld,aa0(i),aa1(i)
                   IF(XK(ne).lt.0.and.xff0.gt.0.)iresultd=1
                   iresulte=max(iresult,iresultd)
                   iresulte=1
                   if(iresulte.eq.1)then
!
!--- special treatment for stability closures
!

                      if(xff0.gt.0.)then
                         xf_ens(i,j,nall+1)=max(0.,-xff_ens3(1)/xk(ne)) &
                                        +massfld
                         xf_ens(i,j,nall+2)=max(0.,-xff_ens3(2)/xk(ne)) &
                                        +massfld
                         xf_ens(i,j,nall+3)=max(0.,-xff_ens3(3)/xk(ne)) &
                                        +massfld
                         xf_ens(i,j,nall+13)=max(0.,-xff_ens3(13)/xk(ne)) &
                                        +massfld
                         xf_ens(i,j,nall+14)=max(0.,-xff_ens3(14)/xk(ne)) &
                                        +massfld
                         xf_ens(i,j,nall+15)=max(0.,-xff_ens3(15)/xk(ne)) &
                                        +massfld
                         xf_ens(i,j,nall+16)=max(0.,-xff_ens3(16)/xk(ne)) &
                                        +massfld
                      else
                         xf_ens(i,j,nall+1)=massfld
                         xf_ens(i,j,nall+2)=massfld
                         xf_ens(i,j,nall+3)=massfld
                         xf_ens(i,j,nall+13)=massfld
                         xf_ens(i,j,nall+14)=massfld
                         xf_ens(i,j,nall+15)=massfld
                         xf_ens(i,j,nall+16)=massfld
                      endif
!
!--- if iresult.eq.1, following independent of xff0
!
                         xf_ens(i,j,nall+4)=max(0.,xff_ens3(4) &
                            +massfld)
                         xf_ens(i,j,nall+5)=max(0.,xff_ens3(5) &
                                        +massfld)
                         xf_ens(i,j,nall+6)=max(0.,xff_ens3(6) &
                                        +massfld)
                         a1=max(1.e-3,pr_ens(i,j,nall+7))
                         xf_ens(i,j,nall+7)=max(0.,xff_ens3(7) &
                                     /a1)
                         a1=max(1.e-3,pr_ens(i,j,nall+8))
                         xf_ens(i,j,nall+8)=max(0.,xff_ens3(8) &
                                     /a1)
                         a1=max(1.e-3,pr_ens(i,j,nall+9))
                         xf_ens(i,j,nall+9)=max(0.,xff_ens3(9) &
                                     /a1)
                         if(XK(ne).lt.0.)then
                            xf_ens(i,j,nall+10)=max(0., &
                                        -xff_ens3(10)/xk(ne)) &
                                        +massfld
                            xf_ens(i,j,nall+11)=max(0., &
                                        -xff_ens3(11)/xk(ne)) &
                                        +massfld
                            xf_ens(i,j,nall+12)=max(0., &
                                        -xff_ens3(12)/xk(ne)) &
                                        +massfld
                         else
                            xf_ens(i,j,nall+10)=massfld
                            xf_ens(i,j,nall+11)=massfld
                            xf_ens(i,j,nall+12)=massfld
                         endif
                      if(icoic.ge.1)then
                      closure_n(i)=0.
                      xf_ens(i,j,nall+1)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+2)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+3)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+4)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+5)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+6)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+7)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+8)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+9)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+10)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+11)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+12)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+13)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+14)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+15)=xf_ens(i,j,nall+icoic)
                      xf_ens(i,j,nall+16)=xf_ens(i,j,nall+icoic)
                      endif
!
! replace 13-16 for now with other stab closures
! (13 gave problems for mass model)
!
!                     xf_ens(i,j,nall+13)=xf_ens(i,j,nall+1)
                      if(icoic.eq.0)xf_ens(i,j,nall+14)=xf_ens(i,j,nall+13)
!                     xf_ens(i,j,nall+15)=xf_ens(i,j,nall+11)
!                     xf_ens(i,j,nall+16)=xf_ens(i,j,nall+11)
!                     xf_ens(i,j,nall+7)=xf_ens(i,j,nall+4)
!                     xf_ens(i,j,nall+8)=xf_ens(i,j,nall+5)
!                     xf_ens(i,j,nall+9)=xf_ens(i,j,nall+6)
!
!--- store new for next time step
!
                      do nens3=1,maxens3
                        massfln(i,j,nall+nens3)=edt(i) &
                                                *xf_ens(i,j,nall+nens3)
                        massfln(i,j,nall+nens3)=max(0., &
                                              massfln(i,j,nall+nens3))
                      enddo
!
!
!--- do some more on the caps!!! ne=1 for 175, ne=2 for 100,....
!
!     do not care for caps here for closure groups 1 and 5,
!     they are fine, do not turn them off here
!
!
                if(ne.eq.2.and.ierr2(i).gt.0)then
                      xf_ens(i,j,nall+1) =0.
                      xf_ens(i,j,nall+2) =0.
                      xf_ens(i,j,nall+3) =0.
                      xf_ens(i,j,nall+4) =0.
                      xf_ens(i,j,nall+5) =0.
                      xf_ens(i,j,nall+6) =0.
                      xf_ens(i,j,nall+7) =0.
                      xf_ens(i,j,nall+8) =0.
                      xf_ens(i,j,nall+9) =0.
                      xf_ens(i,j,nall+10)=0.
                      xf_ens(i,j,nall+11)=0.
                      xf_ens(i,j,nall+12)=0.
                      xf_ens(i,j,nall+13)=0.
                      xf_ens(i,j,nall+14)=0.
                      xf_ens(i,j,nall+15)=0.
                      xf_ens(i,j,nall+16)=0.
                      massfln(i,j,nall+1)=0.
                      massfln(i,j,nall+2)=0.
                      massfln(i,j,nall+3)=0.
                      massfln(i,j,nall+4)=0.
                      massfln(i,j,nall+5)=0.
                      massfln(i,j,nall+6)=0.
                      massfln(i,j,nall+7)=0.
                      massfln(i,j,nall+8)=0.
                      massfln(i,j,nall+9)=0.
                      massfln(i,j,nall+10)=0.
                      massfln(i,j,nall+11)=0.
                      massfln(i,j,nall+12)=0.
                      massfln(i,j,nall+13)=0.
                      massfln(i,j,nall+14)=0.
                      massfln(i,j,nall+15)=0.
                      massfln(i,j,nall+16)=0.
                endif
                if(ne.eq.3.and.ierr3(i).gt.0)then
                      xf_ens(i,j,nall+1) =0.
                      xf_ens(i,j,nall+2) =0.
                      xf_ens(i,j,nall+3) =0.
                      xf_ens(i,j,nall+4) =0.
                      xf_ens(i,j,nall+5) =0.
                      xf_ens(i,j,nall+6) =0.
                      xf_ens(i,j,nall+7) =0.
                      xf_ens(i,j,nall+8) =0.
                      xf_ens(i,j,nall+9) =0.
                      xf_ens(i,j,nall+10)=0.
                      xf_ens(i,j,nall+11)=0.
                      xf_ens(i,j,nall+12)=0.
                      xf_ens(i,j,nall+13)=0.
                      xf_ens(i,j,nall+14)=0.
                      xf_ens(i,j,nall+15)=0.
                      xf_ens(i,j,nall+16)=0.
                      massfln(i,j,nall+1)=0.
                      massfln(i,j,nall+2)=0.
                      massfln(i,j,nall+3)=0.
                      massfln(i,j,nall+4)=0.
                      massfln(i,j,nall+5)=0.
                      massfln(i,j,nall+6)=0.
                      massfln(i,j,nall+7)=0.
                      massfln(i,j,nall+8)=0.
                      massfln(i,j,nall+9)=0.
                      massfln(i,j,nall+10)=0.
                      massfln(i,j,nall+11)=0.
                      massfln(i,j,nall+12)=0.
                      massfln(i,j,nall+13)=0.
                      massfln(i,j,nall+14)=0.
                      massfln(i,j,nall+15)=0.
                      massfln(i,j,nall+16)=0.
                endif

                   endif
 350            continue
! ne=1, cap=175
!
                   nall=(iens-1)*maxens3*maxens*maxens2 &
                        +(iedt-1)*maxens*maxens3
! ne=2, cap=100
!
                   nall2=(iens-1)*maxens3*maxens*maxens2 &
                        +(iedt-1)*maxens*maxens3 &
                        +(2-1)*maxens3
                      xf_ens(i,j,nall+4) = xf_ens(i,j,nall2+4)
                      xf_ens(i,j,nall+5) =xf_ens(i,j,nall2+5)
                      xf_ens(i,j,nall+6) =xf_ens(i,j,nall2+6)
                      xf_ens(i,j,nall+7) =xf_ens(i,j,nall2+7)
                      xf_ens(i,j,nall+8) =xf_ens(i,j,nall2+8)
                      xf_ens(i,j,nall+9) =xf_ens(i,j,nall2+9)
                      xf_ens(i,j,nall+10)=xf_ens(i,j,nall2+10)
                      xf_ens(i,j,nall+11)=xf_ens(i,j,nall2+11)
                      xf_ens(i,j,nall+12)=xf_ens(i,j,nall2+12)
                go to 100
             endif
          elseif(ierr(i).ne.20.and.ierr(i).ne.0)then
             do n=1,ensdim
               xf_ens(i,j,n)=0.
               massfln(i,j,n)=0.
             enddo
          endif
 100   continue

   END SUBROUTINE cup_forcing_ens


   SUBROUTINE cup_kbcon(cap_inc,iloop,k22,kbcon,he_cup,hes_cup, &
              ierr,kbmax,p_cup,cap_max,                         &
              ids,ide, jds,jde, kds,kde,                        &
              ims,ime, jms,jme, kms,kme,                        &
              its,ite, jts,jte, kts,kte                        )

   IMPLICIT NONE
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
        ids,ide, jds,jde, kds,kde,           &
        ims,ime, jms,jme, kms,kme,           &
        its,ite, jts,jte, kts,kte
  ! 
  ! 
  ! 
  ! ierr error value, maybe modified in this routine
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        he_cup,hes_cup,p_cup
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        cap_max,cap_inc
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        kbmax
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        kbcon,k22,ierr
     integer                                                           &
        ,intent (in   )                   ::                           &
        iloop
!
!  local variables in this routine
!

     integer                              ::                           &
        i
     real                                 ::                           &
        pbcdif,plus
     integer :: itf

     itf=ite
!
!--- DETERMINE THE LEVEL OF CONVECTIVE CLOUD BASE  - KBCON
!
       DO 27 i=its,itf
      kbcon(i)=1
      IF(ierr(I).ne.0)GO TO 27
      KBCON(I)=K22(I)
      GO TO 32
 31   CONTINUE
      KBCON(I)=KBCON(I)+1
      IF(KBCON(I).GT.KBMAX(i)+4)THEN
         if(iloop.lt.4)ierr(i)=3
!        if(iloop.lt.4)ierr(i)=997
        GO TO 27
      ENDIF
 32   CONTINUE
      IF(HE_cup(I,K22(I)).LT.HES_cup(I,KBCON(I)))GO TO 31

!     cloud base pressure and max moist static energy pressure
!     i.e., the depth (in mb) of the layer of negative buoyancy
      if(KBCON(I)-K22(I).eq.1)go to 27
      PBCDIF=-P_cup(I,KBCON(I))+P_cup(I,K22(I))
      plus=max(25.,cap_max(i)-float(iloop-1)*cap_inc(i))
      if(iloop.eq.4)plus=cap_max(i)
      IF(PBCDIF.GT.plus)THEN
        K22(I)=K22(I)+1
        KBCON(I)=K22(I)
        GO TO 32
      ENDIF
 27   CONTINUE

   END SUBROUTINE cup_kbcon


   SUBROUTINE cup_kbcon_cin(iloop,k22,kbcon,he_cup,hes_cup,  &
              z,tmean,qes,ierr,kbmax,p_cup,cap_max,xl,cp,    &
              ids,ide, jds,jde, kds,kde,                     &
              ims,ime, jms,jme, kms,kme,                     &
              its,ite, jts,jte, kts,kte                     )

   IMPLICIT NONE
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
        ids,ide, jds,jde, kds,kde,           &
        ims,ime, jms,jme, kms,kme,           &
        its,ite, jts,jte, kts,kte
  ! 
  ! 
  ! 
  ! ierr error value, maybe modified in this routine
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        he_cup,hes_cup,p_cup,z,tmean,qes
     real,    dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        cap_max
     real                                                              &
        ,intent (in   )                   ::                           &
        xl,cp
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        kbmax
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        kbcon,k22,ierr
     integer                                                           &
        ,intent (in   )                   ::                           &
        iloop
!
!  local variables in this routine
!

     integer                              ::                           &
        i,k
     real                                 ::                           &
        cin,cin_max,dh,tprim,gamma
!
     integer :: itf

     itf=ite
!
!
    
!
!--- DETERMINE THE LEVEL OF CONVECTIVE CLOUD BASE  - KBCON
!
       DO 27 i=its,itf
      cin_max=-cap_max(i)
      kbcon(i)=1
      cin = 0.
      IF(ierr(I).ne.0)GO TO 27
      KBCON(I)=K22(I)
      GO TO 32
 31   CONTINUE
      KBCON(I)=KBCON(I)+1
      IF(KBCON(I).GT.KBMAX(i)+2)THEN
         if(iloop.eq.1)ierr(i)=3
!        if(iloop.eq.2)ierr(i)=997
        GO TO 27
      ENDIF
 32   CONTINUE
      dh      = HE_cup(I,K22(I)) - HES_cup(I,KBCON(I))
      if (dh.lt. 0.) then
        GAMMA=(xl/cp)*(xl/(461.525*(Tmean(I,K22(i))**2)))*QES(I,K22(i))
        tprim = dh/(cp*(1.+gamma))

        cin = cin + 9.8066 * tprim &
              *(z(i,k22(i))-z(i,k22(i)-1)) / tmean(i,k22(i))
        go to 31
      end if


!     If negative energy in negatively buoyant layer
!       exceeds convective inhibition (CIN) threshold,
!       then set K22 level one level up and see if that
!       will work.

      IF(cin.lT.cin_max)THEN
!       print *,i,cin,cin_max,k22(i),kbcon(i)
        K22(I)=K22(I)+1
        KBCON(I)=K22(I)
        GO TO 32
      ENDIF
 27   CONTINUE

   END SUBROUTINE cup_kbcon_cin


   SUBROUTINE cup_ktop(ilo,dby,kbcon,ktop,ierr,              &
              ids,ide, jds,jde, kds,kde,                     &
              ims,ime, jms,jme, kms,kme,                     &
              its,ite, jts,jte, kts,kte                     )

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
        ids,ide, jds,jde, kds,kde,           &
        ims,ime, jms,jme, kms,kme,           &
        its,ite, jts,jte, kts,kte
  ! dby = buoancy term
  ! ktop = cloud top (output)
  ! ilo  = flag
  ! ierr error value, maybe modified in this routine
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (inout)                   ::                           &
        dby
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        kbcon
     integer                                                           &
        ,intent (in   )                   ::                           &
        ilo
     integer, dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        ktop
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
!
!  local variables in this routine
!

     integer                              ::                           &
        i,k
!
     integer :: itf, ktf

     itf=ite
     ktf=kte
!
!
        DO 42 i=its,itf
        ktop(i)=1
         IF(ierr(I).EQ.0)then
          DO 40 K=KBCON(I)+1,ktf-1
            IF(DBY(I,K).LE.0.)THEN
                KTOP(I)=K-1
                GO TO 41
             ENDIF
  40      CONTINUE
          if(ilo.eq.1)ierr(i)=5
!         if(ilo.eq.2)ierr(i)=998
          GO TO 42
  41     CONTINUE
         do k=ktop(i)+1,ktf
           dby(i,k)=0.
         enddo
         endif
  42     CONTINUE

   END SUBROUTINE cup_ktop


   SUBROUTINE cup_MAXIMI(ARRAY,KS,KE,MAXX,ierr,              &
              ids,ide, jds,jde, kds,kde,                     &
              ims,ime, jms,jme, kms,kme,                     &
              its,ite, jts,jte, kts,kte                     )

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
         ids,ide, jds,jde, kds,kde,                                    &
         ims,ime, jms,jme, kms,kme,                                    &
         its,ite, jts,jte, kts,kte
  ! array input array
  ! x output array with return values
  ! kt output array of levels
  ! ks,kend  check-range
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
         array
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
         ierr,ke
     integer                                                           &
        ,intent (in   )                   ::                           &
         ks
     integer, dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
         maxx
     real,    dimension (its:ite)         ::                           &
         x
     real                                 ::                           &
         xar
     integer                              ::                           &
         i,k
     integer :: itf

     itf=ite

       DO 200 i=its,itf
       MAXX(I)=KS
       if(ierr(i).eq.0)then
      X(I)=ARRAY(I,KS)
!
       DO 100 K=KS,KE(i)
         XAR=ARRAY(I,K)
         IF(XAR.GE.X(I)) THEN
            X(I)=XAR
            MAXX(I)=K
         ENDIF
 100  CONTINUE
      endif
 200  CONTINUE

   END SUBROUTINE cup_MAXIMI


   SUBROUTINE cup_minimi(ARRAY,KS,KEND,KT,ierr,              &
              ids,ide, jds,jde, kds,kde,                     &
              ims,ime, jms,jme, kms,kme,                     &
              its,ite, jts,jte, kts,kte                     )

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
         ids,ide, jds,jde, kds,kde,                                    &
         ims,ime, jms,jme, kms,kme,                                    &
         its,ite, jts,jte, kts,kte
  ! array input array
  ! x output array with return values
  ! kt output array of levels
  ! ks,kend  check-range
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
         array
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
         ierr,ks,kend
     integer, dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
         kt
     real,    dimension (its:ite)         ::                           &
         x
     integer                              ::                           &
         i,k,kstop

     integer :: itf

     itf=ite

       DO 200 i=its,itf
      KT(I)=KS(I)
      if(ierr(i).eq.0)then
      X(I)=ARRAY(I,KS(I))
       KSTOP=MAX(KS(I)+1,KEND(I))
!
       DO 100 K=KS(I)+1,KSTOP
         IF(ARRAY(I,K).LT.X(I)) THEN
              X(I)=ARRAY(I,K)
              KT(I)=K
         ENDIF
 100  CONTINUE
      endif
 200  CONTINUE

   END SUBROUTINE cup_MINIMI


   SUBROUTINE cup_output_ens(xf_ens,ierr,dellat,dellaq,dellaqc,  &
              outtem,outq,outqc,pre,pw,xmb,ktop,                 &
              j,name,nx,nx2,iens,ierr2,ierr3,pr_ens,             &
              maxens3,ensdim,massfln,xmbmax,                     &
              APR_GR,APR_W,APR_MC,APR_ST,APR_AS,                 &
              APR_CAPMA,APR_CAPME,APR_CAPMI,closure_n,xland1,    &
              ids,ide, jds,jde, kds,kde, &
              ims,ime, jms,jme, kms,kme, &
              its,ite, jts,jte, kts,kte)

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
        ids,ide, jds,jde, kds,kde,           &
        ims,ime, jms,jme, kms,kme,           &
        its,ite, jts,jte, kts,kte
     integer, intent (in   )              ::                           &
        j,ensdim,nx,nx2,iens,maxens3
  ! xf_ens = ensemble mass fluxes
  ! pr_ens = precipitation ensembles
  ! dellat = change of temperature per unit mass flux of cloud ensemble
  ! dellaq = change of q per unit mass flux of cloud ensemble
  ! dellaqc = change of qc per unit mass flux of cloud ensemble
  ! outtem = output temp tendency (per s)
  ! outq   = output q tendency (per s)
  ! outqc  = output qc tendency (per s)
  ! pre    = output precip
  ! xmb    = total base mass flux
  ! xfac1  = correction factor
  ! pw = pw -epsilon*pd (ensemble dependent)
  ! ierr error value, maybe modified in this routine
  !
     real,    dimension (ims:ime,jms:jme,1:ensdim)                     &
        ,intent (inout)                   ::                           &
       xf_ens,pr_ens,massfln
     real,    dimension (ims:ime,jms:jme)                              &
        ,intent (inout)                   ::                           &
               APR_GR,APR_W,APR_MC,APR_ST,APR_AS,APR_CAPMA,            &
               APR_CAPME,APR_CAPMI 

     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
        outtem,outq,outqc
     real,    dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        pre,xmb
     real,    dimension (its:ite)                                      &
        ,intent (inout  )                   ::                           &
        closure_n,xland1
     real,    dimension (its:ite)                                      &
        ,intent (in     )                   ::                           &
        xmbmax
     real,    dimension (its:ite,kts:kte,1:nx)                         &
        ,intent (in   )                   ::                           &
       dellat,dellaqc,dellaq,pw
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        ktop
     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr,ierr2,ierr3
!
!  local variables in this routine
!

     integer                              ::                           &
        i,k,n,ncount
     real                                 ::                           &
        outtes,ddtes,dtt,dtq,dtqc,dtpw,tuning,prerate,clos_wei
     real,    dimension (its:ite)         ::                           &
       xfac1
     real,    dimension (its:ite)::                           &
       xmb_ske,xmb_ave,xmb_std,xmb_cur,xmbweight
     real,    dimension (its:ite)::                           &
       pr_ske,pr_ave,pr_std,pr_cur
     real,    dimension (its:ite,jts:jte)::                           &
               pr_gr,pr_w,pr_mc,pr_st,pr_as,pr_capma,     &
               pr_capme,pr_capmi

!
      character *(*), intent (in)        ::                           &
       name
!
     integer :: itf, ktf

     itf=ite
     ktf=kte
     tuning=0.2
!
!
      DO k=kts,ktf
      do i=its,itf
        outtem(i,k)=0.
        outq(i,k)=0.
        outqc(i,k)=0.
      enddo
      enddo
      do i=its,itf
        pre(i)=0.
        xmb(i)=0.
         xfac1(i)=1.
        xmbweight(i)=1.
      enddo
      do i=its,itf
        IF(ierr(i).eq.0)then
        do n=(iens-1)*nx*nx2*maxens3+1,iens*nx*nx2*maxens3
           if(pr_ens(i,j,n).le.0.)then
             xf_ens(i,j,n)=0.
           endif
        enddo
        endif
      enddo
!
!--- calculate ensemble average mass fluxes
!
       call massflx_stats(xf_ens,ensdim,nx2,nx,maxens3,      &
            xmb_ave,xmb_std,xmb_cur,xmb_ske,j,ierr,1,    &
            APR_GR,APR_W,APR_MC,APR_ST,APR_AS,           &
            APR_CAPMA,APR_CAPME,APR_CAPMI,               &
            pr_gr,pr_w,pr_mc,pr_st,pr_as,                &
            pr_capma,pr_capme,pr_capmi,                  &
            ids,ide, jds,jde, kds,kde,                   &
            ims,ime, jms,jme, kms,kme,                   &
            its,ite, jts,jte, kts,kte                   )
       call massflx_stats(pr_ens,ensdim,nx2,nx,maxens3,  &
            pr_ave,pr_std,pr_cur,pr_ske,j,ierr,2,        &
            APR_GR,APR_W,APR_MC,APR_ST,APR_AS,           &
            APR_CAPMA,APR_CAPME,APR_CAPMI,               &
            pr_gr,pr_w,pr_mc,pr_st,pr_as,                &
            pr_capma,pr_capme,pr_capmi,                  &
            ids,ide, jds,jde, kds,kde,                   &
            ims,ime, jms,jme, kms,kme,                   &
            its,ite, jts,jte, kts,kte                   )
!
!-- now do feedback
!
      ddtes=200.
!     if(name.eq.'shal')ddtes=200.
      do i=its,itf
        if(ierr(i).eq.0)then
         if(xmb_ave(i).le.0.)then
              ierr(i)=13
              xmb_ave(i)=0.
         endif
!        xmb(i)=max(0.,xmb_ave(i))
         xmb(i)=max(.1*xmb_ave(i),xmb_ave(i)-tuning*xmb_std(i))
!        xmb(i)=max(.1*xmb_ave(i),(1.+tuning)*xmb_ave(i))
! --- Now use proper count of how many closures were actually
!       used in cup_forcing_ens (including screening of some
!       closures over water) to properly normalize xmb
           clos_wei=16./max(1.,closure_n(i))
           if (xland1(i).lt.0.5)xmb(i)=xmb(i)*clos_wei
! make sure you take out more mass than what is there
           xmb(i)=min(xmb(i),xmbmax(i))
           if(xmb(i).eq.0.)then
              ierr(i)=19
           endif
           if(xmb(i).gt.100.)then
              ierr(i)=19
           endif
           xfac1(i)=xmb(i)

        endif
        xfac1(i)=xmb_ave(i)
      ENDDO
      DO k=kts,ktf
      do i=its,itf
            dtt=0.
            dtq=0.
            dtqc=0.
            dtpw=0.
        IF(ierr(i).eq.0.and.k.le.ktop(i))then
           do n=1,nx
              dtt=dtt+dellat(i,k,n)
              dtq=dtq+dellaq(i,k,n)
              dtqc=dtqc+dellaqc(i,k,n)
              dtpw=dtpw+pw(i,k,n)
           enddo
           outtes=dtt*XMB(I)*86400./float(nx)
           IF((OUTTES.GT.2.*ddtes.and.k.gt.2))THEN
             XMB(I)= 2.*ddtes/outtes * xmb(i)
             outtes=1.*ddtes
           endif
           if (outtes .lt. -ddtes) then
             XMB(I)= -ddtes/outtes * xmb(i)
             outtes=-ddtes
           endif
           if (outtes .gt. .5*ddtes.and.k.le.2) then
             XMB(I)= ddtes/outtes * xmb(i)
             outtes=.5*ddtes
           endif
           OUTTEM(I,K)=XMB(I)*dtt/float(nx)
           OUTQ(I,K)=XMB(I)*dtq/float(nx)
           OUTQC(I,K)=XMB(I)*dtqc/float(nx)
           PRE(I)=PRE(I)+XMB(I)*dtpw/float(nx)
        endif
      enddo
      enddo
!     do i=its,itf
!       if(ierr(i).eq.0)then
!          prerate=pre(i)*3600.
!          if(prerate.lt.0.1)then
!             if(ierr2(i).gt.0.or.ierr3(i).gt.0)then
!                pre(i)=0.
!                ierr(i)=221
!                do k=kts,ktf
!                   outtem(i,k)=0.
!                   outq(i,k)=0.
!                   outqc(i,k)=0.
!                enddo
!                do k=(iens-1)*nx*nx2*maxens3+1,iens*nx*nx2*maxens3
!                  massfln(i,j,k)=0.
!                  xf_ens(i,j,k)=0.
!                enddo
!              endif
!           endif

!       endif
!     ENDDO

      do i=its,itf
        if(ierr(i).eq.0)then
        xfac1(i)=xmb(i)/xfac1(i)
        do k=(iens-1)*nx*nx2*maxens3+1,iens*nx*nx2*maxens3
          massfln(i,j,k)=massfln(i,j,k)*xfac1(i)
          xf_ens(i,j,k)=xf_ens(i,j,k)*xfac1(i)
        enddo
        endif
      ENDDO

   END SUBROUTINE cup_output_ens


   SUBROUTINE cup_up_aa0(aa0,z,zu,dby,GAMMA_CUP,t_cup,       &
              kbcon,ktop,ierr,                               &
              ids,ide, jds,jde, kds,kde,                     &
              ims,ime, jms,jme, kms,kme,                     &
              its,ite, jts,jte, kts,kte                     )

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
        ids,ide, jds,jde, kds,kde,                                     &
        ims,ime, jms,jme, kms,kme,                                     &
        its,ite, jts,jte, kts,kte
  ! aa0 cloud work function
  ! gamma_cup = gamma on model cloud levels
  ! t_cup = temperature (Kelvin) on model cloud levels
  ! dby = buoancy term
  ! zu= normalized updraft mass flux
  ! z = heights of model levels 
  ! ierr error value, maybe modified in this routine
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        z,zu,gamma_cup,t_cup,dby
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        kbcon,ktop
!
! input and output
!


     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
     real,    dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        aa0
!
!  local variables in this routine
!

     integer                              ::                           &
        i,k
     real                                 ::                           &
        dz,da
!
     integer :: itf, ktf

     itf = ite
     ktf = kte

        do i=its,itf
         aa0(i)=0.
        enddo
        DO 100 k=kts+1,ktf
        DO 100 i=its,itf
         IF(ierr(i).ne.0)GO TO 100
         IF(K.LE.KBCON(I))GO TO 100
         IF(K.Gt.KTOP(I))GO TO 100
         DZ=Z(I,K)-Z(I,K-1)
         da=zu(i,k)*DZ*(9.81/(1004.*( &
                (T_cup(I,K)))))*DBY(I,K-1)/ &
             (1.+GAMMA_CUP(I,K))
         IF(K.eq.KTOP(I).and.da.le.0.)go to 100
         AA0(I)=AA0(I)+da
         if(aa0(i).lt.0.)aa0(i)=0.
100     continue

   END SUBROUTINE cup_up_aa0


   SUBROUTINE cup_up_he(k22,hkb,z_cup,cd,entr,he_cup,hc,     &
              kbcon,ierr,dby,he,hes_cup,                     &
           xlamue, xlamud,            &
              ids,ide, jds,jde, kds,kde,                     &
              ims,ime, jms,jme, kms,kme,                     &
              its,ite, jts,jte, kts,kte                     )

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
                                  ids,ide, jds,jde, kds,kde,           &
                                  ims,ime, jms,jme, kms,kme,           &
                                  its,ite, jts,jte, kts,kte
  ! hc = cloud moist static energy
  ! hkb = moist static energy at originating level
  ! he = moist static energy on model levels
  ! he_cup = moist static energy on model cloud levels
  ! hes_cup = saturation moist static energy on model cloud levels
  ! dby = buoancy term
  ! cd= detrainment function 
  ! z_cup = heights of model cloud levels 
  ! entr = entrainment rate
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        he,he_cup,hes_cup,z_cup,cd
  ! entr= entrainment rate 
     real,    dimension (its:ite,kts:kte)              &
        ,intent (in   )                   ::                           &
        xlamue
     real,    dimension (its:ite)                      &
        ,intent (in   )                   ::                           &
        xlamud                    

     real                                                              &
        ,intent (in   )                   ::                           &
        entr
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        kbcon,k22
!
! input and output
!

   ! ierr error value, maybe modified in this routine

     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr

     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
        hc,dby
     real,    dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        hkb
!
!  local variables in this routine
!

     integer                              ::                           &
        i,k
     real                                 ::           &
        dz,upentr,updetr
!
     integer :: itf, ktf

     itf = ite
     ktf = kte
!
!--- moist static energy inside cloud
!
      do i=its,itf
      if(ierr(i).eq.0.)then
      hkb(i)=he_cup(i,k22(i))
      do k=1,k22(i)
        hc(i,k)=he_cup(i,k)
        DBY(I,K)=0.
      enddo
      do k=k22(i),kbcon(i)-1
        hc(i,k)=hkb(i)
        DBY(I,K)=0.
      enddo
        k=kbcon(i)
        hc(i,k)=hkb(i)
        DBY(I,Kbcon(i))=Hkb(I)-HES_cup(I,K)
      endif
      enddo
      do k=kts+1,ktf
      do i=its,itf
        if(k.gt.kbcon(i).and.ierr(i).eq.0.)then
           DZ=Z_cup(i,K)-Z_cup(i,K-1)
!          upentr = 0.5 * (xlamue(i,k)+xlamue(i,k-1)) * dz
           upentr = xlamue(i,k) * dz
           updetr = 0.5 * xlamud(i) * dz
!          updetr = xlamud(i) * dz
!          HC(i,K)=(HC(i,K-1)*(1.-updetr)+    &
!                  upentr*0.5*(he(i,k)+HE(i,K-1)))/   &
!                  (1. + upentr - updetr)

           HC(i,K)=(HC(i,K-1)*(1.-updetr)+upentr* &
                HE(i,K-1))/(1.+upentr-updetr)
!          HC(i,K)=(HC(i,K-1)*(1.-.5*CD(i,K)*DZ)+entr* &
!               DZ*HE(i,K-1))/(1.+entr*DZ-.5*cd(i,k)*dz)
           DBY(I,K)=HC(I,K)-HES_cup(I,K)
        endif
      enddo

      enddo

   END SUBROUTINE cup_up_he


   SUBROUTINE cup_up_moisture(ierr,z_cup,qc,qrc,pw,pwav,     &
              kbcon,ktop,cd,dby,mentr_rate,clw_all,          &
              q,GAMMA_cup,zu,qes_cup,k22,qe_cup,xl,          &
           xlamue, xlamud,            &
              ids,ide, jds,jde, kds,kde,                     &
              ims,ime, jms,jme, kms,kme,                     &
              its,ite, jts,jte, kts,kte                     )

   IMPLICIT NONE
!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
                                  ids,ide, jds,jde, kds,kde,           &
                                  ims,ime, jms,jme, kms,kme,           &
                                  its,ite, jts,jte, kts,kte
  ! cd= detrainment function 
  ! q = environmental q on model levels
  ! qe_cup = environmental q on model cloud levels
  ! qes_cup = saturation q on model cloud levels
  ! dby = buoancy term
  ! cd= detrainment function 
  ! zu = normalized updraft mass flux
  ! gamma_cup = gamma on model cloud levels
  ! mentr_rate = entrainment rate
  !
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
        q,zu,gamma_cup,qe_cup,dby,qes_cup,z_cup,cd,xlamue
  ! entr= entrainment rate 
     real                                                              &
        ,intent (in   )                   ::                           &
        mentr_rate,xl
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
        kbcon,ktop,k22
     real,    dimension (its:ite)                      &
        ,intent (in   )                   ::                           &
        xlamud

!
! input and output
!

   ! ierr error value, maybe modified in this routine

     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
        ierr
   ! qc = cloud q (including liquid water) after entrainment
   ! qrch = saturation q in cloud
   ! qrc = liquid water content in cloud after rainout
   ! pw = condensate that will fall out at that level
   ! pwav = totan normalized integrated condensate (I1)
   ! c0 = conversion rate (cloud to rain)

     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
        qc,qrc,pw,clw_all
     real,    dimension (its:ite)                                      &
        ,intent (out  )                   ::                           &
        pwav
!
!  local variables in this routine
!

     integer                              ::                           &
        iall,i,k
     real                                 ::                           &
        dh,qrch,c0,dz,radius,tem,tem1,factor
!
     integer :: itf, ktf

     itf = ite
     ktf = kte

        iall=0
        c0=.002
!
!--- no precip for small clouds
!
        if(mentr_rate.gt.0.)then
          radius=.2/mentr_rate
          if(radius.lt.900.)c0=0.
!         if(radius.lt.900.)iall=0
        endif
        do i=its,itf
          pwav(i)=0.
        enddo
        do k=kts,ktf
        do i=its,itf
          pw(i,k)=0.
          if(ierr(i).eq.0)qc(i,k)=qes_cup(i,k)
          clw_all(i,k)=0.
          qrc(i,k)=0.
        enddo
        enddo
      do i=its,itf
      if(ierr(i).eq.0.)then
      do k=k22(i),kbcon(i)
        qc(i,k)=qe_cup(i,k22(i))
      enddo
      endif
      enddo

        DO 100 k=kts+1,ktf
        DO 100 i=its,itf
         IF(ierr(i).ne.0)GO TO 100
         IF(K.Le.KBCON(I))GO TO 100
         IF(K.Gt.KTOP(I))GO TO 100
         DZ=Z_cup(i,K)-Z_cup(i,K-1)
!           tem  = 0.5 * (xlamue(i,k)+xlamue(i,k-1)) * dz
            tem  = xlamue(i,k) * dz
            tem1 = 0.5 * xlamud(i) * dz
!           tem1 = xlamud(i) * dz
            factor = 1. + tem - tem1

!
!------    1. steady state plume equation, for what could
!------       be in cloud without condensation
!
!
!       QC(i,K)=(QC(i,K-1)*(1.-.5*CD(i,K)*DZ)+mentr_rate* &
!               DZ*Q(i,K-1))/(1.+mentr_rate*DZ-.5*cd(i,k)*dz)
!       QC(i,K)=(QC(i,K-1)*(1.-tem1)+tem &
!               *Q(i,K-1))/factor
! q without condensation
            qc(i,k) = ((1.-tem1)*qc(i,k-1) +       &
                      tem*q(i,k-1))   &
                      /factor

!
!--- saturation  in cloud, this is what is allowed to be in it
!
         QRCH=QES_cup(I,K)+(1./XL)*(GAMMA_cup(i,k) &
              /(1.+GAMMA_cup(i,k)))*DBY(I,K)
!
!------- LIQUID WATER CONTENT IN cloud after rainout
!
        clw_all(i,k)=QC(I,K)-QRCH
        QRC(I,K)=(QC(I,K)-QRCH)/(1.+C0*DZ*zu(i,k))
        if(qrc(i,k).lt.0.)then
          qrc(i,k)=0.
        endif
!
!-------   3.Condensation
!
         PW(i,k)=c0*dz*QRC(I,K)*zu(i,k)
        if(iall.eq.1)then
          qrc(i,k)=0.
          pw(i,k)=(QC(I,K)-QRCH)*zu(i,k)
          if(pw(i,k).lt.0.)pw(i,k)=0.
        endif
!
!----- set next level
!
         QC(I,K)=QRC(I,K)+qrch
!
!--- integrated normalized ondensate
!
         PWAV(I)=PWAV(I)+PW(I,K)
 100     CONTINUE

   END SUBROUTINE cup_up_moisture


   SUBROUTINE cup_up_nms(zu,z_cup,entr,cd,kbcon,ktop,ierr,k22,  &
           xlamue, xlamud,            &
              ids,ide, jds,jde, kds,kde,                        &
              ims,ime, jms,jme, kms,kme,                        &
              its,ite, jts,jte, kts,kte                        )

   IMPLICIT NONE

!
!  on input
!

   ! only local wrf dimensions are need as of now in this routine

     integer                                                           &
        ,intent (in   )                   ::                           &
         ids,ide, jds,jde, kds,kde,                                    &
         ims,ime, jms,jme, kms,kme,                                    &
         its,ite, jts,jte, kts,kte
  ! cd= detrainment function 
     real,    dimension (its:ite,kts:kte)                              &
        ,intent (in   )                   ::                           &
         z_cup,cd
  ! entr= entrainment rate 
     real                                                              &
        ,intent (in   )                   ::                           &
         entr
     integer, dimension (its:ite)                                      &
        ,intent (in   )                   ::                           &
         kbcon,ktop,k22
     real,    dimension (its:ite,kts:kte)              &
        ,intent (in   )                   ::                           &
        xlamue
     real,    dimension (its:ite)                      &
        ,intent (in   )                   ::                           &
        xlamud

!
! input and output
!

   ! ierr error value, maybe modified in this routine

     integer, dimension (its:ite)                                      &
        ,intent (inout)                   ::                           &
         ierr
   ! zu is the normalized mass flux

     real,    dimension (its:ite,kts:kte)                              &
        ,intent (out  )                   ::                           &
         zu
!
!  local variables in this routine
!

     integer                              ::                           &
         i,k
     real                                 ::                           &
         dz,ptem
     integer :: itf, ktf

     itf = ite
     ktf = kte
!
!   initialize for this go around
!
       do k=kts,ktf
       do i=its,itf
         zu(i,k)=0.
       enddo
       enddo
!
! do normalized mass budget
!
       do i=its,itf
          IF(ierr(I).eq.0)then
! subcloud levels: work downward from cloudbase
!            zu(i,kbcon(i))=1.
!            do k = kbcon(i)-1,k22(i),-1
!              DZ = Z_cup(i,K+1) - Z_cup(i,K)
!              ptem = 0.5*(xlamue(i,k)+xlamue(i,k+1)) - xlamud(i)
!              ptem = xlamue(i,k) - xlamud(i)
!              ZU(i,K) = ZU(i,K+1) / (1. + ptem*DZ)
!              zu(i,k)=1.
!            enddo
             do k=k22(i),kbcon(i)
               zu(i,k)=1.
             enddo
! levels above cloud base: work upward from cloudbase
!             DO K=KBcon(i)+1,KTOP(i)
             DO K = KBcon(i)+1,ktf-1
               DZ = Z_cup(i,K)-Z_cup(i,K-1)
!              ptem = 0.5*(xlamue(i,k)+xlamue(i,k-1)) - xlamud(i)
               ptem = xlamue(i,k) - xlamud(i)
               ZU(i,K) = ZU(i,K-1)*(1. + ptem*DZ)
             enddo

!            DO K=KBcon(i)+1,KTOP(i)
!              DZ=Z_cup(i,K)-Z_cup(i,K-1)
!              ZU(i,K)=ZU(i,K-1)*(1.+(entr-cd(i,k))*DZ)
!            enddo
          endif
       enddo

   END SUBROUTINE cup_up_nms

!====================================================================
   SUBROUTINE gdinit(RTHCUTEN,RQVCUTEN,RQCCUTEN,RQICUTEN,           &
                        MASS_FLUX,cp,restart,                       &
                        P_QC,P_QI,P_FIRST_SCALAR,                   &
                        RTHFTEN, RQVFTEN,                           &
                        APR_GR,APR_W,APR_MC,APR_ST,APR_AS,          &
                        APR_CAPMA,APR_CAPME,APR_CAPMI,              &
                        allowed_to_read,                            &
                        ids, ide, jds, jde, kds, kde,               &
                        ims, ime, jms, jme, kms, kme,               &
                        its, ite, jts, jte, kts, kte               )
!--------------------------------------------------------------------   
   IMPLICIT NONE
!--------------------------------------------------------------------
   LOGICAL , INTENT(IN)           ::  restart,allowed_to_read
   INTEGER , INTENT(IN)           ::  ids, ide, jds, jde, kds, kde, &
                                      ims, ime, jms, jme, kms, kme, &
                                      its, ite, jts, jte, kts, kte
   INTEGER , INTENT(IN)           ::  P_FIRST_SCALAR, P_QI, P_QC
   REAL,     INTENT(IN)           ::  cp

   REAL,     DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(OUT) ::       &
                                                          RTHCUTEN, &
                                                          RQVCUTEN, &
                                                          RQCCUTEN, &
                                                          RQICUTEN   

   REAL,     DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(OUT) ::       &
                                                          RTHFTEN,  &
                                                          RQVFTEN

   REAL,     DIMENSION( ims:ime , jms:jme ) , INTENT(OUT) ::        &
                                APR_GR,APR_W,APR_MC,APR_ST,APR_AS,  &
                                APR_CAPMA,APR_CAPME,APR_CAPMI,      &
                                MASS_FLUX

   IF(.not.restart)THEN
        RTHCUTEN=0.
        RQVCUTEN=0.
        RTHFTEN=0.
        RQVFTEN=0.

     IF (P_QC .ge. P_FIRST_SCALAR) THEN
           RQCCUTEN=0.
     ENDIF

     IF (P_QI .ge. P_FIRST_SCALAR) THEN
           RQICUTEN=0.
     ENDIF

        mass_flux=0.

   ENDIF
        APR_GR=0.
        APR_ST=0.
        APR_W=0.
        APR_MC=0.
        APR_AS=0.
        APR_CAPMA=0.
        APR_CAPME=0.
        APR_CAPMI=0.

   END SUBROUTINE gdinit


   SUBROUTINE massflx_stats(xf_ens,ensdim,maxens,maxens2,maxens3, &
              xt_ave,xt_std,xt_cur,xt_ske,j,ierr,itest,           &
              APR_GR,APR_W,APR_MC,APR_ST,APR_AS,                  &
              APR_CAPMA,APR_CAPME,APR_CAPMI,                      &
              pr_gr,pr_w,pr_mc,pr_st,pr_as,                       &
              pr_capma,pr_capme,pr_capmi,                         &
              ids,ide, jds,jde, kds,kde,  &
              ims,ime, jms,jme, kms,kme,  &
              its,ite, jts,jte, kts,kte)

   IMPLICIT NONE

   integer, intent (in   )              ::                                    &
                     j,ensdim,maxens3,maxens,maxens2,itest
   INTEGER,      INTENT(IN   ) ::                                             &
                                  ids,ide, jds,jde, kds,kde,                  &
                                  ims,ime, jms,jme, kms,kme,                  &
                                  its,ite, jts,jte, kts,kte


     real, dimension (its:ite)                                                &
         , intent(inout) ::                                                   &
           xt_ave,xt_cur,xt_std,xt_ske
     integer, dimension (its:ite), intent (in) ::                             &
           ierr
     real, dimension (ims:ime,jms:jme,1:ensdim)                               &
         , intent(in   ) ::                                                   &
           xf_ens
     real, dimension (ims:ime,jms:jme)                                        &
         , intent(inout) ::                                                   &
           APR_GR,APR_W,APR_MC,APR_ST,APR_AS,                                 &
           APR_CAPMA,APR_CAPME,APR_CAPMI
     real, dimension (its:ite,jts:jte)                                        &
         , intent(inout) ::                                                   &
           pr_gr,pr_w,pr_mc,pr_st,pr_as,                                      &
           pr_capma,pr_capme,pr_capmi

!
! local stuff
!
     real, dimension (its:ite , 1:maxens3 )       ::                          &
           x_ave,x_cur,x_std,x_ske
     real, dimension (its:ite , 1:maxens  )       ::                          &
           x_ave_cap


      integer, dimension (1:maxens3) :: nc1
      integer :: i,k
      integer :: num,kk,num2,iedt
      real :: a3,a4

      num=ensdim/maxens3
      num2=ensdim/maxens
      if(itest.eq.1)then
      do i=its,ite
       pr_gr(i,j) =  0.
       pr_w(i,j) =  0.
       pr_mc(i,j) = 0.
       pr_st(i,j) = 0.
       pr_as(i,j) = 0.
       pr_capma(i,j) =  0.
       pr_capme(i,j) = 0.
       pr_capmi(i,j) = 0.
      enddo
      endif

      do k=1,maxens
      do i=its,ite
        x_ave_cap(i,k)=0.
      enddo
      enddo
      do k=1,maxens3
      do i=its,ite
        x_ave(i,k)=0.
        x_std(i,k)=0.
        x_ske(i,k)=0.
        x_cur(i,k)=0.
      enddo
      enddo
      do i=its,ite
        xt_ave(i)=0.
        xt_std(i)=0.
        xt_ske(i)=0.
        xt_cur(i)=0.
      enddo
      do kk=1,num
      do k=1,maxens3
      do i=its,ite
        if(ierr(i).eq.0)then
        x_ave(i,k)=x_ave(i,k)+xf_ens(i,j,maxens3*(kk-1)+k)
        endif
      enddo
      enddo
      enddo
      do iedt=1,maxens2
      do k=1,maxens
      do kk=1,maxens3
      do i=its,ite
        if(ierr(i).eq.0)then
        x_ave_cap(i,k)=x_ave_cap(i,k)                               &
            +xf_ens(i,j,maxens3*(k-1)+(iedt-1)*maxens*maxens3+kk)
        endif
      enddo
      enddo
      enddo
      enddo
      do k=1,maxens
      do i=its,ite
        if(ierr(i).eq.0)then
        x_ave_cap(i,k)=x_ave_cap(i,k)/float(num2)
        endif
      enddo
      enddo

      do k=1,maxens3
      do i=its,ite
        if(ierr(i).eq.0)then
        x_ave(i,k)=x_ave(i,k)/float(num)
        endif
      enddo
      enddo
      do k=1,maxens3
      do i=its,ite
        if(ierr(i).eq.0)then
        xt_ave(i)=xt_ave(i)+x_ave(i,k)
        endif
      enddo
      enddo
      do i=its,ite
        if(ierr(i).eq.0)then
        xt_ave(i)=xt_ave(i)/float(maxens3)
        endif
      enddo
!
!--- now do std, skewness,curtosis
!
      do kk=1,num
      do k=1,maxens3
      do i=its,ite
        if(ierr(i).eq.0.and.x_ave(i,k).gt.0.)then
!       print *,i,j,k,kk,x_std(i,k),xf_ens(i,j,maxens3*(kk-1)+k),x_ave(i,k)
        x_std(i,k)=x_std(i,k)+(xf_ens(i,j,maxens3*(kk-1)+k)-x_ave(i,k))**2
        x_ske(i,k)=x_ske(i,k)+(xf_ens(i,j,maxens3*(kk-1)+k)-x_ave(i,k))**3
        x_cur(i,k)=x_cur(i,k)+(xf_ens(i,j,maxens3*(kk-1)+k)-x_ave(i,k))**4
        endif
      enddo
      enddo
      enddo
      do k=1,maxens3
      do i=its,ite
        if(ierr(i).eq.0.and.xt_ave(i).gt.0.)then
        xt_std(i)=xt_std(i)+(x_ave(i,k)-xt_ave(i))**2
        xt_ske(i)=xt_ske(i)+(x_ave(i,k)-xt_ave(i))**3
        xt_cur(i)=xt_cur(i)+(x_ave(i,k)-xt_ave(i))**4
        endif
      enddo
      enddo
      do k=1,maxens3
      do i=its,ite
        if(ierr(i).eq.0.and.x_std(i,k).gt.0.)then
           x_std(i,k)=x_std(i,k)/float(num)
           a3=max(1.e-6,x_std(i,k))
           x_std(i,k)=sqrt(a3)
           a3=max(1.e-6,x_std(i,k)**3)
           a4=max(1.e-6,x_std(i,k)**4)
           x_ske(i,k)=x_ske(i,k)/float(num)/a3
           x_cur(i,k)=x_cur(i,k)/float(num)/a4
        endif
!       print*,'                               '
!       print*,'Some statistics at gridpoint i,j, ierr',i,j,ierr(i)
!       print*,'statistics for closure number ',k
!       print*,'Average= ',x_ave(i,k),'  Std= ',x_std(i,k)
!       print*,'Skewness= ',x_ske(i,k),' Curtosis= ',x_cur(i,k)
!       print*,'                               '

      enddo
      enddo
      do i=its,ite
        if(ierr(i).eq.0.and.xt_std(i).gt.0.)then
           xt_std(i)=xt_std(i)/float(maxens3)
           a3=max(1.e-6,xt_std(i))
           xt_std(i)=sqrt(a3)
           a3=max(1.e-6,xt_std(i)**3)
           a4=max(1.e-6,xt_std(i)**4)
           xt_ske(i)=xt_ske(i)/float(maxens3)/a3
           xt_cur(i)=xt_cur(i)/float(maxens3)/a4
!       print*,'                               '
!       print*,'Total ensemble independent statistics at i =',i
!       print*,'Average= ',xt_ave(i),'  Std= ',xt_std(i)
!       print*,'Skewness= ',xt_ske(i),' Curtosis= ',xt_cur(i)
!       print*,'                               '
!
!  first go around: store massflx for different closures/caps
!
      if(itest.eq.1)then
       pr_gr(i,j) = .333*(x_ave(i,1)+x_ave(i,2)+x_ave(i,3))
       pr_w(i,j) = .333*(x_ave(i,4)+x_ave(i,5)+x_ave(i,6))
       pr_mc(i,j) = .333*(x_ave(i,7)+x_ave(i,8)+x_ave(i,9))
       pr_st(i,j) = .333*(x_ave(i,10)+x_ave(i,11)+x_ave(i,12))
       pr_as(i,j) = .25*(x_ave(i,13)+x_ave(i,14)+x_ave(i,15) &
                     + x_ave(i,16))
       pr_capma(i,j) = x_ave_cap(i,1)
       pr_capme(i,j) = x_ave_cap(i,2)
       pr_capmi(i,j) = x_ave_cap(i,3)
!
!  second go around: store preciprates (mm/hour) for different closures/caps
!
        else if (itest.eq.2)then
       APR_GR(i,j)=.333*(x_ave(i,1)+x_ave(i,2)+x_ave(i,3))*      &
                  3600.*pr_gr(i,j) +APR_GR(i,j)
       APR_W(i,j)=.333*(x_ave(i,4)+x_ave(i,5)+x_ave(i,6))*       &
                  3600.*pr_w(i,j) +APR_W(i,j)
       APR_MC(i,j)=.333*(x_ave(i,7)+x_ave(i,8)+x_ave(i,9))*      &
                  3600.*pr_mc(i,j) +APR_MC(i,j)
       APR_ST(i,j)=.333*(x_ave(i,10)+x_ave(i,11)+x_ave(i,12))*   &
                  3600.*pr_st(i,j) +APR_ST(i,j)
       APR_AS(i,j)=.25*(x_ave(i,13)+x_ave(i,14)+x_ave(i,15)      &
                           + x_ave(i,16))*                       &
                  3600.*pr_as(i,j) +APR_AS(i,j)
       APR_CAPMA(i,j) = x_ave_cap(i,1)*                          &
                  3600.*pr_capma(i,j) +APR_CAPMA(i,j)
       APR_CAPME(i,j) = x_ave_cap(i,2)*                          &
                  3600.*pr_capme(i,j) +APR_CAPME(i,j)
       APR_CAPMI(i,j) = x_ave_cap(i,3)*                          &
                  3600.*pr_capmi(i,j) +APR_CAPMI(i,j)
        endif
        endif
      enddo

   END SUBROUTINE massflx_stats


   SUBROUTINE neg_check(dt,q,outq,outt,outqc,pret,its,ite,kts,kte,itf,ktf)

   INTEGER,      INTENT(IN   ) ::            its,ite,kts,kte,itf,ktf

     real, dimension (its:ite,kts:kte  )                    ,                 &
      intent(inout   ) ::                                                     &
       q,outq,outt,outqc
     real, dimension (its:ite  )                            ,                 &
      intent(inout   ) ::                                                     &
       pret
     real                                                                     &
        ,intent (in  )                   ::                                   &
        dt
     real :: thresh,qmem,qmemf,qmem2,qtest,qmem1
!
! first do check on vertical heating rate
!
      thresh=200.01
      do i=its,itf
      qmemf=1.
      qmem=0.
      do k=kts,ktf
         qmem=outt(i,k)*86400.
         if(qmem.gt.2.*thresh)then
           qmem2=2.*thresh/qmem
           qmemf=min(qmemf,qmem2)
!
!
!          print *,'1',' adjusted massflux by factor ',i,k,qmem,qmem2,qmemf
         endif
         if(qmem.lt.-thresh)then
           qmem2=-thresh/qmem
           qmemf=min(qmemf,qmem2)
!
!
!          print *,'2',' adjusted massflux by factor ',i,k,qmem,qmem2,qmemf
         endif
      enddo
      do k=kts,ktf
         outq(i,k)=outq(i,k)*qmemf
         outt(i,k)=outt(i,k)*qmemf
         outqc(i,k)=outqc(i,k)*qmemf
      enddo
      pret(i)=pret(i)*qmemf 
      enddo
!
! check whether routine produces negative q's. This can happen, since 
! tendencies are calculated based on forced q's. This should have no
! influence on conservation properties, it scales linear through all
! tendencies
!
      thresh=1.e-10
      do i=its,itf
      qmemf=1.
      do k=kts,ktf
         qmem=outq(i,k)
         if(abs(qmem).gt.0.)then
         qtest=q(i,k)+outq(i,k)*dt
         if(qtest.lt.thresh)then
!
! qmem2 would be the maximum allowable tendency
!
           qmem1=outq(i,k)
           qmem2=(thresh-q(i,k))/dt
           qmemf=min(qmemf,qmem2/qmem1)
!          qmemf=max(0.,qmemf)
!          print *,'4 adjusted tendencies ',i,k,qmem,qmem2,qmemf
         endif
         endif
      enddo
      do k=kts,ktf
         outq(i,k)=outq(i,k)*qmemf
         outt(i,k)=outt(i,k)*qmemf
         outqc(i,k)=outqc(i,k)*qmemf
      enddo
      pret(i)=pret(i)*qmemf 
      enddo

   END SUBROUTINE neg_check


!-------------------------------------------------------
END MODULE module_cu_gd
