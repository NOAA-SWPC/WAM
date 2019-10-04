      SUBROUTINE DELDIFS_noq
     &                  (WE,QME,XE,YE,TEME,
     &                   WO,QMO,XO,YO,TEMO,DELTIM,SL,
     &                   LS_NODE,COEF00,K_LEVEL,
     &                   hybrid,gen_coord_hybrid)
!
! program log
! 2011 02 20 : henry juang, diffusion without tracer for ndsl advection
!
      use gfs_dyn_resol_def
      use gfs_dyn_layout1
      use gfs_dyn_coordinate_def				! hmhj
      use gfs_dyn_deldifs_def
      use gfs_dyn_physcons, rerth => con_rerth
     &                    ,  rd => con_rd, cp => con_cp
      IMPLICIT NONE
!
      logical hybrid, gen_coord_hybrid
      REAL(KIND=KIND_EVOD)  WE(LEN_TRIE_LS,2)
     &,                    QME(LEN_TRIE_LS,2)
     &,                     XE(LEN_TRIE_LS,2)
     &,                     YE(LEN_TRIE_LS,2)
     &,                   TEME(LEN_TRIE_LS,2)
     &,                     PE(LEN_TRIE_LS,2)
!
      REAL(KIND=KIND_EVOD)  WO(LEN_TRIO_LS,2)
     &,                    QMO(LEN_TRIO_LS,2)
     &,                     XO(LEN_TRIO_LS,2)
     &,                     YO(LEN_TRIO_LS,2)
     &,                   TEMO(LEN_TRIO_LS,2)
     &,                     PO(LEN_TRIO_LS,2)
!
      REAL(KIND=KIND_EVOD) DELTIM, SL(LEVS)
!
      INTEGER              LS_NODE(LS_DIM,3)
!
!CMR  LS_NODE(1,1) ... LS_NODE(LS_MAX_NODE,1) : VALUES OF L
!CMR  LS_NODE(1,2) ... LS_NODE(LS_MAX_NODE,2) : VALUES OF JBASEV
!CMR  LS_NODE(1,3) ... LS_NODE(LS_MAX_NODE,3) : VALUES OF JBASOD
!
      REAL(KIND=KIND_EVOD) COEF00(LEVS,NTRAC)
!
      INTEGER              K_LEVEL
!
      INTEGER              I,IS,IT,JDEL,JDELH,K,KD,KU
      INTEGER              L,LOCL,N,N0,ND,NP,NPD
!
      INTEGER              INDEV
      INTEGER              INDOD
      integer              indev1,indev2
      integer              indod1,indod2
      real(kind=kind_evod), parameter :: rkappa = cp / rd
      REAL(KIND=KIND_EVOD) DN1,REALVAL,RTNP,DF_DK,FACT
     &,                    SLRD0,FTRD1,RFACT,RFACTRD,RTRD1,FSHK, tem
!
      REAL(KIND=KIND_EVOD), parameter :: CONS0=0.0, CONS1=1.0, CONS2=2.0
!
      INTEGER              INDLSEV,JBASEV
      INTEGER              INDLSOD,JBASOD
!
      INCLUDE 'function2'
!
!     print *,' enter deldifs_fd ' 					! hmhj
!......................................................................
!
      IF(K_LEVEL.EQ.0) THEN
!
        CALL COUNTPERF(0,15,0.)
!!
        allocate(RTRD(LEVS),RTHK(LEVS),sf(levs))
        ALLOCATE ( DNE(LEN_TRIE_LS) )
        ALLOCATE ( DNO(LEN_TRIO_LS) )
        ALLOCATE ( BKLY(levs) )        					! hmhj
        ALLOCATE ( CKLY(levs) )        					! hmhj
        BKLY(:) = 1.0
        CKLY(:) = 0.0
        if (gen_coord_hybrid) then					! hmhj
          DO  k=1,LEVS							! hmhj
! hmhj ak5, bk5, ck5 in gen_coord_hybrid is the same order as model index
            BKLY(k)=0.5*(bk5(k)+bk5(k+1))				! hmhj
            CKLY(k)=0.5*(ck5(k)+ck5(k+1))*rkappa/thref(k)	        ! hmhj
!           if( me.eq.0 )						! hmhj
!    &         print*,'sl bkly ckly  in deldif=',k,sl(k),bkly(k),ckly(k)! hmhj
          enddo								! hmhj
        else if (hybrid) then						! hmhj
          DO  k=1,LEVS
! hmhj   sl(k) go bottom to top but bk(k) go top to bottom
            BKLY(k)=0.5*(bk5(levs-k+1)+bk5(levs-k+2))/SL(k)
!           if( me.eq.0 ) print*,'sl bkly in deldif=',k,sl(k),bkly(k)
          enddo
        endif
!
        IF(JCAP.GT.170) THEN
!         RECIPROCAL OF TIME SCALE OF DIFFUSION AT REFERENCE WAVENUMBER NP
!         RTNP=(JCAP/170.)**4*1.1/3600
!         RTNP=(JCAP/170.)**4*1.1/3600./24.0
          RTNP=(JCAP/170.)**4*1.1/3600./24000.0
          NP=JCAP
          N0=0             ! MAXIMUM WAVENUMBER FOR ZERO DIFFUSION
          JDEL=8           ! ORDER OF DIFFUSION (EVEN POWER TO RAISE DEL)
          FSHK=2.2         ! EXTRA HEIGHT-DEPENDENT DIFFUSION FACTOR PER SCALE HEIGHT
        ELSEIF(JCAP.EQ.170) THEN
!         RECIPROCAL OF TIME SCALE OF DIFFUSION AT REFERENCE WAVENUMBER NP
          RTNP=4*3.E15/(RERTH**4)*FLOAT(80*81)**2
          NP=JCAP
          N0=0.55*JCAP     ! MAXIMUM WAVENUMBER FOR ZERO DIFFUSION
          JDEL=2           ! ORDER OF DIFFUSION (EVEN POWER TO RAISE DEL)
          FSHK=1.0         ! EXTRA HEIGHT-DEPENDENT DIFFUSION FACTOR PER SCALE HEIGHT
        ELSEIF(JCAP.EQ.126) THEN					! hmhj
!         BELOW HAS BEEN TESTED IN SIHMA-THETA FOR 2 YEAR CFS RUN	! hmhj
          RTNP=4*3.E15/(RERTH**4)*FLOAT(80*81)**2			! hmhj
          NP=JCAP							! hmhj
          N0=0.0           						! hmhj
          JDEL=4           						! hmhj
          FSHK=1.0         						! hmhj
        ELSE
!         RECIPROCAL OF TIME SCALE OF DIFFUSION AT REFERENCE WAVENUMBER NP
          RTNP=1*3.E15/(RERTH**4)*FLOAT(80*81)**2
          NP=JCAP
          N0=0.55*JCAP     ! MAXIMUM WAVENUMBER FOR ZERO DIFFUSION
          JDEL=2           ! ORDER OF DIFFUSION (EVEN POWER TO RAISE DEL)
          FSHK=1.0         ! EXTRA HEIGHT-DEPENDENT DIFFUSION FACTOR PER SCALE HEIGHT
        ENDIF
!
        IF (ME.EQ.0) THEN
          PRINT 6,RTNP,NP,N0,JDEL
    6     FORMAT(' IN DELDIF_NOQ: HORIZONTAL DIFFUSION PARAMETERS'/
     &  '   EFFECTIVE ',6PF10.3,' MICROHERTZ AT WAVENUMBER ',I4/
     &  '   MAXIMUM WAVENUMBER FOR ZERO DIFFUSION ',I4/
     &  '   ORDER OF DIFFUSION ',I2)
        ENDIF
!
        SLRD0=0.002        ! SIGMA LEVEL AT WHICH TO BEGIN RAYLEIGH DAMPING
        RTRD1 = 0.0        !Rayleigh friction is now done in gbphys
                           !on the grid using rayleigh_damp.f
        if(adiabatic) RTRD1=1./(5*86400) ! RECIPROCAL OF TIME SCALE PER SCALE HEIGHT
                           !  ABOVE BEGINNING SIGMA LEVEL FOR RAYLEIGH DAMPING
!
        DO K=1,LEVS
          IF(SL(K).LT.SLRD0) THEN
            if (k .gt. levr) then
              RTRD(K)=RTRD1*LOG(SLRD0/SL(K)) ** 2
            else
              RTRD(K)=RTRD1*LOG(SLRD0/SL(K))
            endif
          ELSE
            RTRD(K)=0
          ENDIF
          RTHK(K)=(SL(K))**LOG(1/FSHK)
        ENDDO
!
        JDELH=JDEL/2
        NPD=MAX(NP-N0,0)
        REALVAL=NPD*(NPD+1)
        DN1=CONS2*RTNP/REALVAL**JDELH
!
!......................................................................
!
        DO LOCL=1,LS_MAX_NODE
               L=LS_NODE(LOCL,1)
          JBASEV=LS_NODE(LOCL,2)
          INDEV=INDLSEV(L,L)
          DO N=L,JCAP,2
            ND=MAX(N-N0,0)
            REALVAL=ND*(ND+1)
            DNE(INDEV)=DN1*REALVAL**JDELH
            INDEV=INDEV+1
          ENDDO
        ENDDO
!
!......................................................................
!
        DO LOCL=1,LS_MAX_NODE
               L=LS_NODE(LOCL,1)
          JBASEV=LS_NODE(LOCL,2)
          if (mod(L,2).eq.mod(jcap+1,2)) then
            DNE(INDLSEV(JCAP+1,L))=CONS0 ! SET THE EVEN (N-L) TERMS OF THE TOP ROW TO ZERO
          ENDIF
        ENDDO
!
!......................................................................
!
        DO LOCL=1,LS_MAX_NODE
               L=LS_NODE(LOCL,1)
          JBASOD=LS_NODE(LOCL,3)
          INDOD=INDLSOD(L+1,L)
          DO N=L+1,JCAP,2
            ND=MAX(N-N0,0)
            REALVAL=ND*(ND+1)
            DNO(INDOD)=DN1*REALVAL**JDELH
            INDOD=INDOD+1
          ENDDO
        ENDDO
!
!......................................................................
!
        DO LOCL=1,LS_MAX_NODE
               L=LS_NODE(LOCL,1)
          JBASOD=LS_NODE(LOCL,3)
          if (mod(L,2).ne.mod(jcap+1,2)) then
            DNO(INDLSOD(JCAP+1,L))=CONS0 ! SET THE ODD (N-L) TERMS OF THE TOP ROW TO ZERO
          ENDIF
        ENDDO
!
!......................................................................
!
        DO K=1,LEVS
          KD=MAX(K-1,1)
          KU=MIN(K+1,LEVS)
          SF(K)=SL(K)/(SL(KU)-SL(KD))/SQRT(CONS2)     !CONSTANT
        ENDDO
!
        CALL COUNTPERF(1,15,0.)
!!
        RETURN
      ENDIF
!
!......................................................................
!
      CALL COUNTPERF(0,13,0.)
!!
      K=K_LEVEL
!!
!
!     TEM = COEF00(K,1) * BKLY(K)
!
      DO LOCL=1,LS_MAX_NODE
              L=LS_NODE(LOCL,1)
         JBASEV=LS_NODE(LOCL,2)
         IF (L.EQ.0) THEN
            N0=2
         ELSE
            N0=L
         ENDIF
         indev1 = indlsev(N0,L)
         if (mod(L,2).eq.mod(jcap+1,2)) then
           indev2 = indlsev(jcap+1,L)
         else
           indev2 = indlsev(jcap  ,L)
         endif
!!       DO N = N0, JCAP+1, 2
         DO INDEV = indev1 , indev2
 
           FACT             = DELTIM*DNE(INDEV)*RTHK(K)
           RFACT            = CONS1/(CONS1+FACT)
           RFACTRD          = CONS1/(CONS1+FACT+DELTIM*RTRD(K))
 
           WE(INDEV,1)      = WE(INDEV,1)*RFACTRD
           WE(INDEV,2)      = WE(INDEV,2)*RFACTRD
 
           XE(INDEV,1)      = XE(INDEV,1)*RFACTRD
           XE(INDEV,2)      = XE(INDEV,2)*RFACTRD
 
           PE(INDEV,1)=BKLY(K)*QME(INDEV,1)+CKLY(K)*TEME(INDEV,1)	! hmhj
           PE(INDEV,2)=BKLY(K)*QME(INDEV,2)+CKLY(K)*TEME(INDEV,2)	! hmhj
 
           YE(INDEV,1)      =  ( YE(INDEV,1)+
     X                 FACT*COEF00(K,1)* PE(INDEV,1) )*RFACT		! hmhj
 
           YE(INDEV,2)      = ( YE(INDEV,2) +
     X                 FACT*COEF00(K,1)* PE(INDEV,2) )*RFACT		! hmhj
 
         ENDDO
       ENDDO
!
!......................................................................
!
!      DO L = 0, JCAP
       DO LOCL=1,LS_MAX_NODE
              L=LS_NODE(LOCL,1)
         JBASOD=LS_NODE(LOCL,3)
         indod1 = indlsod(L+1,L)
         if (mod(L,2).eq.mod(jcap+1,2)) then
           indod2 = indlsod(jcap  ,L)
         else
           indod2 = indlsod(jcap+1,L)
         endif
!        DO N = L+1, JCAP+1, 2
         DO INDOD = indod1 , indod2
 
           FACT             = DELTIM*DNO(INDOD)*RTHK(K)
           RFACT            = CONS1/(CONS1+FACT)
           RFACTRD          = CONS1/(CONS1+FACT+DELTIM*RTRD(K))
 
           WO(INDOD,1)      = WO(INDOD,1)*RFACTRD
           WO(INDOD,2)      = WO(INDOD,2)*RFACTRD
 
           XO(INDOD,1)      = XO(INDOD,1)*RFACTRD
           XO(INDOD,2)      = XO(INDOD,2)*RFACTRD

           PO(INDOD,1)=BKLY(K)*QMO(INDOD,1)+CKLY(K)*TEMO(INDOD,1)	! hmhj
           PO(INDOD,2)=BKLY(K)*QMO(INDOD,2)+CKLY(K)*TEMO(INDOD,2)	! hmhj
 
           YO(INDOD,1)      = ( YO(INDOD,1)+
     X                 FACT*COEF00(K,1)* PO(INDOD,1) )*RFACT		! hmhj
 
           YO(INDOD,2)      = ( YO(INDOD,2)+
     X                 FACT*COEF00(K,1)* PO(INDOD,2) )*RFACT		! hmhj
 
         ENDDO
       ENDDO
!
!......................................................................
!
!
      CALL COUNTPERF(1,13,0.)

!     print *,' leave deldifs_fd ' 					! hmhj
!!
      RETURN
      END
