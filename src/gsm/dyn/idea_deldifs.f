!raa********************************************************************
! Starting May 27, 2008: modified by RAA to impose horizontal physical
!     dissipation on top of numerical diffusion in all layers.
! Jan 4, 2007: modified by Rashid Akmaev based on DELDIFS_hyb to do
! horizontal viscosity, thermal conduction, and diffusion of major 
! species (O and O2) with global mean coefficients.
!raa********************************************************************
      SUBROUTINE idea_deldifs(RTE,WE,QME,XE,YE,TEME,
     X                   RTO,WO,QMO,XO,YO,TEMO,DELTIM,
     X                   LS_NODE,COEF00,K_LEVEL,visc,cond,diff)
!
! Apr 06 2012   Henry Juang, copy from deldif and modify for idea
! Jan 16 2013   Jun Wang,    separate init from deldif, and remove module
!                            variables to make this routine thread safe
!
      use gfs_dyn_resol_def
      use gfs_dyn_layout1
      use gfs_dyn_coordinate_def				! hmhj
      use gfs_dyn_deldifs_def
      use gfs_dyn_physcons, rerth => con_rerth
     &                    ,  rd => con_rd, cp => con_cp
      
      IMPLICIT NONE
!
      REAL(KIND=KIND_EVOD) RTE(LEN_TRIE_LS,2,LEVS,NTRAC)
     &,                     WE(LEN_TRIE_LS,2)
     &,                    QME(LEN_TRIE_LS,2)
     &,                     XE(LEN_TRIE_LS,2)
     &,                     YE(LEN_TRIE_LS,2)
     &,                   TEME(LEN_TRIE_LS,2)
     &,                     PE(LEN_TRIE_LS,2)
!
      REAL(KIND=KIND_EVOD) RTO(LEN_TRIO_LS,2,LEVS,NTRAC)
     &,                     WO(LEN_TRIO_LS,2)
     &,                    QMO(LEN_TRIO_LS,2)
     &,                     XO(LEN_TRIO_LS,2)
     &,                     YO(LEN_TRIO_LS,2)
     &,                   TEMO(LEN_TRIO_LS,2)
     &,                     PO(LEN_TRIO_LS,2)
!
      REAL(KIND=KIND_EVOD) DELTIM
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
!raa********************************************************************
! idea change1
! IDEA-related changes
      INTEGER              L,LOCL,N,N0,ND,NP,NPD
!
      INTEGER              INDEV
      INTEGER              INDOD
      integer              indev1,indev2
      integer              indod1,indod2
      real(kind=kind_evod), parameter :: rkappa = cp / rd
      REAL(KIND=KIND_EVOD) DN1,REALVAL,RTNP,DF_DK,FACT
     &,                    SLRD0,FTRD1,RFACT,RFACTRD,RTRD1,FSHK, tem
! INPUT
! - global mean coefficients of viscosity, conductivity, and diffusion
!
      REAL(KIND=KIND_EVOD),intent(in):: visc(levs),cond(levs),diff(levs)
!
! LOCAL VARIABLES
! Physical coefficients DN*, FACT*, RFACT*, etc., are all different for
!     temperature, dynamics, and tracer transport: dn2, dn3, dn4, fact2,
!     fact3, fact4, etc.
!
      REAL(KIND=KIND_EVOD) dn2,dn3,dn4,fact2,fact3,fact4,
     $     rfact2,rfact3,rfact4
!
!C
! idea change1 end
!raa********************************************************************
!
      REAL(KIND=KIND_EVOD), parameter :: CONS0=0.0, CONS1=1.0, CONS2=2.0
!
      INTEGER              INDLSEV,JBASEV
      INTEGER              INDLSOD,JBASOD
!
      INCLUDE 'function2'
!
!     print *,' enter idea_deldifs ' 				! hmhj
!......................................................................
!
craa********************************************************************
craa********************************************************************
!
!......................................................................
!
      CALL COUNTPERF(0,13,0.)
!!
      K=K_LEVEL
!
!!
!raa********************************************************************
! idea change4
!
! Global mean horizontal transport of heat, momentum and tracers is
! incorporated in all layers based on following assumptions
!     1) physical diffusion acts on all scales (N0=0) and
!     2) it is quadratic (JDEL=2)
!     3) diffusion coefficients for all tracers are the same as for
!       pairs O-N2 and O-O2, which are equal (D12=D13), which should be
!       fine for all tracers as of May 2008, because those different
!       from these three are not present where horizontal diffusion is 
!       important
!
! Precalculate arrays of coefficients dneidea dnoide. are set up in 
! in idea_deldif_init, the scalar coefficients are computed every 
! time step
!
! First calculate scalar coefficients (for a given vertical index k)
! dn* for conductivity, viscosity, and diffusion 
!
      dn2=cons2*cond(k)/(rerth**2)
      dn3=cons2*visc(k)/(rerth**2)
      dn4=cons2*diff(k)/(rerth**2)
!
! (I think cons2 is here because the time step is 2.*deltim)
!
!C......................................................................
!C
! The following do-loops originally did only numerical diffusion and
! Rayleigh friction. Now physical diffusion is inserted as well.
!
! -Evens
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
!
! fact* and rfact* are for conductivity, viscosity, and diffusion
!
            FACT             = DELTIM*DNE(INDEV)*RTHK(K)
            fact2=deltim*dneidea(indev)*dn2        
            fact3=deltim*dneidea(indev)*dn3
            fact4=deltim*dneidea(indev)*dn4
!            RFACT            = CONS1/(CONS1+FACT)
!            RFACTRD          = CONS1/(CONS1+FACT+DELTIM*RTRD(K))
            rfact2  =CONS1/(CONS1+FACT+fact2) 
            rfact3  =CONS1/(CONS1+FACT+DELTIM*RTRD(K)+fact3) 
            rfact4  =CONS1/(CONS1+FACT+fact4) 
            
            WE(INDEV,1)      = WE(INDEV,1)*rfact3
            WE(INDEV,2)      = WE(INDEV,2)*rfact3
            XE(INDEV,1)      = XE(INDEV,1)*rfact3
            XE(INDEV,2)      = XE(INDEV,2)*rfact3
!          if(k==71.and.locl==1.and.INDEV==4) 
!     &    print *,'in idea_del,locl=',locl,'l=',l,'indev=',
!     &  indev,indev1,indev2,'dneidea(indev)=',dneidea(indev),'dn3=',dn3,
!     &   dneidea(indev)*dn3,'fact3=',fact3,'rfact3=',rfact3,'fact=',
!     &   fact,'RTRD(K)=',RTRD(K),'we(4,1)=',we(4,1),we(4,2)
            
            RTE(INDEV,1,k,1) = RTE(INDEV,1,k,1)*rfact4
            RTE(INDEV,2,k,1) = RTE(INDEV,2,k,1)*rfact4

            PE(INDEV,1)=BKLY(K)*QME(INDEV,1)+CKLY(K)*TEME(INDEV,1) ! hmhj
            PE(INDEV,2)=BKLY(K)*QME(INDEV,2)+CKLY(K)*TEME(INDEV,2) ! hmhj
            YE(INDEV,1)      =  ( YE(INDEV,1)+
     X           (FACT+fact2)*COEF00(K,1)* PE(INDEV,1) )*rfact2
            YE(INDEV,2)      = ( YE(INDEV,2) +
     X           (FACT+fact2)*COEF00(K,1)* PE(INDEV,2) )*rfact2
 
         ENDDO
      ENDDO
!     print*,'ok5'
!
! -Odds
!
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
            fact2=deltim*dnoidea(indod)*dn2
            fact3=deltim*dnoidea(indod)*dn3
            fact4=deltim*dnoidea(indod)*dn4
            rfact2  =CONS1/(CONS1+FACT+fact2)
            rfact3  =CONS1/(CONS1+FACT+DELTIM*RTRD(K)+fact3)
            rfact4  =CONS1/(CONS1+FACT+fact4)
            
            WO(INDOD,1)      = WO(INDOD,1)*rfact3
            WO(INDOD,2)      = WO(INDOD,2)*rfact3
            XO(INDOD,1)      = XO(INDOD,1)*rfact3
            XO(INDOD,2)      = XO(INDOD,2)*rfact3
            
            RTO(INDOD,1,k,1) = RTO(INDOD,1,k,1)*rfact4
            RTO(INDOD,2,k,1) = RTO(INDOD,2,k,1)*rfact4
            
            PO(INDOD,1)=BKLY(K)*QMO(INDOD,1)+CKLY(K)*TEMO(INDOD,1) ! hmhj
            PO(INDOD,2)=BKLY(K)*QMO(INDOD,2)+CKLY(K)*TEMO(INDOD,2) ! hmhj
            YO(INDOD,1)      = ( YO(INDOD,1)+
     X           (FACT+fact2)*COEF00(K,1)* PO(INDOD,1) )*rfact2
            YO(INDOD,2)      = ( YO(INDOD,2)+
     X           (FACT+fact2)*COEF00(K,1)* PO(INDOD,2) )*rfact2
            
         ENDDO
      ENDDO
!     print*,'ok6'
!
!......................................................................
!
! Now do tracers 2-ntrac. Apply the same diffusion coefficient to all.
! This version also does correction for these tracers similar to the one
! done for enthalpy above. The correction may be gone in the future.
!
      if (ntrac .gt. 1) then
         do it=2,ntrac
!
! -Evens
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
               DO INDEV = indev1 , indev2
 
                  FACT              = DELTIM*DNE(INDEV)*RTHK(K)
                  fact4=deltim*dneidea(indev)*dn4
                  rfact4  =CONS1/(CONS1+FACT+fact4)
 
                  PE(INDEV,1)=BKLY(K)*QME(INDEV,1)+CKLY(K)*TEME(INDEV,1) ! hmhj
                  PE(INDEV,2)=BKLY(K)*QME(INDEV,2)+CKLY(K)*TEME(INDEV,2) ! hmhj

                  RTE(INDEV,1,k,IT) = ( RTE(INDEV,1,k,IT)+
     X                 (FACT+fact4)*COEF00(K,it)*PE(INDEV,1))*rfact4
                  RTE(INDEV,2,k,IT) = ( RTE(INDEV,2,k,IT)+
     X                 (FACT+fact4)*COEF00(K,it)* PE(INDEV,2))*rfact4
 
               ENDDO
            ENDDO
!     print*,'ok7'
!
!......................................................................
!
!
! -Odds
!
            DO LOCL=1,LS_MAX_NODE
                  L=LS_NODE(LOCL,1)
               JBASOD=LS_NODE(LOCL,3)
               indod1 = indlsod(L+1,L)
               if (mod(L,2).eq.mod(jcap+1,2)) then
                  indod2 = indlsod(jcap  ,L)
               else
                  indod2 = indlsod(jcap+1,L)
               endif

               DO INDOD = indod1 , indod2
 
                  FACT              = DELTIM*DNO(INDOD)*RTHK(K)
                  fact4=deltim*dnoidea(indod)*dn4
                  rfact4  =CONS1/(CONS1+FACT+fact4)
 
                  PO(INDOD,1)=BKLY(K)*QMO(INDOD,1)+CKLY(K)*TEMO(INDOD,1) ! hmhj
                  PO(INDOD,2)=BKLY(K)*QMO(INDOD,2)+CKLY(K)*TEMO(INDOD,2) ! hmhj

                  RTO(INDOD,1,k,IT) = ( RTO(INDOD,1,k,IT)+
     X                 (FACT+fact4)*COEF00(K,it)* PO(INDOD,1) )*rfact4
                  RTO(INDOD,2,k,IT) = ( RTO(INDOD,2,k,IT)+
     X                 (FACT+fact4)*COEF00(K,it)* PO(INDOD,2) )*rfact4
 
               ENDDO
            ENDDO
!     print*,'ok8'
         enddo                  ! Tracer do loop end
      endif                     ! If ntrac > 1
! idea change4 end
!raa********************************************************************
!
      CALL COUNTPERF(1,13,0.)

!     print *,' leave idea_deldifs ' 	! hmhj
!!
      RETURN
      END
