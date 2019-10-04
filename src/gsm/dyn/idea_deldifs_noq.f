craa********************************************************************
c Starting May 27, 2008: modified by RAA to impose horizontal physical
C     dissipation on top of numerical diffusion in all layers.
c Jan 4, 2007: modified by Rashid Akmaev based on DELDIFS_hyb to do
c horizontal viscosity, thermal conduction, and diffusion of major 
c species (O and O2) with global mean coefficients.
c add without traver
craa********************************************************************
      SUBROUTINE idea_deldifs_noq(WE,QME,XE,YE,TEME,
     X                   WO,QMO,XO,YO,TEMO,DELTIM,
     X                   LS_NODE,COEF00,K_LEVEL,visc,cond,diff)
!
! Apr 06 2012    Henry Juang, copy from deldif, modify for idea and ndsl
! Jan 16 2013    Jun Wang,    separate init from this file, remove module
!                             variables to make this routine thread safe
!
      use gfs_dyn_resol_def
      use gfs_dyn_layout1
      use gfs_dyn_coordinate_def				! hmhj
      use gfs_dyn_deldifs_def
      use gfs_dyn_physcons, rerth => con_rerth
     &                    ,  rd => con_rd, cp => con_cp
      
      IMPLICIT NONE
!
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
craa********************************************************************
c idea change1
c IDEA-related changes
      INTEGER              L,LOCL,N,N0,ND,NP,NPD
!
      INTEGER              INDEV
      INTEGER              INDOD
      integer              indev1,indev2
      integer              indod1,indod2
      real(kind=kind_evod), parameter :: rkappa = cp / rd
      REAL(KIND=KIND_EVOD) DN1,REALVAL,RTNP,DF_DK,FACT
     &,                    SLRD0,FTRD1,RFACT,RFACTRD,RTRD1,FSHK, tem
c INPUT
c - global mean coefficients of viscosity, conductivity, and diffusion
c
      REAL(KIND=KIND_EVOD),intent(in):: visc(levs),cond(levs),diff(levs)
c
c LOCAL VARIABLES
c Physical coefficients DN*, FACT*, RFACT*, etc., are all different for
C     temperature, dynamics, and tracer transport: dn2, dn3, dn4, fact2,
C     fact3, fact4, etc.
c
      REAL(KIND=KIND_EVOD) dn2,dn3,dn4,fact2,fact3,fact4,
     $     rfact2,rfact3,rfact4
CC
c idea change1 end
craa********************************************************************
!
      REAL(KIND=KIND_EVOD), parameter :: CONS0=0.0, CONS1=1.0, CONS2=2.0
!
      INTEGER              INDLSEV,JBASEV
      INTEGER              INDLSOD,JBASOD
!
      INCLUDE 'function2'
!
!     print *,' enter idea_deldifs_noq ' 				! hmhj
!......................................................................
!
craa********************************************************************
!
!......................................................................
!
      CALL COUNTPERF(0,13,0.)
!!
      K=K_LEVEL
!
!!
craa********************************************************************
c idea change4
c
c Global mean horizontal transport of heat, momentum and tracers is
c incorporated in all layers based on following assumptions
c     1) physical diffusion acts on all scales (N0=0) and
c     2) it is quadratic (JDEL=2)
c     3) diffusion coefficients for all tracers are the same as for
c       pairs O-N2 and O-O2, which are equal (D12=D13), which should be
c       fine for all tracers as of May 2008, because those different
c       from these three are not present where horizontal diffusion is 
c       important
c
c Precalculate arrays of coefficients dnce, dnve, dnde, etc., analogous
c to the coefficients dne and dno for numerical diffusion. Has to be 
c done every timestep. The only reason these coefficients are 
c precalculated is to keep the original structure of this subroutine.
c Otherwise all the relevant calculations could have been done 
c compactly and probably more efficiently in one loop.
c
c First calculate scalar coefficients (for a given vertical index k)
c dn* for conductivity, viscosity, and diffusion 
c
      dn2=cons2*cond(k)/(rerth**2)
      dn3=cons2*visc(k)/(rerth**2)
      dn4=cons2*diff(k)/(rerth**2)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
CC
CC......................................................................
CC
c The following do-loops originally did only numerical diffusion and
c Rayleigh friction. Now physical diffusion is inserted as well.
c
c -Evens
c
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
c
c fact* and rfact* are for conductivity, viscosity, and diffusion
c
            FACT             = DELTIM*DNE(INDEV)*RTHK(K)
            fact2=deltim*dneidea(indev)*dn2        
            fact3=deltim*dneidea(indev)*dn3
            fact4=deltim*dneidea(indev)*dn4
c            RFACT            = CONS1/(CONS1+FACT)
c            RFACTRD          = CONS1/(CONS1+FACT+DELTIM*RTRD(K))
            rfact2  =CONS1/(CONS1+FACT+fact2) 
            rfact3  =CONS1/(CONS1+FACT+DELTIM*RTRD(K)+fact3) 
            rfact4  =CONS1/(CONS1+FACT+fact4) 
            
c            WE(INDEV,1)      = WE(INDEV,1)*RFACTRD
c            WE(INDEV,2)      = WE(INDEV,2)*RFACTRD
c            XE(INDEV,1)      = XE(INDEV,1)*RFACTRD
c            XE(INDEV,2)      = XE(INDEV,2)*RFACTRD
            WE(INDEV,1)      = WE(INDEV,1)*rfact3
            WE(INDEV,2)      = WE(INDEV,2)*rfact3
            XE(INDEV,1)      = XE(INDEV,1)*rfact3
            XE(INDEV,2)      = XE(INDEV,2)*rfact3
            
            PE(INDEV,1)=BKLY(K)*QME(INDEV,1)+CKLY(K)*TEME(INDEV,1) ! hmhj
            PE(INDEV,2)=BKLY(K)*QME(INDEV,2)+CKLY(K)*TEME(INDEV,2) ! hmhj
c            YE(INDEV,1)      =  ( YE(INDEV,1)+
c     X           FACT*COEF00(K,1)* PE(INDEV,1) )*RFACT ! hmhj
c            YE(INDEV,2)      = ( YE(INDEV,2) +
c     X           FACT*COEF00(K,1)* PE(INDEV,2) )*RFACT ! hmhj
            YE(INDEV,1)      =  ( YE(INDEV,1)+
     X           (FACT+fact2)*COEF00(K,1)* PE(INDEV,1) )*rfact2
            YE(INDEV,2)      = ( YE(INDEV,2) +
     X           (FACT+fact2)*COEF00(K,1)* PE(INDEV,2) )*rfact2
 
         ENDDO
      ENDDO
c     print*,'ok5'
c
c -Odds
c
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
c            RFACT            = CONS1/(CONS1+FACT)
c            RFACTRD          = CONS1/(CONS1+FACT+DELTIM*RTRD(K))
            rfact2  =CONS1/(CONS1+FACT+fact2)
            rfact3  =CONS1/(CONS1+FACT+DELTIM*RTRD(K)+fact3)
            rfact4  =CONS1/(CONS1+FACT+fact4)
            
c            WO(INDOD,1)      = WO(INDOD,1)*RFACTRD
c            WO(INDOD,2)      = WO(INDOD,2)*RFACTRD
c            XO(INDOD,1)      = XO(INDOD,1)*RFACTRD
c            XO(INDOD,2)      = XO(INDOD,2)*RFACTRD
            WO(INDOD,1)      = WO(INDOD,1)*rfact3
            WO(INDOD,2)      = WO(INDOD,2)*rfact3
            XO(INDOD,1)      = XO(INDOD,1)*rfact3
            XO(INDOD,2)      = XO(INDOD,2)*rfact3
            
            PO(INDOD,1)=BKLY(K)*QMO(INDOD,1)+CKLY(K)*TEMO(INDOD,1) ! hmhj
            PO(INDOD,2)=BKLY(K)*QMO(INDOD,2)+CKLY(K)*TEMO(INDOD,2) ! hmhj
c            YO(INDOD,1)      = ( YO(INDOD,1)+
c     X           FACT*COEF00(K,1)* PO(INDOD,1) )*RFACT ! hmhj
c            YO(INDOD,2)      = ( YO(INDOD,2)+
c     X           FACT*COEF00(K,1)* PO(INDOD,2) )*RFACT ! hmhj
            YO(INDOD,1)      = ( YO(INDOD,1)+
     X           (FACT+fact2)*COEF00(K,1)* PO(INDOD,1) )*rfact2
            YO(INDOD,2)      = ( YO(INDOD,2)+
     X           (FACT+fact2)*COEF00(K,1)* PO(INDOD,2) )*rfact2
            
         ENDDO
      ENDDO
!
!......................................................................
!
      CALL COUNTPERF(1,13,0.)

!     print *,' leave idea_deldifs_noq ' 	! hmhj
!!
      RETURN
      END
