      SUBROUTINE omega_sig(NJEFF,NSIZE_AR,RCL,
     1                     PS,DPHI,DLAM,DG,UG,VG,VVEL)

C....   CODE LIFTED FROM POST (MCP1840) JUN 88--COMPUTES VVEL (CB/SEC)
C....    INPUT PS IN CB,OUTPUT VVEL IN CB/SEC
C....   DO LOOPS ALTERED FOR BETTER VECTORIZATION POSSIBILITIES..K.A.C.
cc
      USE gfs_dyn_MACHINE , ONLY : kind_grid
      use gfs_dyn_resol_def
      use gfs_dyn_vert_def
      implicit none

!!
      integer              i,k,le,NSIZE_AR,NJEFF
      real(kind=kind_grid) rcl
cc
      real(kind=kind_grid) DPHI(NSIZE_AR),DLAM(NSIZE_AR),PS(NSIZE_AR)
      real(kind=kind_grid) UG(NSIZE_AR,levs),VG(NSIZE_AR,levs)
      real(kind=kind_grid) DG(NSIZE_AR,levs)
      real(kind=kind_grid) VVEL(NSIZE_AR,levs)
C...   VVEL CONTAINS OMEGA IN LAYERS ON RETURN FROM SUBROUTINE...
      real(kind=kind_grid) DB(NJEFF,levs),CB(NJEFF,levs),
     &                     DOT(NJEFF,levs+1),CG(NJEFF,levs)
!!
      DO 1 K=1,levs+1
        DO 49 i=1,NJEFF
          DOT(i,K) = 0.E0
   49 CONTINUE
    1 CONTINUE
C...  COMPUTE C=V(TRUE)*DEL(LN(PS)).DIVIDE BY COS FOR DEL COS FOR V
      DO 48 i=1,NJEFF
        DPHI(i)=DPHI(i)*RCL
        DLAM(i)=DLAM(i)*RCL
   48 CONTINUE
      DO 5 LE=1,levs
        DO 50 i=1,NJEFF
          CG(i,LE)=UG(i,LE)*DLAM(i)+VG(i,LE)*DPHI(i)
   50 CONTINUE
    5 CONTINUE
      DO 51 i=1,NJEFF
        DB(i,1)=DEL(1)*DG(i,1)
        CB(i,1)=DEL(1)*CG(i,1)
   51 CONTINUE
!!
      DO 6 LE=1,levs-1
        DO 52 i=1,NJEFF
          DB(i,LE+1)=DB(i,LE)+DEL(LE+1)*DG(i,LE+1)
          CB(i,LE+1)=CB(i,LE)+DEL(LE+1)*CG(i,LE+1)
   52 CONTINUE
    6 CONTINUE
!!
C...    SIGMA DOT COMPUTED ONLY AT INTERIOR INTERFACES
      DO 7 K=1,levs-1
        DO 53 i=1,NJEFF
          DOT(i,K+1)=DOT(i,K)+DEL(K)
     1               *(DB(i,levs)+CB(i,levs)-DG(i,K)-CG(i,K))
   53 CONTINUE
    7 CONTINUE
!!
      DO 8 K=1,levs
        DO 54 i=1,NJEFF
          VVEL(i,K)=  SL(K)*(CG(i,K)-CB(i,levs)-DB(i,levs))-
     1                0.5*(DOT(i,K+1)+DOT(i,K))
          VVEL(i,K)=VVEL(i,K)*PS(i)
   54 CONTINUE
    8 CONTINUE
!!
      RETURN
      END
