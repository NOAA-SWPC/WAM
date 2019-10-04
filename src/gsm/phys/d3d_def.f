      module d3d_def
      use machine
      implicit none
!
      REAL(KIND=KIND_RAD) ,ALLOCATABLE :: DT3DT(:,:,:,:,:)
      REAL(KIND=KIND_RAD) ,ALLOCATABLE :: DQ3DT(:,:,:,:,:)
      REAL(KIND=KIND_RAD) ,ALLOCATABLE :: DU3DT(:,:,:,:,:)
      REAL(KIND=KIND_RAD) ,ALLOCATABLE :: DV3DT(:,:,:,:,:)
      REAL(KIND=KIND_RAD) ,ALLOCATABLE :: upd_mf(:,:,:,:)
     &,                                   dwn_mf(:,:,:,:)
     &,                                   det_mf(:,:,:,:)
!    &,                                   dkh(:,:,:,:)
!    &,                                   rnp(:,:,:,:)
      REAL(KIND=KIND_RAD) ,ALLOCATABLE :: CLDCOV(:,:,:)
!
      contains
!
      subroutine d3d_init(ngptc,nblck,lonr,lats_node_r,levs,pl_coeff)
      implicit none
      integer ngptc,nblck,lonr,lats_node_r,levs,pl_coeff
!
        allocate (DT3DT(NGPTC,LEVS,6,NBLCK,lats_node_r))
        allocate (DU3DT(NGPTC,LEVS,4,NBLCK,lats_node_r))
        allocate (DV3DT(NGPTC,LEVS,4,NBLCK,lats_node_r))
        allocate (DQ3DT(NGPTC,LEVS,5+pl_coeff,NBLCK,lats_node_r))
        allocate (CLDCOV(levs,LONR,lats_node_r))
        allocate (upd_mf(NGPTC,LEVS,NBLCK,lats_node_r))
        allocate (dwn_mf(NGPTC,LEVS,NBLCK,lats_node_r))
        allocate (det_mf(NGPTC,LEVS,NBLCK,lats_node_r))
!
      end subroutine d3d_init
!
      subroutine d3d_zero
      implicit none
      real, parameter :: zero=0.0
!
         DT3DT  = zero
         DU3DT  = zero
         DV3DT  = zero
         DQ3DT  = zero
         CLDCOV = zero
         upd_mf = zero
         dwn_mf = zero
         det_mf = zero
!
      end subroutine d3d_zero
      end module d3d_def
