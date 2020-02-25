!!------------------------------------------------
! !MODULE:   wam_diag3d_mod: WAM_physics 
!            Define the 3d WAM Diagnostics 
!                                     
!
! 
!
! Apr 2017 VAY: For diagnostics of GW physics nc-files
!------------------------------------------------
      MODULE wam_diag3d_mod

      use machine , only : kind_phys
      IMPLICIT none

      TYPE G3D_WAMD
!----------------------------------------------------
! 3D diag fields
      real(kind=kind_phys),  pointer :: daxz(:,:,:)
      real(kind=kind_phys),  pointer :: dayz(:,:,:)
      real(kind=kind_phys),  pointer :: deps(:,:,:)
      real(kind=kind_phys),  pointer :: dked(:,:,:)
!add 7-arrays
      real(kind=kind_phys),  pointer :: zgkm(:,:,:)
!
      real(kind=kind_phys),  pointer :: hrad(:,:,:)
      real(kind=kind_phys),  pointer :: crad(:,:,:)
      real(kind=kind_phys),  pointer :: joul(:,:,:)
      real(kind=kind_phys),  pointer :: diax(:,:,:)
      real(kind=kind_phys),  pointer :: diay(:,:,:)
      real(kind=kind_phys),  pointer :: kzzt(:,:,:)
!----------------------------------------------------
!  oro-waves 
!    real(kind=kind_phys),  pointer :: maxz(:,:,:)
!    real(kind=kind_phys),  pointer :: mayz(:,:,:)
!    real(kind=kind_phys),  pointer :: meps(:,:,:)
!    real(kind=kind_phys),  pointer :: mked(:,:,:)  
! convective waves 
!    real(kind=kind_phys),  pointer :: caxz(:,:,:)
!    real(kind=kind_phys),  pointer :: cayz(:,:,:)
!    real(kind=kind_phys),  pointer :: ceps(:,:,:)
!    real(kind=kind_phys),  pointer :: cked(:,:,:)
      end type G3D_WAMD


      contains

! subROUTINE: g3dwamd_alloc ---

!---------------------------------------------------------------------------
      subroutine g3dwamd_alloc(dim1, dim2, dim3, g3d_wamfld, iret)

      implicit none
      integer, intent(in)                :: dim1, dim2, dim3
      TYPE(G3D_WAMD), INTENT(out)    :: g3d_wamfld
      integer, intent(out)               :: iret
!

      allocate(                                    
     &      g3d_wamfld%daxz    (dim1,dim2,dim3), 
     &      g3d_wamfld%dayz    (dim1,dim2,dim3), 
     &      g3d_wamfld%deps (dim1,dim2,dim3), 
     &      g3d_wamfld%dked (dim1,dim2,dim3), 
     &                          stat = iret)
      if(iret.ne.0) iret=-3
 
      g3d_wamfld%daxz    (1:dim1, 1:dim2, 1:dim3) = 0.
      g3d_wamfld%dayz    (1:dim1, 1:dim2, 1:dim3) = 0.
      g3d_wamfld%deps    (1:dim1, 1:dim2, 1:dim3) = 0.
      g3d_wamfld%dked    (1:dim1, 1:dim2, 1:dim3) = 0.
!
      allocate(                                    
     &      g3d_wamfld%diax    (dim1,dim2,dim3), 
     &      g3d_wamfld%diay    (dim1,dim2,dim3), 
     &      g3d_wamfld%joul    (dim1,dim2,dim3), 
     &      g3d_wamfld%hrad    (dim1,dim2,dim3), 
     &      g3d_wamfld%crad    (dim1,dim2,dim3), 
     &      g3d_wamfld%kzzt    (dim1,dim2,dim3), 
     &      g3d_wamfld%zgkm (dim1,dim2,dim3), 
     &                          stat = iret)
      g3d_wamfld%diax    (1:dim1, 1:dim2, 1:dim3) = 0.
      g3d_wamfld%diay    (1:dim1, 1:dim2, 1:dim3) = 0.
      g3d_wamfld%joul    (1:dim1, 1:dim2, 1:dim3) = 0.
      g3d_wamfld%hrad    (1:dim1, 1:dim2, 1:dim3) = 0.
      g3d_wamfld%crad    (1:dim1, 1:dim2, 1:dim3) = 0.
      g3d_wamfld%kzzt    (1:dim1, 1:dim2, 1:dim3) = 0.
      g3d_wamfld%zgkm    (1:dim1, 1:dim2, 1:dim3) = 0.
      return
      end subroutine g3dwamd_alloc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      END MODULE  wam_diag3d_mod
!
!      MODULE wam_pass_diag_types
!      use wam_diag3d_mod, only :  G3D_WAMD
!      TYPE(G3D_WAMD)   :: GIS_WAM
!      END MODULE wam_pass_diag_types
