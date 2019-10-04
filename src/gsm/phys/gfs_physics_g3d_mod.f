!
! !MODULE: gfs_physics_g3d_mod  ---      Definition of the 3d atmospheric Diag 
!                                        fields in the ESMF internal state.
!
! !DESCRIPTION: gfs_physics_g3d_mod ---    Define the 3d atmospheric diag 
!                                          states in the ESMF internal state.
!---------------------------------------------------------------------------
! !REVISION HISTORY:
!
!  2009/12/08      Sarah Lu, Initial code.
!  2009/12/14      Sarah Lu, dqdt type changed from kind_rad to kind_phys;
!  2011/12/14      Sarah Lu, initialize fcld and dqdt
!  2014/12/18      Jun Wang, add cnv fields cnv_mfc, cnv_mfd, cnv_qc
!
! !INTERFACE:
!
 MODULE gfs_physics_g3d_mod

 use machine , only : kind_rad, kind_phys
 IMPLICIT none

 TYPE G3D_Var_Data

! 3D diag fields
    real(kind=kind_rad),   pointer :: fcld(:,:,:)
    real(kind=kind_phys),  pointer :: dqdt(:,:,:)
    real(kind=kind_phys),  pointer :: cnv_mfc(:,:,:)
    real(kind=kind_phys),  pointer :: cnv_mfd(:,:,:)
    real(kind=kind_phys),  pointer :: cnv_qc(:,:,:)
   
 end type G3D_Var_Data

 contains

! !IROUTINE: g3d_aldata ---

!---------------------------------------------------------------------------
    subroutine g3d_aldata(dim1, dim2, dim3, g3d_fld, iret)

    implicit none
    integer, intent(in)                :: dim1, dim2, dim3
    TYPE(G3D_Var_Data), INTENT(out)    :: g3d_fld
    integer, intent(out)               :: iret
!

allocate(                                    &
           g3d_fld%fcld    (dim1,dim2,dim3), &
           g3d_fld%dqdt    (dim1,dim2,dim3), &
           g3d_fld%cnv_mfc (dim1,dim2,dim3), &
           g3d_fld%cnv_mfd (dim1,dim2,dim3), &
           g3d_fld%cnv_qc  (dim1,dim2,dim3), &
                               stat = iret)
    if(iret.ne.0) iret=-3

    g3d_fld%fcld     (1:dim1, 1:dim2, 1:dim3) = 0.
    g3d_fld%dqdt     (1:dim1, 1:dim2, 1:dim3) = 0.
    g3d_fld%cnv_mfc  (1:dim1, 1:dim2, 1:dim3) = 0.
    g3d_fld%cnv_mfd  (1:dim1, 1:dim2, 1:dim3) = 0.
    g3d_fld%cnv_qc   (1:dim1, 1:dim2, 1:dim3) = 0.


    return
  end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 END MODULE gfs_physics_g3d_mod
