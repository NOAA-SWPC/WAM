!
! !MODULE: Nst_Var_ESMFMod  ---                Definition of the Nst_Var model
!                                           fields in the ESMF internal state.
!
! !DESCRIPTION: Nst_Var_ESMFMod ---            Define the Nst_Var model  variables
!                                            in the ESMF internal state.
!---------------------------------------------------------------------------
! !REVISION HISTORY:
!
!  May 2008      Shrinivas Moorthi Initial code.
!  Aug 2009      Xu Li for DTM-1p
!
! !INTERFACE:
!
 MODULE gfs_physics_nst_var_mod

 use machine , only : kind_phys

 IMPLICIT none

 TYPE Nst_Var_Data
    real(kind_phys),pointer:: slmsk    (:,:)=>null()
    real(kind_phys),pointer:: xt       (:,:)=>null()
    real(kind_phys),pointer:: xs       (:,:)=>null()
    real(kind_phys),pointer:: xu       (:,:)=>null()
    real(kind_phys),pointer:: xv       (:,:)=>null()
    real(kind_phys),pointer:: xz       (:,:)=>null()
    real(kind_phys),pointer:: zm       (:,:)=>null()
    real(kind_phys),pointer:: xtts     (:,:)=>null()
    real(kind_phys),pointer:: xzts     (:,:)=>null()
    real(kind_phys),pointer:: dt_cool  (:,:)=>null()
    real(kind_phys),pointer:: z_c      (:,:)=>null()
    real(kind_phys),pointer:: c_0      (:,:)=>null()
    real(kind_phys),pointer:: c_d      (:,:)=>null()
    real(kind_phys),pointer:: w_0      (:,:)=>null()
    real(kind_phys),pointer:: w_d      (:,:)=>null()
    real(kind_phys),pointer:: d_conv   (:,:)=>null()
    real(kind_phys),pointer:: ifd      (:,:)=>null()
    real(kind_phys),pointer:: tref     (:,:)=>null()
    real(kind_phys),pointer:: Qrain    (:,:)=>null()
 end type Nst_Var_Data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    contains
    subroutine nstvar_aldata(dim1,dim2,data,iret)
       implicit none
       integer, intent(in)               :: dim1, dim2
       type(nst_var_data),intent(inout)  :: data
       integer, intent(out)              :: iret
!
allocate(                         &
      data%slmsk   (dim1,dim2),   &
      data%xt      (dim1,dim2),   &
      data%xs      (dim1,dim2),   &
      data%xu      (dim1,dim2),   &
      data%xv      (dim1,dim2),   &
      data%xz      (dim1,dim2),   &
      data%zm      (dim1,dim2),   &
      data%xtts    (dim1,dim2),   &
      data%xzts    (dim1,dim2),   &
      data%dt_cool (dim1,dim2),   &
      data%z_c     (dim1,dim2),   &
      data%c_0     (dim1,dim2),   &
      data%c_d     (dim1,dim2),   &
      data%w_0     (dim1,dim2),   &
      data%w_d     (dim1,dim2),   &
      data%d_conv  (dim1,dim2),   &
      data%ifd     (dim1,dim2),   &
      data%tref    (dim1,dim2),   &
      data%Qrain   (dim1,dim2),   &
      stat=iret)
    if(iret.ne.0) iret=-3
    return
  end subroutine nstvar_aldata
 END MODULE gfs_physics_nst_var_mod
