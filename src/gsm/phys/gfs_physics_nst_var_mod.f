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
!  Jun 2016      Fanglin Yang remove pointer for digital filter
!  Jun 2016      Xu Li add nst_init
!
! !INTERFACE:
!
 MODULE gfs_physics_nst_var_mod

 use machine , only : kind_phys

 IMPLICIT none

 TYPE Nst_Var_Data
    real(kind_phys),allocatable, dimension(:,:) ::                &
           slmsk,   xt,   xs,  xu,  xv,  xz,  zm,     xtts, xzts, &
           dt_cool, z_c,  c_0, c_d, w_0, w_d, d_conv, ifd,        &
           tref,    Qrain
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

    subroutine nst_init(nst_fld,iret)
       implicit none
       type(nst_var_data),intent(inout)  :: nst_fld
       integer, intent(out)              :: iret
!
      nst_fld%slmsk   = 0.0 
      nst_fld%xt      = 0.0      
      nst_fld%xs      = 0.0      
      nst_fld%xu      = 0.0      
      nst_fld%xv      = 0.0      
      nst_fld%xz      = 0.0      
      nst_fld%zm      = 0.0      
      nst_fld%xtts    = 0.0    
      nst_fld%xzts    = 0.0    
      nst_fld%dt_cool = 0.0 
      nst_fld%z_c     = 0.0     
      nst_fld%c_0     = 0.0     
      nst_fld%c_d     = 0.0     
      nst_fld%w_0     = 0.0     
      nst_fld%w_d     = 0.0     
      nst_fld%d_conv  = 0.0  
      nst_fld%ifd     = 0.0     
      nst_fld%tref    = 0.0   
      nst_fld%Qrain   = 0.0   
      if(iret.ne.0) iret=-3
      return
    end subroutine nst_init
    subroutine nstvar_axdata(data)
       implicit none
       type(nst_var_data),intent(inout)  :: data

       deallocate(                                                              &
           data%slmsk,  data%xt,   data%xs,   data%xu,      data%xv,   data%xz, &
           data%zm,     data%xtts, data%xzts, data%dt_cool, data%z_c,  data%c_0,&
           data%c_d,    data%w_0,  data%w_d,  data%d_conv,  data%ifd,           &
           data%tref,   data%Qrain )

    return
    end subroutine nstvar_axdata
 END MODULE gfs_physics_nst_var_mod
