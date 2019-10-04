!
! !MODULE: gfs_physics_sfc_flx_mod  ---      Definition of the surface
!                                            fields in the ESMF internal state.
!
! !DESCRIPTION: gfs_physics_sfc_flx_mod ---    Define the surfacee  variables
!                                              in the ESMF internal state.
!---------------------------------------------------------------------------
! !REVISION HISTORY:
!
!  March 2007      Shrinivas Moorthi Initial code.
!  February 2009   Tom Henderson Adapted for FIM from nems r3038.
!
! !INTERFACE:
!
 MODULE gfs_physics_sfc_flx_set_mod
!SMS$IGNORE BEGIN

 use infnan , only: inf

 IMPLICIT none

    contains
    subroutine sfcvar_aldata(dim1s, dim1e, dim2, dim3, sfc_fld, iret)

    USE gfs_physics_sfc_flx_mod, ONLY: Sfc_Var_Data

    implicit none
    TYPE(Sfc_Var_Data), INTENT(inout) :: sfc_fld 
    integer, intent(in)               :: dim1s, dim1e, dim2, dim3

    integer, intent(out)             :: iret
!
allocate(                                  &
           sfc_fld%tsea   (dim1s:dim1e,dim2),     &
           sfc_fld%smc    (dim3,dim1s:dim1e,dim2),&
           sfc_fld%sheleg (dim1s:dim1e,dim2),     &
!TBH:  not used by FIM yet
! jbao new gfs phys
           sfc_fld%sncovr (dim1s:dim1e,dim2),     &
           sfc_fld%stc    (dim3,dim1s:dim1e,dim2),&
           sfc_fld%tg3    (dim1s:dim1e,dim2),     &
           sfc_fld%zorl   (dim1s:dim1e,dim2),     &
           sfc_fld%cv     (dim1s:dim1e,dim2),     &
           sfc_fld%cvb    (dim1s:dim1e,dim2),     &
           sfc_fld%cvt    (dim1s:dim1e,dim2),     &
           sfc_fld%alvsf  (dim1s:dim1e,dim2),     &
           sfc_fld%alvwf  (dim1s:dim1e,dim2),     &
           sfc_fld%alnsf  (dim1s:dim1e,dim2),     &
           sfc_fld%alnwf  (dim1s:dim1e,dim2),     &
           sfc_fld%slmsk  (dim1s:dim1e,dim2),     &
           sfc_fld%vfrac  (dim1s:dim1e,dim2),     &
           sfc_fld%canopy (dim1s:dim1e,dim2),     &
           sfc_fld%f10m   (dim1s:dim1e,dim2),     &
           sfc_fld%t2m    (dim1s:dim1e,dim2),     &
           sfc_fld%q2m    (dim1s:dim1e,dim2),     &
           sfc_fld%vtype  (dim1s:dim1e,dim2),     &
           sfc_fld%stype  (dim1s:dim1e,dim2),     &
           sfc_fld%facsf  (dim1s:dim1e,dim2),     &
           sfc_fld%facwf  (dim1s:dim1e,dim2),     &
           sfc_fld%uustar (dim1s:dim1e,dim2),     &
           sfc_fld%ffmm   (dim1s:dim1e,dim2),     &
           sfc_fld%ffhh   (dim1s:dim1e,dim2),     &
           sfc_fld%hice   (dim1s:dim1e,dim2),     &
           sfc_fld%fice   (dim1s:dim1e,dim2),     &
!TBH:  not used by FIM yet
! jbao new gfs phys
           sfc_fld%tisfc  (dim1s:dim1e,dim2),     &
           sfc_fld%tprcp  (dim1s:dim1e,dim2),     &
           sfc_fld%srflag (dim1s:dim1e,dim2),     &
           sfc_fld%snwdph (dim1s:dim1e,dim2),     &
           sfc_fld%slc    (dim3,dim1s:dim1e,dim2),&
           sfc_fld%shdmin (dim1s:dim1e,dim2),     &
           sfc_fld%shdmax (dim1s:dim1e,dim2),     &
           sfc_fld%slope  (dim1s:dim1e,dim2),     &
           sfc_fld%snoalb (dim1s:dim1e,dim2),     &
!TBH:  not used by FIM yet
!          sfc_fld%oro    (dim1s:dim1e,dim2),     &
           stat=iret)
    if(iret.ne.0) iret=-3
#ifndef LAHEY
    sfc_fld%tsea(:,:) = inf
    sfc_fld%smc(:,:,:) = inf
    sfc_fld%sheleg(:,:) = inf
    sfc_fld%sncovr(:,:) = inf
    sfc_fld%stc(:,:,:) = inf
    sfc_fld%tg3(:,:) = inf
    sfc_fld%zorl(:,:) = inf
    sfc_fld%cv(:,:) = inf
    sfc_fld%cvb(:,:) = inf
    sfc_fld%cvt(:,:) = inf
    sfc_fld%alvsf(:,:) = inf
    sfc_fld%alvwf(:,:) = inf
    sfc_fld%alnsf(:,:) = inf
    sfc_fld%alnwf(:,:) = inf
    sfc_fld%slmsk(:,:) = inf
    sfc_fld%vfrac(:,:) = inf
    sfc_fld%canopy(:,:) = inf
    sfc_fld%f10m(:,:) = inf
    sfc_fld%t2m(:,:) = inf
    sfc_fld%q2m(:,:) = inf
    sfc_fld%vtype(:,:) = inf
    sfc_fld%stype(:,:) = inf
    sfc_fld%facsf(:,:) = inf
    sfc_fld%facwf(:,:) = inf
    sfc_fld%uustar(:,:) = inf
    sfc_fld%ffmm(:,:) = inf
    sfc_fld%ffhh(:,:) = inf
    sfc_fld%hice(:,:) = inf
    sfc_fld%fice(:,:) = inf
#endif

    return
  end subroutine
    subroutine flxvar_aldata(dim1s, dim1e, dim2, flx_fld, iret)

    USE gfs_physics_sfc_flx_mod, ONLY: Flx_Var_Data
    implicit none
    TYPE(Flx_Var_Data), INTENT(inout) :: flx_fld 
    integer, intent(in)               :: dim1s, dim1e, dim2

    integer, intent(out)             :: iret
!
    allocate(                          &
          flx_fld%SFCDSW  (dim1s:dim1e,dim2), &
          flx_fld%COSZEN  (dim1s:dim1e,dim2), &
          flx_fld%TMPMIN  (dim1s:dim1e,dim2), &
          flx_fld%TMPMAX  (dim1s:dim1e,dim2), &
          flx_fld%DUSFC   (dim1s:dim1e,dim2), &
          flx_fld%DVSFC   (dim1s:dim1e,dim2), &
          flx_fld%DTSFC   (dim1s:dim1e,dim2), &
          flx_fld%DQSFC   (dim1s:dim1e,dim2), &
          flx_fld%DLWSFC  (dim1s:dim1e,dim2), &
          flx_fld%ULWSFC  (dim1s:dim1e,dim2), &
          flx_fld%GFLUX   (dim1s:dim1e,dim2), &
          flx_fld%RUNOFF  (dim1s:dim1e,dim2), &
          flx_fld%EP      (dim1s:dim1e,dim2), &
          flx_fld%CLDWRK  (dim1s:dim1e,dim2), &
          flx_fld%DUGWD   (dim1s:dim1e,dim2), &
          flx_fld%DVGWD   (dim1s:dim1e,dim2), &
          flx_fld%PSMEAN  (dim1s:dim1e,dim2), &
          flx_fld%GESHEM  (dim1s:dim1e,dim2), &
          !TBH:  added RAINC, EVAP, HFLX for FIM
          flx_fld%RAINC   (dim1s:dim1e,dim2), &
          flx_fld%EVAP    (dim1s:dim1e,dim2), &
          flx_fld%HFLX    (dim1s:dim1e,dim2), &
          flx_fld%BENGSH  (dim1s:dim1e,dim2), &
          flx_fld%SFCNSW  (dim1s:dim1e,dim2), &
          flx_fld%SFCDLW  (dim1s:dim1e,dim2), &
          flx_fld%TSFLW   (dim1s:dim1e,dim2), &
          flx_fld%PSURF   (dim1s:dim1e,dim2), &
          flx_fld%U10M    (dim1s:dim1e,dim2), &
          flx_fld%V10M    (dim1s:dim1e,dim2), &
          flx_fld%HPBL    (dim1s:dim1e,dim2), &
          flx_fld%PWAT    (dim1s:dim1e,dim2), &
!TBH:  not used by FIM yet
!         flx_fld%CHH     (dim1s:dim1e,dim2), &
!         flx_fld%CMM     (dim1s:dim1e,dim2), &
! jbao new gfs phys
          flx_fld%EPI     (dim1s:dim1e,dim2), &
!         flx_fld%DLWSFCI (dim1s:dim1e,dim2), &
!         flx_fld%ULWSFCI (dim1s:dim1e,dim2), &
!         flx_fld%USWSFCI (dim1s:dim1e,dim2), &
!         flx_fld%DSWSFCI (dim1s:dim1e,dim2), &
!         flx_fld%DTSFCI  (dim1s:dim1e,dim2), &
!         flx_fld%DQSFCI  (dim1s:dim1e,dim2), &
!         flx_fld%GFLUXI  (dim1s:dim1e,dim2), &
!         flx_fld%SRUNOFF (dim1s:dim1e,dim2), &
!         flx_fld%T1      (dim1s:dim1e,dim2), &
!         flx_fld%Q1      (dim1s:dim1e,dim2), &
!         flx_fld%U1      (dim1s:dim1e,dim2), &
!         flx_fld%V1      (dim1s:dim1e,dim2), &
!         flx_fld%ZLVL    (dim1s:dim1e,dim2), &
!         flx_fld%EVBSA   (dim1s:dim1e,dim2), &
!         flx_fld%EVCWA   (dim1s:dim1e,dim2), &
!         flx_fld%TRANSA  (dim1s:dim1e,dim2), &
!         flx_fld%SBSNOA  (dim1s:dim1e,dim2), &
!         flx_fld%SNOWCA  (dim1s:dim1e,dim2), &
!         flx_fld%SOILM   (dim1s:dim1e,dim2), &
          stat=iret)

    if(iret.ne.0) iret=-4
    return
  end subroutine

    subroutine flx_init(flx_fld, iret)

    USE gfs_physics_sfc_flx_mod, ONLY: Flx_Var_Data
    implicit none
    TYPE(Flx_Var_Data), INTENT(inout) :: flx_fld 

    integer, intent(out)             :: iret
!
    flx_fld%TMPMIN  = 1.e4
    flx_fld%TMPMAX  = 0.
    flx_fld%GESHEM  = 0.
    !TBH:  added RAINC for FIM
    flx_fld%RAINC   = 0.
    flx_fld%BENGSH  = 0.
    flx_fld%DUSFC   = 0.
    flx_fld%DVSFC   = 0.
    flx_fld%DTSFC   = 0.
    flx_fld%DQSFC   = 0.
    flx_fld%DLWSFC  = 0.
    flx_fld%ULWSFC  = 0.
    flx_fld%GFLUX   = 0.
!
    flx_fld%RUNOFF  = 0.
    flx_fld%EP      = 0.
    flx_fld%CLDWRK  = 0.
    flx_fld%DUGWD   = 0.
    flx_fld%DVGWD   = 0.
    flx_fld%PSMEAN  = 0.
!
!TBH:  not used by FIM yet
!   flx_fld%EVBSA   = 0.
!   flx_fld%EVCWA   = 0.
!   flx_fld%TRANSA  = 0.
!   flx_fld%SBSNOA  = 0.
!   flx_fld%SNOWCA  = 0.
!   flx_fld%SRUNOFF = 0.

     return
  end subroutine
!SMS$IGNORE END
 END MODULE gfs_physics_sfc_flx_set_mod
