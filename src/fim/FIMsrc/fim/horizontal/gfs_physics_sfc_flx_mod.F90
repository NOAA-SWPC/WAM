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
 MODULE gfs_physics_sfc_flx_mod
!SMS$IGNORE BEGIN

 use machine , only : kind_phys

 IMPLICIT none

 TYPE Sfc_Var_Data
    real(kind=kind_phys),pointer:: tsea(:,:)=>null()
    real(kind=kind_phys),pointer:: smc(:,:,:)=>null()
    real(kind=kind_phys),pointer:: sheleg(:,:)=>null()
!TBH:  not used by FIM yet
! jbao new gfs phys
    real(kind=kind_phys),pointer:: sncovr(:,:)=>null()
    real(kind=kind_phys),pointer:: stc(:,:,:)=>null()
    real(kind=kind_phys),pointer:: tg3(:,:)=>null()
    real(kind=kind_phys),pointer:: zorl(:,:)=>null()
    real(kind=kind_phys),pointer:: cv(:,:)=>null()
    real(kind=kind_phys),pointer:: cvb(:,:)=>null()
    real(kind=kind_phys),pointer:: cvt(:,:)=>null()
    real(kind=kind_phys),pointer:: alvsf(:,:)=>null()
    real(kind=kind_phys),pointer:: alvwf(:,:)=>null()
    real(kind=kind_phys),pointer:: alnsf(:,:)=>null()
    real(kind=kind_phys),pointer:: alnwf(:,:)=>null()
    real(kind=kind_phys),pointer:: slmsk(:,:)=>null()
    real(kind=kind_phys),pointer:: vfrac(:,:)=>null()
    real(kind=kind_phys),pointer:: canopy(:,:)=>null()
    real(kind=kind_phys),pointer:: f10m(:,:)=>null()
    real(kind=kind_phys),pointer:: t2m(:,:)=>null()
    real(kind=kind_phys),pointer:: q2m(:,:)=>null()
    real(kind=kind_phys),pointer:: vtype(:,:)=>null()
    real(kind=kind_phys),pointer:: stype(:,:)=>null()
    real(kind=kind_phys),pointer:: facsf(:,:)=>null()
    real(kind=kind_phys),pointer:: facwf(:,:)=>null()
    real(kind=kind_phys),pointer:: uustar(:,:)=>null()
    real(kind=kind_phys),pointer:: ffmm(:,:)=>null()
    real(kind=kind_phys),pointer:: ffhh(:,:)=>null()
    real(kind=kind_phys),pointer:: hice(:,:)=>null()
    real(kind=kind_phys),pointer:: fice(:,:)=>null()
!TBH:  not used by FIM yet
! jbao new gfs phys
    real(kind=kind_phys),pointer:: tisfc(:,:)=>null()
    real(kind=kind_phys),pointer:: tprcp(:,:)=>null()
    real(kind=kind_phys),pointer:: srflag(:,:)=>null()
    real(kind=kind_phys),pointer:: snwdph(:,:)=>null()
    real(kind=kind_phys),pointer:: slc(:,:,:)=>null()
    real(kind=kind_phys),pointer:: shdmin(:,:)=>null()
    real(kind=kind_phys),pointer:: shdmax(:,:)=>null()
    real(kind=kind_phys),pointer:: slope(:,:)=>null()
    real(kind=kind_phys),pointer:: snoalb(:,:)=>null()
!TBH:  not used by FIM yet
!   real(kind=kind_phys),pointer:: oro(:,:)=>null()
 end type Sfc_Var_Data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 TYPE Flx_Var_Data
    real(kind=kind_phys),pointer:: SFCDSW(:,:)=>null()
    real(kind=kind_phys),pointer:: COSZEN(:,:)=>null()
    real(kind=kind_phys),pointer:: TMPMIN(:,:)=>null()
    real(kind=kind_phys),pointer:: TMPMAX(:,:)=>null()
    real(kind=kind_phys),pointer:: DUSFC(:,:)=>null()
    real(kind=kind_phys),pointer:: DVSFC(:,:)=>null()
    real(kind=kind_phys),pointer:: DTSFC(:,:)=>null()
    real(kind=kind_phys),pointer:: DQSFC(:,:)=>null()
    real(kind=kind_phys),pointer:: DLWSFC(:,:)=>null()
    real(kind=kind_phys),pointer:: ULWSFC(:,:)=>null()
    real(kind=kind_phys),pointer:: GFLUX(:,:)=>null()
    real(kind=kind_phys),pointer:: RUNOFF(:,:)=>null()
    real(kind=kind_phys),pointer:: EP(:,:)=>null()
    real(kind=kind_phys),pointer:: CLDWRK(:,:)=>null()
    real(kind=kind_phys),pointer:: DUGWD(:,:)=>null()
    real(kind=kind_phys),pointer:: DVGWD(:,:)=>null()
    real(kind=kind_phys),pointer:: PSMEAN(:,:)=>null()
    real(kind=kind_phys),pointer:: GESHEM(:,:)=>null()
    !TBH:  added RAINC, EVAP, HFLX for FIM
    real(kind=kind_phys),pointer:: RAINC(:,:)=>null()
    real(kind=kind_phys),pointer:: EVAP(:,:)=>null()
    real(kind=kind_phys),pointer:: HFLX(:,:)=>null()
    real(kind=kind_phys),pointer:: BENGSH(:,:)=>null()
    real(kind=kind_phys),pointer:: SFCNSW(:,:)=>null()
    real(kind=kind_phys),pointer:: SFCDLW(:,:)=>null()
    real(kind=kind_phys),pointer:: TSFLW(:,:)=>null()
    real(kind=kind_phys),pointer:: PSURF(:,:)=>null()
    real(kind=kind_phys),pointer:: U10M(:,:)=>null()
    real(kind=kind_phys),pointer:: V10M(:,:)=>null()
    real(kind=kind_phys),pointer:: HPBL(:,:)=>null()
    real(kind=kind_phys),pointer:: PWAT(:,:)=>null()
!TBH:  not used by FIM yet
!   real(kind=kind_phys),pointer:: CHH(:,:)=>null()
!   real(kind=kind_phys),pointer:: CMM(:,:)=>null()
! jbao new gfs phys
    real(kind=kind_phys),pointer:: EPI(:,:)=>null()
!   real(kind=kind_phys),pointer:: DLWSFCI(:,:)=>null()
!   real(kind=kind_phys),pointer:: ULWSFCI(:,:)=>null()
!   real(kind=kind_phys),pointer:: USWSFCI(:,:)=>null()
!   real(kind=kind_phys),pointer:: DSWSFCI(:,:)=>null()
!   real(kind=kind_phys),pointer:: DTSFCI(:,:)=>null()
!   real(kind=kind_phys),pointer:: DQSFCI(:,:)=>null()
!   real(kind=kind_phys),pointer:: GFLUXI(:,:)=>null()
!   real(kind=kind_phys),pointer:: SRUNOFF(:,:)=>null()
!   real(kind=kind_phys),pointer:: T1(:,:)=>null()
!   real(kind=kind_phys),pointer:: Q1(:,:)=>null()
!   real(kind=kind_phys),pointer:: U1(:,:)=>null()
!   real(kind=kind_phys),pointer:: V1(:,:)=>null()
!   real(kind=kind_phys),pointer:: ZLVL(:,:)=>null()
!   real(kind=kind_phys),pointer:: EVBSA(:,:)=>null()
!   real(kind=kind_phys),pointer:: EVCWA(:,:)=>null()
!   real(kind=kind_phys),pointer:: TRANSA(:,:)=>null()
!   real(kind=kind_phys),pointer:: SBSNOA(:,:)=>null()
!   real(kind=kind_phys),pointer:: SNOWCA(:,:)=>null()
!   real(kind=kind_phys),pointer:: SOILM(:,:)=>null()
 end type Flx_Var_Data
!SMS$IGNORE END
 END MODULE gfs_physics_sfc_flx_mod
