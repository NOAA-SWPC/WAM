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
!  March 2008      Y.-T. Hou         add Sunshine_Duration (suntim) to Flx_Var_Data
!  Jan 2009        Moorthi           add Ho Chun's changes
!  Apr 2009        Y.-T. Hou         add surface lw emissivity (sfcemis)
!  Nov 2009        Sarah Lu,         add rain and rainc
!  Sep 2010        Sarah Lu,         add wet1
!  June 2013       H Chuang          add sr
!  Oct 2013        Xingren Wu, add DUSFCI/DVSFCI
!
! !INTERFACE:
!
 MODULE gfs_physics_sfc_flx_mod

 use machine , only : kind_phys

 IMPLICIT none

 TYPE Sfc_Var_Data
    real(kind=kind_phys),pointer:: tsea(:,:)=>null()
    real(kind=kind_phys),pointer:: smc(:,:,:)=>null()
    real(kind=kind_phys),pointer:: weasd(:,:)=>null()
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
    real(kind=kind_phys),pointer:: tisfc(:,:)=>null()
    real(kind=kind_phys),pointer:: tprcp(:,:)=>null()
    real(kind=kind_phys),pointer:: srflag(:,:)=>null()
    real(kind=kind_phys),pointer:: snwdph(:,:)=>null()
    real(kind=kind_phys),pointer:: slc(:,:,:)=>null()
    real(kind=kind_phys),pointer:: shdmin(:,:)=>null()
    real(kind=kind_phys),pointer:: shdmax(:,:)=>null()
    real(kind=kind_phys),pointer:: slope(:,:)=>null()
    real(kind=kind_phys),pointer:: snoalb(:,:)=>null()
    real(kind=kind_phys),pointer:: oro(:,:)=>null()
    real(kind=kind_phys),pointer:: oro_uf(:,:)=>null()
 end type Sfc_Var_Data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 TYPE Flx_Var_Data
    real(kind=kind_phys),pointer:: SFCDSW(:,:)=>null()
    real(kind=kind_phys),pointer:: COSZEN(:,:)=>null()
    real(kind=kind_phys),pointer:: TMPMIN(:,:)=>null()
    real(kind=kind_phys),pointer:: TMPMAX(:,:)=>null()
!jwang add spfhmax/spfhmin
    real(kind=kind_phys),pointer:: SPFHMIN(:,:)=>null()
    real(kind=kind_phys),pointer:: SPFHMAX(:,:)=>null()
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
    real(kind=kind_phys),pointer:: BENGSH(:,:)=>null()
    real(kind=kind_phys),pointer:: SFCNSW(:,:)=>null()
    real(kind=kind_phys),pointer:: SFCDLW(:,:)=>null()
    real(kind=kind_phys),pointer:: TSFLW(:,:)=>null()
    real(kind=kind_phys),pointer:: PSURF(:,:)=>null()
    real(kind=kind_phys),pointer:: U10M(:,:)=>null()
    real(kind=kind_phys),pointer:: V10M(:,:)=>null()
    real(kind=kind_phys),pointer:: HPBL(:,:)=>null()
    real(kind=kind_phys),pointer:: PWAT(:,:)=>null()
    real(kind=kind_phys),pointer:: CHH(:,:)=>null()
    real(kind=kind_phys),pointer:: CMM(:,:)=>null()
    real(kind=kind_phys),pointer:: EPI(:,:)=>null()
    real(kind=kind_phys),pointer:: DLWSFCI(:,:)=>null()
    real(kind=kind_phys),pointer:: ULWSFCI(:,:)=>null()
    real(kind=kind_phys),pointer:: USWSFCI(:,:)=>null()
    real(kind=kind_phys),pointer:: DSWSFCI(:,:)=>null()
    real(kind=kind_phys),pointer:: DUSFCI(:,:)=>null()
    real(kind=kind_phys),pointer:: DVSFCI(:,:)=>null()
    real(kind=kind_phys),pointer:: DTSFCI(:,:)=>null()
    real(kind=kind_phys),pointer:: DQSFCI(:,:)=>null()
    real(kind=kind_phys),pointer:: GFLUXI(:,:)=>null()
    real(kind=kind_phys),pointer:: SRUNOFF(:,:)=>null()
    real(kind=kind_phys),pointer:: T1(:,:)=>null()
    real(kind=kind_phys),pointer:: Q1(:,:)=>null()
    real(kind=kind_phys),pointer:: U1(:,:)=>null()
    real(kind=kind_phys),pointer:: V1(:,:)=>null()
    real(kind=kind_phys),pointer:: ZLVL(:,:)=>null()
    real(kind=kind_phys),pointer:: EVBSA(:,:)=>null()
    real(kind=kind_phys),pointer:: EVCWA(:,:)=>null()
    real(kind=kind_phys),pointer:: TRANSA(:,:)=>null()
    real(kind=kind_phys),pointer:: SBSNOA(:,:)=>null()
    real(kind=kind_phys),pointer:: SNOWCA(:,:)=>null()
    real(kind=kind_phys),pointer:: SOILM(:,:)=>null()
    real(kind=kind_phys),pointer:: SNOHFA(:,:)=>null()
    real(kind=kind_phys),pointer:: SMCWLT2(:,:)=>null()
    real(kind=kind_phys),pointer:: SMCREF2(:,:)=>null()
    real(kind=kind_phys),pointer:: suntim(:,:)=>null()       ! sunshine durationtime
    real(kind=kind_phys),pointer:: sfcemis(:,:)=>null()      ! surface emissivity
    real(kind=kind_phys),pointer:: RAIN(:,:)=>null()
    real(kind=kind_phys),pointer:: RAINC(:,:)=>null()
    real(kind=kind_phys),pointer:: WET1(:,:)=>null()
    real(kind=kind_phys),pointer:: sr(:,:)=>null()
 end type Flx_Var_Data
 END MODULE gfs_physics_sfc_flx_mod
