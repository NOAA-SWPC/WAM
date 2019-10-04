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
!  Nov 2009        Sarah Lu, add rain and rainc
!  Sep 2010        Sarah Lu, add wet1
!  Nov 2011        Sarah Lu, init wet1
!  Oct 2012        S. Moorthi add oro_uf
!  Aug 2013        S. Moorthi - Adding Huiya's addition - sr variable from gfs
!  Oct 2013        Xingren Wu add flx_fld%DUSFCI/DVSFCI
!  Aug 2015        S. Moorthi - add subroutine sfc_init
!
! !INTERFACE:
!
 MODULE gfs_physics_sfc_flx_set_mod

 use machine , only : kind_phys

 IMPLICIT none

    contains
    subroutine sfcvar_aldata(dim1, dim2, dim3, sfc_fld, iret)

    USE gfs_physics_sfc_flx_mod, ONLY: Sfc_Var_Data
    implicit none

    TYPE(Sfc_Var_Data), INTENT(inout) :: sfc_fld 
    integer, intent(in)               :: dim1, dim2, dim3

    integer, intent(out)             :: iret
!
allocate(                                  &
           sfc_fld%tsea   (dim1,dim2),     &
           sfc_fld%smc    (dim3,dim1,dim2),&
           sfc_fld%weasd  (dim1,dim2),     &
           sfc_fld%sncovr (dim1,dim2),     &
           sfc_fld%stc    (dim3,dim1,dim2),&
           sfc_fld%tg3    (dim1,dim2),     &
           sfc_fld%zorl   (dim1,dim2),     &
           sfc_fld%cv     (dim1,dim2),     &
           sfc_fld%cvb    (dim1,dim2),     &
           sfc_fld%cvt    (dim1,dim2),     &
           sfc_fld%alvsf  (dim1,dim2),     &
           sfc_fld%alvwf  (dim1,dim2),     &
           sfc_fld%alnsf  (dim1,dim2),     &
           sfc_fld%alnwf  (dim1,dim2),     &
           sfc_fld%slmsk  (dim1,dim2),     &
           sfc_fld%vfrac  (dim1,dim2),     &
           sfc_fld%canopy (dim1,dim2),     &
           sfc_fld%f10m   (dim1,dim2),     &
           sfc_fld%t2m    (dim1,dim2),     &
           sfc_fld%q2m    (dim1,dim2),     &
           sfc_fld%vtype  (dim1,dim2),     &
           sfc_fld%stype  (dim1,dim2),     &
           sfc_fld%facsf  (dim1,dim2),     &
           sfc_fld%facwf  (dim1,dim2),     &
           sfc_fld%uustar (dim1,dim2),     &
           sfc_fld%ffmm   (dim1,dim2),     &
           sfc_fld%ffhh   (dim1,dim2),     &
           sfc_fld%hice   (dim1,dim2),     &
           sfc_fld%fice   (dim1,dim2),     &
           sfc_fld%tisfc  (dim1,dim2),     &
           sfc_fld%tprcp  (dim1,dim2),     &
           sfc_fld%srflag (dim1,dim2),     &
           sfc_fld%snwdph (dim1,dim2),     &
           sfc_fld%slc    (dim3,dim1,dim2),&
           sfc_fld%shdmin (dim1,dim2),     &
           sfc_fld%shdmax (dim1,dim2),     &
           sfc_fld%slope  (dim1,dim2),     &
           sfc_fld%snoalb (dim1,dim2),     &
           sfc_fld%oro    (dim1,dim2),     &
           sfc_fld%oro_uf (dim1,dim2),     &
           stat=iret)
    if(iret.ne.0) iret=-3
    return
  end subroutine
    subroutine flxvar_aldata(dim1, dim2, flx_fld, iret)

    USE gfs_physics_sfc_flx_mod, ONLY: Flx_Var_Data
    implicit none
    TYPE(Flx_Var_Data), INTENT(inout) :: flx_fld 
    integer, intent(in)               :: dim1, dim2

    integer, intent(out)             :: iret
!
    allocate(                          &
          flx_fld%SFCDSW  (dim1,dim2), &
          flx_fld%COSZEN  (dim1,dim2), &
          flx_fld%TMPMIN  (dim1,dim2), &
          flx_fld%TMPMAX  (dim1,dim2), &
!jwang add spfhmax/spfhmin
          flx_fld%SPFHMIN (dim1,dim2), &
          flx_fld%SPFHMAX (dim1,dim2), &
          flx_fld%DUSFC   (dim1,dim2), &
          flx_fld%DVSFC   (dim1,dim2), &
          flx_fld%DTSFC   (dim1,dim2), &
          flx_fld%DQSFC   (dim1,dim2), &
          flx_fld%DLWSFC  (dim1,dim2), &
          flx_fld%ULWSFC  (dim1,dim2), &
          flx_fld%GFLUX   (dim1,dim2), &
          flx_fld%RUNOFF  (dim1,dim2), &
          flx_fld%EP      (dim1,dim2), &
          flx_fld%CLDWRK  (dim1,dim2), &
          flx_fld%DUGWD   (dim1,dim2), &
          flx_fld%DVGWD   (dim1,dim2), &
          flx_fld%PSMEAN  (dim1,dim2), &
          flx_fld%GESHEM  (dim1,dim2), &
          flx_fld%BENGSH  (dim1,dim2), &
          flx_fld%SFCNSW  (dim1,dim2), &
          flx_fld%SFCDLW  (dim1,dim2), &
          flx_fld%TSFLW   (dim1,dim2), &
          flx_fld%PSURF   (dim1,dim2), &
          flx_fld%U10M    (dim1,dim2), &
          flx_fld%V10M    (dim1,dim2), &
          flx_fld%HPBL    (dim1,dim2), &
          flx_fld%PWAT    (dim1,dim2), &
          flx_fld%CHH     (dim1,dim2), &
          flx_fld%CMM     (dim1,dim2), &
          flx_fld%EPI     (dim1,dim2), &
          flx_fld%DLWSFCI (dim1,dim2), &
          flx_fld%ULWSFCI (dim1,dim2), &
          flx_fld%USWSFCI (dim1,dim2), &
          flx_fld%DSWSFCI (dim1,dim2), &
          flx_fld%DUSFCI  (dim1,dim2), &
          flx_fld%DVSFCI  (dim1,dim2), &
          flx_fld%DTSFCI  (dim1,dim2), &
          flx_fld%DQSFCI  (dim1,dim2), &
          flx_fld%GFLUXI  (dim1,dim2), &
          flx_fld%SRUNOFF (dim1,dim2), &
          flx_fld%T1      (dim1,dim2), &
          flx_fld%Q1      (dim1,dim2), &
          flx_fld%U1      (dim1,dim2), &
          flx_fld%V1      (dim1,dim2), &
          flx_fld%ZLVL    (dim1,dim2), &
          flx_fld%EVBSA   (dim1,dim2), &
          flx_fld%EVCWA   (dim1,dim2), &
          flx_fld%TRANSA  (dim1,dim2), &
          flx_fld%SBSNOA  (dim1,dim2), &
          flx_fld%SNOWCA  (dim1,dim2), &
          flx_fld%SOILM   (dim1,dim2), &
          flx_fld%SNOHFA  (dim1,dim2), &
          flx_fld%SMCWLT2 (dim1,dim2), &
          flx_fld%SMCREF2 (dim1,dim2), &
          flx_fld%suntim  (dim1,dim2), &                !yth mar/08
          flx_fld%sfcemis (dim1,dim2), &                !yth apr/09
          flx_fld%RAIN    (dim1,dim2), &
          flx_fld%RAINC   (dim1,dim2), &
          flx_fld%WET1    (dim1,dim2), &
          flx_fld%sr      (dim1,dim2), &

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
!jwang add spfhmax/spfhmin
    flx_fld%SPFHMIN = 1.e10
    flx_fld%SPFHMAX = 0.
    flx_fld%GESHEM  = 0.
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
    flx_fld%EVBSA   = 0.
    flx_fld%EVCWA   = 0.
    flx_fld%TRANSA  = 0.
    flx_fld%SBSNOA  = 0.
    flx_fld%SNOWCA  = 0.
    flx_fld%SRUNOFF = 0.
    flx_fld%SNOHFA  = 0.
!jw
    flx_fld%SUNTIM  = 0.
!for gocart
    flx_fld%wet1    = 0.
    flx_fld%sr      = 0.

     iret = 0

     return
  end subroutine
    subroutine sfc_init(sfc_fld, iret)

    USE gfs_physics_sfc_flx_mod, ONLY: Sfc_Var_Data
    implicit none

    TYPE(Sfc_Var_Data), INTENT(inout) :: sfc_fld 

    integer, intent(out)             :: iret
!
    sfc_fld%tsea   = 0.0
    sfc_fld%smc    = 0.0
    sfc_fld%weasd  = 0.0
    sfc_fld%sncovr = 0.0
    sfc_fld%stc    = 0.0
    sfc_fld%tg3    = 0.0
    sfc_fld%zorl   = 0.0
    sfc_fld%cv     = 0.0
    sfc_fld%cvb    = 0.0
    sfc_fld%cvt    = 0.0
    sfc_fld%alvsf  = 0.0
    sfc_fld%alvwf  = 0.0
    sfc_fld%alnsf  = 0.0
    sfc_fld%alnwf  = 0.0
    sfc_fld%slmsk  = 0.0
    sfc_fld%vfrac  = 0.0
    sfc_fld%canopy = 0.0
    sfc_fld%f10m   = 0.0
    sfc_fld%t2m    = 0.0
    sfc_fld%q2m    = 0.0
    sfc_fld%vtype  = 0.0
    sfc_fld%stype  = 0.0
    sfc_fld%facsf  = 0.0
    sfc_fld%facwf  = 0.0
    sfc_fld%uustar = 0.0
    sfc_fld%ffmm   = 0.0
    sfc_fld%ffhh   = 0.0
    sfc_fld%hice   = 0.0
    sfc_fld%fice   = 0.0
    sfc_fld%tisfc  = 0.0
    sfc_fld%tprcp  = 0.0
    sfc_fld%srflag = 0.0
    sfc_fld%snwdph = 0.0
    sfc_fld%slc    = 0.0
    sfc_fld%shdmin = 0.0
    sfc_fld%shdmax = 0.0
    sfc_fld%slope  = 0.0
    sfc_fld%snoalb = 0.0
    sfc_fld%oro    = 0.0
    sfc_fld%oro_uf = 0.0
    if(iret.ne.0) iret=-4
    return
  end subroutine
 END MODULE gfs_physics_sfc_flx_set_mod
