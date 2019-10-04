!
! !MODULE: gfs_physics_aoi_mod  ---        Definition of the A/O/I coupling
!                                          fields in the ESMF internal state.
!
! !DESCRIPTION: gfs_physics_aoi_mod ---    Define the A/O/I coupling variables
!                                          in the ESMF internal state.
!---------------------------------------------------------------------------
! !REVISION HISTORY:
!
!  Feb 2014      Xingren Wu Initial code.
!  May 2015      Xingren Wu Add T/Q/U/V/P/Z from lowest model level
!  Sep 2015      Xingren Wu Add dqsfcin/dtsfcin/dusfcin/dvsfcin/ulwsfcin/
!                               slimskin/tseain/tisfcin/ficein/hicein/hsnoin
!  Oct 2015      Xingren Wu Add snow
!  Jan 2016      Patrick Tripp - added coupled fields for NUOPC/GSM merge
!
! !INTERFACE:
!
 MODULE gfs_physics_aoi_var_mod

 use machine , only : kind_phys

 IMPLICIT none

 TYPE AOI_Var_Data
    real(kind_phys),pointer:: dusfc    (:,:)=>null()
    real(kind_phys),pointer:: dvsfc    (:,:)=>null()
    real(kind_phys),pointer:: dtsfc    (:,:)=>null()
    real(kind_phys),pointer:: dqsfc    (:,:)=>null()
    real(kind_phys),pointer:: dlwsfc   (:,:)=>null()
    real(kind_phys),pointer:: dswsfc   (:,:)=>null()
    real(kind_phys),pointer:: dnirbm   (:,:)=>null()
    real(kind_phys),pointer:: dnirdf   (:,:)=>null()
    real(kind_phys),pointer:: dvisbm   (:,:)=>null()
    real(kind_phys),pointer:: dvisdf   (:,:)=>null()
    real(kind_phys),pointer:: nlwsfc   (:,:)=>null()
    real(kind_phys),pointer:: nswsfc   (:,:)=>null()
    real(kind_phys),pointer:: nnirbm   (:,:)=>null()
    real(kind_phys),pointer:: nnirdf   (:,:)=>null()
    real(kind_phys),pointer:: nvisbm   (:,:)=>null()
    real(kind_phys),pointer:: nvisdf   (:,:)=>null()
    real(kind_phys),pointer:: rain     (:,:)=>null()
    real(kind_phys),pointer:: snow     (:,:)=>null()
    real(kind_phys),pointer:: dusfci   (:,:)=>null()
    real(kind_phys),pointer:: dvsfci   (:,:)=>null()
    real(kind_phys),pointer:: dtsfci   (:,:)=>null()
    real(kind_phys),pointer:: dqsfci   (:,:)=>null()
    real(kind_phys),pointer:: dlwsfci  (:,:)=>null()
    real(kind_phys),pointer:: dswsfci  (:,:)=>null()
    real(kind_phys),pointer:: dnirbmi  (:,:)=>null()
    real(kind_phys),pointer:: dnirdfi  (:,:)=>null()
    real(kind_phys),pointer:: dvisbmi  (:,:)=>null()
    real(kind_phys),pointer:: dvisdfi  (:,:)=>null()
    real(kind_phys),pointer:: nlwsfci  (:,:)=>null()
    real(kind_phys),pointer:: nswsfci  (:,:)=>null()
    real(kind_phys),pointer:: nnirbmi  (:,:)=>null()
    real(kind_phys),pointer:: nnirdfi  (:,:)=>null()
    real(kind_phys),pointer:: nvisbmi  (:,:)=>null()
    real(kind_phys),pointer:: nvisdfi  (:,:)=>null()
    real(kind_phys),pointer:: nirbmi   (:,:)=>null()
    real(kind_phys),pointer:: nirdfi   (:,:)=>null()
    real(kind_phys),pointer:: visbmi   (:,:)=>null()
    real(kind_phys),pointer:: visdfi   (:,:)=>null()
    real(kind_phys),pointer:: t2mi     (:,:)=>null()
    real(kind_phys),pointer:: q2mi     (:,:)=>null()
    real(kind_phys),pointer:: u10mi    (:,:)=>null()
    real(kind_phys),pointer:: v10mi    (:,:)=>null()
    real(kind_phys),pointer:: tseai    (:,:)=>null()
    real(kind_phys),pointer:: psurfi   (:,:)=>null()
    real(kind_phys),pointer:: oro      (:,:)=>null()
    real(kind_phys),pointer:: slimsk   (:,:)=>null()
    real(kind_phys),pointer:: nirbmdi  (:,:)=>null()
    real(kind_phys),pointer:: nirdfdi  (:,:)=>null()
    real(kind_phys),pointer:: visbmdi  (:,:)=>null()
    real(kind_phys),pointer:: visdfdi  (:,:)=>null()
    real(kind_phys),pointer:: nirbmui  (:,:)=>null()
    real(kind_phys),pointer:: nirdfui  (:,:)=>null()
    real(kind_phys),pointer:: visbmui  (:,:)=>null()
    real(kind_phys),pointer:: visdfui  (:,:)=>null()
    real(kind_phys),pointer:: tboti    (:,:)=>null()
    real(kind_phys),pointer:: qboti    (:,:)=>null()
    real(kind_phys),pointer:: uboti    (:,:)=>null()
    real(kind_phys),pointer:: vboti    (:,:)=>null()
    real(kind_phys),pointer:: pboti    (:,:)=>null()
    real(kind_phys),pointer:: zboti    (:,:)=>null()
    real(kind_phys),pointer:: dqsfcin  (:,:)=>null()
    real(kind_phys),pointer:: dtsfcin  (:,:)=>null()
    real(kind_phys),pointer:: dusfcin  (:,:)=>null()
    real(kind_phys),pointer:: dvsfcin  (:,:)=>null()
    real(kind_phys),pointer:: ulwsfcin (:,:)=>null()
    real(kind_phys),pointer:: slimskin (:,:)=>null()
    real(kind_phys),pointer:: tseain   (:,:)=>null()
    real(kind_phys),pointer:: tisfcin  (:,:)=>null()
    real(kind_phys),pointer:: ficein   (:,:)=>null()
    real(kind_phys),pointer:: hicein   (:,:)=>null()
    real(kind_phys),pointer:: hsnoin  (:,:)=>null()
 end type AOI_Var_Data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    contains
    subroutine aoivar_aldata(dim1,dim2,aoi_fld,iret)
       implicit none
       integer, intent(in)               :: dim1, dim2
       type(aoi_var_data),intent(inout)  :: aoi_fld
       integer, intent(out)              :: iret
!
    allocate(                         &
      aoi_fld%dusfc    (dim1,dim2),   &
      aoi_fld%dvsfc    (dim1,dim2),   &
      aoi_fld%dtsfc    (dim1,dim2),   &
      aoi_fld%dqsfc    (dim1,dim2),   &
      aoi_fld%dlwsfc   (dim1,dim2),   &
      aoi_fld%dswsfc   (dim1,dim2),   &
      aoi_fld%dnirbm   (dim1,dim2),   &
      aoi_fld%dnirdf   (dim1,dim2),   &
      aoi_fld%dvisbm   (dim1,dim2),   &
      aoi_fld%dvisdf   (dim1,dim2),   &
      aoi_fld%nlwsfc   (dim1,dim2),   &
      aoi_fld%nswsfc   (dim1,dim2),   &
      aoi_fld%nnirbm   (dim1,dim2),   &
      aoi_fld%nnirdf   (dim1,dim2),   &
      aoi_fld%nvisbm   (dim1,dim2),   &
      aoi_fld%nvisdf   (dim1,dim2),   &
      aoi_fld%rain     (dim1,dim2),   &
      aoi_fld%snow     (dim1,dim2),   &
      aoi_fld%dusfci   (dim1,dim2),   &
      aoi_fld%dvsfci   (dim1,dim2),   &
      aoi_fld%dtsfci   (dim1,dim2),   &
      aoi_fld%dqsfci   (dim1,dim2),   &
      aoi_fld%dlwsfci  (dim1,dim2),   &
      aoi_fld%dswsfci  (dim1,dim2),   &
      aoi_fld%dnirbmi  (dim1,dim2),   &
      aoi_fld%dnirdfi  (dim1,dim2),   &
      aoi_fld%dvisbmi  (dim1,dim2),   &
      aoi_fld%dvisdfi  (dim1,dim2),   &
      aoi_fld%nlwsfci  (dim1,dim2),   &
      aoi_fld%nswsfci  (dim1,dim2),   &
      aoi_fld%nnirbmi  (dim1,dim2),   &
      aoi_fld%nnirdfi  (dim1,dim2),   &
      aoi_fld%nvisbmi  (dim1,dim2),   &
      aoi_fld%nvisdfi  (dim1,dim2),   &
      aoi_fld%nirbmi   (dim1,dim2),   &
      aoi_fld%nirdfi   (dim1,dim2),   &
      aoi_fld%visbmi   (dim1,dim2),   &
      aoi_fld%visdfi   (dim1,dim2),   &
      aoi_fld%t2mi     (dim1,dim2),   &
      aoi_fld%q2mi     (dim1,dim2),   &
      aoi_fld%u10mi    (dim1,dim2),   &
      aoi_fld%v10mi    (dim1,dim2),   &
      aoi_fld%tseai    (dim1,dim2),   &
      aoi_fld%psurfi   (dim1,dim2),   &
      aoi_fld%oro      (dim1,dim2),   &
      aoi_fld%slimsk   (dim1,dim2),   &
      aoi_fld%nirbmdi  (dim1,dim2),   &
      aoi_fld%nirdfdi  (dim1,dim2),   &
      aoi_fld%visbmdi  (dim1,dim2),   &
      aoi_fld%visdfdi  (dim1,dim2),   &
      aoi_fld%nirbmui  (dim1,dim2),   &
      aoi_fld%nirdfui  (dim1,dim2),   &
      aoi_fld%visbmui  (dim1,dim2),   &
      aoi_fld%visdfui  (dim1,dim2),   &
      aoi_fld%tboti    (dim1,dim2),   &
      aoi_fld%qboti    (dim1,dim2),   &
      aoi_fld%uboti    (dim1,dim2),   &
      aoi_fld%vboti    (dim1,dim2),   &
      aoi_fld%pboti    (dim1,dim2),   &
      aoi_fld%zboti    (dim1,dim2),   &
      aoi_fld%dqsfcin  (dim1,dim2),   &
      aoi_fld%dtsfcin  (dim1,dim2),   &
      aoi_fld%dusfcin  (dim1,dim2),   &
      aoi_fld%dvsfcin  (dim1,dim2),   &
      aoi_fld%ulwsfcin (dim1,dim2),   &
      aoi_fld%slimskin (dim1,dim2),   &
      aoi_fld%tseain   (dim1,dim2),   &
      aoi_fld%tisfcin  (dim1,dim2),   &
      aoi_fld%ficein   (dim1,dim2),   &
      aoi_fld%hicein   (dim1,dim2),   &
      aoi_fld%hsnoin   (dim1,dim2),   &
      stat=iret)
    if(iret.ne.0) iret=-3
    return
  end subroutine aoivar_aldata

    subroutine aoivar_init(aoi_fld, iret)

    implicit none

    type(aoi_var_data),intent(inout)  :: aoi_fld
    integer, intent(out)              :: iret
!
      aoi_fld%dusfc    = 0.
      aoi_fld%dvsfc    = 0.
      aoi_fld%dtsfc    = 0.
      aoi_fld%dqsfc    = 0.
      aoi_fld%dlwsfc   = 0.
      aoi_fld%dswsfc   = 0.
      aoi_fld%dnirbm   = 0.
      aoi_fld%dnirdf   = 0.
      aoi_fld%dvisbm   = 0.
      aoi_fld%dvisdf   = 0.
      aoi_fld%nlwsfc   = 0.
      aoi_fld%nswsfc   = 0.
      aoi_fld%nnirbm   = 0.
      aoi_fld%nnirdf   = 0.
      aoi_fld%nvisbm   = 0.
      aoi_fld%nvisdf   = 0.
      aoi_fld%rain     = 0.
      aoi_fld%snow     = 0.
      aoi_fld%dusfci   = 0.
      aoi_fld%dvsfci   = 0.
      aoi_fld%dtsfci   = 0.
      aoi_fld%dqsfci   = 0.
      aoi_fld%dlwsfci  = 0.
      aoi_fld%dswsfci  = 0.
      aoi_fld%dnirbmi  = 0.
      aoi_fld%dnirdfi  = 0.
      aoi_fld%dvisbmi  = 0.
      aoi_fld%dvisdfi  = 0.
      aoi_fld%nlwsfci  = 0.
      aoi_fld%nswsfci  = 0.
      aoi_fld%nnirbmi  = 0.
      aoi_fld%nnirdfi  = 0.
      aoi_fld%nvisbmi  = 0.
      aoi_fld%nvisdfi  = 0.
      aoi_fld%nirbmi   = 0.
      aoi_fld%nirdfi   = 0.
      aoi_fld%visbmi   = 0.
      aoi_fld%visdfi   = 0.
      aoi_fld%t2mi     = 0.
      aoi_fld%q2mi     = 0.
      aoi_fld%u10mi    = 0.
      aoi_fld%v10mi    = 0.
      aoi_fld%tseai    = 0.
      aoi_fld%psurfi   = 0.
      aoi_fld%oro      = 0.
      aoi_fld%slimsk   = 0.
      aoi_fld%nirbmdi  = 0.
      aoi_fld%nirdfdi  = 0.
      aoi_fld%visbmdi  = 0.
      aoi_fld%visdfdi  = 0.
      aoi_fld%nirbmui  = 0.
      aoi_fld%nirdfui  = 0.
      aoi_fld%visbmui  = 0.
      aoi_fld%visdfui  = 0.
      aoi_fld%tboti    = 0.
      aoi_fld%qboti    = 0.
      aoi_fld%uboti    = 0.
      aoi_fld%vboti    = 0.
      aoi_fld%pboti    = 0.
      aoi_fld%zboti    = 0.
      aoi_fld%dqsfcin  = 0.
      aoi_fld%dtsfcin  = 0.
      aoi_fld%dusfcin  = 0.
      aoi_fld%dvsfcin  = 0.
      aoi_fld%ulwsfcin = 0.
      aoi_fld%slimskin = 0.
      aoi_fld%tseain   = 0.
      aoi_fld%tisfcin  = 0.
      aoi_fld%ficein   = 0.
      aoi_fld%hicein   = 0.
      aoi_fld%hsnoin   = 0.

     iret = 0

     return
  end subroutine aoivar_init
 END MODULE gfs_physics_aoi_var_mod
