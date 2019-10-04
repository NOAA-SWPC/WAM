module module_do_physics_one_step

!*********************************************************************
!     do_physics_one_step
!	Calculates column forcing for global fim
!	12/21/2005 - Alexander E. MacDonald     - original version
!	05/01/2006 - Jian-Wen Bao               - modified for GFS physics
!       04/14/2008 - Stan Benjamin, John Brown  - modifications
!                      for introduction of virtual pot temp for temp prog 
!                      variable instead of previous non-virtual pot temp
!       02/26/2009 - Tom Henderson       - moved here from physics.F90 to more 
!                                          closely match new GFS r3038
!       07/21/2009 - Jian-wen Bao        - change to random number generator
!                      for xkt2 (cloud-top height) instead of previous 0.6 constant
!                      This follows NCEP's use of this random number generator
!                      (Mersenne twister) and appears to qualitatively improve
!                      forecasts in tropics.
!*********************************************************************

contains

subroutine do_physics_one_step(dtp, kdt, phour,                       &
                 gfs_ps, gfs_dp, gfs_dpdt, gfs_p, gfs_u, gfs_v,       &
                 gfs_t, gfs_q, gfs_oz, gfs_cld,                       &
                 sfc_fld, flx_fld,                                    &
                 XLON,XLAT,COSZDG,                                    &
                 HPRIME, SWH, HLW, FLUXR, SFALB, SLAG, SDEC, CDEC,    &
                 phy_f3d, phy_f2d, NBLCK,                             &
                 CLDCOV,                                              &
                 nvl, LATS_NODE_R, NMTVR, num_p3d, num_p2d, NFXR,     &
                 nip, lsoil, GravityWaveDrag, CallRadiation,          &
                 yyyymmddhhmm,inv_perm,skip_cu_physics,               &
                 skip_mp_physics,skip_chem,ipn,sscal,ext_cof,asymp,   &
                 extlw_cof)

!SMS$IGNORE BEGIN
USE gfs_physics_sfc_flx_mod, ONLY: Sfc_Var_Data, Flx_Var_Data
use module_initial_chem_namelists,only:chem_opt
use layout1, only: me
!SMS$IGNORE END

USE MACHINE       , ONLY : kind_phys,kind_rad,kind_evod
USE resol_def     , ONLY : levs,levp1
USE module_variables, ONLY:diaga,diagb,u_tdcy_phy,v_tdcy_phy,trc_tdcy_phy
USE module_wrf_variables, ONLY:phys3dwrf,phys2dwrf
!SMS$IGNORE BEGIN
USE mersenne_twister
!SMS$IGNORE END
use module_do_physics_one_step_chem,only:do_physics_one_step_chem

implicit none

! Dimension and type external variables:

real (kind=kind_phys), intent(in   ) :: dtp
integer              , intent(in   ) :: kdt
real (kind=kind_rad) , intent(in   ) :: phour
!SMS$DISTRIBUTE (dh,nip) BEGIN
real (kind=kind_evod), intent(in   ) :: gfs_ps(nip)
real (kind=kind_evod), intent(in   ) :: gfs_dp(nip,nvl)
real (kind=kind_evod), intent(in   ) :: gfs_dpdt(nip,nvl)
real (kind=kind_evod), intent(in   ) :: gfs_p(nip,nvl)
real (kind=kind_evod), intent(inout) :: gfs_u(nip,nvl)
real (kind=kind_evod), intent(inout) :: gfs_v(nip,nvl)
real (kind=kind_evod), intent(inout) :: gfs_t(nip,nvl)
real (kind=kind_evod), intent(inout) :: gfs_q(nip,nvl)
real (kind=kind_evod), intent(inout) :: gfs_oz(nip,nvl)
real (kind=kind_evod), intent(inout) :: gfs_cld(nip,nvl)

real*8 temp_ps(nip)
real*8 temp_dp(nip,nvl)
real*8 temp_dpdt(nip,nvl)
real*8 temp_p(nip,nvl)
real*8 temp_u(nip,nvl)
real*8 temp_v(nip,nvl)
real*8 temp_t(nip,nvl)
real*8 temp_q(nip,nvl)
real*8 temp_oz(nip,nvl)
real*8 temp_cld(nip,nvl)
!SMS$DISTRIBUTE END
TYPE(Sfc_Var_Data)   , intent(inout) :: sfc_fld
TYPE(Flx_Var_Data)   , intent(inout) :: flx_fld
!SMS$DISTRIBUTE (dh,nip) BEGIN
real (kind=kind_rad) , intent(in   ) :: XLON(nip,LATS_NODE_R)
real (kind=kind_rad) , intent(in   ) :: XLAT(nip,LATS_NODE_R)
!TODO: check intents between here and phy_f2d
real (kind=kind_rad) , intent(inout) :: COSZDG(nip,LATS_NODE_R)
real (kind=kind_rad) , intent(inout) :: HPRIME(NMTVR,nip,LATS_NODE_R)
real (kind=kind_rad) , intent(inout) :: SWH(nip,nvl,NBLCK,LATS_NODE_R)
real (kind=kind_rad) , intent(inout) :: HLW(nip,nvl,NBLCK,LATS_NODE_R)
real (kind=kind_rad) , intent(inout) :: FLUXR(NFXR,nip,LATS_NODE_R)
real (kind=kind_rad) , intent(inout) :: SFALB(nip,LATS_NODE_R)
real (kind=kind_evod), intent(inout) :: SLAG(nip,LATS_NODE_R)
real (kind=kind_evod), intent(inout) :: SDEC(nip,LATS_NODE_R)
real (kind=kind_evod), intent(inout) :: CDEC(nip,LATS_NODE_R)
real (kind=kind_rad) , intent(inout) :: phy_f3d(nip,nvl,NBLCK,LATS_NODE_R,num_p3d)
real (kind=kind_rad) , intent(inout) :: phy_f2d(nip,LATS_NODE_R,num_p2d)
integer              , intent(in   ) :: NBLCK
real (kind=kind_rad) , intent(inout) :: CLDCOV(nvl,nip,LATS_NODE_R)
!SMS$DISTRIBUTE END
integer              , intent(in   ) :: nvl
integer              , intent(in   ) :: LATS_NODE_R
integer              , intent(in   ) :: NMTVR
integer              , intent(in   ) :: num_p3d
integer              , intent(in   ) :: num_p2d
integer              , intent(in   ) :: NFXR
integer              , intent(in   ) :: nip
integer              , intent(in   ) :: lsoil
logical              , intent(in   ) :: GravityWaveDrag
integer              , intent(in   ) :: CallRadiation
CHARACTER(len=12)    , intent(in   ) :: yyyymmddhhmm
logical              , intent(in   ) :: skip_cu_physics,skip_mp_physics,skip_chem
integer              , intent(inout) :: ipn
!sms$distribute (dh,2) begin
real,intent(inout)::sscal(:,:,:),ext_cof(:,:,:),asymp(:,:,:),extlw_cof(:,:,:)
!sms$distribute end

! Local variables
!----------------------------------------------------------------------
!

!Parameters and arrays used in the GFS physics
!----------------------------------------------------------------------
! TODO:  move ntrac and nrcm to module resol_def
integer IM, IX, ntrac, nrcm
parameter (IM = 1, IX = 1)
parameter (ntrac = 3, nrcm = 1)
integer  ncld,ntoz,ntcw,lonf,latg, jcap,nlons(im)
parameter (ntcw = 3)

integer levshc(im), levshcm     ! Needed for pry version
LOGICAL lssav,lsfwd
    
real(kind=kind_phys) dtf,FHOUR,solhr, prsshc,cubot,cutop

real(kind=kind_phys) UG   (IX,NVL)  , VG   (IX,NVL)      ,  &
		     TG   (IX,NVL)  , qg   (IX,nvl,ntrac),  &
		                      qg1  (IX,nvl,2)    ,  &
		     GT0  (IX,NVL)  , GU0  (IX,NVL)      ,  &
		     GV0  (IX,NVL)  , gq0  (IX,nvl,ntrac),  &
		     DEL  (IX,NVL)  , PRSI (IX,NVL+1)    ,  &
		     PRSL (IX,NVL)  , PRSLK(IX,NVL)      ,  &
		     PRSIK(IX,NVL+1), PHII (IX,NVL+1)    ,  &
		     PHIL (IX,NVL)  , dkt  (im,NVL-1)    ,  &
		     PGR(IM)        ,  XKT2(IM,nrcm)

! TBH:  These arrays are used to avoid recomputation of values 
! TBH:  between calls to GLOOPR and GBPHYS.  
real(kind=kind_phys) PRSL_S(IX,NVL),  PRSI_S(IX,NVL+1)

! Constants
real(kind=kind_phys) RCS2(IM),clstp

! Local variables
real(kind=kind_phys) SINLAT(IM), COSLAT(IM)

real(kind=kind_phys) acv(IM), acvb(IM), acvt(IM)

!TODO:  Move these to new module d3d_def and (maybe) use d3d_zero to set
real(kind=kind_rad) dt3dt(IX,nvl,6),  dq3dt(IX,nvl,7),  &
                    du3dt(IX,nvl,4),  dv3dt(IX,nvl,4)

logical, save ::  nsst_active=.false.
logical, save ::  lggfs3d=.false.
logical old_monin, cnvgwd
logical sashal,newsas,mom4ice,mstrat,trans_trac,cal_pre 
integer KO3,pl_coeff,ncw(2),lsm,lat
logical, save ::  lssav_cc_dummy=.false.
PARAMETER (KO3=46,pl_coeff=2) !ozone levels in climatology
real(kind=kind_phys) sdiaga(im,nvl),sdiagb(im,nvl)
real(kind=kind_phys) poz(KO3), prdout(IX,ko3,pl_coeff), disout(IX,ko3)
real(kind=kind_phys) flgmin(2), ccwf, ctei_rm,suntim(im),SNCOVR(IM)  &
                     ,SPFHMIN(IM),SPFHMAX(IM)
real(kind=kind_phys) ifd(im),time_old(im),time_ins(im),I_Sw(im),  &
                     I_Q(im),I_Qrain(im),I_M(im),I_Tau(im),  &
                     I_Sw_Zw(im),I_Q_Ts(im),I_M_Ts(im),Tref(im),  &
                     dt_cool(im),z_c(im),dt_warm(im),z_w(im),  &
                     c_0(im),c_d(im),w_0(im),w_d(im), dpshc(IM), crtrh(3),  &
                     CHH(IM),CMM(IM),PI(IM),DLWSFCI(IM),ULWSFCI(IM),USWSFCI(IM),  &
                     DSWSFCI(IM),DTSFCI(IM),DQSFCI(IM),GFLUXI(IM),SRUNOFF(IM),T1(IM),  &
                     Q1(IM),U1(IM),V1(IM),ZLVL(IM), TISFC(IM), &
                     EVBSA(IM),EVCWA(IM),TRANSA(IM),SBSNOA(IM),SNOWCA(IM),SOILM(IM),  &
! jbao new gfs phys
                     SNOHFA(IM),SMCWLT2(IM),SMCREF2(IM),   &
                     gsoil(im), gtmp2m(im), gustar(im), gpblh(im), gu10m(im),  &
                     gv10m(im), gzorl(im), goro(im), dkh(IX,nvl), rnp(ix,nvl),  &
                     upd_mf(ix,nvl), dwn_mf(ix,nvl), det_mf(ix,nvl), oro(im)
      real(kind=kind_phys),dimension(IM):: DLWSFC_cc_dummy,ULWSFC_cc_dummy,SWSFC_cc_dummy,XMU_cc_dummy,  &
       DLW_cc_dummy,DSW_cc_dummy,SNW_cc_dummy,LPREC_cc_dummy,  &
       DUSFC_cc_dummy,DVSFC_cc_dummy,DTSFC_cc_dummy,DQSFC_cc_dummy,  &
       PRECR_cc_dummy

logical RAS,LDIAG3D,pre_rad
logical flipv

integer               :: ivl     ! Index for vertical level
integer               :: global_lats_r(1),lonsperlar(1) !JFM
INTEGER hour
logical :: CallRadiationNow
real(kind=kind_phys), parameter :: cons_1p0d9=1.0E9
integer iseed
real :: wrk(1)
type(random_stat), allocatable, save :: rstat(:)
!SMS$DISTRIBUTE (dh,nip) BEGIN
integer              :: inv_perm(nip)
integer              :: seed0   (nip)
real(kind=kind_phys) :: rannum  (nip)
!SMS$DISTRIBUTE END
INTEGER year,month,day
logical, save :: first_rand=.TRUE.

      integer ipnGlobal,its,mype
      logical DiagPrint

real(kind=kind_rad),allocatable::sscal_rad(:,:)
real(kind=kind_rad),allocatable::asymp_rad(:,:)
real(kind=kind_rad),allocatable::ext_cof_rad(:,:)
real(kind=kind_rad),allocatable::extlw_cof_rad(:,:)

!print *,'DEBUG do_physics_one_step():  SIZE(pb2d) = ',SIZE(pb2d)
LEVS             = nvl    !JFM
levp1            = levs+1 !JFM
global_lats_r(1) = 1      !JFM
lonsperlar   (1) = 1      !JFM

CallRadiationNow = (mod(kdt,CallRadiation) == 0 .or. kdt == 1)

READ(UNIT=yyyymmddhhmm(9:10), FMT='(I2)') hour

me   = 0
!SMS$insert call nnt_me(me)
dtf = dtp
if (first_rand) then
    ! set up seeds for each column, first time only
  ALLOCATE( rstat(LBOUND(rannum,1):UBOUND(rannum,1)) )
  READ(UNIT=yyyymmddhhmm(1:4), FMT='(I4)') year
  READ(UNIT=yyyymmddhhmm(5:6), FMT='(I2)') month
  READ(UNIT=yyyymmddhhmm(7:8), FMT='(I2)') day
!NOTE: For large G-levels this serial could be a memory problem.
!      The serial is used to match the random number to the lat-lon for different curves.
!SMS$SERIAL (<inv_perm,IN>, <seed0,OUT> : default=ignore) BEGIN
  do ipn=1,nip
    seed0(inv_perm(ipn)) = year+month+day+hour+ipn
  enddo
!SMS$SERIAL END
!SMS$PARALLEL (dh,ipn) BEGIN
  do ipn=1,nip
    call random_setseed(seed0(ipn),rstat(ipn))
  enddo
  first_rand = .false.
endif
do ipn=1,nip
  call random_number(wrk,rstat(ipn))
  rannum(ipn) = wrk(1)
enddo

call flush(6)
!SMS$PARALLEL END

allocate(sscal_rad(size(sscal,1),size(sscal,3)))
allocate(asymp_rad(size(asymp,1),size(asymp,3)))
allocate(ext_cof_rad(size(ext_cof,1),size(ext_cof,3)))
allocate(extlw_cof_rad(size(extlw_cof,1),size(extlw_cof,3)))

!----------------------------------------------------------------------
!    Loop begins over all horizontal grid points
!----------------------------------------------------------------------

!SMS$PARALLEL (dh,ipn) BEGIN
do ivl=1,nvl
  do ipn=1,nip
    temp_ps(ipn) = gfs_ps(ipn)
    temp_dp(ipn,ivl) = gfs_dp(ipn,ivl)
    temp_dpdt(ipn,ivl) = gfs_dpdt(ipn,ivl)
    temp_p(ipn,ivl) = gfs_p(ipn,ivl)
    temp_u(ipn,ivl) = gfs_u(ipn,ivl)
    temp_v(ipn,ivl) = gfs_v(ipn,ivl)
    temp_t(ipn,ivl) = gfs_t(ipn,ivl)
    temp_q(ipn,ivl) = gfs_q(ipn,ivl)
    temp_oz(ipn,ivl) = gfs_oz(ipn,ivl)
    temp_cld(ipn,ivl) =  gfs_cld(ipn,ivl)
  enddo
enddo

!sms$compare_var(temp_ps, "do_physics.F90 - gfs_ps0 ")
!sms$compare_var(temp_dp, "do_physics.F90 - gfs_dp0 ")
!sms$compare_var(temp_dpdt, "do_physics.F90 - gfs_dpdt0 ")
!sms$compare_var(temp_p, "do_physics.F90 - gfs_p0 ")
!sms$compare_var(temp_oz, "do_physics.F90 - gfs_oz0 ")
!sms$compare_var(temp_cld, "do_physics.F90 - gfs_cld0 ")
!sms$compare_var(temp_u, "do_physics.F90 - gfs_u0 ")
!sms$compare_var(temp_v, "do_physics.F90 - gfs_v0 ")
!sms$compare_var(temp_t, "do_physics.F90 - gfs_t0 ")
!sms$compare_var(temp_q, "do_physics.F90 - gfs_q0 ")

if(CallRadiationNow) then ! Call radiation
  FLUXR=0.0
endif

do ipn=1,nip
  PRSI_S(1,1)       = gfs_ps(ipn)
  do ivl=1,nvl
    PRSI_S(1,ivl+1) = PRSI_S(1,ivl) - gfs_dp(ipn,ivl)
  enddo
  do ivl=1,nvl
    PRSL_S  (1,ivl) = 0.5*(PRSI_S(1,ivl)+PRSI_S(1,ivl+1))
  enddo
  
  do ivl=1,nvl
      !TODO:  can this copy be avoided, eliminating tg?  
    tg      (1,ivl)   = gfs_t(ipn,ivl)
!jbao      PRSLK   (1,ivl)   = theta_nv(ipn,ivl) /tg(1,ivl)
!jbao      PRSLK   (1,ivl)   = 1./PRSLK(1,ivl)
      !TODO:  can these copies be avoided, eliminating qg?  
    qg      (1,ivl,1) = gfs_q(ipn,ivl)
    qg      (1,ivl,2) = gfs_oz(ipn,ivl)
    qg      (1,ivl,3) = gfs_cld(ipn,ivl)
  enddo

    !
    ! Copy data from FIM arrays to set up for GBPHYS call (below)
    !
  do ivl=1,nvl
      !TODO:  can this copy be avoided, eliminating ug and vg?  
    ug(1,ivl) = gfs_u(ipn,ivl)
    vg(1,ivl) = gfs_v(ipn,ivl)
  enddo
  ncld          = 1
    
!----------------------------------------------------------------------
  if(CallRadiationNow) then ! Call radiation
!----------------------------------------------------------------------
    do ivl=1,nvl
        !TODO:  can these copies be avoided, eliminating qg1?  
      qg1  (1,ivl,1) = gfs_oz(ipn,ivl)
      qg1  (1,ivl,2) = gfs_cld(ipn,ivl)
      PRSL (1,ivl)   = PRSL_S(1,ivl)
    enddo
    do ivl=1,nvl+1
      PRSI (1,ivl)   = PRSI_S(1,ivl)
    enddo
    TISFC = sfc_fld%TSEA(ipn,1)
!          write(6,*)'before gloopr ipn ',ipn
!     if (ipn.eq.8570) then
!         write(6,*)'land point before gloopr'
!         write(6,*)'flx_fld%SFCDLW(ipn,1),flx_fld%TSFLW(ipn,1)'
!         write(6,*)flx_fld%SFCDLW(ipn,1),flx_fld%TSFLW(ipn,1)
!         write(6,*)'flx_fld%SFCNSW(ipn,1),flx_fld%SFCDSW(ipn,1)'
!         write(6,*)flx_fld%SFCNSW(ipn,1),flx_fld%SFCDSW(ipn,1)
!         write(6,*)'SWH(ipn,:,1,1)'
!         write(6,*)SWH(ipn,:,1,1)
!         write(6,*)'HLW(ipn,:,1,1)'
!         write(6,*)HLW(ipn,:,1,1)
!         write(6,*)'qg(1,1,1),qg1(1,1,1)'
!         write(6,*)qg(1,1,1),qg1(1,1,1)
!     endif

    if (chem_opt.gt.0) then
      sscal_rad=sscal(:,ipn,:)
      asymp_rad=asymp(:,ipn,:)
      ext_cof_rad=ext_cof(:,ipn,:)
      extlw_cof_rad=extlw_cof(:,ipn,:)
    else
      sscal_rad=999.
    endif

!----------------------------------------------------------------------
!   Call GFS radiation (longwave and shortwave)
!----------------------------------------------------------------------
    CALL GLOOPR                                           &
        (ncld,                                              &  ! jbao ncld needs to passed into gloopr, add to new gloopr and declare
        1,global_lats_r,                                    &
        lonsperlar,                                         &
!jbao not needed by fim        1.0D0,1.0D0,1.0D0,1.0D0,                            &
!jbao not needed by fim        1.0D0,1.0D0,                                        &  ! jbao new gfs, ndexev, ndexod not used now in new gfs
!jbao old gfs        1.0D0,1.0D0,1,1,                          ! jbao new gfs, ndexev, ndexod not used now in new gfs          &
!jbao old gfs       1.0D0,1.0D0,1.0D0,1.0D0,1.0D0,1.0D0,   &  ! jbao newgfs doesn't need plnew_r, plnow_r
!jbao not needed by fim        1.0D0,1.0D0,1.0D0,1.0D0,                            &
        phour,                                              & ! jbao fcst hour
        XLON(ipn,1),XLAT(ipn,1),COSZDG(ipn,1),              &
        flx_fld%COSZEN(ipn,1),                              & ! COSZEN is output and used in gbphys
        sfc_fld%SLMSK(ipn,1),                               &
        sfc_fld%SHELEG(ipn,1),                              &
        sfc_fld%SNCOVR(ipn,1),                              & ! jbao new gfs needs sncovr
        sfc_fld%SNOALB(ipn,1),                              & ! jbao new gfs nees snoalb
        sfc_fld%ZORL(ipn,1),                                &
        sfc_fld%TSEA(ipn,1),                                &
! jbao old gfs needs stc        sfc_fld%STC(1,ipn,1),HPRIME(1,ipn,1),               &
        HPRIME(1,ipn,1),                                    &
        SFALB(ipn,1),                                       & ! SFALB is set in gloopr but then set to 0 before call gbphys
        sfc_fld%ALVSF(ipn,1),                               &
        sfc_fld%ALNSF(ipn,1),                               &
        sfc_fld%ALVWF(ipn,1),                               &
        sfc_fld%ALNWF(ipn,1),                               &
        sfc_fld%FACSF(ipn,1),                               &
        sfc_fld%FACWF(ipn,1),sfc_fld%CV(ipn,1),             &
        sfc_fld%CVT(ipn,1),                                 &
        sfc_fld%CVB(ipn,1),SWH(ipn,:,1,1),                  &
        HLW(ipn,:,1,1),                                     &
        flx_fld%SFCNSW(ipn,1),                              &
        flx_fld%SFCDLW(ipn,1),                              & ! SWH,HLW,SFCNSW,SFCDLW output and used in gbphys
        sfc_fld%FICE(ipn,1) ,                               &
        TISFC,                                              & ! jbao new gfs needs tisfc
        flx_fld%SFCDSW(ipn,1),                              & ! FOR SEA-ICE - XW Nov04, SFCDSW output and used in gbphys
        flx_fld%TSFLW(ipn,1),FLUXR(:,ipn,1),                & ! jbao new gfs does not need cldcov
        phy_f3d(ipn,:,1,1,:),SLAG(ipn,1) ,SDEC(ipn,1),      & ! jbao new gfs needs phy_f3d 
        CDEC(ipn,1),1,KDT,                                  & 
        0.0D0, prsl(1,1),prsi(1,1),prslk(1,1),tg(1,1),      &
        qg(1,1,1),qg1(1,1,1),sscal_rad,                &
        asymp_rad,ext_cof_rad,extlw_cof_rad,yyyymmddhhmm)

!     sscal(:,ipn,:)=sscal_rad
!     asymp(:,ipn,:)=asymp_rad
!     ext_cof(:,ipn,:)=ext_cof_rad

  endif

  if (ipn.eq.8570) then
!         write(6,*)'land point after gloopr'
!         write(6,*)'flx_fld%SFCDLW(ipn,1),flx_fld%TSFLW(ipn,1)'
!         write(6,*)flx_fld%SFCDLW(ipn,1),flx_fld%TSFLW(ipn,1)
!         write(6,*)'flx_fld%SFCNSW(ipn,1),flx_fld%SFCDSW(ipn,1)'
!         write(6,*)flx_fld%SFCNSW(ipn,1),flx_fld%SFCDSW(ipn,1)
!         write(6,*)'SWH(ipn,:,1,1)'
!         write(6,*)SWH(ipn,:,1,1)
!         write(6,*)'HLW(ipn,:,1,1)'
!         write(6,*)HLW(ipn,:,1,1)
!         write(6,*)'qg(1,1,1),qg1(1,1,1)'
!         write(6,*)qg(1,1,1),qg1(1,1,1)
  endif

    !
    ! Local set up for GBPHYS
    !
  do ivl=1,nvl
    PRSL(1,ivl)  = PRSL_S(1,ivl)
  enddo
  do ivl=1,nvl+1
    PRSI(1,ivl)  = PRSI_S(1,ivl)
  enddo
  PGR             = PRSI_S(1,1)
!jbao    prsik (1,1)     = PRSLK(1,1)
!jbao    prsik (1,nvl+1) = PRSLK(1,nvl)
!jbao    do ivl=2,nvl
!jbao      prsik (1,ivl) = 0.5*(PRSLK(1,ivl)+PRSLK(1,ivl-1))
!jbao    enddo

  ntoz      = 2
  lonf      = 200
  latg      = 94
  jcap      = 126
  ras       = .false.
  sashal = .true.   ! jbao from gfs namelist
  newsas = .true.   ! jbao from gfs namelist
  lsm       = 1
  old_monin = .false.
  mstrat = .false.
  mom4ice = .false.
  trans_trac = .false.
  cal_pre = .false.
  suntim = 0.0
  phil = 0.0
  phii = 0.0
  prsik = 0.0
  dpshc(1) = 0.3 * prsi(1,1) ! jbao new GFS physics as of Feb 2010
  crtrh(:) = 0.85
  cnvgwd = .false.
  TISFC = sfc_fld%TSEA(ipn,1)
  lat       = 1
  nlons     = 200
  xkt2      = 0.6 ! rannum(ipn) ! 0.6
  pre_rad   = .false.
  sinlat    = sin(xlat(ipn,1))
  coslat    = cos(xlat(ipn,1))
  rcs2      = 1.0
  prsshc    = PRSI_S(1,1)
  fhour     = phour
  lssav     = .true.
  solhr     = REAL(mod(REAL(fhour+hour),24.0),kind_phys)
  lsfwd     = .true.
  clstp     = 1110.0
  poz       = 0.0
  prdout    = 0.0
  disout    = 0.0
  flx_fld%PSMEAN(ipn,1) = prsshc
  flx_fld%PSURF(ipn,1)  = prsshc
!TODO:  call flx_init() instead as in GFS r3038 do_physics_one_step.f and 
!TODO:  remove some of these statements...  
  flx_fld%GESHEM(ipn,1) = 0.0
  flx_fld%RAINC(ipn,1)  = 0.0
  cubot                 = 0.0
  cutop                 = 0.0
  if(skip_cu_physics) then
    flx_fld%RAINC(ipn,1)  = phys2dwrf (ipn,6)
    cubot                 = phys2dwrf (ipn,7)
    cutop                 = phys2dwrf (ipn,8)
!       if(flx_fld%RAINC(ipn,1).gt.0)write(6,*)'do_phys',flx_fld%RAINC(ipn,1), &
!                                    cubot,cutop
  endif
  flx_fld%DUSFC(ipn,1)  = 10.0
  flx_fld%DVSFC(ipn,1)  = 10.0
  flx_fld%DTSFC(ipn,1)  = 1.0
  flx_fld%DQSFC(ipn,1)  = 0.0
  flx_fld%GFLUX(ipn,1)  = 0.0
  flx_fld%RUNOFF(ipn,1) = 0.0
  flx_fld%EP(ipn,1)     = 0.0
  flx_fld%CLDWRK(ipn,1) = 0.0
  flx_fld%DUGWD(ipn,1)  = 0.0
  flx_fld%DVGWD(ipn,1)  = 0.0
  flx_fld%BENGSH(ipn,1) = 0.0
  flx_fld%U10M(ipn,1)   = 1.0
  flx_fld%V10M(ipn,1)   = 1.0
  sfc_fld%T2M(ipn,1)    = 300.0
  sfc_fld%Q2M(ipn,1)    = 0.001
  DT3DT     = 0.0
  DQ3DT     = 0.0
  DU3DT     = 0.0
  DV3DT     = 0.0
  LDIAG3D   = .true.
  flipv     = .true.
  EVBSA(:)  = 0.0
  EVCWA(:)  = 0.0
  TRANSA(:) = 0.0
  SBSNOA(:) = 0.0
  SNOWCA(:) = 0.0
  SNOHFA(:) = 0.0
  SPFHMAX(:) = 0.0
  SPFHMIN(:) = 1.e10
  gsoil(:) = 0.0
  gtmp2m(:) = 0.0
  gu10m(:)  =  0.0
  gv10m(:)  =  0.0
  gustar(:) =  0.0
  gzorl(:)  =  0.0
  goro(:) = 0.0
  oro(:) = 0.0
  gpblh(:) = 0.0
  upd_mf(:,:) = 0.0
  dwn_mf(:,:) = 0.0
  det_mf(:,:) = 0.0
  sdiaga(:,:) = 0.0
  sdiagb(:,:) = 0.0
  SRUNOFF(:) = 0.0
   
!----------------------------------------------------------------------
!   Call all other (non-radiation) GFS physics parameterizations
!----------------------------------------------------------------------
  if (ipn.eq.8867) then
!         write(6,*)'water point before call to gbphys'
!         write(6,*)'flx_fld%DLWSFC(ipn,1),flx_fld%ULWSFC(ipn,1)'
!         write(6,*)flx_fld%DLWSFC(ipn,1),flx_fld%ULWSFC(ipn,1)
!         write(6,*)'flx_fld%SFCDLW(ipn,1),flx_fld%TSFLW(ipn,1)'
!         write(6,*)flx_fld%SFCDLW(ipn,1),flx_fld%TSFLW(ipn,1)
  endif

  if (ipn.eq.8570) then
!          write(6,*)'land point before call to gbphys'
!          write(6,*)' flx_fld%HFLX(ipn,1)', flx_fld%HFLX(ipn,1)
!          write(6,*)' flx_fld%EVAP(ipn,1)', flx_fld%EVAP(ipn,1)
!         write(6,*)'flx_fld%DLWSFC(ipn,1),flx_fld%ULWSFC(ipn,1)'
!         write(6,*)flx_fld%DLWSFC(ipn,1),flx_fld%ULWSFC(ipn,1)
!         write(6,*)'flx_fld%SFCDLW(ipn,1),flx_fld%TSFLW(ipn,1)'
!         write(6,*)flx_fld%SFCDLW(ipn,1),flx_fld%TSFLW(ipn,1)
  endif

  call GBPHYS(IM,IX,nvl,lsoil,lsm,ntrac,ncld, &
       ntoz,ntcw,nmtvr,lonf,latg,jcap,ras,nlons,xkt2,nrcm,pre_rad, &
       UG,VG,PGR,TG,QG,gfs_dpdt(ipn,:), &
       GT0,GQ0,GU0,GV0,sinlat,coslat,rcs2,sdiaga,sdiagb, &
       prsi,prsl,prslk,prsik,phii,phil,dpshc,fhour,lssav,solhr, &
       lsfwd,clstp,dtp,dtf,poz,prdout,ko3,pl_coeff, &
       nsst_active,ifd,time_old,time_ins,I_Sw,I_Q,I_Qrain, &
       I_M,I_Tau,I_Sw_Zw,I_Q_Ts,I_M_Ts, &
       Tref,dt_cool,z_c,dt_warm,z_w,c_0,c_d,w_0,w_d, &
       sfc_fld%HICE(ipn,1),sfc_fld%FICE(ipn,1),TISFC,flx_fld%SFCDSW(ipn,1), &           
       sfc_fld%TPRCP(ipn,1), sfc_fld%SRFLAG(ipn,1), &
       sfc_fld%SLC(:,ipn,1),sfc_fld%SNWDPH(ipn,1),sfc_fld%SLOPE(ipn,1),sfc_fld%SHDMIN(ipn,1),sfc_fld%SHDMAX(ipn,1),sfc_fld%SNOALB(ipn,1),SFALB(ipn,1), &
       CHH,CMM,flx_fld%EPI(ipn,1),DLWSFCI,ULWSFCI,USWSFCI,DSWSFCI,DTSFCI, &
       DQSFCI,GFLUXI,SRUNOFF,T1,Q1,U1,V1,ZLVL,EVBSA,EVCWA, &
       TRANSA,SBSNOA,SNOWCA,SOILM,SNOHFA,SMCWLT2,SMCREF2, &
       gsoil,gtmp2m,gustar,gpblh,gu10m,gv10m,gzorl,goro, &
       sfc_fld%TSEA(ipn,1)  ,sfc_fld%SHELEG(ipn,1),sfc_fld%SNCOVR(ipn,1), sfc_fld%TG3(ipn,1),  &
       sfc_fld%ZORL(ipn,1)  ,sfc_fld%CV(ipn,1)    ,sfc_fld%CVB(ipn,1)   ,sfc_fld%CVT(ipn,1)   , &
       sfc_fld%SLMSK(ipn,1) ,sfc_fld%VFRAC(ipn,1) ,sfc_fld%CANOPY(ipn,1),sfc_fld%F10M(ipn,1)  , &
       sfc_fld%VTYPE(ipn,1) ,sfc_fld%STYPE(ipn,1) ,sfc_fld%UUSTAR(ipn,1),sfc_fld%FFMM(ipn,1)  ,sfc_fld%FFHH(ipn,1)  , &
       flx_fld%TMPMIN(ipn,1),flx_fld%TMPMAX(ipn,1), SPFHMIN,SPFHMAX, &
       flx_fld%GESHEM(ipn,1),flx_fld%DUSFC(ipn,1) ,flx_fld%DVSFC(ipn,1) ,flx_fld%DTSFC(ipn,1) , &
       flx_fld%DQSFC(ipn,1) ,flx_fld%DLWSFC(ipn,1),flx_fld%ULWSFC(ipn,1), suntim, &
       flx_fld%GFLUX(ipn,1) ,flx_fld%RUNOFF(ipn,1),flx_fld%EP(ipn,1)    ,flx_fld%CLDWRK(ipn,1), &
       flx_fld%DUGWD(ipn,1) ,flx_fld%DVGWD(ipn,1) ,flx_fld%PSMEAN(ipn,1),flx_fld%RAINC(ipn,1),XLON(ipn,1)  , &
       flx_fld%COSZEN(ipn,1),flx_fld%SFCNSW(ipn,1),XLAT(ipn,1)  , &
       flx_fld%SFCDLW(ipn,1),flx_fld%TSFLW(ipn,1) ,flx_fld%PSURF(ipn,1) ,flx_fld%U10M(ipn,1)  , &
       flx_fld%V10M(ipn,1)  ,sfc_fld%T2M(ipn,1)   ,sfc_fld%Q2M(ipn,1)   , &
       flx_fld%HPBL(ipn,1)  ,flx_fld%PWAT(ipn,1)  ,SWH(ipn,:,1,1),HLW(ipn,:,1,1), &
       sfc_fld%SMC(:,ipn,1),sfc_fld%STC(:,ipn,1),HPRIME(:,ipn,1),slag(ipn,1),sdec(ipn,1),cdec(ipn,1), &
       acv(1),acvb(1),acvt(1), &
       phy_f3d(ipn,:,1,1,:), phy_f2d(ipn,1,:), num_p3d, num_p2d, flgmin, &
       DT3DT, DQ3DT, DU3DT, DV3DT, upd_mf, dwn_mf, det_mf, &
       dkt,dkh,   rnp,  LDIAG3D, lggfs3d, &
       flipv, me,kdt,lat,oro, crtrh, ncw, old_monin,cnvgwd,ccwf,ctei_rm, &
       sashal,newsas,mom4ice,mstrat,trans_trac,cal_pre, &
       flx_fld%HFLX(ipn,1), flx_fld%EVAP(ipn,1), &
       lssav_cc_dummy,DLWSFC_cc_dummy,ULWSFC_cc_dummy,SWSFC_cc_dummy, &
       XMU_cc_dummy, &
       DLW_cc_dummy,DSW_cc_dummy,SNW_cc_dummy,LPREC_cc_dummy, &     
       DUSFC_cc_dummy,DVSFC_cc_dummy,DTSFC_cc_dummy,DQSFC_cc_dummy, &
       PRECR_cc_dummy,skip_cu_physics,skip_mp_physics,cubot,cutop)

  if (ipn.eq.8867) then
!         write(6,*)'water point after call to gbphys'
!         write(6,*)'flx_fld%DLWSFC(ipn,1),flx_fld%ULWSFC(ipn,1)'
!         write(6,*)flx_fld%DLWSFC(ipn,1),flx_fld%ULWSFC(ipn,1)
!         write(6,*)'flx_fld%SFCDLW(ipn,1),flx_fld%TSFLW(ipn,1)'
!         write(6,*)flx_fld%SFCDLW(ipn,1),flx_fld%TSFLW(ipn,1)
  endif

  if (ipn.eq.8570) then
!          write(6,*)'land point after call to gbphys'
!          write(6,*)' flx_fld%HFLX(ipn,1)', flx_fld%HFLX(ipn,1)
!          write(6,*)' flx_fld%EVAP(ipn,1)', flx_fld%EVAP(ipn,1)
!         write(6,*)'land point after call to gbphys'
!          write(6,*)'flx_fld%DLWSFC(ipn,1),flx_fld%ULWSFC(ipn,1)'
!          write(6,*)flx_fld%DLWSFC(ipn,1),flx_fld%ULWSFC(ipn,1)
!          write(6,*)'flx_fld%SFCDLW(ipn,1),flx_fld%TSFLW(ipn,1)'
!          write(6,*)flx_fld%SFCDLW(ipn,1),flx_fld%TSFLW(ipn,1)
  endif

!jbao  old call before Feb 2010
!bao     call GBPHYS(IM,IX,nvl,lsoil,ntrac,ncld,ntoz,ntcw,                &
!bao     nmtvr,lonf,latg,jcap,ras,nlons,xkt2,nrcm,pre_rad,UG,VG,          &
!bao     PGR,TG,QG,gfs_dpdt(ipn,:),GT0,GQ0,GU0,GV0,sinlat,coslat,        &
!bao     rcs2,prsi,prsl,prslk,prsik,phii,phil,prsshc,fhour,lssav,solhr,   &
!bao     lsfwd,clstp,dtp,dtf,poz,prdout,disout,                           &
!bao     sfc_fld%HICE(ipn,1),sfc_fld%FICE(ipn,1),                         &
!bao     flx_fld%SFCDSW(ipn,1),sfc_fld%TPRCP(ipn,1),                      &
!bao     sfc_fld%SRFLAG(ipn,1),sfc_fld%SLC(:,ipn,1),                      &
!bao     sfc_fld%SNWDPH(ipn,1),sfc_fld%SLOPE(ipn,1) ,                     &
!bao     sfc_fld%SHDMIN(ipn,1),sfc_fld%SHDMAX(ipn,1),                     &
!bao     sfc_fld%SNOALB(ipn,1),SFALB(ipn,1),                              &
!bao     sfc_fld%TSEA(ipn,1)  ,sfc_fld%SHELEG(ipn,1),                     &
!bao     sfc_fld%TG3(ipn,1)   ,sfc_fld%ZORL(ipn,1)  ,                     &
!bao     sfc_fld%CV(ipn,1)    ,sfc_fld%CVB(ipn,1)   ,sfc_fld%CVT(ipn,1),  &
!bao     sfc_fld%SLMSK(ipn,1) ,sfc_fld%VFRAC(ipn,1),sfc_fld%CANOPY(ipn,1),&
!bao     sfc_fld%F10M(ipn,1)  ,sfc_fld%VTYPE(ipn,1) ,                     &
!bao     sfc_fld%STYPE(ipn,1) ,sfc_fld%UUSTAR(ipn,1),                     &
!bao     sfc_fld%FFMM(ipn,1)  ,sfc_fld%FFHH(ipn,1)  ,                     &
!bao     flx_fld%TMPMIN(ipn,1),flx_fld%TMPMAX(ipn,1),                     &
!bao     flx_fld%GESHEM(ipn,1),flx_fld%DUSFC(ipn,1) ,                     &
!bao     flx_fld%DVSFC(ipn,1) ,flx_fld%DTSFC(ipn,1) ,                     &
!bao     flx_fld%DQSFC(ipn,1) ,flx_fld%DLWSFC(ipn,1),                     &
!bao     flx_fld%ULWSFC(ipn,1),flx_fld%GFLUX(ipn,1) ,                     &
!bao     flx_fld%RUNOFF(ipn,1),flx_fld%EP(ipn,1)    ,                     &
!bao     flx_fld%CLDWRK(ipn,1),flx_fld%DUGWD(ipn,1) ,                     &
!bao     flx_fld%DVGWD(ipn,1) ,flx_fld%PSMEAN(ipn,1),                     &
!bao     flx_fld%BENGSH(ipn,1),XLON(ipn,1),flx_fld%COSZEN(ipn,1),         &
!bao     flx_fld%SFCNSW(ipn,1),XLAT(ipn,1),flx_fld%SFCDLW(ipn,1),         &
!bao     flx_fld%TSFLW(ipn,1),flx_fld%PSURF(ipn,1),                       &
!bao     flx_fld%U10M(ipn,1)  ,flx_fld%V10M(ipn,1)  ,                     &
!bao     sfc_fld%T2M(ipn,1)   ,sfc_fld%Q2M(ipn,1)   ,                     &
!bao     flx_fld%HPBL(ipn,1)  ,flx_fld%PWAT(ipn,1)  ,                     &
!bao     SWH(ipn,:,1,1),HLW(ipn,:,1,1),sfc_fld%SMC(:,ipn,1),              &
!bao     sfc_fld%STC(:,ipn,1),HPRIME(:,ipn,1),slag(ipn,1),                &
!bao     sdec(ipn,1),cdec(ipn,1),acv(1),                                  &
!bao     acvb(1),acvt(1),phy_f3d(ipn,:,1,1,:),phy_f2d(ipn,1,:),           &
!bao     num_p3d,num_p2d,DT3DT, DQ3DT, DU3DT, DV3DT, LDIAG3D,flipv,       &
!bao     me,kdt,1,flx_fld%HFLX(ipn,1), flx_fld%EVAP(ipn,1),               &
!bao     flx_fld%RAINC(ipn,1),GravityWaveDrag                  )
  if (ipn.eq.8570) then
!          write(6,*)'land point before gbphys t and moist'
!          write(6,*)'u',gfs_u(ipn,:)
!          write(6,*)'v',gfs_v(ipn,:)
!          write(6,*)'t',gfs_t(ipn,:)
!          write(6,*)'q',gfs_q(ipn,:)
!          write(6,*)'cld',gfs_cld(ipn,:)
  endif

  if ( .not. skip_cu_physics  .or. .not. skip_mp_physics) then
!     if(ipn.eq.3071)then
!        do ivl=1,nvl-1
!          write(6,*)'gq0 = ',gq0(1,ivl,1)
!        enddo
!     endif
    do ivl=1,nvl-1
!     diaga(ivl,ipn)    = sdiaga(1,ivl)
!     diaga(ivl,ipn)    = dkt(1,ivl)
      diagb(ivl,ipn)    = sdiagb(1,ivl)
    enddo
  endif

! Set tendency arrays from physics

  do ivl=1,nvl
    u_tdcy_phy(ivl,ipn) = (gu0(1,ivl)-gfs_u(ipn,ivl))/dtp
    v_tdcy_phy(ivl,ipn) = (gv0(1,ivl)-gfs_u(ipn,ivl))/dtp
    trc_tdcy_phy(ivl,ipn,1) = (gt0(1,ivl)-gfs_t(ipn,ivl))/dtp
    trc_tdcy_phy(ivl,ipn,2) = (gq0(1,ivl,1)-gfs_q(ipn,ivl))/dtp
    trc_tdcy_phy(ivl,ipn,3) = (gq0(1,ivl,3)-gfs_cld(ipn,ivl))/dtp
    trc_tdcy_phy(ivl,ipn,4) = (gq0(1,ivl,2)-gfs_oz(ipn,ivl))/dtp
    gfs_u(ipn,ivl)   = gu0(1,ivl)
    gfs_v(ipn,ivl)   = gv0(1,ivl)
    gfs_t(ipn,ivl)   = gt0(1,ivl)
    gfs_q(ipn,ivl)   = gq0(1,ivl,1)
    gfs_oz(ipn,ivl)  = gq0(1,ivl,2)
    gfs_cld(ipn,ivl) = gq0(1,ivl,3)
  enddo
  if (ipn.eq.8570) then
!          write(6,*)'land point after gbphys t and moist'
!          write(6,*)'u',gfs_u(ipn,:)
!          write(6,*)'v',gfs_v(ipn,:)
!          write(6,*)'t',gfs_t(ipn,:)
!          write(6,*)'q',gfs_q(ipn,:)
!          write(6,*)'cld',gfs_cld(ipn,:)
  endif

  if (chem_opt.gt.0 .or. skip_cu_physics  .or. skip_mp_physics) then
    call do_physics_one_step_chem(ipn,skip_chem,skip_cu_physics,skip_mp_physics,dkt,dq3dt,dt3dt)
  endif
enddo !ipn
!jbao     stop

! In the next loop, set temp arrays to new values output from physics

do ivl=1,nvl
  do ipn=1,nip
    temp_ps(ipn) = gfs_ps(ipn)
    temp_dp(ipn,ivl) = gfs_dp(ipn,ivl)
    temp_dpdt(ipn,ivl) = gfs_dpdt(ipn,ivl)
    temp_p(ipn,ivl) = gfs_p(ipn,ivl)
    temp_u(ipn,ivl) = gfs_u(ipn,ivl)
    temp_v(ipn,ivl) = gfs_v(ipn,ivl)
    temp_t(ipn,ivl) = gfs_t(ipn,ivl)
    temp_q(ipn,ivl) = gfs_q(ipn,ivl)
    temp_oz(ipn,ivl) = gfs_oz(ipn,ivl)
    temp_cld(ipn,ivl) =  gfs_cld(ipn,ivl)
    if(ivl.eq.1.and.ipn.eq.1) then 
      call PhysicsGetIpnItsMype(ipnGlobal,its,mype,DiagPrint)
!           print"('do_physics',i5,i9,i5)",mype,ipnGlobal,ivl
    endif
  enddo
enddo
!sms$compare_var(temp_ps, "do_physics.F90 - gfs_ps1 ")
!sms$compare_var(temp_dp, "do_physics.F90 - gfs_dp1 ")
!sms$compare_var(temp_dpdt, "do_physics.F90 - gfs_dpdt1 ")
!sms$compare_var(temp_p, "do_physics.F90 - gfs_p1 ")
!sms$compare_var(temp_oz, "do_physics.F90 - gfs_oz1 ")
!sms$compare_var(temp_cld, "do_physics.F90 - gfs_cld1 ")
!sms$compare_var(temp_u, "do_physics.F90 - gfs_u1 ")
!sms$compare_var(temp_v, "do_physics.F90 - gfs_v1 ")
!sms$compare_var(temp_t, "do_physics.F90 - gfs_t1 ")
!sms$compare_var(temp_q, "do_physics.F90 - gfs_q1 ")
!call PhysicsGetIpnItsMype(ipnGlobal,its,mype,DiagPrint)
!SMS$PARALLEL END

return
end subroutine do_physics_one_step

end module module_do_physics_one_step
