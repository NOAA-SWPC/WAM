       subroutine gloopr
     x    (ncld,
     x     lats_nodes_r,global_lats_r,
     x     lonsperlar,
!jbao not needed by fim     x     epse,epso,epsedn,epsodn,
!jbao not needed by fim     x     snnp1ev,snnp1od, plnev_r,plnod_r,
!jbao not needed by fim     x     pddev_r,pddod_r,
!    x     snnp1ev,snnp1od,ndexev,ndexod,
!    x     plnev_r,plnod_r,pddev_r,pddod_r,plnew_r,plnow_r,
     x     phour,
     &     xlon,xlat,coszdg,COSZEN,
     &     SLMSK,SHELEG,SNCOVR,SNOALB,ZORL,TSEA,HPRIME,
Clu [+1L]: extract snow-free albedo (SFALB)
     +     SFALB,
     &     ALVSF,ALNSF ,ALVWF ,ALNWF,FACSF ,FACWF,CV,CVT ,
     &     CVB  ,SWH,HLW,SFCNSW,SFCDLW,
     &     FICE ,TISFC, SFCDSW,                          ! FOR SEA-ICE - XW Nov04
     &     TSFLW,FLUXR ,       phy_f3d,slag,sdec,cdec,NBLCK,KDT,
     &     global_times_r,prsl,prsi,prslk,gt,gr,gr1,
     &     sscal,asymp,ext_cof,extlw_cof,yyyymmddhhmm)
cc
!jbao new gfs#include "f_hpm.h"
!
      USE MACHINE              ,     ONLY : kind_phys
      USE FUNCPHYS             ,     ONLY : fpkap
      USE PHYSCONS, FV => con_fvirt, rerth => con_rerth 	! hmhj

      use module_radiation_driver,   only : radinit, grrad
      use module_radiation_astronomy,only : astronomy
!
!! ---  for optional spectral band heating outputs
!!    use module_radsw_parameters,   only : NBDSW
!!    use module_radlw_parameters,   only : NBDLW
!
      use resol_def
      use layout1
      use gg_def
      use vert_def
      use date_def
      use namelist_def
      use coordinate_def					! hmhj
      use tracer_const						! hmhj
! jbao old gloopr cldcov passed in, we don't have d3d_def in fim??
!      use d3d_def , only : cldcov
!
      implicit none
! jbao does fim need this???
      include 'mpif.h'
!
! jbao new gfs add ncld , grrad needs it
       integer, intent(in) :: ncld
! jbao newgfs add ncld , grrad needs it

      real (kind=kind_phys), parameter :: QMIN =1.0e-10
      real (kind=kind_evod), parameter :: Typical_pgr = 95.0
      real (kind=kind_evod), parameter :: cons0 = 0.0,  cons2 = 2.0
!
!  --- ...  inputs:
      integer                ls_node, ls_nodes,  max_ls_nodes
      integer, intent(in) :: lats_nodes_r,                              &
     &                       global_lats_r(LATR), lonsperlar(LATR)

      integer, intent(in) :: NBLCK

!jbao not needed by fim      real(kind=kind_evod), dimension(LEN_TRIE_LS), intent(in) ::       &
!jbao not needed by fim     &                       epse, epsedn, snnp1ev

!jbao not needed by fim      real(kind=kind_evod), dimension(LEN_TRIO_LS), intent(in) ::       &
!jbao not needed by fim     &                       epso, epsodn, snnp1od

!jbao not needed by fim      real(kind=kind_evod), intent(in) :: plnev_r(LEN_TRIE_LS, LATR2)
!jbao not needed by fim      real(kind=kind_evod), intent(in) :: plnod_r(LEN_TRIO_LS, LATR2)

      real (kind=kind_phys), dimension(LONR,LATS_NODE_R), intent(in) :: &
     &                       xlon, xlat, slmsk, sheleg, zorl, tsea,     &
     &                       alvsf, alnsf, alvwf, alnwf, facsf, facwf,  &
     &                       cv, cvt, cvb, FICE, tisfc, sncovr, snoalb

      real (kind=kind_phys), intent(in) ::                              &
! jbao orig, in old fim newgfs hprime is defined as 1,lonr..     &                    hprime(NMTVR,LONR,LATS_NODE_R), phour,        &
     &                     HPRIME(    1,lonr,lats_node_r), phour,        &
     &                    phy_f3d(NGPTC,LEVS,NBLCK,LATS_NODE_R,NUM_P3D)
!
!  --- ...  input and output:
!jbao      real(kind=kind_evod), intent(inout) ::                            &
!jbao     &                    trie_ls(LEN_TRIE_LS,2,11*LEVS+3*LEVH+6),      &
!jbao     &                    trio_ls(LEN_TRIO_LS,2,11*LEVS+3*LEVH+6)
      integer              ipt_ls                                       ! hmhj
      real(kind=kind_evod) reall                                        ! hmhj
      real(kind=kind_evod) rlcs2(jcap1)                                 ! hmhj


      real (kind=kind_phys), intent(inout) ::                           &
     &                    fluxr (33,LONR,LATS_NODE_R)
! jbao orig ,nfxr isn't defined til later but=27  and added cldcov not 27, but levs  &                    fluxr (NFXR,LONR,LATS_NODE_R)
      real (kind=kind_phys) CLDCOV(LEVS,lonr,lats_node_r)               !, intent(out) ::                             &
!    &                    cldcov (27,LONR,LATS_NODE_R)
! jbao end new gfs define cldcov

!  --- ...  inputs but not used anymore:
!jbao not needed by fim      real(kind=kind_evod), intent(in) :: pddev_r(LEN_TRIE_LS,LATR2),   &
!jbao not needed by fim     &                                    pddod_r(LEN_TRIO_LS,LATR2)    &
!    &                                    plnew_r(LEN_TRIE_LS,LATR2),   &
!    &                                    plnow_r(LEN_TRIO_LS,LATR2)
!    &                                    syn_ls_r(4*LS_DIM,LOTS,LATR2)

!     integer, intent(in) :: ndexev(LEN_TRIE_LS), ndexod(LEN_TRIO_LS)
      integer, intent(in) :: KDT
!  --- ...  outputs:
! jbao old gloopr had global_times_r defined as the following:
             real(kind=kind_evod) global_times_r
! jbao orig       real(kind=kind_evod), intent(out) ::                              &
! jbao orig in old glooper get rid of dimesnions     &                    global_times_r(LATG,NODES)
      real(kind=kind_evod) ::                                           &
     &                    for_gr_r_1(LONRX*LOTS,LATS_DIM_R),            &
     &                    dyn_gr_r_1(lonrx*lotd,lats_dim_r),            ! hmhj
     &                    for_gr_r_2(LONRX*LOTS,LATS_DIM_R),
     &                    dyn_gr_r_2(lonrx*lotd,lats_dim_r)             ! hmhj

      real (kind=kind_phys), intent(out) ::                             &
     &                    swh(NGPTC,LEVS,NBLCK,LATS_NODE_R),            &
     &                    hlw(NGPTC,LEVS,NBLCK,LATS_NODE_R)

      real (kind=kind_phys),dimension(LONR,LATS_NODE_R), intent(out) :: &
     &                    coszdg, coszen, sfcnsw, sfcdlw, tsflw,        &
     &                    sfcdsw, SFALB

      real (kind=kind_phys), intent(out) :: slag, sdec, cdec

!! --- ...  optional spectral band heating rates
!!    real (kind=kind_phys), optional, intent(out) ::                   &
!!   &                 htrswb(NGPTC,LEVS,NBDSW,NBLCK,LATS_NODE_R),      &
!!   &                 htrlwb(NGPTC,LEVS,NBDLW,NBLCK,LATS_NODE_R)

!  --- ...  locals:
!     real(kind=kind_phys) :: prsl(NGPTC,LEVS), prdel(NGPTC,LEVS),      &
      real(kind=kind_phys) :: prsl(NGPTC,LEVS),  prslk(NGPTC,LEVS),     &
     &                        prsi(NGPTC,LEVP1), prsik(NGPTC,LEVP1)

! jbao newgfs change levr to levs 
      real (kind=kind_phys) :: si_loc(levs+1)
!     real (kind=kind_phys) :: si_loc(levs+1), prslk(NGPTC,levs)

      real (kind=kind_phys) :: gu(NGPTC,LEVS), gv1(NGPTC,LEVS),         &
     &                        gt(NGPTC,levs), gd (NGPTC,LEVS),          &
! jbao orig old has ngptc,levs,2 ntrac -1     &                        gr(NGPTC,levs), gr1(NGPTC,levs,NTRAC-1),  &
     &                        gr(NGPTC,levs), gr1(NGPTC,levs,2),        &
     &                        gphi(NGPTC), glam(NGPTC), gq(NGPTC),      &
     &                        sumq(NGPTC,levs), xcp(NGPTC,levs),        &! hmhj	
     &                        gtv(NGPTC,levs), gtvx(NGPTC,levs),        &! hmhj	
     &                        gtvy(NGPTC,levs)				 ! hmhj

      real (kind=kind_phys) :: f_ice(NGPTC,LEVS), f_rain(NGPTC,LEVS),   &
     &                        r_rime(NGPTC,LEVS)

      real (kind=kind_phys) :: cldcov_v(NGPTC,LEVS), hprime_v(NGPTC),   &
     &                        fluxr_v(NGPTC,33), vvel(NGPTC,LEVS)
! jbao old glooper nfxr not defined yet is 27     &                        fluxr_v(NGPTC,NFXR), vvel(NGPTC,LEVS)
      real (kind=kind_phys) :: flgmin_l(ngptc), work1, work2

      real (kind=kind_phys) :: rinc(5), dtsw, dtlw, solcon

      real (kind=kind_phys), save :: facoz
!-----aerosols Bao-----------------------------------------
       real (kind=kind_phys)  sscal(LEVS,14)
       real (kind=kind_phys)  asymp(LEVS,14)
       real (kind=kind_phys)  ext_cof(LEVS,14)
       real (kind=kind_phys)  extlw_cof(LEVS,16)
!-----end aerosols Bao-----------------------------------------

      integer :: njeff, lon, lan, lat, iblk, lon_dim, lons_lat, istrt
      integer :: idat(8), jdat(8), DAYS(13), iday, imon, midmon, id
      integer :: lmax
      CHARACTER(len=12) :: yyyymmddhhmm
      INTEGER year, month, day, hour, minute

      integer, save :: icwp, k1oz, k2oz, midm, midp

!  ---  number of days in a month
      data DAYS / 31,28,31,30,31,30,31,31,30,31,30,31,30 /

!  --- ...  control parameters: 
!           (some of the them may be moved into model namelist)

!  ---  ICTM=yyyy#, controls time sensitive external data (e.g. CO2, solcon, aerosols, etc)
!     integer, parameter :: ICTM =    0 ! use data at initial cond time, if not
!                                       ! available, use latest, no extrapolation.
!!    integer, parameter :: ICTM =    1 ! use data at the forecast time, if not
!                                       ! available, use latest and extrapolation.
!     integer, parameter :: ICTM =yyyy0 ! use yyyy data for the forecast time,
!                                       ! no further data extrapolation.
!     integer, parameter :: ICTM =yyyy1 ! use yyyy data for the fcst. if needed, do
!                                       ! extrapolation to match the fcst time.

!  ---  ISOL controls solar constant data source
!!    integer, parameter :: ISOL = 0   ! use prescribed solar constant
!     integer, parameter :: ISOL = 1   ! use varying solar const with 11-yr cycle

!  ---  ICO2 controls co2 data source for radiation
!     integer, parameter :: ICO2 = 0   ! prescribed global mean value (old opernl)
!!    integer, parameter :: ICO2 = 1   ! use obs co2 annual mean value only
!     integer, parameter :: ICO2 = 2   ! use obs co2 monthly data with 2-d variation

!  ---  IALB controls surface albedo for sw radiation
!!    integer, parameter :: IALB = 0   ! use climatology alb, based on sfc type
!     integer, parameter :: IALB = 1   ! use modis derived alb (to be developed)

!  ---  IEMS controls surface emissivity for lw radiation
!!    integer, parameter :: IEMS = 0   ! use fixed value of 1.0
!     integer, parameter :: IEMS = 1   ! use varying sfc emiss, based on sfc type
!  ---  IAER  controls aerosols scheme selections
!     Old definition
!     integer, parameter :: IAER  = 1  ! opac climatology, without volc forcing
!     integer, parameter :: IAER  =11  ! opac climatology, with volcanic forcing
!     integer, parameter :: IAER  = 2  ! gocart prognostic, without volc forcing
!     integer, parameter :: IAER  =12  ! gocart prognostic, with volcanic forcing
!     New definition in this code
!  IAER =   0 --> no aerosol effect at all (volc, sw, lw)
!       =   1 --> only tropospheric sw aerosols, no trop-lw and volc
!       =  10 --> only tropospheric lw aerosols, no trop-sw and volc
!       =  11 --> both trop-sw and trop-lw aerosols, no volc
!       = 100 --> only strato-volc aeros, no trop-sw and trop-lw
!       = 101 --> only sw aeros (trop + volc), no lw aeros
!       = 110 --> only lw aeros (trop + volc), no sw aeros
!       = 111 --> both sw and lw aeros (trop + volc) 
!

!  ---  IOVR controls cloud overlapping method in radiation:
!     integer, parameter :: IOVR_SW = 0  ! sw: random overlap clouds
!!    integer, parameter :: IOVR_SW = 1  ! sw: max-random overlap clouds

!     integer, parameter :: IOVR_LW = 0  ! lw: random overlap clouds
!!    integer, parameter :: IOVR_LW = 1  ! lw: max-random overlap clouds

!  ---  iflip indicates model vertical index direction:
!     integer, parameter :: IFLIP = 0    ! virtical profile index from top to bottom
      integer, parameter :: IFLIP = 1    ! virtical profile index from bottom to top
!
!    The following parameters are from gbphys
!
      real (kind=kind_phys), parameter :: dxmax=-16.118095651,          &
     &                dxmin=-9.800790154, dxinv=1.0/(dxmax-dxmin)

      integer :: kr, kt, kd, kq, ku, kv, ierr, dimg, kx, ky
      integer :: i, j, k, n
      integer :: kdtphi,kdtlam,ks                                ! hmhj

      logical :: lslag, change, lprnt
      data  lslag / .false. /,    lprnt / .false. /
      logical, save :: first, sas_shal
      data  first / .true. /

!  ---  timers:
      real*8 :: rtc, timer1, timer2
! jbao new variable for the update as of Feb 2010
      integer nfxr

!
!===> *** ...  begin here
!
!!
cc
!cmr  ls_node(1,1) ... ls_node(ls_max_node,1) : values of L
!cmr  ls_node(1,2) ... ls_node(ls_max_node,2) : values of jbasev
!cmr  ls_node(1,3) ... ls_node(ls_max_node,3) : values of jbasod
cc
cc
c$$$      integer                lots,lotd,lota
c$$$cc
c$$$      parameter            ( lots = 5*levs+1*levh+3 )
c$$$      parameter            ( lotd = 6*levs+2*levh+0 )
c$$$      parameter            ( lota = 3*levs+1*levh+1 )
cc
cc
      integer              kap,kar,kat,kau,kav,kdrlam
      integer              ksd,ksplam,kspphi,ksq,ksr,kst
      integer              ksu,ksv,ksz,node
cc
!     real(kind=kind_evod) spdlat(levs,lats_dim_r)
!Moor real(kind=kind_phys) slk(levs)
!     real(kind=kind_evod) spdmax_node (levs)
!     real(kind=kind_evod) spdmax_nodes(levs,nodes)
cc
cc
ccxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
cc
cc
cc................................................................
cc  syn(1, 0*levs+0*levh+1, lan)  ze
cc  syn(1, 1*levs+0*levh+1, lan)  di
cc  syn(1, 2*levs+0*levh+1, lan)  te
cc  syn(1, 3*levs+0*levh+1, lan)  rq
cc  syn(1, 3*levs+1*levh+1, lan)  q
cc  syn(1, 3*levs+1*levh+2, lan)  dpdlam
cc  syn(1, 3*levs+1*levh+3, lan)  dpdphi
cc  syn(1, 3*levs+1*levh+4, lan)  uln
cc  syn(1, 4*levs+1*levh+4, lan)  vln
cc................................................................
cc  dyn(1, 0*levs+0*levh+1, lan)  d(t)/d(phi)
cc  dyn(1, 1*levs+0*levh+1, lan)  d(rq)/d(phi)
cc  dyn(1, 1*levs+1*levh+1, lan)  d(t)/d(lam)
cc  dyn(1, 2*levs+1*levh+1, lan)  d(rq)/d(lam)
cc  dyn(1, 2*levs+2*levh+1, lan)  d(u)/d(lam)
cc  dyn(1, 3*levs+2*levh+1, lan)  d(v)/d(lam)
cc  dyn(1, 4*levs+2*levh+1, lan)  d(u)/d(phi)
cc  dyn(1, 5*levs+2*levh+1, lan)  d(v)/d(phi)
cc................................................................
cc  anl(1, 0*levs+0*levh+1, lan)  w     dudt
cc  anl(1, 1*levs+0*levh+1, lan)  x     dvdt
cc  anl(1, 2*levs+0*levh+1, lan)  y     dtdt
cc  anl(1, 3*levs+0*levh+1, lan)  rt    drdt
cc  anl(1, 3*levs+1*levh+1, lan)  z     dqdt
cc................................................................
cc
cc
c$$$      parameter(ksz     =0*levs+0*levh+1,
c$$$     x          ksd     =1*levs+0*levh+1,
c$$$     x          kst     =2*levs+0*levh+1,
c$$$     x          ksr     =3*levs+0*levh+1,
c$$$     x          ksq     =3*levs+1*levh+1,
c$$$     x          ksplam  =3*levs+1*levh+2,
c$$$     x          kspphi  =3*levs+1*levh+3,
c$$$     x          ksu     =3*levs+1*levh+4,
c$$$     x          ksv     =4*levs+1*levh+4)
cc
c$$$      parameter(kdtphi  =0*levs+0*levh+1,
c$$$     x          kdrphi  =1*levs+0*levh+1,
c$$$     x          kdtlam  =1*levs+1*levh+1,
c$$$     x          kdrlam  =2*levs+1*levh+1,
c$$$     x          kdulam  =2*levs+2*levh+1,
c$$$     x          kdvlam  =3*levs+2*levh+1,
c$$$     x          kduphi  =4*levs+2*levh+1,
c$$$     x          kdvphi  =5*levs+2*levh+1)
cc
c$$$      parameter(kau     =0*levs+0*levh+1,
c$$$     x          kav     =1*levs+0*levh+1,
c$$$     x          kat     =2*levs+0*levh+1,
c$$$     x          kar     =3*levs+0*levh+1,
c$$$     x          kap     =3*levs+1*levh+1)
cc
cc
c$$$      integer   P_gz,P_zem,P_dim,P_tem,P_rm,P_qm
c$$$      integer   P_ze,P_di,P_te,P_rq,P_q,P_dlam,P_dphi,P_uln,P_vln
c$$$      integer   P_w,P_x,P_y,P_rt,P_zq
c$$$cc
c$$$cc                                               old common /comfspec/
c$$$      parameter(P_gz  = 0*levs+0*levh+1,  !      gze/o(lnte/od,2),
c$$$     x          P_zem = 0*levs+0*levh+2,  !     zeme/o(lnte/od,2,levs),
c$$$     x          P_dim = 1*levs+0*levh+2,  !     dime/o(lnte/od,2,levs),
c$$$     x          P_tem = 2*levs+0*levh+2,  !     teme/o(lnte/od,2,levs),
c$$$     x          P_rm  = 3*levs+0*levh+2,  !      rme/o(lnte/od,2,levh),
c$$$     x          P_qm  = 3*levs+1*levh+2,  !      qme/o(lnte/od,2),
c$$$     x          P_ze  = 3*levs+1*levh+3,  !      zee/o(lnte/od,2,levs),
c$$$     x          P_di  = 4*levs+1*levh+3,  !      die/o(lnte/od,2,levs),
c$$$     x          P_te  = 5*levs+1*levh+3,  !      tee/o(lnte/od,2,levs),
c$$$     x          P_rq  = 6*levs+1*levh+3,  !      rqe/o(lnte/od,2,levh),
c$$$     x          P_q   = 6*levs+2*levh+3,  !       qe/o(lnte/od,2),
c$$$     x          P_dlam= 6*levs+2*levh+4,  !  dpdlame/o(lnte/od,2),
c$$$     x          P_dphi= 6*levs+2*levh+5,  !  dpdphie/o(lnte/od,2),
c$$$     x          P_uln = 6*levs+2*levh+6,  !     ulne/o(lnte/od,2,levs),
c$$$     x          P_vln = 7*levs+2*levh+6,  !     vlne/o(lnte/od,2,levs),
c$$$     x          P_w   = 8*levs+2*levh+6,  !       we/o(lnte/od,2,levs),
c$$$     x          P_x   = 9*levs+2*levh+6,  !       xe/o(lnte/od,2,levs),
c$$$     x          P_y   =10*levs+2*levh+6,  !       ye/o(lnte/od,2,levs),
c$$$     x          P_rt  =11*levs+2*levh+6,  !      rte/o(lnte/od,2,levh),
c$$$     x          P_zq  =11*levs+3*levh+6)  !      zqe/o(lnte/od,2)
cc
cc
!     print *,' in gloopr vertcoord_id =',vertcoord_id
!=================================================================================
!!!   ********* temporary return - Stan B - 14 Aug 2010

!     return
!!!   
!=================================================================================


! jbao new variable for the update as of Feb 2010
      ls_node = 1
      ls_nodes = 1
      max_ls_nodes = 1
      f_ice = 0.0
      f_rain = 0.0
      r_rime =0.0
      cldcov =0.0

!  some of these are declared--jbao do we need to declare them?
      num_p3d=4
      nfxr=33
      sashal=.true. ! jbao true from gfs namelist
      norad_precip=.false. ! jbao false from gfs namelist
      crick_proof=.false. ! jbao false from gfs namelist
      ccnorm=.false. ! jbao false from gfs namelist
      lggfs3d=.false.
      lprnt=.false.
      ras       = .false.  ! jbao needed by new gfs from do onestep

      ntrac = 3 ! new gfs physics 3 !JFM and Bao originally in new gloopr was 2
      ntcw = 3 !JFM and Bao
      ntoz = 2 !JFM and Bao
      vvel = 0.0 !JFM and Bao
      LDIAG3D = .false. !JFM
      lssav = .true. !JFM false means FLUXR is not calculated in grrad, ok because fluxr not used in physics.F90
      fluxr = 0.0
      fluxr_v = 0.0
! jbao end new gfs

!
! jbao not needed by fim      ksz     =0*levs+0*levh+1
! jbao not needed by fim      ksd     =1*levs+0*levh+1
! jbao not needed by fim      kst     =2*levs+0*levh+1
! jbao not needed by fim      ksr     =3*levs+0*levh+1
! jbao not needed by fim      ksq     =3*levs+1*levh+1
! jbao not needed by fim      ksplam  =3*levs+1*levh+2
! jbao not needed by fim      kspphi  =3*levs+1*levh+3
! jbao not needed by fim      ksu     =3*levs+1*levh+4
! jbao not needed by fim      ksv     =4*levs+1*levh+4

! jbao not needed by fim      kdtphi  =0*levs+0*levh+1                          ! hmhj
! jbao not needed by fim      kdtlam  =1*levs+1*levh+1                          ! hmhj
cc
! jbao not needed by fim      lslag=.false.
cc
      idat = 0
! get date info from the date string
      READ(UNIT=yyyymmddhhmm(1:4), FMT='(I4)') year
      READ(UNIT=yyyymmddhhmm(5:6), FMT='(I2)') month
      READ(UNIT=yyyymmddhhmm(7:8), FMT='(I2)') day
      READ(UNIT=yyyymmddhhmm(9:10), FMT='(I2)') hour
      READ(UNIT=yyyymmddhhmm(11:12), FMT='(I2)') minute
      idat(1) = year
      idat(2) = month
      idat(3) = day
      idat(5) = hour

      rinc = 0.
      rinc(2) = phour
      call w3movdat(rinc, idat, jdat)  
!
      if (ntoz .le. 0) then                ! Climatological Ozone!
!
!     if(me .eq. 0) WRITE (6,989) jdat(1),jdat(2),jdat(3),jdat(5)
! 989 FORMAT(' UPDATING OZONE FOR ', I4,I3,I3,I3)
!
        IDAY   = jdat(3)
        IMON   = jdat(2)
        MIDMON = DAYS(IMON)/2 + 1
        CHANGE = FIRST .OR.
     &          ( (IDAY .EQ. MIDMON) .AND. (jdat(5).EQ.0) )
!
        IF (CHANGE) THEN
            IF (IDAY .LT. MIDMON) THEN
               K1OZ = MOD(IMON+10,12) + 1
               MIDM = DAYS(K1OZ)/2 + 1
               K2OZ = IMON
               MIDP = DAYS(K1OZ) + MIDMON
            ELSE
               K1OZ = IMON
               MIDM = MIDMON
               K2OZ = MOD(IMON,12) + 1
               MIDP = DAYS(K2OZ)/2 + 1 + DAYS(K1OZ)
            ENDIF
        ENDIF
!
        IF (IDAY .LT. MIDMON) THEN
           ID = IDAY + DAYS(K1OZ)
        ELSE
           ID = IDAY
        ENDIF
        FACOZ = real (ID-MIDM) / real (MIDP-MIDM)
      endif
! jbao do we need to do this or as in old gloopr goto 11111??
! jbao need ras in here
!
      if (first) then
        sas_shal = sashal .and. (.not. ras)
      goto 1111  ! jbao skip find levels
!
! jbao do we need to bring in these variables? they are defined elsewhere? hybrid not defined 
        if( hybrid.or.gen_coord_hybrid ) then                             ! hmhj

          if( gen_coord_hybrid ) then                                     ! hmhj
            si_loc(levs+1) = si(levp1)                                    ! hmhj
            do k=1,levs                                                   ! hmhj
              si_loc(k) = si(k)                                           ! hmhj
            enddo                                                         ! hmhj
          else                                                            ! hmhj
!  ---  get some sigma distribution for radiation-cloud initialization
!sela   si(k)=(ak5(k)+bk5(k)*Typical_pgr)/Typical_pgr   !ak(k) bk(k) go top to botto
            si_loc(levs+1)= ak5(1)/typical_pgr+bk5(1)
            do k=1,levs
              si_loc(levs+1-k)= ak5(levp1-levs+k)/typical_pgr
     &                        + bk5(levp1-levs+k)
            enddo
          endif
! jbao it will do this comment out??
        else
          do k = 1, levs
            si_loc(k) = si(k)
          enddo
          si_loc(levs+1) = si(levp1)
        endif       ! end_if_hybrid

!  --- determin prognostic/diagnostic cloud scheme
! jbao new gfs continue- here??
1111    continue
        icwp   = 0
        if (NTCW > 0) icwp = 1
           
! jbao new gfs: old gloopr first=false isn't done 111 continue is after this but icwp =0 is done
        first = .false.
           
      endif         ! end_if_first
! jbao end do we need to do this
! jbao new gfs needs sinlat_r, coslat_r, in previous gloopr in fim, did that here
         allocate (sinlat_r(1),coslat_r(1))
         sinlat_r(1) = sin(xlat(1,1)) ! jbao
         coslat_r(1) = cos(xlat(1,1)) ! jbao
         lsswr = .true. ! jbao
         lslwr = .true. ! jbao
         fhswr = 1.0 ! jbao test 0 cause floating point
         fhlwr = 1.0 ! jbao test 0 cause floating point
         isol = 0 ! jbao
         ico2 = 1  ! jbao
!   0 - use default value
!   1 - read in CO2 file
         ialb = 0  ! jbao
         iaer = 1  ! 0 ! jbao test to turn off aerosols
         iovr_lw=1
         iovr_sw=1
         iems = 0 ! jbao
         ictm = 1 ! jbao
         do k=1,levp1
           SI_loc(k)= prsi(1,k)/prsi(1,1)
         enddo


! end jbao new gfs needs sinlat_r, coslat_r, other variables above in previous gloopr in fim, did that here

!
!===> *** ...  radiation initialization
!
      dtsw  = 3600.0 * fhswr
      dtlw  = 3600.0 * fhlwr
                                                                                                            
! jbao is si_loc correct? what about ico2?, isol?,ialb,iaer (before, passed in numbers),ictm not defined
      call radinit                                                      &
!  ---  input:
     &     ( si_loc, levs, IFLIP, NUM_P3D, ICTM,                        &
     &       ISOL, ICO2, ICWP, IALB, IEMS, IAER, idat, jdat, me )
!  ---  output: ( none )
                                                                                                            
!
!===> *** ...  astronomy for sw radiation calculation.
!
      call astronomy                                                    &
!  ---  inputs:
     &     ( lonsperlar, global_lats_r, sinlat_r, coslat_r, xlon,       &
!    &       fhswr, jdat, deltim,                                       &
     &       fhswr, jdat,                                               &
     &       LONR, LATS_NODE_R, LATR, IPT_LATS_NODE_R, lsswr, me,       &
!  ---  outputs:
     &       solcon, slag, sdec, cdec, coszen, coszdg                   &
     &      )
                                                                                                            
!     print *,' returned from astro'
!
!===> *** ...  spectrum to grid transformation for radiation calculation.
!     -----------------------------------
cc
! jbao new gfs not needed by fim      call f_hpmstart(61,"gr delnpe")
! jbao new gfs not needed by fim      call delnpe(trie_ls(1,1,P_q   ),
! jbao new gfs not needed by fim     x            trio_ls(1,1,P_dphi),
! jbao new gfs not needed by fim     x            trie_ls(1,1,P_dlam),
! jbao new gfs not needed by fim     x            epse,epso,ls_node)
! jbao new gfs not needed by fim      call f_hpmstop(61)
cc
! jbao new gfs not needed by fim      call f_hpmstart(62,"gr delnpo")
! jbao new gfs not needed by fim      call delnpo(trio_ls(1,1,P_q   ),
! jbao new gfs not needed by fim     x            trie_ls(1,1,P_dphi),
! jbao new gfs not needed by fim     x            trio_ls(1,1,P_dlam),
! jbao new gfs not needed by fim     x            epse,epso,ls_node)
! jbao new gfs not needed by fim      call f_hpmstop(62)
cc
cc
! jbao new gfs not needed by fim      call f_hpmstart(63,"gr dezouv dozeuv")
!
! jbao new gfs not needed by fim!$omp parallel do shared(trie_ls,trio_ls)
! jbao new gfs not needed by fim!$omp+shared(epsedn,epsodn,snnp1ev,snnp1od,ls_node)
! jbao new gfs not needed by fim!$omp+private(k)
! jbao new gfs not needed by fim      do k=1,levs
! jbao new gfs not needed by fim         call dezouv(trie_ls(1,1,P_di +k-1), trio_ls(1,1,P_ze +k-1),
! jbao new gfs not needed by fim     x               trie_ls(1,1,P_uln+k-1), trio_ls(1,1,P_vln+k-1),
! jbao new gfs not needed by fim     x               epsedn,epsodn,snnp1ev,snnp1od,ls_node)
cc
! jbao new gfs not needed by fim         call dozeuv(trio_ls(1,1,P_di +k-1), trie_ls(1,1,P_ze +k-1),
! jbao new gfs not needed by fim     x               trio_ls(1,1,P_uln+k-1), trie_ls(1,1,P_vln+k-1),
! jbao new gfs not needed by fim     x               epsedn,epsodn,snnp1ev,snnp1od,ls_node)
! jbao new gfs not needed by fim      enddo
! jbao new gfs not needed by fim      call f_hpmstop(63)
cc
! jbao new gfs not needed by fim!sela print*,'completed call to dztouv'
cc
! jbao?? is next line needed
!cmr  call mpi_barrier (mpi_comm_world,ierr)
cc
      CALL countperf(0,5,0.)
      CALL synctime()
      CALL countperf(1,5,0.)
!!
! jbao new gfs dimg is not needed comment out??
! jbao new gfs physics      dimg=0
      CALL countperf(0,1,0.)
cc
!jbao new gfs not needed by fim      call f_hpmstart(67,"gr sumfln")
cc
!jbao new gfs not needed by fim!sela print*,'begining  call to sumfln'
!jbao new gfs not needed by fim      call sumflna_r(trie_ls(1,1,P_ze),
!jbao new gfs not needed by fim     x            trio_ls(1,1,P_ze),
!jbao new gfs not needed by fim     x            lat1s_r,
!jbao new gfs not needed by fim     x            plnev_r,plnod_r,
!jbao new gfs not needed by fim     x            lots,ls_node,latr2,
!jbao new gfs not needed by fim     x            lslag,lats_dim_a,lots,for_gr_r_1,
!jbao new gfs not needed by fim     x            ls_nodes,max_ls_nodes,
!jbao new gfs not needed by fim     x            lats_nodes_r,global_lats_r,
!jbao new gfs not needed by fim     x            lats_node_r,ipt_lats_node_r,lon_dims_r,dimg,
!jbao new gfs not needed by fim     x            lonsperlar,lonrx,latr)
cc
!jbao new gfs not needed by fim!sela print*,'completed call to sumfln'
!jbao new gfs not needed by fim      call f_hpmstop(67)
cc
      CALL countperf(1,1,0.)
cc
! -----------------------------------
!jbao new gfs not needed by fim      if( vertcoord_id.eq.3. ) then
! -----------------------------------
!jbao new gfs not needed by fim      CALL countperf(0,1,0.)                                            ! hmhj
!
!jbao new gfs not needed by fim      call f_hpmstart(68,"gr sumder2")                                  ! hmhj
!
!jbao new gfs not needed by fim      call sumdera_r(trie_ls(1,1,P_te),                                 ! hmhj
!jbao new gfs not needed by fim     x             trio_ls(1,1,P_te),                                   ! hmhj
!jbao new gfs not needed by fim     x             lat1s_r,                                             ! hmhj
!jbao new gfs not needed by fim     x             pddev_r,pddod_r,                                     ! hmhj
!jbao new gfs not needed by fim     x             levs,ls_node,latr2,                                  ! hmhj
!jbao new gfs not needed by fim     x             lslag,lats_dim_r,lotd,                               ! hmhj
!jbao new gfs not needed by fim     x             dyn_gr_r_1,                                          ! hmhj
!jbao new gfs not needed by fim     x             ls_nodes,max_ls_nodes,                               ! hmhj
!jbao new gfs not needed by fim     x             lats_nodes_r,global_lats_r,                          ! hmhj
!jbao new gfs not needed by fim     x             lats_node_r,ipt_lats_node_r,lon_dims_r,dimg,         ! hmhj
!jbao new gfs not needed by fim     x             lonsperlar,lonrx,latr)                               ! hmhj
!
!jbao new gfs not needed by fim      call f_hpmstop(68)                                                ! hmhj
!
!jbao new gfs not needed by fim      CALL countperf(1,1,0.)                                            ! hmhj
! --------------------------------
!jbao new gfs not needed by fim      endif     ! vertcoord_id=3
! --------------------------------
!
 
!cmr  call mpi_barrier (mpi_comm_world,ierr)
 
! jbao not needed by fim      do lan=1,lats_node_r
! jbao new gfs for fim not needed       timer1=rtc()
! jbao not needed by fimcc
! jbao not needed by fim         lat = global_lats_r(ipt_lats_node_r-1+lan)
! jbao not needed by fimcc
! jbao ?? do we need lon_dims_r if so how do we get it in?
! jbao not needed by fim         lon_dim = lon_dims_r(lan)
cc
! jbao not needed by fim         lons_lat = lonsperlar(lat)

! -------------------------------------------------------
! jbao not needed by fim         if( gen_coord_hybrid .and. vertcoord_id.eq.3. ) then
! -------------------------------------------------------
!
! jbao not needed by fim           lmax = min(jcap,lons_lat/2)                                  ! hmhj
! jbao not needed by fim!
! jbao not needed by fim          ipt_ls=min(lat,latr-lat+1)                                    ! hmhj
! jbao not needed by fim
! jbao not needed by fim          do i=1,lmax+1                                                 ! hmhj
! jbao not needed by fim            if ( ipt_ls .ge. lat1s_r(i-1) ) then                        ! hmhj
! jbao not needed by fim               reall=i-1                                                ! hmhj
! jbao not needed by fim               rlcs2(i)=reall*rcs2_r(ipt_ls)/rerth                      ! hmhj
! jbao not needed by fim            else                                                        ! hmhj
! jbao not needed by fim               rlcs2(i)=cons0     !constant                             ! hmhj
! jbao not needed by fim            endif                                                       ! hmhj
! jbao not needed by fim          enddo                                                         ! hmhj
!
! jbao not needed by fim!$omp parallel do private(k,i)
! jbao not needed by fim          do k=1,levs                                                   ! hmhj
! jbao not needed by fim            do i=1,lmax+1                                               ! hmhj
! jbao not needed by fim!
! jbao not needed by fim!           d(t)/d(lam)                                                 ! hmhj
! jbao not needed by fim               dyn_gr_r_1(2*i-1+(kdtlam-2+k)*lon_dim,lan)=              ! hmhj
! jbao not needed by fim     x        -for_gr_r_1(2*i  +(kst   -2+k)*lon_dim,lan)*rlcs2(i)      ! hmhj
! jbao not needed by fim               dyn_gr_r_1(2*i  +(kdtlam-2+k)*lon_dim,lan)=              ! hmhj
! jbao not needed by fim     x         for_gr_r_1(2*i-1+(kst   -2+k)*lon_dim,lan)*rlcs2(i)      ! hmhj
! jbao not needed by fim            enddo                                                       ! hmhj
! jbao not needed by fim          enddo                                                         ! hmhj
! jbao not needed by fim! --------------------
! jbao not needed by fim        endif       ! gc and vertcoord_id=3
! ---------------------
!
cc
! jbao not needed by fim         CALL countperf(0,6,0.)
!jbao newgfs not needed in fim!sela    print*,' beginning call four2grid',lan
!jbao newgfs not needed in fim         CALL FOUR2GRID_thread(for_gr_r_1(1,lan),for_gr_r_2(1,lan),     &
!jbao newgfs not needed in fim     &                  lon_dim,lons_lat,lonrx,5*levs+levh+3,lan,me)

! -------------------------------------------------------
! jbao not needed by fim        if( gen_coord_hybrid.and.vertcoord_id.eq.3. ) then              ! hmhj
! -------------------------------------------------------
! jbao not needed by fim           CALL FOUR2GRID_thread(dyn_gr_r_1(1,lan),dyn_gr_r_2(1,lan),   ! hmhj
! jbao not needed by fim     &                    lon_dim,lons_lat,lonrx,levs,lan,me)           ! hmhj
! jbao not needed by fim           CALL FOUR2GRID_thread(dyn_gr_r_1((kdtlam-1)*lon_dim+1,lan),  ! hmhj
! jbao not needed by fim     &                           dyn_gr_r_2((kdtlam-1)*lon_dim+1,lan),  ! hmhj
! jbao not needed by fim     &                    lon_dim,lons_lat,lonrx,levs,lan,me)           ! hmhj
! -------------------------
! jbao not needed by fim        endif       ! gc and vertcoord_id=3
! -------------------------

! jbao not needed by fim!sela    print*,' completed call four2grid lan=',lan
! jbao ?? new gfs do we need the following?
! jbao not needed by fim         CALL countperf(1,6,0.)
!!
! jbao do we need to this?? gen_coord_hyrbrid is false does it do this?
! jbao not needed by fim        if( .not. gen_coord_hybrid ) then                               ! hmhj
! jbao not needed by fim
! jbao not needed by fim          do k = 1, LEVS
! jbao not needed by fim            kr = (KSR + k - 2) * lon_dim
! jbao not needed by fim            kt = (KST + k - 2) * lon_dim
! jbao not needed by fim            do j = 1, lons_lat
! jbao not needed by fim              if (for_gr_r_2(j+kr,lan) <= 0.0) then
! jbao not needed by fim                for_gr_r_2(j+kr,lan) = QMIN
! jbao not needed by fim              endif
! jbao not needed by fim              for_gr_r_2(j+kt,lan) = for_gr_r_2(j+kt,lan)               &
! jbao not needed by fim     &                             / (1.0 + FV*for_gr_r_2(j+kr,lan))
! jbao not needed by fim            enddo
! jbao not needed by fim          enddo
! jbao not needed by fim          kq = (KSQ - 1)*lon_dim
! jbao not needed by fim          do j = 1, lons_lat
! jbao not needed by fim            for_gr_r_2(j+kq,lan) = exp( for_gr_r_2(j+kq,lan) )
! jbao not needed by fim          enddo
! jbao not needed by fim
! jbao not needed by fim        endif                                                           ! hmhj
c
! jbao don't think next 2 are needed for fim
!        timer2=rtc()
!        global_times_r(lat,me+1)=timer2-timer1
c$$$    print*,'timeloopr',me,timer1,timer2,global_times_r(lat,me+1)
 
!!
! jbao not needed by fim      enddo   !lan
!
! jbao newgfs not needed in fim      call f_hpmstart(69,"gr lat_loop2")
!
!===> *** ...  starting latitude loop
!
      do lan=1,lats_node_r
cc
         lat = global_lats_r(ipt_lats_node_r-1+lan)
cc
         lons_lat = lonsperlar(lat)

!!
!$omp parallel do schedule(dynamic,1) private(lon,j,k,lon_dim)
!$omp+private(istrt,njeff,iblk,ku,kv,kd,kq,kt,kr,kx,ky,ks,n)
!$omp+private(vvel,gu,gv1,gd,gt,gr,gr1,gq,gphi,glam)
!$omp+private(gtv,gtvx,gtvy,sumq,xcp)
!$omp+private(cldcov_v,hprime_v,fluxr_v,f_ice,f_rain,r_rime)
!$omp+private(prslk,prsl,prsik,prsi)

        DO lon=1,lons_lat,NGPTC
!!
!jbao newgfs do we need to define lon_dimsr in do physics?
! jbao no need to use          lon_dim = lon_dims_r(lan)
! jbao newgfs as previous gloopr for fim  set njeff and istrt=1
!          NJEFF   = MIN(NGPTC,lons_lat-lon+1)
!          ISTRT   = lon
          njeff=1
          ISTRT   = 1
          IBLK    = 1
! jbao newgfs iblk is set to 1 before, here it will be 1 but should we just set to 1 ??
! jbao           if (NGPTC.ne.1) then
! jbao             IBLK  = lon/NGPTC+1
! jbao           else
! jbao             IBLK  = lon
! jbao           endif
! jbao not needed by fim          do k = 1, LEVS
! jbao not needed by fim            ku = lon - 1 + (KSU + k - 2)*lon_dim
! jbao not needed by fim            kv = lon - 1 + (KSV + k - 2)*lon_dim
! jbao not needed by fim            kd = lon - 1 + (KSD + k - 2)*lon_dim
! jbao not needed by fim            do j = 1, njeff
! jbao not needed by fim              gu(j,k)  = for_gr_r_2(j+ku,lan)
! jbao not needed by fim              gv1(j,k) = for_gr_r_2(j+kv,lan)
! jbao not needed by fim              gd(j,k)  = for_gr_r_2(j+kd,lan)
! jbao not needed by fim            enddo
! jbao not needed by fim          enddo

! jbao not needed by fim          if( gen_coord_hybrid ) then                                    ! hmhj

! jbao not needed by fim            do k=1,levs                                                  ! hmhj
! jbao not needed by fim              kt = lon - 1 + (KST + k - 2)*lon_dim
! jbao not needed by fim              kr = lon - 1 + (KSR + k - 2)*lon_dim
! jbao not needed by fim              do j=1,njeff                                               ! hmhj
! jbao not needed by fim                gtv(j,k)   = for_gr_r_2(j+kt,lan)
! jbao not needed by fim                gr(j,k)    = max(qmin, for_gr_r_2(j+kr,lan))
! jbao not needed by fim              enddo                                                      ! hmhj
! jbao not needed by fim            enddo                                                        ! hmhj
! --------------------------------------
! jbao not needed by fim            if( vertcoord_id.eq.3. ) then
! --------------------------------------
! jbao not needed by fim              do k=1,levs                                                ! hmhj
! jbao not needed by fim                kx = lon - 1 + (kdtlam + k - 2)*lon_dim
! jbao not needed by fim                ky = lon - 1 + (kdtphi + k - 2)*lon_dim
! jbao not needed by fim                do j=1,njeff                                             ! hmhj
! jbao not needed by fim                  gtvx(j,k)  = dyn_gr_r_2(j+kx,lan)
! jbao not needed by fim                  gtvy(j,k)  = dyn_gr_r_2(j+ky,lan)
! jbao not needed by fim                enddo                                                    ! hmhj
! jbao not needed by fim              enddo							 ! hmhj
! -----------------------------
! jbao not needed by fim            endif
! -----------------------------
! jbao thermodyn_id not defined?
! jbao not needed by fim            if( thermodyn_id.eq.3 ) then
! get dry temperature from enthalpy					! hmhj
! jbao not needed by fim              sumq=0.0							! hmhj
! jbao not needed by fim              xcp=0.0							! hmhj
! jbao not needed by fim              do i=1,ntrac						! hmhj
! jbao not needed by fim                if( cpi(i).ne.0.0 ) then				! hmhj
! jbao not needed by fim                ks=ksr+(i-1)*levs					! hmhj
! jbao not needed by fim                do k=1,levs						! hmhj
! jbao not needed by fim                  kr = lon - 1 + (ks + k - 2)*lon_dim			! hmhj
! jbao not needed by fim                  do j=1,njeff						! hmhj
! jbao not needed by fim                    sumq(j,k)=sumq(j,k)+for_gr_r_2(j+kr,lan)		! hmhj
! jbao not needed by fim                    xcp(j,k)=xcp(j,k)+cpi(i)*for_gr_r_2(j+kr,lan)	! hmhj
! jbao not needed by fim                  enddo							! hmhj
! jbao not needed by fim                enddo							! hmhj
! jbao not needed by fim                endif							! hmhj
! jbao not needed by fim              enddo							! hmhj
! jbao not needed by fim              do k=1,levs						! hmhj
! jbao not needed by fim                do j=1,njeff						! hmhj
! jbao not needed by fim                  xcp(j,k)=(1.-sumq(j,k))*cpi(0)+xcp(j,k)		! hmhj
! jbao not needed by fim                  gt(j,k)=gtv(j,k)/xcp(j,k)				! hmhj
! jbao not needed by fim                enddo							! hmhj
! jbao not needed by fim              enddo							! hmhj
! jbao not needed by fim            else if( thermodyn_id.le.1 ) then				! hmhj
! jbao not needed by fim! get dry temperature from virtual temperature				! hmhj
! jbao not needed by fim             do k=1,levs                                                ! hmhj
! jbao not needed by fim              do j=1,njeff                                              ! hmhj
! jbao not needed by fim                gt(j,k)    = gtv(j,k) / (1.+fv*gr(j,k))                 ! hmhj
! jbao not needed by fim              enddo                                                     ! hmhj
! jbao not needed by fim             enddo							! hmhj
! jbao not needed by fim            else
! jbao not needed by fim! get dry temperature from dry temperature             			! hmhj
! jbao not needed by fim             do k=1,levs                                                ! hmhj
! jbao not needed by fim              do j=1,njeff                                              ! hmhj
! jbao not needed by fim                gt(j,k)    = gtv(j,k)                                   ! hmhj
! jbao not needed by fim              enddo                                                     ! hmhj
! jbao not needed by fim             enddo 
! jbao not needed by fim            endif
! jbao not needed by fim
! jbao not needed by fim          else                                                          ! hmhj
!
! jbao not needed by fim            do k = 1, levs
! jbao not needed by fim              kt = lon - 1 + (KST + k - 2)*lon_dim
! jbao not needed by fim              kr = lon - 1 + (KSR + k - 2)*lon_dim
! jbao not needed by fim              do j = 1, njeff
! jbao not needed by fim                gt(j,k)  = for_gr_r_2(j+kt,lan)
! jbao not needed by fim                gr(j,k)  = for_gr_r_2(j+kr,lan)
! jbao not needed by fim              enddo
! jbao not needed by fim            enddo

! jbao not needed by fim          endif
!
!       Remaining tracers
!
! jbao not needed by fim          do n = 1, NTRAC-1
! jbao not needed by fim            do k = 1, levs
! jbao not needed by fim              kr = lon - 1 + (KSR + n*LEVS + k - 2)*lon_dim
! jbao not needed by fim              do j = 1, njeff
! jbao not needed by fim                gr1(j,k,n) = for_gr_r_2(j+kr,lan)
! jbao not needed by fim              enddo
! jbao not needed by fim            enddo
! jbao not needed by fim          enddo
! jbao not needed by fim          kq = lon - 1 + (KSQ - 1)*lon_dim
! jbao not needed by fim          kt = lon - 1 + (KSPPHI - 1)*lon_dim
! jbao not needed by fim          kr = lon - 1 + (KSPLAM - 1)*lon_dim
! jbao not needed by fim          do j = 1, njeff
! jbao not needed by fim            gq  (j) = for_gr_r_2(j+kq,lan)
! jbao not needed by fim            gphi(j) = for_gr_r_2(j+kt,lan)
! jbao not needed by fim            glam(j) = for_gr_r_2(j+kr,lan)
! jbao not needed by fim          enddo
!!
!  ---  vertical structure variables:   del,si,sl,prslk,prdel
!
! jbao not needed by fim          if( gen_coord_hybrid ) then                                    ! hmhj
! jbao not needed by fim!Moor       call  hyb2press_gc(njeff,ngptc,gq,gtv,prsi,prsl,prdel)       ! hmhj
! jbao not needed by fim            call  hyb2press_gc(njeff,ngptc,gq,gtv,prsi,prsl,prsik,prslk) ! hmhj
! jbao not needed by fim            call omegtes_gc(njeff,ngptc,rcs2_r(min(lat,latr-lat+1)),     ! hmhj
! jbao not needed by fim     &                     gq,gphi,glam,gtv,gtvx,gtvy,gd,gu,gv1,vvel)    ! hmhj
! jbao not needed by fim          elseif (hybrid) then
! jbao not needed by fim !Moor      call  hyb2press(njeff,ngptc,gq, prsi, prsl,prdel)
! jbao not needed by fim            call  hyb2press(njeff,ngptc,gq, prsi, prsl,prsik, prslk)
! jbao not needed by fim            call omegtes(njeff,ngptc,rcs2_r(min(lat,latr-lat+1)),
! jbao not needed by fim     &                   gq,gphi,glam,gd,gu,gv1,vvel)
! jbao not needed by fim          else
! jbao newgfs for fim comment these out??
! jbao not needed by fim !Moor      call  sig2press(njeff,ngptc,gq,sl,del,si, prsi, prsl,prdel)
! jbao not needed by fim            call  sig2press(njeff,ngptc,gq,sl,si,slk,sik,
! jbao not needed by fim     &                                        prsi,prsl,prsik,prslk)
! jbao not needed by fim            CALL countperf(0,12,0.)
! jbao not needed by fim            call omegast3(njeff,ngptc,levs,
! jbao not needed by fim     &                    gphi,glam,gu,gv1,gd,del,
! jbao not needed by fim     &                    rcs2_r(min(lat,latr-lat+1)),vvel,gq,sl)
! jbao not needed by fim          endif
! end jbao newgfs for fim comment these out??
!.....
! jbao this was one of them levr do we need to comment this out
          do k=1,levs
            do j=1,njeff
              prslk(j,k)    = fpkap(prsl(j,k)*1000.0)
              cldcov_v(j,k) = cldcov(k,istrt+j-1,lan)
            enddo
          enddo

!jbao          if (levs .lt. levs) then
!jbao            do j=1,njeff
!jbao              prsi(j,levs+1) = prsi(j,levp1)
!jbao              prsl(j,levs)   = (prsi(j,levp1)+prsi(j,levs)) * 0.5
!jbao              prsik(j,levs+1) = prslk(j,levp1)
!jbao              prslk(j,levs)   = fpkap(prsl(j,levs)*1000.0)
!jbao            enddo
!jbao          endif
! jbao both false doesn't do 
          if (ldiag3d .or. lggfs3d) then
            do k=1,levs
              do j=1,njeff
!Moor           prslk(j,k)    = fpkap(prsl(j,k)*1000.0)
! jbao ?? cldcov_v goes in grrad, but ldiag3d and lggfs3d are false so not defined 
                cldcov_v(j,k) = cldcov(k,istrt+j-1,lan)
              enddo
            enddo
          endif
!
          do j=1,njeff
            hprime_v(j) = hprime(1,istrt+j-1,lan)
          enddo
!
          do k=1, nfxr
            do j=1,njeff
!jbao?? new gfs commented out in previous fim version of gfs but fluxr_v needed in grrad
!              fluxr_v(j,k) = fluxr(k,istrt+j-1,lan)
            enddo
          enddo
          if (NUM_P3D == 3) then
            do k = 1, levs
              do j = 1, njeff
                f_ice (j,k) = phy_f3d(j,k,iblk,lan,1)
                f_rain(j,k) = phy_f3d(j,k,iblk,lan,2)
                r_rime(j,k) = phy_f3d(j,k,iblk,lan,3)
              enddo
            enddo
          endif
! jbao before this is done in gbphys????
! jbao          work1 = (log(coslat_r(lat) / (lons_lat*latg)) - dxmin) * dxinv
! jbao          work1 = max(0.0, min(1.0,work1))
! jbao          work2 = flgmin(1)*work1 + flgmin(2)*(1.0-work1)
! jbao          do j=1,njeff
! jbao            flgmin_l(j) = work2
! jbao          enddo
 
!  *** ...  calling radiation driver
 
!
!     lprnt = me .eq. 0 .and. kdt .ge. 120
!     if (lprnt) then
!     if (kdt .gt. 85) then
!     print *,' calling grrad for me=',me,' lan=',lan,' lat=',lat
!    &,' num_p3d=',num_p3d,' snoalb=',snoalb(lon,lan),' lon=',lon
!    &,' tsea=',tsea(lon,lan),' sncovr=',sncovr(lon,lan),
!    &' sheleg=',sheleg(lon,lan)
!
!jbao old gloopr          call grrad                                                    &
!jbao old gloopr!  ---  inputs:
!jbao old gloopr     &     ( prsi,prsl,prslk,gt,gr,gr1,vvel,slmsk(lon,lan),             &
!jbao old gloopr     &       xlon(lon,lan),xlat(lon,lan),tsea(lon,lan),                 &
!jbao old gloopr     &       sheleg(lon,lan),sncovr(lon,lan),snoalb(lon,lan),           &
!jbao old gloopr     &       zorl(lon,lan),hprime_v,                                    &
!jbao old gloopr     &       alvsf(lon,lan),alnsf(lon,lan),alvwf(lon,lan),              &
!jbao old gloopr     &       alnwf(lon,lan),facsf(lon,lan),facwf(lon,lan),              &
!jbao old gloopr                                          ! fice FOR SEA-ICE XW Nov04
!jbao old gloopr     &       fice(lon,lan),tisfc(lon,lan),                              &
!jbao old gloopr     &       solcon,coszen(lon,lan),coszdg(lon,lan),k1oz,k2oz,facoz,    &
!jbao old gloopr     &       cv(lon,lan),cvt(lon,lan),cvb(lon,lan),                     &
!jbao old gloopr     &       IOVR_SW,IOVR_LW,f_ice,f_rain,r_rime,flgmin_l,              &
!jbao old gloopr     &       NUM_P3D,NTCW-1,NCLD,NTOZ-1,NTRAC-1,NFXR,                   &
!jbao old gloopr     &       dtlw,dtsw, lsswr,lslwr,lssav,ldiag3d,sas_shal,norad_precip,&
!jbao old gloopr     &       crick_proof, ccnorm,lggfs3d,                               &
!jbao old gloopr!    &       dtlw,dtsw, lsswr,lslwr,lssav,ldiag3d,lggfs3d,              &
!jbao old gloopr     &       NGPTC,njeff,levs,IFLIP, me, lprnt,                         &
!jbao old gloopr!  ---  outputs:
!jbao old gloopr     &       swh(1,1,iblk,lan),sfcnsw(lon,lan),sfcdsw(lon,lan),         & ! sfcdsw FOR SEA-ICE XW Nov04
!jbao old gloopr     &       sfalb(lon,lan),                                            & ! lu [+1L]: add sfalb
!jbao old gloopr     &       hlw(1,1,iblk,lan),sfcdlw(lon,lan),tsflw(lon,lan),          &
!jbao old gloopr!  ---  input/output:
!jbao old gloopr     &       fluxr_v,cldcov_v                                           &
!jbao old gloopr!! ---  optional outputs:
!jbao old gloopr!!   &,      HTRSWB=htrswb(1,1,1,iblk,lan),                             &
!jbao old gloopr!!   &,      HTRLWB=htrlwb(1,1,1,iblk,lan)                              &
!jbao old gloopr     &     )
!
! jbao new gloopr call to grrad
          call grrad                                                    &
!  ---  inputs:
     &     ( prsi,prsl,prslk,gt,gr,gr1,vvel,slmsk(lon,lan),             &
     &       xlon(lon,lan),xlat(lon,lan),tsea(lon,lan),                 &
     &       sheleg(lon,lan),sncovr(lon,lan),snoalb(lon,lan),           &
     &       zorl(lon,lan),hprime_v,                                    &
     &       alvsf(lon,lan),alnsf(lon,lan),alvwf(lon,lan),              &
     &       alnwf(lon,lan),facsf(lon,lan),facwf(lon,lan),              &
                                          ! fice FOR SEA-ICE XW Nov04
     &       fice(lon,lan),tisfc(lon,lan),                              &
     &       solcon,coszen(lon,lan),coszdg(lon,lan),k1oz,k2oz,facoz,    &
     &       cv(lon,lan),cvt(lon,lan),cvb(lon,lan),                     &
     &       IOVR_SW,IOVR_LW,f_ice,f_rain,r_rime,flgmin_l,              &
     &       NUM_P3D,NTCW-1,NCLD,NTOZ-1,NTRAC-1,NFXR,                   &
     &       dtlw,dtsw, lsswr,lslwr,lssav,ldiag3d,sas_shal,norad_precip,&
     &       crick_proof, ccnorm,lggfs3d,                               &
     &       sscal,asymp,ext_cof,extlw_cof,                             &
!    &       dtlw,dtsw, lsswr,lslwr,lssav,ldiag3d,lggfs3d,              &
     &       NGPTC,njeff,levs,IFLIP, me, lprnt,                         &
!  ---  outputs:
     &       swh(1,1,iblk,lan),sfcnsw(lon,lan),sfcdsw(lon,lan),         & ! sfcdsw FOR SEA-ICE XW Nov04
     &       sfalb(lon,lan),                                            & ! l
     &       hlw(1,1,iblk,lan),sfcdlw(lon,lan),tsflw(lon,lan),          &
!  ---  input/output:
     &       fluxr_v,cldcov_v                                           &
!! ---  optional outputs:
!!   &,      HTRSWB=htrswb(1,1,1,iblk,lan),                             &
!!   &,      HTRLWB=htrlwb(1,1,1,iblk,lan)                              &
     &     )


!     if (lprnt) print *,' returned from grrad for me=',me,' lan=',
!    &lan,' lat=',lat,' kdt=',kdt
!
!
! jbao not needed in fim          if (ldiag3d .or. lggfs3d) then
! jbao not needed in fim            do k=1,levs
! jbao not needed in fim              do j=1,njeff
! jbao not needed in fim                cldcov(k,istrt+j-1,lan) = cldcov_v(j,k)
! jbao not needed in fim              enddo
! jbao not needed in fim            enddo
! jbao not needed in fim          endif
          do k=1,nfxr
            do j=1,njeff
              fluxr(k,istrt+j-1,lan) = fluxr_v(j,k)
            enddo
          enddo
! jbao original newgfs had levr .lt levs?? keep like this???
! jbao not needed in fim          if (levs .lt. levs) then
! jbao original newgfs had levr+1,levs keep like this???
! jbao not needed in fim            do k=levs+1,levs
! jbao not needed in fim              do j=1,njeff
! jbao not needed in fim                hlw(j,k,iblk,lan) = hlw(j,levs,iblk,lan)
! jbao not needed in fim                swh(j,k,iblk,lan) = swh(j,levs,iblk,lan)
! jbao not needed in fim              enddo
! jbao not needed in fim            enddo
! jbao not needed in fim          endif
 
c$$$          write(2900+lat,*) ' ilon = ',istrt
c$$$          write(2900+lat,'("swh",T16,"hlw")')
c$$$      do k=1,levs
c$$$         write(2900+lat,
c$$$     .         '(e10.3,T16,e10.3,T31,e10.3)')
c$$$     . swh(1,k,iblk,lan),hlw(1,k,iblk,lan)
c$$$       enddo
 
!!
          CALL countperf(1,12,0.)
          ENDDO
!
      enddo
cc
! jbao newgfs not in fim      call f_hpmstop(69)
!!
      CALL countperf(0,5,0.)
      CALL synctime()
      CALL countperf(1,5,0.)
!sela print*,'completed gloopr_v kdt=',kdt
!!
! jbao new gfs for fim ??? this is what it was in previous fim gfs??do we 
! need to allocate / deallocatek
       deallocate(sinlat_r,coslat_r)
      return
      end subroutine gloopr
