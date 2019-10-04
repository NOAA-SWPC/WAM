!#include "../../../ESMFVersionDefine.h"

!
! !module: gfs_dynamics_internal_state_mod 
!                         --- internal state definition of the
!                             esmf gridded component of the gfs dynamics.
!
! !description:  define the gfs dynamics internal state used to
!                                             create the esmf internal state.
!---------------------------------------------------------------------------
! !revision history:
!
!  november 2004      weiyu yang initial code.
!  may      2005      weiyu yang for the updated gfs version.
!  february 2006      shrinivas moorthi updated for the new version of gfs
!  january  2007      hann-ming henry juang for gfs dynamics only
!  may 2009           jun wang for gfs write grid component
!  oct 05   2009      sarah lu, grid_gr unfolded from 2D to 3D
!  oct 16   2009      sarah lu, add gfs_dyn_tracer
!  Nov 2009           weiyu yang, modified for the GEFS ensemble system.
!  Nov 2009           jun wang, add digital filter variables into internal state
!  Feb 2010           jun wang, add logical restart
!  Aug 2010           sarah lu, add glbsum and ptrc
!  Feb 2011           hann-ming henry jaung, add non-iteration dimensional-split
!                                            (NDSL) semi-Lagrangian dynamics
!  Feb 2011           sarah lu, add thermodyn_id, sfcpress_id
!  Sep 2011           weiyu yang, modified for using the ESMF 5.2.0r library.
!  Sep 2012           jun wang, add nemsio_in to specify input files
!  Aug    2016        P Pegion  add iniauinterval
!
! !interface:
!
      module gfs_dynamics_internal_state_mod

!!uses:
!------
      USE ESMF
      USE gfs_dynamics_namelist_mod
      USE gfs_dyn_machine, only: kind_grid, kind_evod
      use gfs_dyn_layout1
      use gfs_dyn_gg_def
      use gfs_dyn_vert_def
      use gfs_dyn_resol_def
      use gfs_dyn_date_def
      use namelist_dynamics_def
      use gfs_dyn_mpi_def
      use gfs_dyn_coordinate_def                                      ! hmhj
      use gfs_dyn_tracer_const                                        ! hmhj
      use gfs_dyn_dfi_mod, only : gfs_dfi_grid_gr                     ! jw
      use gfs_dyn_tracer_config, only: gfs_dyn_tracer_type, glbsum

      implicit none

! -----------------------------------------------
      type gfs_dynamics_internal_state		! start type define
! -----------------------------------------------

      type(nam_gfs_dyn_namelist)   :: nam_gfs_dyn
      type(gfs_dyn_state_namelist) :: esmf_sta_list
!      type(gfs_dyn_tracer_type)    :: gfs_dyn_tracer

      integer                   :: me, nodes
      integer                   :: lnt2_s, llgg_s
      integer                   :: lnt2
      integer                   :: grib_inp

!
      integer nxpt,nypt,jintmx
      integer jcap,levs,lonf,latg,lats_node_a_max
      integer ntrac, ntoz, ntcw, ncld, ntke, levr
      integer thermodyn_id, sfcpress_id
      integer npe_single_member

      character(16)                     ::  cfhour1
!jws
      integer                           ::  num_file
      character(32)        ,allocatable ::  filename_base(:)
      integer                           ::  ipt_lats_node_a
      integer                           ::  lats_node_a
      logical                           ::  adiabatic,iniauinterval
      logical                           ::  slg_flag
!jwe

      integer                           ::  nblck,kdt
      real                              ::  deltim

      integer              ,allocatable ::      lonsperlat (:)
      integer              ,allocatable ::      ls_node    (:)
      integer              ,allocatable ::      ls_nodes   (:, :)
      integer              ,allocatable ::  max_ls_nodes   (:)

      integer              ,allocatable ::  lats_nodes_a   (:)
      integer              ,allocatable ::  global_lats_a  (:)
      integer              ,allocatable ::  lats_nodes_ext (:)
      integer              ,allocatable ::  global_lats_ext(:)

! Add xlon, xlat, lats_nodes_a_fix for mGrid definition
      real (kind=kind_grid),allocatable ::  xlon(:,:),xlat(:,:)
      integer              ,allocatable ::  lats_nodes_a_fix (:)

      real(kind=kind_evod) ,allocatable ::        epse  (:)
      real(kind=kind_evod) ,allocatable ::        epso  (:)
      real(kind=kind_evod) ,allocatable ::        epsedn(:)
      real(kind=kind_evod) ,allocatable ::        epsodn(:)

      real(kind=kind_evod) ,allocatable ::       snnp1ev(:)
      real(kind=kind_evod) ,allocatable ::       snnp1od(:)

      integer              ,allocatable ::        ndexev(:)
      integer              ,allocatable ::        ndexod(:)

      real(kind=kind_evod) ,allocatable ::       plnev_a(:,:)
      real(kind=kind_evod) ,allocatable ::       plnod_a(:,:)
      real(kind=kind_evod) ,allocatable ::       pddev_a(:,:)
      real(kind=kind_evod) ,allocatable ::       pddod_a(:,:)
      real(kind=kind_evod) ,allocatable ::       plnew_a(:,:)
      real(kind=kind_evod) ,allocatable ::       plnow_a(:,:)

      real(kind=kind_evod) ,allocatable ::       trie_ls(:,:,:)
      real(kind=kind_evod) ,allocatable ::       trio_ls(:,:,:)

      INTEGER                               :: TRIEO_TOTAL_SIZE
      INTEGER, ALLOCATABLE, DIMENSION(:)    :: TRIE_LS_SIZE
      INTEGER, ALLOCATABLE, DIMENSION(:)    :: TRIO_LS_SIZE
      INTEGER, ALLOCATABLE, DIMENSION(:)    :: TRIEO_LS_SIZE
      INTEGER, ALLOCATABLE, DIMENSION(:)    :: LS_MAX_NODE_GLOBAL
      INTEGER, ALLOCATABLE, DIMENSION(:, :) :: LS_NODE_GLOBAL

      real(kind=kind_evod) ,allocatable ::      syn_ls_a(:,:,:)
      real(kind=kind_evod) ,allocatable ::      dyn_ls_a(:,:,:)

      real(kind=kind_evod) ,allocatable ::      grid_gr (:, :, :)
      real(kind=kind_evod) ,allocatable ::      grid_gr6(:, :, :)
!* digital filter (jun wang)
      type(gfs_dfi_grid_gr)             ::      grid_gr_dfi
!* stochastic physics addtions(philip pegion)
!
      real(kind=kind_evod) ,allocatable ::    syn_gr_a_1(:,:)
      real(kind=kind_evod) ,allocatable ::    syn_gr_a_2(:,:)
      real(kind=kind_evod) ,allocatable ::    pyn_gr_a_1(:,:)
      real(kind=kind_evod) ,allocatable ::    pyn_gr_a_2(:,:)
      real(kind=kind_evod) ,allocatable ::    dyn_gr_a_1(:,:)
      real(kind=kind_evod) ,allocatable ::    dyn_gr_a_2(:,:)
      real(kind=kind_evod) ,allocatable ::    anl_gr_a_1(:,:)
      real(kind=kind_evod) ,allocatable ::    anl_gr_a_2(:,:)
! stochastic physics weights
      real(kind=kind_evod) ,allocatable ::    shum_wts(:,:,:)
      real(kind=kind_evod) ,allocatable ::    sppt_wts(:,:,:)
      real(kind=kind_evod) ,allocatable ::    skebu_wts(:,:,:)
      real(kind=kind_evod) ,allocatable ::    skebv_wts(:,:,:)
      real(kind=kind_evod) ,allocatable ::    vcu_wts(:,:,:)
      real(kind=kind_evod) ,allocatable ::    vcv_wts(:,:,:)

!!
! for nasa ozon production and distruction rates:(input throu fixio_r)
      integer 	lev,levmax
      real 	phour
      integer 	kfhour
!
      real (kind=kind_grid) pdryini
      real (kind=kind_grid) ,allocatable ::	pwat(:,:)
      real (kind=kind_grid) ,allocatable ::	ptot(:,:)
      real (kind=kind_grid) ,allocatable ::     ptrc(:,:,:)      !glbsum

      integer              init,jcount,jpt,node,ibmsign,lon_dim,ilat

      real(kind=kind_evod) colat1, rone, rlons_lat, scale_ibm

      integer   p_gz,p_zslam,p_zsphi,p_dlam,p_dphi,p_uln,p_vln
      integer   p_zem,p_dim,p_tem,p_rm,p_dpm,p_qm
      integer   p_ze ,p_di ,p_te ,p_rq,p_dp ,p_q
      integer   p_w  ,p_x  ,p_y  ,p_rt,p_dpn,p_zq
      integer   p_zz ,p_dpphi,p_dplam,p_zzphi,p_zzlam

      integer   g_uum,g_vvm,g_ttm,g_rm,g_dpm,g_qm,g_gz,g_zz
      integer   g_uu ,g_vv ,g_tt ,g_rq ,g_dp ,g_q
      integer   g_uup,g_vvp,g_ttp,g_rqp,g_dpp,g_zqp, g_rqtk
      integer   g_u  ,g_v  ,g_t  ,g_rt ,g_dpn,g_zq ,g_p ,g_dpdt

      integer              lotls,lotgr,lots,lots_slg,lotd,lota,lotp

      integer              ibrad,ifges,ihour,ini,j,jdt,ksout,maxstp
      integer              mdt,idt,timetot,timer,time0
      integer              mods,n1,n2,n3,n4,ndgf,ndgi,nfiles,nflps
      integer              n1hyb, n2hyb
      integer              nges,ngpken,niter,nnmod,nradf,nradr
      integer              nsfcf,nsfci,nsfcs,nsigi,nsigs,nstep
      integer              nznlf,nznli,nznls,id,iret,nsout,ndfi

      integer              ierr,iprint,k,l,locl,n
      integer              lan,lat
      integer              spectral_loop

! carry fcst hour and initial date for imp/export
      real(kind=kind_evod) ,allocatable ::  fhour_idate(:,:)
      real(kind=kind_evod) chour
      real(kind=kind_evod) zhour
      integer              nfcstdate7(7)

      logical restart_run
      logical start_step, reset_step, end_step, restart_step

      logical nemsio_in
      logical lsout,ldfi
      logical :: dfiend_step = .false.

      real(kind=kind_evod),allocatable :: tee1(:)

      integer ikey,nrank_all,kcolor

      real(kind=kind_evod) cons0p5,cons1200,cons3600,cons0

      LOGICAL :: Cpl_flag
      LOGICAL ENS
!
! -----------------------------------------------------
      end type gfs_dynamics_internal_state		! end type define
! -----------------------------------------------------

! this state is supported by c pointer not f90 pointer, thus
! need this wrap.
!-----------------------------------------------------------
      type gfs_dyn_wrap		! begin type define
          type (gfs_dynamics_internal_state), pointer :: int_state
      end type gfs_dyn_wrap	! end type define

      end module gfs_dynamics_internal_state_mod
