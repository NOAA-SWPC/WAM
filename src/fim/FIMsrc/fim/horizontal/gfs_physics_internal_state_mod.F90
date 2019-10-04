!
! !module: gfs_physics_internal_state_mod 
!                         --- internal state definition of the
!                             esmf gridded component of the gfs physics.
!
! !description:  define the gfs physics internal state used to
!                                             create the esmf internal state.
!---------------------------------------------------------------------------
! !revision history:
!
!  november 2004      weiyu yang initial code.
!  may 2005           weiyu yang for the updated gfs version.
!  february 2006      shrinivas moorthi updated for the new version of gfs
!  january 2007       hann-ming henry juang for gfs dynamics only
!  july    2007       shrinivas moorthi for gfs physics only
!  november 2007       hann-ming henry juang continue for gfs physics
!  february 2009      tom henderson adapt for FIM from nems r3038
!
! !interface:
!
      module gfs_physics_internal_state_mod
!SMS$IGNORE BEGIN

!!uses:
!------
!TBH remove nam_gfs_phy_namelist for now
!TBH      use gfs_physics_namelist_mod, ONLY: nam_gfs_phy_namelist, gfs_phy_state_namelist
      use gfs_physics_namelist_mod, ONLY: gfs_phy_state_namelist
      use gfs_physics_sfc_flx_mod,  ONLY: Sfc_Var_Data, Flx_Var_Data

      use machine, only: kind_phys, kind_rad, kind_evod
      implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!TBH:  many fields removed for now, restoring them as I go...  

! -----------------------------------------------
      type gfs_physics_internal_state		! start type define
! -----------------------------------------------

!TBH      type(nam_gfs_phy_namelist)   :: nam_gfs_phy
      type(gfs_phy_state_namelist) :: esmf_sta_list

      TYPE(Sfc_Var_Data)        :: sfc_fld
      TYPE(Flx_Var_Data)        :: flx_fld

      integer                   :: me, nodes
      INTEGER                   :: llgg_s, lonr_s, latr_s
!     integer                   :: grib_inp

!
      integer ntrac,nxpt,nypt,jintmx,jcap,levs,lonr,latr,lats_node_r_max
      integer ntoz, ntcw, ncld, lsoil, nmtvr, num_p2d,levr
      integer :: num_p3d = 4
      logical :: ras=.false.
      integer lats_node_r
      integer thermodyn_id, sfcpress_id

      !TODO:  move to module resol_def
      integer nfxr

      character(16)                     ::  cfhour1

      integer                           ::  nblck, kdt
      real(kind=kind_phys)              ::  deltim

      integer              ,allocatable ::      lonsperlar (:)
      integer              ,allocatable ::  lats_nodes_r   (:)
      integer              ,allocatable ::  global_lats_r  (:)
      integer              ,allocatable ::  lats_nodes_ext (:)
      integer              ,allocatable ::  global_lats_ext(:)

      integer   lotgr
      real(kind=kind_evod) ,pointer :: ps(:)
      real(kind=kind_evod) ,pointer :: dp(:,:)
      real(kind=kind_evod) ,pointer :: dpdt(:,:)
      real(kind=kind_evod) ,pointer :: p(:,:)
      real(kind=kind_evod) ,pointer :: u(:,:)
      real(kind=kind_evod) ,pointer :: v(:,:)
      real(kind=kind_evod) ,pointer :: t(:,:)
      real(kind=kind_evod) ,pointer :: q(:,:)
      real(kind=kind_evod) ,pointer :: oz(:,:)
      real(kind=kind_evod) ,pointer :: cld(:,:)

      REAL(KIND=KIND_RAD) ,ALLOCATABLE :: XLON(:,:),XLAT(:,:)
      REAL(KIND=KIND_RAD) ,ALLOCATABLE :: COSZDG(:,:)
      REAL(KIND=KIND_RAD) ,ALLOCATABLE :: sfalb(:,:)
      REAL(KIND=KIND_RAD) ,ALLOCATABLE :: CLDCOV(:,:,:)
      REAL(KIND=KIND_RAD) ,ALLOCATABLE :: HPRIME(:,:,:)
      REAL(KIND=KIND_RAD) ,ALLOCATABLE :: SWH(:,:,:,:),HLW(:,:,:,:)
      REAL(KIND=KIND_RAD) ,ALLOCATABLE :: FLUXR(:,:,:)
!!

      REAL(KIND=KIND_RAD) ,ALLOCATABLE :: phy_f3d(:,:,:,:,:)
      REAL(KIND=KIND_RAD) ,ALLOCATABLE :: phy_f2d(:,:,:)
!
! carry fhour and initial date, may not be necessary later
      real(kind=kind_evod) ,allocatable :: fhour_idate(:,:)
      real(kind=kind_rad) :: phour
      INTEGER :: KFHOUR
      real, allocatable  :: poz(:),ozplin(:,:,:,:)
!     FOR OZON INTERPOLATION:
      INTEGER,ALLOCATABLE:: JINDX1(:),JINDX2(:)
!
      REAL,ALLOCATABLE:: DDY(:)
      !TBH:  promoted to arrays for FIM
      !TODO:  Is this really needed?  
      REAL(KIND=KIND_RAD) ,ALLOCATABLE :: SLAG(:,:),SDEC(:,:),CDEC(:,:)

      !TBH:  added for FIM to avoid SAVE
      REAL(KIND=KIND_PHYS) ,ALLOCATABLE :: acv(:,:),acvb(:,:),acvt(:,:)

      !TBH:  added for FIM to mimic NMMB
      !TBH:  memory bounds
      INTEGER :: ims,ime
      !TBH:  patch (distributed-memory loop) bounds
      INTEGER :: ips,ipe

!!
! for nasa ozon production and distruction rates:(input throu fixio_r)
      integer 	lev,levmax
!
      integer              init,jcount,jpt,node
      integer              ibmsign
      integer              lon_dim,ilat

      real(kind=kind_evod) colat1
!!
!     real(kind=kind_evod) rone
!     real(kind=kind_evod) rlons_lat
!     real(kind=kind_evod) scale_ibm


!     integer              ibrad,ifges,ihour,ini,j,jdt,ksout,maxstp
!     integer              mdt,idt,timetot,timer,time0
!     integer              mods,n1,n2,n3,n4,ndgf,ndgi,nfiles,nflps
!     integer              nges,ngpken,niter,nnmod,nradf,nradr
!     integer              nsfcf,nsfci,nsfcs,nsigi,nsigs,nstep
!     integer              nznlf,nznli,nznls,id,iret,nsout

      integer              iret, n3, n4

      integer              ierr,iprint,k,l,locl,n
      integer              lan,lat

      real(kind=kind_phys) chour
      real(kind=kind_phys) zhour

!     logical start_step
!     logical end_step
      logical lsout

      integer ikey,nrank_all,kcolor

      real(kind=kind_phys) cons0p5,cons1200,cons3600    !constant
      real(kind=kind_phys) cons0                        !constant

!
! -----------------------------------------------------
      end type gfs_physics_internal_state		! end type define
! -----------------------------------------------------

! this state is supported by c pointer not f90 pointer, thus
! need this wrap.
!-----------------------------------------------------------
      type gfs_phy_wrap		! begin type define
          type (gfs_physics_internal_state), pointer :: int_state
      end type gfs_phy_wrap	! end type define

!TBH:  temporary delcaration of int_state here
!TODO:  move this elsewhere...
      TYPE(gfs_physics_internal_state), pointer :: gis_phy

!SMS$IGNORE END
      end module gfs_physics_internal_state_mod
