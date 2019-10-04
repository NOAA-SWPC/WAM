      module resol_def
      
      implicit none
      
      integer   jcap,jcap1,jcap2,latg,latg2,latr,latr2
      integer   levh,levm1,levp1,levs,lnt,lnt2,lnt22,levr
      integer   lnte,lnted,lnto,lntod,lnuv
      integer   lonf,lonfx,lonr,lonrx
      integer   kdt_start
!
!jw      integer   ntrac
      integer,target ::   ntrac,ntke
      integer   nxpt,nypt,jintmx,latrd
      integer   ntoz,ntcw
      integer   lsoil,nmtvr,ncld,num_p3d,num_p2d,nrcm,npdf3d
      integer   ncnvcld3d
      integer   nshoc_3d, nshoc_2d, ntot3d, ntot2d
      integer   ngrids_sfcc, ngrids_flx, nfxr, ngrids_gg
      integer   ngrids_aer                                  ! for g2d_fld
      integer   ngrids_nst,nr_nst,nf_nst
      integer   ngrids_aoi                           ! for A/O/I coupling
!jws      integer   ivsupa, ivssfc, ivssfc_restart, ivsinp
      integer   ivsupa, ivssfc_restart, ivsinp
      integer,target  :: thermodyn_id, sfcpress_id                      ! hmhj
      integer,target  :: ivssfc,ivsnst
      integer   ngrids_sfcc2d,ngrids_sfcc3d
      integer   idvt
!jwe
!
      integer   nlunit
!jw      integer   thermodyn_id, sfcpress_id			! hmhj
!
      integer   g_gz, g_ps, g_t, g_u, g_v, g_q, g_p, g_dp, g_dpdt
      integer   lotgr

      integer kwq,kwte,kwdz,kwrq

!     For Ensemble concurrency run. Weiyu
!     INTEGER :: Ensemble_Id, Total_member

!     The option to add 2d/3d diag fields to physics export state
      logical :: lgocart

!     For GOCART mapping
      integer, allocatable :: scatter_lats(:)
      integer, allocatable :: scatter_lons(:)

      end module resol_def
!
