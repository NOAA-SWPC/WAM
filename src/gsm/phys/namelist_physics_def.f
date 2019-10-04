      module namelist_physics_def

!! Code Revision
!! oct 12 2009     Sarah Lu, add grid_aldata
!! Jan 12 2010     Sarah Lu, add fdaer
!! June   2010     Shrinivas Moorthi - upgrade GFS physics
!! Aug 03 2010     Jun Wang, add fhdfi,ndfi,ldfi
!! Apr 06 2012     Henry Juang, add idea
!! Oct 18 2012     Shrinivas Moorthi add use_ufo, dtphys and remove hdif_fac
!! Dec 06 2012     Jun Wang, add nemsio_in/nemsio_out,sfcio
!! Apr 09 2013     Jun Wang, add ivegsrc and cu_physics 
!! Oct 30 2013     Xingren Wu, add a2oi_out and ngrid_a2oi
!! Mar 18 2014     Sarah Lu, remove iaer_mdl
!! APr 18 2013     Xingren Wu, add cplflx
!! Apr 28 2014     Jun Wang, add prslrd0,cgwf
!! May  2,2014     Philip Pegion, add stochastic physics
!! Aug 01 2014     s moorthi - add fixtrc
!! Feb 05 2014     s moorthi - add stochphys
!! Mar 03 2016     J. Han - add cnvcld for cnv cloudiness enhancement
!! Mar 04 2016     J. Han - change newsas & sashal to imfdeepcnv
!                           & imfshalcnv, respectively
!! Mar 23 2016     s moorthi - add ral_ts - time scale for rayleigh friction
!! Mar 31 2016      Xu  Li   - change nst_fcst to be nstf_name

      use machine, ONLY: kind_evod
      implicit none

      ! PT Moved from gloopb
      logical :: flipv = .true.
      
      integer nszer,nsres,nslwr,nsout,nsswr,nscyc,ndfi,igen,jo3,ngptc
     &,       lsm,ens_mem,ncw(2),lsea,nsout_hf
      integer,target :: ivegsrc,cu_physics
      real(kind=kind_evod) fhswr,fhlwr,fhrot,fhseg,fhmax,fhout,fhres,
     & fhzer,fhini,fhcyc,fhdfi,crtrh(3),flgmin(2),
     & ccwf(2),dlqf(2),ctei_rm(2),fhgoc3d,fhout_hf,fhmax_hf,cdmbgwd(2),
     & bkgd_vdif_m, bkgd_vdif_h, psautco(2), prautco(2), evpco
     &,bkgd_vdif_s,wminco(2),prslrd0,cgwf(2),sup,dtphys,ral_ts

! iau parameters
      logical    :: iau = .false. ! iau forcing included
      integer    :: iaulnp=2 ! 0: no pressure update, 1: log pressure, 2: pressure
      integer    :: iau_delthrs = 6 ! iau time interval (to scale increments)
      character(len=120),   dimension(7) :: iaufiles_fg,iaufiles_anl
      real(kind=kind_evod), dimension(7) :: iaufhrs

      integer imfshalcnv,imfdeepcnv
      logical ldiag3d,ras,zhao_mic,crick_proof,ccnorm
      logical shal_cnv, cscnv, do_shoc, shoc_cld, shocaftcnv
      logical mom4ice,mstrat,trans_trac,moist_adj,lggfs3d,cal_pre
      logical lsfwd,lssav,lscca,lsswr,lslwr,ldfi
! idea hmhj add
      logical lsidea
      logical shuff_lats_r,reshuff_lats_r,reduced_grid
      logical hybrid,gen_coord_hybrid
!     logical hybrid,gen_coord_hybrid,zflxtvd
      logical pre_rad,random_clds,old_monin,cnvgwd 
      logical restart, nemsio_in, nemsio_out, sfcio_out, griboro
      logical climate
      logical use_ufo, semilag, fixtrc(25), gg_tracers
      logical pdfcld,shcnvcw,redrag,hybedmf,dspheat,stochphys,cnvcld
      logical a2oi_out,cplflx
      character*20 ens_nam

      integer ngrid_a2oi

      logical nst_anl    ! nsst analysis control: false (default) = off; true = on
      integer, dimension(5) :: nstf_name

!     nstf_name contains the NSST related parameters
!     nstf_name(1) : 0 = NSSTM off, 1 = NSSTM on but uncoupled, 2 =
!     NSSTM on and coupled
!     nstf_name(2) : 1 = NSSTM spin up on, 0 = NSSTM spin up off
!     nstf_name(3) : 1 = NSST analysis on, 0 = NSSTM analysis off
!     nstf_name(4) : zsea1 in mm
!     nstf_name(5) : zsea2 in mm

!
!     Radiation control parameters
!
      logical norad_precip
      integer isol, ico2, ialb, iems, iaer, iovr_sw, iovr_lw, ictm

      integer isubc_sw, isubc_lw

      real(kind=kind_evod), dimension(5) ::  skeb,sppt,shum,vcamp
      logical                         ::  do_skeb,do_sppt,do_shum,do_vc
      real(kind=kind_evod)               ::  vc
!
!     Chemistry control parameters                       
!
      logical grid_aldata           ! option to allocate grid_fld
      real(kind=kind_evod)  fdaer   ! relaxation time in days to gocart anal/clim
!
      end module namelist_physics_def
