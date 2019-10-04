      module n_namelist_physics_def

!! Code Revision
!! oct 12 2009     Sarah Lu, add grid_aldata
!! Jan 12 2010     Sarah Lu, add fdaer
!! June   2010     Shrinivas Moorthi - upgrade GFS physics
!! Aug 03 2010     Jun Wang, add fhdfi,ndfi,ldfi

      use machine, ONLY: kind_evod
      implicit none
      
      integer nszer,nsres,nslwr,nsout,nsswr,nscyc,ndfi,igen,jo3,ngptc
     &,       lsm,ens_mem,ncw(2),lsea,nsout_hf,num_reduce
      real(kind=kind_evod) fhswr,fhlwr,fhrot,fhseg,fhmax,fhout,fhres,
     & fhzer,fhini,fhcyc,fhdfi,crtrh(3),flgmin(2),
     & ccwf(2),dlqf(2),ctei_rm(2),fhgoc3d,fhout_hf,fhmax_hf,cdmbgwd(2),
     & bkgd_vdif_m, bkgd_vdif_h, hdif_fac, psautco(2), prautco(2), evpco
     &,bkgd_vdif_s,wminco(2)
      logical ldiag3d,ras,zhao_mic,sashal,newsas,crick_proof,ccnorm
      logical shal_cnv
      logical mom4ice,mstrat,trans_trac,moist_adj,lggfs3d,cal_pre
      logical lsfwd,lssav,lscca,lsswr,lslwr,ldfi
      logical shuff_lats_r,reshuff_lats_r,reduced_grid
      logical hybrid,gen_coord_hybrid
!     logical hybrid,gen_coord_hybrid,zflxtvd
      logical pre_rad,random_clds,old_monin,cnvgwd 
      logical restart
!     logical restart, gfsio_in, gfsio_out
      character*20 ens_nam

      integer nst_fcst
      logical nst_spinup
!
!     Radiation control parameters
!
      logical norad_precip
      integer isol, ico2, ialb, iems, iaer, iovr_sw, iovr_lw, ictm

      integer isubc_sw, isubc_lw
!
!     Chemistry control parameters                       
!
      logical grid_aldata           ! option to allocate grid_fld
      real(kind=kind_evod)  fdaer   ! relaxation time in days to gocart anal/clim
!
      end module n_namelist_physics_def
