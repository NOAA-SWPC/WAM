      SUBROUTINE do_physics_one_step(deltim,kdt,PHOUR, 
     &                               grid_fld, sfc_fld,
     &                               flx_fld, nst_fld, g3d_fld,
     &                               g2d_fld, aoi_fld, importData,
     &                               lats_nodes_r,global_lats_r,
     &                               lonsperlar,XLON,XLAT,COSZDG, 
     &                               HPRIME,SWH,swhc,HLW,hlwc,
!    &                               HTRSWB,HTRLWB,      ! idea add
     &                               FLUXR,SFALB, SLAG,SDEC,CDEC,
     &                               OZPLIN,JINDX1,JINDX2, DDY,
     &                               phy_f3d,  phy_f2d, phy_fctd, nctp,
     &                               NBLCK,       zhour_dfi,n3, n4,
!    &                               NBLCK,ZHOUR, zhour_dfi,n3, n4,
     &                               LSOUT,COLAT1,CFHOUR1,restart_step,
     &                               mdl_parm)
!!

!!
!! Code Revision:
!! oct 11 2009     Sarah Lu, grid_gr is replaced by grid_fld
!! dec 01 2009     Sarah Lu, add CLDCOV/FCLD check print
!! dec 08 2009     Sarah Lu, add g3d_fld to gloopr calling argument
!! dec 15 2009     Sarah Lu, add g3d_fld to gloopb calling argument;
!!                           add DQDT check print
!! Feb 05 2010     J. Wang, write out restart file
!! Apr 10 2010     Sarah Lu, debug print removed
!! Jul 07 2010     S. Moorthi Added nst_fld and other changes
!! Jul 21 2010     Sarah Lu, output 2d aerosol diag fields
!! Aug 03 2010     Jun Wang, set llsav through ndfi,ldfi
!! Aug 10 2010     Sarah Lu, zerout g2d_fld if needed
!! Aug 25 2010     J. Wang, add half dfi filtered fields output
!! Sep 11 2010     Sarah Lu, g2d_fld zerout call modified
!! Apr 06 2012     Henry Juang, add idea
!! Oct 18 2012     S. Moorthi Added oro_uf and modifications to nst
!! Mar 08 2013     J. Wang, add restart capibility for idea
!! Mar 26 2014     Xingren Wu, add aoi_fld for A/O/I coupling
!! Mar 31 2014     S Moorthi Add ocn_tmp as the input argument and use
!!                           when it contains valid data - for coupled model
!! Jul -- 2014     S Moorthi - merge with GFS, fix init-micro call etc
!! jun    2014     y-t hou,  revised sw sfc spectral component fluxes
!!                           and ocean albedo (no ice contamination) for coupled mdl
!! Sep 16 2014     Shrinivas Moorthi - cleanup and rearrange argumets
!!                                       for gloopr and glopb
!! Sep 30 2014     Sarah Lu, remove fscav (the option to compute tracer
!!                           scavenging in GFS is disable)
!! Apr 13 2014     Shrinivas Moorthi  - do physics only for comp_task
!! Aug    2015     Xu Li,    change nst_fcst and nst_spinup to be nstf_name and
!!                           introduce the depth mean SST
!! Jun 09 2015     G Theurich    Generalize importData handling
!! Sep 25 2015     Xingren Wu   Connect importData to GSM for A/O/I coupling
!! Jan    2016     P. Tripp     Coupling import/exportFieldsList for NUOPC/GSM merge
!!                              ocn_tmp moved to importData
!! Feb 25 2016     S Moorthi - add kdt_dif and associated chnages
!! March    2016  Hang Lei      Add DDT mdl_parm for physics driver
!! Mar 24 2016     Xingren Wu   Connect Ice/Snow thickness (volume)
!! Aug 17 2016     P Pegion   add logic to know if in iauinterval
!!
      USE machine, ONLY: KIND_GRID, KIND_RAD, kind_phys
      use resol_def
      use layout1
      use vert_def
      use date_def
      use namelist_physics_def
      use physcons, only : rlapse
      use mpi_def
      use ozne_def
      use module_nst_water_prop, only: get_dtzm_2d
      use gfs_physics_sfc_flx_mod
      use gfs_physics_sfc_flx_set_mod
      use gfs_physics_gridgr_mod,   ONLY: Grid_Var_Data
      use gfs_physics_nst_var_mod,  ONLY: Nst_Var_Data
      use gfs_physics_aoi_var_mod,  ONLY: aoi_var_data
      use gfs_physics_g3d_mod,      ONLY: G3D_Var_Data
      use gfs_physics_g2d_mod,      ONLY: G2D_Var_Data, g2d_zerout
!     use gfs_phy_tracer_config,    ONLY: gfs_phy_tracer_type
      use wam_jh_integral,          ONLY: jh_integral_zero,
     &                                    do_jh_integral,
     &                                    write_jh_output,
     &                                    jh_nh_integral

      use d3d_def, ONLY: d3d_zero, CLDCOV
! idea add by hmhj
      use module_radsw_parameters,  only : NBDSW
      use module_radlw_parameters,  only : NBDLW
      use module_CPLFIELDS,         only : NImportFields, queryFieldList
     &,                                    importFieldsList 
     &,                                    importFieldsValid
      use nuopc_physics,            only: model_parameters
      IMPLICIT NONE
!!     
      TYPE(Sfc_Var_Data)        :: sfc_fld
      TYPE(Flx_Var_Data)        :: flx_fld
      TYPE(Grid_Var_Data)       :: grid_fld 
      TYPE(Nst_Var_Data)        :: nst_fld 
      TYPE(G3D_Var_Data)        :: g3d_fld 
      TYPE(G2D_Var_Data)        :: g2d_fld
      type(aoi_var_data)        :: aoi_fld
      logical                   :: iniauinterval
!     type(gfs_phy_tracer_type) :: gfs_phy_tracer
!*    REAL(KIND=KIND_GRID)      GRID_GR(lonr*lats_node_r_max,lotgr)
      CHARACTER(16)             :: CFHOUR1
      logical                   :: restart_step
!!    
! The following dummy array must be an explicit shape dummy array!!!!!
      REAL(KIND=KIND_EVOD),INTENT(IN):: 
     &  importData(lonr,lats_node_r,NImportFields)
! importData(:,:,n) containts the import fields:
! See module_CPLFIELDS.F90 for that list

      REAL(KIND=KIND_EVOD),INTENT(IN)    :: deltim,PHOUR
      REAL(KIND=KIND_EVOD),INTENT(INOUT) :: ZHOUR_DFI
!     REAL(KIND=KIND_EVOD),INTENT(INOUT) :: ZHOUR,ZHOUR_DFI
!!
      REAL(KIND=KIND_EVOD)  :: delt_cpl  ! xw - add for A/O/I coupling
!!     
      INTEGER n3, n4, nblck, nctp, kdt, lan, lons_lat, lat, lon, njeff,
     &        iblk, item
      character(len=128) :: fldname
!!
      INTEGER               LATS_NODES_R(NODES)
      integer, dimension(latr) :: global_lats_r, lonsperlar
      real (kind=kind_rad) dtzm(lonr,lats_node_r)
!!     
      real(kind=kind_evod) zsea1,zsea2
      real(kind=kind_evod) colat1, phyhour, phydt, dtp
      real (kind=kind_phys), dimension(lonr,lats_node_r) :: xlon, xlat,
     &                                                     coszdg, sfalb
      real (kind=kind_phys), dimension(ngptc,levs,nblck,lats_node_r) ::
     &                          swh, swhc, hlw, hlwc, dtjh_global
      REAL (KIND=KIND_RAD) HPRIME(NMTVR,LONR,LATS_NODE_R),
     &                     FLUXR(nfxr,LONR,LATS_NODE_R)
! idea add by hmhj  - commented by moorthi since unused
!    &,                    HTRSWB(NGPTC,LEVS,NBDSW,NBLCK,LATS_NODE_R)
!    &,                    HTRLWB(NGPTC,LEVS,NBDLW,NBLCK,LATS_NODE_R)

      REAL (kind=kind_phys) phy_f3d(NGPTC,LEVS,ntot3d,NBLCK,lats_node_r)
     &,                     phy_f2d(lonr,lats_node_r,ntot2d)
     &,                     phy_fctd(lonr,lats_node_r,nctp)
     &,                     DDY(LATS_NODE_R)

!     real(kind=kind_evod) global_times_r(latr,nodes)

      INTEGER JINDX1(LATS_NODE_R),JINDX2(LATS_NODE_R)
      REAL    OZPLIN(LATSOZP,LEVOZP,pl_coeff,timeoz) !OZONE PL Coeff

      REAL(KIND=KIND_EVOD) SLAG,SDEC,CDEC
      INTEGER              IERR,I,J,K,L,LOCL,N,iprint, findex, kdt_dif
      LOGICAL LSOUT
      real(kind=kind_phys), parameter:: omz1 = 10.0  ! for nst model
      real(kind=kind_phys) :: pf_nh_integral, jh_fac
!
!     real*8 rtc, timer1, timer2

      real (kind=kind_phys) dt_warm, tem1, tem2
      real (kind=kind_phys), save :: zhour_dfin=0.
! NUOPC physics driver types
       type(model_parameters), intent(in)     :: mdl_parm
!
      zsea1=0.001*real(nstf_name(4))
      zsea2=0.001*real(nstf_name(5))

!     SHOUR   = SHOUR + deltim
      shour   = kdt * deltim
      fhour   = shour / 3600.
      lsfwd   = kdt == 1
!jws
      kdt_dif = kdt - kdt_start
!     if (me == 0) write(0,*)' in do_onestep ndfi=',ndfi,' kdt_dif=',kdt_dif
!    &,' ldfi=',ldfi,' me=',me,' kdt=',kdt,' kdt_start=',kdt_start

      lssav = .true.
      iniauinterval = .false.
      if(ndfi > 0 .and. kdt_dif > ndfi/2 .and. 
     &                  kdt_dif <= ndfi .and. ldfi ) then
        lssav = .false.
      endif
      if(.not. ldfi .and. ndfi > 0. and. kdt_dif == ndfi/2+1) then
         zhour = zhour_dfin
      endif
!     if (me == 0) write(0,*)' in do_onestep ndfi=',ndfi,' kdt_dif=',
!    &            kdt_dif,' lssav=',lssav,' kdt=',kdt,' zhour=',zhour
!jwe
      lscca   = mod(KDT ,nsswr) == 0
      lsswr   = mod(KDT ,nsswr) == 1
      lslwr   = mod(KDT ,nslwr) == 1
! test repro
!     phyhour = phour + deltim/3600.
      phyhour = phour                           ! Moorthi No4 24, 2014
      phydt   = deltim

!     if (me == 0) write(0,*)' in do_onestep phyhour=',phyhour

      if(.not. semilag .and. lsfwd) phydt = 0.5*deltim
!
!jw now all the pes are fcst pe
!jw ifnems (.NOT.LIOPE.or.icolor.ne.2) then

      if (comp_task) then

        if (nscyc > 0 .and. mod(kdt,nscyc) == 1) then
          if (me == 0) print*,' calling gcycle at kdt=',kdt
          if ( nst_anl ) then     ! when nst analysis is on
            do j = 1, lats_node_r
              do i = 1, lonr
                if ( sfc_fld%slmsk(i,j) == 0 ) then
                  sfc_fld%tsea(i,j) = nst_fld%tref(i,j)
                endif
              enddo
            enddo

            call gcycle(me,lats_node_r,lonsperlar,global_lats_r,
     &                ipt_lats_node_r,idate,phour,fhcyc,
     &                xlon ,xlat, sfc_fld, ialb)

            do j = 1, lats_node_r
              do i = 1, lonr
                if ( sfc_fld%slmsk(i,j) == 0 ) then
                  nst_fld%tref(i,j) = sfc_fld%tsea(i,j)
                endif
              enddo
            enddo

            call get_dtzm_2d(nst_fld%xt,nst_fld%xz,nst_fld%dt_cool,
     &                       nst_fld%z_c,nst_fld%slmsk,zsea1,zsea2,
     &                       lonr,lats_node_r,dtzm)

            do j = 1, lats_node_r
              do i = 1, lonr
                if ( sfc_fld%slmsk(i,j) == 0 ) then
                  sfc_fld%tsea(i,j) = nst_fld%tref(i,j) + dtzm(i,j)
                endif
              enddo
            enddo

          else            ! when nst analysis is off
            call gcycle(me,LATS_NODE_R,LONSPERLAR,GLOBAL_LATS_R,
     &                ipt_lats_node_r,idate,phour,fhcyc,
     &                xlon ,xlat, sfc_fld, ialb)
          endif    ! if ( nst_anl) then
        endif      ! if (nscyc > 0 .and.
!
        if (num_p3d  ==  3) then        ! Ferrier Microphysics initialization
        dtp = min(phydt,dtphys)
          call INIT_MICRO(dtp,ngptc,levs,ntot3d,
     &                    nblck*lats_node_r, phy_f3d(1,1,1,1,1),
     &                    phour, me)
        endif
!
!-> Coupling insertion

        if (cplflx) then

!  -> Mask
!     ----
          fldname='land_mask'
          findex = QueryFieldList(ImportFieldsList,fldname)
          if (importFieldsValid(findex) .and. 
     &        importData(1,1,findex) > -99999.0) then
!$omp parallel do private(i,j)
            do j = 1, lats_node_r
              do i = 1, lonr
                aoi_fld%slimskin(i,j) = 1.0
                if (importData(i,j,findex) < 0.01) then
                   aoi_fld%FICEIN(i,j) = 0.0
                   aoi_fld%slimskin(i,j) = 3.0
                endif
              enddo
            enddo
          else
            if (me == 0) 
     &        write(0,*) 'do_physics_one_step skip field ',trim(fldname)
          endif

!  -> SST
!     ---
          fldname='sea_surface_temperature'
          findex = QueryFieldList(ImportFieldsList,fldname)
          if (importFieldsValid(findex) .and. 
     &        importData(1,1,findex) > -99999.0) then
            do j = 1, lats_node_r
              do i = 1, lonr
                if (aoi_fld%slimskin(i,j) < 3.1 .and.
     &              aoi_fld%slimskin(i,j) > 2.9) then
                  if (sfc_fld%slmsk(i,j) < 0.1 .or.
     &                sfc_fld%slmsk(i,j) > 1.9) then
                    sfc_fld%TSEA(i,j) = importData(i,j,findex)
                  endif
                endif
              enddo
            enddo
          else
            if (me == 0) 
     &        write(0,*) 'do_physics_one_step skip field ',trim(fldname)
          endif

        if (nstf_name(1) > 1) then     ! update tsea
          if (importFieldsValid(findex) .and. 
     &        importData(1,1,findex) > -99999.0) then
!
!         get tf from T1 (sfc_fld%tsea, OGCM layer 1 temperature) and
!         NSST-Profile
!
            call get_dtzm_2d(nst_fld%xt,nst_fld%xz,nst_fld%dt_cool,
     &                     nst_fld%z_c,sfc_fld%slmsk,
     &                     0.0,omz1,lonr,lats_node_r,dtzm)
!$omp parallel do private(j,i)
            do j = 1, lats_node_r
              do i = 1, lonr
                if ( sfc_fld%slmsk(i,j) == 0 ) then
                  nst_fld%tref(i,j) =  sfc_fld%tsea(i,j) - dtzm(i,j)
     &                 + (sfc_fld%oro(i,j)-sfc_fld%oro_uf(i,j))*rlapse
                endif
              enddo
            enddo
!
!           get tsea from Tf and NSST-Profile
!
            call get_dtzm_2d(nst_fld%xt,nst_fld%xz,nst_fld%dt_cool,
     &                     nst_fld%z_c,sfc_fld%slmsk,
     &                     zsea1,zsea2,lonr,lats_node_r,dtzm)

!$omp parallel do private(j,i)
            do j = 1, lats_node_r
              do i = 1, lonr
                if ( sfc_fld%slmsk(i,j) == 0 ) then
                  sfc_fld%tsea(i,j) = nst_fld%tref(i,j) + dtzm(i,j)
     &                  - (sfc_fld%oro(i,j)-sfc_fld%oro_uf(i,j))*rlapse
                endif
              enddo
            enddo

          else
!
!         get tf from sfc_fld%tsea and NSST-Profile
!
            call get_dtzm_2d(nst_fld%xt,nst_fld%xz,nst_fld%dt_cool,
     &                     nst_fld%z_c,sfc_fld%slmsk,
     &                     zsea1,zsea2,lonr,lats_node_r,dtzm)
!$omp parallel do private(j,i)
            do j = 1, lats_node_r
              do i = 1, lonr
                if ( sfc_fld%slmsk(i,j) == 0 ) then
                  nst_fld%tref(i,j) =  sfc_fld%tsea(i,j) - dtzm(i,j)
     &                 + (sfc_fld%oro(i,j)-sfc_fld%oro_uf(i,j))*rlapse
                endif
              enddo
            enddo
          endif           ! if (importData(1,1,findex) > -99999.0) then
        endif             ! if (nstf_name(1) > 1) then 

!  -> Ts
!     --
          fldname='surface_temperature'
          findex = QueryFieldList(ImportFieldsList,fldname)
          if (importFieldsValid(findex) .and. 
     &        importData(1,1,findex) > -99999.0) then
            do j = 1, lats_node_r
!$omp parallel do private(i)
              do i = 1, lonr
                if (aoi_fld%slimskin(i,j) < 3.1 .and.
     &             aoi_fld%slimskin(i,j) > 2.9) then
                  aoi_fld%TSEAIN(i,j) = importData(i,j,findex)
                endif
              enddo
            enddo
          else
            if (me == 0) 
     &        write(0,*) 'do_physics_one_step skip field ',trim(fldname)
          endif
!  -> SeaIce
!     ------
          fldname='ice_fraction'
          findex = QueryFieldList(ImportFieldsList,fldname)
          if (importFieldsValid(findex) .and. 
     &        importData(1,1,findex) > -99999.0) then
            do j = 1, lats_node_r
!$omp parallel do private(i)
              do i = 1, lonr
                if (importData(i,j,findex) > 0.15) then
                  if (aoi_fld%slimskin(i,j) < 3.1 .and.
     &                aoi_fld%slimskin(i,j) > 2.9) then
                   if (sfc_fld%slmsk(i,j) < 0.1 .or.
     &                 sfc_fld%slmsk(i,j) > 1.9) then
                      aoi_fld%FICEIN(i,j) = importData(i,j,findex)
                      sfc_fld%FICE(i,j) = importData(i,j,findex)
                      aoi_fld%slimskin(i,j) = 4.0
                      sfc_fld%slmsk(i,j) = 2.0
                      sfc_fld%TSEA(i,j) = aoi_fld%TSEAIN(i,j)
                    endif
                  endif
                else
                  if (aoi_fld%slimskin(i,j) > 2.9 .and.
     &                aoi_fld%slimskin(i,j) < 3.1 .and.
     &                sfc_fld%FICE(i,j) > 0.15) then
                     sfc_fld%FICE(i,j) = 0.0
                     sfc_fld%slmsk(i,j) = 0.0
                     sfc_fld%TSEA(i,j) = aoi_fld%TSEAIN(i,j)
                  endif
                endif
              enddo
            enddo
          else
            if (me == 0) 
     &        write(0,*) 'do_physics_one_step skip field ',trim(fldname)
          endif

          fldname='mean_ice_volume'
          findex = QueryFieldList(ImportFieldsList,fldname)
          if (importFieldsValid(findex) .and.
     &        importData(1,1,findex) > -99999.0) then
            do j = 1, lats_node_r
!$omp parallel do private(i)
              do i = 1, lonr
                if (aoi_fld%slimskin(i,j) > 2.5) then
                    aoi_fld%HICEIN(i,j)=importData(i,j,findex)
                    if (sfc_fld%FICE(i,j) > 0.15) then
                       sfc_fld%HICE(i,j)=aoi_fld%HICEIN(i,j)
                    else
                       sfc_fld%HICE(i,j)=0.
                    endif
                endif
              enddo
            enddo
          else
            if (me == 0)
     &        write(0,*) 'do_physics_one_step skip field ',trim(fldname)
          endif

          fldname='mean_snow_volume'
          findex = QueryFieldList(ImportFieldsList,fldname)
          if (importFieldsValid(findex) .and.
     &        importData(1,1,findex) > -99999.0) then
            do j = 1, lats_node_r
!$omp parallel do private(i)
              do i = 1, lonr
                if (aoi_fld%slimskin(i,j) > 2.5) then
                    aoi_fld%HSNOIN(i,j)=importData(i,j,findex)
                    if (sfc_fld%FICE(i,j) > 0.15) then
                       sfc_fld%SNWDPH(i,j)=aoi_fld%HSNOIN(i,j)
                    else
                       sfc_fld%SNWDPH(i,j)=0.
                    endif
                endif
              enddo
            enddo
          else
            if (me == 0)
     &        write(0,*) 'do_physics_one_step skip field ',trim(fldname)
          endif

          fldname='mean_up_lw_flx'
          findex = QueryFieldList(ImportFieldsList,fldname)
          if (importFieldsValid(findex) .and. 
     &        importData(1,1,findex) > -99999.0) then
            do j = 1, lats_node_r
!$omp parallel do private(i)
              do i = 1, lonr
                if (aoi_fld%slimskin(i,j) > 3.9 .and.
     &              aoi_fld%slimskin(i,j) < 4.1) then
                     aoi_fld%ULWSFCIN(i,j)=-importData(i,j,findex)
                endif
              enddo
            enddo
          else
            if (me == 0) 
     &        write(0,*) 'do_physics_one_step skip field ',trim(fldname)
          endif
          fldname='mean_laten_heat_flx'
          findex = QueryFieldList(ImportFieldsList,fldname)
          if (importFieldsValid(findex) .and. 
     &        importData(1,1,findex) > -99999.0) then
            do j = 1, lats_node_r
!$omp parallel do private(i)
              do i = 1, lonr
                if (aoi_fld%slimskin(i,j) > 3.9 .and.
     &              aoi_fld%slimskin(i,j) < 4.1) then
                     aoi_fld%DQSFCIN(i,j) =-importData(i,j,findex)
                endif
              enddo
            enddo
          else
            if (me == 0) 
     &        write(0,*) 'do_physics_one_step skip field ',trim(fldname)
          endif
          fldname='mean_sensi_heat_flx'
          findex = QueryFieldList(ImportFieldsList,fldname)
          if (importFieldsValid(findex) .and. 
     &        importData(1,1,findex) > -99999.0) then
            do j = 1, lats_node_r
!$omp parallel do private(i)
              do i = 1, lonr
                if (aoi_fld%slimskin(i,j) > 3.9 .and.
     &              aoi_fld%slimskin(i,j) < 4.1) then
                     aoi_fld%DTSFCIN(i,j) =-importData(i,j,findex)
                endif
              enddo
            enddo
          else
            if (me == 0) 
     &        write(0,*) 'do_physics_one_step skip field ',trim(fldname)
          endif
          fldname='mean_zonal_moment_flx'
          findex = QueryFieldList(ImportFieldsList,fldname)
          if (importFieldsValid(findex) .and. 
     &        importData(1,1,findex) > -99999.0) then
            do j = 1, lats_node_r
!$omp parallel do private(i)
              do i = 1, lonr
                if (aoi_fld%slimskin(i,j) > 3.9 .and.
     &              aoi_fld%slimskin(i,j) < 4.1) then
                     aoi_fld%DUSFCIN(i,j) =-importData(i,j,findex)
                endif
              enddo
            enddo
          else
            if (me == 0) 
     &        write(0,*) 'do_physics_one_step skip field ',trim(fldname)
          endif
          fldname='mean_merid_moment_flx'
          findex = QueryFieldList(ImportFieldsList,fldname)
          if (importFieldsValid(findex) .and. 
     &        importData(1,1,findex) > -99999.0) then
            do j = 1, lats_node_r
              do i = 1, lonr
                if (aoi_fld%slimskin(i,j) > 3.9 .and.
     &              aoi_fld%slimskin(i,j) < 4.1) then
                    aoi_fld%DVSFCIN(i,j) =-importData(i,j,findex)
                endif
              enddo
            enddo
          else
            if (me == 0) 
     &        write(0,*) 'do_physics_one_step skip field ',trim(fldname)
          endif

        endif   ! if cplflx

        if (lsswr .or. lslwr) then         ! Radiation Call!
! idea add by hmhj
          if(lsidea) then
            if(lsswr) then
              swh    = 0.
!             htrswb = 0.
            endif
            if(lslwr) then
              hlw    = 0.
!             htrlwb = 0.
            endif
          endif
!*        CALL GLOOPR (grid_gr,
          CALL GLOOPR (grid_fld, g3d_fld, aoi_fld, lats_nodes_r
     &,                GLOBAL_LATS_R, LONSPERLAR, phyhour
     &,                deltim, XLON, XLAT, COSZDG, flx_fld%COSZEN
     &,                sfc_fld%SLMSK,  sfc_fld%weasd,  sfc_fld%SNCOVR
     &,                sfc_fld%SNOALB, sfc_fld%ZORL,   sfc_fld%TSEA
     &,                HPRIME, SFALB,  sfc_fld%ALVSF,  sfc_fld%ALNSF
     &,                sfc_fld%ALVWF,  sfc_fld%ALNWF,  sfc_fld%FACSF
     &,                sfc_fld%FACWF,  sfc_fld%CV,     sfc_fld%CVT
     &,                sfc_fld%CVB, SWH, SWHC, HLW, HLWC, flx_fld%SFCNSW
     &,                flx_fld%SFCDLW, sfc_fld%FICE,   sfc_fld%TISFC
     &,                flx_fld%SFCDSW, flx_fld%sfcemis

     &,                flx_fld%TSFLW,  FLUXR, phy_f3d, phy_f2d
     &,                SLAG, SDEC, CDEC, NBLCK, KDT, mdl_parm
!    &,                HTRSWB,HTRLWB       !idea add by hmhj
     &                 )
!          if (iprint .eq. 1) print*,' me = fin gloopr ',me

        endif
!
!!
!*       call gloopb ( grid_gr,
         if (iau) call checkiauforcing(phour*3600,iniauinterval)
         call gloopb (grid_fld,     g3d_fld,       sfc_fld,
     &                flx_fld,      aoi_fld,       nst_fld,
     &                lats_nodes_r, global_lats_r, lonsperlar,
     &                phydt,        phyhour,       sfalb,  xlon,
     &                swh,          swhc,          hlw,    hlwc,
!    &                nbdsw,        nbdlw,         HTRSWB, HTRLWB,          !idea add by hmhj
     &                hprime,       slag,          sdec,   cdec,
     &                ozplin,       jindx1,        jindx2, ddy,
     &                phy_f3d,      phy_f2d,       phy_fctd, nctp,
     &                xlat,         nblck,  kdt,   restart_step,
     &                mdl_parm,     iniauinterval, pf_nh_integral,
     &                dtjh_global)
!
!!
      endif ! if (comp_task) then

!jw   endif !.NOT.LIOPE.or.icolor.ne.2
!--------------------------------------------
!
!      write(0,*)'in do one phys step, lsout=',lsout,'kdt=',kdt, 
!     &   'nszer=',nszer,'ldfi=',ldfi,'ndfi=',ndfi,
!     &   'fhour=',fhour,'zhour=',zhour,'zhour_dfin=',zhour_dfin,
!     &   'zhour_dfi=',zhour_dfi

      if (kdt /= 0) then
         call do_jh_integral(global_lats_r, lonsperlar, nblck,ngptc,me)
         jh_fac = pf_nh_integral/jh_nh_integral - 2.0
         do lan=1,lats_node_r
           lat         = global_lats_r(ipt_lats_node_r-1+lan)
           lons_lat    = lonsperlar(lat)
           do lon=1,lons_lat,ngptc
             njeff = min(ngptc,lons_lat-lon+1)
             iblk  = (lon-1)/ngptc + 1
             do k=1,levs
               do i=1,njeff
                 item = lon+i-1
                 grid_fld%t(item,lan,k) = grid_fld%t(item,lan,k) +
     &                                jh_fac * dtjh_global(i,k,iblk,lan)
               enddo
             enddo
           enddo
         enddo
         if (output_jh_integral .and. me.eq.0) call write_jh_output(kdt,
     &                                              pf_nh_integral)
      endif

      if (lsout .and. kdt /= 0 ) then
!WY bug fix.
!-----------
        IF(.NOT. ALLOCATED(SL)) ALLOCATE(SL(levs))
        IF(.NOT. ALLOCATED(SI)) ALLOCATE(SI(levs + 1))
        CALL WRTOUT_physics(phyhour,FHOUR,ZHOUR,IDATE,
     &                      SL,SI,
     &                      sfc_fld, flx_fld, nst_fld, g2d_fld,
     &                      fluxr,
     &                      global_lats_r,lonsperlar,nblck,
!    &                      lats_nodes_r,global_lats_r,lonsperlar,nblck,
     &                      COLAT1,CFHOUR1,pl_coeff,
     &                     'SFC.F','NST.F','FLX.F','D3D.F')
!
      endif ! if ls_out
!
! A/O/I coupling - fluxes accumulated at every atm time step for coupling
!     Xingren Wu
!
      if (cplflx) then        ! for NUOPC coupling
          delt_cpl = deltim
          call aoicpl_prep(deltim,delt_cpl,phyhour,fhour,idate,
     &                     aoi_fld,global_lats_r,lonsperlar)
      endif   ! if cplflx
!
       IF (kdt > 0 .and. mod(kdt,nsres) == 0) THEN
!           write(0,*)'wrt_restart_physics,kdt=',kdt,'nsres=',nsres
           CALL wrtout_restart_physics(sfc_fld, nst_fld, fhour,idate,
     &                      lats_nodes_r,global_lats_r,lonsperlar,
     &                      phy_f3d, phy_f2d, ngptc, nblck, ens_nam)
       endif
!
      IF (mod(kdt,nszer) == 0 .and. lsout.and.kdt /= 0) THEN
        call flx_init(flx_fld,ierr)
        if(ldfi .and. kdt_dif == ndfi/2) then
         zhour_dfi  = zhour
         zhour_dfin = fhour
        endif
        zhour = fhour
!$omp parallel do private(n,i,j)
        do j=1,lats_node_r
          do i=1,lonr
            do n=1,nfxr
              fluxr(n,i,j) = 0.0
            enddo
          enddo
        enddo
!
        call jh_integral_zero
!
        if (ldiag3d) then
!         if(me==0) print *, 'LU_CLDCOV: zero out d3d fields'
          call d3d_zero

!         if ( gfs_phy_tracer%doing_GOCART ) then
!           call g2d_zerout(gfs_phy_tracer, g2d_fld)
!         endif

        endif
!
        if ( lgocart ) then
          call g2d_zerout(g2d_fld,ierr)
        endif

      ENDIF
!
      if(ldfi .and. kdt_dif == ndfi) then
         zhour = zhour_dfi
      endif
!     print *,'in phys one,kdt=',kdt,'zhour=',zhour,                   &
!     &  'zhour_dfi=',zhour_dfi,'zhour_dfin=',zhour_dfin
 
      if(ndfi > 0 .and. kdt_dif == ndfi .and. ldfi ) then
        ldfi = .false.
      endif
!
      RETURN
      END

      subroutine do_physics_gridcheck(grid_gr,g_pnt,km,
     &                                global_lats_r,lonsperlar,chr)
      use machine
      use resol_def
      use layout1

      real(kind=kind_grid) grid_gr(lonr*lats_node_r_max,lotgr)
      integer,intent(in):: global_lats_r(latr),g_pnt,km
      integer,intent(in):: lonsperlar(latr)
      character*(*) chr

      integer lan,lat,lons_lat,k

      do lan=1,lats_node_r
        lat = global_lats_r(ipt_lats_node_r-1+lan)
        lons_lat = lonsperlar(lat)
        do k=1,km
          call mymaxmin(grid_gr(1,g_pnt+k-1),lons_lat,lonr,1,chr)
        enddo
      enddo
 
      return
      end subroutine do_physics_gridcheck
