       subroutine gloopr (grid_fld, g3d_fld, aoi_fld,                   
     &                    lats_nodes_r,global_lats_r, lonsperlar, phour,
     &                    deltim,xlon,xlat,coszdg,COSZEN,slmsk,weasd,   
     &                    SNCOVR,SNOALB,ZORL,TSEA,HPRIME,SFALB,         
     &                    ALVSF,ALNSF ,ALVWF ,ALNWF,FACSF ,FACWF,CV,CVT,
     &                    CVB  ,SWH,swhc,HLW,hlwc,SFCNSW,SFCDLW,        
     &                    FICE ,TISFC, SFCDSW, sfcemis,                 
     &                    TSFLW,FLUXR, phy_f3d,phy_f2d,                 
     &                    slag,sdec,cdec,NBLCK,KDT                      
!    &,    htrswb,htrlwb                                                &
     &,                   mdl_parm)
!    &     global_times_r)

!! Code Revision:
!! Oct 11 2009     Sarah Lu, grid_gr is replaced by grid_fld
!! Oct 16 2009     Sarah Lu, grid_fld%tracers used
!! Dec 01 2009     Sarah Lu, update fcld (instant cloud cover) in addition
!!                           to cldcov (cumulative cloud cover)
!! Dec 09 2009     Sarah Lu, (1) g3d_fld added to calling argument; (2) grrad
!!                    returns instant cloud cover (cldcov_v); the accumulative 
!!                    and instant cloud cover fields are updated after grrad call
!! Dec 11 2009     Sarah Lu, ldiag3d removed from grrad calling argument
!! Jul/Aug 2009    S. Moorthi Merged with McICA version of radiation from YuTai
!  apr 2012  -     y-t hou  -  modified radiation calling interval
!                     parameters fhswr/fhlwr.  changed units from hour to second
!                     (in multiple of model time step length), thus radiation can
!                     be invoked more frequently than the minimum of 1 hr limit.
!                            -  moved radiation initialization part to
!                     the beginning of model initialization process. in the time
!                     loop, added a subroutine "radupdate" to updae radiation related
!                     external data sets (such as aerosols, gases, astronomy, etc.)
!                     moved computation of coszen to inside subrouine "grrad"
!  nov 2012  -     y-t hou     - modified many physics/radiation control
!                     parameters through module 'physparam' which contains the
!                     use of module 'machine'.
!! Nov 20 2012     Jun Wang  fix the vertical index from levs to levr for gt,gr..etc
!  may 2013        s. mooorthi - removed fpkapx
!  oct 2013        henry Juang  compute prsi from model top with model pressure
!                        thickness to surface for accuracy
!  Mar 2014        Xingren Wu - Add visbm/visdf/nirbm/nirdf to extract for A/O/I coupling
!  Jul 2014        S Moorthi  - merge with GFS and clean up a bit
!  jun 2014        y-t hou     - revised to utilize surface up/down spectral component sw
!                     fluxes for grid composite mean (combined land,ocean,ice,snow), also
!                    estimate ocean values without ice contamination.
!  Aug 2015        S Moorthi - Adding sgs clouds from shoc as an option
!!
!  Mar 2016 J. Han - add ncnvcld3d for conv cloudiness enhancement
!              include shal_cnv, imfdeepcnv, imfshalcnv
!              make new logical lmfshal & lmfdeep2
!  Apr 2016    W.Yang    - bug fix : change  do k = 1, levs  to be  do k = 1, min(levs, levr) 
!                           
!  Mar 2016    Hang Lei  -use physics driver
      use physparam
      use physcons, fv => con_fvirt, rerth => con_rerth,  rk => con_rocp

      use namelist_physics_def,      only : shoc_cld
      use module_radiation_driver,   only : grrad, radupdate
!     use module_radiation_astronomy,only : astronomy
      USE gfs_phy_tracer_config,     only : gfs_phy_tracer

      use module_radsw_parameters,   only : topfsw_type, sfcfsw_type
      use module_radlw_parameters,   only : topflw_type, sfcflw_type
!
!! ---  for optional spectral band heating outputs
      use module_radsw_parameters,   only : nbdsw
      use module_radlw_parameters,   only : nbdlw
!
      use resol_def,             ONLY: levs, levr, latr, lonr, lotgr,   
     &                                 g_t, g_p, g_q, g_dp, g_ps,       
     &                                 ntcw, ntoz, ncld,num_p3d,        
     &                                 nmtvr, ntrac, levp1, nfxr,g_dpdt,
     &                                 lgocart,ntot3d,ntot2d,           
     &                                 npdf3d,ncnvcld3d
      use layout1,               ONLY: me, nodes, lats_node_r,          
     &                                 lats_node_r_max, ipt_lats_node_r
      use gg_def,                ONLY: coslat_r, sinlat_r
      use date_def,              ONLY: idate
      use namelist_physics_def,  ONLY: lsswr,iaer,lslwr,ras,shal_cnv,   
     &                                 lssav, flgmin, ldiag3d,          
     &                                 imfshalcnv, imfdeepcnv,          
     &                                 iovr_lw, iovr_sw, isol, iems,    
     &                                 ialb, fhlwr, fhswr, ico2, ngptc, 
     &                                 crick_proof, norad_precip,ccnorm,
     &                                 ictm, isubc_sw, isubc_lw, fdaer, 
     &                                 sup, ndfi, fhdfi, cplflx
      use d3d_def ,                ONLY: cldcov
      use gfs_physics_gridgr_mod,  ONLY: Grid_Var_Data
      use gfs_physics_g3d_mod,     ONLY: G3D_Var_Data
      use gfs_physics_aoi_var_mod, ONLY: aoi_var_data
      use mersenne_twister,        only: random_setseed, random_index,  
     &                                   random_stat

      use nuopc_physics,
     &   only: state_fields_in, state_fields_out, sfc_properties,
     &         diagnostics, dynamic_parameters,
     &         interface_fields, cloud_properties, radiation_tendencies,
     &         model_parameters, nuopc_rad_run,
     &         nuopc_rad_update,
     &         rad_run_savein, rad_run_readin, rad_run_saveout,
     &         rad_run_readout,
     &         use_nuopc
!
      implicit none
!
      real (kind=kind_phys), parameter :: QMIN =1.0e-10                 
     &,                                   Typical_pgr = 95.0            
     &,                                   cons0 = 0.0,  cons2 = 2.0     
     &,                                   pt00001=1.0e-5
!    &,                                   pt01=0.01
!
!  --- ...  inputs:
      integer, intent(in), dimension(nodes) :: lats_nodes_r
      integer, intent(in), dimension(latr)  :: global_lats_r, lonsperlar

!*    real(kind=kind_grid) grid_gr(lonr*lats_node_r_max,lotgr)

      TYPE(Grid_Var_Data) :: grid_fld 
      TYPE(G3D_Var_Data)  :: g3d_fld 
      type(aoi_var_data)  :: aoi_fld

      integer, intent(in) :: NBLCK


      real (kind=kind_phys), dimension(LONR,LATS_NODE_R), intent(in) :: 
     &                       xlon,  xlat,  slmsk, weasd, zorl,   tsea,  
     &                       alvsf, alnsf, alvwf, alnwf, facsf,  facwf, 
     &                       cv, cvt, cvb, FICE,  tisfc, sncovr, snoalb

      real (kind=kind_phys), intent(inout) ::                           
     &    hprime(NMTVR,LONR,LATS_NODE_R), phour, deltim,                
     &    phy_f3d(NGPTC,LEVS,ntot3d,NBLCK,LATS_NODE_R),                 
     &    phy_f2d(lonr,lats_node_r,ntot2d)
!

      real (kind=kind_phys), intent(inout) ::                           
     &                    fluxr (NFXR,LONR,LATS_NODE_R)

      integer, intent(in) :: KDT
!  --- ...  outputs:
!     real(kind=kind_evod), intent(out) ::                              &
!    &                    global_times_r(latr,NODES)

      real (kind=kind_phys), intent(out), dimension                     
     &          (ngptc,levs,nblck,lats_node_r) :: swh, swhc, hlw, hlwc

      real (kind=kind_phys),dimension(LONR,LATS_NODE_R), intent(out) :: 
     &                    coszdg, coszen, sfcnsw, sfcdlw, tsflw,        
     &                    sfcdsw, SFALB, sfcemis

      real (kind=kind_phys), intent(out) :: slag, sdec, cdec

!! --- ...  optional spectral band heating rates
!     real (kind=kind_phys), optional, intent(out) ::                   &
!    &                 htrswb(NGPTC,LEVS,NBDSW,NBLCK,LATS_NODE_R),      &
!    &                 htrlwb(NGPTC,LEVS,NBDLW,NBLCK,LATS_NODE_R)

!  --- ...  locals:
      real(kind=kind_phys) :: prsl(NGPTC,LEVS),  prslk(NGPTC,LEVS),     
     &                        prsi(NGPTC,LEVP1)

!     real (kind=kind_phys) :: si_loc(LEVR+1)

      real (kind=kind_phys) :: dswcmp(NGPTC,4), uswcmp(NGPTC,4),        
     &                         gt(NGPTC,LEVR),                          
     &                         gr(NGPTC,LEVR), gr1(NGPTC,LEVR,NTRAC-1)

      logical :: lmfshal, lmfdeep2

      real (kind=kind_phys), dimension(ngptc,levs) :: deltaq, cnvw, cnvc
     &,            f_ice, f_rain, r_rime, hlwc_v, swhc_v, cldcov_v, vvel

      real (kind=kind_phys) ::  fluxr_v(NGPTC,NFXR)
      real (kind=kind_phys) ::  work1, work2

      real (kind=kind_phys), dimension(ngptc) :: coszen_v, coszdg_v,    
     &             sinlat_v, coslat_v, hprime_v, flgmin_v

!     real (kind=kind_phys), dimension(LONR,LATS_NODE_R) ::             &
!    &                         sinlat_v, coslat_v

      real (kind=kind_phys) :: rinc(5), dtsw, dtlw, solcon, raddt, solhr
      real (kind=4)         :: rinc4(5)
      integer w3kindreal,w3kindint

      real (kind=kind_phys), save :: facoz

      integer :: njeff, lon, lan, lat, iblk, lons_lat, kk
      integer :: idat(8), jdat(8), DAYS(13), iday, imon, midmon, id
      integer :: nlnsp(LATS_NODE_R)

!  ---  variables of instantaneous calculated toa/sfc radiation fluxes

      type (topfsw_type), dimension(NGPTC) :: topfsw
      type (sfcfsw_type), dimension(NGPTC) :: sfcfsw

      type (topflw_type), dimension(NGPTC) :: topflw
      type (sfcflw_type), dimension(NGPTC) :: sfcflw

!  ---  variables used for random number generator (thread safe mode)
      type (random_stat) :: stat
      integer :: numrdm(LONR*LATR*2), ixseed(LONR,LATS_NODE_R,2)
      integer :: ipseed, icsdlw(NGPTC), icsdsw(NGPTC)
      integer, parameter :: ipsdlim = 1.0e8      ! upper limit for random seeds


!     integer, save :: icwp, k1oz, k2oz, midm, midp, ipsd0, iaerflg

!  ---  number of days in a month
!     data DAYS / 31,28,31,30,31,30,31,31,30,31,30,31,30 /

!  --- ...  control parameters: 
!           (some of the them may be moved into model namelist)

!  ---  ICTM=yyyy#, controls time sensitive external data (e.g. CO2, solcon, aerosols, etc)
!     integer, parameter :: ICTM =   -2 ! same as ICTM=0, but add seasonal cycle
!                                       ! from climatology. no extrapolation.
!     integer, parameter :: ICTM =   -1 ! use user provided external data set for the
!                                       ! forecast time. no extrapolation.
!     integer, parameter :: ICTM =    0 ! use data at initial cond time, if not
!                                       ! available, use latest, no extrapolatio n.
!!    integer, parameter :: ICTM =    1 ! use data at the forecast time, if not
!                                       ! available, use latest and extrapolation.
!     integer, parameter :: ICTM =yyyy0 ! use yyyy data for the forecast time,
!                                       ! no further data extrapolation.
!     integer, parameter :: ICTM =yyyy1 ! use yyyy data for the fcst. if needed, do
!                                       ! extrapolation to match the fcst time.

!  ---  ISOL controls solar constant data source
!!    integer, parameter :: ISOL  = 0  ! use prescribed solar constant
!     integer, parameter :: ISOL  = 1  ! use varying solar const with 11-yr cycle

!  ---  ICO2 controls co2 data source for radiation
!     integer, parameter :: ICO2 = 0   ! prescribed global mean value (old opernl)
!!    integer, parameter :: ICO2 = 1   ! use obs co2 annual mean value only
!     integer, parameter :: ICO2 = 2   ! use obs co2 monthly data with 2-d variation

!  ---  IALB controls surface albedo for sw radiation
!!    integer, parameter :: IALB = 0   ! use climatology alb, based on sfc type
!     integer, parameter :: IALB = 1   ! use modis derived alb (to be developed)

!  ---  IEMS controls surface emissivity and sfc air/ground temp for lw radiation
!        ab: 2-digit control flags. a-for sfc temperature;  b-for emissivity
!!    integer, parameter :: IEMS = 00  ! same air/ground temp; fixed emis = 1.0
!!    integer, parameter :: IEMS = 01  ! same air/ground temp; varying veg typ based emis
!!    integer, parameter :: IEMS = 10  ! diff air/ground temp; fixed emis = 1.0
!!    integer, parameter :: IEMS = 11  ! diff air/ground temp; varying veg typ based emis

!  ---  IAER  controls aerosols scheme selections
!     integer, parameter :: IAER = abc --> abc are three digits integer
!     numbers
!                                          to control aerosol effect
!     a: stratospheric-volcanic forcing;  b: lw radiation;  c: sw
!     radiation.
!      a=0: no stratospheric-volcanic aerosol effect;   =1: include the
!      effect.
!      b=0: no lw tropospheric aerosols; =1: lw compute 1 bnd; =2: lw
!      compute multi bnd.
!      c=0: no sw tropospheric aerosols; =1: sw compute multi band.

!  ---  IOVR controls cloud overlapping method in radiation:
!     integer, parameter :: IOVR_SW = 0  ! sw: random overlap clouds
!!    integer, parameter :: IOVR_SW = 1  ! sw: max-random overlap clouds

!     integer, parameter :: IOVR_LW = 0  ! lw: random overlap clouds
!!    integer, parameter :: IOVR_LW = 1  ! lw: max-random overlap clouds

!  ---  ISUBC controls sub-column cloud approximation in radiation:
!!    integer, parameter :: ISUBC_SW = 0 ! sw: without sub-col clds approx
!     integer, parameter :: ISUBC_SW = 1 ! sw: sub-col clds with
!     prescribed seeds
!     integer, parameter :: ISUBC_SW = 2 ! sw: sub-col clds with random seeds

!!    integer, parameter :: ISUBC_LW = 0 ! lw: without sub-col clds approx
!     integer, parameter :: ISUBC_LW = 1 ! lw: sub-col clds with
!     prescribed seeds
!     integer, parameter :: ISUBC_LW = 2 ! lw: sub-col clds with random seeds

!  ---  iflip indicates model vertical index direction:
!     integer, parameter :: IFLIP = 0    ! virtical profile index from top to bottom
!!    integer, parameter :: IFLIP = 1    ! virtical profile index from bottom to top
!
!    The following parameters are from gbphys
!
!     real (kind=kind_phys), parameter :: dxmax=-16.118095651,          &
!    &                                    dxmin=-9.800790154,           &
!    &                                    dxinv=1.0/(dxmax-dxmin)

!     logical :: change
      logical, save :: first
      data  first / .true. /
!
! NUOPC physics driver types - PT
      type(model_parameters)     :: mdl_parm


! Physics driver containers
      type(state_fields_in)      :: state_fldin
      type(sfc_properties)       :: sfc_prop
      type(diagnostics)          :: diags
      type(interface_fields)     :: intrfc_fld
      type(cloud_properties)     :: cld_prop
      type(radiation_tendencies) :: rad_tend
      type(dynamic_parameters)   :: dyn_parm

      integer :: lonbnd   ! upper lon dimension in lon/lan loop!!
      logical :: savecon  ! Run nuopc save states

!  ---  for debug test use
      real (kind=kind_phys) :: temlon, temlat, alon, alat
      integer               :: ipt
      logical               :: lprnt

      integer   i,j,k,n,item,dbgu

!  ---  timers:
!     real*8 :: rtc, timer1, timer2
!
!===> *** ...  begin here
!!
!!
      idat    = 0
      idat(1) = idate(4)
      idat(2) = idate(2)
      idat(3) = idate(3)
      idat(5) = idate(1)
      rinc    = 0.

! test repro
!     rinc(2) = fhour
      rinc(2) = phour
!     print *,' idate ',idate,' rinc ',rinc
      call w3kind(w3kindreal,w3kindint)
      if(w3kindreal == 4) then
        rinc4 = rinc
        call w3movdat(rinc4, idat, jdat)
      else
        call w3movdat(rinc, idat, jdat)
      endif
!     print *,' jdat ',jdat
!
!===> *** ...  radiation initialization
!
      dtsw  = fhswr                  ! fhswr is in sec
      dtlw  = fhlwr                  ! fhlwr is in sec
      raddt = min(dtsw, dtlw)
!
!  ---  set up parameters for Xu & Randell's cloudiness computation 

      lmfshal = shal_cnv .and. (imfshalcnv > 0)
      lmfdeep2 = (imfdeepcnv == 2)
!
!     if ( fhour > 0.0 ) then
      if ( phour > 0.0 ) then
        solhr = mod(float(jdat(5)),24.0) ! hour after 00z at current fcst time
      else
        solhr = idate(1)                 ! initial time
      endif

!     if (me == 0) write(0,*)' in gloopr solhr=',solhr
      call radupdate                                                    
!  ---  inputs:
     &     ( idat, jdat, dtsw, deltim, lsswr, me,                       
!  ---  outputs:
     &       slag, sdec, cdec, solcon                                   
     &     )

!  --- ...  generate 2-d random seeds array for sub-grid cloud-radiation

      if ( ISUBC_LW == 2 .or. ISUBC_SW == 2 ) then
!       ipseed = mod(nint(100.0*sqrt(fhour*3600)), ipsdlim) + 1 + ipsd0
        ipseed = mod(nint(100.0*sqrt(phour*3600)), ipsdlim) + 1 + ipsd0

        call random_setseed                                             
!  ---  inputs:
     &     ( ipseed,                                                    
!  ---  outputs:
     &       stat                                                       
     &      )
        call random_index                                               
!  ---  inputs:
     &     ( ipsdlim,                                                   
!  ---  outputs:
     &       numrdm, stat                                               
     &     )

        do k = 1, 2
          n = (k-1)*LATR
!$omp parallel do private(i,j,lat,item)
          do j = 1, lats_node_r
            lat  = global_lats_r(ipt_lats_node_r-1+j)
            item = (lat-1)*LONR + n
            do i = 1, LONR
              ixseed(i,j,k) = numrdm(i+item)
            enddo
          enddo
        enddo
      endif
!
!===> *** ...  starting latitude loop
!              ----------------------
!
      do lan=1,lats_node_r
!
         lat      = global_lats_r(ipt_lats_node_r-1+lan)
         lons_lat = lonsperlar(lat)

!$omp parallel do schedule(dynamic,1)
!$omp+private(lon,i,j,k,item,njeff,iblk,n)
!$omp+private(vvel,gt,gr,gr1,work1,work2,flgmin_v)
!$omp+private(cldcov_v,hprime_v,fluxr_v,f_ice,f_rain,r_rime)
!$omp+private(deltaq,cnvw,cnvc)
!$omp+private(coszen_v,coszdg_v,sinlat_v,coslat_v)
!$omp+private(prslk,prsl,prsi,topfsw,sfcfsw,topflw,sfcflw)
!$omp+private(icsdsw,icsdlw)
!$omp+private(hlwc_v,swhc_v)
!$omp+private(lprnt,ipt)
!$omp+private(dswcmp,uswcmp)
!$omp+private(state_fldin,diags,cld_prop,rad_tend)
!$omp+private(intrfc_fld, sfc_prop,dyn_parm)
!$omp+private(lonbnd,savecon)

!!$omp+private(temlon,temlat,lprnt,ipt)
!!$omp+private(lprnt,ipt,dbgu)

!!!$omp+private(temlon,temlat,lprnt,ipt)

        do lon=1,lons_lat,ngptc
!!
          njeff = min(ngptc,lons_lat-lon+1)
          iblk  = (lon-1)/ngptc + 1
!
          lprnt = .false.
          ipt   = 1
!
!  --- ...  for debug test
!         alon = 0.0
!         alat = 2.0
!         alon = 236.25
!         alat = 56.189
!         alon = 22.5
!         alat = -12.381
!         ipt = 0
!         do i = 1, njeff
!           item = lon + i - 1
!           temlon = xlon(item,lan) * 57.29578
!           if (temlon < 0.0) temlon = temlon + 360.0
!           temlat = xlat(item,lan) * 57.29578
!           lprnt = abs(temlon-alon) < 0.5 .and. abs(temlat-alat) < 0.5
!    &          .and. kdt > 0
!           if ( lprnt ) then
!             ipt = i
!             print *,' ipt=',ipt,' lon=',lon,' lan=',lan
!             exit
!           endif
!         enddo
!
          do k = 1, levr
            do i = 1, njeff
              item      = lon+i-1
              gt(i,k)   = grid_fld%t(item,lan,k)
              gr(i,k)   = max(qmin,grid_fld%tracers(1)%flds(item,lan,k))
              prsl(i,k) = grid_fld%p(item,lan,k)
              vvel(i,k) = grid_fld%dpdt(item,lan,k)
            enddo
          enddo

!
!hmhj prsi should be computed from model top for accuray
          do i = 1, njeff
            prsi (i,levs+1) = 0.0
          enddo
          do k = levs, 1, -1
            do i = 1, njeff
              item      = lon+i-1
              prsi(i,k) = prsi(i,k+1) + grid_fld%dp(item,lan,k)
            enddo
          enddo
!     write(1000+me,*)' in gloopr lan=',lan,' prsl=',prsl(njeff,:)
!     write(1000+me,*)' in gloopr lan=',lan,' prsi=',prsi(njeff,:)
!
!       Remaining tracers
!
         do n = 2, ntrac
            item = n - 1
            do k = 1, levr
              do i = 1, njeff
                gr1(i,k,item) = grid_fld%tracers(n)%flds(lon+i-1,lan,k)
              enddo
            enddo
          enddo


          if(num_p3d == 4) then
            if(kdt == 1.or.(kdt == ndfi+1.and.fhdfi == 0.)) then
              do k = 1, min(levs, levr)
                do i = 1, njeff
                  item = lon + i - 1
                  phy_f3d(i,k,1,iblk,lan) = gt(i,k)
                  phy_f3d(i,k,2,iblk,lan) = gr(i,k)
                  phy_f3d(i,k,3,iblk,lan) = gt(i,k)
                  phy_f3d(i,k,4,iblk,lan) = gr(i,k)
                enddo
              enddo
              do i=1,njeff
                item = lon + i - 1
                phy_f2d(item,lan,1) = prsi(i,1)
                phy_f2d(item,lan,2) = prsi(i,1)
              enddo
            endif
          endif
!!
!.....
          if (levr < levs) then
            do i=1,njeff
              prsi(i,levr+1)  = prsi(i,levp1)
              prsl(i,levr)    = (prsi(i,levp1)+prsi(i,levr)) * 0.5
            enddo
          endif
!
          if (ntoz <= 0 .or. iaerflg == 2) then
            do k = 1, levs
              do i = 1, njeff
                prslk(i,k) = (prsl(i,k)*pt00001)**rk
              enddo
            enddo
          endif
          if (shoc_cld) then
            do k = 1, levs
              do i = 1, njeff
                cldcov_v(i,k) = phy_f3d(i,k,ntot3d-2,iblk,lan)
              enddo
!             write(1000+me,*)'sgs_clds=',maxval(cldcov_v(1:njeff,:)),
!    &       ' lan=',lan
            enddo
          endif
!
          do i=1,njeff
            hprime_v(i) = hprime(1,lon+i-1,lan)
          enddo
!
          do k=1,nfxr
            do i=1,njeff
              fluxr_v(i,k) = fluxr(k,lon+i-1,lan)
            enddo
          enddo

          do j = 1, njeff
            sinlat_v(j) = sinlat_r(lat)
            coslat_v(j) = coslat_r(lat)
          enddo

          if (NUM_P3D == 3) then
            do k = 1, LEVR
              do i = 1, njeff
                f_ice (i,k) = phy_f3d(i,k,1,iblk,lan)
                f_rain(i,k) = phy_f3d(i,k,2,iblk,lan)
                r_rime(i,k) = phy_f3d(i,k,3,iblk,lan)
              enddo
            enddo

            work1 = (log(coslat_r(lat)/(lons_lat*latr)) - dxmin) * dxinv
            work1 = max(0.0, min(1.0,work1))
            work2 = flgmin(1)*work1 + flgmin(2)*(1.0-work1)
            do i=1,njeff
              flgmin_v(i) = work2
            enddo
          else
            do i=1,njeff
              flgmin_v(i) = 0.0
            enddo
          endif

          if(num_p3d == 4 .and. npdf3d == 3) then
            do k = 1, levr
              do i = 1, njeff
                deltaq(i,k) = phy_f3d(i,k,5,iblk,lan)
                cnvw(i,k)   = phy_f3d(i,k,6,iblk,lan)
                cnvc(i,k)   = phy_f3d(i,k,7,iblk,lan)
              enddo
            enddo
          else if(npdf3d == 0 .and. ncnvcld3d == 1) then
            kk = num_p3d + 1
            do k = 1, levr
              do i = 1, njeff
                deltaq(i,k) = 0.
                cnvw(i,k)   = phy_f3d(i,k,kk,iblk,lan)
                cnvc(i,k)   = 0.
              enddo
            enddo
          else
            do k = 1, levr
              do i = 1, njeff
                deltaq(i,k) = 0.0
                cnvw(i,k)   = 0.0
                cnvc(i,k)   = 0.0
              enddo
            enddo
          endif

!  --- ...  assign random seeds for sw and lw radiations

          if ( ISUBC_LW==2 .or. ISUBC_SW==2 ) then
            do i = 1, njeff
              icsdsw(i) = ixseed(lon+i-1,lan,1)
              icsdlw(i) = ixseed(lon+i-1,lan,2)
            enddo
          endif
 
!  *** ...  calling radiation driver
!
          if (use_nuopc) then

          lonbnd=lon+njeff-1

          call dyn_parm%setrad(
     &                  xlon(lon:lonbnd,lan), xlat(lon:lonbnd,lan),
     &                  sinlat_v, coslat_v, solhr, NGPTC, njeff, kdt,
     &                  jdat, solcon, icsdsw, icsdlw, dtlw, dtsw,
     &                  lsswr, lslwr, lssav, lmfshal, lmfdeep2, 
     &                  ipt, lprnt, deltim,
     &                  slag, sdec, cdec )

          call state_fldin%setrad(prsi, prsl, prslk,
     &              gt, gr, gr1, vvel )

          call sfc_prop%setrad(
     &      slmsk(lon:lonbnd,lan), tsea(lon:lonbnd,lan),
     &      weasd(lon:lonbnd,lan), sncovr(lon:lonbnd,lan),
     &      snoalb(lon:lonbnd,lan), zorl(lon:lonbnd,lan), hprime_v,
     &      fice(lon:lonbnd,lan), tisfc(lon:lonbnd,lan),
     &      alvsf(lon:lonbnd,lan), alnsf(lon:lonbnd,lan),
     &      alvwf(lon:lonbnd,lan), alnwf(lon:lonbnd,lan),
     &      facsf(lon:lonbnd,lan), facwf(lon:lonbnd,lan)
     &      )

          call diags%setrad(NFXR, fluxr_v, topfsw, topflw,
     &                   dswcmp, uswcmp )

          call cld_prop%setrad(cv(lon:lonbnd,lan), cvt(lon:lonbnd,lan),
     &      cvb(lon:lonbnd,lan), f_ice, f_rain, r_rime, flgmin_v,
     &      cldcov_v, deltaq, sup, cnvw, cnvc )
!!!         cv(lon:lonbnd,lan), cvt(lon:lonbnd,lan),cvb(lon:lonbnd,lan)

          call rad_tend%set(
     &      swh(1:ngptc,1:levr,iblk,lan), sfalb(lon:lonbnd,lan),
     &      coszen_v, hlw(1:ngptc,1:levr,iblk,lan),
     &      tsflw(lon:lonbnd,lan), sfcemis(lon:lonbnd,lan),
     &      coszdg_v )

          call intrfc_fld%setrad(sfcfsw, sfcflw)
!     &       htrlw0, htrsw0, htrswb, htrlwb)

          if (me .eq. 0 .and. kdt .eq. 1) then
            print *, me, ": NUOPC DBG: before nuopc_rad_run in gloopr"
          end if

! PT This is much easier interface - no change needed here
!          if (savecon) then
!            call rad_run_savein(state_fldin, sfc_prop,
!     &                        diags, intrfc_fld,
!     &                        cld_prop, rad_tend, mdl_parm, dyn_parm)
!!! PT TEST - read previously saved inputs
!!            call rad_run_readin(state_fldin, sfc_prop,
!!     &                        diags, intrfc_fld,
!!     &                        cld_prop, rad_tend, mdl_parm, dyn_parm)
!          end if

          call nuopc_rad_run (state_fldin, sfc_prop,
     &                        diags, intrfc_fld,
     &                        cld_prop, rad_tend, mdl_parm, dyn_parm)

          if (me .eq. 0 .and. kdt .eq. 1) then
            print *, me, ": NUOPC DBG: after nuopc_rad_run in gloopr"
          end if

!          if (savecon) then
!            call rad_run_saveout(diags, intrfc_fld,
!     &                        cld_prop, rad_tend)
!          end if

         else

          call grrad
!  ---  inputs:
     &     ( prsi,prsl,prslk,gt,gr,gr1,vvel,slmsk(lon,lan),             
     &       xlon(lon,lan),xlat(lon,lan),tsea(lon,lan),                 
     &       weasd(lon,lan),sncovr(lon,lan),snoalb(lon,lan),            
     &       zorl(lon,lan),hprime_v,                                    
     &       alvsf(lon,lan),alnsf(lon,lan),alvwf(lon,lan),              
     &       alnwf(lon,lan),facsf(lon,lan),facwf(lon,lan),              
     &       fice(lon,lan),tisfc(lon,lan),                              
     &       sinlat_v,coslat_v,solhr,jdat,solcon,                       
     &       cv(lon,lan),cvt(lon,lan),cvb(lon,lan),                     
     &       f_ice,f_rain,r_rime,flgmin_v,                              
     &       icsdsw,icsdlw,NTCW-1,NCLD,NTOZ-1,NTRAC-1,NFXR,             
!    &       dtlw,dtsw, lsswr,lslwr,lssav,                              
     &       dtlw,dtsw, lsswr,lslwr,lssav,shoc_cld,lmfshal,lmfdeep2,    
     &       NGPTC,njeff,LEVR,me, lprnt, ipt, kdt, deltaq,sup,cnvw,cnvc,

!  ---  outputs:
     &       swh(1:ngptc,1:levr,iblk,lan),topfsw,sfcfsw,dswcmp,uswcmp,  
     &       sfalb(lon,lan),coszen_v,coszdg_v,                          
     &       hlw(1:ngptc,1:levr,iblk,lan),topflw,sfcflw,tsflw(lon,lan), 
     &       sfcemis(lon,lan),cldcov_v,                                 
!  ---  input/output:
     &       fluxr_v,                                                    
!    &       fluxr_v, dbgu                                               
!! ---  optional outputs:
     &       htrlw0=hlwc_v,htrsw0=swhc_v                                

!!   &,      HTRSWB=htrswb(1,1,1,iblk,lan),                             
!!   &,      HTRLWB=htrlwb(1,1,1,iblk,lan)                              

     &     )

       end if   ! if use_nuopc

          do j = 1, njeff
            coszen(lon+j-1,lan) = coszen_v(j)
            coszdg(lon+j-1,lan) = coszdg_v(j)
          enddo

          do k=1,nfxr
            do i=1,njeff
              fluxr(k,lon+i-1,lan) = fluxr_v(i,k)
            enddo
          enddo

!  --- ...  radiation fluxes for other physics process or diagnostics

          if (lsswr) then
            do i = 1, njeff
              j = lon + i - 1
              sfcdsw(j,lan) = sfcfsw(i)%dnfxc
              sfcnsw(j,lan) = sfcfsw(i)%dnfxc - sfcfsw(i)%upfxc

!             sfcnirbmd(j,lan)=dswcmp(i,1)
!             sfcnirdfd(j,lan)=dswcmp(i,2)
!             sfcvisbmd(j,lan)=dswcmp(i,3)
!             sfcvisdfd(j,lan)=dswcmp(i,4)
!             sfcnirbmu(j,lan)=uswcmp(i,1)
!             sfcnirdfu(j,lan)=uswcmp(i,2)
!             sfcvisbmu(j,lan)=uswcmp(i,3)
!             sfcvisdfu(j,lan)=uswcmp(i,4)
            enddo
            if (cplflx) then
              do i = 1, njeff
                j = lon + i - 1
                aoi_fld%nirbmdi(j,lan) = dswcmp(i,1)
                aoi_fld%nirdfdi(j,lan) = dswcmp(i,2)
                aoi_fld%visbmdi(j,lan) = dswcmp(i,3)
                aoi_fld%visdfdi(j,lan) = dswcmp(i,4)
                aoi_fld%nirbmui(j,lan) = uswcmp(i,1)
                aoi_fld%nirdfui(j,lan) = uswcmp(i,2)
                aoi_fld%visbmui(j,lan) = uswcmp(i,3)
                aoi_fld%visdfui(j,lan) = uswcmp(i,4)
              enddo

            endif
          endif
          if (lslwr) then
            do i = 1, njeff
              j = lon + i - 1
              sfcdlw(j,lan) = sfcflw(i)%dnfxc
            enddo
            do k=1,levr
              do i=1,njeff
                hlwc(i,k,iblk,lan) = hlwc_v(i,k)
              enddo
            enddo
          endif
          do k=1,levr
            do i=1,njeff
              swhc(i,k,iblk,lan) = swhc_v(i,k)
            enddo
          enddo

          if (levr < levs) then
            do k=levr+1,levs
              do i=1,njeff
                hlw(i,k,iblk,lan)  = hlw(i,levr,iblk,lan)
                hlwc(i,k,iblk,lan) = hlwc(i,levr,iblk,lan)
                swh(i,k,iblk,lan)  = swh(i,levr,iblk,lan)
                swhc(i,k,iblk,lan) = swhc(i,levr,iblk,lan)
              enddo
            enddo
          endif

!     if (lprnt) then
!       write(0,*)' hlw=',hlw(ipt,1:10,iblk,lan),' lan=',lan,' ipt=',ipt,
!    &            ' njeff=',njeff
!       write(0,*)' swh=',swh(ipt,1:10,iblk,lan),' lan=',lan
!     endif
!
! grrad routine computes cldcov_v (instant 3D cloud cover)    -- Sarah Lu
! if ldiag3d is T, update cldcov (accumulative 3D cloud cover)
! if lgocart is T, update fcld (instant 3D cloud cover)

          if (lssav) then
            if (ldiag3d) then
              do k=1,levr
                do i=1,njeff
                  item = lon+i-1
                  cldcov(k,item,lan) = cldcov(k,item,lan)                 
     &                               + cldcov_v(i,k) * raddt
                enddo
              enddo
            endif
          endif
          if (lgocart) then
            do k=1,levr
              do i=1,njeff
                g3d_fld%fcld(lon+i-1,lan,k) = cldcov_v(i,k) 
              enddo
            enddo
          endif

!$$$          write(2900+lat,*) ' ilon = ',lon
!$$$          write(2900+lat,'("swh",T16,"hlw")')
!$$$      do k=1,levs
!$$$         write(2900+lat,
!$$$     .         '(e10.3,T16,e10.3,T31,e10.3)')
!$$$     . swh(1,k,iblk,lan),hlw(1,k,iblk,lan)
!$$$       enddo
!
        enddo                       ! end of lon if loop
!
      enddo                         ! end of lan if loop
!

      return
      end subroutine gloopr
