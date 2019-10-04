      subroutine aoicpl_prep(deltim,delt_cpl,phour,fhour,idate,
     &                  aoi_fld,global_lats_r,lonsperlar)
!
! Feb 2014 Xingren Wu - Collect Fields for Atm/Ocn/Ice coupling
! Apr 2014 Xingren Wu - Modified: Add SW fluxes for Nir/uv+vis
! Jun 2014 Xingren Wu - Modified: Revise Net SW heat fluxes
! Jul 2014 Xingren Wu - Modified: Add Sea/Land mask
! May 2015 Xingren Wu - Modified: Add t/q/u/v/p/zbot
! Oct 2015 Xingren Wu - Modified: add snow and fixed a bug
!                                 for "fluxes" accumulation
! Jan 2016 P. Tripp   - NUOPC/GSM merge - use exportData instead of array
!
      use resol_def,               ONLY: latr, lonr
      use layout1,                 ONLY: me, nodes, lats_node_r,
     &                                   nodes_comp
      use namelist_physics_def,    ONLY: ens_nam,
     &                                   sfcio_out
      use mpi_def,                 ONLY: mpi_comm_all, 
     &                                   mc_comp, mpi_comm_null
      use gfs_physics_aoi_var_mod, ONLY: AOI_Var_Data
      USE machine,                 ONLY: kind_evod

      implicit none
!!
      TYPE(AOI_Var_Data) :: aoi_fld

      integer              idate(4)
      real(kind=kind_evod) deltim,delt_cpl,phour,fhour

      INTEGER              GLOBAL_LATS_R(LATR),   lonsperlar(LATR)
!!
      integer, parameter ::  noa2oi=555
      integer      ierr
!!
!!
      integer   IOPROC
!!
!-------------------------------------------------------------------------
!     print *,' in aoicpl_prep me=',me
      call mpi_barrier(mpi_comm_all,ierr)
      if(nodes_comp < 1 .or. nodes_comp > nodes) then
        print *, '  NODES_COMP UNDEFINED, CANNOT DO I.O '
        call mpi_finalize()
         stop 333
      endif
!
      ioproc = nodes_comp-1
       
      call MPI_BARRIER(mpi_comm_all,ierr)
! 
!!   Data for Atm/Ocn/Ice coupling (average and instantaneous)
!      print *,'---- start aoi_fld collection section -----'
      if(mc_comp /= MPI_COMM_NULL) then
        if(sfcio_out) then
           call cplflx_collect(IOPROC,noa2oi,
     &                         deltim,delt_cpl,FHOUR,IDATE,
     &                         aoi_fld,global_lats_r,lonsperlar)
        endif
      endif                 ! comp node

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE cplflx_collect(IOPROC,noa2oi,
     &                          deltim,delt_cpl,FHOUR,IDATE,
     &                          aoi_fld,global_lats_r,lonsperlar)
!!
!
      use resol_def,               ONLY: lonr, latr
      use layout1,                 ONLY: me, lats_node_r,lats_node_r_max
      use gfs_physics_aoi_var_mod, ONLY: AOI_Var_Data
      USE machine,             ONLY: kind_evod, kind_io8, kind_io4
      use namelist_physics_def,    ONLY: a2oi_out,fhout
      use module_CPLFIELDS,        ONLY: fillExportFields,NExportFields
      use module_CPLFIELDS,        ONLY: queryFieldList,ExportFieldsList
      use namelist_dynamics_def,   ONLY: wam_ipe_coupling

      implicit none
!!
      TYPE(AOI_Var_Data)        :: aoi_fld
      INTEGER              GLOBAL_LATS_R(LATR)
      INTEGER              lonsperlar(LATR)
      integer   IOPROC
      real(kind=kind_io8) slmsk(LONR,LATS_NODE_R)
!!
      integer j,i,k,noa2oi,ngrid2d,iret
      real (kind=kind_io8) rtime
      real (kind=kind_io8) fhour
      real (kind=kind_evod) deltim,delt_cpl
!
      INTEGER     IDATE(4)
      real (kind=kind_io8)   glolal(lonr,LATS_NODE_R)
      real (kind=kind_io8)   buffo(lonr,LATS_NODE_R)
      real (kind=kind_io8)   buffol(lonr,latr)

      real (kind=kind_io8)   exportData(lonr,latr,NExportFields)

      logical     chk_cplflx
!..................................................................
      integer, dimension(lonr,lats_node_r_max) ::  kmsk, kmsk0
      integer iyr,imo,ida,ihr,itmn,ndig
      character*8 s
      character*40 a2oifname,citmn,cform

!
      if(me == ioproc) print *,'begin cplflx_collect'
!
      IYR     = IDATE(4)
      IMO     = IDATE(2)
      IDA     = IDATE(3)
      IHR     = IDATE(1)
      ITMN    = NINT(FHOUR*60)
      NDIG    = MAX(log10(ITMN+0.5)+1.,2.)
!!
      kmsk0    = 0
      RTIME  = 1./delt_cpl
!..........................................................
!
!     MEAN Zonal compt of momentum flux (N/m**2)
!
      glolal  = aoi_fld%DUSFC*RTIME
      ngrid2d = queryfieldlist(exportFieldsList, 
     &          'mean_zonal_moment_flx')
      if (ngrid2d > 0) then
         CALL uninterpred(1,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
         call unsplit2d_phys_r(ioproc,buffol,buffo,global_lats_r)
         exportData(:,:,ngrid2d)=buffol(:,:)
      endif
!..........................................................
!
!     MEAN Merid compt of momentum flux (N/m**2)
!
      glolal  = aoi_fld%DVSFC*RTIME
      ngrid2d = queryfieldlist(exportFieldsList, 
     &          'mean_merid_moment_flx')
      if (ngrid2d > 0) then
         CALL uninterpred(1,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
         call unsplit2d_phys_r(ioproc,buffol,buffo,global_lats_r)
         exportData(:,:,ngrid2d)=buffol(:,:)
      endif
!..........................................................
!
!     MEAN Sensible heat flux (W/m**2)
!
      glolal  = aoi_fld%DTSFC*RTIME
      ngrid2d = queryfieldlist(exportFieldsList, 
     &          'mean_sensi_heat_flx')
      if (ngrid2d > 0) then
         CALL uninterpred(1,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
         call unsplit2d_phys_r(ioproc,buffol,buffo,global_lats_r)
         exportData(:,:,ngrid2d)=buffol(:,:)
      endif
!..........................................................
!
!     MEAN Latent heat flux (W/m**2)
!
      glolal  = aoi_fld%DQSFC*RTIME
      ngrid2d = queryfieldlist(exportFieldsList, 
     &          'mean_laten_heat_flx')
      if (ngrid2d > 0) then
         CALL uninterpred(1,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
         call unsplit2d_phys_r(ioproc,buffol,buffo,global_lats_r)
         exportData(:,:,ngrid2d)=buffol(:,:)
      endif
!..........................................................
!
!     MEAN Downward long wave radiation flux (W/m**2)
!
      glolal  = aoi_fld%DLWSFC*RTIME
      ngrid2d = queryfieldlist(exportFieldsList, 
     &          'mean_down_lw_flx')
      if (ngrid2d > 0) then
         CALL uninterpred(1,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
         call unsplit2d_phys_r(ioproc,buffol,buffo,global_lats_r)
         exportData(:,:,ngrid2d)=buffol(:,:)
      endif
!..........................................................
!
!     MEAN Downward solar radiation flux (W/m**2)
!
      glolal  = aoi_fld%DSWSFC*RTIME
      ngrid2d = queryfieldlist(exportFieldsList, 
     &          'mean_down_sw_flx')
      if (ngrid2d > 0) then
         CALL uninterpred(1,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
         call unsplit2d_phys_r(ioproc,buffol,buffo,global_lats_r)
         exportData(:,:,ngrid2d)=buffol(:,:)
      endif
!..........................................................
!
!     MEAN Precipitation rate (kg/m**2/s)
!
!     glolal = aoi_fld%GESHEM*1.E3*RTIME
      glolal = aoi_fld%rain*1.E3*RTIME
      ngrid2d = queryfieldlist(exportFieldsList, 
     &          'mean_prec_rate')
      if (ngrid2d > 0) then
         CALL uninterpred(1,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
         call unsplit2d_phys_r(ioproc,buffol,buffo,global_lats_r)
         exportData(:,:,ngrid2d)=buffol(:,:)
      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!..........................................................
!
!     Instataneous Zonal compt of momentum flux (N/m**2)
!
      glolal  = aoi_fld%DUSFCI
      ngrid2d = queryfieldlist(exportFieldsList, 
     &          'inst_zonal_moment_flx')
      if (ngrid2d > 0) then
         CALL uninterpred(1,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
         call unsplit2d_phys_r(ioproc,buffol,buffo,global_lats_r)
         exportData(:,:,ngrid2d)=buffol(:,:)
      endif
!..........................................................
!
!     Instataneous Merid compt of momentum flux (N/m**2)
!
      glolal  = aoi_fld%DVSFCI
      ngrid2d = queryfieldlist(exportFieldsList, 
     &          'inst_merid_moment_flx')
      if (ngrid2d > 0) then
         CALL uninterpred(1,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
         call unsplit2d_phys_r(ioproc,buffol,buffo,global_lats_r)
         exportData(:,:,ngrid2d)=buffol(:,:)
      endif
!..........................................................
!
!     Instataneous Sensible heat flux (W/m**2)
!
      glolal = aoi_fld%DTSFCI
      ngrid2d = queryfieldlist(exportFieldsList, 
     &          'inst_sensi_heat_flx')
      if (ngrid2d > 0) then
         CALL uninterpred(1,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
         call unsplit2d_phys_r(ioproc,buffol,buffo,global_lats_r)
         exportData(:,:,ngrid2d)=buffol(:,:)
      endif
!..........................................................
!
!     Instataneous Latent heat flux (W/m**2)
!
      glolal = aoi_fld%DQSFCI
      ngrid2d = queryfieldlist(exportFieldsList, 
     &          'inst_laten_heat_flx')
      if (ngrid2d > 0) then
         CALL uninterpred(1,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
         call unsplit2d_phys_r(ioproc,buffol,buffo,global_lats_r)
         exportData(:,:,ngrid2d)=buffol(:,:)
      endif
!..........................................................
!
!     Instataneous Downward long wave radiation flux (W/m**2)
!
      glolal = aoi_fld%DLWSFCI
      ngrid2d = queryfieldlist(exportFieldsList, 
     &          'inst_down_lw_flx')
      if (ngrid2d > 0) then
         CALL uninterpred(1,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
         call unsplit2d_phys_r(ioproc,buffol,buffo,global_lats_r)
         exportData(:,:,ngrid2d)=buffol(:,:)
      endif
!..........................................................
!
!     Instataneous Downward solar radiation flux (W/m**2)
!
      glolal = aoi_fld%DSWSFCI
      ngrid2d = queryfieldlist(exportFieldsList, 
     &          'inst_down_sw_flx')
      if (ngrid2d > 0) then
         CALL uninterpred(1,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
         call unsplit2d_phys_r(ioproc,buffol,buffo,global_lats_r)
         exportData(:,:,ngrid2d)=buffol(:,:)
      endif
!..........................................................
!
!     Instataneous Temperature (K) 2 m above ground
!
      glolal = aoi_fld%t2mi
      ngrid2d = queryfieldlist(exportFieldsList, 
     &          'inst_temp_height2m')
      if (ngrid2d > 0) then
         CALL uninterpred(1,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
         call unsplit2d_phys_r(ioproc,buffol,buffo,global_lats_r)
         exportData(:,:,ngrid2d)=buffol(:,:)
      endif
!..........................................................
!
!     Instataneous Specific humidity (kg/kg) 2 m above ground
!
      glolal = aoi_fld%q2mi
      ngrid2d = queryfieldlist(exportFieldsList, 
     &          'inst_spec_humid_height2m')
      if (ngrid2d > 0) then
         CALL uninterpred(1,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
         call unsplit2d_phys_r(ioproc,buffol,buffo,global_lats_r)
         exportData(:,:,ngrid2d)=buffol(:,:)
      endif
!..........................................................
!
!     Instataneous u wind (m/s) 10 m above ground
!
      glolal = aoi_fld%u10mi
      ngrid2d = queryfieldlist(exportFieldsList, 
     &          'inst_zonal_wind_height10m')
      if (ngrid2d > 0) then
         CALL uninterpred(1,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
         call unsplit2d_phys_r(ioproc,buffol,buffo,global_lats_r)
         exportData(:,:,ngrid2d)=buffol(:,:)
      endif
!..........................................................
!
!     Instataneous v wind (m/s) 10 m above ground
!
      glolal = aoi_fld%v10mi
      ngrid2d = queryfieldlist(exportFieldsList, 
     &          'inst_merid_wind_height10m')
      if (ngrid2d > 0) then
         CALL uninterpred(1,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
         call unsplit2d_phys_r(ioproc,buffol,buffo,global_lats_r)
         exportData(:,:,ngrid2d)=buffol(:,:)
      endif

!..........................................................
!
!     Instataneous Temperature (K) at surface
!
      glolal  = aoi_fld%tseai
      ngrid2d = queryfieldlist(exportFieldsList, 
     &          'inst_temp_height_surface')
      if (ngrid2d > 0) then
         CALL uninterpred(1,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
         call unsplit2d_phys_r(ioproc,buffol,buffo,global_lats_r)
         exportData(:,:,ngrid2d)=buffol(:,:)
      endif
!..........................................................
!
!     Instataneous Pressure (Pa) land and sea surface
!
      glolal  = aoi_fld%PSURFI
      ngrid2d = queryfieldlist(exportFieldsList, 
     &          'inst_pres_height_surface')
      if (ngrid2d > 0) then
         CALL uninterpred(1,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
         call unsplit2d_phys_r(ioproc,buffol,buffo,global_lats_r)
         exportData(:,:,ngrid2d)=buffol(:,:)
      endif
!..........................................................
!
!     Instataneous Surface height (m)
!
      glolal = aoi_fld%oro
      ngrid2d = queryfieldlist(exportFieldsList, 
     &          'inst_surface_height')
      if (ngrid2d > 0) then
         CALL uninterpred(1,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
         call unsplit2d_phys_r(ioproc,buffol,buffo,global_lats_r)
         exportData(:,:,ngrid2d)=buffol(:,:)
      endif
!..........................................................
!
!     MEAN NET long wave radiation flux (W/m**2)
!
      glolal  = aoi_fld%NLWSFC*RTIME
      ngrid2d = queryfieldlist(exportFieldsList, 
     &          'mean_net_lw_flx')
      if (ngrid2d > 0) then
         CALL uninterpred(1,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
         call unsplit2d_phys_r(ioproc,buffol,buffo,global_lats_r)
         exportData(:,:,ngrid2d)=buffol(:,:)
      endif
!..........................................................
!
!     MEAN NET solar radiation flux over the ocean (W/m**2)
!
      glolal  = aoi_fld%NSWSFC*RTIME
      ngrid2d = queryfieldlist(exportFieldsList, 
     &          'mean_net_sw_flx')
      if (ngrid2d > 0) then
         CALL uninterpred(1,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
         call unsplit2d_phys_r(ioproc,buffol,buffo,global_lats_r)
         exportData(:,:,ngrid2d)=buffol(:,:)
      endif
!..........................................................
!
!     Instataneous NET long wave radiation flux (W/m**2)
!
      glolal = aoi_fld%NLWSFCI
      ngrid2d = queryfieldlist(exportFieldsList, 
     &          'inst_net_lw_flx')
      if (ngrid2d > 0) then
         CALL uninterpred(1,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
         call unsplit2d_phys_r(ioproc,buffol,buffo,global_lats_r)
         exportData(:,:,ngrid2d)=buffol(:,:)
      endif
!..........................................................
!
!     Instataneous NET solar radiation flux over the ocean (W/m**2)
!
      glolal = aoi_fld%NSWSFCI
      ngrid2d = queryfieldlist(exportFieldsList, 
     &          'inst_net_sw_flx')
      if (ngrid2d > 0) then
         CALL uninterpred(1,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
         call unsplit2d_phys_r(ioproc,buffol,buffo,global_lats_r)
         exportData(:,:,ngrid2d)=buffol(:,:)
      endif
!..........................................................
!
!     MEAN sfc downward nir direct flux (W/m**2)
!
      glolal = aoi_fld%dnirbm*rtime
      ngrid2d = queryfieldlist(exportFieldsList, 
     &          'mean_down_sw_ir_dir_flx')
      if (ngrid2d > 0) then
         CALL uninterpred(1,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
         call unsplit2d_phys_r(ioproc,buffol,buffo,global_lats_r)
         exportData(:,:,ngrid2d)=buffol(:,:)
      endif
!..........................................................
!
!     MEAN sfc downward nir diffused flux (W/m**2)
!
      glolal = aoi_fld%dnirdf*rtime
      ngrid2d = queryfieldlist(exportFieldsList, 
     &          'mean_down_sw_ir_dif_flx')
      if (ngrid2d > 0) then
         CALL uninterpred(1,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
         call unsplit2d_phys_r(ioproc,buffol,buffo,global_lats_r)
         exportData(:,:,ngrid2d)=buffol(:,:)
      endif
!..........................................................
!
!     MEAN sfc downward uv+vis direct flux (W/m**2)
!
      glolal = aoi_fld%dvisbm*rtime
      ngrid2d = queryfieldlist(exportFieldsList, 
     &          'mean_down_sw_vis_dir_flx')
      if (ngrid2d > 0) then
         CALL uninterpred(1,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
         call unsplit2d_phys_r(ioproc,buffol,buffo,global_lats_r)
         exportData(:,:,ngrid2d)=buffol(:,:)
      endif
!..........................................................
!
!     MEAN sfc downward uv+vis diffused flux (W/m**2)
!
      glolal = aoi_fld%dvisdf*rtime
      ngrid2d = queryfieldlist(exportFieldsList, 
     &          'mean_down_sw_vis_dif_flx')
      if (ngrid2d > 0) then
         CALL uninterpred(1,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
         call unsplit2d_phys_r(ioproc,buffol,buffo,global_lats_r)
         exportData(:,:,ngrid2d)=buffol(:,:)
      endif
!..........................................................
!
!     Instataneous sfc downward nir direct flux (W/m**2)
!
      glolal = aoi_fld%dnirbmi
      ngrid2d = queryfieldlist(exportFieldsList, 
     &          'inst_down_sw_ir_dir_flx')
      if (ngrid2d > 0) then
         CALL uninterpred(1,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
         call unsplit2d_phys_r(ioproc,buffol,buffo,global_lats_r)
         exportData(:,:,ngrid2d)=buffol(:,:)
      endif
!..........................................................
!
!     Instataneous sfc downward nir diffused flux (W/m**2)
!
      glolal = aoi_fld%dnirdfi
      ngrid2d = queryfieldlist(exportFieldsList, 
     &          'inst_down_sw_ir_dif_flx')
      if (ngrid2d > 0) then
         CALL uninterpred(1,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
         call unsplit2d_phys_r(ioproc,buffol,buffo,global_lats_r)
         exportData(:,:,ngrid2d)=buffol(:,:)
      endif
!..........................................................
!
!     Instataneous sfc downward uv+vis direct flux (W/m**2)
!
      glolal = aoi_fld%dvisbmi
      ngrid2d = queryfieldlist(exportFieldsList, 
     &          'inst_down_sw_vis_dir_flx')
      if (ngrid2d > 0) then
         CALL uninterpred(1,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
         call unsplit2d_phys_r(ioproc,buffol,buffo,global_lats_r)
         exportData(:,:,ngrid2d)=buffol(:,:)
      endif
!..........................................................
!
!     Instataneous sfc downward uv+vis diffused flux (W/m**2)
!
      glolal = aoi_fld%dvisdfi
      ngrid2d = queryfieldlist(exportFieldsList, 
     &          'inst_down_sw_vis_dif_flx')
      if (ngrid2d > 0) then
         CALL uninterpred(1,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
         call unsplit2d_phys_r(ioproc,buffol,buffo,global_lats_r)
         exportData(:,:,ngrid2d)=buffol(:,:)
      endif
!..........................................................
!
!     MEAN NET downward nir direct flux over the ocean (W/m**2)
!
      glolal = aoi_fld%nnirbm*rtime
      ngrid2d = queryfieldlist(exportFieldsList, 
     &          'mean_net_sw_ir_dir_flx')
      if (ngrid2d > 0) then
         CALL uninterpred(1,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
         call unsplit2d_phys_r(ioproc,buffol,buffo,global_lats_r)
         exportData(:,:,ngrid2d)=buffol(:,:)
      endif
!..........................................................
!
!     MEAN NET downward nir diffused flux over the ocean (W/m**2)
!
      glolal = aoi_fld%nnirdf*rtime
      ngrid2d = queryfieldlist(exportFieldsList, 
     &          'mean_net_sw_ir_dif_flx')
      if (ngrid2d > 0) then
         CALL uninterpred(1,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
         call unsplit2d_phys_r(ioproc,buffol,buffo,global_lats_r)
         exportData(:,:,ngrid2d)=buffol(:,:)
      endif
!..........................................................
!
!     MEAN NET downward uv+vis direct flux over the ocean (W/m**2)
!
      glolal = aoi_fld%nvisbm*rtime
      ngrid2d = queryfieldlist(exportFieldsList, 
     &          'mean_net_sw_vis_dir_flx')
      if (ngrid2d > 0) then
         CALL uninterpred(1,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
         call unsplit2d_phys_r(ioproc,buffol,buffo,global_lats_r)
         exportData(:,:,ngrid2d)=buffol(:,:)
      endif
!..........................................................
!
!     MEAN NET downward uv+vis diffused flux over the ocean (W/m**2)
!
      glolal = aoi_fld%nvisdf*rtime
      ngrid2d = queryfieldlist(exportFieldsList, 
     &          'mean_net_sw_vis_dif_flx')
      if (ngrid2d > 0) then
         CALL uninterpred(1,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
         call unsplit2d_phys_r(ioproc,buffol,buffo,global_lats_r)
         exportData(:,:,ngrid2d)=buffol(:,:)
      endif
!..........................................................
!
!     Instataneous NET downward nir direct flux over the ocean (W/m**2)
!
      glolal = aoi_fld%nnirbmi
      ngrid2d = queryfieldlist(exportFieldsList, 
     &          'inst_net_sw_ir_dir_flx')
      if (ngrid2d > 0) then
         CALL uninterpred(1,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
         call unsplit2d_phys_r(ioproc,buffol,buffo,global_lats_r)
         exportData(:,:,ngrid2d)=buffol(:,:)
      endif
!..........................................................
!
!     Instataneous NET downward nir diffused flux over the ocean (W/m**2)
!
      glolal = aoi_fld%nnirdfi
      ngrid2d = queryfieldlist(exportFieldsList, 
     &          'inst_net_sw_ir_dif_flx')
      if (ngrid2d > 0) then
         CALL uninterpred(1,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
         call unsplit2d_phys_r(ioproc,buffol,buffo,global_lats_r)
         exportData(:,:,ngrid2d)=buffol(:,:)
      endif
!..........................................................
!
!     Instataneous NET downward uv+vis direct flux over the ocean (W/m**2)
!
      glolal = aoi_fld%nvisbmi
      ngrid2d = queryfieldlist(exportFieldsList, 
     &          'inst_net_sw_vis_dir_flx')
      if (ngrid2d > 0) then
         CALL uninterpred(1,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
         call unsplit2d_phys_r(ioproc,buffol,buffo,global_lats_r)
         exportData(:,:,ngrid2d)=buffol(:,:)
      endif
!..........................................................
!
!     Instataneous NET downward uv+vis diffused flux over the ocean (W/m**2)
!
      glolal = aoi_fld%nvisdfi
      ngrid2d = queryfieldlist(exportFieldsList, 
     &          'inst_net_sw_vis_dif_flx')
      if (ngrid2d > 0) then
         CALL uninterpred(1,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
         call unsplit2d_phys_r(ioproc,buffol,buffo,global_lats_r)
         exportData(:,:,ngrid2d)=buffol(:,:)
      endif
!..........................................................
!
!     Land-Sea Mask
!
      glolal = aoi_fld%slimsk
      ngrid2d = queryfieldlist(exportFieldsList, 
     &          'inst_land_sea_mask')
      if (ngrid2d > 0) then
         kmsk = nint(glolal)
         CALL uninterpred(1,kmsk,slmsk,glolal,global_lats_r,lonsperlar)
         buffo = mod(slmsk,2._kind_io8)
      endif
      call unsplit2d_phys_r(ioproc,buffol,buffo,global_lats_r)
      exportData(:,:,ngrid2d)=buffol(:,:)
!..........................................................
!
!     Instataneous Temperature (K) at the lowest GSM level
!
      glolal = aoi_fld%tboti
      ngrid2d = queryfieldlist(exportFieldsList, 
     &          'inst_temp_height_lowest')
      if (ngrid2d > 0) then
         CALL uninterpred(1,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
         call unsplit2d_phys_r(ioproc,buffol,buffo,global_lats_r)
         exportData(:,:,ngrid2d)=buffol(:,:)
      endif
!..........................................................
!
!     Instataneous Specific humidity (kg/kg) at the lowest GSM level
!
      glolal = aoi_fld%qboti
      ngrid2d = queryfieldlist(exportFieldsList, 
     &          'inst_spec_humid_height_lowest')
      if (ngrid2d > 0) then
         CALL uninterpred(1,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
         call unsplit2d_phys_r(ioproc,buffol,buffo,global_lats_r)
         exportData(:,:,ngrid2d)=buffol(:,:)
      endif
!..........................................................
!
!     Instataneous u wind (m/s) at the lowest GSM level
!
      glolal = aoi_fld%uboti
      ngrid2d = queryfieldlist(exportFieldsList, 
     &          'inst_zonal_wind_height_lowest')
      if (ngrid2d > 0) then
         CALL uninterpred(1,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
         call unsplit2d_phys_r(ioproc,buffol,buffo,global_lats_r)
         exportData(:,:,ngrid2d)=buffol(:,:)
      endif
!..........................................................
!
!     Instataneous v wind (m/s) at the lowest GSM level
!
      glolal = aoi_fld%vboti
      ngrid2d = queryfieldlist(exportFieldsList, 
     &          'inst_merid_wind_height_lowest')
      if (ngrid2d > 0) then
         CALL uninterpred(1,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
         call unsplit2d_phys_r(ioproc,buffol,buffo,global_lats_r)
         exportData(:,:,ngrid2d)=buffol(:,:)
      endif
!..........................................................
!
!     Instataneous Pressure (Pa) at the lowest GSM level
!
      glolal  = aoi_fld%pboti
      ngrid2d = queryfieldlist(exportFieldsList, 
     &          'inst_pres_height_lowest')
      if (ngrid2d > 0) then
        CALL uninterpred(1,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
        call unsplit2d_phys_r(ioproc,buffol,buffo,global_lats_r)
        exportData(:,:,ngrid2d)=buffol(:,:)
      endif
!..........................................................
!
!     Instataneous height (m) at the lowest GSM level
!
      glolal  = aoi_fld%zboti
      ngrid2d = queryfieldlist(exportFieldsList, 
     &          'inst_height_lowest')
      if (ngrid2d > 0) then
        CALL uninterpred(1,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
        call unsplit2d_phys_r(ioproc,buffol,buffo,global_lats_r)
        exportData(:,:,ngrid2d)=buffol(:,:)
      endif
!..........................................................
!
!     MEAN snow rate (kg/m**2/s)
!
      glolal = aoi_fld%snow*1.E3*RTIME
      ngrid2d = queryfieldlist(exportFieldsList,
     &          'mean_fprec_rate')
      if (ngrid2d > 0) then
         CALL uninterpred(1,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
         call unsplit2d_phys_r(ioproc,buffol,buffo,global_lats_r)
         exportData(:,:,ngrid2d)=buffol(:,:)
      endif
!..........................................................
      if(me.eq.ioproc) then
         print *,'ioproc,ngrid2d,global_lats_r,',
     &            ioproc,ngrid2d,global_lats_r

         if (a2oi_out) then
           write(s,'(i4,i2.2,i2.2)'),iyr,imo,ida
           write(cform,'("(i",i1,".",i1,")")') ndig,ndig
           write(citmn,cform) itmn
           a2oifname='A2OIf.'//trim(s)//'.'//trim(citmn)
           open(noa2oi,file=a2oifname,form='unformatted',
     &          status='unknown',iostat=iret)
           if (iret /= 0) then
             write(0,*)' iret while opening a2oifname ',iret
           else
             do k=1,NExportFields
                write(noa2oi) exportData(:,:,k)
             enddo
             close(noa2oi)
           endif
           PRINT *,'chk_cplflx,FHOUR,DELTIM,DELT_CPL',
     &              FHOUR,DELTIM,DELT_CPL
         endif
      endif
!!
!..........................................................
! Fill the export Fields for ESMF/NUOPC style coupling
      if (.not. wam_ipe_coupling) then
         call fillExportFields(exportData, lonr, latr, rootPet=ioproc)
      endif
!..........................................................
!!
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

      if(me == ioproc) PRINT *,'cplflx_collect completed'
      RETURN
      END

