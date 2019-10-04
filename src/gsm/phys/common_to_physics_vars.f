!      subroutine common_to_physics_vars(psg,ttg,rqg,uug,vvg,
!     subroutine common_to_physics_vars(psg,ttg,uug,vvg,
      subroutine common_to_physics_vars(psg,ttg,
     &                                  ppg,dpg,dpdtg,
     &                                  global_lats_r,lonsperlar)
!!
!! program log:
!!
!! hmhj - this routine change variables from common usage to model
!!        common usage are t=dry temperature (k), p is pascal, real winds
!!        model  usage are t=virtal temperature (k) or enthalpy, 
!!                         p is centibar, mapping winds
!! Mar, 2010 J.WANG - compute mdl phys psg from common variables psg in pascal
!! Jul, 2010 S.Moorthi - removed uug, vvg 
!!
      use resol_def,            ONLY: lonr, latr, levs, levh
      use layout1,              ONLY: lats_node_r_max, ipt_lats_node_r, 
     &                                lats_node_r
      use gg_def,               ONLY: coslat_r
      use namelist_physics_def, ONLY: gen_coord_hybrid
      USE module_gfs_machine,   ONLY: KIND_GRID, kind_evod
      implicit none
!!
!
      integer              global_lats_r(latr)
      integer                 lonsperlar(latr)
!
      REAL(KIND=KIND_GRID) psg    (lonr,lats_node_r_max)
      REAL(KIND=KIND_GRID) ttg    (lonr,lats_node_r_max,levs)
!     REAL(KIND=KIND_GRID) uug    (lonr,lats_node_r_max,levs)
!     REAL(KIND=KIND_GRID) vvg    (lonr,lats_node_r_max,levs)
!      REAL(KIND=KIND_GRID) rqg    (lonr,lats_node_r_max,levh)
      REAL(KIND=KIND_GRID)   ppg  (lonr,lats_node_r_max,levs)
      REAL(KIND=KIND_GRID)   dpg  (lonr,lats_node_r_max,levs)
      REAL(KIND=KIND_GRID) dpdtg  (lonr,lats_node_r_max,levs)
!
      integer              i,j,k,kk, nn, nnl
      integer              l,lan,lat
      integer              lons_lat
!
      real(kind=kind_evod), parameter :: pa2cb=0.001
!--------------------------------------------------------------------
!
      do lan=1,lats_node_r
        lat = global_lats_r(ipt_lats_node_r-1+lan)
        lons_lat = lonsperlar(lat)
! check
!       print *,' phy: me lan lat lons_lat ',me,lan,lat,lons_lat

!       do k=1,levs
!         do i=1,lons_lat
!           uug(i,lan,k) = uug(i,lan,k) * coslat_r(lat)
!           vvg(i,lan,k) = vvg(i,lan,k) * coslat_r(lat)
!         enddo
!       enddo
!
! get model ps (log surface pressure or surface pressure)
!-- psg from nems gfs commom variabales in pascal
        do i=1,lons_lat
          psg (i,lan) =  psg(i,lan) * pa2cb 
        enddo
!
        do k=1,levs
          do i=1,lons_lat
            ppg  (i,lan,k) =  ppg  (i,lan,k) * pa2cb
            dpg  (i,lan,k) =  dpg  (i,lan,k) * pa2cb
            dpdtg(i,lan,k) =  dpdtg(i,lan,k) * pa2cb 
          enddo
        enddo
!
      enddo
!
!     print *,' exit common_to_model_vars '
!!
      return
      end
