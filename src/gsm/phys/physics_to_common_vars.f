!     subroutine physics_to_common_vars (psg,ttg,rqg,uug,vvg,
      subroutine physics_to_common_vars (psg,ttg,rqg,
     &                                   ppg,dpg,dpdtg,
     &                                   global_lats_r,lonsperlar)
!!
!! program log:
!!
!! hmhj - this routine change variables from model usage to common
!!        common usage are t=dry temperature (k), p is pascal, real winds
!!        model  usage are t=virtal temperature (k) or enthalpy, 
!!                         p is centibar, mapping winds
!! Mar, 2010 J.WANG - compute mdl phys psg to common variables psg in pascal
!! Jul, 2010 S.Moorthi - removed uug, vvg
!!
      use resol_def,            ONLY: levh, lonr, latr, levs
      use layout1,              ONLY: ipt_lats_node_r, lats_node_r,
     &                                lats_node_r_max
      use gg_def,               ONLY: coslat_r
      use namelist_physics_def, ONLY: gen_coord_hybrid
      USE machine,              ONLY: KIND_GRID, kind_evod
    
      implicit none
!!
!
      integer              global_lats_r(latr)
      integer                 lonsperlar(latr)
!
      REAL(KIND=KIND_GRID)   psg(lonr,lats_node_r_max)
      REAL(KIND=KIND_GRID)   ttg(lonr,lats_node_r_max,levs)
!     REAL(KIND=KIND_GRID)   uug(lonr,lats_node_r_max,levs)
!     REAL(KIND=KIND_GRID)   vvg(lonr,lats_node_r_max,levs)
      REAL(KIND=KIND_GRID)   rqg(lonr,lats_node_r_max,levh)
      REAL(KIND=KIND_GRID)   ppg(lonr,lats_node_r_max,levs)
      REAL(KIND=KIND_GRID)   dpg(lonr,lats_node_r_max,levs)
      REAL(KIND=KIND_GRID) dpdtg(lonr,lats_node_r_max,levs)
!
      integer              i,j,k, nn, nnl
      integer              l,lan,lat
      integer              lons_lat
!
      real(kind=kind_evod), parameter :: cb2pa=1000.
!
!--------------------------------------------------------------------
!
      do lan=1,lats_node_r
        lat = global_lats_r(ipt_lats_node_r-1+lan)
        lons_lat = lonsperlar(lat)
!
!       do k=1,levs
!         do i=1,lons_lat
!           uug(i,lan,k) = uug(i,lan,k) / coslat_r(lat)
!           vvg(i,lan,k) = vvg(i,lan,k) / coslat_r(lat)
!         enddo
!       enddo
!
!-- psg from nems gfs in pascal
        do i=1,lons_lat
          psg(i,lan) =  psg(i,lan) * cb2pa 
        enddo
!
        do k=1,levs
          do i=1,lons_lat
            ppg  (i,lan,k) =  ppg  (i,lan,k) * cb2pa 
            dpg  (i,lan,k) =  dpg  (i,lan,k) * cb2pa 
            dpdtg(i,lan,k) =  dpdtg(i,lan,k) * cb2pa 
          enddo
        enddo
!
      enddo
!
!     print *,' exit model_to_physics_vars '
!!
      return
      end
