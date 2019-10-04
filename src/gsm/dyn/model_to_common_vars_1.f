      SUBROUTINE model_to_common_vars_1 (uug, vvg, ttg, rqg, psg, 
     &                                 global_lats_a, lonsperlat)
!!
!! hmhj - this routine change variables from model usage to common
!!        common usage are t=dry temperature (k), p is pascal, real winds
!!        model  usage are t=virtal temperature (k) or enthalpy, 
!!                         p is centibar, mapping winds
!!
      USE gfs_dyn_resol_def,     ONLY: lonf, ntrac, thermodyn_id, levs, 
     &    levh, latg
      USE gfs_dyn_layout1,       ONLY: lats_node_a_max, lats_node_a, 
     &    ipt_lats_node_a
      USE gfs_dyn_gg_def,        ONLY: coslat_a
      USE gfs_dyn_machine,       ONLY: KIND_EVOD
      USE gfs_dyn_physcons,      ONLY: fv => con_fvirt
      USE namelist_dynamics_def, ONLY: gen_coord_hybrid
      USE gfs_dyn_tracer_const,  ONLY: cpi
      IMPLICIT none
!!
!
      integer              global_lats_a(latg)
      integer                 lonsperlat(latg)
!
      REAL(KIND = KIND_EVOD) ::  uug(lonf,lats_node_a_max,levs)
      REAL(KIND = KIND_EVOD) ::  vvg(lonf,lats_node_a_max,levs)
      REAL(KIND = KIND_EVOD) ::  ttg(lonf,lats_node_a_max,levs)
      REAL(KIND = KIND_EVOD) ::  rqg(lonf,lats_node_a_max,levh)
      REAL(KIND = KIND_EVOD) ::  psg(lonf,lats_node_a_max)
!
      REAL(KIND = KIND_EVOD) :: tfac(lonf,levs), sumq(lonf,levs)
!
      integer              i, k, nn, nnl
      integer              l,lan,lat
      integer              lons_lat
!
      REAL(KIND = KIND_EVOD), PARAMETER :: qmin=1.0D-10
      REAL(KIND = KIND_EVOD), PARAMETER :: one=1.0D0, cb2pa=1000.0D0
!
!--------------------------------------------------------------------
!
      do lan=1,lats_node_a
        lat = global_lats_a(ipt_lats_node_a-1+lan)
        lons_lat = lonsperlat(lat)

        if (thermodyn_id == 3) then
          do k=1,levs
            do i=1,lons_lat
              tfac(i,k) = 0.0D0
              sumq(i,k) = 0.0D0
            enddo
          enddo
          do nn=1,ntrac
            nnl = (nn-1)*levs
            if (cpi(nn) .ne. 0.0) then
              do k=1,levs
                do i=1,lons_lat
                  sumq(i,k) = sumq(i,k) + rqg(i,lan,nnl+k)
                  tfac(i,k) = tfac(i,k) + 
     &                       cpi(nn)*rqg(i,lan,nnl+k)
                enddo
              enddo
            endif
          enddo
          do k=1,levs
            do i=1,lons_lat
              tfac(i,k) = (one-sumq(i,k))*cpi(0) + tfac(i,k)
            enddo
          enddo
        else
          do k=1,levs
            do i=1,lons_lat
              tfac(i,k) = one + fv*max(rqg(i,lan,k),qmin) 
            enddo
          enddo
        endif
!
!       call mymaxmin(ttg(1,lan,1),lons_lat,lonf,1,' h1 in mdl to com')
!       call mymaxmin(uug(1,lan,1),lons_lat,lonf,1,' u1 in mdl to com')
!       call mymaxmin(vvg(1,lan,1),lons_lat,lonf,1,' v1 in mdl to com')
        do k=1,levs
          do i=1,lons_lat
            ttg(i,lan,k) = ttg(i,lan,k) / tfac(i,k)
            uug(i,lan,k) = uug(i,lan,k) / coslat_a(lat)
            vvg(i,lan,k) = vvg(i,lan,k) / coslat_a(lat)
          enddo
        enddo

!       do k=1,levs
!       call mymaxmin(rqg(1,lan,k),lons_lat,lonf,1,' rq1 ')
!       call mymaxmin(rqg(1,lan,k+levs),lons_lat,lonf,1,' rq2 ')
!       call mymaxmin(rqg(1,lan,k+levs*2),lons_lat,lonf,1,' rq3 ')
!       enddo

!       call mymaxmin(tfac(1,1),lons_lat,lonf,1,' tfac in mdl to com')
!
        if (gen_coord_hybrid) then   ! Ps is the prognostic variable
          do i=1,lons_lat
            psg(i,lan) =  psg(i,lan) * cb2pa 
          enddo
        else                         ! ln(Ps) is the prognostic variable
          do i=1,lons_lat
            psg(i,lan) = exp( psg(i,lan) ) * cb2pa
          enddo
        endif
      enddo
!
!     print *,' exit model_to_common_vars_1 '
!!
      END SUBROUTINE model_to_common_vars_1
