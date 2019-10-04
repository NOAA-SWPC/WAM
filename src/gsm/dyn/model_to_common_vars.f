      subroutine model_to_common_vars (psg,ttg,rqg,uug,vvg,dpg,
     &                                 ppg,dpdtg,
     &                                 global_lats_a,lonsperlat,indxp)
!!
!! hmhj - this routine change variables from model usage to common
!!        common usage are t=dry temperature (k), p is pascal, real winds
!!        model  usage are t=virtal temperature (k) or enthalpy, 
!!                         p is centibar, mapping winds
! program log
! 2011 02 20 : henry juang, update code to have mass_dp options.
!!
      use gfs_dyn_resol_def
      use gfs_dyn_layout1
      use gfs_dyn_gg_def
      use gfs_dyn_vert_def
      use gfs_dyn_date_def
      use namelist_dynamics_def
      use gfs_dyn_coordinate_def 
      use gfs_dyn_tracer_const
      use gfs_dyn_physcons, fv => con_fvirt
      implicit none
!!
!
      integer, dimension(latg) ::  global_lats_a, lonsperlat
      integer                      indxp
!
      real(kind=kind_grid), dimension(lonf,lats_node_a_max) :: psg
      real(kind=kind_grid), dimension(lonf,lats_node_a_max,levs) ::
     &                          ttg, uug, vvg, dpg, ppg, dpdtg
      real(kind=kind_grid), dimension(lonf,lats_node_a_max,levh) ::
     &                          rqg
!
      real(kind=kind_evod), dimension(lonf,levs) :: tfac, sumq
!
      integer              i,j,k, nn, nnl, l, lan, lat, lons_lat
!
      real(kind=kind_evod)  wrk
      real(kind=kind_evod), parameter :: qmin=1.0e-10
     &,                                  one=1.0, cb2pa=1000.
!
!--------------------------------------------------------------------
!
      do lan=1,lats_node_a
        lat      = global_lats_a(ipt_lats_node_a-1+lan)
        lons_lat = lonsperlat(lat)
        wrk      = one / coslat_a(lat)
!
        if (thermodyn_id == 3) then
          do k=1,levs
            do i=1,lons_lat
              tfac(i,k) = 0.0
              sumq(i,k) = 0.0
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
!       call mymaxmin(ttg(1,lan,1),lons_lat,lonf,1,' h1 in mdl to com')
!       call mymaxmin(uug(1,lan,1),lons_lat,lonf,1,' u1 in mdl to com')
!       call mymaxmin(vvg(1,lan,1),lons_lat,lonf,1,' v1 in mdl to com')
        do k=1,levs
          do i=1,lons_lat
            ttg(i,lan,k) = ttg(i,lan,k) / tfac(i,k)
            uug(i,lan,k) = uug(i,lan,k) * wrk
            vvg(i,lan,k) = vvg(i,lan,k) * wrk
            dpg(i,lan,k) = dpg(i,lan,k) * cb2pa 
          enddo
        enddo
!     if (lan == 1) then
!     write(0,*)' dpg=',dpg(11:11,lan,:)
!     endif

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
!         if (semilag) then
!           do i=1,lons_lat
!             psg(i,lan) =  psg(i,lan) * cb2pa 
!           enddo
!         else
            do i=1,lons_lat
              psg(i,lan) = exp( psg(i,lan) ) * cb2pa
            enddo
!         endif
        endif
!
        if( indxp == 1 ) then
          do k=1,levs
            do i=1,lons_lat
              ppg  (i,lan,k) =  ppg  (i,lan,k) * cb2pa 
              dpdtg(i,lan,k) =  dpdtg(i,lan,k) * cb2pa 
            enddo
          enddo
        endif
!
      enddo
!
!     write(0,*)' psg=',psg(1,1),' cb2pa=',cb2pa
!     write(0,*)' ppg=',ppg(1,1,:),' cb2pa=',cb2pa
!     write(150+me,*)' rqg=',rqg(1,1,levh-3:levh),' levh=',levh
!    &,' me=',me
!     print *,' exit model_to_common_vars '
!!
      return
      end
