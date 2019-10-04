!     subroutine getpwatptot (kdt,psg,ttg,rqg, global_lats_a,lonsperlat,
      subroutine getpwatptot (psg,ttg,rqg, global_lats_a,lonsperlat,
     &                        pwat,ptot,ptrc)
!!
!! program log
!!   2007      Henry H. Juang
!!   20100205  J. WANG - this routins is moved from input_fields, it
!!                       computes pwat and ptot
!!   20100825  Sarah Lu - modified to compute tracer global sum
!!
      use gfs_dyn_resol_def
      use gfs_dyn_layout1
      use gfs_dyn_gg_def
      use gfs_dyn_vert_def
      use gfs_dyn_date_def
      use namelist_dynamics_def
      use gfs_dyn_coordinate_def 
      use gfs_dyn_tracer_const
      use gfs_dyn_physcons, fv => con_fvirt, rerth => con_rerth,
     &              grav => con_g,  cp => con_cp , rd => con_rd
      use gfs_dyn_tracer_config, only: glbsum                     !glbsum
      implicit none
!
      integer, dimension(latg), intent(in) :: global_lats_a, lonsperlat
!
      real(kind=kind_evod), intent(in) :: psg(lonf,lats_node_a_max)
      real(kind=kind_evod), intent(in) :: ttg(lonf,lats_node_a_max,levs)
      real(kind=kind_evod), intent(in) :: rqg(lonf,lats_node_a_max,levh)
      real(kind=kind_evod), intent(out):: ptrc(lonf,lats_node_a,ntrac)
      real(kind=kind_evod), dimension(lonf,lats_node_a), intent(out) ::
     &                                                   pwat, ptot

!     integer kdt
!
!local vars
      real(kind=kind_evod)  work(lonf), tkrt0
      real(kind=kind_evod), dimension(lonf,levp1) ::  tki, prsi
      real(kind=kind_evod), dimension(lonf,levs)  ::  tfac, sumq
!
      integer              i, j, k, kk, nn, nnl, l, lan,lat, lons_lat
!
      real(kind=kind_evod), parameter :: qmin=1.0e-10, rkappa = cp / rd
     &,                                  one=1.0, pa2cb=0.001
!
!--------------------------------------------------------------------
!
!     write(1000+me,*)' in getpwatptot kdt=',kdt,' psg=',psg(1,1)
!    &,' ttg=',ttg(1,1,:)
!     if (psg(1,1) < 1 .and. kdt >1) then
!       write(1000+me,*)' calling mpi_quit at kdt=',kdt
!       call mpi_quit(8888)
!     endif

      do lan=1,lats_node_a
        lat      = global_lats_a(ipt_lats_node_a-1+lan)
!       lat      = global_lats_a(ipt_lats_node_a+lats_node_a-lan)
        lons_lat = lonsperlat(lat)
!
!     if (kdt >= 18) then
!     write(1000+me,*)'lan=',lan,' lons_lat=',lons_lat
!    &,' psg=',psg(1:lons_lat,lan),' kdt=',kdt

!     write(1000+me,*)'lan=',lan,' lons_lat=',lons_lat,' psgmin=',
!    &minval(psg(1:lons_lat,lan)),' psgmax=',maxval(psg(1:lons_lat,lan))
!    &,' ipt_lats_node_a=',ipt_lats_node_a
!     endif
! save surface pressure as mass for dry mass adjuctment
! from get model ps (log surface pressure or surface pressure)

        if (gen_coord_hybrid) then   ! Ps is the prognostic variable
          do i=1,lons_lat
            ptot (i,lan) =  psg(i,lan)
          enddo
!                           get pressure at interfaces for pwat 
          tki = 0.0
          do k=2,levs
            do i=1,lons_lat
              tkrt0     = (ttg(i,lan,k-1)+ttg(i,lan,k))
     &                               /(thref(k-1)+thref(k))
              tkrt0     = tkrt0**rkappa
              tki (i,k) = ck5(k)*tkrt0
            enddo
          enddo
          do k=1,levp1
            do i=1,lons_lat
              prsi(i,k)  = ak5(k) + bk5(k)*ptot(i,lan) + tki(i,k) 
            enddo
          enddo
        else if (hybrid) then    ! ln(Ps) is the prognostic variable
          do i=1,lons_lat
!           if (psg(i,lan) < 2.0.or. psg(i,lan) > 5.0)
!    & write(1000+me,*)' kdt=',kdt,'i=',i,' psg=',psg(i,lan),' lan=',lan
            ptot(i,lan) = exp(psg(i,lan))
          enddo
          do k=1,levp1
            kk=levp1+1-k
            do i=1,lons_lat
              prsi(i,k)  = ak5(kk) + bk5(kk)*ptot(i,lan)
            enddo
          enddo
        else
          do i=1,lons_lat
            ptot(i,lan) = exp(psg(i,lan))
          enddo
          do k=1,levp1
            do i=1,lons_lat
              prsi(i,k)  = si(k)*ptot(i,lan)
            enddo
          enddo
        endif                      

!       call mymaxmin(psg(1,lan),lons_lat,lonf,1,' psg in com to mdl')
!
! get pwat (total vertical integrated water)
        do i=1,lons_lat
          pwat(i,lan) = 0.0
        enddo
        do k=1,levs
          do i=1,lons_lat
            work(i) = 0.0
          enddo
          if (ncld > 0) then
            do nn=ntcw,ntcw+ncld-1
              nnl = (nn-1)*levs
              do i=1,lons_lat
                work(i) = work(i) + rqg(i,lan,nnl+k)
              enddo
            enddo
          endif
          do i=1,lons_lat
            pwat(i,lan) = pwat(i,lan) + (prsi(i,k)-prsi(i,k+1))
     &                                * (rqg(i,lan,k) + work(i))
          enddo
        enddo

!
! get ptrc (tracer global sum)                                   !glbsum
!
        if( glbsum ) then                                        !glbsum
          do nn = 1, ntrac                                       !glbsum
            nnl = (nn-1)*levs                                    !glbsum
            do i=1,lons_lat                                      !glbsum
             ptrc(i,lan,nn) = 0.0                                !glbsum
             do k=1,levs                                         !glbsum
               ptrc(i,lan,nn) = ptrc(i,lan,nn) +                 !glbsum
     &         (prsi(i,k)-prsi(i,k+1))*rqg(i,lan,nnl+k)          !glbsum
             enddo                                               !glbsum
            enddo                                                !glbsum
          enddo                                                  !glbsum
        endif                                                    !glbsum

!
!       call mymaxmin(rqg(1,lan,1),lons_lat,lonf,1,' rqg in com to mdl')
!       call mymaxmin(pwat(1,lan),lons_lat,lonf,1,' pwat in com to mdl')

      enddo           ! end of lan loop
!
!!
      return
      end
